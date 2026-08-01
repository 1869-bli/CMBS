"""Round-trip tests for the CMBS codec.

Run directly:  python test_cmbs.py
Or via pytest: pytest test_cmbs.py
"""

import os
import sys
import tempfile

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from cmbs import (
    build_grid,
    decode_payload,
    payload_capacity,
    encode_to_image,
    decode_image,
    DecodeError,
)
from cmbs.codec import _rotate_clockwise, _data_cells, parse_grid, DATA_CELLS
from cmbs.reed_solomon import rs

TESTS = []


def test(fn):
    TESTS.append(fn)
    return fn


@test
def test_payload_sizes():
    assert payload_capacity("L") == 272
    assert payload_capacity("M") == 238
    assert payload_capacity("Q") == 210
    for level in "LMQ":
        assert payload_capacity(level) >= 200


@test
def test_roundtrip_all_levels():
    for level in "LMQ":
        cap = payload_capacity(level)
        for size in (0, 1, 2, 3, 33, 100, cap):
            data = bytes((i * 7 + size) & 0xFF for i in range(size))
            grid = build_grid(data, level)
            out, rot = decode_payload(grid)
            assert out == data, (level, size)
            assert rot == 0


@test
def test_payload_too_big():
    try:
        build_grid(b"x" * (payload_capacity("Q") + 1), "Q")
    except ValueError:
        pass
    else:
        raise AssertionError("expected ValueError for oversized payload")


@test
def test_unicode_roundtrip():
    text = "Hello, CMBS! \u00e9\u00e8\u4f60\u597d\U0001f600\U0001f534\U0001f7e2"
    data = text.encode("utf-8")
    grid = build_grid(data, "M")
    out, _ = decode_payload(grid)
    assert out.decode("utf-8") == text


@test
def test_rotation():
    data = b"which way is up?"
    grid = build_grid(data, "M")
    g = grid
    for rot in (0, 1, 2, 3):
        if rot:
            g = _rotate_clockwise(g)
        out, found = decode_payload(g)
        assert out == data, (rot, found)
        assert found == (4 - rot) % 4


@test
def test_rs_roundtrip():
    msg = list(bytes(range(143)))
    ecc = rs.encode(msg, 25)
    codeword = msg + ecc
    fixed, ok = rs.decode(codeword, 25)
    assert ok and fixed == codeword


@test
def test_rs_corrects_errors():
    msg = list(bytes(range(143)))
    ecc = rs.encode(msg, 25)
    codeword = msg + ecc
    codeword[0] ^= 0xFF
    codeword[5] ^= 0x55
    codeword[100] ^= 0x01
    codeword[167] ^= 0x80
    fixed, ok = rs.decode(codeword, 25)
    assert ok
    assert fixed[:143] == msg


@test
def test_rs_too_many_errors():
    msg = list(bytes(range(126)))
    ecc = rs.encode(msg, 42)
    codeword = msg + ecc
    for i in range(0, 30):
        codeword[i] ^= 0xFF
    fixed, ok = rs.decode(codeword, 42)
    assert not ok or fixed is None


@test
def test_corrupted_cells():
    import random

    rng = random.Random(1234)
    data = bytes(range(100))
    grid = build_grid(data, "Q")
    # flip 20 random cells (up to 20 corrupted bytes spread over 2 RS blocks)
    for _ in range(20):
        r = rng.randrange(32)
        c = rng.randrange(32)
        if (r < 6 and c < 6) or (r < 6 and c >= 26) or (r >= 26 and c < 6):
            continue
        grid[r][c] = (grid[r][c] + 1 + rng.randrange(7)) % 8
    out, _ = decode_payload(grid)
    assert out == data


@test
def test_serpentine_covers_grid():
    from cmbs.codec import _is_reserved

    seen = set(_data_cells)
    assert len(_data_cells) == 864
    assert len(seen) == 864
    reserved = 0
    for r in range(32):
        for c in range(32):
            if _is_reserved(r, c):
                reserved += 1
            else:
                assert (r, c) in seen, (r, c)
    assert reserved == 160


@test
def test_image_roundtrip():
    data = b"image round trip \xf0\x9f\x9a\x80"
    with tempfile.TemporaryDirectory() as tmp:
        path = os.path.join(tmp, "code.png")
        img = encode_to_image(data, "M", cell=20, border=4)
        img.save(path)
        out, rot = decode_image(path)
        assert out == data
        assert rot == 0


@test
def test_image_scaled():
    data = b"scaled code"
    with tempfile.TemporaryDirectory() as tmp:
        path = os.path.join(tmp, "code.png")
        img = encode_to_image(data, "Q", cell=30, border=6)
        img = img.resize((img.width // 2, img.height // 2), __import__("PIL").Image.BILINEAR)
        img.save(path)
        out, _ = decode_image(path)
        assert out == data


@test
def test_image_rotated():
    data = b"rotated on paper"
    with tempfile.TemporaryDirectory() as tmp:
        path = os.path.join(tmp, "code.png")
        encode_to_image(data, "M", cell=24, border=5).save(path)
        from PIL import Image

        img = Image.open(path).rotate(90, expand=True, fillcolor=(255, 255, 255))
        p2 = os.path.join(tmp, "rot.png")
        img.save(p2)
        out, rot = decode_image(p2)
        assert out == data
        assert rot == 1


@test
def test_data_cell_count():
    assert DATA_CELLS == 864
    # 864 cells * 3 bits = 2592 bits = 324 bytes
    grid = build_grid(b"", "M")
    stream = parse_grid(grid)
    assert len(stream) == 324


@test
def test_garbage_grid():
    import random

    rng = random.Random(9)
    grid = [[rng.randrange(8) for _ in range(32)] for _ in range(32)]
    try:
        decode_payload(grid)
    except DecodeError:
        pass
    else:
        raise AssertionError("expected DecodeError for garbage grid")


def _synthetic_photo(data, level="M", rotate=None, perspective=False):
    """Build a fake camera frame: a rendered code on a plain background."""
    import cv2
    import numpy as np
    from PIL import Image

    img = np.array(encode_to_image(data, level, cell=10, border=4).convert("RGB"))
    img = img[:, :, ::-1].copy()  # RGB -> BGR
    h, w = img.shape[:2]
    canvas = np.full((600, 760, 3), 205, np.uint8)
    x0, y0 = 170, 120
    canvas[y0:y0 + h, x0:x0 + w] = img

    if perspective:
        src = np.array([[x0, y0], [x0 + w, y0], [x0 + w, y0 + h], [x0, y0 + h]], dtype=np.float32)
        dst = np.array([
            [x0 + 18, y0 - 14],
            [x0 + w - 26, y0 + 22],
            [x0 + w + 8, y0 + h - 18],
            [x0 - 20, y0 + h + 16],
        ], dtype=np.float32)
        canvas = cv2.warpPerspective(canvas, cv2.getPerspectiveTransform(src, dst),
                                     (canvas.shape[1], canvas.shape[0]),
                                     borderValue=(205, 205, 205))

    if rotate is not None:
        canvas = cv2.rotate(canvas, rotate)
    return canvas


@test
def test_camera_plain():
    try:
        import cv2
        from cmbs.camera import decode_frame
    except ImportError:
        return
    data = b"camera test 123"
    out = decode_frame(_synthetic_photo(data))
    assert out is not None
    assert out[0] == data


@test
def test_camera_rotations():
    try:
        import cv2
        from cmbs.camera import decode_frame
    except ImportError:
        return
    for rotate in (cv2.ROTATE_90_CLOCKWISE, cv2.ROTATE_180,
                   cv2.ROTATE_90_COUNTERCLOCKWISE):
        data = b"rotated camera frame"
        frame = _synthetic_photo(data, rotate=rotate)
        out = decode_frame(frame)
        assert out is not None, "decode failed for rotation %s" % rotate
        assert out[0] == data


@test
def test_camera_perspective():
    try:
        import cv2
        from cmbs.camera import decode_frame
    except ImportError:
        return
    data = b"keystone perspective code"
    out = decode_frame(_synthetic_photo(data, perspective=True))
    assert out is not None, "perspective decode failed"
    assert out[0] == data


@test
def test_camera_empty_frame():
    try:
        import cv2
        import numpy as np
        from cmbs.camera import decode_frame
    except ImportError:
        return
    blank = np.full((480, 640, 3), 200, np.uint8)
    assert decode_frame(blank) is None


def main():
    import traceback

    failed = 0
    for fn in TESTS:
        try:
            fn()
            print("PASS  %s" % fn.__name__)
        except Exception:
            failed += 1
            print("FAIL  %s" % fn.__name__)
            traceback.print_exc()
    print()
    print("%d passed, %d failed" % (len(TESTS) - failed, failed))
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
