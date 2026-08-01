"""Camera / photo decoding.

`decode_frame` takes a raw camera frame (BGR numpy array, as OpenCV returns)
and tries to read a CMBS code from it:

  1. downscale the frame,
  2. locate the four 6x6 finder patterns via multi-scale template matching
     (all four rotations, so any code orientation works),
  3. label them as the code's TL / TR / BR / BL corners by their cyclic order
     around the code centroid (rotation invariant; the exact starting corner
     does not matter because the decoder tries all four rotations),  4. build a homography from the 144 known finder-cell correspondences
     (36 per finder) with RANSAC and warp the code to a square 32x32 grid,
  5. sample each cell's centre colour and run the normal decoder.

If the finder search fails, a plain bounding-box fallback is attempted
(useful when the code fills the frame against a uniform background).
"""

import cv2
import numpy as np

from .palette import PALETTE, nearest_color_index
from .codec import GRID_SIZE, FINDER_SIZE, FINDER_PATTERN, decode_payload, DecodeError

try:
    _HAS_CV2 = True
except ImportError:  # pragma: no cover
    _HAS_CV2 = False


def has_camera_support():
    return _HAS_CV2


# Finders are stored in RGB in the palette; OpenCV works in BGR.
_PALETTE_BGR = [(b, g, r) for (r, g, b) in PALETTE]

_FINDER_SCALES = (3, 4, 5, 6, 8, 10, 12)
_FINDER_MATCH_THRESHOLD = 0.5


def _render_finder(cell):
    img = np.zeros((FINDER_SIZE * cell, FINDER_SIZE * cell, 3), np.uint8)
    for r in range(FINDER_SIZE):
        for c in range(FINDER_SIZE):
            img[r * cell:(r + 1) * cell, c * cell:(c + 1) * cell] = _PALETTE_BGR[FINDER_PATTERN[r][c]]
    return img


_templates = None


def _get_templates():
    """Pre-render the finder pattern at several scales x 4 rotations."""
    global _templates
    if _templates is None:
        _templates = []
        rot_codes = {1: cv2.ROTATE_90_CLOCKWISE,
                     2: cv2.ROTATE_180, 3: cv2.ROTATE_90_COUNTERCLOCKWISE}
        for cell in _FINDER_SCALES:
            base = _render_finder(cell)
            for rot in range(4):
                t = base if rot == 0 else cv2.rotate(base, rot_codes[rot])
                _templates.append((cell, rot, t))
    return _templates


def _sample_grid(warped, size=GRID_SIZE):
    """Sample a warped code image into a grid of palette indices (centre 60% of each cell)."""
    h, w = warped.shape[:2]
    grid = [[0] * size for _ in range(size)]
    for r in range(size):
        for c in range(size):
            x0 = c * w // size
            x1 = (c + 1) * w // size
            y0 = r * h // size
            y1 = (r + 1) * h // size
            ix0 = x0 + (x1 - x0) // 5
            ix1 = x1 - (x1 - x0) // 5
            iy0 = y0 + (y1 - y0) // 5
            iy1 = y1 - (y1 - y0) // 5
            if ix1 <= ix0:
                ix1 = ix0 + 1
            if iy1 <= iy0:
                iy1 = iy0 + 1
            patch = warped[iy0:iy1, ix0:ix1]
            avg = patch.reshape(-1, 3).mean(axis=0)
            grid[r][c] = nearest_color_index((int(avg[2]), int(avg[1]), int(avg[0])))
    return grid


def _detect_finders(img):
    """Locate the 4 finder patterns.  Returns [TL, TR, BR, BL] match tuples
    (score, x, y, cell, rot) or None."""
    h, w = img.shape[:2]
    raw = []
    for cell, rot, t in _get_templates():
        if t.shape[0] >= h or t.shape[1] >= w:
            continue
        res = cv2.matchTemplate(img, t, cv2.TM_CCOEFF_NORMED)
        ys, xs = np.where(res >= _FINDER_MATCH_THRESHOLD)
        for x, y in zip(xs, ys):
            raw.append((float(res[y, x]), int(x), int(y), cell, rot))
    raw.sort(key=lambda m: -m[0])

    kept = []
    for m in raw:
        if all(abs(m[1] - k[1]) > 24 or abs(m[2] - k[2]) > 24 for k in kept):
            kept.append(m)
    if len(kept) < 4:
        return None

    groups = []
    for m in kept:
        for g in groups:
            if abs(m[1] - g[0][1]) < 30 and abs(m[2] - g[0][2]) < 30:
                g.append(m)
                break
        else:
            groups.append([m])
    groups.sort(key=lambda g: -max(x[0] for x in g))
    if len(groups) < 4:
        return None

    top = [max(g, key=lambda x: x[0]) for g in groups[:4]]

    def center(m):
        return m[1] + 2.5 * m[3], m[2] + 2.5 * m[3]

    centers = [center(m) for m in top]
    cx = sum(p[0] for p in centers) / 4
    cy = sum(p[1] for p in centers) / 4

    # Order the four markers cyclically around the code centroid; for screen
    # coordinates (y down) the ascending angle order is TL -> TR -> BR -> BL
    # up to a rotation.  The starting corner does not matter: decode_payload
    # tries all four rotations, and the homography absorbs each marker's own
    # orientation via its image-order cell correspondences.
    order = sorted(range(4), key=lambda i: np.arctan2(centers[i][1] - cy, centers[i][0] - cx))
    tl, tr, br, bl = order
    return [top[tl], top[tr], top[br], top[bl]]


def _warp_from_finders(img, finders):
    """Build a homography from the finder cell correspondences and warp to a square grid."""
    src = []
    dst = []
    for label, (score, x, y, cell, rot) in zip(("TL", "TR", "BR", "BL"), finders):
        fx, fy = {
            "TL": (0, 0),
            "TR": (GRID_SIZE - FINDER_SIZE, 0),
            "BR": (GRID_SIZE - FINDER_SIZE, GRID_SIZE - FINDER_SIZE),
            "BL": (0, GRID_SIZE - FINDER_SIZE),
        }[label]
        for r in range(FINDER_SIZE):
            for c in range(FINDER_SIZE):
                # Image-order cells: the homography (not the template rotation)
                # absorbs each marker's own orientation, so the cyclic label
                # assignment stays a pure rotation for any code orientation.
                src.append((x + c * cell, y + r * cell))
                dst.append(((fx + c) * 8, (fy + r) * 8))
    matrix, _ = cv2.findHomography(
        np.array(src, dtype=np.float32), np.array(dst, dtype=np.float32), cv2.RANSAC, 3.0)
    if matrix is None:
        return None
    return cv2.warpPerspective(img, matrix, (GRID_SIZE * 8, GRID_SIZE * 8))


def _bbox_corners(img):
    """Fallback: locate a roughly square, code-shaped region against a uniform background."""
    h, w = img.shape[:2]
    bg = img[:8, :, :].reshape(-1, 3).mean(axis=0)
    diff = img.astype(int) - bg
    mask = np.sum(diff * diff, axis=2) > 900  # squared distance ~ (30)^2
    ys, xs = np.where(mask)
    if len(xs) < 64:
        return None
    x0, y0, x1, y1 = int(xs.min()), int(ys.min()), int(xs.max()), int(ys.max())
    bw, bh = x1 - x0, y1 - y0
    if bw < 32 or bh < 32 or bh / bw < 0.7 or bh / bw > 1.4:
        return None
    if bw * bh < 0.06 * w * h:
        return None
    return np.array([[x0, y0], [x1, y0], [x1, y1], [x0, y1]], dtype=np.float32)


_BBOX_DST = np.array([[0, 0], [32, 0], [32, 32], [0, 32]], dtype=np.float32) * 8


def _warp_bbox(img, src):
    matrix = cv2.getPerspectiveTransform(src, _BBOX_DST)
    return cv2.warpPerspective(img, matrix, (GRID_SIZE * 8, GRID_SIZE * 8))


def decode_frame(frame_bgr, max_size=640):
    """Try to decode a CMBS code from one BGR camera frame.

    Returns (payload_bytes, rotation) or None if nothing was found.
    """
    if not _HAS_CV2:
        return None
    h, w = frame_bgr.shape[:2]
    scale = min(1.0, max_size / max(h, w))
    if scale < 1.0:
        img = cv2.resize(frame_bgr, (int(w * scale), int(h * scale)),
                         interpolation=cv2.INTER_AREA)
    else:
        img = frame_bgr

    warped = None
    finders = _detect_finders(img)
    if finders is not None:
        warped = _warp_from_finders(img, finders)
    if warped is None:
        src = _bbox_corners(img)
        if src is not None:
            warped = _warp_bbox(img, src)
    if warped is None:
        return None

    grid = _sample_grid(warped)
    try:
        payload, rot = decode_payload(grid)
        return payload, rot
    except DecodeError:
        return None
