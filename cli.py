"""CMBS command-line interface.

Examples:
  python cli.py encode -o hello.png "Hello, CMBS!"
  python cli.py encode -o note.png -e Q -t note.txt
  python cli.py decode hello.png
  python cli.py info
"""

import argparse
import sys

from cmbs import (
    __version__,
    encode_to_image,
    decode_image,
    payload_capacity,
    COLOR_NAMES,
    COLOR_GLYPHS,
    PALETTE,
    GRID_SIZE,
    ECC_LEVELS,
    DecodeError,
)
from cmbs.codec import _RS_PARAMS


def cmd_info(args):
    print("CMBS v%s" % __version__)
    print("Colour Manifested Byte Storage - a colour QR-style code")
    print()
    print("Grid:        %dx%d (%d cells) x 3 bits/cell = %.1f bytes" % (
        GRID_SIZE, GRID_SIZE, GRID_SIZE * GRID_SIZE, GRID_SIZE * GRID_SIZE * 3 / 8))
    print("Reserved:    4x 6x6 finder patterns + 2x 8-cell format blocks")
    print("Data cells:  %d (%d bytes), serpentine scan" % (864, 864 * 3 // 8))
    print("Error corr.: Reed-Solomon over GF(256), 2 blocks, interleaved")
    print()
    print("%-4s %-28s %-10s %s" % ("Lvl", "data/ecc bytes per block", "payload", "corrects / block"))
    for level in ECC_LEVELS:
        p = _RS_PARAMS[level]
        print("%-4s %-28s %-10d %d byte-errors" % (
            level,
            "%d + %d" % (p["data"], p["ecc"]),
            payload_capacity(level),
            p["ecc"] // 2,
        ))
    print()
    print("Palette (3-bit value -> colour):")
    for i, name in enumerate(COLOR_NAMES):
        r, g, b = PALETTE[i]
        print("   %3d  %-8s (%d, %d, %d)" % (i, name, r, g, b))


def cmd_encode(args):
    if args.text:
        data = " ".join(args.text).encode("utf-8")
    else:
        data = open(args.file, "rb").read() if args.file else sys.stdin.buffer.read()
    cap = payload_capacity(args.ecc)
    if len(data) > cap:
        sys.stderr.write(
            "error: %d bytes is too big for ECC level %s (capacity %d bytes)\n"
            % (len(data), args.ecc, cap))
        return 1
    img = encode_to_image(data, args.ecc, cell=args.cell, border=args.border)
    img.save(args.out)
    print("wrote %s: %d bytes at ECC level %s (capacity %d)" % (
        args.out, len(data), args.ecc, cap))
    return 0


def cmd_decode(args):
    try:
        data, rot = decode_image(args.path)
    except DecodeError as exc:
        sys.stderr.write("error: %s\n" % exc)
        return 1
    try:
        text = data.decode("utf-8")
    except UnicodeDecodeError:
        text = None
    if args.bytes:
        sys.stdout.buffer.write(data)
    elif text is not None and "\x00" not in text and text.isprintable():
        try:
            print(text)
        except UnicodeEncodeError:
            print("raw bytes (%d): %s" % (len(data), data.hex()))
    else:
        print("raw bytes (%d): %s" % (len(data), data.hex()))
    print("(rotation %d)" % rot, file=sys.stderr)
    return 0


def main(argv=None):
    parser = argparse.ArgumentParser(prog="cmbs", description="CMBS colour QR code encoder/decoder")
    sub = parser.add_subparsers(dest="cmd")

    e = sub.add_parser("encode", help="encode text/bytes into a CMBS image")
    e.add_argument("-o", "--out", required=True, help="output PNG path")
    e.add_argument("-e", "--ecc", choices=list(ECC_LEVELS), default="M",
                   help="error correction level (default M)")
    e.add_argument("-c", "--cell", type=int, default=20, help="pixels per cell (default 20)")
    e.add_argument("-b", "--border", type=int, default=4, help="quiet zone in cells (default 4)")
    e.add_argument("-t", "--file", help="read payload from this file (binary)")
    e.add_argument("text", nargs="*", help="payload text (omit to read stdin)")
    e.set_defaults(func=cmd_encode)

    d = sub.add_parser("decode", help="decode a CMBS image")
    d.add_argument("path", help="input PNG path")
    d.add_argument("--bytes", action="store_true", help="write raw bytes to stdout")
    d.set_defaults(func=cmd_decode)

    i = sub.add_parser("info", help="print the CMBS format specification")
    i.set_defaults(func=cmd_info)

    args = parser.parse_args(argv)
    if not getattr(args, "func", None):
        parser.print_help()
        return 0
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
