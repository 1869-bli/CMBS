"""CMBS - Colour Manifested Byte Storage.

Encode bytes as a 32x32 grid of 8-colour cells (3 bits each), with finder
patterns for orientation and Reed-Solomon error correction, like a colour QR
code.
"""

from .palette import PALETTE, COLOR_NAMES, COLOR_GLYPHS, nearest_color_index
from .codec import (
    GRID_SIZE,
    ECC_LEVELS,
    payload_capacity,
    build_grid,
    decode_payload,
    parse_grid,
    DecodeError,
)
from .image import render_image, encode_to_image, image_to_grid, decode_image

__version__ = "0.1.0"

__all__ = [
    "PALETTE",
    "COLOR_NAMES",
    "COLOR_GLYPHS",
    "nearest_color_index",
    "GRID_SIZE",
    "ECC_LEVELS",
    "payload_capacity",
    "build_grid",
    "decode_payload",
    "parse_grid",
    "render_image",
    "encode_to_image",
    "image_to_grid",
    "decode_image",
    "DecodeError",
    "__version__",
]
