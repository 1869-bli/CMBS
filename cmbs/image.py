"""CMBS image helpers (Pillow): render a grid to a PNG and sample it back.

Rendering adds a white "quiet zone" border.  Reading crops to the bounding
box of pixels that differ from the corner colour (the finder patterns always
pin all four edges of the grid, so the bounding box lands on the code), then
samples the centre 60% of each 32x32 cell and maps the average colour to the
nearest palette entry.
"""

from .palette import PALETTE, nearest_color_index
from .codec import GRID_SIZE, decode_payload

from PIL import Image, ImageDraw


def render_image(grid, cell=20, border=4):
    """Render a palette-index grid to a Pillow Image with a white border."""
    n = GRID_SIZE
    size = n * cell
    img = Image.new("RGB", (size + 2 * border * cell, size + 2 * border * cell), (255, 255, 255))
    draw = ImageDraw.Draw(img)
    off = border * cell
    for r in range(n):
        for c in range(n):
            x0 = off + c * cell
            y0 = off + r * cell
            draw.rectangle([x0, y0, x0 + cell - 1, y0 + cell - 1], fill=PALETTE[grid[r][c]])
    return img


def encode_to_image(payload, level="M", cell=20, border=4):
    from .codec import build_grid
    return render_image(build_grid(payload, level), cell=cell, border=border)


def _content_bbox(img):
    w, h = img.size
    corner = img.getpixel((0, 0))
    px = img.load()
    thr = 625  # squared distance (25^2) from the corner colour
    minx, miny, maxx, maxy = w, h, -1, -1
    for y in range(h):
        for x in range(w):
            r, g, b = px[x, y]
            dr, dg, db = r - corner[0], g - corner[1], b - corner[2]
            if dr * dr + dg * dg + db * db > thr:
                if x < minx:
                    minx = x
                if x > maxx:
                    maxx = x
                if y < miny:
                    miny = y
                if y > maxy:
                    maxy = y
    if maxx < minx or maxy < miny:
        raise ValueError("no CMBS code found in image")
    return minx, miny, maxx + 1, maxy + 1


def _average_color(region):
    return region.resize((1, 1), Image.BILINEAR).getpixel((0, 0))


def image_to_grid(image):
    """Read a Pillow Image and sample it into a 32x32 grid of palette indices."""
    img = image.convert("RGB")
    box = _content_bbox(img)
    crop = img.crop(box)
    cw, ch = crop.size
    grid = [[0] * GRID_SIZE for _ in range(GRID_SIZE)]
    for r in range(GRID_SIZE):
        for c in range(GRID_SIZE):
            x0 = c * cw // GRID_SIZE
            x1 = (c + 1) * cw // GRID_SIZE
            y0 = r * ch // GRID_SIZE
            y1 = (r + 1) * ch // GRID_SIZE
            ix0 = x0 + (x1 - x0) // 5
            ix1 = x1 - (x1 - x0) // 5
            iy0 = y0 + (y1 - y0) // 5
            iy1 = y1 - (y1 - y0) // 5
            if ix1 <= ix0:
                ix1 = ix0 + 1
            if iy1 <= iy0:
                iy1 = iy0 + 1
            avg = _average_color(crop.crop((ix0, iy0, ix1, iy1)))
            grid[r][c] = nearest_color_index(avg)
    return grid


def decode_image(image_or_path):
    """Open (or accept) an image, sample it and decode the payload."""
    if isinstance(image_or_path, str):
        image = Image.open(image_or_path)
    else:
        image = image_or_path
    grid = image_to_grid(image)
    payload, rot = decode_payload(grid)
    return payload, rot
