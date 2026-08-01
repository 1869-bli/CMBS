"""CMBS codec: grid layout, format blocks, bit/cell mapping, encode/decode.

Physical layout of a 32x32 grid:

  * four 6x6 finder patterns in the four corners (diagonal colour gradient,
    used to validate orientation); the bottom-right pattern is the alignment
    marker that anchors the homography under perspective,
  * two replicated 8-cell format blocks on the left edge (magic, ECC level,
    checksum) that also tell the decoder which way up the code is,
  * everything else is data: 864 cells = 324 bytes.

3 bits fit in one cell and 3 bytes fill exactly 8 cells, so the byte stream
maps to the grid with no fractional cells.  Data is laid down in a serpentine
(boustrophedon) scan so localised damage spreads across the Reed-Solomon
blocks instead of clumping.
"""

from .reed_solomon import rs

GRID_SIZE = 32
BITS_PER_CELL = 3
FINDER_SIZE = 6
FINDER_POS = (
    (0, 0),
    (0, GRID_SIZE - FINDER_SIZE),
    (GRID_SIZE - FINDER_SIZE, 0),
    (GRID_SIZE - FINDER_SIZE, GRID_SIZE - FINDER_SIZE),
)
FORMAT_BLOCKS = ((7, 1), (24, 1))  # (top-left row, top-left col), 8 cells wide
FORMAT_CELLS = 8

ECC_LEVELS = ("L", "M", "Q")

# Usable data cells (324 bytes = 864 cells; the BR finder reserves 36).
DATA_CELLS_TOTAL = 864
DATA_CELLS = 864
DATA_BYTES = DATA_CELLS * BITS_PER_CELL // 8  # 324

# n_blocks x (data bytes, ecc bytes) per ECC level.
_RS_PARAMS = {
    "L": dict(n_blocks=2, data=137, ecc=25),
    "M": dict(n_blocks=2, data=120, ecc=42),
    "Q": dict(n_blocks=2, data=106, ecc=56),
}

FORMAT_MAGIC = 0xDA

# Finder pattern: 6x6 diagonal gradient using all 8 colours.
FINDER_PATTERN = [[(r + c) % 8 for c in range(FINDER_SIZE)] for r in range(FINDER_SIZE)]


class DecodeError(Exception):
    pass


# -- reserved cells -------------------------------------------------------

def _is_reserved(r, c):
    for fr, fc in FINDER_POS:
        if fr <= r < fr + FINDER_SIZE and fc <= c < fc + FINDER_SIZE:
            return True
    for fr, fc in FORMAT_BLOCKS:
        if r == fr and fc <= c < fc + FORMAT_CELLS:
            return True
    return False


_data_cells = []
for _r in range(GRID_SIZE):
    _row = [c for c in range(GRID_SIZE) if not _is_reserved(_r, c)]
    if _r % 2 == 1:
        _row.reverse()
    _data_cells.extend((_r, c) for c in _row)


# -- bit / cell helpers ----------------------------------------------------

def _bits_to_cells(bits):
    return [
        bits[i * 3] * 4 + bits[i * 3 + 1] * 2 + bits[i * 3 + 2]
        for i in range(len(bits) // 3)
    ]


def _cells_to_bits(values):
    bits = []
    for v in values:
        bits.append((v >> 2) & 1)
        bits.append((v >> 1) & 1)
        bits.append(v & 1)
    return bits


def _bytes_to_cells(stream, n_cells):
    bits = []
    for b in stream:
        bits.append((b >> 7) & 1)
        bits.append((b >> 6) & 1)
        bits.append((b >> 5) & 1)
        bits.append((b >> 4) & 1)
        bits.append((b >> 3) & 1)
        bits.append((b >> 2) & 1)
        bits.append((b >> 1) & 1)
        bits.append(b & 1)
    while len(bits) % 3:
        bits.append(0)
    cells = _bits_to_cells(bits)
    if len(cells) < n_cells:
        cells += [7] * (n_cells - len(cells))
    return cells[:n_cells]


def _cells_to_bytes(values):
    bits = _cells_to_bits(values)
    out = bytearray()
    for i in range(0, len(bits), 8):
        b = 0
        for j in range(8):
            b = (b << 1) | bits[i + j]
        out.append(b)
    return bytes(out)


# -- format block -----------------------------------------------------------

def payload_capacity(level):
    level = _check_level(level)
    return _RS_PARAMS[level]["n_blocks"] * _RS_PARAMS[level]["data"] - 2


def _check_level(level):
    level = str(level).upper()
    if level not in ECC_LEVELS:
        raise ValueError("ECC level must be one of %s" % ", ".join(ECC_LEVELS))
    return level


def _format_bits(level):
    bits = [0] * 24
    for i in range(8):
        bits[i] = (FORMAT_MAGIC >> (7 - i)) & 1
    e = ECC_LEVELS.index(level)
    bits[8] = (e >> 1) & 1
    bits[9] = e & 1
    # bits 10..11: version = 0
    chk = sum(bits[:20]) & 0xF
    for i in range(4):
        bits[20 + i] = (chk >> (3 - i)) & 1
    return bits


def _read_format(bits):
    magic = int("".join(str(b) for b in bits[:8]), 2)
    if magic != FORMAT_MAGIC:
        return None
    e = bits[8] * 2 + bits[9]
    if e not in (0, 1, 2):
        return None
    chk = int("".join(str(b) for b in bits[20:24]), 2)
    if sum(bits[:20]) & 0xF != chk:
        return None
    return ECC_LEVELS[e]


def _read_format_block(grid):
    levels = []
    for fr, fc in FORMAT_BLOCKS:
        bits = _cells_to_bits([grid[fr][fc + i] for i in range(FORMAT_CELLS)])
        level = _read_format(bits)
        if level is not None:
            levels.append(level)
    if not levels:
        return None
    if len(set(levels)) == 1:
        return levels[0]
    return levels[0]


def _finders_ok(grid):
    """Check that enough of the finder patterns are intact.

    A warped camera frame rarely samples every finder cell perfectly, so a
    finder counts once at least 75% of its cells match and the grid is
    accepted when at least two finders survive.  Wrong rotations only ever
    match a small fraction of cells, so orientation disambiguation is
    unaffected.
    """
    good = 0
    need = FINDER_SIZE * FINDER_SIZE * 3 // 4  # 27 of 36
    for fr, fc in FINDER_POS:
        matches = sum(1 for r in range(FINDER_SIZE) for c in range(FINDER_SIZE)
                      if grid[fr + r][fc + c] == FINDER_PATTERN[r][c])
        if matches >= need:
            good += 1
    return good >= 2


def _rotate_clockwise(grid):
    return [[grid[GRID_SIZE - 1 - c][r] for c in range(GRID_SIZE)] for r in range(GRID_SIZE)]


# -- encode -----------------------------------------------------------------

def build_grid(payload, level="M"):
    """Encode `payload` bytes into a 32x32 grid of palette indices."""
    level = _check_level(level)
    p = _RS_PARAMS[level]
    total_data = p["n_blocks"] * p["data"]
    if len(payload) > total_data - 2:
        raise ValueError(
            "payload of %d bytes exceeds capacity (%d bytes) for ECC level %s"
            % (len(payload), total_data - 2, level)
        )
    raw = len(payload).to_bytes(2, "big") + payload
    raw += bytes(total_data - len(raw))

    blocks = [raw[i * p["data"]:(i + 1) * p["data"]] for i in range(p["n_blocks"])]
    ecc_blocks = [rs.encode(list(b), p["ecc"]) for b in blocks]

    stream = bytearray()
    for i in range(p["data"]):
        for b in blocks:
            stream.append(b[i])
    for i in range(p["ecc"]):
        for b in ecc_blocks:
            stream.append(b[i])

    cells = _bytes_to_cells(bytes(stream), DATA_CELLS)

    grid = [[7] * GRID_SIZE for _ in range(GRID_SIZE)]

    for fr, fc in FINDER_POS:
        for r in range(FINDER_SIZE):
            for c in range(FINDER_SIZE):
                grid[fr + r][fc + c] = FINDER_PATTERN[r][c]

    fbits = _format_bits(level)
    fcells = _bits_to_cells(fbits)
    for fr, fc in FORMAT_BLOCKS:
        for i in range(FORMAT_CELLS):
            grid[fr][fc + i] = fcells[i]

    for i, (r, c) in enumerate(_data_cells[:DATA_CELLS]):
        grid[r][c] = cells[i]

    return grid


# -- decode -----------------------------------------------------------------

def parse_grid(grid):
    """Read the 324 data bytes back out of a grid (no error correction)."""
    values = [grid[r][c] for (r, c) in _data_cells[:DATA_CELLS]]
    return _cells_to_bytes(values)


def _extract(grid, level):
    p = _RS_PARAMS[level]
    stream = parse_grid(grid)
    total_data = p["n_blocks"] * p["data"]
    data_stream = stream[:total_data]
    ecc_stream = stream[total_data:]

    blocks = []
    for b in range(p["n_blocks"]):
        data_b = data_stream[b::p["n_blocks"]]
        ecc_b = ecc_stream[b::p["n_blocks"]]
        codeword = list(data_b) + list(ecc_b)
        fixed, ok = rs.decode(codeword, p["ecc"])
        if not ok:
            raise DecodeError("too many errors to correct")
        blocks.append(fixed[:p["data"]])

    raw = b"".join(bytes(b) for b in blocks)
    length = int.from_bytes(raw[:2], "big")
    if length > total_data - 2:
        raise DecodeError("bad length prefix")
    return raw[2:2 + length]


def decode_payload(grid):
    """Decode a grid, trying all 4 rotations.  Returns (payload, rotation)."""
    g = grid
    for rot in range(4):
        level = _read_format_block(g)
        if level is not None and _finders_ok(g):
            try:
                return _extract(g, level), rot
            except DecodeError:
                pass
        g = _rotate_clockwise(g)
    raise DecodeError("could not decode: not a valid CMBS code")
