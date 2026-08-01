"""CMBS colour palette.

8 colours, one per 3-bit value (000..111), tuned so the colours are as
separated as possible in RGB space while keeping enough contrast against a
white background for a camera / scanner to disambiguate them in normal
lighting.  Ordering matches the original spec: black, red, green, blue,
yellow, purple(magenta), orange, white.
"""

PALETTE = (
    (0, 0, 0),          # 000 black
    (200, 30, 30),      # 001 red
    (40, 160, 70),      # 010 green
    (20, 90, 220),      # 011 blue
    (250, 200, 30),     # 100 yellow
    (170, 30, 180),     # 101 purple / magenta
    (240, 130, 20),     # 110 orange
    (245, 245, 245),    # 111 white
)

COLOR_NAMES = (
    "black",
    "red",
    "green",
    "blue",
    "yellow",
    "purple",
    "orange",
    "white",
)

COLOR_GLYPHS = (
    "\u26ab",
    "\U0001f534",
    "\U0001f7e2",
    "\U0001f535",
    "\U0001f7e1",
    "\U0001f7e3",
    "\U0001f7e0",
    "\u26aa",
)


def nearest_color_index(rgb):
    """Return palette index whose colour is closest to `rgb` (squared RGB distance)."""
    r, g, b = rgb
    best = 0
    best_d = None
    for i, (pr, pg, pb) in enumerate(PALETTE):
        dr, dg, db = pr - r, pg - g, pb - b
        d = dr * dr + dg * dg + db * db
        if best_d is None or d < best_d:
            best_d = d
            best = i
    return best
