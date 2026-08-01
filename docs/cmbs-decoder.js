/* CMBS decoder in JavaScript.
 *
 * Pure-JS port of the Python cmbs package (palette, Reed-Solomon, grid
 * decode) plus a camera pipeline (finder detection, homography warp,
 * bounding-box fallback) so codes can be read in the browser.
 *
 * Works in a browser (window.CMBS) and in Node (module.exports).
 */

var CMBS = (function () {
  'use strict';

  // ---------------------------------------------------------------- palette
  var PALETTE = [
    [0, 0, 0],          // 000 black
    [200, 30, 30],      // 001 red
    [40, 160, 70],      // 010 green
    [20, 90, 220],      // 011 blue
    [250, 200, 30],     // 100 yellow
    [170, 30, 180],     // 101 purple / magenta
    [240, 130, 20],     // 110 orange
    [245, 245, 245],    // 111 white
  ];

  function nearestIndex(r, g, b) {
    var best = 0, bestD = Infinity;
    for (var i = 0; i < 8; i++) {
      var pr = PALETTE[i][0], pg = PALETTE[i][1], pb = PALETTE[i][2];
      var dr = pr - r, dg = pg - g, db = pb - b;
      var d = dr * dr + dg * dg + db * db;
      if (d < bestD) { bestD = d; best = i; }
    }
    return best;
  }

  // ------------------------------------------------------------- Reed-Solomon
  var RS_EXP = new Uint8Array(512);
  var RS_LOG = new Uint8Array(256);
  (function () {
    var x = 1;
    for (var i = 0; i < 255; i++) {
      RS_EXP[i] = x;
      RS_LOG[x] = i;
      x <<= 1;
      if (x & 0x100) x ^= 0x11D;
    }
    for (var i = 255; i < 512; i++) RS_EXP[i] = RS_EXP[i - 255];
  })();

  function rsMul(a, b) {
    if (a === 0 || b === 0) return 0;
    return RS_EXP[RS_LOG[a] + RS_LOG[b]];
  }
  function rsDiv(a, b) {
    if (b === 0) throw new Error('division by zero in GF(256)');
    if (a === 0) return 0;
    return RS_EXP[((RS_LOG[a] - RS_LOG[b]) % 255 + 255) % 255];
  }
  function rsPow(base, exp) {
    if (base === 0) return 0;
    return RS_EXP[(RS_LOG[base] * exp) % 255];
  }
  function evalPoly(poly, x) {
    var y = 0;
    for (var i = 0; i < poly.length; i++) y = rsMul(y, x) ^ poly[i];
    return y;
  }
  function genPoly(n) {
    var g = [1];
    for (var i = 0; i < n; i++) {
      var root = RS_EXP[i];
      var newg = new Array(g.length + 1).fill(0);
      for (var k = 0; k < g.length; k++) {
        newg[k + 1] ^= rsMul(g[k], root);
        newg[k] ^= g[k];
      }
      g = newg;
    }
    return g;
  }
  function rsEncode(msg, nEcc) {
    var g = genPoly(nEcc);
    var res = Array.prototype.slice.call(msg);
    for (var i = 0; i < nEcc; i++) res.push(0);
    for (var i = 0; i < msg.length; i++) {
      var coef = res[i];
      if (coef) {
        for (var j = 1; j < g.length; j++) res[i + j] ^= rsMul(g[j], coef);
      }
    }
    return res.slice(msg.length);
  }
  function syndromes(msg, nEcc) {
    var out = [];
    for (var i = 0; i < nEcc; i++) out.push(evalPoly(msg, RS_EXP[i]));
    return out;
  }
  function berlekampMassey(synd) {
    var n = synd.length;
    var lamb = [1], b = [1], L = 0, m = 1, bb = 1;
    for (var i = 0; i < n; i++) {
      var d = synd[i];
      for (var j = 1; j <= L; j++) d ^= rsMul(lamb[j], synd[i - j]);
      if (d === 0) {
        m++;
      } else if (2 * L <= i) {
        var t = lamb.slice();
        var coef = rsDiv(d, bb);
        while (lamb.length < m + b.length) lamb.push(0);
        for (var j = 0; j < b.length; j++) lamb[j + m] ^= rsMul(coef, b[j]);
        L = i + 1 - L;
        b = t;
        bb = d;
        m = 1;
      } else {
        var coef2 = rsDiv(d, bb);
        while (lamb.length < m + b.length) lamb.push(0);
        for (var j = 0; j < b.length; j++) lamb[j + m] ^= rsMul(coef2, b[j]);
        m++;
      }
    }
    while (lamb.length > 1 && lamb[lamb.length - 1] === 0) lamb.pop();
    return lamb;
  }
  function chien(lamb, msgLen) {
    var deg = lamb.length - 1;
    var lrev = lamb.slice().reverse();
    var positions = [];
    for (var j = 0; j < msgLen; j++) {
      var x = RS_EXP[(255 - j) % 255];
      if (evalPoly(lrev, x) === 0) positions.push(msgLen - 1 - j);
    }
    if (positions.length !== deg) return null;
    return positions;
  }
  function forney(synd, lamb, positions, nEcc, codeword) {
    var omega = new Array(nEcc).fill(0);
    for (var i = 0; i < synd.length; i++) {
      var s = synd[i];
      if (s === 0) continue;
      for (var j = 0; j < lamb.length; j++) {
        var c = lamb[j];
        if (c === 0) continue;
        var deg = i + j;
        if (deg < nEcc) omega[deg] ^= rsMul(s, c);
      }
    }
    var omegaRev = omega.slice().reverse();
    var msgLen = codeword.length;
    for (var pi = 0; pi < positions.length; pi++) {
      var idx = positions[pi];
      var j = msgLen - 1 - idx;
      var x = RS_EXP[j];
      var xinv = RS_EXP[(255 - j) % 255];
      var omegaX = evalPoly(omegaRev, xinv);
      var deriv = 0;
      for (var k = 1; k < lamb.length; k++) {
        if (k & 1) deriv ^= rsMul(lamb[k], rsPow(xinv, k - 1));
      }
      var magnitude = rsMul(x, rsDiv(omegaX, deriv));
      codeword[idx] ^= magnitude;
    }
  }
  function rsDecode(codeword, nEcc) {
    var msg = Array.prototype.slice.call(codeword);
    var synd = syndromes(msg, nEcc);
    var allZero = true;
    for (var i = 0; i < synd.length; i++) if (synd[i] !== 0) { allZero = false; break; }
    if (allZero) return msg;
    var lamb = berlekampMassey(synd);
    var deg = lamb.length - 1;
    if (deg * 2 > nEcc) return null;
    var positions = chien(lamb, msg.length);
    if (positions === null) return null;
    forney(synd, lamb, positions, nEcc, msg);
    synd = syndromes(msg, nEcc);
    allZero = true;
    for (var i = 0; i < synd.length; i++) if (synd[i] !== 0) { allZero = false; break; }
    if (allZero) return msg;
    return null;
  }

  // ------------------------------------------------------------------ codec
  var GRID_SIZE = 32;
  var FINDER_SIZE = 6;
  var FINDER_POS = [[0, 0], [0, 26], [26, 0], [26, 26]];
  var FORMAT_BLOCKS = [[7, 1], [24, 1]];
  var FORMAT_CELLS = 8;
  var FORMAT_MAGIC = 0xDA;
  var DATA_CELLS = 864;
  var RS_PARAMS = {
    L: { n_blocks: 2, data: 137, ecc: 25 },
    M: { n_blocks: 2, data: 120, ecc: 42 },
    Q: { n_blocks: 2, data: 106, ecc: 56 },
  };

  function isReserved(r, c) {
    for (var i = 0; i < FINDER_POS.length; i++) {
      var fr = FINDER_POS[i][0], fc = FINDER_POS[i][1];
      if (r >= fr && r < fr + FINDER_SIZE && c >= fc && c < fc + FINDER_SIZE) return true;
    }
    for (var i = 0; i < FORMAT_BLOCKS.length; i++) {
      var fr2 = FORMAT_BLOCKS[i][0], fc2 = FORMAT_BLOCKS[i][1];
      if (r === fr2 && c >= fc2 && c < fc2 + FORMAT_CELLS) return true;
    }
    return false;
  }

  var DATA_CELLS_LIST = [];
  (function () {
    for (var r = 0; r < GRID_SIZE; r++) {
      var row = [];
      for (var c = 0; c < GRID_SIZE; c++) if (!isReserved(r, c)) row.push(c);
      if (r % 2 === 1) row.reverse();
      for (var i = 0; i < row.length; i++) DATA_CELLS_LIST.push([r, row[i]]);
    }
  })();

  function cellsToBits(values) {
    var bits = [];
    for (var i = 0; i < values.length; i++) {
      var v = values[i];
      bits.push((v >> 2) & 1, (v >> 1) & 1, v & 1);
    }
    return bits;
  }
  function cellsToBytes(values) {
    var bits = cellsToBits(values);
    var out = new Uint8Array(bits.length >> 3);
    for (var i = 0; i < out.length; i++) {
      var b = 0;
      for (var j = 0; j < 8; j++) b = (b << 1) | bits[i * 8 + j];
      out[i] = b;
    }
    return out;
  }
  function readFormat(bits) {
    var magic = 0;
    for (var i = 0; i < 8; i++) magic = magic * 2 + bits[i];
    if (magic !== FORMAT_MAGIC) return null;
    var e = bits[8] * 2 + bits[9];
    if (e < 0 || e > 2) return null;
    var chk = 0;
    for (var i = 20; i < 24; i++) chk = chk * 2 + bits[i];
    var sum = 0;
    for (var i = 0; i < 20; i++) sum += bits[i];
    if ((sum & 0xF) !== chk) return null;
    return ['L', 'M', 'Q'][e];
  }
  function readFormatBlock(grid) {
    var level = null;
    for (var i = 0; i < FORMAT_BLOCKS.length; i++) {
      var fr = FORMAT_BLOCKS[i][0], fc = FORMAT_BLOCKS[i][1];
      var vals = [];
      for (var k = 0; k < FORMAT_CELLS; k++) vals.push(grid[fr][fc + k]);
      var l = readFormat(cellsToBits(vals));
      if (l !== null) level = l;
    }
    return level;
  }
  function findersOk(grid) {
    var good = 0, need = FINDER_SIZE * FINDER_SIZE * 3 / 4;
    for (var i = 0; i < FINDER_POS.length; i++) {
      var fr = FINDER_POS[i][0], fc = FINDER_POS[i][1];
      var matches = 0;
      for (var r = 0; r < FINDER_SIZE; r++)
        for (var c = 0; c < FINDER_SIZE; c++)
          if (grid[fr + r][fc + c] === ((r + c) % 8)) matches++;
      if (matches >= need) good++;
    }
    return good >= 2;
  }
  function rotateGrid(g) {
    var out = [];
    for (var r = 0; r < GRID_SIZE; r++) {
      out.push(new Array(GRID_SIZE));
      for (var c = 0; c < GRID_SIZE; c++) out[r][c] = g[GRID_SIZE - 1 - c][r];
    }
    return out;
  }
  function parseGrid(grid) {
    var values = [];
    for (var i = 0; i < DATA_CELLS; i++) values.push(grid[DATA_CELLS_LIST[i][0]][DATA_CELLS_LIST[i][1]]);
    return cellsToBytes(values);
  }
  function extract(grid, level) {
    var p = RS_PARAMS[level];
    var stream = parseGrid(grid);
    var total = p.n_blocks * p.data;
    var dataStream = stream.subarray(0, total);
    var eccStream = stream.subarray(total);
    var blocks = [];
    for (var b = 0; b < p.n_blocks; b++) {
      var dataB = [];
      for (var i = b; i < dataStream.length; i += p.n_blocks) dataB.push(dataStream[i]);
      var eccB = [];
      for (var i = b; i < eccStream.length; i += p.n_blocks) eccB.push(eccStream[i]);
      var codeword = dataB.concat(eccB);
      var fixed = rsDecode(codeword, p.ecc);
      if (!fixed) throw new Error('too many errors to correct');
      blocks.push(fixed.slice(0, p.data));
    }
    var raw = new Uint8Array(total);
    for (var b = 0; b < p.n_blocks; b++)
      for (var i = 0; i < p.data; i++) raw[b * p.data + i] = blocks[b][i];
    var length = (raw[0] << 8) | raw[1];
    if (length > total - 2) throw new Error('bad length prefix');
    return raw.slice(2, 2 + length);
  }
  function decodePayload(grid) {
    var g = grid;
    for (var rot = 0; rot < 4; rot++) {
      var level = readFormatBlock(g);
      if (level !== null && findersOk(g)) {
        try {
          return { payload: extract(g, level), rotation: rot };
        } catch (e) { /* try next rotation */ }
      }
      g = rotateGrid(g);
    }
    return null;
  }

  // ------------------------------------------------------------- camera math
  function boxDownscale(img, w, h, scale) {
    var ow = Math.max(1, Math.round(w * scale)), oh = Math.max(1, Math.round(h * scale));
    var out = new Uint8Array(ow * oh * 3);
    for (var y = 0; y < oh; y++) {
      var y0 = Math.floor(y / scale), y1 = Math.max(y0 + 1, Math.floor((y + 1) / scale));
      for (var x = 0; x < ow; x++) {
        var x0 = Math.floor(x / scale), x1 = Math.max(x0 + 1, Math.floor((x + 1) / scale));
        var r = 0, g = 0, b = 0, n = 0;
        for (var yy = y0; yy < y1; yy++)
          for (var xx = x0; xx < x1; xx++) {
            var o = (yy * w + xx) * 3;
            r += img[o]; g += img[o + 1]; b += img[o + 2]; n++;
          }
        var o2 = (y * ow + x) * 3;
        out[o2] = Math.round(r / n); out[o2 + 1] = Math.round(g / n); out[o2 + 2] = Math.round(b / n);
      }
    }
    return [out, ow, oh];
  }

  var FINDER_SCALES = [4, 5, 6, 8, 10, 12, 16, 20, 24];

  function finderPatternAt(rot, r, c) {
    // Colour of the cell that appears at finder-local (r, c) when the code
    // is rotated `rot` times 90 degrees clockwise.
    var r0, c0;
    if (rot === 0) { r0 = r; c0 = c; }
    else if (rot === 1) { r0 = 5 - c; c0 = r; }
    else if (rot === 2) { r0 = 5 - r; c0 = 5 - c; }
    else { r0 = c; c0 = 5 - r; }
    return (r0 + c0) % 8;
  }

  // Coarse finder scoring: single-pixel cell centres from the palette index
  // map, with per-scale lookup tables and early exit.  Fast enough to sweep
  // every position at every scale.
  var _COARSE_CACHE = {};
  var _COARSE_NEED = 17; // 0.45 * 36, rounded up

  function centerScore(idx, w, h, x, y, s) {
    var cache = _COARSE_CACHE[s];
    if (!cache) {
      var half = s >> 1;
      var rx = new Int16Array(36), ry = new Int16Array(36);
      var pv = [[], [], [], []];
      for (var r = 0; r < 6; r++)
        for (var c = 0; c < 6; c++) {
          var i = r * 6 + c;
          rx[i] = c * s + half;
          ry[i] = r * s + half;
          for (var rot = 0; rot < 4; rot++) pv[rot][i] = finderPatternAt(rot, r, c);
        }
      cache = { rx: rx, ry: ry, pv: pv };
      _COARSE_CACHE[s] = cache;
    }
    var best = 0;
    for (var rot = 0; rot < 4; rot++) {
      var m = 0;
      for (var i = 0; i < 36; i++) {
        if (m + (36 - i) < _COARSE_NEED) break;
        var px = x + cache.rx[i], py = y + cache.ry[i];
        if (py < 0 || py >= h || px < 0 || px >= w) continue;
        if (idx[py * w + px] === cache.pv[rot][i]) m++;
      }
      if (m > best) best = m;
      if (best >= _COARSE_NEED && m >= 36 - (36 - _COARSE_NEED)) break;
    }
    return best / 36;
  }

  // Accurate finder scoring: each cell is averaged over its central region via
  // integral images, so only a window aligned to the true cell grid scores
  // high (a half-cell offset smears two colours per cell and loses matches).
  function buildIntegrals(img, w, h) {
    var w1 = w + 1;
    var ir = new Int32Array(w1 * h), ig = new Int32Array(w1 * h), ib = new Int32Array(w1 * h);
    for (var y = 0; y < h; y++) {
      var rowR = 0, rowG = 0, rowB = 0;
      var base = y * w1, prev = base - w1;
      for (var x = 0; x < w; x++) {
        var o = (y * w + x) * 3;
        rowR += img[o]; rowG += img[o + 1]; rowB += img[o + 2];
        var idx = base + x + 1;
        ir[idx] = ir[prev + x + 1] + rowR;
        ig[idx] = ig[prev + x + 1] + rowG;
        ib[idx] = ib[prev + x + 1] + rowB;
      }
    }
    return { ir: ir, ig: ig, ib: ib, w1: w1 };
  }

  function regionSum(integral, w1, x0, y0, x1, y1) {
    return integral[y1 * w1 + x1] - integral[y0 * w1 + x1] - integral[y1 * w1 + x0] + integral[y0 * w1 + x0];
  }

  function avgScore(ints, w, h, x, y, s) {
    var inner = Math.max(1, Math.floor(s * 0.2));
    var best = 0, bestRot = 0;
    for (var rot = 0; rot < 4; rot++) {
      var m = 0;
      for (var r = 0; r < 6; r++)
        for (var c = 0; c < 6; c++) {
          var xa = x + c * s + inner, xb = x + (c + 1) * s - inner;
          var ya = y + r * s + inner, yb = y + (r + 1) * s - inner;
          if (xb <= xa || yb <= ya) continue;
          var area = (xb - xa) * (yb - ya);
          var sr = regionSum(ints.ir, ints.w1, xa, ya, xb, yb);
          var sg = regionSum(ints.ig, ints.w1, xa, ya, xb, yb);
          var sb = regionSum(ints.ib, ints.w1, xa, ya, xb, yb);
          var pix = nearestIndex(Math.round(sr / area), Math.round(sg / area), Math.round(sb / area));
          if (pix === finderPatternAt(rot, r, c)) m++;
        }
      if (m > best) { best = m; bestRot = rot; }
    }
    return { score: best / 36, rot: bestRot };
  }

  // Sum of squared RGB distances of each finder cell's average colour from its
  // expected palette colour.  Near zero only at the exact cell-grid alignment,
  // so it pins the finder position to the true grid (sub-pixel sharp).
  // Sub-pixel refine: fit a parabola to the distance function along each axis
  // and nudge the position to its vertex.  Keeps the warp accurate even when
  // the code is small (a few pixels per cell).
  function refineSubpixel(ints, w, h, x, y, s, rot) {
    function distAt(px, py) { return distScore(ints, w, h, px, py, s, rot, 1); }
    var d0 = distAt(x - 1, y), d1 = distAt(x, y), d2 = distAt(x + 1, y);
    var dx = 0, dy = 0;
    var denom = d0 - 2 * d1 + d2;
    if (denom !== 0) dx = (d0 - d2) / (2 * denom);
    if (dx < -1) dx = -1; else if (dx > 1) dx = 1;
    var e0 = distAt(x, y - 1), e1 = d1, e2 = distAt(x, y + 1);
    denom = e0 - 2 * e1 + e2;
    if (denom !== 0) dy = (e0 - e2) / (2 * denom);
    if (dy < -1) dy = -1; else if (dy > 1) dy = 1;
    return { x: x + dx, y: y + dy };
  }

  function distScore(ints, w, h, x, y, s, rot, inner) {
    if (inner === undefined) inner = Math.max(1, Math.floor(s * 0.2));
    var total = 0;
    for (var r = 0; r < 6; r++)
      for (var c = 0; c < 6; c++) {
        var xa = x + c * s + inner, xb = x + (c + 1) * s - inner;
        var ya = y + r * s + inner, yb = y + (r + 1) * s - inner;
        if (xb <= xa || yb <= ya) continue;
        var area = (xb - xa) * (yb - ya);
        var sr = regionSum(ints.ir, ints.w1, xa, ya, xb, yb) / area;
        var sg = regionSum(ints.ig, ints.w1, xa, ya, xb, yb) / area;
        var sb = regionSum(ints.ib, ints.w1, xa, ya, xb, yb) / area;
        var ec = PALETTE[finderPatternAt(rot, r, c)];
        var dr = sr - ec[0], dg = sg - ec[1], db = sb - ec[2];
        total += dr * dr + dg * dg + db * db;
      }
    return total;
  }

  function coarseSweep(idx, w, h, scales) {
    var candidates = [];
    for (var si = 0; si < scales.length; si++) {
      var s = scales[si];
      if (s * 6 >= w || s * 6 >= h) continue;
      var step = Math.max(2, s >> 1);
      for (var y = 0; y <= h - s * 6; y += step) {
        for (var x = 0; x <= w - s * 6; x += step) {
          var score = centerScore(idx, w, h, x, y, s);
          if (score < 0.45) continue;
          candidates.push({ score: score, x: x, y: y, s: s });
        }
      }
    }
    return candidates;
  }

  function detectFinders(idx, ints, w, h) {
    // Tier 1: large scales only (the common case).  Small scales only get
    // swept when the large ones find nothing (far-away codes).
    var result = refineCandidates(idx, ints, w, h, coarseSweep(idx, w, h, [10, 12, 16, 20, 24]));
    if (!result) result = refineCandidates(idx, ints, w, h, coarseSweep(idx, w, h, FINDER_SCALES));
    return result;
  }

  function refineCandidates(idx, ints, w, h, candidates) {
    if (candidates.length === 0) return null;
    candidates.sort(function (a, b) { return b.score - a.score; });
    var kept = [];
    for (var i = 0; i < candidates.length; i++) {
      var c = candidates[i];
      var cx = c.x + 3 * c.s, cy = c.y + 3 * c.s;
      var tooClose = false;
      for (var k = 0; k < kept.length; k++) {
        var kx = kept[k].x + 3 * kept[k].s, ky = kept[k].y + 3 * kept[k].s;
        var minD = Math.max(kept[k].s, c.s) * 6 * 0.5;
        if (Math.hypot(kx - cx, ky - cy) < minD) { tooClose = true; break; }
      }
      if (!tooClose) kept.push(c);
      if (kept.length >= 12) break;
    }
    var groups = [];
    for (var i = 0; i < kept.length; i++) {
      var c2 = kept[i];
      var c2x = c2.x + 3 * c2.s, c2y = c2.y + 3 * c2.s;
      var placed = false;
      for (var g = 0; g < groups.length; g++) {
        var gx = groups[g][0].x + 3 * groups[g][0].s, gy = groups[g][0].y + 3 * groups[g][0].s;
        if (Math.hypot(gx - c2x, gy - c2y) < 0.06 * w) { groups[g].push(c2); placed = true; break; }
      }
      if (!placed) groups.push([c2]);
    }
    groups.sort(function (a, b) {
      var ma = a.reduce(function (m, x) { return Math.max(m, x.score); }, 0);
      var mb = b.reduce(function (m, x) { return Math.max(m, x.score); }, 0);
      return mb - ma;
    });
    if (groups.length < 4) return null;
    var top = [];
    for (var i = 0; i < 4; i++) {
      var best = groups[i][0];
      for (var j = 1; j < groups[i].length; j++) if (groups[i][j].score > best.score) best = groups[i][j];
      top.push(best);
    }
    // Refine each candidate: position then scale via averaged-cell scoring,
    // then a distance-minimising fine step that pins the exact grid alignment.
    var refined = top.map(function (m) {
      var bx = Math.max(0, m.x - m.s), by = Math.max(0, m.y - m.s);
      var ex = Math.min(w - m.s * 6, m.x + m.s), ey = Math.min(h - m.s * 6, m.y + m.s);
      var bestA = avgScore(ints, w, h, m.x, m.y, m.s);
      var bestX = m.x, bestY = m.y;
      for (var yy = by; yy <= ey; yy++) {
        for (var xx = bx; xx <= ex; xx++) {
          var sc = avgScore(ints, w, h, xx, yy, m.s);
          if (sc.score > bestA.score) { bestA = sc; bestX = xx; bestY = yy; }
        }
      }
      var bestS = m.s;
      for (var ss = Math.max(3, m.s - 2); ss <= m.s + 2; ss++) {
        if (ss * 6 >= w || ss * 6 >= h) continue;
        var sc2 = avgScore(ints, w, h, bestX, bestY, ss);
        if (sc2.score > bestA.score) { bestA = sc2; bestS = ss; }
      }
      var fx = bestX, fy = bestY, fs = bestS;
      var fine = distScore(ints, w, h, bestX, bestY, bestS, bestA.rot, 1);
      var rng = Math.max(4, m.s >> 1);
      for (var ss2 = Math.max(3, bestS - 1); ss2 <= bestS + 1; ss2++) {
        if (ss2 * 6 >= w || ss2 * 6 >= h) continue;
        for (var yy2 = Math.max(0, bestY - rng); yy2 <= Math.min(h - ss2 * 6, bestY + rng); yy2++) {
          for (var xx2 = Math.max(0, bestX - rng); xx2 <= Math.min(w - ss2 * 6, bestX + rng); xx2++) {
            var d = distScore(ints, w, h, xx2, yy2, ss2, bestA.rot, 1);
            if (d < fine) { fine = d; fx = xx2; fy = yy2; fs = ss2; }
          }
        }
      }
      var sub = refineSubpixel(ints, w, h, fx, fy, fs, bestA.rot);
      return { score: bestA.score, x: sub.x, y: sub.y, s: fs, rot: bestA.rot };
    });
    var centers = refined.map(function (m) { return [m.x + 2.5 * m.s, m.y + 2.5 * m.s]; });
    var cx = 0, cy = 0;
    for (var i = 0; i < 4; i++) { cx += centers[i][0]; cy += centers[i][1]; }
    cx /= 4; cy /= 4;
    var order = [0, 1, 2, 3].sort(function (i, j) {
      return Math.atan2(centers[i][1] - cy, centers[i][0] - cx) - Math.atan2(centers[j][1] - cy, centers[j][0] - cx);
    });
    return [refined[order[0]], refined[order[1]], refined[order[2]], refined[order[3]]];
  }

  function addOuter(A, row) {
    for (var i = 0; i < 9; i++)
      for (var j = 0; j < 9; j++)
        A[i][j] += row[i] * row[j];
  }

  function normalize(points) {
    var minx = Infinity, maxx = -Infinity, miny = Infinity, maxy = -Infinity;
    for (var i = 0; i < points.length; i++) {
      if (points[i][0] < minx) minx = points[i][0];
      if (points[i][0] > maxx) maxx = points[i][0];
      if (points[i][1] < miny) miny = points[i][1];
      if (points[i][1] > maxy) maxy = points[i][1];
    }
    var cx = (minx + maxx) / 2, cy = (miny + maxy) / 2;
    var s = Math.max(maxx - minx, maxy - miny) / 2;
    if (s < 1e-9) s = 1;
    var pts = points.map(function (p) { return [(p[0] - cx) / s, (p[1] - cy) / s]; });
    return { pts: pts, cx: cx, cy: cy, s: s };
  }

  // Symmetric Jacobi eigensolver; returns eigenvalues and eigenvectors (columns).
  function jacobiSymmetric(A) {
    var n = A.length;
    var V = [];
    for (var i = 0; i < n; i++) {
      V.push([]);
      for (var j = 0; j < n; j++) V[i].push(i === j ? 1 : 0);
    }
    var B = A.map(function (row) { return row.slice(); });
    for (var it = 0; it < 60; it++) {
      var p = 0, q = 1, maxv = 0;
      for (var i = 0; i < n; i++)
        for (var j = i + 1; j < n; j++)
          if (Math.abs(B[i][j]) > maxv) { maxv = Math.abs(B[i][j]); p = i; q = j; }
      if (maxv < 1e-11) break;
      var theta = 0.5 * Math.atan2(2 * B[p][q], B[q][q] - B[p][p]);
      var c = Math.cos(theta), s = Math.sin(theta);
      var app = B[p][p], aqq = B[q][q], apq = B[p][q];
      B[p][p] = c * c * app - 2 * s * c * apq + s * s * aqq;
      B[q][q] = s * s * app + 2 * s * c * apq + c * c * aqq;
      B[p][q] = 0; B[q][p] = 0;
      for (var k = 0; k < n; k++) {
        if (k !== p && k !== q) {
          var akp = B[k][p], akq = B[k][q];
          B[k][p] = c * akp - s * akq; B[p][k] = B[k][p];
          B[k][q] = s * akp + c * akq; B[q][k] = B[k][q];
        }
      }
      for (var k = 0; k < n; k++) {
        var vkp = V[k][p], vkq = V[k][q];
        V[k][p] = c * vkp - s * vkq;
        V[k][q] = s * vkp + c * vkq;
      }
    }
    var d = [];
    for (var i = 0; i < n; i++) d.push(B[i][i]);
    return { values: d, vectors: V };
  }

  function mul3(a, b) {
    var out = [[0, 0, 0], [0, 0, 0], [0, 0, 0]];
    for (var i = 0; i < 3; i++)
      for (var j = 0; j < 3; j++)
        for (var k = 0; k < 3; k++) out[i][j] += a[i][k] * b[k][j];
    return out;
  }

  function toMatrix(h) {
    return [[h[0], h[1], h[2]], [h[3], h[4], h[5]], [h[6], h[7], h[8]]];
  }
  function toFlat(m) {
    return [m[0][0], m[0][1], m[0][2], m[1][0], m[1][1], m[1][2], m[2][0], m[2][1], m[2][2]];
  }

  // Normalized DLT: homography as the null-space of the design matrix, found
  // with a symmetric Jacobi eigensolver (handles both affine and perspective).
  function homographyFrom(src, dst) {
    var ns = normalize(src), nd = normalize(dst);
    var A = [];
    for (var i = 0; i < 9; i++) A.push(new Array(9).fill(0));
    for (var i = 0; i < ns.pts.length; i++) {
      var x = ns.pts[i][0], y = ns.pts[i][1], u = nd.pts[i][0], v = nd.pts[i][1];
      addOuter(A, [x, y, 1, 0, 0, 0, -u * x, -u * y, -u]);
      addOuter(A, [0, 0, 0, x, y, 1, -v * x, -v * y, -v]);
    }
    var eig = jacobiSymmetric(A);
    var k = 0;
    for (var i = 1; i < 9; i++) if (eig.values[i] < eig.values[k]) k = i;
    var hn = [];
    for (var i = 0; i < 9; i++) hn.push(eig.vectors[i][k]);
    var Ts = [[1 / ns.s, 0, -ns.cx / ns.s], [0, 1 / ns.s, -ns.cy / ns.s], [0, 0, 1]];
    var TdInv = [[nd.s, 0, nd.cx], [0, nd.s, nd.cy], [0, 0, 1]];
    return toFlat(mul3(mul3(TdInv, toMatrix(hn)), Ts));
  }

  function invertH(h) {
    var a = h[0], b = h[1], c = h[2], d = h[3], e = h[4], f = h[5], g = h[6], k = h[7], l = h[8];
    var det = a * (e * l - f * k) - b * (d * l - f * g) + c * (d * k - e * g);
    if (Math.abs(det) < 1e-12) return null;
    return [
      (e * l - f * k) / det, (c * k - b * l) / det, (b * f - c * e) / det,
      (f * g - d * l) / det, (a * l - c * g) / det, (c * d - a * f) / det,
      (d * k - e * g) / det, (b * g - a * k) / det, (a * e - b * d) / det,
    ];
  }

  function warpPerspective(img, w, h, H, ow, oh) {
    var inv = invertH(H);
    if (!inv) return null;
    var out = new Uint8Array(ow * oh * 3);
    for (var y = 0; y < oh; y++) {
      for (var x = 0; x < ow; x++) {
        var z = inv[6] * x + inv[7] * y + inv[8];
        var sx = (inv[0] * x + inv[1] * y + inv[2]) / z;
        var sy = (inv[3] * x + inv[4] * y + inv[5]) / z;
        var X = Math.floor(sx), Y = Math.floor(sy);
        var o = (y * ow + x) * 3;
        if (X < 0 || Y < 0 || X + 1 >= w || Y + 1 >= h) { out[o] = 255; out[o + 1] = 255; out[o + 2] = 255; continue; }
        var fx = sx - X, fy = sy - Y;
        var o0 = (Y * w + X) * 3;
        for (var c = 0; c < 3; c++) {
          var p00 = img[o0 + c], p10 = img[o0 + 3 + c], p01 = img[o0 + w * 3 + c], p11 = img[o0 + (w + 1) * 3 + c];
          var v = p00 * (1 - fx) * (1 - fy) + p10 * fx * (1 - fy) + p01 * (1 - fx) * fy + p11 * fx * fy;
          out[o + c] = Math.round(v);
        }
      }
    }
    return out;
  }

  function warpFromFinders(data, w, h, finders) {
    var src = [], dst = [];
    var label = [[0, 0], [26, 0], [26, 26], [0, 26]];
    for (var i = 0; i < 4; i++) {
      var m = finders[i];
      var fx = label[i][0], fy = label[i][1];
      for (var r = 0; r < 6; r++)
        for (var c = 0; c < 6; c++) {
          src.push([m.x + c * m.s, m.y + r * m.s]);
          dst.push([(fx + c) * 8, (fy + r) * 8]);
        }
    }
    var H = homographyFrom(src, dst);
    if (!H) return null;
    return warpPerspective(data, w, h, H, 256, 256);
  }

  function bboxCorners(data, w, h) {
    var rS = 0, gS = 0, bS = 0, n = 0;
    for (var y = 0; y < 8; y++)
      for (var x = 0; x < w; x++) {
        var o = (y * w + x) * 3;
        rS += data[o]; gS += data[o + 1]; bS += data[o + 2]; n++;
      }
    var br = rS / n, bg = gS / n, bb = bS / n;
    var minx = w, miny = h, maxx = -1, maxy = -1, count = 0;
    for (var y = 0; y < h; y++)
      for (var x = 0; x < w; x++) {
        var o = (y * w + x) * 3;
        var dr = data[o] - br, dg = data[o + 1] - bg, db = data[o + 2] - bb;
        if (dr * dr + dg * dg + db * db > 900) {
          count++;
          if (x < minx) minx = x;
          if (x > maxx) maxx = x;
          if (y < miny) miny = y;
          if (y > maxy) maxy = y;
        }
      }
    if (count < 64) return null;
    var bw = maxx - minx, bh = maxy - miny;
    if (bw < 32 || bh < 32 || bh / bw < 0.7 || bh / bw > 1.4) return null;
    if (bw * bh < 0.06 * w * h) return null;
    return [[minx, miny], [maxx, miny], [maxx, maxy], [minx, maxy]];
  }

  function warpBBox(data, w, h, src) {
    var dst = [[0, 0], [256, 0], [256, 256], [0, 256]];
    var H = homographyFrom(src, dst);
    if (!H) return null;
    return warpPerspective(data, w, h, H, 256, 256);
  }

  function sampleGrid(warped) {
    return sampleGridAt(warped, 0, 0);
  }

  function sampleGridAt(warped, dx, dy) {
    var size = GRID_SIZE, dim = 256;
    var grid = [];
    for (var r = 0; r < size; r++) {
      var row = [];
      for (var c = 0; c < size; c++) {
        var x0 = (c * dim / size) | 0, x1 = ((c + 1) * dim / size) | 0;
        var y0 = (r * dim / size) | 0, y1 = ((r + 1) * dim / size) | 0;
        var ix0 = x0 + ((x1 - x0) / 5 | 0) + dx, ix1 = x1 - ((x1 - x0) / 5 | 0) + dx;
        var iy0 = y0 + ((y1 - y0) / 5 | 0) + dy, iy1 = y1 - ((y1 - y0) / 5 | 0) + dy;
        if (ix1 <= ix0) ix1 = ix0 + 1;
        if (iy1 <= iy0) iy1 = iy0 + 1;
        var r2 = 0, g2 = 0, b2 = 0, n2 = 0;
        for (var yy = Math.max(0, iy0); yy < Math.min(dim, iy1); yy++)
          for (var xx = Math.max(0, ix0); xx < Math.min(dim, ix1); xx++) {
            var o = (yy * dim + xx) * 3;
            r2 += warped[o]; g2 += warped[o + 1]; b2 += warped[o + 2]; n2++;
          }
        row.push(n2 ? nearestIndex(Math.round(r2 / n2), Math.round(g2 / n2), Math.round(b2 / n2)) : 7);
      }
      grid.push(row);
    }
    return grid;
  }

  // Counts mismatched cells across the four finder regions for a grid sampled
  // with the given sub-cell shift.  The finder patterns are known exactly, so
  // this scores grid alignment cheaply without touching the data cells.  The
  // warp grid may carry any of the four rotations, so the minimum mismatch
  // over all rotations is used.
  function finderMismatch(warped, dx, dy) {
    var size = GRID_SIZE, dim = 256;
    var label = [[0, 0], [26, 0], [26, 26], [0, 26]];
    var cells = [];
    for (var fi = 0; fi < 4; fi++) {
      var fr = label[fi][0], fc = label[fi][1];
      for (var r = 0; r < 6; r++)
        for (var c = 0; c < 6; c++) {
          var x0 = ((fc + c) * dim / size) | 0, x1 = ((fc + c + 1) * dim / size) | 0;
          var y0 = ((fr + r) * dim / size) | 0, y1 = ((fr + r + 1) * dim / size) | 0;
          var ix0 = x0 + ((x1 - x0) / 5 | 0) + dx, ix1 = x1 - ((x1 - x0) / 5 | 0) + dx;
          var iy0 = y0 + ((y1 - y0) / 5 | 0) + dy, iy1 = y1 - ((y1 - y0) / 5 | 0) + dy;
          if (ix1 <= ix0) ix1 = ix0 + 1;
          if (iy1 <= iy0) iy1 = iy0 + 1;
          var r2 = 0, g2 = 0, b2 = 0, n2 = 0;
          for (var yy = Math.max(0, iy0); yy < Math.min(dim, iy1); yy++)
            for (var xx = Math.max(0, ix0); xx < Math.min(dim, ix1); xx++) {
              var o = (yy * dim + xx) * 3;
              r2 += warped[o]; g2 += warped[o + 1]; b2 += warped[o + 2]; n2++;
            }
          cells.push(n2 ? nearestIndex(Math.round(r2 / n2), Math.round(g2 / n2), Math.round(b2 / n2)) : 7);
        }
    }
    var best = Infinity;
    for (var rot = 0; rot < 4; rot++) {
      var wrong = 0;
      for (var fi2 = 0; fi2 < 4; fi2++)
        for (var r2 = 0; r2 < 6; r2++)
          for (var c2 = 0; c2 < 6; c2++) {
            var idx = fi2 * 36 + r2 * 6 + c2;
            if (cells[idx] !== finderPatternAt(rot, r2, c2)) wrong++;
          }
      if (wrong < best) best = wrong;
    }
    return best;
  }

  // The finder-based homography can be off by a fraction of a cell on small
  // codes; a tiny uniform shift of the sampling grid then smears cell colours.
  // Search a small shift range for the alignment that best matches the known
  // finder patterns, then re-sample the full grid there.
  function alignGrid(warped) {
    var best = finderMismatch(warped, 0, 0);
    var bx = 0, by = 0;
    for (var dy = -4; dy <= 4; dy++)
      for (var dx = -4; dx <= 4; dx++) {
        if (dx === 0 && dy === 0) continue;
        var m = finderMismatch(warped, dx, dy);
        if (m < best) { best = m; bx = dx; by = dy; }
      }
    for (var fx = -1; fx <= 1; fx += 0.5)
      for (var fy = -1; fy <= 1; fy += 0.5) {
        if (fx === 0 && fy === 0) continue;
        var m2 = finderMismatch(warped, bx + fx, by + fy);
        if (m2 < best) { best = m2; bx += fx; by += fy; }
      }
    return { grid: sampleGridAt(warped, bx, by), dx: bx, dy: by };
  }

  // --------------------------------------------------------- public decode
  function decodeFrame(data, width, height) {
    var w = width, h = height, img = data;
    var scale = Math.min(1, 640 / Math.max(w, h));
    if (scale < 1) {
      var resized = boxDownscale(img, w, h, scale);
      img = resized[0]; w = resized[1]; h = resized[2];
    }
    var idx = new Uint8Array(w * h);
    for (var i = 0; i < w * h; i++) {
      var o = i * 3;
      idx[i] = nearestIndex(img[o], img[o + 1], img[o + 2]);
    }
    var ints = buildIntegrals(img, w, h);
    var finders = detectFinders(idx, ints, w, h);
    var warped = null;
    if (finders) warped = warpFromFinders(img, w, h, finders);
    if (!warped) {
      var src = bboxCorners(img, w, h);
      if (src) warped = warpBBox(img, w, h, src);
    }
    if (!warped) return null;
    var grid = null;
    if (finders) {
      grid = alignGrid(warped).grid;
    } else {
      grid = sampleGrid(warped);
    }
    var res = decodePayload(grid);
    if (!res) return null;
    var corners = null;
    if (finders) {
      corners = [
        [finders[0].x, finders[0].y],
        [finders[1].x + 6 * finders[1].s, finders[1].y],
        [finders[2].x + 6 * finders[2].s, finders[2].y + 6 * finders[2].s],
        [finders[3].x, finders[3].y + 6 * finders[3].s],
      ];
    }
    return { payload: res.payload, rotation: res.rotation, corners: corners };
  }

  return {
    nearestIndex: nearestIndex,
    rsEncode: rsEncode,
    decodePayload: decodePayload,
    decodeFrame: decodeFrame,
    // internals, exposed for testing / debugging
    detectFinders: detectFinders,
    buildIntegrals: buildIntegrals,
    avgScore: avgScore,
    distScore: distScore,
    bboxCorners: bboxCorners,
    homographyFrom: homographyFrom,
    warpFromFinders: warpFromFinders,
    warpBBox: warpBBox,
    sampleGrid: sampleGrid,
    sampleGridAt: sampleGridAt,
    finderMismatch: finderMismatch,
    alignGrid: alignGrid,
    boxDownscale: boxDownscale,
  };
})();

if (typeof module !== 'undefined' && module.exports) {
  module.exports = CMBS;
} else {
  var root = typeof self !== 'undefined' ? self : (typeof window !== 'undefined' ? window : null);
  if (root) root.CMBS = CMBS;
}
