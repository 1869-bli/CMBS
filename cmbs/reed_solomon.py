"""Pure-Python Reed-Solomon over GF(256).

Uses the same field and conventions as QR codes:
  * field generator polynomial x^8 + x^4 + x^3 + x^2 + 1 (0x11D)
  * generator roots alpha^0 .. alpha^(n_ecc-1)
  * systematic encoding: codeword = data bytes followed by parity bytes
  * decoder: syndromes -> Berlekamp-Massey -> Chien search -> Forney

Polynomials are stored highest-degree-first: poly[0] is the coefficient of
the highest power of x (the QR / classic "coders" convention).  No third-party
dependencies.
"""


class ReedSolomon:
    def __init__(self):
        self.exp = [0] * 512
        self.log = [0] * 256
        x = 1
        for i in range(255):
            self.exp[i] = x
            self.log[x] = i
            x <<= 1
            if x & 0x100:
                x ^= 0x11D
        for i in range(255, 512):
            self.exp[i] = self.exp[i - 255]

    # -- field arithmetic ------------------------------------------------

    def mul(self, a, b):
        if a == 0 or b == 0:
            return 0
        return self.exp[self.log[a] + self.log[b]]

    def div(self, a, b):
        if b == 0:
            raise ZeroDivisionError("division by zero in GF(256)")
        if a == 0:
            return 0
        return self.exp[(self.log[a] - self.log[b]) % 255]

    def pow(self, base, exp):
        if base == 0:
            return 0
        return self.exp[(self.log[base] * exp) % 255]

    def eval_poly(self, poly, x):
        """Horner's rule.  `poly` is highest-degree-first."""
        y = 0
        for c in poly:
            y = self.mul(y, x) ^ c
        return y

    # -- encoding --------------------------------------------------------

    def _gen_poly(self, n):
        """Generator polynomial: prod_{i=0}^{n-1} (x - alpha^i), monic, highest-first."""
        g = [1]  # polynomial "1", degree 0
        for i in range(n):
            root = self.exp[i]
            newg = [0] * (len(g) + 1)
            for k in range(len(g)):
                newg[k + 1] ^= self.mul(g[k], root)  # root * g
                newg[k] ^= g[k]                      # x * g
            g = newg
        return g

    def encode(self, msg, n_ecc):
        """Return the n_ecc parity bytes for `msg` (data first, then parity)."""
        g = self._gen_poly(n_ecc)
        res = list(msg) + [0] * n_ecc
        for i in range(len(msg)):
            coef = res[i]
            if coef:
                for j in range(1, len(g)):
                    res[i + j] ^= self.mul(g[j], coef)
        return res[len(msg):]

    # -- decoding --------------------------------------------------------

    def _syndromes(self, msg, n_ecc):
        return [self.eval_poly(msg, self.exp[i]) for i in range(n_ecc)]

    def _berlekamp_massey(self, synd):
        n = len(synd)
        lamb = [1]
        b = [1]
        L = 0
        m = 1
        bb = 1
        for i in range(n):
            d = synd[i]
            for j in range(1, L + 1):
                d ^= self.mul(lamb[j], synd[i - j])
            if d == 0:
                m += 1
            elif 2 * L <= i:
                t = lamb[:]
                coef = self.div(d, bb)
                if len(lamb) < m + len(b):
                    lamb += [0] * (m + len(b) - len(lamb))
                for j in range(len(b)):
                    lamb[j + m] ^= self.mul(coef, b[j])
                L = i + 1 - L
                b = t
                bb = d
                m = 1
            else:
                coef = self.div(d, bb)
                if len(lamb) < m + len(b):
                    lamb += [0] * (m + len(b) - len(lamb))
                for j in range(len(b)):
                    lamb[j + m] ^= self.mul(coef, b[j])
                m += 1
        while len(lamb) > 1 and lamb[-1] == 0:
            lamb.pop()
        return lamb

    def _chien(self, lamb, msg_len):
        """Find error positions (array indices) via the error locator's roots."""
        deg = len(lamb) - 1
        lrev = lamb[::-1]  # BM returns the locator lowest-degree-first
        positions = []
        for j in range(msg_len):
            x = self.exp[(255 - j) % 255]  # alpha^-j
            if self.eval_poly(lrev, x) == 0:
                positions.append(msg_len - 1 - j)
        if len(positions) != deg:
            return None
        return positions

    def _forney(self, synd, lamb, positions, n_ecc, codeword):
        omega = [0] * n_ecc
        for i, s in enumerate(synd):
            if s == 0:
                continue
            for j, c in enumerate(lamb):
                if c == 0:
                    continue
                deg = i + j
                if deg < n_ecc:
                    omega[deg] ^= self.mul(s, c)
        omega_rev = omega[::-1]
        msg_len = len(codeword)
        for idx in positions:
            j = msg_len - 1 - idx
            x = self.exp[j]                    # alpha^j
            xinv = self.exp[(255 - j) % 255]   # alpha^-j
            omega_x = self.eval_poly(omega_rev, xinv)
            deriv = 0
            for k in range(1, len(lamb)):
                if k & 1:
                    deriv ^= self.mul(lamb[k], self.pow(xinv, k - 1))
            magnitude = self.mul(x, self.div(omega_x, deriv))
            codeword[idx] ^= magnitude

    def decode(self, codeword, n_ecc):
        """Correct `codeword` (data + parity).  Returns (corrected, ok).

        On success the returned list is the full corrected codeword; the
        caller slices off the parity.  If the errors exceed the correction
        capacity, returns (None, False).
        """
        msg = list(codeword)
        synd = self._syndromes(msg, n_ecc)
        if all(s == 0 for s in synd):
            return msg, True
        lamb = self._berlekamp_massey(synd)
        deg = len(lamb) - 1
        if deg * 2 > n_ecc:
            return None, False
        positions = self._chien(lamb, len(msg))
        if positions is None:
            return None, False
        self._forney(synd, lamb, positions, n_ecc, msg)
        if all(s == 0 for s in self._syndromes(msg, n_ecc)):
            return msg, True
        return None, False


rs = ReedSolomon()
