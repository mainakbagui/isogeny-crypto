"""csidh.py -- CSIDH non-interactive key exchange in pure Python 3 (no dependencies).

Mainak Bagui. Reference: W. Castryck, T. Lange, C. Martindale, L. Panny, J. Renes,
"CSIDH: An Efficient Post-Quantum Commutative Group Action", ASIACRYPT 2018.

Setting. p = 4 * l_1 * ... * l_n - 1 with l_i distinct odd primes, and
E_A : y^2 = x^3 + A x^2 + x over F_p. Supersingular E_A has #E_A(F_p) = p + 1, so
pi^2 = -p and each l_i splits in Z[sqrt(-p)] as l_i = (l_i, pi - 1)(l_i, pi + 1).
The class group of Z[sqrt(-p)] acts freely and transitively on such curves; the
ideal l_i = (l_i, pi - 1) sends E_A to E_A / <K>, where K is a point of order l_i
in E_A(F_p), and l_i^{-1} uses a point of order l_i with x in F_p but y not in F_p
(a point of the quadratic twist). A secret key is e in [-m, m]^n, representing
l_1^{e_1} ... l_n^{e_n}; the public key is the Montgomery coefficient of [e] E_0.

Everything is x-only: P = (X : Z) in P^1(F_p). Not constant-time (see README).
"""
import secrets


class CSIDH:
    def __init__(self, ells, m):
        self.ells, self.m = list(ells), m
        self.p = 4 * prod(self.ells) - 1
        assert is_probable_prime(self.p), "4 * prod(l_i) - 1 must be prime"

    # -------- x-only Montgomery arithmetic, a24 = (A + 2) / 4 --------
    def xDBL(self, P, a24):
        p, (X, Z) = self.p, P
        t0, t1 = (X + Z) ** 2 % p, (X - Z) ** 2 % p
        t2 = t0 - t1                                         # 4XZ
        return t0 * t1 % p, t2 * (t1 + a24 * t2) % p

    def xADD(self, P, Q, D):                                 # x(P + Q) from x(P), x(Q), x(P - Q)
        p = self.p
        u = (P[0] - P[1]) * (Q[0] + Q[1]) % p
        v = (P[0] + P[1]) * (Q[0] - Q[1]) % p
        return D[1] * (u + v) ** 2 % p, D[0] * (u - v) ** 2 % p

    def xMUL(self, k, P, A):                                 # Montgomery ladder for x([k] P)
        if k == 0 or P[1] == 0:
            return (1, 0)
        if P[0] == 0:                                        # (0, 0) has order 2; xADD degenerates
            return P if k % 2 else (1, 0)
        a24 = (A + 2) * pow(4, -1, self.p) % self.p
        R0, R1 = P, self.xDBL(P, a24)
        for bit in bin(k)[3:]:
            if bit == '1':
                R0, R1 = self.xADD(R1, R0, P), self.xDBL(R1, a24)
            else:
                R0, R1 = self.xDBL(R0, a24), self.xADD(R1, R0, P)
        return R0

    # -------- odd-degree isogeny (Velu, in Montgomery form) --------
    def isogeny(self, A, K, l, Q):
        """Codomain A' of E_A -> E_A / <K> (K of order l) and the image x(phi(Q)).

        With x_i = x([i] K) for 1 <= i <= d = (l - 1)/2,
            phi(x)  = x * prod_i ((x x_i - 1) / (x - x_i))^2,
            A'      = tau * (A - 3 sigma),  tau = prod_i x_i^2,  sigma = 2 sum_i (x_i - 1/x_i).
        """
        p = self.p
        a24 = (A + 2) * pow(4, -1, p) % p
        mults = [K, self.xDBL(K, a24)]
        while len(mults) < (l - 1) // 2:
            mults.append(self.xADD(mults[-1], K, mults[-2]))
        xs = [X * pow(Z, -1, p) % p for X, Z in mults[:(l - 1) // 2]]
        tau, sigma = 1, 0
        X, Z = Q
        Xn, Zn = X, Z
        for xi in xs:
            tau = tau * xi % p
            sigma = (sigma + xi - pow(xi, -1, p)) % p
            Xn = Xn * (X * xi - Z) ** 2 % p
            Zn = Zn * (X - xi * Z) ** 2 % p
        return tau * tau * (A - 6 * sigma) % p, (Xn, Zn)

    # -------- class group action (Algorithm 2 of the CSIDH paper) --------
    def act(self, A, e):
        p, e = self.p, list(e)
        while any(e):
            x = secrets.randbelow(p - 1) + 1
            rhs = (x * x * x + A * x * x + x) % p
            if rhs == 0:
                continue
            s = 1 if legendre(rhs, p) == 1 else -1          # P in E(F_p) or in the twist
            S = [i for i, ei in enumerate(e) if ei != 0 and (ei > 0) == (s > 0)]
            if not S:
                continue
            k = prod(self.ells[i] for i in S)
            Q = self.xMUL((p + 1) // k, (x, 1), A)            # order of Q divides k
            for i in S:
                l = self.ells[i]
                K = self.xMUL(k // l, Q, A)                   # order 1 or l
                if K[1] != 0:
                    A, Q = self.isogeny(A, K, l, Q)
                    e[i] -= s
                k //= l
        return A

    # -------- keys, shared secret, validation --------
    def keygen(self):
        sk = [secrets.randbelow(2 * self.m + 1) - self.m for _ in self.ells]
        return sk, self.act(0, sk)

    def shared(self, sk, pk):
        if not self.validate(pk):
            raise ValueError("public key is not a supersingular Montgomery coefficient")
        return self.act(pk, sk)

    def validate(self, A):
        """Supersingularity test (CSIDH paper, Section 5): find a point whose order d
        divides p + 1 with d > 4 sqrt(p); Hasse then forces #E_A(F_p) = p + 1."""
        p = self.p
        if A % p in (2, p - 2):
            return False                                     # singular curve
        while True:
            P = (secrets.randbelow(p - 1) + 1, 1)
            d = 1
            for l in self.ells:
                Q = self.xMUL((p + 1) // l, P, A)
                if self.xMUL(l, Q, A)[1] != 0:
                    return False                             # [p + 1] P != O
                if Q[1] != 0:
                    d *= l
                    if d * d > 16 * p:
                        return True


# -------- helpers --------
def prod(xs):
    r = 1
    for x in xs:
        r *= x
    return r


def legendre(a, p):
    return pow(a, (p - 1) // 2, p)


def is_probable_prime(n, rounds=40):
    if n < 2:
        return False
    for q in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        if n % q == 0:
            return n == q
    d, s = n - 1, 0
    while d % 2 == 0:
        d, s = d // 2, s + 1
    for _ in range(rounds):
        x = pow(secrets.randbelow(n - 3) + 2, d, n)
        if x in (1, n - 1):
            continue
        for _ in range(s - 1):
            x = x * x % n
            if x == n - 1:
                break
        else:
            return False
    return True


def odd_primes(count):
    ps, n = [], 3
    while len(ps) < count:
        if all(n % q for q in ps if q * q <= n):
            ps.append(n)
        n += 2
    return ps


TOY = CSIDH([3, 5, 7], m=2)                                  # p = 419, the toy example of the paper
CSIDH512 = CSIDH(odd_primes(73) + [587], m=5)                 # 511-bit p, NIST level 1 (classical)


if __name__ == "__main__":
    import sys, time
    C = TOY if "--toy" in sys.argv else CSIDH512
    print(f"p = {C.p}  ({C.p.bit_length()} bits), n = {len(C.ells)}, m = {C.m}")
    t = time.time()
    a, A = C.keygen()
    b, B = C.keygen()
    sa, sb = C.shared(a, B), C.shared(b, A)
    print(f"Alice public key  A = {A:#x}")
    print(f"Bob   public key  B = {B:#x}")
    print(f"shared secret       = {sa:#x}")
    assert sa == sb
    print(f"Alice and Bob agree. ({time.time() - t:.1f} s)")
