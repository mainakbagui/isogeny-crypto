"""test_csidh.py -- checks for csidh.py.  Run: python3 test_csidh.py  (add --512 for CSIDH-512).

Toy prime p = 419 = 4*3*5*7 - 1 (small enough for brute force):
  (i)   every curve reached has #E_A(F_p) = p + 1, counted point by point;
  (ii)  each isogeny maps points of E_A(F_p) to points of E_A'(F_p) (x-only check);
  (iii) the action commutes: [a][b]E_0 = [b][a]E_0, and [-a][a]E_0 = E_0;
  (iv)  the orbit of E_0 has h(-4p) elements, the class number of Z[sqrt(-p)]
        (free and transitive action, CSIDH Theorem 7), computed from reduced forms;
  (v)   validate() accepts exactly the supersingular A in F_p.
"""
import itertools, math, random, sys
from csidh import TOY, CSIDH512, legendre

p = TOY.p


def count_points(A):                                        # #E_A(F_p), including O
    return 1 + sum(1 + (0 if r == 0 else (1 if legendre(r, p) == 1 else -1))
                   for r in ((x ** 3 + A * x * x + x) % p for x in range(p)))


def class_number(D):                                        # reduced primitive forms of discriminant D < 0
    h = 0
    for a in range(1, math.isqrt(-D // 3) + 1):
        for b in range(-a + 1, a + 1):
            if (b * b - D) % (4 * a):
                continue
            c = (b * b - D) // (4 * a)
            if c < a or (c == a and b < 0) or math.gcd(math.gcd(a, b), c) != 1:
                continue
            h += 1
    return h


def test_isogeny_maps_points():
    for A in (0, TOY.act(0, [1, 0, 0]), TOY.act(0, [0, -1, 2])):
        for l in TOY.ells:
            K = next(K for x in range(1, p)
                     for K in [TOY.xMUL((p + 1) // l, (x, 1), A)] if K[1] != 0)
            for x in range(1, p):
                rhs = (x ** 3 + A * x * x + x) % p
                if rhs == 0:
                    continue
                B, (X, Z) = TOY.isogeny(A, K, l, (x, 1))
                if Z == 0:
                    continue                                 # x in the kernel
                u = X * pow(Z, -1, p) % p
                rhs2 = (u ** 3 + B * u * u + u) % p
                assert rhs2 == 0 or legendre(rhs2, p) == legendre(rhs, p)
    print("isogenies map E(F_p) to E'(F_p) and the twist to the twist: ok")


def test_toy():
    orbit = {0}
    for e in itertools.product(range(-3, 4), repeat=3):
        orbit.add(TOY.act(0, e))
    assert all(count_points(A) == p + 1 for A in orbit)
    print(f"all {len(orbit)} curves reached are supersingular: ok")
    h = class_number(-4 * p)
    assert len(orbit) == h, (len(orbit), h)
    print(f"orbit size = h(-4p) = {h}: ok")
    supersingular = {A for A in range(p) if A not in (2, p - 2) and count_points(A) == p + 1}
    assert all(TOY.validate(A) == (A in supersingular) for A in range(p))
    print(f"validate() agrees with point counting on all A in F_{p}: ok "
          f"({len(supersingular)} supersingular A, orbit is {'all' if supersingular == orbit else 'a subset'} of them)")
    for _ in range(50):
        a = [random.randint(-2, 2) for _ in TOY.ells]
        b = [random.randint(-2, 2) for _ in TOY.ells]
        assert TOY.act(TOY.act(0, a), b) == TOY.act(TOY.act(0, b), a)
        assert TOY.act(TOY.act(0, a), [-x for x in a]) == 0
    print("commutativity and inverses: ok")


def test_512():
    C = CSIDH512
    a, A = C.keygen()
    b, B = C.keygen()
    assert C.validate(A) and C.validate(B) and not C.validate(5)
    assert C.shared(a, B) == C.shared(b, A)
    print("CSIDH-512 key exchange agrees, public keys validate, A = 5 rejected: ok")


if __name__ == "__main__":
    test_isogeny_maps_points()
    test_toy()
    if "--512" in sys.argv:
        test_512()
