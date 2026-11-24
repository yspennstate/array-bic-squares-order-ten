"""Build the mutually orthogonal array-bic squares MOAS(7n+3, 3) of the paper
"Mutually orthogonal array-bic squares of order ten and of every large order"
and check them against Definitions 2.16 and 2.17 of Adriaensen, Bamberg and Williford (arXiv:2611.17855).

The group is F_7 x H with H a product of fields F_q, each q a prime at least 5, 8 or 9.
Points: 1, 2, 3 are the positions, 4, ..., 7n+3 the elements of the group.

Usage:
  python verify_moas.py                 runs n = 1, 5, 7, 8, 9, 11, 25 (5x5), 35 (5x7), with a negative control each
  python verify_moas.py 5 8             runs the single case H = F_5 x F_8 (n = 40)
  python verify_moas.py --print10       prints the three squares of order 10
"""
import itertools, random, sys

# ---------- checker straight from Defs 2.16 and 2.17 ----------
def check_square(A, s, i):
    """A: dict (x,y)->symbol for x!=y in 1..s.  Returns list of failures."""
    bad = []
    allowed = set(range(0, s + 1)) - {i}
    for x in range(1, s + 1):
        for y in range(1, s + 1):
            if x == y:
                if (x, y) in A: bad.append(('diag', x))
                continue
            if A.get((x, y)) not in allowed: bad.append(('entry', x, y, A.get((x, y))))
    for x in range(1, s + 1):
        if x != i and A[(x, i)] != 0: bad.append(('col_i', x))
    for y in range(1, s + 1):
        if y == i: continue
        col = sorted(A[(x, y)] for x in range(1, s + 1) if x != y)
        if col != sorted(set(range(1, s + 1)) - {i}): bad.append(('col', y))
    for x in range(1, s + 1):
        row = sorted(A[(x, y)] for y in range(1, s + 1) if y != x)
        want = sorted(set(range(1, s + 1)) - {i, x}) + ([] if x == i else [0])
        if row != sorted(want): bad.append(('row', x))
    return bad

def check_orth(Ai, i, Ak, k, s):
    seen = {}
    for x in range(1, s + 1):
        for y in range(1, s + 1):
            if x != y:
                p = (Ai[(x, y)], Ak[(x, y)])
                seen[p] = seen.get(p, 0) + 1
    want = set(itertools.product(set(range(s + 1)) - {i}, set(range(s + 1)) - {k}))
    want -= {(k, i)} | {(c, c) for c in set(range(s + 1)) - {i, k}}
    return set(seen) == want and all(v == 1 for v in seen.values())

def check_moas(squares, s):
    out = []
    for i, A in squares.items():
        b = check_square(A, s, i)
        if b: out.append(('square', i, b[:3]))
    for i, k in itertools.combinations(sorted(squares), 2):
        if not check_orth(squares[i], i, squares[k], k, s): out.append(('orth', i, k))
    return out

# ---------- additive groups G = F7 x H, H a product of fields ----------
def field(q):
    """(add, mul, elements) for F_q, q in {prime, 8, 9}."""
    if q in (8, 9):
        p, m, poly = (2, 3, (1, 1, 0)) if q == 8 else (3, 2, (1, 0))  # x^3 = x + 1 ; x^2 = -1
        els = list(itertools.product(range(p), repeat=m))
        add = lambda a, b: tuple((u + v) % p for u, v in zip(a, b))
        def mul(a, b):
            prod = [0] * (2 * m - 1)
            for e1, u in enumerate(a):
                for e2, v in enumerate(b):
                    prod[e1 + e2] = (prod[e1 + e2] + u * v) % p
            for e in range(2 * m - 2, m - 1, -1):  # reduce x^e = x^(e-m) * x^m
                c = prod[e]; prod[e] = 0
                if q == 8:   # x^3 = 1 + x  (coefficient tuple is (c0, c1, c2))
                    prod[e - 3] = (prod[e - 3] + c) % p; prod[e - 2] = (prod[e - 2] + c) % p
                else:        # x^2 = -1
                    prod[e - 2] = (prod[e - 2] - c) % p
            return tuple(prod[:m])
        return add, mul, els
    els = [(a,) for a in range(q)]
    return (lambda a, b: ((a[0] + b[0]) % q,)), (lambda a, b: ((a[0] * b[0]) % q,)), els

def build(qs):
    """qs: list of field sizes for H (empty -> n = 1).  Returns (s, squares)."""
    fs = [field(q) for q in qs]
    H = list(itertools.product(*[f[2] for f in fs])) if fs else [()]
    G = [(a, h) for a in range(7) for h in H]
    zeroH = tuple(f[2][0] for f in fs)  # first element is the zero tuple
    gadd = lambda g1, g2: ((g1[0] + g2[0]) % 7, tuple(f[0](u, v) for f, u, v in zip(fs, g1[1], g2[1])))
    s = 3 + len(G)
    label = {g: 4 + idx for idx, g in enumerate(G)}
    sq = {1: {}, 2: {}, 3: {}}
    def put(r, c, outs):
        for pos, v in zip((1, 2, 3), outs):
            key = (r, c)
            assert key not in sq[pos], ('clash', key)
            sq[pos][key] = v
    # inner cells: Example 2.18 MOAS(3,3)
    ex = {1: [[None, 3, 2], [0, None, 3], [0, 2, None]],
          2: [[None, 0, 3], [3, None, 1], [1, 0, None]],
          3: [[None, 2, 0], [1, None, 0], [2, 1, None]]}
    for pos in (1, 2, 3):
        for x in range(3):
            for y in range(3):
                if x != y: sq[pos][(x + 1, y + 1)] = ex[pos][x][y]
    # entries: ('I',k) inner label, ('O',u) outer residue in F7, '*' the symbol 0
    R = [ (('I',1), ('O',0), ('O',0), ('O',3), ('O',4)),
          (('O',0), ('I',1), '*', ('O',1), ('O',4)),
          (('O',0), ('O',1), ('I',2), ('O',6), ('O',3)),
          (('O',0), ('O',5), ('I',3), ('O',4), ('O',6)) ]
    def rot(e):
        if e == '*': return e
        if e[0] == 'I': return ('I', e[1] % 3 + 1)
        return ('O', (2 * e[1]) % 7)
    recs = []
    for r in R:
        cur = r
        for _ in range(3):
            recs.append(cur)
            row, col, a, b, c = cur
            cur = (rot(row), rot(col), rot(c), rot(a), rot(b))
        assert cur == r
    def ev(e, g):
        if e == '*': return 0
        if e[0] == 'I': return e[1]
        return label[gadd(((e[1]) % 7, zeroH), g)]
    for g in G:
        for row, col, a, b, c in recs:
            put(ev(row, g), ev(col, g), (ev(a, g), ev(b, g), ev(c, g)))
    # orbits of (0, d, 2d, 3d, 4d)-type records for d with nonzero H-part
    if fs:
        cks = []
        for f in fs:   # five distinct scalars c_0..c_4 in each field
            els = f[2]
            cks.append([els[k] for k in range(5)])
        for d in G:
            if d[1] == zeroH: continue
            for g in G:
                ent = []
                for k in range(5):
                    a = (k * d[0]) % 7
                    hb = tuple(f[1](ck[k], u) for f, ck, u in zip(fs, cks, d[1]))
                    ent.append(label[gadd((a, hb), g)])
                put(ent[0], ent[1], ent[2:])
    return s, sq

def mutate(sq, s):
    rnd = random.Random(7)
    sq2 = {p: dict(A) for p, A in sq.items()}
    keys = [k for k in sq2[2] if sq2[2][k] != 0]
    k1, k2 = rnd.sample(keys, 2)
    sq2[2][k1], sq2[2][k2] = sq2[2][k2], sq2[2][k1]   # a swap keeps the symbol multiset
    return sq2


def run(qs):
    for q in qs:
        if not (q in (8, 9) or (q >= 5 and all(q % d for d in range(2, int(q ** 0.5) + 1)))):
            raise SystemExit('field size %d not supported: use primes at least 5, 8 or 9' % q)
    s, sq = build(qs)
    cells = sum(len(A) for A in sq.values()) // 3
    bad = check_moas(sq, s)
    neg = check_moas(mutate(sq, s), s)
    print('n=%-4s s=%-4d cells=%d (want %d)  failures=%d  negative control failures=%d  %s'
          % ('x'.join(map(str, qs)) or '1', s, cells, s * (s - 1), len(bad), len(neg),
             'OK' if not bad and neg and cells == s * (s - 1) else 'PROBLEM'))
    return not bad and bool(neg)

if __name__ == '__main__':
    args = sys.argv[1:]
    if args == ['--print10']:
        s, sq = build([])
        for p in (1, 2, 3):
            print('A_%d' % p)
            for x in range(1, 11):
                print(' '.join('.' if x == y else str(sq[p][(x, y)]) for y in range(1, 11)))
    elif args:
        sys.exit(0 if run([int(a) for a in args]) else 1)
    else:
        ok = all([run(qs) for qs in ([], [5], [7], [8], [9], [11], [5, 5], [5, 7])])
        print('ALL CASES PASS' if ok else 'SOME CASE FAILED')
        sys.exit(0 if ok else 1)
