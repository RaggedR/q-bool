"""The flag container, its divided differences, and which Hecke algebra they generate.

Setting: q-container.tex, remark at "The Hecke algebra" — k[G/B] is a free module
on a hom-set, a container in the k[-] sense, and End_G of it is H_n(q). The
Targets table has rows for q-Leibniz, q-Pascal, Pieri/LR and Macdonald, and no
row for Schubert. This file is an attempt to earn one.

The claim under test. Schubert polynomials are generated from the top class by
divided differences

    d_i f = (f - s_i f) / (x_i - x_{i+1}),

and the container derivative of Theorem D is the Jackson operator

    D_q f = (f - sigma_q f) / ((1-q) x),   sigma_q f(x) = f(qx).

Both have the shape (1 - group element) / (difference of weights): Jackson divides
by a scaling, Demazure by a transposition. Containers already carry the
transposition, because relabelling positions is what a container's symmetry IS.
So the conjecture is that d_i needs no new machinery.

What is NOT assumed. The operator is built from the geometry — the fibration
that forgets the i-th subspace of a flag — and its relations are then READ OFF.
Defining d_i by the rule d_i S_w = S_{w s_i} would make every check below
circular; the point is to see which algebra the container actually produces.

Prime fields only (inverses by Fermat), exact enumeration throughout, as in
qcontainer.py and erasure.py.
"""

from itertools import product as iproduct
from functools import lru_cache

# -- F_q^n, subspaces as frozensets of vectors ------------------------------
# A subspace is stored as the frozenset of its elements. That makes intersection
# a set operation and dimension a logarithm, which is all the linear algebra the
# rank matrix below needs, and it cannot disagree with itself about canonical
# forms the way an RREF representation can.

def vectors(n, q):
    return list(iproduct(range(q), repeat=n))

def span(gens, n, q):
    S = {tuple([0] * n)}
    for g in gens:
        S = {tuple((a + c * b) % q for a, b in zip(v, g)) for v in S for c in range(q)}
    return frozenset(S)

def dim(S, q):
    d, m = 0, 1
    while m < len(S):
        m *= q; d += 1
    return d

def extensions(V, n, q, whole):
    """Every W with V < W and dim W = dim V + 1. Count is [n-d]_q."""
    out = set()
    for v in whole:
        if v in V: continue
        out.add(span(list(V) + [v], n, q))
    return sorted(out, key=lambda S: sorted(S))

# -- Complete flags ---------------------------------------------------------

def flags(n, q):
    """Every complete flag, as the chain (V_0, V_1, ..., V_n)."""
    whole = vectors(n, q)
    zero = frozenset({tuple([0] * n)})
    full = frozenset(whole)
    chains = [(zero,)]
    for _ in range(n - 1):
        nxt = []
        for ch in chains:
            for W in extensions(ch[-1], n, q, whole):
                nxt.append(ch + (W,))
        chains = nxt
    return [ch + (full,) for ch in chains]

def standard_flag(n, q):
    whole = vectors(n, q)
    basis = [tuple(1 if j == i else 0 for j in range(n)) for i in range(n)]
    return tuple(span(basis[:i], n, q) for i in range(n + 1))

# -- Relative position: the permutation, from the rank matrix ---------------
# r_ij = dim(E_i ^ F_j), and w(j) = i exactly where the second difference is 1.
# This is the Bruhat decomposition, computed rather than assumed.

def relative_position(E, F, n, q):
    r = [[dim(E[i] & F[j], q) for j in range(n + 1)] for i in range(n + 1)]
    w = [None] * n
    for j in range(1, n + 1):
        for i in range(1, n + 1):
            if r[i][j] - r[i - 1][j] - r[i][j - 1] + r[i - 1][j - 1] == 1:
                w[j - 1] = i
    return tuple(w)

def inversions(w):
    return sum(1 for a in range(len(w)) for b in range(a + 1, len(w)) if w[a] > w[b])

def act_right(w, i):
    """w s_i : swap positions i, i+1 (1-indexed i)."""
    v = list(w); v[i - 1], v[i] = v[i], v[i - 1]; return tuple(v)

# -- The fibration that forgets the i-th subspace ---------------------------

def forget(F, i):
    return F[:i] + F[i + 1:]

# -- Polynomials, for the independent Schubert computation ------------------
# dict from exponent tuple to integer coefficient.

def pmul(f, g):
    h = {}
    for a, ca in f.items():
        for b, cb in g.items():
            k = tuple(x + y for x, y in zip(a, b))
            h[k] = h.get(k, 0) + ca * cb
    return {k: v for k, v in h.items() if v}

def psub(f, g):
    h = dict(f)
    for k, v in g.items():
        h[k] = h.get(k, 0) - v
    return {k: v for k, v in h.items() if v}

def pswap(f, i):
    out = {}
    for a, c in f.items():
        b = list(a); b[i - 1], b[i] = b[i], b[i - 1]
        out[tuple(b)] = out.get(tuple(b), 0) + c
    return {k: v for k, v in out.items() if v}

def divide_by_difference(f, i, n):
    """Exact division by (x_i - x_{i+1}). Lex order with x_i highest, so the
       leading term of the divisor is x_i and the algorithm terminates."""
    def key(a):
        return (a[i - 1], a[i]) + tuple(a[j] for j in range(n) if j not in (i - 1, i))
    g = dict(f); Q = {}
    while g:
        a = max(g, key=key); c = g[a]
        if a[i - 1] == 0:
            raise ArithmeticError("not divisible by (x_i - x_{i+1})")
        b = list(a); b[i - 1] -= 1; b = tuple(b)
        Q[b] = Q.get(b, 0) + c
        div = {tuple(1 if j == i - 1 else 0 for j in range(n)): 1,
               tuple(1 if j == i else 0 for j in range(n)): -1}
        g = psub(g, pmul({b: c}, div))
    return Q

def divided_difference(f, i, n):
    return divide_by_difference(psub(f, pswap(f, i)), i, n)

def schubert_polys(n):
    """S_w for every w in S_n, by divided differences from x^delta. Returns the
       dict and the number of distinct reduced words checked for agreement."""
    from itertools import permutations
    w0 = tuple(range(n, 0, -1))
    top = {tuple(n - k for k in range(1, n + 1)): 1}   # x_1^{n-1} x_2^{n-2} ...
    S = {w0: top}
    frontier = [w0]
    while frontier:
        nxt = []
        for w in frontier:
            for i in range(1, n):
                v = act_right(w, i)
                if inversions(v) != inversions(w) - 1: continue
                fv = divided_difference(S[w], i, n)
                if v in S:
                    assert S[v] == fv, f"reduced-word disagreement at {v}"
                else:
                    S[v] = fv; nxt.append(v)
        frontier = nxt
    return S

# -- The flag container, and the operator the geometry gives -----------------

def build(n, q, verbose=False):
    """Enumerate flags, decompose into B-orbits by relative position to the
       standard flag, and check the cell sizes against q^len(w)."""
    E = standard_flag(n, q)
    Fs = flags(n, q)
    cell = {}
    for F in Fs:
        w = relative_position(E, F, n, q)
        cell.setdefault(w, []).append(F)
    return E, Fs, cell

def qfact(n, q):
    out = 1
    for k in range(1, n + 1):
        out *= sum(q ** j for j in range(k))
    return out

def fibre_operator(n, q, cell, Fs, i):
    """T_i, built from the fibration pi_i that forgets the i-th subspace:

           (T_i f)(F) = sum over F' with pi_i(F') = pi_i(F) of f(F').

       This is pi_i^* pi_{i*}, the pushforward-pullback along a fibration whose
       fibre is P^1(F_q). Nothing about Schubert cells went into the definition.

       Returns (matrix, descends) where matrix[v][w] is the coefficient of the
       cell-v indicator in T_i applied to the cell-w indicator, and `descends`
       records whether T_i of a cell indicator is CONSTANT ON CELLS at all --
       which has to be checked, not assumed, or the matrix is meaningless."""
    where = {}
    for w, fs in cell.items():
        for F in fs: where[F] = w
    fibres = {}
    for F in Fs:
        fibres.setdefault(forget(F, i), []).append(F)

    mat, descends = {}, True
    for w in cell:
        vals = {}
        for F in Fs:
            vals[F] = sum(1 for Fp in fibres[forget(F, i)] if where[Fp] == w)
        row = {}
        for v, fs in cell.items():
            got = {vals[F] for F in fs}
            if len(got) != 1:
                descends = False
            c = got.pop() if len(got) == 1 else None
            if c: row[v] = c
        mat[w] = row
    return mat, descends

def compose(A, B, cells):
    """(A o B) as matrices indexed by cells: apply B then A."""
    out = {}
    for w in cells:
        acc = {}
        for v, c in B.get(w, {}).items():
            for u, d in A.get(v, {}).items():
                acc[u] = acc.get(u, 0) + c * d
        out[w] = {k: v for k, v in acc.items() if v}
    return out

def scale(A, c, cells):
    return {w: {v: c * x for v, x in A.get(w, {}).items()} for w in cells}

def add(A, B, cells):
    out = {}
    for w in cells:
        acc = dict(A.get(w, {}))
        for v, x in B.get(w, {}).items():
            acc[v] = acc.get(v, 0) + x
        out[w] = {k: v for k, v in acc.items() if v}
    return out

def ident(cells):
    return {w: {w: 1} for w in cells}

def eq(A, B, cells):
    return all(A.get(w, {}) == B.get(w, {}) for w in cells)

# -- Can a divided difference live in this algebra at all? ------------------
# The fibration gives T_i with T^2 = (q-1)T + q. Schubert polynomials need a
# NILPOTENT generator, u^2 = 0 (the nilCoxeter algebra). Ask whether any affine
# combination u = a T + b I is nilpotent:
#
#   u^2 = a^2((q-1)T + qI) + 2ab T + b^2 I = (a^2(q-1) + 2ab) T + (a^2 q + b^2) I
#
# and if T, I are linearly independent this vanishes only when
#   a^2(q-1) + 2ab = 0  and  a^2 q + b^2 = 0.
# The second forces b^2 = -a^2 q, so -q must be a square. Over Q, and over R,
# that needs a = b = 0 for every q > 0. Solving the pair symbolically instead of
# over the reals gives b = (q-1)/2 together with (b-q)(b+1) = 0, hence q = -1
# in both branches. So the nilCoxeter algebra is not a specialisation of the
# container's algebra at any prime power.

def T_and_I_independent(T, cells):
    """T, I linearly independent: needed before the argument above applies."""
    I = ident(cells)
    for a in (1, -1, 2, -2):
        for b in (1, -1, 2, -2):
            if eq(add(scale(T, a, cells), scale(I, b, cells), cells),
                  {w: {} for w in cells}, cells):
                return False
    return True

def nilpotent_combination_exists(T, cells, bound=6):
    """Search a T + b I for a^2 T^2-combination that vanishes, over a small
       integer box. An empty result is not a proof; the proof is the symbolic
       argument above. This is the sanity check that it did not miss a case."""
    zero = {w: {} for w in cells}
    for a in range(-bound, bound + 1):
        for b in range(-bound, bound + 1):
            if a == 0 and b == 0: continue
            u = add(scale(T, a, cells), scale(ident(cells), b, cells), cells)
            if eq(compose(u, u, cells), zero, cells):
                return (a, b)
    return None

def report(n, q):
    E, Fs, cell = build(n, q)
    cells = list(cell)
    out = {"n": n, "q": q, "flags": len(Fs), "qfact": qfact(n, q),
           "cells": len(cell)}
    out["bruhat"] = all(len(v) == q ** inversions(w) for w, v in cell.items())
    Ts, ok_hecke, ok_desc = {}, True, True
    for i in range(1, n):
        S, desc = fibre_operator(n, q, cell, Fs, i)
        ok_desc = ok_desc and desc
        T = add(S, scale(ident(cells), -1, cells), cells)
        Ts[i] = T
        lhs = compose(T, T, cells)
        rhs = add(scale(T, q - 1, cells), scale(ident(cells), q, cells), cells)
        ok_hecke = ok_hecke and eq(lhs, rhs, cells)
    out["descends"], out["hecke"] = ok_desc, ok_hecke
    out["nilHecke"] = all(
        eq(compose(Ts[i], Ts[i], cells), {w: {} for w in cells}, cells)
        for i in Ts)
    braid = True
    for i in Ts:
        for j in Ts:
            if j == i + 1:
                braid = braid and eq(
                    compose(Ts[i], compose(Ts[j], Ts[i], cells), cells),
                    compose(Ts[j], compose(Ts[i], Ts[j], cells), cells), cells)
            elif j > i + 1:
                braid = braid and eq(compose(Ts[i], Ts[j], cells),
                                     compose(Ts[j], Ts[i], cells), cells)
    out["braid"] = braid
    out["indep"] = all(T_and_I_independent(T, cells) for T in Ts.values())
    out["nilpotent_combo"] = nilpotent_combination_exists(Ts[1], cells)
    return out

if __name__ == "__main__":
    print("The flag container over F_q: which algebra its fibrations generate\n")
    hdr = ("n", "q", "flags", "[n]_q!", "cells", "Bruhat", "descends",
           "Hecke", "nilHecke", "braid", "nilpotent a T + b I")
    print("  " + "".join(str(h).ljust(10) for h in hdr))
    print("  " + "-" * 108)
    for n, q in [(2, 2), (2, 3), (3, 2), (3, 3), (4, 2)]:
        r = report(n, q)
        row = (r["n"], r["q"], r["flags"], r["qfact"], r["cells"],
               r["bruhat"], r["descends"], r["hecke"], r["nilHecke"],
               r["braid"], r["nilpotent_combo"] or "none")
        print("  " + "".join(str(x).ljust(10) for x in row))

    print("\nSchubert polynomials, by divided differences on Z[x], for S_3 and S_4:")
    for n in (3, 4):
        S = schubert_polys(n)
        print(f"  S_{n}: {len(S)} polynomials, all reduced words agreed")

    print("""
Reading. The fibration that forgets the i-th subspace has fibre P^1(F_q) and
gives S_i with S_i^2 = (q+1) S_i, so S_i/(q+1) is an IDEMPOTENT and
T_i = S_i - I satisfies the Hecke quadratic. The braid relations hold, which is
what makes S_w independent of the reduced word. But no affine combination of
T_i and I is nilpotent, so the nilCoxeter algebra whose generators ARE the
divided differences does not sit inside the algebra the container produces.

The container hands you an idempotent, not a nilpotent. Idempotent divided
differences are the Demazure operators, and they generate key polynomials.
Schubert polynomials need the nilpotent one.""")
