"""
q-containers: containers whose positions are finite-dimensional F_q-vector spaces.

    A q-container is  S <| P   with S a finite set of shapes,
    P : S -> finVect_Fq  assigning a position SPACE to each shape.

    Extension (Set-valued):   [| S <| P |] V  =  Sigma_{s in S} Hom_Fq(P s, V)

    Morphism (S<|P) -> (T<|Q):  u : S -> T   together with
                                f_s : Q(u s) -> P s   LINEAR, backwards.

At q = 1 this should degenerate to ordinary (pointed) containers.

We work in the skeleton of finVect_Fq: objects are dimensions 0,1,2,...
and Hom(m,n) is the set of n-by-m matrices over F_q.

Everything here is exact enumeration; no floating point, no sampling.
"""

from itertools import product
from math import prod

# ---------------------------------------------------------------- field ----

class GF:
    """Prime field F_p.  (Prime powers are not needed for these tests.)"""
    def __init__(self, p):
        self.p = p
    def elements(self):
        return range(self.p)

# --------------------------------------------------------------- linalg ----

def homs(F, m, n):
    """All linear maps F^m -> F^n.

    A map is stored as (n, m, rows) with rows an n-tuple of m-tuples.
    The shape is carried explicitly: a 0-by-d matrix and a 0-by-e matrix are
    different maps, and forgetting that breaks composition through the zero
    object."""
    if m == 0 or n == 0:
        return [(n, m, tuple(() for _ in range(n)))]
    return [(n, m, tuple(rows)) for rows in product(
        list(product(F.elements(), repeat=m)), repeat=n)]

def hom_size(F, m, n):
    return F.p ** (m * n)

def compose(F, g, f):
    """g : F^n -> F^p,  f : F^m -> F^n  |->  g.f : F^m -> F^p."""
    gp, gn, grows = g
    fn, fm, frows = f
    assert gn == fn, f"dimension mismatch: {gn} != {fn}"
    if gp == 0 or fm == 0:
        return (gp, fm, tuple(() for _ in range(gp)))
    out = []
    for i in range(gp):
        row = []
        for j in range(fm):
            acc = 0
            for k in range(gn):
                acc += grows[i][k] * frows[k][j]
            row.append(acc % F.p)
        out.append(tuple(row))
    return (gp, fm, tuple(out))

def identity(F, n):
    if n == 0:
        return (0, 0, ())
    return (n, n, tuple(tuple(1 if i == j else 0 for j in range(n))
                        for i in range(n)))

# ----------------------------------------------------------- q-container ----

# A q-container is represented by the tuple of dimensions of its position
# spaces, one per shape:  C = (d_0, d_1, ..., d_{k-1}).

def extension(F, C, n):
    """[| C |](F^n) as an explicit list of (shape, phi) with phi : P s -> F^n."""
    out = []
    for s, d in enumerate(C):
        for phi in homs(F, d, n):
            out.append((s, phi))
    return out

def extension_size(F, C, n):
    return sum(hom_size(F, d, n) for d in C)

def act(F, f, elem):
    """Functorial action of f : F^n -> F^m on (s, phi) in [|C|](F^n)."""
    s, phi = elem
    return (s, compose(F, f, phi))

# ------------------------------------------------- morphisms vs naturals ----

def container_morphisms(F, C, D):
    """Explicit list of container morphisms C -> D.
       A morphism is (u, (f_s)_s) with u : S -> T and f_s : Q(u s) -> P s."""
    S, T = len(C), len(D)
    out = []
    for u in product(range(T), repeat=S):
        per_shape = [homs(F, D[u[s]], C[s]) for s in range(S)]
        for fs in product(*per_shape):
            out.append((u, fs))
    return out

def container_morphism_count(F, C, D):
    """Closed form:  prod_s  sum_t  |Hom(Q t, P s)|."""
    return prod(sum(hom_size(F, e, d) for e in D) for d in C)

def apply_morphism(F, C, D, mor, elem):
    """Extension of a container morphism: (s, phi) |-> (u s, phi . f_s)."""
    u, fs = mor
    s, phi = elem
    return (u[s], compose(F, phi, fs[s]))

def _search_space(F, C, D, maxdim):
    """Size of the exhaustive search: prod over source elements of |target|."""
    total = 1
    for n in range(maxdim + 1):
        total *= extension_size(F, D, n) ** extension_size(F, C, n)
        if total > 10 ** 8:
            return total
    return total

def _is_natural(F, C, D, alpha, dims, src, index, arrows):
    """Check the naturality square for every arrow and every element."""
    for n in dims:
        for m in dims:
            for f in arrows[(n, m)]:
                for i, e in enumerate(src[n]):
                    lhs = act(F, f, alpha[(n, i)])
                    e2 = act(F, f, e)
                    if lhs != alpha[(m, index[m][e2])]:
                        return False
    return True

def _setup(F, C, D, maxdim):
    dims = list(range(maxdim + 1))
    src = {n: extension(F, C, n) for n in dims}
    tgt = {n: extension(F, D, n) for n in dims}
    arrows = {(n, m): homs(F, n, m) for n in dims for m in dims}
    index = {n: {e: i for i, e in enumerate(src[n])} for n in dims}
    return dims, src, tgt, arrows, index

def brute_force_naturals(F, C, D, maxdim):
    """EXHAUSTIVE: try every function family, keep the natural ones.
       Makes no use of Yoneda.  Only viable for tiny cases."""
    dims, src, tgt, arrows, index = _setup(F, C, D, maxdim)
    slots = [(n, i) for n in dims for i in range(len(src[n]))]
    naturals = []
    for combo in product(*[tgt[n] for (n, _) in slots]):
        alpha = {slot: val for slot, val in zip(slots, combo)}
        if _is_natural(F, C, D, alpha, dims, src, index, arrows):
            naturals.append(combo)
    return naturals

def generic_naturals(F, C, D, maxdim):
    """Enumerate candidate natural families by their values on the GENERIC
       elements (s, id_{P s}), extend by alpha(s, phi) := D(phi)(y_s), and then
       VERIFY naturality of each candidate in full.

       Returns (n_candidates, n_verified_natural, determined) where
       `determined` records that every element of the source is C(phi) applied
       to a generic, so a natural family cannot do anything else -- checked,
       not assumed.  Requires maxdim >= max dim P s."""
    assert maxdim >= max(C), "need the generic objects inside the truncation"
    dims, src, tgt, arrows, index = _setup(F, C, D, maxdim)

    # every source element is C(phi) of a generic?  check it.
    determined = True
    for n in dims:
        for (s, phi) in src[n]:
            gen = (s, identity(F, C[s]))
            if act(F, phi, gen) != (s, phi):
                determined = False
    # candidate values live in [|D|](P s)
    per_shape = [extension(F, D, C[s]) for s in range(len(C))]
    n_cand = 0
    n_nat = 0
    for ys in product(*per_shape):
        n_cand += 1
        alpha = {}
        for n in dims:
            for i, (s, phi) in enumerate(src[n]):
                alpha[(n, i)] = act(F, phi, ys[s])
        if _is_natural(F, C, D, alpha, dims, src, index, arrows):
            n_nat += 1
    return n_cand, n_nat, determined

# --------------------------------------------- the linear (Vect) variant ----

def linear_nat_count(F, C, D):
    """If the extension is Vect-valued, [|C|]V = direct-sum_s Hom(P s, V),
       then Nat = prod_s  direct-sum_t Hom(Q t, P s)  (enriched Yoneda),
       whose cardinality is prod_s prod_t |Hom(Q t, P s)|."""
    return prod(prod(hom_size(F, e, d) for e in D) for d in C)

# ------------------------------------------------------------ derivative ----

def gauss_binom(q, n, k):
    if k < 0 or k > n:
        return 0
    num = den = 1
    for i in range(k):
        num *= q ** (n - i) - 1
        den *= q ** (k - i) - 1
    return num // den

def q_int(q, n):
    """[n]_q = 1 + q + ... + q^{n-1} = number of lines = number of hyperplanes."""
    return sum(q ** i for i in range(n))

def d_sub(q, C):
    """Restrict to a hyperplane: shapes = (s, H) with H a hyperplane of P s."""
    out = []
    for d in C:
        out.extend([d - 1] * gauss_binom(q, d, d - 1) if d >= 1 else [])
    return tuple(out)

def d_quot(q, C):
    """Quotient by a line: shapes = (s, l) with l a line of P s."""
    out = []
    for d in C:
        out.extend([d - 1] * gauss_binom(q, d, 1) if d >= 1 else [])
    return tuple(out)

def counting_series(F, C, n):
    """|[|C|](F^n)| = sum_s x^{dim P s}  with  x = q^n."""
    x = F.p ** n
    return sum(x ** d for d in C)
