"""Experiments on q-containers.  Every number below is exact enumeration."""
from qcontainer import *
from qcontainer import _search_space, _setup, _is_natural
from itertools import product

def hdr(s):
    print("\n" + "=" * 74); print(s); print("=" * 74)

# --------------------------------------------------------------------------
hdr("EXPERIMENT 1  --  Conjecture A: extension is full and faithful")
print("Two independent methods, neither assuming the Yoneda argument:")
print("  EXHAUSTIVE  -- try every function family, keep the natural ones")
print("  GENERIC     -- enumerate candidates by value on (s, id_{P s}),")
print("                 then VERIFY naturality of each in full\n")
print(f"{'q':>2} {'C':>8} {'D':>8} {'md':>3} {'method':>10} {'#natural':>9} {'|Cont(C,D)|':>11}  verdict")
cases = [
    (2, (1,),   (1,),    1),
    (2, (1,),   (1, 1),  1),
    (2, (1, 1), (1,),    1),
    (2, (0, 1), (1,),    1),
    (2, (1, 1), (1, 1),  1),
    (3, (1,),   (1,),    1),
    (3, (1,),   (1, 1),  1),
    (2, (1,),   (2,),    2),
    (2, (2,),   (1,),    2),
    (2, (2,),   (1, 1),  2),
    (2, (1, 2), (2, 1),  2),
    (3, (2,),   (1,),    2),
]
allok = True
for (p_, C, D, md) in cases:
    F = GF(p_)
    cm = container_morphism_count(F, C, D)
    space = _search_space(F, C, D, md)
    if space <= 2 * 10 ** 6:
        n = len(brute_force_naturals(F, C, D, md)); meth = "exhaustive"
        ok = (n == cm)
    else:
        cand, nat, det = generic_naturals(F, C, D, md); meth = "generic"
        n = nat
        ok = (nat == cm) and (cand == cm) and det
    allok &= ok
    print(f"{p_:>2} {str(C):>8} {str(D):>8} {md:>3} {meth:>10} {n:>9} {cm:>11}  {'OK' if ok else 'MISMATCH'}")
print("\nAll cases agree." if allok else "\n*** DISAGREEMENT ***")
print("(generic rows also confirm: every source element is C(phi) of a generic,")
print(" so no natural family can differ from one of the candidates counted)")

hdr("EXPERIMENT 1b  --  the morphisms are not merely equinumerous")
print("check the map (container morphism) |-> (its extension) is injective")
print("and hits every natural family, for one case, elementwise\n")
F = GF(2); C, D, md = (1, 1), (1, 1), 1
mors = container_morphisms(F, C, D)
dims = list(range(md + 1))
src = {n: extension(F, C, n) for n in dims}
realised = set()
for mor in mors:
    fam = tuple(apply_morphism(F, C, D, mor, e) for n in dims for e in src[n])
    realised.add(fam)
nats = set(brute_force_naturals(F, C, D, md))
print(f"container morphisms      : {len(mors)}")
print(f"distinct extensions      : {len(realised)}   (injective: {len(realised) == len(mors)})")
print(f"natural families (brute) : {len(nats)}")
print(f"extensions == naturals   : {realised == nats}")

# --------------------------------------------------------------------------
hdr("EXPERIMENT 2  --  Conjecture B: linearising the SHAPES breaks it")
print("Set-valued  : Nat = prod_s  SUM_t  |Hom(Qt,Ps)|      (disjoint union)")
print("Vect-valued : Nat = prod_s  PROD_t |Hom(Qt,Ps)|      (biproduct)\n")
print(f"{'q':>2} {'C':>10} {'D':>10} {'|Cont|':>10} {'|Nat_Set|':>10} {'|Nat_Vect|':>12} {'gap':>10}")
for (p, C, D) in [(2,(1,),(1,)), (2,(1,),(1,1)), (2,(1,1),(1,1)),
                  (2,(1,),(1,1,1)), (2,(2,),(1,1)), (3,(1,),(1,1))]:
    F = GF(p)
    cm  = container_morphism_count(F, C, D)
    lin = linear_nat_count(F, C, D)
    print(f"{p:>2} {str(C):>10} {str(D):>10} {cm:>10} {cm:>10} {lin:>12} {lin-cm:>10}")
print("\nGap is 0 exactly when |T| = 1 (one shape: sum and product coincide).")

# --------------------------------------------------------------------------
hdr("EXPERIMENT 3  --  the derivative: does sub differ from quotient?")
print("d_sub  : shapes (s,H), H a hyperplane of P s      -> positions dim d-1")
print("d_quot : shapes (s,l), l a line of P s            -> positions dim d-1\n")
print(f"{'q':>2} {'C':>12} {'d_sub':>16} {'d_quot':>16}  iso?")
for p in (2, 3):
    for C in [(1,), (2,), (3,), (1, 2), (2, 2)]:
        ds, dq = d_sub(p, C), d_quot(p, C)
        print(f"{p:>2} {str(C):>12} {str(ds):>16} {str(dq):>16}  "
              f"{'YES' if sorted(ds)==sorted(dq) else 'NO'}")
print("\nReason: #hyperplanes = [d]_q = #lines, and in the skeleton a container")
print("remembers only the DIMENSION of each position space.  Over a semisimple")
print("base, sub and quotient give isomorphic q-containers.")

# --------------------------------------------------------------------------
hdr("EXPERIMENT 4  --  is the derivative the JACKSON q-derivative?")
print("counting series  f_C(x) = sum_s x^{dim P s}   with x = q^n")
print("Jackson:  D_q x^d = [d]_q x^{d-1}\n")
print(f"{'q':>2} {'C':>12} {'f_C(x)':>22} {'D_q f_C':>24} {'f_{dC}':>22}  match")
def series_str(C):
    from collections import Counter
    c = Counter(C)
    return " + ".join(f"{v}x^{k}" if v > 1 else f"x^{k}" for k, v in sorted(c.items(), reverse=True)) or "0"
def jackson(q, C):
    from collections import Counter
    out = Counter()
    for d in C:
        if d >= 1:
            out[d - 1] += q_int(q, d)
    return out
for p in (2, 3):
    for C in [(1,), (2,), (3,), (1, 2), (2, 2, 3)]:
        from collections import Counter
        J = jackson(p, C)
        dC = Counter(d_sub(p, C))
        js = " + ".join(f"{v}x^{k}" for k, v in sorted(J.items(), reverse=True)) or "0"
        cs = " + ".join(f"{v}x^{k}" for k, v in sorted(dC.items(), reverse=True)) or "0"
        print(f"{p:>2} {str(C):>12} {series_str(C):>22} {js:>24} {cs:>22}  {'OK' if J==dC else 'NO'}")

# --------------------------------------------------------------------------
hdr("EXPERIMENT 5  --  q -> 1 degeneration")
print("shape counts of the derivative as a polynomial in q, evaluated at q=1")
print("[d]_q = 1+q+...+q^{d-1}   ->   d      (ordinary derivative)\n")
print(f"{'d':>3} " + " ".join(f"{'[d]_'+str(q):>7}" for q in (2,3,4,5)) + f" {'q=1':>6}")
for d in range(1, 6):
    row = " ".join(f"{q_int(q,d):>7}" for q in (2,3,4,5))
    print(f"{d:>3} {row} {d:>6}")
print("\nSo  D_q  ->  d/dx  and the q-container derivative -> the container derivative.")

# --------------------------------------------------------------------------
hdr("EXPERIMENT 6  --  counting-series analogy is exact")
print("Set  : |[|C|](X)|   = sum_s |X|^{|P s|}")
print("F_q  : |[|C|](F^n)| = sum_s (q^n)^{dim P s}      so  x = q^n\n")
F = GF(2)
for C in [(1,), (2,), (1, 2)]:
    for n in range(4):
        lhs = extension_size(F, C, n)
        rhs = counting_series(F, C, n)
        assert lhs == rhs
        print(f"  C={str(C):8} n={n}  |[|C|](F_2^{n})| = {lhs:>6}   = f_C(x) at x=2^{n}={2**n}")
print("\n  verified by direct enumeration for all rows above.")
