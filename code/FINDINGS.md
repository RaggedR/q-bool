# q-containers: first computational results

Run: `python3 experiments.py` in this directory. All numbers are exact
enumeration over F_2 and F_3; no sampling.

## Definition under test

A **q-container** is `S <| P` with `S` a finite set of shapes and
`P : S -> finVect_Fq`. Extension, Set-valued:

    [| S <| P |] V  =  Sigma_{s in S} Hom_Fq(P s, V)

Morphisms `(S<|P) -> (T<|Q)`: `u : S -> T` plus linear `f_s : Q(u s) -> P s`,
backwards on positions. At q = 1 this degenerates to pointed containers.

## Result 1 — CONFIRMED: the representation theorem holds

The extension functor `Cont_q -> [finVect_q, Set]` is full and faithful.

Verified two ways, neither assuming Yoneda: exhaustive search over all
function families (small cases), and enumeration-by-generic-value with full
naturality verification plus a check that every source element is `C(phi)`
applied to a generic. 12/12 cases agree, q = 2 and 3, up to dim 2.

Closed form, confirmed against the enumeration:

    |Cont_q(C,D)|  =  prod_s sum_t q^{e_t * d_s}        (d = dim P, e = dim Q)

Proof (now clear): `[|C|] = Sigma_s Hom(P s, -)` is a coproduct of
corepresentables, so
`Nat = prod_s [|D|](P s) = prod_s Sigma_t Hom(Q t, P s)` by Yoneda, and
`Sigma_u prod_s Hom(Q(u s), P s)` equals that by **distributivity of prod over
Sigma in Set**. The proof uses precisely the topos-side distributive law that
fails in Vect — which is why Result 2 comes out as it does.

## Result 2 — STRENGTHENED: linearising the shapes destroys everything

Conjecture B said full-and-faithfulness fails for the Vect-valued extension
`(+)_s Hom(P s, V)`. It fails, and maximally:

    |Nat_Vect|  =  prod_s prod_t q^{e_t d_s}  =  q^{(sum_s d_s)(sum_t e_t)}
                =  |Hom( (+)_t Q t , (+)_s P s )|

> **The Vect-valued extension depends only on the total position space
> `(+)_s P s`. All shape information is destroyed.**

Sum vs product: disjoint union keeps the shapes apart, biproduct lets natural
transformations take linear combinations across different target shapes, and
the shape set collapses. Gap is zero exactly when |T| = 1.

## Result 3 — REFUTED: q-containers do not have two derivatives

Predicted: `d_sub` (restrict to a hyperplane) and `d_quot` (quotient by a
line) are distinct, exchanged by the polarity W -> W^perp, merging at q = 1.

**False over Vect_q.** #hyperplanes = [d]_q = #lines, both give position
spaces of dimension d-1, and a q-container remembers only the dimension. So
`d_sub` and `d_quot` are isomorphic for every case tested (q = 2,3; dims 1-3).

The two forms of the Gaussian Pascal recursion are therefore **not** two
derivatives; they are one derivative counted two ways.

Diagnosis: **Vect_q is self-dual** (semisimple, every SES splits, sub and
quotient are interchangeable). The prediction should have been conditional on
that.

### Repaired conjecture (hand-checked, not yet coded)

The two derivatives separate exactly when the base category is **not**
self-dual. Smallest witness, the A_2 quiver `1 -> 2` over F_q, with position
object the projective `P_1 = (F -id-> F)`:

    proper nonzero subobjects of P_1 : only S_2 = (0, F)
    P_1 / S_2                        : S_1 = (F, 0)
    S_1 not iso S_2                  =>  d_sub(P_1) = S_2  /=  S_1 = d_quot(P_1)

So "containers over a quiver" is where the two derivatives are real. Note the
Hall algebra becomes genuinely non-commutative there too, for the same reason
(non-split extensions), so the two repairs are one repair.

## Result 4 — CONFIRMED: the derivative is the Jackson q-derivative

Counting series `f_C(x) = sum_s x^{dim P s}` with `x = q^n`, since
`|[|C|](F_q^n)| = sum_s (q^n)^{dim P s}` (verified by direct enumeration).
Then

    f_{dC}(x)  =  D_q f_C(x),        D_q x^d = [d]_q x^{d-1}

exact for all cases tested, and `[d]_q -> d` as q -> 1, recovering the
ordinary container derivative. `[d]_q` is simultaneously the number of
hyperplanes and the number of lines — which is Result 3 restated.

## Status of the conjectures

| # | claim | prior | verdict |
|---|-------|-------|---------|
| A | representation theorem | 85% | CONFIRMED, and provable by Yoneda |
| B | fails when shapes linearise | 85% | CONFIRMED, and stronger than stated |
| C | two derivatives over F_q | 60% | REFUTED — Vect_q is self-dual |
| C'| two derivatives over a non-self-dual base | new | hand-checked on A_2 |
| D | derivative = Jackson q-derivative | unstated | CONFIRMED |

---

# The flag container and Schubert polynomials

Run: `python3 flag.py`. Exact enumeration of complete flags over F_2 and F_3,
n <= 4; no sampling. Independent divided-difference computation in Z[x].

## The question

`q-container.tex` Targets table has rows for q-Leibniz, q-Pascal, Pieri/LR and
Macdonald, and none for Schubert. The remark at "The Hecke algebra" says
`k[G/B]` is a container and `End_G` of it is `H_n(q)`. Conjecture under test:
the divided difference `d_i f = (f - s_i f)/(x_i - x_{i+1})` that generates
Schubert polynomials is a container operation of the same type as the Jackson
derivative of Theorem D — both being `(1 - group element)/(weight difference)`,
Jackson dividing by a scaling and Demazure by a transposition.

The operator was built from the geometry (the fibration forgetting the i-th
subspace) and its relations read off. Defining it by `d_i S_w = S_{w s_i}`
would have made every check circular.

## Result 5 — CONFIRMED: the Bruhat decomposition, by enumeration

|Fl_n(F_q)| = [n]_q! and |C_w| = q^{l(w)} for every cell, computed from the
rank matrix r_ij = dim(E_i ^ F_j) rather than assumed. n <= 4, q = 2, 3;
5/5 cases. The fibration pi_i has fibre exactly q+1 = |P^1(F_q)| everywhere,
and the fibre-sum operator is constant on cells, so it descends to the cell
module. That it descends is checked, not assumed.

## Result 6 — CONFIRMED: the container gives the braid relations for free

`T_i T_j T_i = T_j T_i T_j` for |i-j| = 1 and `T_i T_j = T_j T_i` for
|i-j| >= 2, on all 24 cells at n = 4. So `S_w` is independent of the reduced
word at container level, which is the half of the conjecture that survives:
no new machinery is needed for well-definedness, and the symmetric group action
on positions is enough.

## Result 7 — REFUTED: the container hands you an idempotent, not a nilpotent

The conjecture said the container operator would satisfy the nilCoxeter
relation `d_i^2 = 0`. It does not, in any case tested.

    fibre sum        S_i^2 = (q+1) S_i        so S_i/(q+1) is IDEMPOTENT
    T_i = S_i - I    T_i^2 = (q-1) T_i + q I  the Hecke quadratic
    nilCoxeter       d_i^2 = 0                NOT satisfied, 5/5 cases

and the failure is structural rather than a bad choice of generator. Writing
`u = a T + b I` and imposing `u^2 = 0` with T, I independent gives
`a^2(q-1) + 2ab = 0` and `a^2 q + b^2 = 0`; the second needs `-q` to be a
square, and solving the pair gives `(b-q)(b+1) = 0` with `b = (q-1)/2`, hence
`q = -1` in both branches. A box search over `|a|,|b| <= 6` finds no nilpotent
combination in any case. **The nilCoxeter algebra is not a specialisation of
the flag container's algebra at any prime power.**

Consequence, and it sharpens the dictionary rather than closing it:

| operator | relation | generates |
|---|---|---|
| S_i/(q+1) | idempotent | Demazure operators, **key polynomials** |
| T_i = S_i - I | Hecke quadratic | H_n(q), Hall-Littlewood, Green |
| d_i | nilpotent | **Schubert polynomials** — not reachable here |

So the container analogue of the Schubert basis is not the Schubert polynomial.
Counting F_q-points sees the Hecke algebra and the Poincare polynomial
`sum_w q^{l(w)} = [n]_q!`; the divided differences live in the Chow ring, which
is a different decategorification of the same cells.

## What the Schubert basis IS, at container level

The flag container is `sum_{w in S_n} X^{l(w)}`: shapes are the permutations,
and the positions over shape w are `F_q^{l(w)}`, the cell. So a Schubert class
is a **shape** of the flag container together with its position space, and the
single-variable shadow of `S_w` is `x^{l(w)}`. The multivariable `S_w(x)` is a
refinement the point count cannot see, because it records the torus weights on
positions and the point count records only their number.

Proposal for the next computation: grade the positions by the T-action and ask
whether the multigraded flag container has `S_w(x_1,...,x_n)` as its shape
polynomial. That is the row the Targets table is missing.

## Status

| # | claim | prior | verdict |
|---|-------|-------|---------|
| 5 | Bruhat decomposition by enumeration | n/a | CONFIRMED, 5/5 |
| 6 | braid relations hold at container level | 70% | CONFIRMED, n <= 4 |
| 7 | d_i^2 = 0 at container level | 60% | **REFUTED** — it is Hecke, not nilCoxeter |
| 8 | multigraded flag container gives S_w(x) | new | open, not attempted |
