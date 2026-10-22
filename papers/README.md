# papers/

External sources for the notes in `../`.

## koekoek-swarttouw-askey-scheme.pdf

Roelof Koekoek and René F. Swarttouw, *The Askey-scheme of hypergeometric
orthogonal polynomials and its q-analogue*, Delft University of Technology
Report 98-17 (dated 20 February 1996 in the preprint). arXiv:math/9602214.
127 pages, freely distributed by the authors — this is the open-access
predecessor of Koekoek–Lesky–Swarttouw, *Hypergeometric Orthogonal Polynomials
and Their q-Analogues* (Springer 2010), cited as [KLS] in
`../andrews-askey-roy-research.tex`.

### The charts (PDF page numbers, not printed page numbers)

| PDF page | printed | what |
|---|---|---|
| **24** | 22 | **the Askey scheme** — one page, Wilson/Racah at the top, Hermite at the bottom, levels labelled by ₄F₃ … ₁F₁ |
| **50** | 48 | **the q-scheme, continuous half** — Askey–Wilson at the top |
| **51** | 49 | **the q-scheme, discrete half** — q-Racah at the top |
| 23 | 21 | §0.8, the q-derivative and q-integral; (0.8.5) is the q-product rule |

Chapter 2 is the limit relations between classical families (the *edges* of the
Askey scheme); chapter 4 is the same for the q-families; chapter 5 is how the
classical scheme is recovered from the q-scheme as q → 1.

### Why it is here

`../andrews-askey-roy-research.tex` Remark 1.1: the chart turns the pile into a
poset. Its Phase C asks what the *edges* mean container-theoretically, so
chapters 2 and 4 are the ones that matter to the programme — and §6 F1 predicts
in advance that those edges are analytic limits and will not survive. Chapter 2
is therefore the evidence to read against that prediction.

Note also (0.8.5), the q-product rule `D_q[fg] = f(qx) D_q g + g D_q f` — the
same q-Leibniz rule worked as an instance in `../q-container.tex`.
