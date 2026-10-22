#!/usr/bin/env python3
"""Symbolic verification for q-dirichlet.tex.

Every numbered claim in the paper's "Summary of claims" table that is marked
"verified" is checked here.  Run: python3 q-dirichlet-verify.py
"""
import sympy as sp

q, u = sp.symbols('q u')

def qint(n):  return sum(q**i for i in range(n))
def qfac(n):
    r = sp.Integer(1)
    for i in range(1, n + 1): r *= qint(i)
    return sp.expand(r)
def qbin(n, k): return sp.cancel(qfac(n) / (qfac(k) * qfac(n - k)))
def GL(n):
    r = sp.Integer(1)
    for i in range(n): r *= (q**n - q**i)
    return sp.expand(r)
def C2(n): return n * (n - 1) // 2

ok = lambda label, cond: print(("  [ok]   " if cond else "  [FAIL] ") + label)

print("1. |GL_n| = q^C(n,2) (q-1)^n [n]_q!   and   |GL_n|/[n]_q! = |U_B||T|")
for n in range(1, 9):
    ok(f"n={n}", sp.simplify(GL(n) - q**C2(n) * (q - 1)**n * qfac(n)) == 0)

print("\n2. additive (Levi) section coefficient = q^{k(n-k)} [n,k]_q")
for n in range(1, 7):
    for k in range(0, n + 1):
        lhs = sp.cancel(GL(n) / (GL(k) * GL(n - k)))
        assert sp.simplify(lhs - sp.expand(q**(k * (n - k)) * qbin(n, k))) == 0, (n, k)
ok("all n <= 6", True)

print("\n3. tensor leg: c_{a,b} = |GL_ab|/(|GL_a||GL_b|); (q-1)c is a polynomial")
for a in range(1, 4):
    for b in range(1, 4):
        c = sp.cancel(GL(a * b) / (GL(a) * GL(b)))
        den = sp.fraction(sp.together(sp.cancel((q - 1) * c)))[1]
        ok(f"(a,b)=({a},{b})  orbit count (q-1)c in Z[q]", sp.degree(den, q) == 0)

print("\n4. unit of the tensor convolution is (q-1) delta^1:  c_{1,n} = 1/(q-1)")
for n in range(1, 6):
    ok(f"n={n}", sp.simplify(sp.cancel(GL(n) / (GL(1) * GL(n))) - 1 / (q - 1)) == 0)

print("\n5. no F_1 limit unless ab = a+b  (Theorem 4.3)")
for (a, b) in [(1, 2), (2, 2), (2, 3), (3, 3)]:
    lim = sp.limit(sp.cancel(GL(a * b) / (GL(a) * GL(b))), q, 1)
    cls = sp.factorial(a * b) / (sp.factorial(a) * sp.factorial(b))
    agree = sp.simplify(lim - cls) == 0
    ok(f"(a,b)=({a},{b}) ab-a-b={a*b-a-b:+d}  limit={lim} classical={cls}",
       agree == (a * b == a + b))

print("\n6. D_{a,b}(q) in N[q], D(1) = (ab)!/(a!b!)   (Theorem 4.6)")
def D(a, b):
    n = a * b
    return sp.cancel(q**(C2(n) - C2(a) - C2(b)) * qfac(n) / (qfac(a) * qfac(b)))
for (a, b) in [(1, 1), (1, 3), (2, 2), (2, 3), (2, 4), (3, 3), (2, 5), (3, 4)]:
    d = sp.expand(D(a, b))
    den = sp.fraction(sp.together(d))[1]
    poly = sp.degree(den, q) == 0
    P = sp.Poly(d, q)
    cls = sp.factorial(a * b) // (sp.factorial(a) * sp.factorial(b))
    ok(f"(a,b)=({a},{b}) deg={P.degree()} nonneg={all(c >= 0 for c in P.all_coeffs())} "
       f"D(1)={sum(P.all_coeffs())} vs {cls}",
       poly and all(c >= 0 for c in P.all_coeffs()) and sum(P.all_coeffs()) == cls)

print("\n7. q-Mobius: mu_q(p) = -1/(q^C(p,2) [p]_q!) for p prime; q->1 gives -1/p!")
N = 8
c = {n: sp.cancel(1 / (q**C2(n) * qfac(n))) for n in range(1, N + 1)}
M = {1: sp.cancel(1 / c[1])}
for n in range(2, N + 1):
    M[n] = sp.cancel(-sum(c[n // d] * M[d] for d in sp.divisors(n) if d < n) / c[1])
for p in [2, 3, 5, 7]:
    ok(f"p={p}  mu_q({p}) = {sp.factor(M[p])}   q->1 = {sp.limit(M[p], q, 1)}",
       sp.simplify(M[p] + c[p]) == 0 and sp.limit(M[p], q, 1) == -sp.Rational(1, sp.factorial(p)))

print("\n8. F_q[t]: Euler product over irreducibles equals 1/(1-qu)")
def I(n): return sp.cancel(sum(sp.mobius(d) * q**(n // d) for d in sp.divisors(n)) / n)
N = 8
logE = sum(-I(n) * sp.series(sp.log(1 - u**n), u, 0, N + 1).removeO() for n in range(1, N + 1))
logZ = sp.series(sp.log(1 / (1 - q * u)), u, 0, N + 1).removeO()
diff = sp.expand(sp.expand(logE) - sp.expand(logZ))
diff = sum(sp.cancel(sp.expand(diff).coeff(u, k)) * u**k for k in range(0, N + 1))
ok(f"log Euler - log 1/(1-qu) = {sp.simplify(diff)} up to u^{N}", sp.simplify(diff) == 0)
print("   #irreducibles of degree n over F_2:",
      [int(I(n).subs(q, 2)) for n in range(1, 8)])
