"""Worked examples for error-correcting.pdf: composition and recoverability."""
from itertools import combinations, product

def rank(rows, p):
    """rank over F_p of a list of row-vectors"""
    M = [list(r) for r in rows]; R = 0
    if not M: return 0
    ncols = len(M[0])
    for c in range(ncols):
        piv = next((r for r in range(R, len(M)) if M[r][c] % p), None)
        if piv is None: continue
        M[R], M[piv] = M[piv], M[R]
        inv = pow(M[R][c], p-2, p)
        M[R] = [(x*inv) % p for x in M[R]]
        for r in range(len(M)):
            if r != R and M[r][c] % p:
                fac = M[r][c]
                M[r] = [(x - fac*y) % p for x, y in zip(M[r], M[R])]
        R += 1
    return R

def cols(G):           # G given as list of rows
    return [tuple(row[j] for row in G) for j in range(len(G[0]))]

def survives(G, p, k):
    """which subsets of surviving columns still span F_p^k"""
    C = cols(G); n = len(C)
    out = {}
    for size in range(1, n+1):
        good = [S for S in combinations(range(n), size)
                if rank([C[j] for j in S], p) == k]
        out[size] = (len(good), len(list(combinations(range(n), size))))
    return out

print("="*70); print("EXAMPLE 1  parity code vs duplication: same size, different protection")
print("="*70)
G_spc = [[1,0,1],[0,1,1]]     # (x,y) -> (x, y, x+y)
G_dup = [[1,0,1],[0,1,0]]     # (x,y) -> (x, y, x)
for name, G in [("parity  (x,y,x+y)", G_spc), ("duplicate (x,y,x)", G_dup)]:
    C = cols(G)
    print(f"\n{name}   generator columns: {C}")
    for size, (good, tot) in survives(G, 2, 2).items():
        if size in (2,3):
            print(f"   any {size} of 3 symbols arrive: {good}/{tot} patterns recover")
    lost = [j for j in range(3) if rank([C[i] for i in range(3) if i != j], 2) < 2]
    print(f"   single losses that are FATAL: {lost if lost else 'none'}")

print()
print("="*70); print("EXAMPLE 2  Reed-Solomon over F_5: any k of n suffice (MDS)")
print("="*70)
p, k, n = 5, 3, 5
pts = list(range(1, n+1))
G_rs = [[pow(a, i, p) for a in pts] for i in range(k)]     # Vandermonde
print(f"\nk={k} message symbols, n={n} transmitted, evaluations at x={pts}")
for r in G_rs: print("   ", r)
res = survives(G_rs, p, k)
for size in range(1, n+1):
    good, tot = res[size]
    mark = "  <-- any k=3 of 5 recover" if size == k and good == tot else ""
    print(f"   {size} symbols arrive: {good}/{tot} patterns recover{mark}")

print()
print("="*70); print("EXAMPLE 3  a pipeline: rank never increases, so the weakest stage decides")
print("="*70)
def matmul(A, B, p):
    return [[sum(A[i][t]*B[t][j] for t in range(len(B))) % p
             for j in range(len(B[0]))] for i in range(len(A))]
# stage 1: encode 3 -> 5  (Reed-Solomon above, a 3x5 matrix)
# stage 2: channel drops symbols 1 and 4  -> select columns {0,2,3}: 5 -> 3
sel = [[1 if j == s else 0 for s in (0,2,3)] for j in range(5)]
# stage 3: a downstream node that only forwards two of the three it holds
fwd = [[1,0],[0,1],[0,0]]
s1 = G_rs
s12 = matmul(s1, sel, p)
s123 = matmul(s12, fwd, p)
print(f"\nafter encode          : rank {rank([list(r) for r in zip(*s1)], p)}  (need 3)")
print(f"after channel loss    : rank {rank([list(r) for r in zip(*s12)], p)}  (need 3)  -> still recoverable")
print(f"after lossy forwarder : rank {rank([list(r) for r in zip(*s123)], p)}  (need 3)  -> LOST, permanently")
print("\nno downstream stage can raise the rank: rank(fg) <= min(rank f, rank g).")
print("redundancy must be added BEFORE the loss, never after.")
