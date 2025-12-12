# ------------------ helpers: gcd / lcm ------------------
def igcd(a, b):
    a = abs(a); b = abs(b)
    while b:
        a, b = b, a % b
    return a

def ilcm(a, b):
    return abs(a // igcd(a, b) * b) if a and b else 0

# ------------------ rationals as (num, den) reduced ------------------
def rnorm(n, d):
    if d == 0:
        raise ZeroDivisionError
    if n == 0:
        return (0, 1)
    if d < 0:
        n, d = -n, -d
    g = igcd(n, d)
    return (n // g, d // g)

def radd(x, y):
    return rnorm(x[0]*y[1] + y[0]*x[1], x[1]*y[1])

def rsub(x, y):
    return rnorm(x[0]*y[1] - y[0]*x[1], x[1]*y[1])

def rmul(x, y):
    return rnorm(x[0]*y[0], x[1]*y[1])

def rdiv(x, y):
    return rnorm(x[0]*y[1], x[1]*y[0])

def riszero(x):
    return x[0] == 0

# ------------------ exact integer determinant (Bareiss) ------------------
def det_bareiss(M, n):
    A = [row[:] for row in M]
    sign = 1
    prev = 1
    if n == 1:
        return A[0][0]

    for k in range(n - 1):
        if A[k][k] == 0:
            swap_row = -1
            for r in range(k + 1, n):
                if A[r][k] != 0:
                    swap_row = r
                    break
            if swap_row == -1:
                return 0
            A[k], A[swap_row] = A[swap_row], A[k]
            sign *= -1

        pivot = A[k][k]

        for i in range(k + 1, n):
            for j in range(k + 1, n):
                A[i][j] = (A[i][j] * pivot - A[i][k] * A[k][j]) // prev

        for i in range(k + 1, n):
            A[i][k] = 0

        prev = pivot

    return sign * A[n - 1][n - 1]

# ------------------ RREF over rationals ------------------
def rref(mat):
    A = [[rnorm(v, 1) for v in row] for row in mat]
    rows = len(A)
    cols = len(A[0])
    r = 0
    pivots = []

    for c in range(cols):
        # find pivot row
        piv = -1
        for rr in range(r, rows):
            if not riszero(A[rr][c]):
                piv = rr
                break
        if piv == -1:
            continue

        # swap
        A[r], A[piv] = A[piv], A[r]

        # scale pivot row to make pivot = 1
        pv = A[r][c]
        for j in range(c, cols):
            A[r][j] = rdiv(A[r][j], pv)

        # eliminate other rows
        for rr in range(rows):
            if rr == r:
                continue
            if not riszero(A[rr][c]):
                factor = A[rr][c]
                for j in range(c, cols):
                    A[rr][j] = rsub(A[rr][j], rmul(factor, A[r][j]))

        pivots.append((r, c))
        r += 1
        if r == rows:
            break

    return A, pivots

# ------------------ nullspace basis of B (Bx=0) ------------------
def nullspace_basis(B):
    R, pivots = rref(B)
    n = len(B)
    pivot_cols = [c for (_, c) in pivots]
    free_cols = [c for c in range(n) if c not in pivot_cols]

    if not free_cols:
        return []  # only zero vector

    basis = []
    # For each free variable = 1, others free = 0, solve pivots
    for fcol in free_cols:
        x = [rnorm(0, 1) for _ in range(n)]
        x[fcol] = rnorm(1, 1)
        # back out pivot variables from RREF rows
        for (r, c) in pivots:
            # x_c + sum_j R[r][j] x_j = 0  => x_c = -sum_j R[r][j] x_j
            s = rnorm(0, 1)
            for j in free_cols:
                if not riszero(R[r][j]):
                    s = radd(s, rmul(R[r][j], x[j]))
            x[c] = rnorm(-s[0], s[1])

        # convert rational vector to integer vector
        den_lcm = 1
        for v in x:
            den_lcm = ilcm(den_lcm, v[1])
        ints = [v[0] * (den_lcm // v[1]) for v in x]

        # reduce by gcd
        g = 0
        for t in ints:
            g = igcd(g, t)
        if g != 0:
            ints = [t // g for t in ints]

        basis.append(ints)

    return basis

# ------------------ main ------------------
n = int(input("Matrix dimension (max 4): "))

A = []
print("Enter matrix values:")
for r in range(n):
    row = []
    for c in range(n):
        row.append(int(input("a{}{}: ".format(r + 1, c + 1))))
    A.append(row)

lo = int(input("Search lambda from: "))
hi = int(input("to: "))

print("\nInteger eigenvalues and eigenvectors (basis):")
any_found = False

for lam in range(lo, hi + 1):
    # B = A - lam I
    B = []
    for r in range(n):
        row = []
        for c in range(n):
            row.append(A[r][c] - lam if r == c else A[r][c])
        B.append(row)

    if det_bareiss(B, n) == 0:
        any_found = True
        print("\nlambda =", lam)

        basis = nullspace_basis(B)
        if not basis:
            print("  (No nonzero eigenvector found — unexpected)")
        else:
            for idx, v in enumerate(basis):
                print("  v{} = {}".format(idx + 1, v))

if not any_found:
    print("No integer eigenvalues in range")
