# --------- Exact integer determinant (Bareiss algorithm) ----------
def det_bareiss(M, n):
    A = [row[:] for row in M]  # copy
    sign = 1
    prev = 1

    if n == 1:
        return A[0][0]

    for k in range(n - 1):
        # pivoting if needed
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

        # zero out below/left (not required but keeps numbers tidy)
        for i in range(k + 1, n):
            A[i][k] = 0

        prev = pivot

    return sign * A[n - 1][n - 1]


# --------- Main program ----------
n = int(input("Matrix dimension (max 4): "))

A = []
print("Enter matrix values:")
for r in range(n):
    row = []
    for c in range(n):
        row.append(int(input("a{}{}: ".format(r + 1, c + 1))))
    A.append(row)

# Let user choose search range (important!)
lo = int(input("Search lambda from: "))
hi = int(input("to: "))

print("\nEigenvalues (integers found):")
found = 0

for lam in range(lo, hi + 1):
    # Build B = A - lam*I
    B = []
    for r in range(n):
        row = []
        for c in range(n):
            if r == c:
                row.append(A[r][c] - lam)
            else:
                row.append(A[r][c])
        B.append(row)

    if det_bareiss(B, n) == 0:
        print(lam)
        found += 1

if found == 0:
    print("No integer eigenvalues in range")
