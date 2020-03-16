def solve_diagonal(lower, diag, upper, b):
    """Solve matrix Ax=b when A is diagonal"""
    from scipy.sparse import linalg, diags
    N = len(diag)
    A = diags(
        diagonals=[lower, diag, upper],
        offsets=[-1, 0, 1],
        shape=(N, N),
        format='csr'
    )
    return linalg.spsolve(A, b)