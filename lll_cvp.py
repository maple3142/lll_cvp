from sage.all import (
    matrix,
    vector,
    matrix,
    block_matrix,
    Sequence,
    ZZ,
    diagonal_matrix,
)
from subprocess import check_output
from re import findall
import shutil, logging

logger = logging.getLogger(__name__)


def build_lattice(mat, lb, ub):
    n = mat.ncols()  # num equations
    m = mat.nrows()  # num variables
    if n != len(ub) or n != len(lb):
        raise ValueError("Number of equations must match number of bounds")
    if any([l > u for l, u in zip(lb, ub)]):
        raise ValueError("All lower bounds must be less than upper bounds")

    L = matrix(ZZ, mat)
    target = vector([(l + u) // 2 for u, l in zip(ub, lb)])

    bounds = [u - l for u, l in zip(ub, lb)]
    K = max(bounds) or L.det()
    Q = matrix.diagonal([K // x if x != 0 else K * n for x in bounds])
    return L, target, Q


def LLL(M):
    logger.debug(f"LLL reduction on matrix of size {M.nrows()}x{M.ncols()}")
    return M.LLL()


def BKZ(M):
    logger.debug(f"BKZ reduction on matrix of size {M.nrows()}x{M.ncols()}")
    return M.BKZ()


def flatter(M):
    logger.debug(f"flatter reduction on matrix of size {M.nrows()}x{M.ncols()}")
    # compile https://github.com/keeganryan/flatter and put it in $PATH
    z = "[[" + "]\n[".join(" ".join(map(str, row)) for row in M) + "]]"
    ret = check_output(["flatter"], input=z.encode())
    return matrix(M.nrows(), M.ncols(), map(int, findall(b"-?\\d+", ret)))


if shutil.which("flatter"):
    has_flatter = True
else:
    has_flatter = False


def auto_choose_reduction(M):
    if not has_flatter:
        return LLL(M)
    if max(M.dimensions()) < 32:
        # prefer LLL for small matrices
        return LLL(M)
    if M.is_square():
        return flatter(M)
    # flatter also works in linear depedent case
    nr, nc = M.dimensions()
    if nr > nc:
        # definitely not linearly independent
        return LLL(M)
    if M.rank() < nc:
        return LLL(M)
    return flatter(M)


default_reduction = auto_choose_reduction


def set_default_reduction(reduction):
    global default_reduction
    default_reduction = reduction


def reduction(M):
    return default_reduction(M)


def babai_cvp(mat, target, reduction=reduction):
    M = reduction(matrix(ZZ, mat))
    G = M.gram_schmidt()[0]
    diff = target
    for i in reversed(range(G.nrows())):
        diff -= M[i] * ((diff * G[i]) / (G[i] * G[i])).round()
    return target - diff


def kannan_cvp(mat, target, reduction=reduction, weight=None):
    if weight is None:
        weight = max(target)
    L = block_matrix([[mat, 0], [-matrix(target), weight]])
    for row in reduction(L):
        if row[-1] < 0:
            row = -row
        if row[-1] == weight:
            return row[:-1] + target


def kannan_cvp_ex(mat, target, reduction=reduction, weight=None):
    # still kannan cvp, but return all possible solutions
    # along with a reduced basis (useful for cvp enumeration)
    if weight is None:
        weight = max(target)
    L = block_matrix([[mat, 0], [-matrix(target), weight]])
    cvps = []
    basis = []
    for row in reduction(L):
        if row[-1] < 0:
            row = -row
        if row[-1] == weight:
            cvps.append(row[:-1] + target)
        elif row[-1] == 0:
            basis.append(row[:-1])
    return matrix(ZZ, cvps), matrix(ZZ, basis)


def solve_inequality(M, lb, ub, cvp=kannan_cvp):
    # find an vector x such that x*M is bounded by lb and ub
    # not checked for correctness
    # note that the returned value is x*M, not x
    L, target, Q = build_lattice(M, lb, ub)
    return Q.solve_left(cvp(L * Q, Q * target))


def solve_inequality_ex(M, lb, ub, cvp_ex=kannan_cvp_ex):
    # find an vector x such that x*M is bounded by lb and ub
    # not checked for correctness
    # note that the returned value is x*M, not x
    L, target, Q = build_lattice(M, lb, ub)
    cvps, basis = cvp_ex(L * Q, Q * target)
    Qi = matrix.diagonal([1 / x for x in Q.diagonal()])
    cvps = (cvps * Qi).change_ring(ZZ)
    basis = (basis * Qi).change_ring(ZZ)
    return cvps, basis


def solve_underconstrained_equations(M, target, lb, ub, cvp=kannan_cvp):
    # find an vector x such that x*M=target and x is bounded by lb and ub
    # not checked for correctness
    n = M.ncols()  # number of equations
    m = M.nrows()  # number of variables
    if n != len(target):
        raise ValueError("number of equations and target mismatch")
    if n >= m:
        raise ValueError("use gauss elimination instead")
    M = block_matrix([[matrix(ZZ, M), 1], [matrix(target), 0]])
    lb = [0] * n + lb
    ub = [0] * n + ub
    sol = solve_inequality(M, lb, ub, cvp=cvp)
    return sol[-m:]


def polynomials_to_matrix(polys):
    # coefficients_monomials is a replacement for coefficient_matrix in sage 10.3
    # and coefficient_matrix is now deprecated
    S = Sequence(polys)
    if hasattr(S, "coefficients_monomials"):
        return S.coefficients_monomials(sparse=False)
    M, monos = S.coefficient_matrix(sparse=False)
    return M, vector(monos)


def solve_multi_modulo_equations(
    eqs, mods, lb, ub, reduction=reduction, cvp=kannan_cvp
):
    # solve a linear system of equations modulo different modulus
    # eqs: a list of equations over ZZ
    # mods: a list of modulus
    if len(eqs) != len(mods):
        raise ValueError("number of equations and modulus mismatch")
    if len(lb) != len(ub):
        raise ValueError("number of lower bounds and upper bounds mismatch")
    M, v = polynomials_to_matrix(eqs)
    assert v.list()[-1] == 1, "only support equations with constant term"
    A, b = M[:, :-1], -M[:, -1]
    M = A.dense_matrix().T
    nr, nc = M.dimensions()
    L = M.stack(diagonal_matrix(mods))
    L = L.augment(matrix.identity(nr).stack(matrix.zero(nc, nr)))
    lbx = b.list() + lb
    ubx = b.list() + ub
    return solve_inequality(L, lbx, ubx, cvp=cvp)[-len(lb) :]


def solve_underconstrained_equations_general(n, eqs, bounds, reduction=reduction):
    # given a underconstrained list of polynomials over Z/nZ (or ZZ if n is None)
    # where the unknown variables are bounded by some bounds
    # bounds should be a dict mapping variable x to an positive integer W, such that |x| < W
    # non-linear monomials will be linearized
    M, monos = polynomials_to_matrix(eqs)
    if n is None:
        L = block_matrix(ZZ, [[M.T, 1]])
    else:
        L = block_matrix(ZZ, [[n, 0], [M.T, 1]])
    bounds = [1] * len(eqs) + [ZZ(m.subs(bounds)) for m in monos.list()]
    K = max(bounds)
    Q = diagonal_matrix([K // b for b in bounds])
    L *= Q
    L = reduction(L)
    L /= Q
    L = L.change_ring(ZZ)
    for row in L:
        if row[: len(eqs)] == 0:
            sol = row[len(eqs) :]
            yield vector(monos), sol


__all__ = [
    "build_lattice",
    "LLL",
    "BKZ",
    "flatter",
    "auto_choose_reduction",
    "set_default_reduction",
    "babai_cvp",
    "kannan_cvp",
    "kannan_cvp_ex",
    "solve_inequality",
    "solve_underconstrained_equations",
    "solve_multi_modulo_equations",
    "polynomials_to_matrix",
    "solve_underconstrained_equations_general",
]
