import contextlib
import contextvars
import itertools
import logging
import shutil
from re import findall
from subprocess import check_output

from sage.all import (
    ZZ,
    Sequence,
    Zmod,
    block_matrix,
    diagonal_matrix,
    matrix,
    pari,
    vector,
)

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
    return M.dense_matrix().LLL()


def BKZ(M):
    logger.debug(f"BKZ reduction on matrix of size {M.nrows()}x{M.ncols()}")
    return M.dense_matrix().BKZ()


def pari_qflll(M):
    logger.debug(f"pari.qflll reduction on matrix of size {M.nrows()}x{M.ncols()}")
    pm = pari(M.T)
    return (pm * pm.qflll()).sage().T


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


def auto_reduction(M):
    """
    Compute a LLL or flatter reduced basis for the lattice M

    :param M: a matrix
    """
    if not has_flatter:
        return LLL(M)
    nr, nc = M.dimensions()
    mx = abs(max(M.list()))
    if max(nr, nc) < 32 and int(mx).bit_length() < 4096:
        # prefer LLL for small matrices
        return LLL(M)
    # flatter only works in linear indepedent case
    if M.is_square() and M.rank() == nc:
        return flatter(M)
    if nr > nc:
        # definitely not linearly independent
        return LLL(M)
    if M.rank() < nc:
        return LLL(M)
    return flatter(M)


_ctx_reduction = contextvars.ContextVar("reduction", default=auto_reduction)


def reduction(M):
    return _ctx_reduction.get()(M)


def cvp(mat, target):
    return _ctx_cvp.get()(mat, target)


def cvp_ex(mat, target):
    return _ctx_cvp_ex.get()(mat, target)


def babai_cvp(mat, target, reduction=reduction):
    M = reduction(matrix(ZZ, mat))
    G = M.gram_schmidt()[0]
    diff = target
    for i in reversed(range(G.nrows())):
        diff -= M[i] * ((diff * G[i]) / (G[i] * G[i])).round()
    return target - diff


def kannan_cvp(mat, target, reduction=reduction, weight=None):
    """
    Solve closest vector problem using Kannan's algorithm

    :param mat: a matrix
    :param target: a vector
    :returns: a solution as a vector
    """
    return kannan_cvp_ex(mat, target, reduction, weight)[0][0]


_ctx_cvp = contextvars.ContextVar("cvp", default=kannan_cvp)


def kannan_cvp_ex(mat, target, reduction=reduction, weight=None):
    """
    Solve closest vector problem using Kannan's algorithm and return all possible solutions and a reduced basis for enumeration

    :param mat: a matrix
    :param target: a vector
    :returns: a pair of (solutions, basis)
    """
    # still kannan cvp, but return all possible solutions
    # along with a reduced basis (useful for cvp enumeration)
    if weight is None:
        mat = reduction(mat)
        weight = mat[-1].norm().round() + 1
        # another reasonable choice is `weight = max(target)` based on my experience
    L = block_matrix(ZZ, [[mat, 0], [-matrix(target), weight]])
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


_ctx_cvp_ex = contextvars.ContextVar("cvp_ex", default=kannan_cvp_ex)


@contextlib.contextmanager
def lattice_context(*, reduction=auto_reduction, cvp=kannan_cvp, cvp_ex=kannan_cvp_ex):
    """
    Set the reduction function for the context
    :param reduction: a function that takes a lattice basis matrix and returns a reduced basis matrix
    """
    token = _ctx_reduction.set(reduction)
    try:
        yield
    finally:
        _ctx_reduction.reset(token)


def solve_inequality(M, lb, ub, cvp=cvp):
    """
    Find an vector x such that x*M is bounded by lb and ub without checking for correctness
    note that the returned vector is x*M, not x

    :param M: a matrix
    :param lb: a list of lower bounds
    :param ub: a list of upper bounds
    :returns: a solution as a vector
    """
    L, target, Q = build_lattice(M, lb, ub)
    return Q.solve_left(cvp(L * Q, Q * target))


def solve_inequality_ex(M, lb, ub, cvp_ex=cvp_ex):
    """
    Find vectors x such that x*M is bounded by lb and ub without checking for correctness along with a reduced basis for enumeration
    note that the returned vector is x*M, not x

    :param M: a matrix
    :param lb: a list of lower bounds
    :param ub: a list of upper bounds
    :returns: a pair of (solutions, basis)
    """
    L, target, Q = build_lattice(M, lb, ub)
    cvps, basis = cvp_ex(L * Q, Q * target)
    Qi = matrix.diagonal([1 / x for x in Q.diagonal()])
    cvps = (cvps * Qi).change_ring(ZZ)
    basis = (basis * Qi).change_ring(ZZ)
    return cvps, basis


def solve_underconstrained_equations(M, target, lb, ub, cvp=cvp):
    """
    Find an vector x such that x*M=target and x is bounded by lb and ub without checking for correctness

    :param M: a matrix
    :param target: a vector
    :param lb: a list of lower bounds
    :param ub: a list of upper bounds
    :returns: a solution as a vector
    """
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
    """
    Convert polynomials to a matrix and a vector of monomials

    :param polys: a list of polynomials
    :returns: a pair of (matrix, vector) that maxtrix * vector = polys
    """
    # coefficients_monomials is a replacement for coefficient_matrix in sage 10.3
    # and coefficient_matrix is now deprecated
    S = Sequence(polys)
    if hasattr(S, "coefficients_monomials"):
        return S.coefficients_monomials(sparse=False)
    M, monos = S.coefficient_matrix(sparse=False)
    return M, vector(monos)


def solve_multi_modulo_equations(eqs, mods, lb, ub, reduction=reduction, cvp=cvp):
    """
    Solve a linear system of equations modulo different modulus

    :param eqs: a list of equations over ZZ
    :param mods: a list of modulus
    :param lb: a list of lower bounds
    :param ub: a list of upper bounds
    :returns: a solution as a vector
    """
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
    """
    Solve an underconstrained polynomial system over Z/nZ (or ZZ if n is None) where the unknown variables are bounded by some bounds

    :param n: the modulus, can be None if the system is over ZZ
    :param eqs: a list of equations over ZZ
    :param bounds: a dict mapping variable x to an positive integer W, such that |x| < W
    :returns: a generator of solutions, each solution is a pair of (monomials, solution)
    """
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


def interval_mult(it1, it2):
    l1, u1 = it1
    l2, u2 = it2
    l = [l1 * l2, l1 * u2, u1 * l2, u1 * u2]
    l.sort()
    return l[0], l[3]


def interval_pow(it, e):
    l, u = it
    return l**e, u**e


def compute_mono_bounds(mono, bounds):
    # compute upper and lower bounds of a monomial using interval arithmetic
    vs = mono.parent().gens()
    exps = mono.exponents(as_ETuples=False)[0]
    ret = (1, 1)
    for v, e in zip(vs, exps):
        if e > 0:
            ret = interval_mult(ret, interval_pow(bounds[v], e))
    return ret


def solve_underconstrained_equations_general_v2(
    n, eqs, bounds, cvp_ex=cvp_ex, mode="default"
):
    """
    Solve an underconstrained polynomial system over Z/nZ (or ZZ if n is None) where the unknown variables are bounded by some bounds
    This version can possibly achieve better bounds when unbalanced bounds are given
    However, the vectors it returns might not be as good as before due to its use of kannan cvp
    So custom enumeration might be needed

    :param n: the modulus, can be None if the system is over ZZ
    :param eqs: a list of equations over ZZ
    :param bounds: a dict mapping variable x to a tuple of (lower bound, upper bound), can be a single integer W to mean (-W, W)
    :returns: a generator of solutions where each solution is a pair of (monomials, solution) if mode is "default"
              a tuple of (monomials, cvps, basis) if mode is "raw"
    """
    M, monos = polynomials_to_matrix(eqs)
    if n is None:
        L = block_matrix(ZZ, [[M.T, 1]])
    else:
        L = block_matrix(ZZ, [[n, 0], [M.T, 1]])
    bounds = {
        x: tuple(w) if isinstance(w, (list, tuple)) else (-w, w)
        for x, w in bounds.items()
    }
    mono_bounds = [compute_mono_bounds(m, bounds) for m in monos.list()]
    monos_lb, monos_ub = zip(*mono_bounds)
    lb = [0] * len(eqs) + list(monos_lb)
    ub = [0] * len(eqs) + list(monos_ub)
    cvps, basis = solve_inequality_ex(L, lb, ub, cvp_ex=cvp_ex)
    if mode == "raw":
        return monos, cvps, basis
    zeros = [b[len(eqs) :] for b in basis if b[: len(eqs)] == 0]
    for sol in cvps:
        if sol[: len(eqs)] == 0:
            yield vector(monos), sol[len(eqs) :]
            for z in zeros:
                # it is also possible to enumerate ... + ? * z
                yield vector(monos), sol[len(eqs) :] + z
                yield vector(monos), sol[len(eqs) :] - z


def enum_brute(base, basis, lb, ub, *, n=5):
    """
    Enumerate solutions v = base + x @ basis, where lb <= v <= ub using bruteforce
    would try at most (2*n+1)^basis.nrows() solutions

    :param base: the base vector, can be None
    :param basis: the basis matrix, should be a reduced basis
    :param lb: the lower bound vector
    :param ub: the upper bound vector
    :param n: the search limit from [-n, n], default 5
    """
    for muls in itertools.product(range(-n, n + 1), repeat=basis.nrows()):
        v = base + vector(muls) * basis
        if all([l <= x <= u for (l, x, u) in enumerate(zip(lb, v, ub))]):
            yield v


def enum_ilp(base, basis, lb, ub, *, lp=10):
    """
    Enumerate solutions v = base + x @ basis, where lb <= v <= ub using ILP
    powered by CPMPy + ORTools (pip install cpmpy ortools)

    :param base: the base vector, can be None
    :param basis: the basis matrix, should be a reduced basis
    :param lb: the lower bound vector
    :param ub: the upper bound vector
    :param lp: MILP search limit from [-lp, lp], default 10
    """
    import multiprocessing as mp
    from multiprocessing.synchronize import Event as EventClass

    from cpmpy import Model, SolverLookup, cpm_array, intvar
    from cpmpy.solvers.ortools import OrtSolutionPrinter
    from ortools.sat.python import cp_model as ort

    if base is None:
        base = vector([0] * basis.ncols())

    tgt = intvar(-lp, lp, int(basis.nrows())) @ cpm_array(basis)
    model = Model([tgt >= cpm_array(lb - base), tgt <= cpm_array(ub - base)])

    solver = SolverLookup.get("ortools", model)

    q = mp.Queue()
    search = mp.Event()

    def do_solve(model: Model, qu: mp.Queue, search: EventClass):
        def display():
            qu.put(tgt.value().tolist())
            # do not search until it is instructed to do so
            search.wait()
            search.clear()

        solver = SolverLookup.get("ortools", model)
        cb = OrtSolutionPrinter(solver, display=display)
        slvr = ort.CpSolver()
        slvr.parameters.enumerate_all_solutions = True
        slvr.Solve(solver.ort_model, cb)

    proc = mp.Process(target=do_solve, args=(model, q, search))
    proc.start()

    try:
        while True:
            yield base + vector(q.get())
            search.set()
    finally:
        proc.kill()
        proc.join()


def find_ortho(mod=None, *vecs, reduction=reduction):
    """
    Find an short orthogonal basis to the given vectors

    :param mod: the modulus, can be None if the system is over ZZ
    :param vecs: a list of vectors
    :returns: a matrix, each row is orthogonal to the input vectors
    """
    assert len(set(len(v) for v in vecs)) == 1, "vectors have different lengths"
    base = [[matrix(vecs).T, matrix.identity(len(vecs[0]))]]
    if mod is not None:
        base += [[ZZ(mod), 0]]
    L = block_matrix(ZZ, base)
    nv = len(vecs)
    L[:, :nv] *= mod if mod is not None else max([max(v) for v in vecs]) * 2**10
    L = reduction(L)
    ret = []
    for row in L:
        if row[:nv] == 0:
            ret.append(row[nv:])
    return matrix(ret)


def reduce_mod_p(M, p, reduction=reduction, is_field=False):
    """
    Find a short basis for the lattice M modulo p
    p should be a prime number as it need to compute echelon_form

    :param M: a matrix, should be able to change ring to Zmod(p)
    :param p: a prime number
    :returns: a matrix, a short basis for the lattice M modulo p
    """
    nr, nc = M.dimensions()
    if nc < nr:
        raise ValueError("number of columns most not be less than number of rows")
    Me = M.change_ring(Zmod(p, is_field=is_field)).echelon_form()
    L = Me.change_ring(ZZ).stack(
        matrix.zero(nc - nr, nr).augment(matrix.identity(nc - nr) * p)
    )
    return reduction(L)


__all__ = [
    "build_lattice",
    "LLL",
    "BKZ",
    "flatter",
    "auto_reduction",
    "lattice_context",
    "babai_cvp",
    "kannan_cvp",
    "kannan_cvp_ex",
    "solve_inequality",
    "solve_inequality_ex",
    "solve_underconstrained_equations",
    "solve_multi_modulo_equations",
    "polynomials_to_matrix",
    "solve_underconstrained_equations_general",
    "solve_underconstrained_equations_general_v2",
    "enum_brute",
    "enum_ilp",
    "find_ortho",
    "reduce_mod_p",
]
