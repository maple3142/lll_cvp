# lll_cvp

A library for solving linear (in)equalities using lattice reduction algorithms. See [example.sage](example.sage) for what it can do. You are also welcome to read and understand the code in [lll_cvp.py](lll_cvp.py).

> This library was previously maintained [here](https://gist.github.com/maple3142/5a88040d4d3cb09c4505991cf0f1fe98).

## Installation

```bash
git clone https://github.com/maple3142/lll_cvp
cd lll_cvp
# ensure you are in an Python environment with sage installed
# venv, conda, etc.
pip install -e .
```

or you can just download [lll_cvp.py](lll_cvp.py) and put it in your project.

### Flatter

It is recommended to have [flatter](https://github.com/keeganryan/flatter) (a faster LLL) installed in `PATH` and this library will automatically use it if available.

## About this library

This is heavily inspired by [Inequality Solving with CVP](https://github.com/rkm0959/Inequality_Solving_with_CVP) at start, but also incorporated a lot of tricks I learned from solving CTF challenges.

This library only focus on solving constrained and underdetermined linear equations, so it doesn't have anything related to Coppersmith's method. Please refer to the following libraries if you need them:

* https://github.com/defund/coppersmith
* https://github.com/josephsurin/lattice-based-cryptanalysis
* https://github.com/kionactf/coppersmith

> Technically, `solve_underconstrained_equations_general` can handle non-linear equations by linearizing them, but Coppermith's method works much better in most cases.

There are also some libraries that have some overlap with this library:

* https://github.com/TheBlupper/linineq
* https://github.com/josephsurin/lattice-based-cryptanalysis
* https://github.com/jvdsn/crypto-attacks
* https://github.com/nneonneo/pwn-stuff/blob/master/math/solvelinmod.py
