# Maintainer: maple3142 <kirby741852963@gmail.com>
pkgname=python-lll_cvp
pkgver=0.0.1
pkgrel=1
pkgdesc="A library for solving linear (in)equalities using lattice reduction algorithms"
arch=('x86_64')
url="https://github.com/maple3142/lll_cvp"
license=('GPL')
depends=('sagemath')
makedepends=(python-build python-installer python-wheel)
source=('pyproject.toml' 'lll_cvp.py')
md5sums=('SKIP' 'SKIP')
_name=${pkgname#python-}

build() {
    python -m build --wheel --no-isolation
}

package() {
    python -m installer --destdir="$pkgdir" dist/*.whl
}
