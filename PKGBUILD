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
source=("lll_cvp::git+https://github.com/maple3142/lll_cvp")
md5sums=('SKIP')
_name=${pkgname#python-}

build() {
    cd $_name
    python -m build --wheel --no-isolation
}

package() {
    cd $_name
    python -m installer --destdir="$pkgdir" dist/*.whl
}
