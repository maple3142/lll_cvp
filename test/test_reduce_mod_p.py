from sage.all import *
from lll_cvp import reduce_mod_p
from unittest import TestCase


class TestReduceModp(TestCase):
    def test_handle_rref_left_not_identity_case(self):
        p = 65537
        F = GF(p)
        v0 = vector(F, [21, 0, 33, 0, 0, 0])
        v1 = vector(F, [0, 55, 0, 0, 0, 23])
        v2 = vector(F, [0, 0, 0, 65, 42, 0])
        M = random_matrix(F, 3, 3) * matrix(F, [v0, v1, v2])
        self.assertEqual(M.pivots(), (0, 1, 3))
        R = reduce_mod_p(M, p)
        i3 = 1 / F(3)
        self.assertIn(R[0], (v0 * i3, v0 * -i3))
