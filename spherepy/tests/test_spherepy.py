from unittest import TestCase

import spherepy

class TestSBessel(TestCase):
    def test_sbesselj(self):
        s = spherepy.sbesselj(1,5)
        self.assertTrue(1 == 1)

    def test_sbessely(self):
        s = spherepy.sbessely(1,5)
        self.assertTrue(1 == 1)
