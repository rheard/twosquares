from pathlib import Path

from pytest import raises
from sympy import primerange

import twosquares

from twosquares import decompose_number, decompose_prime

def test_compiled_tests():
    """Verify that we are running these tests with a compiled version of twosquares"""
    path = Path(twosquares.__file__)
    assert path.suffix.lower() == '.pyd'


class TestPrimeDecomposition:
    """Tests for decompose_prime"""

    def test_primes_below_1000(self):
        """Verify small primes"""

        for i in primerange(1000):
            if i % 4 != 1:  # Primes that are not 1 mod 4 will produce an error (see source for decompose_prime)
                continue

            x, y = decompose_prime(i)
            assert i == x**2 + y**2

    def test_invalid_prime(self):
        """Verify invalid primes"""
        with raises(ValueError, match="Could not decompose"):
            decompose_prime(11)

    def test_high_range(self):
        """Verify large primes"""

        for i in primerange(2**31 + 1, 2**31 + 1001):
            if i % 4 != 1:
                continue

            x, y = decompose_prime(i)
            assert i == x**2 + y**2


class TestNumberDecomposition:
    """Tests for decompose_number"""

    def test_numbers_below_1000(self):
        """Verify small numbers"""

        for i in range(1000):
            for x, y in decompose_number(i):
                assert i == x**2 + y**2

    def test_outside_range(self):
        """Verify large numbers"""

        for i in range(2**31 + 1, 2**31 + 1001):
            for x, y in decompose_number(i):
                assert i == x**2 + y**2
