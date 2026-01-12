from pathlib import Path

from pytest import raises
from sympy import primerange

import twosquares

from twosquares import decompose_number, decompose_prime

def test_compiled_tests():
    """Verify that we are running these tests with a compiled version of twosquares"""
    path = Path(twosquares.__file__)
    assert path.suffix.lower() != '.py'


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

        found_one = False
        for i in primerange(2**31 + 1, 2**31 + 1001):
            if i % 4 != 1:
                continue

            x, y = decompose_prime(i)
            assert i == x**2 + y**2
            found_one = True

        assert found_one

    def test_examples(self):
        """Test some verified examples"""
        examples = {
            19889: (17, 140),
        }

        for example_p, decomposition in examples.items():
            assert decompose_prime(example_p) == decomposition

    def test_d_examples(self):
        """Test some verified examples of higher d values"""
        examples = {
            (41, 2): (3, 4),
            (43, 2): (5, 3),
            (19, 3): (4, 1),
            (37, 3): (5, 2),
            (157, 12): (7, 3),
            (181, 12): (13, 1),
            (2147483929, 10000): (34173, 313),
        }

        for (example_p, example_d), decomposition in examples.items():
            assert decompose_prime(example_p, example_d) == decomposition, \
                f"Not matching for {example_p}, {example_d}"

    def test_two(self):
        """Verify two"""
        assert decompose_prime(2) == (1, 1)
        assert decompose_prime(2, 2) == (0, 1)
        with raises(ValueError, match="Could not decompose"):
            decompose_prime(2, 3)

    def test_three(self):
        """Verify three"""
        with raises(ValueError, match="Could not decompose"):
            decompose_prime(3)
        assert decompose_prime(3, 2) == (1, 1)
        assert decompose_prime(3, 3) == (0, 1)
        with raises(ValueError, match="Could not decompose"):
            decompose_prime(3, 4)



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

    def test_example(self):
        """Test the example from my documentation"""
        answers = decompose_number(19890)

        assert len(answers) == 4
        for x, y in answers:
            assert 19890 == x**2 + y**2

