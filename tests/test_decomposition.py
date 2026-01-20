import math
import os
import random

from functools import cache
from pathlib import Path

from pytest import mark, raises
from sympy import primerange

import twosquares

from twosquares import decompose_number, decompose_prime

def test_compiled_tests():
    """Verify that we are running these tests with a compiled version of twosquares"""
    path = Path(twosquares.__file__)
    assert path.suffix.lower() != '.py'


@cache
def brute_force_twosquares(n: int) -> set[tuple[int, int]]:
    """
    Brute-force all canonical (x,y) with 0 <= x <= y and x^2 + y^2 = n.

    Runs in O(isqrt(n)) time by scanning x and checking whether n - x^2 is a square.

    Returns:
        set: The solutions found to n = x**2 + y**2
    """
    out: set[tuple[int, int]] = set()
    lim = math.isqrt(n)
    for x in range(lim + 1):
        r = n - x * x
        y = math.isqrt(r)
        if y * y == r and x <= y:
            out.add((x, y))
    return out


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

    def test_five(self):
        """Verify five"""
        assert decompose_prime(5) == (1, 2)
        with raises(ValueError, match="Could not decompose"):
            decompose_prime(5, 2)
        with raises(ValueError, match="Could not decompose"):
            decompose_prime(5, 3)
        assert decompose_prime(5, 4) == (1, 1)
        assert decompose_prime(5, 5) == (0, 1)


class TestNumberDecomposition:
    """Tests for decompose_number"""

    def test_small_numbers(self):
        """Verify small numbers"""
        max_n = 10_000 if os.getenv("CI") else 50_000

        for n in range(1, max_n + 1):
            got = decompose_number(n)
            expect = {(x, y) for x, y in brute_force_twosquares(n) if x != 0 and y != 0 and x != y}

            assert got == expect, f"Mismatch for n={n}: missing={expect - got}, extra={got - expect}"

    def test_all_small_numbers(self):
        """Verify all solutions for small numbers"""
        max_n = 10_000 if os.getenv("CI") else 50_000

        for n in range(1, max_n + 1):
            got = decompose_number(n, no_trivial_solutions=False)
            expect = brute_force_twosquares(n)

            assert got == expect, f"Mismatch for n={n}: missing={expect - got}, extra={got - expect}"

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
            assert x**2 + y**2 == 19890

    def test_four(self):
        """Verify four"""
        assert decompose_number(4) == set()
        assert decompose_number(4, no_trivial_solutions=False) == {(0, 2)}

    @mark.parametrize(
        "n",
        [
            # Larger hand-picked composites / squares / near-32bit boundary
            10**6,
            10**6 + 1,
            999_999,
            2**31 - 1,
            2**31,
            2**31 + 1,
            2**31 + 12345,
        ],
    )
    def test_all_examples(self, n: int):
        """Validate completeness for the given examples"""
        got = decompose_number(n, no_trivial_solutions=False)
        expect = brute_force_twosquares(n)
        assert got == expect, f"Mismatch for n={n}: missing={expect - got}, extra={got - expect}"

    def test_fuzzed_large_numbers_match_bruteforce(self):
        """Validate completeness for some random examples"""
        rng = random.Random(0)
        count = 20 if os.getenv("CI") else 200

        for _ in range(count):
            n = rng.randrange(1, 2**31 + 100_000)
            got = decompose_number(n, no_trivial_solutions=False)
            expect = brute_force_twosquares(n)
            assert got == expect, f"Mismatch for n={n}: missing={expect - got}, extra={got - expect}"
