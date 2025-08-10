import math

from functools import cache
from itertools import product
from typing import Optional, Union

from complexint import complexint
from sympy import factorint

def _euclids_algorithm(a: int, b: int, c: int) -> Optional[tuple[int, int]]:
    """Runs Euclid's algorithm and yields remainders"""
    first = None
    while a != 1:
        r = a % b
        a, b = b, r
        if r > c:
            continue
        if not r:
            return None
        if first is not None:
            return r, first
        first = r

    raise ValueError("Euclid's algorithm failed")  # This will never be reached, for ruff


@cache
def decompose_prime(p: int) -> tuple[int, int]:
    """
    Decompose a prime number into a**2 + b**2

    There will be at most 1 solution for primes, but only if the prime is equal 1 mod 4 according
        to Fermat's theorem on sums of two squares.

    This is based on the algorithm described by Stan Wagon (1990),
        based on work by Serret and Hermite (1848), and Cornacchia (1908)

    Returns:
        tuple<int, int>: a and b

    Raises:
        ValueError: If p cannot be decomposed because it is 3 mod 4
    """
    if p % 4 != 1:
        if p == 2:
            return 1, 1

        raise ValueError(f'Could not decompose {p!r}')

    p_sqrt = math.isqrt(p)
    for a in range(1, p):  # a must be co-prime to p
        if pow(a, (p - 1) // 2, p) == p - 1:
            # Found a quadratic non-residue of p! (a)
            #   Now run the Euclidean algorithm with a and p
            res = _euclids_algorithm(p, pow(a, (p - 1) // 4, p), p_sqrt)
            if res:
                return res

    raise ValueError(f'Could not decompose {p!r}')


def decompose_number(n: Union[dict, int],
                     check_count: Optional[int] = None,
                     *,
                     limited_checks: bool = False) -> set[tuple[int, int]]:
    """
    Decompose any number into all the possible x**2 + y**2 solutions

    There may be many solutions.

    Args:
        n (int, dict): The number to decompose. Can be an integer which will be factored,
            or the already factored number.
        check_count (int): If provided, and it is predicted that a number will have fewer than this many solutions,
            that number is skipped and an empty list is returned instead.
        limited_checks (bool): Only run limited checks. Should only be used with prepared input
            or false positive will appear.

    Returns:
        set<tuple<int, int>>: All unique solutions (x, y)
    """

    # Step 1: Factor n. This is the most time consuming step, especially on larger numbers. Avoid if possible
    factors = n if isinstance(n, dict) else factorint(n)

    # Look for shortcuts
    if len(factors) == 1 and sum(factors.values()) == 1:
        # Only 1 factor with a power of 1, this is a prime number
        p = next(iter(factors))
        if check_count and check_count > 1:
            return set()  # There will only be 1 solution. If check_count is greater than that, do nothing
        if p % 4 != 1:
            return set()  # Primes == 1 mod 4 have no solutions

        # Return the single decomposition for this prime...
        return {decompose_prime(p)}

    # Divide the factors into the ones equivalent to 1 mod 4, and 3 mod 4:
    p_1, p_3 = {}, {}
    for p, k in factors.items():
        p_mod_4 = p % 4
        if p_mod_4 == 1:
            p_1[p] = k
        elif p_mod_4 == 3:
            p_3[p] = k

    if not limited_checks and any(k % 2 == 1 for k in p_3.values()):
        # There is a prime == 3 mod 4 and at least one has an odd exponent, so no results
        return set()

    if not p_1:
        # There aren't any primes = 1 mod 4, in which case, escape to no results
        return set()

    if check_count and math.prod(f + 1 for f in p_1.values()) < check_count:
        # Provided a check_count and the expected number of results is less than that
        #   expected number = (f_1 + 1) * (f_2 + 1) * ...
        #   where f is an exponent of a prime = 1 mod 4 in the factorization
        return set()

    two_power = factors.get(2, 0)  # 2 is a special case. Get that exponent
    factors = p_1  # Convenience. These are the only factors we care about going forward...

    # Handle the primes == 3 mod 4: For each one p^k, take -p*j ** max(k // 2, 1),
    #   then multiply that all together
    p_3_coefficients = math.prod(complexint(imag=-p)**max(k // 2, 1) for p, k in p_3.items())

    p_decompositions = {p: decompose_prime(p) for p in factors}
    # Here we create the a+bj and a-bj pairs...
    p_exp_decompositions = {p: (complexint(*d), complexint(d[0], -d[1])) for p, d in p_decompositions.items()}
    first_p = next(iter(factors))
    factors[first_p] -= 1  # subtract 1 for the base item

    # Add the 2's power special case to the p_3_coefficients for later...
    p_3_coefficients *= complexint(1, -1) ** two_power

    # Base item only needs the positive value (a+bj) and doesn't need to vary
    base_item = p_exp_decompositions[first_p][0]
    base_item *= p_3_coefficients  # Add the previously calculated p_3 coefficients (including the 2's power)
    found = set()
    for choices in product([0, 1], repeat=sum(factors.values())):  # This will run ONCE if repeat=0
        # Initial total for this loop with the base item
        total = base_item
        choice_i = 0
        for p, k in factors.items():
            # Get a choice for each factor, on whether to use (a+bj) or (a-bj) here...
            for _ in range(k):
                plus_or_minus = choices[choice_i]
                total *= p_exp_decompositions[p][plus_or_minus]
                choice_i += 1
        sol = tuple(sorted((abs(total.real), abs(total.imag))))
        if sol[0] == sol[1]:
            continue  # Skip symmetrical solutions with repeat numbers (a**2 + a**2)
        if any(x == 0 for x in sol):
            continue  # Skip solutions containing 0
        found.add(sol)
    return found
