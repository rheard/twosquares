import math

from functools import cache
from itertools import product
from typing import Optional, Union

from sympy import factorint as sympy_factorint
from quadint import complexint

try:
    from cypari import pari

    CYPARI = True
except ImportError:
    pari = None
    CYPARI = False


def factorint(n: int) -> dict[int, int]:
    if CYPARI:
        # PARI returns a matrix: [p1 e1; p2 e2; ...]
        f = pari(n).factor()

        # Convert to Python dict[int, int]
        return {
            int(p): int(e) for p, e in zip(f[0], f[1])
        }
    else:
        return sympy_factorint(n)


def _mod_sqrt_prime(n: int, p: int) -> Optional[int]:
    """Return x such that x*x % p == n % p, or None if no sqrt exists. p must be prime."""
    n %= p
    if n == 0:
        return 0

    if p == 2:
        return n

    # Legendre symbol test: residue iff n^((p-1)/2) == 1 (mod p)
    if pow(n, (p - 1) // 2, p) != 1:
        return None

    # Fast path when p ≡ 3 (mod 4)
    if p % 4 == 3:
        return pow(n, (p + 1) // 4, p)

    # Tonelli-Shanks
    q = p - 1
    s = 0
    while q % 2 == 0:
        s += 1
        q //= 2

    z = 2
    while pow(z, (p - 1) // 2, p) != p - 1:
        z += 1

    m = s
    c = pow(z, q, p)
    t = pow(n, q, p)
    r = pow(n, (q + 1) // 2, p)

    while t != 1:
        i = 1
        t2i = (t * t) % p
        while i < m and t2i != 1:
            t2i = (t2i * t2i) % p
            i += 1

        b = pow(c, 1 << (m - i - 1), p)
        r = (r * b) % p
        c = (b * b) % p
        t = (t * c) % p
        m = i

    return r


def _euclids_algorithm(a: int, b: int, c: int) -> Optional[int]:
    """Runs Euclid's algorithm and returns remainder"""
    while b > c:
        r = a % b
        a, b = b, r
        if not b:
            return None

    return b


@cache
def decompose_prime(p: int, d: int = 1) -> tuple[int, int]:
    """
    Decompose a prime number into a**2 + d * b**2

    There will be at most 1 solution for primes.
        If d == 1, this will only be if the prime is equal 1 mod 4 according to Fermat's theorem on sums of two squares.

    This is based on the algorithm described by Stan Wagon (1990),
        based on work by Serret and Hermite (1848), and Cornacchia (1908)

    Returns:
        tuple<int, int>: a and b

    Raises:
        ValueError: If p cannot be decomposed because it is 3 mod 4
    """
    if d < 1:
        raise ValueError(f"d must be >= 1, got {d!r}")

    if p == 2:
        if d == 1:
            return 1, 1
        if d == 2:
            return 0, 1

        raise ValueError(f"Could not decompose {p!r} with d={d!r}")

    # If sqrt(-d) mod p doesn't exist, no solution for this prime
    t = _mod_sqrt_prime((-d) % p, p)  # Note this replaces the powers in the old algorithm
    if t is None:
        raise ValueError(f"Could not decompose {p!r} with d={d!r}")

    def _try_cornacchia_root(p: int, d: int, t: int) -> Optional[tuple[int, int]]:
        p_sqrt = math.isqrt(p)

        x = _euclids_algorithm(p, t, p_sqrt)
        if x is None:
            return None

        rhs = p - x * x
        if rhs < 0 or rhs % d != 0:
            return None

        y2 = rhs // d
        y = math.isqrt(y2)
        if y * y != y2:
            return None

        return abs(x), abs(y)

    res = _try_cornacchia_root(p, d, t)
    if res is None:
        res = _try_cornacchia_root(p, d, (p - t) % p)

    if res is None:
        raise ValueError(f"Could not decompose {p!r} with d={d!r}")

    if res[0] > res[1] and d == 1:
        return res[1], res[0]

    return res


def decompose_number(n: Union[dict, int],
                     check_count: Optional[int] = None,
                     *,
                     limited_checks: bool = False,
                     no_trivial_solutions: bool = True) -> set[tuple[int, int]]:
    """
    Decompose any number into all the possible x**2 + y**2 solutions.

    There may be many solutions.

    Args:
        n (int, dict): The number to decompose. Can be an integer which will be factored,
            or the already factored number.
        check_count (int): If provided, and it is predicted that a number will have fewer than this many solutions,
            that number is skipped and an empty list is returned instead.
        limited_checks (bool): Only run limited checks. Should only be used with prepared input
            or false positive will appear.
        no_trivial_solutions (bool): Exclude trivial solutions? Defined as any symmetrical solution, or any
            solution with 0. Essentially excludes perfect squares and doubles of perfect squares.
            Note that a value of False will make the algorithm quite a bit slower.

    Returns:
        set<tuple<int, int>>: All unique solutions (x, y)
    """

    # Step 1: Factor n. This is the most time consuming step, especially on larger numbers. Avoid if possible
    factors: dict[int, int] = n if isinstance(n, dict) else factorint(n)

    # Look for shortcuts
    if len(factors) == 1 and sum(factors.values()) == 1:
        # Only 1 factor with a power of 1, this is a prime number
        p = next(iter(factors))
        if check_count and check_count > 1:
            return set()  # There will only be 1 solution. If check_count is greater than that, do nothing
        if p % 4 == 3:
            return set()  # Primes == 1 mod 4 have no solutions
        if p == 0:
            return {(0, 0)}
        if p == 2 and no_trivial_solutions:
            # The solution for 2 is the only prime trivial solution (1, 1)
            return set()

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

    if (not limited_checks or no_trivial_solutions) and any(k % 2 == 1 for k in p_3.values()):
        # There is a prime == 3 mod 4 and at least one has an odd exponent, so no results
        return set()

    if not p_1 and no_trivial_solutions:
        # There aren't any primes = 1 mod 4, in which case, escape to no results
        return set()

    if check_count and math.prod(f + 1 for f in p_1.values()) < check_count:
        # Provided a check_count and the expected number of results is less than that
        #   expected number = (f_1 + 1) * (f_2 + 1) * ...
        #   where f is an exponent of a prime = 1 mod 4 in the factorization
        return set()

    two_power = factors.get(2, 0)  # 2 is a special case. Get that exponent
    factors = p_1  # Convenience. These are the only factors we care about going forward...

    # Handle the primes == 3 mod 4: For each one p^k, take p*j ** max(k // 2, 1),
    #   then multiply that all together
    p_3_coefficients = math.prod(p**(k // 2) for p, k in p_3.items())

    p_decompositions = {p: decompose_prime(p) for p in factors}
    # Here we create the complex numbers and their conjugates...
    p_exp_decompositions = {p: (complexint(*d), complexint(d[0], -d[1])) for p, d in p_decompositions.items()}

    # Add the 2's power special case to the p_3_coefficients for later...
    p_3_coefficients *= complexint(1, 1) ** two_power

    base_item = p_3_coefficients  # Add the previously calculated p_3 coefficients (including the 2's power)

    if no_trivial_solutions:
        first_p = next(iter(factors))
        factors[first_p] -= 1  # subtract 1 for the base item

        # Base item only needs the positive value (a+bj) and doesn't need to vary
        base_item *= p_exp_decompositions[first_p][0]

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
        sol: tuple[int, int] = abs(total.real), abs(total.imag)
        if no_trivial_solutions and (sol[0] == sol[1] or sol[0] == 0 or sol[1] == 0):
            continue
        if sol[1] < sol[0]:
            sol = (sol[1], sol[0])
        found.add(sol)
    return found
