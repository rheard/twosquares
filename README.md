# twosquares
Fast, exact decompositions of integers as sums of two squares.

`twosquares` exposes two simple methods for this number-theory task:

* `decompose_prime(p)` — find the unique nonnegative pair (a, b) with a² + b² = p for p = 2 or prime p ≡ 1 (mod 4).

* `decompose_number(n)` — enumerate all nonnegative pairs (a, b) with a² + b² = n (canonical pairs only).
