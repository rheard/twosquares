I've defined this process as "decomposing" a number $n$ into the possible solutions $x^2 + y^2 = n$.

# twosquares

`twosquares` exposes two simple methods for this number-theory task:

* `decompose_prime(p)` — find the unique non-negative pair (a, b) with $a^2 + b^2 = p$ for $p = 2$ or any prime $p \equiv 1 \bmod 4$.

* `decompose_number(n)` — find the unique non-negative pairs (a, b) with $a^2 + b^2 = n$ for any number $n$.

### Algorithm
Here is an explanation for the algorithm in this repo:

This builds on the algorithm described by Stan Wagon (1990),
        based on work by Serret and Hermite (1848), and Cornacchia (1908).

The first thing to know is that there is a deterministic algorithm to quickly find the 1 and only solution for primes $p$ that are $\equiv 1 \bmod 4$. 
    The solution is to simply use Euclid's algorithm. This is the heart of the general algorithm for any $n$, which is described here:

1. Factor $n$. This is the hardest, most time-consuming step. 
    We'll say this factoring has the form $2^t * (p_1 ^ {k_1} * p_2 ^ {k_2} * ...) * (q_1 ^ {j_1} * q_2 ^ {j_2} * ...)$
    Where $t$ is the 2's exponent, all of the primes $p$ are $\equiv 1 \bmod 4$, and all of the primes $q$ are $\equiv 3 \bmod 4$

   1. If any of primes $q$ (which are $\equiv 3 \bmod 4$) have an exponent $j$ that is odd, then there are no solutions.
   2. If there are no primes $p$ (that are $\equiv 1 \bmod 4$), there are no solutions. 

2. Construct a "base number" to use during the combinatorics later. We'll start with the 2's power, we will start the base number off as $(1 - i)^t$

3. Now, for each prime $q^j$ in the factoring (the ones $\equiv 3 \bmod 4$):

   1. Multiply the base number by $(-qi)^{ j / 2 }$ Where $-qi$ is an imaginary number with $0$ real and $-q$ complex part, and it is raised to the power of $j / 2$. 
        Note that from the rules laid out in step 1 above, $j$ is guaranteed to be even so this will always be an integer.

4. Next will begin the combinatorics for the $p$ group, however 1 member of the $p$ group does not need to engage in this combinatorics (reasons why below). 
       So we'll select the first prime $p \equiv 3 \bmod 4$, and create it's "imaginary decomposition". This is a complex number $x+yi$ made from the solution $x^2 + y^2 = p$. 
       Multiply the base number by this number too.

   1. If the exponent for this $p$ (the $k$ value) was 1 then $p$ can be removed entirely from the group. 
   2. If $k$ was greater than $1$, it should be decremented, and the remaining instances of $p$ in the factorization still need to undergo the following combinatorics.

5. Now we need to use combinations of (`True`, `False`) of length $\sum k$ to drive the combinatorics going forward. 
       Python has a product method for this, or you can simply count up using binary numbers to `1<<sum(k)` and look at the bits of this counter. 

   For every possible combination of true/false called "choices":
   1. Start this solution with the base number.
   2. For each factor $p$ left, one time for each exponent $k$:

      1. Get the next "choice" (true/false)
      2. Get the "imaginary decomposition" of the factor, either $x+yi$ if the choice was true or $x-yi$ if the choice was false 
      3. Multiply the total number by this either positive or negative imaginary decomposition

   3. The real and imaginary part of the total number now constitute a solution for $x^2 + y^2 = n$! Amazing!!

      1. The numbers are then sorted so that $x<y$, and this is a solution that may or may not have been found already.

Doing this we can rapidly break any number $n$ up into all of its possible $x^2 + y^2$ solutions!

Note: Remember how we didn't do combinatorics for a single exponent of the $p$ group? 
    If that exponent of that base were included, we would simply get back the same solutions but in reverse order.  
    The problem uses addition and is associative, so we do not care about order. 
    Thus we can effectively halve the compute time with the combinatorics by excluding the other half of that particular exponent.

### Example

As an example, lets look at $n = 19890$. 
    This number's factorization looks like: $2 * 3^2 * 5 * 13 * 17$. 
    The set of primes $p$ (that are $\equiv 1 \bmod 4$) are $5$, $13$ and $17$ (all having an exponent $k$ of $1$). 
    The set of primes $q$ (that are $\equiv 3 \bmod 4$) is $3$, with the only exponent $j$ being 2.

Starting with the rules we can see that all of the primes $q$ have an even exponent ($2$ in this case) and there are primes in the $p$ group, so there must be at least 1 solution for this number!

The 2's exponent is $1$, so we will say our base number is $(1-i)^1$ or just $1 - i$

There is only a single prime in the $q$ group, so we will multiply the base number by $(-3i)^1$ which gives $-3 - 3i$

We'll take the first prime in the $p$ group ($5$) and decompose that number, we find we get $5 = 1^2 + 2^2$. 
    We use this composition to construct a complex number $1 + 2i$. Now multiply the base number by this. 
    This is the real base number, which is $3 - 9i$

Now for combinatorics to produce all the different solutions. The remaining 2 primes $p$ have the following decompositions:

$13 = 2^2 + 3^2$

$17 = 1^2 + 4^2$

1. Using the positive values: Multiply the base number by $2+3i$ and $1+4i$ (from the solutions for both of these primes), 
    and we get $69 + 123i$, which is magically our first solution: $19890 = 69^2 + 123^2$

2. Using the negative values: Multiply the base number by $2-3i$ and $1-4i$, and we get $-129 + 57i$, 
   which using the absolute values of these gives our next solution: $19890 = 129^2 + 57^2$

3. Using positive for 13 and negative for 17: Multiply the base number by $2+3i$ and $1-4i$, and we get $-3 - 141i$, which again gives our next solution: $19890 = 3^2 + 141^2$

4. Lastly it should be obvious: Multiply the base number by $2-3i$ and $1+4i$, we get our final answer: $19890 = 87^2 + 111^2$

So in conclusion, the following are all equal to $19890$: $69^2 + 123^2$, $129^2 + 57^2$, $3^2 + 141^2$, $87^2 + 111^2$

Everything laid out in this example is performed by the following Python code:

```python
from twosquares import decompose_number

print(decompose_number(19890))  # prints {(69, 123), (57, 129), (3, 141), (87, 111)} 
```

## Advanced Usage

It is possible to skip the factoring step and instead build a factored dictionary and pass that to `decompose_number` instead. 
    `decompose_number({2: 1, 3: 2, 5: 1, 13: 1, 17: 1})` will produce the same result as `decompose_number(19890)`.
    If the factoring dictionary has been carefully crafted to comply with the rules laid out in step 1 of the algorithm section above,
    then `limited_checks=True` can be passed as an argument to skip these validations.

It is an interesting fact that the upper bound of solutions can be quickly computed by $\prod (j + 1)$.
    If a minimum number of solutions is required, this can be quickly validated by passing `check_count` as an integer.
