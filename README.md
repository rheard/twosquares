# twosquares

**Deprecated:** `twosquares` has been merged into the `quadint` project as `quadint.sums`.

This package is no longer the recommended way to decompose numbers into sums of two squares. New code should install and import from `quadint` instead:

```bash
python -m pip install quadint
```

```python
from quadint.sums import decompose_prime, decompose_number

print(decompose_prime(19889))
# (17, 140)

print(decompose_number(19890))
# {(69, 123), (57, 129), (3, 141), (87, 111)}
```

## Notes

An explanation of how this package was originally derived can still be found [here](https://github.com/rheard/quadint/tree/main/quadint/sums).