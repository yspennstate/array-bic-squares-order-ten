# Mutually orthogonal array-bic squares of order ten and of every large order

The paper is `main.tex` (source) and `main.pdf`; archived at Zenodo, doi:10.5281/zenodo.24493269.

`verify_moas.py` builds the three mutually orthogonal array-bic squares of order 7n+3 described in the paper
"Mutually orthogonal array-bic squares of order ten and of every large order" (Y. Shmalo, 2026), and checks
them directly against Definitions 2.16 and 2.17 of S. Adriaensen, J. Bamberg and J. Williford, *On configurations
that are nearly generalised polygons*, arXiv:2611.17855.

The group is F_7 x H, where H is a product of fields F_q with each q a prime at least 5, 8 or 9. Points 1, 2, 3
are the positions; points 4, ..., 7n+3 are the elements of the group.

```
python verify_moas.py              # n = 1, 5, 7, 8, 9, 11, 25 (5x5), 35 (5x7), each with a negative control
python verify_moas.py 5 8          # one case, H = F_5 x F_8 (n = 40, order 283)
python verify_moas.py --print10    # the three squares of order 10
```

The negative control swaps two entries of one square and must fail the check. Python 3, standard library only.
