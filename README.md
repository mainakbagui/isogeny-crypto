# Isogeny-based cryptography: code

Mainak Bagui · M.Sc. (Mathematics), Pondicherry University · ORCID [0009-0001-5365-6092](https://orcid.org/0009-0001-5365-6092)

Two self-contained pieces of code on supersingular isogenies.

| folder | what it does | run with |
|---|---|---|
| [`csidh/`](csidh/) | CSIDH key exchange (Castryck–Lange–Martindale–Panny–Renes, 2018): the action of $\mathrm{cl}(\mathbb{Z}[\sqrt{-p}])$ on supersingular Montgomery curves over $\mathbb{F}_p$, with public-key validation. Toy prime $p = 419$ and CSIDH-512. Tests check the action against point counts and the class number $h(-4p)$. | Python 3.8+, no packages |
| [`isogeny-graphs/`](isogeny-graphs/) | Supersingular $\ell$-isogeny graphs over $\mathbb{F}_{p^2}$ and their Brandt matrices $B(\ell)$. Checks that $\operatorname{charpoly} B(\ell)$ equals that of $T_\ell$ on $M_2(\Gamma_0(p))$, the weighted symmetry of $B(\ell)$, and the Ramanujan bound $\lvert\lambda\rvert \le 2\sqrt{\ell}$. An independent PARI/GP version is included. | SageMath ([SageMathCell](https://sagecell.sagemath.org)) or PARI/GP |

```
cd csidh && python3 test_csidh.py && python3 csidh.py
```

Each folder has its own notes on the mathematics and on limitations. The code is for study and research. It is not constant-time and should not be used to protect real data.
