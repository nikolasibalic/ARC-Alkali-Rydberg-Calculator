
ARC (Alkali.ne Rydberg Calculator)
==================================


ARC (Alkali.ne Rydberg Calculator)  is package of routines written in Python, using object-oriented programming (OOP) to make modular, reusable and extendable collection of routines and data for performing useful calculations of single atom and two-atom properties, like level diagrams, interactions and transition strengths for alkali metal and divalent atoms.

Start by installing the latest version of the ARC package calling Python pip ([see prerequisites](https://arc-alkali-rydberg-calculator.readthedocs.io/en/latest/installation.html)) from the command line:

```pip install ARC-Alkali-Rydberg-Calculator ```

For documentation see [online documentation on Read The Docs](http://arc-alkali-rydberg-calculator.readthedocs.io). 

For examples of use check [IPython example notebooks](https://arc-alkali-rydberg-calculator.readthedocs.io/en/latest/getting_started.html#ipython-notebook-with-examples).

For online access to a selection of package functions see [online Atom Calculator](https://atomcalc.org).

If you want to contribute to the project, [check this page](https://arc-alkali-rydberg-calculator.readthedocs.io/en/latest/contribute.html).

![Documentation Status](https://readthedocs.org/projects/arc-alkali-rydberg-calculator/badge/?version=latest) [![PyPI Linux](https://github.com/nikolasibalic/ARC-Alkali-Rydberg-Calculator/actions/workflows/pypi_linux.yaml/badge.svg)](https://github.com/nikolasibalic/ARC-Alkali-Rydberg-Calculator/actions/workflows/pypi_linux.yaml) [![PyPI Windows](https://github.com/nikolasibalic/ARC-Alkali-Rydberg-Calculator/actions/workflows/pypi_windows.yaml/badge.svg)](https://github.com/nikolasibalic/ARC-Alkali-Rydberg-Calculator/actions/workflows/pypi_windows.yaml) [![PyPI macos](https://github.com/nikolasibalic/ARC-Alkali-Rydberg-Calculator/actions/workflows/pypi_macos.yaml/badge.svg)](https://github.com/nikolasibalic/ARC-Alkali-Rydberg-Calculator/actions/workflows/pypi_macos.yaml) [![PyPI version](https://badge.fury.io/py/ARC-Alkali-Rydberg-Calculator.svg)](https://badge.fury.io/py/ARC-Alkali-Rydberg-Calculator)  [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/nikolasibalic/ARC-Alkali-Rydberg-Calculator.git/master?urlpath=lab%2Ftree%2Fdoc%2FRydberg_atoms_a_primer_notebook.ipynb)

-------
Authors
-------

[Nikola Šibalić](https://github.com/nikolasibalic), [Elizabeth J. Robertson](https://www.heibrids.berlin/people/doctoral-researchers/elizabeth-robertson/), [Jonathan D. Pritchard](http://photonics.phys.strath.ac.uk/people/dr-jonathan-pritchard/), [Robert M. Potvliege](https://www.durham.ac.uk/staff/r-m-potvliege/), [Matthew P. A. Jones](https://www.durham.ac.uk/staff/m-p-a-jones/), [Charles S. Adams](https://www.durham.ac.uk/staff/c-s-adams/), [Kevin J. Weatherill](https://www.durham.ac.uk/staff/k-j-weatherill/) and [contributors](https://github.com/nikolasibalic/ARC-Alkali-Rydberg-Calculator/graphs/contributors)

**Please cite as:** 

The simplest way to obtain correct reference(s), given the number of contributions,
is to **call - at the end of your Python script that uses ARC - following function**
```python
from arc import *
# use ARC
print(getCitationForARC())
```
The `getCitationForARC()` will print references that introduced methods
you used into ARC library. Otherwise, you can do manual decision making
based on the logic below:

If you use alkali atoms:
N. Šibalić, J. D. Pritchard, K. J. Weatherill, C. S. Adams,
ARC: An open-source library for calculating properties of alkali Rydberg atoms,
*Computer Physics Communications* **220**, 319 (2017), [https://doi.org/10.1016/j.cpc.2017.06.015](https://doi.org/10.1016/j.cpc.2017.06.015)

If you use divalent atoms (Sr, Ca, Yb ...) or new featutures from ARC 3.0:
E. J. Robertson, N. Šibalić, R. M. Potvliege, M. P. A. Jones,
ARC 3.0: An expanded Python toolbox for atomic physics calculations, *Computer Physics Communications* **261**, 107814 (2021) [https://doi.org/10.1016/j.cpc.2020.107814](https://doi.org/10.1016/j.cpc.2020.107814)

**In addition to main reference above**: If you are using modules from `arc.advanced` please
[paper](arc/advanced/README.md) that introduced relevant ARC extension.

If you are using modules AC Stark calculations `ShirleyMethod` or `RWAStarkShift` please **also cite** [D. H. Meyer, Z. A. Castillo, K. C. Cox, P. D. Kunz, J. Phys. B: At. Mol. Opt. Phys., 53, 034001 (2020) ](https://doi.org/10.1088/1361-6455/ab6051).

**Who are the users of this library?** Check papers that cite us in [ADS](https://ui.adsabs.harvard.edu/abs/2017CoPhC.220..319S/citations) and [Google Scholar](https://scholar.google.com/scholar?cites=3162548955488940394&as_sdt=2005&sciodt=0,5&hl=en).

-------
License
-------

All the files distributed with this program are provided subject to the
BSD-3-Clause license. A copy of the license is provided.
