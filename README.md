# PyTUN

Tunneling Correction to Chemical Kinetics using Eckart Potential 

Equations obtained from "Garrett, B. C., & Truhlar, D. G. (1979). Semiclassical tunneling calculations. Journal of Physical Chemistry, 83(22), 2921-2926."

Reads Gaussian09 output files to generate Tunneling Correction κ (dimensionless), Corrected Rate k<sub>corr</sub> (s<sup>-1</sup>) and Apparent Free Energy Barrier (kJmol<sup>-1</sup>)

Uncorrected rate [k<sub>uncorr</sub>] is obtained using the Eyring-Polanyi Equation, while Apparent Free Energy Barrier is back-transformed from the k<sub>corr</sub> [k<sub>corr</sub>=κk<sub>uncorr</sub>] using the same equation. 

Integration to obtain κ was calculated using 10-point Gauss-Legendre quadrature(Gives exact solution for polynomials of order up till 21)

Zero-point energy corrections are added to each species in the calculation.

User has to provide 6 files placed in the same directory as PyTUN.py:
1. reactant.out (Frequency calculation of reactant)
2. reactant_sp.out (Single-point calculation of reactant usually at a higher level of theory)
3. product.out (Frequency calculation of product)
4. product_sp.out (Single-point calculation of product)
5. ts.out (Frequency calculation of transition state)
6. ts_sp.out (Single-point calculation of transition state)

Run PyTUN with Python 2.7 (Python 3 is currently not supported)

Simply Type:
```Python
Python PyTUN.py
```
Sample Output:
```
Uncorrected delta G dagger is 80.000000
Uncorrected rate is 5.9960817096E-02
The tunneling correction is 6.500000
Corrected rate is 3.8974531112E-01
Apparent delta G dagger is 75.359883
```
