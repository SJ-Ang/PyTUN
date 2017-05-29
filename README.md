# PyTUN
Tunneling Correction to Chemical Kinetics using Eckart Potential 

Equations obtained from "Garrett, B. C., & Truhlar, D. G. (1979). Semiclassical tunneling calculations. Journal of Physical Chemistry, 83(22), 2921-2926."

Reads Gaussian09 output files to generate "Tunneling Correction" κ and Apparent Free Energy Barrier

Integration to obtain κ was calculated using 10-point Gauss-Legendre quadrature(Exact solution for polynomials of order up till 21)

Zero-point energy corrections are added to each species in the calculation.

User has to provide 6 files placed in the same directory as PyTUN.py:
1. reactant.out (Frequency calculation of reactant)
2. reactant_sp.out (Single-point calculation of reactant usually at a higher level of theory)
3. product.out (Frequency calculation of product)
4. product_sp.out (Single-point calculation of product)
5. ts.out (Frequency calculation of transition state)
6. ts_sp.out (Single-point calculation of transition state)

'''Python PyTUN.py'''
