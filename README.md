# PyTUN

Tunneling Correction to Chemical Kinetics 

Equations obtained from: 

Wigner, E. Z., Über das Überschreiten von Potentialschwellen bei chemischen Reaktionen. 1932. Z. Phys. Chem. B, 19, 203. (Wigner's approximation)

Garrett, B. C., & Truhlar, D. G., Semiclassical tunneling calculations. J. of Phys. Chem., 1979, 83(22), 2921-2926. (Eckart Potential)

R. T. Skodje, D. G. Truhlar and B. C. Garrett, A general small-curvature approximation for transition-state-theory transmission coefficients, J. Phys. Chem., 1981, 85, 3019–3023. (Parabolic Approximation)

Read Gaussian09 output files to generate Tunneling Correction κ (dimensionless), Corrected Rate k<sub>corr</sub> (s<sup>-1</sup>) and Apparent Free Energy Barrier (kJmol<sup>-1</sup>)

Uncorrected rate [k<sub>uncorr</sub>] is obtained using the Eyring-Polanyi Equation, while Apparent Free Energy Barrier is back-transformed from the k<sub>corr</sub> [k<sub>corr</sub>=κk<sub>uncorr</sub>] using the same equation. 

To evaluate κ for Eckart Potential, integration to obtain κ was calculated using 10-point Gauss-Legendre quadrature(Gives exact solution for polynomials of order up till 21)

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
```
python PyTUN.py
```
Sample Output:
```
Uncorrected delta G dagger is 148.392236
Uncorrected rate is 6.2526688665E-14
The Wigner tunneling correction is 1.924509
The Skodje tunneling correction is 3.328375
The Eckart tunneling correction is 3.001622
For Wigner Tunneling:
Corrected rate is 1.2033316318E-13
Apparent delta G dagger is 146.769335
For Skodje Tunneling:
Corrected rate is 2.0811226750E-13
Apparent delta G dagger is 145.411330
For Eckart Tunneling:
Corrected rate is 1.8768150989E-13
Apparent delta G dagger is 145.667483

```
