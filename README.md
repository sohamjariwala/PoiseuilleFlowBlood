# Chebyshev Pseudo-spectral solver for Poiseuille blood flow in a microtube
## Model description
The model is developed using a pseudo-spectral scheme where the domain is discretized
using a Chebyshev expansion at finite points across the cross-section of the microtube


No slip
__________________
  .
  .
  .
  .
Chebysehev nodes
  .
  .
  .
  .
  __________________
No slip

Each nodal point is a truncated chebyshev expansion,

v = Sum_i a_i Phi_i(r)

where Phi_(r) is the pseudo-spectral basis function. Different basis functions are used depending on the boundary conditions of the variable being solved.

