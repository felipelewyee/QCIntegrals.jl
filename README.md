# QCIntegrals.jl

QCIntegrals is a simple package written in Julia for computing common Integrals for Quantum Chemistry calculations.

## Basics

The following code defines the geometry of a molecule in Bohrs, and compute one-electron integrals.

```julia

molecule = """
O   0.000000  0.000000  0.000000
H   0.277400  0.892900  0.254400
H   0.606800 -0.238300 -0.716900
"""

shells = get_all_shells_from_xyz(molecule,"aug-cc-pvdz")
Zs,coords = get_Z_xyz(molecule)

S = get_S(shells)
T = get_T(shells)
V = get_V(shells,Zs,coords)
```

Two-electron integrals can also be computed
```julia
I4 = get_I4(shells)

shells_aux = get_all_shells_from_xyz(molecule,"aug-cc-pvdz-jkfit",normalized=false,auxiliar=true)
I2 = get_I2(shells_aux)
I3 = get_I3(shells_aux)
```
