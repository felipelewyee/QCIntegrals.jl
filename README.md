# QCIntegrals.jl

QCIntegrals is a simple package written in Julia for computing common Integrals for Quantum Chemistry calculations.

## Basics

The following code defines the geometry of a molecule in Bohrs, and compute one-electron integrals.

```julia

mol = QCIntegrals.build_molecule("""
O   0.000000  0.000000  0.000000
H   0.277400  0.892900  0.254400
H   0.606800 -0.238300 -0.716900
""")

shells = QCIntegrals.build_basis(molecule,"aug-cc-pvdz")

S = QCIntegrals.build_S(shells)
T = QCIntegrals.build_T(shells)
V = QCIntegrals.build_V(shells,mol.Zs,mol.coords)
```

Two-electron integrals can also be computed
```julia
I4 = QCIntegrals.build_I4(shells)

shells_aux = QCIntegrals.build_basis(mol,"aug-cc-pvdz-jkfit",normalized=false,auxiliar=true)
I2 = QCIntegrals.build_I2(shells_aux)
I3 = QCIntegrals.build_I3(shells,shells_aux)
```
