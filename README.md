# LAPD_dev
Numerical work towards FEM simulation of Large Plasma Device

Firedrake scripts:

1. SOL_1D_DG_upwind_tdep_irksome.py: this evolves simple system of 1D SOL equations (only density and velocity evolve, system is isothermal).  Uses DG for density and CG for velocity.  Uses an upwind flux scheme for DG.  (Note corresponding output animation is from an old slightly different version of the code but the corrected version gives a very similar output.)
