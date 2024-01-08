# Attempt at time-dependent solver of 1D SOL equations, with upwind flux
# based on:
# https://www.firedrakeproject.org/demos/DG_advection.py.html
# irksome implementation follows
# https://www.firedrakeproject.org/Irksome/demos/demo_cahnhilliard.py.html
# and
# https://www.firedrakeproject.org/Irksome/demos/demo_monodomain_FHN.py.html

from firedrake import *
import math
from irksome import Dt, GaussLegendre, MeshConstant, TimeStepper

meshres = 200
mesh = IntervalMesh(meshres, -1, 1)
V1 = FunctionSpace(mesh, "DG", 1)
V2 = VectorFunctionSpace(mesh, "CG", 1)  # IMPORTANT velocity space needs to be continuous
V = V1*V2

# time parameters
T = 4.0
timeres = 100
t = Constant(0.0)
dt = Constant(T/timeres)

# parameters for irksome
butcher_tableau = GaussLegendre(2)

# model parameters
nstar = Constant(1.0)
diffusion = 1.0e0
Temp = 1.0

x = SpatialCoordinate(mesh)
nu = Function(V)
n, u = split(nu)
#n, u = nu.split()  # doesn't work with Irksome
v1, v2 = TestFunctions(V)

# sonic outflow equilibrium init data
#nu.sub(0).interpolate((nstar/sqrt(Temp))*(1+sqrt(1-x[0]*x[0])))
#nu.sub(1).interpolate((sqrt(Temp)/(x[0]))*(1-sqrt(1-x[0]*x[0])))

# Gaussian blob init data
width = 0.1
nu.sub(0).interpolate((1.0+0.2*(1/sqrt(2*math.pi*width**2))*exp(-x[0]**2/(2*width**2))))
nu.sub(1).interpolate(as_vector([0.0+1.0*x[0]]))

# source function
nstarFunc = Function(V)
nstarFunc.sub(0).interpolate(nstar + 0.0*x[0])

# TRIALCODE check init data
File("SOL_1D_DG_upwind_tdep_irksome_init.pvd").write(nu.sub(0), nu.sub(1))
#quit()

norm = FacetNormal(mesh)
u_n = 0.5*(dot(u,norm)+abs(dot(u,norm)))

# outflow BC imposed weakly in here, penalty term in ds (note conditional not really needed here as both are outflow)

F = -((Dt(n)*v1)*dx + (dot(Dt(u),v2))*dx) \
   + (n*dot(u, grad(v1))+v1*nstar)*dx \
   + (n*u[0]*dot(u, grad(v2[0]))+Temp*n*grad(v2[0])[0])*dx \
   - (v1('+') - v1('-'))*(u_n('+')*n('+') - u_n('-')*n('-'))*dS \
   - (v2('+')[0]-v2('-')[0])*(n('+')*u('+')[0]*u[0]('+')-n('-')*u('-')[0]*u[0]('-'))*dS \
   - (v2('+')[0]-v2('-')[0])*(n('+')-n('-'))*Temp*dS \
   - conditional(dot(u, norm) > 0, v1*dot(u, norm)*n, 0.0)*ds \

# params taken from Cahn-Hilliard example cited above
params = {'snes_monitor': None, 'snes_max_it': 100,
          'snes_linesearch_type': 'l2',
          'ksp_type': 'preonly',
          'pc_type': 'lu', 'mat_type': 'aij',
          'pc_factor_mat_solver_type': 'mumps'}

# Dirichlet BCs are needed for boundary velocity

bc_test1 = DirichletBC(V.sub(1),as_vector([-1]),1)
bc_test2 = DirichletBC(V.sub(1),as_vector([1]),2)

stepper = TimeStepper(F, butcher_tableau, t, dt, nu, solver_parameters=params, bcs=[bc_test1, bc_test2])

outfile = File("SOL_1D_DG_upwind_tdep_irksome.pvd")

while float(t) < float(T):
    if (float(t) + float(dt)) >= T:
        dt.assign(T - float(t))
    outfile.write(nu.sub(0), nu.sub(1))
    stepper.advance()
    t.assign(float(t) + float(dt))
    print(float(t), float(dt))

print("done.")
print("\n")

#ns, us = split(nu)
#File("SOL_1D_DG_upwind_tdep_irksome.pvd").write(ns, us)
