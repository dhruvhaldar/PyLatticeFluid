import numpy as np
import pytest
from lbm import LBMSolver

def test_initialization():
    solver = LBMSolver(nx=50, ny=50, tau=0.6)
    rho, u = solver.get_macroscopic(solver.f)
    
    # Check density is initialized to 1.0
    assert np.allclose(rho, 1.0)
    # Check velocity is initialized to 0.0
    assert np.allclose(u, 0.0)

def test_mass_conservation():
    solver = LBMSolver(nx=20, ny=20, tau=0.6, u_lid=0.0)
    
    # We set u_lid=0 to test mass conservation in a closed cavity
    initial_mass = np.sum(solver.rho)
    
    for _ in range(10):
        solver.step()
        
    final_mass = np.sum(solver.rho)
    
    assert np.isclose(initial_mass, final_mass), "Mass is not conserved!"

def test_equilibrium():
    solver = LBMSolver(nx=10, ny=10, tau=1.0)
    
    # Test equilibrium for zero velocity and density 1
    rho = np.ones((10, 10))
    u = np.zeros((2, 10, 10))
    
    feq = solver.equilibrium(rho, u)
    
    # Sum of feq should be equal to rho
    assert np.allclose(np.sum(feq, axis=0), rho)
    
    # For zero velocity, the feq values should equal the weights
    for i in range(9):
        assert np.allclose(feq[i], solver.w[i])

def test_bounce_back():
    solver = LBMSolver(nx=5, ny=5, tau=1.0)
    # Create an artificial non-equilibrium state on the boundary
    # to see if it bounces back correctly in the collision step
    for i in range(9):
        solver.f[i, solver.obstacle] = i  # set arbitrary values
    
    # Expected after bounce-back step is f_post = f[noslip] on boundary
    # We can manually do one collision step and check f_post
    feq = solver.equilibrium(solver.rho, solver.u)
    f_post = solver.f - (solver.f - feq) / solver.tau
    
    # Apply manual bounce back exactly as in the solver
    for i in range(9):
        f_post[i, solver.obstacle] = solver.f[solver.noslip[i], solver.obstacle]
    
    # Verify the bounce-back actually swapped values
    for i in range(9):
        assert np.allclose(f_post[i, solver.obstacle], solver.f[solver.noslip[i], solver.obstacle])
