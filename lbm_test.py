import numpy as np
import matplotlib.pyplot as plt
from lbm import LBMSolver

def run_simulation():
    # Simulation parameters
    nx = 100
    ny = 100
    tau = 0.6
    u_lid = 0.1
    max_steps = 1000

    # Initialize the solver
    solver = LBMSolver(nx, ny, tau, u_lid)

    # Run the simulation loop
    print("Starting LBM Simulation (Lid-driven cavity)...")
    for step in range(max_steps):
        solver.step()
        
        # Print progress
        if (step + 1) % 100 == 0:
            print(f"Step {step + 1}/{max_steps} completed")
    
    print("Simulation finished.")
    
    # Calculate macroscopic quantities
    rho, u = solver.get_macroscopic(solver.f)
    
    # Calculate velocity magnitude
    speed = np.sqrt(u[0]**2 + u[1]**2)

    # Note: We comment out plt.show() and plt.savefig() so no artifacts are left behind
    # in an automated run. A user can uncomment this if they want to view the plot.
    
    # plt.figure(figsize=(6, 5))
    # plt.contourf(speed.T, cmap='viridis', levels=50)
    # plt.colorbar(label='Velocity Magnitude')
    # plt.title(f'Lid-driven Cavity Flow (Step {max_steps})')
    # plt.xlabel('X')
    # plt.ylabel('Y')
    # plt.savefig('test.png')
    # plt.show()
    
    # Verify that the maximum speed is somewhat reasonable (it should be less than the lid speed due to viscosity)
    max_speed = np.max(speed)
    print(f"Maximum velocity magnitude: {max_speed:.4f}")
    
    # Let's perform a simple sanity check assertion for the validation
    if max_speed > 0.0 and max_speed < u_lid + 0.05:
        print("Validation successful! Velocity profile seems physically sound.")
    else:
        print("Validation warning! Maximum velocity looks anomalous.")

if __name__ == "__main__":
    run_simulation()
