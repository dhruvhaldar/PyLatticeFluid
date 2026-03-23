import numpy as np

class LBMSolver:
    def __init__(self, nx, ny, tau, u_lid=0.1):
        self.nx = nx
        self.ny = ny
        self.tau = tau
        self.u_lid = u_lid
        
        # D2Q9 lattice constants
        # Directions:
        # 6 2 5
        # 3 0 1
        # 7 4 8
        self.c = np.array([
            [0,  1,  0, -1,  0,  1, -1, -1,  1], # x-components
            [0,  0,  1,  0, -1,  1,  1, -1, -1]  # y-components
        ])
        
        self.w = np.array([
            4/9, 1/9, 1/9, 1/9, 1/9, 1/36, 1/36, 1/36, 1/36
        ])
        
        # Opposite directions for bounce-back
        self.noslip = [0, 3, 4, 1, 2, 7, 8, 5, 6]
        
        # Initialize macroscopic variables
        self.rho = np.ones((nx, ny))
        self.u = np.zeros((2, nx, ny))
        
        # Solid boundaries (stationary walls)
        self.obstacle = np.zeros((nx, ny), dtype=bool)
        self.obstacle[:, 0] = True   # Bottom wall
        self.obstacle[0, :] = True   # Left wall
        self.obstacle[-1, :] = True  # Right wall
        
        # Moving lid
        self.lid = np.zeros((nx, ny), dtype=bool)
        self.lid[:, -1] = True       # Top wall (lid)
        
        # Initialize distributions to equilibrium
        self.f = self.equilibrium(self.rho, self.u)
        
    def equilibrium(self, rho, u):
        """Calculates equilibrium distribution."""
        usqr = 1.5 * (u[0]**2 + u[1]**2)
        feq = np.zeros((9, self.nx, self.ny))
        for i in range(9):
            cu = 3.0 * (self.c[0,i] * u[0] + self.c[1,i] * u[1])
            feq[i] = rho * self.w[i] * (1.0 + cu + 0.5 * cu**2 - usqr)
        return feq
        
    def get_macroscopic(self, f):
        """Computes density and velocity from distributions."""
        rho = np.sum(f, axis=0)
        u = np.zeros((2, self.nx, self.ny))
        for i in range(9):
            u[0] += self.c[0,i] * f[i]
            u[1] += self.c[1,i] * f[i]
        
        # Avoid division by zero
        rho_safe = np.where(rho > 0, rho, 1.0)
        u[0] /= rho_safe
        u[1] /= rho_safe
        return rho, u
        
    def step(self):
        """Performs one LBM time step."""
        # 1. Macroscopic variables
        self.rho, self.u = self.get_macroscopic(self.f)
        
        # 2. Collision
        feq = self.equilibrium(self.rho, self.u)
        f_post = self.f - (self.f - feq) / self.tau
        
        # Boundary condition: Bounce-back for stationary obstacles
        # To avoid the double-collision bug, we bounce back from self.f (pre-collision)
        for i in range(9):
            f_post[i, self.obstacle] = self.f[self.noslip[i], self.obstacle]
            
        # Boundary condition: Moving lid with Zou/He style momentum injection
        for i in range(9):
            # Applying momentum transfer according to Zou/He style moving wall
            # Bounce back + momentum transfer
            f_post[i, self.lid] = self.f[self.noslip[i], self.lid] - 6.0 * self.w[i] * self.rho[self.lid] * (self.c[0,i] * self.u_lid)
            
        # 3. Streaming
        # Note: roll with negative shifts because if particle goes positive x,
        # its index i needs to be shifted positively to end up at x+1 (so we roll by c).
        # Wait, if we roll f_post[i] by c[0,i], a value at x goes to x+c. That is correct.
        for i in range(9):
            self.f[i] = np.roll(
                np.roll(f_post[i], self.c[0,i], axis=0),
                self.c[1,i], axis=1
            )
