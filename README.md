# TransmissionLossThroughWall
Solution using MATLAB to pressure wave equation with transmission loss of pressure through a wall, solved on a square domain. The purpose of this code is to emulate the effect of sound polution from roads near to someones home. The pressure is initialized as quiescent except the top and right boundaries, with a periodic magnitude condition using a sine function to vary the pressure on the boundary with the hope to make it more akin to having traffic on a road. We use a central-time-central-space discretization of the wave equation for the propogation of pressure in the following way:
$$P_{i,j}^{n+1} = 2P_{i,j}^n - P_{i,j}^{n-1} + C^2 \left( P_{i+1,j}^n + P_{i-1,j}^n + P_{i,j+1}^n + P_{i,j-1}^n - 4P_{i,j}^n \right)$$
We note that we expressed our wave equation in a non-dimensional form $\partial_{tt}\overline P = \overline\nabla^2\overline P$, using dimensionless parameters $\overline P = P/P_0$, $\overline\nabla = L\nabla$, and $\overline t = (c/L) t$. The transmission loss through the wall is dictated by:
$$TL = 10\log_10 \left( 1 + \left( \frac{\rho_w h \omega}{2\rho_0 c} \right)^2 \right)$$
