# TransmissionLossThroughWall
Solution using MATLAB to pressure wave equation with transmission loss of pressure through a wall, solved on a square domain. The purpose of this code is to emulate the effect of sound polution from roads near to someones home. The pressure is initialized as quiescent except the top and right boundaries, with a periodic magnitude condition using a sine function to vary the pressure on the boundary with the hope to make it more akin to having traffic on a road. We use a central-time-central-space discretization of the wave equation for the propogation of pressure.

$$p_{i,j}^{n+1} = 2*p_{i,j}^n - p_{i,j}^{n-1} + C^2*(p_{i+1,j}^n + p_{i-1,j}^n + p_{i,j+1}^n + p_{i,j-1}^n - 4*p_{i,j}^n)$$
