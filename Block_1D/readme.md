This example is a 1D double integrator example. Two calls to fmincon are used, the first to move the block 2 meters in 2 seconds with a target final velocity of 0 m/s and no constraint on the input.
The second takes this original trajectory and minimizes the total control input used.
