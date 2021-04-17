This github contains some of the code and data for the paper "Improved bounds on entropy production in living systems".

EntropyEst.m is the core function that performs the optimization given the measured statistics. There are tolerances in there, such as constraint and function tolerances, that may need to be adjusted. If it can't converge it runs again with a lower constraint tolerance (by increasing tolerances the entropy production bound decreases). Alternatively one could run the optimization multiple times until convergence is reached.

BiasedWalker.m is a good example script to see how to use EntropyEst.m and reproduces some of the results from Fig 2.

The following scripts relate to the experimental figures. They will take a while to run, in particular the bootstrap results were run with a very high n for the figures, running them for lower n will get roughly, but not exactly the same results.

FlagellaMotor.m shows how to get results from Fig 3. The data files are in file_sruct.mat.

MTBootstrap.m shows how to get results from Fig 4. The data files are in MT_traj.mat

HEKOsc.m shows how to get results from Fig 4. The data files are in HEK_data.mat

The remaining functions are helper functions and will be called by the above scripts.
