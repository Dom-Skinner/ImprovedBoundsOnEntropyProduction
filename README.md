This github contains some of the code and data for the paper "Improved bounds on entropy production in living systems", https://doi.org/10.1073/pnas.2024300118.

EntropyEst.m is the core function that performs the optimization given the measured statistics. There are tolerances in there, such as constraint and function tolerances, that may need to be adjusted. If it can't converge it runs again with a lower constraint tolerance (by increasing tolerances the entropy production bound decreases). Alternatively one could run the optimization multiple times until convergence is reached.

BiasedWalker.m is a good example script to see how to use EntropyEst.m and reproduces some of the results from Fig 2.

The following scripts relate to the experimental figures. They will take a while to run, in particular the bootstrap results were run with a very high n for the figures, running them for lower n will get roughly, but not exactly the same results.

FlagellaMotor.m shows how to get results from Fig 3. The data files are in MTB24_10mM.mat, MTB24_85mM.mat, MTB32_85mM.mat.
Experimental data courtesy of Jasmine Nirody, and was obtained similarly to J. A. Nirody, A. L. Nord, and R. M. Berry, “Load-
dependent adaptation near zero load in the bacterial flagellar motor,” J. Royal Soc. Interface 16, 20190300 (2019).

MTBootstrap.m shows how to get results from Fig 4. The data files are in MT_traj.mat.
Experimental data courtesy of Benjamin Lacroix.

HEKOsc.m shows how to get results from Fig 4. The data files are in HEK_data.mat.
Experimental data from K. Thurley, S. C. Tovey, G. Moenke, V. L. Prince, A. Meena, A. P. Thomas, A. Skupin, C. W. Taylor, and M. Falcke, “Reliable encoding of stimulus intensities within random sequences of intracellular ca2+ spikes,” Sci. Signal. 7, ra59 (2014).

The remaining functions are helper functions and will be called by the above scripts.
