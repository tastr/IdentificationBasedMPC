Use main.m to run multiple simulations. 
Every 100 simulations (of one type, aka DD, LSQ, MDL) are saved in a .mat file. To adjust the name and path of the files use <fpath>, <fname_base>. (Note: if a file already exists, it will not be overwritten but adds an _x at the end. Somewhere it works sometimes wrong, and the file names become <file_name.mat_x>. Just rename the file to .mat, for fast fix, if needed. Did not bother me, so sorry for the bug). 
The parameters lambda_list, lambda_alpha_list, lambda_sigma_list, N_list, L_list can be used for multiple simulations (will be run in a loop, saving every 100 simulations in a .mat file)


Used for plots data: the files 
DDsimulations_fixedalpha 
DDsimulations_fixedsigma
LSQ_simulations

Only some parameters are adjusted in the scripts, the simulation is then running in simulation_main script. 

The matrix-plots are contained in 
plot_ddMatrix
plot_lsqMatrix
