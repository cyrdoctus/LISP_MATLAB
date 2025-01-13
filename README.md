# LISP_MATLAB (Lambert Interplanetary Space Plotter)
The following code is adaptation of Lambert Problem Solver based on Howard C. Textbook into a MATLAB code that allows universally simulating transfers and plotting between two chosen planets as well as the launch and arrival dates. 

The program utilizes MATLAB language with Aerospace and Parallel Computing toolbox installed. The Aerospace Toolbox is applied through using the Ephemeris addon which allows for choice of planet and Jualian date to gain the precise position and velocity of the selected planet. This toolbox could be substituted with manual ephemeris entry or through approximations of the plaet orbits, however for higher accuracy the addon is utilized and is necessary for the functionality of the code. 

Ephemeris Addon: https://www.mathworks.com/matlabcentral/fileexchange/46671-ephemeris-data-for-aerospace-toolbox?s_tid=srchtitle_support_results_2_phemeris

Once the above addon is installed the code has segment where a user can enter a name of the desired planets in this case the launch and arrival planets and the specificed trajectory times and orbientation. For time it is meant the year/month/day/hour/minute format for launch and arrival, with trajectory orientation being determined by prograde or retrograde choice, based on the Howard Curtis textbook code for Lambert 1 Transfer.
