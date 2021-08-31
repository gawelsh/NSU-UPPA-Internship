# NSU-UPPA-Internship
Project by Gaël Dufau (UPPA) with Prof. Andrey Marchuk (NSU, Tsunami Lab)

This code simulates the behaviour of a long wave with linear and nonlinear shallow water equations and contains the followings:

-longwaves.sh : shell  lines  that  execute  the  Makefile  and  launch  simulations (From  a  Linuxterminal with the command./longwaves.sh). If you want to simulate a wave thank to linear equation, please add "./include/longwaves.exe" or "./include/longwavesCN.exe". If you want to simulate a wave thank to nonlinear equation, please add "./include/longwavesNL.exe".

-common.c/common.h:  C program that contains initialisation information needed to enhance the program.  You can modify depth function and initialisation function incommon.cand theparameters in common.h.

-python programs: examples of interpretation that can be drawn with the simulated .dat files.

Best regards.
