This is the source code for the deformable particle model (both adaptive and fixed $p_0$ versions).
The code is written in C++. To create the executable, run:

_g++ -o main main.cpp_

To run the executable in the same directory:

./main $\epsilon_a$ $\epsilon_b$ $\epsilon_l$ $\beta$ $Nt$ $f_0$ $\tau$ $\tau_p$

The inputs to the script are the following parameters:

$\epsilon_a$ = cell area conservation energy

$\epsilon_b$ = cell membrane bending energy

$\epsilon_l$ = cell perimeter elastic energy

$\beta$ = drag constant

$Nt$ = number of timesteps for simulation

$f_0$ = active brownian force

$\tau$ = persistence time for active brownian force

$\tau_p$ = perimeter adaptation timescale

