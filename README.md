This is the source code for the deformable particle model (both adaptive and fixed p0 versions).
The code is written in C++. To create the executable, run:

_g++ -o main main.cpp_

To run the executable in the same directory:

./main $\epsilon_a$ eb el B Nt f0 tau tau_p

The inputs to the script are the following parameters:

ea = cell area conservation energy

eb = cell membrane bending energy

el = cell perimeter elastic energy

B = drag constant

Nt = number of timesteps for simulation

f0 = active brownian force

tau = persistence time for active brownian force

tau_p = perimeter adaptation timescale

