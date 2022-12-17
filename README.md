# Numerical Simulation of an Antihydrogen atom in a magnetic trap
This repo contains the source code to perform the numerical integration of the equation of motion of an antihydrogen atom in a magentic trap.
The spatial degrees of motion of the atom are treated classically. Given that this is a neutral particle this trajectory is just rectilinear motion at a fixed speed. The spin degree of freedom, which in this case is that of a spin 1 particle, is treated semiclasically by treating the magnetic field of the trap as a classical field. 
The atom is initialized in one of the three possible spin states and we compute the probability for the spin of the particle to flip as it crosses the magnetic trap, specially near the points where the field goes to zero where there is a possibility for a non adiabatic transition to take place. 
The main file is 
