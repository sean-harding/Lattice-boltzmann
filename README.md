# Lattice-boltzmann
-------------------|SIMULATION OF NON-EQUILIBRIUM FLUID FLOW USING THE LATTICE BOLTZMANN METHOD-------------------|
-> Algorithm details are provided in "Graduate texts in physics: The Lattice Boltzmann Method - Principles and Practice"
-> The basic idea is to work with the DISTRIBUTION FUNCTION f(r,t,v) of the fluid, which gives the probability of a fluid
molecule having a particular velocity, v at a position r and time t.
-> The equation which determines the evolution of the distribution function is the BOLTZMANN EQUATION
-> By discretizing the fluid into boxes/cubes at discrete positions, the distribution function is discretized as
    f(r,t,v)-> f_i(r_k,t), where r_k points to a discrete lattice point. Each f_i goes along with a VELOCITY v_i
    which points between the lattice site r_k and some other lattice site at r_k + t*v_i
-> By deciding on a TRUNCATED SET of velocities, the distribution function is represented simply as a set of numbers.
-> Which velocity set is used depends on the calculation, and how the time interval t= n*dT for integer n and double dT is
chosen

In this simulation, we are working with a fluid in 2D, with 9 discrete velocities. For a given lattice site, the velocity 
vector in position [i] in the array points to adjacent sites. The order of these sites is:
        6 2 5
        3 0 1
        7 4 8
for example, the magnitude of the velocity vector (-1,0) relative to a lattice site is stored in array index [3].

-> The LATTICE BOLTZMANN equation is the discrete version of the Boltzmann equation. It allows us to update the distribution
function in time. This is done in two steps:
1.) At each site at a given time, the system is in a non-equilibrium distribution. If particles are constrained so that
they cannot move between lattice sites, then this distribution is assumed to relax to equilibrium at some characteristic
rate, due to collisions between molecules in the fluid. In practice, the timestep of our simulation is smaller than this
relaxation time, so a COLLISION function has been written to take into account partial relaxation of the distribution
2.) Particles with a given distribution are STREAMED to neighboring sites, representing the advection of matter.

-> The LB algorithm consists of successive COLLISION and STREAMING steps to simulate the dynamics of the fluid.
