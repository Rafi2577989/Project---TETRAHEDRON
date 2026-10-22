# Project---TETRAHEDRON
NUMERICAL SIMULATION OF A SPRING-CONNECTED TETRAHEDRON IN A THREE-DIMENSIONAL BOX
My project simulates a deformable tetrahedron made of four equal-mass particles connected by six harmonic springs.  
Initially, all springs are at their equilibrium length and all particles have the same velocity, so the tetrahedron starts in pure translational motion.  
I studied two cases. In the first case, the tetrahedron collides with the walls of a three-dimensional box. In the second case, it collides with a rigid sphere placed at the centre of the box.  
The main goal was to study how centre-of-mass kinetic energy is transformed into internal vibration and rotation after collisions, while checking that the total mechanical energy remains conserved.  
I implemented the simulation in C++ using the Velocity Verlet method and used Python for plots and animations.
