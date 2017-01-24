# hydropad
Hydropad numerical code developed for studying the formation and evolution of cosmological structures in both baryonic 
and dark components. Collisional matter is treated as a fluid and the corresponding hydrodynamic equations are solved using the PPM scheme on a fixed Eulerian grid. We have described the changes to the basic method required by the cosmological applications. Particular care has been taken in including expansion and gravity in the Riemann solver and in the final integration step. This has required the calculation of the characteristic form of the hydrodynamic equations in expanding coordinates. A double formulation of the energy equation has allowed a proper treatment of the highly supersonicflows common in cosmological simulations. 
The behaviour of the dark matter, the gravitational field and other physical processes are not implemented yet.
Great emphasis is given to performance and scalability to run efficiently on hybrid accelerated HPC systems.
At the moment the code supports only a regular, single resolution mesh, but it is designe to easily support Structured AMR. 
