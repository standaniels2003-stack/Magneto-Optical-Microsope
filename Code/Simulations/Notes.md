# Notes for the simulation codes
In this project i will attempt to make a functiong and physically correct simulation to compare my experimental results to.

## Landau–Lifshitz–Bloch, Two Temperature Model
This simulation is almost completely based on the article Raposo (2020) Micromagnetic Modeling of All Optical Switching of Ferromagnetic Thin Films: The Role of Inverse Faraday Effect and Magnetic Circular Dichroism. It is an attempt to see wether or not i can recreate an already existing simulation using the articles equations and logic.

**Wed Dec 3:**
For recreating the simulation i first created a laser pulse class so that i have a object oreinted systme that makes the simulation easier to understand. Then i expanded it to a object for a sequence of lasers consisting of several laserpulses, this is customizable and thus allows for easy control. After that i added a ferromagnetic material to compute the IFE response to the laserpulse. I also added the two temperature model and used it to simulate how the electron and lattice temperatures change. Finally i made a animation showing the time evolution for the laser power, IFE field, electron temperature and the lattice temperature and automatically make a animations folder in the correct place and save a .mp4 file there with the final animation.