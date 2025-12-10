# Notes for the simulation codes
In this project i will attempt to make a functiong and physically correct simulation to compare my experimental results to.

## Landau–Lifshitz–Bloch, Two Temperature Model
This simulation is almost completely based on the article Raposo (2020) Micromagnetic Modeling of All Optical Switching of Ferromagnetic Thin Films: The Role of Inverse Faraday Effect and Magnetic Circular Dichroism. It is an attempt to see wether or not i can recreate an already existing simulation using the articles equations and logic.

**Wed Dec 3:**
For recreating the simulation i first created a laser pulse class so that i have a object oreinted systme that makes the simulation easier to understand. Then i expanded it to a object for a sequence of lasers consisting of several laserpulses, this is customizable and thus allows for easy control. After that i added a ferromagnetic material to compute the IFE response to the laserpulse. I also added the two temperature model and used it to simulate how the electron and lattice temperatures change. Finally i made a animation showing the time evolution for the laser power, IFE field, electron temperature and the lattice temperature and automatically make a animations folder in the correct place and save a .mp4 file there with the final animation.

**Fri Dec 5:**
Attempt at adding the Landau-Lifshitz-Bloch (LLB) model to allow for actual switching of the mz component, thus far i only magnagezed to demagnetize it not te switch it. I don't think the implementation is necceceraly worng but it might be that more thorough investigation is required and that i should just try more parameter combinations to see if i can get some switching to happen. It might be that something in the simulation framwork is wrong, i have to investigate thet.

**Sat Dec 6:**
Tried to fix the model_llb solver but i think i just made things worse. I might need to just delete it and redo everything properly this time. I don't know why suddenly the simulation goes crazy but it is annoying me. Volgende keer alle onderdelen van H_effective controleren en kijken welke de problemen veroorzaakt, ik ben bijna zeker dat het H_m is maar ik wil het zeker weten.