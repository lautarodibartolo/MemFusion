![Logo](Logo.png)

# MemFusion

To perform this example you will need to:
- Download [Plumed](https://www.plumed.org/).
- Download the **MemFusion.cpp** from this repository and put it in the **/plumed-2.7.2/src/colvar** folder.
- Install Plumed.
- Download and install [Gromacs](https://manual.gromacs.org/documentation/).
- Patch Gromacs with Plumed.
- Reinstall the patched Gromacs.

After that, you can download the **Example** folder from this repository and execute the **Production.sh** script. It will perform a 10ns biased simulation that takes a system of two flat and parallel bilayers with a Synaptotagmine-1 C2B Domain in the cytosolic region and fuse the membranes.
