# MemFusion with PLUMED CV

Installation:
- Download and uncompress [PLUMED](https://www.plumed.org/).
- Download **MemFusion.cpp** and copy it in **/plumed-2.7.2/src/colvar**.
- Install PLUMED.
- Download and install [GROMACS](https://manual.gromacs.org/documentation/).
- Patch GROMACS with PLUMED.
- Reinstall the patched GROMACS.

Example folder contains all files needed for MD to fuse bilayers with one Syt1-C2B domain.

Update 1: If you prefer the new parallel version of the CV, download **MemFusionP.cpp** instead of **MemFusion.cpp** and follow the same process as above. An example folder for the parallel release of the CV is available.
