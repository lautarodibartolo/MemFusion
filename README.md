# MemFusion with PLUMED CV

Installation:
- Download and uncompress [PLUMED](https://www.plumed.org/).
- Download **MemFusion.cpp** and copy it in **/plumed-2.7.2/src/colvar**.
- Install PLUMED.
- Download and install [GROMACS](https://manual.gromacs.org/documentation/).
- Patch GROMACS with PLUMED.
- Reinstall the patched GROMACS.

Example folder contains all files needed for MD to fuse bilayers with one Syt1-C2B domain.

## Cites
If you use this code please cite:  
- Synaptotagmin-1 C2B domains cooperatively stabilize the fusion stalk via a master-servant mechanism, Ary Lautaro Di Bartolo and Diego Masone. https://doi.org/10.1039/D1SC06711G 
- Probing a Continuous Polar Defect: A Reaction Coordinate for Pore Formation in Lipid Membranes, Jochen S. Hub and Neha Awasthi.  https://doi.org/10.1021/acs.jctc.7b00106. The implementation of the CV directly to GROMACS done by Hub's group can be found here: https://gitlab.com/cbjh/gromacs-chain-coordinate

## Updates
* If you prefer the new parallel version of the CV, download **MemFusionP.cpp** instead of **MemFusion.cpp** and follow the same process as above. An example folder for the parallel release of the CV is available.

* Install [PLUMED 2 Development Version](https://github.com/plumed/plumed2.git) with membranefusion module enabled. You can activate it at configure time using the keyword --enable-modules=membranefusion. This will install PLUMED directly with the parallel CV enabled.
