#include "Colvar.h"
#include "ActionRegister.h"
#include <cmath>
#include <omp.h>
#include <iostream>
#include "tools/OpenMP.h"

namespace PLMD
{
    namespace colvar
    {

        class memFusionP : public Colvar
        {
            std::vector<AtomNumber> UMEM, LMEM, TAILS;
            std::vector<double> NSMEM, DSMEM, HMEM, RCYLMEM, ZETAMEM, ONEOVERS2C2CUTOFF, XCYL, YCYL;

        public:
            explicit memFusionP(const ActionOptions &);
            void calculate() override;
            static void registerKeywords(Keywords &keys);
        };

        PLUMED_REGISTER_ACTION(memFusionP, "MEMFUSIONP")

        void memFusionP::registerKeywords(Keywords &keys)
        {
            Colvar::registerKeywords(keys);
            keys.add("atoms", "UMEMBRANE", "all the beads of the upper membrane.");
            keys.add("atoms", "LMEMBRANE", "all the beads of the lower membrane.");
            keys.add("atoms", "TAILS", "all the tail beads of the system.");
            keys.add("optional", "NSMEM", "the number of slices of the membrane fusion cylinder.");
            keys.add("optional", "DSMEM", "thickness of the slices of the membrane fusion cylinder.");
            keys.add("optional", "HMEM", "parameter of the step function θ(x,h) for the membrane fusion.");
            keys.add("optional", "RCYLMEM", "the radius of the membrane fusion cylinder.");
            keys.add("optional", "ZETAMEM", "the radius of the membrane fusion cylinder.");
            keys.add("optional", "ONEOVERS2C2CUTOFF", "cut off large values for the derivative of the atan2 function to avoid violate energy.");
            keys.add("optional", "XCYL", "X coordinate of the fixed cylinder, if not present this will be calculated.");
            keys.add("optional", "YCYL", "X coordinate of the fixed cylinder, if not present this will be calculated.");
        }

        memFusionP::memFusionP(const ActionOptions &ao) : PLUMED_COLVAR_INIT(ao)
        {
            parseAtomList("UMEMBRANE", UMEM);
            if (UMEM.size() == 0)
                error("UMEMBRANE has not any atom specified.");

            parseAtomList("LMEMBRANE", LMEM);
            if (LMEM.size() == 0)
                error("LMEMBRANE has not any atom specified.");
            
            parseAtomList("TAILS", TAILS);
            if (TAILS.size() == 0)
                error("TAILS has not any atom specified.");

            parseVector("NSMEM", NSMEM);
            if (NSMEM.size() > 1)
                error("NSMEM cannot take more than one value.");
            if (NSMEM.size() == 0)
                NSMEM.push_back(86);

            parseVector("DSMEM", DSMEM);
            if (DSMEM.size() > 1)
                error("DSMEM cannot take more than one value.");
            if (DSMEM.size() == 0)
                DSMEM.push_back(0.1);

            parseVector("HMEM", HMEM);
            if (HMEM.size() > 1)
                error("HMEM cannot take more than one value.");
            if (HMEM.size() == 0)
                HMEM.push_back(0.25);

            parseVector("RCYLMEM", RCYLMEM);
            if (RCYLMEM.size() > 1)
                error("RCYLMEM cannot take more than one value.");
            if (RCYLMEM.size() == 0)
                RCYLMEM.push_back(1.75);

            parseVector("ZETAMEM", ZETAMEM);
            if (ZETAMEM.size() > 1)
                error("ZETA cannot take more than one value.");
            if (ZETAMEM.size() == 0)
                ZETAMEM.push_back(0.5);

            parseVector("ONEOVERS2C2CUTOFF", ONEOVERS2C2CUTOFF);
            if (ONEOVERS2C2CUTOFF.size() > 1)
                error("ONEOVERS2C2CUTOFF cannot take more than one value.");
            if (ONEOVERS2C2CUTOFF.size() == 0)
                ONEOVERS2C2CUTOFF.push_back(500);

            parseVector("XCYL", XCYL);
            if (XCYL.size() > 1)
                error("XCYL cannot take more than one value.");
            if (XCYL.size() == 0)
                XCYL.push_back(-1.0);

            parseVector("YCYL", YCYL);
            if (YCYL.size() > 1)
                error("YCYL cannot take more than one value.");
            if (YCYL.size() == 0)
                YCYL.push_back(-1.0);

            checkRead();

            std::vector<AtomNumber> atoms;
            for (unsigned i = 0; i < UMEM.size(); i++)
            {
                atoms.push_back(UMEM[i]);
            }
            for (unsigned i = 0; i < LMEM.size(); i++)
            {
                atoms.push_back(LMEM[i]);
            }
            for (unsigned i = 0; i < TAILS.size(); i++)
            {
                atoms.push_back(TAILS[i]);
            }

            addValueWithDerivatives();
            setNotPeriodic();
            requestAtoms(atoms);
        }

        void memFusionP::calculate()
        {           
            /*************************
            *                        *
            *         System         *
            *                        *
            **************************/
            
            // Box dimensions.
            double Lx = getBox()[0][0], Ly = getBox()[1][1], Lz = getBox()[2][2];

            // Z center of the upper membrane (uMem) and lower membrane (lMem) for systems with PBC: https://en.wikipedia.org/wiki/Center_of_mass#Systems_with_periodic_boundary_conditions .
            double ZuMem, ZlMem, ZuMemcos = 0.0, ZuMemsin = 0.0, ZlMemcos = 0.0, ZlMemsin = 0.0;

            unsigned nt = OpenMP::getNumThreads();
            #pragma omp parallel for num_threads(nt) schedule(static)
            for(unsigned t = 0; t < nt; t++)
            {
                unsigned i_start = (UMEM.size() * t) / nt;
                unsigned i_end = (UMEM.size() * (t+1)) / nt;
                
                // Local variables.
                double uMemAngle_loc, lMemAngle_loc, ZuMemcos_loc = 0.0, ZuMemsin_loc = 0.0, ZlMemcos_loc = 0.0, ZlMemsin_loc = 0.0;

                for(unsigned i = i_start; i < i_end; i++)
                {
                    uMemAngle_loc = 2.0 * M_PI * getPbc().realToScaled(pbcDistance(Vector(0.0, 0.0, 0.0), getPosition(i)))[2];
                    lMemAngle_loc = 2.0 * M_PI * getPbc().realToScaled(pbcDistance(Vector(0.0, 0.0, 0.0), getPosition(i + UMEM.size())))[2];
                    ZuMemcos_loc += cos(uMemAngle_loc);
                    ZuMemsin_loc += sin(uMemAngle_loc);
                    ZlMemcos_loc += cos(lMemAngle_loc);
                    ZlMemsin_loc += sin(lMemAngle_loc);
                }
                
                // Shared variables.
                ZuMemcos += ZuMemcos_loc / UMEM.size();
                ZuMemsin += ZuMemsin_loc / UMEM.size();
                ZlMemcos += ZlMemcos_loc / UMEM.size();
                ZlMemsin += ZlMemsin_loc / UMEM.size();
            }

            ZuMem = Lz * (atan2(-ZuMemsin, -ZuMemcos) + M_PI) / (2.0 * M_PI);
            ZlMem = Lz * (atan2(-ZlMemsin, -ZlMemcos) + M_PI) / (2.0 * M_PI);

            // Z center of the boths membranes (upper and lower).
            double ZMems = (ZuMem + ZlMem) / 2.0;

            /*************************
            *                        *
            *         Xi_Mem         *
            *                        *
            **************************/
           
            // Quantity of beads of the membranes.
            unsigned membraneBeads = UMEM.size() + LMEM.size();

            // Z position of the first slice.
            double firstSliceZ_Mem = ZMems + (0.0 + 0.5 - NSMEM[0] / 2.0) * DSMEM[0];

            // Z distance between the first slice and the Z center of the membrane.
            double firstSliceZDist_Mem = pbcDistance(Vector(0.0, 0.0, firstSliceZ_Mem), Vector(0.0, 0.0, ZMems))[2];

            // Slices to analyze per particle.
            std::vector<unsigned> s1_Mem(TAILS.size()), s2_Mem(TAILS.size());

            // Mark the particles to analyze.
            std::vector<double> analyzeThisParticle_Mem(TAILS.size());

            // Eq. 7 Hub & Awasthi JCTC 2017.
            std::vector<double> faxial_Mem(TAILS.size() * NSMEM[0]);

            // Eq. 16 Hub & Awasthi JCTC 2017.
            std::vector<double> d_faxial_Mem_dz(TAILS.size() * NSMEM[0]);

            // Eq. 10 Hub & Awasthi JCTC 2017.
            std::vector<double> Fs_Mem(NSMEM[0]);

            // Eq. 11 Hub & Awasthi JCTC 2017.
            std::vector<double> ws_Mem(NSMEM[0]);

            // Eq. 10 Hub & Awasthi JCTC 2017.
            double W_Mem = 0.0;

            // Eq. 21 and 22 Hub & Awasthi JCTC 2017.
            std::vector<double> sx_Mem(NSMEM[0]), sy_Mem(NSMEM[0]), cx_Mem(NSMEM[0]), cy_Mem(NSMEM[0]);

            // Eq. 10 Hub & Awasthi JCTC 2017.
            double Xsc_Mem = 0.0, Xcc_Mem = 0.0, Ysc_Mem = 0.0, Ycc_Mem = 0.0;

            #pragma omp parallel for num_threads(nt) schedule(static)
            for(unsigned t = 0; t < nt; t++)
            {
                unsigned i_start = (TAILS.size() * t) / nt;
                unsigned i_end = (TAILS.size() * (t+1)) / nt;

                // Z distance from the lipid tail to the geometric center of both membranes.
                double ZTailDistance_loc;

                // Position in the cylinder.
                double PositionS_Mem_loc;

                // Scaled position of the lipid tail respect the origin of coordinates.
                Vector TailPosition_loc;

                // Aux.
                double x_loc, aux_loc;

                for (unsigned i = i_start; i < i_end; i++)
                {
                    ZTailDistance_loc = pbcDistance(Vector(0.0, 0.0, ZMems), getPosition(i + membraneBeads))[2];
                    PositionS_Mem_loc = (ZTailDistance_loc + firstSliceZDist_Mem) / DSMEM[0];
                    
                    // If the following condition is met the particle is in the Z space of the cylinder.
                    if ((PositionS_Mem_loc >= (-0.5 - HMEM[0])) && (PositionS_Mem_loc <= (NSMEM[0] + 0.5 - 1.0 + HMEM[0])))
                    {
                        analyzeThisParticle_Mem[i] = 1.0;
                        // Defining the slices to analyze each particle.
                        if (PositionS_Mem_loc < 1)
                        {
                            s1_Mem[i] = 0;
                            s2_Mem[i] = 2;
                        }
                        else if (PositionS_Mem_loc <= (NSMEM[0] - 2.0))
                        {
                            s1_Mem[i] = floor(PositionS_Mem_loc) - 1;
                            s2_Mem[i] = floor(PositionS_Mem_loc) + 1;
                        }
                        else
                        {
                            s1_Mem[i] = NSMEM[0] - 3;
                            s2_Mem[i] = NSMEM[0] - 1;
                        }

                        TailPosition_loc = getPbc().realToScaled(pbcDistance(Vector(0.0, 0.0, 0.0), getPosition(i + membraneBeads)));

                        for (unsigned s = s1_Mem[i]; s <= s2_Mem[i]; s++)
                        {
                            x_loc = (ZTailDistance_loc - (s + 0.5 - NSMEM[0] / 2.0) * DSMEM[0]) * 2.0 / DSMEM[0];
                            if (!((x_loc <= -1.0 - HMEM[0]) || (x_loc >= 1.0 + HMEM[0])))
                            {
                                if (((-1.0 + HMEM[0]) <= x_loc) && (x_loc <= (1.0 - HMEM[0])))
                                {
                                    faxial_Mem[i + TAILS.size() * s] = 1.0;
                                    Fs_Mem[s] += 1.0;
                                    sx_Mem[s] += sin(2.0 * M_PI * TailPosition_loc[0]);
                                    sy_Mem[s] += sin(2.0 * M_PI * TailPosition_loc[1]);
                                    cx_Mem[s] += cos(2.0 * M_PI * TailPosition_loc[0]);
                                    cy_Mem[s] += cos(2.0 * M_PI * TailPosition_loc[1]);
                                }
                                else if (((1.0 - HMEM[0]) < x_loc) && (x_loc < (1.0 + HMEM[0])))
                                {
                                    aux_loc = 0.5 - ((3.0 * x_loc - 3.0) / (4.0 * HMEM[0])) + (pow((x_loc - 1.0), 3) / (4.0 * pow(HMEM[0], 3)));
                                    faxial_Mem[i + TAILS.size() * s] = aux_loc;
                                    d_faxial_Mem_dz[i + TAILS.size() * s] = ((-3.0 / (4.0 * HMEM[0])) + ((3.0 * pow((x_loc - 1), 2)) / (4.0 * pow(HMEM[0], 3)))) * 2.0 / DSMEM[0];
                                    Fs_Mem[s] += aux_loc;
                                    sx_Mem[s] += aux_loc * sin(2.0 * M_PI * TailPosition_loc[0]);
                                    sy_Mem[s] += aux_loc * sin(2.0 * M_PI * TailPosition_loc[1]);
                                    cx_Mem[s] += aux_loc * cos(2.0 * M_PI * TailPosition_loc[0]);
                                    cy_Mem[s] += aux_loc * cos(2.0 * M_PI * TailPosition_loc[1]);
                                }
                                else if (((-1.0 - HMEM[0]) < x_loc) && (x_loc < (-1.0 + HMEM[0])))
                                {
                                    aux_loc = 0.5 + ((3.0 * x_loc + 3.0) / (4.0 * HMEM[0])) - (pow((x_loc + 1.0), 3) / (4.0 * pow(HMEM[0], 3)));
                                    faxial_Mem[i + TAILS.size() * s] = aux_loc;
                                    d_faxial_Mem_dz[i + TAILS.size() * s] = ((3.0 / (4.0 * HMEM[0])) - ((3.0 * pow((x_loc + 1), 2)) / (4.0 * pow(HMEM[0], 3)))) * 2.0 / DSMEM[0];
                                    Fs_Mem[s] += aux_loc;
                                    sx_Mem[s] += (aux_loc * sin(2.0 * M_PI * TailPosition_loc[0]));
                                    sy_Mem[s] += (aux_loc * sin(2.0 * M_PI * TailPosition_loc[1]));
                                    cx_Mem[s] += (aux_loc * cos(2.0 * M_PI * TailPosition_loc[0]));
                                    cy_Mem[s] += (aux_loc * cos(2.0 * M_PI * TailPosition_loc[1]));
                                }
                            }
                        }
                    }
                }
            }

            for (unsigned s = 0; s < NSMEM[0]; s++)
            {
                if (Fs_Mem[s] != 0.0)
                {
                    ws_Mem[s] = tanh(Fs_Mem[s]);
                    W_Mem += ws_Mem[s];
                    sx_Mem[s] = sx_Mem[s] / Fs_Mem[s];
                    sy_Mem[s] = sy_Mem[s] / Fs_Mem[s];
                    cx_Mem[s] = cx_Mem[s] / Fs_Mem[s];
                    cy_Mem[s] = cy_Mem[s] / Fs_Mem[s];
                    Xsc_Mem += sx_Mem[s] * ws_Mem[s];
                    Ysc_Mem += sy_Mem[s] * ws_Mem[s];
                    Xcc_Mem += cx_Mem[s] * ws_Mem[s];
                    Ycc_Mem += cy_Mem[s] * ws_Mem[s];
                }
            }

            Xsc_Mem = Xsc_Mem / W_Mem;
            Ysc_Mem = Ysc_Mem / W_Mem;
            Xcc_Mem = Xcc_Mem / W_Mem;
            Ycc_Mem = Ycc_Mem / W_Mem;

            // Eq. 12 Hub & Awasthi JCTC 2017.
            double Xcyl_Mem, Ycyl_Mem;

            if ((XCYL[0] > 0.0) && (YCYL[0] > 0.0))
            {
                Xcyl_Mem = XCYL[0];
                Ycyl_Mem = YCYL[0];
            }
            else
            {
                Xcyl_Mem = (atan2(-Xsc_Mem, -Xcc_Mem) + M_PI) * Lx / (2 * M_PI);
                Ycyl_Mem = (atan2(-Ysc_Mem, -Ycc_Mem) + M_PI) * Ly / (2 * M_PI);
            }

            //CHEQUEADO HASTA ACÁ

            // Eq. 29 Hub & Awasthi JCTC 2017.
            double d_ws_Mem_dz;

            // Center of the cylinder. XY components are calculated (or defined), Z is the Z geometric center of the membranes of the system.
            Vector xyzCyl_Mem = pbcDistance(Vector(0.0, 0.0, 0.0), Vector(Xcyl_Mem, Ycyl_Mem, ZMems));

            // Distances from the lipid tails to center of the cylinder.
            std::vector<Vector> CylDistances_Mem(TAILS.size());

            // Eq. 8 Hub & Awasthi JCTC 2017.
            double fradial_Mem;

            // Eq. 15 Hub & Awasthi JCTC 2017.
            std::vector<double> d_fradial_Mem_dx(TAILS.size()), d_fradial_Mem_dy(TAILS.size());

            // Eq. 35, 36, 37 and 38 Hub & Awasthi JCTC 2017.
            std::vector<double> d_Xcyl_Mem_dx(TAILS.size()), d_Xcyl_Mem_dz(TAILS.size()), d_Ycyl_Mem_dy(TAILS.size()), d_Ycyl_Mem_dz(TAILS.size());

            // To avoid rare instabilities auxX_Mem and auxY_Mem are truncated at a configurable value (default = 500).
            double auxX_Mem = (1 / (pow(Xsc_Mem, 2) + pow(Xcc_Mem, 2))), auxY_Mem = (1 / (pow(Ysc_Mem, 2) + pow(Ycc_Mem, 2)));

            if (auxX_Mem > ONEOVERS2C2CUTOFF[0])
            {
                auxX_Mem = Lx * ONEOVERS2C2CUTOFF[0] / (2 * M_PI);
            }
            else
            {
                auxX_Mem = Lx * auxX_Mem / (2 * M_PI);
            }

            if (auxY_Mem > ONEOVERS2C2CUTOFF[0])
            {
                auxY_Mem = Ly * ONEOVERS2C2CUTOFF[0] / (2 * M_PI);
            }
            else
            {
                auxY_Mem = Ly * auxY_Mem / (2 * M_PI);
            }

            // Number of lipid tails within the slice s of the membranes cylinder.
            std::vector<double> Nsp_Mem(NSMEM[0]), psi_Mem(NSMEM[0]), d_psi_Mem(NSMEM[0]);

            // Eq. 3 Hub & Awasthi JCTC 2017.
            double b_Mem = (ZETAMEM[0] / (1.0 - ZETAMEM[0])), c_Mem = ((1.0 - ZETAMEM[0]) * exp(b_Mem));

            // Eq. 19 Hub & Awasthi JCTC 2017.
            std::vector<double> fradial_Mem_d_faxial_Mem_dz(TAILS.size() * NSMEM[0]);

            // Eq. 20 Hub & Awasthi JCTC 2017.
            std::vector<double> Axs_Mem(NSMEM[0]), Ays_Mem(NSMEM[0]);

            // Eq. 1 Hub & Awasthi JCTC 2017. This is the CV that describes de Pore Nucleation.
            double Xi_Mem = 0.0;

            #pragma omp parallel for num_threads(nt) schedule(static)
            for(unsigned t = 0; t < nt; t++)
            {
                unsigned i_start = (TAILS.size() * t) / nt;
                unsigned i_end = (TAILS.size() * (t+1)) / nt;

                // Scaled position of the lipid tail respect the origin of coordinates.
                Vector TailPosition_loc;

                // Eq. 31, 32 and 33 Hub & Awasthi JCTC 2017
                double d_Xsc_Mem_dx_loc, d_Xsc_Mem_dz_loc, d_Xcc_Mem_dx_loc, d_Xcc_Mem_dz_loc;
                double d_Ysc_Mem_dy_loc, d_Ysc_Mem_dz_loc, d_Ycc_Mem_dy_loc, d_Ycc_Mem_dz_loc;

                // Eq. 25, 26 and 27 Hub & Awasthi JCTC 2017.
                double d_sx_Mem_dx_loc, d_sx_Mem_dz_loc, d_sy_Mem_dy_loc, d_sy_Mem_dz_loc;
                double d_cx_Mem_dx_loc, d_cx_Mem_dz_loc, d_cy_Mem_dy_loc, d_cy_Mem_dz_loc;

                // Aux.
                double x_loc;

                // XY distance from the lipid tails to the center of the cylinder.
                double ri_Mem_loc;

                for (unsigned i = i_start; i < i_end; i++)
                {
                    if (analyzeThisParticle_Mem[i])
                    {
                        TailPosition_loc = getPbc().realToScaled(pbcDistance(Vector(0.0, 0.0, 0.0), getPosition(i + membraneBeads)));
                        d_Xsc_Mem_dx_loc = 0.0;
                        d_Xcc_Mem_dx_loc = 0.0;
                        d_Ysc_Mem_dy_loc = 0.0;
                        d_Ycc_Mem_dy_loc = 0.0;
                        d_Xsc_Mem_dz_loc = 0.0;
                        d_Xcc_Mem_dz_loc = 0.0;
                        d_Ysc_Mem_dz_loc = 0.0;
                        d_Ycc_Mem_dz_loc = 0.0;

                        for (unsigned s = s1_Mem[i]; s <= s2_Mem[i]; s++)
                        {
                            if (Fs_Mem[s] != 0.0)
                            {
                                d_sx_Mem_dx_loc = faxial_Mem[i + TAILS.size() * s] * 2.0 * M_PI * cos(2.0 * M_PI * TailPosition_loc[0]) / (Lx * Fs_Mem[s]);
                                d_sy_Mem_dy_loc = faxial_Mem[i + TAILS.size() * s] * 2.0 * M_PI * cos(2.0 * M_PI * TailPosition_loc[1]) / (Ly * Fs_Mem[s]);
                                d_cx_Mem_dx_loc = -faxial_Mem[i + TAILS.size() * s] * 2.0 * M_PI * sin(2.0 * M_PI * TailPosition_loc[0]) / (Lx * Fs_Mem[s]);
                                d_cy_Mem_dy_loc = -faxial_Mem[i + TAILS.size() * s] * 2.0 * M_PI * sin(2.0 * M_PI * TailPosition_loc[1]) / (Ly * Fs_Mem[s]);
                                d_Xsc_Mem_dx_loc += ws_Mem[s] * d_sx_Mem_dx_loc / W_Mem;
                                d_Xcc_Mem_dx_loc += ws_Mem[s] * d_cx_Mem_dx_loc / W_Mem;
                                d_Ysc_Mem_dy_loc += ws_Mem[s] * d_sy_Mem_dy_loc / W_Mem;
                                d_Ycc_Mem_dy_loc += ws_Mem[s] * d_cy_Mem_dy_loc / W_Mem;

                                d_sx_Mem_dz_loc = d_faxial_Mem_dz[i + TAILS.size() * s] * (sin(2.0 * M_PI * TailPosition_loc[0]) - sx_Mem[s]) / Fs_Mem[s];
                                d_sy_Mem_dz_loc = d_faxial_Mem_dz[i + TAILS.size() * s] * (sin(2.0 * M_PI * TailPosition_loc[1]) - sy_Mem[s]) / Fs_Mem[s];
                                d_cx_Mem_dz_loc = d_faxial_Mem_dz[i + TAILS.size() * s] * (cos(2.0 * M_PI * TailPosition_loc[0]) - cx_Mem[s]) / Fs_Mem[s];
                                d_cy_Mem_dz_loc = d_faxial_Mem_dz[i + TAILS.size() * s] * (cos(2.0 * M_PI * TailPosition_loc[1]) - cy_Mem[s]) / Fs_Mem[s];
                                d_ws_Mem_dz = (1 - pow(ws_Mem[s], 2)) * d_faxial_Mem_dz[i + TAILS.size() * s];
                                d_Xsc_Mem_dz_loc += (ws_Mem[s] * d_sx_Mem_dz_loc + d_ws_Mem_dz * (sx_Mem[s] - Xsc_Mem)) / W_Mem;
                                d_Xcc_Mem_dz_loc += (ws_Mem[s] * d_cx_Mem_dz_loc + d_ws_Mem_dz * (cx_Mem[s] - Xcc_Mem)) / W_Mem;
                                d_Ysc_Mem_dz_loc += (ws_Mem[s] * d_sy_Mem_dz_loc + d_ws_Mem_dz * (sy_Mem[s] - Ysc_Mem)) / W_Mem;
                                d_Ycc_Mem_dz_loc += (ws_Mem[s] * d_cy_Mem_dz_loc + d_ws_Mem_dz * (cy_Mem[s] - Ycc_Mem)) / W_Mem;
                            }
                        }
                        d_Xcyl_Mem_dx[i] = auxX_Mem * (-Xsc_Mem * d_Xcc_Mem_dx_loc + Xcc_Mem * d_Xsc_Mem_dx_loc);
                        d_Xcyl_Mem_dz[i] = auxX_Mem * (-Xsc_Mem * d_Xcc_Mem_dz_loc + Xcc_Mem * d_Xsc_Mem_dz_loc);
                        d_Ycyl_Mem_dy[i] = auxY_Mem * (-Ysc_Mem * d_Ycc_Mem_dy_loc + Ycc_Mem * d_Ysc_Mem_dy_loc);
                        d_Ycyl_Mem_dz[i] = auxY_Mem * (-Ysc_Mem * d_Ycc_Mem_dz_loc + Ycc_Mem * d_Ysc_Mem_dz_loc);

                        CylDistances_Mem[i] = pbcDistance(xyzCyl_Mem, pbcDistance(Vector(0.0, 0.0, 0.0), getPosition(i + membraneBeads)));
                        ri_Mem_loc = sqrt(pow(CylDistances_Mem[i][0], 2) + pow(CylDistances_Mem[i][1], 2));
                        x_loc = ri_Mem_loc / RCYLMEM[0];
                        if (!((x_loc <= -1.0 - HMEM[0]) || (x_loc >= 1.0 + HMEM[0])))
                        {
                            if (((-1.0 + HMEM[0]) <= x_loc) && (x_loc <= (1.0 - HMEM[0])))
                            {
                                fradial_Mem = 1.0;
                            }
                            else if (((1.0 - HMEM[0]) < x_loc) && (x_loc < (1.0 + HMEM[0])))
                            {
                                fradial_Mem = 0.5 - ((3.0 * x_loc - 3.0) / (4.0 * HMEM[0])) + (pow((x_loc - 1.0), 3) / (4.0 * pow(HMEM[0], 3)));
                                d_fradial_Mem_dx[i] = ((-3.0 / (4.0 * HMEM[0])) + ((3.0 * pow((x_loc - 1), 2)) / (4.0 * pow(HMEM[0], 3)))) * CylDistances_Mem[i][0] / (RCYLMEM[0] * ri_Mem_loc);
                                d_fradial_Mem_dy[i] = ((-3.0 / (4.0 * HMEM[0])) + ((3.0 * pow((x_loc - 1), 2)) / (4.0 * pow(HMEM[0], 3)))) * CylDistances_Mem[i][1] / (RCYLMEM[0] * ri_Mem_loc);
                            }
                            else if (((-1.0 - HMEM[0]) < x_loc) && (x_loc < (-1.0 + HMEM[0])))
                            {
                                fradial_Mem = 0.5 + ((3.0 * x_loc + 3.0) / (4.0 * HMEM[0])) - (pow((x_loc + 1.0), 3) / (4.0 * pow(HMEM[0], 3)));
                                d_fradial_Mem_dx[i] = ((3.0 / (4.0 * HMEM[0])) - ((3.0 * pow((x_loc + 1), 2)) / (4.0 * pow(HMEM[0], 3)))) * CylDistances_Mem[i][0] / (RCYLMEM[0] * ri_Mem_loc);
                                d_fradial_Mem_dy[i] = ((3.0 / (4.0 * HMEM[0])) - ((3.0 * pow((x_loc + 1), 2)) / (4.0 * pow(HMEM[0], 3)))) * CylDistances_Mem[i][1] / (RCYLMEM[0] * ri_Mem_loc);
                            }

                            for (unsigned s = s1_Mem[i]; s <= s2_Mem[i]; s++)
                            {
                                Nsp_Mem[s] += fradial_Mem * faxial_Mem[i + TAILS.size() * s];
                                Axs_Mem[s] += faxial_Mem[i + TAILS.size() * s] * d_fradial_Mem_dx[i];
                                Ays_Mem[s] += faxial_Mem[i + TAILS.size() * s] * d_fradial_Mem_dy[i];
                                fradial_Mem_d_faxial_Mem_dz[i + TAILS.size() * s] = fradial_Mem * d_faxial_Mem_dz[i + TAILS.size() * s];
                            }                        
                        }
                    }
                }
            }

            for (unsigned s = 0; s < NSMEM[0]; s++)
            {
                if (Nsp_Mem[s] <= 1.0)
                {
                    psi_Mem[s] = ZETAMEM[0] * Nsp_Mem[s];
                    d_psi_Mem[s] = ZETAMEM[0];
                    Xi_Mem += psi_Mem[s];
                }
                else
                {
                    psi_Mem[s] = 1.0 - c_Mem * exp(-b_Mem * Nsp_Mem[s]);
                    d_psi_Mem[s] = b_Mem * c_Mem * exp(-b_Mem * Nsp_Mem[s]);
                    Xi_Mem += psi_Mem[s];
                }
            }

            Xi_Mem = Xi_Mem / NSMEM[0];

            // Eq. 18 Hub & Awasthi JCTC 2017.
            std::vector<double> faxial_Mem_d_fradial_Mem_dx(TAILS.size() * NSMEM[0]), faxial_Mem_d_fradial_Mem_dy(TAILS.size() * NSMEM[0]), faxial_Mem_d_fradial_Mem_dz(TAILS.size() * NSMEM[0]);

            // Eq. 13 Hub & Awasthi JCTC 2017.
            std::vector<Vector> derivatives_Mem(TAILS.size());


            #pragma omp parallel for num_threads(nt) schedule(static)
            for(unsigned t = 0; t < nt; t++)
            {
                unsigned i_start = (TAILS.size() * t) / nt;
                unsigned i_end = (TAILS.size() * (t+1)) / nt;

                // Aux.
                double aux_loc;

                for (unsigned i = i_start; i < i_end; i++)
                {
                    if (analyzeThisParticle_Mem[i])
                    {
                        for (unsigned s = s1_Mem[i]; s <= s2_Mem[i]; s++)
                        {
                            if (faxial_Mem[i + TAILS.size() * s])
                            {
                                faxial_Mem_d_fradial_Mem_dx[i + TAILS.size() * s] = faxial_Mem[i + TAILS.size() * s] * d_fradial_Mem_dx[i] - d_Xcyl_Mem_dx[i] * Axs_Mem[s];
                                faxial_Mem_d_fradial_Mem_dy[i + TAILS.size() * s] = faxial_Mem[i + TAILS.size() * s] * d_fradial_Mem_dy[i] - d_Ycyl_Mem_dy[i] * Ays_Mem[s];
                                faxial_Mem_d_fradial_Mem_dz[i + TAILS.size() * s] = -d_Xcyl_Mem_dz[i] * Axs_Mem[s] - d_Ycyl_Mem_dz[i] * Ays_Mem[s];
                            }
                        }

                        for (unsigned s = s1_Mem[i]; s <= s2_Mem[i]; s++)
                        {
                            aux_loc = d_psi_Mem[s] / NSMEM[0];
                            derivatives_Mem[i][0] += aux_loc * faxial_Mem_d_fradial_Mem_dx[i + TAILS.size() * s];
                            derivatives_Mem[i][1] += aux_loc * faxial_Mem_d_fradial_Mem_dy[i + TAILS.size() * s];
                            derivatives_Mem[i][2] += aux_loc * (faxial_Mem_d_fradial_Mem_dz[i + TAILS.size() * s] + fradial_Mem_d_faxial_Mem_dz[i + TAILS.size() * s]);
                        }
                    }
                }
            }

            // Derivatives and virial for the Xi_Mem.
            Tensor virial;
            for (unsigned i = 0; i < TAILS.size(); i++)
            {
                setAtomsDerivatives((i + membraneBeads), derivatives_Mem[i]);
                virial -= Tensor(CylDistances_Mem[i], derivatives_Mem[i]);
            }

            setValue(Xi_Mem);
            setBoxDerivatives(virial);
        }
    }
}
