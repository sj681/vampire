//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Andrea Meo and Richard F L Evans 2016. All rights reserved.
//
//------------------------------------------------------------------------------
//

// C++ standard library headers
#include <cmath>
#include <cstdlib>
#include <iostream>

// C library headers
#include <fenv.h>
#include <signal.h>

// Vampire headers
#include "cells.hpp" // needed for cells::cell_id_array but to be removed
#include "dipole.hpp"
#include "vio.hpp"
#include "vutil.hpp"
#include "math.h"
#include "errors.hpp"

// dipole module headers
#include "internal.hpp"

namespace dp = dipole::internal;


int largest(int x, int y, int z){
    int largest = x;
    if (y > largest) largest = y;
    if (z > largest) largest = z;

   return largest;
}


namespace dipole{

   namespace internal{


      void calcualte_cells_levels(const double system_dimensions_x,
                      const double system_dimensions_y,
                      const double system_dimensions_z,
                      const double unit_cell_size_x,
                      const double unit_cell_size_y,
                      const double unit_cell_size_z,
                      const std::vector<double>& atom_coords_x,
                      const std::vector<double>& atom_coords_y,
                      const std::vector<double>& atom_coords_z,
                      const std::vector<int>& atom_type_array,
                      const std::vector<int>& atom_cell_id_array,
                      const int num_total_atoms_for_dipole,
                      const int num_atoms
      ){


      //    // For MPI version, only add local atoms
      //     #ifdef MPICF
      //        int num_local_atoms = vmpi::num_core_atoms+vmpi::num_bdry_atoms;
      //     #else
      //        int num_local_atoms = num_atoms;
      //     #endif
      //
      //     // Determine number of total atoms
      //     #ifdef MPICF
      //      int num_total_atoms=0;
      //      int total_non_mag_removed_atoms=0;
      //      MPI_Reduce(&num_local_atoms,&num_total_atoms, 1,MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
      //      MPI_Reduce(&create::num_total_atoms_non_filler,&total_non_mag_removed_atoms, 1, MPI_INT, MPI_SUM, 0,MPI_COMM_WORLD);
      //      int total_atoms_non_filler = num_total_atoms + total_non_mag_removed_atoms;
      //      MPI_Bcast(&total_atoms_non_filler,1,MPI_INT,0,MPI_COMM_WORLD);
      //     #else
      //      int total_atoms_non_filler = atoms::num_atoms+create::num_total_atoms_non_filler;
      //     #endif
      //
      //
      //    int largest_dimension = largest(system_dimensions_x, system_dimensions_y,system_dimensions_z);
      //
      //    int num_levels = ceil(std::log(2*largest_dimension/cells::macro_cell_size)/std::log(2.0));
      //
      //
      //    //resize the dipole arrays to store the hierarchial cell information.
      //
      //    dp::num_cells_level.resize(num_levels);
      //    dp::num_local_cells_level.resize(num_levels);
      //    dp::num_atoms_magnetic_level.resize(num_levels);
      //
      //
      //    dp::cell_position_array_level.resize(num_levels);
      //    dp::pos_and_mom_array_level.resize(num_levels);
      //
      //    dp::mag_array_x_level.resize(num_levels);
      //    dp::mag_array_y_level.resize(num_levels);
      //    dp::mag_array_z_level.resize(num_levels);
      //
      //    dp::field_array_x_level.resize(num_levels);
      //    dp::field_array_y_level.resize(num_levels);
      //    dp::field_array_z_level.resize(num_levels);
      //
      //    dp::num_atoms_in_cell_level.resize(num_levels);
      //    dp::num_atoms_in_cell_global_level.resize(num_levels);
      //    dp::volume_array_level.resize(num_levels);
      //
      //    dp::total_moment_array_level.resize(num_levels);
      //
      //    dp::total_moment_array_level.resize(num_levels);
      //
      //    dp::local_cell_array_level.resize(num_levels);
      //
      //    for (int level = 0; level < num_levels; level ++){
      //
      //      int cell_size = pow(2,level)*cells::macro_cell_size;
      //
      //      // Calculate number of microcells
      //      // determine number of stacks in x and y (global)
      //      int ncx = static_cast<unsigned int>(ceil((system_dimensions_x+0.01)/cell_size));
      //      int ncy = static_cast<unsigned int>(ceil((system_dimensions_y+0.01)/cell_size));
      //      int ncz = static_cast<unsigned int>(ceil((system_dimensions_z+0.01)/cell_size));
      //
      //      int temp_num_cells = ncx*ncy*ncz;
      //      dp::num_cells_level[level] = temp_num_cells;
      //
      //      dp::cell_position_array_level[level].resize(temp_num_cells*3);
      //
      //      //---------------------------------------------------
      //      // Determine which atoms belong to which cell
      //      //---------------------------------------------------
      //
      //
      //      // Set cell and stack counters
      //      int cell=0;
      //
      //      // Allocate space for 3D supercell array (ST coordinate system)
      //      std::vector<std::vector<std::vector<int> > > supercell_array;
      //
      //      supercell_array.resize(ncx);
      //      for(int i=0;i<ncx;++i){
      //         supercell_array[i].resize(ncy);
      //         for(int j=0;j<ncy;++j){
      //            supercell_array[i][j].resize(ncz);
      //            // store cell coordinates
      //            for(int k=0; k<ncz; ++k){
      //               // associate cell with position i,j,k
      //               supercell_array[i][j][k]=cell;
      //               // increment cell number
      //               cell++;
      //            }
      //         }
      //      }
      //
      //      // Determine number of cells in x,y,z
      //      const int d[3]={ncx,ncy,ncz};
      //      const int cs[3] = {cell_size, cell_size, cell_size}; // cell size
      //
      //      dp::atom_cell_id_array_level[level].resize(num_local_atoms,0.0);
      //      // Assign atoms to cells
      //      for(int atom=0;atom<num_local_atoms;atom++){
      //         // temporary for atom coordinates
      //         double c[3];
      //         // convert atom coordinates to st reference frame
      //         c[0]=atom_coords_x[atom]+0.0001;
      //         c[1]=atom_coords_y[atom]+0.0001;
      //         c[2]=atom_coords_z[atom]+0.0001;
      //         int scc[3]={0,0,0}; // super cell coordinates
      //         // Determine supercell coordinates for atom (rounding down)
      //         scc[0]=int(c[0]/cs[0]);
      //         scc[1]=int(c[1]/cs[1]);
      //         scc[2]=int(c[2]/cs[2]);
      //         for(int i=0;i<3;i++){
      //            // Always check cell in range
      //            if(scc[i]<0 || scc[i]>= d[i]){
      //               terminaltextcolor(RED);
      //               std::cerr << "Error - atom out of supercell range in cell calculation!" << std::endl;
      //               terminaltextcolor(WHITE);
      //               #ifdef MPICF
      //               terminaltextcolor(RED);
      //               std::cerr << "\tCPU Rank: " << vmpi::my_rank << std::endl;
      //               terminaltextcolor(WHITE);
      //               #endif
      //               terminaltextcolor(RED);
      //               std::cerr << "\tAtom number:      " << atom << std::endl;
      //               std::cerr << "\tCell size:        " << cs[0] << "\t" << cs[1] << "\t" << cs[2] << "\t" << std::endl;
      //               std::cerr << "\tAtom coordinates: " << c[0] << "\t" << c[1] << "\t" << c[2] << "\t" << std::endl;
      //               std::cerr << "\tReal coordinates: " << atom_coords_x[atom] << "\t" << atom_coords_y[atom] << "\t" << atom_coords_z[atom] << "\t" << std::endl;
      //               std::cerr << "\tCell coordinates: " << scc[0] << "\t" << scc[1] << "\t" << scc[2] << "\t" << std::endl;
      //               std::cerr << "\tCell maxima:      " << d[0] << "\t" << d[1] << "\t" << d[2] << std::endl;
      //               terminaltextcolor(WHITE);
      //               err::vexit();
      //            }
      //         }
      //         // If no error for range then assign atom to cell
      //         dp::atom_cell_id_array_level[level][atom] = supercell_array[scc[0]][scc[1]][scc[2]];
      //      }
      // //      //-------------------------------------------------------------------------------------
      // //      // Determine number of microcells computed locally
      // //      //-------------------------------------------------------------------------------------
      //   //
      // //      // Resize new cell arrays
      //   //
      // //      dp::pos_and_mom_array_level[level].resize(4*temp_num_cells,0.0);
      //   //
      // //      dp::mag_array_x_level[level].resize(temp_num_cells,0.0);
      // //      dp::mag_array_y_level[level].resize(temp_num_cells,0.0);
      // //      dp::mag_array_z_level[level].resize(temp_num_cells,0.0);
      //   //
      // //      dp::field_array_x_level[level].resize(temp_num_cells,0.0);
      // //      dp::field_array_y_level[level].resize(temp_num_cells,0.0);
      // //      dp::field_array_z_level[level].resize(temp_num_cells,0.0);
      //   //
      // //      dp::num_atoms_in_cell_level[level].resize(temp_num_cells,0);
      // //      dp::num_atoms_in_cell_global_level[level].resize(0);
      // //      dp::volume_array_level[level].resize(temp_num_cells,0.0);
      //   //
      // //      dp::total_moment_array_level[level].resize(temp_num_cells,0.0);
      //   //
      // //      std::vector < double > temp_num_atoms_in_cell(temp_num_cells,0.0);
      // //      std::vector < double > temp_pos_and_mom_array(4*temp_num_cells,0.0);
      // //      std::vector < double >  temp_num_atoms_in_cell_global;
      // //      temp_num_atoms_in_cell_global.resize(0);
      //   //
      // //      int temp_num_atoms_magnetic = 0;
      // //      // Now add atoms to each cell as magnetic 'centre of mass'
      // //      for(int atom=0;atom<num_local_atoms;atom++){
      // //         int local_cell=dp::atom_cell_id_array_level[level][atom];
      // //         //int type = cells::internal::atom_type_array[atom];
      // //         int type = atom_type_array[atom];
      // //         const double mus = mp::material[type].mu_s_SI;
      // //         // Consider only magnetic elements
      // //         if(mp::material[type].non_magnetic==0){
      //   //
      // //            temp_pos_and_mom_array[4*local_cell+0] += atom_coords_x[atom]*mus;
      // //            temp_pos_and_mom_array[4*local_cell+1] += atom_coords_y[atom]*mus;
      // //            temp_pos_and_mom_array[4*local_cell+2] += atom_coords_z[atom]*mus;
      // //            temp_pos_and_mom_array[4*local_cell+3] += mus;
      //   //
      // //            temp_num_atoms_in_cell[local_cell]++;
      // //            temp_num_atoms_magnetic++;
      // //         }
      // //      }
      //   //
      // //      // Save local cells index
      // //      for(int local_cell=0;local_cell<temp_num_cells;local_cell++){
      // //        if(dp::num_atoms_in_cell_level[level][local_cell]>0){
      // //            // add index of cell only if there are atoms inside
      // //            dp::cell_id_array_level[level].push_back(local_cell);
      // //        }
      // //      }
      //   //
      //   //
      //   //
      // //      #ifdef MPICF
      // //        MPI_Allreduce(MPI_IN_PLACE, &temp_num_atoms_in_cell[0],     temp_num_atoms_in_cell.size(),    MPI_INT,    MPI_SUM, MPI_COMM_WORLD);
      // //        MPI_Allreduce(MPI_IN_PLACE, &temp_pos_and_mom_array[0],     temp_pos_and_mom_array.size(),    MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      // //        temp_num_atoms_in_cell_global.resize(temp_num_cells);
      // //        temp_num_atoms_in_cell_global = temp_num_atoms_in_cell;
      // //        MPI_Allreduce(MPI_IN_PLACE, &temp_num_atoms_magnetic, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
      // //      #else
      // //        // copy num_atoms_in_cell to global version
      // //        temp_num_atoms_in_cell_global = temp_num_cells;
      // //      #endif
      //   //
      //   //
      // //      const double factor_for_volume = double(total_atoms_non_filler)/double(temp_num_atoms_magnetic);
      // //      const double atomic_volume =  factor_for_volume * unit_cell_size_x*unit_cell_size_y*unit_cell_size_z/double(cells::num_atoms_in_unit_cell);
      //   //
      // //      // Now find mean coordinates via magnetic 'centre of mass'
      // //      for(int local_cell=0;local_cell<temp_num_cells;local_cell++){
      // //         if(temp_num_atoms_in_cell[local_cell]>0){
      //   //
      // //            dp::pos_and_mom_array_level[level][4*local_cell+0] = temp_pos_and_mom_array[4*local_cell+0]/(temp_pos_and_mom_array[4*local_cell+3]);
      // //            dp::pos_and_mom_array_level[level][4*local_cell+1] = temp_pos_and_mom_array[4*local_cell+1]/(temp_pos_and_mom_array[4*local_cell+3]);
      // //            dp::pos_and_mom_array_level[level][4*local_cell+2] = temp_pos_and_mom_array[4*local_cell+2]/(temp_pos_and_mom_array[4*local_cell+3]);
      // //            dp::pos_and_mom_array_level[level][4*local_cell+3] = temp_pos_and_mom_array[4*local_cell+3];
      //   //
      // //            dp::volume_array_level[level][local_cell] = double(temp_num_atoms_in_cell[local_cell])*atomic_volume;
      // //            dp::num_atoms_in_cell_level[level][local_cell] = temp_num_atoms_in_cell[local_cell];
      // //            dp::num_atoms_magnetic_level[level] = temp_num_atoms_magnetic;
      // //         }
      // //      }
      //   //
      // //      // Calculate number of local cells
      // //      for(int cell=0;cell<temp_num_cells;cell++){
      // //         if(dp::num_atoms_in_cell_level[level][cell]!=0){
      // //            dp::local_cell_array_level[level].push_back(cell);
      // //            dp::num_local_cells_level[level]++;
      // //         }
      // //      }
      //   //
        //
        //
    //     } //end of loop over levels
     } //end of calculate cells function


   } //close internal namespace
} //close dipole namespace
