//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Sarah Jenkins and Richard F L Evans 2016. All rights reserved.
//
//   Email: sj681@york.ac.uk
//
//------------------------------------------------------------------------------
//

// C++ standard library headers
#include <iostream>

// Vampire headers
#include "atoms.hpp"
#include "cells.hpp"
#include "material.hpp"
#include "micromagnetic.hpp"
#include "vio.hpp"
#include "vmpi.hpp"

// micromagnetic module headers
#include "internal.hpp"

// use short version of namespace
namespace mm = micromagnetic::internal;

namespace micromagnetic{

   //----------------------------------------------------------------------------
   // Function to initialize micromagnetic module
   //----------------------------------------------------------------------------
   void initialize(int num_local_cells,
                   int num_cells,
                   int num_atoms,
                   int num_materials,
                   std::vector<int> cell_array,                     //1D array storing which cell each atom is in
                   std::vector<int> neighbour_list_array,           //1D vector listing the nearest neighbours of each atom
                   std::vector<int> neighbour_list_start_index,     //1D vector storing the start index for each atom in the neighbour_list_array
                   std::vector<int> neighbour_list_end_index,       //1D vector storing the end index for each atom in the neighbour_list_array
                   const std::vector<int> type_array,               //1D array storing which material each cell is
                   std::vector <mp::materials_t> material,          //1D vector of type material_t stiring the material properties
                   std::vector <double> x_coord_array,
                   std::vector <double> y_coord_array,
                   std::vector <double> z_coord_array,
                   std::vector <double> volume_array,               //1D vector storing the volume of each cell
                   double Temperature,
                   double num_atoms_in_unit_cell,
                   double system_dimensions_x,
                   double system_dimensions_y,
                   double system_dimensions_z,
                   std::vector<int> local_cell_array){

      //--------------------------------------------------------------------
      // If micromagnetics not enabled, do nothing
      //--------------------------------------------------------------------
      if(micromagnetic::enabled == false ) return;

      zlog << zTs() << "Initialising micromagnetic calculation..." << std::endl;

      //resizes the vectors used to store the cell parameters
      mm::A.resize(num_cells*num_cells,0.0);
      mm::alpha.resize(num_cells,0.0);
      mm::one_o_chi_perp.resize(num_cells,0.0);
      mm::one_o_chi_para.resize(num_cells,0.0);
      mm::gamma.resize(num_cells,0.0);
      mm::ku.resize(num_cells,0.0);
      mm::ms.resize(num_cells,0.0);
      mm::Tc.resize(num_cells,0.0);
      mm::alpha_para.resize(num_cells,0.0);
      mm::alpha_perp.resize(num_cells,0.0);
      mm::m_e.resize(num_cells,1.0); // initialise as fully magnetizated
      mm::macro_neighbour_list_start_index.resize(num_cells,0.0);
      mm::macro_neighbour_list_end_index.resize(num_cells,0.0);
      micromagnetic::cell_discretisation_micromagnetic.resize(num_cells,true);
      mm::ext_field.resize(3,0.0);
      mm::pinning_field_x.resize(num_cells,0.0);
      mm::pinning_field_y.resize(num_cells,0.0);
      mm::pinning_field_z.resize(num_cells,0.0);
      mm::cell_material_array.resize(num_cells,0.0);

      mm::thermal_field_array_x.resize(num_cells,0.0);
      mm::thermal_field_array_y.resize(num_cells,0.0);
      mm::thermal_field_array_z.resize(num_cells,0.0);

      // These functions vectors with the parameters calcualted from the function
      mm::ms =                   mm::calculate_ms(num_local_cells,num_atoms,num_cells, cell_array, type_array,material,local_cell_array);
      mm::alpha =                mm::calculate_alpha(num_local_cells,num_atoms, num_cells, cell_array, type_array, material,local_cell_array);
      mm::Tc =                   mm::calculate_tc(num_local_cells, local_cell_array,num_atoms, num_cells, cell_array,neighbour_list_array,
                                                  neighbour_list_start_index, neighbour_list_end_index, type_array, material);
      mm::ku =                   mm::calculate_ku(num_atoms, num_cells, cell_array, type_array, material);
      mm::gamma =                mm::calculate_gamma(num_atoms, num_cells, cell_array,type_array,material,num_local_cells,local_cell_array);
      mm::one_o_chi_para =       mm::calculate_chi_para(num_local_cells, local_cell_array,num_cells, Temperature);
      mm::one_o_chi_perp =       mm::calculate_chi_perp(num_local_cells, local_cell_array,num_cells, Temperature);
      mm::A =                    mm::calculate_a(num_atoms, num_cells, num_local_cells,cell_array, neighbour_list_array, neighbour_list_start_index,
                                                 neighbour_list_end_index, type_array,  material, volume_array, x_coord_array,
                                                 y_coord_array, z_coord_array, num_atoms_in_unit_cell, local_cell_array);


// for (int cell = 0; cell < num_cells; cell++)
// std::cout << cells::pos_and_mom_array[4*cell+0] << '\t' << cells::pos_and_mom_array[4*cell+1] << '\t' << cells::pos_and_mom_array[4*cell+2] << '\t' << mm::ms[cell] <<std::endl;

      //----------------------------------------------------------
      // Check for antiferromagnetic cells
      //----------------------------------------------------------
      if (discretisation_type == micromagnetics){
         for (int lc = 0; lc < num_local_cells; lc++){
            int cell = local_cell_array[lc];
            // Force atomistic simulation for anti-ferromagnets
            if (mm::Tc[cell] < 0.0) {
               discretisation_type = multiscale;

            }
         }
      }

      //------------------------------------------------------------------------------------------
      // if multiscale simulation work out which cells/atoms are micromagnetic/atomistic
      //------------------------------------------------------------------------------------------
      if (discretisation_type == multiscale){
         //loops over all atoms and if any atom in the cell is atomistic the whole cell becomes atomistic else the cell is micromagnetic.
         for (int atom =0; atom < num_atoms; atom++){
            int cell = cell_array[atom];
            int mat  = type_array[atom];
            micromagnetic::cell_discretisation_micromagnetic[cell] = mp::material[mat].micromagnetic_enabled;
            //unless the cell contains AFM atoms, then it is always atomistic
            if (mm::Tc[cell] < 0.0) micromagnetic::cell_discretisation_micromagnetic[cell] = 0;

         }

         // loops over all atoms saves each atom at micromagnetic or
         // atomistic depending on whether the cell is microamgnetic or atomistic
         for (int atom =0; atom < num_atoms; atom++){
            int cell = cell_array[atom];
            //id atomistic add to numner of atomisic atoms
            if (micromagnetic::cell_discretisation_micromagnetic[cell] == 0) {
               list_of_atomistic_atoms.push_back(atom);
               number_of_atomistic_atoms++;
            }
            // if micromagnetic add to the number of micromagnetic cells
            else {
               list_of_none_atomistic_atoms.push_back(atom);
               number_of_none_atomistic_atoms++;
            }
         }

         // if simulation is micromagetic all cells are made micromagnetic cells
         for (int lc = 0; lc < num_local_cells; lc++){
            int cell = local_cell_array[lc];
            // Check that cell is micromagnetic, has a reasonable moment and reasonable Curie temperature
            if (micromagnetic::cell_discretisation_micromagnetic[cell] == 1 && mm::ms[cell] > 1e-30 && mm::Tc[cell] > 0.1) {
               list_of_micromagnetic_cells.push_back(cell);
               number_of_micromagnetic_cells++;
            }
            // Otherwise list cells as non-magnetic (empty)
            else{
               list_of_empty_micromagnetic_cells.push_back(cell);
            }
         }

      } // end of multiscale initialisation

      //----------------------------------------------------------------------------------------------
      // if not multiscale then micromagnetic simulation:
      // all cells are micromagnetic and all atoms are micromagnetic
      //----------------------------------------------------------------------------------------------
      else {
         // decide which micromagnetic cells to integrate
         for (int lc = 0; lc < num_local_cells; lc++){

            // get cell ID
            int cell = local_cell_array[lc];

            // Check that cell is micromgantic, has a reasonable moment and reasonable Curie temperature
            if(mm::ms[cell] > 1e-30 && mm::Tc[cell] > 0.1){

               // save cell in list
               list_of_micromagnetic_cells.push_back(cell);

               // update number of cells
               number_of_micromagnetic_cells ++;

            }
            // Otherwise list cells as non-magnetic (empty)
            else{
               list_of_empty_micromagnetic_cells.push_back(cell);
            }

         }

         // remove all atomistic atoms from integration in miromagnetic mode
         for (int atom =0; atom < num_atoms; atom++){
            list_of_none_atomistic_atoms.push_back(atom);
            number_of_none_atomistic_atoms++;
         }

      }

      //--------------------------------------------------------------------------------------------
      // for field calcualtions you need to access the atoms in numerically consecutive lists.
      // therefore you need to create lists of consecutive lists
      // loops over all atoms if atom is not one minus the previous atom then create a new list.
      //--------------------------------------------------------------------------------------------
      if (number_of_atomistic_atoms > 0){
         int end = list_of_atomistic_atoms[0];
         int begin = list_of_atomistic_atoms[0];
         for(int atom_list=1;atom_list<number_of_atomistic_atoms;atom_list++){
            int atom = list_of_atomistic_atoms[atom_list];
            int last_atom = list_of_atomistic_atoms[atom_list - 1];
            if ((atom != last_atom +1) || (atom_list == number_of_atomistic_atoms -1)){
               end = atom +1;
               mm::fields_neighbouring_atoms_begin.push_back(begin);
               mm::fields_neighbouring_atoms_end.push_back(end);

               begin = atom + 1;
            }
         }
      }

      // Print informative message to screen if there are no cells integrated and no atoms integrated
      if(number_of_micromagnetic_cells == 0 && number_of_atomistic_atoms == 0){
         terminaltextcolor(YELLOW);
         std::cout << "Warning! No micromagnetic cells or atoms will be integrated (Tc or Ms is too low in all cells)" << std::endl;
         terminaltextcolor(WHITE);
         zlog << zTs() << "Warning! No micromagnetic cells or atoms will be integrated (Tc or Ms is too low in all cells)" << std::endl;
      }

     /*for (int cell = 0; cell < num_cells; cell++){
       std::cerr << '\t' << cells::cell_coords_array_z[cell] << '\t' <<  mm::ms[cell] << '\t' << mm::ku[cell] << '\t' << mm::A[cell] << "\t" << mm::Tc[cell] << "\t" <<micromagnetic::cell_discretisation_micromagnetic[cell] << "\t" << mm::alpha[cell] << '\t' << mm::gamma[cell] << std::endl;
    }*/

    //--------------------------------------------------------------------------------------------
    // What does this do?
    //--------------------------------------------------------------------------------------------

     std::vector < double > temp(num_cells,0);

     //int num_calculations = mm::fields_neighbouring_atoms_begin.size();

     for (int atom = 0; atom < num_atoms; atom ++){
        int mat = type_array[atom];
        int cell = cell_array[atom];
        mm::cell_material_array[cell] = mat;
     }

     for (int cell = 0; cell < num_cells; cell++ ){

        double zi = cells::pos_and_mom_array[4*cell+2];
         int mat = mm::cell_material_array[cell];
        // if ( mp::material[mat].pinning_field_unit_vector[1] != 0) std::cout << mat << "\t" << mp::material[mat].pinning_field_unit_vector[1] << '\t' << mm::pinning_field_height << "\t" << zi << std::endl;
         if (zi < mm::pinning_field_height && mp::material[mat].pinning_field_unit_vector[0]+ mp::material[mat].pinning_field_unit_vector[1] + mp::material[mat].pinning_field_unit_vector[2]!= 0.0){
           double Area = cells::macro_cell_size[0]*cells::macro_cell_size[1];
         //  std::cout << cell << '\t' << mat << '\t' <<  mp::material[mat].pinning_field_unit_vector[0] << '\t' <<  mp::material[mat].pinning_field_unit_vector[1] <<std::endl;
            mm::pinning_field_x[cell] = Area*mp::material[mat].pinning_field_unit_vector[0]/mm::ms[cell];
            mm::pinning_field_y[cell] = Area*mp::material[mat].pinning_field_unit_vector[1]/mm::ms[cell];
            mm::pinning_field_z[cell] = Area*mp::material[mat].pinning_field_unit_vector[2]/mm::ms[cell];
          //  std::cout << mp::material[mat].pinning_field_unit_vector[1] << "\t" << Area << '\t'  << mm::ms[cell] <<  '\t' <<mm::pinning_field_x[cell] << '\t' << mm::pinning_field_y[cell] << '\t' << mm::pinning_field_z[cell] <<std::endl;
        }
     }


   //  std::cout << mm::mm_correction <<std::endl;
  //   if (mm::mm_correction == true){
  //      for (int cell = 0; cell < num_cells; cell++ ){
  //         mm::pinning_field_x[cell] = 2*mm::pinning_field_x[cell]/cells::macro_cell_size[2];
  //         mm::pinning_field_y[cell] = 2*mm::pinning_field_y[cell]/cells::macro_cell_size[2];
  //         mm::pinning_field_z[cell] = 2*mm::pinning_field_z[cell]/cells::macro_cell_size[2];
//           std::cout << cell << '\t' << mm::pinning_field_x[cell] <<'\t' << mm::pinning_field_y[cell] <<'\t' << mm::pinning_field_z[cell] <<std::endl;
    //    }
  //}

      //boltzman stuff
      P.resize(101);
      for (int i = 0; i < 101; i++) P[i].resize(101,0.0);

      //--------------------------------------------
      // print informative message to users
      //--------------------------------------------
      int total_micromagnetic_cells = number_of_micromagnetic_cells;
      #ifdef MPICF
         MPI_Allreduce(&number_of_micromagnetic_cells, &total_micromagnetic_cells, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
      #endif

      std::cout << "Micromagnetic initialisation complete" << std::endl;
      zlog << zTs() << "Micromagnetic initialisation complete" << std::endl;
      std::cout << "Number of micromagnetic cells (all CPUs): " << total_micromagnetic_cells << std::endl;
      zlog << zTs() << "Number of micromagnetic cells: " << total_micromagnetic_cells << std::endl;

      return;

   }

} // end of micromagnetic namespace
