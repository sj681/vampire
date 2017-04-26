
// Vampire headers
#include "micromagnetic.hpp"

// micromagnetic module headers
#include "internal.hpp"
#include "material.hpp"

namespace micromagnetic {

   namespace internal {


      //calculates the curie temperature of each cell
      std::vector<double> calculate_tc(int num_atoms,
                                       int num_cells,
                                       std::vector<int> cell_array,                   //1D array storing which cell each atom is in
                                       std::vector<int> neighbour_list_array,         //1D list of the interactions between atoms,
                                       std::vector<int> neighbour_list_start_index,   //the list for each atom in neighbour list array starts at the start index
                                       std::vector<int> neighbour_list_end_index,     //and ends at the end index - stored in these vectors
                                       const std::vector<int> type_array,            //1D array storing which material each atom is
                                       std::vector <mp::materials_t> material){      //class of material parameters for the atoms

         std::vector<double>  J(num_cells,0.0);         //stores the exchange constant for each cell
         std::vector<double>  N(num_cells,0.0);         //1D vector containg the number of atoms in each cell.
         std::vector<double>  Tc(num_cells,0.0);        //1D vector storing the curie temperature of each cell
         const double e = 1.3;                          //mea field constant
         const double kB = 1.38064852e-23;              //boltzman constant

         //-----------------------------------------------------------------------

         //             TC = sum(Jij)* e/(3kB N)

         //------------------------------------------------------------------------

         //Sums Jij for each cell and calculates the number of atoms in each cell
         for (int atom = 0; atom <num_atoms; atom++){
            N[cell_array[atom]]++;
            for (int neighbour = neighbour_list_start_index[atom]; neighbour < neighbour_list_end_index[atom] +1; neighbour ++){
               J[cell_array[atom]] = J[cell_array[atom]] +  mp::material[type_array[atom]].Jij_matrix_SI[type_array[neighbour_list_array[neighbour]]];
            }
         }


         //calculates Tc for each cell
         for (int cell = 0; cell < num_cells; cell++){
            Tc[cell] = -J[cell]*e/(3*kB*N[cell]);
         }
         return Tc;             //returns a 1D array containing the curiue temepratures
      }
   } //closes the internal namspace
}  //closes the micromagnetic namespace
