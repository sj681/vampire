
// Vampire headers
#include "micromagnetic.hpp"

// micromagnetic module headers
#include "internal.hpp"
#include "material.hpp"

namespace micromagnetic{

   namespace internal{

      // Calculates the intrinsic damping constant of each cell from the atomsitic parameters in each cell.

      //alpha = sum(alpha)/N for each cell

      std::vector<double> calculate_alpha(int num_atoms,
                                          int num_cells,
                                          std::vector<int> cell_array,                 //1D array storing which cell each atom is in
                                          const std::vector<int> type_array,           //1D array storing which material each atom is
                                          std::vector <mp::materials_t> material){     //class of material parameters for the atoms


         std::vector<double>  alpha(num_cells,0.0);                                    //vectors to store alpha
         std::vector<double>  N(num_cells,0.0);                                        //vector stores the number of atoms per micromagnetic cell

         //sums over all atoms to calulcate the sum of alpha and number of atoms for each cell
         for (int atom = 0; atom <num_atoms; atom++){

            alpha[cell_array[atom]] = alpha[cell_array[atom]] + mp::material[type_array[atom]].alpha;
            N[cell_array[atom]]++;
         }

         //calculates the average alpha per cell
         for (int cell = 0; cell < num_cells; cell++){
            alpha[cell] = alpha[cell]/N[cell];
         }
         return alpha;          //return an array of damping constants for each cell
      }
   } //closes the internal namspace
}  //closes the micromagnetic namespace
