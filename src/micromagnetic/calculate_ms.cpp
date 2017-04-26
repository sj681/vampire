
// Vampire headers
#include "micromagnetic.hpp"

// micromagnetic module headers
#include "internal.hpp"

namespace micromagnetic {
   namespace internal {


      //calcualtes the saturation magnetisation of each cell as a sum of the magnetic moment of each atom in the cell.
      //ms = sum muS for each cell
      std::vector<double> calculate_ms(const int num_atoms,
                                       const int num_cells,
                                       std::vector<int> cell_array,                  //1D array storing which cell each atom is in
                                       const std::vector<int> type_array,            //1D array storing which material each atom is
                                       std::vector <mp::materials_t> material){      //class of material parameters for the atoms
         //stores ms for each cell
         std::vector<double> ms(num_cells,0.0);

         //sums over all atoms to sum the muS per cell
         for (int atom = 0; atom < num_atoms; atom++) ms[cell_array[atom]] = ms[cell_array[atom]] + material[type_array[atom]].mu_s_SI;

         return ms;           //returns a 1D array containg the saturation magnetisation of every cell 
      }
   } //closes the internal namspace
}  //closes the micromagnetic namespace
