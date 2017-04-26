// Vampire headers
#include "micromagnetic.hpp"

// micromagnetic module headers
#include "internal.hpp"
#include "material.hpp"

namespace micromagnetic{
   namespace internal{

      //calculates the magnetocrystaline anisotropy for each cell from the individual atomic anisotropies
      //ku = sum (ku) for each cell
      std::vector<double> calculate_ku(const int num_atoms,
                                       const int num_cells,
                                       std::vector<int> cell_array,                  //1D array storing which cell each atom is in
                                       const std::vector<int> type_array,            //1D array storing which material each atom is
                                       std::vector <mp::materials_t> material){      //class of material parameters for the atoms

         //ku = sum ku
         std::vector<double> ku(num_cells,0.0);
         //sums over all atoms to add the ku per cell
         for (int atom = 0; atom < num_atoms; atom++){

            ku[cell_array[atom]] = ku[cell_array[atom]] +  mp::material[type_array[atom]].Ku1_SI;
         }
         return ku;               //returns a 
      }
   } //closes the internal namspace
}  //closes the micromagnetic namespace
