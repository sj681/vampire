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


namespace dipole{

   namespace internal{

      void hierarchical_tensor(){

         // int level = 0;
         //
         // for (int cell = 0; cell < num_cells_level[level]; cell++)
         // std::cout << dp::num_atoms_in_cell_level[level][cell] << '\t' << std::endl;

      }
   }
}
