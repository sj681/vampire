

// Vampire headers
#include "micromagnetic.hpp"

// micromagnetic module headers
#include "internal.hpp"

#include "random.hpp"
#include "sim.hpp"
#include "cells.hpp"
#include <cmath>
#include <iostream>
#include <algorithm>
#include <fstream>


namespace mm = micromagnetic::internal;

namespace micromagnetic{

std::vector<double> mm::calculate_fields(std::vector <double > m, bool LLG, double temperature, int num_cells, int cell, std::vector<double> x_array, std::vector<double> y_array, std::vector<double> z_array){

      std::vector<double> spin_field(3,0.0);

      double m_squared;
      //chi is usually sued as 2/chi
      const double one_o_chi_perp = 1.0/mm::chi_perp[cell];
      const double one_o_2_chi_para = 1.0/(2.0*mm::chi_para[cell]);

      //the temperature is usually used as a reduced temperature.
      const double reduced_temperature = temperature/mm::Tc[cell];
      const double Tc_o_Tc_m_T = mm::Tc[cell]/(temperature - mm::Tc[cell]);

      //calculates me and alpha depending on the temperature.
      if (temperature<=mm::Tc[cell]){
         mm::m_e[cell] = pow((mm::Tc[cell]-temperature)/(mm::Tc[cell]),0.365);
         mm::alpha_para[cell] = (2.0/3.0)*mm::alpha[cell]*reduced_temperature;
         mm::alpha_perp[cell] = mm::alpha[cell]*(1.0-temperature/(3.0*mm::Tc[cell]));
      }
      else{
         if (LLG) mm::m_e[cell] = 0.0000001;
         else mm::m_e[cell] = 0.0;
         mm::alpha_para[cell] = mm::alpha[cell]*(2.0/3.0)*reduced_temperature;
         mm::alpha_perp[cell] = mm::alpha_para[cell];
      }


      const double m_e_squared = mm::m_e[cell]*mm::m_e[cell];
      if (LLG)    m_squared = m_e_squared;
      else        m_squared = m[1]*m[1]+m[2]*m[2]+m[0]*m[0];

      //calculates the final term of the field pf - this is depnedent on temperature.
      double pf = 0;
      if(temperature<=mm::Tc[cell]) {
         if (LLG) pf = 0.0;
         else pf = one_o_2_chi_para*(1.0 - m_squared/m_e_squared);
      }
      else pf = -2.0*one_o_2_chi_para*(1.0 + Tc_o_Tc_m_T*3.0*m_squared/5.0);
   //   std::cout << m_e <<std::endl;

      //calculates the exchage fields as me^1.66 *A*(xi-xj)/m_e^2
      double exchange_field[3]={0.0,0.0,0.0};
      //is T < TC the exchange field = 0
      int cellj;
      double mj;

      if (num_cells > 1){
         int j2 = cell*num_cells;
         //loops over all other cells to sum the interaction
         for(int j = mm::macro_neighbour_list_start_index[cell];j<mm::macro_neighbour_list_end_index[cell] +1;j++){
            // calculate reduced exchange constant factor

            cellj = mm::macro_neighbour_list_array[j];
            if (LLG) mj = mm::m_e[cell];
            else mj = sqrt(x_array[cellj]*x_array[cellj] + y_array[cellj]*y_array[cellj] + z_array[cellj]*z_array[cellj]);
            const double A = mm::A[cellj]*pow(mj,1.66);
            exchange_field[0] -= A*(x_array[cellj] - x_array[cell]);
            exchange_field[1] -= A*(y_array[cellj] - y_array[cell]);
            exchange_field[2] -= A*(z_array[cellj] - z_array[cell]);
               }
      }
      //Sum H = H_exch + H_A +H_exch_grains +H_App + H+dip
      spin_field[0] = pf*m[0] - one_o_chi_perp*m[0] + mm::ext_field[0] + cells::x_field_array[cell]*cells::num_atoms_in_cell[cell] + exchange_field[0] ;
      spin_field[1] = pf*m[1] - one_o_chi_perp*m[1] + mm::ext_field[1] + cells::y_field_array[cell]*cells::num_atoms_in_cell[cell] + exchange_field[1] ;
      spin_field[2] = pf*m[2]                       + mm::ext_field[2] + cells::z_field_array[cell]*cells::num_atoms_in_cell[cell] + exchange_field[2] ;

      return spin_field;
}
}
