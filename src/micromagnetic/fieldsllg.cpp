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

// Vampire headers
#include "cells.hpp"
#include "dipole.hpp"
#include "errors.hpp"
#include "environment.hpp"
#include "random.hpp"
#include "sim.hpp"
#include "vio.hpp"

// micromagnetic module headers
#include "micromagnetic.hpp"
#include "internal.hpp"

namespace micromagnetic{

   namespace internal{

      std::vector<double> calculate_llg_fields(std::vector <double>& m,
                                               double temperature,
                                               int num_cells,
                                               int cell,
                                               std::vector<double>& x_array,
                                               std::vector<double>& y_array,
                                               std::vector<double>& z_array){


      std::vector<double> spin_field(3,0.0);

      const double kB = 1.3806503e-23;

      //chi is usually used as 2/chi
      const double one_o_2_chi_para = (one_o_chi_para[cell]/2.0);

      //the temperature is usually used as a reduced temperature.
      const double reduced_temperature = temperature/Tc[cell];
      const double Tc_o_Tc_m_T = Tc[cell]/(temperature - Tc[cell]);

      //calcualted m_e and alpha temperature dependant
      if(temperature<=Tc[cell]){
         m_e[cell] = pow((Tc[cell]-temperature)/(Tc[cell]),0.365);
         alpha_para[cell] = (2.0/3.0)*alpha[cell]*reduced_temperature;
         alpha_perp[cell] = alpha[cell]*(1.0-temperature/(3.0*Tc[cell]));
      }
      else{
         m_e[cell] = 0.01;
         alpha_para[cell] = alpha[cell]*(2.0/3.0)*reduced_temperature;
         alpha_perp[cell] = alpha_para[cell];
      }

      //saved me_2
      const double m_e_squared = m_e[cell]*m_e[cell];

      //calculates the intercell exchange field (pf) - this is dependent on temperature.
      double pf = 0;
      if(temperature<=Tc[cell]) pf = 0;
      else pf = -2.0*one_o_2_chi_para*(1.0 + Tc_o_Tc_m_T*3.0*m_e_squared/5.0);

      //calculates the exchage fields as me^1.66 *A*(xi-xj)/m_e^2
      double exchange_field[3]={0.0,0.0,0.0};
      //is T < TC the exchange field = 0
      if (num_cells > 1){


         const int start = macro_neighbour_list_start_index[cell];
         const int end = macro_neighbour_list_end_index[cell];
         for(int j = start;j<end;j++){
            // calculate reduced exchange constant factor
            const int cellj = macro_neighbour_list_array[j];
            const double mj = m_e[cellj];
            double Ac = A[j]*pow(mj,1.66);

            int mat  = cell_material_array[cell];
            int matj =cell_material_array[cellj];
           if (mp::material[mat].enable_SAF == true && mp::material[matj].enable_SAF == true){
              if (mat != matj){

                 double Area = cells::macro_cell_size[0]*cells::macro_cell_size[1];
                Ac = -pow(mj,1.66)*Area*mp::material[mat].SAF[matj]/ms[cell];
                //if (mm_correction == true) Ac = 2*Ac/cells::macro_cell_size[2];

              }
           }
      //      std::cout << mat << '\t' << matj << "\t" << Ac <<std::endl;

            exchange_field[0] -= Ac*(x_array[cellj]*m_e[cellj] - x_array[cell]*m_e[cell]);
            exchange_field[1] -= Ac*(y_array[cellj]*m_e[cellj] - y_array[cell]*m_e[cell]);
            exchange_field[2] -= Ac*(z_array[cellj]*m_e[cellj] - z_array[cell]*m_e[cell]);



         }

      }

      //calcualtes thesigma values
      double sigma_para = sqrt(2*kB*temperature*alpha_para[cell]/(ms[cell]*mp::dt));
      double sigma_perp = sqrt(2*kB*temperature*(alpha_perp[cell]-alpha_para[cell])/(mp::dt*ms[cell]*alpha_perp[cell]*alpha_perp[cell]));

      //Sum H = H_exch + H_A +H_exch_grains +H_App + H+dip
      spin_field[0] = ext_field[0] + exchange_field[0] + one_o_chi_perp[cell]*m[0]*m_e[cell] + sigma_para*mtrandom::gaussian() + pinning_field_x[cell];// + sim::track_field_x[cell];
      spin_field[1] =ext_field[1] + exchange_field[1] + one_o_chi_perp[cell]*m[1]*m_e[cell]  + sigma_para*mtrandom::gaussian() + pinning_field_y[cell];// + sim::track_field_y[cell];
      spin_field[2] =ext_field[2] + exchange_field[2]                                        + sigma_para*mtrandom::gaussian() + pinning_field_z[cell];// + sim::track_field_z[cell];

      // Add dipole field if enabled
      if (dipole::activated){
         spin_field[0] += dipole::cells_field_array_x[cell];
         spin_field[1] += dipole::cells_field_array_x[cell];
         spin_field[2] += dipole::cells_field_array_x[cell];
      }

      if (sim::enable_fmr){
         spin_field[0] = spin_field[0] + fmr_H[0];
         spin_field[1] = spin_field[1] + fmr_H[1];
         spin_field[2] = spin_field[2] + fmr_H[2];
      }

      if (environment::enabled){
         spin_field[0] = spin_field[0] + environment::environment_field_x[cell];
         spin_field[1] = spin_field[1] + environment::environment_field_y[cell];
         spin_field[2] = spin_field[2] + environment::environment_field_z[cell];
      }

      if (sim::track_field_x.size() != 0 ){
        spin_field[0] = spin_field[0] + sim::track_field_x[cell];
        spin_field[1] = spin_field[1] + sim::track_field_y[cell];
        spin_field[2] = spin_field[2] + sim::track_field_z[cell];
      }



      return spin_field;

   }

} // end of internal namespace

} // end of micromagnetic namespace
