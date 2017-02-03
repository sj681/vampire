
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


   void mm::stepLLG( int num_cells,
                  double temperature,
                  std::vector<double> x_array,
                  std::vector<double> y_array,
                  std::vector<double> z_array,
                  std::vector<double> ext_field,
                  double dt,
                  std::vector<double>& new_x_array,
                  std::vector<double>& new_y_array,
                  std::vector<double>& new_z_array){


      const double kB = 1.3806503e-23;

      std::vector<double> m(3,0.0);
      std::vector<double> spin_field(3,0.0);


      //calculte chi(T).
      mm::chi_para =  mm::calculate_chi_para(num_cells, temperature);
      mm::chi_perp =  mm::calculate_chi_perp(num_cells, temperature);

      //loops over all the cells to calculate the spin terms per cell - only filled cells where MS>0
      for (int cell =0; cell <num_cells; cell++){
         if (mm::ms[cell] > 1e-100){
            m[0] = x_array[cell];
            m[1] = y_array[cell];
            m[2] = z_array[cell];

            spin_field = mm::calculate_fields(m, true, temperature, num_cells, cell, x_array,y_array,z_array);

            const double one_oneplusalpha_sq = 1/(1+mm::alpha[cell]*mm::alpha[cell]); // material specific alpha and gamma
            const double alpha_oneplusalpha_sq = mm::alpha[cell]/(1+mm::alpha[cell]*mm::alpha[cell]);

            const double S[3] = {m[0],m[1],m[2]};
            const double H[3] = {-spin_field[0], -spin_field[1], -spin_field[2]};

      double xyz[3] = {0.0,0.0,0.0};
            xyz[0]=(one_oneplusalpha_sq)*(S[1]*H[2]-S[2]*H[1]) + (alpha_oneplusalpha_sq)*(S[1]*(S[0]*H[1]-S[1]*H[0])-S[2]*(S[2]*H[0]-S[0]*H[2]));
            xyz[1]=(one_oneplusalpha_sq)*(S[2]*H[0]-S[0]*H[2]) + (alpha_oneplusalpha_sq)*(S[2]*(S[1]*H[2]-S[2]*H[1])-S[0]*(S[0]*H[1]-S[1]*H[0]));
            xyz[2]=(one_oneplusalpha_sq)*(S[0]*H[1]-S[1]*H[0]) + (alpha_oneplusalpha_sq)*(S[0]*(S[2]*H[0]-S[0]*H[2])-S[1]*(S[1]*H[2]-S[2]*H[1]));

            // Store dS in euler array
            new_x_array[cell] = xyz[0];
            new_y_array[cell] = xyz[1];
            new_z_array[cell] = xyz[2];

         }
      }
   }

   void mm::step( int num_cells,
                  double temperature,
                  std::vector<double> x_array,
                  std::vector<double> y_array,
                  std::vector<double> z_array,
                  std::vector<double> ext_field,
                  double dt,
                  std::vector<double>& new_x_array,
                  std::vector<double>& new_y_array,
                  std::vector<double>& new_z_array){

      const double kB = 1.3806503e-23;


      std::vector<double> m(3,0.0);
      std::vector<double> spin_field(3,0.0);

      //6 arrays of gaussian random numbers to store the stochastic noise terms for x,y,z parallel and perperdicular
      std::vector <double> GW1x(num_cells);
      std::vector <double> GW1y(num_cells);
      std::vector <double> GW1z(num_cells);
      std::vector <double> GW2x(num_cells);
      std::vector <double> GW2y(num_cells);
      std::vector <double> GW2z(num_cells);

      //fill the noise terms
      generate (GW1x.begin(),GW1x.end(), mtrandom::gaussian);
      generate (GW1y.begin(),GW1y.end(), mtrandom::gaussian);
      generate (GW1z.begin(),GW1z.end(), mtrandom::gaussian);
      generate (GW2x.begin(),GW2x.end(), mtrandom::gaussian);
      generate (GW2y.begin(),GW2y.end(), mtrandom::gaussian);
      generate (GW2z.begin(),GW2z.end(), mtrandom::gaussian);


   //loops over all the cells to calculate the spin terms per cell - only filled cells where MS>0
   for (int cell =0; cell <num_cells; cell++){
      if (mm::ms[cell] > 1e-100){
         m[0] = x_array[cell];
         m[1] = y_array[cell];
         m[2] = z_array[cell];

         spin_field = mm::calculate_fields(m,false, temperature, num_cells, cell, x_array,y_array,z_array);

         //calculates the stochatic parallel and perpendicular terms
         double a;
         if (micromagnetic::stochastic == true) a = 1.0;
         else if (micromagnetic::stochastic == false) a = 0.0;
         double sigma_para = a*sqrt(2*kB*temperature*mm::alpha_para[cell]/(mm::ms[cell]*dt)); //why 1e-27
         double sigma_perp = a*sqrt(2*kB*temperature*(mm::alpha_perp[cell]-mm::alpha_para[cell])/(dt*mm::ms[cell]*mm::alpha_perp[cell]*mm::alpha_perp[cell]));

         const double H[3] = {spin_field[0], spin_field[1], spin_field[2]};

         //saves the noise terms to an array
         const double GW2t[3] = {GW2x[cell],GW2y[cell],GW2z[cell]};
         const double one_o_m_squared = 1.0/(m[0]*m[0]+m[1]*m[1]+m[2]*m[2]);
         const double SdotH = m[0]*H[0] + m[1]*H[1] + m[2]*H[2];


      double xyz[3] = {0.0,0.0,0.0};
         //calculates the LLB equation
         xyz[0]= 	- (m[1]*H[2]-m[2]*H[1])
                  + mm::alpha_para[cell]*m[0]*SdotH*one_o_m_squared
                  - mm::alpha_perp[cell]*(m[1]*(m[0]*H[1]-m[1]*H[0])-m[2]*(m[2]*H[0]-m[0]*H[2]))*one_o_m_squared
                  + GW1x[cell]*sigma_para
                  - mm::alpha_perp[cell]*(m[1]*(m[0]*GW2t[1]-m[1]*GW2t[0])-m[2]*(m[2]*GW2t[0]-m[0]*GW2t[2]))*one_o_m_squared*sigma_perp;

        xyz[1]= 	- (m[2]*H[0]-m[0]*H[2])
                  + mm::alpha_para[cell]*m[1]*SdotH*one_o_m_squared
                  - mm::alpha_perp[cell]*(m[2]*(m[1]*H[2]-m[2]*H[1])-m[0]*(m[0]*H[1]-m[1]*H[0]))*one_o_m_squared
                  + GW1y[cell]*sigma_para
                  - mm::alpha_perp[cell]*(m[2]*(m[1]*GW2t[2]-m[2]*GW2t[1])-m[0]*(m[0]*GW2t[1]-m[1]*GW2t[0]))*one_o_m_squared*sigma_perp;

        xyz[2]=	- (m[0]*H[1]-m[1]*H[0])
                  + mm::alpha_para[cell]*m[2]*SdotH*one_o_m_squared
                  - mm::alpha_perp[cell]*(m[0]*(m[2]*H[0]-m[0]*H[2])-m[1]*(m[1]*H[2]-m[2]*H[1]))*one_o_m_squared
                  + GW1z[cell]*sigma_para
                  - mm::alpha_perp[cell]*(m[0]*(m[2]*GW2t[0]-m[0]*GW2t[2])-m[1]*(m[1]*GW2t[2]-m[2]*GW2t[1]))*one_o_m_squared*sigma_perp;

                  //returns the values of the LLB as the new
                  new_x_array[cell] = xyz[0];
                  new_y_array[cell] = xyz[1];
                  new_z_array[cell] = xyz[2];

         }

      }
   }
}
