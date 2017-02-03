
// Vampire headers
#include "micromagnetic.hpp"

// micromagnetic module headers
#include "internal.hpp"

#include "random.hpp"
#include "sim.hpp"
#include "vmpi.hpp"
#include "cells.hpp"
#include "atoms.hpp"
#include "vio.hpp"
#include <cmath>
#include <iostream>
#include <algorithm>
#include <fstream>


namespace mm = micromagnetic::internal;


int LLG_serial_heun(int num_steps,int num_cells,double temperature,std::vector<double>& x_mag_array,std::vector<double>& y_mag_array,std::vector<double>& z_mag_array,double Hx,double Hy,double Hz, double H,double dt,std::vector <double> volume_array, int N);

namespace micromagnetic{
/// Master LLB Function - dispatches code path to desired LLB routine
/// \f$ \frac{\partial S}{\partial t} \f$
int LLG(int num_steps,int num_cells,double temperature,std::vector<double>& x_mag_array,std::vector<double>& y_mag_array,std::vector<double>& z_mag_array,double Hx,double Hy,double Hz, double H,double dt, std::vector <double> volume_array, int N){

   //----------------------------------------------------------
   // check calling of routine if error checking is activated
   //----------------------------------------------------------


   //if(err::check==true){std::cout << "LLB has been called" << std::endl;}

   #ifdef MPICF
      //LLB_mpi(num_steps);
   #else
   LLG_serial_heun(num_steps,num_cells,temperature,x_mag_array,y_mag_array,z_mag_array,Hx,Hy,Hz,H,dt, volume_array, N);
   #endif

   return 0;
}
}


int LLG_serial_heun( int num_steps,
                     int num_cells,
                     double temperature,
                     std::vector<double>& x_mag_array,
                     std::vector<double>& y_mag_array,
                     std::vector<double>& z_mag_array,
                     double Hx,
                     double Hy,
                     double Hz,
                     double H,
                     double dt,
                     std::vector <double> volume_array,
                     int N
                  ){

   //The external fields equal the length of the field times the applied field vector.
   //This is saved to an array.
   mm::ext_field[0] = H*Hx;
   mm::ext_field[1] = H*Hy;
   mm::ext_field[2] = H*Hz;

   //loop over all cells and calculate m = M/Ms and save it to x,y,z_array
   //there are lots of blank cells so only the cells where ms!=0 are calcualted.
   //save this new m as the initial value, so it can be saved and used in the final equation.
   for (int cell = 0; cell < num_cells; cell++){
      mm::x_array[cell] = x_mag_array[cell];
      mm::y_array[cell] = y_mag_array[cell];
      mm::z_array[cell] = z_mag_array[cell];
      double m_squared = sqrt(mm::x_array[cell]*mm::x_array[cell] + mm::y_array[cell]*mm::y_array[cell] + mm::z_array[cell]*mm::z_array[cell]);
      mm::x_array[cell] = mm::x_array[cell]/m_squared;
      mm::y_array[cell] = mm::y_array[cell]/m_squared;
      mm::z_array[cell] = mm::z_array[cell]/m_squared;

   }


   mm::stepLLG(num_cells, temperature, mm::x_array,mm::y_array,mm::z_array, mm::ext_field, dt, mm::x_euler_array, mm::y_euler_array, mm::z_euler_array);


   double S_new[3] = {0.0,0.0,0.0};
   double mod_S;
	// Calculate Euler Step
   for (int cell = 0; cell < num_cells; cell++){
      if (mm::ms[cell] > 0){

         S_new[0] = mm::x_array[cell] + mm::x_euler_array[cell]*dt;
         S_new[1] = mm::y_array[cell] + mm::y_euler_array[cell]*dt;
         S_new[2] = mm::z_array[cell] + mm::z_euler_array[cell]*dt;

	      // Normalise Spin Length
	      mod_S = 1.0/sqrt(S_new[0]*S_new[0] + S_new[1]*S_new[1] + S_new[2]*S_new[2]);

	      S_new[0]=S_new[0]*mod_S;
	      S_new[1]=S_new[1]*mod_S;
	      S_new[2]=S_new[2]*mod_S;

	      //Writing of Spin Values to Storage Array
	      mm::mx_store[cell]=S_new[0];
	      mm::my_store[cell]=S_new[1];
	      mm::mz_store[cell]=S_new[2];

      }
   }


mm::stepLLG(num_cells, temperature, mm::mx_store, mm::my_store, mm::mz_store, mm::ext_field, dt, mm::x_heun_array, mm::y_heun_array, mm::z_heun_array);

for (int cell = 0; cell < num_cells; cell++){
   if (mm::ms[cell] > 0){

      S_new[0]=mm::x_array[cell]+(dt/2)*(mm::x_euler_array[cell]+mm::x_heun_array[cell]);
      S_new[1]=mm::y_array[cell]+(dt/2)*(mm::y_euler_array[cell]+mm::y_heun_array[cell]);
      S_new[2]=mm::z_array[cell]+(dt/2)*(mm::z_euler_array[cell]+mm::z_heun_array[cell]);

      // Normalise Spin Length
      mod_S = 1.0/sqrt(S_new[0]*S_new[0] + S_new[1]*S_new[1] + S_new[2]*S_new[2]);

      S_new[0]=S_new[0]*mod_S;
      S_new[1]=S_new[1]*mod_S;
      S_new[2]=S_new[2]*mod_S;

      // Copy new spins to spin array
      mm::x_array[cell]=S_new[0];
      mm::y_array[cell]=S_new[1];
      mm::z_array[cell]=S_new[2];

      cells::x_mag_array[cell] = mm::x_array[cell]*mm::ms[cell];
      cells::y_mag_array[cell] = mm::y_array[cell]*mm::ms[cell];
      cells::z_mag_array[cell] = mm::z_array[cell]*mm::ms[cell];
   }
}



   double x = 0,y = 0,z=0, l = 0, ml= 0;

   for (int cell = 0; cell < num_cells; cell ++)
   {
      x = x + mm::x_array[cell]*mm::m_e[cell];
      y = y + mm::y_array[cell]*mm::m_e[cell];
      z = z + mm::z_array[cell]*mm::m_e[cell];
      l = l + mm::m_e[cell];
      ml= ml + 1;

   }


   if((sim::time +1)%vout::output_rate==0){

      for (int atom = 0; atom < atoms::num_atoms; atom ++)
      {
         int cell = atoms::cell_array[atom];
         atoms::x_spin_array[atom] = mm::x_array[cell]*mm::m_e[cell];
         atoms::y_spin_array[atom] = mm::y_array[cell]*mm::m_e[cell];
         atoms::z_spin_array[atom] = mm::z_array[cell]*mm::m_e[cell];
         atoms::m_spin_array[atom] = mm::m_e[cell];
      }

   }

}
