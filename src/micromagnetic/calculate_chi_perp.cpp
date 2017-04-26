// Vampire headers
#include "micromagnetic.hpp"

// micromagnetic module headers
#include "internal.hpp"
#include <math.h>
#include "material.hpp"

namespace micromagnetic{

   namespace internal{
     //calcuaktes the temperature dependant perpendicular component of the susceptability from the analytical equation

      std::vector<double> calculate_chi_perp(int num_cells, double T){

         // -----------------------------------------------------------------------------------

         //               chi_perp = a0 + sum_i ai . (Tc/1.068Tc)^i        if T < TC
         //               chi_perp = b0 +1/4pi . Tc/(T-TC)                 if T > TC

         // -----------------------------------------------------------------------------------

         std::vector<double>  chi(num_cells,0.0);     //1D vector the store the perpendicular susceptability for each cell

         //Fit Parameters
       	double a0 = 0.00211549427182711;
       	double a1 = 0.110224660661792;
       	double a2 = -0.855153260915204;
       	double a3 = 3.42088365387997;
       	double a4 = -7.84821585896818;
       	double a5 = 10.3247035469514;
       	double a6 = -6.85608273303224;
       	double a7 = 0.797453198330591;
       	double a8 = 1.53787854178089;
       	double a9 = -0.627128148404525;
       	double PI= 3.14159;

         //normalisation when T=0
         double chi_0 = (a0+ a1+ a2+ a3 + a4 + a5 + a6 + a7+ a8+ a9);
         //conversion to SI
         chi_0 = chi_0*9.54393845712027;

         for (int cell = 0; cell < num_cells; cell ++){
            //calculates the perpendicular susceptability from the analytical expression
            if(T<(1.068*Tc[cell])) chi[cell] = a0+ a1*((1.068*Tc[cell]-T)/(1.068*Tc[cell]))+ a2*pow((((1.068*Tc[cell])-T)/(1.068*Tc[cell])),2.)+ a3*pow((((1.068*Tc[cell])-T)/(1.068*Tc[cell])),3.)+ a4*pow((((1.068*Tc[cell])-T)/(1.068*Tc[cell])),4.) + a5*pow((((1.068*Tc[cell])-T)/(1.068*Tc[cell])),5.) + a6*pow((((1.068*Tc[cell])-T)/(1.068*Tc[cell])),6.) + a7*pow((((1.068*Tc[cell])-T)/(1.068*Tc[cell])),7.)+ a8*pow((((1.068*Tc[cell])-T)/(1.068*Tc[cell])),8.) + a9*pow((((1.068*Tc[cell])-T)/(1.068*Tc[cell])),9.);
            else chi[cell] = (0.8*1.4/660.*Tc[cell])/(4*PI)/(T-Tc[cell]);

            //converts to SI
            chi[cell] = (chi[cell]*9.54393845712027); // (Tesla)
            //normalises
            chi[cell] = chi[cell]/chi_0;
            //returns the anisotropy constant
            chi[cell] = -ms[cell]*chi[cell]/ku[cell];
         }
         return chi;            //returns the 1D vector containing the perpencicular susceptability for each cell,
      }
   } //closes the internal namspace
}  //closes the micromagnetic namespace
