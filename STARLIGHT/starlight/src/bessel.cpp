///////////////////////////////////////////////////////////////////////////
//
//    Copyright 2010
//
//    This file is part of starlight.
//
//    starlight is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    starlight is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with starlight. If not, see <http://www.gnu.org/licenses/>.
//
///////////////////////////////////////////////////////////////////////////
//
// File and Version Information:
// $Rev:: 28                          $: revision of last commit
// $Author:: bgrube                   $: author of last commit
// $Date:: 2010-12-10 19:30:01 +0100 #$: date of last commit
//
// Description:
//    Bessel functions taken from ROOT
//
//
///////////////////////////////////////////////////////////////////////////


#include <iostream>
#include <fstream>
#include <cmath>

#include "bessel.h"


using namespace std;


//______________________________________________________________________________
double bessel::besI0(double x)
{
  //FROM ROOT...
  // Parameters of the polynomial approximation
  const double p1=1.0,          p2=3.5156229,    p3=3.0899424,
    p4=1.2067492,    p5=0.2659732,    p6=3.60768e-2,  p7=4.5813e-3;
  
  const double q1= 0.39894228,  q2= 1.328592e-2, q3= 2.25319e-3,
    q4=-1.57565e-3,  q5= 9.16281e-3,  q6=-2.057706e-2,
    q7= 2.635537e-2, q8=-1.647633e-2, q9= 3.92377e-3;
  
  const double k1 = 3.75;
  double ax = fabs(x);
  
  double y=0., result=0.;
  
  if (ax < k1) {
    double xx = x/k1;
    y = xx*xx;
    result = p1+y*(p2+y*(p3+y*(p4+y*(p5+y*(p6+y*p7)))));
  } else {
    y = k1/ax;
    result = (exp(ax)/sqrt(ax))*(q1+y*(q2+y*(q3+y*(q4+y*(q5+y*(q6+y*(q7+y*(q8+y*q9))))))));
  }
  return result;
}


//______________________________________________________________________________
double bessel::dbesk0(double x)
{
   // Compute the modified Bessel function K_0(x) for positive real x.  //should be k0?
   //
   //  M.Abramowitz and I.A.Stegun, Handbook of Mathematical Functions,
   //     Applied Mathematics Series vol. 55 (1964), Washington.
   //
   //--- NvE 12-mar-2000 UU-SAP Utrecht

   // Parameters of the polynomial approximation
  const double p1= -0.57721566,   p2= 0.42278420,  p3=0.23069756,
               p4=3.488590e-2,  p5=2.62698e-3, p6=1.0750e-4,  p7=7.4e-6;
   
   const double q1= 1.25331414,  q2= -7.832358e-2,  q3=2.189568e-2,
                q4= -1.062446e-2, q5=5.87872e-3,  q6= -2.51540e-3,  q7= 5.3208e-4;

   if (x <= 0) {
     cout << "BesselK0 *K0* Invalid argument x = " << x << endl;  //Should be k0?
     return 0;
   }
   
   double y=0.,result=0.;
   
   if (x <= 2) {
     y = x*x/4.;
     result = (-log(x/2.)*bessel::besI0(x))+(p1+y*(p2+y*(p3+y*(p4+y*(p5+y*(p6+y*p7))))));
   } else {
     y = 2./x;
     result = (exp(-x)/sqrt(x))*(q1+y*(q2+y*(q3+y*(q4+y*(q5+y*(q6+y*q7))))));
   }
   return result;	
}


//______________________________________________________________________________
double bessel::besI1(double x)
{
   // Compute the modified Bessel function I_1(x) for any real x.
   //
   //  M.Abramowitz and I.A.Stegun, Handbook of Mathematical Functions,
   //     Applied Mathematics Series vol. 55 (1964), Washington.
   //
   //--- NvE 12-mar-2000 UU-SAP Utrecht

   // Parameters of the polynomial approximation
   const double p1=0.5,          p2=0.87890594,   p3=0.51498869,
                p4=0.15084934,   p5=2.658733e-2,  p6=3.01532e-3,  p7=3.2411e-4;

   const double q1= 0.39894228,  q2=-3.988024e-2, q3=-3.62018e-3,
                  q4= 1.63801e-3,  q5=-1.031555e-2, q6= 2.282967e-2,
                  q7=-2.895312e-2, q8= 1.787654e-2, q9=-4.20059e-3;

   const double k1 = 3.75;
   double ax = fabs(x);

   double y=0., result=0.;
   
   if (ax < k1) {
     double xx = x/k1;
     y = xx*xx;
     result = x*(p1+y*(p2+y*(p3+y*(p4+y*(p5+y*(p6+y*p7))))));
   } else {
     y = k1/ax;
     result = (exp(ax)/sqrt(ax))*(q1+y*(q2+y*(q3+y*(q4+y*(q5+y*(q6+y*(q7+y*(q8+y*q9))))))));
     if (x < 0) result = -result;
   }
   return result;
}


//______________________________________________________________________________
double bessel::dbesk1(double x)
{
   // Compute the modified Bessel function K_1(x) for positive real x.
   //
   //  M.Abramowitz and I.A.Stegun, Handbook of Mathematical Functions,
   //     Applied Mathematics Series vol. 55 (1964), Washington.
   //
   //--- NvE 12-mar-2000 UU-SAP Utrecht

   // Parameters of the polynomial approximation
   const double p1= 1.,          p2= 0.15443144,  p3=-0.67278579,
                p4=-0.18156897,  p5=-1.919402e-2, p6=-1.10404e-3,  p7=-4.686e-5;
   
   const double q1= 1.25331414,  q2= 0.23498619,  q3=-3.655620e-2,
                q4= 1.504268e-2, q5=-7.80353e-3,  q6= 3.25614e-3,  q7=-6.8245e-4;

   if (x <= 0) {
     cout << "bessel:dbesk1 *K1* Invalid argument x = " << x << endl;
     return 0;
   }
   
   double y=0.,result=0.;

   if (x <= 2) {
     y = x*x/4.;
     result = (log(x/2.)*bessel::besI1(x))+(1./x)*(p1+y*(p2+y*(p3+y*(p4+y*(p5+y*(p6+y*p7))))));
   } else {
     y = 2./x;
     result = (exp(-x)/sqrt(x))*(q1+y*(q2+y*(q3+y*(q4+y*(q5+y*(q6+y*q7))))));
   }
   return result;
}
