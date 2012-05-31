//--------------------------------------------------------------------------
//
// Environment:
//      This software is part of the EvtGen package developed jointly
//      for the BaBar and CLEO collaborations.  If you use all or part
//      of it, please give an appropriate acknowledgement.
//
// Copyright Information: See EvtGen/COPYRIGHT
//      Copyright (C) 2000      Caltech, UCSB
//
// Module: 
// Description: Form factors for b->sll according to Ali, Ball et al.
//              hep-ph/9910221v2
//
// Modification history:
//
//    Ryd     January 5, 2000         Module created
//
//    jjhollar October 7 2005         Add Ball & Zwicky '05 LCSR
//    ofte     August 2006            Add some more pi l l FF models
//------------------------------------------------------------------------
#include "EvtGenBase/EvtPatches.hh"

#include "EvtGenBase/EvtPatches.hh"
#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenModels/EvtbTosllBallFF.hh"
#include <math.h>

EvtbTosllBallFF::EvtbTosllBallFF(int ffmodel)
{
  _theFFModel = ffmodel;
}


void EvtbTosllBallFF::getScalarFF(EvtId parent, EvtId daught,
				  double t, double /*mass*/, 
				  double& fp,double& f0,double& ft){

    int model = _theFFModel;

  double m=EvtPDL::getMeanMass(parent);
  double md=EvtPDL::getMeanMass(daught);
  
  double shat=t/(m*m);
  double shat2=shat*shat;
  double shat3=shat2*shat;

  if (daught == EvtPDL::getId(std::string("K+")) ||
      daught == EvtPDL::getId(std::string("K-")) ||
      daught == EvtPDL::getId(std::string("K_S0")) ||
      daught == EvtPDL::getId(std::string("K0")) ||
      daught == EvtPDL::getId(std::string("anti-K0")) ||
      daught == EvtPDL::getId(std::string("K_L0"))
     ) 
 {
   // B --> K form factors
  if (model == 1) {
        //this is Ali-Ball '01 (or really Ali-Ball'99 minimum allowed)
    fp = 0.278*exp(1.568*shat+0.470*shat2+0.885*shat3);
    f0 = 0.278*exp(0.740*shat+0.080*shat2+0.425*shat3);
    ft = 0.300*exp(1.600*shat+0.501*shat2+0.796*shat3);
  }
  if (model == 2) {
        //this is Ali-Ball '99 (central values)
    fp = 0.319*exp(1.465*shat+0.372*shat2+0.782*shat3);
    f0 = 0.319*exp(0.633*shat-0.095*shat2+0.591*shat3);
    ft = 0.355*exp(1.478*shat+0.373*shat2+0.700*shat3);
  }
  if (model == 3) {
        //QCD sum rules (Colangelo et al)
    fp = 0.25/(1.-t/(5.0*5.0));
    f0 = 0.25/(1.-t/(7.0*7.0));
    ft = - 0.14/((1.0 - t/(5.0*5.0))*(1.0 - t/(7.0*7.0)));
  } 
  if (model == 4) {
        // Quark model  (Melikhov et al - hep-ph/9711362)         
    fp = 0.36/(1. - 0.048*t + 0.00063*t*t);
    double fm = -0.30/(1. - 0.050*t + 0.00061*t*t);
    f0 = fp + fm*(t/(m*m - md*md));
    ft  = -(m+md)*0.06/(1 -0.049*t + 0.00064*t*t);
  }
  if (model == 5) {
    fp = 0.341/(1. - 1.41*shat + 0.406*shat*shat);
    f0 = 0.341/(1. - 0.41*shat -0.361*shat*shat);
    ft = 0.374/(1. - 1.42*shat + 0.434*shat*shat);
  }
  if (model == 6) {
    // Ball-Zwicky LCSR '05 "set 2" (mb = 4.8)
    fp = (0.1616 / (1. - (t/29.30))) + 
      (0.1730 / (1. - (t/29.30)) / (1. - (t/29.30)));
    f0 = (0.3302 / (1. - (t/37.46/37.46)));
    ft = (0.1614 / (1. - (t/29.30))) + 
      (0.1981 / (1. - (t/29.30)) / (1. - (t/29.30)));
  }
  if (model == 7){
    // Ball-Zwicky LCSR '05 "set 4" (mb = 4.6)
    fp = (0.1903 / (1. - (t/29.30))) + 
      (0.1478 / (1. - (t/29.30)) / (1. - (t/29.30)));
    f0 = (0.3338 / (1. - (t/38.98/38.98)));
    ft = (0.1851 / (1. - (t/29.30))) + 
      (0.1905 / (1. - (t/29.30)) / (1. - (t/29.30)));
  }
 }
 else if (daught == EvtPDL::getId(std::string("pi+")) ||
	  daught == EvtPDL::getId(std::string("pi-")) ||
	  daught == EvtPDL::getId(std::string("pi0"))
	  ) 
  {
    if(model == 1){
      // B --> pi form factors from Ball-Zwicky'01 (tabulated in hep-ph/0306251)
      fp = 0.261/(1. - 2.03*shat + 1.293*shat*shat);
      f0 = 0.261/(1. - 0.27*shat  -0.752*shat*shat);
      ft = 0.296/(1. - 1.28*shat + 0.193*shat*shat);
    }

    // The following two (2) and (3) should preferably be replaced with
    // something better. ft not provided by the papers, and the ft formula 
    // (from Colangelo'96 (hep-ph/9510403v2) equation (5.1)) seem to be no good. 
    //  if (model == 2) {
    //    // LatticeQCD: Okamoto'04
    //    // ... with f_T from eq. (5.1) of Colangelo'95
    //    fp = 0.23/((1.-shat)*(1.-0.63*shat));
    //    f0 = 0.23/((1.-shat)/1.18);
    //    ft = ( -(m+md)/(2*m) )/ (fp - (m*m-md*md) * ((f0-fp)/t) );
    //  }
    //  if (model == 3) {
    //    // LatticeQCD: Shigemitsu'04
    //    // ... with f_T from eq. (5.1) of Colangelo'95
    //    fp = 0.42*(1.0-0.41)/((1.-shat)*(1.-0.63*shat));
    //    f0 = 0.42*(1.0-0.41)/((1.-shat)/1.18);
    //    ft = ( -(m+md)/(2*m) )/ (fp - (m*m-md*md) * ((f0-fp)/t) );
    //   }

    if (model == 4) {
      // Quark model - B -> pi set 1 (Melikhov-Nikitin'96) 
      fp =  0.29/pow((1.-t/(6.48*6.48)),2.54);
      double fm = -0.26/pow((1.-t/(6.34*6.34)),2.49);
      f0 = fp + fm*(t/(m*m - md*md));
      ft  = -(m+md)*0.05/pow((1.-t/(6.47*6.47)),2.50);
    }
    if (model == 5) {
      // Melikhov-Stech '00. (hep-ph/0001113)
      // relativistic dispersion approach based on constituent quark picture
      fp = 0.29/((1.-shat)*(1.-0.48*shat));
      f0 = 0.29/(1.-0.76*shat +0.28*shat*shat);
      ft = 0.28/((1.-shat)*(1.-0.48*shat));
    }
    if(model == 6){
      // Ball-Zwicky LCSR '05 "set 2" (mb = 4.8)
      fp = (0.744 / (1. - (t/(5.32*5.32)))) + (-0.486 / (1. - (t/40.73)));
      f0 = (0.258 / (1. - (t/33.81)));
      ft = (1.387 / (1. - (t/(5.32*5.32)))) + (-1.134 / (1. - (t/32.22)));
    }
    if(model == 7){
      // Ball-Zwicky LCSR '05 "set 4" (mb = 4.6)
      fp = (0.944 / (1. - (t/(5.32*5.32)))) + (-0.669 / (1. - (t/34.27)));
      f0 = (0.270 / (1. - (t/33.63)));
      ft = (0.152 / (1. - (t/(5.32*5.32)))) + 
	(0.122 / (1. - (t/28.40)) / (1. - (t/28.40)));
    }
  }
 else if (daught == EvtPDL::getId(std::string("eta")) ||
	  daught == EvtPDL::getId(std::string("eta'"))
	  )
 {
   if(model == 1){
     // B --> eta form factors
     fp = 0.261/(1. - 2.03*shat + 1.293*shat*shat);
     f0 = 0.261/(1. - 0.27*shat  -0.752*shat*shat);
     ft = 0.296/(1. - 1.28*shat + 0.193*shat*shat);
   }
   if(model == 6){
    // Ball-Zwicky LCSR '05 "set 2" (mb = 4.8)
    fp = (0.1220 / (1. - (t/28.40))) + 
      (0.1553 / (1. - (t/28.40)) / (1. - (t/28.40)));
    f0 = (0.2734 / (1. - (t/31.03/31.03)));
    ft = (0.1108 / (1. - (t/28.40))) + 
      (0.1752 / (1. - (t/28.40)) / (1. - (t/28.40)));
   }
   if(model == 7){
    // Ball-Zwicky LCSR '05 "set 2" (mb = 4.8)
    fp = (0.1380 / (1. - (t/28.40))) + 
      (0.1462 / (1. - (t/28.40)) / (1. - (t/28.40)));
    f0 = (0.2799 / (1. - (t/30.46/30.46)));
    ft = (0.1160 / (1. - (t/28.40))) + 
      (0.1841 / (1. - (t/28.40)) / (1. - (t/28.40)));
   }
 }

  // cout << "shat "<<shat<<"\t"<<"fp "<<fp<<"\t"<<"f0 "<<f0<<"\t"
  //    <<"ft "<<ft<<endl;

}


void EvtbTosllBallFF::getVectorFF(EvtId parent, EvtId daught,
				  double t, double /*mass*/, 
				  double& a1,double& a2,double& a0, double& v,
				  double& t1, double& t2, double& t3 ){

  int model = _theFFModel;
  
  double m=EvtPDL::getMeanMass(parent);
  double md=EvtPDL::getMeanMass(daught);
  
  double shat=t/(m*m);
  double shat2=shat*shat;

  if (
      daught == EvtPDL::getId(std::string("K*+")) ||
      daught == EvtPDL::getId(std::string("K*-")) ||
      daught == EvtPDL::getId(std::string("K*0")) ||
      daught == EvtPDL::getId(std::string("anti-K*0"))
     ) 
 {
  if (model == 1) {
        //this is Ali-Ball '01
    a1 = 0.294*exp(0.656*shat+0.456*shat2);
    a2 = 0.246*exp(1.237*shat+0.822*shat2);
    a0 = 0.412*exp(1.543*shat+0.954*shat2);
    v  = 0.399*exp(1.537*shat+1.123*shat2);
    
    t1 = 0.334*exp(1.575*shat+1.140*shat2);
    t2 = 0.334*exp(0.562*shat+0.481*shat2);
    t3 = 0.234*exp(1.230*shat+1.089*shat2);
  }
  if (model == 2) {
        //this is Ali-Ball '99
    a1=0.337*exp(0.602*shat+0.258*shat2);
    a2=0.282*exp(1.172*shat+0.567*shat2);
    a0=0.471*exp(1.505*shat+0.710*shat2);
    v=0.457*exp(1.482*shat+1.015*shat2);
    
    t1=0.379*exp(1.519*shat+1.030*shat2);
    t2=0.379*exp(0.517*shat+0.426*shat2);
    t3=0.260*exp(1.129*shat+1.128*shat2);
  }
  if (model == 3) {
    //QCD sum rules (Colangelo et al)
    // JJH - Changed in accordance with erratum (Phys. Rev. D 57, 3186 (1998))
    a1 = 0.37*(1 - 0.023*t);
    a2 = 0.40*(1 + 0.034*t);
    a0 = 0.3/(1.- t/(4.8*4.8));
    v = 0.47/(1.- t/(5.0*5.0));
    
    t1 = 0.19/(1.-t/(5.3*5.3));
    t2 = 0.19*(1. - 0.02*t);
    t3 = 0.3*(1. + 0.01*t); 
  }

  if (model == 4) {
        // Quark model  (Melikhov et al)         
    a1 = 1.6/(1 - 0.0288*t + 0.00028*t*t); a1 = a1/(m+md);
    a2 = (m+md)*0.036/(1. - 0.053*t + 0.00082*t*t);
    double aminus = 0.041/(1. - 0.055*t + 0.00088*t*t);
    double f =  1.60/(1. - 0.0288*t + 0.00028*t*t);
    double aplus = -0.036/(1. - 0.053*t + 0.00082*t*t);
    a0 = (t*aminus + f + (m*m-md*md)*aplus)/(2.0*md);
    v = (m+md)*0.048/(1. - 0.057*t + 0.00085*t*t);

    t1 = 0.28/(1. - 0.058*t + 0.0009*t*t);
    double gplus = -0.28/(1. - 0.058*t + 0.0009*t*t);
    double gminus = 0.24/(1. - 0.059*t + 0.00096*t*t);
    t2 = -gplus - (t*gminus)/(m*m-md*md);
    double h = 0.0037/(1. - 0.075*t + 0.0016*t*t);
    t3 = (m+md)*(m+md)*((gminus/(m*m-md*md) - h/2.));
        
  }
  if (model == 5) {
    a1 = 0.337/(1. - 0.60*shat - 0.023*shat*shat);
    a2 = 0.283/(1. - 1.18*shat + 0.281*shat*shat);
    a0 = 0.470/(1. - 1.55*shat + 0.680*shat*shat);
    v  = 0.458/(1. - 1.55*shat + 0.575*shat*shat);
    t1 = 0.379/(1. - 1.59*shat + 0.615*shat*shat);
    t2 = 0.379/(1. - 0.49*shat - 0.241*shat*shat);
    t3 = 0.261/(1. - 1.20*shat + 0.098*shat*shat);
  }
  if (model == 6) {
    // Ball-Zwicky LCSR '05 (mb = 4.8)
    a1 = 0.290 / (1. - (t/40.38));
    a2 = (-0.084 / (1. - (t/52.00))) + 
      (0.342 / (1. - (t/52.00))/(1. - (t/52.00)));
    a0 = (1.364 / (1. - (t/(5.28*5.28)))) + 
      (-0.990 / (1. - (t/36.78)));
    v = (0.923 / (1. - (t/(5.32*5.32)))) + 
      (-0.511 / (1. - (t/49.40)));
    t1 = (0.823 / (1. - (t/(5.32*5.32)))) + 
      (-0.491 / (1 - (t/46.31)));
    t2 = 0.333 / (1. - (t/41.41));
    t3 = (-0.036 / (1. - (t/48.10))) + 
      (0.368 / (1. - (t/48.10)) / (1. - (t/48.10)));
  }
 }
 else if (daught == EvtPDL::getId(std::string("rho+")) ||
	  daught == EvtPDL::getId(std::string("rho-")) ||
	  daught == EvtPDL::getId(std::string("rho0")) 
	  )
   {
     if(model == 1){
       // B --> rho form factors
       a1 = 0.261/(1. - 0.29*shat - 0.415*shat*shat);
       a2 = 0.223/(1. - 0.93*shat - 0.092*shat*shat);
       a0 = 0.372/(1. - 1.40*shat + 0.437*shat*shat);
       v  = 0.338/(1. - 1.37*shat + 0.315*shat*shat);
       t1 = 0.285/(1. - 1.41*shat + 0.361*shat*shat);
       t2 = 0.285/(1. - 0.28*shat - 0.500*shat*shat);
       t3 = 0.202/(1. - 1.06*shat - 0.076*shat*shat);
     }
     if(model == 6){
       // Ball-Zwicky LCSR '05 (mb = 4.8)
       a1 = 0.240 / (1. - (t/37.51));
       a2 = (0.009 / (1. - (t/40.82))) + 
	 (0.212 / (1. - (t/40.82)) / (1. - (t/40.82)));
       a0 = (1.527 / (1. - (t/(5.28*5.28)))) + 
	 (-1.220 / (1. - (t/33.36)));
       v = (1.045 / (1. - (t/(5.32*5.32)))) + 
	 (-0.721 / (1. - (t/38.34)));
       t1 = (0.897 / (1. - (t/(5.32*5.32)))) + 
	 (-0.629 / (1 - (t/38.04)));
       t2 = 0.267 / (1. - (t/38.59));
       t3 = (0.022 / (1. - (t/40.88))) + 
	 (0.246 / (1. - (t/40.88)) / (1. - (t/40.88)));
     }
   }
  else if (daught == EvtPDL::getId(std::string("omega"))
	   ) 
    {
     if(model == 1){
       // B --> rho form factors
       a1 = 0.261/(1. - 0.29*shat - 0.415*shat*shat);
       a2 = 0.223/(1. - 0.93*shat - 0.092*shat*shat);
       a0 = 0.372/(1. - 1.40*shat + 0.437*shat*shat);
       v  = 0.338/(1. - 1.37*shat + 0.315*shat*shat);
       t1 = 0.285/(1. - 1.41*shat + 0.361*shat*shat);
       t2 = 0.285/(1. - 0.28*shat - 0.500*shat*shat);
       t3 = 0.202/(1. - 1.06*shat - 0.076*shat*shat);
     }
     if(model == 6){
       // Ball-Zwicky LCSR '05 (mb = 4.8)
       a1 = -0.217 / (1. - (t/37.01));
       a2 = (0.006 / (1. - (t/41.24))) + 
	 (0.192 / (1. - (t/41.24)) / (1. - (t/41.24)));
       a0 = (1.321 / (1. - (t/(5.28*5.28)))) + 
	 (-1.040 / (1. - (t/34.47)));
       v = (1.006 / (1. - (t/(5.32*5.32)))) + 
	 (-0.713 / (1. - (t/37.45)));
       t1 = (0.865 / (1. - (t/(5.32*5.32)))) + 
	 (-0.622 / (1 - (t/37.19)));
       t2 = 0.242 / (1. - (t/37.95));
       t3 = (0.023 / (1. - (t/40.87))) + 
	 (0.220 / (1. - (t/40.87)) / (1. - (t/40.87)));
     }
    }
  else if (daught == EvtPDL::getId(std::string("phi"))
	   ) 
    {
     if(model == 6){
       // Ball-Zwicky LCSR '05 (mb = 4.8)
       a1 = 0.308 / (1. - (t/36.54));
       a2 = (-0.054 / (1. - (t/48.94))) + 
	 (0.288 / (1. - (t/48.94)) / (1. - (t/48.94)));
       a0 = (3.310 / (1. - (t/(5.28*5.28)))) + 
	 (-2.835 / (1. - (t/31.57)));
       v = (1.484 / (1. - (t/(5.32*5.32)))) + 
	 (-1.049 / (1. - (t/39.52)));
       t1 = (1.303 / (1. - (t/(5.32*5.32)))) + 
	 (-0.954 / (1 - (t/38.28)));
       t2 = 0.349 / (1. - (t/37.21));
       t3 = (0.027 / (1. - (t/45.56))) + 
	 (0.321 / (1. - (t/45.56)) / (1. - (t/45.56)));
     }
    }

}







