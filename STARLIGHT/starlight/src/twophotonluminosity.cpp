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
// $Rev:: 293                         $: revision of last commit
// $Author:: butter                   $: author of last commit
// $Date:: 2017-11-11 15:46:05 +0100 #$: date of last commit
//
// Description:
//    Added incoherent factor to luminosity table output--Joey
//    Added BRANGE method to integrate from bmin to bmax, irrespective of R_A, etc.  Spencer July, 2017  
//
///////////////////////////////////////////////////////////////////////////


#include <iostream>
#include <fstream>
#include <cmath>

#include "inputParameters.h"
#include "beambeamsystem.h"
#include "beam.h"
#include "starlightconstants.h"
#include "nucleus.h"
#include "bessel.h"
#include "twophotonluminosity.h"

using namespace std;
using namespace starlightConstants;



//______________________________________________________________________________
twoPhotonLuminosity::twoPhotonLuminosity(const inputParameters& inputParametersInstance, beam beam_1,beam beam_2):
beamBeamSystem(inputParametersInstance, beam_1,beam_2)
,_gamma(inputParametersInstance.beamLorentzGamma())
,_nWbins(inputParametersInstance.nmbWBins())
,_nYbins(inputParametersInstance.nmbRapidityBins())
,_wMin(inputParametersInstance.minW())
,_yMin(-inputParametersInstance.maxRapidity())
,_wMax(inputParametersInstance.maxW())
,_yMax(inputParametersInstance.maxRapidity())
,_productionMode(inputParametersInstance.productionMode())
,_beamBreakupMode(inputParametersInstance.beamBreakupMode())
,_interferenceEnabled(inputParametersInstance.interferenceEnabled())
,_interferenceStrength(inputParametersInstance.interferenceStrength())
,_maxPtInterference(inputParametersInstance.maxPtInterference())
,_nmbPtBinsInterference(inputParametersInstance.nmbPtBinsInterference())
,_xsecCalcMethod(inputParametersInstance.xsecCalcMethod())
,_baseFileName(inputParametersInstance.baseFileName())
,_bmin(inputParametersInstance.bmin())
,_bmax(inputParametersInstance.bmax())
{
  //Lets check to see if we need to recalculate the luminosity tables
  double bmina = _bmin;
  double bmaxa = _bmax;
  twoPhotonDifferentialLuminosity(bmina,bmaxa);
}


//______________________________________________________________________________
twoPhotonLuminosity::~twoPhotonLuminosity()
{ }


//______________________________________________________________________________

void twoPhotonLuminosity::twoPhotonDifferentialLuminosity(const double _bmin, const double _bmax)
{
  std::string wyFileName;
  wyFileName = _baseFileName +".txt";

  ofstream wylumfile;
  wylumfile.precision(15);
  wylumfile.open(wyFileName.c_str());
  std::vector<double> w(_nWbins);
  std::vector<double> y(_nYbins);
  double xlum = 0.; 
  double Normalize = 0.,OldNorm;
  double wmev = 0;
 
  Normalize = 1./sqrt(1*(double)_nWbins*_nYbins); //if your grid is very fine, you'll want high accuracy-->small Normalize
  OldNorm   = Normalize;
  
  //Writing out our input parameters+(w,y)grid+diff._lum.
  wylumfile << beam1().Z() <<endl;
  wylumfile << beam1().A() <<endl;
  wylumfile << beam2().Z() <<endl;
  wylumfile << beam2().A() <<endl;
  wylumfile << _gamma <<endl;
  wylumfile << _wMax <<endl;
  wylumfile << _wMin <<endl;
  wylumfile << _nWbins <<endl;
  wylumfile << _yMax <<endl;
  wylumfile << _nYbins <<endl;
  wylumfile << _productionMode <<endl;
  wylumfile << _beamBreakupMode <<endl;
  wylumfile << _interferenceEnabled <<endl;
  wylumfile << _interferenceStrength <<endl;
  wylumfile << _ip->deuteronSlopePar() <<endl;
  wylumfile << _maxPtInterference <<endl;
  wylumfile << _nmbPtBinsInterference <<endl;
  for (unsigned int i = 0; i < _nWbins; i++) {
    w[i] = _wMin + (_wMax-_wMin)/(_nWbins-1)*i;
    wylumfile << w[i] <<endl;
  }
  for (unsigned int i = 0; i < _nYbins; i++) {
    y[i] = -_yMax + 2.*_yMax*i/(_nYbins-1);
    wylumfile << y[i] <<endl;
  }

  if (_beamBreakupMode ==8) { // new method - integrate over fixed impact parameter range, regardless of nuclear breakup SRK 7/2017
    
    cout <<"Calculating two photon luminosity from "<<_bmin<<" fm to "<<_bmax<<" fm."<<endl;
       for (unsigned int i = 0; i < _nWbins; i++) {   //For each (w,y) pair, calculate the diff. _lum
      printf("Calculating cross section: %2.0f %% \r", float(i)/float(_nWbins)*100);
      fflush(stdout);
      for (unsigned int j = 0; j < _nYbins; j++) {
        xlum = w[i] * D2LDMDYBRANGE(w[i],y[j],_bmin,_bmax);   //Convert photon flux dN/dW to Lorentz invariant photon number WdN/dW
        wylumfile << xlum <<endl;
        // cout<<" i: "<<i<<" j: "<<j<<" W*dN/dW: "<<xlum<<endl; 
      }
    }
    return;}  // Just return here to simplify program structure
    

  if(_xsecCalcMethod == 0) {  // faster (analytic) method, only works for symmetric collisions
    
    for (unsigned int i = 0; i < _nWbins; ++i) {   //For each (w,y) pair, calculate the diff. _lum
      for (unsigned int j = 0; j < _nYbins; ++j) {
        wmev = w[i]*1000.;
        xlum = wmev * D2LDMDY(wmev,y[j],Normalize);   //Convert photon flux dN/dW to Lorentz invariant photon number WdN/dW
        if (j==0) OldNorm = Normalize;       //Save value of integral for each new W(i) and Y(i)
        wylumfile << xlum <<endl;
      }
      Normalize = OldNorm;
    }

  }
  else if(_xsecCalcMethod == 1) {  // slower method, always works for any colliding system

    for (unsigned int i = 0; i < _nWbins; i++) {   //For each (w,y) pair, calculate the diff. _lum
      printf("Calculating cross section: %2.0f %% \r", float(i)/float(_nWbins)*100);
      fflush(stdout);
      for (unsigned int j = 0; j < _nYbins; j++) {
        xlum = w[i] * D2LDMDY(w[i],y[j]);   //Convert photon flux dN/dW to Lorentz invariant photon number WdN/dW
        wylumfile << xlum <<endl;
        // cout<<" i: "<<i<<" j: "<<j<<" W*dN/dW: "<<xlum<<endl; 
      }
    }
  }
    return;
}

//______________________________________________________________________________
double twoPhotonLuminosity::D2LDMDY(double M,double Y,double &Normalize)
{
  // double differential luminosity using faster, analytic method (symmetric systems only)

  double D2LDMDYx = 0.;

  _W1    =  M/2.0*exp(Y);
  _W2    =  M/2.0*exp(-Y);
  int Zin1=beam1().Z();
  int Zin2=beam2().Z();
  D2LDMDYx = 2.0/M*Zin1*Zin1*Zin2*Zin2*(starlightConstants::alpha*starlightConstants::alpha)*integral(Normalize); 
  Normalize = D2LDMDYx*M/(2.0*Zin1*Zin1*Zin2*Zin2*
                          starlightConstants::alpha*starlightConstants::alpha);
  return D2LDMDYx;
}



//______________________________________________________________________________
double twoPhotonLuminosity::D2LDMDY(double M, double Y) const 
{
  // double differential luminosity using slower numerical method for standard case

  double D2LDMDYx = 0.;
  double w1    =  M/2.0*exp(Y);
  double w2    =  M/2.0*exp(-Y);
  
  double r_nuc1 = beam1().nuclearRadius();
  double r_nuc2 = beam2().nuclearRadius();
  
  double b1min = r_nuc1;
  double b2min = r_nuc2;
  
  double b1max = max(5.*_gamma*hbarc/w1,5*r_nuc1);
  double b2max = max(5.*_gamma*hbarc/w2,5*r_nuc2);
  
  const int nbins_b1 = 120;
  const int nbins_b2 = 120;
  
  double log_delta_b1 = (log(b1max)-log(b1min))/nbins_b1;
  double log_delta_b2 = (log(b2max)-log(b2min))/nbins_b2;
  double sum = 0;
  for(int i = 0; i < nbins_b1; ++i)
  {
      // Sum from nested integral
      double sum_b2 = 0;
      double b1_low = b1min*exp(i*log_delta_b1);
      double b1_high = b1min*exp((i+1)*log_delta_b1);
      double b1_cent = (b1_high+b1_low)/2.;
      for(int j = 0; j < nbins_b2; ++j)
      {
	// Sum from nested  
	double sum_phi = 0;
	double b2_low = b2min*exp(j*log_delta_b2);
	double b2_high = b2min*exp((j+1)*log_delta_b2);
	double b2_cent = (b2_high+b2_low)/2.;
	
	// Gaussian integration n = 10
	// Since cos is symmetric around 0 we only need 5 of the 
	// points in the gaussian integration.
	// these are points 1, 3, 5, 7, 9 in the standard n=10 Gaussian weights table
	const int ngi = 5;
	double weights[ngi] = 
	{
	  0.2955242247147529,
	  0.2692667193099963,
	  0.2190863625159820,
	  0.1494513491505806,
	  0.0666713443086881,
	};
	
	double abscissas[ngi] =
	{
	  -0.1488743389816312,
	  -0.4333953941292472,
	  -0.6794095682990244,
	  -0.8650633666889845,
	  -0.9739065285171717,
	};
	
	for(int k = 0; k < ngi; ++k)
	{
	    double b_rel = sqrt(b1_cent*b1_cent+b2_cent*b2_cent + 2.*b1_cent*b2_cent*cos(pi*(abscissas[k]+1)));
	    
	    sum_phi += weights[k] * probabilityOfBreakup(b_rel)*2;
	}
	sum_b2 += beam2().photonDensity(b2_cent,w2)*pi*sum_phi*b2_cent*(b2_high-b2_low);
      }
      
      sum += beam1().photonDensity(b1_cent, w1)*sum_b2*b1_cent*(b1_high-b1_low);
      
  }
  D2LDMDYx = 2.*pi*M/2.*sum;
  return D2LDMDYx; 
}

/// new routine goes here

//______________________________________________________________________________
// Version of twoPhoton  Luminosity with a fixed impact parameter range from bmin to bmax  SRK August, 2017
double twoPhotonLuminosity::D2LDMDYBRANGE(double M, double Y, double bmin, double bmax) const 
{
  // double differential luminosity using slower numerical method

  double D2LDMDYx = 0.;
  double w1    =  M/2.0*exp(Y);
  double w2    =  M/2.0*exp(-Y);
  
  double r_nuc1 = beam1().nuclearRadius();
  double r_nuc2 = beam2().nuclearRadius();
  
  double b1min = r_nuc1;
  double b2min = r_nuc2;
  
  double b1max = max(5.*_gamma*hbarc/w1,5*r_nuc1);
  double b2max = max(5.*_gamma*hbarc/w2,5*r_nuc2);
  
  const int nbins_b1 = 120;
  const int nbins_b2 = 120;
  
  double log_delta_b1 = (log(b1max)-log(b1min))/nbins_b1;
  double log_delta_b2 = (log(b2max)-log(b2min))/nbins_b2;
  double sum = 0;
  for(int i = 0; i < nbins_b1; ++i)
  {
      // Sum from nested integral
      double sum_b2 = 0;
      double b1_low = b1min*exp(i*log_delta_b1);
      double b1_high = b1min*exp((i+1)*log_delta_b1);
      double b1_cent = (b1_high+b1_low)/2.;
      for(int j = 0; j < nbins_b2; ++j)
      {
	// Sum from nested  
	double sum_phi = 0;
	double b2_low = b2min*exp(j*log_delta_b2);
	double b2_high = b2min*exp((j+1)*log_delta_b2);
	double b2_cent = (b2_high+b2_low)/2.;

// integrate over ion-ion separation B from bmin to bmax; this corresponds to an angular range.
// This follows Eq. 2.4 of G. Baur et al., Nucl. Phys. A518, 786 (1990)
// cos(theta) = (B^2 -b_1^2-b_2^2)/2b_1b_2
// We don't care about nuclear breakup, so do the delta phi integral analytically; it's just delta phi
// integration from bmin to bmax is equivalent to integrating from thetamin to thetamax
	
// be careful - some choices of b_1,b_2, bmin and/or bmax, do not form a valid triangle
// We take care of this by adjustin the limits, restricting |cos(thetamax|<1, ditto for thetamin.  
	
	    double costhetamin=(b1_cent*b1_cent+b2_cent*b2_cent-bmin*bmin)/(2.*b1_cent*b2_cent);
	    double costhetamax=(b1_cent*b1_cent+b2_cent*b2_cent-bmax*bmax)/(2.*b1_cent*b2_cent);
	    if (costhetamin>1.) {costhetamin=1.;}
	    if (costhetamax>1.) {costhetamax=1.;}
	    if (costhetamin<-1.){costhetamin=-1.;}
            if (costhetamax<-1.){costhetamax=-1.;}

	    double thetamin=acos(costhetamin);
	    double thetamax =acos(costhetamax);
	      
	    double deltatheta=thetamax-thetamin;
	    if (deltatheta < 0. || bmin > bmax || abs(deltatheta)>pi) {
	          cout<< "Error: angle < 0. deltatheta= "<<deltatheta<<"  bmin, bmax="<<bmin<<bmax<<endl;
		  cout <<"thetamin, thetamax="<<thetamin<<" " << thetamax<<endl;
		  cout <<"b1_cent, b2_cent="<<b1_cent<<" "<<b2_cent<<endl;
		  exit(0);
	    }

// Since we do not care about the probability of nuclear breakup here, the result of the integral over phi is just delta phi.
// In the standard version of the code, one uses the Gaussian integration weights, but they just sum to 1.
// the 2 is because the integral only covers half of the phase space

	      sum_phi=2.*deltatheta;
	      sum_b2 += beam2().photonDensity(b2_cent,w2)*sum_phi*b2_cent*(b2_high-b2_low);
	    }
      
	    sum += beam1().photonDensity(b1_cent, w1)*sum_b2*b1_cent*(b1_high-b1_low);
      
      }
  D2LDMDYx = 2.*sum*pi*M/2;

  if (D2LDMDYx == 0.){
    cout<<"error - returning luminosity zero.  bmin,bmax="<<bmin<<bmax<<endl;
  }
  return D2LDMDYx; 
}


void * twoPhotonLuminosity::D2LDMDY_Threaded(void * a)
{
  difflumiargs *args = (difflumiargs*)a;
  double M = args->m;
  double Y = args->y;
  args->res = args->self->D2LDMDY(M, Y);
  
  return NULL;
}


//______________________________________________________________________________ 
double twoPhotonLuminosity::integral(double Normalize)
{
  int NIter = 0;
  int NIterMin = 0;
  double EPS = 0.;
  // double RM = 0.;
  double RM1 = 0.;
  double RM2 = 0.;
  double u1 = 0.;
  double u2 = 0.;
  double B1 = 0.;
  double B2 = 0.;
  double Integrala = 0.;
  double totsummary = 0.;
  double NEval      = 0.;
  double Lower[3];
  double Upper[3];
  double *WK = new double[500000];
  double Result, Summary, ResErr, NFNEVL;

  EPS = .01*Normalize;   //This is EPS for integration, 1% of previous integral value.
  // Change this to the Woods-Saxon radius to be consistent with the older calculations (JN 230710) 
  RM1  = beam1().nuclearRadius()/starlightConstants::hbarcmev;  
  RM2  = beam2().nuclearRadius()/starlightConstants::hbarcmev;  
  // RM  = beam1().woodSaxonRadius()/starlightConstants::hbarcmev;  

  NIter = 10000 + (int)1000000*(int)Normalize; //if integral value is very small, we don't do too many intertions to get precision down to 1%
  NIterMin = 600;
  u1 = 9.*_gamma/_W1; //upper boundary in B1
  u2 = 9.*_gamma/_W2; //upper boundary in B2
  B1 = .4*_gamma/_W1; //intermediate boundary in B1
  B2 = .4*_gamma/_W2; //intermediate boundary in B2
  //The trick is that u1,2 and b1,2 could be less than RM-the lower integration boundary, thus integration area splits into 4,2 or 1 pieces
  
  if (u1 < RM1){
    Integrala = 0;
    totsummary = 0;
    NEval      = 0;
  }
  else if (B1 > RM1){
    if (u2 < RM2){
      Integrala = 0;
      totsummary = 0;
      NEval      = 0;
    }
    else if (B2 > RM2){            //integral has 4 parts
      Integrala = 0;
      totsummary = 40000;
      NEval      = 0;
      Lower[0]   = RM1;       //1
      Lower[1]   = RM2;       //2
      Lower[2]   = 0.;       //3
      Upper[2]   = 2.*starlightConstants::pi;    //3
      Upper[0]   = B1;       //1
      Upper[1]   = B2;       //2
      radmul(3,Lower,Upper,NIterMin,NIter,EPS,WK,NIter,Result,ResErr,NFNEVL,Summary);
      Integrala   = Integrala + Result;
      totsummary = totsummary + 1000*Summary;
      NEval      = NEval + NFNEVL;
      Upper[0]   = u1;       //1
      Upper[1]   = B2;       //2
      Lower[0]   = B1;       //1
      Lower[1]   = RM2;       //2
      radmul(3,Lower,Upper,NIterMin,NIter,EPS,WK,NIter,Result,ResErr,NFNEVL,Summary);
      Integrala   = Integrala + Result;
      totsummary = totsummary + 100*Summary;
      NEval      = NEval + NFNEVL;
      Upper[0]   = B1;       //1
      Upper[1]   = u2;       //2
      Lower[0]   = RM1;       //1
      Lower[1]   = B2;       //2
      radmul(3,Lower,Upper,NIterMin,NIter,EPS,WK,NIter,Result,ResErr,NFNEVL,Summary);
      Integrala   = Integrala + Result;
      totsummary = totsummary + 100*Summary;
      NEval      = NEval + NFNEVL;
      Upper[0]   = u1;       //1
      Upper[1]   = u2;       //2
      Lower[0]   = B1;       //1
      Lower[1]   = B2;       //2
      radmul(3,Lower,Upper,NIterMin,NIter,EPS,WK,NIter,Result,ResErr,NFNEVL,Summary);
      Integrala   = Integrala + Result;
      totsummary = totsummary + Summary;
      NEval      = NEval + NFNEVL;
    }
    else {
      //integral has 2 parts, b2 integral has only 1 component
      Integrala   = 0;
      totsummary = 20000;
      NEval      = 0;
      Lower[0]   = RM1;       //1
      Lower[1]   = RM2;       //2
      Lower[2]   = 0.;       //3
      Upper[2]   = 2.*starlightConstants::pi;    //3
      Upper[0]   = B1;       //1
      Upper[1]   = u2;       //2
      radmul(3,Lower,Upper,NIterMin,NIter,EPS,WK,NIter,Result,ResErr,NFNEVL,Summary);
      Integrala   = Integrala + Result;
      totsummary = totsummary + 100*Summary;
      NEval      = NEval + NFNEVL;
      Upper[0]   = u1;       //1
      Lower[0]   = B1;       //1
      radmul(3,Lower,Upper,NIterMin,NIter,EPS,WK,NIter,Result,ResErr,NFNEVL,Summary);
      Integrala   = Integrala + Result;
      totsummary = totsummary + Summary;
      NEval      = NEval + NFNEVL;
    }
  }
  else{
    if (u2 < RM2 ){
      Integrala   = 0;
      totsummary = 0;
      NEval      = 0;
    }
    else if (B2 > RM2){
      //integral has 2 parts, b1 integral has only 1 component
      Integrala   = 0;
      totsummary = 20000;
      NEval      = 0;
      Lower[0]   = RM1;       //1
      Lower[1]   = RM2;       //2
      Lower[2]   = 0.;       //2
      Upper[2]   = 2.*starlightConstants::pi;    //3
      Upper[0]   = u1;       //1
      Upper[1]   = B2;       //2
      radmul(3,Lower,Upper,NIterMin,NIter,EPS,WK,NIter,Result,ResErr,NFNEVL,Summary);
      Integrala   = Integrala + Result;
      totsummary = totsummary + 100*Summary;
      NEval      = NEval + NFNEVL;
      Upper[1]   = u2;       //2
      Lower[1]   = B2;       //2
      radmul(3,Lower,Upper,NIterMin,NIter,EPS,WK,NIter,Result,ResErr,NFNEVL,Summary);
      Integrala   = Integrala + Result;
      totsummary = totsummary + Summary;
      NEval      = NEval + NFNEVL;
    }
    else{                 //integral has 1 part
      Integrala   = 0;
      totsummary = 10000;
      NEval      = 0;
      Lower[0]   = RM1;       //1
      Lower[1]   = RM2;       //2
      Lower[2]   = 0.;       //3
      Upper[2]   = 2.*starlightConstants::pi;    //3
      Upper[0]   = u1;       //1
      Upper[1]   = u2;       //2
      radmul(3,Lower,Upper,NIterMin,NIter,EPS,WK,NIter,Result,ResErr,NFNEVL,Summary);
      Integrala   = Integrala + Result;
      totsummary = totsummary + Summary;
      NEval      = NEval + NFNEVL;
    }
  }
  Integrala = 2*starlightConstants::pi*Integrala;
  delete [] WK;
  return Integrala;
}

//______________________________________________________________________________ 
double twoPhotonLuminosity::radmul(int N,double *A,double *B,int MINPTS,int MAXPTS,double EPS,double *WK,int IWK,double &RESULT,double &RELERR,double &NFNEVL,double &IFAIL)
{
  double wn1[14] = {       -0.193872885230909911, -0.555606360818980835,
			   -0.876695625666819078, -1.15714067977442459,  -1.39694152314179743,
			   -1.59609815576893754,  -1.75461057765584494,  -1.87247878880251983,
			   -1.94970278920896201,  -1.98628257887517146,  -1.98221815780114818,
			   -1.93750952598689219,  -1.85215668343240347,  -1.72615963013768225};
  
  double wn3[14] = {        0.0518213686937966768,    0.0314992633236803330,
			    0.0111771579535639891, -0.00914494741655235473,  -0.0294670527866686986,
			    -0.0497891581567850424, -0.0701112635269013768,   -0.0904333688970177241,
			    -0.110755474267134071,  -0.131077579637250419,    -0.151399685007366752,
			    -0.171721790377483099,  -0.192043895747599447,    -0.212366001117715794};
  
    double wn5[14] = {        0.871183254585174982e-01,  0.435591627292587508e-01,
			      0.217795813646293754e-01, 0.108897906823146873e-01,  0.544489534115734364e-02,
			      0.272244767057867193e-02, 0.136122383528933596e-02,  0.680611917644667955e-03,
			      0.340305958822333977e-03, 0.170152979411166995e-03,  0.850764897055834977e-04,
			      0.425382448527917472e-04, 0.212691224263958736e-04,  0.106345612131979372e-04};

    double wpn1[14] = {       -1.33196159122085045, -2.29218106995884763,
			      -3.11522633744855959, -3.80109739368998611, -4.34979423868312742,
			      -4.76131687242798352, -5.03566529492455417, -5.17283950617283939,
			      -5.17283950617283939, -5.03566529492455417, -4.76131687242798352,
			      -4.34979423868312742, -3.80109739368998611, -3.11522633744855959};
    
    double wpn3[14] = {        0.0445816186556927292, -0.0240054869684499309,
			       -0.0925925925925925875, -0.161179698216735251,  -0.229766803840877915,
			       -0.298353909465020564,  -0.366941015089163228,  -0.435528120713305891,
			       -0.504115226337448555,  -0.572702331961591218,  -0.641289437585733882,
			       -0.709876543209876532,  -0.778463648834019195,  -0.847050754458161859};

    RESULT = 0;

    double ABSERR = 0.;
    double ctr[15], wth[15], wthl[15], z[15];
    double R1  = 1.;
    double HF  = R1/2.;
    double xl2 = 0.358568582800318073;
    double xl4 = 0.948683298050513796;
    double xl5 = 0.688247201611685289;
    double w2  = 980./6561.;
    double w4  = 200./19683.;
    double wp2 = 245./486.;
    double wp4 = 25./729.;
    int j1 =0;
    double SUM1, SUM2, SUM3, SUM4, SUM5, RGNCMP=0., RGNVAL, RGNERR, F2, F3, DIF, DIFMAX, IDVAXN =0.;
    IFAIL  = 3;
    if (N < 2 || N > 15) return 0;
    if (MINPTS > MAXPTS) return 0;
    int IFNCLS = 0;
    int IDVAX0 = 0;
    bool LDV = false;
    double TWONDM = pow(2.,(double)(N));
    int IRGNST = 2*N+3;
    int IRLCLS = (int)pow(2.,(N))+2*N*(N+1)+1;
    int ISBRGN = IRGNST;
    int ISBTMP, ISBTPP;
    int ISBRGS = IRGNST;

    if ( MAXPTS < IRLCLS ) return 0;
    for ( int j = 0; j < N; j++ ){  //10
      ctr[j]=(B[j] + A[j])*HF;
      wth[j]=(B[j] - A[j])*HF;      
    }
 L20:
    double RGNVOL = TWONDM; //20
    for ( int j = 0; j < N; j++ ){            //30
      RGNVOL = RGNVOL*wth[j];
      z[j] = ctr[j];  
    }
    
    SUM1 = integrand(N,z); 
    
    DIFMAX = 0;
    SUM2   = 0;
    SUM3   = 0;
    for ( int j = 0; j < N; j++ ) {   //40
      z[j]=ctr[j]-xl2*wth[j];
      F2=integrand(N,z);
      
      z[j]=ctr[j]+xl2*wth[j];
      F2=F2+integrand(N,z);
      wthl[j]=xl4*wth[j];
      
      z[j]=ctr[j]-wthl[j];
      F3=integrand(N,z);
      
      z[j]=ctr[j]+wthl[j];
      F3=F3+integrand(N,z);
      SUM2=SUM2+F2;
      SUM3=SUM3+F3;
      DIF=fabs(7.*F2-F3-12.*SUM1);
      DIFMAX=max(DIF,DIFMAX);
      
      if ( DIFMAX == DIF) IDVAXN = j+1;
      z[j]=ctr[j];
      
    }
    
    SUM4   = 0;
    for ( int j = 1; j < N; j++){ //70
      
	j1=j-1;
        for ( int k = j; k < N; k++){  //60
	  for ( int l = 0; l < 2; l++){      //50
	    wthl[j1]=-wthl[j1];
	    z[j1]=ctr[j1]+wthl[j1];
	    for ( int m = 0; m < 2; m++){   //50
	      wthl[k]=-wthl[k];
	      z[k]=ctr[k]+wthl[k];
	      SUM4=SUM4+integrand(N,z);
	    }
	  }
	  z[k]=ctr[k];
	}
        z[j1]=ctr[j1];
    }
    
    SUM5  = 0;
    
    for ( int j = 0; j < N; j++){             //80
      wthl[j]=-xl5*wth[j];
      z[j]=ctr[j]+wthl[j];
    }
 L90:
    SUM5=SUM5+integrand(N,z);   //line 90
    
    for (int j = 0; j < N; j++){ //100
      wthl[j]=-wthl[j];
      z[j]=ctr[j]+wthl[j];  
      if ( wthl[j] > 0. ) goto L90;
    }

    RGNCMP = RGNVOL*(wpn1[N-2]*SUM1+wp2*SUM2+wpn3[N-2]*SUM3+wp4*SUM4);
    RGNVAL = wn1[N-2]*SUM1+w2*SUM2+wn3[N-2]*SUM3+w4*SUM4+wn5[N-2]*SUM5;
   
    RGNVAL = RGNVOL*RGNVAL;
    RGNERR = fabs(RGNVAL-RGNCMP);
    RESULT = RESULT+RGNVAL;
    ABSERR = ABSERR+RGNERR;
    IFNCLS = IFNCLS+IRLCLS;
    
    
    if (LDV){
    L110:

      ISBTMP = 2*ISBRGN;

      if ( ISBTMP > ISBRGS ) goto L160;
      if ( ISBTMP < ISBRGS ){
	ISBTPP = ISBTMP + IRGNST;
	
	if ( WK[ISBTMP-1] < WK[ISBTPP-1] ) ISBTMP = ISBTPP;
      }
      
      if ( RGNERR >= WK[ISBTMP-1] ) goto L160;
      for ( int k = 0; k < IRGNST; k++){
	WK[ISBRGN-k-1] = WK[ISBTMP-k-1];
                }
      ISBRGN = ISBTMP;
      goto L110;
    }
 L140:
    
    ISBTMP = (ISBRGN/(2*IRGNST))*IRGNST;
    
    if ( ISBTMP >= IRGNST && RGNERR > WK[ISBTMP-1] ){
      for ( int k = 0; k < IRGNST; k++){
	WK[ISBRGN-k-1]=WK[ISBTMP-k-1];
      }
      ISBRGN = ISBTMP;
      goto L140;
    }
 L160:
    
    WK[ISBRGN-1] = RGNERR;
    WK[ISBRGN-2] = RGNVAL;
    WK[ISBRGN-3] = IDVAXN;
    
    for ( int j = 0; j < N; j++) {
      
      ISBTMP = ISBRGN-2*j-4;
      WK[ISBTMP]=ctr[j];
      WK[ISBTMP-1]=wth[j];
    }
    if (LDV) {
      LDV = false;
      ctr[IDVAX0-1]=ctr[IDVAX0-1]+2*wth[IDVAX0-1];
      ISBRGS = ISBRGS + IRGNST;
      ISBRGN = ISBRGS;
      goto L20;
    }
    
    RELERR=ABSERR/fabs(RESULT);
    
    
    if ( ISBRGS + IRGNST > IWK ) IFAIL = 2;
    if ( IFNCLS + 2*IRLCLS > MAXPTS ) IFAIL = 1;
    if ( RELERR < EPS && IFNCLS >= MINPTS ) IFAIL = 0;
    
    if ( IFAIL == 3 ) {
      LDV = true;
      ISBRGN = IRGNST;
      ABSERR = ABSERR-WK[ISBRGN-1];
      RESULT = RESULT-WK[ISBRGN-2];
      IDVAX0 = (int)WK[ISBRGN-3];
      
      for ( int j = 0; j < N; j++) {
	ISBTMP = ISBRGN-2*j-4;
	ctr[j] = WK[ISBTMP];
	wth[j] = WK[ISBTMP-1];
      }
      
      wth[IDVAX0-1] = HF*wth[IDVAX0-1];
      ctr[IDVAX0-1] = ctr[IDVAX0-1]-wth[IDVAX0-1];
      goto L20;
    }
    NFNEVL=IFNCLS;
    return 1;
}


//______________________________________________________________________________
double twoPhotonLuminosity::integrand(double ,  // N (unused)
                                      double X[])
{
  double  b1 = X[0];      //1
  double  b2 = X[1];      //2
  double  theta = X[2];   //3
  //breakup effects distances in fermis, so convert to fermis(factor of hbarcmev)
  double  D = sqrt(b1*b1+b2*b2-2*b1*b2*cos(theta))*starlightConstants::hbarcmev;
  double  integrandx = Nphoton(_W1,_gamma,b1)*Nphoton(_W2,_gamma,b2)*b1*b2*probabilityOfBreakup(D); 
  return integrandx;
}


//______________________________________________________________________________
double twoPhotonLuminosity::Nphoton(double W,double gamma,double Rho)
{
 double Nphoton1 =0.;
 double WGamma = W/gamma;
 double WGR    = 1.0*WGamma*Rho;
 //factor of Z^2*alpha is omitted
 double Wphib = WGamma*bessel::dbesk1(WGR);

 Nphoton1 = 1.0/(starlightConstants::pi*starlightConstants::pi)*(Wphib*Wphib);
 return Nphoton1;
}
