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
// $Rev:: 45                          $: revision of last commit
// $Author:: bgrube                   $: author of last commit
// $Date:: 2011-02-27 20:52:35 +0100 #$: date of last commit
//
// Description:
//
//
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
#include "incoherentPhotonNucleusLuminosity.h"


using namespace std;
using namespace starlightConstants;


//______________________________________________________________________________
incoherentPhotonNucleusLuminosity::incoherentPhotonNucleusLuminosity(const inputParameters& inputParametersInstance, beamBeamSystem& bbsystem)
  : photonNucleusCrossSection(inputParametersInstance, bbsystem)
  ,_baseFileName(inputParametersInstance.baseFileName())
  ,_beamLorentzGamma(inputParametersInstance.beamLorentzGamma())
  ,_maxW(inputParametersInstance.maxW())
  ,_minW(inputParametersInstance.minW())
  ,_nmbWBins(inputParametersInstance.nmbWBins())
  ,_maxRapidity(inputParametersInstance.maxRapidity())
  ,_nmbRapidityBins(inputParametersInstance.nmbRapidityBins())
  ,_productionMode(inputParametersInstance.productionMode())
  ,_beamBreakupMode(inputParametersInstance.beamBreakupMode())
  ,_interferenceEnabled(inputParametersInstance.interferenceEnabled())
  ,_interferenceStrength(inputParametersInstance.interferenceStrength())
  ,_maxPtInterference(inputParametersInstance.maxPtInterference())
  ,_nmbPtBinsInterference(inputParametersInstance.nmbPtBinsInterference())
  ,_protonEnergy(inputParametersInstance.protonEnergy())
  ,_parameterValueKey(inputParametersInstance.parameterValueKey())
{
  cout <<"Creating Luminosity Tables for incoherent vector meson production."<<endl;
  incoherentPhotonNucleusDifferentialLuminosity();
  cout <<"Luminosity Tables created."<<endl;
}


//______________________________________________________________________________
incoherentPhotonNucleusLuminosity::~incoherentPhotonNucleusLuminosity()
{ }


//______________________________________________________________________________
void incoherentPhotonNucleusLuminosity::incoherentPhotonNucleusDifferentialLuminosity()
{
  double W,dW,dY;
  double Egamma,Y;
  double testint,dndWdY;
  double C;  
  int beam; 

  std::string wyFileName;
  wyFileName = _baseFileName +".txt";

  ofstream wylumfile;
  wylumfile.precision(15);
  
  double  bwnorm,Eth;

  dW = (_wMax - _wMin)/_nWbins;
  dY  = (_yMax-(-1.0)*_yMax)/_nYbins;
    
  // Write the values of W used in the calculation to slight.txt.  
  wylumfile.open(wyFileName.c_str());
  wylumfile << getbbs().beam1().Z() <<endl;
  wylumfile << getbbs().beam1().A() <<endl;
  wylumfile << getbbs().beam2().Z() <<endl;
  wylumfile << getbbs().beam2().A() <<endl;
  wylumfile << _beamLorentzGamma <<endl;
  wylumfile << _maxW <<endl;
  wylumfile << _minW <<endl;
  wylumfile << _nmbWBins <<endl;
  wylumfile << _maxRapidity <<endl;
  wylumfile << _nmbRapidityBins <<endl;
  wylumfile << _productionMode <<endl;
  wylumfile << _beamBreakupMode <<endl;
  wylumfile << _interferenceEnabled <<endl;
  wylumfile << _interferenceStrength <<endl;
  wylumfile << starlightConstants::deuteronSlopePar <<endl;
  wylumfile << _maxPtInterference <<endl;
  wylumfile << _nmbPtBinsInterference <<endl;
  
  //     Normalize the Breit-Wigner Distribution and write values of W to slight.txt
  testint=0.0;
  //Grabbing default value for C in the breit-wigner calculation
  C=getDefaultC();
  for(unsigned int i = 0; i <= _nWbins - 1; ++i) {
    W = _wMin + double(i)*dW + 0.5*dW;
    testint = testint + breitWigner(W,C)*dW;
    wylumfile << W << endl;
  }
  bwnorm = 1./testint;
  
  //     Write the values of Y used in the calculation to slight.txt.
  for(unsigned int i = 0; i <= _nYbins - 1; ++i) {
    Y = -1.0*_yMax + double(i)*dY + 0.5*dY;
    wylumfile << Y << endl;
  }

  int A_1 = getbbs().beam1().A(); 
  int A_2 = getbbs().beam2().A();

  // Do this first for the case when the first beam is the photon emitter 
  // Treat pA separately with defined beams 
  // The variable beam (=1,2) defines which nucleus is the target 
  for(unsigned int i = 0; i <= _nWbins - 1; ++i) {

    W = _wMin + double(i)*dW + 0.5*dW;

    double Ep = _protonEnergy;

    Eth=0.5*(((W+starlightConstants::protonMass)*(W+starlightConstants::protonMass)-starlightConstants::protonMass*starlightConstants::protonMass)/
	   (Ep + sqrt(Ep*Ep-starlightConstants::protonMass*starlightConstants::protonMass)));
    
    for(unsigned int j = 0; j <= _nYbins - 1; ++j) {

      Y = -1.0*_yMax + double(j)*dY + 0.5*dY;

      if( A_2 == 1 && A_1 != 1 ){
        // pA, first beam is the nucleus and photon emitter
        Egamma = 0.5*W*exp(Y);
        beam = 2; 
      } else if( A_1 ==1 && A_2 != 1){
        // pA, second beam is the nucleus and photon emitter
        Egamma = 0.5*W*exp(-Y); 
        beam = 1; 
      } else {
        Egamma = 0.5*W*exp(Y);        
        beam = 2; 
      }
      
      dndWdY = 0.; 

      if(Egamma > Eth){
	if(Egamma > maxPhotonEnergy())Egamma = maxPhotonEnergy();
        double Wgp = sqrt(2.*Egamma*(Ep+sqrt(Ep*Ep-starlightConstants::protonMass*
                                 starlightConstants::protonMass))+starlightConstants::protonMass*starlightConstants::protonMass);

        double localsig = sigmagp(Wgp); 
        if( A_1 == 1 && A_2 != 1 ){
          dndWdY = Egamma*photonFlux(Egamma,beam)*localsig*breitWigner(W,bwnorm); 
        }else if (A_2 ==1 && A_1 !=1){
          dndWdY = Egamma*photonFlux(Egamma,beam)*localsig*breitWigner(W,bwnorm); 
        }else{ 
          double csVN = sigma_N(Wgp);         
          double csVA = sigma_A(csVN,beam); 
          double csgA= (csVA/csVN)*sigmagp(Wgp); 
          dndWdY = Egamma*photonFlux(Egamma,beam)*csgA*breitWigner(W,bwnorm); 
        }
      }

      wylumfile << dndWdY << endl;

    }
  }

  // Repeat the loop for the case when the second beam is the photon emitter. 
  // Don't repeat for pA
  if( !( (A_2 == 1 && A_1 != 1) || (A_1 == 1 && A_2 != 1) ) ){ 
    for(unsigned int i = 0; i <= _nWbins - 1; ++i) {

      W = _wMin + double(i)*dW + 0.5*dW;

      double Ep = _protonEnergy;

      Eth=0.5*(((W+starlightConstants::protonMass)*(W+starlightConstants::protonMass)-starlightConstants::protonMass*starlightConstants::protonMass)/
	   (Ep + sqrt(Ep*Ep-starlightConstants::protonMass*starlightConstants::protonMass)));
    
      for(unsigned int j = 0; j <= _nYbins - 1; ++j) {

        Y = -1.0*_yMax + double(j)*dY + 0.5*dY;

        beam = 1; 
        Egamma = 0.5*W*exp(-Y);        
      
        dndWdY = 0.; 

        if(Egamma > Eth){
	  if(Egamma > maxPhotonEnergy())Egamma = maxPhotonEnergy();
          double Wgp = sqrt(2.*Egamma*(Ep+sqrt(Ep*Ep-starlightConstants::protonMass*
                                 starlightConstants::protonMass))+starlightConstants::protonMass*starlightConstants::protonMass);

          double csVN = sigma_N(Wgp);         
          double csVA = sigma_A(csVN,beam); 
          double csgA= (csVA/csVN)*sigmagp(Wgp); 
          dndWdY = Egamma*photonFlux(Egamma,beam)*csgA*breitWigner(W,bwnorm); 
        
        }

        wylumfile << dndWdY << endl;

      }
    }
  }

  wylumfile << bwnorm << endl;
  wylumfile << _parameterValueKey << endl;
  wylumfile.close();
  
}


