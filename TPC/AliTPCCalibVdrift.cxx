/**************************************************************************
 * Copyright(c) 2006-07, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/


///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// Class describing the Vdrift dependencies on E,T,P and GasComposition      //
// Authors: Stefan Rossegger, Haavard Helstrup                               //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "TSystem.h"
#include "TObject.h"
#include "TMath.h"
#include "AliTPCTempMap.h"
#include "AliTPCSensorTempArray.h"

#include "AliTPCCalibVdrift.h"

ClassImp(AliTPCCalibVdrift)

namespace paramDefinitions {
    
    // Standard Conditions used as origin in the Magbolz simulations
    // Dimesions E [kV/cm], T [K], P [TORR], Cco2 [%], Cn2 [%]
    const Double_t kstdE = 400;
    const Double_t kstdT = 293;
    const Double_t kstdP = 744;
    const Double_t kstdCco2 = 9.52;
    const Double_t kstdCn2 = 4.76;
    // Driftvelocity at Standardcontitions [cm/microSec]
    const Double_t kstdVdrift = 2.57563;
    
    // Vdrift dependencies simulated with Magbolz [%(Vdrift)/[unit]]
    const Double_t kdvdE = 0.24;
    const Double_t kdvdT = 0.30;
    const Double_t kdvdP = -0.13;
    const Double_t kdvdCco2 = -6.60;
    const Double_t kdvdCn2 = -1.74;
  
    Double_t krho = 0.934246; // density of TPC-Gas [kg/m^3]
                              // method of calculation: weighted average
    Double_t kg = 9.81;
}

using namespace paramDefinitions;

AliTPCCalibVdrift::AliTPCCalibVdrift(AliTPCSensorTempArray *SensTemp, AliDCSSensor *SensPres, TObject *SensGasComp):
  TNamed(),
  fSensTemp(0),
  fSensPres(0),
  fTempMap(0),
  fSensGasComp(0)
{
  //
  //  Standard constructor
  //

  fSensTemp = SensTemp;
  fSensPres = SensPres;
  fTempMap  = new AliTPCTempMap(fSensTemp);
  fSensGasComp = SensGasComp;
}

AliTPCCalibVdrift::AliTPCCalibVdrift(const AliTPCCalibVdrift& source) :
  TNamed(source),
  fSensTemp(source.fSensTemp),
  fSensPres(source.fSensPres),
  fTempMap(source.fTempMap),
  fSensGasComp(source.fSensGasComp)
{
  //
  //  Copy constructor
  //
}

//_____________________________________________________________________________

AliTPCCalibVdrift& AliTPCCalibVdrift::operator=(const AliTPCCalibVdrift& source){
  //
  // assignment operator
  //
  if (&source == this) return *this;
  new (this) AliTPCCalibVdrift(source);
  
  return *this;  
}

//_____________________________________________________________________________
AliTPCCalibVdrift::~AliTPCCalibVdrift()
{
  //
  // AliTPCCalibVdrift destructor
  // 

}

Double_t AliTPCCalibVdrift::GetPTRelative(UInt_t timeSec, Int_t side){
  //
  // Get Relative difference of p/T for given time stamp
  // timeSec - absolute time
  // side    - 0 - A side 1 -C side
  //
  TTimeStamp tstamp(timeSec);
  if (!fSensPres) return 0;
  Double_t pressure    = fSensPres->GetValue(tstamp);
  TLinearFitter * fitter = fTempMap->GetLinearFitter(3,side,tstamp);
  if (!fitter) return 0;
  TVectorD vec;
  fitter->GetParameters(vec);
  delete fitter;
  if (vec[0]<10) return 0;
  Double_t  temperature = vec[0]+273.15;
  Double_t povertMeas = pressure/temperature;
  const Double_t torrTokPascal = 0.75006;
  Double_t povertNom =  kstdP/(torrTokPascal*kstdT);
  return povertMeas/povertNom;
}


//_____________________________________________________________________________
Double_t AliTPCCalibVdrift::VdriftLinearHyperplaneApprox(Double_t dE, Double_t dT, Double_t dP, Double_t dCco2, Double_t dCn2) 
{
  //
  // Returns approximated value for the driftvelocity based on  
  // linear Hyperplane approximation (~ Taylorapproximation of 1st order)
  //

  Double_t vdrift = (dE*kdvdE+dT*kdvdT+dP*kdvdP+dCco2*kdvdCco2+dCn2*kdvdCn2);
  
  return vdrift; 

}
//_____________________________________________________________________________

Double_t AliTPCCalibVdrift::GetVdriftNominal() 
{
  // returns nominal Driftvelocity at StandardConditions
  return kstdVdrift;
}

//_____________________________________________________________________________

Double_t AliTPCCalibVdrift::GetVdriftChange(Double_t x, Double_t y, Double_t z, UInt_t timeSec)
{
  // 
  // Calculates Vdrift change in percent of Vdrift_nominal 
  // (under nominal conditions) at x,y,z,timeSec
  //

  // Get E-field Value --------------------------
  Double_t dE = 0; //FIXME: eventually include Field-Inhomogenities

  // Get Temperature Value ----------------------  
  AliTPCTempMap *tempMap = fTempMap;
  Double_t tempValue = tempMap->GetTemperature(x, y, z, timeSec);
  Double_t dT = tempValue+273.15 - kstdT;
  
  // Get Main Pressure Value ---------------------
  // FIXME: READ REAL PRESSURE SENSOR  
  //        through TObject *fSensPres; 
  //        e.g. Double_t PO = fSensPres->GetValue(timeSec);  
  Double_t p0 = 744;
  // recalculate Pressure according to height in TPC and transform to
  // TORR (with simplified hydrostatic formula)   
  Double_t dP = p0 - krho*kg*y/10000 /1000*760 - kstdP;
   
  // Get GasComposition
  // FIXME: include Goofy values for CO2 and N2 conzentration out of DCS? 
  //   through TObject *fSensGasComp and calculate difference to stdCondit.
  Double_t dCco2 = 0;
  Double_t dCn2 = 0;

  // Calculate change in drift velocity in terms of Vdrift_nominal
  Double_t vdrift = VdriftLinearHyperplaneApprox(dE, dT, dP, dCco2, dCn2); 
  
  return vdrift;    
}

//_____________________________________________________________________________

Double_t AliTPCCalibVdrift::GetMeanZVdriftChange(Double_t x, Double_t y, UInt_t timeSec)
{
  // 
  // Calculates Meanvalue in z direction of Vdrift change in percent 
  // of Vdrift_nominal (under standard conditions) at position x,y,timeSec
  // with help of 'nPopints' base points
  //
  
  Int_t nPoints = 5;
 
  Double_t vdriftSum = 0;

  for (Int_t i = 0; i<nPoints; i++) {
    Double_t z = (Double_t)i/(nPoints-1)*500-250;
    vdriftSum = vdriftSum + GetVdriftChange(x, y, z, timeSec);
  }
  
  Double_t meanZVdrift = vdriftSum/nPoints;

  return meanZVdrift;

}

//_____________________________________________________________________________

TGraph *AliTPCCalibVdrift::MakeGraphMeanZVdriftChange(Double_t x, Double_t y, Int_t nPoints)
{
  //
  // Make graph from start time to end time of Mean Drift Velocity in 
  // Z direction at given x and y position
  //

  UInt_t startTime = fSensTemp->GetStartTime();
  UInt_t endTime = fSensTemp->GetEndTime();
  
  UInt_t stepTime = (endTime - startTime)/nPoints;


  Double_t *xvec = new Double_t[nPoints];
  Double_t *yvec = new Double_t[nPoints];

  for (Int_t ip=0; ip<nPoints; ip++) {
    xvec[ip] = startTime+ip*stepTime;
    yvec[ip] = GetMeanZVdriftChange(x, y, ip*stepTime);
  }

  TGraph *graph = new TGraph(nPoints,xvec,yvec);

  delete [] xvec;
  delete [] yvec;

  graph->GetXaxis()->SetTimeDisplay(1);
  graph->GetXaxis()->SetLabelOffset(0.02);
  graph->GetXaxis()->SetTimeFormat("#splitline{%d/%m}{%H:%M}");

  return graph;
}
