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
  // 2nd order effect Taylor expansion
  const Double_t kdvdE2nd = -0.00107628;
  const Double_t kdvdT2nd = -0.00134441;
  const Double_t kdvdP2nd = 0.000135325;
  const Double_t kdvdCco22nd = 0.328761;
  const Double_t kdvdCn22nd = 0.151605;

  const Double_t torrTokPascal = 0.750061683;
 
  Double_t krho = 0.934246; // density of TPC-Gas [kg/m^3]
                            // method of calculation: weighted average
  Double_t kg = 9.81;

  //
  // Nominal value obtained from 2008 data
  //
  const Double_t kKelvin       =273.15; // degree to Kelvin
  const Double_t kNominalTemp  =19.03;  // mean between A and C side  in degree
  const Double_t kNominalPress =973.9;  // pressure sensor - in mbar- 
                                        // calibDB->GetPressure(tstamp,irun,1)
}


using namespace paramDefinitions;

AliTPCCalibVdrift::AliTPCCalibVdrift():
  TNamed(),
  fSensTemp(0),
  fSensPres(0),
  fTempMap(0),
  fSensGasComp(0),
  fNominalTemp(0),    // nominal temperature in Kelvin
  fNominalPress(0)    // nominal pressure    in mbar 
{
  //
  //  default constructor
  //
}

AliTPCCalibVdrift::AliTPCCalibVdrift(AliTPCSensorTempArray *SensTemp, AliDCSSensor *SensPres, TObject *SensGasComp):
  TNamed(),
  fSensTemp(0),
  fSensPres(0),
  fTempMap(0),
  fSensGasComp(0),
  fNominalTemp(0),    // nominal temperature in Kelvin
  fNominalPress(0)    // nominal pressure    in mbar 
{
  //
  //  Standard constructor
  //

  fSensTemp = SensTemp;
  fSensPres = SensPres;
  if (fSensTemp) {
    fTempMap  = new AliTPCTempMap(fSensTemp);
  } else {
    fTempMap = 0;
  }
  fSensGasComp = SensGasComp;
  fNominalTemp = kNominalTemp;
  fNominalPress= kNominalPress;
}

//_____________________________________________________________________________
AliTPCCalibVdrift::AliTPCCalibVdrift(const AliTPCCalibVdrift& source) :
  TNamed(source),
  fSensTemp(source.fSensTemp),
  fSensPres(source.fSensPres),
  fTempMap(source.fTempMap),
  fSensGasComp(source.fSensGasComp),
  fNominalTemp(source.fNominalTemp),    // nominal temperature in Kelvin
  fNominalPress(source.fNominalPress)    // nominal pressure    in mbar 

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

//_____________________________________________________________________________
Double_t AliTPCCalibVdrift::GetPTRelative(UInt_t absTimeSec, Int_t side){
  //
  // Get Relative difference of p/T for given time stamp
  // absTimeSec - absolute time in secounds
  // side: 0 - A side |  1 - C side
  //

  TTimeStamp tstamp(absTimeSec);

  if (!fSensPres||!fSensTemp) return 0;
  Double_t pressure = fSensPres->GetValue(tstamp);
  TLinearFitter * fitter = fTempMap->GetLinearFitter(3,side,tstamp);
  if (!fitter) return 0;
  TVectorD vec;
  fitter->GetParameters(vec);
  delete fitter;
  if (vec[0]<10) return 0;
  //
  //
  //
  Double_t  temperature = vec[0];  //vec[0] temeperature 
  Double_t  tpnom       = (fNominalTemp+kKelvin)/(fNominalPress);
  Double_t  tpmeasured  = (temperature+kKelvin)/(pressure);
  Double_t  result      = (tpmeasured-tpnom)/tpnom;

  return result;

}


//_____________________________________________________________________________
Double_t AliTPCCalibVdrift::VdriftLinearHyperplaneApprox(Double_t dE, Double_t dT, Double_t dP, Double_t dCco2, Double_t dCn2) 
{
  //
  // Returns approximated value for the driftvelocity change (in percent)
  // based on a Hyperplane approximation (~ Taylorapproximation of 2nd order)
  //

  Double_t termE   = dE*kdvdE + TMath::Power(dE,2)*kdvdE2nd;
  Double_t termT   = dT*kdvdT + TMath::Power(dT,2)*kdvdT2nd;
  Double_t termP   = dP*kdvdP + TMath::Power(dP,2)*kdvdP2nd;
  Double_t termCo2 = dCco2*kdvdCco2 + TMath::Power(dCco2,2)*kdvdCco22nd;
  Double_t termN2  = dCn2*kdvdCn2 + TMath::Power(dCn2,2)*kdvdCn22nd;

  Double_t vdChange = termE+termT+termP+termCo2+termN2;

  return vdChange;

}

//_____________________________________________________________________________

Double_t AliTPCCalibVdrift::GetVdriftNominal() 
{
  // returns nominal Driftvelocity at StandardConditions
  return kstdVdrift;
}

//_____________________________________________________________________________

Double_t AliTPCCalibVdrift::GetVdriftChange(Double_t x, Double_t y, Double_t z, UInt_t absTimeSec)
{
  // 
  // Calculates Vdrift change in percent of Vdrift_nominal 
  // (under nominal conditions) at x,y,z at absolute time (in sec)
  //

  TTimeStamp tstamp(absTimeSec);

  // Get E-field Value --------------------------
  Double_t dE = 0.23; // StandardOffset if CE is set to 100kV

  // Get Temperature Value ----------------------  
  AliTPCTempMap *tempMap = fTempMap;
  Double_t dT = 0;
  if (fTempMap) {
    Double_t tempValue = tempMap->GetTemperature(x, y, z, tstamp);
    dT = tempValue + 273.15 - kstdT;
  }
    
  // Get Main Pressure Value ---------------------
  Double_t dP = 0;
  if (fSensPres==0) {
    // Just the pressure drop over the TPC height
    dP = - krho*kg*y/10000*torrTokPascal;
  } else {
    // pressure sensors plus additional 0.4mbar overpressure within the TPC
    Double_t pressure = fSensPres->GetValue(tstamp) + 0.4; 
    // calculate pressure drop according to height in TPC and transform to
    // TORR (with simplified hydrostatic formula)
    dP = (pressure - krho*kg*y/10000) * torrTokPascal - kstdP;
  }

  // Get GasComposition
  // FIXME: include Goofy values for CO2 and N2 conzentration out of OCDB
  //        Goofy not yet reliable ... 
  Double_t dCco2 = 0;
  Double_t dCn2 = 0;

  // Calculate change in drift velocity in terms of Vdrift_nominal
  Double_t vdChange = VdriftLinearHyperplaneApprox(dE, dT, dP, dCco2, dCn2); 
  
  return vdChange;
    
}

//_____________________________________________________________________________

Double_t AliTPCCalibVdrift::GetMeanZVdriftChange(Double_t x, Double_t y, UInt_t absTimeSec)
{
  // 
  // Calculates Meanvalue in z direction of Vdrift change in percent 
  // of Vdrift_nominal (under standard conditions) at position x,y,absTimeSec
  // with help of 'nPopints' base points
  //
  
  Int_t nPoints = 5;
 
  Double_t vdriftSum = 0;

  for (Int_t i = 0; i<nPoints; i++) {
    Double_t z = (Double_t)i/(nPoints-1)*500-250;
    vdriftSum = vdriftSum + GetVdriftChange(x, y, z, absTimeSec);
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
    yvec[ip] = GetMeanZVdriftChange(x, y, fSensTemp->GetStartTime().GetSec() + ip*stepTime);
  }

  TGraph *graph = new TGraph(nPoints,xvec,yvec);

  delete [] xvec;
  delete [] yvec;

  graph->GetXaxis()->SetTimeDisplay(1);
  graph->GetXaxis()->SetLabelOffset(0.02);
  graph->GetXaxis()->SetTimeFormat("#splitline{%d/%m}{%H:%M}");

  return graph;
}
