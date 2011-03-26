/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
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
//  TPC calibration class for temperature maps and tendencies                //
//  (based on TPC Temperature Sensors and FiniteElement Simulation)          //
//                                                                           //
//  Authors: Stefan Rossegger, Haavard Helstrup                              //
//                                                                           //
//  Note: Obvioulsy some changes by Marian, but when ???                     //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "AliTPCSensorTempArray.h"
#include "TLinearFitter.h"
#include "TString.h"
#include "TGraph2D.h"
#include "TTimeStamp.h"

#include "AliTPCTempMap.h"


ClassImp(AliTPCTempMap)
  
  const char kStringFEsimulation[] = "FEsimulation.txt";

//_____________________________________________________________________________
AliTPCTempMap::AliTPCTempMap(AliTPCSensorTempArray *sensorDCS):
  TNamed(),
  fTempArray(0),
  fStringFEsimulation(kStringFEsimulation)
{
  //
  // AliTPCTempMap default constructor
  //

  fTempArray = sensorDCS;

}

//_____________________________________________________________________________
AliTPCTempMap::AliTPCTempMap(const AliTPCTempMap &c):
  TNamed(c),
  fTempArray(c.fTempArray),
  fStringFEsimulation(c.fStringFEsimulation)
{
  //
  // AliTPCTempMap copy constructor
  //

}

//_____________________________________________________________________________
AliTPCTempMap::~AliTPCTempMap()
{
  //
  // AliTPCTempMap destructor
  //
  
}

//_____________________________________________________________________________
AliTPCTempMap &AliTPCTempMap::operator=(const AliTPCTempMap &c)
{
  //
  // Assignment operator
  //

  if (this != &c) ((AliTPCTempMap &) c).Copy(*this);
  return *this;

}

//_____________________________________________________________________________
void AliTPCTempMap::Copy(TObject &c) const
{
  //
  // Copy function
  //
  
  TObject::Copy(c);
  
}

//_____________________________________________________________________________

Double_t AliTPCTempMap::GetTempGradientY(UInt_t timeSec, Int_t side){
 //
 // Extract Linear Vertical Temperature Gradient [K/cm] within the TPC on 
 // Shaft Side(A): 0
 // Muon  Side(C): 1
 // Values based on TemperatureSensors within the TPC ( type: 3 (TPC) )
 //
 // FIXME: Also return residual-distribution, covariance Matrix
 //        or simply chi2 for validity check? 
 //        -> better use GetLinearFitter - function in this case!
  
 TLinearFitter *fitter = new TLinearFitter(3,"x0++x1++x2");
 TVectorD param(3);
 Int_t i = 0;

 Int_t nsensors = fTempArray->NumSensors();
 for (Int_t isensor=0; isensor<nsensors; isensor++) { // loop over all sensors
   AliTPCSensorTemp *entry = (AliTPCSensorTemp*)fTempArray->GetSensorNum(isensor);
   
   if (entry->GetType()==3 && entry->GetSide()==side) { // take SensorType:TPC 
     Double_t x[3];
     x[0]=1;
     x[1]=entry->GetX();
     x[2]=entry->GetY();    
     Double_t y = fTempArray->GetValue(timeSec,isensor); // get temperature value
     if (IsOK(y)) fitter->AddPoint(x,y,1); // add values to LinearFitter
     i++;
   }

 }  
 fitter->Eval();
 fitter->GetParameters(param);

 fitter->~TLinearFitter();

 return param[2]; // return vertical (Y) tempGradient in [K/cm]
  
}

//_____________________________________________________________________________
TLinearFitter *AliTPCTempMap::GetLinearFitter(Int_t type, Int_t side, TTimeStamp &stamp)
{
  //
  // absolute time stamp used
  // see AliTPCTempMap::GetLinearFitter(Int_t type, Int_t side, UInt_t timeSec) for details
  //
  Int_t timeSec = stamp.GetSec()-fTempArray->GetStartTime().GetSec();
  return GetLinearFitter(type,side,timeSec);
}

//_____________________________________________________________________________
TLinearFitter *AliTPCTempMap::GetLinearFitter(Int_t type, Int_t side, UInt_t timeSec)
{
  // 
  // Creates a TlinearFitter object for the desired region of the TPC 
  // (via choosen type and side of TPC temperature sensors) at a given 
  // timeSec (in secounds) after start time
  // type: 0 ... ReadOutChambers (ROC)
  //       1 ... OuterContainmentVessel (OFC)
  //       2 ... InnerContainmentVessel (IFC) + ThermalScreener (TS)
  //       3 ... Within the TPC (DriftVolume) (TPC)
  //       4 ... Only InnerContainmentVessel (IFC) 
  // side: Can be choosen for type 0 and 3 (otherwise it will be ignored in 
  //       in order to get all temperature sensors of interest)
  //       0 ... Shaft Side (A)
  //       1 ... Muon Side (C)
  // 

  TLinearFitter *fitter = new TLinearFitter(3);
  Double_t x[3]={0};
  Double_t y = 0;
  const Float_t kMaxDelta=0.5;
  
  if (type == 1 || type == 2 || type == 4) {
    fitter->SetFormula("x0++x1++TMath::Sin(x2)"); // returns Z,Y gradient
  } else {
    fitter->SetFormula("x0++x1++x2"); // returns X,Y gradient
  }

  Int_t i = 0;
  Int_t nsensors = fTempArray->NumSensors();

  Float_t temps[1000];
  for (Int_t isensor=0; isensor<nsensors; isensor++) { // loop over all sensors
    AliTPCSensorTemp *entry = (AliTPCSensorTemp*)fTempArray->GetSensorNum(isensor);
    if (entry->GetType()==type && entry->GetSide()==side){
      Float_t temperature= fTempArray->GetValue(timeSec,isensor); // get temperature value
      if (IsOK(temperature)) {temps[i]=temperature; i++;}
    }
  }
  Float_t medianTemp = TMath::Median(i, temps);
  if (i<3) return 0;
  Float_t rmsTemp = TMath::RMS(i, temps);
  
  i=0;
  
  for (Int_t isensor=0; isensor<nsensors; isensor++) { // loop over all sensors
    AliTPCSensorTemp *entry = (AliTPCSensorTemp*)fTempArray->GetSensorNum(isensor);
    
    if (type==0 || type==3) { // 'side' information used
      if (entry->GetType()==type && entry->GetSide()==side) {
	x[0]=1;
	x[1]=entry->GetX();
	x[2]=entry->GetY();    
	y = fTempArray->GetValue(timeSec,isensor); // get temperature value
	if (TMath::Abs(y-medianTemp)>kMaxDelta+4.*rmsTemp) continue;
	if (IsOK(y)) fitter->AddPoint(x,y,1); // add values to LinearFitter
	i++;
      }
    } else if (type==2) { // in case of IFC also usage of TS values
      if ((entry->GetType()==2) || (entry->GetType()==5)) {
	x[0]=1;
	x[1]=entry->GetZ();
	x[2]=entry->GetPhi();    
	y = fTempArray->GetValue(timeSec,isensor);
	if (TMath::Abs(y-medianTemp)>kMaxDelta+4.*rmsTemp) continue;
	if (IsOK(y)) fitter->AddPoint(x,y,1); 
	i++;
      }
    } else if (type==1){
      if (entry->GetType()==type) {
	x[0]=1;
	x[1]=entry->GetZ();
	x[2]=entry->GetPhi();    
	y = fTempArray->GetValue(timeSec,isensor);
	if (TMath::Abs(y-medianTemp)>kMaxDelta+4.*rmsTemp) continue;
	if (IsOK(y)) fitter->AddPoint(x,y,1);
	i++;	
      }
    } else if (type==4) { // ONLY IFC
      if (entry->GetType()==2) {
	x[0]=1;
	x[1]=entry->GetZ();
	x[2]=entry->GetPhi();    
	y = fTempArray->GetValue(timeSec,isensor);
	if (TMath::Abs(y-medianTemp)>kMaxDelta+4.*rmsTemp) continue;
	if (IsOK(y)) fitter->AddPoint(x,y,1); 
	i++;
      }
    }
  }  
  fitter->Eval();
  //fitter->EvalRobust(0.9); // Evaluates fitter
  
  delete [] x;

  return fitter; 

  // returns TLinearFitter object where Chi2, Fitparameters and residuals can 
  // be extracted via usual memberfunctions
  // example: fitter->GetParameters(param)
  // In case of type IFC or OFC, the parameters are the gradients in 
  // Z and Y direction (see fitformula)
  // Caution: Parameters are [K/cm] except Y at IFC,OFC ([K/radius]) 
}

//_____________________________________________________________________________

TGraph2D *AliTPCTempMap::GetTempMapsViaSensors(Int_t type, Int_t side, UInt_t timeSec)
{
  // 
  // Creates a TGraph2D object for the desired region of the TPC 
  // (via choosen type and side of TPC temperature sensors) at a given 
  // timeSec (in secounds) after start time
  // type: 0 ... ReadOutChambers (ROC)
  //       1 ... OuterContainmentVessel (OFC)
  //       2 ... InnerContainmentVessel (IFC) + ThermalScreener (TS)
  //       3 ... Within the TPC (DriftVolume) (TPC)
  // side: Can be choosen for type 0 and 3 (otherwise it will be ignored in 
  //       in order to get all temperature sensors of interest)
  //       0 ... A - side
  //       1 ... C - side
  // 

  TGraph2D *graph2D = new TGraph2D();

  Int_t i = 0;
  
  Int_t nsensors = fTempArray->NumSensors();
  for (Int_t isensor=0; isensor<nsensors; isensor++) { // loop over all sensors
    AliTPCSensorTemp *entry = (AliTPCSensorTemp*)fTempArray->GetSensorNum(isensor);

    Double_t x, y, z, r, phi, tempValue;
    x = entry->GetX();
    y = entry->GetY();
    z = entry->GetZ();
    r = entry->GetR();
    phi = entry->GetPhi();
    tempValue = fTempArray->GetValue(timeSec,isensor);
    //    printf("%d type %d: x=%lf y=%lf temp=%lf\n",isensor,entry->GetType(),x,y, tempValue);
    if (type==0 || type==3) { // 'side' information used
      if (entry->GetType()==type && entry->GetSide()==side) {
	graph2D->SetPoint(i,x,y,tempValue);
	i++;
      }
    } else if (type==2) { // in case of IFC also usage of TS values
      if (entry->GetType()==2 || entry->GetType()==5) {
	graph2D->SetPoint(i,z,phi,tempValue);
	i++;
      }
    } else if (type==1){
      if (entry->GetType()==type) {
	graph2D->SetPoint(i,z,phi,tempValue);
	i++;
      }
    }
  }  
  
  if (type==0 || type==3) {
    graph2D->GetXaxis()->SetTitle("X[cm]");
    graph2D->GetYaxis()->SetTitle("Y[cm]");
    if (type==0 && side==0) {
      graph2D->SetTitle("ROC A side");
    } else if (type==0 && side==1) {
      graph2D->SetTitle("ROC C side");
    } else if (type==3 && side==0) {
      graph2D->SetTitle("TPC A side (Inside the TPC)");
    } else if (type==3 && side==1) {
      graph2D->SetTitle("TPC C side (Inside the TPC)");
    }
  } else if (type==1 || type==2) {
    graph2D->GetXaxis()->SetTitle("Z[cm]");
    graph2D->GetYaxis()->SetTitle("Phi[RAD]");
    if (type==1) {
      graph2D->SetTitle("Outer Containment Vessel");
    } else if (type==2) {
      graph2D->SetTitle("Inner Containment Vessel");
    }
  }

  if (!graph2D->GetN()) {
    printf("Returned TGraph2D is empty: check type and side values\n");
  }

  graph2D->GetXaxis()->SetLabelOffset(0.0);
  graph2D->GetYaxis()->SetLabelOffset(0.005);
  graph2D->GetZaxis()->SetLabelOffset(-0.04);
  

  return graph2D; // returns TGgraph2D object
  
}


//_____________________________________________________________________________

TGraph *AliTPCTempMap::MakeGraphGradient(Int_t axis, Int_t side, Int_t nPoints)
{  
  //
  // Make graph from start time to end time of TempGradient in axis direction
  // axis: 0 ... horizontal Temperature Gradient (X)
  //       1 ... vertical Temperature Gradient (Y)
  //       2 ... longitudenal Temperature Gradient (Z) (side is ignored) 
  //             z gradient value based on OFC temperature sensors
  //             Caution!: better z gradient values through difference between 
  //             param[0] A- and param[0] C-side !
  // side for X and Y gradient: 
  //       0 ... Shaft Side (A)
  //       1 ... Muon Side (C)
  //
  
  TVectorD param(3);
  TLinearFitter *fitter = new TLinearFitter(3);

  UInt_t fStartTime = fTempArray->AliTPCSensorTempArray::GetStartTime();
  UInt_t fEndTime = fTempArray->AliTPCSensorTempArray::GetEndTime();
  
  UInt_t stepTime = (fEndTime-fStartTime)/nPoints;

  Double_t *x = new Double_t[nPoints];
  Double_t *y = new Double_t[nPoints];
  for (Int_t ip=0; ip<nPoints; ip++) {
    x[ip] = fStartTime+ip*stepTime;
    if (axis==2) {// Gradient in Z direction (based on OFC tempSensors)
      fitter = GetLinearFitter(1, side, ip*stepTime);
    } else {// Gradient in X or Y direction (based on TPC tempSensors)
      fitter = GetLinearFitter(3, side, ip*stepTime);
    }
    fitter->GetParameters(param);
    // multiplied by 500 since TempGradient is in [K/cm] 
    // (TPC diameter and length ~500cm)
    if (axis==1) { // Y axis
      y[ip] = param[2]*500;
    } else { // X axis
      y[ip] = param[1]*500;
    }
  }

  TGraph *graph = new TGraph(nPoints,x,y);

  fitter->~TLinearFitter(); 
  delete [] x;
  delete [] y;

  graph->GetXaxis()->SetTimeDisplay(1);
  graph->GetXaxis()->SetLabelOffset(0.02);
  graph->GetXaxis()->SetTimeFormat("#splitline{%d/%m}{%H:%M}");

  return graph;
}


//_____________________________________________________________________________
Double_t AliTPCTempMap::GetTemperature(Double_t x, Double_t y, Double_t z, TTimeStamp &stamp)
{
  //
  // absolute time stamp used
  // see also Double_t AliTPCTempMap::GetTemperature(Double_t x, Double_t y, Double_t z, UInt_t timeSec) for details
  //

  Int_t timeSec = stamp.GetSec()-fTempArray->GetStartTime().GetSec();
  return GetTemperature(x, y, z, timeSec);
}

//_____________________________________________________________________________

Double_t AliTPCTempMap::GetTemperature(Double_t x, Double_t y, Double_t z, UInt_t timeSec)
{  
  //
  // Returns estimated Temperature at given position (x,y,z[cm]) at given time 
  // (timeSec) after starttime
  // Method: so far just a linear interpolation between Linar fits of 
  //         the TPC temperature sensors
  //         FIXME: 'Educated Fit' through FiniteElement Simulation results!
  // FIXXME: Return 0? if x,y,z out of range
  //
  
  TVectorD paramA(3), paramC(3);
  TLinearFitter *fitterA = 0;
  TLinearFitter *fitterC = 0;

  fitterA = GetLinearFitter(3, 0, timeSec);
  fitterA->GetParameters(paramA);
  fitterC = GetLinearFitter(3, 1, timeSec);
  fitterC->GetParameters(paramC);

  Double_t fvalA = paramA[0]+paramA[1]*x+paramA[2]*y;
  Double_t fvalC = paramC[0]+paramC[1]*x+paramC[2]*y;

  Double_t k = (fvalA-fvalC)/(2*247);
  Double_t tempValue = fvalC+(fvalA-fvalC)/2+k*z;

  delete fitterA;
  delete fitterC;

  return tempValue;

}


Bool_t  AliTPCTempMap::IsOK(Float_t value){
  //
  // checks if value is within a certain range
  //
  const Float_t kMinT=15;
  const Float_t kMaxT=25;
  return (value>kMinT && value<kMaxT);
}
