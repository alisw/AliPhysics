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


////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// Class describing TPC temperature sensors (including pointers to graphs/fits//
// Authors: Marian Ivanov, Haavard Helstrup and Martin Siska                  //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////


#include "AliDCSSensor.h"
ClassImp(AliDCSSensor)

const Double_t kSecInHour = 3600.; // seconds in one hour



AliDCSSensor::AliDCSSensor():
  fId(),
  fIdDCS(0),
  fStringID(),
  fStartTime(0),
  fEndTime(0),
  fGraph(0),
  fFit(0),
  fX(0),
  fY(0),
  fZ(0)
{
  //
  //  Standard constructor
  //
}

AliDCSSensor::AliDCSSensor(const AliDCSSensor& source) :
   TNamed(source),
   fId(source.fId),
   fIdDCS(source.fIdDCS),
   fStringID(source.fStringID),
   fStartTime(source.fStartTime),
   fEndTime(source.fEndTime),
   fGraph(source.fGraph),
   fFit(source.fFit),
   fX(source.fX),
   fY(source.fY),
   fZ(source.fZ)
//
//  Copy constructor
//
{ }

AliDCSSensor& AliDCSSensor::operator=(const AliDCSSensor& source){
//
// assignment operator
//
  if (&source == this) return *this;
  new (this) AliDCSSensor(source);

  return *this;
}

//_____________________________________________________________________________
Double_t AliDCSSensor::GetValue(UInt_t timeSec)
{
 //
 // Get temperature value for actual sensor
 //  timeSec given as offset from start-of-run measured in seconds
 //
 Bool_t inside;
 return Eval(TTimeStamp(fStartTime+timeSec),inside);
}
//_____________________________________________________________________________
Double_t AliDCSSensor::GetValue(TTimeStamp time) 
{
 // Get temperature value for actual sensor
 //  time given as absolute TTimeStamp
 //
 Bool_t inside;
 return Eval(time, inside);
}

//_____________________________________________________________________________

Double_t AliDCSSensor::Eval(const TTimeStamp& time, Bool_t inside) const
{
  // 
  // Return temperature at given time
  //  If time < start of map  return value at start of map, inside = false
  //  If time > end of map    return value at end of map, inside = false
  
  UInt_t timeSec = time.GetSec();
  UInt_t diff = timeSec-fStartTime;
  inside = true;
  
  if ( timeSec < fStartTime ) { 
     inside=false;
     diff=0;
  }
  if ( timeSec > fEndTime ) {
     inside=false;
     diff = fEndTime-fStartTime;
  }
 
  Double_t timeHour = diff/kSecInHour;
  if ( fFit ) {
     return fFit->Eval(timeHour); 
  } else {
     return -99;
  }
}

TGraph* AliDCSSensor::MakeGraph(Int_t nPoints) const
{
  //
  // Make graph from start time to end time of DCS values 
  //

  UInt_t stepTime = (fEndTime-fStartTime)/nPoints;
  
  if ( !fFit ) return 0;

  Double_t *x = new Double_t[nPoints+1];
  Double_t *y = new Double_t[nPoints+1];
  for (Int_t ip=0; ip<nPoints; ip++) {
    x[ip] = fStartTime+ip*stepTime;
    y[ip] = fFit->Eval(ip*stepTime/kSecInHour);
  }
  
  TGraph *graph = new TGraph(nPoints,x,y);
  delete [] x;
  delete [] y;
  
  graph->GetXaxis()->SetTimeDisplay(1);
  graph->GetXaxis()->SetLabelOffset(0.02);
  graph->GetXaxis()->SetTimeFormat("#splitline{%d/%m}{%H:%M}");

  return graph;
}

//_____________________________________________________________________________

TClonesArray * AliDCSSensor::ReadTree(TTree* tree) {
  //
  // read values from ascii file
  //

  Int_t nentries = tree->GetEntries();

  char stringId[100];
  Int_t num=0;
  Int_t idDCS=0;
  Double_t x=0;
  Double_t y=0;
  Double_t z=0;

  tree->SetBranchAddress("StringID",&stringId);
  tree->SetBranchAddress("IdDCS",&idDCS);
  tree->SetBranchAddress("Num",&num);
  tree->SetBranchAddress("X",&x);
  tree->SetBranchAddress("Y",&y);
  tree->SetBranchAddress("Z",&z);

  // firstSensor = (Int_t)tree->GetMinimum("ECha");
  // lastSensor = (Int_t)tree->GetMaximum("ECha");

  TClonesArray * array = new TClonesArray("AliDCSSensor",nentries);
   printf ("nentries = %d\n",nentries);

  for (Int_t isensor=0; isensor<nentries; isensor++){
    AliDCSSensor * sens = new ((*array)[isensor])AliDCSSensor;
    tree->GetEntry(isensor);
    sens->SetId(isensor);
    sens->SetIdDCS(idDCS);
    sens->SetStringID(TString(stringId));
    sens->SetX(x);
    sens->SetY(y);
    sens->SetZ(z);

  }
  return array;
}

