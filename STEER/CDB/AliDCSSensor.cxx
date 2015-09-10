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
// Class describing time dependent values read from DCS sensors               //  
// (including pointers to graphs/fits)                                        //
// Authors: Marian Ivanov, Haavard Helstrup and Martin Siska                  //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////


#include "AliDCSSensor.h"
#include "TDatime.h"
#include "TCanvas.h"
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
   fGraph(0),
   fFit(0),
   fX(source.fX),
   fY(source.fY),
   fZ(source.fZ)
//
//  Copy constructor
//
{ 
   if (source.fGraph) fGraph = (TGraph*)source.fGraph->Clone();
   if (source.fFit) fFit = (AliSplineFit*)source.fFit->Clone();
}

AliDCSSensor::~AliDCSSensor(){
  //
  // Destructor
  //
  if(fGraph)
    delete fGraph;
  fGraph=0;

}

AliDCSSensor& AliDCSSensor::operator=(const AliDCSSensor& source){
//
// assignment operator
//
  if (&source == this) return *this;
  new (this) AliDCSSensor(source);

  return *this;
}


void AliDCSSensor::Print(const Option_t* option) const{
  //
  // print function
  //  
  TString opt = option; opt.ToLower();
  printf("%s:%s\n",GetTitle(), GetName());
  printf("%s\n",fStringID.Data());

}

void AliDCSSensor::Draw(Option_t* option) {
  //
  // draw function - to viusalize sensor
  // Unfortuantelly - it  make a memory leak as function Draw does not return the object pointer
  //
  TCanvas * canvas = new TCanvas((fStringID+option).Data(), (fStringID+option).Data()); 
  if (fGraph){
    // transform points to time in s
    Int_t npoints = fGraph->GetN();
    for (Int_t i=0; i<npoints; i++){
      fGraph->GetX()[i]=fGraph->GetX()[i]*3600+fStartTime;
    }
    fGraph->Draw("alp");
    return;
  }
  canvas->cd();
  TGraph * graph = MakeGraph(100);  // memory leak - we can not modify the content - const method
  if (graph){
    graph->Draw(option);              // 
  }
}



//_____________________________________________________________________________
Double_t AliDCSSensor::GetValue(UInt_t timeSec)
{
 // 
 // Get DCS value for actual sensor
 //  timeSec given as offset from start-of-map measured in seconds
 //  *NOTE* In the current TPC setup, start-of-map is defined as the 
 //         first measured point for each sensor. This will be different
 //         for each sensor in the array. If you want to get a value at the 
 //         same absolute time, use AliDCSSensor::GetValue(TTimeStamp time)
 //         or AliDCSSensorArray::GetValue (UInt_t timeSec, Int_t sensor)
 //         which measure offsets with respect to the (global) start-of-run
 //
 Bool_t inside=kTRUE;
 return Eval(TTimeStamp((time_t)(fStartTime+timeSec),0),inside);
}
//_____________________________________________________________________________
Double_t AliDCSSensor::GetValue(TTimeStamp time) 
{
 // Get DCS value for actual sensor
 //  time given as absolute TTimeStamp
 //
 Bool_t inside=kTRUE;
 return Eval(time, inside);
}

//_____________________________________________________________________________

Double_t AliDCSSensor::Eval(const TTimeStamp& time, Bool_t& inside) const
{
  // 
  // Return DCS value at given time
  //  The value is calculated from the AliSplineFit, if a fit is not available 
  //    the most recent reading from the Graph of DCS points is returned (if 
  //    the graph is present)
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
     if ( fGraph ) {
       return EvalGraph(timeHour);
     } else {  
       return -99;
     }
  }
}
//_____________________________________________________________________________

Double_t AliDCSSensor::EvalGraph(const TTimeStamp& time, Bool_t& inside) const
{
  // 
  // Return DCS value from graph of DCS points (i.e return last reading before
  //  the time specified by TTimeStamp
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
  if ( fGraph ) {
     return EvalGraph(timeHour);
  } else {  
     return -99;
  }  
}
//_____________________________________________________________________________
Double_t AliDCSSensor::EvalGraph(const Double_t& timeHour) const 
{
  //
  // Extract last value in graph observed before time given by timeHour
  //

  // return -99 if point specified is before beginning of graph
  Double_t x=0; Double_t y=0;
  fGraph->GetPoint(0,x,y);
  if ( timeHour < x ) return -99;
  
  // return previous point when first time > timeHour is observed
  
  Int_t npoints = fGraph->GetN();
  for (Int_t i=1; i<npoints; i++) {
     fGraph->GetPoint(i,x,y);
     if ( timeHour < x ) {
       fGraph->GetPoint(i-1,x,y);
       return y;
     }
  }
  
  // return last point if all times are < timeHour
  return y;
} 
	

//_____________________________________________________________________________
TGraph* AliDCSSensor::MakeGraph(Int_t nPoints, Bool_t debug) const
{
  //
  // Make graph from start time to end time of DCS values 
  //

 

  UInt_t stepTime = (fEndTime-fStartTime)/nPoints;
  
  if (debug==kTRUE) {
     printf ("Start time %d, End time %d, step time %d\n",
     fStartTime,fEndTime,stepTime);
     TTimeStamp t((time_t)fStartTime,0); t.Print();
     TTimeStamp t2((time_t)fEndTime,0); t2.Print();
  }     
  
  if ( !fFit ) return 0;

  Double_t *x = new Double_t[nPoints+1];
  Double_t *y = new Double_t[nPoints+1];
  for (Int_t ip=0; ip<nPoints; ip++) {
    x[ip] = (time_t)(fStartTime+ip*stepTime);
    y[ip] = fFit->Eval(ip*stepTime/kSecInHour);
    if (debug==kTRUE) {
     TTimeStamp t3((time_t)x[ip],0); 
     printf ("x=%f, y=%f  ",x[ip],y[ip]);
     t3.Print();
    }
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

