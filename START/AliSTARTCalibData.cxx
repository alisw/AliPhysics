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

/* $Id$ */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// class for T0 calibration                       TM-AC-AM_6-02-2006         //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <TCanvas.h>

#include "AliSTARTCalibData.h"
#include "TObjArray.h"
#include "TGraph.h"
#include "TFile.h"
#include "AliLog.h"
#include "TObjString.h"

#include "TAxis.h"
#include "TH2F.h"



ClassImp(AliSTARTCalibData)

//________________________________________________________________
  AliSTARTCalibData::AliSTARTCalibData():   TNamed()

{
  //
}

//________________________________________________________________
AliSTARTCalibData::AliSTARTCalibData(const char* name):TNamed()
{
  TString namst = "Calib_";
  namst += name;
  SetName(namst.Data());
  SetTitle(namst.Data());

}

//________________________________________________________________
AliSTARTCalibData::AliSTARTCalibData(const AliSTARTCalibData& calibda) :
  TNamed(calibda)
{
// copy constructor
  SetName(calibda.GetName());
  SetTitle(calibda.GetName());


}

//________________________________________________________________
AliSTARTCalibData &AliSTARTCalibData::operator =(const AliSTARTCalibData& calibda)
{
// assignment operator
  SetName(calibda.GetName());
  SetTitle(calibda.GetName());
 
  return *this;
}

//________________________________________________________________
AliSTARTCalibData::~AliSTARTCalibData()
{

}
//________________________________________________________________
void AliSTARTCalibData::Reset()
{
    memset(fTimeDelayCFD,1,24*sizeof(Float_t));
    memset(fTimeDelayLED,1,24*sizeof(Float_t));
    memset(fGain,1,24*sizeof(Float_t));
}


//________________________________________________________________
void  AliSTARTCalibData::Print(Option_t*) const
{

  printf("\n	----	PM Arrays	----\n\n");
  printf(" Time delay CFD \n");
  for (Int_t i=0; i<24; i++) printf("  %f",fTimeDelayCFD[i]);
  printf(" \n LED \n");
  for (Int_t i=0; i<24; i++) printf("  %f",fTimeDelayLED[i]);
  printf(" \n Gain \n");
  for (Int_t i=0; i<24; i++) printf("  %f",fGain[i]);
  printf(" \n");
} 



//________________________________________________________________
void AliSTARTCalibData::SetTimeDelayCFD(Float_t* TimeDelay)
{
  if(TimeDelay) for(int t=0; t<24; t++) fTimeDelayCFD[t] = TimeDelay[t];
  //  else for(int t=0; t<24; t++) fTimeDelay[t] = 0.;
}
//________________________________________________________________
void AliSTARTCalibData::SetTimeDelayLED(Float_t* TimeDelay)
{
  if(TimeDelay) for(int t=0; t<24; t++) fTimeDelayLED[t] = TimeDelay[t];
  //  else for(int t=0; t<24; t++) fTimeDelay[t] = 0.;
}

//________________________________________________________________
void AliSTARTCalibData::SetGain(Float_t* Gain)
{
  if(Gain) for(int t=0; t<24; t++) fGain[t] = Gain[t];
  // else for(int t=0; t<24; t++) fGain[t] = 0.;
}


//________________________________________________________________
void AliSTARTCalibData::SetWalk(Int_t ipmt, const Char_t *filename)
{

  TFile *file = new TFile(filename);
  char funcname[256];
  sprintf(funcname,"CFD%i",ipmt+1);
  TF1* gr = (TF1*)file->Get(funcname);
  gr->Print();
  fWalk.AddAtAndExpand(gr,ipmt);
  file->Close();
}


//________________________________________________________________

void AliSTARTCalibData::SetSlewingLED(Int_t ipmt,const Char_t *filename)
{
  Float_t mv, ps; 
  Float_t x[100], y[100];
  cout<<"capacity  "<<fSlewingLED.Capacity()<<endl;
  string buffer;
  
  ifstream inFile(filename);
  if(!inFile) {AliError(Form("Cannot open file %s !",filename));}
  
  inFile >> mv>>ps;
  Int_t i=0;
  
  while(getline(inFile,buffer)){
    x[i]=mv; y[i]=ps;	
    inFile >> mv >> ps;
    i++;
  }
  inFile.close();
  TGraph* gr = new TGraph(i,x,y);
  fSlewingLED.AddAtAndExpand(gr,ipmt);
  cout<<"capacity end "<<fSlewingLED.Capacity()<<endl;

  
}

