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

/* $Id: AliACORDECalibData.cxx,                                            */


#include "AliACORDECalibData.h"
#include "TList.h"
#include "TCanvas.h"

ClassImp(AliACORDECalibData)


//________________________________________________________________
AliACORDECalibData::AliACORDECalibData():
fHits(0),
fTHits(0),
fMultiHits(0),
fTMultiHits(0)
{
  
}

//________________________________________________________________
void AliACORDECalibData::Reset()
{
  
}

//________________________________________________________________
AliACORDECalibData::AliACORDECalibData(const char* name) :
  TNamed(),
  fHits(0),
  fTHits(0),
  fMultiHits(0),
  fTMultiHits(0)
{
  TString namst = "Calib_";
  namst += name;
  SetName(namst.Data());
  SetTitle(namst.Data());

}

//________________________________________________________________
AliACORDECalibData::AliACORDECalibData(const AliACORDECalibData& calibda) :
  TNamed(calibda),
  fHits(0),
  fTHits(0),
  fMultiHits(0),
  fTMultiHits(0)
{
// copy constructor

  SetName(calibda.GetName());
  SetTitle(calibda.GetName());
  
  // there are 60 modules. Note that number of first module is 1 (one)
  for(int t=0; t<60; t++) 
  {
  	fEfficiencies[t] =calibda.GetEfficiency(t+1);
  	fRates[t] = calibda.GetRate(t+1);
	fModulesActivity[t] = calibda.GetModuleActivity(t+1);
  }
}
//_______________________________________________________________
void AliACORDECalibData::Draw(Option_t *)
{
 

  //fHits->Draw();



  TCanvas *ch;
  TString canvasHistoName="Histos";
  ch=new TCanvas(canvasHistoName,canvasHistoName,20,20,600,600);
  ch->Divide(2,2);
  ch->cd(1);
  fHits->Draw();
  ch->cd(2);
  fTHits->Draw();
  ch->cd(3);
  fMultiHits->Draw();
  ch->cd(4);
  fTMultiHits->Draw();

 
}
//________________________________________________________________
AliACORDECalibData &AliACORDECalibData::operator =(const AliACORDECalibData& calibda)
{
// assignment operator

  SetName(calibda.GetName());
  SetTitle(calibda.GetName());
  // there are 60 modules. Note that number of first module is 1 (one)
  for(int t=0; t<60; t++) 
  {
  	fEfficiencies[t] =calibda.GetEfficiency(t+1);
  	fRates[t] = calibda.GetRate(t+1);
	fModulesActivity[t] = calibda.GetModuleActivity(t+1);
  }
  return *this;
}
//_______________________________________________________________
/*void AliACORDECalibData::AddHisto(TH1D *fHist)
{
    


 = (TH1D*)fHist->Clone("hnew");

     
   
 
}
*/

//________________________________________________________________
AliACORDECalibData::~AliACORDECalibData()
{
  
}

                                                                                   

//________________________________________________________________
void AliACORDECalibData::SetEfficiencies(Float_t* Eff)
{
  // there are 60 modules. Note that number of first module is 1 (one)
  if(Eff) for(int t=0; t<60; t++) fEfficiencies[t] = Eff[t];
  else for(int t=0; t<60; t++) fEfficiencies[t] = 0.0;
}

void AliACORDECalibData::SetRates(Float_t* Rt)
{
   if(Rt) for (int t=0;t<60; t++) fRates[t] = Rt[t];
else for (int t=0;t<60; t++) fRates[t] = 0.0;
}

void AliACORDECalibData::SetModulesActivity(Float_t* Mac)
{
	if(Mac) for (int t=0;t<60;t++) fModulesActivity[t] = Mac[t];
	else for (int t=0;t<60;t++) fModulesActivity[t] = 0.0;
}

