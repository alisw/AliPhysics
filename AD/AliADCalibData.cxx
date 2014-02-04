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

/* $Id: AliADCalibData.cxx,                                            */


#include "AliADCalibData.h"
#include "TList.h"
#include "TCanvas.h"

ClassImp(AliADCalibData)


//________________________________________________________________
AliADCalibData::AliADCalibData()
{
 	for (Int_t imod = 0; imod < 60; imod++)
	{
		fEfficiencies[imod]=0.;
		fRates[imod]=0.;
		fModulesActivity[imod]=0.;
	} 
}

//________________________________________________________________
void AliADCalibData::Reset()
{
  
}

//________________________________________________________________
AliADCalibData::AliADCalibData(const char* name) :
  TNamed()
{
  TString namst = "Calib_";
  namst += name;
  SetName(namst.Data());
  SetTitle(namst.Data());

}

//________________________________________________________________
AliADCalibData::AliADCalibData(const AliADCalibData& calibda) :
  TNamed(calibda)
{
// copy constructor

  SetName(calibda.GetName());
  SetTitle(calibda.GetName());
  
  // there are 60 modules. Note that number of first module is 1 (one)
  for(int t=0; t<60; t++) 
  {
  	fEfficiencies[t] =calibda.GetEfficiency(t);
  	fRates[t] = calibda.GetRate(t);
	fModulesActivity[t] = calibda.GetModuleActivity(t);
  }
}
//_______________________________________________________________
void AliADCalibData::Draw(Option_t *)
{
 

  //fHits->Draw();

 
}
//________________________________________________________________
AliADCalibData &AliADCalibData::operator =(const AliADCalibData& calibda)
{
// assignment operator

  SetName(calibda.GetName());
  SetTitle(calibda.GetName());
  // there are 60 modules. Note that number of first module is 1 (one)
  for(int t=0; t<60; t++) 
  {
  	fEfficiencies[t] =calibda.GetEfficiency(t);
  	fRates[t] = calibda.GetRate(t);
	fModulesActivity[t] = calibda.GetModuleActivity(t);
  }
  return *this;
}
//_______________________________________________________________
/*void AliADCalibData::AddHisto(TH1D *fHist)
{
    


 = (TH1D*)fHist->Clone("hnew");

     
   
 
}
*/

//________________________________________________________________
AliADCalibData::~AliADCalibData()
{
  
}

                                                                                   

//________________________________________________________________
void AliADCalibData::SetEfficiencies(Float_t* Eff)
{
  // there are 60 modules. Note that number of first module is 1 (one)
  if(Eff) for(int t=0; t<60; t++) fEfficiencies[t] = Eff[t];
  else for(int t=0; t<60; t++) fEfficiencies[t] = 0.0;
}

void AliADCalibData::SetRates(Float_t* Rt)
{
   if(Rt) for (int t=0;t<60; t++) fRates[t] = Rt[t];
else for (int t=0;t<60; t++) fRates[t] = 0.0;
}

void AliADCalibData::SetModulesActivity(Float_t* Mac)
{
	if(Mac) for (int t=0;t<60;t++) fModulesActivity[t] = Mac[t];
	else for (int t=0;t<60;t++) fModulesActivity[t] = 0.0;
}

