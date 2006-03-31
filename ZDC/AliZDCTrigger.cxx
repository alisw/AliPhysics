/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                         *
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

#include "AliLog.h"
#include "AliRun.h"
#include "AliLoader.h"
#include "AliRunLoader.h"
#include "AliTriggerInput.h"

#include "AliZDC.h"
#include "AliZDCDigit.h"
#include "AliZDCTrigger.h"

//________________________________________________________________
ClassImp(AliZDCTrigger)

//________________________________________________________________
AliZDCTrigger::AliZDCTrigger() : AliTriggerDetector() 
{
   SetName("ZDC");
   CreateInputs();
   //
   SetZNMinCut(0);
   SetZDCLeftMinCut(0);
   SetZDCRightMinCut(0);
   SetZEMMinCut(0);
   SetZDCLeftEMDCuts(0,0);
   SetZDCRightEMDCuts(0,0);
   SetZDCLeftMBCut(0);
   SetZDCRightMBCut(0);
   SetZDCLeftCentrCut(0);
   SetZDCRightCentrCut(0);
   SetZDCLeftSemiCentrCut(0);
   SetZDCRightSemiCentrCut(0);
   SetZEMCentrCut(0);

}

//________________________________________________________________
void AliZDCTrigger::CreateInputs()
{
   // inputs 
   
   // Do not create inputs again!!
   if( fInputs.GetEntriesFast() > 0 ) return;
   
   fInputs.AddLast(new AliTriggerInput("ZDC_1_L1",   "ZDC Minimum Bias", 0x01));
   fInputs.AddLast(new AliTriggerInput("ZDC_2_L1",   "ZDC Central", 0x02));
   fInputs.AddLast(new AliTriggerInput("ZDC_3_L1",   "ZDC Semi-central", 0x04));
   fInputs.AddLast(new AliTriggerInput("ZDC_EMD_L1", "ZDC EMD events", 0x08));
}

//________________________________________________________________
void AliZDCTrigger::Trigger()
{


   AliRunLoader *runLoader = gAlice->GetRunLoader();

   AliLoader *ZDCLoader = runLoader->GetLoader("ZDCLoader");
   ZDCLoader->LoadDigits("READ");
   AliZDCDigit digit;
   AliZDCDigit* pdigit = &digit;
   TTree* TD = ZDCLoader->TreeD();
   if (!TD) cerr<<"AliZDCTrigger: digits tree not found\n";
   TD->SetBranchAddress("ZDC", &pdigit);
   //
   Float_t ZNLeftSignal[2], ZPLeftSignal[2], ZDCLeftSumSignal[2];
   Float_t ZNRightSignal[2], ZPRightSignal[2], ZDCRightSumSignal[2];
   Float_t ZEMSumSignal[2];
   for(Int_t iDigit=0; iDigit<TD->GetEntries(); iDigit++){
      TD->GetEntry(iDigit);
      //
      // *** ZDC LEFT
      if(digit.GetSector(0)==1)
         for(Int_t i=0; i<2; i++){ //0=high range; 1=low range
	    ZNLeftSignal[i] += digit.GetADCValue(i);
	    ZDCLeftSumSignal[i] += digit.GetADCValue(i);
	  }
      else if(digit.GetSector(0)==2)
         for(Int_t i=0; i<2; i++){
	    ZPLeftSignal[i] += digit.GetADCValue(i);
	    ZDCLeftSumSignal[i] += digit.GetADCValue(i);
	  }
      else if(digit.GetSector(0)==3)
         for(Int_t i=0; i<2; i++) ZEMSumSignal[i] += digit.GetADCValue(i);
      // *** ZDC RIGHT
      else if(digit.GetSector(0)==4)
         for(Int_t i=0; i<2; i++){ //0=high range; 1=low range
	    ZNRightSignal[i] += digit.GetADCValue(i);
	    ZDCRightSumSignal[i] += digit.GetADCValue(i);
	  }
      else if(digit.GetSector(0)==5)
         for(Int_t i=0; i<2; i++){
	    ZPRightSignal[i] += digit.GetADCValue(i);
	    ZDCRightSumSignal[i] += digit.GetADCValue(i);
	  }
   }
   // *******************************************************************
   if(ZNLeftSignal[0]>fZDCLeftEMDCuts[0] && ZNLeftSignal[0]<fZDCLeftEMDCuts[1] && 
      ZNRightSignal[0]>fZDCRightEMDCuts[0] && ZNRightSignal[0]<fZDCRightEMDCuts[1] &&
      ZEMSumSignal[1]<fZEMMinCut){ // *** 1n EMD trigger
        SetInput("ZDC_EMD_L1");
   }
   // *******************************************************************
   if(ZDCLeftSumSignal[1]>fZDCLeftMBCut && ZDCRightSumSignal[1]>fZDCRightMBCut) 
       // *** ZDC minimum bias trigger
       SetInput("ZDC_1_L1");
   // *******************************************************************
   if(ZDCLeftSumSignal[1]>fZDCLeftMinCut && ZDCLeftSumSignal[1]<fZDCLeftCentrCut &&
      ZDCRightSumSignal[1]>fZDCRightMinCut && ZDCRightSumSignal[1]<fZDCRightCentrCut &&
      ZEMSumSignal[1]>fZEMCentrCut) 
       // *** ZDC central (0-10%)
       SetInput("ZDC_2_L1");
   // *******************************************************************
   if(ZDCLeftSumSignal[1]>fZDCLeftCentrCut && ZDCLeftSumSignal[1]<fZDCLeftSemiCentrCut &&
      ZDCRightSumSignal[1]>fZDCRightCentrCut && ZDCRightSumSignal[1]<fZDCRightSemiCentrCut
      && ZEMSumSignal[1]>fZEMCentrCut) 
       // *** ZDC semi-central (10-40%)
       SetInput("ZDC_3_L1");
   
}

//________________________________________________________________
void AliZDCTrigger::SetZNMinCut(Float_t ZNMinCut) 
{
  if(ZNMinCut)  fZNMinCut = ZNMinCut;
  else  fZNMinCut = 400.;
}

//________________________________________________________________
void AliZDCTrigger::SetZDCLeftMinCut(Float_t ZDCLeftMinCut) 
{
  if(ZDCLeftMinCut)  fZDCLeftMinCut = ZDCLeftMinCut;
  else  fZDCLeftMinCut = 800.;
}
//________________________________________________________________
void AliZDCTrigger::SetZDCRightMinCut(Float_t ZDCRightMinCut) 
{
  if(ZDCRightMinCut)  fZDCRightMinCut = ZDCRightMinCut;
  else  fZDCRightMinCut = 800.;
}

//________________________________________________________________
void AliZDCTrigger::SetZEMMinCut(Float_t ZEMMinCut) 
{
  if(ZEMMinCut)  fZEMMinCut = ZEMMinCut;
  else  fZEMMinCut = 80.;
}
//________________________________________________________________
void AliZDCTrigger::SetZDCLeftEMDCuts(Float_t* ZDCLeftEMDCuts) 
{
  if(ZDCLeftEMDCuts) for(int j=0; j<2; j++) fZDCLeftEMDCuts[j] = ZDCLeftEMDCuts[j];
  else{
    fZDCLeftEMDCuts[0] = 600.;
    fZDCLeftEMDCuts[1] = 1000.;
  }
}
//________________________________________________________________
void AliZDCTrigger::SetZDCLeftEMDCuts(Float_t ZDCLeftEMDCutInf, 
	Float_t ZDCLeftEMDCutSup) 
{
  if(ZDCLeftEMDCutInf && ZDCLeftEMDCutSup){
    fZDCLeftEMDCuts[0]=ZDCLeftEMDCutInf; 
    fZDCLeftEMDCuts[1]=ZDCLeftEMDCutSup;
  } 	
  else{
    fZDCLeftEMDCuts[0] = 600.;
    fZDCLeftEMDCuts[1] = 1000.;
  }
}
//________________________________________________________________
void AliZDCTrigger::SetZDCRightEMDCuts(Float_t* ZDCRightEMDCuts) 
{
  if(ZDCRightEMDCuts) for(int j=0; j<2; j++) fZDCRightEMDCuts[j] = ZDCRightEMDCuts[j];
  else{
    fZDCRightEMDCuts[0] = 600.;
    fZDCRightEMDCuts[1] = 1000.;
  }
}
//________________________________________________________________
void AliZDCTrigger::SetZDCRightEMDCuts(Float_t ZDCRightEMDCutInf, 
	Float_t ZDCRightEMDCutSup) 
{
  if(ZDCRightEMDCutInf && ZDCRightEMDCutSup){
    fZDCRightEMDCuts[0]=ZDCRightEMDCutInf; 
    fZDCRightEMDCuts[1]=ZDCRightEMDCutSup;
  } 	
  else{
    fZDCRightEMDCuts[0] = 600.;
    fZDCRightEMDCuts[1] = 1000.;
  }
}
//________________________________________________________________
void AliZDCTrigger::SetZDCLeftMBCut(Float_t ZDCLeftMBCut) 
{
  if(ZDCLeftMBCut) fZDCLeftMBCut = ZDCLeftMBCut;
  else fZDCLeftMBCut = 800.;
}
//________________________________________________________________
void AliZDCTrigger::SetZDCRightMBCut(Float_t ZDCRightMBCut) 
{
  if(ZDCRightMBCut) fZDCRightMBCut = ZDCRightMBCut;
  else fZDCRightMBCut = 800.;
}
//________________________________________________________________
void AliZDCTrigger::SetZDCLeftCentrCut(Float_t ZDCLeftCentrCut) 
{
  if(ZDCLeftCentrCut) fZDCLeftCentrCut = ZDCLeftCentrCut;
  else fZDCLeftCentrCut = 10000.;
}
//________________________________________________________________
void AliZDCTrigger::SetZDCRightCentrCut(Float_t ZDCRightCentrCut) 
{
  if(ZDCRightCentrCut) fZDCRightCentrCut = ZDCRightCentrCut;
  else fZDCRightCentrCut = 10000.;
}
//________________________________________________________________
void AliZDCTrigger::SetZDCLeftSemiCentrCut(Float_t ZDCLeftSemiCentrCut) 
{
  if(ZDCLeftSemiCentrCut) fZDCLeftSemiCentrCut = ZDCLeftSemiCentrCut;
  else fZDCLeftSemiCentrCut = 18500.;
}
//________________________________________________________________
void AliZDCTrigger::SetZDCRightSemiCentrCut(Float_t ZDCRightSemiCentrCut) 
{
  if(ZDCRightSemiCentrCut) fZDCRightSemiCentrCut = ZDCRightSemiCentrCut;
  else fZDCRightSemiCentrCut = 18500.;
}
//________________________________________________________________
void AliZDCTrigger::SetZEMCentrCut(Float_t ZEMCentrCut) 
{
  if(ZEMCentrCut) fZEMCentrCut = ZEMCentrCut;
  else fZEMCentrCut = 210.;
}
