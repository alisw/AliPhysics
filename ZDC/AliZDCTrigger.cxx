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
   SetZDCMinCut(0);
   SetZEMMinCut(0);
   SetZDCLeftEMDCuts(0,0);
   SetZDCRightEMDCuts(0,0);
   SetZDCMBCut(0);
   SetZDCCentrCut(0);
   SetZDCSemiCentrCut(0);
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
   Float_t ZNSignal[2], ZPSignal[2], ZDCSumSignal[2], ZEMSumSignal[2];
   for(Int_t iDigit=0; iDigit<TD->GetEntries(); iDigit++){
      TD->GetEntry(iDigit);
      //
      if(digit.GetSector(0)==1)
         for(Int_t i=0; i<2; i++){
	    ZNSignal[i] += digit.GetADCValue(i);
	    ZDCSumSignal[i] += digit.GetADCValue(i);
	  }
      else if(digit.GetSector(0)==2)
         for(Int_t i=0; i<2; i++){
	    ZPSignal[i] += digit.GetADCValue(i);
	    ZDCSumSignal[i] += digit.GetADCValue(i);
	  }
      else if(digit.GetSector(0)==3)
         for(Int_t i=0; i<2; i++) ZEMSumSignal[i] += digit.GetADCValue(i);
   }
   // *******************************************************************
   if(ZNSignal[1]>fZDCLeftEMDCuts[0] && ZNSignal[1]<fZDCLeftEMDCuts[1] && 
   	ZEMSumSignal[1]<fZEMMinCut){ // *** n EMD event
       SetInput("ZDC_EMD_L1");
   }
   // *******************************************************************
   if(ZDCSumSignal[0]>fZDCMBCut) // *** ZDC minimum bias
     SetInput("ZDC_1_L1");
   // *******************************************************************
   if(ZDCSumSignal[0]>fZDCMinCut && ZDCSumSignal[0]<fZDCMBCut) 
     // *** ZDC central (0-10%)
     SetInput("ZDC_2_L1");
   // *******************************************************************
   if(ZDCSumSignal[0]>fZDCCentrCut && ZDCSumSignal[0]<fZDCSemiCentrCut
   	   && ZEMSumSignal[0]>fZEMCentrCut) 
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
void AliZDCTrigger::SetZDCMinCut(Float_t ZDCMinCut) 
{
  if(ZDCMinCut)  fZDCMinCut = ZDCMinCut;
  else  fZDCMinCut = 800.;
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
void AliZDCTrigger::SetZDCMBCut(Float_t ZDCMBCut) 
{
  if(ZDCMBCut) fZDCMBCut = ZDCMBCut;
  else fZDCMBCut = 800.;
}
//________________________________________________________________
void AliZDCTrigger::SetZDCCentrCut(Float_t ZDCCentrCut) 
{
  if(ZDCCentrCut) fZDCCentrCut = ZDCCentrCut;
  else fZDCCentrCut = 10000.;
}
//________________________________________________________________
void AliZDCTrigger::SetZDCSemiCentrCut(Float_t ZDCSemiCentrCut) 
{
  if(ZDCSemiCentrCut) fZDCSemiCentrCut = ZDCSemiCentrCut;
  else fZDCSemiCentrCut = 18500.;
}
//________________________________________________________________
void AliZDCTrigger::SetZEMCentrCut(Float_t ZEMCentrCut) 
{
  if(ZEMCentrCut) fZEMCentrCut = ZEMCentrCut;
  else fZEMCentrCut = 210.;
}
