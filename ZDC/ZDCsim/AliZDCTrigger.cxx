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

// ****************************************************************
//
//	Trigger class for ZDC
//
// ****************************************************************

#include <TTree.h>
#include "AliLog.h"
#include "AliRun.h"
#include "AliLoader.h"
#include "AliRunLoader.h"
#include "AliTriggerInput.h"

#include "AliZDC.h"
#include "AliZDCDigit.h"
#include "AliZDCTrigger.h"
#include <iostream>

//________________________________________________________________
using std::cerr;
ClassImp(AliZDCTrigger)

//________________________________________________________________
AliZDCTrigger::AliZDCTrigger() : 
   AliTriggerDetector(), 
   fZDCLeftMinCut(0),
   fZDCRightMinCut(0),
   fZEMMinCut(0),
   fZDCLeftMBCut(0),
   fZDCRightMBCut(0),
   fZDCLeftCentrCut(0),
   fZDCRightCentrCut(0),
   fZDCLeftSemiCentrCut(0),
   fZDCRightSemiCentrCut(0),
   fZEMCentrCut(0)
{  
   // Constructor
   SetName("ZDC");
   CreateInputs();
   //
   SetZDCLeftEMDCuts(0,0);
   SetZDCRightEMDCuts(0,0);

}

//________________________________________________________________
void AliZDCTrigger::CreateInputs()
{
   // Trigger inputs
   
   // Do not create inputs again!!
   if( fInputs.GetEntriesFast() > 0 ) return;
   
   fInputs.AddLast(new AliTriggerInput("ZDC_1_L1",   "ZDC", 1));
   fInputs.AddLast(new AliTriggerInput("ZDC_2_L1", "ZDC", 1));
   fInputs.AddLast(new AliTriggerInput("ZDC_3_L1", "ZDC", 1));
   fInputs.AddLast(new AliTriggerInput("ZDC_EMD_L1",  "ZDC", 1));
}

//________________________________________________________________
void AliZDCTrigger::Trigger()
{

   // Trigger selection
   //
   AliRunLoader *runLoader = AliRunLoader::Instance();

   AliLoader *aZDCLoader = runLoader->GetLoader("ZDCLoader");
   
   aZDCLoader->LoadDigits("READ");
   AliZDCDigit digit;
   AliZDCDigit* pdigit = &digit;
   TTree* tD = aZDCLoader->TreeD();
   if (!tD) {
     cerr<<"AliZDCTrigger: digits tree not found\n";
     return;
   }
   tD->SetBranchAddress("ZDC", &pdigit);
   //
   Float_t signalZNLeft[]={0,0}, signalZPLeft[]={0,0}, signalZDCLeftSum[]={0,0};
   Float_t signalZNRight[]={0,0}, signalZPRight[]={0,0}, signalZDCRightSum[]={0,0};
   Float_t signalZEMSum[]={0,0};
   for(Int_t iDigit=0; iDigit<tD->GetEntries(); iDigit++){
      tD->GetEntry(iDigit);
      //
      // *** ZDC LEFT
      if(digit.GetSector(0)==1)
         for(Int_t i=0; i<2; i++){ //0=high range; 1=low range
	    signalZNLeft[i] += digit.GetADCValue(i);
	    signalZDCLeftSum[i] += digit.GetADCValue(i);
	  }
      else if(digit.GetSector(0)==2)
         for(Int_t i=0; i<2; i++){
	    signalZPLeft[i] += digit.GetADCValue(i);
	    signalZDCLeftSum[i] += digit.GetADCValue(i);
	  }
      else if(digit.GetSector(0)==3)
         for(Int_t i=0; i<2; i++) signalZEMSum[i] += digit.GetADCValue(i);
      // *** ZDC RIGHT
      else if(digit.GetSector(0)==4)
         for(Int_t i=0; i<2; i++){ //0=high range; 1=low range
	    signalZNRight[i] += digit.GetADCValue(i);
	    signalZDCRightSum[i] += digit.GetADCValue(i);
	  }
      else if(digit.GetSector(0)==5)
         for(Int_t i=0; i<2; i++){
	    signalZPRight[i] += digit.GetADCValue(i);
	    signalZDCRightSum[i] += digit.GetADCValue(i);
	  }
   }
   // *******************************************************************
   if(signalZDCLeftSum[1]>fZDCLeftMBCut && signalZDCRightSum[1]>fZDCRightMBCut) 
       // *** ZDC minimum bias trigger
       SetInput("ZDC_1_L1");
   // *******************************************************************
   if(signalZDCLeftSum[1]>fZDCLeftCentrCut && signalZDCLeftSum[1]<fZDCLeftSemiCentrCut &&
      signalZDCRightSum[1]>fZDCRightCentrCut && signalZDCRightSum[1]<fZDCRightSemiCentrCut
      && signalZEMSum[1]>fZEMCentrCut) 
       // *** ZDC semi-central (10-40%)
       SetInput("ZDC_2_L1");
   // *******************************************************************
   if(signalZDCLeftSum[1]>fZDCLeftMinCut && signalZDCLeftSum[1]<fZDCLeftCentrCut &&
      signalZDCRightSum[1]>fZDCRightMinCut && signalZDCRightSum[1]<fZDCRightCentrCut &&
      signalZEMSum[1]>fZEMCentrCut) 
       // *** ZDC central (0-10%)
       SetInput("ZDC_3_L1");
   // *******************************************************************
   if(signalZNLeft[0]>fZDCLeftEMDCuts[0] && signalZNLeft[0]<fZDCLeftEMDCuts[1] && 
      signalZNRight[0]>fZDCRightEMDCuts[0] && signalZNRight[0]<fZDCRightEMDCuts[1] &&
      signalZEMSum[1]<fZEMMinCut){ // *** 1n EMD trigger
        SetInput("ZDC_EMD_L1");
   }
   
}

//________________________________________________________________
void AliZDCTrigger::SetZDCLeftMinCut(Float_t ZDCLeftMinCut) 
{
  // Set default cut values for ZDC trigger
  //
  if(ZDCLeftMinCut)  fZDCLeftMinCut = ZDCLeftMinCut;
  else  fZDCLeftMinCut = 800.;
}
//________________________________________________________________
void AliZDCTrigger::SetZDCRightMinCut(Float_t ZDCRightMinCut) 
{
  // Set default cut values for ZDC trigger
  //
  if(ZDCRightMinCut)  fZDCRightMinCut = ZDCRightMinCut;
  else  fZDCRightMinCut = 800.;
}

//________________________________________________________________
void AliZDCTrigger::SetZEMMinCut(Float_t ZEMMinCut) 
{
  // Set default cut values for ZDC trigger
  //
  if(ZEMMinCut)  fZEMMinCut = ZEMMinCut;
  else  fZEMMinCut = 80.;
}
//________________________________________________________________
void AliZDCTrigger::SetZDCLeftEMDCuts(Float_t* ZDCLeftEMDCuts) 
{
  // Set default cut values for ZDC trigger
  //
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
  // Set default cut values for ZDC trigger
  //
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
  // Set default cut values for ZDC trigger
  //
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
  // Set default cut values for ZDC trigger
  //
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
  // Set default cut values for ZDC trigger
  //
  if(ZDCLeftMBCut) fZDCLeftMBCut = ZDCLeftMBCut;
  else fZDCLeftMBCut = 800.;
}
//________________________________________________________________
void AliZDCTrigger::SetZDCRightMBCut(Float_t ZDCRightMBCut) 
{
  // Set default cut values for ZDC trigger
  //
  if(ZDCRightMBCut) fZDCRightMBCut = ZDCRightMBCut;
  else fZDCRightMBCut = 800.;
}
//________________________________________________________________
void AliZDCTrigger::SetZDCLeftCentrCut(Float_t ZDCLeftCentrCut) 
{
  // Set default cut values for ZDC trigger
  //
  if(ZDCLeftCentrCut) fZDCLeftCentrCut = ZDCLeftCentrCut;
  else fZDCLeftCentrCut = 10000.;
}
//________________________________________________________________
void AliZDCTrigger::SetZDCRightCentrCut(Float_t ZDCRightCentrCut) 
{
  // Set default cut values for ZDC trigger
  //
  if(ZDCRightCentrCut) fZDCRightCentrCut = ZDCRightCentrCut;
  else fZDCRightCentrCut = 10000.;
}
//________________________________________________________________
void AliZDCTrigger::SetZDCLeftSemiCentrCut(Float_t ZDCLeftSemiCentrCut) 
{
  // Set default cut values for ZDC trigger
  //
  if(ZDCLeftSemiCentrCut) fZDCLeftSemiCentrCut = ZDCLeftSemiCentrCut;
  else fZDCLeftSemiCentrCut = 18500.;
}
//________________________________________________________________
void AliZDCTrigger::SetZDCRightSemiCentrCut(Float_t ZDCRightSemiCentrCut) 
{
  // Set default cut values for ZDC trigger
  //
  if(ZDCRightSemiCentrCut) fZDCRightSemiCentrCut = ZDCRightSemiCentrCut;
  else fZDCRightSemiCentrCut = 18500.;
}
//________________________________________________________________
void AliZDCTrigger::SetZEMCentrCut(Float_t ZEMCentrCut) 
{
  // Set default cut values for ZDC trigger
  //
  if(ZEMCentrCut) fZEMCentrCut = ZEMCentrCut;
  else fZEMCentrCut = 210.;
}
