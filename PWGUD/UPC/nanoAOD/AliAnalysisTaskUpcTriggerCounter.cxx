/*************************************************************************
* Copyright(c) 1998-2008, ALICE Experiment at CERN, All rights reserved. *
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

// c++ headers
#include <iostream>
#include <string.h>

// root headers
#include "TList.h"
#include "TH1I.h"
#include "TH2I.h"
#include "TString.h"

// aliroot headers
#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"
#include "AliAODEvent.h"

// my headers
#include "AliAnalysisTaskUpcTriggerCounter.h"

ClassImp(AliAnalysisTaskUpcTriggerCounter);

TString gTriggerName[ntrig] = {"CCUP4-B","CCUP2-B","CCUP7-B","CINT1-B","CTEST58-B", "CTEST59-B", "CTEST60-B", "CTEST61-B", "CCUP8-B", "CCUP9-B", "CCUP10-B", "CCUP11-B", "CCUP12-B"};

using std::cout;
using std::endl;

//Task counting upc triggers. To be used in nanoAOD train 
// michal.broz@cern.ch

//_____________________________________________________________________________
AliAnalysisTaskUpcTriggerCounter::AliAnalysisTaskUpcTriggerCounter() 
  : AliAnalysisTaskSE(),fListHist(0),fHistTriggerCounter(0),fHistTriggerCounterIR1(0),fHistTriggerCounterIR2(0)

{

//Dummy constructor

}//AliAnalysisTaskUpcTriggerCounter


//_____________________________________________________________________________
AliAnalysisTaskUpcTriggerCounter::AliAnalysisTaskUpcTriggerCounter(const char *name) 
  : AliAnalysisTaskSE(name),fListHist(0),fHistTriggerCounter(0),fHistTriggerCounterIR1(0),fHistTriggerCounterIR2(0)

{

  DefineOutput(1, TList::Class());

}//AliAnalysisTaskUpcTriggerCounter



//_____________________________________________________________________________
AliAnalysisTaskUpcTriggerCounter::~AliAnalysisTaskUpcTriggerCounter() 
{
  // Destructor
  if(fListHist){
     delete fListHist;
     fListHist = 0x0;
  }

}//~AliAnalysisTaskUpcTriggerCounter


//_____________________________________________________________________________
void AliAnalysisTaskUpcTriggerCounter::UserCreateOutputObjects()
{
  
  fListHist = new TList();
  fListHist ->SetOwner();
  
  fHistTriggerCounter = new TH1I("fHistTriggerCounter","Number of recored UPC triggers",ntrig,0,ntrig);
  fHistTriggerCounterIR1 = new TH2I("fHistTriggerCounterIR1","Number of recored UPC triggers",ntrig,0,ntrig,91,-0.5,90.5);
  fHistTriggerCounterIR2 = new TH2I("fHistTriggerCounterIR2","Number of recored UPC triggers",ntrig,0,ntrig,91,-0.5,90.5);

  for (Int_t i = 0; i<ntrig; i++){ 
  	fHistTriggerCounter->GetXaxis()->SetBinLabel(i+1,gTriggerName[i].Data());
	fHistTriggerCounterIR1->GetXaxis()->SetBinLabel(i+1,gTriggerName[i].Data());
	fHistTriggerCounterIR2->GetXaxis()->SetBinLabel(i+1,gTriggerName[i].Data());
	}
  
  fListHist->Add(fHistTriggerCounter); 
  fListHist->Add(fHistTriggerCounterIR1);
  fListHist->Add(fHistTriggerCounterIR2);
  
  PostData(1, fListHist);

}//UserCreateOutputObjects


//_____________________________________________________________________________
void AliAnalysisTaskUpcTriggerCounter::UserExec(Option_t *) 
{

  
  AliAODEvent *aod = (AliAODEvent*) InputEvent();
  if(!aod) return;
  
  TBits fIR1Map = aod->GetHeader()->GetIRInt1InteractionMap();
  TBits fIR2Map = aod->GetHeader()->GetIRInt2InteractionMap();
  Int_t fClosestIR1 = 100;
  Int_t fClosestIR2 = 100;
  for(Int_t item=-1; item>=-90; item--) {
    Int_t bin = 90+item;
    Bool_t isFired = fIR1Map.TestBitNumber(bin);
    if(isFired) {
      fClosestIR1 = TMath::Abs(item);
      break;
    }
  if(fClosestIR1 == 100)fClosestIR1 = 0;
  }
  for(Int_t item=-1; item>=-90; item--) {
    Int_t bin = 90+item;
    Bool_t isFired = fIR2Map.TestBitNumber(bin);
    if(isFired) {
      fClosestIR2 = TMath::Abs(item);
      break;
    }
  }
  if(fClosestIR2 == 100)fClosestIR2 = 0;
	
  TString trigger = aod->GetFiredTriggerClasses();
  for (Int_t i = 0; i<ntrig; i++){
  	if(trigger.Contains(gTriggerName[i].Data())){
		fHistTriggerCounter->Fill(i);
		fHistTriggerCounterIR1->Fill(i,fClosestIR1);
		fHistTriggerCounterIR2->Fill(i,fClosestIR2);
		}
	}
	
  PostData(1, fListHist);

}//UserExec


//_____________________________________________________________________________
void AliAnalysisTaskUpcTriggerCounter::Terminate(Option_t *) 
{

  cout<<"Analysis complete."<<endl;
}//Terminate
