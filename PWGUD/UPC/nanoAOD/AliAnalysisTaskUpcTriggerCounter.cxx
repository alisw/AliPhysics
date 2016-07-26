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
#include "TH1I.h"
#include "TString.h"

// aliroot headers
#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"
#include "AliAODEvent.h"

// my headers
#include "AliAnalysisTaskUpcTriggerCounter.h"

ClassImp(AliAnalysisTaskUpcTriggerCounter);

using std::cout;
using std::endl;

//Task counting upc triggers. To be used in nanoAOD train 
// michal.broz@cern.ch

//_____________________________________________________________________________
AliAnalysisTaskUpcTriggerCounter::AliAnalysisTaskUpcTriggerCounter() 
  : AliAnalysisTaskSE(),fHistTriggerCounter(0)

{

//Dummy constructor

}//AliAnalysisTaskUpcTriggerCounter


//_____________________________________________________________________________
AliAnalysisTaskUpcTriggerCounter::AliAnalysisTaskUpcTriggerCounter(const char *name) 
  : AliAnalysisTaskSE(name),fHistTriggerCounter(0)

{

  DefineOutput(1, TH1I::Class());

}//AliAnalysisTaskUpcTriggerCounter



//_____________________________________________________________________________
AliAnalysisTaskUpcTriggerCounter::~AliAnalysisTaskUpcTriggerCounter() 
{
  // Destructor
  if(fHistTriggerCounter){
     delete fHistTriggerCounter;
     fHistTriggerCounter = 0x0;
  }

}//~AliAnalysisTaskUpcTriggerCounter


//_____________________________________________________________________________
void AliAnalysisTaskUpcTriggerCounter::UserCreateOutputObjects()
{
  
  fHistTriggerCounter = new TH1I("fHistTriggerCounter","Number of recored UPC triggers",13,0.5,13.5);
  TString gTriggerName[13] = {"CCUP4-B","CCUP2-B","CCUP7-B","CINT1-B","CTEST58-B", "CTEST59-B", "CTEST60-B", "CTEST61-B", "CCUP8-B", "CCUP9-B", "CCUP10-B", "CCUP11-B", "CCUP12-B"};
  for (Int_t i = 0; i<13; i++) fHistTriggerCounter->GetXaxis()->SetBinLabel(i+1,gTriggerName[i].Data());
  
  PostData(1, fHistTriggerCounter);

}//UserCreateOutputObjects


//_____________________________________________________________________________
void AliAnalysisTaskUpcTriggerCounter::UserExec(Option_t *) 
{

  
  AliAODEvent *aod = (AliAODEvent*) InputEvent();
  if(!aod) return;
	
  TString trigger = aod->GetFiredTriggerClasses();

  if(trigger.Contains("CCUP4-B")) 	fHistTriggerCounter->Fill(1);
  if(trigger.Contains("CCUP2-B")) 	fHistTriggerCounter->Fill(2);
  if(trigger.Contains("CCUP7-B")) 	fHistTriggerCounter->Fill(3);
  if(trigger.Contains("CINT1-B")) 	fHistTriggerCounter->Fill(4);
  if(trigger.Contains("CTEST58-B")) 	fHistTriggerCounter->Fill(5);
  if(trigger.Contains("CTEST59-B")) 	fHistTriggerCounter->Fill(6);
  if(trigger.Contains("CTEST60-B")) 	fHistTriggerCounter->Fill(7);
  if(trigger.Contains("CTEST61-B")) 	fHistTriggerCounter->Fill(8);
  if(trigger.Contains("CCUP8-B")) 	fHistTriggerCounter->Fill(9);
  if(trigger.Contains("CCUP9-B")) 	fHistTriggerCounter->Fill(10);
  if(trigger.Contains("CCUP10-B")) 	fHistTriggerCounter->Fill(11);
  if(trigger.Contains("CCUP11-B")) 	fHistTriggerCounter->Fill(12);
  if(trigger.Contains("CCUP12-B")) 	fHistTriggerCounter->Fill(13);
  
  PostData(1, fHistTriggerCounter);

}//UserExec


//_____________________________________________________________________________
void AliAnalysisTaskUpcTriggerCounter::Terminate(Option_t *) 
{

  cout<<"Analysis complete."<<endl;
}//Terminate
