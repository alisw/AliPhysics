/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: Satyajit Jena.                                                 *
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


#include "Riostream.h"
#include "TChain.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TList.h"
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisManager.h"
#include "AliStack.h"
#include "AliESDtrackCuts.h"
#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliAODEvent.h"
#include "AliMCEvent.h"
#include "AliMultiplicity.h"
#include "AliEbyEFluctuationAnalysis.h"
#include "AliEbyEEventBase.h"

#include "AliEbyEFluctuationTask.h"

ClassImp(AliEbyEFluctuationTask)

//________________________________________________________________________
AliEbyEFluctuationTask::AliEbyEFluctuationTask() 
:AliAnalysisTaskSE(),
  fListPhy(0),
  fEbyEMFBase(0),  
  fEvtCounter(0),
  fECnt(0)
{
  // Dummy constructor ALWAYS needed for I/O.
}

//________________________________________________________________________
AliEbyEFluctuationTask::AliEbyEFluctuationTask(const char *name) 
  :AliAnalysisTaskSE(name),
   fListPhy(0),
   fEbyEMFBase(0),  
   fEvtCounter(0),
   fECnt(0)
{
  
  DefineOutput(1, TList::Class());                                         
}

//________________________________________________________________________
AliEbyEFluctuationTask::~AliEbyEFluctuationTask()
{
   
    if (fListPhy && !AliAnalysisManager::GetAnalysisManager()->IsProofMode()) {
        delete fListPhy;
    }
}

//________________________________________________________________________
void AliEbyEFluctuationTask::UserCreateOutputObjects()
{
  
    fListPhy = new TList();
    fListPhy->SetOwner();  
    fEvtCounter = new TH1F("hEvtCounter","Event Statistic",50,0,50);
    fListPhy->Add(fEvtCounter);
    
    fListPhy->Add(dynamic_cast<AliEbyEEventBase*>(fEbyEMFBase->GetEbyEEventBaseObject())->GetQA());
    fListPhy->Add(fEbyEMFBase->GetListMeasure());
    
    PostData(1, fListPhy); 
}

//________________________________________________________________________
void AliEbyEFluctuationTask::UserExec(Option_t *) 
{
   
  fEvtCounter->Fill(0);  
  fECnt++;
  
  TString gAnalType = dynamic_cast<AliEbyEEventBase*>(fEbyEMFBase->GetEbyEEventBaseObject())->GetAnalysisLevel();
 
  if(gAnalType == "ESD") 
    {
      AliESDEvent* gESD = dynamic_cast<AliESDEvent*>(InputEvent()); 
      if (!gESD) {
	Printf("ERROR: gESD not available");
	return;    
      }
      const AliESDVertex *vertex = dynamic_cast<AliEbyEEventBase*>(fEbyEMFBase->GetEbyEEventBaseObject())->GetVertex(gESD,dynamic_cast<AliEbyEEventBase*>(fEbyEMFBase->GetEbyEEventBaseObject())->GetAnalysisMode(),dynamic_cast<AliEbyEEventBase*>(fEbyEMFBase->GetEbyEEventBaseObject())->GetVxMax(),dynamic_cast<AliEbyEEventBase*>(fEbyEMFBase->GetEbyEEventBaseObject())->GetVyMax(),dynamic_cast<AliEbyEEventBase*>(fEbyEMFBase->GetEbyEEventBaseObject())->GetVzMax());
      if(vertex){

      
	fEbyEMFBase->Analyze(gESD);
      }

    }

}


//________________________________________________________________________
void AliEbyEFluctuationTask::Terminate(Option_t *) 
{
    // Draw result to screen, or perform fitting, normalizations
    // Called once at the end of the query
        
  fListPhy = dynamic_cast<TList*> (GetOutputData(1));
  if(!fListPhy) { Printf("ERROR: could not retrieve TList fOutput"); return; }
        

}
