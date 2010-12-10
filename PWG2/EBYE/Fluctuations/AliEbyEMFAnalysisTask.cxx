
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

/*-------------------------------------------------------------------------
 *                     AliEbyEMFAnalysisTask Class  
 *       This class deals with runing  Multiplicity Fluctuation Task 
 *                 origin: Satyajit Jena <sjena@cern.ch>
 *       MF: Multiplicity Fluctuation 
 *------------------------------------------------------------------------*/
 
#include "TChain.h"
#include "TString.h"
#include "TList.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3D.h"
#include "TCanvas.h"
#include "AliAnalysisManager.h"
#include "AliVEvent.h"
#include "AliESD.h"
#include "AliESDEvent.h"
#include "AliAODEvent.h"
#include "AliEbyEEventBase.h"
#include "AliEbyEMultiplicityFluctuationAnalysis.h"
#include "AliEbyEMFAnalysisTask.h"

ClassImp(AliEbyEMFAnalysisTask)

//________________________________________________________________________
AliEbyEMFAnalysisTask::AliEbyEMFAnalysisTask(const char *name) 
  : AliAnalysisTaskSE(name),
  fListPhy(0),
  fEbyEMFBase(0),  
  fEvtCounter(0),
  fECnt(0)
  
{
  // DefineInput(0, TChain::Class());
  DefineOutput(1, TList::Class());
}

//________________________________________________________________________
void AliEbyEMFAnalysisTask::UserCreateOutputObjects()
{
  fListPhy = new TList();

  fEvtCounter = new TH1F("hEvtCounter","Event Statistic",50,0,50);
  fListPhy->Add(fEvtCounter);
  fListPhy->Add(dynamic_cast<AliEbyEEventBase*>(fEbyEMFBase->GetEbyEEventBaseObject())->GetQA());
  fListPhy->Add(fEbyEMFBase->GetListMFQA());
  fListPhy->Add(fEbyEMFBase->GetListMeasureMF());

  PostData(1, fListPhy);
}

//________________________________________________________________________
void AliEbyEMFAnalysisTask::UserExec(Option_t *) 
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
  PostData(1, fListPhy); 

}

//______________________________________________________________________//
void AliEbyEMFAnalysisTask::Terminate(Option_t *) {

fListPhy = dynamic_cast<TList*> (GetOutputData(1));
  if (!fListPhy) {
    Error("Terminate","Out Put List not available");
    return;
  }

  Info("AliEbyEMFAnalysisTask"," AliEbyEMFAnalysisTask Successfully finished");
  
}
