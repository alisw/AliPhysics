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
 
#include "AliAnalysisTaskDiJets.h"
#include "AliAODEvent.h"
#include "AliAODJet.h"
#include "AliAODDiJet.h"
#include <TClonesArray.h>
#include <TLorentzVector.h>


ClassImp(AliAnalysisTaskDiJets)

////////////////////////////////////////////////////////////////////////

AliAnalysisTaskDiJets::AliAnalysisTaskDiJets():
    AliAnalysisTaskSE(),
    fDiJets(0),
    fDiJetsIn(0)
{
  // Default constructor
}

AliAnalysisTaskDiJets::AliAnalysisTaskDiJets(const char* name):
    AliAnalysisTaskSE(name),
    fDiJets(0),
    fDiJetsIn(0)
{
  // Default constructor
}

void AliAnalysisTaskDiJets::UserCreateOutputObjects()
{
// Create the output container
//
    if (fDebug > 1) printf("AnalysisTaskDiJets::CreateOutPutData() \n");
    fDiJets = new TClonesArray("AliAODDiJet", 0);
    fDiJets->SetName("Dijets");
    AddAODBranch("TClonesArray", &fDiJets);
}

void AliAnalysisTaskDiJets::Init()
{
    // Initialization
    if (fDebug > 1) printf("AnalysisTaskDiJets::Init() \n");
}

void AliAnalysisTaskDiJets::UserExec(Option_t */*option*/)
{
// Execute analysis for current event
//
    fDiJets->Delete();
    AliAODEvent* aod   = dynamic_cast<AliAODEvent*> (InputEvent());

    if(!aod){
      // We do not have an input AOD, look in the output
      aod = AODEvent();
      if(!aod){
	Printf("%s:%d AODEvent not found in the Output",(char*)__FILE__,__LINE__);
	return;
      }    
    }

    TClonesArray* jets = aod->GetJets();

    // N.B. if we take the aod from the output this is always
    // empty and since it is the same as fDiJets 
    fDiJetsIn = (TClonesArray*) (aod->GetList()->FindObject("Dijets"));
   
    if (fDiJetsIn) {
	printf("Found %d dijets in old list \n", fDiJetsIn->GetEntries());
	AliAODJet* jj1, *jj2;
	AliAODDiJet* testJ;
	
	if (fDiJetsIn->GetEntries() > 0) {
	    testJ = (AliAODDiJet*) (fDiJetsIn->At(0));
	    jj1 = testJ->Jet(0);
	    jj1->Print("");
	    jj2 = testJ->Jet(1);
	    jj2->Print("");
	}
    }
    
    
    Int_t nj = jets->GetEntriesFast();
    printf("There are %5d jets in the event \n", nj);
    if (nj > 1) {
		for (Int_t iJet1=0; iJet1<nj; iJet1++){
			AliAODJet* jet1 = (AliAODJet*) (jets->At(iJet1));
			TLorentzVector v1 = *(jet1->MomentumVector());
			for (Int_t iJet2=iJet1+1; iJet2<nj; iJet2++){
				AliAODJet* jet2 = (AliAODJet*) (jets->At(iJet2));
				TLorentzVector v2 = *(jet2->MomentumVector());
				TLorentzVector v = v1 + v2;
				Int_t ndi = fDiJets->GetEntriesFast();
				TClonesArray &lref = *fDiJets;
				new(lref[ndi]) AliAODDiJet(v);;
				AliAODDiJet* dijet = (AliAODDiJet*) (fDiJets->At(ndi));
				dijet->SetJetRefs(jet1, jet2);
			}
		}
    }
    
}

void AliAnalysisTaskDiJets::Terminate(Option_t */*option*/)
{
// Terminate analysis
//
    if (fDebug > 1) printf("AnalysisDiJets: Terminate() \n");
}

