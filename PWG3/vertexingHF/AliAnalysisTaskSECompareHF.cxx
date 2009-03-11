/**************************************************************************
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

/////////////////////////////////////////////////////////////
//
// AliAnalysisTaskSE for the comparison of heavy flavor
// decay candidates with the MC truth.
//
// Author: A.Dainese, andrea.dainese@lnl.infn.it
/////////////////////////////////////////////////////////////

#include <TClonesArray.h>
#include <TNtuple.h>
#include <TList.h>
#include <TH1F.h>

#include "AliAODEvent.h"
#include "AliAODVertex.h"
#include "AliAODTrack.h"
#include "AliAODMCHeader.h"
#include "AliAODMCParticle.h"
#include "AliAODRecoDecayHF2Prong.h"
#include "AliAODRecoCascadeHF.h"
#include "AliAnalysisVertexingHF.h"
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisTaskSECompareHF.h"

ClassImp(AliAnalysisTaskSECompareHF)


//________________________________________________________________________
AliAnalysisTaskSECompareHF::AliAnalysisTaskSECompareHF():
AliAnalysisTaskSE(),
fOutput(0), 
fNtupleD0Cmp(0),
fHistMass(0),
fVHF(0)
{
  // Default constructor

  // Output slot #1 writes into a TList container
  DefineOutput(1,TList::Class());  //My private output
}

//________________________________________________________________________
AliAnalysisTaskSECompareHF::AliAnalysisTaskSECompareHF(const char *name):
AliAnalysisTaskSE(name),
fOutput(0), 
fNtupleD0Cmp(0),
fHistMass(0),
fVHF(0)
{
  // Default constructor

  // Output slot #1 writes into a TList container
  DefineOutput(1,TList::Class());  //My private output
}

//________________________________________________________________________
AliAnalysisTaskSECompareHF::~AliAnalysisTaskSECompareHF()
{
  // Destructor
  if (fOutput) {
    delete fOutput;
    fOutput = 0;
  }
  if (fVHF) {
    delete fVHF;
    fVHF = 0;
  }
}  

//________________________________________________________________________
void AliAnalysisTaskSECompareHF::Init()
{
  // Initialization

  if(fDebug > 1) printf("AnalysisTaskSECompareHF::Init() \n");

  gROOT->LoadMacro("ConfigVertexingHF.C");

  fVHF = (AliAnalysisVertexingHF*)gROOT->ProcessLine("ConfigVertexingHF()");  
  fVHF->PrintStatus();

  return;
}

//________________________________________________________________________
void AliAnalysisTaskSECompareHF::UserCreateOutputObjects()
{
  // Create the output container
  //
  if(fDebug > 1) printf("AnalysisTaskSECompareHF::UserCreateOutputObjects() \n");

  // Several histograms are more conveniently managed in a TList
  fOutput = new TList();
  fOutput->SetOwner();

  fHistMass = new TH1F("fHistMass", "D^{0} invariant mass; M [GeV]; Entries",200,1.765,1.965);
  fHistMass->Sumw2();
  fHistMass->SetMinimum(0);
  fOutput->Add(fHistMass);

  fNtupleD0Cmp = new TNtuple("fNtupleD0Cmp","D0 comparison","pdg:VxRec:VxTrue:PtRec:PtTrue");
  fOutput->Add(fNtupleD0Cmp);

  return;
}

//________________________________________________________________________
void AliAnalysisTaskSECompareHF::UserExec(Option_t */*option*/)
{
  // Execute analysis for current event:
  // heavy flavor candidates association to MC truth
  
  AliAODEvent *aod = dynamic_cast<AliAODEvent*> (InputEvent());
 
  // load D0->Kpi candidates                                                   
  TClonesArray *inputArrayD0toKpi =
    (TClonesArray*)aod->GetList()->FindObject("D0toKpi");
  if(!inputArrayD0toKpi) {
    printf("AliAnalysisTaskSECompareHF::UserExec: D0toKpi branch not found!\n");
    return;
  }

  /*
  // load D*+ candidates                                                   
  TClonesArray *inputArrayDstar =
    (TClonesArray*)aod->GetList()->FindObject("Dstar");
  if(!inputArrayDstar) {
    printf("AliAnalysisTaskSECompareHF::UserExec: Dstar branch not found!\n");
    return;
  }
  */

  // AOD primary vertex
  AliAODVertex *vtx1 = (AliAODVertex*)aod->GetPrimaryVertex();
  //vtx1->Print();

  // load MC particles
  TClonesArray *mcArray = 
    (TClonesArray*)aod->GetList()->FindObject(AliAODMCParticle::StdBranchName());
  if(!mcArray) {
    printf("AliAnalysisTaskSECompareHF::UserExec: MC particles branch not found!\n");
    return;
  }

  // load MC header
  AliAODMCHeader *mcHeader = 
    (AliAODMCHeader*)aod->GetList()->FindObject(AliAODMCHeader::StdBranchName());
  if(!mcHeader) {
    printf("AliAnalysisTaskSECompareHF::UserExec: MC header branch not found!\n");
    return;
  }
  
    
  // loop over D0->Kpi candidates
  Int_t nInD0toKpi = inputArrayD0toKpi->GetEntriesFast();
  printf("Number of D0->Kpi: %d\n",nInD0toKpi);

  for (Int_t iD0toKpi = 0; iD0toKpi < nInD0toKpi; iD0toKpi++) {
    AliAODRecoDecayHF2Prong *d = (AliAODRecoDecayHF2Prong*)inputArrayD0toKpi->UncheckedAt(iD0toKpi);
    Bool_t unsetvtx=kFALSE;
    if(!d->GetOwnPrimaryVtx()) {
      d->SetOwnPrimaryVtx(vtx1); // needed to compute all variables
      unsetvtx=kTRUE;
    }
    Int_t okD0=0,okD0bar=0; 
    if(d->SelectD0(fVHF->GetD0toKpiCuts(),okD0,okD0bar)) {
      Int_t labD0 = d->MatchToMC(421,mcArray);
      if(labD0>=0) {
	AliAODMCParticle *partD0 = (AliAODMCParticle*)mcArray->At(labD0);
	Int_t pdgD0 = partD0->GetPdgCode();
	Double_t invmass = (pdgD0==421 ? d->InvMassD0() : d->InvMassD0bar());
	fHistMass->Fill(invmass);
	// Post the data already here
	PostData(1,fOutput);

	fNtupleD0Cmp->Fill(pdgD0,d->Xv(),partD0->Xv(),d->Pt(),partD0->Pt());
      }
    }
    if(unsetvtx) d->UnsetOwnPrimaryVtx();
  } // end loop on D0->Kpi


  // SIMPLE EXAMPLE OF D*+ MATCHING
  /*
  // loop over D*+ candidates
  for (Int_t iDstar = 0; iDstar < inputArrayDstar->GetEntries(); iDstar++) {
    AliAODRecoCascadeHF *c = (AliAODRecoCascadeHF*)inputArrayDstar->UncheckedAt(iDstar);
    Int_t labDstar = c->MatchToMC(413,421,mcArray);
    if(labDstar>=0) printf("GOOD MATCH FOR D*+\n"); 
  }
  */

  return;
}

//________________________________________________________________________
void AliAnalysisTaskSECompareHF::Terminate(Option_t */*option*/)
{
  // Terminate analysis
  //
  if(fDebug > 1) printf("AnalysisTaskSECompareHF: Terminate() \n");

  fOutput = dynamic_cast<TList*> (GetOutputData(1));
  if (!fOutput) {     
    printf("ERROR: fOutput not available\n");
    return;
  }

  fNtupleD0Cmp = dynamic_cast<TNtuple*>(fOutput->FindObject("fNtupleD0Cmp"));
  fHistMass = dynamic_cast<TH1F*>(fOutput->FindObject("fHistMass"));

  return;
}

