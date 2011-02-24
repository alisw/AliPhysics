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

/* $Id$ */

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

#include "AliAnalysisManager.h"
#include "AliAODHandler.h"
#include "AliAODEvent.h"
#include "AliAODVertex.h"
#include "AliESDtrack.h"
#include "AliAODTrack.h"
#include "AliAODMCHeader.h"
#include "AliAODMCParticle.h"
#include "AliAODRecoDecayHF2Prong.h"
#include "AliAODRecoDecayHF3Prong.h"
#include "AliAODRecoDecayHF4Prong.h"
#include "AliAODRecoCascadeHF.h"
#include "AliAnalysisVertexingHF.h"
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisTaskSECompareHF.h"

ClassImp(AliAnalysisTaskSECompareHF)


//________________________________________________________________________
AliAnalysisTaskSECompareHF::AliAnalysisTaskSECompareHF():
AliAnalysisTaskSE(),
fOutput(0), 
fNtupleCmp(0),
fHistMass(0),
fHistNEvents(0),
fVHF(0)
{
  // Default constructor

  // NO DefineOutput() HERE (ONLY IN STANDARD CONSTRUCTOR)
}

//________________________________________________________________________
AliAnalysisTaskSECompareHF::AliAnalysisTaskSECompareHF(const char *name):
AliAnalysisTaskSE(name),
fOutput(0), 
fNtupleCmp(0),
fHistMass(0),
fHistNEvents(0),
fVHF(0)
{
  // Standard constructor

  // Output slot #1 writes into a TList container
  DefineOutput(1,TList::Class());  //My private output
  // Output slot #2 writes into a TNtuple container
  DefineOutput(2,TNtuple::Class());  //My private output
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

  fHistNEvents = new TH1F("fHistNEvents", "Number of processed events; ; Events",3,-1.5,1.5);
  fHistNEvents->Sumw2();
  fHistNEvents->SetMinimum(0);
  fOutput->Add(fHistNEvents);

  OpenFile(2); // 2 is the slot number of the ntuple
  fNtupleCmp = new TNtuple("fNtupleCmp","Charm comparison","pdg:nprongs:VxRec:VxTrue:ErrVx:VyRec:VyTrue:ErrVy:VzRec:VzTrue:ErrVz:Chi2toNDF:PtRec:Mrec:CPta:Prodd0");

  return;
}

//________________________________________________________________________
void AliAnalysisTaskSECompareHF::UserExec(Option_t */*option*/)
{
  // Execute analysis for current event:
  // heavy flavor candidates association to MC truth

  
  AliAODEvent *aod = dynamic_cast<AliAODEvent*> (InputEvent());

  TClonesArray *inputArrayVertices = 0;
  TClonesArray *inputArrayD0toKpi = 0;
  TClonesArray *inputArrayDstar = 0;

  if(!aod && AODEvent() && IsStandardAOD()) {
    // In case there is an AOD handler writing a standard AOD, use the AOD 
    // event in memory rather than the input (ESD) event.    
    aod = dynamic_cast<AliAODEvent*> (AODEvent());
    // in this case the braches in the deltaAOD (AliAOD.VertexingHF.root)
    // have to taken from the AOD event hold by the AliAODExtension
    AliAODHandler* aodHandler = (AliAODHandler*) 
      ((AliAnalysisManager::GetAnalysisManager())->GetOutputEventHandler());
    if(aodHandler->GetExtensions()) {
      AliAODExtension *ext = (AliAODExtension*)aodHandler->GetExtensions()->FindObject("AliAOD.VertexingHF.root");
      AliAODEvent *aodFromExt = ext->GetAOD();
      // load HF vertices                
      inputArrayVertices = (TClonesArray*)aodFromExt->GetList()->FindObject("VerticesHF");
      // load D0->Kpi candidates
      inputArrayD0toKpi = (TClonesArray*)aodFromExt->GetList()->FindObject("D0toKpi");
      // load D*+ candidates                                                   
      inputArrayDstar = (TClonesArray*)aodFromExt->GetList()->FindObject("Dstar");
    }
  } else if(aod) {
    // load HF vertices                
    inputArrayVertices = (TClonesArray*)aod->GetList()->FindObject("VerticesHF");
    // load D0->Kpi candidates                                                 
    inputArrayD0toKpi = (TClonesArray*)aod->GetList()->FindObject("D0toKpi");
    // load D*+ candidates                                                   
    inputArrayDstar = (TClonesArray*)aod->GetList()->FindObject("Dstar");
  }


  if(!inputArrayVertices || !aod) {
    printf("AliAnalysisTaskSECompareHF::UserExec: Vertices branch not found!\n");
    return;
  }
  if(!inputArrayD0toKpi) {
    printf("AliAnalysisTaskSECompareHF::UserExec: D0toKpi branch not found!\n");
    return;
  }
  if(!inputArrayDstar) {
    printf("AliAnalysisTaskSECompareHF::UserExec: Dstar branch not found!\n");
    return;
  }
  

  fHistNEvents->Fill(0); // count event
  // Post the data already here
  PostData(1,fOutput);

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
  
  Int_t nprongs,lab,okD0,okD0bar,pdg;
  Bool_t unsetvtx;    
  Double_t invmass,posRec[3],posTrue[3],covRec[6],errx,erry,errz;
  AliAODRecoDecayHF2Prong *d2=0;
  AliAODRecoDecayHF3Prong *d3=0;
  AliAODRecoDecayHF4Prong *d4=0;

  Int_t pdgDgD0toKpi[2]={321,211};
  Int_t pdgDgDplustoKpipi[3]={321,211,211};
  Int_t pdgDgD0toKpipipi[4]={321,211,211,211};

  // loop over vertices
  Int_t nVertices = inputArrayVertices->GetEntriesFast();
  if(fDebug>1) printf("Number of vertices: %d\n",nVertices);

  for (Int_t iVtx = 0; iVtx < nVertices; iVtx++) {
    AliAODVertex *vtx = (AliAODVertex*)inputArrayVertices->UncheckedAt(iVtx);

    vtx->GetXYZ(posRec);  
    vtx->GetCovarianceMatrix(covRec);
    errx=1.; erry=1.; errz=1.;
    if(covRec[0]>0) errx = TMath::Sqrt(covRec[0]);
    if(covRec[2]>0) erry = TMath::Sqrt(covRec[2]);
    if(covRec[5]>0) errz = TMath::Sqrt(covRec[5]);
	  

    // get parent AliAODRecoDecayHF
    nprongs= vtx->GetNDaughters();
    // check that the daughters have kITSrefit and kTPCrefit
    Bool_t allDgOK=kTRUE;
    for(Int_t idg=0; idg<nprongs; idg++) {
      AliAODTrack *track = (AliAODTrack*)vtx->GetDaughter(idg);
      if(!(track->GetStatus()&AliESDtrack::kITSrefit)) allDgOK=kFALSE;
      if(!(track->GetStatus()&AliESDtrack::kTPCrefit)) allDgOK=kFALSE;
    }
    if(!allDgOK) continue;


    switch(nprongs) {
    case 2: // look for D0->Kpi
      d2 = (AliAODRecoDecayHF2Prong*)vtx->GetParent();
      if(d2->IsLikeSign()) continue;
      if(d2->Charge() != 0) continue; // these are D* 
      lab = d2->MatchToMC(421,mcArray,2,pdgDgD0toKpi);
      if(lab>=0) {
	unsetvtx=kFALSE;
	if(!d2->GetOwnPrimaryVtx()) {
	  d2->SetOwnPrimaryVtx(vtx1); // needed to compute all variables
	  unsetvtx=kTRUE;
	}
	okD0=1; okD0bar=1; 
	AliAODMCParticle *dMC = (AliAODMCParticle*)mcArray->At(lab);
	pdg = dMC->GetPdgCode();
	invmass = (pdg==421 ? d2->InvMassD0() : d2->InvMassD0bar());
	// get a daughter for true pos of decay vertex
	AliAODMCParticle *dg0MC = (AliAODMCParticle*)mcArray->At(dMC->GetDaughter(0));
	dg0MC->XvYvZv(posTrue);
	fHistMass->Fill(invmass);
	// Post the data already here
	PostData(1,fOutput);
	Float_t tmp[16]={(Float_t)pdg,(Float_t)nprongs,
			 (Float_t)posRec[0],(Float_t)posTrue[0],(Float_t)errx,
			 (Float_t)posRec[1],(Float_t)posTrue[1],(Float_t)erry,
			 (Float_t)posRec[2],(Float_t)posTrue[2],(Float_t)errz,
			 (Float_t)d2->GetReducedChi2(),(Float_t)d2->Pt(),(Float_t)invmass,
			 (Float_t)d2->CosPointingAngle(),(Float_t)d2->Prodd0d0()};
	fNtupleCmp->Fill(tmp);
	PostData(2,fNtupleCmp);
      
	if(unsetvtx) d2->UnsetOwnPrimaryVtx();
      }
      break;
    case 3: // look for D+
      d3 = (AliAODRecoDecayHF3Prong*)vtx->GetParent();
      if(d3->IsLikeSign()) continue;
      lab = d3->MatchToMC(411,mcArray,3,pdgDgDplustoKpipi);
      if(lab>=0) {
	unsetvtx=kFALSE;
	if(!d3->GetOwnPrimaryVtx()) {
	  d3->SetOwnPrimaryVtx(vtx1); // needed to compute all variables
	  unsetvtx=kTRUE;
	}
	//if(d3->SelectDplus(fVHF->GetDplusCuts())) {
	AliAODMCParticle *dMC = (AliAODMCParticle*)mcArray->At(lab);
	pdg = dMC->GetPdgCode();
	invmass = d3->InvMassDplus();
	// get a daughter for true pos of decay vertex
	AliAODMCParticle *dg0MC = (AliAODMCParticle*)mcArray->At(dMC->GetDaughter(0));
	dg0MC->XvYvZv(posTrue);
	Float_t tmp[16]={(Float_t)pdg,(Float_t)nprongs,
			 (Float_t)posRec[0],(Float_t)posTrue[0],(Float_t)errx,
			 (Float_t)posRec[1],(Float_t)posTrue[1],(Float_t)erry,
			 (Float_t)posRec[2],(Float_t)posTrue[2],(Float_t)errz,
			 (Float_t)d3->GetReducedChi2(),(Float_t)d3->Pt(),(Float_t)invmass,
			 (Float_t)d3->CosPointingAngle(),(Float_t)(d3->Getd0Prong(0)*d3->Getd0Prong(1)*d3->Getd0Prong(2))};
	fNtupleCmp->Fill(tmp);
	PostData(2,fNtupleCmp);
	//}
	if(unsetvtx) d3->UnsetOwnPrimaryVtx();
      }
      break;
    case 4: // look for D0->Kpipipi
      d4 = (AliAODRecoDecayHF4Prong*)vtx->GetParent();
      if(d4->IsLikeSign()) continue;
      lab = d4->MatchToMC(421,mcArray,4,pdgDgD0toKpipipi);
      if(lab>=0) {
	unsetvtx=kFALSE;
	if(!d4->GetOwnPrimaryVtx()) {
	  d4->SetOwnPrimaryVtx(vtx1); // needed to compute all variables
	  unsetvtx=kTRUE;
	}
	okD0=0; okD0bar=0; 
	//if(d4->SelectD0(fVHF->GetD0to4ProngsCuts(),okD0,okD0bar)) {
	AliAODMCParticle *dMC = (AliAODMCParticle*)mcArray->At(lab);
	pdg = dMC->GetPdgCode();
	//invmass = (pdg==421 ? d->InvMassD0() : d->InvMassD0bar());
	invmass = 10.; 	  // dummy
	// get a daughter for true pos of decay vertex
	AliAODMCParticle *dg0MC = (AliAODMCParticle*)mcArray->At(dMC->GetDaughter(0));
	dg0MC->XvYvZv(posTrue);
	Float_t tmp[16]={(Float_t)pdg,(Float_t)nprongs,
			 (Float_t)posRec[0],(Float_t)posTrue[0],(Float_t)errx,
			 (Float_t)posRec[1],(Float_t)posTrue[1],(Float_t)erry,
			 (Float_t)posRec[2],(Float_t)posTrue[2],(Float_t)errz,
			 (Float_t)d4->GetReducedChi2(),(Float_t)d4->Pt(),(Float_t)invmass,
			 (Float_t)d4->CosPointingAngle(),(Float_t)(d4->Getd0Prong(0)*d4->Getd0Prong(1)*d4->Getd0Prong(2)*d4->Getd0Prong(3))};
	fNtupleCmp->Fill(tmp);
	PostData(2,fNtupleCmp);
	//}
	if(unsetvtx) d4->UnsetOwnPrimaryVtx();
      }
      break;
    }

  } // end loop on vertices


  // loop over D*+ candidates
  /*
  for (Int_t iDstar = 0; iDstar < inputArrayDstar->GetEntries(); iDstar++) {
    AliAODRecoCascadeHF *c = (AliAODRecoCascadeHF*)inputArrayDstar->UncheckedAt(iDstar);
    Int_t labDstar = c->MatchToMC(413,421,mcArray);
    if(labDstar>=0) printf("GOOD MATCH FOR D*+\n"); 
  }
  */

  // Post the data already here
  PostData(1,fOutput);
  PostData(2,fNtupleCmp);

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

  fHistMass = dynamic_cast<TH1F*>(fOutput->FindObject("fHistMass"));
  fHistNEvents = dynamic_cast<TH1F*>(fOutput->FindObject("fHistNEvents"));

  //fNtupleCmp = dynamic_cast<TNtuple*> (GetOutputData(2));

  return;
}

