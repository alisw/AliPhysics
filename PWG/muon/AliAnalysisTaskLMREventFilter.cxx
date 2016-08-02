/**************************************************************************
 * Copyright(c) 2004-2006, ALICE Experiment at CERN, All rights reserved. *
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

//========================================================================
//
//     Contact author: boris.teyssier@cern.ch | antonio.uras@cern.ch
//
//=========================================================================


#include "TChain.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "THnSparse.h"
#include "TClonesArray.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliESDVZERO.h"
#include "AliAODEvent.h"
#include "AliVVertex.h"
#include "AliVParticle.h"
#include "AliCentrality.h"
#include "AliLog.h"
#include "TBits.h"

#include "AliOADBContainer.h"
#include "AliOADBMultSelection.h"
#include "AliMultEstimator.h"
#include "AliMultVariable.h"
#include "AliMultInput.h"
#include "AliMultSelection.h"
#include "AliLMRMuon.h"
#include "AliLMREvent.h"
#include "AliAnalysisTaskLMREventFilter.h"

ClassImp(AliAnalysisTaskLMREventFilter)

//====================================================================================================================================================

AliAnalysisTaskLMREventFilter::AliAnalysisTaskLMREventFilter() : 
  AliAnalysisTaskSE(), 
  fMuonTrackCuts(0x0),
  fEventTree(0x0),
  fOutputList(0x0),
  fAliLMREvent(0),  
  fhTriggers(0),
  fhNMu(0)
{
  //
  // Default constructor
  //
  fNTrigClass=7;
  for(Int_t i=0;i<7;i++)
    fTriggerMask[i]=1<<i;
  fTriggerClasses[0]="-B-";
  fTriggerClasses[1]="CMSL7";
  fTriggerClasses[2]="CMSH7";
  fTriggerClasses[3]="CMUL7";
  fTriggerClasses[4]="CMLL7";
  fTriggerClasses[5]="CMSL7 & 0MLL";
  fTriggerClasses[6]="CMSL7 & 0MUL";

}


//====================================================================================================================================================

AliAnalysisTaskLMREventFilter::AliAnalysisTaskLMREventFilter(const Char_t *name, AliMuonTrackCuts *cuts) : 
  AliAnalysisTaskSE(name), 
  fMuonTrackCuts(cuts),
  fEventTree(0x0),
  fOutputList(0x0),
  fAliLMREvent(0),  
  fhTriggers(0),
  fhNMu(0)
{
  // Constructor
  fNTrigClass=7;
  for(Int_t i=0;i<7;i++)
    fTriggerMask[i]=1<<i;
  fTriggerClasses[0]="-B-";
  fTriggerClasses[1]="CMSL7";
  fTriggerClasses[2]="CMSH7";
  fTriggerClasses[3]="CMUL7";
  fTriggerClasses[4]="CMLL7";
  fTriggerClasses[5]="CMSL7 & 0MLL";
  fTriggerClasses[6]="CMSL7 & 0MUL";

  // Define input and output slots here
  DefineOutput(1, TList::Class());
  DefineOutput(2, TTree::Class());

}

//====================================================================================================================================================

AliAnalysisTaskLMREventFilter::~AliAnalysisTaskLMREventFilter() 
{
  delete fMuonTrackCuts; 
  fMuonTrackCuts=NULL;
  delete fEventTree;
  fEventTree=NULL;
  delete fOutputList;
  fOutputList=NULL;
  delete fAliLMREvent;
  fAliLMREvent=NULL;
  delete fhTriggers;
  delete fhNMu;
  fhTriggers=NULL;
  fhNMu=NULL;
}

//====================================================================================================================================================

void AliAnalysisTaskLMREventFilter::UserCreateOutputObjects()
 {
  // Called once

  fEventTree = new TTree("Data","Data");
  fAliLMREvent = new AliLMREvent();
  fEventTree->Branch("fAliLMREvent", &fAliLMREvent);
  
  fOutputList = new TList();
  fOutputList->SetOwner(kTRUE);
  
  fhTriggers = new TH1D("hTriggers","L2 Triggers",fNTrigClass,0,fNTrigClass);
  fOutputList->Add(fhTriggers);	  
  fhTriggers->Sumw2();

  fhNMu = new TH2D("hNMu","Number of Muon",50,0,50,fNTrigClass,0,fNTrigClass);
  fOutputList->Add(fhNMu);
  fhNMu->Sumw2();

  for (Int_t i=0;i<fNTrigClass;i++)
    {
      fhTriggers->GetXaxis()->SetBinLabel(i+1,fTriggerClasses[i]);
      fhNMu->GetYaxis()->SetBinLabel(i+1,fTriggerClasses[i]);
    }

  PostData(1, fOutputList);
  PostData(2, fEventTree);
  printf("End of create Output\n");
}

//====================================================================================================================================================

void AliAnalysisTaskLMREventFilter::UserExec(Option_t *)
 {
  //   Main loop
  //   Called for each event	
  UShort_t evtTrigSelect=0;
  AliAODEvent *fAOD = dynamic_cast<AliAODEvent *>(InputEvent());  
  if (!fAOD) 
    return;

  if (!IsSelectedTrigger(fAOD, kTRUE,evtTrigSelect)) 
    return; 
  
  if (!fAOD->GetNumberOfTracks())
    return;
    
  Int_t nmu = fAOD->GetNumberOfMuonTracks();
  
  if(evtTrigSelect&fTriggerMask[0])
    {
      fhNMu->Fill(nmu,fTriggerClasses[0],1);
      PostData(1, fOutputList); 
      return;
    }

  fhNMu->Fill(nmu,fTriggerClasses[0],1);
  for(Int_t i=1;i<fNTrigClass;i++)
    {
      if(evtTrigSelect&fTriggerMask[i])
	fhNMu->Fill(nmu,fTriggerClasses[i],1);
    }
  
  TString triggerWord(((AliAODHeader*) fAOD->GetHeader())->GetFiredTriggerClasses());
  // ---

  AliMultSelection *MultSelection = (AliMultSelection*)fAOD-> FindListObject("MultSelection");

  // all multiplicity are initialized at 166 and is used for error code of multselection non actived 
  Double_t Multiplicity_V0M          = 166.;
  Double_t Multiplicity_ADM          = 166.;
  Double_t Multiplicity_SPDTracklets = 166.;
  Double_t Multiplicity_SPDClusters  = 166.;
  Double_t Multiplicity_RefMult05    = 166.;
  Double_t Multiplicity_RefMult08    = 166.;

  if (MultSelection) 
    {
      Multiplicity_V0M          = MultSelection->GetMultiplicityPercentile("V0M");
      Multiplicity_ADM          = MultSelection->GetMultiplicityPercentile("ADM");
      Multiplicity_SPDTracklets = MultSelection->GetMultiplicityPercentile("SPDTracklets");
      Multiplicity_SPDClusters  = MultSelection->GetMultiplicityPercentile("SPDClusters");
      Multiplicity_RefMult05    = MultSelection->GetMultiplicityPercentile("RefMult05");
      Multiplicity_RefMult08    = MultSelection->GetMultiplicityPercentile("RefMult08");
    }
  
  
  AliAODVertex *vert = fAOD->GetPrimaryVertex();
  if (!vert) {
    printf ("No vertex found\n");
    return;
  }
  Double_t xvert  = vert->GetX();
  Double_t yvert  = vert->GetY();
  Double_t zvert  = vert->GetZ();
  Int_t vtxcontrib= vert->GetNContributors();
  Double_t evtPlane = fAOD->GetEventplane()->GetEventplane("V0",fAOD,2);
    
  Int_t runNumber = ((AliAODHeader*) fAOD->GetHeader())->GetRunNumber();
  fAliLMREvent->SetRunNumber(runNumber);
  fAliLMREvent->SetEventPlane(evtPlane);
  fAliLMREvent->SetXVertex(xvert);
  fAliLMREvent->SetYVertex(yvert);
  fAliLMREvent->SetZVertex(zvert);
  fAliLMREvent->SetVtxContributors(vtxcontrib);
  fAliLMREvent->SetMultiplicity("V0M",Multiplicity_V0M);
  fAliLMREvent->SetMultiplicity("ADM",Multiplicity_ADM);
  fAliLMREvent->SetMultiplicity("SPDTracklets",Multiplicity_SPDTracklets);
  fAliLMREvent->SetMultiplicity("SPDClusters",Multiplicity_SPDClusters);
  fAliLMREvent->SetMultiplicity("RefMult05",Multiplicity_RefMult05);
  fAliLMREvent->SetMultiplicity("RefMult08",Multiplicity_RefMult08);
  fAliLMREvent->SetTriggerString(triggerWord);

  if (nmu>0)
    {
      Int_t ntotTr = fAOD->GetNumberOfTracks(); 
      AliLMRMuon *trk = NULL;
      for (Int_t itr=0; itr<ntotTr; itr++) 
	{ 
	  AliAODTrack *track = dynamic_cast<AliAODTrack*>(fAOD->GetTrack(itr));
	  if(!track->IsMuonTrack()) continue;
    
	  Double_t p[3];
    
	  Int_t charge = track->Charge(); 
	  track->GetPxPyPz(p);
	  Double_t chi2Match = track->GetChi2MatchTrigger();
	  Int_t match = track->GetMatchTrigger();	
	  Double_t chi2 = track->Chi2perNDF();
	  Double_t rAbs = track->GetRAtAbsorberEnd();
    
	  TVector3 dcaAtVz  = fMuonTrackCuts->GetCorrectedDCA(track);
	  Double_t pTotMean = fMuonTrackCuts->GetAverageMomentum(track);
	  Double_t pDca = pTotMean * dcaAtVz.Mag();
    
	  // Create new Muon
	  trk = fAliLMREvent->AddMuon();
	  trk->SetMomentum(p[0],p[1],p[2]);
	  trk->SetCharge(charge);
	  trk->SetChi2Match(chi2Match);
	  trk->SetChi2(chi2);
	  trk->SetRabs(rAbs);
	  trk->SetpDCA(pDca);
	  trk->SetTriggerMatch(match);
	}
    }
  
 
 fEventTree->Fill();
 fAliLMREvent->Clear("");
  
 PostData(1, fOutputList); 
 PostData(2, fEventTree);

}

//====================================================================================================================================================

void AliAnalysisTaskLMREventFilter::Terminate(Option_t *) 
{
  // Draw result to the screen
  // Called once at the end of the query
  
  fOutputList = dynamic_cast<TList*> (GetOutputData(1));
  if (!fOutputList) {
    printf("ERROR: Output list not available\n");
    return;
  }  
}

//====================================================================================================================================================

Bool_t AliAnalysisTaskLMREventFilter::IsSelectedTrigger(AliAODEvent *fAOD, Bool_t fillHisto,UShort_t &evtTrigSelect)
{
  Bool_t evtToBeProcessed = kFALSE;
  // -- Check Trigger 
  TString trigStr(((AliAODHeader*) fAOD->GetHeader())->GetFiredTriggerClasses());
  TObjArray * tokens = trigStr.Tokenize(" ");
  Int_t ntokens = tokens->GetEntriesFast();
  UShort_t mask=0,tmpMask=0;
  Bool_t goodTrig = kFALSE;
  for (Int_t itoken = 0; itoken < ntokens; ++itoken)
    {
      if ((((TObjString*)tokens->At(itoken))->String()).Contains(fTriggerClasses[0])) 
	{
	  goodTrig = kTRUE;
	  mask=1<<0;
	  break;
	}
    }
  if(!goodTrig)
    return evtToBeProcessed;

  evtToBeProcessed = kTRUE;
  if (fillHisto) 
    fhTriggers->Fill(fTriggerClasses[0].Data(),1);
  
  for(Int_t i=1;i<fNTrigClass-2;i++)
    {
      goodTrig = kFALSE;
      for (Int_t itoken = 0; itoken < ntokens; ++itoken)
	{
	  if ((((TObjString*)tokens->At(itoken))->String()).Contains(fTriggerClasses[i])) 
	    {
	      goodTrig = kTRUE;
	      tmpMask|=fTriggerMask[i];
	      break;
	    }
	}
      if(goodTrig)
	if (fillHisto) 
	  fhTriggers->Fill(fTriggerClasses[i].Data(),1);
    }
  
  if(trigStr.Contains("CMSL7"))
    {
      Int_t inputMLL = 20; // reference to MLL L0 trigger
      Int_t inputMUL = 21; // reference to MUL L0 trigger
      UInt_t inpmask = fAOD->GetHeader()->GetL0TriggerInputs();
      Int_t is0MLLfired = (inpmask & (1<<(inputMLL-1)));
      Int_t is0MULfired = (inpmask & (1<<(inputMUL-1)));
      
      if(is0MLLfired)
	{
	  fhTriggers->Fill(fTriggerClasses[fNTrigClass-2].Data(),1);
	  tmpMask|=fTriggerMask[fNTrigClass-2];
	}
      if(is0MULfired)
	{
	  fhTriggers->Fill(fTriggerClasses[fNTrigClass-1].Data(),1);
	  tmpMask|=fTriggerMask[fNTrigClass-1];
	}
    }
  delete tokens;
  if(tmpMask!=0)
    mask=tmpMask;

  evtTrigSelect=mask;
  return evtToBeProcessed;
  
}
