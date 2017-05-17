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

#include "AliCDBManager.h"
#include "AliCDBEntry.h"
#include "AliGRPObject.h"
#include "AliTriggerInput.h"
#include "AliTriggerConfiguration.h"
#include "AliInputEventHandler.h"

#include "AliAnalysisMuonUtility.h"
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
  fhL0TriggerInputMLL(0),
  fhL0TriggerInputMUL(0),
  fhL0TriggerInputMSL(0),
  fhNMu(0)
{
  //
  // Default constructor
  //
  fNTrigClass=29;
  for(Int_t i=0;i<fNTrigClass;i++)
    fTriggerMask[i]=1<<i;
  fTriggerClasses[0]="-B-";
  fTriggerClasses[1]="CMSL7";
  fTriggerClasses[2]="CMSL8";
  fTriggerClasses[3]="CMSH7";
  fTriggerClasses[4]="CMUL7";
  fTriggerClasses[5]="CMLL7";
  fTriggerClasses[6]="CMSH8";
  fTriggerClasses[7]="CMUL8";
  fTriggerClasses[8]="CMLL8";
  fTriggerClasses[9]="CINT7";
  fTriggerClasses[10]="C0TVX";

  fTriggerClasses[11]="CINT7 (PS)";
  fTriggerClasses[12]="C0TVX (PS)";
  fTriggerClasses[13]="CMSL7 (PS CINT7)";
  fTriggerClasses[14]="CMSL8 (PS C0TVX)";
  fTriggerClasses[15]="CMSH7 (PS CINT7)";
  fTriggerClasses[16]="CMUL7 (PS CINT7)";
  fTriggerClasses[17]="CMLL7 (PS CINT7)";
  fTriggerClasses[18]="CMSH8 (PS C0TVX)";
  fTriggerClasses[19]="CMUL8 (PS C0TVX)";
  fTriggerClasses[20]="CMLL8 (PS C0TVX)";

  fTriggerClasses[21]="CMSL7 & 0MLL";
  fTriggerClasses[22]="CMSL7 & 0MUL";
  fTriggerClasses[23]="CMSL8 & 0MLL";
  fTriggerClasses[24]="CMSL8 & 0MUL";
  fTriggerClasses[25]="CINT7 & 0MUL";
  fTriggerClasses[26]="C0TVX & 0MUL";
  fTriggerClasses[27]="CINT7 & 0MSL";
  fTriggerClasses[28]="C0TVX & 0MSL";

  fL0TriggerInputMLL = 20; // reference to MLL L0 trigger
  fL0TriggerInputMUL = 21; // reference to MUL L0 trigger
  fL0TriggerInputMSL = 18; // reference to MSL L0 trigger
  minContributorsPileUp = 3;
}


//====================================================================================================================================================

AliAnalysisTaskLMREventFilter::AliAnalysisTaskLMREventFilter(const Char_t *name, AliMuonTrackCuts *cuts) : 
  AliAnalysisTaskSE(name), 
  fMuonTrackCuts(cuts),
  fEventTree(0x0),
  fOutputList(0x0),
  fAliLMREvent(0),  
  fhTriggers(0),
  fhL0TriggerInputMLL(0),
  fhL0TriggerInputMUL(0),
  fhL0TriggerInputMSL(0),
  fhNMu(0)
{
  // Constructor
  fNTrigClass=29;
  for(Int_t i=0;i<fNTrigClass;i++)
    fTriggerMask[i]=1<<i;
  fTriggerClasses[0]="-B-";
  fTriggerClasses[1]="CMSL7";
  fTriggerClasses[2]="CMSL8";
  fTriggerClasses[3]="CMSH7";
  fTriggerClasses[4]="CMUL7";
  fTriggerClasses[5]="CMLL7";
  fTriggerClasses[6]="CMSH8";
  fTriggerClasses[7]="CMUL8";
  fTriggerClasses[8]="CMLL8";
  fTriggerClasses[9]="CINT7";
  fTriggerClasses[10]="C0TVX";

  fTriggerClasses[11]="CINT7 (PS)";
  fTriggerClasses[12]="C0TVX (PS)";
  fTriggerClasses[13]="CMSL7 (PS CINT7)";
  fTriggerClasses[14]="CMSL8 (PS C0TVX)";
  fTriggerClasses[15]="CMSH7 (PS CINT7)";
  fTriggerClasses[16]="CMUL7 (PS CINT7)";
  fTriggerClasses[17]="CMLL7 (PS CINT7)";
  fTriggerClasses[18]="CMSH8 (PS C0TVX)";
  fTriggerClasses[19]="CMUL8 (PS C0TVX)";
  fTriggerClasses[20]="CMLL8 (PS C0TVX)";

  fTriggerClasses[21]="CMSL7 & 0MLL";
  fTriggerClasses[22]="CMSL7 & 0MUL";
  fTriggerClasses[23]="CMSL8 & 0MLL";
  fTriggerClasses[24]="CMSL8 & 0MUL";
  fTriggerClasses[25]="CINT7 & 0MUL";
  fTriggerClasses[26]="C0TVX & 0MUL";
  fTriggerClasses[27]="CINT7 & 0MSL";
  fTriggerClasses[28]="C0TVX & 0MSL";

  fL0TriggerInputMLL = 20; // reference to MLL L0 trigger
  fL0TriggerInputMUL = 21; // reference to MUL L0 trigger
  fL0TriggerInputMSL = 18; // reference to MSL L0 trigger
  minContributorsPileUp = 3;

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
  delete fhL0TriggerInputMLL;
  delete fhL0TriggerInputMUL;
  delete fhL0TriggerInputMSL;
  delete fhNMu;
  fhTriggers=NULL;
  fhNMu=NULL;
}

//====================================================================================================================================================

void AliAnalysisTaskLMREventFilter::UserCreateOutputObjects()
 {
  // Called once
  fMuonTrackCuts->SetFilterMask(AliMuonTrackCuts::kMuPdca); 
  fMuonTrackCuts->SetAllowDefaultParams(kTRUE);

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

  fhL0TriggerInputMLL = new TH1D("fhL0TriggerInputMLL","",120,-0.5,119.5);
  fOutputList->Add(fhL0TriggerInputMLL);
  fhL0TriggerInputMLL->Sumw2();

  fhL0TriggerInputMUL = new TH1D("fhL0TriggerInputMUL","",120,-0.5,119.5);
  fOutputList->Add(fhL0TriggerInputMUL);
  fhL0TriggerInputMUL->Sumw2();

  fhL0TriggerInputMSL = new TH1D("fhL0TriggerInputMSL","",120,-0.5,119.5);
  fOutputList->Add(fhL0TriggerInputMSL);
  fhL0TriggerInputMSL->Sumw2();

  for (Int_t i=0;i<fNTrigClass;i++)
    {
      fhTriggers->GetXaxis()->SetBinLabel(i+1,fTriggerClasses[i]);
      fhNMu->GetYaxis()->SetBinLabel(i+1,fTriggerClasses[i]);
    }

  PostData(1, fOutputList);
  PostData(2, fEventTree);
  printf("End of create Output\n");
}


void AliAnalysisTaskLMREventFilter::NotifyRun()
{
  fMuonTrackCuts->SetRun(fInputHandler);
  AliCDBManager *man = AliCDBManager::Instance();
  man->Init();
  man->SetDefaultStorage("raw://"); 
  man->SetRun(fInputHandler->GetEvent()->GetRunNumber()); 

  AliGRPObject* fGRPData = (AliGRPObject*) man->Get("GRP/GRP/Data")->GetObject();
  if(fGRPData->GetBeamType().Contains("p-p"))
    minContributorsPileUp = 3;
  if((fGRPData->GetBeamType().Contains("p-A")) || (fGRPData->GetBeamType().Contains("A-p")))
    minContributorsPileUp = 5;
  if(fGRPData->GetBeamType().Contains("A-A"))
    minContributorsPileUp = 5;
    
  AliTriggerConfiguration *cfg=(AliTriggerConfiguration*)man->Get("GRP/CTP/Config")->GetObject(); 
  TObjArray  inputs = cfg->GetInputs(); 
  for(Int_t i=0;i<inputs.GetEntriesFast();i++)
    {
      AliTriggerInput* inp =  (AliTriggerInput*) inputs[i];
      TString name=inp->GetName();
      if(name.Contains("0MUL"))
	fL0TriggerInputMUL = inp->GetIndexCTP(); // reference to MUL L0 trigger
      if(name.Contains("0MLL"))
	fL0TriggerInputMLL = inp->GetIndexCTP(); // reference to MLL L0 trigger
      if(name.Contains("0MSL"))
	fL0TriggerInputMSL = inp->GetIndexCTP(); // reference to MSL L0 trigger
    }
  fhL0TriggerInputMLL->Fill(fL0TriggerInputMLL);
  fhL0TriggerInputMUL->Fill(fL0TriggerInputMUL);
  fhL0TriggerInputMSL->Fill(fL0TriggerInputMSL);
}
//====================================================================================================================================================

void AliAnalysisTaskLMREventFilter::UserExec(Option_t *)
 {
  //   Main loop
  //   Called for each event
  UInt_t evtTrigSelect=0;
  UShort_t physicsSelectionMask=1<<0;
  AliAODEvent *fAOD = dynamic_cast<AliAODEvent *>(InputEvent());  
  if (!fAOD) 
    return;
  
  if (!IsSelectedTrigger(fAOD, kTRUE,evtTrigSelect,physicsSelectionMask)) 
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

  // All multiplicity are initialized at 166 and is used for error code of multselection non actived
  Float_t Multiplicity_V0M          = 166.;
  Float_t Multiplicity_ADM          = 166.;
  Float_t Multiplicity_SPDTracklets = 166.;
  Float_t Multiplicity_SPDClusters  = 166.;
  Float_t Multiplicity_RefMult05    = 166.;
  Float_t Multiplicity_RefMult08    = 166.;
  Float_t Multiplicity_V0A          = 166.;
  Float_t Multiplicity_V0C          = 166.;
  Float_t Multiplicity_V0EqA        = 166.;
  Float_t Multiplicity_V0EqC        = 166.;
  Float_t Multiplicity_V0EqM        = 166.;
  Float_t Multiplicity_ZNA          = 166.;
  Float_t Multiplicity_ZNC          = 166.;

  if (MultSelection) 
    {
      Multiplicity_V0M          = MultSelection->GetMultiplicityPercentile("V0M");
      Multiplicity_ADM          = MultSelection->GetMultiplicityPercentile("ADM");
      Multiplicity_SPDTracklets = MultSelection->GetMultiplicityPercentile("SPDTracklets");
      Multiplicity_SPDClusters  = MultSelection->GetMultiplicityPercentile("SPDClusters");
      Multiplicity_RefMult05    = MultSelection->GetMultiplicityPercentile("RefMult05");
      Multiplicity_RefMult08    = MultSelection->GetMultiplicityPercentile("RefMult08");
      Multiplicity_RefMult05    = MultSelection->GetMultiplicityPercentile("RefMult05");
      Multiplicity_V0A          = MultSelection->GetMultiplicityPercentile("V0A");
      Multiplicity_V0C          = MultSelection->GetMultiplicityPercentile("V0C");
      Multiplicity_V0EqA        = MultSelection->GetMultiplicityPercentile("V0EqA");
      Multiplicity_V0EqC        = MultSelection->GetMultiplicityPercentile("V0EqC");
      Multiplicity_V0EqM        = MultSelection->GetMultiplicityPercentile("V0EqM");
      Multiplicity_ZNA          = MultSelection->GetMultiplicityPercentile("ZNA");
      Multiplicity_ZNC          = MultSelection->GetMultiplicityPercentile("ZNC");
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
  fAliLMREvent->SetMultiplicity("V0A",Multiplicity_V0A);
  fAliLMREvent->SetMultiplicity("V0C",Multiplicity_V0C);
  fAliLMREvent->SetMultiplicity("V0EqA",Multiplicity_V0EqA);
  fAliLMREvent->SetMultiplicity("V0EqC",Multiplicity_V0EqC);
  fAliLMREvent->SetMultiplicity("V0EqM",Multiplicity_V0EqM);
  fAliLMREvent->SetMultiplicity("ZNA",Multiplicity_ZNA);
  fAliLMREvent->SetMultiplicity("ZNC",Multiplicity_ZNC);
  fAliLMREvent->SetTriggerString(triggerWord);
  fAliLMREvent->SetPhysicsSelectionMask(physicsSelectionMask);
  fAliLMREvent->SetIsPileupFromSPD(fAOD->IsPileupFromSPD(minContributorsPileUp));
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
	  trk->SetSelectionMask(fMuonTrackCuts->GetSelectionMask(track));
	  trk->SetLocalBoard((UShort_t)AliAnalysisMuonUtility::GetLoCircuit(track));
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

Bool_t AliAnalysisTaskLMREventFilter::IsSelectedTrigger(AliAODEvent *fAOD, Bool_t fillHisto,UInt_t &evtTrigSelect,UShort_t &physicsSelectionMask)
{
  Bool_t evtToBeProcessed = kFALSE;
  // -- Check Trigger 
  TString trigStr(((AliAODHeader*) fAOD->GetHeader())->GetFiredTriggerClasses());
  TObjArray * tokens = trigStr.Tokenize(" ");
  Int_t ntokens = tokens->GetEntriesFast();
  UInt_t mask=0,tmpMask=0;
  Bool_t goodTrig = kFALSE;
  for (Int_t itoken = 0; itoken < ntokens; ++itoken)
    {
      if ((((TObjString*)tokens->At(itoken))->String()).Contains(fTriggerClasses[0])) //for -B-
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

  
  for(Int_t i=1;i<=10;i++) // check constructor for fTriggerClasses[*]; 10 is the last basic data trigger name
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

  UInt_t fSelectedMaskMuonINT7=AliVEvent::kMuonSingleLowPt7|AliVEvent::kMuonSingleHighPt7|AliVEvent::kMuonLikeLowPt7|AliVEvent::kMuonUnlikeLowPt7;
  UInt_t fSelectedMaskMuonINT8=AliVEvent::kMuonSingleLowPt8|AliVEvent::kMuonSingleHighPt8|AliVEvent::kMuonLikeLowPt8|AliVEvent::kMuonUnlikeLowPt8;

  UInt_t fSelectMask=fInputHandler->IsEventSelected(); 
  Bool_t isMuonINT7selected  = fSelectMask&fSelectedMaskMuonINT7; 
  Bool_t isMuonC0TVXselected = fSelectMask&fSelectedMaskMuonINT8; 

  if(fSelectMask&AliVEvent::kINT7)
    physicsSelectionMask|=1<<1;
  if(fSelectMask&AliVEvent::kINT8)
    physicsSelectionMask|=1<<2;
  if(isMuonINT7selected)
    physicsSelectionMask|=1<<3;
  if(isMuonC0TVXselected)
    physicsSelectionMask|=1<<4;

  UInt_t inpmask = fAOD->GetHeader()->GetL0TriggerInputs();

  if(trigStr.Contains("CMSL7"))
    {
      if(isMuonINT7selected)
	{
	  fhTriggers->Fill(fTriggerClasses[13].Data(),1);
	  tmpMask|=fTriggerMask[13];
	  Int_t is0MLLfired = (inpmask & (1<<(fL0TriggerInputMLL-1)));
	  Int_t is0MULfired = (inpmask & (1<<(fL0TriggerInputMUL-1)));
	  if(is0MLLfired)
	    {
	      fhTriggers->Fill(fTriggerClasses[21].Data(),1);
	      tmpMask|=fTriggerMask[21];
	    }
	  if(is0MULfired)
	    {
	      fhTriggers->Fill(fTriggerClasses[22].Data(),1);
	      tmpMask|=fTriggerMask[22];
	    }
	}
    }
  if(trigStr.Contains("CMSL8"))
    {
      if(isMuonC0TVXselected)
	{
	  fhTriggers->Fill(fTriggerClasses[9].Data(),1);
	  tmpMask|=fTriggerMask[9];
	  Int_t is0MLLfired = (inpmask & (1<<(fL0TriggerInputMLL-1)));
	  Int_t is0MULfired = (inpmask & (1<<(fL0TriggerInputMUL-1)));
	  if(is0MLLfired)
	    {
	      fhTriggers->Fill(fTriggerClasses[23].Data(),1);
	      tmpMask|=fTriggerMask[23];
	    }
	  if(is0MULfired)
	    {
	      fhTriggers->Fill(fTriggerClasses[24].Data(),1);
	      tmpMask|=fTriggerMask[24];
	    }
	}
    }
  for(Int_t i=3;i<=5;i++)
    {
      if(trigStr.Contains(fTriggerClasses[i].Data()))
	{
	  if(isMuonINT7selected)
	    {
	      fhTriggers->Fill(fTriggerClasses[12+i].Data(),1);
	      tmpMask|=fTriggerMask[12+i]; 
	    }
	}
      if(trigStr.Contains(fTriggerClasses[i+3].Data()))
	if(isMuonC0TVXselected)
	  {
	    fhTriggers->Fill(fTriggerClasses[15+i].Data(),1);
	    tmpMask|=fTriggerMask[15+i]; 
	  }
    }
  if(trigStr.Contains("CINT7"))
    {
      if(fSelectMask&AliVEvent::kINT7)
	{
	  fhTriggers->Fill(fTriggerClasses[11].Data(),1);
	  tmpMask|=fTriggerMask[11];
	  Int_t is0MULfired = (inpmask & (1<<(fL0TriggerInputMUL-1)));
	  Int_t is0MSLfired = (inpmask & (1<<(fL0TriggerInputMSL-1)));
	  if(is0MULfired)
	    {
	      fhTriggers->Fill(fTriggerClasses[25].Data(),1);
	      tmpMask|=fTriggerMask[25];
	    }
	  if(is0MSLfired)
	    {
	      fhTriggers->Fill(fTriggerClasses[27].Data(),1);
	      tmpMask|=fTriggerMask[27];
	    }
	}
    }
  if(trigStr.Contains("C0TVX"))
    {
      if(fSelectMask&AliVEvent::kINT8)
	{
	  fhTriggers->Fill(fTriggerClasses[12].Data(),1);
	  tmpMask|=fTriggerMask[12];
	  Int_t is0MULfired = (inpmask & (1<<(fL0TriggerInputMUL-1)));
	  Int_t is0MSLfired = (inpmask & (1<<(fL0TriggerInputMSL-1)));
	  if(is0MULfired)
	    {
	      fhTriggers->Fill(fTriggerClasses[26].Data(),1);
	      tmpMask|=fTriggerMask[26];
	    }
	  if(is0MSLfired)
	    {
	      fhTriggers->Fill(fTriggerClasses[28].Data(),1);
	      tmpMask|=fTriggerMask[28];
	    }
	}
    }

  delete tokens;
  if(tmpMask!=0)
    mask=tmpMask;

  evtTrigSelect=mask;
  return evtToBeProcessed;
  
}
