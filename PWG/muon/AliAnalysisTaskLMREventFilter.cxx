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
  fhBeamType(0),
  fhL0TriggerInputMLL(0),
  fhL0TriggerInputMUL(0),
  fhL0TriggerInputMSL(0),
  fhL0TriggerInputTVX(0),
  fhL2TriggerCINT7CENTNOTRD(0),
  fhNMu(0),
  fminContributorsPileUp(0)
{
  //
  // Default constructor
  //

  SetTriggerClasses();
 
}

//====================================================================================================================================================

AliAnalysisTaskLMREventFilter::AliAnalysisTaskLMREventFilter(const Char_t *name, AliMuonTrackCuts *cuts) :
  AliAnalysisTaskSE(name),
  fMuonTrackCuts(cuts),
  fEventTree(0x0),
  fOutputList(0x0),
  fAliLMREvent(0), 
  fhTriggers(0),
  fhBeamType(0),
  fhL0TriggerInputMLL(0),
  fhL0TriggerInputMUL(0),
  fhL0TriggerInputMSL(0),
  fhL0TriggerInputTVX(0),
  fhL2TriggerCINT7CENTNOTRD(0),
  fhNMu(0),
  fminContributorsPileUp(0)
{

  // Constructor

  SetTriggerClasses();

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
  delete fhBeamType;
  delete fhL0TriggerInputMLL;
  delete fhL0TriggerInputMUL;
  delete fhL0TriggerInputMSL;
  delete fhL0TriggerInputTVX;
  delete fhL2TriggerCINT7CENTNOTRD;
  delete fhNMu;
  fhTriggers=NULL;
  fhNMu=NULL;
}

//====================================================================================================================================================

void AliAnalysisTaskLMREventFilter::SetTriggerClasses()
{
 
  fNTrigClass=7;
  fTriggerClasses[0]="-B-";
  fTriggerClasses[1]="CMSL7-B-";
  fTriggerClasses[2]="CMSH7-B-";
  fTriggerClasses[3]="CMUL7-B-";
  fTriggerClasses[4]="CMLL7-B-";
  fTriggerClasses[5]="C0TVX-B-NOPF-CENT";
  fTriggerClasses[6]="CINT7-B-NOPF-CENT";
 
  fL0TriggerInputMLL = 20; // reference to MLL L0 trigger
  fL0TriggerInputMUL = 21; // reference to MUL L0 trigger
  fL0TriggerInputMSL = 18; // reference to MSL L0 trigger
  fL0TriggerInputTVX =  3; // reference to TVX L0 trigger
 
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
 
  fhTriggers = new TH1D("hTriggers","L2 Triggers",1,0,1);
  fOutputList->Add(fhTriggers);	 
  fhTriggers->Sumw2();


  fhL0TriggerInputMLL = new TH1D("fhL0TriggerInputMLL","",121,-1.5,119.5);
  fOutputList->Add(fhL0TriggerInputMLL);
  fhL0TriggerInputMLL->Sumw2();

  fhL0TriggerInputMUL = new TH1D("fhL0TriggerInputMUL","",121,-1.5,119.5);
  fOutputList->Add(fhL0TriggerInputMUL);
  fhL0TriggerInputMUL->Sumw2();

  fhL0TriggerInputMSL = new TH1D("fhL0TriggerInputMSL","",121,-1.5,119.5);
  fOutputList->Add(fhL0TriggerInputMSL);
  fhL0TriggerInputMSL->Sumw2();

  fhL0TriggerInputTVX = new TH1D("fhL0TriggerInputTVX","",121,-1.5,119.5);
  fOutputList->Add(fhL0TriggerInputTVX);
  fhL0TriggerInputTVX->Sumw2();

  fhL2TriggerCINT7CENTNOTRD = new TH1D("fhL2TriggerCINT7CENTNOTRD","",121,-1.5,119.5);
  fOutputList->Add(fhL2TriggerCINT7CENTNOTRD);
  fhL2TriggerCINT7CENTNOTRD->Sumw2();

  fhBeamType = new TH1D("fhBeamType","",4,-0.5,3.5);
  fhBeamType->Sumw2();
  fhBeamType->GetXaxis()->SetBinLabel(1,"p-p");
  fhBeamType->GetXaxis()->SetBinLabel(2,"p-A");
  fhBeamType->GetXaxis()->SetBinLabel(3,"A-p");
  fhBeamType->GetXaxis()->SetBinLabel(4,"A-A");
  fOutputList->Add(fhBeamType);

  fhTriggers->GetXaxis()->SetBinLabel(1,fTriggerClasses[0]);

  for (Int_t i=1;i<fNTrigClass;i++)
    {
      fhTriggers->SetBins(fhTriggers->GetNbinsX()+1,0,1);
      fhTriggers->GetXaxis()->SetBinLabel(fhTriggers->GetNbinsX(),fTriggerClasses[i]);
      fhTriggers->SetBins(fhTriggers->GetNbinsX()+1,0,1);
      fhTriggers->GetXaxis()->SetBinLabel(fhTriggers->GetNbinsX(),Form("%s (PS)",fTriggerClasses[i].Data()));   
    }

  fhTriggers->SetBins(fhTriggers->GetNbinsX()+1,0,1);
  fhTriggers->GetXaxis()->SetBinLabel(fhTriggers->GetNbinsX(),"CMSL7-B- &0MUL (PS)");
  fhTriggers->SetBins(fhTriggers->GetNbinsX()+1,0,1);
  fhTriggers->GetXaxis()->SetBinLabel(fhTriggers->GetNbinsX(),"CMSL7-B- &0MLL (PS)");

 
  fhTriggers->SetBins(fhTriggers->GetNbinsX()+1,0,1);
  fhTriggers->GetXaxis()->SetBinLabel(fhTriggers->GetNbinsX(),"(CINT7-CENT || C0TVX-CENT) &0TVX (PS)");
  fhTriggers->SetBins(fhTriggers->GetNbinsX()+1,0,1);
  fhTriggers->GetXaxis()->SetBinLabel(fhTriggers->GetNbinsX(),"(CINT7-CENT || C0TVX-CENT) & (CINT7-CENT &0MSL) (PS)");


  fhTriggers->SetBins(fhTriggers->GetNbinsX()+1,0,1);
  fhTriggers->GetXaxis()->SetBinLabel(fhTriggers->GetNbinsX(),"CINT7-CENT &0MUL (PS)");
  fhTriggers->SetBins(fhTriggers->GetNbinsX()+1,0,1);
  fhTriggers->GetXaxis()->SetBinLabel(fhTriggers->GetNbinsX(),"CINT7-CENT &0MSL (PS)");



  fhNMu = new TH2D("hNMu","Number of Muon",20,0,20,fhTriggers->GetNbinsX(),0,1);
  for(Int_t i=1;i<fhTriggers->GetNbinsX();i++)
    fhNMu->GetYaxis()->SetBinLabel(i,fhTriggers->GetXaxis()->GetBinLabel(i));
  fOutputList->Add(fhNMu);
  fhNMu->Sumw2();

  PostData(1, fOutputList);
  PostData(2, fEventTree);
  printf("End of create Output\n");
}

//====================================================================================================================================================

void AliAnalysisTaskLMREventFilter::NotifyRun()
{
  fMuonTrackCuts->SetRun(fInputHandler);
  AliCDBManager *man = AliCDBManager::Instance();
  man->Init();
  man->SetDefaultStorage("raw://");
  man->SetRun(fInputHandler->GetEvent()->GetRunNumber());

  AliGRPObject* fGRPData = (AliGRPObject*) man->Get("GRP/GRP/Data")->GetObject();
  if(fGRPData->GetBeamType().Contains("p-p"))
    {
      fhBeamType->Fill("p-p",1);
      fminContributorsPileUp = 3;
    }
  else if(fGRPData->GetBeamType().Contains("p-A"))
    {
      fhBeamType->Fill("p-A",1);
      fminContributorsPileUp = 5;
    }
  else if(fGRPData->GetBeamType().Contains("A-p"))
    {
      fhBeamType->Fill("A-p",1);
      fminContributorsPileUp = 5;
    }
  else if(fGRPData->GetBeamType().Contains("A-A"))
    {
      fhBeamType->Fill("A-A",1);
      fminContributorsPileUp = 5;
    }
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
      if(name.Contains("0TVX"))
	fL0TriggerInputTVX = inp->GetIndexCTP(); // reference to TVX L0 trigger
    }
  fhL0TriggerInputMLL->Fill(fL0TriggerInputMLL);
  fhL0TriggerInputMUL->Fill(fL0TriggerInputMUL);
  fhL0TriggerInputMSL->Fill(fL0TriggerInputMSL);
  fhL0TriggerInputTVX->Fill(fL0TriggerInputTVX);

  Int_t fL2TriggerCINT7CENTNOTRD = cfg -> GetClassIndexFromName("CINT7-B-NOPF-CENTNOTRD");
  fhL2TriggerCINT7CENTNOTRD->Fill(fL2TriggerCINT7CENTNOTRD);
}
//====================================================================================================================================================

void AliAnalysisTaskLMREventFilter::UserExec(Option_t *)
 {
  //   Main loop
  //   Called for each event
  AliAODEvent *fAOD = dynamic_cast<AliAODEvent *>(InputEvent()); 
  if (!fAOD)
    return;

  UShort_t physicsSelectionMask=1<<0;
  UShort_t L0TriggerInput=1<<0;
 
  if (!IsSelectedTrigger(fAOD,physicsSelectionMask,L0TriggerInput))
    return;
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
  fAliLMREvent->SetL0TriggerInput(L0TriggerInput);
  fAliLMREvent->SetPhysicsSelectionMask(physicsSelectionMask);
  fAliLMREvent->SetIsPileupFromSPD(fAOD->IsPileupFromSPD(fminContributorsPileUp));
  fAliLMREvent->SetIsPileupFromSPDInMultBins(fAOD->IsPileupFromSPDInMultBins());
  Int_t nmu=0;
  if(fAOD->GetNumberOfTracks())
     nmu= fAOD->GetNumberOfMuonTracks();
  if (nmu>0)
    {
      Int_t ntotTr = fAOD->GetNumberOfTracks();
      AliLMRMuon *trk = NULL;
      for (Int_t itr=0; itr<ntotTr; itr++)
	{
	  AliAODTrack *track = dynamic_cast<AliAODTrack*>(fAOD->GetTrack(itr));
	  if(!track->IsMuonTrack()) continue;
   
	  Double_t p[3];
   
	  Short_t charge = track->Charge();
	  track->GetPxPyPz(p);
	  Double_t chi2Match = track->GetChi2MatchTrigger();
	  Int_t match = track->GetMatchTrigger();	
	  Double_t chi2 = track->Chi2perNDF();
	  Double_t rAbs = track->GetRAtAbsorberEnd();
   
	  TVector3 dcaAtVz  = fMuonTrackCuts->GetCorrectedDCA(track);
	  Double_t pTotMean = fMuonTrackCuts->GetAverageMomentum(track);
	  Double_t pDca = pTotMean * dcaAtVz.Mag();
	  // Create new Muon	
	  trk=fAliLMREvent->AddMuon();
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
 // fAOD=NULL;

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

Bool_t AliAnalysisTaskLMREventFilter::IsSelectedTrigger(AliAODEvent *fAOD, UShort_t &physicsSelectionMask,UShort_t &L0TriggerInput)
{
 
  Int_t nmu=0;
  if (fAOD->GetNumberOfTracks())
    nmu= fAOD->GetNumberOfMuonTracks();
 
  // -- Check Trigger
 
  TString trigStr(((AliAODHeader*) fAOD->GetHeader())->GetFiredTriggerClasses());
 
  if(!trigStr.Contains(fTriggerClasses[0]))    // if trigger word does not contain -B- the event is rejected
    return kFALSE;

  fhNMu->Fill(nmu,fTriggerClasses[0],1);
  fhTriggers->Fill(fTriggerClasses[0].Data(),1);

  UInt_t inpmask = fAOD->GetHeader()->GetL0TriggerInputs();       // L0 trigger mask
  UInt_t fPhysSelMask = fInputHandler->IsEventSelected();         // Physics Selection (PS) mask

  // -------------------- Setting the L0 mask of the event

  Int_t is0TVXfired = (inpmask & (1<<(fL0TriggerInputTVX-1)));
  Int_t is0MSLfired = (inpmask & (1<<(fL0TriggerInputMSL-1)));
  Int_t is0MULfired = (inpmask & (1<<(fL0TriggerInputMUL-1)));
  Int_t is0MLLfired = (inpmask & (1<<(fL0TriggerInputMLL-1)));

  if(is0TVXfired)
    L0TriggerInput|=1<<1;
  if(is0MSLfired)
    L0TriggerInput|=1<<2;
  if(is0MULfired)
    L0TriggerInput|=1<<3;
  if(is0MLLfired)
    L0TriggerInput|=1<<4;

  // -------------------- Setting the PS mask of the event

  Short_t kPSV0     = 1<<1;
  Short_t kPST0     = 1<<2;
  Short_t kPSV0Muon = 1<<3;

  if(fPhysSelMask&(AliVEvent::kINT7|AliVEvent::kINT7inMUON))
    physicsSelectionMask |= kPSV0;
  if(fPhysSelMask&AliVEvent::kAny)
    {
      if(trigStr.Contains(fTriggerClasses[5].Data()))
	physicsSelectionMask |= kPST0;
    }
  if (fPhysSelMask & (AliVEvent::kMuonSingleLowPt7|AliVEvent::kMuonSingleHighPt7|AliVEvent::kMuonLikeLowPt7|AliVEvent::kMuonUnlikeLowPt7))
    physicsSelectionMask |= kPSV0Muon;

  // -------------------- Filling the trigger information: WITHOUT REQUESTING PS

  for(Int_t i=1;i<fNTrigClass;i++)
    {
      if(trigStr.Contains(fTriggerClasses[i]))
	{
	  fhNMu->Fill(nmu, fTriggerClasses[i].Data(),1);
	  fhTriggers->Fill(fTriggerClasses[i].Data(),1);
	}   
    }

  // -------------------- Filling the trigger information: REQUESTING PS

  // The following trigger conditions are mainly used in the evaluation of the luminosity and the downscaling factors
 
  if(physicsSelectionMask & kPSV0Muon)      // kPSMuon is the PS to be used for the CMSL7-B-, CMSH7-B-, CMUL7-B- and CMLL7-B- trigger classes
    {
      for(Int_t i=1;i<5;i++)
	{
	  if(trigStr.Contains(fTriggerClasses[i].Data()))
	    {
	      fhTriggers->Fill(Form("%s (PS)",fTriggerClasses[i].Data()),1);
	      fhNMu ->Fill(nmu,Form("%s (PS)",fTriggerClasses[i].Data()),1);
	    }
	 }
      if(trigStr.Contains(fTriggerClasses[1].Data()))
	{
	  if(is0MULfired) // Used for both luminosity and downscaling
	    {
	      fhTriggers->Fill(Form("%s &0MUL (PS)",fTriggerClasses[1].Data()),1);
	      fhNMu ->Fill(nmu,Form("%s &0MUL (PS)",fTriggerClasses[1].Data()),1);
	    }
	  if(is0MLLfired) // Used for downscaling only
	    {
	      fhTriggers->Fill(Form("%s &0MLL (PS)",fTriggerClasses[1].Data()),1);
	      fhNMu ->Fill(nmu,Form("%s &0MLL (PS)",fTriggerClasses[1].Data()),1);
	    }
	}
    }

  if(physicsSelectionMask & kPSV0)      // kPSV0 is the PS to be used for the V0-based MB trigger classes (CINT7-B-*)
    {
      if(trigStr.Contains(fTriggerClasses[6].Data()))
	{
	  fhTriggers->Fill(Form("%s (PS)",fTriggerClasses[6].Data()),1);
	  fhNMu ->Fill(nmu,Form("%s (PS)",fTriggerClasses[6].Data()),1);
	  if(is0MSLfired) 
	    {	
	      fhTriggers->Fill("CINT7-CENT &0MSL (PS)",1);
	      fhNMu ->Fill(nmu,"CINT7-CENT &0MSL (PS)",1);
	    }
	  if(is0MULfired)
	    {	
	      fhTriggers->Fill("CINT7-CENT &0MUL (PS)",1);
	      fhNMu ->Fill(nmu,"CINT7-CENT &0MUL (PS)",1);
	    }
	}
    }

 
  if( (trigStr.Contains(fTriggerClasses[5].Data()) && (physicsSelectionMask & kPST0)) || ( trigStr.Contains(fTriggerClasses[6].Data()) && (physicsSelectionMask & kPSV0)) )     // PS adapted to include the T0 information for the MB trigger condition
    {
      if(is0TVXfired)
	{
	  fhTriggers->Fill("(CINT7-CENT || C0TVX-CENT) &0TVX (PS)",1);
	  fhNMu ->Fill(nmu,"(CINT7-CENT || C0TVX-CENT) &0TVX (PS)",1);
	}
      if(trigStr.Contains(fTriggerClasses[6].Data()) && (physicsSelectionMask & kPSV0))
	{
	  if(is0MSLfired)
	    {
	      fhTriggers->Fill("(CINT7-CENT || C0TVX-CENT) & (CINT7-CENT &0MSL) (PS)",1);
	      fhNMu ->Fill(nmu,"(CINT7-CENT || C0TVX-CENT) & (CINT7-CENT &0MSL) (PS)",1);
	    }
	}
    }

  return kTRUE;

}

//====================================================================================================================================================
