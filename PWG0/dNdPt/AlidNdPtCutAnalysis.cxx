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
//------------------------------------------------------------------------------
// AlidNdPtCutAnalysis class. 
//
// a. functionality:
// - fills generic cut histograms
// - generates cuts (selection criteria)
//
// b. data members:
// - generic cut histograms
// - control histograms
//
// Author: J.Otwinowski 04/11/2008 
//------------------------------------------------------------------------------
#include "TH1.h"
#include "TH2.h"

#include "AliHeader.h"  
#include "AliGenEventHeader.h"  
#include "AliStack.h"  
#include "AliESDEvent.h"  
#include "AliMCEvent.h"  
#include "AliESDtrackCuts.h"  
#include "AliLog.h" 
#include "AliTracker.h" 

#include "AlidNdPtEventCuts.h"
#include "AlidNdPtAcceptanceCuts.h"
#include "AlidNdPtBackgroundCuts.h"
#include "AlidNdPtAnalysis.h"
#include "AliPhysicsSelection.h"

#include "AliPWG0Helper.h"
#include "AlidNdPtHelper.h"
#include "AlidNdPtCutAnalysis.h"

using namespace std;

ClassImp(AlidNdPtCutAnalysis)

//_____________________________________________________________________________
  AlidNdPtCutAnalysis::AlidNdPtCutAnalysis(): AlidNdPt(),
  fAnalysisFolder(0),
  fEventCount(0),
  fRecEventHist(0),
  fMCEventHist(0),
  fRecMCEventHist(0),
  fRecMCTrackHist(0)
{
  // default constructor
  Init();
}

//_____________________________________________________________________________
AlidNdPtCutAnalysis::AlidNdPtCutAnalysis(Char_t* name, Char_t* title): AlidNdPt(name,title),
  fAnalysisFolder(0),
  fEventCount(0),
  fRecEventHist(0),
  fMCEventHist(0),
  fRecMCEventHist(0),
  fRecMCTrackHist(0)
{
  // constructor
  Init();
}

//_____________________________________________________________________________
AlidNdPtCutAnalysis::~AlidNdPtCutAnalysis() {
  // 
  if(fEventCount) delete fEventCount; fEventCount=0;
  if(fRecEventHist) delete fRecEventHist; fRecEventHist=0;
  if(fMCEventHist) delete fMCEventHist; fMCEventHist=0;
  if(fRecMCEventHist) delete fRecMCEventHist; fRecMCEventHist=0;
  if(fRecMCTrackHist) delete fRecMCTrackHist; fRecMCTrackHist=0;

  if(fAnalysisFolder) delete fAnalysisFolder; fAnalysisFolder=0;
}

//_____________________________________________________________________________
void AlidNdPtCutAnalysis::Init(){
  //
  // Init histograms
  //
  const Int_t ptNbins = 56; 
  const Double_t ptMin = 0.; 
  const Double_t ptMax = 16.; 

  Double_t binsPt[ptNbins+1] = {0.,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.2,2.4,2.6,2.8,3.0,3.2,3.4,3.6,3.8,4.0,4.5,5.0,5.5,6.0,6.5,7.0,7.5,8.0,9.0,10.0,11.0,12.0,13.0,14.0,15.0,16.0};

  // 
  Int_t binsEventCount[2]={2,2};
  Double_t minEventCount[2]={0,0}; 
  Double_t maxEventCount[2]={2,2}; 
  fEventCount = new THnSparseF("fEventCount","trig vs trig+vertex",2,binsEventCount,minEventCount,maxEventCount);
  fEventCount->GetAxis(0)->SetTitle("trig");
  fEventCount->GetAxis(1)->SetTitle("trig+vert");
  fEventCount->Sumw2();

  //Xv:Yv:Zv:ResZv:Mult
  Double_t kFact = 1.0;

  Int_t binsRecEventHist[5]={80,80,100,80,150};
  Double_t minRecEventHist[5]={-3.*kFact,-3.*kFact,-35.,0.,0.}; 
  Double_t maxRecEventHist[5]={3.*kFact,3.*kFact,35.,10.,150.}; 
  fRecEventHist = new THnSparseF("fRecEventHist","Xv:Yv:Zv:ResZv:Mult",5,binsRecEventHist,minRecEventHist,maxRecEventHist);
  fRecEventHist->GetAxis(0)->SetTitle("Xv (cm)");
  fRecEventHist->GetAxis(1)->SetTitle("Yv (cm)");
  fRecEventHist->GetAxis(2)->SetTitle("Zv (cm)");
  fRecEventHist->GetAxis(3)->SetTitle("ResZv (cm)");
  fRecEventHist->GetAxis(4)->SetTitle("Mult");
  fRecEventHist->Sumw2();

  //Xv:Yv:Zv
  Int_t binsMCEventHist[3]={80,80,100};
  Double_t minMCEventHist[3]={-0.1,-0.1,-35.}; 
  Double_t maxMCEventHist[3]={0.1,0.1,35.}; 
  fMCEventHist = new THnSparseF("fMCEventHist","mcXv:mcYv:mcZv",3,binsMCEventHist,minMCEventHist,maxMCEventHist);
  fMCEventHist->GetAxis(0)->SetTitle("mcXv (cm)");
  fMCEventHist->GetAxis(1)->SetTitle("mcYv (cm)");
  fMCEventHist->GetAxis(2)->SetTitle("mcZv (cm)");
  fMCEventHist->Sumw2();

  //Xv-mcXv:Yv-mcYv:Zv-mcZv:Mult
  Int_t binsRecMCEventHist[4]={100,100,100,150};
  Double_t minRecMCEventHist[4]={-1.0*kFact,-1.0*kFact,-1.0*kFact,0.}; 
  Double_t maxRecMCEventHist[4]={1.0*kFact,1.0*kFact,1.0*kFact,150.}; 
  fRecMCEventHist = new THnSparseF("fRecMCEventHist","Xv-mcXv:Yv-mcYv:Zv-mcZv:Mult",4,binsRecMCEventHist,minRecMCEventHist,maxRecMCEventHist);
  fRecMCEventHist->GetAxis(0)->SetTitle("Xv-mcXv (cm)");
  fRecMCEventHist->GetAxis(1)->SetTitle("Yv-mcYv (cm)");
  fRecMCEventHist->GetAxis(2)->SetTitle("Zv-mcZv (cm)");
  fRecMCEventHist->GetAxis(3)->SetTitle("Mult");
  fRecMCEventHist->Sumw2();

  //
  // THnSparse track histograms
  //

  //nClust:chi2PerClust:nClust/nFindableClust:DCAy:DCAz:eta:phi:pt:kinkIdx:isPrim:polarity
  Int_t binsRecMCTrackHist[11]={160,80,80,100,100,30,90,ptNbins, 3, 2, 2};
  Double_t minRecMCTrackHist[11]={0., 0., 0., -5.,-5.,-1.5, 0., ptMin, -1., 0., 0.};
  Double_t maxRecMCTrackHist[11]={160.,10.,1.2, 5.,5.,1.5, 2.*TMath::Pi(), ptMax, 2., 2., 2.};

  fRecMCTrackHist = new THnSparseF("fRecMCTrackHist","nClust:chi2PerClust:nClust/nFindableClust:DCAy:DCAz:eta:phi:pt:kinkIdx:isPrim:polarity",11,binsRecMCTrackHist,minRecMCTrackHist,maxRecMCTrackHist);
  fRecMCTrackHist->SetBinEdges(7,binsPt);

  fRecMCTrackHist->GetAxis(0)->SetTitle("nClust");
  fRecMCTrackHist->GetAxis(1)->SetTitle("chi2PerClust");
  fRecMCTrackHist->GetAxis(2)->SetTitle("nClust/nFindableClust");
  fRecMCTrackHist->GetAxis(3)->SetTitle("DCAy (cm)");
  fRecMCTrackHist->GetAxis(4)->SetTitle("DCAz (cm)");
  fRecMCTrackHist->GetAxis(5)->SetTitle("#eta");
  fRecMCTrackHist->GetAxis(6)->SetTitle("#phi (rad)");
  fRecMCTrackHist->GetAxis(7)->SetTitle("p_{T} (GeV/c)");
  fRecMCTrackHist->GetAxis(8)->SetTitle("kinkIdx"); // 0 - no kink, -1 - kink mother, 1 - kink daugther 
  fRecMCTrackHist->GetAxis(9)->SetTitle("isPrim");
  fRecMCTrackHist->GetAxis(10)->SetTitle("polarity");
  fRecMCTrackHist->Sumw2();

  // init output folder
  fAnalysisFolder = CreateFolder("folderdNdPt","Analysis dNdPt Folder");

}

//_____________________________________________________________________________
void AlidNdPtCutAnalysis::Process(AliESDEvent *const esdEvent, AliMCEvent * const mcEvent)
{
  //
  // Process real and/or simulated events
  //
  if(!esdEvent) {
    AliDebug(AliLog::kError, "esdEvent not available");
    return;
  }

  // get selection cuts
  AlidNdPtEventCuts *evtCuts = GetEventCuts(); 
  AlidNdPtAcceptanceCuts *accCuts = GetAcceptanceCuts(); 
  AliESDtrackCuts *esdTrackCuts = GetTrackCuts(); 

  if(!evtCuts || !accCuts  || !esdTrackCuts) {
    AliDebug(AliLog::kError, "cuts not available");
    return;
  }

  // trigger selection
  Bool_t isEventTriggered = kTRUE;
  AliPhysicsSelection *trigSel = NULL;
  AliTriggerAnalysis *trigAna = NULL; // needed for andV0

  if(evtCuts->IsTriggerRequired())  
  {
    //
    trigSel = GetPhysicsTriggerSelection();
    if(!trigSel) {
      AliDebug(AliLog::kError, "cannot get trigSel");
      return;
    }
    
    if(IsUseMCInfo()) { 
      trigSel->SetAnalyzeMC();
      isEventTriggered = trigSel->IsCollisionCandidate(esdEvent);

      trigAna = trigSel->GetTriggerAnalysis();
      if(!trigAna) 
        return;

      if(GetTrigger() == AliTriggerAnalysis::kV0AND)
        isEventTriggered = trigAna->IsOfflineTriggerFired(esdEvent, GetTrigger());
    }
    else {
      //
      // 0-multiplicity bin for LHC background correction
      //
      if(GetAnalysisMode() == AlidNdPtHelper::kTPCTrackSPDvtx || GetAnalysisMode() == AlidNdPtHelper::kTPCTrackSPDvtxUpdate || 
         GetAnalysisMode() == AlidNdPtHelper::kTPCITSHybridTrackSPDvtx || GetAnalysisMode() == AlidNdPtHelper::kTPCITSHybridTrackSPDvtxDCArPt) 
      {
        trigSel->SetBin0CallbackViaPointer(&AlidNdPtAnalysis::IsBinZeroTrackSPDvtx);
      } else {
        trigSel->SetBin0CallbackViaPointer(&AlidNdPtAnalysis::IsBinZeroSPDvtx);
      }

      isEventTriggered = trigSel->IsCollisionCandidate(esdEvent);

      trigAna = trigSel->GetTriggerAnalysis();
      if(!trigAna) 
        return;

      if(GetTrigger() == AliTriggerAnalysis::kV0AND)
        isEventTriggered = trigAna->IsOfflineTriggerFired(esdEvent, GetTrigger());
    }
  }

  // use MC information
  AliHeader* header = 0;
  AliGenEventHeader* genHeader = 0;
  AliStack* stack = 0;
  TArrayF vtxMC(3);
  AliPWG0Helper::MCProcessType evtType = AliPWG0Helper::kInvalidProcess;

  if(IsUseMCInfo())
  {
    if(!mcEvent) {
      AliDebug(AliLog::kError, "mcEvent not available");
      return;
    }

    // get MC event header
    header = mcEvent->Header();
    if (!header) {
      AliDebug(AliLog::kError, "Header not available");
      return;
    }

    // MC particle stack
    stack = mcEvent->Stack();
    if (!stack) {
      AliDebug(AliLog::kError, "Stack not available");
      return;
    }

    // get event type (ND=0x1, DD=0x2, SD=0x4)
    evtType = AliPWG0Helper::GetEventProcessType(header);
    AliDebug(AliLog::kDebug+1, Form("Found process type %d", evtType));

    // get MC vertex
    genHeader = header->GenEventHeader();
    if (!genHeader) {
      AliDebug(AliLog::kError, "Could not retrieve genHeader from Header");
      return;
    }
    genHeader->PrimaryVertex(vtxMC);

    // Fill MC event histogram
    Double_t vMCEventHist[3]={vtxMC[0],vtxMC[1],vtxMC[2]};
    fMCEventHist->Fill(vMCEventHist);

  } // end bUseMC

  // get reconstructed vertex  
  Bool_t bRedoTPCVertex = evtCuts->IsRedoTPCVertex();
  Bool_t bUseConstraints = evtCuts->IsUseBeamSpotConstraint();
  const AliESDVertex* vtxESD = AlidNdPtHelper::GetVertex(esdEvent,evtCuts,accCuts,esdTrackCuts,GetAnalysisMode(),kFALSE,bRedoTPCVertex,bUseConstraints); 
  Bool_t isRecVertex = AlidNdPtHelper::TestRecVertex(vtxESD, esdEvent->GetPrimaryVertexSPD(), GetAnalysisMode(), kFALSE);
  Bool_t isEventOK = evtCuts->AcceptEvent(esdEvent,mcEvent,vtxESD) && isRecVertex;

  TObjArray *allChargedTracks=0;
  Int_t multAll=0;
  
  //
  // event counter
  // 
  //printf("isEventOK %d, isEventTriggered %d \n",isEventOK,isEventTriggered);

  Bool_t isTrigAndVertex = isEventTriggered && isEventOK;
  Double_t vEventCount[2] = { isEventTriggered, isTrigAndVertex};
  fEventCount->Fill(vEventCount);

  //
  // cosmic background and splitted tracks
  //
  if(GetParticleMode() == AlidNdPtHelper::kBackgroundTrack) 
  {
    AlidNdPtBackgroundCuts *backCuts = GetBackgroundCuts(); 
    if(!backCuts) return;

    for (Int_t iTrack = 0; iTrack < esdEvent->GetNumberOfTracks(); iTrack++) 
    {
      AliESDtrack *track1 = esdEvent->GetTrack(iTrack);
      if(!track1) continue; 
      if(track1->Charge()==0) continue; 

      for (Int_t jTrack = iTrack+1; jTrack < esdEvent->GetNumberOfTracks(); jTrack++) 
      {
        AliESDtrack *track2 = esdEvent->GetTrack(jTrack);
        if(!track2) continue; 
        if(track2->Charge()==0) continue; 

	//printf("track2->Charge() %d\n",track2->Charge());

        backCuts->IsBackgroundTrack(track1, track2);
      }
    }
  }

  // check event cuts
  if(isEventOK && isEventTriggered)
  {
    // get all charged tracks
    allChargedTracks = AlidNdPtHelper::GetAllChargedTracks(esdEvent,GetAnalysisMode());
    if(!allChargedTracks) return;

    Int_t entries = allChargedTracks->GetEntries();
    for(Int_t i=0; i<entries;++i) 
    {
      AliESDtrack *track = (AliESDtrack*)allChargedTracks->At(i);
      if(!track) continue;

      if(!esdTrackCuts->AcceptTrack(track)) continue;

      //
      Bool_t isOK = kFALSE;
      Double_t x[3]; track->GetXYZ(x);
      Double_t b[3]; AliTracker::GetBxByBz(x,b);

      //
      // if TPC-ITS hybrid tracking (kTPCITSHybrid)
      // replace track parameters with TPC-ony track parameters
      //
      if( GetAnalysisMode() == AlidNdPtHelper::kTPCITSHybrid || GetAnalysisMode() == AlidNdPtHelper::kTPCITSHybridTrackSPDvtx || GetAnalysisMode() == AlidNdPtHelper::kTPCITSHybridTrackSPDvtxDCArPt) 
      {
        // Relate TPC-only tracks to SPD vertex
        isOK = track->RelateToVertexTPCBxByBz(vtxESD, b, kVeryBig);
        if(!isOK) continue;

	// replace esd track parameters with TPCinner
        AliExternalTrackParam  *tpcTrack  = new AliExternalTrackParam(*(track->GetTPCInnerParam()));
	if (!tpcTrack) return;
        track->Set(tpcTrack->GetX(),tpcTrack->GetAlpha(),tpcTrack->GetParameter(),tpcTrack->GetCovariance());

        if(tpcTrack) delete tpcTrack; 
      } 

      //
      if (GetAnalysisMode()==AlidNdPtHelper::kTPCSPDvtxUpdate || GetAnalysisMode() == AlidNdPtHelper::kTPCTrackSPDvtxUpdate) 
      {
        //
        // update track parameters
	//
        AliExternalTrackParam cParam;
        isOK = track->RelateToVertexTPCBxByBz(vtxESD, b, kVeryBig, &cParam);
	if(!isOK) continue;

	track->Set(cParam.GetX(),cParam.GetAlpha(),cParam.GetParameter(),cParam.GetCovariance());
      }

      FillHistograms(track, stack);
      multAll++;
    }

    Double_t vRecEventHist[5] = {vtxESD->GetXv(),vtxESD->GetYv(),vtxESD->GetZv(),vtxESD->GetZRes(),multAll};
    fRecEventHist->Fill(vRecEventHist);

    if(IsUseMCInfo()) {
      Double_t vRecMCEventHist[5] = {vtxESD->GetXv()-vtxMC[0],vtxESD->GetYv()-vtxMC[1],vtxESD->GetZv()-vtxMC[2],multAll};
      fRecMCEventHist->Fill(vRecMCEventHist);
    }
  }

  if(allChargedTracks) delete allChargedTracks; allChargedTracks = 0;

}

//_____________________________________________________________________________
void AlidNdPtCutAnalysis::FillHistograms(AliESDtrack *const esdTrack, AliStack *const stack) const
{
  //
  // Fill ESD track and MC histograms 
  //
  if(!esdTrack) return;
  if(esdTrack->Charge() == 0.) return;

  Float_t pt = esdTrack->Pt();
  Float_t eta = esdTrack->Eta();
  Float_t phi = esdTrack->Phi();
  //Int_t nClust = esdTrack->GetTPCclusters(0);
  Int_t nClust = esdTrack->GetTPCNclsIter1();
  Int_t nFindableClust = esdTrack->GetTPCNclsF();

  Float_t chi2PerCluster = 0.;
  //if(nClust>0.) chi2PerCluster = esdTrack->GetTPCchi2()/Float_t(nClust);
  if(nClust>0.) chi2PerCluster = esdTrack->GetTPCchi2Iter1()/Float_t(nClust);

  Float_t clustPerFindClust = 0.;
  if(nFindableClust>0.) clustPerFindClust = Float_t(nClust)/nFindableClust;

  Float_t b[2], bCov[3];
  esdTrack->GetImpactParameters(b,bCov);

  // kink idx
  Int_t kinkIdx = 0;
  //if(esdTrack->GetKinkIndex(0) > 0.)   isKink  = kTRUE;
  if(esdTrack->GetKinkIndex(0) > 0)      kinkIdx = 1;   // kink daughter
  else if(esdTrack->GetKinkIndex(0) < 0) kinkIdx = -1;  // kink mother
  else kinkIdx = 0; // not kink

  //printf("esdTrack->GetKinkIndex(0) %d \n", esdTrack->GetKinkIndex(0));
  //printf("esdTrack->GetKinkIndex(1) %d \n", esdTrack->GetKinkIndex(1));
  //printf("esdTrack->GetKinkIndex(2) %d \n", esdTrack->GetKinkIndex(2));
  //printf("kinkIdx %d \n", kinkIdx);

  //
  // Fill rec vs MC information
  //
  Bool_t isPrim = kTRUE;

  if(IsUseMCInfo()) {
    if(!stack) return;
    Int_t label = TMath::Abs(esdTrack->GetLabel()); 
    TParticle* particle = stack->Particle(label);
    if(!particle) return;
    if(particle->GetPDG() && particle->GetPDG()->Charge()==0.) return;
    isPrim = stack->IsPhysicalPrimary(label);

    //if(isPrim && pt > 1.5 && kinkIdx == -1) printf("nClust  %d \n", nClust);
  }
  
  // fill histo
  Int_t polarity = -2;
  if (esdTrack->Charge() < 0.) polarity = 0; 
  else polarity = 1; 
  Double_t vRecMCTrackHist[11] = { nClust,chi2PerCluster,clustPerFindClust,b[0],b[1],eta,phi,pt,kinkIdx,isPrim, polarity }; 
  fRecMCTrackHist->Fill(vRecMCTrackHist);
}

//_____________________________________________________________________________
Long64_t AlidNdPtCutAnalysis::Merge(TCollection* const list) 
{
  // Merge list of objects (needed by PROOF)

  if (!list)
  return 0;

  if (list->IsEmpty())
  return 1;

  TIterator* iter = list->MakeIterator();
  TObject* obj = 0;

  //TList *collPhysSelection = new TList;

  // collection of generated histograms
  Int_t count=0;
  while((obj = iter->Next()) != 0) {
    AlidNdPtCutAnalysis* entry = dynamic_cast<AlidNdPtCutAnalysis*>(obj);
    if (entry == 0) continue; 
  
    // event histo
    fEventCount->Add(entry->fEventCount);
    fRecEventHist->Add(entry->fRecEventHist);
    fRecMCEventHist->Add(entry->fRecMCEventHist);
    fMCEventHist->Add(entry->fMCEventHist);

    // track histo
    fRecMCTrackHist->Add(entry->fRecMCTrackHist);

    // physics selection
    //collPhysSelection->Add(entry->GetPhysicsTriggerSelection());
    
  count++;
  }

  //AliPhysicsSelection *trigSelection = GetPhysicsTriggerSelection();
  //trigSelection->Merge(collPhysSelection);

  //if(collPhysSelection) delete collPhysSelection;

return count;
}

//_____________________________________________________________________________
void AlidNdPtCutAnalysis::Analyse() 
{
  //
  // Analyse histograms
  //
  TH1::AddDirectory(kFALSE);
  TObjArray *aFolderObj = new TObjArray;
  TH1D *h1D = 0; 
  TH2D *h2D = 0; 


  //
  // get cuts
  //
  AlidNdPtEventCuts *evtCuts = GetEventCuts(); 
  AlidNdPtAcceptanceCuts *accCuts = GetAcceptanceCuts(); 
  AliESDtrackCuts *esdTrackCuts = GetTrackCuts(); 

  if(!evtCuts || !accCuts || !esdTrackCuts) {
    Error("AlidNdPtCutAnalysis::Analyse()", "cuts not available");
    return;
  }

  //
  // set min and max values
  //
  Double_t minPt = accCuts->GetMinPt();
  Double_t maxPt = accCuts->GetMaxPt();
  Double_t minEta = accCuts->GetMinEta();
  Double_t maxEta = accCuts->GetMaxEta()-0.00001;

  Double_t maxDCAr = accCuts->GetMaxDCAr();

  //
  // Event counters
  //
  h2D = (TH2D*)fEventCount->Projection(0,1);
  h2D->SetName("trig_vs_trigANDvertex");
  aFolderObj->Add(h2D);

  fEventCount->GetAxis(0)->SetRange(2,2); // triggered
  h1D = (TH1D*)fEventCount->Projection(1);
  h1D->SetTitle("rec. vertex for triggered events");
  h1D->SetName("trigANDvertex");
  aFolderObj->Add(h1D);

  //
  // Create rec. event histograms
  //
  h1D = (TH1D *)fRecEventHist->Projection(0);
  h1D->SetName("rec_xv");
  aFolderObj->Add(h1D);

  h1D = (TH1D *)fRecEventHist->Projection(1);
  h1D->SetName("rec_yv");
  aFolderObj->Add(h1D);

  h1D = (TH1D *)fRecEventHist->Projection(2);
  h1D->SetName("rec_zv");
  aFolderObj->Add(h1D);

  h2D = (TH2D *)fRecEventHist->Projection(3,4);
  h2D->SetName("rec_resZv_vs_Mult");
  aFolderObj->Add(h2D);

  h2D = (TH2D *)fRecEventHist->Projection(0,1);
  h2D->SetName("rec_xv_vs_yv");
  aFolderObj->Add(h2D);

  h2D = (TH2D *)fRecEventHist->Projection(0,2);
  h2D->SetName("rec_xv_vs_zv");
  aFolderObj->Add(h2D);

  h2D = (TH2D *)fRecEventHist->Projection(3,4);
  h2D->SetName("rec_resZv_vs_Mult");
  aFolderObj->Add(h2D);

  //
  // MC available
  //
  if(IsUseMCInfo()) {

  //
  // Create mc event histograms
  //
  h2D = (TH2D *)fMCEventHist->Projection(0,1);
  h2D->SetName("mc_xv_vs_yv");
  aFolderObj->Add(h2D);

  h2D = (TH2D *)fMCEventHist->Projection(0,2);
  h2D->SetName("mc_xv_vs_zv");
  aFolderObj->Add(h2D);

  //
  // Create rec-mc event histograms
  //
  h2D = (TH2D *)fRecMCEventHist->Projection(0,3);
  h2D->SetName("rec_mc_deltaXv_vs_mult");
  aFolderObj->Add(h2D);

  h2D = (TH2D *)fRecMCEventHist->Projection(1,3);
  h2D->SetName("rec_mc_deltaYv_vs_mult");
  aFolderObj->Add(h2D);

  h2D = (TH2D *)fRecMCEventHist->Projection(2,3);
  h2D->SetName("rec_mc_deltaZv_vs_mult");
  aFolderObj->Add(h2D);

  } // end use MC info 



  //
  // Create rec-mc track track histograms 
  //

  // DCA cuts
  fRecMCTrackHist->GetAxis(3)->SetRangeUser(-maxDCAr,maxDCAr);
  fRecMCTrackHist->GetAxis(4)->SetRangeUser(-maxDCAr,maxDCAr);

  h2D = (TH2D *)fRecMCTrackHist->Projection(7,5);
  h2D->SetName("pt_vs_eta");
  aFolderObj->Add(h2D);

  fRecMCTrackHist->GetAxis(7)->SetRangeUser(minPt,maxPt);  

  h2D = (TH2D *)fRecMCTrackHist->Projection(0,5);
  h2D->SetName("nClust_vs_eta");
  aFolderObj->Add(h2D);

  h2D = (TH2D *)fRecMCTrackHist->Projection(1,5);
  h2D->SetName("chi2PerClust_vs_eta");
  aFolderObj->Add(h2D);

  h2D = (TH2D *)fRecMCTrackHist->Projection(2,5);
  h2D->SetName("ratio_nClust_nFindableClust_vs_eta");
  aFolderObj->Add(h2D);

  h2D = (TH2D *)fRecMCTrackHist->Projection(5,6);
  h2D->SetName("eta_vs_phi");
  aFolderObj->Add(h2D);

  //
  fRecMCTrackHist->GetAxis(5)->SetRangeUser(minEta,maxEta);  

  h2D = (TH2D *)fRecMCTrackHist->Projection(0,6);
  h2D->SetName("nClust_vs_phi");
  aFolderObj->Add(h2D);

  h2D = (TH2D *)fRecMCTrackHist->Projection(1,6);
  h2D->SetName("chi2PerClust_vs_phi");
  aFolderObj->Add(h2D);

  h2D = (TH2D *)fRecMCTrackHist->Projection(2,6);
  h2D->SetName("ratio_nClust_nFindableClust_vs_phi");
  aFolderObj->Add(h2D);

  //
  fRecMCTrackHist->GetAxis(7)->SetRangeUser(0.0,maxPt);  

  h2D = (TH2D *)fRecMCTrackHist->Projection(0,7);
  h2D->SetName("nClust_vs_pt");
  aFolderObj->Add(h2D);

  h2D = (TH2D *)fRecMCTrackHist->Projection(1,7);
  h2D->SetName("chi2PerClust_vs_pt");
  aFolderObj->Add(h2D);

  h2D = (TH2D *)fRecMCTrackHist->Projection(2,7);
  h2D->SetName("ratio_nClust_nFindableClust_vs_pt");
  aFolderObj->Add(h2D);

  h2D = (TH2D *)fRecMCTrackHist->Projection(6,7);
  h2D->SetName("phi_vs_pt");
  aFolderObj->Add(h2D);


  // fiducial volume
  fRecMCTrackHist->GetAxis(5)->SetRangeUser(minEta,maxEta);  
  fRecMCTrackHist->GetAxis(7)->SetRangeUser(minPt,maxPt);  

  // DCA cuts
  fRecMCTrackHist->GetAxis(3)->SetRangeUser(-maxDCAr,maxDCAr);
  fRecMCTrackHist->GetAxis(4)->SetRangeUser(-maxDCAr,maxDCAr);

  h2D = (TH2D *)fRecMCTrackHist->Projection(0,1);
  h2D->SetName("nClust_vs_chi2PerClust");
  aFolderObj->Add(h2D);

  h2D = (TH2D *)fRecMCTrackHist->Projection(0,2);
  h2D->SetName("nClust_vs_ratio_nClust_nFindableClust");
  aFolderObj->Add(h2D);

  //
  // DCAy cuts
  //
  fRecMCTrackHist->GetAxis(0)->SetRange(50,160); // nClust/track > 50
  fRecMCTrackHist->GetAxis(1)->SetRangeUser(0.,3.9999); // chi2/cluster < 4.0
  fRecMCTrackHist->GetAxis(3)->SetRange(1,fRecMCTrackHist->GetAxis(3)->GetNbins());
  //fRecMCTrackHist->GetAxis(4)->SetRangeUser(-1.0,1.0);
  fRecMCTrackHist->GetAxis(4)->SetRange(1,fRecMCTrackHist->GetAxis(4)->GetNbins());

  // sec
  fRecMCTrackHist->GetAxis(9)->SetRange(1,1);
  h1D = (TH1D *)fRecMCTrackHist->Projection(3);
  h1D->SetName("dcay_sec");
  aFolderObj->Add(h1D);

  // prim
  fRecMCTrackHist->GetAxis(9)->SetRange(2,2);
  h1D = (TH1D *)fRecMCTrackHist->Projection(3);
  h1D->SetName("dcay_prim");
  aFolderObj->Add(h1D);

  // DCAz cuts
  //fRecMCTrackHist->GetAxis(3)->SetRangeUser(-1.0,1.0);
  fRecMCTrackHist->GetAxis(4)->SetRange(1,fRecMCTrackHist->GetAxis(4)->GetNbins());

  // sec
  fRecMCTrackHist->GetAxis(9)->SetRange(1,1);
  h1D = (TH1D *)fRecMCTrackHist->Projection(4);
  h1D->SetName("dcaz_sec");
  aFolderObj->Add(h1D);

  // prim
  fRecMCTrackHist->GetAxis(9)->SetRange(2,2);
  h1D = (TH1D *)fRecMCTrackHist->Projection(4);
  h1D->SetName("dcaz_prim");
  aFolderObj->Add(h1D);


  // export objects to analysis folder
  fAnalysisFolder = ExportToFolder(aFolderObj);

  // delete only TObjArray
  if(aFolderObj) delete aFolderObj;
}

//_____________________________________________________________________________
TFolder* AlidNdPtCutAnalysis::ExportToFolder(TObjArray * const array) 
{
  // recreate folder avery time and export objects to new one
  //
  AlidNdPtCutAnalysis * comp=this;
  TFolder *folder = comp->GetAnalysisFolder();

  TString name, title;
  TFolder *newFolder = 0;
  Int_t i = 0;
  Int_t size = array->GetSize();

  if(folder) { 
     // get name and title from old folder
     name = folder->GetName();  
     title = folder->GetTitle();  

	 // delete old one
     delete folder;

	 // create new one
     newFolder = CreateFolder(name.Data(),title.Data());
     newFolder->SetOwner();

	 // add objects to folder
     while(i < size) {
	   newFolder->Add(array->At(i));
	   i++;
	 }
  }

return newFolder;
}

//_____________________________________________________________________________
TFolder* AlidNdPtCutAnalysis::CreateFolder(TString name,TString title) { 
// create folder for analysed histograms
//
TFolder *folder = 0;
  folder = new TFolder(name.Data(),title.Data());

  return folder;
}
