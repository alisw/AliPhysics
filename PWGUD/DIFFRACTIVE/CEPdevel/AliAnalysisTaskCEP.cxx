/*************************************************************************
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
//
// Select events according to gap conditions
//
// Author:
//  Paul Buehler <paul.buehler@oeaw.ac.at>

#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TList.h>
#include <TObject.h>
#include <TRandom3.h>

#include "AliAnalysisManager.h"
#include "AliProdInfo.h"
#include "AliESDInputHandler.h"
#include "AliAODEvent.h"
#include "AliAnalysisUtils.h"
#include "AliRawEventHeaderBase.h"
#include "AliTriggerAnalysis.h"
#include "AliPIDResponse.h"
#include "AliPIDCombined.h"
#include "AliPhysicsSelection.h"
#include "AliMCEvent.h"
#include "AliStack.h"
#include "AliESDtrack.h"
#include "AliESDtrackCuts.h"
#include "AliMultiplicitySelectionCP.h"

#include "AliAnalysisTaskCEP.h"
#include "AliCEPBase.h"
#include "AliCEPUtils.h"

//------------------------------------------------------------------------------
// standard constructor (the one which should be used)
AliAnalysisTaskCEP::AliAnalysisTaskCEP(const char* name,
  Long_t state,
  Int_t rnummin, Int_t rnummax,
  Int_t numTracksMax,
  Double_t fracDG, Double_t fracNDG,
  UInt_t ETmaskDG, UInt_t ETpatternDG,
  UInt_t ETmaskNDG, UInt_t ETpatternNDG,
  UInt_t TTmask, UInt_t TTpattern):
	AliAnalysisTaskSE(name)
	, fAnalysisStatus(state)
  , frnummin(rnummin)
  , frnummax(rnummax)
  , fnumTracksMax(numTracksMax)
  , ffracDG(fracDG)
  , ffracNDG(fracNDG)
  , fETmaskDG(ETmaskDG)
  , fETpatternDG(ETpatternDG)
  , fETmaskNDG(ETmaskNDG)
  , fETpatternNDG(ETpatternNDG)
  , fTTmask(TTmask)
  , fTTpattern(TTpattern)
  , fLHCPeriod(TString(""))
  , fRun(-1)
  , fESDRun(0x0)
  , fESDEvent(0x0)
  , fCEPEvent(0x0)
  , fTracks(0x0)
  , fTrackStatus(0x0)
  , fVtxPos(TVector3(-999.9,-999.9,-999.9))
  , fPIDResponse(0x0)
  , fPIDCombined1(0x0)
  , fPIDCombined2(0x0)
  , fPIDCombined3(0x0)
	, fTrigger(0x0)
	, fPhysicsSelection(0x0)
  , fEventCuts(0x0)
  , fTrackCuts(0x0)
  , fMartinSel(0x0)
  , fMCCEPSystem(TLorentzVector(0,0,0,0))
	, flSPDpileup(0x0)
	, flnClunTra(0x0)
	, flVtx(0x0)
	, fhStatsFlow(0x0)
	, fHist(new TList())
	, fCEPtree(0x0)
  , fCEPUtil(0x0)
{

	// ensures that the histograms are all deleted on exit!
	fHist->SetOwner();

	// slot in TaskSE must start from 1
	Int_t iOutputSlot = 1;
	DefineOutput(iOutputSlot++, TList::Class());
  DefineOutput(iOutputSlot++, TTree::Class());
  
}


//------------------------------------------------------------------------------
AliAnalysisTaskCEP::AliAnalysisTaskCEP():
	AliAnalysisTaskSE()
	, fAnalysisStatus(AliCEPBase::kBitConfigurationSet)
  , frnummin(100000)
  , frnummax(300000)
  , fnumTracksMax(6)
  , ffracDG(1.0)
  , ffracNDG(0.0)
  , fETmaskDG(AliCEPBase::kETBaseLine)
  , fETpatternDG(AliCEPBase::kETBaseLine)
  , fETmaskNDG(AliCEPBase::kETBaseLine)
  , fETpatternNDG(AliCEPBase::kETBaseLine)
  , fTTmask(AliCEPBase::kTTBaseLine)
  , fTTpattern(AliCEPBase::kTTBaseLine)
  , fLHCPeriod(TString(""))
  , fRun(-1)
  , fESDEvent(0x0)
  , fCEPEvent(0x0)
  , fTracks(0x0)
  , fTrackStatus(0x0)
  , fVtxPos(TVector3(-999.9,-999.9,-999.9))
  , fPIDResponse(0x0)
  , fPIDCombined1(0x0)
  , fPIDCombined2(0x0)
  , fPIDCombined3(0x0)
	, fTrigger(0x0)
	, fPhysicsSelection(0x0)
  , fEventCuts(0x0)
	, fTrackCuts(0x0)
  , fMartinSel(0x0)
  , fMCCEPSystem(TLorentzVector(0,0,0,0))
  , flQArnum(0x0)
  , flBBFlag(0x0)
	, flSPDpileup(0x0)
	, flnClunTra(0x0)
	, flVtx(0x0)
  , flV0(0x0)
  , flFMD(0x0)
	, fhStatsFlow(0x0)
	, fHist(new TList())
	, fCEPtree(0x0)
  , fCEPUtil(0x0)
{

}


//------------------------------------------------------------------------------
AliAnalysisTaskCEP::~AliAnalysisTaskCEP()
{
	// Destructor
  
  // delet all objects created with new  
  // delete CEPEventBuffer
  if (fCEPEvent) {
    fCEPEvent->Reset();
    delete fCEPEvent;
    fCEPEvent = 0x0;
  }
  
	// delete array of tracks
  if (fTracks) {
    fTracks->SetOwner(kTRUE);
    fTracks->Clear();
		delete fTracks;
		fTracks = 0x0;
	}

	// delete array of TrackStatus
  if (fTrackStatus) {
    fTrackStatus->Reset();
		delete fTrackStatus;
		fTrackStatus = 0x0;
	}

  // delete fTrigger
  if (fTrigger) {
    delete fTrigger;
    fTrigger = 0x0;
  }
  
  // delete PhysicsSelection
  if (fPhysicsSelection) {
    delete fPhysicsSelection;
    fPhysicsSelection = 0x0;
  }
  
  // delete fEventCuts
  if (fEventCuts) {
    delete fEventCuts;
    fEventCuts = 0x0;
  }

  // delete fTrackCuts
  if (fTrackCuts) {
    delete fTrackCuts;
    fTrackCuts = 0x0;
  }

  // delete fMartinSel
  if (fMartinSel) {
    delete fMartinSel;
    fMartinSel = 0x0;
  }

  // delet PID objects
  if (fPIDCombined1) {
    delete fPIDCombined1;
    fPIDCombined1 = 0x0;
  }
  if (fPIDCombined2) {
    delete fPIDCombined2;
    fPIDCombined2 = 0x0;
  }
  
  // delete lists of QA histograms
  if (flQArnum) {
    delete flQArnum;
    flQArnum = 0x0;
  }
  if (flBBFlag) {
    delete flBBFlag;
    flBBFlag = 0x0;
  }
  if (flSPDpileup) {
    delete flSPDpileup;
    flSPDpileup = 0x0;
  }
  if (flnClunTra) {
    delete flnClunTra;
    flnClunTra = 0x0;
  }
  if (flVtx) {
    delete flVtx;
    flVtx = 0x0;
  }
  if (flV0) {
    delete flV0;
    flV0 = 0x0;
  }
  if (flFMD) {
    delete flFMD;
    flFMD = 0x0;
  }

  // delete fHist and fCEPtree
  if (fHist) {
    delete fHist;
    fHist = 0x0;
  }
  if (fCEPtree) {
    delete fCEPtree;
    fCEPtree = 0x0;
  }

  // CEPUtils
  if (fCEPUtil) {
    delete fCEPUtil;
    fCEPUtil = 0x0;
  }
  
}


//------------------------------------------------------------------------------
void AliAnalysisTaskCEP::UserCreateOutputObjects()
{

  // set the log level
  AliLog::SetGlobalLogLevel(AliLog::kFatal);

  // initialisize some variables
  // fCEPEvent
  fCEPEvent = new CEPEventBuffer();
  // fTracks
  fTracks = new TObjArray();
  // fTrackStatus
  fTrackStatus = new TArrayI();

  // fPhysicsSelection
  //fPhysicsSelection = new AliPhysicsSelection();
  fPhysicsSelection = static_cast<AliPhysicsSelection*> (fInputHandler->GetEventSelection());
  
  // fTrigger
  fTrigger = new AliTriggerAnalysis();
 	fTrigger->SetDoFMD(kTRUE);
	fTrigger->SetFMDThreshold(0.3,0.5);
  fTrigger->ApplyPileupCuts(kTRUE);
 
	// fEventCuts
  fEventCuts = new AliEventCuts();
	
  // fTrackCuts
  Bool_t isRun1 = fCEPUtil->checkstatus(fAnalysisStatus,
    AliCEPBase::kBitisRun1,AliCEPBase::kBitisRun1);
  Int_t clusterCut = 1;                 // 1 = b and c, 0 = d and e
  fTrackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2010(kTRUE,clusterCut);
  if (fCEPUtil->checkstatus(fAnalysisStatus,
    AliCEPBase::kBitTrackCutStudy,AliCEPBase::kBitTrackCutStudy)) {
    printf("Switching TrackCuts Histograms on ...\n");
    fTrackCuts->DefineHistograms();
    fHist->Add(fTrackCuts);
  }
  
  // Martin's selection
	fMartinSel = new AliMultiplicitySelectionCP();
	fMartinSel->InitDefaultTrackCuts(1);          // 1 = b and c, 0 = d and e
  
  // AliCEPUtils
  fCEPUtil = new AliCEPUtils();
  fCEPUtil->SetTPCnclsS(3);             // limit for number of shared clusters
  fCEPUtil->SetTrackDCA(500);           // limit for DCA
  fCEPUtil->SetTrackDCAz(6);            // limit for DCAz
  fCEPUtil->SetTrackEtaRange(-0.9,0.9); // accepted eta range
  
  fCEPUtil->InitTrackCuts(isRun1,clusterCut);  

	// prepare PID Combined
  // three types of Bayes PID
  // 1. with priors
  // 2. without priors
  // 3. with CEP priors
  
  // 1. with priors
	fPIDCombined1 = new AliPIDCombined;
  fPIDCombined1->SetSelectedSpecies(AliPID::kSPECIES);  // This is default
	fPIDCombined1->SetEnablePriors(kTRUE);
	fPIDCombined1->SetDefaultTPCPriors();
  
  // 2. without priors
	fPIDCombined2 = new AliPIDCombined;
  fPIDCombined2->SetSelectedSpecies(AliPID::kSPECIES);  // This is default
	fPIDCombined2->SetEnablePriors(kFALSE);               // priors are set to 1
  
  // 3. set CEP specific priors
  //fPIDCombined3 = new AliPIDCombined;
  //fPIDCombined3->SetSelectedSpecies(AliPID::kSPECIES);  //This is default
	//fPIDCombined3->SetEnablePriors(kTRUE);               // priors are set to 1
  //TString fnameMyPriors = TString("/home/pbuehler/physics/projects/alice/CEP/working/forpass4/res/20160501/MergedPriors.root");
  //TH1F *priordistr[AliPID::kSPECIES];
  //GetMyPriors(fnameMyPriors,priordistr);
  //for (Int_t ii=0; ii<AliPID::kSPECIES; ii++)
  //  fPIDCombined3->SetPriorDistribution((AliPID::EParticleType)ii,priordistr[ii]);
  
  // define the detctors to use for PID ...
  // ... only TPC
  //UInt_t Maskin =
  //UInt_t Maskin = AliPIDResponse::kDetTPC;

  // ... TPC and TOF
  UInt_t Maskin =
    AliPIDResponse::kDetTPC |
    AliPIDResponse::kDetTOF;

  // ... TPC, TOF, and ITS
  //UInt_t Maskin =
  //AliPIDResponse::kDetTPC |
  //  AliPIDResponse::kDetTOF |
  //  AliPIDResponse::kDetITS;

  // ... TPC, TOF, ITS, and TRD
  //UInt_t Maskin =
  //AliPIDResponse::kDetTPC |
  //  AliPIDResponse::kDetTOF |
  //  AliPIDResponse::kDetITS |
  //  AliPIDResponse::kDetTRD;

  fPIDCombined1->SetDetectorMask(Maskin);
  fPIDCombined2->SetDetectorMask(Maskin);
  //fPIDCombined3->SetDetectorMask(Maskin);


  // CreateOutputObjects
  // histograms for various QA tasks
  // histograms for QA versus rnum study
  if (fCEPUtil->checkstatus(fAnalysisStatus,
    AliCEPBase::kBitQArnumStudy,AliCEPBase::kBitQArnumStudy)) {
    
    // get list of histograms
    flQArnum = new TList();
    flQArnum = fCEPUtil->GetQArnumHists(frnummin,frnummax);
    
    // add histograms to the output list
    for (Int_t ii=0; ii<flQArnum->GetEntries(); ii++)
      fHist->Add((TObject*)flQArnum->At(ii));
  }
  
  // histograms BBFlag study
  if (fCEPUtil->checkstatus(fAnalysisStatus,
    AliCEPBase::kBitBBFlagStudy,AliCEPBase::kBitBBFlagStudy)) {
    
    // get list of histograms
    flBBFlag = new TList();
    flBBFlag = fCEPUtil->GetBBFlagQAHists();
    
    // add histograms to the output list
    for (Int_t ii=0; ii<flBBFlag->GetEntries(); ii++)
      fHist->Add((TObject*)flBBFlag->At(ii));
  }
  
  // histograms for SPD pile-up study
  if (fCEPUtil->checkstatus(fAnalysisStatus,
    AliCEPBase::kBitSPDPileupStudy,AliCEPBase::kBitSPDPileupStudy)) {
    
    // get list of histograms
    flSPDpileup = new TList();
    flSPDpileup = fCEPUtil->GetSPDPileupQAHists();
    
    // add histograms to the output list
    for (Int_t ii=0; ii<flSPDpileup->GetEntries(); ii++)
      fHist->Add((TObject*)flSPDpileup->At(ii));
  }
  
  // histograms for nClunTra BG rejection study
  if (fCEPUtil->checkstatus(fAnalysisStatus,
    AliCEPBase::kBitnClunTraStudy,AliCEPBase::kBitnClunTraStudy)) {
    
    // get list of histograms
    flnClunTra = new TList();
    flnClunTra = fCEPUtil->GetnClunTraQAHists();
    
    // add histograms to the output list
    for (Int_t ii=0; ii<flnClunTra->GetEntries(); ii++)
      fHist->Add((TObject*)flnClunTra->At(ii));
  }
  
  // histograms for vertex selection study
  if (fCEPUtil->checkstatus(fAnalysisStatus,
    AliCEPBase::kBitVtxStudy,AliCEPBase::kBitVtxStudy)) {
    
    // get list of histograms
    flVtx = new TList();
    flVtx = fCEPUtil->GetVtxQAHists();
    
    // add histograms to the output list
    for (Int_t ii=0; ii<flVtx->GetEntries(); ii++)
      fHist->Add((TObject*)flVtx->At(ii));
  } else flVtx = NULL;
  
  // histograms for Armenteros-Podolanski plot
  if (fCEPUtil->checkstatus(fAnalysisStatus,
    AliCEPBase::kBitV0Study,AliCEPBase::kBitV0Study)) {
            
    // get list of histograms
    flV0 = new TList();
    flV0 = fCEPUtil->GetV0QAHists();
    
    // add histograms to the output list
    for (Int_t ii=0; ii<flV0->GetEntries(); ii++)
      fHist->Add((TObject*)flV0->At(ii));
  } else flV0 = NULL;
  
  // histograms for FMD study
  if (fCEPUtil->checkstatus(fAnalysisStatus,
    AliCEPBase::kBitFMDStudy,AliCEPBase::kBitFMDStudy)) {
            
    // get list of histograms
    flFMD = new TList();
    flFMD = fCEPUtil->GetFMDQAHists();
    
    // add histograms to the output list
    for (Int_t ii=0; ii<flFMD->GetEntries(); ii++)
      fHist->Add((TObject*)flFMD->At(ii));
  } else flFMD = NULL;
  
  // histogram for event statistics
  fhStatsFlow = fCEPUtil->GetHistStatsFlow();
  fHist->Add(fhStatsFlow);
  
  // add QA histograms for event cut
  fEventCuts->AddQAplotsToList(fHist);

	//CEP tree
	fCEPtree = new TTree("CEP", "CEP");

	// add MonteCarlo truth
  fCEPtree->Branch("MCtruth", &fMCCEPSystem);

  // add branch with CEPEventBuffer
  Int_t split = 2;         // branches are split
  Int_t bsize = 16000; 
  fCEPtree->Branch("CEPEvents","CEPEventBuffer",&fCEPEvent, bsize, split);
  
	PostOutputs();
  
}


//------------------------------------------------------------------------------
void AliAnalysisTaskCEP::UserExec(Option_t *)
{

  // update stats flow
  // events running through UserExec - thus all analyzed events
  fhStatsFlow->Fill(AliCEPBase::kBinTotalInput);

  // get the next event, check ..
  // availability
  // PID
  // magnetic field
  // update the run number fRun
  // get MC event fMCEvent
  if ( !CheckInput() ) {
		PostOutputs();
		return;
	}
  fhStatsFlow->Fill(AliCEPBase::kBinGoodInput);
  ((TH1F*)flQArnum->At(0))->Fill(fRun);
  if (fMCEvent) fhStatsFlow->Fill(AliCEPBase::kBinMCEvent);
  
  AliAODEvent* aodEvent = AODEvent();
  if (aodEvent) printf("Got AODEvent ...\n");
  
    // get event characteristics like ...
  // run number
  fRun = fEvent->GetRunNumber();
  // event number
  Int_t fEventNum = fEvent->GetEventNumberInFile();
  // period
  UInt_t fPeriod = fEvent->GetPeriodNumber();
  // orbit number
  UInt_t fOrbit = fEvent->GetOrbitNumber();
  // bunch cross number
  UShort_t fBunchCross = fEvent->GetBunchCrossNumber();
  // magnetic field strenght
  Double_t fMagField = fEvent->GetMagneticField();
  // number of tracks
  Int_t fNumTracks = fEvent->GetNumberOfTracks();
  // Number of V0s
  Int_t fNumV0s = fEvent->GetNumberOfV0s();
  
  //printf("Run number:         %i\n",fRun);
  //printf("Event number:       %i\n",fEvent);
  //printf("Period number:      %i\n",fPeriod);
  //printf("Orbit number:       %i\n",fOrbit);
  //printf("Bunch cross number: %i\n",fBunchCross);
  //printf("LHC period:         %s\n",fLHCPeriod.Data());
    
  // collission type
  // A:  CINT1A-ABCE-NOPF-ALL         : beam from A-side
  // C:  CINT1C-ABCE-NOPF-ALL         : beam from C-side
  // E:  CINT1-E-NOPF-ALL, CDG5-E     : no beam
  // I:  CINT1B-ABCE-NOPF-ALL, CDG5-I : beam from both sides
  // AC: CDG5-AC                      : beam from one side, side not specified
  Int_t fEventType = fCEPUtil->GetEventType(fEvent);
  //printf("<I - UserExec> EventType: %i\n",fEventType);
  Bool_t isPhyEvent =
    fEvent->GetEventType() == AliRawEventHeaderBase::kPhysicsEvent;
  if (isPhyEvent) fhStatsFlow->Fill(AliCEPBase::kBinPhysEvent);
  
  // does event pass the PhysicsSelection?
  Bool_t isPhys = kFALSE;
  if (fPhysicsSelection) {
    if (fMCEvent) fPhysicsSelection->SetAnalyzeMC(kTRUE);
    isPhys = fPhysicsSelection->IsCollisionCandidate(fEvent)>0;
  }
  if (isPhys) fhStatsFlow->Fill(AliCEPBase::kBinPhysel);

  // pileup and cluster cut
  // see recommendations at
  // https://twiki.cern.ch/twiki/bin/view/ALICE/PWGPPEvSelRun2pp
  Bool_t isPileup = fEvent->IsPileupFromSPD(3,0.8,3.,2.,5.);
  if (isPileup) fhStatsFlow->Fill(AliCEPBase::kBinPileup);
  
  if (flSPDpileup) {
    fCEPUtil->SPDVtxAnalysis(fEvent,3,0.8,3.,2.,5.,flSPDpileup);
  }
  
  AliAnalysisUtils fAnalysisUtils;
  //fAnalysisUtils.SetBSPDCvsTCut(4);
	//fAnalysisUtils.SetASPDCvsTCut(65);
  Bool_t isClusterCut = !fAnalysisUtils.IsSPDClusterVsTrackletBG(fEvent);
  if (isClusterCut) fhStatsFlow->Fill(AliCEPBase::kBinClusterCut);

  if (flnClunTra) {
    fCEPUtil->SPDClusterVsTrackletBGAnalysis(fEvent,flnClunTra);
  }
     
  // EventCuts only works with run2 data
  // this is doing a lot more than the physics selection (by default uses kAny)
  // run number > 225000
  Bool_t isEventCutsel = kFALSE;
  if (fRun >= 225000 && fRun <= 260187) {
    isEventCutsel = fEventCuts->AcceptEvent(fEvent);
  }
  if (isEventCutsel) fhStatsFlow->Fill(AliCEPBase::kBinEventCut);

  // vertex type and position
	// first check vertex from tracks
  // if that does not exists, then try vertex from SPD
  // VertexType:
  // -1: no vertex
  //  1: from SPD
  //  2: from tracks
  Int_t kVertexType = fCEPUtil->GetVtxPos(fEvent,&fVtxPos);
  if (kVertexType!=AliCEPBase::kVtxUnknown)
    fhStatsFlow->Fill(AliCEPBase::kBinVtx);

  if (flVtx) {
    fCEPUtil->VtxAnalysis(fEvent,flVtx);
  }
     
  // did the double-gap trigger (CCUP13-B-SPD1-CENTNOTRD) fire?
  // this is relevant for the LHC16[k,l,o,p] data
  // in case of MC data and data containing no DG trigger
  // the trigger needs to be replaied
  // different triggers are considered
  // CINT11-B-NOPF-CENTNOTRD, DG trigger has to be replaied, LHC16[d,e]
  // CCUP2-B-SPD1-CENTNOTRD, DG trigger has to be replaied, LHC16[h,i,j]
  // CCUP13-B-SPD1-CENTNOTRD = DG trigger, LHC16[k,l,o,p]
  TString firedTriggerClasses = fEvent->GetFiredTriggerClasses();
  Bool_t isReplay = fMCEvent
    || firedTriggerClasses.Contains("CINT11-B-NOPF-CENTNOTRD")
    || firedTriggerClasses.Contains("CCUP2-B-SPD1-CENTNOTRD");
  
  // The following is needed to replay the DG trigger
  // The STG trigger in 2016 required two online tracklets without additional
  // topology.
  TBits foMap = fInputEvent->GetMultiplicity()->GetFastOrFiredChips();
  Bool_t isSTGtriggerFired = IsSTGFired(&foMap);
    
  Bool_t isDGTrigger = kFALSE;
  if (isReplay) {
    
    // this part of the code was proposed by Evgeny Kryshen
    // to replay the DG trigger in MC data
    Bool_t isV0Afired=0;
    Bool_t isV0Cfired=0;
    AliVVZERO* esdV0 = fEvent->GetVZEROData();
    for (Int_t i=32; i<64; i++) isV0Afired |= esdV0->GetBBFlag(i);
    for (Int_t i=0; i<32; i++)  isV0Cfired |= esdV0->GetBBFlag(i);

    isDGTrigger = isSTGtriggerFired && !isV0Afired && !isV0Cfired;
    
    // printf("DG trigger replaied: %i %i %i\n",
    //  isSTGtriggerFired,!isV0Afired,!isV0Cfired);
  
  } else {
  
    if (firedTriggerClasses.Contains("CCUP13-B-SPD1-CENTNOTRD")) {
      isDGTrigger = kTRUE;
      
      fhStatsFlow->Fill(AliCEPBase::kBinDGTrigger);
      ((TH1F*)flQArnum->At(1))->Fill(fRun);

      // past-future trigger protection
      if (flBBFlag)
        fCEPUtil->BBFlagAnalysis(fEvent,flBBFlag);
    }
  }
  
  
  // number of tracklets
  const AliVMultiplicity *mult = fEvent->GetMultiplicity();
  Int_t nTracklets = mult->GetNumberOfTracklets();
  
  // get trigger information using AliTriggerAnalysis.IsOfflineTriggerFired
  // kSPDGFOBits: SPD (any fired chip)
  // kV0A, kV0C: V0
  // kMB1: SPD || V0A || V0B
  // kMB2 = kV0AND: SPD && (V0A || V0B)
  // kFMAA, kFMDC: FMD
  // kADA, kADC: AD
  // kZDC: hit on any side of ZDC
  // kZNA, kZNC: ZN
  // printf("<I - UserExec> Doing trigger analysis ...\n");
	Bool_t isSPD  =
    (fTrigger->IsOfflineTriggerFired(fEvent,AliTriggerAnalysis::kSPDGFO));
	// Bool_t isSPD  =
  //   (fTrigger->IsOfflineTriggerFired(fEvent,AliTriggerAnalysis::kSPDGFOBits));
	Bool_t isMBOR  =
    (fTrigger->IsOfflineTriggerFired(fEvent,AliTriggerAnalysis::kMB1));
	// Bool_t isMBAND  =
  //   (fTrigger->IsOfflineTriggerFired(fEvent,AliTriggerAnalysis::kMB2));
  Bool_t isMBAND =
    (fTrigger->IsOfflineTriggerFired(fEvent,AliTriggerAnalysis::kV0AND));
  Bool_t isV0A =
    (fTrigger->IsOfflineTriggerFired(fEvent,AliTriggerAnalysis::kV0A));
  Bool_t isV0C =
    (fTrigger->IsOfflineTriggerFired(fEvent,AliTriggerAnalysis::kV0C));
  Bool_t isADA =
    (fTrigger->IsOfflineTriggerFired(fEvent,AliTriggerAnalysis::kADA));
  Bool_t isADC =
    (fTrigger->IsOfflineTriggerFired(fEvent,AliTriggerAnalysis::kADC));
  Bool_t isFMDA =
    (fTrigger->IsOfflineTriggerFired(fEvent,AliTriggerAnalysis::kFMDA));
  Bool_t isFMDC =
    (fTrigger->IsOfflineTriggerFired(fEvent,AliTriggerAnalysis::kFMDC)); 
	Bool_t isZDC  =
    (fTrigger->IsOfflineTriggerFired(fEvent,AliTriggerAnalysis::kZDC));
	Bool_t isZDNA  =
    (fTrigger->IsOfflineTriggerFired(fEvent,AliTriggerAnalysis::kZNA));
	Bool_t isZDNC  =
    (fTrigger->IsOfflineTriggerFired(fEvent,AliTriggerAnalysis::kZNC));
  
  Bool_t isV0DG = isSPD && !(isV0A || isV0C);
  Bool_t isADDG = isSPD && !(isADA || isADC);
  Bool_t isFMDDG = isSPD && !(isFMDA || isFMDC);

  if (isMBOR) fhStatsFlow->Fill(AliCEPBase::kBinMBOR);
  if (isMBAND) fhStatsFlow->Fill(AliCEPBase::kBinMBAND);

  if (isMBOR) ((TH1F*)flQArnum->At(2))->Fill(fRun);
  if (isV0DG) ((TH1F*)flQArnum->At(4))->Fill(fRun);
  if (isADDG) ((TH1F*)flQArnum->At(5))->Fill(fRun);
  if (isFMDDG)((TH1F*)flQArnum->At(6))->Fill(fRun);
  
  
  // compre isSPD and isSTGtriggerFired
  // printf("SPD fired: %i %i\n",isSPD,isSTGtriggerFired);
  
  // determine the gap condition using
  // ITS,V0,FMD,AD,ZD
  // printf("<I - UserExec> Creating GapCondition flag ...\n");
  Int_t fEventCondition = AliCEPBase::kETBaseLine
    + isEventCutsel * AliCEPBase::kETEventCut
    + isPhys * AliCEPBase::kETPhyssel
    + isPileup * AliCEPBase::kETPileup
    + isClusterCut * AliCEPBase::kETClusterCut
    + isDGTrigger * AliCEPBase::kETDGTrigger
    + (kVertexType!=AliCEPBase::kVtxUnknown) * AliCEPBase::kETVtx
    + isMBOR * AliCEPBase::kETMBOR + isMBAND * AliCEPBase::kETMBAND
    + isSPD  * AliCEPBase::kETSPDA
    + isV0A  * AliCEPBase::kETV0A  + isV0C   * AliCEPBase::kETV0C
    + isFMDA * AliCEPBase::kETFMDA + isFMDC  * AliCEPBase::kETFMDC
    + isADA  * AliCEPBase::kETADA  + isADC   * AliCEPBase::kETADC
    + isZDC  * AliCEPBase::kETZDCA
    + isZDNA * AliCEPBase::kETZDNA + isZDNC  * AliCEPBase::kETZDNC;
  
  // update fhStatsFlow
  if (!isV0A  && !isV0C)  fhStatsFlow->Fill(AliCEPBase::kBinnoV0);
  if (!isFMDA && !isFMDC) fhStatsFlow->Fill(AliCEPBase::kBinnoFMD);
  if (!isADA  && !isADC)  fhStatsFlow->Fill(AliCEPBase::kBinnoAD);

  // count the number of tracks
  // apply standard cuts to tracks and count accepted tracks
  // and fill the trackcut QA histograms
  if (fCEPUtil->checkstatus(fAnalysisStatus,
    AliCEPBase::kBitTrackCutStudy,AliCEPBase::kBitTrackCutStudy)) {
    Int_t nTracksTotal = fESDEvent->GetNumberOfTracks();
    Int_t nTracksTPCITS = 0;
    for (Int_t ii = 0; ii < nTracksTotal; ii++) {
      AliESDtrack *track = fESDEvent->GetTrack(ii);
	  	if (track) if (fTrackCuts->AcceptTrack(track)) nTracksTPCITS++;
	  }
  }
  
  // get an TObjArray of tracks with the corresponding TrackStatus
  // The TrackStatus is contained in an array of UInt_t
  // The definition of the TrackStatus bits is given in AliCEPBase.h
  Int_t nTracks = fCEPUtil->AnalyzeTracks(fESDEvent,fTracks,fTrackStatus);
  // printf("Total number of charged tracks: %i\n",nTracks);
  
  // V0 study
  if (flV0)  fCEPUtil->V0Analysis(fESDEvent,flV0);
  if (flFMD && !isPileup) fCEPUtil->FMDAnalysis(fESDEvent,fTrigger,flFMD);

  // nbad track with !kTTTPCScluster
  UInt_t mask = AliCEPBase::kTTTPCScluster;
  UInt_t pattern = 0;
  Int_t nbad = fCEPUtil->countstatus(fTrackStatus,mask,pattern);
  if (nbad==0) fhStatsFlow->Fill(AliCEPBase::kBinSharedCluster);
  fEventCondition += (nbad==0) * AliCEPBase::kETSClusterCut;
  
	// Martin's selection
	TArrayI Mindices;
	Int_t nMartinSel = fMartinSel->GetNumberOfITSTPCtracks(fESDEvent,Mindices);
  //printf("<I - UserExec> nMartinSel: %i\n",nMartinSel);
    
  // use the AliCEPUtils::GetCEPTracks method
  TArrayI *Pindices  = new TArrayI();
  Int_t nPaulSel = fCEPUtil->GetCEPTracks(fESDEvent,fTrackStatus,Pindices);
  // if (nMartinSel>0 || nPaulSel>0) printf("%i/%i good CEP tracks\n", nPaulSel,nMartinSel);

  // get the tracks which meet the TT conditions
  // the TT mask and pattern are given as input parameter
  // if TTmask=TTpattern=kTTBaseLine, then all tracks are selected
  TArrayI *TTindices  = new TArrayI();
  Int_t nTracksTT = fCEPUtil->countstatus(fTrackStatus,
    fTTmask, fTTpattern, TTindices);
  // printf("Number of selected tracks: %i\n",nTracksTT);

  // get tracks which pass default TPCITS cuts
  //TArrayI *TPCITSindices  = new TArrayI();
  //Int_t nTracksTPCITS = fCEPUtil->countstatus(fTrackStatus,
  //  AliCEPBase::kTTAccITSTPC, AliCEPBase::kTTAccITSTPC, TPCITSindices);
  
  // for test purposes -------------------------------------------------------
  if (kFALSE) {
    
    // now use the fTrackStatus to scrutinize the event and tracks
    // e.g. Martin's selection
    // 1. no track with !kTTTPCScluster
    // 2. kTTITSpure -> npureITSTracks
    // 3. kTTDCA && !kTTV0 && !kTTITSpure && kTTZv -> nTrackSel
    // 4. 3. && kTTeta && (kTTAccITSTPC || kTTAccITSSA) -> nTrackAccept
    // 5. nTrackSel>=npureITSTracks && nTrackSel>=nTracklets && nTrackAccept==nTrackSel
    // 6. pass fpassedFiredChipsTest
    //  
    // if all criteria are met the event is accepted and the number of
    // good tracks is nTrackSel = nMartinSel
    TArrayI *masks    = new TArrayI();
    TArrayI *patterns = new TArrayI();
    TArrayI *indices  = new TArrayI();

    // the number of ITSpure tracks
    mask = AliCEPBase::kTTITSpure;
    pattern = AliCEPBase::kTTITSpure;
    Int_t npureITSTracks = fCEPUtil->countstatus(fTrackStatus,mask,pattern);
      
    // the number of tracks which fulfill 
    mask = AliCEPBase::kTTDCA+AliCEPBase::kTTV0+AliCEPBase::kTTITSpure+
      AliCEPBase::kTTZv;
    pattern =  AliCEPBase::kTTDCA+AliCEPBase::kTTZv;
    Int_t nTrackSel = fCEPUtil->countstatus(fTrackStatus,mask,pattern,indices);
  
    // mask += AliCEPBase::kTTeta+AliCEPBase::kTTAccITSTPC+AliCEPBase::kTTAccITSSA+
    //   AliCEPBase::kTTFiredChips;
    // pattern += AliCEPBase::kTTeta+AliCEPBase::kTTAccITSTPC+AliCEPBase::kTTAccITSSA+
    //   AliCEPBase::kTTFiredChips;
  
    Int_t nTrackAccept;
    mask += AliCEPBase::kTTeta;
    pattern += AliCEPBase::kTTeta;
    //nTrackAccept = fCEPUtil->countstatus(fTrackStatus,mask,pattern);
    masks->Set(2);
    patterns->Set(2);
    masks->AddAt(mask+AliCEPBase::kTTAccITSTPC,0);
    masks->AddAt(mask+AliCEPBase::kTTAccITSSA, 1);
    patterns->AddAt(pattern+AliCEPBase::kTTAccITSTPC,0);
    patterns->AddAt(pattern+AliCEPBase::kTTAccITSSA, 1);
    nTrackAccept = fCEPUtil->countstatus(fTrackStatus,masks,patterns);
  
    mask = AliCEPBase::kTTFiredChips;
    pattern = AliCEPBase::kTTFiredChips;
    Int_t nTracksFiredChips = fCEPUtil->countstatus(fTrackStatus,mask,pattern);
    Bool_t fpassedFiredChipsTest =
      fCEPUtil->TestFiredChips(fESDEvent,indices);

    //printf("<I - UserExec> nTracksTotal    : %i\n",nTracksTotal);
    printf("\n");
    printf("<I - UserExec> nTracklets      : %i\n",nTracklets);
    printf("<I - UserExec> npureITSTracks  : %i\n",npureITSTracks);
    printf("<I - UserExec> nTracks         : %i\n",nTracks);
    printf("<I - UserExec> nTracksTT       : %i\n",nTracksTT);
    printf("<I - UserExec> nBad            : %i\n",nbad);
    printf("<I - UserExec> nTrackSel       : %i\n",nTrackSel);
    printf("<I - UserExec> nTrackFiredChips: %i - event passed test: %i\n",nTracksFiredChips,fpassedFiredChipsTest);
    printf("<I - UserExec> nTrackAccept    : %i\n",nTrackAccept);
    //printf("<I - UserExec> nTrackAccept    : %i %i %i\n",nTrackAccept,nMartinSel,nPaulSel);

    // clean up
    if (masks) {
      delete masks;
      masks = 0x0;
    }
    if (patterns) {
      delete patterns;
      patterns = 0x0;
    }
    if (indices) {
      delete indices;
      indices = 0x0;
    }
  
  }
  // for test purposes -------------------------------------------------------
  
  // decide which events have to be saved
  // if the number of considered tracks (nTracksTT) is within the limits ...
  Bool_t isToSave = (0 < nTracksTT) && (nTracksTT <= fnumTracksMax);

  // consider in addition that only a fraction (ffracDG, ffracNDG) of the
  // accepted events is saved
  TRandom3 rnd(0);
  
  // AND DG ET conditions are met ...
  Bool_t isToSaveDG =
    fCEPUtil->checkstatus(fEventCondition,fETmaskDG,fETpatternDG);
  if (rnd.Rndm() < (1. - ffracDG)) isToSaveDG = kFALSE;
    
  // AND NDG ET conditions are met ...
  Bool_t isToSaveNDG =
    fCEPUtil->checkstatus(fEventCondition,fETmaskNDG,fETpatternNDG);
  if (rnd.Rndm() < (1. - ffracNDG)) isToSaveNDG = kFALSE;
    
  // OR kBitSaveAllEvents is set
  isToSave = ( isToSave && (isToSaveDG || isToSaveNDG) ) ||
    fCEPUtil->checkstatus(fAnalysisStatus,
    AliCEPBase::kBitSaveAllEvents,AliCEPBase::kBitSaveAllEvents);
  
  if ( isToSave ) {
    
    // update fhStatsFlow
    ((TH1F*)flQArnum->At(3))->Fill(fRun);
    fhStatsFlow->Fill(AliCEPBase::kBinSaved);
    if (isToSaveDG) {
      ((TH1F*)flQArnum->At(7))->Fill(fRun);
      fhStatsFlow->Fill(AliCEPBase::kBinDG);
    }
    if (isToSaveNDG) {
      ((TH1F*)flQArnum->At(8))->Fill(fRun);
      fhStatsFlow->Fill(AliCEPBase::kBinNDG);
    }
    
    // fill the CEPEventBuffer
    fCEPEvent->Reset();
    
    // general event informations
    fCEPEvent->SetRunNumber(fRun);
    fCEPEvent->SetEventNumber(fEventNum);
    fCEPEvent->SetPeriodNumber(fPeriod);
    fCEPEvent->SetOrbitNumber(fOrbit);
    fCEPEvent->SetBunchCrossNumber(fBunchCross);
    fCEPEvent->SetCollissionType(fEventType);
    fCEPEvent->SetMagnField(fMagField);
    fCEPEvent->SetFiredTriggerClasses(firedTriggerClasses);
    
    // FP Flags of V0 and AD
    fCEPEvent->SetPFFlags(fEvent);
  
    // set event status word
    fCEPEvent->SetEventCondition(fEventCondition);

    // number of tracklets and residuals, vertex
    fCEPEvent->SetnTracksTotal(nTracks);
    fCEPEvent->SetnTracklets(nTracklets);
    fCEPEvent->SetnResiduals(fCEPUtil->GetResiduals(fESDEvent));
    fCEPEvent->SetnMSelection(nMartinSel);
    fCEPEvent->SetnV0(fNumV0s);
    fCEPEvent->SetVtxType(kVertexType);
    fCEPEvent->SetVtxPos(fVtxPos);
  
    // if this is a MC event then get the MC true information
    // and save it into the event buffer
    AliStack *stack;
    if (fMCEvent) {
	    
      // MC generator and process type
      TString fMCGenerator;
	    Int_t fMCProcess; 
      fCEPUtil->DetermineMCprocessType(fMCEvent,fMCGenerator,fMCProcess);
      
      // get MC vertex
      TParticle* prot1 = NULL;
      stack = fMCEvent->Stack();
      if (stack) {
        Int_t nPrimaries = stack->GetNprimary();
        //printf("number of tracks: primaries - %i, reconstructed - %i\n",
        //  nPrimaries,nTracks);
        prot1 = stack->Particle(0);
      }

      // update the event buffer
      fCEPEvent->SetMCGenerator(fMCGenerator);
      fCEPEvent->SetMCProcessType(fMCProcess);
      fCEPEvent->SetMCVtxPos(prot1->Vx(),prot1->Vy(),prot1->Vz());
    }
    
    // add tracks
    // all tracks are added, independent of the track status
    Double_t mom[3];
    Double_t stat,nsig,probs[AliPID::kSPECIES];
    for (Int_t ii=0; ii<nTracksTT; ii++) {
    
      // proper poiter into fTracks and fTrackStatus
      Int_t trkIndex = TTindices->At(ii);
      
      // the original track
      AliESDtrack *tmptrk = (AliESDtrack*) fTracks->At(trkIndex);

      // create new CEPTrackBuffer and fill it
      CEPTrackBuffer *trk = new CEPTrackBuffer();
      
      trk->SetTrackStatus(fTrackStatus->At(trkIndex));
      trk->SetTOFBunchCrossing(tmptrk->GetTOFBunchCrossing());
      trk->SetChargeSign((Int_t)tmptrk->Charge());
      trk->SetITSncls(tmptrk->GetNumberOfITSClusters());
      trk->SetTPCncls(tmptrk->GetNumberOfTPCClusters());
      trk->SetTRDncls(tmptrk->GetNumberOfTRDClusters());
      trk->SetTPCnclsS(tmptrk->GetTPCnclsS());
      trk->SetZv(tmptrk->Zv());
      tmptrk->GetPxPyPz(mom);
      trk->SetMomentum(TVector3(mom));
      
      // set PID information
      // ... TPC
      stat = fPIDResponse->ComputePIDProbability(AliPIDResponse::kTPC,tmptrk,AliPID::kSPECIES,probs);
      trk->SetPIDTPCStatus(stat);
      trk->SetPIDTPCSignal(tmptrk->GetTPCsignal());
      for (Int_t jj=0; jj<AliPID::kSPECIES; jj++) {
        stat = fPIDResponse->NumberOfSigmas(
          AliPIDResponse::kTPC,tmptrk,(AliPID::EParticleType)jj,nsig);
        trk->SetPIDTPCnSigma(jj,nsig);
        trk->SetPIDTPCProbability(jj,probs[jj]);
      }
      
      // ... TOF
      stat = fPIDResponse->ComputePIDProbability(AliPIDResponse::kTOF,tmptrk,AliPID::kSPECIES,probs);
      trk->SetPIDTOFStatus(stat);
      trk->SetPIDTOFSignal(tmptrk->GetTOFsignal());
      for (Int_t jj=0; jj<AliPID::kSPECIES; jj++) {
        stat = fPIDResponse->NumberOfSigmas(
          AliPIDResponse::kTOF,tmptrk,(AliPID::EParticleType)jj,nsig);
        trk->SetPIDTOFnSigma(jj,nsig);
        trk->SetPIDTOFProbability(jj,probs[jj]);
      }
      
      // ... Bayes
      stat = fPIDCombined1->ComputeProbabilities(tmptrk, fPIDResponse, probs);
      trk->SetPIDBayesStatus(stat);
      for (Int_t jj=0; jj<AliPID::kSPECIES; jj++)
        trk->SetPIDBayesProbability(jj,probs[jj]);
      
      // get MC truth
      Int_t MCind = tmptrk->GetLabel();
      if (fMCEvent && MCind >= 0) {
        
        TParticle* part = stack->Particle(MCind);
        // printf("MC particle (%i): %f/%f/%f - %f/%f\n",
        //  ii,part->Px(),part->Py(),part->Pz(),part->GetMass(),part->Energy());
        
        // set MC mass and momentum
        TLorentzVector lv;
        part->Momentum(lv);
        
        trk->SetMCPID(part->GetPdgCode());
        trk->SetMCMass(part->GetMass());
        trk->SetMCMomentum(TVector3(lv.Px(),lv.Py(),lv.Pz()));

      }

      // add track to the CEPEventBuffer
      fCEPEvent->AddTrack(trk);
    
    }
    
    // save output
    fCEPtree->Fill();
    PostOutputs();
    
  } else {
    PostData(1, fHist);
  }
    
    
  // clean up
  //if (Pindices) {
  //  delete Pindices;
  //  Pindices = 0x0;
  //}
  //if (TPCITSindices) {
  //  delete TPCITSindices;
  //  TPCITSindices = 0x0;
  //}
  if (TTindices) {
    delete TTindices;
    TTindices = 0x0;
  }

}

//------------------------------------------------------------------------------
// privat methods
//
//------------------------------------------------------------------------------
Bool_t AliAnalysisTaskCEP::CheckInput()
{
  if (!fInputHandler) return kFALSE;

	//general protection
	// are we dealing with an ESD or AOD?
	fEvent = fInputHandler->GetEvent();
  if (fEvent->GetDataLayoutType()==AliVEvent::kESD) fisESD = kTRUE;
  if (fEvent->GetDataLayoutType()==AliVEvent::kAOD) fisAOD = kTRUE;
  
  if (fisESD) fESDEvent = (AliESDEvent*)fEvent;
  
  if (!fisESD) {
		printf("<E - CheckInput> No valid ESD event!\n");
    return kFALSE;
	}
    
  // check incomplete DAQ, see recommendation at
  // https://twiki.cern.ch/twiki/bin/view/ALICE/PWGPPEvSelRun2pp
  if (fEvent->IsIncompleteDAQ()) {
		printf("<E - CheckInput> Incomplete DAQ!\n");
    return kFALSE;
	}

	// get the LHC period name
  //TList *uiList = fInputHandler->GetUserInfo();
  //AliProdInfo prodInfo(uiList);
  //fLHCPeriod = TString(prodInfo.GetLHCPeriod());

  // PID response
	fPIDResponse = (AliPIDResponse*)fInputHandler->GetPIDResponse();
	if (!fPIDResponse){
		printf("<E - CheckInput> No PID information!\n");
		return kFALSE;
	}

	// check magnetic field
  if (TMath::Abs(fEvent->GetMagneticField())<1) {
		printf("<E - CheckInput> Wrong Bfield value! %f\n",
      fEvent->GetMagneticField());
		return kFALSE;
	}

	// update the run number
  Int_t tmprun = fEvent->GetRunNumber();
	if (fRun != tmprun) {
		fRun = tmprun;
		fCEPUtil->SPDLoadGeom(fRun);
	}

	// get MC event
	fMCEvent = MCEvent();
  if (fMCEvent) {  
    // Bug fix 28.05.2016 - do not trust to presence of MC handler, check if the content is valid
    //                    - proper solution (autodetection of MC information) to be implemented 
    if (fMCEvent->Stack()==NULL) fMCEvent=NULL;
  }
  /*
  if (fMCEvent) {
    printf("Is MC event!\n");
  } else {
    printf("Is NOT MC event!\n");
  }
  */
  
	return kTRUE;
}

//------------------------------------------------------------------------------
void AliAnalysisTaskCEP::PostOutputs()
{
	// PostData
	Int_t iOutputSlot = 1;
	PostData(iOutputSlot++, fHist);
  PostData(iOutputSlot++, fCEPtree);
	
}


//------------------------------------------------------------------------------
// code to check if the STG trigger had fired
// code from Evgeny Kryshen
Bool_t AliAnalysisTaskCEP::IsSTGFired(TBits* fFOmap,Int_t dphiMin,Int_t dphiMax)
{

  Int_t hitcnt = 0;
  Bool_t stg = kFALSE;
  
  if (!fFOmap) {
    // printf("<AliAnalysisTaskCEP::IsSTGFired> Problem with F0map - a!\n");
    return stg;
  }

  Int_t n1 = fFOmap->CountBits(400);
  Int_t n0 = fFOmap->CountBits()-n1;
  if (n0<1 || n1<1) {
    // printf("<AliAnalysisTaskCEP::IsSTGFired> Problem with F0map - b!\n");
    return stg;
  }
  
  
  Bool_t l0[20]={0};
  Bool_t l1[40]={0};
  Bool_t phi[20]={0};
  for (Int_t i=0;   i< 400; ++i) if (fFOmap->TestBitNumber(i)) l0[      i/20] = 1;
  for (Int_t i=400; i<1200; ++i) if (fFOmap->TestBitNumber(i)) l1[(i-400)/20] = 1;
  for (Int_t i=0; i<20; ++i) phi[i] = l0[i] & (l1[(2*i)%40] | l1[(2*i+1)%40] | l1[(2*i+2)%40] | l1[(2*i+39)%40]);
  for (Int_t dphi=dphiMin;dphi<=dphiMax;dphi++) {
    for (Int_t i=0; i<20; ++i) {
      stg |= phi[i] & phi[(i+dphi)%20];
      if (phi[i] & phi[(i+dphi)%20]) hitcnt++;
    }
  }
  // printf("hitcnt: %i\n",hitcnt);

  return stg;

} 
