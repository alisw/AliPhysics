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

#include "AliAnalysisManager.h"
#include "AliProdInfo.h"
#include "AliESDInputHandler.h"
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
  Long_t state, Int_t numTracksMax):
	AliAnalysisTaskSE(name)
  , fnumTracksMax(numTracksMax)
  , fLHCPeriod(TString(""))
  , fRun(-1)
	, fAnalysisStatus(state)
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
	, fhStatsFlow(0x0)
	, fHist(0x0)
	, fCEPtree(0x0)
  , fCEPUtil(0x0)
{

	// slot in TaskSE must start from 1
	Int_t iOutputSlot = 1;
	DefineOutput(iOutputSlot++, TList::Class());
  DefineOutput(iOutputSlot++, TTree::Class());

}


//------------------------------------------------------------------------------
AliAnalysisTaskCEP::AliAnalysisTaskCEP():
	AliAnalysisTaskSE()
  , fnumTracksMax(0)
  , fLHCPeriod(TString(""))
  , fRun(-1)
	, fAnalysisStatus(0x0)
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
	, fhStatsFlow(0x0)
	, fHist(0x0)
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
  fTrigger->ApplyPileupCuts(kTRUE);
  
	// fEventCuts
  fEventCuts = new AliEventCuts();
	// fTrackCuts
  Bool_t isRun1 = fCEPUtil->checkstatus(fAnalysisStatus,
    AliCEPBase::kBitisRun1,AliCEPBase::kBitisRun1);
  Int_t clusterCut = 1;                 // 1 = b and c, 0 = d and e
  fTrackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2010(kTRUE,clusterCut);
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
	// TList for Histograms
	fHist = new TList;
	fHist->SetOwner();    // ensures that the histograms are all deleted on exit!

  fhStatsFlow = fCEPUtil->GetHistStatsFlow();
  fHist->Add(fhStatsFlow);

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
  // printf("<I - UserExec> Printing ESD event ...\n");
  // fESDEvent->Print();
  fhStatsFlow->Fill(AliCEPBase::kBinGoodInput);
  if (fMCEvent) fhStatsFlow->Fill(AliCEPBase::kBinMCEvent);
  
  // get event characteristics like ...
  // run number
  fRun = fESDEvent->GetRunNumber();
  // event number
  Int_t fEvent = fESDEvent->GetEventNumberInFile();
  // period
  UInt_t fPeriod = fESDEvent->GetPeriodNumber();
  // orbit number
  UInt_t fOrbit = fESDEvent->GetOrbitNumber();
  // bunch cross number
  UShort_t fBunchCross = fESDEvent->GetBunchCrossNumber();
  // magnetic field strenght
  Double_t fMagField = fESDEvent->GetMagneticField();
  // number of tracks
  Int_t fNumTracks = fESDEvent->GetNumberOfTracks();
  // Number of V0s
  Int_t fNumV0s = fESDEvent->GetNumberOfV0s();
  
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
  Int_t fEventType = fCEPUtil->GetEventType(fESDEvent);
  //printf("<I - UserExec> EventType: %i\n",fEventType);
  Bool_t isPhyEvent =
    fESDEvent->GetEventType() == AliRawEventHeaderBase::kPhysicsEvent;
  if (isPhyEvent) fhStatsFlow->Fill(AliCEPBase::kBinPhysEvent);
  
  // does event pass the PhysicsSelection?
  if (fMCEvent) fPhysicsSelection->SetAnalyzeMC(kTRUE);
  Bool_t isPhys = kFALSE;
  if (fPhysicsSelection) isPhys = fPhysicsSelection->IsCollisionCandidate(fESDEvent)>0;
  if (isPhys) fhStatsFlow->Fill(AliCEPBase::kBinPhySel);

  // pileup and cluster cut
  // see recommendations at
  // https://twiki.cern.ch/twiki/bin/view/ALICE/PWGPPEvSelRun2pp
  AliAnalysisUtils fAnalysisUtils;
  Bool_t isPileup = fESDEvent->IsPileupFromSPD();
  if (isPileup) fhStatsFlow->Fill(AliCEPBase::kBinPileup);
  
  //fAnalysisUtils.SetBSPDCvsTCut(4);
	//fAnalysisUtils.SetASPDCvsTCut(65);
  Bool_t isClusterCut = !fAnalysisUtils.IsSPDClusterVsTrackletBG(fESDEvent);
  if (isClusterCut) fhStatsFlow->Fill(AliCEPBase::kBinClusterCut);

  // EventCuts only works with run2 data
  // this is doing a lot more than the physics selection (by default uses kAny)
  // run number > 225000
  Bool_t isEventCutsel = kFALSE;
  if (fRun > 225000) {
    isEventCutsel = fEventCuts->AcceptEvent((AliVEvent*)fESDEvent);
  }
  if (isEventCutsel) fhStatsFlow->Fill(AliCEPBase::kBinEventCut);

  // vertex type and position
	// first check vertex from tracks
  // if that does not exists, then try vertex from SPD
  // VertexType:
  // -1: no vertex
  //  1: from SPD
  //  2: from tracks
  Int_t kVertexType = fCEPUtil->GetVtxPos(fESDEvent,&fVtxPos);
  printf("VtxPos 2 %f %f %f \n",fVtxPos.X(),fVtxPos.Y(),fVtxPos.Z());
  
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
    (fTrigger->IsOfflineTriggerFired(fESDEvent,AliTriggerAnalysis::kSPDGFO));
	// Bool_t isSPD  =
  //   (fTrigger->IsOfflineTriggerFired(fESDEvent,AliTriggerAnalysis::kSPDGFOBits));
	Bool_t isMBOR  =
    (fTrigger->IsOfflineTriggerFired(fESDEvent,AliTriggerAnalysis::kMB1));
	// Bool_t isMBAND  =
  //   (fTrigger->IsOfflineTriggerFired(fESDEvent,AliTriggerAnalysis::kMB2));
  Bool_t isMBAND =
    (fTrigger->IsOfflineTriggerFired(fESDEvent,AliTriggerAnalysis::kV0AND));
  Bool_t isV0A =
    (fTrigger->IsOfflineTriggerFired(fESDEvent,AliTriggerAnalysis::kV0A));
  Bool_t isV0C =
    (fTrigger->IsOfflineTriggerFired(fESDEvent,AliTriggerAnalysis::kV0C));
  Bool_t isADA =
    (fTrigger->IsOfflineTriggerFired(fESDEvent,AliTriggerAnalysis::kADA));
  Bool_t isADC =
    (fTrigger->IsOfflineTriggerFired(fESDEvent,AliTriggerAnalysis::kADC));
	fTrigger->SetFMDThreshold(0.3,0.5);
  Bool_t isFMDA =
    (fTrigger->IsOfflineTriggerFired(fESDEvent,AliTriggerAnalysis::kFMDA));
  Bool_t isFMDC =
    (fTrigger->IsOfflineTriggerFired(fESDEvent,AliTriggerAnalysis::kFMDC));
	Bool_t isZDC  =
    (fTrigger->IsOfflineTriggerFired(fESDEvent,AliTriggerAnalysis::kZDC));
	Bool_t isZDNA  =
    (fTrigger->IsOfflineTriggerFired(fESDEvent,AliTriggerAnalysis::kZNA));
	Bool_t isZDNC  =
    (fTrigger->IsOfflineTriggerFired(fESDEvent,AliTriggerAnalysis::kZNC));

  if (isMBOR) fhStatsFlow->Fill(AliCEPBase::kBinMBOR);
  if (isMBAND) fhStatsFlow->Fill(AliCEPBase::kBinMBAND);
  
  // determine the gap condition using
  // ITS,V0,FMD,AD,ZD
  // printf("<I - UserExec> Creating GapCondition flag ...\n");
  Int_t fCurrentGapCondition = AliCEPBase::kBitBaseLine
    + isMBOR * AliCEPBase::kBitMBOR + isMBAND * AliCEPBase::kBitMBAND
    + isSPD  * AliCEPBase::kBitSPDA
    + isV0A  * AliCEPBase::kBitV0A  + isV0C   * AliCEPBase::kBitV0C
    + isFMDA * AliCEPBase::kBitFMDA + isFMDC  * AliCEPBase::kBitFMDC
    + isADA  * AliCEPBase::kBitADA  + isADC   * AliCEPBase::kBitADC
    + isZDC  * AliCEPBase::kBitZDCA
    + isZDNA * AliCEPBase::kBitZDNA + isZDNC  * AliCEPBase::kBitZDNC;
  
	// did the double-gap trigger (CCUP13-B-SPD1-CENTNOTRD) fire?
  // this is relevant for the 2016 data
  // compare with (isSPD  && (!isV0A && !isV0C))
  Bool_t isDGTrigger = kFALSE;
  TString firedTriggerClasses = fESDEvent->GetFiredTriggerClasses();
  if (firedTriggerClasses.Contains("CCUP13-B-SPD1-CENTNOTRD")) {
    isDGTrigger = kTRUE;
  }
  //printf("<I - UserExec> firedTriggerClasses: %s\n",firedTriggerClasses.Data());
  if (isDGTrigger) fhStatsFlow->Fill(AliCEPBase::kBinDGTrigger);
  
  // count the number of tracks
  // apply standard cuts to tracks and count accepted tracks
  Int_t nmbTrk_ST = 0;
  for (Int_t ii = 0; ii < fESDEvent->GetNumberOfTracks(); ii++) {
    AliESDtrack *track = fESDEvent->GetTrack(ii);
		if (track) {
      if (fTrackCuts->AcceptTrack(track)) {
		    nmbTrk_ST++;
      }
    }
	}
  
	// Martin's selection
	TArrayI Mindices;
	Int_t nMartinSel = fMartinSel->GetNumberOfITSTPCtracks(fESDEvent,Mindices);
  // printf("<I - UserExec> nMartinSel: %i\n",nMartinSel);
    
  // get an TObjArray of tracks with the corresponding TrackStatus
  // The TrackStatus is contained in an array of UInt_t
  // The definition of the TrackStatus bits is given in AliCEPBase.h
  Int_t nTracks = fCEPUtil->AnalyzeTracks(fESDEvent,fTracks,fTrackStatus);

  // use the AliCEPUtils::GetCEPTracks method
  TArrayI *Pindices  = new TArrayI();
  Int_t nPaulSel = fCEPUtil->GetCEPTracks(fESDEvent,fTrackStatus,Pindices);
  if (nMartinSel>0) printf("%i/%i good CEP tracks\n", nPaulSel,nMartinSel);

  const AliMultiplicity *mult = fESDEvent->GetMultiplicity();
  Int_t nTracklets = mult->GetNumberOfTracklets();
  
  
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
    UInt_t mask, pattern;
    TArrayI *masks    = new TArrayI();
    TArrayI *patterns = new TArrayI();
    TArrayI *indices  = new TArrayI();

    // nbad track with !kTTTPCScluster
    mask = AliCEPBase::kTTTPCScluster;
    pattern = 0;
    Int_t nbad = fCEPUtil->countstatus(fTrackStatus,mask,pattern);
    if (nbad==0) fhStatsFlow->Fill(AliCEPBase::kBinSharedCluster);
  
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

    printf("<I - UserExec> nTracks         : %i\n",nTracks);
    printf("<I - UserExec> nTracklets      : %i\n",nTracklets);
    printf("<I - UserExec> npureITSTracks  : %i\n",npureITSTracks);
    printf("<I - UserExec> nTrackStandard  : %i\n",nmbTrk_ST);
    printf("<I - UserExec> nBad            : %i\n",nbad);
    printf("<I - UserExec> nTrackSel       : %i\n",nTrackSel);
    printf("<I - UserExec> nTrackFiredChips: %i - event passed test: %i\n",nTracksFiredChips,fpassedFiredChipsTest);
    printf("<I - UserExec> nTrackAccept    : %i %i %i\n",nTrackAccept,nMartinSel,nPaulSel);

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
  // if the number of accepted tracks is within the limits ...
  Bool_t isToSave = 0 < nTracks && nTracks <= fnumTracksMax;

  // OR if the .AliCEPBase::kBitSaveAllEvents bit is set ..
  // isToSave = isToSave ||
  //  fCEPUtil->checkstatus(fAnalysisStatus,
  //    AliCEPBase::kBitSaveAllEvents,AliCEPBase::kBitSaveAllEvents);
  
  // AND DGTrigger
  isToSave = isToSave && isDGTrigger;
  
  // OR kBitSaveAllEvents is set
  isToSave = isToSave ||
    fCEPUtil->checkstatus(fAnalysisStatus,
    AliCEPBase::kBitSaveAllEvents,AliCEPBase::kBitSaveAllEvents);
  if ( isToSave ) {
    
    // update fhStatsFlow
    fhStatsFlow->Fill(AliCEPBase::kBinSaved);

    // fill the CEPEventBuffer
    fCEPEvent->Reset();
    
    // general event informations
    fCEPEvent->SetRunNumber(fRun);
    fCEPEvent->SetEventNumber(fEvent);
    fCEPEvent->SetPeriodNumber(fPeriod);
    fCEPEvent->SetOrbitNumber(fOrbit);
    fCEPEvent->SetBunchCrossNumber(fBunchCross);
    fCEPEvent->SetCollissionType(fEventType);
    fCEPEvent->SetMagnField(fMagField);
    fCEPEvent->SetFiredTriggerClasses(firedTriggerClasses);
  
    // selections and cuts
    fCEPEvent->SetPhyssel(isPhys);
    fCEPEvent->SetEventCutsel(isEventCutsel);
    fCEPEvent->SetPileup(isPileup);
    fCEPEvent->SetClusterCut(isClusterCut);
    fCEPEvent->SetDGTrigger(isDGTrigger);
    fCEPEvent->SetGapCondition(fCurrentGapCondition);

    // number of tracklets and residuals, vertex
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
    for (Int_t ii=0; ii<nTracks; ii++) {
    
      // the original track
      AliESDtrack *tmptrk = (AliESDtrack*) fTracks->At(ii);

      // create new CEPTrackBuffer and fill it
      CEPTrackBuffer *trk = new CEPTrackBuffer();
      
      trk->SetTrackStatus(fTrackStatus->At(ii));
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
	if (const AliESDInputHandler *esdH =
    dynamic_cast<AliESDInputHandler*>(fInputHandler)) {
    fESDEvent = (AliESDEvent*)esdH->GetEvent();
	  fESDRun   = (AliESDRun*)fESDEvent->GetESDRun();
    //fESDRun->Print();
    
	} else {
		printf("<E - CheckInput> No valid event!\n");
    return kFALSE;
	}
  if (!fESDEvent) {
		printf("<E - CheckInput> No valid event!\n");
    return kFALSE;
	}
    
  // check incomplete DAQ, see recommendation at
  // https://twiki.cern.ch/twiki/bin/view/ALICE/PWGPPEvSelRun2pp
  if (fESDEvent->IsIncompleteDAQ()) {
		printf("<E - CheckInput> Incomplete DAQ!\n");
    return kFALSE;
	}

	// get the LHC period name
  TList *uiList = fInputHandler->GetUserInfo();
  AliProdInfo prodInfo(uiList);
  fLHCPeriod = TString(prodInfo.GetLHCPeriod());

  // PID response
	fPIDResponse = (AliPIDResponse*)fInputHandler->GetPIDResponse();
	if (!fPIDResponse){
		printf("<E - CheckInput> No PID information!\n");
		return kFALSE;
	}

	// check magnetic field
  if (TMath::Abs(fESDEvent->GetMagneticField())<1) {
		printf("<E - CheckInput> Wrong Bfield value! %f\n",
      fESDEvent->GetMagneticField());
		return kFALSE;
	}

	// update the run number
  Int_t tmprun = fESDEvent->GetRunNumber();
	if (fRun != tmprun) {
		fRun = tmprun;
		fCEPUtil->SPDLoadGeom(fRun);
	}

	// get MC event
	fMCEvent = MCEvent();
  if (fMCEvent){  
    // Bug fix 28.05.2016 - do not trust to presence of MC handler, check if the content is valid
    //                    - proper solution (autodetection of MC information) to be implemented 
    if (fMCEvent->Stack()==NULL) fMCEvent=NULL; 
  }

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
