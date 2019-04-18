#include <iostream>
#include "TChain.h"
#include "TH1F.h"
#include "TList.h"
#include "TGeoGlobalMagField.h"
#include "AliVEvent.h"
#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"
#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliESDtrack.h"
#include "AliAODEvent.h"
#include "AliHeader.h"
#include "AliMCEvent.h"
#include "AliGenEventHeader.h"
#include "AlidNdPtTools.h"
#include "AliESDtrackCuts.h"
#include "AliMultSelection.h"
#include "AliCentrality.h"
#include "AliEventCuts.h"
#include "AliAnalysisTaskMKBase.h"

class AliAnalysisTaskMKBase;

using namespace std;

/// \cond CLASSIMP
ClassImp(AliAnalysisTaskMKBase)
/// \endcond

//_____________________________________________________________________________

AliAnalysisTaskMKBase::AliAnalysisTaskMKBase() 
    : AliAnalysisTaskSE()
    , fAnalysisManager(0)
    , fInputEventHandler(0)
    , fEventSelected(0)
    , fEvent(0)
    , fESD(0)
    , fAOD(0)
    , fMC(0)
    , fVZERO(0)
    , fMCStack(0)
    , fMCHeader(0)
    , fMCGenHeader(0)
    , fNTracksESD(0)
    , fNTracksAcc(0)
    , fIsMC(kFALSE)
    , fRunNumber(0)
    , fRunNumberString("")
    , fFiredTriggerClasses("")
    , fEventSpecie(0)
    , fOldCentPercentileV0M(-1)
    , fMultPercentileV0M(-1)
    , fEventCutsPassed(kFALSE)
    , fZv(0)
    , fMCzv(0)
    , fMultMB(-1)
    , fMultV0M(-1)
    , fMCb(-1)
    , fMCnPrim(-1)
    , fMCnPrimV0M(-1)
    , fMCnTracks(-1)
    , fIsTrigger(kFALSE)
    , fHasVertex(kFALSE)
    , fIsIncompleteDAQ(kFALSE)
    , fIsSPDClusterVsTrackletBG(kFALSE)
    , fIsFirstEventInChunk(kFALSE)
    , fIsPileUpMV(kFALSE)
    , fIsOutOfBunchPileUp(kFALSE)
    , fIsPileUpEvent(kFALSE)
    , fIsPileUpSPD(kFALSE)
    , fNOTIsVertexSelected2013pA(kFALSE)
    , fIsPileupFromSPD508(kFALSE)
    , fESDTrack(0)
    , fPt(0)
    , fEta(0)
    , fPhi(0)
    , fDCA{0,0}
    , fDCACov{0,0,0}
    , fDCAr(0)
    , fDCAz(0)
    , fSigma1Pt2(0)
    , fSigma1Pt(0)
    , fSigned1Pt(0)
    , f1Pt(0)
    , fMCParticle(0)
    , fMCLabel(0)
    , fMCPt(0)
    , fMCEta(0)
    , fMCPhi(0)
    , fMCisPrim(kFALSE)
    , fMCisSec(kFALSE)
    , fMCisSecDecay(kFALSE)
    , fMCisSecMat(kFALSE)
    , fMCPrimSec(-1)
    , fMCParticleType(AlidNdPtTools::kUndefined)
    , fMCProdcutionType(AlidNdPtTools::kUnknown)
    , fMCPDGCode(0)
    , fMCCharge(-9999)
    , fMCQ(-9999)
    , fMCIsCharged(kFALSE)
    , fMCChargeSign(-9999)
    , fInnerP(0)
    , fTPCinnerP(0)
    , fPtInner(0)
    , fEtaInner(0)
    , fPhiInner(0)
    , fPtInnerTPC(0)
    , fEtaInnerTPC(0)
    , fPhiInnerTPC(0)
    , fDCATPC{0,0}
    , fDCACovTPC{0,0,0}
    , fDCArTPC(0)
    , fDCAzTPC(0)
    , fEventCuts(0)
    , fUseEventCuts(kFALSE)
    , fESDtrackCutsM(0)
    , fESDtrackCuts{0,0,0,0,0,0,0,0,0,0}
    , fAcceptTrack{kFALSE,kFALSE,kFALSE,kFALSE,kFALSE,kFALSE,kFALSE,kFALSE,kFALSE,kFALSE}
    , fMultEstimator("")
    , fCentEstimator("")
    , fTriggerMaskRequired(0)
    , fTriggerMaskRejected(0)
    , fOutputList(0)
    , fLogHist(0)
    , fLogErr(0)
    , fLogEvent(0)
    , fRunHist(0)
    , fRunHistSelected(0)
    , fTrigInfo(0)
    , fTrigHist(0)
{
    // default constructor for root    
}

//_____________________________________________________________________________

AliAnalysisTaskMKBase::AliAnalysisTaskMKBase(const char* name) 
    : AliAnalysisTaskSE(name)    
    , fAnalysisManager(0)
    , fInputEventHandler(0)
    , fEventSelected(0)
    , fEvent(0)
    , fESD(0)
    , fAOD(0)
    , fMC(0)
    , fVZERO(0)
    , fMCStack(0)
    , fMCHeader(0)
    , fMCGenHeader(0)
    , fNTracksESD(0)
    , fNTracksAcc(0)
    , fIsMC(kFALSE)
    , fRunNumber(0)
    , fRunNumberString("")
    , fFiredTriggerClasses("")
    , fEventSpecie(0)
    , fOldCentPercentileV0M(-1)
    , fMultPercentileV0M(-1)
    , fEventCutsPassed(kFALSE)
    , fZv(0)
    , fMCzv(0)
    , fMultMB(-1)
    , fMultV0M(-1)
    , fMCb(-1)
    , fMCnPrim(-1)
    , fMCnPrimV0M(-1)
    , fMCnTracks(-1)
    , fIsTrigger(kFALSE)
    , fHasVertex(kFALSE)
    , fIsIncompleteDAQ(kFALSE)
    , fIsSPDClusterVsTrackletBG(kFALSE)
    , fIsFirstEventInChunk(kFALSE)
    , fIsPileUpMV(kFALSE)
    , fIsOutOfBunchPileUp(kFALSE)
    , fIsPileUpEvent(kFALSE)
    , fIsPileUpSPD(kFALSE)
    , fNOTIsVertexSelected2013pA(kFALSE)
    , fIsPileupFromSPD508(kFALSE)
    , fESDTrack(0)
    , fPt(0)
    , fEta(0)
    , fPhi(0)
    , fDCA{0,0}
    , fDCACov{0,0,0}
    , fDCAr(0)
    , fDCAz(0)
    , fSigma1Pt2(0)
    , fSigma1Pt(0)
    , fSigned1Pt(0)
    , f1Pt(0)
    , fMCParticle(0)
    , fMCLabel(0)
    , fMCPt(0)
    , fMCEta(0)
    , fMCPhi(0)
    , fMCisPrim(kFALSE)
    , fMCisSec(kFALSE)
    , fMCisSecDecay(kFALSE)
    , fMCisSecMat(kFALSE)
    , fMCPrimSec(-1)
    , fMCParticleType(AlidNdPtTools::kUndefined)
    , fMCProdcutionType(AlidNdPtTools::kUnknown)    
    , fMCPDGCode(0)
    , fMCCharge(-9999)
    , fMCQ(-9999)
    , fMCIsCharged(kFALSE)
    , fMCChargeSign(-9999)
    , fInnerP(0)
    , fTPCinnerP(0)
    , fPtInner(0)
    , fEtaInner(0)
    , fPhiInner(0)
    , fPtInnerTPC(0)
    , fEtaInnerTPC(0)
    , fPhiInnerTPC(0)
    , fDCATPC{0,0}
    , fDCACovTPC{0,0,0}
    , fDCArTPC(0)
    , fDCAzTPC(0)
    , fEventCuts(0)
    , fUseEventCuts(kFALSE)    
    , fESDtrackCutsM(0)
    , fESDtrackCuts{0,0,0,0,0,0,0,0,0,0}
    , fAcceptTrack{kFALSE,kFALSE,kFALSE,kFALSE,kFALSE,kFALSE,kFALSE,kFALSE,kFALSE,kFALSE}
    , fMultEstimator("")
    , fCentEstimator("")
    , fTriggerMaskRequired(0)
    , fTriggerMaskRejected(0)
    , fOutputList(0)
    , fLogHist(0)
    , fLogErr(0)
    , fLogEvent(0)
    , fRunHist(0)
    , fRunHistSelected(0)
    , fTrigInfo(0)
    , fTrigHist(0)
{    
    // constructor    
    DefineInput(0, TChain::Class()); 
    DefineOutput(1, TList::Class());
}

//_____________________________________________________________________________

AliAnalysisTaskMKBase::~AliAnalysisTaskMKBase()
{
    // destructor
    if (fOutputList) { delete fOutputList; }
}

//_____________________________________________________________________________

void AliAnalysisTaskMKBase::UserCreateOutputObjects()
{
    // create output list
    fOutputList = new TList(); 
    fOutputList->SetOwner(kTRUE); 
    
    // create defualt histograms    
    
    fLogHist = CreateLogHist("fLogHist");
    fOutputList->Add(fLogHist);

    fLogErr = CreateLogHist("fLogErr");
    fOutputList->Add(fLogErr);
    
    fLogEvent = CreateLogHist("fLogEvent");
    fOutputList->Add(fLogEvent);    
    
    fRunHist = CreateLogHist("fRunHist");
    fOutputList->Add(fRunHist);
    
    fRunHistSelected = CreateLogHist("fRunHistSelected");
    fOutputList->Add(fRunHistSelected);
    
    fTrigInfo = CreateLogHist("fTrigInfo");
    fOutputList->Add(fTrigInfo);
    
    fTrigHist = CreateLogHist("fTrigHist");
    fOutputList->Add(fTrigHist);    


    AddOutput();
    // postdata 
    PostData(1, fOutputList);
}

//_____________________________________________________________________________
void AliAnalysisTaskMKBase::CheckEvent()
{

    // incomplete daq events
    if ((fIsIncompleteDAQ = fESD->IsIncompleteDAQ())) {
        LogEvent("event.IsIncompleteDAQ");
    }

    // background rejection etc.
    AliAnalysisUtils utils;
    
    if ((fIsSPDClusterVsTrackletBG = utils.IsSPDClusterVsTrackletBG(fESD))) {    
        LogEvent("utils.IsSPDClusterVsTrackletBG");
    }   
    // first event in chunk 
    if ((fIsFirstEventInChunk = utils.IsFirstEventInChunk(fESD))) {      
        LogEvent("utils.IsFirstEventInChunk");
    }
    if ((fIsPileUpMV = utils.IsPileUpMV(fESD))) {
        LogEvent("utils.IsPileUpMV");
    }
    if ((fIsOutOfBunchPileUp = utils.IsOutOfBunchPileUp(fESD))) {
        LogEvent("utils.IsOutOfBunchPileUp");
    }
    if ((fIsPileUpEvent = utils.IsPileUpEvent(fESD))) {
        LogEvent("utils.IsPileUpEvent");
    }
    if ((fIsPileUpSPD = utils.IsPileUpSPD(fESD))) {
        LogEvent("utils.IsPileUpSPD");
    }
    if ((fNOTIsVertexSelected2013pA = !utils.IsVertexSelected2013pA(fESD))) {
        LogEvent("utils.NOT.IsVertexSelected2013pA");
    }
    if ((fIsPileupFromSPD508 = fESD->IsPileupFromSPD(5,0.8))) {
        LogEvent("event.IsPileupFromSPD(5,0.8)");
    }
}

//_____________________________________________________________________________

/// Read the data Event and sets pointer to 
/// fEvent, fESD
/// also sets fRunNumber and fRunNumberString and fills fRunHist histogram
/// later to included the functionally to run on AOD
///
/// \return kTRUE if ESD event is present

Bool_t AliAnalysisTaskMKBase::ReadEvent()
{
    fAnalysisManager = AliAnalysisManager::GetAnalysisManager();
    if (!fAnalysisManager) {
        Err("noAliAnalysisManager");  
        fInputEventHandler = 0;
    } else {
        fInputEventHandler = dynamic_cast<AliInputEventHandler*>(fAnalysisManager->GetInputEventHandler());
    }
    if (fInputEventHandler) { 
        fEventSelected = fInputEventHandler->IsEventSelected();
        FillTriggerLog();
    } else {
        Err("noInputEventHandler");  
    }
    
    fEvent = InputEvent();
    if (!fEvent) { Err("noEvent"); return kFALSE; }
    
    fFiredTriggerClasses = fEvent->GetFiredTriggerClasses();
    
    fESD = dynamic_cast<AliESDEvent*>(InputEvent());
    fAOD = dynamic_cast<AliAODEvent*>(InputEvent());   

    Log(fTrigInfo,fFiredTriggerClasses.Data());
    
    
    if (!fESD) { Err("noESD"); return kFALSE; } // for now analyse only ESD
    if (fESD)  { LogEvent("ESD"); }
    
    ReadMCEvent();
    
    CheckEvent();
    
    UInt_t fEventSpecie = fESD->GetEventSpecie();
    
    fRunNumber = fESD->GetRunNumber();
    fRunNumberString = "";
    fRunNumberString += fRunNumber;
    
    Log(fRunHist,fRunNumberString.Data());
    
    InitEvent();       
    InitEventMult();
    InitEventCent();
    InitMCEvent(); 
    InitVZERO();
          
    return kTRUE;
        
}

//_____________________________________________________________________________

/// Read the MC Event and sets pointer to 
/// fMC, fMCHeader, fMCGenHeader, fMCStack
/// sets fIsMC to kTRUE if all are present
///
/// \return kTRUE if MCEvent, MCGenHeader and Stack are all present

Bool_t AliAnalysisTaskMKBase::ReadMCEvent()
{
    fIsMC = kFALSE;
    fMC = MCEvent();
    if (!fMC) { 
        LogEvent("noMC"); 
        return kFALSE; 
    } else {
        LogEvent("MC");
    }
    
    fIsMC = kTRUE;
    
    fMCStack = fMC->Stack();
    if (!fMCStack) { 
        Err("noMCstack"); 
        fIsMC=kFALSE; 
    }
    
    fMCGenHeader = 0;
    fMCHeader = fMC->Header();
    if (!fMCHeader) { 
        Err("noMCHeader"); 
        fIsMC=kFALSE; 
    } else {
        fMCGenHeader = fMCHeader->GenEventHeader();
    }
    
    if (!fMCGenHeader) {
        Err("noMCGenHeader");
        fIsMC = kFALSE;
    } else {
        TString s = "mcHeader=";
        s += fMCGenHeader->GetName();
        LogEvent(s.Data());
    }
    
    return fIsMC;
}

//_____________________________________________________________________________

/// Initialize event-related properties, needs to be called from user
/// sets multiplicities, vertex, etc
///
/// \return kTRUE 

Bool_t AliAnalysisTaskMKBase::InitEvent() 
{
    if (!fESD) { return kFALSE;}
      // this is needed for some esd track cuts, to be on the save side we call it here
      if (!TGeoGlobalMagField::Instance()->GetField()) { fESD->InitMagneticField(); }
      
      fNTracksESD = fESD->GetNumberOfTracks();
      
      // loop over all tracks to get the accepted multiplicity
      // only if track cuts are set
      fNTracksAcc = 0;
      if (fESDtrackCutsM) {
          for (Int_t i = 0; i < fNTracksESD; i++) {
              AliESDtrack* track = static_cast<AliESDtrack*>(fESD->GetTrack(i));
              if (!track) continue;
              if (fESDtrackCutsM->AcceptTrack(track) ) { fNTracksAcc++; }              
          }
      } else {
          fNTracksAcc = fNTracksESD;
          Log("noAliESDTrackCutsM");
      } 
     
      return kTRUE;
}

//_____________________________________________________________________________

/// function to check all triggers and log into trigger histogram

void AliAnalysisTaskMKBase::FillTriggerLog()
{    
    if (fEventSelected & AliVEvent::kMB)                { Log(fTrigHist,"kMB"); }                 // Minimum bias trigger in PbPb 2010-11
    if (fEventSelected & AliVEvent::kINT1)              { Log(fTrigHist,"kINT1"); }               // V0A | V0C | SPD minimum bias trigger
    if (fEventSelected & AliVEvent::kINT7)              { Log(fTrigHist,"kINT7"); }               // V0AND minimum bias trigger
    if (fEventSelected & AliVEvent::kMUON)              { Log(fTrigHist,"kMUON"); }               // Single muon trigger in pp2010-11, INT1 suite
    if (fEventSelected & AliVEvent::kHighMult)          { Log(fTrigHist,"kHighMult"); }           // High-multiplicity SPD trigger
    if (fEventSelected & AliVEvent::kHighMultSPD)       { Log(fTrigHist,"kHighMultSPD"); }        // High-multiplicity SPD trigger
    if (fEventSelected & AliVEvent::kEMC1)              { Log(fTrigHist,"kEMC1"); }               // EMCAL trigger in pp2011, INT1 suite
    if (fEventSelected & AliVEvent::kCINT5)             { Log(fTrigHist,"kCINT5"); }              // V0OR minimum bias trigger
    if (fEventSelected & AliVEvent::kINT5)              { Log(fTrigHist,"kINT5"); }               // V0OR minimum bias trigger
    if (fEventSelected & AliVEvent::kCMUS5)             { Log(fTrigHist,"kCMUS5"); }              // Single muon trigger, INT5 suite
    if (fEventSelected & AliVEvent::kMUSPB)             { Log(fTrigHist,"kMUSPB"); }              // Single muon trigger in PbPb 2011
    if (fEventSelected & AliVEvent::kINT7inMUON)        { Log(fTrigHist,"kINT7inMUON"); }         // INT7 in MUON or MUFAST cluster
    if (fEventSelected & AliVEvent::kMuonSingleHighPt7) { Log(fTrigHist,"kMuonSingleHighPt7"); }  // Single muon high-pt, INT7 suite
    if (fEventSelected & AliVEvent::kMUSH7)             { Log(fTrigHist,"kMUSH7"); }              // Single muon high-pt, INT7 suite
    if (fEventSelected & AliVEvent::kMUSHPB)            { Log(fTrigHist,"kMUSHPB"); }             // Single muon high-pt in PbPb 2011
    if (fEventSelected & AliVEvent::kMuonLikeLowPt7)    { Log(fTrigHist,"kMuonLikeLowPt7"); }     // Like-sign dimuon low-pt, INT7 suite
    if (fEventSelected & AliVEvent::kMUL7)              { Log(fTrigHist,"kMUL7"); }               // Like-sign dimuon low-pt, INT7 suite
    if (fEventSelected & AliVEvent::kMuonLikePB)        { Log(fTrigHist,"kMuonLikePB"); }         // Like-sign dimuon low-pt in PbPb 2011
    if (fEventSelected & AliVEvent::kMuonUnlikeLowPt7)  { Log(fTrigHist,"kMuonUnlikeLowPt7"); }   // Unlike-sign dimuon low-pt, INT7 suite
    if (fEventSelected & AliVEvent::kMUU7)              { Log(fTrigHist,"kMUU7"); }               // Unlike-sign dimuon low-pt, INT7 suite
    if (fEventSelected & AliVEvent::kMuonUnlikePB)      { Log(fTrigHist,"kMuonUnlikePB"); }       // Unlike-sign dimuon low-pt in PbPb 2011
    if (fEventSelected & AliVEvent::kEMC7)              { Log(fTrigHist,"kEMC7"); }               // EMCAL/DCAL L0 trigger, INT7 suite
    if (fEventSelected & AliVEvent::kEMC8)              { Log(fTrigHist,"kEMC8"); }               // EMCAL/DCAL L0 trigger, INT8 suite
    if (fEventSelected & AliVEvent::kMUS7)              { Log(fTrigHist,"kMUS7"); }               // Single muon low-pt, INT7 suite
    if (fEventSelected & AliVEvent::kMuonSingleLowPt7)  { Log(fTrigHist,"kMuonSingleLowPt7"); }   // Single muon low-pt, INT7 suite
    if (fEventSelected & AliVEvent::kPHI1)              { Log(fTrigHist,"kPHI1"); }               // PHOS L0 trigger in pp2011, INT1 suite
    if (fEventSelected & AliVEvent::kPHI7)              { Log(fTrigHist,"kPHI7"); }               // PHOS trigger, INT7 suite
    if (fEventSelected & AliVEvent::kPHI8)              { Log(fTrigHist,"kPHI8"); }               // PHOS trigger, INT8 suite
    if (fEventSelected & AliVEvent::kPHOSPb)            { Log(fTrigHist,"kPHOSPb"); }             // PHOS trigger in PbPb 2011
    if (fEventSelected & AliVEvent::kEMCEJE)            { Log(fTrigHist,"kEMCEJE"); }             // EMCAL/DCAL L1 jet trigger
    if (fEventSelected & AliVEvent::kEMCEGA)            { Log(fTrigHist,"kEMCEGA"); }             // EMCAL/DCAL L1 gamma trigger
    if (fEventSelected & AliVEvent::kHighMultV0)        { Log(fTrigHist,"kHighMultV0"); }         // High-multiplicity V0 trigger
    if (fEventSelected & AliVEvent::kCentral)           { Log(fTrigHist,"kCentral"); }            // Central trigger in PbPb 2011
    if (fEventSelected & AliVEvent::kSemiCentral)       { Log(fTrigHist,"kSemiCentral"); }        // Semicentral trigger in PbPb 2011
    if (fEventSelected & AliVEvent::kDG)                { Log(fTrigHist,"kDG"); }                 // Double gap diffractive
    if (fEventSelected & AliVEvent::kDG5)               { Log(fTrigHist,"kDG5"); }                // Double gap diffractive
    if (fEventSelected & AliVEvent::kZED)               { Log(fTrigHist,"kZED"); }                // ZDC electromagnetic dissociation
    if (fEventSelected & AliVEvent::kSPI7)              { Log(fTrigHist,"kSPI7"); }               // Power interaction trigger
    if (fEventSelected & AliVEvent::kSPI)               { Log(fTrigHist,"kSPI"); }                // Power interaction trigger
    if (fEventSelected & AliVEvent::kINT8)              { Log(fTrigHist,"kINT8"); }               // 0TVX trigger
    if (fEventSelected & AliVEvent::kMuonSingleLowPt8)  { Log(fTrigHist,"kMuonSingleLowPt8"); }   // Single muon low-pt, INT8 suite
    if (fEventSelected & AliVEvent::kMuonSingleHighPt8) { Log(fTrigHist,"kMuonSingleHighPt8"); }  // Single muon high-pt, INT8 suite
    if (fEventSelected & AliVEvent::kMuonLikeLowPt8)    { Log(fTrigHist,"kMuonLikeLowPt8"); }     // Like-sign dimuon low-pt, INT8 suite
    if (fEventSelected & AliVEvent::kMuonUnlikeLowPt8)  { Log(fTrigHist,"kMuonUnlikeLowPt8"); }   // Unlike-sign dimuon low-pt, INT8 suite
    if (fEventSelected & AliVEvent::kMuonUnlikeLowPt0)  { Log(fTrigHist,"kMuonUnlikeLowPt0"); }   // Unlike-sign dimuon low-pt, no additional L0 requirement
    if (fEventSelected & AliVEvent::kUserDefined)       { Log(fTrigHist,"kUserDefined"); }        // Set when custom trigger classes are set in AliPhysicsSelection
    if (fEventSelected & AliVEvent::kTRD)               { Log(fTrigHist,"kTRD"); }                // TRD trigger
    if (fEventSelected & AliVEvent::kMuonCalo)          { Log(fTrigHist,"kMuonCalo"); }           // Muon-calo triggers
    if (fEventSelected & AliVEvent::kCaloOnly)          { Log(fTrigHist,"kCaloOnly"); }           // MB, EMCAL and PHOS triggers in CALO or CALOFAST cluster    
    if (fEventSelected & AliVEvent::kFastOnly)          { Log(fTrigHist,"kFastOnly"); }           // The fast cluster fired. This bit is set in to addition another trigger bit, e.g. kMB
    if (fEventSelected & AliVEvent::kAny)               { Log(fTrigHist,"kAny"); }                // to accept any defined trigger
    if (fEventSelected & AliVEvent::kAnyINT)            { Log(fTrigHist,"kAnyINT"); }             // to accept any interaction (aka minimum bias) trigger    
}

//_____________________________________________________________________________

Bool_t AliAnalysisTaskMKBase::InitMCEvent() 
{
    if (!fIsMC) return kFALSE;
    fMCnTracks = fMC->GetNumberOfTracks();
    return kTRUE;
}

//_____________________________________________________________________________

Bool_t AliAnalysisTaskMKBase::InitVZERO()
{
    fVZERO = fESD->GetVZEROData();
    if (!fVZERO) { 
        Log("noVZEROData");
        return kFALSE;
    }
    return kTRUE;
}

//_____________________________________________________________________________

Bool_t AliAnalysisTaskMKBase::InitTrack()
{    
    if (!fESDTrack) { return kFALSE; }
    fPt = fESDTrack->Pt();
    fEta = fESDTrack->Eta();
    fPhi = fESDTrack->Phi();
    fESDTrack->GetImpactParameters(fDCA,fDCACov);
    fDCAr = fDCA[0];
    fDCAz = fDCA[1];
    fSigma1Pt2 = fESDTrack->GetSigma1Pt2();
    if (fSigma1Pt2 < 0) { 
        Err("Sigma1Pt2<0");        
    }
    fSigma1Pt =  TMath::Sqrt(fSigma1Pt2);
    fSigned1Pt = fESDTrack->GetSigned1Pt();
    f1Pt = TMath::Abs(fSigned1Pt);
    return kTRUE;
}

//_____________________________________________________________________________

Bool_t AliAnalysisTaskMKBase::InitTrackCuts()
{        
    fAcceptTrackM = (fESDtrackCutsM)?  fESDtrackCutsM->AcceptTrack(fESDTrack) : kFALSE;
    for (int i=0; i<10; i++) { fAcceptTrack[i] = (fESDtrackCuts[i])? fESDtrackCuts[i]->AcceptTrack(fESDTrack) : kFALSE;}   
    return kTRUE;
}

//_____________________________________________________________________________

Bool_t AliAnalysisTaskMKBase::InitMCTrack()
{ 
    if (!fIsMC) return kFALSE;
    fMCPrimSec = -1; 
    if (!fESDTrack) return kFALSE;
    if (!fMC) return kFALSE;
    
    fMCLabel = TMath::Abs(fESDTrack->GetLabel());
    if (fMCLabel < 0) { Log("tracklabel<0"); }
    fMCParticle  = static_cast<AliMCParticle*>(fMC->GetTrack(fMCLabel));
    if (!fMCParticle) { 
        Err("particleNOTinStack"); 
        return kFALSE;
    }    
    
    return InitMCParticle();
}

//_____________________________________________________________________________

Bool_t AliAnalysisTaskMKBase::InitMCParticle()
{ 
     if (!fIsMC) return kFALSE;
    // set all mc particle related things
     
    fMCPt  = fMCParticle->Pt(); 
    fMCEta = fMCParticle->Eta(); 
    fMCPhi = fMCParticle->Phi(); 
    fMCCharge =  fMCParticle->Charge();
    fMCQ =  fMCCharge/3.0;    
    if (fMCCharge == 0) {
        fMCIsCharged = kFALSE;
        fMCChargeSign = 0;
    } else {
        fMCIsCharged = kTRUE;
        fMCChargeSign = (fMCCharge > 0)? +1 : -1;
    }
    
    fMCisPrim     = fMC->IsPhysicalPrimary(fMCLabel);
    
    Bool_t isPhysPrimStack = fMCStack->IsPhysicalPrimary(fMCLabel);
    if (fMCisPrim != isPhysPrimStack) { Log("PhysPrimMismatch"); }
    
    fMCisSecDecay = fMC->IsSecondaryFromWeakDecay(fMCLabel);
    fMCisSecMat   = fMC->IsSecondaryFromMaterial(fMCLabel);        
    fMCisSec      = fMCisSecMat || fMCisSecDecay;    
    if (fMCisPrim)     { fMCPrimSec = 0; fMCProdcutionType = AlidNdPtTools::kPrim; }
    if (fMCisSecDecay) { fMCPrimSec = 1; fMCProdcutionType = AlidNdPtTools::kSecDecay; }
    if (fMCisSecMat)   { fMCPrimSec = 2; fMCProdcutionType = AlidNdPtTools::kSecMaterial; }
    if (fMCPrimSec == -1)             { Err("NOTprimORsec"); }
    if (fMCisPrim && fMCisSec)        { Err("primANDsec"); }
    if (fMCisSecDecay && fMCisSecMat) { Err("decayANDmat"); }
    fMCPDGCode = fMCParticle->PdgCode();
    fMCParticleType = AlidNdPtTools::ParticleTypeFromPDG(fMCPDGCode);
    
    return kTRUE;
}

//_____________________________________________________________________________

Bool_t AliAnalysisTaskMKBase::InitEventCent()
{ 
    if (!fESD) { Err("MultCentRequiresESD"); return kFALSE; }
    AliCentrality *centrality = fESD->GetCentrality();
    fOldCentPercentileV0M = -1;
    if (!centrality) { 
        LogEvent("noOldCentrality");
        return kFALSE;
    } else {
        fOldCentPercentileV0M = centrality->GetCentralityPercentile("V0M");
        if (fOldCentPercentileV0M < 0.) { LogEvent("fOldCentPercentileV0M<0"); }
        if (fOldCentPercentileV0M > 100.) { LogEvent("fOldCentPercentileV0M>100"); }
        if (fOldCentPercentileV0M == 0.) { LogEvent("fOldCentPercentileV0M==0"); }
        if (fOldCentPercentileV0M == 100.) { LogEvent("fOldCentPercentileV0M==100"); }
    }
    return kTRUE;    
}
//_____________________________________________________________________________

Bool_t AliAnalysisTaskMKBase::InitEventMult()
{
    fEventCutsPassed = kFALSE;
    fMultPercentileV0M = -1;
    AliMultSelection *multSelection = 0; 
    multSelection = (AliMultSelection*) fEvent->FindListObject("MultSelection");
    if(!multSelection) {
        LogEvent("noMultSelection");
        return kFALSE;
    } else {
        fMultPercentileV0M = multSelection->GetMultiplicityPercentile("V0M");
        if (fMultPercentileV0M < 0.) { LogEvent("fMultPercentileV0M<0"); }
        if (fMultPercentileV0M > 100.) { LogEvent("fMultPercentileV0M>100"); }
        if (fMultPercentileV0M == 0.) { LogEvent("fMultPercentileV0M==0"); }
        if (fMultPercentileV0M == 100.) { LogEvent("fMultPercentileV0M==100"); }        
    }
    
    // AliEventCut needs multiplcity task, otherwise it crashes
    // so we put it here  
    fEventCutsPassed = fEventCuts.AcceptEvent(fEvent);
    if (fEventCutsPassed) { LogEvent("AliEventCutsPassed"); }         
    
    return kTRUE;
}

//_____________________________________________________________________________

Bool_t AliAnalysisTaskMKBase::InitTrackIP()
{
    fInnerP = fESDTrack->GetInnerParam();
    if (!fInnerP) { 
        Err("noInnerParam"); 
        fPtInner = 0;
        fEtaInner = 0;
        fPhiInner = 0;
        return kFALSE;
    } else {
        fPtInner = fInnerP->Pt();
        fEtaInner = fInnerP->Eta();
        fPhiInner = fInnerP->Phi();
    }
    return kTRUE;
}

//_____________________________________________________________________________

Bool_t AliAnalysisTaskMKBase::InitTrackTPC()
{
    fESDTrack->GetImpactParametersTPC(fDCATPC,fDCACovTPC);
    fDCArTPC = fDCATPC[0];
    fDCAzTPC = fDCATPC[1];

    fTPCinnerP = fESDTrack->GetTPCInnerParam();
    if (!fTPCinnerP) { 
        Err("noTPCinner"); 
        fPtInnerTPC  = 0;
        fEtaInnerTPC = 0;
        fPhiInnerTPC = 0;
        return kFALSE;
    } else {
        fPtInnerTPC  = fTPCinnerP->Pt();
        fEtaInnerTPC = fTPCinnerP->Eta();
        fPhiInnerTPC = fTPCinnerP->Phi();
    }            
    return kTRUE;
}
        
//_____________________________________________________________________________
        
        
void AliAnalysisTaskMKBase::UserExec(Option_t *)
{   
    // call analysis for data and mc
    ReadEvent();    
    
    AnaEvent();
    if (fIsMC) {
        AnaMCEvent(); 
    } else {
        AnaDATAEvent();
    }
    
    // postdata
    PostData(1, fOutputList);
}

//_____________________________________________________________________________

void AliAnalysisTaskMKBase::LoopOverAllTracks()
{    
    fNTracksESD = fESD->GetNumberOfTracks();
    for (Int_t i = 0; i < fNTracksESD; i++) {
        fESDTrack = static_cast<AliESDtrack*>(fESD->GetTrack(i));
        if (!fESDTrack) { Err("noESDtrack"); continue; }
        InitTrack();
        InitTrackCuts();
        InitMCTrack();
        InitTrackIP();
        InitTrackTPC();
        AnaTrack();
        if (fIsMC) AnaMCTrack();
    }
}


//_____________________________________________________________________________

void AliAnalysisTaskMKBase::LoopOverAllParticles()
{    
    if (!fIsMC) return;
    fMCnTracks = fMC->GetNumberOfTracks();
    for (Int_t i = 0; i < fMCnTracks; i++) {
        fMCParticle  = static_cast<AliMCParticle*>(fMC->GetTrack(i));
        if (!fMCParticle) { Err("noMCParticle"); continue; }         
        fMCLabel = i;
        InitMCParticle();
        AnaMCParticle();        
    }
}


//_____________________________________________________________________________

void AliAnalysisTaskMKBase::Terminate(Option_t *)
{
    // terminate  
    cout<<"AliAnalysisTaskMKBase terminated."<<endl;
}

//_____________________________________________________________________________

void AliAnalysisTaskMKBase::FinishTaskOutput()
{
    // finish task output
}

//_____________________________________________________________________________

/// static factory function that serves as AddTask macro
/// Created and configures a task and attaches it
/// to the AnalysisManager
///
/// \return the created task

AliAnalysisTaskMKBase* AliAnalysisTaskMKBase::AddTaskMKBase(const char* name, const char* outfile) 
{
    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    if (!mgr) {
        ::Error("AddTaskMKBase", "No analysis manager to connect to.");
        return 0;
    }

    // Check the analysis type using the event handlers connected to the analysis manager.
    //==============================================================================
    if (!mgr->GetInputEventHandler()) {
        ::Error("AddTaskMKBase", "This task requires an input event handler");
        return NULL;
    }
    
    // Setup output file
    //===========================================================================
    TString fileName = AliAnalysisManager::GetCommonFileName();        
    fileName += ":";
    fileName += name;  // create a subfolder in the file
    if (outfile) { // if a finename is given, use that one
        fileName = TString(outfile);
    }
    

    // create the task
    //===========================================================================
    AliAnalysisTaskMKBase *task = new AliAnalysisTaskMKBase(name);  
    if (!task) { return 0; }
    
    // configure the task
    //===========================================================================
    task->SelectCollisionCandidates(AliVEvent::kAnyINT);
    task->SetESDtrackCutsM(AlidNdPtTools::CreateESDtrackCuts("default"));
    
    // attach the task to the manager and configure in and ouput
    //===========================================================================
    mgr->AddTask(task);
    mgr->ConnectInput(task,0,mgr->GetCommonInputContainer());
    mgr->ConnectOutput(task,1,mgr->CreateContainer(name, TList::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));
    
  return task;
}

//_____________________________________________________________________________
