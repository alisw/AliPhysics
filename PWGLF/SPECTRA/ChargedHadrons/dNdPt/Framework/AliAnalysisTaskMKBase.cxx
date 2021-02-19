#include "AliAODEvent.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisTask.h"
#include "AliCentrality.h"
#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliESDVertex.h"
#include "AliESDtrack.h"
#include "AliESDtrackCuts.h"
#include "AliEventCuts.h"
#include "AliAnalysisUtils.h"
#include "AliGenDPMjetEventHeader.h"
#include "AliGenEventHeader.h"
#include "AliGenPythiaEventHeader.h"
#include "AliHeader.h"
#include "AliInputEventHandler.h"
#include "AliMCEvent.h"
#include "AliMultEstimator.h"
#include "AliMultSelection.h"
#include "AliVEvent.h"
#include "AlidNdPtTools.h"
#include "TChain.h"
#include "TGeoGlobalMagField.h"
#include "TH1.h"
#include "TList.h"
#include <iostream>
#include "AliExternalTrackParam.h"

#include "AliAnalysisTaskMKBase.h"


//****************************************************************************************
/**
 * Default constructor.
 */
//****************************************************************************************
AliAnalysisTaskMKBase::AliAnalysisTaskMKBase()
: AliAnalysisTaskSE(), fAnalysisManager(0), fInputEventHandler(0),
fEventSelected(0), fEvent(0), fESD(0), fAOD(0), fMC(0), fVZERO(0),
fMCStack(0), fMCHeader(0), fMCGenHeader(0), fMCGenHeaderPythia(0),
fMCGenHeaderDPMjet(0), fMCEventType(AlidNdPtTools::kInvalidProcess),
fMCProcessTypeFlag(0), fMultSelection(0), fCentrality(0), fNTracksESD(0),
fNTracksAcc(0), fIsMC(kFALSE), fRunNumber(0), fRunNumberString(""),
fTimeStamp(0), fEventNumberInFile(0), fFiredTriggerClasses(""),
fEventSpecie(0), fOldCentPercentileV0M(-1), fMultPercentileV0M(-1),
fIsAcceptedAliEventCuts(kFALSE), fVtx(0), fXv(0), fYv(0), fZv(0),
fXvRes(0), fYvRes(0), fZvRes(0), fVtxNContrib(0), fVtxStatus(kFALSE),
fVtxDispersion(0), fUsedVtxTRK(kFALSE), fUsedVtxSPD(kFALSE),
fUsedVtxTPC(kFALSE), fVtxTRK(0), fXvTRK(0), fYvTRK(0), fZvTRK(0),
fXvResTRK(0), fYvResTRK(0), fZvResTRK(0), fVtxNContribTRK(0),
fVtxStatusTRK(kFALSE), fVtxDispersionTRK(0), fVtxSPD(0), fXvSPD(0),
fYvSPD(0), fZvSPD(0), fXvResSPD(0), fYvResSPD(0), fZvResSPD(0),
fVtxNContribSPD(0), fVtxStatusSPD(kFALSE), fVtxDispersionSPD(0),
fVtxTPC(0), fXvTPC(0), fYvTPC(0), fZvTPC(0), fXvResTPC(0), fYvResTPC(0),
fZvResTPC(0), fVtxNContribTPC(0), fVtxStatusTPC(kFALSE),
fVtxDispersionTPC(0), fMCxv(0), fMCyv(0), fMCzv(0), fMultMB(-1),
fMultV0A(-1), fMultV0C(-1), fMultV0M(-1), fMultV0MmultSelection(-1),
fMCb(-1), fMCnPrimPtCut(-1), fMCnPrim10(-1), fMCnPrim08(-1),
fMCnPrim05(-1), fMCnPrimV0M(-1), fMCnTracks(-1), fMCnPrim(-1),
fIsTrigger(kFALSE), fHasVertex(kFALSE), fIsIncompleteDAQ(kFALSE),
fIsSPDClusterVsTrackletBG(kFALSE), fIsFirstEventInChunk(kFALSE),
fIsPileUpMV(kFALSE), fIsOutOfBunchPileUp(kFALSE), fIsPileUpEvent(kFALSE),
fIsPileUpSPD(kFALSE), fIsVertexRejected2013pA(kFALSE),
fIsPileupFromSPD508(kFALSE), fIsEventAccepted(kFALSE), fESDTrack(0),
fPt(0), fP(0), fEta(0), fPhi(0), fDCA{0, 0}, fDCACov{0, 0, 0}, fDCAr(0),
fDCAz(0), fSigma1Pt2(0), fSigma1Pt(0), fSigned1Pt(0), f1Pt(0),
fChargeSign(0), fTPCSignalN(0), fX(0.), fY(0.), fZ(0.), fAlpha(0.), fSnp(0.),
fTgl(0.), fFlags(0), fITSFoundClusters(0.), fITSChi2PerCluster(0.),
fITSClusterMap(0), fTPCFindableClusters(0.), fTPCFoundClusters(0.), fTPCSharedClusters(0.),
fTPCFractionSharedClusters(0.), fTPCCrossedRows(0.),
fTPCCrossedRowsOverFindableClusters(0.), fTPCChi2PerCluster(0.), fTPCGoldenChi2(0.), fTPCGeomLength(0.),
fMCParticle(0), fMCLabel(0), fMCPt(0),
fMCEta(0), fMCPhi(0), fMCisPrim(kFALSE), fMCisSec(kFALSE),
fMCisSecDecay(kFALSE), fMCisSecMat(kFALSE), fMCPrimSec(-1), fMCPileUpTrack(0),
fMCParticleType(AlidNdPtTools::kUndefined),
fMCProdcutionType(AlidNdPtTools::kUnknown), fMCPDGCode(0),
fMCCharge(-9999), fMCQ(-9999), fMCIsCharged(kFALSE), fMCChargeSign(-9999),
fInnerP(0), fTPCinnerP(0), fPtInner(0), fEtaInner(0), fPhiInner(0),
fZInner(0), fPtInnerTPC(0), fEtaInnerTPC(0),
fPhiInnerTPC(0),fZInnerTPC(0), fDCATPC{0, 0}, fDCACovTPC{0, 0, 0}, fDCArTPC(0),
fDCAzTPC(0), fEventCuts(), fUseEventCuts(kFALSE), fESDtrackCutsM(0),
fAcceptTrackM(kFALSE), fESDtrackCuts{0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
fAcceptTrack{kFALSE, kFALSE, kFALSE, kFALSE, kFALSE,
  kFALSE, kFALSE, kFALSE, kFALSE, kFALSE},
fMultEstimator(""), fCentEstimator(""), fTriggerMaskRequired(0),
fTriggerMaskRejected(0), fInternalLoop(kFALSE), fOutputList(0),
fLogHist(0), fLogErr(0), fLogEvent(0), fRunHist(0), fRunHistSelected(0),
fTrigInfo(0), fTrigInfoSelected(0), fTrigHist(0), fTrigHistSelected(0),
fCentralityEstimator(AliAnalysisTaskMKBase::CentralityEstimator::kV0M),
fUseBaseOutput(kFALSE), fNeedEventVertex(kFALSE), fNeedEventCent(kFALSE),
fNeedEventMult(kFALSE), fNeedEventVZERO(kFALSE), fNeedEventQA(kFALSE), fNeedTrackIP(kFALSE),
fNeedTrackTPC(kFALSE), fNeedTrackPID(kFALSE), fNeedTrackQA(kFALSE), fSkipMCtruth(kFALSE)
{
}

//****************************************************************************************
/**
 * ROOT IO constructor.
 */
//****************************************************************************************
AliAnalysisTaskMKBase::AliAnalysisTaskMKBase(const char* name)
: AliAnalysisTaskSE(name), fAnalysisManager(0), fInputEventHandler(0),
fEventSelected(0), fEvent(0), fESD(0), fAOD(0), fMC(0), fVZERO(0),
fMCStack(0), fMCHeader(0), fMCGenHeader(0), fMCGenHeaderPythia(0),
fMCGenHeaderDPMjet(0), fMCEventType(AlidNdPtTools::kInvalidProcess),
fMCProcessTypeFlag(0), fMultSelection(0), fCentrality(0), fNTracksESD(0),
fNTracksAcc(0), fIsMC(kFALSE), fRunNumber(0), fRunNumberString(""),
fTimeStamp(0), fEventNumberInFile(0), fFiredTriggerClasses(""),
fEventSpecie(0), fOldCentPercentileV0M(-1), fMultPercentileV0M(-1),
fIsAcceptedAliEventCuts(kFALSE), fVtx(0), fXv(0), fYv(0), fZv(0),
fXvRes(0), fYvRes(0), fZvRes(0), fVtxNContrib(0), fVtxStatus(kFALSE),
fVtxDispersion(0), fUsedVtxTRK(kFALSE), fUsedVtxSPD(kFALSE),
fUsedVtxTPC(kFALSE), fVtxTRK(0), fXvTRK(0), fYvTRK(0), fZvTRK(0),
fXvResTRK(0), fYvResTRK(0), fZvResTRK(0), fVtxNContribTRK(0),
fVtxStatusTRK(kFALSE), fVtxDispersionTRK(0), fVtxSPD(0), fXvSPD(0),
fYvSPD(0), fZvSPD(0), fXvResSPD(0), fYvResSPD(0), fZvResSPD(0),
fVtxNContribSPD(0), fVtxStatusSPD(kFALSE), fVtxDispersionSPD(0),
fVtxTPC(0), fXvTPC(0), fYvTPC(0), fZvTPC(0), fXvResTPC(0), fYvResTPC(0),
fZvResTPC(0), fVtxNContribTPC(0), fVtxStatusTPC(kFALSE),
fVtxDispersionTPC(0), fMCxv(0), fMCyv(0), fMCzv(0), fMultMB(-1),
fMultV0A(-1), fMultV0C(-1), fMultV0M(-1), fMultV0MmultSelection(-1),
fMCb(-1), fMCnPrimPtCut(-1), fMCnPrim10(-1), fMCnPrim08(-1),
fMCnPrim05(-1), fMCnPrimV0M(-1), fMCnTracks(-1), fMCnPrim(-1),
fIsTrigger(kFALSE), fHasVertex(kFALSE), fIsIncompleteDAQ(kFALSE),
fIsSPDClusterVsTrackletBG(kFALSE), fIsFirstEventInChunk(kFALSE),
fIsPileUpMV(kFALSE), fIsOutOfBunchPileUp(kFALSE), fIsPileUpEvent(kFALSE),
fIsPileUpSPD(kFALSE), fIsVertexRejected2013pA(kFALSE),
fIsPileupFromSPD508(kFALSE), fIsEventAccepted(kFALSE), fESDTrack(0),
fPt(0), fP(0), fEta(0), fPhi(0), fDCA{0, 0}, fDCACov{0, 0, 0}, fDCAr(0),
fDCAz(0), fSigma1Pt2(0), fSigma1Pt(0), fSigned1Pt(0), f1Pt(0),
fChargeSign(0), fTPCSignalN(0),
fX(0.), fY(0.), fZ(0.), fAlpha(0.), fSnp(0.), fTgl(0.), fFlags(0), fITSFoundClusters(0.), fITSChi2PerCluster(0.),
fITSClusterMap(0), fTPCFindableClusters(0.), fTPCFoundClusters(0.), fTPCSharedClusters(0.),
fTPCFractionSharedClusters(0.), fTPCCrossedRows(0.),
fTPCCrossedRowsOverFindableClusters(0.), fTPCChi2PerCluster(0.), fTPCGoldenChi2(0.), fTPCGeomLength(0.),
fMCParticle(0), fMCLabel(0), fMCPt(0),
fMCEta(0), fMCPhi(0), fMCisPrim(kFALSE), fMCisSec(kFALSE),
fMCisSecDecay(kFALSE), fMCisSecMat(kFALSE), fMCPrimSec(-1), fMCPileUpTrack(0),
fMCParticleType(AlidNdPtTools::kUndefined),
fMCProdcutionType(AlidNdPtTools::kUnknown), fMCPDGCode(0),
fMCCharge(-9999), fMCQ(-9999), fMCIsCharged(kFALSE), fMCChargeSign(-9999),
fInnerP(0), fTPCinnerP(0), fPtInner(0), fEtaInner(0), fPhiInner(0),
fZInner(0), fPtInnerTPC(0), fEtaInnerTPC(0),
fPhiInnerTPC(0),fZInnerTPC(0), fDCATPC{0, 0}, fDCACovTPC{0, 0, 0}, fDCArTPC(0),
fDCAzTPC(0), fEventCuts(), fUseEventCuts(kFALSE), fESDtrackCutsM(0),
fAcceptTrackM(kFALSE), fESDtrackCuts{0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
fAcceptTrack{kFALSE, kFALSE, kFALSE, kFALSE, kFALSE,
  kFALSE, kFALSE, kFALSE, kFALSE, kFALSE},
fMultEstimator(""), fCentEstimator(""), fTriggerMaskRequired(0),
fTriggerMaskRejected(0), fInternalLoop(kFALSE), fOutputList(0),
fLogHist(0), fLogErr(0), fLogEvent(0), fRunHist(0), fRunHistSelected(0),
fTrigInfo(0), fTrigInfoSelected(0), fTrigHist(0), fTrigHistSelected(0),
fCentralityEstimator(AliAnalysisTaskMKBase::CentralityEstimator::kV0M),
fUseBaseOutput(kFALSE), fNeedEventVertex(kFALSE), fNeedEventCent(kFALSE),
fNeedEventMult(kFALSE), fNeedEventVZERO(kFALSE), fNeedEventQA(kFALSE), fNeedTrackIP(kFALSE),
fNeedTrackTPC(kFALSE), fNeedTrackPID(kFALSE), fNeedTrackQA(kFALSE), fSkipMCtruth(kFALSE)
{
  DefineInput(0, TChain::Class());
  DefineOutput(1, TList::Class());
}

//****************************************************************************************
/**
 * Destructor.
 */
//****************************************************************************************
AliAnalysisTaskMKBase::~AliAnalysisTaskMKBase()
{
  if (fOutputList) delete fOutputList;
}

//****************************************************************************************
/**
 * Adds some output of general interest. Can be disabled via SetUseBaseOutput(false).
 */
//****************************************************************************************
void AliAnalysisTaskMKBase::BaseAddOutput()
{
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
  
  fTrigInfoSelected = CreateLogHist("fTrigInfoSelected");
  fOutputList->Add(fTrigInfoSelected);
  
  fTrigHist = CreateLogHist("fTrigHist");
  fOutputList->Add(fTrigHist);
  
  fTrigHistSelected = CreateLogHist("fTrigHistSelected");
  fOutputList->Add(fTrigHistSelected);
}

//****************************************************************************************
/**
 * Function to be called per event.
 */
//****************************************************************************************
void AliAnalysisTaskMKBase::UserCreateOutputObjects()
{
  // create output list
  fOutputList = new TList();
  fOutputList->SetOwner(kTRUE);
  
  // add default histograms
  if (fUseBaseOutput) BaseAddOutput();
  
  // add user histograms
  AddOutput();
  // postdata
  PostData(1, fOutputList);
}

//****************************************************************************************
/**
 * Internal function that checks event quality and sets the corresponding flags.
 *
 * REQUIRES:
 * fESD, fMultSelection, vertex infos
 *
 * SETS:
 * fIsIncompleteDAQ, fIsSPDClusterVsTrackletBG, fIsFirstEventInChunk,
 * fIsPileUpMV, fIsOutOfBunchPileUp, fIsPileUpEvent, fIsPileUpSPD,
 * fIsVertexRejected2013pA, fIsPileupFromSPD508
 */
//****************************************************************************************
void AliAnalysisTaskMKBase::InitEventQA()
{
  
  // incomplete daq events
  if ((fIsIncompleteDAQ = fESD->IsIncompleteDAQ())) {
    LogEvent("event.IsIncompleteDAQ");
  }
  
  // background rejection etc. using AliAnalysisUtils
  AliAnalysisUtils utils; // FIXME: this is constructed per event with all its 20 members...
  
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
  if ((fIsVertexRejected2013pA = !utils.IsVertexSelected2013pA(fESD))) {
    LogEvent("utils.NOT.IsVertexSelected2013pA");
  }
  if ((fIsPileupFromSPD508 = fESD->IsPileupFromSPD(5, 0.8))) {
    LogEvent("event.IsPileupFromSPD(5,0.8)");
  }
  
  // for vertex cuts check and log, but do not cut or set flags (for now)
  if (fVtxStatus) {
    if (fZvRes > 0.25) {
      LogEvent("PrimVtx.ResolutionGreater0.25");
    }
  }
  if (fVtxStatusSPD)
  {
    if (fVtxSPD->IsFromVertexerZ() && fVtxDispersionSPD > 0.04) {
      LogEvent("PrimVtxSPD.fromVertexerZ&DispersionGreater0.04");
    }
    if (fVtxDispersionSPD > 0.04) {
      LogEvent("PrimVtxSPD.DispersionGreater0.04");
    }
    if (fVtxDispersionSPD > 0.02) {
      LogEvent("PrimVtxSPD.DispersionGreater0.02");
    }
  }
  if (fVtxStatusTRK && fVtxStatusSPD)
  {
    if (TMath::Abs(fZvTRK - fZvSPD) > 0.5) {
      LogEvent("abs(zvTRK-zvSPD)Greater0.5");
    }
  }
}

//****************************************************************************************
/**
 * Reads the data Event and sets pointer to  fEvent and fESD.
 * Sets fRunNumber and fRunNumberString and fills fRunHist histogram.
 * Functionally to run on AOD to be included later.
 * Sets:  fAnalysisManager, fInputEventHandler, fEventSelected, fEvent, fAOD, fESD
 * \return kTRUE if ESD event is present
 */
//****************************************************************************************
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
  } else {
    Err("noInputEventHandler");
    fEventSelected = 0;
  }
  fEvent = InputEvent();
  if (!fEvent) {
    Err("noEvent");
    return kFALSE;
  }
  
  fFiredTriggerClasses = fEvent->GetFiredTriggerClasses();
  
  fESD = dynamic_cast<AliESDEvent*>(InputEvent());
  fAOD = dynamic_cast<AliAODEvent*>(InputEvent());
  
  // for now analyse only ESDs
  if (!fESD) {
    Err("noESD");
    return kFALSE;
  }
  if (fESD) {
    LogEvent("ESD");
  }
  if (fAOD) {
    LogEvent("AOD");
  }
  
  // read the mc event and set all mc related properties
  if(!fSkipMCtruth) ReadMCEvent();
  
  // set all the event related properties
  InitEvent();
  
  return kTRUE;
}

//****************************************************************************************
/**
 * Read the MC Event and sets pointer to
 * Sets fMC, fMCHeader, fMCGenHeader, fMCStack
 * fIsMC to kTRUE if all are present
 * \return kTRUE if MCEvent, MCGenHeader and Stack are all present
 */
//****************************************************************************************
Bool_t AliAnalysisTaskMKBase::ReadMCEvent()
{
  // FIXME: no need to determine if it is mc for each and every event!!
  fIsMC = kFALSE;
  fMC = MCEvent();
  // FIXME: remove useless logging
  if (!fMC) {
    LogEvent("noMC");
    return kFALSE;
  } else {
    LogEvent("MC");
  }
  
  fIsMC = kTRUE;
  
  // FIXME: using mc stack is outdated and not recommended!
  fMCStack = fMC->Stack();
  if (!fMCStack) {
    Err("noMCstack");
    fIsMC = kFALSE;
  }
  
  fMCGenHeader = nullptr;
  fMCHeader = fMC->Header();
  if (!fMCHeader) {
    Err("noMCHeader");
    fIsMC = kFALSE;
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
  InitMCEvent();
  return fIsMC;
}


//****************************************************************************************
/**
 * Initialise event-related properties, internal function.
 * Sets multiplicities, vertex, fEventSpecie, fRunNumber, fRunNumberString, fNTracksESD, fNTracksAcc, etc
 */
//****************************************************************************************
Bool_t AliAnalysisTaskMKBase::InitEvent()
{
  fMultSelection = static_cast<AliMultSelection*>(fEvent->FindListObject("MultSelection"));
  
  // apply event cuts
  fIsAcceptedAliEventCuts = fEventCuts.AcceptEvent(fEvent);
  fIsAcceptedAliEventCuts = fEventCuts.AcceptEvent(fEvent);
    
  //FIXME: in principle all of the following is not necessary if event will be rejected anyway...
  
  if(fNeedEventVertex)  InitEventVertex();
  if(fNeedEventQA)      InitEventQA(); // needs vertex info!
  if(fNeedEventCent)    InitEventCent();
  if(fNeedEventMult)    InitEventMult(); //FIXME: this requires renaming as it is actually centrality..
  // FIXME: rename fMultPercentileV0M to fCent everywhere as it represents not necessarily v0m centrlaity
  if(fNeedEventVZERO)   InitEventVZERO();

  
  fEventSpecie = fESD->GetEventSpecie();
  
  fRunNumber = fESD->GetRunNumber();
  fRunNumberString = "";
  fRunNumberString += fRunNumber;
  
  fTimeStamp = fESD->GetTimeStamp();
  fEventNumberInFile = fESD->GetEventNumberInFile();
  
  // FIXME: check if this is outdated!
  // this is needed for some esd track cuts, to be on the save side we call it here
  if (!TGeoGlobalMagField::Instance()->GetField()) {
    fESD->InitMagneticField();
  }
  fNTracksESD = fESD->GetNumberOfTracks();
  
  // loop over all tracks to get the accepted multiplicity
  // only if track cuts are set
  fNTracksAcc = 0;
  if (fESDtrackCutsM) {
    // FIXME: this has no effect ?!
    fInternalLoop = kTRUE;
    LoopOverAllTracks();
    fInternalLoop = kFALSE;
  } else {
    fNTracksAcc = fNTracksESD;
    Log("noAliESDTrackCutsM");
  }
  
  // FIXME: what is this??
  fNTracksAcc = 0;
  if (fESDtrackCutsM) {
    for (Int_t i = 0; i < fNTracksESD; i++) {
      AliESDtrack* track = static_cast<AliESDtrack*>(fESD->GetTrack(i));
      if (!track)
        continue;
      if (fESDtrackCutsM->AcceptTrack(track)) {
        ++fNTracksAcc;
      }
    }
  } else {
    fNTracksAcc = fNTracksESD;
    Log("noAliESDTrackCutsM");
  }
  
  return kTRUE;
}

//****************************************************************************************
/**
 * Function to check all triggers and log into trigger histogram.
 */
//****************************************************************************************
void AliAnalysisTaskMKBase::FillTrigHist(TH1D* h)
{
  if (fEventSelected & AliVEvent::kMB) {
    Log(h, "kMB");
  } // Minimum bias trigger in PbPb 2010-11
  if (fEventSelected & AliVEvent::kINT1) {
    Log(h, "kINT1");
  } // V0A | V0C | SPD minimum bias trigger
  if (fEventSelected & AliVEvent::kINT7) {
    Log(h, "kINT7");
  } // V0AND minimum bias trigger
  if (fEventSelected & AliVEvent::kMUON) {
    Log(h, "kMUON");
  } // Single muon trigger in pp2010-11, INT1 suite
  if (fEventSelected & AliVEvent::kHighMult) {
    Log(h, "kHighMult");
  } // High-multiplicity SPD trigger
  if (fEventSelected & AliVEvent::kHighMultSPD) {
    Log(h, "kHighMultSPD");
  } // High-multiplicity SPD trigger
  if (fEventSelected & AliVEvent::kEMC1) {
    Log(h, "kEMC1");
  } // EMCAL trigger in pp2011, INT1 suite
  if (fEventSelected & AliVEvent::kCINT5) {
    Log(h, "kCINT5");
  } // V0OR minimum bias trigger
  if (fEventSelected & AliVEvent::kINT5) {
    Log(h, "kINT5");
  } // V0OR minimum bias trigger
  if (fEventSelected & AliVEvent::kCMUS5) {
    Log(h, "kCMUS5");
  } // Single muon trigger, INT5 suite
  if (fEventSelected & AliVEvent::kMUSPB) {
    Log(h, "kMUSPB");
  } // Single muon trigger in PbPb 2011
  if (fEventSelected & AliVEvent::kINT7inMUON) {
    Log(h, "kINT7inMUON");
  } // INT7 in MUON or MUFAST cluster
  if (fEventSelected & AliVEvent::kMuonSingleHighPt7) {
    Log(h, "kMuonSingleHighPt7");
  } // Single muon high-pt, INT7 suite
  if (fEventSelected & AliVEvent::kMUSH7) {
    Log(h, "kMUSH7");
  } // Single muon high-pt, INT7 suite
  if (fEventSelected & AliVEvent::kMUSHPB) {
    Log(h, "kMUSHPB");
  } // Single muon high-pt in PbPb 2011
  if (fEventSelected & AliVEvent::kMuonLikeLowPt7) {
    Log(h, "kMuonLikeLowPt7");
  } // Like-sign dimuon low-pt, INT7 suite
  if (fEventSelected & AliVEvent::kMUL7) {
    Log(h, "kMUL7");
  } // Like-sign dimuon low-pt, INT7 suite
  if (fEventSelected & AliVEvent::kMuonLikePB) {
    Log(h, "kMuonLikePB");
  } // Like-sign dimuon low-pt in PbPb 2011
  if (fEventSelected & AliVEvent::kMuonUnlikeLowPt7) {
    Log(h, "kMuonUnlikeLowPt7");
  } // Unlike-sign dimuon low-pt, INT7 suite
  if (fEventSelected & AliVEvent::kMUU7) {
    Log(h, "kMUU7");
  } // Unlike-sign dimuon low-pt, INT7 suite
  if (fEventSelected & AliVEvent::kMuonUnlikePB) {
    Log(h, "kMuonUnlikePB");
  } // Unlike-sign dimuon low-pt in PbPb 2011
  if (fEventSelected & AliVEvent::kEMC7) {
    Log(h, "kEMC7");
  } // EMCAL/DCAL L0 trigger, INT7 suite
  if (fEventSelected & AliVEvent::kEMC8) {
    Log(h, "kEMC8");
  } // EMCAL/DCAL L0 trigger, INT8 suite
  if (fEventSelected & AliVEvent::kMUS7) {
    Log(h, "kMUS7");
  } // Single muon low-pt, INT7 suite
  if (fEventSelected & AliVEvent::kMuonSingleLowPt7) {
    Log(h, "kMuonSingleLowPt7");
  } // Single muon low-pt, INT7 suite
  if (fEventSelected & AliVEvent::kPHI1) {
    Log(h, "kPHI1");
  } // PHOS L0 trigger in pp2011, INT1 suite
  if (fEventSelected & AliVEvent::kPHI7) {
    Log(h, "kPHI7");
  } // PHOS trigger, INT7 suite
  if (fEventSelected & AliVEvent::kPHI8) {
    Log(h, "kPHI8");
  } // PHOS trigger, INT8 suite
  if (fEventSelected & AliVEvent::kPHOSPb) {
    Log(h, "kPHOSPb");
  } // PHOS trigger in PbPb 2011
  if (fEventSelected & AliVEvent::kEMCEJE) {
    Log(h, "kEMCEJE");
  } // EMCAL/DCAL L1 jet trigger
  if (fEventSelected & AliVEvent::kEMCEGA) {
    Log(h, "kEMCEGA");
  } // EMCAL/DCAL L1 gamma trigger
  if (fEventSelected & AliVEvent::kHighMultV0) {
    Log(h, "kHighMultV0");
  } // High-multiplicity V0 trigger
  if (fEventSelected & AliVEvent::kCentral) {
    Log(h, "kCentral");
  } // Central trigger in PbPb 2011
  if (fEventSelected & AliVEvent::kSemiCentral) {
    Log(h, "kSemiCentral");
  } // Semicentral trigger in PbPb 2011
  if (fEventSelected & AliVEvent::kDG) {
    Log(h, "kDG");
  } // Double gap diffractive
  if (fEventSelected & AliVEvent::kDG5) {
    Log(h, "kDG5");
  } // Double gap diffractive
  if (fEventSelected & AliVEvent::kZED) {
    Log(h, "kZED");
  } // ZDC electromagnetic dissociation
  if (fEventSelected & AliVEvent::kSPI7) {
    Log(h, "kSPI7");
  } // Power interaction trigger
  if (fEventSelected & AliVEvent::kSPI) {
    Log(h, "kSPI");
  } // Power interaction trigger
  if (fEventSelected & AliVEvent::kINT8) {
    Log(h, "kINT8");
  } // 0TVX trigger
  if (fEventSelected & AliVEvent::kMuonSingleLowPt8) {
    Log(h, "kMuonSingleLowPt8");
  } // Single muon low-pt, INT8 suite
  if (fEventSelected & AliVEvent::kMuonSingleHighPt8) {
    Log(h, "kMuonSingleHighPt8");
  } // Single muon high-pt, INT8 suite
  if (fEventSelected & AliVEvent::kMuonLikeLowPt8) {
    Log(h, "kMuonLikeLowPt8");
  } // Like-sign dimuon low-pt, INT8 suite
  if (fEventSelected & AliVEvent::kMuonUnlikeLowPt8) {
    Log(h, "kMuonUnlikeLowPt8");
  } // Unlike-sign dimuon low-pt, INT8 suite
  if (fEventSelected & AliVEvent::kMuonUnlikeLowPt0) {
    Log(h, "kMuonUnlikeLowPt0");
  } // Unlike-sign dimuon low-pt, no additional L0 requirement
  if (fEventSelected & AliVEvent::kUserDefined) {
    Log(h, "kUserDefined");
  } // Set when custom trigger classes are set in AliPhysicsSelection
  if (fEventSelected & AliVEvent::kTRD) {
    Log(h, "kTRD");
  } // TRD trigger
  if (fEventSelected & AliVEvent::kMuonCalo) {
    Log(h, "kMuonCalo");
  } // Muon-calo triggers
  if (fEventSelected & AliVEvent::kCaloOnly) {
    Log(h, "kCaloOnly");
  } // MB, EMCAL and PHOS triggers in CALO or CALOFAST cluster
  if (fEventSelected & AliVEvent::kFastOnly) {
    Log(h, "kFastOnly");
  } // The fast cluster fired. This bit is set in to addition another trigger
  // bit, e.g. kMB
  if (fEventSelected & AliVEvent::kAny) {
    Log(h, "kAny");
  } // to accept any defined trigger
  if (fEventSelected & AliVEvent::kAnyINT) {
    Log(h, "kAnyINT");
  } // to accept any interaction (aka minimum bias) trigger
}

//****************************************************************************************
/**
 * Base function for analysing MC particles.
 */
//****************************************************************************************
void AliAnalysisTaskMKBase::BaseAnaParticleMC(Int_t flag)
{
  // for non-interal loop call AnaParticleMC
  if (!fInternalLoop) {
    AnaParticleMC(flag);
    return;
  }
  
  // internal loop to get multiplicities
  if (fMCIsCharged && fMCisPrim && !fMCPileUpTrack) {
    fMCnPrim++;
    if (TMath::Abs(fMCEta) < 1.) {
      fMCnPrim10++;
    }
    if (TMath::Abs(fMCEta) < 0.8) {
      fMCnPrim08++;
    }
    if (TMath::Abs(fMCEta) < 0.5) {
      fMCnPrim05++;
    }
    if ((2.8 < fMCEta) && (fMCEta < 5.1)) {
      fMCnPrimV0M++;
    } // V0A
    if ((-3.7 < fMCEta) && (fMCEta < -1.7)) {
      fMCnPrimV0M++;
    } // V0C
    if ((TMath::Abs(fMCEta) < 0.8) && (fMCPt > 0.15)) {
      fMCnPrimPtCut++;
    }
  }
}

//****************************************************************************************
/**
 * Base function for analysing tracks.
 */
//****************************************************************************************
void AliAnalysisTaskMKBase::BaseAnaTrack(Int_t flag) {
  if (!fInternalLoop) {
    // for non-interal loop call AnaTrack
    AnaTrack(flag);
    if (fIsMC) {
      AnaTrackMC(flag);
    } else {
      AnaTrackDATA(flag);
    }
  } else {
    // internal loop to get multiplicities
    if (fAcceptTrackM) {
      fNTracksAcc++;
    }
  }
}

//****************************************************************************************
/**
 * Function for initializing MC events.
 */
//****************************************************************************************
Bool_t AliAnalysisTaskMKBase::InitMCEvent()
{
  // FIXME: this is almost never needed and should be optional!
  InitMCEventType();

  TArrayF vtxMC(3);
  // mc vertex
  fMCGenHeader->PrimaryVertex(vtxMC);
  fMCxv = vtxMC[0];
  fMCyv = vtxMC[1];
  fMCzv = vtxMC[2];
  
  // mc mult
  fMCnTracks = fMC->GetNumberOfTracks();
  
  fMCnPrimPtCut = 0;
  fMCnPrim10 = 0;
  fMCnPrim08 = 0;
  fMCnPrim05 = 0;
  fMCnPrimV0M = 0;
  fMCnPrim = 0;
  
  fInternalLoop = kTRUE;
  LoopOverAllParticles();
  fInternalLoop = kFALSE;
  
  return kTRUE;
}

//****************************************************************************************
/**
 * Determine type of MC event.
 */
//****************************************************************************************
Bool_t AliAnalysisTaskMKBase::InitMCEventType()
{
  AliGenPythiaEventHeader* fMCGenHeaderPythia =
  dynamic_cast<AliGenPythiaEventHeader*>(fMCGenHeader);
  AliGenDPMjetEventHeader* fMCGenHeaderDPMjet =
  dynamic_cast<AliGenDPMjetEventHeader*>(fMCGenHeader);
  
  fMCProcessTypeFlag = 0;
  fMCEventType = AlidNdPtTools::kInvalidProcess;
  
  if (fMCGenHeaderPythia) {
    LogEvent("PythiaGenHeader");
    fMCProcessTypeFlag = fMCGenHeaderPythia->ProcessType();
    if (fMCProcessTypeFlag == 92) {
      fMCEventType = AlidNdPtTools::kSD;
    } else if (fMCProcessTypeFlag == 93) {
      fMCEventType = AlidNdPtTools::kSD;
    } else if (fMCProcessTypeFlag == 94) {
      fMCEventType = AlidNdPtTools::kDD;
    } else if (fMCProcessTypeFlag == 91) {
      fMCEventType = AlidNdPtTools::kElastic;
    } else if (fMCProcessTypeFlag == -2) {
      fMCEventType = AlidNdPtTools::kCD;
    } else if (fMCProcessTypeFlag == -1) {
      fMCEventType = AlidNdPtTools::kND;
    } else {
      Err("UnknownPythiaEventType");
    }
  } else if (fMCGenHeaderDPMjet) {
    LogEvent("DPMjetGenHeader");
    fMCProcessTypeFlag = fMCGenHeaderDPMjet->ProcessType();
    if (fMCProcessTypeFlag == 5) {
      fMCEventType = AlidNdPtTools::kSD;
    } else if (fMCProcessTypeFlag == 6) {
      fMCEventType = AlidNdPtTools::kSD;
    } else if (fMCProcessTypeFlag == 7) {
      fMCEventType = AlidNdPtTools::kDD;
    } else if (fMCProcessTypeFlag == 2) {
      fMCEventType = AlidNdPtTools::kElastic;
    } else if (fMCProcessTypeFlag == 4) {
      fMCEventType = AlidNdPtTools::kCD;
    } else if (fMCProcessTypeFlag == 1) {
      fMCEventType = AlidNdPtTools::kND;
    } else {
      Err("UnknownDPMjetEventType");
    }
  } else {
    LogEvent("unknownGenHeader");
  }
  return kTRUE;
}

//****************************************************************************************
/**
 * Initialize V0 related variables.
 */
//****************************************************************************************
Bool_t AliAnalysisTaskMKBase::InitEventVZERO()
{
  fVZERO = fESD->GetVZEROData();
  if (!fVZERO) {
    Log("noVZEROData");
    fMultV0A = -1;
    fMultV0C = -1;
    fMultV0M = -1;
    return kFALSE;
  }
  fMultV0A = fVZERO->GetMTotV0A();
  fMultV0C = fVZERO->GetMTotV0C();
  fMultV0M = fMultV0A + fMultV0C;
  
  fMultV0MmultSelection = -1;
  if (fMultSelection) {
    AliMultEstimator* estv0m = fMultSelection->GetEstimator("V0M");
    fMultV0MmultSelection = estv0m->GetValue();
    TString defv0m = estv0m->GetDefinition();
    
    //         cout<<endl;
    //         cout<<"V0M DEFINTION"<<endl;
    //         cout<<defv0m<<endl;
    //         cout<<endl;
    //         cout<<"fAmplitude_V0A "<<fMultV0A<<endl;
    //         cout<<"fAmplitude_V0C "<<fMultV0C<<endl;
    //         cout<<"fEvSel_VtxZ    "<<fZv<<endl;
    //         cout<<endl;
    //         cout<<"fMultV0MmultSelection "<<fMultV0MmultSelection<<endl;
    
    defv0m.ReplaceAll("(fAmplitude_V0A)", "[0]");
    defv0m.ReplaceAll("(fAmplitude_V0C)", "[1]");
    defv0m.ReplaceAll("(fEvSel_VtxZ)", "[2]");
    TFormula formv0m("form", defv0m);
    formv0m.SetParameters(fMultV0A, fMultV0C, fZv);
    //         cout<<"my own calulation     "<<formv0m.Eval(0)<<endl;
    //         cout<<endl;
    //         cout<<endl;
    if (TMath::Abs(formv0m.Eval(0) - fMultV0MmultSelection) > 0.01) {
      Err("multSelection.V0Mmismatch");
    }
  }
  
  return kTRUE;
}

//****************************************************************************************
/**
 * Initialize track related variables.
 */
//****************************************************************************************
Bool_t AliAnalysisTaskMKBase::InitTrack()
{
  fPt = fESDTrack->Pt();
  fEta = fESDTrack->Eta();
  fPhi = fESDTrack->Phi();
  fP = fESDTrack->GetP();
  fChargeSign = fESDTrack->Charge();
  fESDTrack->GetImpactParameters(fDCA, fDCACov);
  fDCAr = fDCA[0];
  fDCAz = fDCA[1];
  fSigma1Pt2 = fESDTrack->GetSigma1Pt2();
  fTPCSignalN = fESDTrack->GetTPCsignalN();
  if(fSigma1Pt2 < 0) {
    Err("Sigma1Pt2<0");
  }
  fSigma1Pt = TMath::Sqrt(fSigma1Pt2);
  fSigned1Pt = fESDTrack->GetSigned1Pt();
  f1Pt = TMath::Abs(fSigned1Pt);
  
  InitTrackCuts();
  if(fNeedTrackQA)  InitTrackQA();  // CAVEAT: better apply this after basic cuts to reduce computation effort
  if(fNeedTrackIP)  InitTrackIP();
  if(fNeedTrackTPC) InitTrackTPC();
  if(fNeedTrackPID) InitTrackPID();
  if(fIsMC && !fSkipMCtruth) InitMCTrack();

  return kTRUE;
}

//****************************************************************************************
/**
 * Initialize track qa related variables.
 */
//****************************************************************************************
Bool_t AliAnalysisTaskMKBase::InitTrackQA()
{
  fX = fESDTrack->GetX();
  fY = fESDTrack->GetY();
  fZ = fESDTrack->GetZ();
  fAlpha = fESDTrack->GetAlpha();
  fSnp = fESDTrack->GetSnp();
  fTgl = fESDTrack->GetTgl();
  fITSFoundClusters = fESDTrack->GetITSclusters(0);
  fITSChi2PerCluster = (fITSFoundClusters > 0.) ? (fESDTrack->GetITSchi2()/fITSFoundClusters) : -1.;
  fITSClusterMap = fESDTrack->GetITSClusterMap();
  fFlags = fESDTrack->GetStatus();
  fTPCFindableClusters = fESDTrack->GetTPCNclsF();
  fTPCFoundClusters = fESDTrack->GetTPCclusters(0);
  fTPCChi2PerCluster = (fTPCFoundClusters > 0.) ? fESDTrack->GetTPCchi2()/fTPCFoundClusters : -1.;
  fTPCSharedClusters = fESDTrack->GetTPCnclsS();
  fTPCCrossedRows = fESDTrack->GetTPCCrossedRows();
  fTPCCrossedRowsOverFindableClusters = (fTPCFindableClusters>0) ? fTPCCrossedRows/fTPCFindableClusters : -1.;
  fTPCFractionSharedClusters = (fTPCFoundClusters > 0.) ? fTPCSharedClusters/fTPCFoundClusters : -1.;
  
  if(fESD->GetPrimaryVertex() && fESD->GetPrimaryVertex()->GetStatus())
  {
    fTPCGoldenChi2 = fESDTrack->GetChi2TPCConstrainedVsGlobal(fESD->GetPrimaryVertex());
  }

  const Double_t kDeadZoneWidth = 3.; // this is the default value used in the cut
  const Double_t kMaxZ = 220;
  if(fESDTrack->GetInnerParam()) fTPCGeomLength = fESDTrack->GetLengthInActiveZone(1, kDeadZoneWidth, kMaxZ, fESD->GetMagneticField());
  
  return kTRUE;
}
//****************************************************************************************
/**
 * Initialize track cuts.
 */
//****************************************************************************************
Bool_t AliAnalysisTaskMKBase::InitTrackCuts()
{
  fAcceptTrackM =
  (fESDtrackCutsM) ? fESDtrackCutsM->AcceptTrack(fESDTrack) : kFALSE;
  for (int i = 0; i < 30; i++) {
    fAcceptTrack[i] = (fESDtrackCuts[i])
    ? fESDtrackCuts[i]->AcceptTrack(fESDTrack)
    : kFALSE;
  }
  return kTRUE;
}

//****************************************************************************************
/**
 * Initialize trackPID.
 */
//****************************************************************************************
Bool_t AliAnalysisTaskMKBase::InitTrackPID()
{
  // TODO: implement track PID information
  return kTRUE;
}

//****************************************************************************************
/**
 * Initialize MC track related variables.
 */
//****************************************************************************************
Bool_t AliAnalysisTaskMKBase::InitMCTrack()
{
  fMCPrimSec = -1;
  
  // FIXME: remove thes obsolete checks!
  if (!fESDTrack) return kFALSE;
  if (!fMC) return kFALSE;
  
  fMCLabel = TMath::Abs(fESDTrack->GetLabel());
  if (fMCLabel < 0) {
    Log("tracklabel<0");
  }
  fMCParticle = dynamic_cast<AliMCParticle*>(fMC->GetTrack(fMCLabel));
  if (!fMCParticle) {
    Err("particleNOTinStack");
    return kFALSE;
  }
  InitMCParticle();
  if (fMCPrimSec == -1) {
    Err("TrackNOTprimORsec");
    return kFALSE;
  }
  return kTRUE;
}

//****************************************************************************************
/**
 * Initialize MC particle related variables.
 */
//****************************************************************************************
Bool_t AliAnalysisTaskMKBase::InitMCParticle()
{
  fMCPt = fMCParticle->Pt();
  fMCEta = fMCParticle->Eta();
  fMCPhi = fMCParticle->Phi();
  fMCCharge = fMCParticle->Charge();
  fMCQ = fMCCharge / 3.0;
  if (fMCCharge == 0) {
    fMCIsCharged = kFALSE;
    fMCChargeSign = 0;
  } else {
    fMCIsCharged = kTRUE;
    fMCChargeSign = (fMCCharge > 0) ? +1 : -1;
  }
  
  fMCisPrim = fMC->IsPhysicalPrimary(fMCLabel);
  
  Bool_t isPhysPrimStack = fMCStack->IsPhysicalPrimary(fMCLabel);
  if (fMCisPrim != isPhysPrimStack) {
    Log("PhysPrimMismatch");
  }
  
  fMCisSecDecay = fMC->IsSecondaryFromWeakDecay(fMCLabel);
  fMCisSecMat = fMC->IsSecondaryFromMaterial(fMCLabel);
  fMCisSec = fMCisSecMat || fMCisSecDecay;
  fMCPileUpTrack = AliAnalysisUtils::IsParticleFromOutOfBunchPileupCollision(fMCLabel, fMC);
  if (fMCisPrim) {
    fMCPrimSec = 0;
    fMCProdcutionType = AlidNdPtTools::kPrim;
  }
  else if (fMCisSecDecay) {
    fMCPrimSec = 1;
    fMCProdcutionType = AlidNdPtTools::kSecDecay;
  }
  else if (fMCisSecMat) {
    fMCPrimSec = 2;
    fMCProdcutionType = AlidNdPtTools::kSecMaterial;
  }
  else if(fMCPileUpTrack){
    fMCPrimSec = 3;
    fMCProdcutionType = AlidNdPtTools::ProductionType::kPileUpTrack;
  }
  else {
      fMCPrimSec = -1;
      fMCProdcutionType = AlidNdPtTools::ProductionType::kUnknown;
  }
  // if (fMCPrimSec == -1)             { Err("NOTprimORsec"); }
  if (fMCisPrim && fMCisSec) {
    Err("primANDsec");
  }
  if (fMCisSecDecay && fMCisSecMat) {
    Err("decayANDmat");
  }
  fMCPDGCode = fMCParticle->PdgCode();
  fMCParticleType = AlidNdPtTools::ParticleTypeFromPDG(fMCPDGCode);


  return kTRUE;
}

//****************************************************************************************
/**
 * Initialize event vertex.
 */
//****************************************************************************************
Bool_t AliAnalysisTaskMKBase::InitEventVertex()
{
  fVtxTRK = fESD->GetPrimaryVertexTracks();
  fVtxStatusTRK = kFALSE;
  if (fVtxTRK) {
    fVtxNContribTRK = fVtxTRK->GetNContributors();
    fXvTRK = fVtxTRK->GetX();
    fYvTRK = fVtxTRK->GetY();
    fZvTRK = fVtxTRK->GetZ();
    fXvResTRK = fVtxTRK->GetXRes();
    fYvResTRK = fVtxTRK->GetYRes();
    fZvResTRK = fVtxTRK->GetZRes();
    fVtxStatusTRK = fVtxTRK->GetStatus();
    fVtxDispersionTRK = fVtxTRK->GetDispersion();
    if (!fVtxStatusTRK) {
      LogEvent("noVertexTRK");
    }
  } else {
    LogEvent("PrimVtxTRK==0");
  }
  
  fVtxSPD = fESD->GetPrimaryVertexSPD();
  fVtxStatusSPD = kFALSE;
  if (fVtxSPD) {
    fVtxNContribSPD = fVtxSPD->GetNContributors();
    fXvSPD = fVtxSPD->GetX();
    fYvSPD = fVtxSPD->GetY();
    fZvSPD = fVtxSPD->GetZ();
    fXvResSPD = fVtxSPD->GetXRes();
    fYvResSPD = fVtxSPD->GetYRes();
    fZvResSPD = fVtxSPD->GetZRes();
    fVtxStatusSPD = fVtxSPD->GetStatus();
    fVtxDispersionSPD = fVtxSPD->GetDispersion();
    if (!fVtxStatusSPD) {
      LogEvent("noVertexSPD");
    }
  } else {
    LogEvent("PrimVtxSPD==0");
  }
  
  fVtxTPC = fESD->GetPrimaryVertexTPC();
  fVtxStatusTPC = kFALSE;
  if (fVtxTPC) {
    fVtxNContribTPC = fVtxTPC->GetNContributors();
    fXvTPC = fVtxTPC->GetX();
    fYvTPC = fVtxTPC->GetY();
    fZvTPC = fVtxTPC->GetZ();
    fXvResTPC = fVtxTPC->GetXRes();
    fYvResTPC = fVtxTPC->GetYRes();
    fZvResTPC = fVtxTPC->GetZRes();
    fVtxStatusTPC = fVtxTPC->GetStatus();
    fVtxDispersionTPC = fVtxTPC->GetDispersion();
    if (!fVtxStatusTPC) {
      LogEvent("noVertexTPC");
    }
  } else {
    LogEvent("PrimVtxTPC==0");
  }
  
  // get best available vertex according to AliESDevent
  fVtx = fESD->GetPrimaryVertex();
  fVtxStatus = kFALSE;
  if (fVtx) {
    fVtxNContrib = fVtx->GetNContributors();
    fXv = fVtx->GetX();
    fYv = fVtx->GetY();
    fZv = fVtx->GetZ();
    fXvRes = fVtx->GetXRes();
    fYvRes = fVtx->GetYRes();
    fZvRes = fVtx->GetZRes();
    fVtxStatus = fVtx->GetStatus();
    fVtxDispersion = fVtx->GetDispersion();
    if (!fVtxStatus) {
      LogEvent("noGoodVertex");
    }
  } else {
    LogEvent("PrimVtx==0");
  }
  
  fUsedVtxTRK = (fVtx == fVtxTRK);
  fUsedVtxSPD = (fVtx == fVtxSPD);
  fUsedVtxTPC = (fVtx == fVtxTPC);
  
  if (fUsedVtxTRK) {
    LogEvent("usedVertexTRK");
  }
  if (fUsedVtxSPD) {
    LogEvent("usedVertexSPD");
  }
  if (fUsedVtxTPC) {
    LogEvent("usedVertexTPC");
  }
  
  fMultMB = fVtxNContrib;
  fHasVertex = fVtxStatus; // fHasVertex is only an alias
  return kTRUE;
}

//****************************************************************************************
/**
 * Initialize event centrality using old centrality task (2011 PbPb 2.76 TeV data).
 */
//****************************************************************************************
Bool_t AliAnalysisTaskMKBase::InitEventCent()
{
  if (!fESD) {
    Err("MultCentRequiresESD");
    return kFALSE;
  }
  // FIXME: the centrality should be stored in standard centrality variable since
  // FIXME: both at the same time will never occur and then the code will always work when using fCentralityPercentile
  AliCentrality* fCentrality = fESD->GetCentrality();
  fOldCentPercentileV0M = -1;
  if (!fCentrality) {
    LogEvent("noOldCentrality");
    return kFALSE;
  } else {
    fOldCentPercentileV0M = fCentrality->GetCentralityPercentile("V0M");
    if (fOldCentPercentileV0M < 0.) {
      LogEvent("fOldCentPercentileV0M<0");
    }
    if (fOldCentPercentileV0M > 100.) {
      LogEvent("fOldCentPercentileV0M>100");
    }
    if (fOldCentPercentileV0M == 0.) {
      LogEvent("fOldCentPercentileV0M==0");
    }
    if (fOldCentPercentileV0M == 100.) {
      LogEvent("fOldCentPercentileV0M==100");
    }
  }
  return kTRUE;
}

//****************************************************************************************
/**
 * Initialize event centrality respectively multiplicity percentile in case of pp.
 */
//****************************************************************************************
Bool_t AliAnalysisTaskMKBase::InitEventMult()
{
  fMultPercentileV0M = -1;
  if (!fMultSelection) {
    LogEvent("noMultSelection");
    return kFALSE;
  } else
  { //FIXME: why determine this string per event?
    std::string estimator{"V0M"};
    switch (fCentralityEstimator) {
      case kV0M:
        //estimator = "V0M"; // already initialized whith this value
        break;
      case kCL0:
        estimator = "CL0";
        break;
      case kCL1:
        estimator = "CL1";
        break;
      case kV0Mplus05:
        estimator = "V0Mplus05";
        break;
      case kV0Mplus10:
        estimator = "V0Mplus10";
        break;
      case kV0Mminus05:
        estimator = "V0Mminus05";
        break;
      case kV0Mminus10:
        estimator = "V0Mminus10";
        break;
      case kSPDClustersCorr:
        estimator = "SPDClustersCorr";
        break;
      case kSPDTracklets:
        estimator = "SPDTracklets";
        break;
      default:
        break;
    }
    fMultPercentileV0M = fMultSelection->GetMultiplicityPercentile(estimator.data());
    
    // FIXME: this useless overhead should be avoided...
    /*
    LogEvent(Form("fMultPercentile_%s", estimator.data()));
    if (fMultPercentileV0M < 0.) {
      LogEvent("fMultPercentileV0M<0");
    }
    if (fMultPercentileV0M > 100.) {
      LogEvent("fMultPercentileV0M>100");
    }
    if (fMultPercentileV0M == 0.) {
      LogEvent("fMultPercentileV0M==0");
    }
    if (fMultPercentileV0M == 100.) {
      LogEvent("fMultPercentileV0M==100");
    }
     */
  }
  return kTRUE;
}

//****************************************************************************************
/**
 * Initialize interaction point related variables.
 */
//****************************************************************************************
Bool_t AliAnalysisTaskMKBase::InitTrackIP()
{
  fInnerP = fESDTrack->GetInnerParam();
  if (!fInnerP) {
    Err("noInnerParam");
    fPtInner = 0;
    fEtaInner = 0;
    fPhiInner = 0;
    fZInner = -100;
    return kFALSE;
  } else {
    fPtInner = fInnerP->Pt();
    fEtaInner = fInnerP->Eta();
    fPhiInner = fInnerP->Phi();
    fZInner = fInnerP->GetZ();
  }
  return kTRUE;
}

//****************************************************************************************
/**
 * Initialize TPC track related variables.
 */
//****************************************************************************************
Bool_t AliAnalysisTaskMKBase::InitTrackTPC()
{
  fESDTrack->GetImpactParametersTPC(fDCATPC, fDCACovTPC);
  fDCArTPC = fDCATPC[0];
  fDCAzTPC = fDCATPC[1];
  
  fTPCinnerP = fESDTrack->GetTPCInnerParam();
  if (!fTPCinnerP) {
    Err("noTPCinner");
    fPtInnerTPC = 0;
    fEtaInnerTPC = 0;
    fPhiInnerTPC = 0;
    return kFALSE;
  } else {
    fPtInnerTPC = fTPCinnerP->Pt();
    fEtaInnerTPC = fTPCinnerP->Eta();
    fPhiInnerTPC = fTPCinnerP->Phi();
    fZInnerTPC = fTPCinnerP->GetZ();
  }
  return kTRUE;
}

//****************************************************************************************
/**
 * Event loop called by the AnalysisManager.
 * Calls the event selection and AnaEvent functions of the derived classes.
 */
//****************************************************************************************
void AliAnalysisTaskMKBase::UserExec(Option_t*)
{
  // Read the event and set all event related properties
  if (!ReadEvent()) {
    Err("AliAnalysisTaskMKBase::UserExec:ErrorReadingEvent");
    return;
  }
  // we analyse only events that could be properly read
  if (fUseBaseOutput) FillDefaultHistograms(0);
  fIsEventAccepted = IsEventSelected();
  if (fIsEventAccepted)
  {
    if (fUseBaseOutput) FillDefaultHistograms(1);
    // call user analysis of the event
    AnaEvent();
    // call mc and data anlayses
    if (fIsMC && !fSkipMCtruth) {
      AnaEventMC();
    } else {
      AnaEventDATA();
    }
  }
  // postdata
  PostData(1, fOutputList);
}

//****************************************************************************************
/**
 * Function to fill default event histograms of the base class before and after the user implemented event selection
 * To disable filling of the histograms, overwrite this function in the derived class.
 */
//****************************************************************************************
void AliAnalysisTaskMKBase::FillDefaultHistograms(Int_t step)
{
  if (step == 0) {
    Log(fTrigInfo, fFiredTriggerClasses.Data());
    FillTrigHist(fTrigHist);
    Log(fRunHist, fRunNumberString.Data());
  } else if (step == 1) {
    Log(fTrigInfoSelected, fFiredTriggerClasses.Data());
    FillTrigHist(fTrigHistSelected);
    Log(fRunHistSelected, fRunNumberString.Data());
  } else {
    Err("AliAnalysisTaskMKBase::FillDefaultHistograms:InvalidStep");
  }
}


//****************************************************************************************
/**
 * Loop over all tracks.
 */
//****************************************************************************************
void AliAnalysisTaskMKBase::LoopOverAllTracks(Int_t flag)
{
  fNTracksESD = fESD->GetNumberOfTracks();
  for (Int_t i = 0; i < fNTracksESD; i++) {
    fESDTrack = dynamic_cast<AliESDtrack*>(fESD->GetTrack(i));
    if (!fESDTrack) {
      Err("noESDtrack");
      continue;
    }
    InitTrack();
    BaseAnaTrack(flag);
  }
}

//****************************************************************************************
/**
 * Loop over all particles.
 */
//****************************************************************************************
void AliAnalysisTaskMKBase::LoopOverAllParticles(Int_t flag)
{
  fMCnTracks = fMC->GetNumberOfTracks();
  for (Int_t i = 0; i < fMCnTracks; i++) {
    fMCParticle = dynamic_cast<AliMCParticle*>(fMC->GetTrack(i));
    if (!fMCParticle) {
      Err("noMCParticle");
      continue;
    }
    fMCLabel = i;
    InitMCParticle();
    BaseAnaParticleMC(flag);
  }
}

//****************************************************************************************
/**
 * Terminating function called by the manager after all events were processed.
 */
//****************************************************************************************
void AliAnalysisTaskMKBase::Terminate(Option_t*)
{
}

//****************************************************************************************
/**
 * Static factory function that serves as AddTask macro.
 * Creates and configures a task and attaches it to the AnalysisManager.
 * \return the created task
 */
//****************************************************************************************
AliAnalysisTaskMKBase* AliAnalysisTaskMKBase::AddTaskMKBase(const char* name, const char* outfile)
{
  AliAnalysisManager* mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskMKBase", "No analysis manager to connect to.");
    return 0;
  }

  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTaskMKBase", "This task requires an input event handler");
    return NULL;
  }
  
  // Setup output file
  TString fileName = AliAnalysisManager::GetCommonFileName();
  fileName += ":";
  fileName += name; // create a subfolder in the file
  if (outfile) {    // if a finename is given, use that one
    fileName = TString(outfile);
  }
  
  // create the task
  AliAnalysisTaskMKBase* task = new AliAnalysisTaskMKBase(name);
  if (!task) {
    return 0;
  }
  
  // configure the task
  task->SelectCollisionCandidates(AliVEvent::kAnyINT);
  task->SetESDtrackCutsM(AlidNdPtTools::CreateESDtrackCuts("defaultEta08"));
  
  // attach the task to the manager and configure in and ouput
  mgr->AddTask(task);
  mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(
                     task, 1,
                     mgr->CreateContainer(name, TList::Class(),
                                          AliAnalysisManager::kOutputContainer,
                                          fileName.Data()));
  
  return task;
}
