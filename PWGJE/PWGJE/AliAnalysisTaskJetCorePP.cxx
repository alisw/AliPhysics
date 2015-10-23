
// ******************************************
// This task computes several jet observables like 
// the fraction of energy in inner and outer coronnas,
// jet-track correlations,triggered jet shapes and 
// correlation strength distribution of particles inside jets.    
// Author: lcunquei@cern.ch
// *******************************************


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


#include "TChain.h"
#include "TClonesArray.h"
#include "TTree.h"
#include "TList.h"
#include "TMath.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TH2F.h"
#include "TH3F.h"
#include "THnSparse.h"
#include "TCanvas.h"
#include "TArrayI.h" 
#include "TProfile.h"
#include "TFile.h"
#include "TKey.h"
#include "TRandom3.h"

#include "AliLog.h"

#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"

#include "AliVEvent.h"
#include "AliESDEvent.h"
#include "AliMCEvent.h"
#include "AliESDInputHandler.h"
#include "AliCentrality.h"
#include "AliAnalysisHelperJetTasks.h"
#include "AliInputEventHandler.h"
#include "AliAODJetEventBackground.h"
#include "AliGenPythiaEventHeader.h"
#include "AliAODMCParticle.h"
#include "AliMCParticle.h"
#include "AliAODEvent.h"
#include "AliAODHandler.h"
#include "AliAODJet.h"
#include "AliVVertex.h"
#include "AliVTrack.h"
#include "AliAnalysisTaskJetCorePP.h"
#include "AliHeader.h" //KF//
#include <iostream>

using std::cout;
using std::endl;

ClassImp(AliAnalysisTaskJetCorePP)

//Filip Krizek 1st March 2013

//---------------------------------------------------------------------
AliAnalysisTaskJetCorePP::AliAnalysisTaskJetCorePP() :
AliAnalysisTaskSE(),
fESD(0x0),
fAODIn(0x0),
fAODOut(0x0),
fAODExtension(0x0),
fMcEvent(0x0),
fMcHandler(0x0),
fJetBranchName(""),
fJetBranchNameChargMC(""),
fJetBranchNameKine(""),
fJetBranchNameFullMC(""),
fJetBranchNameBg(""),
fJetBranchNameBgChargMC(""),
fJetBranchNameBgKine(""),
fListJets(0x0),
fListJetsGen(0x0),
fListJetsGenFull(0x0),
fListJetsBg(0x0),
fListJetsBgGen(0x0),
fNonStdFile(""),
fSystem(0), //pp=0  pPb=1
fJetParamR(0.4),
fBgJetParamR(0.3),
fBgMaxJetPt(8.0),
fBgConeR(0.4),
fOfflineTrgMask(AliVEvent::kAny),
fMinContribVtx(1),
fVtxZMin(-10.0),
fVtxZMax(10.0),
fFilterMask(0),
fCentMin(0.0),
fCentMax(100.0),
fJetEtaMin(-0.5),
fJetEtaMax(0.5),
fTriggerEtaCut(0.9),
fTrackEtaCut(0.9),
fTrackLowPtCut(0.15),
fUseExchContainer(0),
fOutputList(0x0),
fHistEvtSelection(0x0),
fh2Ntriggers(0x0),
fHJetSpec(0x0),
fHJetSpecSubUeMedian(0x0),
fHJetSpecSubUeCone(0x0),
fHJetPhiCorr(0x0),
fHJetPhiCorrSubUeMedian(0x0),
fHJetPhiCorrSubUeCone(0x0),
fHJetUeMedian(0x0),
fHJetUeCone(0x0),
fHRhoUeMedianVsCone(0x0),
//fHJetDensity(0x0),
//fHJetDensityA4(0x0),
fhJetPhi(0x0),
fhTriggerPhi(0x0),
fhJetEta(0x0),
fhTriggerEta(0x0),
fhVertexZ(0x0),
fhVertexZAccept(0x0),
fhContribVtx(0x0),
fhContribVtxAccept(0x0),
fhDphiTriggerJet(0x0),
fhDphiTriggerJetAccept(0x0),
fhCentrality(0x0),
fhCentralityAccept(0x0),
fhNofMultipleTriggers(0x0),
fhNofMultipleTriggersCone(0x0),
fhNofMultipleTriggersConeLow(0x0),
fhNofMultipleTriggersConeHigh(0x0),
fhDeltaRMultTriggersLow(0x0),
fhDeltaRMultTriggersHigh(0x0),
fhDeltaPhiMultTriggersLow(0x0),
fhDeltaPhiMultTriggersHigh(0x0),
fhDeltaPhiMultTriggersInclLow(0x0),
fhDeltaPhiMultTriggersInclHigh(0x0),
fhInclTrigCounter(0x0),
fhDeltaPtConeBg(0x0),
fhDeltaPtMedianBg(0x0),
//fHJetPtRaw(0x0),
//fHLeadingJetPtRaw(0x0), 
//fHDphiVsJetPtAll(0x0), 
fhJetPtGenVsJetPtRec(0x0),
fhJetPtGenVsJetPtRecSubUeMedian(0x0),
fhJetPtGenVsJetPtRecSubUeCone(0x0),
fhJetPtGen(0x0), 
fhJetPtSubUeMedianGen(0x0), 
fhJetPtSubUeConeGen(0x0),
fhJetPtResolutionVsPtGen(0x0),
fhJetPtResolutionVsPtConeGen(0x0), 
fhJetPtGenChargVsJetPtGenFull(0x0),
fhJetPtGenFull(0x0),
fh2NtriggersGen(0x0),
fHJetSpecGen(0x0),
fHJetSpecSubUeMedianGen(0x0),
fHJetSpecSubUeConeGen(0x0),
fHJetPhiCorrGen(0x0),
fHJetPhiCorrSubUeMedianGen(0x0),
fHJetPhiCorrSubUeConeGen(0x0),
fHJetUeMedianGen(0x0),
fHJetUeConeGen(0x0),
fhPtTrkTruePrimRec(0x0),
fhPtTrkTruePrimGen(0x0),
fhPtTrkSecOrFakeRec(0x0),
fHRhoUeMedianVsConeGen(0x0),
fhEntriesToMedian(0x0),
fhEntriesToMedianGen(0x0),
fhCellAreaToMedian(0x0),
fhCellAreaToMedianGen(0x0),
fhNofMultipleTriggersGen(0x0),
fhNofMultipleTriggersConeGen(0x0),
fhNofMultipleTriggersConeGenLow(0x0),
fhNofMultipleTriggersConeGenHigh(0x0),
fhDeltaRMultTriggersGenLow(0x0),
fhDeltaRMultTriggersGenHigh(0x0),
fhDeltaPhiMultTriggersGenLow(0x0),
fhDeltaPhiMultTriggersGenHigh(0x0),
fhDeltaPhiMultTriggersInclGenLow(0x0),
fhDeltaPhiMultTriggersInclGenHigh(0x0),
fhInclTrigCounterGen(0x0),
fhNofMultipleTriggersConeGenA(0x0),
fhNofMultipleTriggersConeGenALow(0x0),
fhNofMultipleTriggersConeGenAHigh(0x0),
fhDeltaRMultTriggersGenALow(0x0),
fhDeltaPhiMultTriggersGenALow(0x0),
fhDeltaRMultTriggersGenAHigh(0x0),
fhDeltaPhiMultTriggersGenAHigh(0x0),
fhNofTriggersGenA(0x0),
fhNofMultipleTriggersConeGenB(0x0),
fhNofMultipleTriggersConeGenBLow(0x0),
fhNofMultipleTriggersConeGenBHigh(0x0),
fhDeltaRMultTriggersGenBLow(0x0),
fhDeltaPhiMultTriggersGenBLow(0x0),
fhDeltaRMultTriggersGenBHigh(0x0),
fhDeltaPhiMultTriggersGenBHigh(0x0),
fhNofTriggersGenB(0x0),
fhTriggerCounterGenLevel(0x0),
fhDeltaRMultTriggersGenLevelLow(0x0),
fhDeltaPhiMultTriggersGenLevelLow(0x0),
fhDeltaRMultTriggersGenLevelHigh(0x0),
fhDeltaPhiMultTriggersGenLevelHigh(0x0),
fIsChargedMC(0),
fIsKine(0),
fIsFullMC(0),
faGenIndex(0),
faRecIndex(0),
fkAcceptance(2.0*TMath::Pi()*1.8),
fkDeltaPhiCut(TMath::Pi()-0.8),
fh1Xsec(0x0),
fh1Trials(0x0),
fh1AvgTrials(0x0),
fh1PtHard(0x0),
fh1PtHardNoW(0x0),  
fh1PtHardTrials(0x0),
fAvgTrials(1),
fHardest(0),
fEventNumberRangeLow(0),
fEventNumberRangeHigh(99),
fTriggerPtRangeLow(0.0),
fTriggerPtRangeHigh(10000.0),
fFillRespMx(0),
fRandom(0x0),
fnTrials(1000),
fJetFreeAreaFrac(0.5),
fnPhi(9),
fnEta(2),
fEtaSize(0.7),
fPhiSize(2*TMath::Pi()/fnPhi),
fCellArea(fPhiSize*fEtaSize),
fSafetyMargin(1.1),
fDoubleBinning(kFALSE)
{
   // default Constructor
}

//---------------------------------------------------------------------

AliAnalysisTaskJetCorePP::AliAnalysisTaskJetCorePP(const char *name) :
AliAnalysisTaskSE(name),
fESD(0x0),
fAODIn(0x0),
fAODOut(0x0),
fAODExtension(0x0),
fMcEvent(0x0),
fMcHandler(0x0),
fJetBranchName(""),
fJetBranchNameChargMC(""),
fJetBranchNameKine(""),
fJetBranchNameFullMC(""),
fJetBranchNameBg(""),
fJetBranchNameBgChargMC(""),
fJetBranchNameBgKine(""),
fListJets(0x0),
fListJetsGen(0x0),
fListJetsGenFull(0x0),
fListJetsBg(0x0),
fListJetsBgGen(0x0),
fNonStdFile(""),
fSystem(0),  //pp=0   pPb=1
fJetParamR(0.4),
fBgJetParamR(0.3),
fBgMaxJetPt(8.0),
fBgConeR(0.4),
fOfflineTrgMask(AliVEvent::kAny),
fMinContribVtx(1),
fVtxZMin(-10.0),
fVtxZMax(10.0),
fFilterMask(0),
fCentMin(0.0),
fCentMax(100.0),
fJetEtaMin(-0.5),
fJetEtaMax(0.5),
fTriggerEtaCut(0.9),
fTrackEtaCut(0.9),
fTrackLowPtCut(0.15),
fUseExchContainer(0),
fOutputList(0x0),
fHistEvtSelection(0x0),
fh2Ntriggers(0x0),
fHJetSpec(0x0),
fHJetSpecSubUeMedian(0x0),
fHJetSpecSubUeCone(0x0),
fHJetPhiCorr(0x0),
fHJetPhiCorrSubUeMedian(0x0),
fHJetPhiCorrSubUeCone(0x0),
fHJetUeMedian(0x0),
fHJetUeCone(0x0),
fHRhoUeMedianVsCone(0x0),
//fHJetDensity(0x0),
//fHJetDensityA4(0x0),
fhJetPhi(0x0),
fhTriggerPhi(0x0),
fhJetEta(0x0),
fhTriggerEta(0x0),
fhVertexZ(0x0),
fhVertexZAccept(0x0),
fhContribVtx(0x0),
fhContribVtxAccept(0x0),
fhDphiTriggerJet(0x0),
fhDphiTriggerJetAccept(0x0),
fhCentrality(0x0),
fhCentralityAccept(0x0),
fhNofMultipleTriggers(0x0),
fhNofMultipleTriggersCone(0x0),
fhNofMultipleTriggersConeLow(0x0),
fhNofMultipleTriggersConeHigh(0x0),
fhDeltaRMultTriggersLow(0x0),
fhDeltaRMultTriggersHigh(0x0),
fhDeltaPhiMultTriggersLow(0x0),
fhDeltaPhiMultTriggersHigh(0x0),
fhDeltaPhiMultTriggersInclLow(0x0),
fhDeltaPhiMultTriggersInclHigh(0x0),
fhInclTrigCounter(0x0),
fhDeltaPtConeBg(0x0),
fhDeltaPtMedianBg(0x0),
//fHJetPtRaw(0x0),
//fHLeadingJetPtRaw(0x0), 
//fHDphiVsJetPtAll(0x0), 
fhJetPtGenVsJetPtRec(0x0),
fhJetPtGenVsJetPtRecSubUeMedian(0x0),
fhJetPtGenVsJetPtRecSubUeCone(0x0),
fhJetPtGen(0x0),
fhJetPtSubUeMedianGen(0x0), 
fhJetPtSubUeConeGen(0x0), 
fhJetPtResolutionVsPtGen(0x0),
fhJetPtResolutionVsPtConeGen(0x0),
fhJetPtGenChargVsJetPtGenFull(0x0),
fhJetPtGenFull(0x0),
fh2NtriggersGen(0x0),
fHJetSpecGen(0x0),
fHJetSpecSubUeMedianGen(0x0),
fHJetSpecSubUeConeGen(0x0),
fHJetPhiCorrGen(0x0),
fHJetPhiCorrSubUeMedianGen(0x0),
fHJetPhiCorrSubUeConeGen(0x0),
fHJetUeMedianGen(0x0),
fHJetUeConeGen(0x0),
fhPtTrkTruePrimRec(0x0),
fhPtTrkTruePrimGen(0x0),
fhPtTrkSecOrFakeRec(0x0),
fHRhoUeMedianVsConeGen(0x0),
fhEntriesToMedian(0x0),
fhEntriesToMedianGen(0x0),
fhCellAreaToMedian(0x0),
fhCellAreaToMedianGen(0x0),
fhNofMultipleTriggersGen(0x0),
fhNofMultipleTriggersConeGen(0x0),
fhNofMultipleTriggersConeGenLow(0x0),
fhNofMultipleTriggersConeGenHigh(0x0),
fhDeltaRMultTriggersGenLow(0x0),
fhDeltaRMultTriggersGenHigh(0x0),
fhDeltaPhiMultTriggersGenLow(0x0),
fhDeltaPhiMultTriggersGenHigh(0x0),
fhDeltaPhiMultTriggersInclGenLow(0x0),
fhDeltaPhiMultTriggersInclGenHigh(0x0),
fhInclTrigCounterGen(0x0),
fhNofMultipleTriggersConeGenA(0x0),
fhNofMultipleTriggersConeGenALow(0x0),
fhNofMultipleTriggersConeGenAHigh(0x0),
fhDeltaRMultTriggersGenALow(0x0),
fhDeltaPhiMultTriggersGenALow(0x0),
fhDeltaRMultTriggersGenAHigh(0x0),
fhDeltaPhiMultTriggersGenAHigh(0x0),
fhNofTriggersGenA(0x0),
fhNofMultipleTriggersConeGenB(0x0),
fhNofMultipleTriggersConeGenBLow(0x0),
fhNofMultipleTriggersConeGenBHigh(0x0),
fhDeltaRMultTriggersGenBLow(0x0),
fhDeltaPhiMultTriggersGenBLow(0x0),
fhDeltaRMultTriggersGenBHigh(0x0),
fhDeltaPhiMultTriggersGenBHigh(0x0),
fhNofTriggersGenB(0x0),
fhTriggerCounterGenLevel(0x0),
fhDeltaRMultTriggersGenLevelLow(0x0),
fhDeltaPhiMultTriggersGenLevelLow(0x0),
fhDeltaRMultTriggersGenLevelHigh(0x0),
fhDeltaPhiMultTriggersGenLevelHigh(0x0),
fIsChargedMC(0),
fIsKine(0),
fIsFullMC(0),
faGenIndex(0),
faRecIndex(0),
fkAcceptance(2.0*TMath::Pi()*1.8),
fkDeltaPhiCut(TMath::Pi()-0.8),
fh1Xsec(0x0),
fh1Trials(0x0), 
fh1AvgTrials(0x0),
fh1PtHard(0x0),
fh1PtHardNoW(0x0),  
fh1PtHardTrials(0x0),
fAvgTrials(1),
fHardest(0),
fEventNumberRangeLow(0),
fEventNumberRangeHigh(99),
fTriggerPtRangeLow(0.0),
fTriggerPtRangeHigh(10000.0),
fFillRespMx(0),
fRandom(0x0),
fnTrials(1000),
fJetFreeAreaFrac(0.5),
fnPhi(9),
fnEta(2),
fEtaSize(0.7),
fPhiSize(2*TMath::Pi()/fnPhi),
fCellArea(fPhiSize*fEtaSize),
fSafetyMargin(1.1),
fDoubleBinning(kFALSE)
{
// Constructor
   DefineOutput(1, TList::Class());

   TString dummy(name);
   if(dummy.Contains("KINE")){
      DefineInput(1, TClonesArray::Class());
      DefineInput(2, TClonesArray::Class());
   }
}

//--------------------------------------------------------------
AliAnalysisTaskJetCorePP::AliAnalysisTaskJetCorePP(const AliAnalysisTaskJetCorePP& a):
AliAnalysisTaskSE(a.GetName()),
fESD(a.fESD),
fAODIn(a.fAODIn),
fAODOut(a.fAODOut),
fAODExtension(a.fAODExtension),
fMcEvent(a.fMcEvent),
fMcHandler(a.fMcHandler),
fJetBranchName(a.fJetBranchName),
fJetBranchNameChargMC(a.fJetBranchNameChargMC),
fJetBranchNameKine(a.fJetBranchNameKine),
fJetBranchNameFullMC(a.fJetBranchNameFullMC),
fJetBranchNameBg(a.fJetBranchNameBg),
fJetBranchNameBgChargMC(a.fJetBranchNameBgChargMC),
fJetBranchNameBgKine(a.fJetBranchNameBgKine),
fListJets(a.fListJets),
fListJetsGen(a.fListJetsGen),
fListJetsGenFull(a.fListJetsGenFull),
fListJetsBg(a.fListJetsBg),
fListJetsBgGen(a.fListJetsBgGen),
fNonStdFile(a.fNonStdFile),
fSystem(a.fSystem),  
fJetParamR(a.fJetParamR),
fBgJetParamR(a.fBgJetParamR),
fBgMaxJetPt(a.fBgMaxJetPt),
fBgConeR(a.fBgConeR),
fOfflineTrgMask(a.fOfflineTrgMask),
fMinContribVtx(a.fMinContribVtx),
fVtxZMin(a.fVtxZMin),
fVtxZMax(a.fVtxZMax),
fFilterMask(a.fFilterMask),
fCentMin(a.fCentMin),
fCentMax(a.fCentMax),
fJetEtaMin(a.fJetEtaMin),
fJetEtaMax(a.fJetEtaMax),
fTriggerEtaCut(a.fTriggerEtaCut),
fTrackEtaCut(a.fTrackEtaCut),
fTrackLowPtCut(a.fTrackLowPtCut),
fUseExchContainer(a.fUseExchContainer),
fOutputList(a.fOutputList),
fHistEvtSelection(a.fHistEvtSelection),
fh2Ntriggers(a.fh2Ntriggers),
fHJetSpec(a.fHJetSpec),
fHJetSpecSubUeMedian(a.fHJetSpecSubUeMedian),
fHJetSpecSubUeCone(a.fHJetSpecSubUeCone),
fHJetPhiCorr(a.fHJetPhiCorr),
fHJetPhiCorrSubUeMedian(a.fHJetPhiCorrSubUeMedian),
fHJetPhiCorrSubUeCone(a.fHJetPhiCorrSubUeCone),
fHJetUeMedian(a.fHJetUeMedian),
fHJetUeCone(a.fHJetUeCone),
fHRhoUeMedianVsCone(a.fHRhoUeMedianVsCone), 
//fHJetDensity(a.fHJetDensity),
//fHJetDensityA4(a.fHJetDensityA4),
fhJetPhi(a.fhJetPhi),
fhTriggerPhi(a.fhTriggerPhi),
fhJetEta(a.fhJetEta),
fhTriggerEta(a.fhTriggerEta),
fhVertexZ(a.fhVertexZ),
fhVertexZAccept(a.fhVertexZAccept),
fhContribVtx(a.fhContribVtx),
fhContribVtxAccept(a.fhContribVtxAccept),
fhDphiTriggerJet(a.fhDphiTriggerJet),
fhDphiTriggerJetAccept(a.fhDphiTriggerJetAccept),
fhCentrality(a.fhCentrality),
fhCentralityAccept(a.fhCentralityAccept),
fhNofMultipleTriggers(a.fhNofMultipleTriggers),
fhNofMultipleTriggersCone(a.fhNofMultipleTriggersCone),
fhNofMultipleTriggersConeLow(a.fhNofMultipleTriggersConeLow),
fhNofMultipleTriggersConeHigh(a.fhNofMultipleTriggersConeHigh),
fhDeltaRMultTriggersLow(a.fhDeltaRMultTriggersLow),
fhDeltaRMultTriggersHigh(a.fhDeltaRMultTriggersHigh),
fhDeltaPhiMultTriggersLow(a.fhDeltaPhiMultTriggersLow),
fhDeltaPhiMultTriggersHigh(a.fhDeltaPhiMultTriggersHigh),
fhDeltaPhiMultTriggersInclLow(a.fhDeltaPhiMultTriggersInclLow),
fhDeltaPhiMultTriggersInclHigh(a.fhDeltaPhiMultTriggersInclHigh),
fhInclTrigCounter(a.fhInclTrigCounter),
fhDeltaPtConeBg(a.fhDeltaPtConeBg),
fhDeltaPtMedianBg(a.fhDeltaPtMedianBg),
//fHJetPtRaw(a.fHJetPtRaw),
//fHLeadingJetPtRaw(a.fHLeadingJetPtRaw),
//fHDphiVsJetPtAll(a.fHDphiVsJetPtAll),
fhJetPtGenVsJetPtRec(a.fhJetPtGenVsJetPtRec),
fhJetPtGenVsJetPtRecSubUeMedian(a.fhJetPtGenVsJetPtRecSubUeMedian),
fhJetPtGenVsJetPtRecSubUeCone(a.fhJetPtGenVsJetPtRecSubUeCone),
fhJetPtGen(a.fhJetPtGen),
fhJetPtSubUeMedianGen(a.fhJetPtSubUeMedianGen), 
fhJetPtSubUeConeGen(a.fhJetPtSubUeConeGen), 
fhJetPtResolutionVsPtGen(a.fhJetPtResolutionVsPtGen),
fhJetPtResolutionVsPtConeGen(a.fhJetPtResolutionVsPtConeGen),
fhJetPtGenChargVsJetPtGenFull(a.fhJetPtGenChargVsJetPtGenFull),
fhJetPtGenFull(a.fhJetPtGenFull),
fh2NtriggersGen(a.fh2NtriggersGen),
fHJetSpecGen(a.fHJetSpecGen),
fHJetSpecSubUeMedianGen(a.fHJetSpecSubUeMedianGen),
fHJetSpecSubUeConeGen(a.fHJetSpecSubUeConeGen),
fHJetPhiCorrGen(a.fHJetPhiCorrGen),
fHJetPhiCorrSubUeMedianGen(a.fHJetPhiCorrSubUeMedianGen),
fHJetPhiCorrSubUeConeGen(a.fHJetPhiCorrSubUeConeGen),
fHJetUeMedianGen(a.fHJetUeMedianGen),
fHJetUeConeGen(a.fHJetUeConeGen),
fhPtTrkTruePrimRec(a.fhPtTrkTruePrimRec),
fhPtTrkTruePrimGen(a.fhPtTrkTruePrimGen),
fhPtTrkSecOrFakeRec(a.fhPtTrkSecOrFakeRec),
fHRhoUeMedianVsConeGen(a.fHRhoUeMedianVsConeGen),
fhEntriesToMedian(a.fhEntriesToMedian),
fhEntriesToMedianGen(a.fhEntriesToMedianGen),
fhCellAreaToMedian(a.fhCellAreaToMedian),
fhCellAreaToMedianGen(a.fhCellAreaToMedianGen),
fhNofMultipleTriggersGen(a.fhNofMultipleTriggersGen),
fhNofMultipleTriggersConeGen(a.fhNofMultipleTriggersConeGen),
fhNofMultipleTriggersConeGenLow(a.fhNofMultipleTriggersConeGenLow),
fhNofMultipleTriggersConeGenHigh(a.fhNofMultipleTriggersConeGenHigh),
fhDeltaRMultTriggersGenLow(a.fhDeltaRMultTriggersGenLow),
fhDeltaRMultTriggersGenHigh(a.fhDeltaRMultTriggersGenHigh),
fhDeltaPhiMultTriggersGenLow(a.fhDeltaPhiMultTriggersGenLow),
fhDeltaPhiMultTriggersGenHigh(a.fhDeltaPhiMultTriggersGenHigh),
fhDeltaPhiMultTriggersInclGenLow(a.fhDeltaPhiMultTriggersInclGenLow),
fhDeltaPhiMultTriggersInclGenHigh(a.fhDeltaPhiMultTriggersInclGenHigh),
fhInclTrigCounterGen(a.fhInclTrigCounterGen),
fhNofMultipleTriggersConeGenA(a.fhNofMultipleTriggersConeGenA),
fhNofMultipleTriggersConeGenALow(a.fhNofMultipleTriggersConeGenALow),
fhNofMultipleTriggersConeGenAHigh(a.fhNofMultipleTriggersConeGenAHigh),
fhDeltaRMultTriggersGenALow(a.fhDeltaRMultTriggersGenALow),
fhDeltaPhiMultTriggersGenALow(a.fhDeltaPhiMultTriggersGenALow),
fhDeltaRMultTriggersGenAHigh(a.fhDeltaRMultTriggersGenAHigh),
fhDeltaPhiMultTriggersGenAHigh(a.fhDeltaPhiMultTriggersGenAHigh),
fhNofTriggersGenA(a.fhNofTriggersGenA),
fhNofMultipleTriggersConeGenB(a.fhNofMultipleTriggersConeGenB),
fhNofMultipleTriggersConeGenBLow(a.fhNofMultipleTriggersConeGenBLow),
fhNofMultipleTriggersConeGenBHigh(a.fhNofMultipleTriggersConeGenBHigh),
fhDeltaRMultTriggersGenBLow(a.fhDeltaRMultTriggersGenBLow),
fhDeltaPhiMultTriggersGenBLow(a.fhDeltaPhiMultTriggersGenBLow),
fhDeltaRMultTriggersGenBHigh(a.fhDeltaRMultTriggersGenBHigh),
fhDeltaPhiMultTriggersGenBHigh(a.fhDeltaPhiMultTriggersGenBHigh),
fhNofTriggersGenB(a.fhNofTriggersGenB),
fhTriggerCounterGenLevel(a.fhTriggerCounterGenLevel),
fhDeltaRMultTriggersGenLevelLow(a.fhDeltaRMultTriggersGenLevelLow),
fhDeltaPhiMultTriggersGenLevelLow(a.fhDeltaPhiMultTriggersGenLevelLow),
fhDeltaRMultTriggersGenLevelHigh(a.fhDeltaRMultTriggersGenLevelHigh),
fhDeltaPhiMultTriggersGenLevelHigh(a.fhDeltaPhiMultTriggersGenLevelHigh),
fIsChargedMC(a.fIsChargedMC),
fIsKine(a.fIsKine),
fIsFullMC(a.fIsFullMC),
faGenIndex(a.faGenIndex),
faRecIndex(a.faRecIndex),
fkAcceptance(a.fkAcceptance),
fkDeltaPhiCut(a.fkDeltaPhiCut),
fh1Xsec(a.fh1Xsec),
fh1Trials(a.fh1Trials),
fh1AvgTrials(a.fh1AvgTrials),
fh1PtHard(a.fh1PtHard),
fh1PtHardNoW(a.fh1PtHardNoW),  
fh1PtHardTrials(a.fh1PtHardTrials),
fAvgTrials(a.fAvgTrials),
fHardest(a.fHardest),
fEventNumberRangeLow(a.fEventNumberRangeLow),
fEventNumberRangeHigh(a.fEventNumberRangeHigh),
fTriggerPtRangeLow(a.fTriggerPtRangeLow),
fTriggerPtRangeHigh(a.fTriggerPtRangeHigh),
fFillRespMx(a.fFillRespMx),
fRandom(a.fRandom),
fnTrials(a.fnTrials),
fJetFreeAreaFrac(a.fJetFreeAreaFrac),
fnPhi(a.fnPhi),
fnEta(a.fnEta),
fEtaSize(a.fEtaSize),
fPhiSize(a.fPhiSize),
fCellArea(a.fCellArea),
fSafetyMargin(a.fSafetyMargin),
fDoubleBinning(a.fDoubleBinning)
{
   //Copy Constructor
   fRandom->SetSeed(0);
}
//--------------------------------------------------------------

AliAnalysisTaskJetCorePP& AliAnalysisTaskJetCorePP::operator = (const AliAnalysisTaskJetCorePP& a){
  // assignment operator
  this->~AliAnalysisTaskJetCorePP();
  new(this) AliAnalysisTaskJetCorePP(a);
  return *this;
}
//--------------------------------------------------------------

AliAnalysisTaskJetCorePP::~AliAnalysisTaskJetCorePP()
{
   //Destructor 
   delete fListJets;
   delete fListJetsGen;
   delete fListJetsGenFull;
   delete fListJetsBg;
   delete fListJetsBgGen;
   delete fOutputList; // ????
   delete fRandom;
}

//--------------------------------------------------------------


Bool_t AliAnalysisTaskJetCorePP::Notify()
{
   //Implemented Notify() to read the cross sections
   //and number of trials from pyxsec.root
   //inspired by AliAnalysisTaskJetSpectrum2::Notify()
   if(!(fIsChargedMC || fIsKine)) return kFALSE; 
   Float_t xsection = 0;
   Float_t trials  = 1;
   fAvgTrials = 1;

   if(fIsChargedMC){ 
      fESD = dynamic_cast<AliESDEvent*>(InputEvent());
      if(!fESD){
         if(fDebug>1) AliError("ESD not available");
         fAODIn = dynamic_cast<AliAODEvent*>(InputEvent());
      } 
 
      fAODOut = dynamic_cast<AliAODEvent*>(AODEvent());


      if(fNonStdFile.Length()!=0){
         // case that we have an AOD extension we can fetch the jets from the extended output
         AliAODHandler *aodH = dynamic_cast<AliAODHandler*>(AliAnalysisManager::GetAnalysisManager()->GetOutputEventHandler());
         fAODExtension = aodH ? aodH->GetExtension(fNonStdFile.Data()) : 0;
         if(!fAODExtension){
            if(fDebug>1) Printf("AODExtension found for %s",fNonStdFile.Data());
         } 
      }
 
      TTree *tree = AliAnalysisManager::GetAnalysisManager()->GetTree();

      if(tree){
         TFile *curfile = tree->GetCurrentFile();
         if(!curfile) {
            Error("Notify","No current file");
            return kFALSE;
         }
         if(!fh1Xsec || !fh1Trials){
            Printf("%s%d No Histogram fh1Xsec",(char*)__FILE__,__LINE__);
            return kFALSE;
         }
         AliAnalysisHelperJetTasks::PythiaInfoFromFile(curfile->GetName(),xsection,trials);
         fh1Xsec->Fill("<#sigma>",xsection);
         // construct a poor man average trials
         Float_t nEntries = (Float_t)tree->GetTree()->GetEntries();
         if(trials>=nEntries && nEntries>0.) fAvgTrials = trials/nEntries;
         fh1Trials->Fill("#sum{ntrials}",trials);
      }  

      if(fDebug)Printf("Reading File %s",fInputHandler->GetTree()->GetCurrentFile()->GetName());
   }


   return kTRUE;
}
//--------------------------------------------------------------

void AliAnalysisTaskJetCorePP::Init()
{
   // check for jet branches
   if(fJetBranchNameKine.Length()==0){
      if(!strlen(fJetBranchName.Data())){
         AliError("Jet branch name not set.");
      }
   }
   fMcHandler = dynamic_cast<AliInputEventHandler*> (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler()); 


}

//--------------------------------------------------------------

void AliAnalysisTaskJetCorePP::UserCreateOutputObjects()
{
   // Create histograms and initilize variables
   fSafetyMargin = fBgConeR*fBgConeR /(fBgJetParamR*fBgJetParamR);

   // Called once
   fListJets   = new TList();  //reconstructed level antikt jets
   fListJetsBg = new TList();  //reconstructed jets to be removed from bg

   fIsChargedMC = (fJetBranchNameChargMC.Length()>0) ? kTRUE : kFALSE;
   fIsKine      = (fJetBranchNameKine.Length()>0)    ? kTRUE : kFALSE;
   fIsFullMC    = (fJetBranchNameFullMC.Length()>0)  ? kTRUE : kFALSE;

   fRandom = new TRandom3(0);

   if(fIsChargedMC || fIsKine){   //full MC or pure kine
      fListJetsGen   = new TList(); //generator level charged antikt jets
      fListJetsBgGen = new TList(); //generator level jets to be removed from bg 

      if(fIsFullMC)
         fListJetsGenFull = new TList(); //generator level full jets
   }
   OpenFile(1);
   if(!fOutputList) fOutputList = new TList();
   fOutputList->SetOwner(kTRUE);

   Bool_t oldStatus = TH1::AddDirectoryStatus();
   TH1::AddDirectory(kFALSE);

   fHistEvtSelection = new TH1I("fHistEvtSelection", "event selection", 6, -0.5, 5.5);
   fHistEvtSelection->GetXaxis()->SetBinLabel(1,"ACCEPTED");
   fHistEvtSelection->GetXaxis()->SetBinLabel(2,"events IN");
   fHistEvtSelection->GetXaxis()->SetBinLabel(3,"event selection (rejected)");
   fHistEvtSelection->GetXaxis()->SetBinLabel(4,"vertex cut (rejected)");
   fHistEvtSelection->GetXaxis()->SetBinLabel(5,"centrality (rejected)");
   fHistEvtSelection->GetXaxis()->SetBinLabel(6,"event number (rejected)");
   
   fOutputList->Add(fHistEvtSelection);

   Int_t nBinsCentrality = (fSystem==0) ? 1 : 10; // pp=1 else 10
   //trigger pt spectrum (reconstructed) 
   fh2Ntriggers = new TH2F("fh2Ntriggers","# of triggers",
                             nBinsCentrality,0.0,100.0,50,0.0,50.0);
   if(!fIsKine) fOutputList->Add(fh2Ntriggers);

   Int_t bw = (fDoubleBinning==0) ? 1 : 2; //make larger bin width

   //Centrality, A, pTjet, pTtrigg, dphi
   const Int_t dimSpec   = 5;
   const Int_t nBinsSpec[dimSpec]     = {nBinsCentrality,  50, bw*110,  50, TMath::Nint(10*(TMath::Pi()-fkDeltaPhiCut))};
   const Double_t lowBinSpec[dimSpec] = {0.0,             0.0,    -20, 0.0, fkDeltaPhiCut};
   const Double_t hiBinSpec[dimSpec]  = {100.0,           1.0,  200.0,50.0, TMath::Pi()};
   fHJetSpec = new THnSparseF("fHJetSpec",
                   "Recoil jet spectrum [cent,A,pTjet,pTtrig,dphi]",
                   dimSpec,nBinsSpec,lowBinSpec,hiBinSpec);
   if(!fIsKine) fOutputList->Add(fHJetSpec);  

   //background estimated as  median of kT jets 
   fHJetSpecSubUeMedian = (THnSparseF*) fHJetSpec->Clone("fHJetSpecSubUeMedian");
   fHJetSpecSubUeMedian->SetTitle("Recoil jet spectrum [cent,A,pTjet-pTUe,pTtrig,dphi]");
   if(!fIsKine) fOutputList->Add(fHJetSpecSubUeMedian); 
   //background estimated as weighted  median of kT jets  ala Cone
   fHJetSpecSubUeCone = (THnSparseF*) fHJetSpec->Clone("fHJetSpecSubUeCone");
   fHJetSpecSubUeCone->SetTitle("Recoil jet spectrum [cent,A,pTjet-pTUe,pTtrig,dphi]");
   if(!fIsKine) fOutputList->Add(fHJetSpecSubUeCone); 


   //A, pTjet, pTtrigg, dphi
   const Int_t dimCor = 4;
   const Int_t nBinsCor[dimCor]     = {50,  110,   50,      100};
   const Double_t lowBinCor[dimCor] = {0.0, -20,  0.0, -0.5*TMath::Pi()};
   const Double_t hiBinCor[dimCor]  = {1.0, 200, 50.0, 1.5*TMath::Pi()};
   fHJetPhiCorr = new THnSparseF("fHJetPhiCorr",
                                 "Recoil jet spectrum [A,pTjet,pTtrig,dphi]",
                                 dimCor,nBinsCor,lowBinCor,hiBinCor);

   if(!fIsKine) fOutputList->Add(fHJetPhiCorr); // Dphi distribution jet-triger

   //background estimated as  median of kT jets 
   fHJetPhiCorrSubUeMedian = (THnSparseF*) fHJetPhiCorr->Clone("fHJetPhiCorrSubUeMedian");
   fHJetPhiCorrSubUeMedian->SetTitle("Recoil jet spectrum [A,pTjet-pTUe,pTtrig,dphi]");
   if(!fIsKine) fOutputList->Add(fHJetPhiCorrSubUeMedian); 
   //background estimated as weighted  median of kT jets  ala Cone
   fHJetPhiCorrSubUeCone = (THnSparseF*) fHJetPhiCorr->Clone("fHJetPhiCorrSubUeCone");
   fHJetPhiCorrSubUeCone->SetTitle("Recoil jet spectrum [A,pTjet-pTUe,pTtrig,dphi]");
   if(!fIsKine) fOutputList->Add(fHJetPhiCorrSubUeCone); 


   //------------------- HISTOS FOR DIAGNOSTIC ----------------------
   //A, pTjet, pTjet-pTUe, pTUe, rhoUe     bg estimate from kT median
   const Int_t    dimSpecMed   = 5;
   const Int_t    nBinsSpecMed[dimSpecMed]  = {25,     50,    60,    20,   20};
   const Double_t lowBinSpecMed[dimSpecMed] = {0.0,   0.0, -20.0,   0.0,  0.0};
   const Double_t hiBinSpecMed[dimSpecMed]  = {1.0, 100.0, 100.0,  10.0, 20.0};
   fHJetUeMedian = new THnSparseF("fHJetUeMedian",
                   "Recoil jet spectrum [A,pTjet,pTjet-pTUe, pTUe, rhoUe]",
                   dimSpecMed, nBinsSpecMed, lowBinSpecMed, hiBinSpecMed);
   if(!fIsKine) fOutputList->Add(fHJetUeMedian);  
   
   //A, pTjet, pTjet-pTUe, pTUe, rhoUe     bg estimate from kT median Cone
   fHJetUeCone = (THnSparseF*) fHJetUeMedian->Clone("fHJetUeCone");
   fHJetUeCone->SetTitle("Recoil jet spectrum [A,pTjet,pTjet-pTUe, pTUe, rhoUe]");
   if(!fIsKine) fOutputList->Add(fHJetUeCone); 

   //rho bacground reconstructed data
   const Int_t    dimRho   = 2;
   const Int_t    nBinsRho[dimRho]  = {50  ,   50};
   const Double_t lowBinRho[dimRho] = {0.0  , 0.0};
   const Double_t hiBinRho[dimRho]  = {20.0 , 20.0};

   fHRhoUeMedianVsCone = new THnSparseF("hRhoUeMedianVsCone","[Rho Cone, Rho Median]",  
                                      dimRho, nBinsRho, lowBinRho, hiBinRho);
   if(!fIsKine) fOutputList->Add(fHRhoUeMedianVsCone);

   //Jet number density histos [Trk Mult, jet density, pT trigger]
   /*const Int_t    dimJetDens   = 3;
   const Int_t    nBinsJetDens[dimJetDens]  = {100,   100,   10};
   const Double_t lowBinJetDens[dimJetDens] = {0.0,   0.0,  0.0};
   const Double_t hiBinJetDens[dimJetDens]  = {500.0, 5.0, 50.0};

   fHJetDensity = new THnSparseF("fHJetDensity","Jet dens vs trk mult A>0.07",
                                   dimJetDens,nBinsJetDens,lowBinJetDens,hiBinJetDens);

   fHJetDensityA4 =new THnSparseF("fHJetDensityA4","Jet dens vs trk mult A>0.4",
                                   dimJetDens,nBinsJetDens,lowBinJetDens,hiBinJetDens);

   fOutputList->Add(fHJetDensity);
   fOutputList->Add(fHJetDensityA4);
   */      

   //inclusive azimuthal and pseudorapidity histograms
   fhJetPhi = new TH2D("fhJetPhi","Azim dist jets vs pTjet",
                        50, 0, 100, 50,-TMath::Pi(),TMath::Pi());
   fhTriggerPhi= new TH2D("fhTriggerPhi","azim dist trig had vs pTtrigg",
                        50, 0, 50, 50,-TMath::Pi(),TMath::Pi());
   fhJetEta = new TH2D("fhJetEta","Eta dist jets vs pTjet",
                        50,0, 100, 40,-0.9,0.9);
   fhTriggerEta = new TH2D("fhTriggerEta","Eta dist trig had vs pTtrigg",
                        50, 0, 50, 40,-0.9,0.9);

   fhVertexZ = new TH1D("fhVertexZ","z vertex",40,-20,20);  
   fhVertexZAccept = new TH1D("fhVertexZAccept","z vertex after cut",40,-20,20);  
   fhContribVtx = new TH1D("fhContribVtx","contrib to vtx",200,0,200);   
   fhContribVtxAccept = new TH1D("fhContribVtxAccept","contrib to vtx after cut",200,0,200);   
   fhDphiTriggerJet = new TH1D("fhDphiTriggerJet","Deltaphi trig-jet",50, -TMath::Pi(),TMath::Pi()); 
   fhDphiTriggerJetAccept = new TH1D("fhDphiTriggerJetAccept","Deltaphi trig-jet after cut",50, -TMath::Pi(),TMath::Pi()); 
   fhCentrality = new TH1D("fhCentrality","Centrality",20,0,100);
   fhCentralityAccept = new TH1D("fhCentralityAccept","Centrality after cut",20,0,100);
   fhEntriesToMedian = new TH1D("fhEntriesToMedian","fhEntriesToMedian",30,0,30);
   fhCellAreaToMedian =  new TH1D("fhCellAreaToMedian", "fhCellAreaToMedian", 75,0,1.5);

   fhNofMultipleTriggers = new TH1D("fhNofMultipleTriggers","fhNofMultipleTriggers",20,0,20);
   fhNofMultipleTriggersCone = new TH1D("fhNofMultipleTriggersCone","fhNofMultipleTriggersCone R<0.4",20,0,20);
   fhNofMultipleTriggersConeLow  = (TH1D*) fhNofMultipleTriggersCone->Clone("fhNofMultipleTriggersConepTa15to20");
   fhNofMultipleTriggersConeHigh = (TH1D*) fhNofMultipleTriggersCone->Clone("fhNofMultipleTriggersConepTa20to50");

   fhDeltaRMultTriggersLow = new  TH1D("fhDeltaRMultTriggersLow","fhDeltaRMultTriggersLow", 200,0,4);  
   fhDeltaRMultTriggersHigh = (TH1D*) fhDeltaRMultTriggersLow->Clone("fhDeltaRMultTriggersHigh"); 

   fhDeltaPhiMultTriggersLow  = new  TH1D("fhDeltaPhiMultTriggersLow","fhDeltaPhiRultTriggers", 320,-3.2,3.2); //single incl trigger 
   fhDeltaPhiMultTriggersHigh = (TH1D*) fhDeltaPhiMultTriggersLow->Clone("fhDeltaPhiMultTriggersHigh");  //single incl trigger 

   fhDeltaPhiMultTriggersInclLow = (TH1D*) fhDeltaPhiMultTriggersLow->Clone("fhDeltaPhiMultTriggersInclLow");
   fhDeltaPhiMultTriggersInclHigh = (TH1D*) fhDeltaPhiMultTriggersLow->Clone("fhDeltaPhiMultTriggersInclHigh");
   fhInclTrigCounter = new TH1D("fhInclTrigCounter","fhInclTrigCounter",1,0,1);

   fhDeltaPtConeBg   = new TH1D("fhDeltaPtConeBg","fhDeltaPtConeBg",2000,-20,80);
   fhDeltaPtMedianBg = new TH1D("fhDeltaPtMedianBg","fhDeltaPtMedianBg",2000,-20,80);

   if(!fIsKine){
      fOutputList->Add(fhJetPhi);
      fOutputList->Add(fhTriggerPhi);
      fOutputList->Add(fhJetEta);
      fOutputList->Add(fhTriggerEta);
   }
   fOutputList->Add(fhVertexZ);    
   fOutputList->Add(fhVertexZAccept);   
   if(!fIsKine){
      fOutputList->Add(fhContribVtx); 
      fOutputList->Add(fhContribVtxAccept); 
      fOutputList->Add(fhDphiTriggerJet);
      fOutputList->Add(fhDphiTriggerJetAccept);
      fOutputList->Add(fhCentrality); 
      fOutputList->Add(fhCentralityAccept);
      fOutputList->Add(fhEntriesToMedian);
      fOutputList->Add(fhCellAreaToMedian);
      fOutputList->Add(fhNofMultipleTriggers);
      fOutputList->Add(fhNofMultipleTriggersCone);
      fOutputList->Add(fhNofMultipleTriggersConeLow);
      fOutputList->Add(fhNofMultipleTriggersConeHigh);
      fOutputList->Add(fhDeltaRMultTriggersLow);
      fOutputList->Add(fhDeltaRMultTriggersHigh);
      fOutputList->Add(fhDeltaPhiMultTriggersLow);
      fOutputList->Add(fhDeltaPhiMultTriggersHigh);
      fOutputList->Add(fhDeltaPhiMultTriggersInclLow);
      fOutputList->Add(fhDeltaPhiMultTriggersInclHigh);
      fOutputList->Add(fhInclTrigCounter);
   }
   fOutputList->Add(fhDeltaPtConeBg);   
   fOutputList->Add(fhDeltaPtMedianBg); 


   // raw spectra of INCLUSIVE jets  
   //Centrality, pTjet, A
   /*const Int_t dimRaw   = 3;
   const Int_t nBinsRaw[dimRaw]     = {nBinsCentrality,  50,   100};
   const Double_t lowBinRaw[dimRaw] = {0.0,             0.0,   0.0};
   const Double_t hiBinRaw[dimRaw]  = {100.0,           100,   1.0};
   fHJetPtRaw = new THnSparseF("fHJetPtRaw",
                                "Incl. jet spectrum [cent,pTjet,A]",
                                dimRaw,nBinsRaw,lowBinRaw,hiBinRaw);
   fOutputList->Add(fHJetPtRaw);  

   // raw spectra of LEADING jets  
   //Centrality, pTjet, A
   fHLeadingJetPtRaw = new THnSparseF("fHLeadingJetPtRaw",
                                "Leading jet spectrum [cent,pTjet,A]",
                                dimRaw,nBinsRaw,lowBinRaw,hiBinRaw);
   fOutputList->Add(fHLeadingJetPtRaw);  

   // Dphi versus pT jet 
   //Centrality, Dphi=phiTrig-phiJet, pTjet, pTtrigg 
   const Int_t dimDp   = 4;
   const Int_t nBinsDp[dimDp]     = {nBinsCentrality,  50,     50,    50};
   const Double_t lowBinDp[dimDp] = {0.0,       -TMath::Pi(),   0.0,   0.0};
   const Double_t hiBinDp[dimDp]  = {100.0,      TMath::Pi(), 100.0, 100.0};
   fHDphiVsJetPtAll = new THnSparseF("fHDphiVsJetPtAll",
                                "Dphi vs jet pT [cent,Dphi,pTjet,pTtrigg]",
                                dimDp,nBinsDp,lowBinDp,hiBinDp);
   fOutputList->Add(fHDphiVsJetPtAll);  
   */

   //analyze MC generator level 
   if(fIsChargedMC || fIsKine){   
      if(fFillRespMx){
         //Fill response matrix only once 
         fhJetPtGenVsJetPtRec = new TH2D("fhJetPtGenVsJetPtRec","JetPtGenVsJetPtRec", bw*100,0,200, bw*100,0,200); 
         fOutputList->Add(fhJetPtGenVsJetPtRec); //gen MC charg jet pt spectrum versus rec charged jet pt spectrum
         //....
         fhJetPtGenVsJetPtRecSubUeMedian = new TH2D("fhJetPtGenVsJetPtRecSubUeMedian","fhJetPtGenVsJetPtRecSubUeMedian", bw*110,-20,200, bw*110,-20,200); 
         fOutputList->Add(fhJetPtGenVsJetPtRecSubUeMedian); // with kT median bg subtr 
         //....
         fhJetPtGenVsJetPtRecSubUeCone=(TH2D*)fhJetPtGenVsJetPtRecSubUeMedian->Clone("fhJetPtGenVsJetPtRecSubUeCone");
         fhJetPtGenVsJetPtRecSubUeCone->SetTitle("fhJetPtGenVsJetPtRecSubUeCone");
         fOutputList->Add(fhJetPtGenVsJetPtRecSubUeCone); // with weighted kT median bg subtr 
         //....
         fhJetPtGen = new TH1D("fhJetPtGen","Jet Pt (MC Gen)",bw*100,0,200); //MC generator charged jet pt spectrum
         fOutputList->Add(fhJetPtGen);  
         //....
         fhJetPtSubUeMedianGen = new TH1D("fhJetPtSubUeMedianGen","Jet Pt - UE pT (MC Gen)",bw*110,-20,200); 
         fOutputList->Add(fhJetPtSubUeMedianGen);  // with kT median bg subtr
         //....
         fhJetPtSubUeConeGen = (TH1D*) fhJetPtSubUeMedianGen->Clone("fhJetPtSubUeConeGen");
         fOutputList->Add(fhJetPtSubUeConeGen); // with weighted kT median bg subtr
         //....     
         fhJetPtResolutionVsPtGen = new TH2D("fhJetPtResolutionVsPtGen","fhJetPtResolutionVsPtGen",20,0,100, 35,-1.,0.4);
         fOutputList->Add(fhJetPtResolutionVsPtGen); // with weighted kT median bg subtr
          //....     
         fhJetPtResolutionVsPtConeGen = (TH2D*) fhJetPtResolutionVsPtGen->Clone("fhJetPtResolutionVsPtConeGen");
         fOutputList->Add(fhJetPtResolutionVsPtConeGen); // with weighted kT median bg subtr
 
         //....
         if(fIsFullMC){
            fhJetPtGenChargVsJetPtGenFull = new TH2D("fhJetPtGenChargVsJetPtGenFull","fhJetPtGenChargVsJetPtGenFull", 100,0,200, 100,0,200);
            fOutputList->Add(fhJetPtGenChargVsJetPtGenFull); //gen full MC jet pt versus gen charged jet MC pt
         //....
            fhJetPtGenFull = new TH1D("fhJetPtGenFull","Jet Pt (MC Full jets Gen)",100,0,200); //MC generator full jet pt spectrum
            fOutputList->Add(fhJetPtGenFull); 
         }
      }
      //....
      fh2NtriggersGen = (TH2F*) fh2Ntriggers->Clone("fh2NtriggersGen");
      fh2NtriggersGen->SetTitle(Form("%s Gen MC",fh2Ntriggers->GetTitle()));
      fOutputList->Add(fh2NtriggersGen);

      //Centrality, A, pT, pTtrigg
      fHJetSpecGen = (THnSparseF*) fHJetSpec->Clone("fHJetSpecGen");
      fHJetSpecGen->SetTitle(Form("%s Gen MC",fHJetSpec->GetTitle()));
      fOutputList->Add(fHJetSpecGen); 

      fHJetSpecSubUeMedianGen = (THnSparseF*)  fHJetSpecSubUeMedian->Clone("fHJetSpecSubUeMedianGen");
      fHJetSpecSubUeMedianGen->SetTitle(Form("%s Gen MC",fHJetSpecSubUeMedian->GetTitle()));
      fOutputList->Add(fHJetSpecSubUeMedianGen);  

      fHJetSpecSubUeConeGen =  (THnSparseF*) fHJetSpecSubUeCone->Clone("fHJetSpecSubUeConeGen"); 
      fHJetSpecSubUeConeGen->SetTitle(Form("%s Gen MC",fHJetSpecSubUeCone->GetTitle()));
      fOutputList->Add(fHJetSpecSubUeConeGen);
      //---
      fHJetPhiCorrGen = (THnSparseF*) fHJetPhiCorr->Clone("fHJetPhiCorrGen");
      fHJetPhiCorrGen->SetTitle(Form("%s Gen MC",fHJetPhiCorr->GetTitle()));
      fOutputList->Add(fHJetPhiCorrGen); 

      fHJetPhiCorrSubUeMedianGen = (THnSparseF*) fHJetPhiCorrSubUeMedian->Clone("fHJetPhiCorrSubUeMedianGen");
      fHJetPhiCorrSubUeMedianGen->SetTitle(Form("%s Gen MC",fHJetPhiCorrSubUeMedian->GetTitle()));
      fOutputList->Add(fHJetPhiCorrSubUeMedianGen); 
 
      fHJetPhiCorrSubUeConeGen = (THnSparseF*) fHJetPhiCorrSubUeCone->Clone("fHJetPhiCorrSubUeConeGen");
      fHJetPhiCorrSubUeConeGen->SetTitle(Form("%s Gen MC", fHJetPhiCorrSubUeCone->GetTitle()));
      fOutputList->Add(fHJetPhiCorrSubUeConeGen);

      //---
      fHJetUeMedianGen  =  (THnSparseF*) fHJetUeMedian->Clone("fHJetUeMedianGen");  
      fHJetUeMedianGen->SetTitle(Form("%s Gen MC", fHJetUeMedian->GetTitle()));
      fOutputList->Add(fHJetUeMedianGen);

      fHJetUeConeGen =  (THnSparseF*) fHJetUeCone->Clone("fHJetUeConeGen"); 
      fHJetUeConeGen->SetTitle(Form("%s Gen MC", fHJetUeCone->GetTitle())); 
      fOutputList->Add(fHJetUeConeGen);

      if(fFillRespMx){
         //track efficiency/contamination histograms  eta versus pT
         Double_t bins [] = {0, 0.2,0.4,0.6, 0.8, 1., 1.2, 1.4, 1.6, 1.8, 2., 2.5, 3., 3.5, 4., 5., 6., 8., 10., 20., 50.};
         Int_t nbins = sizeof(bins)/sizeof(Double_t)-1;
         fhPtTrkTruePrimRec = new TH2D("fhPtTrkTruePrimRec","PtTrkTruePrimRec",nbins, bins, 18,-0.9,0.9);
         fOutputList->Add(fhPtTrkTruePrimRec); 
      
         fhPtTrkTruePrimGen = (TH2D*) fhPtTrkTruePrimRec->Clone("fhPtTrkTruePrimGen");
         fhPtTrkTruePrimGen->SetTitle("PtTrkTruePrimGen");    
         fOutputList->Add(fhPtTrkTruePrimGen);
      
         fhPtTrkSecOrFakeRec = (TH2D*) fhPtTrkTruePrimRec->Clone("fhPtTrkSecOrFakeRec");    
         fhPtTrkSecOrFakeRec->SetTitle("PtTrkSecOrFakeRec");    
         fOutputList->Add(fhPtTrkSecOrFakeRec);
      }

      fHRhoUeMedianVsConeGen = (THnSparseF*) fHRhoUeMedianVsCone->Clone("hRhoUeMedianVsConeGen");
      fHRhoUeMedianVsConeGen->SetTitle(Form("%s Gen MC", fHRhoUeMedianVsCone->GetTitle())); 
      fOutputList->Add(fHRhoUeMedianVsConeGen);

      fhEntriesToMedianGen = (TH1D*) fhEntriesToMedian->Clone("fhEntriesToMedianGen");
      fhEntriesToMedianGen->SetTitle(Form("%s Gen MC", fhEntriesToMedian->GetTitle())); 
      fOutputList->Add(fhEntriesToMedianGen);

      fhCellAreaToMedianGen  = (TH1D*) fhCellAreaToMedian->Clone("fhCellAreaToMedianGen");
      fhCellAreaToMedianGen->SetTitle(Form("%s Gen MC", fhCellAreaToMedian->GetTitle())); 
      fOutputList->Add(fhCellAreaToMedianGen);

      fhNofMultipleTriggersGen = (TH1D*) fhNofMultipleTriggers->Clone("fhNofMultipleTriggersGen"); 
      fOutputList->Add(fhNofMultipleTriggersGen);

      fhNofMultipleTriggersConeGen = (TH1D*) fhNofMultipleTriggersCone->Clone("fhNofMultipleTriggersConeGen"); 
      fOutputList->Add(fhNofMultipleTriggersConeGen);

      fhNofMultipleTriggersConeGenLow = (TH1D*) fhNofMultipleTriggersCone->Clone("fhNofMultipleTriggersConeGenpta15to20"); 
      fOutputList->Add(fhNofMultipleTriggersConeGenLow);

      fhNofMultipleTriggersConeGenHigh = (TH1D*) fhNofMultipleTriggersCone->Clone("fhNofMultipleTriggersConeGenpta20to50"); 
      fOutputList->Add(fhNofMultipleTriggersConeGenHigh);

      fhDeltaRMultTriggersGenLow  = (TH1D*) fhDeltaRMultTriggersLow->Clone("fhDeltaRMultTriggersGenLow");
      fOutputList->Add(fhDeltaRMultTriggersGenLow);

      fhDeltaRMultTriggersGenHigh  = (TH1D*) fhDeltaRMultTriggersLow->Clone("fhDeltaRMultTriggersGenHigh");
      fOutputList->Add(fhDeltaRMultTriggersGenHigh);

      fhDeltaPhiMultTriggersGenLow  = (TH1D*) fhDeltaPhiMultTriggersLow->Clone("fhDeltaPhiMultTriggersGenLow");
      fOutputList->Add(fhDeltaPhiMultTriggersGenLow);

      fhDeltaPhiMultTriggersGenHigh  = (TH1D*) fhDeltaPhiMultTriggersLow->Clone("fhDeltaPhiMultTriggersGenHigh");
      fOutputList->Add(fhDeltaPhiMultTriggersGenHigh);

      fhDeltaPhiMultTriggersInclGenLow  = (TH1D*) fhDeltaPhiMultTriggersInclLow->Clone("fhDeltaPhiMultTriggersInclGenLow");
      fOutputList->Add(fhDeltaPhiMultTriggersInclGenLow);

      fhDeltaPhiMultTriggersInclGenHigh = (TH1D*) fhDeltaPhiMultTriggersInclHigh->Clone("fhDeltaPhiMultTriggersInclGenHigh");
      fOutputList->Add(fhDeltaPhiMultTriggersInclGenHigh);

      fhInclTrigCounterGen= (TH1D*) fhInclTrigCounter->Clone("fhInclTrigCounterGen");
      fOutputList->Add(fhInclTrigCounterGen);

      fhNofMultipleTriggersConeGenA = (TH1D*) fhNofMultipleTriggersConeGen->Clone("fhNofMultipleTriggersConeGen10"); 
      fOutputList->Add(fhNofMultipleTriggersConeGenA);

      fhNofMultipleTriggersConeGenALow = (TH1D*) fhNofMultipleTriggersConeGen->Clone("fhNofMultipleTriggersConeGen10pTa15to20"); 
      fOutputList->Add(fhNofMultipleTriggersConeGenALow);

      fhNofMultipleTriggersConeGenAHigh = (TH1D*) fhNofMultipleTriggersConeGen->Clone("fhNofMultipleTriggersConeGen10pTa20to50"); 
      fOutputList->Add(fhNofMultipleTriggersConeGenAHigh);

      fhDeltaRMultTriggersGenALow = (TH1D*) fhDeltaRMultTriggersLow->Clone("fhDeltaRMultTriggersGen10pTa15to20");
      fOutputList->Add(fhDeltaRMultTriggersGenALow);

      fhDeltaPhiMultTriggersGenALow = (TH1D*) fhDeltaPhiMultTriggersLow->Clone("fhDeltaPhiMultTriggersGen10pTa15to20");
      fOutputList->Add(fhDeltaPhiMultTriggersGenALow);

      fhDeltaRMultTriggersGenAHigh = (TH1D*) fhDeltaRMultTriggersLow->Clone("fhDeltaRMultTriggersGen10pTa20to50");
      fOutputList->Add(fhDeltaRMultTriggersGenAHigh);

      fhDeltaPhiMultTriggersGenAHigh =  (TH1D*) fhDeltaPhiMultTriggersLow->Clone("fhDeltaPhiMultTriggersGen10pTa20to50");
      fOutputList->Add(fhDeltaPhiMultTriggersGenAHigh);

      fhNofTriggersGenA = (TH1D*) fhInclTrigCounter->Clone("fhNofTriggersGen10");
      fOutputList->Add(fhNofTriggersGenA);


      fhNofMultipleTriggersConeGenB = (TH1D*) fhNofMultipleTriggersConeGen->Clone("fhNofMultipleTriggersConeGen5"); 
      fOutputList->Add(fhNofMultipleTriggersConeGenB);

      fhNofMultipleTriggersConeGenBLow = (TH1D*) fhNofMultipleTriggersConeGen->Clone("fhNofMultipleTriggersConeGen5pTa15to20"); 
      fOutputList->Add(fhNofMultipleTriggersConeGenBLow);

      fhNofMultipleTriggersConeGenBHigh = (TH1D*) fhNofMultipleTriggersConeGen->Clone("fhNofMultipleTriggersConeGen5pTa20to50"); 
      fOutputList->Add(fhNofMultipleTriggersConeGenBHigh);

      fhDeltaRMultTriggersGenBLow = (TH1D*) fhDeltaRMultTriggersLow->Clone("fhDeltaRMultTriggersGen5pTa15to20");
      fOutputList->Add(fhDeltaRMultTriggersGenBLow);

      fhDeltaPhiMultTriggersGenBLow = (TH1D*) fhDeltaPhiMultTriggersLow->Clone("fhDeltaPhiMultTriggersGen5pTa15to20");
      fOutputList->Add(fhDeltaPhiMultTriggersGenBLow);

      fhDeltaRMultTriggersGenBHigh = (TH1D*) fhDeltaRMultTriggersLow->Clone("fhDeltaRMultTriggersGen5pTa20to50");
      fOutputList->Add(fhDeltaRMultTriggersGenBHigh);

      fhDeltaPhiMultTriggersGenBHigh =  (TH1D*) fhDeltaPhiMultTriggersLow->Clone("fhDeltaPhiMultTriggersGen5pTa20to50");
      fOutputList->Add(fhDeltaPhiMultTriggersGenBHigh);

      fhNofTriggersGenB = (TH1D*) fhInclTrigCounter->Clone("fhNofTriggersGen5");
      fOutputList->Add(fhNofTriggersGenB);




      if(fIsChargedMC){  // Filled from gen level if  well reconstructed trigger exists
         fhTriggerCounterGenLevel = (TH1D*) fhInclTrigCounter->Clone("fhTriggerCounterGenLevel");
         fOutputList->Add(fhTriggerCounterGenLevel);

         fhDeltaRMultTriggersGenLevelLow = (TH1D*) fhDeltaRMultTriggersLow->Clone("fhDeltaRMultTriggersGenLevelLow");
         fOutputList->Add(fhDeltaRMultTriggersGenLevelLow);

         fhDeltaPhiMultTriggersGenLevelLow = (TH1D*) fhDeltaPhiMultTriggersLow->Clone("fhDeltaPhiMultTriggersGenLevelLow");
         fOutputList->Add(fhDeltaPhiMultTriggersGenLevelLow);

         fhDeltaRMultTriggersGenLevelHigh = (TH1D*) fhDeltaRMultTriggersLow->Clone("fhDeltaRMultTriggersGenLevelHigh");
         fOutputList->Add(fhDeltaRMultTriggersGenLevelHigh);
   
         fhDeltaPhiMultTriggersGenLevelHigh = (TH1D*) fhDeltaPhiMultTriggersLow->Clone("fhDeltaPhiMultTriggersGenLevelHigh");
         fOutputList->Add(fhDeltaPhiMultTriggersGenLevelHigh);
      }
   }
   //-------------------------------------
   //     pythia histograms
   const Int_t nBinPt = 150;
   Double_t binLimitsPt[nBinPt+1];
   for(Int_t iPt = 0;iPt <= nBinPt;iPt++){
      if(iPt == 0){
         binLimitsPt[iPt] = -50.0;
      }else{// 1.0
         binLimitsPt[iPt] =  binLimitsPt[iPt-1] + 2.;
      }
   }
   
   fh1Xsec = new TProfile("fh1Xsec","xsec from pyxsec.root",1,0,1);
   fh1Xsec->GetXaxis()->SetBinLabel(1,"<#sigma>");
   fOutputList->Add(fh1Xsec);
   fh1Trials = new TH1F("fh1Trials","trials root file",1,0,1);
   fh1Trials->GetXaxis()->SetBinLabel(1,"#sum{ntrials}");
   fOutputList->Add(fh1Trials);
   fh1AvgTrials = new TH1F("fh1AvgTrials","trials root file",1,0,1);
   fh1AvgTrials->GetXaxis()->SetBinLabel(1,"#sum{avg ntrials}");
   fOutputList->Add(fh1AvgTrials);
   fh1PtHard = new TH1F("fh1PtHard","PYTHIA Pt hard;p_{T,hard}",nBinPt,binLimitsPt);
   fOutputList->Add(fh1PtHard);
   fh1PtHardNoW = new TH1F("fh1PtHardNoW","PYTHIA Pt hard no weight;p_{T,hard}",nBinPt,binLimitsPt);
   fOutputList->Add(fh1PtHardNoW);
   fh1PtHardTrials = new TH1F("fh1PtHardTrials","PYTHIA Pt hard weight with trials;p_{T,hard}",nBinPt,binLimitsPt);
   fOutputList->Add(fh1PtHardTrials);      
   

   // =========== Switch on Sumw2 for all histos ===========
   for(Int_t i=0; i<fOutputList->GetEntries(); i++){
      TH1 *h1 = dynamic_cast<TH1*>(fOutputList->At(i));
      if(h1){
         h1->Sumw2();
         continue;
      }
      THnSparse *hn = dynamic_cast<THnSparse*>(fOutputList->At(i));
      if(hn){
         hn->Sumw2();
      }	  
   }
   TH1::AddDirectory(oldStatus);

   PostData(1, fOutputList);
}

//--------------------------------------------------------------------

void AliAnalysisTaskJetCorePP::UserExec(Option_t *)
{
   //User Exec
   

   //Event loop
   Double_t eventW  = 1.0;
   Double_t ptHard  = 0.0;
   Double_t nTrials = 1.0; // Trials for MC trigger
   if(fIsChargedMC || fIsKine) fh1AvgTrials->Fill("#sum{avg ntrials}",fAvgTrials); 

   if(TMath::Abs((Float_t) fJetParamR)<0.00001){
      AliError("ANTIKT Cone radius is set to zero.");  
      return;
   }

   if(TMath::Abs((Float_t) fBgJetParamR)<0.00001){
      AliError("ANTIKT Cone radius is set to zero.");  
      return;
   }

   if(!fIsKine){ //real data or full MC
      if(!strlen(fJetBranchName.Data())){
         AliError("Name of jet branch not set.");
         return;
      }
      if(!strlen(fJetBranchNameBg.Data())){
         AliError("Name of jet bg branch not set.");
         return;
      }
   }else{  //Kine
       if(!strlen(fJetBranchNameBgKine.Data())){
         AliError("Name of jet bg branch for kine not set.");
         return;
       }

      Init(); 
      if(fMcHandler){
         fMcEvent = fMcHandler->MCEvent(); 
      }else{
         if(fDebug > 1) printf("AliAnalysisTaskJetCorePP::Exec() fMcHandler=NULL\n");
         PostData(1, fOutputList);
         return;
      } 
      if(!fMcEvent){
         if(fDebug > 1) printf("AliAnalysisTaskJetCorePP::Exec() fMcEvent=NULL \n");
         PostData(1, fOutputList);
         return;
      }
 
      Float_t xsection = 0;
      Float_t trials  = 0;

      AliGenPythiaEventHeader *genPH =
         dynamic_cast<AliGenPythiaEventHeader*> (fMcEvent->GenEventHeader()); 
      if(genPH){
         xsection = genPH->GetXsection();
         trials   = genPH->Trials();
         ptHard   = genPH->GetPtHard();
      }
      fh1Xsec->Fill("<#sigma>",xsection);
      fh1Trials->Fill("#sum{ntrials}",trials);
      fh1PtHard->Fill(ptHard,eventW);
      fh1PtHardNoW->Fill(ptHard,1);
      fh1PtHardTrials->Fill(ptHard,trials);
   }


   fESD = dynamic_cast<AliESDEvent*>(InputEvent());
   if(!fESD){
      if(fDebug>1) AliError("ESD not available");
      fAODIn = dynamic_cast<AliAODEvent*>(InputEvent());
   } 
 
   fAODOut = dynamic_cast<AliAODEvent*>(AODEvent());
   AliAODEvent* aod = NULL;
   // take all other information from the aod we take the tracks from
   if(!aod && !fIsKine){
      if(!fESD)         aod = fAODIn; //not ESD and not kine  => input AOD
      else              aod = fAODOut;// ESD or kine
   }


   if(fNonStdFile.Length()!=0){
      // case that we have an AOD extension we can fetch the jets from the extended output
      AliAODHandler *aodH = dynamic_cast<AliAODHandler*>(AliAnalysisManager::GetAnalysisManager()->GetOutputEventHandler());
      fAODExtension = aodH ? aodH->GetExtension(fNonStdFile.Data()) : 0;
      if(!fAODExtension){
         if(fDebug>1) Printf("AODExtension found for %s",fNonStdFile.Data());
      } 
   }
    
   // ----------------- EVENT SELECTION  --------------------------
   fHistEvtSelection->Fill(1); // number of events before event selection

   if(!fIsKine){ 
      // physics selection
      AliInputEventHandler* inputHandler = (AliInputEventHandler*)
           ((AliAnalysisManager::GetAnalysisManager())->GetInputEventHandler());

      if(fOfflineTrgMask > 0){
         if(!(inputHandler->IsEventSelected() & fOfflineTrgMask)){
            if(fDebug) Printf(" Trigger Selection: event REJECTED ... ");
            fHistEvtSelection->Fill(2);
            PostData(1, fOutputList);
            return;
         }
      }
   //check AOD pointer
      if(!aod){
         if(fDebug) Printf("%s:%d No AOD",(char*)__FILE__,__LINE__);
         fHistEvtSelection->Fill(3);
         PostData(1, fOutputList);
         return;
      }

      // vertex selection for reconstructed data
      AliAODVertex* primVtx = aod->GetPrimaryVertex();

      if(!primVtx){
         if(fDebug) Printf("%s:%d No primVtx",(char*)__FILE__,__LINE__);
         fHistEvtSelection->Fill(3);
         PostData(1, fOutputList);
         return;
      }  

      Int_t nTracksPrim = primVtx->GetNContributors();
      Float_t vtxz = primVtx->GetZ();
      //Input events
      fhContribVtx->Fill(nTracksPrim);
      if( nTracksPrim > 0 ) fhVertexZ->Fill(vtxz);

      if((nTracksPrim < fMinContribVtx) ||
         (vtxz < fVtxZMin) ||
         (vtxz > fVtxZMax)){
         if(fDebug) Printf("%s:%d primary vertex z = %f: event REJECTED...",
                         (char*)__FILE__,__LINE__,vtxz);
         fHistEvtSelection->Fill(3);
         PostData(1, fOutputList);
         return;
      }else{
      //Accepted events
         fhContribVtxAccept->Fill(nTracksPrim);
         fhVertexZAccept->Fill(vtxz);
      }
   }else{ //KINE cut on vertex
      const AliVVertex *vtxMC = fMcEvent->GetPrimaryVertex();
      Float_t zVtx = vtxMC->GetZ();
      fhVertexZ->Fill(zVtx);
      if(zVtx < fVtxZMin  || zVtx > fVtxZMax){ //vertex cut
         fHistEvtSelection->Fill(3);
         PostData(1, fOutputList);
         return;
      }else{
         fhVertexZAccept->Fill(zVtx);
      }
   }
   //FK// No event class selection imposed
   // event class selection (from jet helper task)
   //Int_t eventClass = AliAnalysisHelperJetTasks::EventClass();
   //if(fDebug) Printf("Event class %d", eventClass);
   //if(eventClass < fEvtClassMin || eventClass > fEvtClassMax){
   //   fHistEvtSelection->Fill(4);
   //   PostData(1, fOutputList);
   //   return;
   //}

   //------------------ CENTRALITY SELECTION ---------------
   AliCentrality *cent = 0x0;
   Double_t centValue  = 0.0; 
   if(fSystem){  //fSystem=0 for pp,   fSystem=1 for pPb
      if(fESD){
         cent = fESD->GetCentrality();
         if(cent) centValue = cent->GetCentralityPercentile("V0M");
      }else{
         //centValue = aod->GetHeader()->GetCentrality();
         centValue = ((AliVAODHeader*)aod->GetHeader())->GetCentrality();
      }   
      if(fDebug) printf("centrality: %f\n", centValue);
      //Input events
      fhCentrality->Fill(centValue); 

      if(centValue < fCentMin || centValue > fCentMax){
         fHistEvtSelection->Fill(4);
         PostData(1, fOutputList);
         return;
      }else{
         //Accepted events
         fhCentralityAccept->Fill( centValue );
      }
   }
 
   //-----------------select disjunct event subsamples ----------------
   if(!fIsKine){ //reconstructed data
      //Int_t eventnum  = aod->GetHeader()->GetEventNumberESDFile();
      AliAODHeader * header = dynamic_cast<AliAODHeader*>(aod->GetHeader());
      if(!header) AliFatal("Not a standard AOD");

      Int_t eventnum  = header->GetEventNumberESDFile();
      Int_t lastdigit = eventnum % 10;
      if(!(fEventNumberRangeLow<=lastdigit && lastdigit<=fEventNumberRangeHigh)){
         fHistEvtSelection->Fill(5);
         PostData(1, fOutputList);
         return;
       } 
   }  

   if(fDebug) std::cout<<" ACCEPTED EVENT "<<endl;
   fHistEvtSelection->Fill(0); // accepted events 
   // ==================== end event selection ============================ 
 
   Double_t tmpArrayFive[5];
   Double_t tmpArrayFour[4];


   TList particleList; //list of tracks
   Int_t indexTrigg = -1; 
   Double_t rhoFromCellMedian=0.0,    rhoCone=0.0;

   if(!fIsKine){ 
      //=============== Reconstructed antikt jets =============== 
      ReadTClonesArray(fJetBranchName.Data()  , fListJets); 
      ReadTClonesArray(fJetBranchNameBg.Data(), fListJetsBg); 

      //============ Estimate background in reconstructed events ===========

      //Find Hadron trigger and Estimate rho from cone
      indexTrigg = GetListOfTracks(&particleList); //index of trigger hadron in Particle list
      EstimateBgCone(fListJets, &particleList, rhoCone);

      //Estimate rho from cell median minus jets
      EstimateBgRhoMedian(fListJetsBg, &particleList, rhoFromCellMedian,0); //real data

      //Delta Pt
      FillDeltaPt(fListJets, &particleList, rhoFromCellMedian, rhoCone);
   }
   //==============  analyze generator level MC  ================ 
   TList particleListGen; //list of tracks in MC

   if(fIsChargedMC || fIsKine){

      if(fIsKine){  //pure kine

         //================= generated charged antikt jets ================
         ReadTClonesArray(fJetBranchNameKine.Data(),   fListJetsGen); 
         ReadTClonesArray(fJetBranchNameBgKine.Data(), fListJetsBgGen); 
      }else{  
         //================= generated charged antikt jets ================
         ReadTClonesArray(fJetBranchNameChargMC.Data(),   fListJetsGen); 
         ReadTClonesArray(fJetBranchNameBgChargMC.Data(), fListJetsBgGen); 

         if(fIsFullMC){ //generated full jets
            ReadTClonesArray(fJetBranchNameFullMC.Data(),  fListJetsGenFull); 
         }
      }
      //========================================================
      //serarch for charged trigger at the MC generator level

      Int_t    indexTriggGen = -1;
      Double_t ptTriggGen    = -1;
      Int_t    iCounterGen   =  0; //number of entries in particleListGen array
      Int_t    triggersMC[200];//list of trigger candidates
      Int_t    ntriggersMC   = 0; //index in triggers array
      Int_t    triggersMCa[200];   //list of trigger candidates  10%eloss
      Int_t    ntriggersMCa   = 0; //index in triggers array   10%eloss
      Int_t    triggersMCb[200];   //list of trigger candidates  5%eloss
      Int_t    ntriggersMCb   = 0; //index in triggers array     5%eloss


      if(!fIsKine){  
         if(fESD){//ESD input

            AliMCEvent* mcEvent = MCEvent();
            if(!mcEvent){
               PostData(1, fOutputList);
               return;
            }

            AliGenPythiaEventHeader *pythiaGenHeader = AliAnalysisHelperJetTasks::GetPythiaEventHeader(mcEvent);
            if(pythiaGenHeader){
               nTrials = pythiaGenHeader->Trials();
               ptHard  = pythiaGenHeader->GetPtHard();
               fh1PtHard->Fill(ptHard,eventW);
               fh1PtHardNoW->Fill(ptHard,1);
               fh1PtHardTrials->Fill(ptHard,nTrials);
            }
         
            for(Int_t it = 0; it < mcEvent->GetNumberOfTracks(); it++){
               if(!mcEvent->IsPhysicalPrimary(it)) continue;
               AliMCParticle* part = (AliMCParticle*) mcEvent->GetTrack(it);
               if(!part) continue;  
               if(SelectMCGenTracks((AliVParticle*) part, &particleListGen, ptTriggGen, indexTriggGen, iCounterGen)){ 
 
                  if(fHardest==0 && ntriggersMC<200){//single inclusive trigger
                     if(indexTriggGen > -1){//trigger candidate was found
                        triggersMC[ntriggersMC] = indexTriggGen;
                        ntriggersMC++; 
                     }
                  }
 
                  iCounterGen++;//index in particleListGen array 
               }
            }
         }else{  //AOD input
            if(!fAODIn){
               PostData(1, fOutputList);
               return;
            }
            TClonesArray *tca = dynamic_cast<TClonesArray*>(fAODIn->FindListObject(AliAODMCParticle::StdBranchName()));
            if(!tca){ 
               PostData(1, fOutputList);
               return;
            }
            for(Int_t it = 0; it < tca->GetEntriesFast(); it++){
               AliAODMCParticle *part = dynamic_cast<AliAODMCParticle*>(tca->At(it));
               if(!part) continue;  
               if(!part->IsPhysicalPrimary()) continue;
               if(SelectMCGenTracks((AliVParticle*) part, &particleListGen, ptTriggGen, indexTriggGen, iCounterGen)){
 
                  if(fHardest==0 && ntriggersMC<200){//single inclusive trigger
                     if(indexTriggGen > -1){ //trigger candidater was found
                        triggersMC[ntriggersMC] = indexTriggGen;
                        ntriggersMC++; 
                     }
                  }
 
                  iCounterGen++;//number of entries in particleListGen array
               }
            }
         }
      }else{ //analyze kine

         for(Int_t it = 0; it < fMcEvent->GetNumberOfTracks(); it++){
            if(!fMcEvent->IsPhysicalPrimary(it)) continue;
            AliMCParticle* part = (AliMCParticle*) fMcEvent->GetTrack(it);
            if(!part) continue;  
            if(SelectMCGenTracks((AliVParticle*) part, &particleListGen, ptTriggGen, indexTriggGen, iCounterGen)){ 
 
               if(fHardest==0 && ntriggersMC<200){//single inclusive trigger
                  if(indexTriggGen > -1){//trigger candidate was found
                     triggersMC[ntriggersMC] = indexTriggGen;
                     ntriggersMC++; 
                  }
               }
 
               iCounterGen++;//index in particleListGen array 
            }
         }
      }
 
      if(fHardest==0){
         Int_t npar = particleListGen.GetEntries();
         for(Int_t ip=0; ip < npar; ip++){
            AliVParticle *part = (AliVParticle*) particleListGen.At(ip);
            if(!part) continue;
            
            Double_t pta = 0.9 * part->Pt(); //10% energy loss 
            Double_t ptb = 0.95 * part->Pt(); //5% energy loss  
            if(fTriggerPtRangeLow <= pta && pta < fTriggerPtRangeHigh  && 
                TMath::Abs(part->Eta()) < fTriggerEtaCut  &&  ntriggersMCa<200){
               triggersMCa[ntriggersMCa] = ip;
               ntriggersMCa++;
            }

            if(fTriggerPtRangeLow <= ptb && ptb < fTriggerPtRangeHigh  && 
               TMath::Abs(part->Eta()) < fTriggerEtaCut  &&  ntriggersMCb<200){
               triggersMCb[ntriggersMCb] = ip;
               ntriggersMCb++;
            }
         }

         if(ntriggersMCa>0){
            Int_t rnda     = fRandom->Integer(ntriggersMCa); //0 to ntriggers-1
            Int_t indexTriggGena = triggersMCa[rnda];

            Double_t deltaPhia, deltaEtaa, deltaRa;
            Int_t aaLow = 0; 
            Int_t aaHigh = 0; 

            //Correlation with single inclusive  TRIGGER
            AliVParticle* tGenTa = (AliVParticle*) particleListGen.At(indexTriggGena);  
            if(tGenTa){
               fhNofTriggersGenA->Fill(0.5); // 15-50
               //for(Int_t ia=0; ia<ntriggersMCa; ia++)
               for(Int_t ia=0; ia< npar; ia++){
                  //if(indexTriggGena == triggersMCa[ia]) continue;
                  if(indexTriggGena == ia) continue;
               
                  //AliVParticle* tGenTz = (AliVParticle*) particleListGen.At(triggersMCa[ia]);  
                  AliVParticle* tGenTz = (AliVParticle*) particleListGen.At(ia);  
                  if(!tGenTz) continue;
                  if(tGenTz->Pt()*0.9<15.0) continue;
                  if(tGenTz->Pt()*0.9>50.0) continue;
               
                  deltaPhia = RelativePhi(tGenTa->Phi(),tGenTz->Phi());
                  deltaEtaa = tGenTa->Eta()-tGenTz->Eta(); 
                  deltaRa = sqrt(deltaPhia*deltaPhia + deltaEtaa*deltaEtaa);
            
                  if(tGenTz->Pt()*0.9<20.0){
                     fhDeltaRMultTriggersGenALow->Fill(deltaRa);
                     fhDeltaPhiMultTriggersGenALow->Fill(deltaPhia);
                  }else{
                     fhDeltaRMultTriggersGenAHigh->Fill(deltaRa);
                     fhDeltaPhiMultTriggersGenAHigh->Fill(deltaPhia);
                  } 
      
                  if(deltaRa<0.4){ 
                     if(tGenTz->Pt()*0.9<20.0) aaLow++; //15-20
                     else aaHigh++;  //20-50
                  }
               }
            }
            fhNofMultipleTriggersConeGenA->Fill(aaLow+aaHigh); // 15-50
            fhNofMultipleTriggersConeGenALow->Fill(aaLow); //15-20
            fhNofMultipleTriggersConeGenAHigh->Fill(aaHigh);//20-50
         }

         if(ntriggersMCb>0){
            Int_t rndb     = fRandom->Integer(ntriggersMCb); //0 to ntriggers-1
            Int_t indexTriggGenb = triggersMCb[rndb];

            Double_t deltaPhib, deltaEtab, deltaRb;
            Int_t bbLow = 0; 
            Int_t bbHigh = 0; 

            //Correlation with single inclusive  TRIGGER
            AliVParticle* tGenTb = (AliVParticle*) particleListGen.At(indexTriggGenb);  
            if(tGenTb){
               fhNofTriggersGenB->Fill(0.5); // 20-50
               //for(Int_t ib=0; ib<ntriggersMCb; ib++)
               for(Int_t ib=0; ib<npar; ib++){
                  //if(indexTriggGenb == triggersMCb[ib]) continue;
                  if(indexTriggGenb == ib) continue;
               
                  //AliVParticle* tGenTz = (AliVParticle*) particleListGen.At(triggersMCb[ib]);  
                  AliVParticle* tGenTz = (AliVParticle*) particleListGen.At(ib);  
                  if(!tGenTz) continue;
                  if(tGenTz->Pt()*0.95<15.0) continue;
                  if(tGenTz->Pt()*0.95>50.0) continue;
               
                  deltaPhib = RelativePhi(tGenTb->Phi(),tGenTz->Phi());
                  deltaEtab = tGenTb->Eta()-tGenTz->Eta(); 
                  deltaRb = sqrt(deltaPhib*deltaPhib + deltaEtab*deltaEtab);
 
                  if(tGenTz->Pt()*0.95<20.0){
                     fhDeltaRMultTriggersGenBLow->Fill(deltaRb);
                     fhDeltaPhiMultTriggersGenBLow->Fill(deltaPhib);
                  }else{
                     fhDeltaRMultTriggersGenBHigh->Fill(deltaRb);
                     fhDeltaPhiMultTriggersGenBHigh->Fill(deltaPhib);
                  } 
                  
                  if(deltaRb<0.4){
                     if(tGenTz->Pt()*0.95<20.0) bbLow++; //15-20
                     else bbHigh++;
                  }
               }
            }
            fhNofMultipleTriggersConeGenB->Fill(bbLow+bbHigh); //15-50
            fhNofMultipleTriggersConeGenBLow->Fill(bbLow);//15-20
            fhNofMultipleTriggersConeGenBHigh->Fill(bbHigh);//20-50
         }
      }


 
      //==============  Estimate bg in generated events ==============
      Double_t rhoFromCellMedianGen=0.0, rhoConeGen=0.0;

      //Estimate rho from cone
      EstimateBgCone(fListJetsGen, &particleListGen, rhoConeGen);

      //Estimate rho from cell median minus jets
      EstimateBgRhoMedian(fListJetsBgGen, &particleListGen, rhoFromCellMedianGen,1);//mc data

      if(ntriggersMC > 0){ //select random inclusive trigger
         fhInclTrigCounterGen->Fill(0.5,ntriggersMC); //count inclusive triggers

         //Lower pTa bin 15-20
         for(Int_t it1=0; it1<ntriggersMC; it1++){
            AliVParticle* tGent1 = (AliVParticle*) particleListGen.At(triggersMC[it1]);  
            if(!tGent1) continue;
            for(Int_t ia=0; ia<particleListGen.GetEntries(); ia++){
               if(ia == triggersMC[it1]) continue;
               AliVParticle* tGent2 = (AliVParticle*) particleListGen.At(ia);  
               if(!tGent2) continue;
               if(tGent2->Pt()<15.0) continue;
               if(tGent2->Pt()>20.0) continue;
               fhDeltaPhiMultTriggersInclGenLow->Fill(RelativePhi(tGent1->Phi(),tGent2->Phi()));
            } 
         }
         //Higher pTa bin 20-50
         for(Int_t it1=0; it1<ntriggersMC-1; it1++){
            AliVParticle* tGent1 = (AliVParticle*) particleListGen.At(triggersMC[it1]);  
            if(!tGent1) continue;
 
            for(Int_t it2=it1+1; it2<ntriggersMC; it2++){
               AliVParticle* tGent2 = (AliVParticle*) particleListGen.At(triggersMC[it2]);  
               if(!tGent2) continue;
               fhDeltaPhiMultTriggersInclGenHigh->Fill(RelativePhi(tGent1->Phi(),tGent2->Phi()));
            }
         }
      }


      //============  Generator trigger+jet ==================
      if(fHardest==0){ //single inclusive trigger
         if(ntriggersMC>0){ //there is at least one trigger 
            Int_t rnd     = fRandom->Integer(ntriggersMC); //0 to ntriggers-1
            indexTriggGen = triggersMC[rnd];

            fhNofMultipleTriggersGen->Fill(ntriggersMC-1);

            Double_t deltaPhi, deltaEta, deltaR;
            Int_t iLow = 0; 
            Int_t iHigh = 0; 

            //Correlation with single inclusive  TRIGGER
            AliVParticle* tGenT1 = (AliVParticle*) particleListGen.At(indexTriggGen);  
            if(tGenT1){
               //for(Int_t ia=0; ia<ntriggersMC; ia++)
               for(Int_t ia=0; ia< particleListGen.GetEntries(); ia++){
                  //if(indexTriggGen == triggersMC[ia]) continue;
                  if(indexTriggGen == ia) continue;
               
                  //AliVParticle* tGenT2 = (AliVParticle*) particleListGen.At(triggersMC[ia]);  
                  AliVParticle* tGenT2 = (AliVParticle*) particleListGen.At(ia);  
                  if(!tGenT2) continue;
                  if(tGenT2->Pt()<15.0) continue; 
                  if(tGenT2->Pt()>50.0) continue; 
                  deltaPhi = RelativePhi(tGenT1->Phi(),tGenT2->Phi());
                  deltaEta = tGenT1->Eta()-tGenT2->Eta(); 
                  deltaR = sqrt(deltaPhi*deltaPhi + deltaEta*deltaEta);
                  
                  if(tGenT2->Pt()<20.0){  //single inclusive trigger + assoc 15-20
                     fhDeltaRMultTriggersGenLow->Fill(deltaR);
                     fhDeltaPhiMultTriggersGenLow->Fill(deltaPhi);
                  }else{   //single inclusive trigger + assoc 20-50
                     fhDeltaRMultTriggersGenHigh->Fill(deltaR);
                     fhDeltaPhiMultTriggersGenHigh->Fill(deltaPhi);
                  }

                  if(deltaR<0.4){
                     if(tGenT2->Pt()<20.0) iLow++;
                     else iHigh++;
                  }
               }
            }
            fhNofMultipleTriggersConeGen->Fill(iLow+iHigh);
            fhNofMultipleTriggersConeGenLow->Fill(iLow);
            fhNofMultipleTriggersConeGenHigh->Fill(iHigh);

            //Single inclusive trigger within |Delta eta


         }else{
            indexTriggGen = -1; //trigger not found
         }
      }

      //-----------
      Int_t ilowGen  = (fHardest==0 || fHardest==1) ? indexTriggGen : 0;
      Int_t ihighGen = (fHardest==0 || fHardest==1) ? indexTriggGen+1 : particleListGen.GetEntries();
      Bool_t fillOnceGen = kTRUE;
      //-----------
      
      for(Int_t igen= ilowGen; igen< ihighGen; igen++){ //loop over possible trigger
         indexTriggGen = igen; //trigger hadron 

         if(indexTriggGen == -1) continue;  
         AliVParticle* triggerGen=(AliVParticle*) particleListGen.At(indexTriggGen);  
         if(!triggerGen) continue;
 
         if(fHardest >= 2){ 
            if(triggerGen->Pt() < 6.0) continue; //all hadrons pt>6  
         }
         if(TMath::Abs((Float_t) triggerGen->Eta()) > fTriggerEtaCut) continue;

         ptTriggGen = triggerGen->Pt(); //count triggers
         fh2NtriggersGen->Fill(centValue, ptTriggGen);

         //Count jets and trigger-jet pairs at MC  generator level
         for(Int_t ij=0; ij<fListJetsGen->GetEntries(); ij++){
            AliAODJet* jet = (AliAODJet*)(fListJetsGen->At(ij));
            if(!jet) continue;
            Double_t etaJetGen = jet->Eta();
            Double_t areaJetGen = jet->EffectiveAreaCharged();
 
            if((fJetEtaMin<=etaJetGen) && (etaJetGen<=fJetEtaMax)){ 

               Double_t ptUeFromCellMedianGen = rhoFromCellMedianGen*areaJetGen;
               Double_t ptUeConeGen         = rhoConeGen*areaJetGen;

               //Rongrong's analysis
               Double_t dPhiGen = jet->Phi() - triggerGen->Phi();
               if(dPhiGen>2*TMath::Pi())    dPhiGen -= 2*TMath::Pi();
               if(dPhiGen<-2*TMath::Pi())   dPhiGen += 2*TMath::Pi();
               if(dPhiGen<-0.5*TMath::Pi()) dPhiGen += 2*TMath::Pi();
               if(dPhiGen>1.5*TMath::Pi())  dPhiGen -= 2*TMath::Pi();

               //[A,pTjet,pTtrig,dphi]
               tmpArrayFour[0] = areaJetGen;
               tmpArrayFour[1] = jet->Pt();
               tmpArrayFour[2] = ptTriggGen;
               tmpArrayFour[3] = dPhiGen;

               fHJetPhiCorrGen->Fill(tmpArrayFour); // Dphi distribution jet-triger

               //[A,pTjet-pTUe,pTtrig,dphi]
               tmpArrayFour[1] = jet->Pt() - ptUeFromCellMedianGen;
               fHJetPhiCorrSubUeMedianGen->Fill(tmpArrayFour); 

               //[A,pTjet-pTUe,pTtrig,dphi]
               tmpArrayFour[1] = jet->Pt() - ptUeConeGen;
               fHJetPhiCorrSubUeConeGen->Fill(tmpArrayFour);

               //Leticia's analysis
               Double_t dphi = RelativePhi(triggerGen->Phi(), jet->Phi()); 
               if(TMath::Abs((Double_t) dphi) < fkDeltaPhiCut) continue;

               //Centrality, A, pT, pTtrigg
               tmpArrayFive[0] = centValue;
               tmpArrayFive[1] = areaJetGen;
               tmpArrayFive[2] = jet->Pt();
               tmpArrayFive[3] = ptTriggGen;
               tmpArrayFive[4] = TMath::Abs((Double_t) dphi);
               fHJetSpecGen->Fill(tmpArrayFive);

               //Centrality, A, pTjet-pTbgCellMedian, pTtrigg, dphi
               tmpArrayFive[0] = centValue;
               tmpArrayFive[1] = areaJetGen;
               tmpArrayFive[2] = jet->Pt() - ptUeFromCellMedianGen;
               tmpArrayFive[3] = ptTriggGen;
               tmpArrayFive[4] = TMath::Abs((Double_t) dphi);
               fHJetSpecSubUeMedianGen->Fill(tmpArrayFive);

               //Centrality, A, pTjet-pTbgCone, pTtrigg, dphi
               tmpArrayFive[0] = centValue;
               tmpArrayFive[1] = areaJetGen;
               tmpArrayFive[2] = jet->Pt() - ptUeConeGen;
               tmpArrayFive[3] = ptTriggGen;
               tmpArrayFive[4] = TMath::Abs((Double_t) dphi);
               fHJetSpecSubUeConeGen->Fill(tmpArrayFive);

               //Ue diagnostics  "[A,pTjet,pTjet-pTUe, pTUe, rhoUe]", kT median
               tmpArrayFive[0] = areaJetGen;
               tmpArrayFive[1] = jet->Pt();
               tmpArrayFive[2] = jet->Pt() - ptUeFromCellMedianGen;
               tmpArrayFive[3] = ptUeFromCellMedianGen;
               tmpArrayFive[4] = rhoFromCellMedianGen;
               fHJetUeMedianGen->Fill(tmpArrayFive); 

               //Ue diagnostics  "[A,pTjet,pTjet-pTUe, pTUe, rhoUe]", perpendicular Cone
               tmpArrayFive[0] = areaJetGen;
               tmpArrayFive[1] = jet->Pt();
               tmpArrayFive[2] = jet->Pt() - ptUeConeGen;
               tmpArrayFive[3] = ptUeConeGen;
               tmpArrayFive[4] = rhoConeGen;
               fHJetUeConeGen->Fill(tmpArrayFive);

               if(fillOnceGen){
                  Double_t fillRhoGen[] = { rhoConeGen, rhoFromCellMedianGen};
                  fHRhoUeMedianVsConeGen->Fill(fillRhoGen); 
                  fillOnceGen = kFALSE;
               }
            }//back to back jet-trigger pair
         }//jet loop
      }//trigger loop

      if(fFillRespMx && !fIsKine){
         //================ RESPONSE MATRIX ===============
         //Count jets and trigger-jet pairs at MC  generator level
         for(Int_t ij=0; ij<fListJetsGen->GetEntries(); ij++){
            AliAODJet* jet = (AliAODJet*)(fListJetsGen->At(ij));
            if(!jet) continue;
            Double_t etaJetGen = jet->Eta();
            Double_t ptJetGen  = jet->Pt();
       
            if((fJetEtaMin<=etaJetGen) && (etaJetGen<=fJetEtaMax)){ 
               fhJetPtGen->Fill(ptJetGen); // gen level pt spectum of jets response mx normalization

               Double_t areaJetGen = jet->EffectiveAreaCharged();
               Double_t ptUeFromCellMedianGen = rhoFromCellMedianGen*areaJetGen;
               Double_t ptUeConeGen         = rhoConeGen*areaJetGen;

               fhJetPtSubUeMedianGen->Fill(ptJetGen - ptUeFromCellMedianGen); 
               fhJetPtSubUeConeGen->Fill(ptJetGen   - ptUeConeGen);    
            }
         }
         if(fListJets->GetEntries()>0 && fListJetsGen->GetEntries()>0){ //at least some reconstructed jets
         
            Int_t ng = (Int_t) fListJetsGen->GetEntries();
            Int_t nr = (Int_t) fListJets->GetEntries();

            //Find closest MC generator - reconstructed jets
            if(faGenIndex.GetSize()<nr) faGenIndex.Set(nr); //idx of gen jet assoc to rec jet
            if(faRecIndex.GetSize()<ng) faRecIndex.Set(ng); //idx of rec jet assoc to gen jet

            if(fDebug){
               Printf("New Rec List %d gen index Array %d",nr,faGenIndex.GetSize());
               Printf("New Gen List %d rec index Array %d",ng,faRecIndex.GetSize());
            }
            //matching of MC genrator level and reconstructed jets
            AliAnalysisHelperJetTasks::GetClosestJets(fListJetsGen,ng,fListJets,nr,faGenIndex,faRecIndex,fDebug); 

            // Fill response matrix
            for(Int_t ir = 0; ir < nr; ir++){
               AliAODJet *recJet   = (AliAODJet*) fListJets->At(ir);
               Double_t etaJetRec  = recJet->Eta();
               Double_t ptJetRec   = recJet->Pt();
               Double_t areaJetRec = recJet->EffectiveAreaCharged();
               //fill response matrix if generator and reconstructed jets are within |eta|<0.9-fiduc

               if((fJetEtaMin <= etaJetRec) && (etaJetRec <= fJetEtaMax)){ 
                  Int_t ig = faGenIndex[ir]; //associated generator level jet
                  if(ig >= 0 && ig < ng){
                     if(fDebug > 10) Printf("%s:%d ig = %d ir = %d",(char*)__FILE__,__LINE__,ig,ir);
                     AliAODJet *genJet  = (AliAODJet*) fListJetsGen->At(ig);
                     Double_t ptJetGen  = genJet->Pt();
                     Double_t etaJetGen = genJet->Eta();

                     //fill response matrix if generator and reconstructed jets are within |eta|<0.9-fiduc
                     if((fJetEtaMin <= etaJetGen) && (etaJetGen <= fJetEtaMax)){
                        fhJetPtGenVsJetPtRec->Fill(ptJetRec, ptJetGen);
                        if(ptJetGen>0){   
                           fhJetPtResolutionVsPtGen->Fill(ptJetGen,(ptJetRec-ptJetGen)/ptJetGen); 
                        }
                        Double_t areaJetGen  = genJet->EffectiveAreaCharged();
                        Double_t ptUeFromCellMedianGen = rhoFromCellMedianGen*areaJetGen;
                        Double_t ptUeConeGen         = rhoConeGen*areaJetGen;
                        Double_t ptUeFromCellMedianRec = rhoFromCellMedian*areaJetRec;
                        Double_t ptUeConeRec         = rhoCone*areaJetRec;
                        fhJetPtGenVsJetPtRecSubUeMedian->Fill(ptJetRec-ptUeFromCellMedianRec, 
                                                           ptJetGen-ptUeFromCellMedianGen);
                        fhJetPtGenVsJetPtRecSubUeCone->Fill(ptJetRec-ptUeConeRec, ptJetGen-ptUeConeGen);

                        if((ptJetGen-ptUeConeGen)>0){   
                           fhJetPtResolutionVsPtConeGen->Fill(ptJetGen-ptUeConeGen,((ptJetRec-ptUeConeRec)- (ptJetGen-ptUeConeGen))/ (ptJetGen-ptUeConeGen)); 
                        }
                     }
                  }//ig>=0
               }//rec jet in eta acceptance
            }//loop over reconstructed jets
         }// # of  rec jets >0
      
         //=========================== Full jet vs charged jet matrix  ==========================
         if(fIsFullMC){
            //Count full jets and charged-jet pairs at MC  generator level
            for(Int_t ij=0; ij<fListJetsGenFull->GetEntries(); ij++){
               AliAODJet* jetFull = (AliAODJet*)(fListJetsGenFull->At(ij));
               if(!jetFull) continue;
               Double_t etaJetFull = jetFull->Eta();
               Double_t ptJetFull  = jetFull->Pt();

               if((fJetEtaMin<=etaJetFull) && (etaJetFull<=fJetEtaMax)){
                  fhJetPtGenFull->Fill(ptJetFull); // generator level pt spectum of full jets 
               }
            }
            if(fListJetsGen->GetEntries()>0 && fListJetsGenFull->GetEntries()>0){ //at least some reconstructed jets
               Int_t nful = (Int_t) fListJetsGenFull->GetEntries();
               Int_t nchr = (Int_t) fListJetsGen->GetEntries();
 
               //Find closest MC generator full - charged jet
               if(faGenIndex.GetSize()<nchr) faGenIndex.Set(nchr); //idx of gen FULL jet assoc to gen CHARGED jet
               if(faRecIndex.GetSize()<nful) faRecIndex.Set(nful); //idx of gen CHARGED jet assoc to gen FULL jet
 
               if(fDebug){
                  Printf("New Charg List %d Full index Array %d",nchr,faGenIndex.GetSize());
                  Printf("New Full List %d Charg index Array %d",nful,faRecIndex.GetSize());
               }
               //matching of MC genrator level and reconstructed jets
               AliAnalysisHelperJetTasks::GetClosestJets(fListJetsGenFull,nful,fListJetsGen,nchr,faGenIndex,faRecIndex,fDebug);

              // Fill response matrix
               for(Int_t ichr = 0; ichr < nchr; ichr++){ //charged jet loop
                  AliAODJet *chJet  = (AliAODJet*) fListJetsGen->At(ichr);
                  Double_t etaJetCh = chJet->Eta();
                  Double_t ptJetCh  = chJet->Pt();
                  //fill response matrix if generator and reconstructed jets are within |eta|<0.9-fiduc
    
                  if((fJetEtaMin <= etaJetCh) && (etaJetCh <= fJetEtaMax)){
                     Int_t iful = faGenIndex[ichr]; //associated generator level jet
                     if(iful >= 0 && iful < nful){
                        if(fDebug > 10) Printf("%s:%d iful = %d ichr = %d",(char*)__FILE__,__LINE__,iful,ichr);
                        AliAODJet *genJetFull  = (AliAODJet*) fListJetsGenFull->At(iful);
                        Double_t ptJetFull  = genJetFull->Pt();
                        Double_t etaJetFull = genJetFull->Eta();

                        //fill response matrix if generator and reconstructed jets are within |eta|<0.9-fiduc
                        if((fJetEtaMin <= etaJetFull) && (etaJetFull <= fJetEtaMax)){
                           fhJetPtGenChargVsJetPtGenFull->Fill(ptJetFull,ptJetCh);
                        }
                     }//iful>=0
                  }//rec jet in eta acceptance
               }//loop over reconstructed jets
            }// # of  rec jets >0
         }//pointer MC generator jets
      } //fill resp mx only for bin 
   }//analyze generator level MC
    
    
   if(fIsKine){ //skip reconstructed data analysis in case of kine 
      PostData(1, fOutputList);
      return;
   } 
   //=============  RECONSTRUCTED INCLUSIVE JETS ===============

   Double_t etaJet  = 0.0;
   Double_t pTJet   = 0.0;
   Double_t areaJet = 0.0;
   Double_t phiJet  = 0.0;
   //Int_t indexLeadingJet     = -1;
   //Double_t pTLeadingJet     = -10.0; 
   //Double_t areaLeadingJet   = -10.0;
  
   for(Int_t ij=0; ij<fListJets->GetEntries(); ij++){
      AliAODJet* jet = (AliAODJet*)(fListJets->At(ij));
      if(!jet) continue;
      etaJet  = jet->Eta();
      phiJet  = jet->Phi();
      pTJet   = jet->Pt();
      if(pTJet==0) continue; 
     
      if((etaJet<fJetEtaMin) || (etaJet>fJetEtaMax)) continue;
      /*areaJet = jet->EffectiveAreaCharged();*/

      //Jet Diagnostics---------------------------------
      fhJetPhi->Fill(pTJet, RelativePhi(phiJet,0.0)); //phi -pi,pi
      fhJetEta->Fill(pTJet, etaJet);
      //search for leading jet
      /*if(pTJet > pTLeadingJet){
         indexLeadingJet  = ij; 
         pTLeadingJet     = pTJet; 
         areaLeadingJet   = areaJet; 
      } 
 
      // raw spectra of INCLUSIVE jets  
      //Centrality, pTjet, A
      Double_t fillraw[] = { centValue,
                             pTJet,
                             areaJet
                           };
      fHJetPtRaw->Fill(fillraw);*/
   }
   /*
   if(indexLeadingJet > -1){ 
      // raw spectra of LEADING jets  
      //Centrality, pTjet,  A
      Double_t fillleading[] = { centValue,
                                 pTLeadingJet,
                                 areaLeadingJet
                               };
      fHLeadingJetPtRaw->Fill(fillleading);
   } 
   */
   //================  TWO PARTICLE CORRELATION CORRECTION ===============
   if(fIsChargedMC && fHardest==0 ){
      if(indexTrigg>=0){   // Reconstructed trigger
         AliVParticle *triggerHadron = (AliVParticle*) particleList.At(indexTrigg);    
         Int_t genTriggIndex = -1; 
         if(triggerHadron){ 
            Int_t trigLabel = TMath::Abs(triggerHadron->GetLabel());

            Int_t nGen = particleListGen.GetEntries();
            for(Int_t ig=0; ig<nGen; ig++){
	       AliVParticle *trkGen = (AliVParticle*) particleListGen.At(ig);  
               if(!trkGen) continue;

               Int_t genLabel = TMath::Abs(trkGen->GetLabel()); 
               if(trigLabel==genLabel){
                  if(fTriggerPtRangeLow <= trkGen->Pt() && 
                     trkGen->Pt() < fTriggerPtRangeHigh &&
                     TMath::Abs((Float_t) trkGen->Eta()) < fTriggerEtaCut){ 
                     genTriggIndex = ig; //There is a gen level particle track matching to rec trigger
                     break;
                  }
               }
            }
            if(genTriggIndex>-1){  // make generator level TT x Assoc
	       AliVParticle *trkGenTT = (AliVParticle*) particleListGen.At(genTriggIndex);  
               if(trkGenTT){
                  fhTriggerCounterGenLevel->Fill(0.5); // count gen triggers that have reconstructed analog
                  Double_t deltaPhi,deltaEta,deltaR;

                  for(Int_t ig=0; ig<nGen; ig++){
                     if(ig==genTriggIndex) continue; //skip trigger
 	             AliVParticle *trkGenA = (AliVParticle*) particleListGen.At(ig);  
                     if(!trkGenA) continue;
                     if(trkGenA->Pt()<15.0) continue;
                     if(trkGenA->Pt()>50.0) continue;
                     deltaPhi = RelativePhi(trkGenTT->Phi(), trkGenA->Phi());
                     deltaEta = trkGenTT->Eta() - trkGenA->Eta(); 
                     deltaR   = sqrt(deltaPhi*deltaPhi + deltaEta*deltaEta);

                     if(trkGenA->Pt()<20.0){ // gen level assoc particles 
                        fhDeltaRMultTriggersGenLevelLow->Fill(deltaR);
                        fhDeltaPhiMultTriggersGenLevelLow->Fill(deltaPhi);
                     }else{
                        fhDeltaRMultTriggersGenLevelHigh->Fill(deltaR);
                        fhDeltaPhiMultTriggersGenLevelHigh->Fill(deltaPhi);
                     }  
                  }
               }
            }
         }
      }
   }  
 
   // ===============  RECONSTRUCTED TRIGGER-JET PAIRS ================
   if(fIsChargedMC && fFillRespMx){
     FillEffHistos(&particleList, &particleListGen); //Fill efficiency histos
   }
   Bool_t filledOnce = kTRUE; //fill rho histogram only once per event
   //set ranges of the trigger loop
   Int_t ilow  = (fHardest==0 || fHardest==1) ? indexTrigg : 0;
   Int_t ihigh = (fHardest==0 || fHardest==1) ? indexTrigg+1 : particleList.GetEntries();

   for(Int_t itrk= ilow; itrk< ihigh; itrk++){ //loop over possible trigger
      indexTrigg = itrk; //trigger hadron with pT >10 GeV 
 
      if(indexTrigg < 0) continue;

      AliVParticle *triggerHadron = (AliVParticle*) particleList.At(indexTrigg);     
      if(!triggerHadron){  
         PostData(1, fOutputList);
         return;
      }
      if(fHardest >= 2){ 
         if(triggerHadron->Pt() < 6.0) continue; //all hadrons pt>6  
      }
      if(TMath::Abs((Float_t) triggerHadron->Eta())> fTriggerEtaCut) continue;
 
      //Fill trigger histograms
      fh2Ntriggers->Fill(centValue,triggerHadron->Pt()); //trigger pT
      fhTriggerPhi->Fill(triggerHadron->Pt(),RelativePhi(triggerHadron->Phi(),0.0)); //phi -pi,pi
      fhTriggerEta->Fill(triggerHadron->Pt(),triggerHadron->Eta());

  
      //---------- make trigger-jet pairs ---------
      //Int_t injet4     = 0;
      //Int_t injet      = 0; 

      for(Int_t ij=0; ij<fListJets->GetEntries(); ij++){
         AliAODJet* jet = (AliAODJet*)(fListJets->At(ij));
         if(!jet) continue;
         etaJet  = jet->Eta();
         phiJet  = jet->Phi();
         pTJet   = jet->Pt();
         if(pTJet==0) continue; 
     
         if((etaJet<fJetEtaMin) || (etaJet>fJetEtaMax)) continue;
         areaJet = jet->EffectiveAreaCharged();

        //subtract bg using estinates base on median of kT jets
         Double_t ptUeFromCellMedian = rhoFromCellMedian*areaJet;
         Double_t ptUeCone         = rhoCone*areaJet;

         //if(areaJet >= 0.07) injet++; 
         //if(areaJet >= 0.4)  injet4++;

         Double_t dphi = RelativePhi(triggerHadron->Phi(), phiJet); 
         fhDphiTriggerJet->Fill(dphi); //Input

         //Dphi versus jet pT   
         //Centrality, Dphi=phiTrig-phiJet, pTjet, pTtrigg 
         /*Double_t filldp[] = { centValue,
                               dphi,
                               pTJet,
                               triggerHadron->Pt()
                             };
         fHDphiVsJetPtAll->Fill(filldp);
         */ 
         //Rongrong's analysis

         Double_t dPhi = phiJet - triggerHadron->Phi();
         if(dPhi>2*TMath::Pi())    dPhi -= 2*TMath::Pi();
         if(dPhi<-2*TMath::Pi())   dPhi += 2*TMath::Pi();
         if(dPhi<-0.5*TMath::Pi()) dPhi += 2*TMath::Pi();
         if(dPhi>1.5*TMath::Pi())  dPhi -= 2*TMath::Pi();

         //[A,pTjet,pTtrig,dphi]
         tmpArrayFour[0] = areaJet;
         tmpArrayFour[1] = pTJet;
         tmpArrayFour[2] = triggerHadron->Pt();
         tmpArrayFour[3] = dPhi;

         fHJetPhiCorr->Fill(tmpArrayFour); // Dphi distribution jet-triger

         //[A,pTjet-pTUe,pTtrig,dphi]
         tmpArrayFour[1] = pTJet - ptUeFromCellMedian;
         fHJetPhiCorrSubUeMedian->Fill(tmpArrayFour); 

         //[A,pTjet-pTUe,pTtrig,dphi]
         tmpArrayFour[1] = pTJet - ptUeCone;
         fHJetPhiCorrSubUeCone->Fill(tmpArrayFour);

         //Leticia's analysis

         // Select back to back trigger - jet pairs
         if(TMath::Abs((Double_t) dphi) < fkDeltaPhiCut) continue;
         fhDphiTriggerJetAccept->Fill(dphi); //Accepted
          
         //NO bg subtraction
         //Centrality, A, pTjet, pTtrigg, dphi
         tmpArrayFive[0] = centValue;
         tmpArrayFive[1] = areaJet;
         tmpArrayFive[2] = pTJet;
         tmpArrayFive[3] = triggerHadron->Pt();
         tmpArrayFive[4] = TMath::Abs((Double_t) dphi);
         fHJetSpec->Fill(tmpArrayFive);

          //Centrality, A, pTjet-pTbgCellMedian, pTtrigg, dphi
         tmpArrayFive[0] = centValue;
         tmpArrayFive[1] = areaJet;
         tmpArrayFive[2] = pTJet-ptUeFromCellMedian;
         tmpArrayFive[3] = triggerHadron->Pt();
         tmpArrayFive[4] = TMath::Abs((Double_t) dphi);
         fHJetSpecSubUeMedian->Fill(tmpArrayFive);

         //Centrality, A, pTjet-pTbgCone, pTtrigg, dphi
         tmpArrayFive[0] = centValue;
         tmpArrayFive[1] = areaJet;
         tmpArrayFive[2] = pTJet - ptUeCone;
         tmpArrayFive[3] = triggerHadron->Pt();
         tmpArrayFive[4] = TMath::Abs((Double_t) dphi);
         fHJetSpecSubUeCone->Fill(tmpArrayFive);

         //Ue diagnostics  "[A,pTjet,pTjet-pTUe, pTUe, rhoUe]", kT median
         tmpArrayFive[0] = areaJet;
         tmpArrayFive[1] = pTJet;
         tmpArrayFive[2] = pTJet - ptUeFromCellMedian;
         tmpArrayFive[3] = ptUeFromCellMedian;
         tmpArrayFive[4] = rhoFromCellMedian;
         fHJetUeMedian->Fill(tmpArrayFive);
 
         //Ue diagnostics  "[A,pTjet,pTjet-pTUe, pTUe, rhoUe]", Cone median
         tmpArrayFive[0] = areaJet;
         tmpArrayFive[1] = pTJet;
         tmpArrayFive[2] = pTJet - ptUeCone;
         tmpArrayFive[3] = ptUeCone;
         tmpArrayFive[4] = rhoCone;
         fHJetUeCone->Fill(tmpArrayFive);

         if(filledOnce){ //fill for each event only once
            Double_t fillRho[] = { rhoCone,rhoFromCellMedian};
            fHRhoUeMedianVsCone->Fill(fillRho);
            filledOnce = kFALSE;
         }
      }//jet loop
 
      //Fill Jet Density In the Event A>0.07
      /*if(injet>0){
         Double_t filldens[]={ (Double_t) particleList.GetEntries(),
                            injet/fkAcceptance,
                            triggerHadron->Pt()
                          };
         fHJetDensity->Fill(filldens);
      } 

      //Fill Jet Density In the Event A>0.4
      if(injet4>0){ 
         Double_t filldens4[]={ (Double_t) particleList.GetEntries(), 
                                injet4/fkAcceptance,
                                triggerHadron->Pt()
                              };
         fHJetDensityA4->Fill(filldens4);
      }*/
   }


   PostData(1, fOutputList);
}

//----------------------------------------------------------------------------
void AliAnalysisTaskJetCorePP::Terminate(const Option_t *)
{
   // Draw result to the screen
   // Called once at the end of the query

   if(fDebug) printf("AliAnalysisTaskJetCorePP DONE\n");
   if(!GetOutputData(1)) return;
}


//----------------------------------------------------------------------------
Int_t  AliAnalysisTaskJetCorePP::GetListOfTracks(TList *list){
   //Fill the list of accepted tracks (passed track cut)
   //return consecutive index of the hardest ch hadron in the list
   Int_t iCount        = 0;
   AliAODEvent *aodevt = NULL;

   if(!fESD) aodevt = fAODIn;
   else      aodevt = fAODOut;   

   if(!aodevt) return -1;

   Int_t    index = -1; //index of the highest particle in the list
   Double_t ptmax = -10;
   Int_t    triggers[200];
   Int_t    ntriggers = 0; //index in triggers array
 

   for(Int_t it = 0; it < aodevt->GetNumberOfTracks(); it++){
      //AliAODTrack *tr = aodevt->GetTrack(it);
      AliAODTrack *tr = dynamic_cast<AliAODTrack*>(aodevt->GetTrack(it));
      if(!tr) AliFatal("Not a standard AOD");

      if((fFilterMask > 0) && !(tr->TestFilterBit(fFilterMask))) continue;
      //if((fFilterMask > 0) && !(tr->IsHybridGlobalConstrainedGlobal())) continue;
      if(TMath::Abs((Float_t) tr->Eta()) > fTrackEtaCut) continue;
      if(tr->Pt() < fTrackLowPtCut) continue;
      list->Add(tr);

      //Search trigger:
      if(fHardest>0){ //leading particle 
         if(tr->Pt()>ptmax &&
            TMath::Abs((Float_t) tr->Eta()) < fTriggerEtaCut){ 
            ptmax = tr->Pt();	
            index = iCount;
         }
      }
    
      if(fHardest==0 && ntriggers<200){ //single inclusive 
         if(fTriggerPtRangeLow <= tr->Pt() && 
            tr->Pt() < fTriggerPtRangeHigh &&
             TMath::Abs((Float_t) tr->Eta()) < fTriggerEtaCut){ 
            triggers[ntriggers] = iCount;
            ntriggers++;
         }
      }

      iCount++;
   }
   //inclusive trigger
   if(ntriggers>0){ //select random inclusive trigger
      fhInclTrigCounter->Fill(0.5,ntriggers); //count inclusive triggers

      //Lower pTa bin
      for(Int_t it1=0; it1<ntriggers; it1++){
         AliVParticle* tGent1 = (AliVParticle*) list->At(triggers[it1]);  
         if(!tGent1) continue;
         for(Int_t ia=0; ia<list->GetEntries(); ia++){
            if(ia == triggers[it1]) continue;
            AliVParticle* tGent2 = (AliVParticle*) list->At(ia);  
            if(!tGent2) continue;
            if(tGent2->Pt()<15.0) continue;
            if(tGent2->Pt()>20.0) continue;
            fhDeltaPhiMultTriggersInclLow->Fill(RelativePhi(tGent1->Phi(),tGent2->Phi()));
         } 
      }
      //Higher pTa bin
      for(Int_t it1=0; it1<ntriggers-1; it1++){
         AliVParticle* tGent1 = (AliVParticle*) list->At(triggers[it1]);  
         if(!tGent1) continue;
 
         for(Int_t it2=it1+1; it2<ntriggers; it2++){
            AliVParticle* tGent2 = (AliVParticle*) list->At(triggers[it2]);  
            if(!tGent2) continue;
            fhDeltaPhiMultTriggersInclHigh->Fill(RelativePhi(tGent1->Phi(),tGent2->Phi()));
         }
      }
   }

   //single inclusive trigger
   if(fHardest==0 && ntriggers>0){ //select random inclusive trigger
      Int_t rnd = fRandom->Integer(ntriggers); //0 to ntriggers-1
      index     = triggers[rnd];

      fhNofMultipleTriggers->Fill(ntriggers-1);
     
      Double_t deltaPhi, deltaEta, deltaR;
      Int_t iLow=0;
      Int_t iHigh=0;
      //Correlation with single inclusive trigger
      AliVParticle* tGent1 = (AliVParticle*) list->At(index);  
      if(tGent1){
         //for(Int_t ia=0; ia<ntriggers; ia++)
         for(Int_t ia=0; ia<list->GetEntries(); ia++){
            //if(triggers[ia]==index) continue;
            if(ia==index) continue;
            //AliVParticle* tGent2 = (AliVParticle*) list->At(triggers[ia]);  
            AliVParticle* tGent2 = (AliVParticle*) list->At(ia);  
            if(!tGent2) continue;
            if(tGent2->Pt()<15.0) continue;
            if(tGent2->Pt()>50.0) continue;
            deltaPhi = RelativePhi(tGent1->Phi(),tGent2->Phi());
            deltaEta = tGent1->Eta()-tGent2->Eta(); 
            deltaR   = sqrt(deltaPhi*deltaPhi + deltaEta*deltaEta);

            if(tGent2->Pt()<20.0){ 
               fhDeltaRMultTriggersLow->Fill(deltaR);
               fhDeltaPhiMultTriggersLow->Fill(deltaPhi);
            }else{
               fhDeltaRMultTriggersHigh->Fill(deltaR);
               fhDeltaPhiMultTriggersHigh->Fill(deltaPhi);
            }
 
            if(deltaR<0.4){
               if(tGent2->Pt()<20.0) iLow++;
               else iHigh++;
            }
         } 
      }
      fhNofMultipleTriggersCone->Fill(iLow+iHigh);
      fhNofMultipleTriggersConeLow->Fill(iLow);
      fhNofMultipleTriggersConeHigh->Fill(iHigh);

   }

   return index;
}

//----------------------------------------------------------------------------
/*
Double_t AliAnalysisTaskJetCorePP::GetBackgroundInPerpCone(Float_t jetR, Double_t jetPhi, Double_t jetEta, TList* trkList){
   //calculate sum of track pT in the cone perpendicular in phi to the jet 
   //jetR = cone radius
   // jetPhi, jetEta = direction of the jet 
   Int_t numberOfTrks = trkList->GetEntries();
   Double_t pTsum = 0.0;
   Double_t perpConePhi = jetPhi + TMath::Pi()/2;//perp cone w.r.t. jet in phi
   for(Int_t it=0; it<numberOfTrks; it++){
      AliVParticle *trk = (AliVParticle*) trkList->At(it); 
      Double_t dphi = RelativePhi(perpConePhi,trk->Phi());     
      Double_t deta = trk->Eta()-jetEta;     

      if( (dphi*dphi + deta*deta)< (jetR*jetR)){
         pTsum += trk->Pt();
      } 
   }

   return pTsum;
}
*/
//----------------------------------------------------------------------------

Double_t AliAnalysisTaskJetCorePP::RelativePhi(Double_t mphi,Double_t vphi){
   //Get relative azimuthal angle of two particles -pi to pi
   if      (vphi < -TMath::Pi()) vphi += TMath::TwoPi();
   else if (vphi > TMath::Pi())  vphi -= TMath::TwoPi();

   if      (mphi < -TMath::Pi()) mphi += TMath::TwoPi();
   else if (mphi > TMath::Pi())  mphi -= TMath::TwoPi();

   Double_t dphi = mphi - vphi;
   if      (dphi < -TMath::Pi()) dphi += TMath::TwoPi();
   else if (dphi > TMath::Pi())  dphi -= TMath::TwoPi();

   return dphi;//dphi in [-Pi, Pi]
}


//----------------------------------------------------------------------------
Bool_t AliAnalysisTaskJetCorePP::SelectMCGenTracks(AliVParticle *trk, TList *trkList, Double_t &ptLeading, Int_t &index, Int_t counter){
   //fill track efficiency denominator
   if(trk->Pt() < fTrackLowPtCut) return kFALSE;
   if(trk->Charge()==0)        return kFALSE;
   if(TMath::Abs((Float_t) trk->Eta()) > fTrackEtaCut) return kFALSE;
   trkList->Add(trk);
   if(fFillRespMx) fhPtTrkTruePrimGen->Fill(trk->Pt(),trk->Eta()); //Efficiency denominator

   //Search trigger:
   if(fHardest>0){ //leading particle 
      if(ptLeading < trk->Pt() &&
         TMath::Abs((Float_t) trk->Eta()) < fTriggerEtaCut){
         index      = counter;
         ptLeading  = trk->Pt();
      }
   }

   if(fHardest==0){ //single inclusive 
      index = -1;  
      if(fTriggerPtRangeLow <= trk->Pt() && 
            trk->Pt() < fTriggerPtRangeHigh &&
             TMath::Abs((Float_t) trk->Eta()) < fTriggerEtaCut){
         index      = counter;
      }
   }

   return kTRUE;
} 

//----------------------------------------------------------------------------
void AliAnalysisTaskJetCorePP::FillEffHistos(TList *recList, TList *genList){
   //fill track efficiency numerator 
   if(!recList) return;
   if(!genList) return;
   Int_t nRec = recList->GetEntries();
   Int_t nGen = genList->GetEntries();
   for(Int_t ir=0; ir<nRec; ir++){ 
      AliAODTrack *trkRec = (AliAODTrack*) recList->At(ir);
      if(!trkRec) continue; 
      Int_t recLabel = TMath::Abs(trkRec->GetLabel());
      Bool_t matched = kFALSE;
 
      for(Int_t ig=0; ig<nGen; ig++){
	AliVParticle *trkGen = (AliVParticle*) genList->At(ig);  
        if(!trkGen) continue;
        Int_t genLabel = TMath::Abs(trkGen->GetLabel()); 
        if(recLabel==genLabel){
          fhPtTrkTruePrimRec->Fill(trkGen->Pt(),trkGen->Eta());
          matched = kTRUE;
          break;
        }
      }
      if(!matched) fhPtTrkSecOrFakeRec->Fill(trkRec->Pt(),trkRec->Eta()); 
   }

   return;
}
//________________________________________________________________
void AliAnalysisTaskJetCorePP::EstimateBgRhoMedian(TList *listJet, TList* listPart,  Double_t &rhoMedian, Int_t mode){
   //Estimate background rho by means of integrating track pT outside identified jet cones

   rhoMedian = 0;

   //identify leading jet in the full track acceptance
   Int_t idxLeadingJet    = -1;
   Double_t pTleading     = 0.0; 
   AliAODJet* jet = NULL;
 
   for(Int_t ij = 0; ij < listJet->GetEntries(); ij++){ //loop over bg kT jets
      jet = (AliAODJet*)(listJet->At(ij));
      if(!jet) continue;
      if(TMath::Abs(jet->Eta()) > fTrackEtaCut) continue; //track acceptance cut
      if(jet->Pt() > pTleading){
         idxLeadingJet = ij; 
         pTleading     = jet->Pt();
      }
   } 

   Int_t  njets = 0;
   static Double_t jphi[100];
   static Double_t jeta[100];
   static Double_t jRsquared[100];

   for(Int_t ij = 0; ij< listJet->GetEntries(); ij++){ //loop over bg kT jets

      jet = (AliAODJet*)(listJet->At(ij));
      if(!jet) continue;
      if(TMath::Abs(jet->Eta()) > fTrackEtaCut) continue; 

      if((jet->Pt() > fBgMaxJetPt || ij== idxLeadingJet) && njets<100){
         jphi[njets]      = RelativePhi(jet->Phi(),0.0); //-pi,pi
         jeta[njets]      = jet->Eta();
         jRsquared[njets] = fSafetyMargin*jet->EffectiveAreaCharged()/TMath::Pi(); //scale jet radius
         njets++;
      }
   }

   
   static Double_t nOutCone[10][4];
   static Double_t sumPtOutOfCone[10][4];

   for(Int_t ie=0; ie<fnEta; ie++){
      for(Int_t ip=0; ip<fnPhi; ip++){
         nOutCone[ip][ie]       = 0.0;     //initialize counter
         sumPtOutOfCone[ip][ie] = 0.0;
      }
   }

   Double_t rndphi, rndeta;
   Double_t rndphishift, rndetashift;
   Double_t dphi, deta; 
   Bool_t   bIsInCone;

   for(Int_t it=0; it<fnTrials; it++){

      rndphi = fRandom->Uniform(0, fPhiSize);
      rndeta = fRandom->Uniform(0, fEtaSize);
                              
      for(Int_t ip=0; ip<fnPhi; ip++){  //move radom position to each cell
         rndphishift = rndphi + ip*fPhiSize - TMath::Pi();
         for(Int_t ie=0; ie<fnEta; ie++){
            rndetashift = rndeta + ie*fEtaSize - fEtaSize;

            bIsInCone = 0; //tag if trial is in the jet cone
            for(Int_t ij=0; ij<njets; ij++){
               deta = jeta[ij] - rndetashift;
               dphi = RelativePhi(rndphishift,jphi[ij]);
               if((dphi*dphi + deta*deta) < jRsquared[ij]){
                  bIsInCone = 1;
                  break;
               }
            }
            if(!bIsInCone) nOutCone[ip][ie]++;
         }
      }
   }
  

   //Sum particle pT outside jets in cells
   Int_t npart = listPart->GetEntries();
   Int_t phicell,etacell; 
   AliVParticle *part = NULL; 
   for(Int_t ip=0; ip < npart; ip++){ 
         
      part = (AliVParticle*) listPart->At(ip);
      if(!part) continue;
      if( TMath::Abs(part->Eta()) > fEtaSize) continue; //skip part outside +-0.7 in eta
 
      bIsInCone = 0; //init
      for(Int_t ij=0; ij<njets; ij++){
         dphi = RelativePhi(jphi[ij], part->Phi());
         deta = jeta[ij] - part->Eta();
         if((dphi*dphi + deta*deta) < jRsquared[ij]){
            bIsInCone = 1;
            break;
         }
      } 
      if(!bIsInCone){
         phicell = TMath::Nint(TMath::Floor((RelativePhi(part->Phi(),0.0) + TMath::Pi())/fPhiSize));
         etacell = TMath::Nint(TMath::Floor((part->Eta()+fEtaSize)/fEtaSize));
         sumPtOutOfCone[phicell][etacell]+= part->Pt();
      }
   }

   static Double_t rhoInCells[20];
   Double_t  relativeArea;
   Int_t  nCells=0;
   Double_t bufferArea=0.0, bufferPt=0.0; //sum cells where A< fJetFreeAreaFrac
   for(Int_t ip=0; ip<fnPhi; ip++){
      for(Int_t ie=0; ie<fnEta; ie++){
         relativeArea = nOutCone[ip][ie]/fnTrials;
         //cout<<"RA=  "<< relativeArea<<"  BA="<<bufferArea<<"    BPT="<< bufferPt <<endl;

         bufferArea += relativeArea;
         bufferPt   += sumPtOutOfCone[ip][ie];
         if(bufferArea > fJetFreeAreaFrac){ 
            rhoInCells[nCells] = bufferPt/(bufferArea*fCellArea);

            if(mode==1) fhCellAreaToMedianGen->Fill(bufferArea*fCellArea);  
            else fhCellAreaToMedian->Fill(bufferArea*fCellArea);

            bufferArea = 0.0; 
            bufferPt   = 0.0;
            nCells++;
         }
      }
   }


   if(nCells>0){
     rhoMedian = TMath::Median(nCells, rhoInCells);
     if(mode==1){ //mc data 
        fhEntriesToMedianGen->Fill(nCells);
     }else{ //real data
        fhEntriesToMedian->Fill(nCells);
     } 
   } 

}

//__________________________________________________________________
void AliAnalysisTaskJetCorePP::EstimateBgCone(TList *listJet, TList* listPart,  Double_t &rhoPerpCone){

   rhoPerpCone = 0.0; //init

   if(!listJet) return;
   Int_t njet  = listJet->GetEntries();
   Int_t npart = listPart->GetEntries();
   Double_t pTleading  = 0.0;
   Double_t phiLeading = 1000.;
   Double_t etaLeading = 1000.;
   
   for(Int_t jsig=0; jsig < njet; jsig++){ 
      AliAODJet* signaljet = (AliAODJet*)(listJet->At(jsig));
      if(!signaljet) continue;
      if((signaljet->Eta()<fJetEtaMin) && (fJetEtaMax<signaljet->Eta())) continue; //acceptance cut
      if(signaljet->Pt() >= pTleading){ //replace leading and subleading jet
         pTleading  = signaljet->Pt();
         phiLeading = signaljet->Phi();
         etaLeading = signaljet->Eta();
      }
   }  
   if( TMath::Abs(etaLeading) >  fTrackEtaCut - fBgConeR){
      etaLeading  =  (etaLeading > 0) ? (fTrackEtaCut - fBgConeR) :  (-fTrackEtaCut + fBgConeR);
   }
  
   Double_t phileftcone  = phiLeading + TMath::Pi()/2; 
   Double_t phirightcone = phiLeading - TMath::Pi()/2; 
   Double_t dp, de;

   for(Int_t ip=0; ip < npart; ip++){ 
         
      AliVParticle *part = (AliVParticle*) listPart->At(ip);
      if(!part){
         continue;
      }
      
      dp = RelativePhi(phileftcone, part->Phi());
      de = etaLeading - part->Eta();
      if( sqrt(dp*dp + de*de)< fBgConeR) rhoPerpCone += part->Pt();

      dp = RelativePhi(phirightcone, part->Phi());
      if( sqrt(dp*dp + de*de)< fBgConeR) rhoPerpCone += part->Pt();
 
   }

   //normalize total pT by two times cone are 
   rhoPerpCone = rhoPerpCone/(2*TMath::Pi()*fBgConeR*fBgConeR);

}
//__________________________________________________________________

void AliAnalysisTaskJetCorePP::ReadTClonesArray(TString bname, TList *list){
   //Convert TClones array of jets to TList 

   if(!strlen(bname.Data())){
      AliError(Form("Jet branch %s not set.", bname.Data()));
      return;
   }

   TClonesArray *array=0x0;
   if(fUseExchContainer){ //take input from exchange containers
      if(fJetBranchNameKine.Length()>0 && bname.CompareTo(fJetBranchNameKine.Data())==0){
         array = dynamic_cast<TClonesArray*>(GetInputData(1)); //connect exchange container slot 1
      }
      if(fJetBranchNameBgKine.Length()>0 && bname.CompareTo(fJetBranchNameBgKine.Data())==0){
        array = dynamic_cast<TClonesArray*>(GetInputData(2)); //connect exchange container slot 2
      }
   }else{ //take input from AOD
      if(fAODOut&&!array){
         array = dynamic_cast<TClonesArray*>(fAODOut->FindListObject(bname.Data()));
      }
      if(fAODExtension&&!array){
         array = dynamic_cast<TClonesArray*>(fAODExtension->GetAOD()->FindListObject(bname.Data()));
      }
      if(fAODIn&&!array){
         array = dynamic_cast<TClonesArray*>(fAODIn->FindListObject(bname.Data()));
      }
   }

   list->Clear(); //clear list beforehand

   if(array){
      if(fDebug){
         Printf("########## %s: %d jets",bname.Data(), array->GetEntriesFast());
      }

      for(Int_t iJet = 0; iJet < array->GetEntriesFast(); iJet++) {
         AliAODJet *jet = dynamic_cast<AliAODJet*>((*array)[iJet]);
         if (jet) list->Add(jet);
      }
   }

   return;
}

//__________________________________________________________________
void AliAnalysisTaskJetCorePP::FillDeltaPt(TList *jetList, TList *trkList, Double_t rhoMedian, Double_t rhoCones){
   //Estimate magnitude of background fluctuations in givent event

   if(!jetList) return;
   Int_t njet  = jetList->GetEntries();
   Double_t pTleading  = 0.0;
   Double_t phiLeading = 1000.;
   Double_t etaLeading = 1000.;
   Double_t pTSubleading  = 0.0;
   Double_t phiSubLeading = 1000.;
   Double_t etaSubLeading = 1000.;
   
   for(Int_t jsig=0; jsig < njet; jsig++){ 
      AliAODJet* signaljet = (AliAODJet*)(jetList->At(jsig));
      if(!signaljet) continue;
      if((signaljet->Eta()<fJetEtaMin) && (fJetEtaMax<signaljet->Eta())) continue; //acceptance cut
      if(signaljet->Pt() >= pTleading){ // find leading + subleading jet

         pTSubleading  = pTleading;
         phiSubLeading = phiLeading;
         etaSubLeading = etaLeading;

         pTleading  = signaljet->Pt();
         phiLeading = signaljet->Phi();
         etaLeading = signaljet->Eta();
      }else if(signaljet->Pt() >= pTSubleading){
         pTSubleading  = signaljet->Pt();
         phiSubLeading = signaljet->Phi();
         etaSubLeading = signaljet->Eta();
      }
   }  

   //Generate random cone far away from leading jet+subleading jet, at least 2*R
   Double_t rndphi = TMath::Pi() *  fRandom->Uniform(-1,1); 
   Double_t rndeta = fRandom->Uniform(fJetEtaMin,fJetEtaMax);
   Int_t iter = 0; 
   
   while( (sqrt( pow(RelativePhi(rndphi, phiLeading),2)   + pow(etaLeading -rndeta,2)) <  2*fJetParamR) ||
          (sqrt( pow(RelativePhi(rndphi, phiSubLeading),2)+ pow(etaSubLeading -rndeta,2)) <  2*fJetParamR)){
      rndphi = TMath::Pi() *  fRandom->Uniform(-1,1); 
      rndeta = fRandom->Uniform(fJetEtaMin,fJetEtaMax); 
      iter++;

      if(iter>1e6) break;
   }

   //Sum pT of particles in radom cone
   Double_t sumPt = 0., dp, de;
   Int_t ntrk = trkList->GetEntries();
   for(Int_t itrk=0; itrk < ntrk; itrk++){ 
      AliVParticle *track = (AliVParticle*) trkList->At(itrk);
      if(!track){
         continue;
      }
      
      dp = RelativePhi( rndphi, track->Phi());
      de = rndeta - track->Eta();
      if( sqrt(dp*dp + de*de)< fJetParamR ) sumPt += track->Pt();
   }

   fhDeltaPtConeBg->Fill(sumPt - TMath::Pi()*fJetParamR * fJetParamR * rhoCones); 
   fhDeltaPtMedianBg->Fill(sumPt - TMath::Pi()*fJetParamR * fJetParamR * rhoMedian);
}
