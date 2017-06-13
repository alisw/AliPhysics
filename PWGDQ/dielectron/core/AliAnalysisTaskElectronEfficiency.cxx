/*************************************************************************
 * Copyright(c) 1998-2009, ALICE Experiment at CERN, All rights reserved. *
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

///////////////////////////////////////////////////////////////////////////
//                                                                       //
//       Single Electron, Pair and Pair-Prefilter Efficiency Task        //
//                                        (description in .h file)       //
//                                                                       //
///////////////////////////////////////////////////////////////////////////

#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "THnSparse.h"
#include "THStack.h"
#include "TGraph.h"
#include "TChain.h"
#include <math.h>
#include "TObjArray.h"
#include "TList.h"
#include <TLorentzVector.h>
#include "TVectorD.h"

#include "AliAnalysisManager.h"
#include "AliRunLoader.h"
#include "AliLoader.h"
#include "AliESDtrack.h"
#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliKFParticle.h"
#include "AliESDVertex.h"
#include "AliStack.h"
#include "AliVertex.h"
#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliPIDResponse.h"
#include "AliESDtrackCuts.h"
#include "AliCentrality.h"
#include "AliMultSelection.h" // Run2 centrality
#include "AliDielectron.h"
#include "AliDielectronVarManager.h"
#include "AliDielectronHistos.h"
#include "AliDielectronPID.h"
#include "AliDielectronVarCuts.h"
#include "AliDielectronTrackCuts.h"
#include "AliDielectronCutGroup.h"
#include "AliDielectronHelper.h"
#include "AliDielectronMC.h"
#include "AliDielectronSignalMC.h"
// todo: clean up includes...

#include "AliAnalysisTaskElectronEfficiency.h"
#include <bitset>

ClassImp(AliAnalysisTaskElectronEfficiency)

using std::cout;
using std::endl;


//________________________________________________________________________
AliAnalysisTaskElectronEfficiency::AliAnalysisTaskElectronEfficiency() :
AliAnalysisTaskSE(),
fESD(0x0),
mcEvent(0x0),
fPIDResponse(0x0),
fPostPIDCntrdCorrTPC(0x0),
fPostPIDWdthCorrTPC(0x0),
fPostPIDCntrdCorrITS(0x0),
fPostPIDWdthCorrITS(0x0),
fPostPIDCntrdCorrTOF(0x0),
fPostPIDWdthCorrTOF(0x0),
fUsedVars(0x0),
fSignalsMC(0x0),
fDoPairing(kFALSE),
fSelectPhysics(kFALSE),
fTriggerMask(AliVEvent::kAny),
fEventFilter(0x0),
fRequireVtx(kFALSE),
fNminEleInEventForRej(2),
fSupportedCutInstance(0),
//fEventcount(0),
fRandomizeDaughters(kFALSE),
fRandom3(4357),//default seed 4357
fMaxVtxZ(0.),
fCentMin( -1.),
fCentMax(101.),
fUseMultSelection(kFALSE),
fEtaMinGEN(-10.),
fEtaMaxGEN( 10.),
fPtMinGEN( -1.),
fPtMaxGEN(100.),
fNptBins(0),
fNetaBins(0),
fNphiBins(0),
fPtBins(0x0),
fEtaBins(0x0),
fPhiBins(0x0),
fsRunBins(""),
fCalcEfficiencyGen(kTRUE),
fCalcEfficiencyRec(kFALSE),
fCalcEfficiencyPoslabel(kTRUE),
fvTrackCuts(),
fvExtraTrackCuts(),
fvDoPrefilterEff(),
fvRejCutMee(),
fvRejCutTheta(),
fvRejCutPhiV(),
fNgen_Ele(0x0),
fvReco_Ele(),
fvReco_Ele_poslabel(),
fNgen_Pos(0x0),
fvReco_Pos(),
fvReco_Pos_poslabel(),
fNgen1_Rec_Ele(0x0),
fNgen2_Rec_Ele(0x0),
fvReco_Rec_Ele(),
fvReco_Rec_Ele_poslabel(),
fNgen1_Rec_Pos(0x0),
fNgen2_Rec_Pos(0x0),
fvReco_Rec_Pos(),
fvReco_Rec_Pos_poslabel(),
fNmeeBins(0),
fNpteeBins(0),
fMeeBins(0x0),
fPteeBins(0x0),
fNgenPairsResonances(0x0),
fNgenPairsDiffMothers(0x0),
fNgenPairsCharm(0x0),
fNgenPairsBeauty(0x0),
fNgenPairsHF(0x0),
fvRecoPairsResonances(),
fvRecoPairsDiffMothers(),
fvRecoPairsCharm(),
fvRecoPairsBeauty(),
fvRecoPairsHF(),
fNgenPairsRecResonances(0x0),
fNgenPairsRecDiffMothers(0x0),
fNgenPairsRecCharm(0x0),
fNgenPairsRecBeauty(0x0),
fNgenPairsRecHF(0x0),
fvRecoPairsRecResonances(),
fvRecoPairsRecDiffMothers(),
fvRecoPairsRecCharm(),
fvRecoPairsRecBeauty(),
fvRecoPairsRecHF(),
fCalcResolution(kFALSE),
fMakeResolutionSparse(kFALSE),
fTHnResElectrons1(0x0),
fTHnResPositrons1(0x0),
fTHnResElectrons2(0x0),
fTHnResPositrons2(0x0),
fDeltaPhiAll(0x0),
fDeltaPhi(0x0),
fDeltaPhi_alpha(0x0),
fDeltaPhi_pt(0x0),
fDeltaPhi_eta(0x0),
fDeltaPhi_MCcharge(0x0),
fDeltaPhi_charge(0x0),
fDeltaMomNbins(1200),
fDeltaMomMin(-10.),
fDeltaMomMax(2.),
fMomNbins(1000),
fMomMin(0.),
fMomMax(10.),
fRelMomNbins(400),
fRelMomMin(0.),
fRelMomMax(2.),
fDeltaEtaNbins(200),
fDeltaEtaMin(-0.4),
fDeltaEtaMax(0.4),
fDeltaThetaNbins(200),
fDeltaThetaMin(-0.4),
fDeltaThetaMax(0.4),
fDeltaPhiNbins(200),
fDeltaPhiMin(-0.4),
fDeltaPhiMax(0.4),
fDeltaAngleNbins(300),
fDeltaAngleMin(-0.5),
fDeltaAngleMax(1.),
fPGen(0x0),
fPRec(0x0),
fPGen_DeltaP(0x0),
fPtGen_DeltaPt(0x0),
fPGen_PrecOverPGen(0x0),
fPtGen_PtRecOverPtGen(0x0),
fPGen_DeltaEta(0x0),
fPGen_DeltaTheta(0x0),
fPGen_DeltaPhi_Ele(0x0),
fPGen_DeltaPhi_Pos(0x0),
fEtaGen_DeltaEta(0x0),
fThetaGen_DeltaTheta(0x0),
fPhiGen_DeltaPhi(0x0),
fOpeningAngleGen_DeltaOpeningAngleUS(0x0),
fOpeningAngleGen_DeltaOpeningAngleLS(0x0),
fMgen_PtGen_mRes_ptRes(0x0),
fPResArr(0x0),
fUseRelPResolution(kFALSE),
fThetaResArr(0x0),
fEtaResArr(0x0),
fPhiEleResArr(0x0),
fPhiPosResArr(0x0),
fResolutionCuts(0x0),
fPairCutMee(-1),
fPairCutTheta(-1),
fPairCutPhiV(3.2),
fvAllPionsForRej(),
fvPionsRejByAllSigns(),
fvPionsRejByUnlike(),
fOutputList(0x0),
fOutputListSupportHistos(0x0),
fEventStat(0x0),
tracksT(0x0),
fWriteTree(kFALSE),
pxESD(-1.), // tree variables
pyESD(-1.),
pzESD(-1.),
pTPC(-1.),
chargeT(-999),
signalITS(-1.),
signalTPC(-1.),
beta(-1.),
kNclsITS(-1),
kITSchi2Cl(-1.),
kNclsSITS(-1),
kNclsTPC(-1),
kTPCchi2Cl(-1.),
kNclsSTPC(-1),
kNclsTPCdEdx(-1),
kNFclsTPCr(-1.),
kNFclsTPCfCross(-1.),
kNtrkltsTRD(-1),
kNtrkltsTRDPID(-1),
sigmaEleITS(-999.),
sigmaEleTPC(-999.),
sigmaEleTOF(-999.),
probEleTRD(-999.),
sigmaPioITS(-999.),
sigmaPioTPC(-999.),
sigmaPioTOF(-999.),
sigmaKaoITS(-999.),
sigmaKaoTPC(-999.),
sigmaProITS(-999.),
sigmaProTPC(-999.),
isGlobalT(kFALSE),
isGlobalSDD(kFALSE),
labelT(-1),
pdgT(-1),
labelmotherT(-1),
pdgmotherT(-1),
labelgrandmotherT(-1),
pdggrandmotherT(-1),
pxMC(-1.),
pyMC(-1.),
pzMC(-1.),
fSelectedByCut(0),
fSelectedByExtraCut(0)
{
  /// Default Constructor
}


//________________________________________________________________________
AliAnalysisTaskElectronEfficiency::AliAnalysisTaskElectronEfficiency(const char *name) :
AliAnalysisTaskSE(name),
fESD(0x0),
mcEvent(0x0),
fPIDResponse(0x0),
fPostPIDCntrdCorrTPC(0x0),
fPostPIDWdthCorrTPC(0x0),
fPostPIDCntrdCorrITS(0x0),
fPostPIDWdthCorrITS(0x0),
fPostPIDCntrdCorrTOF(0x0),
fPostPIDWdthCorrTOF(0x0),
fUsedVars(0x0),
fSignalsMC(0x0),
fDoPairing(kFALSE),
fSelectPhysics(kFALSE),
fTriggerMask(AliVEvent::kAny),
fEventFilter(0x0),
fRequireVtx(kFALSE),
fNminEleInEventForRej(2),
fSupportedCutInstance(0),
//fEventcount(0),
fRandomizeDaughters(kFALSE),
fRandom3(4357),//default seed 4357
fMaxVtxZ(0.),
fCentMin( -1.),
fCentMax(101.),
fUseMultSelection(kFALSE),
fEtaMinGEN(-10.),
fEtaMaxGEN( 10.),
fPtMinGEN( -1.),
fPtMaxGEN(100.),
fNptBins(0),
fNetaBins(0),
fNphiBins(0),
fPtBins(0x0),
fEtaBins(0x0),
fPhiBins(0x0),
fsRunBins(""),
fCalcEfficiencyGen(kTRUE),
fCalcEfficiencyRec(kFALSE),
fCalcEfficiencyPoslabel(kTRUE),
fvTrackCuts(),
fvExtraTrackCuts(),
fvDoPrefilterEff(),
fvRejCutMee(),
fvRejCutTheta(),
fvRejCutPhiV(),
fNgen_Ele(0x0),
fvReco_Ele(),
fvReco_Ele_poslabel(),
fNgen_Pos(0x0),
fvReco_Pos(),
fvReco_Pos_poslabel(),
fNgen1_Rec_Ele(0x0),
fNgen2_Rec_Ele(0x0),
fvReco_Rec_Ele(),
fvReco_Rec_Ele_poslabel(),
fNgen1_Rec_Pos(0x0),
fNgen2_Rec_Pos(0x0),
fvReco_Rec_Pos(),
fvReco_Rec_Pos_poslabel(),
fNmeeBins(0),
fNpteeBins(0),
fMeeBins(0x0),
fPteeBins(0x0),
fNgenPairsResonances(0x0),
fNgenPairsDiffMothers(0x0),
fNgenPairsCharm(0x0),
fNgenPairsBeauty(0x0),
fNgenPairsHF(0x0),
fvRecoPairsResonances(),
fvRecoPairsDiffMothers(),
fvRecoPairsCharm(),
fvRecoPairsBeauty(),
fvRecoPairsHF(),
fNgenPairsRecResonances(0x0),
fNgenPairsRecDiffMothers(0x0),
fNgenPairsRecCharm(0x0),
fNgenPairsRecBeauty(0x0),
fNgenPairsRecHF(0x0),
fvRecoPairsRecResonances(),
fvRecoPairsRecDiffMothers(),
fvRecoPairsRecCharm(),
fvRecoPairsRecBeauty(),
fvRecoPairsRecHF(),
fCalcResolution(kFALSE),
fMakeResolutionSparse(kFALSE),
fTHnResElectrons1(0x0),
fTHnResPositrons1(0x0),
fTHnResElectrons2(0x0),
fTHnResPositrons2(0x0),
fDeltaPhiAll(0x0),
fDeltaPhi(0x0),
fDeltaPhi_alpha(0x0),
fDeltaPhi_pt(0x0),
fDeltaPhi_eta(0x0),
fDeltaPhi_MCcharge(0x0),
fDeltaPhi_charge(0x0),
fMomNbins(1000),
fMomMin(0.),
fMomMax(10.),
fDeltaMomNbins(1200),
fDeltaMomMin(-10.),
fDeltaMomMax(2.),
fRelMomNbins(400),
fRelMomMin(0.),
fRelMomMax(2.),
fDeltaEtaNbins(200),
fDeltaEtaMin(-0.4),
fDeltaEtaMax(0.4),
fDeltaThetaNbins(200),
fDeltaThetaMin(-0.4),
fDeltaThetaMax(0.4),
fDeltaPhiNbins(200),
fDeltaPhiMin(-0.4),
fDeltaPhiMax(0.4),
fDeltaAngleNbins(300),
fDeltaAngleMin(-0.5),
fDeltaAngleMax(1.),
fPGen(0x0),
fPRec(0x0),
fPGen_DeltaP(0x0),
fPtGen_DeltaPt(0x0),
fPGen_PrecOverPGen(0x0),
fPtGen_PtRecOverPtGen(0x0),
fPGen_DeltaEta(0x0),
fPGen_DeltaTheta(0x0),
fPGen_DeltaPhi_Ele(0x0),
fPGen_DeltaPhi_Pos(0x0),
fEtaGen_DeltaEta(0x0),
fThetaGen_DeltaTheta(0x0),
fPhiGen_DeltaPhi(0x0),
fOpeningAngleGen_DeltaOpeningAngleUS(0x0),
fOpeningAngleGen_DeltaOpeningAngleLS(0x0),
fMgen_PtGen_mRes_ptRes(0x0),
fPResArr(0x0),
fUseRelPResolution(kFALSE),
fThetaResArr(0x0),
fEtaResArr(0x0),
fPhiEleResArr(0x0),
fPhiPosResArr(0x0),
fResolutionCuts(0x0),
fPairCutMee(-1),
fPairCutTheta(-1),
fPairCutPhiV(3.2),
fvAllPionsForRej(),
fvPionsRejByAllSigns(),
fvPionsRejByUnlike(),
fOutputList(0x0),
fOutputListSupportHistos(0x0),
fEventStat(0x0),
tracksT(0x0),
fWriteTree(kFALSE),
pxESD(-1.), // tree variables
pyESD(-1.),
pzESD(-1.),
pTPC(-1.),
chargeT(-999),
signalITS(-1.),
signalTPC(-1.),
beta(-1.),
kNclsITS(-1),
kITSchi2Cl(-1.),
kNclsSITS(-1),
kNclsTPC(-1),
kTPCchi2Cl(-1.),
kNclsSTPC(-1),
kNclsTPCdEdx(-1),
kNFclsTPCr(-1.),
kNFclsTPCfCross(-1.),
kNtrkltsTRD(-1),
kNtrkltsTRDPID(-1),
sigmaEleITS(-999.),
sigmaEleTPC(-999.),
sigmaEleTOF(-999.),
probEleTRD(-999.),
sigmaPioITS(-999.),
sigmaPioTPC(-999.),
sigmaPioTOF(-999.),
sigmaKaoITS(-999.),
sigmaKaoTPC(-999.),
sigmaProITS(-999.),
sigmaProTPC(-999.),
isGlobalT(kFALSE),
isGlobalSDD(kFALSE),
labelT(-1),
pdgT(-1),
labelmotherT(-1),
pdgmotherT(-1),
labelgrandmotherT(-1),
pdggrandmotherT(-1),
pxMC(-1.),
pyMC(-1.),
pzMC(-1.),
fSelectedByCut(0),
fSelectedByExtraCut(0)
{
  /// Constructor
  
  fvTrackCuts.clear();
  fvExtraTrackCuts.clear();
  fvDoPrefilterEff.clear();
  fvRejCutMee.clear();
  fvRejCutTheta.clear();
  fvRejCutPhiV.clear();
  fvReco_Ele.clear();
  fvReco_Ele_poslabel.clear();
  fvReco_Pos.clear();
  fvReco_Pos_poslabel.clear();
  fvReco_Rec_Ele.clear();
  fvReco_Rec_Ele_poslabel.clear();
  fvReco_Rec_Pos.clear();
  fvReco_Rec_Pos_poslabel.clear();
  fvAllPionsForRej.clear();
  fvPionsRejByAllSigns.clear();
  fvPionsRejByUnlike.clear();
  
  
  DefineInput(0, TChain::Class());
  DefineOutput(1, TList::Class());
  DefineOutput(2, TList::Class());
  DefineOutput(3, TTree::Class());
  DefineOutput(4, TH1D::Class());
  
  fUsedVars = new TBits(AliDielectronVarManager::kNMaxValues);
  fUsedVars->SetBitNumber(AliDielectronVarManager::kP, kTRUE);
  fUsedVars->SetBitNumber(AliDielectronVarManager::kPIn, kTRUE);
  fUsedVars->SetBitNumber(AliDielectronVarManager::kITSnSigmaEle, kTRUE);
  fUsedVars->SetBitNumber(AliDielectronVarManager::kTPCnSigmaEle, kTRUE);
  // for ESD tracks, all variables are filled anyhow...
}


//________________________________________________________________________
void AliAnalysisTaskElectronEfficiency::UserCreateOutputObjects()
{
  /// Create histograms
  /// Called once
  AliInfo("Create the output objects. Do other one-time tasks.");
  Printf("Now running: CreateOutputObjects()");
  if(!fCalcEfficiencyGen && !fCalcEfficiencyRec){
    AliWarning("WARNING: set at least one of SetCalcEfficiencyGen() and SetCalcEfficiencyRec() to kTRUE in order to fill histograms.");
  }
  if(fCalcEfficiencyGen){
    Printf("Calculating efficiency in generated binning! ");
  }
  if(fCalcEfficiencyRec){
    Printf("Calculating efficiency %s in reconstructed binning! ", ((fCalcEfficiencyGen)?"ALSO":"ONLY"));
    std::cout << "  pResArr:      " << fPResArr ;
    if(fUseRelPResolution) std::cout << "   relative momentum resolution will be used" << std::endl;
    else std::cout << "   momentum difference will be used" << std::endl;
    std::cout << "  thetaResArr:  " << fThetaResArr << std::endl;
    std::cout << "  etaResArr:    " << fEtaResArr << std::endl;
    std::cout << "  phiEleResArr: " << fPhiEleResArr << std::endl;
    std::cout << "  phiPosResArr: " << fPhiPosResArr << std::endl;
  }
  /// Check if an MC signal was attached
  if(!fSignalsMC) {
    AliFatal("Task needs an AliDielectronSignalMC as basis for the electron selection!"
             "for example, define such function in your Config and call it from the AddTask:\n"
             "void SetupMCSignals(AliAnalysisTaskElectronEfficiency* task){\n"
             "  AliDielectronSignalMC* eleFinalState = new AliDielectronSignalMC(\"eleFinalState\",\"eleFinalState\");\n"
             "  eleFinalState->SetFillPureMCStep(kFALSE);\n"
             "  eleFinalState->SetLegPDGs(11,1);//dummy second leg (never MCtrue)\n"
             "  eleFinalState->SetCheckBothChargesLegs(kTRUE,kTRUE);\n"
             "  eleFinalState->SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);\n"
             "  eleFinalState->SetMotherSources(AliDielectronSignalMC::kDirect, AliDielectronSignalMC::kDirect);//equiv. to IsPrimary();\n"
             "  task->AddSignalMC(eleFinalState);\n}");
  }
  
  
  Printf("  __________________________________________________");
  Printf("  - centrality range:  %f  to  %f", fCentMin, fCentMax);
  Printf("  - use multiplicity selection (Run 2): %s", (fUseMultSelection)?"kTRUE":"kFALSE");
  Printf("  ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^");
  
  OpenFile(1, "RECREATE");
  
  TList *pairEffList(0x0);
  if(fDoPairing){
    pairEffList = new TList();
    pairEffList->SetName("pairEfficiency");
    pairEffList->SetOwner();
    
    if(fCalcEfficiencyGen){
      TList *lGenBinning = new TList();
      lGenBinning->SetName("generatedBinning");
      lGenBinning->SetOwner();
      TList *lReso1 = new TList();
      lReso1->SetName("resonances");
      lReso1->SetOwner();
      lReso1->Add(fNgenPairsResonances);
      for(UInt_t iCut=0; iCut<GetNCutsets(); ++iCut) lReso1->Add(fvRecoPairsResonances.at(iCut));
      lGenBinning->Add(lReso1);
      TList *lDiff1 = new TList();
      lDiff1->SetName("differentMothers");
      lDiff1->SetOwner();
      lDiff1->Add(fNgenPairsDiffMothers);
      for(UInt_t iCut=0; iCut<GetNCutsets(); ++iCut) lDiff1->Add(fvRecoPairsDiffMothers.at(iCut));
      lGenBinning->Add(lDiff1);
      TList *lCharm1 = new TList();
      lCharm1->SetName("charm");
      lCharm1->SetOwner();
      lCharm1->Add(fNgenPairsCharm);
      for(UInt_t iCut=0; iCut<GetNCutsets(); ++iCut) lCharm1->Add(fvRecoPairsCharm.at(iCut));
      lGenBinning->Add(lCharm1);
      TList *lBeauty1 = new TList();
      lBeauty1->SetName("beauty");
      lBeauty1->SetOwner();
      lBeauty1->Add(fNgenPairsBeauty);
      for(UInt_t iCut=0; iCut<GetNCutsets(); ++iCut) lBeauty1->Add(fvRecoPairsBeauty.at(iCut));
      lGenBinning->Add(lBeauty1);
      TList *lHF1 = new TList();
      lHF1->SetName("heavyflavor");
      lHF1->SetOwner();
      lHF1->Add(fNgenPairsHF);
      for(UInt_t iCut=0; iCut<GetNCutsets(); ++iCut) lHF1->Add(fvRecoPairsHF.at(iCut));
      lGenBinning->Add(lHF1);
      
      //lGenBinning->Print();
      pairEffList->Add(lGenBinning);
    }
    
    if(fCalcEfficiencyRec){
      TList *lRecBinning = new TList();
      lRecBinning->SetName("reconstructedBinning");
      lRecBinning->SetOwner();
      TList *lReso = new TList();
      lReso->SetName("resonances");
      lReso->SetOwner();
      lReso->Add(fNgenPairsRecResonances);
      for(UInt_t iCut=0; iCut<GetNCutsets(); ++iCut) lReso->Add(fvRecoPairsRecResonances.at(iCut));
      lRecBinning->Add(lReso);
      TList *lDiff = new TList();
      lDiff->SetName("differentMothers");
      lDiff->SetOwner();
      lDiff->Add(fNgenPairsRecDiffMothers);
      for(UInt_t iCut=0; iCut<GetNCutsets(); ++iCut) lDiff->Add(fvRecoPairsRecDiffMothers.at(iCut));
      lRecBinning->Add(lDiff);
      TList *lCharm = new TList();
      lCharm->SetName("charm");
      lCharm->SetOwner();
      lCharm->Add(fNgenPairsRecCharm);
      for(UInt_t iCut=0; iCut<GetNCutsets(); ++iCut) lCharm->Add(fvRecoPairsRecCharm.at(iCut));
      lRecBinning->Add(lCharm);
      TList *lBeauty = new TList();
      lBeauty->SetName("beauty");
      lBeauty->SetOwner();
      lBeauty->Add(fNgenPairsRecBeauty);
      for(UInt_t iCut=0; iCut<GetNCutsets(); ++iCut) lBeauty->Add(fvRecoPairsRecBeauty.at(iCut));
      lRecBinning->Add(lBeauty);
      TList *lHF = new TList();
      lHF->SetName("heavyflavor");
      lHF->SetOwner();
      lHF->Add(fNgenPairsRecHF);
      for(UInt_t iCut=0; iCut<GetNCutsets(); ++iCut) lHF->Add(fvRecoPairsRecHF.at(iCut));
      lRecBinning->Add(lHF);
      
      //lRecBinning->Print();
      pairEffList->Add(lRecBinning);
    }
    pairEffList->Print();
  }
  
  TList *resolutionList(0x0);
  if(fCalcResolution){
    resolutionList = new TList();
    resolutionList->SetName("resolution");
    resolutionList->SetOwner();
    
    if(fMakeResolutionSparse){
      Int_t THnBins[6] = {fMomNbins, fRelMomNbins, fDeltaMomNbins, fDeltaThetaNbins, fDeltaEtaNbins, fDeltaPhiNbins};
      Double_t THnMin[6] = {fMomMin, fRelMomMin, fDeltaMomMin, fDeltaThetaMin, fDeltaEtaMin, fDeltaPhiMin};
      Double_t THnMax[6] = {fMomMax, fRelMomMax, fDeltaMomMax, fDeltaThetaMax, fDeltaEtaMax, fDeltaPhiMax};
      fTHnResElectrons1 = new THnSparseD("pGen_Res_Electrons1", "pGen_Res_Electrons1", 6, THnBins, THnMin, THnMax);
      fTHnResElectrons1->Sumw2();
      
      fTHnResElectrons1->GetAxis(0)->SetName("pGen");
      fTHnResElectrons1->GetAxis(1)->SetName("pGen_Over_pRec");
      fTHnResElectrons1->GetAxis(2)->SetName("deltaP");
      fTHnResElectrons1->GetAxis(3)->SetName("deltaTheta");
      fTHnResElectrons1->GetAxis(4)->SetName("deltaEta");
      fTHnResElectrons1->GetAxis(5)->SetName("deltaPhi");
      
      fTHnResElectrons1->GetAxis(0)->SetTitle("p^{gen} (GeV/c)");
      fTHnResElectrons1->GetAxis(1)->SetTitle("p^{rec} / p^{gen} (GeV/c)");
      fTHnResElectrons1->GetAxis(2)->SetTitle("p^{rec}_{T} - p^{gen}_{T} (GeV/c)");
      fTHnResElectrons1->GetAxis(3)->SetTitle("#theta^{rec} - #theta^{gen} (rad)");
      fTHnResElectrons1->GetAxis(4)->SetTitle("#eta^{rec} - #eta^{gen}");
      fTHnResElectrons1->GetAxis(5)->SetTitle("#varphi^{rec} - #varphi^{gen} (rad)");
      
      fTHnResPositrons1 = static_cast<THnSparseD*> (fTHnResElectrons1->Clone("pGen_Res_Positrons1"));
      fTHnResPositrons1->SetTitle("pGen_Res_Positrons1");
      fTHnResPositrons1->Sumw2();
      
      fTHnResElectrons2 = static_cast<THnSparseD*> (fTHnResElectrons1->Clone("pGen_Res_Electrons_negativeLabel"));
      fTHnResElectrons2->SetTitle("pGen_Res_Electrons_negativeLabel");
      fTHnResElectrons2->Sumw2();
      
      fTHnResPositrons2 = static_cast<THnSparseD*> (fTHnResElectrons1->Clone("pGen_Res_Positrons_negativeLabel"));
      fTHnResPositrons2->SetTitle("pGen_Res_Positrons_negativeLabel");
      fTHnResPositrons2->Sumw2();
            
      resolutionList->Add(fTHnResElectrons1);
      resolutionList->Add(fTHnResPositrons1);
      resolutionList->Add(fTHnResElectrons2);
      resolutionList->Add(fTHnResPositrons2);
    }
    else{
      fDeltaPhiAll = new TH1D("DeltaPhiAll","",320,-0.1,6.4);
      fDeltaPhiAll->GetXaxis()->SetTitle("#varphi_{gen} - #varphi_{rec}");
      fDeltaPhiAll->Sumw2();
      fDeltaPhiAll->SetMarkerStyle(20);
      fDeltaPhiAll->SetMarkerColor(kBlack);
      fDeltaPhiAll->SetLineColor(kBlack);
      
      fDeltaPhi = new TH1D("DeltaPhi","",320,-0.1,6.4);
      fDeltaPhi->GetXaxis()->SetTitle("#varphi_{gen} - #varphi_{rec}");
      fDeltaPhi->Sumw2();
      fDeltaPhi->SetMarkerStyle(20);
      fDeltaPhi->SetMarkerColor(kBlack);
      fDeltaPhi->SetLineColor(kBlack);
      
      fDeltaPhi_alpha = new TH2D("DeltaPhi_alpha","",320,-0.1,6.4,140,-7.,7.);
      fDeltaPhi_alpha->Sumw2();
      fDeltaPhi_pt = new TH2D("DeltaPhi_pt","",320,-0.1,6.4,500,0.,10.);
      fDeltaPhi_pt->Sumw2();
      fDeltaPhi_eta = new TH2D("DeltaPhi_eta","",320,-0.1,6.4,200,-2.,2.);
      fDeltaPhi_eta->Sumw2();
      fDeltaPhi_MCcharge = new TH2D("DeltaPhi_MCcharge","",320,-0.1,6.4,3,-1.5,1.5);
      fDeltaPhi_MCcharge->Sumw2();
      fDeltaPhi_charge = new TH2D("DeltaPhi_charge","",320,-0.1,6.4,3,-1.5,1.5);
      fDeltaPhi_charge->Sumw2();
      
      fPGen                                = new TH1D("PGen",                               "",500,0., 5.);
      fPRec                                = new TH1D("PRec",                               "",500,0., 5.);
      fPGen_DeltaP                         = new TH2D("PGen_DeltaP",                        "",500,0.,10.,fDeltaMomNbins,fDeltaMomMin,fDeltaMomMax);
      fPtGen_DeltaPt                       = new TH2D("PtGen_DeltaPt",                      "",500,0.,10.,fDeltaMomNbins,fDeltaMomMin,fDeltaMomMax);
      fPGen_PrecOverPGen                   = new TH2D("PGen_PrecOverPGen",                  "",500,0.,10.,fRelMomNbins,fRelMomMin,fRelMomMax);
      fPtGen_PtRecOverPtGen                = new TH2D("PtGen_PtRecOverPtGen",               "",500,0.,10.,fRelMomNbins,fRelMomMin,fRelMomMax);
      fPGen_DeltaEta                       = new TH2D("PGen_DeltaEta",                      "",500,0.,10.,fDeltaEtaNbins,fDeltaEtaMin,fDeltaEtaMax);
      fPGen_DeltaTheta                     = new TH2D("PGen_DeltaTheta",                    "",500,0.,10.,fDeltaThetaNbins,fDeltaThetaMin,fDeltaThetaMax);
      fPGen_DeltaPhi_Ele                   = new TH2D("PGen_DeltaPhi_Ele",                  "",500,0.,10.,fDeltaPhiNbins,fDeltaPhiMin,fDeltaPhiMax);
      fPGen_DeltaPhi_Pos                   = new TH2D("PGen_DeltaPhi_Pos",                  "",500,0.,10.,fDeltaPhiNbins,fDeltaPhiMin,fDeltaPhiMax);
      fEtaGen_DeltaEta                     = new TH2D("EtaGen_DeltaEta",                    "",200,-1.,1.,fDeltaEtaNbins,fDeltaEtaMin,fDeltaEtaMax);
      fThetaGen_DeltaTheta                 = new TH2D("ThetaGen_DeltaTheta",                "",220,-0.1*TMath::Pi(),1.1*TMath::Pi(),fDeltaThetaNbins,fDeltaThetaMin,fDeltaThetaMax);
      fPhiGen_DeltaPhi                     = new TH2D("PhiGen_DeltaPhi",                    "",320,-0.1*TMath::Pi(),2.1*TMath::Pi(),fDeltaPhiNbins,fDeltaPhiMin,fDeltaPhiMax);
      
      fPGen                                ->Sumw2();
      fPRec                                ->Sumw2();
      fPGen_DeltaP                         ->Sumw2();
      fPtGen_DeltaPt                       ->Sumw2();
      fPGen_PrecOverPGen                   ->Sumw2();
      fPtGen_PtRecOverPtGen                ->Sumw2();
      fPGen_DeltaEta                       ->Sumw2();
      fPGen_DeltaTheta                     ->Sumw2();
      fPGen_DeltaPhi_Ele                   ->Sumw2();
      fPGen_DeltaPhi_Pos                   ->Sumw2();
      fEtaGen_DeltaEta                     ->Sumw2();
      fThetaGen_DeltaTheta                 ->Sumw2();
      fPhiGen_DeltaPhi                     ->Sumw2();
      
      fPGen                                ->GetXaxis()->SetTitle("p^{gen} (GeV/c)");
      fPRec                                ->GetXaxis()->SetTitle("p^{rec} (GeV/c)");
      fPGen_DeltaP                         ->GetXaxis()->SetTitle("p^{gen} (GeV/c)");
      fPGen_DeltaP                         ->GetYaxis()->SetTitle("p^{rec} - p^{gen} (GeV/c)");
      fPtGen_DeltaPt                       ->GetXaxis()->SetTitle("p^{gen}_{T} (GeV/c)");
      fPtGen_DeltaPt                       ->GetYaxis()->SetTitle("p^{rec}_{T} - p^{gen}_{T} (GeV/c)");
      fPGen_PrecOverPGen                   ->GetXaxis()->SetTitle("p^{gen} (GeV/c)");
      fPGen_PrecOverPGen                   ->GetYaxis()->SetTitle("p^{rec} / p^{gen} (GeV/c)");
      fPtGen_PtRecOverPtGen                ->GetXaxis()->SetTitle("p^{gen}_{T} (GeV/c)");
      fPtGen_PtRecOverPtGen                ->GetYaxis()->SetTitle("p^{rec}_{T} / p^{gen}_{T} (GeV/c)");
      fPGen_DeltaEta                       ->GetXaxis()->SetTitle("p^{gen} (GeV/c)");
      fPGen_DeltaEta                       ->GetYaxis()->SetTitle("#eta^{rec} - #eta^{gen}");
      fPGen_DeltaTheta                     ->GetXaxis()->SetTitle("p^{gen} (GeV/c)");
      fPGen_DeltaTheta                     ->GetYaxis()->SetTitle("#theta^{rec} - #theta^{gen} (rad)");
      fPGen_DeltaPhi_Ele                   ->GetXaxis()->SetTitle("p^{gen} (GeV/c)");
      fPGen_DeltaPhi_Ele                   ->GetYaxis()->SetTitle("#varphi^{rec} - #varphi^{gen} (rad)");
      fPGen_DeltaPhi_Pos                   ->GetXaxis()->SetTitle("p^{gen} (GeV/c)");
      fPGen_DeltaPhi_Pos                   ->GetYaxis()->SetTitle("#varphi^{rec} - #varphi^{gen} (rad)");
      fEtaGen_DeltaEta                     ->GetXaxis()->SetTitle("#eta^{gen}");
      fEtaGen_DeltaEta                     ->GetYaxis()->SetTitle("#eta^{rec} - #eta^{gen}");
      fThetaGen_DeltaTheta                 ->GetXaxis()->SetTitle("#theta^{gen} (rad)");
      fThetaGen_DeltaTheta                 ->GetYaxis()->SetTitle("#theta^{rec} - #theta^{gen} (rad)");
      fPhiGen_DeltaPhi                     ->GetXaxis()->SetTitle("#varphi^{gen} (rad)");
      fPhiGen_DeltaPhi                     ->GetYaxis()->SetTitle("#varphi^{rec} - #varphi^{gen} (rad)");
      
      resolutionList->Add(fDeltaPhiAll);
      resolutionList->Add(fDeltaPhi_alpha);
      resolutionList->Add(fDeltaPhi_pt);
      resolutionList->Add(fDeltaPhi_eta);
      resolutionList->Add(fDeltaPhi_MCcharge);
      resolutionList->Add(fDeltaPhi_charge);
      resolutionList->Add(fDeltaPhi);
      
      resolutionList->Add(fPGen);
      resolutionList->Add(fPRec);
      resolutionList->Add(fPGen_DeltaP);
      resolutionList->Add(fPtGen_DeltaPt);
      resolutionList->Add(fPGen_PrecOverPGen   );
      resolutionList->Add(fPtGen_PtRecOverPtGen);
      resolutionList->Add(fPGen_DeltaEta       );
      resolutionList->Add(fPGen_DeltaTheta     );
      resolutionList->Add(fPGen_DeltaPhi_Ele   );
      resolutionList->Add(fPGen_DeltaPhi_Pos   );
      resolutionList->Add(fEtaGen_DeltaEta);
      resolutionList->Add(fThetaGen_DeltaTheta);
      resolutionList->Add(fPhiGen_DeltaPhi);
    }
    fOpeningAngleGen_DeltaOpeningAngleUS = new TH2D("OpeningAngleGen_DeltaOpeningAngleUS","",330,-0.1,3.2,fDeltaAngleNbins,fDeltaAngleMin,fDeltaAngleMax);
    fOpeningAngleGen_DeltaOpeningAngleLS = new TH2D("OpeningAngleGen_DeltaOpeningAngleLS","",330,-0.1,3.2,fDeltaAngleNbins,fDeltaAngleMin,fDeltaAngleMax);
    fOpeningAngleGen_DeltaOpeningAngleUS ->Sumw2();
    fOpeningAngleGen_DeltaOpeningAngleLS ->Sumw2();
    fOpeningAngleGen_DeltaOpeningAngleUS ->GetXaxis()->SetTitle("#theta_{ee,US}^{gen} (rad)");
    fOpeningAngleGen_DeltaOpeningAngleUS ->GetYaxis()->SetTitle("#theta_{ee,US}^{rec} - #theta_{ee,US}^{gen} (rad)");
    fOpeningAngleGen_DeltaOpeningAngleLS ->GetXaxis()->SetTitle("#theta_{ee,LS}^{gen} (rad)");
    fOpeningAngleGen_DeltaOpeningAngleLS ->GetYaxis()->SetTitle("#theta_{ee,LS}^{rec} - #theta_{ee,LS}^{gen} (rad)");
        
    resolutionList->Add(fOpeningAngleGen_DeltaOpeningAngleUS);
    resolutionList->Add(fOpeningAngleGen_DeltaOpeningAngleLS);

    
    Int_t bins[4] = {50,60,100,100};
    Double_t min[4] = { 0., 0., 0., 0. };
    Double_t max[4] = { 5., 6., 2., 2. };
    fMgen_PtGen_mRes_ptRes          = new THnSparseF("Mgen_PtGen_mRes_ptRes","",4,bins,min,max);
    fMgen_PtGen_mRes_ptRes->Sumw2();
    fMgen_PtGen_mRes_ptRes->GetAxis(0)->SetTitle("m_{ee}^{gen}");
    fMgen_PtGen_mRes_ptRes->GetAxis(1)->SetTitle("p_{T,ee}^{gen}");
    fMgen_PtGen_mRes_ptRes->GetAxis(2)->SetTitle("m_{ee}^{rec} / m_{ee}^{gen}");
    fMgen_PtGen_mRes_ptRes->GetAxis(3)->SetTitle("p_{T,ee}^{rec} / p_{T,ee}^{gen}");
    
    fMgen_PtGen_mRes_ptRes->GetAxis(0)->SetName("m_gen");
    fMgen_PtGen_mRes_ptRes->GetAxis(1)->SetName("pt_gen");
    fMgen_PtGen_mRes_ptRes->GetAxis(2)->SetName("mResolution");
    fMgen_PtGen_mRes_ptRes->GetAxis(3)->SetName("ptResolution");
    resolutionList->Add(fMgen_PtGen_mRes_ptRes);
  }
  
  TList *singleEffList = new TList();
  singleEffList->SetName("electronEfficiency");
  singleEffList->SetOwner();
  
  if(fCalcEfficiencyGen){
    TList *singleEffGenList = new TList();
    singleEffGenList->SetName("generatedBinning");
    singleEffGenList->SetOwner();
    
    singleEffGenList->Add(fNgen_Ele);
    for (UInt_t iCut=0; iCut<GetNCutsets(); ++iCut)
      singleEffGenList->Add(fvReco_Ele.at(iCut));
    if(fCalcEfficiencyPoslabel)
      for (UInt_t iCut=0; iCut<GetNCutsets(); ++iCut)
        singleEffGenList->Add(fvReco_Ele_poslabel.at(iCut));
    
    singleEffGenList->Add(fNgen_Pos);
    for (UInt_t iCut=0; iCut<GetNCutsets(); ++iCut)
      singleEffGenList->Add(fvReco_Pos.at(iCut));
    if(fCalcEfficiencyPoslabel)
      for (UInt_t iCut=0; iCut<GetNCutsets(); ++iCut)
        singleEffGenList->Add(fvReco_Pos_poslabel.at(iCut));
    
    //singleEffGenList->Print();
    singleEffList->Add(singleEffGenList);
  }
  
  if(fCalcEfficiencyRec){
    TList *singleEffRecList = new TList();
    singleEffRecList->SetName("reconstructedBinning");
    singleEffRecList->SetOwner();
    
    singleEffRecList->Add(fNgen1_Rec_Ele);
    singleEffRecList->Add(fNgen2_Rec_Ele);
    for (UInt_t iCut=0; iCut<GetNCutsets(); ++iCut)
      singleEffRecList->Add(fvReco_Rec_Ele.at(iCut));
    if(fCalcEfficiencyPoslabel)
      for (UInt_t iCut=0; iCut<GetNCutsets(); ++iCut)
        singleEffRecList->Add(fvReco_Rec_Ele_poslabel.at(iCut));
    
    singleEffRecList->Add(fNgen1_Rec_Pos);
    singleEffRecList->Add(fNgen2_Rec_Pos);
    for (UInt_t iCut=0; iCut<GetNCutsets(); ++iCut)
      singleEffRecList->Add(fvReco_Rec_Pos.at(iCut));
    if(fCalcEfficiencyPoslabel)
      for (UInt_t iCut=0; iCut<GetNCutsets(); ++iCut)
        singleEffRecList->Add(fvReco_Rec_Pos_poslabel.at(iCut));
    
    //singleEffRecList->Print();
    singleEffList->Add(singleEffRecList);
  }
  
  for (UInt_t iCut=0; iCut<GetNCutsets(); ++iCut){
    if(fvDoPrefilterEff.at(iCut)){
      singleEffList->Add(fvAllPionsForRej.at(iCut));
      singleEffList->Add(fvPionsRejByAllSigns.at(iCut));
      singleEffList->Add(fvPionsRejByUnlike.at(iCut));
    }
  }
  singleEffList->Print();
  
  
  // be really careful if you need to implement this (see comments in UserExec):
  //    fOutputList->Add(fvReco_Pio.at(iCut));
  //    fOutputList->Add(fvReco_Kao.at(iCut));
  //    fOutputList->Add(fvReco_Pro.at(iCut));
  
  fOutputList = new TList();
  fOutputList->SetName(GetName());
  fOutputList->SetOwner();
  
  if(singleEffList)  fOutputList->Add(singleEffList);
  if(pairEffList)    fOutputList->Add(pairEffList);
  if(resolutionList) fOutputList->Add(resolutionList);
  
  fOutputList->Print();
  
  PostData(1, fOutputList);
  
  OpenFile(2, "RECREATE");
  CreateSupportHistos();
  PostData(2, fOutputListSupportHistos);
  
  OpenFile(3, "RECREATE");
  //Define Tree for Track Values
  tracksT = new TTree("tracksT","tracksT");
  // branches for both data and MC:
  tracksT->Branch("pxESD",              &pxESD);
  tracksT->Branch("pyESD",              &pyESD);
  tracksT->Branch("pzESD",              &pzESD);
  tracksT->Branch("pTPC",               &pTPC);
  tracksT->Branch("chargeT",            &chargeT);
  tracksT->Branch("signalITS",          &signalITS);
  tracksT->Branch("signalTPC",          &signalTPC);
  tracksT->Branch("beta",               &beta);
  tracksT->Branch("kNclsITS",           &kNclsITS);
  tracksT->Branch("kITSchi2Cl",         &kITSchi2Cl);
  tracksT->Branch("kNclsSITS",          &kNclsSITS);
  tracksT->Branch("kNclsTPC",           &kNclsTPC);
  tracksT->Branch("kTPCchi2Cl",         &kTPCchi2Cl);
  tracksT->Branch("kNclsSTPC",          &kNclsSTPC);
  tracksT->Branch("kNclsTPCdEdx",       &kNclsTPCdEdx);
  tracksT->Branch("kNFclsTPCr",         &kNFclsTPCr);
  tracksT->Branch("kNFclsTPCfCross",    &kNFclsTPCfCross);
  tracksT->Branch("kNtrkltsTRD",        &kNtrkltsTRD);
  tracksT->Branch("kNtrkltsTRDPID",     &kNtrkltsTRDPID);
  tracksT->Branch("sigmaEleITS",        &sigmaEleITS);
  tracksT->Branch("sigmaEleTPC",        &sigmaEleTPC);
  tracksT->Branch("sigmaEleTOF",        &sigmaEleTOF);
  //tracksT->Branch("sigmaEleTRD",        &sigmaEleTRD);
  tracksT->Branch("probEleTRD",         &probEleTRD);
  tracksT->Branch("sigmaPioITS",        &sigmaPioITS);
  tracksT->Branch("sigmaPioTPC",        &sigmaPioTPC);
  tracksT->Branch("sigmaPioTOF",        &sigmaPioTOF);
  tracksT->Branch("sigmaKaoITS",        &sigmaKaoITS);
  tracksT->Branch("sigmaKaoTPC",        &sigmaKaoTPC);
  tracksT->Branch("sigmaProITS",        &sigmaProITS);
  tracksT->Branch("sigmaProTPC",        &sigmaProTPC);
  tracksT->Branch("isGlobalT",          &isGlobalT);
  tracksT->Branch("isGlobalSDD",        &isGlobalSDD);
  // additional branches for MC: (not possible to set them inside the "if (fIsMC) {...}")
  tracksT->Branch("labelT",             &labelT);
  tracksT->Branch("pdgT",               &pdgT);
  tracksT->Branch("labelmotherT",       &labelmotherT);
  tracksT->Branch("pdgmotherT",         &pdgmotherT);
  tracksT->Branch("labelgrandmotherT",  &labelgrandmotherT);
  tracksT->Branch("pdggrandmotherT",    &pdggrandmotherT);
  //tracksT->Branch("pMC",                &pMC);
  tracksT->Branch("pxMC",               &pxMC);
  tracksT->Branch("pyMC",               &pyMC);
  tracksT->Branch("pzMC",               &pzMC);
  tracksT->Branch("selectedByCut",      &fSelectedByCut);
  tracksT->Branch("selectedByExtraCut", &fSelectedByExtraCut);
  PostData(3, tracksT);
  
  OpenFile(4, "RECREATE");
  if (!fEventStat){
    fEventStat=new TH1D("hEventStat","Event statistics",kEventStatBins,0,kEventStatBins);
    fEventStat->GetXaxis()->SetBinLabel(1,"Before Phys. Sel.");
    fEventStat->GetXaxis()->SetBinLabel(2,"After Phys. Sel. (data only)");
    fEventStat->GetXaxis()->SetBinLabel(3,"After event filter");
  }
  PostData(4, fEventStat);
}


//________________________________________________________________________
void AliAnalysisTaskElectronEfficiency::UserExec(Option_t *)
{
  
  /// Strategy:  Process one cutInstance (or multiple ones) for analysis tracking&PID efficiency (as usual).
  ///            Process optional, separate cutInstance for prefilter efficiencies: it also produces the usual tracking&PID efficiency
  ///            (but of course for the specified prefilter track sample, so mainly for convenience and curiosity),
  ///            and then computes the pair rejection efficiency, using random rejection of "testparticles" with the selected electrons. (further info in 'CalcPrefilterEff()')
  if(!AliDielectronMC::Instance()->ConnectMCEvent()) return;
  mcEvent = AliDielectronMC::Instance()->GetMCEvent();
  if (!mcEvent) { Printf("ERROR: mcEvent not available"); return; }
  AliDielectronVarManager::SetEvent(mcEvent);
  
  AliESDInputHandler *inputHandlerESD = dynamic_cast<AliESDInputHandler*> (AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
  if (!inputHandlerESD) { Printf("ERROR: Could not get ESDInputHandler\n"); }
  else fESD = (AliESDEvent*)inputHandlerESD->GetEvent();
  if (!fESD) { Printf("ERROR: fESD not available"); return; }
  
  if (!fPIDResponse) SetPIDResponse( ((AliESDInputHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->GetPIDResponse() );
  AliDielectronVarManager::SetPIDResponse(fPIDResponse);
  
  // set pid correction function to var manager
  if(fPostPIDCntrdCorrTPC) AliDielectronPID::SetCentroidCorrFunction(fPostPIDCntrdCorrTPC);
  if(fPostPIDWdthCorrTPC)  AliDielectronPID::SetWidthCorrFunction(fPostPIDWdthCorrTPC);
  if(fPostPIDCntrdCorrITS) AliDielectronPID::SetCentroidCorrFunctionITS(fPostPIDCntrdCorrITS);
  if(fPostPIDWdthCorrITS)  AliDielectronPID::SetWidthCorrFunctionITS(fPostPIDWdthCorrITS);
  if(fPostPIDCntrdCorrTOF) AliDielectronPID::SetCentroidCorrFunctionTOF(fPostPIDCntrdCorrTOF);
  if(fPostPIDWdthCorrTOF)  AliDielectronPID::SetWidthCorrFunctionTOF(fPostPIDWdthCorrTOF);
  
  AliStack *fStack = mcEvent->Stack();
  
  Bool_t isESD=kTRUE;
  Bool_t isAOD=kFALSE; //currently not supported!
  // taken from AliAnalysisTaskMultiDielectron::UserExec():
  // -----
  // Was event selected ?
  ULong64_t isSelected = AliVEvent::kAny;
  //Bool_t isRejected = kFALSE;
  if( fSelectPhysics && inputHandlerESD){
    if((isESD && inputHandlerESD->GetEventSelection()) || isAOD){
      isSelected = inputHandlerESD->IsEventSelected();
      //if (fExcludeTriggerMask && (isSelected&fExcludeTriggerMask)) isRejected=kTRUE;
      //if (fTriggerLogic==kAny) isSelected&=fTriggerMask;
      //else if (fTriggerLogic==kExact) isSelected=((isSelected&fTriggerMask)==fTriggerMask);
      //cout << " isSelected = " << (bitset<32>)isSelected << " \t fTriggerMask = " << (bitset<32>)fTriggerMask << endl;
      isSelected&=fTriggerMask;
      //TString firedTriggerClasses=InputEvent()->GetFiredTriggerClasses();
      //if(!fFiredTrigger.IsNull()) isSelected=(firedTriggerClasses.Contains(fFiredTrigger))^fFiredExclude;
    }
  }
  // -----
  // Before physics selection
  fEventStat->Fill(kAllEvents);
  if (isSelected==0) {
    PostData(4,fEventStat);
    return;
  }
  // After physics selection
  fEventStat->Fill(kPhysicsSelectionEvents);
  
  // Event filter
  if (fEventFilter) {
    if (!fEventFilter->IsSelected(fESD)) return; // fESD instead of 'InputEvent()'
  }
  // centrality
  Double_t centralityF=-1;
  if (fUseMultSelection) { // Run 2
    AliMultSelection *multSelection = (AliMultSelection*)fESD->FindListObject("MultSelection");
    if (multSelection) centralityF  = multSelection->GetMultiplicityPercentile("V0M",kFALSE);
  }
  else {
    AliCentrality *esdCentrality    = fESD->GetCentrality();
    if (esdCentrality) centralityF  = esdCentrality->GetCentralityPercentile("V0M");
  }
  if (centralityF<fCentMin || centralityF>=fCentMax) return;
  // vertex, just for monitoring
  const AliESDVertex* vtxESD = fESD->GetPrimaryVertexTracks();
  Double_t vtxZGlobal = -99.;
  Int_t nCtrb = -1;
  if (vtxESD) {
    vtxZGlobal = vtxESD->GetZ(); // was GetZv(); until Jan 2015
    nCtrb = vtxESD->GetNContributors();
  }
  
  // After event filter
  fEventStat->Fill(kFilteredEvents);
  PostData(4,fEventStat);
  
  
  //  ++fEventcount;
  //  Printf("__________ next Event selected ( %i ) __________", fEventcount);
  Int_t kRunNumber = fESD->GetRunNumber();
  
  // for selected events, fill some histos:
  (dynamic_cast<TH1D *>(fOutputListSupportHistos->At(0)))->Fill(0.);
  (dynamic_cast<TH1D *>(fOutputListSupportHistos->At(1)))->Fill(centralityF);
  (dynamic_cast<TH1D *>(fOutputListSupportHistos->At(4)))->Fill(vtxZGlobal);//hVertexZ
  (dynamic_cast<TH1D *>(fOutputListSupportHistos->At(5)))->Fill(nCtrb);//hNvertexCtrb
  //hNTrksEvent_cent->Fill(fESD->GetNumberOfTracks(),centralityF);
  Int_t Nacc = AliDielectronHelper::GetNacc(fESD);
  Double_t sigmaEleITS_Raw;
  Double_t sigmaEleTPC_Raw;
  //
  // store if there is at least one cutset to be used for Prefilter efficiency determination.
  // store for each cutset if there is at least one electron selected by 'fvExtraTrackCuts'.
  //
  Bool_t              atleastOnePrefilterSetting=kFALSE;
  Bool_t              atleastOneEleExtraSelected=kFALSE;
  std::vector<Bool_t> fvAtleastOneEleExtra_perCut;
  for (UInt_t iCut=0; iCut<GetNCutsets(); ++iCut) {
    if (fvDoPrefilterEff.at(iCut) == kTRUE) atleastOnePrefilterSetting=kTRUE;
    fvAtleastOneEleExtra_perCut.push_back(kFALSE); // initialize vector contents
  }
  
  //
  // store all accepted track IDs in vectors to use them for the Prefilter efficiency determination.
  //
  std::vector<Int_t> vecTrkID;
  vecTrkID.push_back(0); // first element is neccessary to hook track vector into cut vector.
  std::vector< std::vector<Int_t> > vecEleCand_perCut; // one vector per cut setting
  for (UInt_t iCut=0; iCut<GetNCutsets(); ++iCut) { vecEleCand_perCut.push_back(vecTrkID); }
  
  
  Int_t NEleSelected = 0;
  UInt_t ResolutionMask = 1000;
  if(fResolutionCuts) ResolutionMask=(1<<fResolutionCuts->GetCuts()->GetEntries())-1;
  if(mcEvent){
    Int_t nMCtracks = mcEvent->GetNumberOfTracks();
    //AliStack *fStack = mcEvent->Stack();
    
    for(Int_t iMCtrack = 0; iMCtrack < nMCtracks; iMCtrack++){
      // New:
      // Select electrons based on an attached AliDielectronSignalMC.
      // Only the first leg (branch==1) of the first attached signal will be checked.
      //
      Bool_t truth1 = AliDielectronMC::Instance()->IsMCTruth(iMCtrack, (AliDielectronSignalMC*)fSignalsMC->At(0), 1);
      if (!truth1) continue;
      
      AliMCParticle *mctrack = dynamic_cast<AliMCParticle *>(mcEvent->GetTrack(iMCtrack));
      if (!mctrack) continue;
      Double_t mcP     = mctrack->P();
      Double_t mcPt    = mctrack->Pt();
      Double_t mcEta   = mctrack->Eta();
      Double_t mcTheta = mctrack->Theta();
      Double_t mcPhi   = mctrack->Phi();
      if(mcPt  < fPtMinGEN  || mcPt  > fPtMaxGEN)  continue;
      if(mcEta < fEtaMinGEN || mcEta > fEtaMaxGEN) continue;
      if(fCalcEfficiencyGen){
        if(mctrack->Charge() < 0) fNgen_Ele->Fill(mcPt,mcEta,mcPhi);
        else                      fNgen_Pos->Fill(mcPt,mcEta,mcPhi);
      }
      Bool_t bFilled(kFALSE);
      for (Int_t iTracks = 0; iTracks < fESD->GetNumberOfTracks(); iTracks++){
        fSelectedByCut = 0;
        fSelectedByExtraCut = 0;
        AliESDtrack* track = fESD->GetTrack(iTracks);
        if (!track) { Printf("ERROR: Could not receive track %d", iTracks); continue; }
        Int_t label = track->GetLabel();
        Int_t abslabel = TMath::Abs( track->GetLabel() );
        if(abslabel != iMCtrack) continue;
        if(mctrack->Charge()/3 != track->Charge()) continue;
        
        Double_t trackPt  = track->Pt();
        Double_t trackEta = track->Eta();
        Double_t trackPhi = track->Phi();
        UInt_t ResolutionMaskSelected = 999;
        if(fResolutionCuts) ResolutionMaskSelected = fResolutionCuts->IsSelected(track);
        if(fCalcEfficiencyRec && !bFilled && ResolutionMask == ResolutionMaskSelected){
          if(mctrack->Charge() < 0) fNgen2_Rec_Ele ->Fill(trackPt,trackEta,trackPhi);
          else                      fNgen2_Rec_Pos ->Fill(trackPt,trackEta,trackPhi);
          bFilled = kTRUE;
        }
        for (UInt_t iCut=0; iCut<GetNCutsets(); ++iCut){ // loop over all specified cutInstances
          //cutting logic taken from AliDielectron::FillTrackArrays()
          UInt_t selectedMask=(1<<fvTrackCuts.at(iCut)->GetCuts()->GetEntries())-1;
          //apply track cuts
          UInt_t cutMask=fvTrackCuts.at(iCut)->IsSelected(track);
          if (cutMask!=selectedMask) continue;
          vecEleCand_perCut.at(iCut).push_back(iTracks);
          if(track->Charge() < 0){
            if(fCalcEfficiencyGen){
              fvReco_Ele.at(iCut)->Fill(mcPt,mcEta,mcPhi);
              if(fCalcEfficiencyPoslabel && label > 0) fvReco_Ele_poslabel.at(iCut)->Fill(mcPt,mcEta,mcPhi);
            }
            if(fCalcEfficiencyRec){
              fvReco_Rec_Ele.at(iCut)->Fill(trackPt,trackEta,trackPhi);
              if(fCalcEfficiencyPoslabel && label > 0) fvReco_Rec_Ele_poslabel.at(iCut)->Fill(trackPt,trackEta,trackPhi);
            }
          }
          else{
            if(fCalcEfficiencyGen){
              fvReco_Pos.at(iCut)->Fill(mcPt,mcEta,mcPhi);
              if(fCalcEfficiencyPoslabel && label > 0) fvReco_Pos_poslabel.at(iCut)->Fill(mcPt,mcEta,mcPhi);
            }
            if(fCalcEfficiencyRec){
              fvReco_Rec_Pos.at(iCut)->Fill(trackPt,trackEta,trackPhi);
              if(fCalcEfficiencyPoslabel && label > 0) fvReco_Rec_Pos_poslabel.at(iCut)->Fill(trackPt,trackEta,trackPhi);
            }
          }
          fSelectedByCut|=1<<(iCut); // store bitwise which cut settings the track survived.
          // store infos related to prefilter efficiency
          if(fvDoPrefilterEff.at(iCut) == kTRUE){
            // check if one of the prefilter electrons also survives the analysis cuts of this cutset.
            UInt_t selectedMaskExtra = (1<<fvExtraTrackCuts.at(iCut)->GetCuts()->GetEntries())-1;
            UInt_t cutMaskExtra      = fvExtraTrackCuts.at(iCut)->IsSelected(track);
            if(cutMaskExtra==selectedMaskExtra){
              fSelectedByExtraCut|=1<<(iCut); // this is just for the tree.
              atleastOneEleExtraSelected = kTRUE;
              fvAtleastOneEleExtra_perCut.at(iCut) = kTRUE;
            }
          }
        } // cut loop
        if(fSelectedByCut==0) continue;// only go on if the track survived at least 1 of the cut settings!
        // get track information
        // feel free to add more information to the tree, some more variables are already defined as branches...
        pxESD = track->Px();
        pyESD = track->Py();
        pzESD = track->Pz();
        chargeT   = track->Charge();
        signalITS = track->GetITSsignal();
        signalTPC = track->GetTPCsignal();
        kNclsTPC  = track->GetTPCNcls();
        // done below. //kTPCchi2Cl  = track->GetTPCchi2() / kNclsTPC; // only for ESD! // AOD like this: particle->Chi2perNDF()*(tpcNcls-5)/tpcNcls
        kNclsTPCdEdx = track->GetTPCsignalN(); // ("fTPCsignalN") number of points used for dEdx - maybe not in AOD...?
        kNFclsTPCr = track->GetTPCClusterInfo(2,1); // number of findable clusters(crossed rows) in the TPC with more robust definition //ESD & AOD!
        Float_t kNFclsTPC = track->GetTPCNclsF(); //tpcClusFindable // number of findable clusters in the TPC //ESD & AOD!
        kNFclsTPCfCross = (kNFclsTPC>0)?(kNFclsTPCr/kNFclsTPC):0; // fraction crossed rows/findable clusters in the TPC, as done in AliESDtrackCuts  //ESD & AOD!
        //kNtrkltsTRD = 0; // only exists for ESD, see below
        kNtrkltsTRDPID = track->GetTRDntrackletsPID();
        // @TODO: the variables above could also be retreived from the values[] array below.
        // fill the AliDielectronVarManager to get some specific variables, e.g. which were used for track selection.
        Double_t values[AliDielectronVarManager::kNMaxValues]={0.};
        AliDielectronVarManager::SetFillMap(fUsedVars); // currently filled manually in the constructor of this task.
        AliDielectronVarManager::Fill(track, values);
        // retrieve some of these variables from VarManager:
        sigmaEleITS_Raw = values[AliDielectronVarManager::kITSnSigmaEleRaw];
        sigmaEleTPC_Raw = values[AliDielectronVarManager::kTPCnSigmaEleRaw];
        sigmaEleITS = values[AliDielectronVarManager::kITSnSigmaEle];
        sigmaEleTPC = values[AliDielectronVarManager::kTPCnSigmaEle];
        sigmaEleTOF = values[AliDielectronVarManager::kTOFnSigmaEle];
        kNclsITS    = values[AliDielectronVarManager::kNclsITS];
        kITSchi2Cl  = values[AliDielectronVarManager::kITSchi2Cl];
        kNclsSITS   = values[AliDielectronVarManager::kNclsSITS];
        kTPCchi2Cl  = values[AliDielectronVarManager::kTPCchi2Cl];
        kNclsSTPC   = values[AliDielectronVarManager::kNclsSTPC];
        // TODO: which momentum is better?
        const AliExternalTrackParam *innerParam = track->GetTPCInnerParam();
        pTPC        = (innerParam) ? innerParam->P() : -1.; //Track parameters estimated at the inner wall of TPC using the TPC stand-alone
        // get MC track information
        pdgT          = mctrack->PdgCode();
        labelmotherT  = mctrack->Particle()->GetFirstMother();
        AliMCParticle *mother = 0x0;
        AliMCParticle *grandmother = 0x0;
        pdgmotherT            =   0;
        labelgrandmotherT     =  -1;
        pdggrandmotherT       =   0;
        if(labelmotherT>=0)
          mother = dynamic_cast<AliMCParticle *>(mcEvent->GetTrack(labelmotherT));
        if(mother){
          pdgmotherT = mother->PdgCode();
          labelgrandmotherT = mother->Particle()->GetFirstMother();
        }
        if(labelgrandmotherT>=0)
          grandmother = (dynamic_cast<AliMCParticle *>(mcEvent->GetTrack(labelgrandmotherT)));
        if(grandmother)
          pdggrandmotherT = grandmother->PdgCode();
        
        if(fWriteTree) // variable 'fSelectedByCut' is stored in the Tree to distinguish the cut settings...
          tracksT->Fill(); //Fill Track Tree
        
        if(fSelectedByCut & 1<<fSupportedCutInstance){ // only go on if the track survived in the cut setting for which you want to fill the support histograms.
          NEleSelected++;
          Float_t pESD = track->P();
          // for selected tracks, fill some histos:
          (dynamic_cast<TH2D *>(fOutputListSupportHistos->At(6)))->Fill(values[AliDielectronVarManager::kPIn], pTPC);//hPInVarMgr_PInStandAlone
          (dynamic_cast<TH1D *>(fOutputListSupportHistos->At(7)))->Fill(track->Pt());//hPt (reco)
          (dynamic_cast<TH2D *>(fOutputListSupportHistos->At(8)))->Fill(pESD, pTPC);//hP_PIn
          (dynamic_cast<TH2D *>(fOutputListSupportHistos->At(9)))->Fill(pESD, mctrack->P());//hP_Pgen
          // ITS
          (dynamic_cast<TH2D *>(fOutputListSupportHistos->At(10)))->Fill(pESD, signalITS);
          (dynamic_cast<TH2D *>(fOutputListSupportHistos->At(11)))->Fill(pESD, sigmaEleITS);
          //(dynamic_cast<TH2D *>(fOutputListSupportHistos->At(12)))->Fill(pESD, );
          //(dynamic_cast<TH2D *>(fOutputListSupportHistos->At(13)))->Fill(pESD, );
          // TPC
          (dynamic_cast<TH2D *>(fOutputListSupportHistos->At(14)))->Fill(pTPC, signalTPC);
          (dynamic_cast<TH2D *>(fOutputListSupportHistos->At(15)))->Fill(pTPC, sigmaEleTPC);
          (dynamic_cast<TH3D *>(fOutputListSupportHistos->At(16)))->Fill(pTPC, sigmaEleTPC, signalTPC);
          // run dependency
          (dynamic_cast<TH3D *>(fOutputListSupportHistos->At(22)))->Fill(pTPC, signalTPC, kRunNumber);//hTPC_dEdx_P_run
          (dynamic_cast<TH3D *>(fOutputListSupportHistos->At(23)))->Fill(pTPC, sigmaEleTPC, kRunNumber);//hTPCnSigmaEle_P_run
          (dynamic_cast<TH3D *>(fOutputListSupportHistos->At(24)))->Fill(pESD, sigmaEleITS, kRunNumber);//hITSnSigmaEle_P_run
          // TOF
          (dynamic_cast<TH2D *>(fOutputListSupportHistos->At(25)))->Fill(pTPC, sigmaEleTOF);//hTOFnSigmaEle_P
          //(dynamic_cast<TH2D *>(fOutputListSupportHistos->At(26)))->Fill(pTPC, );
          // Eta and Phi
          (dynamic_cast<TH1D *>(fOutputListSupportHistos->At(27)))->Fill(mcEta);
          (dynamic_cast<TH1D *>(fOutputListSupportHistos->At(28)))->Fill(mcPhi);
          (dynamic_cast<TH2D *>(fOutputListSupportHistos->At(29)))->Fill(mcEta, mcPhi);
          (dynamic_cast<TH2D *>(fOutputListSupportHistos->At(30)))->Fill(mcEta, signalTPC);
          (dynamic_cast<TH3D *>(fOutputListSupportHistos->At(31)))->Fill(mcEta, signalTPC, pTPC);
          (dynamic_cast<TH2D *>(fOutputListSupportHistos->At(32)))->Fill(mcEta, sigmaEleTPC);//hTPCnSigmaEle_Eta
          (dynamic_cast<TH3D *>(fOutputListSupportHistos->At(33)))->Fill(mcEta, sigmaEleTPC, pTPC);//hTPCnSigmaEle_Eta_P
          //(dynamic_cast<TH3D *>(fOutputListSupportHistos->At(34)))->Fill(mcEta, sigmaEleTPC, );//hTPCnSigmaEle_Eta_RefMultTPConly
          (dynamic_cast<TH3D *>(fOutputListSupportHistos->At(35)))->Fill(mcEta, sigmaEleTPC, Nacc);//hTPCnSigmaEle_Eta_Nacc
          // DCA
          //(dynamic_cast<TH1D *>(fOutputListSupportHistos->At(36)))->Fill();
          //(dynamic_cast<TH1D *>(fOutputListSupportHistos->At(37)))->Fill();
          //(dynamic_cast<TH2D *>(fOutputListSupportHistos->At(38)))->Fill();
          // Quality
          (dynamic_cast<TH1D *>(fOutputListSupportHistos->At(39)))->Fill(kNFclsTPCfCross);
          (dynamic_cast<TH1D *>(fOutputListSupportHistos->At(40)))->Fill(kNFclsTPCr);
          (dynamic_cast<TH1D *>(fOutputListSupportHistos->At(41)))->Fill(kNclsTPC);
          (dynamic_cast<TH1D *>(fOutputListSupportHistos->At(42)))->Fill(kNclsITS);
          (dynamic_cast<TH1D *>(fOutputListSupportHistos->At(43)))->Fill(kTPCchi2Cl);
          (dynamic_cast<TH1D *>(fOutputListSupportHistos->At(44)))->Fill(kITSchi2Cl);
          (dynamic_cast<TH1D *>(fOutputListSupportHistos->At(45)))->Fill(kNclsSTPC);
          (dynamic_cast<TH1D *>(fOutputListSupportHistos->At(46)))->Fill(kNclsSITS);
          (dynamic_cast<TH1D *>(fOutputListSupportHistos->At(47)))->Fill(values[AliDielectronVarManager::kNclsSFracTPC]);
          (dynamic_cast<TH1D *>(fOutputListSupportHistos->At(48)))->Fill(values[AliDielectronVarManager::kNclsSFracITS]);
          (dynamic_cast<TH1D *>(fOutputListSupportHistos->At(49)))->Fill(values[AliDielectronVarManager::kTPCclsDiff]);
          (dynamic_cast<TH1D *>(fOutputListSupportHistos->At(50)))->Fill(kNclsTPCdEdx);
          (dynamic_cast<TH2D *>(fOutputListSupportHistos->At(51)))->Fill(kNclsTPC, kNFclsTPCr);
          (dynamic_cast<TH2D *>(fOutputListSupportHistos->At(52)))->Fill(track->Pt(), kNFclsTPCr);
          (dynamic_cast<TH2D *>(fOutputListSupportHistos->At(60)))->Fill(mcEta, sigmaEleITS);//hITSnSigmaEle_Eta
          (dynamic_cast<TH3D *>(fOutputListSupportHistos->At(61)))->Fill(mcEta, sigmaEleITS, pESD);//hITSnSigmaEle_Eta_P
          //(dynamic_cast<TH3D *>(fOutputListSupportHistos->At(62)))->Fill(mcEta, sigmaEleITS, );//hITSnSigmaEle_Eta_RefMultTPConly
          (dynamic_cast<TH3D *>(fOutputListSupportHistos->At(63)))->Fill(mcEta, sigmaEleITS, Nacc);//hITSnSigmaEle_Eta_Nacc
          // check of post PID correction. raw values from AliDielectronVarManager
          (dynamic_cast<TH2D *>(fOutputListSupportHistos->At(64)))->Fill(pESD, sigmaEleITS_Raw);//hITSnSigmaEleRaw_P
          (dynamic_cast<TH2D *>(fOutputListSupportHistos->At(65)))->Fill(pTPC, sigmaEleTPC_Raw);//hTPCnSigmaEleRaw_P
          (dynamic_cast<TH2D *>(fOutputListSupportHistos->At(66)))->Fill(mcEta, sigmaEleITS_Raw);//hITSnSigmaEleRaw_Eta
          (dynamic_cast<TH2D *>(fOutputListSupportHistos->At(67)))->Fill(mcEta, sigmaEleTPC_Raw);//hTPCnSigmaEleRaw_Eta
          // PDG codes
          (dynamic_cast<TH1D *>(fOutputListSupportHistos->At(70)))->Fill(pdgT);//hPdgCode
          (dynamic_cast<TH1D *>(fOutputListSupportHistos->At(71)))->Fill(pdgmotherT);//hPdgCodeM
          (dynamic_cast<TH1D *>(fOutputListSupportHistos->At(72)))->Fill(pdggrandmotherT);//hPdgCodeGM
          (dynamic_cast<TH2D *>(fOutputListSupportHistos->At(73)))->Fill(pdgmotherT, pdggrandmotherT);//hPdgCodeM_GM
        } //fSupportedCutInstance
      } // reco track loop
      
      if(fCalcEfficiencyRec){
        // smear generated quantities of non-reconstructed particles with external response matrix,
        // in order to compute the efficiency consistently based on "measurable" quantities.
        Double_t trackP     = mcP;
        Double_t trackPt    = mcPt;
        Double_t trackEta   = mcEta;
        Double_t trackTheta = mcTheta;
        Double_t trackPhi   = mcPhi;
        
        if(fPhiEleResArr && fPhiPosResArr){
          Double_t phiSmearing = 0.;
          if(mctrack->Charge() < 0) phiSmearing =  GetSmearing(fPhiEleResArr,mcP);
          else                      phiSmearing =  GetSmearing(fPhiPosResArr,mcP);
          trackPhi = mcPhi + phiSmearing;
        }
        if(fThetaResArr){
          Double_t thetaSmearing = GetSmearing(fThetaResArr,mcP);
          trackTheta = mcTheta + thetaSmearing;
          trackEta = -TMath::Log(TMath::Tan(trackTheta/2.));
        } else if(fEtaResArr){
          Double_t etaSmearing = GetSmearing(fEtaResArr,mcP);
          trackEta = mcEta + etaSmearing;
        }
        if(fPResArr){
          Double_t pSmearing = GetSmearing(fPResArr,mcP);
          if(fUseRelPResolution) trackP  = mcP * pSmearing;
          else  trackP  = mcP + pSmearing;
          trackPt = TMath::Sin(trackTheta) * trackP;
        }
        
        if(mctrack->Charge() < 0){
          fNgen1_Rec_Ele->Fill(trackPt,trackEta,trackPhi);
          if(!bFilled)
            fNgen2_Rec_Ele->Fill(trackPt,trackEta,trackPhi);
        }
        else{
          fNgen1_Rec_Pos->Fill(trackPt,trackEta,trackPhi);
          if(!bFilled)
            fNgen2_Rec_Pos->Fill(trackPt,trackEta,trackPhi);
        }
      }
      
    } // mc track loop
    
    //(dynamic_cast<TH2D *>(fOutputListSupportHistos->At(2)))->Fill(centralityF,0);//hNTrksEvent_cent
    (dynamic_cast<TH2D *>(fOutputListSupportHistos->At(3)))->Fill(centralityF,NEleSelected);//hNEleEvent_cent
    if(atleastOnePrefilterSetting && atleastOneEleExtraSelected) CalcPrefilterEff(mcEvent, vecEleCand_perCut, fvAtleastOneEleExtra_perCut);
    
    if(fDoPairing){  // calculate pair efficiency from signal pairs
      std::vector<LMEEparticle> LMEEelectrons,LMEEpositrons;
      for(Int_t iMC = 0; iMC < nMCtracks; iMC++){
        AliMCParticle *part = dynamic_cast<AliMCParticle*> (mcEvent->GetTrack(iMC));
        if(!part){ Printf("no MCtrack: %d", iMC); continue; }
        if(!AliDielectronMC::Instance()->IsMCTruth(iMC, (AliDielectronSignalMC*)fSignalsMC->At(0), 1)) continue;
        Int_t mLab = part->GetMother();
        AliMCParticle *mother(0x0);
        if(mLab >= 0) mother = dynamic_cast<AliMCParticle*> (mcEvent->GetTrack(mLab));
        if(!mother){ Printf("no mother for MCtrack: %d", iMC); continue; }
        Int_t mPDG = mother->PdgCode();
        Int_t grmLab = mother->GetMother();
        AliMCParticle *grmother(0x0);
        if(grmLab >= 0) grmother = dynamic_cast<AliMCParticle*> (mcEvent->GetTrack(grmLab));
        Int_t grmPDG(0);
        if(grmother) grmPDG = grmother->PdgCode();
        LMEEparticle lmeeLeg(GetNCutsets());
        lmeeLeg.genP     = part->P();
        lmeeLeg.genPt    = part->Pt();
        lmeeLeg.genTheta = part->Theta();
        lmeeLeg.genEta   = part->Eta();
        lmeeLeg.genPhi   = part->Phi();
        lmeeLeg.mlabel   = mLab;
        lmeeLeg.mPDG     = mPDG;
        lmeeLeg.grmlabel = grmLab;
        lmeeLeg.grmPDG   = grmPDG;
        for(Int_t iESD = 0; iESD < fESD->GetNumberOfTracks(); ++iESD){
          AliESDtrack *track = fESD->GetTrack(iESD);
          if(!track){ Printf("ERROR: Could not receive track %d", iESD); continue; }
          if(TMath::Abs(track->GetLabel()) != iMC) continue;
          for(UInt_t iCut=0; iCut<fvTrackCuts.size(); ++iCut){ // loop over all specified cutInstances
            UInt_t selectedMask=(1<<fvTrackCuts.at(iCut)->GetCuts()->GetEntries())-1;
            UInt_t cutMask=fvTrackCuts.at(iCut)->IsSelected(track);
            if(cutMask==selectedMask)
              lmeeLeg.vbRec.at(iCut) = kTRUE;
          } // cut loop
        } // ESD loop
        if(fCalcEfficiencyRec){ 
          if(fPhiEleResArr && fPhiPosResArr){
            Double_t phiSmearing = 0.;
            if(part->Charge() < 0) phiSmearing =  GetSmearing(fPhiEleResArr,lmeeLeg.genP);
            else                   phiSmearing =  GetSmearing(fPhiPosResArr,lmeeLeg.genP);
            lmeeLeg.recPhi = lmeeLeg.genPhi + phiSmearing;
          }
          if(fThetaResArr){
            Double_t thetaSmearing = GetSmearing(fThetaResArr,lmeeLeg.genP);
            lmeeLeg.recTheta = lmeeLeg.genTheta + thetaSmearing;
            lmeeLeg.recEta = -TMath::Log(TMath::Tan(lmeeLeg.recTheta/2.));
          } else if(fEtaResArr){
            Double_t etaSmearing = GetSmearing(fEtaResArr,lmeeLeg.genP);
            lmeeLeg.recEta = lmeeLeg.genEta + etaSmearing;
          }
          if(fPResArr){
            Double_t pSmearing = GetSmearing(fPResArr,lmeeLeg.genP);
            if(fUseRelPResolution) lmeeLeg.recP  = lmeeLeg.genP * pSmearing;
            else                   lmeeLeg.recP  = lmeeLeg.genP + pSmearing;
            lmeeLeg.recPt = TMath::Sin(lmeeLeg.recTheta) * lmeeLeg.recP;
          }
          lmeeLeg.MakeRecLV();
        }
        lmeeLeg.MakeGenLV();
        if(part->PdgCode() == -11)
          LMEEpositrons.push_back(lmeeLeg);
        else if(part->PdgCode() == 11)
          LMEEelectrons.push_back(lmeeLeg);
      } // MC track loop
      if(LMEEelectrons.size() > 0 && LMEEpositrons.size() > 0){
        for(std::vector<LMEEparticle>::iterator it1 = LMEEelectrons.begin(); it1 != LMEEelectrons.end(); ++it1){
          Bool_t bAccGen1 = it1->genPt > 0.2 && TMath::Abs(it1->genEta) < 0.8;
          Bool_t bAccRec1 = it1->recPt > 0.2 && TMath::Abs(it1->recEta) < 0.8;
          if(!(bAccGen1 || bAccRec1)) continue;
          Bool_t charm1  =  TMath::Abs(it1->mPDG) > 400 && TMath::Abs(it1->mPDG) < 440 && !(TMath::Abs(it1->grmPDG) > 500 && TMath::Abs(it1->grmPDG) < 550);
          Bool_t beauty1 = (TMath::Abs(it1->mPDG) > 400 && TMath::Abs(it1->mPDG) < 440 && TMath::Abs(it1->grmPDG) > 500  && TMath::Abs(it1->grmPDG) < 550) || (TMath::Abs(it1->mPDG) > 500 && TMath::Abs(it1->mPDG) < 550);
          for(std::vector<LMEEparticle>::iterator it2 = LMEEpositrons.begin(); it2 != LMEEpositrons.end(); ++it2){
            Bool_t bAccGen2 = it2->genPt > 0.2 && TMath::Abs(it2->genEta) < 0.8;
            Bool_t bAccRec2 = it2->recPt > 0.2 && TMath::Abs(it2->recEta) < 0.8;
            if(!(bAccGen2 || bAccRec2)) continue;
            Bool_t charm2  =  TMath::Abs(it2->mPDG) > 400 && TMath::Abs(it2->mPDG) < 440 && !(TMath::Abs(it2->grmPDG) > 500 && TMath::Abs(it2->grmPDG) < 550);
            Bool_t beauty2 = (TMath::Abs(it2->mPDG) > 400 && TMath::Abs(it2->mPDG) < 440 && TMath::Abs(it2->grmPDG) > 500  && TMath::Abs(it2->grmPDG) < 550) || (TMath::Abs(it2->mPDG) > 500 && TMath::Abs(it2->mPDG) < 550);
            if(fCalcEfficiencyGen && bAccGen1 && bAccGen2){
              Double_t mee  = (it1->genLv + it2->genLv).M();
              Double_t ptee = (it1->genLv + it2->genLv).Pt();
              if(it1->mlabel == it2->mlabel){
                fNgenPairsResonances->Fill(mee,ptee);
                for(UInt_t iCut=0; iCut<fvTrackCuts.size(); ++iCut)
                  if(it1->vbRec.at(iCut) && it2->vbRec.at(iCut)) fvRecoPairsResonances.at(iCut)->Fill(mee,ptee);
              }
              else{
                fNgenPairsDiffMothers->Fill(mee,ptee);
                for(UInt_t iCut=0; iCut<fvTrackCuts.size(); ++iCut)
                  if(it1->vbRec.at(iCut) && it2->vbRec.at(iCut)) fvRecoPairsDiffMothers.at(iCut)->Fill(mee,ptee);
                if(charm1 && charm2){
                  fNgenPairsCharm ->Fill(mee,ptee);
                  for(UInt_t iCut=0; iCut<fvTrackCuts.size(); ++iCut)
                    if(it1->vbRec.at(iCut) && it2->vbRec.at(iCut)) fvRecoPairsCharm.at(iCut)->Fill(mee,ptee);
                }
                if(beauty1 && beauty2){
                  fNgenPairsBeauty->Fill(mee,ptee);
                  for(UInt_t iCut=0; iCut<fvTrackCuts.size(); ++iCut)
                    if(it1->vbRec.at(iCut) && it2->vbRec.at(iCut)) fvRecoPairsBeauty.at(iCut)->Fill(mee,ptee);
                }
                if((charm1 || beauty1) && (charm2 || beauty2)){
                  fNgenPairsHF->Fill(mee,ptee);
                  for(UInt_t iCut=0; iCut<fvTrackCuts.size(); ++iCut)
                    if(it1->vbRec.at(iCut) && it2->vbRec.at(iCut)) fvRecoPairsHF.at(iCut)->Fill(mee,ptee);
                }
              }
            }
            if(fCalcEfficiencyRec && bAccRec1 && bAccRec2){
              Double_t mee  = (it1->recLv + it2->recLv).M();
              Double_t ptee = (it1->recLv + it2->recLv).Pt();
              if(it1->mlabel == it2->mlabel){
                fNgenPairsRecResonances->Fill(mee,ptee);
                for(UInt_t iCut=0; iCut<fvTrackCuts.size(); ++iCut)
                  if(it1->vbRec.at(iCut) && it2->vbRec.at(iCut)) fvRecoPairsRecResonances.at(iCut)->Fill(mee,ptee);
              }
              else{
                fNgenPairsRecDiffMothers->Fill(mee,ptee);
                for(UInt_t iCut=0; iCut<fvTrackCuts.size(); ++iCut)
                  if(it1->vbRec.at(iCut) && it2->vbRec.at(iCut)) fvRecoPairsRecDiffMothers.at(iCut)->Fill(mee,ptee);
                if(charm1 && charm2){
                  fNgenPairsRecCharm ->Fill(mee,ptee);
                  for(UInt_t iCut=0; iCut<fvTrackCuts.size(); ++iCut)
                    if(it1->vbRec.at(iCut) && it2->vbRec.at(iCut)) fvRecoPairsRecCharm.at(iCut)->Fill(mee,ptee);
                }
                if(beauty1 && beauty2){
                  fNgenPairsRecBeauty->Fill(mee,ptee);
                  for(UInt_t iCut=0; iCut<fvTrackCuts.size(); ++iCut)
                    if(it1->vbRec.at(iCut) && it2->vbRec.at(iCut)) fvRecoPairsRecBeauty.at(iCut)->Fill(mee,ptee);
                }
                if((charm1 || beauty1) && (charm2 || beauty2)){
                  fNgenPairsRecHF->Fill(mee,ptee);
                  for(UInt_t iCut=0; iCut<fvTrackCuts.size(); ++iCut)
                    if(it1->vbRec.at(iCut) && it2->vbRec.at(iCut)) fvRecoPairsRecHF.at(iCut)->Fill(mee,ptee);
                }
              }
            }
          } // positron iteration
        } // electron iteration
      } // pairing
    } // do pairing
    
    if(fCalcResolution && fResolutionCuts){ // calculate electron pt resolution
      Int_t Nprimaries = mcEvent->GetNumberOfPrimaries();
      TLorentzVector l1Gen,l2Gen,l1Rec,l2Rec;
      for(Int_t iTracks = 0; iTracks < fESD->GetNumberOfTracks(); iTracks++){
        AliESDtrack *track = fESD->GetTrack(iTracks);
        if (!track) { Printf("ERROR: Could not receive track %d", iTracks); continue; }
        if(fResolutionCuts->IsSelected(track) != ResolutionMask) continue;
        Int_t label = track->GetLabel();
        
        AliMCParticle *part = dynamic_cast<AliMCParticle *>(mcEvent->GetTrack(TMath::Abs(label)));
        if(!part) { Printf("ERROR: Could not receive mc track %d", TMath::Abs(label)); continue; }
        Int_t mcLabel = part->Label();
        if(!fStack->IsPhysicalPrimary(mcLabel)) continue;
        if(TMath::Abs(part->PdgCode()) != 11) continue;

        Double_t mcPt     = part->Pt();
        Double_t mcP      = part->P();
        Double_t mcTheta  = part->Theta();
        Double_t mcEta    = part->Eta();
        Double_t mcPhi    = part->Phi();
        Double_t recPt    = track->Pt();
        Double_t recP     = track->P();
        Double_t recTheta = track->Theta();
        Double_t recEta   = track->Eta();
        Double_t recPhi   = track->Phi();

        if(fMakeResolutionSparse){
          Double_t thnvals[6] = {mcP,recP / mcP,recP - mcP,recTheta - mcTheta,recEta - mcEta,recPhi - mcPhi};
          if(TMath::Abs(mcEta) < 0.9) {
            if(part->Charge()<0){ 
              if(label > 0)
                fTHnResElectrons1->Fill(thnvals);
              else
                fTHnResElectrons2->Fill(thnvals);
            }
            if(part->Charge()>0){ 
              if(label > 0)
                fTHnResPositrons1->Fill(thnvals);
              else
                fTHnResPositrons2->Fill(thnvals);
            }
          }
        }
        else{
          if(TMath::Abs(mcEta) < 0.9) {
            fPGen                ->Fill(mcP);
            fPRec                ->Fill(recP);
            fPGen_DeltaP         ->Fill(mcP, recP - mcP);
            fPtGen_DeltaPt       ->Fill(mcPt,recPt - mcPt);
            fPGen_PrecOverPGen   ->Fill(mcP, recP / mcP);
            fPtGen_PtRecOverPtGen->Fill(mcPt,recPt / mcPt);
            if (part->Charge()<0) fPGen_DeltaPhi_Ele->Fill(mcP, recPhi - mcPhi);
            else                  fPGen_DeltaPhi_Pos->Fill(mcP, recPhi - mcPhi);
            fPhiGen_DeltaPhi     ->Fill(part->Phi(),  recPhi - mcPhi);
          }
          if(mcPt > 0.2){
            fPGen_DeltaEta       ->Fill(mcP, recEta - mcEta);
            fPGen_DeltaTheta     ->Fill(mcP, recTheta - mcTheta);
            fEtaGen_DeltaEta     ->Fill(part->Eta(),  recEta - mcEta);
            fThetaGen_DeltaTheta ->Fill(part->Theta(),recTheta - mcTheta);
            fDeltaPhi->Fill(recPhi - mcPhi);
          }
        }
        
        if(label < 0) continue;
        l1Gen.SetPtEtaPhiM(part ->Pt(),part ->Eta(),part ->Phi(),AliPID::ParticleMass(AliPID::kElectron));
        l1Rec.SetPtEtaPhiM(track->Pt(),track->Eta(),track->Phi(),AliPID::ParticleMass(AliPID::kElectron));
        for(Int_t iTracks2 = iTracks; iTracks2 < fESD->GetNumberOfTracks(); iTracks2++){
          if(iTracks == iTracks2) continue;
          AliESDtrack *track2 = fESD->GetTrack(iTracks2);
          if (!track2) { Printf("ERROR: Could not receive track %d", iTracks); continue; }
          if(fResolutionCuts->IsSelected(track2) != ResolutionMask) continue;
          Int_t label2 = track2->GetLabel();
          if(label == label2) continue;
          if(label2 < 0) continue;
          AliMCParticle *part2 = dynamic_cast<AliMCParticle *>(mcEvent->GetTrack(TMath::Abs(label2)));
          if(!part2) { Printf("ERROR: Could not receive mc track %d", TMath::Abs(label)); continue; }
          if(track2->Charge() != part2->Charge()/3) continue;
          Int_t mcLabel2 = part2->Label();
          if(!fStack->IsPhysicalPrimary(mcLabel2)) continue;
          if(TMath::Abs(part2->PdgCode()) != 11) continue;
          l2Gen.SetPtEtaPhiM(part2 ->Pt(),part2 ->Eta(),part2 ->Phi(),AliPID::ParticleMass(AliPID::kElectron));
          l2Rec.SetPtEtaPhiM(track2->Pt(),track2->Eta(),track2->Phi(),AliPID::ParticleMass(AliPID::kElectron));
          Double_t OpeningAngleGen = l1Gen.Angle(l2Gen.Vect());
          Double_t OpeningAngleRec = l1Rec.Angle(l2Rec.Vect());
          Double_t OpeningAngleResolution = (OpeningAngleGen > 0.) ? OpeningAngleRec/OpeningAngleGen : -1.;
          
          Double_t vals[4];
          vals[0] = (l1Gen+l2Gen).M();
          vals[2] = (vals[0] > 0.) ? (l1Rec+l2Rec).M()/vals[0] : -1.;
          
          vals[1] = (l1Gen+l2Gen).Pt();
          vals[3] = (vals[1] > 0.) ? (l1Rec+l2Rec).Pt()/vals[1] : -1.;
          if(part->Charge() != part2->Charge()){
            fOpeningAngleGen_DeltaOpeningAngleUS ->Fill(OpeningAngleGen,OpeningAngleRec - OpeningAngleGen);
            fMgen_PtGen_mRes_ptRes->Fill(vals);
          }
          else
            fOpeningAngleGen_DeltaOpeningAngleLS ->Fill(OpeningAngleGen,OpeningAngleRec - OpeningAngleGen);
        } // pairing loop
      } // track loop
    } // resolution calculation
  } //MC loop
  //Printf("__________ end of Event ( %i ) __________", fEventcount);
}

//________________________________________________________________________
void AliAnalysisTaskElectronEfficiency::CalcPrefilterEff(AliMCEvent* mcEventLocal, const std::vector< std::vector<Int_t> > & vvEleCand, const std::vector<Bool_t> & vbEleExtra)
{
  /// Determination of random electron rejection efficiency due to pair-prefiltering (used for photon conversion + Dalitz rejection).
  /// It is estimated by pairing a sample of "testparticles" (using primary, non-injected, charged pions) with the electrons of the prefilter sample (stored in 'vvEleCand').
  /// With this method the pair has no real correlation, just like in the case of random rejection of a signal electron.
  /// The pairing is only done for cut settings for which in the current event at least one final analysis electron was found (stored in 'vbEleExtra').
  /// Probably a second electron of any type should be in the event, otherwise the chance of rejecting an electron-electron pair is zero. Controlled by variable 'fNminEleInEventForRej' (=2 by default).
  //  @TODO: question: how to treat ULS and LS pairing in this context?
  /// The prefilter pair cuts are applied to these random pairs. All and rejected testparticles are stored in 3D histograms.
  /// The advantage of this method is a huge number of random pairs, the disadvantage may be a systematic difference in the implicitly available number of pairs...
  /// (by Patrick)
  
  Int_t nMCtracks  = mcEventLocal->GetNumberOfTracks();
  AliStack *fStack = mcEventLocal->Stack();
  Double_t eleMass = AliPID::ParticleMass(AliPID::kElectron);
  
  for(Int_t iMCtrack = 0; iMCtrack < nMCtracks; iMCtrack++)
  {
    // select only primary, non-injected, charged pions:
    if(!fStack->IsPhysicalPrimary(iMCtrack)) continue;
    AliMCParticle *mcPion = dynamic_cast<AliMCParticle *>(mcEventLocal->GetTrack(iMCtrack));
    if( TMath::Abs(mcPion->PdgCode()) !=  211 ) continue;
    
    Double_t mcPt(-1.),mcEta(-9.),mcPhi(-9.);
    mcPt = mcPion->Pt(); mcEta = mcPion->Eta(); mcPhi = mcPion->Phi();
    if(mcPt  < fPtMinGEN  || mcPt  > fPtMaxGEN)  continue;
    if(mcEta < fEtaMinGEN || mcEta > fEtaMaxGEN) continue;
    
    // need TLorentzVectors for the pairing
    Double_t vMomPi[3]; mcPion->PxPyPz(vMomPi);
    Double_t momPi    = mcPion->P();
    Double_t energyPi = TMath::Sqrt(momPi*momPi + eleMass*eleMass);	// need to use hardcoded electron mass.
    TLorentzVector dau1;
    dau1.SetPxPyPzE(vMomPi[0],vMomPi[1],vMomPi[2],energyPi);
    Int_t chargePion = mcPion->Charge();
    
    // loop over all cut settings
    for (UInt_t iCut=0; iCut<GetNCutsets(); ++iCut)
    {
      if (vbEleExtra.at(iCut) == kFALSE) continue; // if there is no global electron in the event, no prefilter pairing should be done.
      if (vvEleCand.at(iCut).size()-1 < fNminEleInEventForRej) continue; // size-1 because vec[0] is not used (always filled with '0').
      if (fvDoPrefilterEff.at(iCut) == kFALSE) continue; // may not even be needed anymore...
      // flags for the pion
      Bool_t rejByAllSigns = kFALSE;
      Bool_t rejByUnlike   = kFALSE;
      
      // pair the testparticle with all selected electron candidates and check for close pairs.
      // pairing starts at vec[1] because vec[0] is always set to 0, since this is needed to initialize the 2D vector (see UserExec).
      for (UInt_t iEle=1; iEle<vvEleCand.at(iCut).size(); iEle++)
      {
        AliMCParticle *mcEle = dynamic_cast<AliMCParticle *>(mcEventLocal->GetTrack(vvEleCand.at(iCut).at(iEle)));
        Double_t vMomEle[3]; mcEle->PxPyPz(vMomEle);
        Double_t momEle    = mcEle->P();
        Double_t energyEle = TMath::Sqrt(momEle*momEle + eleMass*eleMass);	// need to use hardcoded electron mass.
        TLorentzVector dau2;
        dau2.SetPxPyPzE(vMomEle[0],vMomEle[1],vMomEle[2],energyEle);
        Int_t chargeEle = mcEle->Charge();
        
        Double_t mee   = (dau1+dau2).M();
        Double_t ptee  = (dau1+dau2).Pt();
        Double_t theta = dau1.Angle(dau2.Vect());
        Double_t phiv  = PhivPair(fESD->GetMagneticField(), chargePion, chargeEle, dau1.Vect(), dau2.Vect());
        
        (dynamic_cast<TH2D *>(fOutputListSupportHistos->At(53)))->Fill(mee,ptee);//hMeePtee
        (dynamic_cast<TH3D *>(fOutputListSupportHistos->At(54)))->Fill(mee,phiv,theta);//hMeePhiVOpen
        (dynamic_cast<TH2D *>(fOutputListSupportHistos->At(55)))->Fill(mee,theta);//hMeeOpen
        (dynamic_cast<TH2D *>(fOutputListSupportHistos->At(56)))->Fill(mee,phiv);//hMeePhiV
        (dynamic_cast<TH2D *>(fOutputListSupportHistos->At(57)))->Fill(ptee,theta);//hPteeOpen
        (dynamic_cast<TH2D *>(fOutputListSupportHistos->At(58)))->Fill(ptee,phiv);//hPteePhiV
        (dynamic_cast<TH2D *>(fOutputListSupportHistos->At(59)))->Fill(theta,phiv);//hOpenPhiV
        
        // tag the pion as rejected, if the pair falls within following cut:
        if ( (fvRejCutMee.at(iCut)   > 0. || fvRejCutTheta.at(iCut) > 0. || fvRejCutPhiV.at(iCut) < 3.14) // at least one cut must be enabled
            && (mee   < fvRejCutMee.at(iCut)   || fvRejCutMee.at(iCut)   < 0.) // -> within cut or cut disabled
            && (theta < fvRejCutTheta.at(iCut) || fvRejCutTheta.at(iCut) < 0.)
            && (phiv  > fvRejCutPhiV.at(iCut)  || fvRejCutPhiV.at(iCut)  > 3.14) )
        {
          //cout << " pion rejected! mee = " << mee << "  theta = " << theta << "  phiv = " << phiv << endl;
          rejByAllSigns = kTRUE;
          if (chargePion!=chargeEle) rejByUnlike = kTRUE;
        }
      } //electron loop
      
      //fill histograms per cut setting to determine rejection efficiency
      fvAllPionsForRej.at(iCut)->Fill(mcPt,mcEta,mcPhi); //reference histogram (denominator)
      if (rejByAllSigns) fvPionsRejByAllSigns.at(iCut)->Fill(mcPt,mcEta,mcPhi);
      if (rejByUnlike)   fvPionsRejByUnlike.at(iCut)->Fill(mcPt,mcEta,mcPhi);
    } //cut loop
    
  } //pion loop
}
//______________________________________________
Double_t AliAnalysisTaskElectronEfficiency::PhivPair(Double_t MagField, Int_t charge1, Int_t charge2, TVector3 dau1, TVector3 dau2) //const
{
  /// Following the idea to use opening of collinear pairs in magnetic field from e.g. PHENIX
  /// to identify conversions. Angle between ee plane and magnetic field is calculated (0 to pi).
  /// Due to tracking to the primary vertex, conversions with no intrinsic opening angle
  /// always end up as pair in "cowboy" configuration. The function as defined here then
  /// returns values close to pi.
  /// Correlated Like Sign pairs (from double conversion / dalitz + conversion) may show up
  /// at pi or at 0 depending on which leg has the higher momentum. (not checked yet)
  /// This expected ambiguity is not seen due to sorting of track arrays in this framework.
  /// To reach the same result as for ULS (~pi), the legs are flipped for LS.
  /// from PWGDQ/dielectron/AliDielectronPair.cxx
  
  //Define local buffer variables for leg properties
  Double_t px1=-9999.,py1=-9999.,pz1=-9999.;
  Double_t px2=-9999.,py2=-9999.,pz2=-9999.;
  
  TVector3 fD1=dau1;
  TVector3 fD2=dau2;
  Int_t    d1Q=charge1;
  //Int_t    d2Q=charge2;
  if (fRandomizeDaughters) { // randomize daughters if requested
    if (fRandom3.Rndm()>0.5) {
      fD1=dau2;
      fD2=dau1;
      d1Q=charge2;
      //d2Q=charge1;
    }
  }
  else { // sort particles according to pt, as done by default in AliDielectronPair
    if (dau1.Pt() < dau2.Pt()) {
      fD1=dau2;
      fD2=dau1;
      d1Q=charge2;
      //d2Q=charge1;
    }
  }
  
  if (charge1*charge2 > 0.) { // Like Sign
    if(MagField<0){ // inverted behaviour
      if(d1Q>0){
        px1 = fD1.Px();   py1 = fD1.Py();   pz1 = fD1.Pz();
        px2 = fD2.Px();   py2 = fD2.Py();   pz2 = fD2.Pz();
      }else{
        px1 = fD2.Px();   py1 = fD2.Py();   pz1 = fD2.Pz();
        px2 = fD1.Px();   py2 = fD1.Py();   pz2 = fD1.Pz();
      }
    }else{
      if(d1Q>0){
        px1 = fD2.Px();   py1 = fD2.Py();   pz1 = fD2.Pz();
        px2 = fD1.Px();   py2 = fD1.Py();   pz2 = fD1.Pz();
      }else{
        px1 = fD1.Px();   py1 = fD1.Py();   pz1 = fD1.Pz();
        px2 = fD2.Px();   py2 = fD2.Py();   pz2 = fD2.Pz();
      }
    }
  }
  else { // Unlike Sign
    if(MagField>0){ // regular behaviour
      if(d1Q>0){
        px1 = fD1.Px();
        py1 = fD1.Py();
        pz1 = fD1.Pz();
        
        px2 = fD2.Px();
        py2 = fD2.Py();
        pz2 = fD2.Pz();
      }else{
        px1 = fD2.Px();
        py1 = fD2.Py();
        pz1 = fD2.Pz();
        
        px2 = fD1.Px();
        py2 = fD1.Py();
        pz2 = fD1.Pz();
      }
    }else{
      if(d1Q>0){
        px1 = fD2.Px();
        py1 = fD2.Py();
        pz1 = fD2.Pz();
        
        px2 = fD1.Px();
        py2 = fD1.Py();
        pz2 = fD1.Pz();
      }else{
        px1 = fD1.Px();
        py1 = fD1.Py();
        pz1 = fD1.Pz();
        
        px2 = fD2.Px();
        py2 = fD2.Py();
        pz2 = fD2.Pz();
      }
    }
  }
  
  Double_t px = px1+px2;
  Double_t py = py1+py2;
  Double_t pz = pz1+pz2;
  Double_t dppair = TMath::Sqrt(px*px+py*py+pz*pz);
  
  //unit vector of (pep+pem)
  Double_t pl = dppair;
  Double_t ux = px/pl;
  Double_t uy = py/pl;
  Double_t uz = pz/pl;
  Double_t ax = uy/TMath::Sqrt(ux*ux+uy*uy);
  Double_t ay = -ux/TMath::Sqrt(ux*ux+uy*uy);
  
  //momentum of e+ and e- in (ax,ay,az) axis. Note that az=0 by definition.
  //Double_t ptep = iep->Px()*ax + iep->Py()*ay;
  //Double_t ptem = iem->Px()*ax + iem->Py()*ay;
  
  Double_t pxep = px1;
  Double_t pyep = py1;
  Double_t pzep = pz1;
  Double_t pxem = px2;
  Double_t pyem = py2;
  Double_t pzem = pz2;
  
  //vector product of pep X pem
  Double_t vpx = pyep*pzem - pzep*pyem;
  Double_t vpy = pzep*pxem - pxep*pzem;
  Double_t vpz = pxep*pyem - pyep*pxem;
  Double_t vp = sqrt(vpx*vpx+vpy*vpy+vpz*vpz);
  //Double_t thev = acos(vpz/vp);
  
  //unit vector of pep X pem
  Double_t vx = vpx/vp;
  Double_t vy = vpy/vp;
  Double_t vz = vpz/vp;
  
  //The third axis defined by vector product (ux,uy,uz)X(vx,vy,vz)
  Double_t wx = uy*vz - uz*vy;
  Double_t wy = uz*vx - ux*vz;
  //Double_t wz = ux*vy - uy*vx;
  //Double_t wl = sqrt(wx*wx+wy*wy+wz*wz);
  // by construction, (wx,wy,wz) must be a unit vector.
  // measure angle between (wx,wy,wz) and (ax,ay,0). The angle between them
  // should be small if the pair is conversion.
  // this function then returns values close to pi!
  Double_t cosPhiV = wx*ax + wy*ay;
  Double_t phiv = TMath::ACos(cosPhiV);
  
  return phiv;
}
//________________________________________________________________________
void AliAnalysisTaskElectronEfficiency::Terminate(const Option_t *)
{
  // Draw result to the screen
  // Called once at the end of the query
  Printf(" Now running: Terminate()");
  
  //get output data and draw 'fHistPt'
  //  if (!GetOutputData(0)) return;
  //  TH1D *hist=(TH1D*)(((TObjArray*)GetOutputData(0))->FindObject("fHistPt"));
  //  if (hist) hist->Draw();
}
//________________________________________________________________________
AliAnalysisTaskElectronEfficiency::~AliAnalysisTaskElectronEfficiency()
{
  /// Destructor
  Printf(" Now running: ~Destructor");
  
  if (fSignalsMC) delete fSignalsMC;
  
  Printf("deleting TList");
  if (fOutputList) {
    fOutputList->Clear();
    delete fOutputList;
  }
  Printf("setting pointer to 0x0");
  fOutputList=0x0;
  
  if (fOutputListSupportHistos) {
    fOutputListSupportHistos->Clear();
    delete fOutputListSupportHistos;
  }
  fOutputListSupportHistos=0x0;
  
  // other objects (pointers) may be deleted like this:
  //  if (fTrackCutsITSSA)  delete fTrackCutsITSSA;
  //  fTrackCutsITSSA=0x0;
}
//________________________________________________________________________
void AliAnalysisTaskElectronEfficiency::CreateHistograms(TString names, Int_t cutInstance)
{
  Printf(" Now running: CreateHistograms()");
  
  TObjArray *arrNames=names.Tokenize(";");
  TString name=Form("%02d",cutInstance);
  if (cutInstance<arrNames->GetEntriesFast()){
    name=arrNames->At(cutInstance)->GetName();
  }
  //Printf("%i\t %f\t %i\t %f\t %i\t %f\t ",fNptBins,fPtBins[1],fNetaBins,fEtaBins[1],fNphiBins,fPhiBins[1]);
  if(fCalcEfficiencyGen){
    TH3D *hNreco_Ele = new TH3D(Form("Nreco_Ele_%s",name.Data()),"",fNptBins,fPtBins,fNetaBins,fEtaBins,fNphiBins,fPhiBins);
    hNreco_Ele->Sumw2();
    fvReco_Ele.push_back(hNreco_Ele);
    TH3D *hNreco_Pos = new TH3D(Form("Nreco_Pos_%s",name.Data()),"",fNptBins,fPtBins,fNetaBins,fEtaBins,fNphiBins,fPhiBins);
    hNreco_Pos->Sumw2();
    fvReco_Pos.push_back(hNreco_Pos);
    if(fCalcEfficiencyPoslabel){
      TH3D *hNreco_Ele_poslabel = new TH3D(Form("Nreco_Ele_poslabel_%s",name.Data()),"",fNptBins,fPtBins,fNetaBins,fEtaBins,fNphiBins,fPhiBins);
      hNreco_Ele_poslabel->Sumw2();
      fvReco_Ele_poslabel.push_back(hNreco_Ele_poslabel);
      TH3D *hNreco_Pos_poslabel = new TH3D(Form("Nreco_Pos_poslabel_%s",name.Data()),"",fNptBins,fPtBins,fNetaBins,fEtaBins,fNphiBins,fPhiBins);
      hNreco_Pos_poslabel->Sumw2();
      fvReco_Pos_poslabel.push_back(hNreco_Pos_poslabel);
    }
  }
  
  if(fCalcEfficiencyRec){
    TH3D *hNreco_Rec_Ele = new TH3D(Form("Nreco_Rec_Ele_%s",name.Data()),"",fNptBins,fPtBins,fNetaBins,fEtaBins,fNphiBins,fPhiBins);
    hNreco_Rec_Ele->Sumw2();
    fvReco_Rec_Ele.push_back(hNreco_Rec_Ele);
    TH3D *hNreco_Rec_Pos = new TH3D(Form("Nreco_Rec_Pos_%s",name.Data()),"",fNptBins,fPtBins,fNetaBins,fEtaBins,fNphiBins,fPhiBins);
    hNreco_Rec_Pos->Sumw2();
    fvReco_Rec_Pos.push_back(hNreco_Rec_Pos);
    if(fCalcEfficiencyPoslabel){
      TH3D *hNreco_Rec_Ele_poslabel = new TH3D(Form("Nreco_Rec_Ele_poslabel_%s",name.Data()),"",fNptBins,fPtBins,fNetaBins,fEtaBins,fNphiBins,fPhiBins);
      hNreco_Rec_Ele_poslabel->Sumw2();
      fvReco_Rec_Ele_poslabel.push_back(hNreco_Rec_Ele_poslabel);
      TH3D *hNreco_Rec_Pos_poslabel = new TH3D(Form("Nreco_Rec_Pos_poslabel_%s",name.Data()),"",fNptBins,fPtBins,fNetaBins,fEtaBins,fNphiBins,fPhiBins);
      hNreco_Rec_Pos_poslabel->Sumw2();
      fvReco_Rec_Pos_poslabel.push_back(hNreco_Rec_Pos_poslabel);
    }
  }
  
  // one needs the histogram 'hAllPionsForRej' for each cutInstance independently, because empty events may differ between the cutsets which run together.
  TH3D *hAllPionsForRej = new TH3D(Form("AllPionsForRej_%s",name.Data()),Form("AllPionsForRej_%s",name.Data()),fNptBins,fPtBins,fNetaBins,fEtaBins,fNphiBins,fPhiBins);
  TH3D *hNPionsRejByAllSigns = new TH3D(Form("NPionsRejByAllSigns_%s",name.Data()),Form("NPionsRejByAllSigns_%s",name.Data()),fNptBins,fPtBins,fNetaBins,fEtaBins,fNphiBins,fPhiBins);
  TH3D *hNPionsRejByUnlike = new TH3D(Form("NPionsRejByUnlike_%s",name.Data()),Form("NPionsRejByUnlike_%s",name.Data()),fNptBins,fPtBins,fNetaBins,fEtaBins,fNphiBins,fPhiBins);
  fvAllPionsForRej.push_back(hAllPionsForRej);
  fvPionsRejByAllSigns.push_back(hNPionsRejByAllSigns);
  fvPionsRejByUnlike.push_back(hNPionsRejByUnlike);
  
  
  if(fDoPairing){  
    if(fCalcEfficiencyGen){
      TH2D *hNrecoPairsResonances  = new TH2D(Form("NrecoPairsResonances_%s",      name.Data()),"",fNmeeBins,fMeeBins,fNpteeBins,fPteeBins);
      TH2D *hNrecoPairsDiffMothers = new TH2D(Form("NrecoPairsDifferentMothers_%s",name.Data()),"",fNmeeBins,fMeeBins,fNpteeBins,fPteeBins);
      TH2D *hNrecoPairsCharm       = new TH2D(Form("NrecoPairsCharm_%s",           name.Data()),"",fNmeeBins,fMeeBins,fNpteeBins,fPteeBins);
      TH2D *hNrecoPairsBeauty      = new TH2D(Form("NrecoPairsBeauty_%s",          name.Data()),"",fNmeeBins,fMeeBins,fNpteeBins,fPteeBins);
      TH2D *hNrecoPairsHF          = new TH2D(Form("NrecoPairsHF_%s",              name.Data()),"",fNmeeBins,fMeeBins,fNpteeBins,fPteeBins);
      hNrecoPairsResonances ->Sumw2();
      hNrecoPairsDiffMothers->Sumw2();
      hNrecoPairsCharm      ->Sumw2();
      hNrecoPairsBeauty     ->Sumw2();
      hNrecoPairsHF         ->Sumw2();
      fvRecoPairsResonances  .push_back(hNrecoPairsResonances);
      fvRecoPairsDiffMothers .push_back(hNrecoPairsDiffMothers);
      fvRecoPairsCharm       .push_back(hNrecoPairsCharm);
      fvRecoPairsBeauty      .push_back(hNrecoPairsBeauty);
      fvRecoPairsHF          .push_back(hNrecoPairsHF);
    }
    if(fCalcEfficiencyRec){
      TH2D *hNrecoPairsRecResonances  = new TH2D(Form("NrecoPairsRecResonances_%s",      name.Data()),"",fNmeeBins,fMeeBins,fNpteeBins,fPteeBins);
      TH2D *hNrecoPairsRecDiffMothers = new TH2D(Form("NrecoPairsRecDifferentMothers_%s",name.Data()),"",fNmeeBins,fMeeBins,fNpteeBins,fPteeBins);
      TH2D *hNrecoPairsRecCharm       = new TH2D(Form("NrecoPairsRecCharm_%s",           name.Data()),"",fNmeeBins,fMeeBins,fNpteeBins,fPteeBins);
      TH2D *hNrecoPairsRecBeauty      = new TH2D(Form("NrecoPairsRecBeauty_%s",          name.Data()),"",fNmeeBins,fMeeBins,fNpteeBins,fPteeBins);
      TH2D *hNrecoPairsRecHF          = new TH2D(Form("NrecoPairsRecHF_%s",              name.Data()),"",fNmeeBins,fMeeBins,fNpteeBins,fPteeBins);
      hNrecoPairsRecResonances ->Sumw2();
      hNrecoPairsRecDiffMothers->Sumw2();
      hNrecoPairsRecCharm      ->Sumw2();
      hNrecoPairsRecBeauty     ->Sumw2();
      hNrecoPairsRecHF         ->Sumw2();
      fvRecoPairsRecResonances  .push_back(hNrecoPairsRecResonances);
      fvRecoPairsRecDiffMothers .push_back(hNrecoPairsRecDiffMothers);
      fvRecoPairsRecCharm       .push_back(hNrecoPairsRecCharm);
      fvRecoPairsRecBeauty      .push_back(hNrecoPairsRecBeauty);
      fvRecoPairsRecHF          .push_back(hNrecoPairsRecHF);
    }
  }
  // be really careful if you need to implement this (see comments in UserExec):
  //  TH3D *hNreco_Pio = new TH3D(Form("Nreco_Pio_%s",name.Data()),Form("Nreco_Pio_%s",name.Data()),fNptBins,fPtBins,fNetaBins,fEtaBins,fNphiBins,fPhiBins);
  //  TH3D *hNreco_Kao = new TH3D(Form("Nreco_Kao_%s",name.Data()),Form("Nreco_Kao_%s",name.Data()),fNptBins,fPtBins,fNetaBins,fEtaBins,fNphiBins,fPhiBins);
  //  TH3D *hNreco_Pro = new TH3D(Form("Nreco_Pro_%s",name.Data()),Form("Nreco_Pro_%s",name.Data()),fNptBins,fPtBins,fNetaBins,fEtaBins,fNphiBins,fPhiBins);
  //  fvReco_Pio.push_back(hNreco_Pio);
  //  fvReco_Kao.push_back(hNreco_Kao);
  //  fvReco_Pro.push_back(hNreco_Pro);
}
//________________________________________________________________________
void AliAnalysisTaskElectronEfficiency::CreateHistoGen()
{
  Printf(" Now running: CreateHistoGen()");
  
  Printf("fNptBins=%i\t fPtBins[1]=%f\t fNetaBins=%i\t fEtaBins[1]=%f\t fNphiBins=%i\t fPhiBins[1]=%f\t ",fNptBins,fPtBins[1],fNetaBins,fEtaBins[1],fNphiBins,fPhiBins[1]);
  if(fCalcEfficiencyGen){
    fNgen_Ele       = new TH3D("Ngen_electrons","",fNptBins,fPtBins,fNetaBins,fEtaBins,fNphiBins,fPhiBins);
    fNgen_Pos       = new TH3D("Ngen_positrons","",fNptBins,fPtBins,fNetaBins,fEtaBins,fNphiBins,fPhiBins);
  }
  if(fCalcEfficiencyRec){
    fNgen1_Rec_Ele  = new TH3D("Ngen1_Rec_electrons","",fNptBins,fPtBins,fNetaBins,fEtaBins,fNphiBins,fPhiBins);
    fNgen1_Rec_Pos  = new TH3D("Ngen1_Rec_positrons","",fNptBins,fPtBins,fNetaBins,fEtaBins,fNphiBins,fPhiBins);
    fNgen2_Rec_Ele  = new TH3D("Ngen2_Rec_electrons","",fNptBins,fPtBins,fNetaBins,fEtaBins,fNphiBins,fPhiBins);
    fNgen2_Rec_Pos  = new TH3D("Ngen2_Rec_positrons","",fNptBins,fPtBins,fNetaBins,fEtaBins,fNphiBins,fPhiBins);
  }
  if(fDoPairing){
    Printf("fNmeeBins=%i\t fMeeBins[1]=%f\t fNpteeBins=%i\t fPteeBins[1]=%f\t ",fNmeeBins,fMeeBins[1],fNpteeBins,fPteeBins[1]);    
    if(fCalcEfficiencyGen){
      fNgenPairsResonances     = new TH2D("NgenPairsResonances",      "",fNmeeBins,fMeeBins,fNpteeBins,fPteeBins);
      fNgenPairsDiffMothers    = new TH2D("NgenPairsDifferentMothers","",fNmeeBins,fMeeBins,fNpteeBins,fPteeBins);
      fNgenPairsCharm          = new TH2D("NgenPairsCharm",           "",fNmeeBins,fMeeBins,fNpteeBins,fPteeBins);
      fNgenPairsBeauty         = new TH2D("NgenPairsBeauty",          "",fNmeeBins,fMeeBins,fNpteeBins,fPteeBins);
      fNgenPairsHF             = new TH2D("NgenPairsHF",              "",fNmeeBins,fMeeBins,fNpteeBins,fPteeBins);
    }
    if(fCalcEfficiencyRec){
      fNgenPairsRecResonances  = new TH2D("NgenPairsRecResonances",      "",fNmeeBins,fMeeBins,fNpteeBins,fPteeBins);
      fNgenPairsRecDiffMothers = new TH2D("NgenPairsRecDifferentMothers","",fNmeeBins,fMeeBins,fNpteeBins,fPteeBins);
      fNgenPairsRecCharm       = new TH2D("NgenPairsRecCharm",           "",fNmeeBins,fMeeBins,fNpteeBins,fPteeBins);
      fNgenPairsRecBeauty      = new TH2D("NgenPairsRecBeauty",          "",fNmeeBins,fMeeBins,fNpteeBins,fPteeBins);
      fNgenPairsRecHF          = new TH2D("NgenPairsRecHF",              "",fNmeeBins,fMeeBins,fNpteeBins,fPteeBins);
    }
  }
}
//________________________________________________________________________
void AliAnalysisTaskElectronEfficiency::CreateSupportHistos()
{
  Printf(" Now running: CreateSupportHistos()");
  
  fOutputListSupportHistos = new TList();
  fOutputListSupportHistos->SetName(GetName());
  fOutputListSupportHistos->SetOwner();
  
  // Event variables
  TH1D* hnEvents         = new TH1D("nEvents","Number of processed events after cuts;;N_{events}",1,0.,1.);//,AliDielectronVarManager::kNevents);
  TH1D* hCent            = new TH1D("centrality","N. events vs. centrality;centrality [%];N_{events}",100,0,100);//,AliDielectronVarManager::kCentrality);
  TH2D* hNTrksEvent_cent = new TH2D("hNTrksEvent_cent", "N. tracks vs. centrality (not used);centrality [%];N_{Tracks}", 50,0.,100., 50,-0.5,49.5);
  TH2D* hNEleEvent_cent  = new TH2D("hNEleEvent_cent", "N. selected electrons per event vs. centrality;centrality [%];N_{Tracks}", 50,0.,100., 50,-0.5,49.5);
  TH1D* hVertexZ         = new TH1D("hVertexZ","hVertexZ",300,-15.,15.);
  TH1D* hNvertexCtrb     = new TH1D("hNvertexCtrb","hNvertexCtrb",5000,-0.5,4999.5);
  fOutputListSupportHistos->AddAt(hnEvents, 0);
  fOutputListSupportHistos->AddAt(hCent, 1);
  fOutputListSupportHistos->AddAt(hNTrksEvent_cent, 2);
  fOutputListSupportHistos->AddAt(hNEleEvent_cent, 3);
  fOutputListSupportHistos->AddAt(hVertexZ, 4);
  fOutputListSupportHistos->AddAt(hNvertexCtrb, 5);
  
  TH2D* hPInVarMgr_PInStandAlone = new TH2D("PInVarMgr_PInStandAlone","GetTPCInnerParam()->P() vs GetInnerParam()->P() (VarMgr method); PIn global tracking (VarMgr method) (GeV/c); PIn TPC standalone (GeV/c)",
                                            160,0.,8.,160,0.,8.);
  fOutputListSupportHistos->AddAt(hPInVarMgr_PInStandAlone, 6);
  
  // Track variables
  TH1D* hPt      = new TH1D("Pt","Pt;Pt [GeV];#tracks",200,0,10.);//,AliDielectronVarManager::kPt);
  TH2D* hP_PIn   = new TH2D("P_PIn","TPC inner P vs P;P [GeV/c];PIn (pTPC) [GeV/c]",
                            160,0.,8.,160,0.,8.);//,AliDielectronVarManager::kP,AliDielectronVarManager::kPIn);
  TH2D* hP_Pgen  = new TH2D("P_Pgen","P gen vs P;P [GeV/c];p_{gen} (GeV/c)",
                            160,0.,8.,160,0.,8.);//,AliDielectronVarManager::kP,AliDielectronVarManager::kPIn);
  fOutputListSupportHistos->AddAt(hPt,     7);
  fOutputListSupportHistos->AddAt(hP_PIn,  8);
  fOutputListSupportHistos->AddAt(hP_Pgen, 9);
  
  
  // ITS
  TH2D* hITS_dEdx_P = new TH2D("ITS_dEdx_P","ITS dEdx;P [GeV/c];ITS signal (arb units)",
                               160,0.,8.,700,0.,700.);//.,AliDielectronVarManager::kP,AliDielectronVarManager::kITSsignal,makeLogx);
  TH2D* hITSnSigmaEle_P = new TH2D("ITSnSigmaEle_P","ITS number of sigmas Electrons;P [GeV/c];ITS number of sigmas Electrons",
                                   160,0.,8.,100,-5.,5.);//.,AliDielectronVarManager::kP,AliDielectronVarManager::kITSnSigmaEle,makeLogx);
  TH2D* hITSnSigmaPio_P = new TH2D("ITSnSigmaPio_P","ITS number of sigmas Pions;P [GeV/c];ITS number of sigmas Pions",
                                   160,0.,8.,200,-10.,10.);//.,AliDielectronVarManager::kP,AliDielectronVarManager::kITSnSigmaPio,makeLogx);
  TH2D* hITSnSigmaKao_P = new TH2D("ITSnSigmaKao_P","ITS number of sigmas Kaons;P [GeV/c];ITS number of sigmas Kaons",
                                   160,0.,8.,200,-10.,10.);//.,AliDielectronVarManager::kP,AliDielectronVarManager::kITSnSigmaKao,makeLogx);
  fOutputListSupportHistos->AddAt(hITS_dEdx_P, 10);
  fOutputListSupportHistos->AddAt(hITSnSigmaEle_P, 11);
  fOutputListSupportHistos->AddAt(hITSnSigmaPio_P, 12);
  fOutputListSupportHistos->AddAt(hITSnSigmaKao_P, 13);
  
  // TPC
  TH2D* hTPC_dEdx_P = new TH2D("TPC_dEdx_P","TPC dEdx;PIn (pTPC) [GeV/c];TPC signal (arb units)",
                               160,0.,8.,120,0.,120.);//.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCsignal,makeLogx);
  TH2D* hTPCnSigmaEle_P = new TH2D("TPCnSigmaEle_P","TPC number of sigmas Electrons;PIn (pTPC) [GeV/c];TPC number of sigmas Electrons",
                                   160,0.,8.,100,-5.,5.);//.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCnSigmaEle,makeLogx);
  TH3D* hTPCnSigmaEle_P_dEdx = new TH3D("TPCnSigmaEle_P_dEdx","TPC number of sigmas Electrons;PIn (pTPC) [GeV/c];TPC number of sigmas Electrons;TPC signal (arb units)",
                                        80,0.,4.,80,-4.,4.,50,50.,100.);//.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCnSigmaEle,AliDielectronVarManager::kTPCsignal,makeLogx);
  TH2D* hTPCnSigmaPio_P = new TH2D("TPCnSigmaPio_P","TPC number of sigmas Pions;PIn (pTPC) [GeV/c];TPC number of sigmas Pions",
                                   160,0.,8.,200,-10.,10.);//.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCnSigmaPio,makeLogx);
  TH2D* hTPCnSigmaKao_P = new TH2D("TPCnSigmaKao_P","TPC number of sigmas Kaons;PIn (pTPC) [GeV/c];TPC number of sigmas Kaons",
                                   160,0.,8.,200,-10.,10.);//.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCnSigmaKao,makeLogx);
  TH2D* hTPCnSigmaPro_P = new TH2D("TPCnSigmaPro_P","TPC number of sigmas Protons;PIn (pTPC) [GeV/c];TPC number of sigmas Protons",
                                   160,0.,8.,200,-10.,10.);//.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCnSigmaPro,makeLogx);
  fOutputListSupportHistos->AddAt(hTPC_dEdx_P, 14);
  fOutputListSupportHistos->AddAt(hTPCnSigmaEle_P, 15);
  fOutputListSupportHistos->AddAt(hTPCnSigmaEle_P_dEdx, 16);
  fOutputListSupportHistos->AddAt(hTPCnSigmaPio_P, 17);
  fOutputListSupportHistos->AddAt(hTPCnSigmaKao_P, 18);
  fOutputListSupportHistos->AddAt(hTPCnSigmaPro_P, 19);
  
  TH1D* hUnused20 = new TH1D("hUnused20","hUnused20", 1,0.,1.);
  TH1D* hUnused21 = new TH1D("hUnused21","hUnused21", 1,0.,1.);
  fOutputListSupportHistos->AddAt(hUnused20, 20);
  fOutputListSupportHistos->AddAt(hUnused21, 21);
  
  // run dependency
  // run string "fsRunBins" must be sorted in increasing order!
  TObjArray *objaRuns=fsRunBins.Tokenize(", ");
  const Int_t nRuns=objaRuns->GetEntries();
  if (nRuns<1) { fsRunBins = "-99"; cout << "warning: bins for run dependence not specified!" << endl; }
  fsRunBins.Append(Form(", %i", (Int_t) (atoi(objaRuns->At(nRuns-1)->GetName()) + 10))); // create upper limit for bin of last run!
  
  Double_t* dRunBinning = (AliDielectronHelper::MakeArbitraryBinning(fsRunBins))->GetMatrixArray();
  //for (int i=0; i<nRuns+1; i++) { cout << i << " \t " << dRunBinning[i] << " \t "; }
  
  Int_t nbinsPIn=80, nbinsTPCsig=50;
  TH3D* hTPC_dEdx_P_run = new TH3D("TPC_dEdx_P_run","TPC dEdx;TPC inner PIn (pTPC) [GeV/c];TPC signal (arb units);run number",
                                   nbinsPIn   , (AliDielectronHelper::MakeLinBinning(nbinsPIn   ,0.,4.))->GetMatrixArray(),
                                   nbinsTPCsig, (AliDielectronHelper::MakeLinBinning(nbinsTPCsig,50.,100.))->GetMatrixArray(),
                                   nRuns      , dRunBinning
                                   );//,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCsignal,AliDielectronVarManager::kRunNumber);
  TH3D* hTPCnSigmaEle_P_run = new TH3D("TPCnSigmaEle_P_run","TPC number of sigmas Electrons;TPC inner PIn (pTPC) [GeV/c];TPC number of sigmas Electrons;run number",
                                       nbinsPIn   , (AliDielectronHelper::MakeLinBinning(nbinsPIn   ,0.,4.))->GetMatrixArray(),
                                       50         , (AliDielectronHelper::MakeLinBinning(50         ,-5.,5.))->GetMatrixArray(),
                                       nRuns      , dRunBinning
                                       );//,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCnSigmaEle,AliDielectronVarManager::kRunNumber);
  TH3D* hITSnSigmaEle_P_run = new TH3D("ITSnSigmaEle_P_run","ITS number of sigmas Electrons;P [GeV/c];ITS number of sigmas Electrons;run number",
                                       nbinsPIn   , (AliDielectronHelper::MakeLinBinning(nbinsPIn   ,0.,4.))->GetMatrixArray(),
                                       50         , (AliDielectronHelper::MakeLinBinning(50         ,-5.,5.))->GetMatrixArray(),
                                       nRuns      , dRunBinning
                                       );//,AliDielectronVarManager::kP,AliDielectronVarManager::kITSnSigmaEle,AliDielectronVarManager::kRunNumber);
  fOutputListSupportHistos->AddAt(hTPC_dEdx_P_run, 22);
  fOutputListSupportHistos->AddAt(hTPCnSigmaEle_P_run, 23);
  fOutputListSupportHistos->AddAt(hITSnSigmaEle_P_run, 24);
  
  // TOF
  TH2D* hTOFnSigmaEle_P = new TH2D("TOFnSigmaEle_P","TOF number of sigmas Electrons;PIn (pTPC) [GeV/c];TOF number of sigmas Electrons",
                                   160,0.,8.,100,-5.,5.);//,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTOFnSigmaEle,makeLogx);
  TH2D* hTOFbeta = new TH2D("TOFbeta","TOF beta;PIn (pTPC) [GeV/c];TOF beta",
                            160,0.,8.,120,0.,1.2);//,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTOFbeta,makeLogx);
  fOutputListSupportHistos->AddAt(hTOFnSigmaEle_P, 25);
  fOutputListSupportHistos->AddAt(hTOFbeta, 26);
  
  // Eta and Phi
  TH1D* hEta = new TH1D("Eta","Eta; Eta;#tracks",
                        200,-2,2);//,AliDielectronVarManager::kEta);
  TH1D* hPhi = new TH1D("Phi","Phi; Phi;#tracks",
                        320,0.,6.4);//,AliDielectronVarManager::kPhi);
  TH2D* hEta_Phi = new TH2D("Eta_Phi","Eta Phi Map; Eta; Phi",
                            100,-1,1,320,0,6.4);//,AliDielectronVarManager::kEta,AliDielectronVarManager::kPhi);
  fOutputListSupportHistos->AddAt(hEta, 27);
  fOutputListSupportHistos->AddAt(hPhi, 28);
  fOutputListSupportHistos->AddAt(hEta_Phi, 29);
  
  TH2D* hTPC_dEdx_Eta = new TH2D("TPC_dEdx_Eta","TPC dEdx;Eta;TPC signal (arb units)",
                                 100,-1,1,120,0.,120.);//.,AliDielectronVarManager::kEta,AliDielectronVarManager::kTPCsignal);
  TH3D* hTPC_dEdx_Eta_P = new TH3D("TPC_dEdx_Eta_P","TPC dEdx;Eta;TPC signal (arb units); PIn (pTPC) [GeV/c]",
                                   100,-1,1,60,0.,120.,80,0.,4.);//.,AliDielectronVarManager::kEta,AliDielectronVarManager::kTPCsignal,AliDielectronVarManager::kPIn);
  TH2D* hTPCnSigmaEle_Eta = new TH2D("TPCnSigmaEle_Eta","TPC number of sigmas Electrons; Eta; TPC number of sigmas Electrons",
                                     100,-1,1,100,-5.,5.);//.,AliDielectronVarManager::kEta,AliDielectronVarManager::kTPCnSigmaEle);
  // 3D, may be used for dEdx eta correction
  TH3D* hTPCnSigmaEle_Eta_P = new TH3D("TPCnSigmaEle_Eta_P","TPC number of sigmas Electrons; Eta; TPC number of sigmas Electrons; PIn (pTPC) [GeV/c]",
                                       50,-1,1, 50,-5.,5., 80,0.,4.);//.,AliDielectronVarManager::kEta,AliDielectronVarManager::kTPCnSigmaEle,AliDielectronVarManager::kPIn);
  TH3D* hTPCnSigmaEle_Eta_RefMultTPConly = new TH3D("TPCnSigmaEle_Eta_RefMultTPConly","TPC Ref Mult from AOD header';Eta;n#sigma_{ele}^{TPC};N_{TPC ref}",
                                                    50,-1,1, 50,-5.,5., 100,0.,5000.);//,AliDielectronVarManager::kEta,AliDielectronVarManager::kTPCnSigmaEle,AliDielectronVarManager::kRefMultTPConly);
  TH3D* hTPCnSigmaEle_Eta_Nacc = new TH3D("TPCnSigmaEle_Eta_Nacc","Nacc from 'AliDielectronHelper::GetNacc()';Eta;n#sigma_{ele}^{TPC};N_{acc}",
                                          50,-1,1, 50,-5.,5., 100,0.,5000.);//,AliDielectronVarManager::kEta,AliDielectronVarManager::kTPCnSigmaEle,AliDielectronVarManager::kNacc);
  fOutputListSupportHistos->AddAt(hTPC_dEdx_Eta, 30);
  fOutputListSupportHistos->AddAt(hTPC_dEdx_Eta_P, 31);
  fOutputListSupportHistos->AddAt(hTPCnSigmaEle_Eta, 32);
  fOutputListSupportHistos->AddAt(hTPCnSigmaEle_Eta_P, 33);
  fOutputListSupportHistos->AddAt(hTPCnSigmaEle_Eta_RefMultTPConly, 34);
  fOutputListSupportHistos->AddAt(hTPCnSigmaEle_Eta_Nacc, 35);
  //  TH2D* hTPCnSigmaKao_Eta = new TH2D("TPCnSigmaKao_Eta","TPC number of sigmas Kaons; Eta; TPC number of sigmas Kaons",
  //                        100,-1,1,200,-10.,10.);//.,AliDielectronVarManager::kEta,AliDielectronVarManager::kTPCnSigmaKao);
  //  TH3D* hTPCnSigmaKao_Eta_P = new TH3D("TPCnSigmaKao_Eta_P","TPC number of sigmas Kaons; Eta; TPC number of sigmas Kaons; PIn (pTPC) [GeV/c]",
  //                        50,-1,1,100,-10.,10.,80,0.,4.);//.,AliDielectronVarManager::kEta,AliDielectronVarManager::kTPCnSigmaKao,AliDielectronVarManager::kPIn);
  //  fOutputListSupportHistos->AddAt(hTPCnSigmaKao_Eta, 34);
  //  fOutputListSupportHistos->AddAt(hTPCnSigmaKao_Eta_P, 35);
  
  // DCA
  TH1D* hdXY = new TH1D("dXY","dXY;dXY [cm];#tracks",
                        200,-2.,2.);//.,AliDielectronVarManager::kImpactParXY);
  TH1D* hdZ = new TH1D("dZ","dZ;dZ [cm];#tracks",
                       400,-4.,4.);//.,AliDielectronVarManager::kImpactParZ);
  TH2D* hdXY_dZ = new TH2D("dXY_dZ","dXY dZ Map; dXY; dZ",
                           100,-1.,1.,300,-3.,3.);//.,AliDielectronVarManager::kImpactParXY,AliDielectronVarManager::kImpactParZ);
  fOutputListSupportHistos->AddAt(hdXY, 36);
  fOutputListSupportHistos->AddAt(hdZ, 37);
  fOutputListSupportHistos->AddAt(hdXY_dZ, 38);
  
  // Quality
  TH1D* hTPCcrossedRowsOverFindable = new TH1D("TPCcrossedRowsOverFindable","Number of Crossed Rows TPC over Findable;TPC crossed rows over findable;#tracks",120,0.,1.2);//,AliDielectronVarManager::kNFclsTPCfCross);
  TH1D* hTPCcrossedRows = new TH1D("TPCcrossedRows","Number of Crossed Rows TPC;TPC crossed rows;#tracks",160,-0.5,159.5);//,AliDielectronVarManager::kNFclsTPCr);
  TH1D* hTPCnCls = new TH1D("TPCnCls","Number of Clusters TPC;TPC number clusters;#tracks",160,-0.5,159.5);//,AliDielectronVarManager::kNclsTPC);
  TH1D* hITSnCls = new TH1D("ITSnCls","Number of Clusters ITS;ITS number clusters;#tracks",10,-0.5,9.5);//,AliDielectronVarManager::kNclsITS);
  TH1D* hTPCchi2 = new TH1D("TPCchi2","TPC chi2 per Cluster;TPC chi2/Cl;#tracks",100,0.,10.);//,AliDielectronVarManager::kTPCchi2Cl);
  TH1D* hITSchi2 = new TH1D("ITSchi2","ITS chi2 per Cluster;ITS chi2/Cl;#tracks",100,0.,10.);//,AliDielectronVarManager::kITSchi2Cl);
  TH1D* hTPCnClsS = new TH1D("TPCnClsS",";TPC number of shared clusters;#tracks",160,-0.5,159.5);//,AliDielectronVarManager::kNclsSTPC);
  TH1D* hITSnClsS = new TH1D("ITSnClsS",";ITS number of shared clusters;#tracks",10,-0.5,9.5);//,AliDielectronVarManager::kNclsSITS);
  TH1D* hNclsSFracTPC = new TH1D("NclsSFracTPC","Fraction of shared clusters assigned in the TPC;TPC fraction of shared clusters;#tracks",120,0,1.2);//.,AliDielectronVarManager::kNclsSFracTPC);
  TH1D* hNclsSFracITS = new TH1D("NclsSFracITS","Fraction of shared clusters assigned in the ITS;ITS fraction of shared clusters;#tracks",120,0,1.2);//.,AliDielectronVarManager::kNclsSFracITS);
  TH1D* hTPCclsDiff = new TH1D("TPCclsDiff","TPC cluster difference;N_{d#it{E}/d#it{x} points}^{TPC} - N_{cls}^{TPC};#tracks",100,-80,20);//.,AliDielectronVarManager::kTPCclsDiff);
  TH1D* hTPCsignalN = new TH1D("TPCsignalN","Number of PID Clusters TPC;N_{d#it{E}/d#it{x} points}^{TPC};#tracks",160,-0.5,159.5);//.,AliDielectronVarManager::kTPCsignalN); //kNclsTPCdEdx
  fOutputListSupportHistos->AddAt(hTPCcrossedRowsOverFindable, 39);
  fOutputListSupportHistos->AddAt(hTPCcrossedRows, 40);
  fOutputListSupportHistos->AddAt(hTPCnCls, 41);
  fOutputListSupportHistos->AddAt(hITSnCls, 42);
  fOutputListSupportHistos->AddAt(hTPCchi2, 43);
  fOutputListSupportHistos->AddAt(hITSchi2, 44);
  fOutputListSupportHistos->AddAt(hTPCnClsS, 45);
  fOutputListSupportHistos->AddAt(hITSnClsS, 46);
  fOutputListSupportHistos->AddAt(hNclsSFracTPC, 47);
  fOutputListSupportHistos->AddAt(hNclsSFracITS, 48);
  fOutputListSupportHistos->AddAt(hTPCclsDiff, 49);
  fOutputListSupportHistos->AddAt(hTPCsignalN, 50);
  TH2D* hTPCcrossedRows_TPCnCls = new TH2D("TPCcrossedRows_TPCnCls","TPC crossed rows vs TPC number clusters;TPC number clusters;TPC crossed rows",
                                           160,-0.5,159.5,160,-0.5,159.5);//,AliDielectronVarManager::kNclsTPC,AliDielectronVarManager::kNFclsTPCr);
  TH2D* hTPCcrossedRows_Pt = new TH2D("TPCcrossedRows_Pt","TPC crossed rows vs Pt;Pt [GeV];TPC crossed rows",
                                      160,0.,8.,160,-0.5,159.5);//,AliDielectronVarManager::kPt,AliDielectronVarManager::kNFclsTPCr);
  fOutputListSupportHistos->AddAt(hTPCcrossedRows_TPCnCls, 51);
  fOutputListSupportHistos->AddAt(hTPCcrossedRows_Pt, 52);
  
  
  
  // pair-prefilter histograms
  //2D and 3D histograms
  TH2D* hMeePtee    = new TH2D("InvMass_PairPt",";Inv. Mass [GeV];Pair Pt [GeV];#pairs",
                               100,0.,1., 100,0.,2.);//AliDielectronVarManager::kM, AliDielectronVarManager::kPt);
  TH3D* hMeePhiVOpen= new TH3D("InvMass_PhivPair_OpeningAngle",";Inv. Mass [GeV];PhiV;Opening Angle",
                               100,0.,0.5, 100,0.,TMath::Pi(), 50,0.,TMath::Pi()/2.);//AliDielectronVarManager::kM, AliDielectronVarManager::kPhivPair, AliDielectronVarManager::kOpeningAngle);
  //opening angle and PhiV
  TH2D* hMeeOpen    = new TH2D("InvMass_OpeningAngle",";Inv. Mass [GeV];Opening Angle;#pairs",
                               100,0.,1., 100,0.,TMath::Pi());//AliDielectronVarManager::kM, AliDielectronVarManager::kOpeningAngle);
  TH2D* hMeePhiV    = new TH2D("InvMass_PhivPair",";Inv. Mass [GeV];PhiV;#pairs",
                               100,0.,1., 100,0.,TMath::Pi());//AliDielectronVarManager::kM, AliDielectronVarManager::kPhivPair);
  TH2D* hPteeOpen   = new TH2D("PairPt_OpeningAngle",";Pair Pt [GeV];Opening Angle;#pairs",
                               100,0.,2., 100,0.,TMath::Pi());//AliDielectronVarManager::kPt, AliDielectronVarManager::kOpeningAngle);
  TH2D* hPteePhiV   = new TH2D("PairPt_PhivPair",";Pair Pt [GeV];PhiV;#pairs",
                               100,0.,2., 100,0.,TMath::Pi());//AliDielectronVarManager::kPt, AliDielectronVarManager::kPhivPair);
  TH2D* hOpenPhiV   = new TH2D("OpeningAngle_PhivPair",";Opening Angle;PhiV;#pairs",
                               100,0.,TMath::Pi(), 100,0.,TMath::Pi());//AliDielectronVarManager::kOpeningAngle, AliDielectronVarManager::kPhivPair);
  fOutputListSupportHistos->AddAt(hMeePtee,     53);
  fOutputListSupportHistos->AddAt(hMeePhiVOpen, 54);
  fOutputListSupportHistos->AddAt(hMeeOpen,     55);
  fOutputListSupportHistos->AddAt(hMeePhiV,     56);
  fOutputListSupportHistos->AddAt(hPteeOpen,    57);
  fOutputListSupportHistos->AddAt(hPteePhiV,    58);
  fOutputListSupportHistos->AddAt(hOpenPhiV,    59);
  
  // histograms for ITS eta and dEdx correction
  TH2D* hITSnSigmaEle_Eta = new TH2D("ITSnSigmaEle_Eta","ITS number of sigmas Electrons; Eta; ITS number of sigmas Electrons",
                                     100,-1,1,100,-5.,5.);//.,AliDielectronVarManager::kEta,AliDielectronVarManager::kITSnSigmaEle);
  // 3D, may be used for dEdx eta correction
  TH3D* hITSnSigmaEle_Eta_P = new TH3D("ITSnSigmaEle_Eta_P","ITS number of sigmas Electrons; Eta; ITS number of sigmas Electrons; P [GeV/c]",
                                       50,-1,1, 50,-5.,5., 80,0.,4.);//.,AliDielectronVarManager::kEta,AliDielectronVarManager::kITSnSigmaEle,AliDielectronVarManager::kP);
  TH3D* hITSnSigmaEle_Eta_RefMultTPConly = new TH3D("ITSnSigmaEle_Eta_RefMultTPConly","TPC Ref Mult from AOD header';Eta;n#sigma_{ele}^{ITS};N_{TPC ref}",
                                                    50,-1,1, 50,-5.,5., 100,0.,5000.);//,AliDielectronVarManager::kEta,AliDielectronVarManager::kITSnSigmaEle,AliDielectronVarManager::kRefMultTPConly);
  TH3D* hITSnSigmaEle_Eta_Nacc = new TH3D("ITSnSigmaEle_Eta_Nacc","Nacc from 'AliDielectronHelper::GetNacc()';Eta;n#sigma_{ele}^{ITS};N_{acc}",
                                          50,-1,1, 50,-5.,5., 100,0.,5000.);//,AliDielectronVarManager::kEta,AliDielectronVarManager::kITSnSigmaEle,AliDielectronVarManager::kNacc);
  fOutputListSupportHistos->AddAt(hITSnSigmaEle_Eta,      60);
  fOutputListSupportHistos->AddAt(hITSnSigmaEle_Eta_P,    61);
  fOutputListSupportHistos->AddAt(hITSnSigmaEle_Eta_RefMultTPConly, 62);
  fOutputListSupportHistos->AddAt(hITSnSigmaEle_Eta_Nacc, 63);
  
  // raw sigma Ele in ITS and TPC
  TH2D* hITSnSigmaEleRaw_P = new TH2D("ITSnSigmaEleRaw_P","ITS number of sigmas Electrons (raw);P [GeV/c];ITS number of sigmas Electrons (raw)",
                                      160,0.,8.,100,-5.,5.);//.,AliDielectronVarManager::kP,AliDielectronVarManager::kITSnSigmaEleRaw);
  TH2D* hTPCnSigmaEleRaw_P = new TH2D("TPCnSigmaEleRaw_P","TPC number of sigmas Electrons (raw);PIn (pTPC) [GeV/c];TPC number of sigmas Electrons (raw)",
                                      160,0.,8.,100,-5.,5.);//.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCnSigmaEleRaw);
  TH2D* hITSnSigmaEleRaw_Eta = new TH2D("ITSnSigmaEleRaw_Eta","ITS number of sigmas Electrons (raw); Eta; ITS number of sigmas Electrons (raw)",
                                        100,-1,1,100,-5.,5.);//.,AliDielectronVarManager::kEta,AliDielectronVarManager::kITSnSigmaEleRaw);
  TH2D* hTPCnSigmaEleRaw_Eta = new TH2D("TPCnSigmaEleRaw_Eta","TPC number of sigmas Electrons (raw); Eta; TPC number of sigmas Electrons (raw)",
                                        100,-1,1,100,-5.,5.);//.,AliDielectronVarManager::kEta,AliDielectronVarManager::kTPCnSigmaEleRaw);
  fOutputListSupportHistos->AddAt(hITSnSigmaEleRaw_P,     64);
  fOutputListSupportHistos->AddAt(hTPCnSigmaEleRaw_P,     65);
  fOutputListSupportHistos->AddAt(hITSnSigmaEleRaw_Eta,   66);
  fOutputListSupportHistos->AddAt(hTPCnSigmaEleRaw_Eta,   67);
  
  TH1D* hUnused68 = new TH1D("hUnused","hUnused",1,0,1);
  fOutputListSupportHistos->AddAt(hUnused68, 68);
  TH1D* hUnused69 = new TH1D("hUnused","hUnused",1,0,1);
  fOutputListSupportHistos->AddAt(hUnused69, 69);
  
  
  TH1D* hPdgCode   = new TH1D("hPdgCode",  "hPdgCode",   GetPDGcodes()->GetNrows()-1, GetPDGcodes()->GetMatrixArray());
  TH1D* hPdgCodeM  = new TH1D("hPdgCodeM", "hPdgCodeM",  GetPDGcodes()->GetNrows()-1, GetPDGcodes()->GetMatrixArray());
  TH1D* hPdgCodeGM = new TH1D("hPdgCodeGM","hPdgCodeGM", GetPDGcodes()->GetNrows()-1, GetPDGcodes()->GetMatrixArray());
  TH2D* hPdgCodeM_GM = new TH2D("hPdgCodeM_GM","hPdgCodeM_GM", GetPDGcodes()->GetNrows()-1, GetPDGcodes()->GetMatrixArray(),
                                GetPDGcodes()->GetNrows()-1, GetPDGcodes()->GetMatrixArray());//AliDielectronVarManager::kPdgCodeMother, AliDielectronVarManager::kPdgCodeGrandMother);
  fOutputListSupportHistos->AddAt(hPdgCode,     70);
  fOutputListSupportHistos->AddAt(hPdgCodeM,    71);
  fOutputListSupportHistos->AddAt(hPdgCodeGM,   72);
  fOutputListSupportHistos->AddAt(hPdgCodeM_GM, 73);
  
}
//__________________________________________________________________
void AliAnalysisTaskElectronEfficiency::AddSignalMC(AliDielectronSignalMC* signal)
{
  //
  //  Add an MC signal to the signals list
  //
  if(!fSignalsMC) {
    fSignalsMC = new TObjArray();
    fSignalsMC->SetOwner();
  }
  fSignalsMC->Add(signal);
  cout << "added MCsignal: " << fSignalsMC->At(fSignalsMC->GetEntriesFast()-1)->GetName() << endl;
  if (fSignalsMC->GetEntriesFast()>1) {
    AliWarning("WARNING: only the first MC signal will be used for the electron selection, others are ignored.");
  }
}
//______________________________________________
void AliAnalysisTaskElectronEfficiency::SetCentroidCorrFunction(TObject *fun, UInt_t varx, UInt_t vary, UInt_t varz)
{
  UInt_t valType[20] = {0};
  valType[0]=varx;     valType[1]=vary;     valType[2]=varz;
  TString key = Form("cntrdTPC%d%d%d",varx,vary,varz);
  // StoreVariables() sets the uniqueIDs of the axes to match with the VarManager enum and sets axis titles accordingly.
  // Needed so that AliDielectronPID can extract the needed variables into its local TBits 'fUsedVars'.
  // clone temporary histogram, otherwise it will not be streamed to file!
  if (fun->InheritsFrom(TH1::Class())) {
    AliDielectronHistos::StoreVariables(fun, valType);
    fPostPIDCntrdCorrTPC = (TH1*)fun->Clone(key.Data());
  }
  else if (fun->InheritsFrom(TF1::Class())) {
    AliDielectronHistos::StoreVariables(static_cast<TF1*>(fun)->GetHistogram(), valType);
    fPostPIDCntrdCorrTPC = (TH1*) static_cast<TF1*>(fun)->GetHistogram()->Clone(key.Data());
    if(fPostPIDCntrdCorrTPC) {
      fPostPIDCntrdCorrTPC->GetListOfFunctions()->AddAt(fun,0);
    }
  }
  else {
    AliWarning(Form("WARNING: PID correction object has invalid type: %s!", fun->IsA()->GetName())); return;
  }
  if(fPostPIDCntrdCorrTPC) {
    // check for corrections and add their variables to the fill map
    printf("POST TPC PID CORRECTION added for centroids:  ");
    switch(fPostPIDCntrdCorrTPC->GetDimension()) {
      case 3: printf(" %s, ",fPostPIDCntrdCorrTPC->GetZaxis()->GetName());
      case 2: printf(" %s, ",fPostPIDCntrdCorrTPC->GetYaxis()->GetName());
      case 1: printf(" %s ",fPostPIDCntrdCorrTPC->GetXaxis()->GetName());
    }
    printf("\n");
    fUsedVars->SetBitNumber(varx, kTRUE); // probably not needed within efficiency task...
    fUsedVars->SetBitNumber(vary, kTRUE);
    fUsedVars->SetBitNumber(varz, kTRUE);
  }
}
//______________________________________________
void AliAnalysisTaskElectronEfficiency::SetWidthCorrFunction(TObject *fun, UInt_t varx, UInt_t vary, UInt_t varz)
{
  UInt_t valType[20] = {0};
  valType[0]=varx;     valType[1]=vary;     valType[2]=varz;
  TString key = Form("wdthTPC%d%d%d",varx,vary,varz);
  // StoreVariables() sets the uniqueIDs of the axes to match with the VarManager enum and sets axis titles accordingly.
  // Needed so that AliDielectronPID can extract the needed variables into its local TBits 'fUsedVars'.
  // clone temporary histogram, otherwise it will not be streamed to file!
  if (fun->InheritsFrom(TH1::Class())) {
    AliDielectronHistos::StoreVariables(fun, valType);
    fPostPIDWdthCorrTPC = (TH1*)fun->Clone(key.Data());
  }
  else if (fun->InheritsFrom(TF1::Class())) {
    AliDielectronHistos::StoreVariables(static_cast<TF1*>(fun)->GetHistogram(), valType);
    fPostPIDWdthCorrTPC = (TH1*) static_cast<TF1*>(fun)->GetHistogram()->Clone(key.Data());
    if(fPostPIDWdthCorrTPC) {
      fPostPIDWdthCorrTPC->GetListOfFunctions()->AddAt(fun,0);
    }
  }
  else {
    AliWarning(Form("WARNING: PID correction object has invalid type: %s!", fun->IsA()->GetName())); return;
  }
  if(fPostPIDWdthCorrTPC)  {
    // check for corrections and add their variables to the fill map
    printf("POST TPC PID CORRECTION added for widths:  ");
    switch(fPostPIDWdthCorrTPC->GetDimension()) {
      case 3: printf(" %s, ",fPostPIDWdthCorrTPC->GetZaxis()->GetName());
      case 2: printf(" %s, ",fPostPIDWdthCorrTPC->GetYaxis()->GetName());
      case 1: printf(" %s ",fPostPIDWdthCorrTPC->GetXaxis()->GetName());
    }
    printf("\n");
    fUsedVars->SetBitNumber(varx, kTRUE); // probably not needed within efficiency task...
    fUsedVars->SetBitNumber(vary, kTRUE);
    fUsedVars->SetBitNumber(varz, kTRUE);
  }
}
//______________________________________________
void AliAnalysisTaskElectronEfficiency::SetCentroidCorrFunctionITS(TObject *fun, UInt_t varx, UInt_t vary, UInt_t varz)
{
  UInt_t valType[20] = {0};
  valType[0]=varx;     valType[1]=vary;     valType[2]=varz;
  TString key = Form("cntrdITS%d%d%d",varx,vary,varz);
  // StoreVariables() sets the uniqueIDs of the axes to match with the VarManager enum and sets axis titles accordingly.
  // Needed so that AliDielectronPID can extract the needed variables into its local TBits 'fUsedVars'.
  // clone temporary histogram, otherwise it will not be streamed to file!
  if (fun->InheritsFrom(TH1::Class())) {
    AliDielectronHistos::StoreVariables(fun, valType);
    fPostPIDCntrdCorrITS = (TH1*)fun->Clone(key.Data());
  }
  else if (fun->InheritsFrom(TF1::Class())) {
    AliDielectronHistos::StoreVariables(static_cast<TF1*>(fun)->GetHistogram(), valType);
    fPostPIDCntrdCorrITS = (TH1*) static_cast<TF1*>(fun)->GetHistogram()->Clone(key.Data());
    if(fPostPIDCntrdCorrITS) {
      fPostPIDCntrdCorrITS->GetListOfFunctions()->AddAt(fun,0);
    }
  }
  else {
    AliWarning(Form("WARNING: PID correction object has invalid type: %s!", fun->IsA()->GetName())); return;
  }
  if(fPostPIDCntrdCorrITS)  {
    // check for corrections and add their variables to the fill map
    printf("POST ITS PID CORRECTION added for centroids:  ");
    switch(fPostPIDCntrdCorrITS->GetDimension()) {
      case 3: printf(" %s, ",fPostPIDCntrdCorrITS->GetZaxis()->GetName());
      case 2: printf(" %s, ",fPostPIDCntrdCorrITS->GetYaxis()->GetName());
      case 1: printf(" %s ",fPostPIDCntrdCorrITS->GetXaxis()->GetName());
    }
    printf("\n");
    fUsedVars->SetBitNumber(varx, kTRUE); // probably not needed within efficiency task...
    fUsedVars->SetBitNumber(vary, kTRUE);
    fUsedVars->SetBitNumber(varz, kTRUE);
  }
}
//______________________________________________
void AliAnalysisTaskElectronEfficiency::SetWidthCorrFunctionITS(TObject *fun, UInt_t varx, UInt_t vary, UInt_t varz)
{
  UInt_t valType[20] = {0};
  valType[0]=varx;     valType[1]=vary;     valType[2]=varz;
  TString key = Form("wdthITS%d%d%d",varx,vary,varz);
  // StoreVariables() sets the uniqueIDs of the axes to match with the VarManager enum and sets axis titles accordingly.
  // Needed so that AliDielectronPID can extract the needed variables into its local TBits 'fUsedVars'.
  // clone temporary histogram, otherwise it will not be streamed to file!
  if (fun->InheritsFrom(TH1::Class())) {
    AliDielectronHistos::StoreVariables(fun, valType);
    fPostPIDWdthCorrITS = (TH1*)fun->Clone(key.Data());
  }
  else if (fun->InheritsFrom(TF1::Class())) {
    AliDielectronHistos::StoreVariables(static_cast<TF1*>(fun)->GetHistogram(), valType);
    fPostPIDWdthCorrITS = (TH1*) static_cast<TF1*>(fun)->GetHistogram()->Clone(key.Data());
    if(fPostPIDWdthCorrITS) {
      fPostPIDWdthCorrITS->GetListOfFunctions()->AddAt(fun,0);
    }
  }
  else {
    AliWarning(Form("WARNING: PID correction object has invalid type: %s!", fun->IsA()->GetName())); return;
  }
  if(fPostPIDWdthCorrITS)  {
    // check for corrections and add their variables to the fill map
    printf("POST ITS PID CORRECTION added for widths:  ");
    switch(fPostPIDWdthCorrITS->GetDimension()) {
      case 3: printf(" %s, ",fPostPIDWdthCorrITS->GetZaxis()->GetName());
      case 2: printf(" %s, ",fPostPIDWdthCorrITS->GetYaxis()->GetName());
      case 1: printf(" %s ",fPostPIDWdthCorrITS->GetXaxis()->GetName());
    }
    printf("\n");
    fUsedVars->SetBitNumber(varx, kTRUE); // probably not needed within efficiency task...
    fUsedVars->SetBitNumber(vary, kTRUE);
    fUsedVars->SetBitNumber(varz, kTRUE);
  }
}
//______________________________________________
void AliAnalysisTaskElectronEfficiency::SetCentroidCorrFunctionTOF(TObject *fun, UInt_t varx, UInt_t vary, UInt_t varz)
{
  UInt_t valType[20] = {0};
  valType[0]=varx;     valType[1]=vary;     valType[2]=varz;
  TString key = Form("cntrdTOF%d%d%d",varx,vary,varz);
  // StoreVariables() sets the uniqueIDs of the axes to match with the VarManager enum and sets axis titles accordingly.
  // Needed so that AliDielectronPID can extract the needed variables into its local TBits 'fUsedVars'.
  // clone temporary histogram, otherwise it will not be streamed to file!
  if (fun->InheritsFrom(TH1::Class())) {
    AliDielectronHistos::StoreVariables(fun, valType);
    fPostPIDCntrdCorrTOF = (TH1*)fun->Clone(key.Data());
  }
  else if (fun->InheritsFrom(TF1::Class())) {
    AliDielectronHistos::StoreVariables(static_cast<TF1*>(fun)->GetHistogram(), valType);
    fPostPIDCntrdCorrTOF = (TH1*) static_cast<TF1*>(fun)->GetHistogram()->Clone(key.Data());
    if(fPostPIDCntrdCorrTOF) {
      fPostPIDCntrdCorrTOF->GetListOfFunctions()->AddAt(fun,0);
    }
  }
  else {
    AliWarning(Form("WARNING: PID correction object has invalid type: %s!", fun->IsA()->GetName())); return;
  }
  if(fPostPIDCntrdCorrTOF)  {
    // check for corrections and add their variables to the fill map
    printf("POST TOF PID CORRECTION added for centroids:  ");
    switch(fPostPIDCntrdCorrTOF->GetDimension()) {
      case 3: printf(" %s, ",fPostPIDCntrdCorrTOF->GetZaxis()->GetName());
      case 2: printf(" %s, ",fPostPIDCntrdCorrTOF->GetYaxis()->GetName());
      case 1: printf(" %s ",fPostPIDCntrdCorrTOF->GetXaxis()->GetName());
    }
    printf("\n");
    fUsedVars->SetBitNumber(varx, kTRUE); // probably not needed within efficiency task...
    fUsedVars->SetBitNumber(vary, kTRUE);
    fUsedVars->SetBitNumber(varz, kTRUE);
  }
}
//______________________________________________
void AliAnalysisTaskElectronEfficiency::SetWidthCorrFunctionTOF(TObject *fun, UInt_t varx, UInt_t vary, UInt_t varz)
{
  UInt_t valType[20] = {0};
  valType[0]=varx;     valType[1]=vary;     valType[2]=varz;
  TString key = Form("wdthTOF%d%d%d",varx,vary,varz);
  // StoreVariables() sets the uniqueIDs of the axes to match with the VarManager enum and sets axis titles accordingly.
  // Needed so that AliDielectronPID can extract the needed variables into its local TBits 'fUsedVars'.
  // clone temporary histogram, otherwise it will not be streamed to file!
  if (fun->InheritsFrom(TH1::Class())) {
    AliDielectronHistos::StoreVariables(fun, valType);
    fPostPIDWdthCorrTOF = (TH1*)fun->Clone(key.Data());
  }
  else if (fun->InheritsFrom(TF1::Class())) {
    AliDielectronHistos::StoreVariables(static_cast<TF1*>(fun)->GetHistogram(), valType);
    fPostPIDWdthCorrTOF = (TH1*) static_cast<TF1*>(fun)->GetHistogram()->Clone(key.Data());
    if(fPostPIDWdthCorrTOF) {
      fPostPIDWdthCorrTOF->GetListOfFunctions()->AddAt(fun,0);
    }
  }
  else {
    AliWarning(Form("WARNING: PID correction object has invalid type: %s!", fun->IsA()->GetName())); return;
  }
  if(fPostPIDWdthCorrTOF)  {
    // check for corrections and add their variables to the fill map
    printf("POST TOF PID CORRECTION added for widths:  ");
    switch(fPostPIDWdthCorrTOF->GetDimension()) {
      case 3: printf(" %s, ",fPostPIDWdthCorrTOF->GetZaxis()->GetName());
      case 2: printf(" %s, ",fPostPIDWdthCorrTOF->GetYaxis()->GetName());
      case 1: printf(" %s ",fPostPIDWdthCorrTOF->GetXaxis()->GetName());
    }
    printf("\n");
    fUsedVars->SetBitNumber(varx, kTRUE); // probably not needed within efficiency task...
    fUsedVars->SetBitNumber(vary, kTRUE);
    fUsedVars->SetBitNumber(varz, kTRUE);
  }
}
//______________________________________________________________________________________
TVectorD *AliAnalysisTaskElectronEfficiency::GetPDGcodes()
{
  //
  // array of pdgcodes stored in TDatabasePDG
  //
  printf("AliAnalysisTaskElectronEfficiency::GetPDGcodes()\n");
  //TDatabasePDG *pdg = new TDatabasePDG(); // get warning: W-TDatabasePDG::TDatabasePDG: object already instantiated
  TDatabasePDG *pdg = TDatabasePDG::Instance();
  if (!pdg->ParticleList()) pdg->ReadPDGTable(); // it's not good to call it everytime (tries to redefine the table, many warnings).
  TGraph *gr = new TGraph();
  TIter next(pdg->ParticleList());
  TParticlePDG *p;
  Int_t i=0;
  while ((p = (TParticlePDG *)next())) {
    if(TMath::Abs(p->PdgCode()) < 1e+6) {
      //      printf("%s -> %d \n",p->GetName(),p->PdgCode());
      gr->SetPoint(i++, p->PdgCode(),1.);
    }
  }
  gr->Sort();
  TVectorD *vec = new TVectorD(gr->GetN(), gr->GetX());
  //  vec->Print();
  //delete pdg; // Was used in case of using "new TDatabasePDG()", but this is a singleton...
  // Gives a seg fault in AliPID::Init(). Probably due to "#define M(PID) TDatabasePDG::Instance()->GetParticle(fgkParticleCode[(PID)])->Mass()"
  delete gr;
  return vec;
}
//______________________________________________________________________________________
Double_t AliAnalysisTaskElectronEfficiency::GetSmearing(TObjArray *arr, Double_t x)
{
  TH1D *hisSlice(0x0);
  TH2D *hDeltaXvsXgen = static_cast<TH2D*> (arr->At(0));
  Int_t histIndex = TMath::Min( hDeltaXvsXgen->GetXaxis()->FindBin(x), arr->GetLast() );
  if (histIndex<1) histIndex=1;
  hisSlice = static_cast<TH1D*> (arr->At(histIndex));
  // get smear parameter via random selection from the x slices retreived from the deltax plot
  Double_t smearing(0.);
  if(hisSlice) smearing = hisSlice->GetRandom();
  
  return smearing;
}
