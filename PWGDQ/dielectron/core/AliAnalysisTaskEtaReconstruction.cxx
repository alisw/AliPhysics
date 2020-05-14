
 /*************************************************************************
 * Copyright(c) 1998-2018, ALICE Experiment at CERN, All rights reserved. *
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
//       Eta Reconstruction Task via the Dalitz decay channel            //
//                                        (description in .h file)       //
//                                                                       //
///////////////////////////////////////////////////////////////////////////

#include "AliAnalysisTaskEtaReconstruction.h"
#include "AliVTrack.h"
#include "AliAODInputHandler.h"

#include "AliAODEvent.h"
#include "AliESDEvent.h"
#include "AliMCEvent.h"
#include "AliAnalysisFilter.h"


#include "AliDielectronMC.h"
#include "AliDielectronVarManager.h"
#include "AliDielectronSignalMC.h"
#include "AliDielectronPair.h"
#include "AliDielectronHistos.h"

#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TTree.h"
#include "THnSparse.h"
#include "TLorentzVector.h"

#include "TChain.h"
#include "TSystem.h"

#include <iostream>

AliAnalysisTaskEtaReconstruction::~AliAnalysisTaskEtaReconstruction(){
  delete fPtPion;
  delete fPtEta;
  delete fPtEtaPrime;
  delete fPtRho;
  delete fPtOmega;
  delete fPtPhi;
  delete fPtJPsi;
  delete fCocktailFile;
  delete fOutputList;
  delete fSingleElectronList;
  delete fGeneratedPrimaryList;
  delete fGeneratedSecondaryList;
  delete fGeneratedSmearedPrimaryList;
  delete fGeneratedSmearedSecondaryList;
  delete fRecPrimaryList;
  delete fRecSecondaryList;
  delete fGeneratedPrimaryPairsList;
  delete fGeneratedSecondaryPairsList;
  delete fGeneratedSmearedPrimaryPairsList;
  delete fGeneratedSmearedSecondaryPairsList;
  delete fPairList;
  delete fFourPairList;
  delete fResolutionList;

  delete fOutputListSupportHistos;
}

// ############################################################################
// ############################################################################
AliAnalysisTaskEtaReconstruction::AliAnalysisTaskEtaReconstruction(): AliAnalysisTaskSE(), fEventFilter(0x0)
                                                                              , fdebug(false), run1analysis()
                                                                              , fResoFile(0x0), fResoFilename(""), fResoFilenameFromAlien(""), fArrResoPt(0x0), fArrResoEta(0x0), fArrResoPhi_Pos(0x0), fArrResoPhi_Neg(0x0)
                                                                              , fOutputList(0x0), fSingleElectronList(0x0), fGeneratedPrimaryList(0x0), fGeneratedSecondaryList(0x0), fGeneratedSmearedPrimaryList(0x0), fGeneratedSmearedSecondaryList(0x0), fRecPrimaryList(0x0), fRecSecondaryList(0x0)
                                                                              , fGeneratedPrimaryPairsList(0x0), fGeneratedSecondaryPairsList(0x0), fGeneratedSmearedPrimaryPairsList(0x0), fGeneratedSmearedSecondaryPairsList(0x0), fPairList(0x0), fFourPairList(0x0), fGeneratedFourPairsList(0x0), fGeneratedSmearedFourPairsList(0x0), fResolutionList(0x0)
                                                                              , fPGen_DeltaP(0x0), fPGen_PrecOverPGen(0x0), fPtGen_DeltaPt(0x0), fPtGen_DeltaPtOverPtGen(0x0), fPtGen_PtRecOverPtGen(0x0), fPtGen_DeltaPt_wGenSmeared(0x0), fPtGen_DeltaPtOverPtGen_wGenSmeared(0x0), fPtGen_PtRecOverPtGen_wGenSmeared(0x0)
                                                                              , fPGen_DeltaEta(0x0), fPtGen_DeltaEta(0x0), fPGen_DeltaTheta(0x0), fPGen_DeltaPhi_Ele(0x0), fPGen_DeltaPhi_Pos(0x0), fPtGen_DeltaPhi_Ele(0x0)
                                                                              , fPtGen_DeltaPhi_Pos(0x0), fThetaGen_DeltaTheta(0x0), fPhiGen_DeltaPhi(0x0)
                                                                              , fPtBins(), fEtaBins(), fPhiBins(), fThetaBins()
                                                                              , fResolutionDeltaPtBins(), fResolutionRelPtBins(), fResolutionEtaBins(), fResolutionPhiBins(), fResolutionThetaBins()
                                                                              , fMassBins(), fPairPtBins()
                                                                              , fPtMin(0.), fPtMax(0.), fEtaMin(-99.), fEtaMax(99.)
                                                                              , fPtMinGen(0.), fPtMaxGen(0.), fEtaMinGen(-99.), fEtaMaxGen(99.)
                                                                              , fLowerMassCutPrimaries(), fUpperMassCutPrimaries(), fMassCutSecondaries(), fUpperPrimSecPreFilterMass(), fLowerPrimSecPreFilterMass(), fUpperSecSecPreFilterMass(), fLowerSecSecPreFilterMass(), fV0OnFlyStatus()
                                                                              , fSinglePrimaryLegMCSignal(), fSingleSecondaryLegMCSignal(), fPrimaryPairMCSignal(), fSecondaryPairMCSignal(), fFourPairMCSignal(), fPrimaryDielectronPairNotFromSameMother(), fSecondaryDielectronPairNotFromSameMother()
                                                                              , fGeneratorName(""), fGeneratorMCSignalName(""), fGeneratorULSSignalName(""), fGeneratorHashs(), fGeneratorMCSignalHashs(), fGeneratorULSSignalHashs(), fPIDResponse(0x0), fEvent(0x0), fMC(0x0), fTrack(0x0), isAOD(false), fSelectPhysics(false), fTriggerMask(0)
                                                                              , fTrackCuts_primary_PreFilter(), fTrackCuts_primary_standard(), fTrackCuts_secondary_PreFilter(), fTrackCuts_secondary_standard(), fPairCuts_primary(), fPairCuts_secondary_PreFilter(), fPairCuts_secondary_standard(), fUsedVars(0x0)
                                                                              , fSupportMCSignal(0), fSupportCutsetting(0)
                                                                              , fHistEvents(0x0), fHistEventStat(0x0), fHistCentrality(0x0), fHistVertex(0x0), fHistVertexContibutors(0x0), fHistNTracks(0x0)
                                                                              , fMinCentrality(0.), fMaxCentrality(100), fCentralityFile(0x0), fCentralityFilename(""), fHistCentralityCorrection(0x0)
                                                                              , fOutputListSupportHistos(0x0), fTrackCutListVecPrim(), fTrackCutListVecSec(), fPairCutListVecSec(), fFourPairCutListVec()
                                                                              , fHistGenPrimaryPosPart(), fHistGenPrimaryNegPart(), fHistGenSecondaryPosPart(), fHistGenSecondaryNegPart(), fHistGenSmearedPrimaryPosPart(), fHistGenSmearedPrimaryNegPart(), fHistGenSmearedSecondaryPosPart(), fHistGenSmearedSecondaryNegPart(), fHistRecPrimaryPosPart(), fHistRecPrimaryNegPart(), fHistRecSecondaryPosPart(), fHistRecSecondaryNegPart()
                                                                              , fHistGenPrimaryPair(), fHistGenSecondaryPair(), fHistGenSmearedPrimaryPair(), fHistGenSmearedSecondaryPair(), fHistRecPrimaryPair(), fHistRecSecondaryPair(), fHistGenFourPair(), fHistGenSmearedFourPair(), fHistRecFourPair(), fVecHistPrefilters()
                                                                              , fWriteLegsFromPair(false), fPtMinLegsFromPair(-99.), fPtMaxLegsFromPair(-99.), fEtaMinLegsFromPair(-99.), fEtaMaxLegsFromPair(-99.), fPhiMinLegsFromPair(-99.), fPhiMaxLegsFromPair(-99.), fOpAngleMinLegsFromPair(-99.), fOpAngleMaxLegsFromPair(-99.), fPtNBinsLegsFromPair(-99), fEtaNBinsLegsFromPair(-99), fPhiNBinsLegsFromPair(-99), fOpAngleNBinsLegsFromPair(-99), fTHnSparseGenSmearedLegsFromPrimaryPair(), fTHnSparseGenSmearedLegsFromSecondaryPair(), fTHnSparseRecLegsFromPrimaryPair(), fTHnSparseRecLegsFromSecondaryPair()
                                                                              , fDoPairing(false), fDoFourPairing(false), fUsePreFilter(false), fUseSecPreFilter(false), fDoMassCut(), fPhotonMass()
                                                                              , fGenNegPart_primary(), fGenPosPart_primary(), fGenNegPart_secondary(), fGenPosPart_secondary(), fGenSmearedNegPart_primary(), fGenSmearedPosPart_primary(), fGenSmearedNegPart_secondary(), fGenSmearedPosPart_secondary(), fRecNegPart_primary(), fRecPosPart_primary(), fRecNegPart_secondary(), fRecPosPart_secondary(), fGenPairVec_primary(), fGenPairVec_secondary(), fGenSmearedPairVec_primary(), fGenSmearedPairVec_secondary(), fRecPairVec_primary(), fRecPairVec_secondary(), fRecV0Pair(), fPreFilter_BadTracksLabel_primary()
                                                                              , fDoCocktailWeighting(false), fCocktailFilename(""), fCocktailFilenameFromAlien(""), fCocktailFile(0x0)
                                                                              , fPtPion(0x0), fPtEta(0x0), fPtEtaPrime(0x0), fPtRho(0x0), fPtOmega(0x0), fPtPhi(0x0), fPtJPsi(0x0),
                                                                              fPostPIDCntrdCorrTPC(0x0), fPostPIDWdthCorrTPC(0x0), fPostPIDCntrdCorrITS(0x0), fPostPIDWdthCorrITS(0x0), fPostPIDCntrdCorrTOF(0x0), fPostPIDWdthCorrTOF(0x0)
{
// ROOT IO constructor , don ’t allocate memory here !
}


// ############################################################################
// ############################################################################
AliAnalysisTaskEtaReconstruction::AliAnalysisTaskEtaReconstruction(const char * name) : AliAnalysisTaskSE(name), fEventFilter(0x0)
                                                                              , fdebug(false), run1analysis()
                                                                              , fResoFile(0x0), fResoFilename(""), fResoFilenameFromAlien(""), fArrResoPt(0x0), fArrResoEta(0x0), fArrResoPhi_Pos(0x0), fArrResoPhi_Neg(0x0)
                                                                              , fOutputList(0x0), fSingleElectronList(0x0), fGeneratedPrimaryList(0x0), fGeneratedSecondaryList(0x0), fGeneratedSmearedPrimaryList(0x0), fGeneratedSmearedSecondaryList(0x0), fRecPrimaryList(0x0), fRecSecondaryList(0x0)
                                                                              , fGeneratedPrimaryPairsList(0x0), fGeneratedSecondaryPairsList(0x0), fGeneratedSmearedPrimaryPairsList(0x0), fGeneratedSmearedSecondaryPairsList(0x0), fPairList(0x0), fFourPairList(0x0), fGeneratedFourPairsList(0x0), fGeneratedSmearedFourPairsList(0x0), fResolutionList(0x0)
                                                                              , fPGen_DeltaP(0x0), fPGen_PrecOverPGen(0x0), fPtGen_DeltaPt(0x0), fPtGen_DeltaPtOverPtGen(0x0), fPtGen_PtRecOverPtGen(0x0), fPtGen_DeltaPt_wGenSmeared(0x0), fPtGen_DeltaPtOverPtGen_wGenSmeared(0x0), fPtGen_PtRecOverPtGen_wGenSmeared(0x0)
                                                                              , fPGen_DeltaEta(0x0), fPtGen_DeltaEta(0x0), fPGen_DeltaTheta(0x0), fPGen_DeltaPhi_Ele(0x0), fPGen_DeltaPhi_Pos(0x0), fPtGen_DeltaPhi_Ele(0x0)
                                                                              , fPtGen_DeltaPhi_Pos(0x0), fThetaGen_DeltaTheta(0x0), fPhiGen_DeltaPhi(0x0)
                                                                              , fPtBins(), fEtaBins(), fPhiBins(), fThetaBins()
                                                                              , fResolutionDeltaPtBins(), fResolutionRelPtBins(), fResolutionEtaBins(), fResolutionPhiBins(), fResolutionThetaBins()
                                                                              , fMassBins(), fPairPtBins()
                                                                              , fPtMin(0.), fPtMax(0.), fEtaMin(-99.), fEtaMax(99.)
                                                                              , fPtMinGen(0.), fPtMaxGen(0.), fEtaMinGen(-99.), fEtaMaxGen(99.)
                                                                              , fLowerMassCutPrimaries(), fUpperMassCutPrimaries(), fMassCutSecondaries(), fUpperPrimSecPreFilterMass(), fLowerPrimSecPreFilterMass(), fUpperSecSecPreFilterMass(), fLowerSecSecPreFilterMass(), fV0OnFlyStatus()
                                                                              , fSinglePrimaryLegMCSignal(), fSingleSecondaryLegMCSignal(), fPrimaryPairMCSignal(), fSecondaryPairMCSignal(), fFourPairMCSignal(), fPrimaryDielectronPairNotFromSameMother(), fSecondaryDielectronPairNotFromSameMother()
                                                                              , fGeneratorName(""), fGeneratorMCSignalName(""), fGeneratorULSSignalName(""), fGeneratorHashs(), fGeneratorMCSignalHashs(), fGeneratorULSSignalHashs(), fPIDResponse(0x0), fEvent(0x0), fMC(0x0), fTrack(0x0), isAOD(false), fSelectPhysics(false), fTriggerMask(0)
                                                                              , fTrackCuts_primary_PreFilter(), fTrackCuts_primary_standard(), fTrackCuts_secondary_PreFilter(), fTrackCuts_secondary_standard(), fPairCuts_primary(), fPairCuts_secondary_PreFilter(), fPairCuts_secondary_standard(), fUsedVars(0x0)
                                                                              , fSupportMCSignal(0), fSupportCutsetting(0)
                                                                              , fHistEvents(0x0), fHistEventStat(0x0), fHistCentrality(0x0), fHistVertex(0x0), fHistVertexContibutors(0x0), fHistNTracks(0x0)
                                                                              , fMinCentrality(0.), fMaxCentrality(100), fCentralityFile(0x0), fCentralityFilename(""), fHistCentralityCorrection(0x0)
                                                                              , fOutputListSupportHistos(0x0), fTrackCutListVecPrim(), fTrackCutListVecSec(), fPairCutListVecSec(), fFourPairCutListVec()
                                                                              , fHistGenPrimaryPosPart(), fHistGenPrimaryNegPart(), fHistGenSecondaryPosPart(), fHistGenSecondaryNegPart(), fHistGenSmearedPrimaryPosPart(), fHistGenSmearedPrimaryNegPart(), fHistGenSmearedSecondaryPosPart(), fHistGenSmearedSecondaryNegPart(), fHistRecPrimaryPosPart(), fHistRecPrimaryNegPart(), fHistRecSecondaryPosPart(), fHistRecSecondaryNegPart()
                                                                              , fHistGenPrimaryPair(), fHistGenSecondaryPair(), fHistGenSmearedPrimaryPair(), fHistGenSmearedSecondaryPair(), fHistRecPrimaryPair(), fHistRecSecondaryPair(), fHistGenFourPair(), fHistGenSmearedFourPair(), fHistRecFourPair(), fVecHistPrefilters()
                                                                              , fWriteLegsFromPair(false), fPtMinLegsFromPair(-99.), fPtMaxLegsFromPair(-99.), fEtaMinLegsFromPair(-99.), fEtaMaxLegsFromPair(-99.), fPhiMinLegsFromPair(-99.), fPhiMaxLegsFromPair(-99.), fOpAngleMinLegsFromPair(-99.), fOpAngleMaxLegsFromPair(-99.), fPtNBinsLegsFromPair(-99), fEtaNBinsLegsFromPair(-99), fPhiNBinsLegsFromPair(-99), fOpAngleNBinsLegsFromPair(-99), fTHnSparseGenSmearedLegsFromPrimaryPair(), fTHnSparseGenSmearedLegsFromSecondaryPair(), fTHnSparseRecLegsFromPrimaryPair(), fTHnSparseRecLegsFromSecondaryPair()
                                                                              , fDoPairing(false), fDoFourPairing(false), fUsePreFilter(false), fUseSecPreFilter(false), fDoMassCut(), fPhotonMass()
                                                                              , fGenNegPart_primary(), fGenPosPart_primary(), fGenNegPart_secondary(), fGenPosPart_secondary(), fGenSmearedNegPart_primary(), fGenSmearedPosPart_primary(), fGenSmearedNegPart_secondary(), fGenSmearedPosPart_secondary(), fRecNegPart_primary(), fRecPosPart_primary(), fRecNegPart_secondary(), fRecPosPart_secondary(), fGenPairVec_primary(), fGenPairVec_secondary(), fGenSmearedPairVec_primary(), fGenSmearedPairVec_secondary(), fRecPairVec_primary(), fRecPairVec_secondary(), fRecV0Pair(), fPreFilter_BadTracksLabel_primary()
                                                                              , fDoCocktailWeighting(false), fCocktailFilename(""), fCocktailFilenameFromAlien(""), fCocktailFile(0x0)
                                                                              , fPtPion(0x0), fPtEta(0x0), fPtEtaPrime(0x0), fPtRho(0x0), fPtOmega(0x0), fPtPhi(0x0), fPtJPsi(0x0),
                                                                              fPostPIDCntrdCorrTPC(0x0), fPostPIDWdthCorrTPC(0x0), fPostPIDCntrdCorrITS(0x0), fPostPIDWdthCorrITS(0x0), fPostPIDCntrdCorrTOF(0x0), fPostPIDWdthCorrTOF(0x0)

{
  DefineInput (0, TChain::Class());
  DefineOutput (1, TList::Class());



}


// ############################################################################
// ############################################################################
void AliAnalysisTaskEtaReconstruction::Terminate(Option_t* option){
  // fHistEventStat->SetAxisRange(0., fHistEventStat->GetMaximum() * 1.1, "Y");
}


// ############################################################################
// ############################################################################
void AliAnalysisTaskEtaReconstruction::UserCreateOutputObjects(){
                                                                                if(fdebug) gRandom->SetSeed(122);  // Seed in order to proove consistance of code durig changes. Disables the randomness of smearing. Comment out for normal calculation
  fUsedVars = new TBits(AliDielectronVarManager::kNMaxValues);
  fUsedVars->SetBitNumber(AliDielectronVarManager::kP, kTRUE);
  fUsedVars->SetBitNumber(AliDielectronVarManager::kPIn, kTRUE);
  fUsedVars->SetBitNumber(AliDielectronVarManager::kITSnSigmaEle, kTRUE);
  fUsedVars->SetBitNumber(AliDielectronVarManager::kITSnSigmaMuo, kTRUE);
  fUsedVars->SetBitNumber(AliDielectronVarManager::kITSnSigmaPio, kTRUE);
  fUsedVars->SetBitNumber(AliDielectronVarManager::kITSnSigmaKao, kTRUE);
  fUsedVars->SetBitNumber(AliDielectronVarManager::kITSnSigmaPro, kTRUE);
  fUsedVars->SetBitNumber(AliDielectronVarManager::kTPCnSigmaEle, kTRUE);
  fUsedVars->SetBitNumber(AliDielectronVarManager::kTPCnSigmaMuo, kTRUE);
  fUsedVars->SetBitNumber(AliDielectronVarManager::kTPCnSigmaPio, kTRUE);
  fUsedVars->SetBitNumber(AliDielectronVarManager::kTPCnSigmaKao, kTRUE);
  fUsedVars->SetBitNumber(AliDielectronVarManager::kTPCnSigmaPro, kTRUE);
  fUsedVars->SetBitNumber(AliDielectronVarManager::kTOFnSigmaEle, kTRUE);
  fUsedVars->SetBitNumber(AliDielectronVarManager::kTOFnSigmaMuo, kTRUE);
  fUsedVars->SetBitNumber(AliDielectronVarManager::kTOFnSigmaPio, kTRUE);
  fUsedVars->SetBitNumber(AliDielectronVarManager::kTOFnSigmaKao, kTRUE);
  fUsedVars->SetBitNumber(AliDielectronVarManager::kTOFnSigmaPro, kTRUE);
  fUsedVars->SetBitNumber(AliDielectronVarManager::kNFclsTPCr, kTRUE);
  fUsedVars->SetBitNumber(AliDielectronVarManager::kNFclsTPCfCross, kTRUE);
  fUsedVars->SetBitNumber(AliDielectronVarManager::kNclsTPC, kTRUE);
  fUsedVars->SetBitNumber(AliDielectronVarManager::kNclsITS, kTRUE);
  fUsedVars->SetBitNumber(AliDielectronVarManager::kTPCchi2Cl, kTRUE);
  fUsedVars->SetBitNumber(AliDielectronVarManager::kITSchi2Cl, kTRUE);
  fUsedVars->SetBitNumber(AliDielectronVarManager::kNclsSTPC, kTRUE);
  fUsedVars->SetBitNumber(AliDielectronVarManager::kNclsSITS, kTRUE);
  fUsedVars->SetBitNumber(AliDielectronVarManager::kNclsSFracTPC, kTRUE);
  fUsedVars->SetBitNumber(AliDielectronVarManager::kNclsSFracITS, kTRUE);
  fUsedVars->SetBitNumber(AliDielectronVarManager::kTPCclsDiff, kTRUE);
  fUsedVars->SetBitNumber(AliDielectronVarManager::kTPCsignalN, kTRUE);
  fUsedVars->SetBitNumber(AliDielectronVarManager::kNclsTPC, kTRUE);
  fUsedVars->SetBitNumber(AliDielectronVarManager::kNFclsTPCr, kTRUE);
  fUsedVars->SetBitNumber(AliDielectronVarManager::kImpactParXY, kTRUE);
  fUsedVars->SetBitNumber(AliDielectronVarManager::kImpactParZ, kTRUE);
  fUsedVars->SetBitNumber(AliDielectronVarManager::kCosPointingAngle, kTRUE);
  fUsedVars->SetBitNumber(AliDielectronVarManager::kChi2NDF, kTRUE);
  fUsedVars->SetBitNumber(AliDielectronVarManager::kLegDist, kTRUE);
  fUsedVars->SetBitNumber(AliDielectronVarManager::kR, kTRUE);
  fUsedVars->SetBitNumber(AliDielectronVarManager::kPsiPair, kTRUE);
  fUsedVars->SetBitNumber(AliDielectronVarManager::kM, kTRUE);
  fUsedVars->SetBitNumber(AliDielectronVarManager::kArmPt, kTRUE);
  fUsedVars->SetBitNumber(AliDielectronVarManager::kArmAlpha, kTRUE);
  fUsedVars->SetBitNumber(AliDielectronVarManager::kPt, kTRUE);
  AliDielectronVarManager::SetFillMap(fUsedVars); // currently filled manually in the constructor of this task.



  std::cout << "Starting UserCreateOutputObjects()" << std::endl;

  TObjArray arr = *(fGeneratorName.Tokenize(";"));
  std::cout << "Used Generators: " << std::endl;
  for (int i = 0; i < arr.GetEntries(); ++i){
    TString temp = arr.At(i)->GetName();
    std::cout << "--- " << temp << std::endl;
    fGeneratorHashs.push_back(temp.Hash());
  }

  arr = *(fGeneratorMCSignalName.Tokenize(";"));
  std::cout << "Used Generators for MCSignals: " << std::endl;
  for (int i = 0; i < arr.GetEntries(); ++i){
    TString temp = arr.At(i)->GetName();
    std::cout << "--- " << temp << std::endl;
    fGeneratorMCSignalHashs.push_back(temp.Hash());
  }

  arr = *(fGeneratorULSSignalName.Tokenize(";"));
  std::cout << "Used Generators for ULSSignals: " << std::endl;
  for (int i = 0; i < arr.GetEntries(); ++i){
    TString temp = arr.At(i)->GetName();
    std::cout << "--- " << temp << std::endl;
    fGeneratorULSSignalHashs.push_back(temp.Hash());
  }

  if (fResoFilename != ""){
    fResoFile = TFile::Open(fResoFilename.c_str());
    if (fResoFile == 0x0){
      std::cout << "Location in AliEN: " << fResoFilenameFromAlien << std::endl;
      gSystem->Exec(Form("alien_cp alien://%s .", fResoFilenameFromAlien.c_str()));
      std::cout << "Copy resolution from Alien" << std::endl;
      fResoFile = TFile::Open(fResoFilename.c_str());
    }

    if (!fResoFile->IsOpen()) {
      AliError(Form("Could not open file %s", fResoFilename.c_str()));
    }
    fArrResoPt = (TObjArray *)fResoFile->Get("RelPtResArrCocktail");
    fArrResoEta = (TObjArray *)fResoFile->Get("EtaResArrVsPt");
    fArrResoPhi_Pos = (TObjArray *)fResoFile->Get("PhiPosResArrVsPt");
    fArrResoPhi_Neg = (TObjArray *)fResoFile->Get("PhiEleResArrVsPt");
    std::cout << fArrResoPt << " " << fArrResoEta << " " << fArrResoPhi_Pos << " " << fArrResoPhi_Neg << std::endl;
    if (fArrResoPt == 0x0 ||  fArrResoEta == 0x0 || fArrResoPhi_Pos == 0x0 || fArrResoPhi_Neg == 0x0){
      AliError(Form("Could not extract resolution histograms from file %s", fResoFilename.c_str()));
    }
  }

  if (fDoCocktailWeighting && fCocktailFilename != ""){
    std::cout << "Do Cocktail weighting" << std::endl;
    fCocktailFile = TFile::Open(fCocktailFilename.c_str());
    if (fCocktailFile == 0x0){
      std::cout << "Location in AliEN: " << fCocktailFilenameFromAlien << std::endl;
      gSystem->Exec(Form("alien_cp alien://%s .", fCocktailFilenameFromAlien.c_str()));
      std::cout << "Copy cocktail weighting from Alien" << std::endl;
      fCocktailFile = TFile::Open(fCocktailFilename.c_str());
    }

    if (fCocktailFile){
      fPtPion     = dynamic_cast<TH1F*>(fCocktailFile->Get("Pion"));
      fPtEta      = dynamic_cast<TH1F*>(fCocktailFile->Get("Eta"));
      fPtEtaPrime = dynamic_cast<TH1F*>(fCocktailFile->Get("EtaPrime"));
      fPtRho      = dynamic_cast<TH1F*>(fCocktailFile->Get("Rho"));
      fPtOmega    = dynamic_cast<TH1F*>(fCocktailFile->Get("Omega"));
      fPtPhi      = dynamic_cast<TH1F*>(fCocktailFile->Get("Phi"));
      fPtJPsi     = dynamic_cast<TH1F*>(fCocktailFile->Get("JPsi"));

      if (!fPtPion)     { std::cout << "Pion reweighting not loaded"     << std::endl; }
      if (!fPtEta)      { std::cout << "Eta reweighting not loaded"      << std::endl; }
      if (!fPtEtaPrime) { std::cout << "EtaPrime reweighting not loaded" << std::endl; }
      if (!fPtRho)      { std::cout << "Rho reweighting not loaded"      << std::endl; }
      if (!fPtOmega)    { std::cout << "Omega reweighting not loaded"    << std::endl; }
      if (!fPtPhi)      { std::cout << "Phi reweighting not loaded"      << std::endl; }
      if (!fPtJPsi)     { std::cout << "JPsi reweighting not loaded"     << std::endl; }
    }
    else std::cout << "No cocktail weighting file found" << std::endl;
  }

  if (fCentralityFilename != ""){
    fCentralityFile = TFile::Open(fCentralityFilename.c_str());
    if (!fCentralityFile->IsOpen()) {
      AliError(Form("Could not open file %s", fCentralityFilename.c_str()));
    }
    TList* list_temp = (TList*)fCentralityFile->Get("efficiency");
    fHistCentralityCorrection = (TH1F*) list_temp->FindObject("centrality");
    if (fHistCentralityCorrection == 0x0){
      AliError(Form("Could not extract centrality histogram from file %s", fCentralityFilename.c_str()));
    }
  }

  // Check binning for single electron histograms. All 3 dimension must have >= 1 bin
  const int fNptBins = fPtBins.size()-1;
  const int fNetaBins = fEtaBins.size()-1;
  const int fNphiBins = fPhiBins.size()-1;
  const int fNthetaBins = fThetaBins.size()-1;
  const int fNResolutionDeltaptBins = fResolutionDeltaPtBins.size()-1;
  const int fNResolutionRelptBins = fResolutionRelPtBins.size()-1;
  const int fNResolutionetaBins = fResolutionEtaBins.size()-1;
  const int fNResolutionphiBins = fResolutionPhiBins.size()-1;
  const int fNResolutionthetaBins = fResolutionThetaBins.size()-1;
  const int fNmassBins = fMassBins.size()-1;
  const int fNpairptBins = fPairPtBins.size()-1;

  const int nDim = 7;
  Int_t nBins[nDim] = {fPtNBinsLegsFromPair, fEtaNBinsLegsFromPair, fPhiNBinsLegsFromPair, fPtNBinsLegsFromPair, fEtaNBinsLegsFromPair, fPhiNBinsLegsFromPair, fOpAngleNBinsLegsFromPair};
  Double_t min[nDim] = {fPtMinLegsFromPair, fEtaMinLegsFromPair, fPhiMinLegsFromPair, fPtMinLegsFromPair, fEtaMinLegsFromPair, fPhiMinLegsFromPair, fOpAngleMinLegsFromPair};
  Double_t max[nDim] = {fPtMaxLegsFromPair, fEtaMaxLegsFromPair, fPhiMaxLegsFromPair, fPtMaxLegsFromPair, fEtaMaxLegsFromPair, fPhiMaxLegsFromPair, fOpAngleMaxLegsFromPair};

  if (fNptBins < 2|| fNetaBins < 2 || fNphiBins < 2 || fNthetaBins < 2){
    std::cout << "No Pt, Eta and/or Phi binning given: #ptBins=" << fNptBins << " #etaBins=" << fNetaBins << " #phiBins=" << fNphiBins << " #thetaBins=" << fNthetaBins << std::endl;
    return;
  }

  fOutputList = new TList();
  fOutputList->SetOwner();

  // Initialize all histograms
    fHistEvents             = new TH1F("events", "events", 1, 0., 1.);
    fHistEventStat          = new TH1F("eventStats", "eventStats", kLastBin, -0.5, kLastBin-0.5);
    fHistCentrality         = new TH1F("centrality", "centrality", 100, 0., 100.);
    fHistVertex             = new TH1F("zVertex", "zVertex", 300, -15.0, 15.0);
    fHistVertexContibutors  = new TH1F("vtxContributor", "vtxContributor",5000,-0.5,4999.5);
    fHistNTracks            = new TH1F("nTracks", "nTracks", 4000, 0., 40000.);
    fOutputList->Add(fHistEvents);
    fOutputList->Add(fHistEventStat);
    fOutputList->Add(fHistCentrality);
    fOutputList->Add(fHistVertex);
    // fOutputList->Add(fHistVertexContibutors);
    fOutputList->Add(fHistNTracks);


    // ######################################################
    // ##########  Single Electrons #########################
    // ######################################################
    fSingleElectronList = new TList();
      fSingleElectronList->SetName("SingleElectrons");
      fSingleElectronList->SetOwner();
      // Create List with primary generated particles
      fGeneratedPrimaryList = new TList();
      fGeneratedPrimaryList->SetName("Generated_Primary");
      fGeneratedPrimaryList->SetOwner();
      for (unsigned int i = 0; i < fSinglePrimaryLegMCSignal.size(); ++i){
        TH3D* th3_tmp_pos = new TH3D(Form("Ngen_Pos_%s", fSinglePrimaryLegMCSignal.at(i).GetName()),";p_{T};#eta;#varphi",fNptBins,fPtBins.data(),fNetaBins,fEtaBins.data(),fNphiBins,fPhiBins.data());
        th3_tmp_pos->Sumw2();
        fHistGenPrimaryPosPart.push_back(th3_tmp_pos);
        fGeneratedPrimaryList->Add(th3_tmp_pos);
        TH3D* th3_tmp_neg = new TH3D(Form("Ngen_Neg_%s", fSinglePrimaryLegMCSignal.at(i).GetName()),";p_{T};#eta;#varphi",fNptBins,fPtBins.data(),fNetaBins,fEtaBins.data(),fNphiBins,fPhiBins.data());
        th3_tmp_neg->Sumw2();
        fHistGenPrimaryNegPart.push_back(th3_tmp_neg);
        fGeneratedPrimaryList->Add(th3_tmp_neg);
      }
      // Create List with secondary generated particles
      fGeneratedSecondaryList = new TList();
      fGeneratedSecondaryList->SetName("Generated_Secondary");
      fGeneratedSecondaryList->SetOwner();
      for (unsigned int i = 0; i < fSingleSecondaryLegMCSignal.size(); ++i){
        TH3D* th3_tmp_pos = new TH3D(Form("Ngen_Pos_%s", fSingleSecondaryLegMCSignal.at(i).GetName()),";p_{T};#eta;#varphi",fNptBins,fPtBins.data(),fNetaBins,fEtaBins.data(),fNphiBins,fPhiBins.data());
        th3_tmp_pos->Sumw2();
        fHistGenSecondaryPosPart.push_back(th3_tmp_pos);
        fGeneratedSecondaryList->Add(th3_tmp_pos);
        TH3D* th3_tmp_neg = new TH3D(Form("Ngen_Neg_%s", fSingleSecondaryLegMCSignal.at(i).GetName()),";p_{T};#eta;#varphi",fNptBins,fPtBins.data(),fNetaBins,fEtaBins.data(),fNphiBins,fPhiBins.data());
        th3_tmp_neg->Sumw2();
        fHistGenSecondaryNegPart.push_back(th3_tmp_neg);
        fGeneratedSecondaryList->Add(th3_tmp_neg);
      }

      // Create List with generated+smeared particles
      fGeneratedSmearedPrimaryList = new TList();
      fGeneratedSmearedPrimaryList->SetName("GeneratedSmeared_Primary");
      fGeneratedSmearedPrimaryList->SetOwner();
      for (unsigned int i = 0; i < fSinglePrimaryLegMCSignal.size(); ++i){
        TH3D* th3_tmp_pos = new TH3D(Form("Ngen_Pos_%s", fSinglePrimaryLegMCSignal.at(i).GetName()),";p_{T};#eta;#varphi",fNptBins,fPtBins.data(),fNetaBins,fEtaBins.data(),fNphiBins,fPhiBins.data());
        th3_tmp_pos->Sumw2();
        fHistGenSmearedPrimaryPosPart.push_back(th3_tmp_pos);
        fGeneratedSmearedPrimaryList->Add(th3_tmp_pos);
        TH3D* th3_tmp_neg = new TH3D(Form("Ngen_Neg_%s", fSinglePrimaryLegMCSignal.at(i).GetName()),";p_{T};#eta;#varphi",fNptBins,fPtBins.data(),fNetaBins,fEtaBins.data(),fNphiBins,fPhiBins.data());
        th3_tmp_neg->Sumw2();
        fHistGenSmearedPrimaryNegPart.push_back(th3_tmp_neg);
        fGeneratedSmearedPrimaryList->Add(th3_tmp_neg);
      }
      // Create List with generated+smeared particles
      fGeneratedSmearedSecondaryList = new TList();
      fGeneratedSmearedSecondaryList->SetName("GeneratedSmeared_Secondary");
      fGeneratedSmearedSecondaryList->SetOwner();
      for (unsigned int i = 0; i < fSingleSecondaryLegMCSignal.size(); ++i){
        TH3D* th3_tmp_pos = new TH3D(Form("Ngen_Pos_%s", fSingleSecondaryLegMCSignal.at(i).GetName()),";p_{T};#eta;#varphi",fNptBins,fPtBins.data(),fNetaBins,fEtaBins.data(),fNphiBins,fPhiBins.data());
        th3_tmp_pos->Sumw2();
        fHistGenSmearedSecondaryPosPart.push_back(th3_tmp_pos);
        fGeneratedSmearedSecondaryList->Add(th3_tmp_pos);
        TH3D* th3_tmp_neg = new TH3D(Form("Ngen_Neg_%s", fSingleSecondaryLegMCSignal.at(i).GetName()),";p_{T};#eta;#varphi",fNptBins,fPtBins.data(),fNetaBins,fEtaBins.data(),fNphiBins,fPhiBins.data());
        th3_tmp_neg->Sumw2();
        fHistGenSmearedSecondaryNegPart.push_back(th3_tmp_neg);
        fGeneratedSmearedSecondaryList->Add(th3_tmp_neg);
      }

      fSingleElectronList->Add(fGeneratedPrimaryList);
      fSingleElectronList->Add(fGeneratedSecondaryList);
      fSingleElectronList->Add(fGeneratedSmearedPrimaryList);
      fSingleElectronList->Add(fGeneratedSmearedSecondaryList);

      // Generated reconstructed lists for every cutsetting one list and every MCsignal 2 histograms with pos and neg charge
      for (unsigned int list_i = 0; list_i < fTrackCuts_primary_standard.size(); ++list_i){
        fRecPrimaryList = new TList();
        fRecPrimaryList->SetName(Form("%s_Primary",fTrackCuts_primary_standard.at(list_i)->GetName()));
        fRecPrimaryList->SetOwner();

        for (unsigned int i = 0; i < fSinglePrimaryLegMCSignal.size(); ++i){
          TH3D* th3_tmp_pos = new TH3D(Form("Nrec_Pos_%s", fSinglePrimaryLegMCSignal.at(i).GetName()),";p_{T};#eta;#varphi",fNptBins,fPtBins.data(),fNetaBins,fEtaBins.data(),fNphiBins,fPhiBins.data());
          th3_tmp_pos->Sumw2();
          th3_tmp_pos->SetDirectory(0x0);
          fHistRecPrimaryPosPart.push_back(th3_tmp_pos);
          fRecPrimaryList->Add(th3_tmp_pos);
          TH3D* th3_tmp_neg = new TH3D(Form("Nrec_Neg_%s", fSinglePrimaryLegMCSignal.at(i).GetName()),";p_{T};#eta;#varphi",fNptBins,fPtBins.data(),fNetaBins,fEtaBins.data(),fNphiBins,fPhiBins.data());
          th3_tmp_neg->Sumw2();
          th3_tmp_neg->SetDirectory(0x0);
          fHistRecPrimaryNegPart.push_back(th3_tmp_neg);
          fRecPrimaryList->Add(th3_tmp_neg);

        }
        fSingleElectronList->Add(fRecPrimaryList);
      }
      for (unsigned int list_i = 0; list_i < fPairCuts_secondary_standard.size(); ++list_i){
        fRecSecondaryList = new TList();
        fRecSecondaryList->SetName(Form("%s_Secondary",fPairCuts_secondary_standard.at(list_i)->GetName()));
        fRecSecondaryList->SetOwner();

        for (unsigned int i = 0; i < fSingleSecondaryLegMCSignal.size(); ++i){
          TH3D* th3_tmp_pos = new TH3D(Form("Nrec_Pos_%s", fSingleSecondaryLegMCSignal.at(i).GetName()),";p_{T};#eta;#varphi",fNptBins,fPtBins.data(),fNetaBins,fEtaBins.data(),fNphiBins,fPhiBins.data());
          th3_tmp_pos->Sumw2();
          th3_tmp_pos->SetDirectory(0x0);
          fHistRecSecondaryPosPart.push_back(th3_tmp_pos);
          fRecSecondaryList->Add(th3_tmp_pos);
          TH3D* th3_tmp_neg = new TH3D(Form("Nrec_Neg_%s", fSingleSecondaryLegMCSignal.at(i).GetName()),";p_{T};#eta;#varphi",fNptBins,fPtBins.data(),fNetaBins,fEtaBins.data(),fNphiBins,fPhiBins.data());
          th3_tmp_neg->Sumw2();
          th3_tmp_neg->SetDirectory(0x0);
          fHistRecSecondaryNegPart.push_back(th3_tmp_neg);
          fRecSecondaryList->Add(th3_tmp_neg);

        }
        fSingleElectronList->Add(fRecSecondaryList);
      }



      // ######################################################
      // #####################  PAIRS #########################
      // ######################################################

      if (fDoPairing == true){
        fPairList = new TList();
        fPairList->SetName("Pairs");
        fPairList->SetOwner();

        fGeneratedPrimaryPairsList = new TList();
        fGeneratedPrimaryPairsList->SetName("Generated_Primary");
        fGeneratedPrimaryPairsList->SetOwner();
        for (unsigned int i = 0; i < fPrimaryPairMCSignal.size(); ++i){
          TH2D* th2_tmp = new TH2D(Form("Ngen_%s", fPrimaryPairMCSignal.at(i).GetName()),";m_{ee};p_{T,ee}",fNmassBins,fMassBins.data(),fNpairptBins,fPairPtBins.data());
          th2_tmp->Sumw2();
          fHistGenPrimaryPair.push_back(th2_tmp);
          fGeneratedPrimaryPairsList->Add(th2_tmp);
        }

        fGeneratedSecondaryPairsList = new TList();
        fGeneratedSecondaryPairsList->SetName("Generated_Secondary");
        fGeneratedSecondaryPairsList->SetOwner();
        for (unsigned int i = 0; i < fSecondaryPairMCSignal.size(); ++i){
          TH2D* th2_tmp = new TH2D(Form("Ngen_%s", fSecondaryPairMCSignal.at(i).GetName()),";m_{ee};p_{T,ee}",fNmassBins,fMassBins.data(),fNpairptBins,fPairPtBins.data());
          th2_tmp->Sumw2();
          fHistGenSecondaryPair.push_back(th2_tmp);
          fGeneratedSecondaryPairsList->Add(th2_tmp);
        }


        fGeneratedSmearedPrimaryPairsList = new TList();
        fGeneratedSmearedPrimaryPairsList->SetName("GeneratedSmeared_Primary");
        fGeneratedSmearedPrimaryPairsList->SetOwner();
        for (unsigned int i = 0; i < fPrimaryPairMCSignal.size(); ++i){
          TH2D* th2_tmp = new TH2D(Form("Ngen_%s", fPrimaryPairMCSignal.at(i).GetName()),";m_{ee};p_{T,ee}",fNmassBins,fMassBins.data(),fNpairptBins,fPairPtBins.data());
          th2_tmp->Sumw2();
          fHistGenSmearedPrimaryPair.push_back(th2_tmp);
          fGeneratedSmearedPrimaryPairsList->Add(th2_tmp);

          if (fWriteLegsFromPair){

            THnSparseF* fTHnSparseGenSmearedLegsFromPrimaryPair_tmp= new THnSparseF(Form("fTHnSparseGenSmearedLegsFromPrimaryPair_%s", fPrimaryPairMCSignal.at(i).GetName()),Form("fTHnSparseGenSmearedLegsFromPrimaryPair_%s;p_{t,Pos};#eta_{Pos};#phi_{Pos};p_{t,Neg};#eta_{Neg};#phi_{Neg};opAngle", fPrimaryPairMCSignal.at(i).GetName()), nDim, nBins, min, max);
            fTHnSparseGenSmearedLegsFromPrimaryPair_tmp->GetAxis(0)->SetName("ptPos");
            fTHnSparseGenSmearedLegsFromPrimaryPair_tmp->GetAxis(1)->SetName("etaPos");
            fTHnSparseGenSmearedLegsFromPrimaryPair_tmp->GetAxis(2)->SetName("phiPos");
            fTHnSparseGenSmearedLegsFromPrimaryPair_tmp->GetAxis(3)->SetName("ptNeg");
            fTHnSparseGenSmearedLegsFromPrimaryPair_tmp->GetAxis(4)->SetName("etaNeg");
            fTHnSparseGenSmearedLegsFromPrimaryPair_tmp->GetAxis(5)->SetName("phiNeg");
            fTHnSparseGenSmearedLegsFromPrimaryPair_tmp->GetAxis(6)->SetName("opAngle");
            fTHnSparseGenSmearedLegsFromPrimaryPair.push_back(fTHnSparseGenSmearedLegsFromPrimaryPair_tmp);
            fGeneratedSmearedPrimaryPairsList->Add(fTHnSparseGenSmearedLegsFromPrimaryPair_tmp);
          }
        }

        fGeneratedSmearedSecondaryPairsList = new TList();
        fGeneratedSmearedSecondaryPairsList->SetName("GeneratedSmeared_Secondary");
        fGeneratedSmearedSecondaryPairsList->SetOwner();
        for (unsigned int i = 0; i < fSecondaryPairMCSignal.size(); ++i){
          TH2D* th2_tmp = new TH2D(Form("Ngen_%s", fSecondaryPairMCSignal.at(i).GetName()),";m_{ee};p_{T,ee}",fNmassBins,fMassBins.data(),fNpairptBins,fPairPtBins.data());
          th2_tmp->Sumw2();
          fHistGenSmearedSecondaryPair.push_back(th2_tmp);
          fGeneratedSmearedSecondaryPairsList->Add(th2_tmp);

          if (fWriteLegsFromPair){

            THnSparseF* fTHnSparseGenSmearedLegsFromSecondaryPair_tmp= new THnSparseF(Form("fTHnSparseGenSmearedLegsFromSecondaryPair_%s", fSecondaryPairMCSignal.at(i).GetName()),Form("fTHnSparseGenSmearedLegsFromSecondaryPair_%s;p_{t,Pos};#eta_{Pos};#phi_{Pos};p_{t,Neg};#eta_{Neg};#phi_{Neg};opAngle", fSecondaryPairMCSignal.at(i).GetName()), nDim, nBins, min, max);
            fTHnSparseGenSmearedLegsFromSecondaryPair_tmp->GetAxis(0)->SetName("ptPos");
            fTHnSparseGenSmearedLegsFromSecondaryPair_tmp->GetAxis(1)->SetName("etaPos");
            fTHnSparseGenSmearedLegsFromSecondaryPair_tmp->GetAxis(2)->SetName("phiPos");
            fTHnSparseGenSmearedLegsFromSecondaryPair_tmp->GetAxis(3)->SetName("ptNeg");
            fTHnSparseGenSmearedLegsFromSecondaryPair_tmp->GetAxis(4)->SetName("etaNeg");
            fTHnSparseGenSmearedLegsFromSecondaryPair_tmp->GetAxis(5)->SetName("phiNeg");
            fTHnSparseGenSmearedLegsFromSecondaryPair_tmp->GetAxis(6)->SetName("opAngle");
            fTHnSparseGenSmearedLegsFromSecondaryPair.push_back(fTHnSparseGenSmearedLegsFromSecondaryPair_tmp);
            fGeneratedSmearedSecondaryPairsList->Add(fTHnSparseGenSmearedLegsFromSecondaryPair_tmp);
          }
        }

        fPairList->Add(fGeneratedPrimaryPairsList);
        fPairList->Add(fGeneratedSecondaryPairsList);
        fPairList->Add(fGeneratedSmearedPrimaryPairsList);
        fPairList->Add(fGeneratedSmearedSecondaryPairsList);

        // Generated reconstructed lists for every cutsetting one list and every MCsignal 1 histogram
        for (unsigned int list_i = 0; list_i < fTrackCuts_primary_standard.size(); ++list_i){
          TList* list = new TList();
          list->SetName(Form("%s_Primary",fTrackCuts_primary_standard.at(list_i)->GetName()));
          list->SetOwner();

          for (unsigned int i = 0; i < fPrimaryPairMCSignal.size(); ++i){
            TH2D* th2_tmp = new TH2D(Form("Nrec_%s", fPrimaryPairMCSignal.at(i).GetName()),";m_{ee};p_{T,ee}",fNmassBins,fMassBins.data(),fNpairptBins,fPairPtBins.data());
            th2_tmp->Sumw2();
            fHistRecPrimaryPair.push_back(th2_tmp);
            list->Add(th2_tmp);

            if (fWriteLegsFromPair){
              THnSparseF* fTHnSparseRecLegsFromPrimaryPair_tmp= new THnSparseF(Form("fTHnSparseRecLegsFromPrimaryPair_%s", fPrimaryPairMCSignal.at(i).GetName()),Form("fTHnSparseRecLegsFromPair_%s;p_{t,Pos};#eta_{Pos};#phi_{Pos};p_{t,Neg};#eta_{Neg};#phi_{Neg};opAngle", fPrimaryPairMCSignal.at(i).GetName()), nDim, nBins, min, max);
              fTHnSparseRecLegsFromPrimaryPair_tmp->GetAxis(0)->SetName("ptPos");
              fTHnSparseRecLegsFromPrimaryPair_tmp->GetAxis(1)->SetName("etaPos");
              fTHnSparseRecLegsFromPrimaryPair_tmp->GetAxis(2)->SetName("phiPos");
              fTHnSparseRecLegsFromPrimaryPair_tmp->GetAxis(3)->SetName("ptNeg");
              fTHnSparseRecLegsFromPrimaryPair_tmp->GetAxis(4)->SetName("etaNeg");
              fTHnSparseRecLegsFromPrimaryPair_tmp->GetAxis(5)->SetName("phiNeg");
              fTHnSparseRecLegsFromPrimaryPair_tmp->GetAxis(6)->SetName("opAngle");
              fTHnSparseRecLegsFromPrimaryPair_tmp->SetName(Form("fTHnSparseRecLegsFromPrimaryPairTest_%s;ptPos;etaPos;phiPos;ptNeg;etaNeg;phiNeg;opAngle", fPrimaryPairMCSignal.at(i).GetName()));
              fTHnSparseRecLegsFromPrimaryPair_tmp->SetTitle(Form("fTHnSparseRecLegsFromPrimaryPairTest_%s;ptPos;etaPos;phiPos;ptNeg;etaNeg;phiNeg;opAngle", fPrimaryPairMCSignal.at(i).GetName()));
              fTHnSparseRecLegsFromPrimaryPair.push_back(fTHnSparseRecLegsFromPrimaryPair_tmp);
              list->Add(fTHnSparseRecLegsFromPrimaryPair_tmp);
            }
          }

          fPairList->Add(list);
        }

        // Generated reconstructed lists for every cutsetting one list and every MCsignal 1 histogram
        for (unsigned int list_i = 0; list_i < fPairCuts_secondary_standard.size(); ++list_i){
          TList* list1 = new TList();
          list1->SetName(Form("%s_Secondary",fPairCuts_secondary_standard.at(list_i)->GetName()));
          list1->SetOwner();

          for (unsigned int i = 0; i < fSecondaryPairMCSignal.size(); ++i){
            TH2D* th2_tmp = new TH2D(Form("Nrec_%s", fSecondaryPairMCSignal.at(i).GetName()),";m_{ee};p_{T,ee}",fNmassBins,fMassBins.data(),fNpairptBins,fPairPtBins.data());
            th2_tmp->Sumw2();
            fHistRecSecondaryPair.push_back(th2_tmp);
            list1->Add(th2_tmp);

            if (fWriteLegsFromPair){
              THnSparseF* fTHnSparseRecLegsFromSecondaryPair_tmp= new THnSparseF(Form("fTHnSparseRecLegsFromSecondaryPair_%s", fSecondaryPairMCSignal.at(i).GetName()),Form("fTHnSparseRecLegsFromPair_%s;p_{t,Pos};#eta_{Pos};#phi_{Pos};p_{t,Neg};#eta_{Neg};#phi_{Neg};opAngle", fSecondaryPairMCSignal.at(i).GetName()), nDim, nBins, min, max);
              fTHnSparseRecLegsFromSecondaryPair_tmp->GetAxis(0)->SetName("ptPos");
              fTHnSparseRecLegsFromSecondaryPair_tmp->GetAxis(1)->SetName("etaPos");
              fTHnSparseRecLegsFromSecondaryPair_tmp->GetAxis(2)->SetName("phiPos");
              fTHnSparseRecLegsFromSecondaryPair_tmp->GetAxis(3)->SetName("ptNeg");
              fTHnSparseRecLegsFromSecondaryPair_tmp->GetAxis(4)->SetName("etaNeg");
              fTHnSparseRecLegsFromSecondaryPair_tmp->GetAxis(5)->SetName("phiNeg");
              fTHnSparseRecLegsFromSecondaryPair_tmp->GetAxis(6)->SetName("opAngle");
              fTHnSparseRecLegsFromSecondaryPair_tmp->SetName(Form("fTHnSparseRecLegsFromPairTest_%s;ptPos;etaPos;phiPos;ptNeg;etaNeg;phiNeg;opAngle", fSecondaryPairMCSignal.at(i).GetName()));
              fTHnSparseRecLegsFromSecondaryPair_tmp->SetTitle(Form("fTHnSparseRecLegsFromPairTest_%s;ptPos;etaPos;phiPos;ptNeg;etaNeg;phiNeg;opAngle", fSecondaryPairMCSignal.at(i).GetName()));
              fTHnSparseRecLegsFromSecondaryPair.push_back(fTHnSparseRecLegsFromSecondaryPair_tmp);
              list1->Add(fTHnSparseRecLegsFromSecondaryPair_tmp);
            }
          }

          fPairList->Add(list1);
        }

      }

    // ######################################################
    // #####################  Resolution #########################
    // ######################################################
    fResolutionList = new TList();
    fResolutionList->SetName("Resolution");
    fResolutionList->SetOwner();

      fPGen_DeltaP                         = new TH2D("PGen_DeltaP",                        "", 1000, 0., 20., fNResolutionDeltaptBins, fResolutionDeltaPtBins.data());
      fPGen_PrecOverPGen                   = new TH2D("PGen_PrecOverPGen",                  "", 1000, 0., 20., fNResolutionRelptBins, fResolutionRelPtBins.data());
      fPtGen_DeltaPt                       = new TH2D("PtGen_DeltaPt",                      "", 1000, 0., 20., fNResolutionDeltaptBins, fResolutionDeltaPtBins.data());
      fPtGen_DeltaPtOverPtGen              = new TH2D("PtGen_DeltaPtOverPtGen",             "", 1000, 0., 20., fNResolutionDeltaptBins, -1., +1.);
      fPtGen_PtRecOverPtGen                = new TH2D("PtGen_PtRecOverPtGen",               "", 1000, 0., 20., fNResolutionRelptBins, fResolutionRelPtBins.data());
      fPtGen_DeltaPt_wGenSmeared           = new TH2D("PtGen_DeltaPt_wGenSmeared",          "", 1000, 0., 20., fNResolutionDeltaptBins, fResolutionDeltaPtBins.data());
      fPtGen_DeltaPtOverPtGen_wGenSmeared  = new TH2D("PtGen_DeltaPtOverPtGen_wGenSmeared", "", 1000, 0., 20., fNResolutionDeltaptBins, -1., +1.);
      fPtGen_PtRecOverPtGen_wGenSmeared    = new TH2D("PtGen_PtRecOverPtGen_wGenSmeared",   "", 1000, 0., 20., fNResolutionRelptBins, fResolutionRelPtBins.data());
      fPGen_DeltaEta                       = new TH2D("PGen_DeltaEta",                      "", 1000, 0., 20., fNResolutionetaBins, fResolutionEtaBins.data());
      fPtGen_DeltaEta                      = new TH2D("PtGen_DeltaEta",                     "", 1000, 0., 20., fNResolutionetaBins, fResolutionEtaBins.data());
      fPGen_DeltaTheta                     = new TH2D("PGen_DeltaTheta",                    "", 1000, 0., 20., fNResolutionthetaBins, fResolutionThetaBins.data());
      fPGen_DeltaPhi_Ele                   = new TH2D("PGen_DeltaPhi_Ele",                  "", 1000, 0., 20., fNResolutionphiBins, fResolutionPhiBins.data());
      fPGen_DeltaPhi_Pos                   = new TH2D("PGen_DeltaPhi_Pos",                  "", 1000, 0., 20., fNResolutionphiBins, fResolutionPhiBins.data());
      fPtGen_DeltaPhi_Ele                  = new TH2D("PtGen_DeltaPhi_Ele",                 "", 1000, 0., 20., fNResolutionphiBins, fResolutionPhiBins.data());
      fPtGen_DeltaPhi_Pos                  = new TH2D("PtGen_DeltaPhi_Pos",                 "", 1000, 0., 20., fNResolutionphiBins, fResolutionPhiBins.data());
      fThetaGen_DeltaTheta                 = new TH2D("ThetaGen_DeltaTheta",                "", 220, -0.1*TMath::Pi(), 1.1*TMath::Pi(), fNResolutionthetaBins, fResolutionThetaBins.data());
      fPhiGen_DeltaPhi                     = new TH2D("PhiGen_DeltaPhi",                    "", 320, -0.1*TMath::Pi(), 2.1*TMath::Pi(), fNResolutionphiBins, fResolutionPhiBins.data());

      fPGen_DeltaP                         ->Sumw2();
      fPGen_DeltaEta                       ->Sumw2();
      fPtGen_DeltaPt                       ->Sumw2();
      fPtGen_DeltaPtOverPtGen              ->Sumw2();
      fPtGen_DeltaPt_wGenSmeared           ->Sumw2();
      fPtGen_DeltaPtOverPtGen_wGenSmeared  ->Sumw2();
      fPtGen_PtRecOverPtGen_wGenSmeared    ->Sumw2();
      fPGen_PrecOverPGen                   ->Sumw2();
      fPtGen_PtRecOverPtGen                ->Sumw2();
      fPtGen_DeltaEta                      ->Sumw2();
      fPGen_DeltaTheta                     ->Sumw2();
      fPGen_DeltaPhi_Ele                   ->Sumw2();
      fPGen_DeltaPhi_Pos                   ->Sumw2();
      fPtGen_DeltaPhi_Ele                  ->Sumw2();
      fPtGen_DeltaPhi_Pos                  ->Sumw2();
      fThetaGen_DeltaTheta                 ->Sumw2();
      fPhiGen_DeltaPhi                     ->Sumw2();

      fPGen_DeltaP                         ->GetXaxis()->SetTitle("p^{gen} (GeV/c)");
      fPGen_DeltaP                         ->GetYaxis()->SetTitle("p^{rec} - p^{gen} (GeV/c)");
      fPGen_PrecOverPGen                   ->GetXaxis()->SetTitle("p^{gen} (GeV/c)");
      fPGen_PrecOverPGen                   ->GetYaxis()->SetTitle("p^{rec} / p^{gen} (GeV/c)");
      fPtGen_DeltaPt                       ->GetXaxis()->SetTitle("p^{gen}_{T} (GeV/c)");
      fPtGen_DeltaPt                       ->GetYaxis()->SetTitle("p^{rec}_{T} - p^{gen}_{T} (GeV/c)");
      fPtGen_DeltaPtOverPtGen              ->GetXaxis()->SetTitle("p^{gen}_{T} (GeV/c)");
      fPtGen_DeltaPtOverPtGen              ->GetYaxis()->SetTitle("(p^{gen}_{T} - p^{rec}_{T}) / p^{gen}_{T} (GeV/c)");
      fPtGen_PtRecOverPtGen                ->GetXaxis()->SetTitle("p^{gen}_{T} (GeV/c)");
      fPtGen_PtRecOverPtGen                ->GetYaxis()->SetTitle("p^{rec}_{T} / p^{gen}_{T} (GeV/c)");
      fPtGen_DeltaPt_wGenSmeared           ->GetXaxis()->SetTitle("p^{gen}_{T} (GeV/c)");
      fPtGen_DeltaPt_wGenSmeared           ->GetYaxis()->SetTitle("p^{gen+smeared}_{T} - p^{gen}_{T} (GeV/c)");
      fPtGen_DeltaPtOverPtGen_wGenSmeared  ->GetXaxis()->SetTitle("p^{gen}_{T} (GeV/c)");
      fPtGen_DeltaPtOverPtGen_wGenSmeared  ->GetYaxis()->SetTitle("(p^{gen}_{T} - p^{gen+smeared}_{T}) / p^{gen}_{T} (GeV/c)");
      fPtGen_PtRecOverPtGen_wGenSmeared    ->GetXaxis()->SetTitle("p^{gen}_{T} (GeV/c)");
      fPtGen_PtRecOverPtGen_wGenSmeared    ->GetYaxis()->SetTitle("p^{gen+smeared}_{T} / p^{gen}_{T} (GeV/c)");
      fPGen_DeltaEta                       ->GetXaxis()->SetTitle("p^{gen} (GeV/c)");
      fPGen_DeltaEta                       ->GetYaxis()->SetTitle("#eta^{rec} - #eta^{gen}");
      fPtGen_DeltaEta                      ->GetXaxis()->SetTitle("p^{gen}_{T} (GeV/c)");
      fPtGen_DeltaEta                      ->GetYaxis()->SetTitle("#eta^{rec} - #eta^{gen}");
      fPGen_DeltaTheta                     ->GetXaxis()->SetTitle("p^{gen} (GeV/c)");
      fPGen_DeltaTheta                     ->GetYaxis()->SetTitle("#theta^{rec} - #theta^{gen} (rad)");
      fPGen_DeltaPhi_Ele                   ->GetXaxis()->SetTitle("p^{gen} (GeV/c)");
      fPGen_DeltaPhi_Ele                   ->GetYaxis()->SetTitle("#varphi^{rec} - #varphi^{gen} (rad)");
      fPtGen_DeltaPhi_Ele                  ->GetXaxis()->SetTitle("p^{gen}_{T} (GeV/c)");
      fPtGen_DeltaPhi_Ele                  ->GetYaxis()->SetTitle("#varphi^{gen} - #varphi^{rec} (rad)");
      fPtGen_DeltaPhi_Pos                  ->GetXaxis()->SetTitle("p^{gen}_{T} (GeV/c)");
      fPtGen_DeltaPhi_Pos                  ->GetYaxis()->SetTitle("#varphi^{gen} - #varphi^{rec} (rad)");
      fPGen_DeltaPhi_Pos                   ->GetXaxis()->SetTitle("p^{gen}_{T} (GeV/c)");
      fPGen_DeltaPhi_Pos                   ->GetYaxis()->SetTitle("#varphi^{rec} - #varphi^{gen} (rad)");
      fThetaGen_DeltaTheta                 ->GetXaxis()->SetTitle("#theta^{gen} (rad)");
      fThetaGen_DeltaTheta                 ->GetYaxis()->SetTitle("#theta^{rec} - #theta^{gen} (rad)");
      fPhiGen_DeltaPhi                     ->GetXaxis()->SetTitle("#varphi^{gen} (rad)");
      fPhiGen_DeltaPhi                     ->GetYaxis()->SetTitle("#varphi^{rec} - #varphi^{gen} (rad)");
      fResolutionList->Add(fPGen_DeltaP);
      fResolutionList->Add(fPGen_PrecOverPGen);
      fResolutionList->Add(fPtGen_DeltaPt);
      fResolutionList->Add(fPtGen_DeltaPtOverPtGen);
      fResolutionList->Add(fPtGen_PtRecOverPtGen);
      fResolutionList->Add(fPtGen_DeltaPt_wGenSmeared);
      fResolutionList->Add(fPtGen_DeltaPtOverPtGen_wGenSmeared);
      fResolutionList->Add(fPtGen_PtRecOverPtGen_wGenSmeared);
      fResolutionList->Add(fPGen_DeltaEta);
      fResolutionList->Add(fPtGen_DeltaEta);
      fResolutionList->Add(fPGen_DeltaTheta);
      fResolutionList->Add(fPGen_DeltaPhi_Ele);
      fResolutionList->Add(fPtGen_DeltaPhi_Ele);
      fResolutionList->Add(fPtGen_DeltaPhi_Pos);
      fResolutionList->Add(fPGen_DeltaPhi_Pos);
      fResolutionList->Add(fThetaGen_DeltaTheta);
      fResolutionList->Add(fPhiGen_DeltaPhi);


      // ######################################################
      // #################  Four Electrons ####################
      // ######################################################
                                                                                // if(fdebug) std::cout << __LINE__ << " DEBUG_AnalysisTask: fDoFourPairing " << fDoFourPairing << std::endl;
      if (fDoFourPairing == true){
        fFourPairList = new TList();
        fFourPairList->SetName("4 el. Pairs");
        fFourPairList->SetOwner();
                                                                                // if(fdebug) std::cout << __LINE__ << " DEBUG_AnalysisTask: Line cout " << std::endl;

        fGeneratedFourPairsList = new TList();
        fGeneratedFourPairsList->SetName("Generated");
        fGeneratedFourPairsList->SetOwner();
        for (unsigned int i = 0; i < fFourPairMCSignal.size(); /*i++*/ i+=2){
          TH2D* th2_tmp = new TH2D(Form("Ngen_%s", fFourPairMCSignal.at(i).GetName()),";m_{eeee};p_{T,eeee}",fNmassBins,fMassBins.data(),fNpairptBins,fPairPtBins.data());
          th2_tmp->Sumw2();
          fHistGenFourPair.push_back(th2_tmp);
          fGeneratedFourPairsList->Add(th2_tmp);
        }
        fFourPairList->Add(fGeneratedFourPairsList);

        fGeneratedSmearedFourPairsList = new TList();
        fGeneratedSmearedFourPairsList->SetName("GeneratedSmeared");
        fGeneratedSmearedFourPairsList->SetOwner();
        for (unsigned int i = 0; i < fFourPairMCSignal.size(); i+=2){
         TH2D* th2_tmp = new TH2D(Form("Ngen_%s", fFourPairMCSignal.at(i).GetName()),";m_{eeee};p_{T,eeee}",fNmassBins,fMassBins.data(),fNpairptBins,fPairPtBins.data());
         th2_tmp->Sumw2();
         fHistGenSmearedFourPair.push_back(th2_tmp);
         fGeneratedSmearedFourPairsList->Add(th2_tmp);
        }

        fFourPairList->Add(fGeneratedSmearedFourPairsList);

        // Generated reconstructed lists for every cutsetting one list and every MCsignal 1 histogram
        for (unsigned int list_i = 0; list_i < fTrackCuts_primary_standard.size(); ++list_i){
          TList* list = new TList();
          list->SetName(fTrackCuts_primary_standard.at(list_i)->GetName());
          list->SetOwner();

          for (unsigned int i = 0; i < fFourPairMCSignal.size(); i+=2){
            TH2D* th2_tmp = new TH2D(Form("Nrec_%s", fFourPairMCSignal.at(i).GetName()),";m_{eeee};p_{T,eeee}",fNmassBins,fMassBins.data(),fNpairptBins,fPairPtBins.data());
            th2_tmp->Sumw2();
            fHistRecFourPair.push_back(th2_tmp);
            list->Add(th2_tmp);
          }
          fFourPairList->Add(list);
        }

        TH2D* fHistBeforePrimSecPrefilter = new TH2D ("HistBeforePrimSecPrefilter",";m_{eeee};p_{T,eeee}",fNmassBins,fMassBins.data(),fNpairptBins,fPairPtBins.data());
        TH2D* fHistBeforeSecSecPrefilter = new TH2D ("HistBeforeSecSecPrefilter",";m_{eeee};p_{T,eeee}",fNmassBins,fMassBins.data(),fNpairptBins,fPairPtBins.data());
        TH2D* fHistAfterPrimSecPrefilter = new TH2D ("HistAfterPrimSecPrefilter",";m_{eeee};p_{T,eeee}",fNmassBins,fMassBins.data(),fNpairptBins,fPairPtBins.data());
        TH2D* fHistAfterSecSecPrefilter = new TH2D ("HistAfterSecSecPrefilter",";m_{eeee};p_{T,eeee}",fNmassBins,fMassBins.data(),fNpairptBins,fPairPtBins.data());
        fHistBeforePrimSecPrefilter->Sumw2();
        fHistBeforeSecSecPrefilter->Sumw2();
        fHistAfterPrimSecPrefilter->Sumw2();
        fHistAfterSecSecPrefilter->Sumw2();
        fVecHistPrefilters.push_back(fHistBeforeSecSecPrefilter);
        fVecHistPrefilters.push_back(fHistBeforePrimSecPrefilter);
        fVecHistPrefilters.push_back(fHistAfterSecSecPrefilter);
        fVecHistPrefilters.push_back(fHistAfterPrimSecPrefilter);
        fFourPairList->Add(fHistBeforePrimSecPrefilter);
        fFourPairList->Add(fHistAfterPrimSecPrefilter);
        fFourPairList->Add(fHistBeforeSecSecPrefilter);
        fFourPairList->Add(fHistAfterSecSecPrefilter);


      } // end if fDoFourPairing

    fOutputList->Add(fSingleElectronList);
    fOutputList->Add(fResolutionList);
    if (fDoPairing)     fOutputList->Add(fPairList);
    if (fDoFourPairing) fOutputList->Add(fFourPairList);

    CreateSupportHistos();
    fOutputList->Add(fOutputListSupportHistos);


  PostData(1, fOutputList);
}


// ############################################################################
// ############################################################################
void AliAnalysisTaskEtaReconstruction::UserExec(Option_t* option){
  const double pi = TMath::Pi();
                                                                                // if(fdebug) std::cout << __LINE__ << " DEBUG_AnalysisTask: Line cout " << std::endl;

  // THIS MIGHT BE USED IN THE FUTURE TO SEPARATE DIFFERENT GENERATORS FROM EACH OTHER
  // AliMCEvent* mcEvent = MCEvent();
  // if(!mcEvent)return 0;
  //
  // for(int it = 0;it < mcEvent->GetNumberOfTracks();++it){
  //
  //   //if ESD AliMCParticle* part = (AliMCParticle*)mcEvent->GetTrack(it);
  //   //if AOD AliVParticle* part=(AliVParticle*)mcEvent->GetTrack(it);
  //


  fGenNegPart_primary.clear();
  fGenPosPart_primary.clear();
  fGenNegPart_secondary.clear();
  fGenPosPart_secondary.clear();
  fGenSmearedNegPart_primary.clear();
  fGenSmearedPosPart_primary.clear();
  fGenSmearedNegPart_secondary.clear();
  fGenSmearedPosPart_secondary.clear();
  fRecNegPart_primary.clear();
  fRecPosPart_primary.clear();
  fRecNegPart_secondary.clear();
  fRecPosPart_secondary.clear();
  fGenPairVec_primary.clear();
  fGenPairVec_secondary.clear();
  fGenSmearedPairVec_primary.clear();
  fGenSmearedPairVec_secondary.clear();
  fRecPairVec_primary.clear();
  fRecPairVec_secondary.clear();  // not in use at the moment (used in DoRecTwoPairing for secondary case)
  fRecV0Pair.clear();

  fPreFilter_BadTracksLabel_primary.clear();

  // ##########################################################
  // Set MC event
  if(!AliDielectronMC::Instance()->ConnectMCEvent()) return;

  // ##########################################################
  // Manage AOD&ESD handling and the corresponding events

  isAOD = false;
  AliInputEventHandler *eventHandler = nullptr;
  AliInputEventHandler *eventHandlerMC = nullptr;
                                                                                // if(fdebug) std::cout << __LINE__ << " DEBUG_AnalysisTask: Line cout " << std::endl;
  if ((AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler())->IsA() == AliAODInputHandler::Class()){
                                                                                // if(fdebug) std::cout << __LINE__ << " DEBUG_AnalysisTask: Line cout " << std::endl;
    isAOD = true;
    eventHandler = dynamic_cast<AliAODInputHandler*> (AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
    eventHandlerMC = eventHandler;
  }
  else if ((AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler())->IsA() == AliESDInputHandler::Class()){
                                                                                // if(fdebug) std::cout << __LINE__ << " DEBUG_AnalysisTask: Line cout " << std::endl;
    isAOD = false;
    // eventHandler = dynamic_cast<AliESDInputHandler*> (AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
    eventHandlerMC = dynamic_cast<AliMCEventHandler*> (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler());
    eventHandler   = dynamic_cast<AliESDInputHandler*> (AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
  }
  //
  fMC = eventHandlerMC->MCEvent();
  if (!fMC) { Printf("ERROR: fMC not available"); return; }

  if (!fPIDResponse) SetPIDResponse( eventHandler->GetPIDResponse() );
  AliDielectronVarManager::SetPIDResponse(fPIDResponse);

  if(fPostPIDCntrdCorrTPC) {
    // std::cout << "TPC mean correction applied" << std::endl;
    AliDielectronPID::SetCentroidCorrFunction(fPostPIDCntrdCorrTPC);
  }
  if(fPostPIDWdthCorrTPC)  {
    // std::cout << "TPC width correction applied" << std::endl;
    AliDielectronPID::SetWidthCorrFunction(fPostPIDWdthCorrTPC);
  }
  if(fPostPIDCntrdCorrITS) {
    // std::cout << "ITS mean correction applied" << std::endl;
    AliDielectronPID::SetCentroidCorrFunctionITS(fPostPIDCntrdCorrITS);
  }
  if(fPostPIDWdthCorrITS)  {
    // std::cout << "ITS width correction applied" << std::endl;
    AliDielectronPID::SetWidthCorrFunctionITS(fPostPIDWdthCorrITS);
  }
  if(fPostPIDCntrdCorrTOF) {
    // std::cout << "TOF mean correction applied" << std::endl;
    AliDielectronPID::SetCentroidCorrFunctionTOF(fPostPIDCntrdCorrTOF);
  }
  if(fPostPIDWdthCorrTOF)  {
    // std::cout << "TOF width correction applied" << std::endl;
    AliDielectronPID::SetWidthCorrFunctionTOF(fPostPIDWdthCorrTOF);
  }


  if (isAOD) fEvent = static_cast<AliAODEvent*>(eventHandler->GetEvent());
  else       fEvent = static_cast<AliESDEvent*>(eventHandler->GetEvent());

  AliDielectronVarManager::SetEvent(fEvent);
                                                                                // if(fdebug) std::cout << __LINE__ << " DEBUG_AnalysisTask: Line cout " << std::endl;
  // ##########################################################
  // All events before all cuts
  fHistEventStat->Fill(kAllEvents);

  // ##########################################################
  // calculating physics selection stuff
  ULong64_t isSelected = AliVEvent::kMB;
  if( fSelectPhysics && !isAOD){
    if((!isAOD && eventHandler->GetEventSelection())){
      isSelected = eventHandler->IsEventSelected();

      isSelected&=fTriggerMask;
    }
  }
                                                                                // if(fdebug) std::cout << __LINE__ << " DEBUG_AnalysisTask: Line cout " << std::endl;
  // ##########################################################
  // Apply physics selection
  if (isSelected==0) return;
  fHistEventStat->Fill(kPhysicsSelectionEvents);

  // ##########################################################
  // Apply event filter
  if (fEventFilter) {
    if (!fEventFilter->IsSelected(fEvent)) return;
  }
  fHistEventStat->Fill(kFilteredEvents);
                                                                                // if(fdebug) std::cout << __LINE__ << " DEBUG_AnalysisTask: Line cout " << std::endl;
  // ##########################################################
  // Monitor z-vertex
  const AliVVertex* vtx = fEvent->GetPrimaryVertex();
  Double_t vtxZGlobal = -99.;
  Int_t nCtrb = -1;
  if (vtx) {
    vtxZGlobal = vtx->GetZ();
    nCtrb = vtx->GetNContributors();
  }
  fHistVertex->Fill(vtxZGlobal);
  fHistVertexContibutors->Fill(nCtrb);
                                                                                // if(fdebug) std::cout << __LINE__ << " DEBUG_AnalysisTask: Line cout " << std::endl;
  // ##########################################################
  // Apply centrality selection
  double centralityF = -1;
  AliMultSelection *multSelection = (AliMultSelection*)fEvent->FindListObject("MultSelection");
  if (multSelection) centralityF  = multSelection->GetMultiplicityPercentile("V0M",kFALSE);
  if (centralityF == -1 && fMaxCentrality == -1 && fMinCentrality == -1) {/*do nothing*/} // is used for pp and pPb analysis
  // else if (centralityF > fMaxCentrality || centralityF < fMinCentrality) { return;} // reject event

  fHistEventStat->Fill(kCentralityEvents);
  fHistEvents->Fill(0.5);
  fHistCentrality->Fill(centralityF);

  // Calculating the weight when centrality correction is applied
  double centralityWeight = 1.;
  // if (fHistCentralityCorrection != 0x0){
  //   centralityWeight = (fHistCentralityCorrection->GetEntries() / fHistCentralityCorrection->GetNbinsX()) / fHistCentralityCorrection->FindBin(centralityF) ;
  //   std::cout << "cent: " << centralityF << "  " << "weight: " << centralityWeight << std::endl;
  // }

  // ##########################################################
  // Fill Multiplicity histogram
  int nTracks = fEvent->GetNumberOfTracks();
  fHistNTracks->Fill(nTracks);
                                                                                // if(fdebug) std::cout << __LINE__ << " DEBUG_AnalysisTask: Line cout " << std::endl;
                                                                                // if(fdebug) std::cout << __LINE__ << " DEBUG_AnalysisTask: fSinglePrimaryLegMCSignal[0]: " << fSinglePrimaryLegMCSignal[0].GetName() << " fSingleSecondaryLegMCSignal[0]: " << fSingleSecondaryLegMCSignal[0].GetName() << std::endl;

  // ######################################################
  // ######################################################
  // ######################################################
  // Start particle loop
                                                                                // if (fdebug) std::cout << "Generated and Generated Smeared Single Particle Loop ... " << std::endl;
  for(int iPart = 0; iPart < fMC->GetNumberOfTracks(); iPart++) {
    AliVParticle* mcPart1  = (AliVParticle*)fMC->GetTrack(iPart);
    AliVParticle* mcMPart1  = (AliVParticle*)fMC->GetTrack(TMath::Abs(mcPart1->GetMother()));
    if (!mcPart1) continue;
    // ##########################################################
    // Checking minimum and maximum values for generated particles
    if (mcPart1->Pt()  < fPtMinGen  || mcPart1->Pt()  > fPtMaxGen)  continue;
    if (mcPart1->Eta() < fEtaMinGen || mcPart1->Eta() > fEtaMaxGen) continue;
    // ##########################################################
    // Check MC signals
    std::vector<Bool_t> mcSignal_acc_primary(fSinglePrimaryLegMCSignal.size(), kFALSE); // initialize vector which stores if track is accepted by [i]-th mcsignal
    std::vector<Bool_t> mcSignal_acc_secondary(fSingleSecondaryLegMCSignal.size(), kFALSE); // initialize vector which stores if track is accepted by [i]-th mcsignal
    CheckSinglePrimaryLegMCsignals(mcSignal_acc_primary, iPart);
    CheckSingleSecondaryLegMCsignals(mcSignal_acc_secondary, iPart);

    // ##########################################################
    // check if at least one mc signal is true
    if (CheckIfOneIsTrue(mcSignal_acc_primary) == kFALSE && CheckIfOneIsTrue(mcSignal_acc_secondary) == kFALSE) continue;
                                                                                // if(fdebug) std::cout << __LINE__ << " DEBUG_AnalysisTask: Line cout " << std::endl;
    // ##########################################################
    // check if correct generator used
    bool generatorForMCSignal  = CheckGenerator(iPart, fGeneratorMCSignalHashs);
    if (!generatorForMCSignal) continue;
    bool generatorForULSSignal = CheckGenerator(iPart, fGeneratorULSSignalHashs);
    if (!generatorForMCSignal && !generatorForULSSignal) continue;
    // if (!CheckGenerator(iPart, fGeneratorHashs)) continue;
                                                                                // if(fdebug) std::cout << __LINE__ << " DEBUG_AnalysisTask: Line cout " << std::endl;

    // ##########################################################
    // Creating particles to summarize all the data
    int motherID = TMath::Abs(mcPart1->GetMother());
    int grandmotherID = TMath::Abs(mcMPart1->GetMother());
    Particle part = CreateParticle(mcPart1);
    part.isMCSignal_primary = mcSignal_acc_primary;
    part.isMCSignal_secondary = mcSignal_acc_secondary;
    part.SetTrackID(iPart);
    part.SetMotherID(motherID);
    part.SetGrandMotherID(grandmotherID);
    part.SetULSSignalPair(generatorForULSSignal);
    part.SetMCSignalPair(generatorForMCSignal);


                                                                                // if(fdebug) std::cout << __LINE__ << " DEBUG_AnalysisTask: Line cout " << std::endl;
    // ##########################################################
    // Filling generated particle histograms according to MCSignals and push particles into vectors for later pairing
    if ((part.fPt < fPtMin || part.fPt > fPtMax || part.fEta < fEtaMin || part.fEta > fEtaMax) == kFALSE) {                      // Added kinematic cuts for single generated electrons
      for (unsigned int i = 0; i < part.isMCSignal_primary.size(); ++i){
        if (part.isMCSignal_primary[i]) {
          if      (part.fCharge < 0){
            dynamic_cast<TH3D*>(fHistGenPrimaryNegPart.at(i))->Fill(part.fPt, part.fEta, part.fPhi, centralityWeight);
          }
          else if (part.fCharge > 0) {
            dynamic_cast<TH3D*>(fHistGenPrimaryPosPart.at(i))->Fill(part.fPt, part.fEta, part.fPhi, centralityWeight);
          }
        }
      }
      for (unsigned int i = 0; i < part.isMCSignal_secondary.size(); ++i){
        if (part.isMCSignal_secondary[i]) {
          if      (part.fCharge < 0){
            dynamic_cast<TH3D*>(fHistGenSecondaryNegPart.at(i))->Fill(part.fPt, part.fEta, part.fPhi, centralityWeight);
          }
          else if (part.fCharge > 0) {
            dynamic_cast<TH3D*>(fHistGenSecondaryPosPart.at(i))->Fill(part.fPt, part.fEta, part.fPhi, centralityWeight);
          }
        }
      }

      if (fDoPairing == true || fDoFourPairing == true) {
        if      (part.fCharge < 0 && CheckIfOneIsTrue(mcSignal_acc_primary) == kTRUE /*&& part.isMCSignal_primary[0]*/) fGenNegPart_primary.push_back(part);     // store particles for later pairing
        if      (part.fCharge > 0 && CheckIfOneIsTrue(mcSignal_acc_primary) == kTRUE /*&& part.isMCSignal_primary[0]*/) fGenPosPart_primary.push_back(part);     // store particles for later pairing
        if      (part.fCharge < 0 && CheckIfOneIsTrue(mcSignal_acc_secondary) == kTRUE /*&& part.isMCSignal_secondary[0]*/) fGenNegPart_secondary.push_back(part); // store particles for later pairing
        if      (part.fCharge > 0 && CheckIfOneIsTrue(mcSignal_acc_secondary) == kTRUE /*&& part.isMCSignal_secondary[0]*/) fGenPosPart_secondary.push_back(part); // store particles for later pairing
      }
    }
                                                                                // if(fdebug) std::cout << __LINE__ << " DEBUG_AnalysisTask: Line cout " << std::endl;
    // ##########################################################
    // Filling generated smeared particle histograms according to MCSignals
    // and separated into pos and neg charge as well as primary and secondary
    if (fArrResoPt){ // Smear particles to fill "GeneratedSmeared"
      TLorentzVector smearedVec = ApplyResolution(part.fPt, part.fEta, part.fPhi, part.fCharge);
      part.fPt_smeared  = smearedVec.Pt();
      part.fEta_smeared = smearedVec.Eta();

      if (smearedVec.Phi() < 0) part.fPhi_smeared = smearedVec.Phi()+ 2 * pi;
      else part.fPhi_smeared = smearedVec.Phi();

      if ((part.fPt_smeared < fPtMin || part.fPt_smeared > fPtMax || part.fEta_smeared < fEtaMin || part.fEta_smeared > fEtaMax) == kFALSE) { // Added kinematic cuts for single generated smeared electrons

        for (unsigned int i = 0; i < part.isMCSignal_primary.size(); ++i){
          if (part.isMCSignal_primary[i]) {
            if      (part.fCharge < 0){
              dynamic_cast<TH3D*>(fHistGenSmearedPrimaryNegPart.at(i))->Fill(part.fPt_smeared, part.fEta_smeared, part.fPhi_smeared , centralityWeight);
            }
            else if (part.fCharge > 0) {
              dynamic_cast<TH3D*>(fHistGenSmearedPrimaryPosPart.at(i))->Fill(part.fPt_smeared, part.fEta_smeared, part.fPhi_smeared , centralityWeight);
            }
          }
        }
        for (unsigned int i = 0; i < part.isMCSignal_secondary.size(); ++i){
          if (part.isMCSignal_secondary[i]) {
            if      (part.fCharge < 0){
              dynamic_cast<TH3D*>(fHistGenSmearedSecondaryNegPart.at(i))->Fill(part.fPt_smeared, part.fEta_smeared, part.fPhi_smeared , centralityWeight);
            }
            else if (part.fCharge > 0) {
              dynamic_cast<TH3D*>(fHistGenSmearedSecondaryPosPart.at(i))->Fill(part.fPt_smeared, part.fEta_smeared, part.fPhi_smeared , centralityWeight);
            }
          }
        }

        if (fDoPairing == true || fDoFourPairing == true) {
          if      (part.fCharge < 0 && CheckIfOneIsTrue(mcSignal_acc_primary) == kTRUE /*&& part.isMCSignal_primary[0]*/) fGenSmearedNegPart_primary.push_back(part);      // store particles for later pairing
          if      (part.fCharge > 0 && CheckIfOneIsTrue(mcSignal_acc_primary) == kTRUE /*&& part.isMCSignal_primary[0]*/) fGenSmearedPosPart_primary.push_back(part);      // store particles for later pairing
          if      (part.fCharge < 0 && CheckIfOneIsTrue(mcSignal_acc_secondary) == kTRUE /*&& part.isMCSignal_secondary[0]*/) fGenSmearedNegPart_secondary.push_back(part);  // store particles for later pairing
          if      (part.fCharge > 0 && CheckIfOneIsTrue(mcSignal_acc_secondary) == kTRUE /*&& part.isMCSignal_secondary[0]*/) fGenSmearedPosPart_secondary.push_back(part);  // store particles for later pairing
        }
      }
    }

  }// end of MC track loop


                                                                                // if (fdebug) std::cout << __LINE__ << " DEBUG_AnalysisTask: fGenNegPart_primary = "    << fGenNegPart_primary.size()    << " fGenPosPart_primary = "    << fGenPosPart_primary.size()    << std::endl;
                                                                                // if (fdebug) std::cout << __LINE__ << " DEBUG_AnalysisTask: fGenNegPart_secondary = "  << fGenNegPart_secondary.size()  << " fGenPosPart_secondary = "  << fGenPosPart_secondary.size()  << std::endl << std::endl;
                                                                                // if (fdebug) std::cout << __LINE__ << " DEBUG_AnalysisTask: fGenSmearedNegPart_primary = "    << fGenSmearedNegPart_primary.size()    << " fGenSmearedPosPart_primary = "    << fGenSmearedPosPart_primary.size()    << std::endl;
                                                                                // if (fdebug) std::cout << __LINE__ << " DEBUG_AnalysisTask: fGenSmearedNegPart_secondary = "  << fGenSmearedNegPart_secondary.size()  << " fGenSmearedPosPart_secondary = "  << fGenSmearedPosPart_secondary.size()  << std::endl << std::endl;
  // ##########################################################
  // ##########################################################
  // ##########################################################
  // Start reconstructed track Loop
                                                                                // if (fdebug) std::cout << "Reconstructed Single Particle Loop ... " << std::endl;
  for (Int_t iTracks = 0; iTracks < fEvent->GetNumberOfTracks(); iTracks++){

    // ##########################################################
    // Track handling
    AliVParticle* track = fEvent->GetTrack(iTracks);
    if (!track) { Printf("ERROR: Could not receive track %d", iTracks); continue; }
    if (isAOD) track = static_cast<AliAODTrack*>(track);
    else       track = static_cast<AliESDtrack*>(track);
    int label = track->GetLabel();
    int abslabel = TMath::Abs(label);

    // ##########################################################
    // Apply MC signals
    std::vector<Bool_t> mcSignal_acc_primary(fSinglePrimaryLegMCSignal.size(), kFALSE); // vector which stores if track is accepted by [i]-th mcsignal
    std::vector<Bool_t> mcSignal_acc_secondary(fSingleSecondaryLegMCSignal.size(), kFALSE); // vector which stores if track is accepted by [i]-th mcsignal
    CheckSinglePrimaryLegMCsignals(mcSignal_acc_primary, abslabel);
    CheckSingleSecondaryLegMCsignals(mcSignal_acc_secondary, abslabel);

    // ##########################################################
    // check if at least one mc signal is true otherwise skip this particle
    if (CheckIfOneIsTrue(mcSignal_acc_primary) == kFALSE && CheckIfOneIsTrue(mcSignal_acc_secondary) == kFALSE) continue;

    // ##########################################################
    // check if correct generator used
    bool generatorForMCSignal  = CheckGenerator(label, fGeneratorMCSignalHashs);
    bool generatorForULSSignal = CheckGenerator(label, fGeneratorULSSignalHashs);
    // std::cout << "generatorForMCSignal = " << generatorForMCSignal << std::endl;
    // std::cout << "generatorForULSSignal = " << generatorForULSSignal << std::endl;
    if (!generatorForMCSignal && !generatorForULSSignal) continue;
    // if (!CheckGenerator(label, fGeneratorHashs)) continue;

    // ##########################################################
                                                                                // // changed in order to disable reconstructed cuts
                                                                                // std::vector<bool> selected_primary(fTrackCuts_primary_PreFilter.size(), kTRUE); // vector which stores if track is accepted by [i]-th selection cut
                                                                                // std::vector<bool> selected_secondary(fPairCuts_secondary_standard.size(), kTRUE); // vector which stores if track is accepted by [i]-th selection cut
    // Check if particle is passing primary selection cuts
    std::vector<bool> selected_primary(fTrackCuts_primary_PreFilter.size(), kFALSE); // vector which stores if track is accepted by [i]-th selection cut
    for (UInt_t iCut=0; iCut<fTrackCuts_primary_PreFilter.size(); ++iCut){ // loop over all specified cutInstances
      UInt_t selectedMask_primary=( 1 << fTrackCuts_primary_PreFilter.at(iCut)->GetCuts()->GetEntries())-1;
      // cutting logic taken from AliDielectron::FillTrackArrays()
      // apply track cuts
                                                                                // if(fdebug) std::cout << __LINE__ << " DEBUG_AnalysisTask: Cout line" << std::endl;
      UInt_t cutMask_primary = fTrackCuts_primary_PreFilter.at(iCut)->IsSelected(track);
      if (cutMask_primary == selectedMask_primary) {selected_primary[iCut] = kTRUE; /*std::cout << "prim_reconstructed TRUE" << std::endl;*/}
    }

    // // Check if particle is passing secondary selection cuts
    std::vector<bool> selected_secondary(fPairCuts_secondary_standard.size(), kFALSE); // vector which stores if track is accepted by [i]-th selection cut
    for (UInt_t iCut=0; iCut<fPairCuts_secondary_standard.size(); ++iCut){ // loop over all specified cutInstances
      UInt_t selectedMask_secondary=( 1 << fPairCuts_secondary_standard.at(iCut)->GetCuts()->GetEntries())-1;
                                                                                // if(fdebug) std::cout << __LINE__ << " DEBUG_AnalysisTask: fPairCuts_secondary_standard.at(iCut)->GetCuts()->GetEntries()): " <<  fPairCuts_secondary_standard.at(iCut)->GetCuts()->GetEntries() << ", (1<<X)-1: " <<  ((1 << fPairCuts_secondary_standard.at(iCut)->GetCuts()->GetEntries())-1) << std::endl;
      // cutting logic taken from AliDielectron::FillTrackArrays()
      // apply track cuts
                                                                                // if(fdebug) std::cout << __LINE__ << " DEBUG_AnalysisTask: Cout line" << std::endl;
      UInt_t cutMask_secondary = fPairCuts_secondary_standard.at(iCut)->IsSelected(track);
                                                                                // if(fdebug) std::cout << __LINE__ << " DEBUG_AnalysisTask: cutMask_secondary: " <<  cutMask_secondary << std::endl;
      if (cutMask_secondary == selectedMask_secondary) {selected_secondary[iCut] = kTRUE;/* std::cout << "sec_reconstructed TRUE" << std::endl;*/}
    }

    // ##########################################################
    // check if at least one is selected by cuts otherwise skip this particle
    if (CheckIfOneIsTrue(selected_primary) == kFALSE /*&& CheckIfOneIsTrue(selected_secondary) == kFALSE*/) continue;

    // ##########################################################
    // Create summary particle from track info
    int motherID = TMath::Abs(fMC->GetTrack(abslabel)->GetMother());
    Particle part  = CreateParticle(track);
    part.isMCSignal_primary = mcSignal_acc_primary;
    part.isMCSignal_secondary = mcSignal_acc_secondary;
    part.isReconstructed_primary = selected_primary;
    // part.isReconstructed_secondary = selected_secondary;
    part.SetTrackID(iTracks);
    part.SetMotherID(motherID);
    part.SetULSSignalPair(generatorForULSSignal);
    part.SetMCSignalPair(generatorForMCSignal);


    // ##########################################################
    if (fDoPairing == true || fDoFourPairing == true){
      if      (part.fCharge <  0 && CheckIfOneIsTrue(selected_primary) == kTRUE /*&& CheckIfOneIsTrue(selected_secondary) == kFALSE*/) fRecNegPart_primary.push_back(part);
      if      (part.fCharge >  0 && CheckIfOneIsTrue(selected_primary) == kTRUE /*&& CheckIfOneIsTrue(selected_secondary) == kFALSE*/) fRecPosPart_primary.push_back(part);
      // if      (part.fCharge <  0 && CheckIfOneIsTrue(selected_secondary) == kTRUE /*&& CheckIfOneIsTrue(selected_primary) == kFALSE*/) fRecNegPart_secondary.push_back(part);
      // if      (part.fCharge >  0 && CheckIfOneIsTrue(selected_secondary) == kTRUE /*&& CheckIfOneIsTrue(selected_primary) == kFALSE*/) fRecPosPart_secondary.push_back(part);
    }

    // // Check if generated smeared looks the same like reconstructed
    // if (fArrResoPt){
    //   AliVParticle* genTrack = fMC->GetTrack(abslabel);
    //   double pt_temp  = genTrack->Pt();
    //   double eta_temp = genTrack->Eta();
    //   double phi_temp = genTrack->Phi();
    //
    //   TLorentzVector smearedVec = ApplyResolution(pt_temp, eta_temp, phi_temp, part.fCharge);
    //   part.fPt_smeared  = smearedVec.Pt();
    //   part.fEta_smeared = smearedVec.Eta();
    //   // part.fPhi_smeared = smearedVec.Phi();
    //   if (smearedVec.Phi() < 0) part.fPhi_smeared = smearedVec.Phi() + 2 * pi;
    //   else part.fPhi_smeared = smearedVec.Phi();
    //
    //   // std::cout << "pt_rec: " << part.fPt << "  pt_gen: " << pt_temp << "  phi_rec:" << part.fPhi << "  phi_gen:" << phi_temp << "  phi_gen_smeared:" << part.fPhi_smeared << "  pdgCode:" << genTrack->PdgCode() << std::endl;
    //
    //   // if      (part.fCharge <  0 && CheckIfOneIsTrue(selected_primary) == kTRUE /*&& CheckIfOneIsTrue(selected_secondary) == kFALSE*/) { if(fdebug) std::cout   << __LINE__ << " DEBUG_AnalysisTask: Neg_Prim Particle Label = " << label << ", Particle rec_Pt = " << part.fPt << ", gen smeared Pt = " << part.fPt_smeared << ", Particle rec_Eta = " << part.fEta << ", gen smeared Eta = " << part.fEta_smeared << ", Particle rec_Phi = " << part.fPt << ", gen smeared Phi = " << part.fPt_smeared << std::endl;}
    //   // if      (part.fCharge >  0 && CheckIfOneIsTrue(selected_primary) == kTRUE /*&& CheckIfOneIsTrue(selected_secondary) == kFALSE*/) { if(fdebug) std::cout   << __LINE__ << " DEBUG_AnalysisTask: Pos_Prim Particle Label = " << label << ", Particle rec_Pt = " << part.fPt << ", gen smeared Pt = " << part.fPt_smeared << ", Particle rec_Eta = " << part.fEta << ", gen smeared Eta = " << part.fEta_smeared << ", Particle rec_Phi = " << part.fPt << ", gen smeared Phi = " << part.fPt_smeared << std::endl;}
    //   // if      (part.fCharge <  0 && CheckIfOneIsTrue(selected_secondary) == kTRUE /*&& CheckIfOneIsTrue(selected_primary) == kFALSE*/) { if(fdebug) std::cout   << __LINE__ << " DEBUG_AnalysisTask: Neg_Sec Particle Label = "  << label << ", Particle rec_Pt = " << part.fPt << ", gen smeared Pt = " << part.fPt_smeared << ", Particle rec_Eta = " << part.fEta << ", gen smeared Eta = " << part.fEta_smeared << ", Particle rec_Phi = " << part.fPt << ", gen smeared Phi = " << part.fPt_smeared << std::endl;}
    //   // if      (part.fCharge >  0 && CheckIfOneIsTrue(selected_secondary) == kTRUE /*&& CheckIfOneIsTrue(selected_primary) == kFALSE*/) { if(fdebug) std::cout   << __LINE__ << " DEBUG_AnalysisTask: Pos_Sec Particle Label = "  << label << ", Particle rec_Pt = " << part.fPt << ", gen smeared Pt = " << part.fPt_smeared << ", Particle rec_Eta = " << part.fEta << ", gen smeared Eta = " << part.fEta_smeared << ", Particle rec_Phi = " << part.fPt << ", gen smeared Phi = " << part.fPt_smeared << std::endl;}
    //
    // }

    // // ##########################################################
    // // Filling primary reconstructed particle histograms according to MCSignals
    // for (unsigned int i = 0; i < part.isMCSignal_primary.size(); ++i){
    //   for (unsigned int j = 0; j < part.isReconstructed_primary.size(); ++j){
    //     if (part.isMCSignal_primary[i] == kTRUE) {
    //       if (part.isReconstructed_primary[j] == kTRUE){
    //         if      (part.fCharge < 0) {
    //           dynamic_cast<TH3D*>(fHistRecPrimaryNegPart.at(j * part.isMCSignal_primary.size() + i))->Fill(part.fPt, part.fEta, part.fPhi, centralityWeight);
    //         }
    //         else if (part.fCharge > 0) {
    //           dynamic_cast<TH3D*>(fHistRecPrimaryPosPart.at(j * part.isMCSignal_primary.size() + i))->Fill(part.fPt, part.fEta, part.fPhi, centralityWeight);
    //         }
    //       }// is selected by cutsetting
    //     } // is selected by MC signal
    //   } // end of loop over all cutsettings
    // } // end of loop over all MCsignals
    // // Filling secondary reconstructed particle histograms according to MCSignals
    // for (unsigned int i = 0; i < part.isMCSignal_secondary.size(); ++i){
    //   for (unsigned int j = 0; j < part.isReconstructed_secondary.size(); ++j){
    //     if (part.isMCSignal_secondary[i] == kTRUE) {
    //       if (part.isReconstructed_secondary[j] == kTRUE){
    //         if      (part.fCharge < 0) {
    //           dynamic_cast<TH3D*>(fHistRecSecondaryNegPart.at(j * part.isMCSignal_secondary.size() + i))->Fill(part.fPt, part.fEta, part.fPhi, centralityWeight);
    //         }
    //         else if (part.fCharge > 0) {
    //           dynamic_cast<TH3D*>(fHistRecSecondaryPosPart.at(j * part.isMCSignal_secondary.size() + i))->Fill(part.fPt, part.fEta, part.fPhi, centralityWeight);
    //         }
    //       }// is selected by cutsetting
    //     } // is selected by MC signal
    //   } // end of loop over all cutsettings
    // } // end of loop over all MCsignals
    // // ##########################################################
    // // Fill support histograms with first cutsetting and first mcsignal
    // // if(part.isMCSignal_primary[fSupportMCSignal] == true && part.isReconstructed_primary[fSupportCutsetting] == kTRUE){
    // AliVParticle* mcPart1 = fMC->GetTrack(abslabel);
    // for (unsigned int i = 0; i < part.isReconstructed_primary.size(); i++) {
    //   if(part.isReconstructed_primary[i] == kTRUE){    // check if part is reconstructed in each primary cut setting
    //     for (unsigned int j = 0; j < part.isMCSignal_primary.size(); j++) {
    //       if (part.isMCSignal_primary[j] == kTRUE) {
    //         FillTrackHistograms_Primary(track, mcPart1, i, j); // Fill primary support histograms
    //       }
    //     }
    //   }
    // }
    // // for (unsigned int i = 0; i < part.isReconstructed_secondary.size(); i++) {
    // //   if(part.isReconstructed_secondary[i] == kTRUE){
    // //     for (unsigned int j = 0; j < part.isMCSignal_secondary.size(); j++) {
    // //       if (part.isMCSignal_secondary[j] == kTRUE) {
    // //         FillTrackHistograms_Secondary(track, mcPart1, i, j); // Fill secondary support histograms
    // //       }
    // //     }
    // //   }
    // // }

    if (part.isMCSignal_primary[fSupportMCSignal] == true && part.isReconstructed_primary[fSupportCutsetting] == kTRUE) {
      // ##########################################################
      // Fill resolution histograms
      // ##########################################################
      AliVParticle* mcPart1 = fMC->GetTrack(abslabel);
      double mcP     = mcPart1->P();
      double mcPt    = mcPart1->Pt();
      double mcEta   = mcPart1->Eta();
      double mcPhi   = mcPart1->Phi();
      double mcTheta = mcPart1->Theta();
      double P     = track->P();
      double Pt    = part.fPt;
      double Phi   = part.fPhi;
      double Eta   = part.fEta;
      double Theta = track->Theta();

      if(TMath::Abs(mcEta) < 1.0) {
        // fPGen                ->Fill(mcP);
        // fPRec                ->Fill(recP);
        fPGen_DeltaP           ->Fill(mcP,  P - mcP, centralityWeight);
        fPtGen_DeltaPt         ->Fill(mcPt, Pt - mcPt, centralityWeight);
        fPtGen_DeltaPtOverPtGen->Fill(mcPt, (mcPt - Pt) / mcPt, centralityWeight);
        fPtGen_PtRecOverPtGen  ->Fill(mcPt, Pt / mcPt, centralityWeight);

        if (fArrResoPt){
          double Pt_genSmeared = part.fPt_smeared;
          fPtGen_DeltaPt_wGenSmeared         ->Fill(mcPt, Pt_genSmeared - mcPt, centralityWeight);
          fPtGen_DeltaPtOverPtGen_wGenSmeared->Fill(mcPt, (mcPt - Pt_genSmeared) / mcPt, centralityWeight);
          fPtGen_PtRecOverPtGen_wGenSmeared  ->Fill(mcPt, Pt_genSmeared / mcPt, centralityWeight);
        }

        fPGen_PrecOverPGen     ->Fill(mcP,  P / mcP, centralityWeight);
        if (part.fCharge<0) {
          fPGen_DeltaPhi_Ele   ->Fill(mcP,  Phi - mcPhi, centralityWeight);
          fPtGen_DeltaPhi_Ele  ->Fill(mcPt, mcPhi - Phi, centralityWeight);
        }
        else {
          fPGen_DeltaPhi_Pos   ->Fill(mcP,  Phi - mcPhi, centralityWeight);
          fPtGen_DeltaPhi_Pos  ->Fill(mcPt, mcPhi - Phi, centralityWeight);
        }
        fPhiGen_DeltaPhi       ->Fill(mcPhi, Phi - mcPhi, centralityWeight);
      }

      if(mcPt > fPtMinGen){
        fPGen_DeltaEta       ->Fill(mcP,     Eta - mcEta, centralityWeight);
        fPtGen_DeltaEta      ->Fill(mcPt,    Eta - mcEta, centralityWeight);
        fPGen_DeltaTheta     ->Fill(mcP,     Theta - mcTheta, centralityWeight);
        fThetaGen_DeltaTheta ->Fill(mcTheta, Theta - mcTheta, centralityWeight);
      }
    }
  }

                                                                                // if(fdebug) std::cout << __LINE__ << " DEBUG_AnalysisTask: fTrackCuts_primary_PreFilter.size: " << fTrackCuts_primary_PreFilter.size() << std::endl;
                                                                                // if(fdebug){
                                                                                //   for (size_t i = 0; i < fTrackCuts_primary_PreFilter.size(); i++) {
                                                                                //   std::cout << __LINE__ << " DEBUG_AnalysisTask: " << i+1 << "-th CutSetting: " << fTrackCuts_primary_PreFilter.at(i)->GetName() << std::endl;
                                                                                //   const AliDielectronCutGroup* cgPIDCutsAna = dynamic_cast< const AliDielectronCutGroup*> (fTrackCuts_primary_PreFilter.at(i)->GetCuts()->At(0));
                                                                                //   const AliDielectronCutGroup* cgTrackCutsAnaSPDfirst = dynamic_cast< const AliDielectronCutGroup*> (cgPIDCutsAna->GetCut(3));
                                                                                //   const AliDielectronVarCuts*  trackCutsAOD = dynamic_cast< const AliDielectronVarCuts*> (cgTrackCutsAnaSPDfirst->GetCut(1));
                                                                                //
                                                                                //   // std::cout << trackCutsAOD->GetNCuts() << std::endl;
                                                                                //   for (Int_t j = 0; j < trackCutsAOD->GetNCuts(); j++) {
                                                                                //              std::cout << __LINE__ << " DEBUG_AnalysisTask: Cut "<< j+1 <<" " <<  trackCutsAOD->GetCutName(j) << std::endl;
                                                                                //    }
                                                                                //   }
                                                                                // }
                                                                                // if(fdebug) std::cout << __LINE__ << " DEBUG_AnalysisTask: fPairCuts_secondary_standard.size: " << fPairCuts_secondary_standard.size() << std::endl;
                                                                                // if(fdebug){
                                                                                //   for (size_t i = 0; i < fPairCuts_secondary_standard.size(); i++) {
                                                                                //   std::cout << __LINE__ << " DEBUG_AnalysisTask: " << i+1 << "-th CutSetting: " << fPairCuts_secondary_standard.at(i)->GetName() << std::endl;
                                                                                //   const AliDielectronCutGroup* cgPIDCutsAna = dynamic_cast< const AliDielectronCutGroup*> (fPairCuts_secondary_standard.at(i)->GetCuts()->At(0));
                                                                                //   const AliDielectronCutGroup* cgTrackCutsAnaSPDfirst = dynamic_cast< const AliDielectronCutGroup*> (cgPIDCutsAna->GetCut(3));
                                                                                //   const AliDielectronVarCuts*  trackCutsAOD = dynamic_cast< const AliDielectronVarCuts*> (cgTrackCutsAnaSPDfirst->GetCut(1));
                                                                                //
                                                                                //   // std::cout << trackCutsAOD->GetNCuts() << std::endl;
                                                                                //   for (Int_t j = 0; j < trackCutsAOD->GetNCuts(); j++) {
                                                                                //              std::cout << __LINE__ << " DEBUG_AnalysisTask: Cut "<< j+1 <<" " <<  trackCutsAOD->GetCutName(j) << std::endl;
                                                                                //    }
                                                                                //   }
                                                                                // }
                                                                                // if(fdebug) std::cout << __LINE__ << " DEBUG_AnalysisTask: fRecNegPart_primary:    " << fRecNegPart_primary.size() <<    " fRecPosPart_primary:    " << fRecPosPart_primary.size()   << std::endl;
                                                                                // if(fdebug) std::cout << __LINE__ << " DEBUG_AnalysisTask: fRecNegPart_secondary:  " << fRecNegPart_secondary.size() <<  " fRecPosPart_secondary:  " << fRecPosPart_secondary.size() << std::endl;


  // ##########################################################
  // ##########################################################
  // ##########################################################
  // DO PAIRING
  // ##########################################################

  // float ptPos       = -999;
  // float etaPos      = -999;
  // float phiPos      = -999;
  // float ptNeg       = -999;
  // float etaNeg      = -999;
  // float phiNeg      = -999;
  // float op_angle    = -999;


  if (fDoPairing){
    // if(fdebug) std::cout << "start two pairing" << std::endl;
                                                                                // Debug code in order to check Pdg Codes of Particle, Mother and GrandMother
                                                                                // if (fdebug) {
                                                                                //   if (fdebug) std::cout << __LINE__ << " Primary Particles: " << std::endl;
                                                                                //   for (size_t neg_i = 0; neg_i < fGenNegPart_primary.size(); neg_i++) {
                                                                                //     if (fdebug) std::cout << __LINE__ << " Neg PartID = " << fMC->GetTrack(fGenNegPart_primary[neg_i].GetTrackID())->GetLabel() <<
                                                                                //      ", PDGCode = " << fMC->GetTrack(fGenNegPart_primary[neg_i].GetTrackID())->PdgCode() <<
                                                                                //      ", Mother PdgCode = " << fMC->GetTrack(fGenNegPart_primary[neg_i].GetMotherID())->PdgCode() <<
                                                                                //      /*", GrandMother PdgCode = " << fMC->GetTrack(fGenNegPart_primary[neg_i].GetGrandMotherID())->PdgCode()<<*/ std::endl;
                                                                                //   }
                                                                                //   for (size_t pos_i = 0; pos_i < fGenPosPart_primary.size(); pos_i++) {
                                                                                //     if (fdebug) std::cout << __LINE__ << " Pos PartID = " << fMC->GetTrack(fGenPosPart_primary[pos_i].GetTrackID())->GetLabel() <<
                                                                                //     ", PDGCode = " << fMC->GetTrack(fGenPosPart_primary[pos_i].GetTrackID())->PdgCode() <<
                                                                                //     ", Mother PdgCode = " << fMC->GetTrack(fGenPosPart_primary[pos_i].GetMotherID())->PdgCode() <<
                                                                                //     /*", GrandMother PdgCode = " << fMC->GetTrack(fGenPosPart_primary[pos_i].GetGrandMotherID())->PdgCode()<<*/ std::endl;
                                                                                //   }
                                                                                //
                                                                                //   if (fdebug) std::cout << __LINE__ << " Secondary Particles: " << std::endl;
                                                                                //   for (size_t neg_i = 0; neg_i < fGenNegPart_secondary.size(); neg_i++) {
                                                                                //     if (fdebug) std::cout << __LINE__ << " Neg PartID = " << fMC->GetTrack(fGenNegPart_secondary[neg_i].GetTrackID())->GetLabel() <<
                                                                                //      ", PDGCode = " << fMC->GetTrack(fGenNegPart_secondary[neg_i].GetTrackID())->PdgCode() <<
                                                                                //      ", Mother PdgCode = " << fMC->GetTrack(fGenNegPart_secondary[neg_i].GetMotherID())->PdgCode() <<
                                                                                //      ", GrandMother PdgCode = " << fMC->GetTrack(fGenNegPart_secondary[neg_i].GetGrandMotherID())->PdgCode()<< std::endl;
                                                                                //   }
                                                                                //   for (size_t pos_i = 0; pos_i < fGenPosPart_secondary.size(); pos_i++) {
                                                                                //     if (fdebug) std::cout << __LINE__ << " Pos PartID = " << fMC->GetTrack(fGenPosPart_secondary[pos_i].GetTrackID())->GetLabel() <<
                                                                                //     ", PDGCode = " << fMC->GetTrack(fGenPosPart_secondary[pos_i].GetTrackID())->PdgCode() <<
                                                                                //     ", Mother PdgCode = " << fMC->GetTrack(fGenPosPart_secondary[pos_i].GetMotherID())->PdgCode() <<
                                                                                //     ", GrandMother PdgCode = " << fMC->GetTrack(fGenPosPart_secondary[pos_i].GetGrandMotherID())->PdgCode()<< std::endl;
                                                                                //   }
                                                                                // }

    Bool_t PrimaryPair = kTRUE;
    Bool_t SmearedPair  = kTRUE;
                                                                                // if(fdebug) std::cout << __LINE__ <<  " DEBUG_AnalysisTask: fGenNegPart_primary:   "  << fGenNegPart_primary.size()   << " fGenPosPart_primary:   " <<  fGenPosPart_primary.size()   << std::endl;
                                                                                // if(fdebug) std::cout << __LINE__ <<  " DEBUG_AnalysisTask: fGenNegPart_secondary: "  << fGenNegPart_secondary.size() << " fGenPosPart_secondary: " <<  fGenPosPart_secondary.size() << std::endl;
                                                                                // if (fdebug) std::cout << "Do primary two generated pairing" << std::endl;
    DoGenAndGenSmearTwoPairing(&fGenNegPart_primary, &fGenPosPart_primary, PrimaryPair, !SmearedPair, centralityWeight);
                                                                                // if (fdebug) std::cout << "Do secondary two generated pairing" << std::endl;
    DoGenAndGenSmearTwoPairing(&fGenNegPart_secondary, &fGenPosPart_secondary, !PrimaryPair, !SmearedPair, centralityWeight);

    if(fArrResoPt){
                                                                                // if (fdebug) std::cout << "Do primary two generated smeared pairing" << std::endl;
      DoGenAndGenSmearTwoPairing(&fGenSmearedNegPart_primary, &fGenSmearedPosPart_primary, PrimaryPair, SmearedPair, centralityWeight);
                                                                                // if (fdebug) std::cout << "Do secondary two generated smeared pairing" << std::endl;
      DoGenAndGenSmearTwoPairing(&fGenSmearedNegPart_secondary, &fGenSmearedPosPart_secondary, !PrimaryPair, SmearedPair, centralityWeight);
    }
                                                                                // if (fdebug) std::cout << "Do primary two reconstructed pairing" << std::endl;
    DoRecTwoPairing(fRecNegPart_primary, fRecPosPart_primary, fPrimaryPairMCSignal,  PrimaryPair, centralityWeight); // ! Note: use for secondaries is commented out, since secondaries are selected with V0-Finder
                                                                                // if (fdebug) std::cout << "Do secondary two reconstructed pairing" << std::endl;
    // DoRecTwoPairing(fRecNegPart_secondary, fRecPosPart_secondary, fSecondaryPairMCSignal, !PrimaryPair, centralityWeight);
    DoRecTwoPairingV0(fSecondaryPairMCSignal);

                                                                                // if (fdebug) std::cout << __LINE__ << " DEBUG_AnalysisTask: Size of Vectors after TwoPairing, " << std::endl <<
                                                                                // " fRecNegPart_primary = " << fRecNegPart_primary.size() << std::endl <<
                                                                                // " fRecPosPart_primary = " << fRecPosPart_primary.size() << std::endl <<
                                                                                // /*" fGenPairVec_primary = "   << fGenPairVec_primary.size()   << " fGenSmearedPairVec_primary =   " << fGenSmearedPairVec_primary.size()   << */
                                                                                // " fRecPairVec_primary =   " << fRecPairVec_primary.size()   << std::endl <<
                                                                                // /*" fGenPairVec_secondary = " << fGenPairVec_secondary.size() << " fGenSmearedPairVec_secondary = " << fGenSmearedPairVec_secondary.size() << */
                                                                                // " fRecPairVec_secondary = " << fRecV0Pair.size()/*fRecPairVec_secondary.size()*/ << std::endl;



                                                                                // if (fdebug) {
                                                                                // Cout to look at V0's daughter labels
                                                                                //   for (size_t i = 0; i < fRecV0Pair.size(); i++) {
                                                                                //     std::cout << "V0Status = " << fRecV0Pair.at(i).GetV0Status() << " Tracks: Label1 = " << fRecV0Pair.at(i).GetFirstDaughter() << " Label2 = " <<  fRecV0Pair.at(i).GetSecondDaughter() << std::endl;
                                                                                //   }
                                                                                // }



    // if(fdebug) std::cout << __LINE__ << " fRecV0Pair Vector size: " << fRecV0Pair.size() << std::endl;

  } // End of pairing

  /*  ------ \/ ------ Four Pairing ------ \/ ------  */
  if (fDoFourPairing){
    //##########################################################
    //################### Four Pairing #########################
    //##########################################################
                                                                                // if(fdebug) std::cout << "Doing four pairing..." << std::endl;
    Bool_t SmearedPair  = kTRUE;
    Bool_t ReconstructedPair = kTRUE;
    Bool_t PairPrimary = kTRUE;
    Bool_t TrackCuts = kTRUE;


    // #############################################
    // #     DoFourPairing of Gen & GenSmeared     #
    // #############################################
                                                                                // if (fdebug) std::cout << "Do generated four pairing" << std::endl;
    DoFourPairing(fGenPairVec_primary, fGenPairVec_secondary, !ReconstructedPair, !SmearedPair, centralityWeight);
                                                                                // if(fdebug) std::cout << __LINE__ << " DEBUG_AnalysisTask: Line cout " << std::endl;
    if(fArrResoPt){
      // if (fdebug) std::cout << "Do generated smeared four pairing" << std::endl;
      DoFourPairing(fGenSmearedPairVec_primary, fGenSmearedPairVec_secondary, !ReconstructedPair, SmearedPair, centralityWeight);
    }
                                                                                // if(fdebug) std::cout << __LINE__ << " DEBUG_AnalysisTask: Line cout " << std::endl;

    // #############################################
    // #     DoFourPairing Reconstructed           #
    // #############################################
    // #      Apply PreFitlers        #
    // ################################
                                                                                // if(fdebug) std::cout << __LINE__ << " Start Four PreFilter " << std::endl;
    // PreFilter using Prim and Sec Vector is rejecting tracks out of pos and neg particle vector, thats why the primary pair vector is cleared and primary paring is done again
    if(fUsePreFilter)DoFourPreFilter(&fRecPairVec_primary, &fRecV0Pair);
    if(fUseSecPreFilter)DoFourPreFilter(&fRecV0Pair, &fRecV0Pair);

                                                                                // if (fdebug) std::cout << __LINE__ << " DEBUG_AnalysisTask: Size of Vectors after applying PreFilters, " << std::endl <<
                                                                                // " fRecNegPart_primary = " << fRecNegPart_primary.size() << std::endl <<
                                                                                // " fRecPosPart_primary = " << fRecPosPart_primary.size() << std::endl <<
                                                                                // /*" fGenPairVec_primary = "   << fGenPairVec_primary.size()   << " fGenSmearedPairVec_primary =   " << fGenSmearedPairVec_primary.size()   << */
                                                                                // " fRecPairVec_primary =   " << fRecPairVec_primary.size()   << std::endl <<
                                                                                // /*" fGenPairVec_secondary = " << fGenPairVec_secondary.size() << " fGenSmearedPairVec_secondary = " << fGenSmearedPairVec_secondary.size() << */
                                                                                // " fRecPairVec_secondary = " << fRecV0Pair.size()/*fRecPairVec_secondary.size()*/ << std::endl;

    // Clear Primary Part Vector and do TwoPairing again with pos/neg particle vectors reduced by PreFilter
    fRecPairVec_primary.clear();
    DoRecTwoPairing(fRecNegPart_primary, fRecPosPart_primary, fPrimaryPairMCSignal,  PairPrimary, centralityWeight);

                                                                                // if (fdebug) std::cout << __LINE__ << " DEBUG_AnalysisTask: Size of Vectors after the second primary pairing, " << std::endl <<
                                                                                // " fRecNegPart_primary = " << fRecNegPart_primary.size() << std::endl <<
                                                                                // " fRecPosPart_primary = " << fRecPosPart_primary.size() << std::endl <<
                                                                                // /*" fGenPairVec_primary = "   << fGenPairVec_primary.size()   << " fGenSmearedPairVec_primary =   " << fGenSmearedPairVec_primary.size()   << */
                                                                                // " fRecPairVec_primary =   " << fRecPairVec_primary.size()   << std::endl <<
                                                                                // /*" fGenPairVec_secondary = " << fGenPairVec_secondary.size() << " fGenSmearedPairVec_secondary = " << fGenSmearedPairVec_secondary.size() << */
                                                                                // " fRecPairVec_secondary = " << fRecV0Pair.size()/*fRecPairVec_secondary.size()*/ << std::endl;

    // ################################
    // #      Apply Standard Cut      #
    // ################################
                                                                                // if(fdebug) std::cout << __LINE__ << " Apply StandardCuts for Pairing " << std::endl;
    ApplyStandardCutsAndFillHists(&fRecPairVec_primary, fTrackCuts_primary_standard  ,  TrackCuts,  PairPrimary, centralityWeight);
    ApplyStandardCutsAndFillHists(&fRecPairVec_primary, fPairCuts_primary            , !TrackCuts,  PairPrimary, centralityWeight);
    ApplyStandardCutsAndFillHists(&fRecV0Pair         , fTrackCuts_secondary_standard,  TrackCuts, !PairPrimary, centralityWeight);
    ApplyStandardCutsAndFillHists(&fRecV0Pair         , fPairCuts_secondary_standard , !TrackCuts, !PairPrimary, centralityWeight);


                                                                                // if (fdebug) std::cout << __LINE__ << " DEBUG_AnalysisTask: Size of Vectors after applying standard cuts, " << std::endl <<
                                                                                // " fRecNegPart_primary = " << fRecNegPart_primary.size() << std::endl <<
                                                                                // " fRecPosPart_primary = " << fRecPosPart_primary.size() << std::endl <<
                                                                                // /*" fGenPairVec_primary = "   << fGenPairVec_primary.size()   << " fGenSmearedPairVec_primary =   " << fGenSmearedPairVec_primary.size()   << */
                                                                                // " fRecPairVec_primary =   " << fRecPairVec_primary.size()   << std::endl <<
                                                                                // /*" fGenPairVec_secondary = " << fGenPairVec_secondary.size() << " fGenSmearedPairVec_secondary = " << fGenSmearedPairVec_secondary.size() << */
                                                                                // " fRecPairVec_secondary = " << fRecV0Pair.size()/*fRecPairVec_secondary.size()*/ << std::endl;



    // ################################
    // #      Do actual pairing       #
    // ################################
                                                                                // if(fdebug) std::cout << __LINE__ << " Start Four Reconstructed Pairing " << std::endl;
    DoFourPairing(fRecPairVec_primary, fRecV0Pair, ReconstructedPair, !SmearedPair, centralityWeight);
    // DoFourPairing(fRecPairVec_primary, fRecPairVec_secondary, ReconstructedPair, !SmearedPair, centralityWeight);

										                                                            // if(fdebug) std::cout << __LINE__ << " DEBUG_AnalysisTask: Line cout " << std::endl;
                                                                                // if(fdebug) std::cout << __LINE__ <<  "DEBUG_AnalysisTask: fGenNegPart_primary:   "  << fGenNegPart_primary.size()   << " fGenPosPart_primary:   " <<  fGenPosPart_primary.size()   << std::endl;
                                                                                // if(fdebug) std::cout << __LINE__ <<  "DEBUG_AnalysisTask: fGenNegPart_secondary: "  << fGenNegPart_secondary.size() << " fGenPosPart_secondary: " <<  fGenPosPart_secondary.size() << std::endl;

  } // end of fDoFourPairing
  /*  ------ /\ ------ Four Pairing ------ /\ ------  */


  PostData(1, fOutputList);
}


// ############################################################################
// ############################################################################
void    AliAnalysisTaskEtaReconstruction::FillTrackHistograms_Primary(AliVParticle* track, AliVParticle* mcTrack, int iCutList , int iMCSignal){
  Double_t values[AliDielectronVarManager::kNMaxValues]={0.};
  AliDielectronVarManager::SetFillMap(fUsedVars); // currently filled manually in the constructor of this task.

  AliDielectronVarManager::Fill(track, values);
  // std::cout << "pt var  manager = " << values[AliDielectronVarManager::kPt] << std::endl;
  // std::cout << "SITS    manager = " << values[AliDielectronVarManager::kNclsSITS] << std::endl;
  // std::cout << "TPCnSig manager = " << values[AliDielectronVarManager::kTPCnSigmaEle] << std::endl;
  TString genname;

    // (dynamic_cast<TH1D *>(fMCSigListVecPrim.at(iMCSignal).at(iCutList)->At(0)))->Fill(values[AliDielectronVarManager::kPt]);//hPt (reco)
    // (dynamic_cast<TH1D *>(fTrackCutListVecPrim.at(iCutList)->At(0)))->Fill(values[AliDielectronVarManager::kPt]);//hPt (reco)
    (dynamic_cast<TH1D *>(((TList*)((TList*)fTrackCutListVecPrim.at(iCutList))->At(iMCSignal))->At(0)))->Fill(values[AliDielectronVarManager::kPt]);//hPt (reco)
    (dynamic_cast<TH2D *>(((TList*)((TList*)fTrackCutListVecPrim.at(iCutList))->At(iMCSignal))->At(1)))->Fill(values[AliDielectronVarManager::kP],   values[AliDielectronVarManager::kITSnSigmaEle]);
    (dynamic_cast<TH2D *>(((TList*)((TList*)fTrackCutListVecPrim.at(iCutList))->At(iMCSignal))->At(2)))->Fill(values[AliDielectronVarManager::kPIn], values[AliDielectronVarManager::kTPCnSigmaEle]);
    (dynamic_cast<TH2D *>(((TList*)((TList*)fTrackCutListVecPrim.at(iCutList))->At(iMCSignal))->At(3)))->Fill(values[AliDielectronVarManager::kP],   values[AliDielectronVarManager::kTOFnSigmaEle]);
    (dynamic_cast<TH1D *>(((TList*)((TList*)fTrackCutListVecPrim.at(iCutList))->At(iMCSignal))->At(4)))->Fill(values[AliDielectronVarManager::kEta]);
    (dynamic_cast<TH1D *>(((TList*)((TList*)fTrackCutListVecPrim.at(iCutList))->At(iMCSignal))->At(5)))->Fill(values[AliDielectronVarManager::kPhi]);
    (dynamic_cast<TH2D *>(((TList*)((TList*)fTrackCutListVecPrim.at(iCutList))->At(iMCSignal))->At(6)))->Fill(values[AliDielectronVarManager::kEta], values[AliDielectronVarManager::kPhi]);

    (dynamic_cast<TH1D *>(((TList*)((TList*)fTrackCutListVecPrim.at(iCutList))->At(iMCSignal))->At(7 )))->Fill(values[AliDielectronVarManager::kNFclsTPCfCross]);
    (dynamic_cast<TH1D *>(((TList*)((TList*)fTrackCutListVecPrim.at(iCutList))->At(iMCSignal))->At(8 )))->Fill(values[AliDielectronVarManager::kNFclsTPCr]);
    (dynamic_cast<TH1D *>(((TList*)((TList*)fTrackCutListVecPrim.at(iCutList))->At(iMCSignal))->At(9 )))->Fill(values[AliDielectronVarManager::kNclsTPC]);
    (dynamic_cast<TH1D *>(((TList*)((TList*)fTrackCutListVecPrim.at(iCutList))->At(iMCSignal))->At(10)))->Fill(values[AliDielectronVarManager::kNclsITS]);
    (dynamic_cast<TH1D *>(((TList*)((TList*)fTrackCutListVecPrim.at(iCutList))->At(iMCSignal))->At(11)))->Fill(values[AliDielectronVarManager::kTPCchi2Cl]);
    (dynamic_cast<TH1D *>(((TList*)((TList*)fTrackCutListVecPrim.at(iCutList))->At(iMCSignal))->At(12)))->Fill(values[AliDielectronVarManager::kITSchi2Cl]);
    (dynamic_cast<TH1D *>(((TList*)((TList*)fTrackCutListVecPrim.at(iCutList))->At(iMCSignal))->At(13)))->Fill(values[AliDielectronVarManager::kNclsSTPC]);
    (dynamic_cast<TH1D *>(((TList*)((TList*)fTrackCutListVecPrim.at(iCutList))->At(iMCSignal))->At(14)))->Fill(values[AliDielectronVarManager::kNclsSITS]);
    (dynamic_cast<TH1D *>(((TList*)((TList*)fTrackCutListVecPrim.at(iCutList))->At(iMCSignal))->At(15)))->Fill(values[AliDielectronVarManager::kNclsSFracTPC]);
    (dynamic_cast<TH1D *>(((TList*)((TList*)fTrackCutListVecPrim.at(iCutList))->At(iMCSignal))->At(16)))->Fill(values[AliDielectronVarManager::kNclsSFracITS]);
    (dynamic_cast<TH1D *>(((TList*)((TList*)fTrackCutListVecPrim.at(iCutList))->At(iMCSignal))->At(17)))->Fill(values[AliDielectronVarManager::kTPCclsDiff]);
    (dynamic_cast<TH1D *>(((TList*)((TList*)fTrackCutListVecPrim.at(iCutList))->At(iMCSignal))->At(18)))->Fill(values[AliDielectronVarManager::kTPCsignalN]);
    (dynamic_cast<TH2D *>(((TList*)((TList*)fTrackCutListVecPrim.at(iCutList))->At(iMCSignal))->At(19)))->Fill(values[AliDielectronVarManager::kNclsTPC], values[AliDielectronVarManager::kNFclsTPCr]);
    (dynamic_cast<TH2D *>(((TList*)((TList*)fTrackCutListVecPrim.at(iCutList))->At(iMCSignal))->At(20)))->Fill(values[AliDielectronVarManager::kPt], values[AliDielectronVarManager::kNFclsTPCr]);
    // (dynamic_cast<TH1D *>(((TList*)((TList*)fTrackCutListVecPrim.at(iCutList))->At(iMCSignal))->At(21)))->Fill(values[AliDielectronVarManager::kPdgCode]);
    (dynamic_cast<TH1D *>(((TList*)((TList*)fTrackCutListVecPrim.at(iCutList))->At(iMCSignal))->At(21)))->Fill(mcTrack->PdgCode());
    (dynamic_cast<TH1D *>(((TList*)((TList*)fTrackCutListVecPrim.at(iCutList))->At(iMCSignal))->At(22)))->Fill( (fMC->GetTrack(TMath::Abs(mcTrack->GetMother())))->PdgCode());
    if(fMC->GetCocktailGenerator(TMath::Abs(track->GetLabel()), genname))    (dynamic_cast<TH1D *>(((TList*)((TList*)fTrackCutListVecPrim.at(iCutList))->At(iMCSignal))->At(23)))->Fill( genname,1);
    else (dynamic_cast<TH1D *>(((TList*)((TList*)fTrackCutListVecPrim.at(iCutList))->At(iMCSignal))->At(23)))->Fill( "none",1);
    (dynamic_cast<TH1D *>(((TList*)((TList*)fTrackCutListVecPrim.at(iCutList))->At(iMCSignal))->At(24)))->Fill(values[AliDielectronVarManager::kImpactParXY]);
    (dynamic_cast<TH1D *>(((TList*)((TList*)fTrackCutListVecPrim.at(iCutList))->At(iMCSignal))->At(25)))->Fill(values[AliDielectronVarManager::kImpactParZ]);
    (dynamic_cast<TH1D *>(((TList*)((TList*)fTrackCutListVecPrim.at(iCutList))->At(iMCSignal))->At(26)))->Fill(values[AliDielectronVarManager::kM]);
}

void    AliAnalysisTaskEtaReconstruction::FillTrackHistograms_Secondary(AliVParticle* track, AliVParticle* mcTrack, int iCutList, int iMCSignal){
  Double_t values[AliDielectronVarManager::kNMaxValues]={0.};
  AliDielectronVarManager::SetFillMap(fUsedVars); // currently filled manually in the constructor of this task.

  AliDielectronVarManager::Fill(track, values);
  // std::cout << "pt var  manager = " << values[AliDielectronVarManager::kPt] << std::endl;
  // std::cout << "SITS    manager = " << values[AliDielectronVarManager::kNclsSITS] << std::endl;
  // std::cout << "TPCnSig manager = " << values[AliDielectronVarManager::kTPCnSigmaEle] << std::endl;
  // std::cout << fTrackCutListVecSec.at(iCutList) << std::endl;
  TString genname;
    (dynamic_cast<TH1D *>(((TList*)((TList*)fTrackCutListVecSec.at(iCutList))->At(iMCSignal))->At(0)))->Fill(values[AliDielectronVarManager::kPt]);//hPt (reco)
    (dynamic_cast<TH2D *>(((TList*)((TList*)fTrackCutListVecSec.at(iCutList))->At(iMCSignal))->At(1)))->Fill(values[AliDielectronVarManager::kP],   values[AliDielectronVarManager::kITSnSigmaEle]);
    (dynamic_cast<TH2D *>(((TList*)((TList*)fTrackCutListVecSec.at(iCutList))->At(iMCSignal))->At(2)))->Fill(values[AliDielectronVarManager::kPIn], values[AliDielectronVarManager::kTPCnSigmaEle]);
    (dynamic_cast<TH2D *>(((TList*)((TList*)fTrackCutListVecSec.at(iCutList))->At(iMCSignal))->At(3)))->Fill(values[AliDielectronVarManager::kP],   values[AliDielectronVarManager::kTOFnSigmaEle]);
    (dynamic_cast<TH1D *>(((TList*)((TList*)fTrackCutListVecSec.at(iCutList))->At(iMCSignal))->At(4)))->Fill(values[AliDielectronVarManager::kEta]);
    (dynamic_cast<TH1D *>(((TList*)((TList*)fTrackCutListVecSec.at(iCutList))->At(iMCSignal))->At(5)))->Fill(values[AliDielectronVarManager::kPhi]);
    (dynamic_cast<TH2D *>(((TList*)((TList*)fTrackCutListVecSec.at(iCutList))->At(iMCSignal))->At(6)))->Fill(values[AliDielectronVarManager::kEta], values[AliDielectronVarManager::kPhi]);

    (dynamic_cast<TH1D *>(((TList*)((TList*)fTrackCutListVecSec.at(iCutList))->At(iMCSignal))->At(7 )))->Fill(values[AliDielectronVarManager::kNFclsTPCfCross]);
    (dynamic_cast<TH1D *>(((TList*)((TList*)fTrackCutListVecSec.at(iCutList))->At(iMCSignal))->At(8 )))->Fill(values[AliDielectronVarManager::kNFclsTPCr]);
    (dynamic_cast<TH1D *>(((TList*)((TList*)fTrackCutListVecSec.at(iCutList))->At(iMCSignal))->At(9 )))->Fill(values[AliDielectronVarManager::kNclsTPC]);
    (dynamic_cast<TH1D *>(((TList*)((TList*)fTrackCutListVecSec.at(iCutList))->At(iMCSignal))->At(10)))->Fill(values[AliDielectronVarManager::kNclsITS]);
    (dynamic_cast<TH1D *>(((TList*)((TList*)fTrackCutListVecSec.at(iCutList))->At(iMCSignal))->At(11)))->Fill(values[AliDielectronVarManager::kTPCchi2Cl]);
    (dynamic_cast<TH1D *>(((TList*)((TList*)fTrackCutListVecSec.at(iCutList))->At(iMCSignal))->At(12)))->Fill(values[AliDielectronVarManager::kITSchi2Cl]);
    (dynamic_cast<TH1D *>(((TList*)((TList*)fTrackCutListVecSec.at(iCutList))->At(iMCSignal))->At(13)))->Fill(values[AliDielectronVarManager::kNclsSTPC]);
    (dynamic_cast<TH1D *>(((TList*)((TList*)fTrackCutListVecSec.at(iCutList))->At(iMCSignal))->At(14)))->Fill(values[AliDielectronVarManager::kNclsSITS]);
    (dynamic_cast<TH1D *>(((TList*)((TList*)fTrackCutListVecSec.at(iCutList))->At(iMCSignal))->At(15)))->Fill(values[AliDielectronVarManager::kNclsSFracTPC]);
    (dynamic_cast<TH1D *>(((TList*)((TList*)fTrackCutListVecSec.at(iCutList))->At(iMCSignal))->At(16)))->Fill(values[AliDielectronVarManager::kNclsSFracITS]);
    (dynamic_cast<TH1D *>(((TList*)((TList*)fTrackCutListVecSec.at(iCutList))->At(iMCSignal))->At(17)))->Fill(values[AliDielectronVarManager::kTPCclsDiff]);
    (dynamic_cast<TH1D *>(((TList*)((TList*)fTrackCutListVecSec.at(iCutList))->At(iMCSignal))->At(18)))->Fill(values[AliDielectronVarManager::kTPCsignalN]);
    (dynamic_cast<TH2D *>(((TList*)((TList*)fTrackCutListVecSec.at(iCutList))->At(iMCSignal))->At(19)))->Fill(values[AliDielectronVarManager::kNclsTPC], values[AliDielectronVarManager::kNFclsTPCr]);
    (dynamic_cast<TH2D *>(((TList*)((TList*)fTrackCutListVecSec.at(iCutList))->At(iMCSignal))->At(20)))->Fill(values[AliDielectronVarManager::kPt], values[AliDielectronVarManager::kNFclsTPCr]);
    // (dynamic_cast<TH1D *>(((TList*)((TList*)fTrackCutListVecSec.at(iCutList))->At(iMCSignal))->At(21)))->Fill(values[AliDielectronVarManager::kPdgCode]);
    (dynamic_cast<TH1D *>(((TList*)((TList*)fTrackCutListVecSec.at(iCutList))->At(iMCSignal))->At(21)))->Fill(mcTrack->PdgCode());
    (dynamic_cast<TH1D *>(((TList*)((TList*)fTrackCutListVecSec.at(iCutList))->At(iMCSignal))->At(22)))->Fill( (fMC->GetTrack(TMath::Abs(mcTrack->GetMother())))->PdgCode());
    if(fMC->GetCocktailGenerator(TMath::Abs(track->GetLabel()), genname))    (dynamic_cast<TH1D *>(((TList*)((TList*)fTrackCutListVecSec.at(iCutList))->At(iMCSignal))->At(23)))->Fill( genname,1);
    else (dynamic_cast<TH1D *>(((TList*)((TList*)fTrackCutListVecSec.at(iCutList))->At(iMCSignal))->At(23)))->Fill( "none",1);
    (dynamic_cast<TH1D *>(((TList*)((TList*)fTrackCutListVecSec.at(iCutList))->At(iMCSignal))->At(24)))->Fill(values[AliDielectronVarManager::kImpactParXY]);
    (dynamic_cast<TH1D *>(((TList*)((TList*)fTrackCutListVecSec.at(iCutList))->At(iMCSignal))->At(25)))->Fill(values[AliDielectronVarManager::kImpactParZ]);
    (dynamic_cast<TH1D *>(((TList*)((TList*)fTrackCutListVecSec.at(iCutList))->At(iMCSignal))->At(26)))->Fill(values[AliDielectronVarManager::kM]);

                                                (dynamic_cast<TH2D *>(((TList*)((TList*)((TList*)fTrackCutListVecSec.at(iCutList))->At(iMCSignal))->At(27))->At(0)))->Fill(values[AliDielectronVarManager::kP], values[AliDielectronVarManager::kITSnSigmaEle]);
    if(TMath::Abs(mcTrack->PdgCode()) == 11 )   (dynamic_cast<TH2D *>(((TList*)((TList*)((TList*)fTrackCutListVecSec.at(iCutList))->At(iMCSignal))->At(27))->At(1)))->Fill(values[AliDielectronVarManager::kP], values[AliDielectronVarManager::kITSnSigmaEle]);
    if(TMath::Abs(mcTrack->PdgCode()) == 13 )   (dynamic_cast<TH2D *>(((TList*)((TList*)((TList*)fTrackCutListVecSec.at(iCutList))->At(iMCSignal))->At(27))->At(2)))->Fill(values[AliDielectronVarManager::kP], values[AliDielectronVarManager::kITSnSigmaEle]);
    if(TMath::Abs(mcTrack->PdgCode()) == 211)   (dynamic_cast<TH2D *>(((TList*)((TList*)((TList*)fTrackCutListVecSec.at(iCutList))->At(iMCSignal))->At(27))->At(3)))->Fill(values[AliDielectronVarManager::kP], values[AliDielectronVarManager::kITSnSigmaEle]);
    if(TMath::Abs(mcTrack->PdgCode()) == 321)   (dynamic_cast<TH2D *>(((TList*)((TList*)((TList*)fTrackCutListVecSec.at(iCutList))->At(iMCSignal))->At(27))->At(4)))->Fill(values[AliDielectronVarManager::kP], values[AliDielectronVarManager::kITSnSigmaEle]);
    if(TMath::Abs(mcTrack->PdgCode()) == 2212)  (dynamic_cast<TH2D *>(((TList*)((TList*)((TList*)fTrackCutListVecSec.at(iCutList))->At(iMCSignal))->At(27))->At(5)))->Fill(values[AliDielectronVarManager::kP], values[AliDielectronVarManager::kITSnSigmaEle]);

                                                (dynamic_cast<TH2D *>(((TList*)((TList*)((TList*)fTrackCutListVecSec.at(iCutList))->At(iMCSignal))->At(27))->At(6)))->Fill(values[AliDielectronVarManager::kP], values[AliDielectronVarManager::kITSnSigmaMuo]);
    if(TMath::Abs(mcTrack->PdgCode()) == 11 )   (dynamic_cast<TH2D *>(((TList*)((TList*)((TList*)fTrackCutListVecSec.at(iCutList))->At(iMCSignal))->At(27))->At(7)))->Fill(values[AliDielectronVarManager::kP], values[AliDielectronVarManager::kITSnSigmaMuo]);
    if(TMath::Abs(mcTrack->PdgCode()) == 13 )   (dynamic_cast<TH2D *>(((TList*)((TList*)((TList*)fTrackCutListVecSec.at(iCutList))->At(iMCSignal))->At(27))->At(8)))->Fill(values[AliDielectronVarManager::kP], values[AliDielectronVarManager::kITSnSigmaMuo]);
    if(TMath::Abs(mcTrack->PdgCode()) == 211)   (dynamic_cast<TH2D *>(((TList*)((TList*)((TList*)fTrackCutListVecSec.at(iCutList))->At(iMCSignal))->At(27))->At(9)))->Fill(values[AliDielectronVarManager::kP], values[AliDielectronVarManager::kITSnSigmaMuo]);
    if(TMath::Abs(mcTrack->PdgCode()) == 321)   (dynamic_cast<TH2D *>(((TList*)((TList*)((TList*)fTrackCutListVecSec.at(iCutList))->At(iMCSignal))->At(27))->At(10)))->Fill(values[AliDielectronVarManager::kP], values[AliDielectronVarManager::kITSnSigmaMuo]);
    if(TMath::Abs(mcTrack->PdgCode()) == 2212)  (dynamic_cast<TH2D *>(((TList*)((TList*)((TList*)fTrackCutListVecSec.at(iCutList))->At(iMCSignal))->At(27))->At(11)))->Fill(values[AliDielectronVarManager::kP], values[AliDielectronVarManager::kITSnSigmaMuo]);

                                                (dynamic_cast<TH2D *>(((TList*)((TList*)((TList*)fTrackCutListVecSec.at(iCutList))->At(iMCSignal))->At(27))->At(12)))->Fill(values[AliDielectronVarManager::kP], values[AliDielectronVarManager::kITSnSigmaPio]);
    if(TMath::Abs(mcTrack->PdgCode()) == 11 )   (dynamic_cast<TH2D *>(((TList*)((TList*)((TList*)fTrackCutListVecSec.at(iCutList))->At(iMCSignal))->At(27))->At(13)))->Fill(values[AliDielectronVarManager::kP], values[AliDielectronVarManager::kITSnSigmaPio]);
    if(TMath::Abs(mcTrack->PdgCode()) == 13 )   (dynamic_cast<TH2D *>(((TList*)((TList*)((TList*)fTrackCutListVecSec.at(iCutList))->At(iMCSignal))->At(27))->At(14)))->Fill(values[AliDielectronVarManager::kP], values[AliDielectronVarManager::kITSnSigmaPio]);
    if(TMath::Abs(mcTrack->PdgCode()) == 211)   (dynamic_cast<TH2D *>(((TList*)((TList*)((TList*)fTrackCutListVecSec.at(iCutList))->At(iMCSignal))->At(27))->At(15)))->Fill(values[AliDielectronVarManager::kP], values[AliDielectronVarManager::kITSnSigmaPio]);
    if(TMath::Abs(mcTrack->PdgCode()) == 321)   (dynamic_cast<TH2D *>(((TList*)((TList*)((TList*)fTrackCutListVecSec.at(iCutList))->At(iMCSignal))->At(27))->At(16)))->Fill(values[AliDielectronVarManager::kP], values[AliDielectronVarManager::kITSnSigmaPio]);
    if(TMath::Abs(mcTrack->PdgCode()) == 2212)  (dynamic_cast<TH2D *>(((TList*)((TList*)((TList*)fTrackCutListVecSec.at(iCutList))->At(iMCSignal))->At(27))->At(17)))->Fill(values[AliDielectronVarManager::kP], values[AliDielectronVarManager::kITSnSigmaPio]);

                                                (dynamic_cast<TH2D *>(((TList*)((TList*)((TList*)fTrackCutListVecSec.at(iCutList))->At(iMCSignal))->At(27))->At(18)))->Fill(values[AliDielectronVarManager::kP], values[AliDielectronVarManager::kITSnSigmaKao]);
    if(TMath::Abs(mcTrack->PdgCode()) == 11 )   (dynamic_cast<TH2D *>(((TList*)((TList*)((TList*)fTrackCutListVecSec.at(iCutList))->At(iMCSignal))->At(27))->At(19)))->Fill(values[AliDielectronVarManager::kP], values[AliDielectronVarManager::kITSnSigmaKao]);
    if(TMath::Abs(mcTrack->PdgCode()) == 13 )   (dynamic_cast<TH2D *>(((TList*)((TList*)((TList*)fTrackCutListVecSec.at(iCutList))->At(iMCSignal))->At(27))->At(20)))->Fill(values[AliDielectronVarManager::kP], values[AliDielectronVarManager::kITSnSigmaKao]);
    if(TMath::Abs(mcTrack->PdgCode()) == 211)   (dynamic_cast<TH2D *>(((TList*)((TList*)((TList*)fTrackCutListVecSec.at(iCutList))->At(iMCSignal))->At(27))->At(21)))->Fill(values[AliDielectronVarManager::kP], values[AliDielectronVarManager::kITSnSigmaKao]);
    if(TMath::Abs(mcTrack->PdgCode()) == 321)   (dynamic_cast<TH2D *>(((TList*)((TList*)((TList*)fTrackCutListVecSec.at(iCutList))->At(iMCSignal))->At(27))->At(22)))->Fill(values[AliDielectronVarManager::kP], values[AliDielectronVarManager::kITSnSigmaKao]);
    if(TMath::Abs(mcTrack->PdgCode()) == 2212)  (dynamic_cast<TH2D *>(((TList*)((TList*)((TList*)fTrackCutListVecSec.at(iCutList))->At(iMCSignal))->At(27))->At(23)))->Fill(values[AliDielectronVarManager::kP], values[AliDielectronVarManager::kITSnSigmaKao]);

                                                (dynamic_cast<TH2D *>(((TList*)((TList*)((TList*)fTrackCutListVecSec.at(iCutList))->At(iMCSignal))->At(27))->At(24)))->Fill(values[AliDielectronVarManager::kP], values[AliDielectronVarManager::kITSnSigmaPro]);
    if(TMath::Abs(mcTrack->PdgCode()) == 11 )   (dynamic_cast<TH2D *>(((TList*)((TList*)((TList*)fTrackCutListVecSec.at(iCutList))->At(iMCSignal))->At(27))->At(25)))->Fill(values[AliDielectronVarManager::kP], values[AliDielectronVarManager::kITSnSigmaPro]);
    if(TMath::Abs(mcTrack->PdgCode()) == 13 )   (dynamic_cast<TH2D *>(((TList*)((TList*)((TList*)fTrackCutListVecSec.at(iCutList))->At(iMCSignal))->At(27))->At(26)))->Fill(values[AliDielectronVarManager::kP], values[AliDielectronVarManager::kITSnSigmaPro]);
    if(TMath::Abs(mcTrack->PdgCode()) == 211)   (dynamic_cast<TH2D *>(((TList*)((TList*)((TList*)fTrackCutListVecSec.at(iCutList))->At(iMCSignal))->At(27))->At(27)))->Fill(values[AliDielectronVarManager::kP], values[AliDielectronVarManager::kITSnSigmaPro]);
    if(TMath::Abs(mcTrack->PdgCode()) == 321)   (dynamic_cast<TH2D *>(((TList*)((TList*)((TList*)fTrackCutListVecSec.at(iCutList))->At(iMCSignal))->At(27))->At(28)))->Fill(values[AliDielectronVarManager::kP], values[AliDielectronVarManager::kITSnSigmaPro]);
    if(TMath::Abs(mcTrack->PdgCode()) == 2212)  (dynamic_cast<TH2D *>(((TList*)((TList*)((TList*)fTrackCutListVecSec.at(iCutList))->At(iMCSignal))->At(27))->At(29)))->Fill(values[AliDielectronVarManager::kP], values[AliDielectronVarManager::kITSnSigmaPro]);

                                                (dynamic_cast<TH2D *>(((TList*)((TList*)((TList*)fTrackCutListVecSec.at(iCutList))->At(iMCSignal))->At(28))->At(0)))->Fill(values[AliDielectronVarManager::kPIn], values[AliDielectronVarManager::kTPCnSigmaEle]);
    if(TMath::Abs(mcTrack->PdgCode()) == 11 )   (dynamic_cast<TH2D *>(((TList*)((TList*)((TList*)fTrackCutListVecSec.at(iCutList))->At(iMCSignal))->At(28))->At(1)))->Fill(values[AliDielectronVarManager::kPIn], values[AliDielectronVarManager::kTPCnSigmaEle]);
    if(TMath::Abs(mcTrack->PdgCode()) == 13 )   (dynamic_cast<TH2D *>(((TList*)((TList*)((TList*)fTrackCutListVecSec.at(iCutList))->At(iMCSignal))->At(28))->At(2)))->Fill(values[AliDielectronVarManager::kPIn], values[AliDielectronVarManager::kTPCnSigmaEle]);
    if(TMath::Abs(mcTrack->PdgCode()) == 211)   (dynamic_cast<TH2D *>(((TList*)((TList*)((TList*)fTrackCutListVecSec.at(iCutList))->At(iMCSignal))->At(28))->At(3)))->Fill(values[AliDielectronVarManager::kPIn], values[AliDielectronVarManager::kTPCnSigmaEle]);
    if(TMath::Abs(mcTrack->PdgCode()) == 321)   (dynamic_cast<TH2D *>(((TList*)((TList*)((TList*)fTrackCutListVecSec.at(iCutList))->At(iMCSignal))->At(28))->At(4)))->Fill(values[AliDielectronVarManager::kPIn], values[AliDielectronVarManager::kTPCnSigmaEle]);
    if(TMath::Abs(mcTrack->PdgCode()) == 2212)  (dynamic_cast<TH2D *>(((TList*)((TList*)((TList*)fTrackCutListVecSec.at(iCutList))->At(iMCSignal))->At(28))->At(5)))->Fill(values[AliDielectronVarManager::kPIn], values[AliDielectronVarManager::kTPCnSigmaEle]);

                                                (dynamic_cast<TH2D *>(((TList*)((TList*)((TList*)fTrackCutListVecSec.at(iCutList))->At(iMCSignal))->At(28))->At(6)))->Fill(values[AliDielectronVarManager::kPIn], values[AliDielectronVarManager::kTPCnSigmaMuo]);
    if(TMath::Abs(mcTrack->PdgCode()) == 11 )   (dynamic_cast<TH2D *>(((TList*)((TList*)((TList*)fTrackCutListVecSec.at(iCutList))->At(iMCSignal))->At(28))->At(7)))->Fill(values[AliDielectronVarManager::kPIn], values[AliDielectronVarManager::kTPCnSigmaMuo]);
    if(TMath::Abs(mcTrack->PdgCode()) == 13 )   (dynamic_cast<TH2D *>(((TList*)((TList*)((TList*)fTrackCutListVecSec.at(iCutList))->At(iMCSignal))->At(28))->At(8)))->Fill(values[AliDielectronVarManager::kPIn], values[AliDielectronVarManager::kTPCnSigmaMuo]);
    if(TMath::Abs(mcTrack->PdgCode()) == 211)   (dynamic_cast<TH2D *>(((TList*)((TList*)((TList*)fTrackCutListVecSec.at(iCutList))->At(iMCSignal))->At(28))->At(9)))->Fill(values[AliDielectronVarManager::kPIn], values[AliDielectronVarManager::kTPCnSigmaMuo]);
    if(TMath::Abs(mcTrack->PdgCode()) == 321)   (dynamic_cast<TH2D *>(((TList*)((TList*)((TList*)fTrackCutListVecSec.at(iCutList))->At(iMCSignal))->At(28))->At(10)))->Fill(values[AliDielectronVarManager::kPIn], values[AliDielectronVarManager::kTPCnSigmaMuo]);
    if(TMath::Abs(mcTrack->PdgCode()) == 2212)  (dynamic_cast<TH2D *>(((TList*)((TList*)((TList*)fTrackCutListVecSec.at(iCutList))->At(iMCSignal))->At(28))->At(11)))->Fill(values[AliDielectronVarManager::kPIn], values[AliDielectronVarManager::kTPCnSigmaMuo]);

                                                (dynamic_cast<TH2D *>(((TList*)((TList*)((TList*)fTrackCutListVecSec.at(iCutList))->At(iMCSignal))->At(28))->At(12)))->Fill(values[AliDielectronVarManager::kPIn], values[AliDielectronVarManager::kTPCnSigmaPio]);
    if(TMath::Abs(mcTrack->PdgCode()) == 11 )   (dynamic_cast<TH2D *>(((TList*)((TList*)((TList*)fTrackCutListVecSec.at(iCutList))->At(iMCSignal))->At(28))->At(13)))->Fill(values[AliDielectronVarManager::kPIn], values[AliDielectronVarManager::kTPCnSigmaPio]);
    if(TMath::Abs(mcTrack->PdgCode()) == 13 )   (dynamic_cast<TH2D *>(((TList*)((TList*)((TList*)fTrackCutListVecSec.at(iCutList))->At(iMCSignal))->At(28))->At(14)))->Fill(values[AliDielectronVarManager::kPIn], values[AliDielectronVarManager::kTPCnSigmaPio]);
    if(TMath::Abs(mcTrack->PdgCode()) == 211)   (dynamic_cast<TH2D *>(((TList*)((TList*)((TList*)fTrackCutListVecSec.at(iCutList))->At(iMCSignal))->At(28))->At(15)))->Fill(values[AliDielectronVarManager::kPIn], values[AliDielectronVarManager::kTPCnSigmaPio]);
    if(TMath::Abs(mcTrack->PdgCode()) == 321)   (dynamic_cast<TH2D *>(((TList*)((TList*)((TList*)fTrackCutListVecSec.at(iCutList))->At(iMCSignal))->At(28))->At(16)))->Fill(values[AliDielectronVarManager::kPIn], values[AliDielectronVarManager::kTPCnSigmaPio]);
    if(TMath::Abs(mcTrack->PdgCode()) == 2212)  (dynamic_cast<TH2D *>(((TList*)((TList*)((TList*)fTrackCutListVecSec.at(iCutList))->At(iMCSignal))->At(28))->At(17)))->Fill(values[AliDielectronVarManager::kPIn], values[AliDielectronVarManager::kTPCnSigmaPio]);

                                                (dynamic_cast<TH2D *>(((TList*)((TList*)((TList*)fTrackCutListVecSec.at(iCutList))->At(iMCSignal))->At(28))->At(18)))->Fill(values[AliDielectronVarManager::kPIn], values[AliDielectronVarManager::kTPCnSigmaKao]);
    if(TMath::Abs(mcTrack->PdgCode()) == 11 )   (dynamic_cast<TH2D *>(((TList*)((TList*)((TList*)fTrackCutListVecSec.at(iCutList))->At(iMCSignal))->At(28))->At(19)))->Fill(values[AliDielectronVarManager::kPIn], values[AliDielectronVarManager::kTPCnSigmaKao]);
    if(TMath::Abs(mcTrack->PdgCode()) == 13 )   (dynamic_cast<TH2D *>(((TList*)((TList*)((TList*)fTrackCutListVecSec.at(iCutList))->At(iMCSignal))->At(28))->At(20)))->Fill(values[AliDielectronVarManager::kPIn], values[AliDielectronVarManager::kTPCnSigmaKao]);
    if(TMath::Abs(mcTrack->PdgCode()) == 211)   (dynamic_cast<TH2D *>(((TList*)((TList*)((TList*)fTrackCutListVecSec.at(iCutList))->At(iMCSignal))->At(28))->At(21)))->Fill(values[AliDielectronVarManager::kPIn], values[AliDielectronVarManager::kTPCnSigmaKao]);
    if(TMath::Abs(mcTrack->PdgCode()) == 321)   (dynamic_cast<TH2D *>(((TList*)((TList*)((TList*)fTrackCutListVecSec.at(iCutList))->At(iMCSignal))->At(28))->At(22)))->Fill(values[AliDielectronVarManager::kPIn], values[AliDielectronVarManager::kTPCnSigmaKao]);
    if(TMath::Abs(mcTrack->PdgCode()) == 2212)  (dynamic_cast<TH2D *>(((TList*)((TList*)((TList*)fTrackCutListVecSec.at(iCutList))->At(iMCSignal))->At(28))->At(23)))->Fill(values[AliDielectronVarManager::kPIn], values[AliDielectronVarManager::kTPCnSigmaKao]);

                                                (dynamic_cast<TH2D *>(((TList*)((TList*)((TList*)fTrackCutListVecSec.at(iCutList))->At(iMCSignal))->At(28))->At(24)))->Fill(values[AliDielectronVarManager::kPIn], values[AliDielectronVarManager::kTPCnSigmaPro]);
    if(TMath::Abs(mcTrack->PdgCode()) == 11 )   (dynamic_cast<TH2D *>(((TList*)((TList*)((TList*)fTrackCutListVecSec.at(iCutList))->At(iMCSignal))->At(28))->At(25)))->Fill(values[AliDielectronVarManager::kPIn], values[AliDielectronVarManager::kTPCnSigmaPro]);
    if(TMath::Abs(mcTrack->PdgCode()) == 13 )   (dynamic_cast<TH2D *>(((TList*)((TList*)((TList*)fTrackCutListVecSec.at(iCutList))->At(iMCSignal))->At(28))->At(26)))->Fill(values[AliDielectronVarManager::kPIn], values[AliDielectronVarManager::kTPCnSigmaPro]);
    if(TMath::Abs(mcTrack->PdgCode()) == 211)   (dynamic_cast<TH2D *>(((TList*)((TList*)((TList*)fTrackCutListVecSec.at(iCutList))->At(iMCSignal))->At(28))->At(27)))->Fill(values[AliDielectronVarManager::kPIn], values[AliDielectronVarManager::kTPCnSigmaPro]);
    if(TMath::Abs(mcTrack->PdgCode()) == 321)   (dynamic_cast<TH2D *>(((TList*)((TList*)((TList*)fTrackCutListVecSec.at(iCutList))->At(iMCSignal))->At(28))->At(28)))->Fill(values[AliDielectronVarManager::kPIn], values[AliDielectronVarManager::kTPCnSigmaPro]);
    if(TMath::Abs(mcTrack->PdgCode()) == 2212)  (dynamic_cast<TH2D *>(((TList*)((TList*)((TList*)fTrackCutListVecSec.at(iCutList))->At(iMCSignal))->At(28))->At(29)))->Fill(values[AliDielectronVarManager::kPIn], values[AliDielectronVarManager::kTPCnSigmaPro]);

                                                (dynamic_cast<TH2D *>(((TList*)((TList*)((TList*)fTrackCutListVecSec.at(iCutList))->At(iMCSignal))->At(29))->At(0)))->Fill(values[AliDielectronVarManager::kP], values[AliDielectronVarManager::kTPCnSigmaEle]);
    if(TMath::Abs(mcTrack->PdgCode()) == 11 )   (dynamic_cast<TH2D *>(((TList*)((TList*)((TList*)fTrackCutListVecSec.at(iCutList))->At(iMCSignal))->At(29))->At(1)))->Fill(values[AliDielectronVarManager::kP], values[AliDielectronVarManager::kTOFnSigmaEle]);
    if(TMath::Abs(mcTrack->PdgCode()) == 13 )   (dynamic_cast<TH2D *>(((TList*)((TList*)((TList*)fTrackCutListVecSec.at(iCutList))->At(iMCSignal))->At(29))->At(2)))->Fill(values[AliDielectronVarManager::kP], values[AliDielectronVarManager::kTOFnSigmaEle]);
    if(TMath::Abs(mcTrack->PdgCode()) == 211)   (dynamic_cast<TH2D *>(((TList*)((TList*)((TList*)fTrackCutListVecSec.at(iCutList))->At(iMCSignal))->At(29))->At(3)))->Fill(values[AliDielectronVarManager::kP], values[AliDielectronVarManager::kTOFnSigmaEle]);
    if(TMath::Abs(mcTrack->PdgCode()) == 321)   (dynamic_cast<TH2D *>(((TList*)((TList*)((TList*)fTrackCutListVecSec.at(iCutList))->At(iMCSignal))->At(29))->At(4)))->Fill(values[AliDielectronVarManager::kP], values[AliDielectronVarManager::kTOFnSigmaEle]);
    if(TMath::Abs(mcTrack->PdgCode()) == 2212)  (dynamic_cast<TH2D *>(((TList*)((TList*)((TList*)fTrackCutListVecSec.at(iCutList))->At(iMCSignal))->At(29))->At(5)))->Fill(values[AliDielectronVarManager::kP], values[AliDielectronVarManager::kTOFnSigmaEle]);

                                                (dynamic_cast<TH2D *>(((TList*)((TList*)((TList*)fTrackCutListVecSec.at(iCutList))->At(iMCSignal))->At(29))->At(6)))->Fill(values[AliDielectronVarManager::kP], values[AliDielectronVarManager::kTPCnSigmaMuo]);
    if(TMath::Abs(mcTrack->PdgCode()) == 11 )   (dynamic_cast<TH2D *>(((TList*)((TList*)((TList*)fTrackCutListVecSec.at(iCutList))->At(iMCSignal))->At(29))->At(7)))->Fill(values[AliDielectronVarManager::kP], values[AliDielectronVarManager::kTOFnSigmaMuo]);
    if(TMath::Abs(mcTrack->PdgCode()) == 13 )   (dynamic_cast<TH2D *>(((TList*)((TList*)((TList*)fTrackCutListVecSec.at(iCutList))->At(iMCSignal))->At(29))->At(8)))->Fill(values[AliDielectronVarManager::kP], values[AliDielectronVarManager::kTOFnSigmaMuo]);
    if(TMath::Abs(mcTrack->PdgCode()) == 211)   (dynamic_cast<TH2D *>(((TList*)((TList*)((TList*)fTrackCutListVecSec.at(iCutList))->At(iMCSignal))->At(29))->At(9)))->Fill(values[AliDielectronVarManager::kP], values[AliDielectronVarManager::kTOFnSigmaMuo]);
    if(TMath::Abs(mcTrack->PdgCode()) == 321)   (dynamic_cast<TH2D *>(((TList*)((TList*)((TList*)fTrackCutListVecSec.at(iCutList))->At(iMCSignal))->At(29))->At(10)))->Fill(values[AliDielectronVarManager::kP], values[AliDielectronVarManager::kTOFnSigmaMuo]);
    if(TMath::Abs(mcTrack->PdgCode()) == 2212)  (dynamic_cast<TH2D *>(((TList*)((TList*)((TList*)fTrackCutListVecSec.at(iCutList))->At(iMCSignal))->At(29))->At(11)))->Fill(values[AliDielectronVarManager::kP], values[AliDielectronVarManager::kTOFnSigmaMuo]);

                                                (dynamic_cast<TH2D *>(((TList*)((TList*)((TList*)fTrackCutListVecSec.at(iCutList))->At(iMCSignal))->At(29))->At(12)))->Fill(values[AliDielectronVarManager::kP], values[AliDielectronVarManager::kTPCnSigmaPio]);
    if(TMath::Abs(mcTrack->PdgCode()) == 11 )   (dynamic_cast<TH2D *>(((TList*)((TList*)((TList*)fTrackCutListVecSec.at(iCutList))->At(iMCSignal))->At(29))->At(13)))->Fill(values[AliDielectronVarManager::kP], values[AliDielectronVarManager::kTOFnSigmaPio]);
    if(TMath::Abs(mcTrack->PdgCode()) == 13 )   (dynamic_cast<TH2D *>(((TList*)((TList*)((TList*)fTrackCutListVecSec.at(iCutList))->At(iMCSignal))->At(29))->At(14)))->Fill(values[AliDielectronVarManager::kP], values[AliDielectronVarManager::kTOFnSigmaPio]);
    if(TMath::Abs(mcTrack->PdgCode()) == 211)   (dynamic_cast<TH2D *>(((TList*)((TList*)((TList*)fTrackCutListVecSec.at(iCutList))->At(iMCSignal))->At(29))->At(15)))->Fill(values[AliDielectronVarManager::kP], values[AliDielectronVarManager::kTOFnSigmaPio]);
    if(TMath::Abs(mcTrack->PdgCode()) == 321)   (dynamic_cast<TH2D *>(((TList*)((TList*)((TList*)fTrackCutListVecSec.at(iCutList))->At(iMCSignal))->At(29))->At(16)))->Fill(values[AliDielectronVarManager::kP], values[AliDielectronVarManager::kTOFnSigmaPio]);
    if(TMath::Abs(mcTrack->PdgCode()) == 2212)  (dynamic_cast<TH2D *>(((TList*)((TList*)((TList*)fTrackCutListVecSec.at(iCutList))->At(iMCSignal))->At(29))->At(17)))->Fill(values[AliDielectronVarManager::kP], values[AliDielectronVarManager::kTOFnSigmaPio]);

                                                (dynamic_cast<TH2D *>(((TList*)((TList*)((TList*)fTrackCutListVecSec.at(iCutList))->At(iMCSignal))->At(29))->At(18)))->Fill(values[AliDielectronVarManager::kP], values[AliDielectronVarManager::kTPCnSigmaKao]);
    if(TMath::Abs(mcTrack->PdgCode()) == 11 )   (dynamic_cast<TH2D *>(((TList*)((TList*)((TList*)fTrackCutListVecSec.at(iCutList))->At(iMCSignal))->At(29))->At(19)))->Fill(values[AliDielectronVarManager::kP], values[AliDielectronVarManager::kTOFnSigmaKao]);
    if(TMath::Abs(mcTrack->PdgCode()) == 13 )   (dynamic_cast<TH2D *>(((TList*)((TList*)((TList*)fTrackCutListVecSec.at(iCutList))->At(iMCSignal))->At(29))->At(20)))->Fill(values[AliDielectronVarManager::kP], values[AliDielectronVarManager::kTOFnSigmaKao]);
    if(TMath::Abs(mcTrack->PdgCode()) == 211)   (dynamic_cast<TH2D *>(((TList*)((TList*)((TList*)fTrackCutListVecSec.at(iCutList))->At(iMCSignal))->At(29))->At(21)))->Fill(values[AliDielectronVarManager::kP], values[AliDielectronVarManager::kTOFnSigmaKao]);
    if(TMath::Abs(mcTrack->PdgCode()) == 321)   (dynamic_cast<TH2D *>(((TList*)((TList*)((TList*)fTrackCutListVecSec.at(iCutList))->At(iMCSignal))->At(29))->At(22)))->Fill(values[AliDielectronVarManager::kP], values[AliDielectronVarManager::kTOFnSigmaKao]);
    if(TMath::Abs(mcTrack->PdgCode()) == 2212)  (dynamic_cast<TH2D *>(((TList*)((TList*)((TList*)fTrackCutListVecSec.at(iCutList))->At(iMCSignal))->At(29))->At(23)))->Fill(values[AliDielectronVarManager::kP], values[AliDielectronVarManager::kTOFnSigmaKao]);

                                                (dynamic_cast<TH2D *>(((TList*)((TList*)((TList*)fTrackCutListVecSec.at(iCutList))->At(iMCSignal))->At(29))->At(24)))->Fill(values[AliDielectronVarManager::kP], values[AliDielectronVarManager::kTPCnSigmaPro]);
    if(TMath::Abs(mcTrack->PdgCode()) == 11 )   (dynamic_cast<TH2D *>(((TList*)((TList*)((TList*)fTrackCutListVecSec.at(iCutList))->At(iMCSignal))->At(29))->At(25)))->Fill(values[AliDielectronVarManager::kP], values[AliDielectronVarManager::kTOFnSigmaPro]);
    if(TMath::Abs(mcTrack->PdgCode()) == 13 )   (dynamic_cast<TH2D *>(((TList*)((TList*)((TList*)fTrackCutListVecSec.at(iCutList))->At(iMCSignal))->At(29))->At(26)))->Fill(values[AliDielectronVarManager::kP], values[AliDielectronVarManager::kTOFnSigmaPro]);
    if(TMath::Abs(mcTrack->PdgCode()) == 211)   (dynamic_cast<TH2D *>(((TList*)((TList*)((TList*)fTrackCutListVecSec.at(iCutList))->At(iMCSignal))->At(29))->At(27)))->Fill(values[AliDielectronVarManager::kP], values[AliDielectronVarManager::kTOFnSigmaPro]);
    if(TMath::Abs(mcTrack->PdgCode()) == 321)   (dynamic_cast<TH2D *>(((TList*)((TList*)((TList*)fTrackCutListVecSec.at(iCutList))->At(iMCSignal))->At(29))->At(28)))->Fill(values[AliDielectronVarManager::kP], values[AliDielectronVarManager::kTOFnSigmaPro]);
    if(TMath::Abs(mcTrack->PdgCode()) == 2212)  (dynamic_cast<TH2D *>(((TList*)((TList*)((TList*)fTrackCutListVecSec.at(iCutList))->At(iMCSignal))->At(29))->At(29)))->Fill(values[AliDielectronVarManager::kP], values[AliDielectronVarManager::kTOFnSigmaPro]);

}

void    AliAnalysisTaskEtaReconstruction::FillPairHistograms_Secondary(AliDielectronPair* pair, int iCutList, int iMCSignal){
                                                                                // if(fdebug) std::cout << __LINE__ << " DEBUG_AnalysisTask: Called FillPairHistograms_Secondary" << std::endl;
  Double_t values[AliDielectronVarManager::kNMaxValues]={0.};
  AliDielectronVarManager::SetFillMap(fUsedVars); // currently filled manually in the constructor of this task.

  AliDielectronVarManager::Fill(pair, values);
    (dynamic_cast<TH1D *>(((TList*)((TList*)fPairCutListVecSec.at(iCutList))->At(iMCSignal))->At(0)))->Fill(values[AliDielectronVarManager::kCosPointingAngle]);
    (dynamic_cast<TH1D *>(((TList*)((TList*)fPairCutListVecSec.at(iCutList))->At(iMCSignal))->At(1)))->Fill(values[AliDielectronVarManager::kChi2NDF]);
    (dynamic_cast<TH1D *>(((TList*)((TList*)fPairCutListVecSec.at(iCutList))->At(iMCSignal))->At(2)))->Fill(values[AliDielectronVarManager::kLegDist]);
    (dynamic_cast<TH1D *>(((TList*)((TList*)fPairCutListVecSec.at(iCutList))->At(iMCSignal))->At(3)))->Fill(values[AliDielectronVarManager::kR]);
    (dynamic_cast<TH1D *>(((TList*)((TList*)fPairCutListVecSec.at(iCutList))->At(iMCSignal))->At(4)))->Fill(values[AliDielectronVarManager::kPsiPair]);
    (dynamic_cast<TH1D *>(((TList*)((TList*)fPairCutListVecSec.at(iCutList))->At(iMCSignal))->At(5)))->Fill(values[AliDielectronVarManager::kM]);
    (dynamic_cast<TH1D *>(((TList*)((TList*)fPairCutListVecSec.at(iCutList))->At(iMCSignal))->At(6)))->Fill(values[AliDielectronVarManager::kArmPt]);
    (dynamic_cast<TH1D *>(((TList*)((TList*)fPairCutListVecSec.at(iCutList))->At(iMCSignal))->At(7)))->Fill(values[AliDielectronVarManager::kArmAlpha]);
    (dynamic_cast<TH1D *>(((TList*)((TList*)fPairCutListVecSec.at(iCutList))->At(iMCSignal))->At(8)))->Fill(values[AliDielectronVarManager::kPt]);
    (dynamic_cast<TH2D *>(((TList*)((TList*)fPairCutListVecSec.at(iCutList))->At(iMCSignal))->At(9)))->Fill(values[AliDielectronVarManager::kArmAlpha], values[AliDielectronVarManager::kArmPt]);
}


// ############################################################################
// ############################################################################
AliAnalysisTaskEtaReconstruction::Particle AliAnalysisTaskEtaReconstruction::CreateParticle(AliVParticle* mcPart1){
  double  pt1      = mcPart1->Pt();
  double  eta1     = mcPart1->Eta();
  double  phi1     = mcPart1->Phi();
  short   charge1  = mcPart1->Charge();
  Particle part(pt1, eta1, phi1, charge1);


  return part;
}


// ############################################################################
// ############################################################################
Bool_t AliAnalysisTaskEtaReconstruction::CheckIfOneIsTrue(std::vector<Bool_t>& vec){
  bool min_one_is_true = kFALSE;
  unsigned int size = vec.size();
  for (unsigned int i = 0; i < size; ++i){
    if (vec[i] == kTRUE) {min_one_is_true = kTRUE; break;}
  }
  return min_one_is_true;
}



// ############################################################################
// ############################################################################
void AliAnalysisTaskEtaReconstruction::SetBinsLinear(const std::string var, const double min, const double max, const unsigned int steps){
  if      (var == "pt")     fPtBins.clear();
  else if (var == "eta")    fEtaBins.clear();
  else if (var == "phi")    fPhiBins.clear();
  else if (var == "theta")  fThetaBins.clear();
  else if (var == "ptDelta_reso")fResolutionDeltaPtBins.clear();
  else if (var == "ptRel_reso")  fResolutionRelPtBins.clear();
  else if (var == "eta_reso")    fResolutionEtaBins.clear();
  else if (var == "phi_reso")    fResolutionPhiBins.clear();
  else if (var == "theta_reso")  fResolutionThetaBins.clear();
  else if (var == "mass")   fMassBins.clear();
  else if (var == "pairpt") fPairPtBins.clear();

  const double stepSize = (max - min) / steps;
  for (unsigned int i = 0; i < steps+1; ++i){
    if      (var == "pt")     fPtBins.push_back(i * stepSize + min);
    else if (var == "eta")    fEtaBins.push_back(i * stepSize + min);
    else if (var == "phi")    fPhiBins.push_back(i * stepSize + min);
    else if (var == "theta")  fThetaBins.push_back(i * stepSize + min);
    else if (var == "ptDelta_reso")fResolutionDeltaPtBins.push_back(i * stepSize + min);
    else if (var == "ptRel_reso")  fResolutionRelPtBins.push_back(i * stepSize + min);
    else if (var == "eta_reso")    fResolutionEtaBins.push_back(i * stepSize + min);
    else if (var == "phi_reso")    fResolutionPhiBins.push_back(i * stepSize + min);
    else if (var == "theta_reso")  fResolutionThetaBins.push_back(i * stepSize + min);
    else if (var == "mass")   fMassBins.push_back(i * stepSize + min);
    else if (var == "pairpt") fPairPtBins.push_back(i * stepSize + min);
  }
}


// ############################################################################
// ############################################################################
void AliAnalysisTaskEtaReconstruction::CheckSinglePrimaryLegMCsignals(std::vector<Bool_t>& vec, const int tracklabel){
  for (unsigned int i = 0; i < fSinglePrimaryLegMCSignal.size(); ++i){
    vec.at(i) = AliDielectronMC::Instance()->IsMCTruth(tracklabel, &(fSinglePrimaryLegMCSignal[i]), 1);
  }
}

void AliAnalysisTaskEtaReconstruction::CheckSingleSecondaryLegMCsignals(std::vector<Bool_t>& vec, const int tracklabel){
  for (unsigned int i = 0; i < fSingleSecondaryLegMCSignal.size(); ++i){
    vec.at(i) = AliDielectronMC::Instance()->IsMCTruth(tracklabel, &(fSingleSecondaryLegMCSignal[i]), 1);
  }
}

// ############################################################################
// ############################################################################
void AliAnalysisTaskEtaReconstruction::CreateSupportHistos()
{
  Printf(" Now running: CreateSupportHistos()");

  fOutputListSupportHistos = new TList();
  fOutputListSupportHistos->SetName("Support");
  fOutputListSupportHistos->SetOwner();


  for (unsigned int list_i = 0; list_i < fTrackCuts_primary_standard.size(); ++list_i){
    TList* list1 = new TList();
    list1->SetName(Form("PrimaryTrackCuts: %s",fTrackCuts_primary_standard.at(list_i)->GetName()));
    list1->SetOwner();
    for (unsigned int i = 0; i < fSinglePrimaryLegMCSignal.size(); ++i){
      TList* list_temp = new TList();
      list_temp->SetName(Form("PrimaryTrackCuts: %s", fSinglePrimaryLegMCSignal.at(i).GetName()));
      list_temp->SetOwner();

      // Track variables
      TH1D* hPt_prim     = new TH1D("Pt","Pt ;Pt [GeV];#tracks",160,0.,8.);//,AliDielectronVarManager::kPt);

      // PID
      TH2D* hITSnSigmaEle_P_prim = new TH2D("ITSnSigmaPrimEle_P","ITS number of sigmas, primary electrons ;P [GeV/c];ITS number of sigmas ", 160,0.,8.,100,-5.,5.);//.,AliDielectronVarManager::kP,AliDielectronVarManager::kITSnSigmaEle,makeLogx);
      TH2D* hTPCnSigmaEle_P_prim = new TH2D("TPCnSigmaEle_P","TPC number of sigmas, primary electrons ;PIn (pTPC) [GeV/c];TPC number of sigmas Electrons", 160,0.,8.,100,-5.,5.);//.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCnSigmaEle,makeLogx);
      TH2D* hTOFnSigmaEle_P_prim = new TH2D("TOFnSigmaEle_P","TOF number of sigmas, primary electrons ;PIn (pTPC) [GeV/c];TOF number of sigmas Electrons", 160,0.,8.,100,-5.,5.);//,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTOFnSigmaEle,makeLogx);

      // Track kinematic
      TH1D* hEta_prim = new TH1D("Eta","Eta ; Eta;#tracks", 200, -2, 2);//,AliDielectronVarManager::kEta);
      TH1D* hPhi_prim = new TH1D("Phi","Phi ; Phi;#tracks", 320, 0., 6.4);//,AliDielectronVarManager::kPhi);
      TH2D* hEta_Phi_prim = new TH2D("Eta_Phi","Eta Phi Map, primary electrons ; Eta; Phi", 100, -1, 1, 320, 0, 6.4);//,AliDielectronVarManager::kEta,AliDielectronVarManager::kPhi);

      // Quality
      TH1D* hTPCcrossedRowsOverFindable_prim = new TH1D("TPCcrossedRowsOverFindable","Number of Crossed Rows TPC over Findable ;TPC crossed rows over findable;#tracks",120,0.,1.2);//,AliDielectronVarManager::kNFclsTPCfCross);
      TH1D* hTPCcrossedRows_prim = new TH1D("TPCcrossedRows","Number of Crossed Rows TPC ;TPC crossed rows;#tracks",160,-0.5,159.5);//,AliDielectronVarManager::kNFclsTPCr);
      TH1D* hTPCnCls_prim = new TH1D("TPCnCls","Number of Clusters TPC ;TPC number clusters;#tracks",160,-0.5,159.5);//,AliDielectronVarManager::kNclsTPC);
      TH1D* hITSnCls_prim = new TH1D("ITSnCls","Number of Clusters ITS ;ITS number clusters;#tracks",10,-0.5,9.5);//,AliDielectronVarManager::kNclsITS);
      TH1D* hTPCchi2_prim = new TH1D("TPCchi2","TPC chi2 per Cluster ;TPC chi2/Cl;#tracks",100,0.,10.);//,AliDielectronVarManager::kTPCchi2Cl);
      TH1D* hITSchi2_prim = new TH1D("ITSchi2","ITS chi2 per Cluster ;ITS chi2/Cl;#tracks",100,0.,10.);//,AliDielectronVarManager::kITSchi2Cl);
      TH1D* hTPCnClsS_prim = new TH1D("TPCnClsS",";TPC number of shared clusters ;#tracks",160,-0.5,159.5);//,AliDielectronVarManager::kNclsSTPC);
      TH1D* hITSnClsS_prim = new TH1D("ITSnClsS",";ITS number of shared clusters ;#tracks",7,-0.5,6.5);//,AliDielectronVarManager::kNclsSITS);
      TH1D* hNclsSFracTPC_prim = new TH1D("NclsSFracTPC","Fraction of shared clusters assigned in the TPC ;TPC fraction of shared clusters;#tracks",120,0,1.2);//.,AliDielectronVarManager::kNclsSFracTPC);
      TH1D* hNclsSFracITS_prim = new TH1D("NclsSFracITS","Fraction of shared clusters assigned in the ITS ;ITS fraction of shared clusters;#tracks",120,0,1.2);//.,AliDielectronVarManager::kNclsSFracITS);
      TH1D* hTPCclsDiff_prim = new TH1D("TPCclsDiff","TPC cluster difference ;N_{d#it{E}/d#it{x} points}^{TPC} - N_{cls}^{TPC};#tracks",100,-80,20);//.,AliDielectronVarManager::kTPCclsDiff);
      TH1D* hTPCsignalN_prim = new TH1D("TPCsignalN","Number of PID Clusters TPC ;N_{d#it{E}/d#it{x} points}^{TPC};#tracks",160,-0.5,159.5);//.,AliDielectronVarManager::kTPCsignalN); //kNclsTPCdEdx
      TH2D* hTPCcrossedRows_TPCnCls_prim = new TH2D("TPCcrossedRows_TPCnCls","TPC crossed rows vs TPC number clusters, primary electrons ;TPC number clusters;TPC crossed rows;#tracks",
                                               160,-0.5,159.5,160,-0.5,159.5);//,AliDielectronVarManager::kNclsTPC,AliDielectronVarManager::kNFclsTPCr);
      TH2D* hTPCcrossedRows_Pt_prim = new TH2D("TPCcrossedRows_Pt","TPC crossed rows vs Pt, primary electrons ;Pt [GeV];TPC crossed rows",
                                          160,0.,8.,160,-0.5,159.5);//,AliDielectronVarManager::kPt,AliDielectronVarManager::kNFclsTPCr);

      TH1D* hPDGCode_prim = new TH1D("PDGCode","PDGCode ;PDG Code",10001, -5000, 5000);//.,AliDielectronVarManager::kTPCsignalN); //kNclsTPCdEdx

      TH1D* hPDGCodeMother_prim = new TH1D("PDGCodeMother","PDGCodeMother ;Mother PDG Code",10001, -5000, 5000);//.,AliDielectronVarManager::kTPCsignalN); //kNclsTPCdEdx
      // TH2D* hPDGCode_PDGCodeMother = new TH2D("PDGCode_PDGCodeMother",";PDG code;PDG code Mother",
      // 10001,-5000,5000,10001,-5000,5000);//,AliDielectronVarManager::kPt,AliDielectronVarManager::kNFclsTPCr);
      TH1D* hMCGenCode_prim = new TH1D("MCGenerator","MCGenerator ;MC Generators;#tracks",1, 0, 0);//.,AliDielectronVarManager::kTPCsignalN); //kNclsTPCdEdx
      TH1D* hImpactParXY_prim = new TH1D("ImpactParXY","Impact parameter XY ;Impact parameter XY;#tracks",200, -5., 5.);//.,AliDielectronVarManager::kTPCsignalN); //kImpactParXY_prim
      TH1D* hImpactParZ_prim = new TH1D("ImpactParZ","Impact parameter Z ;Impact parameter Z;#tracks",200, -5., 5.);//.,AliDielectronVarManager::kTPCsignalN); //kImpactParZ_prim
      TH1D* hM_prim = new TH1D("M","Mass ;m ;#tracks",200, 0., 1.); // kM

      // TH1D* hCosPointingAngle_prim = new TH1D("CosPointingAngle","Cosine of the pointing angle;#tracks",200, -0.5, 6.5); // kCosPointingAngle
      // TH1D* hChi2NDF_prim = new TH1D("Chi2NDF","Chi^2/NDF ;Chi^2/NDF;#tracks",200, -5., 5.); // kChi2NDF (NDF = number degree of freedom)
      // TH1D* hLegDist_prim = new TH1D("LegDist","Leg Distance ;Leg Distance;#tracks",200, -0.5, 5.); // kLegDist
      // TH1D* hR_prim = new TH1D("R","Distance to the origin ;r;#tracks",200, -0.5, 5.); // kR
      // TH1D* hPsiPair_prim = new TH1D("PsiPair","phi in mother's rest frame in Collins-Soper picture ;PsiPair;#tracks",200, -5., 5.); // kPsiPair
      // TH1D* hArmPt_prim = new TH1D("ArmPt","Armenteros-Podolanski pt ;Pt [GeV];#tracks",160,0.,1.); // kArmPt
      // TH1D* hArmAlpha_prim = new TH1D("ArmAlpha","Armenteros-Podolanski alpha ;alpha;#tracks",200, -1., 1.); // kArmAlpha

      list_temp->AddAt(hPt_prim,     0);
      list_temp->AddAt(hITSnSigmaEle_P_prim, 1);
      list_temp->AddAt(hTPCnSigmaEle_P_prim, 2);
      list_temp->AddAt(hTOFnSigmaEle_P_prim, 3);
      list_temp->AddAt(hEta_prim, 4);
      list_temp->AddAt(hPhi_prim, 5);
      list_temp->AddAt(hEta_Phi_prim, 6);
      list_temp->AddAt(hTPCcrossedRowsOverFindable_prim, 7);
      list_temp->AddAt(hTPCcrossedRows_prim, 8);
      list_temp->AddAt(hTPCnCls_prim, 9);
      list_temp->AddAt(hITSnCls_prim, 10);
      list_temp->AddAt(hTPCchi2_prim, 11);
      list_temp->AddAt(hITSchi2_prim, 12);
      list_temp->AddAt(hTPCnClsS_prim, 13);
      list_temp->AddAt(hITSnClsS_prim, 14);
      list_temp->AddAt(hNclsSFracTPC_prim, 15);
      list_temp->AddAt(hNclsSFracITS_prim, 16);
      list_temp->AddAt(hTPCclsDiff_prim, 17);
      list_temp->AddAt(hTPCsignalN_prim, 18);
      list_temp->AddAt(hTPCcrossedRows_TPCnCls_prim, 19);
      list_temp->AddAt(hTPCcrossedRows_Pt_prim, 20);
      list_temp->AddAt(hPDGCode_prim, 21);
      // list_temp->AddAt(hPDGCode_PDGCodeMother, 21);
      list_temp->AddAt(hPDGCodeMother_prim, 22);
      list_temp->AddAt(hMCGenCode_prim, 23);
      list_temp->AddAt(hImpactParXY_prim, 24);
      list_temp->AddAt(hImpactParZ_prim, 25);
      list_temp->AddAt(hM_prim,26);
      // list_temp->AddAt(hCosPointingAngle_prim,27);
      // list_temp->AddAt(hChi2NDF_prim,28);
      // list_temp->AddAt(hLegDist_prim,29);
      // list_temp->AddAt(hR_prim,30);
      // list_temp->AddAt(hPsiPair_prim,31);
      // list_temp->AddAt(hArmPt_prim,32);
      // list_temp->AddAt(hArmAlpha_prim,33);

      // fOutputListSupportHistos->Add(list_temp);
      list1->Add(list_temp);
    }
  fTrackCutListVecPrim.push_back(list1);
  fOutputListSupportHistos->Add(list1);
  }

  for (unsigned int list_i = 0; list_i < fPairCuts_secondary_standard.size(); ++list_i){
    TList* list2 = new TList();
    list2->SetName(Form("SecondaryTrackCuts: %s",fPairCuts_secondary_standard.at(list_i)->GetName()));
    list2->SetOwner();
    for (unsigned int i = 0; i < fSecondaryPairMCSignal.size(); ++i){
      TList* list_temp = new TList();
      list_temp->SetName(Form("SecondaryTrackCuts: %s", fSecondaryPairMCSignal.at(i).GetName()));
      list_temp->SetOwner();

      TList* list_nSigmaITS = new TList();
      list_nSigmaITS->SetName("ITS_nSigmas");
      list_nSigmaITS->SetOwner();

      TList* list_nSigmaTPC = new TList();
      list_nSigmaTPC->SetName("TPC_nSigmas");
      list_nSigmaTPC->SetOwner();

      TList* list_nSigmaTOF = new TList();
      list_nSigmaTOF->SetName("TOF_nSigmas");
      list_nSigmaTOF->SetOwner();

    TH1D* hPt_sec      = new TH1D("Pt" ,"Pt ;Pt [GeV];#tracks",640,0.,8.);//,AliDielectronVarManager::kPt);
    TH2D* hITSnSigmaEle_P_sec  = new TH2D("ITSnSigmaSecEle_P","ITS number of sigmas, secondary particles ;P [GeV/c];ITS number of sigmas Electrons", 160,0.,8.,100,-5.,5.);//.,AliDielectronVarManager::kP,AliDielectronVarManager::kITSnSigmaEle,makeLogx);
    TH2D* hTPCnSigmaEle_P_sec = new TH2D("TPCnSigmaEle_P","TPC number of sigmas, secondary particles ;PIn (pTPC) [GeV/c];TPC number of sigmas Electrons", 160,0.,8.,100,-5.,5.);//.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCnSigmaEle,makeLogx);
    TH2D* hTOFnSigmaEle_P_sec = new TH2D("TOFnSigmaEle_P","TOF number of sigmas, secondary particles ;P [GeV/c];TOF number of sigmas Electrons", 160,0.,8.,100,-5.,5.);//,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTOFnSigmaEle,makeLogx);
    TH1D* hEta_sec = new TH1D("Eta","Eta ; Eta;#tracks", 200, -2, 2);//,AliDielectronVarManager::kEta);
    TH1D* hPhi_sec = new TH1D("Phi","Phi ; Phi;#tracks", 320, 0., 6.4);//,AliDielectronVarManager::kPhi);
    TH2D* hEta_Phi_sec = new TH2D("Eta_Phi","Eta Phi Map, secondary electrons ; Eta; Phi", 100, -1, 1, 320, 0, 6.4);//,AliDielectronVarManager::kEta,AliDielectronVarManager::kPhi);
    TH1D* hTPCcrossedRowsOverFindable_sec = new TH1D("TPCcrossedRowsOverFindable","Number of Crossed Rows TPC over Findable ;TPC crossed rows over findable;#tracks",120,0.,1.2);//,AliDielectronVarManager::kNFclsTPCfCross);
    TH1D* hTPCcrossedRows_sec = new TH1D("TPCcrossedRows","Number of Crossed Rows TPC ;TPC crossed rows;#tracks",160,-0.5,159.5);//,AliDielectronVarManager::kNFclsTPCr);
    TH1D* hTPCnCls_sec = new TH1D("TPCnCls","Number of Clusters TPC ;TPC number clusters;#tracks",160,-0.5,159.5);//,AliDielectronVarManager::kNclsTPC);
    TH1D* hITSnCls_sec = new TH1D("ITSnCls","Number of Clusters ITS ;ITS number clusters;#tracks",10,-0.5,9.5);//,AliDielectronVarManager::kNclsITS);
    TH1D* hTPCchi2_sec = new TH1D("TPCchi2","TPC chi2 per Cluster ;TPC chi2/Cl;#tracks",100,0.,10.);//,AliDielectronVarManager::kTPCchi2Cl);
    TH1D* hITSchi2_sec = new TH1D("ITSchi2","ITS chi2 per Cluster ;ITS chi2/Cl;#tracks",100,0.,10.);//,AliDielectronVarManager::kITSchi2Cl);
    TH1D* hTPCnClsS_sec = new TH1D("TPCnClsS",";TPC number of shared clusters ;#tracks",160,-0.5,159.5);//,AliDielectronVarManager::kNclsSTPC);
    TH1D* hITSnClsS_sec = new TH1D("ITSnClsS",";ITS number of shared clusters ;#tracks",7,-0.5,6.5);//,AliDielectronVarManager::kNclsSITS);
    TH1D* hNclsSFracTPC_sec = new TH1D("NclsSFracTPC","Fraction of shared clusters assigned in the TPC ;TPC fraction of shared clusters;#tracks",120,0,1.2);//.,AliDielectronVarManager::kNclsSFracTPC);
    TH1D* hNclsSFracITS_sec = new TH1D("NclsSFracITS","Fraction of shared clusters assigned in the ITS ;ITS fraction of shared clusters;#tracks",120,0,1.2);//.,AliDielectronVarManager::kNclsSFracITS);
    TH1D* hTPCclsDiff_sec = new TH1D("TPCclsDiff","TPC cluster difference ;N_{d#it{E}/d#it{x} points}^{TPC} - N_{cls}^{TPC};#tracks",100,-80,20);//.,AliDielectronVarManager::kTPCclsDiff);
    TH1D* hTPCsignalN_sec = new TH1D("TPCsignalN","Number of PID Clusters TPC ;N_{d#it{E}/d#it{x} points}^{TPC};#tracks",160,-0.5,159.5);//.,AliDielectronVarManager::kTPCsignalN); //kNclsTPCdEdx
    TH2D* hTPCcrossedRows_TPCnCls_sec = new TH2D("TPCcrossedRows_TPCnCls","TPC crossed rows vs TPC number clusters, secondary electrons ;TPC number clusters;TPC crossed rows",
                                             160,-0.5,159.5,160,-0.5,159.5);//,AliDielectronVarManager::kNclsTPC,AliDielectronVarManager::kNFclsTPCr);
    TH2D* hTPCcrossedRows_Pt_sec = new TH2D("TPCcrossedRows_Pt","TPC crossed rows vs Pt, secondary electrons ;Pt [GeV];TPC crossed rows",
                                        160,0.,8.,160,-0.5,159.5);//,AliDielectronVarManager::kPt,AliDielectronVarManager::kNFclsTPCr);
    TH1D* hPDGCode_sec = new TH1D("PDGCode","PDGCode ;PDG Code;#tracks",10001, -5000, 5000);//.,AliDielectronVarManager::kTPCsignalN); //kNclsTPCdEdx
    TH1D* hPDGCodeMother_sec = new TH1D("PDGCodeMother","PDGCodeMother ;Mother PDG Code;#tracks",10001, -5000, 5000);//.,AliDielectronVarManager::kTPCsignalN); //kNclsTPCdEdx
    TH1D* hMCGenCode_sec = new TH1D("MCGenerator","MCGenerator ;MC Generators;#tracks",1, 0, 0);//.,AliDielectronVarManager::kTPCsignalN); //kNclsTPCdEdx
    TH1D* hImpactParXY_sec = new TH1D("ImpactParXY","Impact parameter XY ;Impact parameter XY;#tracks",200, -5., 5.);//.,AliDielectronVarManager::kTPCsignalN); //kImpactParXY_sec
    TH1D* hImpactParZ_sec = new TH1D("ImpactParZ","Impact parameter Z ;Impact parameter Z;#tracks",200, -5., 5.);//.,AliDielectronVarManager::kTPCsignalN); //kImpactParZ_sec
    TH1D* hM_single_sec = new TH1D("SingleMass","single Mass ;mass ;#tracks",200, 0., 1.); // kM

    TH2D* hITSnSigmaEle_P_sec_all  = new TH2D("ITSnSigmaEle_P_true_all ","ITS number of sigmas, filled with all kind  of secondaries ;P [GeV/c];ITS number of sigmas Electrons", 160,0.,8.,100,-5.,5.);
    TH2D* hITSnSigmaEle_P_sec_elec = new TH2D("ITSnSigmaEle_P_true_elec","ITS number of sigmas, filled with secondary elec ;P [GeV/c];ITS number of sigmas Electrons", 160,0.,8.,100,-5.,5.);
    TH2D* hITSnSigmaEle_P_sec_muon = new TH2D("ITSnSigmaEle_P_true_muon","ITS number of sigmas, filled with secondary muon ;P [GeV/c];ITS number of sigmas Electrons", 160,0.,8.,100,-5.,5.);
    TH2D* hITSnSigmaEle_P_sec_pion = new TH2D("ITSnSigmaEle_P_true_pion","ITS number of sigmas, filled with secondary pion ;P [GeV/c];ITS number of sigmas Electrons", 160,0.,8.,100,-5.,5.);
    TH2D* hITSnSigmaEle_P_sec_kaon = new TH2D("ITSnSigmaEle_P_true_kaon","ITS number of sigmas, filled with secondary kaon ;P [GeV/c];ITS number of sigmas Electrons", 160,0.,8.,100,-5.,5.);
    TH2D* hITSnSigmaEle_P_sec_prot = new TH2D("ITSnSigmaEle_P_true_prot","ITS number of sigmas, filled with secondary prot ;P [GeV/c];ITS number of sigmas Electrons", 160,0.,8.,100,-5.,5.);

    TH2D* hITSnSigmaMuon_P_sec_all  = new TH2D("ITSnSigmaMuon_P_true_all ","ITS number of sigmas, filled with all kind  of secondaries ;P [GeV/c];ITS number of sigmas Muon", 160,0.,8.,100,-5.,5.);
    TH2D* hITSnSigmaMuon_P_sec_elec = new TH2D("ITSnSigmaMuon_P_true_elec","ITS number of sigmas, filled with secondary elec ;P [GeV/c];ITS number of sigmas Muon", 160,0.,8.,100,-5.,5.);
    TH2D* hITSnSigmaMuon_P_sec_muon = new TH2D("ITSnSigmaMuon_P_true_muon","ITS number of sigmas, filled with secondary muon ;P [GeV/c];ITS number of sigmas Muon", 160,0.,8.,100,-5.,5.);
    TH2D* hITSnSigmaMuon_P_sec_pion = new TH2D("ITSnSigmaMuon_P_true_pion","ITS number of sigmas, filled with secondary pion ;P [GeV/c];ITS number of sigmas Muon", 160,0.,8.,100,-5.,5.);
    TH2D* hITSnSigmaMuon_P_sec_kaon = new TH2D("ITSnSigmaMuon_P_true_kaon","ITS number of sigmas, filled with secondary kaon ;P [GeV/c];ITS number of sigmas Muon", 160,0.,8.,100,-5.,5.);
    TH2D* hITSnSigmaMuon_P_sec_prot = new TH2D("ITSnSigmaMuon_P_true_prot","ITS number of sigmas, filled with secondary prot ;P [GeV/c];ITS number of sigmas Muon", 160,0.,8.,100,-5.,5.);

    TH2D* hITSnSigmaPion_P_sec_all  = new TH2D("ITSnSigmaPion_P_true_all ","ITS number of sigmas, filled with all kind  of secondaries ;P [GeV/c];ITS number of sigmas Pion", 160,0.,8.,100,-5.,5.);
    TH2D* hITSnSigmaPion_P_sec_elec = new TH2D("ITSnSigmaPion_P_true_elec","ITS number of sigmas, filled with secondary elec ;P [GeV/c];ITS number of sigmas Pions", 160,0.,8.,100,-5.,5.);
    TH2D* hITSnSigmaPion_P_sec_muon = new TH2D("ITSnSigmaPion_P_true_muon","ITS number of sigmas, filled with secondary muon ;P [GeV/c];ITS number of sigmas Pions", 160,0.,8.,100,-5.,5.);
    TH2D* hITSnSigmaPion_P_sec_pion = new TH2D("ITSnSigmaPion_P_true_pion","ITS number of sigmas, filled with secondary pion ;P [GeV/c];ITS number of sigmas Pions", 160,0.,8.,100,-5.,5.);
    TH2D* hITSnSigmaPion_P_sec_kaon = new TH2D("ITSnSigmaPion_P_true_kaon","ITS number of sigmas, filled with secondary kaon ;P [GeV/c];ITS number of sigmas Pions", 160,0.,8.,100,-5.,5.);
    TH2D* hITSnSigmaPion_P_sec_prot = new TH2D("ITSnSigmaPion_P_true_prot","ITS number of sigmas, filled with secondary prot ;P [GeV/c];ITS number of sigmas Pions", 160,0.,8.,100,-5.,5.);

    TH2D* hITSnSigmaKaon_P_sec_all  = new TH2D("ITSnSigmaKaon_P_true_all ","ITS number of sigmas, filled with all kind  of secondaries ;P [GeV/c];ITS number of sigmas Kaon", 160,0.,8.,100,-5.,5.);
    TH2D* hITSnSigmaKaon_P_sec_elec = new TH2D("ITSnSigmaKaon_P_true_elec","ITS number of sigmas, filled with secondary elec ;P [GeV/c];ITS number of sigmas Kaon", 160,0.,8.,100,-5.,5.);
    TH2D* hITSnSigmaKaon_P_sec_muon = new TH2D("ITSnSigmaKaon_P_true_muon","ITS number of sigmas, filled with secondary muon ;P [GeV/c];ITS number of sigmas Kaon", 160,0.,8.,100,-5.,5.);
    TH2D* hITSnSigmaKaon_P_sec_pion = new TH2D("ITSnSigmaKaon_P_true_pion","ITS number of sigmas, filled with secondary pion ;P [GeV/c];ITS number of sigmas Kaon", 160,0.,8.,100,-5.,5.);
    TH2D* hITSnSigmaKaon_P_sec_kaon = new TH2D("ITSnSigmaKaon_P_true_kaon","ITS number of sigmas, filled with secondary kaon ;P [GeV/c];ITS number of sigmas Kaon", 160,0.,8.,100,-5.,5.);
    TH2D* hITSnSigmaKaon_P_sec_prot = new TH2D("ITSnSigmaKaon_P_true_prot","ITS number of sigmas, filled with secondary prot ;P [GeV/c];ITS number of sigmas Kaon", 160,0.,8.,100,-5.,5.);

    TH2D* hITSnSigmaProton_P_sec_all  = new TH2D("ITSnSigmaProton_P_true_all ","ITS number of sigmas, filled with all kind  of secondaries ;P [GeV/c];ITS number of sigmas Proton", 160,0.,8.,100,-5.,5.);
    TH2D* hITSnSigmaProton_P_sec_elec = new TH2D("ITSnSigmaProton_P_true_elec","ITS number of sigmas, filled with secondary elec ;P [GeV/c];ITS number of sigmas Proton", 160,0.,8.,100,-5.,5.);
    TH2D* hITSnSigmaProton_P_sec_muon = new TH2D("ITSnSigmaProton_P_true_muon","ITS number of sigmas, filled with secondary muon ;P [GeV/c];ITS number of sigmas Proton", 160,0.,8.,100,-5.,5.);
    TH2D* hITSnSigmaProton_P_sec_pion = new TH2D("ITSnSigmaProton_P_true_pion","ITS number of sigmas, filled with secondary pion ;P [GeV/c];ITS number of sigmas Proton", 160,0.,8.,100,-5.,5.);
    TH2D* hITSnSigmaProton_P_sec_kaon = new TH2D("ITSnSigmaProton_P_true_kaon","ITS number of sigmas, filled with secondary kaon ;P [GeV/c];ITS number of sigmas Proton", 160,0.,8.,100,-5.,5.);
    TH2D* hITSnSigmaProton_P_sec_prot = new TH2D("ITSnSigmaProton_P_true_prot","ITS number of sigmas, filled with secondary prot ;P [GeV/c];ITS number of sigmas Proton", 160,0.,8.,100,-5.,5.);


    TH2D* hTPCnSigmaEle_P_sec_all  = new TH2D("TPCnSigmaEle_P_true_all ","TPC number of sigmas, filled with all kind  of secondaries ;PIn (pTPC) [GeV/c];TPC number of sigmas Electrons", 160,0.,8.,100,-5.,5.);
    TH2D* hTPCnSigmaEle_P_sec_elec = new TH2D("TPCnSigmaEle_P_true_elec","TPC number of sigmas, filled with secondary elec ;PIn (pTPC) [GeV/c];TPC number of sigmas Electrons", 160,0.,8.,100,-5.,5.);//.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCnSigmaEle,makeLogx);
    TH2D* hTPCnSigmaEle_P_sec_muon = new TH2D("TPCnSigmaEle_P_true_muon","TPC number of sigmas, filled with secondary muon ;PIn (pTPC) [GeV/c];TPC number of sigmas Electrons", 160,0.,8.,100,-5.,5.);//.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCnSigmaEle,makeLogx);
    TH2D* hTPCnSigmaEle_P_sec_pion = new TH2D("TPCnSigmaEle_P_true_pion","TPC number of sigmas, filled with secondary pion ;PIn (pTPC) [GeV/c];TPC number of sigmas Electrons", 160,0.,8.,100,-5.,5.);//.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCnSigmaEle,makeLogx);
    TH2D* hTPCnSigmaEle_P_sec_kaon = new TH2D("TPCnSigmaEle_P_true_kaon","TPC number of sigmas, filled with secondary kaon ;PIn (pTPC) [GeV/c];TPC number of sigmas Electrons", 160,0.,8.,100,-5.,5.);//.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCnSigmaEle,makeLogx);
    TH2D* hTPCnSigmaEle_P_sec_prot = new TH2D("TPCnSigmaEle_P_true_prot","TPC number of sigmas, filled with secondary prot ;PIn (pTPC) [GeV/c];TPC number of sigmas Electrons", 160,0.,8.,100,-5.,5.);//.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCnSigmaEle,makeLogx);

    TH2D* hTPCnSigmaMuon_P_sec_all  = new TH2D("TPCnSigmaMuon_P_true_all ","TPC number of sigmas, filled with all kind  of secondaries ;PIn (pTPC) [GeV/c];TPC number of sigmas Muon", 160,0.,8.,100,-5.,5.);
    TH2D* hTPCnSigmaMuon_P_sec_elec = new TH2D("TPCnSigmaMuon_P_true_elec","TPC number of sigmas, filled with secondary elec ;PIn (pTPC) [GeV/c];TPC number of sigmas Muon", 160,0.,8.,100,-5.,5.);//.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCnSigmaEle,makeLogx);
    TH2D* hTPCnSigmaMuon_P_sec_muon = new TH2D("TPCnSigmaMuon_P_true_muon","TPC number of sigmas, filled with secondary muon ;PIn (pTPC) [GeV/c];TPC number of sigmas Muon", 160,0.,8.,100,-5.,5.);//.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCnSigmaEle,makeLogx);
    TH2D* hTPCnSigmaMuon_P_sec_pion = new TH2D("TPCnSigmaMuon_P_true_pion","TPC number of sigmas, filled with secondary pion ;PIn (pTPC) [GeV/c];TPC number of sigmas Muon", 160,0.,8.,100,-5.,5.);//.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCnSigmaEle,makeLogx);
    TH2D* hTPCnSigmaMuon_P_sec_kaon = new TH2D("TPCnSigmaMuon_P_true_kaon","TPC number of sigmas, filled with secondary kaon ;PIn (pTPC) [GeV/c];TPC number of sigmas Muon", 160,0.,8.,100,-5.,5.);//.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCnSigmaEle,makeLogx);
    TH2D* hTPCnSigmaMuon_P_sec_prot = new TH2D("TPCnSigmaMuon_P_true_prot","TPC number of sigmas, filled with secondary prot ;PIn (pTPC) [GeV/c];TPC number of sigmas Muon", 160,0.,8.,100,-5.,5.);//.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCnSigmaEle,makeLogx);

    TH2D* hTPCnSigmaPion_P_sec_all  = new TH2D("TPCnSigmaPion_P_true_all ","TPC number of sigmas, filled with all kind  of secondaries ;PIn (pTPC) [GeV/c];TPC number of sigmas Pion", 160,0.,8.,100,-5.,5.);
    TH2D* hTPCnSigmaPion_P_sec_elec = new TH2D("TPCnSigmaPion_P_true_elec","TPC number of sigmas, filled with secondary elec ;PIn (pTPC) [GeV/c];TPC number of sigmas Pions", 160,0.,8.,100,-5.,5.);//.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCnSigmaEle,makeLogx);
    TH2D* hTPCnSigmaPion_P_sec_muon = new TH2D("TPCnSigmaPion_P_true_muon","TPC number of sigmas, filled with secondary muon ;PIn (pTPC) [GeV/c];TPC number of sigmas Pions", 160,0.,8.,100,-5.,5.);//.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCnSigmaEle,makeLogx);
    TH2D* hTPCnSigmaPion_P_sec_pion = new TH2D("TPCnSigmaPion_P_true_pion","TPC number of sigmas, filled with secondary pion ;PIn (pTPC) [GeV/c];TPC number of sigmas Pions", 160,0.,8.,100,-5.,5.);//.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCnSigmaEle,makeLogx);
    TH2D* hTPCnSigmaPion_P_sec_kaon = new TH2D("TPCnSigmaPion_P_true_kaon","TPC number of sigmas, filled with secondary kaon ;PIn (pTPC) [GeV/c];TPC number of sigmas Pions", 160,0.,8.,100,-5.,5.);//.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCnSigmaEle,makeLogx);
    TH2D* hTPCnSigmaPion_P_sec_prot = new TH2D("TPCnSigmaPion_P_true_prot","TPC number of sigmas, filled with secondary prot ;PIn (pTPC) [GeV/c];TPC number of sigmas Pions", 160,0.,8.,100,-5.,5.);//.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCnSigmaEle,makeLogx);

    TH2D* hTPCnSigmaKaon_P_sec_all  = new TH2D("TPCnSigmaKaon_P_true_all ","TPC number of sigmas, filled with all kind  of secondaries ;PIn (pTPC) [GeV/c];TPC number of sigmas Kaon", 160,0.,8.,100,-5.,5.);
    TH2D* hTPCnSigmaKaon_P_sec_elec = new TH2D("TPCnSigmaKaon_P_true_elec","TPC number of sigmas, filled with secondary elec ;PIn (pTPC) [GeV/c];TPC number of sigmas Kaon", 160,0.,8.,100,-5.,5.);//.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCnSigmaEle,makeLogx);
    TH2D* hTPCnSigmaKaon_P_sec_muon = new TH2D("TPCnSigmaKaon_P_true_muon","TPC number of sigmas, filled with secondary muon ;PIn (pTPC) [GeV/c];TPC number of sigmas Kaon", 160,0.,8.,100,-5.,5.);//.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCnSigmaEle,makeLogx);
    TH2D* hTPCnSigmaKaon_P_sec_pion = new TH2D("TPCnSigmaKaon_P_true_pion","TPC number of sigmas, filled with secondary pion ;PIn (pTPC) [GeV/c];TPC number of sigmas Kaon", 160,0.,8.,100,-5.,5.);//.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCnSigmaEle,makeLogx);
    TH2D* hTPCnSigmaKaon_P_sec_kaon = new TH2D("TPCnSigmaKaon_P_true_kaon","TPC number of sigmas, filled with secondary kaon ;PIn (pTPC) [GeV/c];TPC number of sigmas Kaon", 160,0.,8.,100,-5.,5.);//.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCnSigmaEle,makeLogx);
    TH2D* hTPCnSigmaKaon_P_sec_prot = new TH2D("TPCnSigmaKaon_P_true_prot","TPC number of sigmas, filled with secondary prot ;PIn (pTPC) [GeV/c];TPC number of sigmas Kaon", 160,0.,8.,100,-5.,5.);//.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCnSigmaEle,makeLogx);

    TH2D* hTPCnSigmaProton_P_sec_all  = new TH2D("TPCnSigmaProton_P_true_all ","TPC number of sigmas, filled with all kind  of secondaries ;PIn (pTPC) [GeV/c];TPC number of sigmas Proton", 160,0.,8.,100,-5.,5.);
    TH2D* hTPCnSigmaProton_P_sec_elec = new TH2D("TPCnSigmaProton_P_true_elec","TPC number of sigmas, filled with secondary elec ;PIn (pTPC) [GeV/c];TPC number of sigmas Proton", 160,0.,8.,100,-5.,5.);//.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCnSigmaEle,makeLogx);
    TH2D* hTPCnSigmaProton_P_sec_muon = new TH2D("TPCnSigmaProton_P_true_muon","TPC number of sigmas, filled with secondary muon ;PIn (pTPC) [GeV/c];TPC number of sigmas Proton", 160,0.,8.,100,-5.,5.);//.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCnSigmaEle,makeLogx);
    TH2D* hTPCnSigmaProton_P_sec_pion = new TH2D("TPCnSigmaProton_P_true_pion","TPC number of sigmas, filled with secondary pion ;PIn (pTPC) [GeV/c];TPC number of sigmas Proton", 160,0.,8.,100,-5.,5.);//.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCnSigmaEle,makeLogx);
    TH2D* hTPCnSigmaProton_P_sec_kaon = new TH2D("TPCnSigmaProton_P_true_kaon","TPC number of sigmas, filled with secondary kaon ;PIn (pTPC) [GeV/c];TPC number of sigmas Proton", 160,0.,8.,100,-5.,5.);//.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCnSigmaEle,makeLogx);
    TH2D* hTPCnSigmaProton_P_sec_prot = new TH2D("TPCnSigmaProton_P_true_prot","TPC number of sigmas, filled with secondary prot ;PIn (pTPC) [GeV/c];TPC number of sigmas Proton", 160,0.,8.,100,-5.,5.);//.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCnSigmaEle,makeLogx);


    TH2D* hTOFnSigmaEle_P_sec_all  = new TH2D("TOFnSigmaEle_P_true_all ","TOF number of sigmas, filled with all kind  of secondaries ;P [GeV/c];TOF number of sigmas Electrons", 160,0.,8.,100,-5.,5.);
    TH2D* hTOFnSigmaEle_P_sec_elec = new TH2D("TOFnSigmaEle_P_true_elec","TOF number of sigmas, filled with secondary elec ;P [GeV/c];TOF number of sigmas Electrons", 160,0.,8.,100,-5.,5.);
    TH2D* hTOFnSigmaEle_P_sec_muon = new TH2D("TOFnSigmaEle_P_true_muon","TOF number of sigmas, filled with secondary muon ;P [GeV/c];TOF number of sigmas Electrons", 160,0.,8.,100,-5.,5.);
    TH2D* hTOFnSigmaEle_P_sec_pion = new TH2D("TOFnSigmaEle_P_true_pion","TOF number of sigmas, filled with secondary pion ;P [GeV/c];TOF number of sigmas Electrons", 160,0.,8.,100,-5.,5.);
    TH2D* hTOFnSigmaEle_P_sec_kaon = new TH2D("TOFnSigmaEle_P_true_kaon","TOF number of sigmas, filled with secondary kaon ;P [GeV/c];TOF number of sigmas Electrons", 160,0.,8.,100,-5.,5.);
    TH2D* hTOFnSigmaEle_P_sec_prot = new TH2D("TOFnSigmaEle_P_true_prot","TOF number of sigmas, filled with secondary prot ;P [GeV/c];TOF number of sigmas Electrons", 160,0.,8.,100,-5.,5.);

    TH2D* hTOFnSigmaMuon_P_sec_all  = new TH2D("TOFnSigmaMuon_P_true_all ","TOF number of sigmas, filled with all kind  of secondaries ;P [GeV/c];TOF number of sigmas Muon", 160,0.,8.,100,-5.,5.);
    TH2D* hTOFnSigmaMuon_P_sec_elec = new TH2D("TOFnSigmaMuon_P_true_elec","TOF number of sigmas, filled with secondary elec ;P [GeV/c];TOF number of sigmas Muon", 160,0.,8.,100,-5.,5.);
    TH2D* hTOFnSigmaMuon_P_sec_muon = new TH2D("TOFnSigmaMuon_P_true_muon","TOF number of sigmas, filled with secondary muon ;P [GeV/c];TOF number of sigmas Muon", 160,0.,8.,100,-5.,5.);
    TH2D* hTOFnSigmaMuon_P_sec_pion = new TH2D("TOFnSigmaMuon_P_true_pion","TOF number of sigmas, filled with secondary pion ;P [GeV/c];TOF number of sigmas Muon", 160,0.,8.,100,-5.,5.);
    TH2D* hTOFnSigmaMuon_P_sec_kaon = new TH2D("TOFnSigmaMuon_P_true_kaon","TOF number of sigmas, filled with secondary kaon ;P [GeV/c];TOF number of sigmas Muon", 160,0.,8.,100,-5.,5.);
    TH2D* hTOFnSigmaMuon_P_sec_prot = new TH2D("TOFnSigmaMuon_P_true_prot","TOF number of sigmas, filled with secondary prot ;P [GeV/c];TOF number of sigmas Muon", 160,0.,8.,100,-5.,5.);

    TH2D* hTOFnSigmaPion_P_sec_all  = new TH2D("TOFnSigmaPion_P_true_all ","TOF number of sigmas, filled with all kind  of secondaries ;P [GeV/c];TOF number of sigmas Pion", 160,0.,8.,100,-5.,5.);
    TH2D* hTOFnSigmaPion_P_sec_elec = new TH2D("TOFnSigmaPion_P_true_elec","TOF number of sigmas, filled with secondary elec ;P [GeV/c];TOF number of sigmas Pions", 160,0.,8.,100,-5.,5.);
    TH2D* hTOFnSigmaPion_P_sec_muon = new TH2D("TOFnSigmaPion_P_true_muon","TOF number of sigmas, filled with secondary muon ;P [GeV/c];TOF number of sigmas Pions", 160,0.,8.,100,-5.,5.);
    TH2D* hTOFnSigmaPion_P_sec_pion = new TH2D("TOFnSigmaPion_P_true_pion","TOF number of sigmas, filled with secondary pion ;P [GeV/c];TOF number of sigmas Pions", 160,0.,8.,100,-5.,5.);
    TH2D* hTOFnSigmaPion_P_sec_kaon = new TH2D("TOFnSigmaPion_P_true_kaon","TOF number of sigmas, filled with secondary kaon ;P [GeV/c];TOF number of sigmas Pions", 160,0.,8.,100,-5.,5.);
    TH2D* hTOFnSigmaPion_P_sec_prot = new TH2D("TOFnSigmaPion_P_true_prot","TOF number of sigmas, filled with secondary prot ;P [GeV/c];TOF number of sigmas Pions", 160,0.,8.,100,-5.,5.);

    TH2D* hTOFnSigmaKaon_P_sec_all  = new TH2D("TOFnSigmaKaon_P_true_all ","TOF number of sigmas, filled with all kind  of secondaries ;P [GeV/c];TOF number of sigmas Kaon", 160,0.,8.,100,-5.,5.);
    TH2D* hTOFnSigmaKaon_P_sec_elec = new TH2D("TOFnSigmaKaon_P_true_elec","TOF number of sigmas, filled with secondary elec ;P [GeV/c];TOF number of sigmas Kaon", 160,0.,8.,100,-5.,5.);
    TH2D* hTOFnSigmaKaon_P_sec_muon = new TH2D("TOFnSigmaKaon_P_true_muon","TOF number of sigmas, filled with secondary muon ;P [GeV/c];TOF number of sigmas Kaon", 160,0.,8.,100,-5.,5.);
    TH2D* hTOFnSigmaKaon_P_sec_pion = new TH2D("TOFnSigmaKaon_P_true_pion","TOF number of sigmas, filled with secondary pion ;P [GeV/c];TOF number of sigmas Kaon", 160,0.,8.,100,-5.,5.);
    TH2D* hTOFnSigmaKaon_P_sec_kaon = new TH2D("TOFnSigmaKaon_P_true_kaon","TOF number of sigmas, filled with secondary kaon ;P [GeV/c];TOF number of sigmas Kaon", 160,0.,8.,100,-5.,5.);
    TH2D* hTOFnSigmaKaon_P_sec_prot = new TH2D("TOFnSigmaKaon_P_true_prot","TOF number of sigmas, filled with secondary prot ;P [GeV/c];TOF number of sigmas Kaon", 160,0.,8.,100,-5.,5.);

    TH2D* hTOFnSigmaProton_P_sec_all  = new TH2D("TOFnSigmaProton_P_true_all ","TOF number of sigmas, filled with all kind  of secondaries ;P [GeV/c];TOF number of sigmas Proton", 160,0.,8.,100,-5.,5.);
    TH2D* hTOFnSigmaProton_P_sec_elec = new TH2D("TOFnSigmaProton_P_true_elec","TOF number of sigmas, filled with secondary elec ;P [GeV/c];TOF number of sigmas Proton", 160,0.,8.,100,-5.,5.);
    TH2D* hTOFnSigmaProton_P_sec_muon = new TH2D("TOFnSigmaProton_P_true_muon","TOF number of sigmas, filled with secondary muon ;P [GeV/c];TOF number of sigmas Proton", 160,0.,8.,100,-5.,5.);
    TH2D* hTOFnSigmaProton_P_sec_pion = new TH2D("TOFnSigmaProton_P_true_pion","TOF number of sigmas, filled with secondary pion ;P [GeV/c];TOF number of sigmas Proton", 160,0.,8.,100,-5.,5.);
    TH2D* hTOFnSigmaProton_P_sec_kaon = new TH2D("TOFnSigmaProton_P_true_kaon","TOF number of sigmas, filled with secondary kaon ;P [GeV/c];TOF number of sigmas Proton", 160,0.,8.,100,-5.,5.);
    TH2D* hTOFnSigmaProton_P_sec_prot = new TH2D("TOFnSigmaProton_P_true_prot","TOF number of sigmas, filled with secondary prot ;P [GeV/c];TOF number of sigmas Proton", 160,0.,8.,100,-5.,5.);

    list_temp->AddAt(hPt_sec,     0);
    list_temp->AddAt(hITSnSigmaEle_P_sec, 1);
    list_temp->AddAt(hTPCnSigmaEle_P_sec, 2);
    list_temp->AddAt(hTOFnSigmaEle_P_sec, 3);
    list_temp->AddAt(hEta_sec, 4);
    list_temp->AddAt(hPhi_sec, 5);
    list_temp->AddAt(hEta_Phi_sec, 6);
    list_temp->AddAt(hTPCcrossedRowsOverFindable_sec, 7);
    list_temp->AddAt(hTPCcrossedRows_sec, 8);
    list_temp->AddAt(hTPCnCls_sec, 9);
    list_temp->AddAt(hITSnCls_sec, 10);
    list_temp->AddAt(hTPCchi2_sec, 11);
    list_temp->AddAt(hITSchi2_sec, 12);
    list_temp->AddAt(hTPCnClsS_sec, 13);
    list_temp->AddAt(hITSnClsS_sec, 14);
    list_temp->AddAt(hNclsSFracTPC_sec, 15);
    list_temp->AddAt(hNclsSFracITS_sec, 16);
    list_temp->AddAt(hTPCclsDiff_sec, 17);
    list_temp->AddAt(hTPCsignalN_sec, 18);
    list_temp->AddAt(hTPCcrossedRows_TPCnCls_sec, 19);
    list_temp->AddAt(hTPCcrossedRows_Pt_sec, 20);
    // list_temp->AddAt(hPDGCode_PDGCodeMother, 21);
    list_temp->AddAt(hPDGCode_sec, 21);
    list_temp->AddAt(hPDGCodeMother_sec, 22);
    list_temp->AddAt(hMCGenCode_sec, 23);
    list_temp->AddAt(hImpactParXY_sec, 24);
    list_temp->AddAt(hImpactParZ_sec, 25);
    list_temp->AddAt(hM_single_sec,26);

    list_nSigmaITS->AddAt(hITSnSigmaEle_P_sec_all ,0);
    list_nSigmaITS->AddAt(hITSnSigmaEle_P_sec_elec,1);
    list_nSigmaITS->AddAt(hITSnSigmaEle_P_sec_muon,2);
    list_nSigmaITS->AddAt(hITSnSigmaEle_P_sec_pion,3);
    list_nSigmaITS->AddAt(hITSnSigmaEle_P_sec_kaon,4);
    list_nSigmaITS->AddAt(hITSnSigmaEle_P_sec_prot,5);
    list_nSigmaITS->AddAt(hITSnSigmaMuon_P_sec_all ,6);
    list_nSigmaITS->AddAt(hITSnSigmaMuon_P_sec_elec,7);
    list_nSigmaITS->AddAt(hITSnSigmaMuon_P_sec_muon,8);
    list_nSigmaITS->AddAt(hITSnSigmaMuon_P_sec_pion,9);
    list_nSigmaITS->AddAt(hITSnSigmaMuon_P_sec_kaon,10);
    list_nSigmaITS->AddAt(hITSnSigmaMuon_P_sec_prot,11);
    list_nSigmaITS->AddAt(hITSnSigmaPion_P_sec_all ,12);
    list_nSigmaITS->AddAt(hITSnSigmaPion_P_sec_elec,13);
    list_nSigmaITS->AddAt(hITSnSigmaPion_P_sec_muon,14);
    list_nSigmaITS->AddAt(hITSnSigmaPion_P_sec_pion,15);
    list_nSigmaITS->AddAt(hITSnSigmaPion_P_sec_kaon,16);
    list_nSigmaITS->AddAt(hITSnSigmaPion_P_sec_prot,17);
    list_nSigmaITS->AddAt(hITSnSigmaKaon_P_sec_all ,18);
    list_nSigmaITS->AddAt(hITSnSigmaKaon_P_sec_elec,19);
    list_nSigmaITS->AddAt(hITSnSigmaKaon_P_sec_muon,20);
    list_nSigmaITS->AddAt(hITSnSigmaKaon_P_sec_pion,21);
    list_nSigmaITS->AddAt(hITSnSigmaKaon_P_sec_kaon,22);
    list_nSigmaITS->AddAt(hITSnSigmaKaon_P_sec_prot,23);
    list_nSigmaITS->AddAt(hITSnSigmaProton_P_sec_all ,24);
    list_nSigmaITS->AddAt(hITSnSigmaProton_P_sec_elec,25);
    list_nSigmaITS->AddAt(hITSnSigmaProton_P_sec_muon,26);
    list_nSigmaITS->AddAt(hITSnSigmaProton_P_sec_pion,27);
    list_nSigmaITS->AddAt(hITSnSigmaProton_P_sec_kaon,28);
    list_nSigmaITS->AddAt(hITSnSigmaProton_P_sec_prot,29);

    list_nSigmaTPC->AddAt(hTPCnSigmaEle_P_sec_all ,0);
    list_nSigmaTPC->AddAt(hTPCnSigmaEle_P_sec_elec,1);
    list_nSigmaTPC->AddAt(hTPCnSigmaEle_P_sec_muon,2);
    list_nSigmaTPC->AddAt(hTPCnSigmaEle_P_sec_pion,3);
    list_nSigmaTPC->AddAt(hTPCnSigmaEle_P_sec_kaon,4);
    list_nSigmaTPC->AddAt(hTPCnSigmaEle_P_sec_prot,5);
    list_nSigmaTPC->AddAt(hTPCnSigmaMuon_P_sec_all ,6);
    list_nSigmaTPC->AddAt(hTPCnSigmaMuon_P_sec_elec,7);
    list_nSigmaTPC->AddAt(hTPCnSigmaMuon_P_sec_muon,8);
    list_nSigmaTPC->AddAt(hTPCnSigmaMuon_P_sec_pion,9);
    list_nSigmaTPC->AddAt(hTPCnSigmaMuon_P_sec_kaon,10);
    list_nSigmaTPC->AddAt(hTPCnSigmaMuon_P_sec_prot,11);
    list_nSigmaTPC->AddAt(hTPCnSigmaPion_P_sec_all ,12);
    list_nSigmaTPC->AddAt(hTPCnSigmaPion_P_sec_elec,13);
    list_nSigmaTPC->AddAt(hTPCnSigmaPion_P_sec_muon,14);
    list_nSigmaTPC->AddAt(hTPCnSigmaPion_P_sec_pion,15);
    list_nSigmaTPC->AddAt(hTPCnSigmaPion_P_sec_kaon,16);
    list_nSigmaTPC->AddAt(hTPCnSigmaPion_P_sec_prot,17);
    list_nSigmaTPC->AddAt(hTPCnSigmaKaon_P_sec_all ,18);
    list_nSigmaTPC->AddAt(hTPCnSigmaKaon_P_sec_elec,19);
    list_nSigmaTPC->AddAt(hTPCnSigmaKaon_P_sec_muon,20);
    list_nSigmaTPC->AddAt(hTPCnSigmaKaon_P_sec_pion,21);
    list_nSigmaTPC->AddAt(hTPCnSigmaKaon_P_sec_kaon,22);
    list_nSigmaTPC->AddAt(hTPCnSigmaKaon_P_sec_prot,23);
    list_nSigmaTPC->AddAt(hTPCnSigmaProton_P_sec_all ,24);
    list_nSigmaTPC->AddAt(hTPCnSigmaProton_P_sec_elec,25);
    list_nSigmaTPC->AddAt(hTPCnSigmaProton_P_sec_muon,26);
    list_nSigmaTPC->AddAt(hTPCnSigmaProton_P_sec_pion,27);
    list_nSigmaTPC->AddAt(hTPCnSigmaProton_P_sec_kaon,28);
    list_nSigmaTPC->AddAt(hTPCnSigmaProton_P_sec_prot,29);

    list_nSigmaTOF->AddAt(hTOFnSigmaEle_P_sec_elec,0);
    list_nSigmaTOF->AddAt(hTOFnSigmaEle_P_sec_all ,1);
    list_nSigmaTOF->AddAt(hTOFnSigmaEle_P_sec_muon,2);
    list_nSigmaTOF->AddAt(hTOFnSigmaEle_P_sec_pion,3);
    list_nSigmaTOF->AddAt(hTOFnSigmaEle_P_sec_kaon,4);
    list_nSigmaTOF->AddAt(hTOFnSigmaEle_P_sec_prot,5);
    list_nSigmaTOF->AddAt(hTOFnSigmaMuon_P_sec_elec,6);
    list_nSigmaTOF->AddAt(hTOFnSigmaMuon_P_sec_all ,7);
    list_nSigmaTOF->AddAt(hTOFnSigmaMuon_P_sec_muon,8);
    list_nSigmaTOF->AddAt(hTOFnSigmaMuon_P_sec_pion,9);
    list_nSigmaTOF->AddAt(hTOFnSigmaMuon_P_sec_kaon,10);
    list_nSigmaTOF->AddAt(hTOFnSigmaMuon_P_sec_prot,11);
    list_nSigmaTOF->AddAt(hTOFnSigmaPion_P_sec_elec,12);
    list_nSigmaTOF->AddAt(hTOFnSigmaPion_P_sec_all ,13);
    list_nSigmaTOF->AddAt(hTOFnSigmaPion_P_sec_muon,14);
    list_nSigmaTOF->AddAt(hTOFnSigmaPion_P_sec_pion,15);
    list_nSigmaTOF->AddAt(hTOFnSigmaPion_P_sec_kaon,16);
    list_nSigmaTOF->AddAt(hTOFnSigmaPion_P_sec_prot,17);
    list_nSigmaTOF->AddAt(hTOFnSigmaKaon_P_sec_elec,18);
    list_nSigmaTOF->AddAt(hTOFnSigmaKaon_P_sec_all ,19);
    list_nSigmaTOF->AddAt(hTOFnSigmaKaon_P_sec_muon,20);
    list_nSigmaTOF->AddAt(hTOFnSigmaKaon_P_sec_pion,21);
    list_nSigmaTOF->AddAt(hTOFnSigmaKaon_P_sec_kaon,22);
    list_nSigmaTOF->AddAt(hTOFnSigmaKaon_P_sec_prot,23);
    list_nSigmaTOF->AddAt(hTOFnSigmaProton_P_sec_elec,24);
    list_nSigmaTOF->AddAt(hTOFnSigmaProton_P_sec_all ,25);
    list_nSigmaTOF->AddAt(hTOFnSigmaProton_P_sec_muon,26);
    list_nSigmaTOF->AddAt(hTOFnSigmaProton_P_sec_pion,27);
    list_nSigmaTOF->AddAt(hTOFnSigmaProton_P_sec_kaon,28);
    list_nSigmaTOF->AddAt(hTOFnSigmaProton_P_sec_prot,29);

    list2->Add(list_temp);
    list_temp->AddAt(list_nSigmaITS,27);
    list_temp->AddAt(list_nSigmaTPC,28);
    list_temp->AddAt(list_nSigmaTOF,29);
    }
  fTrackCutListVecSec.push_back(list2);
  fOutputListSupportHistos->Add(list2);
  }

  for (unsigned int list_i = 0; list_i < fPairCuts_secondary_standard.size(); ++list_i){
    TList* listPairCuts = new TList();
    listPairCuts->SetName(Form("SecondaryPairCuts: %s",fPairCuts_secondary_standard.at(list_i)->GetName()));
    listPairCuts->SetOwner();
    for (unsigned int i = 0; i < fSecondaryPairMCSignal.size(); ++i){
      TList* list_mcSig = new TList();
      list_mcSig->SetName(Form("SecondaryPairCuts: %s", fSecondaryPairMCSignal.at(i).GetName()));
      list_mcSig->SetOwner();

    // pair variable histos
    TH1D* hCosPointingAngle_sec = new TH1D("CosPointingAngle","Cosine of the pointing angle;#tracks",200, 0.8, 1.); // kCosPointingAngle
    TH1D* hChi2NDF_sec = new TH1D("Chi2NDF","Chi^2/NDF ;Chi^2/NDF;#tracks",2000, 0., 12.); // kChi2NDF (NDF = number degree of freedom)
    TH1D* hLegDist_sec = new TH1D("LegDist","Leg Distance ;Leg Distance;#tracks",1000, 0., 4.); // kLegDist
    TH1D* hR_sec = new TH1D("R","Distance to the origin ;r;#tracks",1500, 0., 150.); // kR
    TH1D* hPsiPair_sec = new TH1D("PsiPair","phi in mother's rest frame in Collins-Soper picture ;PsiPair;#tracks",200, 0., 2.); // kPsiPair
    TH1D* hM_pair_sec = new TH1D("PairMass","pair Mass ;mass ;#tracks",200, 0., 1.); // kM
    TH1D* hArmPt_sec = new TH1D("ArmPt","Armenteros-Podolanski pt ;Pt [GeV];#tracks",160,0.,1.); // kArmPt
    TH1D* hArmAlpha_sec = new TH1D("ArmAlpha","Armenteros-Podolanski alpha ;alpha;#tracks",200, -3., 3.); // kArmAlpha
    TH1D* hPairPt_sec = new TH1D("PairPt","Pt of pair ; Pt [GeV];#tracks",640, 0., 8.); // kArmAlpha
    TH2D* hArmAlpha_ArmPt_sec = new TH2D("Armenteros Podolanski Plot","Armenteros Podolanski Plot, secondary electrons ; ArmAlpha; ArmPt", 100, -1, 1, 320, 0, 0.1);//,AliDielectronVarManager::kEta,AliDielectronVarManager::kPhi);
    // 4 pair variable histos
    // TH1D* hFourPairOpeningAngle = new TH1D("Angle between primary and secondary pair", "Angle between primary and secondary pair ;rad ; #counts", 1000, 0. , 6.4);

    // pair support histos
    list_mcSig->AddAt(hCosPointingAngle_sec,0);
    list_mcSig->AddAt(hChi2NDF_sec,1);
    list_mcSig->AddAt(hLegDist_sec,2);
    list_mcSig->AddAt(hR_sec,3);
    list_mcSig->AddAt(hPsiPair_sec,4);
    list_mcSig->AddAt(hM_pair_sec,5);
    list_mcSig->AddAt(hArmPt_sec,6);
    list_mcSig->AddAt(hArmAlpha_sec,7);
    list_mcSig->AddAt(hPairPt_sec,8);
    list_mcSig->AddAt(hArmAlpha_ArmPt_sec,9);
    // 4 pair variable histos
    // list_mcSig->AddAt(hFourPairOpeningAngle,37);
    listPairCuts->Add(list_mcSig);
    }
  fPairCutListVecSec.push_back(listPairCuts);
  fOutputListSupportHistos->Add(listPairCuts);
  }


  for (unsigned int list_i = 0; list_i < fPairCuts_primary.size(); ++list_i){
    TList* listFourPair = new TList();
    listFourPair->SetName(Form("FourPairCuts: %s",fPairCuts_primary.at(list_i)->GetName()));
    listFourPair->SetOwner();
    for (unsigned int i = 0; i < fFourPairMCSignal.size(); i+=2){
      TList* list_temp = new TList();
      list_temp->SetName(Form("FourPairCuts: %s", fFourPairMCSignal.at(i).GetName()));
      list_temp->SetOwner();

    // 4 pair variable histos
    TH1D* hFourPairOpeningAngle = new TH1D("Angle between primary and secondary pair", "Angle between primary and secondary pair ;rad ; #counts", 1000, 0. , 6.4);

    list_temp->AddAt(hFourPairOpeningAngle,0);
    listFourPair->Add(list_temp);
    }
  fFourPairCutListVec.push_back(listFourPair);
  fOutputListSupportHistos->Add(listFourPair);
  }

}

// ############################################################################
// ############################################################################
TLorentzVector AliAnalysisTaskEtaReconstruction::ApplyResolution(double pt, double eta, double phi, short ch) {
  // from Theos LightFlavorGenerator, modified
                                                                                if(fdebug) gRandom->SetSeed(122);  // Seed in order to proove consistance of code durig changes. Disables the randomness of smearing. Comment out for normal calculation

  TLorentzVector resvec;
  // resvec.SetPtEtaPhiM(pt, eta, phi, AliPID::ParticleMass(AliPID::kElectron));

    // smear pt
    Int_t ptbin     = ((TH2D *)(fArrResoPt->At(0)))->GetXaxis()->FindBin(pt);
    Int_t ptbin_max = ((TH2D *)(fArrResoPt->At(0)))->GetXaxis()->GetNbins();
    // make sure that no underflow or overflow bins are used
    if (ptbin < 1)
      ptbin = 1;
    else if (ptbin > ptbin_max)
      ptbin = ptbin_max;
    Double_t smearing = ((TH1D *)(fArrResoPt->At(ptbin)))->GetRandom() * pt;
    const Double_t sPt = pt - smearing;

    // smear eta
    ptbin     = ((TH2D *)(fArrResoEta->At(0)))->GetXaxis()->FindBin(pt);
    ptbin_max = ((TH2D *)(fArrResoEta->At(0)))->GetXaxis()->GetNbins();
    if (ptbin < 1)
      ptbin = 1;
    else if (ptbin > ptbin_max)
      ptbin = ptbin_max;
    smearing = ((TH1D *)(fArrResoEta->At(ptbin)))->GetRandom();
    const Double_t sEta = eta - smearing;

    // smear phi
    ptbin     = ((TH2D *)(fArrResoPhi_Pos->At(0)))->GetXaxis()->FindBin(pt);
    ptbin_max = ((TH2D *)(fArrResoPhi_Pos->At(0)))->GetXaxis()->GetNbins();
    if (ptbin < 1)
      ptbin = 1;
    if (ptbin > ptbin_max)
      ptbin = ptbin_max;
    if (ch > 0) {
      smearing = ((TH1D *)(fArrResoPhi_Pos->At(ptbin)))->GetRandom();
    } else if (ch < 0) {
      smearing = ((TH1D *)(fArrResoPhi_Neg->At(ptbin)))->GetRandom();
    }
    const Double_t sPhi = phi - smearing;

    // printf(" Original Pt = %f Phi %f Eta %f -> final pt = %f Phi %f Eta %f\n",pt,phi,eta,sPt,sPhi,sEta);

    resvec.SetPtEtaPhiM(sPt, sEta, sPhi, AliPID::ParticleMass(AliPID::kElectron));
  return resvec;
}

// ############################################################################
// ############################################################################
Double_t AliAnalysisTaskEtaReconstruction::GetSmearing(TObjArray *arr, Double_t x)
{
  TH1D *hisSlice(0x0);
  TH2D *hDeltaXvsXgen = static_cast<TH2D*> (arr->At(0));
  Int_t histIndex = TMath::Min( hDeltaXvsXgen->GetXaxis()->FindBin(x), arr->GetLast() );
  if (histIndex<1) histIndex=1;
  hisSlice = static_cast<TH1D*> (arr->At(histIndex));
  // get smear parameter via random selection from the x slices retreived from the deltax plot
  Double_t smearing(0.);
  if(hisSlice) smearing = hisSlice->GetRandom();
  delete hisSlice;
  delete hDeltaXvsXgen;
  return smearing;
}

bool AliAnalysisTaskEtaReconstruction::CheckGenerator(int trackID, std::vector<unsigned int> vecHashes){
  if (vecHashes.size() == 0) return true;

  TString genname;
  Bool_t hasGenerator = fMC->GetCocktailGenerator(TMath::Abs(trackID), genname);
  // std::cout << genname << std::endl;
  if(!hasGenerator) {
    Printf("no cocktail header list was found for this track");
    return false;
  }
  else{
    for (unsigned int i = 0; i < vecHashes.size(); ++i){
      // std::cout << genname.Hash() << " " << vecHashes[i] << std::endl;
      if (genname.Hash() == vecHashes[i]) return true;

    }
    return false;
  }
  return false; // should not happen
}

double AliAnalysisTaskEtaReconstruction::GetWeight(Particle part1, Particle part2, double motherpt){
  int pdgMother = 0;
  double weight = 0;
  int bin = -1;

  pdgMother = fMC->GetTrack(part2.GetMotherID())->PdgCode();
  if      (pdgMother == 111) {
    if (fPtPion) {
      bin = fPtPion->FindBin(motherpt);
      weight = fPtPion->GetBinContent(bin);
    }
  }
  else if (pdgMother == 221) {
    if (fPtEta) {
      bin = fPtEta->FindBin(motherpt);
      weight = fPtEta->GetBinContent(bin);
    }
  }
  else if (pdgMother == 331) {
    if (fPtEtaPrime) {
      bin = fPtEtaPrime->FindBin(motherpt);
      weight = fPtEtaPrime->GetBinContent(bin);
    }
  }
  else if (pdgMother == 113) {
    if (fPtRho) {
      bin = fPtRho->FindBin(motherpt);
      weight = fPtRho->GetBinContent(bin);
    }
  }
  else if (pdgMother == 223) {
    if (fPtOmega) {
      bin = fPtOmega->FindBin(motherpt);
      weight = fPtOmega->GetBinContent(bin);
    }
  }
  else if (pdgMother == 333) {
    if (fPtPhi) {
      bin = fPtPhi->FindBin(motherpt);
      weight = fPtPhi->GetBinContent(bin);
    }
  }
  else if (pdgMother == 443) {
    if (fPtJPsi) {
      bin = fPtJPsi->FindBin(motherpt);
      weight = fPtJPsi->GetBinContent(bin);
    }
  }

  // std::cout << "weight from " << pdgMother << " for pt = " << motherpt << "  weight: " << weight << "  bin: " << bin << std::endl;

  return weight;
}

//______________________________________________
void AliAnalysisTaskEtaReconstruction::SetCentroidCorrFunction(Detector det, TObject *fun, UInt_t varx, UInt_t vary, UInt_t varz)
{
  if (fun == 0x0) {
    std::cout << "No correction function chosen" << std::endl;
    return;
  }

  std::string detector = "";
  if      (det == kITS) detector = "ITS";
  else if (det == kTPC) detector = "TPC";
  else if (det == kTOF) detector = "TOF";
  else {
    std::cout << "ERROR: No matching detector, must be kITS, kTPC or kTOF. return;" << std::endl;
    return;
  }
  std::cout << "Do centroid correction with detector " << detector << std::endl;

  TH1* correctionMap = 0x0;
  if      (det == kITS) {
    fPostPIDCntrdCorrITS = dynamic_cast<TH1*>(fun);
    correctionMap = fPostPIDCntrdCorrITS;
  }
  else if (det == kTPC) {
    fPostPIDCntrdCorrTPC = dynamic_cast<TH1*>(fun);
    correctionMap = fPostPIDCntrdCorrTPC;
  }
  else if (det == kTOF) {
    fPostPIDCntrdCorrTOF = dynamic_cast<TH1*>(fun);
    correctionMap = fPostPIDCntrdCorrTOF;
  }

  UInt_t valType[20] = {0};
  valType[0]=varx;     valType[1]=vary;     valType[2]=varz;
  TString key = Form("cntrdTPC%d%d%d",varx,vary,varz);
  // StoreVariables() sets the uniqueIDs of the axes to match with the VarManager enum and sets axis titles accordingly.
  // Needed so that AliDielectronPID can extract the needed variables into its local TBits 'fUsedVars'.
  // clone temporary histogram, otherwise it will not be streamed to file!
  if (fun->InheritsFrom(TH1::Class())) {
    AliDielectronHistos::StoreVariables(fun, valType);
    correctionMap = (TH1*)fun->Clone(key.Data());
  }
  else {
    AliWarning(Form("WARNING: PID correction object has invalid type: %s!", fun->IsA()->GetName())); return;
  }
  if(correctionMap) {
    // check for corrections and add their variables to the fill map
    std::cout << "POST " << detector << "PID CORRECTION added for centroids:  " << std::endl;
    switch(correctionMap->GetDimension()) {
      case 3: printf(" %s, ",correctionMap->GetZaxis()->GetName());
      case 2: printf(" %s, ",correctionMap->GetYaxis()->GetName());
      case 1: printf(" %s ",correctionMap->GetXaxis()->GetName());
    }
    printf("\n");
    // fUsedVars->SetBitNumber(varx, kTRUE); // probably not needed within efficiency task...
    // fUsedVars->SetBitNumber(vary, kTRUE);
    // fUsedVars->SetBitNumber(varz, kTRUE);
  }
}
//______________________________________________
void AliAnalysisTaskEtaReconstruction::SetWidthCorrFunction(Detector det, TObject *fun, UInt_t varx, UInt_t vary, UInt_t varz)
{
  if (fun == 0x0) {
    std::cout << "No correction function chosen" << std::endl;
    return;
  }

  std::string detector = "";
  if      (det == kITS) detector = "ITS";
  else if (det == kTPC) detector = "TPC";
  else if (det == kTOF) detector = "TOF";
  else {
    std::cout << "ERROR: No matching detector, must be kITS, kTPC or kTOF. return;" << std::endl;
    return;
  }
  std::cout << "Do width correction with detector " << detector << std::endl;

  TH1* correctionMap = 0x0;
  if      (det == kITS) {
    fPostPIDWdthCorrITS = dynamic_cast<TH1*>(fun);
    correctionMap = fPostPIDWdthCorrITS;
  }
  else if (det == kTPC) {
    fPostPIDWdthCorrTPC = dynamic_cast<TH1*>(fun);
    correctionMap = fPostPIDWdthCorrTPC;
  }
  else if (det == kTOF) {
    fPostPIDWdthCorrTOF = dynamic_cast<TH1*>(fun);
    correctionMap = fPostPIDWdthCorrTOF;
  }


  UInt_t valType[20] = {0};
  valType[0]=varx;     valType[1]=vary;     valType[2]=varz;
  TString key = Form("wdthTPC%d%d%d",varx,vary,varz);
  // StoreVariables() sets the uniqueIDs of the axes to match with the VarManager enum and sets axis titles accordingly.
  // Needed so that AliDielectronPID can extract the needed variables into its local TBits 'fUsedVars'.
  // clone temporary histogram, otherwise it will not be streamed to file!
  if (fun->InheritsFrom(TH1::Class())) {
    AliDielectronHistos::StoreVariables(fun, valType);
    correctionMap = (TH1*)fun->Clone(key.Data());
  }
  else {
    AliWarning(Form("WARNING: PID correction object has invalid type: %s!", fun->IsA()->GetName())); return;
  }
  if(correctionMap)  {
    // check for corrections and add their variables to the fill map
    std::cout << "POST " << detector << "PID CORRECTION added for widths:  " << std::endl;
    switch(correctionMap->GetDimension()) {
      case 3: printf(" %s, ",correctionMap->GetZaxis()->GetName());
      case 2: printf(" %s, ",correctionMap->GetYaxis()->GetName());
      case 1: printf(" %s ",correctionMap->GetXaxis()->GetName());
    }
    printf("\n");
    // fUsedVars->SetBitNumber(varx, kTRUE); // probably not needed within efficiency task...
    // fUsedVars->SetBitNumber(vary, kTRUE);
    // fUsedVars->SetBitNumber(varz, kTRUE);
  }
}


void AliAnalysisTaskEtaReconstruction::ApplyStandardCutsAndFillHists (std::vector<TwoPair>* fPairVec, std::vector<AliAnalysisFilter*> fCutsetting, Bool_t TrackCuts, Bool_t PairPrimary, double centralityWeight) {
  // If TrackCuts   = true :   Aplly standard TrackCuts to primary particles (there are no track cuts on secondaries yet)
  // If TrackCuts   = false:   Aplly standard PairCuts to either primary pairs or secondary pairs
  // If PairPrimary = true :   fPairVec is the reconstructed primary pair vector
  // If PairPrimary = false:   fPairVec is the reconstructed secondary pair vector

  // ###############################
  // primary part
  // ###############################
  for (unsigned int iPair = 0; iPair < fPairVec->size(); iPair++){
    if (PairPrimary == kTRUE) {
      AliVParticle* negTrack = fEvent->GetTrack(fPairVec->at(iPair).GetFirstDaughter());
      AliVParticle* posTrack = fEvent->GetTrack(fPairVec->at(iPair).GetSecondDaughter());

      int neglabel = negTrack->GetLabel();
      int poslabel = posTrack->GetLabel();
      int absneglabel = TMath::Abs(neglabel);
      int absposlabel = TMath::Abs(poslabel);

      Particle negPart  = CreateParticle(negTrack);
      Particle posPart  = CreateParticle(posTrack);

      // ###############################
      // Track cuts on primary pair
      // ###############################
      if(TrackCuts == kTRUE){
                                                                                // if(fdebug) std::cout << __LINE__ << " Cout Debug: Apply standard priamry track cuts. " << std::endl;
        //  ---------- RECONSTRUCTED Tracks  ---------- //
        // Check if pos and neg Track are passing cut selections
        std::vector<bool> selected_posTrack(fCutsetting.size(), kFALSE); // vector which stores if track is accepted by [i]-th selection cut
        for (UInt_t iCut=0; iCut<fCutsetting.size(); ++iCut){ // loop over all specified cutInstances
          UInt_t selectedMask_posTrack=( 1 << fCutsetting.at(iCut)->GetCuts()->GetEntries())-1;
          // cutting logic taken from AliDielectron::FillTrackArrays()          FillPairArrays()
          // apply track cuts
          UInt_t cutMask_posTrack = fCutsetting.at(iCut)->IsSelected(posTrack);
          if (cutMask_posTrack == selectedMask_posTrack) {selected_posTrack[iCut] = kTRUE; /*std::cout << "reconstructed TRUE" << std::endl;*/}
        }
        std::vector<bool> selected_negTrack(fCutsetting.size(), kFALSE); // vector which stores if track is accepted by [i]-th selection cut
        for (UInt_t iCut=0; iCut<fCutsetting.size(); ++iCut){ // loop over all specified cutInstances
          UInt_t selectedMask_negTrack=( 1 << fCutsetting.at(iCut)->GetCuts()->GetEntries())-1;
          // cutting logic taken from AliDielectron::FillTrackArrays()          FillPairArrays()
          // apply track cuts
          UInt_t cutMask_negTrack = fCutsetting.at(iCut)->IsSelected(negTrack);
          if (cutMask_negTrack == selectedMask_negTrack) {selected_negTrack[iCut] = kTRUE; /*std::cout << "reconstructed TRUE" << std::endl;*/}
        }
        // ##########################################################
        // check if at least one is selected by cuts otherwise remove this pair from vector
        if ((CheckIfOneIsTrue(selected_posTrack) == kFALSE) || (CheckIfOneIsTrue(selected_negTrack) == kFALSE)) {
          fPairVec->erase(fPairVec->begin()+iPair);
          iPair--;
          continue;
        }
        std::vector<Bool_t> mcSignal_acc_posPart(fSinglePrimaryLegMCSignal.size(), kFALSE); // vector which stores if track is accepted by [i]-th mcsignal
        std::vector<Bool_t> mcSignal_acc_negPart(fSinglePrimaryLegMCSignal.size(), kFALSE); // vector which stores if track is accepted by [i]-th mcsignal
        CheckSinglePrimaryLegMCsignals(mcSignal_acc_posPart, absposlabel);
        CheckSinglePrimaryLegMCsignals(mcSignal_acc_negPart, absneglabel);


        posPart.isMCSignal_primary   = mcSignal_acc_posPart;
        negPart.isMCSignal_primary   = mcSignal_acc_negPart;
        posPart.isReconstructed_primary = selected_posTrack;
        negPart.isReconstructed_primary = selected_negTrack;

        fPairVec->at(iPair).SetDauthersAreReconstructed(selected_negTrack,selected_posTrack);

        AliVParticle* mcPosPart = fMC->GetTrack(absposlabel);
        AliVParticle* mcNegPart = fMC->GetTrack(absneglabel);

        // check if MCSignal and isReconstructed for positive primary track
        if (posPart.fCharge > 0) {
          for (unsigned int i = 0; i < posPart.isMCSignal_primary.size(); ++i){
            for (unsigned int j = 0; j < posPart.isReconstructed_primary.size(); ++j){
              if (posPart.isMCSignal_primary[i] == kTRUE) {
                if (posPart.isReconstructed_primary[j] == kTRUE){
                  dynamic_cast<TH3D*>(fHistRecPrimaryPosPart.at(j * posPart.isMCSignal_primary.size() + i))->Fill(posPart.fPt, posPart.fEta, posPart.fPhi, centralityWeight);
                  FillTrackHistograms_Primary(posTrack, mcPosPart, j, i); // Fill primary track support histograms
                }// is selected by cutsetting
              } // is selected by MC signal
            } // end of loop over all cutsettings
          } // end of loop over all MCsignals
        } // end if particle positive

        // check if MCSignal and isReconstructed for negative primaryTrack
        if (negPart.fCharge < 0) {
          for (unsigned int i = 0; i < negPart.isMCSignal_primary.size(); ++i){
            for (unsigned int j = 0; j < negPart.isReconstructed_primary.size(); ++j){
              if (negPart.isMCSignal_primary[i] == kTRUE) {
                if (negPart.isReconstructed_primary[j] == kTRUE){
                  dynamic_cast<TH3D*>(fHistRecPrimaryNegPart.at(j * negPart.isMCSignal_primary.size() + i))->Fill(negPart.fPt, negPart.fEta, negPart.fPhi, centralityWeight);
                  FillTrackHistograms_Primary(negTrack, mcNegPart, j, i); // Fill primary track support histograms
                }// is selected by cutsetting
              } // is selected by MC signal
            } // end of loop over all cutsettings
          } // end of loop over all MCsignals
        } // end if particle negative
      } // end if TrackCuts true



      // ###############################
      // Pair cuts on primary pair
      // ###############################
      else if (TrackCuts == kFALSE) {
                                                                                // if(fdebug) std::cout << __LINE__ << " Cout Debug: Apply standard primary pair cuts. " << std::endl;
        AliDielectronPair *pair = new AliDielectronPair;
        pair->SetKFUsage(false);
        pair->SetTracks(static_cast<AliVTrack*>(negTrack), 11, static_cast<AliVTrack*>(posTrack), -11);

        //  ---------- RECONSTRUCTED Pairs: Apply pair cuts on PAIR  ---------- //
        // Check if pair is passing cut selections
        std::vector<bool> selected(fCutsetting.size(), kFALSE); // vector which stores if pair is accepted by [i]-th selection cut
        for (UInt_t iCut=0; iCut<fCutsetting.size(); ++iCut){ // loop over all specified cutInstances
          if (fPairVec->at(iPair).fFirstPartIsReconstructed[iCut] == kFALSE || fPairVec->at(iPair).fSecondPartIsReconstructed[iCut] == kFALSE) continue;  //check if both daughters are reconstructed
          UInt_t selectedMask=( 1 << fCutsetting.at(iCut)->GetCuts()->GetEntries())-1;
          // cutting logic taken from AliDielectron::FillTrackArrays()          FillPairArrays()
          // apply pair cuts
          UInt_t cutMask = fCutsetting.at(iCut)->IsSelected(pair);
          if (cutMask == selectedMask) {selected[iCut] = kTRUE; /*std::cout << "reconstructed TRUE" << std::endl;*/}
        }
        // ##########################################################
        // check if at least one is selected by cuts otherwise remove this pair from vector
        if (CheckIfOneIsTrue(selected) == kFALSE) {
          fPairVec->erase(fPairVec->begin()+iPair);
          iPair--;
          continue;
        }

        fPairVec->at(iPair).isReconstructed_primary = selected;


        std::vector<Bool_t> mcTwoSignal_acc(fPrimaryPairMCSignal.size(), kFALSE); // vector which stores if pair is accepted by [i]-th mcsignal
        for (unsigned int i = 0; i < fPrimaryPairMCSignal.size(); ++i){
          mcTwoSignal_acc[i] = AliDielectronMC::Instance()->IsMCTruth(pair, &(fPrimaryPairMCSignal[i]));
        }

        // check if at least one mc signal is true
        if (CheckIfOneIsTrue(mcTwoSignal_acc) == kFALSE) {
          fPairVec->erase(fPairVec->begin()+iPair);
          iPair--;
          continue;
        }

        double weight  = 1.;
        double mass    = fPairVec->at(iPair).fMass;
        double pairpt  = fPairVec->at(iPair).fPt;


        // Filling reconstructed primary particle histograms according to MCSignals
        if (PairPrimary == kTRUE) {
          for (unsigned int i = 0; i < fPairVec->at(iPair).isReconstructed_primary.size(); ++i){
            if (fPairVec->at(iPair).isReconstructed_primary[i] == kTRUE){
              for (unsigned int j = 0; j < mcTwoSignal_acc.size(); ++j){
                if (mcTwoSignal_acc[j] == kTRUE){
                  fHistRecPrimaryPair.at(  i *   mcTwoSignal_acc.size() + j)->Fill(mass, pairpt, weight * centralityWeight);
                  // FillPairHistograms_Primary(pair, i, j);                 // Fill pair primary histograms
                }// is selected by cutsetting
              } // end of loop over all cutsettings
            } // is selected by MCSignal
          } // end of loop over all MCsignals
        } // end of if PairPrimary
        delete pair;
      } // end if TrackCuts false
    } // end if Pair primary




    // ###############################
    // secondary part
    // ###############################
    if (PairPrimary == kFALSE) {
      // Get V0 from Event, selected trough V0ID stored in TwoPair and is set in reconstructed secondary two pairing
      Int_t fV0ID = fPairVec->at(iPair).GetV0ID();
      AliAODv0* fV0 = ((AliAODEvent*)fEvent)->GetV0(fV0ID);
      AliAODTrack *TrackD1 = (AliAODTrack *) (fV0->GetSecondaryVtx()->GetDaughter(0));
      AliAODTrack *TrackD2 = (AliAODTrack *) (fV0->GetSecondaryVtx()->GetDaughter(1));

      int labelD1 = TrackD1->GetLabel();
      int labelD2 = TrackD2->GetLabel();
      int absD1label = TMath::Abs(labelD1);
      int absD2label = TMath::Abs(labelD2);

      Particle PartD1  = CreateParticle(TrackD1);
      Particle PartD2  = CreateParticle(TrackD2);

      // ###############################
      // Track cuts on secondary tracks
      // ###############################
      if (TrackCuts == kTRUE) {
                                                                                // if(fdebug) std::cout << __LINE__ << " Cout Debug: Apply standard secondary track cuts. " << std::endl;
        // Check if Daughter1 and Daughter2 Tracks are passing cut selections
        std::vector<bool> selected_TrackD1(fCutsetting.size(), kFALSE); // vector which stores if track is accepted by [i]-th selection cut
        for (UInt_t iCut=0; iCut<fCutsetting.size(); ++iCut){ // loop over all specified cutInstances
          UInt_t selectedMask_TrackD1=( 1 << fCutsetting.at(iCut)->GetCuts()->GetEntries())-1;
          // cutting logic taken from AliDielectron::FillTrackArrays()          FillPairArrays()
          // apply track cuts
          UInt_t cutMask_TrackD1 = fCutsetting.at(iCut)->IsSelected(TrackD1);
          if (cutMask_TrackD1 == selectedMask_TrackD1) {selected_TrackD1[iCut] = kTRUE; /*std::cout << "reconstructed TRUE" << std::endl;*/}
        }
        std::vector<bool> selected_TrackD2(fCutsetting.size(), kFALSE); // vector which stores if track is accepted by [i]-th selection cut
        for (UInt_t iCut=0; iCut<fCutsetting.size(); ++iCut){ // loop over all specified cutInstances
          UInt_t selectedMask_TrackD2=( 1 << fCutsetting.at(iCut)->GetCuts()->GetEntries())-1;
          // cutting logic taken from AliDielectron::FillTrackArrays()          FillPairArrays()
          // apply track cuts
          UInt_t cutMask_TrackD2 = fCutsetting.at(iCut)->IsSelected(TrackD2);
          if (cutMask_TrackD2 == selectedMask_TrackD2) {selected_TrackD2[iCut] = kTRUE; /*std::cout << "reconstructed TRUE" << std::endl;*/}
        }
        // ##########################################################
        // check if at least one is selected by cuts otherwise remove this pair from vector
        if ((CheckIfOneIsTrue(selected_TrackD1) == kFALSE) || (CheckIfOneIsTrue(selected_TrackD2) == kFALSE)) {
          fPairVec->erase(fPairVec->begin()+iPair);
          iPair--;
          continue;
        }

        std::vector<Bool_t> mcSignal_acc_TrackD1(fSingleSecondaryLegMCSignal.size(), kFALSE); // vector which stores if track is accepted by [i]-th mcsignal
        std::vector<Bool_t> mcSignal_acc_TrackD2(fSingleSecondaryLegMCSignal.size(), kFALSE); // vector which stores if track is accepted by [i]-th mcsignal
        CheckSingleSecondaryLegMCsignals(mcSignal_acc_TrackD1, absD1label);
        CheckSingleSecondaryLegMCsignals(mcSignal_acc_TrackD2, absD2label);

        PartD1.isMCSignal_secondary   = mcSignal_acc_TrackD1;
        PartD2.isMCSignal_secondary   = mcSignal_acc_TrackD2;
        PartD1.isReconstructed_secondary = selected_TrackD1;
        PartD2.isReconstructed_secondary = selected_TrackD2;

        fPairVec->at(iPair).SetDauthersAreReconstructed(selected_TrackD1,selected_TrackD2);

        // AliVParticle* mcPartD1 = fMC->GetTrack(absD1label);
        // AliVParticle* mcPartD2 = fMC->GetTrack(absD2label);


        // check if MCSignal and isReconstructed for positive primary track
        for (unsigned int i = 0; i < PartD1.isMCSignal_secondary.size(); ++i){
          for (unsigned int j = 0; j < PartD1.isReconstructed_secondary.size(); ++j){
            if (PartD1.isMCSignal_secondary[i] == kTRUE) {
              if (PartD1.isReconstructed_secondary[j] == kTRUE){
                if (PartD1.fCharge > 0) {
                  dynamic_cast<TH3D*>(fHistRecSecondaryPosPart.at(j * PartD1.isMCSignal_secondary.size() + i))->Fill(PartD1.fPt, PartD1.fEta, PartD1.fPhi, centralityWeight);
                  // FillTrackHistograms_Secondary(TrackD1, mcPartD1, j, i); // Fill secondary track support histograms
                } // end if particle positive
                if (PartD1.fCharge < 0) {
                  dynamic_cast<TH3D*>(fHistRecSecondaryNegPart.at(j * PartD1.isMCSignal_secondary.size() + i))->Fill(PartD1.fPt, PartD1.fEta, PartD1.fPhi, centralityWeight);
                  // FillTrackHistograms_Secondary(TrackD1, mcPartD1, j, i); // Fill secondary track support histograms
                } // end if particle negative
              }// is selected by cutsetting
            } // is selected by MC signal
          } // end of loop over all cutsettings
        } // end of loop over all MCsignals
        // check if MCSignal and isReconstructed for negative primaryTrack
        for (unsigned int i = 0; i < PartD2.isMCSignal_secondary.size(); ++i){
          for (unsigned int j = 0; j < PartD2.isReconstructed_secondary.size(); ++j){
            if (PartD2.isMCSignal_secondary[i] == kTRUE) {
              if (PartD2.isReconstructed_secondary[j] == kTRUE){
                if (PartD2.fCharge > 0) {
                  dynamic_cast<TH3D*>(fHistRecSecondaryPosPart.at(j * PartD2.isMCSignal_secondary.size() + i))->Fill(PartD2.fPt, PartD2.fEta, PartD2.fPhi, centralityWeight);
                  // FillTrackHistograms_Secondary(TrackD2, mcPartD2, j, i); // Fill secondary track support histograms
                } // end if particle positive
                if (PartD2.fCharge < 0) {
                  dynamic_cast<TH3D*>(fHistRecSecondaryNegPart.at(j * PartD2.isMCSignal_secondary.size() + i))->Fill(PartD2.fPt, PartD2.fEta, PartD2.fPhi, centralityWeight);
                  // FillTrackHistograms_Secondary(TrackD2, mcPartD2, j, i); // Fill secondary track support histograms
                } // end if particle negative
              }// is selected by cutsetting
            } // is selected by MC signal
          } // end of loop over all cutsettings
        } // end of loop over all MCsignals
      } // end if TrackCuts true



      // ###############################
      // Pair cuts on secondary pair
      // ###############################
      else if(TrackCuts == kFALSE) {
                                                                                // if(fdebug) std::cout << __LINE__ << " Cout Debug: Apply standard secondary pair cuts. " << std::endl;
        AliDielectronPair *pair = new AliDielectronPair;
        pair->SetKFUsage(false);
        pair->SetTracks(static_cast<AliAODTrack*>(TrackD2), 11, static_cast<AliAODTrack*>(TrackD1), -11);

        //  ---------- RECONSTRUCTED Pairs with PAIR  ---------- //
        // Check if pair is passing cut selections
        std::vector<bool> selected(fCutsetting.size(), kFALSE); // vector which stores if pair is accepted by [i]-th selection cut
        for (UInt_t iCut=0; iCut<fCutsetting.size(); ++iCut){ // loop over all specified cutInstances
          if (fPairVec->at(iPair).fFirstPartIsReconstructed[iCut] == kFALSE || fPairVec->at(iPair).fSecondPartIsReconstructed[iCut] == kFALSE) continue;  //check if both daughters are reconstructed
          UInt_t selectedMask=( 1 << fCutsetting.at(iCut)->GetCuts()->GetEntries())-1;
          // cutting logic taken from AliDielectron::FillTrackArrays()          FillPairArrays()
          // apply pair cuts
          UInt_t cutMask = fCutsetting.at(iCut)->IsSelected(pair);
          if (cutMask == selectedMask) {selected[iCut] = kTRUE; /*std::cout << "reconstructed TRUE" << std::endl;*/}
        }
        // ##########################################################
        // check if at least one is selected by cuts otherwise remove this pair from vector
        if (CheckIfOneIsTrue(selected) == kFALSE) {
          fPairVec->erase(fPairVec->begin()+iPair);
          iPair--;
          continue;
        }

        fPairVec->at(iPair).isReconstructed_secondary = selected;

        std::vector<Bool_t> mcTwoSignal_acc(fSecondaryPairMCSignal.size(), kFALSE); // vector which stores if pair is accepted by [i]-th mcsignal
        for (unsigned int i = 0; i < fSecondaryPairMCSignal.size(); ++i){
          mcTwoSignal_acc[i] = AliDielectronMC::Instance()->IsMCTruth(pair, &(fSecondaryPairMCSignal[i]));
        }

        // check if at least one mc signal is true
        if (CheckIfOneIsTrue(mcTwoSignal_acc) == kFALSE) continue;

        double weight  = 1.;
        double mass    = fPairVec->at(iPair).fMass;
        double pairpt  = fPairVec->at(iPair).fPt;


        // track label, points back to MC track
        Int_t TrackD1Label = TMath::Abs(TrackD1->GetLabel());
        Int_t TrackD2Label = TMath::Abs(TrackD2->GetLabel());

        // Fill secondary support histos
        AliVParticle* mcPosPart = fMC->GetTrack(TrackD1Label);
        AliVParticle* mcNegPart = fMC->GetTrack(TrackD2Label);


        for (unsigned int i = 0; i < fPairVec->at(iPair).isReconstructed_secondary.size(); ++i){
          if (fPairVec->at(iPair).isReconstructed_secondary[i] == kTRUE){
            for (unsigned int j = 0; j < mcTwoSignal_acc.size(); ++j){
              if (mcTwoSignal_acc[j] == kTRUE){
                fHistRecSecondaryPair.at(  i *   mcTwoSignal_acc.size() + j)->Fill(mass, pairpt, weight * centralityWeight);
                FillTrackHistograms_Secondary(TrackD1, mcPosPart, i, j); // Fill track secondary support histograms
                FillTrackHistograms_Secondary(TrackD2, mcNegPart, i, j); // Fill track secondary support histograms
                FillPairHistograms_Secondary(pair, i, j);                 // Fill pair  secondary support histograms
              } // is selected by MCSignal
            } // end of loop over all MCsignals
          } // is selected by cutsetting
        } // end of loop over all cutsettings
        delete pair;
      } // end if TrackCuts false
    } // end if Pair secondary
  } // End of iPair loop
}


// ##########################################################
// ##########################################################
// ##########################################################
// Doing Two Pairing

// generated and generated smeared
void AliAnalysisTaskEtaReconstruction::DoGenAndGenSmearTwoPairing(std::vector<Particle>* vec_negParticle, std::vector<Particle>* vec_posParticle, Bool_t PartPrimary, Bool_t SmearedPair, double centralityWeight){
  for (unsigned int neg_i = 0; neg_i < vec_negParticle->size(); ++neg_i){
    for (unsigned int pos_i = 0; pos_i < vec_posParticle->size(); ++pos_i){
      AliVParticle* mcPart1 = fMC->GetTrack(vec_negParticle->at(neg_i).GetTrackID());
      AliVParticle* mcPart2 = fMC->GetTrack(vec_posParticle->at(pos_i).GetTrackID());


      // Check if electrons are from MCSignal Generator
      if (!vec_posParticle->at(pos_i).GetMCSignalPair() || !vec_negParticle->at(neg_i).GetMCSignalPair()) continue;
      // std::cout << "vec_posParticle->at(pos_i).GetMCSignalPair() = " << vec_posParticle->at(pos_i).GetMCSignalPair() << std::endl;
      // std::cout << "vec_negParticle->at(neg_i).GetMCSignalPair() = " << vec_negParticle->at(neg_i).GetMCSignalPair() << std::endl;
      // std::cout << "#########" << std::endl;
      // Apply MC signals
      std::vector<Bool_t> mcTwoSignal_acc_primary(fPrimaryPairMCSignal.size(), kFALSE); // vector which stores if track is accepted by [i]-th mcsignal
      std::vector<Bool_t> mcTwoSignal_acc_secondary(fSecondaryPairMCSignal.size(), kFALSE); // vector which stores if track is accepted by [i]-th mcsignal

      // Check if it according to mcsignals
      for (unsigned int i = 0; i < fPrimaryPairMCSignal.size(); ++i){
        mcTwoSignal_acc_primary[i] = AliDielectronMC::Instance()->IsMCTruth(mcPart1, mcPart2, &(fPrimaryPairMCSignal[i]));
      }
      for (unsigned int i = 0; i < fSecondaryPairMCSignal.size(); ++i){
        mcTwoSignal_acc_secondary[i] = AliDielectronMC::Instance()->IsMCTruth(mcPart1, mcPart2, &(fSecondaryPairMCSignal[i]));
      }

      // check if at least one mc signal is true
      if      (PartPrimary == kTRUE  && CheckIfOneIsTrue(mcTwoSignal_acc_primary)   == kFALSE) continue;
      else if (PartPrimary == kFALSE && CheckIfOneIsTrue(mcTwoSignal_acc_secondary) == kFALSE) continue;


      // Construct pair variables from LorentzVectors
      TLorentzVector Lvec1;
      TLorentzVector Lvec2;
      if (SmearedPair == kFALSE) {
        Lvec1.SetPtEtaPhiM(mcPart1->Pt(), mcPart1->Eta(), mcPart1->Phi(), AliPID::ParticleMass(AliPID::kElectron));
        Lvec2.SetPtEtaPhiM(mcPart2->Pt(), mcPart2->Eta(), mcPart2->Phi(), AliPID::ParticleMass(AliPID::kElectron));
      }
      else if (SmearedPair == kTRUE) {
        Lvec1.SetPtEtaPhiM(vec_negParticle->at(neg_i).fPt_smeared, vec_negParticle->at(neg_i).fEta_smeared, vec_negParticle->at(neg_i).fPhi_smeared, AliPID::ParticleMass(AliPID::kElectron));
        Lvec2.SetPtEtaPhiM(vec_posParticle->at(pos_i).fPt_smeared, vec_posParticle->at(pos_i).fEta_smeared, vec_posParticle->at(pos_i).fPhi_smeared, AliPID::ParticleMass(AliPID::kElectron));
      }
      TLorentzVector LvecM = Lvec1 + Lvec2;
      double mass = LvecM.M();
      double pairpt = LvecM.Pt();
      double weight = 1.;

      if (fDoMassCut == kTRUE && PartPrimary == kTRUE  && (fLowerMassCutPrimaries >= mass || mass >= fUpperMassCutPrimaries)) continue; // Apply mass cuts for primary pairs
      if (fDoMassCut == kTRUE && PartPrimary == kFALSE && mass >= fMassCutSecondaries)     continue; // Apply mass cuts for  secondary pairs

      if (fCocktailFile) {
        if (vec_negParticle->at(neg_i).GetMotherID() == vec_posParticle->at(pos_i).GetMotherID()){
          double motherpt = fMC->GetTrack(vec_negParticle->at(neg_i).GetMotherID())->Pt();
          weight *= GetWeight(vec_negParticle->at(neg_i), vec_posParticle->at(pos_i), motherpt);
        }
        else{
          weight = 0; // if should not fail by definition. but does in 13 / 10000000 cases
        }
      }

                                                                                // // Check particles mother and grandmother
                                                                                // if (fdebug) {
                                                                                //   int label1 = mcPart1->GetLabel();
                                                                                //   int label2 = mcPart2->GetLabel();
                                                                                //   int abslabel1 = TMath::Abs(label1);
                                                                                //   int abslabel2 = TMath::Abs(label2);
                                                                                //   int mother1 = TMath::Abs(mcPart1->GetMother());
                                                                                //   int mother2 = TMath::Abs(mcPart2->GetMother());
                                                                                //   int grandMotherID1 = TMath::Abs(fMC->GetTrack(mother1)->GetMother());
                                                                                //   int grandMotherID2 = TMath::Abs(fMC->GetTrack(mother2)->GetMother());
                                                                                //   // int nDauthersGM1 =fMC->GetTrack(grandMotherID1)->GetNDaughters();
                                                                                //   // int nDauthersGM2 =fMC->GetTrack(grandMotherID2)->GetNDaughters();
                                                                                //   if (fMC->GetTrack(grandMotherID1)->PdgCode() != 221 || fMC->GetTrack(grandMotherID2)->PdgCode() != 221 ) {
                                                                                //     std::cout << __LINE__ << " DEBUG_AnalysisTask: Neg Part. Label = " << label1 << ", PdgCode = " << mcPart1->PdgCode() << ", PdgCode Mother = " << fMC->GetTrack(mother1)->PdgCode() << ", PdgCode GrandMother = " << fMC->GetTrack(grandMotherID1)->PdgCode() << std::endl;
                                                                                //     std::cout << __LINE__ << " DEBUG_AnalysisTask: Pos Part. Label = " << label2 << ", PdgCode = " << mcPart2->PdgCode() << ", PdgCode Mother = " << fMC->GetTrack(mother2)->PdgCode() << ", PdgCode GrandMother = " << fMC->GetTrack(grandMotherID2)->PdgCode() << std::endl;
                                                                                //   }
                                                                                // }


      // ##########################################################
      // create MotherParticle from TwoPair that is used in 4 pairing and to summarize all the data
      // AliVParticle MotherParticle;
      double  eta     = LvecM.Eta();
      double  phi     = LvecM.Phi();
      short   charge  = mcPart1->Charge() + mcPart2->Charge();
      TwoPair MotherParticle(pairpt, eta, phi, mass, charge);

      MotherParticle.SetDaughterTrackID(vec_negParticle->at(neg_i).GetTrackID(), vec_posParticle->at(pos_i).GetTrackID());  // first: neg track, second: pos track

      if      (PartPrimary == kTRUE)   MotherParticle.SetMCTwoSignal_acc_prim(mcTwoSignal_acc_primary);
      else if (PartPrimary == kFALSE)  MotherParticle.SetMCTwoSignal_acc_sec(mcTwoSignal_acc_secondary);

      // std::cout << __LINE__ << " DEBUG_AnalysisTask: Mother Particlel: Pt = " << MotherParticle.fPt <<
      // ", Eta = " << MotherParticle.fEta <<
      // ", Phi = " << MotherParticle.fPhi <<
      // ", Mass = " << MotherParticle.fMass <<
      // ", charge = " << MotherParticle.fCharge << std::endl;
      if (SmearedPair == kFALSE) {
        if(PartPrimary == kTRUE)  fGenPairVec_primary.push_back(MotherParticle);
        if(PartPrimary == kFALSE) fGenPairVec_secondary.push_back(MotherParticle);
      }
      else if (SmearedPair == kTRUE) {
        if(PartPrimary == kTRUE)  fGenSmearedPairVec_primary.push_back(MotherParticle);
        if(PartPrimary == kFALSE) fGenSmearedPairVec_secondary.push_back(MotherParticle);
      }


      float ptPos       = -999;
      float etaPos      = -999;
      float phiPos      = -999;
      float ptNeg       = -999;
      float etaNeg      = -999;
      float phiNeg      = -999;
      float op_angle    = -999;
      //  Filling primary generated histos
      for (unsigned int i = 0; i < mcTwoSignal_acc_primary.size(); ++i){
        if (mcTwoSignal_acc_primary[i] == kTRUE){
          if (SmearedPair == kFALSE) {
            if (PartPrimary == kTRUE)   fHistGenPrimaryPair.at(i)->Fill(mass, pairpt, weight * centralityWeight);
          }
          else if (SmearedPair == kTRUE) {
            if (PartPrimary == kTRUE)   fHistGenSmearedPrimaryPair.at(i)->Fill(mass, pairpt, weight * centralityWeight);

            if (fWriteLegsFromPair){
              ptNeg  = vec_negParticle->at(neg_i).fPt_smeared;
              etaNeg = vec_negParticle->at(neg_i).fEta_smeared;
              phiNeg = vec_negParticle->at(neg_i).fPhi_smeared;
              ptPos  = vec_posParticle->at(pos_i).fPt_smeared;
              etaPos = vec_posParticle->at(pos_i).fEta_smeared;
              phiPos = vec_posParticle->at(pos_i).fPhi_smeared;
              op_angle = Lvec2.Angle(Lvec1.Vect());

              double tuple[7] = {ptPos,etaPos,phiPos,ptNeg,etaNeg,phiNeg,op_angle};
              fTHnSparseGenSmearedLegsFromPrimaryPair[i]->Fill(tuple);
            }
          }
        }
      } // end of loop over all primary MCsignals


      //  Filling secondary generated histos
      for (unsigned int i = 0; i < mcTwoSignal_acc_secondary.size(); ++i){
        if (mcTwoSignal_acc_secondary[i] == kTRUE){
          if (SmearedPair == kFALSE) {
            if (PartPrimary == kFALSE) fHistGenSecondaryPair.at(i)->Fill(mass, pairpt, weight * centralityWeight);
          }
          else if (SmearedPair == kTRUE) {
            if (PartPrimary == kFALSE) fHistGenSmearedSecondaryPair.at(i)->Fill(mass, pairpt, weight * centralityWeight);

            if (fWriteLegsFromPair){
              ptNeg  = vec_negParticle->at(neg_i).fPt_smeared;
              etaNeg = vec_negParticle->at(neg_i).fEta_smeared;
              phiNeg = vec_negParticle->at(neg_i).fPhi_smeared;
              ptPos  = vec_posParticle->at(pos_i).fPt_smeared;
              etaPos = vec_posParticle->at(pos_i).fEta_smeared;
              phiPos = vec_posParticle->at(pos_i).fPhi_smeared;
              op_angle = Lvec2.Angle(Lvec1.Vect());

              double tuple[7] = {ptPos,etaPos,phiPos,ptNeg,etaNeg,phiNeg,op_angle};
              fTHnSparseGenSmearedLegsFromSecondaryPair[i]->Fill(tuple);
            }
          }
        }
      } // end of loop over all MCsignals
    } // end of loop over all positive particles
  } // end of loop over all negative particles
} // end of DoGenAndGenSmearTwoPairing


// Fill reconstructed  pairs
void AliAnalysisTaskEtaReconstruction::DoRecTwoPairing(std::vector<Particle> fRecNegPart, std::vector<Particle> fRecPosPart, std::vector<AliDielectronSignalMC> fPairMCSignal, Bool_t PartPrimary, double centralityWeight){
  for (unsigned int neg_i = 0; neg_i < fRecNegPart.size(); ++neg_i){
    for (unsigned int pos_i = 0; pos_i < fRecPosPart.size(); ++pos_i){
                                                                                // if(fdebug) std::cout << __LINE__ << " DEBUG_AnalysisTask: fRecNegPart.size = " << fRecNegPart.size() << " fRecPosPart.size = " << fRecPosPart.size() << std::endl;
                                                                                // if(fdebug) {
                                                                                //   for (size_t neg_i = 0; neg_i < fRecNegPart.size(); neg_i++) {
                                                                                //     AliVParticle* track1  = fEvent->GetTrack(fRecNegPart[neg_i].GetTrackID());
                                                                                //     std::cout << __LINE__ << " DEBUG_AnalysisTask: Track1 Label = " << track1->GetLabel() << " and TrackID = " << fRecNegPart[neg_i].GetTrackID() << std::endl;
                                                                                //   }
                                                                                //   for (size_t pos_i = 0; pos_i < fRecPosPart.size(); pos_i++) {
                                                                                //     AliVParticle* track2  = fEvent->GetTrack(fRecPosPart[pos_i].GetTrackID());
                                                                                //     std::cout << __LINE__ << " DEBUG_AnalysisTask: Track2 Label = " << track2->GetLabel() << " and TrackID = " << fRecPosPart[pos_i].GetTrackID() << std::endl;
                                                                                //   }
                                                                                // }
      AliVParticle* track1  = fEvent->GetTrack(fRecNegPart[neg_i].GetTrackID());
      AliVParticle* track2  = fEvent->GetTrack(fRecPosPart[pos_i].GetTrackID());

      // Check if electrons are from MCSignal Generator
      if (!fRecPosPart[pos_i].GetMCSignalPair() || !fRecNegPart[neg_i].GetMCSignalPair()) continue;

      // Apply MC signals
      std::vector<Bool_t> mcTwoSignal_acc(fPairMCSignal.size(), kFALSE); // vector which stores if track is accepted by [i]-th mcsignal

      // Check if it according to mcsignals
      AliDielectronPair *pair = new AliDielectronPair;

      pair->SetKFUsage(false);
      pair->SetTracks(static_cast<AliVTrack*>(track1), 11, static_cast<AliVTrack*>(track2), -11);
      for (unsigned int i = 0; i < fPairMCSignal.size(); ++i){
        mcTwoSignal_acc[i] = AliDielectronMC::Instance()->IsMCTruth(pair, &(fPairMCSignal[i]));
      }
      // check if at least one mc signal is true
      if (CheckIfOneIsTrue(mcTwoSignal_acc) == kFALSE) continue;

      std::vector<bool> selected_primary(fPairCuts_primary.size(), kFALSE); // vector which stores if track is accepted by [i]-th selection cut
      // std::vector<bool> selected_secondary(fPairCuts_secondary_standard.size(), kFALSE); // vector which stores if track is accepted by [i]-th selection cut
      if (PartPrimary == kTRUE) {
        //  ---------- RECONSTRUCTED for PAIR  ---------- //
        // Check if particle is passing primary selection cuts
        for (UInt_t iCut=0; iCut<fPairCuts_primary.size(); ++iCut){ // loop over all specified cutInstances
          UInt_t selectedMask_primary=( 1 << fPairCuts_primary.at(iCut)->GetCuts()->GetEntries())-1;
          // cutting logic taken from AliDielectron::FillTrackArrays()          FillPairArrays()
          // apply track cuts
          UInt_t cutMask_primary = fPairCuts_primary.at(iCut)->IsSelected(pair);
          if (cutMask_primary == selectedMask_primary) {selected_primary[iCut] = kTRUE; /*std::cout << "sec_reconstructed TRUE" << std::endl;*/}
        }                                                                      // if (fdebug) std::cout << __LINE__ << " DEBUG_AnalysisTask: Size of Vectors after applying standard cuts, " << std::endl <<
                                                                                // " fRecNegPart_primary = " << fRecNegPart_primary.size() << std::endl <<
                                                                                // " fRecPosPart_primary = " << fRecPosPart_primary.size() << std::endl <<
                                                                                // /*" fGenPairVec_primary = "   << fGenPairVec_primary.size()   << " fGenSmearedPairVec_primary =   " << fGenSmearedPairVec_primary.size()   << */
                                                                                // " fRecPairVec_primary =   " << fRecPairVec_primary.size()   << std::endl <<
                                                                                // /*" fGenPairVec_secondary = " << fGenPairVec_secondary.size() << " fGenSmearedPairVec_secondary = " << fGenSmearedPairVec_secondary.size() << */
                                                                                // " fRecPairVec_secondary = " << fRecV0Pair.size()/*fRecPairVec_secondary.size()*/ << std::endl;
        // ##########################################################
        // check if at least one is selected by cuts otherwise skip this particle
        if (CheckIfOneIsTrue(selected_primary) == kFALSE) continue;
      }

      // else if (PartPrimary == kFALSE) {
      //   //  ---------- RECONSTRUCTED for PAIR  ---------- //
      //   // Check if particle is passing secondary selection cuts
      //   for (UInt_t iCut=0; iCut<fPairCuts_secondary_standard.size(); ++iCut){ // loop over all specified cutInstances
      //     UInt_t selectedMask_secondary=( 1 << fPairCuts_secondary_standard.at(iCut)->GetCuts()->GetEntries())-1;
      //     // cutting logic taken from AliDielectron::FillTrackArrays()          FillPairArrays()
      //     // apply track cuts
      //     UInt_t cutMask_secondary = fPairCuts_secondary_standard.at(iCut)->IsSelected(pair);
      //     if (cutMask_secondary == selectedMask_secondary) {selected_secondary[iCut] = kTRUE; /*std::cout << "sec_reconstructed TRUE" << std::endl;*/}
      //   }
      //   // ##########################################################
      //   // check if at least one is selected by cuts otherwise skip this particle
      //   if (CheckIfOneIsTrue(selected_secondary) == kFALSE) continue;
      // }
      // Construct pair variables from LorentzVectors
      TLorentzVector Lvec1;
      TLorentzVector Lvec2;
      Lvec1.SetPtEtaPhiM(track1->Pt(), track1->Eta(), track1->Phi(), AliPID::ParticleMass(AliPID::kElectron));
      Lvec2.SetPtEtaPhiM(track2->Pt(), track2->Eta(), track2->Phi(), AliPID::ParticleMass(AliPID::kElectron));
      TLorentzVector LvecM = Lvec1 + Lvec2;
      double mass = LvecM.M();
      if (fDoMassCut == kTRUE && PartPrimary == kTRUE  && (fLowerMassCutPrimaries >= mass || mass >= fUpperMassCutPrimaries)) continue; // mass cut for primaries
      // if (fDoMassCut == kTRUE && PartPrimary == kFALSE && mass >= fMassCutSecondaries   ) continue; // mass cut for primaries
      double pairpt = LvecM.Pt();
      double weight = 1.;
      if (fCocktailFile) {
        if (fRecNegPart[neg_i].GetMotherID() == fRecPosPart[pos_i].GetMotherID()){
          double motherpt = fMC->GetTrack(fRecNegPart[neg_i].GetMotherID())->Pt();
          weight *= GetWeight(fRecNegPart[neg_i], fRecPosPart[pos_i], motherpt);
        }
        else{
          weight = 0; // if should not fail by definition. but does in 13 / 10000000 cases
        }
      }
                                                                                      // if(fdebug) std::cout << __LINE__ << " DEBUG_AnalysisTask: Line cout " << std::endl;
                                                                                // Check particles mother and grandmother if mass > 400 MeV/c^2 and pt < 400 MeV/c
                                                                                // if (fdebug) {
                                                                                  // std::cout << __LINE__ << " DEBUG_AnalysisTask: fEvent->TrackID = " << fEvent->GetTrack(fRecNegPart[neg_i].GetTrackID())->GetLabel() << ", fMC->TrackID = " << fMC->GetTrack(fRecNegPart[neg_i].GetTrackID())->GetLabel() << std::endl;
                                                                                  // if (mass > 0.160 && pairpt < 0.200) {
                                                                                  // int label1 = track1->GetLabel();
                                                                                  // int label2 = track2->GetLabel();
                                                                                  // int abslabel1 = TMath::Abs(label1);
                                                                                  // int abslabel2 = TMath::Abs(label2);
                                                                                  // int motherID1 = TMath::Abs(fMC->GetTrack(abslabel1)->GetMother());
                                                                                  // int motherID2 = TMath::Abs(fMC->GetTrack(abslabel2)->GetMother());
                                                                                  // int grandMotherID1 = TMath::Abs(fMC->GetTrack(motherID1)->GetMother());
                                                                                  // int grandMotherID2 = TMath::Abs(fMC->GetTrack(motherID2)->GetMother());
                                                                                  // int NDauthersGM1 =fMC->GetTrack(grandMotherID1)->GetNDaughters();
                                                                                  // int NDauthersGM2 =fMC->GetTrack(grandMotherID2)->GetNDaughters();
                                                                                  //
                                                                                  //   std::cout << __LINE__ << " DEBUG_AnalysisTask: Cout if particle has mass > 160 MeV/c^2 and pt < 200 MeV/c" << std::endl;
                                                                                  //   std::cout << __LINE__ << "                      ------------------------------------                              " << std::endl;
                                                                                  // if (fMC->GetTrack(grandMotherID1)->PdgCode() != 221 ||fMC->GetTrack(grandMotherID2)->PdgCode() != 221) {
                                                                                  //   std::cout << __LINE__ << " DEBUG_AnalysisTask: Neg Part. fEvent Label = " << label1 << ", PdgCode = " << fMC->GetTrack(abslabel1)->PdgCode() << ", PdgCode Mother1 = " << fMC->GetTrack(motherID1)->PdgCode() << ", PdgCode GrandMother1 = " << fMC->GetTrack(grandMotherID1)->PdgCode() << std::endl;
                                                                                  //   std::cout << __LINE__ << " DEBUG_AnalysisTask: Pos Part. fEvent Label = " << label2 << ", PdgCode = " << fMC->GetTrack(abslabel2)->PdgCode() << ", PdgCode Mother2 = " << fMC->GetTrack(motherID2)->PdgCode() << ", PdgCode GrandMother2 = " << fMC->GetTrack(grandMotherID2)->PdgCode() << std::endl;
                                                                                  // }
                                                                                  //   std::cout << __LINE__ << "                      ------------------------------------                              " << std::endl;
                                                                                  //
                                                                                  //   std::cout << __LINE__ << " DEBUG_AnalysisTask: GrandMother1 :: PdgCode = " << fMC->GetTrack(grandMotherID1)->PdgCode() << ", Number of Daughters of GrandMother1 = " << NDauthersGM1 << std::endl;
                                                                                  //   for (int NGDaughters = 0; NGDaughters < NDauthersGM1; NGDaughters++) {
                                                                                  //     int GDaughterID = fMC->GetTrack(grandMotherID1)->GetDaughterFirst()+NGDaughters;
                                                                                  //     int NDauthersM =fMC->GetTrack(GDaughterID)->GetNDaughters(); // Anzahl der Toechter der Tochter der GrandMother (toechter der mutter)
                                                                                  //     std::cout << __LINE__ << " DEBUG_AnalysisTask: Mother :: " << NGDaughters+1 <<  " Daughter of GrandMother: PdgCode = " << fMC->GetTrack(GDaughterID)->PdgCode() << ", Number of Daughters from Mother = " << NDauthersM << std::endl;
                                                                                  //     for (int NMDaughters = 0; NMDaughters < NDauthersM; NMDaughters++) {
                                                                                  //       int MDaughterID = fMC->GetTrack(GDaughterID)->GetDaughterFirst()+NMDaughters;
                                                                                  //       std::cout << __LINE__ << " DEBUG_AnalysisTask: Particle " << NMDaughters+1 <<  " :: Daughter of Mother: PdgCode = " << fMC->GetTrack(MDaughterID)->PdgCode() << std::endl;
                                                                                  //     }
                                                                                  //   }
                                                                                  //
                                                                                  //   std::cout << __LINE__ << " DEBUG_AnalysisTask: GrandMother2 :: PdgCode = " << fMC->GetTrack(grandMotherID2)->PdgCode() << ", Number of Daughters of GrandMother2 = " << NDauthersGM2 << std::endl;
                                                                                  //   for (int NGDaughters = 0; NGDaughters < NDauthersGM2; NGDaughters++) {
                                                                                  //     int GDaughterID = fMC->GetTrack(grandMotherID2)->GetDaughterFirst()+NGDaughters;
                                                                                  //     int NDauthersM =fMC->GetTrack(GDaughterID)->GetNDaughters();
                                                                                  //     std::cout << __LINE__ << " DEBUG_AnalysisTask: Mother :: " << NGDaughters+1 <<  " Daughter of GrandMother: PdgCode = " << fMC->GetTrack(GDaughterID)->PdgCode() << ", Number of Daughters from Mother = " << NDauthersM << std::endl;
                                                                                  //     for (int NMDaughters = 0; NMDaughters < NDauthersM; NMDaughters++) {
                                                                                  //       int MDaughterID = fMC->GetTrack(GDaughterID)->GetDaughterFirst()+NMDaughters;
                                                                                  //       std::cout << __LINE__ << " DEBUG_AnalysisTask: Particle :: " << NMDaughters+1 <<  " Daughter of Mother: PdgCode = " << fMC->GetTrack(MDaughterID)->PdgCode() << std::endl;
                                                                                  //     }
                                                                                  //   }
                                                                                  // }
                                                                                // }

      // ##########################################################
      // create MotherParticle from TwoPair that is used in 4 pairing and to summarize all the data
      // int motherID = TMath::Abs(fMC->GetTrack(iPart)->GetMother());


      // AliVParticle MotherParticle;
      double  eta     = LvecM.Eta();
      double  phi     = LvecM.Phi();
      short   charge  = track1->Charge() + track2->Charge();
      TwoPair MotherParticle(pairpt, eta, phi, mass, charge);
      MotherParticle.SetDaughterTrackID(fRecNegPart[neg_i].GetTrackID(), fRecPosPart[pos_i].GetTrackID());  // first: neg track, second: pos track

      if (PartPrimary == kTRUE) {
        MotherParticle.isReconstructed_primary   = selected_primary;
        // MotherParticle.SetDauthersAreReconstructed(fRecNegPart[neg_i].isReconstructed_primary, fRecPosPart[pos_i].isReconstructed_primary);
        MotherParticle.SetMCTwoSignal_acc_prim(mcTwoSignal_acc);
        fRecPairVec_primary.push_back(MotherParticle);
      }
      // else if (PartPrimary == kFALSE) {
      //   MotherParticle.isReconstructed_secondary = selected_secondary;
      //   MotherParticle.SetDauthersAreReconstructed(fRecNegPart[neg_i].isReconstructed_secondary, fRecPosPart[pos_i].isReconstructed_secondary);
      //   MotherParticle.SetMCTwoSignal_acc_sec(mcTwoSignal_acc);
      //   fRecPairVec_secondary.push_back(MotherParticle);
      // }

                                                                                // if(fdebug){ std::cout << __LINE__ << " DEBUG_AnalysisTask: Mother Particlel: Pt = " << MotherParticle.fPt <<
                                                                                // ", Eta = " << MotherParticle.fEta <<
                                                                                // ", Phi = " << MotherParticle.fPhi <<
                                                                                // ", Mass = " << MotherParticle.fMass <<
                                                                                // ", charge = " << MotherParticle.fCharge << std::endl;}

      // // ##########################################################
      // float ptPos       = -999;
      // float etaPos      = -999;
      // float phiPos      = -999;
      // float ptNeg       = -999;
      // float etaNeg      = -999;
      // float phiNeg      = -999;
      // float op_angle    = -999;
      // for (unsigned int i = 0; i < mcTwoSignal_acc.size(); ++i){
      //
      //   // Filling reconstructed primary particle histograms according to MCSignals
      //   if (mcTwoSignal_acc[i] == kTRUE && PartPrimary == kTRUE){
      //     for (unsigned int j = 0; j < fRecNegPart[neg_i].isReconstructed_primary.size(); ++j){
      //       if (fRecNegPart[neg_i].isReconstructed_primary[j] == kTRUE && fRecPosPart[pos_i].isReconstructed_primary[j] == kTRUE){
      //         fHistRecPrimaryPair.at(j * mcTwoSignal_acc.size() + i)->Fill(mass, pairpt, weight * centralityWeight);
      //
      //         if (fWriteLegsFromPair){
      //           ptNeg  = fRecNegPart[neg_i].fPt;
      //           etaNeg = fRecNegPart[neg_i].fEta;
      //           phiNeg = fRecNegPart[neg_i].fPhi;
      //           ptPos  = fRecPosPart[pos_i].fPt;
      //           etaPos = fRecPosPart[pos_i].fEta;
      //           phiPos = fRecPosPart[pos_i].fPhi;
      //           op_angle = Lvec2.Angle(Lvec1.Vect());
      //
      //           const double tuple[7] = {ptPos,etaPos,phiPos,ptNeg,etaNeg,phiNeg,op_angle};
      //           fTHnSparseRecLegsFromPrimaryPair.at(j * mcTwoSignal_acc.size() + i)->Fill(tuple);
      //         }
      //       }// is selected by cutsetting
      //     } // end of loop over all cutsettings
      //   } // is selected by MCSignal
      //
      //   // Filling reconstructed secondary particle histograms according to MCSignals
      //   if (mcTwoSignal_acc[i] == kTRUE && PartPrimary == kFALSE){
      //     for (unsigned int j = 0; j < fRecNegPart[neg_i].isReconstructed_secondary.size(); ++j){
      //       if (fRecNegPart[neg_i].isReconstructed_secondary[j] == kTRUE && fRecPosPart[pos_i].isReconstructed_secondary[j] == kTRUE){
      //         fHistRecSecondaryPair.at(j * mcTwoSignal_acc.size() + i)->Fill(mass, pairpt, weight * centralityWeight);
      //
      //         if (fWriteLegsFromPair){
      //           ptNeg  = fRecNegPart[neg_i].fPt;
      //           etaNeg = fRecNegPart[neg_i].fEta;
      //           phiNeg = fRecNegPart[neg_i].fPhi;
      //           ptPos  = fRecPosPart[pos_i].fPt;
      //           etaPos = fRecPosPart[pos_i].fEta;
      //           phiPos = fRecPosPart[pos_i].fPhi;
      //           op_angle = Lvec2.Angle(Lvec1.Vect());
      //
      //           const double tuple[7] = {ptPos,etaPos,phiPos,ptNeg,etaNeg,phiNeg,op_angle};
      //           fTHnSparseRecLegsFromSecondaryPair.at(j * mcTwoSignal_acc.size() + i)->Fill(tuple);
      //         }
      //       }// is selected by cutsetting
      //     } // end of loop over all cutsettings
      //   } // is selected by MCSignal
      // } // end of loop over all MCsignals
      delete pair;
    }// end of positive particle loop
  } // end of negative particle loop
}


void AliAnalysisTaskEtaReconstruction::DoRecTwoPairingV0(std::vector<AliDielectronSignalMC> fPairMCSignal){
  // ######################### V0 Secondaries #####################################
  //
  for (Int_t iV0 = 0; iV0 < fEvent->GetNumberOfV0s(); iV0++){
                                                                                // if (fdebug) std::cout << __LINE__ << " DEBUG_AnalysisTask: GetNumberOfV0s = " << fEvent->GetNumberOfV0s() << std::endl;
    // ##########################################################
    // Track handling
    if (isAOD) {
      AliAODv0* fV0 = ((AliAODEvent*)fEvent)->GetV0(iV0);
      // if (!fV0) { Printf("ERROR: Could not receive V0 %d", iV0); continue; }
      // // if (isAOD) track = static_cast<AliAODTrack*>(track);
      // else       track = static_cast<AliESDtrack*>(track);
      // int v0label = fV0->GetLabel();
      // int absv0label = TMath::Abs(v0label);

      Bool_t fV0Status = fV0 ->GetOnFlyStatus();
      if (fV0Status != fV0OnFlyStatus) continue;

      // new part
      AliAODTrack *posAODTrack = (AliAODTrack *) (fV0->GetSecondaryVtx()->GetDaughter(0));
      AliAODTrack *negAODTrack = (AliAODTrack *) (fV0->GetSecondaryVtx()->GetDaughter(1));


      // Apply MC signals
      std::vector<Bool_t> mcTwoSignal_acc(fPairMCSignal.size(), kFALSE); // vector which stores if track is accepted by [i]-th mcsignal

      // Check if it according to mcsignals
      AliDielectronPair *pair = new AliDielectronPair;
      pair->SetKFUsage(false);
      pair->SetTracks(static_cast<AliAODTrack*>(negAODTrack), 11, static_cast<AliAODTrack*>(posAODTrack), -11);
      for (unsigned int i = 0; i < fPairMCSignal.size(); ++i){
        mcTwoSignal_acc[i] = AliDielectronMC::Instance()->IsMCTruth(pair, &(fPairMCSignal[i]));
      }
      // check if at least one mc signal is true
      if (CheckIfOneIsTrue(mcTwoSignal_acc) == kFALSE) continue;



      //  ---------- RECONSTRUCTED for PAIR  ---------- //
      // Check if particle is passing secondary selection cuts

      std::vector<bool> selected_secondary(fPairCuts_secondary_PreFilter.size(), kFALSE); // vector which stores if pair is accepted by [i]-th selection cut
      for (UInt_t iCut=0; iCut<fPairCuts_secondary_PreFilter.size(); ++iCut){ // loop over all specified cutInstances
        UInt_t selectedMask_secondary=( 1 << fPairCuts_secondary_PreFilter.at(iCut)->GetCuts()->GetEntries())-1;
        // cutting logic taken from AliDielectron::FillTrackArrays()          FillPairArrays()
        // apply pair cuts
        UInt_t cutMask_secondary = fPairCuts_secondary_PreFilter.at(iCut)->IsSelected(pair);
        if (cutMask_secondary == selectedMask_secondary) {selected_secondary[iCut] = kTRUE; /*std::cout << "sec_reconstructed TRUE" << std::endl;*/}
      }



      // // //  ---------- RECONSTRUCTED for TRACKS  ---------- //
      // std::vector<bool> selected_secondary_firstPar(fTrackCuts_secondary_PreFilter.size(), kFALSE); // vector which stores if track is accepted by [i]-th selection cut
      // std::vector<bool> selected_secondary_secPar  (fTrackCuts_secondary_PreFilter.size(), kFALSE); // vector which stores if track is accepted by [i]-th selection cut
      // for (UInt_t iCut=0; iCut<fTrackCuts_secondary_PreFilter.size(); ++iCut){ // loop over all specified cutInstances
      //   UInt_t selectedMask_secondary_firstPar=( 1 << fTrackCuts_secondary_PreFilter.at(iCut)->GetCuts()->GetEntries())-1;
      //   UInt_t selectedMask_secondary_secPar=( 1 << fTrackCuts_secondary_PreFilter.at(iCut)->GetCuts()->GetEntries())-1;
      //   // cutting logic taken from AliDielectron::FillTrackArrays()          FillPairArrays()
      //   // apply track cuts
      //   UInt_t cutMask_secondary_firstPar = fTrackCuts_secondary_PreFilter.at(iCut)->IsSelected(posAODTrack);
      //   UInt_t cutMask_secondary_secPar = fTrackCuts_secondary_PreFilter.at(iCut)->IsSelected(negAODTrack);
      //   if (cutMask_secondary_firstPar == selectedMask_secondary_firstPar) {selected_secondary_firstPar[iCut] = kTRUE; /*std::cout << "sec_reconstructed TRUE" << std::endl;*/}
      //   if (cutMask_secondary_secPar == selectedMask_secondary_secPar) {selected_secondary_secPar[iCut] = kTRUE; /*std::cout << "sec_reconstructed TRUE" << std::endl;*/}
      // }


      // ##########################################################
      // check if at least one is selected by cuts otherwise skip this particle
      if (CheckIfOneIsTrue(selected_secondary) == kFALSE) continue;
      // if (CheckIfOneIsTrue(selected_secondary_firstPar) == kFALSE && CheckIfOneIsTrue(selected_secondary_secPar) == kFALSE) continue;




                                                                              // AliVParticle* track1  = fEvent->GetTrack(fRecNegPart[neg_i].GetTrackID());
                                                                              // int label1 = track1->GetLabel();
                                                                              // int abslabel1 = TMath::Abs(label1);
                                                                              // int motherID1 = TMath::Abs(fMC->GetTrack(abslabel1)->GetMother());
                                                                              // int grandMotherID1 = TMath::Abs(fMC->GetTrack(motherID1)->GetMother());





                                                                              // if (fdebug) std::cout << __LINE__ << " DEBUG_AnalysisTask: fV0 Label: " << v0label << " PosID: " << posTrackID << " PosLabel: " << posTrackLabel << " NegID: " << negTrackID << " NegLabel: " << negTrackLabel << std::endl;

                                                                              // Check particles mother and grandmother
                                                                              // std::cout << __LINE__ << " DEBUG_AnalysisTask: fEvent->TrackID = " << fEvent->GetTrack(fRecNegPart[neg_i].GetTrackID())->GetLabel() << ", fMC->TrackID = " << fMC->GetTrack(fRecNegPart[neg_i].GetTrackID())->GetLabel() << std::endl;
                                                                              // std::cout << __LINE__ << /*" DEBUG_AnalysisTask: fEvent->TrackLabel = " << ((AliAODTrack*)fEvent->GetTrack(negTrackID))->GetID() <<*/ ", fMC->TrackLabel = " << fMC->GetTrack(negTrackID)->GetLabel() << std::endl;



                                                                              // // track label, points back to MC track
                                                                              // Int_t posTrackLabel = TMath::Abs(posAODTrack->GetLabel());
                                                                              // Int_t negTrackLabel = TMath::Abs(negAODTrack->GetLabel());
                                                                              // // unique track ID, points back to the ESD track
                                                                              // // Short_t posTrackID = posAODTrack->GetID();
                                                                              // // Short_t negTrackID = negAODTrack->GetID();
                                                                              // int label1 = fMC->GetTrack(negTrackLabel)->GetLabel();
                                                                              // int label2 = fMC->GetTrack(posTrackLabel)->GetLabel();
                                                                              // // int abslabel1 = TMath::Abs(label1);
                                                                              // // int abslabel2 = TMath::Abs(label2);
                                                                              // int motherLabel1 = TMath::Abs(fMC->GetTrack(negTrackLabel)->GetMother());
                                                                              // int motherLabel2 = TMath::Abs(fMC->GetTrack(posTrackLabel)->GetMother());
                                                                              // int grandMotherLabel1 = TMath::Abs(fMC->GetTrack(motherLabel1)->GetMother());
                                                                              // int grandMotherLabel2 = TMath::Abs(fMC->GetTrack(motherLabel2)->GetMother());
                                                                              // // int NDauthersGM1 =fMC->GetTrack(grandMotherID1)->GetNDaughters();
                                                                              // // int NDauthersGM2 =fMC->GetTrack(grandMotherID2)->GetNDaughters();
                                                                              //
                                                                              // int negPdgCode = fMC->GetTrack(negTrackLabel)->PdgCode();
                                                                              // int posPdgCode = fMC->GetTrack(posTrackLabel)->PdgCode();
                                                                              // int negMPdgCode = fMC->GetTrack(motherLabel1)->PdgCode();
                                                                              // int posMPdgCode = fMC->GetTrack(motherLabel2)->PdgCode();
                                                                              // int negGMPdgCode = fMC->GetTrack(grandMotherLabel1)->PdgCode();
                                                                              // int posGMPdgCode = fMC->GetTrack(grandMotherLabel2)->PdgCode();
                                                                              //
                                                                              //
                                                                              //
                                                                              // if(fdebug) std::cout << __LINE__ << " DEBUG_AnalysisTask: COUT LINE" << std::endl;
                                                                              // if (fdebug) {
                                                                              //   // if ((negGMPdgCode != 111 && negGMPdgCode != 221) || (posGMPdgCode != 111 && posGMPdgCode != 221)) {
                                                                              //   if (grandMotherLabel1 == grandMotherLabel2) {
                                                                              //     std::cout << __LINE__ << "                      ------------------------------------                              " << std::endl;
                                                                              //     std::cout << __LINE__ << " DEBUG_AnalysisTask: Neg Part.   PdgCode = " << negPdgCode << ", PdgCode Mother1 = " << negMPdgCode  << ", PdgCode GrandMother1 = " << negGMPdgCode << std::endl;
                                                                              //     std::cout << __LINE__ << " DEBUG_AnalysisTask: Pos Part.   PdgCode = " << posPdgCode << ", PdgCode Mother2 = " << posMPdgCode  << ", PdgCode GrandMother2 = " << posGMPdgCode << std::endl;
                                                                              //     std::cout << __LINE__ << " DEBUG_AnalysisTask: Neg Part. fMC Label = " << label1     << ", NegMotherLabel  = " << motherLabel1 << ", negGrandMotherLabel  = " << grandMotherLabel1 << std::endl;
                                                                              //     std::cout << __LINE__ << " DEBUG_AnalysisTask: Pos Part. fMC Label = " << label2     << ", PosMotherLabel  = " << motherLabel2 << ", posGrandMotherLabel  = " << grandMotherLabel2 << std::endl;
                                                                              //
                                                                              //     std::cout << __LINE__ << "                      ------------------------------------                              " << std::endl;
                                                                              //     int NDauthersGM1 =fMC->GetTrack(grandMotherLabel1)->GetNDaughters();
                                                                              //     // int NDauthersGM2 =fMC->GetTrack(grandMotherLabel2)->GetNDaughters();
                                                                              //     std::cout << __LINE__ << " DEBUG_AnalysisTask: GrandMother1 :: PdgCode = " << fMC->GetTrack(grandMotherLabel1)->PdgCode() << ", Number of Daughters of GrandMother1 = " << NDauthersGM1 << std::endl;
                                                                              //     for (int NGDaughters = 0; NGDaughters < NDauthersGM1; NGDaughters++) {
                                                                              //       int GDaughterID = fMC->GetTrack(grandMotherLabel1)->GetDaughterFirst()+NGDaughters;
                                                                              //       int NDauthersM =fMC->GetTrack(GDaughterID)->GetNDaughters(); // Anzahl der Toechter von der Tochter der GrandMother (toechter der mutter)
                                                                              //       std::cout << __LINE__ << " DEBUG_AnalysisTask: Mother :: " << NGDaughters+1 <<  " Daughter of GrandMother: PdgCode = " << fMC->GetTrack(GDaughterID)->PdgCode() << ", Number of Daughters from Mother = " << NDauthersM << std::endl;
                                                                              //       for (int NMDaughters = 0; NMDaughters < NDauthersM; NMDaughters++) {
                                                                              //         int MDaughterID = fMC->GetTrack(GDaughterID)->GetDaughterFirst()+NMDaughters;
                                                                              //         std::cout << __LINE__ << " DEBUG_AnalysisTask: Particle " << NMDaughters+1 <<  " :: Daughter of Mother: PdgCode = " << fMC->GetTrack(MDaughterID)->PdgCode() << std::endl;
                                                                              //       }
                                                                              //     }
                                                                              //
                                                                              //     // std::cout << __LINE__ << "                      ------------------------------------                              " << std::endl;
                                                                              //     // std::cout << __LINE__ << " DEBUG_AnalysisTask: GrandMother2 :: PdgCode = " << fMC->GetTrack(grandMotherLabel2)->PdgCode() << ", Number of Daughters of GrandMother2 = " << NDauthersGM2 << std::endl;
                                                                              //     // for (int NGDaughters = 0; NGDaughters < NDauthersGM2; NGDaughters++) {
                                                                              //       // int GDaughterID = fMC->GetTrack(grandMotherLabel2)->GetDaughterFirst()+NGDaughters;
                                                                              //       // int NDauthersM =fMC->GetTrack(GDaughterID)->GetNDaughters();
                                                                              //       // std::cout << __LINE__ << " DEBUG_AnalysisTask: Mother :: " << NGDaughters+1 <<  " Daughter of GrandMother: PdgCode = " << fMC->GetTrack(GDaughterID)->PdgCode() << ", Number of Daughters from Mother = " << NDauthersM << std::endl;
                                                                              //       // for (int NMDaughters = 0; NMDaughters < NDauthersM; NMDaughters++) {
                                                                              //         // int MDaughterID = fMC->GetTrack(GDaughterID)->GetDaughterFirst()+NMDaughters;
                                                                              //         // std::cout << __LINE__ << " DEBUG_AnalysisTask: Particle :: " << NMDaughters+1 <<  " Daughter of Mother: PdgCode = " << fMC->GetTrack(MDaughterID)->PdgCode() << std::endl;
                                                                              //       // }
                                                                              //     // }
                                                                              //
                                                                              //     std::cout << __LINE__ << "                      ====================================                              " << std::endl;
                                                                              //   }
                                                                              // int label1 = negAODTrack->GetLabel();
                                                                              // int label2 = posAODTrack->GetLabel();
                                                                              // int abslabel1 = TMath::Abs(label1);
                                                                              // int abslabel2 = TMath::Abs(label2);
                                                                              // int daughter1pdg = fMC->GetTrack(abslabel1)->PdgCode();
                                                                              // int daughter2pdg = fMC->GetTrack(abslabel2)->PdgCode();
                                                                              // int motherID1 = TMath::Abs(fMC->GetTrack(abslabel1)->GetMother());
                                                                              // int motherID2 = TMath::Abs(fMC->GetTrack(abslabel2)->GetMother());
                                                                              // int mother1pdg = fMC->GetTrack(motherID1)->PdgCode();
                                                                              // int mother2pdg = fMC->GetTrack(motherID2)->PdgCode();
                                                                              // if(fdebug && (TMath::Abs(daughter1pdg) != 11 || TMath::Abs(daughter2pdg) != 11))  std::cout << __LINE__ << " DEBUG_AnalysisTask: Daughters of V0 have PDGCode: Daughter1 = " << daughter1pdg << ", and Daughter2 = " << daughter2pdg << std::endl;
                                                                              // if(fdebug && (mother1pdg != 22 || mother2pdg != 22)) std::cout << __LINE__<< " Mother is no Photon: posMother1 = " << mother1pdg << ", negMother 2 = " << mother2pdg << std::endl;
                                                                              // }




      // Construct pair variables from LorentzVectors
      TLorentzVector Lvec1;
      TLorentzVector Lvec2;
      Lvec1.SetPtEtaPhiM(posAODTrack->Pt(), posAODTrack->Eta(), posAODTrack->Phi(), AliPID::ParticleMass(AliPID::kElectron));
      Lvec2.SetPtEtaPhiM(negAODTrack->Pt(), negAODTrack->Eta(), negAODTrack->Phi(), AliPID::ParticleMass(AliPID::kElectron));
      TLorentzVector LvecM = Lvec1 + Lvec2;
      LvecM.SetE(sqrt(LvecM.P()*LvecM.P() - fPhotonMass*fPhotonMass)); // Set Mass to photon mass = 0
      double mass = LvecM.M();
      if (fDoMassCut == kTRUE && mass >= fMassCutSecondaries   ) continue; // mass cut for primaries
      double pairpt = LvecM.Pt();
      // double weight = 1.;

      // AliVTrack MotherParticle;
      double  eta     = LvecM.Eta();
      double  phi     = LvecM.Phi();
      short   charge  = posAODTrack->Charge() + negAODTrack->Charge();
      TwoPair MotherParticle(pairpt, eta, phi, mass, charge);


      // MotherParticle.isReconstructed_secondary = selected_secondary;
      MotherParticle.SetV0ID(iV0);
      // MotherParticle.SetDaughterTrackID(negAODTrack->GetLabel(), posAODTrack->GetLabel()); // first: neg track, second: pos track
      MotherParticle.SetV0Status(fV0Status);
      // MotherParticle.SetDauthersAreReconstructed(fRecNegPart[neg_i].isReconstructed_secondary, fRecPosPart[pos_i].isReconstructed_secondary);
      // MotherParticle.SetMCTwoSignal_acc_sec(mcTwoSignal_acc);

                                                                              // if(fdebug) std::cout << __LINE__ << " DEBUG_AnalysisTask: COUT LINE" << std::endl;
                                                                              // if (fdebug) {
                                                                              //     std::cout << "V0: Label1 = " << negAODTrack->GetLabel() << " Label2 = " <<  posAODTrack->GetLabel() << std::endl;
                                                                              // }
      fRecV0Pair.push_back(MotherParticle);


      // // Fill secondary support histos
      // AliVParticle* mcPosPart = fMC->GetTrack(posTrackLabel);
      // AliVParticle* mcNegPart = fMC->GetTrack(negTrackLabel);
      //
      // // Filling reconstructed secondary particle histograms according to MCSignals
      //   for (unsigned int i = 0; i < MotherParticle.isReconstructed_secondary.size(); ++i){
      //     if (MotherParticle.isReconstructed_secondary[i] == kTRUE){
      //       for (unsigned int j = 0; j < mcTwoSignal_acc.size(); ++j){
      //         if (mcTwoSignal_acc[j] == kTRUE){
      //         fHistRecSecondaryPair.at(  i *   mcTwoSignal_acc.size() + j)->Fill(mass, pairpt, weight * centralityWeight);
      //         FillTrackHistograms_Secondary(posAODTrack, mcPosPart, i, j); // Fill secondary support histograms
      //         FillTrackHistograms_Secondary(negAODTrack, mcNegPart, i, j); // Fill secondary ssupport histograms
      //         FillPairHistograms_Secondary(pair, i, j);                    // Fill pair secondary histograms
      //       }// is selected by cutsetting
      //     } // end of loop over all cutsettings
      //   } // is selected by MCSignal
      // } // end of loop over all MCsignals

      delete pair;
    }
  }
}



void AliAnalysisTaskEtaReconstruction::DoFourPreFilter(std::vector<TwoPair>* fPairVec_primary, std::vector<TwoPair>* fPairVec_secondary){
  if(fPairVec_primary->size() == 0 || fPairVec_secondary->size() == 0) return;
  // if(fdebug) std::cout << " fPairVec_primary: "    << fPairVec_primary->size()   << std::endl;
  // if(fdebug) std::cout << " fPairVec_secondary: "  << fPairVec_secondary->size() << std::endl;
  if ((fPairVec_primary == fPairVec_secondary) && (fPairVec_secondary->size()==1)) {/*std::cout << "return cause vector has jsut one V0" << std::endl;*/
                                                                                    return;}
  // for (int prim_i = fPairVec_primary->size()-1; prim_i >= 0; --prim_i){
  //                                                                               if(fdebug) std::cout << __LINE__ << " DEBUG_AnalysisTask: COUT Loop parameter prim_i = " << prim_i << std::endl;
  //   for (int sec_i = fPairVec_secondary->size()-1; sec_i >= 0; --sec_i){
  //                                                                               if(fdebug) std::cout << __LINE__ << " DEBUG_AnalysisTask: COUT Loop parameter sec_i = " << sec_i << std::endl;
  //
  //     // Get LorentzVectors from reconstructed first pair (primary) and reconstructed V0 (secondary V0)
  //     TLorentzVector Lvec1;
  //     TLorentzVector Lvec2;
  //                                                                               if(fdebug) std::cout << __LINE__ << " DEBUG_AnalysisTask: COUT LINE "  << std::endl;
  //
  //     Lvec1.SetPtEtaPhiM(fPairVec_primary[prim_i].fPt , fPairVec_primary[prim_i].fEta , fPairVec_primary[prim_i].fPhi , fPairVec_primary[prim_i].fMass);
  //     Lvec2.SetPtEtaPhiM(fPairVec_secondary[sec_i].fPt, fPairVec_secondary[sec_i].fEta, fPairVec_secondary[sec_i].fPhi, fPairVec_secondary[sec_i].fMass);
  //                                                                               if(fdebug) std::cout << __LINE__ << " DEBUG_AnalysisTask: COUT LINE " << std::endl;
  //     TLorentzVector LvecM = Lvec1 + Lvec2;
  //     double mass = LvecM.M();
  //     if(mass > 0.15 || mass < 0.13) continue;
  //     else {
  //                                                                               if(fdebug) std::cout << __LINE__ << " DEBUG_AnalysisTask: COUT LINE " << std::endl;
  //       fPairVec_primary->erase(fPairVec_primary->begin()+prim_i);
  //       fPairVec_secondary->erase(fPairVec_secondary->begin()+sec_i);
  //       break;
  //     }
  //   } // end of loop sec_i
  // } // end of loop prim_i

                                                                                // if(fdebug) std::cout << __LINE__ << " DEBUG_AnalysisTask: COUT LINE" << std::endl;
  unsigned int sec_iStartValue = 0;
  for (unsigned int prim_i = 0; prim_i < fPairVec_primary->size(); ++prim_i){
                                                                                // if(fdebug) std::cout << __LINE__ << " DEBUG_AnalysisTask: COUT Loop parameter prim_i = " << prim_i << std::endl;
    if (fPairVec_primary == fPairVec_secondary) sec_iStartValue = prim_i+1;      // setting sec_i=prim_i+1 is avoiding pairing of same pair (V0)
    for (unsigned int sec_i = sec_iStartValue; sec_i < fPairVec_secondary->size(); ++sec_i){
                                                                                // if(fdebug) std::cout << __LINE__ << " DEBUG_AnalysisTask: COUT Loop parameter sec_i = " << sec_i << std::endl;

      // Use if case, when PreFilter is used with twice the same vector.
      // e.g. Try rejection of pions trough pairing of 2 secondary photon conversion pairs.    (Pi0 -> 2 photons -> e+ e-  e+ e-)
      if (fPairVec_primary == fPairVec_secondary) {
        Int_t fFirstV0ID  = fPairVec_primary->at(prim_i) .GetV0ID();
        Int_t fSecondV0ID = fPairVec_secondary->at(sec_i).GetV0ID();
        AliAODv0* fFirstV0 = ((AliAODEvent*)fEvent)->GetV0(fFirstV0ID);
        AliAODv0* fSecondV0 = ((AliAODEvent*)fEvent)->GetV0(fSecondV0ID);
        AliAODTrack *track1 = (AliAODTrack *) (fFirstV0 ->GetSecondaryVtx()->GetDaughter(0));
        AliAODTrack *track2 = (AliAODTrack *) (fFirstV0 ->GetSecondaryVtx()->GetDaughter(1));
        AliAODTrack *track3 = (AliAODTrack *) (fSecondV0->GetSecondaryVtx()->GetDaughter(0));
        AliAODTrack *track4 = (AliAODTrack *) (fSecondV0->GetSecondaryVtx()->GetDaughter(1));


        // Int_t track1ID = fPairVec_primary->at(prim_i).GetFirstDaughter();    // negative track
        // Int_t track2ID = fPairVec_primary->at(prim_i).GetSecondDaughter();   // positive track
        // Int_t track3ID = fPairVec_secondary->at(sec_i).GetFirstDaughter();   // negative track
        // Int_t track4ID = fPairVec_secondary->at(sec_i).GetSecondDaughter();  // positive track
        if (track1 == track2 || track1 == track3 || track1 == track4 || track2 == track3 || track2 == track4 || track3 == track4)  continue;

        // if (fdebug) {
        //   // int label1 = fMC->GetTrack(track1ID)->GetLabel();
        //   // int label2 = fMC->GetTrack(track2ID)->GetLabel();
        //   // int label3 = fMC->GetTrack(track3ID)->GetLabel();
        //   // int label4 = fMC->GetTrack(track4ID)->GetLabel();
        //   Int_t trackLabel1 = TMath::Abs(track1->GetLabel());                                  int track1PDG = fMC->GetTrack(trackLabel1)->PdgCode();
        //   Int_t trackLabel2 = TMath::Abs(track2->GetLabel());                                  int track2PDG = fMC->GetTrack(trackLabel2)->PdgCode();
        //   Int_t trackLabel3 = TMath::Abs(track3->GetLabel());                                  int track3PDG = fMC->GetTrack(trackLabel3)->PdgCode();
        //   Int_t trackLabel4 = TMath::Abs(track4->GetLabel());                                  int track4PDG = fMC->GetTrack(trackLabel4)->PdgCode();
        //   int motherLabel1 = TMath::Abs(fMC->GetTrack(trackLabel1)->GetMother());              int mother1PDG = fMC->GetTrack(motherLabel1)->PdgCode();
        //   int motherLabel2 = TMath::Abs(fMC->GetTrack(trackLabel2)->GetMother());              int mother2PDG = fMC->GetTrack(motherLabel2)->PdgCode();
        //   int motherLabel3 = TMath::Abs(fMC->GetTrack(trackLabel3)->GetMother());              int mother3PDG = fMC->GetTrack(motherLabel3)->PdgCode();
        //   int motherLabel4 = TMath::Abs(fMC->GetTrack(trackLabel4)->GetMother());              int mother4PDG = fMC->GetTrack(motherLabel4)->PdgCode();
        //   int grandMotherLabel1 = TMath::Abs(fMC->GetTrack(motherLabel1)->GetMother());        int grandMother1PDG = fMC->GetTrack(grandMotherLabel1)->PdgCode();
        //   int grandMotherLabel2 = TMath::Abs(fMC->GetTrack(motherLabel2)->GetMother());        int grandMother2PDG = fMC->GetTrack(grandMotherLabel2)->PdgCode();
        //   int grandMotherLabel3 = TMath::Abs(fMC->GetTrack(motherLabel3)->GetMother());        int grandMother3PDG = fMC->GetTrack(grandMotherLabel3)->PdgCode();
        //   int grandMotherLabel4 = TMath::Abs(fMC->GetTrack(motherLabel4)->GetMother());        int grandMother4PDG = fMC->GetTrack(grandMotherLabel4)->PdgCode();
        //   // int NDauthersGM1 =fMC->GetTrack(grandMotherLabel1)->GetNDaughters();
        //   // int NDauthersGM2 =fMC->GetTrack(grandMotherLabel2)->GetNDaughters();
        //   // int NDauthersGM3 =fMC->GetTrack(grandMotherLabel3)->GetNDaughters();
        //   // int NDauthersGM4 =fMC->GetTrack(grandMotherLabel4)->GetNDaughters();
        //
        //   if (grandMotherLabel1 == grandMotherLabel3 || grandMotherLabel2 == grandMotherLabel4) {
        //     std::cout << __LINE__ << " Decay from same GrandMother: " << " prim_i = " << prim_i << " sec_i = " << sec_i << std::endl;
        //     std::cout << "trackLabel1 = " << trackLabel1 << " track1PDG = " << track1PDG << " mother1Label = " << motherLabel1 << " mother1PDG = " << mother1PDG << " grandMotherLabel1 = "  << grandMotherLabel1 << " grandMother1PDG = " << grandMother1PDG << std::endl;
        //     std::cout << "trackLabel2 = " << trackLabel2 << " track2PDG = " << track2PDG << " mother2Label = " << motherLabel2 << " mother2PDG = " << mother2PDG << " grandMotherLabel2 = "  << grandMotherLabel2 << " grandMother2PDG = " << grandMother2PDG << std::endl;
        //     std::cout << "trackLabel3 = " << trackLabel3 << " track3PDG = " << track3PDG << " mother3Label = " << motherLabel3 << " mother3PDG = " << mother3PDG << " grandMotherLabel3 = "  << grandMotherLabel3 << " grandMother3PDG = " << grandMother3PDG << std::endl;
        //     std::cout << "trackLabel4 = " << trackLabel4 << " track4PDG = " << track4PDG << " mother4Label = " << motherLabel4 << " mother4PDG = " << mother4PDG << " grandMotherLabel4 = "  << grandMotherLabel4 << " grandMother4PDG = " << grandMother4PDG << std::endl;
        //   }
        // }

      }



                                                                                // if(fdebug) std::cout << __LINE__ << " DEBUG_AnalysisTask: COUT LINE" << std::endl;


      // Get LorentzVectors from reconstructed first pair (primary) and reconstructed V0 (secondary V0)
      TLorentzVector Lvec1;
      TLorentzVector Lvec2;
                                                                                // if(fdebug) std::cout << __LINE__ << " DEBUG_AnalysisTask: COUT LINE "  << std::endl;

      Lvec1.SetPtEtaPhiM(fPairVec_primary->at(prim_i).fPt , fPairVec_primary->at(prim_i).fEta , fPairVec_primary->at(prim_i).fPhi , fPairVec_primary->at(prim_i).fMass);
      Lvec2.SetPtEtaPhiM(fPairVec_secondary->at(sec_i).fPt, fPairVec_secondary->at(sec_i).fEta, fPairVec_secondary->at(sec_i).fPhi, fPairVec_secondary->at(sec_i).fMass);
                                                                                // if(fdebug) std::cout << __LINE__ << " DEBUG_AnalysisTask: COUT LINE " << std::endl;
      TLorentzVector LvecM = Lvec1 + Lvec2;
      double mass = LvecM.M();
      double pairpt = LvecM.Pt();
                                                                                // if(fdebug) std::cout << __LINE__ << " DEBUG_AnalysisTask: COUT LINE " << std::endl;

      if (fPairVec_primary != fPairVec_secondary) fVecHistPrefilters.at(0)->Fill(mass, pairpt/*, weight * centralityWeight*/);
      else if (fPairVec_primary == fPairVec_secondary) fVecHistPrefilters.at(2)->Fill(mass, pairpt/*, weight * centralityWeight*/);

      // if(fdebug) std::cout << __LINE__ << " DEBUG_AnalysisTask:Calculated PreFilter mass = " << mass << std::endl;

      if      ( (fPairVec_primary != fPairVec_secondary) && (mass > fUpperPrimSecPreFilterMass || mass < fLowerPrimSecPreFilterMass) ) {/*std::cout << __LINE__ << " Accept 4Pair with mass = " << mass << std::endl;*/
                                                                                                                                        fVecHistPrefilters.at(1)->Fill(mass, pairpt/*, weight * centralityWeight*/);
                                                                                                                                        continue;}
      else if ( (fPairVec_primary == fPairVec_secondary) && (mass > fUpperSecSecPreFilterMass || mass < fLowerSecSecPreFilterMass) )   {/*std::cout << __LINE__ << " Accept 4Pair with mass = " << mass << std::endl;*/
                                                                                                                                        fVecHistPrefilters.at(3)->Fill(mass, pairpt/*, weight * centralityWeight*/);
                                                                                                                                        continue;}
      else {
        AliVParticle* track1 = fEvent->GetTrack(fPairVec_primary->at(prim_i).GetFirstDaughter());   // primary eletrons, neg, first Daughter
        AliVParticle* track2 = fEvent->GetTrack(fPairVec_primary->at(prim_i).GetSecondDaughter());  // primary eletrons, pos, second Daughter

        Int_t labelD1 = -1;
        Int_t labelD2 = -1;
        labelD1 = track1->GetLabel();
        labelD2 = track2->GetLabel();
        fPreFilter_BadTracksLabel_primary.push_back(labelD1);
        fPreFilter_BadTracksLabel_primary.push_back(labelD2);


                                                                                // if(fdebug && (fPairVec_primary != fPairVec_secondary)) std::cout << __LINE__ << " DEBUG_AnalysisTask: PreFilter: erasing pair from primary and secondary vector. mass = " << mass << std::endl;
                                                                                // if(fdebug && (fPairVec_primary == fPairVec_secondary)) std::cout << __LINE__ << " DEBUG_AnalysisTask: PreFilter: erasing V0's from secondary vector. mass = " << mass << std::endl;
        fPairVec_primary->erase(fPairVec_primary->begin()+prim_i);
        fPairVec_secondary->erase(fPairVec_secondary->begin()+sec_i);
        prim_i--;
        sec_i--;
        break;
      }
    } // end of loop sec_i
  } // end of loop prim_i

  // remove tracks from neg and pos primary vector, if prim sec prefilter is applied
  if (fPairVec_primary != fPairVec_secondary) {
    for (unsigned int iNeg = 0; iNeg < fRecNegPart_primary.size(); iNeg++) {
        AliVParticle* track = fEvent->GetTrack(fRecNegPart_primary[iNeg].GetTrackID());
        Int_t tracklabel = track->GetLabel();
        if (std::find(fPreFilter_BadTracksLabel_primary.begin(), fPreFilter_BadTracksLabel_primary.end(), tracklabel) != fPreFilter_BadTracksLabel_primary.end()) fRecNegPart_primary.erase(fRecNegPart_primary.begin()+iNeg);
    }
    for (unsigned int iPos = 0; iPos < fRecPosPart_primary.size(); iPos++) {
        AliVParticle* track = fEvent->GetTrack(fRecPosPart_primary[iPos].GetTrackID());
        Int_t tracklabel = track->GetLabel();
        if (std::find(fPreFilter_BadTracksLabel_primary.begin(), fPreFilter_BadTracksLabel_primary.end(), tracklabel) != fPreFilter_BadTracksLabel_primary.end()) fRecPosPart_primary.erase(fRecPosPart_primary.begin()+iPos);
    }
  }
}


void AliAnalysisTaskEtaReconstruction::DoFourPairing(std::vector<TwoPair> fPairVec_primary, std::vector<TwoPair> fPairVec_secondary, Bool_t CaseRec, Bool_t CaseSmearing, double centralityWeight){
                                                                                // if(fdebug) std::cout << __LINE__ << " DEBUG_AnalysisTask: Line cout " << std::endl;
  								                                                              // if(fdebug) std::cout << __LINE__ << " DEBUG_AnalysisTask: fPairVec_primary = " << fPairVec_primary.size() << ", fPairVec_secondary = " << fPairVec_secondary.size() << ", CaseRec = " << CaseRec << ", CaseSmearing = " << CaseSmearing << std::endl;
                                                                                // if(fdebug) std::cout << " fPairVec_primary: "    << fPairVec_primary.size()   << std::endl;
  								                                                              // if(fdebug) std::cout << " fPairVec_secondary: "  << fPairVec_secondary.size() << std::endl;

  // TH1D* hFourPairOpeningAngle = new TH1D("Angle between primary and secondary pair", "Angle between primary and secondary pair ;rad ; #counts", 1000, 0. , 6.4);
  // hFourPairOpeningAngle->Sumw2();



  // std::cout << "Bool fDoFourPairing is " << fDoFourPairing << std::endl;
  for (unsigned int prim_i = 0; prim_i < fPairVec_primary.size(); ++prim_i){
    // if (fPairVec_primary[prim_i].fMCTwoSignal_acc_prim[0] == false) continue;
    // std::cout << "prim_i = " << prim_i << std::endl;
    for (unsigned int sec_i = 0; sec_i < fPairVec_secondary.size(); ++sec_i){
      // if (fPairVec_secondary[sec_i].fMCTwoSignal_acc_sec[0] == false) continue;
      // std::cout << "sec_i = " << sec_i << std::endl;

      std::vector<Bool_t> mcSignal_acc(fFourPairMCSignal.size()/2, kFALSE); // vector which stores if track is accepted by [i]-th mcsignal


      if (CaseRec == kFALSE) {
        AliVParticle* track1 = fMC->GetTrack(fPairVec_primary[prim_i].GetFirstDaughter());  // primary eletrons, neg, first Daughter
        AliVParticle* track2 = fMC->GetTrack(fPairVec_primary[prim_i].GetSecondDaughter());  // primary eletrons, pos, second Daughter
        AliVParticle* track3 = fMC->GetTrack(fPairVec_secondary[sec_i].GetFirstDaughter());  // secondary electron, neg, first Daughter
        AliVParticle* track4 = fMC->GetTrack(fPairVec_secondary[sec_i].GetSecondDaughter());  // secondary electron, pos, second Daughter

        // Check if electrons are from MCSignal Generator
              // if (!fGenPosPart_primary[prim_i].GetMCSignalPair() || !fGenNegPart_primary[prim_i].GetMCSignalPair() || !fGenNegPart_secondary[sec_i].GetMCSignalPair() || !fGenPosPart_secondary[sec_i].GetMCSignalPair()) continue;

        // Apply MC signals
        // Check MCSignal based on particles (not fully implemented! Check AliDielectronMC & AliDielectronSignalMC)
        for (unsigned int i = 0, j = 0 ; i < fFourPairMCSignal.size()/2 && j < fFourPairMCSignal.size(); ++i, j+=2){
          mcSignal_acc[i] = AliDielectronMC::Instance()->IsMCTruth(track1, track2, track3, track4, &(fFourPairMCSignal[j]), &(fFourPairMCSignal[j+1]));
        }


        // check if at least one mc signal is true
        if (CheckIfOneIsTrue(mcSignal_acc) == kFALSE) continue;
      }

      if (CaseRec == kTRUE) {
        AliVParticle* track1 = fEvent->GetTrack(fPairVec_primary[prim_i].GetFirstDaughter());  // primary eletrons, neg, first Daughter
        AliVParticle* track2 = fEvent->GetTrack(fPairVec_primary[prim_i].GetSecondDaughter());  // primary eletrons, pos, second Daughter

        // GetV0 daughters
        Int_t        iV0    = fPairVec_secondary[sec_i].GetV0ID(); //  i-th V0 in the Event
        AliAODv0*    fV0    = ((AliAODEvent*) fEvent)->GetV0(iV0); // select i-th V0 in the Event
        AliAODTrack* track3 = (AliAODTrack *) (fV0->GetSecondaryVtx()->GetDaughter(0));  // secondary electron, neg, first Daughter from V0
        AliAODTrack* track4 = (AliAODTrack *) (fV0->GetSecondaryVtx()->GetDaughter(1));  // secondary electron, pos, second Daughter from V0

        // Check if electrons are from MCSignal Generator
              // if (!fGenPosPart_primary[prim_i].GetMCSignalPair() || !fGenNegPart_primary[prim_i].GetMCSignalPair() || !fGenNegPart_secondary[sec_i].GetMCSignalPair() || !fGenPosPart_secondary[sec_i].GetMCSignalPair()) continue;
                                                                                // if(fdebug) std::cout << __LINE__ << " DEBUG_AnalysisTask: Line cout " << std::endl;

        // Apply MC signals
        // Check if it is according to mcsignals
        AliDielectronPair firstpair;
        AliDielectronPair secondpair;
        firstpair.SetKFUsage(false);
        secondpair.SetKFUsage(false);
        firstpair.SetTracks(static_cast<AliVTrack*>(track1), 11, static_cast<AliVTrack*>(track2), -11);
        secondpair.SetTracks(static_cast<AliVTrack*>(track3), 11, static_cast<AliVTrack*>(track4), -11);
        for (unsigned int i = 0, j = 0 ; i < fFourPairMCSignal.size()/2 && j < fFourPairMCSignal.size(); ++i, j+=2){
          mcSignal_acc[i] = AliDielectronMC::Instance()->IsMCTruth(&firstpair, &secondpair, &(fFourPairMCSignal[j]), &(fFourPairMCSignal[j+1]));
        }

        // check if at least one mc signal is true
        if (CheckIfOneIsTrue(mcSignal_acc) == kFALSE) continue;
      }
      // Construct pair variables from LorentzVectors
      TLorentzVector Lvec1;
      TLorentzVector Lvec2;
      Lvec1.SetPtEtaPhiM(fPairVec_primary[prim_i].fPt , fPairVec_primary[prim_i].fEta , fPairVec_primary[prim_i].fPhi , fPairVec_primary[prim_i].fMass);
      Lvec2.SetPtEtaPhiM(fPairVec_secondary[sec_i].fPt, fPairVec_secondary[sec_i].fEta, fPairVec_secondary[sec_i].fPhi, fPairVec_secondary[sec_i].fMass);


      TLorentzVector LvecM = Lvec1 + Lvec2;
      double mass = LvecM.M();
      double pairpt = LvecM.Pt();
      double weight = 1.;


      // if(fdebug) {
      //   if(mass<0.13){
      //     if (CaseRec == kFALSE && CaseSmearing == kFALSE) {
      //       AliVParticle* track1 = fMC->GetTrack(fPairVec_primary[prim_i].GetFirstDaughter());  // primary eletrons, neg, first Daughter
      //       AliVParticle* track2 = fMC->GetTrack(fPairVec_primary[prim_i].GetSecondDaughter());  // primary eletrons, pos, second Daughter
      //       AliVParticle* track3 = fMC->GetTrack(fPairVec_secondary[sec_i].GetFirstDaughter());  // secondary electron, neg, first Daughter
      //       AliVParticle* track4 = fMC->GetTrack(fPairVec_secondary[sec_i].GetSecondDaughter());  // secondary electron, pos, second Daughter
      //
      //       Int_t labelD1 = -1;   labelD1 = TMath::Abs(track1->GetLabel());     Int_t pdgD1 = track1->PdgCode();
      //       Int_t labelD2 = -1;   labelD2 = TMath::Abs(track2->GetLabel());     Int_t pdgD2 = track2->PdgCode();
      //       Int_t labelD3 = -1;   labelD3 = TMath::Abs(track3->GetLabel());     Int_t pdgD3 = track3->PdgCode();
      //       Int_t labelD4 = -1;   labelD4 = TMath::Abs(track4->GetLabel());     Int_t pdgD4 = track4->PdgCode();
      //
      //       Int_t labelM1 = -1;   labelM1 = TMath::Abs(track1->GetMother());    AliVParticle* mother1 = fMC->GetTrack(labelM1);       Int_t pdgM1 = mother1->PdgCode();
      //       Int_t labelM2 = -1;   labelM2 = TMath::Abs(track2->GetMother());    AliVParticle* mother2 = fMC->GetTrack(labelM2);       Int_t pdgM2 = mother2->PdgCode();
      //       Int_t labelM3 = -1;   labelM3 = TMath::Abs(track3->GetMother());    AliVParticle* mother3 = fMC->GetTrack(labelM3);       Int_t pdgM3 = mother3->PdgCode();
      //       Int_t labelM4 = -1;   labelM4 = TMath::Abs(track4->GetMother());    AliVParticle* mother4 = fMC->GetTrack(labelM4);       Int_t pdgM4 = mother4->PdgCode();
      //
      //       Int_t labelG1 = -1;   labelG1 = TMath::Abs(mother1->GetMother());   AliVParticle* grandmother1 = fMC->GetTrack(labelG1);     Int_t pdgG1 = grandmother1->PdgCode();       Int_t nDaughtersG1 =fMC->GetTrack(labelG1)->GetNDaughters();
      //       Int_t labelG2 = -1;   labelG2 = TMath::Abs(mother2->GetMother());   AliVParticle* grandmother2 = fMC->GetTrack(labelG2);     Int_t pdgG2 = grandmother2->PdgCode();       Int_t nDaughtersG2 =fMC->GetTrack(labelG2)->GetNDaughters();
      //       Int_t labelG3 = -1;   labelG3 = TMath::Abs(mother3->GetMother());   AliVParticle* grandmother3 = fMC->GetTrack(labelG3);     Int_t pdgG3 = grandmother3->PdgCode();       Int_t nDaughtersG3 =fMC->GetTrack(labelG3)->GetNDaughters();
      //       Int_t labelG4 = -1;   labelG4 = TMath::Abs(mother4->GetMother());   AliVParticle* grandmother4 = fMC->GetTrack(labelG4);     Int_t pdgG4 = grandmother4->PdgCode();       Int_t nDaughtersG4 =fMC->GetTrack(labelG4)->GetNDaughters();
      //
      //       std::cout << "DEBUG_AnalysisTask: Line cout " << std::endl;
      //       std::cout << "DEBUG_AnalysisTask: Mass = " << mass << ", pair_pT = " << pairpt << std::endl;
      //       std::cout << "   List Labels of the four paring particles and their mother as well as grand mother" << std::endl;
      //       std::cout << "    Label D1 = " << labelD1 << " pdgD1 = " << pdgD1 << std::endl;
      //       std::cout << "    Label D2 = " << labelD2 << " pdgD1 = " << pdgD2 << std::endl;
      //       std::cout << "    Label D3 = " << labelD3 << " pdgD1 = " << pdgD3 << std::endl;
      //       std::cout << "    Label D4 = " << labelD4 << " pdgD1 = " << pdgD4 << std::endl;
      //       std::cout << "    -----------" << std::endl;
      //       std::cout << "    Label M1 = " << labelM1 << " pdgM1 = " << pdgM1 << " Mass of particle " << mother1->M() << std::endl;
      //       std::cout << "    Label M2 = " << labelM2 << " pdgM1 = " << pdgM2 << " Mass of particle " << mother2->M() << std::endl;
      //       std::cout << "    Label M3 = " << labelM3 << " pdgM1 = " << pdgM3 << " Mass of particle " << mother3->M() << std::endl;
      //       std::cout << "    Label M4 = " << labelM4 << " pdgM1 = " << pdgM4 << " Mass of particle " << mother4->M() << std::endl;
      //       std::cout << "    -----------" << std::endl;
      //       std::cout << "    Label G1 = " << labelG1 << " pdgG1 = " << pdgG1 << " nDaughtersG1 = " << nDaughtersG1 << " Mass of particle " << grandmother1->M() << std::endl;
      //       std::cout << "    Label G2 = " << labelG2 << " pdgG1 = " << pdgG2 << " nDaughtersG2 = " << nDaughtersG2 << " Mass of particle " << grandmother2->M() << std::endl;
      //       std::cout << "    Label G3 = " << labelG3 << " pdgG1 = " << pdgG3 << " nDaughtersG3 = " << nDaughtersG3 << " Mass of particle " << grandmother3->M() << std::endl;
      //       std::cout << "    Label G4 = " << labelG4 << " pdgG1 = " << pdgG4 << " nDaughtersG4 = " << nDaughtersG4 << " Mass of particle " << grandmother4->M() << std::endl;
      //       std::cout << "    -----------" << std::endl;
      //       TLorentzVector fLVPart1; TLorentzVector fLVPart2; TLorentzVector fLVPart3; TLorentzVector fLVPart4; TLorentzVector fLVMother;
      //       fLVPart1.SetPtEtaPhiM(track1->Pt(), track1->Eta(), track1->Phi(), track1->M());
      //       fLVPart2.SetPtEtaPhiM(track2->Pt(), track2->Eta(), track2->Phi(), track2->M());
      //       fLVPart3.SetPtEtaPhiM(track3->Pt(), track3->Eta(), track3->Phi(), track3->M());
      //       fLVPart4.SetPtEtaPhiM(track4->Pt(), track4->Eta(), track4->Phi(), track4->M());
      //       fLVMother = fLVPart1+fLVPart2+fLVPart3+fLVPart4;
      //       std::cout << "   Recombine 4 particles to mother: mass = " << fLVMother.M() << std::endl;
      //     }
      //   }
      // }



      if (CaseRec == kTRUE) {
        double  pairAngle = Lvec1.Angle(Lvec2.Vect());
        for (unsigned int i = 0; i < fPairVec_secondary[sec_i].isReconstructed_secondary.size(); ++i){
          if (fPairVec_secondary[sec_i].isReconstructed_secondary[i] == kTRUE && fPairVec_primary[prim_i].isReconstructed_primary[i] == kTRUE){
            for (unsigned int j = 0; j < mcSignal_acc.size(); ++j){
              if (mcSignal_acc[j] == kTRUE){
                (dynamic_cast<TH1D *>(((TList*)((TList*)fFourPairCutListVec.at(i))->At(j))->At(0)))->Fill(pairAngle);
              }
            }
          }
        }
      //   hFourPairOpeningAngle->Fill(pairAngle);
      //   fOutputListSupportHistos->Add(hFourPairOpeningAngle);
      }

      if (CaseRec == kFALSE) {
        if (CaseSmearing == kFALSE) {
          for (unsigned int i = 0; i < mcSignal_acc.size(); ++i){
            if (mcSignal_acc[i] == kTRUE){
              fHistGenFourPair.at(i)->Fill(mass, pairpt, weight * centralityWeight);
            }
          } // end of loop over all MCsignals
        }
        if (CaseSmearing == kTRUE) {
          for (unsigned int i = 0; i < mcSignal_acc.size(); ++i){
            if (mcSignal_acc[i] == kTRUE){
              fHistGenSmearedFourPair.at(i)->Fill(mass, pairpt, weight * centralityWeight);
            }
          } // end of loop over all MCsignals
        }
      }
                                                                                // if(fdebug) std::cout << __LINE__ << " DEBUG_AnalysisTask: 4Pair with mass = " << mass << std::endl;
      if (CaseRec == kTRUE) {
        for (unsigned int i = 0; i < mcSignal_acc.size(); ++i){
          if (mcSignal_acc[i] == kTRUE){
            for (unsigned int j = 0; j < fPairVec_primary[prim_i].fFirstPartIsReconstructed.size(); ++j){
              if (fPairVec_primary[prim_i].fFirstPartIsReconstructed[j] == kTRUE && fPairVec_primary[prim_i].fSecondPartIsReconstructed[j] == kTRUE && fPairVec_secondary[sec_i].isReconstructed_secondary[j] == kTRUE ){
                fHistRecFourPair.at(j * mcSignal_acc.size() + i)->Fill(mass, pairpt, weight * centralityWeight);
              }
            }
          } // is selected by MCSignal
        } // end of loop over all MCsignals
      }
    } // end of loop sec_i
  } // end of loop prim_i
}
