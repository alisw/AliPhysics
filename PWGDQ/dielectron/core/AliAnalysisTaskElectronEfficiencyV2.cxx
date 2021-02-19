
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
//       Single Electron and Pair Efficiency Task                        //
//                                        (description in .h file)       //
//                                                                       //
///////////////////////////////////////////////////////////////////////////

#include "AliAnalysisTaskElectronEfficiencyV2.h"
#include "AliVTrack.h"
#include "AliAODInputHandler.h"

#include "AliAODEvent.h"
#include "AliESDEvent.h"
#include "AliMCEvent.h"
#include "AliAnalysisFilter.h"
#include "AliAnalysisUtils.h"

#include "AliGenEventHeader.h"
#include "AliGenCocktailEventHeader.h"
#include "AliGenDPMjetEventHeader.h"
#include "AliGenPythiaEventHeader.h"
#include "AliGenHijingEventHeader.h"

#include "AliDielectronMC.h"
#include "AliDielectronVarManager.h"
#include "AliDielectronSignalMC.h"
#include "AliDielectronPair.h"
#include "AliDielectronHistos.h"

#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TAxis.h"
#include "TTree.h"
#include "THnSparse.h"
#include "TLorentzVector.h"

#include "TChain.h"
#include "TSystem.h"

#include <iostream>

AliAnalysisTaskElectronEfficiencyV2::~AliAnalysisTaskElectronEfficiencyV2(){
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
  delete fPairList;
  delete fResolutionList;
}

// ############################################################################
// ############################################################################
AliAnalysisTaskElectronEfficiencyV2::AliAnalysisTaskElectronEfficiencyV2(): AliAnalysisTaskSE(),
  fEventFilter(0x0),
  fResoFile(0x0),
  fResoFilename(""),
  fResoFilenameFromAlien(""),
  fArrResoPt(0x0),
  fArrResoEta(0x0),
  fArrResoPhi_Pos(0x0),
  fArrResoPhi_Neg(0x0),
  fOutputList(0x0),
  fSingleElectronList(0x0),
  fPairList(0x0),
  fResolutionList(0x0),
  fPGen_DeltaP(0x0),
  fPGen_PrecOverPGen(0x0),
  fPtGen_DeltaPt(0x0),
  fPtGen_DeltaPtOverPtGen(0x0),
  fPtGen_PtRecOverPtGen(0x0),
  fPtGen_DeltaPt_wGenSmeared(0x0),
  fPtGen_DeltaPtOverPtGen_wGenSmeared(0x0),
  fPtGen_PtRecOverPtGen_wGenSmeared(0x0),
  fPGen_DeltaEta(0x0),
  fPtGen_DeltaEta(0x0),
  fPGen_DeltaTheta(0x0),
  fPGen_DeltaPhi_Ele(0x0),
  fPGen_DeltaPhi_Pos(0x0),
  fPtGen_DeltaPhi_Ele(0x0),
  fPtGen_DeltaPhi_Pos(0x0),
  fThetaGen_DeltaTheta(0x0),
  fPhiGen_DeltaPhi(0x0),
  fPtBins(),
  fEtaBins(),
  fPhiBins(),
  fThetaBins(),
  fResolutionDeltaPtBins(),
  fResolutionRelPtBins(),
  fResolutionEtaBins(),
  fResolutionPhiBins(),
  fResolutionThetaBins(),
  fMassBins(),
  fPairPtBins(),
  fPhiVBins(),
  fDoGenSmearing(false),
  fPtMin(0.),
  fPtMax(0.),
  fEtaMin(-99.),
  fEtaMax(99.),
  fPtMinGen(0.),
  fPtMaxGen(0.),
  fEtaMinGen(-99.),
  fEtaMaxGen(99.),
  fSingleLegMCSignal(),
  fPairMCSignal(),
  fDielectronPairNotFromSameMother(),
  fGeneratorName(""),
  fGeneratorMCSignalName(""),
  fGeneratorULSSignalName(""),
  fGeneratorHashs(),
  fGeneratorMCSignalHashs(),
  fGeneratorULSSignalHashs(),
  fCheckGenID(kFALSE),
  fGeneratorIndex(),
  fGeneratorMCSignalIndex(),
  fGeneratorULSSignalIndex(),
  fPIDResponse(0x0),
  fEvent(0x0),
  fMC(0x0),
  fTrack(0x0),
  isAOD(false),
  fSelectPhysics(false),
  fTriggerMask(0),
  fTrackCuts(),
  fUsedVars(0x0),
  fSupportMCSignal(0),
  fSupportCutsetting(0),
  fHistEvents(0x0),
  fHistEventStat(0x0),
  fHistCentralityRaw(0x0),
  fHistCentrality(0x0),
  fHistVertex(0x0),
  fHistVertexContibutors(0x0),
  fHistNTracks(0x0),
  fMinCentrality(0.),
  fMaxCentrality(100),
  fCentralityEst("V0M"),
  fCentralityFile(0x0),
  fCentralityFilename(""),
  fCentralityFilenameFromAlien(""),
  fHistCentralityCorrection(0x0),
  fNBinsCentralityCorr(0.),
  fEntriesCentralityCorr(0.),
  fOutputListSupportHistos(0x0),
  fHistGenPosPart(),
  fHistGenNegPart(),
  fHistGenSmearedPosPart(),
  fHistGenSmearedNegPart(),
  fHistRecPosPart(),
  fHistRecNegPart(),
  fHistGenPair(),
  fHistGenSmearedPair(),
  fHistRecPair(),
  fHistGenPair_ULSandLS(),
  fHistGenSmearedPair_ULSandLS(),
  fHistRecPair_ULSandLS(),
  fWriteLegsFromPair(false),
  fPtMinLegsFromPair(-99.),
  fPtMaxLegsFromPair(-99.),
  fEtaMinLegsFromPair(-99.),
  fEtaMaxLegsFromPair(-99.),
  fPhiMinLegsFromPair(-99.),
  fPhiMaxLegsFromPair(-99.),
  fOpAngleMinLegsFromPair(-99.),
  fOpAngleMaxLegsFromPair(-99.),
  fPtNBinsLegsFromPair(-99),
  fEtaNBinsLegsFromPair(-99),
  fPhiNBinsLegsFromPair(-99),
  fOpAngleNBinsLegsFromPair(-99),
  fTHnSparseGenSmearedLegsFromPair(),
  fTHnSparseRecLegsFromPair(),
  fDoFillPhiV(false),
  fApplyPhivCut(false),
  fMaxMee(-1),
  fMinPhiV(3.2),
  fDoPairing(false),
  fDoULSandLS(false),
  fDeactivateLS(false),
  fGenNegPart(),
  fGenPosPart(),
  fRecNegPart(),
  fRecPosPart(),
  fDoCocktailWeighting(false),
  fCocktailFilename(""),
  fCocktailFilenameFromAlien(""),
  fCocktailFile(0x0),
  fPtPion(0x0),
  fPtEta(0x0),
  fPtEtaPrime(0x0),
  fPtRho(0x0),
  fPtOmega(0x0),
  fPtPhi(0x0),
  fPtJPsi(0x0),
  fPostPIDCntrdCorrTPC(0x0),
  fPostPIDWdthCorrTPC(0x0),
  fPostPIDCntrdCorrITS(0x0),
  fPostPIDWdthCorrITS(0x0),
  fPostPIDCntrdCorrTOF(0x0),
  fPostPIDWdthCorrTOF(0x0)
{
// ROOT IO constructor , don ’t allocate memory here !
}


// ############################################################################
// ############################################################################
AliAnalysisTaskElectronEfficiencyV2::AliAnalysisTaskElectronEfficiencyV2(const char * name) : AliAnalysisTaskSE(name),
  fEventFilter(0x0),
  fResoFile(0x0),
  fResoFilename(""),
  fResoFilenameFromAlien(""),
  fArrResoPt(0x0),
  fArrResoEta(0x0),
  fArrResoPhi_Pos(0x0),
  fArrResoPhi_Neg(0x0),
  fOutputList(0x0),
  fSingleElectronList(0x0),
  fPairList(0x0),
  fResolutionList(0x0),
  fPGen_DeltaP(0x0),
  fPGen_PrecOverPGen(0x0),
  fPtGen_DeltaPt(0x0),
  fPtGen_DeltaPtOverPtGen(0x0),
  fPtGen_PtRecOverPtGen(0x0),
  fPtGen_DeltaPt_wGenSmeared(0x0),
  fPtGen_DeltaPtOverPtGen_wGenSmeared(0x0),
  fPtGen_PtRecOverPtGen_wGenSmeared(0x0),
  fPGen_DeltaEta(0x0),
  fPtGen_DeltaEta(0x0),
  fPGen_DeltaTheta(0x0),
  fPGen_DeltaPhi_Ele(0x0),
  fPGen_DeltaPhi_Pos(0x0),
  fPtGen_DeltaPhi_Ele(0x0),
  fPtGen_DeltaPhi_Pos(0x0),
  fThetaGen_DeltaTheta(0x0),
  fPhiGen_DeltaPhi(0x0),
  fPtBins(),
  fEtaBins(),
  fPhiBins(),
  fThetaBins(),
  fResolutionDeltaPtBins(),
  fResolutionRelPtBins(),
  fResolutionEtaBins(),
  fResolutionPhiBins(),
  fResolutionThetaBins(),
  fMassBins(),
  fPairPtBins(),
  fPhiVBins(),
  fDoGenSmearing(false),
  fPtMin(0.),
  fPtMax(0.),
  fEtaMin(-99.),
  fEtaMax(99.),
  fPtMinGen(0.),
  fPtMaxGen(0.),
  fEtaMinGen(-99.),
  fEtaMaxGen(99.),
  fSingleLegMCSignal(),
  fPairMCSignal(),
  fDielectronPairNotFromSameMother(),
  fGeneratorName(""),
  fGeneratorMCSignalName(""),
  fGeneratorULSSignalName(""),
  fGeneratorHashs(),
  fGeneratorMCSignalHashs(),
  fGeneratorULSSignalHashs(),
  fCheckGenID(kFALSE),
  fGeneratorIndex(),
  fGeneratorMCSignalIndex(),
  fGeneratorULSSignalIndex(),
  fPIDResponse(0x0),
  fEvent(0x0),
  fMC(0x0),
  fTrack(0x0),
  isAOD(false),
  fSelectPhysics(false),
  fTriggerMask(0),
  fTrackCuts(),
  fUsedVars(0x0),
  fSupportMCSignal(0),
  fSupportCutsetting(0),
  fHistEvents(0x0),
  fHistEventStat(0x0),
  fHistCentralityRaw(0x0),
  fHistCentrality(0x0),
  fHistVertex(0x0),
  fHistVertexContibutors(0x0),
  fHistNTracks(0x0),
  fMinCentrality(0.),
  fMaxCentrality(100),
  fCentralityEst("V0M"),
  fCentralityFile(0x0),
  fCentralityFilename(""),
  fCentralityFilenameFromAlien(""),
  fHistCentralityCorrection(0x0),
  fNBinsCentralityCorr(0.),
  fEntriesCentralityCorr(0.),
  fOutputListSupportHistos(0x0),
  fHistGenPosPart(),
  fHistGenNegPart(),
  fHistGenSmearedPosPart(),
  fHistGenSmearedNegPart(),
  fHistRecPosPart(),
  fHistRecNegPart(),
  fHistGenPair(),
  fHistGenSmearedPair(),
  fHistRecPair(),
  fHistGenPair_ULSandLS(),
  fHistGenSmearedPair_ULSandLS(),
  fHistRecPair_ULSandLS(),
  fWriteLegsFromPair(false),
  fPtMinLegsFromPair(-99.),
  fPtMaxLegsFromPair(-99.),
  fEtaMinLegsFromPair(-99.),
  fEtaMaxLegsFromPair(-99.),
  fPhiMinLegsFromPair(-99.),
  fPhiMaxLegsFromPair(-99.),
  fOpAngleMinLegsFromPair(-99.),
  fOpAngleMaxLegsFromPair(-99.),
  fPtNBinsLegsFromPair(-99),
  fEtaNBinsLegsFromPair(-99),
  fPhiNBinsLegsFromPair(-99),
  fOpAngleNBinsLegsFromPair(-99),
  fTHnSparseGenSmearedLegsFromPair(),
  fTHnSparseRecLegsFromPair(),
  fDoFillPhiV(false),
  fApplyPhivCut(false),
  fMaxMee(-1),
  fMinPhiV(3.2),
  fDoPairing(false),
  fDoULSandLS(false),
  fDeactivateLS(false),
  fGenNegPart(),
  fGenPosPart(),
  fRecNegPart(),
  fRecPosPart(),
  fDoCocktailWeighting(false),
  fCocktailFilename(""),
  fCocktailFilenameFromAlien(""),
  fCocktailFile(0x0),
  fPtPion(0x0),
  fPtEta(0x0),
  fPtEtaPrime(0x0),
  fPtRho(0x0),
  fPtOmega(0x0),
  fPtPhi(0x0),
  fPtJPsi(0x0),
  fPostPIDCntrdCorrTPC(0x0),
  fPostPIDWdthCorrTPC(0x0),
  fPostPIDCntrdCorrITS(0x0),
  fPostPIDWdthCorrITS(0x0),
  fPostPIDCntrdCorrTOF(0x0),
  fPostPIDWdthCorrTOF(0x0)
{
  DefineInput (0, TChain::Class());
  DefineOutput (1, TList::Class());



}


// ############################################################################
// ############################################################################
void AliAnalysisTaskElectronEfficiencyV2::Terminate(Option_t* option){
  // fHistEventStat->SetAxisRange(0., fHistEventStat->GetMaximum() * 1.1, "Y");
}


// ############################################################################
// ############################################################################
void AliAnalysisTaskElectronEfficiencyV2::UserCreateOutputObjects(){
  fUsedVars = new TBits(AliDielectronVarManager::kNMaxValues);
  fUsedVars->SetBitNumber(AliDielectronVarManager::kP, kTRUE);
  fUsedVars->SetBitNumber(AliDielectronVarManager::kPIn, kTRUE);
  fUsedVars->SetBitNumber(AliDielectronVarManager::kITSnSigmaEle, kTRUE);
  fUsedVars->SetBitNumber(AliDielectronVarManager::kTPCnSigmaEle, kTRUE);
  fUsedVars->SetBitNumber(AliDielectronVarManager::kTOFnSigmaEle, kTRUE);
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

    if (!fResoFile) {
      AliFatal(Form("Could not open file %s", fResoFilename.c_str()));
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
    if (!fCocktailFile) {
      AliFatal(Form("Could not open file %s", fCocktailFilename.c_str()));
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
    if (fCentralityFile == 0x0){
      std::cout << "Location in AliEN: " <<  fCentralityFilenameFromAlien << std::endl;
      gSystem->Exec(Form("alien_cp alien://%s .", fCentralityFilenameFromAlien.c_str()));
      std::cout << "Copy centrality weighting from Alien" << std::endl;
      fCentralityFile = TFile::Open(fCentralityFilename.c_str());
    }
    if (!fCentralityFile) {
      AliFatal(Form("Could not open file %s", fCentralityFilename.c_str()));
    }
    TList* list_temp = (TList*)fCentralityFile->Get("efficiency");
    fHistCentralityCorrection = (TH1F*) list_temp->FindObject("centrality");
    if (fHistCentralityCorrection == 0x0){
      AliError(Form("Could not extract centrality histogram from file %s", fCentralityFilename.c_str()));
    }
    else {

      if((fMinCentrality<0.)||(fMaxCentrality<0.)) {
	fEntriesCentralityCorr =  fHistCentralityCorrection->Integral();
      }
      else {
	TAxis *xaxis = fHistCentralityCorrection->GetXaxis();
	Int_t bina_mix = xaxis->FindBin(fMinCentrality);
	Int_t binb_mix = xaxis->FindBin(fMaxCentrality);
	Double_t lowedge_mix = xaxis->GetBinLowEdge(bina_mix);
	Double_t upedge_mix = xaxis->GetBinUpEdge(binb_mix);
	if(lowedge_mix > (fMinCentrality+0.0000000001)) bina_mix--;
	if(lowedge_mix < (fMinCentrality-0.0000000001)) bina_mix++;
	if(upedge_mix < (fMaxCentrality-0.0000000001)) binb_mix++;
	if(upedge_mix > (fMaxCentrality+0.0000000001)) binb_mix--;
	fEntriesCentralityCorr =  fHistCentralityCorrection->Integral(bina_mix,binb_mix);
      }
      fNBinsCentralityCorr = CalculateNbins();
      
      std::cout << "Centrality correction On in the range " << fMinCentrality << " " << fMaxCentrality << std::endl;
      std::cout << "nbins: " << fNBinsCentralityCorr << std::endl;
      std::cout << "entries: " << fEntriesCentralityCorr << std::endl;
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
  const int fNphivBins = fPhiVBins.size()-1;

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
    fHistCentralityRaw      = new TH1F("centralityRaw", "centralityRaw", 100, 0., 100.);
    fHistCentrality         = new TH1F("centrality", "centrality", 100, 0., 100.);
    fHistVertex             = new TH1F("zVertex", "zVertex", 300, -15.0, 15.0);
    fHistVertexContibutors  = new TH1F("vtxContributor", "vtxContributor",5000,-0.5,4999.5);
    fHistNTracks            = new TH1F("nTracks", "nTracks", 4000, 0., 40000.);
    fOutputList->Add(fHistEvents);
    fOutputList->Add(fHistEventStat);
    fOutputList->Add(fHistCentralityRaw);
    fOutputList->Add(fHistCentrality);
    fOutputList->Add(fHistVertex);
    // fOutputList->Add(fHistVertexContibutors);
    fOutputList->Add(fHistNTracks);


    // ######################################################
    // ##########  Single Electrons #########################
    // ######################################################
    fSingleElectronList = new TList();
    fSingleElectronList->SetOwner();
      fSingleElectronList->SetName("SingleElectrons");
      // Create List with generated particles
      TList* Generated = new TList();
      Generated->SetOwner();
      Generated->SetName("Generated");
      for (unsigned int i = 0; i < fSingleLegMCSignal.size(); ++i){
        TH3D* th3_tmp_pos = new TH3D(Form("Ngen_Pos_%s", fSingleLegMCSignal.at(i).GetName()),";p_{T};#eta;#varphi",fNptBins,fPtBins.data(),fNetaBins,fEtaBins.data(),fNphiBins,fPhiBins.data());
        th3_tmp_pos->Sumw2();
        fHistGenPosPart.push_back(th3_tmp_pos);
        Generated->Add(th3_tmp_pos);
        TH3D* th3_tmp_neg = new TH3D(Form("Ngen_Neg_%s", fSingleLegMCSignal.at(i).GetName()),";p_{T};#eta;#varphi",fNptBins,fPtBins.data(),fNetaBins,fEtaBins.data(),fNphiBins,fPhiBins.data());
        th3_tmp_neg->Sumw2();
        fHistGenNegPart.push_back(th3_tmp_neg);
        Generated->Add(th3_tmp_neg);
      }

      // Create List with generated+smeared particles
      TList* GeneratedSmeared = new TList();
      GeneratedSmeared->SetName("GeneratedSmeared");
      GeneratedSmeared->SetOwner();
      for (unsigned int i = 0; i < fSingleLegMCSignal.size(); ++i){
        TH3D* th3_tmp_pos = new TH3D(Form("Ngen_Pos_%s", fSingleLegMCSignal.at(i).GetName()),";p_{T};#eta;#varphi",fNptBins,fPtBins.data(),fNetaBins,fEtaBins.data(),fNphiBins,fPhiBins.data());
        th3_tmp_pos->Sumw2();
        fHistGenSmearedPosPart.push_back(th3_tmp_pos);
        GeneratedSmeared->Add(th3_tmp_pos);
        TH3D* th3_tmp_neg = new TH3D(Form("Ngen_Neg_%s", fSingleLegMCSignal.at(i).GetName()),";p_{T};#eta;#varphi",fNptBins,fPtBins.data(),fNetaBins,fEtaBins.data(),fNphiBins,fPhiBins.data());
        th3_tmp_neg->Sumw2();
        fHistGenSmearedNegPart.push_back(th3_tmp_neg);
        GeneratedSmeared->Add(th3_tmp_neg);
      }

      fSingleElectronList->Add(Generated);
      fSingleElectronList->Add(GeneratedSmeared);

      // Generated reconstructed lists for every cutsetting one list and every MCsignal 2 histograms with pos and neg charge
      for (unsigned int list_i = 0; list_i < fTrackCuts.size(); ++list_i){
        TList* list = new TList();
        list->SetName(fTrackCuts.at(list_i)->GetName());
        list->SetOwner();

        for (unsigned int i = 0; i < fSingleLegMCSignal.size(); ++i){
          TH3D* th3_tmp_pos = new TH3D(Form("Nrec_Pos_%s", fSingleLegMCSignal.at(i).GetName()),";p_{T};#eta;#varphi",fNptBins,fPtBins.data(),fNetaBins,fEtaBins.data(),fNphiBins,fPhiBins.data());
          th3_tmp_pos->Sumw2();
          th3_tmp_pos->SetDirectory(0x0);
          fHistRecPosPart.push_back(th3_tmp_pos);
          list->Add(th3_tmp_pos);
          TH3D* th3_tmp_neg = new TH3D(Form("Nrec_Neg_%s", fSingleLegMCSignal.at(i).GetName()),";p_{T};#eta;#varphi",fNptBins,fPtBins.data(),fNetaBins,fEtaBins.data(),fNphiBins,fPhiBins.data());
          th3_tmp_neg->Sumw2();
          th3_tmp_neg->SetDirectory(0x0);
          fHistRecNegPart.push_back(th3_tmp_neg);
          list->Add(th3_tmp_neg);

        }
        fSingleElectronList->Add(list);
      }
      if (fDoGenSmearing == kTRUE){
        for (unsigned int list_i = 0; list_i < fTrackCuts.size(); ++list_i){
          TList* list = new TList();
          std::string gen_smeared_name = fTrackCuts.at(list_i)->GetName();
          gen_smeared_name += "_gen_smeared";
          list->SetName(gen_smeared_name.c_str());
          list->SetOwner();

          for (unsigned int i = 0; i < fSingleLegMCSignal.size(); ++i){
            TH3D* th3_tmp_pos_gen_smeared = new TH3D(Form("Nrec_Pos_%s_gen_smeared", fSingleLegMCSignal.at(i).GetName()),";p_{T};#eta;#varphi",fNptBins,fPtBins.data(),fNetaBins,fEtaBins.data(),fNphiBins,fPhiBins.data());
            th3_tmp_pos_gen_smeared->Sumw2();
            th3_tmp_pos_gen_smeared->SetDirectory(0x0);
            fHistRecPosPart.push_back(th3_tmp_pos_gen_smeared);
            list->Add(th3_tmp_pos_gen_smeared);
            TH3D* th3_tmp_neg_gen_smeared = new TH3D(Form("Nrec_Neg_%s_gen_smeared", fSingleLegMCSignal.at(i).GetName()),";p_{T};#eta;#varphi",fNptBins,fPtBins.data(),fNetaBins,fEtaBins.data(),fNphiBins,fPhiBins.data());
            th3_tmp_neg_gen_smeared->Sumw2();
            th3_tmp_neg_gen_smeared->SetDirectory(0x0);
            fHistRecNegPart.push_back(th3_tmp_neg_gen_smeared);
            list->Add(th3_tmp_neg_gen_smeared);

          }
          fSingleElectronList->Add(list);
        }
      }


      // ######################################################
      // #####################  PAIRS #########################
      // ######################################################

      if (fDoPairing == true){
        fPairList = new TList();
        fPairList->SetName("Pairs");
        fPairList->SetOwner();

        TList* GeneratedPairs = new TList();
        GeneratedPairs->SetName("Generated");
        GeneratedPairs->SetOwner();
        for (unsigned int i = 0; i < fPairMCSignal.size(); ++i){
          TH2D* th2_tmp = new TH2D(Form("Ngen_%s", fPairMCSignal.at(i).GetName()),";m_{ee};p_{T,ee}",fNmassBins,fMassBins.data(),fNpairptBins,fPairPtBins.data());
          th2_tmp->Sumw2();
          fHistGenPair.push_back(th2_tmp);
          GeneratedPairs->Add(th2_tmp);


        }
        if (fDoULSandLS == true){
          for (unsigned int i = 0; i < fSingleLegMCSignal.size(); ++i){
            TH2D* th2_tmp_ULS = new TH2D(Form("Ngen_ULS_%s", fSingleLegMCSignal.at(i).GetName()),";m_{ee};p_{T,ee}",fNmassBins,fMassBins.data(),fNpairptBins,fPairPtBins.data());
            th2_tmp_ULS->Sumw2();
            fHistGenPair_ULSandLS.push_back(th2_tmp_ULS);
            GeneratedPairs->Add(th2_tmp_ULS);
            if (!fDeactivateLS) {
              TH2D* th2_tmp_LSpp = new TH2D(Form("Ngen_LSpp_%s", fSingleLegMCSignal.at(i).GetName()),";m_{ee};p_{T,ee}",fNmassBins,fMassBins.data(),fNpairptBins,fPairPtBins.data());
              th2_tmp_LSpp->Sumw2();
              fHistGenPair_ULSandLS.push_back(th2_tmp_LSpp);
              GeneratedPairs->Add(th2_tmp_LSpp);
              TH2D* th2_tmp_LSnn = new TH2D(Form("Ngen_LSnn_%s", fSingleLegMCSignal.at(i).GetName()),";m_{ee};p_{T,ee}",fNmassBins,fMassBins.data(),fNpairptBins,fPairPtBins.data());
              th2_tmp_LSnn->Sumw2();
              fHistGenPair_ULSandLS.push_back(th2_tmp_LSnn);
              GeneratedPairs->Add(th2_tmp_LSnn);
            }
          }
        }
        fPairList->Add(GeneratedPairs);

        TList* GeneratedSmearedPairs = new TList();
        GeneratedSmearedPairs->SetName("GeneratedSmeared");
        GeneratedSmearedPairs->SetOwner();
        for (unsigned int i = 0; i < fPairMCSignal.size(); ++i){
          TH2D* th2_tmp = new TH2D(Form("Ngen_%s", fPairMCSignal.at(i).GetName()),";m_{ee};p_{T,ee}",fNmassBins,fMassBins.data(),fNpairptBins,fPairPtBins.data());
          th2_tmp->Sumw2();
          fHistGenSmearedPair.push_back(th2_tmp);
          GeneratedSmearedPairs->Add(th2_tmp);

          if (fWriteLegsFromPair){

            THnSparseF* fTHnSparseGenSmearedLegsFromPair_tmp= new THnSparseF(Form("fTHnSparseGenSmearedLegsFromPair_%s", fPairMCSignal.at(i).GetName()),Form("fTHnSparseGenSmearedLegsFromPair_%s;p_{t,Pos};#eta_{Pos};#phi_{Pos};p_{t,Neg};#eta_{Neg};#phi_{Neg};opAngle", fPairMCSignal.at(i).GetName()), nDim, nBins, min, max);
            fTHnSparseGenSmearedLegsFromPair_tmp->GetAxis(0)->SetName("ptPos");
            fTHnSparseGenSmearedLegsFromPair_tmp->GetAxis(1)->SetName("etaPos");
            fTHnSparseGenSmearedLegsFromPair_tmp->GetAxis(2)->SetName("phiPos");
            fTHnSparseGenSmearedLegsFromPair_tmp->GetAxis(3)->SetName("ptNeg");
            fTHnSparseGenSmearedLegsFromPair_tmp->GetAxis(4)->SetName("etaNeg");
            fTHnSparseGenSmearedLegsFromPair_tmp->GetAxis(5)->SetName("phiNeg");
            fTHnSparseGenSmearedLegsFromPair_tmp->GetAxis(6)->SetName("opAngle");
            fTHnSparseGenSmearedLegsFromPair.push_back(fTHnSparseGenSmearedLegsFromPair_tmp);
            GeneratedSmearedPairs->Add(fTHnSparseGenSmearedLegsFromPair_tmp);
          }
        }
        if (fDoULSandLS == true){
          for (unsigned int i = 0; i < fSingleLegMCSignal.size(); ++i){
            TH2D* th2_tmp_ULS = new TH2D(Form("Ngen_ULS_%s", fSingleLegMCSignal.at(i).GetName()),";m_{ee};p_{T,ee}",fNmassBins,fMassBins.data(),fNpairptBins,fPairPtBins.data());
            th2_tmp_ULS->Sumw2();
            fHistGenSmearedPair_ULSandLS.push_back(th2_tmp_ULS);
            GeneratedSmearedPairs->Add(th2_tmp_ULS);
            if (!fDeactivateLS) {
              TH2D* th2_tmp_LSpp = new TH2D(Form("Ngen_LSpp_%s", fSingleLegMCSignal.at(i).GetName()),";m_{ee};p_{T,ee}",fNmassBins,fMassBins.data(),fNpairptBins,fPairPtBins.data());
              th2_tmp_LSpp->Sumw2();
              fHistGenSmearedPair_ULSandLS.push_back(th2_tmp_LSpp);
              GeneratedSmearedPairs->Add(th2_tmp_LSpp);
              TH2D* th2_tmp_LSnn = new TH2D(Form("Ngen_LSnn_%s", fSingleLegMCSignal.at(i).GetName()),";m_{ee};p_{T,ee}",fNmassBins,fMassBins.data(),fNpairptBins,fPairPtBins.data());
              th2_tmp_LSnn->Sumw2();
              fHistGenSmearedPair_ULSandLS.push_back(th2_tmp_LSnn);
              GeneratedSmearedPairs->Add(th2_tmp_LSnn);
            }
          }
        }
        fPairList->Add(GeneratedSmearedPairs);

        // Generated reconstructed lists for every cutsetting one list and every MCsignal 1 histogram
        for (unsigned int list_i = 0; list_i < fTrackCuts.size(); ++list_i){
          TList* list = new TList();
          list->SetName(fTrackCuts.at(list_i)->GetName());
          list->SetOwner();

          for (unsigned int i = 0; i < fPairMCSignal.size(); ++i){

            if(fDoFillPhiV){
              TH3D* th3_tmp_PhiV = new TH3D(Form("Nrec_%s_MPtPhiV", fPairMCSignal.at(i).GetName()),";m_{ee};p_{T,ee};#varphi_{V}",fNmassBins,fMassBins.data(),fNpairptBins,fPairPtBins.data(),fNphivBins,fPhiVBins.data());
              th3_tmp_PhiV->Sumw2();
              fHistRecPair.push_back(th3_tmp_PhiV);
              list->Add(th3_tmp_PhiV);
            }
            else{
              TH2D* th2_tmp = new TH2D(Form("Nrec_%s", fPairMCSignal.at(i).GetName()),";m_{ee};p_{T,ee}",fNmassBins,fMassBins.data(),fNpairptBins,fPairPtBins.data());
              th2_tmp->Sumw2();
              fHistRecPair.push_back(th2_tmp);
              list->Add(th2_tmp);
            }

            if (fWriteLegsFromPair){
              THnSparseF* fTHnSparseRecLegsFromPair_tmp= new THnSparseF(Form("fTHnSparseRecLegsFromPair_%s", fPairMCSignal.at(i).GetName()),Form("fTHnSparseRecLegsFromPair_%s;p_{t,Pos};#eta_{Pos};#phi_{Pos};p_{t,Neg};#eta_{Neg};#phi_{Neg};opAngle", fPairMCSignal.at(i).GetName()), nDim, nBins, min, max);
              fTHnSparseRecLegsFromPair_tmp->GetAxis(0)->SetName("ptPos");
              fTHnSparseRecLegsFromPair_tmp->GetAxis(1)->SetName("etaPos");
              fTHnSparseRecLegsFromPair_tmp->GetAxis(2)->SetName("phiPos");
              fTHnSparseRecLegsFromPair_tmp->GetAxis(3)->SetName("ptNeg");
              fTHnSparseRecLegsFromPair_tmp->GetAxis(4)->SetName("etaNeg");
              fTHnSparseRecLegsFromPair_tmp->GetAxis(5)->SetName("phiNeg");
              fTHnSparseRecLegsFromPair_tmp->GetAxis(6)->SetName("opAngle");
              fTHnSparseRecLegsFromPair_tmp->SetName(Form("fTHnSparseRecLegsFromPairTest_%s;ptPos;etaPos;phiPos;ptNeg;etaNeg;phiNeg;opAngle", fPairMCSignal.at(i).GetName()));
              fTHnSparseRecLegsFromPair_tmp->SetTitle(Form("fTHnSparseRecLegsFromPairTest_%s;ptPos;etaPos;phiPos;ptNeg;etaNeg;phiNeg;opAngle", fPairMCSignal.at(i).GetName()));
              fTHnSparseRecLegsFromPair.push_back(fTHnSparseRecLegsFromPair_tmp);
              list->Add(fTHnSparseRecLegsFromPair_tmp);
            }
          }
          if (fDoULSandLS == true){
            for (unsigned int i = 0; i < fSingleLegMCSignal.size(); ++i){
              TH2D* th2_tmp_ULS = new TH2D(Form("Nrec_ULS_%s", fSingleLegMCSignal.at(i).GetName()),";m_{ee};p_{T,ee}",fNmassBins,fMassBins.data(),fNpairptBins,fPairPtBins.data());
              th2_tmp_ULS->Sumw2();
              fHistRecPair_ULSandLS.push_back(th2_tmp_ULS);
              list->Add(th2_tmp_ULS);
              if (!fDeactivateLS) {
                TH2D* th2_tmp_LSpp = new TH2D(Form("Nrec_LSpp_%s", fSingleLegMCSignal.at(i).GetName()),";m_{ee};p_{T,ee}",fNmassBins,fMassBins.data(),fNpairptBins,fPairPtBins.data());
                th2_tmp_LSpp->Sumw2();
                fHistRecPair_ULSandLS.push_back(th2_tmp_LSpp);
                list->Add(th2_tmp_LSpp);
                TH2D* th2_tmp_LSnn = new TH2D(Form("Nrec_LSnn_%s", fSingleLegMCSignal.at(i).GetName()),";m_{ee};p_{T,ee}",fNmassBins,fMassBins.data(),fNpairptBins,fPairPtBins.data());
                th2_tmp_LSnn->Sumw2();
                fHistRecPair_ULSandLS.push_back(th2_tmp_LSnn);
                list->Add(th2_tmp_LSnn);
              }
            }
          }
          fPairList->Add(list);
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

    fOutputList->Add(fSingleElectronList);
    fOutputList->Add(fResolutionList);
    if (fDoPairing) fOutputList->Add(fPairList);

    CreateSupportHistos();
    fOutputList->Add(fOutputListSupportHistos);


  PostData(1, fOutputList);
}


// ############################################################################
// ############################################################################
void AliAnalysisTaskElectronEfficiencyV2::UserExec(Option_t* option){
  const double pi = TMath::Pi();

  // THIS MIGHT BE USED IN THE FUTURE TO SEPARATE DIFFERENT GENERATORS FROM EACH OTHER
  // AliMCEvent* mcEvent = MCEvent();
  // if(!mcEvent)return 0;
  //
  // for(int it = 0;it < mcEvent->GetNumberOfTracks();++it){
  //
  //   //if ESD AliMCParticle* part = (AliMCParticle*)mcEvent->GetTrack(it);
  //   //if AOD AliVParticle* part=(AliVParticle*)mcEvent->GetTrack(it);
  //


  fGenNegPart.clear();
  fGenPosPart.clear();
  fRecNegPart.clear();
  fRecPosPart.clear();

  // ##########################################################
  // Set MC event
  if(!AliDielectronMC::Instance()->ConnectMCEvent()) return;

  // ##########################################################
  // Manage AOD&ESD handling and the corresponding events

  isAOD = false;
  AliInputEventHandler *eventHandler = nullptr;
  AliInputEventHandler *eventHandlerMC = nullptr;

  if ((AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler())->IsA() == AliAODInputHandler::Class()){
    isAOD = true;
    eventHandler = dynamic_cast<AliAODInputHandler*> (AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
    eventHandlerMC = eventHandler;
  }
  else if ((AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler())->IsA() == AliESDInputHandler::Class()){
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


  // ##########################################################
  // print generator header
  if(isAOD){//for AOD
    AliAODMCHeader* mcHeader = (AliAODMCHeader*)fEvent->GetList()->FindObject(AliAODMCHeader::StdBranchName());
    TList *chl = (TList*)mcHeader->GetCocktailHeaders();//cocktail hearder list
    const Int_t Ngen = chl->GetEntries();
    AliInfo(Form("N generators = %d",Ngen));
    for(Int_t igen=0;igen<Ngen;igen++){
      AliGenEventHeader *gh = (AliGenEventHeader*)chl->At(igen);
      AliInfo(Form("Generator name = %s , NProduced = %d.",gh->GetName(),gh->NProduced()));
    }//end of generator loop
  }
  else{//for ESD
    AliGenEventHeader* genHeader = fMC->GenEventHeader();
    AliGenHijingEventHeader* hijingGenHeader = dynamic_cast<AliGenHijingEventHeader*>(genHeader);
    AliGenPythiaEventHeader* pythiaGenHeader = dynamic_cast<AliGenPythiaEventHeader*>(genHeader);
    AliGenDPMjetEventHeader* dpmjetGenHeader = dynamic_cast<AliGenDPMjetEventHeader*>(genHeader);

    if(hijingGenHeader == NULL && pythiaGenHeader == NULL && dpmjetGenHeader == NULL){
      AliGenCocktailEventHeader* genCocktailHeader = dynamic_cast<AliGenCocktailEventHeader*>(genHeader);
      TList *chl = (TList*)genCocktailHeader->GetHeaders();//cocktail header list
      const Int_t Ngen = chl->GetEntries();
      AliInfo(Form("N generators = %d",Ngen));
      for (Int_t igen=0; igen<Ngen; igen++) {
        AliGenEventHeader *gh = (AliGenEventHeader*)chl->At(igen);
        AliInfo(Form("Cocktail header is found : Generator name = %s , NProduced = %d.",gh->GetName(),gh->NProduced()));
      }
    }
    else if(hijingGenHeader) AliInfo(Form("Hijing header is found : Generator name = %s , NProduced = %d.",hijingGenHeader->GetName(),hijingGenHeader->NProduced()));
    else if(pythiaGenHeader) AliInfo(Form("Pythia header is found : Generator name = %s , NProduced = %d.",pythiaGenHeader->GetName(),pythiaGenHeader->NProduced()));
    else if(dpmjetGenHeader) AliInfo(Form("DPMjet header is found : Generator name = %s , NProduced = %d.",dpmjetGenHeader->GetName(),dpmjetGenHeader->NProduced()));
  }


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

  // ##########################################################
  // Apply centrality selection
  double centralityF = -1;
  AliMultSelection *multSelection = (AliMultSelection*)fEvent->FindListObject("MultSelection");
  if(multSelection){
    centralityF  = multSelection->GetMultiplicityPercentile(fCentralityEst, kFALSE);
  }
  if(centralityF == -1 || (fMaxCentrality == -1 && fMinCentrality == -1)){
    /*Centrality estimation failed or no requirement requested.*/
  }
  else if(centralityF > fMaxCentrality || centralityF < fMinCentrality) {
    return; // reject event
  }

  fHistEventStat->Fill(kCentralityEvents);
  fHistEvents->Fill(0.5);
  fHistCentralityRaw->Fill(centralityF);


  // Calculating the weight when centrality correction is applied
  double centralityWeight = 1.;
  if (fHistCentralityCorrection != 0x0){

    // find bin
    TAxis *xaxis = fHistCentralityCorrection->GetXaxis();
    Int_t bin_cent = xaxis->FindBin(centralityF);
    Double_t lowedge_cent = xaxis->GetBinLowEdge(bin_cent);
    Double_t upedge_cent = xaxis->GetBinUpEdge(bin_cent);
    if(lowedge_cent > (centralityF+0.0000000001)) bin_cent--;
    if(upedge_cent < (centralityF-0.0000000001)) bin_cent++;

    
    if((fNBinsCentralityCorr>0.) && (fHistCentralityCorrection->GetBinContent(bin_cent)>0.)) centralityWeight = fEntriesCentralityCorr/(fNBinsCentralityCorr*fHistCentralityCorrection->GetBinContent(bin_cent));
    
    //centralityWeight = (fHistCentralityCorrection->GetEntries() / fHistCentralityCorrection->GetNbinsX()) / fHistCentralityCorrection->FindBin(centralityF) ;
    //std::cout << "cent: " << centralityF << "  " << "weight: " << centralityWeight << std::endl;
  }
  fHistCentrality->Fill(centralityF,centralityWeight);

  // ##########################################################
  // Fill Multiplicity histogram
  int nTracks = fEvent->GetNumberOfTracks();
  fHistNTracks->Fill(nTracks,centralityWeight);


  // ######################################################
  // ######################################################
  // ######################################################
  // Start particle loop
  for(int iPart = 0; iPart < fMC->GetNumberOfTracks(); iPart++) {
    AliVParticle* mcPart1  = (AliVParticle*)fMC->GetTrack(iPart);
    if (!mcPart1) continue;

    // ##########################################################
    // Checking minimum and maximum values for generated particles
    if (mcPart1->Pt()  < fPtMinGen  || mcPart1->Pt()  > fPtMaxGen)  continue;
    if (mcPart1->Eta() < fEtaMinGen || mcPart1->Eta() > fEtaMaxGen) continue;

    // ##########################################################
    // Check MC signals
    std::vector<Bool_t> mcSignal_acc(fSingleLegMCSignal.size(), kFALSE); // initialize vector which stores if track is accepted by [i]-th mcsignal
    CheckSingleLegMCsignals(mcSignal_acc, iPart);

    // ##########################################################
    // check if at least one mc signal is true
    if (CheckIfOneIsTrue(mcSignal_acc) == kFALSE) continue;

    // ##########################################################
    // check if correct generator used
    bool generatorForMCSignal  = CheckGenerator(iPart, fGeneratorMCSignalHashs);
    bool generatorForULSSignal = CheckGenerator(iPart, fGeneratorULSSignalHashs);

    if(fCheckGenID){
      generatorForMCSignal  = CheckGeneratorIndex(iPart, fGeneratorMCSignalIndex);
      generatorForULSSignal = CheckGeneratorIndex(iPart, fGeneratorULSSignalIndex);
    }

    if (!generatorForMCSignal && !generatorForULSSignal) continue;
    // if (!CheckGenerator(iPart, fGeneratorHashs)) continue;


    // ##########################################################
    // Creating particles to summarize all the data
    int motherID = TMath::Abs(mcPart1->GetMother());
    Particle part = CreateParticle(mcPart1);
    part.isMCSignal = mcSignal_acc;
    part.SetTrackID(iPart);
    part.SetMotherID(motherID);
    part.SetULSSignalPair(generatorForULSSignal);
    part.SetMCSignalPair(generatorForMCSignal);

    // ##########################################################
    // check if electron comes from a mother with ele+pos as daughters
    CheckIfFromMotherWithDielectronAsDaughter(part);

    // ##########################################################
    // Filling generated particle histograms according to MCSignals
    for (unsigned int i = 0; i < part.isMCSignal.size(); ++i){
      if (part.isMCSignal[i]) {
        if      (part.fCharge < 0){
          dynamic_cast<TH3D*>(fHistGenNegPart.at(i))->Fill(part.fPt, part.fEta, part.fPhi, centralityWeight);
        }
        else if (part.fCharge > 0) {
          dynamic_cast<TH3D*>(fHistGenPosPart.at(i))->Fill(part.fPt, part.fEta, part.fPhi, centralityWeight);
        }
      }
    }

    // ##########################################################
    // Filling generated+smeared particle histograms according to MCSignals
    // and separated into pos and neg charge
    if (fArrResoPt){ // Smear particles to fill "GeneratedSmeared"
      TLorentzVector smearedVec = ApplyResolution(part.fPt, part.fEta, part.fPhi, part.fCharge);
      part.fPt_smeared  = smearedVec.Pt();
      part.fEta_smeared = smearedVec.Eta();
      if (smearedVec.Phi() < 0) part.fPhi_smeared = smearedVec.Phi()+ 2 * pi;
      else part.fPhi_smeared = smearedVec.Phi();

      for (unsigned int i = 0; i < part.isMCSignal.size(); ++i){
        if (part.isMCSignal[i]) {
          if      (part.fCharge < 0){
            dynamic_cast<TH3D*>(fHistGenSmearedNegPart.at(i))->Fill(part.fPt_smeared, part.fEta_smeared, part.fPhi_smeared , centralityWeight);
          }
          else if (part.fCharge > 0) {
            dynamic_cast<TH3D*>(fHistGenSmearedPosPart.at(i))->Fill(part.fPt_smeared, part.fEta_smeared, part.fPhi_smeared , centralityWeight);
          }
        }
      }
    }

    if      (fDoPairing == true && part.fCharge < 0) fGenNegPart.push_back(part); // store particles for later pairing
    else if (fDoPairing == true && part.fCharge > 0) fGenPosPart.push_back(part); // store particles for later pairing

  }// end of MC track loop


  // ##########################################################
  // ##########################################################
  // ##########################################################
  // Start reconstructed track Loop

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
    std::vector<Bool_t> mcSignal_acc(fSingleLegMCSignal.size(), kFALSE); // vector which stores if track is accepted by [i]-th mcsignal
    CheckSingleLegMCsignals(mcSignal_acc, abslabel);

    // ##########################################################
    // check if at least one mc signal is true otherwise skip this particle
    if (CheckIfOneIsTrue(mcSignal_acc) == kFALSE) continue;

    // ##########################################################
    // check if correct generator used
    bool generatorForMCSignal  = CheckGenerator(label, fGeneratorMCSignalHashs);
    bool generatorForULSSignal = CheckGenerator(label, fGeneratorULSSignalHashs);

    if(fCheckGenID){
      generatorForMCSignal  = CheckGeneratorIndex(label, fGeneratorMCSignalIndex);
      generatorForULSSignal = CheckGeneratorIndex(label, fGeneratorULSSignalIndex);
    }

    // std::cout << "generatorForMCSignal = " << generatorForMCSignal << std::endl;
    // std::cout << "generatorForULSSignal = " << generatorForULSSignal << std::endl;
    if (!generatorForMCSignal && !generatorForULSSignal) continue;
    // if (!CheckGenerator(label, fGeneratorHashs)) continue;

    // ##########################################################
    // Check if particle is passing selection cuts
    std::vector<bool> selected(fTrackCuts.size(), kFALSE); // vector which stores if track is accepted by [i]-th selection cut
    for (UInt_t iCut=0; iCut<fTrackCuts.size(); ++iCut){ // loop over all specified cutInstances
      UInt_t selectedMask=( 1 << fTrackCuts.at(iCut)->GetCuts()->GetEntries())-1;
      // cutting logic taken from AliDielectron::FillTrackArrays()
      // apply track cuts
      UInt_t cutMask = fTrackCuts.at(iCut)->IsSelected(track);
      if (cutMask == selectedMask) {selected[iCut] = kTRUE; /*std::cout << "reconstructed TRUE" << std::endl;*/}
    }

    // ##########################################################
    // check if at least one is selected by cuts otherwise skip this particle
    if (CheckIfOneIsTrue(selected) == kFALSE) continue;

    // ##########################################################
    // Create summary particle from track info
    int motherID = TMath::Abs(fMC->GetTrack(abslabel)->GetMother());
    Particle part  = CreateParticle(track);
    part.isMCSignal = mcSignal_acc;
    part.isReconstructed = selected;
    part.SetTrackID(iTracks);
    part.SetMotherID(motherID);
    part.SetULSSignalPair(generatorForULSSignal);
    part.SetMCSignalPair(generatorForMCSignal);
    // ##########################################################
    // check if electron comes from a mother with ele+pos as daughters
    CheckIfFromMotherWithDielectronAsDaughter(part);

    // ##########################################################
    if      (fDoPairing == true && part.fCharge < 0) fRecNegPart.push_back(part);
    else if (fDoPairing == true && part.fCharge > 0) fRecPosPart.push_back(part);

    if (fDoGenSmearing == true && fArrResoPt){
      AliVParticle* genTrack = fMC->GetTrack(abslabel);
      double pt_temp  = genTrack->Pt();
      double eta_temp = genTrack->Eta();
      double phi_temp = genTrack->Phi();

      TLorentzVector smearedVec = ApplyResolution(pt_temp, eta_temp, phi_temp, part.fCharge);
      part.fPt_smeared  = smearedVec.Pt();
      part.fEta_smeared = smearedVec.Eta();
      // part.fPhi_smeared = smearedVec.Phi();
      if (smearedVec.Phi() < 0) part.fPhi_smeared = smearedVec.Phi() + 2 * pi;
      else part.fPhi_smeared = smearedVec.Phi();

      // std::cout << "pt_rec: " << part.fPt << "  pt_gen: " << pt_temp << "  phi_rec:" << part.fPhi << "  phi_gen:" << phi_temp << "  phi_gen_smeared:" << part.fPhi_smeared << "  pdgCode:" << genTrack->PdgCode() << std::endl;
    }


    // ##########################################################
    // Filling generated particle histograms according to MCSignals
    for (unsigned int i = 0; i < part.isMCSignal.size(); ++i){
      for (unsigned int j = 0; j < part.isReconstructed.size(); ++j){
        if (part.isMCSignal[i] == kTRUE) {
          if (part.isReconstructed[j] == kTRUE){
            if      (part.fCharge < 0) {
              dynamic_cast<TH3D*>(fHistRecNegPart.at(j * part.isMCSignal.size() + i))->Fill(part.fPt, part.fEta, part.fPhi, centralityWeight);
            }
            else if (part.fCharge > 0) {
              dynamic_cast<TH3D*>(fHistRecPosPart.at(j * part.isMCSignal.size() + i))->Fill(part.fPt, part.fEta, part.fPhi, centralityWeight);
            }
            if (fDoGenSmearing == true && fArrResoPt){
              if      (part.fCharge < 0) {
                dynamic_cast<TH3D*>(fHistRecNegPart.at(j * part.isMCSignal.size() + i + part.isMCSignal.size() * part.isReconstructed.size()))->Fill(part.fPt_smeared, part.fEta_smeared, part.fPhi_smeared, centralityWeight);
              }
              else if (part.fCharge > 0) {
                dynamic_cast<TH3D*>(fHistRecPosPart.at(j * part.isMCSignal.size() + i + part.isMCSignal.size() * part.isReconstructed.size()))->Fill(part.fPt_smeared, part.fEta_smeared, part.fPhi_smeared, centralityWeight);
              }

            }

          }// is selected by cutsetting
        } // is selected by MC signal
      } // end of loop over all cutsettings
    } // end of loop over all MCsignals

    // ##########################################################
    // Fill support histograms with first cutsetting and first mcsignal
    if(part.isMCSignal[fSupportMCSignal] == true && part.isReconstructed[fSupportCutsetting] == kTRUE){

      AliVParticle* mcPart1 = fMC->GetTrack(abslabel);

      FillTrackHistograms(track, mcPart1); // Fill support histograms

      // ##########################################################
      // Fill resolution histograms
      // ##########################################################
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

        if (fDoGenSmearing == true && fArrResoPt){
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


  // ##########################################################
  // ##########################################################
  // ##########################################################
  // DO PAIRING
  // ##########################################################

  float ptPos       = -999;
  float etaPos      = -999;
  float phiPos      = -999;
  float ptNeg       = -999;
  float etaNeg      = -999;
  float phiNeg      = -999;
  float op_angle    = -999;

  if (fDoPairing){
    for (unsigned int neg_i = 0; neg_i < fGenNegPart.size(); ++neg_i){
      for (unsigned int pos_i = 0; pos_i < fGenPosPart.size(); ++pos_i){
        AliVParticle* mcPart1 = fMC->GetTrack(fGenNegPart[neg_i].GetTrackID());
        AliVParticle* mcPart2 = fMC->GetTrack(fGenPosPart[pos_i].GetTrackID());

        if (fDoULSandLS && fGenNegPart[neg_i].GetULSSignalPair() && fGenPosPart[pos_i].GetULSSignalPair()){
          // Calculates for single leg signals unlike sign

          // checks if unsmeared and smeared survive acceptance cuts
          bool selectedByKinematicCuts = true;
          if (fGenNegPart[neg_i].fPt < fPtMin || fGenNegPart[neg_i].fPt > fPtMax || fGenNegPart[neg_i].fEta < fEtaMin || fGenNegPart[neg_i].fEta > fEtaMax) selectedByKinematicCuts = false;
          if (fGenPosPart[pos_i].fPt < fPtMin || fGenPosPart[pos_i].fPt > fPtMax || fGenPosPart[pos_i].fEta < fEtaMin || fGenPosPart[pos_i].fEta > fEtaMax) selectedByKinematicCuts = false;
          bool selectedByKinematicCuts_smeared = true;
          if (fGenNegPart[neg_i].fPt_smeared < fPtMin || fGenNegPart[neg_i].fPt_smeared > fPtMax || fGenNegPart[neg_i].fEta_smeared < fEtaMin || fGenNegPart[neg_i].fEta_smeared > fEtaMax) selectedByKinematicCuts_smeared = false;
          if (fGenPosPart[pos_i].fPt_smeared < fPtMin || fGenPosPart[pos_i].fPt_smeared > fPtMax || fGenPosPart[pos_i].fEta_smeared < fEtaMin || fGenPosPart[pos_i].fEta_smeared > fEtaMax) selectedByKinematicCuts_smeared = false;

          if (selectedByKinematicCuts){
            TLorentzVector Lvec1;
            TLorentzVector Lvec2;
            Lvec1.SetPtEtaPhiM(fGenNegPart[neg_i].fPt, fGenNegPart[neg_i].fEta, fGenNegPart[neg_i].fPhi, AliPID::ParticleMass(AliPID::kElectron));
            Lvec2.SetPtEtaPhiM(fGenPosPart[pos_i].fPt, fGenPosPart[pos_i].fEta, fGenPosPart[pos_i].fPhi, AliPID::ParticleMass(AliPID::kElectron));
            TLorentzVector LvecM = Lvec1 + Lvec2;
            double mass = LvecM.M();
            double pairpt = LvecM.Pt();
            double weight = 1;
            for (unsigned int iMCSignal = 0; iMCSignal < fGenNegPart[neg_i].isMCSignal.size(); ++iMCSignal){
              // if (fGenNegPart[neg_i].isMCSignal[iMCSignal] == true && fGenPosPart[pos_i].isMCSignal[iMCSignal] == true){
              if (fGenNegPart[neg_i].isMCSignal[iMCSignal] == true && fGenPosPart[pos_i].isMCSignal[iMCSignal] == true &&
                  fGenNegPart[neg_i].DielectronPairFromSameMother[iMCSignal] == false && fGenPosPart[pos_i].DielectronPairFromSameMother[iMCSignal] == false){
                if (!fDeactivateLS) {
//                  std::cout << "Deactivate" << std::endl;
                 fHistGenPair_ULSandLS.at(3*iMCSignal)->Fill(mass, pairpt, weight);
                }
                else {
                  fHistGenPair_ULSandLS.at(1*iMCSignal)->Fill(mass, pairpt, weight);
                }
              }
            }
          }
          if (selectedByKinematicCuts_smeared){
            TLorentzVector Lvec1;
            TLorentzVector Lvec2;
            Lvec1.SetPtEtaPhiM(fGenNegPart[neg_i].fPt_smeared, fGenNegPart[neg_i].fEta_smeared, fGenNegPart[neg_i].fPhi_smeared, AliPID::ParticleMass(AliPID::kElectron));
            Lvec2.SetPtEtaPhiM(fGenPosPart[pos_i].fPt_smeared, fGenPosPart[pos_i].fEta_smeared, fGenPosPart[pos_i].fPhi_smeared, AliPID::ParticleMass(AliPID::kElectron));
            TLorentzVector LvecM = Lvec1 + Lvec2;
            double mass = LvecM.M();
            double pairpt = LvecM.Pt();
            double weight = 1;

            for (unsigned int iMCSignal = 0; iMCSignal < fGenNegPart[neg_i].isMCSignal.size(); ++iMCSignal){
              // if (fGenNegPart[neg_i].isMCSignal[iMCSignal] == true && fGenPosPart[pos_i].isMCSignal[iMCSignal] == true)
              if (fGenNegPart[neg_i].isMCSignal[iMCSignal] == true && fGenPosPart[pos_i].isMCSignal[iMCSignal] == true &&
                  fGenNegPart[neg_i].DielectronPairFromSameMother[iMCSignal] == false && fGenPosPart[pos_i].DielectronPairFromSameMother[iMCSignal] == false){
                if (!fDeactivateLS) {
                 fHistGenSmearedPair_ULSandLS.at(3*iMCSignal)->Fill(mass, pairpt, weight);
                }
                else {
                  fHistGenSmearedPair_ULSandLS.at(1*iMCSignal)->Fill(mass, pairpt, weight);
                }
              }
            }
          }
        } // End of ULS

        // Check if electrons are from MCSignal Generator
        if (!fGenPosPart[pos_i].GetMCSignalPair() || !fGenNegPart[neg_i].GetMCSignalPair()) continue;
        // std::cout << "fGenPosPart[pos_i].GetMCSignalPair() = " << fGenPosPart[pos_i].GetMCSignalPair() << std::endl;
        // std::cout << "fGenNegPart[neg_i].GetMCSignalPair() = " << fGenNegPart[neg_i].GetMCSignalPair() << std::endl;
        // std::cout << "#########" << std::endl;
        // Apply MC signals
        std::vector<Bool_t> mcSignal_acc(fPairMCSignal.size(), kFALSE); // vector which stores if track is accepted by [i]-th mcsignal

        // Check if it according to mcsignals
        for (unsigned int i = 0; i < fPairMCSignal.size(); ++i){
          mcSignal_acc[i] = AliDielectronMC::Instance()->IsMCTruth(mcPart1, mcPart2, &(fPairMCSignal[i]));
        }

        // check if at least one mc signal is true
        if (CheckIfOneIsTrue(mcSignal_acc) == kFALSE) continue;

        // This if clause is needed here because smearing can potentially smear tracks into the selected kinematic region
        bool selectedByKinematicCuts = true;
        if (fGenNegPart[neg_i].fPt < fPtMin || fGenNegPart[neg_i].fPt > fPtMax || fGenNegPart[neg_i].fEta < fEtaMin || fGenNegPart[neg_i].fEta > fEtaMax) selectedByKinematicCuts = false;
        if (fGenPosPart[pos_i].fPt < fPtMin || fGenPosPart[pos_i].fPt > fPtMax || fGenPosPart[pos_i].fEta < fEtaMin || fGenPosPart[pos_i].fEta > fEtaMax) selectedByKinematicCuts = false;

        if (selectedByKinematicCuts == true){
          // Construct pair variables from LorentzVectors
          TLorentzVector Lvec1;
          TLorentzVector Lvec2;
          Lvec1.SetPtEtaPhiM(mcPart1->Pt(), mcPart1->Eta(), mcPart1->Phi(), AliPID::ParticleMass(AliPID::kElectron));
          Lvec2.SetPtEtaPhiM(mcPart2->Pt(), mcPart2->Eta(), mcPart2->Phi(), AliPID::ParticleMass(AliPID::kElectron));
          TLorentzVector LvecM = Lvec1 + Lvec2;
          double mass = LvecM.M();
          double pairpt = LvecM.Pt();
          double weight = 1.;
          if (fCocktailFile) {
            if (fGenNegPart[neg_i].GetMotherID() == fGenPosPart[pos_i].GetMotherID()){
              double motherpt = fMC->GetTrack(fGenNegPart[neg_i].GetMotherID())->Pt();
              weight *= GetWeight(fGenNegPart[neg_i], fGenPosPart[pos_i], motherpt);
            }
            else{
              weight = 0; // if should not fail by definition. but does in 13 / 10000000 cases
            }
          }

          for (unsigned int i = 0; i < mcSignal_acc.size(); ++i){
            if (mcSignal_acc[i] == kTRUE){
              fHistGenPair.at(i)->Fill(mass, pairpt, weight * centralityWeight);
            }
          } // end of loop over all MCsignals
        }
        if (fArrResoPt){ // Smear particles to fill "GeneratedSmeared"

          if (fGenNegPart[neg_i].fPt_smeared < fPtMin || fGenNegPart[neg_i].fPt_smeared > fPtMax || fGenNegPart[neg_i].fEta_smeared < fEtaMin || fGenNegPart[neg_i].fEta_smeared > fEtaMax) continue;
          if (fGenPosPart[pos_i].fPt_smeared < fPtMin || fGenPosPart[pos_i].fPt_smeared > fPtMax || fGenPosPart[pos_i].fEta_smeared < fEtaMin || fGenPosPart[pos_i].fEta_smeared > fEtaMax) continue;

          // Construct pair variables from LorentzVectors
          TLorentzVector Lvec1_smeared;
          TLorentzVector Lvec2_smeared;
          Lvec1_smeared.SetPtEtaPhiM(fGenNegPart[neg_i].fPt_smeared, fGenNegPart[neg_i].fEta_smeared, fGenNegPart[neg_i].fPhi_smeared, AliPID::ParticleMass(AliPID::kElectron));
          Lvec2_smeared.SetPtEtaPhiM(fGenPosPart[pos_i].fPt_smeared, fGenPosPart[pos_i].fEta_smeared, fGenPosPart[pos_i].fPhi_smeared, AliPID::ParticleMass(AliPID::kElectron));
          TLorentzVector LvecM_smeared = Lvec1_smeared + Lvec2_smeared;
          double massSmeared = LvecM_smeared.M();
          double pairptSmeared = LvecM_smeared.Pt();
          double weight = 1.;
          if (fCocktailFile) {
            if (fGenNegPart[neg_i].GetMotherID() == fGenPosPart[pos_i].GetMotherID()){
              double motherpt = fMC->GetTrack(fGenNegPart[neg_i].GetMotherID())->Pt();
              weight *= GetWeight(fGenNegPart[neg_i], fGenPosPart[pos_i], motherpt);
            }
            else{
              weight = 0; // if should not fail by definition. but does in 13 / 10000000 cases
            }
          }

          for (unsigned int i = 0; i <  mcSignal_acc.size(); ++i){
            if (mcSignal_acc[i] == kTRUE){
              fHistGenSmearedPair.at(i)->Fill(massSmeared, pairptSmeared, weight * centralityWeight);
              if (fWriteLegsFromPair){
                ptNeg  = fGenNegPart[neg_i].fPt_smeared;
                etaNeg = fGenNegPart[neg_i].fEta_smeared;
                phiNeg = fGenNegPart[neg_i].fPhi_smeared;
                ptPos  = fGenPosPart[pos_i].fPt_smeared;
                etaPos = fGenPosPart[pos_i].fEta_smeared;
                phiPos = fGenPosPart[pos_i].fPhi_smeared;
                op_angle = Lvec2_smeared.Angle(Lvec1_smeared.Vect());

                double tuple[7] = {ptPos,etaPos,phiPos,ptNeg,etaNeg,phiNeg,op_angle};
                fTHnSparseGenSmearedLegsFromPair[i]->Fill(tuple);
              }
            }
          } // end of loop over all MCsignals
        } // end of smearing
      } // end of loop over all positive particles
    } // end of loop over all negative particles

    if (fDoULSandLS && !fDeactivateLS){
      // Calculated for single leg signals LS-- pairs
      for (unsigned int neg_i = 0; neg_i < fGenNegPart.size(); ++neg_i){
        for (unsigned int neg_j = neg_i + 1; neg_j < fGenNegPart.size(); ++neg_j){
          if (!fGenNegPart[neg_i].GetULSSignalPair() || !fGenNegPart[neg_j].GetULSSignalPair()) continue;
          bool selectedByKinematicCuts = true;
          if (fGenNegPart[neg_i].fPt < fPtMin || fGenNegPart[neg_i].fPt > fPtMax || fGenNegPart[neg_i].fEta < fEtaMin || fGenNegPart[neg_i].fEta > fEtaMax) selectedByKinematicCuts = false;
          if (fGenNegPart[neg_j].fPt < fPtMin || fGenNegPart[neg_j].fPt > fPtMax || fGenNegPart[neg_j].fEta < fEtaMin || fGenNegPart[neg_j].fEta > fEtaMax) selectedByKinematicCuts = false;
          bool selectedByKinematicCuts_smeared = true;
          if (fGenNegPart[neg_i].fPt_smeared < fPtMin || fGenNegPart[neg_i].fPt_smeared > fPtMax || fGenNegPart[neg_i].fEta_smeared < fEtaMin || fGenNegPart[neg_i].fEta_smeared > fEtaMax) selectedByKinematicCuts_smeared = false;
          if (fGenNegPart[neg_j].fPt_smeared < fPtMin || fGenNegPart[neg_j].fPt_smeared > fPtMax || fGenNegPart[neg_j].fEta_smeared < fEtaMin || fGenNegPart[neg_j].fEta_smeared > fEtaMax) selectedByKinematicCuts_smeared = false;

          if (selectedByKinematicCuts){
            TLorentzVector Lvec1;
            TLorentzVector Lvec2;
            Lvec1.SetPtEtaPhiM(fGenNegPart[neg_i].fPt, fGenNegPart[neg_i].fEta, fGenNegPart[neg_i].fPhi, AliPID::ParticleMass(AliPID::kElectron));
            Lvec2.SetPtEtaPhiM(fGenNegPart[neg_j].fPt, fGenNegPart[neg_j].fEta, fGenNegPart[neg_j].fPhi, AliPID::ParticleMass(AliPID::kElectron));
            TLorentzVector LvecM = Lvec1 + Lvec2;
            double mass = LvecM.M();
            double pairpt = LvecM.Pt();
            double weight = 1;
            for (unsigned int iMCSignal = 0; iMCSignal < fGenNegPart[neg_i].isMCSignal.size(); ++iMCSignal){
              // if (fGenNegPart[neg_i].isMCSignal[iMCSignal] == true && fGenNegPart[neg_j].isMCSignal[iMCSignal] == true)
              if (fGenNegPart[neg_i].isMCSignal[iMCSignal] == true && fGenNegPart[neg_j].isMCSignal[iMCSignal] == true &&
                  fGenNegPart[neg_i].DielectronPairFromSameMother[iMCSignal] == false && fGenNegPart[neg_j].DielectronPairFromSameMother[iMCSignal] == false)
               fHistGenPair_ULSandLS.at(3*iMCSignal+2)->Fill(mass, pairpt, weight * centralityWeight);
            }
          }
          if (selectedByKinematicCuts_smeared){
            TLorentzVector Lvec1;
            TLorentzVector Lvec2;
            Lvec1.SetPtEtaPhiM(fGenNegPart[neg_i].fPt_smeared, fGenNegPart[neg_i].fEta_smeared, fGenNegPart[neg_i].fPhi_smeared, AliPID::ParticleMass(AliPID::kElectron));
            Lvec2.SetPtEtaPhiM(fGenNegPart[neg_j].fPt_smeared, fGenNegPart[neg_j].fEta_smeared, fGenNegPart[neg_j].fPhi_smeared, AliPID::ParticleMass(AliPID::kElectron));
            TLorentzVector LvecM = Lvec1 + Lvec2;
            double mass = LvecM.M();
            double pairpt = LvecM.Pt();
            double weight = 1;
            for (unsigned int iMCSignal = 0; iMCSignal < fGenNegPart[neg_i].isMCSignal.size(); ++iMCSignal){
              // if (fGenNegPart[neg_i].isMCSignal[iMCSignal] == true && fGenNegPart[neg_j].isMCSignal[iMCSignal] == true)
              if (fGenNegPart[neg_i].isMCSignal[iMCSignal] == true && fGenNegPart[neg_j].isMCSignal[iMCSignal] == true &&
                  fGenNegPart[neg_i].DielectronPairFromSameMother[iMCSignal] == false && fGenNegPart[neg_j].DielectronPairFromSameMother[iMCSignal] == false)
               fHistGenSmearedPair_ULSandLS.at(3*iMCSignal+2)->Fill(mass, pairpt, weight * centralityWeight);

            }
          }
        }
      }

      // Calculated for single leg signals LS++ pairs
      for (unsigned int pos_i = 0; pos_i < fGenPosPart.size(); ++pos_i){
        for (unsigned int pos_j = pos_i + 1; pos_j < fGenPosPart.size(); ++pos_j){
          if (!fGenPosPart[pos_i].GetULSSignalPair() || !fGenPosPart[pos_j].GetULSSignalPair()) continue;
          bool selectedByKinematicCuts = true;
          if (fGenPosPart[pos_i].fPt < fPtMin || fGenPosPart[pos_i].fPt > fPtMax || fGenPosPart[pos_i].fEta < fEtaMin || fGenPosPart[pos_i].fEta > fEtaMax) selectedByKinematicCuts = false;
          if (fGenPosPart[pos_j].fPt < fPtMin || fGenPosPart[pos_j].fPt > fPtMax || fGenPosPart[pos_j].fEta < fEtaMin || fGenPosPart[pos_j].fEta > fEtaMax) selectedByKinematicCuts = false;
          bool selectedByKinematicCuts_smeared = true;
          if (fGenPosPart[pos_i].fPt_smeared < fPtMin || fGenPosPart[pos_i].fPt_smeared > fPtMax || fGenPosPart[pos_i].fEta_smeared < fEtaMin || fGenPosPart[pos_i].fEta_smeared > fEtaMax) selectedByKinematicCuts_smeared = false;
          if (fGenPosPart[pos_j].fPt_smeared < fPtMin || fGenPosPart[pos_j].fPt_smeared > fPtMax || fGenPosPart[pos_j].fEta_smeared < fEtaMin || fGenPosPart[pos_j].fEta_smeared > fEtaMax) selectedByKinematicCuts_smeared = false;

          if (selectedByKinematicCuts){
            TLorentzVector Lvec1;
            TLorentzVector Lvec2;
            Lvec1.SetPtEtaPhiM(fGenPosPart[pos_i].fPt, fGenPosPart[pos_i].fEta, fGenPosPart[pos_i].fPhi, AliPID::ParticleMass(AliPID::kElectron));
            Lvec2.SetPtEtaPhiM(fGenPosPart[pos_j].fPt, fGenPosPart[pos_j].fEta, fGenPosPart[pos_j].fPhi, AliPID::ParticleMass(AliPID::kElectron));
            TLorentzVector LvecM = Lvec1 + Lvec2;
            double mass = LvecM.M();
            double pairpt = LvecM.Pt();
            double weight = 1;
            for (unsigned int iMCSignal = 0; iMCSignal < fGenPosPart[pos_i].isMCSignal.size(); ++iMCSignal){
              // if (fGenPosPart[pos_i].isMCSignal[iMCSignal] == true && fGenPosPart[pos_j].isMCSignal[iMCSignal] == true)
              if (fGenPosPart[pos_i].isMCSignal[iMCSignal] == true && fGenPosPart[pos_j].isMCSignal[iMCSignal] == true &&
                  fGenPosPart[pos_i].DielectronPairFromSameMother[iMCSignal] == false && fGenPosPart[pos_j].DielectronPairFromSameMother[iMCSignal] == false)
               fHistGenPair_ULSandLS.at(3*iMCSignal+1)->Fill(mass, pairpt, weight * centralityWeight);
            }
          }
          if (selectedByKinematicCuts_smeared){
            TLorentzVector Lvec1;
            TLorentzVector Lvec2;
            Lvec1.SetPtEtaPhiM(fGenPosPart[pos_i].fPt_smeared, fGenPosPart[pos_i].fEta_smeared, fGenPosPart[pos_i].fPhi_smeared, AliPID::ParticleMass(AliPID::kElectron));
            Lvec2.SetPtEtaPhiM(fGenPosPart[pos_j].fPt_smeared, fGenPosPart[pos_j].fEta_smeared, fGenPosPart[pos_j].fPhi_smeared, AliPID::ParticleMass(AliPID::kElectron));
            TLorentzVector LvecM = Lvec1 + Lvec2;
            double mass = LvecM.M();
            double pairpt = LvecM.Pt();
            double weight = 1;
            for (unsigned int iMCSignal = 0; iMCSignal < fGenPosPart[pos_i].isMCSignal.size(); ++iMCSignal){
              // if (fGenPosPart[pos_i].isMCSignal[iMCSignal] == true && fGenPosPart[pos_j].isMCSignal[iMCSignal] == true)
              if (fGenPosPart[pos_i].isMCSignal[iMCSignal] == true && fGenPosPart[pos_j].isMCSignal[iMCSignal] == true &&
                  fGenPosPart[pos_i].DielectronPairFromSameMother[iMCSignal] == false && fGenPosPart[pos_j].DielectronPairFromSameMother[iMCSignal] == false)
               fHistGenSmearedPair_ULSandLS.at(3*iMCSignal+1)->Fill(mass, pairpt, weight * centralityWeight);

            }
          }
        }
      }
    }

    // ##########################################################
    // ##########################################################
    // ##########################################################
    // Fill reconstructed pairs
    for (unsigned int neg_i = 0; neg_i < fRecNegPart.size(); ++neg_i){
      for (unsigned int pos_i = 0; pos_i < fRecPosPart.size(); ++pos_i){
        AliVParticle* track1  = fEvent->GetTrack(fRecNegPart[neg_i].GetTrackID());
        AliVParticle* track2  = fEvent->GetTrack(fRecPosPart[pos_i].GetTrackID());

        if (fDoULSandLS && fRecNegPart[neg_i].GetULSSignalPair() && fRecPosPart[pos_i].GetULSSignalPair()){
          // Calculated for single leg signals unlike sign
          TLorentzVector Lvec1;
          TLorentzVector Lvec2;
          Lvec1.SetPtEtaPhiM(fRecNegPart[neg_i].fPt, fRecNegPart[neg_i].fEta, fRecNegPart[neg_i].fPhi, AliPID::ParticleMass(AliPID::kElectron));
          Lvec2.SetPtEtaPhiM(fRecPosPart[pos_i].fPt, fRecPosPart[pos_i].fEta, fRecPosPart[pos_i].fPhi, AliPID::ParticleMass(AliPID::kElectron));
          TLorentzVector LvecM = Lvec1 + Lvec2;
          double mass = LvecM.M();
          double pairpt = LvecM.Pt();

          for (unsigned int i = 0; i < fRecNegPart[neg_i].isMCSignal.size(); ++i){
            // if (fRecNegPart[neg_i].isMCSignal[i] == kTRUE && fRecPosPart[pos_i].isMCSignal[i] == kTRUE){
            if (fRecNegPart[neg_i].isMCSignal[i] == true && fRecPosPart[pos_i].isMCSignal[i] == true &&
                fRecNegPart[neg_i].DielectronPairFromSameMother[i] == false && fRecPosPart[pos_i].DielectronPairFromSameMother[i] == false){
              for (unsigned int j = 0; j < fRecNegPart[neg_i].isReconstructed.size(); ++j){
                if (fRecNegPart[neg_i].isReconstructed[j] == kTRUE && fRecPosPart[pos_i].isReconstructed[j] == kTRUE){
                  if (!fDeactivateLS) {
                    fHistRecPair_ULSandLS[j * 3 * fSingleLegMCSignal.size() + 3 * i]->Fill(mass, pairpt, centralityWeight);
                  }
                  else {
                    fHistRecPair_ULSandLS[j * 1 * fSingleLegMCSignal.size() + 1 * i]->Fill(mass, pairpt, centralityWeight);
                  }
                }// is selected by cutsetting
              } // end of loop over all cutsettings
            } // is selected by MC Signal
          } // end of loop over all MCsignals
        } // end of ULS loops

        // Check if electrons are from MCSignal Generator
        if (!fRecPosPart[pos_i].GetMCSignalPair() || !fRecNegPart[neg_i].GetMCSignalPair()) continue;

        // Apply MC signals
        std::vector<Bool_t> mcSignal_acc(fPairMCSignal.size(), kFALSE); // vector which stores if track is accepted by [i]-th mcsignal

        // Check if it according to mcsignals
        AliDielectronPair pair;
        pair.SetKFUsage(false);
        pair.SetTracks(static_cast<AliVTrack*>(track1), 11, static_cast<AliVTrack*>(track2), -11);
        for (unsigned int i = 0; i < fPairMCSignal.size(); ++i){
          mcSignal_acc[i] = AliDielectronMC::Instance()->IsMCTruth(&pair, &(fPairMCSignal[i]));
        }
        // check if at least one mc signal is true
        if (CheckIfOneIsTrue(mcSignal_acc) == kFALSE) continue;

        // Construct pair variables from LorentzVectors
        TLorentzVector Lvec1;
        TLorentzVector Lvec2;
        Lvec1.SetPtEtaPhiM(track1->Pt(), track1->Eta(), track1->Phi(), AliPID::ParticleMass(AliPID::kElectron));
        Lvec2.SetPtEtaPhiM(track2->Pt(), track2->Eta(), track2->Phi(), AliPID::ParticleMass(AliPID::kElectron));
        TLorentzVector LvecM = Lvec1 + Lvec2;
        double mass = LvecM.M();
        double pairpt = LvecM.Pt();
        double phiv = PhivPair(fEvent->GetMagneticField(),track1->Charge(),track2->Charge(),Lvec1.Vect(),Lvec2.Vect());

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

        // ##########################################################
        // Filling reconstructed particle histograms according to MCSignals
        for (unsigned int i = 0; i < mcSignal_acc.size(); ++i){
          if (mcSignal_acc[i] == kTRUE){
            for (unsigned int j = 0; j < fRecNegPart[neg_i].isReconstructed.size(); ++j){
              if (fRecNegPart[neg_i].isReconstructed[j] == kTRUE && fRecPosPart[pos_i].isReconstructed[j] == kTRUE){

                if(fDoFillPhiV) dynamic_cast<TH3D*>(fHistRecPair.at(j * mcSignal_acc.size() + i))->Fill(mass, pairpt, phiv ,weight * centralityWeight);//3D
                else{
                  if(fApplyPhivCut){
                    if(!(mass < fMaxMee && phiv > fMinPhiV)) dynamic_cast<TH2D*>(fHistRecPair.at(j * mcSignal_acc.size() + i))->Fill(mass, pairpt, weight * centralityWeight);//2D
                  }
                  else dynamic_cast<TH2D*>(fHistRecPair.at(j * mcSignal_acc.size() + i))->Fill(mass, pairpt, weight * centralityWeight);//2D
                }
                if (fWriteLegsFromPair){
                  ptNeg  = fRecNegPart[neg_i].fPt;
                  etaNeg = fRecNegPart[neg_i].fEta;
                  phiNeg = fRecNegPart[neg_i].fPhi;
                  ptPos  = fRecPosPart[pos_i].fPt;
                  etaPos = fRecPosPart[pos_i].fEta;
                  phiPos = fRecPosPart[pos_i].fPhi;
                  op_angle = Lvec2.Angle(Lvec1.Vect());

                  const double tuple[7] = {ptPos,etaPos,phiPos,ptNeg,etaNeg,phiNeg,op_angle};
                  fTHnSparseRecLegsFromPair.at(j * mcSignal_acc.size() + i)->Fill(tuple);
                }
              }// is selected by cutsetting
            } // end of loop over all cutsettings
          } // is selected by MCSignal
        } // end of loop over all MCsignals

      }// end of positive particle loop
    } // end of negative particle loop
  } // End of pairing

  if (fDoULSandLS && !fDeactivateLS){
    //LS--
    for (unsigned int neg_i = 0; neg_i < fRecNegPart.size(); ++neg_i){
      for (unsigned int neg_j = neg_i + 1; neg_j < fRecNegPart.size(); ++neg_j){
        if (!fRecNegPart[neg_i].GetULSSignalPair() || !fRecNegPart[neg_j].GetULSSignalPair()) continue;
        // Calculated for single leg signals like sign
        TLorentzVector Lvec1;
        TLorentzVector Lvec2;
        Lvec1.SetPtEtaPhiM(fRecNegPart[neg_i].fPt, fRecNegPart[neg_i].fEta, fRecNegPart[neg_i].fPhi, AliPID::ParticleMass(AliPID::kElectron));
        Lvec2.SetPtEtaPhiM(fRecNegPart[neg_j].fPt, fRecNegPart[neg_j].fEta, fRecNegPart[neg_j].fPhi, AliPID::ParticleMass(AliPID::kElectron));
        TLorentzVector LvecM = Lvec1 + Lvec2;
        double mass = LvecM.M();
        double pairpt = LvecM.Pt();

        for (unsigned int i = 0; i < fRecNegPart[neg_i].isMCSignal.size(); ++i){
          // if (fRecNegPart[neg_i].isMCSignal[i] == kTRUE && fRecNegPart[neg_j].isMCSignal[i] == kTRUE){
          if (fRecNegPart[neg_i].isMCSignal[i] == true && fRecNegPart[neg_j].isMCSignal[i] == true &&
              fRecNegPart[neg_i].DielectronPairFromSameMother[i] == false && fRecNegPart[neg_j].DielectronPairFromSameMother[i] == false){
            for (unsigned int j = 0; j < fRecNegPart[neg_i].isReconstructed.size(); ++j){
              if (fRecNegPart[neg_i].isReconstructed[j] == kTRUE && fRecNegPart[neg_j].isReconstructed[j] == kTRUE){
                fHistRecPair_ULSandLS[j * 3 * fSingleLegMCSignal.size() + 3 * i + 2]->Fill(mass, pairpt, centralityWeight);
              }// is selected by cutsetting
            } // end of loop over all cutsettings
          } // is selected by MC Signal
        } // end of loop over all MCsignals
      }
    }
    // LS++
    for (unsigned int pos_i = 0; pos_i < fRecPosPart.size(); ++pos_i){
      for (unsigned int pos_j = pos_i + 1; pos_j < fRecPosPart.size(); ++pos_j){
        if (!fRecPosPart[pos_i].GetULSSignalPair() || !fRecPosPart[pos_j].GetULSSignalPair()) continue;
        // Calculated for single leg signals like sign
        TLorentzVector Lvec1;
        TLorentzVector Lvec2;
        Lvec1.SetPtEtaPhiM(fRecPosPart[pos_i].fPt, fRecPosPart[pos_i].fEta, fRecPosPart[pos_i].fPhi, AliPID::ParticleMass(AliPID::kElectron));
        Lvec2.SetPtEtaPhiM(fRecPosPart[pos_j].fPt, fRecPosPart[pos_j].fEta, fRecPosPart[pos_j].fPhi, AliPID::ParticleMass(AliPID::kElectron));
        TLorentzVector LvecM = Lvec1 + Lvec2;
        double mass = LvecM.M();
        double pairpt = LvecM.Pt();

        for (unsigned int i = 0; i < fRecPosPart[pos_i].isMCSignal.size(); ++i){
          // if (fRecPosPart[pos_i].isMCSignal[i] == kTRUE && fRecPosPart[pos_j].isMCSignal[i] == kTRUE){
          if (fRecPosPart[pos_i].isMCSignal[i] == true && fRecPosPart[pos_j].isMCSignal[i] == true &&
              fRecPosPart[pos_i].DielectronPairFromSameMother[i] == false && fRecPosPart[pos_j].DielectronPairFromSameMother[i] == false){
            for (unsigned int j = 0; j < fRecPosPart[pos_i].isReconstructed.size(); ++j){
              if (fRecPosPart[pos_i].isReconstructed[j] == kTRUE && fRecPosPart[pos_j].isReconstructed[j] == kTRUE){
                fHistRecPair_ULSandLS[j * 3 * fSingleLegMCSignal.size() + 3 * i + 1]->Fill(mass, pairpt, centralityWeight);
              }// is selected by cutsetting
            } // end of loop over all cutsettings
          } // is selected by MC Signal
        } // end of loop over all MCsignals
      }
    }


  }

  PostData(1, fOutputList);
}


// ############################################################################
// ############################################################################
void    AliAnalysisTaskElectronEfficiencyV2::FillTrackHistograms(AliVParticle* track, AliVParticle* mcTrack){
  Double_t values[AliDielectronVarManager::kNMaxValues]={0.};
  AliDielectronVarManager::SetFillMap(fUsedVars); // currently filled manually in the constructor of this task.

  AliDielectronVarManager::Fill(track, values);
  // std::cout << "pt var  manager = " << values[AliDielectronVarManager::kPt] << std::endl;
  // std::cout << "SITS    manager = " << values[AliDielectronVarManager::kNclsSITS] << std::endl;
  // std::cout << "TPCnSig manager = " << values[AliDielectronVarManager::kTPCnSigmaEle] << std::endl;
  // std::cout << fOutputListSupportHistos << std::endl;
  TString genname;
  (dynamic_cast<TH1D *>(fOutputListSupportHistos->At(0)))->Fill(values[AliDielectronVarManager::kPt]);//hPt (reco)
  (dynamic_cast<TH2D *>(fOutputListSupportHistos->At(1)))->Fill(values[AliDielectronVarManager::kP],   values[AliDielectronVarManager::kITSnSigmaEle]);
  (dynamic_cast<TH2D *>(fOutputListSupportHistos->At(2)))->Fill(values[AliDielectronVarManager::kPIn], values[AliDielectronVarManager::kTPCnSigmaEle]);
  (dynamic_cast<TH2D *>(fOutputListSupportHistos->At(3)))->Fill(values[AliDielectronVarManager::kP],   values[AliDielectronVarManager::kTOFnSigmaEle]);
  (dynamic_cast<TH1D *>(fOutputListSupportHistos->At(4)))->Fill(values[AliDielectronVarManager::kEta]);
  (dynamic_cast<TH1D *>(fOutputListSupportHistos->At(5)))->Fill(values[AliDielectronVarManager::kPhi]);
  (dynamic_cast<TH2D *>(fOutputListSupportHistos->At(6)))->Fill(values[AliDielectronVarManager::kEta], values[AliDielectronVarManager::kPhi]);

  (dynamic_cast<TH1D *>(fOutputListSupportHistos->At(7 )))->Fill(values[AliDielectronVarManager::kNFclsTPCfCross]);
  (dynamic_cast<TH1D *>(fOutputListSupportHistos->At(8 )))->Fill(values[AliDielectronVarManager::kNFclsTPCr]);
  (dynamic_cast<TH1D *>(fOutputListSupportHistos->At(9 )))->Fill(values[AliDielectronVarManager::kNclsTPC]);
  (dynamic_cast<TH1D *>(fOutputListSupportHistos->At(10)))->Fill(values[AliDielectronVarManager::kNclsITS]);
  (dynamic_cast<TH1D *>(fOutputListSupportHistos->At(11)))->Fill(values[AliDielectronVarManager::kTPCchi2Cl]);
  (dynamic_cast<TH1D *>(fOutputListSupportHistos->At(12)))->Fill(values[AliDielectronVarManager::kITSchi2Cl]);
  (dynamic_cast<TH1D *>(fOutputListSupportHistos->At(13)))->Fill(values[AliDielectronVarManager::kNclsSTPC]);
  (dynamic_cast<TH1D *>(fOutputListSupportHistos->At(14)))->Fill(values[AliDielectronVarManager::kNclsSITS]);
  (dynamic_cast<TH1D *>(fOutputListSupportHistos->At(15)))->Fill(values[AliDielectronVarManager::kNclsSFracTPC]);
  (dynamic_cast<TH1D *>(fOutputListSupportHistos->At(16)))->Fill(values[AliDielectronVarManager::kNclsSFracITS]);
  (dynamic_cast<TH1D *>(fOutputListSupportHistos->At(17)))->Fill(values[AliDielectronVarManager::kTPCclsDiff]);
  (dynamic_cast<TH1D *>(fOutputListSupportHistos->At(18)))->Fill(values[AliDielectronVarManager::kTPCsignalN]);
  (dynamic_cast<TH2D *>(fOutputListSupportHistos->At(19)))->Fill(values[AliDielectronVarManager::kNclsTPC], values[AliDielectronVarManager::kNFclsTPCr]);
  (dynamic_cast<TH2D *>(fOutputListSupportHistos->At(20)))->Fill(values[AliDielectronVarManager::kPt], values[AliDielectronVarManager::kNFclsTPCr]);
  // (dynamic_cast<TH1D *>(fOutputListSupportHistos->At(21)))->Fill(values[AliDielectronVarManager::kPdgCode]);
  (dynamic_cast<TH1D *>(fOutputListSupportHistos->At(21)))->Fill(mcTrack->PdgCode());
  (dynamic_cast<TH1D *>(fOutputListSupportHistos->At(22)))->Fill( (fMC->GetTrack(TMath::Abs(mcTrack->GetMother())))->PdgCode());
  if(fMC->GetCocktailGenerator(TMath::Abs(track->GetLabel()), genname))    (dynamic_cast<TH1D *>(fOutputListSupportHistos->At(23)))->Fill( genname,1);
  else (dynamic_cast<TH1D *>(fOutputListSupportHistos->At(23)))->Fill( "none",1);
}


// ############################################################################
// ############################################################################
AliAnalysisTaskElectronEfficiencyV2::Particle AliAnalysisTaskElectronEfficiencyV2::CreateParticle(AliVParticle* mcPart1){
  double  pt1      = mcPart1->Pt();
  double  eta1     = mcPart1->Eta();
  double  phi1     = mcPart1->Phi();
  short   charge1  = mcPart1->Charge();
  Particle part(pt1, eta1, phi1, charge1);
  part.DielectronPairFromSameMother.resize(fDielectronPairNotFromSameMother.size(), false);


  return part;
}

void AliAnalysisTaskElectronEfficiencyV2::CheckIfFromMotherWithDielectronAsDaughter(Particle& part){
  if (isAOD && fDoULSandLS){

    for (unsigned int k = 0; k < fDielectronPairNotFromSameMother.size(); ++k){
      if (part.isMCSignal[k] == true && fDielectronPairNotFromSameMother[k] == true){
        AliAODMCParticle* mother = dynamic_cast<AliAODMCParticle*> (fMC->GetTrack(part.GetMotherID()));
        // int number_of_daugthers = mother->GetNDaughters() ;
        int LabelFirstDaughter = mother->GetDaughterFirst();
        int LabelLastDaughter = mother->GetDaughterLast();
        // std::cout << "number_of_daughters = " << number_of_daugthers << "  first_daugther = " << LabelFirstDaughter << "  last_daugther = " << LabelLastDaughter << std::endl;

        bool ele_from_same_mother = false;
        bool pos_from_same_mother = false;
        for (int daughter_i = LabelFirstDaughter; daughter_i <= LabelLastDaughter; daughter_i++){
          int pdgCode = fMC->GetTrack(daughter_i)->PdgCode();
          if      (pdgCode == 11)  ele_from_same_mother = true;
          else if (pdgCode == -11) pos_from_same_mother = true;
          // std::cout << "daugther[" << daughter_i << "] with pdgcode: " << pdgCode << std::endl;
        }
        if (ele_from_same_mother == true && pos_from_same_mother == true) {
          part.DielectronPairFromSameMother[k] = true;
          // std::cout << "dielectron pair from same mother" << std::endl;
        }
        else{
          part.DielectronPairFromSameMother[k] = false;
        }
      }
      else{
        part.DielectronPairFromSameMother[k] = false;
      }
    }
  }
  if (part.DielectronPairFromSameMother.size() != fDielectronPairNotFromSameMother.size()){
    std::cout << "ERROR IN SOME PART" << std::endl;
    // vec = std::vector<bool>(fDielectronPairNotFromSameMother.size(), false);
  }
  // part.DielectronPairFromSameMother = vec;

}


// ############################################################################
// ############################################################################
Bool_t AliAnalysisTaskElectronEfficiencyV2::CheckIfOneIsTrue(std::vector<Bool_t>& vec){
  bool min_one_is_true = kFALSE;
  unsigned int size = vec.size();
  for (unsigned int i = 0; i < size; ++i){
    if (vec[i] == kTRUE) {min_one_is_true = kTRUE; break;}
  }
  return min_one_is_true;
}


// ############################################################################
// ############################################################################
void AliAnalysisTaskElectronEfficiencyV2::SetBinsLinear(const std::string var, const double min, const double max, const unsigned int steps){
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
  else if (var == "phiv") fPhiVBins.clear();

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
    else if (var == "phiv") fPhiVBins.push_back(i * stepSize + min);
  }
}


// ############################################################################
// ############################################################################
void AliAnalysisTaskElectronEfficiencyV2::CheckSingleLegMCsignals(std::vector<Bool_t>& vec, const int tracklabel){
  for (unsigned int i = 0; i < fSingleLegMCSignal.size(); ++i){
    vec.at(i) = AliDielectronMC::Instance()->IsMCTruth(tracklabel, &(fSingleLegMCSignal[i]), 1);
  }
}

// ############################################################################
// ############################################################################
void AliAnalysisTaskElectronEfficiencyV2::CreateSupportHistos()
{
  Printf(" Now running: CreateSupportHistos()");

  fOutputListSupportHistos = new TList();
  fOutputListSupportHistos->SetName("Support");
  fOutputListSupportHistos->SetOwner();


  // Track variables
  TH1D* hPt      = new TH1D("Pt","Pt;Pt [GeV];#tracks",200,0.,10.);//,AliDielectronVarManager::kPt);
  fOutputListSupportHistos->AddAt(hPt,     0);

  // PID
  TH2D* hITSnSigmaEle_P = new TH2D("ITSnSigmaEle_P","ITS number of sigmas Electrons;P [GeV/c];ITS number of sigmas Electrons", 200,0.,10.,100,-5.,5.);//.,AliDielectronVarManager::kP,AliDielectronVarManager::kITSnSigmaEle,makeLogx);
  fOutputListSupportHistos->AddAt(hITSnSigmaEle_P, 1);
  TH2D* hTPCnSigmaEle_P = new TH2D("TPCnSigmaEle_P","TPC number of sigmas Electrons;PIn (pTPC) [GeV/c];TPC number of sigmas Electrons", 200,0.,10.,100,-5.,5.);//.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCnSigmaEle,makeLogx);
  fOutputListSupportHistos->AddAt(hTPCnSigmaEle_P, 2);
  TH2D* hTOFnSigmaEle_P = new TH2D("TOFnSigmaEle_P","TOF number of sigmas Electrons;PIn (pTPC) [GeV/c];TOF number of sigmas Electrons", 200,0.,10.,100,-5.,5.);//,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTOFnSigmaEle,makeLogx);
  fOutputListSupportHistos->AddAt(hTOFnSigmaEle_P, 3);

  // Track kinematic
  TH1D* hEta = new TH1D("Eta","Eta; Eta;#tracks", 200, -2, 2);//,AliDielectronVarManager::kEta);
  TH1D* hPhi = new TH1D("Phi","Phi; Phi;#tracks", 320, 0., 6.4);//,AliDielectronVarManager::kPhi);
  TH2D* hEta_Phi = new TH2D("Eta_Phi","Eta Phi Map; Eta; Phi", 100, -1, 1, 320, 0, 6.4);//,AliDielectronVarManager::kEta,AliDielectronVarManager::kPhi);
  fOutputListSupportHistos->AddAt(hEta, 4);
  fOutputListSupportHistos->AddAt(hPhi, 5);
  fOutputListSupportHistos->AddAt(hEta_Phi, 6);

  // Quality
  TH1D* hTPCcrossedRowsOverFindable = new TH1D("TPCcrossedRowsOverFindable","Number of Crossed Rows TPC over Findable;TPC crossed rows over findable;#tracks",120,0.,1.2);//,AliDielectronVarManager::kNFclsTPCfCross);
  TH1D* hTPCcrossedRows = new TH1D("TPCcrossedRows","Number of Crossed Rows TPC;TPC crossed rows;#tracks",160,-0.5,159.5);//,AliDielectronVarManager::kNFclsTPCr);
  TH1D* hTPCnCls = new TH1D("TPCnCls","Number of Clusters TPC;TPC number clusters;#tracks",160,-0.5,159.5);//,AliDielectronVarManager::kNclsTPC);
  TH1D* hITSnCls = new TH1D("ITSnCls","Number of Clusters ITS;ITS number clusters;#tracks",10,-0.5,9.5);//,AliDielectronVarManager::kNclsITS);
  TH1D* hTPCchi2 = new TH1D("TPCchi2","TPC chi2 per Cluster;TPC chi2/Cl;#tracks",100,0.,10.);//,AliDielectronVarManager::kTPCchi2Cl);
  TH1D* hITSchi2 = new TH1D("ITSchi2","ITS chi2 per Cluster;ITS chi2/Cl;#tracks",100,0.,10.);//,AliDielectronVarManager::kITSchi2Cl);
  TH1D* hTPCnClsS = new TH1D("TPCnClsS",";TPC number of shared clusters;#tracks",160,-0.5,159.5);//,AliDielectronVarManager::kNclsSTPC);
  TH1D* hITSnClsS = new TH1D("ITSnClsS",";ITS number of shared clusters;#tracks",7,-0.5,6.5);//,AliDielectronVarManager::kNclsSITS);
  TH1D* hNclsSFracTPC = new TH1D("NclsSFracTPC","Fraction of shared clusters assigned in the TPC;TPC fraction of shared clusters;#tracks",120,0,1.2);//.,AliDielectronVarManager::kNclsSFracTPC);
  TH1D* hNclsSFracITS = new TH1D("NclsSFracITS","Fraction of shared clusters assigned in the ITS;ITS fraction of shared clusters;#tracks",120,0,1.2);//.,AliDielectronVarManager::kNclsSFracITS);
  TH1D* hTPCclsDiff = new TH1D("TPCclsDiff","TPC cluster difference;N_{d#it{E}/d#it{x} points}^{TPC} - N_{cls}^{TPC};#tracks",100,-80,20);//.,AliDielectronVarManager::kTPCclsDiff);
  TH1D* hTPCsignalN = new TH1D("TPCsignalN","Number of PID Clusters TPC;N_{d#it{E}/d#it{x} points}^{TPC};#tracks",160,-0.5,159.5);//.,AliDielectronVarManager::kTPCsignalN); //kNclsTPCdEdx
  fOutputListSupportHistos->AddAt(hTPCcrossedRowsOverFindable, 7);
  fOutputListSupportHistos->AddAt(hTPCcrossedRows, 8);
  fOutputListSupportHistos->AddAt(hTPCnCls, 9);
  fOutputListSupportHistos->AddAt(hITSnCls, 10);
  fOutputListSupportHistos->AddAt(hTPCchi2, 11);
  fOutputListSupportHistos->AddAt(hITSchi2, 12);
  fOutputListSupportHistos->AddAt(hTPCnClsS, 13);
  fOutputListSupportHistos->AddAt(hITSnClsS, 14);
  fOutputListSupportHistos->AddAt(hNclsSFracTPC, 15);
  fOutputListSupportHistos->AddAt(hNclsSFracITS, 16);
  fOutputListSupportHistos->AddAt(hTPCclsDiff, 17);
  fOutputListSupportHistos->AddAt(hTPCsignalN, 18);
  TH2D* hTPCcrossedRows_TPCnCls = new TH2D("TPCcrossedRows_TPCnCls","TPC crossed rows vs TPC number clusters;TPC number clusters;TPC crossed rows",
                                           160,-0.5,159.5,160,-0.5,159.5);//,AliDielectronVarManager::kNclsTPC,AliDielectronVarManager::kNFclsTPCr);
  TH2D* hTPCcrossedRows_Pt = new TH2D("TPCcrossedRows_Pt","TPC crossed rows vs Pt;Pt [GeV];TPC crossed rows",
                                      200,0.,10.,160,-0.5,159.5);//,AliDielectronVarManager::kPt,AliDielectronVarManager::kNFclsTPCr);
  fOutputListSupportHistos->AddAt(hTPCcrossedRows_TPCnCls, 19);
  fOutputListSupportHistos->AddAt(hTPCcrossedRows_Pt, 20);

  TH1D* hPDGCode = new TH1D("PDGCode","PDGCode;#tracks",10001, -5000, 5000);//.,AliDielectronVarManager::kTPCsignalN); //kNclsTPCdEdx
  fOutputListSupportHistos->AddAt(hPDGCode, 21);

  TH1D* hPDGCodeMother = new TH1D("PDGCodeMother","PDGCodeMother;#tracks",10001, -5000, 5000);//.,AliDielectronVarManager::kTPCsignalN); //kNclsTPCdEdx
  fOutputListSupportHistos->AddAt(hPDGCodeMother, 22);
  // TH2D* hPDGCode_PDGCodeMother = new TH2D("PDGCode_PDGCodeMother",";PDG code;PDG code Mother",
  // 10001,-5000,5000,10001,-5000,5000);//,AliDielectronVarManager::kPt,AliDielectronVarManager::kNFclsTPCr);
  // fOutputListSupportHistos->AddAt(hPDGCode_PDGCodeMother, 21);
  TH1D* hMCGenCode = new TH1D("MCGenerator","MCGenerator;#tracks",1, 0, 0);//.,AliDielectronVarManager::kTPCsignalN); //kNclsTPCdEdx
  fOutputListSupportHistos->AddAt(hMCGenCode, 23);

}


// ############################################################################
// ############################################################################
TLorentzVector AliAnalysisTaskElectronEfficiencyV2::ApplyResolution(double pt, double eta, double phi, short ch) {
  // from Theos LightFlavorGenerator, modified

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
Double_t AliAnalysisTaskElectronEfficiencyV2::GetSmearing(TObjArray *arr, Double_t x)
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

bool AliAnalysisTaskElectronEfficiencyV2::CheckGenerator(int trackID, std::vector<unsigned int> vecHashes){
  if (vecHashes.size() == 0) return true;

  if(isAOD){//for AOD
    AliAODMCHeader* mcHeader = (AliAODMCHeader*)fEvent->GetList()->FindObject(AliAODMCHeader::StdBranchName());
    if (!mcHeader) {
      AliError("Could not find MC headr in AOD");
      return false;
    }
    TClonesArray* mcArray = dynamic_cast<TClonesArray*>(fEvent->FindListObject(AliAODMCParticle::StdBranchName()));
    if (!mcArray) {
      AliError("Could not find MC array in AOD");
      return false;
    }
    if(AliAnalysisUtils::IsParticleFromOutOfBunchPileupCollision(TMath::Abs(trackID), mcHeader, mcArray)) return false;//particles from pileup collision should NOT be used.
  }
  else{//for ESD
    if(AliAnalysisUtils::IsParticleFromOutOfBunchPileupCollision(TMath::Abs(trackID), fMC)) return false;//particles from pileup collision should NOT be used.
  }

  TString genname="";
  Bool_t hasGenerator = fMC->GetCocktailGenerator(TMath::Abs(trackID), genname);
  //AliMCParticle* p = (AliMCParticle*)fMC->GetTrack(TMath::Abs(trackID));
  //Int_t genID = p->GetGeneratorIndex();
  //AliInfo(Form("genID = %d , generator name = %s",genID,genname.Data()));

  //std::cout << genname << std::endl;
  if(!hasGenerator) {
    Printf("no cocktail header list was found for this track");
    return false;
  }
  else{
    for (unsigned int i = 0; i < vecHashes.size(); ++i){
      //std::cout << genname.Hash() << " " << vecHashes[i] << std::endl;
      if (genname.Hash() == vecHashes[i]) return true;
    }//end of vecHashes loop
    return false;
  }
  return false; // should not happen
}

bool AliAnalysisTaskElectronEfficiencyV2::CheckGeneratorIndex(int trackID, std::vector<unsigned int> vecGenIDs){
  if (vecGenIDs.size() == 0) return true;

  if(isAOD){//for AOD
    AliAODMCHeader* mcHeader = (AliAODMCHeader*)fEvent->GetList()->FindObject(AliAODMCHeader::StdBranchName());
    if (!mcHeader) {
      AliError("Could not find MC headr in AOD");
      return false;
    }
    TClonesArray* mcArray = dynamic_cast<TClonesArray*>(fEvent->FindListObject(AliAODMCParticle::StdBranchName()));
    if (!mcArray) {
      AliError("Could not find MC array in AOD");
      return false;
    }
    if(AliAnalysisUtils::IsParticleFromOutOfBunchPileupCollision(TMath::Abs(trackID), mcHeader, mcArray)) return false;//particles from pileup collision should NOT be used.
  }
  else{//for ESD
    if(AliAnalysisUtils::IsParticleFromOutOfBunchPileupCollision(TMath::Abs(trackID), fMC)) return false;//particles from pileup collision should NOT be used.
  }

  AliMCParticle* p = (AliMCParticle*)fMC->GetTrack(TMath::Abs(trackID));
  Int_t genID = p->GetGeneratorIndex();

  //TString genname="";
  //Bool_t hasGenerator = fMC->GetCocktailGenerator(TMath::Abs(trackID), genname);
  //AliInfo(Form("genID = %d , generator name = %s",genID,genname.Data()));

  for (unsigned int i = 0; i < vecGenIDs.size(); ++i){
    //std::cout << genID << " " << vecGenIDs[i] << std::endl;
    if (genID == vecGenIDs[i]) return true;
  }//end of vecHashes loop
  return false;
}

double AliAnalysisTaskElectronEfficiencyV2::GetWeight(Particle part1, Particle part2, double motherpt){
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
void AliAnalysisTaskElectronEfficiencyV2::SetCentroidCorrFunction(Detector det, TObject *fun, UInt_t varx, UInt_t vary, UInt_t varz)
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
void AliAnalysisTaskElectronEfficiencyV2::SetWidthCorrFunction(Detector det, TObject *fun, UInt_t varx, UInt_t vary, UInt_t varz)
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
//______________________________________________
Double_t AliAnalysisTaskElectronEfficiencyV2::PhivPair(Double_t MagField, Int_t charge1, Int_t charge2, TVector3 dau1, TVector3 dau2) //const
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
//____________________________
Double_t AliAnalysisTaskElectronEfficiencyV2::CalculateNbins() {

  //
  // Calculate the number of bins not emptied in the correction centrality histo in the centrality range we are interested in
  //

  Double_t nbins=0.;

  if(!fHistCentralityCorrection) return nbins;

  
  for(Int_t k=0; k < fHistCentralityCorrection->GetNbinsX(); k++) {
     
    TAxis *xaxis = fHistCentralityCorrection->GetXaxis();
    Double_t cent = xaxis->GetBinCenter(k+1);

    if((fMinCentrality < 0.) || (fMaxCentrality < 0.)) {
      if( fHistCentralityCorrection->GetBinContent(k+1)>0.) nbins = nbins + 1.;
    }
    else {
      if((cent>=fMinCentrality) && (cent<=fMaxCentrality)) {
	if( fHistCentralityCorrection->GetBinContent(k+1)>0.) nbins = nbins + 1.;
      }
    }
    
  }
    
   
  return nbins;
  
}
//________________________________________________________________________
