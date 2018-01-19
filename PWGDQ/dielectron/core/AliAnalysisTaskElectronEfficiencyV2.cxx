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


#include "AliDielectronMC.h"
#include "AliDielectronVarManager.h"
#include "AliDielectronSignalMC.h"
#include "AliDielectronPair.h"

#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TLorentzVector.h"

#include "TChain.h"
#include "TSystem.h"

#include <iostream>

AliAnalysisTaskElectronEfficiencyV2::~AliAnalysisTaskElectronEfficiencyV2(){
  delete fOutputList;
  delete fSingleElectronList;
  delete fPairList;
  delete fResolutionList;
}

// ############################################################################
// ############################################################################
AliAnalysisTaskElectronEfficiencyV2::AliAnalysisTaskElectronEfficiencyV2(): AliAnalysisTaskSE(), fEventFilter(0x0)
                                                                              , fResoFile(0x0), fResoFilename(""), fResoFilenameFromAlien(""), fArrResoPt(0x0), fArrResoEta(0x0), fArrResoPhi_Pos(0x0), fArrResoPhi_Neg(0x0)
                                                                              , fOutputList(0x0), fSingleElectronList(0x0), fPairList(0x0), fResolutionList(0x0)
                                                                              , fPGen_DeltaP(0x0), fPtGen_DeltaPt(0x0), fPtGen_DeltaPtOverPtGen(0x0), fPGen_PrecOverPGen(0x0), fPtGen_PtRecOverPtGen(0x0)
                                                                              , fPGen_DeltaEta(0x0), fPtGen_DeltaEta(0x0), fPGen_DeltaTheta(0x0), fPGen_DeltaPhi_Ele(0x0), fPGen_DeltaPhi_Pos(0x0), fPtGen_DeltaPhi_Ele(0x0)
                                                                              , fPtGen_DeltaPhi_Pos(0x0), fThetaGen_DeltaTheta(0x0), fPhiGen_DeltaPhi(0x0)
                                                                              , fPtBins(), fEtaBins(), fPhiBins(), fThetaBins()
                                                                              , fResolutionDeltaPtBins(), fResolutionRelPtBins(), fResolutionEtaBins(), fResolutionPhiBins(), fResolutionThetaBins()
                                                                              , fMassBins(), fPairPtBins()
                                                                              , fPtMin(0.), fPtMax(0.), fEtaMin(-99.), fEtaMax(99.)
                                                                              , fPtMinGen(0.), fPtMaxGen(0.), fEtaMinGen(-99.), fEtaMaxGen(99.)
                                                                              , fSingleLegMCSignal(), fPairMCSignal()
                                                                              , fPIDResponse(0x0), fEvent(0x0), fMC(0x0), fTrack(0x0), isAOD(false), fSelectPhysics(false), fTriggerMask(0)
                                                                              , fTrackCuts(), fUsedVars(0x0)
                                                                              , fSupportMCSignal(0), fSupportCutsetting(0)
                                                                              , fHistEvents(0x0), fHistEventStat(0x0), fHistCentrality(0x0), fHistVertex(0x0), fHistVertexContibutors(0x0), fHistNTracks(0x0)
                                                                              , fMinCentrality(0.), fMaxCentrality(100), fCentralityFile(0x0), fCentralityFilename(""), fHistCentralityCorrection(0x0)
                                                                              , fOutputListSupportHistos(0x0)
                                                                              , fHistGenPosPart(), fHistGenNegPart(), fHistGenSmearedPosPart(), fHistGenSmearedNegPart(), fHistRecPosPart(), fHistRecNegPart()
                                                                              , fHistGenPair(), fHistGenSmearedPair(), fHistRecPair(), fHistGenPair_ULSandLS(), fHistGenSmearedPair_ULSandLS(), fHistRecPair_ULSandLS()
                                                                              , fDoPairing(false), fDoULSandLS(false)
                                                                              , fGenNegPart(), fGenPosPart(), fRecNegPart(), fRecPosPart()
{
// ROOT IO constructor , don ’t allocate memory here !
}


// ############################################################################
// ############################################################################
AliAnalysisTaskElectronEfficiencyV2::AliAnalysisTaskElectronEfficiencyV2(const char * name) : AliAnalysisTaskSE(name), fEventFilter(0x0)
                                                                              , fResoFile(0x0), fResoFilename(""), fResoFilenameFromAlien(""), fArrResoPt(0x0), fArrResoEta(0x0), fArrResoPhi_Pos(0x0), fArrResoPhi_Neg(0x0)
                                                                              , fOutputList(0x0), fSingleElectronList(0x0), fPairList(0x0), fResolutionList(0x0)
                                                                              , fPGen_DeltaP(0x0), fPtGen_DeltaPt(0x0), fPtGen_DeltaPtOverPtGen(0x0), fPGen_PrecOverPGen(0x0), fPtGen_PtRecOverPtGen(0x0)
                                                                              , fPGen_DeltaEta(0x0), fPtGen_DeltaEta(0x0), fPGen_DeltaTheta(0x0), fPGen_DeltaPhi_Ele(0x0), fPGen_DeltaPhi_Pos(0x0), fPtGen_DeltaPhi_Ele(0x0)
                                                                              , fPtGen_DeltaPhi_Pos(0x0), fThetaGen_DeltaTheta(0x0), fPhiGen_DeltaPhi(0x0)
                                                                              , fPtBins(), fEtaBins(), fPhiBins(), fThetaBins()
                                                                              , fResolutionDeltaPtBins(), fResolutionRelPtBins(), fResolutionEtaBins(), fResolutionPhiBins(), fResolutionThetaBins()
                                                                              , fMassBins(), fPairPtBins()
                                                                              , fPtMin(0.), fPtMax(0.), fEtaMin(-99.), fEtaMax(99.)
                                                                              , fPtMinGen(0.), fPtMaxGen(0.), fEtaMinGen(-99.), fEtaMaxGen(99.)
                                                                              , fSingleLegMCSignal(), fPairMCSignal()
                                                                              , fPIDResponse(0x0), fEvent(0x0), fMC(0x0), fTrack(0x0), isAOD(false), fSelectPhysics(false), fTriggerMask(0)
                                                                              , fTrackCuts(), fUsedVars(0x0)
                                                                              , fSupportMCSignal(0), fSupportCutsetting(0)
                                                                              , fHistEvents(0x0), fHistEventStat(0x0), fHistCentrality(0x0), fHistVertex(0x0), fHistVertexContibutors(0x0), fHistNTracks(0x0)
                                                                              , fMinCentrality(0.), fMaxCentrality(100), fCentralityFile(0x0), fCentralityFilename(""), fHistCentralityCorrection(0x0)
                                                                              , fOutputListSupportHistos(0x0)
                                                                              , fHistGenPosPart(), fHistGenNegPart(), fHistGenSmearedPosPart(), fHistGenSmearedNegPart(), fHistRecPosPart(), fHistRecNegPart()
                                                                              , fHistGenPair(), fHistGenSmearedPair(), fHistRecPair(), fHistGenPair_ULSandLS(), fHistGenSmearedPair_ULSandLS(), fHistRecPair_ULSandLS()
                                                                              , fDoPairing(false), fDoULSandLS(false)
                                                                              , fGenNegPart(), fGenPosPart(), fRecNegPart(), fRecPosPart()
{
  DefineInput (0, TChain::Class());
  DefineOutput (1, TList::Class());

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

}


// ############################################################################
// ############################################################################
void AliAnalysisTaskElectronEfficiencyV2::Terminate(Option_t* option){
  // fHistEventStat->SetAxisRange(0., fHistEventStat->GetMaximum() * 1.1, "Y");
}


// ############################################################################
// ############################################################################
void AliAnalysisTaskElectronEfficiencyV2::UserCreateOutputObjects(){
  std::cout << "Starting UserCreateOutputObjects()" << std::endl;


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
        fPairList->Add(GeneratedPairs);

        TList* GeneratedSmearedPairs = new TList();
        GeneratedSmearedPairs->SetName("GeneratedSmeared");
        GeneratedSmearedPairs->SetOwner();
        for (unsigned int i = 0; i < fPairMCSignal.size(); ++i){
          TH2D* th2_tmp = new TH2D(Form("Ngen_%s", fPairMCSignal.at(i).GetName()),";m_{ee};p_{T,ee}",fNmassBins,fMassBins.data(),fNpairptBins,fPairPtBins.data());
          th2_tmp->Sumw2();
          fHistGenSmearedPair.push_back(th2_tmp);
          GeneratedSmearedPairs->Add(th2_tmp);
        }
        if (fDoULSandLS == true){
          for (unsigned int i = 0; i < fSingleLegMCSignal.size(); ++i){
            TH2D* th2_tmp_ULS = new TH2D(Form("Ngen_ULS_%s", fSingleLegMCSignal.at(i).GetName()),";m_{ee};p_{T,ee}",fNmassBins,fMassBins.data(),fNpairptBins,fPairPtBins.data());
            th2_tmp_ULS->Sumw2();
            fHistGenSmearedPair_ULSandLS.push_back(th2_tmp_ULS);
            GeneratedSmearedPairs->Add(th2_tmp_ULS);
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
        fPairList->Add(GeneratedSmearedPairs);

        // Generated reconstructed lists for every cutsetting one list and every MCsignal 1 histogram
        for (unsigned int list_i = 0; list_i < fTrackCuts.size(); ++list_i){
          TList* list = new TList();
          list->SetName(fTrackCuts.at(list_i)->GetName());
          list->SetOwner();

          for (unsigned int i = 0; i < fPairMCSignal.size(); ++i){
            TH2D* th2_tmp = new TH2D(Form("Nrec_%s", fPairMCSignal.at(i).GetName()),";m_{ee};p_{T,ee}",fNmassBins,fMassBins.data(),fNpairptBins,fPairPtBins.data());
            th2_tmp->Sumw2();
            fHistRecPair.push_back(th2_tmp);
            list->Add(th2_tmp);
          }
          if (fDoULSandLS == true){
            for (unsigned int i = 0; i < fSingleLegMCSignal.size(); ++i){
              TH2D* th2_tmp_ULS = new TH2D(Form("Nrec_ULS_%s", fSingleLegMCSignal.at(i).GetName()),";m_{ee};p_{T,ee}",fNmassBins,fMassBins.data(),fNpairptBins,fPairPtBins.data());
              th2_tmp_ULS->Sumw2();
              fHistRecPair_ULSandLS.push_back(th2_tmp_ULS);
              list->Add(th2_tmp_ULS);
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
          fPairList->Add(list);
        }
      }

    // ######################################################
    // #####################  Resolution #########################
    // ######################################################
    fResolutionList = new TList();
    fResolutionList->SetName("Resolution");
    fResolutionList->SetOwner();

      fPGen_DeltaP                         = new TH2D("PGen_DeltaP",                        "", 500, 0., 10., fNResolutionDeltaptBins, fResolutionDeltaPtBins.data());
      fPtGen_DeltaPt                       = new TH2D("PtGen_DeltaPt",                      "", 500, 0., 10., fNResolutionDeltaptBins, fResolutionDeltaPtBins.data());
      fPtGen_DeltaPtOverPtGen              = new TH2D("PtGen_DeltaPtOverPtGen",             "", 500, 0., 10., fNResolutionDeltaptBins, -1., +1.);
      fPGen_PrecOverPGen                   = new TH2D("PGen_PrecOverPGen",                  "", 500, 0., 10., fNResolutionRelptBins, fResolutionRelPtBins.data());
      fPtGen_PtRecOverPtGen                = new TH2D("PtGen_PtRecOverPtGen",               "", 500, 0., 10., fNResolutionRelptBins, fResolutionRelPtBins.data());
      fPGen_DeltaEta                       = new TH2D("PGen_DeltaEta",                      "", 500, 0., 10., fNResolutionetaBins, fResolutionEtaBins.data());
      fPtGen_DeltaEta                      = new TH2D("PtGen_DeltaEta",                     "", 500, 0., 10., fNResolutionetaBins, fResolutionEtaBins.data());
      fPGen_DeltaTheta                     = new TH2D("PGen_DeltaTheta",                    "", 500, 0., 10., fNResolutionthetaBins, fResolutionThetaBins.data());
      fPGen_DeltaPhi_Ele                   = new TH2D("PGen_DeltaPhi_Ele",                  "", 500, 0., 10., fNResolutionphiBins, fResolutionPhiBins.data());
      fPGen_DeltaPhi_Pos                   = new TH2D("PGen_DeltaPhi_Pos",                  "", 500, 0., 10., fNResolutionphiBins, fResolutionPhiBins.data());
      fPtGen_DeltaPhi_Ele                  = new TH2D("PtGen_DeltaPhi_Ele",                 "", 500, 0., 10., fNResolutionphiBins, fResolutionPhiBins.data());
      fPtGen_DeltaPhi_Pos                  = new TH2D("PtGen_DeltaPhi_Pos",                 "", 500, 0., 10., fNResolutionphiBins, fResolutionPhiBins.data());
      fThetaGen_DeltaTheta                 = new TH2D("ThetaGen_DeltaTheta",                "", 220, -0.1*TMath::Pi(), 1.1*TMath::Pi(), fNResolutionthetaBins, fResolutionThetaBins.data());
      fPhiGen_DeltaPhi                     = new TH2D("PhiGen_DeltaPhi",                    "", 320, -0.1*TMath::Pi(), 2.1*TMath::Pi(), fNResolutionphiBins, fResolutionPhiBins.data());

      fPGen_DeltaP                         ->Sumw2();
      fPtGen_DeltaPt                       ->Sumw2();
      fPtGen_DeltaPtOverPtGen              ->Sumw2();
      fPGen_PrecOverPGen                   ->Sumw2();
      fPtGen_PtRecOverPtGen                ->Sumw2();
      fPGen_DeltaEta                       ->Sumw2();
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
      fPtGen_DeltaPt                       ->GetXaxis()->SetTitle("p^{gen}_{T} (GeV/c)");
      fPtGen_DeltaPt                       ->GetYaxis()->SetTitle("p^{rec}_{T} - p^{gen}_{T} (GeV/c)");
      fPtGen_DeltaPtOverPtGen              ->GetXaxis()->SetTitle("p^{gen}_{T} (GeV/c)");
      fPtGen_DeltaPtOverPtGen              ->GetYaxis()->SetTitle("(p^{gen}_{T} - p^{rec}_{T}) / p^{gen}_{T} (GeV/c)");
      fPGen_PrecOverPGen                   ->GetXaxis()->SetTitle("p^{gen} (GeV/c)");
      fPGen_PrecOverPGen                   ->GetYaxis()->SetTitle("p^{rec} / p^{gen} (GeV/c)");
      fPtGen_PtRecOverPtGen                ->GetXaxis()->SetTitle("p^{gen}_{T} (GeV/c)");
      fPtGen_PtRecOverPtGen                ->GetYaxis()->SetTitle("p^{rec}_{T} / p^{gen}_{T} (GeV/c)");
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
      fResolutionList->Add(fPtGen_DeltaPt);
      fResolutionList->Add(fPtGen_DeltaPtOverPtGen);
      fResolutionList->Add(fPGen_PrecOverPGen);
      fResolutionList->Add(fPtGen_PtRecOverPtGen);
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
  //  TString genname;
  //   Bool_t yesno=mcEvent->GetCocktailGenerator(it,genname);
  //   if(!yesno) Printf("no cocktail header list was found for this event");
  //   if(yesno) {Printf("cocktail header name is %s", genname.Data());
  //   //you may want to check wether it is Hijing, for example.
  //   if(genname.Contains("Hijing")) Printf("this particle comes from HIJING");}
  //   }

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
  AliMCEvent* fMC = eventHandlerMC->MCEvent();
  if (!fMC) { Printf("ERROR: fMC not available"); return; }

  if (!fPIDResponse) SetPIDResponse( eventHandler->GetPIDResponse() );
  AliDielectronVarManager::SetPIDResponse(fPIDResponse);

  if (isAOD) fEvent = static_cast<AliAODEvent*>(eventHandler->GetEvent());
  else       fEvent = static_cast<AliESDEvent*>(eventHandler->GetEvent());

  AliDielectronVarManager::SetEvent(fEvent);

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
  if (multSelection) centralityF  = multSelection->GetMultiplicityPercentile("V0M",kFALSE);
  if (centralityF == -1 && fMaxCentrality == -1 && fMinCentrality == -1) {/*do nothing*/} // is used for pp and pPb analysis
  else if (centralityF > fMaxCentrality || centralityF < fMinCentrality) { return;} // reject event

  fHistEventStat->Fill(kCentralityEvents);
  fHistEvents->Fill(0.5);
  fHistCentrality->Fill(centralityF);

  // Calculating the weight when centrality correction is applied
  double centralityWeight = 1.;
  if (fHistCentralityCorrection != 0x0){
    centralityWeight = (fHistCentralityCorrection->GetEntries() / fHistCentralityCorrection->GetNbinsX()) / fHistCentralityCorrection->FindBin(centralityF) ;
    std::cout << "cent: " << centralityF << "  " << "weight: " << centralityWeight << std::endl;
  }

  // ##########################################################
  // Fill Multiplicity histogram
  int nTracks = fEvent->GetNumberOfTracks();
  fHistNTracks->Fill(nTracks);


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
    // Creating particles to summarize all the data
    Particle part = CreateParticle(mcPart1);
    part.isMCSignal = mcSignal_acc;
    part.SetTrackID(iPart);

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
      part.fPhi_smeared = smearedVec.Phi()+ pi;

      for (unsigned int i = 0; i < part.isMCSignal.size(); ++i){
        if (part.isMCSignal[i]) {
          if      (part.fCharge < 0){
            dynamic_cast<TH3D*>(fHistGenSmearedNegPart.at(i))->Fill(smearedVec.Pt(), smearedVec.Eta(), smearedVec.Phi() + pi, centralityWeight);
          }
          else if (part.fCharge > 0) {
            dynamic_cast<TH3D*>(fHistGenSmearedPosPart.at(i))->Fill(smearedVec.Pt(), smearedVec.Eta(), smearedVec.Phi() + pi, centralityWeight);
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

    // check if at least one mc signal is true otherwise skip this particle
    if (CheckIfOneIsTrue(mcSignal_acc) == kFALSE) continue;


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
    Particle part = CreateParticle(track);
    part.isMCSignal = mcSignal_acc;
    part.isReconstructed = selected;
    part.SetTrackID(iTracks);


    // ##########################################################
    if      (fDoPairing == true && part.fCharge < 0) fRecNegPart.push_back(part);
    else if (fDoPairing == true && part.fCharge > 0) fRecPosPart.push_back(part);


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

          }// is selected by cutsetting
        } // is selected by MC signal
      } // end of loop over all cutsettings
    } // end of loop over all MCsignals

    // ##########################################################
    // Fill support histograms with first cutsetting and first mcsignal
    if(part.isMCSignal[fSupportMCSignal] == true && part.isReconstructed[fSupportCutsetting] == kTRUE){

      FillTrackHistograms(track); // Fill support histograms


      // ##########################################################
      // Fill resolution histograms
      // ##########################################################
      AliVParticle* mcPart1 = fMC->GetTrack(abslabel);
      double mcP = mcPart1->P();
      double mcPt = mcPart1->Pt();
      double mcEta = mcPart1->Eta();
      double mcPhi = mcPart1->Phi();
      double mcTheta = mcPart1->Theta();

      if(TMath::Abs(mcEta) < 1.0) {
        // fPGen                ->Fill(mcP);
        // fPRec                ->Fill(recP);
        fPGen_DeltaP           ->Fill(mcP,  track->P() - mcP, centralityWeight);
        fPtGen_DeltaPt         ->Fill(mcPt, part.fPt - mcPt, centralityWeight);
        fPtGen_DeltaPtOverPtGen->Fill(mcPt, (mcPt - part.fPt) / mcPt, centralityWeight);
        fPGen_PrecOverPGen     ->Fill(mcP,  track->P() / mcP, centralityWeight);
        fPtGen_PtRecOverPtGen  ->Fill(mcPt, part.fPt / mcPt, centralityWeight);
        if (part.fCharge<0) {
          fPGen_DeltaPhi_Ele   ->Fill(mcP,  part.fPhi - mcPhi, centralityWeight);
          fPtGen_DeltaPhi_Ele  ->Fill(mcPt, mcPhi - part.fPhi, centralityWeight);
        }
        else {
          fPGen_DeltaPhi_Pos   ->Fill(mcP,  part.fPhi - mcPhi, centralityWeight);
          fPtGen_DeltaPhi_Pos  ->Fill(mcPt, mcPhi - part.fPhi, centralityWeight);
        }
        fPhiGen_DeltaPhi       ->Fill(mcPhi, part.fPhi - mcPhi, centralityWeight);
      }

      if(mcPt > fPtMinGen){
        fPGen_DeltaEta       ->Fill(mcP,     part.fEta - mcEta, centralityWeight);
        fPtGen_DeltaEta      ->Fill(mcPt,    part.fEta - mcEta, centralityWeight);
        fPGen_DeltaTheta     ->Fill(mcP,     track->Theta() - mcTheta, centralityWeight);
        fThetaGen_DeltaTheta ->Fill(mcTheta, track->Theta() - mcTheta, centralityWeight);
      }
    }
  }


  // ##########################################################
  // ##########################################################
  // ##########################################################
  // DO PAIRING
  // ##########################################################

  if (fDoPairing){
    for (unsigned int neg_i = 0; neg_i < fGenNegPart.size(); ++neg_i){
      for (unsigned int pos_i = 0; pos_i < fGenPosPart.size(); ++pos_i){
        AliVParticle* mcPart1 = fMC->GetTrack(fGenNegPart[neg_i].GetTrackID());
        AliVParticle* mcPart2 = fMC->GetTrack(fGenPosPart[pos_i].GetTrackID());

        if (fDoULSandLS){
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
              if (fGenNegPart[neg_i].isMCSignal[iMCSignal] == true && fGenPosPart[pos_i].isMCSignal[iMCSignal] == true)
               fHistGenPair_ULSandLS.at(3*iMCSignal)->Fill(mass, pairpt, weight);
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
              if (fGenNegPart[neg_i].isMCSignal[iMCSignal] == true && fGenPosPart[pos_i].isMCSignal[iMCSignal] == true)
               fHistGenSmearedPair_ULSandLS.at(3*iMCSignal)->Fill(mass, pairpt, weight);

            }
          }
        }

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

          for (unsigned int i = 0; i < mcSignal_acc.size(); ++i){
            if (mcSignal_acc[i] == kTRUE){
              fHistGenPair.at(i)->Fill(mass, pairpt, weight);
            }
          } // end of loop over all MCsignals
        }
        if (fArrResoPt){ // Smear particles to fill "GeneratedSmeared"

          if (fGenNegPart[neg_i].fPt_smeared < fPtMin || fGenNegPart[neg_i].fPt_smeared > fPtMax || fGenNegPart[neg_i].fEta_smeared < fEtaMin || fGenNegPart[neg_i].fEta_smeared > fEtaMax) continue;
          if (fGenPosPart[pos_i].fPt_smeared < fPtMin || fGenPosPart[pos_i].fPt_smeared > fPtMax || fGenPosPart[pos_i].fEta_smeared < fEtaMin || fGenPosPart[pos_i].fEta_smeared > fEtaMax) continue;

          // Construct pair variables from LorentzVectors
          TLorentzVector Lvec1;
          TLorentzVector Lvec2;
          Lvec1.SetPtEtaPhiM(fGenNegPart[neg_i].fPt_smeared, fGenNegPart[neg_i].fEta_smeared, fGenNegPart[neg_i].fPhi_smeared, AliPID::ParticleMass(AliPID::kElectron));
          Lvec2.SetPtEtaPhiM(fGenPosPart[pos_i].fPt_smeared, fGenPosPart[pos_i].fEta_smeared, fGenPosPart[pos_i].fPhi_smeared, AliPID::ParticleMass(AliPID::kElectron));
          TLorentzVector LvecM = Lvec1 + Lvec2;
          double mass = LvecM.M();
          double pairpt = LvecM.Pt();

          for (unsigned int i = 0; i < mcSignal_acc.size(); ++i){
            if (mcSignal_acc[i] == kTRUE){
              fHistGenSmearedPair.at(i)->Fill(mass, pairpt, centralityWeight);
            }
          } // end of loop over all MCsignals
        } // end of smearing
      } // end of loop over all positive particles
    } // end of loop over all negative particles

    if (fDoULSandLS){
      // Calculated for single leg signals LS-- pairs
      for (unsigned int neg_i = 0; neg_i < fGenNegPart.size(); ++neg_i){
        for (unsigned int neg_j = neg_i + 1; neg_j < fGenNegPart.size(); ++neg_j){
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
              if (fGenNegPart[neg_i].isMCSignal[iMCSignal] == true && fGenNegPart[neg_j].isMCSignal[iMCSignal] == true)
               fHistGenPair_ULSandLS.at(3*iMCSignal+2)->Fill(mass, pairpt, weight);
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
              if (fGenNegPart[neg_i].isMCSignal[iMCSignal] == true && fGenNegPart[neg_j].isMCSignal[iMCSignal] == true)
               fHistGenSmearedPair_ULSandLS.at(3*iMCSignal+2)->Fill(mass, pairpt, weight);

            }
          }
        }
      }

      // Calculated for single leg signals LS++ pairs
      for (unsigned int pos_i = 0; pos_i < fGenPosPart.size(); ++pos_i){
        for (unsigned int pos_j = pos_i + 1; pos_j < fGenPosPart.size(); ++pos_j){
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
              if (fGenPosPart[pos_i].isMCSignal[iMCSignal] == true && fGenPosPart[pos_j].isMCSignal[iMCSignal] == true)
               fHistGenPair_ULSandLS.at(3*iMCSignal+1)->Fill(mass, pairpt, weight);
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
              if (fGenPosPart[pos_i].isMCSignal[iMCSignal] == true && fGenPosPart[pos_j].isMCSignal[iMCSignal] == true)
               fHistGenSmearedPair_ULSandLS.at(3*iMCSignal+1)->Fill(mass, pairpt, weight);

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

        if (fDoULSandLS){
          // Calculated for single leg signals unlike sign
          TLorentzVector Lvec1;
          TLorentzVector Lvec2;
          Lvec1.SetPtEtaPhiM(fRecNegPart[neg_i].fPt, fRecNegPart[neg_i].fEta, fRecNegPart[neg_i].fPhi, AliPID::ParticleMass(AliPID::kElectron));
          Lvec2.SetPtEtaPhiM(fRecPosPart[pos_i].fPt, fRecPosPart[pos_i].fEta, fRecPosPart[pos_i].fPhi, AliPID::ParticleMass(AliPID::kElectron));
          TLorentzVector LvecM = Lvec1 + Lvec2;
          double mass = LvecM.M();
          double pairpt = LvecM.Pt();

          for (unsigned int i = 0; i < fRecNegPart[neg_i].isMCSignal.size(); ++i){
            if (fRecNegPart[neg_i].isMCSignal[i] == kTRUE && fRecPosPart[pos_i].isMCSignal[i] == kTRUE){
              for (unsigned int j = 0; j < fRecNegPart[neg_i].isReconstructed.size(); ++j){
                if (fRecNegPart[neg_i].isReconstructed[j] == kTRUE && fRecPosPart[pos_i].isReconstructed[j] == kTRUE){
                  fHistRecPair_ULSandLS[j * 3 * fSingleLegMCSignal.size() + 3 * i]->Fill(mass, pairpt, centralityWeight);
                }// is selected by cutsetting
              } // end of loop over all cutsettings
            } // is selected by MC Signal
          } // end of loop over all MCsignals
        }

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

        // ##########################################################
        // Filling reconstructed particle histograms according to MCSignals
        for (unsigned int i = 0; i < mcSignal_acc.size(); ++i){
          if (mcSignal_acc[i] == kTRUE){
            for (unsigned int j = 0; j < fRecNegPart[neg_i].isReconstructed.size(); ++j){
              if (fRecNegPart[neg_i].isReconstructed[j] == kTRUE && fRecPosPart[pos_i].isReconstructed[j] == kTRUE){
                fHistRecPair.at(j * mcSignal_acc.size() + i)->Fill(mass, pairpt, centralityWeight);
              }// is selected by cutsetting
            } // end of loop over all cutsettings
          } // is selected by MCSignal
        } // end of loop over all MCsignals

      }// end of positive particle loop
    } // end of negative particle loop
  } // End of pairing

  if (fDoULSandLS){
    //LS--
    for (unsigned int neg_i = 0; neg_i < fRecNegPart.size(); ++neg_i){
      for (unsigned int neg_j = neg_i + 1; neg_j < fRecNegPart.size(); ++neg_j){
        // Calculated for single leg signals like sign
        TLorentzVector Lvec1;
        TLorentzVector Lvec2;
        Lvec1.SetPtEtaPhiM(fRecNegPart[neg_i].fPt, fRecNegPart[neg_i].fEta, fRecNegPart[neg_i].fPhi, AliPID::ParticleMass(AliPID::kElectron));
        Lvec2.SetPtEtaPhiM(fRecNegPart[neg_j].fPt, fRecNegPart[neg_j].fEta, fRecNegPart[neg_j].fPhi, AliPID::ParticleMass(AliPID::kElectron));
        TLorentzVector LvecM = Lvec1 + Lvec2;
        double mass = LvecM.M();
        double pairpt = LvecM.Pt();

        for (unsigned int i = 0; i < fRecNegPart[neg_i].isMCSignal.size(); ++i){
          if (fRecNegPart[neg_i].isMCSignal[i] == kTRUE && fRecNegPart[neg_j].isMCSignal[i] == kTRUE){
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
        // Calculated for single leg signals like sign
        TLorentzVector Lvec1;
        TLorentzVector Lvec2;
        Lvec1.SetPtEtaPhiM(fRecPosPart[pos_i].fPt, fRecPosPart[pos_i].fEta, fRecPosPart[pos_i].fPhi, AliPID::ParticleMass(AliPID::kElectron));
        Lvec2.SetPtEtaPhiM(fRecPosPart[pos_j].fPt, fRecPosPart[pos_j].fEta, fRecPosPart[pos_j].fPhi, AliPID::ParticleMass(AliPID::kElectron));
        TLorentzVector LvecM = Lvec1 + Lvec2;
        double mass = LvecM.M();
        double pairpt = LvecM.Pt();

        for (unsigned int i = 0; i < fRecPosPart[pos_i].isMCSignal.size(); ++i){
          if (fRecPosPart[pos_i].isMCSignal[i] == kTRUE && fRecPosPart[pos_j].isMCSignal[i] == kTRUE){
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
void    AliAnalysisTaskElectronEfficiencyV2::FillTrackHistograms(AliVParticle* track){
  Double_t values[AliDielectronVarManager::kNMaxValues]={0.};
  AliDielectronVarManager::SetFillMap(fUsedVars); // currently filled manually in the constructor of this task.
  AliDielectronVarManager::Fill(track, values);
  // std::cout << "pt var manager = " << values[AliDielectronVarManager::kPt] << std::endl;
  // std::cout << fOutputListSupportHistos << std::endl;
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
}


// ############################################################################
// ############################################################################
AliAnalysisTaskElectronEfficiencyV2::Particle AliAnalysisTaskElectronEfficiencyV2::CreateParticle(AliVParticle* mcPart1){
  double  pt1      = mcPart1->Pt();
  double  eta1     = mcPart1->Eta();
  double  phi1     = mcPart1->Phi();
  short   charge1  = mcPart1->Charge();
  Particle part(pt1, eta1, phi1, charge1);
  return part;
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

  const double stepSize = (max - min) / steps;
  for (unsigned int i = 0; i < steps; ++i){
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
  TH1D* hPt      = new TH1D("Pt","Pt;Pt [GeV];#tracks",160,0.,8.);//,AliDielectronVarManager::kPt);
  fOutputListSupportHistos->AddAt(hPt,     0);

  // PID
  TH2D* hITSnSigmaEle_P = new TH2D("ITSnSigmaEle_P","ITS number of sigmas Electrons;P [GeV/c];ITS number of sigmas Electrons", 160,0.,8.,100,-5.,5.);//.,AliDielectronVarManager::kP,AliDielectronVarManager::kITSnSigmaEle,makeLogx);
  fOutputListSupportHistos->AddAt(hITSnSigmaEle_P, 1);
  TH2D* hTPCnSigmaEle_P = new TH2D("TPCnSigmaEle_P","TPC number of sigmas Electrons;PIn (pTPC) [GeV/c];TPC number of sigmas Electrons", 160,0.,8.,100,-5.,5.);//.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCnSigmaEle,makeLogx);
  fOutputListSupportHistos->AddAt(hTPCnSigmaEle_P, 2);
  TH2D* hTOFnSigmaEle_P = new TH2D("TOFnSigmaEle_P","TOF number of sigmas Electrons;PIn (pTPC) [GeV/c];TOF number of sigmas Electrons", 160,0.,8.,100,-5.,5.);//,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTOFnSigmaEle,makeLogx);
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
                                      160,0.,8.,160,-0.5,159.5);//,AliDielectronVarManager::kPt,AliDielectronVarManager::kNFclsTPCr);
  fOutputListSupportHistos->AddAt(hTPCcrossedRows_TPCnCls, 19);
  fOutputListSupportHistos->AddAt(hTPCcrossedRows_Pt, 20);


}


// ############################################################################
// ############################################################################
TLorentzVector AliAnalysisTaskElectronEfficiencyV2::ApplyResolution(double pt, double eta, double phi, short ch) {
  // from Theos LightFlavorGenerator, modified

  TLorentzVector resvec;
  resvec.SetPtEtaPhiM(pt, eta, phi, AliPID::ParticleMass(AliPID::kElectron));

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
