/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Authors: Friederike Bock, Daniel Muehlheim                             *
 * Version 1.0                                                            *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

////////////////////////////////////////////////
//---------------------------------------------
// Class handling all kinds of selection cuts for
// Photon from EMCAL clusters
//---------------------------------------------
////////////////////////////////////////////////

#include "AliCaloPhotonCuts.h"
#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"
#include "AliMCEventHandler.h"
#include "AliAODHandler.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "AliStack.h"
#include "AliAODConversionMother.h"
#include "AliAODConversionPhoton.h"
#include "TObjString.h"
#include "AliAODEvent.h"
#include "AliESDEvent.h"
#include "AliCentrality.h"
#include "TList.h"
#include "TFile.h"
#include "AliLog.h"
#include "AliV0ReaderV1.h"
#include "AliAODMCParticle.h"
#include "AliAODMCHeader.h"
#include "AliPicoTrack.h"
#include "AliPHOSGeoUtils.h"
#include "AliTrackerBase.h"
#include "AliVCaloCells.h"
#include "AliVCluster.h"
#include "AliTender.h"
#include "AliTenderSupply.h"
#include "AliEMCALTenderSupply.h"
#include "AliEmcalTenderTask.h"
#include "AliPHOSTenderSupply.h"
#include "AliOADBContainer.h"
#include "AliESDtrackCuts.h"
#include "AliCaloTrackMatcher.h"
#include <vector>

class iostream;

using namespace std;

ClassImp(AliCaloPhotonCuts)


const char* AliCaloPhotonCuts::fgkCutNames[AliCaloPhotonCuts::kNCuts] = {
  "ClusterType",          //0   0: all,    1: EMCAL,   2: PHOS
  "EtaMin",               //1   0: -10,    1: -0.6687, 2: -0,5, 3: -2
  "EtaMax",               //2   0: 10,     1: 0.66465, 2: 0.5,  3: 2
  "PhiMin",               //3   0: -10000, 1: 1.39626
  "PhiMax",               //4   0: 10000, 1: 3.125
  "NonLinearity1"         //5
  "NonLinearity2"         //6
  "DistanceToBadChannel", //7   0: 0,      1: 5
  "Timing",               //8   0: no cut
  "TrackMatching",        //9   0: 0,      1: 5
  "ExoticCluster",        //10   0: no cut
  "MinEnergy",            //11   0: no cut, 1: 0.05,    2: 0.1,  3: 0.15, 4: 0.2, 5: 0.3, 6: 0.5, 7: 0.75, 8: 1, 9: 1.25 (all GeV)
  "MinNCells",            //12  0: no cut, 1: 1,       2: 2,    3: 3,    4: 4,   5: 5,   6: 6
  "MinM02",               //13
  "MaxM02",               //14
  "MinM20",               //15
  "MaxM20",               //16
  "MaximumDispersion",    //17
  "NLM"                   //18
};


//________________________________________________________________________
AliCaloPhotonCuts::AliCaloPhotonCuts(Int_t isMC, const char *name,const char *title) :
  AliAnalysisCuts(name,title),
  fHistograms(NULL),
  fHistExtQA(NULL),
  fCaloTrackMatcher(NULL),
  fGeomEMCAL(NULL),
  fEMCALRecUtils(NULL),
  fEMCALInitialized(kFALSE),
  fGeomPHOS(NULL),
  fPHOSInitialized(kFALSE),
  fPHOSCurrentRun(-1),
  fEMCALBadChannelsMap(NULL),
  fPHOSBadChannelsMap(NULL),
  fBadChannels(NULL),
  fNMaxEMCalModules(12),
  fNMaxPHOSModules(5),
  fHistoModifyAcc(NULL),
  fDoLightOutput(kFALSE),
  fIsMC(0),
  fIsCurrentClusterAcceptedBeforeTM(kFALSE),
  fV0ReaderName("V0ReaderV1"),
  fCaloTrackMatcherName("CaloTrackMatcher_1"),
  fPeriodName(""),
  fCurrentMC(kNoMC),
  fClusterType(0),
  fMinEtaCut(-10),
  fMaxEtaCut(10),
  fUseEtaCut(0),
  fMinPhiCut(-10000),
  fMaxPhiCut(-10000),
  fUsePhiCut(0),
  fMinDistanceToBadChannel(0),
  fUseDistanceToBadChannel(0),
  fMaxTimeDiff(10e10),
  fMinTimeDiff(-10e10),
  fUseTimeDiff(0),
  fMaxDistTrackToClusterEta(0),
  fMinDistTrackToClusterPhi(0),
  fMaxDistTrackToClusterPhi(0),
  fUseDistTrackToCluster(0),
  fUsePtDepTrackToCluster(0),
  fFuncPtDepEta(0),
  fFuncPtDepPhi(0),
  fExtendedMatchAndQA(0),
  fExoticEnergyFracCluster(0),
  fExoticMinEnergyCell(1),
  fUseExoticCluster(0),
  fDoExoticsQA(kFALSE),
  fMinEnergy(0),
  fSeedEnergy(0.1),
  fLocMaxCutEDiff(0.03),
  fUseMinEnergy(0),
  fMinNCells(0),
  fUseNCells(0),
  fMaxM02(1000),
  fMinM02(0),
  fUseM02(0),
  fMaxM02CutNr(0),
  fMinM02CutNr(0),
  fMaxM20(1000),
  fMinM20(0),
  fUseM20(0),
  fMaxDispersion(1000),
  fUseDispersion(0),
  fMinNLM(0),
  fMaxNLM(1000),
  fUseNLM(0),
  fNonLinearity1(0),
  fNonLinearity2(0),
  fSwitchNonLinearity(0),
  fUseNonLinearity(kFALSE),
  fIsPureCalo(0),
  fVectorMatchedClusterIDs(0),
  fCutString(NULL),
  fHistCutIndex(NULL),
  fHistAcceptanceCuts(NULL),
  fHistClusterIdentificationCuts(NULL),
  fHistClusterEtavsPhiBeforeAcc(NULL),
  fHistClusterEtavsPhiAfterAcc(NULL),
  fHistClusterEtavsPhiAfterQA(NULL),
  fHistClusterTimevsEBeforeQA(NULL),
  fHistClusterTimevsEAfterQA(NULL),
  fHistEnergyOfClusterBeforeNL(NULL),
  fHistEnergyOfClusterAfterNL(NULL),
  fHistEnergyOfClusterBeforeQA(NULL),
  fHistEnergyOfClusterAfterQA(NULL),
  fHistNCellsBeforeQA(NULL),
  fHistNCellsAfterQA(NULL),
  fHistNLMVsNCellsAfterQA(NULL),
  fHistNLMVsEAfterQA(NULL),
  fHistM02BeforeQA(NULL),
  fHistM02AfterQA(NULL),
  fHistM20BeforeQA(NULL),
  fHistM20AfterQA(NULL),
  fHistDispersionBeforeQA(NULL),
  fHistDispersionAfterQA(NULL),
  fHistNLMBeforeQA(NULL),
  fHistNLMAfterQA(NULL),
//   fHistNLMAvsNLMBBeforeQA(NULL),
  fHistClusterEnergyvsMod(NULL),
  fHistNCellsBigger100MeVvsMod(NULL),
  fHistNCellsBigger1500MeVvsMod(NULL),
  fHistEnergyOfModvsMod(NULL),
  fHistClusterEnergyvsNCells(NULL),
  fHistCellEnergyvsCellID(NULL),
  fHistCellTimevsCellID(NULL),
  fHistClusterEM02BeforeQA(NULL),
  fHistClusterEM02AfterQA(NULL),
  fHistClusterIncludedCellsBeforeQA(NULL),
  fHistClusterIncludedCellsAfterQA(NULL),
  fHistClusterEnergyFracCellsBeforeQA(NULL),
  fHistClusterEnergyFracCellsAfterQA(NULL),
  fHistClusterIncludedCellsTimingAfterQA(NULL),
  fHistClusterIncludedCellsTimingEnergyAfterQA(NULL),
  fHistClusterDistanceInTimeCut(NULL),
  fHistClusterDistanceOutTimeCut(NULL),
  fHistClusterDistance1DInTimeCut(NULL),
  fHistClusterRBeforeQA(NULL),
  fHistClusterRAfterQA(NULL),
  fHistClusterdEtadPhiBeforeQA(NULL),
  fHistClusterdEtadPhiAfterQA(NULL),
  fHistDistanceTrackToClusterBeforeQA(NULL),
  fHistDistanceTrackToClusterAfterQA(NULL),
  fHistClusterdEtadPhiPosTracksBeforeQA(NULL),
  fHistClusterdEtadPhiNegTracksBeforeQA(NULL),
  fHistClusterdEtadPhiPosTracksAfterQA(NULL),
  fHistClusterdEtadPhiNegTracksAfterQA(NULL),
  fHistClusterdEtadPhiPosTracksP_000_075BeforeQA(NULL),
  fHistClusterdEtadPhiPosTracksP_075_125BeforeQA(NULL),
  fHistClusterdEtadPhiPosTracksP_125_999BeforeQA(NULL),
  fHistClusterdEtadPhiNegTracksP_000_075BeforeQA(NULL),
  fHistClusterdEtadPhiNegTracksP_075_125BeforeQA(NULL),
  fHistClusterdEtadPhiNegTracksP_125_999BeforeQA(NULL),
  fHistClusterdEtadPtBeforeQA(NULL),
  fHistClusterdEtadPtAfterQA(NULL),
  fHistClusterdEtadPtTrueMatched(NULL),
  fHistClusterdPhidPtPosTracksBeforeQA(NULL),
  fHistClusterdPhidPtNegTracksBeforeQA(NULL),
  fHistClusterdPhidPtAfterQA(NULL),
  fHistClusterdPhidPtPosTracksTrueMatched(NULL),
  fHistClusterdPhidPtNegTracksTrueMatched(NULL),
  fHistClusterM20M02BeforeQA(NULL),
  fHistClusterM20M02AfterQA(NULL),
  fHistClusterEtavsPhiExotics(NULL),
  fHistClusterEM02Exotics(NULL),
  fHistClusterEnergyvsNCellsExotics(NULL),
  fHistClusterEEstarExotics(NULL),
  fHistClusterTMEffiInput(NULL),
  fHistClusterEvsTrackECharged(NULL),
  fHistClusterEvsTrackEChargedLead(NULL),
  fHistClusterEvsTrackENeutral(NULL),
  fHistClusterEvsTrackENeutralSubCharged(NULL),
  fHistClusterEvsTrackEGamma(NULL),
  fHistClusterEvsTrackEGammaSubCharged(NULL),
  fHistClusterEvsTrackEConv(NULL),
  fHistClusterENMatchesNeutral(NULL),
  fHistClusterENMatchesCharged(NULL),
  fHistClusterEvsTrackEPrimaryButNoElec(NULL),
  fHistClusterEvsTrackSumEPrimaryButNoElec(NULL),
  fNMaxDCalModules(8),
  fgkDCALCols(32)    
{
  for(Int_t jj=0;jj<kNCuts;jj++){fCuts[jj]=0;}
  fCutString=new TObjString((GetCutNumber()).Data());

  fIsMC = isMC;
}

//________________________________________________________________________
AliCaloPhotonCuts::AliCaloPhotonCuts(const AliCaloPhotonCuts &ref) :
  AliAnalysisCuts(ref),
  fHistograms(NULL),
  fHistExtQA(NULL),
  fCaloTrackMatcher(NULL),
  fGeomEMCAL(NULL),
  fEMCALRecUtils(NULL),
  fEMCALInitialized(kFALSE),
  fGeomPHOS(NULL),
  fPHOSInitialized(kFALSE),
  fPHOSCurrentRun(-1),
  fEMCALBadChannelsMap(NULL),
  fPHOSBadChannelsMap(NULL),
  fBadChannels(NULL),
  fNMaxEMCalModules(ref.fNMaxEMCalModules),
  fNMaxPHOSModules(ref.fNMaxPHOSModules),
  fHistoModifyAcc(NULL),
  fDoLightOutput(ref.fDoLightOutput),
  fIsMC(ref.fIsMC),
  fIsCurrentClusterAcceptedBeforeTM(kFALSE),
  fV0ReaderName(ref.fV0ReaderName),
  fCaloTrackMatcherName(ref.fCaloTrackMatcherName),
  fPeriodName(ref.fPeriodName),
  fCurrentMC(ref.fCurrentMC),
  fClusterType(ref.fClusterType),
  fMinEtaCut(ref.fMinEtaCut),
  fMaxEtaCut(ref.fMaxEtaCut),
  fUseEtaCut(ref.fUseEtaCut),
  fMinPhiCut(ref.fMinPhiCut),
  fMaxPhiCut(ref.fMaxPhiCut),
  fUsePhiCut(ref.fUsePhiCut),
  fMinDistanceToBadChannel(ref.fMinDistanceToBadChannel),
  fUseDistanceToBadChannel(ref.fUseDistanceToBadChannel),
  fMaxTimeDiff(ref.fMaxTimeDiff),
  fMinTimeDiff(ref.fMinTimeDiff),
  fUseTimeDiff(ref.fUseTimeDiff),
  fMaxDistTrackToClusterEta(ref.fMaxDistTrackToClusterEta),
  fMinDistTrackToClusterPhi(ref.fMinDistTrackToClusterPhi),
  fMaxDistTrackToClusterPhi(ref.fMaxDistTrackToClusterPhi),
  fUseDistTrackToCluster(ref.fUseDistTrackToCluster),
  fUsePtDepTrackToCluster(ref.fUsePtDepTrackToCluster),
  fFuncPtDepEta(ref.fFuncPtDepEta),
  fFuncPtDepPhi(ref.fFuncPtDepPhi),
  fExtendedMatchAndQA(ref.fExtendedMatchAndQA),
  fExoticEnergyFracCluster(ref.fExoticEnergyFracCluster),
  fExoticMinEnergyCell(ref.fExoticMinEnergyCell),
  fUseExoticCluster(ref.fUseExoticCluster),
  fDoExoticsQA(ref.fDoExoticsQA),
  fMinEnergy(ref.fMinEnergy),
  fSeedEnergy(ref.fSeedEnergy),
  fLocMaxCutEDiff(ref.fLocMaxCutEDiff),
  fUseMinEnergy(ref.fUseMinEnergy),
  fMinNCells(ref.fMinNCells),
  fUseNCells(ref.fUseNCells),
  fMaxM02(ref.fMaxM02),
  fMinM02(ref.fMinM02),
  fUseM02(ref.fUseM02),
  fMaxM02CutNr(ref.fMaxM02CutNr),
  fMinM02CutNr(ref.fMinM02CutNr),
  fMaxM20(ref.fMaxM20),
  fMinM20(ref.fMinM20),
  fUseM20(ref.fUseDispersion),
  fMaxDispersion(ref.fMaxDispersion),
  fUseDispersion(ref.fUseDispersion),
  fMinNLM(ref.fMinNLM),
  fMaxNLM(ref.fMaxNLM),
  fUseNLM(ref.fUseNLM),
  fNonLinearity1(ref.fNonLinearity1),
  fNonLinearity2(ref.fNonLinearity2),
  fSwitchNonLinearity(ref.fSwitchNonLinearity),
  fUseNonLinearity(ref.fUseNonLinearity),
  fIsPureCalo(ref.fIsPureCalo),
  fVectorMatchedClusterIDs(0),
  fCutString(NULL),
  fHistCutIndex(NULL),
  fHistAcceptanceCuts(NULL),
  fHistClusterIdentificationCuts(NULL),
  fHistClusterEtavsPhiBeforeAcc(NULL),
  fHistClusterEtavsPhiAfterAcc(NULL),
  fHistClusterEtavsPhiAfterQA(NULL),
  fHistClusterTimevsEBeforeQA(NULL),
  fHistClusterTimevsEAfterQA(NULL),
  fHistEnergyOfClusterBeforeNL(NULL),
  fHistEnergyOfClusterAfterNL(NULL),
  fHistEnergyOfClusterBeforeQA(NULL),
  fHistEnergyOfClusterAfterQA(NULL),
  fHistNCellsBeforeQA(NULL),
  fHistNCellsAfterQA(NULL),
  fHistM02BeforeQA(NULL),
  fHistM02AfterQA(NULL),
  fHistM20BeforeQA(NULL),
  fHistM20AfterQA(NULL),
  fHistDispersionBeforeQA(NULL),
  fHistDispersionAfterQA(NULL),
  fHistNLMBeforeQA(NULL),
  fHistNLMAfterQA(NULL),
//   fHistNLMAvsNLMBBeforeQA(NULL),
  fHistNLMVsNCellsAfterQA(NULL),
  fHistNLMVsEAfterQA(NULL),
  fHistClusterEnergyvsMod(NULL),
  fHistNCellsBigger100MeVvsMod(NULL),
  fHistNCellsBigger1500MeVvsMod(NULL),
  fHistEnergyOfModvsMod(NULL),
  fHistClusterEnergyvsNCells(NULL),
  fHistCellEnergyvsCellID(NULL),
  fHistCellTimevsCellID(NULL),
  fHistClusterEM02BeforeQA(NULL),
  fHistClusterEM02AfterQA(NULL),
  fHistClusterIncludedCellsBeforeQA(NULL),
  fHistClusterIncludedCellsAfterQA(NULL),
  fHistClusterEnergyFracCellsBeforeQA(NULL),
  fHistClusterEnergyFracCellsAfterQA(NULL),
  fHistClusterIncludedCellsTimingAfterQA(NULL),
  fHistClusterIncludedCellsTimingEnergyAfterQA(NULL),
  fHistClusterDistanceInTimeCut(NULL),
  fHistClusterDistanceOutTimeCut(NULL),
  fHistClusterDistance1DInTimeCut(NULL),
  fHistClusterRBeforeQA(NULL),
  fHistClusterRAfterQA(NULL),
  fHistClusterdEtadPhiBeforeQA(NULL),
  fHistClusterdEtadPhiAfterQA(NULL),
  fHistDistanceTrackToClusterBeforeQA(NULL),
  fHistDistanceTrackToClusterAfterQA(NULL),
  fHistClusterdEtadPhiPosTracksBeforeQA(NULL),
  fHistClusterdEtadPhiNegTracksBeforeQA(NULL),
  fHistClusterdEtadPhiPosTracksAfterQA(NULL),
  fHistClusterdEtadPhiNegTracksAfterQA(NULL),
  fHistClusterdEtadPhiPosTracksP_000_075BeforeQA(NULL),
  fHistClusterdEtadPhiPosTracksP_075_125BeforeQA(NULL),
  fHistClusterdEtadPhiPosTracksP_125_999BeforeQA(NULL),
  fHistClusterdEtadPhiNegTracksP_000_075BeforeQA(NULL),
  fHistClusterdEtadPhiNegTracksP_075_125BeforeQA(NULL),
  fHistClusterdEtadPhiNegTracksP_125_999BeforeQA(NULL),
  fHistClusterdEtadPtBeforeQA(NULL),
  fHistClusterdEtadPtAfterQA(NULL),
  fHistClusterdEtadPtTrueMatched(NULL),
  fHistClusterdPhidPtPosTracksBeforeQA(NULL),
  fHistClusterdPhidPtNegTracksBeforeQA(NULL),
  fHistClusterdPhidPtAfterQA(NULL),
  fHistClusterdPhidPtPosTracksTrueMatched(NULL),
  fHistClusterdPhidPtNegTracksTrueMatched(NULL),
  fHistClusterM20M02BeforeQA(NULL),
  fHistClusterM20M02AfterQA(NULL),
  fHistClusterEtavsPhiExotics(NULL),
  fHistClusterEM02Exotics(NULL),
  fHistClusterEnergyvsNCellsExotics(NULL),
  fHistClusterEEstarExotics(NULL),
  fHistClusterTMEffiInput(NULL),
  fHistClusterEvsTrackECharged(NULL),
  fHistClusterEvsTrackEChargedLead(NULL),
  fHistClusterEvsTrackENeutral(NULL),
  fHistClusterEvsTrackENeutralSubCharged(NULL),
  fHistClusterEvsTrackEGamma(NULL),
  fHistClusterEvsTrackEGammaSubCharged(NULL),
  fHistClusterEvsTrackEConv(NULL),
  fHistClusterENMatchesNeutral(NULL),
  fHistClusterENMatchesCharged(NULL),
  fHistClusterEvsTrackEPrimaryButNoElec(NULL),
  fHistClusterEvsTrackSumEPrimaryButNoElec(NULL),
  fNMaxDCalModules(ref.fNMaxDCalModules),
  fgkDCALCols(ref.fgkDCALCols)
{
  // Copy Constructor
  for(Int_t jj=0;jj<kNCuts;jj++){fCuts[jj]=ref.fCuts[jj];}
  fCutString=new TObjString((GetCutNumber()).Data());

}


//________________________________________________________________________
AliCaloPhotonCuts::~AliCaloPhotonCuts() {
  // Destructor
  if(fCutString != NULL){
    delete fCutString;
    fCutString = NULL;
  }

  if(fPHOSBadChannelsMap != NULL){
    delete[] fPHOSBadChannelsMap;
    fPHOSBadChannelsMap = NULL;
  }

  if(fFuncPtDepEta) delete fFuncPtDepEta;
  if(fFuncPtDepPhi) delete fFuncPtDepPhi;
}

//________________________________________________________________________
void AliCaloPhotonCuts::InitCutHistograms(TString name){

  // Initialize Cut Histograms for QA (only initialized and filled if function is called)
  TH1::AddDirectory(kFALSE);

  
  
  if (fDoLightOutput) 
    fDoExoticsQA    = kFALSE;
  
  if(fHistograms != NULL){
    delete fHistograms;
    fHistograms=NULL;
  }
  if(fHistograms==NULL){
    fHistograms     = new TList();
    fHistograms->SetOwner(kTRUE);
    if(name=="")fHistograms->SetName(Form("CaloCuts_%s",GetCutNumber().Data()));
    else fHistograms->SetName(Form("%s_%s",name.Data(),GetCutNumber().Data()));
  }

  if(fHistExtQA != NULL){
    delete fHistExtQA;
    fHistExtQA=NULL;
  }
  if(fHistExtQA==NULL){
    fHistExtQA      = new TList();
    fHistExtQA->SetOwner(kTRUE);
    if(name=="")fHistExtQA->SetName(Form("CaloExtQA_%s",GetCutNumber().Data()));
    else fHistExtQA->SetName(Form("%s_ExtQA_%s",name.Data(),GetCutNumber().Data()));
  }

  Int_t nBinsClusterEFine           = 400;
  Int_t nBinsClusterECoarse         = 100;
  Double_t minClusterELog           = 0.5;
  Double_t maxClusterELog           = 100.0;
  Int_t nBinsCellECoarse            = 100;
  Double_t minCellELog              = 0.05;
  Double_t maxCellELog              = 50.0;
  Int_t nBinsModuleECoarse          = 400;
  Double_t minModuleELog            = 0.1;
  Double_t maxModuleELog            = 400.0;
  
  // IsPhotonSelected
  fHistCutIndex                   = new TH1F(Form("IsPhotonSelected %s",GetCutNumber().Data()),"IsPhotonSelected",5,-0.5,4.5);
  fHistCutIndex->GetXaxis()->SetBinLabel(kPhotonIn+1,"in");
  fHistCutIndex->GetXaxis()->SetBinLabel(kDetector+1,"detector");
  fHistCutIndex->GetXaxis()->SetBinLabel(kAcceptance+1,"acceptance");
  fHistCutIndex->GetXaxis()->SetBinLabel(kClusterQuality+1,"cluster QA");
  fHistCutIndex->GetXaxis()->SetBinLabel(kPhotonOut+1,"out");
  fHistograms->Add(fHistCutIndex);

  // Acceptance Cuts
  fHistAcceptanceCuts             = new TH1F(Form("AcceptanceCuts %s",GetCutNumber().Data()),"AcceptanceCuts",5,-0.5,4.5);
  fHistAcceptanceCuts->GetXaxis()->SetBinLabel(1,"in");
  fHistAcceptanceCuts->GetXaxis()->SetBinLabel(2,"eta");
  fHistAcceptanceCuts->GetXaxis()->SetBinLabel(3,"phi");
  fHistAcceptanceCuts->GetXaxis()->SetBinLabel(4,"dist bad channel/acc map");
  fHistAcceptanceCuts->GetXaxis()->SetBinLabel(5,"out");
  fHistograms->Add(fHistAcceptanceCuts);

  // Cluster Cuts
  fHistClusterIdentificationCuts  = new TH2F(Form("ClusterQualityCuts vs E %s",GetCutNumber().Data()),"ClusterQualityCuts",11,-0.5,10.5,nBinsClusterEFine,minClusterELog,maxClusterELog);
  fHistClusterIdentificationCuts->GetXaxis()->SetBinLabel(1,"in");
  fHistClusterIdentificationCuts->GetXaxis()->SetBinLabel(2,"timing");
  fHistClusterIdentificationCuts->GetXaxis()->SetBinLabel(3,"Exotics");
  fHistClusterIdentificationCuts->GetXaxis()->SetBinLabel(4,"minimum energy");
  fHistClusterIdentificationCuts->GetXaxis()->SetBinLabel(5,"minimum NCells");
  fHistClusterIdentificationCuts->GetXaxis()->SetBinLabel(6,"NLM");
  fHistClusterIdentificationCuts->GetXaxis()->SetBinLabel(7,"M02");
  fHistClusterIdentificationCuts->GetXaxis()->SetBinLabel(8,"M20");
  fHistClusterIdentificationCuts->GetXaxis()->SetBinLabel(9,"dispersion");
  fHistClusterIdentificationCuts->GetXaxis()->SetBinLabel(10,"track matching");
  fHistClusterIdentificationCuts->GetXaxis()->SetBinLabel(11,"out");
  SetLogBinningYTH2(fHistClusterIdentificationCuts);
  fHistograms->Add(fHistClusterIdentificationCuts);

  // Acceptance related histogramms
  if( fClusterType == 1 ){ //EMCAL
    const Int_t nEmcalEtaBins             = 96;
    const Int_t nEmcalPhiBins             = 124;
    Float_t EmcalEtaBins[nEmcalEtaBins+1] = {-0.66687,-0.653,-0.63913,-0.62526,-0.61139,-0.59752,-0.58365,-0.56978,-0.55591,-0.54204,-0.52817,-0.5143,-0.50043,-0.48656,-0.47269,-0.45882,-0.44495,-0.43108,-0.41721,-0.40334,-0.38947,-0.3756,-0.36173,-0.34786,-0.33399,-0.32012,-0.30625,-0.29238,-0.27851,-0.26464,-0.25077,-0.2369,-0.22303,-0.20916,-0.19529,-0.18142,-0.16755,-0.15368,-0.13981,-0.12594,-0.11207,-0.0982,-0.08433,-0.07046,-0.05659,-0.04272,-0.02885,-0.01498,-0.00111,0.01276,0.02663,0.0405,0.05437,0.06824,0.08211,0.09598,0.10985,0.12372,0.13759,0.15146,0.16533,0.1792,0.19307,0.20694,0.22081,0.23468,0.24855,0.26242,0.27629,0.29016,0.30403,0.3179,0.33177,0.34564,0.35951,0.37338,0.38725,0.40112,0.41499,0.42886,0.44273,0.4566,0.47047,0.48434,0.49821,0.51208,0.52595,0.53982,0.55369,0.56756,0.58143,0.5953,0.60917,0.62304,0.63691,0.65078,0.66465};
    Float_t EmcalPhiBins[nEmcalPhiBins+1] = {1.408,1.4215,1.435,1.4485,1.462,1.4755,1.489,1.5025,1.516,1.5295,1.543,1.5565,1.57,1.5835,1.597,1.6105,1.624,1.6375,1.651,1.6645,1.678,1.6915,1.705,1.7185,1.732,1.758,1.7715,1.785,1.7985,1.812,1.8255,1.839,1.8525,1.866,1.8795,1.893,1.9065,1.92,1.9335,1.947,1.9605,1.974,1.9875,2.001,2.0145,2.028,2.0415,2.055,2.0685,2.082,2.108,2.1215,2.135,2.1485,2.162,2.1755,2.189,2.2025,2.216,2.2295,2.243,2.2565,2.27,2.2835,2.297,2.3105,2.324,2.3375,2.351,2.3645,2.378,2.3915,2.405,2.4185,2.432,2.456,2.4695,2.483,2.4965,2.51,2.5235,2.537,2.5505,2.564,2.5775,2.591,2.6045,2.618,2.6315,2.645,2.6585,2.672,2.6855,2.699,2.7125,2.726,2.7395,2.753,2.7665,2.78,2.804,2.8175,2.831,2.8445,2.858,2.8715,2.885,2.8985,2.912,2.9255,2.939,2.9525,2.966,2.9795,2.993,3.0065,3.02,3.0335,3.047,3.0605,3.074,3.0875,3.101,3.1145,3.128};

    if(!fDoLightOutput){
      fHistClusterEtavsPhiBeforeAcc = new TH2F(Form("EtaPhi_beforeAcceptance %s",GetCutNumber().Data()),"EtaPhi_beforeAcceptance",nEmcalPhiBins,EmcalPhiBins,nEmcalEtaBins,EmcalEtaBins);
      fHistograms->Add(fHistClusterEtavsPhiBeforeAcc);
      fHistClusterEtavsPhiAfterAcc  = new TH2F(Form("EtaPhi_afterAcceptance %s",GetCutNumber().Data()),"EtaPhi_afterAcceptance",nEmcalPhiBins,EmcalPhiBins,nEmcalEtaBins,EmcalEtaBins);
      fHistograms->Add(fHistClusterEtavsPhiAfterAcc);
    }
    fHistClusterEtavsPhiAfterQA     = new TH2F(Form("EtaPhi_afterClusterQA %s",GetCutNumber().Data()),"EtaPhi_afterClusterQA",nEmcalPhiBins,EmcalPhiBins,nEmcalEtaBins,EmcalEtaBins);
    fHistograms->Add(fHistClusterEtavsPhiAfterQA);
  } else if( fClusterType == 2 ){ //PHOS
    const Int_t nPhosEtaBins        = 56;
    const Int_t nPhosPhiBins        = 192;
    const Float_t PhosEtaRange[2]   = {-0.16, 0.16};
    const Float_t PhosPhiRange[2]   = {4.5, 5.6};

    if(!fDoLightOutput){
      fHistClusterEtavsPhiBeforeAcc = new TH2F(Form("EtaPhi_beforeAcceptance %s",GetCutNumber().Data()),"EtaPhi_beforeAcceptance",nPhosPhiBins,PhosPhiRange[0],PhosPhiRange[1],nPhosEtaBins,PhosEtaRange[0],PhosEtaRange[1]);
      fHistograms->Add(fHistClusterEtavsPhiBeforeAcc);
      fHistClusterEtavsPhiAfterAcc  = new TH2F(Form("EtaPhi_afterAcceptance %s",GetCutNumber().Data()),"EtaPhi_afterAcceptance",nPhosPhiBins,PhosPhiRange[0],PhosPhiRange[1],nPhosEtaBins,PhosEtaRange[0],PhosEtaRange[1]);
      fHistograms->Add(fHistClusterEtavsPhiAfterAcc);
    }
    fHistClusterEtavsPhiAfterQA     = new TH2F(Form("EtaPhi_afterClusterQA %s",GetCutNumber().Data()),"EtaPhi_afterClusterQA",nPhosPhiBins,PhosPhiRange[0],PhosPhiRange[1],nPhosEtaBins,PhosEtaRange[0],PhosEtaRange[1]);
    fHistograms->Add(fHistClusterEtavsPhiAfterQA);
  } else if( fClusterType == 3){ //DCAL
    const Int_t nDcalEtaBins             = 96;
    const Int_t nDcalPhiBins             = 124;
    if(!fDoLightOutput){
      fHistClusterEtavsPhiBeforeAcc = new TH2F(Form("EtaPhi_beforeAcceptance %s",GetCutNumber().Data()),"EtaPhi_beforeAcceptance",nDcalPhiBins,4.5,5.7,nDcalEtaBins,-0.66687,0.66465);
      fHistograms->Add(fHistClusterEtavsPhiBeforeAcc);
      fHistClusterEtavsPhiAfterAcc  = new TH2F(Form("EtaPhi_afterAcceptance %s",GetCutNumber().Data()),"EtaPhi_afterAcceptance",nDcalPhiBins,4.5,5.7,nDcalEtaBins,-0.66687,0.66465);
      fHistograms->Add(fHistClusterEtavsPhiAfterAcc);
}
      fHistClusterEtavsPhiAfterQA     = new TH2F(Form("EtaPhi_afterClusterQA %s",GetCutNumber().Data()),"EtaPhi_afterClusterQA",nDcalPhiBins,4.5,5.7,nDcalEtaBins,-0.66687,0.66465);
      fHistograms->Add(fHistClusterEtavsPhiAfterQA);
  } else if( fClusterType == 0 ){ //all
    if(!fDoLightOutput){
      fHistClusterEtavsPhiBeforeAcc = new TH2F(Form("EtaPhi_beforeAcceptance %s",GetCutNumber().Data()),"EtaPhi_beforeAcceptance",462,0,2*TMath::Pi(),110,-0.7,0.7);
      fHistograms->Add(fHistClusterEtavsPhiBeforeAcc);
      fHistClusterEtavsPhiAfterAcc  = new TH2F(Form("EtaPhi_afterAcceptance %s",GetCutNumber().Data()),"EtaPhi_afterAcceptance",462,0,2*TMath::Pi(),110,-0.7,0.7);
      fHistograms->Add(fHistClusterEtavsPhiAfterAcc);
    }
    fHistClusterEtavsPhiAfterQA     = new TH2F(Form("EtaPhi_afterClusterQA %s",GetCutNumber().Data()),"EtaPhi_afterClusterQA",462,0,2*TMath::Pi(),110,-0.7,0.7);
    fHistograms->Add(fHistClusterEtavsPhiAfterQA);
  } else{AliError(Form("Cluster Type is not EMCAL nor PHOS nor all: %i",fClusterType));}

  if(fIsMC > 1){
    if(!fDoLightOutput){
      fHistClusterEtavsPhiBeforeAcc->Sumw2();
      fHistClusterEtavsPhiAfterAcc->Sumw2();
    }
    fHistClusterEtavsPhiAfterQA->Sumw2();
  }
  
  // Cluster quality related histograms
  Double_t timeMin                  = -2e-6;
  Double_t timeMax                  = 8e-6;
  if( fClusterType == 1 || fClusterType == 3){
    timeMin                         = -2e-7;
    timeMax                         = 12e-7;
  }
  
  if(!fDoLightOutput){
    fHistClusterTimevsEBeforeQA     = new TH2F(Form("ClusterTimeVsE_beforeClusterQA %s",GetCutNumber().Data()),"ClusterTimeVsE_beforeClusterQA",800,timeMin,timeMax,nBinsClusterECoarse,minClusterELog,maxClusterELog);
    SetLogBinningYTH2(fHistClusterTimevsEBeforeQA);
    fHistograms->Add(fHistClusterTimevsEBeforeQA);
  }
  fHistClusterTimevsEAfterQA        = new TH2F(Form("ClusterTimeVsE_afterClusterQA %s",GetCutNumber().Data()),"ClusterTimeVsE_afterClusterQA",800,timeMin,timeMax,nBinsClusterECoarse,minClusterELog,maxClusterELog);
  SetLogBinningYTH2(fHistClusterTimevsEAfterQA);
  fHistograms->Add(fHistClusterTimevsEAfterQA);
  if(fUseNonLinearity && !fDoLightOutput){
    fHistEnergyOfClusterBeforeNL    = new TH1F(Form("EnergyOfCluster_beforeNonLinearity %s",GetCutNumber().Data()),"EnergyOfCluster_beforeNonLinearity",nBinsClusterEFine,minClusterELog,maxClusterELog);
    SetLogBinningXTH1(fHistEnergyOfClusterBeforeNL);
    fHistograms->Add(fHistEnergyOfClusterBeforeNL);
    fHistEnergyOfClusterAfterNL     = new TH1F(Form("EnergyOfCluster_afterNonLinearity %s",GetCutNumber().Data()),"EnergyOfCluster_afterNonLinearity",nBinsClusterEFine,minClusterELog,maxClusterELog);
    SetLogBinningXTH1(fHistEnergyOfClusterAfterNL);
    fHistograms->Add(fHistEnergyOfClusterAfterNL);
  }
  fHistEnergyOfClusterBeforeQA      = new TH1F(Form("EnergyOfCluster_beforeClusterQA %s",GetCutNumber().Data()),"EnergyOfCluster_beforeClusterQA",nBinsClusterEFine,minClusterELog,maxClusterELog);
  SetLogBinningXTH1(fHistEnergyOfClusterBeforeQA);
  fHistograms->Add(fHistEnergyOfClusterBeforeQA);
  fHistEnergyOfClusterAfterQA       = new TH1F(Form("EnergyOfCluster_afterClusterQA %s",GetCutNumber().Data()),"EnergyOfCluster_afterClusterQA",nBinsClusterEFine,minClusterELog,maxClusterELog);
  SetLogBinningXTH1(fHistEnergyOfClusterAfterQA);
  fHistograms->Add(fHistEnergyOfClusterAfterQA);
  if(!fDoLightOutput){
    fHistNCellsBeforeQA             = new TH1F(Form("NCellPerCluster_beforeClusterQA %s",GetCutNumber().Data()),"NCellPerCluster_beforeClusterQA",50,0,50);
    fHistograms->Add(fHistNCellsBeforeQA);
    fHistNCellsAfterQA              = new TH1F(Form("NCellPerCluster_afterClusterQA %s",GetCutNumber().Data()),"NCellPerCluster_afterClusterQA",50,0,50);
    fHistograms->Add(fHistNCellsAfterQA);
    fHistM02BeforeQA                = new TH1F(Form("M02_beforeClusterQA %s",GetCutNumber().Data()),"M02_beforeClusterQA",400,0,5);
    fHistograms->Add(fHistM02BeforeQA);
    fHistM02AfterQA                 = new TH1F(Form("M02_afterClusterQA %s",GetCutNumber().Data()),"M02_afterClusterQA",400,0,5);
    fHistograms->Add(fHistM02AfterQA);
    fHistM20BeforeQA                = new TH1F(Form("M20_beforeClusterQA %s",GetCutNumber().Data()),"M20_beforeClusterQA",400,0,2.5);
    fHistograms->Add(fHistM20BeforeQA);
    fHistM20AfterQA                 = new TH1F(Form("M20_afterClusterQA %s",GetCutNumber().Data()),"M20_afterClusterQA",400,0,2.5);
    fHistograms->Add(fHistM20AfterQA);
    fHistDispersionBeforeQA         = new TH1F(Form("Dispersion_beforeClusterQA %s",GetCutNumber().Data()),"Dispersion_beforeClusterQA",100,0,4);
    fHistograms->Add(fHistDispersionBeforeQA);
    fHistDispersionAfterQA          = new TH1F(Form("Dispersion_afterClusterQA %s",GetCutNumber().Data()),"Dispersion_afterClusterQA",100,0,4);
    fHistograms->Add(fHistDispersionAfterQA);
    fHistNLMBeforeQA                = new TH1F(Form("NLM_beforeClusterQA %s",GetCutNumber().Data()),"NLM_beforeClusterQA",10,0,10);
    fHistograms->Add(fHistNLMBeforeQA);
    fHistNLMAfterQA                 = new TH1F(Form("NLM_afterClusterQA %s",GetCutNumber().Data()),"NLM_afterClusterQA",10,0,10);
    fHistograms->Add(fHistNLMAfterQA);
//     fHistNLMAvsNLMBBeforeQA         = new TH2F(Form("NLMAvsNLMB_beforeClusterQA %s",GetCutNumber().Data()),"NLMAvsNLMB_beforeClusterQA",10,0,10,10,0,10);
//     fHistograms->Add(fHistNLMAvsNLMBBeforeQA);
    fHistNLMVsNCellsAfterQA         = new TH2F(Form("NLM_NCells_afterClusterQA %s",GetCutNumber().Data()),"NLM_NCells_afterClusterQA",10,0,10,50,0,50);
    fHistograms->Add(fHistNLMVsNCellsAfterQA);
    fHistNLMVsEAfterQA              = new TH2F(Form("NLM_E_afterClusterQA %s",GetCutNumber().Data()),"NLM_E_afterClusterQA",10,0,10,nBinsClusterECoarse,minClusterELog,maxClusterELog);
    SetLogBinningYTH2(fHistNLMVsEAfterQA);
    fHistograms->Add(fHistNLMVsEAfterQA);
    if(fExtendedMatchAndQA > 1 || fIsPureCalo > 0){
      fHistClusterEM02AfterQA       = new TH2F(Form("EVsM02_afterClusterQA %s",GetCutNumber().Data()),"EVsM02_afterClusterQA",nBinsClusterEFine,minClusterELog,maxClusterELog,400,0,5);
      SetLogBinningXTH2(fHistClusterEM02AfterQA);
      fHistograms->Add(fHistClusterEM02AfterQA);
      fHistClusterEM02BeforeQA      = new TH2F(Form("EVsM02_beforeClusterQA %s",GetCutNumber().Data()),"EVsM02_beforeClusterQA",nBinsClusterEFine,minClusterELog,maxClusterELog,400,0,5);
      SetLogBinningXTH2(fHistClusterEM02BeforeQA);
      fHistograms->Add(fHistClusterEM02BeforeQA);
    }
  }
//----------------
  if(fIsMC > 1){
    if(!fDoLightOutput) fHistClusterTimevsEBeforeQA->Sumw2();
    fHistClusterTimevsEAfterQA->Sumw2();
    if(fUseNonLinearity && !fDoLightOutput){
      fHistEnergyOfClusterBeforeNL->Sumw2();
      fHistEnergyOfClusterAfterNL->Sumw2();
    }
    fHistEnergyOfClusterBeforeQA->Sumw2();
    fHistEnergyOfClusterAfterQA->Sumw2();
    if(!fDoLightOutput){
      fHistNCellsBeforeQA->Sumw2();
      fHistNCellsAfterQA->Sumw2();
      fHistM02BeforeQA->Sumw2();
      fHistM02AfterQA->Sumw2();
      fHistM20BeforeQA->Sumw2();
      fHistM20AfterQA->Sumw2();
      fHistDispersionBeforeQA->Sumw2();
      fHistDispersionAfterQA->Sumw2();
      fHistNLMBeforeQA->Sumw2();
      fHistNLMAfterQA->Sumw2();
//       fHistNLMAvsNLMBBeforeQA->Sumw2();
      fHistNLMVsNCellsAfterQA->Sumw2();
      fHistNLMVsEAfterQA->Sumw2();
      if(fExtendedMatchAndQA > 1 || fIsPureCalo > 0){
        fHistClusterEM02AfterQA->Sumw2();
        fHistClusterEM02BeforeQA->Sumw2();
      }
    }
  }
//----------------
  if(!fDoLightOutput){
    if( fClusterType == 1 ){
      Int_t nMaxCellsEMCAL      = fNMaxEMCalModules*48*24;
      fBadChannels              = new TProfile("EMCal - Bad Channels","EMCal - Bad Channels",nMaxCellsEMCAL,0,nMaxCellsEMCAL);
      fHistExtQA->Add(fBadChannels);
    } else if( fClusterType == 2 ){
      Int_t nMaxCellsPHOS       = fNMaxPHOSModules*56*64;
      fBadChannels              = new TProfile("PHOS - Bad Channels","PHOS - Bad Channels",nMaxCellsPHOS,0,nMaxCellsPHOS);
      fHistExtQA->Add(fBadChannels);
    } else if( fClusterType == 3 ){
      Int_t nStartCellDCAL      = 12288;
      Int_t nMaxCellsDCAL       = fNMaxDCalModules*32*24;
      fBadChannels              = new TProfile("DCAL - Bad Channels","DCAL - Bad Channels",nMaxCellsDCAL,nStartCellDCAL,nStartCellDCAL+nMaxCellsDCAL);
      fHistExtQA->Add(fBadChannels);
    }
  }
  
  if(fExtendedMatchAndQA > 1){
    if( fClusterType == 1 ){ //EMCAL
      // detailed cluster QA histos for EMCAL
      Int_t nMaxCellsEMCAL            = fNMaxEMCalModules*48*24;
      fHistClusterEnergyvsMod         = new TH2F(Form("ClusterEnergyVsModule_afterClusterQA %s",GetCutNumber().Data()),"ClusterEnergyVsModule_afterClusterQA", nBinsClusterEFine, minClusterELog, maxClusterELog, fNMaxEMCalModules, 0, fNMaxEMCalModules);
      SetLogBinningXTH2(fHistClusterEnergyvsMod);
      fHistExtQA->Add(fHistClusterEnergyvsMod);
      fHistNCellsBigger100MeVvsMod    = new TH2F(Form("NCellsAbove100VsModule %s",GetCutNumber().Data()),"NCellsAbove100VsModule",200,0,200,fNMaxEMCalModules,0,fNMaxEMCalModules);
      fHistExtQA->Add(fHistNCellsBigger100MeVvsMod);
      fHistNCellsBigger1500MeVvsMod   = new TH2F(Form("NCellsAbove1500VsModule %s",GetCutNumber().Data()),"NCellsAbove1500VsModule",100,0,100,fNMaxEMCalModules,0,fNMaxEMCalModules);
      fHistExtQA->Add(fHistNCellsBigger1500MeVvsMod);
      fHistEnergyOfModvsMod           = new TH2F(Form("ModuleEnergyVsModule %s",GetCutNumber().Data()),"ModuleEnergyVsModule",nBinsModuleECoarse, minModuleELog, maxModuleELog, 
                                                 fNMaxEMCalModules,0,fNMaxEMCalModules);
      SetLogBinningXTH2(fHistEnergyOfModvsMod);
      fHistExtQA->Add(fHistEnergyOfModvsMod);
      fHistClusterEnergyvsNCells      = new TH2F(Form("ClusterEnergyVsNCells_afterQA %s",GetCutNumber().Data()),"ClusterEnergyVsNCells_afterQA",nBinsClusterECoarse, minClusterELog, maxClusterELog, 50, 0, 50);
      SetLogBinningXTH2(fHistClusterEnergyvsNCells);
      fHistExtQA->Add(fHistClusterEnergyvsNCells);
      fHistClusterDistanceInTimeCut   = new TH2F(Form("ClusterDistanceTo_withinTimingCut %s",GetCutNumber().Data()),"ClusterDistanceTo_withinTimingCut; dist row; dist col",20,-10,10,20,-10,10);
      fHistClusterDistanceInTimeCut->Sumw2();
      fHistExtQA->Add(fHistClusterDistanceInTimeCut);
      fHistClusterDistanceOutTimeCut  = new TH2F(Form("ClusterDistanceTo_outsideTimingCut %s",GetCutNumber().Data()),"ClusterDistanceTo_outsideTimingCut; dist row; dist col",20,-10,10,20,-10,10);
      fHistClusterDistanceOutTimeCut->Sumw2();
      fHistExtQA->Add(fHistClusterDistanceOutTimeCut);
      fHistClusterDistance1DInTimeCut   = new TH1F(Form("Cluster1D_DistanceTo_withinTimingCut %s",GetCutNumber().Data()),"Cluster1D_DistanceTo_withinTimingCut; dist 1D; #entries",200,0.,0.5);
      fHistClusterDistance1DInTimeCut->Sumw2();
      fHistExtQA->Add(fHistClusterDistance1DInTimeCut);

      // detailed cell QA histos for EMCAL
      if(fExtendedMatchAndQA > 3){
        fHistCellEnergyvsCellID                 = new TH2F(Form("CellEnergyVsCellID %s",GetCutNumber().Data()),"CellEnergyVsCellID",nBinsCellECoarse, minCellELog, maxCellELog, nMaxCellsEMCAL, 0, nMaxCellsEMCAL);
        SetLogBinningXTH2(fHistCellEnergyvsCellID);
        fHistExtQA->Add(fHistCellEnergyvsCellID);
        fHistCellTimevsCellID                   = new TH2F(Form("CellTimeVsCellID %s",GetCutNumber().Data()),"CellTimeVsCellID",600,-timeMax,timeMax,nMaxCellsEMCAL,0,nMaxCellsEMCAL);
        fHistExtQA->Add(fHistCellTimevsCellID);
        fHistClusterIncludedCellsBeforeQA       = new TH1F(Form("ClusterIncludedCells_beforeClusterQA %s",GetCutNumber().Data()),"ClusterIncludedCells_beforeClusterQA",nMaxCellsEMCAL,0,nMaxCellsEMCAL);
        fHistExtQA->Add(fHistClusterIncludedCellsBeforeQA);
        fHistClusterIncludedCellsAfterQA        = new TH1F(Form("ClusterIncludedCells_afterClusterQA %s",GetCutNumber().Data()),"ClusterIncludedCells_afterClusterQA",nMaxCellsEMCAL,0,nMaxCellsEMCAL);
        fHistExtQA->Add(fHistClusterIncludedCellsAfterQA);
        fHistClusterEnergyFracCellsBeforeQA     = new TH1F(Form("ClusterEnergyFracCells_beforeClusterQA %s",GetCutNumber().Data()),"ClusterEnergyFracCells_beforeClusterQA",nMaxCellsEMCAL,0,nMaxCellsEMCAL);
        fHistClusterEnergyFracCellsBeforeQA->Sumw2();
        fHistExtQA->Add(fHistClusterEnergyFracCellsBeforeQA);
        fHistClusterEnergyFracCellsAfterQA      = new TH1F(Form("ClusterEnergyFracCells_afterClusterQA %s",GetCutNumber().Data()),"ClusterEnergyFracCells_afterClusterQA",nMaxCellsEMCAL,0,nMaxCellsEMCAL);
        fHistClusterEnergyFracCellsAfterQA->Sumw2();
        fHistExtQA->Add(fHistClusterEnergyFracCellsAfterQA);
        fHistClusterIncludedCellsTimingAfterQA  = new TH2F(Form("ClusterIncludedCellsTiming_afterClusterQA %s",GetCutNumber().Data()),"ClusterIncludedCellsTiming_afterClusterQA; cluster E (GeV);t (ns)", 
                                                           nBinsClusterECoarse, minClusterELog, maxClusterELog, 200, -500, 500);
        SetLogBinningXTH2(fHistClusterIncludedCellsTimingAfterQA);
        fHistClusterIncludedCellsTimingAfterQA->Sumw2();
        fHistExtQA->Add(fHistClusterIncludedCellsTimingAfterQA);
        fHistClusterIncludedCellsTimingEnergyAfterQA  = new TH2F(Form("ClusterIncludedCellsTimingEnergy_afterClusterQA %s",GetCutNumber().Data()),"ClusterIncludedCellsTimingEnergy_afterClusterQA; cell E (GeV);t (ns)",
                                                                 nBinsCellECoarse, minCellELog, maxCellELog, 200, -500, 500);
        SetLogBinningXTH2(fHistClusterIncludedCellsTimingEnergyAfterQA);
        fHistClusterIncludedCellsTimingEnergyAfterQA->Sumw2();
        fHistExtQA->Add(fHistClusterIncludedCellsTimingEnergyAfterQA);
      }
      
    } else if( fClusterType == 2 ){ //PHOS
      // detailed cluster QA histos for PHOS
      Int_t nMaxCellsPHOS             = fNMaxPHOSModules*56*64;
      fHistClusterEnergyvsMod         = new TH2F(Form("ClusterEnergyVsModule_afterClusterQA %s",GetCutNumber().Data()),"ClusterEnergyVsModule_afterClusterQA",nBinsClusterEFine, minClusterELog, maxClusterELog,
                                                 fNMaxPHOSModules, 0, fNMaxPHOSModules);
      SetLogBinningXTH2(fHistClusterEnergyvsMod);
      fHistExtQA->Add(fHistClusterEnergyvsMod);
      fHistNCellsBigger100MeVvsMod    = new TH2F(Form("NCellsAbove100VsModule %s",GetCutNumber().Data()),"NCellsAbove100VsModule",200,0,200,fNMaxPHOSModules,0,fNMaxPHOSModules);
      fHistExtQA->Add(fHistNCellsBigger100MeVvsMod);
      fHistNCellsBigger1500MeVvsMod   = new TH2F(Form("NCellsAbove1500VsModule %s",GetCutNumber().Data()),"NCellsAbove1500VsModule",100,0,100,fNMaxPHOSModules,0,fNMaxPHOSModules);
      fHistExtQA->Add(fHistNCellsBigger1500MeVvsMod);
      fHistEnergyOfModvsMod           = new TH2F(Form("ModuleEnergyVsModule %s",GetCutNumber().Data()),"ModuleEnergyVsModule",nBinsModuleECoarse, minModuleELog, maxModuleELog,fNMaxPHOSModules,0,fNMaxPHOSModules);
      SetLogBinningXTH2(fHistEnergyOfModvsMod);
      fHistExtQA->Add(fHistEnergyOfModvsMod);
      fHistClusterEnergyvsNCells      = new TH2F(Form("ClusterEnergyVsNCells_afterQA %s",GetCutNumber().Data()),"ClusterEnergyVsNCells_afterQA",nBinsClusterECoarse, minClusterELog, maxClusterELog, 50, 0, 50);
      SetLogBinningXTH2(fHistClusterEnergyvsNCells);
      fHistExtQA->Add(fHistClusterEnergyvsNCells);
      fHistClusterDistanceInTimeCut   = new TH2F(Form("ClusterDistanceTo_withinTimingCut %s",GetCutNumber().Data()),"ClusterDistanceTo_withinTimingCut; dist row; dist col",20,-10,10,20,-10,10);
      fHistClusterDistanceInTimeCut->Sumw2();
      fHistExtQA->Add(fHistClusterDistanceInTimeCut);
      fHistClusterDistanceOutTimeCut  = new TH2F(Form("ClusterDistanceTo_outsideTimingCut %s",GetCutNumber().Data()),"ClusterDistanceTo_outsideTimingCut; dist row; dist col",20,-10,10,20,-10,10);
      fHistClusterDistanceOutTimeCut->Sumw2();
      fHistExtQA->Add(fHistClusterDistanceOutTimeCut);
      fHistClusterDistance1DInTimeCut   = new TH1F(Form("Cluster1D_DistanceTo_withinTimingCut %s",GetCutNumber().Data()),"Cluster1D_DistanceTo_withinTimingCut; dist 1D; #entries",200,0.,0.5);
      fHistClusterDistance1DInTimeCut->Sumw2();
      fHistExtQA->Add(fHistClusterDistance1DInTimeCut);
      
      // detailed cell QA histos for PHOS
      if(fExtendedMatchAndQA > 3){
        fHistCellEnergyvsCellID         = new TH2F(Form("CellEnergyVsCellID %s",GetCutNumber().Data()),"CellEnergyVsCellID",nBinsCellECoarse, minCellELog, maxCellELog, nMaxCellsPHOS,0,nMaxCellsPHOS);
        SetLogBinningXTH2(fHistCellEnergyvsCellID);
        fHistExtQA->Add(fHistCellEnergyvsCellID);
        fHistCellTimevsCellID           = new TH2F(Form("CellTimeVsCellID %s",GetCutNumber().Data()),"CellTimeVsCellID",600,-timeMax,timeMax,nMaxCellsPHOS,0,nMaxCellsPHOS);
        fHistExtQA->Add(fHistCellTimevsCellID);
        fHistClusterIncludedCellsBeforeQA   = new TH1F(Form("ClusterIncludedCells_beforeClusterQA %s",GetCutNumber().Data()),"ClusterIncludedCells_beforeClusterQA",nMaxCellsPHOS,0,nMaxCellsPHOS);
        fHistExtQA->Add(fHistClusterIncludedCellsBeforeQA);
        fHistClusterIncludedCellsAfterQA    = new TH1F(Form("ClusterIncludedCells_afterClusterQA %s",GetCutNumber().Data()),"ClusterIncludedCells_afterClusterQA",nMaxCellsPHOS,0,nMaxCellsPHOS);
        fHistExtQA->Add(fHistClusterIncludedCellsAfterQA);
        fHistClusterEnergyFracCellsBeforeQA = new TH1F(Form("ClusterEnergyFracCells_beforeClusterQA %s",GetCutNumber().Data()),"ClusterEnergyFracCells_beforeClusterQA",nMaxCellsPHOS,0,nMaxCellsPHOS);
        fHistClusterEnergyFracCellsBeforeQA->Sumw2();
        fHistExtQA->Add(fHistClusterEnergyFracCellsBeforeQA);
        fHistClusterEnergyFracCellsAfterQA  = new TH1F(Form("ClusterEnergyFracCells_afterClusterQA %s",GetCutNumber().Data()),"ClusterEnergyFracCells_afterClusterQA",nMaxCellsPHOS,0,nMaxCellsPHOS);
        fHistClusterEnergyFracCellsAfterQA->Sumw2();
        fHistExtQA->Add(fHistClusterEnergyFracCellsAfterQA);
        fHistClusterIncludedCellsTimingAfterQA        = new TH2F(Form("ClusterIncludedCellsTiming_afterClusterQA %s",GetCutNumber().Data()),"ClusterIncludedCellsTiming_afterClusterQA; E (GeV);t (ns)", 
                                                                 nBinsClusterECoarse, minClusterELog, maxClusterELog, 200, -500, 500);
        SetLogBinningXTH2(fHistClusterIncludedCellsTimingAfterQA);
        fHistClusterIncludedCellsTimingAfterQA->Sumw2();
        fHistExtQA->Add(fHistClusterIncludedCellsTimingAfterQA);
        fHistClusterIncludedCellsTimingEnergyAfterQA  = new TH2F(Form("ClusterIncludedCellsTimingEnergy_afterClusterQA %s",GetCutNumber().Data()),"ClusterIncludedCellsTimingEnergy_afterClusterQA; cell E (GeV);t (ns)",
                                                                 nBinsCellECoarse, minCellELog, maxCellELog, 200, -500, 500);
        SetLogBinningXTH2(fHistClusterIncludedCellsTimingEnergyAfterQA);
        fHistClusterIncludedCellsTimingEnergyAfterQA->Sumw2();
        fHistExtQA->Add(fHistClusterIncludedCellsTimingEnergyAfterQA);
      }
    } else if( fClusterType == 3 ){ //DCAL
      Int_t nModulesStart = 12;
      Int_t nCellsStart = 12288;
      Int_t nMaxCellsDCAL            = fNMaxDCalModules*32*24;

      fHistClusterEnergyvsMod         = new TH2F(Form("ClusterEnergyVsModule_afterClusterQA %s",GetCutNumber().Data()),"ClusterEnergyVsModule_afterClusterQA",nBinsClusterEFine,minClusterELog, maxClusterELog, fNMaxDCalModules,nModulesStart,fNMaxDCalModules+nModulesStart);
      SetLogBinningXTH2(fHistClusterEnergyvsMod);
      fHistExtQA->Add(fHistClusterEnergyvsMod);
      fHistNCellsBigger100MeVvsMod    = new TH2F(Form("NCellsAbove100VsModule %s",GetCutNumber().Data()),"NCellsAbove100VsModule",200,0,200,fNMaxDCalModules,nModulesStart,fNMaxDCalModules+nModulesStart);
      fHistExtQA->Add(fHistNCellsBigger100MeVvsMod);
      fHistNCellsBigger1500MeVvsMod   = new TH2F(Form("NCellsAbove1500VsModule %s",GetCutNumber().Data()),"NCellsAbove1500VsModule",100,0,100,fNMaxDCalModules,nModulesStart,nModulesStart+fNMaxDCalModules);
      fHistExtQA->Add(fHistNCellsBigger1500MeVvsMod);
      fHistEnergyOfModvsMod           = new TH2F(Form("ModuleEnergyVsModule %s",GetCutNumber().Data()),"ModuleEnergyVsModule",1000,0,100,fNMaxDCalModules,nModulesStart,fNMaxDCalModules+nModulesStart);
      //SetLogBinningXTH2(fHistEnergyOfModvsMod);
      fHistExtQA->Add(fHistEnergyOfModvsMod);
      fHistClusterEnergyvsNCells      = new TH2F(Form("ClusterEnergyVsNCells_afterQA %s",GetCutNumber().Data()),"ClusterEnergyVsNCells_afterQA",nBinsClusterECoarse, minClusterELog, maxClusterELog,50,0,50);
      SetLogBinningXTH2(fHistClusterEnergyvsNCells);
      fHistExtQA->Add(fHistClusterEnergyvsNCells);
      fHistClusterDistanceInTimeCut   = new TH2F(Form("ClusterDistanceTo_withinTimingCut %s",GetCutNumber().Data()),"ClusterDistanceTo_withinTimingCut; dist row; dist col",20,-10,10,20,-10,10);
      fHistClusterDistanceInTimeCut->Sumw2();
      fHistExtQA->Add(fHistClusterDistanceInTimeCut);
      fHistClusterDistanceOutTimeCut  = new TH2F(Form("ClusterDistanceTo_outsideTimingCut %s",GetCutNumber().Data()),"ClusterDistanceTo_outsideTimingCut; dist row; dist col",20,-10,10,20,-10,10);
      fHistClusterDistanceOutTimeCut->Sumw2();
      fHistExtQA->Add(fHistClusterDistanceOutTimeCut);
      fHistClusterDistance1DInTimeCut   = new TH1F(Form("Cluster1D_DistanceTo_withinTimingCut %s",GetCutNumber().Data()),"Cluster1D_DistanceTo_withinTimingCut; dist 1D; #entries",200,0.,0.5);
      fHistClusterDistance1DInTimeCut->Sumw2();
      fHistExtQA->Add(fHistClusterDistance1DInTimeCut);

      // detailed cell QA histos for DCAL
      if(fExtendedMatchAndQA > 3){
        fHistCellEnergyvsCellID         = new TH2F(Form("CellEnergyVsCellID %s",GetCutNumber().Data()),"CellEnergyVsCellID",nBinsCellECoarse, minCellELog, maxCellELog,nMaxCellsDCAL,nCellsStart,nMaxCellsDCAL+nCellsStart);
	SetLogBinningXTH2(fHistCellEnergyvsCellID);
        fHistExtQA->Add(fHistCellEnergyvsCellID);
        fHistCellTimevsCellID           = new TH2F(Form("CellTimeVsCellID %s",GetCutNumber().Data()),"CellTimeVsCellID",600,-timeMax,timeMax,nMaxCellsDCAL,nCellsStart,nMaxCellsDCAL+nCellsStart);
        fHistExtQA->Add(fHistCellTimevsCellID);
        fHistClusterIncludedCellsBeforeQA       = new TH1F(Form("ClusterIncludedCells_beforeClusterQA %s",GetCutNumber().Data()),"ClusterIncludedCells_beforeClusterQA",nMaxCellsDCAL,nCellsStart,nMaxCellsDCAL+nCellsStart);
        fHistExtQA->Add(fHistClusterIncludedCellsBeforeQA);
        fHistClusterIncludedCellsAfterQA        = new TH1F(Form("ClusterIncludedCells_afterClusterQA %s",GetCutNumber().Data()),"ClusterIncludedCells_afterClusterQA",nMaxCellsDCAL,nCellsStart,nMaxCellsDCAL+nCellsStart);
        fHistExtQA->Add(fHistClusterIncludedCellsAfterQA);
        fHistClusterEnergyFracCellsBeforeQA     = new TH1F(Form("ClusterEnergyFracCells_beforeClusterQA %s",GetCutNumber().Data()),"ClusterEnergyFracCells_beforeClusterQA",nMaxCellsDCAL,nCellsStart,nMaxCellsDCAL+nCellsStart);
        fHistClusterEnergyFracCellsBeforeQA->Sumw2();
        fHistExtQA->Add(fHistClusterEnergyFracCellsBeforeQA);
        fHistClusterEnergyFracCellsAfterQA      = new TH1F(Form("ClusterEnergyFracCells_afterClusterQA %s",GetCutNumber().Data()),"ClusterEnergyFracCells_afterClusterQA",nMaxCellsDCAL,nCellsStart,nMaxCellsDCAL+nCellsStart);
        fHistClusterEnergyFracCellsAfterQA->Sumw2();
        fHistExtQA->Add(fHistClusterEnergyFracCellsAfterQA);
        fHistClusterIncludedCellsTimingAfterQA  = new TH2F(Form("ClusterIncludedCellsTiming_afterClusterQA %s",GetCutNumber().Data()),"ClusterIncludedCellsTiming_afterClusterQA; cluster E (GeV);t (ns)",nBinsClusterECoarse, minClusterELog, maxClusterELog,200,-500,500);
        SetLogBinningXTH2(fHistClusterIncludedCellsTimingAfterQA);
        fHistClusterIncludedCellsTimingAfterQA->Sumw2();
        fHistExtQA->Add(fHistClusterIncludedCellsTimingAfterQA);
        fHistClusterIncludedCellsTimingEnergyAfterQA  = new TH2F(Form("ClusterIncludedCellsTimingEnergy_afterClusterQA %s",GetCutNumber().Data()),"ClusterIncludedCellsTimingEnergy_afterClusterQA; cell E (GeV);t (ns)",nBinsCellECoarse, minCellELog, maxCellELog,200,-500,500);
        SetLogBinningXTH2(fHistClusterIncludedCellsTimingEnergyAfterQA);
        fHistClusterIncludedCellsTimingEnergyAfterQA->Sumw2();
        fHistExtQA->Add(fHistClusterIncludedCellsTimingEnergyAfterQA);
      }
    } else{AliError(Form("fExtendedMatchAndQA (%i) not (yet) defined for cluster type (%i)",fExtendedMatchAndQA,fClusterType));}
  }

  //TrackMatching histograms
  if(fUseDistTrackToCluster && !fDoLightOutput) {
    const Int_t nEtaBins        = 300;
    const Int_t nPhiBins        = 300;
    const Float_t EtaRange[2]   = {-0.3, 0.3};
    const Float_t PhiRange[2]   = {-0.3, 0.3};

    fHistClusterdEtadPhiBeforeQA    = new TH2F(Form("dEtaVsdPhi_beforeClusterQA %s",GetCutNumber().Data()),"dEtaVsdPhi_beforeClusterQA", nEtaBins, EtaRange[0], EtaRange[1], nPhiBins, PhiRange[0], PhiRange[1]);
    fHistograms->Add(fHistClusterdEtadPhiBeforeQA);
    fHistClusterdEtadPhiAfterQA     = new TH2F(Form("dEtaVsdPhi_afterClusterQA %s",GetCutNumber().Data()),"dEtaVsdPhi_afterClusterQA", nEtaBins, EtaRange[0], EtaRange[1], nPhiBins, PhiRange[0], PhiRange[1]);
    fHistograms->Add(fHistClusterdEtadPhiAfterQA);
    //----------------
    if(fIsMC > 1){
      fHistClusterdEtadPhiBeforeQA->Sumw2();
      fHistClusterdEtadPhiAfterQA->Sumw2();
    }
    //----------------
    if(fExtendedMatchAndQA == 1 || fExtendedMatchAndQA == 3 || fExtendedMatchAndQA == 5 ){
        
      // QA histograms for track matching
      fHistClusterRBeforeQA                           = new TH1F(Form("R_Cluster_beforeClusterQA %s",GetCutNumber().Data()),"R of cluster",200,400,500);
      fHistograms->Add(fHistClusterRBeforeQA);
      fHistClusterRAfterQA                            = new TH1F(Form("R_Cluster_afterClusterQA %s",GetCutNumber().Data()),"R of cluster_matched",200,400,500);
      fHistograms->Add(fHistClusterRAfterQA);
      fHistDistanceTrackToClusterBeforeQA             = new TH1F(Form("DistanceToTrack_beforeClusterQA %s",GetCutNumber().Data()),"DistanceToTrack_beforeClusterQA",200,0,2);
      fHistograms->Add(fHistDistanceTrackToClusterBeforeQA);
      fHistDistanceTrackToClusterAfterQA              = new TH1F(Form("DistanceToTrack_afterClusterQA %s",GetCutNumber().Data()),"DistanceToTrack_afterClusterQA",200,0,2);
      fHistograms->Add(fHistDistanceTrackToClusterAfterQA);
      fHistClusterdEtadPhiPosTracksBeforeQA           = new TH2F(Form("dEtaVsdPhi_posTracks_beforeClusterQA %s",GetCutNumber().Data()),"dEtaVsdPhi_posTracks_beforeClusterQA",
                                                                 nEtaBins, EtaRange[0], EtaRange[1], nPhiBins, PhiRange[0], PhiRange[1]);
      fHistograms->Add(fHistClusterdEtadPhiPosTracksBeforeQA);
      fHistClusterdEtadPhiNegTracksBeforeQA           = new TH2F(Form("dEtaVsdPhi_negTracks_beforeClusterQA %s",GetCutNumber().Data()),"dEtaVsdPhi_negTracks_beforeClusterQA", 
                                                                 nEtaBins, EtaRange[0], EtaRange[1], nPhiBins, PhiRange[0], PhiRange[1]);
      fHistograms->Add(fHistClusterdEtadPhiNegTracksBeforeQA);
      fHistClusterdEtadPhiPosTracksAfterQA            = new TH2F(Form("dEtaVsdPhi_posTracks_afterClusterQA %s",GetCutNumber().Data()),"dEtaVsdPhi_posTracks_afterClusterQA", 
                                                                 nEtaBins, EtaRange[0], EtaRange[1], nPhiBins, PhiRange[0], PhiRange[1]);
      fHistograms->Add(fHistClusterdEtadPhiPosTracksAfterQA);
      fHistClusterdEtadPhiNegTracksAfterQA            = new TH2F(Form("dEtaVsdPhi_negTracks_afterClusterQA %s",GetCutNumber().Data()),"dEtaVsdPhi_negTracks_afterClusterQA", 
                                                                 nEtaBins, EtaRange[0], EtaRange[1], nPhiBins, PhiRange[0], PhiRange[1]);
      fHistograms->Add(fHistClusterdEtadPhiNegTracksAfterQA);
      fHistClusterdEtadPhiPosTracksP_000_075BeforeQA  = new TH2F(Form("dEtaVsdPhi_posTracks_P<0.75_beforeClusterQA %s",GetCutNumber().Data()),"dEtaVsdPhi_posTracks_P<0.75_beforeClusterQA", 
                                                                 nEtaBins, EtaRange[0], EtaRange[1], nPhiBins, PhiRange[0], PhiRange[1]);
      fHistograms->Add(fHistClusterdEtadPhiPosTracksP_000_075BeforeQA);
      fHistClusterdEtadPhiPosTracksP_075_125BeforeQA  = new TH2F(Form("dEtaVsdPhi_posTracks_0.75<P<1.25_beforeClusterQA %s",GetCutNumber().Data()),"dEtaVsdPhi_posTracks_0.75<P<1.25_beforeClusterQA", 
                                                                 nEtaBins, EtaRange[0], EtaRange[1], nPhiBins, PhiRange[0], PhiRange[1]);
      fHistograms->Add(fHistClusterdEtadPhiPosTracksP_075_125BeforeQA);
      fHistClusterdEtadPhiPosTracksP_125_999BeforeQA  = new TH2F(Form("dEtaVsdPhi_posTracks_P>1.25_beforeClusterQA %s",GetCutNumber().Data()),"dEtaVsdPhi_posTracks_P>1.25_beforeClusterQA",
                                                                 nEtaBins, EtaRange[0], EtaRange[1], nPhiBins, PhiRange[0], PhiRange[1]);
      fHistograms->Add(fHistClusterdEtadPhiPosTracksP_125_999BeforeQA);
      fHistClusterdEtadPhiNegTracksP_000_075BeforeQA  = new TH2F(Form("dEtaVsdPhi_negTracks_P<0.75_beforeClusterQA %s",GetCutNumber().Data()),"dEtaVsdPhi_negTracks_P<0.75_beforeClusterQA", 
                                                                 nEtaBins, EtaRange[0], EtaRange[1], nPhiBins, PhiRange[0], PhiRange[1]);
      fHistograms->Add(fHistClusterdEtadPhiNegTracksP_000_075BeforeQA);
      fHistClusterdEtadPhiNegTracksP_075_125BeforeQA  = new TH2F(Form("dEtaVsdPhi_negTracks_0.75<P<1.25_beforeClusterQA %s",GetCutNumber().Data()),"dEtaVsdPhi_negTracks_0.75<P<1.25_beforeClusterQA", 
                                                                 nEtaBins, EtaRange[0], EtaRange[1], nPhiBins, PhiRange[0], PhiRange[1]);
      fHistograms->Add(fHistClusterdEtadPhiNegTracksP_075_125BeforeQA);
      fHistClusterdEtadPhiNegTracksP_125_999BeforeQA  = new TH2F(Form("dEtaVsdPhi_negTracks_P>1.25_beforeClusterQA %s",GetCutNumber().Data()),"dEtaVsdPhi_negTracks_P>1.25_beforeClusterQA", 
                                                                 nEtaBins, EtaRange[0], EtaRange[1], nPhiBins, PhiRange[0], PhiRange[1]);
      fHistograms->Add(fHistClusterdEtadPhiNegTracksP_125_999BeforeQA);
      fHistClusterdEtadPtBeforeQA                     = new TH2F(Form("dEtaVsPt_beforeClusterQA %s",GetCutNumber().Data()),"dEtaVsPt_beforeClusterQA",nEtaBins,EtaRange[0],EtaRange[1], 
                                                                 nBinsModuleECoarse,  minClusterELog, maxClusterELog);
      SetLogBinningYTH2(fHistClusterdEtadPtBeforeQA);
      fHistograms->Add(fHistClusterdEtadPtBeforeQA);
      fHistClusterdEtadPtAfterQA                      = new TH2F(Form("dEtaVsPt_afterClusterQA %s",GetCutNumber().Data()),"dEtaVsPt_afterClusterQA",nEtaBins,EtaRange[0],EtaRange[1], 
                                                                 nBinsModuleECoarse,  minClusterELog, maxClusterELog);
      SetLogBinningYTH2(fHistClusterdEtadPtAfterQA);
      fHistograms->Add(fHistClusterdEtadPtAfterQA);
      fHistClusterdPhidPtPosTracksBeforeQA            = new TH2F(Form("dPhiVsPt_posTracks_beforeClusterQA %s",GetCutNumber().Data()),"dPhiVsPt_posTracks_beforeClusterQA", 
                                                                 2*nPhiBins, 2*PhiRange[0], 2*PhiRange[1], nBinsModuleECoarse,  minClusterELog, maxClusterELog);
      SetLogBinningYTH2(fHistClusterdPhidPtPosTracksBeforeQA);
      fHistograms->Add(fHistClusterdPhidPtPosTracksBeforeQA);
      fHistClusterdPhidPtNegTracksBeforeQA            = new TH2F(Form("dPhiVsPt_negTracks_beforeClusterQA %s",GetCutNumber().Data()),"dPhiVsPt_negTracks_beforeClusterQA", 
                                                                 2*nPhiBins, 2*PhiRange[0], 2*PhiRange[1], nBinsModuleECoarse,  minClusterELog, maxClusterELog);
      SetLogBinningYTH2(fHistClusterdPhidPtNegTracksBeforeQA);
      fHistograms->Add(fHistClusterdPhidPtNegTracksBeforeQA);
      fHistClusterdPhidPtAfterQA                      = new TH2F(Form("dPhiVsPt_afterClusterQA %s",GetCutNumber().Data()),"dPhiVsPt_afterClusterQA",2*nPhiBins,2*PhiRange[0],2*PhiRange[1], 
                                                                 nBinsModuleECoarse,  minClusterELog, maxClusterELog);
      SetLogBinningYTH2(fHistClusterdPhidPtAfterQA);
      fHistograms->Add(fHistClusterdPhidPtAfterQA);
      
      if(fIsMC > 0){
        fHistClusterdEtadPtTrueMatched                = new TH2F(Form("dEtaVsPt_TrueMatched %s",GetCutNumber().Data()),"dEtaVsPt_TrueMatched",nEtaBins,EtaRange[0],EtaRange[1],
                                                                  nBinsModuleECoarse,  minClusterELog, maxClusterELog);
        SetLogBinningYTH2(fHistClusterdEtadPtTrueMatched);
        fHistograms->Add(fHistClusterdEtadPtTrueMatched);
        fHistClusterdPhidPtPosTracksTrueMatched       = new TH2F(Form("dPhiVsPt_posTracks_TrueMatched %s",GetCutNumber().Data()),"dPhiVsPt_posTracks_TrueMatched",2*nPhiBins,2*PhiRange[0],2*PhiRange[1], 
                                                                  nBinsModuleECoarse,  minClusterELog, maxClusterELog);
        SetLogBinningYTH2(fHistClusterdPhidPtPosTracksTrueMatched);
        fHistograms->Add(fHistClusterdPhidPtPosTracksTrueMatched);
        fHistClusterdPhidPtNegTracksTrueMatched       = new TH2F(Form("dPhiVsPt_negTracks_TrueMatched %s",GetCutNumber().Data()),"dPhiVsPt_negTracks_TrueMatched",2*nPhiBins,2*PhiRange[0],2*PhiRange[1],
                                                                  nBinsModuleECoarse,  minClusterELog, maxClusterELog);
        SetLogBinningYTH2(fHistClusterdPhidPtNegTracksTrueMatched);
        fHistograms->Add(fHistClusterdPhidPtNegTracksTrueMatched);
      }

      // QA histos for shower shape correlations
      fHistClusterM20M02BeforeQA                      = new TH2F(Form("M20VsM02_beforeClusterQA %s",GetCutNumber().Data()),"M20VsM02_beforeClusterQA",200,0,2.5,400,0,5);
      fHistograms->Add(fHistClusterM20M02BeforeQA);
      fHistClusterM20M02AfterQA                       = new TH2F(Form("M20VsM02_afterClusterQA %s",GetCutNumber().Data()),"M20VsM02_afterClusterQA",200,0,2.5,400,0,5);
      fHistograms->Add(fHistClusterM20M02AfterQA);

      //----------------
      if(fIsMC > 1){
        fHistClusterRBeforeQA->Sumw2();
        fHistClusterRAfterQA->Sumw2();
        fHistDistanceTrackToClusterBeforeQA->Sumw2();
        fHistDistanceTrackToClusterAfterQA->Sumw2();
        fHistClusterdEtadPhiPosTracksBeforeQA->Sumw2();
        fHistClusterdEtadPhiNegTracksBeforeQA->Sumw2();
        fHistClusterdEtadPhiPosTracksAfterQA->Sumw2();
        fHistClusterdEtadPhiNegTracksAfterQA->Sumw2();
        fHistClusterdEtadPhiPosTracksP_000_075BeforeQA->Sumw2();
        fHistClusterdEtadPhiPosTracksP_075_125BeforeQA->Sumw2();
        fHistClusterdEtadPhiPosTracksP_125_999BeforeQA->Sumw2();
        fHistClusterdEtadPhiNegTracksP_000_075BeforeQA->Sumw2();
        fHistClusterdEtadPhiNegTracksP_075_125BeforeQA->Sumw2();
        fHistClusterdEtadPhiNegTracksP_125_999BeforeQA->Sumw2();
        fHistClusterdEtadPtBeforeQA->Sumw2();
        fHistClusterdEtadPtAfterQA->Sumw2();
        fHistClusterdPhidPtPosTracksBeforeQA->Sumw2();
        fHistClusterdPhidPtNegTracksBeforeQA->Sumw2();
        fHistClusterdPhidPtAfterQA->Sumw2();
        fHistClusterdEtadPtTrueMatched->Sumw2();
        fHistClusterdPhidPtPosTracksTrueMatched->Sumw2();
        fHistClusterdPhidPtNegTracksTrueMatched->Sumw2();
        fHistClusterM20M02BeforeQA->Sumw2();
        fHistClusterM20M02AfterQA->Sumw2();
      }
      //----------------
    }
  }
  if( fUseDistTrackToCluster && fIsMC && (fExtendedMatchAndQA == 1 || fExtendedMatchAndQA == 3 || fExtendedMatchAndQA == 5 )){
    // TM efficiency histograms
    fHistClusterTMEffiInput                       = new TH2F(Form("TMEffiInputHisto %s",GetCutNumber().Data()),"TMEffiInputHisto",nBinsClusterEFine, minClusterELog, maxClusterELog, 22, -0.5, 21.5);
    SetLogBinningXTH2(fHistClusterTMEffiInput);
    fHistClusterTMEffiInput->GetYaxis()->SetBinLabel(1,"All cl");
    fHistClusterTMEffiInput->GetYaxis()->SetBinLabel(2,"Ch cl");
    fHistClusterTMEffiInput->GetYaxis()->SetBinLabel(3,"Ne cl");
    fHistClusterTMEffiInput->GetYaxis()->SetBinLabel(4,"Ne cl sub ch");
    fHistClusterTMEffiInput->GetYaxis()->SetBinLabel(5,"Ga cl");
    fHistClusterTMEffiInput->GetYaxis()->SetBinLabel(6,"Ga cl sub ch");
    fHistClusterTMEffiInput->GetYaxis()->SetBinLabel(7,"conv cl");
    fHistClusterTMEffiInput->GetYaxis()->SetBinLabel(8,"Ch cl prim");
    fHistClusterTMEffiInput->GetYaxis()->SetBinLabel(9,"El cl");
    fHistClusterTMEffiInput->GetYaxis()->SetBinLabel(10,"All cl match");
    fHistClusterTMEffiInput->GetYaxis()->SetBinLabel(11,"Ch cl match");
    fHistClusterTMEffiInput->GetYaxis()->SetBinLabel(12,"Ch cl match w lead");
    fHistClusterTMEffiInput->GetYaxis()->SetBinLabel(13,"Ne cl match");
    fHistClusterTMEffiInput->GetYaxis()->SetBinLabel(14,"Ne cl sub ch match");
    fHistClusterTMEffiInput->GetYaxis()->SetBinLabel(15,"Ga cl match");
    fHistClusterTMEffiInput->GetYaxis()->SetBinLabel(16,"Ga cl sub ch match");
    fHistClusterTMEffiInput->GetYaxis()->SetBinLabel(17,"conv cl match");
    fHistClusterTMEffiInput->GetYaxis()->SetBinLabel(18,"conv cl match w lead");
    fHistClusterTMEffiInput->GetYaxis()->SetBinLabel(19,"Ch cl match");
    fHistClusterTMEffiInput->GetYaxis()->SetBinLabel(20,"Ch cl match w lead");
    fHistClusterTMEffiInput->GetYaxis()->SetBinLabel(21,"El cl match");
    fHistClusterTMEffiInput->GetYaxis()->SetBinLabel(22,"El cl match w lead");
    fHistClusterTMEffiInput->GetXaxis()->SetTitle("#it{E}_{cl} (GeV)");
    fHistograms->Add(fHistClusterTMEffiInput);

    fHistClusterEvsTrackECharged                  = new TH2F(Form("ClusterE_TrackE_ChargedCluster %s",GetCutNumber().Data()),"ClusterE TrackE ChargedCluster",
                                                             nBinsClusterEFine, minClusterELog, maxClusterELog, nBinsClusterEFine, minClusterELog, maxClusterELog);
    SetLogBinningXTH2(fHistClusterEvsTrackECharged);
    SetLogBinningYTH2(fHistClusterEvsTrackECharged);
    fHistClusterEvsTrackECharged->GetXaxis()->SetTitle("#it{E}_{cl} (GeV)");
    fHistClusterEvsTrackECharged->GetYaxis()->SetTitle("#it{E}_{tr} (GeV)");
    fHistograms->Add(fHistClusterEvsTrackECharged);
    fHistClusterEvsTrackEChargedLead              = new TH2F(Form("ClusterE_TrackE_ChargedCluster_LeadMatched %s",GetCutNumber().Data()),"ClusterE TrackE ChargedCluster LeadMatched",
                                                             nBinsClusterEFine, minClusterELog, maxClusterELog, nBinsClusterEFine, minClusterELog, maxClusterELog);
    SetLogBinningXTH2(fHistClusterEvsTrackEChargedLead);
    SetLogBinningYTH2(fHistClusterEvsTrackEChargedLead);
    fHistClusterEvsTrackEChargedLead->GetXaxis()->SetTitle("#it{E}_{cl} (GeV)");
    fHistClusterEvsTrackEChargedLead->GetYaxis()->SetTitle("#it{E}_{tr} (GeV)");
    fHistograms->Add(fHistClusterEvsTrackEChargedLead);
    fHistClusterEvsTrackENeutral                  = new TH2F(Form("ClusterE_TrackE_NeutralCluster %s",GetCutNumber().Data()),"ClusterE TrackE NeutralCluster",
                                                             nBinsClusterEFine, minClusterELog, maxClusterELog, nBinsClusterEFine, minClusterELog, maxClusterELog);
    SetLogBinningXTH2(fHistClusterEvsTrackENeutral);
    SetLogBinningYTH2(fHistClusterEvsTrackENeutral);
    fHistClusterEvsTrackENeutral->GetXaxis()->SetTitle("#it{E}_{cl} (GeV)");
    fHistClusterEvsTrackENeutral->GetYaxis()->SetTitle("#it{E}_{tr} (GeV)");
    fHistograms->Add(fHistClusterEvsTrackENeutral);
    fHistClusterEvsTrackENeutralSubCharged        = new TH2F(Form("ClusterE_TrackE_NeutralClusterSubCharged %s",GetCutNumber().Data()),"ClusterE TrackE NeutralCluster SubCharged",
                                                             nBinsClusterEFine, minClusterELog, maxClusterELog, nBinsClusterEFine, minClusterELog, maxClusterELog);
    SetLogBinningXTH2(fHistClusterEvsTrackENeutralSubCharged);
    SetLogBinningYTH2(fHistClusterEvsTrackENeutralSubCharged);
    fHistClusterEvsTrackENeutralSubCharged->GetXaxis()->SetTitle("#it{E}_{cl} (GeV)");
    fHistClusterEvsTrackENeutralSubCharged->GetYaxis()->SetTitle("#it{E}_{tr} (GeV)");
    fHistograms->Add(fHistClusterEvsTrackENeutralSubCharged);
    fHistClusterEvsTrackEGamma                    = new TH2F(Form("ClusterE_TrackE_GammaCluster %s",GetCutNumber().Data()),"ClusterE TrackE GammaCluster",
                                                             nBinsClusterEFine, minClusterELog, maxClusterELog, nBinsClusterEFine, minClusterELog, maxClusterELog);
    SetLogBinningXTH2(fHistClusterEvsTrackEGamma);
    SetLogBinningYTH2(fHistClusterEvsTrackEGamma);
    fHistClusterEvsTrackEGamma->GetXaxis()->SetTitle("#it{E}_{cl} (GeV)");
    fHistClusterEvsTrackEGamma->GetYaxis()->SetTitle("#it{E}_{tr} (GeV)");
    fHistograms->Add(fHistClusterEvsTrackEGamma);
    fHistClusterEvsTrackEGammaSubCharged          = new TH2F(Form("ClusterE_TrackE_GammaClusterSubCharged %s",GetCutNumber().Data()),"ClusterE TrackE GammaCluster SubCharged",
                                                             nBinsClusterEFine, minClusterELog, maxClusterELog, nBinsClusterEFine, minClusterELog, maxClusterELog);
    SetLogBinningXTH2(fHistClusterEvsTrackEGammaSubCharged);
    SetLogBinningYTH2(fHistClusterEvsTrackEGammaSubCharged);
    fHistClusterEvsTrackEGammaSubCharged->GetXaxis()->SetTitle("#it{E}_{cl} (GeV)");
    fHistClusterEvsTrackEGammaSubCharged->GetYaxis()->SetTitle("#it{E}_{tr} (GeV)");
    fHistograms->Add(fHistClusterEvsTrackEGammaSubCharged);
    fHistClusterEvsTrackEConv                    = new TH2F(Form("ClusterE_TrackE_ConvCluster %s",GetCutNumber().Data()),"ClusterE TrackE ConvCluster",
                                                             nBinsClusterEFine, minClusterELog, maxClusterELog, nBinsClusterEFine, minClusterELog, maxClusterELog);
    SetLogBinningXTH2(fHistClusterEvsTrackEConv);
    SetLogBinningYTH2(fHistClusterEvsTrackEConv);
    fHistClusterEvsTrackEConv->GetXaxis()->SetTitle("#it{E}_{cl} (GeV)");
    fHistClusterEvsTrackEConv->GetYaxis()->SetTitle("#it{E}_{tr} (GeV)");
    fHistograms->Add(fHistClusterEvsTrackEConv);

    fHistClusterENMatchesNeutral                  = new TH2F(Form("ClusterE_NMatches_NeutralCluster %s",GetCutNumber().Data()),"ClusterE NMatches NeutralCluster",
                                                             nBinsClusterEFine, minClusterELog, maxClusterELog, 20, -0.5, 19.5);
    SetLogBinningXTH2(fHistClusterENMatchesNeutral);
    fHistClusterENMatchesNeutral->GetXaxis()->SetTitle("#it{E}_{cl} (GeV)");
    fHistClusterENMatchesNeutral->GetYaxis()->SetTitle("#it{N}_{matches}");
    fHistograms->Add(fHistClusterENMatchesNeutral);
    fHistClusterENMatchesCharged                  = new TH2F(Form("ClusterE_NMatches_ChargedCluster %s",GetCutNumber().Data()),"ClusterE NMatches ChargedCluster",
                                                             nBinsClusterEFine, minClusterELog, maxClusterELog, 20, -0.5, 19.5);
    SetLogBinningXTH2(fHistClusterENMatchesCharged);
    fHistClusterENMatchesCharged->GetXaxis()->SetTitle("#it{E}_{cl} (GeV)");
    fHistClusterENMatchesCharged->GetYaxis()->SetTitle("#it{N}_{matches}");
    fHistograms->Add(fHistClusterENMatchesCharged);

    fHistClusterEvsTrackEPrimaryButNoElec        = new TH2F(Form("ClusterE_TrackE_ChargedClusterNoElec %s",GetCutNumber().Data()),"ClusterE TrackE ChargedClusterNoElec",
                                                             nBinsClusterEFine, minClusterELog, maxClusterELog, nBinsClusterEFine, minClusterELog, maxClusterELog);
    SetLogBinningXTH2(fHistClusterEvsTrackEPrimaryButNoElec);
    SetLogBinningYTH2(fHistClusterEvsTrackEPrimaryButNoElec);
    fHistClusterEvsTrackEPrimaryButNoElec->GetXaxis()->SetTitle("#it{E}_{cl} (GeV)");
    fHistClusterEvsTrackEPrimaryButNoElec->GetYaxis()->SetTitle("#it{E}_{tr} (GeV)");
    fHistograms->Add(fHistClusterEvsTrackEPrimaryButNoElec);
    fHistClusterEvsTrackSumEPrimaryButNoElec        = new TH2F(Form("ClusterE_TrackSumE_ChargedClusterNoElec %s",GetCutNumber().Data()),"ClusterE TrackSumE ChargedClusterNoElec",
                                                             nBinsClusterEFine, minClusterELog, maxClusterELog, nBinsClusterEFine, minClusterELog, maxClusterELog);
    SetLogBinningXTH2(fHistClusterEvsTrackSumEPrimaryButNoElec);
    SetLogBinningYTH2(fHistClusterEvsTrackSumEPrimaryButNoElec);
    fHistClusterEvsTrackSumEPrimaryButNoElec->GetXaxis()->SetTitle("#it{E}_{cl} (GeV)");
    fHistClusterEvsTrackSumEPrimaryButNoElec->GetYaxis()->SetTitle("#sum#it{E}_{tr} (GeV)");
    fHistograms->Add(fHistClusterEvsTrackSumEPrimaryButNoElec);

    if(fIsMC > 1){
      fHistClusterTMEffiInput->Sumw2();
      fHistClusterEvsTrackECharged->Sumw2();
      fHistClusterEvsTrackEChargedLead->Sumw2();
      fHistClusterEvsTrackENeutral->Sumw2();
      fHistClusterEvsTrackENeutralSubCharged->Sumw2();
      fHistClusterEvsTrackEGamma->Sumw2();
      fHistClusterEvsTrackEGammaSubCharged->Sumw2();
      fHistClusterEvsTrackEConv->Sumw2();
      fHistClusterENMatchesNeutral->Sumw2();
      fHistClusterENMatchesCharged->Sumw2();
      fHistClusterEvsTrackEPrimaryButNoElec->Sumw2();
      fHistClusterEvsTrackSumEPrimaryButNoElec->Sumw2();
    }
  }

  if (fDoExoticsQA){
    if( fClusterType == 1 ){ //EMCAL
      const Int_t nEmcalEtaBins             = 96;
      const Int_t nEmcalPhiBins             = 124;
      Float_t EmcalEtaBins[nEmcalEtaBins+1] = {-0.66687,-0.653,-0.63913,-0.62526,-0.61139,-0.59752,-0.58365,-0.56978,-0.55591,-0.54204,-0.52817,-0.5143,-0.50043,-0.48656,-0.47269,-0.45882,-0.44495,-0.43108,-0.41721,-0.40334,-0.38947,-0.3756,-0.36173,-0.34786,-0.33399,-0.32012,-0.30625,-0.29238,-0.27851,-0.26464,-0.25077,-0.2369,-0.22303,-0.20916,-0.19529,-0.18142,-0.16755,-0.15368,-0.13981,-0.12594,-0.11207,-0.0982,-0.08433,-0.07046,-0.05659,-0.04272,-0.02885,-0.01498,-0.00111,0.01276,0.02663,0.0405,0.05437,0.06824,0.08211,0.09598,0.10985,0.12372,0.13759,0.15146,0.16533,0.1792,0.19307,0.20694,0.22081,0.23468,0.24855,0.26242,0.27629,0.29016,0.30403,0.3179,0.33177,0.34564,0.35951,0.37338,0.38725,0.40112,0.41499,0.42886,0.44273,0.4566,0.47047,0.48434,0.49821,0.51208,0.52595,0.53982,0.55369,0.56756,0.58143,0.5953,0.60917,0.62304,0.63691,0.65078,0.66465};
      Float_t EmcalPhiBins[nEmcalPhiBins+1] = {1.408,1.4215,1.435,1.4485,1.462,1.4755,1.489,1.5025,1.516,1.5295,1.543,1.5565,1.57,1.5835,1.597,1.6105,1.624,1.6375,1.651,1.6645,1.678,1.6915,1.705,1.7185,1.732,1.758,1.7715,1.785,1.7985,1.812,1.8255,1.839,1.8525,1.866,1.8795,1.893,1.9065,1.92,1.9335,1.947,1.9605,1.974,1.9875,2.001,2.0145,2.028,2.0415,2.055,2.0685,2.082,2.108,2.1215,2.135,2.1485,2.162,2.1755,2.189,2.2025,2.216,2.2295,2.243,2.2565,2.27,2.2835,2.297,2.3105,2.324,2.3375,2.351,2.3645,2.378,2.3915,2.405,2.4185,2.432,2.456,2.4695,2.483,2.4965,2.51,2.5235,2.537,2.5505,2.564,2.5775,2.591,2.6045,2.618,2.6315,2.645,2.6585,2.672,2.6855,2.699,2.7125,2.726,2.7395,2.753,2.7665,2.78,2.804,2.8175,2.831,2.8445,2.858,2.8715,2.885,2.8985,2.912,2.9255,2.939,2.9525,2.966,2.9795,2.993,3.0065,3.02,3.0335,3.047,3.0605,3.074,3.0875,3.101,3.1145,3.128};

      fHistClusterEtavsPhiExotics     = new TH2F(Form("EtaPhi_Exotics %s",GetCutNumber().Data()),"EtaPhi_Exotics",nEmcalPhiBins,EmcalPhiBins,nEmcalEtaBins,EmcalEtaBins);
      fHistograms->Add(fHistClusterEtavsPhiExotics);
    } else if( fClusterType == 2 ){ //PHOS
      const Int_t nPhosEtaBins        = 56;
      const Int_t nPhosPhiBins        = 192;
      const Float_t PhosEtaRange[2]   = {-0.16, 0.16};
      const Float_t PhosPhiRange[2]   = {4.5, 5.6};

      fHistClusterEtavsPhiExotics     = new TH2F(Form("EtaPhi_Exotics %s",GetCutNumber().Data()),"EtaPhi_Exotics",nPhosPhiBins,PhosPhiRange[0],PhosPhiRange[1],nPhosEtaBins,PhosEtaRange[0],PhosEtaRange[1]);
      fHistograms->Add(fHistClusterEtavsPhiExotics);
    }
    if( fClusterType == 3 ){ //DCAL
      const Int_t nDcalEtaBins             = 96;
      const Int_t nDcalPhiBins             = 124;

      fHistClusterEtavsPhiExotics     = new TH2F(Form("EtaPhi_Exotics %s",GetCutNumber().Data()),"EtaPhi_Exotics",nDcalPhiBins,4.5,5.7,nDcalEtaBins,-0.66687,0.66465);
      fHistograms->Add(fHistClusterEtavsPhiExotics);
    } 
    fHistClusterEM02Exotics           = new TH2F(Form("EVsM02_Exotics %s",GetCutNumber().Data()),"EVsM02_afterClusterQA",nBinsClusterEFine,minClusterELog,maxClusterELog,200,0,2);
    SetLogBinningXTH2(fHistClusterEM02Exotics);
    fHistograms->Add(fHistClusterEM02Exotics);
    fHistClusterEnergyvsNCellsExotics = new TH2F(Form("ClusterEnergyVsNCells_Exotics %s",GetCutNumber().Data()),"ClusterEnergyVsNCells_Exotics",nBinsClusterECoarse,minClusterELog,maxClusterELog,50,0,50);
    SetLogBinningXTH2(fHistClusterEnergyvsNCellsExotics);
    fHistograms->Add(fHistClusterEnergyvsNCellsExotics);
    fHistClusterEEstarExotics         = new TH2F(Form("ClusterEnergyVsEnergystar_Exotics %s",GetCutNumber().Data()),"ClusterEnergyVsEnergystar_Exotics", nBinsClusterECoarse, minClusterELog, maxClusterELog,
                                                 nBinsClusterECoarse, minClusterELog, maxClusterELog);
    SetLogBinningXTH2(fHistClusterEEstarExotics);
    SetLogBinningYTH2(fHistClusterEEstarExotics);
    fHistograms->Add(fHistClusterEEstarExotics);
    
    if (fIsMC > 1){
      fHistClusterEtavsPhiExotics->Sumw2();
      fHistClusterEM02Exotics->Sumw2();
      fHistClusterEnergyvsNCellsExotics->Sumw2();
      fHistClusterEEstarExotics->Sumw2();
    }
  }
  
  fVectorMatchedClusterIDs.clear();

  TH1::AddDirectory(kTRUE);
  return;
}

//________________________________________________________________________
void AliCaloPhotonCuts::InitializeEMCAL(AliVEvent *event){

  if (fClusterType == 1 || fClusterType == 3){
    AliTender* alitender=0x0;
    AliEmcalTenderTask* emcaltender=0x0;
  
    if(event->IsA()==AliESDEvent::Class()) 
      alitender         = (AliTender*) AliAnalysisManager::GetAnalysisManager()->GetTopTasks()->FindObject("AliTender");
    else if( event->IsA()==AliAODEvent::Class()) 
      emcaltender       = (AliEmcalTenderTask*) AliAnalysisManager::GetAnalysisManager()->GetTopTasks()->FindObject("AliEmcalTenderTask");
    
    if(alitender){
      TIter next(alitender->GetSupplies());
      AliTenderSupply *supply;
      while ((supply=(AliTenderSupply*)next())) if(supply->IsA()==AliEMCALTenderSupply::Class()) break;
      fEMCALRecUtils        = ((AliEMCALTenderSupply*)supply)->GetRecoUtils();
      fEMCALBadChannelsMap  = fEMCALRecUtils->GetEMCALBadChannelStatusMapArray();
    } else if(emcaltender){
      fEMCALRecUtils        = ((AliEMCALTenderSupply*)emcaltender->GetEMCALTenderSupply())->GetRecoUtils();
      fEMCALBadChannelsMap  = fEMCALRecUtils->GetEMCALBadChannelStatusMapArray();
    }
    if (fEMCALRecUtils) fEMCALInitialized = kTRUE;

    fGeomEMCAL              = AliEMCALGeometry::GetInstance();
    if(!fGeomEMCAL){ AliFatal("EMCal geometry not initialized!");}

    //retrieve pointer to trackMatcher Instance
    if(fUseDistTrackToCluster) fCaloTrackMatcher = (AliCaloTrackMatcher*)AliAnalysisManager::GetAnalysisManager()->GetTask(fCaloTrackMatcherName.Data());
    if(!fCaloTrackMatcher && fUseDistTrackToCluster ){ AliFatal("CaloTrackMatcher instance could not be initialized!");}

    if(!fDoLightOutput){
      Int_t nMaxCellsEMCAL  = fNMaxEMCalModules*48*24;
      Int_t nMinCellsDCAL = 12288;
      Int_t nMaxCellsDCAL = nMinCellsDCAL+fNMaxDCalModules*32*24;
      Int_t nMaxCells;
      Int_t nMinCells;
      if(fClusterType == 1){
        nMaxCells = nMaxCellsEMCAL;
        nMinCells = 0;
      } else if(fClusterType == 3){
        nMaxCells = nMaxCellsDCAL;
        nMinCells = nMinCellsDCAL;
      }


      Int_t imod = -1;Int_t iTower = -1, iIphi = -1, iIeta = -1;
      Int_t icol = -1;Int_t irow = -1;

      for(Int_t iCell=nMinCells;iCell<nMaxCells;iCell++){
        fGeomEMCAL->GetCellIndex(iCell,imod,iTower,iIphi,iIeta);
        if (fEMCALBadChannelsMap->GetEntries() <= imod) continue;
        fGeomEMCAL->GetCellPhiEtaIndexInSModule(imod,iTower,iIphi,iIeta,irow,icol);
        Int_t iBadCell      = (Int_t) ((TH2I*)fEMCALBadChannelsMap->At(imod))->GetBinContent(icol,irow);
        if(iBadCell > 0) fBadChannels->Fill(iCell,1);
        else fBadChannels->Fill(iCell,0);
      }
    }
  }
  return;
}

//________________________________________________________________________
void AliCaloPhotonCuts::InitializePHOS (AliVEvent *event){

  if (fClusterType == 2){
    if(!fDoLightOutput){
      fGeomPHOS = AliPHOSGeometry::GetInstance();
      if(!fGeomPHOS) AliFatal("PHOS geometry not initialized!");
      Int_t nModules = fGeomPHOS->GetNModules();
      //cout << nModules << endl;

      fPHOSBadChannelsMap = new TH2I*[nModules];

      AliOADBContainer badmapContainer(Form("phosBadMap"));
      badmapContainer.InitFromFile("$ALICE_PHYSICS/OADB/PHOS/PHOSBadMaps.root","phosBadMap");
      TObjArray *maps = (TObjArray*)badmapContainer.GetObject(event->GetRunNumber(),"phosBadMap");

      if(!maps){
        AliError(Form("Can not read Bad map for run %d. \n You may choose to use your map with ForceUsingBadMap()\n",event->GetRunNumber())) ;
        for(Int_t mod=0;mod<nModules;mod++) fPHOSBadChannelsMap[mod] = NULL;
      }else{
        AliInfo(Form("Setting PHOS bad map with name %s \n",maps->GetName())) ;
        for(Int_t mod=0;mod<nModules;mod++){
          TH2I * h = (TH2I*)maps->At(mod);
          //cout << mod << ", " << h << ", " << __LINE__ << endl;
          if(h) fPHOSBadChannelsMap[mod] = new TH2I(*h);
          else fPHOSBadChannelsMap[mod] = NULL;
        }

        Int_t nMaxCellsPHOS = fNMaxPHOSModules*56*64;
        Int_t relid[4];

        for(Int_t iCell=0;iCell<nMaxCellsPHOS;iCell++){
          fGeomPHOS->AbsToRelNumbering(iCell,relid);
          //cout << relid[0] << ", "  << relid[1] << ", "  << relid[2] << ", "  << relid[3] << ", " << __LINE__ << endl;
          if(relid[1]!=0) AliFatal("PHOS CPV in PHOS cell array?");
          if(relid[0]>=nModules || relid[0]<0 || !fPHOSBadChannelsMap[relid[0]]) continue;
          Int_t iBadCell = (Int_t) ((TH2I*)fPHOSBadChannelsMap[relid[0]])->GetBinContent(relid[2],relid[3]);
          if(iBadCell > 0) fBadChannels->Fill(iCell,1);
          else fBadChannels->Fill(iCell,0);
        }
      }
    }

    //retrieve pointer to trackMatcher Instance
    if(fUseDistTrackToCluster) fCaloTrackMatcher = (AliCaloTrackMatcher*)AliAnalysisManager::GetAnalysisManager()->GetTask(fCaloTrackMatcherName.Data());
    if(!fCaloTrackMatcher && fUseDistTrackToCluster){ AliFatal("CaloTrackMatcher instance could not be initialized!");}

    fPHOSInitialized = kTRUE;
    fPHOSCurrentRun = event->GetRunNumber();
  }
  return;
}

//________________________________________________________________________
Bool_t AliCaloPhotonCuts::ClusterIsSelectedMC(TParticle *particle,AliStack *fMCStack){
   // MonteCarlo Photon Selection

  if(!fMCStack)return kFALSE;
  if(!particle) return kFALSE;

  if (particle->GetPdgCode() == 22){

    if ( particle->Eta() < fMinEtaCut || particle->Eta() > fMaxEtaCut ) return kFALSE;
    if ( particle->Phi() < fMinPhiCut || particle->Phi() > fMaxPhiCut ) return kFALSE;
    
    if(particle->GetMother(0) >-1 && fMCStack->Particle(particle->GetMother(0))->GetPdgCode() == 22){
      return kFALSE;// no photon as mothers!
    }
    return kTRUE;
  }
  return kFALSE;
}

//________________________________________________________________________
Bool_t AliCaloPhotonCuts::ClusterIsSelectedElecMC(TParticle *particle,AliStack *fMCStack){
   // MonteCarlo Photon Selection

  if(!fMCStack)return kFALSE;
  if(!particle) return kFALSE;

  if (TMath::Abs(particle->GetPdgCode()) == 11){

    if ( particle->Eta() < fMinEtaCut || particle->Eta() > fMaxEtaCut ) return kFALSE;
    if ( particle->Phi() < fMinPhiCut || particle->Phi() > fMaxPhiCut ) return kFALSE;
    
    if(particle->GetMother(0) >-1 && fMCStack->Particle(particle->GetMother(0))->GetPdgCode() == 11){
      return kFALSE;// no photon as mothers!
    }
    return kTRUE;
  }
  return kFALSE;
}

//________________________________________________________________________
Bool_t AliCaloPhotonCuts::ClusterIsSelectedElecAODMC(AliAODMCParticle *particle,TClonesArray *aodmcArray){
   // MonteCarlo Photon Selection

  if(!aodmcArray)return kFALSE;
  if(!particle) return kFALSE;

  if (TMath::Abs(particle->GetPdgCode()) == 11){

    if ( particle->Eta() < fMinEtaCut || particle->Eta() > fMaxEtaCut ) return kFALSE;
    if ( particle->Phi() < fMinPhiCut || particle->Phi() > fMaxPhiCut ) return kFALSE;
    if(particle->GetMother() >-1 && (static_cast<AliAODMCParticle*>(aodmcArray->At(particle->GetMother())))->GetPdgCode() == 11){
      return kFALSE;// no photon as mothers!
    }
    return kTRUE;
  }
  return kFALSE;
}

//________________________________________________________________________
Bool_t AliCaloPhotonCuts::ClusterIsSelectedAODMC(AliAODMCParticle *particle,TClonesArray *aodmcArray){
  // MonteCarlo Photon Selection

  if(!aodmcArray)return kFALSE;
  if(!particle) return kFALSE;

  if (particle->GetPdgCode() == 22){
    if ( particle->Eta() < fMinEtaCut || particle->Eta() > fMaxEtaCut ) return kFALSE;
    if ( particle->Phi() < fMinPhiCut || particle->Phi() > fMaxPhiCut ) return kFALSE;
    if(particle->GetMother() > -1 && (static_cast<AliAODMCParticle*>(aodmcArray->At(particle->GetMother())))->GetPdgCode() == 22){
        return kFALSE;// no photon as mothers!
    }
    return kTRUE;// return in case of accepted gamma
  }
  return kFALSE;
}

//________________________________________________________________________
// This function selects the clusters based on their quality criteria
//________________________________________________________________________
Bool_t AliCaloPhotonCuts::ClusterQualityCuts(AliVCluster* cluster, AliVEvent *event, AliVEvent* mcEvent, Int_t isMC, Double_t weight, Long_t clusterID)
{   // Specific Photon Cuts

  // Initialize EMCAL rec utils if not initialized
  if(!fEMCALInitialized && (fClusterType == 1 || fClusterType == 3)) InitializeEMCAL(event);

  fIsCurrentClusterAcceptedBeforeTM = kFALSE;
  Int_t cutIndex = 0;
  if(fHistClusterIdentificationCuts)fHistClusterIdentificationCuts->Fill(cutIndex,cluster->E());
  cutIndex++;

  // cluster position defintion
  Float_t clusPos[3]={0,0,0};
  cluster->GetPosition(clusPos);
  TVector3 clusterVector(clusPos[0],clusPos[1],clusPos[2]);
  Double_t etaCluster = clusterVector.Eta();
  Double_t phiCluster = clusterVector.Phi();
  if (phiCluster < 0) phiCluster += 2*TMath::Pi();

  
  Int_t nLM = GetNumberOfLocalMaxima(cluster, event);
//   Int_t nLMGustavo = fEMCALCaloUtils->GetNumberOfLocalMaxima(cluster, event->GetEMCALCells()) ;
//   cout << "mine: " << nLM << "\t Gustavo: " << nLMGustavo << endl;
  
  // Fill Histos before Cuts
  if(fHistClusterTimevsEBeforeQA) fHistClusterTimevsEBeforeQA->Fill(cluster->GetTOF(), cluster->E(), weight);
  if(fHistEnergyOfClusterBeforeQA) fHistEnergyOfClusterBeforeQA->Fill(cluster->E(), weight);
  if(fHistNCellsBeforeQA) fHistNCellsBeforeQA->Fill(cluster->GetNCells(), weight);
  if(fHistM02BeforeQA) fHistM02BeforeQA->Fill(cluster->GetM02(), weight);
  if(fHistM20BeforeQA) fHistM20BeforeQA->Fill(cluster->GetM20(), weight);
  if(fHistDispersionBeforeQA) fHistDispersionBeforeQA->Fill(cluster->GetDispersion(), weight);
  if(fHistNLMBeforeQA) fHistNLMBeforeQA->Fill(nLM, weight);
//   if(fHistNLMAvsNLMBBeforeQA) fHistNLMAvsNLMBBeforeQA->Fill(nLM, nLMGustavo, weight);
  if(fHistClusterEM02BeforeQA) fHistClusterEM02BeforeQA->Fill(cluster->E(),cluster->GetM02(), weight);

  AliVCaloCells* cells = NULL;
  if(fExtendedMatchAndQA > 1){
    if(cluster->IsEMCAL()){ //EMCAL
      cells = event->GetEMCALCells();
    }else if(cluster->IsPHOS()){ //PHOS
      cells = event->GetPHOSCells();
    }
    if(fHistClusterIncludedCellsBeforeQA){
      Int_t nCellCluster = cluster->GetNCells();
      for(Int_t iCell=0;iCell<nCellCluster;iCell++){
        fHistClusterIncludedCellsBeforeQA->Fill(cluster->GetCellAbsId(iCell));
        if(cluster->E()>0.) fHistClusterEnergyFracCellsBeforeQA->Fill(cluster->GetCellAbsId(iCell),cells->GetCellAmplitude(cluster->GetCellAbsId(iCell))/cluster->E());
      }
    }
  }

  // Check wether timing is ok
  if (fUseTimeDiff){
    if( (cluster->GetTOF() < fMinTimeDiff || cluster->GetTOF() > fMaxTimeDiff) && !(isMC>0)){
      if(fHistClusterIdentificationCuts)fHistClusterIdentificationCuts->Fill(cutIndex, cluster->E());//1
      return kFALSE;
    }
  }  
  cutIndex++;//2, next cut

  // exotic cluster cut 
  Float_t energyStar      = 0;
  if(fUseExoticCluster && IsExoticCluster(cluster, event, energyStar)){
    if(fHistClusterIdentificationCuts)fHistClusterIdentificationCuts->Fill(cutIndex, cluster->E());//3
    if (fDoExoticsQA){
      // replay cuts
      Bool_t failed     = kFALSE;
      Bool_t failedM02  = kFALSE;
      if (fUseMinEnergy)
        if(cluster->E() < fMinEnergy)
          failed = kTRUE;
      if (fUseNCells)
        if(cluster->GetNCells() < fMinNCells) 
          failed = kTRUE;
      if (fUseNLM)
        if( nLM < fMinNLM || nLM > fMaxNLM ) 
          failed = kTRUE;
      if (fUseM02 == 1){
        if( cluster->GetM02()< fMinM02 || cluster->GetM02() > fMaxM02 ) 
          failedM02  = kTRUE;
      } else if (fUseM02 ==2 ) {
        if( cluster->GetM02()< CalculateMinM02(fMinM02CutNr, cluster->E()) || 
            cluster->GetM02() > CalculateMaxM02(fMaxM02CutNr, cluster->E()) ) 
          failedM02  = kTRUE;
      }  
      if (fUseM20)
        if( cluster->GetM20()< fMinM20 || cluster->GetM20() > fMaxM20 ) 
          failed = kTRUE;
      if (fUseDispersion)
        if( cluster->GetDispersion()> fMaxDispersion) 
          failed = kTRUE;
      if (fVectorMatchedClusterIDs.size()>0 && fUseDistTrackToCluster)
        if( CheckClusterForTrackMatch(cluster) )
          failed = kTRUE;
      if ( !( failed || failedM02 ) ){  
        if(fHistClusterEtavsPhiExotics) fHistClusterEtavsPhiExotics->Fill(phiCluster, etaCluster, weight);
        if(fHistClusterEnergyvsNCellsExotics) fHistClusterEnergyvsNCellsExotics->Fill(cluster->E(), cluster->GetNCells(), weight);
        if(fHistClusterEEstarExotics) fHistClusterEEstarExotics->Fill(cluster->E(),energyStar, weight);
      }
      if ( !failed ){
        if(fHistClusterEM02Exotics) fHistClusterEM02Exotics->Fill(cluster->E(), cluster->GetM02(), weight);
      }  
    }  
    return kFALSE;
  }
  cutIndex++;//3, next cut
  
  // minimum cell energy cut
  if (fUseMinEnergy){
    if(cluster->E() < fMinEnergy){
      if(fHistClusterIdentificationCuts)fHistClusterIdentificationCuts->Fill(cutIndex, cluster->E());//4
      return kFALSE;
    }
  }  
  cutIndex++;//4, next cut

  // minimum number of cells
  if (fUseNCells){
    if(cluster->GetNCells() < fMinNCells) {
      if(fHistClusterIdentificationCuts)fHistClusterIdentificationCuts->Fill(cutIndex, cluster->E());//5
      return kFALSE;
    }
  }  
  cutIndex++;//5, next cut

  // NLM cut
  if (fUseNLM){
    if( nLM < fMinNLM || nLM > fMaxNLM ) {
      if(fHistClusterIdentificationCuts)fHistClusterIdentificationCuts->Fill(cutIndex, cluster->E());//9
      return kFALSE;
    }
  }  
  cutIndex++;//6, next cut
  
  // M02 cut
  if (fUseM02 == 1){
    if( cluster->GetM02()< fMinM02 || cluster->GetM02() > fMaxM02 ) {
      if(fHistClusterIdentificationCuts)fHistClusterIdentificationCuts->Fill(cutIndex, cluster->E());//6
      return kFALSE;
    }
  } else if (fUseM02 ==2 ) {
    if( cluster->GetM02()< CalculateMinM02(fMinM02CutNr, cluster->E()) || 
      cluster->GetM02() > CalculateMaxM02(fMaxM02CutNr, cluster->E()) ) {
      if(fHistClusterIdentificationCuts)fHistClusterIdentificationCuts->Fill(cutIndex, cluster->E());//6
      return kFALSE;
    }
  }
  cutIndex++;//7, next cut

  // M20 cut
  if (fUseM20){
    if( cluster->GetM20()< fMinM20 || cluster->GetM20() > fMaxM20 ) {
      if(fHistClusterIdentificationCuts)fHistClusterIdentificationCuts->Fill(cutIndex, cluster->E());//7
      return kFALSE;
    }
  }  
  cutIndex++;//8, next cut
  
  // dispersion cut
  if (fUseDispersion){
    if( cluster->GetDispersion()> fMaxDispersion) {
      if(fHistClusterIdentificationCuts)fHistClusterIdentificationCuts->Fill(cutIndex, cluster->E());//8
      return kFALSE;
    }
  }  
  cutIndex++;//9, next cut


  // Classify events
  AliESDEvent *esdev  = dynamic_cast<AliESDEvent*>(event);
  AliAODEvent *aodev  = 0;
  Bool_t isESD        = kTRUE;
  if (!esdev) {
    isESD             = kFALSE;
    aodev             = dynamic_cast<AliAODEvent*>(event);
    if (!aodev) {
      AliError("Task needs AOD or ESD event, returning");
      return kFALSE;
    }
  }

  fIsCurrentClusterAcceptedBeforeTM = kTRUE;
  
  // classification of clusters for TM efficiency
  // 0: Neutral cluster
  // 1: Neutral cluster sub charged
  // 2: Gamma cluster
  // 3: Gamma cluster sub charged
  // 4: Gamma conv cluster
  // 5: Charged cluster
  // 6: Electron
  // 7: prim charged cluster
  Int_t classification  = -1;
  Long_t leadMCLabel    = -1;
  if (isESD) 
    leadMCLabel         = ((AliESDCaloCluster*)cluster)->GetLabel();
  else 
    leadMCLabel         = ((AliAODCaloCluster*)cluster)->GetLabel();
  
  // TM efficiency histograms before TM
  if (fIsMC && isMC && (fExtendedMatchAndQA == 1 || fExtendedMatchAndQA == 3 || fExtendedMatchAndQA == 5) && fUseDistTrackToCluster){
    classification    = ClassifyClusterForTMEffi(cluster, event, mcEvent, isESD); 
    fHistClusterTMEffiInput->Fill(cluster->E(), 0., weight); //All cl 
    if (classification == 5 )
      fHistClusterTMEffiInput->Fill(cluster->E(), 1., weight); //Ch cl 
    if (classification == 7 )
      fHistClusterTMEffiInput->Fill(cluster->E(), 7., weight); //Ch cl 
    if (classification == 4)
      fHistClusterTMEffiInput->Fill(cluster->E(), 6., weight); //conv electron cl 
    if (classification == 6)
      fHistClusterTMEffiInput->Fill(cluster->E(), 8., weight); // electron cl 
    if (classification == 0 || classification == 1)
      fHistClusterTMEffiInput->Fill(cluster->E(), 2., weight); // Ne cl match
    if (classification == 1)
      fHistClusterTMEffiInput->Fill(cluster->E(), 3., weight); // Ne cl sub ch match
    if (classification == 2 || classification == 3)
      fHistClusterTMEffiInput->Fill(cluster->E(), 4., weight); // Ga cl match
    if ( classification == 3)
      fHistClusterTMEffiInput->Fill(cluster->E(), 5., weight); // Ga cl sub ch match      

    Int_t nlabelsMatchedTracks      = 0;
    if (!fUsePtDepTrackToCluster)
      nlabelsMatchedTracks          = fCaloTrackMatcher->GetNMatchedTrackIDsForCluster(event, cluster->GetID(), fMaxDistTrackToClusterEta, -fMaxDistTrackToClusterEta, 
                                                                                      fMaxDistTrackToClusterPhi, fMinDistTrackToClusterPhi);
    else 
      nlabelsMatchedTracks          = fCaloTrackMatcher->GetNMatchedTrackIDsForCluster(event, cluster->GetID(), fFuncPtDepEta, fFuncPtDepPhi);
    if (classification < 4 && classification > -1)
      fHistClusterENMatchesNeutral->Fill(cluster->E(), nlabelsMatchedTracks);
    else 
      fHistClusterENMatchesCharged->Fill(cluster->E(), nlabelsMatchedTracks);
  }  
  
  
  if (fVectorMatchedClusterIDs.size()>0 && fUseDistTrackToCluster){
    if( CheckClusterForTrackMatch(cluster) ){
      if(fHistClusterIdentificationCuts)fHistClusterIdentificationCuts->Fill(cutIndex, cluster->E());//2
      // TM efficiency histos after TM
      if (fIsMC && isMC && (fExtendedMatchAndQA == 1 || fExtendedMatchAndQA == 3 || fExtendedMatchAndQA == 5)){

        fHistClusterTMEffiInput->Fill(cluster->E(), 9., weight);
        if (classification == 5 )
          fHistClusterTMEffiInput->Fill(cluster->E(), 10., weight); //Ch cl match
        if (classification == 4)
          fHistClusterTMEffiInput->Fill(cluster->E(), 16., weight); //conv cl match
        if (classification == 0 || classification == 1)
          fHistClusterTMEffiInput->Fill(cluster->E(), 12., weight); // Ne cl match
        if ( classification == 1)
          fHistClusterTMEffiInput->Fill(cluster->E(), 13., weight); // Ne cl sub ch match
        if (classification == 2 || classification == 3)
          fHistClusterTMEffiInput->Fill(cluster->E(), 14., weight); // Ga cl match
        if ( classification == 3)
          fHistClusterTMEffiInput->Fill(cluster->E(), 15., weight); // Ga cl sub ch match
        if ( classification == 7)
          fHistClusterTMEffiInput->Fill(cluster->E(), 18., weight); // Ch prim cl match
        if ( classification == 6)
          fHistClusterTMEffiInput->Fill(cluster->E(), 20., weight); // El cl match
          
        vector<Int_t> labelsMatchedTracks;
        if (!fUsePtDepTrackToCluster)
          labelsMatchedTracks           = fCaloTrackMatcher->GetMatchedTrackIDsForCluster(event, cluster->GetID(), fMaxDistTrackToClusterEta, -fMaxDistTrackToClusterEta, 
                                                                                          fMaxDistTrackToClusterPhi, fMinDistTrackToClusterPhi);
        else 
          labelsMatchedTracks           = fCaloTrackMatcher->GetMatchedTrackIDsForCluster(event, cluster->GetID(), fFuncPtDepEta, fFuncPtDepPhi);
        
        Int_t idHighestPt = -1;
        Double_t ptMax    = -1;
        Double_t eMax     = -1;
        Double_t eSum     = 0;
        Bool_t foundLead  = kFALSE;
        Double_t eLead    = -1;
        Int_t idLead      = -1;
        for (Int_t i = 0; i < (Int_t)labelsMatchedTracks.size(); i++){
          AliVTrack* currTrack  = dynamic_cast<AliVTrack*>(event->GetTrack(labelsMatchedTracks.at(i)));
          eSum += currTrack->E();
          if (ptMax < currTrack->Pt()){
            ptMax               = currTrack->Pt();
            eMax                = currTrack->E();
            idHighestPt         = labelsMatchedTracks.at(i);
          }
          if (classification == 4 || classification == 5 || classification == 6 || classification == 7){
            Long_t mcLabelTrack = -1;
            if (isESD)
              mcLabelTrack      = TMath::Abs(((AliESDtrack*)currTrack)->GetLabel());
            else  
              mcLabelTrack      = TMath::Abs(((AliAODTrack*)currTrack)->GetLabel());
            if (mcLabelTrack!= -1 && mcLabelTrack == leadMCLabel){
              foundLead         = kTRUE;
              eLead             = currTrack->E();
              idLead            = labelsMatchedTracks.at(i);
            }  
          }  
        }
        if (classification == 5 || classification == 7 || classification == 6){
          fHistClusterEvsTrackECharged->Fill(cluster->E(), eMax, weight);
          if (classification == 5 || classification == 7){
            fHistClusterEvsTrackEPrimaryButNoElec->Fill(cluster->E(), eMax, weight);
            fHistClusterEvsTrackSumEPrimaryButNoElec->Fill(cluster->E(), eSum, weight);
          }

          if (foundLead ){
            if (classification == 5)
              fHistClusterTMEffiInput->Fill(cluster->E(), 11., weight); //Ch cl match w lead
            if (classification == 7)
              fHistClusterTMEffiInput->Fill(cluster->E(), 19., weight); //Ch prim cl match w lead
            if (classification == 6)
              fHistClusterTMEffiInput->Fill(cluster->E(), 21., weight); //El cl match w lead
            fHistClusterEvsTrackEChargedLead->Fill(cluster->E(), eLead, weight);
          }
        }    
        if (classification == 4){
          fHistClusterEvsTrackEConv->Fill(cluster->E(), eMax, weight);
          if (foundLead)
          fHistClusterTMEffiInput->Fill(cluster->E(), 17., weight); //conv cl match w lead
        }
        if (classification == 0 )
          fHistClusterEvsTrackENeutral->Fill(cluster->E(), eMax, weight);    
        if (classification == 1)  
          fHistClusterEvsTrackENeutralSubCharged->Fill(cluster->E(), eMax, weight);
        if (classification == 2)
          fHistClusterEvsTrackEGamma->Fill(cluster->E(), eMax, weight);    
        if (classification == 3)
          fHistClusterEvsTrackEGammaSubCharged->Fill(cluster->E(), eMax, weight);    
        
        labelsMatchedTracks.clear();
      }

      return kFALSE;
    }
  }


  cutIndex++;//10, next cut
  
  // DONE with selecting photons
  if(fHistClusterIdentificationCuts)fHistClusterIdentificationCuts->Fill(cutIndex, cluster->E());//10

  // Histos after Cuts

  if(fHistClusterEtavsPhiAfterQA) fHistClusterEtavsPhiAfterQA->Fill(phiCluster, etaCluster, weight);
  if(fHistClusterTimevsEAfterQA) fHistClusterTimevsEAfterQA->Fill(cluster->GetTOF(), cluster->E(), weight);
  if(fHistEnergyOfClusterAfterQA) fHistEnergyOfClusterAfterQA->Fill(cluster->E(), weight);
  if(fHistNCellsAfterQA) fHistNCellsAfterQA->Fill(cluster->GetNCells(), weight);
  if(fHistM02AfterQA) fHistM02AfterQA->Fill(cluster->GetM02(), weight);
  if(fHistM20AfterQA) fHistM20AfterQA->Fill(cluster->GetM20(), weight);
  if(fHistDispersionAfterQA) fHistDispersionAfterQA->Fill(cluster->GetDispersion(), weight);
  if(fHistNLMAfterQA) fHistNLMAfterQA->Fill(nLM, weight);
  if(fHistNLMVsNCellsAfterQA) fHistNLMVsNCellsAfterQA->Fill(nLM,cluster->GetNCells(), weight);
  if(fHistNLMVsEAfterQA) fHistNLMVsEAfterQA->Fill(nLM, cluster->E(), weight);
  if(fHistClusterEM02AfterQA) fHistClusterEM02AfterQA->Fill(cluster->E(), cluster->GetM02(), weight);
  
  if(fExtendedMatchAndQA > 1){
    if(fHistClusterIncludedCellsAfterQA){
      Int_t nCellCluster = cluster->GetNCells();
      for(Int_t iCell=0;iCell<nCellCluster;iCell++){
        Int_t cellID = cluster->GetCellAbsId(iCell);
        Double_t cellAmp = cells->GetCellAmplitude(cellID);
        Double_t cellTime = cells->GetCellTime(cellID);
        fHistClusterIncludedCellsAfterQA->Fill(cellID);
        if(cluster->E()>0.) fHistClusterEnergyFracCellsAfterQA->Fill(cellID,cellAmp/cluster->E());
        if(isMC==0){
          fHistClusterIncludedCellsTimingAfterQA->Fill(cluster->E(),cellTime*1E9);
          fHistClusterIncludedCellsTimingEnergyAfterQA->Fill(cellAmp,cellTime*1E9);
        }else{
          fHistClusterIncludedCellsTimingAfterQA->Fill(cluster->E(),cellTime*1E8);
          fHistClusterIncludedCellsTimingEnergyAfterQA->Fill(cellAmp,cellTime*1E8);
        }
      }
    }
    
    if(fHistClusterEnergyvsNCells) fHistClusterEnergyvsNCells->Fill(cluster->E(),cluster->GetNCells());
    if(cluster->IsEMCAL()){
      Int_t iSuperModule = -1;
      fGeomEMCAL = AliEMCALGeometry::GetInstance();
      if(!fGeomEMCAL){ AliFatal("EMCal geometry not initialized!");}
      if(fHistClusterEnergyvsMod && fGeomEMCAL->SuperModuleNumberFromEtaPhi(clusterVector.Eta(),clusterVector.Phi(),iSuperModule)){
        fHistClusterEnergyvsMod->Fill(cluster->E(),iSuperModule);
      }
    }else if(cluster->IsPHOS()){
      Int_t relId[4] = {0,0,0,0};
      fGeomPHOS = AliPHOSGeometry::GetInstance();
      if(!fGeomPHOS){ AliFatal("PHOS geometry not initialized!");}
      if(fHistClusterEnergyvsMod && fGeomPHOS->GlobalPos2RelId(clusterVector,relId)){
        fHistClusterEnergyvsMod->Fill(cluster->E(),relId[0]);
      }
    }
  }

  return kTRUE;
}



//________________________________________________________________________
void AliCaloPhotonCuts::FillHistogramsExtendedQA(AliVEvent *event, Int_t isMC)
{
  if(fExtendedMatchAndQA < 2) return;

  AliVCaloCells* cells;

  Int_t nModules = 0;
  Int_t* nCellsBigger100MeV;
  Int_t* nCellsBigger1500MeV;
  Double_t* EnergyOfMod;

  if( (fClusterType == 1 || fClusterType == 3) && !fEMCALInitialized ) InitializeEMCAL(event);
  if( fClusterType == 2 && ( !fPHOSInitialized || (fPHOSCurrentRun != event->GetRunNumber()) ) ) InitializePHOS(event);
    
  Int_t nModulesStart = 0;
  if( fClusterType == 1 || fClusterType == 3){ //EMCAL & DCAL
    cells = event->GetEMCALCells();
    fGeomEMCAL = AliEMCALGeometry::GetInstance();
    if(!fGeomEMCAL) AliFatal("EMCal geometry not initialized!");
    if(!fEMCALBadChannelsMap) AliFatal("EMCal bad channels map not initialized!");
    nModules = fGeomEMCAL->GetNumberOfSuperModules();
    if( fClusterType == 3) {nModules = 8; nModulesStart = 12;}
  } else if( fClusterType == 2 ){ //PHOS
    cells = event->GetPHOSCells();
    fGeomPHOS = AliPHOSGeometry::GetInstance();
    if(!fGeomPHOS) AliFatal("PHOS geometry not initialized!");
    nModules = fGeomPHOS->GetNModules();
  } else{
    AliError(Form("fExtendedMatchAndQA(%i):FillHistogramsExtendedMatchAndQA() not (yet) defined for cluster type (%i)",fExtendedMatchAndQA,fClusterType));
  }
  
  nCellsBigger100MeV = new Int_t[nModules];
  nCellsBigger1500MeV = new Int_t[nModules];
  EnergyOfMod = new Double_t[nModules];

  for(Int_t iModule=0;iModule<nModules;iModule++){nCellsBigger100MeV[iModule]=0;nCellsBigger1500MeV[iModule]=0;EnergyOfMod[iModule]=0;}

  for(Int_t iCell=0;iCell<cells->GetNumberOfCells();iCell++){
    Short_t cellNumber=0;
    Double_t cellAmplitude=0;
    Double_t cellTime=0;
    Double_t cellEFrac=0;
    Int_t cellMCLabel=0;
    Int_t nMod = -1;

    cells->GetCell(iCell,cellNumber,cellAmplitude,cellTime,cellMCLabel,cellEFrac);
    if( fClusterType == 3 && cellNumber < 12288){continue;}
    if( fClusterType == 2 && cellNumber < 0){continue;} //Scip CPV cells in PHOS case
    Int_t imod = -1;Int_t iTower = -1, iIphi = -1, iIeta = -1;
    Int_t icol = -1;Int_t irow = -1;
    Int_t relid[4];// for PHOS

    Bool_t doBadCell = kTRUE;
    if( fClusterType == 1 || fClusterType == 3){
      nMod = fGeomEMCAL->GetSuperModuleNumber(cellNumber);
      fGeomEMCAL->GetCellIndex(cellNumber,imod,iTower,iIphi,iIeta);
      if (fEMCALBadChannelsMap->GetEntries() <= imod) doBadCell=kFALSE;
      fGeomEMCAL->GetCellPhiEtaIndexInSModule(imod,iTower,iIphi,iIeta,irow,icol);
    }else if( fClusterType == 2 ){
      fGeomPHOS->AbsToRelNumbering(cellNumber,relid);
      if(relid[1]!=0) AliFatal("PHOS CPV in PHOS cell array?");
      nMod = relid[0];//(Int_t) (1 + (cellNumber-1)/3584);
      if(nMod>=nModules || nMod<0 || !fPHOSBadChannelsMap[nMod]) doBadCell=kFALSE;
    }

    Int_t iBadCell = 0;
    if( (fClusterType == 1 || fClusterType == 3) && doBadCell){
      iBadCell = (Int_t) ((TH2I*)fEMCALBadChannelsMap->At(imod))->GetBinContent(icol,irow);
    }else if( fClusterType == 2 && doBadCell){
      iBadCell = (Int_t) ((TH2I*)fPHOSBadChannelsMap[nMod])->GetBinContent(relid[2],relid[3]);
    }

    if(iBadCell > 0) continue;
	// nModulesStart == 0 for EMCAL and PHOS
    if(cellAmplitude > 0.1) nCellsBigger100MeV[nMod-nModulesStart]++; 
    if(cellAmplitude > 1.5) nCellsBigger1500MeV[nMod-nModulesStart]++;
    if(cellAmplitude > 0.05) EnergyOfMod[nMod-nModulesStart]+=cellAmplitude;
      
    if(fExtendedMatchAndQA > 3){
      if(fHistCellEnergyvsCellID && (cellAmplitude > 0.05)) fHistCellEnergyvsCellID->Fill(cellAmplitude,cellNumber);
      if(fHistCellTimevsCellID && (cellAmplitude > 0.2)) fHistCellTimevsCellID->Fill(cellTime,cellNumber);
    }
  }
  
  for(Int_t iModule=0;iModule<nModules;iModule++){
    if(fHistNCellsBigger100MeVvsMod) fHistNCellsBigger100MeVvsMod->Fill(nCellsBigger100MeV[iModule],iModule+nModulesStart);
    if(fHistNCellsBigger1500MeVvsMod) fHistNCellsBigger1500MeVvsMod->Fill(nCellsBigger1500MeV[iModule],iModule+nModulesStart);
    if(fHistEnergyOfModvsMod) fHistEnergyOfModvsMod->Fill(EnergyOfMod[iModule],iModule+nModulesStart);
  }

  delete[] nCellsBigger100MeV;nCellsBigger100MeV=0x0;
  delete[] nCellsBigger1500MeV;nCellsBigger1500MeV=0x0;
  delete[] EnergyOfMod;EnergyOfMod=0x0;

  //fill distClusterTo_withinTiming/outsideTiming
  Int_t nclus = event->GetNumberOfCaloClusters();
  AliVCluster* cluster = 0x0;
  AliVCluster* clusterMatched = 0x0;
  for(Int_t iClus=0; iClus<nclus ; iClus++){
    if(event->IsA()==AliESDEvent::Class()) cluster = new AliESDCaloCluster(*(AliESDCaloCluster*)event->GetCaloCluster(iClus));
    else if(event->IsA()==AliAODEvent::Class()) cluster = new AliAODCaloCluster(*(AliAODCaloCluster*)event->GetCaloCluster(iClus));

    if( (fClusterType == 1 || fClusterType == 3) && !cluster->IsEMCAL()){delete cluster; continue;}
    if( fClusterType == 2 && cluster->GetType() !=AliVCluster::kPHOSNeutral){delete cluster; continue;}

    Float_t clusPos[3]={0,0,0};
    cluster->GetPosition(clusPos);
    TVector3 clusterVector(clusPos[0],clusPos[1],clusPos[2]);
    Double_t etaCluster = clusterVector.Eta();
    Double_t phiCluster = clusterVector.Phi();
    if (phiCluster < 0) phiCluster += 2*TMath::Pi();
    Int_t nLM = GetNumberOfLocalMaxima(cluster, event);

    //acceptance cuts
    if (fUseEtaCut && (etaCluster < fMinEtaCut || etaCluster > fMaxEtaCut)){delete cluster; continue;}
    if (fUsePhiCut && (phiCluster < fMinPhiCut || phiCluster > fMaxPhiCut)){delete cluster; continue;}
    if (fUseDistanceToBadChannel>0 && CheckDistanceToBadChannel(cluster,event)){delete cluster; continue;}
    //cluster quality cuts
    if (fVectorMatchedClusterIDs.size()>0 && CheckClusterForTrackMatch(cluster)){delete cluster; continue;}
    if (fUseMinEnergy && (cluster->E() < fMinEnergy)){delete cluster; continue;}
    if (fUseNCells && (cluster->GetNCells() < fMinNCells)){delete cluster; continue;}
    if (fUseNLM && (nLM < fMinNLM || nLM > fMaxNLM)){delete cluster; continue;}
    if (fUseM02 == 1 && (cluster->GetM02() < fMinM02 || cluster->GetM02() > fMaxM02)){delete cluster; continue;}
    if (fUseM02 == 2 && (cluster->GetM02() < CalculateMinM02(fMinM02CutNr, cluster->E()) || cluster->GetM02() > CalculateMaxM02(fMaxM02CutNr, cluster->E()))){delete cluster; continue;}
    if (fUseM20 && (cluster->GetM20() < fMinM20 || cluster->GetM20() > fMaxM20)){delete cluster; continue;}
    if (fUseDispersion && (cluster->GetDispersion() > fMaxDispersion)){delete cluster; continue;}
    //cluster within timing cut
    if (!(isMC>0) && (cluster->GetTOF() < fMinTimeDiff || cluster->GetTOF() > fMaxTimeDiff)){delete cluster; continue;}

    Int_t largestCellicol = -1, largestCellirow = -1;
    Int_t largestCellID = FindLargestCellInCluster(cluster,event);
    if(largestCellID==-1) AliFatal("FillHistogramsExtendedQA: FindLargestCellInCluster found cluster with NCells<1?");
    Int_t largestCelliMod = GetModuleNumberAndCellPosition(largestCellID, largestCellicol, largestCellirow);
    if(largestCelliMod < 0) AliFatal("FillHistogramsExtendedQA: GetModuleNumberAndCellPosition found SM with ID<0?");

    for(Int_t iClus2=iClus+1; iClus2<nclus; iClus2++){
      if(event->IsA()==AliESDEvent::Class()) clusterMatched = new AliESDCaloCluster(*(AliESDCaloCluster*)event->GetCaloCluster(iClus2));
      else if(event->IsA()==AliAODEvent::Class()) clusterMatched = new AliAODCaloCluster(*(AliAODCaloCluster*)event->GetCaloCluster(iClus2));

      if( (fClusterType == 1 || fClusterType == 3) && !clusterMatched->IsEMCAL()){delete clusterMatched; continue;}
      if( fClusterType == 2 && clusterMatched->GetType() !=AliVCluster::kPHOSNeutral){delete clusterMatched; continue;}

      Float_t clusPos2[3]={0,0,0};
      clusterMatched->GetPosition(clusPos2);
      TVector3 clusterMatchedVector(clusPos2[0],clusPos2[1],clusPos2[2]);
      Double_t etaclusterMatched = clusterMatchedVector.Eta();
      Double_t phiclusterMatched = clusterMatchedVector.Phi();
      if (phiclusterMatched < 0) phiclusterMatched += 2*TMath::Pi();
      Int_t nLMMatched = GetNumberOfLocalMaxima(clusterMatched, event);

      //acceptance cuts
      if (fUseEtaCut && (etaclusterMatched < fMinEtaCut || etaclusterMatched > fMaxEtaCut)){delete clusterMatched; continue;}
      if (fUsePhiCut && (phiclusterMatched < fMinPhiCut || phiclusterMatched > fMaxPhiCut)){delete clusterMatched; continue;}
      if (fUseDistanceToBadChannel>0 && CheckDistanceToBadChannel(clusterMatched,event)){delete clusterMatched; continue;}
      //cluster quality cuts
      if (fVectorMatchedClusterIDs.size()>0 && CheckClusterForTrackMatch(clusterMatched)){delete clusterMatched; continue;}
      if (fUseMinEnergy && (clusterMatched->E() < fMinEnergy)){delete clusterMatched; continue;}
      if (fUseNCells && (clusterMatched->GetNCells() < fMinNCells)){delete clusterMatched; continue;}
      if (fUseNLM && (nLMMatched < fMinNLM || nLMMatched > fMaxNLM)){delete clusterMatched; continue;}
      if (fUseM02 == 1 && (clusterMatched->GetM02() < fMinM02 || clusterMatched->GetM02() > fMaxM02)){delete clusterMatched; continue;}
      if (fUseM02 == 2 && (clusterMatched->GetM02() < CalculateMinM02(fMinM02CutNr, clusterMatched->E()) || cluster->GetM02() > CalculateMaxM02(fMaxM02CutNr, clusterMatched->E()))){delete clusterMatched; continue;}
      if (fUseM20 && (clusterMatched->GetM20() < fMinM20 || clusterMatched->GetM20() > fMaxM20)){delete clusterMatched; continue;}
      if (fUseDispersion && (clusterMatched->GetDispersion() > fMaxDispersion)){delete clusterMatched; continue;}

      // Get rowdiff and coldiff

      Int_t matched_largestCellicol = -1, matched_largestCellirow = -1;
      Int_t matched_largestCellID = FindLargestCellInCluster(clusterMatched,event);
      if(matched_largestCellID==-1) AliFatal("FillHistogramsExtendedQA: FindLargestCellInCluster found cluster with NCells<1?");
      Int_t matched_largestCelliMod = GetModuleNumberAndCellPosition(matched_largestCellID, matched_largestCellicol, matched_largestCellirow);
      if(matched_largestCelliMod < 0) AliFatal("FillHistogramsExtendedQA: GetModuleNumberAndCellPosition found SM with ID<0?");

//      cout << "\n+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
//      cout << "Cluster: " << largestCelliMod << ", " << largestCellirow << ", " << largestCellicol << " , time: " << cluster->GetTOF() << endl;
//      cout << "Matched: " << matched_largestCelliMod << ", " << matched_largestCellirow << ", " << matched_largestCellicol << " , time: " << clusterMatched->GetTOF() << endl;

      Int_t rowdiff = -100;
      Int_t coldiff = -100;
      Bool_t calculatedDiff = kFALSE;

      Int_t ClusID = largestCelliMod/2;
      Int_t matchClusID = matched_largestCelliMod/2;

      if( matched_largestCelliMod == largestCelliMod){
        rowdiff = largestCellirow - matched_largestCellirow;
        coldiff = largestCellicol - matched_largestCellicol;
        calculatedDiff = kTRUE;
      }else if( TMath::Abs(matched_largestCelliMod - largestCelliMod) == 1 && (ClusID == matchClusID) ){
        if(matched_largestCelliMod%2){
          matched_largestCelliMod -= 1;
          matched_largestCellicol += AliEMCALGeoParams::fgkEMCALCols;
        }else{
          matched_largestCelliMod += 1;
          matched_largestCellicol -= AliEMCALGeoParams::fgkEMCALCols;
        }

        if( matched_largestCelliMod == largestCelliMod ){
          rowdiff = largestCellirow - matched_largestCellirow;
          coldiff = largestCellicol - matched_largestCellicol;
          calculatedDiff = kTRUE;
        }
      }
//      cout << "\t\t ROWDIFF: " << rowdiff << endl;
//      cout << "\t\t COLDIFF: " << coldiff << endl;
//      cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n" << endl;
      //cluster outside timing cut
      if( calculatedDiff ){
        Float_t dist1D = TMath::Sqrt(TMath::Power(etaCluster-etaclusterMatched,2)+TMath::Power(phiCluster-phiclusterMatched,2));
        if( !(isMC>0) ){
          if( (clusterMatched->GetTOF() > fMinTimeDiff && clusterMatched->GetTOF() < fMaxTimeDiff) ){
            fHistClusterDistanceInTimeCut->Fill(rowdiff,coldiff);
            fHistClusterDistance1DInTimeCut->Fill(dist1D);
          }
          else fHistClusterDistanceOutTimeCut->Fill(rowdiff,coldiff);
        }else{
          fHistClusterDistanceInTimeCut->Fill(rowdiff,coldiff);
          fHistClusterDistance1DInTimeCut->Fill(dist1D);
        }
      }
      delete clusterMatched;
    }

    delete cluster;
  }

  return;
}

//________________________________________________________________________
//************** Find number of local maxima in cluster ******************
//* derived from G. Conesa Balbastre's AliCalorimeterUtils *******************
//************************************************************************
Int_t AliCaloPhotonCuts::GetNumberOfLocalMaxima(AliVCluster* cluster, AliVEvent * event){

  const Int_t   nc = cluster->GetNCells();
  
  Int_t   absCellIdList[nc];
  Float_t maxEList[nc];

  Int_t nMax = GetNumberOfLocalMaxima(cluster, event, absCellIdList, maxEList);
  
  return nMax;
}  

//________________________________________________________________________
Int_t AliCaloPhotonCuts::FindSecondLargestCellInCluster(AliVCluster* cluster, AliVEvent* event){

  const Int_t nCells      = cluster->GetNCells();
  AliVCaloCells* cells    = NULL;
  
  if (fClusterType == 1 || fClusterType == 3) 
    cells                 = event->GetEMCALCells();
  else if (fClusterType ==2 ) 
    cells                 = event->GetPHOSCells();
  
//   cout << "NCells: "<< nCells<< " cluster energy: " << cluster->E() << endl;
  Float_t eMax            = 0.;
  Int_t idMax             = -1;
  Int_t idMax2            = -1;
  Int_t iCellMax          = -1;
  
  if (nCells < 2) return idMax;
  for (Int_t iCell = 1;iCell < nCells;iCell++){
    if (cells->GetCellAmplitude(cluster->GetCellsAbsId()[iCell])> eMax){
      eMax                = cells->GetCellAmplitude(cluster->GetCellsAbsId()[iCell]);
      idMax               = cluster->GetCellsAbsId()[iCell];
      iCellMax            = iCell;
    }  
  }  
  
  eMax                    = 0.;
  for (Int_t iCell = 1;iCell < nCells;iCell++){
    if (iCell == iCellMax) continue;
    if (cells->GetCellAmplitude(cluster->GetCellsAbsId()[iCell])> eMax){
      eMax                = cells->GetCellAmplitude(cluster->GetCellsAbsId()[iCell]);
      idMax2              = cluster->GetCellsAbsId()[iCell];
    }  
  }  

  return idMax2;
}

//________________________________________________________________________
Int_t AliCaloPhotonCuts::FindLargestCellInCluster(AliVCluster* cluster, AliVEvent* event){

  const Int_t nCells      = cluster->GetNCells();
  AliVCaloCells* cells    = NULL;

  if (fClusterType == 1 || fClusterType == 3) 
    cells                 = event->GetEMCALCells();
  else if (fClusterType ==2 ) 
    cells                 = event->GetPHOSCells();
  
//   cout << "NCells: "<< nCells<< " cluster energy: " << cluster->E() << endl;
  Float_t eMax            = 0.;
  Int_t idMax             = -1;
  
  if (nCells < 1) return idMax;
  for (Int_t iCell = 0;iCell < nCells;iCell++){
    Int_t cellAbsID       = cluster->GetCellsAbsId()[iCell];
    if (cells->GetCellAmplitude(cellAbsID)> eMax){
      eMax                = cells->GetCellAmplitude(cellAbsID);
      idMax               = cellAbsID;
    }
  }
  return idMax;
  
}


//________________________________________________________________________
//************** Find number of local maxima in cluster ******************
//* derived from G. Conesa Balbastre's AliCalorimeterUtils ***************
//************************************************************************
Int_t AliCaloPhotonCuts::GetNumberOfLocalMaxima(AliVCluster* cluster, AliVEvent * event, Int_t *absCellIdList, Float_t* maxEList){

  Int_t absCellId1        = -1;
  Int_t absCellId2        = -1;
  const Int_t nCells      = cluster->GetNCells();
  AliVCaloCells* cells    = NULL;
  
  if (fClusterType == 1 || fClusterType == 3) 
    cells                 = event->GetEMCALCells();
  else if (fClusterType ==2 ) 
    cells                 = event->GetPHOSCells();
  
//   cout << "NCells: "<< nCells<< " cluster energy: " << cluster->E() << endl;
  Float_t eMax            = 0.;
  Int_t idMax             = -1;
  
  for (Int_t iCell = 0;iCell < nCells;iCell++){
    absCellIdList[iCell]  = cluster->GetCellsAbsId()[iCell];
//     Int_t imod = -1, icol = -1, irow = -1;
//     imod = GetModuleNumberAndCellPosition(absCellIdList[iCell], icol, irow);
//     cout << absCellIdList[iCell] <<"\t" << cells->GetCellAmplitude(absCellIdList[iCell]) << "\t"<< imod << "\t" << icol << "\t" << irow << endl;
    if (cells->GetCellAmplitude(absCellIdList[iCell])> eMax){
      eMax                = cells->GetCellAmplitude(absCellIdList[iCell]);
      idMax               = absCellIdList[iCell];
    }  
  }  

  // find the largest separated cells in cluster
  for (Int_t iCell = 0;iCell < nCells;iCell++){
    // check whether valid cell number is selected
    if (absCellIdList[iCell] >= 0){
      // store current energy and cell id
      absCellId1          = cluster->GetCellsAbsId()[iCell];
      Float_t en1         = cells->GetCellAmplitude(absCellId1);
      
      // loop over other cells in cluster
      for (Int_t iCellN = 0;iCellN < nCells;iCellN++){
        // jump out if array has changed in the meantime
        if (absCellIdList[iCell] == -1) continue;
        // get cell id & check whether its valid
        absCellId2        = cluster->GetCellsAbsId()[iCellN];

        // don't compare to yourself
        if (absCellId2 == -1) continue;
        if (absCellId1 == absCellId2) continue;
        
        // get cell energy of second cell
        Float_t en2       = cells->GetCellAmplitude(absCellId2);
        
        // check if cells are Neighbours
        if (AreNeighbours(absCellId1, absCellId2)){
          // determine which cell has larger energy, mask the other
//           cout << "found neighbour: " << absCellId1 << "\t" << absCellId2 << endl;
//           cout << "energies: " << en1 << "\t" << en2 << endl;
          if (en1 > en2 ){
            absCellIdList[iCellN]       = -1;
            if (en1 < en2 + fLocMaxCutEDiff)
                absCellIdList[iCell]    = -1;
          } else {
            absCellIdList[iCell]        = -1;
            if (en1 > en2 - fLocMaxCutEDiff)
                absCellIdList[iCellN]   = -1;
          }  
        }
      }
    }    
  }  

  // shrink list of cells to only maxima
  Int_t nMaximaNew        = 0;
  for (Int_t iCell = 0;iCell < nCells;iCell++){
//     cout << iCell << "\t" << absCellIdList[iCell] << endl;
    if (absCellIdList[iCell] > -1){
      Float_t en          = cells->GetCellAmplitude(absCellIdList[iCell]);
      // check whether cell energy is larger than required seed
      if (en < fSeedEnergy) continue;
      absCellIdList[nMaximaNew]   = absCellIdList[iCell];
      maxEList[nMaximaNew]        = en;
      nMaximaNew++;
    }  
  }  

  // check whether a local maximum was found
  // if no maximum was found use highest cell as maximum
  if (nMaximaNew == 0){
    nMaximaNew            = 1;
    maxEList[0]           = eMax;
    absCellIdList[0]      = idMax;
  }  

  return nMaximaNew;
}  

//________________________________________________________________________
//************** Function to determine neighbours of cells ***************
//* derived from G. Conesa Balbastre's AliCalorimeterUtils ***************
//************************************************************************
Bool_t AliCaloPhotonCuts::AreNeighbours(Int_t absCellId1, Int_t absCellId2){
  Bool_t areNeighbours = kFALSE ;
  
  Int_t irow1 = -1, icol1 = -1;
  Int_t irow2 = -1, icol2 = -1;
  
  Int_t rowdiff =  0, coldiff =  0;
  
  Int_t nSupMod1          = GetModuleNumberAndCellPosition(absCellId1, icol1, irow1);
  Int_t nSupMod2          = GetModuleNumberAndCellPosition(absCellId2, icol2, irow2);
    
  // check if super modules are correct
  if (nSupMod1== -1 || nSupMod2 == -1) return areNeighbours;

  if(fClusterType==1 && nSupMod1!=nSupMod2) {
    // In case of a shared cluster, index of SM in C side, columns start at 48 and ends at 48*2-1
    // C Side impair SM, nSupMod%2=1;A side pair SM nSupMod%2=0
    if(nSupMod1%2) icol1+=AliEMCALGeoParams::fgkEMCALCols;
    else           icol2+=AliEMCALGeoParams::fgkEMCALCols;
  }
  
  rowdiff = TMath::Abs( irow1 - irow2 ) ;
  coldiff = TMath::Abs( icol1 - icol2 ) ;
  
//   if (( coldiff <= 1 )  && ( rowdiff <= 1 ) && (coldiff + rowdiff <= 2))
  if ((coldiff + rowdiff == 1 ))
    areNeighbours         = kTRUE ;
  
  return areNeighbours;
}


//________________________________________________________________________
//************** Function to obtain module number, row and column ********
//* derived from G. Conesa Balbastre's AliCalorimeterUtils ***************
//************************************************************************
Int_t AliCaloPhotonCuts::GetModuleNumberAndCellPosition(Int_t absCellId, Int_t & icol, Int_t & irow){
  if( fClusterType == 1 || fClusterType == 3){ //EMCAL & DCAL
    fGeomEMCAL = AliEMCALGeometry::GetInstance();
    if(!fGeomEMCAL) AliFatal("EMCal geometry not initialized!");
  } else if( fClusterType == 2 ){ //PHOS
    fGeomPHOS = AliPHOSGeometry::GetInstance();
    if(!fGeomPHOS) AliFatal("PHOS geometry not initialized!");
  }
  
  Int_t imod = -1;Int_t iTower = -1, iIphi = -1, iIeta = -1;
  if( fClusterType == 1 || fClusterType == 3){
    fGeomEMCAL->GetCellIndex(absCellId,imod,iTower,iIphi,iIeta);
    fGeomEMCAL->GetCellPhiEtaIndexInSModule(imod,iTower,iIphi,iIeta,irow,icol);
  } else if ( fClusterType == 2 ){
    Int_t relId[4];
    fGeomPHOS->AbsToRelNumbering(absCellId,relId);
    irow                  = relId[2];
    icol                  = relId[3];
    imod                  = relId[0]-1;
  }
  return imod;
}

//___________________________________________________________________________
// Split energy of cluster between the 2 local maxima, sum energy on 3x3, and if the 2 
// maxima are too close and have common cells, split the energy between the 2.
//* derived from G. Conesa Balbastre's AliCalorimeterUtils *******************
//___________________________________________________________________________
void AliCaloPhotonCuts::SplitEnergy(Int_t absCellId1, Int_t absCellId2,
                                    AliVCluster* cluster,
                                    AliVEvent* event, 
                                    Int_t isMC,
                                    AliAODCaloCluster* cluster1,
                                    AliAODCaloCluster* cluster2){
                    
  const Int_t ncells      = cluster->GetNCells();
  Int_t absCellIdList[ncells];

  AliVCaloCells* cells    = NULL;
  if (fClusterType == 1 || fClusterType == 3) 
    cells                 = event->GetEMCALCells();
  else if (fClusterType ==2 ) 
    cells                 = event->GetPHOSCells();

  Float_t e1              = 0;
  Float_t e2              = 0;
  Float_t eCluster        = 0;
  
  for(Int_t iCell    = 0;iCell < ncells;iCell++ ) {
    absCellIdList[iCell]  = cluster->GetCellsAbsId()[iCell];
    Float_t ec            = cells->GetCellAmplitude(absCellIdList[iCell]);
    eCluster+=ec;
  }

  UShort_t absCellIdList1 [12];
  Double_t fracList1      [12];
  UShort_t absCellIdList2 [12];
  Double_t fracList2      [12];

  // Init counters and variables
  Int_t ncells1         = 1 ;
  absCellIdList1[0]     = absCellId1 ;
  fracList1 [0]         = 1. ;
  
  Float_t ecell1        = cells->GetCellAmplitude(absCellId1);
  e1                    = ecell1;
  
  Int_t ncells2         = 1 ;
  absCellIdList2[0]     = absCellId2 ;
  fracList2 [0]         = 1. ;
  
  Float_t ecell2        = cells->GetCellAmplitude(absCellId2);
  e2                    = ecell2;
    
//   cout << "Cluster: " << eCluster << "\t cell1: " << absCellId1 << "\t" << e1 << "\t cell2: " << absCellId2 << "\t" << e2 << endl;
  // Very rough way to share the cluster energy
  Float_t eRemain           = (eCluster-ecell1-ecell2)/2;
  Float_t shareFraction1    = (ecell1+eRemain)/eCluster;
  Float_t shareFraction2    = (ecell2+eRemain)/eCluster;

//   cout << eRemain << "\t" << shareFraction1<< "\t" << shareFraction2 << endl;
  
  for(Int_t iCell = 0;iCell < ncells;iCell++){
    
    Int_t absId         = absCellIdList[iCell];
    if ( absId==absCellId1 || absId==absCellId2 || absId < 0 ) continue;
    
    Float_t ecell = cells->GetCellAmplitude(absId);
    if(AreNeighbours(absCellId1,absId )){ 
      absCellIdList1[ncells1] = absId;
      if(AreNeighbours(absCellId2,absId )){ 
        fracList1[ncells1] = shareFraction1;
        e1 += ecell*shareFraction1;
      } else {
        fracList1[ncells1] = 1.;
        e1 += ecell;
      }    
      ncells1++;
    } // neigbour to cell1
    
    if(AreNeighbours(absCellId2,absId )) { 
      absCellIdList2[ncells2]= absId;
    
      if(AreNeighbours(absCellId1,absId )){ 
        fracList2[ncells2] = shareFraction2;
        e2 += ecell*shareFraction2;
      } else { 
        fracList2[ncells2] = 1.;
        e2 += ecell;
      }
      ncells2++;
    } // neigbour to cell2  
  }
//   cout << "Cluster: " << eCluster << "\t cell1: " << absCellId1 << "\t" << e1 << "\t cell2: " << absCellId2 << "\t" << e2 << endl;
        
  cluster1->SetE(e1);
  cluster2->SetE(e2);
  
  cluster1->SetNCells(ncells1);
  cluster2->SetNCells(ncells2);
  
  cluster1->SetCellsAbsId(absCellIdList1);
  cluster2->SetCellsAbsId(absCellIdList2);
  
  cluster1->SetCellsAmplitudeFraction(fracList1);
  cluster2->SetCellsAmplitudeFraction(fracList2);
  
  // Correct linearity
  if (fClusterType == 1 || fClusterType == 3){
    CorrectEMCalNonLinearity(cluster1, isMC) ;
    CorrectEMCalNonLinearity(cluster2, isMC) ;
  }

  // Initialize EMCAL rec utils if not initialized
  if(!fEMCALInitialized && (fClusterType == 1 || fClusterType == 3) ) InitializeEMCAL(event);
  
  if(fEMCALInitialized && (fClusterType == 1 || fClusterType == 3) ){
    fEMCALRecUtils->RecalculateClusterPosition(fGeomEMCAL, cells, cluster1);
    fEMCALRecUtils->RecalculateClusterPosition(fGeomEMCAL, cells, cluster2);
  }
}

//________________________________________________________________________
Bool_t AliCaloPhotonCuts::CheckDistanceToBadChannel(AliVCluster* cluster, AliVEvent* event)
{
  if(fUseDistanceToBadChannel != 1 && fUseDistanceToBadChannel != 2) return kFALSE;

  //not yet fully implemented for PHOS:
  if( fClusterType == 2 ) return kFALSE;

  if( (fClusterType == 1 || fClusterType == 3) && !fEMCALInitialized ) InitializeEMCAL(event);
  if( fClusterType == 2 && ( !fPHOSInitialized || (fPHOSCurrentRun != event->GetRunNumber()) ) ) InitializePHOS(event);

  Int_t largestCellID = FindLargestCellInCluster(cluster,event);
  if(largestCellID==-1) AliFatal("CheckDistanceToBadChannel: FindLargestCellInCluster found cluster with NCells<1?");

  Int_t largestCellicol = -1, largestCellirow = -1;
  Int_t rowdiff =  0, coldiff =  0;

  Int_t largestCelliMod = GetModuleNumberAndCellPosition(largestCellID, largestCellicol, largestCellirow);
  if(largestCelliMod < 0) AliFatal("CheckDistanceToBadChannel: GetModuleNumberAndCellPosition found SM with ID<0?");

  Int_t nMinRows = 0, nMaxRows = 0;
  Int_t nMinCols = 0, nMaxCols = 0;

  Bool_t checkNextSM = kFALSE;
  Int_t distanceForLoop = fMinDistanceToBadChannel+1;
  if( fClusterType == 1 ){
    nMinRows = largestCellirow - distanceForLoop;
    nMaxRows = largestCellirow + distanceForLoop;
    if(nMinRows < 0) nMinRows = 0;
    if(nMaxRows > AliEMCALGeoParams::fgkEMCALRows) nMaxRows = AliEMCALGeoParams::fgkEMCALRows;

    nMinCols = largestCellicol - distanceForLoop;
    nMaxCols = largestCellicol + distanceForLoop;

    if(largestCelliMod%2){
      if(nMinCols < 0){
        nMinCols = 0;
        checkNextSM = kTRUE;
      }
      if(nMaxCols > AliEMCALGeoParams::fgkEMCALCols) nMaxCols = AliEMCALGeoParams::fgkEMCALCols;
    }else{
      if(nMinCols < 0) nMinCols = 0;
      if(nMaxCols > AliEMCALGeoParams::fgkEMCALCols){
        nMaxCols = AliEMCALGeoParams::fgkEMCALCols;
        checkNextSM = kTRUE;
      }
    }
  }else if( fClusterType == 3 ){
    nMinRows = largestCellirow - distanceForLoop;
    nMaxRows = largestCellirow + distanceForLoop;
    if(nMinRows < 0) nMinRows = 0;
    if(nMaxRows > AliEMCALGeoParams::fgkEMCALCols) nMaxRows = AliEMCALGeoParams::fgkEMCALCols; //AliEMCALGeoParams::fgkDCALRows; <- doesnt exist yet (DCAl = EMCAL here)

    nMinCols = largestCellicol - distanceForLoop;
    nMaxCols = largestCellicol + distanceForLoop;
    if(nMinCols < 0) nMinCols = 0;
    if(nMaxCols > fgkDCALCols) nMaxCols = fgkDCALCols; // AliEMCALGeoParams::fgkDCALCols; <- doesnt exist yet

  }else if( fClusterType == 2 ){
//    nMaxRows = 64;
//    nMaxCols = 56;
  }

//  cout << "Cluster: " << fClusterType << ",checkNextSM: " << checkNextSM << endl;
//  cout << "largestCell: " << largestCellID << ",mod: " << largestCelliMod << ",col: " << largestCellicol << ",row: " << largestCellirow << endl;
//  cout << "distanceForLoop: " << distanceForLoop << ",nMinRows: " << nMinRows << ",nMaxRows: " << nMaxRows << ",nMinCols: " << nMinCols << ",nMaxCols: " << nMaxCols << endl;

  //check bad cells within respective SM
  for (Int_t irow = nMinRows;irow < nMaxRows;irow++)
  {
    for (Int_t icol = nMinCols;icol < nMaxCols;icol++)
    {
      if(irow == largestCellirow && icol == largestCellicol) continue;

      Int_t iBadCell = 0;
      if( (fClusterType == 1 || fClusterType == 3) && largestCelliMod<fEMCALBadChannelsMap->GetEntries()){
        iBadCell = (Int_t) ((TH2I*)fEMCALBadChannelsMap->At(largestCelliMod))->GetBinContent(icol,irow);
      }else if( fClusterType == 2 && fPHOSBadChannelsMap[largestCelliMod+1]){
        iBadCell = (Int_t) ((TH2I*)fPHOSBadChannelsMap[largestCelliMod+1])->GetBinContent(icol,irow);
      }
      //cout << "largestCelliMod: " << largestCelliMod << ",iBadCell: " << iBadCell << ",icol: " << icol << ",irow: " << irow << endl;
      if(iBadCell==0) continue;

      rowdiff = TMath::Abs( largestCellirow - irow ) ;
      coldiff = TMath::Abs( largestCellicol - icol ) ;
      //cout << "rowdiff: " << rowdiff << ",coldiff: " << coldiff << endl;
      if(fUseDistanceToBadChannel==1){
        if ((coldiff + rowdiff <= fMinDistanceToBadChannel )) return kTRUE;
      }else if(fUseDistanceToBadChannel==2){
        if (( coldiff <= fMinDistanceToBadChannel )  && ( rowdiff <= fMinDistanceToBadChannel ) && (coldiff + rowdiff <= fMinDistanceToBadChannel*2)) return kTRUE;
      }
      //cout << "not within distanceToBadChannel!" << endl;
    }
  }

  //check bad cells in neighboring SM only if within chosen distanceToBadChannel from maxEnergyCell the next SM could be reached
  if(checkNextSM) {
    // In case of a shared cluster, index of SM in C side, columns start at 48 and ends at 48*2-1
    // C Side impair SM, nSupMod%2=1;A side pair SM nSupMod%2=0
    if( fClusterType == 1 ){
      if(largestCelliMod%2){
        nMinCols = largestCellicol - distanceForLoop + AliEMCALGeoParams::fgkEMCALCols;
        nMaxCols = AliEMCALGeoParams::fgkEMCALCols;

        largestCelliMod -= 1;
        largestCellicol += AliEMCALGeoParams::fgkEMCALCols;
      }else{
        nMinCols = 0;
        nMaxCols = largestCellicol + distanceForLoop - AliEMCALGeoParams::fgkEMCALCols;

        largestCelliMod += 1;
        largestCellicol -= AliEMCALGeoParams::fgkEMCALCols;
      }
    }else if( fClusterType == 2 ){
//      nMaxRows = 64;
//      nMaxCols = 56;
    }
    //cout << "largestCell: " << largestCellID << ",mod: " << largestCelliMod << ",col: " << largestCellicol << ",row: " << largestCellirow << endl;
    //cout << "distanceForLoop: " << distanceForLoop << ",nMinRows: " << nMinRows << ",nMaxRows: " << nMaxRows << ",nMinCols: " << nMinCols << ",nMaxCols: " << nMaxCols << endl;
    for (Int_t irow = nMinRows;irow < nMaxRows;irow++)
    {
      for (Int_t icol = nMinCols;icol < nMaxCols;icol++)
      {
        Int_t iBadCell = 0;
        if( fClusterType == 1 && largestCelliMod<fEMCALBadChannelsMap->GetEntries()){
          iBadCell = (Int_t) ((TH2I*)fEMCALBadChannelsMap->At(largestCelliMod))->GetBinContent(icol,irow);
        }else if( fClusterType == 2 && fPHOSBadChannelsMap[largestCelliMod+1]){
          iBadCell = (Int_t) ((TH2I*)fPHOSBadChannelsMap[largestCelliMod+1])->GetBinContent(icol,irow);
        }
        //cout << "largestCelliMod: " << largestCelliMod << ",iBadCell: " << iBadCell << ",icol: " << icol << ",irow: " << irow << endl;
        if(iBadCell==0) continue;

        rowdiff = TMath::Abs( largestCellirow - irow ) ;
        coldiff = TMath::Abs( largestCellicol - icol ) ;
        //cout << "rowdiff: " << rowdiff << ",coldiff: " << coldiff << endl;
        if(fUseDistanceToBadChannel==1){
          if ((coldiff + rowdiff <= fMinDistanceToBadChannel )) return kTRUE;
        }else if(fUseDistanceToBadChannel==2){
          if (( coldiff <= fMinDistanceToBadChannel )  && ( rowdiff <= fMinDistanceToBadChannel ) && (coldiff + rowdiff <= fMinDistanceToBadChannel*2)) return kTRUE;
        }
        //cout << "not within distanceToBadChannel!" << endl;
      }
    }
  }

  return kFALSE;
}


//________________________________________________________________________
Bool_t AliCaloPhotonCuts::ClusterIsSelected(AliVCluster *cluster, AliVEvent * event, AliVEvent * mcEvent, Int_t isMC, Double_t weight, Long_t clusterID)
{
  //Selection of Reconstructed photon clusters with Calorimeters

  FillClusterCutIndex(kPhotonIn);

  // do NonLinearity if switched on
  if(fUseNonLinearity && cluster->IsEMCAL()){
    if(fHistEnergyOfClusterBeforeNL) fHistEnergyOfClusterBeforeNL->Fill(cluster->E(),weight);
    CorrectEMCalNonLinearity(cluster,isMC);
    if(fHistEnergyOfClusterAfterNL) fHistEnergyOfClusterAfterNL->Fill(cluster->E(),weight);
  }

//  Double_t vertex[3] = {0,0,0};
//  event->GetPrimaryVertex()->GetXYZ(vertex);
    // TLorentzvector with cluster
//  TLorentzVector clusterVector;
//  cluster->GetMomentum(clusterVector,vertex);

  Float_t clusPos[3]={0,0,0};
  cluster->GetPosition(clusPos);
  TVector3 clusterVector(clusPos[0],clusPos[1],clusPos[2]);
  Double_t etaCluster = clusterVector.Eta();
  Double_t phiCluster = clusterVector.Phi();
  if (phiCluster < 0) phiCluster += 2*TMath::Pi();

  // Histos before cuts
  if(fHistClusterEtavsPhiBeforeAcc) fHistClusterEtavsPhiBeforeAcc->Fill(phiCluster,etaCluster,weight);
  
  // Cluster Selection - 0= accept any calo cluster
  if (fClusterType > 0){
    //Select EMCAL cluster
    if ( (fClusterType == 1 || fClusterType == 3) && !cluster->IsEMCAL()){
      FillClusterCutIndex(kDetector);
      return kFALSE;
    }
    //Select PHOS cluster
    if (fClusterType == 2 && !cluster->IsPHOS()){
      FillClusterCutIndex(kDetector);
      // Apply nonlinearity on PHOS cluster
      if(fHistEnergyOfClusterBeforeNL) fHistEnergyOfClusterBeforeNL->Fill(cluster->E(),weight);
      CorrectEMCalNonLinearity(cluster,isMC);
      if(fHistEnergyOfClusterAfterNL) fHistEnergyOfClusterAfterNL->Fill(cluster->E(),weight);
      return kFALSE;
    }
  }
  
  // Acceptance Cuts
  if(!AcceptanceCuts(cluster,event,weight)){
    FillClusterCutIndex(kAcceptance);
    return kFALSE;
  }
  // Cluster Quality Cuts
  if(!ClusterQualityCuts(cluster,event,mcEvent,isMC,weight,clusterID)){
    FillClusterCutIndex(kClusterQuality);
    return kFALSE;
  }

  // Photon passed cuts
  FillClusterCutIndex(kPhotonOut);
  return kTRUE;
}


//________________________________________________________________________
Bool_t AliCaloPhotonCuts::AcceptanceCuts(AliVCluster *cluster, AliVEvent* event, Double_t weight)
{
   // Exclude certain areas for photon reconstruction

  Int_t cutIndex=0;
  if(fHistAcceptanceCuts)fHistAcceptanceCuts->Fill(cutIndex);
  cutIndex++;

  Float_t clusPos[3]={0,0,0};
  cluster->GetPosition(clusPos);
  TVector3 clusterVector(clusPos[0],clusPos[1],clusPos[2]);
  Double_t etaCluster = clusterVector.Eta();
  Double_t phiCluster = clusterVector.Phi();
  if (phiCluster < 0) phiCluster += 2*TMath::Pi();
  
  // check eta range
  if (fUseEtaCut){
    if (etaCluster < fMinEtaCut || etaCluster > fMaxEtaCut){
      if(fHistAcceptanceCuts)fHistAcceptanceCuts->Fill(cutIndex);
      return kFALSE;
    }
  }
  cutIndex++;
  
  // check phi range
  if (fUsePhiCut ){
    if (phiCluster < fMinPhiCut || phiCluster > fMaxPhiCut){
      if(fHistAcceptanceCuts)fHistAcceptanceCuts->Fill(cutIndex);
      return kFALSE;
    }
  }
  cutIndex++;
  
  // check distance to bad channel
  if (fUseDistanceToBadChannel>0){
    if (CheckDistanceToBadChannel(cluster,event)){
      if(fHistAcceptanceCuts)fHistAcceptanceCuts->Fill(cutIndex);
      return kFALSE;
    }  
  }
  //alternatively check histogram fHistoModifyAcc if cluster should be rejected
  if(fHistoModifyAcc){
    if(fHistoModifyAcc->GetBinContent(FindLargestCellInCluster(cluster,event)) < 1){
      if(fHistAcceptanceCuts)fHistAcceptanceCuts->Fill(cutIndex);
      return kFALSE;
    }
  }
  cutIndex++;
  if(fHistAcceptanceCuts)fHistAcceptanceCuts->Fill(cutIndex);

  // Histos after cuts
  if(fHistClusterEtavsPhiAfterAcc) fHistClusterEtavsPhiAfterAcc->Fill(phiCluster,etaCluster,weight);
  
  return kTRUE;
}

//________________________________________________________________________
Bool_t AliCaloPhotonCuts::MatchConvPhotonToCluster(AliAODConversionPhoton* convPhoton, AliVCluster* cluster, AliVEvent* event, Double_t weight){

  if (!fUseDistTrackToCluster) return kFALSE;
  if( (fClusterType == 1 || fClusterType == 3) && !fEMCALInitialized ) InitializeEMCAL(event);
  if( fClusterType == 2 && ( !fPHOSInitialized || (fPHOSCurrentRun != event->GetRunNumber()) ) ) InitializePHOS(event);
  
  AliESDEvent *esdev = dynamic_cast<AliESDEvent*>(event);
  AliAODEvent *aodev = 0;
  if (!esdev) {
    aodev = dynamic_cast<AliAODEvent*>(event);
    if (!aodev) {
      AliError("Task needs AOD or ESD event, returning");
      return kFALSE;
    }
  }

  if(!cluster->IsEMCAL() && !cluster->IsPHOS()){AliError("Cluster is neither EMCAL nor PHOS, returning");return kFALSE;}

  Float_t clusterPosition[3] = {0,0,0};
  cluster->GetPosition(clusterPosition);
  Double_t clusterR = TMath::Sqrt( clusterPosition[0]*clusterPosition[0] + clusterPosition[1]*clusterPosition[1] );
  if(fHistClusterRBeforeQA) fHistClusterRBeforeQA->Fill(clusterR,weight);

//cout << "+++++++++ Cluster: x, y, z, R" << clusterPosition[0] << ", " << clusterPosition[1] << ", " << clusterPosition[2] << ", " << clusterR << "+++++++++" << endl;

  Bool_t matched = kFALSE;
  for (Int_t i = 0;i < 2;i++){
    Int_t tracklabel = convPhoton->GetLabel(i);
    AliVTrack *inTrack = 0x0;
    if(esdev) {
      if(tracklabel > event->GetNumberOfTracks() ) continue;
      inTrack = esdev->GetTrack(tracklabel);
    } else {
      if(((AliV0ReaderV1*)AliAnalysisManager::GetAnalysisManager()->GetTask(fV0ReaderName.Data()))->AreAODsRelabeled()){
        inTrack = dynamic_cast<AliVTrack*>(event->GetTrack(tracklabel));
      } else {
        for(Int_t ii=0;ii<event->GetNumberOfTracks();ii++) {
          inTrack = dynamic_cast<AliVTrack*>(event->GetTrack(ii));
          if(inTrack){
            if(inTrack->GetID() == tracklabel) {
              break;
            }
          }
        }
      }
    }

    Float_t dEta = 0;
    Float_t dPhi = 0;
    Bool_t propagated = fCaloTrackMatcher->PropagateV0TrackToClusterAndGetMatchingResidual(inTrack,cluster,event,dEta,dPhi);
    if (propagated){
      Float_t dR2 = dPhi*dPhi + dEta*dEta;
      if(fHistDistanceTrackToClusterBeforeQA)fHistDistanceTrackToClusterBeforeQA->Fill(TMath::Sqrt(dR2), weight);
      if(fHistClusterdEtadPhiBeforeQA) fHistClusterdEtadPhiBeforeQA->Fill(dEta, dPhi, weight);

      Float_t clusM02 = (Float_t) cluster->GetM02();
      Float_t clusM20 = (Float_t) cluster->GetM20();
      if(!fDoLightOutput && (fExtendedMatchAndQA == 1 || fExtendedMatchAndQA == 3 || fExtendedMatchAndQA == 5 )){
        if(inTrack->Charge() > 0) {
          fHistClusterdEtadPhiPosTracksBeforeQA->Fill(dEta, dPhi, weight);
          fHistClusterdPhidPtPosTracksBeforeQA->Fill(dPhi, inTrack->Pt(), weight);
          if(inTrack->P() < 0.75) fHistClusterdEtadPhiPosTracksP_000_075BeforeQA->Fill(dEta, dPhi, weight);
          else if(inTrack->P() < 1.25) fHistClusterdEtadPhiPosTracksP_075_125BeforeQA->Fill(dEta, dPhi, weight);
          else fHistClusterdEtadPhiPosTracksP_125_999BeforeQA->Fill(dEta, dPhi, weight);
        } else {
          fHistClusterdEtadPhiNegTracksBeforeQA->Fill(dEta, dPhi, weight);
          fHistClusterdPhidPtNegTracksBeforeQA->Fill(dPhi, inTrack->Pt(), weight);
          if(inTrack->P() < 0.75) fHistClusterdEtadPhiNegTracksP_000_075BeforeQA->Fill(dEta, dPhi, weight);
          else if(inTrack->P() < 1.25) fHistClusterdEtadPhiNegTracksP_075_125BeforeQA->Fill(dEta, dPhi, weight);
          else fHistClusterdEtadPhiNegTracksP_125_999BeforeQA->Fill(dEta, dPhi, weight);
        }
        fHistClusterdEtadPtBeforeQA->Fill(dEta, inTrack->Pt(), weight);
        fHistClusterM20M02BeforeQA->Fill(clusM20, clusM02, weight);
        if(fCurrentMC != kNoMC && fIsMC > 0){
          Int_t clusterMCLabel = cluster->GetLabel();
          Int_t convPhotonDaughterLabel = -1;
          if(inTrack->Charge() > 0) convPhotonDaughterLabel = convPhoton->GetMCLabelPositive();
          else convPhotonDaughterLabel = convPhoton->GetMCLabelNegative();
          if( (convPhotonDaughterLabel != -1) && (clusterMCLabel != -1) && (convPhotonDaughterLabel == clusterMCLabel)){ //real match
            fHistClusterdEtadPtTrueMatched->Fill(dEta, inTrack->Pt(), weight);
            if(inTrack->Charge() > 0) fHistClusterdPhidPtPosTracksTrueMatched->Fill(dPhi, inTrack->Pt(), weight);
            else fHistClusterdPhidPtNegTracksTrueMatched->Fill(dPhi, inTrack->Pt(), weight);
          }
        }
      }

      Bool_t match_dEta = (TMath::Abs(dEta) < fMaxDistTrackToClusterEta) ? kTRUE : kFALSE;
      Bool_t match_dPhi = kFALSE;
      if( (inTrack->Charge() > 0) && (dPhi > fMinDistTrackToClusterPhi) && (dPhi < fMaxDistTrackToClusterPhi) ) match_dPhi = kTRUE;
      else if( (inTrack->Charge() < 0) && (dPhi < -fMinDistTrackToClusterPhi) && (dPhi > -fMaxDistTrackToClusterPhi) ) match_dPhi = kTRUE;

      if(fUsePtDepTrackToCluster){
        if( TMath::Abs(dEta) < fFuncPtDepEta->Eval(inTrack->Pt())) match_dEta = kTRUE;
        else match_dEta = kFALSE;

        if( TMath::Abs(dPhi) < fFuncPtDepPhi->Eval(inTrack->Pt())) match_dPhi = kTRUE;
        else match_dPhi = kFALSE;
      }

      if(match_dEta && match_dPhi){
            //if(dR2 < fMinDistTrackToCluster*fMinDistTrackToCluster){
        matched = kTRUE;
        if(fHistClusterdEtadPtAfterQA) fHistClusterdEtadPtAfterQA->Fill(dEta,inTrack->Pt());
        if(fHistClusterdPhidPtAfterQA) fHistClusterdPhidPtAfterQA->Fill(dPhi,inTrack->Pt());
      } else {
        if(fHistDistanceTrackToClusterAfterQA)fHistDistanceTrackToClusterAfterQA->Fill(TMath::Sqrt(dR2), weight);
        if(fHistClusterdEtadPhiAfterQA) fHistClusterdEtadPhiAfterQA->Fill(dEta, dPhi, weight);
        if(fHistClusterRAfterQA) fHistClusterRAfterQA->Fill(clusterR, weight);
        if(!fDoLightOutput && (fExtendedMatchAndQA == 1 || fExtendedMatchAndQA == 3 || fExtendedMatchAndQA == 5 )){
          if(inTrack->Charge() > 0) fHistClusterdEtadPhiPosTracksAfterQA->Fill(dEta, dPhi, weight);
          else fHistClusterdEtadPhiNegTracksAfterQA->Fill(dEta, dPhi, weight);
          fHistClusterM20M02AfterQA->Fill(clusM20, clusM02, weight);
        }
      }  
    }
  }

  return matched;

}

//________________________________________________________________________
void AliCaloPhotonCuts::MatchTracksToClusters(AliVEvent* event, Double_t weight, Bool_t isEMCalOnly){
  if( !fUseDistTrackToCluster ) return;
  if( (fClusterType == 1 || fClusterType == 3) && !fEMCALInitialized ) InitializeEMCAL(event);
  if( fClusterType == 2 && ( !fPHOSInitialized || (fPHOSCurrentRun != event->GetRunNumber()) ) ) InitializePHOS(event);

  // not yet fully implemented + tested for PHOS
  // if( fClusterType != 1 && fClusterType != 3) return;

  fVectorMatchedClusterIDs.clear();

  Int_t nClus = event->GetNumberOfCaloClusters();
  Int_t nModules = 0;

  if(fClusterType == 1 || fClusterType == 3){
    fGeomEMCAL = AliEMCALGeometry::GetInstance();
    if(!fGeomEMCAL){ AliFatal("EMCal geometry not initialized!");}
    nModules = fGeomEMCAL->GetNumberOfSuperModules();
  }else if(fClusterType == 2){
    fGeomPHOS = AliPHOSGeometry::GetInstance();
    if(!fGeomPHOS) AliFatal("PHOS geometry not initialized!");
    nModules = fGeomPHOS->GetNModules();
  }

  AliESDEvent *esdev = dynamic_cast<AliESDEvent*>(event);
  AliAODEvent *aodev = 0;
  if (!esdev) {
    aodev = dynamic_cast<AliAODEvent*>(event);
    if (!aodev) {
      AliError("Task needs AOD or ESD event, returning");
      return;
    }
  }

  // if not EMCal only reconstruction (-> hybrid PCM+EMCal), use only primary tracks for basic track matching procedure
  AliESDtrackCuts *EsdTrackCuts = 0x0;
  if(!isEMCalOnly && esdev){
    // Using standard function for setting Cuts
    Int_t runNumber = event->GetRunNumber();
    // if LHC11a or earlier or if LHC13g or if LHC12a-i -> use 2010 cuts
    if( (runNumber<=146860) || (runNumber>=197470 && runNumber<=197692) || (runNumber>=172440 && runNumber<=193766) ){
      EsdTrackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2010();
    // else if run2 data use 2015 PbPb cuts
    }else if (runNumber>=209122){
      EsdTrackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2015PbPb();
    // else use 2011 version of track cuts
    }else{
      EsdTrackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011();
    }
    EsdTrackCuts->SetMaxDCAToVertexZ(2);
    EsdTrackCuts->SetEtaRange(-0.8, 0.8);
    EsdTrackCuts->SetPtRange(0.15);
  }

//  cout << "MatchTracksToClusters: " << event->GetNumberOfTracks() << ", " << fIsPureCalo << ", " << fUseDistTrackToCluster << endl;

  for (Int_t itr=0;itr<event->GetNumberOfTracks();itr++){
    AliExternalTrackParam *trackParam = 0;
    AliVTrack *inTrack = 0x0;
    if(esdev){
      inTrack = esdev->GetTrack(itr);
      if(!inTrack) continue;
      AliESDtrack *esdt = dynamic_cast<AliESDtrack*>(inTrack);
      if(!isEMCalOnly){ //match only primaries for hybrid reconstruction schemes
         if(!EsdTrackCuts->AcceptTrack(esdt)) continue;
      }

      const AliExternalTrackParam *in = esdt->GetInnerParam();
      if (!in){AliDebug(2, "Could not get InnerParam of Track, continue");continue;}
      trackParam = new AliExternalTrackParam(*in);
    } else if(aodev) {
      inTrack = dynamic_cast<AliVTrack*>(aodev->GetTrack(itr));
      if(!inTrack) continue;
      AliAODTrack *aodt = dynamic_cast<AliAODTrack*>(inTrack);
      if(!isEMCalOnly){ //match only primaries for hybrid reconstruction schemes
        if(!aodt->IsHybridGlobalConstrainedGlobal()) continue;
        if(TMath::Abs(aodt->Eta())>0.8) continue;
        if(aodt->Pt()<0.15) continue;
      }

    }

    Float_t clsPos[3] = {0.,0.,0.};
    for(Int_t iclus=0;iclus < nClus;iclus++){
      AliVCluster* cluster = event->GetCaloCluster(iclus);
      if (!cluster) continue;
      Float_t dEta, dPhi;
      if(!fCaloTrackMatcher->GetTrackClusterMatchingResidual(inTrack->GetID(),cluster->GetID(),dEta,dPhi)) continue;
      cluster->GetPosition(clsPos);
      Float_t clusterR = TMath::Sqrt( clsPos[0]*clsPos[0] + clsPos[1]*clsPos[1] );
      Float_t dR2 = dPhi*dPhi + dEta*dEta;
//      cout << "dEta/dPhi: " << dEta << ", " << dPhi << " - ";
//      cout << dR2 << endl;
      if(isEMCalOnly && fHistDistanceTrackToClusterBeforeQA)fHistDistanceTrackToClusterBeforeQA->Fill(TMath::Sqrt(dR2), weight);
      if(isEMCalOnly && fHistClusterdEtadPhiBeforeQA) fHistClusterdEtadPhiBeforeQA->Fill(dEta, dPhi, weight);
      if(isEMCalOnly && fHistClusterRBeforeQA) fHistClusterRBeforeQA->Fill(clusterR, weight);

      Float_t clusM02 = (Float_t) cluster->GetM02();
      Float_t clusM20 = (Float_t) cluster->GetM20();
      if(isEMCalOnly && !fDoLightOutput && (fExtendedMatchAndQA == 1 || fExtendedMatchAndQA == 3 || fExtendedMatchAndQA == 5 )){
        if(inTrack->Charge() > 0) {
          fHistClusterdEtadPhiPosTracksBeforeQA->Fill(dEta, dPhi, weight);
          fHistClusterdPhidPtPosTracksBeforeQA->Fill(dPhi, inTrack->Pt(), weight);
          if(inTrack->P() < 0.75) fHistClusterdEtadPhiPosTracksP_000_075BeforeQA->Fill(dEta, dPhi, weight);
          else if(inTrack->P() < 1.25) fHistClusterdEtadPhiPosTracksP_075_125BeforeQA->Fill(dEta, dPhi, weight);
          else fHistClusterdEtadPhiPosTracksP_125_999BeforeQA->Fill(dEta, dPhi, weight);
        }
        else{
          fHistClusterdEtadPhiNegTracksBeforeQA->Fill(dEta, dPhi, weight);
          fHistClusterdPhidPtNegTracksBeforeQA->Fill(dPhi, inTrack->Pt(), weight);
          if(inTrack->P() < 0.75) fHistClusterdEtadPhiNegTracksP_000_075BeforeQA->Fill(dEta, dPhi, weight);
          else if(inTrack->P() < 1.25) fHistClusterdEtadPhiNegTracksP_075_125BeforeQA->Fill(dEta, dPhi, weight);
          else fHistClusterdEtadPhiNegTracksP_125_999BeforeQA->Fill(dEta, dPhi, weight);
        }
        fHistClusterdEtadPtBeforeQA->Fill(dEta, inTrack->Pt(), weight);
        fHistClusterM20M02BeforeQA->Fill(clusM20, clusM02, weight);
      }

      Bool_t match_dEta = (TMath::Abs(dEta) < fMaxDistTrackToClusterEta) ? kTRUE : kFALSE;
      Bool_t match_dPhi = kFALSE;
      if( (inTrack->Charge() > 0) && (dPhi > fMinDistTrackToClusterPhi) && (dPhi < fMaxDistTrackToClusterPhi) ) match_dPhi = kTRUE;
      else if( (inTrack->Charge() < 0) && (dPhi < -fMinDistTrackToClusterPhi) && (dPhi > -fMaxDistTrackToClusterPhi) ) match_dPhi = kTRUE;

      if(fUsePtDepTrackToCluster){
        if( TMath::Abs(dEta) < fFuncPtDepEta->Eval(inTrack->Pt())) match_dEta = kTRUE;
        else match_dEta = kFALSE;

        if( TMath::Abs(dPhi) < fFuncPtDepPhi->Eval(inTrack->Pt())) match_dPhi = kTRUE;
        else match_dPhi = kFALSE;
      }

      if(match_dEta && match_dPhi){
        fVectorMatchedClusterIDs.push_back(cluster->GetID());
//        cout << "MATCHED!!!!!!!!!!!!!!!!!!!!!!!!! - " << cluster->GetID() << endl;
        if(isEMCalOnly){
          if(fHistClusterdEtadPtAfterQA) fHistClusterdEtadPtAfterQA->Fill(dEta,inTrack->Pt());
          if(fHistClusterdPhidPtAfterQA) fHistClusterdPhidPtAfterQA->Fill(dPhi,inTrack->Pt());
        }
        break;
      } else if(isEMCalOnly){
        if(fHistDistanceTrackToClusterAfterQA)fHistDistanceTrackToClusterAfterQA->Fill(TMath::Sqrt(dR2), weight);
        if(fHistClusterdEtadPhiAfterQA) fHistClusterdEtadPhiAfterQA->Fill(dEta, dPhi, weight);
        if(fHistClusterRAfterQA) fHistClusterRAfterQA->Fill(clusterR, weight);
        if(!fDoLightOutput && (fExtendedMatchAndQA == 1 || fExtendedMatchAndQA == 3 || fExtendedMatchAndQA == 5 )){
          if(inTrack->Charge() > 0) fHistClusterdEtadPhiPosTracksAfterQA->Fill(dEta, dPhi, weight);
          else fHistClusterdEtadPhiNegTracksAfterQA->Fill(dEta, dPhi, weight);
          fHistClusterM20M02AfterQA->Fill(clusM20, clusM02, weight);
        }
//        cout << "no match" << endl;
      }
    }

    delete trackParam;
  }

  if(EsdTrackCuts){
    delete EsdTrackCuts;
    EsdTrackCuts=0x0;
  }

  return;
}

//________________________________________________________________________
Bool_t AliCaloPhotonCuts::CheckClusterForTrackMatch(AliVCluster* cluster){
  vector<Int_t>::iterator it;
  it = find (fVectorMatchedClusterIDs.begin(), fVectorMatchedClusterIDs.end(), cluster->GetID());
  if (it != fVectorMatchedClusterIDs.end()) return kTRUE;
  else return kFALSE;
}

//________________________________________________________________________
Bool_t AliCaloPhotonCuts::UpdateCutString() {
   ///Update the cut string (if it has been created yet)

   if(fCutString && fCutString->GetString().Length() == kNCuts) {
      fCutString->SetString(GetCutNumber());
   } else {
      return kFALSE;
   }
   return kTRUE;
}

//________________________________________________________________________
Bool_t AliCaloPhotonCuts::InitializeCutsFromCutString(const TString analysisCutSelection ) {
  // Initialize Cuts from a given Cut string
  AliInfo(Form("Set CaloCut Number: %s",analysisCutSelection.Data()));
  if(analysisCutSelection.Length()!=kNCuts) {
    AliError(Form("Cut selection has the wrong length! size is %d, number of cuts is %d", analysisCutSelection.Length(), kNCuts));
    return kFALSE;
  }
  if(!analysisCutSelection.IsDigit()){
    AliError("Cut selection contains characters");
    return kFALSE;
  }

  const char *cutSelection = analysisCutSelection.Data();
  #define ASSIGNARRAY(i)  fCuts[i] = cutSelection[i] - '0'
  for(Int_t ii=0;ii<kNCuts;ii++){
    ASSIGNARRAY(ii);
  }

  // Set Individual Cuts
  for(Int_t ii=0;ii<kNCuts;ii++){
    if(!SetCut(cutIds(ii),fCuts[ii]))return kFALSE;
  }
  PrintCutsWithValues();
  return kTRUE;
}

//________________________________________________________________________
Bool_t AliCaloPhotonCuts::SetCut(cutIds cutID, const Int_t value) {
  ///Set individual cut ID

  switch (cutID) {    
    
    case kClusterType:
      if( SetClusterTypeCut(value)) {
        fCuts[kClusterType] = value;
        UpdateCutString();
        return kTRUE;
      } else return kFALSE;
    
    case kEtaMin:
      if( SetMinEtaCut(value)) {
        fCuts[kEtaMin] = value;
        UpdateCutString();
        return kTRUE;
      } else return kFALSE;

    case kEtaMax:
      if( SetMaxEtaCut(value)) {
        fCuts[kEtaMax] = value;
        UpdateCutString();
        return kTRUE;
      } else return kFALSE;

    case kPhiMin:
      if( SetMinPhiCut(value)) {
        fCuts[kPhiMin] = value;
        UpdateCutString();
        return kTRUE;
      } else return kFALSE;

    case kPhiMax:
      if( SetMaxPhiCut(value)) {
        fCuts[kPhiMax] = value;
        UpdateCutString();
        return kTRUE;
      } else return kFALSE;

    case kDistanceToBadChannel:
      if( SetDistanceToBadChannelCut(value)) {
        fCuts[kDistanceToBadChannel] = value;
        UpdateCutString();
        return kTRUE;
      } else return kFALSE;

    case kTiming:
      if( SetTimingCut(value)) {
        fCuts[kTiming] = value;
        UpdateCutString();
        return kTRUE;
      } else return kFALSE;

    case kTrackMatching:
      if( SetTrackMatchingCut(value)) {
        fCuts[kTrackMatching] = value;
        UpdateCutString();
        return kTRUE;
      } else return kFALSE;

    case kExoticCluster:
      if( SetExoticClusterCut(value)) {
        fCuts[kExoticCluster] = value;
        UpdateCutString();
        return kTRUE;
      } else return kFALSE;

    case kMinEnergy:
      if( SetMinEnergyCut(value)) {
        fCuts[kMinEnergy] = value;
        UpdateCutString();
        return kTRUE;
      } else return kFALSE;

    case kNMinCells:
      if( SetMinNCellsCut(value)) {
        fCuts[kNMinCells] = value;
        UpdateCutString();
        return kTRUE;
      } else return kFALSE;
      
    case kMinM02:
      if( SetMinM02(value)) {
        fCuts[kMinM02] = value;
        UpdateCutString();
        return kTRUE;
      } else return kFALSE;

    case kMaxM02:
      if( SetMaxM02(value)) {
        fCuts[kMaxM02] = value;
        UpdateCutString();
        return kTRUE;
      } else return kFALSE;
    
    case kMinM20:
      if( SetMinM20(value)) {
        fCuts[kMinM20] = value;
        UpdateCutString();
        return kTRUE;
      } else return kFALSE;

    case kMaxM20:
      if( SetMaxM20(value)) {
        fCuts[kMaxM20] = value;
        UpdateCutString();
        return kTRUE;
      } else return kFALSE;

    case kDispersion:
      if( SetDispersion(value)) {
        fCuts[kDispersion] = value;
        UpdateCutString();
        return kTRUE;
      } else return kFALSE;

    case kNLM:
      if( SetNLM(value)) {
        fCuts[kNLM] = value;
        UpdateCutString();
        return kTRUE;
      } else return kFALSE;

    case kNonLinearity1:
      if( SetNonLinearity1(value)) {
        fCuts[kNonLinearity1] = value;
        UpdateCutString();
        return kTRUE;
      } else return kFALSE;

    case kNonLinearity2:
      if( SetNonLinearity2(value)) {
        fCuts[kNonLinearity2] = value;
        UpdateCutString();
        return kTRUE;
      } else return kFALSE;

    case kNCuts:
      AliError("Cut id out of range");
      return kFALSE;
  }

  AliError("Cut id %d not recognized");
  return kFALSE;


}

//________________________________________________________________________
void AliCaloPhotonCuts::PrintCuts() {
   // Print out current Cut Selection
   for(Int_t ic = 0;ic < kNCuts;ic++) {
      printf("%-30s : %d \n", fgkCutNames[ic], fCuts[ic]);
   }
}

//________________________________________________________________________
void AliCaloPhotonCuts::PrintCutsWithValues() {
  // Print out current Cut Selection with value
  printf("\nCluster cutnumber \n");
  for(Int_t ic = 0;ic < kNCuts;ic++) {
    printf("%d",fCuts[ic]);
  }
  printf("\n\n");
  if (fIsPureCalo>0) printf("Merged cluster analysis was specified, mode: '%i'\n", fIsPureCalo);
  
  printf("Acceptance cuts: \n");
  if (fClusterType == 0) printf("\tall calorimeter clusters are used\n");
  if (fClusterType == 1) printf("\tEMCAL calorimeter clusters are used\n");
  if (fClusterType == 2) printf("\tPHOS calorimeter clusters are used\n");
  if (fClusterType == 3) printf("\tDCAL calorimeter clusters are used\n");
  if (fUseEtaCut) printf("\t%3.2f < eta_{cluster} < %3.2f\n", fMinEtaCut, fMaxEtaCut );
  if (fUsePhiCut) printf("\t%3.2f < phi_{cluster} < %3.2f\n", fMinPhiCut, fMaxPhiCut );
  if (fUseDistanceToBadChannel>0) printf("\tdistance to bad channel used in mode '%i', distance in cells: %f \n",fUseDistanceToBadChannel, fMinDistanceToBadChannel);
  
  printf("Cluster Quality cuts: \n");
  if (fUseTimeDiff) printf("\t %6.2f ns < time difference < %6.2f ns\n", fMinTimeDiff*1e9, fMaxTimeDiff*1e9 );
  if (fUseDistTrackToCluster) printf("\tmin distance to track in eta > %3.2f, min phi < %3.2f and max phi > %3.2f\n", fMaxDistTrackToClusterEta, fMinDistTrackToClusterPhi, fMaxDistTrackToClusterPhi );
  if (fUseExoticCluster)printf("\t exotic cluster: %3.2f\n", fExoticEnergyFracCluster );
  if (fUseMinEnergy)printf("\t E_{cluster} > %3.2f\n", fMinEnergy );
  if (fUseNCells) printf("\t number of cells per cluster >= %d\n", fMinNCells );
  if (fUseM02 == 1) printf("\t %3.2f < M02 < %3.2f\n", fMinM02, fMaxM02 );
  if (fUseM02 == 2) printf("\t energy dependent M02 cut used with cutnumber min: %d  max: %d \n", fMinM02CutNr, fMaxM02CutNr );
  if (fUseM20) printf("\t %3.2f < M20 < %3.2f\n", fMinM20, fMaxM20 );
  if (fUseDispersion) printf("\t dispersion < %3.2f\n", fMaxDispersion );
  if (fUseNLM) printf("\t %d < NLM < %d\n", fMinNLM, fMaxNLM );

  printf("NonLinearity Correction: \n");
  printf("VO Reader name: %s \n",fV0ReaderName.Data());
  TString periodName = ((AliV0ReaderV1*)AliAnalysisManager::GetAnalysisManager()->GetTask(fV0ReaderName.Data()))->GetPeriodName();
  if (periodName.CompareTo("") != 0) fCurrentMC = FindEnumForMCSet(periodName);
  if (fUseNonLinearity) printf("\t Chose NonLinearity cut '%i', Period name: %s, period-enum: %o \n", fSwitchNonLinearity, periodName.Data(), fCurrentMC );
  else printf("\t No NonLinearity Correction on AnalysisTask level has been chosen\n");
  
}

// EMCAL acceptance 2011
// 1.39626, 3.125 (phi)
// -0.66687,,0.66465


//________________________________________________________________________
Bool_t AliCaloPhotonCuts::SetClusterTypeCut(Int_t clusterType)
{   // Set Cut
  switch(clusterType){
  case 0: // all clusters
    fClusterType=0;
    break;
  case 1: // EMCAL clusters
    fClusterType=1;
    break;
  case 2: // PHOS clusters
    fClusterType=2;
    break;
  case 3: // DCAL clusters
    fClusterType=3;
    break;
  default:
    AliError(Form("ClusterTypeCut not defined %d",clusterType));
    return kFALSE;
  }
  return kTRUE;
}

//___________________________________________________________________
Bool_t AliCaloPhotonCuts::SetMinEtaCut(Int_t minEta)
{
  switch(minEta){
  case 0:
    if (!fUseEtaCut) fUseEtaCut=0;
    fMinEtaCut=-10.;
    break;
  case 1:
    if (!fUseEtaCut) fUseEtaCut=1;
    fMinEtaCut=-0.6687;
    break;
  case 2: 
    if (!fUseEtaCut) fUseEtaCut=1;
    fMinEtaCut=-0.5;
    break;
  case 3: 
    if (!fUseEtaCut) fUseEtaCut=1;
    fMinEtaCut=-2;
    break;
  case 4: 
    if (!fUseEtaCut) fUseEtaCut=1;
    fMinEtaCut = -0.13;
    break;
  case 5:
    if (!fUseEtaCut) fUseEtaCut=1;
    fMinEtaCut=-0.7;
    break;
  case 6:
    if (!fUseEtaCut) fUseEtaCut=1;
    fMinEtaCut=-0.3;
    break;
  case 7:
    if (!fUseEtaCut) fUseEtaCut=1;
    fMinEtaCut=-0.4;
    break;
    
  default:
    AliError(Form("MinEta Cut not defined %d",minEta));
    return kFALSE;
  }
  return kTRUE;
}


//___________________________________________________________________
Bool_t AliCaloPhotonCuts::SetMaxEtaCut(Int_t maxEta)
{
  switch(maxEta){
  case 0: 
    if (!fUseEtaCut) fUseEtaCut=0;
    fMaxEtaCut=10;
    break;
  case 1:
    if (!fUseEtaCut) fUseEtaCut=1;
    fMaxEtaCut=0.66465;
    break;
  case 2: 
    if (!fUseEtaCut) fUseEtaCut=1;
    fMaxEtaCut=0.5;
    break;
  case 3: 
    if (!fUseEtaCut) fUseEtaCut=1;
    fMaxEtaCut=2;
    break;
  case 4: 
    if (!fUseEtaCut) fUseEtaCut=1;
    fMaxEtaCut= 0.13;
    break;    
  case 5:
    if (!fUseEtaCut) fUseEtaCut=1;
    fMaxEtaCut=0.7;
    break;
  case 6:
    if (!fUseEtaCut) fUseEtaCut=1;
    fMaxEtaCut=0.3;
    break;
  case 7:
    if (!fUseEtaCut) fUseEtaCut=1;
    fMaxEtaCut=0.4;
    break;
  default:
    AliError(Form("MaxEta Cut not defined %d",maxEta));
    return kFALSE;
  }
  return kTRUE;
}

//___________________________________________________________________
Bool_t AliCaloPhotonCuts::SetMinPhiCut(Int_t minPhi)
{
  switch(minPhi){
  case 0: 
    if (!fUsePhiCut) fUsePhiCut=0;
    fMinPhiCut=-10000;
    break;
  case 1: // min EMCAL
    if (!fUsePhiCut) fUsePhiCut=1;
    fMinPhiCut=1.39626;
    break;
  case 2: // min EMCAL with TRD 2012 
    if (!fUsePhiCut) fUsePhiCut=1;
    fMinPhiCut=2.10;
    break;
  case 3: // min EMCAL with TRD 2011 
    if (!fUsePhiCut) fUsePhiCut=1;
    fMinPhiCut=2.45;
    break;
  case 4: 
    if( !fUsePhiCut ) fUsePhiCut=1;
    fMinPhiCut = 4.54;//PHOS acceptance
    break;
  case 5:
    if( !fUsePhiCut ) fUsePhiCut=1;
    fMinPhiCut = 4.54;//DCal acceptance
    break;

  default:
    AliError(Form("MinPhi Cut not defined %d",minPhi));
    return kFALSE;
  }
  return kTRUE;
}

//___________________________________________________________________
Bool_t AliCaloPhotonCuts::SetMaxPhiCut(Int_t maxPhi)
{
  switch(maxPhi){
  case 0: 
    if (!fUsePhiCut) fUsePhiCut=0;
    fMaxPhiCut=10000;
    break;
  case 1: // max EMCAL
    if (!fUsePhiCut) fUsePhiCut=1;
    fMaxPhiCut=3.15;
    break;
  case 2: // max EMCAL with TRD 2011
    if (!fUsePhiCut) fUsePhiCut=1;
    fMaxPhiCut=2.45;
    break;
  case 3: // max EMCAL with TRD 2012
    if (!fUsePhiCut) fUsePhiCut=1;
    fMaxPhiCut=2.10;
    break;
  case 4: 
    if( !fUsePhiCut ) fUsePhiCut=1;
    fMaxPhiCut = 5.59;//PHOS acceptance
    break;
  case 5:
    if( !fUsePhiCut ) fUsePhiCut=1;
    fMaxPhiCut = 5.59;//DCal acceptance
    break;
    
  default:
    AliError(Form("Max Phi Cut not defined %d",maxPhi));
    return kFALSE;
  }
  return kTRUE;
}

//___________________________________________________________________
Bool_t AliCaloPhotonCuts::SetDistanceToBadChannelCut(Int_t distanceToBadChannel)
{
  switch(distanceToBadChannel){
  case 0: 
    fUseDistanceToBadChannel=0;
    fMinDistanceToBadChannel=0;
    break;
  case 1: 
    if(fUseDistanceToBadChannel!=1) fUseDistanceToBadChannel=1;
    fMinDistanceToBadChannel=1;
    break;
  case 2:
    if(fUseDistanceToBadChannel!=1) fUseDistanceToBadChannel=1;
    fMinDistanceToBadChannel=2;
    break;
  case 3:
    if(fUseDistanceToBadChannel!=1) fUseDistanceToBadChannel=1;
    fMinDistanceToBadChannel=3;
    break;
  case 4:
    if(fUseDistanceToBadChannel!=1) fUseDistanceToBadChannel=1;
    fMinDistanceToBadChannel=4;
    break;
  case 5:
    if(fUseDistanceToBadChannel!=2) fUseDistanceToBadChannel=2;
    fMinDistanceToBadChannel=1;
    break;
  case 6:
    if(fUseDistanceToBadChannel!=2) fUseDistanceToBadChannel=2;
    fMinDistanceToBadChannel=2;
    break;
  case 7:
    if(fUseDistanceToBadChannel!=2) fUseDistanceToBadChannel=2;
    fMinDistanceToBadChannel=3;
    break;
  case 8:
    if(fUseDistanceToBadChannel!=2) fUseDistanceToBadChannel=2;
    fMinDistanceToBadChannel=4;
    break;
  default:
    AliError(Form("minimum distance to bad channel Cut not defined %d",distanceToBadChannel));
    return kFALSE;
  }
  return kTRUE;
}

//___________________________________________________________________
Bool_t AliCaloPhotonCuts::SetTimingCut(Int_t timing)
{
  switch(timing){
  case 0: 
    fUseTimeDiff=0;
    fMinTimeDiff=-500;
    fMaxTimeDiff=500;
    break;
  case 1: 
    if (!fUseTimeDiff) fUseTimeDiff=1;
    fMinTimeDiff=-10e-7;
    fMaxTimeDiff=10e-7;//1000ns
    break;
  case 2: 
    if (!fUseTimeDiff) fUseTimeDiff=1;
    fMinTimeDiff=-50e-8;
    fMaxTimeDiff=50e-8;//500ns
    break;
  case 3: 
    if (!fUseTimeDiff) fUseTimeDiff=1;
    fMinTimeDiff=-20e-8;
    fMaxTimeDiff=20e-8;//200ns
    break;
  case 4: 
    if (!fUseTimeDiff) fUseTimeDiff=1;
    fMinTimeDiff=-10e-8;
    fMaxTimeDiff=10e-8;//100ns
    break;
  case 5: 
    if (!fUseTimeDiff) fUseTimeDiff=1;
    fMinTimeDiff=-50e-9;
    fMaxTimeDiff=50e-9;//50ns
    break;
  case 6:
    if (!fUseTimeDiff) fUseTimeDiff=1;
    fMinTimeDiff=-30e-9;
    fMaxTimeDiff=35e-9;
    break;
  case 7:
    if (!fUseTimeDiff) fUseTimeDiff=1;
    fMinTimeDiff=-30e-9;
    fMaxTimeDiff=30e-9;//30ns
    break;
  case 8:
    if (!fUseTimeDiff) fUseTimeDiff=1;
    fMinTimeDiff=-20e-9;
    fMaxTimeDiff=30e-9;
    break;
  case 9:
    if (!fUseTimeDiff) fUseTimeDiff=1;
    fMinTimeDiff=-20e-9;
    fMaxTimeDiff=25e-9;
    break;

  default:
    AliError(Form("Timing Cut not defined %d",timing));
    return kFALSE;
  }
  return kTRUE;
}

//___________________________________________________________________
Bool_t AliCaloPhotonCuts::SetTrackMatchingCut(Int_t trackMatching)
{
  // matching parameters for EMCal clusters
  if(fClusterType == 1 || fClusterType == 3){
    switch(trackMatching){
    case 0:
      fUseDistTrackToCluster = kFALSE;
      fMaxDistTrackToClusterEta = 0;
      fMinDistTrackToClusterPhi = 0;
      fMaxDistTrackToClusterPhi = 0;
      break;
    case 1:
      if (!fUseDistTrackToCluster) fUseDistTrackToCluster=kTRUE;
      fMaxDistTrackToClusterEta = 0.008;//0.015;
      fMinDistTrackToClusterPhi = -0.03;//-0.01;
      fMaxDistTrackToClusterPhi = 0.03;//0.03;//0.04;
      break;
    case 2:
      if (!fUseDistTrackToCluster) fUseDistTrackToCluster=kTRUE;
      fMaxDistTrackToClusterEta = 0.012;//0.015;
      fMinDistTrackToClusterPhi = -0.05;//-0.01;
      fMaxDistTrackToClusterPhi = 0.04;//0.035;//0.05;
      break;
    case 3:
      if (!fUseDistTrackToCluster) fUseDistTrackToCluster=kTRUE;
      fMaxDistTrackToClusterEta = 0.016;//0.015;
      fMinDistTrackToClusterPhi = -0.09;//-0.015;
      fMaxDistTrackToClusterPhi = 0.06;//0.04;//0.1;
      break;
    case 4:
      if (!fUseDistTrackToCluster) fUseDistTrackToCluster=kTRUE;
      fMaxDistTrackToClusterEta = 0.018;//0.015;
      fMinDistTrackToClusterPhi = -0.11;//-0.015;
      fMaxDistTrackToClusterPhi = 0.07;//0.045;//0.13;
      break;
    case 5:
      if (!fUseDistTrackToCluster) fUseDistTrackToCluster=kTRUE;
      fMaxDistTrackToClusterEta = 0.02;//0.015;
      fMinDistTrackToClusterPhi = -0.13;//-0.02;
      fMaxDistTrackToClusterPhi = 0.08;//0.05;//0.15
      break;
//    case 6:
//      if (!fUseDistTrackToCluster) fUseDistTrackToCluster=kTRUE;
//      fMaxDistTrackToClusterEta = 0.022;//0.015;
//      fMinDistTrackToClusterPhi = -0.15;//-0.02;
//      fMaxDistTrackToClusterPhi = 0.10;//0.055;//0.2;
//      break;
    // pT dependent matching parameters
    case 6:
      if (!fUseDistTrackToCluster) fUseDistTrackToCluster=kTRUE;
      fUsePtDepTrackToCluster = kTRUE;
      fFuncPtDepEta = new TF1("func", "[1] + 1 / pow(x + pow(1 / ([0] - [1]), 1 / [2]), [2])");
      fFuncPtDepEta->SetParameters(0.03, 0.010, 2.5);

      fFuncPtDepPhi = new TF1("func", "[1] + 1 / pow(x + pow(1 / ([0] - [1]), 1 / [2]), [2])");
      fFuncPtDepPhi->SetParameters(0.08, 0.015, 2.);
      break;
    case 7:
      if (!fUseDistTrackToCluster) fUseDistTrackToCluster=kTRUE;
      fUsePtDepTrackToCluster = kTRUE;
      fFuncPtDepEta = new TF1("func", "[1] + 1 / pow(x + pow(1 / ([0] - [1]), 1 / [2]), [2])");
      fFuncPtDepEta->SetParameters(0.04, 0.010, 2.5);

      fFuncPtDepPhi = new TF1("func", "[1] + 1 / pow(x + pow(1 / ([0] - [1]), 1 / [2]), [2])");
      fFuncPtDepPhi->SetParameters(0.09, 0.015, 2.);
      break;
    case 8:
      if (!fUseDistTrackToCluster) fUseDistTrackToCluster=kTRUE;
      fUsePtDepTrackToCluster = kTRUE;
      fFuncPtDepEta = new TF1("func", "[1] + 1 / pow(x + pow(1 / ([0] - [1]), 1 / [2]), [2])");
      fFuncPtDepEta->SetParameters(0.05, 0.010, 2.5);

      fFuncPtDepPhi = new TF1("func", "[1] + 1 / pow(x + pow(1 / ([0] - [1]), 1 / [2]), [2])");
      fFuncPtDepPhi->SetParameters(0.10, 0.015, 1.75);
      break;
    case 9:
      if (!fUseDistTrackToCluster) fUseDistTrackToCluster=kTRUE;
      fUsePtDepTrackToCluster = kTRUE;
      fFuncPtDepEta = new TF1("func", "[1] + 1 / pow(x + pow(1 / ([0] - [1]), 1 / [2]), [2])");
      fFuncPtDepEta->SetParameters(0.06, 0.015, 2.5);

      fFuncPtDepPhi = new TF1("func", "[1] + 1 / pow(x + pow(1 / ([0] - [1]), 1 / [2]), [2])");
      fFuncPtDepPhi->SetParameters(0.12, 0.020, 1.75);
      break;

    default:
      AliError(Form("Track Matching Cut not defined %d",trackMatching));
      return kFALSE;
    }
    return kTRUE;
  // matching parameters for PHOS clusters
  }else if(fClusterType == 2) {
    switch(trackMatching){
      case 0:
        fUseDistTrackToCluster = kFALSE;
        fMaxDistTrackToClusterEta = 0;
        fMinDistTrackToClusterPhi = 0;
        fMaxDistTrackToClusterPhi = 0;
        break;
      case 1:
        if (!fUseDistTrackToCluster) fUseDistTrackToCluster=kTRUE;
        fMaxDistTrackToClusterEta = 0.005;//0.015;
        fMinDistTrackToClusterPhi = -0.03;//-0.025;
        fMaxDistTrackToClusterPhi = 0.03;//0.06;//0.3;
        break;
      case 2:
        if (!fUseDistTrackToCluster) fUseDistTrackToCluster=kTRUE;
        fMaxDistTrackToClusterEta = 0.01;//0.015;
        fMinDistTrackToClusterPhi = -0.09;//-0.025;
        fMaxDistTrackToClusterPhi = 0.07;//0.07;//0.4;
        break;
      case 3:
        if (!fUseDistTrackToCluster) fUseDistTrackToCluster=kTRUE;
        fMaxDistTrackToClusterEta = 0.015;//0.02;
        fMinDistTrackToClusterPhi = -0.15;//-0.03;
        fMaxDistTrackToClusterPhi = 0.11;//0.1;//0.5;
        break;

      default:
        AliError(Form("Track Matching Cut not defined %d",trackMatching));
        return kFALSE;
    }
    return kTRUE;
  }
  return kTRUE;
}

//___________________________________________________________________
Bool_t AliCaloPhotonCuts::SetExoticClusterCut(Int_t exoticCell)
{
  switch(exoticCell){
  case 0: 
    fUseExoticCluster         = 0;
    fExoticEnergyFracCluster  = 0;
    break;
  case 1: 
    if (!fUseExoticCluster) 
      fUseExoticCluster       = 1;
    fExoticEnergyFracCluster  = 0.995;
    break;
  case 2: 
    if (!fUseExoticCluster) 
      fUseExoticCluster       = 1;
    fExoticEnergyFracCluster  = 0.99;
    break;
  case 3: 
    if (!fUseExoticCluster) 
      fUseExoticCluster       = 1;
    fExoticEnergyFracCluster  = 0.98;
    break;
  case 4: 
    if (!fUseExoticCluster) 
      fUseExoticCluster       = 1;
    fExoticEnergyFracCluster  = 0.975;
    break;
  case 5: 
    if (!fUseExoticCluster) 
      fUseExoticCluster       = 1;
    fExoticEnergyFracCluster  = 0.97;
    break;
  case 6: 
    if (!fUseExoticCluster) 
      fUseExoticCluster       = 1;
    fExoticEnergyFracCluster  = 0.965;
    break;
  case 7: 
    if (!fUseExoticCluster) 
      fUseExoticCluster       = 1;
    fExoticEnergyFracCluster  = 0.96;
    break;
  case 8: 
    if (!fUseExoticCluster) 
      fUseExoticCluster       = 1;
    fExoticEnergyFracCluster  = 0.95;
    break;
  case 9: 
    if (!fUseExoticCluster) 
      fUseExoticCluster       = 1;
    fExoticEnergyFracCluster  = 0.94;
    break;
  default:
    AliError(Form("Exotic cell Cut not defined %d",exoticCell));
    return kFALSE;
  }
  return kTRUE;
}
    
//___________________________________________________________________
Bool_t AliCaloPhotonCuts::SetMinEnergyCut(Int_t minEnergy)
{
  if (fIsPureCalo != 1){
    if (fClusterType!=2) {
      switch(minEnergy){
        case 0:
          if (!fUseMinEnergy) fUseMinEnergy=0;
          fMinEnergy=0.1;
          break;
        case 1:
          if (!fUseMinEnergy) fUseMinEnergy=1;
          fMinEnergy=0.5;
          break;
        case 2:
          if (!fUseMinEnergy) fUseMinEnergy=1;
          fMinEnergy=0.6;
          break;
        case 3:
          if (!fUseMinEnergy) fUseMinEnergy=1;
          fMinEnergy=0.7;
          break;
        case 4:
          if (!fUseMinEnergy) fUseMinEnergy=1;
          fMinEnergy=0.8;
          break;
        case 5:
          if (!fUseMinEnergy) fUseMinEnergy=1;
          fMinEnergy=0.9;
          break;
        case 6:
          if (!fUseMinEnergy) fUseMinEnergy=1;
          fMinEnergy=4.5;
          break;
        case 7:
          if (!fUseMinEnergy) fUseMinEnergy=1;
          fMinEnergy=5.0;
          break;
        case 8:
          if (!fUseMinEnergy) fUseMinEnergy=1;
          fMinEnergy=5.5;
          break;
        case 9:
          if (!fUseMinEnergy) fUseMinEnergy=1;
          fMinEnergy=6.0;
          break;
        default:
          AliError(Form("Minimum Energy Cut not defined %d",minEnergy));
          return kFALSE;
      }
    } else {
      switch(minEnergy){
        case 0:
          if (!fUseMinEnergy) fUseMinEnergy=0;
          fMinEnergy=0.1;
          break;
        case 1:
          if (!fUseMinEnergy) fUseMinEnergy=1;
          fMinEnergy=0.3;
          break;
        case 2:
          if (!fUseMinEnergy) fUseMinEnergy=1;
          fMinEnergy=0.5;
          break;
        case 3:
          if (!fUseMinEnergy) fUseMinEnergy=1;
          fMinEnergy=0.6;
          break;
        case 4:
          if (!fUseMinEnergy) fUseMinEnergy=1;
          fMinEnergy=0.7;
          break;
        case 5:
          if (!fUseMinEnergy) fUseMinEnergy=1;
          fMinEnergy=0.8;
          break;
        case 6:
          if (!fUseMinEnergy) fUseMinEnergy=1;
          fMinEnergy=0.9;
          break;
        default:
          AliError(Form("Minimum Energy Cut not defined %d",minEnergy));
          return kFALSE;
      }
    }
    return kTRUE;
  } else {
    switch(minEnergy){
    case 0:
      if (!fUseMinEnergy) fUseMinEnergy=0;
      fMinEnergy=0.1;
      break;
    case 1: 
      if (!fUseMinEnergy) fUseMinEnergy=1;
      fMinEnergy=4.;
      break;
    case 2: 
      if (!fUseMinEnergy) fUseMinEnergy=1;
      fMinEnergy=5.;
      break;
    case 3: 
      if (!fUseMinEnergy) fUseMinEnergy=1;
      fMinEnergy=6.;
      break;
    case 4: 
      if (!fUseMinEnergy) fUseMinEnergy=1;
      fMinEnergy=7.;
      break;
    case 5: 
      if (!fUseMinEnergy) fUseMinEnergy=1;
      fMinEnergy=7.5;
      break;
    case 6: 
      if (!fUseMinEnergy) fUseMinEnergy=1;
      fMinEnergy=8.;
      break;
    case 7: 
      if (!fUseMinEnergy) fUseMinEnergy=1;
      fMinEnergy=8.5;
      break;
    case 8: 
      if (!fUseMinEnergy) fUseMinEnergy=1;
      fMinEnergy=9.;
      break;
    case 9: 
      if (!fUseMinEnergy) fUseMinEnergy=1;
      fMinEnergy=9.5;
      break;
    default:
      AliError(Form("Minimum Energy Cut not defined %d",minEnergy));
      return kFALSE;
    }
    return kTRUE;
  }
  return kTRUE;
}
    
//___________________________________________________________________
Bool_t AliCaloPhotonCuts::SetMinNCellsCut(Int_t minNCells)
{
  switch(minNCells){
  case 0:
    if (!fUseNCells) fUseNCells=0;
    fMinNCells=0;
    break;
  case 1: 
    if (!fUseNCells) fUseNCells=1;
    fMinNCells=1;
    break;
  case 2: 
    if (!fUseNCells) fUseNCells=1;
    fMinNCells=2;
    break;
  case 3: 
    if (!fUseNCells) fUseNCells=1;
    fMinNCells=3;
    break;
  case 4: 
    if (!fUseNCells) fUseNCells=1;
    fMinNCells=4;
    break;
  case 5: 
    if (!fUseNCells) fUseNCells=1;
    fMinNCells=5;
    break;
  case 6: 
    if (!fUseNCells) fUseNCells=1;
    fMinNCells=6;
    break;

  default:
    AliError(Form("Min N cells Cut not defined %d",minNCells));
    return kFALSE;
  }
  return kTRUE;
}

//___________________________________________________________________
Bool_t AliCaloPhotonCuts::SetMaxM02(Int_t maxM02)
{
  fMaxM02CutNr = maxM02;
  if (fIsPureCalo == 1){
    fUseM02 = 2;
    return kTRUE;
  }
  
  switch(maxM02){
  case 0: 
    if (!fUseM02) fUseM02=0;
    fMaxM02=100;
    break;
  case 1: 
    if (!fUseM02) fUseM02=1;
    fMaxM02=1.;
    break;
  case 2: 
    if (!fUseM02) fUseM02=1;
    fMaxM02=0.7;
    break;
  case 3: 
    if (!fUseM02) fUseM02=1;
    fMaxM02=0.5;
    break;
  case 4: 
    if (!fUseM02) fUseM02=1;
    fMaxM02=0.4;
    break;
  case 5: 
    if (!fUseM02) fUseM02=1;
    fMaxM02=0.3;
    break;
  case 6: 
    if (!fUseM02) fUseM02=1;
    fMaxM02=0.27;
    break;
  default:
    AliError(Form("Max M02 Cut not defined %d",maxM02));
    return kFALSE;
  }
  return kTRUE;
}

//___________________________________________________________________
Float_t AliCaloPhotonCuts::CalculateMaxM02 (Int_t maxM02, Float_t clusEnergy){
  switch (maxM02){
    case 0: 
      return 10;
    case 1:
      if (fMinNLM == 1 && fMaxNLM == 1 ){
        return FunctionM02(clusEnergy, 0.0662, -0.0201, -0.0955, 1.86e-3, 9.91 );
      } else if (fMinNLM == 2 && fMaxNLM == 2 ){
        return FunctionM02(clusEnergy, 0.353, -0.0264, -0.524, 5.59e-3, 21.9 );
      } else {
        return 10;
      }  
    case 2:
      if (fMinNLM == 1 && fMaxNLM == 1 ){
        return FunctionM02(clusEnergy, 0.0662, -0.0201, -0.0, 1.86e-3, 9.91 );
      } else if (fMinNLM == 2 && fMaxNLM == 2 ){
        return FunctionM02(clusEnergy, 0.353, -0.0264, -0.424, 5.59e-3, 21.9 );
      } else {
        return 10;
      }  
    case 3:
      if (fMinNLM == 1 && fMaxNLM == 1 ){
        return FunctionM02(clusEnergy, 0.0662, -0.0201, -0.2, 1.86e-3, 9.91 );
      } else if (fMinNLM == 2 && fMaxNLM == 2 ){
        return FunctionM02(clusEnergy, 0.353, -0.0264, -0.624, 5.59e-3, 21.9 );
      } else {
        return 10;
      }  
      
    default:  
      AliError(Form("Max M02 for merged cluster Cut not defined %d",maxM02));
      return 10;
  }
  return 10;
  
}  

//___________________________________________________________________
Float_t AliCaloPhotonCuts::CalculateMinM02 (Int_t minM02, Float_t clusEnergy){
  switch (minM02){
    case 0: 
      return 0.;
    case 1:
      if (FunctionM02(clusEnergy, 2.135, -0.245, 0., 0., 0. ) > 0.3)
        return FunctionM02(clusEnergy, 2.135, -0.245, 0., 0., 0. );
      else 
        return 0.3;
    case 2:
      if (FunctionM02(clusEnergy, 2.135, -0.245, 0., 0., 0. ) > 0.27)
        return FunctionM02(clusEnergy, 2.135, -0.245, 0., 0., 0. );
      else 
        return 0.27;
    case 3:
      if (FunctionM02(clusEnergy, 2.135, -0.245, 0., 0., 0. ) > 0.25)
        return FunctionM02(clusEnergy, 2.135, -0.245, 0., 0., 0. );
      else 
        return 0.25;
    case 4:
      if (FunctionM02(clusEnergy, 2.135, -0.245, 0.1, 0., 0. ) > 0.27)
        return FunctionM02(clusEnergy, 2.135, -0.245, 0.1, 0., 0. );
      else 
        return 0.27;
    case 5:
      if (FunctionM02(clusEnergy, 2.135, -0.245, -0.1, 0., 0. ) > 0.27)
        return FunctionM02(clusEnergy, 2.135, -0.245, -0.1, 0., 0. );
      else 
        return 0.27;
    case 6: 
      return 0.3;
    case 7:   
      return 0.27;
    case 8: 
      return 0.25;
      
    default:  
      AliError(Form("Min M02 for merged cluster Cut not defined %d",minM02));
      return -1;
  }
  return -1;
}  


//___________________________________________________________________
Bool_t AliCaloPhotonCuts::SetMinM02(Int_t minM02)
{
  fMinM02CutNr = minM02;
  if (fIsPureCalo == 1){
    fUseM02 = 2;
    return kTRUE;
  }

  switch(minM02){
  case 0: 
    if (!fUseM02) fUseM02=0;
    fMinM02=0;
    break;
  case 1: 
    if (!fUseM02) fUseM02=1;
    fMinM02=0.002;
    break;
  case 2:
    if (!fUseM02) fUseM02=1;
    fMinM02=0.1;
    break;
  case 3:
    if (!fUseM02) fUseM02=1;
    fMinM02=0.2;
    break;

  default:
    AliError(Form("Min M02 not defined %d",minM02));
    return kFALSE;
  }
  return kTRUE;
}

//___________________________________________________________________
Bool_t AliCaloPhotonCuts::SetMaxM20(Int_t maxM20)
{
  switch(maxM20){
  case 0: 
    if (!fUseM20) fUseM20=0;
    fMaxM20=100;
    break;
  case 1: 
    if (!fUseM20) fUseM20=1;
    fMaxM20=0.5;
    break;
  default:
    AliError(Form("Max M20 Cut not defined %d",maxM20));
    return kFALSE;
  }
  return kTRUE;
}

//___________________________________________________________________
Bool_t AliCaloPhotonCuts::SetMinM20(Int_t minM20)
{
  switch(minM20){
  case 0: 
    if (!fUseM20) fUseM20=0;
    fMinM20=0;
    break;
  case 1: 
    if (!fUseM20) fUseM20=1;
    fMinM20=0.002;
    break;
  default:
    AliError(Form("Min M20 Cut not defined %d",minM20));
    return kFALSE;
  }
  return kTRUE;
}

//___________________________________________________________________
Bool_t AliCaloPhotonCuts::SetDispersion(Int_t dispersion)
{
  switch(dispersion){
  case 0: 
    if (!fUseDispersion) fUseDispersion=0;
    fMaxDispersion =100;
    break;
  case 1: 
    if (!fUseDispersion) fUseDispersion=1;
    fMaxDispersion=2.;
    break;
  default:
    AliError(Form("Maximum Dispersion Cut not defined %d",dispersion));
    return kFALSE;
  }
  return kTRUE;
}

//___________________________________________________________________
Bool_t AliCaloPhotonCuts::SetNLM(Int_t nlm)
{
  switch(nlm){
  case 0: 
    if (!fUseNLM) fUseNLM=0;
    fMinNLM =0;
    fMaxNLM =100;
    break;
  case 1: 
    if (!fUseNLM) fUseNLM=1;
    fMinNLM =1;
    fMaxNLM =1;
    break;
  case 2: 
    if (!fUseNLM) fUseNLM=1;
    fMinNLM =2;
    fMaxNLM =2;
    break;

  default:
    AliError(Form("NLM Cut not defined %d",nlm));
    return kFALSE;
  }
  return kTRUE;
}

//___________________________________________________________________
Bool_t AliCaloPhotonCuts::SetNonLinearity1(Int_t nl1)
{
  if( nl1 >= 0 && nl1 <=9){
    fNonLinearity1 = nl1;
  }
  else{
    AliError(Form("NonLinearity Correction (part1) not defined %d",nl1));
    return kFALSE;
  }
  return kTRUE;
}

//___________________________________________________________________
Bool_t AliCaloPhotonCuts::SetNonLinearity2(Int_t nl2)
{
  if( nl2 >= 0 && nl2 <=9){
    fNonLinearity2 = nl2;
    if(nl2 == 0) fUseNonLinearity = kFALSE;
    else if(nl2 > 0) fUseNonLinearity = kTRUE;
    fSwitchNonLinearity = fNonLinearity1*10 + fNonLinearity2;
  }
  else{
    AliError(Form("NonLinearity Correction (part2) not defined %d",nl2));
    return kFALSE;
  }
  return kTRUE;
}

//________________________________________________________________________
void AliCaloPhotonCuts::CorrectEMCalNonLinearity(AliVCluster* cluster, Int_t isMC)
{
  if(!fUseNonLinearity) return;

  if (!cluster) {
    AliInfo("Cluster pointer null!");
    return;
  }

  Float_t energy = cluster->E();

  if (energy < 0.05) {
    // Clusters with less than 50 MeV or negative are not possible
    AliInfo(Form("Too Low Cluster energy!, E = %f < 0.05 GeV",energy));
    return;
  }

  if(fCurrentMC==kNoMC){
    AliV0ReaderV1* V0Reader = (AliV0ReaderV1*) AliAnalysisManager::GetAnalysisManager()->GetTask(fV0ReaderName.Data());
    if( V0Reader == NULL ){
      AliFatal(Form("No V0Reader called '%s' could be found within AliCaloPhotonCuts::CorrectEMCalNonLinearity",fV0ReaderName.Data()));
      return;
    }
    fPeriodName = V0Reader->GetPeriodName();
    fCurrentMC = FindEnumForMCSet(fPeriodName);
    
    printf("AliCaloPhotonCuts:Period name has been set to %s, period-enum: %o\n",fPeriodName.Data(),fCurrentMC ) ;
  }
  
    
  Bool_t fPeriodNameAvailable = kTRUE;

  switch(fSwitchNonLinearity){

    // Standard NonLinearity - standard kPi0MCv5 for MC and kSDMv5 for data from Jason
    case 1:
      if( fClusterType == 1|| fClusterType == 3){
        energy *= FunctionNL_kPi0MCv5(energy);
        if(isMC == 0) energy *= FunctionNL_kSDMv5(energy);
      }
      else if ( fClusterType == 2 ){
	// NonLinearity correction from PHOS group for LHC10b-f
        if(isMC>0)
          energy = FunctionNL_PHOS(energy, 1.008, 0.015, 0.4);
      }
      break;

    // kPi0MCv3 for MC and kTestBeamv3 for data
    case 2:
      if(isMC == 0) energy *= FunctionNL_kTestBeamv3(energy);
      else energy *= FunctionNL_kPi0MCv3(energy);
      break;
    // kPi0MCv3 for MC and kTestBeamv2 for data
    case 3:
      if(isMC == 0) energy *= FunctionNL_kTestBeamv2(energy);
      else energy *= FunctionNL_kPi0MCv3(energy);
      break;

    // kPi0MCv2 for MC and kTestBeamv3 for data
    case 4:
      if(isMC == 0) energy *= FunctionNL_kTestBeamv3(energy);
      else energy *= FunctionNL_kPi0MCv2(energy);
      break;
    // kPi0MCv2 for MC and kTestBeamv2 for data
    case 5:
      if(isMC == 0) energy *= FunctionNL_kTestBeamv2(energy);
      else energy *= FunctionNL_kPi0MCv2(energy);
      break;

    // kPi0MCv1 for MC and kTestBeamv3 for data
    case 6:
      if(isMC == 0) energy *= FunctionNL_kTestBeamv3(energy);
      else energy *= FunctionNL_kPi0MCv1(energy);
      break;
    // kPi0MCv1 for MC and kTestBeamv2 for data
    case 7:
      if(isMC == 0) energy *= FunctionNL_kTestBeamv2(energy);
      else energy *= FunctionNL_kPi0MCv1(energy);
      break;

    // kPi0MCv6 for MC and kSDMv6 for data
    case 8:
      if(isMC == 0) energy *= FunctionNL_kSDMv6(energy);
      else energy *= FunctionNL_kPi0MCv6(energy);
      break;

//----------------------------------------------------------------------------------------------------------

// *************** 10 + x **** default tender settings - pp 
      
    // NonLinearity pp ConvCalo - only shifting MC - no timing cut
    case 11:
      label_case_11:
      if(isMC>0){
        // 8TeV LHC12x
        //pass1
        if( fCurrentMC==k14e2a || fCurrentMC==k14e2b )
          energy /= FunctionNL_kSDM(energy, 0.983251, -3.44339, -1.70998);

        else if( fCurrentMC==k14e2c )
          energy /= FunctionNL_kSDM(energy, 0.984462, -3.00363, -2.63773);

        //pass2
        else if( fCurrentMC == k15h1 )
          energy /= FunctionNL_kSDM(energy, 0.96874*0.991*0.9958*0.999, -3.76064, -0.193181);

        else if( fCurrentMC == k15h2 )
          energy /= FunctionNL_kSDM(energy, 0.969703*0.989*0.9969*0.9991, -3.80387, -0.200546);

        else if( fCurrentMC == k16c2 || fCurrentMC == k16c2_plus )
          energy /= FunctionNL_kSDM(energy, 0.974859*0.987*0.996, -3.85842, -0.405277);

        // 2.76TeV LHC11a/LHC13g
        else if( fCurrentMC==k12f1a || fCurrentMC==k12i3 || fCurrentMC==k15g2 )
          energy /= FunctionNL_kSDM(energy, 0.984889*0.995*0.9970, -3.65456, -1.12744);

        else if(fCurrentMC==k12f1b)
          energy /= FunctionNL_kSDM(energy, 0.984384*0.995*0.9970, -3.30287, -1.48516);

        else if( fCurrentMC==k15g1a || fCurrentMC==k15g1b || fCurrentMC==k15a3a || fCurrentMC==k15a3a_plus || fCurrentMC==k15a3b )
          energy /= FunctionNL_kSDM(energy, 0.981892*0.995*0.9970, -5.43438, -1.05468);

        // 7 TeV LHC10x
//        else if( fCurrentMC==k14j4 ) //v1
//          energy /= FunctionNL_kSDM(energy, 0.973866*0.99*0.996*0.999, -4.06436, -0.379);
        else if( fCurrentMC==k14j4 ) //v2
          energy /= FunctionNL_kSDM(energy, 0.974525*0.986*0.999, -4.00247, -0.453046);
        
        // pp 5.02 TeV LHC15n
        else if( fCurrentMC==k16k5a ) {
          if(fClusterType==3) energy /= FunctionNL_kSDM(energy, 0.980211, -4.374598, -0.171988);
          if(fClusterType==1) energy /= FunctionNL_kSDM(energy, 0.984876, -9.999609, -4.999891);
        }
        else if( fCurrentMC==k16k5b ) {
          if(fClusterType==3) energy /= FunctionNL_kSDM(energy, 0.981417, -2.772002, -0.955616);
          if(fClusterType==1) energy /= FunctionNL_kSDM(energy, 0.980275, -3.172374, -0.730326);
        }

        else fPeriodNameAvailable = kFALSE;
      }
      break;

    // NonLinearity pp Calo - only shifting MC - no timing cut
    case 12:
      label_case_12:
      if(isMC>0){
        // 8TeV LHC12x
        //pass1
        if( fCurrentMC==k14e2a || fCurrentMC==k14e2b )
          energy /= FunctionNL_kSDM(2.0*energy, 0.967301, -3.1683, -0.653058);

        else if( fCurrentMC==k14e2c )
          energy /= FunctionNL_kSDM(2.0*energy, 0.96728, -2.96279, -0.903677);

        //pass2
        else if( fCurrentMC == k15h1 )
          energy /= FunctionNL_kSDM(energy, 0.963379*0.9985*0.9992, -3.61217, -0.614043);

        else if( fCurrentMC == k15h2 )
          energy /= FunctionNL_kSDM(energy, 0.96105*0.999*0.9996, -3.62239, -0.556256);

        else if( fCurrentMC == k16c2 || fCurrentMC == k16c2_plus )
          energy /= FunctionNL_kSDM(energy, 0.960596*0.999*0.999, -3.48444, -0.766862);

        // 2.76TeV LHC11a/LHC13g
        else if(  fCurrentMC==k12f1a || fCurrentMC==k12i3 || fCurrentMC==k15g2 )
          energy /= FunctionNL_kSDM(2.0*energy, 0.966151*0.995*0.9981, -2.97974, -0.29463);

        else if( fCurrentMC==k12f1b )
          energy /= FunctionNL_kSDM(2.0*energy, 0.988814*0.995*0.9981, 0.335011, -4.30322);

        else if( fCurrentMC==k15g1a || fCurrentMC==k15g1b || fCurrentMC==k15a3a || fCurrentMC==k15a3a_plus || fCurrentMC==k15a3b )
          energy /= FunctionNL_kSDM(2.0*energy, 0.979994*0.995*0.9981, -3.24431, -0.760205);

        // 7TeV LHC10x
//        else if(  fCurrentMC==k14j4 ) //v1
//          energy /= FunctionNL_kSDM(energy, 0.955095*0.9991, -3.44162, -0.486573);
        else if(  fCurrentMC==k14j4 ) //v2
          energy /= FunctionNL_kSDM(energy, 0.962095*0.9991*0.9993, -3.63967, -0.747825);

        else fPeriodNameAvailable = kFALSE;
      }
      break;

    // NonLinearity ConvCalo - kTestBeamv3 + shifting MC
    case 13:
      energy *= FunctionNL_kTestBeamv3(energy);
      goto label_case_11;// goto previous case for shifting MC
      break;

    // NonLinearity Calo - kTestBeamv3 + shifting MC
    case 14:
      energy *= FunctionNL_kTestBeamv3(energy);
      goto label_case_12;// goto previous case for shifting MC
      break;

    // NonLinearity ConvCalo - kPi0MC + kSDM
    case 15:
      // 8TeV LHC12x
      if ( fCurrentMC==k14e2a || fCurrentMC==k14e2b  || fCurrentMC==k14e2c || fCurrentMC == k15h1 || fCurrentMC == k15h2  || fCurrentMC == k12pp8TeV || fCurrentMC == k16c2 || fCurrentMC == k16c2_plus ){
        energy *= FunctionNL_kPi0MC(energy, 1.0, 0.04979, 1.3, 0.0967998, 219.381, 63.1604, 1.011);
        if(isMC == 0) energy *= FunctionNL_kSDM(energy, 0.9846, -3.319, -2.033);
      
      // 2.76TeV LHC11a/LHC13g
      } else if ( fCurrentMC == k12f1a || fCurrentMC == k12i3 || fCurrentMC == k15g2 || fCurrentMC == k12f1b || 
                  fCurrentMC == k15g1a || fCurrentMC == k15g1b || fCurrentMC == k15a3a || fCurrentMC == k15a3a_plus || fCurrentMC == k15a3b ||
                  fCurrentMC == k11pp2760GeV || fCurrentMC == k13pp2760GeV
                ) {
        energy *= FunctionNL_kPi0MC(energy, 1.0, 0.04123, 1.045, 0.0967998, 219.381, 63.1604, 1.014);
        if(isMC == 0) energy *= FunctionNL_kSDM(energy, 0.9807*0.995*0.9970, -3.377, -0.8535);        
      }
      else fPeriodNameAvailable = kFALSE;
      break;

    // NonLinearity Calo - kPi0MC + kSDM
    case 16:
      // 8TeV LHC12x
      if ( fCurrentMC==k14e2a || fCurrentMC==k14e2b  || fCurrentMC==k14e2c || fCurrentMC == k15h1 || fCurrentMC == k15h2  || fCurrentMC == k12pp8TeV || fCurrentMC == k16c2 || fCurrentMC == k16c2_plus ){
        energy *= FunctionNL_kPi0MC(energy, 1.0, 0.06539, 1.121, 0.0967998, 219.381, 63.1604, 1.011);
        if(isMC == 0) energy *= FunctionNL_kSDM(2.0*energy, 0.9676, -3.216, -0.6828);
      
      // 2.76TeV LHC11a/LHC13g
      } else if ( fCurrentMC == k12f1a || fCurrentMC == k12i3 || fCurrentMC == k15g2 || fCurrentMC == k12f1b || 
                  fCurrentMC == k15g1a || fCurrentMC == k15g1b || fCurrentMC == k15a3a || fCurrentMC == k15a3a_plus || fCurrentMC == k15a3b ||
                  fCurrentMC == k11pp2760GeV || fCurrentMC == k13pp2760GeV
                ) {
        energy *= FunctionNL_kPi0MC(energy, 1.0, 0.06115, 0.9535, 0.0967998, 219.381, 63.1604, 1.013);
        if(isMC == 0) energy *= FunctionNL_kSDM(2.0*energy, 0.9772*0.995*0.9981, -3.256, -0.4449);      
      }
      else fPeriodNameAvailable = kFALSE;
      break;

// *************** 20 + x **** modified tender Settings 1 - pp
    // NonLinearity pp ConvCalo - only shifting MC - no timing cut
    case 21:
      label_case_21:
      if(isMC>0){
        // 2.76TeV LHC11a/LHC13g
        if(  fCurrentMC==k12f1a || fCurrentMC==k12i3 )
          energy /= (FunctionNL_DPOW(energy, 1.0443938253, -0.0691830812, -0.1247555443, 1.1673716264, -0.1853095466, -0.0848801702) - 0.0055);
        else if(fCurrentMC==k15g2)
          energy /= (FunctionNL_DPOW(energy, 1.1716155406, -0.1962930603, -0.0193959829, 1.0336659741, -0.0467778485, -0.4407662248) - 0.0055);
        else if(fCurrentMC==k12f1b)
          energy /= (FunctionNL_DPOW(energy, 1.0166321784, -0.0440799552, -0.2611899222, 1.0636538464, -0.0816662488, -0.2173961316) - 0.007);
        else if( fCurrentMC==k15g1a || fCurrentMC==k15g1b )
          energy /= (FunctionNL_DPOW(energy, 1.1100193881, -0.1389194936, -0.0800000242, 1.1673716264, -0.1853095466, -0.0848801702) - 0.017);
        else if( fCurrentMC==k15a3a || fCurrentMC==k15a3a_plus || fCurrentMC==k15a3b )
          energy /= (FunctionNL_DPOW(energy, 1.0520183153, -0.0806102847, -0.1450415920, 1.0336724056, -0.0467844121, -0.4406992764) - 0.016);
        // 8TeV
        else if( fCurrentMC == k15h1 )
          energy /= (FunctionNL_DPOW(energy, 1.0654169768, -0.0935785719, -0.1137883054, 1.1814766150, -0.1980098061, -0.0854569214) - 0.0138);
        else if( fCurrentMC == k15h2 )
          energy /= (FunctionNL_DPOW(energy, 1.0652493513, -0.0929276101, -0.1113762695, 1.1837801885, -0.1999914832, -0.0854569214) - 0.0145);
        else if( fCurrentMC == k16c2 || fCurrentMC == k16c2_plus )
          energy /= (FunctionNL_DPOW(energy, 1.0489259285, -0.0759079646, -0.1239772934, 1.1835846739, -0.1998987993, -0.0854186691) - 0.014);
        // 7 TeV
//        else if( fCurrentMC == k14j4 ) //v1
//          energy /= (FunctionNL_DPOW(energy, 1.1086453117, -0.1373335557, -0.0800000000, 1.1855482066, -0.1999999504, -0.0830177063) - 0.014);
        else if( fCurrentMC == k14j4 ) //v2
          energy /= (FunctionNL_DPOW(energy, 1.1082846035, -0.1369968318, -0.0800000002, 1.1850179319, -0.1999999950, -0.0863054172) - 0.015);
        else fPeriodNameAvailable = kFALSE;
      }
      break;
      
    // NonLinearity pp Calo - only shifting MC - no timing cut
    case 22:
      label_case_22:
      if(isMC>0){
        // 2.76TeV LHC11a/LHC13g
        if(  fCurrentMC==k12f1a || fCurrentMC==k12i3 )
          energy /= (FunctionNL_DPOW(energy, 0.9980625418, -0.0564782662, -0.5, 1.0383412435, -0.0851830429, -0.4999999996) - 0.00175);
        else if( fCurrentMC==k15g2 )
          energy /= (FunctionNL_DPOW(energy, 1.0795372569, -0.1347324732, -0.1630736190, 1.1614181498, -0.199995361, -0.1711378093) - 0.0035);
        else if( fCurrentMC==k12f1b )
          energy /= (FunctionNL_DPOW(energy, 1.0232969083, -0.090409434, -0.3592406513, 1.0383412435, -0.0851830429, -0.4999999996) + 0.0007);
        else if( fCurrentMC==k15g1a || fCurrentMC==k15g1b )
          energy /= (FunctionNL_DPOW(energy, 1.0106037132, -0.0748250591, -0.4999999996, 1.0383412435, -0.0851830429, -0.4999999996) - 0.014);
        else if( fCurrentMC==k15a3a || fCurrentMC==k15a3a_plus || fCurrentMC==k15a3b ) 
          energy /= (FunctionNL_DPOW(energy, 1.0119417393, -0.0755250741, -0.4999999996, 1.1614181498, -0.1999995361, -0.1711378093) - 0.006);
        //8TeV
        else if( fCurrentMC == k15h1 )
          energy /= (FunctionNL_DPOW(energy, 1.1389201636, -0.1999994717, -0.1622237979, 1.1603460704, -0.1999999989, -0.2194447313) - 0.0025);
        else if( fCurrentMC == k15h2 )
          energy /= (FunctionNL_DPOW(energy, 1.0105301622, -0.0732424689, -0.5000000000, 1.0689250170, -0.1082682369, -0.4388156470) - 0.001);
        else if( fCurrentMC == k16c2 || fCurrentMC == k16c2_plus )
          energy /= (FunctionNL_DPOW(energy, 0.9922456908, -0.0551212559, -0.5000000000, 1.0513459039, -0.0894163252, -0.5000000000) + 0.002);
        // 7 TeV
//        else if( fCurrentMC == k14j4 ) //v1
//          energy /= (FunctionNL_DPOW(energy, 0.9994789138, -0.0601419399, -0.4999999999, 1.1635744933, -0.1999999978, -0.1985578372) - 0.005);
        else if( fCurrentMC == k14j4 ) //v2
          energy /= (FunctionNL_DPOW(energy, 1.0074002842, -0.0682543971, -0.4509341085, 1.1224162203, -0.1586806096, -0.2458351112) - 0.003);
        else fPeriodNameAvailable = kFALSE;
      }    
      break;
    // NonLinearity ConvCalo - kTestBeamv3 + shifting MC  
    case 23:
      energy *= FunctionNL_kTestBeamv3(energy);
      goto label_case_21;// goto previous case for shifting MC
      break;

    // NonLinearity Calo - kTestBeamv3 + shifting MC
    case 24:
      energy *= FunctionNL_kTestBeamv3(energy);
      goto label_case_22;// goto previous case for shifting MC
      break;

      
// *************** 30 + x **** modified tender Settings 2 - pp
    case 31:
      label_case_31:
      if(isMC>0){
        // 2.76TeV LHC11a/LHC13g
        if(  fCurrentMC==k12f1a || fCurrentMC==k12i3 )
          energy /= FunctionNL_kSDM(energy, 0.983176*0.9945, -3.91107, -0.697613);
        else if(fCurrentMC==k15g2)
          energy /= FunctionNL_kSDM(energy, 0.972574*0.9942, -3.19191, -0.946239);
        else if(fCurrentMC==k12f1b)
          energy /= FunctionNL_kSDM(energy, 0.981893*0.9930, -4.05476, -0.710661);
        else if( fCurrentMC==k15g1a || fCurrentMC==k15g1b )
          energy /= FunctionNL_kSDM(energy, 0.983176*0.993*0.99, -1.85546, -3.37696);
        else if( fCurrentMC==k15a3a || fCurrentMC==k15a3a_plus || fCurrentMC==k15a3b )
          energy /= FunctionNL_kSDM(energy, 0.977035*0.9835, -3.82187, -1.04332);               
        else fPeriodNameAvailable = kFALSE;
      }    
      break;
    // NonLinearity pp Calo - only shifting MC - no timing cut
    case 32:
      label_case_32:
      if(isMC>0){
        // 2.76TeV LHC11a/LHC13g
        if(  fCurrentMC==k12f1a || fCurrentMC==k12i3 )
          energy /= FunctionNL_kSDM(energy, 0.974358*0.9987, -2.18037, -1.91622);
        else if( fCurrentMC==k15g2 )
          energy /= FunctionNL_kSDM(energy, 0.963307*0.9962, -3.27998, -0.589806);
        else if( fCurrentMC==k12f1b )
          energy /= FunctionNL_kSDM(energy, 0.97499*0.9995, -0.180148, -4.78066);
        else if( fCurrentMC==k15g1a || fCurrentMC==k15g1b )
          energy /= FunctionNL_kSDM(energy, 0.974424*0.998*0.992, -0.533785, -4.06374);
        else if ( fCurrentMC==k15a3a || fCurrentMC==k15a3a_plus || fCurrentMC==k15a3b ) 
          energy /= FunctionNL_kSDM(energy, 0.963307*0.995, -4.01949, -0.38667);    
        else fPeriodNameAvailable = kFALSE;
      }    
      break;
    // NonLinearity ConvCalo - kTestBeamv3 + shifting MC  
    case 33:
      energy *= FunctionNL_kTestBeamv3(energy);
      goto label_case_31;// goto previous case for shifting MC
      break;

    // NonLinearity Calo - kTestBeamv3 + shifting MC
    case 34:
      energy *= FunctionNL_kTestBeamv3(energy);
      goto label_case_32;// goto previous case for shifting MC
      break;

      
// *************** 40 + x **** default tender Settings - pPb
    // NonLinearity LHC13 pPb ConvCalo  - only shifting MC
    case 41:
      label_case_41:
      if(isMC>0){
        if( fCurrentMC==k13b2_efix || fCurrentMC == k16c3a || fCurrentMC == k16c3b || fCurrentMC == k16c3c ){
          energy /= FunctionNL_kSDM(energy, 0.995*0.978578, -3.80517, -0.581197);//v4
        }
        else if( fCurrentMC==k13e7 ) energy /= FunctionNL_kSDM(energy, 0.979813, -3.53445, -0.733067);//v0

        else fPeriodNameAvailable = kFALSE;
      }
      break;

    // NonLinearity LHC13 pPb Calo  - only shifting MC
    case 42:
      label_case_42:
      if(isMC>0){
        if( fCurrentMC==k13b2_efix || fCurrentMC == k16c3a || fCurrentMC == k16c3b || fCurrentMC == k16c3c ) 
          energy /= FunctionNL_kSDM(energy, 1.002*0.970383, -3.65936, -0.721139);//v4
        else if( fCurrentMC==k13e7 ) 
          energy /= FunctionNL_kSDM(energy, 0.970537, -3.36675, -0.958747);//v0
        else fPeriodNameAvailable = kFALSE;
      }
      break;

    // NonLinearity LHC13 pPb ConvCalo  - kTestBeamv3 + shifting MC
    case 43:
      energy *= FunctionNL_kTestBeamv3(energy);
      goto label_case_41;// goto previous case for shifting MC
      break;

    // NonLinearity LHC13 pPb Calo  - kTestBeamv3 + shifting MC
    case 44:
      energy *= FunctionNL_kTestBeamv3(energy);
      goto label_case_42;// goto previous case for shifting MC
      break;

    // NonLinearity LHC13 pPb Calo - excluding the two lowest pT points
    case 49:
      if(isMC>0){
        if( fCurrentMC==k13b2_efix || fCurrentMC==k13e7 || fCurrentMC == k16c3a || fCurrentMC == k16c3b || fCurrentMC == k16c3c ) 
          energy /= FunctionNL_kSDM(energy, 0.973302, -3.12524, -1.13546);
        else fPeriodNameAvailable = kFALSE;
      }
      break;

// *************** 50 + x **** modified tender Settings 1 - pPb
    // NonLinearity LHC13 pPb ConvCalo  - only shifting MC
    case 51:
      label_case_51:
      if(isMC>0){
        if( fCurrentMC==k13b2_efix || fCurrentMC == k16c3a || fCurrentMC == k16c3b || fCurrentMC == k16c3c ) 
          energy /= (FunctionNL_DPOW(energy, 1.0754004911, -0.0992327361, -0.0802161499, 1.1849304274, -0.1999999986, -0.0828138864) - 0.005);//v4
        else if( fCurrentMC==k13e7 ) 
          energy /= FunctionNL_DPOW(energy, 1.0546114304, -0.0758513555, -0.0800000002, 1.1849400584, -0.1999999970, -0.0826417756);//v4
        else fPeriodNameAvailable = kFALSE;
      }
      break;

    // NonLinearity LHC13 pPb Calo  - only shifting MC
    case 52:
      label_case_52:
      if(isMC>0){
        if( fCurrentMC==k13b2_efix || fCurrentMC == k16c3a || fCurrentMC == k16c3b || fCurrentMC == k16c3c ) 
          energy /= (FunctionNL_DPOW(energy, 1.0358207569, -0.0914347267, -0.3683743201, 1.1549558754, -0.1942615277, -0.2216281109) + 0.002);//v4
        else if( fCurrentMC==k13e7 ) 
          energy /= FunctionNL_DPOW(energy, 1.0149551972, -0.0697288693, -0.4586527438, 1.1549558754, -0.1942615277, -0.2216281109);//v4
        else fPeriodNameAvailable = kFALSE;
      }
      break;

      
// *************** 60 + x **** modified tender Settings 2 - pPb


// *************** 70 + x **** default tender Settings - PbPb
      

// *************** 80 + x **** modified tender Settings 1 - PbPb
      
      
// *************** 90 + x **** modified tender Settings 2 - PbPb
      
      
      
//----------------------------------------------------------------------------------------------------------

    default:
      AliFatal(Form("NonLinearity correction not defined for cut: '%d' ! Returning...",fSwitchNonLinearity));
      return;

  }

  if(!fPeriodNameAvailable){
    AliFatal(Form("NonLinearity correction not defined for fPeriodName: '%s'! Please check cut number (%d) as well as function AliCaloPhotonCuts::CorrectEMCalNonLinearity. Correction failed, returning...",fPeriodName.Data(),fSwitchNonLinearity));
    return;
  }

  cluster->SetE(energy);

  return;
}

//________________________________________________________________________
Float_t AliCaloPhotonCuts::FunctionNL_kPi0MC(Float_t e, Float_t p0, Float_t p1, Float_t p2, Float_t p3, Float_t p4, Float_t p5, Float_t p6){
  return ( p6 / ( p0 * ( 1. / ( 1. + p1 * exp( -e / p2 ) ) * 1. / ( 1. + p3 * exp( ( e - p4 ) / p5 ) ) ) ) );
}

//________________________________________________________________________
Float_t AliCaloPhotonCuts::FunctionNL_kSDM(Float_t e, Float_t p0, Float_t p1, Float_t p2){
  return ( p0 + exp( p1 + ( p2 * e ) ) );
}

//________________________________________________________________________
Float_t AliCaloPhotonCuts::FunctionNL_DPOW(Float_t e, Float_t p0, Float_t p1, Float_t p2, Float_t p3, Float_t p4, Float_t p5){
  Float_t ret = 1;
  if ((p3 +  p4 * TMath::Power(e,p5 ) ) != 0)
    ret = ( (p0 +  p1 * TMath::Power(e,p2 ) )/(p3 +  p4 * TMath::Power(e,p5 ) ) );
  if (ret != 0.)
    return ret;
  else 
    return 1.;
}

//________________________________________________________________________
Float_t AliCaloPhotonCuts::FunctionNL_PHOS(Float_t e, Float_t p0, Float_t p1, Float_t p2){
  return (0.0241+1.0504*e+0.000249*e*e)*p0*(1+p1/(1.+e*e/p2/p2)) ;
}

//************************************************************************
// predefined functions:
//________________________________________________________________________
Float_t AliCaloPhotonCuts::FunctionNL_kPi0MCv1(Float_t e){
  return ( 1.014 * exp( 0.03329 / e ) ) + ( ( -0.3853 / ( 0.5423 * 2. * TMath::Pi() ) * exp( -( e + 0.4335 ) * ( e + 0.4335 ) / (2. * 0.5423 * 0.5423 ) ) ) );
}

//________________________________________________________________________
Float_t AliCaloPhotonCuts::FunctionNL_kPi0MCv2(Float_t e){
  return ( 0.311111 / TMath::Power( e - 0.571666, 0.567995 ) + 1 );
}

//________________________________________________________________________
Float_t AliCaloPhotonCuts::FunctionNL_kPi0MCv3(Float_t e){
  return ( 1.0 / ( 0.981039 * ( 1. / ( 1. + 0.113508 * exp( -e / 1.00173 ) ) * 1. / ( 1. + 0.0967998 * exp( ( e - 219.381 ) / 63.1604 ) ) ) ) );
}

//________________________________________________________________________
Float_t AliCaloPhotonCuts::FunctionNL_kPi0MCv5(Float_t e){
  return ( 1.01286 / ( 1.0 * ( 1. / ( 1. + 0.0664778 * exp( -e / 1.57 ) ) * 1. / ( 1. + 0.0967998 * exp( ( e - 219.381 ) / 63.1604 ) ) ) ) );
}

//________________________________________________________________________
Float_t AliCaloPhotonCuts::FunctionNL_kPi0MCv6(Float_t e){
  return ( 1.00437 / ( 1.0 * ( 1. / ( 1. + 0.0797873 * exp( -e / 1.68322 ) ) * 1. / ( 1. + 0.0806098 * exp( ( e - 244.586 ) / 116.938 ) ) ) ) );
}

// only shifting data, to be used with kPi0MCv5 before
//________________________________________________________________________
Float_t AliCaloPhotonCuts::FunctionNL_kSDMv5(Float_t e){
  return ( 0.964 + exp( -3.132 + ( -0.435 * 2.0 * e ) ) );
}

// be careful: different definition than kSDMv5
//________________________________________________________________________
Float_t AliCaloPhotonCuts::FunctionNL_kSDMv6(Float_t e){
  return ( 0.987054 / ( 1.0 * ( 1. / ( 1. + 0.237767 * exp( -e / 0.651203 ) ) * 1. / ( 1. + 0.183741 * exp( ( e - 155.427 ) / 17.0335 ) ) ) ) );
}

//________________________________________________________________________
Float_t AliCaloPhotonCuts::FunctionNL_kTestBeamv2(Float_t e){
  return ( 0.968 / ( 0.983504 *( 1. / ( 1. + 0.210106 * exp( -e / 0.897274 ) ) * 1. / ( 1. + 0.0829064 * exp( ( e - 152.299 ) / 31.5028 ) ) ) ) );
}

//________________________________________________________________________
Float_t AliCaloPhotonCuts::FunctionNL_kTestBeamv3(Float_t e){
  return ( 0.9615 / ( 0.976941 *( 1. / ( 1. + 0.162310 * exp( -e / 1.08689 ) ) * 1. / ( 1. + 0.0819592 * exp( ( e - 152.338 ) / 30.9594 ) ) ) ) );
}

//************************************************************************
//________________________________________________________________________
Float_t AliCaloPhotonCuts::FunctionM02(Float_t E, Float_t a, Float_t b, Float_t c, Float_t d, Float_t e){
  return ( exp( a+ b*E ) + c + d*E + e/E);
}

//________________________________________________________________________
AliCaloPhotonCuts::MCSet AliCaloPhotonCuts::FindEnumForMCSet(TString namePeriod){
  if(       namePeriod.CompareTo("LHC14e2a")==0)        return k14e2a;
  else if(  namePeriod.CompareTo("LHC14e2b")==0)        return k14e2b;
  else if(  namePeriod.CompareTo("LHC14e2c")==0)        return k14e2c;
  else if(  namePeriod.CompareTo("LHC12f1a")==0)        return k12f1a;
  else if(  namePeriod.CompareTo("LHC12f1b")==0)        return k12f1b;
  else if(  namePeriod.CompareTo("LHC12i3")==0)         return k12i3;
  else if(  namePeriod.CompareTo("LHC15g1a")==0)        return k15g1a;
  else if(  namePeriod.CompareTo("LHC15g1b")==0)        return k15g1b;
  else if(  namePeriod.CompareTo("LHC15g2")==0)         return k15g2;
  else if(  namePeriod.CompareTo("LHC15a3a")==0)        return k15a3a;
  else if(  namePeriod.CompareTo("LHC15a3a_plus")==0)   return k15a3a_plus;
  else if(  namePeriod.CompareTo("LHC15a3b")==0)        return k15a3b;
  else if(  namePeriod.Contains("LHC13b2_efix"))        return k13b2_efix;
  else if(  namePeriod.Contains("LHC13e7"))             return k13e7;
  else if(  namePeriod.CompareTo("LHC15h1a1")==0 ||
            namePeriod.CompareTo("LHC15h1b")==0 ||
            namePeriod.CompareTo("LHC15h1c")==0 ||
            namePeriod.CompareTo("LHC15h1d")==0 ||
            namePeriod.CompareTo("LHC15h1f")==0 ||
            namePeriod.CompareTo("LHC15h1g")==0 ||
            namePeriod.CompareTo("LHC15h1h")==0 ||
            namePeriod.CompareTo("LHC15h1i")==0)        return k15h1;
  else if(  namePeriod.CompareTo("LHC15h2a")==0 ||
            namePeriod.CompareTo("LHC15h2b")==0 ||
            namePeriod.CompareTo("LHC15h2c")==0 ||
            namePeriod.CompareTo("LHC15h2d")==0 ||
            namePeriod.CompareTo("LHC15h2f")==0 ||
            namePeriod.CompareTo("LHC15h2g")==0 ||
            namePeriod.CompareTo("LHC15h2h")==0 ||
            namePeriod.CompareTo("LHC15h2i")==0)        return k15h2;
  else if(  namePeriod.CompareTo("LHC14j4b")==0 ||
            namePeriod.CompareTo("LHC14j4c")==0 ||
            namePeriod.CompareTo("LHC14j4d")==0 ||
            namePeriod.CompareTo("LHC14j4e")==0 ||
            namePeriod.CompareTo("LHC14j4f")==0)        return k14j4;
  else if ( namePeriod.CompareTo("LHC16c2") == 0 )      return k16c2;
  else if ( namePeriod.CompareTo("LHC16c2_plus") == 0 ) return k16c2_plus;
  else if ( namePeriod.CompareTo("LHC16c3a") == 0 )     return k16c3a;
  else if ( namePeriod.CompareTo("LHC16c3b") == 0 )     return k16c3b;
  else if ( namePeriod.CompareTo("LHC16c3c") == 0 )     return k16c3c;
  else if ( namePeriod.CompareTo("LHC16h3") == 0 )      return k16h3;
  else if ( namePeriod.CompareTo("LHC16h3b") == 0 )     return k16h3b;
  else if ( namePeriod.CompareTo("LHC16h8a") == 0 )     return k16h8a;
  else if ( namePeriod.CompareTo("LHC16h8b") == 0 )     return k16h8b;
  else if ( namePeriod.CompareTo("LHC16k3a") == 0 || 
            namePeriod.CompareTo("LHC16k3a2") == 0 )     return k16k3a;
  else if ( namePeriod.CompareTo("LHC16k3b") == 0 || 
            namePeriod.CompareTo("LHC16k3b2") == 0 )     return k16k3b;
  else if ( namePeriod.CompareTo("LHC16k5a") == 0  )     return k16k5a;
  else if ( namePeriod.CompareTo("LHC16k5b") == 0  )     return k16k5b;
  else if ( namePeriod.CompareTo("LHC10b") == 0 ||
            namePeriod.CompareTo("LHC10c") == 0 ||
            namePeriod.CompareTo("LHC10d") == 0 ||
            namePeriod.CompareTo("LHC10e") == 0 ||
            namePeriod.CompareTo("LHC10f") == 0 ||
            namePeriod.CompareTo("LHC10g") == 0 )       return k10pp7TeV;
  else if ( namePeriod.CompareTo("LHC10h") == 0 )       return k10PbPb2760GeV;
  else if ( namePeriod.CompareTo("LHC11a") == 0 )       return k11pp2760GeV;
  else if ( namePeriod.CompareTo("LHC11b") == 0 ||
            namePeriod.CompareTo("LHC11c") == 0 ||
            namePeriod.CompareTo("LHC11d") == 0 ||
            namePeriod.CompareTo("LHC11e") == 0 ||
            namePeriod.CompareTo("LHC11f") == 0 ||
            namePeriod.CompareTo("LHC11g") == 0 )       return k11pp7TeV;
  else if ( namePeriod.CompareTo("LHC11h") == 0 )       return k11PbPb2760GeV;
  else if ( namePeriod.CompareTo("LHC12a") == 0 ||
            namePeriod.CompareTo("LHC12b") == 0 ||
            namePeriod.CompareTo("LHC12c") == 0 ||
            namePeriod.CompareTo("LHC12d") == 0 ||
            namePeriod.CompareTo("LHC12e") == 0 ||
            namePeriod.CompareTo("LHC12f") == 0 ||
            namePeriod.CompareTo("LHC12g") == 0 ||
            namePeriod.CompareTo("LHC12h") == 0 ||
            namePeriod.CompareTo("LHC12i") == 0 )       return k12pp8TeV;
  else if ( namePeriod.CompareTo("LHC13b") == 0 ||
            namePeriod.CompareTo("LHC13c") == 0 ||
            namePeriod.CompareTo("LHC13d") == 0 ||
            namePeriod.CompareTo("LHC13e") == 0 ||
            namePeriod.CompareTo("LHC13f") == 0 )       return k13pPb5023GeV;
  else if ( namePeriod.CompareTo("LHC13g") == 0 )       return k13pp2760GeV;
  else if ( namePeriod.CompareTo("LHC15f") == 0 ||
            namePeriod.CompareTo("LHC15g") == 0 ||
            namePeriod.CompareTo("LHC15h") == 0 ||
            namePeriod.CompareTo("LHC15i") == 0 ||
            namePeriod.CompareTo("LHC15j") == 0 ||
            namePeriod.CompareTo("LHC15k") == 0 ||
            namePeriod.CompareTo("LHC15l") == 0 ||
            namePeriod.CompareTo("LHC15m") == 0 )       return k15pp13TeV;
  else if ( namePeriod.CompareTo("LHC15n") == 0 )       return k15pp5TeV;
  else if ( namePeriod.CompareTo("LHC15o") == 0 )       return k15PbPb5TeV;
  else return kNoMC;
}

//________________________________________________________________________
void AliCaloPhotonCuts::SetLogBinningXTH1(TH1* histoRebin){
  TAxis *axisafter  = histoRebin->GetXaxis();
  Int_t bins        = axisafter->GetNbins();
  Double_t from     = axisafter->GetXmin();
  Double_t to       = axisafter->GetXmax();
  Double_t *newbins = new Double_t[bins+1];
  newbins[0]        = from;
  Double_t factor   = TMath::Power(to/from, 1./bins);
  for(Int_t i=1; i<=bins; ++i) newbins[i] = factor * newbins[i-1];
  axisafter->Set(bins, newbins);
  delete [] newbins;
  return;
}


//________________________________________________________________________
void AliCaloPhotonCuts::SetLogBinningXTH2(TH2* histoRebin){
  TAxis *axisafter  = histoRebin->GetXaxis();
  Int_t bins        = axisafter->GetNbins();
  Double_t from     = axisafter->GetXmin();
  Double_t to       = axisafter->GetXmax();
  Double_t *newbins = new Double_t[bins+1];
  newbins[0]        = from;
  Double_t factor   = TMath::Power(to/from, 1./bins);
  for(Int_t i=1; i<=bins; ++i) newbins[i] = factor * newbins[i-1];
  axisafter->Set(bins, newbins);
  delete [] newbins;
  return;
}

//________________________________________________________________________
void AliCaloPhotonCuts::SetLogBinningYTH2(TH2* histoRebin){
  TAxis *axisafter  = histoRebin->GetYaxis();
  Int_t bins        = axisafter->GetNbins();
  Double_t from     = axisafter->GetXmin();
  Double_t to       = axisafter->GetXmax();
  Double_t *newbins = new Double_t[bins+1];
  newbins[0]        = from;
  Double_t factor   = TMath::Power(to/from, 1./bins);
  for(Int_t i=1; i<=bins; ++i) newbins[i] = factor * newbins[i-1];
  axisafter->Set(bins, newbins);
  delete [] newbins;
  return;
}

//________________________________________________________________________
Double_t AliCaloPhotonCuts::GetDistanceBetweenClusters(AliVCluster* cluster1, AliVCluster* cluster2){
  
  Float_t clusPos1[3]   = {0,0,0};
  cluster1->GetPosition(clusPos1);
  TVector3 clusterVector1(clusPos1[0],clusPos1[1],clusPos1[2]);
  Double_t etaCluster1  = clusterVector1.Eta();
  Double_t phiCluster1  = clusterVector1.Phi();
  if (phiCluster1 < 0) phiCluster1 += 2*TMath::Pi();
 
  Float_t clusPos2[3]   = {0,0,0};
  cluster2->GetPosition(clusPos2);
  TVector3 clusterVector2(clusPos2[0],clusPos2[1],clusPos2[2]);
  Double_t etaCluster2  = clusterVector2.Eta();
  Double_t phiCluster2  = clusterVector2.Phi();
  if (phiCluster2 < 0) phiCluster2 += 2*TMath::Pi();
  
  Double_t deltaEta     = TMath::Abs(etaCluster1-etaCluster2);
  Double_t deltaPhi     = 0;
  if (phiCluster1 > phiCluster2){
    deltaPhi            = phiCluster1-phiCluster2;
    if (deltaPhi > TMath::Pi()) 
      deltaPhi          = 2*TMath::Pi()-deltaPhi;
  } else {
    deltaPhi            = phiCluster2-phiCluster1;
    if (deltaPhi > TMath::Pi()) 
      deltaPhi          = 2*TMath::Pi()-deltaPhi;
  }  
  Double_t r            = TMath::Sqrt(deltaPhi*deltaPhi+deltaEta*deltaEta);
  return r;
}  

//________________________________________________________________________
TString AliCaloPhotonCuts::GetCutNumber(){
   // returns TString with current cut number
   TString a(kNCuts);
   for(Int_t ii=0;ii<kNCuts;ii++){
      a.Append(Form("%d",fCuts[ii]));
   }
   return a;
}


//___________________________________________________________________
// Check if the cluster highest energy tower is exotic.
//___________________________________________________________________
Bool_t AliCaloPhotonCuts::IsExoticCluster( AliVCluster *cluster, AliVEvent *event, Float_t &energyStar ) {  
  
  if (!cluster) {
    AliInfo("Cluster pointer null!");
    return kFALSE;
  }
  energyStar              = 0; 
  
  AliVCaloCells* cells    = NULL;
  if (fClusterType == 1 || fClusterType == 3) 
    cells                 = event->GetEMCALCells();
  else if (fClusterType ==2 ) 
    return kFALSE;
//     cells                 = event->GetPHOSCells();


  Int_t largestCellID     = FindLargestCellInCluster(cluster,event);  
  Float_t ecell1          = cells->GetCellAmplitude(largestCellID); ; 
  Float_t eCross          = GetECross(largestCellID,cells);
  energyStar              = ecell1+eCross;
  
  if (ecell1 < fExoticMinEnergyCell)
    return kFALSE;
  
  if (1-eCross/ecell1 > fExoticEnergyFracCluster) {
    return kTRUE;
  }  else { 
    return kFALSE;
  }  
  return kFALSE;
}

//___________________________________________________________________________
// Calculate the energy in the cross around the energy of a given cell.
// Used in exotic clusters/cells rejection.
//___________________________________________________________________________
Float_t AliCaloPhotonCuts::GetECross( Int_t absID, AliVCaloCells* cells ){  
  
  Int_t imod = -1, iphi =-1, ieta=-1,iTower = -1, iIphi = -1, iIeta = -1; 
  fGeomEMCAL->GetCellIndex(absID,imod,iTower,iIphi,iIeta); 
  fGeomEMCAL->GetCellPhiEtaIndexInSModule(imod,iTower,iIphi, iIeta,iphi,ieta);  
  
  // Get close cells index, energy and time, not in corners  
  Int_t absID1 = -1;
  Int_t absID2 = -1;
  
  if ( iphi < AliEMCALGeoParams::fgkEMCALRows-1) absID1 = fGeomEMCAL->GetAbsCellIdFromCellIndexes(imod, iphi+1, ieta);
  if ( iphi > 0 )                                absID2 = fGeomEMCAL->GetAbsCellIdFromCellIndexes(imod, iphi-1, ieta);
  
  // In case of cell in eta = 0 border, depending on SM shift the cross cell index
  Int_t absID3 = -1;
  Int_t absID4 = -1;
  if ( fClusterType == 1 && ieta == AliEMCALGeoParams::fgkEMCALCols-1 && !(imod%2) ){
    absID3 = fGeomEMCAL-> GetAbsCellIdFromCellIndexes(imod+1, iphi, 0);
    absID4 = fGeomEMCAL-> GetAbsCellIdFromCellIndexes(imod,   iphi, ieta-1); 
  } else if ( fClusterType == 1 && ieta == 0 && imod%2 ) {
    absID3 = fGeomEMCAL-> GetAbsCellIdFromCellIndexes(imod,   iphi, ieta+1);
    absID4 = fGeomEMCAL-> GetAbsCellIdFromCellIndexes(imod-1, iphi, AliEMCALGeoParams::fgkEMCALCols-1); 
  } else {
    if ( ieta < AliEMCALGeoParams::fgkEMCALCols-1 ) 
      absID3 = fGeomEMCAL-> GetAbsCellIdFromCellIndexes(imod, iphi, ieta+1);
    if ( ieta > 0 )                                 
      absID4 = fGeomEMCAL-> GetAbsCellIdFromCellIndexes(imod, iphi, ieta-1); 
  }
    
  Float_t  ecell1  = 0, ecell2  = 0, ecell3  = 0, ecell4  = 0;

  
  // Do not include bad channels found in analysis,
  if (AcceptCellByBadChannelMap(absID1)) ecell1 = cells->GetCellAmplitude(absID1); ; 
  if (AcceptCellByBadChannelMap(absID2)) ecell2 = cells->GetCellAmplitude(absID2); ; 
  if (AcceptCellByBadChannelMap(absID3)) ecell3 = cells->GetCellAmplitude(absID3); ; 
  if (AcceptCellByBadChannelMap(absID4)) ecell4 = cells->GetCellAmplitude(absID4); ; 
  
  return ecell1+ecell2+ecell3+ecell4;
}

//_______________________________________________________________________________
// Determine whether cell would be accepted according to bad channel map
//_______________________________________________________________________________
Bool_t AliCaloPhotonCuts::AcceptCellByBadChannelMap(Int_t absID ){
  if (fClusterType == 1 || fClusterType == 3){
    if(!fGeomEMCAL) {
      AliFatal("No instance of the geometry is available");
      return kFALSE;
    }
    if ( absID < 0 || absID >= 24*48*fGeomEMCAL->GetNumberOfSuperModules() ) 
      return kFALSE;

    Int_t imod = -1, iphi =-1, ieta=-1,iTower = -1, iIphi = -1, iIeta = -1; 
    fGeomEMCAL->GetCellIndex(absID,imod,iTower,iIphi,iIeta);
    fGeomEMCAL->GetCellPhiEtaIndexInSModule(imod,iTower,iIphi, iIeta,iphi,ieta);  

    // Do not include bad channels found in analysis,
    if (fEMCALRecUtils->GetEMCALChannelStatus(imod, ieta, iphi) == 0 ) 
      return kTRUE;
    else 
      return kFALSE;
  } 
  return kTRUE;
}

//_______________________________________________________________________________
// Classify clusters for TM effi
// 0: Neutral cluster
// 1: Neutral cluster sub charged
// 2: Gamma cluster
// 3: Gamma cluster sub charged
// 4: Gamma conv cluster
// 5: Charged cluster
// 6: electron
// 7: primary charged
//_______________________________________________________________________________
Int_t AliCaloPhotonCuts::ClassifyClusterForTMEffi(AliVCluster* cluster, AliVEvent* event, AliVEvent* mcEvent, Bool_t isESD){
  Int_t* mclabelsCluster  = cluster->GetLabels();
  Int_t classification    = -1;
  
  const AliVVertex* primVtxMC   = mcEvent->GetPrimaryVertex();
  Double_t mcProdVtxX   = primVtxMC->GetX();
  Double_t mcProdVtxY   = primVtxMC->GetY();
  Double_t mcProdVtxZ   = primVtxMC->GetZ();

  if (isESD){
    AliStack* fMCStack    = NULL;
    if(mcEvent){
      fMCStack          = dynamic_cast<AliStack*>(dynamic_cast<AliMCEvent*>(mcEvent)->Stack());
    }  
    if (fMCStack == NULL){
      AliError("No stack available in ESD");
      return -1;
    }  
    
    if (cluster->GetNLabels()>0){
      TParticle* particleLead   = (TParticle*)fMCStack->Particle(mclabelsCluster[0]);
      Double_t charge           = ((TParticlePDG*)particleLead->GetPDG())->Charge();
      if (charge == 0 || charge == -0){
        classification        = 0;
        if (particleLead->GetPdgCode() == 22) 
          classification      = 2;
      } else {
        classification        = 5;
        if (particleLead->GetPdgCode() == 11 || particleLead->GetPdgCode() == -11){
          classification      = 6;
          if (particleLead->GetMother(0) > -1)
            if (((TParticle*)fMCStack->Particle(particleLead->GetMother(0)))->GetPdgCode() == 22)
              classification  = 4;
        } else {
            Double_t deltaX = particleLead->Vx() - mcProdVtxX;
            Double_t deltaY = particleLead->Vy() - mcProdVtxY;
            Double_t deltaZ = particleLead->Vz() - mcProdVtxZ;

            Double_t realRadius3D = TMath::Sqrt(deltaX*deltaX+deltaY*deltaY+deltaZ*deltaZ);
            if (realRadius3D < 5.0)
              classification  = 7;
        }  
      }
      if ((classification == 0 || classification == 2) && cluster->GetNLabels() > 1){
        Bool_t goOut = kFALSE;
        for (Int_t i = 1; (i< cluster->GetNLabels() && !goOut); i++){
          TParticle* particleSub    = (TParticle*)fMCStack->Particle(mclabelsCluster[i]);
          Double_t charge           = ((TParticlePDG*)particleSub->GetPDG())->Charge();
          if (!(charge == 0 || charge == -0)){
            classification++;
            goOut = kTRUE;
          }  
        }
      }
    }
  } else {
    TClonesArray *AODMCTrackArray = dynamic_cast<TClonesArray*>(event->FindListObject(AliAODMCParticle::StdBranchName()));
    if (AODMCTrackArray == NULL){
      AliError("No MC particle list available in AOD");
      return -1;
    }
    
    if (cluster->GetNLabels()>0){
      AliAODMCParticle* particleLead    = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(mclabelsCluster[0]));
      Short_t charge                    = particleLead->Charge();
      if (charge == 0 ){
        classification        = 0;
        if (particleLead->GetPdgCode() == 22) 
          classification      = 2;
      } else {
        classification        = 5;
        if (particleLead->GetPdgCode() == 11 || particleLead->GetPdgCode() == -11) {
          classification      = 6;
          if (particleLead->GetMother() > -1){
            AliAODMCParticle* mother    = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(particleLead->GetMother()));
            if (mother->GetPdgCode() == 22)
              classification  = 4;
          }
        } else {
          Double_t deltaX = particleLead->Xv() - mcProdVtxX;
          Double_t deltaY = particleLead->Yv() - mcProdVtxY;
          Double_t deltaZ = particleLead->Zv() - mcProdVtxZ;

          Double_t realRadius3D = TMath::Sqrt(deltaX*deltaX+deltaY*deltaY+deltaZ*deltaZ);
          if (realRadius3D < 5.0)
            classification    = 7;
        }  
      }
      if ((classification == 0 || classification == 2) && cluster->GetNLabels() > 1){
        Bool_t goOut = kFALSE;
        for (Int_t i = 1; (i< cluster->GetNLabels() && !goOut); i++){
          AliAODMCParticle* particleSub = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(mclabelsCluster[i]));
          Double_t charge               = particleSub->Charge();
          if ( charge != 0 ){
            classification++;
            goOut = kTRUE;
          }  
        }
      }
    }
  }  
  return classification;
}  

//_______________________________________________________________________________
std::vector<Int_t> AliCaloPhotonCuts::GetVectorMatchedTracksToCluster(AliVEvent* event, AliVCluster* cluster){
  vector<Int_t> labelsMatched(0);
  if(!fUseDistTrackToCluster) return labelsMatched;

  if (!fUsePtDepTrackToCluster)
    labelsMatched = fCaloTrackMatcher->GetMatchedTrackIDsForCluster(event, cluster->GetID(), fMaxDistTrackToClusterEta, -fMaxDistTrackToClusterEta,
                                                                    fMaxDistTrackToClusterPhi, fMinDistTrackToClusterPhi);
  else
    labelsMatched = fCaloTrackMatcher->GetMatchedTrackIDsForCluster(event, cluster->GetID(), fFuncPtDepEta, fFuncPtDepPhi);

  return labelsMatched;
}

//_______________________________________________________________________________
Bool_t AliCaloPhotonCuts::GetClosestMatchedTrackToCluster(AliVEvent* event, AliVCluster* cluster, Int_t &trackLabel){
  if(!fUseDistTrackToCluster) return kFALSE;
  vector<Int_t> labelsMatched = GetVectorMatchedTracksToCluster(event,cluster);

  if((Int_t) labelsMatched.size()<1) return kFALSE;

  Float_t dEta = -100;
  Float_t dPhi = -100;
  Int_t smallestDistTrack = -1;
  Float_t smallestDist = 100;
  Bool_t isSucessful = kFALSE;
  for (Int_t i = 0; i < (Int_t)labelsMatched.size(); i++){
    AliVTrack* currTrack  = dynamic_cast<AliVTrack*>(event->GetTrack(labelsMatched.at(i)));
    if(!fCaloTrackMatcher->GetTrackClusterMatchingResidual(currTrack->GetID(), cluster->GetID(), dEta, dPhi)) continue;

    Float_t tempDist = TMath::Sqrt((dEta*dEta)+(dPhi*dPhi));
    if(tempDist < smallestDist){
      smallestDist = tempDist;
      smallestDistTrack = labelsMatched.at(i);
      isSucessful = kTRUE;
    }
  }
  if(!isSucessful) return kFALSE;
  trackLabel = smallestDistTrack;
  return kTRUE;
}

//_______________________________________________________________________________
Bool_t AliCaloPhotonCuts::GetHighestPtMatchedTrackToCluster(AliVEvent* event, AliVCluster* cluster, Int_t &trackLabel){
  if(!fUseDistTrackToCluster) return kFALSE;
  vector<Int_t> labelsMatched = GetVectorMatchedTracksToCluster(event,cluster);

  if((Int_t) labelsMatched.size()<1) return kFALSE;

  Int_t highestPtTrack = -1;
  Float_t highestPt = 0;
  Bool_t isSucessful = kFALSE;
  for (Int_t i = 0; i < (Int_t)labelsMatched.size(); i++){
    AliVTrack* currTrack  = dynamic_cast<AliVTrack*>(event->GetTrack(labelsMatched.at(i)));
    if(currTrack->Pt() > highestPt){
      highestPt = currTrack->Pt();
      highestPtTrack = labelsMatched.at(i);
      isSucessful = kTRUE;
    }
  }
  if(!isSucessful) return kFALSE;
  trackLabel = highestPtTrack;
  return kTRUE;
}
