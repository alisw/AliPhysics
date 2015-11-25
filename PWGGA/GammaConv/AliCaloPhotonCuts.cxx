/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                         *
 * Authors: Friederike Bock, Baldo Sahlmueller                           *
 * Version 1.0                                                 *
 *                                                    *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is        *
 * provided "as is" without express or implied warranty.               *
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
  "ExoticCell",           //10   0: no cut
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
AliCaloPhotonCuts::AliCaloPhotonCuts(Bool_t isJetJet, const char *name,const char *title) :
  AliAnalysisCuts(name,title),
  fHistograms(NULL),
  fHistExtQA(NULL),
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
  fIsJetJet(kFALSE),
  fV0ReaderName("V0ReaderV1"),
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
  fExtendedMatchAndQA(0),
  fExoticCell(0),
  fUseExoticCell(0),
  fMinEnergy(0),
  fSeedEnergy(0.3),
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
  //fHistDistanceToBadChannelBeforeAcc(NULL),
  //fHistDistanceToBadChannelAfterAcc(NULL),
  fHistClusterTimevsEBeforeQA(NULL),
  fHistClusterTimevsEAfterQA(NULL),
  //fHistExoticCellBeforeQA(NULL),
  //fHistExoticCellAfterQA(NULL),
  //fHistNMatchedTracks(NULL),
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
  fHistClusterdPhidPtBeforeQA(NULL),
  fHistClusterM20M02BeforeQA(NULL),
  fHistClusterM20M02AfterQA(NULL)
{
  for(Int_t jj=0;jj<kNCuts;jj++){fCuts[jj]=0;}
  fCutString=new TObjString((GetCutNumber()).Data());

  if(isJetJet) fIsJetJet = kTRUE;
}

//________________________________________________________________________
AliCaloPhotonCuts::AliCaloPhotonCuts(const AliCaloPhotonCuts &ref) :
  AliAnalysisCuts(ref),
  fHistograms(NULL),
  fHistExtQA(NULL),
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
  fIsJetJet(ref.fIsJetJet),
  fV0ReaderName(ref.fV0ReaderName),
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
  fExtendedMatchAndQA(ref.fExtendedMatchAndQA),
  fExoticCell(ref.fExoticCell),
  fUseExoticCell(ref.fUseExoticCell),
  fMinEnergy(ref.fMinEnergy),
  fSeedEnergy(ref.fSeedEnergy),
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
  //fHistDistanceToBadChannelBeforeAcc(NULL),
  //fHistDistanceToBadChannelAfterAcc(NULL),
  fHistClusterTimevsEBeforeQA(NULL),
  fHistClusterTimevsEAfterQA(NULL),
  //fHistExoticCellBeforeQA(NULL),
  //fHistExoticCellAfterQA(NULL),
  //fHistNMatchedTracks(NULL),
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
  fHistClusterdPhidPtBeforeQA(NULL),
  fHistClusterM20M02BeforeQA(NULL),
  fHistClusterM20M02AfterQA(NULL)
{
  // Copy Constructor
  for(Int_t jj=0;jj<kNCuts;jj++){fCuts[jj]=ref.fCuts[jj];}
  fCutString=new TObjString((GetCutNumber()).Data());

}


//________________________________________________________________________
AliCaloPhotonCuts::~AliCaloPhotonCuts() {
  // Destructor
  //Deleting fHistograms leads to seg fault it it's added to output collection of a task
  // if(fHistograms)
  //    delete fHistograms;
  // fHistograms = NULL;
  if(fCutString != NULL){
    delete fCutString;
    fCutString = NULL;
  }

  if(fPHOSBadChannelsMap != NULL){
    delete[] fPHOSBadChannelsMap;
    fPHOSBadChannelsMap = NULL;
  }
}

//________________________________________________________________________
void AliCaloPhotonCuts::InitCutHistograms(TString name){

  // Initialize Cut Histograms for QA (only initialized and filled if function is called)
  TH1::AddDirectory(kFALSE);

  if(fHistograms != NULL){
    delete fHistograms;
    fHistograms=NULL;
  }
  if(fHistograms==NULL){
    fHistograms=new TList();
    fHistograms->SetOwner(kTRUE);
    if(name=="")fHistograms->SetName(Form("CaloCuts_%s",GetCutNumber().Data()));
    else fHistograms->SetName(Form("%s_%s",name.Data(),GetCutNumber().Data()));
  }

  if(fHistExtQA != NULL){
    delete fHistExtQA;
    fHistExtQA=NULL;
  }
  if(fHistExtQA==NULL){
    fHistExtQA=new TList();
    fHistExtQA->SetOwner(kTRUE);
    if(name=="")fHistExtQA->SetName(Form("CaloExtQA_%s",GetCutNumber().Data()));
    else fHistExtQA->SetName(Form("%s_ExtQA_%s",name.Data(),GetCutNumber().Data()));
  }

  // IsPhotonSelected
  fHistCutIndex=new TH1F(Form("IsPhotonSelected %s",GetCutNumber().Data()),"IsPhotonSelected",5,-0.5,4.5);
  fHistCutIndex->GetXaxis()->SetBinLabel(kPhotonIn+1,"in");
  fHistCutIndex->GetXaxis()->SetBinLabel(kDetector+1,"detector");
  fHistCutIndex->GetXaxis()->SetBinLabel(kAcceptance+1,"acceptance");
  fHistCutIndex->GetXaxis()->SetBinLabel(kClusterQuality+1,"cluster QA");
  fHistCutIndex->GetXaxis()->SetBinLabel(kPhotonOut+1,"out");
  fHistograms->Add(fHistCutIndex);

  // Acceptance Cuts
  fHistAcceptanceCuts=new TH1F(Form("AcceptanceCuts %s",GetCutNumber().Data()),"AcceptanceCuts",5,-0.5,4.5);
  fHistAcceptanceCuts->GetXaxis()->SetBinLabel(1,"in");
  fHistAcceptanceCuts->GetXaxis()->SetBinLabel(2,"eta");
  fHistAcceptanceCuts->GetXaxis()->SetBinLabel(3,"phi");
  fHistAcceptanceCuts->GetXaxis()->SetBinLabel(4,"distance to bad channel");
  fHistAcceptanceCuts->GetXaxis()->SetBinLabel(5,"out");
  fHistograms->Add(fHistAcceptanceCuts);

  // Cluster Cuts
  fHistClusterIdentificationCuts =new TH1F(Form("ClusterQualityCuts %s",GetCutNumber().Data()),"ClusterQualityCuts",11,-0.5,10.5);
  fHistClusterIdentificationCuts->GetXaxis()->SetBinLabel(1,"in");
  fHistClusterIdentificationCuts->GetXaxis()->SetBinLabel(2,"timing");
  fHistClusterIdentificationCuts->GetXaxis()->SetBinLabel(3,"track matching");
  fHistClusterIdentificationCuts->GetXaxis()->SetBinLabel(4,"Exotics");
  fHistClusterIdentificationCuts->GetXaxis()->SetBinLabel(5,"minimum energy");
  fHistClusterIdentificationCuts->GetXaxis()->SetBinLabel(6,"minimum NCells");
  fHistClusterIdentificationCuts->GetXaxis()->SetBinLabel(7,"NLM");
  fHistClusterIdentificationCuts->GetXaxis()->SetBinLabel(8,"M02");
  fHistClusterIdentificationCuts->GetXaxis()->SetBinLabel(9,"M20");
  fHistClusterIdentificationCuts->GetXaxis()->SetBinLabel(10,"dispersion");
  fHistClusterIdentificationCuts->GetXaxis()->SetBinLabel(11,"out");
  fHistograms->Add(fHistClusterIdentificationCuts);

  // Acceptance related histogramms
  if( fClusterType == 1 ){ //EMCAL
    const Int_t nEmcalEtaBins = 96;
    const Int_t nEmcalPhiBins = 124;
    Float_t EmcalEtaBins[nEmcalEtaBins+1] = {-0.66687,-0.653,-0.63913,-0.62526,-0.61139,-0.59752,-0.58365,-0.56978,-0.55591,-0.54204,-0.52817,-0.5143,-0.50043,-0.48656,-0.47269,-0.45882,-0.44495,-0.43108,-0.41721,-0.40334,-0.38947,-0.3756,-0.36173,-0.34786,-0.33399,-0.32012,-0.30625,-0.29238,-0.27851,-0.26464,-0.25077,-0.2369,-0.22303,-0.20916,-0.19529,-0.18142,-0.16755,-0.15368,-0.13981,-0.12594,-0.11207,-0.0982,-0.08433,-0.07046,-0.05659,-0.04272,-0.02885,-0.01498,-0.00111,0.01276,0.02663,0.0405,0.05437,0.06824,0.08211,0.09598,0.10985,0.12372,0.13759,0.15146,0.16533,0.1792,0.19307,0.20694,0.22081,0.23468,0.24855,0.26242,0.27629,0.29016,0.30403,0.3179,0.33177,0.34564,0.35951,0.37338,0.38725,0.40112,0.41499,0.42886,0.44273,0.4566,0.47047,0.48434,0.49821,0.51208,0.52595,0.53982,0.55369,0.56756,0.58143,0.5953,0.60917,0.62304,0.63691,0.65078,0.66465};
    Float_t EmcalPhiBins[nEmcalPhiBins+1] = {1.408,1.4215,1.435,1.4485,1.462,1.4755,1.489,1.5025,1.516,1.5295,1.543,1.5565,1.57,1.5835,1.597,1.6105,1.624,1.6375,1.651,1.6645,1.678,1.6915,1.705,1.7185,1.732,1.758,1.7715,1.785,1.7985,1.812,1.8255,1.839,1.8525,1.866,1.8795,1.893,1.9065,1.92,1.9335,1.947,1.9605,1.974,1.9875,2.001,2.0145,2.028,2.0415,2.055,2.0685,2.082,2.108,2.1215,2.135,2.1485,2.162,2.1755,2.189,2.2025,2.216,2.2295,2.243,2.2565,2.27,2.2835,2.297,2.3105,2.324,2.3375,2.351,2.3645,2.378,2.3915,2.405,2.4185,2.432,2.456,2.4695,2.483,2.4965,2.51,2.5235,2.537,2.5505,2.564,2.5775,2.591,2.6045,2.618,2.6315,2.645,2.6585,2.672,2.6855,2.699,2.7125,2.726,2.7395,2.753,2.7665,2.78,2.804,2.8175,2.831,2.8445,2.858,2.8715,2.885,2.8985,2.912,2.9255,2.939,2.9525,2.966,2.9795,2.993,3.0065,3.02,3.0335,3.047,3.0605,3.074,3.0875,3.101,3.1145,3.128};

    fHistClusterEtavsPhiBeforeAcc=new TH2F(Form("EtaPhi_beforeAcceptance %s",GetCutNumber().Data()),"EtaPhi_beforeAcceptance",nEmcalPhiBins,EmcalPhiBins,nEmcalEtaBins,EmcalEtaBins);
    fHistograms->Add(fHistClusterEtavsPhiBeforeAcc);
    fHistClusterEtavsPhiAfterAcc=new TH2F(Form("EtaPhi_afterAcceptance %s",GetCutNumber().Data()),"EtaPhi_afterAcceptance",nEmcalPhiBins,EmcalPhiBins,nEmcalEtaBins,EmcalEtaBins);
    fHistograms->Add(fHistClusterEtavsPhiAfterAcc);
    fHistClusterEtavsPhiAfterQA=new TH2F(Form("EtaPhi_afterClusterQA %s",GetCutNumber().Data()),"EtaPhi_afterClusterQA",nEmcalPhiBins,EmcalPhiBins,nEmcalEtaBins,EmcalEtaBins);
    fHistograms->Add(fHistClusterEtavsPhiAfterQA);
  } else if( fClusterType == 2 ){ //PHOS
    const Int_t nPhosEtaBins = 56;
    const Int_t nPhosPhiBins = 192;
    const Float_t PhosEtaRange[2] = {-0.16, 0.16};
    const Float_t PhosPhiRange[2] = {4.5, 5.6};

    fHistClusterEtavsPhiBeforeAcc=new TH2F(Form("EtaPhi_beforeAcceptance %s",GetCutNumber().Data()),"EtaPhi_beforeAcceptance",nPhosPhiBins,PhosPhiRange[0],PhosPhiRange[1],nPhosEtaBins,PhosEtaRange[0],PhosEtaRange[1]);
    fHistograms->Add(fHistClusterEtavsPhiBeforeAcc);
    fHistClusterEtavsPhiAfterAcc=new TH2F(Form("EtaPhi_afterAcceptance %s",GetCutNumber().Data()),"EtaPhi_afterAcceptance",nPhosPhiBins,PhosPhiRange[0],PhosPhiRange[1],nPhosEtaBins,PhosEtaRange[0],PhosEtaRange[1]);
    fHistograms->Add(fHistClusterEtavsPhiAfterAcc);
    fHistClusterEtavsPhiAfterQA=new TH2F(Form("EtaPhi_afterClusterQA %s",GetCutNumber().Data()),"EtaPhi_afterClusterQA",nPhosPhiBins,PhosPhiRange[0],PhosPhiRange[1],nPhosEtaBins,PhosEtaRange[0],PhosEtaRange[1]);
    fHistograms->Add(fHistClusterEtavsPhiAfterQA);
  }  else if( fClusterType == 0 ){ //all
        fHistClusterEtavsPhiBeforeAcc=new TH2F(Form("EtaPhi_beforeAcceptance %s",GetCutNumber().Data()),"EtaPhi_beforeAcceptance",462,0,2*TMath::Pi(),110,-0.7,0.7);
        fHistograms->Add(fHistClusterEtavsPhiBeforeAcc);
        fHistClusterEtavsPhiAfterAcc=new TH2F(Form("EtaPhi_afterAcceptance %s",GetCutNumber().Data()),"EtaPhi_afterAcceptance",462,0,2*TMath::Pi(),110,-0.7,0.7);
        fHistograms->Add(fHistClusterEtavsPhiAfterAcc);
        fHistClusterEtavsPhiAfterQA=new TH2F(Form("EtaPhi_afterClusterQA %s",GetCutNumber().Data()),"EtaPhi_afterClusterQA",462,0,2*TMath::Pi(),110,-0.7,0.7);
        fHistograms->Add(fHistClusterEtavsPhiAfterQA);
  } else{AliError(Form("Cluster Type is not EMCAL nor PHOS nor all: %i",fClusterType));}
//----------------
  if(fIsJetJet){
    fHistClusterEtavsPhiBeforeAcc->Sumw2();
    fHistClusterEtavsPhiAfterAcc->Sumw2();
    fHistClusterEtavsPhiAfterQA->Sumw2();
  }
//----------------
  //fHistDistanceToBadChannelBeforeAcc = new TH1F(Form("DistanceToBadChannel_beforeAcceptance %s",GetCutNumber().Data()),"DistanceToBadChannel_beforeAcceptance",200,0,40);
  //fHistograms->Add(fHistDistanceToBadChannelBeforeAcc);
  //fHistDistanceToBadChannelAfterAcc = new TH1F(Form("DistanceToBadChannel_afterAcceptance %s",GetCutNumber().Data()),"DistanceToBadChannel_afterAcceptance",200,0,40);
  //fHistograms->Add(fHistDistanceToBadChannelAfterAcc);
  
  // Cluster quality related histograms
  Double_t timeMin = -2e-6;
  Double_t timeMax = 8e-6;
  if( fClusterType == 1 ){
    timeMin = -2e-7;
    timeMax = 12e-7;
  }

  fHistClusterTimevsEBeforeQA=new TH2F(Form("ClusterTimeVsE_beforeClusterQA %s",GetCutNumber().Data()),"ClusterTimeVsE_beforeClusterQA",800,timeMin,timeMax,100,0,40);
  fHistograms->Add(fHistClusterTimevsEBeforeQA);
  fHistClusterTimevsEAfterQA=new TH2F(Form("ClusterTimeVsE_afterClusterQA %s",GetCutNumber().Data()),"ClusterTimeVsE_afterClusterQA",800,timeMin,timeMax,100,0,40);
  fHistograms->Add(fHistClusterTimevsEAfterQA);
    //fHistExoticCellBeforeQA=new TH2F(Form("ExoticCell_beforeClusterQA %s",GetCutNumber().Data()),"ExoticCell_beforeClusterQA",400,0,40,50,0.75,1);
    //fHistograms->Add(fHistExoticCellBeforeQA);
    //fHistExoticCellAfterQA=new TH2F(Form("ExoticCell_afterClusterQA %s",GetCutNumber().Data()),"ExoticCell_afterClusterQA",400,0,40,50,0.75,1);
    //fHistograms->Add(fHistExoticCellAfterQA);
    //fHistNMatchedTracks = new TH1F(Form("NMatchedTracks_%s",GetCutNumber().Data()),"NMatchedTracks",22,-1.5,20.5);
    //fHistograms->Add(fHistNMatchedTracks);
  if(fUseNonLinearity){
    fHistEnergyOfClusterBeforeNL = new TH1F(Form("EnergyOfCluster_beforeNonLinearity %s",GetCutNumber().Data()),"EnergyOfCluster_beforeNonLinearity",500,0,50);
    fHistograms->Add(fHistEnergyOfClusterBeforeNL);
    fHistEnergyOfClusterAfterNL = new TH1F(Form("EnergyOfCluster_afterNonLinearity %s",GetCutNumber().Data()),"EnergyOfCluster_afterNonLinearity",500,0,50);
    fHistograms->Add(fHistEnergyOfClusterAfterNL);
  }
  fHistEnergyOfClusterBeforeQA = new TH1F(Form("EnergyOfCluster_beforeClusterQA %s",GetCutNumber().Data()),"EnergyOfCluster_beforeClusterQA",500,0,50);
  fHistograms->Add(fHistEnergyOfClusterBeforeQA);
  fHistEnergyOfClusterAfterQA = new TH1F(Form("EnergyOfCluster_afterClusterQA %s",GetCutNumber().Data()),"EnergyOfCluster_afterClusterQA",500,0,50);
  fHistograms->Add(fHistEnergyOfClusterAfterQA);
  fHistNCellsBeforeQA = new TH1F(Form("NCellPerCluster_beforeClusterQA %s",GetCutNumber().Data()),"NCellPerCluster_beforeClusterQA",50,0,50);
  fHistograms->Add(fHistNCellsBeforeQA);
  fHistNCellsAfterQA = new TH1F(Form("NCellPerCluster_afterClusterQA %s",GetCutNumber().Data()),"NCellPerCluster_afterClusterQA",50,0,50);
  fHistograms->Add(fHistNCellsAfterQA);
  fHistM02BeforeQA = new TH1F(Form("M02_beforeClusterQA %s",GetCutNumber().Data()),"M02_beforeClusterQA",400,0,5);
  fHistograms->Add(fHistM02BeforeQA);
  fHistM02AfterQA = new TH1F(Form("M02_afterClusterQA %s",GetCutNumber().Data()),"M02_afterClusterQA",400,0,5);
  fHistograms->Add(fHistM02AfterQA);
  fHistM20BeforeQA = new TH1F(Form("M20_beforeClusterQA %s",GetCutNumber().Data()),"M20_beforeClusterQA",400,0,2.5);
  fHistograms->Add(fHistM20BeforeQA);
  fHistM20AfterQA = new TH1F(Form("M20_afterClusterQA %s",GetCutNumber().Data()),"M20_afterClusterQA",400,0,2.5);
  fHistograms->Add(fHistM20AfterQA);
  fHistDispersionBeforeQA = new TH1F(Form("Dispersion_beforeClusterQA %s",GetCutNumber().Data()),"Dispersion_beforeClusterQA",100,0,4);
  fHistograms->Add(fHistDispersionBeforeQA);
  fHistDispersionAfterQA = new TH1F(Form("Dispersion_afterClusterQA %s",GetCutNumber().Data()),"Dispersion_afterClusterQA",100,0,4);
  fHistograms->Add(fHistDispersionAfterQA);
  fHistNLMBeforeQA = new TH1F(Form("NLM_beforeClusterQA %s",GetCutNumber().Data()),"NLM_beforeClusterQA",10,0,10);
  fHistograms->Add(fHistNLMBeforeQA);
  fHistNLMAfterQA = new TH1F(Form("NLM_afterClusterQA %s",GetCutNumber().Data()),"NLM_afterClusterQA",10,0,10);
  fHistograms->Add(fHistNLMAfterQA);
  fHistNLMVsNCellsAfterQA = new TH2F(Form("NLM_NCells_afterClusterQA %s",GetCutNumber().Data()),"NLM_NCells_afterClusterQA",10,0,10,50,0,50);
  fHistograms->Add(fHistNLMVsNCellsAfterQA);
  fHistNLMVsEAfterQA = new TH2F(Form("NLM_E_afterClusterQA %s",GetCutNumber().Data()),"NLM_E_afterClusterQA",10,0,10,500,0,50);
  fHistograms->Add(fHistNLMVsEAfterQA);
  if(fExtendedMatchAndQA > 1 || fIsPureCalo > 0){
    fHistClusterEM02AfterQA = new TH2F(Form("EVsM02_afterClusterQA %s",GetCutNumber().Data()),"EVsM02_afterClusterQA",500,0,50,400,0,5);
    fHistograms->Add(fHistClusterEM02AfterQA);
    fHistClusterEM02BeforeQA = new TH2F(Form("EVsM02_beforeClusterQA %s",GetCutNumber().Data()),"EVsM02_beforeClusterQA",500,0,50,400,0,5);
    fHistograms->Add(fHistClusterEM02BeforeQA);
  }
//----------------
  if(fIsJetJet){
    fHistClusterTimevsEBeforeQA->Sumw2();
    fHistClusterTimevsEAfterQA->Sumw2();
    if(fUseNonLinearity){
      fHistEnergyOfClusterBeforeNL->Sumw2();
      fHistEnergyOfClusterAfterNL->Sumw2();
    }
    fHistEnergyOfClusterBeforeQA->Sumw2();
    fHistEnergyOfClusterAfterQA->Sumw2();
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
    fHistNLMVsNCellsAfterQA->Sumw2();
    fHistNLMVsEAfterQA->Sumw2();
    if(fExtendedMatchAndQA > 1 || fIsPureCalo > 0){
      fHistClusterEM02AfterQA->Sumw2();
      fHistClusterEM02BeforeQA->Sumw2();
    }
  }
//----------------
  if( fClusterType == 1 ){
    Int_t nMaxCellsEMCAL = fNMaxEMCalModules*48*24;
    fBadChannels = new TProfile("EMCal - Bad Channels","EMCal - Bad Channels",nMaxCellsEMCAL,0,nMaxCellsEMCAL);
    fHistExtQA->Add(fBadChannels);
  } else if( fClusterType == 2 ){
    Int_t nMaxCellsPHOS = fNMaxPHOSModules*56*64;
    fBadChannels = new TProfile("PHOS - Bad Channels","PHOS - Bad Channels",nMaxCellsPHOS,0,nMaxCellsPHOS);
    fHistExtQA->Add(fBadChannels);
  }
  
  if(fExtendedMatchAndQA > 1){
    if( fClusterType == 1 ){ //EMCAL
      Int_t nMaxCellsEMCAL = fNMaxEMCalModules*48*24;
      fHistClusterEnergyvsMod = new TH2F(Form("ClusterEnergyVsModule_afterClusterQA %s",GetCutNumber().Data()),"ClusterEnergyVsModule_afterClusterQA",500,0,50,fNMaxEMCalModules,0,fNMaxEMCalModules);
      fHistExtQA->Add(fHistClusterEnergyvsMod);
      fHistNCellsBigger100MeVvsMod = new TH2F(Form("NCellsAbove100VsModule %s",GetCutNumber().Data()),"NCellsAbove100VsModule",200,0,200,fNMaxEMCalModules,0,fNMaxEMCalModules);
      fHistExtQA->Add(fHistNCellsBigger100MeVvsMod);
      fHistNCellsBigger1500MeVvsMod = new TH2F(Form("NCellsAbove1500VsModule %s",GetCutNumber().Data()),"NCellsAbove1500VsModule",100,0,100,fNMaxEMCalModules,0,fNMaxEMCalModules);
      fHistExtQA->Add(fHistNCellsBigger1500MeVvsMod);
      fHistEnergyOfModvsMod = new TH2F(Form("ModuleEnergyVsModule %s",GetCutNumber().Data()),"ModuleEnergyVsModule",1000,0,100,fNMaxEMCalModules,0,fNMaxEMCalModules);
      fHistExtQA->Add(fHistEnergyOfModvsMod);
      fHistClusterEnergyvsNCells = new TH2F(Form("ClusterEnergyVsNCells_afterQA %s",GetCutNumber().Data()),"ClusterEnergyVsNCells_afterQA",300,0,30,50,0,50);
      fHistExtQA->Add(fHistClusterEnergyvsNCells);
      fHistCellEnergyvsCellID = new TH2F(Form("CellEnergyVsCellID %s",GetCutNumber().Data()),"CellEnergyVsCellID",300,0,30,nMaxCellsEMCAL,0,nMaxCellsEMCAL);
      fHistExtQA->Add(fHistCellEnergyvsCellID);
      fHistCellTimevsCellID = new TH2F(Form("CellTimeVsCellID %s",GetCutNumber().Data()),"CellTimeVsCellID",1200,-timeMax,timeMax,nMaxCellsEMCAL,0,nMaxCellsEMCAL);
      fHistExtQA->Add(fHistCellTimevsCellID);
      fHistClusterIncludedCellsBeforeQA = new TH1F(Form("ClusterIncludedCells_beforeClusterQA %s",GetCutNumber().Data()),"ClusterIncludedCells_beforeClusterQA",nMaxCellsEMCAL,0,nMaxCellsEMCAL);
      fHistExtQA->Add(fHistClusterIncludedCellsBeforeQA);
      fHistClusterIncludedCellsAfterQA = new TH1F(Form("ClusterIncludedCells_afterClusterQA %s",GetCutNumber().Data()),"ClusterIncludedCells_afterClusterQA",nMaxCellsEMCAL,0,nMaxCellsEMCAL);
      fHistExtQA->Add(fHistClusterIncludedCellsAfterQA);
      fHistClusterEnergyFracCellsBeforeQA = new TH1F(Form("ClusterEnergyFracCells_beforeClusterQA %s",GetCutNumber().Data()),"ClusterEnergyFracCells_beforeClusterQA",nMaxCellsEMCAL,0,nMaxCellsEMCAL);
      fHistClusterEnergyFracCellsBeforeQA->Sumw2();
      fHistExtQA->Add(fHistClusterEnergyFracCellsBeforeQA);
      fHistClusterEnergyFracCellsAfterQA = new TH1F(Form("ClusterEnergyFracCells_afterClusterQA %s",GetCutNumber().Data()),"ClusterEnergyFracCells_afterClusterQA",nMaxCellsEMCAL,0,nMaxCellsEMCAL);
      fHistClusterEnergyFracCellsAfterQA->Sumw2();
      fHistExtQA->Add(fHistClusterEnergyFracCellsAfterQA);
    }
    else if( fClusterType == 2 ){ //PHOS
      Int_t nMaxCellsPHOS = fNMaxPHOSModules*56*64;
      fHistClusterEnergyvsMod = new TH2F(Form("ClusterEnergyVsModule_afterClusterQA %s",GetCutNumber().Data()),"ClusterEnergyVsModule_afterClusterQA",500,0,50,fNMaxPHOSModules,0,fNMaxPHOSModules);
      fHistExtQA->Add(fHistClusterEnergyvsMod);
      fHistNCellsBigger100MeVvsMod = new TH2F(Form("NCellsAbove100VsModule %s",GetCutNumber().Data()),"NCellsAbove100VsModule",200,0,200,fNMaxPHOSModules,0,fNMaxPHOSModules);
      fHistExtQA->Add(fHistNCellsBigger100MeVvsMod);
      fHistNCellsBigger1500MeVvsMod = new TH2F(Form("NCellsAbove1500VsModule %s",GetCutNumber().Data()),"NCellsAbove1500VsModule",100,0,100,fNMaxPHOSModules,0,fNMaxPHOSModules);
      fHistExtQA->Add(fHistNCellsBigger1500MeVvsMod);
      fHistEnergyOfModvsMod = new TH2F(Form("ModuleEnergyVsModule %s",GetCutNumber().Data()),"ModuleEnergyVsModule",1000,0,100,fNMaxPHOSModules,0,fNMaxPHOSModules);
      fHistExtQA->Add(fHistEnergyOfModvsMod);
      fHistClusterEnergyvsNCells = new TH2F(Form("ClusterEnergyVsNCells_afterQA %s",GetCutNumber().Data()),"ClusterEnergyVsNCells_afterQA",300,0,30,50,0,50);
      fHistExtQA->Add(fHistClusterEnergyvsNCells);
      fHistCellEnergyvsCellID = new TH2F(Form("CellEnergyVsCellID %s",GetCutNumber().Data()),"CellEnergyVsCellID",300,0,30,nMaxCellsPHOS,0,nMaxCellsPHOS);
      fHistExtQA->Add(fHistCellEnergyvsCellID);
      fHistCellTimevsCellID = new TH2F(Form("CellTimeVsCellID %s",GetCutNumber().Data()),"CellTimeVsCellID",1200,-timeMax,timeMax,nMaxCellsPHOS,0,nMaxCellsPHOS);
      fHistExtQA->Add(fHistCellTimevsCellID);
      fHistClusterIncludedCellsBeforeQA = new TH1F(Form("ClusterIncludedCells_beforeClusterQA %s",GetCutNumber().Data()),"ClusterIncludedCells_beforeClusterQA",nMaxCellsPHOS,0,nMaxCellsPHOS);
      fHistExtQA->Add(fHistClusterIncludedCellsBeforeQA);
      fHistClusterIncludedCellsAfterQA = new TH1F(Form("ClusterIncludedCells_afterClusterQA %s",GetCutNumber().Data()),"ClusterIncludedCells_afterClusterQA",nMaxCellsPHOS,0,nMaxCellsPHOS);
      fHistExtQA->Add(fHistClusterIncludedCellsAfterQA);
      fHistClusterEnergyFracCellsBeforeQA = new TH1F(Form("ClusterEnergyFracCells_beforeClusterQA %s",GetCutNumber().Data()),"ClusterEnergyFracCells_beforeClusterQA",nMaxCellsPHOS,0,nMaxCellsPHOS);
      fHistClusterEnergyFracCellsBeforeQA->Sumw2();
      fHistExtQA->Add(fHistClusterEnergyFracCellsBeforeQA);
      fHistClusterEnergyFracCellsAfterQA = new TH1F(Form("ClusterEnergyFracCells_afterClusterQA %s",GetCutNumber().Data()),"ClusterEnergyFracCells_afterClusterQA",nMaxCellsPHOS,0,nMaxCellsPHOS);
      fHistClusterEnergyFracCellsAfterQA->Sumw2();
      fHistExtQA->Add(fHistClusterEnergyFracCellsAfterQA);
    }
    else{AliError(Form("fExtendedMatchAndQA (%i) not (yet) defined for cluster type (%i)",fExtendedMatchAndQA,fClusterType));}
  }

  //TrackMatching histograms
  if(fUseDistTrackToCluster) {
    const Int_t nEtaBins = 300;
    const Int_t nPhiBins = 300;
    const Float_t EtaRange[2] = {-0.3, 0.3};
    const Float_t PhiRange[2] = {-0.3, 0.3};

    fHistClusterRBeforeQA = new TH1F(Form("R_Cluster_beforeClusterQA %s",GetCutNumber().Data()),"R of cluster",200,400,500);
    fHistograms->Add(fHistClusterRBeforeQA);
    fHistClusterRAfterQA = new TH1F(Form("R_Cluster_afterClusterQA %s",GetCutNumber().Data()),"R of cluster_matched",200,400,500);
    fHistograms->Add(fHistClusterRAfterQA);
    fHistClusterdEtadPhiBeforeQA=new TH2F(Form("dEtaVsdPhi_beforeClusterQA %s",GetCutNumber().Data()),"dEtaVsdPhi_beforeClusterQA",nEtaBins,EtaRange[0],EtaRange[1],nPhiBins,PhiRange[0],PhiRange[1]);
    fHistograms->Add(fHistClusterdEtadPhiBeforeQA);
    fHistClusterdEtadPhiAfterQA=new TH2F(Form("dEtaVsdPhi_afterClusterQA %s",GetCutNumber().Data()),"dEtaVsdPhi_afterClusterQA",nEtaBins,EtaRange[0],EtaRange[1],nPhiBins,PhiRange[0],PhiRange[1]);
    fHistograms->Add(fHistClusterdEtadPhiAfterQA);
    fHistDistanceTrackToClusterBeforeQA = new TH1F(Form("DistanceToTrack_beforeClusterQA %s",GetCutNumber().Data()),"DistanceToTrack_beforeClusterQA",200,0,2);
    fHistograms->Add(fHistDistanceTrackToClusterBeforeQA);
    fHistDistanceTrackToClusterAfterQA = new TH1F(Form("DistanceToTrack_afterClusterQA %s",GetCutNumber().Data()),"DistanceToTrack_afterClusterQA",200,0,2);
    fHistograms->Add(fHistDistanceTrackToClusterAfterQA);
    //----------------
    if(fIsJetJet){
      fHistClusterRBeforeQA->Sumw2();
      fHistClusterRAfterQA->Sumw2();
      fHistClusterdEtadPhiBeforeQA->Sumw2();
      fHistClusterdEtadPhiAfterQA->Sumw2();
      fHistDistanceTrackToClusterBeforeQA->Sumw2();
      fHistDistanceTrackToClusterAfterQA->Sumw2();
    }
    //----------------
    if(fExtendedMatchAndQA > 0 && fExtendedMatchAndQA < 3){
      fHistClusterdEtadPhiPosTracksBeforeQA = new TH2F(Form("dEtaVsdPhi_posTracks_beforeClusterQA %s",GetCutNumber().Data()),"dEtaVsdPhi_posTracks_beforeClusterQA",nEtaBins,EtaRange[0],EtaRange[1],nPhiBins,PhiRange[0],PhiRange[1]);
      fHistograms->Add(fHistClusterdEtadPhiPosTracksBeforeQA);
      fHistClusterdEtadPhiNegTracksBeforeQA = new TH2F(Form("dEtaVsdPhi_negTracks_beforeClusterQA %s",GetCutNumber().Data()),"dEtaVsdPhi_negTracks_beforeClusterQA",nEtaBins,EtaRange[0],EtaRange[1],nPhiBins,PhiRange[0],PhiRange[1]);
      fHistograms->Add(fHistClusterdEtadPhiNegTracksBeforeQA);
      fHistClusterdEtadPhiPosTracksAfterQA = new TH2F(Form("dEtaVsdPhi_posTracks_afterClusterQA %s",GetCutNumber().Data()),"dEtaVsdPhi_posTracks_afterClusterQA",nEtaBins,EtaRange[0],EtaRange[1],nPhiBins,PhiRange[0],PhiRange[1]);
      fHistograms->Add(fHistClusterdEtadPhiPosTracksAfterQA);
      fHistClusterdEtadPhiNegTracksAfterQA = new TH2F(Form("dEtaVsdPhi_negTracks_afterClusterQA %s",GetCutNumber().Data()),"dEtaVsdPhi_negTracks_afterClusterQA",nEtaBins,EtaRange[0],EtaRange[1],nPhiBins,PhiRange[0],PhiRange[1]);
      fHistograms->Add(fHistClusterdEtadPhiNegTracksAfterQA);
      fHistClusterdEtadPhiPosTracksP_000_075BeforeQA = new TH2F(Form("dEtaVsdPhi_posTracks_P<0.75_beforeClusterQA %s",GetCutNumber().Data()),"dEtaVsdPhi_posTracks_P<0.75_beforeClusterQA",nEtaBins,EtaRange[0],EtaRange[1],nPhiBins,PhiRange[0],PhiRange[1]);
      fHistograms->Add(fHistClusterdEtadPhiPosTracksP_000_075BeforeQA);
      fHistClusterdEtadPhiPosTracksP_075_125BeforeQA = new TH2F(Form("dEtaVsdPhi_posTracks_0.75<P<1.25_beforeClusterQA %s",GetCutNumber().Data()),"dEtaVsdPhi_posTracks_0.75<P<1.25_beforeClusterQA",nEtaBins,EtaRange[0],EtaRange[1],nPhiBins,PhiRange[0],PhiRange[1]);
      fHistograms->Add(fHistClusterdEtadPhiPosTracksP_075_125BeforeQA);
      fHistClusterdEtadPhiPosTracksP_125_999BeforeQA = new TH2F(Form("dEtaVsdPhi_posTracks_P>1.25_beforeClusterQA %s",GetCutNumber().Data()),"dEtaVsdPhi_posTracks_P>1.25_beforeClusterQA",nEtaBins,EtaRange[0],EtaRange[1],nPhiBins,PhiRange[0],PhiRange[1]);
      fHistograms->Add(fHistClusterdEtadPhiPosTracksP_125_999BeforeQA);
      fHistClusterdEtadPhiNegTracksP_000_075BeforeQA= new TH2F(Form("dEtaVsdPhi_negTracks_P<0.75_beforeClusterQA %s",GetCutNumber().Data()),"dEtaVsdPhi_negTracks_P<0.75_beforeClusterQA",nEtaBins,EtaRange[0],EtaRange[1],nPhiBins,PhiRange[0],PhiRange[1]);
      fHistograms->Add(fHistClusterdEtadPhiNegTracksP_000_075BeforeQA);
      fHistClusterdEtadPhiNegTracksP_075_125BeforeQA = new TH2F(Form("dEtaVsdPhi_negTracks_0.75<P<1.25_beforeClusterQA %s",GetCutNumber().Data()),"dEtaVsdPhi_negTracks_0.75<P<1.25_beforeClusterQA",nEtaBins,EtaRange[0],EtaRange[1],nPhiBins,PhiRange[0],PhiRange[1]);
      fHistograms->Add(fHistClusterdEtadPhiNegTracksP_075_125BeforeQA);
      fHistClusterdEtadPhiNegTracksP_125_999BeforeQA = new TH2F(Form("dEtaVsdPhi_negTracks_P>1.25_beforeClusterQA %s",GetCutNumber().Data()),"dEtaVsdPhi_negTracks_P>1.25_beforeClusterQA",nEtaBins,EtaRange[0],EtaRange[1],nPhiBins,PhiRange[0],PhiRange[1]);
      fHistograms->Add(fHistClusterdEtadPhiNegTracksP_125_999BeforeQA);
      fHistClusterdEtadPtBeforeQA = new TH2F(Form("dEtaVsPt_beforeClusterQA %s",GetCutNumber().Data()),"dEtaVsPt_beforeClusterQA",nEtaBins,EtaRange[0],EtaRange[1],300,0,30);
      fHistograms->Add(fHistClusterdEtadPtBeforeQA);
      fHistClusterdPhidPtBeforeQA = new TH2F(Form("dPhiVsPt_beforeClusterQA %s",GetCutNumber().Data()),"dPhiVsPt_beforeClusterQA",2*nPhiBins,2*PhiRange[0],2*PhiRange[1],300,0,30);
      fHistograms->Add(fHistClusterdPhidPtBeforeQA);
      fHistClusterM20M02BeforeQA = new TH2F(Form("M20VsM02_beforeClusterQA %s",GetCutNumber().Data()),"M20VsM02_beforeClusterQA",200,0,2.5,400,0,5);
      fHistograms->Add(fHistClusterM20M02BeforeQA);
      fHistClusterM20M02AfterQA = new TH2F(Form("M20VsM02_afterClusterQA %s",GetCutNumber().Data()),"M20VsM02_afterClusterQA",200,0,2.5,400,0,5);
      fHistograms->Add(fHistClusterM20M02AfterQA);

      //----------------
      if(fIsJetJet){
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
        fHistClusterdPhidPtBeforeQA->Sumw2();
        fHistClusterM20M02BeforeQA->Sumw2();
        fHistClusterM20M02AfterQA->Sumw2();
      }
//----------------
    }
  }

  fVectorMatchedClusterIDs.clear();

  TH1::AddDirectory(kTRUE);
  return;
}

//________________________________________________________________________
void AliCaloPhotonCuts::InitializeEMCAL(AliVEvent *event){

  if (fClusterType == 1){
    AliTender* alitender=0x0;
    AliEmcalTenderTask* emcaltender=0x0;
  
    if(event->IsA()==AliESDEvent::Class()) alitender = (AliTender*) AliAnalysisManager::GetAnalysisManager()->GetTopTasks()->FindObject("AliTender");
    else if( event->IsA()==AliAODEvent::Class()) emcaltender = (AliEmcalTenderTask*) AliAnalysisManager::GetAnalysisManager()->GetTopTasks()->FindObject("AliEmcalTenderTask");
    
    if(alitender){
      TIter next(alitender->GetSupplies());
      AliTenderSupply *supply;
      while ((supply=(AliTenderSupply*)next())) if(supply->IsA()==AliEMCALTenderSupply::Class()) break;
      fEMCALRecUtils = ((AliEMCALTenderSupply*)supply)->GetRecoUtils();
      fEMCALBadChannelsMap = fEMCALRecUtils->GetEMCALBadChannelStatusMapArray();
    } else if(emcaltender){
      fEMCALRecUtils = ((AliEMCALTenderSupply*)emcaltender->GetEMCALTenderSupply())->GetRecoUtils();
      fEMCALBadChannelsMap = fEMCALRecUtils->GetEMCALBadChannelStatusMapArray();
    }
    if (fEMCALRecUtils) fEMCALInitialized = kTRUE;


    Int_t nMaxCellsEMCAL = fNMaxEMCalModules*48*24;
    Int_t imod = -1;Int_t iTower = -1, iIphi = -1, iIeta = -1;
    Int_t icol = -1;Int_t irow = -1;

    fGeomEMCAL = AliEMCALGeometry::GetInstance();
    if(!fGeomEMCAL){ AliFatal("EMCal geometry not initialized!");}

    for(Int_t iCell=0;iCell<nMaxCellsEMCAL;iCell++){
      fGeomEMCAL->GetCellIndex(iCell,imod,iTower,iIphi,iIeta);
      if (fEMCALBadChannelsMap->GetEntries() <= imod) continue;
      fGeomEMCAL->GetCellPhiEtaIndexInSModule(imod,iTower,iIphi,iIeta,irow,icol);
      Int_t iBadCell = (Int_t) ((TH2I*)fEMCALBadChannelsMap->At(imod))->GetBinContent(icol,irow);
      if(iBadCell > 0) fBadChannels->Fill(iCell,1);
      else fBadChannels->Fill(iCell,0);
    }
  }
  return;
}

//________________________________________________________________________
void AliCaloPhotonCuts::InitializePHOS (AliVEvent *event){

  if (fClusterType == 2){
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

    fPHOSInitialized = kTRUE;
    fPHOSCurrentRun = event->GetRunNumber();
  }
  return;
}

//________________________________________________________________________
Bool_t AliCaloPhotonCuts::ClusterIsSelectedMC(TParticle *particle,AliStack *fMCStack){
   // MonteCarlo Photon Selection

  if(!fMCStack)return kFALSE;

  if (particle->GetPdgCode() == 22){

    if ( particle->Eta() < fMinEtaCut || particle->Eta() > fMaxEtaCut ) return kFALSE;
    if ( particle->Phi() < fMinPhiCut || particle->Phi() > fMaxPhiCut ) return kFALSE;
    
    if(particle->GetMother(0) >-1 && fMCStack->Particle(particle->GetMother(0))->GetPdgCode() == 22){
      return kFALSE;// no photon as mothers!
    }
    // decision on prim/sec moved to main task
//     if(particle->GetMother(0) >= fMCStack->GetNprimary()){
//       return kFALSE;// the gamma has a mother, and it is not a primary particle
//     }
    return kTRUE;
  }
  return kFALSE;
}

//________________________________________________________________________
Bool_t AliCaloPhotonCuts::ClusterIsSelectedAODMC(AliAODMCParticle *particle,TClonesArray *aodmcArray){
  // MonteCarlo Photon Selection

  if(!aodmcArray)return kFALSE;

  if (particle->GetPdgCode() == 22){
    if ( particle->Eta() < fMinEtaCut || particle->Eta() > fMaxEtaCut ) return kFALSE;
    if ( particle->Phi() < fMinPhiCut || particle->Phi() > fMaxPhiCut ) return kFALSE;
    if(particle->GetMother() > -1 && (static_cast<AliAODMCParticle*>(aodmcArray->At(particle->GetMother())))->GetPdgCode() == 22){
        return kFALSE;// no photon as mothers!
    }
    // decision on prim/sec moved to main task
//    Bool_t isPrimary = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsConversionPrimaryAOD(fInputEvent, static_cast<AliAODMCParticle*>(aodmcArray->At(particle->GetMother())), mcProdVtxX, mcProdVtxY, mcProdVtxZ);
//    if(!isPrimary){
//      return kFALSE;// the gamma has a mother, and it is not a primary particle
//    }
    return kTRUE;// return in case of accepted gamma
  }
  return kFALSE;
}

//________________________________________________________________________
// This function selects the clusters based on their quality criteria
//________________________________________________________________________
Bool_t AliCaloPhotonCuts::ClusterQualityCuts(AliVCluster* cluster, AliVEvent *event, Int_t isMC, Double_t weight)
{   // Specific Photon Cuts

  Int_t cutIndex = 0;
  if(fHistClusterIdentificationCuts)fHistClusterIdentificationCuts->Fill(cutIndex);
  cutIndex++;

//   Double_t minR = 999.0;
  // get the minimum radius of tracks to cluster
//   if(fHistDistanceTrackToClusterBeforeQA || fHistDistanceTrackToClusterAfterQA){
//     Float_t pos[3];
//     cluster->GetPosition(pos);// Get cluster position
//     
//     TVector3 cp(pos);
//     int NtrMatched = 0;
//     NtrMatched = cluster->GetNTracksMatched();
//     fHistNMatchedTracks->Fill(NtrMatched);
//     //loop over tracks for Jet QA
//     TList *l = event->GetList();
//     TClonesArray *tracks = dynamic_cast<TClonesArray*>(l->FindObject("Tracks"));
//     for(int itrack = 0;itrack < NtrMatched;itrack++){
//       AliVTrack *trackcluster = static_cast<AliVTrack*>(tracks->At(itrack));
//       if (! trackcluster) {
//         AliError(Form("Couldn't get ESD track %d\n", itrack));
//         continue;
//       }
//       Double_t dphi = -999.0;
//       Double_t deta = -999.0;
//       AliPicoTrack::GetEtaPhiDiff(trackcluster, cluster, dphi, deta);
//       cout << "here" << endl;
//       Double_t dr = sqrt(dphi*dphi + deta+deta);
//       if(dr < minR)
//         minR = dr;
//     }//loop over tracks
//   }

  Int_t nLM = GetNumberOfLocalMaxima(cluster, event);

  // Fill Histos before Cuts
  if(fHistClusterTimevsEBeforeQA) fHistClusterTimevsEBeforeQA->Fill(cluster->GetTOF(), cluster->E(), weight);
//   if(fHistExoticCellBeforeQA) fHistExoticCellBeforeQA->Fill(cluster->E(), weight);
//   if(fHistDistanceTrackToClusterBeforeQA) fHistDistanceTrackToClusterBeforeQA->Fill(minR, weight);
  if(fHistEnergyOfClusterBeforeQA) fHistEnergyOfClusterBeforeQA->Fill(cluster->E(), weight);
  if(fHistNCellsBeforeQA) fHistNCellsBeforeQA->Fill(cluster->GetNCells(), weight);
  if(fHistM02BeforeQA) fHistM02BeforeQA->Fill(cluster->GetM02(), weight);
  if(fHistM20BeforeQA) fHistM20BeforeQA->Fill(cluster->GetM20(), weight);
  if(fHistDispersionBeforeQA) fHistDispersionBeforeQA->Fill(cluster->GetDispersion(), weight);
  if(fHistNLMBeforeQA) fHistNLMBeforeQA->Fill(nLM, weight);
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
      if(fHistClusterIdentificationCuts)fHistClusterIdentificationCuts->Fill(cutIndex);//1
      return kFALSE;
    }
  }  
  cutIndex++;//2, next cut

  // Minimum distance to track
//   if (fUseDistTrackToCluster){
//      Float_t pos[3];
//      cluster->GetPosition(pos);// Get cluster position
//      TVector3 cp(pos);
//      int NtrMatched = 0;
//      NtrMatched = cluster->GetNTracksMatched();
//     fHistNMatchedTracks->Fill(NtrMatched);
//      if(NtrMatched > 0){
//       //loop over tracks for QA
//       TList *l = event->GetList();
//       TClonesArray *tracks = dynamic_cast<TClonesArray*>(l->FindObject("Tracks"));
//       
//       Double_t dphi = 999.0;
//       Double_t deta = 999.0;
//       Double_t dr2 = 999.0;
// 
//       for(int itrack = 0;itrack < NtrMatched;itrack++){
//         AliVTrack *trackcluster = NULL;
//         trackcluster = static_cast<AliVTrack*>(tracks->At(itrack));
//         if (! trackcluster) {
//         AliError(Form("Couldn't get ESD track %d\n", itrack));
//         continue;
//         }
//         AliPicoTrack::GetEtaPhiDiff(trackcluster, cluster, dphi, deta);
//         dr2 = dphi*dphi + deta+deta;
//         //cout << dr << endl;
//         if(dr2 < fMinDistTrackToCluster*fMinDistTrackToCluster){
//     //        if(dphi < fMinDistTrackToCluster || deta < fMinDistTrackToCluster){
//         if(fHistClusterIdentificationCuts)fHistClusterIdentificationCuts->Fill(cutIndex);//2
//         return kFALSE;
//         }
//         
//       }//loop over tracks
//     }
// //    if(cluster->GetEmcCpvDistance() < fMinDistTrackToCluster){
// //      if(fHistClusterIdentificationCuts)fHistClusterIdentificationCuts->Fill(cutIndex);//2
// //      return kFALSE;
// //    }
//   }
  if (fIsPureCalo>0 && fVectorMatchedClusterIDs.size()>0){
    if( CheckClusterForTrackMatch(cluster) ){
      if(fHistClusterIdentificationCuts)fHistClusterIdentificationCuts->Fill(cutIndex);//2
      return kFALSE;
    }
  }
  cutIndex++;//3, next cut

  // exotic cell cut --IMPLEMENT LATER---
//   if(!AcceptanceCuts(photon)){
//     if(fHistClusterIdentificationCuts)fHistClusterIdentificationCuts->Fill(cutIndex);//3
//     return kFALSE;
//   }
  cutIndex++;//4, next cut
  
  // minimum cell energy cut
  if (fUseMinEnergy){
    if(cluster->E() < fMinEnergy){
      if(fHistClusterIdentificationCuts)fHistClusterIdentificationCuts->Fill(cutIndex);//4
      return kFALSE;
    }
  }  
  cutIndex++;//5, next cut

  // minimum number of cells
  if (fUseNCells){
    if(cluster->GetNCells() < fMinNCells) {
      if(fHistClusterIdentificationCuts)fHistClusterIdentificationCuts->Fill(cutIndex);//5
      return kFALSE;
    }
  }  
  cutIndex++;//6, next cut

  // NLM cut
  if (fUseNLM){
    if( nLM < fMinNLM || nLM > fMaxNLM ) {
      if(fHistClusterIdentificationCuts)fHistClusterIdentificationCuts->Fill(cutIndex);//9
      return kFALSE;
    }
  }  
  cutIndex++;//7, next cut
  

  // M02 cut
  if (fUseM02 == 1){
    if( cluster->GetM02()< fMinM02 || cluster->GetM02() > fMaxM02 ) {
      if(fHistClusterIdentificationCuts)fHistClusterIdentificationCuts->Fill(cutIndex);//6
      return kFALSE;
    }
  } else if (fUseM02 ==2 ) {
    if( cluster->GetM02()< CalculateMinM02(fMinM02CutNr, cluster->E()) || 
      cluster->GetM02() > CalculateMaxM02(fMaxM02CutNr, cluster->E()) ) {
      if(fHistClusterIdentificationCuts)fHistClusterIdentificationCuts->Fill(cutIndex);//6
      return kFALSE;
    }
  }
  cutIndex++;//8, next cut

  // M20 cut
  if (fUseM20){
    if( cluster->GetM20()< fMinM20 || cluster->GetM20() > fMaxM20 ) {
      if(fHistClusterIdentificationCuts)fHistClusterIdentificationCuts->Fill(cutIndex);//7
      return kFALSE;
    }
  }  
  cutIndex++;//9, next cut
  
  // dispersion cut
  if (fUseDispersion){
    if( cluster->GetDispersion()> fMaxDispersion) {
      if(fHistClusterIdentificationCuts)fHistClusterIdentificationCuts->Fill(cutIndex);//8
      return kFALSE;
    }
  }  
  cutIndex++;//10, next cut
  
  
  // DONE with selecting photons
  if(fHistClusterIdentificationCuts)fHistClusterIdentificationCuts->Fill(cutIndex);//10

  // Histos after Cuts
    Float_t clusPos[3]={0,0,0};
    cluster->GetPosition(clusPos);
    TVector3 clusterVector(clusPos[0],clusPos[1],clusPos[2]);
  Double_t etaCluster = clusterVector.Eta();
  Double_t phiCluster = clusterVector.Phi();
  if (phiCluster < 0) phiCluster += 2*TMath::Pi();

  if(fHistClusterEtavsPhiAfterQA) fHistClusterEtavsPhiAfterQA->Fill(phiCluster, etaCluster, weight);
  if(fHistClusterTimevsEAfterQA) fHistClusterTimevsEAfterQA->Fill(cluster->GetTOF(), cluster->E(), weight);
//   if(fHistExoticCellAfterQA) fHistExoticCellAfterQA->Fill(cluster->E(), weight);
//   if(fHistDistanceTrackToClusterAfterQA) fHistDistanceTrackToClusterAfterQA->Fill(minR, weight);
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
        fHistClusterIncludedCellsAfterQA->Fill(cluster->GetCellAbsId(iCell));
        if(cluster->E()>0.) fHistClusterEnergyFracCellsAfterQA->Fill(cluster->GetCellAbsId(iCell),cells->GetCellAmplitude(cluster->GetCellAbsId(iCell))/cluster->E());
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
void AliCaloPhotonCuts::FillHistogramsExtendedQA(AliVEvent *event)
{
  if(fExtendedMatchAndQA < 2) return;

  AliVCaloCells* cells;

  Int_t nModules = 0;
  Int_t* nCellsBigger100MeV;
  Int_t* nCellsBigger1500MeV;
  Double_t* EnergyOfMod;

  if( fClusterType == 1 && !fEMCALInitialized ) InitializeEMCAL(event);
  if( fClusterType == 2 && ( !fPHOSInitialized || (fPHOSCurrentRun != event->GetRunNumber()) ) ) InitializePHOS(event);
  
  if( fClusterType == 1 ){ //EMCAL
    cells = event->GetEMCALCells();
    fGeomEMCAL = AliEMCALGeometry::GetInstance();
    if(!fGeomEMCAL) AliFatal("EMCal geometry not initialized!");
    if(!fEMCALBadChannelsMap) AliFatal("EMCal bad channels map not initialized!");
    nModules = fGeomEMCAL->GetNumberOfSuperModules();
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

    Int_t imod = -1;Int_t iTower = -1, iIphi = -1, iIeta = -1;
    Int_t icol = -1;Int_t irow = -1;
    Int_t relid[4];// for PHOS

    Bool_t doBadCell = kTRUE;
    if( fClusterType == 1 ){
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
    if( fClusterType == 1 && doBadCell){
      iBadCell = (Int_t) ((TH2I*)fEMCALBadChannelsMap->At(imod))->GetBinContent(icol,irow);
    }else if( fClusterType == 2 && doBadCell){
      iBadCell = (Int_t) ((TH2I*)fPHOSBadChannelsMap[nMod])->GetBinContent(relid[2],relid[3]);
    }

    if(iBadCell > 0) continue;

    if(cellAmplitude > 0.1) nCellsBigger100MeV[nMod]++;
    if(cellAmplitude > 1.5) nCellsBigger1500MeV[nMod]++;
    if(cellAmplitude > 0.05) EnergyOfMod[nMod]+=cellAmplitude;
      
    if(fHistCellEnergyvsCellID && (cellAmplitude > 0.05)) fHistCellEnergyvsCellID->Fill(cellAmplitude,cellNumber);
    if(fHistCellTimevsCellID && (cellAmplitude > 0.2)) fHistCellTimevsCellID->Fill(cellTime,cellNumber);
  }

  for(Int_t iModule=0;iModule<nModules;iModule++){
    if(fHistNCellsBigger100MeVvsMod) fHistNCellsBigger100MeVvsMod->Fill(nCellsBigger100MeV[iModule],iModule);
    if(fHistNCellsBigger1500MeVvsMod) fHistNCellsBigger1500MeVvsMod->Fill(nCellsBigger1500MeV[iModule],iModule);
    if(fHistEnergyOfModvsMod) fHistEnergyOfModvsMod->Fill(EnergyOfMod[iModule],iModule);
  }

  delete[] nCellsBigger100MeV;nCellsBigger100MeV=0x0;
  delete[] nCellsBigger1500MeV;nCellsBigger1500MeV=0x0;
  delete[] EnergyOfMod;EnergyOfMod=0x0;

  return;
}

//________________________________________________________________________
//************** Find number of local maxima in cluster ******************
//* derived from G. C. Balbastre's AliCalorimeterUtils *******************
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
  
  if (fClusterType == 1 ) 
    cells                 = event->GetEMCALCells();
  else if (fClusterType ==2 ) 
    cells                 = event->GetPHOSCells();
  
//   cout << "NCells: "<< nCells<< " cluster energy: " << cluster->E() << endl;
  Float_t eMax            = 0.;
  Int_t idMax             = -1;
  Int_t idMax2            = -1;
  
  if (nCells < 2) return idMax;
  for (Int_t iCell = 1;iCell < nCells;iCell++){
    if (cells->GetCellAmplitude(cluster->GetCellAbsId(iCell))> eMax){
      eMax                = cells->GetCellAmplitude(cluster->GetCellAbsId(iCell));
      idMax               = cluster->GetCellAbsId(iCell);
    }  
  }  
  
  eMax                    = 0.;
  for (Int_t iCell = 1;iCell < nCells;iCell++){
    if (iCell == idMax) continue;
    if (cells->GetCellAmplitude(cluster->GetCellAbsId(iCell))> eMax){
      eMax                = cells->GetCellAmplitude(cluster->GetCellAbsId(iCell));
      idMax2              = cluster->GetCellAbsId(iCell);
    }  
  }  

  return idMax2;
}

//________________________________________________________________________
Int_t AliCaloPhotonCuts::FindLargestCellInCluster(AliVCluster* cluster, AliVEvent* event){

  const Int_t nCells      = cluster->GetNCells();
  AliVCaloCells* cells    = NULL;

  if (fClusterType == 1 ) 
    cells                 = event->GetEMCALCells();
  else if (fClusterType ==2 ) 
    cells                 = event->GetPHOSCells();
  
//   cout << "NCells: "<< nCells<< " cluster energy: " << cluster->E() << endl;
  Float_t eMax            = 0.;
  Int_t idMax             = -1;
  
  if (nCells < 1) return idMax;
  for (Int_t iCell = 0;iCell < nCells;iCell++){
    Int_t cellAbsID       = cluster->GetCellAbsId(iCell);
    if (cells->GetCellAmplitude(cellAbsID)> eMax){
      eMax                = cells->GetCellAmplitude(cellAbsID);
      idMax               = cellAbsID;
    }
  }
  return idMax;
  
}


//________________________________________________________________________
//************** Find number of local maxima in cluster ******************
//* derived from G. C. Balbastre's AliCalorimeterUtils *******************
//************************************************************************
Int_t AliCaloPhotonCuts::GetNumberOfLocalMaxima(AliVCluster* cluster, AliVEvent * event, Int_t *absCellIdList, Float_t* maxEList){

  Int_t absCellId1        = -1;
  Int_t absCellId2        = -1;
  const Int_t nCells      = cluster->GetNCells();
  AliVCaloCells* cells    = NULL;
  
  if (fClusterType == 1 ) 
    cells                 = event->GetEMCALCells();
  else if (fClusterType ==2 ) 
    cells                 = event->GetPHOSCells();
  
//   cout << "NCells: "<< nCells<< " cluster energy: " << cluster->E() << endl;
  Float_t eMax            = 0.;
  Int_t idMax             = -1;
  
  for (Int_t iCell = 0;iCell < nCells;iCell++){
    absCellIdList[iCell]  = cluster->GetCellAbsId(iCell);
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
      absCellId1          = absCellIdList[iCell];
      Float_t en1         = cells->GetCellAmplitude(absCellId1);
      if (en1 < fSeedEnergy) {
        absCellIdList[iCell] = -1;
        continue;
      }  
      // loop over other cells in cluster
      for (Int_t iCellN = 0;iCellN < nCells;iCellN++){
        // jump out if array has changed in the meantime
        if (absCellIdList[iCell] == -1) continue;
        // don't compare to yourself
        if (iCell == iCellN) continue;
        
        // get cell id & check whether its valid
        absCellId2        = absCellIdList[iCellN];
        if (absCellId2 == -1) continue;
        
        // get cell energy of second cell
        Float_t en2       = cells->GetCellAmplitude(absCellId2);
        
        // check if cells are Neighbours
        if (AreNeighbours(absCellId1, absCellId2)){
          // determine which cell has larger energy, mask the other
//           cout << "found neighbour: " << absCellId1 << "\t" << absCellId2 << endl;
//           cout << "energies: " << en1 << "\t" << en2 << endl;
          if (en1 > en2){
            absCellIdList[iCellN] = -1;
          } else {
            absCellIdList[iCell]  = -1;
            continue;
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
      absCellIdList[nMaximaNew]   = absCellIdList[iCell];
      maxEList[nMaximaNew]        = cells->GetCellAmplitude(absCellIdList[iCell]);
      nMaximaNew++;
    }  
  }  

  // check whether a local maximum was found
  // if no maximum was found use highest cell as maximum
  if (nMaximaNew == 0){
    nMaximaNew            =  1;
    maxEList[0]           = eMax;
    absCellIdList[0]      = idMax;
  }  

  return nMaximaNew;
}  

//________________________________________________________________________
//************** Function to determine neighbours of cells ***************
//* derived from G. C. Balbastre's AliCalorimeterUtils *******************
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
//* derived from G. C. Balbastre's AliCalorimeterUtils *******************
//************************************************************************
Int_t AliCaloPhotonCuts::GetModuleNumberAndCellPosition(Int_t absCellId, Int_t & icol, Int_t & irow){
  if( fClusterType == 1 ){ //EMCAL
    fGeomEMCAL = AliEMCALGeometry::GetInstance();
    if(!fGeomEMCAL) AliFatal("EMCal geometry not initialized!");
  } else if( fClusterType == 2 ){ //PHOS
    fGeomPHOS = AliPHOSGeometry::GetInstance();
    if(!fGeomPHOS) AliFatal("PHOS geometry not initialized!");
  }
  
  Int_t imod = -1;Int_t iTower = -1, iIphi = -1, iIeta = -1;
  if( fClusterType == 1 ){
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
//* derived from G. C. Balbastre's AliCalorimeterUtils *******************
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
  if (fClusterType == 1 ) 
    cells                 = event->GetEMCALCells();
  else if (fClusterType ==2 ) 
    cells                 = event->GetPHOSCells();

  Float_t e1              = 0;
  Float_t e2              = 0;
  Float_t eCluster        = 0;
  
  for(Int_t iCell    = 0;iCell < ncells;iCell++ ) {
    absCellIdList[iCell]  = cluster->GetCellAbsId(iCell);
    Float_t ec            = cells->GetCellAmplitude(absCellIdList[iCell]);
    eCluster+=ec;
  }

  UShort_t absCellIdList1  [12];
  Double_t fracList1     [12];
  UShort_t absCellIdList2  [12];
  Double_t fracList2     [12];

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
  if (fClusterType == 1){
    CorrectEMCalNonLinearity(cluster1, isMC) ;
    CorrectEMCalNonLinearity(cluster2, isMC) ;
  }

  // Initialize EMCAL rec utils if not initialized
  if(!fEMCALInitialized && fClusterType == 1) InitializeEMCAL(event);
  
  if(fEMCALInitialized && fClusterType == 1){
    fEMCALRecUtils->RecalculateClusterPosition(fGeomEMCAL, cells, cluster1);
    fEMCALRecUtils->RecalculateClusterPosition(fGeomEMCAL, cells, cluster2);
  }
}

//________________________________________________________________________
Bool_t AliCaloPhotonCuts::CheckDistanceToBadChannel(AliVCluster* cluster, AliVEvent* event)
{
  if(fUseDistanceToBadChannel != 1 && fUseDistanceToBadChannel != 2) return kFALSE;

  //not yet fully implemented for PHOS:
  if( fClusterType == 2) return kFALSE;

  if( fClusterType == 1 && !fEMCALInitialized ) InitializeEMCAL(event);
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
Bool_t AliCaloPhotonCuts::ClusterIsSelected(AliVCluster *cluster, AliVEvent * event, Int_t isMC, Double_t weight)
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
    if (fClusterType == 1 && !cluster->IsEMCAL()){
      FillClusterCutIndex(kDetector);
      return kFALSE;
    }
    //Select PHOS cluster
    if (fClusterType == 2 && !cluster->IsPHOS()){
      FillClusterCutIndex(kDetector);
      return kFALSE;
    }
  }
  
  // Acceptance Cuts
  if(!AcceptanceCuts(cluster,event,weight)){
    FillClusterCutIndex(kAcceptance);
    return kFALSE;
  }
  // Cluster Quality Cuts
  if(!ClusterQualityCuts(cluster,event,isMC,weight)){
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

//  Double_t vertex[3] = {0,0,0};
//   event->GetPrimaryVertex()->GetXYZ(vertex);
  // TLorentzvector with cluster
//  TLorentzVector clusterVector;
//  cluster->GetMomentum(clusterVector,vertex);

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
  cutIndex++;
  if(fHistAcceptanceCuts)fHistAcceptanceCuts->Fill(cutIndex);

  // Histos after cuts
  if(fHistClusterEtavsPhiAfterAcc) fHistClusterEtavsPhiAfterAcc->Fill(phiCluster,etaCluster,weight);
  
  return kTRUE;
}

//________________________________________________________________________
Bool_t AliCaloPhotonCuts::MatchConvPhotonToCluster(AliAODConversionPhoton* convPhoton, AliVCluster* cluster, AliVEvent* event, Double_t weight){

  if (!fUseDistTrackToCluster) return kFALSE;
  
  AliESDEvent *esdev = dynamic_cast<AliESDEvent*>(event);
  AliAODEvent *aodev = 0;
  if (!esdev) {
    aodev = dynamic_cast<AliAODEvent*>(event);
    if (!aodev) {
      AliError("Task needs AOD or ESD event, returning");
      return kFALSE;
    }
  }

//    Double_t vertex[3] = {0,0,0};
//    event->GetPrimaryVertex()->GetXYZ(vertex);

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
      if(((AliV0ReaderV1*)AliAnalysisManager::GetAnalysisManager()->GetTask("V0ReaderV1"))->AreAODsRelabeled()){
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

    AliESDtrack *esdt = dynamic_cast<AliESDtrack*>(inTrack);
    AliAODTrack *aodt = 0;
    if (!esdt) {
      aodt = dynamic_cast<AliAODTrack*>(inTrack);
      if (!aodt){AliError("Track is neither ESD nor AOD, continue");continue;}
    }

    AliExternalTrackParam *trackParam = 0;
    if (esdt) {
      const AliExternalTrackParam *in = esdt->GetInnerParam();
      if (!in){AliError("Could not get InnerParam of Track, continue");continue;}
      trackParam = new AliExternalTrackParam(*in);
    } else {
      Double_t xyz[3] = {0}, pxpypz[3] = {0}, cv[21] = {0};
      aodt->PxPyPz(pxpypz);
      aodt->XvYvZv(xyz);
      aodt->GetCovarianceXYZPxPyPz(cv);
      trackParam = new AliExternalTrackParam(xyz,pxpypz,cv,aodt->Charge());
    }
    if (!trackParam){AliError("Could not get TrackParameters, continue");continue;}
    
    Bool_t propagated = kFALSE;
    AliExternalTrackParam emcParam(*trackParam);
    Float_t dPhi = 0;
    Float_t dEta = 0;

    if(cluster->IsEMCAL()){
      Float_t eta = 0;Float_t phi = 0;Float_t pt = 0;
      propagated = AliEMCALRecoUtils::ExtrapolateTrackToEMCalSurface(&emcParam, 430, 0.000510999, 20, eta, phi, pt);
      if(propagated){
        propagated = AliEMCALRecoUtils::ExtrapolateTrackToCluster(&emcParam, cluster, 0.000510999, 5, dEta, dPhi);
      }
    }
    if(cluster->IsPHOS()){
      propagated = AliTrackerBase::PropagateTrackToBxByBz(&emcParam, clusterR, 0.000510999, 20, kTRUE, 0.8, -1);
      if (propagated){
        Double_t trkPos[3] = {0,0,0};
        emcParam.GetXYZ(trkPos);
        TVector3 trkPosVec(trkPos[0],trkPos[1],trkPos[2]);
        TVector3 clsPosVec(clusterPosition);
        dPhi = clsPosVec.DeltaPhi(trkPosVec);
        dEta = clsPosVec.Eta()-trkPosVec.Eta();
      }
    }

    if (propagated){
      Float_t dR2 = dPhi*dPhi + dEta*dEta;
      if(fHistDistanceTrackToClusterBeforeQA)fHistDistanceTrackToClusterBeforeQA->Fill(TMath::Sqrt(dR2), weight);
      if(fHistClusterdEtadPhiBeforeQA) fHistClusterdEtadPhiBeforeQA->Fill(dEta, dPhi, weight);

            Float_t clusM02 = (Float_t) cluster->GetM02();
            Float_t clusM20 = (Float_t) cluster->GetM20();
      if(fExtendedMatchAndQA > 0 && fExtendedMatchAndQA < 3){
        if(inTrack->Charge() > 0) {
          fHistClusterdEtadPhiPosTracksBeforeQA->Fill(dEta, dPhi, weight);
          if(inTrack->P() < 0.75) fHistClusterdEtadPhiPosTracksP_000_075BeforeQA->Fill(dEta, dPhi, weight);
          else if(inTrack->P() < 1.25) fHistClusterdEtadPhiPosTracksP_075_125BeforeQA->Fill(dEta, dPhi, weight);
          else fHistClusterdEtadPhiPosTracksP_125_999BeforeQA->Fill(dEta, dPhi, weight);
        } else {
          fHistClusterdEtadPhiNegTracksBeforeQA->Fill(dEta, dPhi, weight);
          if(inTrack->P() < 0.75) fHistClusterdEtadPhiNegTracksP_000_075BeforeQA->Fill(dEta, dPhi, weight);
          else if(inTrack->P() < 1.25) fHistClusterdEtadPhiNegTracksP_075_125BeforeQA->Fill(dEta, dPhi, weight);
          else fHistClusterdEtadPhiNegTracksP_125_999BeforeQA->Fill(dEta, dPhi, weight);
        }
        fHistClusterdEtadPtBeforeQA->Fill(dEta, inTrack->Pt(), weight);
        fHistClusterdPhidPtBeforeQA->Fill(dPhi, inTrack->Pt(), weight);
        fHistClusterM20M02BeforeQA->Fill(clusM20, clusM02, weight);
      }

      Bool_t match_dEta = (abs(dEta) < fMaxDistTrackToClusterEta) ? kTRUE : kFALSE;
      Bool_t match_dPhi = kFALSE;
      if( (inTrack->Charge() > 0) && (dPhi > fMinDistTrackToClusterPhi) && (dPhi < fMaxDistTrackToClusterPhi) ) match_dPhi = kTRUE;
      else if( (inTrack->Charge() < 0) && (dPhi < -fMinDistTrackToClusterPhi) && (dPhi > -fMaxDistTrackToClusterPhi) ) match_dPhi = kTRUE;

      if(match_dEta && match_dPhi){
            //if(dR2 < fMinDistTrackToCluster*fMinDistTrackToCluster){
        matched = kTRUE;
      } else {
        if(fHistDistanceTrackToClusterAfterQA)fHistDistanceTrackToClusterAfterQA->Fill(TMath::Sqrt(dR2), weight);
        if(fHistClusterdEtadPhiAfterQA) fHistClusterdEtadPhiAfterQA->Fill(dEta, dPhi, weight);
        if(fHistClusterRAfterQA) fHistClusterRAfterQA->Fill(clusterR, weight);
        if(fExtendedMatchAndQA > 0 && fExtendedMatchAndQA < 3){
          if(inTrack->Charge() > 0) fHistClusterdEtadPhiPosTracksAfterQA->Fill(dEta, dPhi, weight);
          else fHistClusterdEtadPhiNegTracksAfterQA->Fill(dEta, dPhi, weight);
          fHistClusterM20M02AfterQA->Fill(clusM20, clusM02, weight);
        }
      }  
    }
    delete trackParam;
  }

  return matched;

}

//________________________________________________________________________
void AliCaloPhotonCuts::MatchTracksToClusters(AliVEvent* event, Double_t weight){

  if( fIsPureCalo == 0 || !fUseDistTrackToCluster ) return;

  // not yet fully implemented + tested for PHOS
  if( fClusterType != 1) return;

  fVectorMatchedClusterIDs.clear();

  Int_t nClus = event->GetNumberOfCaloClusters();
  Int_t nModules = 0;

  if(fClusterType == 1){
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

//  cout << "MatchTracksToClusters: " << event->GetNumberOfTracks() << ", " << fIsPureCalo << ", " << fUseDistTrackToCluster << endl;

  for (Int_t itr=0;itr<event->GetNumberOfTracks();itr++){
    AliExternalTrackParam *trackParam = 0;
    AliVTrack *inTrack = 0x0;
//    cout << "-------------------------LOOPING: " << itr << endl;
    if(esdev){
      inTrack = esdev->GetTrack(itr);
      AliESDtrack *esdt = dynamic_cast<AliESDtrack*>(inTrack);

      const AliExternalTrackParam *in = esdt->GetInnerParam();
      if (!in){AliError("Could not get InnerParam of Track, continue");continue;}
      trackParam = new AliExternalTrackParam(*in);
    } else if(aodev) {
      if(((AliV0ReaderV1*)AliAnalysisManager::GetAnalysisManager()->GetTask("V0ReaderV1"))->AreAODsRelabeled()){
        inTrack = dynamic_cast<AliVTrack*>(aodev->GetTrack(itr));
      } else {
        for(Int_t ii=0;ii<aodev->GetNumberOfTracks();ii++) {
          inTrack = dynamic_cast<AliVTrack*>(aodev->GetTrack(ii));
          if(inTrack){
            if(inTrack->GetID() == itr) {
              break;
            }
          }
        }
      }

      AliAODTrack *aodt = dynamic_cast<AliAODTrack*>(inTrack);
      Double_t xyz[3] = {0}, pxpypz[3] = {0}, cv[21] = {0};
      aodt->PxPyPz(pxpypz);
      aodt->XvYvZv(xyz);
      aodt->GetCovarianceXYZPxPyPz(cv);
      trackParam = new AliExternalTrackParam(xyz,pxpypz,cv,aodt->Charge());
    }

    if (!trackParam) {AliError("Could not get TrackParameters, continue");continue;}
    AliExternalTrackParam emcParam(*trackParam);
    Float_t eta, phi, pt;

    //propagate tracks to emc surfaces
    if(fClusterType == 1){
      if (!AliEMCALRecoUtils::ExtrapolateTrackToEMCalSurface(&emcParam, 440., 0.139, 20., eta, phi, pt)) {
        delete trackParam;
        continue;
      }

      if( TMath::Abs(eta) > 0.75 ) {
        delete trackParam;
        continue;
      }
      // Save some time and memory in case of no DCal present
      if( nModules < 13 && ( phi < 70*TMath::DegToRad() || phi > 190*TMath::DegToRad())){
        delete trackParam;
        continue;
      }
    }else if(fClusterType == 2){
      if( !AliTrackerBase::PropagateTrackToBxByBz(&emcParam, 460., 0.139, 20, kTRUE, 0.8, -1)){
        delete trackParam;
        continue;
      }
//to do: implement of distance checks
    }

    Float_t dEta=-999, dPhi=-999;
    Float_t clsPos[3] = {0.,0.,0.};
    Double_t exPos[3] = {0.,0.,0.};
    if (!emcParam.GetXYZ(exPos)) continue;

//    cout << "-: " << trackParam << endl;
//    cout << "eta/phi: " << eta << ", " << phi << endl;
//    cout << "nClus: " << nClus << endl;
    for(Int_t iclus=0;iclus < nClus;iclus++){
      AliVCluster* cluster = event->GetCaloCluster(iclus);
      if (!cluster) continue;
//      cout << "-------------------------LOOPING: " << iclus << ", " << cluster->GetID() << endl;
      cluster->GetPosition(clsPos);
      Double_t dR = TMath::Sqrt(TMath::Power(exPos[0]-clsPos[0],2)+TMath::Power(exPos[1]-clsPos[1],2)+TMath::Power(exPos[2]-clsPos[2],2));
      if (dR > 100) continue;
      Double_t clusterR = TMath::Sqrt( clsPos[0]*clsPos[0] + clsPos[1]*clsPos[1] );
//      TVector3 vecC(clsPos);
//      cout << "eta/phi: " << vecC.Eta() << ", " << vecC.Phi() << endl;
      AliExternalTrackParam trackParamTmp(emcParam);//Retrieve the starting point every time before the extrapolation
      if(fClusterType == 1){
        if (!cluster->IsEMCAL()) continue;
        if(!AliEMCALRecoUtils::ExtrapolateTrackToCluster(&trackParamTmp, cluster, 0.139, 5., dEta, dPhi)) continue;
      }else if(fClusterType == 2){
        if (!cluster->IsPHOS()) continue;
        if(!AliTrackerBase::PropagateTrackToBxByBz(&trackParamTmp, clusterR, 0.139, 5., kTRUE, 0.8, -1)) continue;
        Double_t trkPos[3] = {0,0,0};
        trackParamTmp.GetXYZ(trkPos);
        TVector3 trkPosVec(trkPos[0],trkPos[1],trkPos[2]);
        TVector3 clsPosVec(clsPos);
        dPhi = clsPosVec.DeltaPhi(trkPosVec);
        dEta = clsPosVec.Eta()-trkPosVec.Eta();
      }

      Float_t dR2 = dPhi*dPhi + dEta*dEta;
//      cout << "dEta/dPhi: " << dEta << ", " << dPhi << " - ";
//      cout << dR2 << endl;
      if(fHistDistanceTrackToClusterBeforeQA)fHistDistanceTrackToClusterBeforeQA->Fill(TMath::Sqrt(dR2), weight);
      if(fHistClusterdEtadPhiBeforeQA) fHistClusterdEtadPhiBeforeQA->Fill(dEta, dPhi, weight);

      Float_t clusM02 = (Float_t) cluster->GetM02();
      Float_t clusM20 = (Float_t) cluster->GetM20();
      if(fExtendedMatchAndQA > 0 && fExtendedMatchAndQA < 3){
        if(inTrack->Charge() > 0) {
          fHistClusterdEtadPhiPosTracksBeforeQA->Fill(dEta, dPhi, weight);
          if(inTrack->P() < 0.75) fHistClusterdEtadPhiPosTracksP_000_075BeforeQA->Fill(dEta, dPhi, weight);
          else if(inTrack->P() < 1.25) fHistClusterdEtadPhiPosTracksP_075_125BeforeQA->Fill(dEta, dPhi, weight);
          else fHistClusterdEtadPhiPosTracksP_125_999BeforeQA->Fill(dEta, dPhi, weight);
        }
        else{
          fHistClusterdEtadPhiNegTracksBeforeQA->Fill(dEta, dPhi, weight);
          if(inTrack->P() < 0.75) fHistClusterdEtadPhiNegTracksP_000_075BeforeQA->Fill(dEta, dPhi, weight);
          else if(inTrack->P() < 1.25) fHistClusterdEtadPhiNegTracksP_075_125BeforeQA->Fill(dEta, dPhi, weight);
          else fHistClusterdEtadPhiNegTracksP_125_999BeforeQA->Fill(dEta, dPhi, weight);
        }
        fHistClusterdEtadPtBeforeQA->Fill(dEta, inTrack->Pt(), weight);
        fHistClusterdPhidPtBeforeQA->Fill(dPhi, inTrack->Pt(), weight);
        fHistClusterM20M02BeforeQA->Fill(clusM20, clusM02, weight);
      }

      Bool_t match_dEta = (abs(dEta) < fMaxDistTrackToClusterEta) ? kTRUE : kFALSE;
      Bool_t match_dPhi = kFALSE;
      if( (inTrack->Charge() > 0) && (dPhi > fMinDistTrackToClusterPhi) && (dPhi < fMaxDistTrackToClusterPhi) ) match_dPhi = kTRUE;
      else if( (inTrack->Charge() < 0) && (dPhi < -fMinDistTrackToClusterPhi) && (dPhi > -fMaxDistTrackToClusterPhi) ) match_dPhi = kTRUE;

      if(match_dEta && match_dPhi){
        fVectorMatchedClusterIDs.push_back(cluster->GetID());
//        cout << "MATCHED!!!!!!!!!!!!!!!!!!!!!!!!! - " << cluster->GetID() << endl;
        break;
      } else {
        if(fHistDistanceTrackToClusterAfterQA)fHistDistanceTrackToClusterAfterQA->Fill(TMath::Sqrt(dR2), weight);
        if(fHistClusterdEtadPhiAfterQA) fHistClusterdEtadPhiAfterQA->Fill(dEta, dPhi, weight);
        if(fHistClusterRAfterQA) fHistClusterRAfterQA->Fill(clusterR, weight);
        if(fExtendedMatchAndQA > 0 && fExtendedMatchAndQA < 3){
          if(inTrack->Charge() > 0) fHistClusterdEtadPhiPosTracksAfterQA->Fill(dEta, dPhi, weight);
          else fHistClusterdEtadPhiNegTracksAfterQA->Fill(dEta, dPhi, weight);
          fHistClusterM20M02AfterQA->Fill(clusM20, clusM02, weight);
        }
//        cout << "no match" << endl;
      }
    }

    delete trackParam;
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

    case kExoticCell:
      if( SetExoticCellCut(value)) {
        fCuts[kExoticCell] = value;
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
  if (fUseEtaCut) printf("\t%3.2f < eta_{cluster} < %3.2f\n", fMinEtaCut, fMaxEtaCut );
  if (fUsePhiCut) printf("\t%3.2f < phi_{cluster} < %3.2f\n", fMinPhiCut, fMaxPhiCut );
  if (fUseDistanceToBadChannel>0) printf("\tdistance to bad channel used in mode '%i', distance in cells: %f \n",fUseDistanceToBadChannel, fMinDistanceToBadChannel);
  
  printf("Cluster Quality cuts: \n");
  if (fUseTimeDiff) printf("\t %3.2f < time difference < %3.2f\n", fMinTimeDiff, fMaxTimeDiff );
  if (fUseDistTrackToCluster) printf("\tmin distance to track in eta > %3.2f, min phi < %3.2f and max phi > %3.2f\n", fMaxDistTrackToClusterEta, fMinDistTrackToClusterPhi, fMaxDistTrackToClusterPhi );
  if (fUseExoticCell)printf("\t exotic cell: %3.2f\n", fExoticCell );
  if (fUseMinEnergy)printf("\t E_{cluster} > %3.2f\n", fMinEnergy );
  if (fUseNCells) printf("\t number of cells per cluster >= %d\n", fMinNCells );
  if (fUseM02 == 1) printf("\t %3.2f < M02 < %3.2f\n", fMinM02, fMaxM02 );
  if (fUseM02 == 2) printf("\t energy dependent M02 cut used with cutnumber min: %d  max: %d \n", fMinM02CutNr, fMaxM02CutNr );
  if (fUseM20) printf("\t %3.2f < M20 < %3.2f\n", fMinM20, fMaxM20 );
  if (fUseDispersion) printf("\t dispersion < %3.2f\n", fMaxDispersion );
  if (fUseNLM) printf("\t %d < NLM < %d\n", fMinNLM, fMaxNLM );

  printf("NonLinearity Correction: \n");
  if (fUseNonLinearity) printf("\t Chose NonLinearity cut '%i'\n", fSwitchNonLinearity);
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
  // standard: cluster - V0-track matching
  if(fIsPureCalo == 0){
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
    case 6:
      if (!fUseDistTrackToCluster) fUseDistTrackToCluster=kTRUE;
      fMaxDistTrackToClusterEta = 0.022;//0.015;
      fMinDistTrackToClusterPhi = -0.15;//-0.02;
      fMaxDistTrackToClusterPhi = 0.10;//0.055;//0.2;
      break;
    case 7: //PHOS
      if (!fUseDistTrackToCluster) fUseDistTrackToCluster=kTRUE;
      fMaxDistTrackToClusterEta = 0.005;//0.015;
      fMinDistTrackToClusterPhi = -0.03;//-0.025;
      fMaxDistTrackToClusterPhi = 0.03;//0.06;//0.3;
      break;
    case 8: //PHOS
      if (!fUseDistTrackToCluster) fUseDistTrackToCluster=kTRUE;
      fMaxDistTrackToClusterEta = 0.01;//0.015;
      fMinDistTrackToClusterPhi = -0.09;//-0.025;
      fMaxDistTrackToClusterPhi = 0.07;//0.07;//0.4;
      break;
    case 9: //PHOS
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
  // if merged cluster analysis running: use cut for general track - cluster matching
  }else{
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
    case 6:
      if (!fUseDistTrackToCluster) fUseDistTrackToCluster=kTRUE;
      fMaxDistTrackToClusterEta = 0.022;//0.015;
      fMinDistTrackToClusterPhi = -0.15;//-0.02;
      fMaxDistTrackToClusterPhi = 0.10;//0.055;//0.2;
      break;
    case 7: //PHOS
      if (!fUseDistTrackToCluster) fUseDistTrackToCluster=kTRUE;
      fMaxDistTrackToClusterEta = 0.005;//0.015;
      fMinDistTrackToClusterPhi = -0.03;//-0.025;
      fMaxDistTrackToClusterPhi = 0.03;//0.06;//0.3;
      break;
    case 8: //PHOS
      if (!fUseDistTrackToCluster) fUseDistTrackToCluster=kTRUE;
      fMaxDistTrackToClusterEta = 0.01;//0.015;
      fMinDistTrackToClusterPhi = -0.09;//-0.025;
      fMaxDistTrackToClusterPhi = 0.07;//0.07;//0.4;
      break;
    case 9: //PHOS
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
Bool_t AliCaloPhotonCuts::SetExoticCellCut(Int_t exoticCell)
{
  switch(exoticCell){
  case 0: 
    fUseExoticCell=0;
    fExoticCell=0;
    break;
  case 1: 
    if (!fUseExoticCell) fUseExoticCell=1;
    fExoticCell=5;
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
    switch(minEnergy){
    case 0: 
      if (!fUseMinEnergy) fUseMinEnergy=0;
      fMinEnergy=0.1;
      break;
    case 1: 
      if (!fUseMinEnergy) fUseMinEnergy=1;
      fMinEnergy=0.2;
      break;
    case 2: 
      if (!fUseMinEnergy) fUseMinEnergy=1;
      fMinEnergy=0.3;
      break;
    case 3: 
      if (!fUseMinEnergy) fUseMinEnergy=1;
      fMinEnergy=0.4;
      break;
    case 4: 
      if (!fUseMinEnergy) fUseMinEnergy=1;
      fMinEnergy=0.5;
      break;
    case 5: 
      if (!fUseMinEnergy) fUseMinEnergy=1;
      fMinEnergy=0.6;
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
    return kTRUE;
  } else   {
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

  if(isMC>0 && fCurrentMC==kNoMC){
    AliV0ReaderV1* V0Reader = (AliV0ReaderV1*) AliAnalysisManager::GetAnalysisManager()->GetTask(fV0ReaderName.Data());
    if( V0Reader == NULL ){
      AliFatal(Form("No V0Reader called '%s' could be found within AliCaloPhotonCuts::CorrectEMCalNonLinearity",fV0ReaderName.Data()));
      return;
    }
    fPeriodName = V0Reader->GetPeriodName();
    fCurrentMC = FindEnumForMCSet(fPeriodName);
  }
  Bool_t fPeriodNameAvailable = kTRUE;

  switch(fSwitchNonLinearity){

    // Standard NonLinearity - standard kPi0MCv5 for MC and kSDMv5 for data from Jason
    case 1:
      energy *= FunctionNL_kPi0MCv5(energy);
      if(isMC == 0) energy *= FunctionNL_kSDMv5(energy);
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

    // NonLinearity LHC12 ConvCalo - only shifting MC
    case 11:
      label_case_11:
      if(isMC>0){
        //pass1
        if( fCurrentMC==k14e2a || fCurrentMC==k14e2b )
          energy /= FunctionNL_kSDM(energy, 0.983251, -3.44339, -1.70998);

        else if( fCurrentMC==k14e2c )
          energy /= FunctionNL_kSDM(energy, 0.984462, -3.00363, -2.63773);

        //pass2
        else if( fCurrentMC == k15h1 )
          energy /= FunctionNL_kSDM(energy, 0.983251, -3.44339, -1.70998);

        else if( fCurrentMC == k15h2 )
          energy /= FunctionNL_kSDM(energy, 0.984462, -3.00363, -2.63773);

        else fPeriodNameAvailable = kFALSE;
      }
      break;

    // NonLinearity LHC12 Calo - only shifting MC
    case 12:
      label_case_12:
      if(isMC>0){
        //pass1
        if( fCurrentMC==k14e2a || fCurrentMC==k14e2b )
          energy /= FunctionNL_kSDM(2.0*energy, 0.967301, -3.1683, -0.653058);

        else if( fCurrentMC==k14e2c )
          energy /= FunctionNL_kSDM(2.0*energy, 0.96728, -2.96279, -0.903677);

        //pass2
        else if( fCurrentMC == k15h1 )
          energy /= FunctionNL_kSDM(2.0*energy, 0.967301, -3.1683, -0.653058);

        else if( fCurrentMC == k15h2 )
          energy /= FunctionNL_kSDM(2.0*energy, 0.96728, -2.96279, -0.903677);

        else fPeriodNameAvailable = kFALSE;
      }
      break;

    // NonLinearity LHC12 ConvCalo - kTestBeamv3 + shifting MC
    case 13:
      energy *= FunctionNL_kTestBeamv3(energy);
      goto label_case_11;// goto previous case for shifting MC
      break;

    // NonLinearity LHC12 Calo - kTestBeamv3 + shifting MC
    case 14:
      energy *= FunctionNL_kTestBeamv3(energy);
      goto label_case_12;// goto previous case for shifting MC
      break;

    // NonLinearity LHC12 ConvCalo - kPi0MC + kSDM
    case 15:
      energy *= FunctionNL_kPi0MC(energy, 1.0, 0.04979, 1.3, 0.0967998, 219.381, 63.1604, 1.011);
      if(isMC == 0) energy *= FunctionNL_kSDM(energy, 0.9846, -3.319, -2.033);
      break;

    // NonLinearity LHC12 Calo - kPi0MC + kSDM
    case 16:
      energy *= FunctionNL_kPi0MC(energy, 1.0, 0.06539, 1.121, 0.0967998, 219.381, 63.1604, 1.011);
      if(isMC == 0) energy *= FunctionNL_kSDM(2.0*energy, 0.9676, -3.216, -0.6828);
      break;

//----------------------------------------------------------------------------------------------------------

    // NonLinearity LHC11a ConvCalo - only shifting MC
    case 21:
      label_case_21:
      if(isMC>0){
        if( fCurrentMC==k12f1a || fCurrentMC==k12i3 || fCurrentMC==k15g2 )
          energy /= FunctionNL_kSDM(energy, 0.984889*0.995*0.9970, -3.65456, -1.12744);

        else if(fCurrentMC==k12f1b)
          energy /= FunctionNL_kSDM(energy, 0.984384*0.995*0.9970, -3.30287, -1.48516);

        else if( fCurrentMC==k15g1a || fCurrentMC==k15g1b || fCurrentMC==k15a3a || fCurrentMC==k15a3a_plus || fCurrentMC==k15a3b )
          energy /= FunctionNL_kSDM(energy, 0.981892*0.995*0.9970, -5.43438, -1.05468);

        else fPeriodNameAvailable = kFALSE;
      }
      break;

    // NonLinearity LHC11a Calo - only shifting MC
    case 22:
      label_case_22:
      if(isMC>0){
        if(  fCurrentMC==k12f1a || fCurrentMC==k12i3 || fCurrentMC==k15g2 )
          energy /= FunctionNL_kSDM(2.0*energy, 0.966151*0.995*0.9981, -2.97974, -0.29463);

        else if( fCurrentMC==k12f1b )
          energy /= FunctionNL_kSDM(2.0*energy, 0.988814*0.995*0.9981, 0.335011, -4.30322);

        else if( fCurrentMC==k15g1a || fCurrentMC==k15g1b || fCurrentMC==k15a3a || fCurrentMC==k15a3a_plus || fCurrentMC==k15a3b )
          energy /= FunctionNL_kSDM(2.0*energy, 0.979994*0.995*0.9981, -3.24431, -0.760205);

        else fPeriodNameAvailable = kFALSE;
      }
      break;

    // NonLinearity LHC11a ConvCalo - kTestBeamv3 + shifting MC
    case 23:
      energy *= FunctionNL_kTestBeamv3(energy);
      goto label_case_21;// goto previous case for shifting MC
      break;

    // NonLinearity LHC11a Calo - kTestBeamv3 + shifting MC
    case 24:
      energy *= FunctionNL_kTestBeamv3(energy);
      goto label_case_22;// goto previous case for shifting MC
      break;

    // NonLinearity LHC11a ConvCalo - kPi0MC + kSDM
    case 25:
      energy *= FunctionNL_kPi0MC(energy, 1.0, 0.04123, 1.045, 0.0967998, 219.381, 63.1604, 1.014);
      if(isMC == 0) energy *= FunctionNL_kSDM(energy, 0.9807*0.995*0.9970, -3.377, -0.8535);
      break;

    // NonLinearity LHC11a Calo - kPi0MC + kSDM
    case 26:
      energy *= FunctionNL_kPi0MC(energy, 1.0, 0.06115, 0.9535, 0.0967998, 219.381, 63.1604, 1.013);
      if(isMC == 0) energy *= FunctionNL_kSDM(2.0*energy, 0.9772*0.995*0.9981, -3.256, -0.4449);
      break;

//----------------------------------------------------------------------------------------------------------

    // NonLinearity LHC10b-f ConvCalo - only shifting MC
    case 31:
    label_case_31:
      if(isMC>0){
        if( fCurrentMC==k14j4 )
          energy /= FunctionNL_kSDM(energy, 0.975218, -3.90409, -0.783633);

        else fPeriodNameAvailable = kFALSE;
      }
      break;

    // NonLinearity LHC10b-f Calo - only shifting MC
    case 32:
    label_case_32:
      if(isMC>0){
        if(  fCurrentMC==k14j4 )
          energy /= FunctionNL_kSDM(energy, 0.951724, -3.2596, -0.755328);

        else fPeriodNameAvailable = kFALSE;
      }
      break;

    // NonLinearity LHC10b-f ConvCalo - kTestBeamv3 + shifting MC
    case 33:
      energy *= FunctionNL_kTestBeamv3(energy);
      goto label_case_31;// goto previous case for shifting MC
      break;

    // NonLinearity LHC10b-f Calo - kTestBeamv3 + shifting MC
    case 34:
     energy *= FunctionNL_kTestBeamv3(energy);
      goto label_case_32;// goto previous case for shifting MC
      break;

//----------------------------------------------------------------------------------------------------------

    // NonLinearity LHC13 pPb ConvCalo  - only shifting MC
    case 81:
      label_case_81:
      if(isMC>0){
        if( fCurrentMC==k13b2_efix ) energy /= FunctionNL_kSDM(energy, 0.9936*0.976721, -3.60967, -0.43353);//v2
          //energy /= FunctionNL_kSDM(energy, 0.977118, -3.46238, -0.575729);v0
          //energy /= FunctionNL_kSDM(energy, 0.975357, -3.54572, -0.398501);v1
        if( fCurrentMC==k13e7 ) energy /= FunctionNL_kSDM(energy, 0.9936*0.972984, -3.4499, -0.363516);//v0
        else fPeriodNameAvailable = kFALSE;
      }
      break;

    // NonLinearity LHC13 pPb Calo  - only shifting MC
    case 82:
      label_case_82:
      if(isMC>0){
        if( fCurrentMC==k13b2_efix ) energy /= FunctionNL_kSDM(2.0*energy, 0.9970*0.974951, -2.56938, -0.863324);//v2
          //energy /= FunctionNL_kSDM(2.0*energy, 0.975467, -1.9989, -1.23208);v0
          //energy /= FunctionNL_kSDM(2.0*energy, 0.974716, -2.56403, -0.85898);v1
        if( fCurrentMC==k13e7 ) energy /= FunctionNL_kSDM(2.0*energy, 0.9970*0.975679, -3.10412, -0.606205);//v0
        else fPeriodNameAvailable = kFALSE;
      }
      break;

    // NonLinearity LHC13 pPb ConvCalo  - kTestBeamv3 + shifting MC
    case 83:
      energy *= FunctionNL_kTestBeamv3(energy);
      goto label_case_81;// goto previous case for shifting MC
      break;

    // NonLinearity LHC13 pPb Calo  - kTestBeamv3 + shifting MC
    case 84:
      energy *= FunctionNL_kTestBeamv3(energy);
      goto label_case_82;// goto previous case for shifting MC
      break;

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
AliCaloPhotonCuts::MCSet AliCaloPhotonCuts::FindEnumForMCSet(TString nameMC){
  if(nameMC.CompareTo("LHC14e2a")==0)           return k14e2a;
  else if(nameMC.CompareTo("LHC14e2b")==0)      return k14e2b;
  else if(nameMC.CompareTo("LHC14e2c")==0)      return k14e2c;
  else if(nameMC.CompareTo("LHC12f1a")==0)      return k12f1a;
  else if(nameMC.CompareTo("LHC12f1b")==0)      return k12f1b;
  else if(nameMC.CompareTo("LHC12i3")==0)       return k12i3;
  else if(nameMC.CompareTo("LHC15g1a")==0)      return k15g1a;
  else if(nameMC.CompareTo("LHC15g1b")==0)      return k15g1b;
  else if(nameMC.CompareTo("LHC15g2")==0)       return k15g2;
  else if(nameMC.CompareTo("LHC15a3a")==0)      return k15a3a;
  else if(nameMC.CompareTo("LHC15a3a_plus")==0) return k15a3a_plus;
  else if(nameMC.CompareTo("LHC15a3b")==0)      return k15a3b;
  else if(nameMC.Contains("LHC13b2_efix"))      return k13b2_efix;
  else if(nameMC.Contains("LHC13e7"))           return k13e7;
  else if(nameMC.CompareTo("LHC15h1a1")==0 ||
      nameMC.CompareTo("LHC15h1b")==0 ||
      nameMC.CompareTo("LHC15h1c")==0 ||
      nameMC.CompareTo("LHC15h1d")==0 ||
      nameMC.CompareTo("LHC15h1f")==0 ||
      nameMC.CompareTo("LHC15h1g")==0 ||
      nameMC.CompareTo("LHC15h1h")==0 ||
      nameMC.CompareTo("LHC15h1i")==0)          return k15h1;
  else if(nameMC.CompareTo("LHC15h2a")==0 ||
      nameMC.CompareTo("LHC15h2b")==0 ||
      nameMC.CompareTo("LHC15h2c")==0 ||
      nameMC.CompareTo("LHC15h2d")==0 ||
      nameMC.CompareTo("LHC15h2f")==0 ||
      nameMC.CompareTo("LHC15h2g")==0 ||
      nameMC.CompareTo("LHC15h2h")==0 ||
      nameMC.CompareTo("LHC15h2i")==0)         return k15h2;
  else if(nameMC.CompareTo("LHC14j4b")==0 ||
      nameMC.CompareTo("LHC14j4c")==0 ||
      nameMC.CompareTo("LHC14j4d")==0 ||
      nameMC.CompareTo("LHC14j4e")==0 ||
      nameMC.CompareTo("LHC14j4f")==0)         return k14j4;
  else return kNoMC;
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

