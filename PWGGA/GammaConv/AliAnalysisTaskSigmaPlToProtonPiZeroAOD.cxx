/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: Steven Merkel                                                  *
 * Supervisor: Adrian Mechler                                             *
 * Version 1.0                                                            *
 *                                                                        *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its    *
 * documentation strictly for non-commercial purposes is hereby granted    *
 * without fee, provided that the above copyright notice appears in all    *
 * copies and that both the copyright notice and this permission notice    *
 * appear in the supporting documentation. The authors make no claims    *
 * about the suitability of this software for any purpose. It is      *
 * provided "as is" without express or implied warranty.               *
 **************************************************************************/
#include "AliAODCaloCluster.h"
#include "AliPIDResponse.h"
#include "AliAnalysisTaskSigmaPlToProtonPiZeroAOD.h"
#include "AliAODVertex.h"
#include "AliVEvent.h"
#include "TList.h"
#include "TChain.h"
#include "TTree.h"
#include "TBranch.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "THn.h"
#include "THnSparse.h"
#include "TCanvas.h"
#include "TNtuple.h"
#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliAODEvent.h"
#include "AliAODInputHandler.h"
#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliMCParticle.h"
#include "AliCentrality.h"
#include "AliAODVZERO.h"
#include "AliAODPid.h"
#include "AliVParticle.h"
#include "AliAODTrack.h"
#include "AliKFVertex.h"
#include "AliV0ReaderV1.h"
#include "AliGenCocktailEventHeader.h"
#include "AliAODMCParticle.h"
#include "AliAODMCHeader.h"
#include "AliEventplane.h"
#include "AliAnalysisTaskEMCALClusterizeFast.h"
#include "AliAODEvent.h"
#include "AliAODInputHandler.h"
#include "AliAODEvent.h"
#include "AliAODInputHandler.h"
#include "AliInputEventHandler.h"
#include "AliCaloTrackMatcher.h"
#include <vector>
#include <map>
#include <fstream>
#include <iostream>
#include "TVector3.h"
#include "TLorentzVector.h"
#include "AliAODConversionPhoton.h"
#include "AliAODConversionMother.h"



class AliAnalysisTaskSigmaPlToProtonPiZeroAOD;    // your analysis class


using namespace std;            // std namespace: so you can do things like 'cout'

ClassImp(AliAnalysisTaskSigmaPlToProtonPiZeroAOD) // classimp: necessary for root


	AliAnalysisTaskSigmaPlToProtonPiZeroAOD::AliAnalysisTaskSigmaPlToProtonPiZeroAOD() : AliAnalysisTaskSE(),
	fEvent(0),
	fAOD(0),
	fMCEvent(0),
	fPIDResponse(0),
	fAODMCTrackArray(0),
	fOutputList(0),
	fAODList(NULL),
	fHistSigmaPlus(0),
	fHistSigmaPlusMC(0),
	fHistSigmaPlusMCTrueProtonGamma(0),
	fHistSigmaPlusMCTrueProton(0),
	fHistSigmaPlusMCTruePion(0),
	fHistDoubleCountTrueSigmaInvMassPt(0),
	fHistGenSigmaPt(0),
	fHistGenSigmaPerEvent(0),
	fHistGenProtonPt(0),
	fHistGenPiZeroPt(0),
	fHistSigmaPtEta(0),
	fHistProtonPtEta(0),
	fHistPi0PtEta(0),
	fHistGenAngleProtonPiZero(0),
	fHistReconstructedMassPi0(0),
	fHistReconstructedMassPi0MC(0),
	fHistReconstructedMassPi0wCut(0),
	fHistReconstructedMassPi0MCwCut(0),
	fHistPodolanski(0),
	fHistPodolanskiWCut(0),
	fHistPodolanskiWCutTrue(0),
	fHistPodolanskiGenTrue(0),
	fHistAllTracksPt(0),
	fHistProtonPt(0),
	fHistTrueProtonPt(0),
	fHistThetaPhiTrueSigmaPl(0),
	fHistThetaPhi(0),
	fHistThetaPhiProton(0),
	fHistClusterE(0),
	fHistClusterEWOCuts(0),
	fHistNClusWoCuts(0),
	fHistNClusWCuts(0),
	fHistNProtonsPerEvent(0),
	fHistDecayangle(0),
	fHistDecayangleTrue(0),
	fHistDecayanglewCut(0),
	fHistDecayangleTruewCut(0),
	fHistTrackDCAXY(0),
	fHistTrackDCAZ(0),
	fHistTrackDCAXYTrue(0),
	fHistTrackDCAZTrue(0),
	fHistTrackDCAXYwCuts(0),
	fHistTrackDCAZwCuts(0),
	fHistTrackDCAXYTruewCuts(0),
	fHistTrackDCAZTruewCuts(0),
	fHistDEDx(0),
	fHistTOFBeta(0),
	fHistTPCSignal(0),
	fHistTPCCluster(0),
    fHistTPCClusterTrue(0),
    fHistTPCchi2(0),
    fHistTPCchi2True(0),
    fHistITSCluster(0),
    fHistITSClusterTrue(0),
    fHistITSchi2(0),
    fHistITSchi2True(0),
    fHistTPCClusterwCut(0),
    fHistTPCClusterTruewCut(0),
    fHistTPCchi2wCut(0),
    fHistTPCchi2TruewCut(0),
    fHistITSClusterwCut(0),
    fHistITSClusterTruewCut(0),
    fHistITSchi2wCut(0),
    fHistITSchi2TruewCut(0),
	fHistRotationWGammaGamma(0),
	fHistRotationWProtonPion(0),
	fHistSigmaMassPtWoPodCut(0),
	fHistSigmaMassPtWoPodCutMC(0),
	fHistNLoopsProton(0),
	fHistNLoopsGamma(0),
	fHistXi0MC(0),
	fGenPhaseSpace(),
	fV0Reader(NULL),
	fV0ReaderName("V0ReaderV1"),
	fCorrTaskSetting(""),
	fEventCutArray(NULL),
	fClusterCutArray(NULL),                                   // List with Cluster Cuts
	fMesonCutArray(NULL),                                  // ConversionPhotonCutObject                                     // If a jet is near the EMCal in the current event
	fSigmaCutArray(NULL),
	fHistNEvents(NULL),                                        // array of histos with event information
	fHistNEventsWOWeight(0),                                     // array of histos with event information
	fnCuts(0),                                               // number of cuts to be analysed in parallel
	fIsHeavyIon(0),                                          // switch for pp = 0, PbPb = 1, pPb = 2
	fDoLightOutput(kFALSE),                                       // switch for running light output, kFALSE -> normal mode, kTRUE -> light mode
	fIsMC(0),
	fDoInOutTimingCluster(kFALSE),                                // manual timing cut for cluster to combine cluster within timing cut and without
	fMinTimingCluster(0),                                    // corresponding ranges, min
	fMaxTimingCluster(0),                                   // corresponding ranges, max
	fTrackMatcherRunningMode(0),
	fQTCutUpper(0),
	fQTCutLower(0),
	fWeightJetJetMC(1),
  	fOutlierJetReader(NULL)



		      // fAOD(0), fOutputList(0), fHistAllTracksPt(0), fHistProtonPt(0), fHistVertexzPos(0), fHistVertexzPosTrue(0), fHistThetaPhi(0), fHistThetaPhiProton(0)
{
	// default constructor, don't allocate memory here!
	// this is used by root for IO purposes, it needs to remain empty
}
//_____________________________________________________________________________
AliAnalysisTaskSigmaPlToProtonPiZeroAOD::AliAnalysisTaskSigmaPlToProtonPiZeroAOD(const char* name) : AliAnalysisTaskSE(name),
	fEvent(0),
	fAOD(0),
	fMCEvent(0),
	fPIDResponse(0),
	fAODMCTrackArray(0),
	fOutputList(0),
	fAODList(NULL),
	fHistSigmaPlus(0),
	fHistSigmaPlusMC(0),
	fHistSigmaPlusMCTrueProtonGamma(0),
	fHistSigmaPlusMCTrueProton(0),
	fHistSigmaPlusMCTruePion(0),
	fHistDoubleCountTrueSigmaInvMassPt(0),
	fHistGenSigmaPt(0),
	fHistGenSigmaPerEvent(0),
	fHistGenProtonPt(0),
	fHistGenPiZeroPt(0),
	fHistSigmaPtEta(0),
	fHistProtonPtEta(0),
	fHistPi0PtEta(0),
	fHistGenAngleProtonPiZero(0),
	fHistReconstructedMassPi0(0),
	fHistReconstructedMassPi0MC(0),
	fHistReconstructedMassPi0wCut(0),
	fHistReconstructedMassPi0MCwCut(0),
	fHistPodolanski(0),
	fHistPodolanskiWCut(0),
	fHistPodolanskiWCutTrue(0),
	fHistPodolanskiGenTrue(0),
	fHistAllTracksPt(0),
	fHistProtonPt(0),
	fHistTrueProtonPt(0),
	fHistThetaPhiTrueSigmaPl(0),
	fHistThetaPhi(0),
	fHistThetaPhiProton(0),
	fHistClusterE(0),
	fHistClusterEWOCuts(0),
	fHistNClusWoCuts(0),
	fHistNClusWCuts(0),
	fHistNProtonsPerEvent(0),
	fHistDecayangle(0),
	fHistDecayangleTrue(0),
	fHistDecayanglewCut(0),
	fHistDecayangleTruewCut(0),
	fHistTrackDCAXY(0),
	fHistTrackDCAZ(0),
	fHistTrackDCAXYTrue(0),
	fHistTrackDCAZTrue(0),
	fHistTrackDCAXYwCuts(0),
	fHistTrackDCAZwCuts(0),
	fHistTrackDCAXYTruewCuts(0),
	fHistTrackDCAZTruewCuts(0),
	fHistDEDx(0),
	fHistTOFBeta(0),
	fHistTPCSignal(0),
	fHistTPCCluster(0),
    fHistTPCClusterTrue(0),
    fHistTPCchi2(0),
    fHistTPCchi2True(0),
    fHistITSCluster(0),
    fHistITSClusterTrue(0),
    fHistITSchi2(0),
    fHistITSchi2True(0),
    fHistTPCClusterwCut(0),
    fHistTPCClusterTruewCut(0),
    fHistTPCchi2wCut(0),
    fHistTPCchi2TruewCut(0),
    fHistITSClusterwCut(0),
    fHistITSClusterTruewCut(0),
    fHistITSchi2wCut(0),
    fHistITSchi2TruewCut(0),
	fHistRotationWGammaGamma(0),
	fHistRotationWProtonPion(0),
	fHistSigmaMassPtWoPodCut(0),
	fHistSigmaMassPtWoPodCutMC(0),
	fHistNLoopsProton(0),
	fHistNLoopsGamma(0),
	fHistXi0MC(0),
	fGenPhaseSpace(),
	fV0Reader(NULL),
	fV0ReaderName("V0ReaderV1"),
	fCorrTaskSetting(""),
	fEventCutArray(NULL),
	fClusterCutArray(NULL),                                   // List with Cluster Cuts
	fMesonCutArray(NULL),
	fSigmaCutArray(NULL),                                // ConversionPhotonCutObject                                     // If a jet is near the EMCal in the current event
	fHistNEvents(NULL),  
	fHistNEventsWOWeight(0),                                     // array of histos with event information
	fnCuts(0),                                               // number of cuts to be analysed in parallel
	fIsHeavyIon(0),                                          // switch for pp = 0, PbPb = 1, pPb = 2
	fDoLightOutput(kFALSE),                                       // switch for running light output, kFALSE -> normal mode, kTRUE -> light mode
	fIsMC(0),
	fDoInOutTimingCluster(kFALSE),                                // manual timing cut for cluster to combine cluster within timing cut and without
	fMinTimingCluster(0),                                    // corresponding ranges, min
	fMaxTimingCluster(0),                                   // corresponding ranges, max
	fTrackMatcherRunningMode(0),
	fQTCutUpper(0),
	fQTCutLower(0),
	fWeightJetJetMC(1),
  	fOutlierJetReader(NULL)



{
	DefineOutput(1, TList::Class());
	DefineOutput(2, TTree::Class());
}
//_____________________________________________________________________________
AliAnalysisTaskSigmaPlToProtonPiZeroAOD::~AliAnalysisTaskSigmaPlToProtonPiZeroAOD()
{
	// destructor
	// if(fOutputList) {
	// 	delete fOutputList;
	// }
}
//_____________________________________________________________________________
void AliAnalysisTaskSigmaPlToProtonPiZeroAOD::UserCreateOutputObjects()
{
	fV0Reader=(AliV0ReaderV1*)AliAnalysisManager::GetAnalysisManager()->GetTask(fV0ReaderName.Data());
	if(!fV0Reader){printf("Error: No V0 Reader");return;} // GetV0Reader
	if(((AliConvEventCuts*)fEventCutArray->At(0))->GetUseJetFinderForOutliers()){
    fOutlierJetReader=(AliAnalysisTaskJetOutlierRemoval*)AliAnalysisManager::GetAnalysisManager()->GetTask("AliAnalysisTaskJetOutlierRemoval");
    if(!fOutlierJetReader){AliFatal("Error: No AliAnalysisTaskJetOutlierRemoval");} // Get jet outlier task
    else{printf("Found AliAnalysisTaskJetOutlierRemoval used for outlier removal!\n");}
  }

	// Create histograms
	if(fOutputList != NULL){
		delete fOutputList;
		fOutputList  = NULL;
	}
	if(fOutputList == NULL){
		fOutputList  = new TList();
		fOutputList->SetOwner(kTRUE);
	}
	fAODList            = new TList*[fnCuts];

	fHistSigmaPlus = new TH2F*[fnCuts];
	fHistReconstructedMassPi0 = new TH2F*[fnCuts];
	fHistReconstructedMassPi0wCut = new TH2F*[fnCuts];
	fHistPodolanski = new TH2F*[fnCuts];
	fHistPodolanskiWCut = new TH2F*[fnCuts];
	fHistPodolanskiWCutTrue = new TH2F*[fnCuts];
	fHistAllTracksPt = new TH1F*[fnCuts];
	fHistProtonPt = new TH1F*[fnCuts];
	fHistTrueProtonPt = new TH1F*[fnCuts];
	fHistThetaPhi = new TH2F*[fnCuts];
	fHistThetaPhiProton = new TH2F*[fnCuts];
	fHistClusterE = new TH1F*[fnCuts];
	fHistClusterEWOCuts = new TH1F*[fnCuts];
	fHistNClusWoCuts = new TH1F*[fnCuts];
	fHistNClusWCuts = new TH1F*[fnCuts];
	fHistNProtonsPerEvent = new TH1F*[fnCuts];
	fHistTrackDCAXY = new TH1F*[fnCuts];
	fHistTrackDCAZ = new TH1F*[fnCuts];
	fHistTrackDCAXYwCuts = new TH1F*[fnCuts];
	fHistTrackDCAZwCuts = new TH1F*[fnCuts];
	fHistDecayangle = new TH2F*[fnCuts];
	fHistDecayanglewCut = new TH2F*[fnCuts];
	fHistDEDx = new TH2F*[fnCuts];
	fHistTOFBeta = new TH2F*[fnCuts];
	fHistTPCSignal = new TH2F*[fnCuts];
	fHistTPCCluster = new TH1F*[fnCuts];
	fHistTPCchi2 = new TH1F*[fnCuts];
	fHistITSCluster = new TH1F*[fnCuts];
	fHistITSchi2 = new TH1F*[fnCuts];
	fHistTPCClusterwCut = new TH1F*[fnCuts];
	fHistTPCchi2wCut = new TH1F*[fnCuts];
	fHistITSClusterwCut = new TH1F*[fnCuts];
	fHistITSchi2wCut = new TH1F*[fnCuts];
	fHistRotationWProtonPion = new TH2F*[fnCuts];
	fHistRotationWGammaGamma = new TH2F*[fnCuts];
	fHistSigmaMassPtWoPodCut = new TH2F*[fnCuts];
	fHistNEvents = new TH1F*[fnCuts];
	fHistNEventsWOWeight = new TH1F*[fnCuts];
	
	if(fIsMC > 0){
		fHistSigmaPlusMC = new TH2F*[fnCuts];
		fHistSigmaPlusMCTrueProtonGamma = new TH2F*[fnCuts];
		fHistSigmaPlusMCTrueProton = new TH2F*[fnCuts];
		fHistSigmaPlusMCTruePion = new TH2F*[fnCuts];
		fHistDoubleCountTrueSigmaInvMassPt = new TH2F*[fnCuts];
		fHistReconstructedMassPi0MC = new TH2F*[fnCuts];
		fHistReconstructedMassPi0MCwCut = new TH2F*[fnCuts];
		fHistSigmaMassPtWoPodCutMC = new TH2F*[fnCuts];
		fHistTrackDCAXYTrue = new TH1F*[fnCuts];
		fHistTrackDCAZTrue = new TH1F*[fnCuts];
		fHistTrackDCAXYTruewCuts = new TH1F*[fnCuts];
		fHistTrackDCAZTruewCuts = new TH1F*[fnCuts];
		fHistDecayangleTrue = new TH2F*[fnCuts];
		fHistDecayangleTruewCut = new TH2F*[fnCuts];
		fHistTPCClusterTrue = new TH1F*[fnCuts];
		fHistTPCchi2True = new TH1F*[fnCuts];
		fHistITSClusterTrue = new TH1F*[fnCuts];
		fHistITSchi2True = new TH1F*[fnCuts];
		fHistTPCClusterTruewCut = new TH1F*[fnCuts];
		fHistTPCchi2TruewCut = new TH1F*[fnCuts];
		fHistITSClusterTruewCut = new TH1F*[fnCuts];
		fHistITSchi2TruewCut = new TH1F*[fnCuts];
	
		fHistThetaPhiTrueSigmaPl = new TH2F*[fnCuts];
		fHistGenSigmaPt = new TH1F*[fnCuts];
		fHistGenSigmaPerEvent = new TH1F*[fnCuts];
		fHistGenProtonPt = new TH1F*[fnCuts];
		fHistGenPiZeroPt = new TH1F*[fnCuts];
		fHistSigmaPtEta = new TH2F*[fnCuts];
		fHistProtonPtEta = new TH2F*[fnCuts];
		fHistPi0PtEta = new TH2F*[fnCuts];
		fHistGenAngleProtonPiZero = new TH1F*[fnCuts];
		fHistPodolanskiGenTrue = new TH2F*[fnCuts];
		fHistNLoopsProton = new TH1F*[fnCuts];
		fHistNLoopsGamma = new TH1F*[fnCuts];
		fHistXi0MC = new TH2F*[fnCuts];
	}

	for(Int_t iCut = 0; iCut<fnCuts;iCut++){

		TString cutstringEvent    = ((AliConvEventCuts*)fEventCutArray->At(iCut))->GetCutNumber();
		TString cutstringCalo     = ((AliCaloPhotonCuts*)fClusterCutArray->At(iCut))->GetCutNumber();
		TString cutstringMeson    = ((AliConversionMesonCuts*)fMesonCutArray->At(iCut))->GetCutNumber();
		TString cutstringSigma    = ((AliCaloSigmaCuts*)fSigmaCutArray->At(iCut))->GetCutNumber();

		fAODList[iCut]          = new TList();
		fAODList[iCut]->SetOwner(kTRUE);
		fAODList[iCut]->SetName(Form("%s_%s_%s_%s AOD histograms", cutstringEvent.Data(), cutstringCalo.Data(), cutstringMeson.Data(), cutstringSigma.Data()));
		fOutputList->Add(fAODList[iCut]);

		fHistSigmaPlus[iCut] = new TH2F("fHistSigmaPlus", ";#it{m}_{inv} (GeV/#it{c^{2}});#it{p}_{T} (GeV/#it{c})", 100, 1.1, 1.6, 40, 0., 10.);
		fAODList[iCut]->Add(fHistSigmaPlus[iCut]);
		fHistReconstructedMassPi0[iCut] = new TH2F("fHistReconstructedMassPi0",";#it{m}_{inv} (GeV/#it{c^{2}});#it{p}_{T} (GeV/#it{c})", 100, 0., 0.3, 30, 0., 15.);
		fAODList[iCut]->Add(fHistReconstructedMassPi0[iCut]);
		fHistReconstructedMassPi0wCut[iCut] = new TH2F("fHistReconstructedMassPi0wCut",";#it{m}_{inv} (GeV/#it{c^{2}});#it{p}_{T} (GeV/#it{c})", 100, 0., 0.3, 30, 0., 15.);
		fAODList[iCut]->Add(fHistReconstructedMassPi0wCut[iCut]);
		fHistPodolanski[iCut] = new TH2F("fHistPodolanski","", 100, 0., 1., 100, 0., 1.);
		fAODList[iCut]->Add(fHistPodolanski[iCut]);
		fHistPodolanskiWCut[iCut] = new TH2F("fHistPodolanskiWCut","", 100, 0., 1., 100, 0., 1.);
		fAODList[iCut]->Add(fHistPodolanskiWCut[iCut]);
		fHistAllTracksPt[iCut] = new TH1F("fHistAllTracksPt", "fHistAllTracksPt;#it{p}_{T} (GeV/#it{c});Yield", 100, 0, 10.);
		fAODList[iCut]->Add(fHistAllTracksPt[iCut]);
		fHistDecayangle[iCut] = new TH2F("fHistDecayangle", ";#vartheta (rad);#it{p}_{T} (GeV/#it{c})", 100, 0, TMath::Pi(), 40, 0., 10.);
		fAODList[iCut]->Add(fHistDecayangle[iCut]);
		fHistDecayanglewCut[iCut] = new TH2F("fHistDecayanglewCut", ";#vartheta (rad);#it{p}_{T} (GeV/#it{c})", 100, 0, TMath::Pi(), 40, 0., 10.);
		fAODList[iCut]->Add(fHistDecayanglewCut[iCut]);
		fHistProtonPt[iCut] = new TH1F("fHistProtonPt", "fHistProtonPt;#it{p}_{T,Proton} (GeV/#it{c});Yield", 100, 0, 10);
		fAODList[iCut]->Add(fHistProtonPt[iCut]);
		fHistThetaPhi[iCut] = new TH2F("fHistThetaPhi", "fHistThetaPhi;#theta ; #phi", 50, 0., TMath::Pi() ,100, 0., 2*TMath::Pi());
		fAODList[iCut]->Add(fHistThetaPhi[iCut]);
		fHistThetaPhiProton[iCut] = new TH2F("fHistThetaPhiProton", "fHistThetaPhiProton;#theta_{p} ; #phi_{p}", 50, 0., TMath::Pi() ,100, 0., 2* TMath::Pi());
		fAODList[iCut]->Add(fHistThetaPhiProton[iCut]);
		fHistClusterE[iCut] = new TH1F("fHistClusterE", "fHistClusterE;#it{E} (GeV);Yield", 100, 0, 30);
		fAODList[iCut]->Add(fHistClusterE[iCut]);
		fHistClusterEWOCuts[iCut] = new TH1F("fHistClusterEWOCuts", "fHistClusterEWOCuts;#it{E}_{Test} (GeV);Yield", 100, 0, 30);
		fAODList[iCut]->Add(fHistClusterEWOCuts[iCut]);
		fHistNClusWoCuts[iCut] = new TH1F("fHistNClusWoCuts", "fHistNClusWoCuts;#it{N}_{Cluster, wo. Cuts};Yield", 30, -0.5, 29.5);
		fAODList[iCut]->Add(fHistNClusWoCuts[iCut]);
		fHistNClusWCuts[iCut] = new TH1F("fHistNClusWCuts", "fHistNClusWCuts;#it{N}_{Cluster,w. Cuts};Yield", 30, -0.5, 29.5);
		fAODList[iCut]->Add(fHistNClusWCuts[iCut]);
		fHistNProtonsPerEvent[iCut] = new TH1F("fHistNProtonsPerEvent", "fHistNProtonsPerEvent;#it{N}_{Protons per Event};Yield", 10, -0.5, 9.5);
		fAODList[iCut]->Add(fHistNProtonsPerEvent[iCut]);
		fHistTrackDCAXY[iCut] = new TH1F("fHistTrackDCAXY", "fHistTrackDCAXY;#it{N}_{Protons per Event};", 400, 0., 0.5);
		fAODList[iCut]->Add(fHistTrackDCAXY[iCut]);
		fHistTrackDCAZ[iCut] = new TH1F("fHistTrackDCAZ", "fHistTrackDCAZ;#it{N}_{Protons per Event};", 400, 0., 0.5);
		fAODList[iCut]->Add(fHistTrackDCAZ[iCut]);
		fHistTrackDCAXYwCuts[iCut] = new TH1F("fHistTrackDCAXYwCuts", "fHistTrackDCAXYwCuts;#it{N}_{Protons per Event};", 400, 0., 0.5);
		fAODList[iCut]->Add(fHistTrackDCAXYwCuts[iCut]);
		fHistTrackDCAZwCuts[iCut] = new TH1F("fHistTrackDCAZwCuts", "fHistTrackDCAZwCuts;#it{N}_{Protons per Event};", 400, 0., 0.5);
		fAODList[iCut]->Add(fHistTrackDCAZwCuts[iCut]);
		fHistDEDx[iCut] = new TH2F("fHistDEDx", "fHistDEDx;#it{p};d#it{E}/d#it{x}", 200,0.,10.,200,1.,201.);
		fAODList[iCut]->Add(fHistDEDx[iCut]);
		fHistTOFBeta[iCut] = new TH2F("fHistTOFBeta", "fHistTOFBeta;#it{p};#beta", 200,0.,10.,130,0.1,1.3);
		fAODList[iCut]->Add(fHistTOFBeta[iCut]);
		fHistTPCSignal[iCut] = new TH2F("fHistTPCSignal", "fHistTPCSignal;#it{p};#sigma_{TPC}", 200, 0., 10., 60, -3., 3.);
		fAODList[iCut]->Add(fHistTPCSignal[iCut]);
		fHistTPCCluster[iCut] = new TH1F("fHistTPCCluster", "fHistTPCCluster;#it{N}_{Cluster TPC};", 50, 50., 100.);
		fAODList[iCut]->Add(fHistTPCCluster[iCut]);
		fHistTPCchi2[iCut] = new TH1F("fHistTPCchi2", "fHistTPCchi2;#chi^{2};", 100, 0., 10.);
		fAODList[iCut]->Add(fHistTPCchi2[iCut]);
		fHistITSCluster[iCut] = new TH1F("fHistITSCluster", "fHistITSCluster;#it{N}_{Cluster TPC};", 10, 0., 10.);
		fAODList[iCut]->Add(fHistITSCluster[iCut]);
		fHistITSchi2[iCut] = new TH1F("fHistITSchi2", "fHistITSchi2;#chi^{2};", 100, 0., 10.);
		fAODList[iCut]->Add(fHistITSchi2[iCut]);
		fHistTPCClusterwCut[iCut] = new TH1F("fHistTPCClusterwCut", "fHistTPCClusterwCut;#it{N}_{Cluster TPC};", 50, 50., 100.);
		fAODList[iCut]->Add(fHistTPCClusterwCut[iCut]);
		fHistTPCchi2wCut[iCut] = new TH1F("fHistTPCchi2wCut", "fHistTPCchi2wCut;#chi^{2};", 100, 0., 10.);
		fAODList[iCut]->Add(fHistTPCchi2wCut[iCut]);
		fHistITSClusterwCut[iCut] = new TH1F("fHistITSClusterwCut", "fHistITSClusterwCut;#it{N}_{Cluster TPC};", 10, 0., 10.);
		fAODList[iCut]->Add(fHistITSClusterwCut[iCut]);
		fHistITSchi2wCut[iCut] = new TH1F("fHistITSchi2wCut", "fHistITSchi2wCut;#chi^{2};", 100, 0., 10.);
		fAODList[iCut]->Add(fHistITSchi2wCut[iCut]);
		fHistSigmaMassPtWoPodCut[iCut] = new TH2F("fHistSigmaMassPtWoPodCut", ";#it{m}_{inv} (GeV/#it{c^{2}});#it{p}_{T} (GeV/#it{c})", 100, 1.1, 1.6, 40, 0., 10.);
		fAODList[iCut]->Add(fHistSigmaMassPtWoPodCut[iCut]);
		fHistRotationWGammaGamma[iCut] = new TH2F("fHistRotationWGammaGamma", ";#it{m}_{inv} (GeV/#it{c^{2}});#it{p}_{T} (GeV/#it{c})", 100, 1.1, 1.6, 40, 0., 10.);
		fHistRotationWGammaGamma[iCut]->Sumw2();
		fAODList[iCut]->Add(fHistRotationWGammaGamma[iCut]);
		fHistRotationWProtonPion[iCut] = new TH2F("fHistRotationWProtonPion", ";#it{m}_{inv} (GeV/#it{c^{2}});#it{p}_{T} (GeV/#it{c})", 100, 1.1, 1.6, 40, 0., 10.);
		fHistRotationWProtonPion[iCut]->Sumw2();
		fAODList[iCut]->Add(fHistRotationWProtonPion[iCut]);
		
		if(fIsMC > 1) {
			fHistSigmaPlus[iCut]->Sumw2();
			fHistReconstructedMassPi0[iCut]->Sumw2();
			fHistReconstructedMassPi0wCut[iCut]->Sumw2();
			fHistPodolanski[iCut]->Sumw2();
			fHistPodolanskiWCut[iCut]->Sumw2();
			fHistAllTracksPt[iCut]->Sumw2();
			fHistDecayangle[iCut]->Sumw2();
			fHistDecayanglewCut[iCut]->Sumw2();
			fHistProtonPt[iCut]->Sumw2();
			fHistThetaPhi[iCut]->Sumw2();
			fHistThetaPhiProton[iCut]->Sumw2();
			fHistClusterE[iCut]->Sumw2();
			fHistClusterEWOCuts[iCut]->Sumw2();
			fHistNClusWoCuts[iCut]->Sumw2();
			fHistNClusWCuts[iCut]->Sumw2();
			fHistNProtonsPerEvent[iCut]->Sumw2();
			fHistTrackDCAXY[iCut]->Sumw2();
			fHistTrackDCAZ[iCut]->Sumw2();
			fHistTrackDCAXYwCuts[iCut]->Sumw2();
			fHistTrackDCAZwCuts[iCut]->Sumw2();
			fHistDEDx[iCut]->Sumw2();
			fHistTOFBeta[iCut]->Sumw2();
			fHistTPCSignal[iCut]->Sumw2();
			fHistTPCCluster[iCut]->Sumw2();
			fHistTPCchi2[iCut]->Sumw2();
			fHistITSCluster[iCut]->Sumw2();
			fHistITSchi2[iCut]->Sumw2();
			fHistTPCClusterwCut[iCut]->Sumw2();
			fHistTPCchi2wCut[iCut]->Sumw2();
			fHistITSClusterwCut[iCut]->Sumw2();
			fHistITSchi2wCut[iCut]->Sumw2();
			fHistSigmaMassPtWoPodCut[iCut]->Sumw2();
		}


		fHistNEvents[iCut]     = new TH1F("NEvents", "NEvents", 15, -0.5, 14.5);
		fHistNEvents[iCut]->GetXaxis()->SetBinLabel(1,"Accepted");
		fHistNEvents[iCut]->GetXaxis()->SetBinLabel(2,"Centrality");
		fHistNEvents[iCut]->GetXaxis()->SetBinLabel(3,"Miss. MC or inc. ev.");
		fHistNEvents[iCut]->GetXaxis()->SetBinLabel(4,"Trigger");
		fHistNEvents[iCut]->GetXaxis()->SetBinLabel(5,"Vertex Z");
		fHistNEvents[iCut]->GetXaxis()->SetBinLabel(6,"Cont. Vertex");
		fHistNEvents[iCut]->GetXaxis()->SetBinLabel(7,"SPD Pile-Up");
		fHistNEvents[iCut]->GetXaxis()->SetBinLabel(8,"no SDD");
		fHistNEvents[iCut]->GetXaxis()->SetBinLabel(9,"no V0AND");
		fHistNEvents[iCut]->GetXaxis()->SetBinLabel(10,"EMCAL problems");
		fHistNEvents[iCut]->GetXaxis()->SetBinLabel(11,"rejectedForJetJetMC");
		fHistNEvents[iCut]->GetXaxis()->SetBinLabel(12,"SPD hits vs tracklet");
		fHistNEvents[iCut]->GetXaxis()->SetBinLabel(13,"Out-of-Bunch pileup Past-Future");
		fHistNEvents[iCut]->GetXaxis()->SetBinLabel(14,"Pileup V0M-TPCout Tracks");
		fHistNEvents[iCut]->GetXaxis()->SetBinLabel(15,"Sphericity");
		fHistNEvents[iCut]->Sumw2();

		fHistNEventsWOWeight[iCut]     = new TH1F("fHistNEventsWOWeight", "fHistNEventsWOWeight", 15, -0.5, 14.5);
		fHistNEventsWOWeight[iCut]->GetXaxis()->SetBinLabel(1,"Accepted");
		fHistNEventsWOWeight[iCut]->GetXaxis()->SetBinLabel(2,"Centrality");
		fHistNEventsWOWeight[iCut]->GetXaxis()->SetBinLabel(3,"Miss. MC or inc. ev.");
		fHistNEventsWOWeight[iCut]->GetXaxis()->SetBinLabel(4,"Trigger");
		fHistNEventsWOWeight[iCut]->GetXaxis()->SetBinLabel(5,"Vertex Z");
		fHistNEventsWOWeight[iCut]->GetXaxis()->SetBinLabel(6,"Cont. Vertex");
		fHistNEventsWOWeight[iCut]->GetXaxis()->SetBinLabel(7,"SPD Pile-Up");
		fHistNEventsWOWeight[iCut]->GetXaxis()->SetBinLabel(8,"no SDD");
		fHistNEventsWOWeight[iCut]->GetXaxis()->SetBinLabel(9,"no V0AND");
		fHistNEventsWOWeight[iCut]->GetXaxis()->SetBinLabel(10,"EMCAL problems");
		fHistNEventsWOWeight[iCut]->GetXaxis()->SetBinLabel(11,"rejectedForJetJetMC");
		fHistNEventsWOWeight[iCut]->GetXaxis()->SetBinLabel(12,"SPD hits vs tracklet");
		fHistNEventsWOWeight[iCut]->GetXaxis()->SetBinLabel(13,"Out-of-Bunch pileup Past-Future");
		fHistNEventsWOWeight[iCut]->GetXaxis()->SetBinLabel(14,"Pileup V0M-TPCout Tracks");
		fHistNEventsWOWeight[iCut]->GetXaxis()->SetBinLabel(15,"Sphericity");

		fAODList[iCut]->Add(fHistNEvents[iCut]);          // don't forget to add it to the list! the list will be written to file, so if you want
		fAODList[iCut]->Add(fHistNEventsWOWeight[iCut]);          // don't forget to add it to the list! the list will be written to file, so if you want

		if(!((AliConvEventCuts*)fEventCutArray->At(iCut))) continue;
		if(((AliConvEventCuts*)fEventCutArray->At(iCut))->GetCutHistograms()){
			fAODList[iCut]->Add(((AliConvEventCuts*)fEventCutArray->At(iCut))->GetCutHistograms());
		}
		if(!((AliCaloPhotonCuts*)fClusterCutArray->At(iCut))) continue;
		if(((AliCaloPhotonCuts*)fClusterCutArray->At(iCut))->GetCutHistograms()){
			fAODList[iCut]->Add(((AliCaloPhotonCuts*)fClusterCutArray->At(iCut))->GetCutHistograms());
		}
		if(!((AliConversionMesonCuts*)fMesonCutArray->At(iCut))) continue;
		if(((AliConversionMesonCuts*)fMesonCutArray->At(iCut))->GetCutHistograms()){
			fAODList[iCut]->Add(((AliConversionMesonCuts*)fMesonCutArray->At(iCut))->GetCutHistograms());
		}
		if(!((AliCaloSigmaCuts*)fSigmaCutArray->At(iCut))) continue;
    	if(((AliCaloSigmaCuts*)fSigmaCutArray->At(iCut))->GetCutHistograms()){
      		fAODList[iCut]->Add(((AliCaloSigmaCuts*)fSigmaCutArray->At(iCut))->GetCutHistograms());
    	}

		if(fIsMC > 0){
			fHistPodolanskiWCutTrue[iCut] = new TH2F("fHistPodolanskiWCutTrue","", 100, 0., 1., 100, 0., 1.);
			fHistTrueProtonPt[iCut] = new TH1F("fHistTrueProtonPt", "fHistTrueProtonPt;#it{p}_{T,Proton} (GeV/#it{c});Yield", 100, 0, 10);
			fHistSigmaMassPtWoPodCutMC[iCut] = new TH2F("fHistSigmaMassPtWoPodCutMC", ";#it{m}_{inv} (GeV/#it{c^{2}});#it{p}_{T} (GeV/#it{c})", 100, 1.1, 1.6, 40, 0., 10.);
			fHistTrackDCAXYTrue[iCut] = new TH1F("fHistTrackDCAXYTrue", "fHistTrackDCAXYTrue;#it{N}_{Protons per Event};", 400, 0., 0.5);
			fHistTrackDCAZTrue[iCut] = new TH1F("fHistTrackDCAZTrue", "fHistTrackDCAZTrue;#it{N}_{Protons per Event};", 400, 0., 0.5);
			fHistTrackDCAXYTruewCuts[iCut] = new TH1F("fHistTrackDCAXYTruewCuts", "fHistTrackDCAXYTruewCuts;#it{N}_{Protons per Event};", 400, 0., 0.5);
			fHistTrackDCAZTruewCuts[iCut] = new TH1F("fHistTrackDCAZTruewCuts", "fHistTrackDCAZTruewCuts;#it{N}_{Protons per Event};", 400, 0., 0.5);
			fHistDecayangleTrue[iCut] = new TH2F("fHistDecayangleTrue", ";#vartheta (rad);#it{p}_{T,Proton} (GeV/#it{c})", 100, 0, TMath::Pi(), 40, 0., 10.);
			fHistDecayangleTruewCut[iCut] = new TH2F("fHistDecayangleTruewCut", ";#vartheta (rad);#it{p}_{T,Proton} (GeV/#it{c})", 100, 0, TMath::Pi(), 40, 0., 10.);
			fHistSigmaPlusMC[iCut] = new TH2F("fHistSigmaPlusMC", ";#it{m}_{inv} (GeV/#it{c^{2}});#it{p}_{T} (GeV/#it{c})", 100, 1.1, 1.6, 40, 0., 10.);
			fHistSigmaPlusMCTrueProtonGamma[iCut] = new TH2F("fHistSigmaPlusMCTrueProtonGamma", ";#it{m}_{inv} (GeV/#it{c^{2}});#it{p}_{T} (GeV/#it{c})", 100, 1.1, 1.6, 40, 0., 10.);
			fHistSigmaPlusMCTrueProton[iCut] = new TH2F("fHistSigmaPlusMCTrueProton", ";#it{m}_{inv} (GeV/#it{c^{2}});#it{p}_{T} (GeV/#it{c})", 100, 1.1, 1.6, 40, 0., 10.);
			fHistSigmaPlusMCTruePion[iCut] = new TH2F("fHistSigmaPlusMCTruePion", ";#it{m}_{inv} (GeV/#it{c^{2}});#it{p}_{T} (GeV/#it{c})", 100, 1.1, 1.6, 40, 0., 10.);
			fHistDoubleCountTrueSigmaInvMassPt[iCut] = new TH2F("fHistDoubleCountTrueSigmaInvMassPt", ";#it{m}_{inv} (GeV/#it{c^{2}});#it{p}_{T} (GeV/#it{c})", 100, 1.1, 1.6, 40, 0., 10.);
			fHistReconstructedMassPi0MC[iCut] = new TH2F("fHistReconstructedMassPi0MC",";#it{m}_{inv} (GeV/#it{c^{2}});#it{p}_{T} (GeV/#it{c})", 100, 0., 0.3, 30, 0., 15.);
			fHistReconstructedMassPi0MCwCut[iCut] = new TH2F("fHistReconstructedMassPi0MCwCut",";#it{m}_{inv} (GeV/#it{c^{2}});#it{p}_{T} (GeV/#it{c})", 100, 0., 0.3, 30, 0., 15.);
			fHistTPCClusterTrue[iCut] = new TH1F("fHistTPCClusterTrue", "fHistTPCClusterTrue;#it{N}_{Cluster TPC};", 50, 50., 100.);
			fHistTPCchi2True[iCut] = new TH1F("fHistTPCchi2True", "fHistTPCchi2True;#chi^{2};", 100, 0., 10.);
			fHistITSClusterTrue[iCut] = new TH1F("fHistITSClusterTrue", "fHistITSClusterTrue;#it{N}_{Cluster TPC};", 10, 0., 10.);
			fHistITSchi2True[iCut] = new TH1F("fHistITSchi2True", "fHistITSchi2True;#chi^{2};", 100, 0., 10.);
			fHistTPCClusterTruewCut[iCut] = new TH1F("fHistTPCClusterTruewCut", "fHistTPCClusterTruewCut;#it{N}_{Cluster TPC};", 50, 50., 100.);
			fHistTPCchi2TruewCut[iCut] = new TH1F("fHistTPCchi2TruewCut", "fHistTPCchi2TruewCut;#chi^{2};", 100, 0., 10.);
			fHistITSClusterTruewCut[iCut] = new TH1F("fHistITSClusterTruewCut", "fHistITSClusterTruewCut;#it{N}_{Cluster TPC};", 10, 0., 10.);
			fHistITSchi2TruewCut[iCut] = new TH1F("fHistITSchi2TruewCut", "fHistITSchi2TruewCut;#chi^{2};", 100, 0., 10.);
			fHistThetaPhiTrueSigmaPl[iCut] = new TH2F("fHistThetaPhiTrueSigmaPl", "fHistThetaPhiTrueSigmaPl;#theta ; #phi", 50, -1., 4. ,100, 0., 2*TMath::Pi());
			fHistGenSigmaPt[iCut] = new TH1F("fHistGenSigmaPt", "fHistGenSigmaPt;#it{p}_{T,Sigma} (GeV/#it{c});Yield", 100, 0, 30);
			fHistGenSigmaPerEvent[iCut] = new TH1F("fHistGenSigmaPerEvent", "fHistGenSigmaPerEvent;N;Yield", 11, -0.5, 10.5);
			fHistGenSigmaPerEvent[iCut] -> Sumw2();
			fHistGenProtonPt[iCut] = new TH1F("fHistGenProtonPt", "fHistGenProtonPt;#it{p}_{T,Proton} (GeV/#it{c});Yield", 100, 0, 30);
			fHistGenPiZeroPt[iCut] = new TH1F("fHistGenPiZeroPt", "fHistGenPiZeroPt;#it{p}_{T,Pion} (GeV/#it{c});Yield", 100, 0, 30);
			fHistSigmaPtEta[iCut] = new TH2F("fHistSigmaPtEta", "fHistSigmaPtEta;#it{p}_{T,Pion} (GeV/#it{c});Yield", 100, 0, 30, 50, -5, 5);
			fHistProtonPtEta[iCut] = new TH2F("fHistProtonPtEta", "fHistProtonPtEta;#it{p}_{T,Pion} (GeV/#it{c});Yield", 100, 0, 30, 50, -5, 5);
			fHistPi0PtEta[iCut] = new TH2F("fHistPi0PtEta", "fHistPi0PtEta;#it{p}_{T,Pion} (GeV/#it{c});Yield", 100, 0, 10, 50, -5, 5);
			fHistGenAngleProtonPiZero[iCut] = new TH1F("fHistGenAngleProtonPiZero", "fHistGenAngleProtonPiZero;#it{beta}_{Proton,Pion} (rad);Yield", 100, 0, TMath::Pi());
			fHistPodolanskiGenTrue[iCut] = new TH2F("fHistPodolanskiGenTrue","", 100, 0., 1., 100, 0., 1.);
			fHistNLoopsProton[iCut] = new TH1F("fHistNLoopsProton", "fHistNLoopsProton;#chi^{2};", 100, 0., 100.);
			fHistNLoopsGamma[iCut] = new TH1F("fHistNLoopsGamma", "fHistNLoopsGamma;#chi^{2};", 100, 0., 100.);
			fHistXi0MC[iCut] = new TH2F("fHistXi0MC", ";#it{m}_{inv} (GeV/#it{c^{2}});#it{p}_{T} (GeV/#it{c})", 100, 1.1, 1.6, 40, 0., 10.);
			if(fIsMC > 1){
				fHistPodolanskiWCutTrue[iCut]->Sumw2();
				fHistTrueProtonPt[iCut]->Sumw2();
				fHistSigmaMassPtWoPodCutMC[iCut]->Sumw2();
				fHistTrackDCAXYTrue[iCut]->Sumw2();
				fHistTrackDCAZTrue[iCut]->Sumw2();
				fHistTrackDCAXYTruewCuts[iCut]->Sumw2();
				fHistTrackDCAZTruewCuts[iCut]->Sumw2();
				fHistDecayangleTrue[iCut]->Sumw2();
				fHistDecayangleTruewCut[iCut]->Sumw2();
				fHistReconstructedMassPi0MCwCut[iCut]->Sumw2();
				fHistReconstructedMassPi0MC[iCut]->Sumw2();
				fHistDoubleCountTrueSigmaInvMassPt[iCut]->Sumw2();
				fHistSigmaPlusMCTruePion[iCut]->Sumw2();
				fHistSigmaPlusMCTrueProton[iCut]->Sumw2();
				fHistSigmaPlusMCTrueProtonGamma[iCut]->Sumw2();
				fHistSigmaPlusMC[iCut]->Sumw2();
				fHistThetaPhiTrueSigmaPl[iCut]->Sumw2();
				fHistGenSigmaPt[iCut]->Sumw2();
				fHistGenProtonPt[iCut]->Sumw2();
				fHistGenPiZeroPt[iCut]->Sumw2();
				fHistSigmaPtEta[iCut]->Sumw2();
				fHistProtonPtEta[iCut]->Sumw2();
				fHistPi0PtEta[iCut]->Sumw2();
				fHistGenAngleProtonPiZero[iCut]->Sumw2();
				fHistPodolanskiGenTrue[iCut]->Sumw2();
				fHistTPCClusterTrue[iCut]->Sumw2();
				fHistTPCchi2True[iCut]->Sumw2();
				fHistITSClusterTrue[iCut]->Sumw2();
				fHistITSchi2True[iCut]->Sumw2();
				fHistTPCClusterTruewCut[iCut]->Sumw2();
				fHistTPCchi2TruewCut[iCut]->Sumw2();
				fHistITSClusterTruewCut[iCut]->Sumw2();
				fHistITSchi2TruewCut[iCut]->Sumw2();
				fHistNLoopsProton[iCut]->Sumw2();
				fHistNLoopsGamma[iCut]->Sumw2();
				fHistXi0MC[iCut]->Sumw2();
			}
			
			fAODList[iCut]->Add(fHistSigmaPlusMC[iCut]);
			fAODList[iCut]->Add(fHistSigmaPlusMCTrueProtonGamma[iCut]);
			fAODList[iCut]->Add(fHistSigmaPlusMCTrueProton[iCut]);
			fAODList[iCut]->Add(fHistSigmaPlusMCTruePion[iCut]);
			fAODList[iCut]->Add(fHistDoubleCountTrueSigmaInvMassPt[iCut]);
			fAODList[iCut]->Add(fHistPodolanskiGenTrue[iCut]);
			fAODList[iCut]->Add(fHistPodolanskiWCutTrue[iCut]);
			fAODList[iCut]->Add(fHistTrueProtonPt[iCut]);
			fAODList[iCut]->Add(fHistTrackDCAXYTrue[iCut]);
			fAODList[iCut]->Add(fHistTrackDCAZTrue[iCut]);
			fAODList[iCut]->Add(fHistTrackDCAXYTruewCuts[iCut]);
			fAODList[iCut]->Add(fHistTrackDCAZTruewCuts[iCut]);
			fAODList[iCut]->Add(fHistDecayangleTrue[iCut]);
			fAODList[iCut]->Add(fHistDecayangleTruewCut[iCut]);
			fAODList[iCut]->Add(fHistTPCClusterTrue[iCut]);
			fAODList[iCut]->Add(fHistTPCchi2True[iCut]);
			fAODList[iCut]->Add(fHistITSClusterTrue[iCut]);
			fAODList[iCut]->Add(fHistITSchi2True[iCut]);
			fAODList[iCut]->Add(fHistTPCClusterTruewCut[iCut]);
			fAODList[iCut]->Add(fHistTPCchi2TruewCut[iCut]);
			fAODList[iCut]->Add(fHistITSClusterTruewCut[iCut]);
			fAODList[iCut]->Add(fHistITSchi2TruewCut[iCut]);
			fAODList[iCut]->Add(fHistThetaPhiTrueSigmaPl[iCut]);
			fAODList[iCut]->Add(fHistGenSigmaPerEvent[iCut]);
			fAODList[iCut]->Add(fHistGenSigmaPt[iCut]);
			fAODList[iCut]->Add(fHistGenProtonPt[iCut]);
			fAODList[iCut]->Add(fHistGenPiZeroPt[iCut]);
			fAODList[iCut]->Add(fHistSigmaPtEta[iCut]);
			fAODList[iCut]->Add(fHistProtonPtEta[iCut]);
			fAODList[iCut]->Add(fHistPi0PtEta[iCut]);
			fAODList[iCut]->Add(fHistGenAngleProtonPiZero[iCut]);
			fAODList[iCut]->Add(fHistReconstructedMassPi0MC[iCut]);
			fAODList[iCut]->Add(fHistReconstructedMassPi0MCwCut[iCut]);
			fAODList[iCut]->Add(fHistSigmaMassPtWoPodCutMC[iCut]);
			fAODList[iCut]->Add(fHistNLoopsProton[iCut]);
			fAODList[iCut]->Add(fHistNLoopsGamma[iCut]);
			fAODList[iCut]->Add(fHistXi0MC[iCut]);
		}
	}
	if(fV0Reader)
		if((AliConvEventCuts*)fV0Reader->GetEventCuts())
			if(((AliConvEventCuts*)fV0Reader->GetEventCuts())->GetCutHistograms())
				fOutputList->Add(((AliConvEventCuts*)fV0Reader->GetEventCuts())->GetCutHistograms());
	PostData(1, fOutputList);
}
//_____________________________________________________________________________
void AliAnalysisTaskSigmaPlToProtonPiZeroAOD::UserExec(Option_t *)
{

	fEvent = InputEvent();
	if (!fEvent) {
		Printf("ERROR: Could not retrieve event");
		return;
	}
	fAOD = dynamic_cast<AliAODEvent*>(fEvent);
	if (!fAOD) {
		Printf("ERROR: Could not retrieve fAOD");
		return;
	}
	//Eventqualitäts check
	if(fIsMC> 0){
		// fEventHandler = dynamic_cast<AliAODInputHandler*> (AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
		// fMCEvent = fEventHandler->MCEvent();
		fMCEvent = MCEvent();
		if(!fMCEvent)  {
			Printf("ERROR: Could not retrieve fMCEvent");
			return;
		}
	}

	Int_t eventQuality = ((AliConvEventCuts*)fV0Reader->GetEventCuts())->GetEventQuality();
	if(fEvent->IsIncompleteDAQ()==kTRUE) eventQuality = 2;  // incomplete event
	if(eventQuality == 2 || eventQuality == 3){// Event Not Accepted due to MC event missing or wrong trigger for V0ReaderV1 or because it is incomplete
		for(Int_t iCut = 0; iCut<fnCuts; iCut++){
			FillfHistNEvents(iCut,eventQuality, fWeightJetJetMC);
			if (fIsMC>1) fHistNEventsWOWeight[iCut]->Fill(eventQuality);
		}
		return;
	}

	//Find EventVertex
	const AliAODVertex *vertex=fAOD->GetPrimaryVertex();    // 22/8
	if(!vertex) {
		Printf("ERROR: Could not retrieve vertex");
		return;
	}
	//Primär Vertex Bestimmen
	Double_t vpos[3];
	vertex->GetXYZ(vpos);
	//Find Protons and secundary Vertex
	Int_t iTracks=0;
	if(!fAOD->GetNumberOfTracks()) {
		// Printf("ERROR: Could not retrieve fAOD->GetNumberOfTracks()");
		return;
	}
	iTracks=fAOD->GetNumberOfTracks();  // see how many tracks there are in the event

	if(!fPIDResponse) {
		AliAnalysisManager *man = AliAnalysisManager::GetAnalysisManager();
		if(man){
			AliInputEventHandler* inputHandler = (AliInputEventHandler*)(man->GetInputEventHandler());
			if(!inputHandler){
				Printf("ERROR: Could not retrieve inputHandler");
				return;
			}
			fPIDResponse = (AliPIDResponse*)inputHandler->GetPIDResponse();
			if (!fPIDResponse) {
				Printf("ERROR: Could not retrieve fPIDResponse");
				return;
			}
		} else {
			Printf("ERROR: Could not retrieve AliAnalysisManager");
			return;
		}
	}

	for(Int_t iCut = 0; iCut<fnCuts; iCut++){
		//Vektoren für die Einzelnen Teilchen und Vertices initialisieren
		vector < AliAODTrack* > proton;
		vector < AliAODTrack* > tracks;
		vector < AliVCluster* > photon;
		vector < TLorentzVector > pions;

		Bool_t isRunningEMCALrelAna = kFALSE;
		if (((AliCaloPhotonCuts*)fClusterCutArray->At(iCut))->GetClusterType() == 1
				|| ((AliCaloPhotonCuts*)fClusterCutArray->At(iCut))->GetClusterType() == 3
				|| ((AliCaloPhotonCuts*)fClusterCutArray->At(iCut))->GetClusterType() == 4
		   ) isRunningEMCALrelAna = kTRUE;


		if(fIsMC>0){
	      fWeightJetJetMC       = 1;
	  //     cout << fMCEvent << endl;
	      Float_t maxjetpt      = -1.;
	      Float_t pthard = -1;
	      if(((AliConvEventCuts*)fEventCutArray->At(iCut))->GetUseJetFinderForOutliers()) maxjetpt = fOutlierJetReader->GetMaxJetPt();
	      Bool_t isMCJet        = ((AliConvEventCuts*)fEventCutArray->At(iCut))->IsJetJetMCEventAccepted( fMCEvent, fWeightJetJetMC, pthard,fEvent, maxjetpt);
	      if(fIsMC==1) fWeightJetJetMC = 1;
	      if (!isMCJet){
	        fHistNEvents[iCut]->Fill(10,fWeightJetJetMC);
	        if (fIsMC>1) fHistNEventsWOWeight[iCut]->Fill(10);
	        continue;
	      }
	    }
	

		Int_t eventNotAccepted = ((AliConvEventCuts*)fEventCutArray->At(iCut))->IsEventAcceptedByCut(fV0Reader->GetEventCuts(),fEvent,fMCEvent,fIsHeavyIon, isRunningEMCALrelAna);

		if(eventNotAccepted != 0){ // Check Centrality, PileUp, SDD and V0AND --> Not Accepted => eventQuality = 1
			if (eventNotAccepted==3 && fIsMC == 1){
				FillfHistNEvents(iCut,0,fWeightJetJetMC);
				if (fIsMC>1) fHistNEventsWOWeight[iCut]->Fill(0);
			} else {
				FillfHistNEvents(iCut,eventNotAccepted,fWeightJetJetMC);
				if (fIsMC>1) fHistNEventsWOWeight[iCut]->Fill(eventNotAccepted);
				continue;
			}
		} else {
			FillfHistNEvents(iCut,eventNotAccepted,fWeightJetJetMC);
			if (fIsMC>1) fHistNEventsWOWeight[iCut]->Fill(eventNotAccepted);
		}

		// Generiertes Spektrum der Sigma+
		if(fIsMC > 0){
			fAODMCTrackArray = dynamic_cast<TClonesArray*>(fInputEvent->FindListObject(AliAODMCParticle::StdBranchName()));
			if (fAODMCTrackArray == NULL) return;
			Int_t nSigmaperEvent =0;
			for(Long_t i = 1; i < fAODMCTrackArray->GetEntriesFast(); ++i) {
				AliAODMCParticle* sigma = nullptr;
				sigma = static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(i));
				if (!sigma) {
					Printf("ERROR: Could not retrieve AliAODMCParticle");
					continue;
				}
				if(sigma->GetPdgCode() == 3222){
					if(fHistThetaPhiTrueSigmaPl[iCut]) fHistThetaPhiTrueSigmaPl[iCut]->Fill(sigma->Theta(), sigma->Phi(),fWeightJetJetMC);
					Int_t nDaughters = sigma->GetNDaughters();
					TLorentzVector protonTrueVec;
					TLorentzVector pionTrueVec;
					TLorentzVector sigmaTrueVec;
					if(nDaughters == 2){
						if((sigma->GetDaughterLabel(0) < 0) || (sigma->GetDaughterLabel(1) < 0)){
							continue;
						}
						sigmaTrueVec.SetPtEtaPhiM(sigma->Pt(), sigma->Eta(), sigma->Phi(), sigma->M());
						AliAODMCParticle* daughter1 = static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(sigma->GetDaughterLabel(0)));
						AliAODMCParticle* daughter2 = static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(sigma->GetDaughterLabel(1)));
						if((daughter1->GetPdgCode() == 2212) && (daughter2->GetPdgCode() == 111)){
							protonTrueVec.SetPtEtaPhiM(daughter1->Pt(), daughter1->Eta(), daughter1->Phi(), daughter1->M());
							pionTrueVec.SetPtEtaPhiM(daughter2->Pt(), daughter2->Eta(), daughter2->Phi(), daughter2->M());
							if(fHistGenSigmaPt[iCut]) fHistGenSigmaPt[iCut]->Fill(sigma->Pt(),fWeightJetJetMC);
							if(fHistGenProtonPt[iCut]) fHistGenProtonPt[iCut]->Fill(protonTrueVec.Pt(),fWeightJetJetMC);
							if(fHistGenPiZeroPt[iCut]) fHistGenPiZeroPt[iCut]->Fill(pionTrueVec.Pt(),fWeightJetJetMC);
							if(fHistSigmaPtEta[iCut]) fHistSigmaPtEta[iCut]->Fill(sigma->Pt(), sigma->Eta(),fWeightJetJetMC);
							if(fHistProtonPtEta[iCut]) fHistProtonPtEta[iCut]->Fill(protonTrueVec.Pt(), protonTrueVec.Eta(),fWeightJetJetMC);
							if(fHistPi0PtEta[iCut]) fHistPi0PtEta[iCut]->Fill(pionTrueVec.Pt(), pionTrueVec.Eta(),fWeightJetJetMC);
							if(fHistGenAngleProtonPiZero[iCut]) fHistGenAngleProtonPiZero[iCut]->Fill(pionTrueVec.Angle(protonTrueVec.Vect()),fWeightJetJetMC);
							if(fHistPodolanskiGenTrue[iCut]) fHistPodolanskiGenTrue[iCut]->Fill(GetPodAlpha(sigmaTrueVec, protonTrueVec, pionTrueVec),GetQT(sigmaTrueVec, pionTrueVec),fWeightJetJetMC);
							nSigmaperEvent = nSigmaperEvent+1;
						}
						if((daughter1->GetPdgCode() == 111) && (daughter2->GetPdgCode() == 2212)){
							pionTrueVec.SetPtEtaPhiM(daughter1->Pt(), daughter1->Eta(), daughter1->Phi(), daughter1->M());
							protonTrueVec.SetPtEtaPhiM(daughter2->Pt(), daughter2->Eta(), daughter2->Phi(), daughter2->M());
							if(fHistGenSigmaPt[iCut]) fHistGenSigmaPt[iCut]->Fill(sigma->Pt(),fWeightJetJetMC);
							if(fHistGenProtonPt[iCut]) fHistGenProtonPt[iCut]->Fill(protonTrueVec.Pt(),fWeightJetJetMC);
							if(fHistGenPiZeroPt[iCut]) fHistGenPiZeroPt[iCut]->Fill(pionTrueVec.Pt(),fWeightJetJetMC);
							if(fHistSigmaPtEta[iCut]) fHistSigmaPtEta[iCut]->Fill(sigma->Pt(), sigma->Eta(),fWeightJetJetMC);
							if(fHistProtonPtEta[iCut]) fHistProtonPtEta[iCut]->Fill(protonTrueVec.Pt(), protonTrueVec.Eta(),fWeightJetJetMC);
							if(fHistPi0PtEta[iCut]) fHistPi0PtEta[iCut]->Fill(pionTrueVec.Pt(), pionTrueVec.Eta(),fWeightJetJetMC);
							if(fHistGenAngleProtonPiZero[iCut]) fHistGenAngleProtonPiZero[iCut]->Fill(pionTrueVec.Angle(protonTrueVec.Vect()),fWeightJetJetMC);
							if(fHistPodolanskiGenTrue[iCut]) fHistPodolanskiGenTrue[iCut]->Fill(GetPodAlpha(sigmaTrueVec, protonTrueVec, pionTrueVec),GetQT(sigmaTrueVec, pionTrueVec),fWeightJetJetMC);
							nSigmaperEvent = nSigmaperEvent+1;
						}
					}
				}
			}
			if(fHistGenSigmaPerEvent[iCut]) fHistGenSigmaPerEvent[iCut]->Fill(nSigmaperEvent,fWeightJetJetMC);
		}

		AliAODTrack* track;
		for(Int_t i = 0; i < iTracks; i++) {
			AliAODTrack* tracktmp;
			tracktmp = static_cast<AliAODTrack*>(fAOD->GetTrack(i));
			if(!tracktmp){
				Printf("ERROR: Could not retrieve track %i", i);
				continue;
			}
			track = new AliAODTrack(*tracktmp);
			if(!track) continue;
			Float_t trackDCAXY = 0.0, trackDCAZ = 0.0;
			track->GetImpactParameters(trackDCAXY,trackDCAZ);
			if(fHistTPCCluster[iCut] && track->GetTPCNcls())fHistTPCCluster[iCut]->Fill(track->GetTPCNcls(),fWeightJetJetMC);
			if(fHistTPCchi2[iCut] && track->GetTPCchi2perCluster())fHistTPCchi2[iCut]->Fill(track->GetTPCchi2perCluster(),fWeightJetJetMC);
			if(fHistITSCluster[iCut] && track->GetITSNcls())fHistITSCluster[iCut]->Fill(track->GetITSNcls(),fWeightJetJetMC);
			if(fHistITSchi2[iCut] && track->GetITSchi2())fHistITSchi2[iCut]->Fill(track->GetITSchi2(),fWeightJetJetMC);
			if(fIsMC > 0){
				if((IsRealProton(track, fAODMCTrackArray, iCut, fWeightJetJetMC, 0)) > 0){
					if(fHistTPCClusterTrue[iCut] && track->GetTPCNcls())fHistTPCClusterTrue[iCut]->Fill(track->GetTPCNcls(),fWeightJetJetMC);
					if(fHistTPCchi2True[iCut] && track->GetTPCchi2perCluster())fHistTPCchi2True[iCut]->Fill(track->GetTPCchi2perCluster(),fWeightJetJetMC);
					if(fHistITSClusterTrue[iCut] && track->GetITSNcls())fHistITSClusterTrue[iCut]->Fill(track->GetITSNcls(),fWeightJetJetMC);
					if(fHistITSchi2True[iCut] && track->GetITSchi2())fHistITSchi2True[iCut]->Fill(track->GetITSchi2(),fWeightJetJetMC);
				}	
			}	
			if(((((AliCaloSigmaCuts*)fSigmaCutArray->At(iCut))->TrackIsSelected(track, fPIDResponse)) == kTRUE)){
				if(fHistTrackDCAXY[iCut])fHistTrackDCAXY[iCut]->Fill(TMath::Abs(trackDCAXY),fWeightJetJetMC);
				if(fHistTrackDCAZ[iCut])fHistTrackDCAZ[iCut]->Fill(TMath::Abs(trackDCAZ),fWeightJetJetMC);
				if(fHistTPCSignal[iCut]) fHistTPCSignal[iCut]->Fill(track->P(), fPIDResponse->NumberOfSigmasTPC(track, AliPID::kProton),fWeightJetJetMC);
				if(fHistTPCClusterwCut[iCut] && track->GetTPCNcls())fHistTPCClusterwCut[iCut]->Fill(track->GetTPCNcls(),fWeightJetJetMC);
				if(fHistTPCchi2wCut[iCut] && track->GetTPCchi2perCluster())fHistTPCchi2wCut[iCut]->Fill(track->GetTPCchi2perCluster(),fWeightJetJetMC);
				if(fHistITSClusterwCut[iCut] && track->GetITSNcls())fHistITSClusterwCut[iCut]->Fill(track->GetITSNcls(),fWeightJetJetMC);
				if(fHistITSchi2wCut[iCut] && track->GetITSchi2())fHistITSchi2wCut[iCut]->Fill(track->GetITSchi2(),fWeightJetJetMC);
				if(fIsMC > 0){
					if((IsRealProton(track, fAODMCTrackArray, iCut, fWeightJetJetMC, 0)) > 0){
						if(fHistTrackDCAXYTrue[iCut])fHistTrackDCAXYTrue[iCut]->Fill(TMath::Abs(trackDCAXY),fWeightJetJetMC);
						if(fHistTrackDCAZTrue[iCut])fHistTrackDCAZTrue[iCut]->Fill(TMath::Abs(trackDCAZ),fWeightJetJetMC);
						if(fHistTPCClusterTruewCut[iCut] && track->GetTPCNcls())fHistTPCClusterTruewCut[iCut]->Fill(track->GetTPCNcls(),fWeightJetJetMC);
						if(fHistTPCchi2TruewCut[iCut] && track->GetTPCchi2perCluster())fHistTPCchi2TruewCut[iCut]->Fill(track->GetTPCchi2perCluster(),fWeightJetJetMC);
						if(fHistITSClusterTruewCut[iCut] && track->GetITSNcls())fHistITSClusterTruewCut[iCut]->Fill(track->GetITSNcls(),fWeightJetJetMC);
						if(fHistITSchi2TruewCut[iCut] && track->GetITSchi2())fHistITSchi2TruewCut[iCut]->Fill(track->GetITSchi2(),fWeightJetJetMC);
					}	
				}
				if(((((AliCaloSigmaCuts*)fSigmaCutArray->At(iCut))->TrackIsSelectedByDCACut(track)) == kTRUE)){
					if(fHistTrackDCAXYwCuts[iCut])fHistTrackDCAXYwCuts[iCut]->Fill(TMath::Abs(trackDCAXY),fWeightJetJetMC);
					if(fHistTrackDCAZwCuts[iCut])fHistTrackDCAZwCuts[iCut]->Fill(TMath::Abs(trackDCAZ),fWeightJetJetMC);
					if(fHistProtonPt[iCut]) fHistProtonPt[iCut]->Fill(track->Pt(),fWeightJetJetMC);
					if(fHistThetaPhiProton[iCut]) fHistThetaPhiProton[iCut]->Fill(track->Theta(), track->Phi(),fWeightJetJetMC);
					if(fIsMC > 0){
						if((IsRealProton(track, fAODMCTrackArray, iCut, fWeightJetJetMC, 0)) > 0){
							if(fHistTrackDCAXYTruewCuts[iCut])fHistTrackDCAXYTruewCuts[iCut]->Fill(TMath::Abs(trackDCAXY),fWeightJetJetMC);
							if(fHistTrackDCAZTruewCuts[iCut])fHistTrackDCAZTruewCuts[iCut]->Fill(TMath::Abs(trackDCAZ),fWeightJetJetMC);
						}
					}	
					proton.push_back(track);
				}
				else {
					
					tracks.push_back(track);
				}			
			}
			else {
				tracks.push_back(track);
			}
			if(fHistAllTracksPt[iCut]) fHistAllTracksPt[iCut]->Fill(track->Pt(),fWeightJetJetMC);
			if(fHistThetaPhi[iCut]) fHistThetaPhi[iCut]->Fill(track->Theta(), track->Phi(),fWeightJetJetMC);
			//dE/dx Plot
			Double_t ptot = track->P();
			Double_t tpcSignal = track->GetTPCsignal();
			if(fHistDEDx[iCut]) fHistDEDx[iCut]->Fill(ptot, tpcSignal,fWeightJetJetMC);
			if (track->GetTOFsignal()){
  				const float len = track->GetIntegratedLength();
				const float tim = track->GetTOFsignal() - fPIDResponse->GetTOFResponse().GetStartTime(track->GetTPCmomentum());
				const float beta = len / (tim * (2.99792457999999984e-02));
				if(fHistTOFBeta[iCut]) fHistTOFBeta[iCut]->Fill(ptot, beta,fWeightJetJetMC);
			}

		}
		fHistNProtonsPerEvent[iCut]->Fill(proton.size(),fWeightJetJetMC);
		if(proton.size() > 0){
			//Find gammas
			Int_t nclus                     = 0;
			Int_t nClusWCuts                = 0;
			nclus = fAOD->GetNumberOfCaloClusters();
			if(fHistNClusWoCuts[iCut]) fHistNClusWoCuts[iCut]->Fill(nclus,fWeightJetJetMC);
			if(nclus == 0){
				for(unsigned int i = 0; i < photon.size(); ++i) {
					if(photon[i]) delete photon[i];
				}
				photon.clear();
				for(unsigned int i = 0; i < proton.size(); ++i) {
					if(proton[i]) delete proton[i];
				}
				proton.clear();
				for(unsigned int i = 0; i < tracks.size(); ++i) {
					if(tracks[i]) delete tracks[i];
				}
				tracks.clear();
				// for(unsigned int i = 0; i < pions.size(); ++i) {
				// 	if(pions[i]) delete pions[i];
				// }
				pions.clear();
				continue;
			}  
			AliVCluster* clus = NULL;
			for(Int_t i=0; i < nclus; ++i){
				clus = new AliAODCaloCluster(*(AliAODCaloCluster*)fEvent->GetCaloCluster(i));
				if(!clus) {
					Printf("ERROR: Could not find clus %i",i);
					continue;
				}
				if(fHistClusterEWOCuts[iCut]) fHistClusterEWOCuts[iCut]-> Fill(clus->E(),fWeightJetJetMC);
				if((((AliCaloPhotonCuts*)fClusterCutArray->At(iCut))->ClusterIsSelected(clus,fAOD,fMCEvent,fIsMC, 1.,i)) == kTRUE){
					photon.push_back(clus);
					if(fHistClusterE[iCut]) fHistClusterE[iCut]-> Fill(clus->E(),fWeightJetJetMC);
					nClusWCuts+= 1;
				} else {
					delete clus;
				}
				clus = NULL;
			}
			if(fHistNClusWCuts[iCut]) fHistNClusWCuts[iCut]->Fill(nClusWCuts,fWeightJetJetMC);
			TLorentzVector protonVektor;
			TLorentzVector sigmaVektor;
			Int_t trueProtonMotherID = -1;
			Int_t truePhotonMotherID1 = -1;
			Int_t truePhotonMotherID2 = -1;

			Int_t trueProtonFromXi0ID = -1;
			Int_t truePhotonFromXi0ID1 = -1;
			Int_t truePhotonFromXi0ID2 = -1;

			if( photon.size() > 1){
				for(unsigned int iProton = 0; iProton < proton.size(); ++iProton){
					if(!proton[iProton]) {
						Printf("ERROR: Could not find proton[%i][%i]",iCut,iProton);
						continue;
					}
					vector<Int_t>         fVectorDoubleCountTrueSigmas;
					trueProtonFromXi0ID = -1;
					trueProtonMotherID = -1;
					AliAODTrack* protonCandidate = proton[iProton];
					protonVektor.SetPtEtaPhiM(protonCandidate->Pt(),protonCandidate->Eta(),protonCandidate->Phi(), 0.938272);
					if(fIsMC > 0){
						trueProtonMotherID = IsRealProton(protonCandidate, fAODMCTrackArray, iCut, fWeightJetJetMC, 1);
						trueProtonFromXi0ID = IsProtonFromXi0(protonCandidate, fAODMCTrackArray, iCut, fWeightJetJetMC, 1);
						if(trueProtonMotherID > 0){
							if(fHistTrueProtonPt[iCut]) fHistTrueProtonPt[iCut]->Fill(protonCandidate->Pt(),fWeightJetJetMC);
						}
					}
					for(unsigned int iPhoton1 = 0; iPhoton1 < photon.size(); ++iPhoton1) {
						if(!photon[iPhoton1]) {
							Printf("ERROR: Could not find photon[%i][%i]",iCut,iPhoton1);
							continue;
						}
						truePhotonFromXi0ID1 = -1;
						truePhotonMotherID1 = -1;
						AliVCluster* gamma1 = photon[iPhoton1];
						TLorentzVector clusterVector1;
						gamma1->GetMomentum(clusterVector1,vpos);
						AliAODConversionPhoton PhotonCandidate1 = AliAODConversionPhoton(&clusterVector1);
						// Flag Photon as CaloPhoton
						PhotonCandidate1.SetIsCaloPhoton(2);
						PhotonCandidate1.SetCaloClusterRef(iPhoton1);
						PhotonCandidate1.SetLeadingCellID(((AliCaloPhotonCuts*)fClusterCutArray->At(iCut))->FindLargestCellInCluster(gamma1,fEvent));
						// get MC label
						if(fIsMC> 0){
							Int_t* mclabelsCluster = gamma1->GetLabels();
							PhotonCandidate1.SetNCaloPhotonMCLabels(gamma1->GetNLabels());
							// cout << gamma1->GetNLabels() << endl;
							if (gamma1->GetNLabels()>0){
								for (Int_t k =0; k<(Int_t)gamma1->GetNLabels(); k++){
									PhotonCandidate1.SetCaloPhotonMCLabel(k,mclabelsCluster[k]);
									// Int_t pdgCode = fMCEvent->Particle(mclabelsCluster[k])->GetPdgCode();
									// cout << "label " << k << "\t" << mclabelsCluster[k] << " pdg code: " << pdgCode << endl;
								}
							}
							truePhotonMotherID1 = IsRealPhoton(&PhotonCandidate1, fAODMCTrackArray, iCut, fWeightJetJetMC);
							truePhotonFromXi0ID1 = IsPhotonFromXi0(&PhotonCandidate1, fAODMCTrackArray, iCut, fWeightJetJetMC);
						}
						for(unsigned int iPhoton2 = 0; iPhoton2 < photon.size(); ++iPhoton2) {
							if( iPhoton2 > iPhoton1){
								if(!photon[iPhoton2]) {
									Printf("ERROR: Could not find photon[%i][%i]",iCut,iPhoton2);
									continue;
								}
								truePhotonMotherID2 = -1;
								truePhotonFromXi0ID2 = -1;
								AliVCluster* gamma2 = photon[iPhoton2];
								TLorentzVector clusterVector2;
								gamma2->GetMomentum(clusterVector2,vpos);
								AliAODConversionPhoton PhotonCandidate2 = AliAODConversionPhoton(&clusterVector2);
								// Flag Photon as CaloPhoton
								PhotonCandidate2.SetIsCaloPhoton(2);
								PhotonCandidate2.SetCaloClusterRef(iPhoton2);
								PhotonCandidate2.SetLeadingCellID(((AliCaloPhotonCuts*)fClusterCutArray->At(iCut))->FindLargestCellInCluster(gamma2,fEvent));
								// get MC label
								if(fIsMC> 0){
									Int_t* mclabelsCluster = gamma2->GetLabels();
									PhotonCandidate2.SetNCaloPhotonMCLabels(gamma2->GetNLabels());
									// cout << gamma2->GetNLabels() << endl;
									if (gamma2->GetNLabels()>0){
										for (Int_t k =0; k<(Int_t)gamma2->GetNLabels(); k++){
											PhotonCandidate2.SetCaloPhotonMCLabel(k,mclabelsCluster[k]);
											// Int_t pdgCode = fMCEvent->Particle(mclabelsCluster[k])->GetPdgCode();
											// cout << "label " << k << "\t" << mclabelsCluster[k] << " pdg code: " << pdgCode << endl;
										}
									}
									truePhotonMotherID2 = IsRealPhoton(&PhotonCandidate2, fAODMCTrackArray, iCut, fWeightJetJetMC);
									truePhotonFromXi0ID2 = IsPhotonFromXi0(&PhotonCandidate2, fAODMCTrackArray, iCut, fWeightJetJetMC);
								}
								AliAODConversionMother pi0cand = AliAODConversionMother(&PhotonCandidate1,&PhotonCandidate2);
								pi0cand.SetLabels(iPhoton1,iPhoton2);
								if((((AliConversionMesonCuts*)fMesonCutArray->At(iCut))->MesonIsSelected(&pi0cand,kTRUE,((AliConvEventCuts*)fEventCutArray->At(iCut))->GetEtaShift(),PhotonCandidate1.GetLeadingCellID(),PhotonCandidate2.GetLeadingCellID()))){
									if( truePhotonMotherID1 == truePhotonMotherID2 && truePhotonMotherID1 > 0 && (fIsMC > 0)) {
										if((fHistReconstructedMassPi0MC[iCut]) && (iProton == 0)) fHistReconstructedMassPi0MC[iCut]->Fill(pi0cand.M(), pi0cand.Pt(),fWeightJetJetMC);
									}
									if((fHistReconstructedMassPi0[iCut]) && (iProton == 0)) fHistReconstructedMassPi0[iCut]->Fill(pi0cand.M(), pi0cand.Pt(),fWeightJetJetMC);
									if((((AliCaloSigmaCuts*)fSigmaCutArray->At(iCut))->PionIsSelectedByMassCut(pi0cand.M())) == kTRUE){
										TLorentzVector rekombinatedPi0;
										rekombinatedPi0.SetPtEtaPhiM(pi0cand.Pt(), pi0cand.Eta(), pi0cand.Phi(), 0.135);
										if((iProton == 0)) pions.push_back(rekombinatedPi0);
										sigmaVektor = protonVektor + rekombinatedPi0;
										if((fHistReconstructedMassPi0wCut[iCut]) && (iProton == 0)) fHistReconstructedMassPi0wCut[iCut]->Fill(pi0cand.M(), pi0cand.Pt(),fWeightJetJetMC);
										if( truePhotonMotherID1 == truePhotonMotherID2 && truePhotonMotherID1 > 0 && (fIsMC > 0)) {
											if((fHistReconstructedMassPi0MCwCut[iCut]) && (iProton == 0)) fHistReconstructedMassPi0MCwCut[iCut]->Fill(pi0cand.M(), pi0cand.Pt(),fWeightJetJetMC);
										}
										if(fHistPodolanski[iCut]) fHistPodolanski[iCut]->Fill(GetPodAlpha(sigmaVektor, protonVektor, rekombinatedPi0),GetQT(sigmaVektor, rekombinatedPi0),fWeightJetJetMC);
										if(trueProtonMotherID == truePhotonMotherID1 && trueProtonMotherID == truePhotonMotherID2 && trueProtonMotherID > 0 && (fIsMC > 0)) {
											if(fHistPodolanskiWCutTrue[iCut]) fHistPodolanskiWCutTrue[iCut]->Fill(GetPodAlpha(sigmaVektor, protonVektor, rekombinatedPi0),GetQT(sigmaVektor, rekombinatedPi0),fWeightJetJetMC);
											if(fHistSigmaMassPtWoPodCutMC[iCut]) fHistSigmaMassPtWoPodCutMC[iCut]->Fill(sigmaVektor.M(), sigmaVektor.Pt(),fWeightJetJetMC);

										}
										if(fHistSigmaMassPtWoPodCut[iCut])fHistSigmaMassPtWoPodCut[iCut]->Fill(sigmaVektor.M(), sigmaVektor.Pt(),fWeightJetJetMC);
										if((((AliCaloSigmaCuts*)fSigmaCutArray->At(iCut))->ArmenterosLikeQtCut(GetPodAlpha(sigmaVektor, protonVektor, rekombinatedPi0), GetQT(sigmaVektor, protonVektor))) == kFALSE) continue;
										if(fHistPodolanskiWCut[iCut]) fHistPodolanskiWCut[iCut]->Fill(GetPodAlpha(sigmaVektor, protonVektor, rekombinatedPi0),GetQT(sigmaVektor, rekombinatedPi0),fWeightJetJetMC);
										if(fHistDecayangle[iCut])fHistDecayangle[iCut]->Fill(rekombinatedPi0.Angle(protonVektor.Vect()), sigmaVektor.Pt(),fWeightJetJetMC);
										if(trueProtonMotherID == truePhotonMotherID1 && trueProtonMotherID == truePhotonMotherID2 && trueProtonMotherID > 0 && (fIsMC > 0)) {
											if(fHistPodolanskiWCutTrue[iCut]) fHistPodolanskiWCutTrue[iCut]->Fill(GetPodAlpha(sigmaVektor, protonVektor, rekombinatedPi0),GetQT(sigmaVektor, rekombinatedPi0),fWeightJetJetMC);
											if(fHistDecayangleTrue[iCut])fHistDecayangleTrue[iCut]->Fill(rekombinatedPi0.Angle(protonVektor.Vect()), sigmaVektor.Pt(),fWeightJetJetMC);
										}
										if((((AliCaloSigmaCuts*)fSigmaCutArray->At(iCut))->SigmaDaughtersOpeningangleCut(rekombinatedPi0.Angle(protonVektor.Vect()))) == kFALSE) continue;
										if(fHistDecayanglewCut[iCut])fHistDecayanglewCut[iCut]->Fill(rekombinatedPi0.Angle(protonVektor.Vect()), sigmaVektor.Pt(),fWeightJetJetMC);
										if(fHistSigmaPlus[iCut]) fHistSigmaPlus[iCut]-> Fill(sigmaVektor.M(), sigmaVektor.Pt(),fWeightJetJetMC);
										if(trueProtonMotherID == truePhotonMotherID1 && trueProtonMotherID == truePhotonMotherID2 && trueProtonMotherID > 0 && (fIsMC > 0)) {
											if(fHistDecayangleTruewCut[iCut])fHistDecayangleTruewCut[iCut]->Fill(rekombinatedPi0.Angle(protonVektor.Vect()), sigmaVektor.Pt(),fWeightJetJetMC);
											if(fHistSigmaPlusMC[iCut]) fHistSigmaPlusMC [iCut]-> Fill(sigmaVektor.M(), sigmaVektor.Pt(),fWeightJetJetMC);
											if (CheckVectorForDoubleCount(fVectorDoubleCountTrueSigmas,trueProtonMotherID)) fHistDoubleCountTrueSigmaInvMassPt[iCut]->Fill(sigmaVektor.M(), sigmaVektor.Pt(),fWeightJetJetMC);
										}
										if(((trueProtonMotherID == truePhotonMotherID1) && (truePhotonMotherID1 != truePhotonMotherID2) && (fIsMC > 0) && trueProtonMotherID > 0) || ((trueProtonMotherID == truePhotonMotherID2) && (truePhotonMotherID1 != truePhotonMotherID2) && (fIsMC > 0) && trueProtonMotherID > 0)){
											if(fHistSigmaPlusMCTrueProtonGamma[iCut]) fHistSigmaPlusMCTrueProtonGamma [iCut]-> Fill(sigmaVektor.M(), sigmaVektor.Pt(),fWeightJetJetMC);
										}
										if(((fIsMC > 0) && trueProtonMotherID > 0) && (trueProtonMotherID != truePhotonMotherID2)  && (trueProtonMotherID != truePhotonMotherID1)){
											if(fHistSigmaPlusMCTrueProton[iCut]) fHistSigmaPlusMCTrueProton [iCut]-> Fill(sigmaVektor.M(), sigmaVektor.Pt(),fWeightJetJetMC);
										}
										if( (truePhotonMotherID1 == truePhotonMotherID2) && truePhotonMotherID1 > 0 && fIsMC > 0 && (trueProtonMotherID != truePhotonMotherID2)){
											if(fHistSigmaPlusMCTruePion[iCut]) fHistSigmaPlusMCTruePion [iCut]-> Fill(sigmaVektor.M(), sigmaVektor.Pt(),fWeightJetJetMC);
										}
										if(trueProtonFromXi0ID == truePhotonFromXi0ID1 && trueProtonFromXi0ID == truePhotonFromXi0ID2 && trueProtonFromXi0ID > 0 && (fIsMC > 0)) {
											if(fHistXi0MC[iCut]) fHistXi0MC [iCut]-> Fill(sigmaVektor.M(), sigmaVektor.Pt(),fWeightJetJetMC);
										}
									}
								}
							}
						}
					}
					fVectorDoubleCountTrueSigmas.clear();				
					fVectorDoubleCountTrueSigmas.resize(0);				
				}
			}
		}
		if((((AliCaloSigmaCuts*)fSigmaCutArray->At(iCut))->UseRotationmethod()) == 0){
			CalculateBackgroundSwappWGammaGamma(photon, proton, vpos, iCut,fWeightJetJetMC);
		}
		if((((AliCaloSigmaCuts*)fSigmaCutArray->At(iCut))->UseRotationmethod()) == 1){
			CalculateBackgroundSwappWProtonPion(pions, proton, iCut,fWeightJetJetMC);
		}

		for(unsigned int i = 0; i < photon.size(); ++i) {
			if(photon[i]) delete photon[i];
		}
		photon.clear();
		for(unsigned int i = 0; i < proton.size(); ++i) {
			if(proton[i]) delete proton[i];
		}
		proton.clear();
		for(unsigned int i = 0; i < tracks.size(); ++i) {
			if(tracks[i]) delete tracks[i];
		}
		tracks.clear();
		// for(unsigned int i = 0; i < pions.size(); ++i) {
		// 	if(pions[i]) delete pions[i];
		// }
		pions.clear();
	}
	// continue until all the tracks are processed
	PostData(1, fOutputList);                           // stream the results the analysis of this event to
	// the output manager which will take care of writing
	// it to a file
}
//_________________________________________________________________________________
Bool_t AliAnalysisTaskSigmaPlToProtonPiZeroAOD::CheckVectorForDoubleCount(vector<Int_t> &vec, Int_t tobechecked)
{
  if(tobechecked > -1)
  {
    vector<Int_t>::iterator it;
    it = find (vec.begin(), vec.end(), tobechecked);
    if (it != vec.end()) return true;
    else{
      vec.push_back(tobechecked);
      return false;
    }
  }
  return false;
}
//_____________________________________________________________________________
Double_t AliAnalysisTaskSigmaPlToProtonPiZeroAOD::GetPodAlpha(TLorentzVector sigmaVektor, TLorentzVector protonVektor, TLorentzVector rekombinatedPi0)
{
	return (cos(sigmaVektor.Angle(protonVektor.Vect())) * protonVektor.P() - cos(sigmaVektor.Angle(rekombinatedPi0.Vect())) * rekombinatedPi0.P()) / (cos(sigmaVektor.Angle(protonVektor.Vect())) * protonVektor.P() + cos(sigmaVektor.Angle(rekombinatedPi0.Vect())) * rekombinatedPi0.P());

}
//_____________________________________________________________________________
Double_t AliAnalysisTaskSigmaPlToProtonPiZeroAOD::GetQT(TLorentzVector sigmaVektor, TLorentzVector rekombinatedPi0)
{
	return  sin(sigmaVektor.Angle(rekombinatedPi0.Vect())) * rekombinatedPi0.P();
}
//________________________________________________________________________
Int_t AliAnalysisTaskSigmaPlToProtonPiZeroAOD::IsRealProton(AliAODTrack* track, TClonesArray* fAODMCTrackArray, Int_t iCut, Double_t fWeightJetJetMC, Int_t fill)
{//Check if proton comes frome a sigma+
	Int_t mcID = track->GetLabel();
	AliAODMCParticle* mcTrack = static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(TMath::Abs(mcID)));
	if(!mcTrack) return -1;
	if (mcTrack->GetMother() < 0) return -1;
	AliAODMCParticle *TrackMother = static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(mcTrack->GetMother())); //mother MC particle object
	if(!TrackMother) return -1;
	Int_t codeTrack = mcTrack->GetPdgCode();
	Int_t codeMother = TrackMother->GetPdgCode();
	if( (codeTrack==2212) && (codeMother==3222) ){
		return (TrackMother->GetLabel());
	} else if((codeTrack==2212) && (codeMother==2212) ){
		Int_t codeGrandMother = codeMother;
		Int_t labelGrandMother = mcTrack->GetMother();
		Int_t counter = 0;
		while(codeGrandMother == 2212 && counter < 10 && labelGrandMother < fAODMCTrackArray->GetEntriesFast()){
			AliAODMCParticle *TrackGrandMother = static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(labelGrandMother)); //mother MC particle object
			if(!TrackGrandMother) return -1;
			codeGrandMother = TrackGrandMother->GetPdgCode();
			if(codeGrandMother == 3222){
				if(fill == 1) fHistNLoopsProton[iCut]->Fill(counter, fWeightJetJetMC);
				return TrackGrandMother->GetLabel();
			}	
			labelGrandMother = TrackGrandMother->GetMother();
			counter += 1;
		}
		return -1;
	} else {
		return -1;
	}
	return -1;
}
//________________________________________________________________________
Int_t AliAnalysisTaskSigmaPlToProtonPiZeroAOD::IsRealPhoton(AliAODConversionPhoton *PhotonCandidate, TClonesArray* fAODMCTrackArray, Int_t iCut, Double_t fWeightJetJetMC)
{ //checks if a reconstructed photon from sigma plus
	AliAODMCParticle *Photon = NULL;
	if (PhotonCandidate->GetNCaloPhotonMCLabels() > 0){
		// Photon = PhotonCandidate->GetMCParticle(fMCEvent);
		if (PhotonCandidate->GetCaloPhotonMCLabel(0) < 0) return -1;
		Photon = static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(PhotonCandidate->GetCaloPhotonMCLabel(0)));
		if (Photon) {
			if(Photon->GetPdgCode() == 22){
				if (Photon->GetMother() < 0) return -1;
				AliAODMCParticle* motherPart2 = static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(Photon->GetMother()));
				if (motherPart2->GetMother() < 0) return -1;
				AliAODMCParticle* grandmotherPart2 = static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(motherPart2->GetMother()));
				if(motherPart2->GetPdgCode() == 111 && grandmotherPart2->GetPdgCode() == 3222){
					return (grandmotherPart2->GetLabel());	
				} else if(motherPart2->GetPdgCode() == 11 || motherPart2->GetPdgCode() == -11 || motherPart2->GetPdgCode() == 22){
					Int_t codeGrandMother = motherPart2->GetPdgCode();
					Int_t codePotentialMother = Photon->GetPdgCode();
					Int_t labelGrandMother = Photon->GetMother();
					Int_t counter = 0;
					while((codeGrandMother == 11 || codeGrandMother == -11 || codeGrandMother == 22 || codeGrandMother == 111) && counter < 10 && labelGrandMother < fAODMCTrackArray->GetEntriesFast()){
						AliAODMCParticle *TrackGrandMother = static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(labelGrandMother)); //mother MC particle object
						if(!TrackGrandMother) return -1;
						codeGrandMother = TrackGrandMother->GetPdgCode();
						if(codePotentialMother == 111 && codeGrandMother == 3222){
							fHistNLoopsGamma[iCut]->Fill(counter, fWeightJetJetMC);
							return TrackGrandMother->GetLabel();
						}
						labelGrandMother = TrackGrandMother->GetMother();
						codePotentialMother = codeGrandMother;
						counter += 1;
					}
					return -1;
				}else {
					return -1;
				}
			} else if(Photon->GetPdgCode() == 11 || Photon->GetPdgCode() == -11){
				if (Photon->GetMother() < 0) return -1;
				AliAODMCParticle* motherPart2 = static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(Photon->GetMother()));
				if (motherPart2->GetMother() < 0) return -1;
				AliAODMCParticle* grandmotherPart2 = static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(motherPart2->GetMother()));
				if (grandmotherPart2->GetMother() < 0) return -1;
				AliAODMCParticle* grandgrandmotherPart2 = static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(grandmotherPart2->GetMother()));
				if(motherPart2->GetPdgCode() == 22 && grandmotherPart2->GetPdgCode() == 111  && grandgrandmotherPart2->GetPdgCode() == 3222){
					return (grandgrandmotherPart2->GetLabel());
				} else if((motherPart2->GetPdgCode() == 11 || motherPart2->GetPdgCode() == -11 || motherPart2->GetPdgCode() == 22) && (grandmotherPart2->GetPdgCode() == 11 || grandmotherPart2->GetPdgCode() == -11 || grandmotherPart2->GetPdgCode() == 22)){
					Int_t codeGrandMother = motherPart2->GetPdgCode();
					Int_t codePotentialMother = Photon->GetPdgCode();
					Int_t labelGrandMother = Photon->GetMother();
					Int_t counter = 0;
					while((codeGrandMother == 11 || codeGrandMother == -11 || codeGrandMother == 22 || codeGrandMother == 111) && counter < 10 && labelGrandMother < fAODMCTrackArray->GetEntries()){
						AliAODMCParticle *TrackGrandMother = static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(labelGrandMother)); //mother MC particle object
						if(!TrackGrandMother) return -1;
						codeGrandMother = TrackGrandMother->GetPdgCode();
						if(codePotentialMother == 111 && codeGrandMother == 3222){
							fHistNLoopsGamma[iCut]->Fill(counter, fWeightJetJetMC);
							return TrackGrandMother->GetLabel();
						}
						labelGrandMother = TrackGrandMother->GetMother();
						codePotentialMother = codeGrandMother;
						counter += 1;
					}
					return -1;
				}else {
					return -1;
				}
			} else {
				return -1;
			}
		}
	}
	return -1;

}

//________________________________________________________________________
void AliAnalysisTaskSigmaPlToProtonPiZeroAOD::CalculateBackgroundSwappWGammaGamma(vector < AliVCluster* > photon, vector < AliAODTrack* > proton, Double_t vpos[3], Int_t iCut, Double_t fWeightJetJetMC){

    // Double_t rotationAngle = TMath::Pi()/8.0; //0.78539816339; // rotaion angle 45°
    // Double_t rotationAngle = TMath::Pi()/2.0; //0.78539816339; // rotaion angle 90°

    TLorentzVector lvRotationPhoton1;   // photon candidates which get rotated
    TLorentzVector lvRotationPhoton;   // photon candidates which get rotated
    TLorentzVector lvRotationBGPhoton;   // photon candidates which get rotated
    TLorentzVector lvRotationBGProton;   // photon candidates which get rotated
    TVector3 vRotationPion;            // reconstructed mother particle from the two photons
    TLorentzVector lvRotationBGPion;            // reconstructed mother particle from the two photons
    TLorentzVector lvRotationBGPion1;            // reconstructed mother particle from the two photons
    TLorentzVector lvRotationBGSigma;            // reconstructed mother particle from the two photons
    TLorentzVector lvRotationBGSigma1;            // reconstructed mother particle from the two photons
    Int_t cellIDRotatedPhoton = -1; // cell ID of the cluster after rotation
    Int_t cellIDRotatedPhoton1 = -1; // cell ID of the cluster after rotation
    Bool_t bPhotonAccepted = kFALSE;
    Bool_t bPhoton1Accepted = kFALSE;

    TVector3 tvEtaPhigamma1, tvEtaPhigamma2, tvEtaPhigamma1Decay, tvEtaPhigamma2Decay, tvNormBeforeDecay, tvNormAfterDecay;
    Float_t asymBeforeDecay = 0.;
    Float_t asymAfterDecay = 0.;
    Double_t massGamma[2] = {0,0};

    std::vector<std::array<Double_t, 2>> vSwappingInvMassPT;
    vSwappingInvMassPT.clear();
    vSwappingInvMassPT.resize(0);
    Double_t tempMultWeightSwapping = 1.; // weight taking multiplicity of event into account

    // curcial requierment is that the event has at least 3 cluster candidates
    if(photon.size() > 1 && proton.size() > 0 ){
    	for(unsigned int iProtonBG=0;iProtonBG<proton.size();iProtonBG++){
    		AliAODTrack* protonBGCandidate = proton[iProtonBG];
    		lvRotationBGProton.SetX(protonBGCandidate->Px());
            lvRotationBGProton.SetY(protonBGCandidate->Py());
            lvRotationBGProton.SetZ(protonBGCandidate->Pz());
            lvRotationBGProton.SetE(protonBGCandidate->E());

	    	for(unsigned int iCurrent1=0;iCurrent1<photon.size();iCurrent1++){
		        AliVCluster* photonCandidate = photon[iCurrent1];

	        	for(unsigned int iCurrent2=iCurrent1+1;iCurrent2<photon.size();iCurrent2++){
			          AliVCluster* photon1Candidate = photon[iCurrent2];

			        for(int iSwapp = 0; iSwapp < 20; ++iSwapp){
						photonCandidate->GetMomentum(lvRotationPhoton,vpos);
						photon1Candidate->GetMomentum(lvRotationPhoton1,vpos);

			            vRotationPion = (lvRotationPhoton + lvRotationPhoton1).Vect();

			            // //Rotation um festen Winkel
			            // lvRotationPhoton.Rotate(rotationAngle, vRotationPion);
			            // lvRotationPhoton1.Rotate(rotationAngle, vRotationPion);

			            //Rotation TGenPhasespace
		                tvEtaPhigamma1 = lvRotationPhoton.Vect();
		                tvEtaPhigamma2 = lvRotationPhoton1.Vect();
		                tvNormBeforeDecay = tvEtaPhigamma1.Cross(tvEtaPhigamma2);
		                asymBeforeDecay = fabs((lvRotationPhoton.E()-lvRotationPhoton1.E())/(lvRotationPhoton.E()+lvRotationPhoton1.E()));

			            TLorentzVector lvRotationMother = lvRotationPhoton + lvRotationPhoton1;
			            fGenPhaseSpace.SetDecay(lvRotationMother, 2, massGamma);
			            fGenPhaseSpace.Generate();
			            lvRotationPhoton = *fGenPhaseSpace.GetDecay(0);
			            lvRotationPhoton1 = *fGenPhaseSpace.GetDecay(1);

		                tvEtaPhigamma1Decay = lvRotationPhoton.Vect();
		                tvEtaPhigamma2Decay = lvRotationPhoton1.Vect();
		                tvNormAfterDecay = tvEtaPhigamma1Decay.Cross(tvEtaPhigamma2Decay);  // norm vector to decay plane
		                asymAfterDecay = fabs((lvRotationPhoton.E()-lvRotationPhoton1.E())/(lvRotationPhoton.E()+lvRotationPhoton1.E()));
		                // check if decay is nearly the same as original decay: if yes continue with next decay
		                if((tvNormAfterDecay.Angle(tvNormBeforeDecay) < 20*TMath::Pi()/180. || tvNormAfterDecay.Angle(tvNormBeforeDecay) > 340*TMath::Pi()/180.) && ( fabs(asymBeforeDecay - asymAfterDecay) < 0.05 )   ) continue;


						if(lvRotationPhoton1.Phi() < 0.){
		            		if((fabs(lvRotationPhoton1.Eta()) < 0.13) && ((lvRotationPhoton1.Phi()+2.*TMath::Pi()) > 250.f*(TMath::Pi()/180.f)) && ((lvRotationPhoton1.Phi()+2.*TMath::Pi()) < 320.f*(TMath::Pi()/180.f))){
		            		bPhoton1Accepted = kTRUE;
		            		}
		            	}
		            	else{
			            	if( (fabs(lvRotationPhoton1.Eta()) < 0.13) && (lvRotationPhoton1.Phi() > 250.f*(TMath::Pi()/180.f)) && (lvRotationPhoton1.Phi() < 320.f*(TMath::Pi()/180.f))){
			            		bPhoton1Accepted = kTRUE;
			            	}
		            	}

		            	if(lvRotationPhoton.Phi() < 0.){
		            		if( (fabs(lvRotationPhoton.Eta()) < 0.13) && ((lvRotationPhoton.Phi()+2.*TMath::Pi()) > 250.f*(TMath::Pi()/180.f)) && ((lvRotationPhoton.Phi()+2.*TMath::Pi()) < 320.f*(TMath::Pi()/180.f))){
		            			bPhotonAccepted = kTRUE;
		            		}
		            	}
		            	else{
			            	if( (fabs(lvRotationPhoton.Eta()) < 0.13) && (lvRotationPhoton.Phi() > 250.f*(TMath::Pi()/180.f)) && (lvRotationPhoton.Phi() < 320.f*(TMath::Pi()/180.f))){
			            		bPhotonAccepted = kTRUE;
			            	}
		            	}

		            	if(bPhotonAccepted == kFALSE && bPhoton1Accepted == kFALSE) continue;

		            	cellIDRotatedPhoton = ((AliCaloPhotonCuts*)fClusterCutArray->At(iCut))->GetCaloCellIdFromEtaPhi(lvRotationPhoton.Eta(), static_cast<double>((lvRotationPhoton.Phi()<0) ? lvRotationPhoton.Phi() + TMath::Pi()*2. : lvRotationPhoton.Phi()));
	           			cellIDRotatedPhoton1 = ((AliCaloPhotonCuts*)fClusterCutArray->At(iCut))->GetCaloCellIdFromEtaPhi(lvRotationPhoton1.Eta(), static_cast<double>((lvRotationPhoton1.Phi()<0) ? lvRotationPhoton1.Phi() + TMath::Pi()*2. : lvRotationPhoton1.Phi()));
	           			std::unique_ptr<AliAODConversionPhoton> currentEventRotation (new AliAODConversionPhoton(&lvRotationPhoton));
	           			std::unique_ptr<AliAODConversionPhoton> currentEventRotation1 (new AliAODConversionPhoton(&lvRotationPhoton1));

		            	for(unsigned int iPionBG=0;iPionBG<photon.size();iPionBG++){
		            		if(iPionBG == iCurrent2 || iPionBG == iCurrent1) continue;
			          		AliVCluster* photonBGCandidate = photon[iPionBG];
							photonBGCandidate->GetMomentum(lvRotationBGPhoton,vpos);
		           			std::unique_ptr<AliAODConversionPhoton> currentEventGoodBGPhoton (new AliAODConversionPhoton(&lvRotationBGPhoton));
				          	if(bPhotonAccepted == kTRUE){
				          		//First Swapped Gamma
			            		lvRotationBGPion = (lvRotationPhoton + lvRotationBGPhoton);
			            		if(lvRotationBGPion.M() < 0.118 || lvRotationBGPion.M() > 0.148 ) continue;
				          		TLorentzVector rekombinatedPi0BG;
								rekombinatedPi0BG.SetPtEtaPhiM(lvRotationBGPion.Pt(), lvRotationBGPion.Eta(), lvRotationBGPion.Phi(), 0.135);
								std::unique_ptr<AliAODConversionMother> backgroundCandidate(new AliAODConversionMother(currentEventRotation.get(), currentEventGoodBGPhoton.get()));
								if(!(((AliCaloPhotonCuts*)fClusterCutArray->At(iCut))->CheckDistanceToBadChannelSwapping(cellIDRotatedPhoton, lvRotationPhoton.Phi(), fInputEvent)) && lvRotationPhoton.E() > ((AliCaloPhotonCuts*)fClusterCutArray->At(iCut))->GetMinClusterEnergy())
		              			{
		               				if(((AliConversionMesonCuts*) fMesonCutArray->At(iCut))->MesonIsSelected(backgroundCandidate.get(),kFALSE,((AliConvEventCuts*)fEventCutArray->At(iCut))->GetEtaShift(), cellIDRotatedPhoton, currentEventGoodBGPhoton.get()->GetLeadingCellID()))
		               				{
					          		//First Pion
							            lvRotationBGSigma = (rekombinatedPi0BG + lvRotationBGProton);
										if((((AliCaloSigmaCuts*)fSigmaCutArray->At(iCut))->ArmenterosLikeQtCut(GetPodAlpha(lvRotationBGSigma, lvRotationBGProton, rekombinatedPi0BG), GetQT(lvRotationBGSigma, lvRotationBGProton))) == kFALSE) continue;
										if((((AliCaloSigmaCuts*)fSigmaCutArray->At(iCut))->SigmaDaughtersOpeningangleCut(rekombinatedPi0BG.Angle(lvRotationBGProton.Vect()))) == kFALSE) continue;
					                 	vSwappingInvMassPT.push_back({lvRotationBGSigma.M(),lvRotationBGSigma.Pt()});
			                   		}
		              			}
		              		}

	              			if(bPhoton1Accepted == kTRUE){
		              			//Second Swapped Gamma
			                  	lvRotationBGPion1 = (lvRotationPhoton1 + lvRotationBGPhoton);
			            		if(lvRotationBGPion1.M() < 0.118 || lvRotationBGPion1.M() > 0.148 ) continue;
								TLorentzVector rekombinatedPi0BG1;
								rekombinatedPi0BG1.SetPtEtaPhiM(lvRotationBGPion1.Pt(), lvRotationBGPion1.Eta(), lvRotationBGPion1.Phi(), 0.135);
								std::unique_ptr<AliAODConversionMother> backgroundCandidate1(new AliAODConversionMother(currentEventRotation1.get(), currentEventGoodBGPhoton.get()));
		              			if(!(((AliCaloPhotonCuts*)fClusterCutArray->At(iCut))->CheckDistanceToBadChannelSwapping(cellIDRotatedPhoton1, lvRotationPhoton1.Phi(), fInputEvent)) && lvRotationPhoton1.E() > ((AliCaloPhotonCuts*)fClusterCutArray->At(iCut))->GetMinClusterEnergy())
		              			{
		               				if(((AliConversionMesonCuts*) fMesonCutArray->At(iCut))->MesonIsSelected(backgroundCandidate1.get(),kFALSE,((AliConvEventCuts*)fEventCutArray->At(iCut))->GetEtaShift(), cellIDRotatedPhoton1, currentEventGoodBGPhoton.get()->GetLeadingCellID()))
		               				{
					          		//Second Pion
					            		lvRotationBGSigma1 = (rekombinatedPi0BG1 + lvRotationBGProton);
										if((((AliCaloSigmaCuts*)fSigmaCutArray->At(iCut))->ArmenterosLikeQtCut(GetPodAlpha(lvRotationBGSigma1, lvRotationBGProton, rekombinatedPi0BG1), GetQT(lvRotationBGSigma1, lvRotationBGProton))) == kFALSE) continue;
										if(!(((AliCaloSigmaCuts*)fSigmaCutArray->At(iCut))->SigmaDaughtersOpeningangleCut(rekombinatedPi0BG1.Angle(lvRotationBGProton.Vect())))) continue;
					                 	vSwappingInvMassPT.push_back({lvRotationBGSigma1.M(),lvRotationBGSigma1.Pt()});

			                   		}
		              			}
		              		}
	          			}
	        		}
				}
			}
		}
      // Fill the histograms
	    // if(vSwappingInvMassPT.size() > 0){
    	//     tempMultWeightSwapping = (0.5*(fClusterCandidates->GetEntries()*fClusterCandidates->GetEntries() - fClusterCandidates->GetEntries()))/(vSwappingInvMassPT.size());
     //  	}
        for(Int_t i = 0; i < (Int_t)vSwappingInvMassPT.size(); i++){
        	fHistRotationWGammaGamma[iCut]->Fill(vSwappingInvMassPT.at(i)[0], vSwappingInvMassPT.at(i)[1], tempMultWeightSwapping*fWeightJetJetMC);
        }
    }
}


//________________________________________________________________________
void AliAnalysisTaskSigmaPlToProtonPiZeroAOD::CalculateBackgroundSwappWProtonPion(vector < TLorentzVector > pions, vector < AliAODTrack* > proton, Int_t iCut, Double_t fWeightJetJetMC){

    TLorentzVector lvRotationProton;            // reconstructed mother particle from the two photons
    TLorentzVector lvRotationBGProton;            // reconstructed mother particle from the two photons
    TLorentzVector lvRotationBGSigma;            // reconstructed mother particle from the two photons
    TLorentzVector lvRotationBGSigma1;            // reconstructed mother particle from the two photons
    Bool_t bPionAccepted = kFALSE;

    TVector3 tvEtaPhiPion, tvEtaPhiProton, tvEtaPhiPionDecay, tvEtaPhiProtonDecay, tvNormBeforeDecay, tvNormAfterDecay;
    Float_t asymBeforeDecay = 0.;
    Float_t asymAfterDecay = 0.;
    Double_t massDaughters[2] = {0.135,0.938};

    std::vector<std::array<Double_t, 2>> vSwappingInvMassPT;
    vSwappingInvMassPT.clear();
    vSwappingInvMassPT.resize(0);
    Double_t tempMultWeightSwapping = 1.; // weight taking multiplicity of event into account

    // curcial requierment is that the event has at least 3 cluster candidates
    if(pions.size() > 0 && proton.size() > 0 ){
    	for(unsigned int iCurrentProton=0;iCurrentProton<proton.size();iCurrentProton++){
    		AliAODTrack* protonCandidate = proton[iCurrentProton];
    		lvRotationProton.SetX(protonCandidate->Px());
            lvRotationProton.SetY(protonCandidate->Py());
            lvRotationProton.SetZ(protonCandidate->Pz());
            lvRotationProton.SetE(protonCandidate->E());
	    	for(unsigned int iCurrent=0;iCurrent<pions.size();iCurrent++){
		        TLorentzVector lvRotationPion = pions[iCurrent];
		        for(int iSwapp = 0; iSwapp < 4; ++iSwapp){
		            //Rotation TGenPhasespace
	                tvEtaPhiPion = lvRotationPion.Vect();
	                tvEtaPhiProton = lvRotationProton.Vect();
	                tvNormBeforeDecay = tvEtaPhiPion.Cross(tvEtaPhiProton);
	                asymBeforeDecay = fabs((lvRotationPion.E()-lvRotationProton.E())/(lvRotationPion.E()+lvRotationProton.E()));

		            TLorentzVector lvRotationMother = lvRotationPion + lvRotationProton;
		            fGenPhaseSpace.SetDecay(lvRotationMother, 2, massDaughters);
		            fGenPhaseSpace.Generate();
		            lvRotationPion = *fGenPhaseSpace.GetDecay(0);
		            lvRotationProton = *fGenPhaseSpace.GetDecay(1);

	                tvEtaPhiPionDecay = lvRotationPion.Vect();
	                tvEtaPhiProtonDecay = lvRotationProton.Vect();
	                tvNormAfterDecay = tvEtaPhiPionDecay.Cross(tvEtaPhiProtonDecay);  // norm vector to decay plane
	                asymAfterDecay = fabs((lvRotationPion.E()-lvRotationProton.E())/(lvRotationPion.E()+lvRotationProton.E()));
	                // check if decay is nearly the same as original decay: if yes continue with next decay
	                if((tvNormAfterDecay.Angle(tvNormBeforeDecay) < 20*TMath::Pi()/180. || tvNormAfterDecay.Angle(tvNormBeforeDecay) > 340*TMath::Pi()/180.) && ( fabs(asymBeforeDecay - asymAfterDecay) < 0.05 )   ) continue;

	            	if(lvRotationPion.Phi() < 0.){
	            		if( (fabs(lvRotationPion.Eta()) < 0.13) && ((lvRotationPion.Phi()+2.*TMath::Pi()) > 250.f*(TMath::Pi()/180.f)) && ((lvRotationPion.Phi()+2.*TMath::Pi()) < 320.f*(TMath::Pi()/180.f))){
	            			bPionAccepted = kTRUE;
	            		}
	            	}
	            	else{
		            	if( (fabs(lvRotationPion.Eta()) < 0.13) && (lvRotationPion.Phi() > 250.f*(TMath::Pi()/180.f)) && (lvRotationPion.Phi() < 320.f*(TMath::Pi()/180.f))){
		            		bPionAccepted = kTRUE;
		            	}
	            	}
	            	for(unsigned int iPionBG=0;iPionBG<pions.size();iPionBG++){
	            		if(iPionBG == iCurrent) continue;
		          		TLorentzVector lvRotationBGPion = pions[iPionBG];
			            lvRotationBGSigma = (lvRotationBGPion + lvRotationProton);
						if((((AliCaloSigmaCuts*)fSigmaCutArray->At(iCut))->ArmenterosLikeQtCut(GetPodAlpha(lvRotationBGSigma, lvRotationProton, lvRotationBGPion), GetQT(lvRotationBGSigma, lvRotationProton))) == kFALSE) continue;
						if(!(((AliCaloSigmaCuts*)fSigmaCutArray->At(iCut))->SigmaDaughtersOpeningangleCut(lvRotationBGPion.Angle(lvRotationProton.Vect())))) continue;
	                 	vSwappingInvMassPT.push_back({lvRotationBGSigma.M(),lvRotationBGSigma.Pt()});

	              	}
              		if(bPionAccepted == kTRUE){
						for(unsigned int iProtonBG=0;iProtonBG<proton.size();iProtonBG++){
	            			if(iProtonBG == iCurrentProton) continue;
    						AliAODTrack* protonBGCandidate = proton[iCurrentProton];
				    		lvRotationBGProton.SetX(protonBGCandidate->Px());
				            lvRotationBGProton.SetY(protonBGCandidate->Py());
				            lvRotationBGProton.SetZ(protonBGCandidate->Pz());
				            lvRotationBGProton.SetE(protonBGCandidate->E());

				            lvRotationBGSigma1 = (lvRotationPion + lvRotationBGProton);
							if((((AliCaloSigmaCuts*)fSigmaCutArray->At(iCut))->ArmenterosLikeQtCut(GetPodAlpha(lvRotationBGSigma1, lvRotationBGProton, lvRotationPion), GetQT(lvRotationBGSigma1, lvRotationBGProton))) == kFALSE) continue;
							if(!(((AliCaloSigmaCuts*)fSigmaCutArray->At(iCut))->SigmaDaughtersOpeningangleCut(lvRotationPion.Angle(lvRotationBGProton.Vect())))) continue;
		                 	vSwappingInvMassPT.push_back({lvRotationBGSigma1.M(),lvRotationBGSigma1.Pt()});
	              		}
          			}
				}
			}
		}
      // Fill the histograms
	    // if(vSwappingInvMassPT.size() > 0){
    	//     tempMultWeightSwapping = (0.5*(fClusterCandidates->GetEntries()*fClusterCandidates->GetEntries() - fClusterCandidates->GetEntries()))/(vSwappingInvMassPT.size());
     //  	}
        for(Int_t i = 0; i < (Int_t)vSwappingInvMassPT.size(); i++){
        	fHistRotationWProtonPion[iCut]->Fill(vSwappingInvMassPT.at(i)[0], vSwappingInvMassPT.at(i)[1], tempMultWeightSwapping*fWeightJetJetMC);
        }
    }
}
//________________________________________________________________________
Int_t AliAnalysisTaskSigmaPlToProtonPiZeroAOD::IsProtonFromXi0(AliAODTrack* track, TClonesArray* fAODMCTrackArray, Int_t iCut, Double_t fWeightJetJetMC, Int_t fill)
{//Check if proton comes frome a sigma+
	Int_t mcIDProtonCandidate = track->GetLabel();
	AliAODMCParticle* mcTrackProtonFromXi0Candidate = static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(TMath::Abs(mcIDProtonCandidate)));
	if(!mcTrackProtonFromXi0Candidate) return -1;
	if (mcTrackProtonFromXi0Candidate->GetMother() < 0) return -1;
	AliAODMCParticle *TrackLamdaCandidate = static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(mcTrackProtonFromXi0Candidate->GetMother())); //mother MC particle object
	if(!TrackLamdaCandidate) return -1;
	AliAODMCParticle *TrackXi0Candidate = static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(TrackLamdaCandidate->GetMother())); //mother MC particle object
	if(!TrackXi0Candidate) return -1;
	Int_t codeProtonCandidate = mcTrackProtonFromXi0Candidate->GetPdgCode();
	Int_t codeLamdaCandidate = TrackLamdaCandidate->GetPdgCode();
	Int_t codeXi0Candidate = TrackXi0Candidate->GetPdgCode();
	if( (codeProtonCandidate==2212) && (codeLamdaCandidate==3122) && (codeXi0Candidate==3322) ){
		return (TrackXi0Candidate->GetLabel());
	} else if((codeProtonCandidate==2212) && (codeLamdaCandidate==2212) ){
		Int_t codePotentialGrandMother = codeLamdaCandidate;
		Int_t labelPotentialGrandMother = mcTrackProtonFromXi0Candidate->GetMother();
		Int_t codePotentialGrandGrandMother = TrackXi0Candidate->GetLabel();
		Int_t counter = 0;
		while(codePotentialGrandMother == 2212 && counter < 10 && labelPotentialGrandMother < fAODMCTrackArray->GetEntriesFast()){
			AliAODMCParticle *TrackPotentialGrandMother = static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(labelPotentialGrandMother)); //mother MC particle object
			if(!TrackPotentialGrandMother) return -1;
			codePotentialGrandMother = TrackPotentialGrandMother->GetPdgCode();
			AliAODMCParticle *TrackPotentialGrandGrandMother = static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(TrackPotentialGrandMother->GetMother())); //mother MC particle object
			if(!TrackPotentialGrandGrandMother) return -1;
			codePotentialGrandGrandMother = TrackPotentialGrandGrandMother->GetPdgCode();
			if(codePotentialGrandMother == 3122 && codePotentialGrandGrandMother == 3322){
				return TrackPotentialGrandGrandMother->GetLabel();
			}	
			labelPotentialGrandMother = TrackPotentialGrandMother->GetMother();
			counter += 1;
		}
		return -1;
	} else {
		return -1;
	}
	return -1;
}
//________________________________________________________________________
Int_t AliAnalysisTaskSigmaPlToProtonPiZeroAOD::IsPhotonFromXi0(AliAODConversionPhoton *PhotonFromXi0Candidate, TClonesArray* fAODMCTrackArray, Int_t iCut, Double_t fWeightJetJetMC)
{ //checks if a reconstructed photon from sigma plus
	AliAODMCParticle *PhotonFromXi0 = NULL;
	if (PhotonFromXi0Candidate->GetNCaloPhotonMCLabels() > 0){
		// PhotonFromXi0 = PhotonFromXi0Candidate->GetMCParticle(fMCEvent);
		if (PhotonFromXi0Candidate->GetCaloPhotonMCLabel(0) < 0) return -1;
		PhotonFromXi0 = static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(PhotonFromXi0Candidate->GetCaloPhotonMCLabel(0)));
		if (PhotonFromXi0) {
			if(PhotonFromXi0->GetPdgCode() == 22){
				if (PhotonFromXi0->GetMother() < 0) return -1;
				AliAODMCParticle* pi0nFromXi0Candidate = static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(PhotonFromXi0->GetMother()));
				if (pi0nFromXi0Candidate->GetMother() < 0) return -1;
				AliAODMCParticle* Xi0Candidate = static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(pi0nFromXi0Candidate->GetMother()));
				if(pi0nFromXi0Candidate->GetPdgCode() == 111 && Xi0Candidate->GetPdgCode() == 3322){
					return (Xi0Candidate->GetLabel());	
				} else if(pi0nFromXi0Candidate->GetPdgCode() == 11 || pi0nFromXi0Candidate->GetPdgCode() == -11 || pi0nFromXi0Candidate->GetPdgCode() == 22){
					Int_t codeXi0 = pi0nFromXi0Candidate->GetPdgCode();
					Int_t codepionFromXi0 = PhotonFromXi0->GetPdgCode();
					Int_t labelXi0 = PhotonFromXi0->GetMother();
					Int_t counter = 0;
					while((codeXi0 == 11 || codeXi0 == -11 || codeXi0 == 22 || codeXi0 == 111) && counter < 10 && labelXi0 < fAODMCTrackArray->GetEntriesFast()){
						AliAODMCParticle *TrackGrandMotherXi0 = static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(labelXi0)); //mother MC particle object
						if(!TrackGrandMotherXi0) return -1;
						codeXi0 = TrackGrandMotherXi0->GetPdgCode();
						if(codepionFromXi0 == 111 && codeXi0 == 3322){
							return TrackGrandMotherXi0->GetLabel();
						}
						labelXi0 = TrackGrandMotherXi0->GetMother();
						codepionFromXi0 = codeXi0;
						counter += 1;
					}
					return -1;
				}else {
					return -1;
				}	
			} else if(PhotonFromXi0->GetPdgCode() == 11 || PhotonFromXi0->GetPdgCode() == -11){
				if (PhotonFromXi0->GetMother() < 0) return -1;
				AliAODMCParticle* photonFromXi0Candidate = static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(PhotonFromXi0->GetMother()));
				if (photonFromXi0Candidate->GetMother() < 0) return -1;
				AliAODMCParticle* pi0nFromXi0Candidate = static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(photonFromXi0Candidate->GetMother()));
				if (pi0nFromXi0Candidate->GetMother() < 0) return -1;
				AliAODMCParticle* Xi0Candidate = static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(pi0nFromXi0Candidate->GetMother()));
				if(photonFromXi0Candidate->GetPdgCode() == 22 && pi0nFromXi0Candidate->GetPdgCode() == 111  && Xi0Candidate->GetPdgCode() == 3322){
					return (Xi0Candidate->GetLabel());
				} else if((photonFromXi0Candidate->GetPdgCode() == 11 || photonFromXi0Candidate->GetPdgCode() == -11 || photonFromXi0Candidate->GetPdgCode() == 22) && (pi0nFromXi0Candidate->GetPdgCode() == 11 || pi0nFromXi0Candidate->GetPdgCode() == -11 || pi0nFromXi0Candidate->GetPdgCode() == 22)){
					Int_t codeGrandMother = photonFromXi0Candidate->GetPdgCode();
					Int_t codePotentialMother = PhotonFromXi0->GetPdgCode();
					Int_t labelGrandMother = PhotonFromXi0->GetMother();
					Int_t counter = 0;
					while((codeGrandMother == 11 || codeGrandMother == -11 || codeGrandMother == 22 || codeGrandMother == 111) && counter < 10 && labelGrandMother < fAODMCTrackArray->GetEntries()){
						AliAODMCParticle *TrackGrandMother = static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(labelGrandMother)); //mother MC particle object
						if(!TrackGrandMother) return -1;
						codeGrandMother = TrackGrandMother->GetPdgCode();
						if(codePotentialMother == 111 && codeGrandMother == 3322){
							return TrackGrandMother->GetLabel();
						}
						labelGrandMother = TrackGrandMother->GetMother();
						codePotentialMother = codeGrandMother;
						counter += 1;
					}
					return -1;
				}else {
					return -1;
				}
			} else {
				return -1;
			}
		}
	}
	return -1;

}
//_____________________________________________________________________________
void AliAnalysisTaskSigmaPlToProtonPiZeroAOD::Terminate(Option_t *)
{
	// terminate
	// called at the END of the analysis (when all events are processed)
}
//_____________________________________________________________________________
