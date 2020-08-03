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
	fEventHandler(0),
	fAODMCTrackArray(0),
	fOutputList(0),
	fAODList(NULL),
	fHistSigmaPlus(0),
	fHistSigmaPlusMC(0),
	fHistSigmaPlusMCTrueProtonGamma(0),
	fHistSigmaPlusMCTrueProton(0),
	fHistSigmaPlusMCTruePion(0),
	fHistSigmaPlusMCTrueProtonGammaDoubleCounting(0),
	fHistSigmaPlusMCTrueProtonDoubleCounting(0),
	fHistSigmaPlusMCTruePionDoubleCounting(0),
	fHistGenSigmaPt(0),
	fHistGenProtonPt(0),
	fHistGenPiZeroPt(0),
	fHistSigmaPtEta(0),
	fHistProtonPtEta(0),
	fHistPi0PtEta(0),
	fHistGenAngleProtonPiZero(0),
	fHistReconstructedMassPi0(0),
	fHistReconstructedMassPi0MC(0),
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
	fHistTrackDCAXY(0),
	fHistTrackDCAZ(0),
	fHistTrackDCAXYTrue(0),
	fHistTrackDCAZTrue(0),
	fHistDEDx(0),
	fHistTPCSignal(0),
	fHistoMotherBackInvMassPt(0),
	fHistSigmaMassPtWoPodCut(0),
	fHistSigmaMassPtWoPodCutMC(0),
	fFitPodolanskiUpperCut(0),
	fFitPi0MassDataLowPt(0),
	fFitPi0MassDataHighPt(0),
	fFitPi0MassMCLowPt(0),
	fFitPi0MassMCHighPt(0),
	fFitWidthData(0),
	fFitWidthMC(0),
	fV0Reader(NULL),
	fV0ReaderName("V0ReaderV1"),
	fCorrTaskSetting(""),
	fEventCutArray(NULL),
	fEventCuts(NULL),
	fClusterCutArray(NULL),
	fCaloPhotonCuts(NULL),                                    // List with Cluster Cuts
	fMesonCutArray(NULL),
	fMesonCuts(NULL),                                       // MesonCutObject
	fConversionCuts(NULL),                                  // ConversionPhotonCutObject                                     // If a jet is near the EMCal in the current event
	fHistNEvents(NULL),                                        //! array of histos with event information
	fnCuts(0),                                               // number of cuts to be analysed in parallel
	fIsHeavyIon(0),                                          // switch for pp = 0, PbPb = 1, pPb = 2
	fDoLightOutput(kFALSE),                                       // switch for running light output, kFALSE -> normal mode, kTRUE -> light mode
	fIsMC(0),
	fDoInOutTimingCluster(kFALSE),                                // manual timing cut for cluster to combine cluster within timing cut and without
	fMinTimingCluster(0),                                    // corresponding ranges, min
	fMaxTimingCluster(0),                                   // corresponding ranges, max
	fTrackMatcherRunningMode(0),
	fQTCutUpper(0),
	fQTCutLower(0)


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
	fEventHandler(0),
	fAODMCTrackArray(0),
	fOutputList(0),
	fAODList(NULL),
	fHistSigmaPlus(0),
	fHistSigmaPlusMC(0),
	fHistSigmaPlusMCTrueProtonGamma(0),
	fHistSigmaPlusMCTrueProton(0),
	fHistSigmaPlusMCTruePion(0),
	fHistSigmaPlusMCTrueProtonGammaDoubleCounting(0),
	fHistSigmaPlusMCTrueProtonDoubleCounting(0),
	fHistSigmaPlusMCTruePionDoubleCounting(0),
	fHistGenSigmaPt(0),
	fHistGenProtonPt(0),
	fHistGenPiZeroPt(0),
	fHistSigmaPtEta(0),
	fHistProtonPtEta(0),
	fHistPi0PtEta(0),
	fHistGenAngleProtonPiZero(0),
	fHistReconstructedMassPi0(0),
	fHistReconstructedMassPi0MC(0),
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
	fHistTrackDCAXY(0),
	fHistTrackDCAZ(0),
	fHistTrackDCAXYTrue(0),
	fHistTrackDCAZTrue(0),
	fHistDEDx(0),
	fHistTPCSignal(0),
	fHistoMotherBackInvMassPt(0),
	fHistSigmaMassPtWoPodCut(0),
	fHistSigmaMassPtWoPodCutMC(0),
	fFitPodolanskiUpperCut(0),
	fFitPi0MassDataLowPt(0),
	fFitPi0MassDataHighPt(0),
	fFitPi0MassMCLowPt(0),
	fFitPi0MassMCHighPt(0),
	fFitWidthData(0),
	fFitWidthMC(0),
	fV0Reader(NULL),
	fV0ReaderName("V0ReaderV1"),
	fCorrTaskSetting(""),
	fEventCutArray(NULL),
	fEventCuts(NULL),
	fClusterCutArray(NULL),
	fCaloPhotonCuts(NULL),                                    // List with Cluster Cuts
	fMesonCutArray(NULL),
	fMesonCuts(NULL),                                       // MesonCutObject
	fConversionCuts(NULL),                                  // ConversionPhotonCutObject                                     // If a jet is near the EMCal in the current event
	fHistNEvents(NULL),                                        //! array of histos with event information
	fnCuts(0),                                               // number of cuts to be analysed in parallel
	fIsHeavyIon(0),                                          // switch for pp = 0, PbPb = 1, pPb = 2
	fDoLightOutput(kFALSE),                                       // switch for running light output, kFALSE -> normal mode, kTRUE -> light mode
	fIsMC(0),
	fDoInOutTimingCluster(kFALSE),                                // manual timing cut for cluster to combine cluster within timing cut and without
	fMinTimingCluster(0),                                    // corresponding ranges, min
	fMaxTimingCluster(0),                                   // corresponding ranges, max
	fTrackMatcherRunningMode(0),
	fQTCutUpper(0),
	fQTCutLower(0)


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
	fHistTrackDCAXY = new TH2F*[fnCuts];
	fHistTrackDCAZ = new TH2F*[fnCuts];
	fHistDEDx = new TH2F*[fnCuts];
	fHistTPCSignal = new TH2F*[fnCuts];
	fHistoMotherBackInvMassPt = new TH2F*[fnCuts];
	fHistSigmaMassPtWoPodCut = new TH2F*[fnCuts];
	fHistNEvents = new TH1F*[fnCuts];
	fFitPodolanskiUpperCut = new TF1*[fnCuts];
	fFitPi0MassDataLowPt = new TF1*[fnCuts];
	fFitPi0MassDataHighPt = new TF1*[fnCuts];
	fFitWidthData = new TF1*[fnCuts];

	if(fIsMC > 0){
		fHistSigmaPlusMC = new TH2F*[fnCuts];
		fHistSigmaPlusMCTrueProtonGamma = new TH2F*[fnCuts];
		fHistSigmaPlusMCTrueProton = new TH2F*[fnCuts];
		fHistSigmaPlusMCTruePion = new TH2F*[fnCuts];
		fHistSigmaPlusMCTrueProtonGammaDoubleCounting = new TH2F*[fnCuts];
		fHistSigmaPlusMCTrueProtonDoubleCounting = new TH2F*[fnCuts];
		fHistSigmaPlusMCTruePionDoubleCounting = new TH2F*[fnCuts];
		fHistReconstructedMassPi0MC = new TH2F*[fnCuts];
		fHistSigmaMassPtWoPodCutMC = new TH2F*[fnCuts];
		fHistTrackDCAXYTrue = new TH2F*[fnCuts];
		fHistTrackDCAZTrue = new TH2F*[fnCuts];
		fFitPi0MassMCLowPt = new TF1*[fnCuts];
		fFitPi0MassMCHighPt = new TF1*[fnCuts];
		fFitWidthMC = new TF1*[fnCuts];

		fHistThetaPhiTrueSigmaPl = new TH2F("fHistThetaPhiTrueSigmaPl", "fHistThetaPhiTrueSigmaPl;#theta ; #phi", 50, -1., 4. ,100, 0., 2*TMath::Pi());
		fHistGenSigmaPt = new TH1F("fHistGenSigmaPt", "fHistGenSigmaPt;#it{p}_{T,Sigma} (GeV/#it{c});Yield", 100, 0, 30);
		fHistGenProtonPt = new TH1F("fHistGenProtonPt", "fHistGenProtonPt;#it{p}_{T,Proton} (GeV/#it{c});Yield", 100, 0, 30);
		fHistGenPiZeroPt = new TH1F("fHistGenPiZeroPt", "fHistGenPiZeroPt;#it{p}_{T,Pion} (GeV/#it{c});Yield", 100, 0, 10);
		fHistSigmaPtEta = new TH2F("fHistSigmaPtEta", "fHistSigmaPtEta;#it{p}_{T,Pion} (GeV/#it{c});Yield", 100, 0, 30, 10, -5, 5);
		fHistProtonPtEta = new TH2F("fHistProtonPtEta", "fHistProtonPtEta;#it{p}_{T,Pion} (GeV/#it{c});Yield", 100, 0, 30, 10, -5, 5);
		fHistPi0PtEta = new TH2F("fHistPi0PtEta", "fHistPi0PtEta;#it{p}_{T,Pion} (GeV/#it{c});Yield", 100, 0, 10, 10, -5, 5);
		fHistGenAngleProtonPiZero = new TH1F("fHistGenAngleProtonPiZero", "fHistGenAngleProtonPiZero;#it{beta}_{Proton,Pion} (rad);Yield", 10, 0, TMath::Pi());
		fHistPodolanskiGenTrue = new TH2F("fHistPodolanskiGenTrue","", 100, 0., 1., 100, 0., 1.);
	}

	for(Int_t iCut = 0; iCut<fnCuts;iCut++){

		TString cutstringEvent    = ((AliConvEventCuts*)fEventCutArray->At(iCut))->GetCutNumber();
		TString cutstringCalo     = ((AliCaloPhotonCuts*)fClusterCutArray->At(iCut))->GetCutNumber();
		TString cutstringMeson    = "NoMesonCut";

		fAODList[iCut]          = new TList();
		fAODList[iCut]->SetOwner(kTRUE);
		fAODList[iCut]->SetName(Form("%s_%s_%s AOD histograms", cutstringEvent.Data(), cutstringCalo.Data(), cutstringMeson.Data()));
		fOutputList->Add(fAODList[iCut]);

		fHistSigmaPlus[iCut] = new TH2F("fHistSigmaPlus", ";#it{m}_{inv} (GeV/#it{c^{2}});#it{p}_{T} (GeV/#it{c})", 40, 1.1, 1.3, 40, 0., 20.);
		fAODList[iCut]->Add(fHistSigmaPlus[iCut]);
		fHistReconstructedMassPi0[iCut] = new TH2F("fHistReconstructedMassPi0",";#it{m}_{inv} (GeV/#it{c^{2}});#it{p}_{T} (GeV/#it{c})", 60, 0., 0.3, 30, 0., 15.);
		fAODList[iCut]->Add(fHistReconstructedMassPi0[iCut]);
		fHistPodolanski[iCut] = new TH2F("fHistPodolanski","", 100, 0., 1., 100, 0., 1.);
		fAODList[iCut]->Add(fHistPodolanski[iCut]);
		fHistPodolanskiWCut[iCut] = new TH2F("fHistPodolanskiWCut","", 100, 0., 1., 100, 0., 1.);
		fAODList[iCut]->Add(fHistPodolanskiWCut[iCut]);
		fHistAllTracksPt[iCut] = new TH1F("fHistAllTracksPt", "fHistAllTracksPt;#it{p}_{T} (GeV/#it{c});Yield", 100, 0, 10);
		fAODList[iCut]->Add(fHistAllTracksPt[iCut]);
		fHistProtonPt[iCut] = new TH1F("fHistProtonPt", "fHistProtonPt;#it{p}_{T,Proton} (GeV/#it{c});Yield", 100, 0, 10);
		fAODList[iCut]->Add(fHistProtonPt[iCut]);
		fHistThetaPhi[iCut] = new TH2F("fHistThetaPhi", "fHistThetaPhi;#theta ; #phi", 50, 0., TMath::Pi() ,100, 0., 2*TMath::Pi());
		fAODList[iCut]->Add(fHistThetaPhi[iCut]);
		fHistThetaPhiProton[iCut] = new TH2F("fHistThetaPhiProton", "fHistThetaPhiProton;#theta_{p} ; #phi_{p}", 50, 0., TMath::Pi() ,100, 0., 2*TMath::Pi());
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
		fHistTrackDCAXY[iCut] = new TH2F("fHistTrackDCAXY", "fHistTrackDCAXY;#it{N}_{Protons per Event};#it{p}_{T} (GeV/#it{c})", 400, 0., 10.,20, 0., 10.);
		fAODList[iCut]->Add(fHistTrackDCAXY[iCut]);
		fHistTrackDCAZ[iCut] = new TH2F("fHistTrackDCAZ", "fHistTrackDCAZ;#it{N}_{Protons per Event};#it{p}_{T} (GeV/#it{c})", 400, 0., 10.,20, 0., 10.);
		fAODList[iCut]->Add(fHistTrackDCAZ[iCut]);
		fHistDEDx[iCut] = new TH2F("fHistDEDx", "fHistDEDx;#it{p};#d it{E}/d it{x}", 100,0.01,10.,100,1.,200.);
		fAODList[iCut]->Add(fHistDEDx[iCut]);
		fHistTPCSignal[iCut] = new TH2F("fHistTPCSignal", "fHistTPCSignal;#it{p}_{T};#sigma_{TPC}", 100, 0., 10., 20, -3., 3.);
		fAODList[iCut]->Add(fHistTPCSignal[iCut]);
		fHistSigmaMassPtWoPodCut[iCut] = new TH2F("fHistSigmaMassPtWoPodCut", ";#it{m}_{inv} (GeV/#it{c^{2}});#it{p}_{T} (GeV/#it{c})", 20, 1.1, 1.3, 20, 0., 20.);
		fAODList[iCut]->Add(fHistSigmaMassPtWoPodCut[iCut]);
		fHistoMotherBackInvMassPt[iCut] = new TH2F("fHistoMotherBackInvMassPt", ";#it{m}_{inv} (GeV/#it{c^{2}});#it{p}_{T} (GeV/#it{c})", 20, 1.1, 1.3, 20, 0., 20.);
		fHistoMotherBackInvMassPt[iCut]->Sumw2();
		fAODList[iCut]->Add(fHistoMotherBackInvMassPt[iCut]);


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

		fAODList[iCut]->Add(fHistNEvents[iCut]);          // don't forget to add it to the list! the list will be written to file, so if you want

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


		fFitPi0MassDataHighPt[iCut] = new TF1("fFitPi0MassDataHighPt","[0]*1/x*log(x*x)+[1]*1/(x*x)+[2]",1.1,25.);
		fFitPi0MassDataLowPt[iCut] = new TF1("fFitPi0MassDataLowPt","[0]*x*x+[1]*x+[2]",0.3,1.2);
		fFitWidthData[iCut] = new TF1("fFitWidthData","[0]*exp(x*[1]+[2])+[3]*x+ [4] + [5]",0.3,25.);

		if(((AliCaloPhotonCuts*)fClusterCutArray->At(iCut))->GetClusterType() == 2){
			fFitPi0MassDataHighPt[iCut] -> SetParameters(-7.82053e-04, -1.44012e-03, 1.35284e-01);
			fFitPi0MassDataLowPt[iCut] -> SetParameters(1.38489e-02, -3.04366e-02, 1.50638e-01);
			fFitWidthData[iCut] -> SetParameters(1.23693e+01,-7.18506e-01, -7.61979e+00, 7.47276e-05, 2.14652e-03, 2.14652e-03);
		}
		if(((AliCaloPhotonCuts*)fClusterCutArray->At(iCut))->GetClusterType() == 1 || ((AliCaloPhotonCuts*)fClusterCutArray->At(iCut))->GetClusterType() == 3 || ((AliCaloPhotonCuts*)fClusterCutArray->At(iCut))->GetClusterType() == 4){
			fFitPi0MassDataLowPt[iCut] -> SetParameters(1.22791e-02, -2.92736e-02, 1.52764e-01);
			fFitPi0MassDataHighPt[iCut] -> SetParameters(-3.44079e-03, -2.72693e-03, 1.38021e-01);
			fFitWidthData[iCut] -> SetParameters(2.52112e+00,-8.52979e-01, -5.58870e+00, 4.57756e-04, 1.91358e-01, -1.85394e-01);
		}

		fFitPodolanskiUpperCut[iCut] = new TF1("fFitPodolanskiUpperCut","sqrt([0]*[0]*(1-(((x-0.61)*(x-0.61))/(0.41*0.41))))",0.2,1.02);
		fFitPodolanskiUpperCut[iCut]->SetParameter(0, fQTCutUpper);
		
		fAODList[iCut]->Add(fFitPodolanskiUpperCut[iCut]);
		fAODList[iCut]->Add(fFitPi0MassDataHighPt[iCut]);
		fAODList[iCut]->Add(fFitWidthData[iCut]);
		fAODList[iCut]->Add(fFitPi0MassDataLowPt[iCut]);

		if(fIsMC > 0){
			fHistPodolanskiWCutTrue[iCut] = new TH2F("fHistPodolanskiWCutTrue","", 100, 0., 1., 100, 0., 1.);

			fHistTrueProtonPt[iCut] = new TH1F("fHistTrueProtonPt", "fHistTrueProtonPt;#it{p}_{T,Proton} (GeV/#it{c});Yield", 100, 0, 10);
			fHistSigmaMassPtWoPodCutMC[iCut] = new TH2F("fHistSigmaMassPtWoPodCutMC", ";#it{m}_{inv} (GeV/#it{c^{2}});#it{p}_{T} (GeV/#it{c})", 20, 1.1, 1.3, 20, 0., 20.);
			fHistTrackDCAXYTrue[iCut] = new TH2F("fHistTrackDCAXYTrue", "fHistTrackDCAXYTrue;#it{N}_{Protons per Event};#it{p}_{T} (GeV/#it{c})", 400, 0., 10.,50, 0., 10.);
			fHistTrackDCAZTrue[iCut] = new TH2F("fHistTrackDCAZTrue", "fHistTrackDCAZTrue;#it{N}_{Protons per Event};#it{p}_{T} (GeV/#it{c})", 400, 0., 10.,50, 0., 10.);

			fFitPi0MassMCHighPt[iCut] = new TF1("fFitPi0MassMCHighPt","[0]*1/x*log(x*x)+[1]*1/(x*x)+[2]",1.,25.);
			fFitPi0MassMCLowPt[iCut] = new TF1("fFitPi0MassMCLowPt","[0]*x*x+[1]*x+[2]",0.3,1.05);
			fFitWidthMC[iCut] = new TF1("fFitWidthMC","[0]*exp(x*[1]+[2])+[3]*x+[4] +[5]",0.3,25.);

			if(((AliCaloPhotonCuts*)fClusterCutArray->At(iCut))->GetClusterType() == 2){
				fFitPi0MassMCHighPt[iCut] -> SetParameters(-1.17878e-03, -1.90910e-03,  1.35525e-01);
				fFitPi0MassMCLowPt[iCut] -> SetParameters(3.12765e-02, -6.48323e-02, 1.67346e-01);
				fFitWidthMC[iCut] -> SetParameters(2.00567e+01,-9.40838e-01, -8.18957e+00, 3.31007e-05, 2.25791e-03, 2.79006e-03);
			}
			if(((AliCaloPhotonCuts*)fClusterCutArray->At(iCut))->GetClusterType() == 1 || ((AliCaloPhotonCuts*)fClusterCutArray->At(iCut))->GetClusterType() == 3 || ((AliCaloPhotonCuts*)fClusterCutArray->At(iCut))->GetClusterType() == 4){
				fFitPi0MassMCLowPt[iCut] -> SetParameters(5.96171e-03, -1.52107e-02, 1.43861e-01);
				fFitPi0MassMCHighPt[iCut] -> SetParameters(-5.33046e-03, -4.49031e-03, 1.38791e-01);
				fFitWidthMC[iCut] -> SetParameters(3.45251e-01,-8.54696e-01, -3.44592e+00, 5.65500e-04, 2.87614e-03, 3.28479e-03);
			}
			fAODList[iCut]->Add(fFitPi0MassMCHighPt[iCut]);
			fAODList[iCut]->Add(fFitPi0MassMCLowPt[iCut]);
			fAODList[iCut]->Add(fFitWidthMC[iCut]);

			fHistSigmaPlusMC[iCut] = new TH2F("fHistSigmaPlusMC", ";#it{m}_{inv} (GeV/#it{c^{2}});#it{p}_{T} (GeV/#it{c})", 40, 1.1, 1.3, 40, 0., 20.);
			fHistSigmaPlusMCTrueProtonGamma[iCut] = new TH2F("fHistSigmaPlusMCTrueProtonGamma", ";#it{m}_{inv} (GeV/#it{c^{2}});#it{p}_{T} (GeV/#it{c})", 40, 1.1, 1.3, 40, 0., 20.);
			fHistSigmaPlusMCTrueProton[iCut] = new TH2F("fHistSigmaPlusMCTrueProton", ";#it{m}_{inv} (GeV/#it{c^{2}});#it{p}_{T} (GeV/#it{c})", 40, 1.1, 1.3, 40, 0., 20.);
			fHistSigmaPlusMCTruePion[iCut] = new TH2F("fHistSigmaPlusMCTruePion", ";#it{m}_{inv} (GeV/#it{c^{2}});#it{p}_{T} (GeV/#it{c})", 40, 1.1, 1.3, 40, 0., 20.);
			fHistSigmaPlusMCTrueProtonGammaDoubleCounting[iCut] = new TH2F("fHistSigmaPlusMCTrueProtonGammaDoubleCounting", ";#it{m}_{inv} (GeV/#it{c^{2}});#it{p}_{T} (GeV/#it{c})", 40, 1.1, 1.3, 40, 0., 20.);
			fHistSigmaPlusMCTrueProtonDoubleCounting[iCut] = new TH2F("fHistSigmaPlusMCTrueProtonDoubleCounting", ";#it{m}_{inv} (GeV/#it{c^{2}});#it{p}_{T} (GeV/#it{c})", 40, 1.1, 1.3, 40, 0., 20.);
			fHistSigmaPlusMCTruePionDoubleCounting[iCut] = new TH2F("fHistSigmaPlusMCTruePionDoubleCounting", ";#it{m}_{inv} (GeV/#it{c^{2}});#it{p}_{T} (GeV/#it{c})", 40, 1.1, 1.3, 40, 0., 20.);
			fHistReconstructedMassPi0MC[iCut] = new TH2F("fHistReconstructedMassPi0MC",";#it{m}_{inv} (GeV/#it{c^{2}});#it{p}_{T} (GeV/#it{c})", 60, 0., 0.3, 30, 0., 15.);
			fAODList[iCut]->Add(fHistSigmaPlusMC[iCut]);
			fAODList[iCut]->Add(fHistSigmaPlusMCTrueProtonGamma[iCut]);
			fAODList[iCut]->Add(fHistSigmaPlusMCTrueProton[iCut]);
			fAODList[iCut]->Add(fHistSigmaPlusMCTruePion[iCut]);
			fAODList[iCut]->Add(fHistSigmaPlusMCTrueProtonGammaDoubleCounting[iCut]);
			fAODList[iCut]->Add(fHistSigmaPlusMCTrueProtonDoubleCounting[iCut]);
			fAODList[iCut]->Add(fHistSigmaPlusMCTruePionDoubleCounting[iCut]);
			fAODList[iCut]->Add(fHistPodolanskiGenTrue);
			fAODList[iCut]->Add(fHistPodolanskiWCutTrue[iCut]);
			fAODList[iCut]->Add(fHistTrueProtonPt[iCut]);
			fAODList[iCut]->Add(fHistTrackDCAXYTrue[iCut]);
			fAODList[iCut]->Add(fHistTrackDCAZTrue[iCut]);
			fAODList[iCut]->Add(fHistThetaPhiTrueSigmaPl);
			fAODList[iCut]->Add(fHistGenSigmaPt);
			fAODList[iCut]->Add(fHistGenProtonPt);
			fAODList[iCut]->Add(fHistGenPiZeroPt);
			fAODList[iCut]->Add(fHistSigmaPtEta);
			fAODList[iCut]->Add(fHistProtonPtEta);
			fAODList[iCut]->Add(fHistPi0PtEta);
			fAODList[iCut]->Add(fHistGenAngleProtonPiZero);
			fAODList[iCut]->Add(fHistReconstructedMassPi0MC[iCut]);
			fAODList[iCut]->Add(fHistSigmaMassPtWoPodCutMC[iCut]);

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
			FillfHistNEvents(iCut,eventQuality);
			// if (fIsMC>1) fHistoNEventsWOWeight[iCut]->Fill(eventQuality);
		}
		return;
	}

	// Generiertes Spektrum der Sigma+
	if(fIsMC > 0){

		fAODMCTrackArray = dynamic_cast<TClonesArray*>(fInputEvent->FindListObject(AliAODMCParticle::StdBranchName()));
		if (fAODMCTrackArray == NULL) return;

		for(Long_t i = 1; i < fAODMCTrackArray->GetEntriesFast(); ++i) {
			AliAODMCParticle* sigma = nullptr;
			sigma = static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(i));
			if (!sigma) {
				Printf("ERROR: Could not retrieve AliAODMCParticle");
				continue;
			}
			if(sigma->GetPdgCode() == 3222){
				if(fHistThetaPhiTrueSigmaPl) fHistThetaPhiTrueSigmaPl->Fill(sigma->Theta(), sigma->Phi());
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
						if(fHistGenSigmaPt) fHistGenSigmaPt->Fill(sigma->Pt());
						if(fHistGenProtonPt) fHistGenProtonPt->Fill(protonTrueVec.Pt());
						if(fHistGenPiZeroPt) fHistGenPiZeroPt->Fill(pionTrueVec.Pt());
						if(fHistSigmaPtEta) fHistSigmaPtEta->Fill(sigma->Pt(), sigma->Eta());
						if(fHistProtonPtEta) fHistProtonPtEta->Fill(protonTrueVec.Pt(), protonTrueVec.Eta());
						if(fHistPi0PtEta) fHistPi0PtEta->Fill(pionTrueVec.Pt(), pionTrueVec.Eta());
						if(fHistGenAngleProtonPiZero) fHistGenAngleProtonPiZero->Fill(pionTrueVec.Angle(protonTrueVec.Vect()));
						if(fHistPodolanskiGenTrue) fHistPodolanskiGenTrue->Fill(GetPodAlpha(sigmaTrueVec, protonTrueVec, pionTrueVec),GetQT(sigmaTrueVec, pionTrueVec));
					}
					if((daughter1->GetPdgCode() == 111) && (daughter2->GetPdgCode() == 2212)){
						pionTrueVec.SetPtEtaPhiM(daughter1->Pt(), daughter1->Eta(), daughter1->Phi(), daughter1->M());
						protonTrueVec.SetPtEtaPhiM(daughter2->Pt(), daughter2->Eta(), daughter2->Phi(), daughter2->M());
						if(fHistGenSigmaPt) fHistGenSigmaPt->Fill(sigma->Pt());
						if(fHistGenProtonPt) fHistGenProtonPt->Fill(protonTrueVec.Pt());
						if(fHistGenPiZeroPt) fHistGenPiZeroPt->Fill(pionTrueVec.Pt());
						if(fHistSigmaPtEta) fHistSigmaPtEta->Fill(sigma->Pt(), sigma->Eta());
						if(fHistProtonPtEta) fHistProtonPtEta->Fill(protonTrueVec.Pt(), protonTrueVec.Eta());
						if(fHistPi0PtEta) fHistPi0PtEta->Fill(pionTrueVec.Pt(), pionTrueVec.Eta());
						if(fHistGenAngleProtonPiZero) fHistGenAngleProtonPiZero->Fill(pionTrueVec.Angle(protonTrueVec.Vect()));
						if(fHistPodolanskiGenTrue) fHistPodolanskiGenTrue->Fill(GetPodAlpha(sigmaTrueVec, protonTrueVec, pionTrueVec),GetQT(sigmaTrueVec, pionTrueVec));
					}
				}
			}
		}
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
		// vector < TLorentzVector > pions;

		Bool_t isRunningEMCALrelAna = kFALSE;
		if (((AliCaloPhotonCuts*)fClusterCutArray->At(iCut))->GetClusterType() == 1
				|| ((AliCaloPhotonCuts*)fClusterCutArray->At(iCut))->GetClusterType() == 3
				|| ((AliCaloPhotonCuts*)fClusterCutArray->At(iCut))->GetClusterType() == 4
		   ) isRunningEMCALrelAna = kTRUE;

		Int_t eventNotAccepted = ((AliConvEventCuts*)fEventCutArray->At(iCut))->IsEventAcceptedByCut(fV0Reader->GetEventCuts(),fEvent,fMCEvent,fIsHeavyIon, isRunningEMCALrelAna);

		if(eventNotAccepted != 0){ // Check Centrality, PileUp, SDD and V0AND --> Not Accepted => eventQuality = 1
			// if (fIsMC>1) if(fHistoNEventsWOWeight[iCut]) fHistoNEventsWOWeight[iCut]->Fill(eventNotAccepted);
			// 	// cout << "event rejected due to wrong trigger: " <<eventNotAccepted << endl;
			// if (eventNotAccepted==3 && fIsMC > 0){
			//   triggered = kFALSE;
			// }else {
			//   continue;
			// }
			if (eventNotAccepted==3 && fIsMC > 0){
				FillfHistNEvents(iCut,0);
			} else {
				FillfHistNEvents(iCut,eventNotAccepted);
				continue;
			}
		} else {
			FillfHistNEvents(iCut,eventNotAccepted);
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
			if (!track->P()) continue;
			double protonSignalTPC = 0.; //numbers of sigmas of TPC signal
			double protonSignalTOF = 0.; //numbers of sigmas of TOF signal
			if((fPIDResponse->NumberOfSigmasTPC(track, AliPID::kProton))) protonSignalTPC = (fPIDResponse->NumberOfSigmasTPC(track, AliPID::kProton));
			if((fPIDResponse->NumberOfSigmasTOF(track, AliPID::kProton))) protonSignalTOF = (fPIDResponse->NumberOfSigmasTOF(track, AliPID::kProton)); //numbers of sigmas of TOF signal
			Float_t dcaXY = 0.0, dcaZ = 0.0;
			track->GetImpactParameters(dcaXY,dcaZ);
			// Select Proton Candidates
			// if ((TMath::Abs(protonSignalTPC) - 3. < 0.)) {
			if ((TMath::Abs(protonSignalTPC) - 3. < 0.) && (TMath::Abs(protonSignalTOF) - 3. < 0.) && (dcaXY > 0.05) && (dcaZ > 0.05)) {
				if(fHistProtonPt[iCut]) fHistProtonPt[iCut]->Fill(track->Pt());
				if(fHistThetaPhiProton[iCut]) fHistThetaPhiProton[iCut]->Fill(track->Theta(), track->Phi());
				if(fHistTPCSignal[iCut]) fHistTPCSignal[iCut]->Fill(track->Pt(), protonSignalTPC);
				// cout <<"Track ID \t" << i << "\t TPC Sigma \t" << protonSignalTPC << "\t P \t" <<  track->P() << "\t TPC Signal \t" << track->GetTPCsignal() << endl;
				proton.push_back(track);
			} else {
				tracks.push_back(track);
			}
			if(fHistAllTracksPt[iCut]) fHistAllTracksPt[iCut]->Fill(track->Pt());
			if(fHistThetaPhi[iCut]) fHistThetaPhi[iCut]->Fill(track->Theta(), track->Phi());
			//dE/dx Plot
			Double_t ptot = track->P();
			Double_t tpcSignal = track->GetTPCsignal();
			if(fHistDEDx[iCut]) fHistDEDx[iCut]->Fill(ptot, tpcSignal);
		}
		fHistNProtonsPerEvent[iCut]->Fill(proton.size());
		if(proton.size() > 0){
			//Find gammas
			Int_t nclus                     = 0;
			Int_t nClusWCuts                = 0;
			nclus = fAOD->GetNumberOfCaloClusters();
			if(fHistNClusWoCuts[iCut]) fHistNClusWoCuts[iCut]->Fill(nclus);
			if(nclus == 0)  continue;
			AliVCluster* clus = NULL;
			for(Int_t i=0; i < nclus; ++i){
				clus = new AliAODCaloCluster(*(AliAODCaloCluster*)fEvent->GetCaloCluster(i));
				if(!clus) {
					Printf("ERROR: Could not find clus %i",i);
					continue;
				}
				if(fHistClusterEWOCuts[iCut]) fHistClusterEWOCuts[iCut]-> Fill(clus->E());
				if((((AliCaloPhotonCuts*)fClusterCutArray->At(iCut))->ClusterIsSelected(clus,fAOD,fMCEvent,fIsMC, 1.,i)) == kTRUE){
					photon.push_back(clus);
					if(fHistClusterE[iCut]) fHistClusterE[iCut]-> Fill(clus->E());
					nClusWCuts+= 1;
				} else {
					delete clus;
				}
				clus = NULL;
			}
			if(fHistNClusWCuts[iCut]) fHistNClusWCuts[iCut]->Fill(nClusWCuts);
			//Paaring
			//			TLorentzVector photonCandidate1;
			//			TLorentzVector photonCandidate2;
			TLorentzVector protonVektor;
			TLorentzVector sigmaVektor;
			Bool_t trueSigmaProton = kFALSE;
			Bool_t trueSigmaPhoton1 = kFALSE;
			Bool_t trueSigmaPhoton2 = kFALSE;
			Bool_t mesonSelected = kFALSE;

			if( photon.size() > 1){
				for(unsigned int iProton = 0; iProton < proton.size(); ++iProton){
					trueSigmaProton = kFALSE;
					if(!proton[iProton]) {
						Printf("ERROR: Could not find proton[%i][%i]",iCut,iProton);
						continue;
					}
					AliAODTrack* protonCandidate = proton[iProton];
					protonVektor.SetPtEtaPhiM(protonCandidate->Pt(),protonCandidate->Eta(),protonCandidate->Phi(), 0.938272);
					Float_t dcaXYTrue = 0.0, dcaZTrue = 0.0;
					protonCandidate->GetImpactParameters(dcaXYTrue,dcaZTrue);
					if(fIsMC > 0){
						trueSigmaProton = IsRealProton(protonCandidate, fAODMCTrackArray);
						if(trueSigmaProton == kTRUE){
							if(fHistTrueProtonPt[iCut]) fHistTrueProtonPt[iCut]->Fill(protonCandidate->Pt());
						}
					}
					for(unsigned int iPhoton1 = 0; iPhoton1 < photon.size(); ++iPhoton1) {
						trueSigmaPhoton1 = kFALSE;
						if(!photon[iPhoton1]) {
							Printf("ERROR: Could not find photon[%i][%i]",iCut,iPhoton1);
							continue;
						}
						vector < vector < Double_t > > vDoubleCounting;
						vector < Double_t > vDoubleCountingCandidates;
						vector < Double_t > vCandidates;
						Bool_t bAllDaughtersMeasured = kFALSE;
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
									if (k<50)PhotonCandidate1.SetCaloPhotonMCLabel(k,mclabelsCluster[k]);
									// Int_t pdgCode = fMCEvent->Particle(mclabelsCluster[k])->GetPdgCode();
									// cout << "label " << k << "\t" << mclabelsCluster[k] << " pdg code: " << pdgCode << endl;
								}
							}
							trueSigmaPhoton1 = IsRealPhoton(&PhotonCandidate1, fAODMCTrackArray);
						}
						for(unsigned int iPhoton2 = 0; iPhoton2 < photon.size(); ++iPhoton2) {
							if( iPhoton2 > iPhoton1){
								trueSigmaPhoton2 = kFALSE;
								if(!photon[iPhoton2]) {
									Printf("ERROR: Could not find photon[%i][%i]",iCut,iPhoton2);
									continue;
								}
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
											if (k<50)PhotonCandidate2.SetCaloPhotonMCLabel(k,mclabelsCluster[k]);
											// Int_t pdgCode = fMCEvent->Particle(mclabelsCluster[k])->GetPdgCode();
											// cout << "label " << k << "\t" << mclabelsCluster[k] << " pdg code: " << pdgCode << endl;
										}
									}
									trueSigmaPhoton2 = IsRealPhoton(&PhotonCandidate2, fAODMCTrackArray);
								}
								AliAODConversionMother pi0cand = AliAODConversionMother(&PhotonCandidate1,&PhotonCandidate2);
								pi0cand.SetLabels(iPhoton1,iPhoton2);
								if((((AliConversionMesonCuts*)fMesonCutArray->At(iCut))->MesonIsSelected(&pi0cand,kTRUE,((AliConvEventCuts*)fEventCutArray->At(iCut))->GetEtaShift(),PhotonCandidate1.GetLeadingCellID(),PhotonCandidate2.GetLeadingCellID()))){
									if((fHistReconstructedMassPi0[iCut]) && (iProton == 0)) fHistReconstructedMassPi0[iCut]->Fill(pi0cand.M(), pi0cand.Pt());
									if((trueSigmaPhoton1 == kTRUE) && (trueSigmaPhoton2 ==kTRUE) && (fIsMC > 0)){
										if((fHistReconstructedMassPi0MC[iCut]) && (iProton == 0)) fHistReconstructedMassPi0MC[iCut]->Fill(pi0cand.M(), pi0cand.Pt());
									}
									mesonSelected = kFALSE;
									if(fIsMC > 0 ) {
										mesonSelected = IsPi0SelectedMC(&pi0cand, fFitPi0MassMCLowPt[iCut], fFitPi0MassMCHighPt[iCut], fFitWidthData[iCut]);
									}
									if(fIsMC == 0 ) {
										mesonSelected = IsPi0Selected(&pi0cand, fFitPi0MassDataLowPt[iCut], fFitPi0MassDataHighPt[iCut], fFitWidthData[iCut]);
									}
									if(mesonSelected == kTRUE){
										TLorentzVector rekombinatedPi0;
										rekombinatedPi0.SetPtEtaPhiM(pi0cand.Pt(), pi0cand.Eta(), pi0cand.Phi(), 0.135);
										// if((iProton == 0)) pions.push_back(rekombinatedPi0);
										sigmaVektor = protonVektor + rekombinatedPi0;
										if(fHistPodolanski[iCut]) fHistPodolanski[iCut]->Fill(GetPodAlpha(sigmaVektor, protonVektor, rekombinatedPi0),GetQT(sigmaVektor, rekombinatedPi0));
										if((trueSigmaProton == kTRUE) && (trueSigmaPhoton1 == kTRUE) && (trueSigmaPhoton2 ==kTRUE) && (fIsMC > 0)){
											if(fHistPodolanskiWCutTrue[iCut]){
												fHistPodolanskiWCutTrue[iCut]->Fill(GetPodAlpha(sigmaVektor, protonVektor, rekombinatedPi0),GetQT(sigmaVektor, rekombinatedPi0));
												fHistSigmaMassPtWoPodCutMC[iCut]->Fill(sigmaVektor.M(), sigmaVektor.Pt());
											}
										}
										fHistSigmaMassPtWoPodCut[iCut]->Fill(sigmaVektor.M(), sigmaVektor.Pt());
										if(GetPodAlpha(sigmaVektor, protonVektor, rekombinatedPi0) > 1.02 || GetPodAlpha(sigmaVektor, protonVektor, rekombinatedPi0) < 0.2) continue;
										if(GetQT(sigmaVektor, protonVektor) > fFitPodolanskiUpperCut[iCut]->Eval(GetPodAlpha(sigmaVektor, protonVektor, rekombinatedPi0))) continue;
										if(rekombinatedPi0.Angle(protonVektor.Vect()) > 0.35 || rekombinatedPi0.Angle(protonVektor.Vect()) < 0.05) continue;
										if(fHistPodolanskiWCut[iCut]) fHistPodolanskiWCut[iCut]->Fill(GetPodAlpha(sigmaVektor, protonVektor, rekombinatedPi0),GetQT(sigmaVektor, rekombinatedPi0));
										if(fHistSigmaPlus[iCut]) fHistSigmaPlus[iCut]-> Fill(sigmaVektor.M(), sigmaVektor.Pt());
										if(fHistTrackDCAXY[iCut])fHistTrackDCAXY[iCut]->Fill(TMath::Abs(dcaXYTrue),protonVektor.Pt());
										if(fHistTrackDCAZ[iCut])fHistTrackDCAZ[iCut]->Fill(TMath::Abs(dcaZTrue),protonVektor.Pt());
										if((trueSigmaProton == kTRUE) && (trueSigmaPhoton1 == kTRUE) && (trueSigmaPhoton2 ==kTRUE) && (fIsMC > 0)){
											if(fHistSigmaPlusMC[iCut]) fHistSigmaPlusMC [iCut]-> Fill(sigmaVektor.M(), sigmaVektor.Pt());
											if(fHistTrackDCAXYTrue[iCut]) fHistTrackDCAXYTrue[iCut]->Fill(std::abs(dcaXYTrue),protonCandidate->Pt());
											if(fHistTrackDCAZTrue[iCut]) fHistTrackDCAZTrue[iCut]->Fill(std::abs(dcaZTrue),protonCandidate->Pt());
											bAllDaughtersMeasured = kTRUE;
										}
										if(((trueSigmaProton == kTRUE) && (trueSigmaPhoton1 == kTRUE) && (trueSigmaPhoton2 ==kFALSE) && (fIsMC > 0)) || ((trueSigmaProton == kTRUE) && (trueSigmaPhoton1 == kFALSE) && (trueSigmaPhoton2 ==kTRUE) && (fIsMC > 0))){
											if(fHistSigmaPlusMCTrueProtonGamma[iCut]) fHistSigmaPlusMCTrueProtonGamma [iCut]-> Fill(sigmaVektor.M(), sigmaVektor.Pt());
											vDoubleCountingCandidates.push_back(sigmaVektor.M());
											vDoubleCountingCandidates.push_back(sigmaVektor.Pt());
											vDoubleCountingCandidates.push_back(1.);
											vDoubleCounting.push_back(vDoubleCountingCandidates);
										}
										if((trueSigmaProton == kTRUE) && (trueSigmaPhoton1 == kFALSE) && (trueSigmaPhoton2 ==kFALSE) && (fIsMC > 0)){
											if(fHistSigmaPlusMCTrueProton[iCut]) fHistSigmaPlusMCTrueProton [iCut]-> Fill(sigmaVektor.M(), sigmaVektor.Pt());
											vDoubleCountingCandidates.push_back(sigmaVektor.M());
											vDoubleCountingCandidates.push_back(sigmaVektor.Pt());
											vDoubleCountingCandidates.push_back(2.);
											vDoubleCounting.push_back(vDoubleCountingCandidates);
										}
										if((trueSigmaProton == kFALSE) && (trueSigmaPhoton1 == kTRUE) && (trueSigmaPhoton2 ==kTRUE) && (fIsMC > 0)){
											if(fHistSigmaPlusMCTruePion[iCut]) fHistSigmaPlusMCTruePion [iCut]-> Fill(sigmaVektor.M(), sigmaVektor.Pt());
											vDoubleCountingCandidates.push_back(sigmaVektor.M());
											vDoubleCountingCandidates.push_back(sigmaVektor.Pt());
											vDoubleCountingCandidates.push_back(3.);
											vDoubleCounting.push_back(vDoubleCountingCandidates);
										}
									}
								}
							}
						}
						if(bAllDaughtersMeasured == kTRUE){
							for(unsigned int i=0 ; i < vDoubleCounting.size(); ++i){
								vCandidates = vDoubleCounting[i];
								if(vCandidates[2] == 1.){
									if(fHistSigmaPlusMCTrueProtonGammaDoubleCounting[iCut]) fHistSigmaPlusMCTrueProtonGammaDoubleCounting[iCut]-> Fill(vCandidates[0], vCandidates[1]);
								}
								if(vCandidates[2] == 2.){
									if(fHistSigmaPlusMCTrueProtonDoubleCounting[iCut]) fHistSigmaPlusMCTrueProtonDoubleCounting[iCut]-> Fill(vCandidates[0], vCandidates[1]);
								}
								if(vCandidates[2] == 3.){
									if(fHistSigmaPlusMCTruePionDoubleCounting[iCut]) fHistSigmaPlusMCTruePionDoubleCounting[iCut]-> Fill(vCandidates[0], vCandidates[1]);
								}
							}
						}
						vCandidates.clear();
						vDoubleCounting.clear();
						vDoubleCountingCandidates.clear();
					}
				}
			}
		}

		CalculateBackgroundSwapp(photon, proton, vpos, iCut);

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
		// pions.clear();
	}
	// continue until all the tracks are processed
	PostData(1, fOutputList);                           // stream the results the analysis of this event to
	// the output manager which will take care of writing
	// it to a file
}
//_____________________________________________________________________________
Bool_t AliAnalysisTaskSigmaPlToProtonPiZeroAOD::IsPi0Selected(AliAODConversionMother* pi0CandTmp, TF1* fitPi0MassDataLowPt, TF1* fitPi0MassDataHighPt, TF1* fitWidthData)
{
	Double_t minPt = 0.3;
	// Double_t maxPtMassDataFirst = 1.2;
	// Double_t minPtMassDataFirst = 1.1;
	Double_t maxPt = 25.;

	Double_t minPi0Mass = 0;
	Double_t maxPi0Mass = 0;
	// Double_t width = 0;
	// Double_t meanMass = 0;

	if ( ((pi0CandTmp->Pt()) > minPt) && ((pi0CandTmp->Pt()) < maxPt)) {
		// 	width = fitWidthData->Eval(pi0CandTmp->Pt());
		// 	if((pi0CandTmp->Pt()) <= minPtMassDataFirst ){
		// 		meanMass = fitPi0MassDataLowPt->Eval(pi0CandTmp->Pt());
		// 	} else if((pi0CandTmp->Pt()) > minPtMassDataFirst){
		// 		meanMass = fitPi0MassDataHighPt->Eval(pi0CandTmp->Pt());
		// 	} else {
		// 		return kFALSE;
		// 	}
		minPi0Mass = 0.118;
		maxPi0Mass = 0.148;
		  // minPi0Mass = 0.
		// maxPi0Mass = 0.3;
		if((pi0CandTmp->M() < maxPi0Mass) && (pi0CandTmp->M() > minPi0Mass)){
			return kTRUE;
		}
	}
	return kFALSE;
}
//_____________________________________________________________________________
Bool_t AliAnalysisTaskSigmaPlToProtonPiZeroAOD::IsPi0SelectedMC(AliAODConversionMother *pi0CandTmp, TF1* fitPi0MassMCLowPt, TF1* fitPi0MassMCHighPt, TF1* fitWidthMC)
{
	Double_t minPt = 0.3;
	// Double_t maxPtMassMCFirst = 1.05;
	// Double_t minPtMassMCSecond = 1.;
	Double_t maxPt = 25.;

	Double_t minPi0Mass = 0;
	Double_t maxPi0Mass = 0;
	// Double_t width = 0;
	// Double_t meanMass = 0;

	if ( ((pi0CandTmp->Pt()) > minPt) && ((pi0CandTmp->Pt()) < maxPt)) {
		// width = fitWidthMC->Eval(pi0CandTmp->Pt());
		// if((pi0CandTmp->Pt()) <= minPtMassMCSecond ){
		// 	meanMass = fitPi0MassMCLowPt->Eval(pi0CandTmp->Pt());
		// } else if((pi0CandTmp->Pt()) > minPtMassMCSecond ){
		// 	meanMass = fitPi0MassMCHighPt->Eval(pi0CandTmp->Pt());
		// } else {
		// 	return kFALSE;
		// }
		minPi0Mass = 0.118;
		maxPi0Mass = 0.148;
	 // 	 minPi0Mass = 0.
		// maxPi0Mass = 0.3;
		// minPi0Mass = meanMass - 3*width;
		// maxPi0Mass = meanMass + 5*width;
		if(pi0CandTmp->M() < maxPi0Mass && pi0CandTmp->M() > minPi0Mass){
			return kTRUE;
		}
	}

	return kFALSE;
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
Bool_t AliAnalysisTaskSigmaPlToProtonPiZeroAOD::IsRealProton(AliAODTrack* track, TClonesArray* fAODMCTrackArray)
{//Check if proton comes frome a sigma+
	Int_t mcID = track->GetLabel();
	AliAODMCParticle* mcTrack = static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(TMath::Abs(mcID)));
	if(!mcTrack) return kFALSE;
	if (mcTrack->GetMother() < 0) return kFALSE;
	AliAODMCParticle *TrackMother = static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(mcTrack->GetMother())); //mother MC particle object
	if(!TrackMother) return kFALSE;
	Int_t codeTrack = mcTrack->GetPdgCode(); 
	Int_t codeMother = TrackMother->GetPdgCode(); 
	if( (codeTrack==2212) && (codeMother==3222) ){
		return kTRUE;
	} else {
		return kFALSE;
	}
}
//________________________________________________________________________
Bool_t AliAnalysisTaskSigmaPlToProtonPiZeroAOD::IsRealPhoton(AliAODConversionPhoton *PhotonCandidate, TClonesArray* fAODMCTrackArray)
{ //checks if a reconstructed photon from sigma plus
	AliAODMCParticle *Photon = NULL;
	if (PhotonCandidate->GetNCaloPhotonMCLabels() > 0){
		// Photon = PhotonCandidate->GetMCParticle(fMCEvent);
		if (PhotonCandidate->GetCaloPhotonMCLabel(0) < 0) return kFALSE;
		Photon = static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(PhotonCandidate->GetCaloPhotonMCLabel(0)));
		if (Photon) {
			if(Photon->GetPdgCode() == 22){
				if (Photon->GetMother() < 0) return kFALSE;
				AliAODMCParticle* motherPart2 = static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(Photon->GetMother()));
				if (motherPart2->GetMother() < 0) return kFALSE;
				AliAODMCParticle* grandmotherPart2 = static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(motherPart2->GetMother()));
				if(motherPart2->GetPdgCode() == 111 && grandmotherPart2->GetPdgCode() == 3222){
					return kTRUE;
				} else {
					return kFALSE;
				}
			} else if(Photon->GetPdgCode() == 11 || Photon->GetPdgCode() == -11){
				if (Photon->GetMother() < 0) return kFALSE;
				AliAODMCParticle* motherPart2 = static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(Photon->GetMother()));
				if (motherPart2->GetMother() < 0) return kFALSE;
				AliAODMCParticle* grandmotherPart2 = static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(motherPart2->GetMother()));
				if (grandmotherPart2->GetMother() < 0) return kFALSE;
				AliAODMCParticle* grandgrandmotherPart2 = static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(grandmotherPart2->GetMother()));
				if(motherPart2->GetPdgCode() == 22 && grandmotherPart2->GetPdgCode() == 111  && grandgrandmotherPart2->GetPdgCode() == 3222){
					return kTRUE;
				} else {
					return kFALSE;
				}
			} else {
				return kFALSE;
			}
		}
	}
	return kFALSE;

}

//________________________________________________________________________
void AliAnalysisTaskSigmaPlToProtonPiZeroAOD::CalculateBackgroundSwapp(vector < AliVCluster* > photon, vector < AliAODTrack* > proton, Double_t vpos[3], Int_t iCut){

    Double_t rotationAngle = TMath::Pi()/8.0; //0.78539816339; // rotaion angle 45°
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

    std::vector<std::array<Double_t, 2>> vSwappingInvMassPT;
    vSwappingInvMassPT.clear();
    vSwappingInvMassPT.resize(0);
    Double_t tempMultWeightSwapping = 1.; // weight taking multiplicity of event into account

    // curcial requierment is that the event has at least 3 cluster candidates
    if(photon.size() > 2 || proton.size() > 0 ){
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

			        for(int iSwapp = 0; iSwapp < 1; ++iSwapp){
			         // for(int iSwapp = 0; iSwapp < ((AliConversionMesonCuts*)fMesonCutArray->At(iCut))->GetNumberOfSwappsForBg(); ++iSwapp){
						photonCandidate->GetMomentum(lvRotationPhoton,vpos);
						photon1Candidate->GetMomentum(lvRotationPhoton1,vpos);
			         
			            vRotationPion = (lvRotationPhoton + lvRotationPhoton1).Vect();
			            lvRotationPhoton.Rotate(rotationAngle, vRotationPion);
			            lvRotationPhoton1.Rotate(rotationAngle, vRotationPion);

						if(lvRotationPhoton1.Phi() < 0.){
		            		if( (fabs(lvRotationPhoton1.Eta()) > 0.13) || ((lvRotationPhoton1.Phi()+2.*TMath::Pi()) < 250.f*(TMath::Pi()/180.f)) || ((lvRotationPhoton1.Phi()+2.*TMath::Pi()) > 320.f*(TMath::Pi()/180.f))){
		            		continue;
		            		} 
		            	}
		            	else{
			            	if( (fabs(lvRotationPhoton1.Eta()) > 0.13) || (lvRotationPhoton1.Phi() < 250.f*(TMath::Pi()/180.f)) || (lvRotationPhoton1.Phi() > 320.f*(TMath::Pi()/180.f))){
			            		continue;
			            	} 
		            	}

		            	if(lvRotationPhoton.Phi() < 0.){
		            		if( (fabs(lvRotationPhoton.Eta()) > 0.13) || ((lvRotationPhoton.Phi()+2.*TMath::Pi()) < 250.f*(TMath::Pi()/180.f)) || ((lvRotationPhoton.Phi()+2.*TMath::Pi()) > 320.f*(TMath::Pi()/180.f))){
		            		continue;
		            		} 
		            	}
		            	else{
			            	if( (fabs(lvRotationPhoton.Eta()) > 0.13) || (lvRotationPhoton.Phi() < 250.f*(TMath::Pi()/180.f)) || (lvRotationPhoton.Phi() > 320.f*(TMath::Pi()/180.f))){
			            		continue;
			            	} 
		            	}

		            	cellIDRotatedPhoton = ((AliCaloPhotonCuts*)fClusterCutArray->At(iCut))->GetCaloCellIdFromEtaPhi(lvRotationPhoton.Eta(), static_cast<double>((lvRotationPhoton.Phi()<0) ? lvRotationPhoton.Phi() + TMath::Pi()*2. : lvRotationPhoton.Phi()));
	           			cellIDRotatedPhoton1 = ((AliCaloPhotonCuts*)fClusterCutArray->At(iCut))->GetCaloCellIdFromEtaPhi(lvRotationPhoton1.Eta(), static_cast<double>((lvRotationPhoton1.Phi()<0) ? lvRotationPhoton1.Phi() + TMath::Pi()*2. : lvRotationPhoton1.Phi()));
	           			std::unique_ptr<AliAODConversionPhoton> currentEventRotation (new AliAODConversionPhoton(&lvRotationPhoton));
	           			std::unique_ptr<AliAODConversionPhoton> currentEventRotation1 (new AliAODConversionPhoton(&lvRotationPhoton1));

		            	for(unsigned int iPionBG=0;iPionBG<photon.size();iPionBG++){
		            		if(iPionBG == iCurrent2 || iPionBG == iCurrent1) continue;
			          		 AliVCluster* photonBGCandidate = photon[iPionBG];
							photonBGCandidate->GetMomentum(lvRotationBGPhoton,vpos);
			          		//First Swapped Gamma
		            		lvRotationBGPion = (lvRotationPhoton + lvRotationBGPhoton);
		            		if(lvRotationBGPion.M() < 0.118 || lvRotationBGPion.M() > 0.148 ) continue;
			          		//Second Swapped Gamma
		                  	lvRotationBGPion1 = (lvRotationPhoton1 + lvRotationBGPhoton);
		            		if(lvRotationBGPion1.M() < 0.118 || lvRotationBGPion1.M() > 0.148 ) continue;
		            		TLorentzVector rekombinatedPi0BG;
							rekombinatedPi0BG.SetPtEtaPhiM(lvRotationBGPion.Pt(), lvRotationBGPion.Eta(), lvRotationBGPion.Phi(), 0.135);
							TLorentzVector rekombinatedPi0BG1;
							rekombinatedPi0BG1.SetPtEtaPhiM(lvRotationBGPion1.Pt(), lvRotationBGPion1.Eta(), lvRotationBGPion1.Phi(), 0.135);
	           				std::unique_ptr<AliAODConversionPhoton> currentEventGoodBGPhoton (new AliAODConversionPhoton(&lvRotationBGPhoton));
							std::unique_ptr<AliAODConversionMother> backgroundCandidate(new AliAODConversionMother(currentEventRotation.get(), currentEventGoodBGPhoton.get()));
							std::unique_ptr<AliAODConversionMother> backgroundCandidate1(new AliAODConversionMother(currentEventRotation1.get(), currentEventGoodBGPhoton.get()));
	              
	 						if(!(((AliCaloPhotonCuts*)fClusterCutArray->At(iCut))->CheckDistanceToBadChannelSwapping(cellIDRotatedPhoton, lvRotationPhoton.Phi(), fInputEvent)) && lvRotationPhoton.E() > ((AliCaloPhotonCuts*)fClusterCutArray->At(iCut))->GetMinClusterEnergy())
	              			{
	               				if(((AliConversionMesonCuts*) fMesonCutArray->At(iCut))->MesonIsSelected(backgroundCandidate.get(),kFALSE,((AliConvEventCuts*)fEventCutArray->At(iCut))->GetEtaShift(), cellIDRotatedPhoton, currentEventGoodBGPhoton.get()->GetLeadingCellID()))
	               				{
				          		//First Pion
						            lvRotationBGSigma = (rekombinatedPi0BG + lvRotationBGProton);
				            		if(GetPodAlpha(lvRotationBGSigma, lvRotationBGProton, rekombinatedPi0BG) > 1.02 || GetPodAlpha(lvRotationBGSigma, lvRotationBGProton, rekombinatedPi0BG) < 0.2) continue;
									if(GetQT(lvRotationBGSigma, lvRotationBGProton) > fFitPodolanskiUpperCut[iCut]->Eval(GetPodAlpha(lvRotationBGSigma, lvRotationBGProton, rekombinatedPi0BG))) continue;
									if(rekombinatedPi0BG.Angle(lvRotationBGProton.Vect()) > 0.35 || rekombinatedPi0BG.Angle(lvRotationBGProton.Vect()) < 0.05) continue;
				                 	vSwappingInvMassPT.push_back({lvRotationBGSigma.M(),lvRotationBGSigma.Pt()});
		                   		}
	              			}
	              			if(!(((AliCaloPhotonCuts*)fClusterCutArray->At(iCut))->CheckDistanceToBadChannelSwapping(cellIDRotatedPhoton1, lvRotationPhoton1.Phi(), fInputEvent)) && lvRotationPhoton1.E() > ((AliCaloPhotonCuts*)fClusterCutArray->At(iCut))->GetMinClusterEnergy())
	              			{
	               				if(((AliConversionMesonCuts*) fMesonCutArray->At(iCut))->MesonIsSelected(backgroundCandidate1.get(),kFALSE,((AliConvEventCuts*)fEventCutArray->At(iCut))->GetEtaShift(), cellIDRotatedPhoton1, currentEventGoodBGPhoton.get()->GetLeadingCellID()))
	               				{
				          		//Second Pion
				            		lvRotationBGSigma1 = (rekombinatedPi0BG1 + lvRotationBGProton);
				            		if(GetPodAlpha(lvRotationBGSigma1, lvRotationBGProton, rekombinatedPi0BG1) > 1.02 || GetPodAlpha(lvRotationBGSigma1, lvRotationBGProton, rekombinatedPi0BG1) < 0.2) continue;
									if(GetQT(lvRotationBGSigma1, lvRotationBGProton) > fFitPodolanskiUpperCut[iCut]->Eval(GetPodAlpha(lvRotationBGSigma1, lvRotationBGProton, rekombinatedPi0BG1))) continue;
									if(rekombinatedPi0BG1.Angle(lvRotationBGProton.Vect()) > 0.35 || rekombinatedPi0BG1.Angle(lvRotationBGProton.Vect()) < 0.05) continue;
				                 	vSwappingInvMassPT.push_back({lvRotationBGSigma1.M(),lvRotationBGSigma1.Pt()});
						        	
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
        	fHistoMotherBackInvMassPt[iCut]->Fill(vSwappingInvMassPT.at(i)[0], vSwappingInvMassPT.at(i)[1], tempMultWeightSwapping);
        }
    }
}

//_____________________________________________________________________________
void AliAnalysisTaskSigmaPlToProtonPiZeroAOD::Terminate(Option_t *)
{
	// terminate
	// called at the END of the analysis (when all events are processed)
}
//_____________________________________________________________________________
