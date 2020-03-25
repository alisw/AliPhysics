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
#include "AliESDCaloCluster.h"
#include "AliPIDResponse.h"
#include "AliAnalysisTaskSigmaPlToProtonPiZero.h"
#include "AliESDVertex.h"
#include "AliVEvent.h"
#include "AliESDkink.h"
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
#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliMCParticle.h"
#include "AliCentrality.h"
#include "AliESDVZERO.h"
#include "AliESDpid.h"
#include "AliVParticle.h"
#include "AliESDtrack.h"
#include "AliESDtrackCuts.h"
#include "AliKFVertex.h"
#include "AliV0ReaderV1.h"
#include "AliGenCocktailEventHeader.h"
#include "AliAODMCParticle.h"
#include "AliAODMCHeader.h"
#include "AliEventplane.h"
#include "AliAnalysisTaskEMCALClusterizeFast.h"
#include "AliAODEvent.h"
#include "AliAODInputHandler.h"
#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliInputEventHandler.h"
#include "AliCaloTrackMatcher.h"
#include <vector>
#include <map>
#include <fstream>
#include "TVector3.h"
#include "TLorentzVector.h"
#include "AliAODConversionPhoton.h"
#include "AliAODConversionMother.h"



class AliAnalysisTaskSigmaPlToProtonPiZero;    // your analysis class


using namespace std;            // std namespace: so you can do things like 'cout'

ClassImp(AliAnalysisTaskSigmaPlToProtonPiZero) // classimp: necessary for root

double pi = TMath::Pi();


AliAnalysisTaskSigmaPlToProtonPiZero::AliAnalysisTaskSigmaPlToProtonPiZero() : AliAnalysisTaskSE(),
	fEvent(0),
	fESD(0),
  fMCEvent(0),
  fMCStack(0),
  fPIDResponse(0),
  fOutputList(0),
  fESDList(NULL),
  fHistSigmaPlus(0),
  fHistSigmaPlusMC(0),
  fHistSigmaPlusMCGen(0),
  fHistReconstructedMassPi0(0),
  fHistReconstructedMassPi0MC(0),
  fHistPodolanski(0),
  fHistAllTracksPt(0),
  fHistProtonPt(0),
  fHistThetaPhi(0),
  fHistThetaPhiProton(0),
  fHistClusterE(0),
  fHistClusterEWOCuts(0),
  fHistNClusWoCuts(0),
  fHistNClusWCuts(0),
  fHistDEDx(0),
  fHistTPCSignal(0),
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
  fTrackMatcherRunningMode(0)


							// fESD(0), fOutputList(0), fHistAllTracksPt(0), fHistProtonPt(0), fHistVertexzPos(0), fHistVertexzPosTrue(0), fHistThetaPhi(0), fHistThetaPhiProton(0)
{
	// default constructor, don't allocate memory here!
	// this is used by root for IO purposes, it needs to remain empty
}
//_____________________________________________________________________________
AliAnalysisTaskSigmaPlToProtonPiZero::AliAnalysisTaskSigmaPlToProtonPiZero(const char* name) : AliAnalysisTaskSE(name),
	fEvent(0),
	fESD(0),
  fMCEvent(0),
  fMCStack(0),
  fPIDResponse(0),
  fOutputList(0),
  fESDList(NULL),
  fHistSigmaPlus(0),
  fHistSigmaPlusMC(0),
  fHistSigmaPlusMCGen(0),
  fHistReconstructedMassPi0(0),
  fHistReconstructedMassPi0MC(0),
  fHistPodolanski(0),
  fHistAllTracksPt(0),
  fHistProtonPt(0),
  fHistThetaPhi(0),
  fHistThetaPhiProton(0),
  fHistClusterE(0),
  fHistClusterEWOCuts(0),
  fHistNClusWoCuts(0),
  fHistNClusWCuts(0),
  fHistDEDx(0),
  fHistTPCSignal(0),
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
  fMaxTimingCluster(0),                                    // corresponding ranges, max
  fTrackMatcherRunningMode(0)


{
	DefineOutput(1, TList::Class());
	DefineOutput(2, TTree::Class());
}
//_____________________________________________________________________________
AliAnalysisTaskSigmaPlToProtonPiZero::~AliAnalysisTaskSigmaPlToProtonPiZero()
{
	// destructor
	// if(fOutputList) {
	// 	delete fOutputList;
	// }
}
//_____________________________________________________________________________
void AliAnalysisTaskSigmaPlToProtonPiZero::UserCreateOutputObjects()
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
	// example of a histogram
	const Int_t Anzahl_Bins[5] = {20, 30, 15, 10, 10};
	const Double_t Lower_bin_edges[5] = {1.1 ,0. ,0. ,0. ,0.};
	const Double_t upper_bin_edges[5] = {1.3, 30., 15., 250., pi};

	fESDList            = new TList*[fnCuts];

	fHistSigmaPlus = new THnD*[fnCuts];
	fHistReconstructedMassPi0 = new TH2F*[fnCuts];
	fHistPodolanski = new TH2F*[fnCuts];
	fHistAllTracksPt = new TH1F*[fnCuts];
	fHistProtonPt = new TH1F*[fnCuts];
	fHistThetaPhi = new TH2F*[fnCuts];
	fHistThetaPhiProton = new TH2F*[fnCuts];
	fHistClusterE = new TH1F*[fnCuts];
	fHistClusterEWOCuts = new TH1F*[fnCuts];
	fHistNClusWoCuts = new TH1F*[fnCuts];
	fHistNClusWCuts = new TH1F*[fnCuts];
	fHistDEDx = new TH2F*[fnCuts];
	fHistTPCSignal = new TH2F*[fnCuts];
	fHistNEvents = new TH1F*[fnCuts];
	fFitPi0MassDataLowPt = new TF1*[fnCuts];
	fFitPi0MassDataHighPt = new TF1*[fnCuts];
	fFitWidthData = new TF1*[fnCuts];

	if(fIsMC > 0){
		fHistSigmaPlusMC = new THnD*[fnCuts];
		fHistReconstructedMassPi0MC = new TH2F*[fnCuts];
		fFitPi0MassMCLowPt = new TF1*[fnCuts];
		fFitPi0MassMCHighPt = new TF1*[fnCuts];
		fFitWidthMC = new TF1*[fnCuts];

		fHistSigmaPlusMCGen = new TH1F("fHistSigmaPlusMCGen",";#it{p}_{T} (GeV/#it{c});Yield", 30, 0., 30.);
		fHistSigmaPlusMCGen->Sumw2();
	}

  	for(Int_t iCut = 0; iCut<fnCuts;iCut++){

  		TString cutstringEvent    = ((AliConvEventCuts*)fEventCutArray->At(iCut))->GetCutNumber();
    	TString cutstringCalo     = ((AliCaloPhotonCuts*)fClusterCutArray->At(iCut))->GetCutNumber();
    	TString cutstringMeson    = "NoMesonCut";

  		fESDList[iCut]          = new TList();
    	fESDList[iCut]->SetOwner(kTRUE);
    	fESDList[iCut]->SetName(Form("%s_%s_%s ESD histograms", cutstringEvent.Data(), cutstringCalo.Data(), cutstringMeson.Data()));
    	fOutputList->Add(fESDList[iCut]);

	  	fHistSigmaPlus[iCut] = new THnD("fHistSigmaPlus", "", 5, Anzahl_Bins, Lower_bin_edges, upper_bin_edges );
		fHistSigmaPlus[iCut]->Sumw2();
		fESDList[iCut]->Add(fHistSigmaPlus[iCut]);
		fHistReconstructedMassPi0[iCut] = new TH2F("fHistReconstructedMassPi0",";#it{m}_{inv} (GeV/#it{c^{2}});#it{p}_{T} (GeV/#it{c})", 60, 0., 0.3, 30, 0., 15.);
		fHistReconstructedMassPi0[iCut]->Sumw2();
		fESDList[iCut]->Add(fHistReconstructedMassPi0[iCut]);
		fHistPodolanski[iCut] = new TH2F("fHistPodolanski","", 20, 0., 1., 20, 0., 1.);
		fHistPodolanski[iCut]->Sumw2();
		fESDList[iCut]->Add(fHistPodolanski[iCut]);
		fHistAllTracksPt[iCut] = new TH1F("fHistAllTracksPt", "fHistAllTracksPt;#it{p}_{T} (GeV/#it{c});Yield", 100, 0, 10);
		fHistAllTracksPt[iCut]->Sumw2();
		fESDList[iCut]->Add(fHistAllTracksPt[iCut]);
		fHistProtonPt[iCut] = new TH1F("fHistProtonPt", "fHistProtonPt;#it{p}_{T,Proton} (GeV/#it{c});Yield", 100, 0, 10);
		fHistProtonPt[iCut]->Sumw2();
		fESDList[iCut]->Add(fHistProtonPt[iCut]);
		fHistThetaPhi[iCut] = new TH2F("fHistThetaPhi", "fHistThetaPhi;#eta ; #phi", 50, 0., pi ,100, 0., 2*pi);
		fHistThetaPhi[iCut]->Sumw2();
		fESDList[iCut]->Add(fHistThetaPhi[iCut]);
		fHistThetaPhiProton[iCut] = new TH2F("fHistThetaPhiProton", "fHistThetaPhiProton;#eta_{p} ; #phi_{p}", 50, 0., pi ,100, 0., 2*pi);
		fHistThetaPhiProton[iCut]->Sumw2();
		fESDList[iCut]->Add(fHistThetaPhiProton[iCut]);
		fHistClusterE[iCut] = new TH1F("fHistClusterE", "fHistClusterE;#it{E} (GeV);Yield", 100, 0, 30);
		fHistClusterE[iCut]->Sumw2();
		fESDList[iCut]->Add(fHistClusterE[iCut]);
		fHistClusterEWOCuts[iCut] = new TH1F("fHistClusterEWOCuts", "fHistClusterEWOCuts;#it{E}_{Test} (GeV);Yield", 100, 0, 30);
		fHistClusterEWOCuts[iCut]->Sumw2();
		fESDList[iCut]->Add(fHistClusterEWOCuts[iCut]);
		fHistNClusWoCuts[iCut] = new TH1F("fHistNClusWoCuts", "fHistNClusWoCuts;#it{N}_{Cluster, wo. Cuts};Yield", 30, -0.5, 29.5);
		fHistNClusWoCuts[iCut]->Sumw2();
		fESDList[iCut]->Add(fHistNClusWoCuts[iCut]);
		fHistNClusWCuts[iCut] = new TH1F("fHistNClusWCuts", "fHistNClusWCuts;#it{N}_{Cluster,w. Cuts};Yield", 30, -0.5, 29.5);
		fHistNClusWCuts[iCut]->Sumw2();
		fESDList[iCut]->Add(fHistNClusWCuts[iCut]);
		fHistDEDx[iCut] = new TH2F("fHistDEDx", "fHistDEDx;#it{p};#d it{E}/d it{x}", 1000,0.01,100.,1000,1.,2000.);
		fHistDEDx[iCut]->Sumw2();
		fESDList[iCut]->Add(fHistDEDx[iCut]);
		fHistTPCSignal[iCut] = new TH2F("fHistTPCSignal", "fHistTPCSignal;#it{p}_{T};#sigma_{TPC}", 100, 0., 10., 10, -5., 5.);
		fHistTPCSignal[iCut]->Sumw2();
		fESDList[iCut]->Add(fHistTPCSignal[iCut]);

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

		fESDList[iCut]->Add(fHistNEvents[iCut]);          // don't forget to add it to the list! the list will be written to file, so if you want

		if(!((AliConvEventCuts*)fEventCutArray->At(iCut))) continue;
    	if(((AliConvEventCuts*)fEventCutArray->At(iCut))->GetCutHistograms()){
      	fESDList[iCut]->Add(((AliConvEventCuts*)fEventCutArray->At(iCut))->GetCutHistograms());
    	}
    	if(!((AliCaloPhotonCuts*)fClusterCutArray->At(iCut))) continue;
    	if(((AliCaloPhotonCuts*)fClusterCutArray->At(iCut))->GetCutHistograms()){
        fESDList[iCut]->Add(((AliCaloPhotonCuts*)fClusterCutArray->At(iCut))->GetCutHistograms());
    	}
      	if(!((AliConversionMesonCuts*)fMesonCutArray->At(iCut))) continue;
     	if(((AliConversionMesonCuts*)fMesonCutArray->At(iCut))->GetCutHistograms()){
        fESDList[iCut]->Add(((AliConversionMesonCuts*)fMesonCutArray->At(iCut))->GetCutHistograms());
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
		fESDList[iCut]->Add(fFitPi0MassDataHighPt[iCut]);
		fESDList[iCut]->Add(fFitWidthData[iCut]);
		fESDList[iCut]->Add(fFitPi0MassDataLowPt[iCut]);

		if(fIsMC > 0){

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
			fESDList[iCut]->Add(fFitPi0MassMCHighPt[iCut]);
			fESDList[iCut]->Add(fFitPi0MassMCLowPt[iCut]);
			fESDList[iCut]->Add(fFitWidthMC[iCut]);

			fHistSigmaPlusMC[iCut] = new THnD("fHistSigmaPlusMC", "", 5, Anzahl_Bins, Lower_bin_edges, upper_bin_edges );
			fHistSigmaPlusMC[iCut]->Sumw2();
			fESDList[iCut]->Add(fHistSigmaPlusMC[iCut]);

			fESDList[iCut]->Add(fHistSigmaPlusMCGen);
			fHistReconstructedMassPi0MC[iCut] = new TH2F("fHistReconstructedMassPi0MC",";#it{m}_{inv} (GeV/#it{c^{2}});#it{p}_{T} (GeV/#it{c})", 60, 0., 0.3, 30, 0., 15.);
			fHistReconstructedMassPi0MC[iCut]->Sumw2();
			fESDList[iCut]->Add(fHistReconstructedMassPi0MC[iCut]);
		}
	}
	if(fV0Reader)
	  if((AliConvEventCuts*)fV0Reader->GetEventCuts())
		if(((AliConvEventCuts*)fV0Reader->GetEventCuts())->GetCutHistograms())
		  fOutputList->Add(((AliConvEventCuts*)fV0Reader->GetEventCuts())->GetCutHistograms());
	PostData(1, fOutputList);
}
//_____________________________________________________________________________
void AliAnalysisTaskSigmaPlToProtonPiZero::UserExec(Option_t *)
{

	fEvent = InputEvent();
	if (!fEvent) {
		Printf("ERROR: Could not retrieve event");
		return;
	}
	fESD = dynamic_cast<AliESDEvent*>(fEvent);
	if (!fESD) {
		Printf("ERROR: Could not retrieve fESD");
		return;
	}
	//Eventqualitäts check
	if(fIsMC> 0){
		fMCEvent = MCEvent();
		if(!fMCEvent)  {
			Printf("ERROR: Could not retrieve fMCEvent");
			return;
		}
		fMCStack = fMCEvent->Stack(); //stack of MC events
		if(!fMCStack)  {
			Printf("ERROR: Could not retrieve fMCStack");
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
		for(Int_t i = 0; i < fMCEvent->GetNumberOfTracks(); i++) {
			TParticle* sigma = (TParticle *)fMCEvent->Particle(i);
			if (!sigma) continue;
			if(sigma->GetPdgCode() == 3222){
				Double_t sigmaPt = sigma->Pt();
				if(fHistSigmaPlusMCGen) fHistSigmaPlusMCGen->Fill(sigmaPt);
			}
		}
	}


	//Find EventVertex
	const AliESDVertex *vertex=fESD->GetPrimaryVertex();    // 22/8
	if(!vertex) {
		Printf("ERROR: Could not retrieve vertex");
		return;
	}
	//Primär Vertex Bestimmen
	Double_t vpos[3];
	vertex->GetXYZ(vpos);
	//Find Protons and secundary Vertex
	Int_t iTracks=0;
	if(!fESD->GetNumberOfTracks()) {
		Printf("ERROR: Could not retrieve fESD->GetNumberOfTracks()");
		return;
	}
	iTracks=fESD->GetNumberOfTracks();  // see how many tracks there are in the event


	if(!fPIDResponse) {
		AliAnalysisManager *man = AliAnalysisManager::GetAnalysisManager();
		if(man){
			AliInputEventHandler* inputHandler = (AliInputEventHandler*)(man->GetInputEventHandler());
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
		vector < AliESDtrack* > proton;
		vector < AliESDtrack* > tracks;
		vector < AliVCluster* > photon;

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

		AliESDtrack* track;
		for(Int_t i = 0; i < iTracks; i++) {

			AliESDtrack* tracktmp;
			tracktmp = static_cast<AliESDtrack*>(fESD->GetTrack(i));
			if(!tracktmp){
				Printf("ERROR: Could not retrieve track %i", i);
				continue;
			}
			track = new AliESDtrack(*tracktmp);
			//Kinks of doughter
			if(track->GetKinkIndex(0) > 0){
				double protonSignal = (fPIDResponse->NumberOfSigmasTPC(track, AliPID::kProton)); //numbers of sigmas of TPC signal
				// Select Proton Candidates
				if (std::abs(protonSignal) < 3 ) {
					if(fHistProtonPt[iCut]) fHistProtonPt[iCut]->Fill(track->Pt());
					if(fHistThetaPhiProton[iCut]) fHistThetaPhiProton[iCut]->Fill(track->Theta(), track->Phi());
					if(fHistTPCSignal[iCut]) fHistTPCSignal[iCut]->Fill(track->Pt(), protonSignal);
					proton.push_back(track);
				} else {
					tracks.push_back(track);
				}
			} else {
				tracks.push_back(track);
			}
			if(fHistAllTracksPt[iCut]) fHistAllTracksPt[iCut]->Fill(track->Pt());
			if(fHistThetaPhi[iCut]) fHistThetaPhi[iCut]->Fill(track->Theta(), track->Phi());
			//dE/dx Plot
			if (!track->GetInnerParam()) continue;
			Double_t ptot = track->GetInnerParam()->P();
			Double_t tpcSignal = track->GetTPCsignal();
			if(fHistDEDx[iCut]) fHistDEDx[iCut]->Fill(ptot, tpcSignal);
		}

		if(proton.size() > 0){
			//Find gammas
			Int_t nclus                     = 0;
			Int_t nClusWCuts                = 0;
			nclus = fESD->GetNumberOfCaloClusters();
			if(fHistNClusWoCuts[iCut]) fHistNClusWoCuts[iCut]->Fill(nclus);
			if(nclus == 0)  continue;
			AliVCluster* clus = NULL;
			for(Int_t i=0; i < nclus; ++i){
				clus = new AliESDCaloCluster(*(AliESDCaloCluster*)fEvent->GetCaloCluster(i));
				if(!clus) {
					Printf("ERROR: Could not find clus %i",i);
					continue;
				}
				if(fHistClusterEWOCuts[iCut]) fHistClusterEWOCuts[iCut]-> Fill(clus->E());
				if((((AliCaloPhotonCuts*)fClusterCutArray->At(iCut))->ClusterIsSelected(clus,fESD,fMCEvent,fIsMC, 1.,i)) == kTRUE){
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
			TLorentzVector rekombinatedPi0;
			TLorentzVector photonCandidate1;
			TLorentzVector photonCandidate2;
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
					AliESDtrack* protonCandidate = proton[iProton];
					protonVektor.SetPtEtaPhiM(protonCandidate->Pt(),protonCandidate->Eta(),protonCandidate->Phi(), 0.938272);
					AliESDkink *kink=fESD->GetKink(TMath::Abs(protonCandidate->GetKinkIndex(0))-1);
					if(!kink){
						Printf("ERROR: Could not find kink[%i][%i]",iCut,iProton);
						continue;
					}
					const TVector3 vposKink(kink->GetPosition());
					Double_t trackKink[3] = {vposKink[0], vposKink[1], vposKink[2]};
					Double_t kinkDistanceToVertex = sqrt(vposKink[0] * vposKink[0] + vposKink[1] * vposKink[1] + vposKink[2] * vposKink[2]);
					if(fIsMC > 0){
						trueSigmaProton = IsRealProtonKink(kink, fMCStack);
					}
					for(unsigned int iPhoton1 = 0; iPhoton1 < photon.size(); ++iPhoton1) {
						trueSigmaPhoton1 = kFALSE;
						if(!photon[iPhoton1]) {
							Printf("ERROR: Could not find photon[%i][%i]",iCut,iPhoton1);
							continue;
						}
						AliVCluster* gamma1 = photon[iPhoton1];
						TLorentzVector clusterVector1;
					    gamma1->GetMomentum(clusterVector1,trackKink);
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
							trueSigmaPhoton1 = IsRealPhoton(PhotonCandidate1, fMCEvent);
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
							    gamma2->GetMomentum(clusterVector2,trackKink);
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
									trueSigmaPhoton2 = IsRealPhoton(PhotonCandidate2, fMCEvent);
								}
							    AliAODConversionMother pi0cand = AliAODConversionMother(&PhotonCandidate1,&PhotonCandidate2);
								pi0cand.SetLabels(iPhoton1,iPhoton2);
								if((((AliConversionMesonCuts*)fMesonCutArray->At(iCut))->MesonIsSelected(&pi0cand,kTRUE,((AliConvEventCuts*)fEventCutArray->At(iCut))->GetEtaShift(),PhotonCandidate1.GetLeadingCellID(),PhotonCandidate2.GetLeadingCellID()))){
									if(fHistReconstructedMassPi0[iCut]) fHistReconstructedMassPi0[iCut]->Fill(pi0cand.M(), pi0cand.Pt());
									if(trueSigmaPhoton1 == kTRUE && trueSigmaPhoton2 ==kTRUE && fIsMC > 0){
										if(fHistReconstructedMassPi0MC[iCut]) fHistReconstructedMassPi0MC[iCut]->Fill(pi0cand.M(), pi0cand.Pt());
									}
									mesonSelected = kFALSE;
									if(fIsMC > 0 ) {
										mesonSelected = IsPi0SelectedMC(&pi0cand, fFitPi0MassMCLowPt[iCut], fFitPi0MassMCHighPt[iCut], fFitWidthData[iCut]);
									}
									if(fIsMC == 0 ) {
										mesonSelected = IsPi0Selected(&pi0cand, fFitPi0MassDataLowPt[iCut], fFitPi0MassDataHighPt[iCut], fFitWidthData[iCut]);
									}
									if(mesonSelected == kTRUE){
										rekombinatedPi0.SetPtEtaPhiM(pi0cand.Pt(), pi0cand.Eta(), pi0cand.Phi(), 0.135);
										sigmaVektor = protonVektor + rekombinatedPi0;
										if(fHistPodolanski[iCut]) fHistPodolanski[iCut]->Fill(GetPodAlpha(sigmaVektor, protonVektor, rekombinatedPi0),GetQT(sigmaVektor, rekombinatedPi0));
										if(fHistPodolanski[iCut]) fHistPodolanski[iCut]->Fill(GetPodAlpha(sigmaVektor, protonVektor, rekombinatedPi0),GetQT(sigmaVektor, protonVektor));
										Double_t farrSigma[5] = {sigmaVektor.M(), sigmaVektor.Pt(), protonVektor.Pt(), kinkDistanceToVertex, rekombinatedPi0.Angle(protonVektor.Vect())};
										if(fHistSigmaPlus[iCut]) fHistSigmaPlus[iCut]-> Fill(farrSigma);
										if(trueSigmaProton == kTRUE && trueSigmaPhoton1 == kTRUE && trueSigmaPhoton2 ==kTRUE && fIsMC > 0){
											if(fHistSigmaPlusMC [iCut]) fHistSigmaPlusMC [iCut]-> Fill(farrSigma);
										}
									}
								}
							}
						}
					}
				}
			}
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
	}


	// continue until all the tracks are processed
	PostData(1, fOutputList);                           // stream the results the analysis of this event to
	// the output manager which will take care of writing
	// it to a file
}
//_____________________________________________________________________________
Bool_t AliAnalysisTaskSigmaPlToProtonPiZero::IsPi0Selected(AliAODConversionMother* pi0CandTmp, TF1* fitPi0MassDataLowPt, TF1* fitPi0MassDataHighPt, TF1* fitWidthData)
{
	Double_t minPt = 0.3;
	// Double_t maxPtMassDataFirst = 1.2;
	Double_t minPtMassDataFirst = 1.1;
	Double_t maxPt = 25.;

	Double_t minPi0Mass = 0;
	Double_t maxPi0Mass = 0;
	Double_t width = 0;
	Double_t meanMass = 0;

	if ( (pi0CandTmp->Pt()) > minPt && (pi0CandTmp->Pt()) < maxPt) {
		width = fitWidthData->Eval(pi0CandTmp->Pt());
		if((pi0CandTmp->Pt()) <= minPtMassDataFirst ){
			meanMass = fitPi0MassDataLowPt->Eval(pi0CandTmp->Pt());
		} else if((pi0CandTmp->Pt()) > minPtMassDataFirst){
			meanMass = fitPi0MassDataHighPt->Eval(pi0CandTmp->Pt());
		} else {
			return kFALSE;
		}
		minPi0Mass = meanMass - 3*width;
		maxPi0Mass = meanMass + 3*width;
		if(pi0CandTmp->M() < maxPi0Mass && pi0CandTmp->M() > minPi0Mass){
			return kTRUE;
		}
	}
	return kFALSE;
}
//_____________________________________________________________________________
Bool_t AliAnalysisTaskSigmaPlToProtonPiZero::IsPi0SelectedMC(AliAODConversionMother *pi0CandTmp, TF1* fitPi0MassMCLowPt, TF1* fitPi0MassMCHighPt, TF1* fitWidthMC)
{
	Double_t minPt = 0.3;
	// Double_t maxPtMassMCFirst = 1.05;
	Double_t minPtMassMCSecond = 1.;
	Double_t maxPt = 25.;

	Double_t minPi0Mass = 0;
	Double_t maxPi0Mass = 0;
	Double_t width = 0;
	Double_t meanMass = 0;

	if ( (pi0CandTmp->Pt()) > minPt && (pi0CandTmp->Pt()) < maxPt) {
		width = fitWidthMC->Eval(pi0CandTmp->Pt());
		if((pi0CandTmp->Pt()) <= minPtMassMCSecond ){
			meanMass = fitPi0MassMCLowPt->Eval(pi0CandTmp->Pt());
		} else if((pi0CandTmp->Pt()) > minPtMassMCSecond ){
			meanMass = fitPi0MassMCHighPt->Eval(pi0CandTmp->Pt());
		} else {
			return kFALSE;
		}
		minPi0Mass = meanMass - 3*width;
		maxPi0Mass = meanMass + 3*width;
		if(pi0CandTmp->M() < maxPi0Mass && pi0CandTmp->M() > minPi0Mass){
			return kTRUE;
		}
	}

	return kFALSE;
}
//_____________________________________________________________________________
Double_t AliAnalysisTaskSigmaPlToProtonPiZero::GetPodAlpha(TLorentzVector sigmaVektor, TLorentzVector protonVektor, TLorentzVector rekombinatedPi0)
{
	return (cos(sigmaVektor.Angle(protonVektor.Vect())) * protonVektor.P() - cos(sigmaVektor.Angle(rekombinatedPi0.Vect())) * rekombinatedPi0.P()) / (cos(sigmaVektor.Angle(protonVektor.Vect())) * protonVektor.P() + cos(sigmaVektor.Angle(rekombinatedPi0.Vect())) * rekombinatedPi0.P());

}
//_____________________________________________________________________________
Double_t AliAnalysisTaskSigmaPlToProtonPiZero::GetQT(TLorentzVector sigmaVektor, TLorentzVector rekombinatedPi0)
{
	return  sin(sigmaVektor.Angle(rekombinatedPi0.Vect())) * rekombinatedPi0.P();
}
//________________________________________________________________________
Bool_t AliAnalysisTaskSigmaPlToProtonPiZero::IsRealProtonKink(AliESDkink* kink,AliStack* fMCStack)
{ //checks if a reconstructed kink is a pion kink
	Int_t Label1 = kink->GetLabel(0); //mother's track MC label (first component is mother)
	Int_t Label2 = kink->GetLabel(1); //daughter's track MC label (second component is daughter)
	//if (Label1>MCNTracks) return kFALSE;
	//if (Label2>MCNTracks) return kFALSE;

	TParticle *KinkSigmaCandidate = fMCStack->Particle(TMath::Abs(Label1)); //mother MC particle object
	TParticle *KinkProtonCandidate = fMCStack->Particle(TMath::Abs(Label2)); //daughter MC particle object

	Int_t code1 = KinkSigmaCandidate->GetPdgCode(); //mother's pdg code obtained from MC mother object
	Int_t code2 = KinkProtonCandidate->GetPdgCode(); //daughter's pdg code obtained from MC daughter object

	if( (code1==3222) && (code2==2212) ){
		return kTRUE;
	} else {
		return kFALSE;
	}
}
//________________________________________________________________________
Bool_t AliAnalysisTaskSigmaPlToProtonPiZero::IsRealPhoton(AliAODConversionPhoton PhotonCandidate, AliMCEvent* fMCEvent)
{ //checks if a reconstructed photon from sigma plus
	TParticle *Photon = NULL;
	if (PhotonCandidate.GetNCaloPhotonMCLabels() > 0){
		// Photon = PhotonCandidate.GetMCParticle(fMCEvent);
		Photon = fMCEvent->Particle(PhotonCandidate.GetCaloPhotonMCLabel(0));
		if (Photon) {
			if(Photon->GetPdgCode() == 22){
				TParticle* motherPart2 = (TParticle*)fMCEvent->Particle(Photon->GetMother(0));
				TParticle* grandmotherPart2 = (TParticle*)fMCEvent->Particle(motherPart2->GetMother(0));
				if(motherPart2->GetPdgCode() == 111 && grandmotherPart2->GetPdgCode() == 3222){
					return kTRUE;
				} else {
					return kFALSE;
				}
			} else if(Photon->GetPdgCode() == 11 || Photon->GetPdgCode() == -11){
				TParticle* motherPart2 = (TParticle*)fMCEvent->Particle(Photon->GetMother(0));
				TParticle* grandmotherPart2 = (TParticle*)fMCEvent->Particle(motherPart2->GetMother(0));
				TParticle* grandgrandmotherPart2 = (TParticle*)fMCEvent->Particle(grandmotherPart2->GetMother(0));
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
//_____________________________________________________________________________
void AliAnalysisTaskSigmaPlToProtonPiZero::Terminate(Option_t *)
{
	// terminate
	// called at the END of the analysis (when all events are processed)
}
//_____________________________________________________________________________
