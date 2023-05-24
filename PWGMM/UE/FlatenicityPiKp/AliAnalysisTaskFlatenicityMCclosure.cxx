/****************************************************************************
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
 *                                                                        *
 * Author: Antonio Ortiz (antonio.ortiz@nucleares.unam.mx)                *
 * Anatask to compute flatenicity (arXiv:2204.13733)                      *
 **************************************************************************/

class TTree;

class AliESDtrackCuts;
class AliESDAD; // AD

#include "AliAnalysisManager.h"
#include "AliAnalysisTask.h"
#include "AliAnalysisTaskESDfilter.h"
#include "AliAnalysisTaskSE.h"
#include "AliCentrality.h"
#include "AliESDAD.h" //AD
#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliESDUtils.h"
#include "AliESDVZERO.h"
#include "AliESDtrack.h"
#include "AliESDtrackCuts.h"
#include "AliEventCuts.h"
#include "AliGenCocktailEventHeader.h"
#include "AliGenEventHeader.h"
#include "AliInputEventHandler.h"
#include "AliLog.h"
#include "AliMCEvent.h"
#include "AliMCEventHandler.h"
#include "AliMCParticle.h"
#include "AliVZDC.h" //AD
#include <AliVAD.h>  //AD

#include "AliMultEstimator.h"
#include "AliMultInput.h"
#include "AliMultSelection.h"
#include "AliMultVariable.h"
#include "AliMultiplicity.h"
#include "AliOADBContainer.h"
#include "AliOADBMultSelection.h"
#include "AliStack.h"
#include "AliVEvent.h"
#include "AliVTrack.h"
#include "AliVVertex.h"
#include "TCanvas.h"
#include "TChain.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "THnSparse.h"
#include "TLegend.h"
#include "TList.h"
#include "TMath.h"
#include "TParticle.h"
#include "TProfile.h"
#include "TVector3.h"
#include <AliAnalysisFilter.h>
#include <AliESDVertex.h>
#include <AliHeader.h>
#include <AliMultiplicity.h>
#include <Riostream.h>
#include <TBits.h>
#include <TDirectory.h>
#include <TMath.h>
#include <TRandom.h>
#include <TTree.h>
using std::cout;
using std::endl;

#include "AliAnalysisTaskFlatenicityMCclosure.h"

/* static const Int_t nPtbins = 36; */
/* static const Double_t Ptbins[nPtbins + 1] = { */
/* 	0.0, 0.1, 0.15, 0.2,  0.25, 0.3,  0.35, 0.4,  0.45, 0.5,  0.6, 0.7, 0.8, */
/* 	0.9, 1.0, 1.25, 1.5,  2.0,  2.5,  3.0,  3.5,  4.0,  4.5,  5.0, 6.0, 7.0, */
/* 	8.0, 9.0, 10.0, 12.0, 14.0, 16.0, 18.0, 20.0, 30.0, 40.0, 50.0}; */

static const Int_t nPtbins = 13;
static const Double_t Ptbins[nPtbins + 1] = { 0.0, 0.15, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 6.0, 8.0, 10.0, 20.0};

static const Int_t nCent = 9;
static const Double_t centClass[nCent + 1] = {0.0,  1.0,  5.0,  10.0, 20.0, 30.0, 40.0, 50.0, 70.0, 100.0};

// const Int_t nCent = 2;
// Double_t centClass[nCent + 1] = {0.0, 50.0, 100.0};
static const Int_t nDet = 4;
static const Char_t *DetName[nDet] = {"ADC", "V0C", "V0A", "ADA"};
static const Int_t nComb = 3;
static const Char_t *CombName[nComb] = {"V0C_V0A", "ADC_ADA", "V0C_V0A_ADC_ADA"};

using namespace std; // std namespace: so you can do things like 'cout' etc

ClassImp(AliAnalysisTaskFlatenicityMCclosure) // classimp: necessary for root

AliAnalysisTaskFlatenicityMCclosure::AliAnalysisTaskFlatenicityMCclosure()
	: AliAnalysisTaskSE(), fESD(0), fEventCuts(0x0), fMCStack(0), fMC(0),
	fUseMC(kFALSE), fIsCalib(kFALSE), fIsEqualALICE(kTRUE), fVtxz(-1),
	fParVtx(0x0), fV0Mindex(-1), fmultTPC(-1), fmultV0A(-1),
	fmultV0C(-1), fmultADA(-1), fmultADC(-1), fmultTPCmc(-1), fmultV0Amc(-1),
	fmultV0Cmc(-1), fmultADAmc(-1), fmultADCmc(-1), fDetFlat("V0"),
	fRemoveTrivialScaling(kFALSE), fTrackFilter(0x0), fOutputList(0), fEtaCut(0.8),
	fPtMin(0.5), fv0mamplitude(0), fv0mpercentile(0), fFlat(-1), fFlatMC(-1), 
	fMultSelection(0x0), hFlatV0vsFlatTPC(0), hFlatenicityMC(0), hFlatResponse(0), 	
	hActivityV0DataSectBefore(0), hActivityV0DataSect(0), hV0vsVtxz(0),
	hActivityV0McSect(0), hFlatVsNchMC(0), hFlatVsV0M(0), hFlatMCVsV0M(0),
	hEtamc(0), hMultMCmVsV0M(0), hMultMCaVsV0M(0), hMultMCcVsV0M(0), hMultmVsV0M(0),
	hMultmVsV0Malice(0), hMultaVsV0M(0), hMultcVsV0M(0), hV0MBadruns(0),
	hMultV0AV0CvsFlat_BFTrigSel(0), hAmpV0AV0CvsFlat_BFTrigSel(0), hPercentileV0MvsFlat_BFTrigSel(0),
	hMultV0AV0CvsFlat_AFTrigSel(0), hAmpV0AV0CvsFlat_AFTrigSel(0), hPercentileV0MvsFlat_AFTrigSel(0),
	hMultV0AV0CvsFlatvspT_BFTrigSel(0), hAmpV0AV0CvsFlatvspT_BFTrigSel(0), hPercentileV0MvsFlatvspT_BFTrigSel(0),
	hMultV0AV0CvsFlatvspT_AFTrigSel(0), hAmpV0AV0CvsFlatvspT_AFTrigSel(0), hPercentileV0MvsFlatvspT_AFTrigSel(0),
	hAmpV0AV0CvsFlat(0), hPercentileV0MvsFlat(0), hAmpV0AV0CvsFlatvspT(0), hPercentileV0MvsFlatvspT(0), 
	hPtAll(0), hPtPrimaries(0), hPtSecondaries(0), hAmpV0AV0CalicevsFlat(0), hAmpV0AV0CalicevsFlatvspT(0),
	hAmpV0AV0CvsFlatvspT_pi(0), hAmpV0AV0CvsFlatvspT_k(0), hAmpV0AV0CvsFlatvspT_p(0),
	hMultV0AV0CvsFlatvspT_pi_BFTrigSel(0), hMultV0AV0CvsFlatvspT_k_BFTrigSel(0), hMultV0AV0CvsFlatvspT_p_BFTrigSel(0)
{
	for (Int_t i_d = 0; i_d < nDet; ++i_d) {
		hComponentsMult[i_d] = 0;
	}
	for (Int_t i_c = 0; i_c < nComb; ++i_c) {
		hCombinedMult[i_c] = 0;
	}
	for (Int_t i_d = 0; i_d < nDet; ++i_d) {
		hComponentsMultmc[i_d] = 0;
	}
	for (Int_t i_c = 0; i_c < nComb; ++i_c) {
		hCombinedMultmc[i_c] = 0;
	}
	for (Int_t i_c = 0; i_c < nComb; ++i_c) {
		hRmCombinedMult[i_c] = 0;
	}
	for (Int_t i_c = 0; i_c < nCent; ++i_c) {
		hMultMCmVsFlat[i_c] = 0;
	}
	for (Int_t i_c = 0; i_c < nCent; ++i_c) {
		hMultmVsFlat[i_c] = 0;
	}
}

//_____________________________________________________________________________
AliAnalysisTaskFlatenicityMCclosure::AliAnalysisTaskFlatenicityMCclosure(const char *name)
	: AliAnalysisTaskSE(name), fESD(0), fEventCuts(0x0), fMCStack(0), fMC(0),
	fUseMC(kFALSE), fIsCalib(kFALSE), fIsEqualALICE(kTRUE), fVtxz(-1),
	fParVtx(0x0), fV0Mindex(-1), fmultTPC(-1), fmultV0A(-1),
	fmultV0C(-1), fmultADA(-1), fmultADC(-1), fmultTPCmc(-1), fmultV0Amc(-1),
	fmultV0Cmc(-1), fmultADAmc(-1), fmultADCmc(-1), fDetFlat("V0"),
	fRemoveTrivialScaling(kFALSE), fTrackFilter(0x0), fOutputList(0), fEtaCut(0.8),
	fPtMin(0.5), fv0mamplitude(0), fv0mpercentile(0), fFlat(-1), fFlatMC(-1), 
	fMultSelection(0x0), hFlatV0vsFlatTPC(0), hFlatenicityMC(0), hFlatResponse(0), 
	hActivityV0DataSectBefore(0), hActivityV0DataSect(0), hV0vsVtxz(0),
	hActivityV0McSect(0), hFlatVsNchMC(0), hFlatVsV0M(0), hFlatMCVsV0M(0),
	hEtamc(0), hMultMCmVsV0M(0), hMultMCaVsV0M(0), hMultMCcVsV0M(0), hMultmVsV0M(0),
	hMultmVsV0Malice(0), hMultaVsV0M(0), hMultcVsV0M(0), hV0MBadruns(0),
	hMultV0AV0CvsFlat_BFTrigSel(0), hAmpV0AV0CvsFlat_BFTrigSel(0), hPercentileV0MvsFlat_BFTrigSel(0),
	hMultV0AV0CvsFlat_AFTrigSel(0), hAmpV0AV0CvsFlat_AFTrigSel(0), hPercentileV0MvsFlat_AFTrigSel(0),
	hMultV0AV0CvsFlatvspT_BFTrigSel(0), hAmpV0AV0CvsFlatvspT_BFTrigSel(0), hPercentileV0MvsFlatvspT_BFTrigSel(0), 
	hMultV0AV0CvsFlatvspT_AFTrigSel(0), hAmpV0AV0CvsFlatvspT_AFTrigSel(0), hPercentileV0MvsFlatvspT_AFTrigSel(0),
	hAmpV0AV0CvsFlat(0), hPercentileV0MvsFlat(0), hAmpV0AV0CvsFlatvspT(0), hPercentileV0MvsFlatvspT(0), 
	hPtAll(0), hPtPrimaries(0), hPtSecondaries(0), hAmpV0AV0CalicevsFlat(0), hAmpV0AV0CalicevsFlatvspT(0),
	hAmpV0AV0CvsFlatvspT_pi(0), hAmpV0AV0CvsFlatvspT_k(0), hAmpV0AV0CvsFlatvspT_p(0),
	hMultV0AV0CvsFlatvspT_pi_BFTrigSel(0), hMultV0AV0CvsFlatvspT_k_BFTrigSel(0), hMultV0AV0CvsFlatvspT_p_BFTrigSel(0)
{
	for (Int_t i_d = 0; i_d < nDet; ++i_d) {
		hComponentsMult[i_d] = 0;
	}
	for (Int_t i_c = 0; i_c < nComb; ++i_c) {
		hCombinedMult[i_c] = 0;
	}
	for (Int_t i_d = 0; i_d < nDet; ++i_d) {
		hComponentsMultmc[i_d] = 0;
	}
	for (Int_t i_c = 0; i_c < nComb; ++i_c) {
		hCombinedMultmc[i_c] = 0;
	}
	for (Int_t i_c = 0; i_c < nComb; ++i_c) {
		hRmCombinedMult[i_c] = 0;
	}
	for (Int_t i_c = 0; i_c < nCent; ++i_c) {
		hMultMCmVsFlat[i_c] = 0;
	}
	for (Int_t i_c = 0; i_c < nCent; ++i_c) {
		hMultmVsFlat[i_c] = 0;
	}

	DefineInput(0, TChain::Class()); // define the input of the analysis: in this
					 // case you take a 'chain' of events
	DefineOutput(1, TList::Class()); // define the ouptut of the analysis: in this
					 // case it's a list of histograms
}

//_____________________________________________________________________________
AliAnalysisTaskFlatenicityMCclosure::~AliAnalysisTaskFlatenicityMCclosure() {
	// destructor
	if (fOutputList) {
		delete fOutputList; // at the end of your task, it is deleted from memory by
				    // calling this function
		fOutputList = 0x0;
	}
}

//_____________________________________________________________________________
void AliAnalysisTaskFlatenicityMCclosure::UserCreateOutputObjects() {

	// create track filters
	fTrackFilter = new AliAnalysisFilter("trackFilter");

	AliESDtrackCuts* fCuts = new AliESDtrackCuts();
	fCuts->SetMinNCrossedRowsTPC(70);
	fCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);
	fCuts->SetMaxChi2PerClusterTPC(4);
	fCuts->SetAcceptKinkDaughters(kFALSE);
	fCuts->SetRequireTPCRefit(kTRUE);
	fCuts->SetRequireITSRefit(kTRUE);
	fCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kAny);
	fCuts->SetMaxDCAToVertexZ(2);
	fCuts->SetDCAToVertex2D(kFALSE);
	fCuts->SetRequireSigmaToVertex(kFALSE);
	fCuts->SetMaxChi2PerClusterITS(36);
	fCuts->SetMaxDCAToVertexXYPtDep("0.0105+0.0350/pt^1.1");
	fCuts->SetEtaRange(-0.8, 0.8);
	fTrackFilter->AddCuts(fCuts);

	float min_flat = -0.01;
	float max_flat = 1.01;
	int nbins_flat = 1020;
	if (fRemoveTrivialScaling) {
		min_flat = -0.1;
		max_flat = 9.9;
		nbins_flat = 2000;
	}

	const int nFlatbins = 1020;
	double Flatbins[nFlatbins+1] = {0.0};
	for (int i = 0; i <= nFlatbins; ++i) {
		Flatbins[i] = -0.01 + (double)i * 0.001;
	}

	const int nMultbinsForward = 801;
	double MultbinsForward[nMultbinsForward+1] = {0.0};
	for (int i = 0; i <= nMultbinsForward; ++i) {
		MultbinsForward[i] = -0.5 + (double)i;
	}

	const int nV0MAmplitudebins = 701;
	double V0MAmplitudebins[nV0MAmplitudebins+1] = {0.0};
	for (int i = 0; i <= nV0MAmplitudebins; ++i) {
		V0MAmplitudebins[i] = -0.5 + (double)i;
	}

	// create output objects
	fParVtx = new TF1("vtxpar", "pol2", -15, 15);
	fParVtx->SetParameters(89.8737, 0.127185, 0.00572492);

	OpenFile(1);
	fOutputList = new TList(); // this is a list which will contain all of your histograms
	fOutputList->SetOwner(kTRUE); // memory stuff: the list is owner of all

	hFlatV0vsFlatTPC = new TH2D("hFlatV0vsFlatTPC", "counter", nbins_flat, min_flat, max_flat, nbins_flat, min_flat, max_flat);
	fOutputList->Add(hFlatV0vsFlatTPC);

	hAmpV0AV0CvsFlat = new TH2F("hAmpV0AV0CvsFlat", "; V0M measured amplitude; Flattenicity",nV0MAmplitudebins,V0MAmplitudebins,nFlatbins,Flatbins);
	fOutputList->Add(hAmpV0AV0CvsFlat);

	hAmpV0AV0CalicevsFlat = new TH2F("hAmpV0AV0CalicevsFlat", "; V0M measured amplitude alice; Flattenicity",nV0MAmplitudebins,V0MAmplitudebins,nFlatbins,Flatbins);
	fOutputList->Add(hAmpV0AV0CalicevsFlat);

	hPercentileV0MvsFlat = new TH2F("hPercentileV0MvsFlat", "; V0M percentile; Flattenicity",nCent,centClass,nFlatbins,Flatbins);
	fOutputList->Add(hPercentileV0MvsFlat);

	hAmpV0AV0CvsFlatvspT = new TH3F("hAmpV0AV0CvsFlatvspT", "; V0M measured amplitude; Flattenicity; #it{p}_{T} (GeV/#it{c})",nV0MAmplitudebins,V0MAmplitudebins,nFlatbins,Flatbins,nPtbins,Ptbins);
	fOutputList->Add(hAmpV0AV0CvsFlatvspT);

	hAmpV0AV0CalicevsFlatvspT = new TH3F("hAmpV0AV0CalicevsFlatvspT", "; V0M measured amplitude alice; Flattenicity; #it{p}_{T} (GeV/#it{c})",nV0MAmplitudebins,V0MAmplitudebins,nFlatbins,Flatbins,nPtbins,Ptbins);
	fOutputList->Add(hAmpV0AV0CalicevsFlatvspT);

	hPercentileV0MvsFlatvspT = new TH3F("hPercentileV0MvsFlatvspT", "; V0M percentile; Flattenicity; #it{p}_{T} (GeV/#it{c})",nCent,centClass,nFlatbins,Flatbins,nPtbins,Ptbins);
	fOutputList->Add(hPercentileV0MvsFlatvspT);

	hAmpV0AV0CvsFlatvspT_pi = new TH3F("hAmpV0AV0CvsFlatvspT_pi", "; V0M measured amplitude; Flattenicity; #it{p}_{T} (GeV/#it{c})",nV0MAmplitudebins,V0MAmplitudebins,nFlatbins,Flatbins,nPtbins,Ptbins);
	fOutputList->Add(hAmpV0AV0CvsFlatvspT_pi);

	hAmpV0AV0CvsFlatvspT_k = new TH3F("hAmpV0AV0CvsFlatvspT_k", "; V0M measured amplitude; Flattenicity; #it{p}_{T} (GeV/#it{c})",nV0MAmplitudebins,V0MAmplitudebins,nFlatbins,Flatbins,nPtbins,Ptbins);
	fOutputList->Add(hAmpV0AV0CvsFlatvspT_k);

	hAmpV0AV0CvsFlatvspT_p = new TH3F("hAmpV0AV0CvsFlatvspT_p", "; V0M measured amplitude; Flattenicity; #it{p}_{T} (GeV/#it{c})",nV0MAmplitudebins,V0MAmplitudebins,nFlatbins,Flatbins,nPtbins,Ptbins);
	fOutputList->Add(hAmpV0AV0CvsFlatvspT_p);

	hPtAll = new TH1D("hPtAll", "All charged particles; #it{p}_{T} (GeV/#it{c}); Entries",nPtbins,Ptbins);
	fOutputList->Add(hPtAll);

	hPtPrimaries = new TH1D("hPtPrimaries", "Only primary charged particles; #it{p}_{T} (GeV/#it{c}); Entries",nPtbins,Ptbins);
	fOutputList->Add(hPtPrimaries);

	hPtSecondaries = new TH1D("hPtSecondaries", "Only secondary charged particles; #it{p}_{T} (GeV/#it{c}); Entries",nPtbins,Ptbins);
	fOutputList->Add(hPtSecondaries);

	for (Int_t i_d = 0; i_d < nDet; ++i_d) {
		hComponentsMult[i_d] = new TH2D(Form("hAmpl_%s", DetName[i_d]), "", 5000,
				-0.5, 20000.0, 200, -0.5, 199.5);
		fOutputList->Add(hComponentsMult[i_d]);
	}
	for (Int_t i_c = 0; i_c < nComb; ++i_c) {
		hCombinedMult[i_c] = new TH2D(Form("hCombined_%s", CombName[i_c]), "", 500,
				-0.5, 499.5, 200, -0.5, 199.5);
		fOutputList->Add(hCombinedMult[i_c]);
	}

	if (fUseMC) {

		hEtamc = new TH1D("hEtamc", "", 140, -7.0, 7.0);
		fOutputList->Add(hEtamc);

		hFlatenicityMC = new TH1D("hFlatenicityMC", "counter", nbins_flat, min_flat, max_flat);
		fOutputList->Add(hFlatenicityMC);
		hFlatResponse = new TH2D("hFlatResponse", "; true flat; measured flat", nbins_flat, min_flat, max_flat, nbins_flat, min_flat, max_flat);
		fOutputList->Add(hFlatResponse);

		hFlatVsNchMC = new TH2D("hFlatVsNchMC", "; true flat; true Nch", nbins_flat,
				min_flat, max_flat, 100, -0.5, 99.5);
		fOutputList->Add(hFlatVsNchMC);

		for (Int_t i_d = 0; i_d < nDet; ++i_d) {
			hComponentsMultmc[i_d] = new TH2D(Form("hTrueMult_%s", DetName[i_d]), "",
					600, -0.5, 599.0, 200, -0.5, 199.5);
			fOutputList->Add(hComponentsMultmc[i_d]);
		}
		for (Int_t i_c = 0; i_c < nComb; ++i_c) {
			hCombinedMultmc[i_c] = new TH2D(Form("hTrueCombined_%s", CombName[i_c]),
					"", 500, -0.5, 499.5, 200, -0.5, 199.5);
			fOutputList->Add(hCombinedMultmc[i_c]);
		}
		for (Int_t i_c = 0; i_c < nComb; ++i_c) {
			hRmCombinedMult[i_c] =
				new TH2D(Form("hRmCombined_%s", CombName[i_c]),
						"; measured combined mult.; true combined mult.", 500, -0.5,
						499.5, 500, -0.5, 499.5);
			fOutputList->Add(hRmCombinedMult[i_c]);
		}

		hFlatMCVsV0M = new TH2D("hFlatMCVsV0M", "", nCent, centClass, nbins_flat, min_flat, max_flat);
		fOutputList->Add(hFlatMCVsV0M);

		for (Int_t i_c = 0; i_c < nCent; ++i_c) {
			hMultMCmVsFlat[i_c] = new TH2D(Form("hMultMCmVsFlat_c%d", i_c), "", 100,
					0.0, 1.0, 1000, -0.5, 999.5);
			fOutputList->Add(hMultMCmVsFlat[i_c]);
		}
		hMultMCmVsV0M = new TH2D("hMultMCmVsV0M", "", nCent, centClass, 1000, -0.5, 999.5);
		fOutputList->Add(hMultMCmVsV0M);

		hMultMCaVsV0M = new TH2D("hMultMCaVsV0M", "", nCent, centClass, 1000, -0.5, 999.5);
		fOutputList->Add(hMultMCaVsV0M);

		hMultMCcVsV0M = new TH2D("hMultMCcVsV0M", "", nCent, centClass, 1000, -0.5, 999.5);
		fOutputList->Add(hMultMCcVsV0M);

		hMultV0AV0CvsFlat_BFTrigSel = new TH2F("hMultV0AV0CvsFlat_BFTrigSel","Before Trigger selection; True #it{N}_{ch} in the V0A+V0C acceptance; True Flattenicity",nMultbinsForward,MultbinsForward,nFlatbins,Flatbins);
		fOutputList->Add(hMultV0AV0CvsFlat_BFTrigSel);

		hAmpV0AV0CvsFlat_BFTrigSel = new TH2F("hAmpV0AV0CvsFlat_BFTrigSel","Before Trigger selection; V0M measured amplitude; True Flattenicity",nV0MAmplitudebins,V0MAmplitudebins,nFlatbins,Flatbins);
		fOutputList->Add(hAmpV0AV0CvsFlat_BFTrigSel);

		hPercentileV0MvsFlat_BFTrigSel = new TH2F("hPercentileV0MvsFlat_BFTrigSel","Before Trigger selection; V0M percentile; True Flattenicity",nCent,centClass,nFlatbins,Flatbins);
		fOutputList->Add(hPercentileV0MvsFlat_BFTrigSel);

		hMultV0AV0CvsFlat_AFTrigSel = new TH2F("hMultV0AV0CvsFlat_AFTrigSel","After Trigger selection; True #it{N}_{ch} in the V0A+V0C acceptance; True Flattenicity",nMultbinsForward,MultbinsForward,nFlatbins,Flatbins);
		fOutputList->Add(hMultV0AV0CvsFlat_AFTrigSel);

		hAmpV0AV0CvsFlat_AFTrigSel = new TH2F("hAmpV0AV0CvsFlat_AFTrigSel","After Trigger selection; V0M measured amplitude; True Flattenicity",nV0MAmplitudebins,V0MAmplitudebins,nFlatbins,Flatbins);
		fOutputList->Add(hAmpV0AV0CvsFlat_AFTrigSel);

		hPercentileV0MvsFlat_AFTrigSel = new TH2F("hPercentileV0MvsFlat_AFTrigSel","After Trigger selection; V0M percentile; True Flattenicity",nCent,centClass,nFlatbins,Flatbins);
		fOutputList->Add(hPercentileV0MvsFlat_AFTrigSel);

		hMultV0AV0CvsFlatvspT_BFTrigSel = new TH3F("hMultV0AV0CvsFlatvspT_BFTrigSel","Before Trigger selection; True #it{N}_{ch} in the V0A+V0C acceptance; True Flattenicity; True #it{p}_{T} (GeV/#it{c})",nMultbinsForward,MultbinsForward,nFlatbins,Flatbins,nPtbins,Ptbins);
		fOutputList->Add(hMultV0AV0CvsFlatvspT_BFTrigSel);

		hAmpV0AV0CvsFlatvspT_BFTrigSel = new TH3F("hAmpV0AV0CvsFlatvspT_BFTrigSel","Before Trigger selection; V0M measured amplitude; True Flattenicity; True #it{p}_{T} (GeV/#it{c})",nV0MAmplitudebins,V0MAmplitudebins,nFlatbins,Flatbins,nPtbins,Ptbins);
		fOutputList->Add(hAmpV0AV0CvsFlatvspT_BFTrigSel);

		hPercentileV0MvsFlatvspT_BFTrigSel = new TH3F("hPercentileV0MvsFlatvspT_BFTrigSel","Before Trigger selection; V0M percentile; True Flattenicity; True #it{p}_{T} (GeV/#it{c})",nCent,centClass,nFlatbins,Flatbins,nPtbins,Ptbins);
		fOutputList->Add(hPercentileV0MvsFlatvspT_BFTrigSel);

		hMultV0AV0CvsFlatvspT_AFTrigSel = new TH3F("hMultV0AV0CvsFlatvspT_AFTrigSel","After Trigger selection; True #it{N}_{ch} in the V0A+V0C acceptance; True Flattenicity; True #it{p}_{T} (GeV/#it{c})",nMultbinsForward,MultbinsForward,nFlatbins,Flatbins,nPtbins,Ptbins);
		fOutputList->Add(hMultV0AV0CvsFlatvspT_AFTrigSel);

		hAmpV0AV0CvsFlatvspT_AFTrigSel = new TH3F("hAmpV0AV0CvsFlatvspT_AFTrigSel","After Trigger selection; V0M measured amplitude; True Flattenicity; True #it{p}_{T} (GeV/#it{c})",nV0MAmplitudebins,V0MAmplitudebins,nFlatbins,Flatbins,nPtbins,Ptbins);
		fOutputList->Add(hAmpV0AV0CvsFlatvspT_AFTrigSel);

		hPercentileV0MvsFlatvspT_AFTrigSel = new TH3F("hPercentileV0MvsFlatvspT_AFTrigSel","After Trigger selection; V0M percentile; True Flattenicity; True #it{p}_{T} (GeV/#it{c})",nCent,centClass,nFlatbins,Flatbins,nPtbins,Ptbins);
		fOutputList->Add(hPercentileV0MvsFlatvspT_AFTrigSel);

		hMultV0AV0CvsFlatvspT_pi_BFTrigSel = new TH3F("hMultV0AV0CvsFlatvspT_pi_BFTrigSel","Before Trigger selection; True #it{N}_{ch} in the V0A+V0C acceptance; True Flattenicity; True #it{p}_{T} (GeV/#it{c})",nMultbinsForward,MultbinsForward,nFlatbins,Flatbins,nPtbins,Ptbins);
		fOutputList->Add(hMultV0AV0CvsFlatvspT_pi_BFTrigSel);

		hMultV0AV0CvsFlatvspT_k_BFTrigSel = new TH3F("hMultV0AV0CvsFlatvspT_k_BFTrigSel","Before Trigger selection; True #it{N}_{ch} in the V0A+V0C acceptance; True Flattenicity; True #it{p}_{T} (GeV/#it{c})",nMultbinsForward,MultbinsForward,nFlatbins,Flatbins,nPtbins,Ptbins);
		fOutputList->Add(hMultV0AV0CvsFlatvspT_k_BFTrigSel);

		hMultV0AV0CvsFlatvspT_p_BFTrigSel = new TH3F("hMultV0AV0CvsFlatvspT_p_BFTrigSel","Before Trigger selection; True #it{N}_{ch} in the V0A+V0C acceptance; True Flattenicity; True #it{p}_{T} (GeV/#it{c})",nMultbinsForward,MultbinsForward,nFlatbins,Flatbins,nPtbins,Ptbins);
		fOutputList->Add(hMultV0AV0CvsFlatvspT_p_BFTrigSel);
	}

	hActivityV0DataSectBefore = new TProfile(
			"hActivityV0DataSectBefore",
			"rec; V0 sector; (before calib) #LTmultiplicity#GT", 64, -0.5, 63.5);
	fOutputList->Add(hActivityV0DataSectBefore);

	hActivityV0DataSect = new TProfile(
			"hActivityV0DataSect", "rec; V0 sector; (after calib) #LTmultiplicity#GT",
			64, -0.5, 63.5);
	fOutputList->Add(hActivityV0DataSect);

	hV0vsVtxz = new TProfile("hV0vsVtxz", ";total amplitude; vtx_z", 30, -15, 15);
	fOutputList->Add(hV0vsVtxz);

	if (fUseMC) {
		hActivityV0McSect =
			new TProfile("hActivityV0McSect", "true; V0 sector; #LTmultiplicity#GT",
					64, -0.5, 63.5);
		fOutputList->Add(hActivityV0McSect);
	}

	hFlatVsV0M = new TH2D("hFlatVsV0M", "", nCent, centClass, nbins_flat,
			min_flat, max_flat);
	fOutputList->Add(hFlatVsV0M);

	for (Int_t i_c = 0; i_c < nCent; ++i_c) {
		hMultmVsFlat[i_c] = new TH2D(Form("hMultmVsFlat_c%d", i_c), "",nV0MAmplitudebins,V0MAmplitudebins,nFlatbins,Flatbins);
		fOutputList->Add(hMultmVsFlat[i_c]);
	}

	hMultmVsV0M =
		new TH2D("hMultmVsV0M", "", nCent, centClass, 1000, -0.5, 999.5);
	fOutputList->Add(hMultmVsV0M);

	hMultmVsV0Malice =
		new TH2D("hMultmVsV0Malice", "", nCent, centClass, 1000, -0.5, 999.5);
	fOutputList->Add(hMultmVsV0Malice);

	hMultaVsV0M =
		new TH2D("hMultaVsV0M", "", nCent, centClass, 1000, -0.5, 999.5);
	fOutputList->Add(hMultaVsV0M);

	hMultcVsV0M =
		new TH2D("hMultcVsV0M", "", nCent, centClass, 1000, -0.5, 999.5);
	fOutputList->Add(hMultcVsV0M);

	hV0MBadruns = new TH1D("hV0MBadruns", "", 1000, -0.5, 999.5);
	fOutputList->Add(hV0MBadruns);

	fEventCuts.AddQAplotsToList(fOutputList);
	PostData(1, fOutputList); // postdata will notify the analysis manager of
				  // changes / updates to the
}

//_____________________________________________________________________________
void AliAnalysisTaskFlatenicityMCclosure::UserExec(Option_t *) {

	AliVEvent *event = InputEvent();
	if (!event) {
		Error("UserExec", "Could not retrieve event");
		return;
	}

	fESD = dynamic_cast<AliESDEvent *>(event);

	if (!fESD) {
		Printf("%s:%d ESDEvent not found in Input Manager", (char *)__FILE__,
				__LINE__);
		this->Dump();
		return;
	}

	if (fUseMC) {
		//      E S D
		fMC = dynamic_cast<AliMCEvent *>(MCEvent());
		if (!fMC) {
			Printf("%s:%d MCEvent not found in Input Manager", (char *)__FILE__,
					__LINE__);
			this->Dump();
			return;
		}
		fMCStack = fMC->Stack();
	}

	AliHeader *headerMC;
	Bool_t isGoodVtxPosMC = kFALSE;

	if (fUseMC) {
		headerMC = fMC->Header();
		AliGenEventHeader *genHeader = headerMC->GenEventHeader();
		TArrayF vtxMC(3); // primary vertex  MC
		vtxMC[0] = 9999;
		vtxMC[1] = 9999;
		vtxMC[2] = 9999; // initialize with dummy
		if (genHeader)
			genHeader->PrimaryVertex(vtxMC);

		if (TMath::Abs(vtxMC[2]) <= 10)
			isGoodVtxPosMC = kTRUE;
	}

	fMultSelection = (AliMultSelection *)fESD->FindListObject("MultSelection");
	if (!fMultSelection)
		cout << "------- No AliMultSelection Object Found --------"
			<< fMultSelection << endl;

	fv0mpercentile = -999;
	fv0mpercentile = fMultSelection->GetMultiplicityPercentile("V0M");
	fv0mamplitude   = fMultSelection->GetEstimator("V0M")->GetValue();
	int v0multalice = fMultSelection->GetEstimator("V0M")->GetValue();

	if (fUseMC) {
		fFlatMC = -1;
		fFlatMC = GetFlatenicityMC();
	}

	if ((fUseMC) && (isGoodVtxPosMC)) {
		ExtractMultiplicitiesMC();
		MCkinematics();

		if ((fFlatMC >= 0) && (fmultV0Cmc) > 0 && (fmultV0Amc > 0) && (fmultTPCmc > 0)) { 
			// This condition ensures that at least one charged particle 
			// went into the acceptance of the V0A or V0C and there is 
			// at least one charged particle in |eta|<0.8
			hMultV0AV0CvsFlat_BFTrigSel->Fill(fmultV0Amc + fmultV0Cmc,fFlatMC);
			hAmpV0AV0CvsFlat_BFTrigSel->Fill(fv0mamplitude,fFlatMC);
			hPercentileV0MvsFlat_BFTrigSel->Fill(fv0mpercentile,fFlatMC);
			GetMCchargedTrueDists(); 
		}

	}

	// Trigger selection
	UInt_t fSelectMask = fInputHandler->IsEventSelected();
	Bool_t isINT7selected = fSelectMask & AliVEvent::kINT7;
	if (!isINT7selected)
		return;

	// Good events
	if (!fEventCuts.AcceptEvent(event)) {
		PostData(1, fOutputList);
		return;
	}

	// Good vertex
	Bool_t hasRecVertex = kFALSE;
	hasRecVertex = HasRecVertex();
	if (!hasRecVertex)
		return;

	for (Int_t i_c = 0; i_c < nCent; ++i_c) {
		if (fv0mpercentile >= centClass[i_c] &&
				fv0mpercentile < centClass[i_c + 1]) {
			fV0Mindex = i_c;
		} else {
			continue;
		}
	}

	Double_t flatenicity_v0 = -1;
	Double_t flatenicity_tpc = GetFlatenicityTPC();
	if (fIsEqualALICE) {
		flatenicity_v0 = GetFlatenicityV0EqualALICE();
		ExtractMultiplicitiesEqualALICE();
	} else {
		flatenicity_v0 = GetFlatenicityV0();
		ExtractMultiplicities();
	}

	float com1mc = 0;
	float com2mc = 0;
	float com3mc = 0;
	if (fUseMC) {
		if (isGoodVtxPosMC) {
			ExtractMultiplicitiesMC();

			float activityMC[4] = {0, 0, 0, 0};
			activityMC[0] = fmultADCmc;
			activityMC[1] = fmultV0Cmc;
			activityMC[2] = fmultV0Amc;
			activityMC[3] = fmultADAmc;
			for (int i_a = 0; i_a < 4; ++i_a) {
				hComponentsMultmc[i_a]->Fill(activityMC[i_a], fmultTPCmc);
			}
			com1mc = fmultV0Amc + fmultV0Cmc;
			com2mc = fmultADAmc + fmultADCmc;
			com3mc = com1mc + com2mc;
			hCombinedMultmc[0]->Fill(com1mc, fmultTPCmc);
			hCombinedMultmc[1]->Fill(com2mc, fmultTPCmc);
			hCombinedMultmc[2]->Fill(com3mc, fmultTPCmc);
		}
	}

	// these values were obtained from LHC16l pass 2 (and MC)
	float avData[4] = {1819.91, 55.6384, 27.6564, 449.373};
	float avExpect[4] = {9.18629, 15.1672, 11.934, 7.47469};
	if (fUseMC) {
		avData[0] = 1380.66;
		avData[1] = 60.009;
		avData[2] = 35.1942;
		avData[3] = 202.348;
	}
	float weigths[4] = {1.0, 1.0, 1.0, 1.0};
	for (int i_a = 0; i_a < 4; ++i_a) {
		weigths[i_a] = avExpect[i_a] / avData[i_a];
	}

	float activity[4] = {0, 0, 0, 0};
	activity[0] = fmultADC;
	activity[1] = fmultV0C;
	activity[2] = fmultV0A;
	activity[3] = fmultADA;
	for (int i_a = 0; i_a < 4; ++i_a) {
		hComponentsMult[i_a]->Fill(activity[i_a], fmultTPC);
	}
	float com1 = weigths[1] * activity[1] + activity[2] * weigths[2];
	float com2 = weigths[0] * activity[0] + activity[3] * weigths[3];
	float com3 = com1 + com2;
	hCombinedMult[0]->Fill(com1, fmultTPC);
	hCombinedMult[1]->Fill(com2, fmultTPC);
	hCombinedMult[2]->Fill(com3, fmultTPC);

	if (fUseMC) {
		hRmCombinedMult[0]->Fill(com1, com1mc);
		hRmCombinedMult[1]->Fill(com2, com2mc);
		hRmCombinedMult[2]->Fill(com3, com3mc);
	}

	fFlat = flatenicity_v0; // default V0
	if (fDetFlat == "VO_TPC") {
		fFlat = (flatenicity_v0 + flatenicity_tpc) / 2.0;
	}
	if (fDetFlat == "TPC") {
		fFlat = flatenicity_tpc;
	}
	if (fDetFlat == "V0") {
		fFlat = flatenicity_v0;
	}

	// Here the True Flattenicity is calculated
	// However, this has the effects of the Trigger selection
	fFlatMC = -1;
	if ((fUseMC) && (fmultV0Cmc) > 0 && (fmultV0Amc > 0)) {
		fFlatMC = GetFlatenicityMC();
		if (fFlatMC >= 0) {
			hFlatenicityMC->Fill(fFlatMC);
			hFlatResponse->Fill(fFlatMC, fFlat);
			hFlatMCVsV0M->Fill(fv0mpercentile, fFlatMC);
			if (fV0Mindex >= 0) {
				hMultMCmVsFlat[fV0Mindex]->Fill(fFlatMC, fmultV0Cmc + fmultV0Amc);
			}
			hMultMCmVsV0M->Fill(fv0mpercentile, fmultV0Cmc + fmultV0Amc);
			hMultMCcVsV0M->Fill(fv0mpercentile, fmultV0Cmc);
			hMultMCaVsV0M->Fill(fv0mpercentile, fmultV0Amc);

			hMultV0AV0CvsFlat_AFTrigSel->Fill(fmultV0Amc + fmultV0Cmc,fFlatMC); 
			hAmpV0AV0CvsFlat_AFTrigSel->Fill(fv0mamplitude,fFlatMC); 
			hPercentileV0MvsFlat_AFTrigSel->Fill(fv0mpercentile,fFlatMC);
			MakeMCanalysis();
		}
	}

	hFlatV0vsFlatTPC->Fill(flatenicity_tpc, flatenicity_v0);

	if ((fFlat >= 0) && (fmultV0C) > 0 && (fmultV0A > 0)) {
		if (fV0Mindex >= 0) {
			if ((fV0Mindex == nCent - 1) && (v0multalice > 400)) {
				hV0MBadruns->Fill(v0multalice);
			} else {
				hFlatVsV0M->Fill(fv0mpercentile, fFlat);
				hMultmVsFlat[fV0Mindex]->Fill(fmultV0C + fmultV0A,fFlat);
				hMultmVsV0M->Fill(fv0mpercentile, fmultV0C + fmultV0A);
				hMultmVsV0Malice->Fill(fv0mpercentile, v0multalice);
				hMultcVsV0M->Fill(fv0mpercentile, fmultV0C);
				hMultaVsV0M->Fill(fv0mpercentile, fmultV0A);
			
				hAmpV0AV0CvsFlat->Fill(fmultV0A + fmultV0C, fFlat);
				hAmpV0AV0CalicevsFlat->Fill(fv0mamplitude, fFlat);
				hPercentileV0MvsFlat->Fill(fv0mpercentile, fFlat);
				MakeDataanalysis();
			}
		}
	}

	PostData(1, fOutputList); // stream the result of this event to the output
				  // manager which will write it to a file
}

//______________________________________________________________________________
void AliAnalysisTaskFlatenicityMCclosure::Terminate(Option_t *) {}

//______________________________________________________________________________
void AliAnalysisTaskFlatenicityMCclosure::MakeDataanalysis() {

	Int_t nTracks = fESD->GetNumberOfTracks();
	for (Int_t iT = 0; iT < nTracks; ++iT) {

		AliESDtrack *esdtrack = static_cast<AliESDtrack *>(
				fESD->GetTrack(iT)); // get a track (type AliesdTrack)
		if (!esdtrack)
			continue;
		if (!fTrackFilter->IsSelected(esdtrack))
			continue;
		if (TMath::Abs(esdtrack->Eta()) > fEtaCut)
			continue;
		if (esdtrack->Pt() < fPtMin)
			continue;

		Int_t mcLabel = -1;
		mcLabel = TMath::Abs(esdtrack->GetLabel());

		hPtAll->Fill(esdtrack->Pt());
		if(fMC->IsPhysicalPrimary(mcLabel)) { hPtPrimaries->Fill(esdtrack->Pt()); }
		else { hPtSecondaries->Fill(esdtrack->Pt()); }

		hAmpV0AV0CalicevsFlatvspT->Fill(fv0mamplitude, fFlat, esdtrack->Pt());
		hAmpV0AV0CvsFlatvspT->Fill(fmultV0A + fmultV0C, fFlat, esdtrack->Pt());
		hPercentileV0MvsFlatvspT->Fill(fv0mpercentile, fFlat, esdtrack->Pt());

		if (TMath::Abs(esdtrack->Charge())==0)
			continue;
		
		AliMCParticle *particle = nullptr;
		particle = (AliMCParticle*)fMC->GetTrack(mcLabel);

		if (!particle)
			continue;

		int pdgCode = particle->PdgCode();
		int pidCode = GetPidCode(pdgCode);

		if (pidCode==1) { hAmpV0AV0CvsFlatvspT_pi->Fill(fmultV0A + fmultV0C, fFlat, esdtrack->Pt()); }
		if (pidCode==2) { hAmpV0AV0CvsFlatvspT_k->Fill(fmultV0A + fmultV0C, fFlat, esdtrack->Pt()); }
		if ((pidCode==3) && (fMC->IsPhysicalPrimary(mcLabel))) { hAmpV0AV0CvsFlatvspT_p->Fill(fmultV0A + fmultV0C, fFlat, esdtrack->Pt()); }

	}
}

//______________________________________________________________________________
void AliAnalysisTaskFlatenicityMCclosure::MakeMCanalysis() {

	for (Int_t i = 0; i < fMC->GetNumberOfTracks(); ++i) {

		AliMCParticle *particle = (AliMCParticle *)fMC->GetTrack(i);
		if (!particle)
			continue;
		if (!fMC->IsPhysicalPrimary(i))
			continue;
		if (TMath::Abs(particle->Eta()) > fEtaCut)
			continue;
		if (particle->Pt() < fPtMin)
			continue;
		if (TMath::Abs(particle->Charge()) < 0.1)
			continue;
		if (AliAnalysisUtils::IsParticleFromOutOfBunchPileupCollision(i, fMC))
			continue;

		hMultV0AV0CvsFlatvspT_AFTrigSel->Fill(fmultV0Amc + fmultV0Cmc,fFlatMC,particle->Pt());
		hAmpV0AV0CvsFlatvspT_AFTrigSel->Fill(fv0mamplitude,fFlatMC,particle->Pt());
		hPercentileV0MvsFlatvspT_AFTrigSel->Fill(fv0mpercentile,fFlatMC,particle->Pt());
	}
}
//______________________________________________________________________________
Double_t AliAnalysisTaskFlatenicityMCclosure::GetFlatenicityTPC() {

	const int nRings2 = 4;
	const int nSectors2 = 8;
	const int nCells2 = nRings2 * nSectors2;
	float maxEta2[nRings2] = {-0.4, 0.0, +0.4, +0.8};
	float minEta2[nRings2] = {-0.8, -0.4, +0.0, +0.4};
	float maxPhi2[nSectors2] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0};
	float minPhi2[nSectors2] = {0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0};
	float RhoLattice2[nCells2];
	for (int iCh = 0; iCh < nCells2; iCh++) {
		RhoLattice2[iCh] = 0.0;
	}
	int mult_glob = 0;
	Int_t nTracks = fESD->GetNumberOfTracks();
	for (Int_t iT = 0; iT < nTracks; ++iT) {

		AliESDtrack *esdtrack = static_cast<AliESDtrack *>(
				fESD->GetTrack(iT)); // get a track (type AliesdTrack)
		if (!esdtrack)
			continue;
		if (!fTrackFilter->IsSelected(esdtrack))
			continue;
		float eta_a = esdtrack->Eta();
		float phi_a = esdtrack->Phi();

		if (TMath::Abs(eta_a) > fEtaCut)
			continue;
		if (esdtrack->Pt() < fPtMin)
			continue;
		int i_ch = 0;
		for (int ir = 0; ir < nRings2; ir++) {
			for (int is = 0; is < nSectors2; is++) {
				if (eta_a >= minEta2[ir] && eta_a < maxEta2[ir] &&
						phi_a >= minPhi2[is] * 2.0 * M_PI / (1.0 * nSectors2) &&
						phi_a < maxPhi2[is] * 2.0 * M_PI / (1.0 * nSectors2)) {
					RhoLattice2[i_ch]++;
					mult_glob++;
				}
				i_ch++;
			}
		}
	}

	double mRho_glob = 0;
	for (int iCell = 0; iCell < nCells2; ++iCell) {
		mRho_glob += 1.0 * RhoLattice2[iCell];
	}
	// average activity per cell
	mRho_glob /= (1.0 * nCells2);
	// get sigma
	double sRho_glob_tmp = 0;
	for (int iCell = 0; iCell < nCells2; ++iCell) {
		sRho_glob_tmp += TMath::Power(1.0 * RhoLattice2[iCell] - mRho_glob, 2);
	}
	sRho_glob_tmp /= (1.0 * nCells2 * nCells2);
	double sRho_glob = TMath::Sqrt(sRho_glob_tmp);
	float flatenicity_glob = 9999;
	if (mRho_glob > 0) {
		if (fRemoveTrivialScaling) {
			flatenicity_glob = TMath::Sqrt(mult_glob) * sRho_glob / mRho_glob;
		} else {
			flatenicity_glob = sRho_glob / mRho_glob;
		}
	}

	return flatenicity_glob;
}

//______________________________________________________________________________
Double_t AliAnalysisTaskFlatenicityMCclosure::GetFlatenicityV0EqualALICE() {

	AliVVZERO *lVV0 = 0x0;
	AliVEvent *lVevent = 0x0;
	lVevent = dynamic_cast<AliVEvent *>(InputEvent());
	if (!lVevent) {
		AliWarning("ERROR: ESD / AOD event not available \n");
		return -1;
	}
	// Get VZERO Information for multiplicity later
	lVV0 = lVevent->GetVZEROData();
	if (!lVV0) {
		AliError("AliVVZERO not available");
		return 9999;
	}
	// Flatenicity calculation
	const Int_t nRings = 4;
	const Int_t nSectors = 8;
	Float_t minEtaV0C[nRings] = {-3.7, -3.2, -2.7, -2.2};
	Float_t maxEtaV0C[nRings] = {-3.2, -2.7, -2.2, -1.7};
	Float_t maxEtaV0A[nRings] = {5.1, 4.5, 3.9, 3.4};
	Float_t minEtaV0A[nRings] = {4.5, 3.9, 3.4, 2.8};
	// Grid
	const Int_t nCells = nRings * 2 * nSectors;
	float RhoLattice[nCells];
	for (Int_t iCh = 0; iCh < nCells; iCh++) {
		RhoLattice[iCh] = 0.0;
	}

	// before calibration
	for (Int_t iCh = 0; iCh < nCells; iCh++) {
		float mult = 0;
		// only corrected for vertex
		if (iCh < 32) { // V0C
			mult = AliESDUtils::GetCorrV0C(lVV0->GetMultiplicity(iCh), fVtxz);
		} else { // V0A
			mult = AliESDUtils::GetCorrV0A(lVV0->GetMultiplicity(iCh), fVtxz);
		}

		RhoLattice[iCh] = mult;
		hActivityV0DataSectBefore->Fill(iCh, lVV0->GetMultiplicity(iCh));
	}

	// Filling histos with mult info
	float total_v0_tmp = 0;
	for (Int_t iCh = 0; iCh < nCells; iCh++) {
		hActivityV0DataSect->Fill(iCh, RhoLattice[iCh]);
		total_v0_tmp += RhoLattice[iCh];
	}
	int total_v0 = total_v0_tmp;

	hV0vsVtxz->Fill(fVtxz, total_v0);

	Int_t nringA = 0;
	Int_t nringC = 0;
	for (Int_t iCh = 0; iCh < nCells; iCh++) {
		Float_t detaV0 = -1;
		// Float_t mult = lVV0->GetMultiplicity(iCh);
		if (iCh < 32) { // V0C
			if (iCh < 8) {
				nringC = 0;
			} else if (iCh >= 8 && iCh < 16) {
				nringC = 1;
			} else if (iCh >= 16 && iCh < 24) {
				nringC = 2;
			} else {
				nringC = 3;
			}
			detaV0 = maxEtaV0C[nringC] - minEtaV0C[nringC];
		} else { // V0A
			if (iCh < 40) {
				nringA = 0;
			} else if (iCh >= 40 && iCh < 48) {
				nringA = 1;
			} else if (iCh >= 48 && iCh < 56) {
				nringA = 2;
			} else {
				nringA = 3;
			}
			detaV0 = maxEtaV0A[nringA] - minEtaV0A[nringA];
		}
		// consider the different eta coverage
		RhoLattice[iCh] /= detaV0;
	}
	Float_t mRho = 0;
	Float_t flatenicity = -1;
	for (Int_t iCh = 0; iCh < nCells; iCh++) {
		mRho += RhoLattice[iCh];
	}
	Float_t multiplicityV0M = mRho;
	// average activity per cell
	mRho /= (1.0 * nCells);
	// get sigma
	Double_t sRho_tmp = 0;
	for (Int_t iCh = 0; iCh < nCells; iCh++) {
		sRho_tmp += TMath::Power(1.0 * RhoLattice[iCh] - mRho, 2);
	}
	sRho_tmp /= (1.0 * nCells * nCells);
	Float_t sRho = TMath::Sqrt(sRho_tmp);
	if (mRho > 0) {
		if (fRemoveTrivialScaling) {
			flatenicity = TMath::Sqrt(multiplicityV0M) * sRho / mRho;
		} else {
			flatenicity = sRho / mRho;
		}
	} else {
		flatenicity = 9999;
	}
	return flatenicity;
}

//______________________________________________________________________________
Double_t AliAnalysisTaskFlatenicityMCclosure::GetFlatenicityV0() {

	AliVVZERO *lVV0 = 0x0;
	AliVEvent *lVevent = 0x0;
	lVevent = dynamic_cast<AliVEvent *>(InputEvent());
	if (!lVevent) {
		AliWarning("ERROR: ESD / AOD event not available \n");
		return -1;
	}
	// Get VZERO Information for multiplicity later
	lVV0 = lVevent->GetVZEROData();
	if (!lVV0) {
		AliError("AliVVZERO not available");
		return 9999;
	}

	// Flatenicity calculation
	const Int_t nRings = 4;
	const Int_t nSectors = 8;
	Float_t minEtaV0C[nRings] = {-3.7, -3.2, -2.7, -2.2};
	Float_t maxEtaV0C[nRings] = {-3.2, -2.7, -2.2, -1.7};
	Float_t maxEtaV0A[nRings] = {5.1, 4.5, 3.9, 3.4};
	Float_t minEtaV0A[nRings] = {4.5, 3.9, 3.4, 2.8};
	// Grid
	const Int_t nCells = nRings * 2 * nSectors;
	float RhoLattice[nCells];
	for (Int_t iCh = 0; iCh < nCells; iCh++) {
		RhoLattice[iCh] = 0.0;
	}

	// before calibration
	for (Int_t iCh = 0; iCh < nCells; iCh++) {
		Float_t mult = lVV0->GetMultiplicity(iCh);
		RhoLattice[iCh] = mult;
		hActivityV0DataSectBefore->Fill(iCh, RhoLattice[iCh]);
	}
	// after calibration
	if (fIsCalib) {
		for (Int_t iCh = 0; iCh < nCells; iCh++) {
			RhoLattice[iCh] *= fParVtx->Eval(0.0) / fParVtx->Eval(fVtxz);
		}
	}

	// Filling histos with mult info
	float total_v0_tmp = 0;
	for (Int_t iCh = 0; iCh < nCells; iCh++) {
		hActivityV0DataSect->Fill(iCh, RhoLattice[iCh]);
		total_v0_tmp += RhoLattice[iCh];
	}
	int total_v0 = total_v0_tmp;
	hV0vsVtxz->Fill(fVtxz, total_v0);

	Int_t nringA = 0;
	Int_t nringC = 0;
	for (Int_t iCh = 0; iCh < nCells; iCh++) {
		Float_t detaV0 = -1;
		// Float_t mult = lVV0->GetMultiplicity(iCh);
		if (iCh < 32) { // V0C
			if (iCh < 8) {
				nringC = 0;
			} else if (iCh >= 8 && iCh < 16) {
				nringC = 1;
			} else if (iCh >= 16 && iCh < 24) {
				nringC = 2;
			} else {
				nringC = 3;
			}
			detaV0 = maxEtaV0C[nringC] - minEtaV0C[nringC];
		} else { // V0A
			if (iCh < 40) {
				nringA = 0;
			} else if (iCh >= 40 && iCh < 48) {
				nringA = 1;
			} else if (iCh >= 48 && iCh < 56) {
				nringA = 2;
			} else {
				nringA = 3;
			}
			detaV0 = maxEtaV0A[nringA] - minEtaV0A[nringA];
		}
		// consider the different eta coverage
		RhoLattice[iCh] /= detaV0;
	}
	Float_t mRho = 0;
	Float_t flatenicity = -1;
	for (Int_t iCh = 0; iCh < nCells; iCh++) {
		mRho += RhoLattice[iCh];
	}
	Float_t multiplicityV0M = mRho;
	// average activity per cell
	mRho /= (1.0 * nCells);
	// get sigma
	Double_t sRho_tmp = 0;
	for (Int_t iCh = 0; iCh < nCells; iCh++) {
		sRho_tmp += TMath::Power(1.0 * RhoLattice[iCh] - mRho, 2);
	}
	sRho_tmp /= (1.0 * nCells * nCells);
	Float_t sRho = TMath::Sqrt(sRho_tmp);
	if (mRho > 0) {
		if (fRemoveTrivialScaling) {
			flatenicity = TMath::Sqrt(multiplicityV0M) * sRho / mRho;
		} else {
			flatenicity = sRho / mRho;
		}
	} else {
		flatenicity = 9999;
	}
	return flatenicity;
}

//______________________________________________________________________________
Double_t AliAnalysisTaskFlatenicityMCclosure::GetFlatenicityMC() {

	// Flatenicity calculation
	const Int_t nRings = 8;
	Float_t maxEta[nRings] = {-3.2, -2.7, -2.2, -1.7, 5.1, 4.5, 3.9, 3.4};
	Float_t minEta[nRings] = {-3.7, -3.2, -2.7, -2.2, 4.5, 3.9, 3.4, 2.8};

	const Int_t nSectors = 8;
	Float_t PhiBins[nSectors + 1];
	Float_t deltaPhi = (2.0 * TMath::Pi()) / (1.0 * nSectors);
	for (int i_phi = 0; i_phi < nSectors + 1; ++i_phi) {
		PhiBins[i_phi] = 0;
		if (i_phi < nSectors) {
			PhiBins[i_phi] = i_phi * deltaPhi;
		} else {
			PhiBins[i_phi] = 2.0 * TMath::Pi();
		}
	}

	// Grid
	const Int_t nCells = nRings * nSectors;
	Float_t RhoLattice[nCells];
	for (Int_t iCh = 0; iCh < nCells; iCh++) {
		RhoLattice[iCh] = 0.0;
	}

	Int_t nMult = 0;
	for (Int_t i = 0; i < fMC->GetNumberOfTracks(); ++i) {

		AliMCParticle *particle = (AliMCParticle *)fMC->GetTrack(i);
		if (!particle)
			continue;
		if (!fMC->IsPhysicalPrimary(i))
			continue;
		if (particle->Pt() <= 0.0)
			continue;
		if (TMath::Abs(particle->Charge()) < 0.1)
			continue;
		if (AliAnalysisUtils::IsParticleFromOutOfBunchPileupCollision(i, fMC))
			continue;

		Double_t phi = particle->Phi();
		Double_t eta = particle->Eta();

		Int_t i_segment = 0;
		for (int i_eta = 0; i_eta < nRings; ++i_eta) {

			for (int i_phi = 0; i_phi < nSectors; ++i_phi) {

				if (eta >= minEta[i_eta] && eta < maxEta[i_eta] &&
						phi >= PhiBins[i_phi] && phi < PhiBins[i_phi + 1]) {
					nMult++;
					RhoLattice[i_segment] += 1.0;
				}
				i_segment++;
			}
		}
	}

	Int_t i_segment = 0;
	for (int i_eta = 0; i_eta < nRings; ++i_eta) {
		for (int i_phi = 0; i_phi < nSectors; ++i_phi) {
			Float_t deltaEta = TMath::Abs(maxEta[i_eta] - minEta[i_eta]);
			hActivityV0McSect->Fill(i_segment, RhoLattice[i_segment]);
			RhoLattice[i_segment] /= deltaEta;
			// Filling histos with mult info
			i_segment++;
		}
	}

	Float_t mRho = 0;
	Float_t flatenicity = -1;
	for (Int_t iCh = 0; iCh < nCells; iCh++) {
		mRho += RhoLattice[iCh];
	}
	// average activity per cell
	mRho /= (1.0 * nCells);
	// get sigma
	Float_t sRho_tmp = 0;
	for (Int_t iCh = 0; iCh < nCells; iCh++) {
		sRho_tmp += TMath::Power(1.0 * RhoLattice[iCh] - mRho, 2);
	}
	sRho_tmp /= (1.0 * nCells * nCells);
	Float_t sRho = TMath::Sqrt(sRho_tmp);

	if (mRho > 0) {
		if (fRemoveTrivialScaling) {
			flatenicity = TMath::Sqrt(1.0 * nMult) * sRho / mRho;
		} else {
			flatenicity = sRho / mRho;
		}
	} else {
		sRho = 9999;
	}
	hFlatVsNchMC->Fill(flatenicity, nMult);
	return flatenicity;
}
//______________________________________________________________________________
void AliAnalysisTaskFlatenicityMCclosure::ExtractMultiplicities() {

	fmultTPC = 0;
	Int_t nTracks = fESD->GetNumberOfTracks();
	for (Int_t iT = 0; iT < nTracks; ++iT) {

		AliESDtrack *esdtrack = static_cast<AliESDtrack *>(fESD->GetTrack(iT)); // get a track (type AliesdTrack)
		if (!esdtrack)
			continue;
		if (!fTrackFilter->IsSelected(esdtrack))
			continue;
		if (TMath::Abs(esdtrack->Eta()) > fEtaCut)
			continue;
		if (esdtrack->Pt() < fPtMin)
			continue;

		fmultTPC++;
	}

	AliVVZERO *lVV0 = 0x0;
	AliVEvent *lVevent = 0x0;
	lVevent = dynamic_cast<AliVEvent *>(InputEvent());
	if (!lVevent) {
		AliWarning("ERROR: ESD / AOD event not available \n");
		return;
	}
	// Get VZERO Information for multiplicity later
	lVV0 = lVevent->GetVZEROData();
	if (!lVV0) {
		AliError("AliVVZERO not available");
		return;
	}

	const Int_t nChannels = 64;
	float fmultV0C_tmp = 0;
	float fmultV0A_tmp = 0;
	for (Int_t iCh = 0; iCh < nChannels; iCh++) {
		float mult = lVV0->GetMultiplicity(iCh);
		if (fIsCalib) {
			mult *= fParVtx->Eval(0.0) / fParVtx->Eval(fVtxz);
		}
		if (iCh < 32) { // V0C
			fmultV0C_tmp += mult;
		} else { // V0A
			fmultV0A_tmp += mult;
		}
	}

	fmultV0C = fmultV0C_tmp;
	fmultV0A = fmultV0A_tmp;

	AliVAD *lVAD = 0x0;
	lVAD = lVevent->GetADData();
	if (!lVAD) {
		AliError("AliVAD not available");
		return;
	}
	fmultADA = 0;
	fmultADC = 0;
	// Get Multiplicity info per AD 16 channel: C-side : 0-7, A-side 8-15
	for (Int_t i = 0; i < 8; i++) {
		fmultADA += lVAD->GetMultiplicityADA(i);
	}
	for (Int_t i = 8; i < 16; i++) {
		fmultADC += lVAD->GetMultiplicityADC(i - 8);
	}
}
//______________________________________________________________________________
void AliAnalysisTaskFlatenicityMCclosure::ExtractMultiplicitiesEqualALICE() {

	fmultTPC = 0;
	Int_t nTracks = fESD->GetNumberOfTracks();
	for (Int_t iT = 0; iT < nTracks; ++iT) {

		AliESDtrack *esdtrack = static_cast<AliESDtrack *>(
				fESD->GetTrack(iT)); // get a track (type AliesdTrack)
		if (!esdtrack)
			continue;
		if (!fTrackFilter->IsSelected(esdtrack))
			continue;
		if (TMath::Abs(esdtrack->Eta()) > fEtaCut)
			continue;
		if (esdtrack->Pt() < fPtMin)
			continue;
		fmultTPC++;
	}

	AliVVZERO *lVV0 = 0x0;
	AliVEvent *lVevent = 0x0;
	lVevent = dynamic_cast<AliVEvent *>(InputEvent());
	if (!lVevent) {
		AliWarning("ERROR: ESD / AOD event not available \n");
		return;
	}
	// Get VZERO Information for multiplicity later
	lVV0 = lVevent->GetVZEROData();
	if (!lVV0) {
		AliError("AliVVZERO not available");
		return;
	}
	float fmultV0A_tmp = 0;
	float fmultV0C_tmp = 0;
	fmultV0A_tmp = AliESDUtils::GetCorrV0A(lVV0->GetMTotV0A(), fVtxz);
	fmultV0C_tmp = AliESDUtils::GetCorrV0C(lVV0->GetMTotV0C(), fVtxz);

	fmultV0A = fmultV0A_tmp;
	fmultV0C = fmultV0C_tmp;

	AliVAD *lVAD = 0x0;
	lVAD = lVevent->GetADData();
	if (!lVAD) {
		AliError("AliVAD not available");
		return;
	}
	fmultADA = 0;
	fmultADC = 0;
	// Get Multiplicity info per AD 16 channel: C-side : 0-7, A-side 8-15
	for (Int_t i = 0; i < 8; i++) {
		fmultADA += lVAD->GetMultiplicityADA(i);
	}
	for (Int_t i = 8; i < 16; i++) {
		fmultADC += lVAD->GetMultiplicityADC(i - 8);
	}
}

//______________________________________________________________________________
void AliAnalysisTaskFlatenicityMCclosure::ExtractMultiplicitiesMC() {

	fmultV0Amc = 0;
	fmultV0Cmc = 0;
	fmultADAmc = 0;
	fmultADCmc = 0;
	fmultTPCmc = 0;

	for (Int_t i = 0; i < fMC->GetNumberOfTracks(); ++i) {

		AliMCParticle *particle = (AliMCParticle *)fMC->GetTrack(i);
		if (!particle)
			continue;
		if (!fMC->IsPhysicalPrimary(i))
			continue;
		if (particle->Pt() <= 0.0)
			continue;
		if (TMath::Abs(particle->Charge()) < 0.1)
			continue;
		if (AliAnalysisUtils::IsParticleFromOutOfBunchPileupCollision(i, fMC))
			continue;

		Double_t eta_a = particle->Eta();
		if (eta_a >= 2.8 && eta_a < 5.1) { // v0a acceptance (excluding first ring)
			fmultV0Amc++;
		}
		if (eta_a >= 4.8 && eta_a < 6.3) { // ada acceptance
			fmultADAmc++;
		}
		if (eta_a >= -3.7 && eta_a < -1.7) { // v0c
			fmultV0Cmc++;
		}
		if (eta_a >= -7.0 && eta_a < -4.9) { // adc
			fmultADCmc++;
		}
		if (TMath::Abs(eta_a) < 0.8) { // adc
			fmultTPCmc++;
		}
	}
}
//______________________________________________________________________________
void AliAnalysisTaskFlatenicityMCclosure::MCkinematics() {

	for (Int_t i = 0; i < fMC->GetNumberOfTracks(); ++i) {

		AliMCParticle *particle = (AliMCParticle *)fMC->GetTrack(i);
		if (!particle)
			continue;
		if (!fMC->IsPhysicalPrimary(i))
			continue;
		if (particle->Pt() <= 0.0)
			continue;
		if (TMath::Abs(particle->Charge()) < 0.1)
			continue;
		if (AliAnalysisUtils::IsParticleFromOutOfBunchPileupCollision(i, fMC))
			continue;

		hEtamc->Fill(particle->Eta());
	}
}
//______________________________________________________________________
Bool_t AliAnalysisTaskFlatenicityMCclosure::HasRecVertex() {

	float fMaxDeltaSpdTrackAbsolute = 0.5f;
	float fMaxDeltaSpdTrackNsigmaSPD = 1.e14f;
	float fMaxDeltaSpdTrackNsigmaTrack = 1.e14;
	float fMaxResolutionSPDvertex = 0.25f;
	float fMaxDispersionSPDvertex = 1.e14f;

	Bool_t fRequireTrackVertex = true;
	unsigned long fFlag;
	fFlag = BIT(AliEventCuts::kNoCuts);

	const AliVVertex *vtTrc = fESD->GetPrimaryVertex();
	bool isTrackV = true;
	if (vtTrc->IsFromVertexer3D() || vtTrc->IsFromVertexerZ())
		isTrackV = false;
	const AliVVertex *vtSPD = fESD->GetPrimaryVertexSPD();

	if (vtSPD->GetNContributors() > 0)
		fFlag |= BIT(AliEventCuts::kVertexSPD);

	if (vtTrc->GetNContributors() > 1 && isTrackV)
		fFlag |= BIT(AliEventCuts::kVertexTracks);

	if (((fFlag & BIT(AliEventCuts::kVertexTracks)) || !fRequireTrackVertex) &&
			(fFlag & BIT(AliEventCuts::kVertexSPD)))
		fFlag |= BIT(AliEventCuts::kVertex);

	const AliVVertex *&vtx =
		bool(fFlag & BIT(AliEventCuts::kVertexTracks)) ? vtTrc : vtSPD;
	AliVVertex *fPrimaryVertex = const_cast<AliVVertex *>(vtx);
	if (!fPrimaryVertex)
		return kFALSE;

	/// Vertex quality cuts
	double covTrc[6], covSPD[6];
	vtTrc->GetCovarianceMatrix(covTrc);
	vtSPD->GetCovarianceMatrix(covSPD);
	double dz = bool(fFlag & AliEventCuts::kVertexSPD) &&
		bool(fFlag & AliEventCuts::kVertexTracks)
		? vtTrc->GetZ() - vtSPD->GetZ()
		: 0.; /// If one of the two vertices is not available this cut
		      /// is always passed.
	double errTot = TMath::Sqrt(covTrc[5] + covSPD[5]);
	double errTrc =
		bool(fFlag & AliEventCuts::kVertexTracks) ? TMath::Sqrt(covTrc[5]) : 1.;
	double nsigTot = TMath::Abs(dz) / errTot, nsigTrc = TMath::Abs(dz) / errTrc;
	/// vertex dispersion for run1, only for ESD, AOD code to be added here
	const AliESDVertex *vtSPDESD = dynamic_cast<const AliESDVertex *>(vtSPD);
	double vtSPDdispersion = vtSPDESD ? vtSPDESD->GetDispersion() : 0;
	if ((TMath::Abs(dz) <= fMaxDeltaSpdTrackAbsolute &&
				nsigTot <= fMaxDeltaSpdTrackNsigmaSPD &&
				nsigTrc <=
				fMaxDeltaSpdTrackNsigmaTrack) && // discrepancy track-SPD vertex
			(!vtSPD->IsFromVertexerZ() ||
			 TMath::Sqrt(covSPD[5]) <= fMaxResolutionSPDvertex) &&
			(!vtSPD->IsFromVertexerZ() ||
			 vtSPDdispersion <= fMaxDispersionSPDvertex) /// vertex dispersion cut for
								     /// run1, only for ESD
	   ) // quality cut on vertexer SPD z
		fFlag |= BIT(AliEventCuts::kVertexQuality);

	Bool_t hasVtx = (TESTBIT(fFlag, AliEventCuts::kVertex)) &&
		(TESTBIT(fFlag, AliEventCuts::kVertexQuality));
	fVtxz = vtSPD->GetZ();
	return hasVtx;
}
//______________________________________________________________________________
void AliAnalysisTaskFlatenicityMCclosure::GetMCchargedTrueDists() {

	for (Int_t i = 0; i < fMC->GetNumberOfTracks(); ++i) {
		AliMCParticle *particle = (AliMCParticle *)fMC->GetTrack(i);
		if (!particle)
			continue;
		if (!fMC->IsPhysicalPrimary(i))
			continue;
		if (TMath::Abs(particle->Eta()) > fEtaCut)
			continue;
		if (particle->Pt() < fPtMin)
			continue;
		if (TMath::Abs(particle->Charge()) < 0.1)
			continue;
		if (AliAnalysisUtils::IsParticleFromOutOfBunchPileupCollision(i, fMC))
			continue;

		int pdgCode = particle->PdgCode();
		int pidCode = GetPidCode(pdgCode);

		hMultV0AV0CvsFlatvspT_BFTrigSel->Fill(fmultV0Amc + fmultV0Cmc,fFlatMC,particle->Pt());
		hAmpV0AV0CvsFlatvspT_BFTrigSel->Fill(fv0mamplitude,fFlatMC,particle->Pt());
		hPercentileV0MvsFlatvspT_BFTrigSel->Fill(fv0mpercentile,fFlatMC,particle->Pt());

		if (pidCode==1) { hMultV0AV0CvsFlatvspT_pi_BFTrigSel->Fill(fmultV0Amc + fmultV0Cmc,fFlatMC,particle->Pt()); }
		if (pidCode==2) { hMultV0AV0CvsFlatvspT_k_BFTrigSel->Fill(fmultV0Amc + fmultV0Cmc,fFlatMC,particle->Pt()); }
		if (pidCode==3) { hMultV0AV0CvsFlatvspT_p_BFTrigSel->Fill(fmultV0Amc + fmultV0Cmc,fFlatMC,particle->Pt()); }
	}
}
//______________________________________________________________________________
int AliAnalysisTaskFlatenicityMCclosure::GetPidCode(int pdgCode) 
{
	// return our internal code for pions, kaons, and protons

	Int_t pidCode = 6;

	switch (TMath::Abs(pdgCode)) {
		case 211:
			pidCode = 1; // pion
			break;
		case 321:
			pidCode = 2; // kaon
			break;
		case 2212:
			pidCode = 3; // proton
			break;
		case 310:
			pidCode = 4; // K0s
			break;
		case 3122:
			pidCode = 5; // lambda
			break;
		case 333:
			pidCode = 6; // phi
			break;
		case 11:
			pidCode = 7; // electron
			break;
		default:
			pidCode = 8;  // something else?
	};

	return pidCode;
}
