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

class AliPPVsMultUtils;
class AliESDtrackCuts;

#include "AliAnalysisManager.h"
#include "AliAnalysisTask.h"
#include "AliAnalysisTaskESDfilter.h"
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisUtils.h"
#include "AliCentrality.h"
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
#include "AliMultEstimator.h"
#include "AliMultInput.h"
#include "AliMultSelection.h"
#include "AliMultVariable.h"
#include "AliMultiplicity.h"
#include "AliOADBContainer.h"
#include "AliOADBMultSelection.h"
#include "AliPPVsMultUtils.h"
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

#include "AliAnalysisTaskFlatenicityPiKp.h"

static const Int_t nCent = 9;
static const Int_t nEta = 4;

static Double_t centClass[nCent + 1] = {0.0,  1.0,  5.0,  10.0, 20.0,
	30.0, 40.0, 50.0, 70.0, 100.0};

static const Char_t* V0MClass[nCent] = {"0_1","1_5","5_10","10_20",
	"20_30","30_40","40_50","50_70","70_100"};

static const Char_t* etaClass[nEta] = {"02","24","46","68"};
static const Char_t* ParticleType[3] = {"Primaries","MaterialInt","WeakDecays"};

static const double C_Value = TMath::C()*(1.e2/1.e12); // cm/ps

using namespace std; // std namespace: so you can do things like 'cout' etc

ClassImp(AliAnalysisTaskFlatenicityPiKp) // classimp: necessary for root

AliAnalysisTaskFlatenicityPiKp::AliAnalysisTaskFlatenicityPiKp()
	: AliAnalysisTaskSE(), fESD(0), fEventCuts(0x0), fMCStack(0), fMC(0),
	fUseMC(kFALSE), fV0Mindex(-1), fV0MMultiplicity(-1.0), fDetFlat("V0"), fIsMCclosure(kFALSE),
	fDeltaV0(kTRUE), fRemoveTrivialScaling(kFALSE), fnGen(-1), fPIDResponse(0x0),
	fTrackFilter(0x0), fTrackFilterPID(0x0), fOutputList(0), fEtaCut(0.8), fPtMin(0.5), fNcl(70), fV0MEqualisation(kTRUE) ,fdEdxCalibrated(kTRUE),
	fSaveDCAxyHistograms(kFALSE), 
	fEtaCalibrationPos(0x0), fEtaCalibrationNeg(0x0), fV0CCalibration(0x0), fV0ACalibration(0x0), 
	fcutLow(0x0), fcutHigh(0x0), fcutDCAxy(0x0), fPeriod("16l"),
	ftrackmult08(0), fv0mpercentile(0), fMidRapidityMult(0), fFlat(-1), fFlatTPC(-1.), fFlatMC(-1),
	fMultSelection(0x0), hPtPrimIn(0), hPtPrimOut(0), hPtSecOut(0), hPtOut(0),
	hFlatenicityMC(0), hFlatResponse(0), hFlatVsPtMC(0), hActivityV0CV0A(0),
	hActivityV0DataSect(0), hActivityV0McSect(0), hFlatVsNchMC(0),
	hFlatVsV0MVsMult(0),
	hMCPtPionPos(0),hMCPtKaonPos(0),hMCPtProtonPos(0),
	hMCPtPionNeg(0),hMCPtKaonNeg(0),hMCPtProtonNeg(0),
	hTPCRecTracksPionPos(0), hTPCRecTracksKaonPos(0), hTPCRecTracksProtonPos(0),
	hTPCRecTracksPionNeg(0), hTPCRecTracksKaonNeg(0), hTPCRecTracksProtonNeg(0),
	hTOFRecTracksPionPos(0), hTOFRecTracksKaonPos(0), hTOFRecTracksProtonPos(0),
	hTOFRecTracksPionNeg(0), hTOFRecTracksKaonNeg(0), hTOFRecTracksProtonNeg(0),
	hrTPCRecTracksPion(0), hrTPCRecTracksKaon(0), hrTPCRecTracksProton(0),
	hPionTPCDCAxyNegData(0), hPionTPCDCAxyPosData(0), hProtonTPCDCAxyNegData(0), hProtonTPCDCAxyPosData(0),
	hPionTOFDCAxyNegData(0), hPionTOFDCAxyPosData(0), hProtonTOFDCAxyNegData(0), hProtonTOFDCAxyPosData(0),
	hMIPVsEta(0), pMIPVsEta(0), hPlateauVsEta(0), pPlateauVsEta(0)


{
	for (Int_t i_c = 0; i_c < nCent; ++i_c) {
		hFlatVsPtV0M[i_c] = 0;

		for (Int_t i_eta = 0; i_eta < nEta; ++i_eta){

			hNsigmaPiPos[i_c][i_eta] = 0;
			hNsigmaKPos[i_c][i_eta] = 0;
			hNsigmaPPos[i_c][i_eta] = 0;
			hPtTPCEtaPos[i_c][i_eta] = 0;

			hNsigmaPiNeg[i_c][i_eta] = 0;
			hNsigmaKNeg[i_c][i_eta] = 0;
			hNsigmaPNeg[i_c][i_eta] = 0;
			hPtTPCEtaNeg[i_c][i_eta] = 0;

			hBetaPos[i_c][i_eta] = 0;
			hMomentumTOFEtaPos[i_c][i_eta] = 0;
			hPtTOFEtaPos[i_c][i_eta] = 0;

			hBetaNeg[i_c][i_eta] = 0;
			hMomentumTOFEtaNeg[i_c][i_eta] = 0;
			hPtTOFEtaNeg[i_c][i_eta] = 0;

			hdEdx[i_c][i_eta] = 0;
			hPtrTPC[i_c][i_eta] = 0;
			if (i_c==0){
				hPtVsP[i_eta] = 0;	
				pion_cont_in_kaon_h[i_eta] = 0;
				electron_cont_in_kaon_h[i_eta] = 0;
				nsigma_kaon_h[i_eta] = 0;
				pion_cont_in_proton_h[i_eta] = 0;
				electron_cont_in_proton_h[i_eta] = 0;
				nsigma_proton_h[i_eta] = 0;
				kaon_cont_in_pion_h[i_eta] = 0;
				electron_cont_in_pion_h[i_eta] = 0;
				nsigma_pion_h[i_eta] = 0;
			}
		}
	}

	for(Int_t i = 0; i < 3; ++i){
		hPionTOFDCAxyNeg[i] = 0;
		hProtonTOFDCAxyNeg[i] = 0;
		hPionTOFDCAxyPos[i] = 0;
		hProtonTOFDCAxyPos[i] = 0;
		hPionTPCDCAxyNeg[i] = 0;
		hProtonTPCDCAxyNeg[i] = 0;	
		hPionTPCDCAxyPos[i] = 0;	
		hProtonTPCDCAxyPos[i] = 0;

	}

}
//_____________________________________________________________________________
AliAnalysisTaskFlatenicityPiKp::AliAnalysisTaskFlatenicityPiKp(const char *name)
	: AliAnalysisTaskSE(name), fESD(0), fEventCuts(0x0), fMCStack(0), fMC(0),
	fUseMC(kFALSE), fV0Mindex(-1), fV0MMultiplicity(-1.0), fDetFlat("V0"), fIsMCclosure(kFALSE),
	fDeltaV0(kTRUE), fRemoveTrivialScaling(kFALSE), fnGen(-1), fPIDResponse(0x0),
	fTrackFilter(0x0), fTrackFilterPID(0x0), fOutputList(0), fEtaCut(0.8), fPtMin(0.5), fNcl(70), fV0MEqualisation(kTRUE), fdEdxCalibrated(kTRUE), 
	fSaveDCAxyHistograms(kFALSE), 
	fEtaCalibrationPos(0x0), fEtaCalibrationNeg(0x0), fV0CCalibration(0x0), fV0ACalibration(0x0), 
	fcutLow(0x0), fcutHigh(0x0), fcutDCAxy(0x0), fPeriod("16l"),
	ftrackmult08(0), fv0mpercentile(0), fMidRapidityMult(0), fFlat(-1), fFlatTPC(-1.), fFlatMC(-1),
	fMultSelection(0x0), hPtPrimIn(0), hPtPrimOut(0), hPtSecOut(0), hPtOut(0),
	hFlatenicityMC(0), hFlatResponse(0), hFlatVsPtMC(0), hActivityV0CV0A(0),
	hActivityV0DataSect(0), hActivityV0McSect(0), hFlatVsNchMC(0),
	hFlatVsV0MVsMult(0), 
	hMCPtPionPos(0),hMCPtKaonPos(0),hMCPtProtonPos(0),
	hMCPtPionNeg(0),hMCPtKaonNeg(0),hMCPtProtonNeg(0),
	hTPCRecTracksPionPos(0), hTPCRecTracksKaonPos(0), hTPCRecTracksProtonPos(0),
	hTPCRecTracksPionNeg(0), hTPCRecTracksKaonNeg(0), hTPCRecTracksProtonNeg(0),
	hTOFRecTracksPionPos(0), hTOFRecTracksKaonPos(0), hTOFRecTracksProtonPos(0),
	hTOFRecTracksPionNeg(0), hTOFRecTracksKaonNeg(0), hTOFRecTracksProtonNeg(0),
	hrTPCRecTracksPion(0), hrTPCRecTracksKaon(0), hrTPCRecTracksProton(0),
	hPionTPCDCAxyNegData(0), hPionTPCDCAxyPosData(0), hProtonTPCDCAxyNegData(0), hProtonTPCDCAxyPosData(0),
	hPionTOFDCAxyNegData(0), hPionTOFDCAxyPosData(0), hProtonTOFDCAxyNegData(0), hProtonTOFDCAxyPosData(0),
	hMIPVsEta(0), pMIPVsEta(0), hPlateauVsEta(0), pPlateauVsEta(0)

{
	for (Int_t i_c = 0; i_c < nCent; ++i_c) {
		hFlatVsPtV0M[i_c] = 0;

		for (Int_t i_eta = 0; i_eta < nEta; ++i_eta){

			hNsigmaPiPos[i_c][i_eta] = 0;
			hNsigmaKPos[i_c][i_eta] = 0;
			hNsigmaPPos[i_c][i_eta] = 0;
			hPtTPCEtaPos[i_c][i_eta] = 0;

			hNsigmaPiNeg[i_c][i_eta] = 0;
			hNsigmaKNeg[i_c][i_eta] = 0;
			hNsigmaPNeg[i_c][i_eta] = 0;
			hPtTPCEtaNeg[i_c][i_eta] = 0;

			hBetaPos[i_c][i_eta] = 0;
			hMomentumTOFEtaPos[i_c][i_eta] = 0;
			hPtTOFEtaPos[i_c][i_eta] = 0;

			hBetaNeg[i_c][i_eta] = 0;
			hMomentumTOFEtaNeg[i_c][i_eta] = 0;
			hPtTOFEtaNeg[i_c][i_eta] = 0;

			hdEdx[i_c][i_eta] = 0;
			hPtrTPC[i_c][i_eta] = 0;
			if (i_c==0){
				hPtVsP[i_eta] = 0;	
				pion_cont_in_kaon_h[i_eta] = 0;
				electron_cont_in_kaon_h[i_eta] = 0;
				nsigma_kaon_h[i_eta] = 0;
				pion_cont_in_proton_h[i_eta] = 0;
				electron_cont_in_proton_h[i_eta] = 0;
				nsigma_proton_h[i_eta] = 0;
				kaon_cont_in_pion_h[i_eta] = 0;
				electron_cont_in_pion_h[i_eta] = 0;
				nsigma_pion_h[i_eta] = 0;
			}
		}
	}

	for(Int_t i = 0; i < 3; ++i){
		hPionTOFDCAxyNeg[i] = 0;
		hProtonTOFDCAxyNeg[i] = 0;
		hPionTOFDCAxyPos[i] = 0;
		hProtonTOFDCAxyPos[i] = 0;
		hPionTPCDCAxyNeg[i] = 0;
		hProtonTPCDCAxyNeg[i] = 0;	
		hPionTPCDCAxyPos[i] = 0;	
		hProtonTPCDCAxyPos[i] = 0;

	}

	DefineInput(0, TChain::Class()); // define the input of the analysis: in this
					 // case you take a 'chain' of events
	DefineOutput(1, TList::Class()); // define the ouptut of the analysis: in this
					 // case it's a list of histograms
}
//_____________________________________________________________________________
AliAnalysisTaskFlatenicityPiKp::~AliAnalysisTaskFlatenicityPiKp() {
	// destructor
	if (fOutputList) {
		delete fOutputList; // at the end of your task, it is deleted from memory by
				    // calling this function
		fOutputList = 0x0;
	}
}
//_____________________________________________________________________________
void AliAnalysisTaskFlatenicityPiKp::UserCreateOutputObjects() {

	AliAnalysisManager *man = AliAnalysisManager::GetAnalysisManager();
	if(man){
		AliInputEventHandler* inputHandler = (AliInputEventHandler*)(man->GetInputEventHandler());
		if(inputHandler)fPIDResponse = inputHandler->GetPIDResponse();
	}

	// create track filters
	fTrackFilter = new AliAnalysisFilter("trackFilter");
	AliESDtrackCuts *fCuts = new AliESDtrackCuts();
	fCuts->SetAcceptKinkDaughters(kFALSE);
	fCuts->SetRequireTPCRefit(kTRUE);
	fCuts->SetRequireITSRefit(kTRUE);
	fCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kAny);
	fCuts->SetDCAToVertex2D(kFALSE);
	fCuts->SetRequireSigmaToVertex(kFALSE);
	fCuts->SetEtaRange(-0.8, 0.8);
	fCuts->SetMinNCrossedRowsTPC(70);
	fCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);
	fCuts->SetMaxChi2PerClusterTPC(4);
	fCuts->SetMaxDCAToVertexZ(2);
	fCuts->SetCutGeoNcrNcl(3., 130., 1.5, 0.85, 0.7);
	fCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);
	fCuts->SetMaxChi2PerClusterTPC(4);
	fCuts->SetMaxDCAToVertexZ(2);
	fCuts->SetMaxChi2PerClusterITS(36);
	fCuts->SetMaxDCAToVertexXYPtDep("0.0105+0.0350/pt^1.1");
	fCuts->SetMaxChi2PerClusterITS(36);
	fTrackFilter->AddCuts(fCuts);

	fTrackFilterPID = new AliAnalysisFilter("trackFilterPID");
	AliESDtrackCuts* fCutsPID = new AliESDtrackCuts;
	fCutsPID->SetMinNCrossedRowsTPC(70);
	fCutsPID->SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);
	fCutsPID->SetMaxChi2PerClusterTPC(4);
	fCutsPID->SetAcceptKinkDaughters(kFALSE);
	fCutsPID->SetRequireTPCRefit(kTRUE);
	fCutsPID->SetRequireITSRefit(kTRUE);
	fCutsPID->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kAny);
	fCutsPID->SetMaxDCAToVertexZ(2);
	fCutsPID->SetDCAToVertex2D(kFALSE);
	fCutsPID->SetRequireSigmaToVertex(kFALSE);
	fCutsPID->SetMaxChi2PerClusterITS(36);
	fTrackFilterPID->AddCuts(fCutsPID);

	fEtaCalibrationPos = new TF1("fDeDxVsEtaPos", "pol7", 0.0, 1.0);
	fEtaCalibrationNeg = new TF1("fDeDxVsEtaNeg", "pol7", -1.0, 0.0);

	fV0CCalibration = new TF1("fV0CCalibration","(x>-0.5 && x<0.5)*[0]+(x>0.5 && x<1.5)*[1]+(x>1.5 && x<2.5)*[2]+(x>2.5 && x<3.5)*[3]+(x>3.5 && x<4.5)*[4]+(x>4.5 && x<5.5)*[5]+(x>5.5 && x<6.5)*[6]+(x>6.5 && x<7.5)*[7]+(x>7.5 && x<8.5)*[8]+(x>8.5 && x<9.5)*[9]+(x>9.5 && x<10.5)*[10]+(x>10.5 && x<11.5)*[11]+(x>11.5 && x<12.5)*[12]+(x>12.5 && x<13.5)*[13]+(x>13.5 && x<14.5)*[14]+(x>14.5 && x<15.5)*[15]+(x>15.5 && x<16.5)*[16]+(x>16.5 && x<17.5)*[17]+(x>17.5 && x<18.5)*[18]+(x>18.5 && x<19.5)*[19]+(x>19.5 && x<20.5)*[20]+(x>20.5 && x<21.5)*[21]+(x>21.5 && x<22.5)*[22]+(x>22.5 && x<23.5)*[23]+(x>23.5 && x<24.5)*[24]+(x>24.5 && x<25.5)*[25]+(x>25.5 && x<26.5)*[26]+(x>26.5 && x<27.5)*[27]+(x>27.5 && x<28.5)*[28]+(x>28.5 && x<29.5)*[29]+(x>29.5 && x<30.5)*[30]+(x>30.5 && x<31.5)*[31]",-1.0,64.0); 

	fV0ACalibration = new TF1("fV0ACalibration","(x>31.5 && x<32.5)*[0]+(x>32.5 && x<33.5)*[1]+(x>33.5 && x<34.5)*[2]+(x>34.5 && x<35.5)*[3]+(x>35.5 && x<36.5)*[4]+(x>36.5 && x<37.5)*[5]+(x>37.5 && x<38.5)*[6]+(x>38.5 && x<39.5)*[7]+(x>39.5 && x<40.5)*[8]+(x>40.5 && x<41.5)*[9]+(x>41.5 && x<42.5)*[10]+(x>42.5 && x<43.5)*[11]+(x>43.5 && x<44.5)*[12]+(x>44.5 && x<45.5)*[13]+(x>45.5 && x<46.5)*[14]+(x>46.5 && x<47.5)*[15]+(x>47.5 && x<48.5)*[16]+(x>48.5 && x<49.5)*[17]+(x>49.5 && x<50.5)*[18]+(x>50.5 && x<51.5)*[19]+(x>51.5 && x<52.5)*[20]+(x>52.5 && x<53.5)*[21]+(x>53.5 && x<54.5)*[22]+(x>54.5 && x<55.5)*[23]+(x>55.5 && x<56.5)*[24]+(x>56.5 && x<57.5)*[25]+(x>57.5 && x<58.5)*[26]+(x>58.5 && x<59.5)*[27]+(x>59.5 && x<60.5)*[28]+(x>60.5 && x<61.5)*[29]+(x>61.5 && x<62.5)*[30]+(x>62.5 && x<63.5)*[31]",-1.0,64.0); 

	const int nParsV0 = 32;
	const double ParsV0C[nParsV0] = { 1.647294, 1.837619, 1.761105, 1.555255, 1.726714, 1.303948, 1.687114, 1.878224, 1.546167, 2.003266, 1.712996, 1.549461, 1.727268, 1.554480, 1.706153, 1.363934, 1.872374, 1.727050, 2.199335, 0.787233, 1.562607, 1.868105, 2.120761, 1.903412, 2.133782, 2.181180, 1.799161, 1.951639, 1.703616, 1.711802, 1.762529, 1.704079 };

	const double ParsV0A[nParsV0] = { 0.729977, 0.700893, 0.763448, 0.743796, 0.669337, 0.694408, 0.620627, 0.525617, 1.035450, 0.655640, 1.137670, 0.923248, 0.859047, 0.781449, 0.398952, 0.835587, 1.110981, 1.328355, 1.060960, 1.338867, 1.145839, 1.169695, 1.125137, 1.252237, 0.894590, 1.053655, 1.171697, 0.986297, 1.192937, 1.107533, 1.018933, 1.009477 };

	for (Int_t i = 1; i <= nParsV0; ++i){
		fV0CCalibration->SetParameter(i-1,ParsV0C[i-1]);
		fV0ACalibration->SetParameter(i-1,ParsV0A[i-1]);
	}

	fcutLow = new TF1("StandardPhiCutLow",  "0.1/x/x+pi/18.0-0.025", 0, 50);
	fcutHigh = new TF1("StandardPhiCutHigh", "0.12/x+pi/18.0+0.035", 0, 50);

	fcutDCAxy = new TF1("fMaxDCAxy","[0]+[1]/(x^[2])",0,1e10);
	fcutDCAxy->SetParameter(0,0.0105);
	fcutDCAxy->SetParameter(1,0.0350);
	fcutDCAxy->SetParameter(2,1.1);

	const int nPtbins = 56;
	double Ptbins[nPtbins+1] = {
		0.25, 0.30, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 
		0.8, 0.85, 0.9, 0.95, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7,
		1.8, 1.9, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4,
		3.6, 3.8, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 8.0,
		9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 18.0,
		20.0,22.0,24.0,26.0,30.0};

	/*
	   const Int_t nPtbins = 31;
	   Double_t Ptbins[nPtbins+1] = {
	   0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.2, 
	   1.4, 1.6, 1.8, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 
	   5.5, 6.0, 6.5, 7.0, 8.0, 9.0, 10.0, 12.0, 14.0, 16.0, 
	   18.0, 20.0 };
	   */

	// create output objects
	float min_flat = -0.01;
	float max_flat = 1.01;
	int nbins_flat = 1020;
	if (fRemoveTrivialScaling) {
		min_flat = -0.1;
		max_flat = 9.9;
		nbins_flat = 2000;
	}

	const int nFlatbins = 204;
	double Flatbins[nFlatbins+1] = {0.0};
	for (int i = 0; i <= nFlatbins; ++i) {
		Flatbins[i] = -0.01 + (double)i * 0.005;
	}

	const int nMultbins = 100;
	double Multbins[nMultbins+1] = {0.0};
	for (int i = 0; i <= nMultbins; ++i) {
		Multbins[i] = -0.5 + (double)i;
	}

	const int nnSigmabins = 80;
	double nSigmabins[nnSigmabins+1] = {0.0};

	for (int i = 0; i <= nnSigmabins; ++i){
		nSigmabins[i] = -4.0 + i*0.1;
	}

	const int nBetabins   = 151;
	double Betabins[nBetabins+1] = { 0.0 };
	for(int i = 0; i <= nBetabins; ++i){
		Betabins[i] = 0.6 + 0.003333 * ((double)i);
	}

	const int dEdxHigh = 100;
	const int dEdxLow = 40;

	const int ndEdxbins = dEdxHigh-dEdxLow;
	double dEdxbins[ndEdxbins+1];

	for(int i = dEdxLow; i <= dEdxHigh; ++i){
		dEdxbins[i-dEdxLow] = i;
	}

	const int nV0Multbins = 51;
	double V0Multbins[nV0Multbins+1] = { 0.0 };

	for (Int_t i = 0; i <= nV0Multbins; ++i){
		V0Multbins[i] = -0.5 + (double)i;
	}

	const int nV0Sectorsbins = 64;
	double V0Sectorsbins[nV0Sectorsbins+1] = { 0.0 };

	for (int i = 0; i <= nV0Sectorsbins; ++i){
		V0Sectorsbins[i] = -0.5 + (double)i;
	}

	OpenFile(1);
	fOutputList =
		new TList(); // this is a list which will contain all of your histograms
	fOutputList->SetOwner(kTRUE); // memory stuff: the list is owner of all

	hPionTPCDCAxyNegData = new TH2F("hPionTPCDCAxyNeg","; #it{p}_{T} (GeV/#it{c}); DCA_{xy}", nPtbins, Ptbins, 1000, -3.5, 3.5 );	
	hPionTPCDCAxyPosData = new TH2F("hPionTPCDCAxyPos","; #it{p}_{T} (GeV/#it{c}); DCA_{xy}", nPtbins, Ptbins, 1000, -3.5, 3.5 );	
	hProtonTPCDCAxyNegData = new TH2F("hProtonTPCDCAxyNeg","; #it{p}_{T} (GeV/#it{c}); DCA_{xy}", nPtbins, Ptbins, 1000, -3.5, 3.5 );	
	hProtonTPCDCAxyPosData = new TH2F("hProtonTPCDCAxyPos","; #it{p}_{T} (GeV/#it{c}); DCA_{xy}", nPtbins, Ptbins, 1000, -3.5, 3.5 );	

	hPionTOFDCAxyNegData = new TH2F("hPionTOFDCAxyNeg","; #it{p}_{T} (GeV/#it{c}); DCA_{xy}", nPtbins, Ptbins, 1000, -3.5, 3.5 );	
	hPionTOFDCAxyPosData = new TH2F("hPionTOFDCAxyPos","; #it{p}_{T} (GeV/#it{c}); DCA_{xy}", nPtbins, Ptbins, 1000, -3.5, 3.5 );	
	hProtonTOFDCAxyNegData = new TH2F("hProtonTOFDCAxyNeg","; #it{p}_{T} (GeV/#it{c}); DCA_{xy}", nPtbins, Ptbins, 1000, -3.5, 3.5 );	
	hProtonTOFDCAxyPosData = new TH2F("hProtonTOFDCAxyPos","; #it{p}_{T} (GeV/#it{c}); DCA_{xy}", nPtbins, Ptbins, 1000, -3.5, 3.5 );	

	hMIPVsEta = new TH2F("hMIPVsEta","; #eta; dE/dx_{MIP, primary tracks}", 50, -0.8, 0.8, 60-40, 40.0, 60.0);
	pMIPVsEta = new TProfile("pMIPVsEta","; #eta; #LT dE/dx #GT_{MIP, primary tracks}", 50, -0.8, 0.8, 40.0, 60.0);
	hPlateauVsEta = new TH2F("hPlateauVsEta","; #eta; dE/dx_{Plateau, primary tracks}",50, -0.8, 0.8, 50, 60.0, 110.0);
	pPlateauVsEta = new TProfile("pPlateauVsEta","; #eta; #LT dE/dx #GT_{Plateau, primary tracks}", 50, -0.8, 0.8, 60.0, 110.0);

	for (int i_c = 0; i_c < nCent; ++i_c) 
	{
		hFlatVsPtV0M[i_c] = new TH2F(Form("hPtVsFlatV0M_c%d", i_c),"; Flat. V0M; #it{p}_{T} (GeV/#it{c})", nFlatbins, Flatbins, nPtbins, Ptbins);
		//fOutputList->Add(hFlatVsPtV0M[i_c]);

		for (Int_t i_eta = 0; i_eta < nEta; ++i_eta) {

			hNsigmaPiPos[i_c][i_eta] = new TH3F(Form("hNsigmaPiPos_c%s_eta_%s",V0MClass[i_c],etaClass[i_eta]),";#it{p}_{T} (GeV/#it{c}); n#sigma; Flatenicity", nPtbins, Ptbins, nnSigmabins, nSigmabins, nFlatbins, Flatbins);
			hNsigmaKPos[i_c][i_eta] = new TH3F(Form("hNsigmaKPos_c%s_eta_%s",V0MClass[i_c],etaClass[i_eta]),";#it{p}_{T} (GeV/#it{c}); n#sigma; Flatenicity", nPtbins, Ptbins, nnSigmabins, nSigmabins, nFlatbins, Flatbins);
			hNsigmaPPos[i_c][i_eta] = new TH3F(Form("hNsigmaPPos_c%s_eta_%s",V0MClass[i_c],etaClass[i_eta]),";#it{p}_{T} (GeV/#it{c}); n#sigma; Flatenicity", nPtbins, Ptbins, nnSigmabins, nSigmabins, nFlatbins, Flatbins);
			hPtTPCEtaPos[i_c][i_eta] = new TH2F(Form("hPtTPCEtaPos_c%s_eta_%s",V0MClass[i_c],etaClass[i_eta]),"; #it{p}_{T} (GeV/#it{c}; Flatenicity;)", nPtbins, Ptbins, nFlatbins, Flatbins);

			hNsigmaPiNeg[i_c][i_eta] = new TH3F(Form("hNsigmaPiNeg_c%s_eta_%s",V0MClass[i_c],etaClass[i_eta]),";#it{p}_{T} (GeV/#it{c}); n#sigma; Flatenicity", nPtbins, Ptbins, nnSigmabins, nSigmabins, nFlatbins, Flatbins);
			hNsigmaKNeg[i_c][i_eta] = new TH3F(Form("hNsigmaKNeg_c%s_eta_%s",V0MClass[i_c],etaClass[i_eta]),";#it{p}_{T} (GeV/#it{c}); n#sigma; Flatenicity", nPtbins, Ptbins, nnSigmabins, nSigmabins, nFlatbins, Flatbins);
			hNsigmaPNeg[i_c][i_eta] = new TH3F(Form("hNsigmaPNeg_c%s_eta_%s",V0MClass[i_c],etaClass[i_eta]),";#it{p}_{T} (GeV/#it{c}); n#sigma; Flatenicity", nPtbins, Ptbins, nnSigmabins, nSigmabins, nFlatbins, Flatbins);
			hPtTPCEtaNeg[i_c][i_eta] = new TH2F(Form("hPtTPCEtaNeg_c%s_eta_%s",V0MClass[i_c],etaClass[i_eta]),"; #it{p}_{T} (GeV/#it{c}); Flatenicity;)", nPtbins, Ptbins, nFlatbins, Flatbins);

			hBetaPos[i_c][i_eta] = new TH3F(Form("hBetaPos_c%s_eta_%s",V0MClass[i_c],etaClass[i_eta]),";#it{p} (GeV/#it{c}); #beta; Flatenicity", nPtbins, Ptbins, nBetabins,Betabins, nFlatbins, Flatbins);
			hMomentumTOFEtaPos[i_c][i_eta] = new TH2F(Form("hMomentumTOFEtaPos_c%s_eta_%s",V0MClass[i_c],etaClass[i_eta]),";#it{p} (GeV/#it{c}); Flatenicity", nPtbins, Ptbins, nFlatbins, Flatbins);
			hPtTOFEtaPos[i_c][i_eta] = new TH2F(Form("hPtTOFEtaPos_c%s_eta_%s",V0MClass[i_c],etaClass[i_eta]),";#it{p}_{T} (GeV/#it{c}); Flatenicity", nPtbins, Ptbins, nFlatbins, Flatbins);

			hBetaNeg[i_c][i_eta] = new TH3F(Form("hBetaNeg_c%s_eta_%s",V0MClass[i_c],etaClass[i_eta]),";#it{p}; #beta; Flatenicity", nPtbins, Ptbins, nBetabins,Betabins, nFlatbins, Flatbins);
			hMomentumTOFEtaNeg[i_c][i_eta] = new TH2F(Form("hMomentumTOFEtaNeg_c%s_eta_%s",V0MClass[i_c],etaClass[i_eta]),";#it{p} (GeV/#it{c}); Flatenicity", nPtbins, Ptbins, nFlatbins, Flatbins);
			hPtTOFEtaNeg[i_c][i_eta] = new TH2F(Form("hPtTOFEtaNeg_c%s_eta_%s",V0MClass[i_c],etaClass[i_eta]),";#it{p}_{T} (GeV/#it{c}); Flatenicity", nPtbins, Ptbins, nFlatbins, Flatbins);

			hdEdx[i_c][i_eta] = new TH3F(Form("hdEdx_c%s_eta_%s",V0MClass[i_c],etaClass[i_eta]), ";#it{p} (GeV/#it{c}); dE/dx; Flatenicity", nPtbins, Ptbins, ndEdxbins, dEdxbins, nFlatbins, Flatbins);
			hPtrTPC[i_c][i_eta] = new TH2F(Form("hPtrTPC_c%s_eta_%s",V0MClass[i_c],etaClass[i_eta]),";#it{p}_{T} (GeV/#it{c}); Flatenicity", nPtbins, Ptbins, nFlatbins, Flatbins);

			if (i_c==0){
				hPtVsP[i_eta] = new TH2F(Form("hPtVsP_eta_%s",etaClass[i_eta]), ";#it{p} (GeV/#it{c}); #it{p}_{T} (GeV/#it{c})", nPtbins, Ptbins, nPtbins, Ptbins);
			}

			if (!fUseMC && fV0MEqualisation){ 

				fOutputList->Add(hNsigmaPiPos[i_c][i_eta]);
				fOutputList->Add(hNsigmaKPos[i_c][i_eta]);
				fOutputList->Add(hNsigmaPPos[i_c][i_eta]);
				fOutputList->Add(hPtTPCEtaPos[i_c][i_eta]);
				fOutputList->Add(hNsigmaPiNeg[i_c][i_eta]);
				fOutputList->Add(hNsigmaKNeg[i_c][i_eta]);
				fOutputList->Add(hNsigmaPNeg[i_c][i_eta]);
				fOutputList->Add(hPtTPCEtaNeg[i_c][i_eta]);

				fOutputList->Add(hBetaPos[i_c][i_eta]);
				fOutputList->Add(hMomentumTOFEtaPos[i_c][i_eta]);
				fOutputList->Add(hPtTOFEtaPos[i_c][i_eta]);
				fOutputList->Add(hBetaNeg[i_c][i_eta]);
				fOutputList->Add(hMomentumTOFEtaNeg[i_c][i_eta]);
				fOutputList->Add(hPtTOFEtaNeg[i_c][i_eta]);
				fOutputList->Add(hdEdx[i_c][i_eta]);
				fOutputList->Add(hPtrTPC[i_c][i_eta]);

				fOutputList->Add(hPtVsP[i_eta]);

				if (i_eta == 0){
					if (fSaveDCAxyHistograms) {
						fOutputList->Add(hPionTPCDCAxyNegData);
						fOutputList->Add(hPionTPCDCAxyPosData);
						fOutputList->Add(hProtonTPCDCAxyNegData);
						fOutputList->Add(hProtonTPCDCAxyPosData);
						fOutputList->Add(hPionTOFDCAxyNegData);
						fOutputList->Add(hPionTOFDCAxyPosData);
						fOutputList->Add(hProtonTOFDCAxyNegData);
						fOutputList->Add(hProtonTOFDCAxyPosData);
					}
					fOutputList->Add(hMIPVsEta);
					fOutputList->Add(pMIPVsEta);
					fOutputList->Add(hPlateauVsEta);
					fOutputList->Add(pPlateauVsEta);
				}
			} 
		}
	}

	if (fUseMC) {
		hPtPrimIn =
			new TH1D("hPtPrimIn", "Prim In; #it{p}_{T} (GeV/#it{c}; counts)",
					nPtbins, Ptbins);
		fOutputList->Add(hPtPrimIn);

		hPtPrimOut =
			new TH1D("hPtPrimOut", "Prim Out; #it{p}_{T} (GeV/#it{c}; counts)",
					nPtbins, Ptbins);
		fOutputList->Add(hPtPrimOut);

		hPtSecOut =
			new TH1D("hPtSecOut", "Sec Out; #it{p}_{T} (GeV/#it{c}; counts)",
					nPtbins, Ptbins);
		fOutputList->Add(hPtSecOut);

		hPtOut = new TH1D("hPtOut", "all Out; #it{p}_{T} (GeV/#it{c}; counts)",
				nPtbins, Ptbins);
		fOutputList->Add(hPtOut);

		hFlatenicityMC =
			new TH1D("hFlatenicityMC", "counter", nbins_flat, min_flat, max_flat);
		fOutputList->Add(hFlatenicityMC);
		hFlatResponse =
			new TH2D("hFlatResponse", "; true flat; measured flat", nbins_flat,
					min_flat, max_flat, nbins_flat, min_flat, max_flat);
		fOutputList->Add(hFlatResponse);
		hFlatVsPtMC =
			new TH2D("hFlatVsPtMC", "MC true; Flatenicity; #it{p}_{T} (GeV/#it{c})",
					nbins_flat, min_flat, max_flat, nPtbins, Ptbins);
		fOutputList->Add(hFlatVsPtMC);

		hFlatVsNchMC = new TH2D("hFlatVsNchMC", "; true flat; true Nch", nbins_flat,
				min_flat, max_flat, 100, -0.5, 99.5);
		fOutputList->Add(hFlatVsNchMC);

		hMCPtPionPos = new TH1F("hMCPtPionPos",";#it{p}_{T} (GeV/#it{c}); Counts;", nPtbins, Ptbins);
		fOutputList->Add(hMCPtPionPos);
		hMCPtKaonPos = new TH1F("hMCPtKaonPos",";#it{p}_{T} (GeV/#it{c}); Counts;", nPtbins, Ptbins);
		fOutputList->Add(hMCPtKaonPos);
		hMCPtProtonPos = new TH1F("hMCPtProtonPos",";#it{p}_{T} (GeV/#it{c}); Counts;", nPtbins, Ptbins);
		fOutputList->Add(hMCPtProtonPos);

		hMCPtPionNeg = new TH1F("hMCPtPionNeg",";#it{p}_{T} (GeV/#it{c}); Counts;", nPtbins, Ptbins);
		fOutputList->Add(hMCPtPionNeg);
		hMCPtKaonNeg = new TH1F("hMCPtKaonNeg",";#it{p}_{T} (GeV/#it{c}); Counts;", nPtbins, Ptbins);
		fOutputList->Add(hMCPtKaonNeg);
		hMCPtProtonNeg = new TH1F("hMCPtProtonNeg",";#it{p}_{T} (GeV/#it{c}); Counts;", nPtbins, Ptbins);
		fOutputList->Add(hMCPtProtonNeg);

		hTPCRecTracksPionPos = new TH1F("hTPCRecTracksPionPos",";#it{p}_{T} (GeV/#it{c}); Counts;", nPtbins, Ptbins);
		fOutputList->Add(hTPCRecTracksPionPos);
		hTPCRecTracksKaonPos = new TH1F("hTPCRecTracksKaonPos",";#it{p}_{T} (GeV/#it{c}); Counts;", nPtbins, Ptbins);
		fOutputList->Add(hTPCRecTracksKaonPos);
		hTPCRecTracksProtonPos = new TH1F("hTPCRecTracksProtonPos",";#it{p}_{T} (GeV/#it{c}); Counts;", nPtbins, Ptbins);
		fOutputList->Add(hTPCRecTracksProtonPos);

		hTPCRecTracksPionNeg = new TH1F("hTPCRecTracksPionNeg",";#it{p}_{T} (GeV/#it{c}); Counts;", nPtbins, Ptbins);
		fOutputList->Add(hTPCRecTracksPionNeg);
		hTPCRecTracksKaonNeg = new TH1F("hTPCRecTracksKaonNeg",";#it{p}_{T} (GeV/#it{c}); Counts;", nPtbins, Ptbins);
		fOutputList->Add(hTPCRecTracksKaonNeg);
		hTPCRecTracksProtonNeg = new TH1F("hTPCRecTracksProtonNeg",";#it{p}_{T} (GeV/#it{c}); Counts;", nPtbins, Ptbins);
		fOutputList->Add(hTPCRecTracksProtonNeg);

		hTOFRecTracksPionPos = new TH1F("hTOFRecTracksPionPos",";#it{p}_{T} (GeV/#it{c}); Counts;", nPtbins, Ptbins);
		fOutputList->Add(hTOFRecTracksPionPos);
		hTOFRecTracksKaonPos = new TH1F("hTOFRecTracksKaonPos",";#it{p}_{T} (GeV/#it{c}); Counts;", nPtbins, Ptbins);
		fOutputList->Add(hTOFRecTracksKaonPos);
		hTOFRecTracksProtonPos = new TH1F("hTOFRecTracksProtonPos",";#it{p}_{T} (GeV/#it{c}); Counts;", nPtbins, Ptbins);
		fOutputList->Add(hTOFRecTracksProtonPos);

		hTOFRecTracksPionNeg = new TH1F("hTOFRecTracksPionNeg",";#it{p}_{T} (GeV/#it{c}); Counts;", nPtbins, Ptbins);
		fOutputList->Add(hTOFRecTracksPionNeg);
		hTOFRecTracksKaonNeg = new TH1F("hTOFRecTracksKaonNeg",";#it{p}_{T} (GeV/#it{c}); Counts;", nPtbins, Ptbins);
		fOutputList->Add(hTOFRecTracksKaonNeg);
		hTOFRecTracksProtonNeg = new TH1F("hTOFRecTracksProtonNeg",";#it{p}_{T} (GeV/#it{c}); Counts;", nPtbins, Ptbins);
		fOutputList->Add(hTOFRecTracksProtonNeg);

		hrTPCRecTracksPion = new TH1F("hrTPCRecTracksPion",";#it{p}_{T} (GeV/#it{c}); Counts;", nPtbins, Ptbins);
		fOutputList->Add(hrTPCRecTracksPion);
		hrTPCRecTracksKaon = new TH1F("hrTPCRecTracksKaon",";#it{p}_{T} (GeV/#it{c}); Counts;", nPtbins, Ptbins);
		fOutputList->Add(hrTPCRecTracksKaon);
		hrTPCRecTracksProton = new TH1F("hrTPCRecTracksProton",";#it{p}_{T} (GeV/#it{c}); Counts;", nPtbins, Ptbins);
		fOutputList->Add(hrTPCRecTracksProton);

		for(Int_t i = 0; i < 3; ++i){
			hPionTOFDCAxyNeg[i] = new TH2F(Form("hPion%sTOFDCAxyNeg",ParticleType[i]),"; #it{p}_{T} (GeV/#it{c}); DCA_{xy}", nPtbins, Ptbins, 1000, -3.5, 3.5 );	
			fOutputList->Add(hPionTOFDCAxyNeg[i]);
			hProtonTOFDCAxyNeg[i] = new TH2F(Form("hProton%sTOFDCAxyNeg",ParticleType[i]),"; #it{p}_{T} (GeV/#it{c}); DCA_{xy}", nPtbins, Ptbins, 1000, -3.5, 3.5 );	
			fOutputList->Add(hProtonTOFDCAxyNeg[i]);
			hPionTOFDCAxyPos[i] = new TH2F(Form("hPion%sTOFDCAxyPos",ParticleType[i]),"; #it{p}_{T} (GeV/#it{c}); DCA_{xy}", nPtbins, Ptbins, 1000, -3.5, 3.5 );	
			fOutputList->Add(hPionTOFDCAxyPos[i]);
			hProtonTOFDCAxyPos[i] = new TH2F(Form("hProton%sTOFDCAxyPos",ParticleType[i]),"; #it{p}_{T} (GeV/#it{c}); DCA_{xy}", nPtbins, Ptbins, 1000, -3.5, 3.5 );	
			fOutputList->Add(hProtonTOFDCAxyPos[i]);

			hPionTPCDCAxyNeg[i] = new TH2F(Form("hPion%sTPCDCAxyNeg",ParticleType[i]),"; #it{p}_{T} (GeV/#it{c}); DCA_{xy}", nPtbins, Ptbins, 1000, -3.5, 3.5 );	
			fOutputList->Add(hPionTPCDCAxyNeg[i]);
			hProtonTPCDCAxyNeg[i] = new TH2F(Form("hProton%sTPCDCAxyNeg",ParticleType[i]),"; #it{p}_{T} (GeV/#it{c}); DCA_{xy}", nPtbins, Ptbins, 1000, -3.5, 3.5 );	
			fOutputList->Add(hProtonTPCDCAxyNeg[i]);
			hPionTPCDCAxyPos[i] = new TH2F(Form("hPion%sTPCDCAxyPos",ParticleType[i]),"; #it{p}_{T} (GeV/#it{c}); DCA_{xy}", nPtbins, Ptbins, 1000, -3.5, 3.5 );	
			fOutputList->Add(hPionTPCDCAxyPos[i]);
			hProtonTPCDCAxyPos[i] = new TH2F(Form("hProton%sTPCDCAxyPos",ParticleType[i]),"; #it{p}_{T} (GeV/#it{c}); DCA_{xy}", nPtbins, Ptbins, 1000, -3.5, 3.5 );	
			fOutputList->Add(hProtonTPCDCAxyPos[i]);

		}

		for (Int_t i_eta = 0; i_eta < nEta; ++i_eta) 
		{	
			pion_cont_in_kaon_h[i_eta] = new TH2F(Form("pion_cont_in_kaon_h_%s",etaClass[i_eta]),"; #it{p}_{T} (GeV/#it{c}); n#sigma",nPtbins,Ptbins,nnSigmabins,nSigmabins);
			fOutputList->Add(pion_cont_in_kaon_h[i_eta]);
			nsigma_kaon_h[i_eta] = new TH2F(Form("nsigma_kaon_h_%s",etaClass[i_eta]),"; #it{p}_{T} (GeV/#it{c}); n#sigma",nPtbins,Ptbins,nnSigmabins,nSigmabins);
			fOutputList->Add(nsigma_kaon_h[i_eta]);
			electron_cont_in_kaon_h[i_eta] = new TH2F(Form("electron_cont_in_kaon_h_%s",etaClass[i_eta]),"; #it{p}_{T} (GeV/#it{c}); n#sigma",nPtbins,Ptbins,nnSigmabins,nSigmabins);
			fOutputList->Add(electron_cont_in_kaon_h[i_eta]);

			pion_cont_in_proton_h[i_eta] = new TH2F(Form("pion_cont_in_proton_h_%s",etaClass[i_eta]),"; #it{p}_{T} (GeV/#it{c}); n#sigma",nPtbins,Ptbins,nnSigmabins,nSigmabins);
			fOutputList->Add(pion_cont_in_proton_h[i_eta]);
			nsigma_proton_h[i_eta] = new TH2F(Form("nsigma_proton_h_%s",etaClass[i_eta]),"; #it{p}_{T} (GeV/#it{c}); n#sigma",nPtbins,Ptbins,nnSigmabins,nSigmabins);
			fOutputList->Add(nsigma_proton_h[i_eta]);
			electron_cont_in_proton_h[i_eta] = new TH2F(Form("electron_cont_in_proton_h_%s",etaClass[i_eta]),"; #it{p}_{T} (GeV/#it{c}); n#sigma",nPtbins,Ptbins,nnSigmabins,nSigmabins);
			fOutputList->Add(electron_cont_in_proton_h[i_eta]);

			kaon_cont_in_pion_h[i_eta] = new TH2F(Form("kaon_cont_in_pion_h_%s",etaClass[i_eta]),"; #it{p}_{T} (GeV/#it{c}); n#sigma",nPtbins,Ptbins,nnSigmabins,nSigmabins);
			fOutputList->Add(kaon_cont_in_pion_h[i_eta]);
			nsigma_pion_h[i_eta] = new TH2F(Form("nsigma_pion_h_%s",etaClass[i_eta]),"; #it{p}_{T} (GeV/#it{c}); n#sigma",nPtbins,Ptbins,nnSigmabins,nSigmabins);
			fOutputList->Add(nsigma_pion_h[i_eta]);
			electron_cont_in_pion_h[i_eta] = new TH2F(Form("electron_cont_in_pion_h_%s",etaClass[i_eta]),"; #it{p}_{T} (GeV/#it{c}); n#sigma",nPtbins,Ptbins,nnSigmabins,nSigmabins);
			fOutputList->Add(electron_cont_in_pion_h[i_eta]);
		}
	}

	// x: Sector y: Multiplicity z: V0M Multiplicity
	hActivityV0CV0A = new TH3F("hActivityV0CV0A", "; VZERO channel; #it{N}_{ch} per VZERO cahnnel; V0M quantile", nV0Sectorsbins, V0Sectorsbins, nV0Multbins, V0Multbins, nCent, centClass );
	fOutputList->Add(hActivityV0CV0A);

	hActivityV0DataSect = new TProfile("hActivityV0DataSect", "rec; V0 sector; #LTmultiplicity#GT", 64, -0.5, 63.5);
	fOutputList->Add(hActivityV0DataSect);
	if (fUseMC) {
		hActivityV0McSect =
			new TProfile("hActivityV0McSect", "true; V0 sector; #LTmultiplicity#GT",
					64, -0.5, 63.5);
		fOutputList->Add(hActivityV0McSect);
	}

	hFlatVsV0MVsMult = new TH3F("hFlatVsV0MVsMult", "; Flatenicity; #it{N}_{ch}; V0M Percentile", 
			nFlatbins, Flatbins, nMultbins, Multbins, nCent, centClass);
	fOutputList->Add(hFlatVsV0MVsMult);

	fEventCuts.AddQAplotsToList(fOutputList);
	PostData(1, fOutputList); // postdata will notify the analysis manager of
				  // changes / updates to the
}
//_____________________________________________________________________________
void AliAnalysisTaskFlatenicityPiKp::UserExec(Option_t *) {

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
		if (genHeader) {
			genHeader->PrimaryVertex(vtxMC);
		}
		if (TMath::Abs(vtxMC[2]) <= 10)
			isGoodVtxPosMC = kTRUE;
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

	// Multiplicity Estimation
	fv0mpercentile = -999;

	fMultSelection = (AliMultSelection *)fESD->FindListObject("MultSelection");
	if (!fMultSelection)
		cout << "------- No AliMultSelection Object Found --------"
			<< fMultSelection << endl;
	fv0mpercentile = fMultSelection->GetMultiplicityPercentile("V0M");

	for (Int_t i_c = 0; i_c < nCent; ++i_c) {
		if (fv0mpercentile >= centClass[i_c] &&
				fv0mpercentile < centClass[i_c + 1]) {
			fV0Mindex = i_c;
		} else {
			continue;
		}
	}

	fMidRapidityMult = GetMidRapidityMultiplicity();
	//fFlatTPC = GetFlatenicityTPC(); 
	fFlat = GetFlatenicityV0();

	fFlatMC = -1;
	if (fUseMC) {
		fFlatMC = GetFlatenicityMC();
		if (fFlatMC >= 0) {
			hFlatenicityMC->Fill(fFlatMC);
			hFlatResponse->Fill(fFlatMC, fFlat);
			MakeMCanalysis();
			MakeMCanalysisPID();
			nSigmaContamination();
		}
	}

	if (fFlat > 0 && fV0Mindex >= 0) {

		hFlatVsV0MVsMult->Fill(fFlat, fMidRapidityMult, fv0mpercentile);
		//MakeDataanalysis();
		MakePIDanalysis();
	}

	if (fIsMCclosure) {
		Double_t randomUE = -1;
		gRandom->SetSeed(0);
		randomUE = gRandom->Uniform(0.0, 1.0);
		if (randomUE < 0.5) { // corrections (50% stat.)
			if (isGoodVtxPosMC) {
			}
		} else { // for testing the method
		}
	} else {
		if (fUseMC) {
			if (isGoodVtxPosMC) {
			}
		} else {
		}
	}

	PostData(1, fOutputList); // stream the result of this event to the output
				  // manager which will write it to a file
}
//______________________________________________________________________________
void AliAnalysisTaskFlatenicityPiKp::Terminate(Option_t *) {}
//______________________________________________________________________________
void AliAnalysisTaskFlatenicityPiKp::MakeDataanalysis() {

	// rec
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
		hFlatVsPtV0M[fV0Mindex]->Fill(fFlat, esdtrack->Pt());
	}
}

//______________________________________________________________________________
void AliAnalysisTaskFlatenicityPiKp::MakePIDanalysis() {


	Int_t nESDTracks = fESD->GetNumberOfTracks();
	for (Int_t iT = 0; iT < nESDTracks; iT++) {

		AliESDtrack* esdTrack = static_cast<AliESDtrack *>(fESD->GetTrack(iT));

		if (!esdTrack)
			continue;

		if (!fTrackFilterPID->IsSelected(esdTrack))
			continue;

		if (TMath::Abs(esdTrack->Eta()) > fEtaCut )
			continue;

		if (esdTrack->Pt() < fPtMin)
			continue;

		Double_t momentum = esdTrack->P();
		Double_t pt = esdTrack->Pt();
		Double_t eta = esdTrack->Eta();
		Double_t phi = esdTrack->Phi();
		Float_t dEdx = esdTrack->GetTPCsignal();
		UShort_t Ncl = esdTrack->GetTPCsignalN();
		Short_t charge = esdTrack->Charge();
		Float_t DCAxy = 0.0;
		Float_t dcaz = 0.0;
		esdTrack->GetImpactParameters(DCAxy,dcaz);

		Float_t nSigmaPi = fPIDResponse->NumberOfSigmasTPC(esdTrack,AliPID::kPion);
		Float_t nSigmaK = fPIDResponse->NumberOfSigmasTPC(esdTrack,AliPID::kKaon);
		Float_t nSigmaP = fPIDResponse->NumberOfSigmasTPC(esdTrack,AliPID::kProton);
		Float_t nSigmaPiTOF = fPIDResponse->NumberOfSigmasTOF(esdTrack,AliPID::kPion);
		Float_t nSigmaPTOF = fPIDResponse->NumberOfSigmasTOF(esdTrack,AliPID::kProton);

		Int_t nh = -1;
		if (TMath::Abs(eta)<0.2)
			nh = 0;
		else if (TMath::Abs(eta)>=0.2 && TMath::Abs(eta)<0.4)
			nh = 1;
		else if (TMath::Abs(eta)>=0.4 && TMath::Abs(eta)<0.6)
			nh = 2;
		else if (TMath::Abs(eta)>=0.6 && TMath::Abs(eta)<0.8)
			nh = 3;

		if (nh<0)
			continue;

		//! These are used to measure the Feed-Down correction

		if( charge < 0 ){
			if( (nSigmaPi >= -2.0) && (nSigmaPi <= 2.0) ) { hPionTPCDCAxyNegData->Fill(pt,DCAxy); }
			if( (nSigmaP >= -2.0) && (nSigmaP <= 2.0) ) { hProtonTPCDCAxyNegData->Fill(pt,DCAxy); }

			if( TOFPID(esdTrack) ){
				if( TMath::Sqrt(nSigmaPi*nSigmaPi + nSigmaPiTOF*nSigmaPiTOF) < 2.0) { hPionTOFDCAxyNegData->Fill(pt,DCAxy); }
				if( TMath::Sqrt(nSigmaP*nSigmaP + nSigmaPTOF*nSigmaPTOF) < 2.0) { hProtonTOFDCAxyNegData->Fill(pt,DCAxy); }
			}

		}else{
			if( (nSigmaPi >= -2.0) && (nSigmaPi <= 2.0) ) { hPionTPCDCAxyPosData->Fill(pt,DCAxy); }
			if( (nSigmaP >= -2.0) && (nSigmaP <= 2.0)) { hProtonTPCDCAxyPosData->Fill(pt,DCAxy); }

			if( TOFPID(esdTrack) ){
				if(TMath::Sqrt(nSigmaPi*nSigmaPi + nSigmaPiTOF*nSigmaPiTOF) < 2.0) { hPionTOFDCAxyPosData->Fill(pt,DCAxy); }
				if(TMath::Sqrt(nSigmaP*nSigmaP + nSigmaPTOF*nSigmaPTOF) < 2.0) { hProtonTOFDCAxyPosData->Fill(pt,DCAxy); }
			}
		}


		Double_t maxDCAxy = 10.0;
		maxDCAxy = fcutDCAxy->Eval(pt);

		if (TMath::Abs(DCAxy) > maxDCAxy )
			continue;

		if( charge < 0.0){
			hNsigmaPiNeg[fV0Mindex][nh]->Fill(pt,nSigmaPi,fFlat);
			hNsigmaKNeg[fV0Mindex][nh]->Fill(pt,nSigmaK,fFlat);
			hNsigmaPNeg[fV0Mindex][nh]->Fill(pt,nSigmaP,fFlat);
			hPtTPCEtaNeg[fV0Mindex][nh]->Fill(pt,fFlat);
		}else{
			hNsigmaPiPos[fV0Mindex][nh]->Fill(pt,nSigmaPi,fFlat);
			hNsigmaKPos[fV0Mindex][nh]->Fill(pt,nSigmaK,fFlat);
			hNsigmaPPos[fV0Mindex][nh]->Fill(pt,nSigmaP,fFlat);
			hPtTPCEtaPos[fV0Mindex][nh]->Fill(pt,fFlat);
		} 

		if( TOFPID(esdTrack) ){

			Double_t trkLength = esdTrack->GetIntegratedLength();
			Double_t beta = trkLength/((esdTrack->GetTOFsignal()-fPIDResponse->GetTOFResponse().GetStartTime(esdTrack->P()))*C_Value);

			if(esdTrack->Charge() < 0.0){
				hBetaNeg[fV0Mindex][nh]->Fill(momentum,beta,fFlat);
				hMomentumTOFEtaNeg[fV0Mindex][nh]->Fill(momentum,fFlat);
				hPtTOFEtaNeg[fV0Mindex][nh]->Fill(pt,fFlat);
			}else{ 
				hBetaPos[fV0Mindex][nh]->Fill(momentum,beta,fFlat);
				hMomentumTOFEtaPos[fV0Mindex][nh]->Fill(momentum,fFlat);
				hPtTOFEtaPos[fV0Mindex][nh]->Fill(pt,fFlat);
			}
		}

		hPtVsP[nh]->Fill(momentum,pt);

		if(!PhiCut(esdTrack->Pt(), phi, esdTrack->Charge(), 1.0, fcutLow, fcutHigh))
			continue;

		if (Ncl < fNcl)
			continue;

		if (fdEdxCalibrated) dEdx *= 50.0/EtaCalibration(eta);

		hdEdx[fV0Mindex][nh]->Fill(momentum,dEdx,fFlat);
		hPtrTPC[fV0Mindex][nh]->Fill(pt,fFlat);

		Bool_t IsTOFout=kFALSE;
		if ((esdTrack->GetStatus()&AliESDtrack::kTOFout)==0)
			IsTOFout=kTRUE;

		Float_t lengthtrack = esdTrack->GetIntegratedLength();
		Float_t timeTOF = esdTrack->GetTOFsignal();
		Double_t inttime[5]={0,0,0,0,0};
		esdTrack->GetIntegratedTimes(inttime);// Returns the array with integrated times for each particle hypothesis
		Float_t beta = -99;
		if ( !IsTOFout ){
			if ( ( lengthtrack != 0 ) && ( timeTOF != 0) )
				beta = inttime[0] / timeTOF;
		}

		if ( momentum <= 0.6 && momentum >= 0.4 ){ //only p:0.4-0.6 GeV, pion MIP
			if( dEdx < 60.0 && dEdx > 40.0 ){
				hMIPVsEta->Fill(eta,dEdx);
				pMIPVsEta->Fill(eta,dEdx);
			}
			if( dEdx > 70.0 && dEdx < 90.0 ){
				if(TMath::Abs(beta-1)<0.1){
					hPlateauVsEta->Fill(eta,dEdx);
					pPlateauVsEta->Fill(eta,dEdx);
				}
			}
		}


	}//end of track loop


}

//______________________________________________________________________________
void AliAnalysisTaskFlatenicityPiKp::MakeMCanalysis() {

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
		hFlatVsPtMC->Fill(fFlatMC, particle->Pt());
		hPtPrimIn->Fill(particle->Pt());
	}
	// rec
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
		hPtOut->Fill(esdtrack->Pt());
		Int_t mcLabel = -1;
		mcLabel = TMath::Abs(esdtrack->GetLabel());
		if (fMC->IsPhysicalPrimary(mcLabel)) {
			hPtPrimOut->Fill(esdtrack->Pt());
		} else {
			hPtSecOut->Fill(esdtrack->Pt());
		}
	}
}

//______________________________________________________________________________

void AliAnalysisTaskFlatenicityPiKp::MakeMCanalysisPID() {

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

		hFlatVsPtMC->Fill(fFlatMC, particle->Pt());
		hPtPrimIn->Fill(particle->Pt());

		Int_t pdgCode = particle->PdgCode();
		Int_t pidCode = GetPidCode(pdgCode);

		if(particle->Charge() < 0.0){
			if (pidCode==1) hMCPtPionNeg->Fill(particle->Pt());	
			if (pidCode==2) hMCPtKaonNeg->Fill(particle->Pt());	
			if (pidCode==3) hMCPtProtonNeg->Fill(particle->Pt());	
		}else{
			if (pidCode==1) hMCPtPionPos->Fill(particle->Pt());	
			if (pidCode==2) hMCPtKaonPos->Fill(particle->Pt());	
			if (pidCode==3) hMCPtProtonPos->Fill(particle->Pt());	

		}
	}

	// rec
	Int_t nTracks = fESD->GetNumberOfTracks();
	for (Int_t iT = 0; iT < nTracks; iT++){

		AliESDtrack *esdTrack = static_cast<AliESDtrack *>(
				fESD->GetTrack(iT)); // get a track (type AliesdTrack)
		if (!esdTrack)
			continue;

		if (!fTrackFilterPID->IsSelected(esdTrack))
			continue;

		if (TMath::Abs(esdTrack->Eta()) > fEtaCut)
			continue;

		if (esdTrack->Pt() < fPtMin)
			continue;

		Int_t mcLabel = -1;
		mcLabel = TMath::Abs(esdTrack->GetLabel());

		AliMCParticle *mcTrack = 0;
		mcTrack = (AliMCParticle*)fMC->GetTrack(mcLabel);

		if (!mcTrack) 
			continue;

		if ( TMath::Abs(esdTrack->Charge())==0 )
			continue;

		Int_t pdgCode = mcTrack->PdgCode();
		Int_t pidCode = GetPidCode(pdgCode);
		Double_t pt = esdTrack->Pt();

		Double_t maxDCAxy = 10.0;
		maxDCAxy = fcutDCAxy->Eval(pt);
		Float_t DCAxy = 0.0;
		Float_t dcaz = 0.0;
		esdTrack->GetImpactParameters(DCAxy,dcaz);

		if( fMC->IsPhysicalPrimary(mcLabel) ){

			if( TOFPID(esdTrack) ){

				if(esdTrack->Charge() < 0.0){
					if (pidCode==1) hPionTOFDCAxyNeg[0]->Fill(pt,DCAxy);	
					if (pidCode==3) hProtonTOFDCAxyNeg[0]->Fill(pt,DCAxy);	
				}else{ 
					if (pidCode==1) hPionTOFDCAxyPos[0]->Fill(pt,DCAxy);	
					if (pidCode==3) hProtonTOFDCAxyPos[0]->Fill(pt,DCAxy);	
				}
			}

			if(esdTrack->Charge() < 0.0){
				if (pidCode==1) hPionTPCDCAxyNeg[0]->Fill(pt,DCAxy);	
				if (pidCode==3) hProtonTPCDCAxyNeg[0]->Fill(pt,DCAxy);	
			}else{
				if (pidCode==1) hPionTPCDCAxyPos[0]->Fill(pt,DCAxy);	
				if (pidCode==3) hProtonTPCDCAxyPos[0]->Fill(pt,DCAxy);	
			}
		}

		if( fMC->IsSecondaryFromMaterial(mcLabel) ){

			if( TOFPID(esdTrack) ){

				if(esdTrack->Charge() < 0.0){
					if (pidCode==1) hPionTOFDCAxyNeg[1]->Fill(pt,DCAxy);	
					if (pidCode==3) hProtonTOFDCAxyNeg[1]->Fill(pt,DCAxy);	
				}else{ 
					if (pidCode==1) hPionTOFDCAxyPos[1]->Fill(pt,DCAxy);	
					if (pidCode==3) hProtonTOFDCAxyPos[1]->Fill(pt,DCAxy);	
				}
			}

			if(esdTrack->Charge() < 0.0){
				if (pidCode==1) hPionTPCDCAxyNeg[1]->Fill(pt,DCAxy);	
				if (pidCode==3) hProtonTPCDCAxyNeg[1]->Fill(pt,DCAxy);	
			}else{
				if (pidCode==1) hPionTPCDCAxyPos[1]->Fill(pt,DCAxy);	
				if (pidCode==3) hProtonTPCDCAxyPos[1]->Fill(pt,DCAxy);	
			}
		}

		if( fMC->IsSecondaryFromWeakDecay(mcLabel) ){

			if( TOFPID(esdTrack) ){

				if(esdTrack->Charge() < 0.0){
					if (pidCode==1) hPionTOFDCAxyNeg[2]->Fill(pt,DCAxy);	
					if (pidCode==3) hProtonTOFDCAxyNeg[2]->Fill(pt,DCAxy);	
				}else{ 
					if (pidCode==1) hPionTOFDCAxyPos[2]->Fill(pt,DCAxy);	
					if (pidCode==3) hProtonTOFDCAxyPos[2]->Fill(pt,DCAxy);	
				}
			}

			if(esdTrack->Charge() < 0.0){
				if (pidCode==1) hPionTPCDCAxyNeg[2]->Fill(pt,DCAxy);	
				if (pidCode==3) hProtonTPCDCAxyNeg[2]->Fill(pt,DCAxy);	
			}else{
				if (pidCode==1) hPionTPCDCAxyPos[2]->Fill(pt,DCAxy);	
				if (pidCode==3) hProtonTPCDCAxyPos[2]->Fill(pt,DCAxy);	
			}
		}


		//! DCAxy cut to select primaries
		if (TMath::Abs(DCAxy) > maxDCAxy )
			continue;

		if( fMC->IsPhysicalPrimary(mcLabel) ){

			if( TOFPID(esdTrack) ){

				if(esdTrack->Charge() < 0.0){
					if (pidCode==1) hTOFRecTracksPionNeg->Fill(esdTrack->Pt());	
					if (pidCode==2)	hTOFRecTracksKaonNeg->Fill(esdTrack->Pt());	
					if (pidCode==3) hTOFRecTracksProtonNeg->Fill(esdTrack->Pt());	
				}else{ 
					if (pidCode==1) hTOFRecTracksPionPos->Fill(esdTrack->Pt());	
					if (pidCode==2)	hTOFRecTracksKaonPos->Fill(esdTrack->Pt());	
					if (pidCode==3) hTOFRecTracksProtonPos->Fill(esdTrack->Pt());	
				}
			}

			if(esdTrack->Charge() < 0.0){
				if (pidCode==1) hTPCRecTracksPionNeg->Fill(esdTrack->Pt());	
				if (pidCode==2) hTPCRecTracksKaonNeg->Fill(esdTrack->Pt());	
				if (pidCode==3) hTPCRecTracksProtonNeg->Fill(esdTrack->Pt());	
			}else{
				if (pidCode==1) hTPCRecTracksPionPos->Fill(esdTrack->Pt());	
				if (pidCode==2) hTPCRecTracksKaonPos->Fill(esdTrack->Pt());	
				if (pidCode==3) hTPCRecTracksProtonPos->Fill(esdTrack->Pt());	
			}

			if (esdTrack->GetTPCsignalN() < fNcl)
				continue;

			if (!PhiCut(esdTrack->Pt(), esdTrack->Phi(), esdTrack->Charge(), 1, fcutLow, fcutHigh))
				continue;

			if (pidCode==1) hrTPCRecTracksPion->Fill(esdTrack->Pt());	
			if (pidCode==2) hrTPCRecTracksKaon->Fill(esdTrack->Pt());	
			if (pidCode==3) hrTPCRecTracksProton->Fill(esdTrack->Pt());	
		}
	}

}

//______________________________________________________________________________
Double_t AliAnalysisTaskFlatenicityPiKp::GetFlatenicityTPC() {

	const int nRings2 = 8;
	const int nSectors2 = 8;
	const int nCells2 = nRings2 * nSectors2;
	//	float maxEta2[nRings2] = {-0.4, 0.0, +0.4, +0.8};
	//	float minEta2[nRings2] = {-0.8, -0.4, +0.0, +0.4};
	float maxEta2[nRings2] = {-0.6, -0.4, -0.2, +0.0, +0.2, +0.4, +0.6, +0.8};
	float minEta2[nRings2] = {-0.8, -0.6, -0.4, -0.2, +0.0, +0.2, +0.4, +0.6};
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

	// only events with tracks in the central region
	if (mRho_glob <= 0.0) { return 9999; }

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
Double_t AliAnalysisTaskFlatenicityPiKp::GetFlatenicityV0() {

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

	fV0MMultiplicity = -1.0;
	// Flatenicity calculation
	const Int_t nRings = 4;
	const Int_t nSectors = 8;
	Float_t minEtaV0C[nRings] = {-3.7, -3.2, -2.7, -2.2};
	Float_t maxEtaV0C[nRings] = {-3.2, -2.7, -2.2, -1.7};
	Float_t maxEtaV0A[nRings] = {5.1, 4.5, 3.9, 3.4};
	Float_t minEtaV0A[nRings] = {4.5, 3.9, 3.4, 2.8};
	// Grid
	const Int_t nCells = nRings * 2 * nSectors;
	Float_t RhoLattice[nCells];
	for (Int_t iCh = 0; iCh < nCells; iCh++) {
		RhoLattice[iCh] = 0.0;
	}

	Int_t nringA = 0;
	Int_t nringC = 0;
	for (Int_t iCh = 0; iCh < nCells; iCh++) {
		Float_t detaV0 = -1;
		Float_t mult = lVV0->GetMultiplicity(iCh);
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

		if (fDeltaV0) { mult /= detaV0; }

		// Removing the detaV0 before
		// performing the equalization 
		if (!fV0MEqualisation) { 
			RhoLattice[iCh] = mult; 
		} 
		else{ 
			// When doing the equalization
			// one should take into account
			// the Delta Eta
			if (iCh < 32) { RhoLattice[iCh]  = mult * 1.405 / fV0CCalibration->Eval(iCh); }
			if (iCh >= 32) { RhoLattice[iCh] = mult * 1.405 / fV0ACalibration->Eval(iCh); }
		}
	}

	float mRho = 0.0;
	float flatenicity = -1;

	for (Int_t iCh = 0; iCh < nCells; iCh++)
	{
		mRho += RhoLattice[iCh];
	}

	float multiplicityV0M = mRho;

	// If the Multiplicity in the V0C + V0A 
	// is equal to zero, the event is discarded
	if (mRho == 0.0) { return -1; }

	// Filling histos with mult info
	for (Int_t iCh = 0; iCh < nCells; iCh++) {
		hActivityV0CV0A->Fill(iCh,RhoLattice[iCh],fv0mpercentile);
		hActivityV0DataSect->Fill(iCh, RhoLattice[iCh]);
	}

	fV0MMultiplicity = mRho;
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
		flatenicity = -1;
	}
	return flatenicity;
}
//______________________________________________________________________________
Double_t AliAnalysisTaskFlatenicityPiKp::GetMidRapidityMultiplicity() {

	int mult_glob = 0;
	int nTracks = fESD->GetNumberOfTracks();
	for (Int_t iT = 0; iT < nTracks; ++iT) 
	{

		AliESDtrack *esdtrack = static_cast<AliESDtrack *>(
				fESD->GetTrack(iT)); // get a track (type AliesdTrack)
		if (!esdtrack)
			continue;
		if (!fTrackFilter->IsSelected(esdtrack))
			continue;
		float eta_a = esdtrack->Eta();

		if (TMath::Abs(eta_a) > fEtaCut)
			continue;
		if (esdtrack->Pt() < fPtMin)
			continue;

		mult_glob++;
	}

	return (double)mult_glob;
}
//______________________________________________________________________________
Double_t AliAnalysisTaskFlatenicityPiKp::GetFlatenicityMC() {

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
			RhoLattice[i_segment] /= deltaEta;
			// Filling histos with mult info
			hActivityV0McSect->Fill(i_segment, RhoLattice[i_segment]);
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

Bool_t AliAnalysisTaskFlatenicityPiKp::HasRecVertex() {

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

	return hasVtx;
}

Bool_t AliAnalysisTaskFlatenicityPiKp::TOFPID(AliESDtrack * track) {
	UInt_t status;
	status=track->GetStatus();

	if (!(status & AliESDtrack::kTOFout) || !(status & AliESDtrack::kTIME))
		return kFALSE;

	if (track->GetIntegratedLength() < 350.)
		return kFALSE;

	if (TMath::Abs(track->GetTOFsignalDx()) > 10.0 || TMath::Abs(track->GetTOFsignalDz()) > 10.0)
		return kFALSE;

	return kTRUE;
}

Double_t AliAnalysisTaskFlatenicityPiKp::EtaCalibration(const Double_t &eta) {

	double aPos = 0.0;
	double bPos = 0.0;
	double cPos = 0.0;
	double dPos = 0.0;
	double ePos = 0.0;
	double fPos = 0.0;
	double gPos = 0.0;
	double hPos = 0.0;

	double aNeg = 0.0;
	double bNeg = 0.0;
	double cNeg = 0.0;
	double dNeg = 0.0;
	double eNeg = 0.0;
	double fNeg = 0.0;
	double gNeg = 0.0;
	double hNeg = 0.0;

	if (fPeriod == "16l"){
		aPos = 49.9216; bPos = 3.07252; cPos = -42.8044; dPos = 259.666; ePos = -910.432; fPos = 1776.09; gPos = -1740.65; hPos = 662.232;
		aNeg = 49.9732; bNeg = 4.03575; cNeg = 65.6189;  dNeg = 374.429; eNeg = 951.459;  fNeg = 1153.75; gNeg = 618.493;  hNeg = 100.499;
	} else if(fPeriod == "16k"){
		aPos = 49.9421; bPos = 2.3446; cPos = -41.2765; dPos = 279.695; ePos = -1027.73; fPos = 2022.84; gPos = -1967.79; hPos = 738.823;
		aNeg = 50.0477; bNeg = 8.27344; cNeg = 125.29;  dNeg = 736.8;   eNeg = 2057.75;  fNeg = 2935.38; gNeg = 2064.03;  hNeg = 565.983;
	} else if(fPeriod == "16deghijop"){
		aPos = 49.9743; bPos = 2.3388; cPos = -44.1496; dPos = 296.029; ePos = -1056.56; fPos = 2031.44; gPos = -1946.51; hPos = 723.89;
		aNeg = 50.0329; bNeg = 6.99747; cNeg = 107.168;  dNeg = 649.001; eNeg = 1875.17;  fNeg = 2785.78; gNeg = 2063.77;  hNeg = 606.868;
	} else if(fPeriod == "17data"){
		aPos = 49.6097; bPos = 0.922856; cPos = -6.57484; dPos = 65.3117; ePos = -372.142; fPos = 950.451; gPos = -1085.27; hPos = 458.144;
		aNeg = 49.6555; bNeg = 6.98696; cNeg = 102.734;  dNeg = 566.424; eNeg = 1513.64;  fNeg = 2092.01; gNeg = 1429.32;  hNeg = 375.642;
	} else{
		aPos = 49.6975; bPos = 2.32535; cPos = -42.6516; dPos = 283.058; ePos = -1009.58; fPos = 1945.89; gPos = -1871.23; hPos = 698.552;
		aNeg = 49.8071; bNeg = 9.78466; cNeg = 120.018;  dNeg = 603.325; eNeg = 1470.92;  fNeg = 1819.63; gNeg = 1073.82;  hNeg = 230.142;
	}

	for (Int_t i = 0; i < 8; ++i){
		fEtaCalibrationPos->SetParameter(i,0);
		fEtaCalibrationNeg->SetParameter(i,0);
	}

	Double_t Calibration = 0.0;

	if (eta<=0.0){
		fEtaCalibrationNeg->SetParameter(0,aNeg);
		fEtaCalibrationNeg->SetParameter(1,bNeg);
		fEtaCalibrationNeg->SetParameter(2,cNeg);
		fEtaCalibrationNeg->SetParameter(3,dNeg);
		fEtaCalibrationNeg->SetParameter(4,eNeg);
		fEtaCalibrationNeg->SetParameter(5,fNeg);
		fEtaCalibrationNeg->SetParameter(6,gNeg);
		fEtaCalibrationNeg->SetParameter(7,hNeg);

		Calibration = fEtaCalibrationNeg->Eval(eta);
	} else {
		fEtaCalibrationPos->SetParameter(0,aPos);
		fEtaCalibrationPos->SetParameter(1,bPos);
		fEtaCalibrationPos->SetParameter(2,cPos);
		fEtaCalibrationPos->SetParameter(3,dPos);
		fEtaCalibrationPos->SetParameter(4,ePos);
		fEtaCalibrationPos->SetParameter(5,fPos);
		fEtaCalibrationPos->SetParameter(6,gPos);
		fEtaCalibrationPos->SetParameter(7,hPos);

		Calibration = fEtaCalibrationPos->Eval(eta);
	}

	return Calibration;
}

void AliAnalysisTaskFlatenicityPiKp::nSigmaContamination() {

	Int_t iTracks(fESD->GetNumberOfTracks());          

	for (Int_t iT = 0; iT < iTracks; ++iT) {                

		AliESDtrack *esdTrack = static_cast<AliESDtrack *>(
				fESD->GetTrack(iT)); // get a track (type AliesdTrack)
		if (!esdTrack)
			continue;

		if (!fTrackFilterPID->IsSelected(esdTrack))
			continue;

		if (TMath::Abs(esdTrack->Eta()) > fEtaCut)
			continue;

		if (esdTrack->Pt() < fPtMin)
			continue;

		Double_t maxDCAxy = 10.0;
		maxDCAxy = fcutDCAxy->Eval(esdTrack->Pt());
		Float_t DCAxy = 0.0;
		Float_t dcaz = 0.0;
		esdTrack->GetImpactParameters(DCAxy,dcaz);

		//! DCAxy cut to select primaries
		if (TMath::Abs(DCAxy) > maxDCAxy )
			continue;

		Int_t mcLabel = -1;
		mcLabel = TMath::Abs(esdTrack->GetLabel());

		if (!fMC->IsPhysicalPrimary(mcLabel) )
			continue;

		AliMCParticle *mcTrack = 0;
		mcTrack = (AliMCParticle*)fMC->GetTrack(mcLabel);

		if (!mcTrack) 
			continue;

		if (TMath::Abs(esdTrack->Charge())==0)
			continue;

		Int_t pdgCode = mcTrack->PdgCode();
		Int_t pidCode = GetPidCode(pdgCode);

		Int_t nh = -1;
		Double_t eta = esdTrack->Eta();
		if (TMath::Abs(eta)<0.2)
			nh = 0;
		else if (TMath::Abs(eta)>=0.2 && TMath::Abs(eta)<0.4)
			nh = 1;
		else if (TMath::Abs(eta)>=0.4 && TMath::Abs(eta)<0.6)
			nh = 2;
		else
			nh = 3;

		if ( nh < 0 )
			continue;

		Float_t nsigma_pi = fPIDResponse->NumberOfSigmasTPC(esdTrack,AliPID::kPion);
		Float_t nsigma_k = fPIDResponse->NumberOfSigmasTPC(esdTrack,AliPID::kKaon);
		Float_t nsigma_p = fPIDResponse->NumberOfSigmasTPC(esdTrack,AliPID::kProton);
		Double_t pt = esdTrack->Pt();

		if (TMath::Abs(nsigma_k) <= 3.0){
			if (pidCode==2) { nsigma_kaon_h[nh]->Fill(pt,nsigma_k); }
			if (pidCode==1) { pion_cont_in_kaon_h[nh]->Fill(pt,nsigma_k); }
			if (pidCode==7) { electron_cont_in_kaon_h[nh]->Fill(pt,nsigma_k); }
		}

		if (TMath::Abs(nsigma_p) <= 3.0){	
			if (pidCode==3) { nsigma_proton_h[nh]->Fill(pt,nsigma_p); }
			if (pidCode==1) { pion_cont_in_proton_h[nh]->Fill(pt,nsigma_p); }
			if (pidCode==7) { electron_cont_in_proton_h[nh]->Fill(pt,nsigma_p); }
		}

		if (TMath::Abs(nsigma_pi) <= 3.0){	
			if (pidCode==1) { nsigma_pion_h[nh]->Fill(pt,nsigma_pi); }
			if (pidCode==2) { kaon_cont_in_pion_h[nh]->Fill(pt,nsigma_pi); }
			if (pidCode==7) { electron_cont_in_pion_h[nh]->Fill(pt,nsigma_pi); }
		}
	} 

}


Bool_t AliAnalysisTaskFlatenicityPiKp::PhiCut(Double_t pt, Double_t phi, Double_t q, Float_t mag, TF1* phiCutLow, TF1* phiCutHigh) {

	if(pt < 2.0)
		return kTRUE;

	if(fESD->GetMagneticField() < 0)    // for negatve polarity field
		phi = TMath::TwoPi() - phi;
	if(q < 0) // for negatve charge
		phi = TMath::TwoPi()-phi;

	phi += TMath::Pi()/18.0; // to center gap in the middle
	phi = fmod(phi, TMath::Pi()/9.0);

	if(phi<phiCutHigh->Eval(pt)
			&& phi>phiCutLow->Eval(pt))
		return kFALSE; // reject track

	//	hPhi[fCentClass]->Fill(pt, phi);

	return kTRUE;
}

Int_t AliAnalysisTaskFlatenicityPiKp::GetPidCode(Int_t pdgCode) {
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

