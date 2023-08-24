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

// This task is intended to measure the Multiplicity uncorrelated systematic uncertainty
// For this purposes only the periods LHC16k and LHC16l are used

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
#include <AliKFVertex.h>
#include <AliKFParticle.h>
using std::cout;
using std::endl;

#include "AliAnalysisTaskSystematicsV0M.h"

static const Int_t nCent = 9;
static const Int_t nEta = 4;

static double centClass[nCent + 1] = {0.0,1.0,5.0,10.0,20.0,30.0,40.0,50.0,70.0,100.0};
static const Char_t* etaClass[nEta] = {"02","24","46","68"};
static const double C_Value = TMath::C()*(1.e2/1.e12); // cm/ps

using namespace std; // std namespace: so you can do things like 'cout' etc

ClassImp(AliAnalysisTaskSystematicsV0M) // classimp: necessary for root

AliAnalysisTaskSystematicsV0M::AliAnalysisTaskSystematicsV0M()
	: AliAnalysisTaskSE(), fESD(0), fEventCuts(0x0), fMCStack(0), fMC(0),
	fUseMC(kFALSE), fV0MMultiplicity(-1.0), fDetFlat("V0"), fIsMCclosure(kFALSE),
	fDeltaV0(kTRUE), fRemoveTrivialScaling(kFALSE), fnGen(-1), fPIDResponse(0x0),
	fTrackFilter(0x0), fTrackFilterPID(0x0), fOutputList(0), fEtaCut(0.8), fPtMin(0.5), fNcl(70), fdEdxCalibrated(kTRUE),
	fEtaCalibrationPos(0), fEtaCalibrationNeg(0), felededxfitPos(0), felededxfitNeg(0),
	fcutLow(0x0), fcutHigh(0x0), fcutDCAxy(0x0), fPeriod("16l"), fIsSystematicVariation(kTRUE), fSystVarTrkCuts(0),
	ftrackmult08(0), fv0mpercentile(0), fFlat(-1), fFlatTPC(-1.),
	fMultSelection(0x0), hFlat(0), pMIPVsEta(0), pPlateauVsEta(0), pMIPVsEtaV0s(0)
{

	for (int i_eta = 0; i_eta < nEta; ++i_eta)
	{
		hNsigmaPiPos[i_eta] = 0;
		hNsigmaKPos[i_eta] = 0;
		hNsigmaPPos[i_eta] = 0;
		hNsigmaPiNeg[i_eta] = 0;
		hNsigmaKNeg[i_eta] = 0;
		hNsigmaPNeg[i_eta] = 0;
		hPtTPCEtaPos[i_eta] = 0;
		hPtTPCEtaNeg[i_eta] = 0;
		hBetaPos[i_eta] = 0;
		hMomentumTOFEtaPos[i_eta] = 0;
		hPtTOFEtaPos[i_eta] = 0;
		hBetaNeg[i_eta] = 0;
		hMomentumTOFEtaNeg[i_eta] = 0;
		hPtTOFEtaNeg[i_eta] = 0;
		hdEdx[i_eta] = 0;
		hPtrTPC[i_eta] = 0;
		hPtVsP[i_eta] = 0;
		histPiV0[i_eta] = 0;
		histPV0[i_eta] = 0;
		histEV0[i_eta] = 0;
		histPiTof[i_eta] = 0;
		histPiTof2[i_eta] = 0;

	}
}
//_____________________________________________________________________________
AliAnalysisTaskSystematicsV0M::AliAnalysisTaskSystematicsV0M(const char *name)
	: AliAnalysisTaskSE(name), fESD(0), fEventCuts(0x0), fMCStack(0), fMC(0),
	fUseMC(kFALSE), fV0MMultiplicity(-1.0), fDetFlat("V0"), fIsMCclosure(kFALSE),
	fDeltaV0(kTRUE), fRemoveTrivialScaling(kFALSE), fnGen(-1), fPIDResponse(0x0),
	fTrackFilter(0x0), fTrackFilterPID(0x0), fOutputList(0), fEtaCut(0.8), fPtMin(0.5), fNcl(70), fdEdxCalibrated(kTRUE), 
	fEtaCalibrationPos(0), fEtaCalibrationNeg(0), felededxfitPos(0), felededxfitNeg(0),
	fcutLow(0x0), fcutHigh(0x0), fcutDCAxy(0x0), fPeriod("16l"), fIsSystematicVariation(kTRUE), fSystVarTrkCuts(0),
	ftrackmult08(0), fv0mpercentile(0), fFlat(-1), fFlatTPC(-1.),
	fMultSelection(0x0), hFlat(0),
	pMIPVsEta(0), pPlateauVsEta(0), pMIPVsEtaV0s(0)
{
	for (int i_eta = 0; i_eta < nEta; ++i_eta)
	{
		hNsigmaPiPos[i_eta] = 0;
		hNsigmaKPos[i_eta] = 0;
		hNsigmaPPos[i_eta] = 0;
		hNsigmaPiNeg[i_eta] = 0;
		hNsigmaKNeg[i_eta] = 0;
		hNsigmaPNeg[i_eta] = 0;
		hPtTPCEtaPos[i_eta] = 0;
		hPtTPCEtaNeg[i_eta] = 0;
		hBetaPos[i_eta] = 0;
		hMomentumTOFEtaPos[i_eta] = 0;
		hPtTOFEtaPos[i_eta] = 0;
		hBetaNeg[i_eta] = 0;
		hMomentumTOFEtaNeg[i_eta] = 0;
		hPtTOFEtaNeg[i_eta] = 0;
		hdEdx[i_eta] = 0;
		hPtrTPC[i_eta] = 0;
		hPtVsP[i_eta] = 0;	
		histPiV0[i_eta] = 0;
		histPV0[i_eta] = 0;
		histEV0[i_eta] = 0;
		histPiTof[i_eta] = 0;
		histPiTof2[i_eta] = 0;
	}
	DefineInput(0, TChain::Class()); // define the input of the analysis: in this
					 // case you take a 'chain' of events
	DefineOutput(1, TList::Class()); // define the ouptut of the analysis: in this
					 // case it's a list of histograms
}
//_____________________________________________________________________________
AliAnalysisTaskSystematicsV0M::~AliAnalysisTaskSystematicsV0M() {
	// destructor
	if (fOutputList) {
		delete fOutputList; // at the end of your task, it is deleted from memory by
				    // calling this function
		fOutputList = 0x0;
	}
}
//_____________________________________________________________________________
void AliAnalysisTaskSystematicsV0M::UserCreateOutputObjects() {

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

	// This is just to ensure that nominal track cuts
	// are used when fIsSystematicVariation = kFALSE
	if (!fIsSystematicVariation) { fSystVarTrkCuts = 0; }

	if(fSystVarTrkCuts==1){ //! Lower: SetMinNCrossedRowsTPC(60)
		fCutsPID->SetMinNCrossedRowsTPC(60);
		fCutsPID->SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);
		fCutsPID->SetMaxChi2PerClusterTPC(4);
		fCutsPID->SetMaxDCAToVertexZ(2);
		fCutsPID->SetMaxChi2PerClusterITS(36);
	}
	else if(fSystVarTrkCuts==2){ //! Higher: SetMinNCrossedRowsTPC(100)
		fCutsPID->SetMinNCrossedRowsTPC(100);
		fCutsPID->SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);
		fCutsPID->SetMaxChi2PerClusterTPC(4);
		fCutsPID->SetMaxDCAToVertexZ(2);
		fCutsPID->SetMaxChi2PerClusterITS(36);
	}
	else if(fSystVarTrkCuts==3){ //! Lower: SetMinRatioCrossedRowsOverFindableClustersTPC(0.7)
		fCutsPID->SetMinNCrossedRowsTPC(70);
		fCutsPID->SetMinRatioCrossedRowsOverFindableClustersTPC(0.7);
		fCutsPID->SetMaxChi2PerClusterTPC(4);
		fCutsPID->SetMaxDCAToVertexZ(2);
		fCutsPID->SetMaxChi2PerClusterITS(36);
	}
	else if(fSystVarTrkCuts==4){ //! Higher: SetMinRatioCrossedRowsOverFindableClustersTPC(0.9)
		fCutsPID->SetMinNCrossedRowsTPC(70);
		fCutsPID->SetMinRatioCrossedRowsOverFindableClustersTPC(0.9);
		fCutsPID->SetMaxChi2PerClusterTPC(4);
		fCutsPID->SetMaxDCAToVertexZ(2);
		fCutsPID->SetMaxChi2PerClusterITS(36);
	}
	else if(fSystVarTrkCuts==5){ //! Lower: SetMaxChi2PerClusterTPC(3)
		fCutsPID->SetMinNCrossedRowsTPC(70);
		fCutsPID->SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);
		fCutsPID->SetMaxChi2PerClusterTPC(3);
		fCutsPID->SetMaxDCAToVertexZ(2);
		fCutsPID->SetMaxChi2PerClusterITS(36);
	}
	else if(fSystVarTrkCuts==6){ //! Higher: SetMaxChi2PerClusterTPC(5)
		fCutsPID->SetMinNCrossedRowsTPC(70);
		fCutsPID->SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);
		fCutsPID->SetMaxChi2PerClusterTPC(5);
		fCutsPID->SetMaxDCAToVertexZ(2);
		fCutsPID->SetMaxChi2PerClusterITS(36);
	}
	else if(fSystVarTrkCuts==7){ //! Lower: SetMaxChi2PerClusterITS(25)
		fCutsPID->SetMinNCrossedRowsTPC(70);
		fCutsPID->SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);
		fCutsPID->SetMaxChi2PerClusterTPC(4);
		fCutsPID->SetMaxDCAToVertexZ(2);
		fCutsPID->SetMaxChi2PerClusterITS(25);
	}
	else if(fSystVarTrkCuts==8){ //! Higher: SetMaxChi2PerClusterITS(49)
		fCutsPID->SetMinNCrossedRowsTPC(70);
		fCutsPID->SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);
		fCutsPID->SetMaxChi2PerClusterTPC(4);
		fCutsPID->SetMaxDCAToVertexZ(2);
		fCutsPID->SetMaxChi2PerClusterITS(49);
	}
	else if(fSystVarTrkCuts==9){ //! Lower: SetMaxDCAToVertexZ(1)
		fCutsPID->SetMinNCrossedRowsTPC(70);
		fCutsPID->SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);
		fCutsPID->SetMaxChi2PerClusterTPC(4);
		fCutsPID->SetMaxDCAToVertexZ(1);
		fCutsPID->SetMaxChi2PerClusterITS(36);
	}
	else if(fSystVarTrkCuts==10){ //! Lower: SetMaxDCAToVertexZ(5)
		fCutsPID->SetMinNCrossedRowsTPC(70);
		fCutsPID->SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);
		fCutsPID->SetMaxChi2PerClusterTPC(4);
		fCutsPID->SetMaxDCAToVertexZ(5);
		fCutsPID->SetMaxChi2PerClusterITS(36);
	}
	else{ //! Nominal values
		fCutsPID->SetMinNCrossedRowsTPC(70);
		fCutsPID->SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);
		fCutsPID->SetMaxChi2PerClusterTPC(4);
		fCutsPID->SetMaxDCAToVertexZ(2);
		fCutsPID->SetMaxChi2PerClusterITS(36);
	}

	std::cout << "GetMinNCrossedRowsTPC = " << fCutsPID->GetMinNCrossedRowsTPC() << '\n';
	std::cout << "GetMinRatioCrossedRowsOverFindableClustersTPC = " << fCutsPID->GetMinRatioCrossedRowsOverFindableClustersTPC() << '\n';
	std::cout << "GetMaxChi2PerClusterTPC = " << fCutsPID->GetMaxChi2PerClusterTPC() << '\n';
	std::cout << "GetMaxDCAToVertexZ = " << fCutsPID->GetMaxDCAToVertexZ() << '\n';
	std::cout << "GetMaxChi2PerClusterITS = " << fCutsPID->GetMaxChi2PerClusterITS() << '\n';

	fTrackFilterPID->AddCuts(fCutsPID);

	const std::map<std::string, std::array<double, 9>> neg_eta_mip{
		{"18b",{49.8668, 12.2863, 259.435, 2060.93, 8446.22, 19854.3, 27051.4, 19835.5, 6041.59}},
			{"18d",{49.9868, 20.5379, 372.253, 2754.23, 10730.7, 24146.8, 31618.4, 22361.2, 6593.44}},
			{"18e",{49.9502, 19.4593, 362.489, 2717.93, 10675.4, 24154.7, 31761.2, 22543.5, 6669.84}},
			{"18f",{49.9843, 19.3543, 355.345, 2659.38, 10440.5, 23610.6, 31019.1, 21989.3, 6495.44}},
			{"18g",{2.41851e-314, 2.41851e-314, 2.41852e-314, 2.41852e-314, 2.41853e-314, 2.41853e-314, 2.41853e-314, 3.02765e-314, 3.02765e-314}},
			{"18h",{2.41851e-314, 2.41851e-314, 2.41852e-314, 2.41852e-314, 2.41853e-314, 2.41853e-314, 2.41853e-314, 3.02765e-314, 3.02765e-314}},
			{"18i",{49.8488, 15.5784, 317.835, 2415.84, 9475.32, 21359.3, 28019.5, 19886.4, 5895.86}},
			{"18j",{2.41851e-314, 2.41851e-314, 2.41852e-314, 2.41852e-314, 2.41853e-314, 2.41853e-314, 2.41853e-314, 3.02765e-314, 3.02765e-314}},
			{"18k",{49.927, 18.6596, 360.176, 2692.77, 10441.7, 23259.7, 30111.7, 21068.5, 6154.48}},
			{"18l",{49.7486, 20.5767, 394.317, 2987.78, 11792.4, 26741.7, 35191.3, 24979.6, 7389.03}},
			{"18m",{49.7772, 20.4881, 371.619, 2757.56, 10763.5, 24233, 31728.5, 22440.1, 6619.96}},
			{"18o",{49.7787, 18.0261, 338.343, 2531.96, 9940.15, 22541.3, 29769.8, 21249, 6326.01}},
			{"16d",{49.9734, 8.26292, 182.129, 1441.79, 5838.95, 13595.8, 18467.1, 13585.3, 4169.83}},
			{"16e",{49.9699, 11.34, 251.41, 2013.77, 8226.55, 19144.7, 25735.7, 18602, 5588.55}},
			{"16g",{50.0053, 18.9008, 385.551, 2988.53, 11902.3, 26999.7, 35344.2, 24877.4, 7286.33}},
			{"16h",{50.0607, 18.6326, 392.199, 3144.56, 12846.9, 29706.4, 39457.8, 28084.3, 8297.16}},
			{"16i",{50.07, 17.8159, 373.451, 2988.92, 12147.5, 27904.1, 36824.4, 26063.5, 7665.41}},
			{"16j",{50.0813, 19.7903, 412.257, 3289.53, 13358, 30679.1, 40470.6, 28615.6, 8401.52}},
			{"16k",{50.1712, 15.966, 294.311, 2286.7, 9311.34, 21743.3, 29397.8, 21402.8, 6483.97}},
			{"16l",{50.2007, 15.9971, 291.383, 2254.17, 9185.31, 21530.5, 29268.9, 21438.3, 6533.7}},
			{"16o",{50.0424, 15.9223, 343.917, 2772.75, 11372.6, 26452.9, 35405.7, 25413.5, 7571.52}},
			{"16p",{50.0558, 17.5436, 368.288, 2900.47, 11622, 26442.1, 34706.4, 24506, 7204.09}}
	};

	const std::map<std::string, std::array<double, 9>> pos_eta_mip{
		{"18b",{49.804, -5.41725, 150.569, -1281.28, 5360.52, -12658.9, 17229.5, -12588.2, 3812.03}},
			{"18d",{49.8368, -0.700672, 92.1073, -1014.41, 4789.02, -12182.8, 17440.7, -13200.9, 4099.71}},
			{"18e",{49.8207, -2.64836, 119.148, -1162.13, 5188.24, -12712, 17682.9, -13082, 3987.94}},
			{"18f",{49.8689, -4.23401, 141.049, -1333.73, 5942.55, -14612.2, 20405.2, -15142.9, 4627.71}},
			{"18g",{49.8785, -7.15882, 174.404, -1459.63, 6097.56, -14398.4, 19570.9, -14263.2, 4306.47}},
			{"18h",{2.41851e-314, 2.41851e-314, 2.41852e-314, 2.41852e-314, 2.41853e-314, 2.41853e-314, 2.41853e-314, 3.02765e-314, 3.02765e-314}},
			{"18i",{49.7489, -4.13738, 140.656, -1249.28, 5305.89, -12612.7, 17225.4, -12611.7, 3825.18}},
			{"18j",{2.41851e-314, 2.41851e-314, 2.41852e-314, 2.41852e-314, 2.41853e-314, 2.41853e-314, 2.41853e-314, 3.02765e-314, 3.02765e-314}},
			{"18k",{49.7975, -2.9696, 127.625, -1212.21, 5296.91, -12736, 17441.6, -12748.8, 3852.52}},
			{"18l",{49.6025, -1.46048, 105.497, -1094.46, 4999.54, -12406, 17422, -13000.2, 3996.79}},
			{"18m",{49.6393, -4.74755, 143.275, -1307.48, 5671.19, -13622.4, 18657.4, -13639.1, 4121.86}},
			{"18o",{49.6328, -0.941546, 94.8448, -1025.84, 4800.12, -12136.6, 17304.4, -13070, 4057.81}},
			{"16d",{49.9598, -6.27035, 165.621, -1388.86, 5852.57, -14055.9, 19514.1, -14540.6, 4485.2}},
			{"16e",{49.9218, -5.5296, 164.462, -1420.32, 6042.03, -14514.7, 20054.1, -14830.4, 4534.49}},
			{"16g",{49.8685, -3.28044, 137.568, -1273.78, 5562.06, -13506.7, 18736.3, -13864.5, 4234.32}},
			{"16h",{49.9075, -2.49821, 111.517, -1093.5, 4958.3, -12387.3, 17582, -13258.9, 4113.83}},
			{"16i",{49.9301, -3.05879, 121.567, -1200.96, 5509.67, -13884.4, 19811.4, -14984.3, 4656.93}},
			{"16j",{49.9267, -1.58978, 96.2437, -986.505, 4519.85, -11306.1, 16031.9, -12075.5, 3742.29}},
			{"16k",{50.028, -1.65448, 46.8913, -523.373, 2705.39, -7593.22, 11921.5, -9793.56, 3264.29}},
			{"16l",{50.0438, -1.71624, 49.5086, -542.689, 2766.84, -7689.99, 11992.4, -9809.14, 3261.24}},
			{"16o",{49.9187, -3.68235, 132.208, -1237.58, 5477.21, -13453.5, 18847.6, -14068.6, 4330.22}},
			{"16p",{49.9181, -3.37394, 130.832, -1251.02, 5578.24, -13730, 19228.6, -14332.4, 4402.7}}
	};

	const auto map_neg_eta_mip = neg_eta_mip.find(fPeriod);
	const std::array<double,9> pars_neg_eta_mip = map_neg_eta_mip->second;

	const auto map_pos_eta_mip = pos_eta_mip.find(fPeriod);
	const std::array<double,9> pars_pos_eta_mip = map_pos_eta_mip->second;

	fEtaCalibrationPos = new TF1("fDeDxVsEtaPos", "pol8", 0.0, 1.0);
	fEtaCalibrationNeg = new TF1("fDeDxVsEtaNeg", "pol8", -1.0, 0.0);

	for (int par = 0; par < 9; ++par)
	{
		fEtaCalibrationNeg->SetParameter(par,0.0);
		fEtaCalibrationPos->SetParameter(par,0.0);
		fEtaCalibrationNeg->SetParameter(par,pars_neg_eta_mip[par]);
		fEtaCalibrationPos->SetParameter(par,pars_pos_eta_mip[par]);
	}

	std::cout << "fPeriod = " << fPeriod << '\n';
	std::cout << "MIP negative eta parameters: " << '\n';
	for (int par = 0; par < 9; ++par)
	{
		std::cout << "par = " << pars_neg_eta_mip[par] << '\n';
	}

	std::cout << "MIP positive eta parameters: " << '\n';
	for (int par = 0; par < 9; ++par)
	{
		std::cout << "par = " << pars_pos_eta_mip[par] << '\n';
	}

	// These are the Plateau calibration parameters
	const std::map<std::string, std::array<double, 5>> neg_eta{
		{"18b",{78.0584, -9.92534, -43.4941, -74.9382, -43.5307}},
			{"18d",{78.2565, -10.3295, -43.4784, -68.978, -36.1958}},
			{"18e",{78.4238, -8.72674, -39.6735, -68.3232, -39.2252}},
			{"18f",{78.3869, -10.1457, -46.4504, -79.7451, -45.6895}},
			{"18g",{78.2767, -8.64306, -33.9562, -51.945, -26.884}},
			{"18h",{78.2488, -13.6366, -62.1845, -106.582, -61.3608}},
			{"18i",{78.1277, -9.55994, -40.8765, -70.2554, -41.2268}},
			{"18j",{75.4719, -74.3785, -361.481, -602.121, -323.142}},
			{"18k",{78.2423, -10.1552, -43.4763, -73.0076, -41.5902}},
			{"18l",{77.8836, -13.14, -57.3167, -94.4722, -52.2078}},
			{"18m",{78.0229, -10.2569, -45.1534, -77.6259, -45.2149}},
			{"18o",{77.9943, -11.072, -46.4452, -74.317, -39.5181}},
			{"16d",{78.3978, -9.39313, -45.2739, -85.8433, -52.9036}},
			{"16e",{78.5115, -10.575, -48.2195, -84.9406, -48.9549}},
			{"16g",{78.3259, -11.7637, -49.5636, -82.687, -45.3562}},
			{"16h",{78.5328, -12.2705, -57.7146, -101.312, -57.5771}},
			{"16i",{78.7061, -10.1328, -50.0004, -90.9223, -53.1454}},
			{"16j",{78.5559, -11.504, -56.6088, -101.704, -58.4783}},
			{"16k",{79.7335, -8.64085, -44.1529, -84.0726, -51.4506}},
			{"16l",{79.8256, -8.85733, -46.6177, -89.6549, -55.4802}},
			{"16o",{78.4585, -12.6594, -61.3969, -110.308, -63.8506}},
			{"16p",{78.4266, -11.1231, -50.8061, -88.0376, -48.9562}}
	};

	const std::map<std::string, std::array<double, 5>> pos_eta{
		{"18b",{78.1418, 12.3115, -58.9967, 106.955, -63.9435}},
			{"18d",{78.7388, 11.0431, -61.3579, 116.231, -70.7895}},
			{"18e",{78.6326, 11.8068, -61.526, 112.217, -66.4458}},
			{"18f",{78.6414, 12.1778, -64.3007, 117.807, -70.0741}},
			{"18g",{78.2821, 13.5853, -68.2264, 125.483, -75.7091}},
			{"18h",{78.7728, 9.77365, -50.9454, 93.2052, -55.6025}},
			{"18i",{78.2028, 13.0052, -62.5986, 112.765, -67.2408}},
			{"18j",{78.1961, 41.8608, -270.187, 573.874, -384.229}},
			{"18k",{78.5361, 12.6516, -65.5564, 120.462, -72.2361}},
			{"18l",{78.3832, 11.1527, -59.8539, 112.631, -68.757}},
			{"18m",{78.3155, 12.5582, -66.5439, 124.659, -76.06}},
			{"18o",{78.4926, 8.86145, -48.2071, 88.5975, -51.5732}},
			{"16d",{78.689, 9.42784, -45.7495, 85.2206, -51.8471}},
			{"16e",{78.668, 10.7527, -48.0636, 83.0091, -47.2603}},
			{"16g",{78.6351, 14.9114, -73.4522, 132.114, -77.3137}},
			{"16h",{78.8393, 11.6098, -61.0424, 114.416, -68.4445}},
			{"16i",{78.8459, 12.3258, -65.2223, 122.571, -73.7102}},
			{"16j",{78.7301, 12.7389, -66.9754, 124.931, -74.3142}},
			{"16k",{80.0008, 6.93231, -42.3245, 87.4325, -55.7795}},
			{"16l",{80.1504, 4.98, -31.3886, 66.3694, -43.2082}},
			{"16o",{78.7704, 10.0829, -53.0814, 100.967, -61.0128}},
			{"16p",{78.7037, 11.2722, -57.0117, 103.02, -58.7785}}
	};

	const auto map_neg_eta = neg_eta.find(fPeriod);
	const std::array<double,5> pars_neg_eta = map_neg_eta->second;

	const auto map_pos_eta = pos_eta.find(fPeriod);
	const std::array<double,5> pars_pos_eta = map_pos_eta->second;

	felededxfitPos = new TF1("felededxfitPos", "pol4", 0.0, 1.0);
	felededxfitNeg = new TF1("felededxfitNeg", "pol4", -1.0, 0.0);

	for(int i = 0; i < 5; ++i)
	{
		felededxfitNeg->SetParameter(i,0.0);
		felededxfitPos->SetParameter(i,0.0);
		felededxfitNeg->SetParameter(i,pars_neg_eta[i]);
		felededxfitPos->SetParameter(i,pars_pos_eta[i]);
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
		20.0,22.0,24.0,26.0,30.0 };

	const int nPtbins_low_pT = 17;
	double Ptbins_low_pT[nPtbins_low_pT+1] = {
		0.25, 0.30, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 
		0.8, 0.85, 0.9, 0.95, 1.0, 1.1, 1.2 };

	const int nPtbins_intermediate_pT = 30;
	double Ptbins_intermediate_pT[nPtbins_intermediate_pT+1] = {
		0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0, 
		1.1,  1.2, 1.3,  1.4, 1.5,  1.6, 1.7,  1.8, 1.9, 
		2.0,  2.2, 2.4,  2.6, 2.8,  3.0, 3.2,  3.4, 3.6, 
		3.8, 4.0, 4.5, 5.0 };

	const int nPtbins_high_pT = 31;
	double Ptbins_high_pT[nPtbins_high_pT+1] = {
		2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 
		4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 8.0, 9.0, 10.0, 
		11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 18.0, 20.0, 22.0, 24.0,
		26.0, 30.0 };

	const int nPtBinsV0s = 25;
	double ptBinsV0s[nPtBinsV0s+1] = { 0.0 , 0.1 , 0.2 , 0.3 , 0.4 , 0.5 , 0.6 , 0.7 , 0.8 , 0.9 , 1.0 ,
		1.2 , 1.4 , 1.6 , 1.8 , 2.0 , 2.5 , 3.0 , 3.5 , 4.0 , 5.0 , 7.0 ,
		9.0 , 12.0, 15.0, 20.0 };

	const int nDeltaPiBins = 80;
	const double deltaPiLow  = 20;
	const double deltaPiHigh = 100;

	int nFlatbins = 0;
	if (!fIsSystematicVariation) {
		nFlatbins = 1020;
	}
	else {
		nFlatbins = 8;
	}

	double Flatbins[nFlatbins+1];

	if (!fIsSystematicVariation) {
		for (int i = 0; i <= nFlatbins; ++i) 
		{
			Flatbins[i] = -0.01 + (double)i * 0.001;
		}
	}
	else {
		const double Flatbins_16kl[9] = { -0.01, 0.102, 0.12, 0.132, 0.151, 0.168, 0.186, 0.205, 1.01 };
		for (int i = 0; i < 9; ++i )
		{
			Flatbins[i] = Flatbins_16kl[i];
		}
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
	fOutputList = new TList(); // this is a list which will contain all of your histograms
	fOutputList->SetOwner(kTRUE); // memory stuff: the list is owner of all

	pMIPVsEta = new TProfile("pMIPVsEta","; #eta; #LT dE/dx #GT_{MIP, primary tracks}", 50, -0.8, 0.8, 40.0, 60.0);
	pPlateauVsEta = new TProfile("pPlateauVsEta","; #eta; #LT dE/dx #GT_{Plateau, primary tracks}", 50, -0.8, 0.8, 60.0, 110.0);
	pMIPVsEtaV0s = new TProfile("pMIPVsEtaV0s","; #eta; #LT dE/dx #GT_{MIP, secondary tracks}", 50, -0.8, 0.8, 40.0, 60.0);
	hFlat = new TH2F("hFlat_vs_V0M", "; V0M; Flatenicity;",nCent,centClass,1020,-0.01,1.01);

	if (!fUseMC) { 
		fOutputList->Add(hFlat);
		fOutputList->Add(pMIPVsEta);
		fOutputList->Add(pPlateauVsEta);
		fOutputList->Add(pMIPVsEtaV0s);
	}

	for (int i_eta = 0; i_eta < nEta; ++i_eta) 
	{
		hNsigmaPiPos[i_eta] = new TH3F(Form("hNsigmaPiPos_eta_%s",etaClass[i_eta]),";#it{p}_{T} (GeV/#it{c}); n#sigma; V0M",nPtbins_low_pT,Ptbins_low_pT,nnSigmabins,nSigmabins,nCent,centClass);
		hNsigmaKPos[i_eta] = new TH3F(Form("hNsigmaKPos_eta_%s",etaClass[i_eta]),";#it{p}_{T} (GeV/#it{c}); n#sigma; V0M",nPtbins_low_pT,Ptbins_low_pT,nnSigmabins,nSigmabins,nCent,centClass);
		hNsigmaPPos[i_eta] = new TH3F(Form("hNsigmaPPos_eta_%s",etaClass[i_eta]),";#it{p}_{T} (GeV/#it{c}); n#sigma; V0M",nPtbins_low_pT,Ptbins_low_pT,nnSigmabins,nSigmabins,nCent,centClass);
		hNsigmaPiNeg[i_eta] = new TH3F(Form("hNsigmaPiNeg_eta_%s",etaClass[i_eta]),";#it{p}_{T} (GeV/#it{c}); n#sigma; V0M",nPtbins_low_pT,Ptbins_low_pT,nnSigmabins,nSigmabins,nCent,centClass);
		hNsigmaKNeg[i_eta] = new TH3F(Form("hNsigmaKNeg_eta_%s",etaClass[i_eta]),";#it{p}_{T} (GeV/#it{c}); n#sigma; V0M",nPtbins_low_pT,Ptbins_low_pT,nnSigmabins,nSigmabins,nCent,centClass);
		hNsigmaPNeg[i_eta] = new TH3F(Form("hNsigmaPNeg_eta_%s",etaClass[i_eta]),";#it{p}_{T} (GeV/#it{c}); n#sigma; V0M",nPtbins_low_pT,Ptbins_low_pT,nnSigmabins,nSigmabins,nCent,centClass);
		hPtTPCEtaPos[i_eta] = new TH2F(Form("hPtTPCEtaPos_eta_%s",etaClass[i_eta]),"; #it{p}_{T} (GeV/#it{c}; V0M",nPtbins_low_pT,Ptbins_low_pT,nCent,centClass);
		hPtTPCEtaNeg[i_eta] = new TH2F(Form("hPtTPCEtaNeg_eta_%s",etaClass[i_eta]),"; #it{p}_{T} (GeV/#it{c}); V0M",nPtbins_low_pT,Ptbins_low_pT,nCent,centClass);
		hBetaPos[i_eta] = new TH3F(Form("hBetaPos_eta_%s",etaClass[i_eta]),";#it{p} (GeV/#it{c}); #beta; V0M",nPtbins_intermediate_pT,Ptbins_intermediate_pT,nBetabins,Betabins,nCent,centClass);
		hBetaNeg[i_eta] = new TH3F(Form("hBetaNeg_eta_%s",etaClass[i_eta]),";#it{p} (GeV/#it{c}); #beta; V0M", nPtbins_intermediate_pT, Ptbins_intermediate_pT, nBetabins,Betabins,nCent,centClass);
		hMomentumTOFEtaPos[i_eta] = new TH2F(Form("hMomentumTOFEtaPos_eta_%s",etaClass[i_eta]),";#it{p} (GeV/#it{c}); V0M", nPtbins_intermediate_pT, Ptbins_intermediate_pT,nCent,centClass);
		hMomentumTOFEtaNeg[i_eta] = new TH2F(Form("hMomentumTOFEtaNeg_eta_%s",etaClass[i_eta]),";#it{p} (GeV/#it{c}); V0M", nPtbins_intermediate_pT, Ptbins_intermediate_pT,nCent,centClass);
		hPtTOFEtaPos[i_eta] = new TH2F(Form("hPtTOFEtaPos_eta_%s",etaClass[i_eta]),";#it{p}_{T} (GeV/#it{c}); V0M", nPtbins_intermediate_pT, Ptbins_intermediate_pT,nCent,centClass);
		hPtTOFEtaNeg[i_eta] = new TH2F(Form("hPtTOFEtaNeg_eta_%s",etaClass[i_eta]),";#it{p}_{T} (GeV/#it{c}); V0M", nPtbins_intermediate_pT, Ptbins_intermediate_pT,nCent,centClass);
		hdEdx[i_eta] = new TH3F(Form("hdEdx_eta_%s",etaClass[i_eta]), ";#it{p} (GeV/#it{c}); dE/dx; V0M", nPtbins_high_pT, Ptbins_high_pT, ndEdxbins, dEdxbins,nCent,centClass);
		hPtrTPC[i_eta] = new TH2F(Form("hPtrTPC_eta_%s",etaClass[i_eta]),";#it{p}_{T} (GeV/#it{c}); V0M", nPtbins_high_pT, Ptbins_high_pT,nCent,centClass);
		hPtVsP[i_eta] = new TH2F(Form("hPtVsP_eta_%s",etaClass[i_eta]), ";#it{p} (GeV/#it{c}); #it{p}_{T} (GeV/#it{c})", nPtbins, Ptbins, nPtbins, Ptbins);
		histPiV0[i_eta] = new TH2F(Form("hPiV0_%s",etaClass[i_eta]), "Pions id by V0; #it{p} (GeV/#it{c}); d#it{e}d#it{x}", nPtBinsV0s, ptBinsV0s, nDeltaPiBins, deltaPiLow, deltaPiHigh);
		histPV0[i_eta] = new TH2F(Form("hPV0_%s",etaClass[i_eta]), "Protons id by V0; #it{p} (GeV/#it{c}); d#it{e}d#it{x}", nPtBinsV0s, ptBinsV0s, nDeltaPiBins, deltaPiLow, deltaPiHigh);
		histPiTof[i_eta] = new TH2F(Form("hPiTOF_%s",etaClass[i_eta]), "Primary Pions from TOF; #it{p} (GeV/#it{c}); d#it{e}d#it{x}", nPtBinsV0s, ptBinsV0s, nDeltaPiBins, deltaPiLow, deltaPiHigh);
		histPiTof2[i_eta] = new TH2F(Form("hPiTOF2_%s",etaClass[i_eta]), "Primary Pions from TOF; #it{p} (GeV/#it{c}); d#it{e}d#it{x}", nPtBinsV0s, ptBinsV0s, nDeltaPiBins, deltaPiLow, deltaPiHigh);
		histEV0[i_eta] = new TH2F(Form("hEV0_%s",etaClass[i_eta]), "Electrons id by V0; #it{p} (GeV/#it{c}); d#it{e}d#it{x}", nPtBinsV0s, ptBinsV0s, nDeltaPiBins, deltaPiLow, deltaPiHigh);

		if (!fUseMC) { 

			fOutputList->Add(hNsigmaPiPos[i_eta]);
			fOutputList->Add(hNsigmaKPos[i_eta]);
			fOutputList->Add(hNsigmaPPos[i_eta]);
			fOutputList->Add(hNsigmaPiNeg[i_eta]);
			fOutputList->Add(hNsigmaKNeg[i_eta]);
			fOutputList->Add(hNsigmaPNeg[i_eta]);
			fOutputList->Add(hPtTPCEtaPos[i_eta]);
			fOutputList->Add(hPtTPCEtaNeg[i_eta]);
			fOutputList->Add(hBetaPos[i_eta]);
			fOutputList->Add(hMomentumTOFEtaPos[i_eta]);
			fOutputList->Add(hPtTOFEtaPos[i_eta]);
			fOutputList->Add(hBetaNeg[i_eta]);
			fOutputList->Add(hMomentumTOFEtaNeg[i_eta]);
			fOutputList->Add(hPtTOFEtaNeg[i_eta]);
			fOutputList->Add(hdEdx[i_eta]);
			fOutputList->Add(hPtrTPC[i_eta]);
			fOutputList->Add(hPtVsP[i_eta]);
			fOutputList->Add(histPiV0[i_eta]);
			fOutputList->Add(histPV0[i_eta]);
			fOutputList->Add(histPiTof[i_eta]);
			fOutputList->Add(histPiTof2[i_eta]);
			fOutputList->Add(histEV0[i_eta]);
		}
	}

	/* fEventCuts.AddQAplotsToList(fOutputList); */
	PostData(1, fOutputList); // postdata will notify the analysis manager of
				  // changes / updates to the
}
//_____________________________________________________________________________
void AliAnalysisTaskSystematicsV0M::UserExec(Option_t *) {

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
		/* PostData(1, fOutputList); */
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
	float v0multalice = fMultSelection->GetEstimator("V0M")->GetValue();

	// Analyze V0s for the MB sample
	AnalyzeV0s();

	// This condition is added to reject the very few low-multiplicity events V0M Percentile 70-100%
	// that have a very large V0M Amplitude
	if (fv0mpercentile >= 70.0 && fv0mpercentile < 100.0) {
		if (v0multalice >= 400.0) { return; }
	}

	/* fMidRapidityMult = GetMidRapidityMultiplicity(); */
	/* fFlatTPC = GetFlatenicityTPC(); */ 
	fFlat = GetFlatenicityV0();

	if (fFlat >= 0.0) {

		hFlat->Fill(fv0mpercentile,fFlat);
		// piKp as a function of Flattenicity
		MakePIDanalysis();
		// Charged particle spectra as a function of Flattenicity
		/* MakeDataanalysis(); */
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
void AliAnalysisTaskSystematicsV0M::Terminate(Option_t *) {}
//______________________________________________________________________________
void AliAnalysisTaskSystematicsV0M::MakePIDanalysis() 
{
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

		float nSigmaPi = fPIDResponse->NumberOfSigmasTPC(esdTrack,AliPID::kPion);
		float nSigmaK = fPIDResponse->NumberOfSigmasTPC(esdTrack,AliPID::kKaon);
		float nSigmaP = fPIDResponse->NumberOfSigmasTPC(esdTrack,AliPID::kProton);
		float nSigmaPiTOF = fPIDResponse->NumberOfSigmasTOF(esdTrack,AliPID::kPion);

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

		double maxDCAxy = 10.0;
		maxDCAxy = fcutDCAxy->Eval(pt);

		if (TMath::Abs(DCAxy) > maxDCAxy )
			continue;

		if( charge < 0.0){
			hNsigmaPiNeg[nh]->Fill(pt,nSigmaPi,fv0mpercentile);
			hNsigmaKNeg[nh]->Fill(pt,nSigmaK,fv0mpercentile);
			hNsigmaPNeg[nh]->Fill(pt,nSigmaP,fv0mpercentile);
			hPtTPCEtaNeg[nh]->Fill(pt,fv0mpercentile);
		}else{
			hNsigmaPiPos[nh]->Fill(pt,nSigmaPi,fv0mpercentile);
			hNsigmaKPos[nh]->Fill(pt,nSigmaK,fv0mpercentile);
			hNsigmaPPos[nh]->Fill(pt,nSigmaP,fv0mpercentile);
			hPtTPCEtaPos[nh]->Fill(pt,fv0mpercentile);
		} 

		if( TOFPID(esdTrack) ){

			double trkLength = esdTrack->GetIntegratedLength();
			double beta = trkLength/((esdTrack->GetTOFsignal()-fPIDResponse->GetTOFResponse().GetStartTime(esdTrack->P()))*C_Value);

			if(esdTrack->Charge() < 0.0){
				hBetaNeg[nh]->Fill(momentum,beta,fv0mpercentile);
				hMomentumTOFEtaNeg[nh]->Fill(momentum,fv0mpercentile);
				hPtTOFEtaNeg[nh]->Fill(pt,fv0mpercentile);
			}else{ 
				hBetaPos[nh]->Fill(momentum,beta,fv0mpercentile);
				hMomentumTOFEtaPos[nh]->Fill(momentum,fv0mpercentile);
				hPtTOFEtaPos[nh]->Fill(pt,fv0mpercentile);
			}
		}

		hPtVsP[nh]->Fill(momentum,pt);

		if(!PhiCut(esdTrack->Pt(), phi, esdTrack->Charge(), 1.0, fcutLow, fcutHigh))
			continue;

		if (Ncl < fNcl)
			continue;

		if (fdEdxCalibrated) { dEdx *= (50.0/EtaCalibration(eta)); }

		hdEdx[nh]->Fill(momentum,dEdx,fv0mpercentile);
		hPtrTPC[nh]->Fill(pt,fv0mpercentile);

		if ( TOFPID(esdTrack) && (TMath::Abs(nSigmaPiTOF) < 2.0)) { histPiTof[nh]->Fill(momentum,dEdx); }

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

		if(beta>1.0){
			histPiTof2[nh]->Fill(momentum,dEdx);
		}

		if ( momentum <= 0.6 && momentum >= 0.4 ){ //only p:0.4-0.6 GeV, pion MIP
			if( dEdx < 60.0 && dEdx > 40.0 ){
				/* hMIPVsEta->Fill(eta,dEdx); */
				pMIPVsEta->Fill(eta,dEdx);
			}
			if( dEdx > 70.0 && dEdx < 90.0 ){
				if(TMath::Abs(beta-1)<0.1){
					/* hPlateauVsEta->Fill(eta,dEdx); */
					pPlateauVsEta->Fill(eta,dEdx);
				}
			}
		}


	}//end of track loop


}
//________________________________________________________________________
void AliAnalysisTaskSystematicsV0M::AnalyzeV0s()
{

	Int_t nv0s = fESD->GetNumberOfV0s();

	const AliESDVertex *myBestPrimaryVertex = fESD->GetPrimaryVertex();

	if ( !myBestPrimaryVertex ) 
		return;
	if ( !(myBestPrimaryVertex->GetStatus()) )
		return;

	Double_t  lPrimaryVtxPosition[3];
	myBestPrimaryVertex->GetXYZ(lPrimaryVtxPosition);
	Double_t  lPrimaryVtxCov[6];
	myBestPrimaryVertex->GetCovMatrix(lPrimaryVtxCov);
	Double_t  lPrimaryVtxChi2 = myBestPrimaryVertex->GetChi2toNDF();

	AliAODVertex* myPrimaryVertex = new AliAODVertex(lPrimaryVtxPosition, lPrimaryVtxCov, lPrimaryVtxChi2, NULL, -1, AliAODVertex::kPrimary);

	//
	// LOOP OVER V0s, K0s, L, AL
	//

	for (Int_t iV0=0; iV0<nv0s; iV0++)
	{

		AliESDv0 *esdV0 = fESD->GetV0(iV0);
		if (!esdV0) continue;

		//check onfly status
		//		if( !esdV0->GetOnFlyStatus() )
		//			continue;

		if (esdV0->GetOnFlyStatus()!=0)
			continue;

		// AliESDTrack (V0 Daughters)
		UInt_t lKeyPos = (UInt_t)TMath::Abs(esdV0->GetPindex());
		UInt_t lKeyNeg = (UInt_t)TMath::Abs(esdV0->GetNindex());

		AliESDtrack *pTrack = fESD->GetTrack(lKeyPos);
		AliESDtrack *nTrack = fESD->GetTrack(lKeyNeg);

		if (!pTrack || !nTrack) {
			Printf("ERROR: Could not retreive one of the daughter track");
			continue;
		}

		// Remove like-sign
		if (pTrack->GetSign() == nTrack->GetSign())
			continue;

		// Eta cut on decay products
		if(TMath::Abs(pTrack->Eta()) > fEtaCut || TMath::Abs(nTrack->Eta()) > fEtaCut)
			continue;

		if (!fTrackFilterPID->IsSelected(pTrack)) { continue; }
		if (!fTrackFilterPID->IsSelected(nTrack)) { continue; }

		// Check if switch does anything!
		Bool_t isSwitched = kFALSE;
		if (pTrack->GetSign() < 0) { // switch
			isSwitched = kTRUE;
			AliESDtrack* helpTrack = nTrack;
			nTrack = pTrack;
			pTrack = helpTrack;
		}

		AliKFVertex primaryVtxKF( *myPrimaryVertex );
		AliKFParticle::SetField(fESD->GetMagneticField());

		// Also implement switch here!!!!!!
		AliKFParticle* negEKF  = 0; // e-
		AliKFParticle* posEKF  = 0; // e+
		AliKFParticle* negPiKF = 0; // pi -
		AliKFParticle* posPiKF = 0; // pi +
		AliKFParticle* posPKF  = 0; // p
		AliKFParticle* negAPKF = 0; // p-bar

		if(!isSwitched) {
			negEKF  = new AliKFParticle( *(esdV0->GetParamN()) , 11);
			posEKF  = new AliKFParticle( *(esdV0->GetParamP()) ,-11);
			negPiKF = new AliKFParticle( *(esdV0->GetParamN()) ,-211);
			posPiKF = new AliKFParticle( *(esdV0->GetParamP()) , 211);
			posPKF  = new AliKFParticle( *(esdV0->GetParamP()) , 2212);
			negAPKF = new AliKFParticle( *(esdV0->GetParamN()) ,-2212);
		} else { // switch + and -
			negEKF  = new AliKFParticle( *(esdV0->GetParamP()) , 11);
			posEKF  = new AliKFParticle( *(esdV0->GetParamN()) ,-11);
			negPiKF = new AliKFParticle( *(esdV0->GetParamP()) ,-211);
			posPiKF = new AliKFParticle( *(esdV0->GetParamN()) , 211);
			posPKF  = new AliKFParticle( *(esdV0->GetParamN()) , 2212);
			negAPKF = new AliKFParticle( *(esdV0->GetParamP()) ,-2212);
		}

		AliKFParticle v0GKF;  // Gamma e.g. from pi0
		v0GKF+=(*negEKF);
		v0GKF+=(*posEKF);
		v0GKF.SetProductionVertex(primaryVtxKF);

		AliKFParticle v0K0sKF; // K0 short
		v0K0sKF+=(*negPiKF);
		v0K0sKF+=(*posPiKF);
		v0K0sKF.SetProductionVertex(primaryVtxKF);

		AliKFParticle v0LambdaKF; // Lambda
		v0LambdaKF+=(*negPiKF);
		v0LambdaKF+=(*posPKF);
		v0LambdaKF.SetProductionVertex(primaryVtxKF);

		AliKFParticle v0AntiLambdaKF; // Lambda-bar
		v0AntiLambdaKF+=(*posPiKF);
		v0AntiLambdaKF+=(*negAPKF);
		v0AntiLambdaKF.SetProductionVertex(primaryVtxKF);

		Double_t dmassG     = TMath::Abs(v0GKF.GetMass());
		Double_t dmassK     = TMath::Abs(v0K0sKF.GetMass()-0.498);
		Double_t dmassL     = TMath::Abs(v0LambdaKF.GetMass()-1.116);
		Double_t dmassAL    = TMath::Abs(v0AntiLambdaKF.GetMass()-1.116);

		if( dmassG  > 0.1 &&   
				dmassK  > 0.1 &&   
				dmassL  > 0.1 &&   
				dmassAL > 0.1
		  )    
			continue;

		for( Int_t case_v0 = 0; case_v0 < 2; ++case_v0 ){

			switch(case_v0){
				case 0:{

					       Bool_t fillPos = kFALSE;
					       Bool_t fillNeg = kFALSE;

					       if(dmassG < 0.1)
						       continue;

					       if(dmassK>0.01 && dmassL>0.01 && dmassAL>0.01){
						       continue;
					       }

					       if(dmassL<0.01){
						       fillPos = kTRUE;
					       }

					       if(dmassAL<0.01) {
						       if(fillPos)
							       continue;
						       fillNeg = kTRUE;
					       }

					       if(dmassK<0.01) {
						       if(fillPos||fillNeg)
							       continue;
						       fillPos = kTRUE;
						       fillNeg = kTRUE;
					       }

					       for(Int_t j = 0; j < 2; j++) {

						       AliESDtrack* track = 0;

						       if(j==0) {

							       if(fillNeg)
								       track = nTrack;
							       else
								       continue;
						       } else {

							       if(fillPos)
								       track = pTrack;
							       else
								       continue;
						       }

						       if(track->GetTPCsignalN()<fNcl) { continue; }
						       Double_t phi = track->Phi();

						       /* if(!PhiCut(track->Pt(), phi, track->Charge(), Magf, fcutLow, fcutHigh)) */
						       if(!PhiCut(track->Pt(), phi, track->Charge(), 1.0, fcutLow, fcutHigh)) { continue; }

						       Double_t eta      = track->Eta();
						       Double_t momentum = track->P();
						       Double_t dedx     = track->GetTPCsignal();
						       Double_t dedxUnc  = track->GetTPCsignal();

						       if(fdEdxCalibrated) { dedx *= (50.0/EtaCalibration(eta)); }

						       if(fillPos&&fillNeg){
							       if(dedxUnc < 60.0 && dedxUnc > 40.0){
								       if(momentum<0.6&&momentum>0.4){
									       pMIPVsEtaV0s->Fill(eta,dedx);
								       }
							       }
						       }

						       Int_t nh = -1;

						       if(TMath::Abs(eta)<0.2)
							       nh = 0; 
						       else if(TMath::Abs(eta)>=0.2 && TMath::Abs(eta)<0.4)
							       nh = 1; 
						       else if(TMath::Abs(eta)>=0.4 && TMath::Abs(eta)<0.6)
							       nh = 2; 
						       else if(TMath::Abs(eta)>=0.6 && TMath::Abs(eta)<0.8)
							       nh = 3; 

						       if(nh<0) { continue; }

						       if(fillPos&&fillNeg){
							       float nSigmaPiTPC = fPIDResponse->NumberOfSigmasTPC(track,AliPID::kPion);
							       float nSigmaPiTOF = fPIDResponse->NumberOfSigmasTOF(track,AliPID::kPion);
							       if ((momentum < 1.0) && (TMath::Abs(nSigmaPiTPC) < 2.0)) { histPiV0[nh]->Fill(momentum, dedx); }
							       if (!TOFPID(track)) { continue; }
							       if ((momentum >= 1.0 ) && (TMath::Abs(nSigmaPiTOF) < 2.0)) { histPiV0[nh]->Fill(momentum, dedx); }
						       }
						       else{
							       float nSigmaPTPC = fPIDResponse->NumberOfSigmasTPC(track,AliPID::kProton);
							       float nSigmaPTOF = fPIDResponse->NumberOfSigmasTOF(track,AliPID::kProton);
							       if ((momentum < 1.0) && (TMath::Abs(nSigmaPTPC) < 2.5)) { histPV0[nh]->Fill(momentum, dedx); }
							       if (!TOFPID(track)) { continue; }
							       if ((momentum >= 1.0) && (TMath::Abs(nSigmaPTOF) < 2.5)) { histPV0[nh]->Fill(momentum, dedx); }
						       }
					       }//end loop over two tracks

				       };
				       break;

				case 1:{//gammas

					       Bool_t fillPos = kFALSE;
					       Bool_t fillNeg = kFALSE;

					       if( dmassK>0.01 && dmassL>0.01 && dmassAL>0.01 ) {
						       if( dmassG<0.01 && dmassG>0.0001 ) {
							       if( TMath::Abs(nTrack->GetTPCsignal()-EtaCalibrationEl(nTrack->Eta())) < 5.0)
								       fillPos = kTRUE;
						       } else {
							       continue;
						       }
					       }

					       if(fillPos == kTRUE && fillNeg == kTRUE)
						       continue;

					       AliESDtrack* track = 0;
					       if(fillNeg)
						       track = nTrack;
					       else if(fillPos)
						       track = pTrack;
					       else
						       continue;

					       Double_t dedx     = track->GetTPCsignal();
					       Double_t eta      = track->Eta();
					       Double_t phi      = track->Phi();
					       Double_t momentum = track->P();

					       if(fdEdxCalibrated) { dedx *= (50.0/EtaCalibration(eta)); }

					       if(track->GetTPCsignalN()<fNcl) { continue; }

					       if(!PhiCut(track->Pt(), phi, track->Charge(), 1.0, fcutLow, fcutHigh)) { continue; }

					       Int_t nh = -1;

					       if(TMath::Abs(eta)<0.2)
						       nh = 0; 
					       else if(TMath::Abs(eta)>=0.2 && TMath::Abs(eta)<0.4)
						       nh = 1; 
					       else if(TMath::Abs(eta)>=0.4 && TMath::Abs(eta)<0.6)
						       nh = 2; 
					       else if(TMath::Abs(eta)>=0.6 && TMath::Abs(eta)<0.8)
						       nh = 3; 

					       if(nh<0) { continue; }

					       histEV0[nh]->Fill(momentum, dedx);

				       };
				       break;


			}//end switch

		}//end loop over case V0

		// clean up loop over v0

		delete negPiKF;
		delete posPiKF;
		delete posPKF;
		delete negAPKF;



	}

	delete myPrimaryVertex;

}
//______________________________________________________________________________
Double_t AliAnalysisTaskSystematicsV0M::GetFlatenicityTPC() {

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
Double_t AliAnalysisTaskSystematicsV0M::GetFlatenicityV0() {

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
		/* hActivityV0DataSectBefore->Fill(iCh, RhoLattice[iCh]); */
	}
	// after calibration
	/* if (fIsCalib) { */
	/* 	for (Int_t iCh = 0; iCh < nCells; iCh++) { */
	/* 		RhoLattice[iCh] *= fParVtx->Eval(0.0) / fParVtx->Eval(fVtxz); */
	/* 	} */
	/* } */

	// Filling histos with mult info
	/* float total_v0_tmp = 0.0; */
	/* for (Int_t iCh = 0; iCh < nCells; iCh++) { */
	/* 	hActivityV0DataSect->Fill(iCh, RhoLattice[iCh]); */
	/* total_v0_tmp += RhoLattice[iCh]; */
	/* } */
	/* float total_v0 = total_v0_tmp; */
	/* cout << "total_v0 = " << total_v0 << endl; */
	/* hV0vsVtxz->Fill(fVtxz, total_v0); */

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
Double_t AliAnalysisTaskSystematicsV0M::GetMidRapidityMultiplicity() {

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
Bool_t AliAnalysisTaskSystematicsV0M::HasRecVertex() {

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
//______________________________________________________________________________
Bool_t AliAnalysisTaskSystematicsV0M::TOFPID(AliESDtrack * track) 
{
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
//______________________________________________________________________________
Double_t AliAnalysisTaskSystematicsV0M::EtaCalibration(const Double_t &eta) 
{
	double calibration = 999.0;
	if (eta<0.0){ calibration = fEtaCalibrationNeg->Eval(eta); } 
	else{ calibration = fEtaCalibrationPos->Eval(eta); }

	return calibration;
}
//________________________________________________________________________
double AliAnalysisTaskSystematicsV0M::EtaCalibrationEl(const double &eta)
{
	double calibration = 999.0;
	if (eta<0.0){ calibration = felededxfitNeg->Eval(eta); }
	else{ calibration = felededxfitPos->Eval(eta); }

	return calibration;
}
//______________________________________________________________________________
Bool_t AliAnalysisTaskSystematicsV0M::PhiCut(Double_t pt, Double_t phi, Double_t q, Float_t mag, TF1* phiCutLow, TF1* phiCutHigh) {

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
//______________________________________________________________________________
Int_t AliAnalysisTaskSystematicsV0M::GetPidCode(Int_t pdgCode) {
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
