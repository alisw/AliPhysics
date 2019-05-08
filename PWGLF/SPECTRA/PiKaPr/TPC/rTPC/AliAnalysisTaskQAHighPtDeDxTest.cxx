/*************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, All rights reserved. *
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


#include "AliAnalysisTaskQAHighPtDeDxTest.h"

// ROOT includes
#include <TList.h>
#include <TTree.h>
#include <TMath.h>
#include <TH1.h>
#include <TF1.h>
#include <TProfile.h>
#include <TParticle.h>
#include <TFile.h>

// AliRoot includes
#include "AliAnalysisTask.h"
#include <AliAnalysisManager.h>
#include <AliAnalysisFilter.h>
#include <AliESDInputHandler.h>
#include <AliESDEvent.h>
#include <AliEventCuts.h>
#include <AliESDVertex.h>
#include <AliLog.h>
#include <AliExternalTrackParam.h>
#include <AliESDtrackCuts.h>
#include <AliESDVZERO.h>
#include <AliAODVZERO.h>
#include "AliMultSelection.h"
#include <AliPIDResponse.h>
#include "AliTPCPIDResponse.h"

#include <AliMCEventHandler.h>
#include <AliMCEvent.h>
#include <AliStack.h>

#include <TTreeStream.h>

#include <AliHeader.h>
#include <AliGenPythiaEventHeader.h>
#include <AliGenDPMjetEventHeader.h>

#include <AliCentrality.h>
#include <AliESDv0.h>
#include <AliKFVertex.h>
#include <AliAODVertex.h>

#include <AliAODTrack.h>
#include <AliVParticle.h>
#include <AliPID.h>
#include <AliAODPid.h>
#include <AliAODMCHeader.h>


// STL includes
#include <iostream>
using namespace std;


//
// Responsible:
// Antonio Ortiz (Lund)
// Peter Christiansen (Lund)
//




const Double_t AliAnalysisTaskQAHighPtDeDxTest::fgkClight = 2.99792458e-2;
namespace {
	Float_t magf = -1;
	TF1 cutLow ("StandardPhiCutLow",  "0.1/x/x+TMath::Pi()/18.0-0.025", 0, 50);
	TF1 cutHigh("StandardPhiCutHigh", "0.12/x+TMath::Pi()/18.0+0.035", 0, 50);
	Double_t DeDxMIPMin  = 40;
	Double_t DeDxMIPMax  = 60;
	const Int_t nHists = 4;
	Float_t centralityGlobal = -10;
	Int_t etaLow[nHists+5]  = {-8, -8, -6, -4, -2, 0, 2, 4, 6};
	Int_t etaHigh[nHists+5] = { 8, -6, -4, -2,  0, 2, 4, 6, 8};

	Int_t EtaLow[nHists]  = {6, 4, 2, 0};
	Int_t EtaHigh[nHists] = {8, 6, 4, 2};
	const Int_t nCent = 11;
	const Double_t CentMin[nCent] = {0.0,5.0,10.0,20.0,30.0,40.0,50.0,60.0,70.0,80.0,0.0};
	const Double_t CentMax[nCent] = {5.0,10.0,20.0,30.0,40.0,50.0,60.0,70.0,80.0,90.0,90.0};
	const Char_t *CentName[nCent] = {"0.0-5.0","5.0-10.0","10.0-20.0","20.0-30.0","30.0-40.0","40.0-50.0","50.0-60.0","60.0-70.0","70.0-80.0","80.0-90.0","0.0-90.0"};

	const Double_t aPos[nCent] = {49.3854 ,49.2214 ,49.1335 ,49.1671 ,49.3608 ,49.6613 ,49.9794 ,50.2623 ,50.4512 ,50.5573 ,49.3285 };
	const Double_t bPos[nCent] = {8.50827 ,8.20077 ,8.33807 ,9.36101 ,10.857  ,12.5605 ,14.549  ,15.7691 ,16.6262 ,17.6903 ,9.49632 };
	const Double_t cPos[nCent] = {-150.134,-129.652,-109.397,-90.7946,-74.823 ,-63.4573,-61.5329,-59.7291,-58.2066,-68.7638,-107.865};
	const Double_t dPos[nCent] = {838.211 ,715.728 ,587.286 ,457.231 ,327.686 ,220.404 ,166.382 ,127.903 ,92.6146 ,149.173 ,556.093 };
	const Double_t ePos[nCent] = {-2277.88,-1949.1 ,-1603.4 ,-1256.84,-880.111,-553.732,-366.853,-236.343,-100.579,-269.8  ,-1504.6 };
	const Double_t fPos[nCent] = {3318.5  ,2846.65 ,2353.97 ,1879.3  ,1321.22 ,825.927 ,527.742 ,321.84  ,81.0856 ,358.708 ,2209.12 };
	const Double_t gPos[nCent] = {-2470.88,-2119.51,-1754.71,-1422.31,-1001.68,-623.129,-390.347,-231.087,-23.963 ,-258.275,-1651.18};
	const Double_t hPos[nCent] = {734.774 ,628.065 ,517.362 ,422.993 ,295.064 ,179.067 ,107.452 ,58.8568 ,-10.8729,68.8145 ,488.407 };

	const Double_t aNeg[nCent] = {49.3578 ,49.2119 ,49.1413  ,49.1946  ,49.4028 ,49.7114 ,50.037  ,50.3247 ,50.517  ,50.6306 ,49.3393 };
	const Double_t bNeg[nCent] = {-13.9589,-11.045 ,-8.12876 ,-5.51988 ,-4.00106,-3.40983,-3.89816,-4.07962,-4.46431,-4.66132,-8.31485};
	const Double_t cNeg[nCent] = {-202.352,-153.935,-100.674 ,-42.8637 ,4.86577 ,40.8121 ,57.3274 ,70.8767 ,75.2846 ,78.9678 ,-88.9968};
	const Double_t dNeg[nCent] = {-1100.28,-833.932,-532.883 ,-197.343 ,98.1673 ,334.917 ,459.629 ,563.165 ,602.338 ,640.548 ,-448.818};
	const Double_t eNeg[nCent] = {-2986.55,-2265.31,-1439.34 ,-519.437 ,322.639 ,1013.08 ,1385.47 ,1709.42 ,1825.07 ,1965.66 ,-1191.73};
	const Double_t fNeg[nCent] = {-4394.3 ,-3334.72,-2108.53 ,-756.468 ,518.829 ,1576.09 ,2144.36 ,2665.59 ,2829.14 ,3084.37 ,-1733   };
	const Double_t gNeg[nCent] = {-3344.02,-2535.33,-1589.49 ,-561.256 ,432.513 ,1260.8  ,1699.07 ,2123.11 ,2233.53 ,2463    ,-1300.55};
	const Double_t hNeg[nCent] = {-1030.35,-779.155,-481.989 ,-164.124 ,149.287 ,411.113 ,546.43  ,684.361 ,712.202 ,793.691 ,-392.685};

	const Double_t aPosEl[nCent]={78.3895 ,78.3782 ,78.3037 ,78.2321 ,78.2172 ,78.2288 ,78.2487 ,78.4117 ,78.5567,78.5476  ,78.3182 };
	const Double_t bPosEl[nCent]={2.22721 ,3.22037 ,4.70227 ,7.30353 ,10.3975 ,14.1058 ,17.8151 ,19.8804 ,20.5765,24.0961  ,5.69795 };
	const Double_t cPosEl[nCent]={-13.5734,-16.5663,-19.8482,-26.3319,-34.337 ,-43.3906,-53.6702,-57.9954,-59.2553,-77.6615,-22.2188};
	const Double_t dPosEl[nCent]={28.8237 ,30.1572 ,30.0755 ,32.5052 ,38.7437 ,44.5118 ,55.5141 ,55.75   ,54.6671 ,89.0958 ,31.7212 };
	const Double_t ePosEl[nCent]={-20.1399,-19.4987,-17.3879,-15.3834,-16.6006,-16.3666,-20.9339,-18.1111,-16.1267,-37.2048,-17.3347};

	const Double_t aNegEl[nCent]={78.3209 ,78.2647 ,78.1615 ,78.0598 ,78.0252 ,77.975  ,77.9915 ,78.0635 ,78.2788 ,78.3643 ,78.18   };
	const Double_t bNegEl[nCent]={-1.59833,-2.26827,-3.40339,-5.38275,-7.50926,-11.231 ,-14.6979,-17.5359,-17.1721,-17.7016,-4.24817};
	const Double_t cNegEl[nCent]={-7.94146,-9.39369,-11.9299,-16.8957,-21.2631,-30.3458,-40.7906,-49.9978,-46.2792,-48.5651,-14.0596};
	const Double_t dNegEl[nCent]={-18.6273,-16.6219,-15.7565,-17.7392,-17.8171,-23.2   ,-34.6427,-45.1412,-34.3996,-39.208 ,-18.0062};
	const Double_t eNegEl[nCent]={-14.633 ,-11.5233,-9.08855,-7.92521,-5.32888,-4.58932,-8.95339,-12.9314,-4.43948,-8.11577,-9.89118};

	/*
	const Double_t aPos[nCent] = {49.5068 ,49.3267 ,49.2128 ,49.2071 ,49.3462 ,49.5824 ,49.8366 ,50.0635 ,50.2113 ,50.2987 ,0};
	const Double_t bPos[nCent] = {7.99301 ,7.63037 ,7.68758 ,8.50097 ,9.76739 ,11.1737 ,12.8348 ,13.871  ,14.8113 ,15.4976 ,0};
	const Double_t cPos[nCent] = {-148.922,-129.292,-109.707,-91.5511,-76.4792,-65.543 ,-62.9378,-60.5559,-61.2258,-66.5928,0};
	const Double_t dPos[nCent] = {850.273 ,734.612 ,611.057 ,486.88  ,367.204 ,268.673 ,217.937 ,179.637 ,158.53  ,184.901 ,0};
	const Double_t ePos[nCent] = {-2342.24,-2033.14,-1697.72,-1366.61,-1019.28,-721.825,-553.408,-426.339,-331.665,-410.548,0};
	const Double_t fPos[nCent] = {3450.53 ,3007.15 ,2522.94 ,2066.81 ,1551.02 ,1099.95 ,838.608 ,638.621 ,458.864,590.143  ,0};
	const Double_t gPos[nCent] = {-2597.23,-2266.78,-1903.17,-1581.39,-1190.71,-845.194,-646.339,-490.554,-328.57,-440.673 ,0};
	const Double_t hPos[nCent] = {781.335 ,680.872 ,569.019 ,476.95  ,357.339 ,251.067 ,191.488 ,143.235 ,86.6676,125.1    ,0};

	const Double_t aNeg[nCent] = {49.4762 ,49.3119 ,49.2138 ,49.2249 ,49.379  ,49.6208,49.8805,50.1069 ,50.2655 ,50.3576 ,0};
	const Double_t bNeg[nCent] = {-14.4745,-11.7702,-9.06904,-6.63424,-5.1604 ,-4.5837,-4.8944,-5.20679,-5.26467,-5.56416,0};
	const Double_t cNeg[nCent] = {-213.498,-168.72 ,-119.345,-66.0901,-21.8998,10.161 ,26.3479,36.9805,44.5138  ,44.9557 ,0};
	const Double_t dNeg[nCent] = {-1178.99,-933.383,-653.758,-344.876,-71.3742,137.254,255.881,341.299,395.363  ,408.37  ,0};
	const Double_t eNeg[nCent] = {-3241.31,-2577.39,-1807.84,-959.262,-177.258,424.91 ,776.464,1046.49,1202.08  ,1255.38 ,0};
	const Double_t fNeg[nCent] = {-4823.81,-3848.86,-2702.57,-1450.74,-261.322,651.884,1187.16,1622.13,1846.67  ,1947.08 ,0};
	const Double_t gNeg[nCent] = {-3710.7 ,-2966.13,-2079.21,-1122.78,-191.756,517.224,930.416,1283.28,1443.46  ,1535.22 ,0};
	const Double_t hNeg[nCent] = {-1155.44,-923.834,-644.487,-347.314,-52.3973,169.789,297.806,412.004,456.702  ,489.474 ,0};

	const Double_t aPosEl[nCent]={78.321  ,78.3324 ,78.2748 ,78.2127 ,78.2036 ,78.1972 ,78.2077 ,78.2999 ,78.3799 ,78.3387 ,0};
	const Double_t bPosEl[nCent]={2.33678 ,2.9658  ,4.21018 ,6.49506 ,8.85515 ,12.0224 ,15.2225 ,17.0048 ,19.0335 ,21.0828 ,0};
	const Double_t cPosEl[nCent]={-12.4591,-13.8368,-16.3764,-21.9963,-27.3262,-34.4519,-44.0844,-47.555 ,-56.3171,-63.9511,0};
	const Double_t dPosEl[nCent]={27.4761 ,25.9919 ,25.2571 ,26.9013 ,29.3614 ,32.0544 ,43.7582 ,44.0543 ,56.5737 ,68.2497 ,0};
	const Double_t ePosEl[nCent]={-19.8433,-17.523 ,-15.3877,-13.1872,-12.4892,-10.3975,-16.0792,-14.1679,-19.9037,-26.3362,0};

	const Double_t aNegEl[nCent]={78.2558 ,78.2234 ,78.1366 ,78.0618 ,78.0306 ,77.9744 ,77.9784 ,78.002  ,78.1778 ,78.0184 ,0};
	const Double_t bNegEl[nCent]={-1.9557 ,-2.45987,-3.49541,-5.00473,-6.78782,-9.91944,-12.4015,-15.338 ,-14.1248,-18.3535,0};
	const Double_t cNegEl[nCent]={-7.73037,-8.50241,-11.232 ,-14.524 ,-18.3252,-26.0661,-31.9118,-42.5509,-33.8971,-52.0504,0};
	const Double_t dNegEl[nCent]={-18.288 ,-15.12  ,-15.3748,-14.806 ,-15.1126,-19.5903,-23.1564,-37.1449,-17.7372,-47.6752,0};
	const Double_t eNegEl[nCent]={-14.647 ,-10.9019,-9.44581,-6.85271,-4.90317,-4.08212,-3.96587,-10.2503,2.92615 ,-14.0647,0};
	*/

	Int_t nDeltaPiBins   = 80;
	Double_t deltaPiLow  = 20;
	Double_t deltaPiHigh = 100;
	const Double_t dEdxHigh = 200;
	const Double_t dEdxLow  = 40;
	const Char_t *Pid[7]={"Ch","Pion","Kaon","Proton","Electron","Muon","Oher"};
	const Char_t *Q[3]={"", "Neg", "Pos"};
}

ClassImp(AliAnalysisTaskQAHighPtDeDxTest)
//_____________________________________________________________________________
AliAnalysisTaskQAHighPtDeDxTest::AliAnalysisTaskQAHighPtDeDxTest():
		AliAnalysisTaskSE(),
		fESD(0x0),
		fAOD(0x0),
		fEventCuts(0x0),
		fMC(0x0),
		fMCStack(0x0),
		fMCArray(0x0),
		fPIDResponse(0x0),
		fTrackFilter2015PbPb(0x0),
		fTrackFilterGolden(0x0),
		fTrackFilterTPC(0x0),
		fCentEst("V0M"),
		fAnalysisType("ESD"),
		fAnalysisMC(kFALSE),
		fAnalysisPbPb(kFALSE),
		ftrigBit(0x0),
		fRandom(0x0),
		fPileUpRej(kFALSE),
		fEtaCut(0.9),
		cent(3),
		fMinCent(0.0),
		fMaxCent(100.0),
		fStoreMcIn(kFALSE),//
		fMcProcessType(-999),
		fTriggeredEventMB(-999),
		fVtxStatus(-999),
		fZvtx(-999),
		fZvtxMC(-999),
		fRun(-999),
		fEventId(-999),
		fListOfObjects(0),
		fVtxMC(0x0),
		fdEdxCalibrated(0x0),
		fMakePid(0x0),
		fcent(0x0),
		fcentAfterPrimaries(0x0),
		fcentAfterV0s(0x0),
		fEtaCalibrationNeg(0x0),
		fEtaCalibration(0x0),
		felededxfitPos(0x0),
		felededxfitNeg(0x0),
		fcutDCAxy(0x0)


{


	for(Int_t i = 0; i<nCent; ++i){


		// Histograms for PreCalibration

		hMIPVsEta[i]=0;
		pMIPVsEta[i]=0;
		hMIPVsEtaV0s[i]=0;
		pMIPVsEtaV0s[i]=0;
		hPlateauVsEta[i]=0;
		pPlateauVsEta[i]=0;
		hPhi[i]=0;

		// Histograms for PostCalibration

		hPtAll[i]=0;
		hPtAllPos[i]=0;
		hPtAllNeg[i]=0;

		hDCAxyVsPtPiNeg[i]=0;
		hDCAxyVsPtKNeg[i]=0;
		hDCAxyVsPtPNeg[i]=0;
		hDCAxyVsPtPiNegC[i]=0;
		hDCAxyVsPtKNegC[i]=0;
		hDCAxyVsPtPNegC[i]=0;
		hDCAxyVsPtPiPos[i]=0;
		hDCAxyVsPtKPos[i]=0;
		hDCAxyVsPtPPos[i]=0;
		hDCAxyVsPtPiPosC[i]=0;
		hDCAxyVsPtKPosC[i]=0;
		hDCAxyVsPtPPosC[i]=0;

		for(Int_t j=0; j<nHists; ++j){

			hPtPos[i][j]=0;//TH1D, Transverse momentum distribution  positive-charged hadrons
			hPtNeg[i][j]=0;//TH1D, Transverse momentum distribution  negative-charged hadrons
			hPtVsP[i][j]=0;//TH2D, Transverse momentum Vs momentum

			hMIPVsNch[i][j]=0;//TH2D, MIP vs Nch for different eta intervals
			pMIPVsNch[i][j]=0;//TProfile, MIP vs Nch for different eta intervals
			hMIPVsPhi[i][j]=0;//TH2D, MIP vs phi for different eta intervals
			pMIPVsPhi[i][j]=0;//TProfile, MIP vs phi for different eta intervals
			hPlateauVsPhi[i][j]=0;//TH2D, dE/dx vs Phi, electrons 0.4<p<0.6 GeV/c
			pPlateauVsPhi[i][j]=0;//TProfile, dE/dx vs Phi, electrons 0.4<p<0.6 GeV/c

			hDeDxVsP[i][j]=0;//TH2D, DeDx vs P

			hnSigmaPiPos[i][j]=0;//TH2D, nSigmas vs Pt Pions
			hnSigmaKPos[i][j]=0;//TH2D, nSigmas vs Pt Kaons
			hnSigmaPPos[i][j]=0;//TH2D, nSigmas vs Pt Protons

			hnSigmaPiNeg[i][j]=0;//TH2D, nSigmas vs Pt AntiPions
			hnSigmaKNeg[i][j]=0;//TH2D, nSigmas vs Pt AntiKaons
			hnSigmaPNeg[i][j]=0;//TH2D, nSigmas vs Pt AntiProtons

			histPiV0[i][j]=0;//TH2D, dE/dx vs p, pi id by V0s
			histpPiV0[i][j]=0;//TH1D, pi id by V0s
			histPV0[i][j]=0;// TH2D, dE/dx vs p, p id by V0s
			histpPV0[i][j]=0;// TH1D, p id by V0s
			histPiTof[i][j]=0;//TH2D, dE/dx vs p for a "clean" sample of pions, beta>1
			histpPiTof[i][j]=0;//TH1D, for a "clean" sample of pions, beta>1
			histEV0[i][j]=0;

		}

	}

	//default constructor
	for(Int_t cent=0;cent<nCent;++cent){
		for(Int_t pid=0;pid<7;++pid){
			hMcIn[cent][pid]=0;
			hMcOut[cent][pid]=0;
			hMcInNeg[cent][pid]=0;
			hMcInPos[cent][pid]=0;
			hMcOutNeg[cent][pid]=0;
			hMcOutPos[cent][pid]=0;
		}
	}

	for(Int_t cent=0;cent<nCent;++cent){
		for(Int_t pid=0;pid<7;++pid){
			for(Int_t q=0;q<3;++q){
				hDCApTPrim[cent][pid][q]=0;
				hDCApTWDec[cent][pid][q]=0;
				hDCApTMate[cent][pid][q]=0;

				hDCApTPrim2[cent][pid][q]=0;
				hDCApTWDec2[cent][pid][q]=0;
				hDCApTMate2[cent][pid][q]=0;
			}
		}
	}




}


AliAnalysisTaskQAHighPtDeDxTest::AliAnalysisTaskQAHighPtDeDxTest(const char *name):
	AliAnalysisTaskSE(name),
	fESD(0x0),
	fAOD(0x0),
	fEventCuts(0x0),
	fMC(0x0),
	fMCStack(0x0),
	fMCArray(0x0),
	fPIDResponse(0x0),
	fTrackFilter2015PbPb(0x0),
	fTrackFilterGolden(0x0),
	fTrackFilterTPC(0x0),
	fCentEst("V0M"),
	fAnalysisType("ESD"),
	fAnalysisMC(kFALSE),
	fAnalysisPbPb(kFALSE),
	ftrigBit(0x0),
	fRandom(0x0),
	fPileUpRej(kFALSE),
	fEtaCut(0.9),
	cent(3),
	fMinCent(0.0),
	fMaxCent(100.0),
	fStoreMcIn(kFALSE),
	fMcProcessType(-999),
	fTriggeredEventMB(-999),
	fVtxStatus(-999),
	fZvtx(-999),
	fZvtxMC(-999),
	fRun(-999),
	fEventId(-999),
	fListOfObjects(0), 
	fVtxMC(0x0),
	fdEdxCalibrated(0x0),
	fMakePid(0x0),
	fcent(0x0),
	fcentAfterPrimaries(0x0),
	fcentAfterV0s(0x0),
	fEtaCalibrationNeg(0x0),
	fEtaCalibration(0x0),
	felededxfitPos(0x0),
	felededxfitNeg(0x0),
	fcutDCAxy(0x0)

{


	for(Int_t i = 0; i<nCent; ++i){

		hMIPVsEta[i] = 0;
		hMIPVsEta[i]=0;
		pMIPVsEta[i]=0;
		hMIPVsEtaV0s[i]=0;
		pMIPVsEtaV0s[i]=0;
		hPlateauVsEta[i]=0;
		pPlateauVsEta[i]=0;
		hPhi[i]=0;

		hPtAll[i]=0;
		hPtAllPos[i]=0;
		hPtAllNeg[i]=0;

		hDCAxyVsPtPiNeg[i]=0;
		hDCAxyVsPtKNeg[i]=0;
		hDCAxyVsPtPNeg[i]=0;
		hDCAxyVsPtPiNegC[i]=0;
		hDCAxyVsPtKNegC[i]=0;
		hDCAxyVsPtPNegC[i]=0;
		hDCAxyVsPtPiPos[i]=0;
		hDCAxyVsPtKPos[i]=0;
		hDCAxyVsPtPPos[i]=0;
		hDCAxyVsPtPiPosC[i]=0;
		hDCAxyVsPtKPosC[i]=0;
		hDCAxyVsPtPPosC[i]=0;

		for(Int_t j=0; j<nHists; ++j){

			hPtPos[i][j]=0;//TH1D, Transverse momentum distribution  positive-charged hadrons
			hPtNeg[i][j]=0;//TH1D, Transverse momentum distribution  negative-charged hadrons
			hPtVsP[i][j]=0;//TH2D, Transverse momentum Vs momentum

			hMIPVsNch[i][j]=0;//TH2D, MIP vs Nch for different eta intervals
			pMIPVsNch[i][j]=0;//TProfile, MIP vs Nch for different eta intervals
			hMIPVsPhi[i][j]=0;//TH2D, MIP vs phi for different eta intervals
			pMIPVsPhi[i][j]=0;//TProfile, MIP vs phi for different eta intervals
			hPlateauVsPhi[i][j]=0;//TH2D, dE/dx vs Phi, electrons 0.4<p<0.6 GeV/c
			pPlateauVsPhi[i][j]=0;//TProfile, dE/dx vs Phi, electrons 0.4<p<0.6 GeV/c

			hDeDxVsP[i][j]=0;//TH2D, DeDx vs P

			hnSigmaPiPos[i][j]=0;//TH2D, nSigmas vs Pt Pions
			hnSigmaKPos[i][j]=0;//TH2D, nSigmas vs Pt Kaons
			hnSigmaPPos[i][j]=0;//TH2D, nSigmas vs Pt Protons

			hnSigmaPiNeg[i][j]=0;//TH2D, nSigmas vs Pt AntiPions
			hnSigmaKNeg[i][j]=0;//TH2D, nSigmas vs Pt AntiKaons
			hnSigmaPNeg[i][j]=0;//TH2D, nSigmas vs Pt AntiProtons

			histPiV0[i][j]=0;//TH2D, dE/dx vs p, pi id by V0s
			histpPiV0[i][j]=0;//TH1D, pi id by V0s
			histPV0[i][j]=0;// TH2D, dE/dx vs p, p id by V0s
			histpPV0[i][j]=0;// TH1D, p id by V0s
			histPiTof[i][j]=0;//TH2D, dE/dx vs p for a "clean" sample of pions, beta>1
			histpPiTof[i][j]=0;//TH1D, for a "clean" sample of pions, beta>1
			histEV0[i][j]=0;

		}

	}

	// Default constructor (should not be used)
	for(Int_t cent=0; cent<nCent; ++cent){
		for(Int_t pid=0; pid<7; ++pid){
			hMcIn[cent][pid]=0;
			hMcOut[cent][pid]=0;
			hMcInNeg[cent][pid]=0;
			hMcInPos[cent][pid]=0;
			hMcOutNeg[cent][pid]=0;
			hMcOutPos[cent][pid]=0;
		}
	}

	for(Int_t cent=0;cent<nCent;++cent){
		for(Int_t pid=0;pid<7;++pid){
			for(Int_t q=0;q<3;++q){
				hDCApTPrim[cent][pid][q]=0;
				hDCApTWDec[cent][pid][q]=0;
				hDCApTMate[cent][pid][q]=0;

				hDCApTPrim2[cent][pid][q]=0;
				hDCApTWDec2[cent][pid][q]=0;
				hDCApTMate2[cent][pid][q]=0;
			}
		}
	}


	DefineOutput(1, TList::Class());//esto es nuevo
}




AliAnalysisTaskQAHighPtDeDxTest::~AliAnalysisTaskQAHighPtDeDxTest() {
	//
	// Destructor
	//

}



//______________________________________________________________________________



void AliAnalysisTaskQAHighPtDeDxTest::UserCreateOutputObjects()
{
	// This method is called once per worker node
	// Here we define the output: histograms and debug tree if requested
	// We also create the random generator here so it might get different seeds...
	fRandom = new TRandom(0); // 0 means random seed


	AliAnalysisManager *man = AliAnalysisManager::GetAnalysisManager();
	if(man){
		AliInputEventHandler* inputHandler = (AliInputEventHandler*)(man->GetInputEventHandler());
		if(inputHandler)fPIDResponse = inputHandler->GetPIDResponse();
	}


	//OpenFile(1);
	fListOfObjects = new TList();
	fListOfObjects->SetOwner();

	//
	// Histograms
	//

	fcent=new TH1F("fcent","fcent",13,0,13);
	fcentAfterPrimaries =new TH1F("fcentAfterPrimaries","fcentAfterPrimaries",13,0,13);
	fcentAfterV0s =new TH1F("fcentAfterV0s","fcentAfterV0s",13,0,13);
	fListOfObjects->Add(fcent);
	fListOfObjects->Add(fcentAfterPrimaries);
	fListOfObjects->Add(fcentAfterV0s);

	const Int_t nPtBins = 63;
	Double_t ptBins[nPtBins+1] = {

		0.01,0.1,0.12,0.14,0.16,0.18,0.2,0.25,0.3,0.35,0.4,0.45,0.5,
		0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1,1.1,1.2,1.3,1.4,
		1.5,1.6,1.7,1.8,1.9,2,2.2,2.4,2.6,2.8,3,3.2,3.4,3.6,3.8,4,
		4.5,5,5.5,6,6.5,7,8,9,10,11,12,13,14,15,16,18,20,22,24,26,30

	};


	const Int_t nPtBinsV0s = 25;
	Double_t ptBinsV0s[nPtBinsV0s+1] = { 0.0 , 0.1 , 0.2 , 0.3 , 0.4 , 0.5 , 0.6 , 0.7 , 0.8 , 0.9 , 1.0 ,
		1.2 , 1.4 , 1.6 , 1.8 , 2.0 , 2.5 , 3.0 , 3.5 , 4.0 , 5.0 , 7.0 ,
		9.0 , 12.0, 15.0, 20.0 };

	const Int_t ndcaBins = 100;
	Double_t dcaBins[ndcaBins+1] = {
		-4.0, -3.9, -3.8, -3.7, -3.6, -3.5, -3.4, -3.3, -3.2, -3.1,
		-3.0, -2.9, -2.8, -2.7, -2.6, -2.5, -2.4, -2.3, -2.2, -2.1,
		-2.0, -1.9, -1.8, -1.7, -1.6, -1.5, -1.4, -1.3, -1.2, -1.1,
		-1.0, -0.95, -0.9, -0.85, -0.8, -0.75, -0.7, -0.65, -0.6,
		-0.55, -0.5, -0.45, -0.4, -0.35, -0.3, -0.25, -0.2, -0.15,
		-0.1, -0.05, 0.0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4,
		0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95,
		1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2,
		2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3.0, 3.1, 3.2, 3.3, 3.4,
		3.5, 3.6, 3.7, 3.8, 3.9, 4.0};


	const Char_t* ending[nHists] = {"02", "24", "46", "68"};
	const Char_t* LatexEta[nHists] = {"|#eta|<0.2", "0.2<|#eta|<0.4", "0.4<|#eta|<0.6", "0.6<|#eta|<0.8" };

	fcutDCAxy = new TF1("fMaxDCAxy","[0]+[1]/(x^[2])",0,1e10);
	fcutDCAxy->SetParameter(0,0.0105);
	fcutDCAxy->SetParameter(1,0.0350);
	fcutDCAxy->SetParameter(2,1.1);


	fEtaCalibrationNeg = new TF1("fDeDxVsEtaNeg", "pol7", -1.0, 0.0);
	fEtaCalibration    = new TF1("fDeDxVsEtaPos", "pol7", 0.0, 1.0);

	felededxfitPos     = new TF1("felededxfitPos", "pol4", 0.0, 1.0);
	felededxfitNeg     = new TF1("felededxfitNeg", "pol4", -1.0, 0.0);


	Int_t nPhiBins = 36;

	for(Int_t i = 0; i<nCent; ++i){

		hMIPVsEta[i] = new TH2D(Form("hMIPVsEta%.2f-%.2f",CentMin[i],CentMax[i]),"; #eta; dE/dx_{MIP, primary tracks}",50,-0.8,0.8,DeDxMIPMax-DeDxMIPMin, DeDxMIPMin, DeDxMIPMax);
		pMIPVsEta[i] = new TProfile(Form("pMIPVsEta%.2f-%.2f",CentMin[i],CentMax[i]),"; #eta; #LT dE/dx #GT_{MIP, primary tracks}",50,-0.8,0.8, DeDxMIPMin, DeDxMIPMax);
		hMIPVsEtaV0s[i] = new TH2D(Form("hMIPVsEtaV0s%.2f-%.2f",CentMin[i],CentMax[i]),"; #eta; dE/dx_{MIP, secondary tracks}",50,-0.8,0.8,DeDxMIPMax-DeDxMIPMin, DeDxMIPMin, DeDxMIPMax);
		pMIPVsEtaV0s[i] = new TProfile(Form("pMIPVsEtaV0s%.2f-%.2f",CentMin[i],CentMax[i]),"; #eta; #LT dE/dx #GT_{MIP, secondary tracks}",50,-0.8,0.8,DeDxMIPMin, DeDxMIPMax);
		hPlateauVsEta[i] = new TH2D(Form("hPlateauVsEta%.2f-%.2f",CentMin[i],CentMax[i]),"; #eta; dE/dx_{Plateau, primary tracks}",50,-0.8,0.8,50, 60, 110);
		pPlateauVsEta[i] = new TProfile(Form("pPlateauVsEta%.2f-%.2f",CentMin[i],CentMax[i]),"; #eta; #LT dE/dx #GT_{Plateau, primary tracks}",50,-0.8,0.8, 60, 110);
		hPhi[i] = new TH2D(Form("histPhi%.2f-%.2f",CentMin[i],CentMax[i]), ";pt; #phi'", nPtBinsV0s, ptBinsV0s, 90, -0.05, 0.4);

		hPtAll[i] = new TH1D(Form("hPt%.2f-%.2f",CentMin[i],CentMax[i]),";#it{p}_{T};Counts",nPtBins,ptBins);
		hPtAll[i]->Sumw2();
		hPtAllNeg[i] = new TH1D(Form("hPtNeg%.2f-%.2f",CentMin[i],CentMax[i]),";#it{p}_{T};Counts",nPtBins,ptBins);
		hPtAllNeg[i]->Sumw2();
		hPtAllPos[i] = new TH1D(Form("hPtPos%.2f-%.2f",CentMin[i],CentMax[i]),";#it{p}_{T};Counts",nPtBins,ptBins);
		hPtAllPos[i]->Sumw2();

		hDCAxyVsPtPiNeg[i] = new TH2D(Form("hDCAxyVsPtPiNeg%.2f-%.2f",CentMin[i],CentMax[i]), "hDCAxyVsPtPiNeg; #it{p}_{T} (GeV/#it{c}); DCA_{xy} (cm)", nPtBins, ptBins, ndcaBins, dcaBins );
		hDCAxyVsPtPiNegC[i] = new TH2D(Form("hDCAxyVsPtPiNegC%.2f-%.2f",CentMin[i],CentMax[i]), "hDCAxyVsPtPiNeg Close; #it{p}_{T} (GeV/#it{c}); DCA_{xy} (cm)", nPtBins, ptBins, 800, -4.0, 4.0 );
		hDCAxyVsPtPiNeg[i]->Sumw2();
		hDCAxyVsPtPiNegC[i]->Sumw2();
		hDCAxyVsPtKNeg[i] = new TH2D(Form("hDCAxyVsPtKNeg%.2f-%.2f",CentMin[i],CentMax[i]), "hDCAxyVsPtKNeg; #it{p}_{T} (GeV/#it{c}); DCA_{xy} (cm)", nPtBins, ptBins, ndcaBins, dcaBins );
		hDCAxyVsPtKNegC[i] = new TH2D(Form("hDCAxyVsPtKNegC%.2f-%.2f",CentMin[i],CentMax[i]), "hDCAxyVsPtKNeg Close; #it{p}_{T} (GeV/#it{c}); DCA_{xy} (cm)", nPtBins, ptBins, 800, -4.0, 4.0 );
		hDCAxyVsPtKNeg[i]->Sumw2();
		hDCAxyVsPtKNegC[i]->Sumw2();
		hDCAxyVsPtPNeg[i] = new TH2D(Form("hDCAxyVsPtPNeg%.2f-%.2f",CentMin[i],CentMax[i]), "hDCAxyVsPtPNeg; #it{p}_{T} (GeV/#it{c}); DCA_{xy} (cm)", nPtBins, ptBins, ndcaBins, dcaBins );
		hDCAxyVsPtPNegC[i] = new TH2D(Form("hDCAxyVsPtPNegC%.2f-%.2f",CentMin[i],CentMax[i]), "hDCAxyVsPtPNeg Close; #it{p}_{T} (GeV/#it{c}); DCA_{xy} (cm)", nPtBins, ptBins, 800, -4.0, 4.0 );
		hDCAxyVsPtPNeg[i]->Sumw2();
		hDCAxyVsPtPNegC[i]->Sumw2();
		hDCAxyVsPtPiPos[i] = new TH2D(Form("hDCAxyVsPtPiPos%.2f-%.2f",CentMin[i],CentMax[i]), "hDCAxyVsPtPiPos; #it{p}_{T} (GeV/#it{c}); DCA_{xy} (cm)", nPtBins, ptBins, ndcaBins, dcaBins );
		hDCAxyVsPtPiPosC[i] = new TH2D(Form("hDCAxyVsPtPiPosC%.2f-%.2f",CentMin[i],CentMax[i]), "hDCAxyVsPtPiPos Close; #it{p}_{T} (GeV/#it{c}); DCA_{xy} (cm)", nPtBins, ptBins, 800, -4.0, 4.0 );
		hDCAxyVsPtPiPos[i]->Sumw2();
		hDCAxyVsPtPiPosC[i]->Sumw2();
		hDCAxyVsPtKPos[i] = new TH2D(Form("hDCAxyVsPtKPos%.2f-%.2f",CentMin[i],CentMax[i]), "hDCAxyVsPtKPos; #it{p}_{T} (GeV/#it{c}); DCA_{xy} (cm)", nPtBins, ptBins, ndcaBins, dcaBins );
		hDCAxyVsPtKPosC[i] = new TH2D(Form("hDCAxyVsPtKPosC%.2f-%.2f",CentMin[i],CentMax[i]), "hDCAxyVsPtKPos Close; #it{p}_{T} (GeV/#it{c}); DCA_{xy} (cm)", nPtBins, ptBins, 800, -4.0, 4.0 );
		hDCAxyVsPtKPos[i]->Sumw2();
		hDCAxyVsPtKPosC[i]->Sumw2();
		hDCAxyVsPtPPos[i] = new TH2D(Form("hDCAxyVsPtPPos%.2f-%.2f",CentMin[i],CentMax[i]), "hDCAxyVsPtPPos; #it{p}_{T} (GeV/#it{c}); DCA_{xy} (cm)", nPtBins, ptBins, ndcaBins, dcaBins );
		hDCAxyVsPtPPosC[i] = new TH2D(Form("hDCAxyVsPtPPosC%.2f-%.2f",CentMin[i],CentMax[i]), "hDCAxyVsPtPPos Close; #it{p}_{T} (GeV/#it{c}); DCA_{xy} (cm)", nPtBins, ptBins, 800, -4.0, 4.0 );
		hDCAxyVsPtPPos[i]->Sumw2();
		hDCAxyVsPtPPosC[i]->Sumw2();

		for(Int_t j=0; j<nHists; j++) {

			hDeDxVsP[i][j] = new TH2D( Form("hDeDxVsP%.2f-%.2f-%s",CentMin[i],CentMax[i],ending[j]), ";#it{p} [GeV/c]; dE/dx", nPtBins, ptBins, dEdxHigh-dEdxLow, dEdxLow, dEdxHigh);
			hDeDxVsP[i][j]->Sumw2();

			hnSigmaPiPos[i][j] = new TH2D(Form("hnSigmaPiPos%.2f-%.2f-%s",CentMin[i],CentMax[i],ending[j]), ";#it{p}_{T};nSigmaPiPos", nPtBins, ptBins, 40, -10, 10);
			hnSigmaPiPos[i][j]->Sumw2();

			hnSigmaKPos[i][j] = new TH2D(Form("hnSigmaKPos%.2f-%.2f-%s",CentMin[i],CentMax[i],ending[j]), ";#it{p}_{T};nSigmaKPos", nPtBins, ptBins, 40, -10, 10);
			hnSigmaKPos[i][j]->Sumw2();

			hnSigmaPPos[i][j] = new TH2D(Form("hnSigmaPPos%.2f-%.2f-%s",CentMin[i],CentMax[i],ending[j]), ";#it{p}_{T};nSigmaPPos", nPtBins, ptBins, 40, -10, 10);
			hnSigmaPPos[i][j]->Sumw2();

			hnSigmaPiNeg[i][j] = new TH2D(Form("hnSigmaPiNeg%.2f-%.2f-%s",CentMin[i],CentMax[i],ending[j]), ";#it{p}_{T};nSigmaPiNeg", nPtBins, ptBins, 40, -10, 10);
			hnSigmaPiNeg[i][j]->Sumw2();

			hnSigmaKNeg[i][j] = new TH2D(Form("hnSigmaKNeg%.2f-%.2f-%s",CentMin[i],CentMax[i],ending[j]), ";#it{p}_{T};nSigmaKNeg", nPtBins, ptBins, 40, -10, 10);
			hnSigmaKNeg[i][j]->Sumw2();

			hnSigmaPNeg[i][j] = new TH2D(Form("hnSigmaPNeg%.2f-%.2f-%s",CentMin[i],CentMax[i],ending[j]), ";#it{p}_{T};nSigmaPNeg", nPtBins, ptBins, 40, -10, 10);
			hnSigmaPNeg[i][j]->Sumw2();

			hMIPVsPhi[i][j] = new TH2D(Form("hMIPVsPhi%.2f-%.2f-%s",CentMin[i],CentMax[i],ending[j]), Form("%s; #phi (rad); dE/dx MIP",LatexEta[j]), nPhiBins, 0, 2*TMath::Pi(),DeDxMIPMax-DeDxMIPMin, DeDxMIPMin, DeDxMIPMax);
			hMIPVsPhi[i][j]->Sumw2();

			pMIPVsPhi[i][j] = new TProfile(Form("pMIPVsPhi%.2f-%.2f-%s",CentMin[i],CentMax[i],ending[j]), Form("%s; #phi (rad); dE/dx MIP",LatexEta[j]),  nPhiBins, 0, 2*TMath::Pi(),DeDxMIPMin, DeDxMIPMax);
			pMIPVsPhi[i][j]->Sumw2();

			hPlateauVsPhi[i][j]  = new TH2D(Form("hPlateauVsPhi%.2f-%.2f-%s",CentMin[i],CentMax[i],ending[j]), Form("%s; #phi (rad); dE/dx Plateau",LatexEta[j]),  nPhiBins, 0, 2*TMath::Pi(),20, 70, 90);
			hPlateauVsPhi[i][j]->Sumw2();

			pPlateauVsPhi[i][j] = new TProfile(Form("pPlateauVsPhi%.2f-%.2f-%s",CentMin[i],CentMax[i],ending[j]), Form("%s; #phi (rad); dE/dx Plateau",LatexEta[j]), nPhiBins, 0, 2*TMath::Pi(),DeDxMIPMax, 95);
			pPlateauVsPhi[i][j]->Sumw2();

			hMIPVsNch[i][j] = new TH2D(Form("hMIPVsNch%.2f-%.2f-%s",CentMin[i],CentMax[i],ending[j]), Form("%s; TPC track mult. |#eta|<0.8; dE/dx MIP",LatexEta[j]), 800, 1, 4001, DeDxMIPMax-DeDxMIPMin, DeDxMIPMin, DeDxMIPMax);
			hMIPVsNch[i][j]->Sumw2();

			pMIPVsNch[i][j] = new TProfile(Form("pMIPVsNch%.2f-%.2f-%s",CentMin[i],CentMax[i],ending[j]), Form("%s; TPC track mult. |#eta|<0.8; dE/dx MIP",LatexEta[j]), 800, 1, 4001, DeDxMIPMin, DeDxMIPMax);
			pMIPVsNch[i][j]->Sumw2();

			hPtPos[i][j] = new TH1D(Form("hPtPos%.2f-%.2f-%s",CentMin[i],CentMax[i],ending[j]),";#it{p}_{T};Counts",nPtBins,ptBins);
			hPtPos[i][j]->Sumw2();

			hPtNeg[i][j] = new TH1D(Form("hPtNeg%.2f-%.2f-%s",CentMin[i],CentMax[i],ending[j]),";#it{p}_{T};Counts",nPtBins,ptBins);
			hPtNeg[i][j]->Sumw2();

			hPtVsP[i][j] = new TH2D(Form("hPtVsP%.2f-%.2f-%s",CentMin[i],CentMax[i],ending[j]), ";#it{p} [GeV/c]; #it{p}_{T}", nPtBins, ptBins, nPtBins, ptBins);
			hPtVsP[i][j]->Sumw2();

			histPiV0[i][j]  = new TH2D(Form("hPiV0%.2f-%.2f-%s",CentMin[i],CentMax[i],ending[j]), "Pions id by V0", nPtBinsV0s, ptBinsV0s, nDeltaPiBins, deltaPiLow, deltaPiHigh);
			histPiV0[i][j]->Sumw2();

			histpPiV0[i][j]  = new TH1D(Form("hPiV01D%.2f-%.2f-%s",CentMin[i],CentMax[i],ending[j]), "Pions id by V0; #it{p} (GeV/#it{c}); counts", 200, 0, 20);
			histpPiV0[i][j]->Sumw2();

			histPV0[i][j]   = new TH2D(Form("hPV0%.2f-%.2f-%s",CentMin[i],CentMax[i],ending[j]), "Protons id by V0", nPtBinsV0s, ptBinsV0s, nDeltaPiBins, deltaPiLow, deltaPiHigh);
			histPV0[i][j]->Sumw2();

			histpPV0[i][j]  = new TH1D(Form("hPV01D%.2f-%.2f-%s",CentMin[i],CentMax[i],ending[j]), "Protons id by V0; #it{p} (GeV/#it{c}); counts", 200, 0, 20);
			histpPV0[i][j]->Sumw2();

			histPiTof[i][j] = new TH2D(Form("hPiTOF%.2f-%.2f-%s",CentMin[i],CentMax[i],ending[j]), "Primary Pions from TOF; #it{p} (GeV/#it{c}); d#it{e}d#it{x}", nPtBinsV0s, ptBinsV0s, nDeltaPiBins, deltaPiLow, deltaPiHigh);
			histPiTof[i][j]->Sumw2();

			histpPiTof[i][j]  = new TH1D(Form("hPTOF%.2f-%.2f-%s",CentMin[i],CentMax[i],ending[j]), "Primary Pions from TOF ; #it{p} (GeV/#it{c}); counts", 200, 0, 20);
			histpPiTof[i][j]->Sumw2();

			histEV0[i][j]   = new TH2D(Form("hEV0%.2f-%.2f-%s",CentMin[i],CentMax[i],ending[j]), "Electrons id by V0", nPtBinsV0s, ptBinsV0s, nDeltaPiBins, deltaPiLow, deltaPiHigh);
			histEV0[i][j]->Sumw2();

		}// eta loop
	} // centrality loop

	if(!fAnalysisMC){

		for(Int_t i=0; i<nCent; ++i ){

			fListOfObjects->Add(hMIPVsEta[i]);
			fListOfObjects->Add(pMIPVsEta[i]);
			fListOfObjects->Add(hMIPVsEtaV0s[i]);
			fListOfObjects->Add(pMIPVsEtaV0s[i]);
			fListOfObjects->Add(hPlateauVsEta[i]);
			fListOfObjects->Add(pPlateauVsEta[i]);
			fListOfObjects->Add(hPhi[i]);

			for(Int_t j=0; j<nHists; ++j){

			fListOfObjects->Add(hMIPVsNch[i][j]);
			fListOfObjects->Add(pMIPVsNch[i][j]);
			fListOfObjects->Add(hMIPVsPhi[i][j]);
			fListOfObjects->Add(pMIPVsPhi[i][j]);
			fListOfObjects->Add(hPlateauVsPhi[i][j]);
			fListOfObjects->Add(pPlateauVsPhi[i][j]);

			}

			if(fMakePid){

				fListOfObjects->Add(hPtAll[i]);
				/*fListOfObjects->Add(hPtAllNeg[i]);
				fListOfObjects->Add(hPtAllPos[i]);

				fListOfObjects->Add(hDCAxyVsPtPiNeg[i]);
				fListOfObjects->Add(hDCAxyVsPtPiPos[i]);
				//				fListOfObjects->Add(hDCAxyVsPtKNeg[i]);
				//				fListOfObjects->Add(hDCAxyVsPtKPos[i]);
				fListOfObjects->Add(hDCAxyVsPtPNeg[i]);
				fListOfObjects->Add(hDCAxyVsPtPPos[i]);

				fListOfObjects->Add(hDCAxyVsPtPiNegC[i]);
				fListOfObjects->Add(hDCAxyVsPtPiPosC[i]);
				//				fListOfObjects->Add(hDCAxyVsPtKNegC[i]);
				//				fListOfObjects->Add(hDCAxyVsPtKPosC[i]);
				fListOfObjects->Add(hDCAxyVsPtPNegC[i]);
				fListOfObjects->Add(hDCAxyVsPtPPosC[i]);
				*/

				for(Int_t j=0; j<nHists; ++j){

				/*
				fListOfObjects->Add(hnSigmaPiPos[i][j]);
				fListOfObjects->Add(hnSigmaPiNeg[i][j]);	
				fListOfObjects->Add(hnSigmaKPos[i][j]);
				fListOfObjects->Add(hnSigmaKNeg[i][j]);
				fListOfObjects->Add(hnSigmaPPos[i][j]);
				fListOfObjects->Add(hnSigmaPNeg[i][j]);
				fListOfObjects->Add(hPtPos[i][j]);
				fListOfObjects->Add(hPtNeg[i][j]);
				*/

					//					fListOfObjects->Add(histpPiV0[i][j]);
					//					fListOfObjects->Add(histpPV0[i][j]);
					//					fListOfObjects->Add(histpPiTof[i][j]);
				fListOfObjects->Add(hDeDxVsP[i][j]);
				fListOfObjects->Add(hPtVsP[i][j]);
				fListOfObjects->Add(histPiV0[i][j]);
				fListOfObjects->Add(histPiTof[i][j]);
				fListOfObjects->Add(histEV0[i][j]);
				fListOfObjects->Add(histPV0[i][j]);
				}
			}//	if(MakePID) 
		}//	Cent
	}//	!fAnalysisMC


	else{
		for(Int_t cent=0; cent<nCent; cent++) {
			for(Int_t pid=0; pid<7; pid++) {
				hMcIn[cent][pid]=new TH1D(Form("hIn%s-%s",CentName[cent],Pid[pid]), Form("MC in (pid %s)", Pid[pid]),nPtBins,ptBins);
				hMcIn[cent][pid]->Sumw2();
				hMcInNeg[cent][pid]=new TH1D(Form("hInNeg%s-%s",CentName[cent],Pid[pid]),Form("MC in (pid %s)",Pid[pid]),nPtBins,ptBins);
				hMcInNeg[cent][pid]->Sumw2();
				hMcInPos[cent][pid]=new TH1D(Form("hInPos%s-%s",CentName[cent],Pid[pid]),Form("MC in (pid %s)",Pid[pid]),nPtBins,ptBins);
				hMcInPos[cent][pid]->Sumw2();
				hMcOut[cent][pid]=new TH1D(Form("hMcOut%s-%s",CentName[cent],Pid[pid]),Form("MC out (pid %s)",Pid[pid]),nPtBins,ptBins);
				hMcOut[cent][pid]->Sumw2();
				hMcOutNeg[cent][pid]=new TH1D(Form("hMcOutNeg%s-%s",CentName[cent],Pid[pid]),Form("MC out (pid %s)",Pid[pid]),nPtBins,ptBins);
				hMcOutNeg[cent][pid]->Sumw2();
				hMcOutPos[cent][pid]=new TH1D(Form("hMcOutPos%s-%s",CentName[cent],Pid[pid]),Form("MC out (pid %s)",Pid[pid]),nPtBins,ptBins);
				hMcOutPos[cent][pid]->Sumw2();

				fListOfObjects->Add(hMcIn[cent][pid]);
				fListOfObjects->Add(hMcInNeg[cent][pid]);
				fListOfObjects->Add(hMcInPos[cent][pid]);
				fListOfObjects->Add(hMcOut[cent][pid]);
				fListOfObjects->Add(hMcOutNeg[cent][pid]);
				fListOfObjects->Add(hMcOutPos[cent][pid]);

			}	// pid Eff
		}	// cent Eff

		for(Int_t i_cent=0; i_cent<nCent; ++i_cent){
			for(Int_t pid=0; pid<7; ++pid){
				for(Int_t q=0; q<3; ++q){
					hDCApTPrim[i_cent][pid][q] = 0;
					hDCApTPrim[i_cent][pid][q] = new TH2D(Form("hDCA_%s%sPrimcent%d",Pid[pid],Q[q],i_cent),"primaries; #it{p}_{T} (GeV/#it{c}); DCA_{xy} (cm)", nPtBins, ptBins, ndcaBins, dcaBins );
					hDCApTPrim[i_cent][pid][q]->Sumw2();
					hDCApTWDec[i_cent][pid][q] = 0;
					hDCApTWDec[i_cent][pid][q] = new TH2D(Form("hDCA_%s%sWDeccent%d",Pid[pid],Q[q],i_cent),"from weak decays; #it{p}_{T} (GeV/#it{c}); DCA_{xy} (cm)", nPtBins, ptBins, ndcaBins, dcaBins );
					hDCApTWDec[i_cent][pid][q]->Sumw2();
					hDCApTMate[i_cent][pid][q] = 0;
					hDCApTMate[i_cent][pid][q] = new TH2D(Form("hDCA_%s%sMatecent%d",Pid[pid],Q[q],i_cent),"from material; #it{p}_{T} (GeV/#it{c}); DCA_{xy} (cm)", nPtBins, ptBins, ndcaBins, dcaBins );
					hDCApTMate[i_cent][pid][q]->Sumw2();

					// narrower bin size
					hDCApTPrim2[i_cent][pid][q] = 0;
					hDCApTPrim2[i_cent][pid][q] = new TH2D(Form("hDCA2_%s%sPrimcent%d",Pid[pid],Q[q],i_cent),"primaries; #it{p}_{T} (GeV/#it{c}); DCA_{xy} (cm)", nPtBins, ptBins, 800, -4.0, 4.0 );
					hDCApTPrim2[i_cent][pid][q]->Sumw2();
					hDCApTWDec2[i_cent][pid][q] = 0;
					hDCApTWDec2[i_cent][pid][q] = new TH2D(Form("hDCA2_%s%sWDeccent%d",Pid[pid],Q[q],i_cent),"from weak decays; #it{p}_{T} (GeV/#it{c}); DCA_{xy} (cm)", nPtBins, ptBins, 800, -4.0, 4.0 );
					hDCApTWDec2[i_cent][pid][q]->Sumw2();
					hDCApTMate2[i_cent][pid][q] = 0;
					hDCApTMate2[i_cent][pid][q] = new TH2D(Form("hDCA2_%s%sMatecent%d",Pid[pid],Q[q],i_cent),"from material; #it{p}_{T} (GeV/#it{c}); DCA_{xy} (cm)", nPtBins, ptBins, 800, -4.0, 4.0 );
					hDCApTMate2[i_cent][pid][q]->Sumw2();

					fListOfObjects->Add(hDCApTPrim[i_cent][pid][q]);
					fListOfObjects->Add(hDCApTPrim2[i_cent][pid][q]);
					fListOfObjects->Add(hDCApTWDec[i_cent][pid][q]);
					fListOfObjects->Add(hDCApTWDec2[i_cent][pid][q]);
					fListOfObjects->Add(hDCApTMate[i_cent][pid][q]);
					fListOfObjects->Add(hDCApTMate2[i_cent][pid][q]);
				}	// charge MC
			}	// pid MC
		}	//cent DCA MC

		fVtxMC = new TH1F("fVtxMC","mc vtx", 120, -30, 30);
		fListOfObjects->Add(fVtxMC);

	}


	fEventCuts.AddQAplotsToList(fListOfObjects);
	// Post output data.
	PostData(1, fListOfObjects);
}

//______________________________________________________________________________
void AliAnalysisTaskQAHighPtDeDxTest::UserExec(Option_t *)
{
	// Main loop

	//
	// First we make sure that we have valid input(s)!
	//

	AliVEvent *event = InputEvent();

	//	hEvents->Fill(0);
	if (!event) {
		Error("UserExec", "Could not retrieve event");
		return;
	}

	if (fAnalysisType == "ESD"){
		fESD = dynamic_cast<AliESDEvent*>(event);
		if(!fESD){
			Printf("%s:%d ESDEvent not found in Input Manager",(char*)__FILE__,__LINE__);
			this->Dump();
			return;
		}

	} else{
		fAOD = dynamic_cast<AliAODEvent*>(event);
		if(!fAOD){
			Printf("%s:%d AODEvent not found in Input Manager",(char*)__FILE__,__LINE__);
			this->Dump();
			return;
		}
	}

	if (fAnalysisMC){
		if (fAnalysisType == "ESD"){
			fMC = dynamic_cast<AliMCEvent*>(MCEvent());
			if(!fMC){
				Printf("%s:%d MCEvent not found in Input Manager",(char*)__FILE__,__LINE__);
				this->Dump();
				return;
			}

			fMCStack = fMC->Stack();

			if(!fMCStack){
				Printf("%s:%d MCStack not found in Input Manager",(char*)__FILE__,__LINE__);
				this->Dump();
				return;
			}
		} else { // AOD

			fMC = dynamic_cast<AliMCEvent*>(MCEvent());
			if(fMC)
				fMC->Dump();

			fMCArray = (TClonesArray*)fAOD->FindListObject("mcparticles");
			if(!fMCArray){
				Printf("%s:%d AOD MC array not found in Input Manager",(char*)__FILE__,__LINE__);
				this->Dump();
				return;
			}
		}
	}


	//--------------- Event Selection --------------------


	if (fAnalysisType == "ESD")
		AnalyzeESD(fESD);



	if (fAnalysisMC) {
		if (fAnalysisType == "ESD"){

			AliHeader* headerMC = fMC->Header();
			if (headerMC) {

				AliGenEventHeader* genHeader = headerMC->GenEventHeader();
				TArrayF vtxMC(3); // primary vertex  MC
				vtxMC[0]=9999; vtxMC[1]=9999;  vtxMC[2]=9999; //initialize with dummy
				if (genHeader) {
					genHeader->PrimaryVertex(vtxMC);
				}
				fZvtxMC = vtxMC[2];

				// PYTHIA:
				AliGenPythiaEventHeader* pythiaGenHeader =
					dynamic_cast<AliGenPythiaEventHeader*>(headerMC->GenEventHeader());
				if (pythiaGenHeader) {  //works only for pythia
					fMcProcessType =  GetPythiaEventProcessType(pythiaGenHeader->ProcessType());
				}
				// PHOJET:
				AliGenDPMjetEventHeader* dpmJetGenHeader =
					dynamic_cast<AliGenDPMjetEventHeader*>(headerMC->GenEventHeader());
				if (dpmJetGenHeader) {
					fMcProcessType = GetDPMjetEventProcessType(dpmJetGenHeader->ProcessType());
				}
			}
		} else { // AOD
			AliAODMCHeader* mcHeader = dynamic_cast<AliAODMCHeader*>(fAOD->FindListObject("mcHeader"));
			if(mcHeader) {
				fZvtxMC = mcHeader->GetVtxZ();
				if(strstr(mcHeader->GetGeneratorName(), "Pythia")) {
					fMcProcessType =  GetPythiaEventProcessType(mcHeader->GetEventType());
				} else {
					fMcProcessType =  GetDPMjetEventProcessType(mcHeader->GetEventType());
				}
			}
		}
	}


	PostData(1, fListOfObjects);
}

//________________________________________________________________________
void AliAnalysisTaskQAHighPtDeDxTest::AnalyzeESD(AliESDEvent* esdEvent)
{

 	fEventCuts.SetupLHC15o();
        fEventCuts.SetManualMode();

	fEventCuts.fUseStrongVarCorrelationCut = true;
        fEventCuts.fUseVariablesCorrelationCuts = true;

	if (!fEventCuts.AcceptEvent(esdEvent)){
		PostData(1, fListOfObjects);
		return;
	}



	UInt_t maskPhysSel = ((AliInputEventHandler *)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected();
	maskPhysSel &= AliVEvent::kINT7;
	if (!maskPhysSel) {
		return;
	}

	Float_t centrality = -10;

	if(fAnalysisPbPb){
		centrality = fEventCuts.GetCentrality(); /// Centrality calculated with the default estimator (V0M for LHC15o)

		if((centrality>fMaxCent)||(centrality<fMinCent))return;
		for(Int_t icent = 0; icent<(nCent-1); ++icent){
			if(centrality > CentMin[icent] && centrality <= CentMax[icent]){
				cent = icent;
				fcent->Fill(icent+1);
				ProduceArrayTrksESD( esdEvent, cent );
				ProduceArrayV0ESD( esdEvent, cent );

				if(fAnalysisMC)
				ProcessMCTruthESD(cent);
			}
		}

		fcent->Fill(11);

	}



}

//________________________________________________________________________
void AliAnalysisTaskQAHighPtDeDxTest::AnalyzeAOD(AliAODEvent* aodEvent)
{
	fRun  = aodEvent->GetRunNumber();
	fEventId = 0;
	if(aodEvent->GetHeader())
		fEventId = GetEventIdAsLong(aodEvent->GetHeader());

	//UInt_t    time      = 0; // Missing AOD info? aodEvent->GetTimeStamp();
	magf      = aodEvent->GetMagneticField();

	//Int_t     trackmult = 0; // no pt cuts
	//Int_t     nadded    = 0;

	Bool_t isPileup = aodEvent->IsPileupFromSPD();
	if(fPileUpRej)
		if(isPileup)
			return;
	//fn1->Fill(4);



	if(fTriggeredEventMB) {// Only MC case can we have not triggered events

		// accepted event
		//fEvents->Fill(0);

		//if(fVtxStatus!=1) return; // accepted vertex
		//Int_t nAODTracks = aodEvent->GetNumberOfTracks();

		ProduceArrayTrksAOD( aodEvent );
		ProduceArrayV0AOD( aodEvent );

		//fEvents->Fill(1);




	} // end if triggered

}

//_____________________________________________________________________________
Float_t AliAnalysisTaskQAHighPtDeDxTest::GetVertex(const AliVEvent* event) const
{
	Float_t zvtx = -999;

	const AliVVertex* primaryVertex = event->GetPrimaryVertex();

	if(primaryVertex->GetNContributors()>0)
		zvtx = primaryVertex->GetZ();

	return zvtx;
}

//_____________________________________________________________________________
Short_t AliAnalysisTaskQAHighPtDeDxTest::GetPidCode(Int_t pdgCode) const
{
	// return our internal code for pions, kaons, and protons

	Short_t pidCode = 6;

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
		case 11:
			pidCode = 4; // electron
			break;
		case 13:
			pidCode = 5; // muon
			break;
		default:
			pidCode = 6;  // something else?
	};

	return pidCode;
}

//_____________________________________________________________________________
void AliAnalysisTaskQAHighPtDeDxTest::ProcessMCTruthESD(const Int_t Cent)
{
	// Fill the special MC histogram with the MC truth info

	cout<<"Cent Inside ProcessMCTruth ::: "<<Cent<<endl;
	const Int_t nTracksMC = fMCStack->GetNtrack();

	for (Int_t iTracks = 0; iTracks < nTracksMC; iTracks++) {

		TParticle* trackMC = fMCStack->Particle(iTracks);
		if( !trackMC )
			continue;

		if( !(fMCStack->IsPhysicalPrimary(iTracks)) )
			continue;

		TParticlePDG* pdgPart = trackMC ->GetPDG();
		Double_t chargeMC = pdgPart->Charge();

		if(chargeMC==0)
			continue;

		if ( TMath::Abs(trackMC->Eta()) > fEtaCut )
			continue;

		//		if ( TMath::Abs(trackMC->Y()) > 0.5 )
		//			continue;

		Int_t pdgCode = trackMC->GetPdgCode();
		Short_t pidCodeMC = 0;
		pidCodeMC = GetPidCode(pdgCode);

		hMcIn[Cent][0]->Fill(trackMC->Pt());
		hMcIn[Cent][pidCodeMC]->Fill(trackMC->Pt());

		hMcIn[10][0]->Fill(trackMC->Pt());
		hMcIn[10][pidCodeMC]->Fill(trackMC->Pt());

		if( chargeMC < 0 ){
			hMcInNeg[Cent][0]->Fill(trackMC->Pt());
			hMcInNeg[Cent][pidCodeMC]->Fill(trackMC->Pt());

			hMcInNeg[10][0]->Fill(trackMC->Pt());
			hMcInNeg[10][pidCodeMC]->Fill(trackMC->Pt());
		}
		else{
			hMcInPos[Cent][0]->Fill(trackMC->Pt());
			hMcInPos[Cent][pidCodeMC]->Fill(trackMC->Pt());

			hMcInPos[10][0]->Fill(trackMC->Pt());
			hMcInPos[10][pidCodeMC]->Fill(trackMC->Pt());
		}


	}//MC track loop



}

//_____________________________________________________________________________
void AliAnalysisTaskQAHighPtDeDxTest::ProcessMCTruthAOD()
{
	// Fill the special MC histogram with the MC truth info


	const Int_t nTracksMC = fMCArray->GetEntriesFast();

	for (Int_t iTracks = 0; iTracks < nTracksMC; iTracks++) {

		AliAODMCParticle* trackMC = dynamic_cast<AliAODMCParticle*>(fMCArray->At(iTracks));

		if(!trackMC){
			AliError("Cannot get MC particle");
			continue;
		}

		//Cuts
		if(!(trackMC->IsPhysicalPrimary()))
			continue;


		Double_t chargeMC = trackMC->Charge();
		if(chargeMC==0)
			continue;


		if (TMath::Abs(trackMC->Eta()) > fEtaCut )
			continue;

		//		Double_t etaMC = trackMC->Eta();
		Int_t pdgCode = trackMC->GetPdgCode();
		Short_t pidCodeMC = 0;
		pidCodeMC = GetPidCode(pdgCode);

		////cout<<"pidcode="<<pidCodeMC<<endl;
		/*		for(Int_t nh = 0; nh < 9; nh++) {

				if( etaMC > etaHigh[nh]/10.0 || etaMC < etaLow[nh]/10.0 )
				continue;

				hMcIn[0][nh]->Fill(trackMC->Pt());
				hMcIn[pidCodeMC][nh]->Fill(trackMC->Pt());


				}
		 */

	}//MC track loop


}


//_____________________________________________________________________________
Short_t AliAnalysisTaskQAHighPtDeDxTest::GetPythiaEventProcessType(Int_t pythiaType) {
	//
	// Get the process type of the event.  PYTHIA
	//
	// source PWG0   dNdpt

	Short_t globalType = -1; //init

	if(pythiaType==92||pythiaType==93){
		globalType = 2; //single diffractive
	}
	else if(pythiaType==94){
		globalType = 3; //double diffractive
	}
	//else if(pythiaType != 91){ // also exclude elastic to be sure... CKB??
	else {
		globalType = 1;  //non diffractive
	}
	return globalType;
}

//_____________________________________________________________________________
Short_t AliAnalysisTaskQAHighPtDeDxTest::GetDPMjetEventProcessType(Int_t dpmJetType) {
	//
	// get the process type of the event.  PHOJET
	//
	//source PWG0   dNdpt
	// can only read pythia headers, either directly or from cocktalil header
	Short_t globalType = -1;

	if (dpmJetType == 1 || dpmJetType == 4) { // explicitly inelastic plus central diffraction
		globalType = 1;
	}
	else if (dpmJetType==5 || dpmJetType==6) {
		globalType = 2;
	}
	else if (dpmJetType==7) {
		globalType = 3;
	}
	return globalType;
}

//_____________________________________________________________________________
ULong64_t AliAnalysisTaskQAHighPtDeDxTest::GetEventIdAsLong(AliVHeader* header) const
{
	// To have a unique id for each event in a run!
	// Modified from AliRawReader.h
	return ((ULong64_t)header->GetBunchCrossNumber()+
			(ULong64_t)header->GetOrbitNumber()*3564+
			(ULong64_t)header->GetPeriodNumber()*16777215*3564);
}


//____________________________________________________________________
TParticle* AliAnalysisTaskQAHighPtDeDxTest::FindPrimaryMother(AliStack* stack, Int_t label)
{
	//
	// Finds the first mother among the primary particles of the particle identified by <label>,
	// i.e. the primary that "caused" this particle
	//
	// Taken from AliPWG0Helper class
	//

	Int_t motherLabel = FindPrimaryMotherLabel(stack, label);
	if (motherLabel < 0)
		return 0;

	return stack->Particle(motherLabel);
}

//____________________________________________________________________
Int_t AliAnalysisTaskQAHighPtDeDxTest::FindPrimaryMotherLabel(AliStack* stack, Int_t label)
{
	//
	// Finds the first mother among the primary particles of the particle identified by <label>,
	// i.e. the primary that "caused" this particle
	//
	// returns its label
	//
	// Taken from AliPWG0Helper class
	//
	const Int_t nPrim  = stack->GetNprimary();

	while (label >= nPrim) {

		//printf("Particle %d (pdg %d) is not a primary. Let's check its mother %d\n", label, mother->GetPdgCode(), mother->GetMother(0));

		TParticle* particle = stack->Particle(label);
		if (!particle) {

			AliDebugGeneral("FindPrimaryMotherLabel", AliLog::kError, Form("UNEXPECTED: particle with label %d not found in stack.", label));
			return -1;
		}

		// find mother
		if (particle->GetMother(0) < 0) {

			AliDebugGeneral("FindPrimaryMotherLabel", AliLog::kError, Form("UNEXPECTED: Could not find mother of secondary particle %d.", label));
			return -1;
		}

		label = particle->GetMother(0);
	}

	return label;
}

//____________________________________________________________________
AliAODMCParticle* AliAnalysisTaskQAHighPtDeDxTest::FindPrimaryMotherAOD(AliAODMCParticle* startParticle)
{
	//
	// Finds the first mother among the primary particles of the particle identified by <label>,
	// i.e. the primary that "caused" this particle
	//
	// Taken from AliPWG0Helper class
	//

	AliAODMCParticle* mcPart = startParticle;

	while (mcPart)
	{

		if(mcPart->IsPrimary())
			return mcPart;

		Int_t mother = mcPart->GetMother();

		mcPart = dynamic_cast<AliAODMCParticle*>(fMCArray->At(mother));
	}

	return 0;
}


//V0______________________________________
//____________________________________________________________________
TParticle* AliAnalysisTaskQAHighPtDeDxTest::FindPrimaryMotherV0(AliStack* stack, Int_t label)
{
	//
	// Finds the first mother among the primary particles of the particle identified by <label>,
	// i.e. the primary that "caused" this particle
	//
	// Taken from AliPWG0Helper class
	//

	Int_t nSteps = 0;

	Int_t motherLabel = FindPrimaryMotherLabelV0(stack, label, nSteps);
	if (motherLabel < 0)
		return 0;

	return stack->Particle(motherLabel);
}

//____________________________________________________________________
Int_t AliAnalysisTaskQAHighPtDeDxTest::FindPrimaryMotherLabelV0(AliStack* stack, Int_t label, Int_t& nSteps)
{
	//
	// Finds the first mother among the primary particles of the particle identified by <label>,
	// i.e. the primary that "caused" this particle
	//
	// returns its label
	//
	// Taken from AliPWG0Helper class
	//
	nSteps = 0;
	const Int_t nPrim  = stack->GetNprimary();

	while (label >= nPrim) {

		//printf("Particle %d (pdg %d) is not a primary. Let's check its mother %d\n", label, mother->GetPdgCode(), mother->GetMother(0));

		nSteps++; // 1 level down

		TParticle* particle = stack->Particle(label);
		if (!particle) {

			AliDebugGeneral("FindPrimaryMotherLabelV0", AliLog::kError, Form("UNEXPECTED: particle with label %d not found in stack.", label));
			return -1;
		}

		// find mother
		if (particle->GetMother(0) < 0) {

			AliDebugGeneral("FindPrimaryMotherLabelV0", AliLog::kError, Form("UNEXPECTED: Could not find mother of secondary particle %d.", label));
			return -1;
		}

		label = particle->GetMother(0);
	}

	return label;
}

//____________________________________________________________________
AliAODMCParticle* AliAnalysisTaskQAHighPtDeDxTest::FindPrimaryMotherAODV0(AliAODMCParticle* startParticle, Int_t& nSteps)
{
	//
	// Finds the first mother among the primary particles of the particle identified by <label>,
	// i.e. the primary that "caused" this particle
	//
	// Taken from AliPWG0Helper class
	//

	nSteps = 0;

	AliAODMCParticle* mcPart = startParticle;

	while (mcPart)
	{

		if(mcPart->IsPrimary())
			return mcPart;

		Int_t mother = mcPart->GetMother();

		mcPart = dynamic_cast<AliAODMCParticle*>(fMCArray->At(mother));
		nSteps++; // 1 level down
	}

	return 0;
}



//__________________________________________________________________
void AliAnalysisTaskQAHighPtDeDxTest::ProduceArrayTrksESD( AliESDEvent *ESDevent, const Int_t Cent ){

	const Int_t nESDTracks = ESDevent->GetNumberOfTracks();

	fcentAfterPrimaries->Fill(Cent+1);
	fcentAfterPrimaries->Fill(11);

	Int_t multTPC = 0;
	for(Int_t iT = 0; iT < nESDTracks; iT++) {

		AliESDtrack* esdTrack = ESDevent->GetTrack(iT);

		if( TMath::Abs(esdTrack->Eta()) > fEtaCut )
			continue;

		UInt_t selectDebug = 0;
		if (fTrackFilterTPC) {
			selectDebug = fTrackFilterTPC->IsSelected(esdTrack);
			if (!selectDebug) {
				continue;
			}
		}

		multTPC++;

	}

	for(Int_t iT = 0; iT < nESDTracks; iT++) {

		AliESDtrack* esdTrack = ESDevent->GetTrack(iT);
		if(TMath::Abs(esdTrack->Eta()) > fEtaCut)
			continue;

		// TracCuts 2015 Pb-Pb Open DCA
		UInt_t selectDebug = 0;
		if (fTrackFilter2015PbPb) {
			selectDebug = fTrackFilter2015PbPb->IsSelected(esdTrack);
			if (!selectDebug) {
				continue;
			}
		}

		Short_t ncl = esdTrack->GetTPCsignalN();
		if(ncl<fNcl)
		continue;
		
		if(ncl < 70)cout<<" ¡¡¡¡¡¡¡¡¡¡¡     ncl > 70"<<endl;
		

		Double_t eta      = esdTrack->Eta();
		Double_t phi      = esdTrack->Phi();
		Double_t momentum = esdTrack->P();
		Double_t pt       = esdTrack->Pt();
		Double_t dedx     = esdTrack->GetTPCsignal();
		Double_t dedxmb   = esdTrack->GetTPCsignal();
		Double_t dedxUnc  = esdTrack->GetTPCsignal();

		Float_t dcaxy = 0.;
		Float_t dcaz = 0.;
		esdTrack->GetImpactParameters(dcaxy,dcaz);

		if(!PhiCut(esdTrack->Pt(), phi, esdTrack->Charge(), magf, &cutLow, &cutHigh))
			continue;

		if(fdEdxCalibrated){
			if(eta < 0){
				dedx   *= 50/EtaCalibrationNeg(Cent,eta);
				dedxmb *= 50/EtaCalibrationNeg(10,eta);
			}
			else{
				dedx   *= 50/EtaCalibrationPos(Cent,eta);
				dedxmb *= 50/EtaCalibrationPos(10,eta);
			}
		}

		/*if(esdTrack->Charge() < 0.0){
			if(TMath::Sqrt(TMath::Power(fPIDResponse->NumberOfSigmasTPC(esdTrack,AliPID::kPion),2)+TMath::Power(fPIDResponse->NumberOfSigmasTOF(esdTrack,AliPID::kPion),2))<2.0){
				hDCAxyVsPtPiNeg[Cent]->Fill(pt,dcaxy);
				hDCAxyVsPtPiNegC[Cent]->Fill(pt,dcaxy);

			}
			if(TMath::Sqrt(TMath::Power(fPIDResponse->NumberOfSigmasTPC(esdTrack,AliPID::kKaon),2)+TMath::Power(fPIDResponse->NumberOfSigmasTOF(esdTrack,AliPID::kKaon),2))<2.0){
				hDCAxyVsPtKNeg[Cent]->Fill(pt,dcaxy);
				hDCAxyVsPtKNegC[Cent]->Fill(pt,dcaxy);

			}
			if(TMath::Sqrt(TMath::Power(fPIDResponse->NumberOfSigmasTPC(esdTrack,AliPID::kProton),2)+TMath::Power(fPIDResponse->NumberOfSigmasTOF(esdTrack,AliPID::kProton),2))<2.0){
				hDCAxyVsPtPNeg[Cent]->Fill(pt,dcaxy);
				hDCAxyVsPtPNegC[Cent]->Fill(pt,dcaxy);
			}
		}*/
		/*else{
			if(TMath::Sqrt(TMath::Power(fPIDResponse->NumberOfSigmasTPC(esdTrack,AliPID::kPion),2)+TMath::Power(fPIDResponse->NumberOfSigmasTOF(esdTrack,AliPID::kPion),2))<2.0){
				hDCAxyVsPtPiPos[Cent]->Fill(pt,dcaxy);
				hDCAxyVsPtPiPosC[Cent]->Fill(pt,dcaxy);
			}
			if(TMath::Sqrt(TMath::Power(fPIDResponse->NumberOfSigmasTPC(esdTrack,AliPID::kKaon),2)+TMath::Power(fPIDResponse->NumberOfSigmasTOF(esdTrack,AliPID::kKaon),2))<2.0){
				hDCAxyVsPtKPos[Cent]->Fill(pt,dcaxy);
				hDCAxyVsPtKPosC[Cent]->Fill(pt,dcaxy);
			}
			if(TMath::Sqrt(TMath::Power(fPIDResponse->NumberOfSigmasTPC(esdTrack,AliPID::kProton),2)+TMath::Power(fPIDResponse->NumberOfSigmasTOF(esdTrack,AliPID::kProton),2))<2.0){
				hDCAxyVsPtPPos[Cent]->Fill(pt,dcaxy);
				hDCAxyVsPtPPosC[Cent]->Fill(pt,dcaxy);
			}
		}*/


		Short_t pidCode     = 0;
		if(fAnalysisMC) {

			const Int_t label = TMath::Abs(esdTrack->GetLabel());
			TParticle* mcTrack = 0;
			mcTrack = fMCStack->Particle(label);

			if (mcTrack){

				if( esdTrack->Charge()==0 )
					continue;

				Int_t pdgCode = mcTrack->GetPdgCode();
				pidCode = GetPidCode(pdgCode);

				if( fMCStack->IsPhysicalPrimary(label) ){

					if( TMath::Abs(dcaxy) < GetMaxDCApTDep(fcutDCAxy,pt) ){
						hMcOut[Cent][0]->Fill(esdTrack->Pt());
						hMcOut[Cent][pidCode]->Fill(esdTrack->Pt());

						hMcOut[10][0]->Fill(esdTrack->Pt());
						hMcOut[10][pidCode]->Fill(esdTrack->Pt());
					}

					hDCApTPrim[Cent][0][0]->Fill(pt,dcaxy);
					hDCApTPrim[Cent][pidCode][0]->Fill(pt,dcaxy);

					hDCApTPrim2[Cent][0][0]->Fill(pt,dcaxy);
					hDCApTPrim2[Cent][pidCode][0]->Fill(pt,dcaxy);

					if( esdTrack->Charge() < 0.0 ){

						if( TMath::Abs(dcaxy) < GetMaxDCApTDep(fcutDCAxy,pt) )  {
							hMcOutNeg[Cent][0]->Fill(esdTrack->Pt());
							hMcOutNeg[Cent][pidCode]->Fill(esdTrack->Pt());

							hMcOutNeg[10][0]->Fill(esdTrack->Pt());
							hMcOutNeg[10][pidCode]->Fill(esdTrack->Pt());
						}

						hDCApTPrim[Cent][0][1]->Fill(pt,dcaxy);
						hDCApTPrim[Cent][pidCode][1]->Fill(pt,dcaxy);

						hDCApTPrim2[Cent][0][1]->Fill(pt,dcaxy);
						hDCApTPrim2[Cent][pidCode][1]->Fill(pt,dcaxy);
					}
					else{

						if( TMath::Abs(dcaxy) < GetMaxDCApTDep(fcutDCAxy,pt) ){
							hMcOutPos[Cent][0]->Fill(esdTrack->Pt());
							hMcOutPos[Cent][pidCode]->Fill(esdTrack->Pt());

							hMcOutPos[10][0]  ->Fill(esdTrack->Pt());
							hMcOutPos[10][pidCode]  ->Fill(esdTrack->Pt());
						}

						hDCApTPrim[Cent][0][2]->Fill(pt,dcaxy);
						hDCApTPrim[Cent][pidCode][2]->Fill(pt,dcaxy);

						hDCApTPrim2[Cent][0][2]->Fill(pt,dcaxy);
						hDCApTPrim2[Cent][pidCode][2]->Fill(pt,dcaxy);
					}
				}	// Primary particles MC

				if( fMCStack->IsSecondaryFromWeakDecay(label) ){

					hDCApTWDec[Cent][0][0]->Fill(pt,dcaxy);
					hDCApTWDec[Cent][pidCode][0]->Fill(pt,dcaxy);

					hDCApTWDec2[Cent][0][0]->Fill(pt,dcaxy);
					hDCApTWDec2[Cent][pidCode][0]->Fill(pt,dcaxy);

					if( esdTrack->Charge() < 0.0 ){

						hDCApTWDec[Cent][0][1]->Fill(pt,dcaxy);
						hDCApTWDec[Cent][pidCode][1]->Fill(pt,dcaxy);

						hDCApTWDec2[Cent][0][1]->Fill(pt,dcaxy);
						hDCApTWDec2[Cent][pidCode][1]->Fill(pt,dcaxy);

					}
					else{

						hDCApTWDec[Cent][0][2]->Fill(pt,dcaxy);
						hDCApTWDec[Cent][pidCode][2]->Fill(pt,dcaxy);

						hDCApTWDec2[Cent][0][2]->Fill(pt,dcaxy);
						hDCApTWDec2[Cent][pidCode][2]->Fill(pt,dcaxy);

					}
				}	// Weak Decay MC

				if( fMCStack->IsSecondaryFromMaterial(label) ){

					hDCApTMate[Cent][0][0]->Fill(pt,dcaxy);
					hDCApTMate[Cent][pidCode][0]->Fill(pt,dcaxy);

					hDCApTMate2[Cent][0][0]->Fill(pt,dcaxy);
					hDCApTMate2[Cent][pidCode][0]->Fill(pt,dcaxy);

					if( esdTrack->Charge() < 0.0 ){

						hDCApTMate[Cent][0][1]->Fill(pt,dcaxy);
						hDCApTMate[Cent][pidCode][1]->Fill(pt,dcaxy);

						hDCApTMate2[Cent][0][1]->Fill(pt,dcaxy);
						hDCApTMate2[Cent][pidCode][1]->Fill(pt,dcaxy);

					}
					else{

						hDCApTMate[Cent][0][2]->Fill(pt,dcaxy);
						hDCApTMate[Cent][pidCode][2]->Fill(pt,dcaxy);

						hDCApTMate2[Cent][0][2]->Fill(pt,dcaxy);
						hDCApTMate2[Cent][pidCode][2]->Fill(pt,dcaxy);

					}
				}	// Material Inte MC
			}	//mcTrack
		}	//fAnalysis MC

		// DCAxy cut
		if( TMath::Abs(dcaxy) > GetMaxDCApTDep(fcutDCAxy,pt) )
			continue;

		Bool_t IsTOFout=kFALSE;
		if ((esdTrack->GetStatus()&AliESDtrack::kTOFout)==0)
			IsTOFout=kTRUE;
		Float_t lengthtrack=esdTrack->GetIntegratedLength();
		Float_t timeTOF=esdTrack->GetTOFsignal();
		Double_t inttime[5]={0,0,0,0,0};
		esdTrack->GetIntegratedTimes(inttime);// Returns the array with integrated times for each particle hypothesis
		Float_t beta = -99;
		if ( !IsTOFout ){
			if ( ( lengthtrack != 0 ) && ( timeTOF != 0) )
				beta = inttime[0] / timeTOF;
		}

		if(!fdEdxCalibrated){
			if( momentum <= 0.6 && momentum >= 0.4 ){//only p:0.4-0.6 GeV, pion MIP
				if( dedxUnc < DeDxMIPMax && dedxUnc > DeDxMIPMin ){
					hMIPVsEta[Cent]->Fill(eta,dedxUnc);
					pMIPVsEta[Cent]->Fill(eta,dedxUnc);
					hMIPVsEta[10]->Fill(eta,dedxUnc);
					pMIPVsEta[10]->Fill(eta,dedxUnc);
				}
				if( dedxUnc > 70 && dedxUnc < 90 ){
					if(TMath::Abs(beta-1)<0.1){
						hPlateauVsEta[Cent]->Fill(eta,dedxUnc);
						pPlateauVsEta[Cent]->Fill(eta,dedxUnc);
						hPlateauVsEta[10]->Fill(eta,dedxUnc);
						pPlateauVsEta[10]->Fill(eta,dedxUnc);
					}
				}
			}
		}
		else{
			if( momentum <= 0.6 && momentum >= 0.4 ){//only p:0.4-0.6 GeV, pion MIP
				if( dedxUnc < DeDxMIPMax && dedxUnc > DeDxMIPMin ){
					hMIPVsEta[Cent]->Fill(eta,dedx);
					pMIPVsEta[Cent]->Fill(eta,dedx);
					hMIPVsEta[10]->Fill(eta,dedxmb);
					pMIPVsEta[10]->Fill(eta,dedxmb);
				}
				if( dedxUnc > 70 && dedxUnc < 90 ){
					if(TMath::Abs(beta-1)<0.1){
						hPlateauVsEta[Cent]->Fill(eta,dedx);
						pPlateauVsEta[Cent]->Fill(eta,dedx);
						hPlateauVsEta[10]->Fill(eta,dedxmb);
						pPlateauVsEta[10]->Fill(eta,dedxmb);
					}
				}
			}
		}

		hPtAll[Cent]->Fill(pt);
		hPtAll[10]->Fill(pt);

		if(esdTrack->Charge() < 0.){
			hPtAllNeg[Cent]->Fill(pt);
			//			if(TMath::Abs(fPIDResponse->NumberOfSigmasTPC(esdTrack,AliPID::kPion))<3.0)
			//				hDCAxyVsPtPiNegC[Cent]->Fill(pt,dcaxy);
			//			if(TMath::Abs(fPIDResponse->NumberOfSigmasTPC(esdTrack,AliPID::kKaon))<3.0)
			//				hDCAxyVsPtKNegC[Cent]->Fill(pt,dcaxy);
			//			if(TMath::Abs(fPIDResponse->NumberOfSigmasTPC(esdTrack,AliPID::kProton))<3.0)
			//				hDCAxyVsPtPNegC[Cent]->Fill(pt,dcaxy);
		}
		else{
			hPtAllPos[Cent]->Fill(pt);
			//			if(TMath::Abs(fPIDResponse->NumberOfSigmasTPC(esdTrack,AliPID::kPion))<3.0)
			//				hDCAxyVsPtPiPosC[Cent]->Fill(pt,dcaxy);
			//			if(TMath::Abs(fPIDResponse->NumberOfSigmasTPC(esdTrack,AliPID::kKaon))<3.0)
			//				hDCAxyVsPtKPosC[Cent]->Fill(pt,dcaxy);
			//			if(TMath::Abs(fPIDResponse->NumberOfSigmasTPC(esdTrack,AliPID::kProton))<3.0)
			//				hDCAxyVsPtPPosC[Cent]->Fill(pt,dcaxy);
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

		if(nh<0)
			continue;

		/*
		   if(fAnalysisMC){

		   if( !(fMCStack->IsPhysicalPrimary(iT)) )
		   continue;

		   if( esdTrack->Charge()==0 )
		   continue;

		   if ( TMath::Abs(esdTrack->Y()) > 0.5 )
		   continue;

		   hMcOut[Cent][0]->Fill(esdTrack->Pt());
		   hMcOut[Cent][pidCode]->Fill(esdTrack->Pt());

		   hMcOut[10][0]->Fill(esdTrack->Pt());
		   hMcOut[10][pidCode]->Fill(esdTrack->Pt());

		   if( esdTrack->Charge() < 0.0 ){
		   hMcOutNeg[Cent][0]->Fill(esdTrack->Pt());
		   hMcOutNeg[Cent][pidCode]->Fill(esdTrack->Pt());

		   hMcOutNeg[10][0]->Fill(esdTrack->Pt());
		   hMcOutNeg[10][pidCode]->Fill(esdTrack->Pt());
		   }
		   else{
		   hMcOutPos[Cent][0]->Fill(esdTrack->Pt());
		   hMcOutPos[Cent][pidCode]->Fill(esdTrack->Pt());

		   hMcOutPos[10][0]  ->Fill(esdTrack->Pt());
		   hMcOutPos[10][pidCode]  ->Fill(esdTrack->Pt());
		   }
		   }
		 */

		if(beta>1){
			histPiTof[Cent][nh]->Fill(momentum, dedx);
			histPiTof[10][nh]->Fill(momentum, dedxmb);
		}

		if( momentum <= 0.6 && momentum >= 0.4  ){

			if( dedx < DeDxMIPMax && dedx > DeDxMIPMin ){
				hMIPVsPhi[Cent][nh]->Fill(phi,dedx);
				pMIPVsPhi[Cent][nh]->Fill(phi,dedx);

				hMIPVsPhi[10][nh]->Fill(phi,dedxmb);
				pMIPVsPhi[10][nh]->Fill(phi,dedxmb);

				hMIPVsNch[Cent][nh]->Fill(multTPC,dedx);
				pMIPVsNch[Cent][nh]->Fill(multTPC,dedx);

				hMIPVsNch[10][nh]->Fill(multTPC,dedxmb);
				pMIPVsNch[10][nh]->Fill(multTPC,dedxmb);

			}
			if( dedx > 70 && dedx < 90 ){
				if(TMath::Abs(beta-1)<0.1){
					hPlateauVsPhi[Cent][nh]->Fill(phi,dedx);
					pPlateauVsPhi[Cent][nh]->Fill(phi,dedx);

					hPlateauVsPhi[10][nh]->Fill(phi,dedxmb);
					pPlateauVsPhi[10][nh]->Fill(phi,dedxmb);
				}
			}
		}

		hPtVsP[Cent][nh]->Fill(momentum,pt);
		hPtVsP[10][nh]->Fill(momentum,pt);
		hDeDxVsP[Cent][nh]->Fill(momentum,dedx);
		hDeDxVsP[10][nh]->Fill(momentum,dedxmb);

		if(esdTrack->Charge() < 0.){
			hPtNeg[Cent][nh]->Fill(pt);

			hnSigmaPiNeg[Cent][nh]->Fill(pt,fPIDResponse->NumberOfSigmasTPC(esdTrack,AliPID::kPion));
			hnSigmaKNeg[Cent][nh]->Fill(pt,fPIDResponse->NumberOfSigmasTPC(esdTrack,AliPID::kKaon));
			hnSigmaPNeg[Cent][nh]->Fill(pt,fPIDResponse->NumberOfSigmasTPC(esdTrack,AliPID::kProton));
		}
		else{
			hPtPos[Cent][nh]->Fill(pt);

			hnSigmaPiPos[Cent][nh]->Fill(pt,fPIDResponse->NumberOfSigmasTPC(esdTrack,AliPID::kPion));
			hnSigmaKPos[Cent][nh]->Fill(pt,fPIDResponse->NumberOfSigmasTPC(esdTrack,AliPID::kKaon));
			hnSigmaPPos[Cent][nh]->Fill(pt,fPIDResponse->NumberOfSigmasTPC(esdTrack,AliPID::kProton));
		}



	}//end of track loop




}
//__________________________________________________________________
void AliAnalysisTaskQAHighPtDeDxTest::ProduceArrayTrksAOD( AliAODEvent *AODevent ){

	//Bool_t OfficialCalibration = kTRUE;
	Int_t nAODTracks = AODevent->GetNumberOfTracks();
	Int_t multTPC = 0;

	//get multiplicity tpc only track cuts
	for(Int_t iT = 0; iT < nAODTracks; iT++) {

		//AliAODTrack* aodTrack = AODevent->GetTrack(iT);
		AliAODTrack* aodTrack = static_cast<AliAODTrack*>(AODevent->GetTrack(iT));
		if(!aodTrack)printf("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");


		if(TMath::Abs(aodTrack->Eta()) > fEtaCut)
			continue;


		if (fTrackFilterTPC) {
			// TPC only cuts is bit 1
			if(!aodTrack->TestFilterBit(1))
				continue;
		}

		multTPC++;

	}


	for(Int_t iT = 0; iT < nAODTracks; iT++) {

		AliAODTrack* aodTrack = static_cast<AliAODTrack*>(AODevent->GetTrack(iT));
		AliVTrack* const trackPID = dynamic_cast<AliVTrack*>(aodTrack);
		if(!trackPID)cout << "NO trackPID    " << endl;

		if (fTrackFilterGolden) {
			// "Global track RAA analysis QM2011 + Chi2ITS<36"; bit 1024
			if(!aodTrack->TestFilterBit(1024))
				continue;
		}


		if(TMath::Abs(aodTrack->Eta()) > fEtaCut)
			continue;


		Double_t eta  = aodTrack->Eta();
		Double_t phi  = aodTrack->Phi();
		Double_t momentum = aodTrack->P();
		Double_t pt = aodTrack->Pt();


		if(!PhiCut(aodTrack->Pt(), phi, aodTrack->Charge(), magf, &cutLow, &cutHigh))
			continue;


		AliAODPid* aodPid = aodTrack->GetDetPid();
		Short_t ncl     = -10;
		Float_t dedx    = -10;

		//TOF
		Float_t beta = -99;


		if(aodPid) {
			ncl     = aodPid->GetTPCsignalN();
			//			if(!OfficialCalibration)
			dedx    = aodPid->GetTPCsignal();
			//			else
			// 			dedx = fPIDResponse->GetTPCResponse().GetEtaCorrectedTrackdEdx( aodTrack, AliPID::kPion,      AliTPCPIDResponse::kdEdxDefault );

			//TOF
			Bool_t IsTOFout=kFALSE;
			Float_t lengthtrack=-999;//in aod we do not have that information, beta must be: beta=inttime/timeTOF
			Float_t timeTOF=-999;

			if ((aodTrack->GetStatus()&AliESDtrack::kTOFout)==0)
				IsTOFout=kTRUE;

			lengthtrack=-999;//in aod we do not have that information, beta must be: beta=inttime/timeTOF

			timeTOF=aodPid->GetTOFsignal();

			Double_t inttime[5]={0,0,0,0,0};
			aodPid->GetIntegratedTimes(inttime);// Returns the array with integrated times for each particle hypothesis


			if ( !IsTOFout ){
				if ( ( lengthtrack != 0 ) && ( timeTOF != 0) )
					beta = inttime[0] / timeTOF;
			}

		}


		if(ncl<fNcl)
			continue;


		Short_t pidCode     = 0;

		if(fAnalysisMC) {

			const Int_t label = TMath::Abs(aodTrack->GetLabel());
			AliAODMCParticle* mcTrack = dynamic_cast<AliAODMCParticle*>(fMCArray->At(label));

			if (mcTrack){

				Int_t pdgCode = mcTrack->GetPdgCode();
				pidCode = GetPidCode(pdgCode);

			}

		}

		if( momentum <= 0.6 && momentum >= 0.4  ){//only p:0.4-0.6 GeV, pion MIP
			if( dedx < DeDxMIPMax && dedx > DeDxMIPMin ){
				if(momentum<0.6&&momentum>0.4){
					hMIPVsEta[cent]->Fill(eta,dedx);
					pMIPVsEta[cent]->Fill(eta,dedx);
				}
			}
			if( dedx > DeDxMIPMax+1 && dedx < 95 ){
				if(TMath::Abs(beta-1)<0.1){
					hPlateauVsEta[cent]->Fill(eta,dedx);
					pPlateauVsEta[cent]->Fill(eta,dedx);
				}
			}
		}


		for(Int_t nh = 0; nh < 9; nh++){

			if( eta > etaHigh[nh]/10.0 || eta < etaLow[nh]/10.0 )
				continue;

			/*			if(fAnalysisMC){
						hMcOut[0][nh]->Fill(aodTrack->Pt());
						hMcOut[pidCode][nh]->Fill(aodTrack->Pt());
						}
			 */
			if(beta>1){
				histPiTof[cent][nh]->Fill(momentum,dedx);
				histpPiTof[cent][nh]->Fill(momentum);
			}

			if( momentum <= 0.6 && momentum >= 0.4  ){//only p:0.4-0.6 GeV, pion MIP
				//Fill  pion MIP, before calibration
				if( dedx < DeDxMIPMax && dedx > DeDxMIPMin ){
					hMIPVsPhi[cent][nh]->Fill(phi,dedx);
					pMIPVsPhi[cent][nh]->Fill(phi,dedx);

					hMIPVsNch[cent][nh]->Fill(multTPC,dedx);
					pMIPVsNch[cent][nh]->Fill(multTPC,dedx);

				}

				//Fill electrons, before calibration
				if( dedx > DeDxMIPMax+1 && dedx < 95 ){
					if(TMath::Abs(beta-1)<0.1){
						hPlateauVsPhi[cent][nh]->Fill(phi,dedx);
						pPlateauVsPhi[cent][nh]->Fill(phi,dedx);
					}
				}
			}

			//				hPt[nh]->Fill(pt);
			hPtVsP[cent][nh]->Fill(momentum,pt);
			//				histAllCh[nh]->Fill(momentum,dedx);
			hDeDxVsP[cent][nh]->Fill(momentum,dedx);


			if(aodTrack->Charge() < 0.){
				hPtNeg[cent][nh]->Fill(pt);
				hnSigmaPiNeg[cent][nh]->Fill(pt,fPIDResponse->NumberOfSigmasTPC(aodTrack,AliPID::kPion));
				hnSigmaKNeg[cent][nh]->Fill(pt,fPIDResponse->NumberOfSigmasTPC(aodTrack,AliPID::kKaon));
				hnSigmaPNeg[cent][nh]->Fill(pt,fPIDResponse->NumberOfSigmasTPC(aodTrack,AliPID::kProton));
			}

			else{
				hPtPos[cent][nh]->Fill(pt);
				hnSigmaPiPos[cent][nh]->Fill(pt,fPIDResponse->NumberOfSigmasTPC(aodTrack,AliPID::kPion));
				hnSigmaKPos[cent][nh]->Fill(pt,fPIDResponse->NumberOfSigmasTPC(aodTrack,AliPID::kKaon));
				hnSigmaPPos[cent][nh]->Fill(pt,fPIDResponse->NumberOfSigmasTPC(aodTrack,AliPID::kProton));
			}

		}//end loop over eta intervals





	}//end of track loop




}



//----------------------------------------------------------------------------------



void AliAnalysisTaskQAHighPtDeDxTest::ProduceArrayV0ESD( AliESDEvent *ESDevent, const Int_t Cent ){

	Int_t nv0s = ESDevent->GetNumberOfV0s();

	fcentAfterV0s->Fill(Cent+1);
	fcentAfterV0s->Fill(11);

	const AliESDVertex *myBestPrimaryVertex = ESDevent->GetPrimaryVertex();

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

	for (Int_t iV0=0; iV0<nv0s; iV0++) {

		AliESDv0 *esdV0 = ESDevent->GetV0(iV0);
		if ( !esdV0 ) continue;

		if( esdV0->GetOnFlyStatus()!=0 )
			continue;

		// AliESDTrack (V0 Daughters)
		UInt_t lKeyPos = (UInt_t)TMath::Abs(esdV0->GetPindex());
		UInt_t lKeyNeg = (UInt_t)TMath::Abs(esdV0->GetNindex());

		AliESDtrack *pTrack = ESDevent->GetTrack(lKeyPos);
		AliESDtrack *nTrack = ESDevent->GetTrack(lKeyNeg);

		if (!pTrack || !nTrack) {
			Printf("ERROR: Could not retreive one of the daughter track");
			continue;
		}

		// Remove like-sign
		if (pTrack->GetSign() == nTrack->GetSign())
			continue;

		UInt_t selectDebug_p = 0;
		if (fTrackFilter2015PbPb) {
			selectDebug_p = fTrackFilter2015PbPb->IsSelected(pTrack);
			if (!selectDebug_p) {
				continue;
			}
		}

		UInt_t selectDebug_n = 0;
		if (fTrackFilter2015PbPb) {
			selectDebug_n = fTrackFilter2015PbPb->IsSelected(nTrack);
			if (!selectDebug_n) {
				continue;
			}
		}

		// Eta cut on decay products
		if(TMath::Abs(pTrack->Eta()) > fEtaCut || TMath::Abs(nTrack->Eta()) > fEtaCut)
			continue;

		// Check if switch does anything!
		Bool_t isSwitched = kFALSE;
		if (pTrack->GetSign() < 0) { // switch
			isSwitched = kTRUE;
			AliESDtrack* helpTrack = nTrack;
			nTrack = pTrack;
			pTrack = helpTrack;
		}

		AliKFVertex primaryVtxKF( *myPrimaryVertex );
		AliKFParticle::SetField(ESDevent->GetMagneticField());

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


						        if(track->GetTPCsignalN()<fNcl)
								continue;
						       
							Double_t phi     = track->Phi();
						    	if(!PhiCut(track->Pt(), phi, track->Charge(), magf, &cutLow, &cutHigh))
							       	continue;

						       Double_t eta      = track->Eta();
						       Double_t momentum = track->P();
						       Double_t dedx     = track->GetTPCsignal();
						       Double_t dedxmb   = track->GetTPCsignal();
						       Double_t dedxUnc  = track->GetTPCsignal();


						       if(fdEdxCalibrated){
							       if(eta < 0){
								       dedx     *= 50.0/EtaCalibrationNeg(Cent,eta);
								       dedxmb   *= 50.0/EtaCalibrationNeg(10,eta);
							       }
							       else{
								       dedx     *= 50.0/EtaCalibrationPos(Cent,eta);
								       dedxmb   *= 50.0/EtaCalibrationPos(10,eta);
							       }
						       }


						       if(fillPos&&fillNeg){
							       if( dedxUnc < DeDxMIPMax && dedxUnc > DeDxMIPMin ){
								       if(momentum<0.6&&momentum>0.4){
									       hMIPVsEtaV0s[Cent]->Fill(eta,dedx);
									       pMIPVsEtaV0s[Cent]->Fill(eta,dedx);

									       hMIPVsEtaV0s[10]->Fill(eta,dedxmb);
									       pMIPVsEtaV0s[10]->Fill(eta,dedxmb);
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

						       if(nh<0)
							       continue;


						       if(fillPos&&fillNeg){

							       histPiV0[Cent][nh]->Fill(momentum, dedx);
							       histpPiV0[Cent][nh]->Fill(momentum);

							       histPiV0[10][nh]->Fill(momentum, dedxmb);
							       histpPiV0[10][nh]->Fill(momentum);

						       }
						       else{
							       histPV0[Cent][nh]->Fill(momentum, dedx);
							       histpPV0[Cent][nh]->Fill(momentum);

							       histPV0[10][nh]->Fill(momentum, dedxmb);
							       histpPV0[10][nh]->Fill(momentum);

						       }
					       }//end loop over two tracks

				       };
				       break;

				case 1:{//gammas

					       Bool_t fillPos = kFALSE;
					       Bool_t fillNeg = kFALSE;


					       if( dmassK>0.01 && dmassL>0.01 && dmassAL>0.01 ) {
						       if( dmassG<0.01 && dmassG>0.0001 ) {

							       if(nTrack->Eta() > 0)
								       if( TMath::Abs(nTrack->GetTPCsignal()-EtaCalibrationPosEl(Cent,nTrack->Eta())) < 5)
									       fillPos = kTRUE;

							       if(nTrack->Eta() < 0)
								       if( TMath::Abs(nTrack->GetTPCsignal()-EtaCalibrationNegEl(Cent,nTrack->Eta())) < 5)
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
					       Double_t dedxmb   = track->GetTPCsignal();
					       Double_t eta      = track->Eta();
					       Double_t phi      = track->Phi();
					       Double_t momentum = track->P();

					       if(fdEdxCalibrated){
						       if(eta < 0){
							       dedx   *= 50/EtaCalibrationNeg(Cent,eta);
							       dedxmb *= 50/EtaCalibrationNeg(10,eta);
									}
						       else{
							       dedx   *= 50/EtaCalibrationPos(Cent,eta);
							       dedxmb *= 50/EtaCalibrationPos(10,eta);
									}
					       }


					       if(track->GetTPCsignalN()<=70)continue;

					       if(!PhiCut(track->Pt(), phi, track->Charge(), magf, &cutLow, &cutHigh))
						       continue;

					       Int_t nh = -1;

					       if(TMath::Abs(eta)<0.2)
						       nh = 0; 
					       else if(TMath::Abs(eta)>=0.2 && TMath::Abs(eta)<0.4)
						       nh = 1; 
					       else if(TMath::Abs(eta)>=0.4 && TMath::Abs(eta)<0.6)
						       nh = 2; 
					       else if(TMath::Abs(eta)>=0.6 && TMath::Abs(eta)<0.8)
						       nh = 3; 

					       if(nh<0)
						       continue;

					       histEV0[Cent][nh]->Fill(momentum, dedx);
					       histEV0[10][nh]->Fill(momentum, dedxmb);

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



//__________________________________________________________________________



void AliAnalysisTaskQAHighPtDeDxTest::ProduceArrayV0AOD( AliAODEvent *AODevent ){
	Int_t nv0s = AODevent->GetNumberOfV0s();
	/*
	   if(nv0s<1){
	   return;
	   }*/

	AliAODVertex *myBestPrimaryVertex = AODevent->GetPrimaryVertex();
	if (!myBestPrimaryVertex) return;



	// ################################
	// #### BEGINNING OF V0 CODE ######
	// ################################
	// This is the begining of the V0 loop
	for (Int_t iV0 = 0; iV0 < nv0s; iV0++) {
		AliAODv0 *aodV0 = AODevent->GetV0(iV0);
		if (!aodV0) continue;


		//check onfly status
		if(aodV0->GetOnFlyStatus()!=0)
			continue;

		// AliAODTrack (V0 Daughters)
		AliAODVertex* vertex = aodV0->GetSecondaryVtx();
		if (!vertex) {
			Printf("ERROR: Could not retrieve vertex");
			continue;
		}

		AliAODTrack *pTrack = (AliAODTrack*)vertex->GetDaughter(0);
		AliAODTrack *nTrack = (AliAODTrack*)vertex->GetDaughter(1);
		if (!pTrack || !nTrack) {
			Printf("ERROR: Could not retrieve one of the daughter track");
			continue;
		}

		// Remove like-sign
		if (pTrack->Charge() == nTrack->Charge()) {
			//cout<< "like sign, continue"<< endl;
			continue;
		}

		// Make sure charge ordering is ok
		if (pTrack->Charge() < 0) {
			AliAODTrack* helpTrack = pTrack;
			pTrack = nTrack;
			nTrack = helpTrack;
		}

		// Eta cut on decay products
		if(TMath::Abs(pTrack->Eta()) > fEtaCut || TMath::Abs(nTrack->Eta()) > fEtaCut)
			continue;


		Double_t dmassG  = aodV0->InvMass2Prongs(0,1,11,11);
		Double_t dmassK  = aodV0->MassK0Short()-0.498;
		Double_t dmassL  = aodV0->MassLambda()-1.116;
		Double_t dmassAL = aodV0->MassAntiLambda()-1.116;

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

						       AliAODTrack* track = 0;

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

						       if(track->GetTPCsignalN()<=70)continue;

						       Double_t phi     = track->Phi();

						       if(!PhiCut(track->Pt(), phi, track->Charge(), magf, &cutLow, &cutHigh))
							       continue;


						       //if(!PhiCut(pt, phi, charge, magf, &cutLow, &cutHigh, hPhi))
						       //	continue;

						       Double_t eta  = track->Eta();
						       Double_t momentum = track->Pt();
						       Double_t dedx = track->GetTPCsignal();

						       if(fillPos&&fillNeg){


							       if( dedx < DeDxMIPMax && dedx > DeDxMIPMin ){
								       if(momentum<0.6&&momentum>0.4){
									       hMIPVsEtaV0s[cent]->Fill(eta,dedx);
									       pMIPVsEtaV0s[cent]->Fill(eta,dedx);
								       }
							       }


						       }

						       for(Int_t nh = 0; nh < nHists; nh++) {



							       if( eta < etaLow[nh]/10.0 || eta > etaHigh[nh]/10.0 )
								       continue;

							       if(fillPos&&fillNeg){

								       histPiV0[cent][nh]->Fill(momentum, dedx);
								       histpPiV0[cent][nh]->Fill(momentum);

							       }
							       else{

								       histPV0[cent][nh]->Fill(momentum, dedx);
								       histpPV0[cent][nh]->Fill(momentum);

							       }

						       }


					       }//end loop over two tracks
				       };
				       break;

				case 1:{//gammas

					       Bool_t fillPos = kFALSE;
					       Bool_t fillNeg = kFALSE;

					       if(dmassK>0.01 && dmassL>0.01 && dmassAL>0.01) {
						       if(dmassG<0.01 && dmassG>0.0001) {

							       if( TMath::Abs(nTrack->GetTPCsignal() - 85.0) < 5)
								       fillPos = kTRUE;
							       if( TMath::Abs(pTrack->GetTPCsignal() - 85.0) < 5)
								       fillNeg = kTRUE;

						       } else {
							       continue;
						       }
					       }


					       if(fillPos == kTRUE && fillNeg == kTRUE)
						       continue;


					       AliAODTrack* track = 0;
					       if(fillNeg)
						       track = nTrack;
					       else if(fillPos)
						       track = pTrack;
					       else
						       continue;

					       Double_t dedx  = track->GetTPCsignal();
					       Double_t eta  = track->Eta();
					       Double_t phi  = track->Phi();
					       Double_t momentum = track->P();

					       if(track->GetTPCsignalN()<=70)continue;

					       if(!PhiCut(track->Pt(), phi, track->Charge(), magf, &cutLow, &cutHigh))
						       continue;

					       for(Int_t nh = 0; nh < nHists; nh++) {

						       if( eta < etaLow[nh]/10.0 || eta > etaHigh[nh]/10.0 )
							       continue;

						       histEV0[cent][nh]->Fill(momentum, dedx);

					       }

				       };
				       break;


			}//end switch
		}//end loop over V0s cases

	}//end loop over v0's




}

//-------------------------------------------------------------------------
Bool_t AliAnalysisTaskQAHighPtDeDxTest::PhiCut(Double_t pt, Double_t phi, Double_t q, Float_t   mag, TF1* phiCutLow, TF1* phiCutHigh)

{
	if(pt < 2.0)
		return kTRUE;

	//Double_t phi = track->Phi();
	if(mag < 0)    // for negatve polarity field
		phi = TMath::TwoPi() - phi;
	if(q < 0) // for negatve charge
		phi = TMath::TwoPi()-phi;

	phi += TMath::Pi()/18.0; // to center gap in the middle
	phi = fmod(phi, TMath::Pi()/9.0);

	if(phi<phiCutHigh->Eval(pt)
			&& phi>phiCutLow->Eval(pt))
		return kFALSE; // reject track

	hPhi[cent]->Fill(pt, phi);
	hPhi[10]->Fill(pt, phi);

	return kTRUE;
}

//-------------------------------------------------------------------------

Float_t AliAnalysisTaskQAHighPtDeDxTest::GetMaxDCApTDep( TF1 *fMaxDCAxy, Double_t ptI){

	Double_t maxDCAxy = 10;
	maxDCAxy = fMaxDCAxy->Eval(ptI);
	return maxDCAxy;

}


//-------------------------------------------------------------------------

Double_t AliAnalysisTaskQAHighPtDeDxTest::EtaCalibrationNeg( const Int_t Cent, const Double_t eta){


	for(Int_t i=0; i<8; ++i)
		fEtaCalibrationNeg->SetParameter(i,0);

	fEtaCalibrationNeg->SetParameter(0,aNeg[Cent]);
	fEtaCalibrationNeg->SetParameter(1,bNeg[Cent]);
	fEtaCalibrationNeg->SetParameter(2,cNeg[Cent]);
	fEtaCalibrationNeg->SetParameter(3,dNeg[Cent]);
	fEtaCalibrationNeg->SetParameter(4,eNeg[Cent]);
	fEtaCalibrationNeg->SetParameter(5,fNeg[Cent]);
	fEtaCalibrationNeg->SetParameter(6,gNeg[Cent]);
	fEtaCalibrationNeg->SetParameter(7,hNeg[Cent]);

	/*
	   cout<<"----------------------------------"<<endl;
	   cout<<"Values InFunction ::  Cent  "<<Cent<<"Eta   "<<eta<<endl;  
	   printf("Par0 = %f \n",aNeg[Cent]);
	   printf("Par1 = %f \n",bNeg[Cent]);
	   printf("Par2 = %f \n",cNeg[Cent]);
	   printf("Par3 = %f \n",dNeg[Cent]);
	   printf("Par4 = %f \n",eNeg[Cent]);
	   printf("Par5 = %f \n",fNeg[Cent]);
	   printf("Par6 = %f \n",gNeg[Cent]);
	   printf("Par7 = %f \n",hNeg[Cent]);
	   printf("f(eta) = %f \n",fEtaCalibrationNeg->Eval(eta));
	 */
	return fEtaCalibrationNeg->Eval(eta);


}


//-------------------------------------------------------------------------

Double_t AliAnalysisTaskQAHighPtDeDxTest::EtaCalibrationPos( const Int_t Cent, const Double_t eta){


	for(Int_t i=0; i<8; ++i)
		fEtaCalibration->SetParameter(i,0);

	fEtaCalibration->SetParameter(0,aPos[Cent]);
	fEtaCalibration->SetParameter(1,bPos[Cent]);
	fEtaCalibration->SetParameter(2,cPos[Cent]);
	fEtaCalibration->SetParameter(3,dPos[Cent]);
	fEtaCalibration->SetParameter(4,ePos[Cent]);
	fEtaCalibration->SetParameter(5,fPos[Cent]);
	fEtaCalibration->SetParameter(6,gPos[Cent]);
	fEtaCalibration->SetParameter(7,hPos[Cent]);


	/*	cout<<"----------------------------------"<<endl;
		cout<<"Values InFunction ::  Cent  "<<Cent<<"Eta   "<<eta<<endl;  
		printf("Par0 = %f \n",aPos[Cent]);
		printf("Par1 = %f \n",bPos[Cent]);
		printf("Par2 = %f \n",cPos[Cent]);
		printf("Par3 = %f \n",dPos[Cent]);
		printf("Par4 = %f \n",ePos[Cent]);
		printf("Par5 = %f \n",fPos[Cent]);
		printf("Par6 = %f \n",gPos[Cent]);
		printf("Par7 = %f \n",hPos[Cent]);
		printf("f(eta) = %f \n",fEtaCalibration->Eval(eta));
	 */
	return fEtaCalibration->Eval(eta);

}



//-------------------------------------------------------------------------

Double_t AliAnalysisTaskQAHighPtDeDxTest::EtaCalibrationNegEl(const Int_t Cent, const Double_t eta){


	for(Int_t i=0; i<5; ++i)
		felededxfitNeg->SetParameter(i,0);


	felededxfitNeg->SetParameter(0,aNegEl[Cent]);
	felededxfitNeg->SetParameter(1,bNegEl[Cent]);
	felededxfitNeg->SetParameter(2,cNegEl[Cent]);
	felededxfitNeg->SetParameter(3,dNegEl[Cent]);
	felededxfitNeg->SetParameter(4,eNegEl[Cent]);


	return felededxfitNeg->Eval(eta);

}


//-------------------------------------------------------------------------

Double_t AliAnalysisTaskQAHighPtDeDxTest::EtaCalibrationPosEl(const Int_t Cent, const Double_t eta){


	for(Int_t i=0; i<5; ++i)
		felededxfitPos->SetParameter(i,0);

	felededxfitPos->SetParameter(0,aPosEl[Cent]);
	felededxfitPos->SetParameter(1,bPosEl[Cent]);
	felededxfitPos->SetParameter(2,cPosEl[Cent]);
	felededxfitPos->SetParameter(3,dPosEl[Cent]);
	felededxfitPos->SetParameter(4,ePosEl[Cent]);

	return felededxfitPos->Eval(eta);

}

