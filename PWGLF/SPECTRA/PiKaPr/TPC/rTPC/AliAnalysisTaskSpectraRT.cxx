/**************************************************************************
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
 *************************************************************************/

#include "AliAnalysisTaskSpectraRT.h"

// ROOT includes
#include <TList.h>
#include <TChain.h>
#include <TTree.h>
#include <TMath.h>
#include <TRandom.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
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
#include <AliMultiplicity.h>
#include <AliMultSelection.h>
#include <AliPIDResponse.h>
#include "AliTPCPIDResponse.h"
#include "AliAnalysisUtils.h"

#include <AliMCEventHandler.h>
#include <AliMCEvent.h>
#include "AliMCParticle.h"
#include <AliStack.h>

#include <TTreeStream.h>

#include <AliHeader.h>
#include <AliGenPythiaEventHeader.h>
#include <AliGenDPMjetEventHeader.h>
#include "AliGenEventHeader.h"
#include "AliGenCocktailEventHeader.h"

#include <AliCentrality.h>
#include <AliESDv0.h>
#include <AliKFVertex.h>
#include <AliAODVertex.h>

#include <AliAODTrack.h>
#include <AliVParticle.h>
#include <AliPID.h>
#include <AliAODPid.h>
#include <AliAODMCHeader.h>

#include <iostream>
#include <vector>
////class AliAnalysisTaskSpectraRT;    // your analysis class
using namespace std;


//
// Responsible:
// Omar Vazquez (Lund)


static float Magf                = 1;
static const int nPid              = 4;
static const int nRegion           = 4;
static const int nHists            = 4;
const char* Region[4] = {"Toward","Away","Transverse","FullAzimuth"};
const char* Pid[nPid] = {"Charged","Pion","Kaon","Proton"};
const char* Charge[2] = {"Pos","Neg"};
static const double C_Value = TMath::C()*(1.e2/1.e12); // cm/ps

ClassImp(AliAnalysisTaskSpectraRT)
	AliAnalysisTaskSpectraRT::AliAnalysisTaskSpectraRT():
		AliAnalysisTaskSE(),
		fESD(0x0),
		fEventCuts(0x0),
		fMC(0x0),
		fMCStack(0x0),
		fMCArray(0x0),
		fPIDResponse(0x0),
		fGeometricalCut(0x0),
		fTrackFilterDaughters(0x0),
		fTrackFilter(0x0),
		fHybridTrackCuts1(0x0),
		fHybridTrackCuts2(0x0),
		utils(0x0),
		fAnalysisType("ESD"),
		fAnalysisMC(kFALSE),
		fIsMCclosure(kTRUE),
		fRandom(0x0),
		fNcl(70),
		fTrackID(0),
		fEtaCut(0.9),
		fdEdxCalibrated(kTRUE),
		fDeDxMIPMin(40),
		fDeDxMIPMax(60),
		fdEdxHigh(200),
		fdEdxLow(40),
		fPeriod("l"),
		///		fSetTPConlyTrkCuts(kFALSE),
		fSelectHybridTracks(kTRUE),
		fLeadPtCutMin(5.0),
		fLeadPtCutMax(40.0),
		fGenLeadPhi(0.0),
		fGenLeadPt(0.0),
		fGenLeadIn(0.0),
		fRecLeadEta(0.0),
		fRecLeadPhi(0.0),
		fRecLeadPt(0.0),
		fRecLeadIn(0.0),
		fPtMin(0.15),
		fListOfObjects(0),
		fEvents(0x0),
		hNchTSData(0x0),
		hPhiTotal(0x0),
		hPhiStandard(0x0),
		hPhiHybrid1(0x0),
		hPhiHybrid2(0x0),
		hPhiLeading(0x0),
		hDeltaPhiDeltaEta(0x0),
		fEtaCalibrationPos(0x0),
		fEtaCalibrationNeg(0x0),
		fEtaCalibrationPosEl(0x0),
		fEtaCalibrationNegEl(0x0),
		fcutDCAxy(0x0),
		fcutLow(0x0),
		fcutHigh(0x0),
		hMIPVsEta(0x0),
		pMIPVsEta(0x0),
		hPlateauVsEta(0x0),
		pPlateauVsEta(0x0),
		hMIPVsEtaV0s(0x0),
		pMIPVsEtaV0s(0x0),
		hPhirTPC(0x0)

{


	for(int r = 0; r < nRegion; r++){

		if(r < (nRegion-1)){
			hPhiData[r] = 0;
		}

		for(int j = 0; j < nHists; j++){

			if(r < (nRegion-1)){
				hNchVsPtDataPosTOF[r][j] = 0;
				hNchVsPtDataNegTOF[r][j] = 0;
				hNchVsPtPosTPC[r][j] = 0;
				hNchVsPtNegTPC[r][j] = 0;
				hNchVsPtPosTOF[r][j] = 0;
				hNchVsPtNegTOF[r][j] = 0;
				hNchVsPPosTOF[r][j] = 0;
				hNchVsPNegTOF[r][j] = 0;
				hDeDxVsP[r][j] = 0;
				hNchVsPrTPC[r][j] = 0;
				hNchVsPtrTPC[r][j] = 0;
			}

		}	// ending

	}	// region        


	for(int r = 0; r < (nRegion-1); r++){

		for(int j = 0; j < nHists; j++){

			hNchVsPtDataPosPionTPC[r][j] = 0;
			hNchVsPtDataNegPionTPC[r][j] = 0;
			hNchVsPtDataPosKaonTPC[r][j] = 0;
			hNchVsPtDataNegKaonTPC[r][j] = 0;
			hNchVsPtDataPosProtonTPC[r][j] = 0;
			hNchVsPtDataNegProtonTPC[r][j] = 0;
		}	// ending

	}	// region        


	for(int j = 0; j < nHists; ++j){
		histPiTof[j] = 0;
		histEV0[j] = 0;
		histPV0[j] = 0;
		histPiV0[j] = 0;
		hMIPVsPhi[j] = 0;
		pMIPVsPhi[j] = 0;
		hPlateauVsPhi[j] = 0;
		pPlateauVsPhi[j] = 0;
		hPtVsP[j] = 0;
		hnSigmaElectrons[j] = 0;

	}

}


AliAnalysisTaskSpectraRT::AliAnalysisTaskSpectraRT(const char *name):
	AliAnalysisTaskSE(name),
	fESD(0x0),
	fEventCuts(0x0),
	fMC(0x0),
	fMCStack(0x0),
	fMCArray(0x0),
	fPIDResponse(0x0),
	fGeometricalCut(0x0),
	fTrackFilterDaughters(0x0),
	fTrackFilter(0x0),
	fHybridTrackCuts1(0x0),
	fHybridTrackCuts2(0x0),
	utils(0x0),
	fAnalysisType("ESD"),
	fAnalysisMC(kFALSE),
	fIsMCclosure(kTRUE),
	fRandom(0x0),
	fNcl(70),
	fTrackID(0),
	fEtaCut(0.9),
	fdEdxCalibrated(kTRUE),
	fDeDxMIPMin(40),
	fDeDxMIPMax(60),
	fdEdxHigh(200),
	fdEdxLow(40),
	fPeriod("l"),
	///	fSetTPConlyTrkCuts(kFALSE),
	fSelectHybridTracks(kTRUE),
	fLeadPtCutMin(5.0),
	fLeadPtCutMax(40.0),
	fGenLeadPhi(0.0),
	fGenLeadPt(0.0),
	fGenLeadIn(0.0),
	fRecLeadEta(0.0),
	fRecLeadPhi(0.0),
	fRecLeadPt(0.0),
	fRecLeadIn(0.0),
	fPtMin(0.15),
	fListOfObjects(0),
	fEvents(0x0),
	hNchTSData(0x0),
	hPhiTotal(0x0),
	hPhiStandard(0x0),
	hPhiHybrid1(0x0),
	hPhiHybrid2(0x0),
	hPhiLeading(0x0),
	hDeltaPhiDeltaEta(0x0),
	fEtaCalibrationPos(0x0),
	fEtaCalibrationNeg(0x0),
	fEtaCalibrationPosEl(0x0),
	fEtaCalibrationNegEl(0x0),
	fcutDCAxy(0x0),
	fcutLow(0x0),
	fcutHigh(0x0),
	hMIPVsEta(0x0),
	pMIPVsEta(0x0),
	hPlateauVsEta(0x0),
	pPlateauVsEta(0x0),
	hMIPVsEtaV0s(0x0),
	pMIPVsEtaV0s(0x0),
	hPhirTPC(0x0)

{

	for(int r = 0; r < nRegion; r++){

		if(r < (nRegion-1)){
			hPhiData[r] = 0;
		}

		for(int j = 0; j < nHists; j++){

			if(r < (nRegion-1)){
				hNchVsPtDataPosTOF[r][j] = 0;
				hNchVsPtDataNegTOF[r][j] = 0;
				hNchVsPtPosTPC[r][j] = 0;
				hNchVsPtNegTPC[r][j] = 0;
				hNchVsPtPosTOF[r][j] = 0;
				hNchVsPtNegTOF[r][j] = 0;
				hNchVsPPosTOF[r][j] = 0;
				hNchVsPNegTOF[r][j] = 0;
				hDeDxVsP[r][j] = 0;
				hNchVsPrTPC[r][j] = 0;
				hNchVsPtrTPC[r][j] = 0;
			}
		}

	}        

	for(int r = 0; r < (nRegion-1); r++){

		for(int j = 0; j < nHists; j++){

			hNchVsPtDataPosPionTPC[r][j] = 0;
			hNchVsPtDataNegPionTPC[r][j] = 0;
			hNchVsPtDataPosKaonTPC[r][j] = 0;
			hNchVsPtDataNegKaonTPC[r][j] = 0;
			hNchVsPtDataPosProtonTPC[r][j] = 0;
			hNchVsPtDataNegProtonTPC[r][j] = 0;
		}	// ending

	}	// region        


	for(int j = 0; j < nHists; ++j){
		histPiTof[j] = 0;
		histEV0[j] = 0;
		histPV0[j] = 0;
		histPiV0[j] = 0;
		hMIPVsPhi[j] = 0;
		pMIPVsPhi[j] = 0;
		hPlateauVsPhi[j] = 0;
		pPlateauVsPhi[j] = 0;
		hPtVsP[j] = 0;
		hnSigmaElectrons[j] = 0;

	}


	DefineInput(0, TChain::Class());
	DefineOutput(1, TList::Class());//esto es nuevo
}

AliAnalysisTaskSpectraRT::~AliAnalysisTaskSpectraRT() {
	//
	// Destructor
	//

}
//______________________________________________________________________________
void AliAnalysisTaskSpectraRT::UserCreateOutputObjects()
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

	// Quality cuts for selecting the leading particle
	// Hybrid tracks + Geometrical cut
	if(!fGeometricalCut){

		fGeometricalCut = new AliESDtrackCuts("fGeometricalCut");	

		if(fTrackID==11) { fGeometricalCut->SetCutGeoNcrNcl(2, 130, 1.5, 0.85, 0.7); printf("fTrackID = %d\n",fTrackID);}
		else if(fTrackID==12) { fGeometricalCut->SetCutGeoNcrNcl(4, 130, 1.5, 0.85, 0.7); printf("fTrackID = %d\n",fTrackID);}
		else if(fTrackID==13) { fGeometricalCut->SetCutGeoNcrNcl(3, 120, 1.5, 0.85, 0.7); printf("fTrackID = %d\n",fTrackID);}
		else if(fTrackID==14) { fGeometricalCut->SetCutGeoNcrNcl(3, 140, 1.5, 0.85, 0.7); printf("fTrackID = %d\n",fTrackID);}
		else{ 
			printf("Nominal setting for the GeometricalCut - fTrackID = %d\n",fTrackID);
			fGeometricalCut->SetCutGeoNcrNcl(3, 130, 1.5, 0.85, 0.7); 
		}

	}

	// Track Cuts for Nch in the Transverse region and pT spectra
	// Hybrid tracks
	if(!fTrackFilter){

		fTrackFilter = new AliESDtrackCuts("fTrackFilter");	
		//fTrackFilter->SetMinNCrossedRowsTPC(70); //! Variated in track cuts systematics
		//fTrackFilter->SetMinRatioCrossedRowsOverFindableClustersTPC(0.8); //! Variated in track cuts systematics
		//fTrackFilter->SetMaxChi2PerClusterTPC(4); //! Variated in track cuts systematics
		fTrackFilter->SetAcceptKinkDaughters(kFALSE);
		fTrackFilter->SetRequireTPCRefit(kTRUE);
		fTrackFilter->SetRequireITSRefit(kTRUE);
		fTrackFilter->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kAny);
		fTrackFilter->SetMaxDCAToVertexXYPtDep("0.0105+0.0350/pt^1.1");
		//fTrackFilter->SetMaxChi2TPCConstrainedGlobal(36); //! This cut is excluded 
		//fTrackFilter->SetMaxDCAToVertexZ(2); //! Variated in track cuts systematics
		fTrackFilter->SetDCAToVertex2D(kFALSE);
		fTrackFilter->SetRequireSigmaToVertex(kFALSE);
		//fTrackFilter->SetMaxChi2PerClusterITS(36); //! Variated in track cuts systematics
		fTrackFilter->SetEtaRange(-0.8,0.8);

		if(fTrackID==0){ //! Nominal values
			fTrackFilter->SetMinNCrossedRowsTPC(70);
			fTrackFilter->SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);
			fTrackFilter->SetMaxChi2PerClusterTPC(4);
			fTrackFilter->SetMaxDCAToVertexZ(2);
			fTrackFilter->SetMaxChi2PerClusterITS(36);
		}
		if(fTrackID==1){ //! Lower: SetMinNCrossedRowsTPC(60)
			fTrackFilter->SetMinNCrossedRowsTPC(60);
			fTrackFilter->SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);
			fTrackFilter->SetMaxChi2PerClusterTPC(4);
			fTrackFilter->SetMaxDCAToVertexZ(2);
			fTrackFilter->SetMaxChi2PerClusterITS(36);
		}
		if(fTrackID==2){ //! Higher: SetMinNCrossedRowsTPC(100)
			printf("fTrackFilter for fTrackID = %d\n",fTrackID);
			fTrackFilter->SetMinNCrossedRowsTPC(100);
			fTrackFilter->SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);
			fTrackFilter->SetMaxChi2PerClusterTPC(4);
			fTrackFilter->SetMaxDCAToVertexZ(2);
			fTrackFilter->SetMaxChi2PerClusterITS(36);
		}
		if(fTrackID==3){ //! Lower: SetMinRatioCrossedRowsOverFindableClustersTPC(0.7) 
			fTrackFilter->SetMinNCrossedRowsTPC(70);
			fTrackFilter->SetMinRatioCrossedRowsOverFindableClustersTPC(0.7);
			fTrackFilter->SetMaxChi2PerClusterTPC(4);
			fTrackFilter->SetMaxDCAToVertexZ(2);
			fTrackFilter->SetMaxChi2PerClusterITS(36);
		}
		if(fTrackID==4){ //! Higher: SetMinRatioCrossedRowsOverFindableClustersTPC(0.9)
			fTrackFilter->SetMinNCrossedRowsTPC(70);
			fTrackFilter->SetMinRatioCrossedRowsOverFindableClustersTPC(0.9);
			fTrackFilter->SetMaxChi2PerClusterTPC(4);
			fTrackFilter->SetMaxDCAToVertexZ(2);
			fTrackFilter->SetMaxChi2PerClusterITS(36);
		}
		if(fTrackID==5){ //! Lower: SetMaxChi2PerClusterTPC(3)  
			fTrackFilter->SetMinNCrossedRowsTPC(70);
			fTrackFilter->SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);
			fTrackFilter->SetMaxChi2PerClusterTPC(3);
			fTrackFilter->SetMaxDCAToVertexZ(2);
			fTrackFilter->SetMaxChi2PerClusterITS(36);
		}
		if(fTrackID==6){ //! Higher: SetMaxChi2PerClusterTPC(5)
			fTrackFilter->SetMinNCrossedRowsTPC(70);
			fTrackFilter->SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);
			fTrackFilter->SetMaxChi2PerClusterTPC(5);
			fTrackFilter->SetMaxDCAToVertexZ(2);
			fTrackFilter->SetMaxChi2PerClusterITS(36);
		}
		if(fTrackID==7){ //! Lower: SetMaxChi2PerClusterITS(25)
			fTrackFilter->SetMinNCrossedRowsTPC(70);
			fTrackFilter->SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);
			fTrackFilter->SetMaxChi2PerClusterTPC(4);
			fTrackFilter->SetMaxDCAToVertexZ(2);
			fTrackFilter->SetMaxChi2PerClusterITS(25);
		}
		if(fTrackID==8){ //! Higher: SetMaxChi2PerClusterITS(49)
			fTrackFilter->SetMinNCrossedRowsTPC(70);
			fTrackFilter->SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);
			fTrackFilter->SetMaxChi2PerClusterTPC(4);
			fTrackFilter->SetMaxDCAToVertexZ(2);
			fTrackFilter->SetMaxChi2PerClusterITS(49);
		}
		if(fTrackID==9){ //! Lower: SetMaxDCAToVertexZ(1)
			fTrackFilter->SetMinNCrossedRowsTPC(70);
			fTrackFilter->SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);
			fTrackFilter->SetMaxChi2PerClusterTPC(4);
			fTrackFilter->SetMaxDCAToVertexZ(1);
			fTrackFilter->SetMaxChi2PerClusterITS(36);
		}
		if(fTrackID==10){ //! Lower: SetMaxDCAToVertexZ(5)
			fTrackFilter->SetMinNCrossedRowsTPC(70);
			fTrackFilter->SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);
			fTrackFilter->SetMaxChi2PerClusterTPC(4);
			fTrackFilter->SetMaxDCAToVertexZ(5);
			fTrackFilter->SetMaxChi2PerClusterITS(36);
		}
		if((fTrackID==11)||(fTrackID==12)||(fTrackID==13)||(fTrackID==14)){ //! Nominal values
			printf("fTrackFilter for fTrackID = %d\n",fTrackID);
			fTrackFilter->SetMinNCrossedRowsTPC(70);
			fTrackFilter->SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);
			fTrackFilter->SetMaxChi2PerClusterTPC(4);
			fTrackFilter->SetMaxDCAToVertexZ(2);
			fTrackFilter->SetMaxChi2PerClusterITS(36);
		}
	}

	if(!fHybridTrackCuts1){
		fHybridTrackCuts1 = new AliESDtrackCuts("fHybridTrackCuts1");	
		//fHybridTrackCuts1->SetMinNCrossedRowsTPC(70); //! Variated in track cuts systematics
		fHybridTrackCuts1->SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);
		//fHybridTrackCuts1->SetMaxChi2PerClusterTPC(4); //! Variated in track cuts systematics
		fHybridTrackCuts1->SetAcceptKinkDaughters(kFALSE);
		fHybridTrackCuts1->SetRequireTPCRefit(kTRUE);
		fHybridTrackCuts1->SetRequireITSRefit(kTRUE);
		fHybridTrackCuts1->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kOff);
		fHybridTrackCuts1->SetMaxDCAToVertexXYPtDep("0.0105+0.0350/pt^1.1");
		//fHybridTrackCuts1->SetMaxChi2TPCConstrainedGlobal(36); //! This cut is excluded 
		//fHybridTrackCuts1->SetMaxDCAToVertexZ(2); //! Variated in track cuts systematics
		fHybridTrackCuts1->SetDCAToVertex2D(kFALSE);
		fHybridTrackCuts1->SetRequireSigmaToVertex(kFALSE);
		fHybridTrackCuts1->SetMaxChi2PerClusterITS(36);
		fHybridTrackCuts1->SetEtaRange(-0.8,0.8);

		if(fTrackID==0){ //! Nominal values
			fHybridTrackCuts1->SetMinNCrossedRowsTPC(70);
			fHybridTrackCuts1->SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);
			fHybridTrackCuts1->SetMaxChi2PerClusterTPC(4);
			fHybridTrackCuts1->SetMaxDCAToVertexZ(2);
			fHybridTrackCuts1->SetMaxChi2PerClusterITS(36);
		}
		if(fTrackID==1){ //! Lower: SetMinNCrossedRowsTPC(60)
			fHybridTrackCuts1->SetMinNCrossedRowsTPC(60);
			fHybridTrackCuts1->SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);
			fHybridTrackCuts1->SetMaxChi2PerClusterTPC(4);
			fHybridTrackCuts1->SetMaxDCAToVertexZ(2);
			fHybridTrackCuts1->SetMaxChi2PerClusterITS(36);
		}
		if(fTrackID==2){ //! Higher: SetMinNCrossedRowsTPC(100)
			fHybridTrackCuts1->SetMinNCrossedRowsTPC(100);
			fHybridTrackCuts1->SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);
			fHybridTrackCuts1->SetMaxChi2PerClusterTPC(4);
			fHybridTrackCuts1->SetMaxDCAToVertexZ(2);
			fHybridTrackCuts1->SetMaxChi2PerClusterITS(36);
		}
		if(fTrackID==3){ //! Lower: SetMinRatioCrossedRowsOverFindableClustersTPC(0.7) 
			fHybridTrackCuts1->SetMinNCrossedRowsTPC(70);
			fHybridTrackCuts1->SetMinRatioCrossedRowsOverFindableClustersTPC(0.7);
			fHybridTrackCuts1->SetMaxChi2PerClusterTPC(4);
			fHybridTrackCuts1->SetMaxDCAToVertexZ(2);
			fHybridTrackCuts1->SetMaxChi2PerClusterITS(36);
		}
		if(fTrackID==4){ //! Higher: SetMinRatioCrossedRowsOverFindableClustersTPC(0.9)
			fHybridTrackCuts1->SetMinNCrossedRowsTPC(70);
			fHybridTrackCuts1->SetMinRatioCrossedRowsOverFindableClustersTPC(0.9);
			fHybridTrackCuts1->SetMaxChi2PerClusterTPC(4);
			fHybridTrackCuts1->SetMaxDCAToVertexZ(2);
			fHybridTrackCuts1->SetMaxChi2PerClusterITS(36);
		}
		if(fTrackID==5){ //! Lower: SetMaxChi2PerClusterTPC(3)  
			fHybridTrackCuts1->SetMinNCrossedRowsTPC(70);
			fHybridTrackCuts1->SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);
			fHybridTrackCuts1->SetMaxChi2PerClusterTPC(3);
			fHybridTrackCuts1->SetMaxDCAToVertexZ(2);
			fHybridTrackCuts1->SetMaxChi2PerClusterITS(36);
		}
		if(fTrackID==6){ //! Higher: SetMaxChi2PerClusterTPC(5)
			fHybridTrackCuts1->SetMinNCrossedRowsTPC(70);
			fHybridTrackCuts1->SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);
			fHybridTrackCuts1->SetMaxChi2PerClusterTPC(5);
			fHybridTrackCuts1->SetMaxDCAToVertexZ(2);
			fHybridTrackCuts1->SetMaxChi2PerClusterITS(36);
		}
		if(fTrackID==7){ //! Lower: SetMaxChi2PerClusterITS(25)
			fHybridTrackCuts1->SetMinNCrossedRowsTPC(70);
			fHybridTrackCuts1->SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);
			fHybridTrackCuts1->SetMaxChi2PerClusterTPC(4);
			fHybridTrackCuts1->SetMaxDCAToVertexZ(2);
			fHybridTrackCuts1->SetMaxChi2PerClusterITS(25);
		}
		if(fTrackID==8){ //! Higher: SetMaxChi2PerClusterITS(49)
			fHybridTrackCuts1->SetMinNCrossedRowsTPC(70);
			fHybridTrackCuts1->SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);
			fHybridTrackCuts1->SetMaxChi2PerClusterTPC(4);
			fHybridTrackCuts1->SetMaxDCAToVertexZ(2);
			fHybridTrackCuts1->SetMaxChi2PerClusterITS(49);
		}
		if(fTrackID==9){ //! Lower: SetMaxDCAToVertexZ(1)
			fHybridTrackCuts1->SetMinNCrossedRowsTPC(70);
			fHybridTrackCuts1->SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);
			fHybridTrackCuts1->SetMaxChi2PerClusterTPC(4);
			fHybridTrackCuts1->SetMaxDCAToVertexZ(1);
			fHybridTrackCuts1->SetMaxChi2PerClusterITS(36);
		}
		if(fTrackID==10){ //! Lower: SetMaxDCAToVertexZ(5)
			fHybridTrackCuts1->SetMinNCrossedRowsTPC(70);
			fHybridTrackCuts1->SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);
			fHybridTrackCuts1->SetMaxChi2PerClusterTPC(4);
			fHybridTrackCuts1->SetMaxDCAToVertexZ(5);
			fHybridTrackCuts1->SetMaxChi2PerClusterITS(36);
		}
		if((fTrackID==11)||(fTrackID==12)||(fTrackID==13)||(fTrackID==14)){ //! Nominal values
			printf("fHybridTrackCuts1 for fTrackID = %d\n",fTrackID);
			fHybridTrackCuts1->SetMinNCrossedRowsTPC(70);
			fHybridTrackCuts1->SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);
			fHybridTrackCuts1->SetMaxChi2PerClusterTPC(4);
			fHybridTrackCuts1->SetMaxDCAToVertexZ(2);
			fHybridTrackCuts1->SetMaxChi2PerClusterITS(36);
		}
	} 

	if(!fHybridTrackCuts2){
		fHybridTrackCuts2 = new AliESDtrackCuts("fHybridTrackCuts2");	
		//fHybridTrackCuts2->SetMinNCrossedRowsTPC(70); //! Variated in track cuts systematics
		fHybridTrackCuts2->SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);
		//fHybridTrackCuts2->SetMaxChi2PerClusterTPC(4); //! Variated in track cuts systematics
		fHybridTrackCuts2->SetAcceptKinkDaughters(kFALSE);
		fHybridTrackCuts2->SetRequireTPCRefit(kTRUE);
		fHybridTrackCuts2->SetRequireITSRefit(kFALSE);
		fHybridTrackCuts2->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kNone);
		fHybridTrackCuts2->SetMaxDCAToVertexXYPtDep("0.0105+0.0350/pt^1.1");
		//fHybridTrackCuts2->SetMaxChi2TPCConstrainedGlobal(36); //! This cut is excluded 
		//fHybridTrackCuts2->SetMaxDCAToVertexZ(2); //! Variated in track cuts systematics
		fHybridTrackCuts2->SetDCAToVertex2D(kFALSE);
		fHybridTrackCuts2->SetRequireSigmaToVertex(kFALSE);
		fHybridTrackCuts2->SetMaxChi2PerClusterITS(36);
		fHybridTrackCuts2->SetEtaRange(-0.8,0.8);

		if(fTrackID==0){ //! Nominal values
			fHybridTrackCuts2->SetMinNCrossedRowsTPC(70);
			fHybridTrackCuts2->SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);
			fHybridTrackCuts2->SetMaxChi2PerClusterTPC(4);
			fHybridTrackCuts2->SetMaxDCAToVertexZ(2);
			fHybridTrackCuts2->SetMaxChi2PerClusterITS(36);
		}
		if(fTrackID==1){ //! Lower: SetMinNCrossedRowsTPC(60)
			fHybridTrackCuts2->SetMinNCrossedRowsTPC(60);
			fHybridTrackCuts2->SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);
			fHybridTrackCuts2->SetMaxChi2PerClusterTPC(4);
			fHybridTrackCuts2->SetMaxDCAToVertexZ(2);
			fHybridTrackCuts2->SetMaxChi2PerClusterITS(36);
		}
		if(fTrackID==2){ //! Higher: SetMinNCrossedRowsTPC(100)
			fHybridTrackCuts2->SetMinNCrossedRowsTPC(100);
			fHybridTrackCuts2->SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);
			fHybridTrackCuts2->SetMaxChi2PerClusterTPC(4);
			fHybridTrackCuts2->SetMaxDCAToVertexZ(2);
			fHybridTrackCuts2->SetMaxChi2PerClusterITS(36);
		}
		if(fTrackID==3){ //! Lower: SetMinRatioCrossedRowsOverFindableClustersTPC(0.7) 
			fHybridTrackCuts2->SetMinNCrossedRowsTPC(70);
			fHybridTrackCuts2->SetMinRatioCrossedRowsOverFindableClustersTPC(0.7);
			fHybridTrackCuts2->SetMaxChi2PerClusterTPC(4);
			fHybridTrackCuts2->SetMaxDCAToVertexZ(2);
			fHybridTrackCuts2->SetMaxChi2PerClusterITS(36);
		}
		if(fTrackID==4){ //! Higher: SetMinRatioCrossedRowsOverFindableClustersTPC(0.9)
			fHybridTrackCuts2->SetMinNCrossedRowsTPC(70);
			fHybridTrackCuts2->SetMinRatioCrossedRowsOverFindableClustersTPC(0.9);
			fHybridTrackCuts2->SetMaxChi2PerClusterTPC(4);
			fHybridTrackCuts2->SetMaxDCAToVertexZ(2);
			fHybridTrackCuts2->SetMaxChi2PerClusterITS(36);
		}
		if(fTrackID==5){ //! Lower: SetMaxChi2PerClusterTPC(3)  
			fHybridTrackCuts2->SetMinNCrossedRowsTPC(70);
			fHybridTrackCuts2->SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);
			fHybridTrackCuts2->SetMaxChi2PerClusterTPC(3);
			fHybridTrackCuts2->SetMaxDCAToVertexZ(2);
			fHybridTrackCuts2->SetMaxChi2PerClusterITS(36);
		}
		if(fTrackID==6){ //! Higher: SetMaxChi2PerClusterTPC(5)
			fHybridTrackCuts2->SetMinNCrossedRowsTPC(70);
			fHybridTrackCuts2->SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);
			fHybridTrackCuts2->SetMaxChi2PerClusterTPC(5);
			fHybridTrackCuts2->SetMaxDCAToVertexZ(2);
			fHybridTrackCuts2->SetMaxChi2PerClusterITS(36);
		}
		if(fTrackID==7){ //! Lower: SetMaxChi2PerClusterITS(25)
			fHybridTrackCuts2->SetMinNCrossedRowsTPC(70);
			fHybridTrackCuts2->SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);
			fHybridTrackCuts2->SetMaxChi2PerClusterTPC(4);
			fHybridTrackCuts2->SetMaxDCAToVertexZ(2);
			fHybridTrackCuts2->SetMaxChi2PerClusterITS(25);
		}
		if(fTrackID==8){ //! Higher: SetMaxChi2PerClusterITS(49)
			fHybridTrackCuts2->SetMinNCrossedRowsTPC(70);
			fHybridTrackCuts2->SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);
			fHybridTrackCuts2->SetMaxChi2PerClusterTPC(4);
			fHybridTrackCuts2->SetMaxDCAToVertexZ(2);
			fHybridTrackCuts2->SetMaxChi2PerClusterITS(49);
		}
		if(fTrackID==9){ //! Lower: SetMaxDCAToVertexZ(1)
			fHybridTrackCuts2->SetMinNCrossedRowsTPC(70);
			fHybridTrackCuts2->SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);
			fHybridTrackCuts2->SetMaxChi2PerClusterTPC(4);
			fHybridTrackCuts2->SetMaxDCAToVertexZ(1);
			fHybridTrackCuts2->SetMaxChi2PerClusterITS(36);
		}
		if(fTrackID==10){ //! Lower: SetMaxDCAToVertexZ(5)
			fHybridTrackCuts2->SetMinNCrossedRowsTPC(70);
			fHybridTrackCuts2->SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);
			fHybridTrackCuts2->SetMaxChi2PerClusterTPC(4);
			fHybridTrackCuts2->SetMaxDCAToVertexZ(5);
			fHybridTrackCuts2->SetMaxChi2PerClusterITS(36);
		}
		if((fTrackID==11)||(fTrackID==12)||(fTrackID==13)||(fTrackID==14)){ //! Nominal values
			printf("fHybridTrackCuts2 for fTrackID = %d\n",fTrackID);
			fHybridTrackCuts2->SetMinNCrossedRowsTPC(70);
			fHybridTrackCuts2->SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);
			fHybridTrackCuts2->SetMaxChi2PerClusterTPC(4);
			fHybridTrackCuts2->SetMaxDCAToVertexZ(2);
			fHybridTrackCuts2->SetMaxChi2PerClusterITS(36);
		}
	} 

	// Quality cuts for selecting daughters of V0s
	if(!fTrackFilterDaughters){

		fTrackFilterDaughters = new AliESDtrackCuts("fTrackFilterDaughters");

		// TPC
		//fTrackFilterDaughters->SetMinNCrossedRowsTPC(70); //! Variated in track cuts systematics
		//fTrackFilterDaughters->SetMinRatioCrossedRowsOverFindableClustersTPC(0.8); //! Variated in track cuts systematics
		//fTrackFilterDaughters->SetMaxChi2PerClusterTPC(4); //! Variated in track cuts systematics
		fTrackFilterDaughters->SetAcceptKinkDaughters(kFALSE);
		fTrackFilterDaughters->SetRequireTPCRefit(kTRUE);
		// ITS
		fTrackFilterDaughters->SetRequireITSRefit(kTRUE);
		fTrackFilterDaughters->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kAny);
		//fTrackFilterDaughters->SetMaxChi2TPCConstrainedGlobal(36); //! This cut is excluded 
		//fTrackFilterDaughters->SetMaxDCAToVertexZ(2); //! Variated in track cuts systematics
		fTrackFilterDaughters->SetDCAToVertex2D(kFALSE);
		fTrackFilterDaughters->SetRequireSigmaToVertex(kFALSE);
		//fTrackFilterDaughters->SetMaxChi2PerClusterITS(36); //! Variated in track cuts systematics
		fTrackFilterDaughters->SetEtaRange(-0.8,0.8);

		if(fTrackID==0){ //! Nominal values
			fTrackFilterDaughters->SetMinNCrossedRowsTPC(70);
			fTrackFilterDaughters->SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);
			fTrackFilterDaughters->SetMaxChi2PerClusterTPC(4);
			fTrackFilterDaughters->SetMaxDCAToVertexZ(2);
			fTrackFilterDaughters->SetMaxChi2PerClusterITS(36);
		}
		if(fTrackID==1){ //! Lower: SetMinNCrossedRowsTPC(60)
			fTrackFilterDaughters->SetMinNCrossedRowsTPC(60);
			fTrackFilterDaughters->SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);
			fTrackFilterDaughters->SetMaxChi2PerClusterTPC(4);
			fTrackFilterDaughters->SetMaxDCAToVertexZ(2);
			fTrackFilterDaughters->SetMaxChi2PerClusterITS(36);
		}
		if(fTrackID==2){ //! Higher: SetMinNCrossedRowsTPC(100)
			fTrackFilterDaughters->SetMinNCrossedRowsTPC(100);
			fTrackFilterDaughters->SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);
			fTrackFilterDaughters->SetMaxChi2PerClusterTPC(4);
			fTrackFilterDaughters->SetMaxDCAToVertexZ(2);
			fTrackFilterDaughters->SetMaxChi2PerClusterITS(36);
		}
		if(fTrackID==3){ //! Lower: SetMinRatioCrossedRowsOverFindableClustersTPC(0.7) 
			fTrackFilterDaughters->SetMinNCrossedRowsTPC(70);
			fTrackFilterDaughters->SetMinRatioCrossedRowsOverFindableClustersTPC(0.7);
			fTrackFilterDaughters->SetMaxChi2PerClusterTPC(4);
			fTrackFilterDaughters->SetMaxDCAToVertexZ(2);
			fTrackFilterDaughters->SetMaxChi2PerClusterITS(36);
		}
		if(fTrackID==4){ //! Higher: SetMinRatioCrossedRowsOverFindableClustersTPC(0.9)
			fTrackFilterDaughters->SetMinNCrossedRowsTPC(70);
			fTrackFilterDaughters->SetMinRatioCrossedRowsOverFindableClustersTPC(0.9);
			fTrackFilterDaughters->SetMaxChi2PerClusterTPC(4);
			fTrackFilterDaughters->SetMaxDCAToVertexZ(2);
			fTrackFilterDaughters->SetMaxChi2PerClusterITS(36);
		}
		if(fTrackID==5){ //! Lower: SetMaxChi2PerClusterTPC(3)  
			fTrackFilterDaughters->SetMinNCrossedRowsTPC(70);
			fTrackFilterDaughters->SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);
			fTrackFilterDaughters->SetMaxChi2PerClusterTPC(3);
			fTrackFilterDaughters->SetMaxDCAToVertexZ(2);
			fTrackFilterDaughters->SetMaxChi2PerClusterITS(36);
		}
		if(fTrackID==6){ //! Higher: SetMaxChi2PerClusterTPC(5)
			fTrackFilterDaughters->SetMinNCrossedRowsTPC(70);
			fTrackFilterDaughters->SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);
			fTrackFilterDaughters->SetMaxChi2PerClusterTPC(5);
			fTrackFilterDaughters->SetMaxDCAToVertexZ(2);
			fTrackFilterDaughters->SetMaxChi2PerClusterITS(36);
		}
		if(fTrackID==7){ //! Lower: SetMaxChi2PerClusterITS(25)
			fTrackFilterDaughters->SetMinNCrossedRowsTPC(70);
			fTrackFilterDaughters->SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);
			fTrackFilterDaughters->SetMaxChi2PerClusterTPC(4);
			fTrackFilterDaughters->SetMaxDCAToVertexZ(2);
			fTrackFilterDaughters->SetMaxChi2PerClusterITS(25);
		}
		if(fTrackID==8){ //! Higher: SetMaxChi2PerClusterITS(49)
			fTrackFilterDaughters->SetMinNCrossedRowsTPC(70);
			fTrackFilterDaughters->SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);
			fTrackFilterDaughters->SetMaxChi2PerClusterTPC(4);
			fTrackFilterDaughters->SetMaxDCAToVertexZ(2);
			fTrackFilterDaughters->SetMaxChi2PerClusterITS(49);
		}
		if(fTrackID==9){ //! Lower: SetMaxDCAToVertexZ(1)
			fTrackFilterDaughters->SetMinNCrossedRowsTPC(70);
			fTrackFilterDaughters->SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);
			fTrackFilterDaughters->SetMaxChi2PerClusterTPC(4);
			fTrackFilterDaughters->SetMaxDCAToVertexZ(1);
			fTrackFilterDaughters->SetMaxChi2PerClusterITS(36);
		}
		if(fTrackID==10){ //! Lower: SetMaxDCAToVertexZ(5)
			fTrackFilterDaughters->SetMinNCrossedRowsTPC(70);
			fTrackFilterDaughters->SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);
			fTrackFilterDaughters->SetMaxChi2PerClusterTPC(4);
			fTrackFilterDaughters->SetMaxDCAToVertexZ(5);
			fTrackFilterDaughters->SetMaxChi2PerClusterITS(36);
		}
		if((fTrackID==11)||(fTrackID==12)||(fTrackID==13)||(fTrackID==14)){ //! Nominal values
			printf("fTrackFilterDaughters for fTrackID = %d\n",fTrackID);
			fTrackFilterDaughters->SetMinNCrossedRowsTPC(70);
			fTrackFilterDaughters->SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);
			fTrackFilterDaughters->SetMaxChi2PerClusterTPC(4);
			fTrackFilterDaughters->SetMaxDCAToVertexZ(2);
			fTrackFilterDaughters->SetMaxChi2PerClusterITS(36);
		}

	}

	//OpenFile(1);
	fListOfObjects = new TList();
	fListOfObjects->SetOwner(kTRUE);

	//
	// Histograms
	//

	fEvents = new TH1F( "fEvents", "; Evt. Sel.",12,0,12);
	fEvents->GetXaxis()->SetBinLabel(1, "Processed");
	fEvents->GetXaxis()->SetBinLabel(2, "PhysSel+Trigger");
	fEvents->GetXaxis()->SetBinLabel(3, "INEL>0");
	fEvents->GetXaxis()->SetBinLabel(4, "BG");//NotinVertexcut");
	fEvents->GetXaxis()->SetBinLabel(5, "IsPileUpFromSPDinMultBins");//NotinVertexcut");
	fEvents->GetXaxis()->SetBinLabel(6, "Incom DAQ");//NotinVertexcut");
	fEvents->GetXaxis()->SetBinLabel(7, "Res&Proximity");//NotinVertexcut");
	fEvents->GetXaxis()->SetBinLabel(8, "|Vtz|<10cm");//NotinVertexcut");
	fListOfObjects->Add(fEvents);

	const int nPtBinsV0s = 21;
	double ptBinsV0s[nPtBinsV0s+1] = { 
		0.1, 0.2, 0.3, 0.4,
		0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.2, 1.4, 1.6, 
		1.8, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0, 7.0, 9.0};

	const int nPtBins = 45;
	double ptBins[nPtBins+1] = {
		0.0, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 
		0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 
		0.90, 0.95, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7,
		1.80, 1.90, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4,
		3.60, 3.80, 4.0, 4.5, 5.0, 5.5, 6.0, 7.0, 8.0, 10.0};

	const int nBinsNsigma = 100;
	double binsNsigma[nBinsNsigma+1] = {0};

	for(int i = 0; i <= nBinsNsigma; ++i){
		binsNsigma[i] = -10.0+i*0.2;
	}

	const int nDeltaPiBins   = 80;
	double DeltaPiBins[nDeltaPiBins+1] = { 0 };
	for(int i = 0; i <= nDeltaPiBins; ++i){
		DeltaPiBins[i] = 20.0+i*1.0;
	}

	const int ndEdxBins   = 200;
	double dEdxBins[ndEdxBins+1] = { 0 };
	for(int i = 0; i <= ndEdxBins; ++i){
		dEdxBins[i] = fdEdxLow+i*1.0;
	}

	const int nBetaBins   = 100;
	double BetaBins[nBetaBins+1] = { 0 };
	for(int i = 0; i <= nBetaBins; ++i){
		BetaBins[i] = 0.2+((double)i)/100.0;
	}

	const int nBinsRT = 50;
	double binsRT[nBinsRT+1] = {0};

	for(int i = 0; i <= nBinsRT; ++i)
		binsRT[i] = (((double)i)-0.5)/1.0;

	const char* ending[nHists] = {"02", "24", "46", "68"};

	fcutDCAxy = new TF1("fMaxDCAxy","[0]+[1]/(x^[2])",0,1e10);
	fcutDCAxy->SetParameter(0,0.0105);
	fcutDCAxy->SetParameter(1,0.0350);
	fcutDCAxy->SetParameter(2,1.1);

	fcutLow = new TF1("StandardPhiCutLow",  "0.1/x/x+TMath::Pi()/18.0-0.025", 0, 50);
	fcutHigh = new TF1("StandardPhiCutHigh", "0.12/x+TMath::Pi()/18.0+0.035", 0, 50);

	fEtaCalibrationPos = new TF1("fDeDxVsEtaPos", "pol7", 0.0, 1.0);
	fEtaCalibrationNeg = new TF1("fDeDxVsEtaNeg", "pol7", -1.0, 0.0);
	fEtaCalibrationPosEl = new TF1("fDeDxVsEtaPosEl", "pol4", 0.0, 1.0);
	fEtaCalibrationNegEl = new TF1("fDeDxVsEtaNegEl", "pol4", -1.0, 0.0);

	hNchTSData = new TH1F("hMultTSData",";#it{N}_{acc}^{TS}; Entries",nBinsRT,binsRT);
	fListOfObjects->Add(hNchTSData);

	hPhiTotal = new TH2F("hPhiSum","; #eta; #varphi",50,-0.8,0.8,100,-TMath::Pi()/2.0,5.0*TMath::Pi()/2.0);
	fListOfObjects->Add(hPhiTotal);

	hPhiStandard = new TH2F("hPhiStandard","; #eta; #varphi",50,-0.8,0.8,100,-TMath::Pi()/2.0,5.0*TMath::Pi()/2.0);
	fListOfObjects->Add(hPhiStandard);

	hPhiHybrid1 = new TH2F("hPhiHybrid1","; #eta; #varphi",50,-0.8,0.8,100,-TMath::Pi()/2.0,5.0*TMath::Pi()/2.0);
	fListOfObjects->Add(hPhiHybrid1);

	hPhiHybrid2 = new TH2F("hPhiHybrid2","; #eta; #varphi",50,-0.8,0.8,100,-TMath::Pi()/2.0,5.0*TMath::Pi()/2.0);
	fListOfObjects->Add(hPhiHybrid2);

	hPhiLeading = new TH1F("hPhiLeading","; #varphi",100,-TMath::Pi()/2.0,5.0*TMath::Pi()/2.0);
	fListOfObjects->Add(hPhiLeading);

	hDeltaPhiDeltaEta = new TH2F("hDeltaPhiDeltaEta","; #Delta #varphi; #Delta #eta",100,-TMath::Pi(),5.0*TMath::Pi()/2.0,100,-2.0,2.0); 
	fListOfObjects->Add(hDeltaPhiDeltaEta);

	// Histos rTPC

	hMIPVsEta = new TH2F("hMIPVsEta","; #eta; dE/dx_{MIP, primary tracks}",50,-0.8,0.8,fDeDxMIPMax-fDeDxMIPMin,fDeDxMIPMin,fDeDxMIPMax);
	pMIPVsEta = new TProfile("pMIPVsEta","; #eta; #LT dE/dx #GT_{MIP, primary tracks}",50,-0.8,0.8,fDeDxMIPMin,fDeDxMIPMax);

	hPlateauVsEta = new TH2F("hPlateauVsEta","; #eta; dE/dx_{Plateau, primary tracks}",50,-0.8,0.8,50, 60, 110);
	pPlateauVsEta = new TProfile("pPlateauVsEta","; #eta; #LT dE/dx #GT_{Plateau, primary tracks}",50,-0.8,0.8, 60, 110);

	hMIPVsEtaV0s = new TH2F("hMIPVsEtaV0s","; #eta; dE/dx_{MIP, primary tracks}",50,-0.8,0.8,fDeDxMIPMax-fDeDxMIPMin,fDeDxMIPMin,fDeDxMIPMax);
	pMIPVsEtaV0s = new TProfile("pMIPVsEtaV0s","; #eta; #LT dE/dx #GT_{MIP, primary tracks}",50,-0.8,0.8,fDeDxMIPMin,fDeDxMIPMax);

	hPhirTPC = new TH2F("hPhirTPC", ";pt; #phi'", nPtBinsV0s, ptBinsV0s, 90, -0.05, 0.4);

	fListOfObjects->Add(hMIPVsEta);
	fListOfObjects->Add(pMIPVsEta);
	fListOfObjects->Add(hPlateauVsEta);
	fListOfObjects->Add(pPlateauVsEta);
	fListOfObjects->Add(hMIPVsEtaV0s);
	fListOfObjects->Add(pMIPVsEtaV0s);
	fListOfObjects->Add(hPhirTPC);

	for( int j = 0; j < nHists; j++ ){

		histEV0[j] = new TH2F(Form("histEV0_%s",ending[j]),"Electronsons from V0s; #it{p} (GeV/#it{c}); d#it{e}d#it{x}",nPtBinsV0s,ptBinsV0s,nDeltaPiBins,DeltaPiBins);
		histPV0[j] = new TH2F(Form("histPV0_%s",ending[j]),"Electronsons from V0s; #it{p} (GeV/#it{c}); d#it{e}d#it{x}",nPtBinsV0s,ptBinsV0s,nDeltaPiBins,DeltaPiBins);
		histPiV0[j] = new TH2F(Form("histPiV0_%s",ending[j]),"Electronsons from V0s; #it{p} (GeV/#it{c}); d#it{e}d#it{x}",nPtBinsV0s,ptBinsV0s,nDeltaPiBins,DeltaPiBins);
		histPiTof[j] = new TH2F(Form("hPiTOF_%s",ending[j]),"Pions from TOF;#it{p} (GeV/#it{c});d#it{e}d#it{x}",nPtBinsV0s,ptBinsV0s,nDeltaPiBins,DeltaPiBins);

		hMIPVsPhi[j] = new TH2F(Form("hMIPVsPhi_%s",ending[j]),";#phi (rad); dE/dx MIP",100,-TMath::Pi()/2.0,5.0*TMath::Pi()/2.0,fDeDxMIPMax-fDeDxMIPMin,fDeDxMIPMin,fDeDxMIPMax);

		pMIPVsPhi[j] = new TProfile(Form("pMIPVsPhi_%s",ending[j]),";#phi (rad); dE/dx MIP",100,-TMath::Pi()/2.0,5.0*TMath::Pi()/2.0,fDeDxMIPMin, fDeDxMIPMax);

		hPlateauVsPhi[j] = new TH2F(Form("hPlateauVsPhi_%s",ending[j]),";#phi (rad); dE/dx Plateau",100,-TMath::Pi()/2.0,5.0*TMath::Pi()/2.0,20, 70, 90);

		pPlateauVsPhi[j] = new TProfile(Form("pPlateauVsPhi_%s",ending[j]),";#phi (rad); dE/dx Plateau",100,-TMath::Pi()/2.0,5.0*TMath::Pi()/2.0,fDeDxMIPMax, 95);
		pPlateauVsPhi[j]->Sumw2();

		if(!fAnalysisMC){
			fListOfObjects->Add(histEV0[j]);
			fListOfObjects->Add(histPV0[j]);
			fListOfObjects->Add(histPiV0[j]);
			fListOfObjects->Add(histPiTof[j]);
			fListOfObjects->Add(hMIPVsPhi[j]);
			fListOfObjects->Add(pMIPVsPhi[j]);
			fListOfObjects->Add(hPlateauVsPhi[j]);
			fListOfObjects->Add(pPlateauVsPhi[j]);
		}
	}	// Only ending


	for( int j = 0; j < nHists; j++ ){

		hPtVsP[j] = new TH2F(Form("hPtVsP_%s",ending[j]),";#it{p} (GeV/#it{c}); #it{p}_{T} (GeV/#it{c})",nPtBins,ptBins,nPtBins,ptBins);
		fListOfObjects->Add(hPtVsP[j]);

		hnSigmaElectrons[j] = new TH2F(Form("hnSigmaElectrons_%s",ending[j]),";#it{p}_{T}^{rec};n#sigma;#it{N}_{acc}",nPtBins,ptBins,nBinsNsigma,binsNsigma);
		fListOfObjects->Add(hnSigmaElectrons[j]);


	}	// Only ending


	for(int r = 0; r < nRegion; ++r){

		if(r<(nRegion-1)){

			hPhiData[r] = new TH1F(Form("hPhiData_%s",Region[r]),"",64,-TMath::Pi()/2.0,3.0*TMath::Pi()/2.0);


		}

		for( int j = 0; j < nHists; j++ ){

			if(r<(nRegion-1)){

				hNchVsPtPosTPC[r][j] = new TH2F(Form("hNchVsPtPosTPC_%s_%s",Region[r],ending[j]),";#it{p}_{T} (GeV/#it{c}); #it{N}_{acc}^{TS}",nPtBins,ptBins,nBinsRT,binsRT);

				hNchVsPtNegTPC[r][j] = new TH2F(Form("hNchVsPtNegTPC_%s_%s",Region[r],ending[j]),";#it{p}_{T} (GeV/#it{c}); #it{N}_{acc}^{TS}",nPtBins,ptBins,nBinsRT,binsRT);

				hNchVsPtPosTOF[r][j] = new TH2F(Form("hNchVsPtPosTOF_%s_%s",Region[r],ending[j]),";#it{p}_{T} (GeV/#it{c}); #it{N}_{acc}^{TS}",nPtBins,ptBins,nBinsRT,binsRT);

				hNchVsPtNegTOF[r][j] = new TH2F(Form("hNchVsPtNegTOF_%s_%s",Region[r],ending[j]),";#it{p}_{T} (GeV/#it{c}); #it{N}_{acc}^{TS}",nPtBins,ptBins,nBinsRT,binsRT);

				hNchVsPPosTOF[r][j] = new TH2F(Form("hNchVsPPosTOF_%s_%s",Region[r],ending[j]),";#it{p} (GeV/#it{c}); #it{N}_{acc}^{TS}",nPtBins,ptBins,nBinsRT,binsRT);

				hNchVsPNegTOF[r][j] = new TH2F(Form("hNchVsPNegTOF_%s_%s",Region[r],ending[j]),";#it{p} (GeV/#it{c}); #it{N}_{acc}^{TS}",nPtBins,ptBins,nBinsRT,binsRT);

				hNchVsPrTPC[r][j] = new TH2F(Form("hNchVsPrTPC_%s_%s",Region[r],ending[j]),";#it{p} (GeV/#it{c}); #it{N}_{acc}^{TS}",nPtBins,ptBins,nBinsRT,binsRT);
				hNchVsPtrTPC[r][j] = new TH2F(Form("hNchVsPtrTPC_%s_%s",Region[r],ending[j]),";#it{p}_{T} (GeV/#it{c}); #it{N}_{acc}^{TS}",nPtBins,ptBins,nBinsRT,binsRT);


				hNchVsPtDataPosPionTPC[r][j] = new TH3F(Form("hNchVsPtData_Pos_Pion_TPC_%s_%s",Region[r],ending[j]),";#it{p}_{T}^{rec};n#sigma;#it{N}_{acc}",nPtBins,ptBins,nBinsNsigma,binsNsigma,nBinsRT,binsRT);

				hNchVsPtDataNegPionTPC[r][j] = new TH3F(Form("hNchVsPtData_Neg_Pion_TPC_%s_%s",Region[r],ending[j]),";#it{p}_{T}^{rec};n#sigma;#it{N}_{acc}",nPtBins,ptBins,nBinsNsigma,binsNsigma,nBinsRT,binsRT);

				hNchVsPtDataPosKaonTPC[r][j] = new TH3F(Form("hNchVsPtData_Pos_Kaon_TPC_%s_%s",Region[r],ending[j]),";#it{p}_{T}^{rec};n#sigma;#it{N}_{acc}",nPtBins,ptBins,nBinsNsigma,binsNsigma,nBinsRT,binsRT);

				hNchVsPtDataNegKaonTPC[r][j] = new TH3F(Form("hNchVsPtData_Neg_Kaon_TPC_%s_%s",Region[r],ending[j]),";#it{p}_{T}^{rec};n#sigma;#it{N}_{acc}",nPtBins,ptBins,nBinsNsigma,binsNsigma,nBinsRT,binsRT);

				hNchVsPtDataPosProtonTPC[r][j] = new TH3F(Form("hNchVsPtData_Pos_Proton_TPC_%s_%s",Region[r],ending[j]),";#it{p}_{T}^{rec};n#sigma;#it{N}_{acc}",nPtBins,ptBins,nBinsNsigma,binsNsigma,nBinsRT,binsRT);

				hNchVsPtDataNegProtonTPC[r][j] = new TH3F(Form("hNchVsPtData_Neg_Proton_TPC_%s_%s",Region[r],ending[j]),";#it{p}_{T}^{rec};n#sigma;#it{N}_{acc}",nPtBins,ptBins,nBinsNsigma,binsNsigma,nBinsRT,binsRT);

				hNchVsPtDataPosTOF[r][j] = new TH3F(Form("hNchVsPtData_Pos_TOF_%s_%s",Region[r],ending[j]),";#it{p}^{rec};#beta;#it{N}_{acc}",nPtBins,ptBins,nBetaBins,BetaBins,nBinsRT,binsRT);

				hNchVsPtDataNegTOF[r][j] = new TH3F(Form("hNchVsPtData_Neg_TOF_%s_%s",Region[r],ending[j]),";#it{p}^{rec};#beta;#it{N}_{acc}",nPtBins,ptBins,nBetaBins,BetaBins,nBinsRT,binsRT);

				hDeDxVsP[r][j] = new TH3F(Form("hDeDxVsP_%s_%s",Region[r],ending[j]),";#it{p} (GeV/#it{c});#it{d}E/#it{d}x;#it{N}_{acc}^{TS}",nPtBins,ptBins,ndEdxBins,dEdxBins,nBinsRT,binsRT);

			}

			if(!fAnalysisMC){

				if(r<(nRegion-1)){

					fListOfObjects->Add(hNchVsPtPosTPC[r][j]);
					fListOfObjects->Add(hNchVsPtNegTPC[r][j]);
					fListOfObjects->Add(hNchVsPtPosTOF[r][j]);
					fListOfObjects->Add(hNchVsPtNegTOF[r][j]);
					fListOfObjects->Add(hNchVsPPosTOF[r][j]);
					fListOfObjects->Add(hNchVsPNegTOF[r][j]);
					fListOfObjects->Add(hNchVsPrTPC[r][j]);
					fListOfObjects->Add(hNchVsPtrTPC[r][j]);

					fListOfObjects->Add(hNchVsPtDataPosPionTPC[r][j]);
					fListOfObjects->Add(hNchVsPtDataNegPionTPC[r][j]);
					fListOfObjects->Add(hNchVsPtDataPosKaonTPC[r][j]);
					fListOfObjects->Add(hNchVsPtDataNegKaonTPC[r][j]);
					fListOfObjects->Add(hNchVsPtDataPosProtonTPC[r][j]);
					fListOfObjects->Add(hNchVsPtDataNegProtonTPC[r][j]);

					fListOfObjects->Add(hNchVsPtDataPosTOF[r][j]);
					fListOfObjects->Add(hNchVsPtDataNegTOF[r][j]);

					fListOfObjects->Add(hDeDxVsP[r][j]);

				}

			}

		}	// ending 

		if(!fAnalysisMC){
			if(r<(nRegion-1)){
				fListOfObjects->Add(hPhiData[r]);	
			}
		}

	}

	////	fEventCuts.AddQAplotsToList(fListOfObjects);
	PostData(1, fListOfObjects);

}
//______________________________________________________________________________
void AliAnalysisTaskSpectraRT::UserExec(Option_t *)
{
	// Main loop

	//
	// First we make sure that we have valid input(s)!
	//

	AliVEvent *event = InputEvent();
	if (!event) {
		Error("UserExec", "Could not retrieve event");
		return;
	}

	fESD = dynamic_cast<AliESDEvent*>(event);
	if(!fESD){
		Printf("%s:%d ESDEvent not found in Input Manager",(char*)__FILE__,__LINE__);
		this->Dump();
		return;
	}

	if (fAnalysisMC){

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

	}

	utils = new AliAnalysisUtils();
	if (!utils)
	{
		cout<<"------- No AnalysisUtils Object Found --------"<<utils<<endl;
		return;
	}

	fEvents->Fill(0.5);

	AliHeader* headerMC = fMC->Header();
	bool isGoodVtxPosMC = kFALSE;
	if (fAnalysisMC){
		AliGenEventHeader* genHeader = headerMC->GenEventHeader();
		TArrayF vtxMC(3); // primary vertex  MC
		vtxMC[0]=9999; vtxMC[1]=9999;  vtxMC[2]=9999; //initialize with dummy
		if (genHeader) {
			genHeader->PrimaryVertex(vtxMC);
		}

		if(TMath::Abs(vtxMC[2])<=10){
			isGoodVtxPosMC = kTRUE;
		}

		// Before trigger selection
		GetLeadingObject(kTRUE);// leading particle at gen level
	}

	UInt_t fSelectMask= fInputHandler->IsEventSelected();
	bool isINT7selected = fSelectMask&AliVEvent::kINT7;
	if(!isINT7selected)
		return;
	fEvents->Fill(1.5);

	int INEL = -1;
	INEL = AliESDtrackCuts::GetReferenceMultiplicity(fESD, AliESDtrackCuts::kTracklets, 1.0);
	if( INEL < 1 )
		return;
	fEvents->Fill(2.5);

	if( utils->IsSPDClusterVsTrackletBG(fESD) )
		return;
	fEvents->Fill(3.5);

	if( fESD->IsPileupFromSPDInMultBins() )
		return;
	fEvents->Fill(4.5);

	if( fESD->IsIncompleteDAQ())
		return;
	fEvents->Fill(5.5);

	if( !selectVertex2015pp(fESD,kTRUE,kFALSE,kTRUE) )
		return;
	fEvents->Fill(6.5);

	if( !IsGoodZvertexPos(fESD) )
		return;
	fEvents->Fill(7.5);

	GetLeadingObject(kFALSE);


	if(fIsMCclosure){
		double randomUE = gRandom->Uniform(0.0,1.0);
		if(randomUE<0.5){// corrections (50% stat.)
			if(isGoodVtxPosMC){
				if( ( fGenLeadPt>=fLeadPtCutMin && fGenLeadPt<fLeadPtCutMax ) && ( fRecLeadPt>=fLeadPtCutMin && fRecLeadPt<fLeadPtCutMax )){
					GetDetectorResponse();
					GetMCCorrections();
				}
			}
		}
		else{// for testing the method
			if( ( fGenLeadPt>=fLeadPtCutMin && fGenLeadPt<fLeadPtCutMax ) && ( fRecLeadPt>=fLeadPtCutMin && fRecLeadPt<fLeadPtCutMax )){
				GetMultiplicityDistributions();
			}
		}
	}
	else{
		if(fAnalysisMC){
			if(isGoodVtxPosMC){
				if( (fGenLeadPt>=fLeadPtCutMin && fGenLeadPt<fLeadPtCutMax)&&(fRecLeadPt>=fLeadPtCutMin && fRecLeadPt<fLeadPtCutMax)){
					GetMultiplicityDistributions();
					GetDetectorResponse();
					GetMCCorrections();
				}
			}
		}
		else{
			if(( fRecLeadPt>=fLeadPtCutMin && fRecLeadPt<fLeadPtCutMax )){
				ProduceArrayTrksESD();
			}

			ProduceArrayV0ESD();

		}
	}

	PostData(1, fListOfObjects);
}
//_____________________________________________________________________________
void AliAnalysisTaskSpectraRT::GetLeadingObject(bool isMC) {

	Double_t flPt = 0.0;// leading pT
	Double_t flEta = 0.0;
	Double_t flPhi = 0.0;
	Int_t flIndex = 0;

	if(isMC){
		for (Int_t i = 0; i < fMC->GetNumberOfTracks(); i++) {

			AliMCParticle* particle = (AliMCParticle*)fMC->GetTrack(i);
			if (!particle) continue;

			if (!fMC->IsPhysicalPrimary(i)) continue;  
			if (particle->Charge() == 0) continue;
			if ( TMath::Abs(particle->Eta()) > fEtaCut )continue;
			if( particle->Pt() < fPtMin)continue;

			if (flPt<particle->Pt()){
				flPt = particle->Pt();
				flPhi = particle->Phi();
				flIndex = i;
			}
		}

		fGenLeadPhi = flPhi;
		fGenLeadPt  = flPt;
		fGenLeadIn  = flIndex;
	}
	else{

		int iTracks(fESD->GetNumberOfTracks());          
		for(int i=0; i < iTracks; i++) {                

			AliESDtrack* track = static_cast<AliESDtrack*>(fESD->GetTrack(i)); 
			if(!track) continue;
			if(TMath::Abs(track->Eta()) > fEtaCut) continue;
			if(track->Pt() < fPtMin) continue;

			AliESDtrack* track_hybrid = 0x0;
			if(!fSelectHybridTracks){
				if(!fTrackFilter->AcceptTrack(track)) { continue; } 
				else{ track_hybrid = new AliESDtrack(*track); }
			}else{
				track_hybrid = SetHybridTrackCuts(track,kFALSE,kFALSE,kFALSE);
				if(!track_hybrid) { continue; }
			}

			if(!fGeometricalCut->AcceptTrack(track_hybrid)) continue;

			if (flPt<track_hybrid->Pt()){
				flPt  = track_hybrid->Pt();
				flEta = track_hybrid->Eta();
				flPhi = track_hybrid->Phi();
				flIndex = i;
			}

			delete track_hybrid;

		}

		fRecLeadEta = flEta;
		fRecLeadPhi = flPhi;
		fRecLeadPt  = flPt;
		fRecLeadIn  = flIndex;

	}

}
//_____________________________________________________________________________
void AliAnalysisTaskSpectraRT::GetMultiplicityDistributions(){

	int multTSgen=0;
	int multTSrec=0;
	const double pi = TMath::Pi();

	for (int i = 0; i < fMC->GetNumberOfTracks(); i++) {

		AliMCParticle* particle = (AliMCParticle*)fMC->GetTrack(i);
		if (!particle) continue;
		if (i==fGenLeadIn) continue;
		if (!fMC->IsPhysicalPrimary(i)) continue;
		if (particle->Charge() == 0) continue;
		if (TMath::Abs(particle->Eta()) > fEtaCut)continue;
		if (particle->Pt() < fPtMin)continue;

		double DPhi = DeltaPhi(particle->Phi(), fGenLeadPhi);

		if(TMath::Abs(DPhi)<pi/3.0){
			continue;
		}
		else if(TMath::Abs(DPhi-pi)<pi/3.0){
			continue;
		}
		else{
			multTSgen++;
		}
	}

	int iTracks = 0;           
	iTracks = fESD->GetNumberOfTracks();           

	for(int i = 0; i < iTracks; i++) {                 

		if(i==fRecLeadIn) continue;
		AliESDtrack* esdtrack = static_cast<AliESDtrack*>(fESD->GetTrack(i));  
		if(!esdtrack) continue;
		if(esdtrack->Charge() == 0 ) continue;
		if(TMath::Abs(esdtrack->Eta()) > fEtaCut) continue;
		if(esdtrack->Pt() < fPtMin) continue;

		AliESDtrack* track = 0x0;
		if(!fSelectHybridTracks){
			if(!fTrackFilter->AcceptTrack(esdtrack)) { continue; } 
			else{ track = esdtrack; }
		}else{
			track = SetHybridTrackCuts(esdtrack,kFALSE,kFALSE,kFALSE);
			if(!track) { continue; }
		}

		hPhiTotal->Fill(track->Eta(),track->Phi());
		double DPhi = DeltaPhi(track->Phi(), fRecLeadPhi);

		if(TMath::Abs(DPhi)<pi/3.0){
			continue;
		}
		else if(TMath::Abs(DPhi-pi)<pi/3.0){
			continue;
		}
		else{
			multTSrec++;
		}

		delete track;

	}
}
//_____________________________________________________________________________
void AliAnalysisTaskSpectraRT::GetDetectorResponse() {

	int multTSgen = 0;
	int multTSrec = 0;
	double pi = TMath::Pi();

	for (int i = 0; i < fMC->GetNumberOfTracks(); i++) {

		if(i==fGenLeadIn) continue;
		AliMCParticle* particle = (AliMCParticle*)fMC->GetTrack(i);
		if (!particle) continue;
		if (!fMC->IsPhysicalPrimary(i)) continue;
		if (particle->Charge() == 0) continue;
		if (TMath::Abs(particle->Eta()) > fEtaCut)continue;
		if (particle->Pt() < fPtMin)continue;

		double DPhi = DeltaPhi(particle->Phi(), fGenLeadPhi);

		if(TMath::Abs(DPhi)<pi/3.0){
			continue;		
		}
		else if(TMath::Abs(DPhi-pi)<pi/3.0){
			continue;		
		}
		else{
			multTSgen++;
		}
	}

	int iTracks = 0;          
	iTracks = fESD->GetNumberOfTracks();          

	for(int i = 0; i < iTracks; i++){              

		if(i==fRecLeadIn) continue;
		AliESDtrack* esdtrack = static_cast<AliESDtrack*>(fESD->GetTrack(i)); 
		if(!esdtrack) continue;
		if(esdtrack->Charge() == 0 ) continue;
		if(TMath::Abs(esdtrack->Eta()) > fEtaCut) continue;
		if(esdtrack->Pt() < fPtMin) continue;

		AliESDtrack* track = 0x0;
		if(!fSelectHybridTracks){
			if(!fTrackFilter->AcceptTrack(esdtrack)) { continue; } 
			else{ track = esdtrack; }
		}else{
			track = SetHybridTrackCuts(esdtrack,kFALSE,kFALSE,kFALSE);
			if(!track) { continue; }
		}

		double DPhi = DeltaPhi(track->Phi(), fRecLeadPhi);

		if(TMath::Abs(DPhi)<pi/3.0){
			continue;		
		}
		else if(TMath::Abs(DPhi-pi)<pi/3.0){
			continue;		
		}
		else{
			multTSrec++;
		}

		delete track;
	}

}
//_____________________________________________________________________________
void AliAnalysisTaskSpectraRT::GetMCCorrections(){

	int multTSgen = 0;
	int multTSrec = 0;
	double pi = TMath::Pi();

	for (int i = 0; i < fMC->GetNumberOfTracks(); i++){

		if (i==fGenLeadIn) continue;
		AliMCParticle* particle = (AliMCParticle*)fMC->GetTrack(i);
		if (!particle) continue;
		if (!fMC->IsPhysicalPrimary(i)) continue;
		if (particle->Charge() == 0) continue;
		if (TMath::Abs(particle->Eta()) > fEtaCut)continue;
		if (particle->Pt() < fPtMin) continue;

		double DPhi = DeltaPhi(particle->Phi(), fGenLeadPhi);

		if(TMath::Abs(DPhi)<pi/3.0){
			continue;
		}
		else if(TMath::Abs(DPhi-pi)<pi/3.0){
			continue;
		}
		else{
			multTSgen++;
		}
	}

	int iTracks(fESD->GetNumberOfTracks());          
	for(int i = 0; i < iTracks; i++){              

		if(i==fRecLeadIn) continue;
		AliESDtrack* esdtrack = static_cast<AliESDtrack*>(fESD->GetTrack(i)); 
		if(!esdtrack) continue;
		if(TMath::Abs(esdtrack->Eta()) > fEtaCut) continue;
		if(esdtrack->Pt() < fPtMin) continue;
		if(esdtrack->Charge()==0 ) continue;

		AliESDtrack* track = 0x0;
		if(!fSelectHybridTracks){
			if(!fTrackFilter->AcceptTrack(esdtrack)) { continue; } 
			else{ track = esdtrack; }
		}else{
			track = SetHybridTrackCuts(esdtrack,kTRUE,kTRUE,kTRUE);
			if(!track) { continue; }
		}

		double DPhi = DeltaPhi(track->Phi(), fRecLeadPhi);
		if(TMath::Abs(DPhi)<pi/3.0){
			continue;
		}
		else if(TMath::Abs(DPhi-pi)<pi/3.0){
			continue;
		}
		else{
			multTSrec++;
		}

		delete track;

	}
}
//_____________________________________________________________________________
double AliAnalysisTaskSpectraRT::DeltaPhi(Double_t phi, Double_t Lphi,
		Double_t rangeMin, Double_t rangeMax)
{

	double dphi = -999;
	double pi = TMath::Pi();
	//	if(Lphi > 2*pi || Lphi < 0)cout << "Lphi :: " << Lphi << endl;
	//	if(phi  > 2*pi || phi < 0)cout << "phi = " << phi << endl;

	if(phi < 0)          phi += 2*pi;
	else if(phi > 2*pi)  phi -= 2*pi;
	if(Lphi < 0)         Lphi += 2*pi;
	else if(Lphi > 2*pi) Lphi -= 2*pi;
	dphi = Lphi - phi;
	if (dphi < rangeMin)      dphi += 2*pi;
	else if (dphi > rangeMax) dphi -= 2*pi;

	return dphi;
}
//_____________________________________________________________________________
double AliAnalysisTaskSpectraRT::DeltaEta(Double_t eta, Double_t Leta)
{

	double deta = -999.0;
	deta = Leta - eta;
	return deta;
}
//_____________________________________________________________________________
short AliAnalysisTaskSpectraRT::GetPidCode(Int_t pdgCode) const
{
	// return our internal code for pions, kaons, and protons

	short pidCode = 6;

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
void AliAnalysisTaskSpectraRT::ProduceArrayTrksESD(){

	const double pi = TMath::Pi();
	int multTSdata = 0;

	int iTracks(fESD->GetNumberOfTracks());          
	for(int i = 0; i < iTracks; i++){              

		if(i==fRecLeadIn) continue;
		AliESDtrack* esdtrack = static_cast<AliESDtrack*>(fESD->GetTrack(i)); 
		if(!esdtrack) continue;

		AliESDtrack* track = 0x0;
		if(!fSelectHybridTracks){
			if(!fTrackFilter->AcceptTrack(esdtrack)) { continue; } 
			else{ track = new AliESDtrack(*esdtrack); }
		}else{
			track = SetHybridTrackCuts(esdtrack,kTRUE,kTRUE,kTRUE);
			if(!track) { continue; }
		}

		if(TMath::Abs(track->Eta()) > fEtaCut) continue;
		if(track->Pt() < fPtMin) continue;

		hPhiTotal->Fill(track->Eta(),track->Phi());
		double DPhi = DeltaPhi(track->Phi(), fRecLeadPhi);

		if(TMath::Abs(DPhi)<pi/3.0){
			hPhiData[0]->Fill(DPhi);
		}
		else if(TMath::Abs(DPhi-pi)<pi/3.0){
			hPhiData[1]->Fill(DPhi);
		}
		else{
			hPhiData[2]->Fill(DPhi);
			multTSdata++;
		}

		delete track;

	}

	hNchTSData->Fill(multTSdata);
	hPhiLeading->Fill(fRecLeadPhi);

	for(int iT = 0; iT < iTracks; iT++) {

		if(iT==fRecLeadIn) continue;
		AliESDtrack* track = (AliESDtrack*)fESD->GetTrack(iT);
		if(!track) continue;

		AliESDtrack* esdTrack = 0x0;
		if(!fSelectHybridTracks){
			if(!fTrackFilter->AcceptTrack(track)) { continue; } 
			else{ esdTrack = new AliESDtrack(*track); }
		}else{
			esdTrack = SetHybridTrackCuts(track,kFALSE,kFALSE,kFALSE);
			if(!esdTrack) { continue; }
		}

		if(TMath::Abs(esdTrack->Eta()) > fEtaCut) continue;
		if(esdTrack->GetTPCsignalN() < fNcl) continue;
		if(esdTrack->Pt() < fPtMin) continue;

		double DPhi = DeltaPhi(esdTrack->Phi(), fRecLeadPhi);
		double DEta = DeltaEta(esdTrack->Eta(), fRecLeadEta);
		hDeltaPhiDeltaEta->Fill(DPhi,DEta);

		int nh = -1;
		double eta = esdTrack->Eta();
		if(TMath::Abs(eta)<0.2)
			nh = 0;
		else if(TMath::Abs(eta)>=0.2 && TMath::Abs(eta)<0.4)
			nh = 1;
		else if(TMath::Abs(eta)>=0.4 && TMath::Abs(eta)<0.6)
			nh = 2;
		else
			nh = 3;

		if( nh < 0 )
			continue;

		if(TMath::Abs(DPhi)<pi/3.0){
			if(esdTrack->Charge() > 0){
				hNchVsPtPosTPC[0][nh]->Fill(esdTrack->Pt(),multTSdata);
				hNchVsPtDataPosPionTPC[0][nh]->Fill(esdTrack->Pt(),fPIDResponse->NumberOfSigmasTPC(esdTrack,AliPID::kPion),multTSdata);	
				hNchVsPtDataPosKaonTPC[0][nh]->Fill(esdTrack->Pt(),fPIDResponse->NumberOfSigmasTPC(esdTrack,AliPID::kKaon),multTSdata);	
				hNchVsPtDataPosProtonTPC[0][nh]->Fill(esdTrack->Pt(),fPIDResponse->NumberOfSigmasTPC(esdTrack,AliPID::kProton),multTSdata);	
			}
			if(esdTrack->Charge() < 0){
				hNchVsPtNegTPC[0][nh]->Fill(esdTrack->Pt(),multTSdata);
				hNchVsPtDataNegPionTPC[0][nh]->Fill(esdTrack->Pt(),fPIDResponse->NumberOfSigmasTPC(esdTrack,AliPID::kPion),multTSdata);	
				hNchVsPtDataNegKaonTPC[0][nh]->Fill(esdTrack->Pt(),fPIDResponse->NumberOfSigmasTPC(esdTrack,AliPID::kKaon),multTSdata);	
				hNchVsPtDataNegProtonTPC[0][nh]->Fill(esdTrack->Pt(),fPIDResponse->NumberOfSigmasTPC(esdTrack,AliPID::kProton),multTSdata);	
			}
		}
		else if(TMath::Abs(DPhi-pi)<pi/3.0){
			if(esdTrack->Charge() > 0){
				hNchVsPtPosTPC[1][nh]->Fill(esdTrack->Pt(),multTSdata);
				hNchVsPtDataPosPionTPC[1][nh]->Fill(esdTrack->Pt(),fPIDResponse->NumberOfSigmasTPC(esdTrack,AliPID::kPion),multTSdata);	
				hNchVsPtDataPosKaonTPC[1][nh]->Fill(esdTrack->Pt(),fPIDResponse->NumberOfSigmasTPC(esdTrack,AliPID::kKaon),multTSdata);	
				hNchVsPtDataPosProtonTPC[1][nh]->Fill(esdTrack->Pt(),fPIDResponse->NumberOfSigmasTPC(esdTrack,AliPID::kProton),multTSdata);	
			}
			if(esdTrack->Charge() < 0){
				hNchVsPtNegTPC[1][nh]->Fill(esdTrack->Pt(),multTSdata);
				hNchVsPtDataNegPionTPC[1][nh]->Fill(esdTrack->Pt(),fPIDResponse->NumberOfSigmasTPC(esdTrack,AliPID::kPion),multTSdata);	
				hNchVsPtDataNegKaonTPC[1][nh]->Fill(esdTrack->Pt(),fPIDResponse->NumberOfSigmasTPC(esdTrack,AliPID::kKaon),multTSdata);	
				hNchVsPtDataNegProtonTPC[1][nh]->Fill(esdTrack->Pt(),fPIDResponse->NumberOfSigmasTPC(esdTrack,AliPID::kProton),multTSdata);	
			}
		}
		else{
			if(esdTrack->Charge() > 0){
				hNchVsPtPosTPC[2][nh]->Fill(esdTrack->Pt(),multTSdata);
				hNchVsPtDataPosPionTPC[2][nh]->Fill(esdTrack->Pt(),fPIDResponse->NumberOfSigmasTPC(esdTrack,AliPID::kPion),multTSdata);	
				hNchVsPtDataPosKaonTPC[2][nh]->Fill(esdTrack->Pt(),fPIDResponse->NumberOfSigmasTPC(esdTrack,AliPID::kKaon),multTSdata);	
				hNchVsPtDataPosProtonTPC[2][nh]->Fill(esdTrack->Pt(),fPIDResponse->NumberOfSigmasTPC(esdTrack,AliPID::kProton),multTSdata);	
			}
			if(esdTrack->Charge() < 0){
				hNchVsPtNegTPC[2][nh]->Fill(esdTrack->Pt(),multTSdata);
				hNchVsPtDataNegPionTPC[2][nh]->Fill(esdTrack->Pt(),fPIDResponse->NumberOfSigmasTPC(esdTrack,AliPID::kPion),multTSdata);	
				hNchVsPtDataNegKaonTPC[2][nh]->Fill(esdTrack->Pt(),fPIDResponse->NumberOfSigmasTPC(esdTrack,AliPID::kKaon),multTSdata);	
				hNchVsPtDataNegProtonTPC[2][nh]->Fill(esdTrack->Pt(),fPIDResponse->NumberOfSigmasTPC(esdTrack,AliPID::kProton),multTSdata);	
			}
		}

		hnSigmaElectrons[nh]->Fill(esdTrack->Pt(),fPIDResponse->NumberOfSigmasTPC(esdTrack,AliPID::kElectron));

		hPtVsP[nh]->Fill(esdTrack->P(),esdTrack->Pt());

		//
		//_______________________________ TOF PID
		//

		bool IsTOFout = kFALSE;
		IsTOFout = TOFPID(esdTrack);
		if(IsTOFout){

			double trkLength = esdTrack->GetIntegratedLength();
			double beta = trkLength/((esdTrack->GetTOFsignal()-fPIDResponse->GetTOFResponse().GetStartTime(esdTrack->P()))*C_Value);
			double inttime[5]={0,0,0,0,0};
			esdTrack->GetIntegratedTimes(inttime);// Returns the array with integrated times for each particle hypothesis

			if(TMath::Abs(DPhi)<pi/3.0){
				if(esdTrack->Charge() > 0){
					hNchVsPPosTOF[0][nh]->Fill(esdTrack->P(),multTSdata);
					hNchVsPtPosTOF[0][nh]->Fill(esdTrack->Pt(),multTSdata);
					hNchVsPtDataPosTOF[0][nh]->Fill(esdTrack->P(),beta,multTSdata);
				}
				if(esdTrack->Charge() < 0){
					hNchVsPNegTOF[0][nh]->Fill(esdTrack->P(),multTSdata);
					hNchVsPtNegTOF[0][nh]->Fill(esdTrack->Pt(),multTSdata);
					hNchVsPtDataNegTOF[0][nh]->Fill(esdTrack->P(),beta,multTSdata);	
				}
			}
			else if(TMath::Abs(DPhi-pi)<pi/3.0){
				if(esdTrack->Charge() > 0){
					hNchVsPPosTOF[1][nh]->Fill(esdTrack->P(),multTSdata);
					hNchVsPtPosTOF[1][nh]->Fill(esdTrack->Pt(),multTSdata);
					hNchVsPtDataPosTOF[1][nh]->Fill(esdTrack->P(),beta,multTSdata);
				}
				if(esdTrack->Charge() < 0){
					hNchVsPNegTOF[1][nh]->Fill(esdTrack->P(),multTSdata);
					hNchVsPtNegTOF[1][nh]->Fill(esdTrack->Pt(),multTSdata);
					hNchVsPtDataNegTOF[1][nh]->Fill(esdTrack->P(),beta,multTSdata);	
				}
			}
			else{
				if(esdTrack->Charge() > 0){
					hNchVsPPosTOF[2][nh]->Fill(esdTrack->P(),multTSdata);
					hNchVsPtPosTOF[2][nh]->Fill(esdTrack->Pt(),multTSdata);
					hNchVsPtDataPosTOF[2][nh]->Fill(esdTrack->P(),beta,multTSdata);
				}
				if(esdTrack->Charge() < 0){
					hNchVsPNegTOF[2][nh]->Fill(esdTrack->P(),multTSdata);
					hNchVsPtNegTOF[2][nh]->Fill(esdTrack->Pt(),multTSdata);
					hNchVsPtDataNegTOF[2][nh]->Fill(esdTrack->P(),beta,multTSdata);	
				}
			}

			if(TMath::Abs(fPIDResponse->NumberOfSigmasTOF(esdTrack,AliPID::kPion))<2.0)
				histPiTof[nh]->Fill(esdTrack->P(),esdTrack->GetTPCsignal());


		}	// TOF 

		//
		//_______________________________ rTPC PID
		//

		double momentum = esdTrack->P();
		float  dedx     = esdTrack->GetTPCsignal();

		if(!PhiCut(esdTrack->Pt(), esdTrack->Phi(), esdTrack->Charge(), Magf, fcutLow, fcutHigh))
			continue;

		if(fdEdxCalibrated)
			dedx *= 50/EtaCalibration(eta);

		if( (momentum <= 0.6)&&(momentum >= 0.4) ){//only p:0.4-0.6 GeV, pion MIP
			if( (esdTrack->GetTPCsignal() < fDeDxMIPMax) && (esdTrack->GetTPCsignal() > fDeDxMIPMin) ){
				hMIPVsEta->Fill(eta,dedx);
				pMIPVsEta->Fill(eta,dedx);
			}
			if( (esdTrack->GetTPCsignal() > 70.0) && (esdTrack->GetTPCsignal() < 90.0) ){
				if( TMath::Abs(fPIDResponse->NumberOfSigmasTOF(esdTrack,AliPID::kElectron))<2.0 ){
					hPlateauVsEta->Fill(eta,dedx);
					pPlateauVsEta->Fill(eta,dedx);
				}
			}
		}

		if( momentum <= 0.6 && momentum >= 0.4 ){
			if( dedx < fDeDxMIPMax && dedx > fDeDxMIPMin ){
				hMIPVsPhi[nh]->Fill(esdTrack->Phi(),dedx);
				pMIPVsPhi[nh]->Fill(esdTrack->Phi(),dedx);
			}
			if( dedx > 70 && dedx < 90 ){
				if( TMath::Abs(fPIDResponse->NumberOfSigmasTOF(esdTrack,AliPID::kElectron))<2.0 ){
					hPlateauVsPhi[nh]->Fill(esdTrack->Phi(),dedx);
					pPlateauVsPhi[nh]->Fill(esdTrack->Phi(),dedx);
				}
			}
		}

		if(TMath::Abs(DPhi)<pi/3.0){
			hDeDxVsP[0][nh]->Fill(momentum,dedx,multTSdata);
			hNchVsPrTPC[0][nh]->Fill(momentum,multTSdata);
			hNchVsPtrTPC[0][nh]->Fill(esdTrack->Pt(),multTSdata);
		}
		else if(TMath::Abs(DPhi-pi)<pi/3.0){
			hDeDxVsP[1][nh]->Fill(momentum,dedx,multTSdata);
			hNchVsPrTPC[1][nh]->Fill(momentum,multTSdata);
			hNchVsPtrTPC[1][nh]->Fill(esdTrack->Pt(),multTSdata);
		}
		else{
			hDeDxVsP[2][nh]->Fill(momentum,dedx,multTSdata);
			hNchVsPrTPC[2][nh]->Fill(momentum,multTSdata);
			hNchVsPtrTPC[2][nh]->Fill(esdTrack->Pt(),multTSdata);
		}

		delete esdTrack;

	}//end of track loop
}
//________________________________________________________________________
void AliAnalysisTaskSpectraRT::ProduceArrayV0ESD(){

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

	for (Int_t iV0=0; iV0<nv0s; iV0++) {

		AliESDv0 *esdV0 = fESD->GetV0(iV0);
		if ( !esdV0 ) continue;

		if(esdV0->GetOnFlyStatus()!=0)
			continue;

		// AliESDTrack (V0 Daughters)
		UInt_t lKeyPos = (UInt_t)TMath::Abs(esdV0->GetPindex());
		UInt_t lKeyNeg = (UInt_t)TMath::Abs(esdV0->GetNindex());

		AliESDtrack *pTrack = fESD->GetTrack(lKeyPos);
		AliESDtrack *nTrack = fESD->GetTrack(lKeyNeg);

		if(!pTrack || !nTrack) {
			Printf("ERROR: Could not retreive one of the daughter track");
			continue;
		}

		// Remove like-sign
		if(pTrack->GetSign() == nTrack->GetSign())
			continue;

		// Eta cut on decay products
		if(TMath::Abs(pTrack->Eta()) > fEtaCut || TMath::Abs(nTrack->Eta()) > fEtaCut)
			continue;

		if ( fTrackFilterDaughters ){

			if (!fTrackFilterDaughters->AcceptTrack(pTrack)) continue;
			if (!fTrackFilterDaughters->AcceptTrack(nTrack)) continue;

		}

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

					       bool fillPos = kFALSE;
					       bool fillNeg = kFALSE;

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

						       if(track->GetTPCsignalN()< fNcl )
							       continue;

						       if(!PhiCut(track->Pt(),track->Phi(), track->Charge(), Magf, fcutLow, fcutHigh))
							       continue;

						       double eta      = track->Eta();
						       double momentum = track->P();
						       double dedx     = track->GetTPCsignal();

						       if(fdEdxCalibrated)
							       dedx *= 50/EtaCalibration(eta);						      

						       if(fillPos&&fillNeg){
							       if( (track->GetTPCsignal() < fDeDxMIPMax) && (track->GetTPCsignal() > fDeDxMIPMin) ){
								       if(momentum<0.6&&momentum>0.4){
									       hMIPVsEtaV0s->Fill(eta,dedx);
									       pMIPVsEtaV0s->Fill(eta,dedx);
								       }
							       }
						       }

						       int nh = -1;

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
							       histPiV0[nh]->Fill(track->P(), dedx);
						       }
						       else{
							       histPV0[nh]->Fill(track->P(), dedx);
						       }

					       }//end loop over two tracks

				       };
				       break;

				case 1:{//gammas

					       bool fillPos = kFALSE;
					       bool fillNeg = kFALSE;

					       if( dmassK>0.01 && dmassL>0.01 && dmassAL>0.01 ) {
						       if( dmassG<0.01 && dmassG>0.0001 ) {
							       if( TMath::Abs(nTrack->GetTPCsignal()-EtaCalibrationEl(nTrack->Eta())) < 5)
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

					       double dedx = track->GetTPCsignal();
					       double eta  = track->Eta();

					       if(fdEdxCalibrated)
						       dedx *= 50/EtaCalibration(track->Eta());						      

					       if(track->GetTPCsignalN() <= fNcl)
						       continue;

					       if(!PhiCut(track->Pt(), track->Phi(), track->Charge(), Magf, fcutLow, fcutHigh))
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

					       histEV0[nh]->Fill(track->P(),track->GetTPCsignal());

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
//________________________________________________________________________
bool AliAnalysisTaskSpectraRT::selectVertex2015pp(AliESDEvent *esd,
		Bool_t checkSPDres, //enable check on vtx resolution
		Bool_t requireSPDandTrk, //ask for both trk and SPD vertex
		Bool_t checkProximity) //apply cut on relative position of spd and trk verteces
{

	if (!esd) return kFALSE;

	const AliESDVertex * trkVertex = esd->GetPrimaryVertexTracks();
	const AliESDVertex * spdVertex = esd->GetPrimaryVertexSPD();
	bool hasSPD = spdVertex->GetStatus();
	bool hasTrk = trkVertex->GetStatus();

	//Note that AliVertex::GetStatus checks that N_contributors is > 0
	//reject events if both are explicitly requested and none is available
	if (requireSPDandTrk && !(hasSPD && hasTrk)) return kFALSE;

	//reject events if none between the SPD or track verteces are available
	//if no trk vertex, try to fall back to SPD vertex;
	if (!hasTrk) {
		if (!hasSPD) return kFALSE;
		//on demand check the spd vertex resolution and reject if not satisfied
		if (checkSPDres && !IsGoodSPDvertexRes(spdVertex)) return kFALSE;
	} else {
		if (hasSPD) {
			//if enabled check the spd vertex resolution and reject if not satisfied
			//if enabled, check the proximity between the spd vertex and trak vertex, and reject if not satisfied
			if (checkSPDres && !IsGoodSPDvertexRes(spdVertex)) return kFALSE;
			if ((checkProximity && TMath::Abs(spdVertex->GetZ() - trkVertex->GetZ())>0.5)) return kFALSE;
		}
	}

	//Cut on the vertex z position
	//const AliESDVertex * vertex = esd->GetPrimaryVertex();
	//if (TMath::Abs(vertex->GetZ())>10) return kFALSE;
	return kTRUE;
}
//________________________________________________________________________
bool AliAnalysisTaskSpectraRT::IsGoodSPDvertexRes(const AliESDVertex* spdVertex)
{

	if( !spdVertex ) return kFALSE;
	if( spdVertex->IsFromVertexerZ() && !(spdVertex->GetDispersion()<0.04 && spdVertex->GetZRes()<0.25) ) return kFALSE;
	return kTRUE;
}
//________________________________________________________________________
bool AliAnalysisTaskSpectraRT::IsGoodZvertexPos(AliESDEvent *esd)
{

	if( !esd ) return kFALSE;
	//Cut on the vertex z position
	const AliESDVertex * vertex = esd->GetPrimaryVertex();
	if (TMath::Abs(vertex->GetZ())>10) return kFALSE;
	return kTRUE;
}
//________________________________________________________________________
bool AliAnalysisTaskSpectraRT::PhiCut(const double& pt, double phi, const double& q, const float& mag, TF1* phiCutLow, TF1* phiCutHigh)
{
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

	hPhirTPC->Fill(pt, phi);

	return kTRUE;
}
//________________________________________________________________________
float AliAnalysisTaskSpectraRT::GetMaxDCApTDep( TF1 *fMaxDCAxy, Double_t ptI){

	double maxDCAxy = 10;
	maxDCAxy = fMaxDCAxy->Eval(ptI);
	return maxDCAxy;

}
//________________________________________________________________________
/*void AliAnalysisTaskSpectraRT::SetTrackCuts(AliAnalysisFilter* fTrackFilter){

  AliESDtrackCuts* esdTrackCuts = 0x0;
  if(fSetTPConlyTrkCuts){
  esdTrackCuts = AliESDtrackCuts::GetStandardTPCOnlyTrackCuts();
  esdTrackCuts->SetRequireTPCRefit(kTRUE);
  esdTrackCuts->SetRequireITSRefit(kTRUE);
  esdTrackCuts->SetEtaRange(-0.8,0.8);
  }
  else{
  esdTrackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011(kFALSE,1);
  esdTrackCuts->SetMaxDCAToVertexXYPtDep("0.0105+0.0350/pt^1.1");
  esdTrackCuts->SetEtaRange(-0.8,0.8);
  }

  fTrackFilter->AddCuts(esdTrackCuts);
  }*/
//________________________________________________________________________
AliESDtrack* AliAnalysisTaskSpectraRT::SetHybridTrackCuts(AliESDtrack *esdtrack, const bool fillPhiStand, const bool fillPhHyb1, const bool fillPhHyb2){

	// 
	// 	Get the Hybrid Tracks 
	// 	

	AliESDtrack *newTrack = 0x0;

	if(fTrackFilter->AcceptTrack(esdtrack)){
		newTrack = new AliESDtrack(*esdtrack);
		if(fillPhiStand) hPhiStandard->Fill(newTrack->Eta(),newTrack->Phi());
		////			newTrack->SetTRDQuality(0);
	}
	else if(fHybridTrackCuts1->AcceptTrack(esdtrack)){
		if(esdtrack->GetConstrainedParam()){
			newTrack = new AliESDtrack(*esdtrack);
			const AliExternalTrackParam* constrainParam = esdtrack->GetConstrainedParam();
			newTrack->Set(constrainParam->GetX(),constrainParam->GetAlpha(),constrainParam->GetParameter(),constrainParam->GetCovariance());
			////				newTrack->SetTRDQuality(1);
			if(fillPhHyb1) hPhiHybrid1->Fill(newTrack->Eta(),newTrack->Phi());
		}
		else{ return 0x0; }
	}
	/*else if(fHybridTrackCuts2->AcceptTrack(esdtrack)){
	  if(esdtrack->GetConstrainedParam()){
	  newTrack = new AliESDtrack(*esdtrack);
	  const AliExternalTrackParam* constrainParam = esdtrack->GetConstrainedParam();
	  newTrack->Set(constrainParam->GetX(),constrainParam->GetAlpha(),constrainParam->GetParameter(),constrainParam->GetCovariance());
	/////				newTrack->SetTRDQuality(2);
	if(fillPhHyb2) hPhiHybrid2->Fill(newTrack->Eta(),newTrack->Phi());
	}
	else{ return 0x0; }
	}*/
	  else{
		  return 0x0;
	  }

	  return newTrack;

}
//________________________________________________________________________
double AliAnalysisTaskSpectraRT::EtaCalibration(const double &eta){

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

	if(strcmp(fPeriod,"16l")==0){
		aPos = 49.9216; bPos = 3.07252; cPos = -42.8044; dPos = 259.666; ePos = -910.432; fPos = 1776.09; gPos = -1740.65; hPos = 662.232;
		aNeg = 49.9732; bNeg = 4.03575; cNeg = 65.6189;  dNeg = 374.429; eNeg = 951.459;  fNeg = 1153.75; gNeg = 618.493;  hNeg = 100.499;
	}else if(strcmp(fPeriod,"16k")==0){
		aPos = 49.9421; bPos = 2.3446; cPos = -41.2765; dPos = 279.695; ePos = -1027.73; fPos = 2022.84; gPos = -1967.79; hPos = 738.823;
		aNeg = 50.0477; bNeg = 8.27344; cNeg = 125.29;  dNeg = 736.8;   eNeg = 2057.75;  fNeg = 2935.38; gNeg = 2064.03;  hNeg = 565.983;
	}else if(strcmp(fPeriod,"16deghijop")==0){
		aPos = 49.9743; bPos = 2.3388; cPos = -44.1496; dPos = 296.029; ePos = -1056.56; fPos = 2031.44; gPos = -1946.51; hPos = 723.89;
		aNeg = 50.0329; bNeg = 6.99747; cNeg = 107.168;  dNeg = 649.001; eNeg = 1875.17;  fNeg = 2785.78; gNeg = 2063.77;  hNeg = 606.868;
	}else if(strcmp(fPeriod,"17data")==0){
		aPos = 49.6097; bPos = 0.922856; cPos = -6.57484; dPos = 65.3117; ePos = -372.142; fPos = 950.451; gPos = -1085.27; hPos = 458.144;
		aNeg = 49.6555; bNeg = 6.98696; cNeg = 102.734;  dNeg = 566.424; eNeg = 1513.64;  fNeg = 2092.01; gNeg = 1429.32;  hNeg = 375.642;
	}else{
		aPos = 49.6975; bPos = 2.32535; cPos = -42.6516; dPos = 283.058; ePos = -1009.58; fPos = 1945.89; gPos = -1871.23; hPos = 698.552;
		aNeg = 49.8071; bNeg = 9.78466; cNeg = 120.018;  dNeg = 603.325; eNeg = 1470.92;  fNeg = 1819.63; gNeg = 1073.82;  hNeg = 230.142;
	}

	for(int i=0; i<8; ++i){
		fEtaCalibrationPos->SetParameter(i,0);
		fEtaCalibrationNeg->SetParameter(i,0);
	}

	if(eta<0.0){
		fEtaCalibrationNeg->SetParameter(0,aNeg);
		fEtaCalibrationNeg->SetParameter(1,bNeg);
		fEtaCalibrationNeg->SetParameter(2,cNeg);
		fEtaCalibrationNeg->SetParameter(3,dNeg);
		fEtaCalibrationNeg->SetParameter(4,eNeg);
		fEtaCalibrationNeg->SetParameter(5,fNeg);
		fEtaCalibrationNeg->SetParameter(6,gNeg);
		fEtaCalibrationNeg->SetParameter(7,hNeg);

		return fEtaCalibrationNeg->Eval(eta);
	}
	else{
		fEtaCalibrationPos->SetParameter(0,aPos);
		fEtaCalibrationPos->SetParameter(1,bPos);
		fEtaCalibrationPos->SetParameter(2,cPos);
		fEtaCalibrationPos->SetParameter(3,dPos);
		fEtaCalibrationPos->SetParameter(4,ePos);
		fEtaCalibrationPos->SetParameter(5,fPos);
		fEtaCalibrationPos->SetParameter(6,gPos);
		fEtaCalibrationPos->SetParameter(7,hPos);

		return fEtaCalibrationPos->Eval(eta);
	}


}
//________________________________________________________________________
double AliAnalysisTaskSpectraRT::EtaCalibrationEl(const double &eta){

	double aPos = 0.0;
	double bPos = 0.0;
	double cPos = 0.0;
	double dPos = 0.0;
	double ePos = 0.0;

	double aNeg = 0.0;
	double bNeg = 0.0;
	double cNeg = 0.0;
	double dNeg = 0.0;
	double eNeg = 0.0;

	if(strcmp(fPeriod,"16l")==0){
		aPos = 79.4195; bPos = 7.82459; cPos = -23.3466; dPos = 26.5577; ePos = -8.27151;
		aNeg = 79.8571; bNeg = -14.2921; cNeg = -66.6972; dNeg = -103.794; eNeg = -50.5771;
	}else if(strcmp(fPeriod,"16k")==0){
		aPos = 80.254; bPos = 6.37076; cPos = -50.9878; dPos = 116.611; ePos = -79.0483;
		aNeg = 79.8728; bNeg = -3.08265; cNeg = -11.3778; dNeg = -20.6605; eNeg = -12.3861;
	}else if(strcmp(fPeriod,"16deghijop")==0){
		aPos = 80.0719; bPos = 7.10053; cPos = -42.4788; dPos = 86.1074; ePos = -54.0891;
		aNeg = 79.6155; bNeg = -12.1254; cNeg = -66.2488; dNeg = -132.426; eNeg = -85.0155;
	}else if(strcmp(fPeriod,"17data")==0){
		aPos = 82.4621; bPos = 5.20353; cPos = -32.2608; dPos = 63.4788; ePos = -39.3277;
		aNeg = 82.306; bNeg = -4.04076; cNeg = -22.133; dNeg = -40.5782; eNeg = -23.8157;
	}else{
		aPos = 79.7726; bPos = 6.83744; cPos = -40.0469; dPos = 78.987; ePos = -50.1373;
		aNeg = 79.4863; bNeg = -5.00403; cNeg = -21.6184;  dNeg = -39.1295; eNeg = -24.8757;
	}

	for(int i=0; i<5; ++i){
		fEtaCalibrationPosEl->SetParameter(i,0);
		fEtaCalibrationNegEl->SetParameter(i,0);
	}

	if(eta<0.0){
		fEtaCalibrationNegEl->SetParameter(0,aNeg);
		fEtaCalibrationNegEl->SetParameter(1,bNeg);
		fEtaCalibrationNegEl->SetParameter(2,cNeg);
		fEtaCalibrationNegEl->SetParameter(3,dNeg);
		fEtaCalibrationNegEl->SetParameter(4,eNeg);

		return fEtaCalibrationNegEl->Eval(eta);
	}
	else{
		fEtaCalibrationPosEl->SetParameter(0,aPos);
		fEtaCalibrationPosEl->SetParameter(1,bPos);
		fEtaCalibrationPosEl->SetParameter(2,cPos);
		fEtaCalibrationPosEl->SetParameter(3,dPos);
		fEtaCalibrationPosEl->SetParameter(4,ePos);

		return fEtaCalibrationPosEl->Eval(eta);
	}


}
//________________________________________________________________________
bool AliAnalysisTaskSpectraRT::TOFPID(AliESDtrack * track)
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

