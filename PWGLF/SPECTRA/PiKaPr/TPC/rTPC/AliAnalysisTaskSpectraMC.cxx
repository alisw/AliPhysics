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

#include "AliAnalysisTaskSpectraMC.h"

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
////class AliAnalysisTaskSpectraMC;    // your analysis class
using namespace std;

//
// Responsible:
// Omar Vazquez (Lund)

static float Magf = 1;
static const int nPid = 4;
static const int nRegion = 4;
static const int nRt = 5;
static const char* Region[4] = {"Toward","Away","Transverse","FullAzimuth"};
static const char* Pid[nPid+1] = {"Charged","Pion","Kaon","Proton","Phi"};

ClassImp(AliAnalysisTaskSpectraMC)
	AliAnalysisTaskSpectraMC::AliAnalysisTaskSpectraMC():
		AliAnalysisTaskSE(),
		fESD(0x0),
		fEventCuts(0x0),
		fMC(0x0),
		fMCStack(0x0),
		fMCArray(0x0),
		fPIDResponse(0x0),
		fGeometricalCut(0x0),
		fTrackFilter(0x0),
		fHybridTrackCuts1(0x0),
		fHybridTrackCuts2(0x0),
		utils(0x0),
		fAnalysisType("ESD"),
		fAnalysisMC(kFALSE),
		fIsMCclosure(kTRUE),
		fRandom(0x0),
		fTrackID(0),
		fNcl(70),
		fEtaCut(0.9),
		fDeDxMIPMin(40),
		fDeDxMIPMax(60),
		fdEdxHigh(200),
		fdEdxLow(40),
		//fSetTPConlyTrkCuts(kFALSE),
		fSelectHybridTracks(kTRUE),
		fLeadPtCutMin(5.0),
		fLeadPtCutMax(40.0),
		fGenLeadPhi(0.0),
		fGenLeadPt(0.0),
		fGenLeadEta(0x0),
		fRecLeadEta(0x0), 
		fGenLeadIn(0.0),
		fRecLeadPhi(0.0),
		fRecLeadPt(0.0),
		fRecLeadIn(0.0),
		fPtMin(0.15),
		fListOfObjects(0),
		fEvents(0x0),
		hPtRec(0x0),
		hPtRecPion(0x0),
		hPtRecProton(0x0),
		hPtPri(0x0),
		hPtSec(0x0),
		hInvMassPhi(0x0),
		fPtLVsNchRec(0x0),
		hDeltaPtLeading(0x0),
		hPtRecEtaRecLeading(0x0),
		hEtaGenEtaRecLeading(0x0), 
		hPhiEtaRecLeading(0x0),
		hPhiLeading(0x0),
		hMultTSGen(0x0),
		hMultTSRec(0x0),
		hTrigger(0x0),
		hNchResponse(0x0),
		hPhiTotal(0x0),
		hPhiResPhi(0x0),
		hPhiStandard(0x0),
		hPhiHybrid1(0x0),
		hPhiHybrid2(0x0),
		fEtaCalibration(0x0),
		fEtaCalibrationEl(0x0),
		fcutDCAxy(0x0),
		fcutLow(0x0),
		fcutHigh(0x0)

{

	for(int pid = 0; pid < nPid; ++pid){
		hPtGenPosIn[pid] = 0;
		hPtRecPosIn[pid] = 0;
		hPtGenNegIn[pid] = 0;
		hPtRecNegIn[pid] = 0;
		hPtRecInTOF[pid] = 0;
		hPtRecPosInTOF[pid] = 0;
		hPtRecNegInTOF[pid] = 0;
		hPtrTPCRecIn[pid] = 0;
	}

	for(int pid = 0; pid < (nPid+1); ++pid){
		hPtGenIn[pid] = 0;
		hPtRecIn[pid] = 0;
	}

	for(int r = 0; r < nRegion; r++){

		hPhiGen[r] = 0;
		hPhiRec[r] = 0;

		for(int pid = 0; pid < (nPid+1); ++pid){

			hNchGenVsPtGenPID[r][pid] = 0;
			hNchGenVsPtRec[r][pid] = 0;
			hNchVsPtDataTPC[r][pid] = 0;


		}

	}	// region        


}


AliAnalysisTaskSpectraMC::AliAnalysisTaskSpectraMC(const char *name):
	AliAnalysisTaskSE(name),
	fESD(0x0),
	fEventCuts(0x0),
	fMC(0x0),
	fMCStack(0x0),
	fMCArray(0x0),
	fPIDResponse(0x0),
	fGeometricalCut(0x0),
	fTrackFilter(0x0),
	fHybridTrackCuts1(0x0),
	fHybridTrackCuts2(0x0),
	utils(0x0),
	fAnalysisType("ESD"),
	fAnalysisMC(kFALSE),
	fIsMCclosure(kTRUE),
	fRandom(0x0),
	fTrackID(0),
	fNcl(70),
	fEtaCut(0.9),
	fDeDxMIPMin(40),
	fDeDxMIPMax(60),
	fdEdxHigh(200),
	fdEdxLow(40),
	//fSetTPConlyTrkCuts(kFALSE),
	fSelectHybridTracks(kTRUE),
	fLeadPtCutMin(5.0),
	fLeadPtCutMax(40.0),
	fGenLeadPhi(0.0),
	fGenLeadPt(0.0),
	fGenLeadEta(0x0),
	fRecLeadEta(0x0), 
	fGenLeadIn(0.0),
	fRecLeadPhi(0.0),
	fRecLeadPt(0.0),
	fRecLeadIn(0.0),
	fPtMin(0.15),
	fListOfObjects(0),
	fEvents(0x0),
	hPtRec(0x0),
	hPtRecPion(0x0),
	hPtRecProton(0x0),
	hPtPri(0x0),
	hPtSec(0x0),
	hInvMassPhi(0x0),
	fPtLVsNchRec(0x0),
	hDeltaPtLeading(0x0),
	hPtRecEtaRecLeading(0x0),
	hEtaGenEtaRecLeading(0x0), 
	hPhiEtaRecLeading(0x0),
	hPhiLeading(0x0),
	hMultTSGen(0x0),
	hMultTSRec(0x0),
	hTrigger(0x0),
	hNchResponse(0x0),
	hPhiTotal(0x0),
	hPhiResPhi(0x0),
	hPhiStandard(0x0),
	hPhiHybrid1(0x0),
	hPhiHybrid2(0x0),
	fEtaCalibration(0x0),
	fEtaCalibrationEl(0x0),
	fcutDCAxy(0x0),
	fcutLow(0x0),
	fcutHigh(0x0)

{

	for(int pid = 0; pid < nPid; ++pid){
		hPtGenPosIn[pid] = 0;
		hPtRecPosIn[pid] = 0;
		hPtGenNegIn[pid] = 0;
		hPtRecNegIn[pid] = 0;
		hPtRecInTOF[pid] = 0;
		hPtRecPosInTOF[pid] = 0;
		hPtRecNegInTOF[pid] = 0;
		hPtrTPCRecIn[pid] = 0;
	}

	for(int pid = 0; pid < (nPid+1); ++pid){
		hPtGenIn[pid] = 0;
		hPtRecIn[pid] = 0;
	}

	for(int r = 0; r < nRegion; r++){

		hPhiGen[r] = 0;
		hPhiRec[r] = 0;

		for(int pid = 0; pid < (nPid+1); ++pid){

			hNchGenVsPtGenPID[r][pid] = 0;
			hNchGenVsPtRec[r][pid] = 0;
			hNchVsPtDataTPC[r][pid] = 0;

		}

	}        

	DefineInput(0, TChain::Class());
	DefineOutput(1, TList::Class());//esto es nuevo
}

AliAnalysisTaskSpectraMC::~AliAnalysisTaskSpectraMC() {
	//
	// Destructor
	//

}
//______________________________________________________________________________
void AliAnalysisTaskSpectraMC::UserCreateOutputObjects()
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
		//fHybridTrackCuts1->SetMinRatioCrossedRowsOverFindableClustersTPC(0.8); //! Variated in track cuts systematics
		//fHybridTrackCuts1->SetMaxChi2PerClusterTPC(4); //! Variated in track cuts systematics
		fHybridTrackCuts1->SetAcceptKinkDaughters(kFALSE);
		fHybridTrackCuts1->SetRequireTPCRefit(kTRUE);
		fHybridTrackCuts1->SetRequireITSRefit(kTRUE);
		fHybridTrackCuts1->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kOff);
		fHybridTrackCuts1->SetMaxDCAToVertexXYPtDep("0.0105+0.0350/pt^1.1");
		//fHybridTrackCuts1->SetMaxChi2TPCConstrainedGlobal(36); //! This cut is excluded
		//fHybridTrackCuts1->SetMaxDCAToVertexZ(2); //! Variated in track cuts systematics
		fHybridTrackCuts1->SetDCAToVertex2D(kFALSE);
		fHybridTrackCuts1->SetRequireSigmaToVertex(kFALSE);
		//fHybridTrackCuts1->SetMaxChi2PerClusterITS(36); //! Variated in track cuts systematics
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
		//fHybridTrackCuts2->SetMinRatioCrossedRowsOverFindableClustersTPC(0.8); //! Variated in track cuts systematics
		//fHybridTrackCuts2->SetMaxChi2PerClusterTPC(4); //! Variated in track cuts systematics
		fHybridTrackCuts2->SetAcceptKinkDaughters(kFALSE);
		fHybridTrackCuts2->SetRequireTPCRefit(kTRUE);
		fHybridTrackCuts2->SetRequireITSRefit(kFALSE);
		fHybridTrackCuts2->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kOff);
		fHybridTrackCuts2->SetMaxDCAToVertexXYPtDep("0.0105+0.0350/pt^1.1");
		//fHybridTrackCuts2->SetMaxChi2TPCConstrainedGlobal(36); //! This cut is excluded
		//fHybridTrackCuts2->SetMaxDCAToVertexZ(2); //! Variated in track cuts systematics
		fHybridTrackCuts2->SetDCAToVertex2D(kFALSE);
		fHybridTrackCuts2->SetRequireSigmaToVertex(kFALSE);
		//fHybridTrackCuts2->SetMaxChi2PerClusterITS(36); //! Variated in track cuts systematics
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


	const int nPtBins = 45;
	double ptBins[nPtBins+1] = {
		0.0, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45,
		0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85,
		0.90, 0.95, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7,
		1.80, 1.90, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4,
		3.60, 3.80, 4.0, 4.5, 5.0, 5.5, 6.0, 7.0, 8.0, 10.0};


	// Phi binning

	/*	const int nPtBins = 17;
		double ptBins[nPtBins+1] = { 0.0 , 0.15, 0.3, 0.5, 0.7, 0.9, 1.2, 1.4, 1.6,
		1.8, 2.0, 2.2, 2.6, 3.0, 3.5, 4.0, 5.0, 8.0};
		*/

	const int nBinsNsigma = 50;
	double binsNsigma[nBinsNsigma+1] = {0};

	for(int i = 0; i <= nBinsNsigma; ++i){
		binsNsigma[i] = -10.0+i*0.4;
	}

	const int nDeltaPiBins   = 55;
	double DeltaPiBins[nDeltaPiBins+1] = { 0 };
	for(int i = 0; i <= nDeltaPiBins; ++i){
		DeltaPiBins[i] = 20.0+i*1.0;
	}

	const int ndEdxBins   = 100;
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
		binsRT[i] = (((double)i)-0.5);

	fcutDCAxy = new TF1("fMaxDCAxy","[0]+[1]/(x^[2])",0,1e10);
	fcutDCAxy->SetParameter(0,0.0105);
	fcutDCAxy->SetParameter(1,0.0350);
	fcutDCAxy->SetParameter(2,1.1);

	fcutLow = new TF1("StandardPhiCutLow",  "0.1/x/x+TMath::Pi()/18.0-0.025", 0, 50);
	fcutHigh = new TF1("StandardPhiCutHigh", "0.12/x+TMath::Pi()/18.0+0.035", 0, 50);

	fEtaCalibration   = new TF1("fDeDxVsEtaPos", "pol7", 0.0, 1.0);
	fEtaCalibrationEl = new TF1("fDeDxVsEtaEl", "pol4", 0.0, 1.0);

	hPhiTotal = new TH2F("hPhiSum","; #eta; #varphi",50,-0.8,0.8,100,-TMath::Pi()/2.0,5.0*TMath::Pi()/2.0);
	fListOfObjects->Add(hPhiTotal);

	hPhiResPhi = new TH1F("hPhiResPhi","; #varphi",100,-TMath::Pi()/2.0,5.0*TMath::Pi()/2.0);
	fListOfObjects->Add(hPhiResPhi);

	hPhiStandard = new TH2F("hPhiStandard","; #eta; #varphi",50,-0.8,0.8,100,-TMath::Pi()/2.0,5.0*TMath::Pi()/2.0);
	fListOfObjects->Add(hPhiStandard);

	hPhiHybrid1 = new TH2F("hPhiHybrid1","; #eta; #varphi",50,-0.8,0.8,100,-TMath::Pi()/2.0,5.0*TMath::Pi()/2.0);
	fListOfObjects->Add(hPhiHybrid1);

	hPhiHybrid2 = new TH2F("hPhiHybrid2","; #eta; #varphi",50,-0.8,0.8,100,-TMath::Pi()/2.0,5.0*TMath::Pi()/2.0);
	fListOfObjects->Add(hPhiHybrid2);

	if(fAnalysisMC){

		hMultTSRec = new TH1F("hMultTSRec",";#it{N}_{acc}; Entries",nBinsRT,binsRT);
		fListOfObjects->Add(hMultTSRec);

		hTrigger = new TH1F("hTrigger",";#it{N}_{ev}; Entries",7,0.0,7.0);
		hTrigger->GetXaxis()->SetBinLabel(1,"Good Gen");
		hTrigger->GetXaxis()->SetBinLabel(2,"Good Rec");
		hTrigger->GetXaxis()->SetBinLabel(3,"Good Both");
		hTrigger->GetXaxis()->SetBinLabel(4,"Good Gen/Bad Rec");
		hTrigger->GetXaxis()->SetBinLabel(5,"Good Rec/Bad Gen");
		fListOfObjects->Add(hTrigger);

		hMultTSGen = new TH1F("hMultTSGen",";#it{N}_{acc}; Entries",nBinsRT,binsRT);
		fListOfObjects->Add(hMultTSGen);

		fPtLVsNchRec = new TH2F("fPtLVsNchRec", ";#it{p}^{L}_{rec} (GeV/#it{c});#it{N}_{acc}",nPtBins,ptBins,nBinsRT,binsRT);
		fListOfObjects->Add(fPtLVsNchRec);

		hDeltaPtLeading = new TH1F("hDeltaPtLeading",";#Delta#it{p}_{T} (GeV/#it{c}); Entries",80,-10.0,10.0);
		fListOfObjects->Add(hDeltaPtLeading);

		hPtRecEtaRecLeading = new TH2F("hPtRecEtaRecLeading",";#it{p}_{T}^{r};#eta^{r}",nPtBins,ptBins,100,-1.0,1.0); 
		fListOfObjects->Add(hPtRecEtaRecLeading);

		hPhiLeading = new TH1F("hPhiLeading",";#varphi^{leading}; Entries",100,-TMath::Pi()/2.0,5.0*TMath::Pi()/2.0);
		fListOfObjects->Add(hPhiLeading);

		hEtaGenEtaRecLeading = new TH2F("hDEtahDPt",";#eta^{g}-#eta^{r}; #it{p}_{T}^{g}-#it{p}_{T}^{r}",100,-2.0,2.0,80,-10.0,10.0); 
		fListOfObjects->Add(hEtaGenEtaRecLeading);

		hPhiEtaRecLeading = new TH2F("hDEtahDPhi",";#eta^{g}-#eta^{r}; #varphi^{g}-#varphi^{r}",100,-2.0,2.0,100,-TMath::Pi()/2.0,5.0*TMath::Pi()/2.0); 
		fListOfObjects->Add(hPhiEtaRecLeading);

		for( int i = 0; i < nRegion; ++i ){

			hPhiGen[i]= new TH1F(Form("hPhiGen_%s",Region[i]),"",64,-TMath::Pi()/2.0,3.0*TMath::Pi()/2.0);
			fListOfObjects->Add(hPhiGen[i]);

			hPhiRec[i] = new TH1F(Form("hPhiRec_%s",Region[i]),"",64,-TMath::Pi()/2.0,3.0*TMath::Pi()/2.0);
			fListOfObjects->Add(hPhiRec[i]);

			for(int pid = 0; pid < (nPid+1); ++pid){

				hNchVsPtDataTPC[i][pid] = new TH2F(Form("hNchVsPt_%s_%s",Region[i],Pid[pid]),";#it{N}_{acc};#it{p}_{T}^{rec}",nBinsRT,binsRT,nPtBins,ptBins);
				fListOfObjects->Add(hNchVsPtDataTPC[i][pid]);

				hNchGenVsPtGenPID[i][pid] = new TH2F(Form("hNchGenVsPtGen_%s_%s",Region[i],Pid[pid]),";#it{N}_{ch};#it{p}_{T}^{gen};",nBinsRT,binsRT,nPtBins,ptBins);
				fListOfObjects->Add(hNchGenVsPtGenPID[i][pid]);

				hNchGenVsPtRec[i][pid] = new TH2F(Form("hNchGenVsPtRec_%s_%s",Region[i],Pid[pid]),";#it{N}_{ch};#it{p}_{T}^{rec};",nBinsRT,binsRT,nPtBins,ptBins);
				fListOfObjects->Add(hNchGenVsPtRec[i][pid]);

			}
		}

		for(int pid = 0; pid < nPid; pid++){

			hPtGenPosIn[pid] = new TH1F(Form("hPtGenPosIn_%s",Pid[pid]),";#it{p}_{T}^{gen}",nPtBins,ptBins);
			fListOfObjects->Add(hPtGenPosIn[pid]);

			hPtGenNegIn[pid] = new TH1F(Form("hPtGenNegIn_%s",Pid[pid]),";#it{p}_{T}^{gen}",nPtBins,ptBins);
			fListOfObjects->Add(hPtGenNegIn[pid]);

			hPtRecPosIn[pid] = new TH1F(Form("hPtRecPosIn_%s",Pid[pid]),";#it{p}_{T}^{rec}",nPtBins,ptBins);
			fListOfObjects->Add(hPtRecPosIn[pid]);

			hPtRecNegIn[pid] = new TH1F(Form("hPtRecNegIn_%s",Pid[pid]),";#it{p}_{T}^{rec}",nPtBins,ptBins);
			fListOfObjects->Add(hPtRecNegIn[pid]);

			hPtRecInTOF[pid] = new TH1F(Form("hPtRecInTOF_%s",Pid[pid]),";#it{p}_{T}^{rec}",nPtBins,ptBins);
			fListOfObjects->Add(hPtRecInTOF[pid]);

			hPtRecPosInTOF[pid] = new TH1F(Form("hPtRecPosInTOF_%s",Pid[pid]),";#it{p}_{T}^{rec}",nPtBins,ptBins);
			fListOfObjects->Add(hPtRecPosInTOF[pid]);

			hPtRecNegInTOF[pid] = new TH1F(Form("hPtRecNegInTOF_%s",Pid[pid]),";#it{p}_{T}^{rec}",nPtBins,ptBins);
			fListOfObjects->Add(hPtRecNegInTOF[pid]);

			hPtrTPCRecIn[pid] = new TH1F(Form("hPtrTPCRecIn_%s",Pid[pid]),";#it{p}_{T}^{rec}",nPtBins,ptBins);
			fListOfObjects->Add(hPtrTPCRecIn[pid]);

		}

		for(int pid = 0; pid < (nPid+1); pid++){

			hPtGenIn[pid] = new TH1F(Form("hPtGenIn_%s",Pid[pid]),";#it{p}_{T}^{gen}",nPtBins,ptBins);
			fListOfObjects->Add(hPtGenIn[pid]);

			hPtRecIn[pid] = new TH1F(Form("hPtRecIn_%s",Pid[pid]),";#it{p}_{T}^{rec}",nPtBins,ptBins);
			fListOfObjects->Add(hPtRecIn[pid]);

		}

		hNchResponse = new TH2F("hNchResponse","; #it{N}_{ch}; #it{N}_{acc}",nBinsRT,binsRT,nBinsRT,binsRT);
		fListOfObjects->Add(hNchResponse);

		hPtRec = new TH1F("hPtRec","; #it{p}_{T}; Entries",nPtBins,ptBins);
		fListOfObjects->Add(hPtRec);

		hPtRecPion = new TH1F("hPtRecPion","; #it{p}_{T}; Entries",nPtBins,ptBins);
		fListOfObjects->Add(hPtRecPion);

		hPtRecProton = new TH1F("hPtRecProton","; #it{p}_{T}; Entries",nPtBins,ptBins);
		fListOfObjects->Add(hPtRecProton);

		hPtPri = new TH1F("hPtPri","; #it{p}_{T}^{rec}; Entries",nPtBins,ptBins);
		fListOfObjects->Add(hPtPri);

		hPtSec = new TH1F("hPtSec","; #it{p}_{T}^{rec}; Entries",nPtBins,ptBins);
		fListOfObjects->Add(hPtSec);

		hInvMassPhi = new TH1F("hInvMassPhi","; InvMass_{#phi}; Entries",50,0.95,1.35);
		fListOfObjects->Add(hInvMassPhi);

	}

	////	fEventCuts.AddQAplotsToList(fListOfObjects);
	PostData(1, fListOfObjects);

}
//______________________________________________________________________________
void AliAnalysisTaskSpectraMC::UserExec(Option_t *)
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
				}

				GetMCCorrections();
				GetMCCorrectionsPhi();

			}
		}
		else{// for testing the method
			if( ( fGenLeadPt>=fLeadPtCutMin && fGenLeadPt<fLeadPtCutMax ) && ( fRecLeadPt>=fLeadPtCutMin && fRecLeadPt<fLeadPtCutMax )){
				GetMultiplicityDistributions();	
				//GetMultiplicityDistributionsPhi();
			}

			if(( fGenLeadPt>=fLeadPtCutMin && fGenLeadPt<fLeadPtCutMax ) ) { hTrigger->Fill(0.5); }
			if(( fRecLeadPt>=fLeadPtCutMin && fRecLeadPt<fLeadPtCutMax ) ) { hTrigger->Fill(1.5); }
			if( ( fGenLeadPt>=fLeadPtCutMin && fGenLeadPt<fLeadPtCutMax ) && ( fRecLeadPt>=fLeadPtCutMin && fRecLeadPt<fLeadPtCutMax )) { hTrigger->Fill(2.5); }
			if( ( fGenLeadPt>=fLeadPtCutMin && fGenLeadPt<fLeadPtCutMax ) && !( fRecLeadPt>=fLeadPtCutMin && fRecLeadPt<fLeadPtCutMax )){ hTrigger->Fill(3.5); }
			if( !( fGenLeadPt>=fLeadPtCutMin && fGenLeadPt<fLeadPtCutMax ) && ( fRecLeadPt>=fLeadPtCutMin && fRecLeadPt<fLeadPtCutMax )) 
			{ hTrigger->Fill(4.5); 
				hPtRecEtaRecLeading->Fill(fRecLeadPt,fRecLeadEta); 
				hDeltaPtLeading->Fill(fGenLeadPt-fRecLeadPt); 
				hEtaGenEtaRecLeading->Fill(fGenLeadEta-fRecLeadEta,fGenLeadPt-fRecLeadPt); 
				hPhiEtaRecLeading->Fill(fGenLeadEta-fRecLeadEta,fGenLeadPhi-fRecLeadPhi); }

		}
	}
	else{
		if(fAnalysisMC){
			if(isGoodVtxPosMC){
				if( (fGenLeadPt>=fLeadPtCutMin && fGenLeadPt<fLeadPtCutMax)&&(fRecLeadPt>=fLeadPtCutMin && fRecLeadPt<fLeadPtCutMax)){
					GetMultiplicityDistributions();
					//GetMultiplicityDistributionsPhi();
					GetDetectorResponse();
				}

				if(( fGenLeadPt>=fLeadPtCutMin && fGenLeadPt<fLeadPtCutMax ) ) { hTrigger->Fill(0.5); }
				if(( fRecLeadPt>=fLeadPtCutMin && fRecLeadPt<fLeadPtCutMax ) ) { hTrigger->Fill(1.5); }
				if( ( fGenLeadPt>=fLeadPtCutMin && fGenLeadPt<fLeadPtCutMax ) && ( fRecLeadPt>=fLeadPtCutMin && fRecLeadPt<fLeadPtCutMax )) { hTrigger->Fill(2.5); }
				if( ( fGenLeadPt>=fLeadPtCutMin && fGenLeadPt<fLeadPtCutMax ) && !( fRecLeadPt>=fLeadPtCutMin && fRecLeadPt<fLeadPtCutMax )){ hTrigger->Fill(3.5); }
				if( !( fGenLeadPt>=fLeadPtCutMin && fGenLeadPt<fLeadPtCutMax ) && ( fRecLeadPt>=fLeadPtCutMin && fRecLeadPt<fLeadPtCutMax ))
				{ hTrigger->Fill(4.5); 
					hDeltaPtLeading->Fill(fGenLeadPt-fRecLeadPt); 
					hPtRecEtaRecLeading->Fill(fRecLeadPt,fRecLeadEta); 
					hEtaGenEtaRecLeading->Fill(fGenLeadEta-fRecLeadEta,fGenLeadPt-fRecLeadPt); 
					hPhiEtaRecLeading->Fill(fGenLeadEta-fRecLeadEta,fGenLeadPhi-fRecLeadPhi); }

				GetMCCorrections();
				//GetMCCorrectionsPhi();

			}
		}

	}

	PostData(1, fListOfObjects);
}
//_____________________________________________________________________________
void AliAnalysisTaskSpectraMC::GetLeadingObject(bool isMC) {

	double flPt = 0.0;
	double flPhi = 0.0;
	double flEta = 0.0;
	int flIndex = 0;

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
				flEta = particle->Eta();
				flIndex = i;
			}
		}

		fGenLeadPhi = flPhi;
		fGenLeadEta = flEta;
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
				else{ track_hybrid = track; }
			}else{
				track_hybrid = SetHybridTrackCuts(track,kFALSE,kFALSE,kFALSE);
				if(!track_hybrid) { continue; }
			}

			if(!fGeometricalCut->AcceptTrack(track_hybrid)) continue;

			if (flPt<track_hybrid->Pt()){
				flPt  = track_hybrid->Pt();
				flPhi = track_hybrid->Phi();
				flEta = track_hybrid->Eta();
				flIndex = i;
			}

			delete track_hybrid;

		}

		fRecLeadPhi = flPhi;
		fRecLeadEta = flEta;
		fRecLeadPt  = flPt;
		fRecLeadIn  = flIndex;

	}

}
//_____________________________________________________________________________
void AliAnalysisTaskSpectraMC::GetMultiplicityDistributions(){

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

		AliESDtrack* esdtrack = static_cast<AliESDtrack*>(fESD->GetTrack(i));  
		if(!esdtrack) continue;
		if(i==fRecLeadIn) continue;
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

	hMultTSGen->Fill(multTSgen);
	hMultTSRec->Fill(multTSrec);
	fPtLVsNchRec->Fill(fRecLeadPt,multTSrec);

	//	Filling Nch vs pT True
	for (int i = 0; i < fMC->GetNumberOfTracks(); i++) {

		AliMCParticle* particle = (AliMCParticle*)fMC->GetTrack(i);
		if(!particle) continue;
		if(i==fGenLeadIn) continue;
		if(TMath::Abs(particle->Eta()) > fEtaCut) continue;
		if(particle->Pt() < fPtMin) continue;

		short pidCodeMC = 0;
		pidCodeMC = GetPidCode(particle->PdgCode());

		double DPhi = DeltaPhi(particle->Phi(), fGenLeadPhi);

		if(particle->Charge() != 0.0){ // Charge != 0.0 particles

			if(!fMC->IsPhysicalPrimary(i)) // Here, otherwise the Phi will be cut off 
				continue;

			hNchGenVsPtGenPID[3][0]->Fill(multTSgen,particle->Pt());
			if(pidCodeMC==1)
				hNchGenVsPtGenPID[3][1]->Fill(multTSgen,particle->Pt());
			if(pidCodeMC==2)
				hNchGenVsPtGenPID[3][2]->Fill(multTSgen,particle->Pt());
			if(pidCodeMC==3)
				hNchGenVsPtGenPID[3][3]->Fill(multTSgen,particle->Pt());

			if(TMath::Abs(DPhi)<pi/3.0){

				hNchGenVsPtGenPID[0][0]->Fill(multTSgen,particle->Pt());

				if(pidCodeMC==1)
					hNchGenVsPtGenPID[0][1]->Fill(multTSgen,particle->Pt());
				else if(pidCodeMC==2)
					hNchGenVsPtGenPID[0][2]->Fill(multTSgen,particle->Pt());
				else if(pidCodeMC==3)
					hNchGenVsPtGenPID[0][3]->Fill(multTSgen,particle->Pt());
				else
					continue;
			}
			else if(TMath::Abs(DPhi-pi)<pi/3.0){

				hNchGenVsPtGenPID[1][0]->Fill(multTSgen,particle->Pt());

				if(pidCodeMC==1)
					hNchGenVsPtGenPID[1][1]->Fill(multTSgen,particle->Pt());
				else if(pidCodeMC==2)
					hNchGenVsPtGenPID[1][2]->Fill(multTSgen,particle->Pt());
				else if(pidCodeMC==3)
					hNchGenVsPtGenPID[1][3]->Fill(multTSgen,particle->Pt());
				else
					continue;
			}
			else{

				hNchGenVsPtGenPID[2][0]->Fill(multTSgen,particle->Pt());

				if(pidCodeMC==1)
					hNchGenVsPtGenPID[2][1]->Fill(multTSgen,particle->Pt());
				else if(pidCodeMC==2)
					hNchGenVsPtGenPID[2][2]->Fill(multTSgen,particle->Pt());
				else if(pidCodeMC==3)
					hNchGenVsPtGenPID[2][3]->Fill(multTSgen,particle->Pt());
				else
					continue;
			}
		} // Charge != 0.0

		else{ // Charge == 0.0 particles

			if(pidCodeMC==6)
				hNchGenVsPtGenPID[3][4]->Fill(multTSgen,particle->Pt());

			if(TMath::Abs(DPhi)<pi/3.0){

				if(pidCodeMC==6)
					hNchGenVsPtGenPID[0][4]->Fill(multTSgen,particle->Pt());
				else
					continue;
			}
			else if(TMath::Abs(DPhi-pi)<pi/3.0){

				if(pidCodeMC==6)
					hNchGenVsPtGenPID[1][4]->Fill(multTSgen,particle->Pt());
				else
					continue;
			}
			else{
				if(pidCodeMC==6)
					hNchGenVsPtGenPID[2][4]->Fill(multTSgen,particle->Pt());
				else
					continue;
			}

		} // Charge == 0.0

	}

	hPhiLeading->Fill(fRecLeadPhi);

	// Filling rec pT vs UE (for pT I use 2015 track cuts, UE uses TPC-only)
	for(int i=0; i < iTracks; i++){  

		AliESDtrack* esdtrack = static_cast<AliESDtrack*>(fESD->GetTrack(i));
		if(!esdtrack) continue;
		if(i==fRecLeadIn) continue;
		if(esdtrack->Charge() == 0.0 ) continue;
		if(TMath::Abs(esdtrack->Eta()) > fEtaCut) continue;
		if(esdtrack->GetTPCsignalN() < fNcl) continue;
		if(esdtrack->Pt() < fPtMin) continue;

		AliESDtrack* track = 0x0;
		if(!fSelectHybridTracks){
			if(!fTrackFilter->AcceptTrack(esdtrack)) { continue; }
			else{ track = esdtrack; }
		}else{
			track = SetHybridTrackCuts(esdtrack,kFALSE,kFALSE,kFALSE);
			if(!track) { continue; }
		}

		const int label = TMath::Abs(track->GetLabel());
		AliMCParticle *trackMC = (AliMCParticle*)fMC->GetTrack(label);

		short pidCodeMC = 0;
		pidCodeMC = GetPidCode(trackMC->PdgCode());

		double DPhi = DeltaPhi(track->Phi(), fRecLeadPhi);

		hNchGenVsPtRec[3][0]->Fill(multTSgen,track->Pt());
		hNchVsPtDataTPC[3][0]->Fill(multTSrec,track->Pt());

		if(pidCodeMC==1){
			hNchGenVsPtRec[3][1]->Fill(multTSgen,track->Pt());
			hNchVsPtDataTPC[3][1]->Fill(multTSrec,track->Pt());
		}
		if(pidCodeMC==2){
			hNchGenVsPtRec[3][2]->Fill(multTSgen,track->Pt());
			hNchVsPtDataTPC[3][2]->Fill(multTSrec,track->Pt());
		}
		if(pidCodeMC==3){
			hNchGenVsPtRec[3][3]->Fill(multTSgen,track->Pt());
			hNchVsPtDataTPC[3][3]->Fill(multTSrec,track->Pt());
		}

		if(TMath::Abs(DPhi)<pi/3.0){

			hNchGenVsPtRec[0][0]->Fill(multTSgen,track->Pt());
			hNchVsPtDataTPC[0][0]->Fill(multTSrec,track->Pt());

			if(pidCodeMC==1){
				hNchGenVsPtRec[0][1]->Fill(multTSgen,track->Pt());
				hNchVsPtDataTPC[0][1]->Fill(multTSrec,track->Pt());
			}
			else if(pidCodeMC==2){
				hNchGenVsPtRec[0][2]->Fill(multTSgen,track->Pt());
				hNchVsPtDataTPC[0][2]->Fill(multTSrec,track->Pt());
			}
			else if(pidCodeMC==3){
				hNchGenVsPtRec[0][3]->Fill(multTSgen,track->Pt());
				hNchVsPtDataTPC[0][3]->Fill(multTSrec,track->Pt());
			}
			else 
				continue;
		}

		else if(TMath::Abs(DPhi-pi)<pi/3.0){

			hNchGenVsPtRec[1][0]->Fill(multTSgen,track->Pt());
			hNchVsPtDataTPC[1][0]->Fill(multTSrec,track->Pt());

			if(pidCodeMC==1){
				hNchGenVsPtRec[1][1]->Fill(multTSgen,track->Pt());
				hNchVsPtDataTPC[1][1]->Fill(multTSrec,track->Pt());
			}
			else if(pidCodeMC==2){
				hNchGenVsPtRec[1][2]->Fill(multTSgen,track->Pt());
				hNchVsPtDataTPC[1][2]->Fill(multTSrec,track->Pt());
			}
			else if(pidCodeMC==3){
				hNchGenVsPtRec[1][3]->Fill(multTSgen,track->Pt());
				hNchVsPtDataTPC[1][3]->Fill(multTSrec,track->Pt());
			}
			else 
				continue;
		}

		else{
			hNchGenVsPtRec[2][0]->Fill(multTSgen,track->Pt());
			hNchVsPtDataTPC[2][0]->Fill(multTSrec,track->Pt());

			if(pidCodeMC==1){
				hNchGenVsPtRec[2][1]->Fill(multTSgen,track->Pt());
				hNchVsPtDataTPC[2][1]->Fill(multTSrec,track->Pt());
			}
			else if(pidCodeMC==2){
				hNchGenVsPtRec[2][2]->Fill(multTSgen,track->Pt());
				hNchVsPtDataTPC[2][2]->Fill(multTSrec,track->Pt());
			}
			else if(pidCodeMC==3){
				hNchGenVsPtRec[2][3]->Fill(multTSgen,track->Pt());
				hNchVsPtDataTPC[2][3]->Fill(multTSrec,track->Pt());
			}
			else 
				continue;
		}

		delete track;

	}
}
//_____________________________________________________________________________
void AliAnalysisTaskSpectraMC::GetDetectorResponse() {

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
			hPhiGen[0]->Fill(DPhi);
		}
		else if(TMath::Abs(DPhi-pi)<pi/3.0){
			hPhiGen[1]->Fill(DPhi);
		}
		else{
			multTSgen++;
			hPhiGen[2]->Fill(DPhi);
		}

		hPhiGen[3]->Fill(DPhi);
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
			track = SetHybridTrackCuts(esdtrack,kTRUE,kTRUE,kTRUE);
			if(!track) { continue; }
		}

		hPhiTotal->Fill(track->Eta(),track->Phi());
		double DPhi = DeltaPhi(track->Phi(), fRecLeadPhi);

		if(TMath::Abs(DPhi)<pi/3.0){
			hPhiRec[0]->Fill(DPhi);
		}
		else if(TMath::Abs(DPhi-pi)<pi/3.0){
			hPhiRec[1]->Fill(DPhi);
		}
		else{
			multTSrec++;
			hPhiRec[2]->Fill(DPhi);
		}

		hPhiRec[3]->Fill(DPhi);

		delete track;

	}

	hNchResponse->Fill(multTSrec,multTSgen);

}
//_____________________________________________________________________________
void AliAnalysisTaskSpectraMC::GetMCCorrections(){

	int iTracks(fESD->GetNumberOfTracks());          

	for (int i = 0; i < fMC->GetNumberOfTracks(); i++){

		AliMCParticle* particle = (AliMCParticle*)fMC->GetTrack(i);
		if (!particle) continue;
		if (TMath::Abs(particle->Eta()) > fEtaCut)continue;
		if (particle->Pt() < fPtMin)continue;
		short pidCodeMC = GetPidCode(particle->PdgCode());

		if(pidCodeMC==6) // Phi
			hPtGenIn[4]->Fill(particle->Pt());

		if (!fMC->IsPhysicalPrimary(i)) continue;
		if (particle->Charge() == 0.0) continue;

		hPtGenIn[0]->Fill(particle->Pt());
		if(particle->Charge() > 0.0)
			hPtGenPosIn[0]->Fill(particle->Pt());
		if(particle->Charge() < 0.0)
			hPtGenNegIn[0]->Fill(particle->Pt());

		if(pidCodeMC==1){ // Pion
			hPtGenIn[1]->Fill(particle->Pt());
			if(particle->Charge() > 0.0)
				hPtGenPosIn[1]->Fill(particle->Pt());
			if(particle->Charge() < 0.0)
				hPtGenNegIn[1]->Fill(particle->Pt());
		}
		else if(pidCodeMC==2){ // Kaon
			hPtGenIn[2]->Fill(particle->Pt());
			if(particle->Charge() > 0.0)
				hPtGenPosIn[2]->Fill(particle->Pt());
			if(particle->Charge() < 0.0)
				hPtGenNegIn[2]->Fill(particle->Pt());
		}
		else if(pidCodeMC==3){ // Proton
			hPtGenIn[3]->Fill(particle->Pt());
			if(particle->Charge() > 0.0)
				hPtGenPosIn[3]->Fill(particle->Pt());
			if(particle->Charge() < 0.0)
				hPtGenNegIn[3]->Fill(particle->Pt());
		}
		else{
			continue;
		}

	} // loop of true-MC particles

	for(int i = 0; i < iTracks; i++){                

		AliESDtrack* esdtrack = static_cast<AliESDtrack*>(fESD->GetTrack(i)); 
		if(!esdtrack) continue;        
		if(esdtrack->Charge() == 0.0 ) continue;
		if(TMath::Abs(esdtrack->Eta()) > fEtaCut) continue;        
		if(esdtrack->Pt() < fPtMin) continue;
		if(esdtrack->GetTPCsignalN() < fNcl) continue;

		AliESDtrack* track = 0x0;
		if(!fSelectHybridTracks){
			if(!fTrackFilter->AcceptTrack(esdtrack)) { continue; }
			else{ track = esdtrack; }
		}else{
			track = SetHybridTrackCuts(esdtrack,kFALSE,kFALSE,kFALSE);
			if(!track) { continue; }
		}

		const int label = TMath::Abs(track->GetLabel());
		AliMCParticle *trackMC = (AliMCParticle*)fMC->GetTrack(label);
		short pidCodeMC = 0;
		pidCodeMC = GetPidCode(trackMC->PdgCode());

		hPtRec->Fill(track->Pt());

		if(pidCodeMC==1)
			hPtRecPion->Fill(track->Pt());
		if(pidCodeMC==3)
			hPtRecProton->Fill(track->Pt());

		if(fMC->IsSecondaryFromWeakDecay(label) || fMC->IsSecondaryFromMaterial(label))
			hPtSec->Fill(track->Pt());		

		if( fMC->IsPhysicalPrimary(label) ){

			hPtPri->Fill(track->Pt());
			hPtRecIn[0]->Fill(track->Pt());
			if(track->Charge() > 0.0)
				hPtRecPosIn[0]->Fill(track->Pt());
			if(track->Charge() < 0.0)
				hPtRecNegIn[0]->Fill(track->Pt());

			if(TOFPID(track)){
				hPtRecInTOF[0]->Fill(track->Pt());
				if(track->Charge() > 0.0)
					hPtRecPosInTOF[0]->Fill(track->Pt());
				if(track->Charge() < 0.0)
					hPtRecNegInTOF[0]->Fill(track->Pt());

			}	// TOF

			// rTPC efficiency
			if(PhiCut(track->Pt(), track->Phi(), track->Charge(), Magf, fcutLow, fcutHigh)){
				hPtrTPCRecIn[0]->Fill(track->Pt());
			}

			if(pidCodeMC==1){
				hPtRecIn[1]->Fill(track->Pt());
				if(track->Charge() > 0.0)
					hPtRecPosIn[1]->Fill(track->Pt());
				if(track->Charge() < 0.0)
					hPtRecNegIn[1]->Fill(track->Pt());

				if(TOFPID(track)){
					hPtRecInTOF[1]->Fill(track->Pt());
					if(track->Charge() > 0.0)
						hPtRecPosInTOF[1]->Fill(track->Pt());
					if(track->Charge() < 0.0)
						hPtRecNegInTOF[1]->Fill(track->Pt());
				}	// TOF

				// rTPC efficiency
				if(PhiCut(track->Pt(), track->Phi(), track->Charge(), Magf, fcutLow, fcutHigh)){
					hPtrTPCRecIn[1]->Fill(track->Pt());
				}
			}

			else if(pidCodeMC==2){
				hPtRecIn[2]->Fill(track->Pt());
				if(track->Charge() > 0.0)
					hPtRecPosIn[2]->Fill(track->Pt());
				if(track->Charge() < 0.0)
					hPtRecNegIn[2]->Fill(track->Pt());

				if(TOFPID(track)){
					hPtRecInTOF[2]->Fill(track->Pt());
					if(track->Charge() > 0.0)
						hPtRecPosInTOF[2]->Fill(track->Pt());
					if(track->Charge() < 0.0)
						hPtRecNegInTOF[2]->Fill(track->Pt());
				}	// TOF

				// rTPC efficiency
				if(PhiCut(track->Pt(), track->Phi(), track->Charge(), Magf, fcutLow, fcutHigh)){
					hPtrTPCRecIn[2]->Fill(track->Pt());
				}
			}

			else if(pidCodeMC==3){
				hPtRecIn[3]->Fill(track->Pt());
				if(track->Charge() > 0.0)
					hPtRecPosIn[3]->Fill(track->Pt());
				if(track->Charge() < 0.0)
					hPtRecNegIn[3]->Fill(track->Pt());

				if(TOFPID(track)){
					hPtRecInTOF[3]->Fill(track->Pt());
					if(track->Charge() > 0.0)
						hPtRecPosInTOF[3]->Fill(track->Pt());
					if(track->Charge() < 0.0)
						hPtRecNegInTOF[3]->Fill(track->Pt());
				}	// TOF

				// rTPC efficiency
				if(PhiCut(track->Pt(), track->Phi(), track->Charge(), Magf, fcutLow, fcutHigh)){
					hPtrTPCRecIn[3]->Fill(track->Pt());
				}
			}

			else{
				continue;
			}

		} // Primary particles loop 

		delete track;

	} // Tracks loop

}
//_____________________________________________________________________________
void AliAnalysisTaskSpectraMC::GetMultiplicityDistributionsPhi(){

	int multTSgen=0;
	int multTSrec=0;
	const double pi = TMath::Pi();

	for (int i = 0; i < fMC->GetNumberOfTracks(); i++) {

		AliMCParticle* particle = (AliMCParticle*)fMC->GetTrack(i);
		if (!particle) continue;
		if (i==fGenLeadIn) continue;
		if (!fMC->IsPhysicalPrimary(i)) continue;
		if (particle->Charge() == 0.0) continue;
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
		if(esdtrack->Charge() == 0.0 ) continue;
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

	// Filling rec pT vs UE (for pT I use 2015 track cuts, UE uses TPC-only)
	for(int i=0; i < iTracks; i++){  

		AliESDtrack* esdtrack = static_cast<AliESDtrack*>(fESD->GetTrack(i)); 
		if(!esdtrack) continue;        
		if(esdtrack->Charge() == 0.0 ) continue;
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

		const int label = TMath::Abs(track->GetLabel());
		TParticle *trackMC = (TParticle*)fMC->Particle(label);

		int indexm = trackMC->GetFirstMother();
		if (indexm < 0) continue;

		TParticle *trackMC_mother = (TParticle*)fMC->Particle(indexm);

		int index_track1 = i;

		if(trackMC_mother->GetPdgCode() != 333)
			continue;

		for(int i2 = 0; i2 < iTracks; i2++){

			AliESDtrack* esdtrack = static_cast<AliESDtrack*>(fESD->GetTrack(i2)); 
			if(!esdtrack) continue;        
			if(esdtrack->Charge() == 0.0 ) continue;
			if(TMath::Abs(esdtrack->Eta()) > fEtaCut) continue;        
			if(esdtrack->Pt() < fPtMin) continue;
			if(index_track1==i2) continue; // skipping same track

			AliESDtrack* track2 = 0x0;
			if(!fSelectHybridTracks){
				if(!fTrackFilter->AcceptTrack(esdtrack)) { continue; } 
				else{ track2 = esdtrack; }
			}else{
				track2 = SetHybridTrackCuts(esdtrack,kFALSE,kFALSE,kFALSE);
				if(!track2) { continue; }
			}

			const int label2 = TMath::Abs(track2->GetLabel());
			TParticle *trackMC2 = (TParticle*)fMC->Particle(label2);

			int indexm2 = trackMC2->GetFirstMother();
			if (indexm2 < 0) continue;

			TParticle *trackMC_mother2 = (TParticle*)fMC->Particle(indexm2);

			if(trackMC_mother2->GetPdgCode() != 333) // mother not phi
				continue;

			if(TMath::Abs(trackMC->GetPdgCode()) != 321) // kaon
				continue;

			if(TMath::Abs(trackMC2->GetPdgCode()) != 321) // kaon
				continue;

			if(track->Charge() == track2->Charge()) // same sign
				continue;

			if(TMath::Abs(indexm) != TMath::Abs(indexm2)) // different mother
				continue;

			double px1 = track->Px();
			double py1 = track->Py();
			double px2 = track2->Px();
			double py2 = track2->Py();
			double pxm = px1 + px2;
			double pym = py1 + py2;
			double ptm = TMath::Sqrt(pxm*pxm + pym*pym);
			double phim = TMath::ATan2(pym,pxm);
			double DPhi = DeltaPhi(phim, fRecLeadPhi);

			hPhiResPhi->Fill(phim);

			hNchGenVsPtRec[3][4]->Fill(multTSgen,ptm);
			hNchVsPtDataTPC[3][4]->Fill(multTSrec,ptm);

			if(TMath::Abs(DPhi)<pi/3.0){
				hNchGenVsPtRec[0][4]->Fill(multTSgen,ptm);
				hNchVsPtDataTPC[0][4]->Fill(multTSrec,ptm);
			}
			else if(TMath::Abs(DPhi-pi)<pi/3.0){
				hNchGenVsPtRec[1][4]->Fill(multTSgen,ptm);
				hNchVsPtDataTPC[1][4]->Fill(multTSrec,ptm);
			}
			else{
				hNchGenVsPtRec[2][4]->Fill(multTSgen,ptm);
				hNchVsPtDataTPC[2][4]->Fill(multTSrec,ptm);
			}

			delete track2;

		} // Track2 loop

		delete track;

	} // Track1 loop

}

//_____________________________________________________________________________
void AliAnalysisTaskSpectraMC::GetMCCorrectionsPhi()
{

	int iTracks(fESD->GetNumberOfTracks());          

	for(int i = 0; i < iTracks; i++){                

		AliESDtrack* esdtrack = static_cast<AliESDtrack*>(fESD->GetTrack(i)); 
		if(!esdtrack) continue;        
		if(esdtrack->Charge() == 0.0 ) continue;
		if(TMath::Abs(esdtrack->Eta()) > fEtaCut) continue;        
		if(esdtrack->Pt() < fPtMin) continue;

		AliESDtrack* track1 = 0x0;
		if(!fSelectHybridTracks){
			if(!fTrackFilter->AcceptTrack(esdtrack)) { continue; }
			else{ track1 = esdtrack; }
		}else{
			track1 = SetHybridTrackCuts(esdtrack,kFALSE,kFALSE,kFALSE);
			if(!track1) { continue; }
		}

		const int label1 = TMath::Abs(track1->GetLabel());
		TParticle *trackMC1 = (TParticle*)fMC->Particle(label1);

		int indexm = trackMC1->GetFirstMother();
		if (indexm < 0) continue;

		TParticle *trackMC_mother = (TParticle*)fMC->Particle(indexm);

		int index_track1 = i;

		if(trackMC_mother->GetPdgCode() != 333)
			continue;

		for(int i2 = 0; i2 < iTracks; i2++){

			AliESDtrack* esdtrack = static_cast<AliESDtrack*>(fESD->GetTrack(i2)); 
			if(!esdtrack) continue;        
			if(esdtrack->Charge() == 0.0 ) continue;
			if(TMath::Abs(esdtrack->Eta()) > fEtaCut) continue;        
			if(esdtrack->Pt() < fPtMin) continue;
			if(index_track1==i2) continue; // skipping same track

			AliESDtrack* track2 = 0x0;
			if(!fSelectHybridTracks){
				if(!fTrackFilter->AcceptTrack(esdtrack)) { continue; }
				else{ track2 = esdtrack; }
			}else{
				track2 = SetHybridTrackCuts(esdtrack,kFALSE,kFALSE,kFALSE);
				if(!track2) { continue; }
			}

			const int label2 = TMath::Abs(track2->GetLabel());
			TParticle *trackMC2 = (TParticle*)fMC->Particle(label2);

			int indexm2 = trackMC2->GetFirstMother();
			if (indexm2 < 0) continue;

			TParticle *trackMC_mother2 = (TParticle*)fMC->Particle(indexm2);

			if(trackMC_mother2->GetPdgCode() != 333) // mother not phi
				continue;

			if(TMath::Abs(trackMC1->GetPdgCode()) != 321) // kaon
				continue;

			if(TMath::Abs(trackMC2->GetPdgCode()) != 321) // kaon
				continue;

			if(track1->Charge() == track2->Charge()) // same sign
				continue;

			if(TMath::Abs(indexm) != TMath::Abs(indexm2)) // different mother
				continue;

			double p1 = track1->P();
			double px1 = track1->Px();
			double py1 = track1->Py();
			double pz1 = track1->Pz();
			double E1  = TMath::Sqrt(p1*p1 + 0.497648*0.497648);

			double p2  = track2->P();
			double px2 = track2->Px();
			double py2 = track2->Py();
			double pz2 = track2->Pz();
			double E2  = TMath::Sqrt(p2*p2 + 0.497648*0.497648);

			double pxm = px1 + px2;
			double pym = py1 + py2;
			double pzm = pz1 + pz2;
			double Em  = E1 + E2;
			double ptm = TMath::Sqrt(pxm*pxm + pym*pym);
			double invMm = TMath::Sqrt(Em*Em - pxm*pxm - pym*pym - pzm*pzm); // Mother Mass

			hInvMassPhi->Fill(invMm);
			hPtRecIn[4]->Fill(ptm);

			delete track2;

		} // Track2 loop

		delete track1;

	} // Tracks loop

}
//_____________________________________________________________________________
double AliAnalysisTaskSpectraMC::DeltaPhi(Double_t phi, Double_t Lphi,
		Double_t rangeMin, Double_t rangeMax)
{

	Double_t dphi = -999;
	Double_t pi = TMath::Pi();
	if(Lphi > 2*pi || Lphi < 0)cout << "Lphi :: " << Lphi << endl;
	if(phi  > 2*pi || phi < 0)cout << "phi = " << phi << endl;

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
short AliAnalysisTaskSpectraMC::GetPidCode(Int_t pdgCode) const
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

//________________________________________________________________________

bool AliAnalysisTaskSpectraMC::selectVertex2015pp(AliESDEvent *esd,
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
bool AliAnalysisTaskSpectraMC::IsGoodSPDvertexRes(const AliESDVertex* spdVertex)
{

	if( !spdVertex ) return kFALSE;
	if( spdVertex->IsFromVertexerZ() && !(spdVertex->GetDispersion()<0.04 && spdVertex->GetZRes()<0.25) ) return kFALSE;
	return kTRUE;
}
//________________________________________________________________________
bool AliAnalysisTaskSpectraMC::IsGoodZvertexPos(AliESDEvent *esd)
{

	if( !esd ) return kFALSE;
	//Cut on the vertex z position
	const AliESDVertex * vertex = esd->GetPrimaryVertex();
	if (TMath::Abs(vertex->GetZ())>10) return kFALSE;
	return kTRUE;
}
//________________________________________________________________________
bool AliAnalysisTaskSpectraMC::PhiCut(const double& pt, double phi, const double& q, const float& mag, TF1* phiCutLow, TF1* phiCutHigh)
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

	//    hPhi[4]->Fill(pt, phi);

	return kTRUE;
}
//________________________________________________________________________
float AliAnalysisTaskSpectraMC::GetMaxDCApTDep( TF1 *fMaxDCAxy, Double_t ptI){

	double maxDCAxy = 10;
	maxDCAxy = fMaxDCAxy->Eval(ptI);
	return maxDCAxy;

}
//________________________________________________________________________
/*void AliAnalysisTaskSpectraMC::SetTrackCuts(AliAnalysisFilter* fTrackFilter){

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
AliESDtrack* AliAnalysisTaskSpectraMC::SetHybridTrackCuts(AliESDtrack *esdtrack, const bool fillPhiStand, const bool fillPhHyb1, const bool fillPhHyb2){

	// 
	// 	Get the Hybrid Tracks 
	// 	

	AliESDtrack *newTrack = 0x0;

	if(fTrackFilter->AcceptTrack(esdtrack))
	{
		newTrack = new AliESDtrack(*esdtrack);
		if(fillPhiStand) hPhiStandard->Fill(newTrack->Eta(),newTrack->Phi());
		////			newTrack->SetTRDQuality(0);
	}
	else if(fHybridTrackCuts1->AcceptTrack(esdtrack))
	{
		if(esdtrack->GetConstrainedParam())
		{
			newTrack = new AliESDtrack(*esdtrack);
			const AliExternalTrackParam* constrainParam = esdtrack->GetConstrainedParam();
			newTrack->Set(constrainParam->GetX(),constrainParam->GetAlpha(),constrainParam->GetParameter(),constrainParam->GetCovariance());
			////				newTrack->SetTRDQuality(1);
			if(fillPhHyb1) hPhiHybrid1->Fill(newTrack->Eta(),newTrack->Phi());
		}
		else
			return 0x0;
	}
	/*else if(fHybridTrackCuts2->AcceptTrack(esdtrack))
	  {
	  if(esdtrack->GetConstrainedParam())
	  {
	  newTrack = new AliESDtrack(*esdtrack);
	  const AliExternalTrackParam* constrainParam = esdtrack->GetConstrainedParam();
	  newTrack->Set(constrainParam->GetX(),constrainParam->GetAlpha(),constrainParam->GetParameter(),constrainParam->GetCovariance());
	/////				newTrack->SetTRDQuality(2);
	if(fillPhHyb2) hPhiHybrid2->Fill(newTrack->Eta(),newTrack->Phi());
	}
	else
	return 0x0;
	}*/
	  else
	  {
		  return 0x0;
	  }

	  return newTrack;

}
//________________________________________________________________________
double AliAnalysisTaskSpectraMC::EtaCalibration( const int &indx, const double &eta){
	//    h        i         j        l      k        o          p
	const Double_t aPos[nRt+2]      = {49.9044 ,50.0841  ,49.8419 ,49.9799 ,49.9659 ,50.0535 ,50.0649};
	const Double_t bPos[nRt+2]      = {4.05075 ,-0.743724,10.3952 ,2.99619 ,2.91366 ,-2.87404,-4.3589};
	const Double_t cPos[nRt+2]      = {-58.1027,-13.8508 ,-151.227,-45.718 ,-45.5994,31.159  ,66.268 };
	const Double_t dPos[nRt+2]      = {342.297 ,141.269  ,900.83  ,290.013 ,290.042 ,-151.257,-435.715};
	const Double_t ePos[nRt+2]      = {-1098.19,-567.054 ,-2833.25,-1018.42,-1014.49,282.703 ,1325.24};
	const Double_t fPos[nRt+2]      = {1944.32 ,1105.63  ,4884.59 ,1948.68 ,1931.84 ,-98.8756,-2030.25};
	const Double_t gPos[nRt+2]      = {-1749.14,-1022.19 ,-4318.99,-1864.06,-1839.36,-230.114,1542.09};
	const Double_t hPos[nRt+2]      = {617.929 ,355.158  ,1521.66 ,692.752 ,680.421 ,172.854 ,-468.577};

	const Double_t aNeg[nRt+2]      = {49.9261,50.0561,49.9583,50.078 ,50.046 ,49.9496,50.1258};
	const Double_t bNeg[nRt+2]      = {2.92422,4.68965,3.38038,6.67199,6.79992,2.45301,12.7977};
	const Double_t cNeg[nRt+2]      = {61.6661,65.891 ,53.1256,103.662,109.86 ,53.654 ,190.076};
	const Double_t dNeg[nRt+2]      = {421.545,394.542,314.489,611.034,668.241,363.689,1144.11};
	const Double_t eNeg[nRt+2]      = {1283.04,1100.56,825.296,1695.63,1916.44,1115.13,3411.98};
	const Double_t fNeg[nRt+2]      = {1944.85,1516.21,1021.01,2395.88,2815.04,1762.18,5402.99};
	const Double_t gNeg[nRt+2]      = {1442.98,989.24 ,548.44 ,1669.22,2057.21,1421.46,4379.16};
	const Double_t hNeg[nRt+2]      = {419.491,238.333,84.7945,455.362,595.391,469.45 ,1436.76};


	for(Int_t i=0; i<8; ++i)
		fEtaCalibration->SetParameter(i,0);

	if(eta<0){
		fEtaCalibration->SetParameter(0,aNeg[indx]);
		fEtaCalibration->SetParameter(1,bNeg[indx]);
		fEtaCalibration->SetParameter(2,cNeg[indx]);
		fEtaCalibration->SetParameter(3,dNeg[indx]);
		fEtaCalibration->SetParameter(4,eNeg[indx]);
		fEtaCalibration->SetParameter(5,fNeg[indx]);
		fEtaCalibration->SetParameter(6,gNeg[indx]);
		fEtaCalibration->SetParameter(7,hNeg[indx]);
	}
	else{
		fEtaCalibration->SetParameter(0,aPos[indx]);
		fEtaCalibration->SetParameter(1,bPos[indx]);
		fEtaCalibration->SetParameter(2,cPos[indx]);
		fEtaCalibration->SetParameter(3,dPos[indx]);
		fEtaCalibration->SetParameter(4,ePos[indx]);
		fEtaCalibration->SetParameter(5,fPos[indx]);
		fEtaCalibration->SetParameter(6,gPos[indx]);
		fEtaCalibration->SetParameter(7,hPos[indx]);
	}

	return fEtaCalibration->Eval(eta);

}
//________________________________________________________________________
double AliAnalysisTaskSpectraMC::EtaCalibrationEl(const int &indx, const double &eta){

	const Double_t aPosEl[nRt+2]    = {79.8647 ,79.6737 ,80.3915 ,80.1263 ,79.9957 ,79.6537 ,80.6434 };
	const Double_t bPosEl[nRt+2]    = {6.50512 ,16.0745 ,9.53925 ,5.28525 ,7.03079 ,15.0221 ,0.40293 };
	const Double_t cPosEl[nRt+2]    = {-35.9277,-80.5639,-69.3773,-32.7731,-42.9098,-83.6391,-21.8162};
	const Double_t dPosEl[nRt+2]    = {73.1535 ,148.866 ,143.956 ,68.4524 ,88.7057 ,168.5   ,61.9147 };
	const Double_t ePosEl[nRt+2]    = {-47.1041,-90.3376,-89.5518,-44.1566,-56.6554,-107.999,-44.6593};

	const Double_t aNegEl[nRt+2]    = {79.6366 ,80.0767 ,79.6157 ,79.8351 ,79.7387 ,79.3638 ,79.9111 };
	const Double_t bNegEl[nRt+2]    = {-11.3437,-2.51009,-16.2468,-8.46921,-8.60021,-17.1977,-1.66066};
	const Double_t cNegEl[nRt+2]    = {-65.1353,-23.6188,-92.0783,-44.5947,-44.1718,-82.7998,-6.96109};
	const Double_t dNegEl[nRt+2]    = {-134.447,-65.5053,-180.753,-86.2242,-84.4984,-143.394,-16.0465};
	const Double_t eNegEl[nRt+2]    = {-87.7848,-51.1463,-112.997,-53.6285,-51.945 ,-81.3439,-10.3587};


	for(Int_t i=0; i<5; ++i)
		fEtaCalibrationEl->SetParameter(i,0);

	if(eta<0){
		fEtaCalibrationEl->SetParameter(0,aNegEl[indx]);
		fEtaCalibrationEl->SetParameter(1,bNegEl[indx]);
		fEtaCalibrationEl->SetParameter(2,cNegEl[indx]);
		fEtaCalibrationEl->SetParameter(3,dNegEl[indx]);
		fEtaCalibrationEl->SetParameter(4,eNegEl[indx]);
	}
	else{
		fEtaCalibrationEl->SetParameter(0,aPosEl[indx]);
		fEtaCalibrationEl->SetParameter(1,bPosEl[indx]);
		fEtaCalibrationEl->SetParameter(2,cPosEl[indx]);
		fEtaCalibrationEl->SetParameter(3,dPosEl[indx]);
		fEtaCalibrationEl->SetParameter(4,ePosEl[indx]);
	}

	return fEtaCalibrationEl->Eval(eta);

}
//________________________________________________________________________
bool AliAnalysisTaskSpectraMC::TOFPID(AliESDtrack * track)
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


