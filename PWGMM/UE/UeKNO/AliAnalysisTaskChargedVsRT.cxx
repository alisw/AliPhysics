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
 * provided "as is" without express or implied warranty.     
 *                                                                        *
 * Authors:       Luz Tiscareño (luz.elena.tiscareno.montoya@cern.ch)     *
 *                Antonio Ortiz (antonio.ortiz@nucleares.unam.mx)         * 
 *                                                                        *
 **************************************************************************/

/* AliAnaysisTaskChargedVsRT source code
 * The analysis task produce all the histos needed for MC closure test studies
 * Results include only the KNO properties
 */

class TTree;

class AliPPVsMultUtils;
class AliESDtrackCuts;


#include <Riostream.h>
#include "TChain.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TProfile.h"
#include "THnSparse.h"
#include "TVector3.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TLegend.h"
#include "TList.h"
#include "AliLog.h"
#include "AliVEvent.h"
#include "AliVVertex.h"
#include "AliVTrack.h"
//#include "AliV0vertexer.h"
#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliESDEvent.h"
#include "AliESDVZERO.h"
#include "AliESDInputHandler.h"
#include "AliMultiplicity.h"
#include "AliPPVsMultUtils.h"
#include "AliESDtrack.h"
#include "AliESDtrackCuts.h"
#include "AliAnalysisTaskESDfilter.h"
#include "AliMCEventHandler.h"
#include "AliInputEventHandler.h"
#include "AliMCEvent.h"
#include "AliMCParticle.h"
#include "AliStack.h"
#include "TParticle.h"
#include <AliHeader.h>
#include "AliOADBContainer.h"
#include "AliOADBMultSelection.h"
#include "AliMultSelection.h"
#include "AliMultEstimator.h"
#include "AliMultInput.h"
#include "AliMultVariable.h"
#include "AliCentrality.h"
#include "AliMultiplicity.h"
#include "AliESDUtils.h"
#include "AliGenEventHeader.h"
#include "AliGenCocktailEventHeader.h"
//#include "AliGenPythiaEventHeader.h"
#include "AliEventCuts.h"
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisUtils.h"
#include "AliPPVsMultUtils.h"
#include <AliESDVertex.h>
#include <AliMultiplicity.h>
#include <TTree.h>
#include <TMath.h>
#include <TRandom.h>
#include <TDirectory.h>
#include <TBits.h>
#include <AliAnalysisFilter.h>

using std::cout;
using std::endl;

#include "AliAnalysisTaskChargedVsRT.h"


const Char_t * NameReg_3[3]={"NS","AS","TS"};


const Int_t ptNbins = 36;
Double_t ptbins1_3[ptNbins+1] = {
	0.0,  0.1,  0.15,  0.2,  0.25,  0.3,  0.35,  0.4,  0.45,  0.5,  0.6,  0.7,  0.8,  0.9,  1.0, 1.25, 1.5,  2.0,  2.5,  3.0,  3.5,  4.0,  4.5, 5.0, 6.0, 7.0,  8.0,  9.0,  10.0,  12.0,  14.0,  16.0,  18.0,  20.0,  30.0,  40.0,  50.0};

const int nBinsDCAxy = 121;
double binsDCAxy3[nBinsDCAxy+1] = {-3.025,-2.975,-2.925,-2.875,-2.825,-2.775,-2.725,-2.675,-2.625,-2.575,-2.525,-2.475,-2.425,-2.375,-2.325,-2.275,-2.225, -2.175,-2.125,-2.075,-2.025,-1.975,-1.925,-1.875,-1.825,-1.775,-1.725,-1.675,-1.625,-1.575,-1.525,-1.475,-1.425,-1.375,-1.325,-1.275,-1.225,-1.175,-1.125,-1.075,-1.025,-0.975,-0.925,-0.875,-0.825,-0.775,-0.725,-0.675,-0.625,-0.575,-0.525,-0.475,-0.425,-0.375,-0.325,-0.275,-0.225,-0.175,-0.125,-0.075,-0.025,0.025,0.075,0.125,0.175,0.225,0.275,0.325,0.375,0.425,0.475,0.525,0.575,0.625,0.675,0.725,0.775,0.825,0.875,0.925,0.975,1.025,1.075,1.125,1.175,1.225,1.275,1.325,1.375,1.425,1.475,1.525,1.575,1.625,1.675,1.725,1.775,1.825,1.875,1.925,1.975,2.025,2.075,2.125,2.175,2.225,2.275,2.325,2.375,2.425,2.475,2.525,2.575,2.625,2.675,2.725,2.775,2.825,2.875,2.925,2.975,3.025};


const Double_t pi = 3.1415926535897932384626433832795028841971693993751058209749445;
Float_t MultV0M3, MultRef3;
class AliAnalysisTaskChargedVsRT;    // your analysis class

using namespace std;            // std namespace: so you can do things like 'cout' etc

ClassImp(AliAnalysisTaskChargedVsRT) // classimp: necessary for root

	AliAnalysisTaskChargedVsRT::AliAnalysisTaskChargedVsRT() : AliAnalysisTaskSE(),
	fESD(0),
	fEventCuts(0x0),
	fMCStack(0),
	fMC(0),
	fUseMC(kFALSE),
	fIsMCclosure(kFALSE),
	fIsHybAna(kFALSE),
	fnRecHy(-1),
	fnRecHyWoDCA(-1),
	fnGen(-1),
	fTPCclustersVar1(kFALSE),
	fTPCclustersVar2(kFALSE),
	fNcrVar1(kFALSE),
	fNcrVar2(kFALSE),
	fGeoTPCVar1(kFALSE),
	fGeoTPCVar2(kFALSE),
	fGeoTPCVar3(kFALSE),
	fGeoTPCVar4(kFALSE),
	fChisqTPCVar1(kFALSE),
	fChisqTPCVar2(kFALSE),
	fChisqITSVar1(kFALSE),
	fChisqITSVar2(kFALSE),
	fChisqITSmTPCVar1(kFALSE),
	fChisqITSmTPCVar2(kFALSE),
	fDcazVar1(kFALSE),
	fDcazVar2(kFALSE),
	fSPDreqVar1(kFALSE),
	fPIDResponse(0x0),
	fTrackFilterTPC(0x0),
	fTrackFilter2015(0x0),
	fTrackFilter2015woDCA(0x0),
	fTrackFilterHybrid0(0x0),
	fTrackFilterHybrid1(0x0),
	fTrackFilterHybrid2(0x0),
	fTrackFilterHybrid0woDCA(0x0),
	fTrackFilterHybrid1woDCA(0x0),
	fTrackFilterHybrid2woDCA(0x0),
	fOutputList(0),
	fEtaCut(0.8),
	fPtMin(0.5),
	fLeadPtCutMin(5.0),
	fLeadPtCutMax(40.0),
	fGenLeadPhi(0),
	fGenLeadPt(0),
	fGenLeadIn(0),
	fRecLeadPhi(0),
	fRecLeadPt(0),
	fRecLeadIn(0),
	ftrackmult08(0),
	fv0mpercentile(0),
	fdcaxy(-999),
	fdcaz(-999),
	fMultSelection(0x0),
	fMultSelectionbefvtx(0x0),
	hNchTSGen(0),
	hNchTSGenTest(0),
	hNchGen(0),
	hNchGenTest(0),
	hNchTSRec(0),
	hNchTSRecTest(0),
	hNchData(0),
	hNchTSData(0),
	//Only with Hybrid Track Cuts
	hPhiTotal(0), //Sum of all the contributions
	hPhiStandard(0), //Distribution of phi without corrections -w/ SPD & ITS-
	hPhiHybrid1(0), // Correction of phi distribution -w/o SPD & w/ ITS-
	hPhiHybrid2(0), // Correction of phi distribution -no SPD req. & w/o ITS-
	//
	hNchResponse(0),
	hNchRec(0),
	hNchRecTest(0),
	hPtInPrim(0),
	hPtInPrim_pion(0),
	hPtInPrim_kaon(0),
	hPtInPrim_proton(0),
	hPtInPrim_sigmap(0),
	hPtInPrim_sigmam(0),
	hPtInPrim_omega(0),
	hPtInPrim_xi(0),
	hPtInPrim_rest(0),
	hPtOut(0),
	hPtOutPrim(0),
	hPtOutPrim_pion(0),
	hPtOutPrim_kaon(0),
	hPtOutPrim_proton(0),
	hPtOutPrim_sigmap(0),
	hPtOutPrim_sigmam(0),
	hPtOutPrim_omega(0),
	hPtOutPrim_xi(0),
	hPtOutPrim_rest(0),
	hPtOutSec(0),
	hCounter(0),
	hPTVsDCAData(0),
	hptvsdcaPrim(0),
	hptvsdcaDecs(0),
	hptvsdcaMatl(0),
	hptvsdcaAll(0)
{

	for(Int_t i=0;i<3;++i){
		hPtVsUEGenTest[i]=0;
		hPtVsUERecTest[i]=0;
		hPtVsUEData[i]=0;
		hPtVsNchGenTest[i]=0;
		hPtVsNchRecTest[i]=0;
		hPtVsNchData[i]=0;
		hPhiGen[i]=0;
		hPhiRec[i]=0;
	}


	// default constructor, don't allocate memory here!  this is used by root for IO purposes, it needs to remain empty
}
//_____________________________________________________________________________
AliAnalysisTaskChargedVsRT::AliAnalysisTaskChargedVsRT(const char* name) : AliAnalysisTaskSE(name),
	fESD(0),
	fEventCuts(0x0),
	fMCStack(0),
	fMC(0),
	fUseMC(kFALSE),
	fIsMCclosure(kFALSE),
	fIsHybAna(kFALSE),
	fnRecHy(-1),
	fnRecHyWoDCA(-1),
	fnGen(-1),
	fTPCclustersVar1(kFALSE),
	fTPCclustersVar2(kFALSE),
	fNcrVar1(kFALSE),
	fNcrVar2(kFALSE),
	fGeoTPCVar1(kFALSE),
	fGeoTPCVar2(kFALSE),
	fGeoTPCVar3(kFALSE),
	fGeoTPCVar4(kFALSE),
	fChisqTPCVar1(kFALSE),
	fChisqTPCVar2(kFALSE),
	fChisqITSVar1(kFALSE),
	fChisqITSVar2(kFALSE),
	fChisqITSmTPCVar1(kFALSE),
	fChisqITSmTPCVar2(kFALSE),
	fDcazVar1(kFALSE),
	fDcazVar2(kFALSE),
	fSPDreqVar1(kFALSE),
	fPIDResponse(0x0),
	fTrackFilterTPC(0x0),
	fTrackFilter2015(0x0),
	fTrackFilter2015woDCA(0x0),
	fTrackFilterHybrid0(0x0),
	fTrackFilterHybrid1(0x0),
	fTrackFilterHybrid2(0x0),
	fTrackFilterHybrid0woDCA(0x0),
	fTrackFilterHybrid1woDCA(0x0),
	fTrackFilterHybrid2woDCA(0x0),
	fOutputList(0),
	fEtaCut(0.8),
	fPtMin(0.5),
	fLeadPtCutMin(5.0),
	fLeadPtCutMax(40.0),
	fGenLeadPhi(0),
	fGenLeadPt(0),
	fGenLeadIn(0),
	fRecLeadPhi(0),
	fRecLeadPt(0),
	fRecLeadIn(0),
	ftrackmult08(0),
	fv0mpercentile(0),
	fdcaxy(-999),
	fdcaz(-999),
	fMultSelection(0x0),
	fMultSelectionbefvtx(0x0),
	hNchTSGen(0),
	hNchTSGenTest(0),
	hNchGen(0),
	hNchGenTest(0),
	hNchTSRec(0),
	hNchTSRecTest(0),
	hNchData(0),
	hNchTSData(0),
	//Only with Hybrid Track Cuts
	hPhiTotal(0), //Sum of all the contributions
	hPhiStandard(0), //Distribution of phi without corrections -w/ SPD & ITS-
	hPhiHybrid1(0), // Correction of phi distribution -w/o SPD & w/ ITS-
	hPhiHybrid2(0), // Correction of phi distribution -no SPD req. & w/o ITS-
	//
	hNchResponse(0),
	hNchRec(0),
	hNchRecTest(0),
	hPtInPrim(0),
	hPtInPrim_pion(0),
	hPtInPrim_kaon(0),
	hPtInPrim_proton(0),
	hPtInPrim_sigmap(0),
	hPtInPrim_sigmam(0),
	hPtInPrim_omega(0),
	hPtInPrim_xi(0),
	hPtInPrim_rest(0),
	hPtOut(0),
	hPtOutPrim(0),
	hPtOutPrim_pion(0),
	hPtOutPrim_kaon(0),
	hPtOutPrim_proton(0),
	hPtOutPrim_sigmap(0),
	hPtOutPrim_sigmam(0),
	hPtOutPrim_omega(0),
	hPtOutPrim_xi(0),
	hPtOutPrim_rest(0),
	hPtOutSec(0),
	hCounter(0),
	hPTVsDCAData(0),
	hptvsdcaPrim(0),
	hptvsdcaDecs(0),
	hptvsdcaMatl(0),
	hptvsdcaAll(0)

{

	for(Int_t i=0;i<3;++i){
		hPtVsUEGenTest[i]=0;
		hPtVsUERecTest[i]=0;
		hPtVsUEData[i]=0;
		hPtVsNchGenTest[i]=0;
		hPtVsNchRecTest[i]=0;
		hPtVsNchData[i]=0;
		hPhiGen[i]=0;
		hPhiRec[i]=0;
	}

	DefineInput(0, TChain::Class());    // define the input of the analysis: in this case you take a 'chain' of events
	// this chain is created by the analysis manager, so no need to worry about it, does its work automatically
	DefineOutput(1, TList::Class());    // define the ouptut of the analysis: in this case it's a list of histograms

}
//_____________________________________________________________________________
AliAnalysisTaskChargedVsRT::~AliAnalysisTaskChargedVsRT()
{
	// destructor
	if(fOutputList) {
		delete fOutputList;     // at the end of your task, it is deleted from memory by calling this function
		fOutputList = 0x0;
	}

}
//_____________________________________________________________________________
void AliAnalysisTaskChargedVsRT::UserCreateOutputObjects()
{

	
	AliAnalysisManager *man = AliAnalysisManager::GetAnalysisManager();
	if(man){
		AliInputEventHandler* inputHandler = (AliInputEventHandler*)(man->GetInputEventHandler());
		if(inputHandler)fPIDResponse = inputHandler->GetPIDResponse();
	}

	//Int_t  fTrackCutsNchT, fTrackCutsPtSpectra, fTrackCutsPtLeading
	// 0: TPConlyTrackCuts, 1: StandardTPCITS2015, 2: Hybrid

	// Define TPConlyTrackCuts
	fTrackFilterTPC = new AliAnalysisFilter("fTrackFilterTPC");
	AliESDtrackCuts * fCutsTPC = AliESDtrackCuts::GetStandardTPCOnlyTrackCuts();
	fCutsTPC->SetRequireTPCRefit(kTRUE);
	fCutsTPC->SetRequireITSRefit(kTRUE);
	fCutsTPC->SetEtaRange(-0.8,0.8);
	fTrackFilterTPC->AddCuts(fCutsTPC);

	// Define TPCITS2015TrackCuts
	fTrackFilter2015 = new AliAnalysisFilter("trackFilter2015");
	fTrackFilter2015woDCA = new AliAnalysisFilter("trackFilter2015woDCA");// wo DCA cut
	AliESDtrackCuts * fCuts2015 = new AliESDtrackCuts();
	fCuts2015->SetMaxFractionSharedTPCClusters(0.4);//
	fCuts2015->SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);//
	fCuts2015->SetCutGeoNcrNcl(3., 130., 1.5, 0.85, 0.7);//
	fCuts2015->SetMaxChi2PerClusterTPC(4);//
	fCuts2015->SetAcceptKinkDaughters(kFALSE);//
	fCuts2015->SetRequireTPCRefit(kTRUE);//
	fCuts2015->SetRequireITSRefit(kTRUE);//
	fCuts2015->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kAny);//
	fCuts2015->SetMaxDCAToVertexZ(2);//
	fCuts2015->SetDCAToVertex2D(kFALSE);//
	fCuts2015->SetRequireSigmaToVertex(kFALSE);//
	fCuts2015->SetMaxChi2PerClusterITS(36);//
	fCuts2015->SetEtaRange(-0.8,0.8);
	// In case of variations for syst. uncertainties
	if (fTPCclustersVar1) {fCuts2015->SetMaxFractionSharedTPCClusters(0.2);}
	else if (fTPCclustersVar2) {fCuts2015->SetMaxFractionSharedTPCClusters(1.);}
	else {fCuts2015->SetMaxFractionSharedTPCClusters(0.4);}// Default

	if (fNcrVar1) {fCuts2015->SetMinRatioCrossedRowsOverFindableClustersTPC(0.7);}//
	else if (fNcrVar2) {fCuts2015->SetMinRatioCrossedRowsOverFindableClustersTPC(0.9);}//
	else {fCuts2015->SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);}// Default

	if (fGeoTPCVar1) {fCuts2015->SetCutGeoNcrNcl(2., 130., 1.5, 0.85, 0.7);}//
	else if (fGeoTPCVar2) {fCuts2015->SetCutGeoNcrNcl(4., 130., 1.5, 0.85, 0.7);}//
	else if (fGeoTPCVar3) {fCuts2015->SetCutGeoNcrNcl(3., 120., 1.5, 0.85, 0.7);}//
	else if (fGeoTPCVar4) {fCuts2015->SetCutGeoNcrNcl(3., 140., 1.5, 0.85, 0.7);}//
	else {fCuts2015->SetCutGeoNcrNcl(3., 130., 1.5, 0.85, 0.7);}// Default

	if (fChisqTPCVar1) {fCuts2015->SetMaxChi2PerClusterTPC(3);}
	else if (fChisqTPCVar2) {fCuts2015->SetMaxChi2PerClusterTPC(5);}
	else {fCuts2015->SetMaxChi2PerClusterTPC(4);}// Default

	if (!fSPDreqVar1) {fCuts2015->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kAny);}// Default

	if (fDcazVar1) {fCuts2015->SetMaxDCAToVertexZ(1);} // DCAz = 1 cm
	else if (fDcazVar2) {fCuts2015->SetMaxDCAToVertexZ(5);} // DCAz = 5 cm
	else {fCuts2015->SetMaxDCAToVertexZ(2);}// Default

	if (fChisqITSVar1) {fCuts2015->SetMaxChi2PerClusterITS(25);}//
	else if (fChisqITSVar2) {fCuts2015->SetMaxChi2PerClusterITS(49);}//
	else {fCuts2015->SetMaxChi2PerClusterITS(36);}// Default

	fTrackFilter2015woDCA->AddCuts(fCuts2015);

	fCuts2015->SetMaxDCAToVertexXYPtDep("0.0182+0.0350/pt^1.01");//
	fCuts2015->SetMaxChi2TPCConstrainedGlobal(36);
	if (fChisqITSmTPCVar1) {fCuts2015->SetMaxChi2TPCConstrainedGlobal(25);}//
	else if (fChisqITSmTPCVar2) {fCuts2015->SetMaxChi2TPCConstrainedGlobal(49);}//

	fTrackFilter2015->AddCuts(fCuts2015);

	// Define Hybrid 0 (global, 2011 track cuts)
	fTrackFilterHybrid0 = new AliAnalysisFilter("trackFilterHybrid0");
	fTrackFilterHybrid0woDCA = new AliAnalysisFilter("trackFilterHybrid0woDCA");
	AliESDtrackCuts * fCutsHybrid0 = new AliESDtrackCuts();
	fCutsHybrid0 = new AliESDtrackCuts("fCutsHybrid0");
	fCutsHybrid0->SetMinNCrossedRowsTPC(70); //Default
	fCutsHybrid0->SetMinRatioCrossedRowsOverFindableClustersTPC(0.8); //
	fCutsHybrid0->SetMaxChi2PerClusterTPC(4); //
	fCutsHybrid0->SetAcceptKinkDaughters(kFALSE); //
	fCutsHybrid0->SetRequireTPCRefit(kTRUE); //
	fCutsHybrid0->SetRequireITSRefit(kTRUE); //
	fCutsHybrid0->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kAny);
	fCutsHybrid0->SetMaxDCAToVertexZ(2); //
	fCutsHybrid0->SetDCAToVertex2D(kFALSE); //
	fCutsHybrid0->SetRequireSigmaToVertex(kFALSE);
	fCutsHybrid0->SetMaxChi2PerClusterITS(36); //
	fCutsHybrid0->SetMaxDCAToVertexXY(2.4);//
	fCutsHybrid0->SetMaxDCAToVertexZ(3.2); //
	fTrackFilterHybrid0woDCA->AddCuts(fCutsHybrid0);
	fCutsHybrid0->SetMaxDCAToVertexXYPtDep("0.0105+0.0350/pt^1.1"); //
	fTrackFilterHybrid0->AddCuts(fCutsHybrid0);

	// Define Hybrid 1 (w/o SPD + ITSrefit)
	fTrackFilterHybrid1 = new AliAnalysisFilter("trackFilterHybrid1");
	fTrackFilterHybrid1woDCA = new AliAnalysisFilter("trackFilterHybrid1woDCA");
	AliESDtrackCuts * fCutsHybrid1 = new AliESDtrackCuts();
	fCutsHybrid1 = new AliESDtrackCuts("fCutsHybrid1");
	fCutsHybrid1->SetMinNCrossedRowsTPC(70); //Default
	fCutsHybrid1->SetMinRatioCrossedRowsOverFindableClustersTPC(0.8); //
	fCutsHybrid1->SetMaxChi2PerClusterTPC(4); //
	fCutsHybrid1->SetAcceptKinkDaughters(kFALSE); //
	fCutsHybrid1->SetRequireTPCRefit(kTRUE); //
	fCutsHybrid1->SetRequireITSRefit(kTRUE); //
	fCutsHybrid1->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kOff);
	fCutsHybrid1->SetMaxDCAToVertexZ(2); //
	fCutsHybrid1->SetDCAToVertex2D(kFALSE); //
	fCutsHybrid1->SetRequireSigmaToVertex(kFALSE);
	fCutsHybrid1->SetMaxChi2PerClusterITS(36); //
	fCutsHybrid1->SetMaxDCAToVertexXY(2.4);//
	fCutsHybrid1->SetMaxDCAToVertexZ(3.2); //
	fTrackFilterHybrid1woDCA->AddCuts(fCutsHybrid1);
	fCutsHybrid1->SetMaxDCAToVertexXYPtDep("0.0105+0.0350/pt^1.1"); //
	fTrackFilterHybrid1->AddCuts(fCutsHybrid1);

	// Define Hybrid 2 (w/o SPD + NoneSPD)
	fTrackFilterHybrid2 = new AliAnalysisFilter("trackFilterHybrid2");
	fTrackFilterHybrid2woDCA = new AliAnalysisFilter("trackFilterHybrid2woDCA");
	AliESDtrackCuts * fCutsHybrid2 = new AliESDtrackCuts();
	fCutsHybrid2 = new AliESDtrackCuts("fCutsHybrid2");
	fCutsHybrid2->SetMinNCrossedRowsTPC(70); //Default
	fCutsHybrid2->SetMinRatioCrossedRowsOverFindableClustersTPC(0.8); //
	fCutsHybrid2->SetMaxChi2PerClusterTPC(4); //
	fCutsHybrid2->SetAcceptKinkDaughters(kFALSE); //
	fCutsHybrid2->SetRequireTPCRefit(kTRUE); //
	fCutsHybrid2->SetRequireITSRefit(kFALSE); //
	fCutsHybrid2->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kNone);
	fCutsHybrid2->SetMaxDCAToVertexZ(2); //
	fCutsHybrid2->SetDCAToVertex2D(kFALSE); //
	fCutsHybrid2->SetRequireSigmaToVertex(kFALSE);
	fCutsHybrid2->SetMaxChi2PerClusterITS(36); //
	fCutsHybrid2->SetMaxDCAToVertexXY(2.4);//
	fCutsHybrid2->SetMaxDCAToVertexZ(3.2); //
	fTrackFilterHybrid2woDCA->AddCuts(fCutsHybrid2);
	fCutsHybrid2->SetMaxDCAToVertexXYPtDep("0.0105+0.0350/pt^1.1"); //
	fTrackFilterHybrid2->AddCuts(fCutsHybrid2);


	// create output objects

	OpenFile(1);
	fOutputList = new TList();          // this is a list which will contain all of your histograms
	// at the end of the analysis, the contents of this list are written  to the output file
	fOutputList->SetOwner(kTRUE);       // memory stuff: the list is owner of all objects and will delete them if requested

	if(fUseMC)
	{
		hNchTSGen = new TH1D("hNchTSGen","",3000,-0.5,2999.5);
		fOutputList->Add(hNchTSGen);

		hNchTSGenTest = new TH1D("hNchTSGenTest","",3000,-0.5,2999.5); 
		fOutputList->Add(hNchTSGenTest);

		hNchTSRecTest = new TH1D("hNchTSRecTest","",3000,-0.5,2999.5); 
		fOutputList->Add(hNchTSRecTest);

		hNchGen = new TH1D("hNchGen","",3000,-0.5,2999.5);
		fOutputList->Add(hNchGen);

		hNchGenTest = new TH1D("hNchGenTest","",3000,-0.5,2999.5); 
		fOutputList->Add(hNchGenTest);

		hNchRecTest = new TH1D("hNchRecTest","",3000,-0.5,2999.5); 
		fOutputList->Add(hNchRecTest);

		for(Int_t i=0;i<3;++i){
			hPhiGen[i]= new TH1D(Form("hPhiGen_%s",NameReg_3[i]),"",64,-TMath::Pi()/2.0,3.0*TMath::Pi()/2.0);
			fOutputList->Add(hPhiGen[i]);
		}

		hNchResponse = new TH2D("hNchResponse","Detector response; rec mult; gen mult",3000,-0.5,2999.5,3000,-0.5,2999.5);
		fOutputList->Add(hNchResponse);

		hPtInPrim = new TH1D("hPtInPrim","pT prim true; pT; Nch",ptNbins,ptbins1_3);
		fOutputList->Add(hPtInPrim);

		hPtInPrim_pion = new TH1D("hPtInPrim_pion","pT prim true; pT; Nch",ptNbins,ptbins1_3);
		fOutputList->Add(hPtInPrim_pion);

		hPtInPrim_kaon = new TH1D("hPtInPrim_kaon","pT prim true; pT; Nch",ptNbins,ptbins1_3);
		fOutputList->Add(hPtInPrim_kaon);

		hPtInPrim_proton = new TH1D("hPtInPrim_proton","pT prim true; pT; Nch",ptNbins,ptbins1_3);
		fOutputList->Add(hPtInPrim_proton);

		hPtInPrim_sigmap = new TH1D("hPtInPrim_sigmap","pT prim true; pT; Nch",ptNbins,ptbins1_3);
		fOutputList->Add(hPtInPrim_sigmap);

		hPtInPrim_sigmam = new TH1D("hPtInPrim_sigmam","pT prim true; pT; Nch",ptNbins,ptbins1_3);
		fOutputList->Add(hPtInPrim_sigmam);

		hPtInPrim_omega = new TH1D("hPtInPrim_omega","pT prim true; pT; Nch",ptNbins,ptbins1_3);
		fOutputList->Add(hPtInPrim_omega);

		hPtInPrim_xi = new TH1D("hPtInPrim_xi","pT prim true; pT; Nch",ptNbins,ptbins1_3);
		fOutputList->Add(hPtInPrim_xi);

		hPtInPrim_rest = new TH1D("hPtInPrim_rest","pT prim true; pT; Nch",ptNbins,ptbins1_3);
		fOutputList->Add(hPtInPrim_rest);

		hPtOut = new TH1D("hPtOut","pT all rec; pT; Nch",ptNbins,ptbins1_3);
		fOutputList->Add(hPtOut);

		hPtOutPrim = new TH1D("hPtOutPrim","pT prim rec; pT; Nch",ptNbins,ptbins1_3);
		fOutputList->Add(hPtOutPrim);

		hPtOutPrim_pion = new TH1D("hPtOutPrim_pion","pT prim rec; pT; Nch",ptNbins,ptbins1_3);
		fOutputList->Add(hPtOutPrim_pion);

		hPtOutPrim_kaon = new TH1D("hPtOutPrim_kaon","pT prim rec; pT; Nch",ptNbins,ptbins1_3);
		fOutputList->Add(hPtOutPrim_kaon);

		hPtOutPrim_proton = new TH1D("hPtOutPrim_proton","pT prim rec; pT; Nch",ptNbins,ptbins1_3);
		fOutputList->Add(hPtOutPrim_proton);

		hPtOutPrim_sigmap = new TH1D("hPtOutPrim_sigmap","pT prim rec; pT; Nch",ptNbins,ptbins1_3);
		fOutputList->Add(hPtOutPrim_sigmap);

		hPtOutPrim_sigmam = new TH1D("hPtOutPrim_sigmam","pT prim rec; pT; Nch",ptNbins,ptbins1_3);
		fOutputList->Add(hPtOutPrim_sigmam);

		hPtOutPrim_omega = new TH1D("hPtOutPrim_omega","pT prim rec; pT; Nch",ptNbins,ptbins1_3);
		fOutputList->Add(hPtOutPrim_omega);

		hPtOutPrim_xi = new TH1D("hPtOutPrim_xi","pT prim rec; pT; Nch",ptNbins,ptbins1_3);
		fOutputList->Add(hPtOutPrim_xi);

		hPtOutPrim_rest = new TH1D("hPtOutPrim_rest","pT prim rec; pT; Nch",ptNbins,ptbins1_3);
		fOutputList->Add(hPtOutPrim_rest);

		hPtOutSec = new TH1D("hPtOutSec","pT sec rec; pT; Nch",ptNbins,ptbins1_3);
		fOutputList->Add(hPtOutSec);
	}

	hNchTSRec = new TH1D("hNchTSRec","",3000,-0.5,2999.5);
	fOutputList->Add(hNchTSRec);

	hNchRec = new TH1D("hNchRec","",3000,-0.5,2999.5);
	fOutputList->Add(hNchTSRec);

	hNchTSData = new TH1D("hNchTSData","",3000,-0.5,2999.5); 
	fOutputList->Add(hNchTSData);

	hPhiTotal = new TH1F("hPhiSum","#varphi",50, -TMath::Pi()/2.0,5.0*TMath::Pi()/2.0);
	fOutputList->Add(hPhiTotal);

	hPhiStandard = new TH1F("hPhiSPD&ITS","#varphi",50, -TMath::Pi()/2.0,5.0*TMath::Pi()/2.0);
	fOutputList->Add(hPhiStandard);

	hPhiHybrid1 = new TH1F("hPhiITS","#varphi",50, -TMath::Pi()/2.0,5.0*TMath::Pi()/2.0);
	fOutputList->Add(hPhiHybrid1);

	hPhiHybrid2 = new TH1F("hPhiNITS","#varphi",50, -TMath::Pi()/2.0,5.0*TMath::Pi()/2.0);
	fOutputList->Add(hPhiHybrid2);

	hNchData = new TH1D("hNchData","",3000,-0.5,2999.5); 
	fOutputList->Add(hNchData);

	for(Int_t i=0;i<3;++i){
		hPhiRec[i]= new TH1D(Form("hPhiRec_%s",NameReg_3[i]),"",64,-TMath::Pi()/2.0,3.0*TMath::Pi()/2.0);
		fOutputList->Add(hPhiRec[i]);
	}

	for(Int_t i=0;i<3;++i){

		hPtVsUEGenTest[i] = new TH2D(Form("hPtVsUEGenTest_%s",NameReg_3[i]),"gen pT vs nch_transverse",3000,-0.5,2999.5,ptNbins,ptbins1_3);
		fOutputList->Add(hPtVsUEGenTest[i]);

		hPtVsUERecTest[i] = new TH2D(Form("hPtVsUERecTest_%s",NameReg_3[i]),"rec pT vs nch_transverse",3000,-0.5,2999.5,ptNbins,ptbins1_3);
		fOutputList->Add(hPtVsUERecTest[i]);	  

		hPtVsNchGenTest[i] = new TH2D(Form("hPtVsNchGenTest_%s",NameReg_3[i]),"gen pT vs nch_transverse",3000,-0.5,2999.5,ptNbins,ptbins1_3);
		fOutputList->Add(hPtVsNchGenTest[i]);

		hPtVsNchRecTest[i] = new TH2D(Form("hPtVsNchRecTest_%s",NameReg_3[i]),"rec pT vs nch_transverse",3000,-0.5,2999.5,ptNbins,ptbins1_3);
		fOutputList->Add(hPtVsNchRecTest[i]);

	}


	for(Int_t i=0;i<3;++i){

		hPtVsUEData[i] = new TH2D(Form("hPtVsUEData_%s",NameReg_3[i]),"data pT vs nch_transverse",3000,-0.5,2999.5,ptNbins,ptbins1_3);
		fOutputList->Add(hPtVsUEData[i]);

		hPtVsNchData[i] = new TH2D(Form("hPtVsNchData_%s",NameReg_3[i]),"data pT vs nch_transverse",3000,-0.5,2999.5,ptNbins,ptbins1_3);
		fOutputList->Add(hPtVsNchData[i]);

	}

	hPTVsDCAData = new TH2D("hPTVsDCAData","pT vs dcaxy",ptNbins,ptbins1_3,nBinsDCAxy,binsDCAxy3);
	fOutputList->Add(hPTVsDCAData);

	hptvsdcaPrim = new TH2F("hptvsdcaPrim","pt vs dca Primaries",ptNbins,ptbins1_3,nBinsDCAxy,binsDCAxy3);
	fOutputList->Add(hptvsdcaPrim);

	hptvsdcaDecs = new TH2F("hptvsdcaDecs","pt vs dca Decays",ptNbins,ptbins1_3,nBinsDCAxy,binsDCAxy3);
	fOutputList->Add(hptvsdcaDecs);

	hptvsdcaMatl = new TH2F("hptvsdcaMatl","pt vs dca Material",ptNbins,ptbins1_3,nBinsDCAxy,binsDCAxy3);
	fOutputList->Add(hptvsdcaMatl);

	hptvsdcaAll = new TH2F("hptvsdcaAll","pt vs dca all",ptNbins,ptbins1_3,nBinsDCAxy,binsDCAxy3);
	fOutputList->Add(hptvsdcaAll);

	fEventCuts.AddQAplotsToList(fOutputList);
	PostData(1, fOutputList);           // postdata will notify the analysis manager of changes / updates to the

}
//_____________________________________________________________________________
void AliAnalysisTaskChargedVsRT::UserExec(Option_t *)
{

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
	
	if (fUseMC) {

		//      E S D
		fMC = dynamic_cast<AliMCEvent*>(MCEvent());
		if(!fMC){
			Printf("%s:%d MCEvent not found in Input Manager",(char*)__FILE__,__LINE__);
			this->Dump();
			return;
		}
		fMCStack = fMC->Stack();
	}

	AliHeader* headerMC = fMC->Header();
	Bool_t isGoodVtxPosMC = kFALSE;

	vector<Float_t> ptMc;
	vector<Float_t> phiMc;
	vector<Int_t> idMc;
	if (fUseMC){
		AliGenEventHeader* genHeader = headerMC->GenEventHeader();
		TArrayF vtxMC(3); // primary vertex  MC 
		vtxMC[0]=9999; vtxMC[1]=9999;  vtxMC[2]=9999; //initialize with dummy
		if (genHeader) {
			genHeader->PrimaryVertex(vtxMC);
		}
		if(TMath::Abs(vtxMC[2])<=10)
			isGoodVtxPosMC = kTRUE;

		fnGen = FillArrayMC( ptMc,  phiMc, idMc );

		// Before trigger selection
		GetLeadingObjectFromArray(ptMc,phiMc,fnGen, kTRUE);
		// Filling histos for observables at generator level
		if( fGenLeadPt>=fLeadPtCutMin && fGenLeadPt<fLeadPtCutMax )
			GetMultiplicityDistributionsTrue(phiMc,ptMc,fnGen);


	}



	// Define the arrays for hybrid, standard (w DCA cut)
	vector<Float_t> ptHy;
	vector<Float_t> phiHy;
	vector<Float_t> dcaxyHy;
	vector<Float_t> dcazHy;
	vector<Int_t>   isprimHy;
	vector<Int_t>   idHy;
	fnRecHy = FillArray( ptHy, phiHy, dcaxyHy, dcazHy, isprimHy, idHy, kTRUE, fIsHybAna );

	// Before trigger selection
	GetLeadingObjectFromArray(ptHy,phiHy,fnRecHy,kFALSE);
	// Now we get the leading particle pT

	// Define the arrays for hybrid, standard (wo DCA cut)
	vector<Float_t> ptHyWoDCA;
	vector<Float_t> phiHyWoDCA;
	vector<Float_t> dcaxyHyWoDCA;
	vector<Float_t> dcazHyWoDCA;
	vector<Int_t>   isprimHyWoDCA;
	vector<Int_t>   idHyWoDCA;
	fnRecHyWoDCA = FillArray( ptHyWoDCA, phiHyWoDCA, dcaxyHyWoDCA, dcazHyWoDCA, isprimHyWoDCA, idHyWoDCA, kFALSE, fIsHybAna );

	// Trigger selection
	UInt_t fSelectMask= fInputHandler->IsEventSelected();
	Bool_t isINT7selected = fSelectMask&AliVEvent::kINT7;
	if(!isINT7selected)
		return;

	// Good events
	if (!fEventCuts.AcceptEvent(event)) {
		PostData(1, fOutputList);
		return;
	}

	fMultSelectionbefvtx = (AliMultSelection*) fESD->FindListObject("MultSelection");

	// Good vertex
	Bool_t hasRecVertex = kFALSE;
	hasRecVertex=HasRecVertex();
	if(!hasRecVertex)return;

	// Multiplicity Estimation
	fv0mpercentile = -999;

	fMultSelection = (AliMultSelection*) fESD->FindListObject("MultSelection");
	if (!fMultSelection)
		cout<<"------- No AliMultSelection Object Found --------"<<fMultSelection<<endl;
	fv0mpercentile = fMultSelection->GetMultiplicityPercentile("V0M");
	
	if(fIsMCclosure){
		Double_t randomUE = -1;
		gRandom->SetSeed(0);
		randomUE = gRandom->Uniform(0.0,1.0);
		if(randomUE<0.5){                   // corrections (50% stat.)
			if(isGoodVtxPosMC){
				// KNO scaling
				if( ( fGenLeadPt>=fLeadPtCutMin && fGenLeadPt<fLeadPtCutMax ) && ( fRecLeadPt>=fLeadPtCutMin && fRecLeadPt<fLeadPtCutMax )){
					GetDetectorResponse(phiMc,fnGen,phiHy,fnRecHy);
					GetBinByBinCorrections(fnGen,fnRecHy,ptMc,ptHy,idMc,idHy,isprimHy);
				}
			}
		}
		else{// for testing the method
			// KNO scaling
			if(( fRecLeadPt>=fLeadPtCutMin && fRecLeadPt<fLeadPtCutMax ))
			{
				GetMultiplicityDistributions(phiHy,ptHy,fnRecHy,ptHyWoDCA,dcaxyHyWoDCA,isprimHyWoDCA, fnRecHyWoDCA);
			}
		}
	}
	else{
		if(fUseMC){
			if(isGoodVtxPosMC){
				// KNO scaling
				if( ( fGenLeadPt>=fLeadPtCutMin && fGenLeadPt<fLeadPtCutMax ) && ( fRecLeadPt>=fLeadPtCutMin && fRecLeadPt<fLeadPtCutMax ))
				{
					GetDetectorResponse(phiMc,fnGen,phiHy,fnRecHy);
					GetBinByBinCorrections(fnGen,fnRecHy,ptMc,ptHy,idMc,idHy,isprimHy);
					GetMultiplicityDistributions(phiHy,ptHy,fnRecHy,ptHyWoDCA,dcaxyHyWoDCA,isprimHyWoDCA, fnRecHyWoDCA);
				}
			}
		}
		else{
			// KNO scaling
			if(( fRecLeadPt>=fLeadPtCutMin && fRecLeadPt<fLeadPtCutMax ))
				GetMultiplicityDistributionsData(phiHy,ptHy,fnRecHy,ptHyWoDCA,dcaxyHyWoDCA,fnRecHyWoDCA);
		}
	}
	//}

	PostData(1, fOutputList); // stream the result of this event to the output manager which will write it to a file

}


//______________________________________________________________________________
void AliAnalysisTaskChargedVsRT::Terminate(Option_t *)
{

}

//_____________________________________________________________________________
void AliAnalysisTaskChargedVsRT::GetLeadingObjectFromArray(const vector<Float_t> &pt,const vector<Float_t> &phi, Int_t multPart, Bool_t isMC) {

	Float_t flPt = 0;// leading pT
	Float_t flPhi = 0;
	Int_t flIndex = 0;


	for(Int_t i1 = 0; i1 < multPart; ++i1){
		if ( flPt < pt[i1] ){
			flPt  = pt[i1];
			flPhi = phi[i1];
			flIndex = i1;
		}

	}
	if(isMC){
		fGenLeadPhi = flPhi;
		fGenLeadPt  = flPt;
		fGenLeadIn  = flIndex;
	}else{
		fRecLeadPhi = flPhi;
		fRecLeadPt  = flPt;
		fRecLeadIn  = flIndex;
	}
}

//_____________________________________________________________________________
void AliAnalysisTaskChargedVsRT::GetBinByBinCorrections( Int_t multGen, Int_t multRec, const vector<Float_t> &ptGen, const vector<Float_t> &ptRec, const vector<Int_t> &idGen, const vector<Int_t> &idRec, const vector<Int_t> &isprimRec){


	// Histos for efficiencyxacceptance
	for (Int_t i = 0; i < multGen; ++i) {
		hPtInPrim->Fill(ptGen[i]);// inital pT distribution (MC gen)
		if (idGen[i]==0) hPtInPrim_pion->Fill(ptGen[i]); //pions
		else if (idGen[i]==1) hPtInPrim_kaon->Fill(ptGen[i]); //kaons
		else if (idGen[i]==2) hPtInPrim_proton->Fill(ptGen[i]); //protons
		else if (idGen[i]==3) hPtInPrim_sigmap->Fill(ptGen[i]); //sigma plus
		else if (idGen[i]==4) hPtInPrim_sigmam->Fill(ptGen[i]); //sigma minus
		else if (idGen[i]==5) hPtInPrim_omega->Fill(ptGen[i]); //Omega
		else if (idGen[i]==6) hPtInPrim_xi->Fill(ptGen[i]); //Xi
		else hPtInPrim_rest->Fill(ptGen[i]); //rest of the charged particles
	}


	for(Int_t i=0; i < multRec; ++i) {  // loop over all these tracks
		hPtOut->Fill(ptRec[i]);
		if( isprimRec[i] == 0 ){
			hPtOutPrim->Fill(ptRec[i]);
			if (idRec[i]==0) hPtOutPrim_pion->Fill(ptRec[i]); //pions
			else if (idRec[i]==1) hPtOutPrim_kaon->Fill(ptRec[i]); //kaons
			else if (idRec[i]==2) hPtOutPrim_proton->Fill(ptRec[i]); //protons
			else if (idRec[i]==3) hPtOutPrim_sigmap->Fill(ptRec[i]); //sigma plus
			else if (idRec[i]==4) hPtOutPrim_sigmam->Fill(ptRec[i]); //sigma minus
			else if (idRec[i]==5) hPtOutPrim_omega->Fill(ptRec[i]); //Omega
			else if (idRec[i]==6) hPtOutPrim_xi->Fill(ptRec[i]); //Xi
			else hPtOutPrim_rest->Fill(ptRec[i]); //rest of the charged particles
		}
		if( isprimRec[i] == 1 || isprimRec[i] == 2 ){
			hPtOutSec->Fill(ptRec[i]);
		}

	}
}

void AliAnalysisTaskChargedVsRT::GetDetectorResponse(const vector<Float_t> &phiGen, Int_t multGen, const vector<Float_t> &phiRec, Int_t multRec) {

	Int_t multTSgen=0;
	Int_t multTSrec=0;
	for (Int_t i = 0; i < multGen; ++i) {

		if(i==fGenLeadIn)continue;
		Double_t DPhi = DeltaPhi(phiGen[i], fGenLeadPhi);

		// definition of the topological regions
		if(TMath::Abs(DPhi)<pi/3.0){// near side
			hPhiGen[0]->Fill(DPhi);
		}
		else if(TMath::Abs(DPhi-pi)<pi/3.0){// away side
			hPhiGen[1]->Fill(DPhi);
		}
		else{// transverse side
			multTSgen++;
			hPhiGen[2]->Fill(DPhi);
		}
	}

	hNchTSGen->Fill(multTSgen);

	for(Int_t i=0; i < multRec; ++i) {  // loop over all these tracks

		if(i==fRecLeadIn)continue;
		Double_t DPhi = DeltaPhi(phiRec[i], fRecLeadPhi);

		// definition of the topological regions
		if(TMath::Abs(DPhi)<pi/3.0){// near side
			hPhiRec[0]->Fill(DPhi);
		}
		else if(TMath::Abs(DPhi-pi)<pi/3.0){// away side
			hPhiRec[1]->Fill(DPhi);
		}
		else{// transverse side
			multTSrec++;
			hPhiRec[2]->Fill(DPhi);
		}

	}

	hNchTSRec->Fill(multTSrec);
	hNchResponse->Fill(multTSrec,multTSgen);

}
//______________________________________________________________
void AliAnalysisTaskChargedVsRT::GetMultiplicityDistributionsData(const vector<Float_t> &phiRec, const vector<Float_t> &ptRec, Int_t multRec,  const vector<Float_t> &ptRecWoDCA, const vector<Float_t> &dcaxyRecWoDCA, Int_t multRecWoDCA){

	Int_t multTSrec=0;

	for(Int_t i=0; i < multRec; ++i) {                 // loop over all these tracks

		if(i==fRecLeadIn) continue;
		Double_t DPhi = DeltaPhi(phiRec[i], fRecLeadPhi);

		// definition of the topological regions
		if(TMath::Abs(DPhi)<pi/3.0){// near side
			continue;
		}
		else if(TMath::Abs(DPhi-pi)<pi/3.0){// away side
			continue;
		}
		else{// transverse side
			multTSrec++;
		}

	}

	hNchTSData->Fill(multTSrec);
	hNchData->Fill(multRec);

	// Filling rec pT vs UE (for pT *** considering the hybrid track cuts or the 2015 track cuts ***)
	for(Int_t i=0; i < multRec; ++i) {                 // loop over all these tracks

		if(i==fRecLeadIn)continue;
		Double_t DPhi = DeltaPhi(phiRec[i], fRecLeadPhi);

		// definition of the topological regions
		if(TMath::Abs(DPhi)<pi/3.0){// near side
			hPtVsUEData[0]->Fill(multTSrec,ptRec[i]);
			hPtVsNchData[0]->Fill(multRec,ptRec[i]);
		}
		else if(TMath::Abs(DPhi-pi)<pi/3.0){// away side
			hPtVsUEData[1]->Fill(multTSrec,ptRec[i]);
			hPtVsNchData[1]->Fill(multRec,ptRec[i]);
		}
		else{// transverse side
			hPtVsUEData[2]->Fill(multTSrec,ptRec[i]);
			hPtVsNchData[2]->Fill(multRec,ptRec[i]);
		}

	}

	// Filling the histos wo DCA cuts
	for(Int_t i=0; i < multRecWoDCA; ++i) {

		hPTVsDCAData->Fill(ptRecWoDCA[i],dcaxyRecWoDCA[i]);

	}


}
//_____________________________________________________________
void AliAnalysisTaskChargedVsRT::GetMultiplicityDistributionsTrue(const vector<Float_t> &phiGen, const vector<Float_t> &ptGen, Int_t multGen){


	Int_t multTSgen=0;
	//Int_t multTSrec=0;

	for (Int_t i = 0; i < multGen; ++i) {

		if(i==fGenLeadIn)continue;

		Double_t DPhi = DeltaPhi(phiGen[i], fGenLeadPhi);

		// definition of the topological regions
		if(TMath::Abs(DPhi)<pi/3.0){// near side
			continue;
		}
		else if(TMath::Abs(DPhi-pi)<pi/3.0){// away side
			continue;
		}
		else{// transverse side
			multTSgen++;
		}

	}

	hNchTSGenTest->Fill(multTSgen);
	hNchGenTest->Fill(multGen);

	for (Int_t i = 0; i < multGen; ++i) {

		if(i==fGenLeadIn)continue;
		Double_t DPhi = DeltaPhi(phiGen[i], fGenLeadPhi);

		// definition of the topological regions
		if(TMath::Abs(DPhi)<pi/3.0){// near side
			hPtVsUEGenTest[0]->Fill(multTSgen,ptGen[i]);
			hPtVsNchGenTest[0]->Fill(multGen,ptGen[i]);
		}
		else if(TMath::Abs(DPhi-pi)<pi/3.0){// away side
			hPtVsUEGenTest[1]->Fill(multTSgen,ptGen[i]);
			hPtVsNchGenTest[1]->Fill(multGen,ptGen[i]);
		}
		else{// transverse side
			hPtVsUEGenTest[2]->Fill(multTSgen,ptGen[i]);
			hPtVsNchGenTest[2]->Fill(multGen,ptGen[i]);
		}

	}

}

//____________________________________________________________
void AliAnalysisTaskChargedVsRT::GetMultiplicityDistributions(const vector<Float_t> &phiRec, const vector<Float_t> &ptRec, Int_t multRec,  const vector<Float_t> &ptRecWoDCA, const vector<Float_t> &dcaxyRecWoDCA, const vector<Int_t> &isprimRecWoDCA, Int_t multRecWoDCA){
	Int_t multTSrec=0;
	// see how many tracks there are in the event
	for(Int_t i=0; i < multRec; ++i) {                 // loop over all these tracks

		if(i==fRecLeadIn)continue;
		Double_t DPhi = DeltaPhi(phiRec[i], fRecLeadPhi);

		// definition of the topological regions
		if(TMath::Abs(DPhi)<pi/3.0){// near side
			continue;
		}
		else if(TMath::Abs(DPhi-pi)<pi/3.0){// away side
			continue;
		}
		else{// transverse side
			multTSrec++;
		}

	}

	hNchTSRecTest->Fill(multTSrec);
	hNchRecTest->Fill(multRec);

	// Filling rec pT vs UE (for pT *** considering the hybrid track cuts***)

	for(Int_t i=0; i < multRec; ++i) {                 // loop over all these tracks

		if(i==fRecLeadIn)continue;
		Double_t DPhi = DeltaPhi(phiRec[i], fRecLeadPhi);

		// definition of the topological regions
		if(TMath::Abs(DPhi)<pi/3.0){// near side
			hPtVsUERecTest[0]->Fill(multTSrec,ptRec[i]);
			hPtVsNchRecTest[0]->Fill(multRec,ptRec[i]);
		}
		else if(TMath::Abs(DPhi-pi)<pi/3.0){// away side
			hPtVsUERecTest[1]->Fill(multTSrec,ptRec[i]);
			hPtVsNchRecTest[1]->Fill(multRec,ptRec[i]);
		}
		else{// transverse side
			hPtVsUERecTest[2]->Fill(multTSrec,ptRec[i]);
			hPtVsNchRecTest[2]->Fill(multRec,ptRec[i]);
		}

	}

	// Auxiliar distributions to calculate the contamination from secondary particles, it runs over tracks wo DCA cut
	for(Int_t i=0; i < multRecWoDCA; ++i) {                 // loop over all these tracks

		hptvsdcaAll->Fill(ptRecWoDCA[i],dcaxyRecWoDCA[i]);


		if( fUseMC && ( fGenLeadPt>=fLeadPtCutMin && fGenLeadPt<fLeadPtCutMax )){
			if(isprimRecWoDCA[i]==0){
				hptvsdcaPrim->Fill(ptRecWoDCA[i],dcaxyRecWoDCA[i]);
			}
			else if(isprimRecWoDCA[i]==1){
				hptvsdcaDecs->Fill(ptRecWoDCA[i],dcaxyRecWoDCA[i]);
			}
			else if(isprimRecWoDCA[i]==2){
				hptvsdcaMatl->Fill(ptRecWoDCA[i],dcaxyRecWoDCA[i]);
			}

		}
	}



}

Double_t AliAnalysisTaskChargedVsRT::DeltaPhi(Double_t phia, Double_t phib,
		Double_t rangeMin, Double_t rangeMax)
{
	Double_t dphi = -999;
	Double_t pi = TMath::Pi();

	if (phia < 0)         phia += 2*pi;
	else if (phia > 2*pi) phia -= 2*pi;
	if (phib < 0)         phib += 2*pi;
	else if (phib > 2*pi) phib -= 2*pi;
	dphi = phib - phia;
	if (dphi < rangeMin)      dphi += 2*pi;
	else if (dphi > rangeMax) dphi -= 2*pi;

	return dphi;
}
Bool_t AliAnalysisTaskChargedVsRT::HasRecVertex(){


	float fMaxDeltaSpdTrackAbsolute = 0.5f;
	float fMaxDeltaSpdTrackNsigmaSPD = 1.e14f;
	float fMaxDeltaSpdTrackNsigmaTrack = 1.e14;
	float fMaxResolutionSPDvertex = 0.25f;
	float fMaxDispersionSPDvertex = 1.e14f;

	Bool_t fRequireTrackVertex = true;
	unsigned long fFlag;
	fFlag =BIT(AliEventCuts::kNoCuts);

	const AliVVertex* vtTrc = fESD->GetPrimaryVertex();
	bool isTrackV = true;
	if(vtTrc->IsFromVertexer3D() || vtTrc->IsFromVertexerZ()) isTrackV=false;
	const AliVVertex* vtSPD = fESD->GetPrimaryVertexSPD();


	if (vtSPD->GetNContributors() > 0) fFlag |= BIT(AliEventCuts::kVertexSPD);

	if (vtTrc->GetNContributors() > 1 && isTrackV ) fFlag |= BIT(AliEventCuts::kVertexTracks);

	if (((fFlag & BIT(AliEventCuts::kVertexTracks)) ||  !fRequireTrackVertex) && (fFlag & BIT(AliEventCuts::kVertexSPD))) fFlag |= BIT(AliEventCuts::kVertex);

	const AliVVertex* &vtx = bool(fFlag & BIT(AliEventCuts::kVertexTracks)) ? vtTrc : vtSPD;
	AliVVertex   *fPrimaryVertex = const_cast<AliVVertex*>(vtx);
	if(!fPrimaryVertex)return kFALSE;

	/// Vertex quality cuts
	double covTrc[6],covSPD[6];
	vtTrc->GetCovarianceMatrix(covTrc);
	vtSPD->GetCovarianceMatrix(covSPD);
	double dz = bool(fFlag & AliEventCuts::kVertexSPD) && bool(fFlag & AliEventCuts::kVertexTracks) ? vtTrc->GetZ() - vtSPD->GetZ() : 0.; /// If one of the two vertices is not available this cut is always passed.
	double errTot = TMath::Sqrt(covTrc[5]+covSPD[5]);
	double errTrc = bool(fFlag & AliEventCuts::kVertexTracks) ? TMath::Sqrt(covTrc[5]) : 1.;
	double nsigTot = TMath::Abs(dz) / errTot, nsigTrc = TMath::Abs(dz) / errTrc;
	/// vertex dispersion for run1, only for ESD, AOD code to be added here
	const AliESDVertex* vtSPDESD = dynamic_cast<const AliESDVertex*>(vtSPD);
	double vtSPDdispersion = vtSPDESD ? vtSPDESD->GetDispersion() : 0;
	if (
			(TMath::Abs(dz) <= fMaxDeltaSpdTrackAbsolute && nsigTot <= fMaxDeltaSpdTrackNsigmaSPD && nsigTrc <= fMaxDeltaSpdTrackNsigmaTrack) && // discrepancy track-SPD vertex
			(!vtSPD->IsFromVertexerZ() || TMath::Sqrt(covSPD[5]) <= fMaxResolutionSPDvertex) &&
			(!vtSPD->IsFromVertexerZ() || vtSPDdispersion <= fMaxDispersionSPDvertex) /// vertex dispersion cut for run1, only for ESD
	   ) // quality cut on vertexer SPD z
		fFlag |= BIT(AliEventCuts::kVertexQuality);  

	Bool_t hasVtx = (TESTBIT(fFlag,AliEventCuts::kVertex))&&(TESTBIT(fFlag,AliEventCuts::kVertexQuality));

	return hasVtx;

}
//______________________________________________________________________
Int_t AliAnalysisTaskChargedVsRT::FillArrayMC( vector<Float_t> &ptArray, vector<Float_t> &phiArray, vector<Int_t> &idArray ){
	/*
	   id 0: pion, 1: kaon, 2: proton, 3: sigma plus, 4: sigma minus, 5: Omega, 6: Xi, 7: other charged
	 */

	ptArray.clear();
	phiArray.clear();
	idArray.clear();
	Int_t nNchGen    = 0;

	for (Int_t i = 0; i < fMC->GetNumberOfTracks(); ++i) {

		AliMCParticle* particle = (AliMCParticle*)fMC->GetTrack(i);
		if (!particle) continue;
		if (!fMC->IsPhysicalPrimary(i)) continue;
		if (particle->Charge() == 0) continue;
		if ( TMath::Abs(particle->Eta()) > fEtaCut )continue;
		if( particle->Pt() < fPtMin)continue;

		Int_t idPart = -1;
		Int_t partPDG = TMath::Abs(particle->PdgCode());
		if (partPDG==211) idPart = 0; //pions
		else if (partPDG==321) idPart = 1; //kaons
		else if (partPDG==2212) idPart = 2; //protons
		else if (partPDG==3222) idPart = 3; //sigma plus
		else if (partPDG==3112) idPart = 4; //sigma minus
		else if (partPDG==3334) idPart = 5; //Omega
		else if (partPDG==3312) idPart = 6; //Xi
		else idPart = 7; //rest of the charged particles

		ptArray.push_back(particle->Pt());
		phiArray.push_back(particle->Phi());
		idArray.push_back(idPart);
		nNchGen++;

	}
	return nNchGen;


}
//_____________________________
Int_t AliAnalysisTaskChargedVsRT::FillArray( vector<Float_t> &ptArray, vector<Float_t> &phiArray, vector<Float_t> &dcaxyArray, vector<Float_t> &dcazArray, vector<Int_t> &isprimArray, vector<Int_t> &idArray, const bool wDcaCut, const bool useHy ){
	/*
	   id 0: pion, 1: kaon, 2: proton, 3: sigma plus, 4: sigma minus, 5: Omega, 6: Xi, 7: other charged
	 */
	ptArray.clear();
	phiArray.clear();
	dcaxyArray.clear();
	dcazArray.clear();
	isprimArray.clear();
	idArray.clear();

	Int_t nTracks = fESD->GetNumberOfTracks();
	Int_t nNchRec    = 0;
	if(wDcaCut){// with DCA cut
		for(Int_t iT = 0; iT < nTracks; ++iT) {

			AliESDtrack* esdtrack = static_cast<AliESDtrack*>(fESD->GetTrack(iT));  // get a track (type AliesdTrack)
			if(!esdtrack) continue;
			fdcaxy = -999;
			fdcaz = -999;
			if(TMath::Abs(esdtrack->Eta()) > fEtaCut) continue;
			if( esdtrack->Pt() < fPtMin)continue;
			AliESDtrack *newTrack = 0x0;
			Int_t isPrim = -1;
			Int_t idTrack = -1;
			Int_t mcLabel = -1;
			if(fUseMC){// get label: 0: prim, 1: weak decays, 2: material
				mcLabel = TMath::Abs(esdtrack->GetLabel());
				TParticle *mcParticle = fMC->GetTrack(mcLabel)->Particle();
				if(!mcParticle){
					printf("----ERROR: mcParticle not available------------------\n");
					continue;
				}
				Int_t partPDG_rec = TMath::Abs(mcParticle->GetPdgCode());
				if (partPDG_rec==211) idTrack = 0; //pions
				else if (partPDG_rec==321) idTrack = 1; //kaons
				else if (partPDG_rec==2212) idTrack = 2; //protons
				else if (partPDG_rec==3222) idTrack = 3; //sigma plus
				else if (partPDG_rec==3112) idTrack = 4; //sigma minus
				else if (partPDG_rec==3334) idTrack = 5; //Omega
				else if (partPDG_rec==3312) idTrack = 6; //Xi
				else idTrack = 7; //rest of the charged particles
			}
			Bool_t isHy0=kFALSE;
			Bool_t isHy1=kFALSE;
			Bool_t isHy2=kFALSE;
			if(useHy){
				if(fTrackFilterHybrid0->IsSelected(esdtrack))
					isHy0 = kTRUE;
				else if(fTrackFilterHybrid1->IsSelected(esdtrack))
					isHy1 = kTRUE;
				else if(fTrackFilterHybrid2->IsSelected(esdtrack))
					isHy2 = kTRUE;
			}
			else{
				if(fTrackFilter2015->IsSelected(esdtrack))
					isHy0 = kTRUE;
			}

			if(isHy0){
				newTrack = new AliESDtrack(*esdtrack);
				hPhiTotal->Fill(esdtrack->Phi());
				hPhiStandard->Fill(newTrack->Phi());
				newTrack->GetImpactParameters(fdcaxy,fdcaz);
				ptArray.push_back(newTrack->Pt());
				phiArray.push_back(newTrack->Phi());
				dcaxyArray.push_back(fdcaxy);
				dcazArray.push_back(fdcaz);
				if(fUseMC){// get label: 0: prim, 1: weak decays, 2: material
					if(fMC->IsPhysicalPrimary(mcLabel))
						isPrim = 0;
					else if (fMC->IsSecondaryFromWeakDecay(mcLabel))
						isPrim = 1;
					else if(fMC->IsSecondaryFromMaterial(mcLabel))
						isPrim = 2;
					else
						continue;
				}
				isprimArray.push_back(isPrim);
				idArray.push_back(idTrack);
				nNchRec++;
			}
			else if(isHy1){
				newTrack = new AliESDtrack(*esdtrack);
				if(esdtrack->GetConstrainedParam()){
					const AliExternalTrackParam* constrainParam = esdtrack->GetConstrainedParam();
					newTrack->Set(constrainParam->GetX(),constrainParam->GetAlpha(),constrainParam->GetParameter(),constrainParam->GetCovariance());
					hPhiTotal->Fill(esdtrack->Phi());
					hPhiHybrid1->Fill(newTrack->Phi());
					newTrack->GetImpactParameters(fdcaxy,fdcaz);
					ptArray.push_back(newTrack->Pt());
					phiArray.push_back(newTrack->Phi());
					dcaxyArray.push_back(fdcaxy);
					dcazArray.push_back(fdcaz);
					if(fUseMC){// get label: 0: prim, 1: weak decays, 2: material
						if(fMC->IsPhysicalPrimary(mcLabel))
							isPrim = 0;
						else if (fMC->IsSecondaryFromWeakDecay(mcLabel))
							isPrim = 1;
						else if(fMC->IsSecondaryFromMaterial(mcLabel))
							isPrim = 2;
						else
							continue;
					}
					isprimArray.push_back(isPrim);
					idArray.push_back(idTrack);
					nNchRec++;
				}
			}
			else if(isHy2){
				newTrack = new AliESDtrack(*esdtrack);
				if(esdtrack->GetConstrainedParam()){
					const AliExternalTrackParam* constrainParam = esdtrack->GetConstrainedParam();
					newTrack->Set(constrainParam->GetX(),constrainParam->GetAlpha(),constrainParam->GetParameter(),constrainParam->GetCovariance());
					hPhiTotal->Fill(esdtrack->Phi());
					hPhiHybrid2->Fill(newTrack->Phi());
					newTrack->GetImpactParameters(fdcaxy,fdcaz);
					ptArray.push_back(newTrack->Pt());
					phiArray.push_back(newTrack->Phi());
					dcaxyArray.push_back(fdcaxy);
					dcazArray.push_back(fdcaz);
					if(fUseMC){// get label: 0: prim, 1: weak decays, 2: material
						if(fMC->IsPhysicalPrimary(mcLabel))
							isPrim = 0;
						else if (fMC->IsSecondaryFromWeakDecay(mcLabel))
							isPrim = 1;
						else if(fMC->IsSecondaryFromMaterial(mcLabel))
							isPrim = 2;
						else
							continue;
					}
					isprimArray.push_back(isPrim);
					idArray.push_back(idTrack);
					nNchRec++;

				}
			}

		}
	}else{
		for(Int_t iT = 0; iT < nTracks; ++iT) {

			AliESDtrack* esdtrack = static_cast<AliESDtrack*>(fESD->GetTrack(iT));  // get a track (type AliesdTrack)
			if(!esdtrack) continue;
			fdcaxy = -999;
			fdcaz = -999;
			if(TMath::Abs(esdtrack->Eta()) > fEtaCut) continue;
			if( esdtrack->Pt() < fPtMin)continue;
			AliESDtrack *newTrack = 0x0;
			Int_t isPrim = -1;
			Int_t idTrack = -1;
			Int_t mcLabel = -1;
			TParticle *mcParticle = 0;
			if(fUseMC){// get label: 0: prim, 1: weak decays, 2: material
				mcLabel = TMath::Abs(esdtrack->GetLabel());
				mcParticle = fMC->GetTrack(mcLabel)->Particle();
				if(!mcParticle){
					printf("----ERROR: mcParticle not available------------------\n");
					continue;
				}

				Int_t partPDG_rec = TMath::Abs(mcParticle->GetPdgCode());
				if (partPDG_rec==211) idTrack = 0; //pions
				else if (partPDG_rec==321) idTrack = 1; //kaons
				else if (partPDG_rec==2212) idTrack = 2; //protons
				else if (partPDG_rec==3222) idTrack = 3; //sigma plus
				else if (partPDG_rec==3112) idTrack = 4; //sigma minus
				else if (partPDG_rec==3334) idTrack = 5; //Omega
				else if (partPDG_rec==3312) idTrack = 6; //Xi
				else idTrack = 7; //rest of the charged particles
			}
			Bool_t isHy0=kFALSE;
			Bool_t isHy1=kFALSE;
			Bool_t isHy2=kFALSE;
			if(useHy){
				if(fTrackFilterHybrid0woDCA->IsSelected(esdtrack))
					isHy0 = kTRUE;
				else if(fTrackFilterHybrid1woDCA->IsSelected(esdtrack))
					isHy1 = kTRUE;
				else if(fTrackFilterHybrid2woDCA->IsSelected(esdtrack))
					isHy2 = kTRUE;
			}
			else{
				if(fTrackFilter2015woDCA->IsSelected(esdtrack))
					isHy0 = kTRUE;
			}
			if(isHy0){
				newTrack = new AliESDtrack(*esdtrack);
				newTrack->GetImpactParameters(fdcaxy,fdcaz);
				ptArray.push_back(newTrack->Pt());
				phiArray.push_back(newTrack->Phi());
				dcaxyArray.push_back(fdcaxy);
				dcazArray.push_back(fdcaz);
				if(fUseMC){// get label: 0: prim, 1: weak decays, 2: material
					if(fMC->IsPhysicalPrimary(mcLabel))
						isPrim = 0;
					else if (fMC->IsSecondaryFromWeakDecay(mcLabel))
						isPrim = 1;
					else if(fMC->IsSecondaryFromMaterial(mcLabel))
						isPrim = 2;
					else
						continue;
				}
				isprimArray.push_back(isPrim);
				idArray.push_back(idTrack);
				nNchRec++;
			}
			else if(isHy1){
				newTrack = new AliESDtrack(*esdtrack);
				if(esdtrack->GetConstrainedParam()){
					const AliExternalTrackParam* constrainParam = esdtrack->GetConstrainedParam();
					newTrack->Set(constrainParam->GetX(),constrainParam->GetAlpha(),constrainParam->GetParameter(),constrainParam->GetCovariance());
					newTrack->GetImpactParameters(fdcaxy,fdcaz);
					ptArray.push_back(newTrack->Pt());
					phiArray.push_back(newTrack->Phi());
					dcaxyArray.push_back(fdcaxy);
					dcazArray.push_back(fdcaz);
					if(fUseMC){// get label: 0: prim, 1: weak decays, 2: material
						if(fMC->IsPhysicalPrimary(mcLabel))
							isPrim = 0;
						else if (fMC->IsSecondaryFromWeakDecay(mcLabel))
							isPrim = 1;
						else if(fMC->IsSecondaryFromMaterial(mcLabel))
							isPrim = 2;
						else    
							continue;
					}
					isprimArray.push_back(isPrim);
					idArray.push_back(idTrack);
					nNchRec++;

				}
			}
			else if(isHy2){
				newTrack = new AliESDtrack(*esdtrack);
				if(esdtrack->GetConstrainedParam()){
					const AliExternalTrackParam* constrainParam = esdtrack->GetConstrainedParam();
					newTrack->Set(constrainParam->GetX(),constrainParam->GetAlpha(),constrainParam->GetParameter(),constrainParam->GetCovariance());
					newTrack->GetImpactParameters(fdcaxy,fdcaz);
					ptArray.push_back(newTrack->Pt());
					phiArray.push_back(newTrack->Phi());
					dcaxyArray.push_back(fdcaxy);
					dcazArray.push_back(fdcaz);
					if(fUseMC){// get label: 0: prim, 1: weak decays, 2: material
						if(fMC->IsPhysicalPrimary(mcLabel))
							isPrim = 0;
						else if (fMC->IsSecondaryFromWeakDecay(mcLabel))
							isPrim = 1;
						else if(fMC->IsSecondaryFromMaterial(mcLabel))
							isPrim = 2;
						else    
							continue;
					}
					isprimArray.push_back(isPrim);
					idArray.push_back(idTrack);
					nNchRec++;

				}
			}

		}
	}
	return nNchRec;

}
//________________________________________________________________________
