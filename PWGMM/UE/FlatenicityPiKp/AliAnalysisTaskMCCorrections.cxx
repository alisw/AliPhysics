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
#include <AliKFVertex.h>
#include <AliKFParticle.h>
using std::cout;
using std::endl;

#include "AliAnalysisTaskMCCorrections.h"

static const int nCent = 9;
static const int nEta = 4;
/* static const int nFlatbins = 8; */
static const double centClass[nCent+1] = {0.0,1.0,5.0,10.0,20.0,30.0,40.0,50.0,70.0,100.0};
/* static const double Flatbins_16kl[nFlatbins+1] = { 0.0, 0.102, 0.12, 0.132, 0.151, 0.168, 0.186, 0.205, 1.01 }; */
/* static const double Flatbins_lhc16deghijp[nFlatbins+1] = { 0.0, 0.1, 0.117, 0.129, 0.148, 0.165, 0.183, 0.202, 1.01 }; */
/* static const double Flatbins_lhc18bdefghijklmo[nFlatbins+1] = { 0.0, 0.099, 0.117, 0.129, 0.148, 0.165, 0.183, 0.202, 1.01 }; */
static const char* etaClass[nEta] = {"02","24","46","68"};
static const char* ParticleType[3] = {"Primaries","MaterialInt","WeakDecays"};

using namespace std; // std namespace: so you can do things like 'cout' etc

ClassImp(AliAnalysisTaskMCCorrections) // classimp: necessary for root

AliAnalysisTaskMCCorrections::AliAnalysisTaskMCCorrections()
	: AliAnalysisTaskSE(), fESD(0), fEventCuts(0x0), fMCStack(0), fMC(0),
	fUseMC(kFALSE), fV0MMultiplicity(-1.0),
	fDeltaV0(kTRUE), fRemoveTrivialScaling(kFALSE), fPIDResponse(0x0),
	fTrackFilter(0x0), fTrackFilterPID(0x0), fOutputList(0), fEtaCut(0.8), fPtMin(0.5), fNcl(70), 
	fcutLow(0x0), fcutHigh(0x0), fcutDCAxy(0x0), fSystVarTrkCuts(0),
	fv0mpercentile(0), fFlat(-1), fFlatTPC(-1.), fFlatMC(-1),
	fMultSelection(0x0), fPeriod("16k"),
	hFlatenicityMC(0), hFlatenicityMCRec(0), hFlatResponse(0),
	hActivityV0McSect(0), hFlatVsNchMC(0), hTrueINEL_vtx(0), hAccINEL_vtx(0),
	hTrueINELWithFlat_evts(0), hAccINELWithFlat_evts(0),
	hMCPtPionPos(0),hMCPtKaonPos(0),hMCPtProtonPos(0),
	hMCPtPionNeg(0),hMCPtKaonNeg(0),hMCPtProtonNeg(0),
	hTPCRecTracksPionPos(0), hTPCRecTracksKaonPos(0), hTPCRecTracksProtonPos(0),
	hTPCRecTracksPionNeg(0), hTPCRecTracksKaonNeg(0), hTPCRecTracksProtonNeg(0),
	hTOFRecTracksPionPos(0), hTOFRecTracksKaonPos(0), hTOFRecTracksProtonPos(0),
	hTOFRecTracksPionNeg(0), hTOFRecTracksKaonNeg(0), hTOFRecTracksProtonNeg(0),
	hrTPCRecTracksPion(0), hrTPCRecTracksKaon(0), hrTPCRecTracksProton(0)
{

	for (int i_eta = 0; i_eta < nEta; ++i_eta)
	{
		nsigma_kaon_h[i_eta] = 0;
		random_cont_in_kaon_h[i_eta] = 0;
		nsigma_proton_h[i_eta] = 0;
		random_cont_in_proton_h[i_eta] = 0;
		nsigma_pion_h[i_eta] = 0;
		random_cont_in_pion_h[i_eta] = 0;
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

	for(Int_t i = 0; i < 4; ++i){
		hTrueINELWithFlat_pT[i] = 0;
		hAccINELWithFlat_pT[i] = 0;
	}

}
//_____________________________________________________________________________
AliAnalysisTaskMCCorrections::AliAnalysisTaskMCCorrections(const char *name)
	: AliAnalysisTaskSE(name), fESD(0), fEventCuts(0x0), fMCStack(0), fMC(0),
	fUseMC(kFALSE), fV0MMultiplicity(-1.0),
	fDeltaV0(kTRUE), fRemoveTrivialScaling(kFALSE), fPIDResponse(0x0),
	fTrackFilter(0x0), fTrackFilterPID(0x0), fOutputList(0), fEtaCut(0.8), fPtMin(0.5), fNcl(70),
	fcutLow(0x0), fcutHigh(0x0), fcutDCAxy(0x0), fSystVarTrkCuts(0),
	fv0mpercentile(0), fFlat(-1), fFlatTPC(-1.), fFlatMC(-1),
	fMultSelection(0x0), fPeriod("16k"),
	hFlatenicityMC(0), hFlatenicityMCRec(0), hFlatResponse(0),
	hActivityV0McSect(0), hFlatVsNchMC(0), hTrueINEL_vtx(0), hAccINEL_vtx(0),
	hTrueINELWithFlat_evts(0), hAccINELWithFlat_evts(0),
	hMCPtPionPos(0),hMCPtKaonPos(0),hMCPtProtonPos(0),
	hMCPtPionNeg(0),hMCPtKaonNeg(0),hMCPtProtonNeg(0),
	hTPCRecTracksPionPos(0), hTPCRecTracksKaonPos(0), hTPCRecTracksProtonPos(0),
	hTPCRecTracksPionNeg(0), hTPCRecTracksKaonNeg(0), hTPCRecTracksProtonNeg(0),
	hTOFRecTracksPionPos(0), hTOFRecTracksKaonPos(0), hTOFRecTracksProtonPos(0),
	hTOFRecTracksPionNeg(0), hTOFRecTracksKaonNeg(0), hTOFRecTracksProtonNeg(0),
	hrTPCRecTracksPion(0), hrTPCRecTracksKaon(0), hrTPCRecTracksProton(0)
{
	for (int i_eta = 0; i_eta < nEta; ++i_eta)
	{
		nsigma_kaon_h[i_eta] = 0;
		random_cont_in_kaon_h[i_eta] = 0;
		nsigma_proton_h[i_eta] = 0;
		random_cont_in_proton_h[i_eta] = 0;
		nsigma_pion_h[i_eta] = 0;
		random_cont_in_pion_h[i_eta] = 0;
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

	for(Int_t i = 0; i < 4; ++i){
		hTrueINELWithFlat_pT[i] = 0;
		hAccINELWithFlat_pT[i] = 0;
	}

	DefineInput(0, TChain::Class()); // define the input of the analysis: in this
					 // case you take a 'chain' of events
	DefineOutput(1, TList::Class()); // define the ouptut of the analysis: in this
					 // case it's a list of histograms
}
//_____________________________________________________________________________
AliAnalysisTaskMCCorrections::~AliAnalysisTaskMCCorrections() {
	// destructor
	if (fOutputList) {
		delete fOutputList; // at the end of your task, it is deleted from memory by
				    // calling this function
		fOutputList = 0x0;
	}
}
//_____________________________________________________________________________
void AliAnalysisTaskMCCorrections::UserCreateOutputObjects() {

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

	fcutLow = new TF1("StandardPhiCutLow",  "0.1/x/x+pi/18.0-0.025", 0, 50);
	fcutHigh = new TF1("StandardPhiCutHigh", "0.12/x+pi/18.0+0.035", 0, 50);

	fcutDCAxy = new TF1("fMaxDCAxy","[0]+[1]/(x^[2])",0,1e10);
	fcutDCAxy->SetParameter(0,0.0105);
	fcutDCAxy->SetParameter(1,0.0350);
	fcutDCAxy->SetParameter(2,1.1);

	const int nPtbins = 41;
	double Ptbins[nPtbins+1] = {
		0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7,
		0.75, 0.8, 0.85, 0.9, 0.95, 1.0,  1.1, 1.2,  1.3, 1.4,
		1.5, 1.6,  1.7, 1.8,  1.9, 2.0,  2.2, 2.4,  2.6, 2.8,
		3.0, 3.2,  3.4, 3.6,  3.8, 4.0, 4.5, 5.0, 6.0,
		8.0, 10.0, 20.0 };

	// create output objects
	float min_flat = -0.01;
	float max_flat = 1.01;
	int nbins_flat = 1020;
	if (fRemoveTrivialScaling) {
		min_flat = -0.1;
		max_flat = 9.9;
		nbins_flat = 2000;
	}

	const int nnSigmabins = 80;
	double nSigmabins[nnSigmabins+1] = {0.0};

	for (int i = 0; i <= nnSigmabins; ++i){
		nSigmabins[i] = -4.0 + i*0.1;
	}

	const int nDCAbins = 200;
	double DCAbins[nDCAbins+1] = {0.0};

	for (int i = 0; i <= nDCAbins; ++i){
		DCAbins[i] = -3.0 + i*0.03;
	}

	/* double Flatbins[nFlatbins + 1]; */
	/* if (fPeriod=="16k" || fPeriod=="16l") { */
	/* 	for (int i = 0; i <= nFlatbins; ++i) */ 
	/* 	{ */
	/* 		Flatbins[i] = Flatbins_16kl[i]; */
	/* 	} */
	/* } */
	/* else if (fPeriod=="16d" || fPeriod=="16e" || fPeriod=="16g" || fPeriod=="16h" || fPeriod=="16i" || fPeriod=="16j" || fPeriod=="16o" || fPeriod=="16p") { */
	/* 	for (int i = 0; i <= nFlatbins; ++i) */ 
	/* 	{ */
	/* 		Flatbins[i] = Flatbins_lhc16deghijp[i]; */
	/* 	} */
	/* } */
	/* else{ */
	/* 	for (int i = 0; i <= nFlatbins; ++i) */ 
	/* 	{ */
	/* 		Flatbins[i] = Flatbins_lhc18bdefghijklmo[i]; */
	/* 	} */
	/* } */

	const int nFlatbins = 1020;
	double Flatbins[nFlatbins+1] = {0.0};
	for (int i = 0; i <= nFlatbins; ++i) {
		Flatbins[i] = -0.01 + (double)i * 0.001;
	}

	OpenFile(1);
	fOutputList = new TList(); // this is a list which will contain all of your histograms
	fOutputList->SetOwner(kTRUE); // memory stuff: the list is owner of all

	hFlatenicityMC = new TH2F("hFlatenicityMC", ";True Flatenicity; V0M Percentile;", nbins_flat, min_flat, max_flat, nCent, centClass );
	fOutputList->Add(hFlatenicityMC);

	hFlatenicityMCRec = new TH2F("hFlatenicityMCRec",";rec Flatenicity;V0M Percentile;", nbins_flat, min_flat, max_flat, nCent, centClass );
	fOutputList->Add(hFlatenicityMCRec);

	hFlatResponse = new TH2F("hFlatResponse", "; true flat; measured flat", nbins_flat, min_flat, max_flat, nbins_flat, min_flat, max_flat);
	fOutputList->Add(hFlatResponse);

	hFlatVsNchMC = new TH2F("hFlatVsNchMC", "; true flat; true Nch", nbins_flat, min_flat, max_flat, 100, -0.5, 99.5);
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

	for(int i = 0; i < 3; ++i)
	{
		hPionTOFDCAxyNeg[i] = new TH3F(Form("hPion%sTOFDCAxyNeg",ParticleType[i]),"; #it{p}_{T} (GeV/#it{c}); DCA_{xy}; Flattenicity",nPtbins,Ptbins,nDCAbins,DCAbins,nFlatbins,Flatbins);
		fOutputList->Add(hPionTOFDCAxyNeg[i]);
		hProtonTOFDCAxyNeg[i] = new TH3F(Form("hProton%sTOFDCAxyNeg",ParticleType[i]),"; #it{p}_{T} (GeV/#it{c}); DCA_{xy}; Flattenicity",nPtbins,Ptbins,nDCAbins,DCAbins,nFlatbins,Flatbins);	
		fOutputList->Add(hProtonTOFDCAxyNeg[i]);
		hPionTOFDCAxyPos[i] = new TH3F(Form("hPion%sTOFDCAxyPos",ParticleType[i]),"; #it{p}_{T} (GeV/#it{c}); DCA_{xy}; Flattenicity",nPtbins,Ptbins,nDCAbins,DCAbins,nFlatbins,Flatbins);
		fOutputList->Add(hPionTOFDCAxyPos[i]);
		hProtonTOFDCAxyPos[i] = new TH3F(Form("hProton%sTOFDCAxyPos",ParticleType[i]),"; #it{p}_{T} (GeV/#it{c}); DCA_{xy}; Flattenicity",nPtbins,Ptbins,nDCAbins,DCAbins,nFlatbins,Flatbins);	
		fOutputList->Add(hProtonTOFDCAxyPos[i]);

		hPionTPCDCAxyNeg[i] = new TH3F(Form("hPion%sTPCDCAxyNeg",ParticleType[i]),"; #it{p}_{T} (GeV/#it{c}); DCA_{xy}; Flattenicity",nPtbins,Ptbins,nDCAbins,DCAbins,nFlatbins,Flatbins);	
		fOutputList->Add(hPionTPCDCAxyNeg[i]);
		hProtonTPCDCAxyNeg[i] = new TH3F(Form("hProton%sTPCDCAxyNeg",ParticleType[i]),"; #it{p}_{T} (GeV/#it{c}); DCA_{xy}; Flattenicity",nPtbins,Ptbins,nDCAbins,DCAbins,nFlatbins,Flatbins);	
		fOutputList->Add(hProtonTPCDCAxyNeg[i]);
		hPionTPCDCAxyPos[i] = new TH3F(Form("hPion%sTPCDCAxyPos",ParticleType[i]),"; #it{p}_{T} (GeV/#it{c}); DCA_{xy}; Flattenicity",nPtbins,Ptbins,nDCAbins,DCAbins,nFlatbins,Flatbins);	
		fOutputList->Add(hPionTPCDCAxyPos[i]);
		hProtonTPCDCAxyPos[i] = new TH3F(Form("hProton%sTPCDCAxyPos",ParticleType[i]),"; #it{p}_{T} (GeV/#it{c}); DCA_{xy}; Flattenicity",nPtbins,Ptbins,nDCAbins,DCAbins,nFlatbins,Flatbins);	
		fOutputList->Add(hProtonTPCDCAxyPos[i]);
	}

	for (int i_eta = 0; i_eta < nEta; ++i_eta) 
	{	
		nsigma_kaon_h[i_eta] = new TH2F(Form("nsigma_kaon_h_%s",etaClass[i_eta]),"; #it{p}_{T} (GeV/#it{c}); n#sigma",nPtbins,Ptbins,nnSigmabins,nSigmabins);
		fOutputList->Add(nsigma_kaon_h[i_eta]);
		random_cont_in_kaon_h[i_eta] = new TH2F(Form("random_cont_in_kaon_h_%s",etaClass[i_eta]),"; #it{p}_{T} (GeV/#it{c}); n#sigma",nPtbins,Ptbins,nnSigmabins,nSigmabins);
		fOutputList->Add(random_cont_in_kaon_h[i_eta]);

		nsigma_proton_h[i_eta] = new TH2F(Form("nsigma_proton_h_%s",etaClass[i_eta]),"; #it{p}_{T} (GeV/#it{c}); n#sigma",nPtbins,Ptbins,nnSigmabins,nSigmabins);
		fOutputList->Add(nsigma_proton_h[i_eta]);
		random_cont_in_proton_h[i_eta] = new TH2F(Form("random_cont_in_proton_h_%s",etaClass[i_eta]),"; #it{p}_{T} (GeV/#it{c}); n#sigma",nPtbins,Ptbins,nnSigmabins,nSigmabins);
		fOutputList->Add(random_cont_in_proton_h[i_eta]);

		nsigma_pion_h[i_eta] = new TH2F(Form("nsigma_pion_h_%s",etaClass[i_eta]),"; #it{p}_{T} (GeV/#it{c}); n#sigma",nPtbins,Ptbins,nnSigmabins,nSigmabins);
		fOutputList->Add(nsigma_pion_h[i_eta]);
		random_cont_in_pion_h[i_eta] = new TH2F(Form("random_cont_in_pion_h_%s",etaClass[i_eta]),"; #it{p}_{T} (GeV/#it{c}); n#sigma",nPtbins,Ptbins,nnSigmabins,nSigmabins);
		fOutputList->Add(random_cont_in_pion_h[i_eta]);
	}

	hTrueINEL_vtx = new TH1F("hTrueINEL_VtxZ","; Vertex Z (cm)",200, -12.0,12.0);
	fOutputList->Add(hTrueINEL_vtx);
	hAccINEL_vtx = new TH1F("hAccINEL_VtxZ","; Vertex Z (cm)",200, -12.0,12.0);
	fOutputList->Add(hAccINEL_vtx);
	hTrueINELWithFlat_evts = new TH2F("hTrueINEL_vs_V0M_Flat_evts","; Flatenicity; V0M; Counts",nFlatbins,Flatbins,nCent,centClass);
	fOutputList->Add(hTrueINELWithFlat_evts);
	hAccINELWithFlat_evts = new TH2F("hAccINEL_vs_V0M_Flat_evts","; Flatenicity; V0M; Counts",nFlatbins,Flatbins,nCent,centClass);
	fOutputList->Add(hAccINELWithFlat_evts);

	const char* PIDnames[] = {"Charged","Pion","Kaon","Proton"};
	for(int i = 0; i < 4; ++i)
	{
		hTrueINELWithFlat_pT[i] = new TH3F(Form("hTrueINEL_vs_V0M_Flat_%s",PIDnames[i]),"; Flatenicity; V0M; #it{p}_{T} (GeV/#it{c})",nFlatbins,Flatbins,nCent,centClass,nPtbins,Ptbins);
		fOutputList->Add(hTrueINELWithFlat_pT[i]);
		hAccINELWithFlat_pT[i] = new TH3F(Form("hAccINEL_vs_V0M_Flat_%s",PIDnames[i]),"; Flatenicity; V0M; #it{p}_{T} (GeV/#it{c})",nFlatbins,Flatbins,nCent,centClass,nPtbins,Ptbins);
		fOutputList->Add(hAccINELWithFlat_pT[i]);
	}

	hActivityV0McSect = new TProfile("hActivityV0McSect", "true; V0 sector; #LTmultiplicity#GT", 64, -0.5, 63.5);
	fOutputList->Add(hActivityV0McSect);

	/* fEventCuts.AddQAplotsToList(fOutputList); */
	PostData(1, fOutputList); // postdata will notify the analysis manager of
				  // changes / updates to the
}
//_____________________________________________________________________________
void AliAnalysisTaskMCCorrections::UserExec(Option_t *) {

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
		if (TMath::Abs(vtxMC[2]) <= 10) { isGoodVtxPosMC = kTRUE; }

		if (isGoodVtxPosMC) { hTrueINEL_vtx->Fill(vtxMC[2]); }
	}

	// Multiplicity Estimation
	fMultSelection = (AliMultSelection *)fESD->FindListObject("MultSelection");
	if (!fMultSelection) { std::cout << "------- No AliMultSelection Object Found --------" << fMultSelection << '\n'; }
	fv0mpercentile = -999;
	fv0mpercentile = fMultSelection->GetMultiplicityPercentile("V0M");

	// 
	fFlat = GetFlatenicityV0();

	// True INEL > 0 
	// If this condition is false: |z|<= 10 cm
	// then the event is ignored
	if (isGoodVtxPosMC) { TrueINEL(); }
	else { return; }

	// Trigger selection
	UInt_t fSelectMask = fInputHandler->IsEventSelected();
	Bool_t isINT7selected = fSelectMask & AliVEvent::kINT7;
	if (!isINT7selected) { return; }

	// Good events
	if (!fEventCuts.AcceptEvent(event)) {
		PostData(1, fOutputList);
		return;
	}

	// Good vertex
	Bool_t hasRecVertex = kFALSE;
	hasRecVertex = HasRecVertex();
	if (!hasRecVertex) { return; }

	/* const AliVVertex *vtTrc = fMC->GetPrimaryVertex(); */
	/* if (!vtTrc) { return; } */
	/* double mc_vtxZ_pos = vtTrc->GetZ(); */
	const AliESDVertex *vt_rec = fESD->GetPrimaryVertex();
	if (!vt_rec) { return; }
	double rec_vtxZ_pos = vt_rec->GetZ();
	hAccINEL_vtx->Fill(rec_vtxZ_pos);

	// Measured pT spectra with Accepted INEL > 0
	AccINEL();

	/* fMidRapidityMult = GetMidRapidityMultiplicity(); */
	/* fFlatTPC = GetFlatenicityTPC(); */ 
	fFlatMC = -1;
	if (fUseMC) {
		if (!isGoodVtxPosMC) { return; }
		MakeMCanalysisPID();
		nSigmaContamination();
		fFlatMC = GetFlatenicityMC();
		hFlatenicityMC->Fill(fFlatMC,fv0mpercentile);
		hFlatenicityMCRec->Fill(fFlat,fv0mpercentile);
		hFlatResponse->Fill(fFlatMC, fFlat);
	}

	PostData(1, fOutputList); // stream the result of this event to the output
				  // manager which will write it to a file
}
//______________________________________________________________________________
void AliAnalysisTaskMCCorrections::Terminate(Option_t *) {}
//______________________________________________________________________________
void AliAnalysisTaskMCCorrections::TrueINEL() {

	int particles = 0;
	for (Int_t i = 0; i < fMC->GetNumberOfTracks(); ++i) 
	{
		AliMCParticle *particle = (AliMCParticle *)fMC->GetTrack(i);
		if (!particle)
			continue;
		if (!fMC->IsPhysicalPrimary(i))
			continue;
		if (TMath::Abs(particle->Eta()) > 1.0)
			continue;
		if (particle->Pt() <= 0.0)
			continue;
		if (TMath::Abs(particle->Charge()) < 0.1)
			continue;

		if (AliAnalysisUtils::IsParticleFromOutOfBunchPileupCollision(i,fMC))
			continue;

		particles++;
	}

	// This is the INEL > 0 condition
	if (particles < 1) { return; }
	/* hTrueINEL_evts->Fill(fv0mpercentile); */
	hTrueINELWithFlat_evts->Fill(fFlat,fv0mpercentile);

	for (Int_t i = 0; i < fMC->GetNumberOfTracks(); ++i) 
	{

		AliMCParticle *particle = (AliMCParticle *)fMC->GetTrack(i);
		if (!particle)
			continue;
		if (!fMC->IsPhysicalPrimary(i))
			continue;
		if (TMath::Abs(particle->Eta()) > fEtaCut)
			continue;
		if (particle->Pt() < fPtMin)
			continue;
		// To allow K0s and Lambda
		/* if (TMath::Abs(particle->Charge()) < 0.1) */
		/* 	continue; */

		if (AliAnalysisUtils::IsParticleFromOutOfBunchPileupCollision(i,fMC))
			continue;

		int pdgCode = particle->PdgCode();
		int pidCode = GetPidCode(pdgCode);
		double pt = particle->Pt();

		// K0s and Lambda
		/* if ((TMath::Abs(particle->Charge()) < 0.1) && (pidCode==4 || pidCode==5)) { */ 
		/* 	hTrueINELWithFlat_pT[pidCode]->Fill(fFlat,fv0mpercentile,pt); */ 
		/* } */
		if (TMath::Abs(particle->Charge()) < 0.1) { continue; }
		hTrueINELWithFlat_pT[0]->Fill(fFlat,fv0mpercentile,pt);
		// Pions: pidCode = 1
		// Kaons: pidCode = 2
		// Protons: pidCode = 3
		if (pidCode > 3) { continue; }
		hTrueINELWithFlat_pT[pidCode]->Fill(fFlat,fv0mpercentile,pt);
	}
}
//______________________________________________________________________________
void AliAnalysisTaskMCCorrections::AccINEL() {

	int particles = 0;
	for (Int_t i = 0; i < fMC->GetNumberOfTracks(); ++i) 
	{
		AliMCParticle *particle = (AliMCParticle *)fMC->GetTrack(i);
		if (!particle)
			continue;
		if (!fMC->IsPhysicalPrimary(i))
			continue;
		if (TMath::Abs(particle->Eta()) > 1.0)
			continue;
		if (particle->Pt() < fPtMin)
			continue;
		if (TMath::Abs(particle->Charge()) < 0.1)
			continue;

		if (AliAnalysisUtils::IsParticleFromOutOfBunchPileupCollision(i,fMC))
			continue;

		particles++;
	}

	// INEL > 0 condition
	if (particles < 1) { return; }
	/* hAccINEL_evts->Fill(fv0mpercentile); */
	hAccINELWithFlat_evts->Fill(fFlat,fv0mpercentile);

	for (Int_t i = 0; i < fMC->GetNumberOfTracks(); ++i) 
	{

		AliMCParticle *particle = (AliMCParticle *)fMC->GetTrack(i);
		if (!particle)
			continue;
		if (!fMC->IsPhysicalPrimary(i))
			continue;
		if (TMath::Abs(particle->Eta()) > fEtaCut)
			continue;
		if (particle->Pt() < fPtMin)
			continue;
		// To allow K0s and Lambda
		/* if (TMath::Abs(particle->Charge()) < 0.1) */
		/* 	continue; */

		if (AliAnalysisUtils::IsParticleFromOutOfBunchPileupCollision(i,fMC))
			continue;

		int pdgCode = particle->PdgCode();
		int pidCode = GetPidCode(pdgCode);
		double pt = particle->Pt();

		// K0s and Lambda
		/* if ((TMath::Abs(particle->Charge()) < 0.1) && (pidCode==4 || pidCode==5)) { */ 
		/* 	hAccINELWithFlat_pT[pidCode]->Fill(fFlat,fv0mpercentile,pt); */ 
		/* } */
		if (TMath::Abs(particle->Charge()) < 0.1) { continue; }
		hAccINELWithFlat_pT[0]->Fill(fFlat,fv0mpercentile,pt);
		// Pions: pidCode = 1
		// Kaons: pidCode = 2
		// Protons: pidCode = 3
		if (pidCode > 3) { continue; }
		hAccINELWithFlat_pT[pidCode]->Fill(fFlat,fv0mpercentile,pt);
	}
}
//______________________________________________________________________________
void AliAnalysisTaskMCCorrections::MakeMCanalysisPID() {

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
		if (AliAnalysisUtils::IsParticleFromOutOfBunchPileupCollision(i,fMC))
			continue;

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
					if (pidCode==1) hPionTOFDCAxyNeg[0]->Fill(pt,DCAxy,fFlat);	
					if (pidCode==3) hProtonTOFDCAxyNeg[0]->Fill(pt,DCAxy,fFlat);	
				}else{ 
					if (pidCode==1) hPionTOFDCAxyPos[0]->Fill(pt,DCAxy,fFlat);	
					if (pidCode==3) hProtonTOFDCAxyPos[0]->Fill(pt,DCAxy,fFlat);	
				}
			}

			if(esdTrack->Charge() < 0.0){
				if (pidCode==1) hPionTPCDCAxyNeg[0]->Fill(pt,DCAxy,fFlat);	
				if (pidCode==3) hProtonTPCDCAxyNeg[0]->Fill(pt,DCAxy,fFlat);	
			}else{
				if (pidCode==1) hPionTPCDCAxyPos[0]->Fill(pt,DCAxy,fFlat);	
				if (pidCode==3) hProtonTPCDCAxyPos[0]->Fill(pt,DCAxy,fFlat);	
			}
		}

		if( fMC->IsSecondaryFromMaterial(mcLabel) ){

			if( TOFPID(esdTrack) ){

				if(esdTrack->Charge() < 0.0){
					if (pidCode==1) hPionTOFDCAxyNeg[1]->Fill(pt,DCAxy,fFlat);	
					if (pidCode==3) hProtonTOFDCAxyNeg[1]->Fill(pt,DCAxy,fFlat);	
				}else{ 
					if (pidCode==1) hPionTOFDCAxyPos[1]->Fill(pt,DCAxy,fFlat);	
					if (pidCode==3) hProtonTOFDCAxyPos[1]->Fill(pt,DCAxy,fFlat);	
				}
			}

			if(esdTrack->Charge() < 0.0){
				if (pidCode==1) hPionTPCDCAxyNeg[1]->Fill(pt,DCAxy,fFlat);	
				if (pidCode==3) hProtonTPCDCAxyNeg[1]->Fill(pt,DCAxy,fFlat);	
			}else{
				if (pidCode==1) hPionTPCDCAxyPos[1]->Fill(pt,DCAxy,fFlat);	
				if (pidCode==3) hProtonTPCDCAxyPos[1]->Fill(pt,DCAxy,fFlat);	
			}
		}

		if( fMC->IsSecondaryFromWeakDecay(mcLabel) ){

			if( TOFPID(esdTrack) ){

				if(esdTrack->Charge() < 0.0){
					if (pidCode==1) hPionTOFDCAxyNeg[2]->Fill(pt,DCAxy,fFlat);	
					if (pidCode==3) hProtonTOFDCAxyNeg[2]->Fill(pt,DCAxy,fFlat);	
				}else{ 
					if (pidCode==1) hPionTOFDCAxyPos[2]->Fill(pt,DCAxy,fFlat);	
					if (pidCode==3) hProtonTOFDCAxyPos[2]->Fill(pt,DCAxy,fFlat);	
				}
			}

			if(esdTrack->Charge() < 0.0){
				if (pidCode==1) hPionTPCDCAxyNeg[2]->Fill(pt,DCAxy,fFlat);	
				if (pidCode==3) hProtonTPCDCAxyNeg[2]->Fill(pt,DCAxy,fFlat);	
			}else{
				if (pidCode==1) hPionTPCDCAxyPos[2]->Fill(pt,DCAxy,fFlat);	
				if (pidCode==3) hProtonTPCDCAxyPos[2]->Fill(pt,DCAxy,fFlat);	
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
Double_t AliAnalysisTaskMCCorrections::GetFlatenicityV0() {

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
Double_t AliAnalysisTaskMCCorrections::GetMidRapidityMultiplicity() {

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
Double_t AliAnalysisTaskMCCorrections::GetFlatenicityMC() {

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
		if (AliAnalysisUtils::IsParticleFromOutOfBunchPileupCollision(i,fMC))
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
			float mult = RhoLattice[i_segment];
			RhoLattice[i_segment] /= deltaEta;
			// Filling histos with mult info
			if (fDeltaV0) { mult /= deltaEta; }
			else { mult *= 1.0; }

			hActivityV0McSect->Fill(i_segment,mult);
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
	if (mRho > 0.0) {
		if (fRemoveTrivialScaling) {
			flatenicity = TMath::Sqrt(1.0 * nMult) * sRho / mRho;
		} else {
			flatenicity = sRho / mRho;
		}
		hFlatVsNchMC->Fill(flatenicity, nMult);
	} else {
		flatenicity = 999;
	}
	return flatenicity;
}
//______________________________________________________________________________
Bool_t AliAnalysisTaskMCCorrections::HasRecVertex() {

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
Bool_t AliAnalysisTaskMCCorrections::TOFPID(AliESDtrack * track) 
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
void AliAnalysisTaskMCCorrections::nSigmaContamination() {

	Int_t iTracks(fESD->GetNumberOfTracks());          

	for (Int_t iT = 0; iT < iTracks; ++iT) 
	{

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

		/* if (!fMC->IsPhysicalPrimary(mcLabel) ) */
		/* 	continue; */

		AliMCParticle *mc_particle = nullptr;
		mc_particle = (AliMCParticle*)fMC->GetTrack(mcLabel);

		if (!mc_particle) 
			continue;

		if (TMath::Abs(esdTrack->Charge())==0)
			continue;

		Int_t pdgCode = mc_particle->PdgCode();
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

		//	if (TMath::Abs(nsigma_k) <= 4.0){
		if (pidCode==2) { nsigma_kaon_h[nh]->Fill(pt,nsigma_k); }
		if (!(pidCode==2)) { random_cont_in_kaon_h[nh]->Fill(pt,nsigma_k); }
		/* else if (pidCode==1) { pion_cont_in_kaon_h[nh]->Fill(pt,nsigma_k); } */
		/* else if (pidCode==7) { electron_cont_in_kaon_h[nh]->Fill(pt,nsigma_k); } */
		/* else { random_cont_in_kaon_h[nh]->Fill(pt,nsigma_k); } */
		//	}

		//	if (TMath::Abs(nsigma_p) <= 4.0){	
		if (pidCode==3) { nsigma_proton_h[nh]->Fill(pt,nsigma_p); }
		if (!(pidCode==3)) { random_cont_in_proton_h[nh]->Fill(pt,nsigma_p); }
		/* else if (pidCode==1) { pion_cont_in_proton_h[nh]->Fill(pt,nsigma_p); } */
		/* else if (pidCode==7) { electron_cont_in_proton_h[nh]->Fill(pt,nsigma_p); } */
		/* else { random_cont_in_proton_h[nh]->Fill(pt,nsigma_p); } */
		/* //	} */

		//	if (TMath::Abs(nsigma_pi) <= 3.0){	
		if (pidCode==1) { nsigma_pion_h[nh]->Fill(pt,nsigma_pi); }
		if (!(pidCode==1)) { random_cont_in_pion_h[nh]->Fill(pt,nsigma_pi); }
		/* else if (pidCode==2) { kaon_cont_in_pion_h[nh]->Fill(pt,nsigma_pi); } */
		/* else if (pidCode==7) { electron_cont_in_pion_h[nh]->Fill(pt,nsigma_pi); } */
		/* else { random_cont_in_pion_h[nh]->Fill(pt,nsigma_pi); } */
		/* //	} */
	} 

}
//______________________________________________________________________________
Bool_t AliAnalysisTaskMCCorrections::PhiCut(Double_t pt, Double_t phi, Double_t q, Float_t mag, TF1* phiCutLow, TF1* phiCutHigh) {

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
Int_t AliAnalysisTaskMCCorrections::GetPidCode(Int_t pdgCode) {
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


