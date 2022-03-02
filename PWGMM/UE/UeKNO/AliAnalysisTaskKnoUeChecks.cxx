/**************************************************************************
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
 * Author: Feng Fan (feng.fan@cern.ch)                                    *
 *         Antonio Ortiz (antonio.ortiz@nucleares.unam.mx)                *
 **************************************************************************/

/* AliAnalysisTaskKnoUeChecks source code
 * The analysis task produce all the histos needed for MC closure test studies
 * Results include both KNO properties and traditional UE analysis
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

#include "AliAnalysisTaskKnoUeChecks.h"


TF1* f_Eff3;// efficiency for TPC only track

const Char_t * nameReg3[3]={"NS","AS","TS"};
const Int_t ptNbins = 36;
Double_t ptbins3[ptNbins+1] = {
    0.0,  0.1,  0.15,  0.2,  0.25,  0.3,  0.35,  0.4,  0.45,  0.5,  0.6,  0.7,  0.8,  0.9,  1.0,  1.5,  2.0,  2.5,  3.0,  3.5,  4.0,  4.5, 5.0, 6.0, 7.0,  8.0,  9.0,  10.0,  12.0,  14.0,  16.0,  18.0,  20.0,  25.0,  30.0,  40.0,  50.0 };

const Double_t pi = 3.1415926535897932384626433832795028841971693993751058209749445;

class AliAnalysisTaskKnoUeChecks;    // your analysis class

using namespace std;                 // std namespace: so you can do things like 'cout' etc

ClassImp(AliAnalysisTaskKnoUeChecks) // classimp: necessary for root

AliAnalysisTaskKnoUeChecks::AliAnalysisTaskKnoUeChecks() : AliAnalysisTaskSE(),
    fESD(0), fMC(0), fMCStack(0), fUseMC(kFALSE), fIsMCclosure(kFALSE), fLeadingTrackFilter(0x0), fTrackFilter(0x0), fOutputList(0), 
    fEtaCut(0.8), fPtMin(0.5), fLeadPtCutMin(5.0), fLeadPtCutMax(40.0), fGenLeadPhi(0), fGenLeadPt(0), fGenLeadIn(0), fRecLeadPhi(0),
    fRecLeadPt(0), fRecLeadIn(0), hCounter(0), hZvtxAllMeasured(0), hZvtxTrigMeasured(0), hZvtxGoodVtxMeasured(0),
    hZvtxCutAccMeasured(0), hZvtxAllGen(0), hZvtxCutGen(0), hZvtxTrigGen(0), hZvtxGoodVtxGen(0), hZvtxCutAccGen(0), hNchTSData(0),
    hNchTSminData(0), hNchTSmaxData(0), hPhiData_TS1(0),hPhiData_TS2(0), hNchTSGen(0), hNchTSRec(0), hNchTSResponse(0), hNchTSGenTest(0),
    hNchTSRecTest(0), hPhiGen_TS1(0), hPhiGen_TS2(0), hPhiRec_TS1(0), hPhiRec_TS2(0), hPhiGenTest_TS1(0), hPhiGenTest_TS2(0), hPhiRecTest_TS1(0),
    hPhiRecTest_TS2(0), hNchTSminGen(0), hNchTSminRec(0), hNchTSminResponse(0), hNchTSminGenTest(0), hNchTSminRecTest(0), hNchTSmaxGen(0),
    hNchTSmaxRec(0), hNchTSmaxResponse(0), hNchTSmaxGenTest(0), hNchTSmaxRecTest(0), hPtPrimGen(0), hPtPrimRec(0), hNchTSDataTrue(0), 
    hNchTSDataMeasured(0), hNchTSDataResponse(0), hNchTSDataTrueTest(0), hNchTSDataMeasuredTest(0), hNchTSminDataTrue(0), hNchTSminDataMeasured(0), 
    hNchTSminDataResponse(0), hNchTSminDataTrueTest(0), hNchTSminDataMeasuredTest(0), hNchTSmaxDataTrue(0), hNchTSmaxDataMeasured(0), 
    hNchTSmaxDataResponse(0), hNchTSmaxDataTrueTest(0), hNchTSmaxDataMeasuredTest(0)
{
    for (Int_t i = 0;i < 3; ++i) {
	hPhiGen[i] = 0;
	hPhiRec[i] = 0;
    } 

    // default constructor, don't allocate memory here!  this is used by root for IO purposes, it needs to remain empty
}
//_____________________________________________________________________________
AliAnalysisTaskKnoUeChecks::AliAnalysisTaskKnoUeChecks(const char* name) : AliAnalysisTaskSE(name),
    fESD(0), fMC(0), fMCStack(0), fUseMC(kFALSE), fIsMCclosure(kFALSE), fLeadingTrackFilter(0x0), fTrackFilter(0x0), fOutputList(0), 
    fEtaCut(0.8), fPtMin(0.5), fLeadPtCutMin(5.0), fLeadPtCutMax(40.0), fGenLeadPhi(0), fGenLeadPt(0), fGenLeadIn(0), fRecLeadPhi(0),
    fRecLeadPt(0), fRecLeadIn(0), hCounter(0), hZvtxAllMeasured(0), hZvtxTrigMeasured(0), hZvtxGoodVtxMeasured(0),
    hZvtxCutAccMeasured(0), hZvtxAllGen(0), hZvtxCutGen(0), hZvtxTrigGen(0), hZvtxGoodVtxGen(0), hZvtxCutAccGen(0), hNchTSData(0),
    hNchTSminData(0), hNchTSmaxData(0), hPhiData_TS1(0),hPhiData_TS2(0), hNchTSGen(0), hNchTSRec(0), hNchTSResponse(0), hNchTSGenTest(0),
    hNchTSRecTest(0), hPhiGen_TS1(0), hPhiGen_TS2(0), hPhiRec_TS1(0), hPhiRec_TS2(0), hPhiGenTest_TS1(0), hPhiGenTest_TS2(0), hPhiRecTest_TS1(0),
    hPhiRecTest_TS2(0), hNchTSminGen(0), hNchTSminRec(0), hNchTSminResponse(0), hNchTSminGenTest(0), hNchTSminRecTest(0), hNchTSmaxGen(0),
    hNchTSmaxRec(0), hNchTSmaxResponse(0), hNchTSmaxGenTest(0), hNchTSmaxRecTest(0), hPtPrimGen(0), hPtPrimRec(0), hNchTSDataTrue(0), 
    hNchTSDataMeasured(0), hNchTSDataResponse(0), hNchTSDataTrueTest(0), hNchTSDataMeasuredTest(0), hNchTSminDataTrue(0), hNchTSminDataMeasured(0), 
    hNchTSminDataResponse(0), hNchTSminDataTrueTest(0), hNchTSminDataMeasuredTest(0), hNchTSmaxDataTrue(0), hNchTSmaxDataMeasured(0), 
    hNchTSmaxDataResponse(0), hNchTSmaxDataTrueTest(0), hNchTSmaxDataMeasuredTest(0)
{
    for (Int_t i = 0;i < 3; ++i) {
	hPhiGen[i] = 0;
	hPhiRec[i] = 0;
    }

    // constructor
    DefineInput(0, TChain::Class());    // define the input of the analysis: in this case you take a 'chain' of events
    // this chain is created by the analysis manager, so no need to worry about it, does its work automatically
    DefineOutput(1, TList::Class());    // define the ouptut of the analysis: in this case it's a list of histograms

}
//_____________________________________________________________________________
AliAnalysisTaskKnoUeChecks::~AliAnalysisTaskKnoUeChecks()
{

    // destructor
    if(fOutputList) {
	delete fOutputList; // at the end of your task, it is deleted from memory by calling this function
	fOutputList = 0x0;
    }

}
//_____________________________________________________________________________
void AliAnalysisTaskKnoUeChecks::UserCreateOutputObjects()
{
    f_Eff3 = new TF1("fpara","([0]+[1]*x+[2]*x*x)*(x>=0.5&&x<3.5)+([3]+[4]*x+[5]*x*x)*(x>=3.5&&x<=10.0)+([6])*(x>10.0)",0.5,50);
    f_Eff3->SetParameters(8.74533e-01, 1.26958e-02, -5.38522e-03, 8.04741e-01, 1.67890e-02, -8.71589e-04, 8.90427e-01);

    // fCuts1 *** leading particle ***
    if (!fLeadingTrackFilter) {
	fLeadingTrackFilter = new AliAnalysisFilter("trackFilter2015");
	AliESDtrackCuts * fCuts1 = new AliESDtrackCuts();
	// 1 additional cut
	fCuts1->SetMaxFractionSharedTPCClusters(0.4);
	// 13 cuts from StadardITSTPCTrackCuts2015PbPb
	// fCuts1->SetMinNCrossedRowsTPC(70);
	fCuts1->SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);
	fCuts1->SetCutGeoNcrNcl(3., 130., 1.5, 0.85, 0.7);
	fCuts1->SetMaxChi2PerClusterTPC(4);
	fCuts1->SetAcceptKinkDaughters(kFALSE);
	fCuts1->SetRequireTPCRefit(kTRUE);
	fCuts1->SetRequireITSRefit(kTRUE);
	fCuts1->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kAny);
	fCuts1->SetMaxDCAToVertexXYPtDep("0.0182+0.0350/pt^1.01");
	fCuts1->SetMaxChi2TPCConstrainedGlobal(36);
	fCuts1->SetMaxDCAToVertexZ(2);
	fCuts1->SetDCAToVertex2D(kFALSE);
	fCuts1->SetRequireSigmaToVertex(kFALSE);
	fCuts1->SetMaxChi2PerClusterITS(36);
	fLeadingTrackFilter->AddCuts(fCuts1);
    }

    // fCuts2 ***  multiplicity in the transverse, trans-max and trans-min regions ***
    if (!fTrackFilter) {
	fTrackFilter = new AliAnalysisFilter("trackFilterTPConly");
	// 6 cuts from StandardTPConlyTrackCuts
	AliESDtrackCuts * fCuts2 = AliESDtrackCuts::GetStandardTPCOnlyTrackCuts();
	// 2 additional cuts
	fCuts2->SetRequireTPCRefit(kTRUE);
	fCuts2->SetRequireITSRefit(kTRUE);
	fTrackFilter->AddCuts(fCuts2);
    }




    // create output objects
    OpenFile(1);
    // this is a list which will contain all of your histograms, at the end of the analysis, the contents of this list are written  to the output file
    fOutputList = new TList();
    // memory stuff: the list is owner of all objects and will delete them if requested
    fOutputList->SetOwner(kTRUE);

    hCounter = new TH1I("hCounter","Counter; sel; Nev",4,0,4);
    hCounter->GetXaxis()->SetBinLabel(1, "Total processed Events");
    hCounter->GetXaxis()->SetBinLabel(2, "Events selected with trigger");
    hCounter->GetXaxis()->SetBinLabel(3, "Events selected with good vertex");
    hCounter->GetXaxis()->SetBinLabel(4, "Events selected for analysis");
    fOutputList->Add(hCounter);

    hZvtxAllMeasured = new TH1D("hZvtxAllMeasured", "Zvtx before any selection; Zvtx ; Events", 400, -20, 20);
    fOutputList->Add(hZvtxAllMeasured);
    hZvtxTrigMeasured = new TH1D("hZvtxTrigMeasured", "Zvtx after trigger selection; Zvtx ; Events", 400, -20, 20);
    fOutputList->Add(hZvtxTrigMeasured);
    hZvtxGoodVtxMeasured = new TH1D("hZvtxGoodVtxMeasured", "Zvtx after good vertex selection ; Zvtx ; Events", 400, -20, 20);
    fOutputList->Add(hZvtxGoodVtxMeasured);
    hZvtxCutAccMeasured = new TH1D("hZvtxCutAccMeasured", "Zvtx after accpted event; Zvtx ; Events", 400, -20, 20);
    fOutputList->Add(hZvtxCutAccMeasured);

    if (!fUseMC) {// data
	hNchTSData = new TH1D("hNchTSData","; rec mult; Events", 100, -0.5, 99.5);
	fOutputList->Add(hNchTSData);
	hNchTSminData = new TH1D("hNchTSminData","; rec mult; Events", 100, -0.5, 99.5);
	fOutputList->Add(hNchTSminData);
	hNchTSmaxData = new TH1D("hNchTSmaxData","; rec mult; Events",100, -0.5, 99.5);
	fOutputList->Add(hNchTSmaxData);
	hPhiData_TS1= new TH1D("hPhiData_TS1","; #Delta#varphi; Entries", 64, -TMath::Pi()/2.0, 3.0*TMath::Pi()/2.0);
	fOutputList->Add(hPhiData_TS1);
	hPhiData_TS2= new TH1D("hPhiData_TS2","; #Delta#varphi; Entries", 64, -TMath::Pi()/2.0, 3.0*TMath::Pi()/2.0);
	fOutputList->Add(hPhiData_TS2);

	// data driven test
	hNchTSDataTrue = new TH1D("hNchTSDataTrue","detected N_{ch}; true mult; Events", 100, -0.5, 99.5);
	fOutputList->Add(hNchTSDataTrue);
	hNchTSDataMeasured = new TH1D("hNchTSDataMeasured","N_{acc} selecting by efficiency; measured mult; Events", 100, -0.5, 99.5);
	fOutputList->Add(hNchTSDataMeasured);
	hNchTSDataResponse = new TH2D("hNchTSDataResponse","response matrix; measured mult; true mult", 100, -0.5, 99.5, 100, -0.5, 99.5);
	fOutputList->Add(hNchTSDataResponse);
	hNchTSDataTrueTest = new TH1D("hNchTSDataTrueTest","; true mult; Events", 100, -0.5, 99.5);
	fOutputList->Add(hNchTSDataTrueTest);
	hNchTSDataMeasuredTest = new TH1D("hNchTSDataMeasuredTest","; measured mult; Events", 100, -0.5, 99.5);
	fOutputList->Add(hNchTSDataMeasuredTest);

        hNchTSminDataTrue = new TH1D("hNchTSminDataTrue","detected N_{ch}; true mult; Events", 100, -0.5, 99.5);
	fOutputList->Add(hNchTSminDataTrue);
	hNchTSminDataMeasured = new TH1D("hNchTSminDataMeasured","N_{acc} selecting by efficiency; measured mult; Events", 100, -0.5, 99.5);
	fOutputList->Add(hNchTSminDataMeasured);
	hNchTSminDataResponse = new TH2D("hNchTSminDataResponse","response matrix; measured mult; true mult", 100, -0.5, 99.5, 100, -0.5, 99.5);
	fOutputList->Add(hNchTSminDataResponse);
	hNchTSminDataTrueTest = new TH1D("hNchTSminDataTrueTest","; true mult; Events", 100, -0.5, 99.5);
	fOutputList->Add(hNchTSminDataTrueTest);
	hNchTSminDataMeasuredTest = new TH1D("hNchTSminDataMeasuredTest","; measured mult; Events", 100, -0.5, 99.5);
	fOutputList->Add(hNchTSminDataMeasuredTest);
        
	hNchTSmaxDataTrue = new TH1D("hNchTSmaxDataTrue","detected N_{ch}; true mult; Events", 100, -0.5, 99.5);
	fOutputList->Add(hNchTSmaxDataTrue);
	hNchTSmaxDataMeasured = new TH1D("hNchTSmaxDataMeasured","N_{acc} selecting by efficiency; measured mult; Events", 100, -0.5, 99.5);
	fOutputList->Add(hNchTSmaxDataMeasured);
	hNchTSmaxDataResponse = new TH2D("hNchTSmaxDataResponse","response matrix; measured mult; true mult", 100, -0.5, 99.5, 100, -0.5, 99.5);
	fOutputList->Add(hNchTSmaxDataResponse);
	hNchTSmaxDataTrueTest = new TH1D("hNchTSmaxDataTrueTest","; true mult; Events", 100, -0.5, 99.5);
	fOutputList->Add(hNchTSmaxDataTrueTest);
	hNchTSmaxDataMeasuredTest = new TH1D("hNchTSmaxDataMeasuredTest","; measured mult; Events", 100, -0.5, 99.5);
	fOutputList->Add(hNchTSmaxDataMeasuredTest);

    }
    else {// MC
	hZvtxAllGen = new TH1D("hZvtxAllGen", "Zvtx before any event selection; Zvtx ; Events", 400, -20, 20);
	fOutputList->Add(hZvtxAllGen);
	hZvtxCutGen = new TH1D("hZvtxCutGen", "ZvtxCut; Zvtx after vertex position cut; Events", 400, -20, 20);
	fOutputList->Add(hZvtxCutGen);
	hZvtxTrigGen = new TH1D("hZvtxTrigGen", "Zvtx after Trigger selection; Zvtx ; Events", 400, -20, 20);
	fOutputList->Add(hZvtxTrigGen);
	hZvtxGoodVtxGen = new TH1D("hZvtxGoodVtxGen", "Zvtx after good vertex selection; Zvtx ; Events", 400, -20, 20);
	fOutputList->Add(hZvtxGoodVtxGen);
        hZvtxCutAccGen = new TH1D("hZvtxCutAccGen", "Zvtx after vertex cut and accept events; Zvtx ; Events", 400, -20, 20);
	fOutputList->Add(hZvtxCutAccGen);

	for (Int_t i=0;i<3;++i) {
	    hPhiGen[i]= new TH1D(Form("hPhiGen_%s",nameReg3[i]),"; #Delta#varphi; Entries", 64, -TMath::Pi()/2.0, 3.0*TMath::Pi()/2.0);
	    fOutputList->Add(hPhiGen[i]);

	    hPhiRec[i]= new TH1D(Form("hPhiRec_%s",nameReg3[i]),"; #Delta#varphi; Entries", 64, -TMath::Pi()/2.0, 3.0*TMath::Pi()/2.0);
	    fOutputList->Add(hPhiRec[i]);
	}

	hNchTSGen = new TH1D("hNchTSGen","; gen mult; Events", 100, -0.5, 99.5);
	fOutputList->Add(hNchTSGen);
	hNchTSRec = new TH1D("hNchTSRec","; rec mult; Events", 100, -0.5, 99.5);
	fOutputList->Add(hNchTSRec);
	hNchTSResponse = new TH2D("hNchTSResponse","; rec mult; gen mult", 100, -0.5, 99.5, 100, -0.5, 99.5);
	fOutputList->Add(hNchTSResponse);
	hNchTSGenTest = new TH1D("hNchTSGenTest","; gen mult; Events", 100, -0.5, 99.5);
	fOutputList->Add(hNchTSGenTest);
	hNchTSRecTest = new TH1D("hNchTSRecTest","; rec mult; Events", 100, -0.5, 99.5);
	fOutputList->Add(hNchTSRecTest);

	hPhiGen_TS1 = new TH1D("hPhiGen_TS1","; #Delta#varphi; Entries", 64, -TMath::Pi()/2.0, 3.0*TMath::Pi()/2.0);
	fOutputList->Add(hPhiGen_TS1);
	hPhiGen_TS2 = new TH1D("hPhiGen_TS2","; #Delta#varphi; Entries", 64, -TMath::Pi()/2.0, 3.0*TMath::Pi()/2.0);
	fOutputList->Add(hPhiGen_TS2);
	hPhiRec_TS1 = new TH1D("hPhiRec_TS1","; #Delta#varphi; Entries", 64, -TMath::Pi()/2.0, 3.0*TMath::Pi()/2.0);
	fOutputList->Add(hPhiRec_TS1);
	hPhiRec_TS2 = new TH1D("hPhiRec_TS2","; #Delta#varphi; Entries", 64, -TMath::Pi()/2.0, 3.0*TMath::Pi()/2.0);
	fOutputList->Add(hPhiRec_TS2);
	hPhiGenTest_TS1 = new TH1D("hPhiGenTest_TS1","; #Delta#varphi; Entries", 64, -TMath::Pi()/2.0, 3.0*TMath::Pi()/2.0);
	fOutputList->Add(hPhiGenTest_TS1);
	hPhiGenTest_TS2 = new TH1D("hPhiGenTest_TS2","; #Delta#varphi; Entries", 64, -TMath::Pi()/2.0, 3.0*TMath::Pi()/2.0);
	fOutputList->Add(hPhiGenTest_TS2);
	hPhiRecTest_TS1 = new TH1D("hPhiRecTest_TS1","; #Delta#varphi; Entries", 64, -TMath::Pi()/2.0, 3.0*TMath::Pi()/2.0);
	fOutputList->Add(hPhiRecTest_TS1);
	hPhiRecTest_TS2 = new TH1D("hPhiRecTest_TS2","; #Delta#varphi; Entries", 64, -TMath::Pi()/2.0, 3.0*TMath::Pi()/2.0);
	fOutputList->Add(hPhiRecTest_TS2);

	hNchTSminGen = new TH1D("hNchTSminGen","; gen mult; Events", 100, -0.5, 99.5);
	fOutputList->Add(hNchTSminGen);
	hNchTSminRec = new TH1D("hNchTSminRec","; rec mult; Events", 100, -0.5, 99.5);
	fOutputList->Add(hNchTSminRec);
	hNchTSminResponse = new TH2D("hNchTSminResponse","Detector response; rec mult; gen mult", 100, -0.5, 99.5, 100, -0.5, 99.5);
	fOutputList->Add(hNchTSminResponse);
	hNchTSminGenTest = new TH1D("hNchTSminGenTest","; gen mult; Events", 100, -0.5, 99.5);
	fOutputList->Add(hNchTSminGenTest);
	hNchTSminRecTest = new TH1D("hNchTSminRecTest","; rec mult; Events", 100, -0.5, 99.5);
	fOutputList->Add(hNchTSminRecTest);

	hNchTSmaxGen = new TH1D("hNchTSmaxGen","; gen mult; Events", 100, -0.5, 99.5);
	fOutputList->Add(hNchTSmaxGen);
	hNchTSmaxRec = new TH1D("hNchTSmaxRec","; rec mult; Events", 100, -0.5, 99.5);
	fOutputList->Add(hNchTSmaxRec);
	hNchTSmaxResponse = new TH2D("hNchTSmaxResponse","Detector response; rec mult; gen mult", 100, -0.5, 99.5, 100, -0.5, 99.5);
	fOutputList->Add(hNchTSmaxResponse);
	hNchTSmaxGenTest = new TH1D("hNchTSmaxGenTest","; gen mult; Events", 100, -0.5, 99.5);
	fOutputList->Add(hNchTSmaxGenTest);
	hNchTSmaxRecTest = new TH1D("hNchTSmaxRecTest","; rec mult; Events", 100, -0.5, 99.5);
	fOutputList->Add(hNchTSmaxRecTest);

	// tracking efficiency for Nch in transverse, trans-max and trans-min regions
	hPtPrimGen = new TH1D("hPtPrimGen","p_T prim true; p_T; Entries", ptNbins, ptbins3);
	fOutputList->Add(hPtPrimGen);
	hPtPrimRec = new TH1D("hPtPrimRec","p_T prim rec; p_T; Entries", ptNbins, ptbins3);
	fOutputList->Add(hPtPrimRec);
    }


    fEventCuts.AddQAplotsToList(fOutputList);
    PostData(1, fOutputList);// postdata will notify the analysis manager of changes / updates to the

}
//_____________________________________________________________________________
void AliAnalysisTaskKnoUeChecks::UserExec(Option_t *)
{
    AliVEvent *event = InputEvent();
    if (!event) {
	Error("UserExec", "Could not retrieve event");
	return;
    }
    fESD = dynamic_cast<AliESDEvent*>(event);
    if (!fESD) {
	Printf("%s:%d ESDEvent not found in Input Manager",(char*)__FILE__,__LINE__);
	this->Dump();
	return;
    }


    if (fUseMC) {
	// E S D
	fMC = dynamic_cast<AliMCEvent*>(MCEvent());
	if (!fMC) {
	    Printf("%s:%d MCEvent not found in Input Manager",(char*)__FILE__,__LINE__);
	    this->Dump();
	    return;
	}
	fMCStack = fMC->Stack();
    }

    hCounter->Fill(0);

    
    
    Bool_t isGoodVtxPosMC = kFALSE;
    TArrayF vtxMC(3); // primary vertex  MC
    vtxMC[0]=9999; vtxMC[1]=9999;  vtxMC[2]=9999; //initialize with dummy
    if (fUseMC) {
	AliHeader* headerMC = fMC->Header();
	AliGenEventHeader* genHeader = headerMC->GenEventHeader();
	if (genHeader) {
	    genHeader->PrimaryVertex(vtxMC);
	}
	hZvtxAllGen->Fill(vtxMC[2]);

	if(TMath::Abs(vtxMC[2]) <= 10){
	    isGoodVtxPosMC = kTRUE;
            hZvtxCutGen->Fill(vtxMC[2]);
	} 
	// Before trigger selection
	GetLeadingObject(kTRUE);// leading particle at gen level
    }
    // Before trigger selection
    GetLeadingObject(kFALSE);// leading particle at rec level


    const AliVVertex *fPrimaryVtx  = fESD->GetPrimaryVertex();
    Double_t Zvertex = fPrimaryVtx->GetZ();
    hZvtxAllMeasured->Fill(Zvertex);

  
    // === Trigger selection ===
    UInt_t fSelectMask= fInputHandler->IsEventSelected();
    Bool_t isINT7selected = fSelectMask&AliVEvent::kINT7||fSelectMask&AliVEvent::kMB;
    if (!isINT7selected)  return;
    //
    hCounter->Fill(1);
    hZvtxTrigMeasured->Fill(Zvertex);
    if (fUseMC) hZvtxTrigGen->Fill(vtxMC[2]);

    // === Good vertex ===
    Bool_t hasRecVertex = kFALSE;
    hasRecVertex = HasRecVertex();
    if (!hasRecVertex) return;
    //
    hCounter->Fill(2);
    hZvtxGoodVtxMeasured->Fill(Zvertex);
    if (fUseMC) hZvtxGoodVtxGen->Fill(vtxMC[2]);
    
    // === vertex position cut && Good events ===
    fEventCuts.SetMaxVertexZposition(10.0); //setting vertex position cut implemented in Acceptevent()
    if (!fEventCuts.AcceptEvent(event)) {
	PostData(1, fOutputList);
	return;
    }
    //
    hCounter->Fill(3);
    hZvtxCutAccMeasured->Fill(Zvertex);
    if (isGoodVtxPosMC) hZvtxCutAccGen->Fill(vtxMC[2]);



    if (fIsMCclosure) {
	Double_t randomUE = gRandom->Uniform(0.0, 1.0);
	if (randomUE < 0.5) {// corrections (50% stat.)
	    if (isGoodVtxPosMC) {
		// KNO scaling
		if ( (fGenLeadPt>=fLeadPtCutMin && fGenLeadPt<fLeadPtCutMax) && (fRecLeadPt>=fLeadPtCutMin && fRecLeadPt<fLeadPtCutMax) )
		    GetDetectorResponse();
	    }
	}
	else {// for testing the method
	      // KNO scaling
	      if ( (fGenLeadPt>=fLeadPtCutMin && fGenLeadPt<fLeadPtCutMax) && (fRecLeadPt>=fLeadPtCutMin && fRecLeadPt<fLeadPtCutMax) )
		  GetMultiplicityDistributions();
	}
    }
    else {
	if (fUseMC) {// MC
	    if (isGoodVtxPosMC) {
		// KNO scaling
		if ( (fGenLeadPt>=fLeadPtCutMin && fGenLeadPt<fLeadPtCutMax) && (fRecLeadPt>=fLeadPtCutMin && fRecLeadPt<fLeadPtCutMax) )
		    GetDetectorResponse();
		
                // tracking efficiency
		if (fGenLeadPt>=fPtMin && fRecLeadPt>=fPtMin) 
		    GetTrackingEfficiencyTPConly();
	    }
	}
	else {// data
	    Double_t randomUE = gRandom->Uniform(0.0, 1.0);
	    // data driven test
            if (randomUE < 0.5) {// corrections (50% stat.)
	        if (fRecLeadPt>=fLeadPtCutMin && fRecLeadPt<fLeadPtCutMax) {
		    GetDetectorResponseDataDriven();
                    GetMultiplicityDistributionsData();
                }
	    }
	    else {// for testing the method
		if (fRecLeadPt>=fLeadPtCutMin && fRecLeadPt<fLeadPtCutMax)
                    GetMultiplicityDistributionsDataDriven(); 
	    }
	}
    }


    PostData(1, fOutputList); // stream the result of this event to the output manager which will write it to a file
}

//______________________________________________________________________________
void AliAnalysisTaskKnoUeChecks::Terminate(Option_t *)
{

}

//_____________________________________________________________________________
void AliAnalysisTaskKnoUeChecks::GetLeadingObject(Bool_t isMC) 
{
    Double_t flPt = 0;
    Double_t flPhi = 0;
    Int_t flIndex = 0;

    if (isMC) {// loop over generated particles
	for (Int_t i = 0; i < fMC->GetNumberOfTracks(); i++) {
	    AliMCParticle* particle = (AliMCParticle*)fMC->GetTrack(i);
	    if (!particle) continue;

	    if (!fMC->IsPhysicalPrimary(i)) continue;
	    if (particle->Charge() == 0) continue;
	    if (TMath::Abs(particle->Eta()) > fEtaCut) continue;
	    if (particle->Pt() < fPtMin) continue;

	    if (flPt<particle->Pt()) {
		flPt = particle->Pt();
		flPhi = particle->Phi();
		flIndex = i;
	    }
	}

	fGenLeadPhi = flPhi;
	fGenLeadPt  = flPt;
	fGenLeadIn  = flIndex;
    }
    else {// loop over reconstructed tracks
	for (Int_t i=0; i < fESD->GetNumberOfTracks(); i++) {                 
	    AliESDtrack* track = static_cast<AliESDtrack*>(fESD->GetTrack(i)); 
	    if (!track) continue;
	    if (!fLeadingTrackFilter->IsSelected(track)) continue;
	    if (TMath::Abs(track->Eta()) > fEtaCut) continue;
	    if (track->Pt() < fPtMin) continue;

	    if (flPt<track->Pt()){
		flPt  = track->Pt();
		flPhi = track->Phi();
		flIndex = i;
	    }
	} 
	fRecLeadPhi = flPhi;
	fRecLeadPt  = flPt;
	fRecLeadIn  = flIndex;
    }
}

//______________________________________________________________
void AliAnalysisTaskKnoUeChecks::GetTrackingEfficiencyTPConly()
{
    // loop over generated particles
    for (Int_t i = 0; i < fMC->GetNumberOfTracks(); i++) {
	AliMCParticle* particle = (AliMCParticle*)fMC->GetTrack(i);
	if (!particle) continue;
	if (!fMC->IsPhysicalPrimary(i)) continue;
	if (particle->Charge() == 0) continue;
	if (TMath::Abs(particle->Eta()) > fEtaCut) continue;
	if (particle->Pt() < fPtMin) continue;

	hPtPrimGen->Fill(particle->Pt()); // inital pT distribution (MC gen)
    }

    // loop over reconstructed tracks
    for (Int_t i=0; i < fESD->GetNumberOfTracks(); i++) {                
	AliESDtrack* track = static_cast<AliESDtrack*>(fESD->GetTrack(i)); 
	if (!track) continue;
	if (!fTrackFilter->IsSelected(track)) continue; // for KNO analysis we consider TPConly
	if (TMath::Abs(track->Eta()) > fEtaCut) continue;
	if (track->Pt() < fPtMin) continue;

	const Int_t label = TMath::Abs(track->GetLabel());

	if (fMC->IsPhysicalPrimary(label)){
	    hPtPrimRec->Fill(track->Pt());
	}
    }
}

//______________________________________________________________
void AliAnalysisTaskKnoUeChecks::GetDetectorResponseDataDriven()
{
    Int_t mult_true[3]     = {0, 0, 0}; // 0-TS, 1-TSmix, 2-TSmax
    Int_t mult_measured[3] = {0, 0, 0}; // 0-TS, 1-TSmix, 2-TSmax
    Int_t multTS1_true = 0;
    Int_t multTS2_true = 0;
    Int_t multTS1_measured = 0;
    Int_t multTS2_measured = 0;

    
    // loop over reconstructed tracks
    for (Int_t i = 0; i < fESD->GetNumberOfTracks(); i++) {
        if (i == fRecLeadIn) continue;
        
        AliESDtrack* track = static_cast<AliESDtrack*>(fESD->GetTrack(i));
        if (!track) continue;
        if (!fTrackFilter->IsSelected(track)) continue;
        if (TMath::Abs(track->Eta()) > fEtaCut) continue;
        if (track->Pt() < fPtMin) continue;

        Double_t DPhi = DeltaPhi(track->Phi(), fRecLeadPhi);

	if (TMath::Abs(DPhi)< pi/3.0) {// near side
	    continue;
	}
	else if (TMath::Abs(DPhi-pi)< pi/3.0) {// away side
	    continue;
	}
	else {// transverse side
	    mult_true[0]++;
	    if (DPhi>pi/3.0 && DPhi<2*pi/3.0) {
		multTS1_true++;

	    } else {
		multTS2_true++;
	    }
	}
    
        // second track selection following the efficiency
        if (f_Eff3->Eval(track->Pt()) < gRandom->Uniform(0,1)) continue;
        
        if (TMath::Abs(DPhi) < pi/3.0) {// near side
            continue;
        }
        else if (TMath::Abs(DPhi-pi) < pi/3.0) {// away side
            continue;
        }
        else {// transverse side
            mult_measured[0]++;
            
            if (DPhi>pi/3.0 && DPhi<2*pi/3.0) {
                multTS1_measured++;

            } else {
                multTS2_measured++;
            }
        }
    }

    hNchTSDataTrue->Fill(mult_true[0]);
    
    if (multTS2_true >= multTS1_true) {
        mult_true[1] = multTS1_true;
        mult_true[2] = multTS2_true;
    } else {
        mult_true[1] = multTS2_true;
        mult_true[2] = multTS1_true;
    }
    hNchTSminDataTrue->Fill(mult_true[1]);
    hNchTSmaxDataTrue->Fill(mult_true[2]);

    
    hNchTSDataMeasured->Fill(mult_measured[0]);

    if (multTS2_measured >= multTS1_measured) {
        mult_measured[1] = multTS1_measured;
        mult_measured[2] = multTS2_measured;
    } else {
        mult_measured[1] = multTS2_measured;
        mult_measured[2] = multTS1_measured;
    }
    hNchTSminDataMeasured->Fill(mult_measured[1]);
    hNchTSmaxDataMeasured->Fill(mult_measured[2]);

    hNchTSDataResponse->Fill(mult_measured[0], mult_true[0]);
    hNchTSminDataResponse->Fill(mult_measured[1], mult_true[1]);
    hNchTSmaxDataResponse->Fill(mult_measured[2], mult_true[2]);

}

//______________________________________________________________ 
void AliAnalysisTaskKnoUeChecks::GetMultiplicityDistributionsDataDriven()
{
    Int_t mult_true[3]     = {0, 0, 0}; // 0-TS, 1-TSmix, 2-TSmax
    Int_t mult_measured[3] = {0, 0, 0}; // 0-TS, 1-TSmix, 2-TSmax
    Int_t multTS1_true = 0;
    Int_t multTS2_true = 0;
    Int_t multTS1_measured = 0;
    Int_t multTS2_measured = 0;
    

    // loop over reconstructed tracks
    for (Int_t i = 0; i < fESD->GetNumberOfTracks(); i++) {
        if (i == fRecLeadIn) continue;
        
        AliESDtrack* track = static_cast<AliESDtrack*>(fESD->GetTrack(i));
        if (!track) continue;
        if (!fTrackFilter->IsSelected(track)) continue;
        if (TMath::Abs(track->Eta()) > fEtaCut) continue;
        if (track->Pt() < fPtMin) continue;

        
        Double_t DPhi = DeltaPhi(track->Phi(), fRecLeadPhi);
        if (TMath::Abs(DPhi)< pi/3.0) {// near side
            continue;
        }
        else if (TMath::Abs(DPhi-pi)< pi/3.0) {// away side
            continue;
        }
        else {// transverse side
            mult_true[0]++;
            if (DPhi>pi/3.0 && DPhi<2*pi/3.0) {
                multTS1_true++;
                
            } else {
                multTS2_true++;
            }
        }
        
        // second track selection following the efficiency
        if (f_Eff3->Eval(track->Pt()) < gRandom->Uniform(0,1)) continue;
        
        if (TMath::Abs(DPhi) < pi/3.0) {// near side
            continue;
        }
        else if (TMath::Abs(DPhi-pi) < pi/3.0) {// away side
            continue;
        }
        else {// transverse side
            mult_measured[0]++;
            
            if (DPhi>pi/3.0 && DPhi<2*pi/3.0) {
                multTS1_measured++;

            } else {
                multTS2_measured++;
            }
        }
    }

    hNchTSDataTrueTest->Fill(mult_true[0]);

    if (multTS2_true >= multTS1_true) {
        mult_true[1] = multTS1_true;
        mult_true[2] = multTS2_true;
    } else {
        mult_true[1] = multTS2_true;
        mult_true[2] = multTS1_true;
    }
    hNchTSminDataTrueTest->Fill(mult_true[1]);
    hNchTSmaxDataTrueTest->Fill(mult_true[2]);


    hNchTSDataMeasuredTest->Fill(mult_measured[0]);

    if (multTS2_measured >= multTS1_measured) {
        mult_measured[1] = multTS1_measured;
        mult_measured[2] = multTS2_measured;
    } else {
        mult_measured[1] = multTS2_measured;
        mult_measured[2] = multTS1_measured;
    }
    hNchTSminDataMeasuredTest->Fill(mult_measured[1]);
    hNchTSmaxDataMeasuredTest->Fill(mult_measured[2]);
}

//______________________________________________________________
void AliAnalysisTaskKnoUeChecks::GetDetectorResponse()
{
    Int_t multTSgen = 0, multTS1gen = 0, multTS2gen = 0, multTSmin_gen = 0, multTSmax_gen = 0;
    Int_t multTSrec = 0, multTS1rec = 0, multTS2rec = 0, multTSmin_rec = 0, multTSmax_rec = 0;

    // loop over generated particles
    for (Int_t i = 0; i < fMC->GetNumberOfTracks(); i++) {
	if (i==fGenLeadIn) continue;

	AliMCParticle* particle = (AliMCParticle*)fMC->GetTrack(i);
	if (!particle) continue;
	if (!fMC->IsPhysicalPrimary(i)) continue; 
	if (particle->Charge() == 0) continue;
	if (TMath::Abs(particle->Eta()) > fEtaCut) continue;
	if (particle->Pt() < fPtMin) continue;

	Double_t DPhi = DeltaPhi(particle->Phi(), fGenLeadPhi);
	if (TMath::Abs(DPhi) < pi/3.0) {// near side
	    hPhiGen[0]->Fill(DPhi);
	}
	else if (TMath::Abs(DPhi-pi) < pi/3.0) {// away side
	    hPhiGen[1]->Fill(DPhi);
	} 
	else {// transverse side
	    multTSgen++;
	    hPhiGen[2]->Fill(DPhi);

	    if (DPhi>pi/3.0 && DPhi<2*pi/3.0) {
		multTS1gen++;
		hPhiGen_TS1->Fill(DPhi);
	    } else {
		multTS2gen++;
		hPhiGen_TS2->Fill(DPhi);
	    }
	}
    }
    hNchTSGen->Fill(multTSgen);

    if (multTS2gen>=multTS1gen) {
	multTSmin_gen = multTS1gen;
	multTSmax_gen = multTS2gen;
    } else {
	multTSmin_gen = multTS2gen;
	multTSmax_gen = multTS1gen;
    }
    hNchTSminGen->Fill(multTSmin_gen);
    hNchTSmaxGen->Fill(multTSmax_gen);


    // loop over reconstructed tracks
    for (Int_t i = 0; i < fESD->GetNumberOfTracks(); i++) {
	if (i==fRecLeadIn) continue;

	AliESDtrack* track = static_cast<AliESDtrack*>(fESD->GetTrack(i));
	if (!track) continue;
	if (!fTrackFilter->IsSelected(track)) continue;
	if (TMath::Abs(track->Eta()) > fEtaCut) continue;
	if (track->Pt() < fPtMin) continue;

	Double_t DPhi = DeltaPhi(track->Phi(), fRecLeadPhi);
	if (TMath::Abs(DPhi) < pi/3.0) {// near side
	    hPhiRec[0]->Fill(DPhi);
	}
	else if (TMath::Abs(DPhi-pi) < pi/3.0) {// away side
	    hPhiRec[1]->Fill(DPhi);
	}
	else {// transverse side
	    multTSrec++;
	    hPhiRec[2]->Fill(DPhi);

	    if (DPhi>pi/3.0 && DPhi<2*pi/3.0) {
		multTS1rec++;
		hPhiRec_TS1->Fill(DPhi);
	    } else {
		multTS2rec++;
		hPhiRec_TS2->Fill(DPhi);
	    }
	}

    }
    hNchTSRec->Fill(multTSrec);
    
    if (multTS2rec>=multTS1rec) {
	multTSmin_rec = multTS1rec;
	multTSmax_rec = multTS2rec;
    } else {
	multTSmin_rec = multTS2rec;
	multTSmax_rec = multTS1rec;
    }
    hNchTSminRec->Fill(multTSmin_rec);
    hNchTSmaxRec->Fill(multTSmax_rec);

    hNchTSResponse->Fill(multTSrec,multTSgen);
    hNchTSminResponse->Fill(multTSmin_rec,multTSmin_gen);
    hNchTSmaxResponse->Fill(multTSmax_rec,multTSmax_gen);
}

//______________________________________________________________
void AliAnalysisTaskKnoUeChecks::GetMultiplicityDistributionsData()
{
    Int_t multTSrec = 0, multTS1rec = 0, multTS2rec = 0, multTSmin_rec = 0, multTSmax_rec = 0;

    // loop over reconstructed tracks
    for (Int_t i = 0; i < fESD->GetNumberOfTracks(); i++) {
	if (i==fRecLeadIn) continue;

	AliESDtrack* track = static_cast<AliESDtrack*>(fESD->GetTrack(i));
	if (!track) continue;
	if (!fTrackFilter->IsSelected(track)) continue;
	if (TMath::Abs(track->Eta()) > fEtaCut) continue;
	if (track->Pt() < fPtMin) continue;

	Double_t DPhi = DeltaPhi(track->Phi(), fRecLeadPhi);
	if (TMath::Abs(DPhi)<pi/3.0) {// near side
	    continue;
	}
	else if (TMath::Abs(DPhi-pi)<pi/3.0) {// away side
	    continue;
	}
	else {// transverse side
	    multTSrec++;

	    if (DPhi>pi/3.0 && DPhi<2*pi/3.0) {
		multTS1rec++;
		hPhiData_TS1->Fill(DPhi);
	    } else {
		multTS2rec++;
		hPhiData_TS2->Fill(DPhi);
	    }
	}
    }
    hNchTSData->Fill(multTSrec);

    if (multTS2rec >= multTS1rec) {
	multTSmin_rec = multTS1rec;
	multTSmax_rec = multTS2rec;
    } else {
	multTSmin_rec = multTS2rec;
	multTSmax_rec = multTS1rec;
    }
    hNchTSminData->Fill(multTSmin_rec);
    hNchTSmaxData->Fill(multTSmax_rec);

}

//____________________________________________________________
void AliAnalysisTaskKnoUeChecks::GetMultiplicityDistributions(){

    Int_t multTSgen = 0, multTS1gen = 0, multTS2gen = 0, multTSmin_gen = 0, multTSmax_gen = 0;
    Int_t multTSrec = 0, multTS1rec = 0, multTS2rec = 0, multTSmin_rec = 0, multTSmax_rec = 0;

    // loop over generated particles
    for (Int_t i = 0; i < fMC->GetNumberOfTracks(); i++) {
	if (i==fGenLeadIn) continue;

	AliMCParticle* particle = (AliMCParticle*)fMC->GetTrack(i);
	if (!particle) continue;
	if (!fMC->IsPhysicalPrimary(i)) continue; 
	if (particle->Charge() == 0) continue;
	if (TMath::Abs(particle->Eta()) > fEtaCut) continue;
	if (particle->Pt() < fPtMin) continue;

	Double_t DPhi = DeltaPhi(particle->Phi(), fGenLeadPhi);
	if (TMath::Abs(DPhi) < pi/3.0){// near side
	    continue;
	}
	else if (TMath::Abs(DPhi-pi) < pi/3.0){// away side
	    continue;
	}
	else {// transverse side
	    multTSgen++;

	    if (DPhi>pi/3.0 && DPhi<2*pi/3.0) {
		multTS1gen++;
		hPhiGenTest_TS1->Fill(DPhi);
	    } else {
		multTS2gen++;
		hPhiGenTest_TS2->Fill(DPhi);
	    }
	}
    }
    hNchTSGenTest->Fill(multTSgen);

    if (multTS2gen >= multTS1gen) {
	multTSmin_gen = multTS1gen;
	multTSmax_gen = multTS2gen;
    } else {
	multTSmin_gen = multTS2gen;
	multTSmax_gen = multTS1gen;
    }
    hNchTSminGenTest->Fill(multTSmin_gen);
    hNchTSmaxGenTest->Fill(multTSmax_gen);


    // loop over reconstructed tracks
    for (Int_t i = 0; i < fESD->GetNumberOfTracks(); i++) {
	if (i==fRecLeadIn) continue;

	AliESDtrack* track = static_cast<AliESDtrack*>(fESD->GetTrack(i));
	if (!track) continue;
	if (!fTrackFilter->IsSelected(track)) continue;
	if (TMath::Abs(track->Eta()) > fEtaCut) continue;
	if (track->Pt() < fPtMin) continue;

	Double_t DPhi = DeltaPhi(track->Phi(), fRecLeadPhi);
	if (TMath::Abs(DPhi) < pi/3.0) {// near side
	    continue;
	}
	else if (TMath::Abs(DPhi-pi) < pi/3.0) {// away side
	    continue;
	}
	else {// transverse side
	    multTSrec++;

	    if (DPhi>pi/3.0 && DPhi<2*pi/3.0) {
		multTS1rec++;
		hPhiRecTest_TS1->Fill(DPhi);
	    } else {
		multTS2rec++;
		hPhiRecTest_TS2->Fill(DPhi);
	    }
	}
    }
    hNchTSRecTest->Fill(multTSrec);

    if (multTS2rec >= multTS1rec) {
	multTSmin_rec = multTS1rec;
	multTSmax_rec = multTS2rec;
    } else {
	multTSmin_rec = multTS2rec;
	multTSmax_rec = multTS1rec;
    }
    hNchTSminRecTest->Fill(multTSmin_rec);
    hNchTSmaxRecTest->Fill(multTSmax_rec);
}

//______________________________________________________________
Double_t AliAnalysisTaskKnoUeChecks::DeltaPhi(Double_t phia, Double_t phib, Double_t rangeMin, Double_t rangeMax)
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

//______________________________________________________________
Bool_t AliAnalysisTaskKnoUeChecks::HasRecVertex(){

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
    if (vtTrc->IsFromVertexer3D() || vtTrc->IsFromVertexerZ()) isTrackV = false;
    const AliVVertex* vtSPD = fESD->GetPrimaryVertexSPD();


    if (vtSPD->GetNContributors() > 0) fFlag |= BIT(AliEventCuts::kVertexSPD);

    if (vtTrc->GetNContributors() > 1 && isTrackV ) fFlag |= BIT(AliEventCuts::kVertexTracks);

    if (((fFlag & BIT(AliEventCuts::kVertexTracks)) ||  !fRequireTrackVertex) && (fFlag & BIT(AliEventCuts::kVertexSPD))) fFlag |= BIT(AliEventCuts::kVertex);

    const AliVVertex* &vtx = bool(fFlag & BIT(AliEventCuts::kVertexTracks)) ? vtTrc : vtSPD;
    AliVVertex   *fPrimaryVertex = const_cast<AliVVertex*>(vtx);
    if (!fPrimaryVertex) return kFALSE;

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
