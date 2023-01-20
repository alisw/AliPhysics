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
 * Authors: Luz Tiscareño (luz.elena.tiscareno.montoya@cern.ch)           *
 *          Paola Vargas (paola.vargas.torres@cern.ch)                    *
 *          Gyula Bencedi (Gyula.Bencedi@cern.ch)                         *
 *          Antonio Ortiz (antonio.ortiz@nucleares.unam.mx)               *
 *                                                                        *
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

#include "AliAnalysisTaskChargedVsRT.h"

const Char_t *NameReg_3[3] = {"NS", "AS", "TS"};

const Int_t ptNbins = 36;
Double_t ptbins1_3[ptNbins + 1] = {
    0.0, 0.1, 0.15, 0.2,  0.25, 0.3,  0.35, 0.4,  0.45, 0.5,  0.6, 0.7, 0.8,
    0.9, 1.0, 1.25, 1.5,  2.0,  2.5,  3.0,  3.5,  4.0,  4.5,  5.0, 6.0, 7.0,
    8.0, 9.0, 10.0, 12.0, 14.0, 16.0, 18.0, 20.0, 30.0, 40.0, 50.0};

const int nBinsDCAxy = 121;
double binsDCAxy3[nBinsDCAxy + 1] = {
    -3.025, -2.975, -2.925, -2.875, -2.825, -2.775, -2.725, -2.675, -2.625,
    -2.575, -2.525, -2.475, -2.425, -2.375, -2.325, -2.275, -2.225, -2.175,
    -2.125, -2.075, -2.025, -1.975, -1.925, -1.875, -1.825, -1.775, -1.725,
    -1.675, -1.625, -1.575, -1.525, -1.475, -1.425, -1.375, -1.325, -1.275,
    -1.225, -1.175, -1.125, -1.075, -1.025, -0.975, -0.925, -0.875, -0.825,
    -0.775, -0.725, -0.675, -0.625, -0.575, -0.525, -0.475, -0.425, -0.375,
    -0.325, -0.275, -0.225, -0.175, -0.125, -0.075, -0.025, 0.025,  0.075,
    0.125,  0.175,  0.225,  0.275,  0.325,  0.375,  0.425,  0.475,  0.525,
    0.575,  0.625,  0.675,  0.725,  0.775,  0.825,  0.875,  0.925,  0.975,
    1.025,  1.075,  1.125,  1.175,  1.225,  1.275,  1.325,  1.375,  1.425,
    1.475,  1.525,  1.575,  1.625,  1.675,  1.725,  1.775,  1.825,  1.875,
    1.925,  1.975,  2.025,  2.075,  2.125,  2.175,  2.225,  2.275,  2.325,
    2.375,  2.425,  2.475,  2.525,  2.575,  2.625,  2.675,  2.725,  2.775,
    2.825,  2.875,  2.925,  2.975,  3.025};

const Double_t pi = 3.141592653589793238;
// Bin Nch
const Double_t minbinNch = -0.5;

using namespace std; // std namespace: so you can do things like 'cout' etc

ClassImp(AliAnalysisTaskChargedVsRT) // classimp: necessary for root

    AliAnalysisTaskChargedVsRT::AliAnalysisTaskChargedVsRT()
    : AliAnalysisTaskSE(), fESD(0), fEventCuts(0x0), fMCStack(0), fMC(0),
      fUseMC(kFALSE), fIsMCclosure(kFALSE), fIsHybAna(kFALSE),
      fMultPercenV0(kFALSE), fnRecHy(-1), fnRecHyWoDCA(-1), fnGen(-1),
      fNcrVar1(kFALSE), fNcrVar2(kFALSE), fTPCclustersVar1(kFALSE),
      fTPCclustersVar2(kFALSE), fGeoTPCVar1(kFALSE), fGeoTPCVar2(kFALSE),
      fGeoTPCVar3(kFALSE), fGeoTPCVar4(kFALSE), fChisqTPCVar1(kFALSE),
      fChisqTPCVar2(kFALSE), fChisqITSVar1(kFALSE), fChisqITSVar2(kFALSE),
      fDcazVar1(kFALSE), fDcazVar2(kFALSE), fPIDResponse(0x0),
      fTrackFilter(0x0), fTrackFilterwoDCA(0x0), fTrackFilterHybrid0(0x0),
      fTrackFilterHybrid1(0x0), fTrackFilterHybrid0woDCA(0x0),
      fTrackFilterHybrid1woDCA(0x0), fOutputList(0), fEtaCut(0.8), fPtMin(0.5),
      fLeadPtCutMin(5.0), fLeadPtCutMax(40.0), fNchNbin(200), fNchBinMax(199.5),
      fGenLeadPhi(0), fGenLeadPt(0), fGenLeadIn(0), fRecLeadPhi(0),
      fRecLeadPt(0), fRecLeadIn(0), ftrackmult08(0), fv0mpercentile(0),
      fdcaxy(-999), fdcaz(-999), fMultSelection(0x0), fMultSelectionbefvtx(0x0),
      hNchTSGen(0), hNchTSGenTest(0), hNchGen(0), hNchGenTest(0), hNchTSRec(0),
      hNchTSRecTest(0), hNchData(0), hNchTSData(0),
      // Only with Hybrid Track Cuts
      hPhiTotal(0),    // Sum of all the contributions
      hPhiStandard(0), // Distribution of phi without corrections -w/ SPD & ITS-
      hPhiHybrid1(0),  // Correction of phi distribution -w/o SPD & w/ ITS-
      hNchResponse(0), hNchRec(0), hNchRecTest(0), hPtInPrim(0),
      hPtInPrim_lambda(0), hPtInPrim_pion(0), hPtInPrim_kaon(0),
      hPtInPrim_proton(0), hPtInPrim_sigmap(0), hPtInPrim_sigmam(0),
      hPtInPrim_omega(0), hPtInPrim_xi(0), hPtInPrim_rest(0), hPtOut(0),
      hPtOutPrim(0), hPtOutPrim_lambda(0), hPtOutPrim_pion(0),
      hPtOutPrim_kaon(0), hPtOutPrim_proton(0), hPtOutPrim_sigmap(0),
      hPtOutPrim_sigmam(0), hPtOutPrim_omega(0), hPtOutPrim_xi(0),
      hPtOutPrim_rest(0), hPtOutSec(0), hCounter(0), hPTVsDCAData(0),
      hptvsdcaPrim(0), hptvsdcaDecs(0), hptvsdcaMatl(0), hptvsdcaAll(0),
      hV0MvsNchT(0) {

  for (Int_t i = 0; i < 3; ++i) {
    hPtVsUEGenTest[i] = 0;
    hPtVsUERecTest[i] = 0;
    hPtVsUEData[i] = 0;
    hPtVsNchGenTest[i] = 0;
    hPtVsNchRecTest[i] = 0;
    hPtVsNchData[i] = 0;
    hPhiGen[i] = 0;
    hPhiRec[i] = 0;
  }
}
//_____________________________________________________________________________
AliAnalysisTaskChargedVsRT::AliAnalysisTaskChargedVsRT(const char *name)
    : AliAnalysisTaskSE(name), fESD(0), fEventCuts(0x0), fMCStack(0), fMC(0),
      fUseMC(kFALSE), fIsMCclosure(kFALSE), fIsHybAna(kFALSE),
      fMultPercenV0(kFALSE), fnRecHy(-1), fnRecHyWoDCA(-1), fnGen(-1),
      fNcrVar1(kFALSE), fNcrVar2(kFALSE), fTPCclustersVar1(kFALSE),
      fTPCclustersVar2(kFALSE), fGeoTPCVar1(kFALSE), fGeoTPCVar2(kFALSE),
      fGeoTPCVar3(kFALSE), fGeoTPCVar4(kFALSE), fChisqTPCVar1(kFALSE),
      fChisqTPCVar2(kFALSE), fChisqITSVar1(kFALSE), fChisqITSVar2(kFALSE),
      fDcazVar1(kFALSE), fDcazVar2(kFALSE), fPIDResponse(0x0),
      fTrackFilter(0x0), fTrackFilterwoDCA(0x0), fTrackFilterHybrid0(0x0),
      fTrackFilterHybrid1(0x0), fTrackFilterHybrid0woDCA(0x0),
      fTrackFilterHybrid1woDCA(0x0), fOutputList(0), fEtaCut(0.8), fPtMin(0.5),
      fLeadPtCutMin(5.0), fLeadPtCutMax(40.0), fNchNbin(200), fNchBinMax(199.5),
      fGenLeadPhi(0), fGenLeadPt(0), fGenLeadIn(0), fRecLeadPhi(0),
      fRecLeadPt(0), fRecLeadIn(0), ftrackmult08(0), fv0mpercentile(0),
      fdcaxy(-999), fdcaz(-999), fMultSelection(0x0), fMultSelectionbefvtx(0x0),
      hNchTSGen(0), hNchTSGenTest(0), hNchGen(0), hNchGenTest(0), hNchTSRec(0),
      hNchTSRecTest(0), hNchData(0), hNchTSData(0),
      // Only with Hybrid Track Cuts
      hPhiTotal(0),    // Sum of all the contributions
      hPhiStandard(0), // Distribution of phi without corrections -w/ SPD & ITS-
      hPhiHybrid1(0),  // Correction of phi distribution -w/o SPD & w/ ITS-
      hNchResponse(0), hNchRec(0), hNchRecTest(0), hPtInPrim(0),
      hPtInPrim_lambda(0), hPtInPrim_pion(0), hPtInPrim_kaon(0),
      hPtInPrim_proton(0), hPtInPrim_sigmap(0), hPtInPrim_sigmam(0),
      hPtInPrim_omega(0), hPtInPrim_xi(0), hPtInPrim_rest(0), hPtOut(0),
      hPtOutPrim(0), hPtOutPrim_lambda(0), hPtOutPrim_pion(0),
      hPtOutPrim_kaon(0), hPtOutPrim_proton(0), hPtOutPrim_sigmap(0),
      hPtOutPrim_sigmam(0), hPtOutPrim_omega(0), hPtOutPrim_xi(0),
      hPtOutPrim_rest(0), hPtOutSec(0), hCounter(0), hPTVsDCAData(0),
      hptvsdcaPrim(0), hptvsdcaDecs(0), hptvsdcaMatl(0), hptvsdcaAll(0),
      hV0MvsNchT(0)

{

  for (Int_t i = 0; i < 3; ++i) {
    hPtVsUEGenTest[i] = 0;
    hPtVsUERecTest[i] = 0;
    hPtVsUEData[i] = 0;
    hPtVsNchGenTest[i] = 0;
    hPtVsNchRecTest[i] = 0;
    hPtVsNchData[i] = 0;
    hPhiGen[i] = 0;
    hPhiRec[i] = 0;
  }

  DefineInput(0, TChain::Class()); // define the input of the analysis: in this
                                   // case you take a 'chain' of events
  DefineOutput(1, TList::Class()); // define the ouptut of the analysis: in this
                                   // case it's a list of histograms
}
//_____________________________________________________________________________
AliAnalysisTaskChargedVsRT::~AliAnalysisTaskChargedVsRT() {
  // destructor
  if (fOutputList) {
    delete fOutputList; // at the end of your task, it is deleted from memory by
                        // calling this function
    fOutputList = 0x0;
  }
}
//_____________________________________________________________________________
void AliAnalysisTaskChargedVsRT::UserCreateOutputObjects() {
  // AliAnalysisManager *man = AliAnalysisManager::GetAnalysisManager();
  // if(man){
  //     AliInputEventHandler* inputHandler =
  //     (AliInputEventHandler*)(man->GetInputEventHandler());
  //     if(inputHandler)fPIDResponse = inputHandler->GetPIDResponse();
  // }

  // Define Track Cuts
  fTrackFilter = new AliAnalysisFilter("trackFilter");
  fTrackFilterwoDCA = new AliAnalysisFilter("trackFilterwoDCA"); // wo DCA cut
  AliESDtrackCuts *fCuts = new AliESDtrackCuts();

  AliESDtrackCuts *fCutswoDCA = new AliESDtrackCuts();

  SetCutsFilterWoDCA(fCutswoDCA);

  fCuts->SetAcceptKinkDaughters(kFALSE); //
  fCuts->SetRequireTPCRefit(kTRUE);      //
  fCuts->SetRequireITSRefit(kTRUE);      //
  fCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,
                                  AliESDtrackCuts::kAny); //
  fCuts->SetDCAToVertex2D(kFALSE);                        //
  fCuts->SetRequireSigmaToVertex(kFALSE);                 //
  fCuts->SetEtaRange(-0.8, 0.8);

  // to be considered for systematics...
  fCuts->SetMinNCrossedRowsTPC(70); // TBC
  fCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);
  fCuts->SetMaxChi2PerClusterTPC(4);
  fCuts->SetMaxDCAToVertexZ(2);

  // if (fGeoTPCVar1) {fCuts->SetCutGeoNcrNcl(2., 130., 1.5, 0.85, 0.7);}//
  // else if (fGeoTPCVar2) {fCuts->SetCutGeoNcrNcl(4., 130., 1.5, 0.85, 0.7);}//
  // else if (fGeoTPCVar3) {fCuts->SetCutGeoNcrNcl(3., 120., 1.5, 0.85, 0.7);}//
  // else if (fGeoTPCVar4) {fCuts->SetCutGeoNcrNcl(3., 140., 1.5, 0.85, 0.7);}//
  // else {fCuts->SetCutGeoNcrNcl(3., 130., 1.5, 0.85, 0.7);}// Default
  fCuts->SetCutGeoNcrNcl(3., 130., 1.5, 0.85, 0.7); // Default
  // if (fNcrVar1)
  // {fCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(0.7);}// else if
  // (fNcrVar2) {fCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(0.9);}//
  // else {fCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);}// Default
  fCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);
  // if (fChisqTPCVar1) {fCuts->SetMaxChi2PerClusterTPC(3);}
  // else if (fChisqTPCVar2) {fCuts->SetMaxChi2PerClusterTPC(5);}
  // else {fCuts->SetMaxChi2PerClusterTPC(4);}// Default
  fCuts->SetMaxChi2PerClusterTPC(4);
  // if (fDcazVar1) {fCuts->SetMaxDCAToVertexZ(1);} // DCAz = 1 cm
  // else if (fDcazVar2) {fCuts->SetMaxDCAToVertexZ(5);} // DCAz = 5 cm
  // else {fCuts->SetMaxDCAToVertexZ(2);}// Default
  fCuts->SetMaxDCAToVertexZ(2);
  // if (fChisqITSVar1) {fCuts->SetMaxChi2PerClusterITS(25);}//
  // else if (fChisqITSVar2) {fCuts->SetMaxChi2PerClusterITS(49);}//
  // else {fCuts->SetMaxChi2PerClusterITS(36);}// Default
  fCuts->SetMaxChi2PerClusterITS(36);
  fCuts->SetMaxDCAToVertexXYPtDep("0.0105+0.0350/pt^1.1");
  fCuts->SetMaxChi2PerClusterITS(36);

  fTrackFilter->AddCuts(fCuts);

  // Define Hybrid 0 (global, 2011 track cuts)
  //
  fTrackFilterHybrid0 = new AliAnalysisFilter("trackFilterHybrid0");
  fTrackFilterHybrid0woDCA = new AliAnalysisFilter("trackFilterHybrid0woDCA");
  AliESDtrackCuts *fCutsHybrid0 = new AliESDtrackCuts();
  fCutsHybrid0 = new AliESDtrackCuts("fCutsHybrid0");

  AliESDtrackCuts *fCutsHybrid0woDCA = new AliESDtrackCuts();
  fCutsHybrid0woDCA = new AliESDtrackCuts("fCutsHybrid0woDCA");

  SetCutsHybrid0WoDCA(fCutsHybrid0woDCA);

  // to be considered for systematics...
  if (fNcrVar1)
    fCutsHybrid0->SetMinNCrossedRowsTPC(60);
  else if (fNcrVar2)
    fCutsHybrid0->SetMinNCrossedRowsTPC(100);
  else
    fCutsHybrid0->SetMinNCrossedRowsTPC(70); // Default

  if (fTPCclustersVar1)
    fCutsHybrid0->SetMinRatioCrossedRowsOverFindableClustersTPC(0.7);
  else if (fTPCclustersVar2)
    fCutsHybrid0->SetMinRatioCrossedRowsOverFindableClustersTPC(0.9);
  else
    fCutsHybrid0->SetMinRatioCrossedRowsOverFindableClustersTPC(0.8); // Default

  if (fChisqTPCVar1)
    fCutsHybrid0->SetMaxChi2PerClusterTPC(3);
  else if (fChisqTPCVar2)
    fCutsHybrid0->SetMaxChi2PerClusterTPC(5);
  else
    fCutsHybrid0->SetMaxChi2PerClusterTPC(4); // Default

  if (fChisqITSVar1)
    fCutsHybrid0->SetMaxChi2PerClusterITS(25);
  else if (fChisqITSVar2)
    fCutsHybrid0->SetMaxChi2PerClusterITS(49);
  else
    fCutsHybrid0->SetMaxChi2PerClusterITS(36); // Default

  if (fDcazVar1)
    fCutsHybrid0->SetMaxDCAToVertexZ(1);
  else if (fDcazVar2)
    fCutsHybrid0->SetMaxDCAToVertexZ(3);
  else
    fCutsHybrid0->SetMaxDCAToVertexZ(2); // Default

  fCutsHybrid0->SetAcceptKinkDaughters(kFALSE);
  fCutsHybrid0->SetRequireTPCRefit(kTRUE);
  fCutsHybrid0->SetRequireITSRefit(kTRUE);
  fCutsHybrid0->SetClusterRequirementITS(AliESDtrackCuts::kSPD,
                                         AliESDtrackCuts::kOff);
  fCutsHybrid0->SetDCAToVertex2D(kFALSE);
  fCutsHybrid0->SetRequireSigmaToVertex(kFALSE);
  fCutsHybrid0->SetEtaRange(-0.8, 0.8);
  fCutsHybrid0->SetMaxDCAToVertexXYPtDep("0.0105+0.0350/pt^1.1");
  fTrackFilterHybrid0->AddCuts(fCutsHybrid0);

  // Define Hybrid 1
  //
  fTrackFilterHybrid1 = new AliAnalysisFilter("trackFilterHybrid1");
  fTrackFilterHybrid1woDCA = new AliAnalysisFilter("trackFilterHybrid1woDCA");
  AliESDtrackCuts *fCutsHybrid1 = new AliESDtrackCuts();
  fCutsHybrid1 = new AliESDtrackCuts("fCutsHybrid1");

  AliESDtrackCuts *fCutsHybrid1woDCA = new AliESDtrackCuts();
  fCutsHybrid1woDCA = new AliESDtrackCuts("fCutsHybrid1woDCA");

  SetCutsHybrid1WoDCA(fCutsHybrid1woDCA);

  // to be considered for systematics...
  if (fNcrVar1)
    fCutsHybrid1->SetMinNCrossedRowsTPC(60);
  else if (fNcrVar2)
    fCutsHybrid1->SetMinNCrossedRowsTPC(100);
  else
    fCutsHybrid1->SetMinNCrossedRowsTPC(70); // Default

  if (fTPCclustersVar1)
    fCutsHybrid1->SetMinRatioCrossedRowsOverFindableClustersTPC(0.7);
  else if (fTPCclustersVar2)
    fCutsHybrid1->SetMinRatioCrossedRowsOverFindableClustersTPC(0.9);
  else
    fCutsHybrid1->SetMinRatioCrossedRowsOverFindableClustersTPC(0.8); // Default

  if (fChisqTPCVar1)
    fCutsHybrid1->SetMaxChi2PerClusterTPC(3);
  else if (fChisqTPCVar2)
    fCutsHybrid1->SetMaxChi2PerClusterTPC(5);
  else
    fCutsHybrid1->SetMaxChi2PerClusterTPC(4); // Default

  if (fChisqITSVar1)
    fCutsHybrid1->SetMaxChi2PerClusterITS(25);
  else if (fChisqITSVar2)
    fCutsHybrid1->SetMaxChi2PerClusterITS(49);
  else
    fCutsHybrid1->SetMaxChi2PerClusterITS(36); // Default

  if (fDcazVar1)
    fCutsHybrid1->SetMaxDCAToVertexZ(1);
  else if (fDcazVar2)
    fCutsHybrid1->SetMaxDCAToVertexZ(3);
  else
    fCutsHybrid1->SetMaxDCAToVertexZ(2); // Default

  fCutsHybrid1->SetAcceptKinkDaughters(kFALSE);
  fCutsHybrid1->SetRequireTPCRefit(kTRUE);
  fCutsHybrid1->SetRequireITSRefit(kFALSE);
  fCutsHybrid1->SetClusterRequirementITS(AliESDtrackCuts::kSPD,
                                         AliESDtrackCuts::kNone);
  fCutsHybrid1->SetDCAToVertex2D(kFALSE);
  fCutsHybrid1->SetRequireSigmaToVertex(kFALSE);
  fCutsHybrid1->SetEtaRange(-0.8, 0.8);
  fCutsHybrid1->SetMaxDCAToVertexXYPtDep("0.0105+0.0350/pt^1.1");

  fTrackFilterHybrid1->AddCuts(fCutsHybrid1);

  // create output objects

  OpenFile(1);
  fOutputList =
      new TList(); // this is a list which will contain all of your histograms
  // at the end of the analysis, the contents of this list are written  to the
  // output file
  fOutputList->SetOwner(kTRUE); // memory stuff: the list is owner of all
                                // objects and will delete them if requested

  if (fUseMC) {

    hNchTSGen = new TH1D("hNchTSGen", "", fNchNbin, minbinNch, fNchBinMax);
    fOutputList->Add(hNchTSGen);

    hNchTSGenTest =
        new TH1D("hNchTSGenTest", "", fNchNbin, minbinNch, fNchBinMax);
    fOutputList->Add(hNchTSGenTest);

    hNchTSRecTest =
        new TH1D("hNchTSRecTest", "", fNchNbin, minbinNch, fNchBinMax);
    fOutputList->Add(hNchTSRecTest);

    hNchGen = new TH1D("hNchGen", "", fNchNbin, minbinNch, fNchBinMax);
    fOutputList->Add(hNchGen);

    hNchGenTest = new TH1D("hNchGenTest", "", fNchNbin, minbinNch, fNchBinMax);
    fOutputList->Add(hNchGenTest);

    hNchRecTest = new TH1D("hNchRecTest", "", fNchNbin, minbinNch, fNchBinMax);
    fOutputList->Add(hNchRecTest);

    for (Int_t i = 0; i < 3; ++i) {
      hPhiGen[i] = new TH1D(Form("hPhiGen_%s", NameReg_3[i]), "", 64,
                            -TMath::Pi() / 2.0, 3.0 * TMath::Pi() / 2.0);
      fOutputList->Add(hPhiGen[i]);
    }

    hNchResponse = new TH2D(
        "hNchResponse", "Detector response; rec mult; gen mult", fNchNbin,
        minbinNch, fNchBinMax, fNchNbin, minbinNch, fNchBinMax);
    fOutputList->Add(hNchResponse);

    hPtInPrim =
        new TH1D("hPtInPrim", "pT prim true; pT; Nch", ptNbins, ptbins1_3);
    fOutputList->Add(hPtInPrim);

    hPtOut = new TH1D("hPtOut", "pT all rec; pT; Nch", ptNbins, ptbins1_3);
    fOutputList->Add(hPtOut);

    hPtOutPrim =
        new TH1D("hPtOutPrim", "pT prim rec; pT; Nch", ptNbins, ptbins1_3);
    fOutputList->Add(hPtOutPrim);

    hPtOutSec =
        new TH1D("hPtOutSec", "pT sec rec; pT; Nch", ptNbins, ptbins1_3);
    fOutputList->Add(hPtOutSec);

    if (!fIsMCclosure) {
      hPtInPrim_lambda = new TH1D("hPtInPrim_lambda", "pT prim true; pT; Nch",
                                  ptNbins, ptbins1_3);
      fOutputList->Add(hPtInPrim_lambda);

      hPtInPrim_pion = new TH1D("hPtInPrim_pion", "pT prim true; pT; Nch",
                                ptNbins, ptbins1_3);
      fOutputList->Add(hPtInPrim_pion);

      hPtInPrim_kaon = new TH1D("hPtInPrim_kaon", "pT prim true; pT; Nch",
                                ptNbins, ptbins1_3);
      fOutputList->Add(hPtInPrim_kaon);

      hPtInPrim_proton = new TH1D("hPtInPrim_proton", "pT prim true; pT; Nch",
                                  ptNbins, ptbins1_3);
      fOutputList->Add(hPtInPrim_proton);

      hPtInPrim_sigmap = new TH1D("hPtInPrim_sigmap", "pT prim true; pT; Nch",
                                  ptNbins, ptbins1_3);
      fOutputList->Add(hPtInPrim_sigmap);

      hPtInPrim_sigmam = new TH1D("hPtInPrim_sigmam", "pT prim true; pT; Nch",
                                  ptNbins, ptbins1_3);
      fOutputList->Add(hPtInPrim_sigmam);

      hPtInPrim_omega = new TH1D("hPtInPrim_omega", "pT prim true; pT; Nch",
                                 ptNbins, ptbins1_3);
      fOutputList->Add(hPtInPrim_omega);

      hPtInPrim_xi =
          new TH1D("hPtInPrim_xi", "pT prim true; pT; Nch", ptNbins, ptbins1_3);
      fOutputList->Add(hPtInPrim_xi);

      hPtInPrim_rest = new TH1D("hPtInPrim_rest", "pT prim true; pT; Nch",
                                ptNbins, ptbins1_3);
      fOutputList->Add(hPtInPrim_rest);

      hPtOutPrim_lambda = new TH1D("hPtOutPrim_lambda", "pT prim true; pT; Nch",
                                   ptNbins, ptbins1_3);
      fOutputList->Add(hPtOutPrim_lambda);

      hPtOutPrim_pion = new TH1D("hPtOutPrim_pion", "pT prim rec; pT; Nch",
                                 ptNbins, ptbins1_3);
      fOutputList->Add(hPtOutPrim_pion);

      hPtOutPrim_kaon = new TH1D("hPtOutPrim_kaon", "pT prim rec; pT; Nch",
                                 ptNbins, ptbins1_3);
      fOutputList->Add(hPtOutPrim_kaon);

      hPtOutPrim_proton = new TH1D("hPtOutPrim_proton", "pT prim rec; pT; Nch",
                                   ptNbins, ptbins1_3);
      fOutputList->Add(hPtOutPrim_proton);

      hPtOutPrim_sigmap = new TH1D("hPtOutPrim_sigmap", "pT prim rec; pT; Nch",
                                   ptNbins, ptbins1_3);
      fOutputList->Add(hPtOutPrim_sigmap);

      hPtOutPrim_sigmam = new TH1D("hPtOutPrim_sigmam", "pT prim rec; pT; Nch",
                                   ptNbins, ptbins1_3);
      fOutputList->Add(hPtOutPrim_sigmam);

      hPtOutPrim_omega = new TH1D("hPtOutPrim_omega", "pT prim rec; pT; Nch",
                                  ptNbins, ptbins1_3);
      fOutputList->Add(hPtOutPrim_omega);

      hPtOutPrim_xi =
          new TH1D("hPtOutPrim_xi", "pT prim rec; pT; Nch", ptNbins, ptbins1_3);
      fOutputList->Add(hPtOutPrim_xi);

      hPtOutPrim_rest = new TH1D("hPtOutPrim_rest", "pT prim rec; pT; Nch",
                                 ptNbins, ptbins1_3);
      fOutputList->Add(hPtOutPrim_rest);
    }

    for (Int_t i = 0; i < 3; ++i) {

      hPtVsUEGenTest[i] = new TH2D(Form("hPtVsUEGenTest_%s", NameReg_3[i]),
                                   "gen pT vs nch_transverse", fNchNbin,
                                   minbinNch, fNchBinMax, ptNbins, ptbins1_3);
      fOutputList->Add(hPtVsUEGenTest[i]);

      hPtVsNchGenTest[i] = new TH2D(Form("hPtVsNchGenTest_%s", NameReg_3[i]),
                                    "gen pT vs nch_transverse", fNchNbin,
                                    minbinNch, fNchBinMax, ptNbins, ptbins1_3);
      fOutputList->Add(hPtVsNchGenTest[i]);
    }

    hptvsdcaPrim = new TH2F("hptvsdcaPrim", "pt vs dca Primaries", ptNbins,
                            ptbins1_3, nBinsDCAxy, binsDCAxy3);
    fOutputList->Add(hptvsdcaPrim);

    hptvsdcaDecs = new TH2F("hptvsdcaDecs", "pt vs dca Decays", ptNbins,
                            ptbins1_3, nBinsDCAxy, binsDCAxy3);
    fOutputList->Add(hptvsdcaDecs);

    hptvsdcaMatl = new TH2F("hptvsdcaMatl", "pt vs dca Material", ptNbins,
                            ptbins1_3, nBinsDCAxy, binsDCAxy3);
    fOutputList->Add(hptvsdcaMatl);

    hptvsdcaAll = new TH2F("hptvsdcaAll", "pt vs dca all", ptNbins, ptbins1_3,
                           nBinsDCAxy, binsDCAxy3);
    fOutputList->Add(hptvsdcaAll);
  }

  hNchTSRec = new TH1D("hNchTSRec", "", fNchNbin, minbinNch, fNchBinMax);
  fOutputList->Add(hNchTSRec);

  hNchRec = new TH1D("hNchRec", "", fNchNbin, minbinNch, fNchBinMax);
  fOutputList->Add(hNchTSRec);
  hPhiTotal = new TH1F("hPhiSum", "#varphi", 50, -TMath::Pi() / 2.0,
                       5.0 * TMath::Pi() / 2.0);
  fOutputList->Add(hPhiTotal);

  hPhiStandard = new TH1F("hPhiSPD&ITS", "#varphi", 50, -TMath::Pi() / 2.0,
                          5.0 * TMath::Pi() / 2.0);
  fOutputList->Add(hPhiStandard);

  hPhiHybrid1 = new TH1F("hPhiITS", "#varphi", 50, -TMath::Pi() / 2.0,
                         5.0 * TMath::Pi() / 2.0);
  fOutputList->Add(hPhiHybrid1);

  hV0MvsNchT = new TProfile("hV0MvsNchT", "", fNchNbin, minbinNch, fNchBinMax);
  fOutputList->Add(hV0MvsNchT);

  for (Int_t i = 0; i < 3; ++i) {
    hPhiRec[i] = new TH1D(Form("hPhiRec_%s", NameReg_3[i]), "", 64,
                          -TMath::Pi() / 2.0, 3.0 * TMath::Pi() / 2.0);
    fOutputList->Add(hPhiRec[i]);
  }

  for (Int_t i = 0; i < 3; ++i) {

    hPtVsUERecTest[i] = new TH2D(Form("hPtVsUERecTest_%s", NameReg_3[i]),
                                 "rec pT vs nch_transverse", fNchNbin,
                                 minbinNch, fNchBinMax, ptNbins, ptbins1_3);
    fOutputList->Add(hPtVsUERecTest[i]);

    hPtVsNchRecTest[i] = new TH2D(Form("hPtVsNchRecTest_%s", NameReg_3[i]),
                                  "rec pT vs nch_transverse", fNchNbin,
                                  minbinNch, fNchBinMax, ptNbins, ptbins1_3);
    fOutputList->Add(hPtVsNchRecTest[i]);
  }

  if (!fUseMC) {
    hNchTSData = new TH1D("hNchTSData", "", fNchNbin, minbinNch, fNchBinMax);
    fOutputList->Add(hNchTSData);

    hNchData = new TH1D("hNchData", "", fNchNbin, minbinNch, fNchBinMax);
    fOutputList->Add(hNchData);

    for (Int_t i = 0; i < 3; ++i) {

      hPtVsUEData[i] = new TH2D(Form("hPtVsUEData_%s", NameReg_3[i]),
                                "data pT vs nch_transverse", fNchNbin,
                                minbinNch, fNchBinMax, ptNbins, ptbins1_3);
      fOutputList->Add(hPtVsUEData[i]);

      hPtVsNchData[i] = new TH2D(Form("hPtVsNchData_%s", NameReg_3[i]),
                                 "data pT vs nch_transverse", fNchNbin,
                                 minbinNch, fNchBinMax, ptNbins, ptbins1_3);
      fOutputList->Add(hPtVsNchData[i]);
    }

    hPTVsDCAData = new TH2D("hPTVsDCAData", "pT vs dcaxy", ptNbins, ptbins1_3,
                            nBinsDCAxy, binsDCAxy3);
    fOutputList->Add(hPTVsDCAData);
  }

  fEventCuts.AddQAplotsToList(fOutputList);
  PostData(1, fOutputList); // postdata will notify the analysis manager of
                            // changes / updates to the
}
//_____________________________________________________________________________
void AliAnalysisTaskChargedVsRT::UserExec(Option_t *) {

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

  vector<Float_t> ptMc;
  vector<Float_t> phiMc;
  vector<Int_t> idMc;
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

    fnGen = FillArrayMC(ptMc, phiMc, idMc);

    // Before trigger selection
    GetLeadingObjectFromArray(ptMc, phiMc, idMc, fnGen, kTRUE);
    // Filling histos for observables at generator level
    if (fGenLeadPt >= fLeadPtCutMin && fGenLeadPt < fLeadPtCutMax)
      GetMultiplicityDistributionsTrue(phiMc, ptMc, fnGen, idMc);
  }

  // Define the arrays for hybrid, standard (w DCA cut)
  vector<Float_t> ptHy;
  vector<Float_t> phiHy;
  vector<Float_t> dcaxyHy;
  vector<Float_t> dcazHy;
  vector<Int_t> isprimHy;
  vector<Int_t> idHy;
  fnRecHy =
      FillArray(ptHy, phiHy, dcaxyHy, dcazHy, isprimHy, idHy, kTRUE, fIsHybAna);

  // Before trigger selection
  GetLeadingObjectFromArray(ptHy, phiHy, idHy, fnRecHy, kFALSE);
  // Now we get the leading particle pT

  // Define the arrays for hybrid, standard (wo DCA cut)
  vector<Float_t> ptHyWoDCA;
  vector<Float_t> phiHyWoDCA;
  vector<Float_t> dcaxyHyWoDCA;
  vector<Float_t> dcazHyWoDCA;
  vector<Int_t> isprimHyWoDCA;
  vector<Int_t> idHyWoDCA;
  fnRecHyWoDCA = FillArray(ptHyWoDCA, phiHyWoDCA, dcaxyHyWoDCA, dcazHyWoDCA,
                           isprimHyWoDCA, idHyWoDCA, kFALSE, fIsHybAna);

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

  fMultSelectionbefvtx =
      (AliMultSelection *)fESD->FindListObject("MultSelection");

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
  if (fMultPercenV0) {
    fv0mpercentile = fMultSelection->GetMultiplicityPercentile("V0M");
  } else {
    fv0mpercentile = fMultSelection->GetMultiplicityPercentile("V0A");
  }

  if (fIsMCclosure) {
    Double_t randomUE = -1;
    gRandom->SetSeed(0);
    randomUE = gRandom->Uniform(0.0, 1.0);
    if (randomUE < 0.5) { // corrections (50% stat.)
      if (isGoodVtxPosMC) {
        // KNO scaling
        if ((fGenLeadPt >= fLeadPtCutMin && fGenLeadPt < fLeadPtCutMax) &&
            (fRecLeadPt >= fLeadPtCutMin && fRecLeadPt < fLeadPtCutMax)) {
          GetDetectorResponse(phiMc, fnGen, phiHy, fnRecHy, idMc);
          GetBinByBinCorrections(fnGen, fnRecHy, ptMc, ptHy, idMc, idHy,
                                 isprimHy);
        }
      }
    } else { // for testing the method
      // KNO scaling
      if ((fRecLeadPt >= fLeadPtCutMin && fRecLeadPt < fLeadPtCutMax)) {
        GetMultiplicityDistributions(phiHy, ptHy, fnRecHy, ptHyWoDCA,
                                     dcaxyHyWoDCA, isprimHyWoDCA, fnRecHyWoDCA);
      }
    }
  } else {
    if (fUseMC) {
      if (isGoodVtxPosMC) {
        // KNO scaling
        if ((fGenLeadPt >= fLeadPtCutMin && fGenLeadPt < fLeadPtCutMax) &&
            (fRecLeadPt >= fLeadPtCutMin && fRecLeadPt < fLeadPtCutMax)) {
          GetDetectorResponse(phiMc, fnGen, phiHy, fnRecHy, idMc);
          GetBinByBinCorrections(fnGen, fnRecHy, ptMc, ptHy, idMc, idHy,
                                 isprimHy);
          GetMultiplicityDistributions(phiHy, ptHy, fnRecHy, ptHyWoDCA,
                                       dcaxyHyWoDCA, isprimHyWoDCA,
                                       fnRecHyWoDCA);
        }
      }
    } else {
      // KNO scaling
      if ((fRecLeadPt >= fLeadPtCutMin && fRecLeadPt < fLeadPtCutMax))

        GetMultiplicityDistributionsData(phiHy, ptHy, fnRecHy, ptHyWoDCA,
                                         dcaxyHyWoDCA, fnRecHyWoDCA);
    }
  }

  PostData(1, fOutputList); // stream the result of this event to the output
                            // manager which will write it to a file
}

//______________________________________________________________________________
void AliAnalysisTaskChargedVsRT::Terminate(Option_t *) {}

//_____________________________________________________________________________
void AliAnalysisTaskChargedVsRT::GetLeadingObjectFromArray(
    const vector<Float_t> &pt, const vector<Float_t> &phi,
    const vector<Int_t> &id, Int_t multPart, Bool_t isMC) {

  Float_t flPt = 0; // leading pT
  Float_t flPhi = 0;
  Int_t flIndex = 0;

  for (Int_t i1 = 0; i1 < multPart; ++i1) {

    if (isMC) {
      if (id[i1] <= 0 || id[i1] > 8) { // only charged particles
        continue;
      }
    }

    if (flPt < pt[i1]) {
      flPt = pt[i1];
      flPhi = phi[i1];
      flIndex = i1;
    }
  }
  if (isMC) {
    fGenLeadPhi = flPhi;
    fGenLeadPt = flPt;
    fGenLeadIn = flIndex;
  } else {
    fRecLeadPhi = flPhi;
    fRecLeadPt = flPt;
    fRecLeadIn = flIndex;
  }
}

//_____________________________________________________________________________
void AliAnalysisTaskChargedVsRT::GetBinByBinCorrections(
    Int_t multGen, Int_t multRec, const vector<Float_t> &ptGen,
    const vector<Float_t> &ptRec, const vector<Int_t> &idGen,
    const vector<Int_t> &idRec, const vector<Int_t> &isprimRec) {

  // Histos for efficiencyxacceptance
  for (Int_t i = 0; i < multGen; ++i) {
    if (idGen[i] >= 0 && idGen[i] <= 8) {
      if (idGen[i] > 0) {
        hPtInPrim->Fill(
            ptGen[i]); // inital pT distribution (MC gen only charged particles)
      }
      if (!fIsMCclosure) {
        if (idGen[i] == 0)
          hPtInPrim_lambda->Fill(ptGen[i]); // lambdas
        else if (idGen[i] == 1)
          hPtInPrim_pion->Fill(ptGen[i]); // pions
        else if (idGen[i] == 2)
          hPtInPrim_kaon->Fill(ptGen[i]); // kaons
        else if (idGen[i] == 3)
          hPtInPrim_proton->Fill(ptGen[i]); // protons
        else if (idGen[i] == 4)
          hPtInPrim_sigmap->Fill(ptGen[i]); // sigma plus
        else if (idGen[i] == 5)
          hPtInPrim_sigmam->Fill(ptGen[i]); // sigma minus
        else
          hPtInPrim_rest->Fill(ptGen[i]); // rest of the charged particles
        if (idGen[i] == 6)
          hPtInPrim_omega->Fill(ptGen[i]); // Omega
        if (idGen[i] == 7)
          hPtInPrim_xi->Fill(ptGen[i]); // Xi
      }
    }
  }

  for (Int_t i = 0; i < multRec; ++i) { // loop over all these tracks
    hPtOut->Fill(ptRec[i]);
    if (idRec[i] >= 0 && idRec[i] <= 8) {
      if (isprimRec[i] == 0) {
        hPtOutPrim->Fill(ptRec[i]);
        if (!fIsMCclosure) {
          if (idRec[i] == 0)
            hPtOutPrim_lambda->Fill(ptRec[i]); // lambdas
          else if (idRec[i] == 1)
            hPtOutPrim_pion->Fill(ptRec[i]); // pions
          else if (idRec[i] == 2)
            hPtOutPrim_kaon->Fill(ptRec[i]); // kaons
          else if (idRec[i] == 3)
            hPtOutPrim_proton->Fill(ptRec[i]); // protons
          else if (idRec[i] == 4)
            hPtOutPrim_sigmap->Fill(ptRec[i]); // sigma plus
          else if (idRec[i] == 5)
            hPtOutPrim_sigmam->Fill(ptRec[i]); // sigma minus
          else
            hPtOutPrim_rest->Fill(ptRec[i]); // rest of the charged particles
          if (idRec[i] == 6)
            hPtOutPrim_omega->Fill(ptRec[i]); // Omega
          if (idRec[i] == 7)
            hPtOutPrim_xi->Fill(ptRec[i]); // Xi
        }
      }
    }
    if (isprimRec[i] == 1 || isprimRec[i] == 2) {
      hPtOutSec->Fill(ptRec[i]);
    }
  }
}

void AliAnalysisTaskChargedVsRT::GetDetectorResponse(
    const vector<Float_t> &phiGen, Int_t multGen, const vector<Float_t> &phiRec,
    Int_t multRec, const vector<Int_t> &idGen) {

  Int_t multTSgen = 0;
  Int_t multTSrec = 0;
  for (Int_t i = 0; i < multGen; ++i) {

    if (i == fGenLeadIn)
      continue;
    Double_t DPhi = DeltaPhi(phiGen[i], fGenLeadPhi);

    // definition of the topological regions
    if (TMath::Abs(DPhi) < pi / 3.0) { // near side
      hPhiGen[0]->Fill(DPhi);
    } else if (TMath::Abs(DPhi - pi) < pi / 3.0) { // away side
      hPhiGen[1]->Fill(DPhi);
    } else { // transverse side
      if (idGen[i] > 0 && idGen[i] <= 8) {
        multTSgen++;
      }
      hPhiGen[2]->Fill(DPhi);
    }
  }

  hNchTSGen->Fill(multTSgen);

  for (Int_t i = 0; i < multRec; ++i) { // loop over all these tracks

    if (i == fRecLeadIn)
      continue;
    Double_t DPhi = DeltaPhi(phiRec[i], fRecLeadPhi);

    // definition of the topological regions
    if (TMath::Abs(DPhi) < pi / 3.0) { // near side
      hPhiRec[0]->Fill(DPhi);
    } else if (TMath::Abs(DPhi - pi) < pi / 3.0) { // away side
      hPhiRec[1]->Fill(DPhi);
    } else { // transverse side
      multTSrec++;
      hPhiRec[2]->Fill(DPhi);
    }
  }

  hNchTSRec->Fill(multTSrec);
  hNchResponse->Fill(multTSrec, multTSgen);
}
//______________________________________________________________
void AliAnalysisTaskChargedVsRT::GetMultiplicityDistributionsData(
    const vector<Float_t> &phiRec, const vector<Float_t> &ptRec, Int_t multRec,
    const vector<Float_t> &ptRecWoDCA, const vector<Float_t> &dcaxyRecWoDCA,
    Int_t multRecWoDCA) {

  Int_t multTSrec = 0;

  for (Int_t i = 0; i < multRec; ++i) { // loop over all these tracks

    if (i == fRecLeadIn)
      continue;
    Double_t DPhi = DeltaPhi(phiRec[i], fRecLeadPhi);

    // definition of the topological regions
    if (TMath::Abs(DPhi) < pi / 3.0) { // near side
      continue;
    } else if (TMath::Abs(DPhi - pi) < pi / 3.0) { // away side
      continue;
    } else { // transverse side
      multTSrec++;
    }
  }

  hNchTSData->Fill(multTSrec);
  hNchData->Fill(multRec);
  hV0MvsNchT->Fill(multRec, fv0mpercentile);
  // Filling rec pT vs UE (for pT *** considering the hybrid track cuts or the
  // 2011 track cuts ***)
  for (Int_t i = 0; i < multRec; ++i) { // loop over all these tracks

    if (i == fRecLeadIn)
      continue;
    Double_t DPhi = DeltaPhi(phiRec[i], fRecLeadPhi);

    // definition of the topological regions
    if (TMath::Abs(DPhi) < pi / 3.0) { // near side
      hPtVsUEData[0]->Fill(multTSrec, ptRec[i]);
      hPtVsNchData[0]->Fill(multRec, ptRec[i]);
    } else if (TMath::Abs(DPhi - pi) < pi / 3.0) { // away side
      hPtVsUEData[1]->Fill(multTSrec, ptRec[i]);
      hPtVsNchData[1]->Fill(multRec, ptRec[i]);
    } else { // transverse side
      hPtVsUEData[2]->Fill(multTSrec, ptRec[i]);
      hPtVsNchData[2]->Fill(multRec, ptRec[i]);
    }
  }

  // Filling the histos wo DCA cuts
  for (Int_t i = 0; i < multRecWoDCA; ++i) {

    hPTVsDCAData->Fill(ptRecWoDCA[i], dcaxyRecWoDCA[i]);
  }
}
//_____________________________________________________________
void AliAnalysisTaskChargedVsRT::GetMultiplicityDistributionsTrue(
    const vector<Float_t> &phiGen, const vector<Float_t> &ptGen, Int_t multGen,
    const vector<Int_t> &idGen) {

  Int_t multTSgen = 0;

  for (Int_t i = 0; i < multGen; ++i) {

    if (i == fGenLeadIn)
      continue;

    Double_t DPhi = DeltaPhi(phiGen[i], fGenLeadPhi);

    // definition of the topological regions
    if (TMath::Abs(DPhi) < pi / 3.0) { // near side
      continue;
    } else if (TMath::Abs(DPhi - pi) < pi / 3.0) { // away side
      continue;
    } else { // transverse side
      if (idGen[i] > 0 && idGen[i] <= 8) {
        multTSgen++;
      }
    }
  }

  hNchGen->Fill(multGen);
  hNchTSGenTest->Fill(multTSgen);
  hNchGenTest->Fill(multGen);

  for (Int_t i = 0; i < multGen; ++i) {

    if (idGen[i] <= 0 || idGen[i] > 8) {
      continue;
    }

    if (i == fGenLeadIn)
      continue;
    Double_t DPhi = DeltaPhi(phiGen[i], fGenLeadPhi);

    // definition of the topological regions
    if (TMath::Abs(DPhi) < pi / 3.0) { // near side
      hPtVsUEGenTest[0]->Fill(multTSgen, ptGen[i]);
      hPtVsNchGenTest[0]->Fill(multGen, ptGen[i]);
    } else if (TMath::Abs(DPhi - pi) < pi / 3.0) { // away side
      hPtVsUEGenTest[1]->Fill(multTSgen, ptGen[i]);
      hPtVsNchGenTest[1]->Fill(multGen, ptGen[i]);
    } else { // transverse side
      hPtVsUEGenTest[2]->Fill(multTSgen, ptGen[i]);
      hPtVsNchGenTest[2]->Fill(multGen, ptGen[i]);
    }
  }
}

//____________________________________________________________
void AliAnalysisTaskChargedVsRT::GetMultiplicityDistributions(
    const vector<Float_t> &phiRec, const vector<Float_t> &ptRec, Int_t multRec,
    const vector<Float_t> &ptRecWoDCA, const vector<Float_t> &dcaxyRecWoDCA,
    const vector<Int_t> &isprimRecWoDCA, Int_t multRecWoDCA) {
  Int_t multTSrec = 0;
  // see how many tracks there are in the event
  for (Int_t i = 0; i < multRec; ++i) { // loop over all these tracks

    if (i == fRecLeadIn)
      continue;
    Double_t DPhi = DeltaPhi(phiRec[i], fRecLeadPhi);

    // definition of the topological regions
    if (TMath::Abs(DPhi) < pi / 3.0) { // near side
      continue;
    } else if (TMath::Abs(DPhi - pi) < pi / 3.0) { // away side
      continue;
    } else { // transverse side
      multTSrec++;
    }
  }

  hNchTSRecTest->Fill(multTSrec);
  hNchRecTest->Fill(multRec);

  // Filling rec pT vs UE (for pT *** considering the hybrid track cuts***)

  for (Int_t i = 0; i < multRec; ++i) { // loop over all these tracks

    if (i == fRecLeadIn)
      continue;
    Double_t DPhi = DeltaPhi(phiRec[i], fRecLeadPhi);

    // definition of the topological regions
    if (TMath::Abs(DPhi) < pi / 3.0) { // near side
      hPtVsUERecTest[0]->Fill(multTSrec, ptRec[i]);
      hPtVsNchRecTest[0]->Fill(multRec, ptRec[i]);
    } else if (TMath::Abs(DPhi - pi) < pi / 3.0) { // away side
      hPtVsUERecTest[1]->Fill(multTSrec, ptRec[i]);
      hPtVsNchRecTest[1]->Fill(multRec, ptRec[i]);
    } else { // transverse side
      hPtVsUERecTest[2]->Fill(multTSrec, ptRec[i]);
      hPtVsNchRecTest[2]->Fill(multRec, ptRec[i]);
    }
  }

  // Auxiliar distributions to calculate the contamination from secondary
  // particles, it runs over tracks wo DCA cut
  for (Int_t i = 0; i < multRecWoDCA; ++i) { // loop over all these tracks

    hptvsdcaAll->Fill(ptRecWoDCA[i], dcaxyRecWoDCA[i]);

    if (fUseMC && (fGenLeadPt >= fLeadPtCutMin && fGenLeadPt < fLeadPtCutMax)) {
      if (isprimRecWoDCA[i] == 0) {
        hptvsdcaPrim->Fill(ptRecWoDCA[i], dcaxyRecWoDCA[i]);
      } else if (isprimRecWoDCA[i] == 1) {
        hptvsdcaDecs->Fill(ptRecWoDCA[i], dcaxyRecWoDCA[i]);
      } else if (isprimRecWoDCA[i] == 2) {
        hptvsdcaMatl->Fill(ptRecWoDCA[i], dcaxyRecWoDCA[i]);
      }
    }
  }
}

Double_t AliAnalysisTaskChargedVsRT::DeltaPhi(Double_t phia, Double_t phib,
                                              Double_t rangeMin,
                                              Double_t rangeMax) {
  Double_t dphi = -999;
  Double_t pi = TMath::Pi();

  if (phia < 0)
    phia += 2 * pi;
  else if (phia > 2 * pi)
    phia -= 2 * pi;
  if (phib < 0)
    phib += 2 * pi;
  else if (phib > 2 * pi)
    phib -= 2 * pi;
  dphi = phib - phia;
  if (dphi < rangeMin)
    dphi += 2 * pi;
  else if (dphi > rangeMax)
    dphi -= 2 * pi;

  return dphi;
}
Bool_t AliAnalysisTaskChargedVsRT::HasRecVertex() {

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
//______________________________________________________________________
Int_t AliAnalysisTaskChargedVsRT::FillArrayMC(vector<Float_t> &ptArray,
                                              vector<Float_t> &phiArray,
                                              vector<Int_t> &idArray) {
  /*
     id 0: pion, 1: kaon, 2: proton, 3: sigma plus, 4: sigma minus, 5: Omega, 6:
     Xi, 7: other charged
   */

  ptArray.clear();
  phiArray.clear();
  idArray.clear();
  Int_t nNchGen = 0;

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

    Int_t idPart = -1;
    Int_t partPDG = TMath::Abs(particle->PdgCode());
    if (partPDG == 3122)
      idPart = 0; // lambda
    if (particle->Charge() != 0) {
      if (partPDG == 211)
        idPart = 1; // pions
      else if (partPDG == 321)
        idPart = 2; // kaons
      else if (partPDG == 2212)
        idPart = 3; // protons
      else if (partPDG == 3222)
        idPart = 4; // sigma plus
      else if (partPDG == 3112)
        idPart = 5; // sigma minus
      else if (partPDG == 3334)
        idPart = 6; // Omega
      else if (partPDG == 3312)
        idPart = 7; // Xi
      else
        idPart = 8; // rest of the charged particles
    }
    ptArray.push_back(particle->Pt());
    phiArray.push_back(particle->Phi());
    idArray.push_back(idPart);
    nNchGen++;
  }
  return nNchGen;
}
//_____________________________
Int_t AliAnalysisTaskChargedVsRT::FillArray(
    vector<Float_t> &ptArray, vector<Float_t> &phiArray,
    vector<Float_t> &dcaxyArray, vector<Float_t> &dcazArray,
    vector<Int_t> &isprimArray, vector<Int_t> &idArray, const bool wDcaCut,
    const bool useHy) {
  /*
     id 0: pion, 1: kaon, 2: proton, 3: sigma plus, 4: sigma minus, 5: Omega, 6:
     Xi, 7: other charged
   */
  ptArray.clear();
  phiArray.clear();
  dcaxyArray.clear();
  dcazArray.clear();
  isprimArray.clear();
  idArray.clear();

  Int_t nTracks = fESD->GetNumberOfTracks();
  Int_t nNchRec = 0;
  if (wDcaCut) { // with DCA cut
    for (Int_t iT = 0; iT < nTracks; ++iT) {

      AliESDtrack *esdtrack = static_cast<AliESDtrack *>(
          fESD->GetTrack(iT)); // get a track (type AliesdTrack)
      if (!esdtrack)
        continue;
      fdcaxy = -999;
      fdcaz = -999;
      if (TMath::Abs(esdtrack->Eta()) > fEtaCut)
        continue;
      if (esdtrack->Pt() < fPtMin)
        continue;
      AliESDtrack *newTrack = 0x0;
      Int_t isPrim = -1;
      Int_t idTrack = -1;
      Int_t mcLabel = -1;
      if (fUseMC) { // get label: 0: prim, 1: weak decays, 2: material
        mcLabel = TMath::Abs(esdtrack->GetLabel());
        TParticle *mcParticle = fMC->GetTrack(mcLabel)->Particle();
        if (!mcParticle) {
          printf("----ERROR: mcParticle not available------------------\n");
          continue;
        }
        Int_t partPDG_rec = TMath::Abs(mcParticle->GetPdgCode());
        if (partPDG_rec == 3122)
          idTrack = 0; // lambdas
        else if (partPDG_rec == 211)
          idTrack = 1; // pions
        else if (partPDG_rec == 321)
          idTrack = 2; // kaons
        else if (partPDG_rec == 2212)
          idTrack = 3; // protons
        else if (partPDG_rec == 3222)
          idTrack = 4; // sigma plus
        else if (partPDG_rec == 3112)
          idTrack = 5; // sigma minus
        else if (partPDG_rec == 3334)
          idTrack = 6; // Omega
        else if (partPDG_rec == 3312)
          idTrack = 7; // Xi
        else
          idTrack = 8; // rest of the charged particles
      }
      Bool_t isHy0 = kFALSE;
      Bool_t isHy1 = kFALSE;
      if (useHy) {
        if (fTrackFilterHybrid0->IsSelected(esdtrack))
          isHy0 = kTRUE;
        else if (fTrackFilterHybrid1->IsSelected(esdtrack))
          isHy1 = kTRUE;
      } else {
        if (fTrackFilter->IsSelected(esdtrack))
          isHy0 = kTRUE;
      }

      if (isHy0) {
        newTrack = new AliESDtrack(*esdtrack);
        hPhiTotal->Fill(esdtrack->Phi());
        hPhiStandard->Fill(newTrack->Phi());
        newTrack->GetImpactParameters(fdcaxy, fdcaz);
        ptArray.push_back(newTrack->Pt());
        phiArray.push_back(newTrack->Phi());
        dcaxyArray.push_back(fdcaxy);
        dcazArray.push_back(fdcaz);
        if (fUseMC) { // get label: 0: prim, 1: weak decays, 2: material
          if (fMC->IsPhysicalPrimary(mcLabel))
            isPrim = 0;
          else if (fMC->IsSecondaryFromWeakDecay(mcLabel))
            isPrim = 1;
          else if (fMC->IsSecondaryFromMaterial(mcLabel))
            isPrim = 2;
          else
            continue;
        }
        isprimArray.push_back(isPrim);
        idArray.push_back(idTrack);
        nNchRec++;
      } else if (isHy1) {
        newTrack = new AliESDtrack(*esdtrack);
        if (esdtrack->GetConstrainedParam()) {
          const AliExternalTrackParam *constrainParam =
              esdtrack->GetConstrainedParam();
          newTrack->Set(constrainParam->GetX(), constrainParam->GetAlpha(),
                        constrainParam->GetParameter(),
                        constrainParam->GetCovariance());
          hPhiTotal->Fill(esdtrack->Phi());
          hPhiHybrid1->Fill(newTrack->Phi());
          newTrack->GetImpactParameters(fdcaxy, fdcaz);
          ptArray.push_back(newTrack->Pt());
          phiArray.push_back(newTrack->Phi());
          dcaxyArray.push_back(fdcaxy);
          dcazArray.push_back(fdcaz);
          if (fUseMC) { // get label: 0: prim, 1: weak decays, 2: material
            if (fMC->IsPhysicalPrimary(mcLabel))
              isPrim = 0;
            else if (fMC->IsSecondaryFromWeakDecay(mcLabel))
              isPrim = 1;
            else if (fMC->IsSecondaryFromMaterial(mcLabel))
              isPrim = 2;
            else
              continue;
          }
          isprimArray.push_back(isPrim);
          idArray.push_back(idTrack);
          nNchRec++;
        }
      }

      delete newTrack;
    }
  } else {
    for (Int_t iT = 0; iT < nTracks; ++iT) {

      AliESDtrack *esdtrack = static_cast<AliESDtrack *>(
          fESD->GetTrack(iT)); // get a track (type AliesdTrack)
      if (!esdtrack)
        continue;
      fdcaxy = -999;
      fdcaz = -999;
      if (TMath::Abs(esdtrack->Eta()) > fEtaCut)
        continue;
      if (esdtrack->Pt() < fPtMin)
        continue;
      AliESDtrack *newTrack = 0x0;
      Int_t isPrim = -1;
      Int_t idTrack = -1;
      Int_t mcLabel = -1;
      TParticle *mcParticle = 0;
      if (fUseMC) { // get label: 0: prim, 1: weak decays, 2: material
        mcLabel = TMath::Abs(esdtrack->GetLabel());
        mcParticle = fMC->GetTrack(mcLabel)->Particle();
        if (!mcParticle) {
          printf("----ERROR: mcParticle not available------------------\n");
          continue;
        }

        Int_t partPDG_rec = TMath::Abs(mcParticle->GetPdgCode());
        if (partPDG_rec == 3122)
          idTrack = 0; // lambdas
        else if (partPDG_rec == 211)
          idTrack = 1; // pions
        else if (partPDG_rec == 321)
          idTrack = 2; // kaons
        else if (partPDG_rec == 2212)
          idTrack = 3; // protons
        else if (partPDG_rec == 3222)
          idTrack = 4; // sigma plus
        else if (partPDG_rec == 3112)
          idTrack = 5; // sigma minus
        else if (partPDG_rec == 3334)
          idTrack = 6; // Omega
        else if (partPDG_rec == 3312)
          idTrack = 7; // Xi
        else
          idTrack = 8; // rest of the charged particles
      }
      Bool_t isHy0 = kFALSE;
      Bool_t isHy1 = kFALSE;
      if (useHy) {
        if (fTrackFilterHybrid0woDCA->IsSelected(esdtrack))
          isHy0 = kTRUE;
        else if (fTrackFilterHybrid1woDCA->IsSelected(esdtrack))
          isHy1 = kTRUE;
      } else {
        if (fTrackFilterwoDCA->IsSelected(esdtrack))
          isHy0 = kTRUE;
      }
      if (isHy0) {
        newTrack = new AliESDtrack(*esdtrack);
        newTrack->GetImpactParameters(fdcaxy, fdcaz);
        ptArray.push_back(newTrack->Pt());
        phiArray.push_back(newTrack->Phi());
        dcaxyArray.push_back(fdcaxy);
        dcazArray.push_back(fdcaz);
        if (fUseMC) { // get label: 0: prim, 1: weak decays, 2: material
          if (fMC->IsPhysicalPrimary(mcLabel))
            isPrim = 0;
          else if (fMC->IsSecondaryFromWeakDecay(mcLabel))
            isPrim = 1;
          else if (fMC->IsSecondaryFromMaterial(mcLabel))
            isPrim = 2;
          else
            continue;
        }
        isprimArray.push_back(isPrim);
        idArray.push_back(idTrack);
        nNchRec++;
      } else if (isHy1) {
        newTrack = new AliESDtrack(*esdtrack);
        if (esdtrack->GetConstrainedParam()) {
          const AliExternalTrackParam *constrainParam =
              esdtrack->GetConstrainedParam();
          newTrack->Set(constrainParam->GetX(), constrainParam->GetAlpha(),
                        constrainParam->GetParameter(),
                        constrainParam->GetCovariance());
          newTrack->GetImpactParameters(fdcaxy, fdcaz);
          ptArray.push_back(newTrack->Pt());
          phiArray.push_back(newTrack->Phi());
          dcaxyArray.push_back(fdcaxy);
          dcazArray.push_back(fdcaz);
          if (fUseMC) { // get label: 0: prim, 1: weak decays, 2: material
            if (fMC->IsPhysicalPrimary(mcLabel))
              isPrim = 0;
            else if (fMC->IsSecondaryFromWeakDecay(mcLabel))
              isPrim = 1;
            else if (fMC->IsSecondaryFromMaterial(mcLabel))
              isPrim = 2;
            else
              continue;
          }
          isprimArray.push_back(isPrim);
          idArray.push_back(idTrack);
          nNchRec++;
        }
      }

      delete newTrack;
    }
  }
  return nNchRec;
}
//________________________________________________________________________
void AliAnalysisTaskChargedVsRT::SetCutsHybrid0WoDCA(AliESDtrackCuts *cIn) {
  cIn->SetMinNCrossedRowsTPC(70); // TBC
  cIn->SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);
  cIn->SetMaxChi2PerClusterTPC(4);
  cIn->SetMaxDCAToVertexZ(2);

  cIn->SetAcceptKinkDaughters(kFALSE);
  cIn->SetRequireTPCRefit(kTRUE);
  cIn->SetRequireITSRefit(kTRUE);
  cIn->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kOff);
  cIn->SetDCAToVertex2D(kFALSE);
  cIn->SetRequireSigmaToVertex(kFALSE);
  cIn->SetEtaRange(-0.8, 0.8);

  fTrackFilterHybrid0woDCA->AddCuts(cIn);
  // fTrackFilterHybrid1woDCA->AddCuts(cIn);
}
//________________________________________________________________________
void AliAnalysisTaskChargedVsRT::SetCutsHybrid1WoDCA(AliESDtrackCuts *cIn1) {
  cIn1->SetMinNCrossedRowsTPC(70); // TBC
  cIn1->SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);
  cIn1->SetMaxChi2PerClusterTPC(4);
  cIn1->SetMaxDCAToVertexZ(2);

  cIn1->SetAcceptKinkDaughters(kFALSE);
  cIn1->SetRequireTPCRefit(kTRUE);
  cIn1->SetRequireITSRefit(kFALSE);
  cIn1->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kNone);
  cIn1->SetDCAToVertex2D(kFALSE);
  cIn1->SetRequireSigmaToVertex(kFALSE);
  cIn1->SetEtaRange(-0.8, 0.8);

  // fTrackFilterHybrid0woDCA->AddCuts(cIn);
  fTrackFilterHybrid1woDCA->AddCuts(cIn1);
}
//________________________________________________________________________
void AliAnalysisTaskChargedVsRT::SetCutsFilterWoDCA(
    AliESDtrackCuts *cFilterIn) {
  cFilterIn->SetMinNCrossedRowsTPC(70); // TBC
  cFilterIn->SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);
  cFilterIn->SetMaxChi2PerClusterTPC(4);
  cFilterIn->SetMaxDCAToVertexZ(2);

  cFilterIn->SetAcceptKinkDaughters(kFALSE);
  cFilterIn->SetRequireTPCRefit(kTRUE);
  cFilterIn->SetRequireITSRefit(kTRUE);
  cFilterIn->SetClusterRequirementITS(AliESDtrackCuts::kSPD,
                                      AliESDtrackCuts::kAny);
  cFilterIn->SetDCAToVertex2D(kFALSE);
  cFilterIn->SetRequireSigmaToVertex(kFALSE);
  cFilterIn->SetEtaRange(-0.8, 0.8);

  if (fGeoTPCVar1) {
    cFilterIn->SetCutGeoNcrNcl(2., 130., 1.5, 0.85, 0.7);
  } //
  else if (fGeoTPCVar2) {
    cFilterIn->SetCutGeoNcrNcl(4., 130., 1.5, 0.85, 0.7);
  } //
  else if (fGeoTPCVar3) {
    cFilterIn->SetCutGeoNcrNcl(3., 120., 1.5, 0.85, 0.7);
  } //
  else if (fGeoTPCVar4) {
    cFilterIn->SetCutGeoNcrNcl(3., 140., 1.5, 0.85, 0.7);
  } //
  else {
    cFilterIn->SetCutGeoNcrNcl(3., 130., 1.5, 0.85, 0.7);
  } // Default

  if (fNcrVar1) {
    cFilterIn->SetMinRatioCrossedRowsOverFindableClustersTPC(0.7);
  } //
  else if (fNcrVar2) {
    cFilterIn->SetMinRatioCrossedRowsOverFindableClustersTPC(0.9);
  } //
  else {
    cFilterIn->SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);
  } // Default

  if (fChisqTPCVar1) {
    cFilterIn->SetMaxChi2PerClusterTPC(3);
  } else if (fChisqTPCVar2) {
    cFilterIn->SetMaxChi2PerClusterTPC(5);
  } else {
    cFilterIn->SetMaxChi2PerClusterTPC(4);
  } // Default

  if (fDcazVar1) {
    cFilterIn->SetMaxDCAToVertexZ(1);
  } // DCAz = 1 cm
  else if (fDcazVar2) {
    cFilterIn->SetMaxDCAToVertexZ(5);
  } // DCAz = 5 cm
  else {
    cFilterIn->SetMaxDCAToVertexZ(2);
  } // Default

  fTrackFilterwoDCA->AddCuts(cFilterIn);
}
