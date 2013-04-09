/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, proviyaded that the above copyright notice appears in all *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purapose. It is         *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

///////////////////////////////////////////////////////////////////////////
//                                                                       //
//                                                                       //
// Analysis for identified charged hadron spectra in pPb.                //
//                                                                       //
//                                                                       //
///////////////////////////////////////////////////////////////////////////

#include "Riostream.h"
#include "TChain.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2D.h"
#include "TH3F.h"
#include "TList.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TObjArray.h"
#include "TF1.h"
#include "TFile.h"

#include "AliAnalysisTaskSE.h"
#include "AliAnalysisManager.h"

#include "AliHeader.h"
#include "AliGenPythiaEventHeader.h"
#include "AliGenCocktailEventHeader.h"

#include "AliPID.h"
#include "AliESDtrackCuts.h"
#include "AliESDVertex.h"
#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliESDtrack.h"
#include "AliESDpid.h"
#include "AliCentrality.h"
#include "AliESDUtils.h"
#include "AliMultiplicity.h"

#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliStack.h"

#include "AliLog.h"                 

#include "AliAnalysisTPCTOFpA.h"


ClassImp(AliAnalysisTPCTOFpA)

//________________________________________________________________________
AliAnalysisTPCTOFpA::AliAnalysisTPCTOFpA() 
  : AliAnalysisTaskSE("TaskChargedHadron"), fESD(0), fListHist(0), fESDtrackCuts(0),fESDTrackCutsMult(0),fESDpid(0),
  fMCtrue(0),
  fOnlyQA(0),
  fUseHBTmultiplicity(0),
  fUseTPConlyTracks(0),
  fSaveMotherPDG(0),
  fSmallTHnSparse(0),
  fIspA(0),
  fRapCMS(0),
  fCentEst(0),
  fTOFmisMatch(0),
  fTRDinReject(0),
  fTOFwindow(0),
  fDCAzCut(0),
  fCrossedRows(0),
  fRatioRowsClusters(0),
  fTPCnSigmaCutLow(0),
  fTPCnSigmaCutHigh(0),
  fRapidityCutLow(0),
  fRapidityCutHigh(0),
  fEvenDCAbinning(0),
  fAlephParameters(),
  fHistRealTracks(0),
  fHistMCparticles(0),
  fHistPidQA(0),
  fHistMult(0),
  fHistCentrality(0),
  fHistTOFwindow(0)
{
  // default Constructor
  /* fast compilation test
     gSystem->Load("libANALYSIS");
     gSystem->Load("libANALYSISalice");
     .L /d/alice09/akalweit/train/trunk/akalweit_hadronSpectra/AliAnalysisTPCTOFpA.cxx++
   */
}


//________________________________________________________________________
AliAnalysisTPCTOFpA::AliAnalysisTPCTOFpA(const char *name) 
  : AliAnalysisTaskSE(name), fESD(0), fListHist(0), fESDtrackCuts(0),fESDTrackCutsMult(0),fESDpid(0),
    fMCtrue(0),
    fOnlyQA(0),
    fUseHBTmultiplicity(0),
    fUseTPConlyTracks(0),
    fSaveMotherPDG(0),
    fSmallTHnSparse(0),
    fIspA(0),
    fRapCMS(0),
    fCentEst(0),
    fTOFmisMatch(0),
    fTRDinReject(0),
    fTOFwindow(0),
    fDCAzCut(0),
    fCrossedRows(0),
    fRatioRowsClusters(0),
    fTPCnSigmaCutLow(0),
    fTPCnSigmaCutHigh(0),
    fRapidityCutLow(0),
    fRapidityCutHigh(0),
    fEvenDCAbinning(0),
    fAlephParameters(),
    fHistRealTracks(0),
    fHistMCparticles(0),
    fHistPidQA(0),
    fHistMult(0),
    fHistCentrality(0),
    fHistTOFwindow(0)
    {
  //
  // standard constructur which should be used
  //
  Printf("*** CONSTRUCTOR CALLED ****");
  //
  fMCtrue = kTRUE; 
  fOnlyQA = kFALSE;
  fUseHBTmultiplicity = kTRUE;
  fUseTPConlyTracks = kFALSE;


  fUseTPConlyTracks =kFALSE;
  fSaveMotherPDG = kFALSE;
  fSmallTHnSparse = kTRUE;
  fTPCnSigmaCutLow = -3.;
  fTPCnSigmaCutHigh = 3.;
  fRapidityCutLow = 0.;
  fRapidityCutHigh = 0.5;
  fEvenDCAbinning = kFALSE;
  fIspA = kTRUE;
  fRapCMS = kTRUE;
  fCentEst = "V0A";
  fEvenDCAbinning = kFALSE;
  fTOFmisMatch = 2;
  fTRDinReject = kFALSE;
  fTOFwindow = 10.;
  fDCAzCut = 2.;
  fCrossedRows = 70;
  fRatioRowsClusters = 0.8;



  /* real */
  fAlephParameters[0] = 0.0283086;
  fAlephParameters[1] = 2.63394e+01;
  fAlephParameters[2] = 5.04114e-11;
  fAlephParameters[3] = 2.12543e+00;
  fAlephParameters[4] = 4.88663e+00;
  //
  // initialize PID object
  //
  //fESDpid = new AliESDpid();
  //
  // create track cuts
  //
  fESDtrackCuts = new AliESDtrackCuts("AliESDtrackCuts","AliESDtrackCuts");
  //
  //Initialize();
  // Output slot #0 writes into a TList container
  DefineOutput(1, TList::Class());

}


//________________________________________________________________________
void AliAnalysisTPCTOFpA::Initialize()
{
  //
  // updating parameters in case of changes
  //
  // 1. pid parameters
  //
  //fESDpid->GetTPCResponse().SetBetheBlochParameters(fAlephParameters[0],fAlephParameters[1],fAlephParameters[2],fAlephParameters[3],fAlephParameters[4]);
  //
  // 2. track cuts
  //
  /*
  fESDtrackCuts->SetMaxCovDiagonalElements(2, 2, 0.5, 0.5, 2);  // BEWARE STANDARD VALUES ARE: 2, 2, 0.5, 0.5, 2
  fESDtrackCuts->SetMaxNsigmaToVertex(3);
  fESDtrackCuts->SetRequireSigmaToVertex(kTRUE);
  fESDtrackCuts->SetAcceptKinkDaughters(kFALSE);
  fESDtrackCuts->SetMinNClustersTPC(70);
  fESDtrackCuts->SetMaxChi2PerClusterTPC(4);
  fESDtrackCuts->SetMaxDCAToVertexXY(3);
  fESDtrackCuts->SetMaxDCAToVertexZ(3);
  fESDtrackCuts->SetRequireTPCRefit(kTRUE);
  fESDtrackCuts->SetRequireITSRefit(kTRUE);
  fESDtrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kAny); //TEMPORARY <-> REMOVE
  fESDtrackCuts->SetMinNClustersITS(3);
  */
  //fESDtrackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2010(kTRUE); // kTRUE = sel. primaries --> patch for the moment, do TFractionFitter later


  if (!fUseTPConlyTracks) {
    fESDtrackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011(kFALSE,kTRUE);
	  fESDtrackCuts->SetMaxDCAToVertexXY(3);
	  //fESDtrackCuts->SetMaxDCAToVertexZ(2);
	  fESDtrackCuts->SetEtaRange(-0.8,0.8);
  }
  else {
	  //fESDtrackCuts = AliESDtrackCuts::GetStandardTPCOnlyTrackCuts();
	  fESDtrackCuts->SetMinNClustersTPC(70);
	  fESDtrackCuts->SetMaxChi2PerClusterTPC(4);
	  fESDtrackCuts->SetAcceptKinkDaughters(kFALSE);
	  fESDtrackCuts->SetRequireTPCRefit(kFALSE);

	  fESDtrackCuts->SetMaxDCAToVertexXY(15);
	  fESDtrackCuts->SetMaxDCAToVertexZ(6);
	  fESDtrackCuts->SetDCAToVertex2D(kFALSE);
	  fESDtrackCuts->SetRequireSigmaToVertex(kFALSE);

	  fESDtrackCuts->SetEtaRange(-0.9,0.9);
  }


  //change DCA z cut here with flags set before:
  fESDtrackCuts->SetMaxDCAToVertexZ(fDCAzCut);
  fESDtrackCuts->SetMinNCrossedRowsTPC(fCrossedRows);
  fESDtrackCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(fRatioRowsClusters);

  //
  //
  //
  //
  fESDTrackCutsMult = AliESDtrackCuts::GetStandardITSTPCTrackCuts2010(kTRUE);
  fESDTrackCutsMult->SetEtaRange(-1.2,+1.2);
  fESDTrackCutsMult->SetPtRange(0.15,1e10);

}


//________________________________________________________________________
void AliAnalysisTPCTOFpA::UserCreateOutputObjects() 
{
  // Create histograms
  // Called once
  fListHist = new TList();
  fListHist->SetOwner(kTRUE);
  //
  const Int_t kPtBins = 35;
  const Int_t kMultBins = 11;

  Int_t kDcaBinsTemp = 76;
  if (fEvenDCAbinning) kDcaBinsTemp = 150;
  const Int_t kDcaBins = (const Int_t) kDcaBinsTemp;

  const Float_t kDcaBinsTPConlyFactor = 5; //need to change binning of DCA plot for tpconly
  // sort pT-bins ..
  Double_t binsPt[77] = {0., 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0};

  Double_t binsDca[77] =  {-3,-2.85,-2.7,-2.55,-2.4,-2.25,-2.1,-1.95,-1.8,-1.65,-1.5,-1.35,-1.2,-1.05,-0.9,-0.75,-0.6,-0.45,-0.3,-0.285,-0.27,-0.255,-0.24,-0.225,-0.21,-0.195,-0.18,-0.165,-0.15,-0.135,-0.12,-0.105,-0.09,-0.075,-0.06,-0.045,-0.03,-0.015,0,0.015,0.03,0.045,0.06,0.075,0.09,0.105,0.12,0.135,0.15,0.165,0.18,0.195,0.21,0.225,0.24,0.255,0.27,0.285,0.3,0.45,0.6,0.75,0.9,1.05,1.2,1.35,1.5,1.65,1.8,1.95,2.1,2.25,2.4,2.55,2.7,2.85,3};

  // DCA bins borders get multiplied by constant factor for TPConlyTracks
  Double_t binsDcaTPConly[77];
  for (Int_t i = 0; i< 77; i++) {
  	binsDcaTPConly[i] = kDcaBinsTPConlyFactor * binsDca[i];
  }

  //
  // create the histograms with all necessary information --> it is filled 4x for each particle assumption
  //
  // (0.) assumed particle: 0. pion, 1. kaon, 2. proton, 3. deuteron
  // (1.) multiplicity or centrality -- number of accepted ESD tracks per events (deprecated), but now classes from 1 to 10, 0: Min. Bias
  // (2.) pT
  // (3.) sign
  // (4.) rapidity --> filled 4x
  // (5.)  pull TPC dEx --> filled 4x
  // (6.) has valid TOF pid signal
  // (7.) nsigma TOF --> filled 4x
  // (8..) dca_xy
  // (9.) CODE -- only MC 0-generated, 1-true rec. primaries, 2-misident prim, 3-second weak, 4-second material, 5-misident sec, 6-sec. K0, 7-sec. lambda, 8-sec sigma+
  //
  
  //dimensions of standard THnSparse
  //                              0,           1,           2,  3,       4,    5,     6,   7,     8
  Int_t    binsHistReal[9] = {   3,   kMultBins,     kPtBins,  2,       20,    50,    2,  16,    kDcaBins};
  Double_t xminHistReal[9] = {-0.5,        -1.5,           0, -2,      -1.0,   -5, -0.5,  -8,       -3};
  Double_t xmaxHistReal[9] = { 2.5,         9.5,           3,  2,       1.0,    5,  1.5,   8,        3};

  //dimensions of small THnSparse
  //                              0,           1,           2,  3,                        4,   5,       6
  Int_t    binsHistRealSm[7] = {   3,   kMultBins,     kPtBins,  2,  /*    10,   50,*/    2,  80, kDcaBins};
  Double_t xminHistRealSm[7] = {-0.5,        -1.5,           0, -2,  /*  -0.5,   -5,*/ -0.5,  -8,       -3};
  Double_t xmaxHistRealSm[7] = { 2.5,         9.5,           3,  2,  /*   0.5,    5,*/  1.5,   8,        3};

  if (!fSmallTHnSparse) fHistRealTracks = new THnSparseF("fHistRealTracks","real tracks",9,binsHistReal,xminHistReal,xmaxHistReal);
  else                  fHistRealTracks = new THnSparseF("fHistRealTracks","real tracks",7,binsHistRealSm,xminHistRealSm,xmaxHistRealSm);
  //
  fHistRealTracks->GetAxis(2)->Set(kPtBins, binsPt);

  //different DCAxy binning for TPConlyTracks

  Int_t dcaAxisNumber = 8;
  if (fSmallTHnSparse) dcaAxisNumber = 6;

  if (!fUseTPConlyTracks){
    if (fEvenDCAbinning) fHistRealTracks->GetAxis(dcaAxisNumber)->Set(kDcaBins,-3.,3.);
    else fHistRealTracks->GetAxis(dcaAxisNumber)->Set(kDcaBins, binsDca);
  }
  else {
    if (fEvenDCAbinning) fHistRealTracks->GetAxis(dcaAxisNumber)->Set(kDcaBins,-15.,15.);
    else fHistRealTracks->GetAxis(dcaAxisNumber)->Set(kDcaBins, binsDcaTPConly);
   }
  fListHist->Add(fHistRealTracks);
  //
  //                      0.ptot,1.tpcSig,2.hasTOF, 3. assumed part., 4. nclDedx, 5. nSigmaTPC (4x), 6. nSigmaTOF (4x), 7. centrality
  fHistPidQA = new TH3D("fHistPidQA","PID QA",500,0.1,10,1000,0,1000,2,-2,2);
  BinLogAxis(fHistPidQA);
  fListHist->Add(fHistPidQA);
  
  // dimensions of standard THnSparse
  //                            0,            1,           2,  3,      4,   5,    6,   7,       8.,    9.
  Int_t    binsHistMC[10] = {   3,    kMultBins,     kPtBins,  2,     20,  50,    2,  16, kDcaBins,    6};
  Double_t xminHistMC[10] = {-0.5,         -1.5,           0, -2,   -1.0,  -5, -0.5,  -8,       -3, -0.5};
  Double_t xmaxHistMC[10] = { 2.5,          9.5,           3,  2,    1.0,   5,  1.5,   8,        3,  5.5};

  //dimensions of small THnSparse
  //                              0,           1,           2,  3,                     4,   5,       6.,    7.
  Int_t    binsHistMCSm[8] = {   3,    kMultBins,     kPtBins,  2,  /*   10,  50,*/    2,  80, kDcaBins,    6};
  Double_t xminHistMCSm[8] = {-0.5,         -1.5,           0, -2,  /* -0.5,  -5,*/ -0.5,  -8,       -3, -0.5};
  Double_t xmaxHistMCSm[8] = { 2.5,          9.5,           3,  2,  /*  0.5,   5,*/  1.5,   8,        3,  5.5};

  //different binning for CODE axis, if we want to save motherPDG
  if (fSaveMotherPDG) {
    binsHistMC[9] = 9;
    xmaxHistMC[9] = 8.5;
    binsHistMCSm[7] = 9;
    xmaxHistMCSm[7] = 8.5;
  }


  if (!fSmallTHnSparse) fHistMCparticles = new THnSparseF("fHistMCparticles","MC histogram",10,binsHistMC,xminHistMC,xmaxHistMC);
  else                  fHistMCparticles = new THnSparseF("fHistMCparticles","MC histogram",8,binsHistMCSm,xminHistMCSm,xmaxHistMCSm);


  fHistMCparticles->GetAxis(2)->Set(kPtBins, binsPt);

  /*
  //different DCAxy binning for TPConlyTracks
  if (!fEvenDCAbinning){
    if (!fUseTPConlyTracks) fHistMCparticles->GetAxis(dcaAxisNumber)->Set(kDcaBins, binsDca);
    else fHistMCparticles->GetAxis(dcaAxisNumber)->Set(kDcaBins, binsDcaTPConly);
  }
  */

  if (!fUseTPConlyTracks){
    if (fEvenDCAbinning) fHistMCparticles->GetAxis(dcaAxisNumber)->Set(kDcaBins,-3.,3.);
    else fHistMCparticles->GetAxis(dcaAxisNumber)->Set(kDcaBins, binsDca);
  }
  else {
    if (fEvenDCAbinning) fHistMCparticles->GetAxis(dcaAxisNumber)->Set(kDcaBins,-15.,15.);
    else fHistMCparticles->GetAxis(dcaAxisNumber)->Set(kDcaBins, binsDcaTPConly);
   }



  fListHist->Add(fHistMCparticles);
  //
  fHistMult = new TH2D("fHistMult", "control histogram to count number of events", 502, -2.5, 499.5,4,-0.5,3.5);
  fHistCentrality = new TH1D("fHistCentrality", "control histogram to count number of events", 22, -1.5, 20.5);
  fHistTOFwindow = new TH2D("fHistTOFwindow", "control hisogram for TOF window",160,-10.,10.,160,-10.,10.);
  fHistTOFwindow->GetXaxis()->SetTitle("dx");
  fHistTOFwindow->GetYaxis()->SetTitle("dz");
  fListHist->Add(fHistMult);
  fListHist->Add(fHistCentrality);
  fListHist->Add(fHistTOFwindow);

  PostData(1, fListHist);

}

//________________________________________________________________________
void AliAnalysisTPCTOFpA::UserExec(Option_t *) 
{
  //
  // main event loop
  //
  if (!fESDpid) fESDpid = ((AliESDInputHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->GetESDpid();
  if (!fESDpid) {
    fESDpid = new AliESDpid(); // HACK FOR MC PBPB --> PLEASE REMOVE AS SOON AS POSSIBLE
    fESDpid->GetTPCResponse().SetBetheBlochParameters(1.28778e+00/50., 3.13539e+01, TMath::Exp(-3.16327e+01), 1.87901e+00, 6.41583e+00);
  }
  //AliLog::SetGlobalLogLevel(AliLog::kError);
  //
  // Check Monte Carlo information and other access first:
  //
  AliMCEventHandler* eventHandler = dynamic_cast<AliMCEventHandler*> (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler());
  if (!eventHandler) {
    //Printf("ERROR: Could not retrieve MC event handler");
    fMCtrue = kFALSE;
  }
  //
  AliMCEvent* mcEvent = 0x0;
  AliStack* stack = 0x0;
  if (eventHandler) mcEvent = eventHandler->MCEvent();
  if (!mcEvent) {
    //Printf("ERROR: Could not retrieve MC event");
    if (fMCtrue) return;
  }
  if (fMCtrue) {
    stack = mcEvent->Stack();
    if (!stack) return;
  }
  //
  fESD = dynamic_cast<AliESDEvent*>( InputEvent() );
  if (!fESD) {
    //Printf("ERROR: fESD not available");
    return;
  }
  
  if (!fESDtrackCuts) {
    Printf("ERROR: fESDtrackCuts not available");
    return;
  }
  //
  // check if event is selected by physics selection class
  //
  Bool_t isSelected = kTRUE; // for reasons of backward compatibility --> check is now in AddTask macro
  /*
  Bool_t isSelected = kFALSE;
  isSelected = (((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & AliVEvent::kMB) == AliVEvent::kMB;
  */
  //
  // monitor vertex position
  //
  Bool_t isVertexOk = kTRUE;
  const AliESDVertex *vertex = fESD->GetPrimaryVertexTracks();
  if(vertex->GetNContributors()<1) {
    // SPD vertex
    vertex = fESD->GetPrimaryVertexSPD();
    /* quality checks on SPD-vertex */ 
    TString vertexType = vertex->GetTitle();
    if (vertexType.Contains("vertexer: Z") && (vertex->GetDispersion() > 0.04 || vertex->GetZRes() > 0.25))  isVertexOk = kFALSE; //vertex = 0x0; //
    if (vertex->GetNContributors()<1)  isVertexOk = kFALSE; //vertex = 0x0; //
  }  
  //
  // small track loop to determine trigger Pt, multiplicity or centrality
  //
  Double_t triggerPt = 0;
  for (Int_t i=0;i<fESD->GetNumberOfTracks();++i) {
    AliESDtrack *track =fESD->GetTrack(i);
    if (!fESDTrackCutsMult->AcceptTrack(track)) continue;
    if (track->Pt() > triggerPt) triggerPt = track->Pt();
  }
  Int_t trackCounter = fESDTrackCutsMult->CountAcceptedTracks(fESD);
  //
  // 2nd multiplicity estimator SPD hits; replaces multiplicity from HBT
  //
  Float_t spdCorr = -1;
  const AliMultiplicity *mult = fESD->GetMultiplicity();
  Float_t nClusters[6]={0.0,0.0,0.0,0.0,0.0,0.0};
  for(Int_t ilay=0; ilay<6; ilay++) nClusters[ilay] = (Float_t)mult->GetNumberOfITSClusters(ilay);
  if (vertex) spdCorr = AliESDUtils::GetCorrSPD2(nClusters[1],vertex->GetZ());
  //
  Float_t centrality = -1;
  //
  // IMPORTANT CENTRALITY DEFINITION FOR pp
  //
  if (!fUseHBTmultiplicity) {
    // test binning to check if multiplicity-dependence is due to diffractive events being s
    if (trackCounter >= 0  && trackCounter <= 0)  centrality = 0;
    if (trackCounter >= 1  && trackCounter <= 1)  centrality = 1;
    if (trackCounter >= 2  && trackCounter <= 2)  centrality = 2;
    if (trackCounter >= 3  && trackCounter <= 3)  centrality = 3;
    if (trackCounter >= 4  && trackCounter <= 4)  centrality = 4;
    if (trackCounter >= 5  && trackCounter <= 5)  centrality = 5;
    // this was all the original bin 1 being [0..5]
    if (trackCounter >= 6  && trackCounter <= 9)  centrality = 6;
    if (trackCounter >= 10 && trackCounter <= 14) centrality = 7;
    if (trackCounter >= 15 && trackCounter <= 22) centrality = 8;
    if (trackCounter >= 23 && trackCounter <= 32) centrality = 9;
    if (trackCounter >= 33 && trackCounter <= 42) centrality = 10;
    /*
    if (trackCounter >= 43 && trackCounter <= 52) centrality = 7
    if (trackCounter >= 53 && trackCounter <= 62) centrality = 8;
    if (trackCounter >= 63 && trackCounter <= 72) centrality = 9;
    if (trackCounter >= 73 && trackCounter <= 82) centrality = 10;
    */
  } else {
    if (spdCorr >= 0  && spdCorr <=  2)  centrality  = 0;
    if (spdCorr >= 3  && spdCorr <=  5)  centrality  = 1;
    if (spdCorr >= 6  && spdCorr <=  8)  centrality  = 2;
    if (spdCorr >= 9  && spdCorr <= 11)  centrality  = 3;
    if (spdCorr >= 12 && spdCorr <= 14)  centrality  = 4;
    if (spdCorr >= 15 && spdCorr <= 16)  centrality  = 5;
    // this was all the original bin 1 being [0..16]
    if (spdCorr >= 17 && spdCorr <= 30)  centrality =  6;
    if (spdCorr >= 31 && spdCorr <= 45)  centrality =  7;
    if (spdCorr >= 46 && spdCorr <= 68)  centrality =  8;
    if (spdCorr >= 69 && spdCorr <= 97)  centrality =  9;
    if (spdCorr >= 98)                   centrality = 10;
    /*
    if (spdCorr >= 17 && spdCorr <= 30)  centrality = 2;
    if (spdCorr >= 31 && spdCorr <= 45)  centrality = 3;
    if (spdCorr >= 46 && spdCorr <= 68)  centrality = 4;
    if (spdCorr >= 69 && spdCorr <= 97)  centrality = 5;
    if (spdCorr >= 98)                   centrality = 6;
    */
  }
  //
  Int_t rootS = fESD->GetBeamEnergy() < 1000 ? 0 : 1;
  if (fESD->GetEventSpecie() == 4) { // PbPb
    rootS = 2;
    AliCentrality *esdCentrality = fESD->GetCentrality();
    centrality = esdCentrality->GetCentralityClass10("V0M") + 1; // centrality percentile determined with V0
    if (TMath::Abs(centrality - 1) < 1e-5) {
      centrality = esdCentrality->GetCentralityClass5("V0M");
    }
  }

  if (fIspA) {
    AliCentrality *esdCentrality = fESD->GetCentrality();
    Float_t pApercentile = esdCentrality->GetCentralityPercentile(fCentEst.Data()); // centrality percentile determined with V0M
    if (pApercentile >=  0. && pApercentile <  5.) centrality = -1; 
    if (pApercentile >=  5. && pApercentile < 10.) centrality = 0; 
    if (pApercentile >= 10. && pApercentile < 20.) centrality = 1;
    if (pApercentile >= 20. && pApercentile < 30.) centrality = 2;
    if (pApercentile >= 30. && pApercentile < 40.) centrality = 3;
    if (pApercentile >= 40. && pApercentile < 50.) centrality = 4;
    if (pApercentile >= 50. && pApercentile < 60.) centrality = 5; 
    if (pApercentile >= 60. && pApercentile < 70.) centrality = 6;
    if (pApercentile >= 70. && pApercentile < 80.) centrality = 7;
    if (pApercentile >= 80. && pApercentile < 90.) centrality = 8;
    if (pApercentile >= 90. && pApercentile <= 100.) centrality = 9;

    /*
    cout << "*****************ispA switch works***************************" << endl;
    cout << "centrality estimator  is:  " << fCentEst.Data() << endl; 
    cout << "centrality percentile is:  " << pApercentile << endl;
    cout << "*************************************************************" << endl;
    */
  }



  
  Int_t nContributors = 0;
  if (fESD->GetPrimaryVertexTPC()) nContributors = fESD->GetPrimaryVertexTPC()->GetNContributors();
  //
  
  Int_t processtype = 0;
  Int_t processCode = 0;
  //
  // important change: fill generated only after vertex cut in case of heavy-ions
  //

  
  if (!vertex || !isVertexOk) {
    fHistMult->Fill(-1, processCode);
    PostData(1, fListHist);
    return;
  } else {
    if (TMath::Abs(vertex->GetZv()) > 10) {
      fHistMult->Fill(-1, processCode);
      PostData(1, fListHist);
      return;
    }
  }
  //
  
  

  if (fMCtrue) {
    //
    //
    //
    /*
    AliHeader * header = mcEvent->Header();
    processtype = GetPythiaEventProcessType(header);
    // non diffractive
    if (processtype !=92 && processtype !=93 && processtype != 94) processCode = 1;
    // single diffractive
    if ((processtype == 92 || processtype == 93)) processCode = 2;
    // double diffractive
    if (processtype == 94) processCode = 3;
    //
    */

    for(Int_t i = 0; i < stack->GetNtrack(); i++) {
      TParticle * trackMC = stack->Particle(i);
      Int_t pdg = trackMC->GetPdgCode();
      //
      Double_t xv = trackMC->Vx();
      Double_t yv = trackMC->Vy();
      Double_t zv = trackMC->Vz();
      Double_t dxy = 0;
      dxy = TMath::Sqrt(xv*xv + yv*yv); // so stupid to avoid warnings
      Double_t dz = 0;
      dz = TMath::Abs(zv); // so stupid to avoid warnings
      //
      // vertex cut - selection of primaries
      //
      //if (dxy > 3 || dz > 10) continue; // fixed cut at 3cm in r
      //
      if (!stack->IsPhysicalPrimary(i)) continue;
      //
      // fill MC histograms here...
      // 
      Double_t rap = trackMC->Y();
      if (fRapCMS) rap = rap + 0.465;
      Double_t pT  = trackMC->Pt();
      Int_t sign = pdg < 0 ? -1 : 1; // only works for charged pi,K,p !!
//      Double_t transMass = TMath::Sqrt(trackMC->Pt()*trackMC->Pt() + trackMC->GetMass()*trackMC->GetMass()) - trackMC->GetMass();
      //
      Int_t iPart = -1;
      if (TMath::Abs(pdg) == 211)  iPart = 0; // select Pi+/Pi- only
      if (TMath::Abs(pdg) == 321)  iPart = 1; // select K+/K- only
      if (TMath::Abs(pdg) == 2212) iPart = 2; // select p+/p- only
      if (iPart == -1) continue;
      //

      if (!fSmallTHnSparse){
	Double_t vecHistMC[10] = {iPart, centrality,  pT, sign, rap, 0, 1, 0, dxy, 0};
	if (!fOnlyQA) fHistMCparticles->Fill(vecHistMC);
      }
      else{
	if (rap>fRapidityCutLow && rap<fRapidityCutHigh){
	  Double_t vecHistMC[8] = {iPart, centrality,  pT, sign, 1, 0, dxy, 0};
	  if (!fOnlyQA) fHistMCparticles->Fill(vecHistMC);
	}
      }
    }
  }
  //
  if (!isSelected && !fOnlyQA) {
    PostData(1, fListHist);
    return;
  }
  //
  if (!vertex || !isVertexOk) {
    fHistMult->Fill(-1, processCode);
    PostData(1, fListHist);
    return;
  } else {
    if (TMath::Abs(vertex->GetZv()) > 10) {
      fHistMult->Fill(-1, processCode);
      PostData(1, fListHist);
      return;
    }
  }
  //
  // count events after physics selection and after vertex selection
  //
  //cout << "MULTIPLICITY " << trackCounter << " " << fESD->GetEventNumberInFile() <<endl;
  fHistMult->Fill(trackCounter, processCode);
  fHistCentrality->Fill(centrality);


  //***************************************************
  // track loop
  //***************************************************
  //const Float_t kNsigmaCut = 3;
  //const Float_t k2sigmaCorr = 1/(0.5*(TMath::Erf(kNsigmaCut/sqrt(2))-TMath::Erf(-kNsigmaCut/sqrt(2))))/*1/0.9545*/;
  //
  Float_t dca[2], cov[3]; // dca_xy, dca_z, sigma_xy, sigma_xy_z, sigma_z for the vertex cut
  //
  for (Int_t i=0;i<fESD->GetNumberOfTracks();++i) {
	
	AliESDtrack *track = 0;
	AliESDtrack *trackForTOF = 0; //normal track for all TOF information needed when using tpconly-tracks
	
	//normal tracks, if tpconly flag is set, use tpconlytracks
	if (!fUseTPConlyTracks){
    	track =fESD->GetTrack(i); 
	}
	else {
   		track = fESDtrackCuts->GetTPCOnlyTrack(fESD,i);
		if (!track) continue;
		trackForTOF = fESD->GetTrack(i);
	}
    //
    if (!track->GetInnerParam()) {
		if (fUseTPConlyTracks) {delete track; track = 0;} //need to delete tpconlytrack
		continue;
	}
    Double_t ptot = track->GetInnerParam()->GetP(); // momentum for dEdx determination
    Double_t pT = track->Pt();
    track->GetImpactParameters(dca, cov);
    //
    //
    // cut for dead regions in the detector
    // if (track->Eta() > 0.1 && (track->Eta() < 0.2 && track->Phi() > 0.1 && track->Phi() < 0.1) continue;
    //
    // 2.a) apply some standard track cuts according to general recommendations
    //
    if (!fESDtrackCuts->AcceptTrack(track)) {
		if (fUseTPConlyTracks) {delete track; track = 0;} //need to delete tpconlytrack
		continue;
	}

	UInt_t status = 0;
    if (!fUseTPConlyTracks) status = track->GetStatus();
	else status = trackForTOF->GetStatus();
    Bool_t hasTOFout  = status&AliESDtrack::kTOFout; 
    Bool_t hasTOFtime = status&AliESDtrack::kTIME;
    Bool_t hasTOFpid  = status&AliESDtrack::kTOFpid;
    Bool_t hasTOFmismatch  = status&AliESDtrack::kTOFmismatch;
    Bool_t hasTOF     = kFALSE;
    if (hasTOFout && hasTOFtime && hasTOFpid) hasTOF = kTRUE;


    //check if TRD contributed to tracking and throw track away if  fTRDinReject flag is set
    Bool_t hasTRDin = status&AliESDtrack::kTRDin; 
    if (hasTRDin && fTRDinReject) {
      //hasTOF = kFALSE;
      if (fUseTPConlyTracks) {delete track; track = 0;} //need to delete tpconlytrack
      continue;
    }


    //check TOF window
    Float_t dxTOF = track->GetTOFsignalDx();
    Float_t dzTOF = track->GetTOFsignalDz();

    if (hasTOF) fHistTOFwindow->Fill(dxTOF,dzTOF);

    //******************************************
    //*******NEEDS PROPER CUT SETTER************
    //******************************************
    //cut on TOF window here
    if (TMath::Abs(dxTOF) > fTOFwindow || TMath::Abs(dzTOF) > fTOFwindow) hasTOF = kFALSE;

    

    Float_t length = 0.;
    if (!fUseTPConlyTracks) length = track->GetIntegratedLength(); 
    else length = trackForTOF->GetIntegratedLength();

    if (length < 350.) hasTOF = kFALSE;

    //
    // calculate rapidities and kinematics
    // 
    //
    Double_t pvec[3];
    track->GetPxPyPz(pvec);
    Double_t energyPion = TMath::Sqrt(track->GetP()*track->GetP() + AliPID::ParticleMass(AliPID::kPion)*AliPID::ParticleMass(AliPID::kPion));
    Double_t energyKaon = TMath::Sqrt(track->GetP()*track->GetP() + AliPID::ParticleMass(AliPID::kKaon)*AliPID::ParticleMass(AliPID::kKaon));
    Double_t energyProton = TMath::Sqrt(track->GetP()*track->GetP() + AliPID::ParticleMass(AliPID::kProton)*AliPID::ParticleMass(AliPID::kProton));
    Double_t energyDeuteron = TMath::Sqrt(track->GetP()*track->GetP() + 4*AliPID::ParticleMass(AliPID::kProton)*AliPID::ParticleMass(AliPID::kProton));
    //
    Double_t rapPion = 0.5*TMath::Log((energyPion + pvec[2])/(energyPion - pvec[2]));
    Double_t rapKaon = 0.5*TMath::Log((energyKaon + pvec[2])/(energyKaon - pvec[2]));
    Double_t rapProton = 0.5*TMath::Log((energyProton + pvec[2])/(energyProton - pvec[2]));
    Double_t rapDeuteron = 0.5*TMath::Log((energyDeuteron + pvec[2])/(energyDeuteron - pvec[2]));

    if (fRapCMS) {
      rapPion = rapPion + 0.465;
      rapKaon = rapKaon + 0.465;
      rapProton = rapProton + 0.465;
      rapDeuteron = rapDeuteron + 0.465;
    }


    //
//    Double_t transMassPion = TMath::Sqrt(track->Pt()*track->Pt() + AliPID::ParticleMass(AliPID::kPion)*AliPID::ParticleMass(AliPID::kPion))      -  AliPID::ParticleMass(AliPID::kPion);
//    Double_t transMassKaon = TMath::Sqrt(track->Pt()*track->Pt() + AliPID::ParticleMass(AliPID::kKaon)*AliPID::ParticleMass(AliPID::kKaon))     -  AliPID::ParticleMass(AliPID::kKaon);
 //   Double_t transMassProton = TMath::Sqrt(track->Pt()*track->Pt() + AliPID::ParticleMass(AliPID::kProton)*AliPID::ParticleMass(AliPID::kProton))    -  AliPID::ParticleMass(AliPID::kProton);
//    Double_t transMassDeuteron = TMath::Sqrt(track->Pt()*track->Pt() + 4*AliPID::ParticleMass(AliPID::kProton)*AliPID::ParticleMass(AliPID::kProton))    -  2*AliPID::ParticleMass(AliPID::kProton);
    //
    // 3. make the PID
    //
    Double_t sign = track->GetSign();   
    Double_t tpcSignal = track->GetTPCsignal();
    //
    //
    // 3.a. calculate expected signals in nsigma
    //
    //  
    // (0.) assumed particle: 0. pion, 1. kaon, 2. proton, 3. deuteron
    // (1.) multiplicity or centrality -- number of accepted ESD tracks per events (deprecated), but now classes from 1 to 10, 0: Min. Bias
    // (2.) pT
    // (3.) sign
    // (4.) rapidity --> filled 4x
    // (5.)  pull TPC dEx --> filled 4x
    // (6.) has valid TOF pid signal
    // (7.) nsigma TOF --> filled 4x
    // (8..) dca_xy
    // (9.) CODE -- only MC 0-generated, 1-true rec. primaries, 2-misident, 3-second weak, 4-second material
    //
//    Double_t transMass[4] = {transMassPion,transMassKaon,transMassProton,transMassDeuteron};
    Double_t rap[4] = {rapPion,rapKaon,rapProton,rapDeuteron};
    Double_t pullsTPC[4] = {fESDpid->NumberOfSigmasTPC(track,AliPID::kPion),
			    fESDpid->NumberOfSigmasTPC(track,AliPID::kKaon),
			    fESDpid->NumberOfSigmasTPC(track,AliPID::kProton),
			    0}; // ASK FOR PUTTING THE DEUTERON TO AliPID !!!!!!!!!!!!!!


    Float_t time0 = fESDpid->GetTOFResponse().GetStartTime(track->P());
    //Float_t time0 = fESDpid->GetTOFResponse().GetTimeZero(); //old way of getting time0
    //fESDpid->GetTOFResponse().SetTimeResolution(130.);
    Double_t pullsTOF[4] ={0.,0.,0.,0.};
    if (!fUseTPConlyTracks) {
      pullsTOF[0] = fESDpid->NumberOfSigmasTOF(track,AliPID::kPion, time0);
      pullsTOF[1] = fESDpid->NumberOfSigmasTOF(track,AliPID::kKaon, time0);
      pullsTOF[2] = fESDpid->NumberOfSigmasTOF(track,AliPID::kProton, time0);
      pullsTOF[3] = 0; // ASK FOR PUTTING THE DEUTERON TO AliPID !!!!!!!!!!!!!!;
    }
    else {
      pullsTOF[0] = fESDpid->NumberOfSigmasTOF(trackForTOF,AliPID::kPion, time0);
      pullsTOF[1] = fESDpid->NumberOfSigmasTOF(trackForTOF,AliPID::kKaon, time0);
      pullsTOF[2] = fESDpid->NumberOfSigmasTOF(trackForTOF,AliPID::kProton, time0);
      pullsTOF[3] = 0; // ASK FOR PUTTING THE DEUTERON TO AliPID !!!!!!!!!!!!!!;
    }

    //
//    Double_t tpcQA[4] = {fESDpid->NumberOfSigmasTPC(track,AliPID::kElectron),
//			 fESDpid->NumberOfSigmasTPC(track,AliPID::kPion),
//			 fESDpid->NumberOfSigmasTPC(track,AliPID::kKaon),
//			 fESDpid->NumberOfSigmasTPC(track,AliPID::kProton)};

    Double_t tofQA[4] = {0.,0.,0.,0.}; 
    if (!fUseTPConlyTracks) {
      tofQA[0] = fESDpid->NumberOfSigmasTOF(track,AliPID::kElectron, time0);
      tofQA[0] = fESDpid->NumberOfSigmasTOF(track,AliPID::kPion, time0);
      tofQA[0] = fESDpid->NumberOfSigmasTOF(track,AliPID::kKaon, time0);
      tofQA[0] = fESDpid->NumberOfSigmasTOF(track,AliPID::kProton, time0);
    }
    else{
      tofQA[0] = fESDpid->NumberOfSigmasTOF(trackForTOF,AliPID::kElectron, time0);
      tofQA[0] = fESDpid->NumberOfSigmasTOF(trackForTOF,AliPID::kPion, time0);
      tofQA[0] = fESDpid->NumberOfSigmasTOF(trackForTOF,AliPID::kKaon, time0);
      tofQA[0] = fESDpid->NumberOfSigmasTOF(trackForTOF,AliPID::kProton, time0);
    }


    //save information for every particle type  // loop over assumed particle type
    for(Int_t iPart = 0; iPart < 3; iPart++) {

      if (!fSmallTHnSparse) {
	Double_t vecHistReal[9]  = {iPart,  centrality,   pT, sign,  rap[iPart], pullsTPC[iPart], hasTOF, pullsTOF[iPart], dca[0]};
	if (!fOnlyQA) fHistRealTracks->Fill(vecHistReal);
      }
      else {
	if (pullsTPC[iPart]>fTPCnSigmaCutLow && pullsTPC[iPart]<fTPCnSigmaCutHigh && rap[iPart]>fRapidityCutLow && rap[iPart]<fRapidityCutHigh) {
	  Double_t vecHistReal[7]  = {iPart,  centrality,   pT, sign, hasTOF, pullsTOF[iPart], dca[0]};
	  if (!fOnlyQA) fHistRealTracks->Fill(vecHistReal);
	}
      }


      // using MC truth for precise efficiencies...
      //
      if (fMCtrue && !fOnlyQA) {
	Int_t code = 9; // code: 0-generated, 1-true rec. primaries, 2-misident, 3-second weak, 4-second material
	Int_t assumedPdg = 0;//2212(proton); 321(Kaon); 211(pion);
	Int_t motherCode = -1;
	if (iPart == 0) assumedPdg = 211;
	if (iPart == 1) assumedPdg = 321;
	if (iPart == 2) assumedPdg = 2212;
	//
	//
	TParticle *trackMC = stack->Particle(TMath::Abs(track->GetLabel()));
	Int_t pdg = TMath::Abs(trackMC->GetPdgCode());
	//
	if (pdg != assumedPdg && stack->IsPhysicalPrimary(TMath::Abs(track->GetLabel()))) code = 2;
	if (pdg != assumedPdg && stack->IsSecondaryFromWeakDecay(TMath::Abs(track->GetLabel()))) code = 5;
	if (pdg == assumedPdg && stack->IsPhysicalPrimary(TMath::Abs(track->GetLabel()))) code = 1;
	if (pdg == assumedPdg && stack->IsSecondaryFromWeakDecay(TMath::Abs(track->GetLabel()))) {
	  code = 3;
	  if (fSaveMotherPDG){
	    TParticle *trackMother =  stack->Particle(TMath::Abs(trackMC->GetFirstMother()));
	    if (trackMother->GetPdgCode() == 310) motherCode = 6; //K0
	    if (trackMother->GetPdgCode() == 3122) motherCode = 7; //Lambda
	    if (trackMother->GetPdgCode() == 3222) motherCode = 8; //Sigma+
	  }
	}
	if (pdg == assumedPdg && stack->IsSecondaryFromMaterial(TMath::Abs(track->GetLabel()))) code = 4;
	
	//
	// muons need special treatment, because they are indistinguishable from pions
	//
	if (iPart == 0 && pdg == 13  && stack->IsPhysicalPrimary(TMath::Abs(track->GetLabel()))) code = 1;
	if (iPart == 0 && pdg == 13  && stack->IsSecondaryFromWeakDecay(TMath::Abs(track->GetLabel()))) code = 3;
	//
	// check TOF mismatch on MC basis with TOF label
	//
	Int_t tofLabel[3];
	if (!fUseTPConlyTracks) track->GetTOFLabel(tofLabel);
	else trackForTOF->GetTOFLabel(tofLabel);


	//three options:
	//0: do NOT check at all
	//1: do check
	//2: in case of decays, check if mother label matches --> if yes, assign hasTOF = kTRUE

	if (fTOFmisMatch == 0) {
	  //if (TMath::Abs(track->GetLabel()) != TMath::Abs(tofLabel[0]) || tofLabel[1] > 0) hasTOF = kFALSE;
	}
	if (fTOFmisMatch == 1) {
	  if (TMath::Abs(track->GetLabel()) != TMath::Abs(tofLabel[0]) || tofLabel[1] > 0) hasTOF = kFALSE;
	}
	if (fTOFmisMatch == 2) {
	  if (TMath::Abs(track->GetLabel()) != TMath::Abs(tofLabel[0]) || tofLabel[1] > 0) hasTOF = kFALSE;
 	  TParticle *matchedTrack = stack->Particle(TMath::Abs(tofLabel[0]));
	  if (TMath::Abs(matchedTrack->GetFirstMother()) == TMath::Abs(track->GetLabel())) hasTOF = kTRUE;
	}

	  

	  //
	// IMPORTANT BIG PROBLEM HERE THE PROBABLILITY TO HAVE A PID SIGNAL MUST BE IN !!!!!!!!!!!!
	//
	if (!fSmallTHnSparse){
	  Double_t vectorHistMC[10] = {iPart,  centrality,  pT, sign,  rap[iPart], pullsTPC[iPart], hasTOF, pullsTOF[iPart], dca[0], code};
	  if (!fOnlyQA) { 
	    fHistMCparticles->Fill(vectorHistMC);
	    if (motherCode != -1 && fSaveMotherPDG) { //if mother of weak decay is K0, lambda or sigma+ add track again with this information
	      Double_t vectorHistMCmother[10] = {iPart,  centrality,  pT, sign,  rap[iPart], pullsTPC[iPart], hasTOF, pullsTOF[iPart], dca[0], motherCode};
	      fHistMCparticles->Fill(vectorHistMCmother);
	    }
	  }
	}
	else{
	  if (pullsTPC[iPart]>fTPCnSigmaCutLow && pullsTPC[iPart]<fTPCnSigmaCutHigh && rap[iPart]>fRapidityCutLow && rap[iPart]<fRapidityCutHigh) {
	    //                              0,           1,   2,    3,           4,               5,      6,               7,      8,   9
	    Double_t vectorHistMC[8] = {iPart,  centrality,  pT, sign, hasTOF, pullsTOF[iPart], dca[0], code};
	    if (!fOnlyQA) { 
	      fHistMCparticles->Fill(vectorHistMC);
	      if (motherCode != -1 && fSaveMotherPDG) { //if mother of weak decay is K0, lambda or sigma+ add track again with this information
		Double_t vectorHistMCmother[8] = {iPart,  centrality,  pT, sign, hasTOF, pullsTOF[iPart], dca[0], motherCode};
		fHistMCparticles->Fill(vectorHistMCmother);
	      }
	    }
	  }
	}
      }
      //
      //
      Int_t tpcShared = track->GetTPCnclsS();
      if (TMath::Abs(track->Eta()) < 0.8 && iPart == 0 && tpcShared < 4) fHistPidQA->Fill(ptot,tpcSignal,sign);
    } // end loop over assumed particle type

	  //need to delete tpconlytrack
	  if (fUseTPConlyTracks){
 		delete track; 
  		track = 0;     
  	  }

  } // end of track loop
  
  // Post output data  
  PostData(1, fListHist);
  
}      


//________________________________________________________________________
void AliAnalysisTPCTOFpA::Terminate(Option_t *) 
{
  // Draw result to the screen
  // Called once at the end of the query
  Printf("*** CONSTRUCTOR CALLED ****");

}


//________________________________________________________________________
Bool_t AliAnalysisTPCTOFpA::SelectOnImpPar(AliESDtrack* t) {
  //
  // cut on transverse impact parameter // DEPRECATED
  //
  Float_t d0z0[2],covd0z0[3];
  t->GetImpactParameters(d0z0,covd0z0);
  Float_t sigma= 0.0050+0.0060/TMath::Power(t->Pt(),0.9);
  Float_t d0max = 7.*sigma;
  //
  Float_t sigmaZ = 0.0146+0.0070/TMath::Power(t->Pt(),1.114758);
  if (t->Pt() > 1) sigmaZ = 0.0216;
  Float_t d0maxZ = 5.*sigmaZ;
  //
  if(TMath::Abs(d0z0[0]) < d0max && TMath::Abs(d0z0[1]) < d0maxZ) return kTRUE;
  return kFALSE;
}


//________________________________________________________________________
void AliAnalysisTPCTOFpA::BinLogAxis(const TH1 *h) {
  //
  // Method for the correct logarithmic binning of histograms
  //
  TAxis *axis = h->GetXaxis();
  int bins = axis->GetNbins();

  Double_t from = axis->GetXmin();
  Double_t to = axis->GetXmax();
  Double_t *newBins = new Double_t[bins + 1];
   
  newBins[0] = from;
  Double_t factor = pow(to/from, 1./bins);
  
  for (int i = 1; i <= bins; i++) {
   newBins[i] = factor * newBins[i-1];
  }
  axis->Set(bins, newBins);
  delete [] newBins;
  
}


//________________________________________________________________________
Int_t AliAnalysisTPCTOFpA::GetPythiaEventProcessType(const AliHeader* aHeader, const Bool_t adebug) const {
  //
  // get the process type of the event.
  //

  // can only read pythia headers, either directly or from cocktalil header
  AliGenPythiaEventHeader* pythiaGenHeader = dynamic_cast<AliGenPythiaEventHeader*>(aHeader->GenEventHeader());

  if (!pythiaGenHeader) {

    AliGenCocktailEventHeader* genCocktailHeader = dynamic_cast<AliGenCocktailEventHeader*>(aHeader->GenEventHeader());
    if (!genCocktailHeader) {
      //printf("AliAnalysisTPCTOFpA::GetProcessType : Unknown header type (not Pythia or Cocktail). \n");
      return -1;
    }

    TList* headerList = genCocktailHeader->GetHeaders();
    if (!headerList) {
      return -1;
    }

    for (Int_t i=0; i<headerList->GetEntries(); i++) {
      pythiaGenHeader = dynamic_cast<AliGenPythiaEventHeader*>(headerList->At(i));
      if (pythiaGenHeader)
        break;
    }

    if (!pythiaGenHeader) {
      //printf("AliAnalysisTPCTOFpA::GetProcessType : Could not find Pythia header. \n");
      return -1;
    }
  }

  if (adebug) {
    //printf("AliAnalysisTPCTOFpA::GetProcessType : Pythia process type found: %d \n",pythiaGenHeader->ProcessType());
  }

  return pythiaGenHeader->ProcessType();
}
