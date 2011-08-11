/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, proviyaded that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purapose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

///////////////////////////////////////////////////////////////////////////
//                                                                       //
//                                                                       //
// Analysis for identified charged hadron spectra.                       //
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

#include "AliAnalysisCombinedHadronSpectra.h"


ClassImp(AliAnalysisCombinedHadronSpectra)

//________________________________________________________________________
AliAnalysisCombinedHadronSpectra::AliAnalysisCombinedHadronSpectra() 
  : AliAnalysisTaskSE("TaskChargedHadron"), fESD(0), fListHist(0), fESDtrackCuts(0),fESDTrackCutsMult(0),fESDpid(0),
    fMCtrue(0),
    fOnlyQA(0),
    fUseHBTmultiplicity(0),
    fAlephParameters(),
    fHistRealTracks(0),
    fHistMCparticles(0),
    fHistPidQA(0),
    fHistMult(0),
    fHistCentrality(0)
{
  // default Constructor
  /* fast compilation test
     gSystem->Load("libANALYSIS");
     gSystem->Load("libANALYSISalice");
     .L /d/alice09/akalweit/train/trunk/akalweit_hadronSpectra/AliAnalysisCombinedHadronSpectra.cxx++
   */
}


//________________________________________________________________________
AliAnalysisCombinedHadronSpectra::AliAnalysisCombinedHadronSpectra(const char *name) 
  : AliAnalysisTaskSE(name), fESD(0), fListHist(0), fESDtrackCuts(0),fESDTrackCutsMult(0),fESDpid(0),
    fMCtrue(0),
    fOnlyQA(0),
    fUseHBTmultiplicity(0),
    fAlephParameters(),
    fHistRealTracks(0),
    fHistMCparticles(0),
    fHistPidQA(0),
    fHistMult(0),
    fHistCentrality(0)
{
  //
  // standard constructur which should be used
  //
  Printf("*** CONSTRUCTOR CALLED ****");
  //
  fMCtrue = kTRUE; 
  fOnlyQA = kFALSE;
  fUseHBTmultiplicity = kTRUE;
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
  Initialize();
  // Output slot #0 writes into a TList container
  DefineOutput(1, TList::Class());

}


//________________________________________________________________________
void AliAnalysisCombinedHadronSpectra::Initialize()
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
  fESDtrackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2010(kFALSE);
  fESDtrackCuts->SetMaxDCAToVertexXY(3);
  fESDtrackCuts->SetMaxDCAToVertexZ(3);
  //
  //
  //
  //
  fESDTrackCutsMult = AliESDtrackCuts::GetStandardITSTPCTrackCuts2010(kTRUE);
  fESDTrackCutsMult->SetEtaRange(-1.2,+1.2);
  fESDTrackCutsMult->SetPtRange(0.15,1e10);

}


//________________________________________________________________________
void AliAnalysisCombinedHadronSpectra::UserCreateOutputObjects() 
{
  // Create histograms
  // Called once
  fListHist = new TList();
  fListHist->SetOwner(kTRUE);
  //
  const Int_t kPtBins = 35;
  const Int_t kMultBins = 11;
  const Int_t kDcaBins = 76;
  // sort pT-bins ..
  Double_t binsPt[kPtBins+1] = {0., 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0};
  Double_t binsDca[kDcaBins+1] = {-3,-2.85,-2.7,-2.55,-2.4,-2.25,-2.1,-1.95,-1.8,-1.65,-1.5,-1.35,-1.2,-1.05,-0.9,-0.75,-0.6,-0.45,-0.3,-0.285,-0.27,-0.255,-0.24,-0.225,-0.21,-0.195,-0.18,-0.165,-0.15,-0.135,-0.12,-0.105,-0.09,-0.075,-0.06,-0.045,-0.03,-0.015,0,0.015,0.03,0.045,0.06,0.075,0.09,0.105,0.12,0.135,0.15,0.165,0.18,0.195,0.21,0.225,0.24,0.255,0.27,0.285,0.3,0.45,0.6,0.75,0.9,1.05,1.2,1.35,1.5,1.65,1.8,1.95,2.1,2.25,2.4,2.55,2.7,2.85,3};
  //
  // create the histograms with all necessary information --> it is filled 4x for each particle assumption
  //
  // (0.) assumed particle: 0. pion, 1. kaon, 2. proton, 3. deuteron
  // (1.) trigger event class: MB, high-mult, etc. / we could also put here a trigger pT-particle
  // (2.) multiplicity or centrality -- number of accepted ESD tracks per events (deprecated), but now classes from 1 to 10, 0: Min. Bias
  // (3.) sqrt(s) -- center of mass energy, 0. 900 GeV or 1. 7 TeV, 2, PbPb
  // (4.) pT
  // (5.) sign
  // (6.) mT - m0 (transverse kinetic energy) --> filled 4x // SHADOWED FOR THE MOMENT !!!!
  // (7.) rapidity --> filled 4x
  // (8)  pull TPC dEx --> filled 4x
  // (9.) has valid TOF pid signal
  // (10.) nsigma TOF --> filled 4x
  // (11.) dca_xy
  // (12.) dca_z
  // (13.) CODE -- only MC 0-generated, 1-true rec. primaries, 2-misident, 3-second weak, 4-second material 5-unknown
  //
  //                              0,    1,         2,    3,        4,  5,   6,    7,   8,    9,  10,       11, 12
  Int_t    binsHistReal[13] = {   3,    1, kMultBins,    1,  kPtBins,  2,   1,   10,  50,    2,  50, kDcaBins,  5};
  Double_t xminHistReal[13] = {-0.5,    2,      -0.5, -0.5,        0, -2,   0, -0.5,  -5,- 0.5,  -5,       -3, -2};
  Double_t xmaxHistReal[13] = { 2.5,   10,      10.5,  2.5,        3,  2,   3,  0.5,   5,  1.5,   5,        3,  2};
  fHistRealTracks = new THnSparseL("fHistRealTracks","real tracks",13,binsHistReal,xminHistReal,xmaxHistReal);
  //
  fHistRealTracks->GetAxis(4)->Set(kPtBins, binsPt);
  //fHistRealTracks->GetAxis(2)->Set(kMultBins, multBins);
  fHistRealTracks->GetAxis(11)->Set(kDcaBins, binsDca);
  fListHist->Add(fHistRealTracks);
  //
  //                      0.ptot,1.tpcSig,2.hasTOF, 3. assumed part., 4. nclDedx, 5. nSigmaTPC (4x), 6. nSigmaTOF (4x), 7. centrality
  Int_t    binsHistQA[8] = {100,     600,       2,                4,  50,  50, 50,   20};
  Double_t xminHistQA[8] = {0.1,      30,    -0.5,             -0.5,  60,  -5, -5,   0.};
  Double_t xmaxHistQA[8] = {  4,     500,     1.5,              3.5, 160,   5,  5, 4000};
  fHistPidQA = new THnSparseL("fHistPidQA","PID QA",8,binsHistQA,xminHistQA,xmaxHistQA);
  BinLogAxis(fHistPidQA, 0);
  fListHist->Add(fHistPidQA);
  //                            0,    1,         2,    3,        4,  5,   6,    7,   8,    9,  10,       11, 12,   13
  Int_t    binsHistMC[14] = {   3,    1, kMultBins,    1,  kPtBins,  2,   1,   10,  50,    2,  50, kDcaBins,  5,    6};
  Double_t xminHistMC[14] = {-0.5,    2,      -0.5, -0.5,        0, -2,   0, -0.5,  -5,- 0.5,  -5,       -3, -2, -0.5};
  Double_t xmaxHistMC[14] = { 2.5,   10,      10.5,  2.5,        3,  2,   3,  0.5,   5,  1.5,   5,        3,  2,  5.5};
  fHistMCparticles = new THnSparseL("fHistMCparticles","MC histogram",14,binsHistMC,xminHistMC,xmaxHistMC);
  fHistMCparticles->GetAxis(4)->Set(kPtBins, binsPt);
  fHistMCparticles->GetAxis(11)->Set(kDcaBins, binsDca);
  //fHistMCparticles->GetAxis(2)->Set(kMultBins, multBins);
  //fHistMCparticles->GetAxis(6)->Set(kPtBins, binsPt);
  fListHist->Add(fHistMCparticles);
  //
  fHistMult = new TH2D("fHistMult", "control histogram to count number of events", 502, -2.5, 499.5,4,-0.5,3.5);
  fHistCentrality = new TH1D("fHistCentrality", "control histogram to count number of events", 22, -1.5, 20.5);
  fListHist->Add(fHistMult);
  fListHist->Add(fHistCentrality);
  
}

//________________________________________________________________________
void AliAnalysisCombinedHadronSpectra::UserExec(Option_t *) 
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
  const AliESDVertex *vertex = fESD->GetPrimaryVertexTracks();
  if(vertex->GetNContributors()<1) {
    // SPD vertex
    vertex = fESD->GetPrimaryVertexSPD();
    if(vertex->GetNContributors()<1) vertex = 0x0;
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
    if (centrality == 1) {
      centrality = esdCentrality->GetCentralityClass5("V0M");
    }
  }
  Int_t nContributors = 0;
  if (fESD->GetPrimaryVertexTPC()) nContributors = fESD->GetPrimaryVertexTPC()->GetNContributors();
  //
  Int_t processtype = 0;
  Int_t processCode = 0;
  if (fMCtrue) {
    //
    //
    //
    AliHeader * header = mcEvent->Header();
    processtype = GetPythiaEventProcessType(header);
    //
    // HARD DEBUG OF MULTIPLICITY DEPENDENT EFFICIENCY -> CALCULATE EFF. ONLY FOR NON-DIFFRACTIVE EVENTS -- PLEASE REMOVE AS SOON AS POSSIBLE
    //
    // DEBUG DEBUG  DEBUG DEBUG  DEBUG DEBUG  DEBUG DEBUG  DEBUG DEBUG  DEBUG DEBUG  DEBUG DEBUG  DEBUG DEBUG  DEBUG DEBUG
    if (processtype == 92 || processtype ==93 || processtype == 94) return;
    // DEBUG DEBUG  DEBUG DEBUG  DEBUG DEBUG  DEBUG DEBUG  DEBUG DEBUG  DEBUG DEBUG  DEBUG DEBUG  DEBUG DEBUG  DEBUG DEBUG 
    // non diffractive
    if (processtype !=92 && processtype !=93 && processtype != 94) processCode = 1;
    // single diffractive
    if ((processtype == 92 || processtype == 93)) processCode = 2;
    // double diffractive
    if (processtype == 94) processCode = 3;
    //
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
      Double_t pT  = trackMC->Pt();
      Int_t sign = pdg < 0 ? -1 : 1; // only works for charged pi,K,p !!
      Double_t transMass = TMath::Sqrt(trackMC->Pt()*trackMC->Pt() + trackMC->GetMass()*trackMC->GetMass()) 
	- trackMC->GetMass();
      //
      Int_t iPart = -1;
      if (TMath::Abs(pdg) == 211)  iPart = 0; // select Pi+/Pi- only
      if (TMath::Abs(pdg) == 321)  iPart = 1; // select K+/K- only
      if (TMath::Abs(pdg) == 2212) iPart = 2; // select p+/p- only
      if (iPart == -1) continue;
      //
      Double_t vecHistMC[14] = {iPart, triggerPt, centrality, rootS,  pT, sign, transMass, rap, 0, 1, 0, dxy, 0, 0};
      if (!fOnlyQA) fHistMCparticles->Fill(vecHistMC);
    }
  }
  //
  if (!isSelected && !fOnlyQA) {
    PostData(1, fListHist);
    return;
  }
  //
  if (!vertex) {
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
  //
  // track loop
  //
  //const Float_t kNsigmaCut = 3;
  //const Float_t k2sigmaCorr = 1/(0.5*(TMath::Erf(kNsigmaCut/sqrt(2))-TMath::Erf(-kNsigmaCut/sqrt(2))))/*1/0.9545*/;
  //
  Float_t dca[2], cov[3]; // dca_xy, dca_z, sigma_xy, sigma_xy_z, sigma_z for the vertex cut
  //
  for (Int_t i=0;i<fESD->GetNumberOfTracks();++i) {
    AliESDtrack *track =fESD->GetTrack(i); 
    //
    if (!track->GetInnerParam()) continue;
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
    if (!fESDtrackCuts->AcceptTrack(track)) continue;
    UInt_t status = track->GetStatus();
    Bool_t hasTOFout  = status&AliESDtrack::kTOFout; 
    Bool_t hasTOFtime = status&AliESDtrack::kTIME;
    Bool_t hasTOF     = hasTOFout && hasTOFtime;
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
    //
    Double_t transMassPion = TMath::Sqrt(track->Pt()*track->Pt() + AliPID::ParticleMass(AliPID::kPion)*AliPID::ParticleMass(AliPID::kPion)) 
      -  AliPID::ParticleMass(AliPID::kPion);
    Double_t transMassKaon = TMath::Sqrt(track->Pt()*track->Pt() + AliPID::ParticleMass(AliPID::kKaon)*AliPID::ParticleMass(AliPID::kKaon)) 
      -  AliPID::ParticleMass(AliPID::kKaon);
    Double_t transMassProton = TMath::Sqrt(track->Pt()*track->Pt() + AliPID::ParticleMass(AliPID::kProton)*AliPID::ParticleMass(AliPID::kProton)) 
      -  AliPID::ParticleMass(AliPID::kProton);
    Double_t transMassDeuteron = TMath::Sqrt(track->Pt()*track->Pt() + 4*AliPID::ParticleMass(AliPID::kProton)*AliPID::ParticleMass(AliPID::kProton))
      -  2*AliPID::ParticleMass(AliPID::kProton);
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
    //
    // fill the histogram
    // (0.) assumed particle: 0. pion, 1. kaon, 2. proton, 3. deuteron
    // (1.) trigger event class: MB, high-mult, etc.
    // (2.) multiplicity -- number of accepted ESD tracks per events
    // (3.) sqrt(s) -- center of mass energy, 0. 900 GeV or 1. 7 TeV
    // (4.) pT
    // (5.) sign
    // (6.) mT - m0 (transverse kinetic energy) --> filled 4x
    // (7.) rapidity --> filled 4x
    // (8) pull TPC dEx --> filled 4x
    // (9.) has valid TOF pid signal
    // (10.) nsigma TOF --> filled 4x
    // (11.) dca_xy
    // (12.) dca_z
    //
    Double_t transMass[4] = {transMassPion,transMassKaon,transMassProton,transMassDeuteron};
    Double_t rap[4] = {rapPion,rapKaon,rapProton,rapDeuteron};
    Double_t pullsTPC[4] = {fESDpid->NumberOfSigmasTPC(track,AliPID::kPion),
			    fESDpid->NumberOfSigmasTPC(track,AliPID::kKaon),
			    fESDpid->NumberOfSigmasTPC(track,AliPID::kProton),
			    0}; // ASK FOR PUTTING THE DEUTERON TO AliPID !!!!!!!!!!!!!!
    Float_t time0 = fESDpid->GetTOFResponse().GetTimeZero();
    //fESDpid->GetTOFResponse().SetTimeResolution(130.);
    Double_t pullsTOF[4] =  {fESDpid->NumberOfSigmasTOF(track,AliPID::kPion, time0),
			     fESDpid->NumberOfSigmasTOF(track,AliPID::kKaon, time0),
			     fESDpid->NumberOfSigmasTOF(track,AliPID::kProton, time0),
			     0}; // ASK FOR PUTTING THE DEUTERON TO AliPID !!!!!!!!!!!!!!;
    //
    Double_t tpcQA[4] = {fESDpid->NumberOfSigmasTPC(track,AliPID::kElectron),
			 fESDpid->NumberOfSigmasTPC(track,AliPID::kPion),
			 fESDpid->NumberOfSigmasTPC(track,AliPID::kKaon),
			 fESDpid->NumberOfSigmasTPC(track,AliPID::kProton)};
    Double_t tofQA[4] = {fESDpid->NumberOfSigmasTOF(track,AliPID::kElectron, time0),
			 fESDpid->NumberOfSigmasTOF(track,AliPID::kPion, time0),
			 fESDpid->NumberOfSigmasTOF(track,AliPID::kKaon, time0),
			 fESDpid->NumberOfSigmasTOF(track,AliPID::kProton, time0)};
    //
    for(Int_t iPart = 0; iPart < 3; iPart++) { // loop over assumed particle type
      //                              0,         1,          2,     3,   4,    5,                6,          7,               8,      9,              10,     11,     12
      Double_t vecHistReal[13] = {iPart, triggerPt, centrality, rootS,  pT, sign, transMass[iPart], rap[iPart], pullsTPC[iPart], hasTOF, pullsTOF[iPart], dca[0], dca[1]};
      if (!fOnlyQA) fHistRealTracks->Fill(vecHistReal);
      //
      // using MC truth for precise efficiencies...
      //
      if (fMCtrue && !fOnlyQA) {
	Int_t code = 5; // code: 0-generated, 1-true rec. primaries, 2-misident, 3-second weak, 4-second material
	Int_t assumedPdg = 0;//2212(proton); 321(Kaon); 211(pion);
	if (iPart == 0) assumedPdg = 211;
	if (iPart == 1) assumedPdg = 321;
	if (iPart == 2) assumedPdg = 2212;
	//
	//
	TParticle *trackMC = stack->Particle(TMath::Abs(track->GetLabel()));
	Int_t pdg = TMath::Abs(trackMC->GetPdgCode());
	//
	if (pdg != assumedPdg) code = 2;
	if (pdg == assumedPdg && stack->IsPhysicalPrimary(TMath::Abs(track->GetLabel()))) code = 1;
	if (pdg == assumedPdg && stack->IsSecondaryFromWeakDecay(TMath::Abs(track->GetLabel()))) code = 3;
	if (pdg == assumedPdg && stack->IsSecondaryFromMaterial(TMath::Abs(track->GetLabel()))) code = 4;
	/*
	if (pdg == assumedPdg && !stack->IsPhysicalPrimary(TMath::Abs(track->GetLabel())) && trackMC->GetFirstMother() > 0) {
	  TParticle *trackMother = stack->Particle(trackMC->GetFirstMother());
	  Int_t mPdg = TMath::Abs(trackMother->GetPdgCode());
	  if (mPdg==3122||mPdg==3222||mPdg==3212||mPdg==3112||mPdg==3322||mPdg==3312||mPdg==3332||mPdg==130||mPdg==310) code = 3;
	  if (mPdg!=3122&&mPdg!=3222&&mPdg!=3212&&mPdg!=3112&&mPdg!=3322&&mPdg!=3312&&mPdg!=3332&&mPdg!=130&&mPdg!=310) code = 4;
	}
	*/
	//                               0,         1,          2,     3,   4,    5,               6,           7, 8,      9,10,     11,     12,   13
	//
	// IMPORTANT BIG PROBLEM HERE THE PROBABLILITY TO HAVE A PID SIGNAL MUST BE IN !!!!!!!!!!!!
	//
	Double_t vectorHistMC[14] = {iPart, triggerPt, centrality, rootS,  pT, sign, transMass[iPart], rap[iPart], pullsTPC[iPart], hasTOF, pullsTOF[iPart], dca[0], dca[1], code};
	if (!fOnlyQA) fHistMCparticles->Fill(vectorHistMC);
      }
      //
      //
      //                      0.ptot,1.tpcSig, 2.hasTOF, 3. assumed part., 4, nclDedx, 5. nSigmaTPC (4x), 6. nSigmaTOF (4x), 7. evt. mult.
      Double_t vecHistQA[8] = {ptot, tpcSignal, hasTOF, iPart, track->GetTPCsignalN(), tpcQA[iPart], tofQA[iPart], nContributors};
      if (TMath::Abs(track->Eta()) < 0.8 && fHistPidQA->GetEntries() < 1e6 && fOnlyQA) fHistPidQA->Fill(vecHistQA);
    } // end loop over assumed particle type
    
    
  } // end of track loop
  
  // Post output data  
  PostData(1, fListHist);
  
}      


//________________________________________________________________________
void AliAnalysisCombinedHadronSpectra::Terminate(Option_t *) 
{
  // Draw result to the screen
  // Called once at the end of the query
  Printf("*** CONSTRUCTOR CALLED ****");

}


//________________________________________________________________________
Bool_t AliAnalysisCombinedHadronSpectra::SelectOnImpPar(AliESDtrack* t) {
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
void AliAnalysisCombinedHadronSpectra::BinLogAxis(const THnSparse *h, Int_t axisNumber) {
  //
  // Method for the correct logarithmic binning of histograms
  //
  TAxis *axis = h->GetAxis(axisNumber);
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
Int_t AliAnalysisCombinedHadronSpectra::GetPythiaEventProcessType(const AliHeader* aHeader, const Bool_t adebug) const {
  //
  // get the process type of the event.
  //

  // can only read pythia headers, either directly or from cocktalil header
  AliGenPythiaEventHeader* pythiaGenHeader = dynamic_cast<AliGenPythiaEventHeader*>(aHeader->GenEventHeader());

  if (!pythiaGenHeader) {

    AliGenCocktailEventHeader* genCocktailHeader = dynamic_cast<AliGenCocktailEventHeader*>(aHeader->GenEventHeader());
    if (!genCocktailHeader) {
      //printf("AliAnalysisCombinedHadronSpectra::GetProcessType : Unknown header type (not Pythia or Cocktail). \n");
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
      //printf("AliAnalysisCombinedHadronSpectra::GetProcessType : Could not find Pythia header. \n");
      return -1;
    }
  }

  if (adebug) {
    //printf("AliAnalysisCombinedHadronSpectra::GetProcessType : Pythia process type found: %d \n",pythiaGenHeader->ProcessType());
  }

  return pythiaGenHeader->ProcessType();
}
