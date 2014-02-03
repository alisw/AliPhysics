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
// Analysis for identified charged hadron spectra.                       //
//                                                                       //
//                                                                       //
///////////////////////////////////////////////////////////////////////////

#include "Riostream.h"
#include "TChain.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
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

#include "AliAnalysisDeuteronpA.h"

#include "AliAnalysisUtils.h"


ClassImp(AliAnalysisDeuteronpA)

//________________________________________________________________________
AliAnalysisDeuteronpA::AliAnalysisDeuteronpA() 
: AliAnalysisTaskSE("TaskChargedHadron"), fESD(0), fListHist(0), fESDtrackCuts(0),fESDTrackCutsMult(0),fESDpid(0),
  fMCtrue(0),
  fRapCMSpA(0),
  fAlephParameters(),
  fUtils(0),
  fHistRealTracks(0),
  fHistMCparticles(0),
  fHistPidQA(0),
  fHistTofQA(0),
  fHistMult(0),
  fHistCentrality(0),
  fHistMomCorr(0),
  fHistEtaPtGen(0),
  fHistVertex(0),
  fHistVertexRes(0),
  fHistVertexResTracks(0)
{
  // default Constructor
  /* fast compilation test
     gSystem->Load("libANALYSIS");
     gSystem->Load("libANALYSISalice");
     .L AliAnalysisDeuteronpA.cxx++
   */
}


//________________________________________________________________________
AliAnalysisDeuteronpA::AliAnalysisDeuteronpA(const char *name) 
  : AliAnalysisTaskSE(name), fESD(0), fListHist(0), fESDtrackCuts(0),fESDTrackCutsMult(0),fESDpid(0),
    fMCtrue(0),
    fRapCMSpA(0),
    fAlephParameters(),
    fUtils(0),
    fHistRealTracks(0),
    fHistMCparticles(0),
    fHistPidQA(0),
    fHistTofQA(0),
    fHistMult(0),
    fHistCentrality(0),
    fHistMomCorr(0),
    fHistEtaPtGen(0),
    fHistVertex(0),
    fHistVertexRes(0),
    fHistVertexResTracks(0)
{
  //
  // standard constructur which should be used
  //
  Printf("*** CONSTRUCTOR CALLED ****");
  //
  fMCtrue = kTRUE; 
  fRapCMSpA = kTRUE;
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
void AliAnalysisDeuteronpA::Initialize()
{
  //
  // updating parameters in case of changes
  //


  fESDtrackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011(kFALSE,kTRUE);
  fESDtrackCuts->SetMaxDCAToVertexXY(3);
  fESDtrackCuts->SetMaxDCAToVertexZ(2);
  fESDtrackCuts->SetEtaRange(-0.9,0.9);
  //
  //

}


//________________________________________________________________________
void AliAnalysisDeuteronpA::UserCreateOutputObjects() 
{

  //Create analysis utils for event selection and pileup rejection
  fUtils = new AliAnalysisUtils();
  fUtils->SetCutOnZVertexSPD(kFALSE);
  //
  //
  // Create histograms
  // Called once
  fListHist = new TList();
  fListHist->SetOwner(kTRUE);
  //
  const Int_t kPtBins = 28;
  const Int_t kMultBins = 11;
  const Int_t kDcaBins = 38;


  //different binning for YCMS in pA to cover detector acceptance
  Double_t kYlowBorder = -0.6;
  Double_t kYhighBorder = 0.6;
  Double_t kYBins = 12;

  if (fRapCMSpA){
    kYlowBorder = -0.1;
    kYhighBorder = 1.1;
    kYBins = 12;
  }

  //
  // sort pT-bins ..
  //
  Double_t binsPt[kPtBins+1] = {0., 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.2, 4.4, 4.6, 4.8, 5.0, 5.2, 5.4, 5.6};
  //
  // provide a reasonable dca/binning
  //
  Double_t binsDca[kDcaBins+1];
  for(Int_t i = 0; i< kDcaBins+1; i++) {
    binsDca[i] = 1.5*(1./TMath::Exp(-TMath::Abs((i - kDcaBins/2.))/(kDcaBins/2)) -1);
    if ((i - kDcaBins/2.) < 0) binsDca[i] *= -1;
    //cout <<  binsDca[i] << endl;
  }
  //
  // create the histograms with all necessary information --> it is filled 4x for each particle assumption
  //
  // (0.) assumed particle: 0. deuteron, 1. triton, 2. He-3
  // (1.) multiplicity or centrality -- number of accepted ESD tracks per events (deprecated), but now classes from 1 to 10, 0: Min. Bias
  // (2.) pT
  // (3.) sign
  // (4.) rapidity --> filled 4x
  // (5.)  pull TPC dEx --> filled 4x
  // (6.) has valid TOF pid signal
  // (7.) nsigma TOF --> filled 4x XXXXXXXXX no mass*mass
  // (8..) dca_xy
  // (9.) CODE -- only MC 0-generated, 1-true rec. primaries, 2-misident prim, 3-second weak, 4-second material, 5-misident sec
  //
  //                              0,           1,           2,  3,      4,              5,    6,    7,       8
  Int_t    binsHistReal[9] = {   3,   kMultBins,     kPtBins,  2,      kYBins      ,   50,    2,  100, kDcaBins};
  Double_t xminHistReal[9] = {-0.5,        -0.5,           0, -2,      kYlowBorder ,   -5,- 0.5, -2.5,       -3};
  Double_t xmaxHistReal[9] = { 2.5,        10.5,           3,  2,      kYhighBorder,    5,  1.5,  2.5,        3};
  fHistRealTracks = new THnSparseF("fHistRealTracks","real tracks",9,binsHistReal,xminHistReal,xmaxHistReal);
  //
  fHistRealTracks->GetAxis(2)->Set(kPtBins, binsPt);
  fHistRealTracks->GetAxis(8)->Set(kDcaBins, binsDca);
  fListHist->Add(fHistRealTracks);
  //
  //
  fHistPidQA = new TH3F("fHistPidQA","PID QA",500,0.1,10,1000,0,1000,2,-2,2);
  BinLogAxis(fHistPidQA);
  fListHist->Add(fHistPidQA);
  //
  fHistTofQA = new TH2F("fHistTofQA","TOF-QA",200,-4.,4.,400,0,1.1);
  fListHist->Add(fHistTofQA);
  //                            0,            1,           2,  3,      4,            5,    6,    7,        8,    9
  Int_t    binsHistMC[10] = {   3,    kMultBins,     kPtBins,  2,    kYBins      ,  50,    2,  100, kDcaBins,    6};
  Double_t xminHistMC[10] = {-0.5,         -0.5,           0, -2,    kYlowBorder ,  -5,- 0.5, -2.5,       -3, -0.5};
  Double_t xmaxHistMC[10] = { 2.5,         10.5,           3,  2,    kYhighBorder,   5,  1.5,  2.5,        3,  5.5};
  //
  // different binning for CODE axis, if we want to save motherPDG
  //
  fHistMCparticles = new THnSparseF("fHistMCparticles","MC histogram",10,binsHistMC,xminHistMC,xmaxHistMC);
  fHistMCparticles->GetAxis(2)->Set(kPtBins, binsPt);
  fHistMCparticles->GetAxis(8)->Set(kDcaBins, binsDca);
  fListHist->Add(fHistMCparticles);
  //
  fHistMult = new TH2F("fHistMult", "control histogram to count number of events", 502, -2.5, 499.5,4,-0.5,3.5);
  fHistCentrality = new TH1F("fHistCentrality", "control histogram to count number of events", 22, -1.5, 20.5);
  fListHist->Add(fHistMult);
  fListHist->Add(fHistCentrality);
  //
  fHistMomCorr = new TH2F("fHistMomCorr","momentum correction; pT_{MCtrue}; rel. difference",200,0,3,200,-0.2,0.2);
  fListHist->Add(fHistMomCorr);
  //
  fHistEtaPtGen = new TH3F("fHistPtEtaGen","rapidity correction; p_{T} (GeV/c); rapidity y; pid",200,0,12, 200,-1.2,1.2, 3,-0.5,2.5);
  fListHist->Add(fHistEtaPtGen);
  //
  fHistVertex = new TH2F("fHistVertex","vertex position; V_{XY} (cm); V_{Z}",2000,-1.,1.,400,-20.,20.);
  fListHist->Add(fHistVertex);
  //
  fHistVertexRes = new TH1F("fHistVertexRes","vertex position difference MC truth and rec; V_{Z,rec} - V_{Z,truth} (cm)",2000,-0.5,0.5);
  fListHist->Add(fHistVertexRes);
  //
  fHistVertexResTracks = new TH2F("fHistVertexResTracks","vertex res vs no of tracks; no. of ch. primaries |#eta| < 0.9; V_{Z,rec} - V_{Z,truth} (cm)",50,-0.5,499.5,200,-0.05,0.05);
  fListHist->Add(fHistVertexResTracks);

  PostData(1, fListHist);

}

//________________________________________________________________________
void AliAnalysisDeuteronpA::UserExec(Option_t *) 
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
  
  //check with analysis utils
  //first event in chunck --> continue
  if(fUtils->IsFirstEventInChunk(fESD)) {
    PostData(1, fListHist);
    return;
  }
  //check for pileup
  if (fUtils->IsPileUpEvent(fESD)){
    PostData(1, fListHist);
    return;
  }
  //vetex cuts
  Bool_t isVtxOk = kTRUE;
  if(!fUtils->IsVertexSelected2013pA(fESD)) isVtxOk = kFALSE;
  //NEED TO PUT THIS IN ACTION SOMEWHERE


  if (!fESDtrackCuts) {
    Printf("ERROR: fESDtrackCuts not available");
    return;
  }
  //
  // check if event is selected by physics selection class
  //
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
  // fill histogram to monitor vertex position
  if (vertex && isVtxOk) fHistVertex->Fill(TMath::Sqrt(vertex->GetXv()*vertex->GetXv() + vertex->GetYv()*vertex->GetYv()),vertex->GetZv());

  Float_t vertexRes = -100;
  Int_t nPrimaries = 0;


  if (fMCtrue && vertex && isVtxOk) {
    //
    //
    for(Int_t i = 0; i < stack->GetNtrack(); i++) {
      TParticle * trackMC = stack->Particle(i);
      
      if (!trackMC->IsPrimary()) continue;
      //
      //Double_t xv = trackMC->Vx();
      //Double_t yv = trackMC->Vy();
      Double_t zv = trackMC->Vz();

      vertexRes = vertex->GetZv() - zv;
      fHistVertexRes->Fill(vertex->GetZv() - zv);

      break;
    }

    for(Int_t i = 0; i < stack->GetNtrack(); i++) {
      TParticle * trackMC = stack->Particle(i);
      
      //      if (trackMC->IsPrimary()) {
      if (stack->IsPhysicalPrimary(i)) {
	if (trackMC->Eta() > -0.9 && trackMC->Eta() < 0.9 && trackMC->GetPDG()->Charge() != 0){
	  nPrimaries++;
	  //cout << nPrimaries << endl;
	}
      }
    }
  }

  fHistVertexResTracks->Fill(nPrimaries, vertexRes);




  //
  Float_t centrality = -1;
  //
  /*
  Int_t rootS = fESD->GetBeamEnergy() < 1000 ? 0 : 1;
  if (fESD->GetEventSpecie() == 4) { // PbPb
    rootS = 2;
    AliCentrality *esdCentrality = fESD->GetCentrality();
    centrality = esdCentrality->GetCentralityClass10("V0M") + 1; // centrality percentile determined with V0
    if (TMath::Abs(centrality - 1) < 1e-5) {
      centrality = esdCentrality->GetCentralityClass5("V0M");
    }
  }
  */
  AliCentrality *esdCentrality = fESD->GetCentrality();
  centrality = esdCentrality->GetCentralityClass10("V0A") + 1; // centrality percentile determined with V0A
  if (TMath::Abs(centrality - 1) < 1e-5) {
    centrality = esdCentrality->GetCentralityClass5("V0A");
  }


  //
  Int_t processCode = 0;
  //
  // important change: fill generated only after vertex cut in case of heavy-ions
  //
  if ( fESD->GetEventSpecie() == 4) {
    if (!vertex || !isVtxOk) {
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
  }
  //
  if (fMCtrue) {
    //
    //
    for(Int_t i = 0; i < stack->GetNtrack(); i++) {
      TParticle * trackMC = stack->Particle(i);
      Int_t pdg = trackMC->GetPdgCode();
      //
      Double_t xv = trackMC->Vx();
      Double_t yv = trackMC->Vy();
      //Double_t zv = trackMC->Vz();
      Double_t dxy = 0;
      dxy = TMath::Sqrt(xv*xv + yv*yv); // so stupid to avoid warnings
      //Double_t dz = 0;
      //dz = TMath::Abs(zv); // so stupid to avoid warnings
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
      if (fRapCMSpA) rap = rap + 0.465;
      Double_t pT  = trackMC->Pt();
      Int_t sign = pdg < 0 ? -1 : 1; // only works for charged pi,K,p !!
      //      Double_t transMass = TMath::Sqrt(trackMC->Pt()*trackMC->Pt() + trackMC->GetMass()*trackMC->GetMass()) - trackMC->GetMass();
      //
      Int_t iPart = -1;
      if (TMath::Abs(pdg) == 1000010020)  iPart = 0; // select d+/d- only
      if (TMath::Abs(pdg) == 1000010030)  iPart = 1; // select t+/t- only
      if (TMath::Abs(pdg) == 1000020030)  iPart = 2; // select He+/He- only
      if (iPart == -1) continue;
      //
      Double_t vecHistMC[10] = {iPart, centrality,  pT, sign, rap, 0, 1, 0, dxy, 0};
      fHistMCparticles->Fill(vecHistMC);
      //
      if (iPart == 0) fHistEtaPtGen->Fill(pT, trackMC->Y(), 0.);
      if (iPart == 1) fHistEtaPtGen->Fill(pT, trackMC->Y(), 1.);
      if (iPart == 2) fHistEtaPtGen->Fill(pT, trackMC->Y(), 2.);
    }
  }
  //
  // check vertex
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
  fHistMult->Fill(fESDtrackCuts->CountAcceptedTracks(fESD), 0);
  fHistMult->Fill(fESD->GetNumberOfTracks(), 1);
  fHistCentrality->Fill(centrality);
  //
  //***************************************************
  // track loop
  //***************************************************
  //const Float_t kNsigmaCut = 3;
  //const Float_t k2sigmaCorr = 1/(0.5*(TMath::Erf(kNsigmaCut/sqrt(2))-TMath::Erf(-kNsigmaCut/sqrt(2))))/*1/0.9545*/;
  //
  Float_t dca[2], cov[3]; // dca_xy, dca_z, sigma_xy, sigma_xy_z, sigma_z for the vertex cut
  //
  for (Int_t i=0;i<fESD->GetNumberOfTracks();++i) {
    //
    AliESDtrack *track  =fESD->GetTrack(i); 
    if (!track->GetInnerParam()) continue;
    //
    Double_t ptot = track->GetInnerParam()->GetP(); // momentum for dEdx determination
    //
    // momentum correction for different mass assumption in tracking
    //
    Double_t pT = track->Pt()/(1 - 0.333303/TMath::Power(track->Pt() + 0.651111, 5.27268));
    // -[0]/TMath::Power(x - [1],[2])
    //
    track->GetImpactParameters(dca, cov);
    //
    // cut for dead regions in the detector
    // if (track->Eta() > 0.1 && (track->Eta() < 0.2 && track->Phi() > 0.1 && track->Phi() < 0.1) continue;
    //
    // 2.a) apply some standard track cuts according to general recommendations
    //
    if (!fESDtrackCuts->AcceptTrack(track)) continue;
    //
    UInt_t status = track->GetStatus();
    //
    Bool_t hasTOFout  = status&AliESDtrack::kTOFout; 
    Bool_t hasTOFtime = status&AliESDtrack::kTIME;
    Bool_t hasTOF     = kFALSE;
    if (hasTOFout) hasTOF = kTRUE;
    Float_t length = track->GetIntegratedLength(); 
    if (length < 350.) hasTOF = kFALSE;
    //
    // calculate rapidities and kinematics
    // 
    //
    Double_t pvec[3];
    track->GetPxPyPz(pvec);
    Double_t energyDeuteron = TMath::Sqrt(pT*pT + pvec[2]*pvec[2] +
					  AliPID::ParticleMass(AliPID::kDeuteron)*AliPID::ParticleMass(AliPID::kDeuteron));
    Double_t energyTriton = TMath::Sqrt(track->GetP()*track->GetP() + 
					AliPID::ParticleMass(AliPID::kTriton)*AliPID::ParticleMass(AliPID::kTriton));
    Double_t energyHe3 = TMath::Sqrt(track->GetP()*2*track->GetP()*2 + 
					  AliPID::ParticleMass(AliPID::kHe3)*AliPID::ParticleMass(AliPID::kHe3));
    //
    Double_t rapDeuteron = 0.5*TMath::Log((energyDeuteron + pvec[2])/(energyDeuteron - pvec[2]));
    Double_t rapTriton = 0.5*TMath::Log((energyTriton + pvec[2])/(energyTriton - pvec[2]));
    Double_t rapHe3 = 0.5*TMath::Log((energyHe3 + pvec[2]*2)/(energyHe3 - pvec[2]*2));
    //
    if (fRapCMSpA) {
      rapDeuteron = rapDeuteron + 0.465;
      rapTriton = rapTriton + 0.465;
      rapHe3 = rapHe3 + 0.465;
    }
    //
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
    //
    Float_t deutExp = AliExternalTrackParam::BetheBlochAleph(ptot/(0.938*2),1.45802,27.4992,4.00313e-15,2.48485,8.31768);
    Float_t tritExp = AliExternalTrackParam::BetheBlochAleph(ptot/(0.938*3),1.45802,27.4992,4.00313e-15,2.48485,8.31768);
    Float_t hel3Exp = 4*AliExternalTrackParam::BetheBlochAleph(2*ptot/(0.938*3),1.74962,27.4992,4.00313e-15,2.42485,8.31768);
    if (fMCtrue && ptot > 0.3) {
      //Double_t parMC[5] = {1.17329, 27.4992, 4.00313e-15, 2.35563, 9.47569}; // OLD FOR LHC11b9_1 !!
      //Double_t parMC[5] = {1.17329, 27.4992, 4.00313e-15, 2.1204316, 4.1373729};  // NEW!!! PbPb
      Double_t parMC[5] = {20.1533, 2.58127, 0.00114169, 2.0373, 0.502123}; //for LHC13b2efix enhanced
      deutExp = AliExternalTrackParam::BetheBlochAleph(ptot/(0.938*2),parMC[0],parMC[1],parMC[2],parMC[3],parMC[4]);
      tritExp = AliExternalTrackParam::BetheBlochAleph(ptot/(0.938*3),parMC[0],parMC[1],parMC[2],parMC[3],parMC[4]);
      hel3Exp = 4*AliExternalTrackParam::BetheBlochAleph(2*ptot/(0.938*3),parMC[0],parMC[1],parMC[2],parMC[3],parMC[4]);
    }
    //
    //
    Double_t rap[4] = {rapDeuteron, rapTriton, rapHe3,0};
    Double_t pullsTPC[4] = {(tpcSignal - deutExp)/(0.07*deutExp),
			    (tpcSignal - tritExp)/(0.07*tritExp),
			    (tpcSignal - hel3Exp)/(0.07*hel3Exp),
			    0};
    //
    // Process TOF information
    //
    //Float_t time0 = fESDpid->GetTOFResponse().GetTimeZero();
    Float_t time0 = fESDpid->GetTOFResponse().GetStartTime(track->P());//fESDpid->GetTOFResponse().GetTimeZero();
    Float_t mass = 0;
    Float_t time = -1; 
    Float_t beta = 0;
    //
    if (length > 350.) {
      time = track->GetTOFsignal() - time0;
      if (time > 0) {
	beta = length / (2.99792457999999984e-02 * time);
	Float_t gamma = 1/TMath::Sqrt(1 - beta*beta);
	mass = ptot/TMath::Sqrt(gamma*gamma - 1); // using inner TPC mom. as approx.
      }
    }
    //
    /*
    Double_t timeDeuteron = (length/2.99792457999999984e-02) * TMath::Sqrt(1 + 1.876*1.876/(ptot*ptot));
    Double_t timeTriton   = (length/2.99792457999999984e-02) * TMath::Sqrt(1 + 2.809*2.809/(ptot*ptot));
    Double_t timeHe3      = (length/2.99792457999999984e-02) * TMath::Sqrt(1 + 2.809*2.809/(2*ptot*2*ptot));
    */
    //
    Double_t pullsTOF[4] ={0.,0.,0.,0.};
    pullsTOF[0] = mass*mass - 1.876*1.876; // assuming 130.ps time resolution
    pullsTOF[1] = mass*mass - 2.809*2.809;
    pullsTOF[2] = mass*mass - 2.809*2.809;
    pullsTOF[3] = 0;
    //
    //
    if (hasTOFout && hasTOFtime && TMath::Abs(pullsTPC[1]) < 3) fHistTofQA->Fill(ptot*sign, beta);
    //
    //
    for(Int_t iPart = 0; iPart < 3; iPart++) { // loop over assumed particle type
      //
      // temporary measure to reduce memory: fill only deuterons
      //
      if (iPart != 0) continue;
      //                              0,           1,    2,    3,           4,               5,      6,              7,     8
      Double_t vecHistReal[9]  = {iPart,  centrality,   pT, sign,  rap[iPart], pullsTPC[iPart], hasTOF, pullsTOF[iPart], dca[0]};
      fHistRealTracks->Fill(vecHistReal);
      if (TMath::Abs(rap[1]) < 0.5) fHistPidQA->Fill(ptot,tpcSignal,sign); // has to be used for the very difficult deuterons.
      //
      // using MC truth for precise efficiencies...
      //
      if (fMCtrue) {
	Int_t code = 9; // code: 0-generated, 1-true rec. primaries, 2-misident, 3-second weak, 4-second material
	Int_t assumedPdg = 0;//2212(proton); 321(Kaon); 211(pion);
	//
	if (iPart == 0) assumedPdg = 1000010020;
	if (iPart == 1) assumedPdg = 1000010030;
	if (iPart == 2) assumedPdg = 1000020030;
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
	}
	if (pdg == assumedPdg && stack->IsSecondaryFromMaterial(TMath::Abs(track->GetLabel()))) code = 4;
	//
	// fill momentum correction histogram
	//
	if (code == 1) fHistMomCorr->Fill(trackMC->Pt(), (pT - trackMC->Pt())/trackMC->Pt());
	//
	//here cout rap[ipart] and trackMC->Y() and pT
	//if (trackMC->Pt()>1.5 && TMath::Abs(rap[iPart] - trackMC->Y()) > 0.1) {
	//    printf("Y_rec: %4.3f    Y_MC: %4.3f    p_T: %5.3f \n", rap[iPart], trackMC->Y(), trackMC->Pt());
	//}
	 //
	// check TOF mismatch on MC basis with TOF label
	//
	Int_t tofLabel[3] = {0,0,0};
	track->GetTOFLabel(tofLabel);
	//	
	//
	if (TMath::Abs(track->GetLabel()) != TMath::Abs(tofLabel[0]) || tofLabel[1] > 0) {
	  hasTOF = kFALSE;
	  if (TMath::Abs(trackMC->GetFirstMother()) == TMath::Abs(tofLabel[0])) hasTOF = kTRUE;
	}
	//
	// IMPORTANT BIG PROBLEM HERE THE PROBABLILITY TO HAVE A PID SIGNAL MUST BE IN !!!!!!!!!!!!
	//
	//                              0,           1,   2,    3,           4,               5,      6,               7,      8,   9
	Double_t vectorHistMC[10] = {iPart,  centrality,  pT, sign,  rap[iPart], pullsTPC[iPart], hasTOF, pullsTOF[iPart], dca[0], code};
	fHistMCparticles->Fill(vectorHistMC);

      }
      //
      //
    } // end loop over assumed particle type


  } // end of track loop
  
  // Post output data  
  PostData(1, fListHist);
  
}      


//________________________________________________________________________
void AliAnalysisDeuteronpA::Terminate(Option_t *) 
{
  // Draw result to the screen
  // Called once at the end of the query
  Printf("*** CONSTRUCTOR CALLED ****");

}


//________________________________________________________________________
void AliAnalysisDeuteronpA::BinLogAxis(const TH1 *h) {
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

