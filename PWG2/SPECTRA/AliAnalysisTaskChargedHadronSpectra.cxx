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
 **************************************************************************/

///////////////////////////////////////////////////////////////////////////
//                                                                       //
//                                                                       //
// Analysis for identified particle spectra measured with TPC dE/dx.     //
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

#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliStack.h"

#include "AliLog.h"

#include "AliAnalysisTaskChargedHadronSpectra.h"


ClassImp(AliAnalysisTaskChargedHadronSpectra)

//________________________________________________________________________
AliAnalysisTaskChargedHadronSpectra::AliAnalysisTaskChargedHadronSpectra() 
: AliAnalysisTaskSE("TaskChargedHadron"), fESD(0), fListHist(0), fESDtrackCuts(0),fESDpid(0),
  fMCtrue(0),
  fTrackingMode(0),
  fAlephParameters(),
  fHistPtMCKaon(0),
  fHistPtMCProton(0),
  fHistPtMCPion(0),
  fHistPtEtaKaon(0),
  fHistPtEtaKaonNoKink(0),
  fHistPtEtaProton(0),
  fHistPtEtaProtonDCA(0),
  fHistPtEtaPion(0),
  fHistClassicKaon(0),
  fHistClassicProton(0),
  fHistClassicPion(0),
  fDeDx(0),
  fHistTrackPerEvent(0),
  fHistTrackPerEventMC(0),
  fSecProtons(0),
  fVertexZ(0),
  fHistEtaNcls(0),
  fHistEtaPhi(0),
  fHistEffProton(0),
  fHistEffProtonDCA(0),
  fHistEffPion(0),
  fHistEffKaon(0),
  fHistRealTracks(0),
  fHistMCparticles(0)  
{
  // default Constructor
}


//________________________________________________________________________
AliAnalysisTaskChargedHadronSpectra::AliAnalysisTaskChargedHadronSpectra(const char *name) 
  : AliAnalysisTaskSE(name), fESD(0), fListHist(0), fESDtrackCuts(0),fESDpid(0),
    fMCtrue(0),
    fTrackingMode(0),
    fAlephParameters(),
    fHistPtMCKaon(0),
    fHistPtMCProton(0),
    fHistPtMCPion(0),
    fHistPtEtaKaon(0),
    fHistPtEtaKaonNoKink(0),
    fHistPtEtaProton(0),
    fHistPtEtaProtonDCA(0),
    fHistPtEtaPion(0),
    fHistClassicKaon(0),
    fHistClassicProton(0),
    fHistClassicPion(0),
    fDeDx(0),
    fHistTrackPerEvent(0),
    fHistTrackPerEventMC(0),
    fSecProtons(0),
    fVertexZ(0),
    fHistEtaNcls(0),
    fHistEtaPhi(0),
    fHistEffProton(0),
    fHistEffProtonDCA(0),
    fHistEffPion(0),
    fHistEffKaon(0),
    fHistRealTracks(0),
    fHistMCparticles(0)
{
  //
  // standard constructur which should be used
  //
  Printf("*** CONSTRUCTOR CALLED ****");
  //
  fMCtrue = kFALSE; // change three things for processing real data: 1. set this to false! / 2. change ALEPH parameters / 3. change parameters in the Exec()
  fTrackingMode = 0;
  /* real */
  fAlephParameters[0] = 0.0283086;
  fAlephParameters[1] = 2.63394e+01;
  fAlephParameters[2] = 5.04114e-11;
  fAlephParameters[3] = 2.12543e+00;
  fAlephParameters[4] = 4.88663e+00;
  //
  // initialize PID object
  //
  fESDpid = new AliESDpid();
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
void AliAnalysisTaskChargedHadronSpectra::Initialize()
{
  //
  // updating parameters in case of changes
  //
  // 1. pid parameters
  //
  fESDpid->GetTPCResponse().SetBetheBlochParameters(fAlephParameters[0],fAlephParameters[1],fAlephParameters[2],fAlephParameters[3],fAlephParameters[4]);
  //
  // 2. track cuts
  //
  //fESDtrackCuts->SetMaxCovDiagonalElements(2, 2, 0.5, 0.5, 2);  // BEWARE STANDARD VALUES ARE: 2, 2, 0.5, 0.5, 2
  //fESDtrackCuts->SetMaxNsigmaToVertex(3);
  //fESDtrackCuts->SetRequireSigmaToVertex(kTRUE);
  fESDtrackCuts->SetAcceptKinkDaughters(kFALSE);
  fESDtrackCuts->SetMinNClustersTPC(80);
  fESDtrackCuts->SetMaxChi2PerClusterTPC(4);
  //fESDtrackCuts->SetMaxDCAToVertexXY(3);
  //fESDtrackCuts->SetMaxDCAToVertexZ(3);
  //
  if (fTrackingMode == 1) {//for GLOBAL running IN ADDITION
    fESDtrackCuts->SetRequireTPCRefit(kTRUE);
    fESDtrackCuts->SetRequireITSRefit(kTRUE);
    fESDtrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kAny); //TEMPORARY <-> REMOVE
    //
    //fESDtrackCuts->SetMinNClustersITS(3);
    //
    //fESDtrackCuts->SetMaxDCAToVertexXY(0.5);
    //fESDtrackCuts->SetMaxDCAToVertexZ(0.05); // this is now done in this special function...
  }


}


//________________________________________________________________________
void AliAnalysisTaskChargedHadronSpectra::UserCreateOutputObjects() 
{
  // Create histograms
  // Called once
  fListHist = new TList();
  //fListHist->SetOwner(); // Whoever knows how the data handling is ...?

  const Int_t kPtBins = 2*30;
  const Double_t kPtMax  = 2.0;
  const Int_t kEtaBins = 4;
  const Double_t kEtaMax = 1;
  const Int_t kDeDxBins = 100;
  const Double_t kDeDxMax = 1;
  const Int_t kMultBins = 10;
  const Int_t kMultMax = 300;

  Double_t binsPt[kPtBins+1] = {-2.00, -1.90, -1.80, -1.70, -1.60, -1.50, -1.40, -1.30, -1.20, -1.10,
				-1.00, -0.95, -0.90, -0.85, -0.80, -0.75, -0.70, -0.65, -0.60, -0.55,
				-0.50, -0.45, -0.40, -0.35, -0.30, -0.25, -0.20, -0.15, -0.10, -0.05,
				 0.00,  0.05,  0.10,  0.15,  0.20,  0.25,  0.30,  0.35,  0.40,  0.45,
				 0.50,  0.55,  0.60,  0.65,  0.70,  0.75,  0.80,  0.85,  0.90,  0.95,
				 1.00,  1.10,  1.20,  1.30,  1.40,  1.50,  1.60,  1.70,  1.80,  1.90, 2.00};

  // MC histograms
  fHistPtMCKaon = new TH3F("HistPtMCKaon", "PtEtaKaon; mult; #eta; p_{T} (GeV)",kMultBins,-0.5,kMultMax,kEtaBins,0,kEtaMax,kPtBins,-kPtMax,kPtMax);
  fHistPtMCProton = new TH3F("HistPtMCProton", "PtEtaProton;mult; #eta; p_{T} (GeV)",kMultBins,-0.5,kMultMax,kEtaBins,0,kEtaMax,kPtBins,-kPtMax,kPtMax);
  fHistPtMCPion = new TH3F("HistPtMCPion", "PtEtaPion;mult; #eta; p_{T} (GeV)",kMultBins,-0.5,kMultMax,kEtaBins,0,kEtaMax,kPtBins,-kPtMax,kPtMax);
  //
  fHistPtMCKaon->GetZaxis()->Set(kPtBins, binsPt);
  fHistPtMCProton->GetZaxis()->Set(kPtBins, binsPt);
  fHistPtMCPion->GetZaxis()->Set(kPtBins, binsPt);
  //
  fListHist->Add(fHistPtMCKaon);
  fListHist->Add(fHistPtMCProton);
  fListHist->Add(fHistPtMCPion);
  //
  // reconstructed particle histograms
  fHistPtEtaKaon = new TH3F("HistPtEtaKaon", "PtEtaKaon;mult; #eta; p_{T} (GeV)",kMultBins,-0.5,kMultMax,kEtaBins,0,kEtaMax,kPtBins,-kPtMax,kPtMax);
  fHistPtEtaKaonNoKink = new TH3F("HistPtEtaKaonNoKink", "PtEtaKaonNoKink;mult; #eta; p_{T} (GeV)",kMultBins,-0.5,kMultMax,kEtaBins,0,kEtaMax,kPtBins,-kPtMax,kPtMax);
  fHistPtEtaProton = new TH3F("HistPtEtaProton", "PtEtaProton;mult; #eta; p_{T} (GeV)",kMultBins,-0.5,kMultMax,kEtaBins,0,kEtaMax,kPtBins,-kPtMax,kPtMax);
  fHistPtEtaProtonDCA = new TH3F("HistPtEtaProtonDCA", "PtEtaProton;DCA (cm); #eta; p_{T} (GeV)",1000,0,50,kEtaBins,0,kEtaMax,kPtBins,-kPtMax,kPtMax);
  fHistPtEtaPion = new TH3F("HistPtEtaPion", "PtEtaPion;mult; #eta; p_{T} (GeV)",kMultBins,-0.5,kMultMax,kEtaBins,0,kEtaMax,kPtBins,-kPtMax,kPtMax);
  //
  fHistPtEtaKaon->GetZaxis()->Set(kPtBins, binsPt);
  fHistPtEtaKaonNoKink->GetZaxis()->Set(kPtBins, binsPt);
  fHistPtEtaProton->GetZaxis()->Set(kPtBins, binsPt);  
  fHistPtEtaProtonDCA->GetZaxis()->Set(kPtBins, binsPt);
  fHistPtEtaPion->GetZaxis()->Set(kPtBins, binsPt);
  //
  fListHist->Add(fHistPtEtaKaon);
  fListHist->Add(fHistPtEtaKaonNoKink);
  fListHist->Add(fHistPtEtaProton);
  fListHist->Add(fHistPtEtaProtonDCA);
  fListHist->Add(fHistPtEtaPion);
 
  // histograms for the classical analysis
  fHistClassicKaon = new TH3F("HistClassicKaon", "PtEtaKaon;p_{T} (GeV);#eta;delta dEdx",kPtBins,-kPtMax,kPtMax,kEtaBins,0,kEtaMax,kDeDxBins,-kDeDxMax,kDeDxMax);
  fHistClassicProton = new TH3F("HistClassicProton", "PtEtaProton;p_{T} (GeV);#eta;delta dEdx",kPtBins,-kPtMax,kPtMax,kEtaBins,0,kEtaMax,kDeDxBins,-kDeDxMax,kDeDxMax);
  fHistClassicPion = new TH3F("HistClassicPion", "PtEtaPion;p_{T} (GeV);#eta;delta dEdx",kPtBins,-kPtMax,kPtMax,kEtaBins,0,kEtaMax,kDeDxBins,-kDeDxMax,kDeDxMax);
  //
  fHistClassicKaon->GetXaxis()->Set(kPtBins, binsPt);
  fHistClassicProton->GetXaxis()->Set(kPtBins, binsPt);
  fHistClassicPion->GetXaxis()->Set(kPtBins, binsPt);
  //
  fListHist->Add(fHistClassicKaon);
  fListHist->Add(fHistClassicProton);
  fListHist->Add(fHistClassicPion);
  // histograms of general interest
  fDeDx = new TH3F("DeDx","dEdx; momentum p (GeV); TPC signal (a.u.)",500,0.1,100.,1000,0.,1000,2,-2,2);
  BinLogX(fDeDx);
  fListHist->Add(fDeDx);
  fHistTrackPerEvent = new TH2F("HistTrackPerEvent", "Tracks per event;Number of Tracks;Number of Events; code",101,-0.5,100,10,-0.5,9.5);
  fListHist->Add(fHistTrackPerEvent);
  fHistTrackPerEventMC  = new TH3F("HistTrackPerEventMC", "Tracks per event;code;TrackPerEvent;IsSelected",10,-0.5,9,101,-0.5,100,2,-0.5,1.5);
  fListHist->Add(fHistTrackPerEventMC);
  fSecProtons = new TH2F("SecProtons", "xy;x;y",100,-10,10,100,-10,10);
  fListHist->Add(fSecProtons);
  fVertexZ = new TH3F("VertexZ", "vertex position;x (cm); y (cm); z (cm); counts",100,-5,5,100,-5,5,100,-25,25);
  fListHist->Add(fVertexZ);
  //
  fHistEtaNcls = new TH2F("HistEtaNcls", "EtaNcls;#eta;Ncls",100,-1.5,1.5,80,0,160);
  fListHist->Add(fHistEtaNcls);
  fHistEtaPhi = new TH2F("HistEtaPhi", "EtaNcls;#eta;#phi",100,-1.5,1.5,100,0,2*TMath::Pi());
  fListHist->Add(fHistEtaPhi);

  // histograms for a refined efficiency study
  fHistEffProton = new TH3F("HistEffProton", "Eff;p_{T} (GeV); code",10,-0.5,9.5,kEtaBins,0,kEtaMax,kPtBins,-kPtMax,kPtMax);
  fListHist->Add(fHistEffProton);
  fHistEffProtonDCA  = new TH3F("HistEffProtonDCA", "Eff;p_{T} (GeV); code",10,-0.5,9.5,200,0,15,kPtBins,-kPtMax,kPtMax);
  fListHist->Add(fHistEffProtonDCA);
  fHistEffPion = new TH3F("HistEffPion", "Eff;p_{T} (GeV); code",10,-0.5,9.5,kEtaBins,0,kEtaMax,kPtBins,-kPtMax,kPtMax);
  fListHist->Add(fHistEffPion);
  fHistEffKaon = new TH3F("HistEffKaon", "Eff;p_{T} (GeV); code",10,-0.5,9.5,kEtaBins,0,kEtaMax,kPtBins,-kPtMax,kPtMax);
  fListHist->Add(fHistEffKaon);
  //
  fHistEffProton->GetZaxis()->Set(kPtBins, binsPt);
  fHistEffProtonDCA->GetZaxis()->Set(kPtBins, binsPt);
  fHistEffPion->GetZaxis()->Set(kPtBins, binsPt);
  fHistEffKaon->GetZaxis()->Set(kPtBins, binsPt);
  //
  // create the histograms with all necessary information
  //
  //  (0.) multiplicity, (1.) eta, (2.) pT, (3.) sign, (4.) rapPion, (5.) rapKaon, (6.) rapProton, (7.) dca, (8.) nclIter1, (9.) pullPion, (10.) pullKaon, (11.) pullProton
  //                               0,   1, 2, 3,  4 ,   5,   6, 7,  8,  9, 10, 11
  Int_t    binsHistReal[12] = {  301,  16,20, 2,  16,  16,  16,60, 10,100,100,100};
  Double_t xminHistReal[12] = { -0.5,-0.8, 0,-2,-0.8,-0.8,-0.8, 0, 60,-10,-10,-10};
  Double_t xmaxHistReal[12] = {300.5, 0.8, 2, 2, 0.8, 0.8, 0.8,15,160, 10, 10, 10};
  fHistRealTracks = new THnSparseS("fHistRealTracks","real tracks; (0.) multiplicity, (1.) eta, (2.) pT, (3.) sign, (4.) rapPion, (5.) rapKaon, (6.) rapProton, (7.) dca, (8.) nclIter1, (9.) pullPion, (10.) pullKaon, (11.) pullProton",12,binsHistReal,xminHistReal,xmaxHistReal);
  fListHist->Add(fHistRealTracks);
  //
  // (0.) multiplicity, (1.) eta, (2.) pT, (3.) sign, (4.) rap, (5.) dca, (6.) ID, (7.) code
  //                            0,   1, 2, 3,   4, 5,   6,  7
  Int_t    binsHistMC[8] = {  301,  16,20, 2,  16,60,   5,  10};
  Double_t xminHistMC[8] = { -0.5,-0.8, 0,-2,-0.8, 0,-0.5,-0.5};
  Double_t xmaxHistMC[8] = {300.5, 0.8, 2, 2, 0.8,15, 4.5, 9.5};
  fHistMCparticles = new THnSparseS("fHistMCparticles"," (0.) multiplicity, (1.) eta, (2.) pT, (3.) sign, (4.) rap, (5.) dca, (6.) ID, (7.) code",8,binsHistMC,xminHistMC,xmaxHistMC);
  fListHist->Add(fHistMCparticles);
  //
  //
  //
 
  
}

//________________________________________________________________________
void AliAnalysisTaskChargedHadronSpectra::UserExec(Option_t *) 
{
  //
  // main event loop
  //
  AliLog::SetGlobalLogLevel(AliLog::kError);
  //
  // Check Monte Carlo information and other access first:
  //
  AliMCEventHandler* eventHandler = dynamic_cast<AliMCEventHandler*> (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler());
  if (!eventHandler) {
    //Printf("ERROR: Could not retrieve MC event handler");
    if (fMCtrue) return;
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
  
  Int_t trackCounter = 0;
  //
  // check if event is selected by physics selection class
  //
  fHistTrackPerEvent->Fill(fESD->GetNumberOfTracks(),0);
  Bool_t isSelected = kFALSE;
  isSelected = ((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected();
  if (isSelected) fHistTrackPerEvent->Fill(fESD->GetNumberOfTracks(),1); 
  //
  // monitor vertex position
  //
  const AliESDVertex *vertex = fESD->GetPrimaryVertexTracks();
  if(vertex->GetNContributors()<1) {
    // SPD vertex
    vertex = fESD->GetPrimaryVertexSPD();
    if(vertex->GetNContributors()<1) vertex = 0x0;
  }  
  
  // 1st track loop to determine multiplicities
  for (Int_t i=0;i<fESD->GetNumberOfTracks();++i) {
    AliESDtrack *track = 0x0;
    if (fTrackingMode == 0) track = AliESDtrackCuts::GetTPCOnlyTrack(fESD,i);
    if (fTrackingMode == 1) track = new AliESDtrack(*fESD->GetTrack(i));  // alternative: global track GLOBAL GLOBAL  GLOBAL GLOBAL  GLOBAL GLOBAL  GLOBAL
    if (!track) continue;
    if (fESDtrackCuts->AcceptTrack(track) &&  TMath::Abs(track->Eta())< 0.9) {
      trackCounter++;
      //
    }
    delete track;
  }
  
  fHistTrackPerEventMC->Fill(0., trackCounter, isSelected);
  fHistTrackPerEvent->Fill(trackCounter,2);  
  
  if (fMCtrue) {
    //
    //
    //
    AliHeader * header = mcEvent->Header();
    Int_t processtype = GetPythiaEventProcessType(header);
    // non diffractive
    if (processtype !=92 && processtype !=93 && processtype != 94) fHistTrackPerEventMC->Fill(1., trackCounter, isSelected);
    // single diffractive
    if ((processtype == 92 || processtype == 93)) fHistTrackPerEventMC->Fill(2., trackCounter, isSelected);
    // double diffractive
    if (processtype == 94) fHistTrackPerEventMC->Fill(3., trackCounter, isSelected);
    //
    Int_t mult = trackCounter;//stack->GetNprimary();
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
      if (pdg == 321) fHistPtMCKaon->Fill(mult, TMath::Abs(trackMC->Y()), trackMC->Pt()); // select K+ only
      if (pdg == 211) fHistPtMCPion->Fill(mult, TMath::Abs(trackMC->Y()), trackMC->Pt()); // select Pi+ only
      if (pdg == 2212) fHistPtMCProton->Fill(mult, TMath::Abs(trackMC->Y()), trackMC->Pt()); // select p+ only
      //
      if (pdg == -321) fHistPtMCKaon->Fill(mult, TMath::Abs(trackMC->Y()), -trackMC->Pt()); // select K- only
      if (pdg == -211) fHistPtMCPion->Fill(mult, TMath::Abs(trackMC->Y()), -trackMC->Pt()); // select Pi- only
      if (pdg == -2212) fHistPtMCProton->Fill(mult, TMath::Abs(trackMC->Y()), -trackMC->Pt()); // select p- only
      //
      // special Kaon checks - those which decay before entering the TPC fiducial volume
      // <-> should be revisited by looping from first to last daugther and check if UniqueID is 4
      //
      if (TMath::Abs(pdg)==321) {
	TParticle * trackDaughterMC = stack->Particle(TMath::Abs(trackMC->GetFirstDaughter()));
	Int_t pdgDaughter = trackDaughterMC->GetPdgCode();
	if (TMath::Abs(pdgDaughter) == 13 && pdg == 321) fHistEffKaon->Fill(5, TMath::Abs(trackMC->Y()), trackMC->Pt()); // Kaon control hist
	if (TMath::Abs(pdgDaughter) == 13 && pdg == -321) fHistEffKaon->Fill(5, TMath::Abs(trackMC->Y()), -trackMC->Pt()); // Kaon control hist
      }
    }
  }
  //
  //end MC tracks loop <-> now check event selection criteria
  //
  if (!isSelected) {
    PostData(1, fListHist);
    return;
  }
  //
  if (!vertex) {
    fHistTrackPerEventMC->Fill(0., trackCounter, isSelected);
    PostData(1, fListHist);
    return;
  } else {
    fVertexZ->Fill(vertex->GetXv(),vertex->GetYv(),vertex->GetZv());
    if (TMath::Abs(vertex->GetZv()) > 10) {
      fHistTrackPerEventMC->Fill(0., trackCounter, isSelected);
      PostData(1, fListHist);
      return;
    }
  }
  //
  // tarck loop
  //
  const Float_t kNsigmaCut = 3;
  const Float_t k2sigmaCorr = 1/(0.5*(TMath::Erf(kNsigmaCut/sqrt(2))-TMath::Erf(-kNsigmaCut/sqrt(2))))/*1/0.9545*/;
  
  Float_t dca[2], cov[3]; // dca_xy, dca_z, sigma_xy, sigma_xy_z, sigma_z for the vertex cut
  
  for (Int_t i=0;i<fESD->GetNumberOfTracks();++i) {
    //
    // 1. select TPC tracks only and propagate them to the primary vertex determined with the TPC standalone
    //
    AliESDtrack *track = 0x0;
    if (fTrackingMode == 0) track = AliESDtrackCuts::GetTPCOnlyTrack(fESD,i);
    if (fTrackingMode == 1) track = new AliESDtrack(*fESD->GetTrack(i)); // alternative: global track GLOBAL GLOBAL  GLOBAL GLOBAL 
    if (!track) continue;
    AliESDtrack *trackDummy = fESD->GetTrack(i);
    if (!trackDummy->GetInnerParam()) { // DEBUG DEBUG
      delete track;
      continue;
    }
    //
    if (track->GetTPCNclsIter1() < 80) {
      delete track;
      continue;
    }
    //
    AliExternalTrackParam trackIn(*trackDummy->GetInnerParam()); 
    Double_t ptot = trackIn.GetP(); // momentum for dEdx determination
    Double_t pT = track->Pt();
    Double_t eta = TMath::Abs(track->Eta());
    //
    fHistEtaNcls->Fill(track->Eta(),track->GetTPCNcls());
    fHistEtaPhi->Fill(track->Eta(),track->Phi());
    //
    // cut for dead regions in the detector
    // if (track->Eta() > 0.1 && (track->Eta() < 0.2 && track->Phi() > 0.1 && track->Phi() < 0.1) continue;
    //
    // 2.a) apply some standard track cuts according to general recommendations
    //
    if (!fESDtrackCuts->AcceptTrack(track)) {
      //
      if (fMCtrue) {
	TParticle *trackMC = stack->Particle(TMath::Abs(track->GetLabel()));
	Int_t pdg = TMath::Abs(trackMC->GetPdgCode());     
	if (pdg == 321 && stack->IsPhysicalPrimary(TMath::Abs(track->GetLabel()))) fHistEffKaon->Fill(6,eta,track->GetSign()*pT); // Kaon control hist
      }
      //
      delete track;
      continue;
    }
    
    if (fTrackingMode == 1) {
      /*if (!SelectOnImpPar(track)) {
	delete track;
	continue;
	}*/
    }
    
    //
    // calculate rapidities and kinematics
    // 
    //
    Double_t pvec[3];
    track->GetPxPyPz(pvec);
    Double_t energyPion = TMath::Sqrt(track->GetP()*track->GetP() + AliPID::ParticleMass(AliPID::kPion)*AliPID::ParticleMass(AliPID::kPion));
    Double_t energyKaon = TMath::Sqrt(track->GetP()*track->GetP() + AliPID::ParticleMass(AliPID::kKaon)*AliPID::ParticleMass(AliPID::kKaon));
    Double_t energyProton = TMath::Sqrt(track->GetP()*track->GetP() + AliPID::ParticleMass(AliPID::kProton)*AliPID::ParticleMass(AliPID::kProton));
    //
    Double_t rapPion = TMath::Abs(0.5*TMath::Log((energyPion + pvec[2])/(energyPion - pvec[2])));
    Double_t rapKaon = TMath::Abs(0.5*TMath::Log((energyKaon + pvec[2])/(energyKaon - pvec[2])));
    Double_t rapProton = TMath::Abs(0.5*TMath::Log((energyProton + pvec[2])/(energyProton - pvec[2])));
    //
    // 3. make the PID
    //
    Double_t sign = track->GetSign();   
    Double_t tpcSignal = track->GetTPCsignal();
    //
    if (eta < 0.8 && track->GetTPCsignalN() > 70) fDeDx->Fill(ptot, tpcSignal, sign);
    //
    // 3.a. calculate expected signals
    //
    Float_t sigKaon     = fESDpid->GetTPCResponse().GetExpectedSignal(ptot, AliPID::kKaon);
    Float_t sigProton   = fESDpid->GetTPCResponse().GetExpectedSignal(ptot, AliPID::kProton);
    Float_t sigPion     = fESDpid->GetTPCResponse().GetExpectedSignal(ptot, AliPID::kPion);
    //
    // 3.b. fill histograms
    //
    fHistClassicKaon->Fill(sign*pT,rapKaon,(tpcSignal-sigKaon)/sigKaon);
    fHistClassicProton->Fill(sign*pT,rapProton,(tpcSignal-sigProton)/sigProton);
    fHistClassicPion->Fill(sign*pT,rapPion,(tpcSignal-sigPion)/sigPion);
    //
    /* 2sigma PID with 2sigma eff correction! */
    // PION
    if (TMath::Abs(fESDpid->NumberOfSigmasTPC(track,AliPID::kPion)) < kNsigmaCut) {
      fHistPtEtaPion->Fill(trackCounter, rapPion, sign*pT, k2sigmaCorr);
      //
      if (trackCounter < 300 && fMCtrue) {
	TParticle *trackMC = stack->Particle(TMath::Abs(track->GetLabel()));
	Int_t pdg = TMath::Abs(trackMC->GetPdgCode());
	if (pdg == 211 && stack->IsPhysicalPrimary(TMath::Abs(track->GetLabel()))) fHistEffPion->Fill(0,rapPion,sign*pT,k2sigmaCorr);
	if (pdg == 211 && !stack->IsPhysicalPrimary(TMath::Abs(track->GetLabel()))) {
	  fHistEffPion->Fill(1,rapPion,sign*pT,k2sigmaCorr);
	  if (trackMC->GetFirstMother() > 0) {
	    TParticle *trackMother = stack->Particle(trackMC->GetFirstMother());
	    Int_t code = TMath::Abs(trackMother->GetPdgCode());
	    if (code==3122||code==3222||code==3212||code==3112||code==3322||code==3312||code==3332||code==130||code==310) fHistEffPion->Fill(3,rapPion,sign*pT,k2sigmaCorr);
	  }
	}
	if (TMath::Abs(trackMC->GetPdgCode()) != 211) fHistEffPion->Fill(2,rapPion,sign*pT,k2sigmaCorr);
	if (TMath::Abs(trackMC->GetPdgCode()) == 13)  fHistEffPion->Fill(4,rapPion,sign*pT,k2sigmaCorr);
      }
    }
    // KAON
    if (TMath::Abs(fESDpid->NumberOfSigmasTPC(track,AliPID::kKaon)) < kNsigmaCut) {
      fHistPtEtaKaon->Fill(trackCounter, rapKaon, sign*pT, k2sigmaCorr);
      //
      if (trackCounter < 300 && fMCtrue) {
	TParticle *trackMC = stack->Particle(TMath::Abs(track->GetLabel()));
	Int_t pdg = TMath::Abs(trackMC->GetPdgCode());
	if (pdg == 321 && stack->IsPhysicalPrimary(TMath::Abs(track->GetLabel()))) fHistEffKaon->Fill(0,rapKaon,sign*pT,k2sigmaCorr);
	if (pdg == 321 && !stack->IsPhysicalPrimary(TMath::Abs(track->GetLabel()))) {
	  fHistEffKaon->Fill(1,rapKaon,sign*pT,k2sigmaCorr);
	  if (trackMC->GetFirstMother() > 0) {
	    TParticle *trackMother = stack->Particle(trackMC->GetFirstMother());
	    Int_t code = TMath::Abs(trackMother->GetPdgCode());
	    if (code==3122||code==3222||code==3212||code==3112||code==3322||code==3312||code==3332||code==130||code==310) fHistEffKaon->Fill(3,rapKaon,sign*pT,k2sigmaCorr);
	  }
	}
	if (TMath::Abs(trackMC->GetPdgCode()) != 321) fHistEffKaon->Fill(2,rapKaon,sign*pT,k2sigmaCorr);
      }
    }
    // KAON NO KINK
    if (TMath::Abs(fESDpid->NumberOfSigmasTPC(track,AliPID::kKaon)) < kNsigmaCut && track->GetKinkIndex(0)==0) fHistPtEtaKaonNoKink->Fill(trackCounter, rapKaon, sign*pT, k2sigmaCorr);
    // PROTON
    if (TMath::Abs(fESDpid->NumberOfSigmasTPC(track,AliPID::kProton)) < kNsigmaCut) {
      fHistPtEtaProton->Fill(trackCounter, rapProton, sign*pT, k2sigmaCorr);
      //
      track->GetImpactParameters(dca,cov);
      Float_t primVtxDCA = TMath::Sqrt(dca[0]*dca[0] + dca[1]*dca[1]);
      fHistPtEtaProtonDCA->Fill(primVtxDCA, rapProton, sign*pT, k2sigmaCorr);
      //
      if (trackCounter < 300 && fMCtrue) {
	TParticle *trackMC = stack->Particle(TMath::Abs(track->GetLabel()));
	Int_t pdg = TMath::Abs(trackMC->GetPdgCode());
	if (pdg == 2212 && stack->IsPhysicalPrimary(TMath::Abs(track->GetLabel()))) {
	  fHistEffProton->Fill(0.,rapProton,sign*pT,k2sigmaCorr);
	  if (rapProton < 0.15) fHistEffProtonDCA->Fill(0.,primVtxDCA,sign*pT,k2sigmaCorr);
	}
	if (pdg == 2212 && !stack->IsPhysicalPrimary(TMath::Abs(track->GetLabel()))) {
	  fHistEffProton->Fill(1,rapProton,sign*pT,k2sigmaCorr);
	  if (rapProton < 0.15) fHistEffProtonDCA->Fill(1,primVtxDCA,sign*pT,k2sigmaCorr);
	  if (trackMC->GetFirstMother() > 0) {
	    TParticle *trackMother = stack->Particle(trackMC->GetFirstMother());
	    Int_t code = TMath::Abs(trackMother->GetPdgCode());
	    if (code==3122||code==3222||code==3212||code==3112||code==3322||code==3312||code==3332||code==130||code==310) {
	      fHistEffProton->Fill(3,rapProton,sign*pT,k2sigmaCorr);
	      if (rapProton < 0.15) fHistEffProtonDCA->Fill(3,primVtxDCA,sign*pT,k2sigmaCorr);
	    }
	    if (code!=3122&&code!=3222&&code!=3212&&code!=3112&&code!=3322&&code!=3312&&code!=3332&&code!=130&&code!=310) fSecProtons->Fill(trackMC->Vx(),trackMC->Vy());
	  }
	}
	if (TMath::Abs(trackMC->GetPdgCode()) != 2212) fHistEffProton->Fill(2,rapProton,sign*pT,k2sigmaCorr);
      }
    }
    
    delete track;
    
  } // end of track loop
  
  // Post output data  
  PostData(1, fListHist);
  
}      


//________________________________________________________________________
void AliAnalysisTaskChargedHadronSpectra::Terminate(Option_t *) 
{
  // Draw result to the screen
  // Called once at the end of the query
  Printf("*** CONSTRUCTOR CALLED ****");

}


//________________________________________________________________________
Bool_t AliAnalysisTaskChargedHadronSpectra::SelectOnImpPar(AliESDtrack* t) {
  //
  // cut on transverse impact parameter
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
  if (TMath::Abs(d0z0[0]) < d0max && TMath::Abs(d0z0[1]) < d0maxZ) return kTRUE;
  return kFALSE;
}


//________________________________________________________________________
void AliAnalysisTaskChargedHadronSpectra::BinLogX(const TH1 *h) {

  // Method for the correct logarithmic binning of histograms

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
Int_t AliAnalysisTaskChargedHadronSpectra::GetPythiaEventProcessType(const AliHeader* aHeader, const Bool_t adebug) const {
  //
  // get the process type of the event.
  //

  // can only read pythia headers, either directly or from cocktalil header
  AliGenPythiaEventHeader* pythiaGenHeader = dynamic_cast<AliGenPythiaEventHeader*>(aHeader->GenEventHeader());

  if (!pythiaGenHeader) {

    AliGenCocktailEventHeader* genCocktailHeader = dynamic_cast<AliGenCocktailEventHeader*>(aHeader->GenEventHeader());
    if (!genCocktailHeader) {
      //printf("AliAnalysisTaskChargedHadronSpectra::GetProcessType : Unknown header type (not Pythia or Cocktail). \n");
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
      //printf("AliAnalysisTaskChargedHadronSpectra::GetProcessType : Could not find Pythia header. \n");
      return -1;
    }
  }

  if (adebug) {
    //printf("AliAnalysisTaskChargedHadronSpectra::GetProcessType : Pythia process type found: %d \n",pythiaGenHeader->ProcessType());
  }

  return pythiaGenHeader->ProcessType();
} 


//________________________________________________________________________
TH1D * AliAnalysisTaskChargedHadronSpectra::AnalyseClassicProton(const TH3F * input, Int_t EtaBinLow, Int_t EtaBinHigh,const  Double_t * AlephParams) {
  //
  // The multiple Gauss fitting is happening here...
  // input histogram is of the form (Pt, eta, difference in dEdx)
  //

  TF1 *foProton = new TF1("foProton", "AliExternalTrackParam::BetheBlochAleph(x/0.93827,[0],[1],[2],[3],[4])",0.05,20); 
  TF1 *foPion = new TF1("foPion", "AliExternalTrackParam::BetheBlochAleph(x/0.13957,[0],[1],[2],[3],[4])",0.05,20);
  TF1 *foElec = new TF1("foElec", "AliExternalTrackParam::BetheBlochAleph(x/0.000511,[0],[1],[2],[3],[4])",0.05,20);
  TF1 *foKaon = new TF1("foKaon", "AliExternalTrackParam::BetheBlochAleph(x/0.493677,[0],[1],[2],[3],[4])",0.05,20);
  TF1 *foMuon = new TF1("foMuon", "AliExternalTrackParam::BetheBlochAleph(x/0.105658,[0],[1],[2],[3],[4])",0.05,20);
  //
  foProton->SetParameters(AlephParams);  
  foPion->SetParameters(AlephParams);
  foElec->SetParameters(AlephParams);  
  foKaon->SetParameters(AlephParams);  
  foMuon->SetParameters(AlephParams);

  
  const Double_t kSigma = 0.06; // expected resolution (roughly)

  const Int_t kPtBins = 2*30/*input->GetXaxis()->GetNbins()*/;
  const Double_t kPtMax  = input->GetXaxis()->GetXmax();
  const Int_t kEtaBins = input->GetYaxis()->GetNbins();

  TH1D * histPt = new TH1D("histPt", "Pt; pT (Gev); dN",kPtBins,-kPtMax,kPtMax);

  Double_t binsPt[kPtBins+1] = {-2.00, -1.90, -1.80, -1.70, -1.60, -1.50, -1.40, -1.30, -1.20, -1.10,
				-1.00, -0.95, -0.90, -0.85, -0.80, -0.75, -0.70, -0.65, -0.60, -0.55,
				-0.50, -0.45, -0.40, -0.35, -0.30, -0.25, -0.20, -0.15, -0.10, -0.05,
				 0.00,  0.05,  0.10,  0.15,  0.20,  0.25,  0.30,  0.35,  0.40,  0.45,
				 0.50,  0.55,  0.60,  0.65,  0.70,  0.75,  0.80,  0.85,  0.90,  0.95,
				 1.00,  1.10,  1.20,  1.30,  1.40,  1.50,  1.60,  1.70,  1.80,  1.90, 2.00};

  histPt->GetXaxis()->Set(kPtBins, binsPt);

  TCanvas * canvMany = new TCanvas("canvManyProton","canvManyProton");
  canvMany->Print("canvManyProton.ps[");

  for(Int_t x=1; x < kPtBins+1; x++) {
    Double_t pT = input->GetXaxis()->GetBinCenter(x);
    cout << "Now fitting: " << pT << endl;
    for(Int_t y=1; y < kEtaBins+1; y++) {
      Double_t rapidity =  input->GetYaxis()->GetBinCenter(y);
      Double_t eta = TMath::ASinH(TMath::Sqrt(1 + (0.938*0.938)/(pT*pT))*TMath::SinH(rapidity)); 
      Double_t theta = 2*TMath::ATan(TMath::Exp(-eta));
      Double_t ptot = TMath::Abs(pT/TMath::Sin(theta));
      TH1D *histDeDx = input->ProjectionZ(Form("dedx_%i_%i",x,y),x,x,y,y + EtaBinHigh - EtaBinLow);
      histDeDx->SetTitle(Form("pt_%f",pT));
      Float_t maxYield = histDeDx->GetEntries()/(TMath::Sqrt(2*TMath::Pi())*0.04/histDeDx->GetBinWidth(1));
      //
      TF1 * gausFit = new TF1("gausFit", "gaus(0) + gaus(3) + gaus(6) + gaus(9)",-1/*-7*kSigma*/,1/*7*kSigma*/);
      Double_t parameters[12] = {maxYield/2.,0,kSigma,
				 maxYield/2.,(foPion->Eval(ptot)-foProton->Eval(ptot))/foProton->Eval(ptot),kSigma*(foPion->Eval(ptot)/foProton->Eval(ptot)),
				 maxYield/2.,(foElec->Eval(ptot)-foProton->Eval(ptot))/foProton->Eval(ptot),kSigma*(foElec->Eval(ptot)/foProton->Eval(ptot)),
				 maxYield/2.,(foKaon->Eval(ptot)-foProton->Eval(ptot))/foProton->Eval(ptot),kSigma*(foKaon->Eval(ptot)/foProton->Eval(ptot))};
      gausFit->SetParameters(parameters);
      if (TMath::Abs(pT) < 0.7) {
	gausFit = new TF1("gausFit", "gaus(0)",-1/*-7*kSigma*/,1/*7*kSigma*/);
	gausFit->SetRange(-6*kSigma,6*kSigma);
	gausFit->SetParameters(parameters);
	gausFit->SetParLimits(0,0.,maxYield);
	//for(Int_t ipar = 3; ipar < 12; ipar++) gausFit->FixParameter(ipar,0.);
	gausFit->SetParLimits(1,-0.15,0.15);
	gausFit->SetParLimits(2,0.02,0.18);
	//
      } else {
	Double_t pionPosition = (foPion->Eval(ptot)-foProton->Eval(ptot))/foProton->Eval(ptot);
	gausFit->SetParLimits(3,0.,maxYield);
	gausFit->SetParLimits(4, 0.95*pionPosition, 1.05*pionPosition);
	gausFit->SetParLimits(5,0.8*kSigma*(foPion->Eval(ptot)/foProton->Eval(ptot)),1.2*kSigma*(foPion->Eval(ptot)/foProton->Eval(ptot)));
	//
	Double_t elecPosition = (foElec->Eval(ptot)-foProton->Eval(ptot))/foProton->Eval(ptot);
	gausFit->SetParLimits(6,0.,maxYield);
	gausFit->SetParLimits(7, 0.95*elecPosition, 1.05*elecPosition);
	gausFit->SetParLimits(8,0.8*kSigma*(foElec->Eval(ptot)/foProton->Eval(ptot)),1.2*kSigma*(foElec->Eval(ptot)/foProton->Eval(ptot)));
	//
	Double_t kaonPosition = (foKaon->Eval(ptot)-foProton->Eval(ptot))/foProton->Eval(ptot);
	gausFit->SetParLimits(9,0.,maxYield);
	gausFit->SetParLimits(10, 0.95*kaonPosition, 1.05*kaonPosition);
	gausFit->SetParLimits(11,0.8*kSigma*(foKaon->Eval(ptot)/foProton->Eval(ptot)),1.2*kSigma*(foKaon->Eval(ptot)/foProton->Eval(ptot)));
      }
      histDeDx->Fit("gausFit","QNR");
      gausFit->GetParameters(parameters);
      // visualization      
      if (y == EtaBinLow) {
	canvMany->cd(x);
	gPad->SetLogy();  
	histDeDx->SetMarkerStyle(21);
	histDeDx->SetMarkerSize(0.7);
	histDeDx->Draw("E");
	gausFit->SetLineWidth(2);
	gausFit->Draw("same");
	//
	TF1 gausFit0("gausFit0", "gaus(0)",-1,1);
	TF1 gausFit1("gausFit1", "gaus(0)",-1,1);
	TF1 gausFit2("gausFit2", "gaus(0)",-1,1);
	TF1 gausFit3("gausFit3", "gaus(0)",-1,1);
	gausFit0.SetParameters(parameters[0],parameters[1],parameters[2]);
	gausFit1.SetParameters(parameters[3],parameters[4],parameters[5]);
	gausFit2.SetParameters(parameters[6],parameters[7],parameters[8]);
	gausFit3.SetParameters(parameters[9],parameters[10],parameters[11]);
	//
	gausFit0.SetLineColor(4);
	gausFit1.SetLineColor(2);
	gausFit2.SetLineColor(6);
	gausFit3.SetLineColor(8);
	//
	gausFit0.SetLineWidth(1);
	gausFit1.SetLineWidth(1);
	gausFit2.SetLineWidth(1);
	gausFit3.SetLineWidth(1);
	//
	gausFit0.Draw("same");
	gausFit1.Draw("same"); 
	gausFit2.Draw("same"); 
	gausFit3.Draw("same");
	if (TMath::Abs(pT) < 2) canvMany->Print("canvManyProton.ps");
      }
      
      Double_t yield = gausFit->GetParameter(0)*TMath::Sqrt(2*TMath::Pi())*gausFit->GetParameter(2)/histDeDx->GetBinWidth(1); // area of the gaus fit
      delete gausFit;
      //
      if (y == EtaBinLow && yield > 0) histPt->Fill(pT,yield);
      //delete gausYield;
    }
  }
  canvMany->Print("canvManyProton.ps]");

  TCanvas * canvPt = new TCanvas();
  canvPt->cd();
  histPt->Draw();
  
  return histPt;
}



//________________________________________________________________________
TH1D * AliAnalysisTaskChargedHadronSpectra::AnalyseClassicPion(const TH3F * input, Int_t EtaBinLow, Int_t EtaBinHigh,const  Double_t * AlephParams) {
  //
  // The multiple Gauss fitting is happening here...
  // input histogram is of the form (Pt, eta, difference in dEdx)
  //

  TF1 *foProton = new TF1("foProton", "AliExternalTrackParam::BetheBlochAleph(x/0.93827,[0],[1],[2],[3],[4])",0.05,20); 
  TF1 *foPion = new TF1("foPion", "AliExternalTrackParam::BetheBlochAleph(x/0.13957,[0],[1],[2],[3],[4])",0.05,20);
  TF1 *foElec = new TF1("foElec", "AliExternalTrackParam::BetheBlochAleph(x/0.000511,[0],[1],[2],[3],[4])",0.05,20);
  TF1 *foKaon = new TF1("foKaon", "AliExternalTrackParam::BetheBlochAleph(x/0.493677,[0],[1],[2],[3],[4])",0.05,20);
  TF1 *foMuon = new TF1("foMuon", "AliExternalTrackParam::BetheBlochAleph(x/0.105658,[0],[1],[2],[3],[4])",0.05,20);
  //
  foProton->SetParameters(AlephParams);  
  foPion->SetParameters(AlephParams);
  foElec->SetParameters(AlephParams);  
  foKaon->SetParameters(AlephParams);  
  foMuon->SetParameters(AlephParams);

  
  const Double_t kSigma = 0.06; // expected resolution (roughly)

  const Int_t kPtBins = 2*30/*input->GetXaxis()->GetNbins()*/;
  const Double_t kPtMax  = input->GetXaxis()->GetXmax();
  const Int_t kEtaBins = input->GetYaxis()->GetNbins();

  TH1D * histPt = new TH1D("histPt", "Pt; pT (Gev); dN",kPtBins,-kPtMax,kPtMax);

  Double_t binsPt[kPtBins+1] = {-2.00, -1.90, -1.80, -1.70, -1.60, -1.50, -1.40, -1.30, -1.20, -1.10,
				-1.00, -0.95, -0.90, -0.85, -0.80, -0.75, -0.70, -0.65, -0.60, -0.55,
				-0.50, -0.45, -0.40, -0.35, -0.30, -0.25, -0.20, -0.15, -0.10, -0.05,
				 0.00,  0.05,  0.10,  0.15,  0.20,  0.25,  0.30,  0.35,  0.40,  0.45,
				 0.50,  0.55,  0.60,  0.65,  0.70,  0.75,  0.80,  0.85,  0.90,  0.95,
				 1.00,  1.10,  1.20,  1.30,  1.40,  1.50,  1.60,  1.70,  1.80,  1.90, 2.00};

  histPt->GetXaxis()->Set(kPtBins, binsPt);
  //

  TCanvas * canvMany = new TCanvas("canvManyPion","canvManyPion");
  canvMany->Print("canvManyPion.ps[");

  for(Int_t x=1; x < kPtBins+1; x++) {
    Double_t pT = input->GetXaxis()->GetBinCenter(x);
    for(Int_t y=1; y < kEtaBins+1; y++) {
      Double_t rapidity =  input->GetYaxis()->GetBinCenter(y);
      Double_t eta = TMath::ASinH(TMath::Sqrt(1 + (0.1395*0.1395)/(pT*pT))*TMath::SinH(rapidity));  
      Double_t theta = 2*TMath::ATan(TMath::Exp(-eta));
      Double_t ptot = TMath::Abs(pT/TMath::Sin(theta));
      TH1D *histDeDx = input->ProjectionZ(Form("dedx_%i_%i",x,y),x,x,y,y + EtaBinHigh - EtaBinLow);
      histDeDx->SetTitle(Form("pt_%f",pT));
      Float_t maxYield = histDeDx->GetEntries()/(TMath::Sqrt(2*TMath::Pi())*0.04/histDeDx->GetBinWidth(1));
      //
      TF1 * gausFit = new TF1("gausFit", "gaus(0) + gaus(3) + gaus(6) + gaus(9)",-1/*-7*kSigma*/,1/*7*kSigma*/);
      Double_t parameters[12] = {maxYield/2.,0,kSigma,
				 maxYield/2.,(foProton->Eval(ptot)-foPion->Eval(ptot))/foPion->Eval(ptot),kSigma*(foPion->Eval(ptot)/foPion->Eval(ptot)),
				 maxYield/2.,(foElec->Eval(ptot)-foPion->Eval(ptot))/foPion->Eval(ptot),kSigma*(foElec->Eval(ptot)/foPion->Eval(ptot)),
				 maxYield/2.,(foKaon->Eval(ptot)-foPion->Eval(ptot))/foPion->Eval(ptot),kSigma*(foKaon->Eval(ptot)/foPion->Eval(ptot))};
      gausFit->SetParameters(parameters);
      gausFit->SetParLimits(0,0.,maxYield);
      gausFit->SetParLimits(1,-0.05,0.05);
      gausFit->SetParLimits(2,0.04,0.08);
      //
      Double_t protonPosition = (foProton->Eval(ptot)-foPion->Eval(ptot))/foPion->Eval(ptot);
      gausFit->SetParLimits(3,0.,maxYield);
      gausFit->SetParLimits(4, 0.95*protonPosition, 1.05*protonPosition);
      gausFit->SetParLimits(5,0.8*kSigma*(foProton->Eval(ptot)/foPion->Eval(ptot)),1.2*kSigma*(foProton->Eval(ptot)/foPion->Eval(ptot)));
      //
      Double_t elecPosition = (foElec->Eval(ptot)-foPion->Eval(ptot))/foPion->Eval(ptot);
      gausFit->SetParLimits(6,0.,maxYield);
      gausFit->SetParLimits(7, 0.95*elecPosition, 1.05*elecPosition);
      gausFit->SetParLimits(8,0.8*kSigma*(foElec->Eval(ptot)/foPion->Eval(ptot)),1.2*kSigma*(foElec->Eval(ptot)/foPion->Eval(ptot)));
      //
      Double_t kaonPosition = (foKaon->Eval(ptot)-foPion->Eval(ptot))/foPion->Eval(ptot);
      gausFit->SetParLimits(9,0.,maxYield);
      gausFit->SetParLimits(10, 0.95*kaonPosition, 1.05*kaonPosition);
      gausFit->SetParLimits(11,0.8*kSigma*(foKaon->Eval(ptot)/foPion->Eval(ptot)),1.2*kSigma*(foKaon->Eval(ptot)/foPion->Eval(ptot)));
      histDeDx->Fit("gausFit","QNR");
      gausFit->GetParameters(parameters);
      // visualization      
      if (y == EtaBinLow) {
	canvMany->cd(x);
	gPad->SetLogy();  
	histDeDx->SetMarkerStyle(21);
	histDeDx->SetMarkerSize(0.7);
	histDeDx->Draw("E");
	gausFit->SetLineWidth(2);
	gausFit->Draw("same");
	//
	TF1 gausFit0("gausFit0", "gaus(0)",-1,1);
	TF1 gausFit1("gausFit1", "gaus(0)",-1,1);
	TF1 gausFit2("gausFit2", "gaus(0)",-1,1);
	TF1 gausFit3("gausFit3", "gaus(0)",-1,1);
	gausFit0.SetParameters(parameters[0],parameters[1],parameters[2]);
	gausFit1.SetParameters(parameters[3],parameters[4],parameters[5]);
	gausFit2.SetParameters(parameters[6],parameters[7],parameters[8]);
	gausFit3.SetParameters(parameters[9],parameters[10],parameters[11]);
	//
	gausFit0.SetLineColor(4);
	gausFit1.SetLineColor(2);
	gausFit2.SetLineColor(6);
	gausFit3.SetLineColor(8);
	//
	gausFit0.SetLineWidth(1);
	gausFit1.SetLineWidth(1);
	gausFit2.SetLineWidth(1);
	gausFit3.SetLineWidth(1);
	//
	gausFit0.Draw("same");
	gausFit1.Draw("same"); 
	gausFit2.Draw("same"); 
	gausFit3.Draw("same");
	if (TMath::Abs(pT) < 2) canvMany->Print("canvManyPion.ps");
      }
      
      Double_t yield = gausFit->GetParameter(0)*TMath::Sqrt(2*TMath::Pi())*gausFit->GetParameter(2)/histDeDx->GetBinWidth(1); // area of the gaus fit
      delete gausFit;
      //
      if (y == EtaBinLow && yield > 0) histPt->Fill(pT,yield);
      //delete gausYield;
    }
  }
  canvMany->Print("canvManyPion.ps]");

  TCanvas * canvPt = new TCanvas();
  canvPt->cd();
  histPt->Draw();
  
  return histPt;
}


//________________________________________________________________________
TH1D * AliAnalysisTaskChargedHadronSpectra::AnalyseClassicKaon(const TH3F * input, Int_t EtaBinLow, Int_t EtaBinHigh,const  Double_t * AlephParams) {
  //
  // The multiple Gauss fitting is happening here...
  // input histogram is of the form (Pt, eta, difference in dEdx)
  //

  TF1 *foProton = new TF1("foProton", "AliExternalTrackParam::BetheBlochAleph(x/0.93827,[0],[1],[2],[3],[4])",0.05,20); 
  TF1 *foPion = new TF1("foPion", "AliExternalTrackParam::BetheBlochAleph(x/0.13957,[0],[1],[2],[3],[4])",0.05,20);
  TF1 *foElec = new TF1("foElec", "AliExternalTrackParam::BetheBlochAleph(x/0.000511,[0],[1],[2],[3],[4])",0.05,20);
  TF1 *foKaon = new TF1("foKaon", "AliExternalTrackParam::BetheBlochAleph(x/0.493677,[0],[1],[2],[3],[4])",0.05,20);
  TF1 *foMuon = new TF1("foMuon", "AliExternalTrackParam::BetheBlochAleph(x/0.105658,[0],[1],[2],[3],[4])",0.05,20);
  //
  foProton->SetParameters(AlephParams);  
  foPion->SetParameters(AlephParams);
  foElec->SetParameters(AlephParams);  
  foKaon->SetParameters(AlephParams);  
  foMuon->SetParameters(AlephParams);

  
  const Double_t kSigma = 0.055; // expected resolution (roughly)

  const Int_t kPtBins = 2*30/*input->GetXaxis()->GetNbins()*/;
  const Double_t kPtMax  = input->GetXaxis()->GetXmax();
  //const Int_t kEtaBins = input->GetYaxis()->GetNbins();

  TH1D * histPt = new TH1D("histPt", "Pt; pT (Gev); dN",kPtBins,-kPtMax,kPtMax);

  Double_t binsPt[kPtBins+1] = {-2.00, -1.90, -1.80, -1.70, -1.60, -1.50, -1.40, -1.30, -1.20, -1.10,
				-1.00, -0.95, -0.90, -0.85, -0.80, -0.75, -0.70, -0.65, -0.60, -0.55,
				-0.50, -0.45, -0.40, -0.35, -0.30, -0.25, -0.20, -0.15, -0.10, -0.05,
				 0.00,  0.05,  0.10,  0.15,  0.20,  0.25,  0.30,  0.35,  0.40,  0.45,
				 0.50,  0.55,  0.60,  0.65,  0.70,  0.75,  0.80,  0.85,  0.90,  0.95,
				 1.00,  1.10,  1.20,  1.30,  1.40,  1.50,  1.60,  1.70,  1.80,  1.90, 2.00};

  histPt->GetXaxis()->Set(kPtBins, binsPt);

  TCanvas * canvMany = new TCanvas("canvManyKaon","canvManyKaon");
  canvMany->Print("canvManyKaon.ps[");

  for(Int_t x=1; x < kPtBins+1; x++) {
    Double_t pT = input->GetXaxis()->GetBinCenter(x);
    for(Int_t y=1; y < 2/*kEtaBins+1*/; y++) {
      Double_t rapidity =  input->GetYaxis()->GetBinCenter(y);
      Double_t eta = TMath::ASinH(TMath::Sqrt(1 + (0.4936*0.4936)/(pT*pT))*TMath::SinH(rapidity));  
      Double_t theta = 2*TMath::ATan(TMath::Exp(-eta));
      Double_t ptot = TMath::Abs(pT/TMath::Sin(theta));
      TH1D *histDeDx = input->ProjectionZ(Form("dedx_%i_%i",x,y),x,x,y,y + EtaBinHigh - EtaBinLow);
      histDeDx->SetTitle(Form("pt_%f",pT));
      Float_t maxYield = histDeDx->GetEntries()/(TMath::Sqrt(2*TMath::Pi())*0.04/histDeDx->GetBinWidth(1));
      //
      TF1 * gausFit = new TF1("gausFit", "gaus(0) + gaus(3) + gaus(6) + gaus(9)",-1 /*-5*kSigma*/, 1/*5*kSigma*/);
      Double_t parameters[12] = {maxYield*0.1,0,kSigma,
				 maxYield*0.8,(foPion->Eval(ptot)-foKaon->Eval(ptot))/foKaon->Eval(ptot),kSigma*(foPion->Eval(ptot)/foKaon->Eval(ptot)),
				 maxYield*0.01,(foElec->Eval(ptot)-foKaon->Eval(ptot))/foKaon->Eval(ptot),kSigma*(foElec->Eval(ptot)/foKaon->Eval(ptot)),
				 maxYield*0.1,(foProton->Eval(ptot)-foKaon->Eval(ptot))/foKaon->Eval(ptot),kSigma*(foProton->Eval(ptot)/foKaon->Eval(ptot))};
      gausFit->SetParameters(parameters);
      gausFit->SetParLimits(0,0.,maxYield);
      gausFit->SetParLimits(1,-0.15,0.15);
      gausFit->SetParLimits(2,0.04,0.08);
      //
      Double_t pionPosition = (foPion->Eval(ptot)-foKaon->Eval(ptot))/foKaon->Eval(ptot);
      //Double_t elecPosition = (foElec->Eval(ptot)-foKaon->Eval(ptot))/foKaon->Eval(ptot);
      Double_t protonPosition = (foProton->Eval(ptot)-foKaon->Eval(ptot))/foKaon->Eval(ptot);
      //
      if (TMath::Abs(pT) < 0.4) {
	cout << " fixing the parameters " << endl;
	gausFit->SetParameter(3,0.); gausFit->SetParameter(6,0.); gausFit->SetParameter(9,0.);
	gausFit->FixParameter(3,0); gausFit->FixParameter(6,0);	gausFit->FixParameter(9,0);
	//gausFit->SetParLimits(0,0.,maxYield);
	//gausFit->SetParLimits(1,-0.15,0.15);
	//gausFit->SetParLimits(2,0.03,0.18);
      } else {
	//
	gausFit->SetParLimits(0,0.,0.5*maxYield);
	gausFit->SetParLimits(1,-0.05,0.05);
	gausFit->SetParLimits(2,0.04,0.08);
	//
	gausFit->SetParLimits(3,0.,maxYield);
	gausFit->SetParLimits(4, 0.5*pionPosition, 1.5*pionPosition);
	gausFit->SetParLimits(5,0.6*kSigma*(foPion->Eval(ptot)/foKaon->Eval(ptot)),1.4*kSigma*(foPion->Eval(ptot)/foKaon->Eval(ptot)));
	//
	gausFit->FixParameter(6,0);
	gausFit->FixParameter(7,0);
	gausFit->FixParameter(8,0);
	/*
	  gausFit->SetParLimits(6,0.,0.2*maxYield);
	  gausFit->SetParLimits(7, 0.9*elecPosition, 1.1*elecPosition);
	  gausFit->SetParLimits(8,0.6*kSigma*(foElec->Eval(ptot)/foKaon->Eval(ptot)),1.4*kSigma*(foElec->Eval(ptot)/foKaon->Eval(ptot)));
	*/
	//
	gausFit->SetParLimits(9,0.,0.2*maxYield);
	gausFit->SetParLimits(10, 0.9*protonPosition, 1.1*protonPosition);
	gausFit->SetParLimits(11,0.6*kSigma*(foProton->Eval(ptot)/foKaon->Eval(ptot)),1.4*kSigma*(foProton->Eval(ptot)/foKaon->Eval(ptot)));
      }
      histDeDx->Fit("gausFit","QNR");
      gausFit->GetParameters(parameters);
      // visualization      
      if (y == EtaBinLow) {
	canvMany->cd(x);
	//gPad->SetLogy();  
	histDeDx->SetMarkerStyle(21);
	histDeDx->SetMarkerSize(0.7);
	histDeDx->Draw("E");
	gausFit->SetLineWidth(2);
	gausFit->Draw("same");
	//
	TF1 gausFit0("gausFit0", "gaus(0)",-1,1);
	TF1 gausFit1("gausFit1", "gaus(0)",-1,1);
	TF1 gausFit2("gausFit2", "gaus(0)",-1,1);
	TF1 gausFit3("gausFit3", "gaus(0)",-1,1);
	gausFit0.SetParameters(parameters[0],parameters[1],parameters[2]);
	gausFit1.SetParameters(parameters[3],parameters[4],parameters[5]);
	gausFit2.SetParameters(parameters[6],parameters[7],parameters[8]);
	gausFit3.SetParameters(parameters[9],parameters[10],parameters[11]);
	//
	gausFit0.SetLineColor(4);
	gausFit1.SetLineColor(2);
	gausFit2.SetLineColor(6);
	gausFit3.SetLineColor(8);
	//
	gausFit0.SetLineWidth(1);
	gausFit1.SetLineWidth(1);
	gausFit2.SetLineWidth(1);
	gausFit3.SetLineWidth(1);
	//
	gausFit0.Draw("same");
	gausFit1.Draw("same"); 
	gausFit2.Draw("same"); 
	gausFit3.Draw("same");
	if (TMath::Abs(pT) < 2) canvMany->Print("canvManyKaon.ps");
	//
	// propaganda slots for the paper
	//
	if (pT > 0.374 && pT < 0.376) {
	  TFile f1("slice1.root","RECREATE");
	  canvMany->Write();
	  f1.Close();
	}
	if (pT > 0.524 && pT < 0.526) {
	  TFile f1("slice2.root","RECREATE");
	  canvMany->Write();
	  f1.Close();	  
	}
	if (pT > 0.674 && pT < 0.676) {
	  TFile f1("slice3.root","RECREATE");
	  canvMany->Write();
	  f1.Close();
	}
	if (pT > 0.624 && pT < 0.626) {
	  TFile f1("slice4.root","RECREATE");
	  canvMany->Write();
	  f1.Close();
	}
      }
      
      Double_t yield = gausFit->GetParameter(0)*TMath::Sqrt(2*TMath::Pi())*gausFit->GetParameter(2)/histDeDx->GetBinWidth(1); // area of the gaus fit
      cout << pT << " res.: " << parameters[2] << " " << " yield: " << yield << endl;
      delete gausFit;
      //
      if (y == EtaBinLow && yield > 0) histPt->Fill(pT,yield);
      //delete gausYield;
    }
  }
  canvMany->Print("canvManyKaon.ps]");

  TCanvas * canvPt = new TCanvas();
  canvPt->cd();
  histPt->Draw();
  
  return histPt;
}


//________________________________________________________________________
void AliAnalysisTaskChargedHadronSpectra::Postprocess(const TList * ListOfHistogramsMC,const  TList * ListOfHistogramsData, const Char_t *filename) {

  // Define ranges
  const Int_t kMultLow  = 1;
  const Int_t kMultHigh = 10;

  const Int_t kEtaHigh  = 4;

  // Extract histograms for the MC efficiency study
  TH3F* histPtMCKaonEff = ( TH3F *)ListOfHistogramsMC->FindObject("HistPtMCKaon");
  TH3F* histPtMCProtonEff = ( TH3F *)ListOfHistogramsMC->FindObject("HistPtMCProton");
  TH3F* histPtMCPionEff = ( TH3F *)ListOfHistogramsMC->FindObject("HistPtMCPion");
  
  TH3F* histPtEtaKaonEff = ( TH3F *)ListOfHistogramsMC->FindObject("HistPtEtaKaon");
  TH3F* histPtEtaProtonEff = ( TH3F *)ListOfHistogramsMC->FindObject("HistPtEtaProton");
  TH3F* histPtEtaPionEff = ( TH3F *)ListOfHistogramsMC->FindObject("HistPtEtaPion");
  
  TH3F* histEffProton = ( TH3F *)ListOfHistogramsMC->FindObject("HistEffProton");
  TH3F* histEffKaon = ( TH3F *)ListOfHistogramsMC->FindObject("HistEffKaon");
  TH3F* histEffPion = ( TH3F *)ListOfHistogramsMC->FindObject("HistEffPion");

  // Extract data for the final spectra  
  TH3F* histPtEtaKaon = ( TH3F *)ListOfHistogramsData->FindObject("HistPtEtaKaon");
  TH3F* histPtEtaProton = ( TH3F *)ListOfHistogramsData->FindObject("HistPtEtaProton");
  TH3F* histPtEtaPion = ( TH3F *)ListOfHistogramsData->FindObject("HistPtEtaPion");

  // "MC" of the real data for consistency checks
  TH3F* histPtMCKaonData = ( TH3F *)ListOfHistogramsData->FindObject("HistPtMCKaon");
  TH3F* histPtMCProtonData = ( TH3F *)ListOfHistogramsData->FindObject("HistPtMCProton");
  TH3F* histPtMCPionData = ( TH3F *)ListOfHistogramsData->FindObject("HistPtMCPion");


  TFile spectraFile(filename,"RECREATE");

  // 1. Protons  
  for(Int_t iEta = 1; iEta < kEtaHigh+1; iEta++) {
    TH1D *protonRaw     = histPtEtaProtonEff->ProjectionZ(Form("ProtonRawTot_%i",iEta),kMultLow,kMultHigh,iEta,iEta); // red crosses
    protonRaw->Sumw2();
    TH1D *protonMC      = histPtMCProtonEff->ProjectionZ(Form("McProton_%i",iEta),kMultLow,kMultHigh,iEta,iEta); // solid black line
    protonMC->Sumw2();
    TH1D *protonRawPrim = histEffProton->ProjectionZ(Form("ProtonRawPrim_%i",iEta),1,1,iEta,iEta); // thin red line
    protonRawPrim->Sumw2();
    TH1D *protonSecReac = histEffProton->ProjectionZ(Form("ProtonSecReac_%i",iEta),2,2,iEta,iEta); // blue line
    protonSecReac->Sumw2();
    TH1D *protonSecWeak = histEffProton->ProjectionZ(Form("ProtonSecWeak_%i",iEta),4,4,iEta,iEta); // green line
    protonSecWeak->Sumw2();
    TH1D *protonMis     = histEffProton->ProjectionZ(Form("ProtonMis_%i",iEta),3,3,iEta,iEta); // thin black line
    protonMis->Sumw2();
    protonSecReac->Add(protonSecWeak,-1);
    //
    // 1. total (physical) efficiency
    //
    TH1D *protonEffTot = new TH1D(*protonRaw);
    protonEffTot->Divide(protonRaw,protonMC);
    protonEffTot->SetName(Form("ProtonEffTot_%i",iEta));
    protonEffTot->Write();
    //
    // 2. contributions from secondaries created by interactions
    //
    TH1D *protonEffSecReac = new TH1D(*protonRaw);
    TH1D *dummy = new TH1D(*protonRaw);
    dummy->Add(protonSecReac,-1);
    protonEffSecReac->Divide(protonRaw,dummy);
    protonEffSecReac->SetName(Form("ProtonEffSecReac_%i",iEta));
    protonEffSecReac->Write();
    delete dummy;
    //
    // 3. contributions from secondaries from weak decays
    //
    TH1D *protonEffSecWeak = new TH1D(*protonRaw);
    TH1D *dummy2 = new TH1D(*protonRaw);
    dummy2->Add(protonSecWeak,-1);
    protonEffSecWeak->Divide(protonRaw,dummy2);
    protonEffSecWeak->SetName(Form("ProtonEffSecWeak_%i",iEta));
    protonEffSecWeak->Write();
    delete dummy2;
    //
    // 4. decay (Kaon,pion) and detector absorption (energy loss and multiple scattering)
    //
    TH1D *protonEffAbs = new TH1D(*protonRaw);
    protonEffAbs->Divide(protonRawPrim,protonMC);
    protonEffAbs->SetName(Form("ProtonEffAbs_%i",iEta));
    protonEffAbs->Write();
    //
    // 5. contamination
    //
    TH1D *protonEffMis = new TH1D(*protonMis);
    protonEffMis->Divide(protonMis,protonRaw);
    protonEffMis->SetName(Form("ProtonContam_%i",iEta));
    protonEffMis->Write();

    //
    // PROCESS THE REAL DATA FROM HERE ON
    //
    TH1D *histPtProtonEta = histPtEtaProton->ProjectionZ(Form("HistPtProtonEtaRaw_%i",iEta),kMultLow,kMultHigh,iEta,iEta);
    histPtProtonEta->Write();
    //
    // Apply corrections as you would like them to have...
    //
    TH1D *protonFinal = new TH1D(*histPtProtonEta);
    protonFinal->SetName(Form("HistPtProtonEtaFinal_%i",iEta));
    protonFinal->Sumw2();
    protonFinal->Divide(protonEffTot); // enable or disable specific corrections here...
    protonFinal->Write();
    //
    // if available save also the MC truth for consistency checks
    //
    TH1D *histPtProtonEtaMcData = histPtMCProtonData->ProjectionZ(Form("HistPtProtonEtaMcData_%i",iEta),kMultLow,kMultHigh,iEta,iEta);
    histPtProtonEtaMcData->Write();
  }
  
  // 2. Kaons
  for(Int_t iEta = 1; iEta < kEtaHigh+1; iEta++) {
    TH1D *kaonRaw     = histPtEtaKaonEff->ProjectionZ(Form("KaonRawTot_%i",iEta),kMultLow,kMultHigh,iEta,iEta); // red crosses
    kaonRaw->Sumw2();
    TH1D *kaonMC      = histPtMCKaonEff->ProjectionZ(Form("McKaon_%i",iEta),kMultLow,kMultHigh,iEta,iEta); // solid black line
    kaonMC->Sumw2();
    TH1D *kaonRawPrim = histEffKaon->ProjectionZ(Form("KaonRawPrim_%i",iEta),1,1,iEta,iEta); // thin red line
    kaonRawPrim->Sumw2();
    TH1D *kaonSecReac = histEffKaon->ProjectionZ(Form("KaonSecReac_%i",iEta),2,2,iEta,iEta); // blue line
    kaonSecReac->Sumw2();
    TH1D *kaonSecWeak = histEffKaon->ProjectionZ(Form("KaonSecWeak_%i",iEta),4,4,iEta,iEta); // green line
    kaonSecWeak->Sumw2();
    TH1D *kaonMis     = histEffKaon->ProjectionZ(Form("KaonMis_%i",iEta),3,3,iEta,iEta); // thin black line
    kaonMis->Sumw2();
    kaonSecReac->Add(kaonSecWeak,-1);
    //
    // 1. total (physical) efficiency
    //
    TH1D *kaonEffTot = new TH1D(*kaonRaw);
    kaonEffTot->Divide(kaonRaw,kaonMC);
    kaonEffTot->SetName(Form("KaonEffTot_%i",iEta));
    kaonEffTot->Write();
    //
    // 2. contributions from secondaries created by interactions
    //
    TH1D *kaonEffSecReac = new TH1D(*kaonRaw);
    TH1D *dummy = new TH1D(*kaonRaw);
    dummy->Add(kaonSecReac,-1);
    kaonEffSecReac->Divide(kaonRaw,dummy);
    kaonEffSecReac->SetName(Form("KaonEffSecReac_%i",iEta));
    kaonEffSecReac->Write();
    delete dummy;
    //
    // 3. contributions from secondaries from weak decays
    //
    TH1D *kaonEffSecWeak = new TH1D(*kaonRaw);
    TH1D *dummy2 = new TH1D(*kaonRaw);
    dummy2->Add(kaonSecWeak,-1);
    kaonEffSecWeak->Divide(kaonRaw,dummy2);
    kaonEffSecWeak->SetName(Form("KaonEffSecWeak_%i",iEta));
    kaonEffSecWeak->Write();
    delete dummy2;
    //
    // 4. decay (Kaon,pion) and detector absorption (energy loss and multiple scattering)
    //
    TH1D *kaonEffAbs = new TH1D(*kaonRaw);
    kaonEffAbs->Divide(kaonRawPrim,kaonMC);
    kaonEffAbs->SetName(Form("KaonEffAbs_%i",iEta));
    kaonEffAbs->Write();
    //
    // 5. contamination
    //
    TH1D *kaonEffMis = new TH1D(*kaonMis);
    kaonEffMis->Divide(kaonMis,kaonRaw);
    kaonEffMis->SetName(Form("KaonContam_%i",iEta));
    kaonEffMis->Write();
    
    //
    // PROCESS THE REAL DATA FROM HERE ON
    //
    TH1D *histPtKaonEta = histPtEtaKaon->ProjectionZ(Form("HistPtKaonEtaRaw_%i",iEta),kMultLow,kMultHigh,iEta,iEta);
    histPtKaonEta->Write();
    //
    // Apply corrections as you would like them to have...
    //
    TH1D *kaonFinal = new TH1D(*histPtKaonEta);
    kaonFinal->SetName(Form("HistPtKaonEtaFinal_%i",iEta));
    kaonFinal->Sumw2();
    kaonFinal->Divide(kaonEffTot); // enable or disable specific corrections here...
    kaonFinal->Write();
    //
    // if available save also the MC truth for consistency checks
    //
    TH1D *histPtKaonEtaMcData = histPtMCKaonData->ProjectionZ(Form("HistPtKaonEtaMcData_%i",iEta),kMultLow,kMultHigh,iEta,iEta);
    histPtKaonEtaMcData->Write();
  }


  // 3. Pions
  for(Int_t iEta = 1; iEta < kEtaHigh+1; iEta++) {
    TH1D *pionRaw     = histPtEtaPionEff->ProjectionZ(Form("PionRawTot_%i",iEta),kMultLow,kMultHigh,iEta,iEta); // red crosses
    pionRaw->Sumw2();
    TH1D *pionMC      = histPtMCPionEff->ProjectionZ(Form("McPion_%i",iEta),kMultLow,kMultHigh,iEta,iEta); // solid black line
    pionMC->Sumw2();
    TH1D *pionRawPrim = histEffPion->ProjectionZ(Form("PionRawPrim_%i",iEta),1,1,iEta,iEta); // thin red line
    pionRawPrim->Sumw2();
    TH1D *pionSecReac = histEffPion->ProjectionZ(Form("PionSecReac_%i",iEta),2,2,iEta,iEta); // blue line
    pionSecReac->Sumw2();
    TH1D *pionSecWeak = histEffPion->ProjectionZ(Form("PionSecWeak_%i",iEta),4,4,iEta,iEta); // green line
    pionSecWeak->Sumw2();
    TH1D *pionMis     = histEffPion->ProjectionZ(Form("PionMis_%i",iEta),3,3,iEta,iEta); // thin black line
    pionMis->Sumw2();
    pionSecReac->Add(pionSecWeak,-1);
    //
    // 1. total (physical) efficiency
    //
    TH1D *pionEffTot = new TH1D(*pionRaw);
    pionEffTot->Divide(pionRaw,pionMC);
    pionEffTot->SetName(Form("PionEffTot_%i",iEta));
    pionEffTot->Write();
    //
    // 2. contributions from secondaries created by interactions
    //
    TH1D *pionEffSecReac = new TH1D(*pionRaw);
    TH1D *dummy = new TH1D(*pionRaw);
    dummy->Add(pionSecReac,-1);
    pionEffSecReac->Divide(pionRaw,dummy);
    pionEffSecReac->SetName(Form("PionEffSecReac_%i",iEta));
    pionEffSecReac->Write();
    delete dummy;
    //
    // 3. contributions from secondaries from weak decays
    //
    TH1D *pionEffSecWeak = new TH1D(*pionRaw);
    TH1D *dummy2 = new TH1D(*pionRaw);
    dummy2->Add(pionSecWeak,-1);
    pionEffSecWeak->Divide(pionRaw,dummy2);
    pionEffSecWeak->SetName(Form("PionEffSecWeak_%i",iEta));
    pionEffSecWeak->Write();
    delete dummy2;
    //
    // 4. decay (Kaon,pion) and detector absorption (energy loss and multiple scattering)
    //
    TH1D *pionEffAbs = new TH1D(*pionRaw);
    pionEffAbs->Divide(pionRawPrim,pionMC);
    pionEffAbs->SetName(Form("PionEffAbs_%i",iEta));
    pionEffAbs->Write();
    //
    // 5. contamination
    //
    TH1D *pionEffMis = new TH1D(*pionMis);
    pionEffMis->Divide(pionMis,pionRaw);
    pionEffMis->SetName(Form("PionContam_%i",iEta));
    pionEffMis->Write();
    
    //
    // PROCESS THE REAL DATA FROM HERE ON
    //
    TH1D *histPtPionEta = histPtEtaPion->ProjectionZ(Form("HistPtPionEtaRaw_%i",iEta),kMultLow,kMultHigh,iEta,iEta);
    histPtPionEta->Write();
    //
    // Apply corrections as you would like them to have...
    //
    TH1D *pionFinal = new TH1D(*histPtPionEta);
    pionFinal->SetName(Form("HistPtPionEtaFinal_%i",iEta));
    pionFinal->Sumw2();
    pionFinal->Divide(pionEffTot); // enable or disable specific corrections here...
    pionFinal->Write();
    //
    // if available save also the MC truth for consistency checks
    //
    TH1D *histPtPionEtaMcData = histPtMCPionData->ProjectionZ(Form("HistPtPionEtaMcData_%i",iEta),kMultLow,kMultHigh,iEta,iEta);
    histPtPionEtaMcData->Write();

  }

  //
  // Close the file after everthing is done
  //
  spectraFile.Close();

}




