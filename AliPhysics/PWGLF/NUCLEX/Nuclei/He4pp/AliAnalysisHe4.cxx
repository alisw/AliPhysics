
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

#include "AliAnalysisHe4.h"

#include "AliAnalysisUtils.h"
#include "TProfile2D.h"


ClassImp(AliAnalysisHe4)

//________________________________________________________________________
AliAnalysisHe4::AliAnalysisHe4() 
: AliAnalysisTaskSE("TaskChargedHadron"), 
  fESD(0), 
  fListHist(0), 
  fESDtrackCuts(0),
  fESDTrackCutsMult(0),
  fESDpid(0),
  fMCtrue(0),
  fAlephParameters(),
  fUtils(0),
  fHistRealTracks(0),
  fHistMCparticles(0),
  fHistPidQA(0),
  fHistTofQA(0),
  fHistMult(0),
  fHistMomCorr(0),
  fHistEtaPtGen(0),
  fHistVertex(0),
  fHistVertexRes(0),
  fHistVertexResTracks(0),
  fHistDeDx(0),

  fHistHe3(0),
  fHistHe4(0),
  fHistAntiHe3(0),
  fHistAntiHe4(0),
  fHistTPCsigHe3(0),
  fHistTPCsigHe4(0)
{
  // default Constructor
}


//________________________________________________________________________
AliAnalysisHe4::AliAnalysisHe4(const char *name) 
  : AliAnalysisTaskSE(name), fESD(0), fListHist(0), fESDtrackCuts(0),fESDTrackCutsMult(0),fESDpid(0),
    fMCtrue(0),
    fAlephParameters(),
    fUtils(0),
    fHistRealTracks(0),
    fHistMCparticles(0),
    fHistPidQA(0),
    fHistTofQA(0),
    fHistMult(0),
    fHistMomCorr(0),
    fHistEtaPtGen(0),
    fHistVertex(0),
    fHistVertexRes(0),
    fHistVertexResTracks(0),
    fHistDeDx(0),
    fHistHe3(0),
    fHistHe4(0),  
    fHistAntiHe3(0),
    fHistAntiHe4(0),
    fHistTPCsigHe3(0),
    fHistTPCsigHe4(0)
    
{
  //
  // standard constructur which should be used
  //
  // He3
  fAlephParameters[0] = 0.0283086;
  fAlephParameters[1] = 2.63394e+01;
  fAlephParameters[2] = 5.04114e-11;
  fAlephParameters[3] = 2.12543e+00;
  fAlephParameters[4] = 4.88663e+00;
  
  // He4
//   fAlephParameters[0]  = 1.69155;
//   fAlephParameters[1]  = 27.4992;
//   fAlephParameters[2]  = 4.00313e-15;
//   fAlephParameters[3]  = 2.48485;
//   fAlephParameters[4]  = 8.31768;
//     
    
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
void AliAnalysisHe4::Initialize()
{
  //
  // updating parameters in case of changes
  //
  fESDtrackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011(kFALSE,kTRUE);
  fESDtrackCuts->SetMaxDCAToVertexXY(3);
  fESDtrackCuts->SetMaxDCAToVertexZ(2);
  fESDtrackCuts->SetEtaRange(-0.9,0.9);
  //
  fESDtrackCuts->SetAcceptKinkDaughters(kFALSE);
  fESDtrackCuts->SetMinNClustersTPC(70);
  fESDtrackCuts->SetMaxChi2PerClusterTPC(6);
  fESDtrackCuts->SetRequireTPCRefit(kTRUE);
//   fESDtrackCuts->SetRequireITSRefit(kTRUE);
  fESDtrackCuts->SetMinNClustersITS(1);
  


}


//________________________________________________________________________
void AliAnalysisHe4::UserCreateOutputObjects() 
{

  //Create analysis utils for event selection and pileup rejection
  fUtils = new AliAnalysisUtils();
  fUtils->SetCutOnZVertexSPD(kFALSE);


  // Create histograms
  // Called once
  fListHist = new TList();
  fListHist->SetOwner(kTRUE);
  
  // binning 
  const Int_t kPtBins = 28;
  // const Int_t kMultBins = 11;
  const Int_t kDcaBins = 38;


  // sort pT-bins ..
  Double_t binsPt[kPtBins+1] = {0., 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.2, 4.4, 4.6, 4.8, 5.0, 5.2, 5.4, 5.6};
 
  
  // provide a reasonable dca/binning
  Double_t binsDca[kDcaBins+1];
  for(Int_t i = 0; i< kDcaBins+1; i++) {
    binsDca[i] = 1.5*(1./TMath::Exp(-TMath::Abs((i - kDcaBins/2.))/(kDcaBins/2)) -1);
    if ((i - kDcaBins/2.) < 0) binsDca[i] *= -1;
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
  //                           0      2,  3,      4,     5,    6,    7,       8
//   Int_t    binsHistReal[8] = { 4,   kPtBins,  2,      12,   50,    2,  100, kDcaBins};
//   Double_t xminHistReal[8] = {-0.5,    0, -2,    -0.6,   -5,- 0.5, -2.5,       -3};
//   Double_t xmaxHistReal[8] = {3.5,      3,  2,     0.6,    5,  1.5,  2.5,        3};
//   fHistRealTracks = new THnSparseF("fHistRealTracks","real tracks",8,binsHistReal,xminHistReal,xmaxHistReal);
//   //
//   fHistRealTracks->GetAxis(1)->Set(kPtBins, binsPt);
//   fHistRealTracks->GetAxis(7)->Set(kDcaBins, binsDca);
//   
//   fHistRealTracks->GetAxis(0)->SetTitle("assumed particle");
//   fHistRealTracks->GetAxis(1)->SetTitle("p_{T} (GeV/c)");
//   fHistRealTracks->GetAxis(2)->SetTitle("sign");
//   fHistRealTracks->GetAxis(3)->SetTitle("y");
//   fHistRealTracks->GetAxis(4)->SetTitle("TPC dEdx");
//   fHistRealTracks->GetAxis(5)->SetTitle("has valid TOF");
//   fHistRealTracks->GetAxis(6)->SetTitle("nsigma TOF");
//   fHistRealTracks->GetAxis(7)->SetTitle("DCA_{xy}");
  
//   fListHist->Add(fHistRealTracks);

  
//   fHistPidQA = new TH3F("fHistPidQA","PID QA",500,0.1,10,1000,0,1000,2,-2,2);
//   BinLogAxis(fHistPidQA);
//   fListHist->Add(fHistPidQA);
// 
//   fHistTofQA = new TH2F("fHistTofQA","TOF-QA",200,-4.,4.,400,0,1.1);
//   fListHist->Add(fHistTofQA);

  
  
  
// MC hist
//   //                            0,            1,           2,  3,      4,            5,    6,    7,        8,    9
//   Int_t    binsHistMC[10] = {   3,    kMultBins,     kPtBins,  2,    12      ,  50,    2,  100, kDcaBins,    6};
//   Double_t xminHistMC[10] = {-0.5,         -0.5,           0, -2,    -0.6 ,  -5,- 0.5, -2.5,       -3, -0.5};
//   Double_t xmaxHistMC[10] = { 2.5,         10.5,           3,  2,    0.6,   5,  1.5,  2.5,        3,  5.5};
//   //
//   // different binning for CODE axis, if we want to save motherPDG
//   //
//   fHistMCparticles = new THnSparseF("fHistMCparticles","MC histogram",10,binsHistMC,xminHistMC,xmaxHistMC);
//   fHistMCparticles->GetAxis(2)->Set(kPtBins, binsPt);
//   fHistMCparticles->GetAxis(8)->Set(kDcaBins, binsDca);
//   fListHist->Add(fHistMCparticles);
  //
  

  fHistMult = new TH2F("fHistMult", "control histogram to count number of events", 502, -2.5, 499.5,4,-0.5,3.5);
  fListHist->Add(fHistMult);

  //
//   fHistMomCorr = new TH2F("fHistMomCorr","momentum correction; pT_{MCtrue}; rel. difference",200,0,3,200,-0.2,0.2);
//   fListHist->Add(fHistMomCorr);
  //
//   fHistEtaPtGen = new TH3F("fHistPtEtaGen","rapidity correction; p_{T} (GeV/c); rapidity y; pid",200,0,12, 200,-1.2,1.2, 3,-0.5,2.5);
//   fListHist->Add(fHistEtaPtGen);
  //
  fHistVertex = new TH2F("fHistVertex","vertex position; V_{XY} (cm); V_{Z}",2000,-1.,1.,400,-20.,20.);
  fListHist->Add(fHistVertex);
  //
//   fHistVertexRes = new TH1F("fHistVertexRes","vertex position difference MC truth and rec; V_{Z,rec} - V_{Z,truth} (cm)",2000,-0.5,0.5);
//   fListHist->Add(fHistVertexRes);
  //
//   fHistVertexResTracks = new TH2F("fHistVertexResTracks","vertex res vs no of tracks; no. of ch. primaries |#eta| < 0.9; V_{Z,rec} - V_{Z,truth} (cm)",50,-0.5,499.5,200,-0.05,0.05);
//   fListHist->Add(fHistVertexResTracks);

  
  
  fHistDeDx = new TH2F("fHistDeDx", "dE/dx", 100, 0.1, 6.0, 1000, 0.0, 1000);
  fHistDeDx->GetYaxis()->SetTitle("TPC dE/dx signal (a.u.)");
  fHistDeDx->GetXaxis()->SetTitle("#frac{p}{z} (GeV/c)");
  BinLogAxis(fHistDeDx);
  fListHist->Add(fHistDeDx);

  
//    const Int_t ndims = 2;
//    Int_t bins[ndims] = {100, 1000};
//    Double_t xmin[ndims] = {0.1, 0.,};
//    Double_t xmax[ndims] = {6. ,1000.};
//    THnSparse* fHistDeDx = new THnSparseD("fHistDeDx", "dE/dx", ndims, bins, xmin, xmax);
//     fListHist->Add(fHistDeDx);

  
  fHistHe3 = new TH1F("fHistHe3", "pT hist he3 canidates", kPtBins,binsPt);
  fHistHe3->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  fListHist->Add(fHistHe3);
  
  fHistHe4 = new TH1F("fHistHe4", "pT hist he4 canidates",kPtBins,binsPt);
  fHistHe4->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  fListHist->Add(fHistHe4);

  fHistAntiHe3 = new TH1F("fHistAntiHe3", "pT hist anti he3 canidates", kPtBins,binsPt);
  fHistAntiHe3->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  fListHist->Add(fHistAntiHe3);
  
  fHistAntiHe4 = new TH1F("fHistAntiHe4", "pT hist anit he4 canidates",kPtBins,binsPt);
  fHistAntiHe4->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  fListHist->Add(fHistAntiHe4);
  
  

  fHistTPCsigHe3 = new TProfile2D("fHistTPCsigHe3", "hist to check cut on TPC sigma for he3", 100, 0.1, 6.0, 1000, 0.0, 1000.);
  BinLogAxis(fHistTPCsigHe3);
  fListHist->Add(fHistTPCsigHe3);
  
  fHistTPCsigHe4 = new TProfile2D("fHistTPCsigHe4", "hist to check cut on TPC sigma for he4", 100, 0.1, 6.0, 1000, 0.0, 1000.);
  BinLogAxis(fHistTPCsigHe4);
  fListHist->Add(fHistTPCsigHe4);  
  
  PostData(1, fListHist);

}

//________________________________________________________________________
void AliAnalysisHe4::UserExec(Option_t *) 
{

  // main event loop
  if (!fESDpid) fESDpid = ((AliESDInputHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->GetESDpid();
  
  if (!fESDpid){
    fESDpid = new AliESDpid();     
    fESDpid->GetTPCResponse().SetBetheBlochParameters(fAlephParameters[0],fAlephParameters[1], fAlephParameters[2], fAlephParameters[3],fAlephParameters[4]);
  }
  
  
  AliLog::SetGlobalLogLevel(AliLog::kError);
  
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
  
  fESD = dynamic_cast<AliESDEvent*>( InputEvent() );
  if (!fESD) {
    Printf("ERROR: fESD not available");
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
  if (vertex && isVtxOk) fHistVertex->Fill(TMath::Sqrt(vertex->GetX()*vertex->GetX() + vertex->GetY()*vertex->GetY()),vertex->GetZ());

  
  //
  // count events after physics selection and after vertex selection
  //
  fHistMult->Fill(fESDtrackCuts->CountAcceptedTracks(fESD), 0);
  fHistMult->Fill(fESD->GetNumberOfTracks(), 1);

  
  
  //
  //***************************************************
  // track loop
  //***************************************************
  //const Double_t kNsigmaCut = 3;
  //const Double_t k2sigmaCorr = 1/(0.5*(TMath::Erf(kNsigmaCut/sqrt(2))-TMath::Erf(-kNsigmaCut/sqrt(2))))/*1/0.9545*/;
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
    
    
    UInt_t status = track->GetStatus();
    //
    Bool_t hasTOFout  = status&AliESDtrack::kTOFout; 
    // Bool_t hasTOFtime = status&AliESDtrack::kTIME;
    
    Double_t length = track->GetIntegratedLength(); 
    Bool_t hasTOF     = kFALSE;
    if (hasTOFout && length < 350.) hasTOF = kTRUE;

    
    /*
    //
    // calculate rapidities and kinematics
    //
    Double_t massHe3=AliPID::ParticleMass(AliPID::kHe3);
    Double_t massHe4=AliPID::ParticleMass(AliPID::kAlpha);

    
    Double_t energyTri = TMath::Sqrt(track->GetP()*2*track->GetP()*2 + massTri*massTri);
    Double_t energyHe3 = TMath::Sqrt(track->GetP()*2*track->GetP()*2 + massHe3*massHe3);
    Double_t energyHe4 = TMath::Sqrt(track->GetP()*2*track->GetP()*2 + massHe4*massHe4);
    
    Double_t pvec[3];
    track->GetPxPyPz(pvec);    
    
    Double_t rapTri = 0.5*TMath::Log((energyTri + pvec[2]*2)/(energyTri - pvec[2]*2));
    Double_t rapHe3 = 0.5*TMath::Log((energyHe3 + pvec[2]*2)/(energyHe3 - pvec[2]*2));
    Double_t rapHe4 = 0.5*TMath::Log((energyHe4 + pvec[2]*2)/(energyHe4 - pvec[2]*2));
    
    */
    
    
    
    // 3. make the PID
    Double_t sign = track->GetSign();   
    Double_t tpcSignal = track->GetTPCsignal();
    
//     if ( tpcSignal > 1000 || tpcSignal < 120) continue; // do i need that??? very strict cut
//     Double_t tpcVec[3] = {ptot,tpcSignal, sign};
//     fHistDeDx->Fill(tpcVec);
    fHistDeDx->Fill(ptot,tpcSignal, sign);

    
    // fill two hist with he3 and he4 
    // exclude everything with sigma triton < 3
    // he3 if |sigma| < 3 
    if (fESDpid->NumberOfSigmasTPC(track, AliPID::kTriton) < 3) continue;
    
    Bool_t isHe3=kFALSE;
    Bool_t isHe4=kFALSE;
    
    if (TMath::Abs(fESDpid->NumberOfSigmasTPC(track, AliPID::kHe3)) < 3) isHe3=kTRUE;
//     if (fESDpid->NumberOfSigmasTPC(track, AliPID::kHe3) > 3) isHe4=kTRUE;
    if (TMath::Abs(fESDpid->NumberOfSigmasTPC(track, AliPID::kAlpha)) < 3) isHe4=kTRUE;

    
    if(isHe3){
//       Double_t vecTPCHistHe3[3] = {pT,   tpcSignal, fESDpid->NumberOfSigmasTPC(track, AliPID::kHe3)};
      fHistTPCsigHe3->Fill(pT,   tpcSignal, fESDpid->NumberOfSigmasTPC(track, AliPID::kHe3));
      if (sign > 0) fHistHe3->Fill(pT);  
      if (sign < 0) fHistAntiHe3->Fill(pT);  
    }

    if(isHe4){
//       Double_t vecTPCHistHe4[3] = {pT,   tpcSignal, fESDpid->NumberOfSigmasTPC(track, AliPID::kAlpha)};
      fHistTPCsigHe4->Fill(pT,   tpcSignal, fESDpid->NumberOfSigmasTPC(track, AliPID::kAlpha));
      if (sign > 0) fHistHe4->Fill(pT);  
      if (sign < 0) fHistAntiHe4->Fill(pT);  
    }
    
    

  } // end of track loop  
  // Post output data  
  PostData(1, fListHist);

}      


//________________________________________________________________________
void AliAnalysisHe4::Terminate(Option_t *) 
{
  // Draw result to the screen
  // Called once at the end of the query
  Printf("*** CONSTRUCTOR CALLED ****");
}


//________________________________________________________________________
void AliAnalysisHe4::BinLogAxis(const TH1 *h) {
  //
  // Method for the correct logarithmic binning of histograms
  //
  TAxis *axis = const_cast<TAxis*>(h->GetXaxis());
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