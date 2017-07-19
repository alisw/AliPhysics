#include <Riostream.h>
#include "TChain.h"
#include "TTree.h"
#include "TF1.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "THn.h"
#include "TCanvas.h"
#include "AliESDtrackCuts.h"
#include "AliESDpid.h"
#include "AliESDVertex.h"
#include "TFile.h"

#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliCentrality.h"
#include "AliESDEvent.h"
#include "AliESDInputHandler.h"

#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliStack.h"

#include "AliAnalysisTaskAntiHe4.h"

using std::cout;
using std::endl;

///////////////////////////////////////////////////////////
//                                                       //
// Analysis for the observation of anti-alpha particles. //
//                                                       //
///////////////////////////////////////////////////////////


ClassImp(AliAnalysisTaskAntiHe4)

//________________________________________________________________________
AliAnalysisTaskAntiHe4::AliAnalysisTaskAntiHe4()
  : AliAnalysisTaskSE(), 
  fEventHandler(0),
  fESD(0), 
  fHistCentralityClass10(0),
  fHistCentralityPercentile(0),
  fHistTriggerStat(0),
  fHistTriggerStatAfterEventSelection(0),
  fHistDEDx(0), 
  fHistTOF3D(0),
  fHistAlpha(0),
  fHistAlphaSignal(0),
  fGraphAlphaSignal(0),
  fNCounter(0),
  fHistDeDx(0),
  fHistDeDxRegion(0),
  fHistDeDxSharp(0),
  fHistTOF2D(0),
  fHistTOFnuclei(0x0),
  fNTriggers(5),
  fBBParametersLightParticles(),
  fBBParametersNuclei(),
  fMCtrue(0),
  fTriggerFired(),
  fESDtrackCuts(0),
  fESDtrackCutsSharp(0),
  fESDpid(0), 
  fAntiAlpha(0),
  fHistHelium4PtGen(0),
  fHistHelium4PtGenPrim(0),
  fHistHelium4PtGenSec(0),
  fHistHelium4PtGenEta(0),
  fHistHelium4PtGenPrimEta(0),
  fHistAntiHelium4PtGen(0),
  fHistAntiHelium4PtGenPrim(0),
  fHistAntiHelium4PtGenSec(0),
  fHistAntiHelium4PtGenEta(0),
  fHistHelium4PtAso(0),
  fHistHelium4PtAsoPrim(0),
  fHistHelium4PtAsoSec(0),
  fHistAntiHelium4PtAso(0),
  fTree(0), 
  fOutputContainer(0),
  fEvnt(0),
  fItrk(0)
{
  // default Constructor
  
  
  // Define input and output slots here
}

//________________________________________________________________________
AliAnalysisTaskAntiHe4::AliAnalysisTaskAntiHe4(const char *name)
  : AliAnalysisTaskSE(name),
    fEventHandler(0), 
    fESD(0),
    fHistCentralityClass10(0),
    fHistCentralityPercentile(0), 
    fHistTriggerStat(0),
    fHistTriggerStatAfterEventSelection(0),
    fHistDEDx(0), 
    fHistTOF3D(0),
    fHistAlpha(0),
    fHistAlphaSignal(0),
    fGraphAlphaSignal(0),
    fNCounter(0),
    fHistDeDx(0),
    fHistDeDxRegion(0),
    fHistDeDxSharp(0),
    fHistTOF2D(0),
    fHistTOFnuclei(0x0),
    fNTriggers(5),
    fBBParametersLightParticles(),
    fBBParametersNuclei(),
    fMCtrue(0),
    fTriggerFired(),
    fESDtrackCuts(0),
    fESDtrackCutsSharp(0),
    fESDpid(0), 
    fAntiAlpha(0),
    fHistHelium4PtGen(0),
    fHistHelium4PtGenPrim(0),
    fHistHelium4PtGenSec(0),
    fHistHelium4PtGenEta(0),
    fHistHelium4PtGenPrimEta(0),
    fHistAntiHelium4PtGen(0),
    fHistAntiHelium4PtGenPrim(0),
    fHistAntiHelium4PtGenSec(0),
    fHistAntiHelium4PtGenEta(0),
    fHistHelium4PtAso(0),
    fHistHelium4PtAsoPrim(0),
    fHistHelium4PtAsoSec(0),
    fHistAntiHelium4PtAso(0),
    fTree(0), 
    fOutputContainer(0),
    fEvnt(0),
    fItrk(0)
{
  // Constructor

  // Define input and output slots here
  // Input slot #0 works with a TChain
  DefineInput(0, TChain::Class());
  // Output slot #0 writes into a TH1 container
  DefineOutput(1, TObjArray::Class());
  DefineOutput(2, TTree::Class());
  //
  // cuts for candidates
  //
  fESDtrackCuts = new AliESDtrackCuts("AliESDtrackCuts","AliESDtrackCuts");
  //
  fESDtrackCuts->SetAcceptKinkDaughters(kFALSE);
  fESDtrackCuts->SetMinNClustersTPC(70);
  fESDtrackCuts->SetMaxChi2PerClusterTPC(6);
  fESDtrackCuts->SetMaxDCAToVertexXY(3);
  fESDtrackCuts->SetMaxDCAToVertexZ(2);
  fESDtrackCuts->SetRequireTPCRefit(kTRUE);
  //fESDtrackCuts->SetRequireITSRefit(kTRUE);
  fESDtrackCuts->SetMinNClustersITS(1);
  fESDtrackCuts->SetEtaRange(-1.0,1.0);
  //
  // cuts for final plots
  //
  fESDtrackCutsSharp  = new AliESDtrackCuts("AliESDtrackCuts","AliESDtrackCuts");
  fESDtrackCutsSharp->SetAcceptKinkDaughters(kFALSE);
  fESDtrackCutsSharp->SetMinNClustersTPC(80);
  fESDtrackCutsSharp->SetMaxChi2PerClusterITS(10);// TO BE INVESTIGATED !!!!!!!!!!!!!!
  fESDtrackCutsSharp->SetMaxChi2PerClusterTPC(5);
  fESDtrackCutsSharp->SetRequireTPCRefit(kTRUE);
  fESDtrackCutsSharp->SetRequireITSRefit(kTRUE);
  fESDtrackCutsSharp->SetMinNClustersITS(2);
  fESDtrackCutsSharp->SetMaxDCAToVertexXY(0.1);
  fESDtrackCutsSharp->SetMaxDCAToVertexZ(0.5);
  fESDtrackCutsSharp->SetEtaRange(-0.8,0.8);

  //ESD Track cuts  from TestFilterRawTask
  /*  fESDtrackCuts = new AliESDtrackCuts("AliESDtrackCuts","AliESDtrackCuts");
  fESDtrackCuts->SetMinNClustersTPC(80);
  fESDtrackCuts->SetMinNClustersITS(2);
  fESDtrackCuts->SetAcceptKinkDaughters(kFALSE);
  fESDtrackCuts->SetMaxChi2PerClusterTPC(5);
  fESDtrackCuts->SetRequireTPCRefit(kTRUE);
  fESDtrackCuts->SetRequireITSRefit(kTRUE);
  fESDtrackCuts->SetEtaRange(-1,1);
  fESDtrackCuts->SetMaxDCAToVertexXY(1);
  fESDtrackCuts->SetMaxDCAToVertexZ(2);
  //test strickter cuts                  
  //fESDtrackCuts->SetMaxDCAToVertexXY(0.1);
  //fESDtrackCuts->SetMaxDCAToVertexZ(0.5);
  //fESDtrackCuts->SetEtaRange(-0.8,0.8);     
  */

  fMCtrue = kTRUE;

}

//________________________________________________________________________
void AliAnalysisTaskAntiHe4::UserCreateOutputObjects()
{
  // Create histograms
  // Called once
  //
  //histogram to count number of events                                 
  //
  fHistCentralityClass10  = new TH1F("fHistCentralityClass10", "centrality with class10", 11, -0.5, 10.5);
  fHistCentralityClass10->GetXaxis()->SetTitle("Centrality");
  fHistCentralityClass10->GetYaxis()->SetTitle("Entries");
  //
  fHistCentralityPercentile  = new TH1F("fHistCentralityPercentile", "centrality with percentile", 101, -0.1, 100.1);
  fHistCentralityPercentile->GetXaxis()->SetTitle("Centrality");
  fHistCentralityPercentile->GetYaxis()->SetTitle("Entries");
  //
  //trigger statitics histogram
  //
  fHistTriggerStat = new TH1F("fHistTriggerStat","Trigger statistics", fNTriggers,-0.5,fNTriggers-0.5);
  const Char_t* aTriggerNames[] = { "kMB", "kCentral", "kSemiCentral", "kEMCEJE", "kEMCEGA" };
  for ( Int_t ii=0; ii < fNTriggers; ii++ )
    fHistTriggerStat->GetXaxis()->SetBinLabel(ii+1, aTriggerNames[ii]);

  fHistTriggerStatAfterEventSelection = new TH1F("fHistTriggerStatAfterEventSelection","Trigger statistics after event selection", fNTriggers,-0.5,fNTriggers-0.5);
  for ( Int_t ii=0; ii < fNTriggers; ii++ )
    fHistTriggerStatAfterEventSelection->GetXaxis()->SetBinLabel(ii+1, aTriggerNames[ii]);

  //dE/dx performance 
  fHistDEDx= new TH3F("fHistDEDx","DEDx",1000,0.01,100,1000,1,2000,2,-2,2);
  BinLogAxis(fHistDEDx, 0);
  BinLogAxis(fHistDEDx, 1);
  fHistDEDx->GetXaxis()->SetTitle("p_{tot}/sign");
  fHistDEDx->GetYaxis()->SetTitle("TPC signal");
  
  fHistDeDx = new TH2F("fHistDeDx", "dE/dx", 1000, 0.1, 6.0, 1000, 0.0, 1000);
  fHistDeDx->GetYaxis()->SetTitle("TPC dE/dx signal (a.u.)");
  fHistDeDx->GetXaxis()->SetTitle("#frac{p}{z} (GeV/c)");
  BinLogAxis(fHistDeDx);

  fHistDeDxRegion = new TH3F("fHistDeDxRegion", "dE/dx", 400, 0., 6.0, 300, 0., 3., 4, -0.5, 3.5);
  fHistDeDxRegion->GetYaxis()->SetTitle("TPC dE/dx signal (a.u.)");
  fHistDeDxRegion->GetXaxis()->SetTitle("#frac{p}{z} (GeV/c)");

  fHistDeDxSharp = new TH2F("fHistDeDxSharp", "dE/dx", 1000, 0.1, 6.0, 1000, 0.0, 1000);
  fHistDeDxSharp->GetYaxis()->SetTitle("TPC dE/dx signal (a.u.)");
  fHistDeDxSharp->GetXaxis()->SetTitle("#frac{p}{z} (GeV/c)");
  BinLogAxis(fHistDeDxSharp);

  //TOF performance
  fHistTOF3D = new TH3F("fHistTOF3D","mass determination from TOF",400,0.25,4.5,2,-0.5,1.5,2,-1,1);

  fHistTOF2D = new TH2F("fHistTOF2D", "TOF2D", 500, 0.0, 10., 2250, 0.2, 1.1);
  fHistTOF2D->GetYaxis()->SetTitle("#beta");
  fHistTOF2D->GetXaxis()->SetTitle("p (GeV/c)");

  fHistTOFnuclei = new TH2F("fHistTOFnuclei", "TOF", 500, 0.0, 10., 2250, 0.7, 1.);
  fHistTOFnuclei->GetYaxis()->SetTitle("#beta");
  fHistTOFnuclei->GetXaxis()->SetTitle("#frac{p}{z} (GeV/c)");

  //alphas
  fHistAlpha = new TH1F("fHistAlpha", "Anti-Alpha", 22, 1.12, 4.31);
  fHistAlpha->GetYaxis()->SetTitle("Counts");
  fHistAlpha->GetXaxis()->SetTitle("#frac{m^{2}}{z^{2}} (GeV^{2}/c^{4})");

  fHistAlphaSignal  = new TH1F("fHistAlphaSignal", "Anti-Alpha", 22, 1.12, 4.31);
  fHistAlphaSignal->GetYaxis()->SetTitle("Counts");
  fHistAlphaSignal->GetXaxis()->SetTitle("#frac{m^{2}}{z^{2}} (GeV^{2}/c^{4})");

  fGraphAlphaSignal = new TGraph(20);
  fNCounter = 0;
  //
  //  (0.) dca, (1.) sign, (2.) particle Type, (3.) p_tot
  //
  TString axisNameMult[4]={"dca","sign","particleType","ptot"};
  TString axisTitleMult[4]={"dca","sign","particleType","ptot"};
  const Int_t kDcaBins = 76;
  Double_t binsDca[kDcaBins+1] = {-3,-2.85,-2.7,-2.55,-2.4,-2.25,-2.1,-1.95,-1.8,-1.65,-1.5,-1.35,-1.2,-1.05,-0.9,-0.75,-0.6,-0.45,-0.3,-0.285,-0.27,-0.255,-0.24,-0.225,-0.21,-0.195,-0.18,-0.165,-0.15,-0.135,-0.12,-0.105,-0.09,-0.075,-0.06,-0.045,-0.03,-0.015,0,0.015,0.03,0.045,0.06,0.075,0.09,0.105,0.12,0.135,0.15,0.165,0.18,0.195,0.21,0.225,0.24,0.255,0.27,0.285,0.3,0.45,0.6,0.75,0.9,1.05,1.2,1.35,1.5,1.65,1.8,1.95,2.1,2.25,2.4,2.55,2.7,2.85,3};
  //
  //                     (0.) dca, (1.) sign, (2.) particle Type, (3.) p_tot
  Int_t    binsAntiAlpha[4] = {77,         2,                  4,        100};
  Double_t xminAntiAlpha[4] = {-3,        -2,               -0.5,          0};
  Double_t xmaxAntiAlpha[4] = { 3,         2,                3.5,          6};
  //
  fAntiAlpha = new THnF("fAntiAlpha", "AntiAlpha; (0.) dca, (1.) sign, (2.) particle Type, (3.) p_tot", 4, binsAntiAlpha, xminAntiAlpha, xmaxAntiAlpha);
  //
  fAntiAlpha->GetAxis(0)->Set(kDcaBins, binsDca);
  for (Int_t iaxis=0; iaxis<4;iaxis++){
    fAntiAlpha->GetAxis(iaxis)->SetName(axisNameMult[iaxis]);
    fAntiAlpha->GetAxis(iaxis)->SetTitle(axisTitleMult[iaxis]);
  }  
  //
  //  
  //Create histograms for MC
  //generated histogramms
  fHistHelium4PtGen = new TH1F("fHistHelium4PtGen", "Generated  ^{4}He", 200, 0.0, 10.0);
  fHistHelium4PtGen->GetYaxis()->SetTitle("Counts");
  fHistHelium4PtGen->GetXaxis()->SetTitle("p_{T}  (GeV/c)");

  fHistHelium4PtGenPrim = new TH1F("fHistHelium4PtGenPrim", "Primary generated  ^{4}He", 200, 0.0, 10.0);
  fHistHelium4PtGenPrim->GetYaxis()->SetTitle("Counts");
  fHistHelium4PtGenPrim->GetXaxis()->SetTitle("p_{T}  (GeV/c)");

  fHistHelium4PtGenSec = new TH1F("fHistHelium4PtGenSec", "Secondary generated  ^{4}He", 200, 0.0, 10.0);
  fHistHelium4PtGenSec->GetYaxis()->SetTitle("Counts");
  fHistHelium4PtGenSec->GetXaxis()->SetTitle("p_{T}  (GeV/c)");

  fHistHelium4PtGenEta = new TH1F("fHistHelium4PtGenEta", "Generated  ^{4}He in #eta < 0.8", 200, 0.0, 10.0);
  fHistHelium4PtGenEta->GetYaxis()->SetTitle("Counts");
  fHistHelium4PtGenEta->GetXaxis()->SetTitle("p_{T}  (GeV/c)");

  fHistHelium4PtGenPrimEta = new TH1F("fHistHelium4PtGenPrimEta", "Primary generated ^{4}He in #eta < 0.8", 200, 0.0, 10.0);
  fHistHelium4PtGenPrimEta->GetYaxis()->SetTitle("Counts");
  fHistHelium4PtGenPrimEta->GetXaxis()->SetTitle("p_{T}  (GeV/c)");

  fHistAntiHelium4PtGen = new TH1F("fHistAntiHelium4PtGen", "Generated  #bar{^{4}He}", 200, 0.0, 10.0);
  fHistAntiHelium4PtGen->GetYaxis()->SetTitle("Counts");
  fHistAntiHelium4PtGen->GetXaxis()->SetTitle("p_{T}  (GeV/c)");

  fHistAntiHelium4PtGenPrim = new TH1F("fHistAntiHelium4PtGenPrim", "Primary generated  #bar{^{4}He}", 200, 0.0, 10.0);
  fHistAntiHelium4PtGenPrim->GetYaxis()->SetTitle("Counts");
  fHistAntiHelium4PtGenPrim->GetXaxis()->SetTitle("p_{T}  (GeV/c)");

  fHistAntiHelium4PtGenSec = new TH1F("fHistAntiHelium4PtGenSec", "Secondary generated  #bar{^{4}He}", 200, 0.0, 10.0);
  fHistAntiHelium4PtGenSec->GetYaxis()->SetTitle("Counts");
  fHistAntiHelium4PtGenSec->GetXaxis()->SetTitle("p_{T}  (GeV/c)");

  fHistAntiHelium4PtGenEta = new TH1F("fHistAntiHelium4PtGenEta", "Generated  #bar{^{4}He} in #eta < 0.8", 200, 0.0, 10.0);
  fHistAntiHelium4PtGenEta->GetYaxis()->SetTitle("Counts");
  fHistAntiHelium4PtGenEta->GetXaxis()->SetTitle("p_{T}  (GeV/c)");

  //associated histograms
  fHistHelium4PtAso = new TH1F("fHistHelium4PtAso", "Associated  ^{4}He", 200, 0.0, 10.0);
  fHistHelium4PtAso->GetYaxis()->SetTitle("Counts");
  fHistHelium4PtAso->GetXaxis()->SetTitle("p_{T}  (GeV/c)");

  fHistHelium4PtAsoPrim = new TH1F("fHistHelium4PtAsoPrim", "Associated prim ^{4}He", 200, 0.0, 10.0);
  fHistHelium4PtAsoPrim->GetYaxis()->SetTitle("Counts");
  fHistHelium4PtAsoPrim->GetXaxis()->SetTitle("p_{T}  (GeV/c)");

  fHistHelium4PtAsoSec = new TH1F("fHistHelium4PtAsoSec", "Associated sec  ^{4}He", 200, 0.0, 10.0);
  fHistHelium4PtAsoSec->GetYaxis()->SetTitle("Counts");
  fHistHelium4PtAsoSec->GetXaxis()->SetTitle("p_{T}  (GeV/c)");

  fHistAntiHelium4PtAso = new TH1F("fHistAntiHelium4PtAso", "Associated  #bar{^{4}He}", 200, 0.0, 10.0);
  fHistAntiHelium4PtAso->GetYaxis()->SetTitle("Counts");
  fHistAntiHelium4PtAso->GetXaxis()->SetTitle("p_{T}  (GeV/c)");

  //
  //
  fOutputContainer = new TObjArray(1);
  fOutputContainer->SetOwner(kTRUE);
  fOutputContainer->SetName(GetName());
  fOutputContainer->Add(fHistCentralityClass10);
  fOutputContainer->Add(fHistCentralityPercentile);
  fOutputContainer->Add(fHistTriggerStat);
  fOutputContainer->Add(fHistTriggerStatAfterEventSelection);
  fOutputContainer->Add(fHistDEDx);
  fOutputContainer->Add(fHistDeDx);
  fOutputContainer->Add(fHistDeDxRegion);
  fOutputContainer->Add(fHistDeDxSharp);
  fOutputContainer->Add(fAntiAlpha);  
  fOutputContainer->Add(fHistAlpha);
  fOutputContainer->Add(fHistAlphaSignal);
  fOutputContainer->Add(fGraphAlphaSignal);
  fOutputContainer->Add(fHistTOF3D);
  fOutputContainer->Add(fHistTOF2D);
  fOutputContainer->Add(fHistTOFnuclei);
  fOutputContainer->Add(fHistHelium4PtGen);
  fOutputContainer->Add(fHistHelium4PtGenPrim);
  fOutputContainer->Add(fHistHelium4PtGenSec);
  fOutputContainer->Add(fHistHelium4PtGenEta);
  fOutputContainer->Add(fHistHelium4PtGenPrimEta);
  fOutputContainer->Add(fHistAntiHelium4PtGen);
  fOutputContainer->Add(fHistAntiHelium4PtGenPrim);
  fOutputContainer->Add(fHistAntiHelium4PtGenSec);
  fOutputContainer->Add(fHistAntiHelium4PtGenEta);
  fOutputContainer->Add(fHistHelium4PtAso);
  fOutputContainer->Add(fHistHelium4PtAsoPrim);
  fOutputContainer->Add(fHistHelium4PtAsoSec);
  fOutputContainer->Add(fHistAntiHelium4PtAso);


  //------------ Tree and branch definitions ----------------//  
  fTree = new TTree("tree", "alpha tree");     
  //------------ Event variables ------------//
  fTree->Branch("fName",fName,"fName/C");
  fTree->Branch("fEvnt",&fEvnt, "fEvnt/I"); 
  fTree->Branch("fFileName",fFileName,"fFileName/C");
  fTree->Branch("fEventNumber",fEventNumber,"fEventNumber/I");
  fTree->Branch("fItrk",&fItrk, "fItrk/I"); 
  //-------------------------------------------//  
  //----------- Track variables --------------//  
  fTree->Branch("fEta",fEta,"fEta[fItrk]/D");
  fTree->Branch("fKinkIndex",fKinkIndex,"fKinkIndex[fItrk]/I");
  fTree->Branch("fCentrality",fCentrality,"fCentrality[fItrk]/F");
  //
  fTree->Branch("fTPCnCluster",fTPCnCluster,"fTPCnCluster[fItrk]/s");
  fTree->Branch("fTPCNsignal",fTPCNsignal,"fTPCNsignal[fItrk]/s");
  fTree->Branch("fChi2PerClusterTPC",fChi2PerClusterTPC,"fChi2PerClusterTPC[fItrk]/D");
  fTree->Branch("fTPCRefit",fTPCRefit,"fTPCRefit[fItrk]/O");
  fTree->Branch("fTPCsignal0",fTPCsignal0,"fTPCsignal0[fItrk]/D");
  fTree->Branch("fTPCsignal1",fTPCsignal1,"fTPCsignal1[fItrk]/D");
  fTree->Branch("fTPCsignal2",fTPCsignal2,"fTPCsignal2[fItrk]/D");
  fTree->Branch("fTPCsignal3",fTPCsignal3,"fTPCsignal3[fItrk]/D");
  fTree->Branch("fTPCSharedClusters",fTPCSharedClusters,"fTPCSharedClusters[fItrk]/I");
  fTree->Branch("fTPCNclsIter1",fTPCNclsIter1,"fTPCNclsIter1[fItrk]/s");
  //
  fTree->Branch("fITSsignal",fITSsignal,"fITSsignal[fItrk]/D");
  fTree->Branch("fITSnCluster",fITSnCluster,"fITSnCluster[fItrk]/I");
  fTree->Branch("fChi2PerClusterITS",fChi2PerClusterITS,"fChi2PerClusterITS[fItrk]/D");
  fTree->Branch("fITSRefit",fITSRefit,"fITSRefit[fItrk]/O");
  //
  fTree->Branch("fTOFtime",fTOFtime,"fTOFtime[fItrk]/O");
  fTree->Branch("fTOFRefit",fTOFRefit,"fTOFRefit[fItrk]/O");
  fTree->Branch("fTOFout",fTOFout,"fTOFout[fItrk]/O");
  fTree->Branch("fTOFsignalDz",fTOFsignalDz,"fTOFsignalDz[fItrk]/D");
  fTree->Branch("fTOFsignalDx",fTOFsignalDx,"fTOFsignalDx[fItrk]/D");
  //
  fTree->Branch("fTRDin",fTRDin,"fTRDin[fItrk]/O");
  //
  fTree->Branch("fDCAXY",fDCAXY,"fDCAXY[fItrk]/F");
  fTree->Branch("fDCAZ",fDCAZ,"fDCAZ[fItrk]/F");
  //
  fTree->Branch("fTrkPtot",fTrkPtot,"fTrkPtot[fItrk]/D"); 
  fTree->Branch("fTPCPtot",fTPCPtot,"fTPCPtot[fItrk]/D"); 
  fTree->Branch("fTrackPt",fTrackPt,"fTrackPt[fItrk]/D");
  fTree->Branch("fDeDx",fDeDx,"fDeDx[fItrk]/D");  
  fTree->Branch("fSign",fSign,"fSign[fItrk]/D");  
  fTree->Branch("fMass",fMass,"fMass[fItrk]/F");
  fTree->Branch("fTime",fTime,"fTime[fItrk]/F");
  fTree->Branch("fLength",fLength,"fLength[fItrk]/F");
  fTree->Branch("fSigmaQP",fSigmaQP,"fSigmaQP[fItrk]/D");
  //
  fTree->Branch("fAssociated",fAssociated,"fAssociated[fItrk]/O");

}

//________________________________________________________________________
void AliAnalysisTaskAntiHe4::UserExec(Option_t *)
{
  // Main loop
  // Called for each event

  //get Event-Handler for the trigger information
  fEventHandler= dynamic_cast<AliInputEventHandler*> (AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
  if (!fEventHandler) {
    AliError("Could not get InputHandler");
    //return -1;
    return;
  }

  // Monte Carlo
  AliMCEventHandler* eventHandler = dynamic_cast<AliMCEventHandler*> (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler());
  if (!eventHandler){ 
    //Printf("ERROR: Could not retrieve MC event handler");
    fMCtrue = kFALSE;
  }

  AliMCEvent* mcEvent = 0x0;
  AliStack* stack = 0x0;
  if (eventHandler) mcEvent = eventHandler->MCEvent();
  if (!mcEvent){
    //Printf("ERROR: Could not retrieve MC event");
    if (fMCtrue) return;
  }
  
  if (fMCtrue){
    stack = mcEvent->Stack();
    if (!stack) return;
  }
  
  //look for the generated particles
  MCGenerated(stack);
  
  fESD=dynamic_cast<AliESDEvent*>(InputEvent());
  if (!fESD) {
    //Printf("ERROR: fESD not available");
    return;
  }

  if (SetupEvent() < 0) {
    PostData(1, fOutputContainer);
    return;
  }

  const AliESDVertex *vertex = fESD->GetPrimaryVertexTracks();
  if(vertex->GetNContributors()<1) {
    // SPD vertex
    vertex = fESD->GetPrimaryVertexSPD();
    if(vertex->GetNContributors()<1) {
      PostData(1, fOutputContainer);
      return;
    }
    
  }  
  // check if event is selected by physics selection class
  //
  Bool_t isSelected = kFALSE;
  isSelected = ((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected();
  if (!isSelected || TMath::Abs(vertex->GetZ()) > 10) {
    PostData(1, fOutputContainer);
    return;
  }

  //
  //centrality selection in PbPb 
  //
  Int_t centralityClass10 = -1;
  Float_t centralityPercentile = -1;
  //
  if (fESD->GetEventSpecie() == 4) { //species == 4 == PbPb
    // PbPb
    AliCentrality *esdCentrality = fESD->GetCentrality();
    centralityClass10 = esdCentrality->GetCentralityClass10("V0M"); // centrality percentile determined with V0                           
    centralityPercentile = esdCentrality->GetCentralityPercentile("V0M"); // centrality percentile determined with V0                        
    if(!fMCtrue){
      if (centralityPercentile < 0. || centralityPercentile > 80. ) return; //select only events with centralities between 0 and 80 %  
    }
  }
  //
//  if (!fTriggerFired[0] && !fTriggerFired[1] && !fTriggerFired[2]) return; // select only events which pass kMB, kCentral, kSemiCentral
  //
  fHistCentralityClass10->Fill(centralityClass10);
  fHistCentralityPercentile->Fill(centralityPercentile);
  //
  if(fTriggerFired[0] == kTRUE) fHistTriggerStatAfterEventSelection->Fill(0);
  if(fTriggerFired[1] == kTRUE) fHistTriggerStatAfterEventSelection->Fill(1);
  if(fTriggerFired[2] == kTRUE) fHistTriggerStatAfterEventSelection->Fill(2);
  if(fTriggerFired[3] == kTRUE) fHistTriggerStatAfterEventSelection->Fill(3);
  if(fTriggerFired[4] == kTRUE) fHistTriggerStatAfterEventSelection->Fill(4);
  //
  //
  if (!fESDpid) fESDpid = ((AliESDInputHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->GetESDpid();
  if (!fESDpid) {
    fESDpid = new AliESDpid(); // HACK FOR MC PBPB --> PLEASE REMOVE AS SOON AS POSSIBLE
    fESDpid->GetTPCResponse().SetBetheBlochParameters(1.28778e+00/50., 3.13539e+01, TMath::Exp(-3.16327e+01), 1.87901e+00, 6.41583e+00);
  }
  //
  Float_t dca[2], cov[3]; // dca_xy, dca_z, sigma_xy, sigma_xy_z, sigma_z for th // for Anti-Alpha
  fEvnt =  fESD->GetEventNumberInFile();
  sscanf(fInputHandler->GetTree()->GetCurrentFile()->GetName(),"%s", fName);
  fItrk = 0;
  //
  Int_t runNumber = 0;
  runNumber = fESD->GetRunNumber();
  //
  Bool_t fillTree = kFALSE;
  // Track loop to fill the spectram
  for (Int_t iTracks = 0; iTracks < fESD->GetNumberOfTracks(); iTracks++) {

    fEventNumber[fItrk] = -1;

    fEta[fItrk] = -2;
    fCentrality[fItrk] = -1;
    fTPCNsignal[fItrk] = -1;
    fTPCnCluster[fItrk] = -1;
    fChi2PerClusterTPC[fItrk] = -1;
    fTPCRefit[fItrk] = kFALSE;
    fTPCsignal0[fItrk] = -1;
    fTPCsignal1[fItrk] = -1;
    fTPCsignal2[fItrk] = -1;
    fTPCsignal3[fItrk] = -1;
    fTPCSharedClusters[fItrk] = -1;
    fTPCNclsIter1[fItrk] = -1;

    fITSsignal[fItrk] = -1;
    fITSnCluster[fItrk] = -1;
    fChi2PerClusterITS[fItrk] = -1;
    fITSRefit[fItrk] = kFALSE;

    fTOFRefit[fItrk] = kFALSE;
    fTOFtime[fItrk] = kFALSE;
    fTOFout[fItrk] = kFALSE;
    fTOFsignalDz[fItrk] = -1;
    fTOFsignalDx[fItrk] = -1;

    fTRDin[fItrk] = kFALSE;

    fDCAZ[fItrk] = -1;
    fDCAXY[fItrk] = -1;

    fTrkPtot[fItrk] = -1;
    fTPCPtot[fItrk] = -1;
    fTrackPt[fItrk] = -1;
    fDeDx[fItrk] = -1;
    fSign[fItrk] = -2;
    fMass[fItrk] = -1;
    fTime[fItrk] = -1;
    fLength[fItrk] = -1;
    fSigmaQP[fItrk] = -1;    

    fAssociated[fItrk] = kFALSE;

    AliESDtrack* track = dynamic_cast<AliESDtrack*>(fESD->GetTrack(iTracks));
    if (!fESDtrackCuts->AcceptTrack(track)) continue;
    //
    Double_t nClustersTPCPID=track->GetTPCsignalN();
    if(nClustersTPCPID < 60) continue;
    //
    if (!track->GetInnerParam()) continue;
    //
    Double32_t signal[4] = {0,0,0,0};
    Char_t ncl[3];
    Char_t nrows[3];
    if (track->GetTPCdEdxInfo()) {
      track->GetTPCdEdxInfo()->GetTPCSignalRegionInfo(signal,ncl,nrows);
    }
    //
    Double_t ptot = track->GetInnerParam()->GetP();
    Double_t ptotInc = track->GetP(); // total momentum of the incoming particle
    Double_t sign = track->GetSign();
    //
    Double_t eta = TMath::Abs(track->Eta());
    TBits shared = track->GetTPCSharedMap();

    track->GetImpactParameters(dca, cov);
    Float_t dcaXY = TMath::Sqrt(dca[0]*dca[0]);
    Float_t dcaXYsign = dca[0];
    Float_t dcaZ  = TMath::Sqrt(dca[1]*dca[1]);
    //
    Double_t cov1[15];
    track->GetExternalCovariance(cov1);//->GetSigmaQoverP();
    //
    Double_t tpcSignal = track->GetTPCsignal();
    //
    //PID via specific energy loss in the TPC
    //define the arrays for the Bethe-Bloch-Parameters
    //
    
    SetBBParameters(runNumber);
    
    //define expected signals for the various species
    Double_t expSignalProton = AliExternalTrackParam::BetheBlochAleph(ptot/0.93827,fBBParametersLightParticles[0],fBBParametersLightParticles[1],fBBParametersLightParticles[2],fBBParametersLightParticles[3],fBBParametersLightParticles[4]);
    Double_t expSignalDeuteron = AliExternalTrackParam::BetheBlochAleph(ptot/1.875612,fBBParametersLightParticles[0],fBBParametersLightParticles[1],fBBParametersLightParticles[2],fBBParametersLightParticles[3],fBBParametersLightParticles[4]);
    Double_t expSignalTriton = AliExternalTrackParam::BetheBlochAleph(ptot/2.808921,fBBParametersLightParticles[0],fBBParametersLightParticles[1],fBBParametersLightParticles[2],fBBParametersLightParticles[3],fBBParametersLightParticles[4]);

    Double_t expSignalHelium3 = 4*AliExternalTrackParam::BetheBlochAleph(2*ptot/2.80941,fBBParametersNuclei[0],fBBParametersNuclei[1],fBBParametersNuclei[2],fBBParametersNuclei[3],fBBParametersNuclei[4]);
    Double_t expSignalHelium4 = 4*AliExternalTrackParam::BetheBlochAleph(2*ptot/3.728401,fBBParametersNuclei[0],fBBParametersNuclei[1],fBBParametersNuclei[2],fBBParametersNuclei[3],fBBParametersNuclei[4]);

    //
    Float_t mass = 0;
    Float_t time = -1; 
    Float_t beta = 0;
    Float_t gamma = -1;
    Bool_t  hasTOFout  = kFALSE;
    Bool_t  hasTOFtime = kFALSE;
    Float_t length = track->GetIntegratedLength();
    UInt_t  status = track->GetStatus();
    Bool_t  isAssociated = kFALSE;

    if (length > 350) {
      time = track->GetTOFsignal() - fESDpid->GetTOFResponse().GetStartTime(track->P());
      if (time > 0) {
	beta = length / (2.99792457999999984e-02 * time);
	gamma = 1/TMath::Sqrt(1 - beta*beta);
	mass = ptot/TMath::Sqrt(gamma*gamma - 1); // using inner TPC mom. as approx.
      }
    }
    //
    // fill tree and print candidates (for short list)
    //
    Float_t cut = 4*AliExternalTrackParam::BetheBlochAleph(2*ptot/(0.938*3),1.1931,
							   31.9806,
							   5.04114e-11,
							   2.13096,
							   2.38541);
    if (fMCtrue) cut = 4*AliExternalTrackParam::BetheBlochAleph(2*ptot/(0.938*3),0.9931,
								31.9806,
								5.04114e-11,
								2.13096,
								2.38541);
    Bool_t IsDeuteron = kFALSE;
    Bool_t IsTriton = kFALSE;

    Double_t DeuteronSigma = TMath::Abs(tpcSignal - expSignalDeuteron)/expSignalDeuteron;
    Double_t TritonSigma = TMath::Abs(tpcSignal - expSignalTriton)/expSignalTriton;
    

    if(DeuteronSigma < 0.3 && runNumber < 166500) IsDeuteron = kTRUE;
    if(TritonSigma < 0.3 && runNumber < 166500) IsTriton = kTRUE;

    //cout << "isDeuteron: " << IsDeuteron << endl;
    
    //if(tpcSignal > cut || IsDeuteron == kTRUE || IsTriton == kTRUE){
    if(tpcSignal > cut ){
 
      //cout << "isDeuteron: " << IsDeuteron << endl;
      if (eta < 1.0 && tpcSignal < 1000 && dcaZ < 15 && dcaXY < 15 && ptot < 20 && ptot > 1.0 && tpcSignal > 120){// && track->GetTPCsignalN() > 60) { // && ptot > 1.0 &6 tpcSignal > 120
	//
	cout << "AntiAlphaEvent" << " " 
	     << AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()->GetTree()->GetCurrentFile()->GetName() << " " 
	     << "event number in file: " << fESD->GetEventNumberInFile() 
	     << " Index " << iTracks 
	     << " ptot: " << ptot 
	     << " sig: " << tpcSignal <<endl;
	//
	fillTree = kTRUE;
	//
	
	sscanf(fInputHandler->GetTree()->GetCurrentFile()->GetName(),"%s", fFileName);
	fEventNumber[fItrk] = fESD->GetEventNumberInFile();
	
	fEta[fItrk] = eta;
	fKinkIndex[fItrk] = track->GetKinkIndex(0);
	fCentrality[fItrk] = centralityPercentile;
	
	fTPCNsignal[fItrk] = track->GetTPCsignalN();
	fTPCnCluster[fItrk] = track->GetTPCNcls();
	fChi2PerClusterTPC[fItrk] = track->GetTPCchi2()/fTPCnCluster[fItrk];
	if(status&AliESDtrack::kTPCrefit)
	  fTPCRefit[fItrk] = kTRUE;
	else fTPCRefit[fItrk] = kFALSE;
	fTPCsignal0[fItrk] = signal[0];
	fTPCsignal1[fItrk] = signal[1];
	fTPCsignal2[fItrk] = signal[2];
	fTPCsignal3[fItrk] = signal[3];
	fTPCSharedClusters[fItrk] = shared.CountBits();
	fTPCNclsIter1[fItrk] = track->GetTPCNclsIter1();
	
	fITSsignal[fItrk] = track->GetITSsignal();
	fITSnCluster[fItrk] = track->GetNcls(0);
	fChi2PerClusterITS[fItrk] = track->GetITSchi2()/fITSnCluster[fItrk];
	if(status&AliESDtrack::kITSrefit)
	  fITSRefit[fItrk] = kTRUE;
	else fITSRefit[fItrk] = kFALSE;
	
	
	if(status&AliESDtrack::kITSrefit)
	  fITSRefit[fItrk] = kTRUE;
	else fITSRefit[fItrk] = kFALSE;
	hasTOFout = status&AliESDtrack::kTOFout;
	hasTOFtime  = status&AliESDtrack::kTIME;
	fTOFtime[fItrk] = hasTOFtime;
	fTOFout[fItrk]  = hasTOFout;
	fTOFsignalDz[fItrk] = track->GetTOFsignalDz();
	fTOFsignalDx[fItrk] = track->GetTOFsignalDx();
	
	fTRDin[fItrk] = status&AliESDtrack::kTRDin;
	
	fDCAZ[fItrk] = dcaXY;
	fDCAXY[fItrk] = dcaZ;
	
	fTrkPtot[fItrk] = track->P();
	fTPCPtot[fItrk] = ptot;
	fTrackPt[fItrk] = track->Pt();
	fDeDx[fItrk] = tpcSignal;
	fSign[fItrk] = sign;
	fMass[fItrk] = mass;
	fTime[fItrk] = time;
	fLength[fItrk] = length;
	fSigmaQP[fItrk] = cov1[14];

	if (fMCtrue){ //associated
	  
	  Int_t label  = track->GetLabel();
	  TParticle *tparticle = stack->Particle(TMath::Abs(label));
	  
	  Bool_t isPrimary = stack->IsPhysicalPrimary(TMath::Abs(label));
	  Bool_t isSecondary = stack->IsSecondaryFromMaterial(TMath::Abs(label));
	  
	  Long_t pdgCode = tparticle->GetPdgCode();
	  Double_t pT =(track->Pt())*2;
	  
	  if(pdgCode == 1000020040)
	    {
	      fHistHelium4PtAso->Fill(pT);
	      if(isPrimary) fHistHelium4PtAsoPrim->Fill(pT);
	      if(isSecondary)  fHistHelium4PtAsoSec->Fill(pT);
	      isAssociated = kTRUE;
	    }
	  
	  if(pdgCode == -1000020040)
	    {
	      fHistAntiHelium4PtAso->Fill(pT);
	      isAssociated = kTRUE;
	    }
	  
	}
	
	fAssociated[fItrk] = isAssociated;
	
	fItrk++;
      }
    }
    //
    // do pid fill histogram for raw ratios
    //
    //                       (0.) dca, (1.) sign, (2.) particle Type, (3.) p_tot
    Int_t id = -1;
    if (ptot > 0.2 && TMath::Abs(tpcSignal - expSignalProton)/expSignalProton < 0.2) id = 0;
    if (ptot > 0.3 && TMath::Abs(tpcSignal - expSignalDeuteron)/expSignalDeuteron < 0.2) id = 1;
    if (ptot > 0.7 && TMath::Abs(tpcSignal - expSignalTriton)/expSignalTriton < 0.2) id = 2;
    if (ptot > 0.5 && (tpcSignal - expSignalHelium3)/expSignalHelium3 > -0.1 &&  (tpcSignal - expSignalHelium3)/expSignalHelium3 < 0.2) id = 3;
    //
    Double_t vecAntiAlpha[4] = {dcaXYsign, sign, static_cast<Double_t>(id), ptotInc};
    if (id != -1 && tpcSignal > 120) fAntiAlpha->Fill(vecAntiAlpha);
    //
    // fill final histograms
    //
    if (!fESDtrackCutsSharp->AcceptTrack(track) || shared.CountBits() > 1 || track->GetTPCsignalN() < 80 || track->GetKinkIndex(0) != 0 || track->GetTPCNclsIter1() < 80) continue;
    //
    fHistDEDx->Fill(ptot,track->GetTPCsignal(), sign);
    if (TMath::Abs(tpcSignal - expSignalHelium4)/expSignalHelium4 < 0.2 && time < 99998) fHistTOF3D->Fill(mass,1,sign);  
    //
    // remove overlapping tracks with special dE/dx cut
    //
    //if (signal[1] < 6) continue;
    //if (signal[0]/signal[1] > 1.6 || signal[2]/signal[1] > 1.6 || signal[3]/signal[1] > 1.3) continue;
    //
    if(sign<0) {
      fHistDeDx->Fill(ptot, track->GetTPCsignal());
      if (track->GetTPCsignalN() > 100 &&
          TMath::Abs(track->Eta()) < 1.0 &&
          signal[3]/signal[1] > 0.6 &&
          signal[0]/signal[1] > 0.5 &&
          signal[3]/signal[1] < 1.2 &&
          track->GetTPCNclsIter1() > 70 &&
          track->GetTPCnclsS() < 10) {
        fHistDeDxSharp->Fill(ptot, track->GetTPCsignal());
      }
    }
    //
    // dE/dx for specific regions
    //
    for(Int_t iSig = 0; iSig < 4; iSig++) {
      if (signal[1] > 6) fHistDeDxRegion->Fill(ptot,signal[iSig]/signal[1],iSig);
    }
    //
    // alpha TOF plot
    //
    if((track->GetTPCsignal()-expSignalHelium4)/expSignalHelium4 > -0.15 && (track->GetTPCsignal()-expSignalHelium4)/expSignalHelium4 < 0.15) {
      //TOF
      hasTOFout  = status&AliESDtrack::kTOFout;
      hasTOFtime = status&AliESDtrack::kTIME;
      Bool_t hasTOF     = kFALSE;
      
      if (hasTOFout && hasTOFtime) hasTOF = kTRUE;
      //
      if (length < 350.) hasTOF = kFALSE;
      
      Float_t time0 = fESDpid->GetTOFResponse().GetStartTime(track->P());//fESDpid->GetTOFResponse().GetTimeZero();
      //                              
      if (length > 350. &&  hasTOF == kTRUE && ptot < 4) {
	time = track->GetTOFsignal() - time0;
	if (time > 0) {
	  beta = length / (2.99792457999999984e-02 * time);
	  if (beta < 0.975) {
	    gamma = 1/TMath::Sqrt(1 - beta*beta);
	    mass = ptot/TMath::Sqrt(gamma*gamma - 1); // using inner TPC mom. as approx.
	    if (TMath::Sqrt(track->GetTOFsignalDz()*track->GetTOFsignalDz() + track->GetTOFsignalDx()*track->GetTOFsignalDx()) < 5.) {
	      fHistAlpha->Fill(mass*mass);
	      if (mass*mass > 3. && mass*mass < 4.) {
		fHistAlphaSignal->Fill(mass*mass);
		fGraphAlphaSignal->SetPoint(fNCounter, ptot, track->GetTPCsignal());
		fNCounter++;
	      }
	      fHistTOFnuclei->Fill(ptot,beta);
	    }
	  }
	}
      }
    }

  }//end loop over tracks

  if (fillTree) fTree->Fill();

  // Post output data.
  PostData(1, fOutputContainer);
  PostData(2, fTree);
}      

//________________________________________________________________________
void AliAnalysisTaskAntiHe4::Terminate(const Option_t *)
{
  // Draw result to the screen
  // Called once at the end of the query

  //get output data and darw 'fHistPt'
  if (!GetOutputData(0)) return;

}


//________________________________________________________________________
void AliAnalysisTaskAntiHe4::BinLogAxis(const THn *h, Int_t axisNumber) {
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
void AliAnalysisTaskAntiHe4::BinLogAxis(const TH3 *h, Int_t axisNumber) {
  //
  // Method for the correct logarithmic binning of histograms
  //
  TAxis *axis = const_cast<TAxis*>(h->GetXaxis());
  if (axisNumber == 1) axis = const_cast<TAxis*>(h->GetYaxis());
  if (axisNumber == 2) axis = const_cast<TAxis*>(h->GetZaxis());
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
void AliAnalysisTaskAntiHe4::BinLogAxis(const TH1 *h) {
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
//_____________________________________________________________________________
Int_t AliAnalysisTaskAntiHe4::Initialize() {


  // -- Reset Event
  // ----------------
  ResetEvent();

  return 0;
}
//________________________________________________________________________
Int_t AliAnalysisTaskAntiHe4::SetupEvent() {
  // Setup Reading of event

  ResetEvent();
  // -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
  // -- Get Trigger
  // -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --

  //Bool_t isTriggered = IsTriggered();
  IsTriggered();
  return 0;
}
//_____________________________________________________________________________
void AliAnalysisTaskAntiHe4::ResetEvent() {
  // Reset event
  // -- Reset QA
  for (Int_t ii = 0; ii < fNTriggers; ++ii)
    fTriggerFired[ii] = kFALSE;

  return;
}
//________________________________________________________________________
Bool_t AliAnalysisTaskAntiHe4::IsTriggered() {
  // Check if Event is triggered and fill Trigger Histogram

  if ((fEventHandler->IsEventSelected() & AliVEvent::kMB))          fTriggerFired[0] = kTRUE;
  if ((fEventHandler->IsEventSelected() & AliVEvent::kCentral))     fTriggerFired[1] = kTRUE;
  if ((fEventHandler->IsEventSelected() & AliVEvent::kSemiCentral)) fTriggerFired[2] = kTRUE;
  if ((fEventHandler->IsEventSelected() & AliVEvent::kEMCEJE))      fTriggerFired[3] = kTRUE;
  if ((fEventHandler->IsEventSelected() & AliVEvent::kEMCEGA))      fTriggerFired[4] = kTRUE;

  Bool_t isTriggered = kFALSE;

  for (Int_t ii=0; ii<fNTriggers; ++ii) {
    if(fTriggerFired[ii]) {
      isTriggered = kTRUE;
      fHistTriggerStat->Fill(ii);
    }
  }

  return isTriggered;
}
//________________________________________________________________________
void AliAnalysisTaskAntiHe4::SetBBParameters(Int_t runNumber){

  if(runNumber < 166500){ //LHC10h
    
    fBBParametersLightParticles[0]   = 1.45802; 
    fBBParametersLightParticles[1]   = 27.4992;
    fBBParametersLightParticles[2]   = 4.00313e-15;
    fBBParametersLightParticles[3]   = 2.48485;
    fBBParametersLightParticles[4]   = 8.31768;

    fBBParametersNuclei[0]  = 1.69155;
    fBBParametersNuclei[1]  = 27.4992;
    fBBParametersNuclei[2]  = 4.00313e-15;
    fBBParametersNuclei[3]  = 2.48485;
    fBBParametersNuclei[4]  = 8.31768;

  }

  if(runNumber > 166500 && runNumber <= 215150){ //LHC11h

    fBBParametersLightParticles[0]   = 4.69637;//1.11243;
    fBBParametersLightParticles[1]   = 7.51827;//26.1144;
    fBBParametersLightParticles[2]   = 0.0183746;//4.00313e-15;
    fBBParametersLightParticles[3]   = 2.60;//2.72969;
    fBBParametersLightParticles[4]   = 2.7;//9.15038;

    fBBParametersNuclei[0]  = 1.4906;
    fBBParametersNuclei[1]  = 27.9758;
    fBBParametersNuclei[2]  = 4.00313e-15;
    fBBParametersNuclei[3]  = 2.50804;
    fBBParametersNuclei[4]  = 8.31768;

  }
    if(runNumber > 215150){ //LHC RUN 2
        
        fBBParametersLightParticles[0]   = 1.45802;
        fBBParametersLightParticles[1]   = 27.4992;
        fBBParametersLightParticles[2]   = 4.00313e-15;
        fBBParametersLightParticles[3]   = 2.48485;
        fBBParametersLightParticles[4]   = 8.31768;
        
        fBBParametersNuclei[0]  = 1.69155;
        fBBParametersNuclei[1]  = 27.4992;
        fBBParametersNuclei[2]  = 4.00313e-15;
        fBBParametersNuclei[3]  = 2.48485;
        fBBParametersNuclei[4]  = 8.31768;
        
    }

    
}
//______________________________________________________________________________
void AliAnalysisTaskAntiHe4::MCGenerated(AliStack* stack) 
{ 

  // Monte Carlo for genenerated particles
  if (fMCtrue) //MC loop  
    {
 
      Int_t stackN = 0;

      for(stackN = 0; stackN < stack->GetNtrack(); stackN++) //loop over stack
	{

	  Bool_t isPrimary = stack->IsPhysicalPrimary(stackN);
	  Bool_t isSecondary = stack->IsSecondaryFromMaterial(stackN);

	  const TParticle *tparticle = stack->Particle(stackN);
	  Long_t pdgCode = tparticle->GetPdgCode();

	  Double_t pTGen = tparticle->Pt();
	  Double_t eta = tparticle->Eta();
	  //check which particle it is 

	  //Alpha
	  if(pdgCode == 1000020040)
	    {
	      fHistHelium4PtGen->Fill(pTGen);
	      if(isPrimary) fHistHelium4PtGenPrim->Fill(pTGen);
	      if(isSecondary) fHistHelium4PtGenSec->Fill(pTGen);
	      if(TMath::Abs(eta) < 0.8)fHistHelium4PtGenEta->Fill(pTGen);
	      if(isPrimary && TMath::Abs(eta) < 0.8)fHistHelium4PtGenPrimEta->Fill(pTGen);
	    }

	  //Anti-Alpha
	  if(pdgCode == -1000020040) 
	    {
	      fHistAntiHelium4PtGen->Fill(pTGen);
	      if(isPrimary) fHistAntiHelium4PtGenPrim->Fill(pTGen);
	      if(isSecondary) fHistAntiHelium4PtGenSec->Fill(pTGen);
	      if(TMath::Abs(eta) < 0.8)fHistAntiHelium4PtGenEta->Fill(pTGen);
	    }

  	      
	}//end loop over stack
      

    }//end MC
}
