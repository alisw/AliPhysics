#include "TChain.h"
#include "TList.h"
#include "TCanvas.h"
#include "TParticle.h"
#include "TLorentzVector.h"
#include "TGraphErrors.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TArrayF.h"
#include "TF1.h"
#include "TRandom.h"
#include "TFile.h"

#include "AliAnalysisManager.h"
#include "AliLog.h"

#include "AliEventPoolManager.h"           

#include "AliAnalysisTaskToyModel.h"
#include "AliBalance.h"
#include "AliBalancePsi.h"
#include "AliAnalysisTaskTriggeredBF.h"


// Analysis task for the toy model analysis
// Authors: Panos.Christakoglou@nikhef.nl

using std::cout;
using std::endl;

ClassImp(AliAnalysisTaskToyModel)

//________________________________________________________________________
AliAnalysisTaskToyModel::AliAnalysisTaskToyModel() 
: TObject(),
  fUseDebug(kFALSE),
  fBalance(0),
  fRunShuffling(kFALSE), fShuffledBalance(0),
  fRunMixing(kFALSE), fMixedBalance(0), fPoolMgr(0),
  fList(0), fListBF(0), fListBFS(0), fListBFM(0),
  fHistEventStats(0),
  fHistNumberOfAcceptedParticles(0),
  fHistReactionPlane(0),
  fHistEtaTotal(0), fHistEta(0),
  fHistEtaPhiPos(0), fHistEtaPhiNeg(0),
  fHistRapidity(0),
  fHistRapidityPions(0), fHistRapidityKaons(0), fHistRapidityProtons(0),
  fHistPhi(0),
  fHistPhiPions(0), fHistPhiKaons(0), fHistPhiProtons(0),
  fHistPt(0),
  fHistPtPions(0), fHistPtKaons(0), fHistPtProtons(0),
  fTotalMultiplicityMean(100.), fTotalMultiplicitySigma(0.0),
  fNetChargeMean(0.0), fNetChargeSigma(0.0),
  fPtMin(0.0), fPtMax(100.0),
  fEtaMin(-1.0), fEtaMax(1.0),
  fUseAcceptanceParameterization(kFALSE), fAcceptanceParameterization(0),
  fSimulateDetectorEffects(kFALSE),
  fNumberOfInefficientSectors(0),
  fInefficiencyFactorInPhi(1.0),
  fNumberOfDeadSectors(0),
  fEfficiencyDropNearEtaEdges(kFALSE),
  fEfficiencyMatrix(0),
  fUseAllCharges(kFALSE), fParticleMass(0.0),
  fPtSpectraAllCharges(0), fTemperatureAllCharges(100.),
  fReactionPlane(0.0),
  fAzimuthalAngleAllCharges(0), fDirectedFlowAllCharges(0.0), 
  fEllipticFlowAllCharges(0.0), fTriangularFlowAllCharges(0.0),
  fQuandrangularFlowAllCharges(0.0), fPentangularFlowAllCharges(0.0),
  fPionPercentage(0.8), fPionMass(0.0),
  fPtSpectraPions(0), fTemperaturePions(100.),
  fAzimuthalAnglePions(0), fDirectedFlowPions(0.0), 
  fEllipticFlowPions(0.0), fTriangularFlowPions(0.0), 
  fQuandrangularFlowPions(0.0), fPentangularFlowPions(0.0),
  fKaonPercentage(0.8), fKaonMass(0.0),
  fPtSpectraKaons(0), fTemperatureKaons(100.),
  fAzimuthalAngleKaons(0), fDirectedFlowKaons(0.0), 
  fEllipticFlowKaons(0.0), fTriangularFlowKaons(0.0),
  fQuandrangularFlowKaons(0.0), fPentangularFlowKaons(0.0),
  fProtonPercentage(0.8), fProtonMass(0.0),
  fPtSpectraProtons(0), fTemperatureProtons(100.),
  fAzimuthalAngleProtons(0), fDirectedFlowProtons(0.0), 
  fEllipticFlowProtons(0.0), fTriangularFlowProtons(0.0),
  fQuandrangularFlowProtons(0.0), fPentangularFlowProtons(0.0),
  fUseDynamicalCorrelations(kFALSE), fDynamicalCorrelationsPercentage(0.1),
  fUseJets(kFALSE), fPtAssoc(0) {
  // Constructor
}

//________________________________________________________________________
AliAnalysisTaskToyModel::~AliAnalysisTaskToyModel() {
  //Destructor
  if(fUseAllCharges) {
    delete fPtSpectraAllCharges;
    delete fAzimuthalAngleAllCharges;
  }
  else {
    delete fPtSpectraPions;
    delete fAzimuthalAnglePions;
    delete fPtSpectraKaons;
    delete fAzimuthalAngleKaons;
    delete fPtSpectraProtons;
    delete fAzimuthalAngleProtons;
  }
  if(fUseJets) delete fPtAssoc;
}

//________________________________________________________________________
void AliAnalysisTaskToyModel::Init() {
  //Initialize objects
  //==============gRandom Seed=======================//
  gRandom->SetSeed(0);  //seed is set to a random value (depending on machine clock)
  //==============gRandom Seed=======================//

  //==============Particles and spectra==============//
  TParticle *pion = new TParticle();
  pion->SetPdgCode(211);
  fPionMass = pion->GetMass();

  TParticle *kaon = new TParticle();
  kaon->SetPdgCode(321);
  fKaonMass = kaon->GetMass();
  
  TParticle *proton = new TParticle();
  proton->SetPdgCode(2212);
  fProtonMass = proton->GetMass();

  if(fUseAllCharges) {
    fParticleMass = fPionMass;
    fPtSpectraAllCharges = new TF1("fPtSpectraAllCharges","x*TMath::Exp(-TMath::Power([0]*[0]+x*x,0.5)/[1])",0.,5.);
    fPtSpectraAllCharges->SetParName(0,"Mass");
    fPtSpectraAllCharges->SetParName(1,"Temperature");
    //fPtSpectraAllCharges = new TF1("fPtSpectraAllCharges","(x^2/TMath::Sqrt(TMath::Power(x,2) + TMath::Power(0.139,2)))*TMath::Power((1. + x/[0]),-[1])",0.,20.);
    //fPtSpectraAllCharges->SetParName(0,"pt0");
    //fPtSpectraAllCharges->SetParName(1,"b");
  }
  else {
    fPtSpectraPions = new TF1("fPtSpectraPions","x*TMath::Exp(-TMath::Power([0]*[0]+x*x,0.5)/[1])",0.,5.);
    fPtSpectraPions->SetParName(0,"Mass");
    fPtSpectraPions->SetParName(1,"Temperature");
    
    fPtSpectraKaons = new TF1("fPtSpectraKaons","x*TMath::Exp(-TMath::Power([0]*[0]+x*x,0.5)/[1])",0.,5.);
    fPtSpectraKaons->SetParName(0,"Mass");
    fPtSpectraKaons->SetParName(1,"Temperature");
    
    fPtSpectraProtons = new TF1("fPtSpectraProtons","x*TMath::Exp(-TMath::Power([0]*[0]+x*x,0.5)/[1])",0.,5.);
    fPtSpectraProtons->SetParName(0,"Mass");
    fPtSpectraProtons->SetParName(1,"Temperature");
  }
  //==============Particles and spectra==============//

  //==============Flow values==============//
  if(fUseAllCharges) {
    if(fUseDebug) {
      Printf("Directed flow: %lf",fDirectedFlowAllCharges);
      Printf("Elliptic flow: %lf",fEllipticFlowAllCharges);
      Printf("Triangular flow: %lf",fTriangularFlowAllCharges);
      Printf("Quandrangular flow: %lf",fQuandrangularFlowAllCharges);
      Printf("Pentangular flow: %lf",fPentangularFlowAllCharges);
    }

    fAzimuthalAngleAllCharges = new TF1("fAzimuthalAngleAllCharges","1+2.*[1]*TMath::Cos(x-[0])+2.*[2]*TMath::Cos(2*(x-[0]))+2.*[3]*TMath::Cos(3*(x-[0]))+2.*[4]*TMath::Cos(4*(x-[0]))+2.*[5]*TMath::Cos(5*(x-[0]))",0.,2.*TMath::Pi());
    fAzimuthalAngleAllCharges->SetParName(0,"Reaction Plane");
    fAzimuthalAngleAllCharges->SetParName(1,"Directed flow");
    fAzimuthalAngleAllCharges->SetParName(2,"Elliptic flow"); 
    fAzimuthalAngleAllCharges->SetParName(3,"Triangular flow");
    fAzimuthalAngleAllCharges->SetParName(4,"Quandrangular flow");
    fAzimuthalAngleAllCharges->SetParName(5,"Pentangular flow");
  }
  else {
    fAzimuthalAnglePions = new TF1("fAzimuthalAnglePions","1+2.*[1]*TMath::Cos(x-[0])+2.*[2]*TMath::Cos(2*(x-[0]))+2.*[3]*TMath::Cos(3*(x-[0]))+2.*[4]*TMath::Cos(4*(x-[0]))+2.*[5]*TMath::Cos(5*(x-[0]))",0.,2.*TMath::Pi());
    fAzimuthalAnglePions->SetParName(0,"Reaction Plane");
    fAzimuthalAnglePions->SetParName(1,"Directed flow");
    fAzimuthalAnglePions->SetParName(2,"Elliptic flow"); 
    fAzimuthalAnglePions->SetParName(3,"Triangular flow");
    fAzimuthalAnglePions->SetParName(4,"Quandrangular flow");
    fAzimuthalAnglePions->SetParName(5,"Pentangular flow");
    
    fAzimuthalAngleKaons = new TF1("fAzimuthalAngleKaons","1+2.*[1]*TMath::Cos(x-[0])+2.*[2]*TMath::Cos(2*(x-[0]))+2.*[3]*TMath::Cos(3*(x-[0]))+2.*[4]*TMath::Cos(4*(x-[0]))+2.*[5]*TMath::Cos(5*(x-[0]))",0.,2.*TMath::Pi());
    fAzimuthalAngleKaons->SetParName(0,"Reaction Plane");
    fAzimuthalAngleKaons->SetParName(1,"Directed flow");
    fAzimuthalAngleKaons->SetParName(2,"Elliptic flow"); 
    fAzimuthalAngleKaons->SetParName(3,"Triangular flow");
    fAzimuthalAngleKaons->SetParName(4,"Quandrangular flow");
    fAzimuthalAngleKaons->SetParName(5,"Pentangular flow");
    
    fAzimuthalAngleProtons = new TF1("fAzimuthalAngleProtons","1+2.*[1]*TMath::Cos(x-[0])+2.*[2]*TMath::Cos(2*(x-[0]))+2.*[3]*TMath::Cos(3*(x-[0]))+2.*[4]*TMath::Cos(4*(x-[0]))+2.*[5]*TMath::Cos(5*(x-[0]))",0.,2.*TMath::Pi());
    fAzimuthalAngleProtons->SetParName(0,"Reaction Plane");
    fAzimuthalAngleProtons->SetParName(1,"Directed flow");
    fAzimuthalAngleProtons->SetParName(2,"Elliptic flow"); 
    fAzimuthalAngleProtons->SetParName(3,"Triangular flow");
    fAzimuthalAngleProtons->SetParName(4,"Quandrangular flow");
    fAzimuthalAngleProtons->SetParName(5,"Pentangular flow");
  }
  //==============Flow values==============//

  //===================Jets===================//
  if(fUseJets) {
    fPtAssoc = new TF1("fPtAssoc","x*TMath::Exp(-TMath::Power([0]*[0]+x*x,0.5)/[1])",0.,20.);
    fPtAssoc->SetParName(0,"pt0");
    fPtAssoc->SetParName(1,"b");
    fPtAssoc->SetParameter(0,0.139);
    fPtAssoc->SetParameter(1,0.5);
    fPtAssoc->SetLineColor(1);
  }

  //===================Jets===================//

  //==============Efficiency matrix==============//
  if(fSimulateDetectorEffects) SetupEfficiencyMatrix();
  //==============Efficiency matrix==============//
}

//________________________________________________________________________
void AliAnalysisTaskToyModel::SetupEfficiencyMatrix() {
  //Setup the efficiency matrix
  TH1F *hPt = new TH1F("hPt","",200,0.1,20.1);
  TH1F *hEta = new TH1F("hEta","",20,-0.95,0.95);
  TH1F *hPhi = new TH1F("hPhi","",72,0.,2.*TMath::Pi());
  fEfficiencyMatrix = new TH3F("fEfficiencyMatrix","",
			       hEta->GetNbinsX(),
			       hEta->GetXaxis()->GetXmin(),
			       hEta->GetXaxis()->GetXmax(),
			       hPt->GetNbinsX(),
			       hPt->GetXaxis()->GetXmin(),
			       hPt->GetXaxis()->GetXmax(),
			       hPhi->GetNbinsX(),
			       hPhi->GetXaxis()->GetXmin(),
			       hPhi->GetXaxis()->GetXmax());

  //Efficiency in pt
  Double_t epsilon[20] = {0.3,0.6,0.77,0.79,0.80,0.80,0.80,0.80,0.80,0.80,0.80,0.80,0.80,0.80,0.80,0.80,0.80,0.80,0.80,0.80};
  for(Int_t i=1;i<=20;i++) {
    hPt->SetBinContent(i,epsilon[i-1]);
    hPt->SetBinError(i,0.01);
  }
  for(Int_t i=21;i<=200;i++) {
    hPt->SetBinContent(i,epsilon[19]);
    hPt->SetBinError(i,0.01);
  }

  //Efficiency in eta
  for(Int_t i=1;i<=hEta->GetNbinsX();i++) {
    hEta->SetBinContent(i,1.0);
    hEta->SetBinError(i,0.01);
  }
  if(fEfficiencyDropNearEtaEdges) {
    hEta->SetBinContent(1,0.7); hEta->SetBinContent(2,0.8);
    hEta->SetBinContent(3,0.9);
    hEta->SetBinContent(18,0.9); hEta->SetBinContent(19,0.8);
    hEta->SetBinContent(20,0.7);
  }

  //Efficiency in phi
  for(Int_t i=1;i<=hPhi->GetNbinsX();i++) {
    hPhi->SetBinContent(i,1.0);
    hPhi->SetBinError(i,0.01);
  }
  for(Int_t i=1;i<=fNumberOfInefficientSectors;i++)
    hPhi->SetBinContent(hPhi->FindBin(hPhi->GetRandom()),fInefficiencyFactorInPhi);
  for(Int_t i=1;i<=fNumberOfDeadSectors;i++)
    hPhi->SetBinContent(hPhi->FindBin(hPhi->GetRandom()),0.0);
  
  //Fill the 3D efficiency map
  for(Int_t iBinX = 1; iBinX<=fEfficiencyMatrix->GetXaxis()->GetNbins();iBinX++) {
    //cout<<"==================================="<<endl;
    for(Int_t iBinY = 1; iBinY<=fEfficiencyMatrix->GetYaxis()->GetNbins();iBinY++) {
      //cout<<"==================================="<<endl;
      for(Int_t iBinZ = 1; iBinZ<=fEfficiencyMatrix->GetZaxis()->GetNbins();iBinZ++) {
	fEfficiencyMatrix->SetBinContent(iBinX,iBinY,iBinZ,hEta->GetBinContent(iBinX)*hPt->GetBinContent(iBinY)*hPhi->GetBinContent(iBinZ));
	//cout<<"Eta: "<<hEta->GetBinCenter(iBinX)<<" - Pt: "<<hPt->GetBinCenter(iBinY)<<" - Phi: "<<hPhi->GetBinCenter(iBinZ)<<" - "<<hEta->GetBinContent(iBinX)<<" , "<<hPt->GetBinContent(iBinY)<<" , "<<hPhi->GetBinContent(iBinZ)<<" - Efficiency: "<<hEta->GetBinContent(iBinX)*hPt->GetBinContent(iBinY)*hPhi->GetBinContent(iBinZ)<<endl;
      }
    }
  }
}

//________________________________________________________________________
void AliAnalysisTaskToyModel::CreateOutputObjects() {
  // Create histograms
  // Called once

  // global switch disabling the reference 
  // (to avoid "Replacing existing TH1" if several wagons are created in train)
  Bool_t oldStatus = TH1::AddDirectoryStatus();
  TH1::AddDirectory(kFALSE);

  if(!fBalance) {
    fBalance = new AliBalancePsi();
    fBalance->SetDeltaEtaMax(2.0);
  }
  if(fRunShuffling) {
    if(!fShuffledBalance) {
      fShuffledBalance = new AliBalancePsi();
      fShuffledBalance->SetDeltaEtaMax(2.0);
    }
  }
  if(fRunMixing) {
    if(!fMixedBalance) {
      fMixedBalance = new AliBalancePsi();
      fMixedBalance->SetDeltaEtaMax(2.0);
    }
  }

  //QA list
  fList = new TList();
  fList->SetName("listQA");
  fList->SetOwner();

  //Balance Function list
  fListBF = new TList();
  fListBF->SetName("listBF");
  fListBF->SetOwner();

  if(fRunShuffling) {
    fListBFS = new TList();
    fListBFS->SetName("listBFShuffled");
    fListBFS->SetOwner();
  }

  if(fRunMixing) {
    fListBFM = new TList();
    fListBFM->SetName("listBFMixed");
    fListBFM->SetOwner();
  }

  //==============QA================//
  //Event stats.
  TString gCutName[4] = {"Total","Offline trigger",
                         "Vertex","Analyzed"};
  fHistEventStats = new TH1F("fHistEventStats",
                             "Event statistics;;N_{events}",
                             4,0.5,4.5);
  for(Int_t i = 1; i <= 4; i++)
    fHistEventStats->GetXaxis()->SetBinLabel(i,gCutName[i-1].Data());
  fList->Add(fHistEventStats);

  fHistNumberOfAcceptedParticles = new TH1F("fHistNumberOfAcceptedParticles",";N_{acc.};Entries",10000,0,10000);
  fList->Add(fHistNumberOfAcceptedParticles);
  
  fHistReactionPlane = new TH1F("fHistReactionPlane","Reaction plane angle;#Psi [rad];Entries",1000,0.,2.*TMath::Pi());
  fList->Add(fHistReactionPlane);

  //Pseudo-rapidity
  fHistEtaTotal = new TH1F("fHistEtaTotal","Pseudo-rapidity (full phase space);#eta;Entries",1000,-15.,15.); 
  fList->Add(fHistEtaTotal);

  fHistEta = new TH1F("fHistEta","Pseudo-rapidity (acceptance);#eta;Entries",1000,-1.5,1.5); 
  fList->Add(fHistEta);

  fHistEtaPhiPos = new TH2F("fHistEtaPhiPos","#eta-#phi distribution (+);#eta;#varphi (rad)",80,-2.,2.,72,-TMath::Pi()/2.,3.*TMath::Pi()/2.);
  fList->Add(fHistEtaPhiPos);
  fHistEtaPhiNeg = new TH2F("fHistEtaPhiNeg","#eta-#phi distribution (-);#eta;#varphi (rad)",80,-2.,2.,72,-TMath::Pi()/2.,3.*TMath::Pi()/2.);
  fList->Add(fHistEtaPhiNeg);

  //Rapidity
  fHistRapidity = new TH1F("fHistRapidity","Rapidity (acceptance);y;Entries",1000,-1.5,1.5); 
  fList->Add(fHistRapidity);
  fHistRapidityPions = new TH1F("fHistRapidityPions","Rapidity (acceptance - pions);y;Entries",1000,-1.5,1.5); 
  fList->Add(fHistRapidityPions);
  fHistRapidityKaons = new TH1F("fHistRapidityKaons","Rapidity (acceptance - kaons);y;Entries",1000,-1.5,1.5); 
  fList->Add(fHistRapidityKaons);
  fHistRapidityProtons = new TH1F("fHistRapidityProtons","Rapidity (acceptance - protons);y;Entries",1000,-1.5,1.5); 
  fList->Add(fHistRapidityProtons);

  //Phi
  fHistPhi = new TH1F("fHistPhi","Phi (acceptance);#phi (rad);Entries",1000,0.,2*TMath::Pi());
  fList->Add(fHistPhi);

  fHistPhiPions = new TH1F("fHistPhiPions","Phi (acceptance - pions);#phi (rad);Entries",1000,0.,2*TMath::Pi());
  fList->Add(fHistPhiPions);
  fHistPhiKaons = new TH1F("fHistPhiKaons","Phi (acceptance - kaons);#phi (rad);Entries",1000,0.,2*TMath::Pi());
  fList->Add(fHistPhiKaons);
  fHistPhiProtons = new TH1F("fHistPhiProtons","Phi (acceptance - protons);#phi (rad);Entries",1000,0.,2*TMath::Pi());
  fList->Add(fHistPhiProtons);

  //Pt
  fHistPt = new TH1F("fHistPt","Pt (acceptance);p_{t} (GeV/c);Entries",1000,0.,10.);
  fList->Add(fHistPt);

  fHistPtPions = new TH1F("fHistPtPions","Pt (acceptance - pions);p_{t} (GeV/c);Entries",1000,0.,10.);
  fList->Add(fHistPtPions);
  fHistPtKaons = new TH1F("fHistPtKaons","Pt (acceptance - kaons);p_{t} (GeV/c);Entries",1000,0.,10.);
  fList->Add(fHistPtKaons);
  fHistPtProtons = new TH1F("fHistPtProtons","Pt (acceptance - protons);p_{t} (GeV/c);Entries",1000,0.,10.);
  fList->Add(fHistPtProtons);

  if(fEfficiencyMatrix) fList->Add(fEfficiencyMatrix);

  //==============Balance function histograms================//
  // Initialize histograms if not done yet
  if(!fBalance->GetHistNp()){
    AliWarning("Histograms not yet initialized! --> Will be done now");
    AliWarning("--> Add 'gBalance->InitHistograms()' in your configBalanceFunction");
    fBalance->InitHistograms();
  }

  if(fRunShuffling) {
    if(!fShuffledBalance->GetHistNp()) {
      AliWarning("Histograms (shuffling) not yet initialized! --> Will be done now");
      AliWarning("--> Add 'gBalance->InitHistograms()' in your configBalanceFunction");
      fShuffledBalance->InitHistograms();
    }
  }

  fListBF->Add(fBalance->GetHistNp());
  fListBF->Add(fBalance->GetHistNn());
  fListBF->Add(fBalance->GetHistNpn());
  fListBF->Add(fBalance->GetHistNnn());
  fListBF->Add(fBalance->GetHistNpp());
  fListBF->Add(fBalance->GetHistNnp());

  if(fRunShuffling) {
    fListBFS->Add(fShuffledBalance->GetHistNp());
    fListBFS->Add(fShuffledBalance->GetHistNn());
    fListBFS->Add(fShuffledBalance->GetHistNpn());
    fListBFS->Add(fShuffledBalance->GetHistNnn());
    fListBFS->Add(fShuffledBalance->GetHistNpp());
    fListBFS->Add(fShuffledBalance->GetHistNnp());
  }

  if(fRunMixing) {
    fListBFM->Add(fMixedBalance->GetHistNp());
    fListBFM->Add(fMixedBalance->GetHistNn());
    fListBFM->Add(fMixedBalance->GetHistNpn());
    fListBFM->Add(fMixedBalance->GetHistNnn());
    fListBFM->Add(fMixedBalance->GetHistNpp());
    fListBFM->Add(fMixedBalance->GetHistNnp());
  }

  // Event Mixing
  if(fRunMixing){
    Int_t fMixingTracks = 2000;
    Int_t trackDepth = fMixingTracks; 
    Int_t poolsize   = 1000;  // Maximum number of events, ignored in the present implemented of AliEventPoolManager
    
    // centrality bins
    // Double_t centralityBins[] = {0.,1.,2.,3.,4.,5.,7.,10.,20.,30.,40.,50.,60.,70.,80.,100.}; // SHOULD BE DEDUCED FROM CREATED ALITHN!!!
    // Double_t* centbins        = centralityBins;
    // Int_t nCentralityBins     = sizeof(centralityBins) / sizeof(Double_t) - 1;

    // multiplicity bins
    Double_t multiplicityBins[] = {0,10,20,30,40,50,60,70,80,100,100000}; // SHOULD BE DEDUCED FROM CREATED ALITHN!!!
    Double_t* multbins        = multiplicityBins;
    Int_t nMultiplicityBins     = sizeof(multiplicityBins) / sizeof(Double_t) - 1;
    
    // Zvtx bins
    Double_t vertexBins[] = {-10., -7., -5., -3., -1., 1., 3., 5., 7., 10.}; // SHOULD BE DEDUCED FROM CREATED ALITHN!!!
    Double_t* vtxbins     = vertexBins;
    Int_t nVertexBins     = sizeof(vertexBins) / sizeof(Double_t) - 1;
    
    // Event plane angle (Psi) bins
    // Double_t psiBins[] = {0.,45.,135.,215.,305.,360.}; // SHOULD BE DEDUCED FROM CREATED ALITHN!!!
    // Double_t* psibins     = psiBins;
    // Int_t nPsiBins     = sizeof(psiBins) / sizeof(Double_t) - 1;
    
    // // run the event mixing also in bins of event plane (statistics!)
    // if(fRunMixingEventPlane){
    //   if(fEventClass=="Multiplicity"){
    // 	fPoolMgr = new AliEventPoolManager(poolsize, trackDepth, nMultiplicityBins, multbins, nVertexBins, vtxbins, nPsiBins, psibins);
    //   }
    //   else{
    // 	fPoolMgr = new AliEventPoolManager(poolsize, trackDepth, nCentralityBins, centbins, nVertexBins, vtxbins, nPsiBins, psibins);
    //   }
    // }
    // else{
    //if(fEventClass=="Multiplicity"){
    fPoolMgr = new AliEventPoolManager(poolsize, trackDepth, nMultiplicityBins, multbins, nVertexBins, vtxbins);
    //}
    //else{
    //fPoolMgr = new AliEventPoolManager(poolsize, trackDepth, nCentralityBins, centbins, nVertexBins, vtxbins);
    //}
  }

  TH1::AddDirectory(oldStatus);
}

//________________________________________________________________________
void AliAnalysisTaskToyModel::Run(Int_t nEvents) {
  // Main loop
  // Called for each event
  Short_t vCharge = 0;
  Float_t vY = 0.0;
  Float_t vEta = 0.0;
  Float_t vPhi = 0.0;
  Float_t vP[3] = {0.,0.,0.};
  Float_t vPt = 0.0;
  Float_t vE = 0.0;
  Bool_t isPion = kFALSE, isKaon = kFALSE, isProton = kFALSE;

  Double_t gDecideCharge = 0.;

  if(fUseAllCharges) {
    //fPtSpectraAllCharges->SetParameter(0,fParticleMass);
    fPtSpectraAllCharges->SetParameter(0,1.05);
    fPtSpectraAllCharges->SetParameter(1,fTemperatureAllCharges);

    fAzimuthalAngleAllCharges->SetParameter(1,fDirectedFlowAllCharges);
    fAzimuthalAngleAllCharges->SetParameter(2,fEllipticFlowAllCharges);
    fAzimuthalAngleAllCharges->SetParameter(3,fTriangularFlowAllCharges);
    fAzimuthalAngleAllCharges->SetParameter(4,fQuandrangularFlowAllCharges);
    fAzimuthalAngleAllCharges->SetParameter(5,fPentangularFlowAllCharges);
  }
  else {
    fPtSpectraPions->SetParameter(0,fPionMass);
    fPtSpectraPions->SetParameter(1,fTemperaturePions);
    fPtSpectraKaons->SetParameter(0,fKaonMass);
    fPtSpectraKaons->SetParameter(1,fTemperatureKaons);
    fPtSpectraProtons->SetParameter(0,fProtonMass);
    fPtSpectraProtons->SetParameter(1,fTemperatureProtons);

    fAzimuthalAnglePions->SetParameter(1,fDirectedFlowPions);
    fAzimuthalAnglePions->SetParameter(2,fEllipticFlowPions);
    fAzimuthalAnglePions->SetParameter(3,fTriangularFlowPions);
    fAzimuthalAnglePions->SetParameter(4,fQuandrangularFlowPions);
    fAzimuthalAnglePions->SetParameter(5,fPentangularFlowPions);

    fAzimuthalAngleKaons->SetParameter(1,fDirectedFlowKaons);
    fAzimuthalAngleKaons->SetParameter(2,fEllipticFlowKaons);
    fAzimuthalAngleKaons->SetParameter(3,fTriangularFlowKaons);
    fAzimuthalAngleKaons->SetParameter(4,fQuandrangularFlowKaons);
    fAzimuthalAngleKaons->SetParameter(5,fPentangularFlowKaons);

    fAzimuthalAngleProtons->SetParameter(1,fDirectedFlowProtons);
    fAzimuthalAngleProtons->SetParameter(2,fEllipticFlowProtons);
    fAzimuthalAngleProtons->SetParameter(3,fTriangularFlowProtons);
    fAzimuthalAngleProtons->SetParameter(4,fQuandrangularFlowProtons);
    fAzimuthalAngleProtons->SetParameter(5,fPentangularFlowProtons);
  }


  for(Int_t iEvent = 0; iEvent < nEvents; iEvent++) {
    // vector holding the charges/kinematics of all tracks (charge,y,eta,phi,p0,p1,p2,pt,E)
    //vector<Double_t> *chargeVectorShuffle[9];   // this will be shuffled
    //vector<Double_t> *chargeVector[9];          // original charge
    //for(Int_t i = 0; i < 9; i++){
    //chargeVectorShuffle[i] = new vector<Double_t>;
    //chargeVector[i]        = new vector<Double_t>;
    //}

    // TObjArray for the accepted particles 
    // (has to be done here, otherwise mxing with event pool does not work, overwriting pointers!)
    TObjArray *tracksMain = new TObjArray();
    tracksMain->SetOwner(kTRUE);
    TObjArray *tracksMixing = 0x0;
    if(fRunMixing) {
      tracksMixing = new TObjArray();
      tracksMixing->SetOwner(kTRUE);
    }

    tracksMain->Clear();
    if(fRunMixing) tracksMixing->Clear();
    
    fHistEventStats->Fill(1);
    fHistEventStats->Fill(2);
    fHistEventStats->Fill(3);

    if((iEvent%1000) == 0) 
      cout<<"Event: "<<iEvent<<"/"<<nEvents<<endl;
    
    //Multiplicities
    Int_t nMultiplicity = (Int_t)(gRandom->Gaus(fTotalMultiplicityMean,fTotalMultiplicitySigma));
    Int_t nNetCharge = (Int_t)(gRandom->Gaus(fNetChargeMean,fNetChargeSigma));
    
    Int_t nGeneratedPositive = (Int_t)((nMultiplicity/2) + nNetCharge);
    Int_t nGeneratedNegative = nMultiplicity - nGeneratedPositive;
    if(fUseDebug) 
      Printf("Total multiplicity: %d - Generated positive: %d - Generated negative: %d",nMultiplicity,nGeneratedPositive,nGeneratedNegative);
    
    //Randomization of the reaction plane
    fReactionPlane = 2.0*TMath::Pi()*gRandom->Rndm();
    //fReactionPlane = 0.0;
    if(fUseAllCharges) 
      fAzimuthalAngleAllCharges->SetParameter(0,fReactionPlane);
    else {
      fAzimuthalAnglePions->SetParameter(0,fReactionPlane);
      fAzimuthalAngleKaons->SetParameter(0,fReactionPlane);
      fAzimuthalAngleProtons->SetParameter(0,fReactionPlane);
    }
    
    Int_t gNumberOfAcceptedParticles = 0;
    Int_t gNumberOfAcceptedPositiveParticles = 0;
    Int_t gNumberOfAcceptedNegativeParticles = 0;
    Int_t nGeneratedPions = 0, nGeneratedKaons = 0, nGeneratedProtons = 0;
 
    //Generate particles
    for(Int_t iParticleCount = 0; iParticleCount < nMultiplicity; iParticleCount++) {
      isPion = kFALSE; isKaon = kFALSE; isProton = kFALSE;
      if(fUseDebug) 
	Printf("Generating positive: %d(%d)",iParticleCount+1,nGeneratedPositive);

      //Pseudo-rapidity sampled from a Gaussian centered @ 0
      vEta = gRandom->Gaus(0.0,4.0);

      //Fill QA histograms (full phase space)
      fHistEtaTotal->Fill(vEta);

      //Decide the charge
      gDecideCharge = gRandom->Rndm();
      if(gDecideCharge <= 1.*nGeneratedPositive/nMultiplicity)       
	vCharge = 1;
      else 
	vCharge = -1;
      
      //Acceptance
      if((vEta < fEtaMin) || (vEta > fEtaMax)) continue;

      if(!fUseAllCharges) {
	//Decide the specie
	Double_t randomNumberSpecies = gRandom->Rndm();
	if((randomNumberSpecies >= 0.0)&&(randomNumberSpecies < fPionPercentage)) {
	  nGeneratedPions += 1;
	  vPt = fPtSpectraPions->GetRandom();
	  vPhi = fAzimuthalAnglePions->GetRandom();
	  fParticleMass = fPionMass;
	  isPion = kTRUE;
	}
	else if((randomNumberSpecies >= fPionPercentage)&&(randomNumberSpecies < fPionPercentage + fKaonPercentage)) {
	  nGeneratedKaons += 1;
	  vPt = fPtSpectraKaons->GetRandom();
	  vPhi = fAzimuthalAngleKaons->GetRandom();
	  fParticleMass = fKaonMass;
	  isKaon = kTRUE;
	}
	else if((randomNumberSpecies >= fPionPercentage + fKaonPercentage)&&(randomNumberSpecies < fPionPercentage + fKaonPercentage + fProtonPercentage)) {
	  nGeneratedProtons += 1;
	  vPt = fPtSpectraProtons->GetRandom();
	  vPhi = fAzimuthalAngleProtons->GetRandom();
	  fParticleMass = fProtonMass;
	  isProton = kTRUE;
	}
      }
      else {
	vPt = fPtSpectraAllCharges->GetRandom();
	vPhi = fAzimuthalAngleAllCharges->GetRandom();
      }
      
      vP[0] = vPt*TMath::Cos(vPhi);
      vP[1] = vPt*TMath::Sin(vPhi);
      vP[2] = vPt*TMath::SinH(vEta);
      vE = TMath::Sqrt(TMath::Power(fParticleMass,2) +
			TMath::Power(vP[0],2) +
			TMath::Power(vP[1],2) +
			TMath::Power(vP[2],2));
      
      vY = 0.5*TMath::Log((vE + vP[2])/(vE - vP[2]));
      
      //pt coverage
      if((vPt < fPtMin) || (vPt > fPtMax)) continue;
      //Printf("pt: %lf - mins: %lf - max: %lf",vPt,fPtMin,fPtMax);

      //acceptance filter
      if(fUseAcceptanceParameterization) {
	Double_t gRandomNumberForAcceptance = gRandom->Rndm();
	if(gRandomNumberForAcceptance > fAcceptanceParameterization->Eval(vPt)) 
	  continue;
      }

      //Detector effects
      if(fSimulateDetectorEffects) {
	Double_t randomNumber = gRandom->Rndm();
	if(randomNumber > fEfficiencyMatrix->GetBinContent(fEfficiencyMatrix->FindBin(vEta,vPt,vPhi)))
	  continue;
      }

      //Fill QA histograms (acceptance);
      if(vCharge > 0) {
	gNumberOfAcceptedPositiveParticles += 1;
	if(vPhi > 3.*TMath::Pi()/2.)
	  fHistEtaPhiPos->Fill(vEta,vPhi-2.*TMath::Pi());
	else
	  fHistEtaPhiPos->Fill(vEta,vPhi);
      }
      else {
	gNumberOfAcceptedNegativeParticles += 1;
	if(vPhi > 3.*TMath::Pi()/2.)
	  fHistEtaPhiNeg->Fill(vEta,vPhi-2.*TMath::Pi());
	else
	  fHistEtaPhiNeg->Fill(vEta,vPhi);
      }

      fHistEta->Fill(vEta);
      fHistRapidity->Fill(vY);
      fHistPhi->Fill(vPhi);
      fHistPt->Fill(vPt);
      if(isPion) {
	fHistRapidityPions->Fill(vY);
	fHistPhiPions->Fill(vPhi);
	fHistPtPions->Fill(vPt);
      }
      else if(isKaon) {
	fHistRapidityKaons->Fill(vY);
	fHistPhiKaons->Fill(vPhi);
	fHistPtKaons->Fill(vPt);
      }
      else if(isProton) {
	fHistRapidityProtons->Fill(vY);
	fHistPhiProtons->Fill(vPhi);
	fHistPtProtons->Fill(vPt);
      }
      gNumberOfAcceptedParticles += 1;

      // add the track to the TObjArray
      tracksMain->Add(new AliBFBasicParticle(vEta, vPhi, vPt, vCharge, 1.0));  
      if(fRunMixing) 
	tracksMixing->Add(new AliBFBasicParticle(vEta, vPhi, vPt, vCharge, 1.0));  
    }//generated positive particle loop

    //Jets
    if(fUseJets) {
      const Int_t nAssociated = 3;

      Double_t gPtTrig1 = 0., gPtTrig2 = 0., gPtAssoc = 0.;
      Double_t gPhiTrig1 = 0., gPhiTrig2 = 0., gPhiAssoc = 0.;
      Double_t gEtaTrig1 = 0., gEtaTrig2 = 0., gEtaAssoc = 0.;
      Short_t gChargeTrig1 = 0, gChargeTrig2 = 0, gChargeAssoc = 0;
 
      Double_t gJetCone = 0.2;

      //First leading particle
      gPtTrig1 = gRandom->Uniform(3.,5.);
      gEtaTrig1 = gRandom->Uniform(-0.8,0.8);
      gPhiTrig1 = gRandom->Uniform(0.,TMath::TwoPi());

      //Decide the charge
      gDecideCharge = gRandom->Rndm();
      if(gDecideCharge <= 0.5)
	gChargeTrig1 = 1;
      else 
	gChargeTrig1 = -1;
      
      //Acceptance
      if((gEtaTrig1 < fEtaMin) || (gEtaTrig1 > fEtaMax)) continue;
      //pt coverage
      if((gPtTrig1 < fPtMin) || (gPtTrig1 > fPtMax)) continue;
      //Printf("pt: %lf - mins: %lf - max: %lf",vPt,fPtMin,fPtMax);

      //acceptance filter
      if(fUseAcceptanceParameterization) {
	Double_t gRandomNumberForAcceptance = gRandom->Rndm();
	if(gRandomNumberForAcceptance > fAcceptanceParameterization->Eval(gPtTrig1)) 
	  continue;
      }

      //Detector effects
      if(fSimulateDetectorEffects) {
	Double_t randomNumber = gRandom->Rndm();
	if(randomNumber > fEfficiencyMatrix->GetBinContent(fEfficiencyMatrix->FindBin(gEtaTrig1,gPtTrig1,gPhiTrig1)))
	  continue;
      }

      gNumberOfAcceptedParticles += 1;

      // add the track to the TObjArray
      tracksMain->Add(new AliBFBasicParticle(gEtaTrig1, gPhiTrig1, gPtTrig1, gChargeTrig1, 1.0));  
      if(fRunMixing) 
	tracksMixing->Add(new AliBFBasicParticle(gEtaTrig1, gPhiTrig1, gPtTrig1, gChargeTrig1, 1.0));  

      Int_t iAssociated = 0; 
      while(iAssociated < nAssociated) {
	gPtAssoc = fPtAssoc->GetRandom();
	if(gPtAssoc < gPtTrig1) {
	  gEtaAssoc = gRandom->Uniform(gEtaTrig1 - gJetCone/2.,gEtaTrig1 + gJetCone/2.);
	  gPhiAssoc = gRandom->Uniform(gPhiTrig1 - gJetCone/2.,gPhiTrig1 + gJetCone/2.);
	  if(gPhiAssoc < 0.) gPhiAssoc += TMath::TwoPi();
	  else if(gPhiAssoc > TMath::TwoPi()) gPhiAssoc -= TMath::TwoPi();
	  
	  iAssociated += 1;

	  gDecideCharge = gRandom->Rndm();
	  if(gDecideCharge <= 0.5)
	    gChargeAssoc = 1;
	  else 
	    gChargeAssoc = -1;
	  
	  //Acceptance
	  if((gEtaAssoc < fEtaMin) || (gEtaAssoc > fEtaMax)) continue;
	  //pt coverage
	  if((gPtAssoc < fPtMin) || (gPtAssoc > fPtMax)) continue;
	  //Printf("pt: %lf - mins: %lf - max: %lf",vPt,fPtMin,fPtMax);

	  //acceptance filter
	  if(fUseAcceptanceParameterization) {
	    Double_t gRandomNumberForAcceptance = gRandom->Rndm();
	    if(gRandomNumberForAcceptance > fAcceptanceParameterization->Eval(gPtAssoc)) 
	      continue;
	  }

	  //Detector effects
	  if(fSimulateDetectorEffects) {
	    Double_t randomNumber = gRandom->Rndm();
	    if(randomNumber > fEfficiencyMatrix->GetBinContent(fEfficiencyMatrix->FindBin(gEtaAssoc,gPtAssoc,gPhiAssoc)))
	      continue;
	  }

	  gNumberOfAcceptedParticles += 1;

	  // add the track to the TObjArray
	  tracksMain->Add(new AliBFBasicParticle(gEtaAssoc, gPhiAssoc, gPtAssoc, gChargeAssoc, 1.0));  
	  if(fRunMixing) 
	    tracksMixing->Add(new AliBFBasicParticle(gEtaAssoc, gPhiAssoc, gPtAssoc, gChargeAssoc, 1.0));  
	}//pt,assoc < pt,trig
      }//associated
      
      //back2back
      gPtTrig2 = gPtTrig1;
      gEtaTrig2 = -gEtaTrig1;
      gPhiTrig2 = TMath::Pi() + gPhiTrig1;
      if(gPhiTrig2 < 0.) gPhiTrig2 += TMath::TwoPi();
      else if(gPhiTrig2 > TMath::TwoPi()) gPhiTrig2 -= TMath::TwoPi();

      //Decide the charge
      gDecideCharge = gRandom->Rndm();
      if(gDecideCharge <= 0.5)
	gChargeTrig2 = 1;
      else 
	gChargeTrig2 = -1;
      
      //Acceptance
      if((gEtaTrig2 < fEtaMin) || (gEtaTrig2 > fEtaMax)) continue;
      //pt coverage
      if((gPtTrig2 < fPtMin) || (gPtTrig2 > fPtMax)) continue;
      //Printf("pt: %lf - mins: %lf - max: %lf",vPt,fPtMin,fPtMax);

      //acceptance filter
      if(fUseAcceptanceParameterization) {
	Double_t gRandomNumberForAcceptance = gRandom->Rndm();
	if(gRandomNumberForAcceptance > fAcceptanceParameterization->Eval(gPtTrig2)) 
	  continue;
      }

      //Detector effects
      if(fSimulateDetectorEffects) {
	Double_t randomNumber = gRandom->Rndm();
	if(randomNumber > fEfficiencyMatrix->GetBinContent(fEfficiencyMatrix->FindBin(gEtaTrig2,gPtTrig2,gPhiTrig2)))
	  continue;
      }

      gNumberOfAcceptedParticles += 1;

      // add the track to the TObjArray
      tracksMain->Add(new AliBFBasicParticle(gEtaTrig2, gPhiTrig2, gPtTrig2, gChargeTrig2, 1.0));  
      if(fRunMixing) 
	tracksMixing->Add(new AliBFBasicParticle(gEtaTrig2, gPhiTrig2, gPtTrig2, gChargeTrig2, 1.0));  

      iAssociated = 0; 
      while(iAssociated < nAssociated) {
	gPtAssoc = fPtAssoc->GetRandom();
	if(gPtAssoc < gPtTrig2) {
	  gEtaAssoc = gRandom->Uniform(gEtaTrig2 - gJetCone/2.,gEtaTrig2 + gJetCone/2.);
	  gPhiAssoc = gRandom->Uniform(gPhiTrig2 - gJetCone/2.,gPhiTrig2 + gJetCone/2.);
	  if(gPhiAssoc < 0.) gPhiAssoc += TMath::TwoPi();
	  else if(gPhiAssoc > TMath::TwoPi()) gPhiAssoc -= TMath::TwoPi();
	  
	  iAssociated += 1;

	  gDecideCharge = gRandom->Rndm();
	  if(gDecideCharge <= 0.5)
	    gChargeAssoc = 1;
	  else 
	    gChargeAssoc = -1;
	  
	  //Acceptance
	  if((gEtaAssoc < fEtaMin) || (gEtaAssoc > fEtaMax)) continue;
	  //pt coverage
	  if((gPtAssoc < fPtMin) || (gPtAssoc > fPtMax)) continue;
	  //Printf("pt: %lf - mins: %lf - max: %lf",vPt,fPtMin,fPtMax);

	  //acceptance filter
	  if(fUseAcceptanceParameterization) {
	    Double_t gRandomNumberForAcceptance = gRandom->Rndm();
	    if(gRandomNumberForAcceptance > fAcceptanceParameterization->Eval(gPtAssoc)) 
	      continue;
	  }

	  //Detector effects
	  if(fSimulateDetectorEffects) {
	    Double_t randomNumber = gRandom->Rndm();
	    if(randomNumber > fEfficiencyMatrix->GetBinContent(fEfficiencyMatrix->FindBin(gEtaAssoc,gPtAssoc,gPhiAssoc)))
	      continue;
	  }

	  gNumberOfAcceptedParticles += 1;

	  // add the track to the TObjArray
	  tracksMain->Add(new AliBFBasicParticle(gEtaAssoc, gPhiAssoc, gPtAssoc, gChargeAssoc, 1.0));  
	  if(fRunMixing) 
	    tracksMixing->Add(new AliBFBasicParticle(gEtaAssoc, gPhiAssoc, gPtAssoc, gChargeAssoc, 1.0));  
	}//pt,assoc < pt,trig
      }//associated
    }//Jet usage








    
    //Dynamical correlations
    Int_t nGeneratedPositiveDynamicalCorrelations = 0;
    Int_t nGeneratedNegativeDynamicalCorrelations = 0;
    /*Double_t vChargePrime = 0;
    Double_t vYPrime = 0.0;
    Double_t vEtaPrime = 0.0;
    Double_t vPhiPrime = 0.0;
    Double_t vPPrime[3] = {0.,0.,0.};
    Double_t vPtPrime = 0.0;
    Double_t vEPrime = 0.0;
    //Generate "correlated" particles 
    if(fUseDynamicalCorrelations) {
      Int_t gNumberOfDynamicalCorrelations = (Int_t)(0.5*gNumberOfAcceptedParticles*fDynamicalCorrelationsPercentage);
      for(Int_t iDynamicalCorrelations = 0; iDynamicalCorrelations < gNumberOfDynamicalCorrelations; iDynamicalCorrelations++) {
	isPion = kFALSE; isKaon = kFALSE; isProton = kFALSE;
	
	//Pseudo-rapidity sampled from a Gaussian centered @ 0
	vEta = gRandom->Gaus(0.0,0.1);
	vCharge = 1.0;
	nGeneratedPositiveDynamicalCorrelations += 1;
	
	vEtaPrime = -vEta;
	vChargePrime = -1.0;
	nGeneratedNegativeDynamicalCorrelations += 1;
	  
	//Acceptance
	if((vEta < fEtaMin) || (vEta > fEtaMax)) continue;
	if((vEtaPrime < fEtaMin) || (vEtaPrime > fEtaMax)) continue;
      
	if(!fUseAllCharges) {
	  //Decide the specie
	  Double_t randomNumberSpecies = gRandom->Rndm();
	  if((randomNumberSpecies >= 0.0)&&(randomNumberSpecies < fPionPercentage)) {
	    nGeneratedPions += 1;
	    vPt = fPtSpectraPions->GetRandom();
	    vPhi = fAzimuthalAnglePions->GetRandom();
	    fParticleMass = fPionMass;
	    isPion = kTRUE;
	  }
	  else if((randomNumberSpecies >= fPionPercentage)&&(randomNumberSpecies < fPionPercentage + fKaonPercentage)) {
	    nGeneratedKaons += 1;
	    vPt = fPtSpectraKaons->GetRandom();
	    vPhi = fAzimuthalAngleKaons->GetRandom();
	    fParticleMass = fKaonMass;
	    isKaon = kTRUE;
	  }
	  else if((randomNumberSpecies >= fPionPercentage + fKaonPercentage)&&(randomNumberSpecies < fPionPercentage + fKaonPercentage + fProtonPercentage)) {
	    nGeneratedProtons += 1;
	    vPt = fPtSpectraProtons->GetRandom();
	    vPtPrime = vPt;
	    vPhi = fAzimuthalAngleProtons->GetRandom();
	    fParticleMass = fProtonMass;
	    isProton = kTRUE;
	  }
	}
	else {
	  vPt = fPtSpectraAllCharges->GetRandom();
	  vPhi = fAzimuthalAngleAllCharges->GetRandom();
	}
	vPtPrime = vPt;
	vPhiPrime = vPhi;

	vP[0] = vPt*TMath::Cos(vPhi);
	vP[1] = vPt*TMath::Sin(vPhi);
	vP[2] = vPt*TMath::SinH(vEta);
	vE = TMath::Sqrt(TMath::Power(fParticleMass,2) +
			  TMath::Power(vP[0],2) +
			  TMath::Power(vP[1],2) +
			  TMath::Power(vP[2],2));
	
	vY = 0.5*TMath::Log((vE + vP[2])/(vE - vP[2]));

	vPPrime[0] = vPtPrime*TMath::Cos(vPhiPrime);
	vPPrime[1] = vPtPrime*TMath::Sin(vPhiPrime);
	vPPrime[2] = vPtPrime*TMath::SinH(vEtaPrime);
	vEPrime = TMath::Sqrt(TMath::Power(fParticleMass,2) +
			  TMath::Power(vPPrime[0],2) +
			  TMath::Power(vPPrime[1],2) +
			  TMath::Power(vPPrime[2],2));
	
	vYPrime = 0.5*TMath::Log((vEPrime + vPPrime[2])/(vEPrime - vPPrime[2]));
      
	//pt coverage
	if((vPt < fPtMin) || (vPt > fPtMax)) continue;
	if((vPtPrime < fPtMin) || (vPtPrime > fPtMax)) continue;

	//acceptance filter
	if(fUseAcceptanceParameterization) {
	  Double_t gRandomNumberForAcceptance = gRandom->Rndm();
	  if(gRandomNumberForAcceptance > fAcceptanceParameterization->Eval(vPt)) 
	    continue;
	  
	  Double_t gRandomNumberForAcceptancePrime = gRandom->Rndm();
	  if(gRandomNumberForAcceptancePrime > fAcceptanceParameterization->Eval(vPtPrime)) 
	    continue;
	}
      
	// fill charge vector (positive)
	chargeVector[0]->push_back(vCharge);
	chargeVector[1]->push_back(vY);
	chargeVector[2]->push_back(vEta);
	chargeVector[3]->push_back(vPhi);
	chargeVector[4]->push_back(vP[0]);
	chargeVector[5]->push_back(vP[1]);
	chargeVector[6]->push_back(vP[2]);
	chargeVector[7]->push_back(vPt);
	chargeVector[8]->push_back(vE);
	
	if(fRunShuffling) {
	  chargeVectorShuffle[0]->push_back(vCharge);
	  chargeVectorShuffle[1]->push_back(vY);
	  chargeVectorShuffle[2]->push_back(vEta);
	  chargeVectorShuffle[3]->push_back(vPhi);
	  chargeVectorShuffle[4]->push_back(vP[0]);
	  chargeVectorShuffle[5]->push_back(vP[1]);
	  chargeVectorShuffle[6]->push_back(vP[2]);
	  chargeVectorShuffle[7]->push_back(vPt);
	  chargeVectorShuffle[8]->push_back(vE);
	}

	// fill charge vector (negative)
	chargeVector[0]->push_back(vChargePrime);
	chargeVector[1]->push_back(vYPrime);
	chargeVector[2]->push_back(vEtaPrime);
	chargeVector[3]->push_back(vPhiPrime);
	chargeVector[4]->push_back(vPPrime[0]);
	chargeVector[5]->push_back(vPPrime[1]);
	chargeVector[6]->push_back(vPPrime[2]);
	chargeVector[7]->push_back(vPtPrime);
	chargeVector[8]->push_back(vEPrime);
	
	if(fRunShuffling) {
	  chargeVectorShuffle[0]->push_back(vChargePrime);
	  chargeVectorShuffle[1]->push_back(vYPrime);
	  chargeVectorShuffle[2]->push_back(vEtaPrime);
	  chargeVectorShuffle[3]->push_back(vPhiPrime);
	  chargeVectorShuffle[4]->push_back(vPPrime[0]);
	  chargeVectorShuffle[5]->push_back(vPPrime[1]);
	  chargeVectorShuffle[6]->push_back(vPPrime[2]);
	  chargeVectorShuffle[7]->push_back(vPtPrime);
	  chargeVectorShuffle[8]->push_back(vEPrime);
	}

	gNumberOfAcceptedParticles += 2;
	}//loop over the dynamical correlations
	}*/// usage of dynamical correlations

    if(fUseDebug) {
      Printf("=======================================================");
      Printf("Total: %d - Total positive: %d - Total negative: %d",nMultiplicity,nGeneratedPositive,nGeneratedNegative);
      Printf("Accepted positive: %d - Accepted negative: %d",gNumberOfAcceptedPositiveParticles,gNumberOfAcceptedNegativeParticles);
      Printf("Correlations: %d - Correlations positive: %d - Correlations negative: %d",nGeneratedPositiveDynamicalCorrelations+nGeneratedNegativeDynamicalCorrelations,nGeneratedPositiveDynamicalCorrelations,nGeneratedNegativeDynamicalCorrelations);
      Printf("Number of accepted particles: %d",gNumberOfAcceptedParticles);
      if(!fUseAllCharges)
	Printf("Pions: %lf - Kaons: %lf - Protons: %lf",1.*nGeneratedPions/nMultiplicity,1.*nGeneratedKaons/nMultiplicity,1.*nGeneratedProtons/nMultiplicity);
      //Printf("Calculating the balance function for %d particles",chargeVector[0]->size());
    }

    fHistEventStats->Fill(4);
    fHistNumberOfAcceptedParticles->Fill(gNumberOfAcceptedParticles);
    fHistReactionPlane->Fill(fReactionPlane);

    // Event mixing 
    if (fRunMixing)
      {
	// 1. First get an event pool corresponding in mult (cent) and
	//    zvertex to the current event. Once initialized, the pool
	//    should contain nMix (reduced) events. This routine does not
	//    pre-scan the chain. The first several events of every chain
	//    will be skipped until the needed pools are filled to the
	//    specified depth. If the pool categories are not too rare, this
	//    should not be a problem. If they are rare, you could lose`
	//    statistics.
	
	// 2. Collect the whole pool's content of tracks into one TObjArray
	//    (bgTracks), which is effectively a single background super-event.
	
	// 3. The reduced and bgTracks arrays must both be passed into
	//    FillCorrelations(). Also nMix should be passed in, so a weight
	//    of 1./nMix can be applied.
	
	AliEventPool* pool = fPoolMgr->GetEventPool(1., 0.,fReactionPlane);
      
	if (!pool){
	  AliFatal(Form("No pool found for centrality = %f, zVtx = %f, psi = %f", 1., 0.,fReactionPlane));
	}
	else{
	  
	  //pool->SetDebug(1);
	  
	  if (pool->IsReady() || pool->NTracksInPool() > 50000 / 10 || pool->GetCurrentNEvents() >= 5){ 
	    
	    
	    Int_t nMix = pool->GetCurrentNEvents();
	    //cout << "nMix = " << nMix << " tracks in pool = " << pool->NTracksInPool() << ", tracks in this event = "<<gNumberOfAcceptedParticles <<", event Plane = "<<fReactionPlane<<endl;
	    
	    //((TH1F*) fListOfHistos->FindObject("eventStat"))->Fill(2);
	    //((TH2F*) fListOfHistos->FindObject("mixedDist"))->Fill(centrality, pool->NTracksInPool());
	    //if (pool->IsReady())
	    //((TH1F*) fListOfHistos->FindObject("eventStat"))->Fill(3);
	    
	    // Fill mixed-event histos here  
	    for (Int_t jMix=0; jMix<nMix; jMix++) 
	      {
		TObjArray* tracksMixed = pool->GetEvent(jMix);
		fMixedBalance->CalculateBalance(fReactionPlane,tracksMain,tracksMixed,1,1.,0.);
	      }
	  }
	  
	  // Update the Event pool
	  pool->UpdatePool(tracksMain);
	  //pool->PrintInfo();
	  
	}//pool NULL check  
      }//run mixing
    
    //Calculate the balance function
    fBalance->CalculateBalance(fReactionPlane,tracksMain,NULL,1,1.,0.);
    //if(fRunMixing)
    //  fMixedBalance->CalculateBalance(fReactionPlane,tracksMixing,NULL,1,1.,0.);	     

  }//event loop
}      

//________________________________________________________________________
void  AliAnalysisTaskToyModel::FinishOutput() {
  //Printf("END BF");

  TFile *gOutput = TFile::Open("outputToyModel.root","recreate");
  TDirectoryFile *dir = new TDirectoryFile("PWGCFEbyE.outputBalanceFunctionPsiAnalysis","PWGCFEbyE.outputBalanceFunctionPsiAnalysis");
  
  fList->SetName("listQA");
  fList->SetOwner(kTRUE);
  dir->Add(fList);
  //fList->Write("listQA",TObject::kSingleKey);

  if (!fBalance) {
    Printf("ERROR: fBalance not available");
    return;
  }  
  fListBF->SetName("listBF");
  fListBF->SetOwner(kTRUE);
  dir->Add(fListBF);
  //fListBF->Write("listBF",TObject::kSingleKey);

  if(fRunShuffling) {
    if (!fShuffledBalance) {
      Printf("ERROR: fShuffledBalance not available");
      return;
    }
    fListBFS->SetName("listBFShuffled");
    fListBFS->SetOwner(kTRUE);
    dir->Add(fListBFS);
    //fListBFS->Write("listBFShuffled",TObject::kSingleKey);
  }

  if(fRunMixing) {
    if (!fMixedBalance) {
      Printf("ERROR: fMixedBalance not available");
      return;
    }
    fListBFM->SetName("listBFMixed");
    fListBFM->SetOwner(kTRUE);
    dir->Add(fListBFM);
    //fListBFM->Write("listBFMixed",TObject::kSingleKey);
  }

  dir->Write(dir->GetName(),TObject::kSingleKey);
  gOutput->Close();
}

