#include "TChain.h"
#include "TList.h"
#include "TCanvas.h"
#include "TParticle.h"
#include "TLorentzVector.h"
#include "TGraphErrors.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TArrayF.h"
#include "TF1.h"
#include "TRandom.h"
#include "TFile.h"

#include "AliAnalysisManager.h"
#include "AliLog.h"

#include "AliAnalysisTaskToyModel.h"
#include "AliBalance.h"


// Analysis task for the toy model analysis
// Authors: Panos.Christakoglou@nikhef.nl

ClassImp(AliAnalysisTaskToyModel)

//________________________________________________________________________
AliAnalysisTaskToyModel::AliAnalysisTaskToyModel() 
: TObject(),
  fUseDebug(kFALSE),
  fBalance(0),
  fRunShuffling(kFALSE), fShuffledBalance(0),
  fList(0), fListBF(0), fListBFS(0),
  fHistEventStats(0),
  fHistNumberOfAcceptedParticles(0),
  fHistReactionPlane(0),
  fHistEtaTotal(0), fHistEta(0),
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
  fUseDynamicalCorrelations(kFALSE), fDynamicalCorrelationsPercentage(0.1) {
  // Constructor
}

//________________________________________________________________________
AliAnalysisTaskToyModel::~AliAnalysisTaskToyModel() {
  //Destructor
  delete fPtSpectraAllCharges;
  delete fAzimuthalAngleAllCharges;
  delete fPtSpectraPions;
  delete fAzimuthalAnglePions;
  delete fPtSpectraKaons;
  delete fAzimuthalAngleKaons;
  delete fPtSpectraProtons;
  delete fAzimuthalAngleProtons;
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
}

//________________________________________________________________________
void AliAnalysisTaskToyModel::CreateOutputObjects() {
  // Create histograms
  // Called once
  if(!fBalance) {
    fBalance = new AliBalance();
    fBalance->SetInterval(-1,-0.8,0.8,16,0.,1.6);
  }
  if(fRunShuffling) {
    if(!fShuffledBalance) {
      fShuffledBalance = new AliBalance();
      fShuffledBalance->SetInterval(-1,-0.8,0.8,16,0.,1.6);
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

  //==============Balance function histograms================//
  // Initialize histograms if not done yet
  if(!fBalance->GetHistNp(0)){
    AliWarning("Histograms not yet initialized! --> Will be done now");
    AliWarning("--> Add 'gBalance->InitHistograms()' in your configBalanceFunction");
    fBalance->InitHistograms();
  }

  if(fRunShuffling) {
    if(!fShuffledBalance->GetHistNp(0)) {
      AliWarning("Histograms (shuffling) not yet initialized! --> Will be done now");
      AliWarning("--> Add 'gBalance->InitHistograms()' in your configBalanceFunction");
      fShuffledBalance->InitHistograms();
    }
  }

  for(Int_t a = 0; a < ANALYSIS_TYPES; a++){
    fListBF->Add(fBalance->GetHistNp(a));
    fListBF->Add(fBalance->GetHistNn(a));
    fListBF->Add(fBalance->GetHistNpn(a));
    fListBF->Add(fBalance->GetHistNnn(a));
    fListBF->Add(fBalance->GetHistNpp(a));
    fListBF->Add(fBalance->GetHistNnp(a));

    if(fRunShuffling) {
      fListBFS->Add(fShuffledBalance->GetHistNp(a));
      fListBFS->Add(fShuffledBalance->GetHistNn(a));
      fListBFS->Add(fShuffledBalance->GetHistNpn(a));
      fListBFS->Add(fShuffledBalance->GetHistNnn(a));
      fListBFS->Add(fShuffledBalance->GetHistNpp(a));
      fListBFS->Add(fShuffledBalance->GetHistNnp(a));
    }  
  }

  // Post output data.
  //PostData(1, fList);
  //PostData(2, fListBF);
  //if(fRunShuffling) PostData(3, fListBFS);
}

//________________________________________________________________________
void AliAnalysisTaskToyModel::Run(Int_t nEvents) {
  // Main loop
  // Called for each event
  Double_t v_charge = 0;
  Double_t v_y = 0.0;
  Double_t v_eta = 0.0;
  Double_t v_phi = 0.0;
  Double_t v_p[3] = {0.,0.,0.};
  Double_t v_pt = 0.0;
  Double_t v_E = 0.0;
  Bool_t isPion = kFALSE, isKaon = kFALSE, isProton = kFALSE;

  if(fUseAllCharges) {
    fPtSpectraAllCharges->SetParameter(0,fParticleMass);
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
    vector<Double_t> *chargeVectorShuffle[9];   // this will be shuffled
    vector<Double_t> *chargeVector[9];          // original charge
    for(Int_t i = 0; i < 9; i++){
      chargeVectorShuffle[i] = new vector<Double_t>;
      chargeVector[i]        = new vector<Double_t>;
    }

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

    //Int_t nGeneratedPositive = 0, nGeneratedNegative = 0;
    Int_t nGeneratedPions = 0, nGeneratedKaons = 0, nGeneratedProtons = 0;
    
    //Randomization of the reaction plane
    //fReactionPlane = 2.0*TMath::Pi()*gRandom->Rndm();
    fReactionPlane = 0.0;
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
 
    //Generate positive particles
    for(Int_t iParticleCount = 0; iParticleCount < nGeneratedPositive; iParticleCount++) {
      isPion = kFALSE; isKaon = kFALSE; isProton = kFALSE;
      if(fUseDebug) 
	Printf("Generating positive: %d(%d)",iParticleCount+1,nGeneratedPositive);

      //Pseudo-rapidity sampled from a Gaussian centered @ 0
      v_eta = gRandom->Gaus(0.0,4.0);

      //Fill QA histograms (full phase space)
      fHistEtaTotal->Fill(v_eta);

      v_charge = 1.0;
      //nGeneratedPositive += 1;
      
      //Acceptance
      if((v_eta < fEtaMin) || (v_eta > fEtaMax)) continue;

      if(!fUseAllCharges) {
	//Decide the specie
	Double_t randomNumberSpecies = gRandom->Rndm();
	if((randomNumberSpecies >= 0.0)&&(randomNumberSpecies < fPionPercentage)) {
	  nGeneratedPions += 1;
	  v_pt = fPtSpectraPions->GetRandom();
	  v_phi = fAzimuthalAnglePions->GetRandom();
	  fParticleMass = fPionMass;
	  isPion = kTRUE;
	}
	else if((randomNumberSpecies >= fPionPercentage)&&(randomNumberSpecies < fPionPercentage + fKaonPercentage)) {
	  nGeneratedKaons += 1;
	  v_pt = fPtSpectraKaons->GetRandom();
	  v_phi = fAzimuthalAngleKaons->GetRandom();
	  fParticleMass = fKaonMass;
	  isKaon = kTRUE;
	}
	else if((randomNumberSpecies >= fPionPercentage + fKaonPercentage)&&(randomNumberSpecies < fPionPercentage + fKaonPercentage + fProtonPercentage)) {
	  nGeneratedProtons += 1;
	  v_pt = fPtSpectraProtons->GetRandom();
	  v_phi = fAzimuthalAngleProtons->GetRandom();
	  fParticleMass = fProtonMass;
	  isProton = kTRUE;
	}
      }
      else {
	v_pt = fPtSpectraAllCharges->GetRandom();
	v_phi = fAzimuthalAngleAllCharges->GetRandom();
      }
      
      v_p[0] = v_pt*TMath::Cos(v_phi);
      v_p[1] = v_pt*TMath::Sin(v_phi);
      v_p[2] = v_pt*TMath::SinH(v_eta);
      v_E = TMath::Sqrt(TMath::Power(fParticleMass,2) +
			TMath::Power(v_p[0],2) +
			TMath::Power(v_p[1],2) +
			TMath::Power(v_p[2],2));
      
      v_y = 0.5*TMath::Log((v_E + v_p[2])/(v_E - v_p[2]));
      
      //pt coverage
      if((v_pt < fPtMin) || (v_pt > fPtMax)) continue;
      //Printf("pt: %lf - mins: %lf - max: %lf",v_pt,fPtMin,fPtMax);

      //acceptance filter
      if(fUseAcceptanceParameterization) {
	Double_t gRandomNumberForAcceptance = gRandom->Rndm();
	if(gRandomNumberForAcceptance > fAcceptanceParameterization->Eval(v_pt)) 
	  continue;
      }

      gNumberOfAcceptedPositiveParticles += 1;

      //Fill QA histograms (acceptance);
      fHistEta->Fill(v_eta);
      fHistRapidity->Fill(v_y);
      fHistPhi->Fill(v_phi);
      fHistPt->Fill(v_pt);
      if(isPion) {
	fHistRapidityPions->Fill(v_y);
	fHistPhiPions->Fill(v_phi);
	fHistPtPions->Fill(v_pt);
      }
      else if(isKaon) {
	fHistRapidityKaons->Fill(v_y);
	fHistPhiKaons->Fill(v_phi);
	fHistPtKaons->Fill(v_pt);
      }
      else if(isProton) {
	fHistRapidityProtons->Fill(v_y);
	fHistPhiProtons->Fill(v_phi);
	fHistPtProtons->Fill(v_pt);
      }

      // fill charge vector
      chargeVector[0]->push_back(v_charge);
      chargeVector[1]->push_back(v_y);
      chargeVector[2]->push_back(v_eta);
      chargeVector[3]->push_back(TMath::RadToDeg()*v_phi);
      chargeVector[4]->push_back(v_p[0]);
      chargeVector[5]->push_back(v_p[1]);
      chargeVector[6]->push_back(v_p[2]);
      chargeVector[7]->push_back(v_pt);
      chargeVector[8]->push_back(v_E);
      
      if(fRunShuffling) {
	chargeVectorShuffle[0]->push_back(v_charge);
	chargeVectorShuffle[1]->push_back(v_y);
	chargeVectorShuffle[2]->push_back(v_eta);
	chargeVectorShuffle[3]->push_back(TMath::RadToDeg()*v_phi);
	chargeVectorShuffle[4]->push_back(v_p[0]);
	chargeVectorShuffle[5]->push_back(v_p[1]);
	chargeVectorShuffle[6]->push_back(v_p[2]);
	chargeVectorShuffle[7]->push_back(v_pt);
	chargeVectorShuffle[8]->push_back(v_E);
      }
      gNumberOfAcceptedParticles += 1;
    }//generated positive particle loop
 
    //Generate negative particles
    for(Int_t iParticleCount = 0; iParticleCount < nGeneratedNegative; iParticleCount++) {
      isPion = kFALSE; isKaon = kFALSE; isProton = kFALSE;
      if(fUseDebug) 
	Printf("Generating negative: %d(%d)",iParticleCount+1,nGeneratedNegative);

      //Pseudo-rapidity sampled from a Gaussian centered @ 0
      v_eta = gRandom->Gaus(0.0,4.0);

      //Fill QA histograms (full phase space)
      fHistEtaTotal->Fill(v_eta);

      v_charge = -1.0;
      //nGeneratedNegative += 1;
      
      //Acceptance
      if((v_eta < fEtaMin) || (v_eta > fEtaMax)) continue;

      if(!fUseAllCharges) {
	//Decide the specie
	Double_t randomNumberSpecies = gRandom->Rndm();
	if((randomNumberSpecies >= 0.0)&&(randomNumberSpecies < fPionPercentage)) {
	  nGeneratedPions += 1;
	  v_pt = fPtSpectraPions->GetRandom();
	  v_phi = fAzimuthalAnglePions->GetRandom();
	  fParticleMass = fPionMass;
	  isPion = kTRUE;
	}
	else if((randomNumberSpecies >= fPionPercentage)&&(randomNumberSpecies < fPionPercentage + fKaonPercentage)) {
	  nGeneratedKaons += 1;
	  v_pt = fPtSpectraKaons->GetRandom();
	  v_phi = fAzimuthalAngleKaons->GetRandom();
	  fParticleMass = fKaonMass;
	  isKaon = kTRUE;
	}
	else if((randomNumberSpecies >= fPionPercentage + fKaonPercentage)&&(randomNumberSpecies < fPionPercentage + fKaonPercentage + fProtonPercentage)) {
	  nGeneratedProtons += 1;
	  v_pt = fPtSpectraProtons->GetRandom();
	  v_phi = fAzimuthalAngleProtons->GetRandom();
	  fParticleMass = fProtonMass;
	  isProton = kTRUE;
	}
      }
      else {
	v_pt = fPtSpectraAllCharges->GetRandom();
	v_phi = fAzimuthalAngleAllCharges->GetRandom();
      }
      
      v_p[0] = v_pt*TMath::Cos(v_phi);
      v_p[1] = v_pt*TMath::Sin(v_phi);
      v_p[2] = v_pt*TMath::SinH(v_eta);
      v_E = TMath::Sqrt(TMath::Power(fParticleMass,2) +
			TMath::Power(v_p[0],2) +
			TMath::Power(v_p[1],2) +
			TMath::Power(v_p[2],2));
      
      v_y = 0.5*TMath::Log((v_E + v_p[2])/(v_E - v_p[2]));
      
      //pt coverage
      if((v_pt < fPtMin) || (v_pt > fPtMax)) continue;
      //Printf("pt: %lf - mins: %lf - max: %lf",v_pt,fPtMin,fPtMax);

     //acceptance filter
      if(fUseAcceptanceParameterization) {
	Double_t gRandomNumberForAcceptance = gRandom->Rndm();
	if(gRandomNumberForAcceptance > fAcceptanceParameterization->Eval(v_pt)) 
	  continue;
      }

      gNumberOfAcceptedNegativeParticles += 1;

      //Fill QA histograms (acceptance);
      fHistEta->Fill(v_eta);
      fHistRapidity->Fill(v_y);
      fHistPhi->Fill(v_phi);
      fHistPt->Fill(v_pt);
      if(isPion) {
	fHistRapidityPions->Fill(v_y);
	fHistPhiPions->Fill(v_phi);
	fHistPtPions->Fill(v_pt);
      }
      else if(isKaon) {
	fHistRapidityKaons->Fill(v_y);
	fHistPhiKaons->Fill(v_phi);
	fHistPtKaons->Fill(v_pt);
      }
      else if(isProton) {
	fHistRapidityProtons->Fill(v_y);
	fHistPhiProtons->Fill(v_phi);
	fHistPtProtons->Fill(v_pt);
      }

      // fill charge vector
      chargeVector[0]->push_back(v_charge);
      chargeVector[1]->push_back(v_y);
      chargeVector[2]->push_back(v_eta);
      chargeVector[3]->push_back(TMath::RadToDeg()*v_phi);
      chargeVector[4]->push_back(v_p[0]);
      chargeVector[5]->push_back(v_p[1]);
      chargeVector[6]->push_back(v_p[2]);
      chargeVector[7]->push_back(v_pt);
      chargeVector[8]->push_back(v_E);
      
      if(fRunShuffling) {
	chargeVectorShuffle[0]->push_back(v_charge);
	chargeVectorShuffle[1]->push_back(v_y);
	chargeVectorShuffle[2]->push_back(v_eta);
	chargeVectorShuffle[3]->push_back(TMath::RadToDeg()*v_phi);
	chargeVectorShuffle[4]->push_back(v_p[0]);
	chargeVectorShuffle[5]->push_back(v_p[1]);
	chargeVectorShuffle[6]->push_back(v_p[2]);
	chargeVectorShuffle[7]->push_back(v_pt);
	chargeVectorShuffle[8]->push_back(v_E);
      }
      gNumberOfAcceptedParticles += 1;
    }//generated negative particle loop
    
    //Dynamical correlations
    Double_t v_chargePrime = 0;
    Double_t v_yPrime = 0.0;
    Double_t v_etaPrime = 0.0;
    Double_t v_phiPrime = 0.0;
    Double_t v_pPrime[3] = {0.,0.,0.};
    Double_t v_ptPrime = 0.0;
    Double_t v_EPrime = 0.0;
    Int_t nGeneratedPositiveDynamicalCorrelations = 0;
    Int_t nGeneratedNegativeDynamicalCorrelations = 0;
    //Generate "correlated" particles 
    if(fUseDynamicalCorrelations) {
      Int_t gNumberOfDynamicalCorrelations = (Int_t)(0.5*gNumberOfAcceptedParticles*fDynamicalCorrelationsPercentage);
      for(Int_t iDynamicalCorrelations = 0; iDynamicalCorrelations < gNumberOfDynamicalCorrelations; iDynamicalCorrelations++) {
	isPion = kFALSE; isKaon = kFALSE; isProton = kFALSE;
	
	//Pseudo-rapidity sampled from a Gaussian centered @ 0
	v_eta = gRandom->Gaus(0.0,0.1);
	v_charge = 1.0;
	nGeneratedPositiveDynamicalCorrelations += 1;
	
	v_etaPrime = -v_eta;
	v_chargePrime = -1.0;
	nGeneratedNegativeDynamicalCorrelations += 1;
	  
	//Acceptance
	if((v_eta < fEtaMin) || (v_eta > fEtaMax)) continue;
	if((v_etaPrime < fEtaMin) || (v_etaPrime > fEtaMax)) continue;
      
	if(!fUseAllCharges) {
	  //Decide the specie
	  Double_t randomNumberSpecies = gRandom->Rndm();
	  if((randomNumberSpecies >= 0.0)&&(randomNumberSpecies < fPionPercentage)) {
	    nGeneratedPions += 1;
	    v_pt = fPtSpectraPions->GetRandom();
	    v_phi = fAzimuthalAnglePions->GetRandom();
	    fParticleMass = fPionMass;
	    isPion = kTRUE;
	  }
	  else if((randomNumberSpecies >= fPionPercentage)&&(randomNumberSpecies < fPionPercentage + fKaonPercentage)) {
	    nGeneratedKaons += 1;
	    v_pt = fPtSpectraKaons->GetRandom();
	    v_phi = fAzimuthalAngleKaons->GetRandom();
	    fParticleMass = fKaonMass;
	    isKaon = kTRUE;
	  }
	  else if((randomNumberSpecies >= fPionPercentage + fKaonPercentage)&&(randomNumberSpecies < fPionPercentage + fKaonPercentage + fProtonPercentage)) {
	    nGeneratedProtons += 1;
	    v_pt = fPtSpectraProtons->GetRandom();
	    v_ptPrime = v_pt;
	    v_phi = fAzimuthalAngleProtons->GetRandom();
	    fParticleMass = fProtonMass;
	    isProton = kTRUE;
	  }
	}
	else {
	  v_pt = fPtSpectraAllCharges->GetRandom();
	  v_phi = fAzimuthalAngleAllCharges->GetRandom();
	}
	v_ptPrime = v_pt;
	v_phiPrime = v_phi;

	v_p[0] = v_pt*TMath::Cos(v_phi);
	v_p[1] = v_pt*TMath::Sin(v_phi);
	v_p[2] = v_pt*TMath::SinH(v_eta);
	v_E = TMath::Sqrt(TMath::Power(fParticleMass,2) +
			  TMath::Power(v_p[0],2) +
			  TMath::Power(v_p[1],2) +
			  TMath::Power(v_p[2],2));
	
	v_y = 0.5*TMath::Log((v_E + v_p[2])/(v_E - v_p[2]));

	v_pPrime[0] = v_ptPrime*TMath::Cos(v_phiPrime);
	v_pPrime[1] = v_ptPrime*TMath::Sin(v_phiPrime);
	v_pPrime[2] = v_ptPrime*TMath::SinH(v_etaPrime);
	v_EPrime = TMath::Sqrt(TMath::Power(fParticleMass,2) +
			  TMath::Power(v_pPrime[0],2) +
			  TMath::Power(v_pPrime[1],2) +
			  TMath::Power(v_pPrime[2],2));
	
	v_yPrime = 0.5*TMath::Log((v_EPrime + v_pPrime[2])/(v_EPrime - v_pPrime[2]));
      
	//pt coverage
	if((v_pt < fPtMin) || (v_pt > fPtMax)) continue;
	if((v_ptPrime < fPtMin) || (v_ptPrime > fPtMax)) continue;

	//acceptance filter
	if(fUseAcceptanceParameterization) {
	  Double_t gRandomNumberForAcceptance = gRandom->Rndm();
	  if(gRandomNumberForAcceptance > fAcceptanceParameterization->Eval(v_pt)) 
	    continue;
	  
	  Double_t gRandomNumberForAcceptancePrime = gRandom->Rndm();
	  if(gRandomNumberForAcceptancePrime > fAcceptanceParameterization->Eval(v_ptPrime)) 
	    continue;
	}
      
	// fill charge vector (positive)
	chargeVector[0]->push_back(v_charge);
	chargeVector[1]->push_back(v_y);
	chargeVector[2]->push_back(v_eta);
	chargeVector[3]->push_back(TMath::RadToDeg()*v_phi);
	chargeVector[4]->push_back(v_p[0]);
	chargeVector[5]->push_back(v_p[1]);
	chargeVector[6]->push_back(v_p[2]);
	chargeVector[7]->push_back(v_pt);
	chargeVector[8]->push_back(v_E);
	
	if(fRunShuffling) {
	  chargeVectorShuffle[0]->push_back(v_charge);
	  chargeVectorShuffle[1]->push_back(v_y);
	  chargeVectorShuffle[2]->push_back(v_eta);
	  chargeVectorShuffle[3]->push_back(TMath::RadToDeg()*v_phi);
	  chargeVectorShuffle[4]->push_back(v_p[0]);
	  chargeVectorShuffle[5]->push_back(v_p[1]);
	  chargeVectorShuffle[6]->push_back(v_p[2]);
	  chargeVectorShuffle[7]->push_back(v_pt);
	  chargeVectorShuffle[8]->push_back(v_E);
	}

	// fill charge vector (negative)
	chargeVector[0]->push_back(v_chargePrime);
	chargeVector[1]->push_back(v_yPrime);
	chargeVector[2]->push_back(v_etaPrime);
	chargeVector[3]->push_back(TMath::RadToDeg()*v_phiPrime);
	chargeVector[4]->push_back(v_pPrime[0]);
	chargeVector[5]->push_back(v_pPrime[1]);
	chargeVector[6]->push_back(v_pPrime[2]);
	chargeVector[7]->push_back(v_ptPrime);
	chargeVector[8]->push_back(v_EPrime);
	
	if(fRunShuffling) {
	  chargeVectorShuffle[0]->push_back(v_chargePrime);
	  chargeVectorShuffle[1]->push_back(v_yPrime);
	  chargeVectorShuffle[2]->push_back(v_etaPrime);
	  chargeVectorShuffle[3]->push_back(TMath::RadToDeg()*v_phiPrime);
	  chargeVectorShuffle[4]->push_back(v_pPrime[0]);
	  chargeVectorShuffle[5]->push_back(v_pPrime[1]);
	  chargeVectorShuffle[6]->push_back(v_pPrime[2]);
	  chargeVectorShuffle[7]->push_back(v_ptPrime);
	  chargeVectorShuffle[8]->push_back(v_EPrime);
	}

	gNumberOfAcceptedParticles += 2;
      }//loop over the dynamical correlations
    }// usage of dynamical correlations

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

    //Calculate the balance function
    fBalance->CalculateBalance(gNumberOfAcceptedParticles,chargeVector);
    
    if(fRunShuffling) {
      // shuffle charges
      random_shuffle( chargeVectorShuffle[0]->begin(), chargeVectorShuffle[0]->end() );
      fShuffledBalance->CalculateBalance(gNumberOfAcceptedParticles,chargeVectorShuffle);
    }
  }
}      

//________________________________________________________________________
void  AliAnalysisTaskToyModel::FinishOutput() {
  //Printf("END BF");

  TFile *gOutput = TFile::Open("outputToyModel.root","recreate");
  fList->Write();

  if (!fBalance) {
    Printf("ERROR: fBalance not available");
    return;
  }  
  fListBF->Write();

  if(fRunShuffling) {
    if (!fShuffledBalance) {
      Printf("ERROR: fShuffledBalance not available");
      return;
    }
    fListBFS->Write();
  }
  gOutput->Close();
}

