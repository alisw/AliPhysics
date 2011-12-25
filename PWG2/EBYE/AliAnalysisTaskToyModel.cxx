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

#include "AliAnalysisTaskSE.h"
#include "AliAnalysisManager.h"
#include "AliLog.h"

#include "AliAnalysisTaskToyModel.h"
#include "AliBalance.h"


// Analysis task for the toy model analysis
// Authors: Panos.Christakoglou@nikhef.nl

ClassImp(AliAnalysisTaskToyModel)

//________________________________________________________________________
AliAnalysisTaskToyModel::AliAnalysisTaskToyModel(const char *name) 
: AliAnalysisTaskSE(name), 
  fBalance(0),
  fRunShuffling(kFALSE), fShuffledBalance(0),
  fList(0), fListBF(0), fListBFS(0),
  fHistEventStats(0),
  fTotalMultiplicityMean(100.), fTotalMultiplicitySigma(0.0),
  fNetChargeMean(0.0), fNetChargeSigma(0.0),
  fPtMin(0.0), fPtMax(100.0),
  fEtaMin(-1.0), fEtaMax(1.0),
  fUseAcceptanceParameterization(kFALSE), fAcceptanceParameterization(0),
  fPtSpectraAllCharges(0), fTemperatureAllCharges(100.),
  fReactionPlane(0.0),
  fAzimuthalAngleAllCharges(0), fDirectedFlowAllCharges(0.0), 
  fEllipticFlowAllCharges(0.0), fTriangularFlowAllCharges(0.0),
  fQuandrangularFlowAllCharges(0.0), fPentangularFlowAllCharges(0.0),
  fPionPercentage(0.8),
  fPtSpectraPions(0), fTemperaturePions(100.),
  fAzimuthalAnglePions(0), fDirectedFlowPions(0.0), 
  fEllipticFlowPions(0.0), fTriangularFlowPions(0.0), 
  fQuandrangularFlowPions(0.0), fPentangularFlowPions(0.0),
  fKaonPercentage(0.8),
  fPtSpectraKaons(0), fTemperatureKaons(100.),
  fAzimuthalAngleKaons(0), fDirectedFlowKaons(0.0), 
  fEllipticFlowKaons(0.0), fTriangularFlowKaons(0.0),
  fQuandrangularFlowKaons(0.0), fPentangularFlowKaons(0.0),
  fProtonPercentage(0.8),
  fPtSpectraProtons(0), fTemperatureProtons(100.),
  fAzimuthalAngleProtons(0), fDirectedFlowProtons(0.0), 
  fEllipticFlowProtons(0.0), fTriangularFlowProtons(0.0),
  fQuandrangularFlowProtons(0.0), fPentangularFlowProtons(0.0) {
  // Constructor

  // Define input and output slots here
  // Input slot #0 works with a TChain
  DefineInput(0, TChain::Class());
  // Output slot #0 writes into a TH1 container
  DefineOutput(1, TList::Class());
  DefineOutput(2, TList::Class());
  DefineOutput(3, TList::Class());
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
  //==============Particles and spectra==============//
  TParticle *pion = new TParticle();
  pion->SetPdgCode(211);
  Double_t gPionMass = pion->GetMass();

  fPtSpectraAllCharges = new TF1("fPtSpectraAllCharges","x*TMath::Exp(-TMath::Power([0]*[0]+x*x,0.5)/[1])",0.,100.);
  fPtSpectraAllCharges->SetParName(0,"Mass");
  fPtSpectraAllCharges->SetParameter(0,gPionMass);
  fPtSpectraAllCharges->SetParName(1,"Temperature");
  fPtSpectraAllCharges->SetParameter(1,fTemperatureAllCharges);

  fPtSpectraPions = new TF1("fPtSpectraPions","x*TMath::Exp(-TMath::Power([0]*[0]+x*x,0.5)/[1])",0.,100.);
  fPtSpectraPions->SetParName(0,"Mass");
  fPtSpectraPions->SetParameter(0,gPionMass);
  fPtSpectraPions->SetParName(1,"Temperature");
  fPtSpectraPions->SetParameter(1,fTemperaturePions);

  TParticle *kaon = new TParticle();
  kaon->SetPdgCode(321);
  Double_t gKaonMass = kaon->GetMass();
  fPtSpectraKaons = new TF1("fPtSpectraKaons","x*TMath::Exp(-TMath::Power([0]*[0]+x*x,0.5)/[1])",0.,100.);
  fPtSpectraKaons->SetParName(0,"Mass");
  fPtSpectraKaons->SetParameter(0,gKaonMass);
  fPtSpectraKaons->SetParName(1,"Temperature");
  fPtSpectraKaons->SetParameter(1,fTemperatureKaons);

  TParticle *proton = new TParticle();
  proton->SetPdgCode(2212);
  Double_t gProtonMass = proton->GetMass();
  fPtSpectraProtons = new TF1("fPtSpectraProtons","x*TMath::Exp(-TMath::Power([0]*[0]+x*x,0.5)/[1])",0.,5.);
  fPtSpectraProtons->SetParName(0,"Mass");
  fPtSpectraProtons->SetParameter(0,gProtonMass);
  fPtSpectraProtons->SetParName(1,"Temperature");
  fPtSpectraProtons->SetParameter(1,fTemperatureProtons);
  //==============Particles and spectra==============//

  //==============Flow values==============//
  fAzimuthalAngleAllCharges = new TF1("fAzimuthalAngleAllCharges","1+2.*[1]*TMath::Cos(x-[0])+2.*[2]*TMath::Cos(2*(x-[0]))+2.*[3]*TMath::Cos(3*(x-[0]))+2.*[4]*TMath::Cos(4*(x-[0]))+2.*[5]*TMath::Cos(5*(x-[0]))",0.,2.*TMath::Pi());
  fAzimuthalAngleAllCharges->SetParName(0,"Reaction Plane");
  fAzimuthalAngleAllCharges->SetParameter(0,fReactionPlane);
  fAzimuthalAngleAllCharges->SetParName(1,"Directed flow");
  fAzimuthalAngleAllCharges->SetParameter(1,fDirectedFlowAllCharges);
  fAzimuthalAngleAllCharges->SetParName(2,"Elliptic flow"); 
  fAzimuthalAngleAllCharges->SetParameter(2,fEllipticFlowAllCharges);
  fAzimuthalAngleAllCharges->SetParName(3,"Triangular flow");
  fAzimuthalAngleAllCharges->SetParameter(3,fTriangularFlowAllCharges);
  fAzimuthalAngleAllCharges->SetParName(4,"Quandrangular flow");
  fAzimuthalAngleAllCharges->SetParameter(4,fQuandrangularFlowAllCharges);
  fAzimuthalAngleAllCharges->SetParName(5,"Pentangular flow");
  fAzimuthalAngleAllCharges->SetParameter(5,fPentangularFlowAllCharges);

  fAzimuthalAnglePions = new TF1("fAzimuthalAnglePions","1+2.*[1]*TMath::Cos(x-[0])+2.*[2]*TMath::Cos(2*(x-[0]))+2.*[3]*TMath::Cos(3*(x-[0]))+2.*[4]*TMath::Cos(4*(x-[0]))+2.*[5]*TMath::Cos(5*(x-[0]))",0.,2.*TMath::Pi());
  fAzimuthalAnglePions->SetParName(0,"Reaction Plane");
  fAzimuthalAnglePions->SetParameter(0,fReactionPlane);
  fAzimuthalAnglePions->SetParName(1,"Directed flow");
  fAzimuthalAnglePions->SetParameter(1,fDirectedFlowPions);
  fAzimuthalAnglePions->SetParName(2,"Elliptic flow"); 
  fAzimuthalAnglePions->SetParameter(2,fEllipticFlowPions);
  fAzimuthalAnglePions->SetParName(3,"Triangular flow");
  fAzimuthalAnglePions->SetParameter(3,fTriangularFlowPions);
  fAzimuthalAnglePions->SetParName(4,"Quandrangular flow");
  fAzimuthalAnglePions->SetParameter(4,fQuandrangularFlowPions);
  fAzimuthalAnglePions->SetParName(5,"Pentangular flow");
  fAzimuthalAnglePions->SetParameter(5,fPentangularFlowPions);

  fAzimuthalAngleKaons = new TF1("fAzimuthalAngleKaons","1+2.*[1]*TMath::Cos(x-[0])+2.*[2]*TMath::Cos(2*(x-[0]))+2.*[3]*TMath::Cos(3*(x-[0]))+2.*[4]*TMath::Cos(4*(x-[0]))+2.*[5]*TMath::Cos(5*(x-[0]))",0.,2.*TMath::Pi());
  fAzimuthalAngleKaons->SetParName(0,"Reaction Plane");
  fAzimuthalAngleKaons->SetParameter(0,fReactionPlane);
  fAzimuthalAngleKaons->SetParName(1,"Directed flow");
  fAzimuthalAngleKaons->SetParameter(1,fDirectedFlowKaons);
  fAzimuthalAngleKaons->SetParName(2,"Elliptic flow"); 
  fAzimuthalAngleKaons->SetParameter(2,fEllipticFlowKaons);
  fAzimuthalAngleKaons->SetParName(3,"Triangular flow");
  fAzimuthalAngleKaons->SetParameter(3,fTriangularFlowKaons);
  fAzimuthalAngleKaons->SetParName(4,"Quandrangular flow");
  fAzimuthalAngleKaons->SetParameter(4,fQuandrangularFlowKaons);
  fAzimuthalAngleKaons->SetParName(5,"Pentangular flow");
  fAzimuthalAngleKaons->SetParameter(5,fPentangularFlowKaons);

  fAzimuthalAngleProtons = new TF1("fAzimuthalAngleProtons","1+2.*[1]*TMath::Cos(x-[0])+2.*[2]*TMath::Cos(2*(x-[0]))+2.*[3]*TMath::Cos(3*(x-[0]))+2.*[4]*TMath::Cos(4*(x-[0]))+2.*[5]*TMath::Cos(5*(x-[0]))",0.,2.*TMath::Pi());
  fAzimuthalAngleProtons->SetParName(0,"Reaction Plane");
  fAzimuthalAngleProtons->SetParameter(0,fReactionPlane);
  fAzimuthalAngleProtons->SetParName(1,"Directed flow");
  fAzimuthalAngleProtons->SetParameter(1,fDirectedFlowProtons);
  fAzimuthalAngleProtons->SetParName(2,"Elliptic flow"); 
  fAzimuthalAngleProtons->SetParameter(2,fEllipticFlowProtons);
  fAzimuthalAngleProtons->SetParName(3,"Triangular flow");
  fAzimuthalAngleProtons->SetParameter(3,fTriangularFlowProtons);
  fAzimuthalAngleProtons->SetParName(4,"Quandrangular flow");
  fAzimuthalAngleProtons->SetParameter(4,fQuandrangularFlowProtons);
  fAzimuthalAngleProtons->SetParName(5,"Pentangular flow");
  fAzimuthalAngleProtons->SetParameter(5,fPentangularFlowProtons);
  //==============Flow values==============//
}

//________________________________________________________________________
void AliAnalysisTaskToyModel::UserCreateOutputObjects() {
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

  //Event stats.
  TString gCutName[4] = {"Total","Offline trigger",
                         "Vertex","Analyzed"};
  fHistEventStats = new TH1F("fHistEventStats",
                             "Event statistics;;N_{events}",
                             4,0.5,4.5);
  for(Int_t i = 1; i <= 4; i++)
    fHistEventStats->GetXaxis()->SetBinLabel(i,gCutName[i-1].Data());
  fList->Add(fHistEventStats);

  // Balance function histograms
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
  PostData(1, fList);
  PostData(2, fListBF);
  if(fRunShuffling) PostData(3, fListBFS);
}

//________________________________________________________________________
void AliAnalysisTaskToyModel::UserExec(Option_t *) {
  // Main loop
  // Called for each event

  // vector holding the charges/kinematics of all tracks (charge,y,eta,phi,p0,p1,p2,pt,E)
  vector<Double_t> *chargeVectorShuffle[9];   // this will be shuffled
  vector<Double_t> *chargeVector[9];          // original charge
  for(Int_t i = 0; i < 9; i++){
    chargeVectorShuffle[i] = new vector<Double_t>;
    chargeVector[i]        = new vector<Double_t>;
  }

  Double_t v_charge;
  Double_t v_y;
  Double_t v_eta;
  Double_t v_phi;
  Double_t v_p[3];
  Double_t v_pt;
  Double_t v_E;

  //Multiplicities
  Int_t nMultiplicity = (Int_t)(gRandom->Gaus(fTotalMultiplicityMean,fTotalMultiplicitySigma));
  Int_t nNetCharge = (Int_t)(gRandom->Gaus(fNetChargeMean,fNetChargeSigma));
  Int_t nGeneratedPositive = 0.5*(nMultiplicity + nNetCharge);
  Int_t nGeneratedNegative = nMultiplicity - nGeneratedPositive;
  Int_t nGeneratedPositivePions = (Int_t)(fPionPercentage*nGeneratedPositive);
  Int_t nGeneratedNegativePions = (Int_t)(fPionPercentage*nGeneratedNegative);
  Int_t nGeneratedPositiveKaons = (Int_t)(fKaonPercentage*nGeneratedPositive);
  Int_t nGeneratedNegativeKaons = (Int_t)(fKaonPercentage*nGeneratedNegative);
  Int_t nGeneratedPositiveProtons = (Int_t)(fProtonPercentage*nGeneratedPositive);
  Int_t nGeneratedNegativeProtons = (Int_t)(fProtonPercentage*nGeneratedNegative);  

  Printf("Total positive: %d - Total negative: %d",nGeneratedPositive,nGeneratedNegative);
  Printf("Positive pions: %d - Negative pions: %d",nGeneratedPositivePions,nGeneratedNegativePions);
  Printf("Positive kaons: %d - Negative kaons: %d",nGeneratedPositiveKaons,nGeneratedNegativeKaons);
  Printf("Positive protons: %d - Negative protons: %d",nGeneratedPositiveProtons,nGeneratedNegativeProtons);

  Int_t gNumberOfAcceptedParticles = 0;
  //positive particles
  for(Int_t iPosCount = 0; iPosCount < nGeneratedPositive; iPosCount++) {
    v_charge = 1.0;

    v_pt = fPtSpectraAllCharges->GetRandom();
    v_phi = fAzimuthalAngleAllCharges->GetRandom();
    v_eta = gRandom->Gaus(0.0,4.0);

    v_p[0] = v_pt*TMath::Cos(v_phi);
    v_p[1] = v_pt*TMath::Sin(v_phi);
    v_p[2] = v_pt*TMath::SinH(v_eta);
    v_E = TMath::Sqrt(TMath::Power(0.139,2) +
		      TMath::Power(v_p[0],2) +
		      TMath::Power(v_p[1],2) +
		      TMath::Power(v_p[2],2));

    v_y = 0.5*TMath::Log((v_E + v_p[2])/(v_E - v_p[2]));

    //Acceptance
    if((v_eta < fEtaMin) || (v_eta > fEtaMax)) continue;
    if((v_pt < fPtMin) || (v_pt > fPtMax)) continue;

    // fill charge vector
    chargeVector[0]->push_back(v_charge);
    chargeVector[1]->push_back(v_y);
    chargeVector[2]->push_back(v_eta);
    chargeVector[3]->push_back(v_phi);
    chargeVector[4]->push_back(v_p[0]);
    chargeVector[5]->push_back(v_p[1]);
    chargeVector[6]->push_back(v_p[2]);
    chargeVector[7]->push_back(v_pt);
    chargeVector[8]->push_back(v_E);
    
    if(fRunShuffling) {
      chargeVectorShuffle[0]->push_back(v_charge);
      chargeVectorShuffle[1]->push_back(v_y);
      chargeVectorShuffle[2]->push_back(v_eta);
      chargeVectorShuffle[3]->push_back(v_phi);
      chargeVectorShuffle[4]->push_back(v_p[0]);
      chargeVectorShuffle[5]->push_back(v_p[1]);
      chargeVectorShuffle[6]->push_back(v_p[2]);
      chargeVectorShuffle[7]->push_back(v_pt);
      chargeVectorShuffle[8]->push_back(v_E);
    }
    gNumberOfAcceptedParticles += 1;
	    
  }//positive particle loop

  //negative particles
  for(Int_t iNegCount = 0; iNegCount < nGeneratedNegative; iNegCount++) {
    v_charge = -1.0;

    v_pt = fPtSpectraAllCharges->GetRandom();
    v_phi = fAzimuthalAngleAllCharges->GetRandom();
    v_eta = gRandom->Gaus(0.0,4.0);

    v_p[0] = v_pt*TMath::Cos(v_phi);
    v_p[1] = v_pt*TMath::Sin(v_phi);
    v_p[2] = v_pt*TMath::SinH(v_eta);
    v_E = TMath::Sqrt(TMath::Power(0.139,2) +
		      TMath::Power(v_p[0],2) +
		      TMath::Power(v_p[1],2) +
		      TMath::Power(v_p[2],2));

    v_y = 0.5*TMath::Log((v_E + v_p[2])/(v_E - v_p[2]));

    //Acceptance
    if((v_eta < fEtaMin) || (v_eta > fEtaMax)) continue;
    if((v_pt < fPtMin) || (v_pt > fPtMax)) continue;

    // fill charge vector
    chargeVector[0]->push_back(v_charge);
    chargeVector[1]->push_back(v_y);
    chargeVector[2]->push_back(v_eta);
    chargeVector[3]->push_back(v_phi);
    chargeVector[4]->push_back(v_p[0]);
    chargeVector[5]->push_back(v_p[1]);
    chargeVector[6]->push_back(v_p[2]);
    chargeVector[7]->push_back(v_pt);
    chargeVector[8]->push_back(v_E);
    
    if(fRunShuffling) {
      chargeVectorShuffle[0]->push_back(v_charge);
      chargeVectorShuffle[1]->push_back(v_y);
      chargeVectorShuffle[2]->push_back(v_eta);
      chargeVectorShuffle[3]->push_back(v_phi);
      chargeVectorShuffle[4]->push_back(v_p[0]);
      chargeVectorShuffle[5]->push_back(v_p[1]);
      chargeVectorShuffle[6]->push_back(v_p[2]);
      chargeVectorShuffle[7]->push_back(v_pt);
      chargeVectorShuffle[8]->push_back(v_E);
    }
    gNumberOfAcceptedParticles += 1;
	    
  }//negative particle loop

  fBalance->CalculateBalance(gNumberOfAcceptedParticles,chargeVector);

  if(fRunShuffling) {
    // shuffle charges
    random_shuffle( chargeVectorShuffle[0]->begin(), chargeVectorShuffle[0]->end() );
    fShuffledBalance->CalculateBalance(gNumberOfAcceptedParticles,chargeVectorShuffle);
  }
}      

//________________________________________________________________________
void  AliAnalysisTaskToyModel::FinishTaskOutput() {
  //Printf("END BF");

  if (!fBalance) {
    Printf("ERROR: fBalance not available");
    return;
  }  
  if(fRunShuffling) {
    if (!fShuffledBalance) {
      Printf("ERROR: fShuffledBalance not available");
      return;
    }
  }
}

//________________________________________________________________________
void AliAnalysisTaskToyModel::Terminate(Option_t *) {
  // Draw result to the screen
  // Called once at the end of the query

  // not implemented ...

}
