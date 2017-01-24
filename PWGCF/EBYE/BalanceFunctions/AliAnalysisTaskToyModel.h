#ifndef ALIANALYSISTASKTOYMODEL_CXX
#define ALIANALYSISTASKTOYMODEL_CXX

// Analysis task for a simple toy model, currently used for the 
// balance function
// Authors: Panos Christakoglou@nikhef.nl

class TList;
class TH1F;
class TH2F;
class TH3F;
class TF1;
class TLorentzVector;

class AliBalancePsi;
class AliEventPoolManager;

class AliAnalysisTaskToyModel : public TObject {
 public:
  AliAnalysisTaskToyModel();
  virtual ~AliAnalysisTaskToyModel(); 
  
  virtual void   Init();
  virtual void   CreateOutputObjects();
  virtual void   Run(Int_t nEvents);
  virtual void   FinishOutput();
  
  void SetDebugFlag() {fUseDebug = kTRUE;}

  void SetAnalysisObject(AliBalancePsi *const analysis) {
    fBalance         = analysis;
  }
  void SetShufflingObject(AliBalancePsi *const analysisShuffled) {
    fRunShuffling = kTRUE;
    fShuffledBalance = analysisShuffled;
  }
  void SetMixingObject(AliBalancePsi *const analysisMixed) {
    fRunMixing = kTRUE;
    fMixedBalance = analysisMixed;
  }

  //============Toy model: List of setters============//
  void SetTotalMultiplicity(Double_t mean, Double_t sigma) {
    fTotalMultiplicityMean = mean;
    fTotalMultiplicitySigma = sigma;}
  void SetNetCharge(Double_t mean, Double_t sigma) {
    fNetChargeMean = mean;
    fNetChargeSigma = sigma;}

  //Acceptance
  void SetKinematicsCutsMC(Double_t ptmin, Double_t ptmax,
                           Double_t etamin, Double_t etamax){
    fPtMin  = ptmin; fPtMax  = ptmax;
    fEtaMin = etamin; fEtaMax = etamax;
  }

  // Additional options (for integral studies)
  void SetSigmaGaussEta(Double_t sigmaGaussEta){
    fSigmaGaussEta = sigmaGaussEta;
  }
  void SetConstantEta(Double_t constantEta){
    fConstantEta = constantEta;
  }
  void SetFixPt(Double_t fixPt){
    fFixPt = fixPt;
  } 
  void SetFixedPositiveRatio(Bool_t fixedPositiveRatio=kTRUE){
    fFixedPositiveRatio = fixedPositiveRatio;
  }

  //Acceptance filter
  void SetAcceptanceParameterization(TF1 *parameterization) {
    fUseAcceptanceParameterization = kTRUE;
    fAcceptanceParameterization = parameterization;}

  //Acceptance - simulate detector effects/inefficiencies
  void SimulateDetectorEffects() {fSimulateDetectorEffects = kTRUE;}
  void SetNumberOfInefficientSectorsInPhi(Int_t numberOfInefficientSectors) {
    fNumberOfInefficientSectors = numberOfInefficientSectors;
    fInefficiencyFactorInPhi = 0.5;}
  void SetInefficiencyFactor(Double_t gInefficiencyFactorInPhi) {
    fInefficiencyFactorInPhi = gInefficiencyFactorInPhi;}
  void SetNumberOfDeadSectorsInPhi(Int_t numberOfDeadSectors) {
    fNumberOfDeadSectors = numberOfDeadSectors;}
  void EnableEfficiencyDropNearEtaEdges() {
    fEfficiencyDropNearEtaEdges = kTRUE;}
  void SimulateDetectorEffectsCorrection(TString filename, 
					 Int_t nCentralityBins, 
					 Double_t *centralityArrayForCorrections);

  //All charges
  void SetSpectraTemperatureForAllCharges(Double_t temperature) {
    fUseAllCharges = kTRUE;
    fTemperatureAllCharges = temperature;}
  void SetDirectedFlowForAllCharges(Double_t v1) {
    fUseAllCharges = kTRUE;
    fDirectedFlowAllCharges = v1;}
  void SetEllipticFlowForAllCharges(Double_t v2) {
    fUseAllCharges = kTRUE;
    fEllipticFlowAllCharges = v2;}
  void SetTriangularFlowForAllCharges(Double_t v3) {
    fUseAllCharges = kTRUE;
    fTriangularFlowAllCharges = v3;}
  void SetQuandrangularFlowForAllCharges(Double_t v4) {
    fUseAllCharges = kTRUE;
    fQuandrangularFlowAllCharges = v4;}
  void SetPentangularFlowForAllCharges(Double_t v5) {
    fUseAllCharges = kTRUE;
    fPentangularFlowAllCharges = v5;}

  //Pions
  void SetPionPercentage(Double_t percentage) {
    fPionPercentage = percentage;}
  void SetSpectraTemperatureForPions(Double_t temperature) {
    fTemperaturePions = temperature;}
  void SetDirectedFlowForPions(Double_t v1) {
    fDirectedFlowPions = v1;}
  void SetEllipticFlowForPions(Double_t v2) {
    fEllipticFlowPions = v2;}
  void SetTriangularFlowForPions(Double_t v3) {
    fTriangularFlowPions = v3;}
  void SetQuandrangularFlowForPions(Double_t v4) {
    fQuandrangularFlowPions = v4;}
  void SetPentangularFlowForPions(Double_t v5) {
    fPentangularFlowPions = v5;}
  
  //Kaons
  void SetKaonPercentage(Double_t percentage) {
    fKaonPercentage = percentage;}
  void SetSpectraTemperatureForKaons(Double_t temperature) {
    fTemperatureKaons = temperature;}
  void SetDirectedFlowForKaons(Double_t v1) {
    fDirectedFlowKaons = v1;}
  void SetEllipticFlowForKaons(Double_t v2) {
    fEllipticFlowKaons = v2;}
  void SetTriangularFlowForKaons(Double_t v3) {
    fTriangularFlowKaons = v3;}
  void SetQuandrangularFlowForKaons(Double_t v4) {
    fQuandrangularFlowKaons = v4;}
  void SetPentangularFlowForKaons(Double_t v5) {
    fPentangularFlowKaons = v5;}

  //Protons
  void SetProtonPercentage(Double_t percentage) {
    fProtonPercentage = percentage;}
  void SetSpectraTemperatureForProtons(Double_t temperature) {
    fTemperatureProtons = temperature;}
  void SetDirectedFlowForProtons(Double_t v1) {
    fDirectedFlowProtons = v1;}
  void SetEllipticFlowForProtons(Double_t v2) {
    fEllipticFlowProtons = v2;}
  void SetTriangularFlowForProtons(Double_t v3) {
    fTriangularFlowProtons = v3;}
  void SetQuandrangularFlowForProtons(Double_t v4) {
    fQuandrangularFlowProtons = v4;}
  void SetPentangularFlowForProtons(Double_t v5) {
    fPentangularFlowProtons = v5;}

  // particle selection: -1 no selection, 0 pions, 1 kaons, 2 protons
  void SelectParticleSpecies(Int_t iParticle){
    fSelectParticle = iParticle;
  }

  //Dynamical correlations
  void SetCorrelationPercentage(Double_t percentage) {
    fUseDynamicalCorrelations = kTRUE; 
    fDynamicalCorrelationsPercentage = percentage;
  }
  void SetCorrelationDeltaEta(Double_t deltaEta) {
    fUseDynamicalCorrelations = kTRUE; 
    fDynamicalCorrelationsDeltaEta = deltaEta;
  }
  void SetCorrelationDeltaPhi(Double_t deltaPhi) {
    fUseDynamicalCorrelations = kTRUE; 
    fDynamicalCorrelationsDeltaPhi = deltaPhi;
  }

  //Rapidity shift
  void SetRapidityShift() {
    fUseRapidityShift = kTRUE;
  }

  //Use rapidity instead of eta in correlation functions
  void SetUseRapidity() {
    fUseRapidity = kTRUE;
  }
  
  //Jet-like structures
  void SetUseJets() {fUseJets = kTRUE;}

  // Local charge conservation (LCC)
  void SetUseLCC() {fUseLCC = kTRUE;}//default values for sigma (pt=0.1,eta=0.5,phi=0.5)

  void SetUseLCC(Double_t sigmaPt, Double_t sigmaEta, Double_t sigmaPhi) {
    fUseLCC   = kTRUE;
    fSigmaPt  = sigmaPt; 
    fSigmaEta = sigmaEta; 
    fSigmaPhi = sigmaPhi; 
  }
  
  //============Toy model: List of setters============//

 private:
  void SetupEfficiencyMatrix();//setup the efficiency matrix

  Double_t GetTrackbyTrackCorrectionMatrix(Double_t vEta, 
					      Double_t vPhi, 
					      Double_t vPt, 
					      Short_t vCharge, 
					      Double_t gCentrality);//efficiency factor for the track

  Bool_t fUseDebug; //Debug flag

  AliBalancePsi *fBalance; //BF object
  Bool_t fRunShuffling;//run shuffling or not
  AliBalancePsi *fShuffledBalance; //BF object (shuffled)
  Bool_t fRunMixing;//run mixing or not
  AliBalancePsi *fMixedBalance; //BF object (mixed)
  AliEventPoolManager*     fPoolMgr;         //! event pool manager
  TList *fList; //fList object
  TList *fListBF; //fList object
  TList *fListBFS; //fList object (shuffling)
  TList *fListBFM; //fList object (mixing)

  TH1F *fHistEventStats; //event stats
  TH1F *fHistNumberOfAcceptedParticles; //number of accepted particles
  TH1F *fHistReactionPlane; //reaction plane angle
  TH1F *fHistEtaTotal; //pseudo-rapidity (full phase space)
  TH1F *fHistEta; //pseudo-rapidity (acceptance)
  TH2F *fHistEtaPhiPos; //eta-phi pos
  TH2F *fHistEtaPhiNeg; //eta-phi neg
  TH1F *fHistRapidity; //rapidity (acceptance)
  TH1F *fHistRapidityPions; //rapidity (acceptance)
  TH1F *fHistRapidityKaons; //rapidity (acceptance)
  TH1F *fHistRapidityProtons; //rapidity (acceptance)
  TH1F *fHistPhi; //phi (acceptance)
  TH1F *fHistPhiPions; //phi (acceptance)
  TH1F *fHistPhiKaons; //phi (acceptance)
  TH1F *fHistPhiProtons; //phi (acceptance)
  TH1F *fHistPt; //pt (acceptance)
  TH1F *fHistPtPions; //pt (acceptance)
  TH1F *fHistPtKaons; //pt (acceptance)
  TH1F *fHistPtProtons; //pt (acceptance)

  //Toy model input
  Double_t fTotalMultiplicityMean; //mean for the total multiplicity
  Double_t fTotalMultiplicitySigma; //sigma for the total multiplicity
  Double_t fNetChargeMean; //mean for the net charge
  Double_t fNetChargeSigma; //sigma for the net charge
  Double_t fPtMin; //pt min for acceptance
  Double_t fPtMax; //pt max for acceptance
  Double_t fEtaMin; //eta min for acceptance
  Double_t fEtaMax; //eta max for acceptance
  Double_t fSigmaGaussEta; //sigma for the Gaussian distribution of randomly produced particles (default = 4.0)
  Double_t fConstantEta; //eta value for a constant distribution of particles, if -1. then Gauss is used (default = -1)
  Double_t fFixPt; //fixed pT for unidentified particles (default = -1., i.e. randomly produced pT from fPtSpectraAllCharges)
  Bool_t fFixedPositiveRatio; //fix ratio of produced positive to negative particles (dafult = kFALSE, i.e. randomly produced ratio)
  

  //Acceptance parameterization
  Bool_t fUseAcceptanceParameterization; //flag acceptance parameterization
  TF1 *fAcceptanceParameterization; //acceptance parameterization

  //Simulate detector effects
  Bool_t fSimulateDetectorEffects;//simulate detector effects in pT
  Int_t fNumberOfInefficientSectors;//inefficient secotrs in phi
  Double_t fInefficiencyFactorInPhi;//efficiency factor < 1
  Int_t fNumberOfDeadSectors;//number of dead sectors
  Bool_t fEfficiencyDropNearEtaEdges;//efficiency drop in eta edges
  TH3F *fEfficiencyMatrix; //efficiency matrix in eta-pt-phi
  
  Bool_t fSimulateDetectorEffectsCorrection;//simulate detector effects as used for correction of data
  TH3F *fHistCorrectionPlus[101]; //correction matrix Plus
  TH3F *fHistCorrectionMinus[101]; //correction matrix minus
  Double_t fCentralityArrayForCorrections[101];//centrality array for correction
  Int_t fCentralityArrayBinsForCorrections;//number of centralitry bins
  Double_t fPtMinForCorrections;//only used for AODs
  Double_t fPtMaxForCorrections;//only used for AODs
  Double_t fPtBinForCorrections; //=================================correction
  Double_t fEtaMinForCorrections;//only used for AODs
  Double_t fEtaMaxForCorrections;//only used for AODs
  Double_t fEtaBinForCorrections; //=================================correction
  Double_t fPhiMinForCorrections;//only used for AODs
  Double_t fPhiMaxForCorrections;//only used for AODs
  Double_t fPhiBinForCorrections; //=================================correction

  //Kinematics
  Bool_t   fUseAllCharges; //use all charges
  Int_t    fSelectParticle; //select particle species (no selection -1 [default], pions 0, kaons 1, protons 2)
  Double_t fParticleMass; //particle mass
  TF1     *fPtSpectraAllCharges; //spectra for all charges
  Double_t fTemperatureAllCharges; //temperature for pt spectra
  Double_t fReactionPlane; //reaction plane angle
  TF1     *fAzimuthalAngleAllCharges; //azimuthal angle
  Double_t fDirectedFlowAllCharges; //directed flow value
  Double_t fEllipticFlowAllCharges; //elliptic flow value
  Double_t fTriangularFlowAllCharges; //triangular flow value
  Double_t fQuandrangularFlowAllCharges; //quadrangular flow value
  Double_t fPentangularFlowAllCharges; //pentangular flow value

  Double_t fPionPercentage; //percentage of pions
  Double_t fPionMass; //pion mass
  TF1     *fPtSpectraPions; //spectra for pions
  Double_t fTemperaturePions; //temperature for pt spectra
  TF1     *fAzimuthalAnglePions; //azimuthal angle for pions
  Double_t fDirectedFlowPions; //directed flow value
  Double_t fEllipticFlowPions; //elliptic flow value
  Double_t fTriangularFlowPions; //triangular flow value
  Double_t fQuandrangularFlowPions; //quadrangular flow value
  Double_t fPentangularFlowPions; //pentangular flow value

  Double_t fKaonPercentage; //percentage of kaons
  Double_t fKaonMass; //kaon mass
  TF1     *fPtSpectraKaons; //spectra for kaons
  Double_t fTemperatureKaons; //temperature for pt spectra
  TF1     *fAzimuthalAngleKaons; //azimuthal angle for kaons
  Double_t fDirectedFlowKaons; //directed flow value
  Double_t fEllipticFlowKaons; //elliptic flow value
  Double_t fTriangularFlowKaons; //triangular flow value
  Double_t fQuandrangularFlowKaons; //quadrangular flow value
  Double_t fPentangularFlowKaons; //pentangular flow value

  Double_t fProtonPercentage; //percentage of protons
  Double_t fProtonMass; //proton mass
  TF1     *fPtSpectraProtons; //spectra for protons
  Double_t fTemperatureProtons; //temperature for pt spectra
  TF1     *fAzimuthalAngleProtons; //azimuthal angle for protons
  Double_t fDirectedFlowProtons; //directed flow value
  Double_t fEllipticFlowProtons; //elliptic flow value
  Double_t fTriangularFlowProtons; //triangular flow value
  Double_t fQuandrangularFlowProtons; //quadrangular flow value
  Double_t fPentangularFlowProtons; //pentangular flow value

  Bool_t fUseDynamicalCorrelations; //Usage of dynamical correlations
  Double_t fDynamicalCorrelationsPercentage; //Percentage of correlations
  Double_t fDynamicalCorrelationsDeltaEta; //Gaussian Width in Eta for correlations
  Double_t fDynamicalCorrelationsDeltaPhi; //Gaussian Width in Eta correlations

  Bool_t fUseRapidityShift; //Usage of rapidity shift
  Double_t fRapidityShift;//shift in rapidity (only applicable if fUseAllCharges==kFALSE)
  Bool_t fUseRapidity;  //Usage of rapidity instead of eta in correlation functions
  TLorentzVector *vTmp;// temporary vector for calculating shifted observables
  TLorentzVector *vShift;// rapidity shift
  TLorentzVector *vBeam_p;// proton beam
  TLorentzVector *vBeam_Pb;// Pb beam

  Bool_t fUseJets;//Usage of jet-like structures
  TF1 *fPtAssoc;//pt of associated

  Bool_t fUseLCC;//Usage of Local Charge Conservation
  Double_t fSigmaPt;//sigma for LCC spread in pT
  Double_t fSigmaEta;//sigma for LCC spread in Eta
  Double_t fSigmaPhi;//sigma for LCC spread in Phi

  AliAnalysisTaskToyModel(const AliAnalysisTaskToyModel&); // not implemented
  AliAnalysisTaskToyModel& operator=(const AliAnalysisTaskToyModel&); // not implemented
  
  ClassDef(AliAnalysisTaskToyModel, 1); // example of analysis
};

#endif
