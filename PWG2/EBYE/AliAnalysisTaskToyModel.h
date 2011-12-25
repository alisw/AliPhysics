#ifndef ALIANALYSISTASKTOYMODEL_CXX
#define ALIANALYSISTASKTOYMODEL_CXX

// Analysis task for a simple toy model, currently used for the 
// balance function
// Authors: Panos Christakoglou@nikhef.nl

class TList;
class TH1F;
class TH2F;
class TF1;

class AliBalance;

#include "AliAnalysisTaskSE.h"
#include "AliBalance.h"


class AliAnalysisTaskToyModel : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskToyModel(const char *name = "AliAnalysisTaskToyModel");
  virtual ~AliAnalysisTaskToyModel(); 
  
  virtual void   Init();
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   FinishTaskOutput();
  virtual void   Terminate(Option_t *);
  
  void SetAnalysisObject(AliBalance *const analysis) {
    fBalance         = analysis;
  }
  void SetShufflingObject(AliBalance *const analysisShuffled) {
    fRunShuffling = kTRUE;
    fShuffledBalance = analysisShuffled;
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

  //Acceptance filter
  void SetAcceptanceParameterization(TF1 *parameterization) {
    fAcceptanceParameterization = parameterization;}

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
  void SetQuadrangularFlowForAllCharges(Double_t v4) {
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
  void SetQuadrangularFlowForPions(Double_t v4) {
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
  void SetQuadrangularFlowForKaons(Double_t v4) {
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
  void SetQuadrangularFlowForProtons(Double_t v4) {
    fQuandrangularFlowProtons = v4;}
  void SetPentangularFlowForProtons(Double_t v5) {
    fPentangularFlowProtons = v5;}
  //============Toy model: List of setters============//

 private:
  AliBalance *fBalance; //BF object
  Bool_t fRunShuffling;//run shuffling or not
  AliBalance *fShuffledBalance; //BF object (shuffled)
  TList *fList; //fList object
  TList *fListBF; //fList object
  TList *fListBFS; //fList object

  TH1F *fHistEventStats; //event stats

  //Toy model input
  Double_t fTotalMultiplicityMean; //mean for the total multiplicity
  Double_t fTotalMultiplicitySigma; //sigma for the total multiplicity
  Double_t fNetChargeMean; //mean for the net charge
  Double_t fNetChargeSigma; //sigma for the net charge
  Double_t fPtMin; //pt min for acceptance
  Double_t fPtMax; //pt max for acceptance
  Double_t fEtaMin; //eta min for acceptance
  Double_t fEtaMax; //eta max for acceptance

  Bool_t fUseAcceptanceParameterization; //flag acceptance parameterization
  TF1 *fAcceptanceParameterization; //acceptance parameterization

  Bool_t   fUseAllCharges; //use all charges
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

  AliAnalysisTaskToyModel(const AliAnalysisTaskToyModel&); // not implemented
  AliAnalysisTaskToyModel& operator=(const AliAnalysisTaskToyModel&); // not implemented
  
  ClassDef(AliAnalysisTaskToyModel, 1); // example of analysis
};

#endif
