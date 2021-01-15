#ifndef AliUEHistograms_H
#define AliUEHistograms_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: AliUEHistograms.h 20164 2007-08-14 15:31:50Z morsch $ */

// encapsulates several AliUEHist objects for a full UE analysis plus additional control histograms

#include "TNamed.h"
#include "AliUEHist.h"
#include "TMath.h"
#include "THn.h" // in cxx file causes .../THn.h:257: error: conflicting declaration ‘typedef class THnT<float> THnF’

class AliVParticle;

class TList;
class TSeqCollection;
class TObjArray;
class TH1F;
class TH2F;
class TH3F;

class AliUEHistograms : public TNamed
{
 public:
  AliUEHistograms(const char* name = "AliUEHistograms", const char* histograms = "", const char* binning = 0);
  virtual ~AliUEHistograms();
  
  void Fill(Int_t eventType, Float_t zVtx, AliUEHist::CFStep step, AliVParticle* leading, TList* toward, TList* away, TList* min, TList* max);
  void FillCorrelations(Double_t centrality, Float_t zVtx, AliUEHist::CFStep step, TObjArray* particles, TObjArray* mixed = 0, Float_t weight = 1, Bool_t firstTime = kTRUE, Bool_t twoTrackCuts = kTRUE, Float_t bSign = 0, Float_t twoTrackEfficiencyCutValue = -1, Bool_t applyEfficiency = kFALSE);
  void Fill(AliVParticle* leadingMC, AliVParticle* leadingReco);
  void FillEvent(Int_t eventType, Int_t step);
  void FillEvent(Double_t centrality, Int_t step);
  void FillTrackingEfficiency(TObjArray* mc, TObjArray* recoPrim, TObjArray* recoAll, TObjArray* recoPrimPID, TObjArray* recoAllPID, TObjArray* fake, Int_t particleType, Double_t centrality = 0, Double_t zVtx = 0);
  void FillFakePt(TObjArray* fake, Double_t centrality);
 
  void CopyReconstructedData(AliUEHistograms* from);
  void DeepCopy(AliUEHistograms* from);
  
  AliUEHist* GetUEHist(Int_t id);
  
  AliUEHist* GetNumberDensitypT() { return fNumberDensitypT; }
  AliUEHist* GetSumpT() { return fSumpT; }
  AliUEHist* GetNumberDensityPhi() { return fNumberDensityPhi; }
  
  void SetNumberDensitypT(AliUEHist* obj) { fNumberDensitypT = obj; }
  void SetSumpT(AliUEHist* obj) { fSumpT = obj; }
  void SetNumberDensityPhi(AliUEHist* obj) { fNumberDensityPhi = obj; }
  
  void SetRunNumber(Long64_t runNumber) { fRunNumber = runNumber; }
  
  void SetEfficiencyCorrectionTriggers(THnF* hist)   { fEfficiencyCorrectionTriggers = hist;   }
  void SetEfficiencyCorrectionAssociated(THnF* hist) { fEfficiencyCorrectionAssociated = hist; }
  
  TH2F* GetCorrelationpT()  { return fCorrelationpT; }
  TH2F* GetCorrelationEta() { return fCorrelationEta; }
  TH2F* GetCorrelationPhi() { return fCorrelationPhi; }
  TH2F* GetCorrelationR()   { return fCorrelationR; }
  TH2F* GetCorrelationLeading2Phi() { return fCorrelationLeading2Phi; }
  TH2F* GetCorrelationMultiplicity() { return fCorrelationMultiplicity; }
  TH3F* GetYield() { return fYields; }
  TH2F* GetInvYield() { return fInvYield2; }
  TH3F* GetYieldEtaPhiPT() { return fYieldsEtaPhiPT; }
  
  TH2F* GetEventCount()     { return fEventCount; }
  TH3F* GetEventCountDifferential() { return fEventCountDifferential; }
  TH1F* GetVertexContributors() { return fVertexContributors; }
  TH1F* GetCentralityDistribution() { return fCentralityDistribution; }
  TH2F* GetCentralityCorrelation() { return fCentralityCorrelation; }
  Long64_t GetRunNumber() { return fRunNumber; }
  Int_t GetMergeCount() { return fMergeCount; }
  TH3F* GetTwoTrackDistance(Int_t i) { return fTwoTrackDistancePt[i]; }
  TH2F* GetControlConvResoncances() { return fControlConvResoncances; }
  Bool_t GetWeightPerEvent() { return fWeightPerEvent; }
  
  void Correct(AliUEHistograms* corrections);
  
  void SetEtaRange(Float_t etaMin, Float_t etaMax);
  void SetPtRange(Float_t ptMin, Float_t ptMax);
  void SetPartSpecies(Int_t species);
  void SetZVtxRange(Float_t min, Float_t max);
  void SetContaminationEnhancement(TH1F* hist);
  void SetCombineMinMax(Bool_t flag);
  void SetTrackEtaCut(Float_t value);
  void SetWeightPerEvent(Bool_t flag);
  void SetSelectCharge(Int_t selectCharge) { fSelectCharge = selectCharge; }
  void SetSelectTriggerCharge(Int_t selectCharge) { fTriggerSelectCharge = selectCharge; }
  void SetSelectAssociatedCharge(Int_t selectCharge) { fAssociatedSelectCharge = selectCharge; }
  void SetTriggerRestrictEta(Float_t eta) { fTriggerRestrictEta = eta; }
  void SetEtaOrdering(Bool_t flag) { fEtaOrdering = flag; }
  void SetPairCuts(Float_t conversions, Float_t resonances) { fCutConversionsV = conversions; fCutK0sV = resonances; fCutLambdaV = resonances; }
  void SetCutOnPhi(bool cutOnPhi) { if (cutOnPhi) fCutPhiV = 0.005; }
  void SetCutOnRho(bool cutOnRho) { if (cutOnRho) fCutRhoV = 0.005; }
  void SetCutOnK0s(Float_t cutOnK0sV) { fCutK0sV = cutOnK0sV; }
  void SetCutOnLambda(Float_t cutOnLambdaV) { fCutLambdaV = cutOnLambdaV; }
  void SetCutOnPhi(Float_t cutOnPhiV) { fCutPhiV = cutOnPhiV; }
  void SetCutOnRho(Float_t cutOnRhoV) { fCutRhoV = cutOnRhoV; }
  void SetCustomCut(Float_t cutCustomMass, Float_t cutCustomFirst, Float_t cutCustomSecond, Float_t cutCustomV) { fCutCustomMass = cutCustomMass; fCutCustomFirst = cutCustomFirst; fCutCustomSecond = cutCustomSecond; fCutCustomV = cutCustomV; }
  
  void SetRejectResonanceDaughters(Int_t value) { fRejectResonanceDaughters = value; }
  void SetOnlyOneEtaSide(Int_t flag)    { fOnlyOneEtaSide = flag; }
  void SetOnlyOneAssocEtaSide(Int_t flag)    { fOnlyOneAssocEtaSide = flag; }
  void SetPtOrder(Bool_t flag) { fPtOrder = flag; }
  void SetTwoTrackCutMinRadius(Float_t min) { fTwoTrackCutMinRadius = min; }

  void SetCheckEventNumberInCorrelation(Bool_t val) { fCheckEventNumberInCorrelation = val; }
  void ExtendTrackingEfficiency(Bool_t verbose = kFALSE);
  void Reset();

  AliUEHistograms(const AliUEHistograms &c);
  AliUEHistograms& operator=(const AliUEHistograms& c);
  virtual void Copy(TObject& c) const;

  virtual Long64_t Merge(TCollection* list);
  void Scale(Double_t factor);
  
protected:
  void FillRegion(AliUEHist::Region region, Float_t zVtx, AliUEHist::CFStep step, AliVParticle* leading, TList* list, Int_t multiplicity);
  Int_t CountParticles(TList* list, Float_t ptMin);
  void DeleteContainers();
  inline Float_t GetInvMassSquared(Float_t pt1, Float_t eta1, Float_t phi1, Float_t pt2, Float_t eta2, Float_t phi2, Float_t m0_1, Float_t m0_2);
  inline Float_t GetInvMassSquaredCheap(Float_t pt1, Float_t eta1, Float_t phi1, Float_t pt2, Float_t eta2, Float_t phi2, Float_t m0_1, Float_t m0_2);
  inline Float_t GetDPhiStar(Float_t phi1, Float_t pt1, Float_t charge1, Float_t phi2, Float_t pt2, Float_t charge2, Float_t radius, Float_t bSign);
  
  static const Int_t fgkUEHists; // number of histograms

  AliUEHist* fNumberDensitypT;   // d^2N/dphideta vs pT,lead
  AliUEHist* fSumpT;             // d^2 sum(pT)/dphideta vs pT,lead
  AliUEHist* fNumberDensityPhi;  // d^2N/dphideta vs delta phi,lead (in pT,lead bins)
  
  TH2F* fCorrelationpT;         // pT,lead: true vs reco
  TH2F* fCorrelationEta;        // #eta,lead; true vs reco
  TH2F* fCorrelationPhi;        // #phi,lead; true vs reco
  TH2F* fCorrelationR;          // R = sqrt(delta eta^2 + delta phi^2) (true vs reco) vs pT,lead,MC
  TH2F* fCorrelationLeading2Phi;// delta phi (true vs reco) vs pT,lead,MC
  TH2F* fCorrelationMultiplicity; // number of mc particls vs reco particles (for pT > 0.5 GeV/c)
  TH3F* fYields;                // centrality vs pT vs eta
  TH3F* fYieldsEtaPhiPT;        // pT vs eta vs phi
  TH2F* fInvYield2; 		// invariant yield as cross check of tracking
  
  TH2F* fEventCount;            // event count as function of step, (for pp: event type (plus additional step -1 for all events without vertex range even in MC)) (for PbPb: centrality)
  TH3F* fEventCountDifferential;// event count as function of leading pT, step, event type
  
  TH1F* fVertexContributors;    // number of contributors to the vertex
  TH1F* fCentralityDistribution; // distribution of the variable used for centrality selection
  TH2F* fCentralityCorrelation;  // centrality vs multiplicity
  
  TH3F* fITSClusterMap;          // its cluster map vs centrality vs pT
  
  TH3F* fTwoTrackDistancePt[2];    // control histograms for two-track efficiency study: dphi*_min vs deta (0 = before cut, 1 = after cut)
  TH2F* fControlConvResoncances; // control histograms for cuts on conversions and resonances
  
  THnF* fEfficiencyCorrectionTriggers;   // if non-0 this efficiency correction is applied on the fly to the filling for trigger particles. The factor is multiplicative, i.e. should contain 1/efficiency
  THnF* fEfficiencyCorrectionAssociated;   // if non-0 this efficiency correction is applied on the fly to the filling for associated particles. The factor is multiplicative, i.e. should contain 1/efficiency
  
  Int_t fSelectCharge;           // (un)like sign selection when building correlations: 0: no selection; 1: unlike sign; 2: like sign
  Int_t fTriggerSelectCharge;    // select charge of trigger particle
  Int_t fAssociatedSelectCharge; // select charge of associated particle
  Float_t fTriggerRestrictEta;   // restrict eta range for trigger particle (default: -1 [off])
  Bool_t fEtaOrdering;           // activate eta ordering to prevent shape distortions. see FillCorrelation for the details
  Float_t fCutConversionsV;        // cut on conversions (inv mass)
  Float_t fCutK0sV;              // cut on K0s (inv mass)
  Float_t fCutLambdaV;           // cut on Lambda (inv mass)
  Float_t fCutPhiV;              // cut on Phi (inv mass)
  Float_t fCutRhoV;              // cut on Rho (inv mass)
  Float_t fCutCustomMass;       // user-defined inv mass value
  Float_t fCutCustomFirst;      // user-defined mass of the 1st particle
  Float_t fCutCustomSecond;     // user-defined mass of the 2nd particle
  Float_t fCutCustomV;          // cut on user-defined value (inv mass)
  Int_t fRejectResonanceDaughters; // reject all daughters of all resonance candidates (1: test method (cut at m_inv=0.9); 2: k0; 3: lambda)
  Int_t fOnlyOneEtaSide;       // decides that only trigger particle from one eta side are considered (0 = all; -1 = negative, 1 = positive)
  Int_t fOnlyOneAssocEtaSide;       // decides that only associated particle from one eta side are considered (0 = all; -1 = negative, 1 = positive)
  Bool_t fWeightPerEvent;	// weight with the number of trigger particles per event
  Bool_t fPtOrder;		// apply pT,a < pT,t condition
  Float_t fTwoTrackCutMinRadius; // min radius for TTR cut

  Bool_t fCheckEventNumberInCorrelation; // do not correlate two particles from the same event (only works for AliBasicParticles)

  Long64_t fRunNumber;           // run number that has been processed
  
  Int_t fMergeCount;		// counts how many objects have been merged together
  
  ClassDef(AliUEHistograms, 33)  // underlying event histogram container
};

Float_t AliUEHistograms::GetDPhiStar(Float_t phi1, Float_t pt1, Float_t charge1, Float_t phi2, Float_t pt2, Float_t charge2, Float_t radius, Float_t bSign)
{ 
  //
  // calculates dphistar
  //
  
  Float_t dphistar = phi1 - phi2 - charge1 * bSign * TMath::ASin(0.075 * radius / pt1) + charge2 * bSign * TMath::ASin(0.075 * radius / pt2);
  
  static const Double_t kPi = TMath::Pi();
  
  // circularity
//   if (dphistar > 2 * kPi)
//     dphistar -= 2 * kPi;
//   if (dphistar < -2 * kPi)
//     dphistar += 2 * kPi;
  
  if (dphistar > kPi)
    dphistar = kPi * 2 - dphistar;
  if (dphistar < -kPi)
    dphistar = -kPi * 2 - dphistar;
  if (dphistar > kPi) // might look funny but is needed
    dphistar = kPi * 2 - dphistar;
  
  return dphistar;
}

Float_t AliUEHistograms::GetInvMassSquared(Float_t pt1, Float_t eta1, Float_t phi1, Float_t pt2, Float_t eta2, Float_t phi2, Float_t m0_1, Float_t m0_2)
{
  // calculate inv mass squared
  // same can be achieved, but with more computing time with
  /*TLorentzVector photon, p1, p2;
  p1.SetPtEtaPhiM(triggerParticle->Pt(), triggerEta, triggerParticle->Phi(), 0.510e-3);
  p2.SetPtEtaPhiM(particle->Pt(), eta[j], particle->Phi(), 0.510e-3);
  photon = p1+p2;
  photon.M()*/
  
  Float_t tantheta1 = 1e10;
  
  if (eta1 < -1e-10 || eta1 > 1e-10)
  {
    Float_t expTmp = TMath::Exp(-eta1);
    tantheta1 = 2.0 * expTmp / ( 1.0 - expTmp*expTmp);
  }
  
  Float_t tantheta2 = 1e10;
  if (eta2 < -1e-10 || eta2 > 1e-10)
  {
    Float_t expTmp = TMath::Exp(-eta2);
    tantheta2 = 2.0 * expTmp / ( 1.0 - expTmp*expTmp);
  }
  
  Float_t e1squ = m0_1 * m0_1 + pt1 * pt1 * (1.0 + 1.0 / tantheta1 / tantheta1);
  Float_t e2squ = m0_2 * m0_2 + pt2 * pt2 * (1.0 + 1.0 / tantheta2 / tantheta2);
  
  Float_t mass2 = m0_1 * m0_1 + m0_2 * m0_2 + 2 * ( TMath::Sqrt(e1squ * e2squ) - ( pt1 * pt2 * ( TMath::Cos(phi1 - phi2) + 1.0 / tantheta1 / tantheta2 ) ) );
  
//   Printf(Form("%f %f %f %f %f %f %f %f %f", pt1, eta1, phi1, pt2, eta2, phi2, m0_1, m0_2, mass2));
  
  return mass2;
}

Float_t AliUEHistograms::GetInvMassSquaredCheap(Float_t pt1, Float_t eta1, Float_t phi1, Float_t pt2, Float_t eta2, Float_t phi2, Float_t m0_1, Float_t m0_2)
{
  // calculate inv mass squared approximately
  
  Float_t tantheta1 = 1e10;
  
  if (eta1 < -1e-10 || eta1 > 1e-10)
  {
    Float_t expTmp = 1.0-eta1+eta1*eta1/2-eta1*eta1*eta1/6+eta1*eta1*eta1*eta1/24;
    tantheta1 = 2.0 * expTmp / ( 1.0 - expTmp*expTmp);
  }
  
  Float_t tantheta2 = 1e10;
  if (eta2 < -1e-10 || eta2 > 1e-10)
  {
    Float_t expTmp = 1.0-eta2+eta2*eta2/2-eta2*eta2*eta2/6+eta2*eta2*eta2*eta2/24;
    tantheta2 = 2.0 * expTmp / ( 1.0 - expTmp*expTmp);
  }
  
  Float_t e1squ = m0_1 * m0_1 + pt1 * pt1 * (1.0 + 1.0 / tantheta1 / tantheta1);
  Float_t e2squ = m0_2 * m0_2 + pt2 * pt2 * (1.0 + 1.0 / tantheta2 / tantheta2);
  
  // fold onto 0...pi
  Float_t deltaPhi = TMath::Abs(phi1 - phi2);
  while (deltaPhi > TMath::TwoPi())
    deltaPhi -= TMath::TwoPi();
  if (deltaPhi > TMath::Pi())
    deltaPhi = TMath::TwoPi() - deltaPhi;
  
  Float_t cosDeltaPhi = 0;
  if (deltaPhi < TMath::Pi()/3)
    cosDeltaPhi = 1.0 - deltaPhi*deltaPhi/2 + deltaPhi*deltaPhi*deltaPhi*deltaPhi/24;
  else if (deltaPhi < 2*TMath::Pi()/3)
    cosDeltaPhi = -(deltaPhi - TMath::Pi()/2) + 1.0/6 * TMath::Power((deltaPhi - TMath::Pi()/2), 3);
  else
    cosDeltaPhi = -1.0 + 1.0/2.0*(deltaPhi - TMath::Pi())*(deltaPhi - TMath::Pi()) - 1.0/24.0 * TMath::Power(deltaPhi - TMath::Pi(), 4);
  
  Float_t mass2 = m0_1 * m0_1 + m0_2 * m0_2 + 2 * ( TMath::Sqrt(e1squ * e2squ) - ( pt1 * pt2 * ( cosDeltaPhi + 1.0 / tantheta1 / tantheta2 ) ) );
  
//   Printf(Form("%f %f %f %f %f %f %f %f %f", pt1, eta1, phi1, pt2, eta2, phi2, m0_1, m0_2, mass2));
  
  return mass2;
}

#endif
