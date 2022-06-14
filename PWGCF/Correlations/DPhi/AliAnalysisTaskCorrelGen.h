#ifndef AliAnalysisTaskCorrelGen_H
#define AliAnalysisTaskCorrelGen_H

/**
 * \file AliAnalysisTaskCorrelGen.h
 * \brief Declaration of class AliAnalysisTaskCorrelGen
 *
 * The task selects candidates for K0s, Lambdas, AntiLambdas, Xi+ and Xi-
 * and calculates correlations with charged unidentified particles in phi and eta.
 * The charged unidentified particles are also taken as trigger particles to have a check.
 * The task works with AOD or ESD (with MC truth) events only
 * and contains also mixing for acceptance corrections.
 * Adopted from AliAnalysisTaskDiHadCorrelHighPt.
 */

/* Copyright(c) 1998-2016, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include "AliAnalysisTaskSE.h"
#include "AliEventPoolManager.h"
#include "AliEventCuts.h"

class AliEventPoolManager;
class THnSparse;
class AliESDtrack;
class AliMCParticle;
class THistManager;

class AliAnalysisTaskCorrelGen : public AliAnalysisTaskSE
{
public:
  AliAnalysisTaskCorrelGen();
  AliAnalysisTaskCorrelGen(const char *name);
  virtual ~AliAnalysisTaskCorrelGen();

  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *option);
  virtual void Terminate(Option_t *option);

  

  void Correlation(TObjArray *triggers, TObjArray *associated, Bool_t hh, Bool_t V0h, Float_t perc, Bool_t hV0);
  void CorrelationMixing(TObjArray *triggers, TObjArray *bgTracks, Float_t perc);
  void CorrelationMixinghV0(TObjArray *bgTracks, TObjArray *assocArray, Float_t perc);

  Double_t GetDeltaPhi(Double_t trigphi, Double_t assocphi) const;

  // Configuration
  void SetNmixedEvts(Int_t nmix) { fNMixEvts = nmix; }
  void SetPtTrigMin(Double_t var) { fPtTrigMin = var; }
  void SetPtAsocMin(Double_t var) { fPtAsocMin = var; }
  void SetNbinsPtTrig(Int_t nbins) { fNbinsPtTrig = nbins; }
  void SetNbinsPtAssoc(Int_t nbins) { fNbinsPtAssoc = nbins; }
  void SetNbinsDeltaPhi(Int_t bins) { fNbinsDPhi = bins; }
  void SetNbinsDeltaEta(Int_t bins) { fNbinsDEta = bins; }
  void SetNbinsEta(Int_t bins) { fNbinsEta = bins; }
  void SetNbinsPhi(Int_t nbins) { fNbinsPhi = nbins; }
  void SetNbinsVz(Int_t nbins) { fNbinsVz = nbins; }
  void SetDoCorrelV0h(Bool_t doUse) { fV0hCorr = doUse; }
  void SetDoCorrelhh(Bool_t doUse) { fhhCorr = doUse; }
  void SetDoCorrelhV0(Bool_t doUse) { fhV0Corr = doUse; }
  void SetFilterBit(Int_t filter) { fFilterBit = filter; }
  void SetDoMixingGen(Bool_t mix) { fMixGen = mix; }
  void SetPVCut(Int_t cut) { fPrimVtxCut = cut; }
  void SetESDtrackCuts(AliESDtrackCuts *const trcuts) { fESDTrkCuts = trcuts; }
  void SetIsAOD(Bool_t aodAnalysis) { fIsAOD = aodAnalysis; }
  void SetTestPions(Bool_t pions) { fTestPions = pions; }
  void SetMaxPtAssoc(Double_t max) { fPtAssocMax = max; }
  void SetMaxPtTrig(Double_t max) { fPtTrigMax = max; }
  void SetMaxDCAToVertexZ(Float_t max) { fMaxDCAtoVtxZ = max; }
  void SetMaxDCAToVertexXY(Float_t max) { fMaxDCAtoVtxXY = max; }
  void SetMinNCrossedRowsTPCprimtracks(Float_t min) { fMinNCrossedRowsTPCprimtracks = min; }
  void SetHighMultTriggerV0(Bool_t tr) { fhighMult = tr; }
  void SetHighMultTriggerSPD(Bool_t tr) { fhighMultSPD = tr; }

  void SetEventQAPlots(Bool_t pl) { fEvtCutsQAPlots = pl; }
  void SetIspp(Bool_t pp) { fpp = pp; }
  void SetNbinsMult(Int_t nBins) { fNbinsMult = nBins; }
  Int_t GetOrigParton(const AliMCParticle *mcTrack) const;
  void SetOnTheFlyMCTrain(Bool_t tr) { fOnFlyMC = tr; }
  AliAODTrack *SetAliAODTrack(Double_t theta, Double_t phi, Double_t pt, Short_t charge);

  AliEventCuts *fAliEventCuts; //!

private:
  AliMCEvent *fmcEvent; //! input MC event
  TList *fOutputList;   //! output list

  THnSparse *fHistMCCorr;   //!
  THnSparse *fHistMCMixGen; //!
  THnSparse *fHistNtrigGen; //!
  TH1D *fHistGenMult;       //!

  Int_t fMixTrk;                 // size of track buffer for event mixing
  AliEventPoolManager *fPoolMgr; //! event pool manager
  AliEventPool *fPoolMCGen;      //!
  Double_t fPtTrigMin;           //
  Double_t fPtAsocMin;           //
  Int_t fNMixEvts;               // number of minimum mixed events

  Int_t fNbinsPtTrig;  //
  Int_t fNbinsPtAssoc; //

  Bool_t fV0hCorr;   // enable to run V0-h correlations separately
  Bool_t fhhCorr;    // enable to run h-h correlations separately
  Bool_t fhV0Corr;   // enable to run h-V0 correlations separately
  Int_t fFilterBit;  // enable to vary filter bit for systematic studies
  Bool_t fIsAOD;     // anable to change between ESDs and AODs
  Int_t fNbinsDPhi;  // Number of DeltaPhi Bins in correlation function and in mixing
  Int_t fNbinsDEta;  // Number of DeltaEta Bins in correlation functionand in mixing
  Int_t fNbinsEta;   // Number of Eta bins for efficiency correction
  Int_t fNbinsPhi;   // Number of Phi bins for efficiency correction
  Bool_t fMixGen;    // enable to compute the correlation function from generated MC particles
  Int_t fNbinsVz;    // number of PV bins for mixing
  Int_t fPrimVtxCut; // PV position acceptance

  AliESDtrackCuts *fESDTrkCuts; //!
  Bool_t fTestPions;            // for testing MC trains

  Double_t fPtAssocMax; // maximal pt of associated Particle
  Double_t fPtTrigMax;  // maximal pt of trigger Particle

  Double_t fPV[3];

  Float_t fMaxDCAtoVtxZ;                 // DCA selection criterium for primary tracks in Z direction
  Float_t fMaxDCAtoVtxXY;                // DCA selection criterium for primary tracks in XY plane
  Float_t fMinNCrossedRowsTPCprimtracks; // TPC quality selection criterium for primary tracks

  TObjArray *fMCTrkSel;          //!  generated associated particles
  TObjArray *fMCGenTrkMix;       //! generated associated particles for Mixing
  TObjArray *fMCTrkTrigSel;      //! generated trigger charged hadrons
  TObjArray *fMCTrkV0Sel;        //! Generated V0 triggers
  TObjArray *fMCV0AssocSel;      //! Generated V0 assoc
  TObjArray *fMCAssocSel;        //! all reconstructed associated particles, with reconstructed pt,phi,eta values - for raw correlation function
  TObjArray *fselectedMCV0assoc; //! all reconstructed V0 as associated particles, with reconstructed pt,phi,eta values - for raw correlation function
  TObjArray *fMCTrigSel;         //! all reconstructed trigger particles, with reconstructed pt,phi,eta values - for raw correlation function
  TObjArray *fMCV0RecTrigSel;    //! All reconstructed V0 candidates for triggers with reconstructed pt,phi,eta values - for raw correlation function
  TObjArray *fTrkSel;            //! tracks for mixing
  TObjArray *fTrkAssocSel;       //!
  TObjArray *fTrigTrkSel;        //!
  TObjArray *fV0TrigSel;         //!
  TObjArray *fV0AssocSel;        //!
  Bool_t fhighMult;              // enable V0 high multiplicity trigger
  Bool_t fhighMultSPD;           // enable SPD high multiplicity trigger

  Bool_t fEvtCutsQAPlots; // enable to save event QA plots
  Bool_t fpp;             // flag pp / PbPb
  Int_t fNbinsMult;       // number of multiplicity bins

  Bool_t fOnFlyMC; // on the fly MC train usage

  AliAnalysisTaskCorrelGen(const AliAnalysisTaskCorrelGen &);            // not implemented
  AliAnalysisTaskCorrelGen &operator=(const AliAnalysisTaskCorrelGen &); // not implemented

  ClassDef(AliAnalysisTaskCorrelGen, 39);
};

/**
 * \class AliV0ChParticleGen
 * \brief Implementation of generated particles.
 *
 * This class is required for correlations calculation and event mixing
 */

class AliV0ChParticleGen : public AliVParticle
{
public:
  AliV0ChParticleGen(Float_t eta, Float_t phi, Float_t pt, Short_t candidate, Double_t mass)
      : fEta(eta), fPhi(phi), fpT(pt), fCandidate(candidate), fLabel(-1), fIDh(-1), fIDpos(-1), fIDneg(-1), fIDbach(0), fMass(mass), fCharge(0), fPz(0), fEnergie(0), fPhi1(0), fPt1(0), fEta1(0), fChar1(0), fPhi2(0), fPt2(0), fEta2(0), fChar2(0), fOriginalParton(-30)
  {
  }
  AliV0ChParticleGen(Float_t eta, Float_t phi, Float_t pt, Short_t candidate, Int_t label)
      : fEta(eta), fPhi(phi), fpT(pt), fCandidate(candidate), fLabel(label), fIDh(-1), fIDpos(-1), fIDneg(-1), fIDbach(0), fMass(0), fCharge(0), fPz(0), fEnergie(0), fPhi1(0), fPt1(0), fEta1(0), fChar1(0), fPhi2(0), fPt2(0), fEta2(0), fChar2(0), fOriginalParton(-30)
  {
  }
  AliV0ChParticleGen(Float_t eta, Float_t phi, Float_t pt, Short_t candidate, Int_t label, Int_t iDh, Short_t charge)
      : fEta(eta), fPhi(phi), fpT(pt), fCandidate(candidate), fLabel(label), fIDh(iDh), fIDpos(-1), fIDneg(-1), fIDbach(0), fMass(0), fCharge(charge), fPz(0), fEnergie(0), fPhi1(0), fPt1(0), fEta1(0), fChar1(0), fPhi2(0), fPt2(0), fEta2(0), fChar2(0), fOriginalParton(-30)
  {
  }
  AliV0ChParticleGen(Float_t eta, Float_t phi, Float_t pt, Short_t candidate, Int_t label, Int_t iDh, Int_t partontype)
      : fEta(eta), fPhi(phi), fpT(pt), fCandidate(candidate), fLabel(label), fIDh(iDh), fIDpos(-1), fIDneg(-1), fIDbach(0), fMass(0), fCharge(0), fPz(0), fEnergie(0), fPhi1(0), fPt1(0), fEta1(0), fChar1(0), fPhi2(0), fPt2(0), fEta2(0), fChar2(0), fOriginalParton(partontype)
  {
  }
  AliV0ChParticleGen(Float_t eta, Float_t phi, Float_t pt, Short_t candidate, Int_t label, Short_t charge, Double_t pz, Double_t energ)
      : fEta(eta), fPhi(phi), fpT(pt), fCandidate(candidate), fLabel(label), fIDh(-1), fIDpos(-1), fIDneg(-1), fIDbach(0), fMass(0), fCharge(charge), fPz(pz), fEnergie(energ), fPhi1(0), fPt1(0), fEta1(0), fChar1(0), fPhi2(0), fPt2(0), fEta2(0), fChar2(0), fOriginalParton(-30)
  {
  }
  AliV0ChParticleGen(Float_t eta, Float_t phi, Float_t pt, Short_t candidate, Int_t label, Int_t iDh, Short_t charge, Double_t pz, Double_t energ)
      : fEta(eta), fPhi(phi), fpT(pt), fCandidate(candidate), fLabel(label), fIDh(iDh), fIDpos(-1), fIDneg(-1), fIDbach(0), fMass(0), fCharge(charge), fPz(pz), fEnergie(energ), fPhi1(0), fPt1(0), fEta1(0), fChar1(0), fPhi2(0), fPt2(0), fEta2(0), fChar2(0), fOriginalParton(-30)
  {
  }
  AliV0ChParticleGen(Float_t eta, Float_t phi, Float_t pt, Short_t candidate, Int_t label, Int_t idpos, Int_t idneg, Double_t mass, Double_t phi1, Double_t pt1, Double_t eta1, Double_t char1, Double_t phi2, Double_t pt2, Double_t eta2, Double_t char2)
      : fEta(eta), fPhi(phi), fpT(pt), fCandidate(candidate), fLabel(label), fIDh(-1), fIDpos(idpos), fIDneg(idneg), fIDbach(0), fMass(mass), fCharge(0), fPz(0), fEnergie(0), fPhi1(phi1), fPt1(pt1), fEta1(eta1), fChar1(char1), fPhi2(phi2), fPt2(pt2), fEta2(eta2), fChar2(char2), fOriginalParton(-30)
  {
  }
  AliV0ChParticleGen(Float_t eta, Float_t phi, Float_t pt, Short_t candidate, Int_t label, Int_t idpos, Int_t idneg, Double_t mass)
      : fEta(eta), fPhi(phi), fpT(pt), fCandidate(candidate), fLabel(label), fIDh(-1), fIDpos(idpos), fIDneg(idneg), fIDbach(0), fMass(mass), fCharge(0), fPz(0), fEnergie(0), fPhi1(0), fPt1(0), fEta1(0), fChar1(0), fPhi2(0), fPt2(0), fEta2(0), fChar2(0), fOriginalParton(-30)
  {
  }
  AliV0ChParticleGen(Float_t eta, Float_t phi, Float_t pt, Short_t candidate, Int_t label, Int_t idpos, Int_t idneg, Double_t mass, Int_t partontype)
      : fEta(eta), fPhi(phi), fpT(pt), fCandidate(candidate), fLabel(label), fIDh(-1), fIDpos(idpos), fIDneg(idneg), fIDbach(0), fMass(mass), fCharge(0), fPz(0), fEnergie(0), fPhi1(0), fPt1(0), fEta1(0), fChar1(0), fPhi2(0), fPt2(0), fEta2(0), fChar2(0), fOriginalParton(partontype)
  {
  }
  AliV0ChParticleGen(Float_t eta, Float_t phi, Float_t pt, Short_t candidate, Double_t mass, Int_t idpos, Int_t idneg, Int_t idbach)
      : fEta(eta), fPhi(phi), fpT(pt), fCandidate(candidate), fLabel(0), fIDh(-1), fIDpos(idpos), fIDneg(idneg), fIDbach(idbach), fMass(mass), fCharge(0), fPz(0), fEnergie(0), fPhi1(0), fPt1(0), fEta1(0), fChar1(0), fPhi2(0), fPt2(0), fEta2(0), fChar2(0), fOriginalParton(-30)
  {
  }
  virtual ~AliV0ChParticleGen() {}

  // kinematics
  virtual Double_t Px() const
  {
    AliFatal("Not implemented");
    return 0;
  }
  virtual Double_t Py() const
  {
    AliFatal("Not implemented");
    return 0;
  }
  virtual Double_t Pz() const { return fPz; }
  virtual Double_t Pt() const { return fpT; }
  virtual Double_t P() const
  {
    AliFatal("Not implemented");
    return 0;
  }
  virtual Bool_t PxPyPz(Double_t[3]) const
  {
    AliFatal("Not implemented");
    return 0;
  }

  virtual Double_t Xv() const
  {
    AliFatal("Not implemented");
    return 0;
  }
  virtual Double_t Yv() const
  {
    AliFatal("Not implemented");
    return 0;
  }
  virtual Double_t Zv() const
  {
    AliFatal("Not implemented");
    return 0;
  }
  virtual Bool_t XvYvZv(Double_t[3]) const
  {
    AliFatal("Not implemented");
    return 0;
  }

  virtual Double_t OneOverPt() const
  {
    AliFatal("Not implemented");
    return 0;
  }
  virtual Double_t Phi() const { return fPhi; }
  virtual Double_t Theta() const
  {
    AliFatal("Not implemented");
    return 0;
  }

  virtual Double_t E() const { return fEnergie; }
  virtual Double_t M() const { return fMass; }

  virtual Double_t Eta() const { return fEta; }
  virtual Double_t Y() const
  {
    AliFatal("Not implemented");
    return 0;
  }

  virtual Short_t Charge() const { return fCharge; }
  virtual Int_t GetLabel() const { return fLabel; }
  // PID
  virtual Int_t PdgCode() const
  {
    AliFatal("Not implemented");
    return 0;
  }
  virtual const Double_t *PID() const
  {
    AliFatal("Not implemented");
    return 0;
  }

  virtual Short_t WhichCandidate() const { return fCandidate; }
  Int_t GetIDCh() const { return fIDh; }
  Int_t GetIDPos() const { return fIDpos; }
  Int_t GetIDNeg() const { return fIDneg; }
  Int_t GetIDBach() const { return fIDbach; }
  Short_t GetCharge1() const { return fChar1; }
  Short_t GetCharge2() const { return fChar2; }
  Double_t GetPhi1() const { return fPhi1; }
  Double_t GetPhi2() const { return fPhi2; }
  Double_t GetPt1() const { return fPt1; }
  Double_t GetPt2() const { return fPt2; }
  Double_t GetEta1() const { return fEta1; }
  Double_t GetEta2() const { return fEta2; }
  Int_t GetOrigalPartonType() const { return fOriginalParton; }

private:
  Double_t fEta;         // eta
  Double_t fPhi;         // phi
  Double_t fpT;          // pT
  Short_t fCandidate;    // V0 candidate: 1 - K0, 2 - Lambda, 3 - Antilambda, 4 - ChTrack
  Int_t fLabel;          // Label of MC particles
  Int_t fIDh;            // Label
  Int_t fIDpos;          // Label of possitive charged daughter
  Int_t fIDneg;          // Label of negative charged daughter
  Int_t fIDbach;         // Label of bachelor track
  Double_t fMass;        // mass
  Short_t fCharge;       // charge of the track
  Double_t fPz;          // pZ
  Double_t fEnergie;     // E
  Double_t fPhi1;        // phi of first daughter
  Double_t fPt1;         // pt of first daughter
  Double_t fEta1;        // eta of first daughter
  Double_t fChar1;       // charge of first daughter
  Double_t fPhi2;        // phi of second daughter
  Double_t fPt2;         // pt of second daughter
  Double_t fEta2;        // eta of second daughter
  Double_t fChar2;       // charge of second daughter
  Int_t fOriginalParton; // 0- quark; 1-gluon, for MC truth information

  ClassDef(AliV0ChParticleGen, 8)
};

#endif
