/// \class AliReducedHypTritEvent
/// \brief Reduced Event class
///
/// Reduced Event class for saving hypertriton information in a TTree for later analysis.
/// Each Event contains an array of v0, which all have a positive and a negative track.
///
/// \author Lukas Kreis <l.kreis@gsi.de>, GSI
/// \date Feb 17, 2016


#ifndef ALIREDUCEDHYPTRITEVENT_H
#define ALIREDUCEDHYPTRITEVENT_H

#ifndef ROOT_Object
#include "TObject.h"
#include "TLorentzVector.h"
#include "TLorentzVector.h"
#include <TClonesArray.h>
#endif

class AliReducedHypTritTrack : public TObject {

  friend class AliAnalysisTaskHypTritEventTree;

 public:
  AliReducedHypTritTrack();
  ~AliReducedHypTritTrack();

  // Getters
  Double_t Pt() const {return fPt;}
  TLorentzVector Momentum() const {return fMomentum;}
  Double_t DCAtoPrim() const {return fDCAtoPrim;}
  Double_t P() const {return fP;}
  Double_t Phi() const {return fPhi;}
  Double_t Eta() const {return fEta;}
  Double_t Charge() const {return fSign;}
  Double_t TPCsignal() const {return fDedx;}
  Double_t GetTPCNcls() const {return fTpcNClusters;}

 private:
  Double_t fPt;                // transverse momentum of track
  TLorentzVector fMomentum;   // momentum of track
  Double_t fDCAtoPrim;        // DCA to prim vertex
  Double_t fP;                // total momentum of  track
  Double_t fDedx;             // specific energyloss in TPC of track
  Double_t fSign;             // sign of  track
  Double_t fEta;              // eta of track
  Double_t fPhi;              // phi of track
  Double_t fTpcNClusters;     // number of clusters


AliReducedHypTritTrack(const AliReducedHypTritTrack&);
AliReducedHypTritTrack &operator = (const AliReducedHypTritTrack&);
ClassDef(AliReducedHypTritTrack, 1)
};

class AliReducedHypTritV0 : public TObject {

  friend class AliAnalysisTaskHypTritEventTree;

 public:
  AliReducedHypTritV0();
  ~AliReducedHypTritV0();

  // Getters
  TLorentzVector Position() const {return fSecVertexPos;}
  Double_t Mass() const {return fMotherInvMass;}
  Double_t Pt() const {return fMotherPt;}
  Double_t P() const {return fMotherP;}
  Double_t SigmaD0() const {return fSigmaD0;}
  Float_t GetDcaV0Daughters() const {return fDCAv0;}
  Float_t GetD() const {return fDCAtoPrim;}
  Float_t GetRr() const {return fDecayRadius;}
  Float_t GetV0CosineOfPointingAngle() const {return fCosPointingAngle;}
  Bool_t GetOnFlyStatus() const {return fOnFlyStatus;}
  Double_t GetChi2V0() const {return fChi2;}
  Bool_t MCTruth() const {return fMCTruth;}
  Double_t GetRapidity() const {return fRapidity;}
  Bool_t IsAntiParticle() const {return fAntiParticle;}
  AliReducedHypTritTrack *GetPosTrack() const {return fPosTrack;}
  AliReducedHypTritTrack *GetNegTrack() const {return fNegTrack;}

 private:
  TLorentzVector fSecVertexPos;         // Lorentzvector of v0
  Double_t fDCAtoPrim;       // distance of v0 to primary vertex
  Double_t fMotherP;       // momentum of mother
  Double_t fMotherPt;       //  transverse momentum of mother
  Double_t fMotherInvMass;  // reconstructed invariant mass of mother
  Float_t fDCAv0;           // DCA of daughters
  Double_t fSigmaD0;         // Sigma of tracks
  Float_t fCosPointingAngle;  // cosine of pointing angle of vertex
  Float_t fDecayRadius;       // decay radius of mother particle
  Bool_t fOnFlyStatus;          // reconstruction of v0
  Double_t fChi2;               // chi2 of V0
  Bool_t fMCTruth;           // Monte Carlo truth of mother type
  Double_t fRapidity;           // Rapidity of V0
  Bool_t fAntiParticle;         // anti or particle
  AliReducedHypTritTrack *fPosTrack;        // positive daughter of v0
  AliReducedHypTritTrack *fNegTrack;        // negative daughter of v0

  AliReducedHypTritV0(const AliReducedHypTritV0&);
  AliReducedHypTritV0 &operator = (const AliReducedHypTritV0&);
  ClassDef(AliReducedHypTritV0, 1);
};

class AliReducedHypTritEvent : public TObject {

  friend class AliAnalysisTaskHypTritEventTree;

 public:
  AliReducedHypTritEvent();
  ~AliReducedHypTritEvent();

  AliReducedHypTritV0 *GetV0(Int_t i) const 
      {return (i < fNV0s ? (AliReducedHypTritV0*) fV0s->At(i) : 0x0);}
  Float_t GetCentrality() const {return fCentrality;}
  Int_t GetNV0s() const {return fNV0s;}
  Bool_t GetTrigger(Int_t i) const {return fTrigger[i];}
  Int_t GetRunNumber() const {return fRunNumber;}
  TLorentzVector VertexPosition() const {return fPrimVertexPos;}

  void ClearEvent();

 private:

  TLorentzVector fPrimVertexPos; // position of primary vertex
  Int_t fNV0s;            // number of v0s in event
  Float_t fCentrality;    // centrality of event
  Int_t fRunNumber;       // number of run
  TClonesArray *fV0s;     // array of v0s in event
  Bool_t fTrigger[3];     // array of Triggers

  static TClonesArray *fgV0s;

  AliReducedHypTritEvent(const AliReducedHypTritEvent&);
  AliReducedHypTritEvent &operator = (const AliReducedHypTritEvent&);
  ClassDef(AliReducedHypTritEvent, 1);
};

#endif


