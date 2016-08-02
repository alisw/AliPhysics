/// \class AliReducedHypTritEvent
/// \brief Reduced Event class
///
/// Reduced Event class for saving hypertriton information in a TTree for later analysis.
/// Each Event contains an array of v0, which all have a positive and a negative track.
///
/// \author Lukas Kreis <l.kreis@gsi.de>, GSI
/// \date Feb 17, 2016


#ifndef AliReducedHypTritEVENT_H
#define AliReducedHypTritEVENT_H

#ifndef ROOT_Object
#include "TObject.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include <TClonesArray.h>
#endif

class AliReducedHypTritTrack : public TObject {

  friend class AliAnalysisTaskHypTritEventTree;

 public:
  AliReducedHypTritTrack();
  ~AliReducedHypTritTrack();

  // Getters
  TLorentzVector P() const {return fP;}
  Double_t Dca() const {return fDca;}
  Double_t Phi() const {return fPhi;}
  Double_t Eta() const {return fEta;}
  Double_t Charge() const {return fCharge;}
  Double_t Dedx() const {return fDedx;}
  Double_t DedxSigma() const {return fDedxSigma;}
  Double_t TPCNcls() const {return fTpcNClusters;}

 private:
  TLorentzVector fP;   // 4 momentum of track
  Double_t fDca;        // DCA to prim vertex
  Double_t fDedx;             // specific energyloss in TPC of track
  Double_t fDedxSigma;        // dEdx sigma
  Double_t fCharge;           // sign of  track
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
  TVector3 Position() const {return fPosition;}
  Double_t M() const {return fM;}
  Double_t Pt() const {return fPt;}
  Double_t P() const {return fP;}
  Double_t SigmaD0() const {return fSigmaD0;}
  Float_t Dca() const {return fDCAv0;}
  Double_t Ct() const {return fDecayLength;}
  Float_t Cos() const {return fCosPointingAngle;}
  Double_t Chi2() const {return fChi2;}
  Bool_t MC() const {return fMCTruth;}
  Double_t Y() const {return fRapidity;}
  Int_t Charge() const {return fCharge;}
  AliReducedHypTritTrack *Pi() const {return fPiTrack;}
  AliReducedHypTritTrack *He() const {return fHeTrack;}

 private:
  TVector3 fPosition;         // Lorentzvector of v0
  Double_t fP;       // momentum of mother
  Double_t fPt;       //  transverse momentum of mother
  Double_t fM;  // reconstructed invariant mass of mother
  Float_t fDCAv0;           // DCA of daughters
  Double_t fSigmaD0;         // Sigma of tracks
  Float_t fCosPointingAngle;  // cosine of pointing angle of vertex
  Float_t fDecayLength;       // decay radius of mother particle
  Double_t fChi2;               // chi2 of V0
  Bool_t fMCTruth;              // Monte Carlo truth of mother type
  Double_t fRapidity;           // Rapidity of V0
  Int_t fCharge;                // anti or particle
  AliReducedHypTritTrack *fPiTrack;        // positive daughter of v0
  AliReducedHypTritTrack *fHeTrack;        // negative daughter of v0

  AliReducedHypTritV0(const AliReducedHypTritV0&);
  AliReducedHypTritV0 &operator = (const AliReducedHypTritV0&);
  ClassDef(AliReducedHypTritV0, 1);
};

class AliReducedHypTritEvent : public TObject {

  friend class AliAnalysisTaskHypTritEventTree;

 public:
  AliReducedHypTritEvent();
  ~AliReducedHypTritEvent();

  AliReducedHypTritV0 *V0(Int_t i) const
      {return (i < fNumberV0s ? (AliReducedHypTritV0*) fV0s->At(i) : 0x0);}
  Float_t Centrality() const {return fCentrality;}
  Int_t NumberV0s() const {return fNumberV0s;}
  Bool_t Trigger() const {return fTrigger;}
  Int_t RunNumber() const {return fRunNumber;}
  TVector3 VertexPosition() const {return fVertexPosition;}
  Int_t EventId() const {return fEventId;}

  void ClearEvent();

 private:
  Int_t fEventId;
  TVector3 fVertexPosition; // position of primary vertex
  Int_t fNumberV0s;            // number of v0s in event
  Float_t fCentrality;    // centrality of event
  Int_t fRunNumber;       // number of run
  TClonesArray *fV0s;     // array of v0s in event
  Int_t fTrigger;     // array of Triggers

  static TClonesArray *fgV0s;

  AliReducedHypTritEvent(const AliReducedHypTritEvent&);
  AliReducedHypTritEvent &operator = (const AliReducedHypTritEvent&);
  ClassDef(AliReducedHypTritEvent, 1);
};

#endif
