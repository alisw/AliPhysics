#ifndef ALIUPCMUONTRACK_H
#define ALIUPCMUONTRACK_H

//_____________________________________________________________________________
//    Class for UPC data
//    Author: Jaroslav Adam
//_____________________________________________________________________________

#include "TObject.h"

class TLorentzVector;

class AliUPCMuonTrack : public TObject
{
public:
  AliUPCMuonTrack();

  virtual ~AliUPCMuonTrack();

  void Clear(Option_t * /*option*/ ="");

  //Setters
  void SetPtEtaPhi(Double_t pt, Double_t eta, Double_t phi) {fPt=pt; fEta=eta; fPhi=phi;}

  void SetCharge(Short_t charge) {fCharge=charge;}

  void SetMatchTrigger(Int_t match) {fMatchTrigger=match;}
  void SetRAtAbsorberEnd(Double_t rabs) {fRabs=rabs;}
  void SetChi2perNDF(Double_t chi2) {fChi2perNDF=chi2;}
  void SetDCA(Double_t dca) {fDca=dca;}
  void SetPxDCA(Bool_t pdca) {fPdca=pdca;}

  Int_t MakeArrayInt(Int_t size);
  Int_t MakeArrayD(Int_t size);

  //Getters
  //virtual void GetPtEtaPhi(Double_t pt, Double_t eta, Double_t phi) const {pt=fPt; eta=fEta; phi=fPhi;}
  Double_t GetPt(void) const { return fPt; }
  Double_t GetEta(void) const { return fEta; }
  Double_t GetPhi(void) const { return fPhi; }
  void GetMomentum(TLorentzVector *v) const;

  Short_t GetCharge(void) const { return fCharge; }

  Int_t GetMatchTrigger(void) const { return fMatchTrigger; }
  Double_t GetRAtAbsorberEnd(void) const { return fRabs; }
  Double_t GetChi2perNDF(void) const { return fChi2perNDF; }
  Double_t GetDCA(void) const { return fDca; }
  Bool_t GetPxDCA(void) const { return fPdca; }

  TArrayI *GetArrayInt(void) const { return fArrayInt; }
  TArrayD *GetArrayD(void) const { return fArrayD; }

protected:
  AliUPCMuonTrack(const AliUPCMuonTrack &o);
  AliUPCMuonTrack &operator=(const AliUPCMuonTrack &o);

  Double_t fPt; // transversal momentum
  Double_t fEta; // pseudorapidity
  Double_t fPhi; // azimutal angle
  Short_t fCharge; // track charge
  Int_t fMatchTrigger; // muon trigger match
  Double_t fRabs; // transverse position r of the track at the end of the absorber
  Double_t fChi2perNDF; // chi2/NDF of momentum fit
  Double_t fDca; // Distance of Closest Approach in the vertex plane
  Bool_t fPdca; // pDCA by AliMuonTrackCuts
  TArrayI *fArrayInt; // extension of the muon track for other integer parameters
  TArrayD *fArrayD; // extension of the muon track for other double parameters

  const Double_t fkMuonMass; // mass of muon

  ClassDef(AliUPCMuonTrack,1);
};

#endif



















