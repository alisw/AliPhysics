/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
/*
 * plesae see source file for more details
 */
#ifndef ALIESDV0KINECUTS_H
#define ALIESDV0KINECUTS_H

#include <TObject.h>

class AliESDv0;
class AliESDEvent;
class AliVEvent;
class AliESDtrack;
class AliVTrack;
class AliKFParticle;
class AliKFVertex;

class AliESDv0KineCuts : public TObject{
 public:
  enum{ // Reconstructed V0
    kUndef = -1,
      kGamma = 0,
      kK0 = 1,
      kLambda = 2,
      kALambda = 3
      };
  
  AliESDv0KineCuts();
  virtual ~AliESDv0KineCuts();

  AliESDv0KineCuts(const AliESDv0KineCuts &ref);
  AliESDv0KineCuts &operator=(const AliESDv0KineCuts &ref);

  Bool_t ProcessV0(AliESDv0* const v0, Int_t &pdgV0, Int_t &pdgP, Int_t &pdgN);
  Bool_t ProcessV0(AliESDv0* const v0, Int_t &pdgP, Int_t &pdgN);

  Int_t  PreselectV0(AliESDv0* const v0);

  Bool_t CaseGamma(AliESDv0* const v0, Int_t &pdgV0, Int_t &pdgP, Int_t &pdgN);
  Bool_t CaseK0(AliESDv0* const v0, Int_t &pdgV0, Int_t &pdgP, Int_t &pdgN);
  Bool_t CaseLambda(AliESDv0* const v0, Int_t &pdgV0, Int_t &pdgP, Int_t &pdgN, Int_t id);

  Bool_t V0CutsCommon(const AliESDv0 * const v0);
  Bool_t SingleTrackCuts(AliESDv0 * const v0);
  void   Armenteros(AliESDv0* const v0, Float_t val[2]);
  Bool_t CheckSigns(AliESDv0* const v0);

  void   SetEvent(AliESDEvent* const event);
  void   SetEvent(AliVEvent* const event);
  void   SetPrimaryVertex(AliKFVertex* const v) { fPrimaryVertex = v; };

  // setter functions for cut values
  void   SetGammaCutChi2NDF(Float_t val)  { fGcutChi2NDF = val; };
  void   SetGammaCutCosPoint(Float_t *val) { 
    fGcutCosPoint[0] = val[0];
    fGcutCosPoint[1] = val[1];
  };
  void   SetGammaCutDCA(Float_t *val){
    fGcutDCA[0] = val[0];
    fGcutDCA[1] = val[1];
  };
  void   SetGammaCutVertexR(Float_t *val){
    fGcutVertexR[0] = val[0];
    fGcutVertexR[1] = val[1];
  };
  void   SetGammaCutPsiPair(Float_t *val){
    fGcutPsiPair[0] = val[0];
    fGcutPsiPair[1] = val[1];
  };
  void   SetGammaCutInvMass(Float_t val){
    fGcutInvMass = val;
  };

  Double_t PsiPair(AliESDv0* const v0);

  Bool_t GetConvPosXY(AliESDtrack * const ptrack, AliESDtrack * const ntrack, Double_t convpos[2]);
  Bool_t GetHelixCenter(AliESDtrack * const track, Double_t b, Int_t charge, Double_t center[2]);

 protected:
  void Copy(TObject &ref) const;

 private:

  AliKFParticle *CreateMotherParticle(const AliVTrack* const pdaughter, const AliVTrack* const ndaughter, Int_t pspec, Int_t nspec);

 private:
  AliESDv0              *fV0;             // current V0 candidate
  AliESDEvent           *fEvent;          // current event
  AliKFVertex           *fPrimaryVertex;  // primary vertex

  // gamma cut values
  Float_t               fGcutChi2NDF;      // Chi2NF cut value for the AliKFparticle gamma
  Float_t               fGcutCosPoint[2];  // cos of the pointing angle [min, max]
  Float_t               fGcutDCA[2];       // DCA between the daughter tracks [min, max]
  Float_t               fGcutVertexR[2];   // radius of the conversion point [min, max]
  Float_t               fGcutPsiPair[2];   // value of the psi pair cut [min, max]
  Float_t               fGcutInvMass;      // upper value on the gamma invariant mass
  // K0 cut values
  Float_t               fK0cutChi2NDF;     // Chi2NF cut value for the AliKFparticle K0
  Float_t               fK0cutCosPoint[2]; // cos of the pointing angle [min, max]
  Float_t               fK0cutDCA[2];      // DCA between the daughter tracks [min, max]
  Float_t               fK0cutVertexR[2];  // radius of the decay point [min, max]
  Float_t               fK0cutInvMass[2];  // invariant mass window
  // Lambda & anti-Lambda cut values
  Float_t               fLcutChi2NDF;      // Chi2NF cut value for the AliKFparticle K0
  Float_t               fLcutCosPoint[2];  // cos of the pointing angle [min, max]
  Float_t               fLcutDCA[2];       // DCA between the daughter tracks [min, max]
  Float_t               fLcutVertexR[2];   // radius of the decay point [min, max]
  Float_t               fLcutInvMass[2];   // invariant mass window
  
  
  ClassDef(AliESDv0KineCuts, 0);

};

#endif
