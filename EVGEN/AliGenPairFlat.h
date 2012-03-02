#ifndef ALIGENPAIRFLAT_H
#define ALIGENPAIRFLAT_H


// Generator for particle pairs in a preset
// kinematic range 
// ranges can be set for invariant mass, pair pT, pair rapidity
// and pair azimuth

// Comments and suggestions: markus.konrad.kohler@cern.ch


#include "AliGenerator.h"

class TF1;

class AliGenPairFlat : public AliGenerator
{
 public:
  AliGenPairFlat();
  virtual ~AliGenPairFlat();
  virtual void Generate();
  virtual void Init();

  virtual void SetPairNPart(Int_t npart) {fPairNpart=npart;}
  virtual void SetPairYRange(Float_t Ymin, Float_t Ymax)
      {fPairYMin = Ymin; fPairYMax = Ymax;}
  virtual void SetPairPhiRange(Float_t phimin, Float_t phimax)
      {fPairPhiMin = phimin; fPairPhiMax = phimax;}
  virtual void SetPairPtRange(Float_t ptmin, Float_t ptmax)
      {fPairPtMin = ptmin; fPairPtMax = ptmax;}
  virtual void SetPairMassRange(Float_t massmin, Float_t massmax)
      {fPairMassMin = massmin; fPairMassMax = massmax;}
  virtual void SetLegPdg(Int_t pdg1, Int_t pdg2)
      {fLegPdg1 = pdg1; fLegPdg2 = pdg2;}
  virtual void SetPolarization(Float_t Alpha)
      {fAlpha = Alpha;}
  void SetDebug(Int_t debug) {fDebug=debug;}

protected:

  Int_t   fPairNpart;  		//Number of generated pairs
  Float_t fPairYMin;  		// Minimum eta 
  Float_t fPairYMax;  		// Maximum eta
  Float_t fPairPhiMin;  	// Minimum phi 
  Float_t fPairPhiMax;  	// Maximum phi
  Float_t fPairPtMin;  		// Minimum pt 
  Float_t fPairPtMax;  		// Maximum pt
  Float_t fPairMassMin;		// Minimum mass 
  Float_t fPairMassMax;   	// Maximum mass
  Int_t   fLegPdg1;		// pdg code of first daughter
  Int_t   fLegPdg2;		// pdg code of second daughter
  Float_t fAlpha;		// Polarization factor
  Int_t   fDebug;   		// debug level
  TF1    *fPol;			// Polarization function

  Bool_t Decay(const TLorentzVector& mother, TLorentzVector &dau1, TLorentzVector &dau2 , TF1* fPol);


  private:
  AliGenPairFlat(const AliGenPairFlat &pair);
  AliGenPairFlat & operator=(const AliGenPairFlat & pair);

  ClassDef(AliGenPairFlat,1) // Flat random pair generator
};

#endif
