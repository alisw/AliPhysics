#ifndef ALIGENCRT_H
#define ALIGENCRT_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#include "AliGenerator.h"
#include "AliCRTConstants.h"

class TClonesArray;
class TF1;

class AliGenCRT : public AliGenerator {
 public:
  AliGenCRT();
  AliGenCRT(Int_t npart);
  AliGenCRT(const AliGenCRT& gen);
  virtual ~AliGenCRT();

  AliGenCRT& operator= (const AliGenCRT& gen);

  virtual void Init();
  virtual void Generate();
  virtual void SetPart(Int_t part) {fIpart = part;}

  void SetMode(ECRMode mode) {fCRMode = mode;}
  const TString*  GetMode() const {return fCRModeName;}

  void SetZenithalAngleRange(Float_t min,Float_t max=0) {fZenithMin=min;fZenithMax=max;}
  void SetAzimuthalAngleRange(Float_t min, Float_t max=0) {fAzimuthMin=min;fAzimuthMax=max;}

  void SetGridRange(Int_t nx,Float_t xwidth, Int_t nz, Float_t zwidth);
  const Float_t GetMomentumResolution() const {return fPResolution;}

  void SetMomentumDistrubutionFunction(TF1 *func) {fMomentumDist=func;}
  void SetZenithalDistributionFunction(TF1 *func) {fZenithDist = func;}
  void SetMomentumResolution(Float_t res=1.) {fPResolution=res;}

  const Float_t GetMomentum() const;
  const Float_t GetZenithAngle(Float_t mom) const;

  // The following methods are for testing pourpuses
  TF1* GetMomentumDistibution() const {return fMomentumDist;}
  TF1* GetUnfoldedDistribution() const {return fUnfoldedMomentumDist;}

  TClonesArray* GetArray() const {return fPDist;}

 protected:
  void InitApWeightFactors();
  void InitMomentumGeneration();
  void InitZenithalAngleGeneration();
  void GenerateOneMuonBundle();
  void GenerateOneSingleMuon(Bool_t withFlatMomentum=kFALSE);

 private:
  Int_t    fIpart;              //! Particle type.
  ECRMode  fCRMode;             //! Cosmic muons generation method flag
  TString* fCRModeName;         //! Cosmic muons generation mode name

  Float_t  fXwidth;             //! X width of the grid
  Int_t    fNx;                 //! Number of divisions in x
  Float_t  fZwidth;             //! Z widht of the  grid
  Int_t    fNz;                 //! Number of divisions in z
  Bool_t   fMuonGrid;           //! Flag for method (Muon-bundles) checkout

  Float_t  fZenithMin;          //! Minimum zenithal angle.
  Float_t  fZenithMax;          //! Maximum zenithal angle.

  Float_t  fAzimuthMin;         //! Minimum azimuthal angle.
  Float_t  fAzimuthMax;         //! Maximum azimuthal angle.

  Float_t  fPRange;             //! Cosmic muon momentum range width in GeVs.
  Float_t  fPResolution;        //! Momentum resolution in GeVs.

  TArrayF* fAp;                 //! a(p) correction factors for the ang. dist.

  TF1*     fMomentumDist;       //! Function to generate the momentum dist.
  TF1*     fUnfoldedMomentumDist; //!
  TF1*     fZenithDist;         //! Function to generate the zenith angle dist.

  TClonesArray* fPDist;         //! Array of fZenithDist, to be used by a(p).

  ClassDef(AliGenCRT, 1) // Generator for AliCRT class
};
#endif // ALIGENCRT_H
