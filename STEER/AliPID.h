#ifndef ALIPID_H
#define ALIPID_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

///////////////////////////////////////////////////////////////////////////////
///                                                                          //
/// particle id probability densities                                        //
///                                                                          //
///////////////////////////////////////////////////////////////////////////////

/* $Id$ */


#include <TObject.h>


class AliPID : public TObject {
 public:
  enum {
    kSPECIES = 5,    // Number of particle species recognized by the PID
    kSPECIESN = 10   // Number of charged+neutral particle species recognized by the PHOS/EMCAL PID
  };
  enum EParticleType {
    kElectron = 0, 
    kMuon = 1, 
    kPion = 2, 
    kKaon = 3, 
    kProton = 4, 
    kPhoton = 5, 
    kPi0 = 6, 
    kNeutron = 7, 
    kKaon0 = 8, 
    kEleCon = 9,
    kUnknown = 10
  };
  static Float_t       ParticleMass(Int_t iType) 
    {return fgkParticleMass[iType];};
  static const char*   ParticleName(Int_t iType) 
    {return fgkParticleName[iType];};
  static const char*   ParticleShortName(Int_t iType) 
    {return fgkParticleShortName[iType];};
  static const char*   ParticleLatexName(Int_t iType) 
    {return fgkParticleLatexName[iType];};
  static Int_t         ParticleCode(Int_t iType) 
    {return fgkParticleCode[iType];};

  AliPID();
  AliPID(const Double_t* probDensity, Bool_t charged = kTRUE);
  AliPID(const Float_t* probDensity, Bool_t charged = kTRUE);
  AliPID(const AliPID& pid);
  AliPID& operator = (const AliPID& pid);

  Double_t             GetProbability(EParticleType iType,
				      const Double_t* prior) const;
  Double_t             GetProbability(EParticleType iType) const;
  void                 GetProbabilities(Double_t* probabilities,
					const Double_t* prior) const;
  void                 GetProbabilities(Double_t* probabilities) const;
  EParticleType        GetMostProbable(const Double_t* prior) const;
  EParticleType        GetMostProbable() const;
  
  void                 SetProbabilities(const Double_t* probabilities,
                                        Bool_t charged = kTRUE);

  static void          SetPriors(const Double_t* prior,
				 Bool_t charged = kTRUE);
  static void          SetPrior(EParticleType iType, Double_t prior);

  AliPID&              operator *= (const AliPID& pid);

 private:

  void                 Init();

  Bool_t               fCharged;                   // flag for charged/neutral
  Double_t             fProbDensity[kSPECIESN];    // probability densities
  static Double_t      fgPrior[kSPECIESN];         // a priori probabilities

  static /*const*/ Float_t fgkParticleMass[kSPECIESN+1]; // particle masses
  static const char*   fgkParticleName[kSPECIESN+1]; // particle names
  static const char*   fgkParticleShortName[kSPECIESN+1]; // particle names
  static const char*   fgkParticleLatexName[kSPECIESN+1]; // particle names
  static const Int_t   fgkParticleCode[kSPECIESN+1]; // particle codes

  ClassDef(AliPID, 2)    // particle id probability densities
};


AliPID operator * (const AliPID& pid1, const AliPID& pid2);


#endif
