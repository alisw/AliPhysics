#ifndef ROOT_DRGEN
#define ROOT_DRGEN
/* Copyright(c) 1998-2013, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include "TGenerator.h"

class TDRgen : public TGenerator {

public:
  TDRgen();
  virtual ~TDRgen();

  void Initialize    ();
  void GenerateEvent ();
  void Finish        ();

  // Setters for COMMON /PP2INI/
  //all these setters can be declared  const
  //because they don't change the _object_
  void SetIPROC    (Int_t   iproc1, Int_t iproc2   ) const;
  void SetSqrtS   (Double_t sqrts                  ) const;
  void SetProcXsec (Double_t *procXsec             ) const;
  void SetAMIN     (Double_t amin                  ) const;
  void SetAMAX     (Double_t amax                  ) const;
  void SetF2Polarization (Double_t *polarization   ) const;    //f2 polarization array: |D0|^2, |D-|^2, |D+|^2, |D--|^2, |D++|^2, phase(D-,D0), phase(D--,D0), phase(D++,D+)

  // Getters for COMMON /PP2EVNT/
  Int_t GetIevent  () const; // current event number
  Int_t GetIproc  () const;    //current ecent process type
  Int_t GetNpart () const;     // number of final-state particles
  Int_t GetIParticleNumber () const;      // number of particles in stack 

  void ImportParticles (TClonesArray * particles);
 
  ClassDef(TDRgen,1)  //Interface to DRgen Event Generator
};

#endif
