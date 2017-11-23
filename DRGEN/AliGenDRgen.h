#ifndef ALIGENDRGEN_H
#define ALIGENDRGEN_H
/* Copyright(c) 1998-2013, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */


#include "AliGenMC.h"
class TDRgen;

//-------------------------------------------------------------
class AliGenDRgen : public AliGenMC
{
public:
  AliGenDRgen();
  virtual ~AliGenDRgen();
  void Generate();
  void Init();
  void SetDebug(Int_t debug) {fDebug=debug;}
  void SetEventListRange(Int_t eventFirst=-1, Int_t eventLast=-1);

  // Setters for initial parameters

  void SetProcess     (Int_t   proc1  =2, Int_t   proc2  =2    );
  void SetBeamEnergy  (Double_t energy=3500.);
  void SetF2Polarization  (Double_t *polarization);  // Set polarization of f2(1270) meson: f2 polarization array: |D0|^2, |D-|^2, |D+|^2, |D--|^2, |D++|^2, phase(D-,D0), phase(D--,D0), phase(D++,D+). Valid only for processes 9 and 10.
  void SetMassRange(Float_t min, Float_t max);
  // Getters of output parameters

  TClonesArray*  GetParticleList ();

 protected:
  TDRgen     *fTDRgen;          //!generator DRgen
  TClonesArray *fParticles;         // Particle  List
  Int_t         fEvent;             //!internal event number
  Int_t         fDebug;             //!debug level
  Int_t         fDebugEventFirst;   //!First event to debug
  Int_t         fDebugEventLast;    //!Last  event to debug
  Bool_t        kMassRange;         //!cut on mass of centrally produced system - flag
  Float_t     fMassMin;        //Minimum mass of centrally produced system
  Float_t     fMassMax;        //Maximum mass of centrally produced system
   

 private:
  AliGenDRgen(const AliGenDRgen & gen);
  AliGenDRgen & operator=(const AliGenDRgen & gen);

  ClassDef(AliGenDRgen,1)     // Generator of 2-pomeron processes in HE pp collisions
};
#endif
