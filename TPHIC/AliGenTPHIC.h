#ifndef ALIGENTPHIC_H
#define ALIGENTPHIC_H
/* Copyright(c) 1998-2003, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */
// Event generator of two-photon processes
// in ultra-peripheral ion collisions
// Author: Yuri.Kharlov@cern.ch
// 15 April 2003

#include "AliGenMC.h"
class TPHICgen;
class AliPythia;

//-------------------------------------------------------------
class AliGenTPHIC : public AliGenMC
{
public:
  AliGenTPHIC();
  virtual ~AliGenTPHIC();
  void Generate();
  void Init();
  void SetDebug(Int_t debug) {fDebug=debug;}
  void SetEventListRange(Int_t eventFirst=-1, Int_t eventLast=-1);

  // Setters for TPHIC initial parameters

  void SetProcess     (Int_t   proc  =2    );
  void SetBeamEnergy  (Float_t energy=3500.);
  void SetBeamZ       (Int_t   z     =20   );
  void SetBeamA       (Int_t   a     =40   );
  void SetYggRange    (Float_t ymin  =-7., Float_t ymax=7.);
  void SetMggRange    (Float_t mmin  = 2., Float_t mmax=20.);
  void SetNgridY      (Int_t   ny    = 20  );
  void SetNgridM      (Int_t   nm    =100  );
  void SetLumFunName  (TString name  ="lum_ca_2_20.dat"  );
  void SetLumFunFlag  (Int_t   flag  =-1   );
  void SetKfFermion   (Int_t   kf    = 13  );
  void SetKfOnium     (Int_t   kf    =441  );
  void SetMassOnium   (Float_t mass        );
  void SetGGwidthOnium(Float_t width       );
  void SetKfVmesons   (Int_t kf=113, Int_t kf2=113);

  // Getters of TPHIC output parameters

  Float_t GetGGmass              ();
  Float_t GetGGrapidity          ();
  Float_t GetG1mass              ();
  Float_t GetG2mass              ();
  TClonesArray*  GetParticleList ();
  TLorentzVector MomentumRecNucl1();
  TLorentzVector MomentumRecNucl2();
  Float_t GetXSectionCurrent     ();
  Float_t GetXSection            ();
  Float_t GetXSectionError       ();
 protected:
  TPHICgen     *fTPHICgen;          //!generator TPHIC17
  AliPythia    *fPythia;            //!generator PYTHIA6
  TClonesArray *fParticles;         // Particle  List
  Int_t         fEvent;             //!internal event number
  Int_t         fDebug;             //!debug level
  Int_t         fDebugEventFirst;   //!First event to debug
  Int_t         fDebugEventLast;    //!Last  event to debug

 private:
  AliGenTPHIC(const AliGenTPHIC & gen);
  AliGenTPHIC & operator=(const AliGenTPHIC & gen);

  ClassDef(AliGenTPHIC,1)     // Generator of 2-photon processes in ultra-peripheral collisions
};
#endif
