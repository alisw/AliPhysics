#ifndef ALIPHOS_H
#define ALIPHOS_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

////////////////////////////////////////////////
//  Manager and hits classes for PHOS     //
//  Version 2                                 //
//  Author M. Volkov, RRC KI                  //
////////////////////////////////////////////////

// --- ROOT system ---
#include <TArray.h> 
#include <TRandom.h> 

// --- galice header files ---
#include "AliDetector.h"
#include "AliHit.h"

class AliPHOSv2 : public AliDetector{

protected:
  Float_t fXtlSize[3]; // PWO crystal dimensions
  Float_t fWrapThickness; // Thickness of Tyvek wrapper
  Float_t fPINSize[3]; // PIN diode dimensions
  Float_t fCPVThickness; // CPV thickness
  Float_t fPHOSFoam[3]; // Outer foam cover dimensions
  Float_t fPHOStxwall[3]; // Textolit box dimensions
  Float_t fPHOSAir[3]; // Inner air filled volume dimensions
  Float_t fRadius[2]; // Distances from IP to outer cover and to Xtal surface
  Float_t fPHOSextra[10]; // Assorted geometrical parameters
  Float_t fNphi; // Number of crystal units in X (phi) direction
  Float_t fNz; // Number of crystal units in Z direction
  Float_t fNModules; // Number of modules constituing PHOS
  Float_t fPHOSAngle[4]; // Position angles of modules
  
public:
  AliPHOSv2(void);
  AliPHOSv2(const char *name, const char *title="");
  virtual ~AliPHOSv2(void);

  virtual Int_t IsVersion(void) const {return 2;}

  virtual void DefPars(void);

  virtual void BuildGeometry(void);
  virtual void CreateGeometry(void);
  virtual void CreateMaterials(void);

  virtual void Init(void);

  virtual void StepManager(void);

  virtual void AddHit(Int_t track, Int_t *vol, Float_t *hits);

  virtual Float_t GetCrystalSize(Int_t n) const {return fXtlSize[n];}
  virtual Float_t GetWrapThickness(void) const {return fWrapThickness;}
  virtual Float_t GetPINSize(Int_t n) const {return fPINSize[n];}
  virtual Float_t GetCPVThickness(void) const {return fCPVThickness;}
  virtual Float_t GetPHOSFoam(Int_t n) const {return fPHOSFoam[n];}
  virtual Float_t GetPHOStxwall(Int_t n) const {return fPHOStxwall[n];}
  virtual Float_t GetPHOSAir(Int_t n) const {return fPHOSAir[n];}
  virtual Float_t GetRadius(Int_t n) const {return fRadius[n];}
  virtual Float_t GetPHOSextra(Int_t n) const {return fPHOSextra[n];}
  virtual Float_t GetNphi(void) const {return fNphi;}
  virtual Float_t GetNz(void) const {return fNz;}
  virtual Float_t GetNModules(void) const {return fNModules;}
  virtual Float_t &GetModuleAngle(Int_t n) {return fPHOSAngle[n];}

ClassDef(AliPHOSv2,1)  // Hits manager for PHOS, version 2

};

//////////////////////////////////////////////////////////////////////////////

class AliPHOShitv2 : public AliHit{

protected:
  Int_t     fVolume[4];  //array of volumes. This is not GEANT NUMBV(), it is (BOX,LAYER,ROW,COLUMN) array.
  Float_t   fELOS;       //Energy deposited
  
public:
  AliPHOShitv2(void) {;}
  AliPHOShitv2(Int_t shunt, Int_t track, Int_t *vol, Float_t *hits);
  AliPHOShitv2(AliPHOShitv2 const &rValue){*this=rValue;}
  virtual ~AliPHOShitv2(void) {;}
  
  Int_t GetVolume(Int_t i) const {return fVolume[i];}
  Float_t GetEnergy(void) const {return fELOS;}
  
  Bool_t operator==(AliPHOShitv2 const &rValue) const;
  AliPHOShitv2 const operator+(AliPHOShitv2 const &rValue) const;

ClassDef(AliPHOShitv2,1)  // Hits object for PHOS

};

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef ALIPHOS_H
