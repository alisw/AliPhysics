#ifndef ALIVZEROv5_H
#define ALIVZEROv5_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */


///////////////////////////////////////////////////
//                                               //
//  Manager and hits classes for set : VZERO     //
//                                   version 5   //
//                                   june 2004   //
//                                               //
///////////////////////////////////////////////////

#include "TLorentzVector.h" 
#include "AliVZERO.h"

class AliVZEROv5 : public AliVZERO {
  
public:
  AliVZEROv5();
  AliVZEROv5(const char *name, const char *title);
  virtual       ~AliVZEROv5() {}
  virtual void   AddHit(Int_t track, Int_t *vol, Float_t *hits); 
  virtual void   AddDigits(Int_t *tracks, Int_t *digits);
  virtual void   CreateGeometry();
  virtual void   CreateMaterials();
  virtual void   DrawModule() const;
  virtual void   Init();
  virtual void   MakeBranch(Option_t *option);
  virtual Int_t  IsVersion() const {return 3;}
  virtual void   StepManager();
  Int_t          GetCellId(Int_t *vol, Float_t *hits);
  
protected:
  Int_t          fCellId;        // Scintillator cell number from 0 to 95 
  TLorentzVector fTrackPosition; // Position of particle entering cell
  TLorentzVector fTrackMomentum; // Momentum of particle entering cell
  
private:  
  Float_t fLightYield;       // Lightyield in BC408   (93.75 eV per photon)
  Float_t fLightAttenuation; // LightAttenuation in fibers (0.05 per meter)
  Float_t fnMeters;          // Number of meters of fibers to PM
  Float_t fFibToPhot;        // Loss in Fibers - Photocathode Connection 
  
  ClassDef(AliVZEROv5,1)  //Class for VZERO version 5
};

#endif


