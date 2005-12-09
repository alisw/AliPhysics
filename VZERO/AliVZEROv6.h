#ifndef ALIVZEROv6_H
#define ALIVZEROv6_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */


/////////////////////////////////////////////////////
//                                                 //
//  Manager and hits classes for set :   VZERO     //
//                                     version 6   //
//                                  september 2005 //
//                                                 //
/////////////////////////////////////////////////////

#include "TLorentzVector.h" 
#include "AliVZERO.h"

class AliVZEROv6 : public AliVZERO {
  
public:
  AliVZEROv6();
  AliVZEROv6(const char *name, const char *title);
  virtual       ~AliVZEROv6() {};
  virtual void   AddHit(Int_t track, Int_t *vol, Float_t *hits); 
  virtual void   AddDigits(Int_t *tracks, Int_t *digits);
  virtual void   CreateGeometry();
  virtual void   BuildGeometry();
  virtual void   CreateMaterials();
  virtual void   DrawModule() const;
  virtual void   DrawGeometry();
  virtual void   Init();
  virtual void   MakeBranch(Option_t *option);
  virtual Int_t  IsVersion() const {return fVersion;};
  virtual void   StepManager();
  Int_t          GetCellId(Int_t *vol, Float_t *hits);
  
protected:
  Int_t          fCellId;        // Scintillator cell number from 0 to 95 
  TLorentzVector fTrackPosition; // Position of particle entering cell
  TLorentzVector fTrackMomentum; // Momentum of particle entering cell
  
private: 

// Parameters related to geometry :
// V0 part in front of muon arm absorber 

  Float_t  fV0CHeight1, fV0CHeight2, fV0CHeight3, fV0CHeight4; 
  Float_t  fV0CRMin, fV0CRBox;
  Float_t  fV0CLidThickness;
  Float_t  fV0CCellThickness;
  Float_t  fV0CBoxThickness; 
  Float_t  fV0COffsetFibers;

// V0 part on the other side with respect to IP

  Float_t  fV0AHeight1, fV0AHeight2, fV0AHeight3, fV0AHeight4; 
  Float_t  fV0ARMin;
  Float_t  fV0ACellThickness;
  
// Parameters related to light production :
 
  Float_t fLightYield;       // Lightyield in BC408   (93.75 eV per photon)
  Float_t fLightAttenuation; // LightAttenuation in fibers (0.05 per meter)
  Float_t fnMeters;          // Number of meters of fibers to PM
  Float_t fFibToPhot;        // Loss in Fibers - Photocathode Connection 
  Int_t   fVersion;          // Version number == IsVersion
  
  ClassDef(AliVZEROv6,1)  //Class for VZERO version 6
};

#endif


