#ifndef AliVZEROv7_H
#define AliVZEROv7_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */


/////////////////////////////////////////////////////
//                                                 //
//  Manager and hits classes for set :   VZERO     //
//                                     version 7   //
//                                     April 2006  //
//                                                 //
/////////////////////////////////////////////////////

#include "TLorentzVector.h" 
#include "AliVZERO.h"

class AliVZEROv7 : public AliVZERO {
  
public:
  AliVZEROv7();
  AliVZEROv7(const char *name, const char *title);
  virtual       ~AliVZEROv7() {};
  virtual void   AddHit(Int_t track, Int_t *vol, Float_t *hits); 
  virtual void   AddDigits(Int_t *tracks, Int_t *digits);
  virtual void   CreateGeometry();
  virtual void   AddAlignableVolumes() const;
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

// V0C part in front of muon arm absorber 
// V0C Parameters related to geometry: 
  Double_t  fV0CHeight1, fV0CHeight2, fV0CHeight3, fV0CHeight4; // Heights of V0C elements
  Double_t  fV0CRMin, fV0CRBox;  // Min and max radii of V0C box
  Double_t  fV0CLidThickness;    // Thickness of V0C box lid
  Double_t  fV0CCellThickness;   // Thickness of V0C cell
  Double_t  fV0CBoxThickness;    // Thickness of V0C box
  Double_t  fV0COffsetFibers;    // Z offsets to output fibers

// V0C Parameters related to light production:
  Float_t fV0CLightYield;        // Lightyield in BC408   (93.75 eV per photon)
  Float_t fV0CLightAttenuation;  // LightAttenuation in fibers (0.05 per meter)
  Float_t fV0CnMeters;           // Number of meters of fibers to PM
  Float_t fV0CFibToPhot;         // Loss in Fibers - Photocathode Connection 

// V0A Parameters related to geometry:
  Double_t fV0AR0, fV0AR1, fV0AR2, fV0AR3, fV0AR4, fV0AR5, fV0AR6;	// Radius of V0A
  Double_t fV0ASciWd, fV0APlaWd, fV0APlaAl, fV0AOctWd, fV0AFraWd; 	// Thickness of elements
  Double_t fV0AOctH1, fV0AOctH2, fV0ABasHt;				// Height of elements
  Double_t fV0AFibRd;							// Radius of Fiber
  Double_t fV0APlaEx;                                                   // Extension of plates to basis
  Double_t fV0APMBWd, fV0APMBHt, fV0APMBTh, fV0APMBWdW, fV0APMBHtW;     // Parameters for Photo-Multiplier
  Double_t fV0APMBAng, fV0APMBThW, fV0APMTR1, fV0APMTR2, fV0APMTR3;     // Parameters for Photo-Multiplier
  Double_t fV0APMTR4, fV0APMTH, fV0APMTB;                               // Parameters for Photo-Multiplier
  Float_t fV0AnMeters;                                                  // Must be calculated depending on each ring
  
// V0A Parameters related to light production:
  Double_t fV0ALightYield;       // Lightyield in BC404
  Double_t fV0ALightAttenuation; // LightAttenuation in fibers
  Double_t fV0AFibToPhot;        // Loss in Fibers - Photocathode Connection

  Int_t   fVersion;              // Version number == IsVersion


  ClassDef(AliVZEROv7,1)  // Class for VZERO version 7
};

#endif
