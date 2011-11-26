#ifndef AliIonPDGCodes_h
#define AliIonPDGCodes_h 
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include <TArrayI.h>
#include <TParticlePDG.h>
#include <TParticleClassPDG.h>

class TArrayI;
class TParticle;
class TDatabasePDG;

class AliIonPDGCodes : public TObject
{

public:
  
  AliIonPDGCodes();
  virtual ~AliIonPDGCodes() {;}
  AliIonPDGCodes(const AliIonPDGCodes &PDGCodes);
  AliIonPDGCodes& operator=(const AliIonPDGCodes& pdg);
  virtual void AddParticlesToPdgDataBase();
//  virtual void MapPDGGEant3Codes();
  virtual void SetPDGCode(Int_t i, Int_t val) {fPDGCode[i]=val;}
  Int_t IdFromPDG(Int_t pdg) const;
  Int_t PDGFromId(Int_t pdg) const;
  
protected:
  
  Int_t fNIon;
  Int_t fPDGCode[200]; 

  ClassDef(AliIonPDGCodes,1) 

};
#endif
