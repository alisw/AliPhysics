#ifndef ALIITSUHIT_H
#define ALIITSUHIT_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

////////////////////////////////////////////////////////////////////////
//
// At the moment the same functionality/data-members as parent AliITShit 
// except the geometry transformation uses AliITSgeomTGeoUp 
//
////////////////////////////////////////////////////////////////////////

#include "TLorentzVector.h"

#include "AliITSMFTHit.h" 

class AliITSUGeomTGeo;

class AliITSUHit : public AliITSMFTHit {

public:
  AliITSUHit();
  AliITSUHit(Int_t shunt, Int_t track, Int_t *vol, Float_t *hits);
  AliITSUHit(Int_t, Int_t, Int_t *,Float_t,Float_t,
	     TLorentzVector &x,TLorentzVector &x0,TLorentzVector &p);
  AliITSUHit(const AliITSUHit &h);
  virtual ~AliITSUHit() {}

  void
  GetChipID(Int_t &lay,Int_t &stav,Int_t &sstav, Int_t &mod, Int_t &det, const AliITSUGeomTGeo *) const;

protected:
  AliITSUHit& operator=(const AliITSUHit &h);
  ClassDef(AliITSUHit,3)  //Hits object
	 
}; 

#endif
