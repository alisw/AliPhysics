/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/


/////////////////////////////////////////////////////
// Class to handle pairs of tracks of opposite charge
// Useful for resonance analysis
// Derives from AliVParticle => 
// usable in Correction Framework
/////////////////////////////////////////////////////
// author : renaud.vernet@cern.ch
/////////////////////////////////////////////////////


#ifndef ALICFPAIR_H
#define ALICFPAIR_H

#include "AliVParticle.h"

class AliESDtrack ;
class AliESDv0;
class AliESDEvent;
class AliAODv0;

class AliCFPair : public AliVParticle {

 public:
  AliCFPair(AliVParticle* t1, AliVParticle* t2);
  AliCFPair(AliESDv0* v0, AliESDEvent* esd);
  AliCFPair(AliAODv0* v0);
  AliCFPair(const AliCFPair& c);
  AliCFPair& operator=(const AliCFPair& c);
  virtual ~AliCFPair(){};

  AliVParticle* GetNeg() const {return fTrackNeg;}
  AliVParticle* GetPos() const {return fTrackPos;}
  AliESDv0*    GetESDV0()  const {return fESDV0;}
  AliAODv0*    GetAODV0()  const {return fAODV0;}
  void         SetV0PDG(Int_t pdg) {fV0PDG=pdg;}
  virtual Bool_t       PxPyPz(Double_t p[3]) const ;
  virtual Double32_t   P()  const ;
  virtual Double32_t   Pt() const ;
  virtual Double32_t   Px() const ;
  virtual Double32_t   Py() const ;
  virtual Double32_t   Pz() const ;
  virtual Double32_t   E () const ;
  virtual Double32_t   Xv() const ;
  virtual Double32_t   Yv() const ;
  virtual Double32_t   Zv() const ;
  virtual Bool_t       XvYvZv(Double_t x[3]) const ;

  virtual Double32_t OneOverPt() const {return 1/Pt();}
  virtual Double32_t Phi()   const ;
  virtual Double32_t Theta() const ;
  virtual Double32_t M() const ;
  virtual Double32_t Eta() const ;
  virtual Double32_t Y() const ;
  virtual Short_t    Charge() const {return 0;} // returns 0 because opposite charge tracks... maybe to extend to all kinds of pairs
  virtual Int_t      GetLabel() const {return fLabel;}
  virtual void       SetLabel(Int_t label) {fLabel=label;}
  // PID
  virtual const Double_t *PID() const {return 0;} // return PID object (to be defined, still)


 private:
  Bool_t fIsV0;            // true if V0 passed to the constructor
  AliVParticle* fTrackNeg; // pointer to the negative track 
  AliVParticle* fTrackPos; // pointer to the positive track 
  AliESDv0*    fESDV0;     // pointer to the ESD V0 if AliESDv0 is passed to the constructor
  AliAODv0*    fAODV0;     // pointer to the AOD V0 if AliAODv0 is passed to the constructor
  Int_t        fLabel;     // associated MC label
  Int_t        fV0PDG;     // assumed V0 PDG
  
  ClassDef(AliCFPair,0);
};

#endif
