#ifndef __CORRELRECOPARENT_H__
#define __CORRELRECOPARENT_H__
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */
/* $Id:  $ */

//______________________________________________________________________________
// Container class for reconstructed parents. Reconstruction method uses
// AliKFParticle or TLorentzVector as indicated by kUseAliKF set in CorrelDefs.h
//-- Author: Paul Constantin

#include "CorrelParticle.h"
#include "CorrelKFTrack.h"

namespace JetCorrelHD {
   
  class CorrelRecoParent_t : public CorrelParticle_t {      
  public:
    CorrelRecoParent_t();
    virtual ~CorrelRecoParent_t() {;}
    CorrelRecoParent_t* operator=(const CorrelRecoParent_t& rhs);
    virtual CorrelRecoParent_t* Copy();

    Float_t Assym()   const {return fAssym;}
    Float_t OpenAng() const {return fOpenAng;}
    AliVEvent* Evt()  const {return fEVT;}
    Bool_t Reconstruct(CorrelParticle_t* p1, CorrelParticle_t* p2);
    void SetEvent(AliVEvent * const v) {fEVT=v;}

    virtual void Show();
    
  private:
    Float_t fAssym;   // children energy assymetry
    Float_t fOpenAng; // children opening angle
    AliVEvent* fEVT;  // input event (ESD or AOD)

    // disable (make private) the copy constructor
    CorrelRecoParent_t(const CorrelRecoParent_t &p);

    Bool_t NotInMass(PartType_t ID, Float_t mass);
  };
  
} // namespace declaration

#endif
