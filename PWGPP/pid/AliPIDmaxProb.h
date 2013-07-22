#ifndef ALIPIDMAXPROB_H
#define ALIPIDMAXPROB_H


/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: AliPIDmaxProb.h 49869 2013-07-10 04:49:51Z fnoferin $ */

/////////////////////////////////////////////////
//                                             //
//             PID cut max prob                //
//           noferini@bo.infn.it               //
/////////////////////////////////////////////////

#include"AliPIDperfCut.h"
#include"AliPIDCombined.h"

class AliPIDmaxProb : public AliPIDperfCut
{
 public:
  AliPIDmaxProb(const char *name);
  AliPIDmaxProb();

  void SetPIDMask(Int_t mask){fMaskPID=mask;}; // set the PID mask to be used

  Bool_t IsSelected(AliVTrack *track,AliPID::EParticleType type) const;

 private:
  AliPIDmaxProb(const AliPIDmaxProb &old);
  AliPIDmaxProb& operator=(const AliPIDmaxProb &source);

  AliPIDCombined *fPIDCombined;         //! PID combined object
  Int_t fMaskPID;                       // PID mask

  ClassDef(AliPIDmaxProb,1)             // PID cut virtual class
};
#endif
