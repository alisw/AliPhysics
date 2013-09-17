#ifndef ALIPIDPERFCUT_H
#define ALIPIDPERFCUT_H


/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: AliPIDperfCut.h 49869 2013-07-10 04:49:51Z fnoferin $ */

/////////////////////////////////////////////////
//                                             //
//         PID cut for performanve             //
//           noferini@bo.infn.it               //
/////////////////////////////////////////////////

#include"AliVTrack.h"
#include"AliPID.h"
#include"TNamed.h"

class AliPIDperfCut : public TNamed
{
 public:
  AliPIDperfCut(const char *name);
  AliPIDperfCut();

  virtual Bool_t IsSelected(AliVTrack *track,AliPID::EParticleType type) const;

 private:

  ClassDef(AliPIDperfCut,1)  // PID cut virtual class
};
#endif
