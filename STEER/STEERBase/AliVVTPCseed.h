#ifndef ALIVVTPCSEED_H
#define ALIVVTPCSEED_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               *
 * Primary Author: Mikolaj Krzewicki, mkrzewic@cern.ch
 */
#include "Rtypes.h"
class AliTPCseed;

class AliVVTPCseed {
  public:
  AliVVTPCseed() {}
  virtual ~AliVVTPCseed() {}
  virtual void CopyToTPCseed( AliTPCseed &) const = 0;
  ClassDef(AliVVTPCseed, 1);
};

#endif
