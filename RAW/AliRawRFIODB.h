#ifndef ALIRAWRFIODB_H
#define ALIRAWRFIODB_H
// @(#)alimdc:$Name$:$Id$
// Author: Fons Rademakers  26/11/99

/* Copyright(c) 1998-2003, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// AliRawRFIODB                                                         //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "AliRawDB.h"


class AliRawRFIODB : public AliRawDB {

public:
   AliRawRFIODB(AliRawEvent *event, Double_t maxsize, Int_t compress);
   ~AliRawRFIODB() { Close(); }

   void Close();

private:
   const char *GetFileName() const;

   ClassDef(AliRawRFIODB,0)  // Raw DB via RFIO
};

#endif
