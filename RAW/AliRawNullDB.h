#ifndef ALIRAWNULLDB_H
#define ALIRAWNULLDB_H
// @(#)alimdc:$Name$:$Id$
// Author: Fons Rademakers  26/11/99

/* Copyright(c) 1998-2003, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// AliRawNullDB                                                         //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "AliRawDB.h"


class AliRawNullDB : public AliRawDB {

public:
   AliRawNullDB(AliRawEvent *event, Double_t maxsize, Int_t compress);
   ~AliRawNullDB() { Close(); }

   void Close();

private:
   const char *GetFileName() const;

   ClassDef(AliRawNullDB,0)  // Raw DB to /dev/null
};

#endif
