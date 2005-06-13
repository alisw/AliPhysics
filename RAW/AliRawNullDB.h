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
   AliRawNullDB(AliRawEvent *event,
		AliESD *esd,
		Int_t compress,
		const char* fileName);
   ~AliRawNullDB() { Close(); }

   Int_t Close();

private:
   const char *GetFileName() const;

   ClassDef(AliRawNullDB,0)  // Raw DB to /dev/null
};

#endif
