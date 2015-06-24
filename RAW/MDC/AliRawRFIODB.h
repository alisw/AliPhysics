#ifndef ALIRAWRFIODB_H
#define ALIRAWRFIODB_H
// @(#) $Id$
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
   AliRawRFIODB(AliRawEventV2 *event,
		AliESDEvent *esd,
		Int_t compress,
		const char* fileName = NULL,
		Int_t basketsize = 32000,
		Long64_t autoflush = -5000000LL);
   ~AliRawRFIODB() { Close(); }

   Long64_t Close();

private:
   const char *GetFileName() const;

   ClassDef(AliRawRFIODB,0)  // Raw DB via RFIO
};

#endif
