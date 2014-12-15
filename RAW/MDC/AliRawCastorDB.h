#ifndef ALIRAWCASTORDB_H
#define ALIRAWCASTORDB_H
// @(#) $Id$
// Author: Fons Rademakers  26/11/99

/* Copyright(c) 1998-2003, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// AliRawCastorDB                                                       //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "AliRawDB.h"


class AliRawCastorDB : public AliRawDB {

public:
   AliRawCastorDB(AliRawEventV2 *event,
		  AliESDEvent *esd,
		  Int_t compress,
		  const char* fileName = NULL,
		  Int_t basketsize = 32000);
   ~AliRawCastorDB() { Close(); }

   const char *GetOpenOption() const { return "-RECREATE"; }
   Int_t       GetNetopt() const { return 0; }
   Long64_t    Close();

private:
   const char *GetFileName() const;

   ClassDef(AliRawCastorDB,0)  // Raw DB via CASTOR and rootd
};

#endif
