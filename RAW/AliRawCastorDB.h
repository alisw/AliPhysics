#ifndef ALIRAWCASTORDB_H
#define ALIRAWCASTORDB_H
// @(#)alimdc:$Name$:$Id$
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
   AliRawCastorDB(AliRawEvent *event,
		  AliESD *esd,
		  Int_t compress,
		  const char* fileName = NULL);
   ~AliRawCastorDB() { Close(); }

   const char *GetOpenOption() const { return "-RECREATE"; }
   Int_t       GetNetopt() const { return 0; }
   void        Close();

private:
   const char *GetFileName() const;

   ClassDef(AliRawCastorDB,0)  // Raw DB via CASTOR and rootd
};

#endif
