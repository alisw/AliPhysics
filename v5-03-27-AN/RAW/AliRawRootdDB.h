#ifndef ALIRAWROOTDDB_H
#define ALIRAWROOTDDB_H
// @(#) $Id$
// Author: Fons Rademakers  26/11/99

/* Copyright(c) 1998-2003, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// AliRawRootdDB                                                        //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "AliRawDB.h"


class AliRawRootdDB : public AliRawDB {

public:
   AliRawRootdDB(AliRawEventV2 *event,
		 AliESDEvent *esd,
		 Int_t compress,
		 const char* fileName = NULL,
		 Int_t basketsize = 32000);
   ~AliRawRootdDB() { Close(); }

   Long64_t Close();

private:
   const char *GetFileName() const;

   ClassDef(AliRawRootdDB,0)  // Raw DB via rootd
};

#endif
