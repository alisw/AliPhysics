#ifndef ALIRAWROOTDDB_H
#define ALIRAWROOTDDB_H
// @(#)alimdc:$Name$:$Id$
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
   AliRawRootdDB(AliRawEvent *event,
#ifdef USE_HLT
		 AliESD *esd,
#endif
		 Double_t maxsize, Int_t compress);
   ~AliRawRootdDB() { Close(); }

   void Close();

private:
   const char *GetFileName() const;

   ClassDef(AliRawRootdDB,0)  // Raw DB via rootd
};

#endif
