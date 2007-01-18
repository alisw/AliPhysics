#ifndef ALITAGNULLDB_H
#define ALITAGNULLDB_H
// @(#)alimdc:$Name$:$Id$
// Author: Fons Rademakers  26/11/99

/* Copyright(c) 1998-2003, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// AliTagNullDB                                                         //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "AliTagDB.h"


class AliTagNullDB : public AliTagDB {

public:
   AliTagNullDB(AliRawEventTag *eventTag);
   ~AliTagNullDB() { Close(); }

   void Close();

private:
   const char *GetFileName() const;

   ClassDef(AliTagNullDB,0)   // Tag DB to /dev/null
};

#endif
