#ifndef ALIRUNDB_H
#define ALIRUNDB_H
// @(#)alimdc:$Name$:$Id$
// Author: Fons Rademakers  26/11/99

/* Copyright(c) 1998-2003, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// AliRunDB                                                             //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#ifndef ROOT_TObject
#include <TObject.h>
#endif


// Forward class declarations
class AliStats;
class TFile;


class AliRunDB : public TObject {

public:
   AliRunDB(Bool_t noLocalDB = kFALSE);
   ~AliRunDB() { Close(); }

   void Update(AliStats *stats);
   void UpdateRDBMS(AliStats *stats);
   void UpdateAliEn(AliStats *stats);
   void Close();

private:
   TFile  *fRunDB;     // run database

   AliRunDB(const AliRunDB& runDB);
   AliRunDB& operator = (const AliRunDB& runDB);

   ClassDef(AliRunDB,0)  // Run (bookkeeping) DB
};

#endif
