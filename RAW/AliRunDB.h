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

#ifndef ROOT_TString
#include <TString.h>
#endif


// Forward class declarations
class AliStats;
class TFile;


class AliRunDB : public TObject {

public:
   AliRunDB(const char* localFS, Bool_t rdbms = kFALSE, 
	    const char* alienHost = NULL, const char* alienDir = NULL);
   ~AliRunDB() { Close(); }

   void Update(AliStats *stats);
   void UpdateLocal(AliStats *stats);
   void UpdateRDBMS(AliStats *stats);
   void UpdateAliEn(AliStats *stats);
   void Close();

private:
   TFile  *fRunDB;     // run database
   Bool_t  fRDBMS;     // flag for usage of central MySQL DB
   TString fAlienHost; // alien host name
   TString fAlienDir;  // alien directory

   AliRunDB(const AliRunDB& runDB);
   AliRunDB& operator = (const AliRunDB& runDB);

   ClassDef(AliRunDB,0)  // Run (bookkeeping) DB
};

#endif
