#ifndef ALITAGDB_H
#define ALITAGDB_H
// @(#)alimdc:$Name$:$Id$
// Author: Fons Rademakers  26/11/99

/* Copyright(c) 1998-2003, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// AliTagDB                                                             //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#ifndef ROOT_TObject
#include <TObject.h>
#endif

#ifndef ROOT_TFile
#include <TFile.h>
#endif

#ifndef ROOT_TTree
#include <TTree.h>
#endif


// Forward class declarations
class AliRawEventHeader;


class AliTagDB : public TObject {

public:
   AliTagDB(AliRawEventHeader *header, Double_t maxsize, Bool_t create = kTRUE);
   virtual ~AliTagDB() { Close(); }

   Bool_t          Create();
   virtual void    Close();
   void            Fill() { fTree->Fill(); }
   Bool_t          FileFull()
            { return (fTagDB->GetBytesWritten() > fMaxSize) ? kTRUE : kFALSE; }

   Bool_t          NextFile();

   Double_t           GetBytesWritten() const { return fTagDB->GetBytesWritten(); }
   TFile             *GetDB() const { return fTagDB; }
   const char        *GetDBName() const { return fTagDB->GetName(); }
   AliRawEventHeader *GetHeader() const { return fHeader; }
   Int_t              GetEvents() const { return (Int_t) fTree->GetEntries(); }
   Float_t            GetCompressionFactor() const;

protected:
   TFile             *fTagDB;     // DB to store header information only (tag)
   TTree             *fTree;      // tree use to store header
   AliRawEventHeader *fHeader;    // header via which data is stored
   Double_t           fMaxSize;   // maximum size in bytes of tag DB

   virtual const char *GetFileName() const;

private:
   AliTagDB(const AliTagDB& tagDB);
   AliTagDB& operator = (const AliTagDB& tagDB);

   ClassDef(AliTagDB,0)  // Tag DB
};

#endif
