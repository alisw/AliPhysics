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

#ifndef ROOT_TString
#include <TString.h>
#endif



// Forward class declarations
class AliRawEventTag;


class AliTagDB : public TObject {

public:
   AliTagDB(AliRawEventTag *eventTag, const char* fileName = NULL);
   virtual ~AliTagDB() { Close(); }

   Bool_t          Create(const char* fileName = NULL);
   virtual void    Close();
   void            Fill() { fTree->Fill(); }
   Bool_t          FileFull()
            { return (fMaxSize >= 0) ? ((fTagDB->GetBytesWritten() > fMaxSize) ? kTRUE : kFALSE) : kFALSE; }

   Bool_t          NextFile(const char* fileName = NULL);

   void            SetMaxSize(Double_t maxSize) { fMaxSize = maxSize; }
   void            SetFS(const char* fs);

   Double_t           GetBytesWritten() const { return fTagDB->GetBytesWritten(); }
   TFile             *GetDB() const { return fTagDB; }
   const char        *GetDBName() const { return fTagDB->GetName(); }
   AliRawEventTag    *GetEventTag() const { return fEventTag; }
   Int_t              GetEvents() const { return (Int_t) fTree->GetEntries(); }
   Float_t            GetCompressionFactor() const;

protected:
   TFile             *fTagDB;     // DB to store header information only (tag)
   TTree             *fTree;      // tree use to store header
   AliRawEventTag    *fEventTag;  // pointer to event tag object via which data is stored
   Double_t           fMaxSize;   // maximum size in bytes of tag DB
   TString            fFS;        // tag DB file system location
   Bool_t             fDeleteFiles; // flag for deletion of files

   virtual const char *GetFileName() const;

private:
   AliTagDB(const AliTagDB& tagDB);
   AliTagDB& operator = (const AliTagDB& tagDB);

   ClassDef(AliTagDB,0)  // Tag DB
};

#endif
