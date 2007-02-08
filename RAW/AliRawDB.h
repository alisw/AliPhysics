#ifndef ALIRAWDB_H
#define ALIRAWDB_H
// @(#)alimdc:$Name$:$Id$
// Author: Fons Rademakers  26/11/99

/* Copyright(c) 1998-2003, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// AliRawDB                                                             //
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
class AliRawEvent;
class AliStats;
class TFile;
class AliESD;

class AliRawDB : public TObject {

public:
   AliRawDB(AliRawEvent *event,
	    AliESD *esd,
	    Int_t compress,
            const char* fileName = NULL);
   virtual ~AliRawDB() { if(Close()==-1) Error("~AliRawDB", "cannot close output file!"); }

   virtual const char *GetOpenOption() const { return "RECREATE"; }
   virtual Int_t       GetNetopt() const { return 0; }
   virtual Bool_t      Create(const char* fileName = NULL);
   virtual Int_t       Close();
   Int_t               Fill();
   Int_t               GetTotalSize();

   void         WriteStats(AliStats* stats);

   void         SetMaxSize(Double_t maxSize) { fMaxSize = maxSize; }
   void         SetFS(const char* fs1, const char* fs2 = NULL);
   void         SetDeleteFiles(Bool_t deleteFiles = kTRUE) { fDeleteFiles = deleteFiles; }
   void         SetGuidFileFolder(const char* guidFileFolder)
     { fGuidFileFolder = guidFileFolder; }

   Bool_t       NextFile(const char* fileName = NULL);

   Double_t     GetBytesWritten() const { return fRawDB->GetBytesWritten(); }
   TFile       *GetDB() const { return fRawDB; }
   const char  *GetDBName() const { return fRawDB->GetName(); }
   Int_t        GetEvents() const { return (Int_t) fTree->GetEntries(); }
   AliRawEvent *GetEvent() const { return fEvent; }
   Float_t      GetCompressionFactor() const;
   Int_t        GetCompressionMode() const { return fRawDB->GetCompressionLevel(); }
   void         Stop() { fStop = kTRUE; }
   static const char *GetAliRootTag();

protected:
   TFile         *fRawDB;         // DB to store raw data
   TTree         *fTree;          // tree used to store raw data
   AliRawEvent   *fEvent;         // AliRawEvent via which data is stored
   TTree         *fESDTree;       // tree for storing HLT ESD information
   AliESD        *fESD;           // pointer to HLT ESD object
   Int_t          fCompress;      // compression mode (1 default)
   Double_t       fMaxSize;       // maximum size in bytes of the raw DB
   TString        fFS1;           // first raw DB file system location
   TString        fFS2;           // second raw DB file system location
   Bool_t         fDeleteFiles;   // flag for deletion of files
   Bool_t         fStop;          // stop execution (triggered by SIGUSR1)
   const char *fGuidFileFolder;   // folder which contains guid files
   static const char  *fgkAliRootTag; // string with the aliroot tag id

   virtual const char *GetFileName() const;
   virtual Bool_t      FSHasSpace(const char *fs) const;
   virtual void        MakeTree();

   void WriteGuidFile();

private:
   AliRawDB(const AliRawDB& rawDB);
   AliRawDB& operator = (const AliRawDB& rawDB);

   ClassDef(AliRawDB,0)  // Raw DB
};

#endif
