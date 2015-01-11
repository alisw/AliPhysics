#ifndef ALIRAWDB_H
#define ALIRAWDB_H
// @(#)alimdc:$Name:  $:$Id$
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

#include "AliDAQ.h"

// Forward class declarations
class AliRawEventV2;
class AliRawDataArrayV2;
class TFile;
class AliESDEvent;

class AliRawDB : public TObject {

public:
   AliRawDB(AliRawEventV2 *event,
	    AliESDEvent *esd,
	    Int_t compress,
            const char* fileName = NULL,
	    Int_t basketsize = 32000);
   virtual ~AliRawDB();

   virtual const char *GetOpenOption() const { return "RECREATE"; }
   virtual Int_t       GetNetopt() const { return 0; }
   virtual Bool_t      Create(const char* fileName = NULL);
   virtual Long64_t    Close();
   Int_t               Fill();
   Long64_t            GetTotalSize();
   Long64_t            AutoSave();

   void         SetMaxSize(Double_t maxSize) { fMaxSize = maxSize; }
   void         SetFS(const char* fs1, const char* fs2 = NULL);
   void         SetDeleteFiles(Bool_t deleteFiles = kTRUE) { fDeleteFiles = deleteFiles; }

   Bool_t       NextFile(const char* fileName = NULL);

   Double_t     GetBytesWritten() const { return fRawDB->GetBytesWritten(); }
   TFile       *GetDB() const { return fRawDB; }
   const char  *GetDBName() const { return fRawDB->GetName(); }
   Int_t        GetEvents() const { return (fTree) ? (Int_t) fTree->GetEntries() : 0; }
   AliRawEventV2 *GetEvent() const { return fEvent; }
   Float_t      GetCompressionFactor() const;
   Int_t        GetCompressionMode() const { return fRawDB->GetCompressionLevel(); }
   Int_t        GetBasketSize() const { return fBasketSize; }
   void         Stop() { fStop = kTRUE; }
   static const char *GetAliRootTag();
   Bool_t       WriteGuidFile(TString &guidFileFolder);
   void         Reset();
   AliRawDataArrayV2 *GetRawDataArray(UInt_t eqSize, UInt_t eqId) const;

protected:
   TFile         *fRawDB;         // DB to store raw data
   TTree         *fTree;          // tree used to store raw data
   AliRawEventV2   *fEvent;         // AliRawEvent via which data is stored
   AliRawDataArrayV2  **fDetRawData[AliDAQ::kNDetectors+1]; // Detectors raw-data payload
   TTree         *fESDTree;       // tree for storing HLT ESD information
   AliESDEvent        *fESD;           // pointer to HLT ESD object
   Int_t          fCompress;      // compression mode (1 default)
   Int_t          fBasketSize;    // root i/o basket size (default = 32000)
   Double_t       fMaxSize;       // maximum size in bytes of the raw DB
   TString        fFS1;           // first raw DB file system location
   TString        fFS2;           // second raw DB file system location
   Bool_t         fDeleteFiles;   // flag for deletion of files
   Bool_t         fStop;          // stop execution (triggered by SIGUSR1)
   static const char  *fgkAliRootTag; // string with the aliroot tag id

   static Int_t   fgkDetBranches[AliDAQ::kNDetectors+1]; // Number of branches in each detector

   virtual const char *GetFileName() const;
   virtual Bool_t      FSHasSpace(const char *fs) const;
   virtual void        MakeTree();

private:
   AliRawDB(const AliRawDB& rawDB);
   AliRawDB& operator = (const AliRawDB& rawDB);

   ClassDef(AliRawDB,5)  // Raw DB
};

#endif
