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


// Forward class declarations
class AliRawEvent;
class AliStats;
class TFile;

#ifdef USE_HLT
class AliESD;
#endif

class AliRawDB : public TObject {

public:
   AliRawDB(AliRawEvent *event,
#ifdef USE_HLT
	    AliESD *esd,
#endif
	    Double_t maxsize, Int_t compress,
            Bool_t create = kTRUE);
   virtual ~AliRawDB() { Close(); }

   virtual const char *GetOpenOption() const { return "RECREATE"; }
   virtual Int_t       GetNetopt() const { return 0; }
   virtual Bool_t      Create();
   virtual void        Close();
   void                Fill() { fTree->Fill();
#ifdef USE_HLT
                                fESDTree->Fill();
#endif
                              }
   Bool_t              FileFull() { return (fRawDB->GetBytesWritten() > fMaxSize) ?
                                    kTRUE : kFALSE; }

   void         WriteStats(AliStats* stats);

   Bool_t       NextFile();

   Double_t     GetBytesWritten() const { return fRawDB->GetBytesWritten(); }
   TFile       *GetDB() const { return fRawDB; }
   const char  *GetDBName() const { return fRawDB->GetName(); }
   Int_t        GetEvents() const { return (Int_t) fTree->GetEntries(); }
   AliRawEvent *GetEvent() const { return fEvent; }
   Float_t      GetCompressionFactor() const;
   Int_t        GetCompressionMode() const { return fRawDB->GetCompressionLevel(); }

protected:
   TFile         *fRawDB;         // DB to store raw data
   TTree         *fTree;          // tree used to store raw data
   AliRawEvent   *fEvent;         // AliRawEvent via which data is stored
#ifdef USE_HLT
   TTree         *fESDTree;       // tree for storing HLT ESD information
   AliESD        *fESD;           // pointer to HLT ESD object
#endif
   Int_t          fCompress;      // compression mode (1 default)
   Double_t       fMaxSize;       // maximum size in bytes of the raw DB

   virtual const char *GetFileName() const;
   virtual Bool_t      FSHasSpace(const char *fs) const;
   virtual void        MakeTree();

private:
   AliRawDB(const AliRawDB& rawDB);
   AliRawDB& operator = (const AliRawDB& rawDB);

   ClassDef(AliRawDB,0)  // Raw DB
};

#endif
