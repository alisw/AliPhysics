// @(#)alimdc:$Name$:$Id$
// Author: Fons Rademakers  26/11/99

/**************************************************************************
 * Copyright(c) 1998-2003, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// AliRawDB                                                             //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include <errno.h>

#include <TSystem.h>
#include <TKey.h>

#include "AliESD.h"
#include "AliRawEvent.h"
#include "AliRawEventHeaderBase.h"
#include "AliStats.h"

#include "AliRawDB.h"


ClassImp(AliRawDB)

const char *AliRawDB::fgkAliRootTag = "$Name$";

//______________________________________________________________________________
AliRawDB::AliRawDB(AliRawEvent *event,
		   AliESD *esd, 
		   Int_t compress,
                   const char* fileName) :
  fRawDB(NULL),
  fTree(NULL),
  fEvent(event),
  fESDTree(NULL),
  fESD(esd),
  fCompress(compress),
  fMaxSize(-1),
  fFS1(""),
  fFS2(""),
  fDeleteFiles(kFALSE),
  fStop(kFALSE)
{
   // Create a new raw DB

   if (fileName) {
      if (!Create(fileName))
         MakeZombie();
   }
}

//______________________________________________________________________________
Bool_t AliRawDB::FSHasSpace(const char *fs) const
{
   // Check for at least fMaxSize bytes of free space on the file system.
   // If the space is not available return kFALSE, kTRUE otherwise.

   Long_t id, bsize, blocks, bfree;

   if (gSystem->GetFsInfo(fs, &id, &bsize, &blocks, &bfree) == 1) {
      Error("FSHasSpace", "could not stat file system %s", fs);
      return kFALSE;
   }

   // Leave 5 percent of diskspace free
   Double_t avail = Double_t(bfree) * 0.95;
   if (avail*bsize > fMaxSize)
      return kTRUE;

   Warning("FSHasSpace", "no space on file system %s", fs);
   return kFALSE;
}

//______________________________________________________________________________
const char *AliRawDB::GetFileName() const
{
   // Return filename based on hostname and date and time. This will make
   // each file unique. Also makes sure (via FSHasSpace()) that there is
   // enough space on the file system to store the file. Returns 0 in
   // case of error or interrupt signal.

   static TString fname;
   static Bool_t  fstoggle = kFALSE;

   TString fs = fstoggle ? fFS2 : fFS1;
   TDatime dt;

   TString hostname = gSystem->HostName();
   Int_t pos;
   if ((pos = hostname.Index(".")) != kNPOS)
      hostname.Remove(pos);

   if (!FSHasSpace(fs)) {
      while (1) {
         fstoggle = !fstoggle;
         fs = fstoggle ? fFS2 : fFS1;
         if (FSHasSpace(fs)) break;
         Info("GetFileName", "sleeping 30 seconds before retrying...");
         gSystem->Sleep(30000);   // sleep for 30 seconds
         if (fStop) return 0;
      }
   }

   fname = fs + "/" + hostname + "_";
   fname += dt.GetDate();
   fname += "_";
   fname += dt.GetTime();
   fname += ".root";

   fstoggle = !fstoggle;

   return fname;
}

//______________________________________________________________________________
void AliRawDB::SetFS(const char* fs1, const char* fs2)
{
// set the file system location

  fFS1 = fs1;
  if (fs1 && !fFS1.Contains(":")) {
    gSystem->ResetErrno();
    gSystem->MakeDirectory(fs1);
    if (gSystem->GetErrno() && gSystem->GetErrno() != EEXIST) {
      SysError("SetFS", "mkdir %s", fs1);
    }
  }

  fFS2 = fs2;
  if (fs2) {
    gSystem->ResetErrno();
    gSystem->MakeDirectory(fs2);
    if (gSystem->GetErrno() && gSystem->GetErrno() != EEXIST) {
      SysError("SetFS", "mkdir %s", fs2);
    }
  }
}

//______________________________________________________________________________
Bool_t AliRawDB::Create(const char* fileName)
{
   // Create a new raw DB.

   const Int_t kMaxRetry = 1;
   const Int_t kMaxSleep = 1;      // seconds
   const Int_t kMaxSleepLong = 10; // seconds
   Int_t retry = 0;

again:
   if (fStop) return kFALSE;

   const char *fname = fileName;
   if (!fname) fname = GetFileName();
   if (!fname) {
      Error("Create", "error getting raw DB file name");
      return kFALSE;
   }

   retry++;

   fRawDB = TFile::Open(fname, GetOpenOption(),
			Form("ALICE raw-data file (%s)", GetAliRootTag()), fCompress,
			GetNetopt());
   if (!fRawDB) {
      if (retry < kMaxRetry) {
         Warning("Create", "failure to open file, sleeping %d %s before retrying...",
                 kMaxSleep, kMaxSleep==1 ? "second" : "seconds");
         gSystem->Sleep(kMaxSleep*1000);
         goto again;
      }
      Error("Create", "failure to open file %s after %d tries", fname, kMaxRetry);
      return kFALSE;
   }
   if (retry > 1)
      Warning("Create", "succeeded to open file after %d retries", retry);

   if (fRawDB->IsZombie()) {
      if (fRawDB->GetErrno() == ENOSPC ||
          fRawDB->GetErrno() == 1018   ||   // SECOMERR
          fRawDB->GetErrno() == 1027) {     // SESYSERR
         fRawDB->ResetErrno();
         delete fRawDB;
         Warning("Create", "file is a zombie (no space), sleeping %d %s before retrying...",
                 kMaxSleepLong, kMaxSleepLong==1 ? "second" : "seconds");
         gSystem->Sleep(kMaxSleepLong*1000);   // sleep 10 seconds before retrying
         goto again;
      }
      Error("Create", "file %s is zombie", fname);
      fRawDB->ResetErrno();
      delete fRawDB;
      fRawDB = 0;
      if (retry < kMaxRetry) {
         Warning("Create", "file is a zombie, sleeping %d %s before retrying...",
                 kMaxSleep, kMaxSleep==1 ? "second" : "seconds");
         gSystem->Sleep(kMaxSleep*1000);
         goto again;
      }
      Error("Create", "failure to open file %s after %d tries", fname, kMaxRetry);
      return kFALSE;
   }

   // Create raw data TTree
   MakeTree();

   return kTRUE;
}

//______________________________________________________________________________
void AliRawDB::MakeTree()
{
   // Create ROOT Tree object container.

   fTree = new TTree("RAW", Form("ALICE raw-data tree (%s)", GetAliRootTag()));
   fTree->SetAutoSave(2000000000);  // autosave when 2 Gbyte written

   Int_t bufsize = 256000;
   // splitting 29.6 MB/s, no splitting 35.3 MB/s on P4 2GHz 15k SCSI
   //Int_t split   = 1;
   Int_t split   = 0;
   fTree->Branch("rawevent", "AliRawEvent", &fEvent, bufsize, split);

   // Create tree which will contain the HLT ESD information

   if (fESD) {
     fESDTree = new TTree("esdTree", Form("ALICE HLT ESD tree (%s)", GetAliRootTag()));
     fESDTree->SetAutoSave(2000000000);  // autosave when 2 Gbyte written
     split   = 0;
     fESDTree->Branch("ESD", "AliESD", &fESD, bufsize, split);
   }

}

//______________________________________________________________________________
Int_t AliRawDB::Close()
{
   // Close raw DB.
   if (!fRawDB) return 0;

   if (!fRawDB->IsOpen()) return 0;

   fRawDB->cd();

   // Write the tree.
   Bool_t error = kFALSE;
   if (fTree->Write() == 0)
     error = kTRUE;
   if (fESDTree)
     if (fESDTree->Write() == 0)
       error = kTRUE;

   // Close DB, this also deletes the fTree
   fRawDB->Close();

   Int_t filesize = fRawDB->GetEND();

   if (fDeleteFiles) {
      gSystem->Unlink(fRawDB->GetName());
      delete fRawDB;
      fRawDB = 0;
      if(!error)
	return filesize;
      else
	return -1;
   }

   delete fRawDB;
   fRawDB = 0;
   if(!error)
     return filesize;
   else
     return -1;
}

//______________________________________________________________________________
Int_t AliRawDB::Fill()
{
   // Fill the trees and return the number of written bytes

   Double_t bytes = fRawDB->GetBytesWritten();
   Bool_t error = kFALSE;
   if (fTree->Fill() == -1)
     error = kTRUE;
   if (fESDTree) 
     if (fESDTree->Fill() == -1)
       error = kTRUE;
   if(!error)
     return Int_t(fRawDB->GetBytesWritten() - bytes);
   else
     return -1;
}

//______________________________________________________________________________
Int_t AliRawDB::GetTotalSize()
{
   // Return the total size of the trees
  Int_t total = 0;

  {
    Int_t skey = 0;
    TDirectory *dir = fTree->GetDirectory();
    if (dir) {
      TKey *key = dir->GetKey(fTree->GetName());
      if (key) skey = key->GetKeylen();
    }
    total += skey;
    if (fTree->GetZipBytes() > 0) total += fTree->GetTotBytes();
    TBuffer b(TBuffer::kWrite,10000);
    TTree::Class()->WriteBuffer(b,fTree);
    total += b.Length();
  }

  if(fESDTree)
    {
      Int_t skey = 0;
      TDirectory *dir = fESDTree->GetDirectory();
      if (dir) {
	TKey *key = dir->GetKey(fESDTree->GetName());
	if (key) skey = key->GetKeylen();
      }
      total += skey;
      if (fESDTree->GetZipBytes() > 0) total += fESDTree->GetTotBytes();
      TBuffer b(TBuffer::kWrite,10000);
      TTree::Class()->WriteBuffer(b,fESDTree);
      total += b.Length();
    }

  return total;
}

//______________________________________________________________________________
void AliRawDB::WriteStats(AliStats* stats)
{
   // Write stats to raw DB, local run DB and global MySQL DB.

   AliRawEventHeaderBase &header = *GetEvent()->GetHeader();

   // Write stats into RawDB
   TDirectory *ds = gDirectory;
   GetDB()->cd();
   stats->SetEvents(GetEvents());
   stats->SetLastId(header.GetP("Id")[0]);
   stats->SetFileSize(GetBytesWritten());
   stats->SetCompressionFactor(GetCompressionFactor());
   stats->SetEndTime();
   stats->Write("stats");
   ds->cd();
}

//______________________________________________________________________________
Bool_t AliRawDB::NextFile(const char* fileName)
{
   // Close te current file and open a new one.
   // Returns kFALSE in case opening failed.

   Close();

   if (!Create(fileName)) return kFALSE;
   return kTRUE;
}

//______________________________________________________________________________
Float_t AliRawDB::GetCompressionFactor() const
{
   // Return compression factor.

   if (fTree->GetZipBytes() == 0.)
      return 1.0;
   else
      return fTree->GetTotBytes()/fTree->GetZipBytes();
}

//______________________________________________________________________________
const char *AliRawDB::GetAliRootTag()
{
  // Return the aliroot tag (version)
  // used to generate the raw data file.
  // Stored in the raw-data file title.

  TString version = fgkAliRootTag;
  version.Remove(TString::kBoth,'$');
  version.ReplaceAll("Name","AliRoot version");

  return version.Data();
}
