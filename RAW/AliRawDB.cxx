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

#ifdef ALI_DATE
#include "event.h"
#endif

#include "AliRawEvent.h"
#include "AliRawEventHeader.h"
#include "AliStats.h"

#include "AliRawDB.h"


ClassImp(AliRawDB)


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

   // Consistency check with DATE header file
#ifdef ALI_DATE
   if (fEvent->GetHeader()->HeaderSize() != EVENT_HEAD_BASE_SIZE) {
      Error("AliRawDB", "inconsistency between DATE and AliRawEvent headers");
      MakeZombie();
      return;
   }
#endif

   if (fileName) {
      if (!Create(fileName))
         MakeZombie();
   }
}

//______________________________________________________________________________
AliRawDB::AliRawDB(const AliRawDB& rawDB): TObject(rawDB)
{
// copy constructor

  Fatal("AliRawDB", "copy constructor not implemented");
}

//______________________________________________________________________________
AliRawDB& AliRawDB::operator = (const AliRawDB& /*rawDB*/)
{
// assignment operator

  Fatal("operator =", "assignment operator not implemented");
  return *this;
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

   const Int_t kMaxRetry = 200;
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
                        Form("ALICE MDC%d raw DB", kMDC), fCompress,
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

   fTree = new TTree("RAW", Form("ALICE MDC%d raw data tree", kMDC));
   fTree->SetAutoSave(2000000000);  // autosave when 2 Gbyte written

   Int_t bufsize = 256000;
   // splitting 29.6 MB/s, no splitting 35.3 MB/s on P4 2GHz 15k SCSI
   //Int_t split   = 1;
   Int_t split   = 0;
   fTree->Branch("rawevent", "AliRawEvent", &fEvent, bufsize, split);

   // Create tree which will contain the HLT ESD information

   if (fESD) {
     fESDTree = new TTree("esdTree", Form("ALICE MDC%d HLT ESD tree", kMDC));
     fESDTree->SetAutoSave(2000000000);  // autosave when 2 Gbyte written
     split   = 99;
     fESDTree->Branch("ESD", "AliESD", &fESD, bufsize, split);
   }

}

//______________________________________________________________________________
void AliRawDB::Close()
{
   // Close raw DB.

   if (!fRawDB) return;

   fRawDB->cd();

   // Write the tree.
   fTree->Write();
   if (fESDTree) fESDTree->Write();

   // Close DB, this also deletes the fTree
   fRawDB->Close();

   if (fDeleteFiles) {
      gSystem->Unlink(fRawDB->GetName());
      delete fRawDB;
      fRawDB = 0;
      return;
   }

   // Create semaphore to say this file is finished
   Int_t tfd = ::creat(Form("%s.done", fRawDB->GetName()), 0644);
   close(tfd);

   delete fRawDB;
   fRawDB = 0;
}

//______________________________________________________________________________
Int_t AliRawDB::Fill()
{
   // Fill the trees and return the number of written bytes

   Double_t bytes = fRawDB->GetBytesWritten();
   fTree->Fill();
   if (fESDTree) fESDTree->Fill();
   return Int_t(fRawDB->GetBytesWritten() - bytes);
}

//______________________________________________________________________________
void AliRawDB::WriteStats(AliStats* stats)
{
   // Write stats to raw DB, local run DB and global MySQL DB.

   AliRawEventHeader &header = *GetEvent()->GetHeader();

   // Write stats into RawDB
   TDirectory *ds = gDirectory;
   GetDB()->cd();
   stats->SetEvents(GetEvents());
   stats->SetLastId(header.GetRunNumber(), header.GetEventInRun());
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
