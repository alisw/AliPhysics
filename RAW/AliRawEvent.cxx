// @(#)alimdc:$Name$:$Id$
// Author: Fons Rademakers  26/11/99
// Updated: Dario Favretto  15/04/2003

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
// AliRawEvent                                                          //
//                                                                      //
// Set of classes defining the ALICE RAW event format. The AliRawEvent  //
// class defines a RAW event. It consists of an AliEventHeader object   //
// an AliEquipmentHeader object, an AliRawData object and an array of   //
// sub-events, themselves also being AliRawEvents. The number of        //
// sub-events depends on the number of DATE LDC's.                      //
// The AliRawEvent objects are written to a ROOT file using different   //
// technologies, i.e. to local disk via AliRawDB or via rfiod using     //
// AliRawRFIODB or via rootd using AliRawRootdDB or to CASTOR via       //
// rootd using AliRawCastorDB (and for performance testing there is     //
// also AliRawNullDB).                                                  //
// The AliRunDB class provides the interface to the run and file        //
// catalogues (AliEn or plain MySQL).                                   //
// The AliStats class provides statics information that is added as     //
// a single keyed object to each raw file.                              //
// The AliTagDB provides an interface to a TAG database.                //
// The AliMDC class is usid by the "alimdc" stand-alone program         //
// that reads data directly from DATE.                                  //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include <sys/types.h>
#include <sys/stat.h>
#if defined(__DECCXX)
#include <sys/statvfs.h>
#else
#if defined(__APPLE__)
#include <sys/param.h>
#include <sys/mount.h>
#else
#include <sys/vfs.h>
#endif
#endif
#include <unistd.h>
#include <stdlib.h>
#include <fcntl.h>
#include <errno.h>

#include <TBuffer.h>
#include <TSystem.h>
#include <TError.h>
#include <TStopwatch.h>
#include <TSQLServer.h>
#include <TSQLResult.h>
#include <TUrl.h>
#include <TGrid.h>
#include <TH1.h>

#if defined(__APPLE__)
#undef Free
#endif
#include "AliRawEvent.h"

// Good for Linux
#define long32 int
#include "DateEvent.h"

#ifdef USE_EB
#include "libDateEb.h"
#endif


ClassImp(AliRawEvent)
ClassImp(AliRawEventHeader)
ClassImp(AliRawEquipmentHeader)
ClassImp(AliRawData)
ClassImp(AliStats)
ClassImp(AliRawDB)
ClassImp(AliRawRFIODB)
ClassImp(AliRawCastorDB)
ClassImp(AliRawRootdDB)
ClassImp(AliRawNullDB)
ClassImp(AliTagDB)
ClassImp(AliTagNullDB)
ClassImp(AliRunDB)
ClassImp(AliMDC)

// Which MDC is this...
const Int_t kMDC = 5;

// Fixed file system locations for the different DB's
#ifdef USE_RDM
const char* const kFifo       = "/tmp/alimdc.fifo";
const char* const kRawDBFS[2] = { "/tmp/mdc1", "/tmp/mdc2" };
const char* const kTagDBFS    = "/tmp/mdc1/tags";
const char* const kRunDBFS    = "/tmp/mdc1/meta";
const char* const kRFIOFS     = "rfio:/castor/cern.ch/user/r/rdm";
const char* const kCastorFS   = "castor:/castor/cern.ch/user/r/rdm";
const char* const kRootdFS    = "root://localhost//tmp/mdc1";
const char* const kAlienHost  = "alien://aliens7.cern.ch:15000/?direct";
const char* const kAlienDir   = "/alice_mdc/DC";
#else
const char* const kFifo       = "/tmp/alimdc.fifo";
const char* const kRawDBFS[2] = { "/data1/mdc", "/data2/mdc" };
const char* const kTagDBFS    = "/data1/mdc/tags";
const char* const kRunDBFS    = "/data1/mdc/meta";
const char* const kRFIOFS     = "rfio:/castor/cern.ch/lcg/dc5";
const char* const kCastorFS   = "castor:/castor/cern.ch/lcg/dc5";
const char* const kRootdFS    = "root://localhost//tmp/mdc1";
const char* const kAlienHost  = "alien://aliens7.cern.ch:15000/?direct";
const char* const kAlienDir   = "/alice_mdc/DC";
#endif

// Maximum size of tag db files
const Double_t kMaxTagFileSize = 2.5e8;    // 250MB

Bool_t AliMDC::fgDeleteFiles = kFALSE;

AliMDC *gAliMDC = 0;


//______________________________________________________________________________
Bool_t AliRawEventHeader::DataIsSwapped() const
{
   // Returns true if event data is swapped.

   if (TEST_SYSTEM_ATTRIBUTE(fTypeAttribute, ATTR_EVENT_SWAPPED))
      return kTRUE;
   return kFALSE;
}

//______________________________________________________________________________
void AliRawEventHeader::Swap()
{
   // Swap header data.

   if (IsSwapped()) {
      fSize         = net2host(fSize);
      fMagic        = net2host(fMagic);
      fHeadLen      = net2host(fHeadLen);
      fVersion      = net2host(fVersion);
      fType         = net2host(fType);
      fRunNb        = net2host(fRunNb);
      for (int i = 0; i < kIdWords; i++)
         fId[i] = net2host(fId[i]);
      for (int i = 0; i < kTriggerWords; i++)
         fTriggerPattern[i] = net2host(fTriggerPattern[i]);
      for (int i = 0; i < kDetectorWords; i++)
         fDetectorPattern[i] = net2host(fDetectorPattern[i]);
      for (int i = 0; i < kAttributeWords; i++)
         fTypeAttribute[i] = net2host(fTypeAttribute[i]);
      fLDCId        = net2host(fLDCId);
      fGDCId        = net2host(fGDCId);
   }
}

//______________________________________________________________________________
UInt_t AliRawEventHeader::GetEventInRun() const
{
   // Get event number in run. Correct for fixed target mode which is used
   // in the Data Challenge Setup.

   if (!TEST_SYSTEM_ATTRIBUTE(fTypeAttribute, ATTR_ORBIT_BC)) {
      return EVENT_ID_GET_NB_IN_RUN(fId);
   }
   return 0;
}

//______________________________________________________________________________
const char *AliRawEventHeader::GetTypeName() const
{
   // Get event type as a string.

   switch (GetType()) {
      case kStartOfRun:
         return "START_OF_RUN";
         break;
      case kEndOfRun:
         return "END_OF_RUN";
         break;
      case kStartOfRunFiles:
         return "START_OF_RUN_FILES";
         break;
      case kEndOfRunFiles:
         return "END_OF_RUN_FILES";
         break;
      case kStartOfBurst:
         return "START_OF_BURST";
         break;
      case kEndOfBurst:
         return "END_OF_BURST";
         break;
      case kPhysicsEvent:
         return "PHYSICS_EVENT";
         break;
      case kCalibrationEvent:
         return "CALIBRATION_EVENT";
         break;
      case kFormatError:
         return "EVENT_FORMAT_ERROR";
         break;
      default:
         return "*** UNKNOWN EVENT TYPE ***";
         break;
   }
}


//______________________________________________________________________________
void AliRawEquipmentHeader::Swap()
{
   // Swap equipment header data. There is no way to see if the data
   // has already been swapped. This method is only called when the
   // header is read from the DATE event builder (GDC).

   fSize                 = net2host(fSize);
   fEquipmentType        = net2host(fEquipmentType);
   fEquipmentID          = net2host(fEquipmentID);
   fBasicElementSizeType = net2host(fBasicElementSizeType);
   for (int i = 0; i < kAttributeWords; i++)
      fTypeAttribute[i] = net2host(fTypeAttribute[i]);
}


//______________________________________________________________________________
AliRawEvent::AliRawEvent()
{
   // Create ALICE event object. If ownData is kFALSE we will use a static
   // raw data object, otherwise a private copy will be made.

   fNSubEvents = 0;
   fEvtHdr     = 0;
   fEqpHdr     = 0;
   fRawData    = 0;
   fSubEvents  = 0;
}

//______________________________________________________________________________
AliRawEventHeader *AliRawEvent::GetHeader()
{
   // Get event header part of AliRawEvent.

   if (!fEvtHdr)
      fEvtHdr = new AliRawEventHeader;

   return fEvtHdr;
}

//______________________________________________________________________________
AliRawEquipmentHeader *AliRawEvent::GetEquipmentHeader()
{
   // Get equipment header part of AliRawEvent.

   if (!fEqpHdr)
      fEqpHdr = new AliRawEquipmentHeader;

   return fEqpHdr;
}

//______________________________________________________________________________
AliRawData *AliRawEvent::GetRawData()
{
   // Get raw data part of AliRawEvent.

   if (!fRawData)
      fRawData = new AliRawData;

   return fRawData;
}

//______________________________________________________________________________
AliRawEvent *AliRawEvent::NextSubEvent()
{
   // Returns next sub-event object.

   if (!fSubEvents)
      fSubEvents = new TObjArray(100); // arbitrary, probably enough to prevent resizing

   if (fSubEvents->GetSize() <= fNSubEvents) {
      fSubEvents->Expand(fNSubEvents+10);
      Warning("NextSubEvent", "expanded fSubEvents by 10 to %d",
              fSubEvents->GetSize());
   }

   AliRawEvent *ev;
   if (!(ev = (AliRawEvent *)fSubEvents->At(fNSubEvents))) {
      ev = new AliRawEvent;
      fSubEvents->AddAt(ev, fNSubEvents);
   }

   fNSubEvents++;

   return ev;
}

//______________________________________________________________________________
AliRawEvent *AliRawEvent::GetSubEvent(Int_t index) const
{
   // Get specified sub event. Returns 0 if sub event does not exist.

   if (!fSubEvents)
      return 0;

   return (AliRawEvent *) fSubEvents->At(index);
}

//______________________________________________________________________________
void AliRawEvent::Reset()
{
   // Reset the event in case it needs to be re-used (avoiding costly
   // new/delete cycle). We reset the size marker for the AliRawData
   // objects and the sub event counter.

   for (int i = 0; i < fNSubEvents; i++) {
      AliRawEvent *ev = (AliRawEvent *)fSubEvents->At(i);
      ev->GetRawData()->SetSize(0);
   }
   fNSubEvents = 0;
}

//______________________________________________________________________________
AliRawEvent::~AliRawEvent()
{
   // Clean up event object. Delete also, possible, private raw data.

   delete fEvtHdr;
   delete fEqpHdr;
   delete fRawData;
   if (fSubEvents)
      fSubEvents->Delete();
   delete fSubEvents;
}

//______________________________________________________________________________
AliStats::AliStats(const char *filename, Int_t compmode, Bool_t filter)
{
   // Create statistics object.

   fEvents     = 0;
   fFirstRun   = 0;
   fFirstEvent = 0;
   fLastRun    = 0;
   fLastEvent  = 0;
   fChunk      = -0.5;
   fFileName   = filename;
   fCompMode   = compmode;
   fFilter     = filter;
   fRTHist     = 0;
}

//______________________________________________________________________________
AliStats::~AliStats()
{
   // Cleanup stats object.

   delete fRTHist;
}

//______________________________________________________________________________
AliStats &AliStats::operator=(const AliStats &rhs)
{
   // AliStats assignment operator.

   if (this != &rhs) {
      TObject::operator=(rhs);
      fEvents     = rhs.fEvents;
      fFirstRun   = rhs.fFirstRun;
      fFirstEvent = rhs.fFirstEvent;
      fLastRun    = rhs.fLastRun;
      fLastEvent  = rhs.fLastEvent;
      fBegin      = rhs.fBegin;
      fEnd        = rhs.fEnd;
      fFileName   = rhs.fFileName;
      fFileSize   = rhs.fFileSize;
      fCompFactor = rhs.fCompFactor;
      fCompMode   = rhs.fCompMode;
      fFilter     = rhs.fFilter;
      fRTHist     = rhs.fRTHist ? (TH1F*) rhs.fRTHist->Clone() : 0;
      fChunk      = rhs.fChunk;
   }
   return *this;
}

//______________________________________________________________________________
void AliStats::Fill(Float_t time)
{
   // Fill histogram. This histogram shows the (hopefully constant) time
   // it takes to fill the ROOT DB.
   // Expects to be called 100 times for each file.

   if (!fRTHist) {
      fRTHist = new TH1F("rtime","Real-time to write data chunk", 100, 0, 100);
      fRTHist->SetDirectory(0);
   }

   fRTHist->Fill(fChunk, time);
   fChunk += 1.0;
}

//______________________________________________________________________________
void AliStats::WriteToDB(AliRawDB *rawdb)
{
   // Write stats to raw DB, local run DB and global MySQL DB.

   AliRawEventHeader &header = *rawdb->GetEvent()->GetHeader();

   // Write stats into RawDB
   TDirectory *ds = gDirectory;
   rawdb->GetDB()->cd();
   SetEvents(rawdb->GetEvents());
   SetLastId(header.GetRunNumber(), header.GetEventInRun());
   SetFileSize(rawdb->GetBytesWritten());
   SetCompressionFactor(rawdb->GetCompressionFactor());
   SetEndTime();
   Write("stats");
   ds->cd();

   // Write stats also in the bookkeeping RunDB
   AliRunDB *rundb = new AliRunDB(kTRUE);
   rundb->Update(this);
   rundb->UpdateRDBMS(this);
   rundb->UpdateAliEn(this);
   delete rundb;
}

//______________________________________________________________________________
AliRawDB::AliRawDB(AliRawEvent *event, Double_t maxsize, Int_t compress,
                   Bool_t create)
{
   // Create a new raw DB containing at most maxsize bytes.

   fEvent    = event;
   fMaxSize  = maxsize;
   fCompress = compress;

   // Consistency check with DATE header file
   if (fEvent->GetHeader()->HeaderSize() != EVENT_HEAD_BASE_SIZE) {
      Error("AliRawDB", "inconsistency between DATE and AliRawEvent headers");
      MakeZombie();
      return;
   }

   if (create) {
      if (!Create())
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

   TString fs = fstoggle ? kRawDBFS[1] : kRawDBFS[0];
   TDatime dt;

   TString hostname = gSystem->HostName();
   Int_t pos;
   if ((pos = hostname.Index(".")) != kNPOS)
      hostname.Remove(pos);

   if (!FSHasSpace(fs)) {
      while (1) {
         fstoggle = !fstoggle;
         fs = fstoggle ? kRawDBFS[1] : kRawDBFS[0];
         if (FSHasSpace(fs)) break;
         Info("GetFileName", "sleeping 30 seconds before retrying...");
         gSystem->Sleep(30000);   // sleep for 30 seconds
         if (gAliMDC && gAliMDC->StopLoop())
            return 0;
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
Bool_t AliRawDB::Create()
{
   // Create a new raw DB.

   const Int_t kMaxRetry = 200;
   const Int_t kMaxSleep = 1;      // seconds
   const Int_t kMaxSleepLong = 10; // seconds
   Int_t retry = 0;

again:
   if (gAliMDC && gAliMDC->StopLoop())
      return kFALSE;

   const char *fname = GetFileName();
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
}

//______________________________________________________________________________
void AliRawDB::Close()
{
   // Close raw DB.

   if (!fRawDB) return;

   fRawDB->cd();

   // Write the tree.
   fTree->Write();

   // Close DB, this also deletes the fTree
   fRawDB->Close();

   if (AliMDC::DeleteFiles()) {
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
Bool_t AliRawDB::NextFile()
{
   // Close te current file and open a new one.
   // Returns kFALSE in case opening failed.

   Close();

   if (!Create()) return kFALSE;
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
AliRawRFIODB::AliRawRFIODB(AliRawEvent *event, Double_t maxsize, Int_t compress)
   : AliRawDB(event, maxsize, compress, kFALSE)
{
   // Create a new raw DB that will be accessed via RFIO.

#ifndef USE_RDM
   static int init = 0;
   // Set STAGE_POOL environment variable to current host
   if (!init) {
      // THESE ENVIRONMENT SYMBOLS ARE NOW DEFINED BY THE ALICE DATE SETUP
      // THEREFORE WE SHALL NOT USE ANY HARDCODED VALUES BUT RATHER USE
      // WHATEVER HAS BEEN SET IN THE DATE SITE
      //gSystem->Setenv("STAGE_POOL", "lcg00");
      //gSystem->Setenv("STAGE_HOST", "stage013");

      // however for sanity we check if they are really set
      if (!gSystem->Getenv("STAGE_POOL"))
         Error("AliRawRFIODB", "STAGE_POOL not set");
      if (!gSystem->Getenv("STAGE_HOST"))
         Error("AliRawRFIODB", "STAGE_HOST not set");
      init = 1;
   }
#endif

   if (!Create())
      MakeZombie();
   else
      fRawDB->UseCache(50, 0x200000);  //0x100000 = 1MB)
}

//______________________________________________________________________________
const char *AliRawRFIODB::GetFileName() const
{
   // Return filename based on hostname and date and time. This will make
   // each file unique. Also the directory will be made unique for each
   // day by adding the date to the fs. Assumes there is always enough
   // space on the device.

   static TString fname;

   TString fs = kRFIOFS;
   TDatime dt;

   // make a new subdirectory for each day
   fs += "/adc-";
   fs += dt.GetDate();

   Long_t id, size, flags, time;
   if (gSystem->GetPathInfo(fs, &id, &size, &flags, &time) == 1) {
      // directory does not exist, create it
      if (gSystem->mkdir(fs, kTRUE) == -1) {
         Error("GetFileName", "cannot create dir %s, using %s", fs.Data(),
               kRFIOFS);
         fs = kRFIOFS;
      }
   }
   // FIXME: should check if fs is a directory

   TString hostname = gSystem->HostName();
   Int_t pos;
   if ((pos = hostname.Index(".")) != kNPOS)
      hostname.Remove(pos);

   fname = fs + "/" + hostname + "_";
   fname += dt.GetDate();
   fname += "_";
   fname += dt.GetTime();
   fname += ".root";

   return fname;
}

//______________________________________________________________________________
void AliRawRFIODB::Close()
{
   // Close raw RFIO DB.

   if (!fRawDB) return;

   fRawDB->cd();

   // Write the tree.
   fTree->Write();

   // Close DB, this also deletes the fTree
   fRawDB->Close();

   if (AliMDC::DeleteFiles()) {
      TUrl u(fRawDB->GetName());
      gSystem->Exec(Form("rfrm %s", u.GetFile()));
   }

   delete fRawDB;
   fRawDB = 0;
}


//______________________________________________________________________________
AliRawCastorDB::AliRawCastorDB(AliRawEvent *event, Double_t maxsize, Int_t compress)
   : AliRawDB(event, maxsize, compress, kFALSE)
{
   // Create a new raw DB that will be accessed via CASTOR and rootd.

#ifndef USE_RDM
   static int init = 0;
   // Set STAGE_POOL environment variable to current host
   if (!init) {
      // THESE ENVIRONMENT SYMBOLS ARE NOW DEFINED BY THE ALICE DATE SETUP
      // THEREFORE WE SHALL NOT USE ANY HARDCODED VALUES BUT RATHER USE
      // WHATEVER HAS BEEN SET IN THE DATE SITE
      //gSystem->Setenv("STAGE_POOL", "lcg00");
      //gSystem->Setenv("STAGE_HOST", "stage013");

      // however for sanity we check if they are really set
      if (!gSystem->Getenv("STAGE_POOL"))
         Error("AliRawRFIODB", "STAGE_POOL not set");
      if (!gSystem->Getenv("STAGE_HOST"))
         Error("AliRawRFIODB", "STAGE_HOST not set");
      init = 1;
   }
#endif

   if (!Create())
      MakeZombie();
   else
      fRawDB->UseCache(50, 0x200000);  //0x100000 = 1MB)
}

//______________________________________________________________________________
const char *AliRawCastorDB::GetFileName() const
{
   // Return filename based on hostname and date and time. This will make
   // each file unique. Also the directory will be made unique for each
   // day by adding the date to the fs. Assumes there is always enough
   // space on the device.

   static TString fname;

   TString fs  = kCastorFS;
   TString fsr = kRFIOFS;
   TDatime dt;

   // make a new subdirectory for each day
   fs += "/adc-";
   fs += dt.GetDate();

   fsr += "/adc-";
   fsr += dt.GetDate();

   Long_t id, size, flags, time;
   if (gSystem->GetPathInfo(fsr, &id, &size, &flags, &time) == 1) {
      // directory does not exist, create it
      if (gSystem->mkdir(fsr, kTRUE) == -1) {
         Error("GetFileName", "cannot create dir %s, using %s", fsr.Data(),
               kRFIOFS);
         fs = kCastorFS;
      }
   }
   // FIXME: should check if fs is a directory

   TString hostname = gSystem->HostName();
   Int_t pos;
   if ((pos = hostname.Index(".")) != kNPOS)
      hostname.Remove(pos);

   fname = fs + "/" + hostname + "_";
   fname += dt.GetDate();
   fname += "_";
   fname += dt.GetTime();
   fname += ".root";

   return fname;
}

//______________________________________________________________________________
void AliRawCastorDB::Close()
{
   // Close raw CASTOR/rootd DB.

   if (!fRawDB) return;

   fRawDB->cd();

   // Write the tree.
   fTree->Write();

   // Close DB, this also deletes the fTree
   fRawDB->Close();

   if (AliMDC::DeleteFiles()) {
      TUrl u(fRawDB->GetName());
      gSystem->Exec(Form("rfrm %s", u.GetFile()));
   }

   delete fRawDB;
   fRawDB = 0;
}


//______________________________________________________________________________
AliRawRootdDB::AliRawRootdDB(AliRawEvent *event, Double_t maxsize, Int_t compress)
   : AliRawDB(event, maxsize, compress, kFALSE)
{
   // Create a new raw DB that will be accessed via rootd daemon.

   if (!Create())
      MakeZombie();
   else
      fRawDB->UseCache(50, 0x200000);  //0x100000 = 1MB)
}

//______________________________________________________________________________
const char *AliRawRootdDB::GetFileName() const
{
   // Return filename based on hostname and date and time. This will make
   // each file unique. Also the directory will be made unique for each
   // day by adding the date to the fs. Assumes there is always enough
   // space on the device.

   static TString fname;

   TString fs = kRootdFS;
   TDatime dt;

#if 0
   // make a new subdirectory for each day
   fs += "/adc-";
   fs += dt.GetDate();

   Long_t id, size, flags, time;
   if (gSystem->GetPathInfo(fs, &id, &size, &flags, &time) == 1) {
      // directory does not exist, create it
      if (gSystem->mkdir(fs, kTRUE) == -1) {
         Error("GetFileName", "cannot create dir %s, using %s", fs.Data(),
               kRootdFS);
         fs = kRootdFS;
      }
   }
   // FIXME: should check if fs is a directory
#endif

   TString hostname = gSystem->HostName();
   Int_t pos;
   if ((pos = hostname.Index(".")) != kNPOS)
      hostname.Remove(pos);

   fname = fs + "/" + hostname + "_";
   fname += dt.GetDate();
   fname += "_";
   fname += dt.GetTime();
   fname += ".root";

   return fname;
}

//______________________________________________________________________________
void AliRawRootdDB::Close()
{
   // Close raw rootd DB.

   if (!fRawDB) return;

   fRawDB->cd();

   // Write the tree.
   fTree->Write();

   // Close DB, this also deletes the fTree
   fRawDB->Close();

#if 0
   // can use services of TFTP
   if (AliMDC::DeleteFiles())
      gSystem->Exec(Form("rfrm %s", fRawDB->GetName()));
#endif

   delete fRawDB;
   fRawDB = 0;
}


//______________________________________________________________________________
AliRawNullDB::AliRawNullDB(AliRawEvent *event, Double_t maxsize, Int_t compress)
   : AliRawDB(event, maxsize, compress, kFALSE)
{
   // Create a new raw DB that will wrtie to /dev/null.

   if (!Create())
      MakeZombie();
}

//______________________________________________________________________________
const char *AliRawNullDB::GetFileName() const
{
   // Return /dev/null as filename.

   return "/dev/null";
}

//______________________________________________________________________________
void AliRawNullDB::Close()
{
   // Close raw RFIO DB.

   if (!fRawDB) return;

   fRawDB->cd();

   // Write the tree.
   fTree->Write();

   // Close DB, this also deletes the fTree
   fRawDB->Close();

   delete fRawDB;
   fRawDB = 0;
}


//______________________________________________________________________________
AliTagDB::AliTagDB(AliRawEventHeader *header, Double_t maxsize, Bool_t create)
{
   // Create tag DB.

   fHeader   = header;
   fMaxSize  = maxsize;

   if (create) {
      if (!Create())
         MakeZombie();
   }
}

//______________________________________________________________________________
Bool_t AliTagDB::Create()
{
   // Create a new tag DB.

   fTagDB = new TFile(GetFileName(), "RECREATE",
                      Form("ALICE MDC%d tag DB", kMDC), 1);
   if (fTagDB->IsZombie()) {
      Error("Create", "error opening tag DB");
      fTagDB = 0;
      return kFALSE;
   }

   // Create ROOT Tree object container
   fTree = new TTree("TAG", Form("ALICE MDC%d header data tree", kMDC));
   fTree->SetAutoSave(100000000);  // autosave when 100 Mbyte written

   Int_t bufsize = 32000;
   Int_t split   = 1;
   fTree->Branch("header", "AliRawEventHeader", &fHeader, bufsize, split);

   return kTRUE;
}

//______________________________________________________________________________
void AliTagDB::Close()
{
   // Close tag DB.

   if (!fTagDB) return;

   fTagDB->cd();

   // Write the tree.
   fTree->Write();

   // Close DB, this also deletes the fTree
   fTagDB->Close();

   if (AliMDC::DeleteFiles())
      gSystem->Unlink(fTagDB->GetName());

   delete fTagDB;
   fTagDB = 0;
}

//______________________________________________________________________________
Bool_t AliTagDB::NextFile()
{
   // Close te current file and open a new one.
   // Returns kFALSE in case opening failed.

   Close();

   if (!Create()) return kFALSE;
   return kTRUE;
}

//______________________________________________________________________________
Float_t AliTagDB::GetCompressionFactor() const
{
   // Return compression factor.

   if (fTree->GetZipBytes() == 0.)
      return 1.0;
   else
      return fTree->GetTotBytes()/fTree->GetZipBytes();
}

//______________________________________________________________________________
const char *AliTagDB::GetFileName() const
{
   // Return filename based on hostname and date and time. This will make
   // each file unique. The tags will be stored in the /data1/tags directory.

   static char fname[64];
   const char *fs = kTagDBFS;

   // check that fs exists (crude check fails if fs is a file)
   gSystem->MakeDirectory(fs);

   char hostname[64];

   strcpy(hostname, gSystem->HostName());

   char *s;
   if ((s = strchr(hostname, '.')))
      *s = 0;

   TDatime dt;

   sprintf(fname, "%s/%s_%d_%d.root", fs, hostname, dt.GetDate(), dt.GetTime());

   return fname;
}


//______________________________________________________________________________
AliTagNullDB::AliTagNullDB(AliRawEventHeader *header, Double_t maxsize) :
   AliTagDB(header, maxsize, kFALSE)
{
   // Create tag db writing to /dev/null.

   if (!Create())
      MakeZombie();
}

//______________________________________________________________________________
const char *AliTagNullDB::GetFileName() const
{
   // Return /dev/null as filename.

   return "/dev/null";
}

//______________________________________________________________________________
void AliTagNullDB::Close()
{
   // Close null tag DB.

   if (!fTagDB) return;

   fTagDB->cd();

   // Write the tree.
   fTree->Write();

   // Close DB, this also deletes the fTree
   fTagDB->Close();

   delete fTagDB;
   fTagDB = 0;
}


//______________________________________________________________________________
AliRunDB::AliRunDB(Bool_t noLocalDB)
{
   // Open run database, and get or create tree.

   fRunDB = 0;

   if (noLocalDB) return;

   // Get hostname
   char hostname[64], filename[64];
   const char *fs = kRunDBFS;

   // check that fs exists (crude check fails if fs is a file)
   gSystem->MakeDirectory(fs);

   strcpy(hostname, gSystem->HostName());

   char *s;
   if ((s = strchr(hostname, '.')))
      *s = 0;

   sprintf(filename, "%s/%s_rundb.root", fs, hostname);

   if (!gSystem->AccessPathName(filename, kFileExists))
      fRunDB = new TFile(filename, "UPDATE");
   else
      fRunDB = new TFile(filename, "CREATE", Form("ALICE MDC%d Run DB", kMDC));
}

//______________________________________________________________________________
void AliRunDB::Update(AliStats *stats)
{
   // Add stats object to database.

   if (!stats || !fRunDB) return;

   TDirectory *ds = gDirectory;
   fRunDB->cd();

   char sname[64];
   char *s = (char*)strrchr(stats->GetFileName(), '/');
   if (s) {
      s++;
      strcpy(sname, s);
   } else
      strcpy(sname, stats->GetFileName());
   s = strchr(sname, '.');
   if (s) *s = 0;

   stats->Write(sname);

   ds->cd();
}

//______________________________________________________________________________
void AliRunDB::UpdateRDBMS(AliStats *stats)
{
   // Add stats object to central MySQL DB.

   if (!stats) return;

   char sql[4096];
   char bt[25], et[25];

   strcpy(bt, stats->GetBeginTime().AsSQLString());
   strcpy(et, stats->GetEndTime().AsSQLString());

   sprintf(sql, "INSERT INTO mdc%dcatalog VALUES (0, '%s', %d, "
           "%d, %d, %d, %d, %d, %d, %.2f, '%s', '%s', '%s')", kMDC,
           stats->GetFileName(), (int)stats->GetFileSize(), stats->GetEvents(),
           stats->GetFirstRun(), stats->GetFirstEvent(), stats->GetLastRun(),
           stats->GetLastEvent(), stats->GetCompressionMode(),
           stats->GetCompressionFactor(), stats->GetFilterState() ? "on" : "off",
           bt, et);

   // open connection to MySQL server on pcsalo
   TSQLServer *db = TSQLServer::Connect("mysql://pcsalo.cern.ch/mdc", "alice", "amdc");

   if (!db || db->IsZombie()) {
      Error("UpdateRDBMS", "failed to connect to MySQL server on pcsalo");
      printf("%s\n", sql);
      delete db;
      return;
   }

   TSQLResult *res = db->Query(sql);

   if (!res) {
      Error("UpdateRDBMS", Form("insert into mdc%dcatalog failed", kMDC));
      printf("%s\n", sql);
   }

   delete res;
   delete db;
}

//______________________________________________________________________________
void AliRunDB::UpdateAliEn(AliStats *stats)
{
   // Record file in AliEn catalog.

   if (!stats) return;

   TGrid *g = TGrid::Connect(kAlienHost, "");

   TString lfn = kAlienDir;
   TDatime dt;

   // make a subdirectory for each day
   lfn += "/adc-";
   lfn += dt.GetDate();

   // check if directory exists, if not create it
   Grid_ResultHandle_t res = 0;
   if (!(res = g->OpenDir(lfn))) {
      // directory does not exist, create it
      if (g->Mkdir(lfn) == -1) {
         Error("UpdateAliEn", "cannot create directory %s", lfn.Data());
         lfn = kAlienDir;
      }
   }
   if (res) g->CloseResult(res);

   lfn += "/";
   lfn += gSystem->BaseName(stats->GetFileName());

   Int_t result = g->AddFile(lfn, stats->GetFileName(),
			     (int)stats->GetFileSize());

   if (result == -1) {
      Error("UpdateAliEn", "error adding file to AliEn catalog");
      printf("AliEn: AddFile(%s, %s, %d)\n", lfn.Data(), stats->GetFileName(),
             (int)stats->GetFileSize());
   }

   delete g;
}

//______________________________________________________________________________
void AliRunDB::Close()
{
   // Close run database.

   if (fRunDB) fRunDB->Close();
   delete fRunDB;
}

//----------------- Use SIGUSR1 to interupt endless loop -----------------------
class AliMDCInterruptHandler : public TSignalHandler {
public:
   AliMDCInterruptHandler(AliMDC *mdc) : TSignalHandler(kSigUser1, kFALSE), fMDC(mdc) { }
   Bool_t Notify() {
      Info("Notify", "received a SIGUSR1 signal");
      fMDC->SetStopLoop();
      return kTRUE;
   }
private:
   AliMDC *fMDC;   // alimdc to signal

   AliMDCInterruptHandler(const AliMDCInterruptHandler &);
   AliMDCInterruptHandler &operator=(const AliMDCInterruptHandler &);
};

//______________________________________________________________________________
AliMDC::AliMDC(Int_t fd, Int_t compress, Double_t maxFileSize, Bool_t useFilter,
               EWriteMode mode, Bool_t useLoop, Bool_t delFiles)
{
   // Create MDC processor object.

   fFd           = fd;
   fCompress     = compress;
   fMaxFileSize  = maxFileSize;
   fUseFilter    = useFilter;
   fWriteMode    = mode;
   fUseLoop      = useLoop;
   fUseFifo      = kFALSE;
   fUseEb        = kFALSE;
   fStopLoop     = kFALSE;
   fNumEvents    = 0;
   fDebugLevel   = 0;
   fgDeleteFiles = delFiles;

   if (fFd == -1) {
#ifdef USE_EB
     if (!ebRegister()) {
        Error("AliMDC", "cannot register with the event builder (%s)",
              ebGetLastError());
        return;
     }
     fUseEb = kTRUE;
#else
     if ((mkfifo(kFifo, 0644) < 0) && (errno != EEXIST)) {
         Error("AliMDC", "cannot create fifo %s", kFifo);
         return;
      }
      if ((chmod(kFifo, 0666) == -1) && (errno != EPERM)) {
         Error("AliMDC", "cannot change permission of fifo %s", kFifo);
         return;
      }
      if ((fFd = open(kFifo, O_RDONLY)) == -1) {
         Error("AliMDC", "cannot open input file %s", kFifo);
         return;
      }
      fUseFifo = kTRUE;
#endif
      fUseLoop = kFALSE;
   }

   printf("<AliMDC::AliMDC>: input = %s, rawdb size = %f, filter = %s, "
          "looping = %s, compression = %d, delete files = %s",
          fUseFifo ? "fifo" : (fUseEb ? "eb" : "file"), fMaxFileSize,
          fUseFilter ? "on" : "off", fUseLoop ? "yes" : "no", fCompress,
          fgDeleteFiles ? "yes" : "no");
   if (fWriteMode == kRFIO)
      printf(", use RFIO\n");
   else if (fWriteMode == kROOTD)
      printf(", use rootd\n");
   else if (fWriteMode == kCASTOR)
      printf(", use CASTOR/rootd\n");
   else if (fWriteMode == kDEVNULL)
      printf(", write raw data to /dev/null\n");
   else
      printf("\n");

   // install SIGUSR1 handler to allow clean interrupts
   gSystem->AddSignalHandler(new AliMDCInterruptHandler(this));

   gAliMDC = this;
}

//______________________________________________________________________________
Int_t AliMDC::Run()
{
   // Run the MDC processor. Read from the input stream and only return
   // when the input gave and EOF or a fatal error occured. On success 0
   // is returned, 1 in case of a fatality.

   TStopwatch timer;
   Int_t status;

   // Make sure needed directories exist
   const char *dirs[4];
   dirs[0] = kRawDBFS[0];
   dirs[1] = kRawDBFS[1];
   dirs[2] = kTagDBFS;
   dirs[3] = kRunDBFS;
   for (int idir = 0; idir < 4; idir++) {
      gSystem->ResetErrno();
      gSystem->MakeDirectory(dirs[idir]);
      if (gSystem->GetErrno() && gSystem->GetErrno() != EEXIST) {
         SysError("Run", "mkdir %s", dirs[idir]);
         return 1;
      }
   }

   // Used for statistics
   timer.Start();
   Double_t told = 0, tnew = 0;
   Float_t  chunkSize = fMaxFileSize/100, nextChunk = chunkSize;

   // Event object used to store event data.
   AliRawEvent *event = new AliRawEvent;

   // Create new raw DB.
   AliRawDB *rawdb;
   if (fWriteMode == kRFIO)
      rawdb = new AliRawRFIODB(event, fMaxFileSize, fCompress);
   else if (fWriteMode == kROOTD)
      rawdb = new AliRawRootdDB(event, fMaxFileSize, fCompress);
   else if (fWriteMode == kCASTOR)
      rawdb = new AliRawCastorDB(event, fMaxFileSize, fCompress);
   else if (fWriteMode == kDEVNULL)
      rawdb = new AliRawNullDB(event, fMaxFileSize, fCompress);
   else
      rawdb = new AliRawDB(event, fMaxFileSize, fCompress);

   if (rawdb->IsZombie()) return 1;
   printf("Filling raw DB %s\n", rawdb->GetDBName());

   // Create new tag DB.
   AliTagDB *tagdb = 0;
#if 0
   // no tagdb for the time being to get maximum speed
   if (fWriteMode == kDEVNULL)
      tagdb = new AliTagNullDB(event->GetHeader(), kMaxTagFileSize);
   else
      tagdb = new AliTagDB(event->GetHeader(), kMaxTagFileSize);
   if (tagdb->IsZombie())
      tagdb = 0;
   else
      printf("Filling tag DB %s\n", tagdb->GetDBName());
#endif

   // Create AliStats object
   AliStats *stats = new AliStats(rawdb->GetDBName(), fCompress, fUseFilter);

   // Shortcut for easy header access
   AliRawEventHeader &header = *event->GetHeader();

   // Process input stream
#ifdef USE_EB
   Int_t eorFlag = 0;
   while (!(eorFlag = ebEor())) {
      struct iovec *ebvec;
      if ((ebvec = ebGetNextEvent()) == (void *)-1) {
         Error("Run", "error getting next event (%s)", ebGetLastError());
         break;
      }
      if (ebvec == 0) {
         // no event, sleep for 1 second and try again
         gSystem->Sleep(1000);
         continue;
      }
      char *ebdata = (char *) ebvec[0].iov_base;
#else
   while (1) {
      char *ebdata = 0;
#endif

      // Read event header
      if ((status = ReadHeader(header, ebdata)) != header.HeaderSize()) {
         if (status == 0) {
            if (fUseLoop) {
#ifndef USE_EB
               ::lseek(fFd, 0, SEEK_SET);
#endif
               continue;
            }
            printf("<AliMDC::Run>: EOF, processed %d events\n", fNumEvents);
            break;
         }
         return 1;
      }
      ALIDEBUG(3)
         header.Dump();

      // If we were in looping mode stop directly after a SIGUSR1 signal
      if (StopLoop()) {
         Info("Run", "Stopping loop, processed %d events", fNumEvents);
         break;
      }

      // Check if event has any hard track flagged
      Bool_t callFilter = kFALSE;
      // This needs to be re-engineered for the next ADC...
      //if (fUseFilter && TEST_USER_ATTRIBUTE(header.GetTypeAttribute(), 0))
      //   callFilter = kTRUE;

      // Check event type and skip "Start of Run", "End of Run",
      // "Start of Run Files" and "End of Run Files"
      switch (header.GetType()) {
         case kStartOfRun:
         case kEndOfRun:
         case kStartOfRunFiles:
         case kEndOfRunFiles:
            {
               Int_t skip = header.GetEventSize() - header.HeaderSize();
#ifndef USE_EB
               ::lseek(fFd, skip, SEEK_CUR);
#endif
               ALIDEBUG(1)
                  Info("Run", "Skipping %s (%d bytes)", header.GetTypeName(), skip);
               continue;
            }
         default:
            ALIDEBUG(1) {
               Int_t s = header.GetEventSize() - header.HeaderSize();
               Info("Run", "Processing %s (%d bytes)", header.GetTypeName(), s);
            }
      }

      // Amount of data left to read for this event
      Int_t toRead = header.GetEventSize() - header.HeaderSize();

      // If there is less data for this event than the next sub-event
      // header, something is wrong. Skip to next event...
      if (toRead < header.HeaderSize()) {
	 ALIDEBUG(1) {
            Warning("Run",
		    "header size (%d) exceeds number of bytes to read (%d)\n",
		    header.HeaderSize(), toRead);
	    header.Dump();
         }
         if ((status = DumpEvent(toRead)) != toRead) {
            if (status == 0)
               break;
            return 1;
         }
         Error("Run", "discarding event %d (too little data for header)", fNumEvents);
         continue;
      }

      // Loop over all sub-events... (LDCs)
      Int_t nsub = 1;
      while (toRead > 0) {
#ifdef USE_EB
         ebdata = (char *)ebvec[nsub].iov_base;
#endif

         ALIDEBUG(1)
            Info("Run", "reading LDC %d", nsub);

         AliRawEvent *subEvent = event->NextSubEvent();

         // Read sub-event header
         AliRawEventHeader &subHeader = *subEvent->GetHeader();
         if ((status = ReadHeader(subHeader, ebdata)) != subHeader.HeaderSize()) {
            if (status == 0) {
               Error("Run", "unexpected EOF reading sub-event header");
               break;
            }
            return 1;
         }

         ALIDEBUG(3)
            subHeader.Dump();

         toRead -= subHeader.HeaderSize();

#ifdef USE_EB
         ebdata = (char *)(ebvec[nsub].iov_base) + subHeader.HeaderSize();
#endif

         Int_t rawSize = subHeader.GetEventSize() - subHeader.HeaderSize();

         // Read Equipment Header (in case of physics or calibration event)
         if (header.GetType() == kPhysicsEvent ||
             header.GetType() == kCalibrationEvent) {
            AliRawEquipmentHeader &equipment = *subEvent->GetEquipmentHeader();
            Int_t equipHeaderSize = equipment.HeaderSize();
            if ((status = ReadEquipmentHeader(equipment, header.DataIsSwapped(),
                                              ebdata)) != equipHeaderSize) {
               if (status == 0) {
                  Error("Run", "unexpected EOF reading equipment-header");
                  break;
               }
               return 1;
            }
            toRead  -= equipHeaderSize;
            rawSize -= equipHeaderSize;
#ifdef USE_EB
            ebdata = (char *)(ebvec[nsub].iov_base) + subHeader.HeaderSize() +
                     equipHeaderSize;
#endif
         }

         // Make sure raw data less than left over bytes for current event
         if (rawSize > toRead) {
            ALIDEBUG(1) {
               Warning("Run", "raw data size (%d) exceeds number of bytes "
		       "to read (%d)\n", rawSize, toRead);
	       subHeader.Dump();
            }
            if ((status = DumpEvent(toRead)) != toRead) {
               if (status == 0)
                  break;
               return 1;
            }
            Error("Run", "discarding event %d (too much data)", fNumEvents);
            continue;
         }

         // Read sub-event raw data
         AliRawData &subRaw = *subEvent->GetRawData();
         if ((status = ReadRawData(subRaw, rawSize, ebdata)) != rawSize) {
            if (status == 0) {
               Error("Run", "unexpected EOF reading sub-event raw data");
               break;
            }
            return 1;
         }

         if (callFilter) {
            if (TEST_USER_ATTRIBUTE(subHeader.GetTypeAttribute(), 0))
               Filter(subRaw);
            else {
               // set size of all sectors without hard track flag to 0
               subRaw.SetSize(0);
            }
         }

         toRead -= rawSize;
         nsub++;
      }

      // Set stat info for first event of this file
      if (rawdb->GetEvents() == 0)
         stats->SetFirstId(header.GetRunNumber(), header.GetEventInRun());

      // Store raw event in tree
      rawdb->Fill();

      // Store header in tree
      if (tagdb) tagdb->Fill();

      fNumEvents++;

      if (!(fNumEvents%10))
         printf("Processed event %d (%d)\n", fNumEvents, rawdb->GetEvents());

      // Filling time statistics
      if (rawdb->GetBytesWritten() > nextChunk) {
         tnew = timer.RealTime();
         stats->Fill(tnew-told);
         told = tnew;
         timer.Continue();
         nextChunk += chunkSize;
      }

      // Check size of raw db. If bigger than maxFileSize, close file
      // and continue with new file.
      if (rawdb->FileFull()) {

         printf("Written raw DB at a rate of %.1f MB/s\n",
                rawdb->GetBytesWritten() / timer.RealTime() / 1000000.);

         // Write stats object to raw db, run db, MySQL and AliEn
         stats->WriteToDB(rawdb);
         delete stats;

         if (!rawdb->NextFile()) {
            Error("Run", "error opening next raw data file");
            return 1;
         }

         printf("Filling raw DB %s\n", rawdb->GetDBName());
         stats = new AliStats(rawdb->GetDBName(), fCompress, fUseFilter);

         timer.Start();
         told = 0, tnew = 0;
         nextChunk = chunkSize;
      }

      // Check size of tag db
      if (tagdb && tagdb->FileFull()) {
         if (!tagdb->NextFile())
            tagdb = 0;
         else
            printf("Filling tag DB %s\n", tagdb->GetDBName());
      }

      // Make top event object ready for next event data
      //printf("Event %d has %d sub-events\n", fNumEvents, event->GetNSubEvents());
      event->Reset();

#ifdef USE_EB
      if (!ebReleaseEvent(ebvec)) {
         Error("Run", "problem releasing event (%s)", ebGetLastError());
         break;
      }
#endif
   }

   printf("Written raw DB at a rate of %.1f MB/s\n",
          rawdb->GetBytesWritten() / timer.RealTime() / 1000000.);

   // Write stats to raw db and run db and delete stats object
   stats->WriteToDB(rawdb);
   delete stats;

   // Close the raw DB
   delete rawdb;

   // Close the tag DB
   delete tagdb;

   // Close input source
   close(fFd);

#if 0
   // Cleanup fifo
   if (fUseFifo && ::unlink(kFifo) == -1) {
      SysError("Run", "unlink");
      return 1;
   }
#endif

#ifdef USE_EB
   // Print eor flag
   if (eorFlag) {
      Info("Run", "event builder reported end of run (%d)", eorFlag);
   }
#endif

   return 0;
}

//______________________________________________________________________________
Int_t AliMDC::Read(void *buffer, Int_t length)
{
   // Read exactly length bytes into buffer. Returns number of bytes
   // received, returns -1 in case of error and 0 for EOF.

   errno = 0;

   if (fFd < 0) return -1;

   Int_t n, nrecv = 0;
   char *buf = (char *)buffer;

   for (n = 0; n < length; n += nrecv) {
      if ((nrecv = read(fFd, buf+n, length-n)) <= 0) {
         if (nrecv == 0)
            break;        // EOF
         if (errno != EINTR)
            SysError("Read", "read");
         return -1;
      }
   }
   return n;
}

//______________________________________________________________________________
Int_t AliMDC::ReadHeader(AliRawEventHeader &header, void *eb)
{
   // Read header info from DATE data stream. Returns bytes read (i.e.
   // AliRawEventHeader::HeaderSize()), -1 in case of error and 0 for EOF.

   Int_t nrecv;

   if (eb) {
      // read from event builder memory area
      memcpy(header.HeaderBegin(), eb, header.HeaderSize());
      nrecv = header.HeaderSize();
   } else {
      // read from fifo or file
      if ((nrecv = Read(header.HeaderBegin(), header.HeaderSize())) !=
           header.HeaderSize()) {
         if (nrecv == 0)
            return 0;
         return -1;
      }
   }

   // Swap header data if needed
   if (header.IsSwapped())
      header.Swap();

   // Is header valid...
   if (!header.IsValid()) {
      Error("ReadHeader", "invalid header format");
      // try recovery... how?
      return -1;
   }
   if (header.GetEventSize() < (UInt_t)header.HeaderSize()) {
      Error("ReadHeader", "invalid header size");
      // try recovery... how?
      return -1;
   }

   return nrecv;
}

//______________________________________________________________________________
Int_t AliMDC::ReadEquipmentHeader(AliRawEquipmentHeader &header,
                                  Bool_t isSwapped, void *eb)
{
   // Read equipment header info from DATE data stream. Returns bytes read
   // (i.e. AliRawEquipmentHeader::HeaderSize()), -1 in case of error and
   // 0 for EOF. If isSwapped is kTRUE the event data is byte swapped
   // and we will swap the header to host format.

   Int_t nrecv;

   if (eb) {
      // read from event builder memory area
      memcpy(header.HeaderBegin(), eb, header.HeaderSize());
      nrecv = header.HeaderSize();
   } else {
      // read from fifo or file
      if ((nrecv = Read(header.HeaderBegin(), header.HeaderSize())) !=
           header.HeaderSize()) {
         if (nrecv == 0)
            return 0;
         return -1;
      }
   }

   // Swap equipment header data if needed
   if (isSwapped)
      header.Swap();

   if (header.GetEquipmentSize() < (UInt_t)header.HeaderSize()) {
      Error("ReadEquipmentHeader", "invalid equipment header size");
      // try recovery... how?
      return -1;
   }

   return nrecv;
}

//______________________________________________________________________________
Int_t AliMDC::ReadRawData(AliRawData &raw, Int_t size, void *eb)
{
   // Read raw data from DATE data stream. Returns bytes read (i.e.
   // AliRawEventHeader::HeaderSize()), -1 in case of error and 0 for EOF.

   Int_t nrecv;

   if (eb) {
      // read from event builder memory area
      raw.SetBuffer(eb, size);
      nrecv = size;
   } else {
      // read from fifo or file
      raw.SetSize(size);
      if ((nrecv = Read(raw.GetBuffer(), size)) != size) {
         if (nrecv == 0) {
            Error("ReadRawData", "unexpected EOF");
            return 0;
         }
         return -1;
      }
   }

   return nrecv;
}

//______________________________________________________________________________
Int_t AliMDC::DumpEvent(Int_t toRead)
{
   // This case should not happen, but if it does try to handle it
   // gracefully by reading the rest of the event and discarding it.
   // Returns bytes read, -1 in case of fatal error and 0 for EOF.

   Error("DumpEvent", "dumping %d bytes of event %d", toRead, fNumEvents);

   Int_t nrecv;
   char *tbuf = new char[toRead];
   if ((nrecv = Read(tbuf, toRead)) != toRead) {
      if (nrecv == 0) {
         Error("DumpEvent", "unexpected EOF");
         return 0;
      }
      return -1;
   }
   delete [] tbuf;

   return nrecv;
}

//______________________________________________________________________________
Int_t AliMDC::Filter(AliRawData &raw)
{
   // Call 3rd level filter for this raw data segment.

#ifdef USE_HLT

   // Add HLT code here

#else

   raw.GetSize();
   printf("Filter called for event %d\n", fNumEvents);

#endif

   return 0;
}
