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

/* $Id$ */

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// AliMDC                                                               //
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

#include <errno.h>

#include <TSystem.h>
#include <TError.h>
#include <TStopwatch.h>

#ifdef ALI_DATE
#include "event.h"
#endif
#ifdef USE_EB
#include "libDateEb.h"
#endif

#include "AliRawEvent.h"
#include "AliRawEventHeader.h"
#include "AliRawEquipmentHeader.h"
#include "AliRawData.h"
#include "AliStats.h"
#include "AliRawDB.h"
#include "AliRawRFIODB.h"
#include "AliRawCastorDB.h"
#include "AliRawRootdDB.h"
#include "AliRawNullDB.h"
#include "AliTagDB.h"

#include "AliMDC.h"


ClassImp(AliMDC)


#define ALIDEBUG(level) \
   if (AliMDC::Instance() && (AliMDC::Instance()->GetDebugLevel() >= (level)))


// Fixed file system locations for the different DB's
#ifdef USE_RDM
const char* const AliMDC::fgkFifo       = "/tmp/alimdc.fifo";
const char* const AliMDC::fgkRawDBFS[2] = { "/tmp/mdc1", "/tmp/mdc2" };
const char* const AliMDC::fgkTagDBFS    = "/tmp/mdc1/tags";
const char* const AliMDC::fgkRunDBFS    = "/tmp/mdc1/meta";
const char* const AliMDC::fgkRFIOFS     = "rfio:/castor/cern.ch/user/r/rdm";
const char* const AliMDC::fgkCastorFS   = "castor:/castor/cern.ch/user/r/rdm";
const char* const AliMDC::fgkRootdFS    = "root://localhost//tmp/mdc1";
const char* const AliMDC::fgkAlienHost  = "alien://aliens7.cern.ch:15000/?direct";
const char* const AliMDC::fgkAlienDir   = "/alice_mdc/DC";
#else
const char* const AliMDC::fgkFifo       = "/tmp/alimdc.fifo";
const char* const AliMDC::fgkRawDBFS[2] = { "/data1/mdc", "/data2/mdc" };
const char* const AliMDC::fgkTagDBFS    = "/data1/mdc/tags";
const char* const AliMDC::fgkRunDBFS    = "/data1/mdc/meta";
const char* const AliMDC::fgkRFIOFS     = "rfio:/castor/cern.ch/lcg/dc5";
const char* const AliMDC::fgkCastorFS   = "castor:/castor/cern.ch/lcg/dc5";
const char* const AliMDC::fgkRootdFS    = "root://localhost//tmp/mdc1";
const char* const AliMDC::fgkAlienHost  = "alien://aliens7.cern.ch:15000/?direct";
const char* const AliMDC::fgkAlienDir   = "/alice_mdc/DC";
#endif

// Maximum size of tag db files
const Double_t AliMDC::fgkMaxTagFileSize = 2.5e8;    // 250MB

Bool_t AliMDC::fgDeleteFiles = kFALSE;
AliMDC* AliMDC::fgInstance = NULL;


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
     if ((mkfifo(fgkFifo, 0644) < 0) && (errno != EEXIST)) {
         Error("AliMDC", "cannot create fifo %s", fgkFifo);
         return;
      }
      if ((chmod(fgkFifo, 0666) == -1) && (errno != EPERM)) {
         Error("AliMDC", "cannot change permission of fifo %s", fgkFifo);
         return;
      }
      if ((fFd = open(fgkFifo, O_RDONLY)) == -1) {
         Error("AliMDC", "cannot open input file %s", fgkFifo);
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

   fgInstance = this;
}

//______________________________________________________________________________
AliMDC::AliMDC(const AliMDC& mdc): TObject(mdc)
{
// copy constructor

  Fatal("AliMDC", "copy constructor not implemented");
}

//______________________________________________________________________________
AliMDC& AliMDC::operator = (const AliMDC& /*mdc*/)
{
// assignment operator

  Fatal("operator =", "assignment operator not implemented");
  return *this;
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
   dirs[0] = fgkRawDBFS[0];
   dirs[1] = fgkRawDBFS[1];
   dirs[2] = fgkTagDBFS;
   dirs[3] = fgkRunDBFS;
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
   if (fWriteMode == fgkDEVNULL)
      tagdb = new AliTagNullDB(event->GetHeader(), fgkMaxTagFileSize);
   else
      tagdb = new AliTagDB(event->GetHeader(), fgkMaxTagFileSize);
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
         case AliRawEventHeader::kStartOfRun:
         case AliRawEventHeader::kEndOfRun:
         case AliRawEventHeader::kStartOfRunFiles:
         case AliRawEventHeader::kEndOfRunFiles:
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
         if (header.GetType() == AliRawEventHeader::kPhysicsEvent ||
             header.GetType() == AliRawEventHeader::kCalibrationEvent) {
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
#ifdef ALI_DATE
            if (TEST_USER_ATTRIBUTE(subHeader.GetTypeAttribute(), 0))
               Filter(subRaw);
            else {
               // set size of all sectors without hard track flag to 0
               subRaw.SetSize(0);
            }
#endif
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
   if (fUseFifo && ::unlink(fgkFifo) == -1) {
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

//______________________________________________________________________________
AliMDC::AliMDCInterruptHandler::AliMDCInterruptHandler(const 
						       AliMDCInterruptHandler&
						       handler): 
  TSignalHandler(handler) 
{
// copy constructor

  Fatal("AliMDCInterruptHandler", "copy constructor not implemented");
}

//______________________________________________________________________________
AliMDC::AliMDCInterruptHandler& 
  AliMDC::AliMDCInterruptHandler::operator = (const AliMDCInterruptHandler& 
					      /*handler*/)
{
// assignment operator

  Fatal("operator =", "assignment operator not implemented");
  return *this;
}
