#ifndef ALIMDC_H
#define ALIMDC_H
// @(#)alimdc:$Name$:$Id$
// Author: Fons Rademakers  26/11/99

/* Copyright(c) 1998-2003, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// AliMDC                                                               //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#ifndef ROOT_TObject
#include <TObject.h>
#endif

#ifndef ROOT_TSysEvtHandler
#include <TSysEvtHandler.h>
#endif

// Forward class declarations
class AliRawEvent;
class AliRawEventHeader;
class AliRawEquipmentHeader;
class AliRawData;
#ifdef USE_HLT
class AliESD;
#endif

class AliMDC : public TObject {

public:
   enum EWriteMode { kLOCAL, kRFIO, kROOTD, kCASTOR, kDEVNULL };

   AliMDC(Int_t fd, Int_t compress, Double_t maxFileSize, Bool_t useFilter,
          EWriteMode mode, Bool_t useLoop, Bool_t delFiles);
   ~AliMDC() { fgInstance = NULL; }

   static AliMDC* Instance() {return fgInstance;}

   Int_t  Run();
   void   SetStopLoop() { fStopLoop = kTRUE; }
   Bool_t StopLoop() const { return fStopLoop; }

   void   SetDebugLevel(Int_t level) { fDebugLevel = level; }
   Int_t  GetDebugLevel() const { return fDebugLevel; }

   static Bool_t DeleteFiles() { return fgDeleteFiles; }

   enum {kMDC = 6};  // Which MDC is this...

   static const char* Fifo() {return fgkFifo;}
   static const char* RawDBFS(Int_t i) {return fgkRawDBFS[i];}
   static const char* TagDBFS() {return fgkTagDBFS;}
   static const char* RunDBFS() {return fgkRunDBFS;}
   static const char* RFIOFS() {return fgkRFIOFS;}
   static const char* CastorFS() {return fgkCastorFS;}
   static const char* RootdFS() {return fgkRootdFS;}
   static const char* AlienHost() {return fgkAlienHost;}
   static const char* AlienDir() {return fgkAlienDir;}

private:
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

     AliMDCInterruptHandler(const AliMDCInterruptHandler& handler);
     AliMDCInterruptHandler& operator=(const AliMDCInterruptHandler& handler);
   };

   static AliMDC* fgInstance;  // singleton instance

   Int_t      fFd;          // DATE input stream
   Int_t      fCompress;    // compression factor used for raw output DB
   Int_t      fNumEvents;   // number of events processed
   Int_t      fDebugLevel;  // controls debug print-out
   Double_t   fMaxFileSize; // maximum size of raw output DB
   EWriteMode fWriteMode;   // write mode (local, rfio, rootd, castor, /dev/null)
   Bool_t     fUseFifo;     // read from fifo, file otherwise
   Bool_t     fUseEb;       // use event builder API instead of fifo
   Bool_t     fUseFilter;   // use 3rd level trigger filter
   Bool_t     fUseLoop;     // loop on input source (must be file)
   Bool_t     fStopLoop;    // break from endless loop (triggered by SIGUSR1)

   static Bool_t fgDeleteFiles;  // flag for deletion of files

   static const Double_t fgkMaxTagFileSize;  // maximal size of tag DB

   // Fixed file system locations for the different DB's
   static const char* const fgkFifo;       // fifo
   static const char* const fgkRawDBFS[2]; // raw DB
   static const char* const fgkTagDBFS;    // tag DB
   static const char* const fgkRunDBFS;    // run DB
   static const char* const fgkRFIOFS;     // rfio
   static const char* const fgkCastorFS;   // castor
   static const char* const fgkRootdFS;    // rootd
   static const char* const fgkAlienHost;  // alien host name
   static const char* const fgkAlienDir;   // alien directory

   AliMDC(const AliMDC& mdc);
   AliMDC& operator = (const AliMDC& mdc);

   Int_t     Read(const char *name) { return TObject::Read(name); }
   Int_t     Read(void *buffer, Int_t length);
   Int_t     ReadHeader(AliRawEventHeader &header, void *eb = 0);
   Int_t     ReadEquipmentHeader(AliRawEquipmentHeader &header,
                                 Bool_t isSwapped, void *eb = 0);
   Int_t     ReadRawData(AliRawData &raw, Int_t size, void *eb = 0);
   Int_t     DumpEvent(Int_t toRead);
   Int_t     Filter(
#ifdef USE_HLT
		    AliRawEvent *event,AliESD *esd
#endif
		    );

   ClassDef(AliMDC,0)  // MDC processor
};

#endif
