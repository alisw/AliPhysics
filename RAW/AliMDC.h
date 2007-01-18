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

#ifndef ROOT_TObjArray
#include <TObjArray.h>
#endif

#ifndef ROOT_TSysEvtHandler
#include <TSysEvtHandler.h>
#endif

// Forward class declarations
class AliRawEvent;
class AliRawEventHeaderBase;
class AliRawEquipmentHeader;
class AliRawData;
class AliRawDB;
class AliTagDB;
class AliRawEventTag;
class AliStats;
class AliESD;

class AliMDC : public TObject {

public:
   enum EWriteMode { kLOCAL, kRFIO, kROOTD, kCASTOR, kDEVNULL };
   enum EFilterMode { kFilterOff, kFilterTransparent, kFilterOn };
   enum EErrorCode { kErrStartEndRun = -1, 
		     kErrHeader = -2, 
		     kErrHeaderSize = -3, 
		     kErrSubHeader = -4, 
		     kErrDataSize = -5, 
		     kErrEquipmentHeader = -6, 
		     kErrEquipment = -7,
                     kErrFileSize = -8,
		     kFilterReject = -9,
                     kErrWriting = -10 };

   AliMDC(Int_t compress, Bool_t deleteFiles, 
	  EFilterMode filterMode = kFilterTransparent, 
	  Double_t maxSizeTagDB = -1, const char* fileNameTagDB = NULL);
   virtual ~AliMDC();

   Int_t      Open(EWriteMode mode, const char* fileName,
		   Double_t maxFileSize = 0,
		   const char* fs1 = NULL, const char* fs2 = NULL);
   Int_t      ProcessEvent(void* event, Bool_t isIovecArray = kFALSE);
   Int_t      GetTotalSize();
   Int_t      Close();

   Int_t      Run(const char* inputFile, Bool_t loop,
	          EWriteMode mode, Double_t maxFileSize, 
		  const char* fs1 = NULL, const char* fs2 = NULL);
   void       Stop();

private:
   class AliMDCInterruptHandler : public TSignalHandler {
   public:
     AliMDCInterruptHandler(AliMDC *mdc) : TSignalHandler(kSigUser1, kFALSE), fMDC(mdc) { }
     Bool_t Notify() {
       Info("Notify", "received a SIGUSR1 signal");
       fMDC->Stop();
       return kTRUE;
     }
   private:
     AliMDC *fMDC;   // alimdc to signal

     AliMDCInterruptHandler(const AliMDCInterruptHandler& handler); // Not implemented
     AliMDCInterruptHandler& operator=(const AliMDCInterruptHandler& handler); // Not implemented
   };

   AliRawEvent *fEvent;       // produced AliRawEvent
   AliESD      *fESD;         // pointer to HLT ESD object
   AliStats    *fStats;       // statistics
   AliRawDB    *fRawDB;       // raw data DB
   AliTagDB    *fTagDB;       // tag DB
   AliRawEventTag *fEventTag; // raw-data event tag object
   Int_t        fCompress;    // compression factor used for raw output DB
   Bool_t       fDeleteFiles; // flag for deletion of files
   EFilterMode  fFilterMode;  // high level filter mode
   TObjArray    fFilters;     // filter algorithms
   Bool_t       fStop;        // stop execution (triggered by SIGUSR1)
   Bool_t       fIsTagDBCreated; // is tag db already created
   Double_t     fMaxSizeTagDB;// max size of the tag DB
   const char*  fFileNameTagDB;// tag DB file name

   // Filter names
   enum {kNFilters = 1};
   static const char* const fgkFilterName[kNFilters];

   AliMDC(const AliMDC& mdc);
   AliMDC& operator = (const AliMDC& mdc);

   Int_t     Read(const char *name) { return TObject::Read(name); }
   Int_t     Read(Int_t fd, void *buffer, Int_t length);
   Int_t     ReadEquipmentHeader(AliRawEquipmentHeader &header,
                                 Bool_t isSwapped, char*& data);
   Int_t     ReadRawData(AliRawData &raw, Int_t size, char*& data);

   ClassDef(AliMDC,0)  // MDC processor
};

#endif
