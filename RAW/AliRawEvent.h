#ifndef ALIRAWEVENT_H
#define ALIRAWEVENT_H
// @(#)alimdc:$Name$:$Id$
// Author: Fons Rademakers  26/11/99
// Updated: Dario Favretto  15/04/2003

/* Copyright(c) 1998-2003, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

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

#ifndef ROOT_TObject
#include <TObject.h>
#endif

#ifndef ROOT_TDatime
#include <TDatime.h>
#endif

#ifndef ROOT_TString
#include <TString.h>
#endif

#ifndef ROOT_TFile
#include <TFile.h>
#endif

#ifndef ROOT_TTree
#include <TTree.h>
#endif


// Forward class declarations
class TH1F;
class AliRawDB;


// The following enumeration can be used once the kEventTypeMask has been
// applied to the raw event type
enum EAliRawEventType {
   kStartOfRun =       1,    // START_OF_RUN
   kEndOfRun =         2,    // END_OF_RUN
   kStartOfRunFiles =  3,    // START_OF_RUN_FILES
   kEndOfRunFiles =    4,    // END_OF_RUN_FILES
   kStartOfBurst =     5,    // START_OF_BURST
   kEndOfBurst =       6,    // END_OF_BURST
   kPhysicsEvent =     7,    // PHYSICS_EVENT
   kCalibrationEvent = 8,    // CALIBRATION_EVENT
   kFormatError =      9     // EVENT_FORMAT_ERROR
};

const Int_t kEventTypeMin = kStartOfRun;
const Int_t kEventTypeMax = kFormatError;

const UInt_t kEventMagicNumber        = 0xDA1E5AFE;
const UInt_t kEventMagicNumberSwapped = 0xFE5A1EDA;

// Type sizes
const Int_t kIdWords        = 2;
const Int_t kTriggerWords   = 2;
const Int_t kDetectorWords  = 1;
const Int_t kAttributeWords = 3;


class AliRawEventHeader : public TObject {

public:
   AliRawEventHeader() { fSize = 0; }
   virtual ~AliRawEventHeader() { }

   void         *HeaderBegin() { return (void *) &fSize; }
   Int_t         HeaderSize() const { return (Long_t) &fGDCId - (Long_t) &fSize + sizeof(fGDCId); }
   Bool_t        DataIsSwapped() const;
   Bool_t        IsSwapped() const { return (fMagic == kEventMagicNumberSwapped) ? kTRUE : kFALSE; }
   Bool_t        IsValid() const { return IsSwapped() ? kTRUE : ((fMagic == kEventMagicNumber) ? kTRUE : kFALSE); }
   void          Swap();

   UInt_t        GetEventSize() const { return fSize; }
   UInt_t        GetMagic() const { return fMagic; }
   UInt_t        GetHeaderLength() const { return fHeadLen; }
   UInt_t        GetVersion() const { return fVersion; }
   UInt_t        GetType() const { return fType; }
   const char   *GetTypeName() const;
   UInt_t        GetRunNumber() const { return fRunNb; }
   UInt_t        GetEventInRun() const;
   const UInt_t *GetId() const { return fId; }
   const UInt_t *GetTriggerPattern() const { return fTriggerPattern; }
   const UInt_t *GetDetectorPattern() const { return fDetectorPattern; }
   const UInt_t *GetTypeAttribute() const { return fTypeAttribute; }
   UInt_t        GetLDCId() const { return fLDCId; }
   UInt_t        GetGDCId() const { return fGDCId; }

private:
   UInt_t fSize;          // size of event in bytes
   UInt_t fMagic;         // magic number used for consistency check
   UInt_t fHeadLen;       // size of header in bytes
   UInt_t fVersion;       // unique version identifier
   UInt_t fType;          // event type
   UInt_t fRunNb;         // run number
   UInt_t fId[kIdWords];  // id field
   UInt_t fTriggerPattern[kTriggerWords];   // trigger pattern
   UInt_t fDetectorPattern[kDetectorWords]; // detector pattern
   UInt_t fTypeAttribute[kAttributeWords];  // system (0,1) and user (2) attributes
   UInt_t fLDCId;         // LDC id
   UInt_t fGDCId;         // GDC id

   ClassDef(AliRawEventHeader,1)  // Alice raw event header
};


class AliRawEquipmentHeader : public TObject {

public:
   AliRawEquipmentHeader() { fSize = 0; }
   ~AliRawEquipmentHeader() { }

   void         *HeaderBegin() { return (void *) &fSize; }
   Int_t         HeaderSize() const { return (Long_t) &fBasicElementSizeType - (Long_t) &fSize + sizeof(fBasicElementSizeType); }
   void          Swap();

   UInt_t        GetEquipmentSize() const { return fSize; }
   UInt_t        GetEquipmentType() const { return fEquipmentType; }
   UInt_t        GetId() const { return fEquipmentID; }
   const UInt_t *GetTypeAttribute() const { return fTypeAttribute; }
   UInt_t        GetBasicSizeType() const { return fBasicElementSizeType; }

private:
   UInt_t fSize;                            // number of raw data bytes
   UInt_t fEquipmentType;                   // equipment type
   UInt_t fEquipmentID;                     // equipment ID
   UInt_t fTypeAttribute[kAttributeWords];  // system (0,1) and user (2) attributes
   UInt_t fBasicElementSizeType;            // basic element size type

   ClassDef(AliRawEquipmentHeader,1) //Alice equipment header
};


class AliRawData : public TObject {

public:
   AliRawData() { fSize = fBufSize = 0; fRawData = 0; fOwner = kTRUE; }
   virtual ~AliRawData() { if (fOwner) delete [] fRawData; }

   inline void SetSize(Int_t size);
   inline void SetBuffer(void *buf, Int_t size);
   Int_t       GetSize() const { return fSize; }
   void       *GetBuffer() { return fRawData; }

private:
   Int_t   fSize;         // number of raw data bytes
   Int_t   fBufSize;      //!actual size of fRawData
   char   *fRawData;      //[fSize] raw event data
   Bool_t  fOwner;        //!if true object owns fRawData buffer

   AliRawData(const AliRawData &);      // not implemented, usage causes
   AliRawData &operator=(const AliRawData &);  // link time error

   ClassDef(AliRawData,1)  // Alice raw event buffer
};

void AliRawData::SetSize(Int_t size)
{
   if (size > fBufSize) {
      if (fOwner) delete [] fRawData;
      fRawData = new char [size];
      fBufSize = size;
      fOwner   = kTRUE;
   }
   fSize = size;
}

void AliRawData::SetBuffer(void *buf, Int_t size)
{
   if (fOwner) delete [] fRawData;
   fRawData = (char *) buf;
   fBufSize = size;
   fSize    = size;
   fOwner   = kFALSE;
}


class AliRawEvent : public TObject {

public:
   AliRawEvent();
   virtual ~AliRawEvent();

   AliRawEventHeader     *GetHeader();
   AliRawEquipmentHeader *GetEquipmentHeader();
   AliRawData            *GetRawData();
   Int_t                  GetNSubEvents() const { return fNSubEvents; }
   AliRawEvent           *NextSubEvent();
   AliRawEvent           *GetSubEvent(Int_t index) const;
   void                   Reset();

private:
   Int_t                  fNSubEvents;  // number of valid sub-events
   AliRawEventHeader     *fEvtHdr;      // event header object
   AliRawEquipmentHeader *fEqpHdr;      // equipment header
   AliRawData            *fRawData;     // raw data container
   TObjArray             *fSubEvents;   // sub AliRawEvent's

   AliRawEvent(const AliRawEvent &);    // not implemented, usage causes
   AliRawEvent &operator=(const AliRawEvent &); // link time error

   ClassDef(AliRawEvent,1)  // ALICE raw event object
};


class AliStats : public TObject {

public:
   AliStats(const char *filename = "", Int_t compmode = 0, Bool_t filter = kFALSE);
   virtual ~AliStats();
   AliStats &operator=(const AliStats &rhs);

   void SetEvents(Int_t events) { fEvents = events; }
   void SetFirstId(Int_t run, Int_t event) { fFirstRun = run; fFirstEvent = event; }
   void SetLastId(Int_t run, Int_t event) { fLastRun = run; fLastEvent = event; }
   void SetBeginTime() { fBegin.Set(); }
   void SetEndTime() { fEnd.Set(); }
   void SetFileSize(Double_t size) { fFileSize = size; }
   void SetCompressionFactor(Float_t comp) { fCompFactor = comp; }
   void Fill(Float_t time);
   void WriteToDB(AliRawDB *rawdb);

   Int_t       GetEvents() const { return fEvents; }
   Int_t       GetFirstRun() const { return fFirstRun; }
   Int_t       GetFirstEvent() const { return fFirstEvent; }
   Int_t       GetLastRun() const { return fLastRun; }
   Int_t       GetLastEvent() const { return fLastEvent; }
   TDatime    &GetBeginTime() { return fBegin; }
   TDatime    &GetEndTime() { return fEnd; }
   Double_t    GetFileSize() const { return fFileSize; }
   Int_t       GetCompressionMode() const { return fCompMode; }
   Float_t     GetCompressionFactor() const { return fCompFactor; }
   Bool_t      GetFilterState() const { return fFilter; }
   const char *GetFileName() const { return fFileName; }
   TH1F       *GetRTHist() const { return fRTHist; }

private:
   Int_t    fEvents;     // number of events in this file
   Int_t    fFirstRun;   // run number of first event in file
   Int_t    fFirstEvent; // event number of first event in file
   Int_t    fLastRun;    // run number of last event in file
   Int_t    fLastEvent;  // event number of last event in file
   TDatime  fBegin;      // begin of filling time
   TDatime  fEnd;        // end of filling time
   TString  fFileName;   // name of file containing this data
   Double_t fFileSize;   // size of file
   Float_t  fCompFactor; // tree compression factor
   Int_t    fCompMode;   // compression mode
   Bool_t   fFilter;     // 3rd level filter on/off
   TH1F    *fRTHist;     // histogram of real-time to process chunck of data
   Float_t  fChunk;      //!chunk to be histogrammed

   AliStats(const AliStats &);  // not implemented

   ClassDef(AliStats,1)  // Statistics object
};


class AliRawDB : public TObject {

public:
   AliRawDB(AliRawEvent *event, Double_t maxsize, Int_t compress,
            Bool_t create = kTRUE);
   virtual ~AliRawDB() { Close(); }

   virtual const char *GetOpenOption() const { return "RECREATE"; }
   virtual Int_t       GetNetopt() const { return 0; }
   virtual Bool_t      Create();
   virtual void        Close();
   void                Fill() { fTree->Fill(); }
   Bool_t              FileFull() { return (fRawDB->GetBytesWritten() > fMaxSize) ?
                                    kTRUE : kFALSE; }

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
   Int_t          fCompress;      // compression mode (1 default)
   Double_t       fMaxSize;       // maximum size in bytes of the raw DB

   virtual const char *GetFileName() const;
   virtual Bool_t      FSHasSpace(const char *fs) const;
   virtual void        MakeTree();

private:
   AliRawDB(const AliRawDB &);       // not implemented, usage causes
   AliRawDB &operator=(const AliRawDB &); // link time error

   ClassDef(AliRawDB,0)  // Raw DB
};


class AliRawRFIODB : public AliRawDB {

public:
   AliRawRFIODB(AliRawEvent *event, Double_t maxsize, Int_t compress);
   ~AliRawRFIODB() { Close(); }

   void Close();

private:
   const char *GetFileName() const;

   ClassDef(AliRawRFIODB,0)  // Raw DB via RFIO
};


class AliRawCastorDB : public AliRawDB {

public:
   AliRawCastorDB(AliRawEvent *event, Double_t maxsize, Int_t compress);
   ~AliRawCastorDB() { Close(); }

   const char *GetOpenOption() const { return "-RECREATE"; }
   Int_t       GetNetopt() const { return 0; }
   void        Close();

private:
   const char *GetFileName() const;

   ClassDef(AliRawCastorDB,0)  // Raw DB via CASTOR and rootd
};


class AliRawRootdDB : public AliRawDB {

public:
   AliRawRootdDB(AliRawEvent *event, Double_t maxsize, Int_t compress);
   ~AliRawRootdDB() { Close(); }

   void Close();

private:
   const char *GetFileName() const;

   ClassDef(AliRawRootdDB,0)  // Raw DB via rootd
};


class AliRawNullDB : public AliRawDB {

public:
   AliRawNullDB(AliRawEvent *event, Double_t maxsize, Int_t compress);
   ~AliRawNullDB() { Close(); }

   void Close();

private:
   const char *GetFileName() const;

   ClassDef(AliRawNullDB,0)  // Raw DB to /dev/null
};


class AliTagDB : public TObject {

public:
   AliTagDB(AliRawEventHeader *header, Double_t maxsize, Bool_t create = kTRUE);
   virtual ~AliTagDB() { Close(); }

   Bool_t          Create();
   virtual void    Close();
   void            Fill() { fTree->Fill(); }
   Bool_t          FileFull()
            { return (fTagDB->GetBytesWritten() > fMaxSize) ? kTRUE : kFALSE; }

   Bool_t          NextFile();

   Double_t           GetBytesWritten() const { return fTagDB->GetBytesWritten(); }
   TFile             *GetDB() const { return fTagDB; }
   const char        *GetDBName() const { return fTagDB->GetName(); }
   AliRawEventHeader *GetHeader() const { return fHeader; }
   Int_t              GetEvents() const { return (Int_t) fTree->GetEntries(); }
   Float_t            GetCompressionFactor() const;

protected:
   TFile             *fTagDB;     // DB to store header information only (tag)
   TTree             *fTree;      // tree use to store header
   AliRawEventHeader *fHeader;    // header via which data is stored
   Double_t           fMaxSize;   // maximum size in bytes of tag DB

   virtual const char *GetFileName() const;

private:
   AliTagDB(const AliTagDB &);       // not implemented, usage causes
   AliTagDB &operator=(const AliTagDB &); // link time error

   ClassDef(AliTagDB,0)  // Tag DB
};


class AliTagNullDB : public AliTagDB {

public:
   AliTagNullDB(AliRawEventHeader *header, Double_t maxsize);
   ~AliTagNullDB() { Close(); }

   void Close();

private:
   const char *GetFileName() const;

   ClassDef(AliTagNullDB,0)   // Tag DB to /dev/null
};


class AliRunDB : public TObject {

public:
   AliRunDB(Bool_t noLocalDB = kFALSE);
   ~AliRunDB() { Close(); }

   void Update(AliStats *stats);
   void UpdateRDBMS(AliStats *stats);
   void UpdateAliEn(AliStats *stats);
   void Close();

private:
   TFile  *fRunDB;     // run database

   AliRunDB(const AliRunDB &);       // not implemented, usage causes
   AliRunDB &operator=(const AliRunDB &); // link time error

   ClassDef(AliRunDB,0)  // Run (bookkeeping) DB
};


class AliMDC : public TObject {

public:
   enum EWriteMode { kLOCAL, kRFIO, kROOTD, kCASTOR, kDEVNULL };

   AliMDC(Int_t fd, Int_t compress, Double_t maxFileSize, Bool_t useFilter,
          EWriteMode mode, Bool_t useLoop, Bool_t delFiles);
   ~AliMDC() { }

   Int_t  Run();
   void   SetStopLoop() { fStopLoop = kTRUE; }
   Bool_t StopLoop() const { return fStopLoop; }

   void   SetDebugLevel(Int_t level) { fDebugLevel = level; }
   Int_t  GetDebugLevel() const { return fDebugLevel; }

   static Bool_t DeleteFiles() { return fgDeleteFiles; }

private:
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

   Int_t     Read(const char *name) { return TObject::Read(name); }
   Int_t     Read(void *buffer, Int_t length);
   Int_t     ReadHeader(AliRawEventHeader &header, void *eb = 0);
   Int_t     ReadEquipmentHeader(AliRawEquipmentHeader &header,
                                 Bool_t isSwapped, void *eb = 0);
   Int_t     ReadRawData(AliRawData &raw, Int_t size, void *eb = 0);
   Int_t     DumpEvent(Int_t toRead);
   Int_t     Filter(AliRawData &raw);

   ClassDef(AliMDC,0)  // MDC processor
};

R__EXTERN AliMDC *gAliMDC;

#define ALIDEBUG(level) \
   if (gAliMDC && (gAliMDC->GetDebugLevel() >= (level)))

#endif
