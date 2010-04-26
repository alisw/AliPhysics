#ifndef ALIITSRAWSTREAMSPD_H
#define ALIITSRAWSTREAMSPD_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

///////////////////////////////////////////////////////////////////////////////
///
/// This class provides access to ITS SPD digits in raw data.
///
///////////////////////////////////////////////////////////////////////////////

#include "AliRawReader.h"
#include "AliITSRawStream.h"
#include "AliITSRawStreamSPDErrorLog.h"

class AliITSRawStreamSPD: public AliITSRawStream {
  public :
    AliITSRawStreamSPD(AliRawReader* rawReader);
    AliITSRawStreamSPD(const AliITSRawStreamSPD& rstream);
    AliITSRawStreamSPD& operator=(const AliITSRawStreamSPD& rstream);
    virtual ~AliITSRawStreamSPD() {};

    virtual Bool_t  Next();
    virtual Int_t   ReadCalibHeader();

    // the 2 methods below are equivalent to AliITSRawStream::GetCoord1 and GetCoord2
    // together with the AliITSRawStream::GetModuleID these are the "offline" coordinates
    Int_t    GetColumn() const {return fCoord1;};
    Int_t    GetRow() const {return fCoord2;};

    // together with the AliRawReader::GetDDLID() these are the "online" coordinates
    UShort_t GetHalfStaveNr() const {return fHalfStaveNr;}
    UShort_t GetChipAddr() const {return fChipAddr;}
    Int_t    GetChipCol() const {return fCol;};
    Int_t    GetChipRow() const {return fRow;};

    // module mapping
    static Int_t GetModuleNumber(UInt_t iDDL, UInt_t iModule);
    static Int_t GetModuleNumber(UInt_t iDDL, UInt_t iHS, UInt_t iChip) 
      {return GetOfflineModuleFromOnline(iDDL,iHS,iChip);}

    // coordinate conversions:
    static Bool_t OfflineToOnline(UInt_t module, UInt_t colM, UInt_t RowM, UInt_t& eq, UInt_t& hs, UInt_t& chip, UInt_t& col, UInt_t& row);
    static Bool_t OnlineToOffline(UInt_t eq, UInt_t hs, UInt_t chip, UInt_t col, UInt_t row, UInt_t& module, UInt_t& colM, UInt_t& rowM);
    // coordinate conversions - offline->online
    static UInt_t GetOnlineEqIdFromOffline(UInt_t module);
    static UInt_t GetOnlineHSFromOffline(UInt_t module);
    static UInt_t GetOnlineChipFromOffline(UInt_t module, UInt_t colM);
    static UInt_t GetOnlineColFromOffline(UInt_t module, UInt_t colM);
    static UInt_t GetOnlineRowFromOffline(UInt_t module, UInt_t rowM);
    // coordinate conversions - online->offline
    static UInt_t GetOfflineModuleFromOnline(UInt_t eqId, UInt_t hs, UInt_t chip);
    static UInt_t GetOfflineChipKeyFromOnline(UInt_t eqId, UInt_t hs, UInt_t chip);
    static UInt_t GetOfflineColFromOnline(UInt_t eqId, UInt_t hs, UInt_t chip, UInt_t col);
    static UInt_t GetOfflineRowFromOnline(UInt_t eqId, UInt_t hs, UInt_t chip, UInt_t row);

    Bool_t  GetFastOrSignal(UInt_t eq, UInt_t hs, UInt_t chip);

    Int_t   GetEventCounter() const {return fEventCounter;}    // get last read event counter value
    Short_t GetEventCounterFullEq(UInt_t eq) const;
    Short_t GetEventCounterFullHS(UInt_t eq, UInt_t hs) const;
    Short_t GetEventCounterFullChip(UInt_t eq, UInt_t hs, UInt_t chip) const;
    Bool_t  IsEventCounterFullConsistent() const;

    Bool_t IsActiveEq(UInt_t eq) const;
    Bool_t IsActiveHS(UInt_t eq, UInt_t hs) const;
    Bool_t IsActiveChip(UInt_t eq, UInt_t hs, UInt_t chip) const;

    Bool_t GetHalfStavePresent(UInt_t hs);

    Int_t  GetHword(UInt_t index);
    // use the methods below to extract the information from the scan calibration header:
    UInt_t GetHrouterNr() const {return (fCalHeadWord[0] & 0x0000003f);}
    Bool_t GetHhalfStaveScanned(UInt_t hs) const;
    UInt_t GetHtype() const {return (Int_t)((fCalHeadWord[1]) & 0x000000ff);}
    Bool_t GetHdataFormat() const {return (Bool_t)(((fCalHeadWord[1]) & 0x00000100)>>8);}
    UInt_t GetHtriggers() const {return fCalHeadWord[2];}
    Bool_t GetHchipPresent(UInt_t hs, UInt_t chip) const;
    UInt_t GetHdacStart() const {return ((fCalHeadWord[5]>>24) & 0x000000ff);}
    UInt_t GetHdacEnd() const {return ((fCalHeadWord[5]>>16) & 0x000000ff);}
    UInt_t GetHdacStep() const {return ((fCalHeadWord[5]>>8) & 0x000000ff);}
    UInt_t GetHdacId() const {return ((fCalHeadWord[5]) & 0x000000ff);}
    UInt_t GetHrowStart() const {return (UInt_t) ((fCalHeadWord[6]>>24) & 0x000000ff);}
    UInt_t GetHrowEnd() const {return (UInt_t) ((fCalHeadWord[6]>>16) & 0x000000ff);}
    UInt_t GetHrowValue() const {return (UInt_t) ((fCalHeadWord[6]>> 8) & 0x000000ff);}
    UInt_t GetHrowSpan() const {return (UInt_t) ((fCalHeadWord[6]) & 0x000000ff);}    
    UInt_t GetHdacValue() const {return (Int_t) ((fCalHeadWord[6]) & 0x000000ff);}
    UInt_t GetHdacHigh(UInt_t hs) const;
    UInt_t GetHdacLow(UInt_t hs) const;
    UInt_t GetHTPAmp(UInt_t hs) const;
    Bool_t GetHminTHchipPresent(UInt_t chip) const;
    UInt_t GetHglobalDBversion() const {return fCalHeadWord[16];}
    // use the methods below to extract the information from the fo calibration header:
    UInt_t GetFOHrouterNr() const {return GetHrouterNr();}
    UInt_t GetFOHtype() const {return GetHtype();}
    UInt_t GetFOHtriggers() const {return GetHtriggers();}
    Bool_t GetFOHchipPresent(UInt_t hs, UInt_t chip) const {return GetHchipPresent(hs,chip);}
    UInt_t GetFOHglobalDBversion() const {return fCalHeadWord[5];}
    UInt_t GetFOHMatrixID() const {return fCalHeadWord[6] & 0x0000000f;}
    UInt_t GetFOHpixelCol() const {return (fCalHeadWord[6] >> 20) & 0x0000001f;}
    UInt_t GetFOHpixelRow() const {return (fCalHeadWord[6] >> 10) & 0x000000ff;}
    UInt_t GetFOHnumDacs() const;
    UInt_t GetFOHdacIndex(UInt_t index) const;
    UInt_t GetFOHdacValue(UInt_t index) const;
    UInt_t GetFOHchipCount(UInt_t hs, UInt_t chip) const;

    void   ActivateAdvancedErrorLog(Bool_t activate, AliITSRawStreamSPDErrorLog* advLogger = NULL);
    AliITSRawStreamSPDErrorLog* GetAdvancedErrorLog() {return fAdvLogger;}

    static const Char_t* GetErrorName(UInt_t errorCode);

    enum {kDDLsNumber = 20};      // number of DDLs in SPD
    enum {kModulesPerDDL = 12};   // number of modules in each DDL
    enum {kCalHeadLenMax = 64};   // maximum number of calib header words
    enum ESPDRawStreamError {
      kTotal,
      kHeaderMissingErr,
      kTrailerMissingErr,
      kTrailerWithoutHeaderErr,
      kHeaderCountErr,
      kTrailerCountErr,
      kFillUnexpectErr,
      kFillMissingErr,
      kWrongFillWordErr,
      kNumberHitsErr,
      kEventCounterErr,
      kDDLNumberErr,
      kHSNumberErr,
      kChipAddrErr,
      kCalHeaderLengthErr,
      kAdvEventCounterErr,     // used by SPDmood
      kAdvEventCounterOrderErr,// used by SPDmood
      kTrailerErrorBitErr,
      kLinkRxDetectorFatalErr,
      kTSMtriggerErr,
      kHighMultiplicityFlag
    };

  private :
    static const Int_t fgkDDLModuleMap[kDDLsNumber][kModulesPerDDL];  // mapping DDL/module -> module number

    Bool_t      ReadNextShort();
    Bool_t      ReadNextInt();
    void        NewEvent();
    void	CheckHeaderAndTrailerCount(Int_t ddlID);

    Int_t       fEventCounter;                // chip event counter
    Short_t     fEventCounterFull[20][6][10]; // chip event counter

    UShort_t    fChipAddr;                    // chip nr
    UShort_t    fHalfStaveNr;                 // half stave nr
    UInt_t      fCol;                         // chip column nr
    UInt_t      fRow;                         // chip row nr
    UInt_t      fCalHeadWord[kCalHeadLenMax]; // calibration header words
    UInt_t      fCalHeadLen;                  // calibration header length

    UShort_t    fData;            // 16 bit data word read
    UInt_t      fOffset;          // offset for cell column
    UInt_t      fHitCount;        // counter of hits
    UChar_t     fDataChar1, fDataChar2, fDataChar3, fDataChar4; // temps part of a 32bit word
    Bool_t      fFirstWord;       // keeps track of which of the two 16bit words out of the 32bit word to read when ReadNextShort is called
    Bool_t      fCalHeadRead[20];            // calibration header read flags (reset at new event)
    UInt_t      fPrevEventId;                // previous event id (needed to know if there is a new event)

    UInt_t      fEqPLBytesRead;              // bytes read for current equipment payload
    UInt_t      fEqPLChipHeadersRead;        // number of chip headers found in current equipment payload
    UInt_t      fEqPLChipTrailersRead;       // number of chip trailers found in current equipment payload
    Bool_t      fHeaderOrTrailerReadLast;    // in previous words, was a header (true) or a trailer (false) read last
    UInt_t      fExpectedHeaderTrailerCount; // expected number of headers and trailers for the current equipment payload
    Bool_t      fFillOutOfSynch;             // has a fill word been put in the wrong place?
    Int_t       fDDLID;                      // ddl id 
    Int_t       fLastDDLID;                  // ddl id for equipment read in previous words

    Bool_t      fAdvancedErrorLog;           // is the advanced error logging activated?
    AliITSRawStreamSPDErrorLog *fAdvLogger;  // pointer to special error logger object

    Bool_t      fFastOrSignal[20][6][10];    // fastor signal bit (one per chip)

    Bool_t      fActiveEq[20];               // which equipments are active (found in data)
    Bool_t      fActiveHS[20][6];            // which half-staves are active (blockattribute bits)
    Bool_t      fActiveChip[20][6][10];      // which chips are active (found in data)

    ClassDef(AliITSRawStreamSPD, 0) // class for reading ITS SPD raw digits
};

#endif
