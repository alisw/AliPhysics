/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//-----------------------------------
//
// decoding of TRD raw data stream 
// and translation into digits
//
//----------------------------------

#ifndef ALITRDRAWSTREAM_H
#define ALITRDRAWSTREAM_H

#include "TObject.h"

#include "AliTRDrawStreamBase.h"
//#include "AliRawReader.h"

class TObjArray;
class TString;
class TTree;
class TBranch;

class AliRawReader;
class AliTRDdigitsManager;
class AliTRDdigitsParam;
class AliTRDarrayADC;
class AliTRDSignalIndex;
class AliTRDtrackletContainer;

class AliTRDrawStream : public AliTRDrawStreamBase
{
 public:
  AliTRDrawStream(AliRawReader *rawReader = 0x0);
  ~AliTRDrawStream();

  Bool_t SetReader(AliRawReader *rawReader) { fRawReader = rawReader; return kTRUE; }
  void SetDigitsManager(AliTRDdigitsManager *digMgr) { fDigitsManager = digMgr; }

  Bool_t ReadEvent(TTree *trackletTree = 0x0);

  Bool_t NextDDL();
  Int_t NextChamber(AliTRDdigitsManager *digMgr, 
		    UInt_t ** /* trackletContainer */, UShort_t ** /* errorContainer */);

  Bool_t ConnectTracklets(TTree *trklTree);

  // legacy code, to be removed
  Bool_t SetRawVersion(Int_t) { return kTRUE; }
  void SetSharedPadReadout(Bool_t) {}
  void SetNoErrorWarning() {}

  // error handling
  enum ErrorCode_t { 
    kUnknown = 0, 
    kLinkMonitor, 
    kPtrgCntMismatch, 
    kNonTrdEq,
    kStackHeaderInvalid,
    kInvalidDetector,
    kNoDigits,
    kHCmismatch,
    kHCcheckFailed,
    kPosUnexp,
    kTPmodeInvalid,
    kTPmismatch,
    kNtimebinsChanged,
    kAdcMaskInconsistent,
    kAdcCheckInvalid,
    kAdcDataAbort,
    kAdcChannelsMiss,
    kMissMcmHeaders,
    kLastErrorCode
  }; 

  enum ErrorBehav_t {
    kTolerate = 0,
    kAbort = 1, 
    kDiscardMCM = 2,
    kDiscardHC = 4
  };

  TTree* GetErrorTree() const { return fErrors; }
  static const char* GetErrorMessage(ErrorCode_t errCode);

  // event statistics
  class AliTRDrawStats : public TObject {
  public:
    AliTRDrawStats() : TObject(), fBytesRead(0) {}
    void ClearStats();

    class AliTRDrawStatsSector : public TObject {
    public:
      AliTRDrawStatsSector() : TObject(), fBytes(0), fBytesRead(0), fNTracklets(0), fNMCMs(0), fNChannels(0) {}
      void ClearStats();

      class AliTRDrawStatsHC : public TObject {
      public:
        AliTRDrawStatsHC() : TObject(), fBytes(0), fBytesRead(0), fNTracklets(0), fNMCMs(0), fNChannels(0) {}
	void ClearStats();

	Int_t fBytes;             // number of bytes (not necessarily read)
 	Int_t fBytesRead;	  // number of bytes read
	Int_t fNTracklets;	  // number of tracklets
	Int_t fNMCMs;		  // number of MCMs (from MCM headers)
	Int_t fNChannels;	  // number of channels
	ClassDef(AliTRDrawStatsHC, 1);
      };

      Int_t fBytes;		     // number of bytes (not necessarily read)
      Int_t fBytesRead;		     // number of bytes read
      Int_t fNTracklets;	     // number of tracklets
      Int_t fNMCMs;		     // number of MCMs (from MCM headers)
      Int_t fNChannels;		     // number of channels
      AliTRDrawStatsHC fStatsHC[60]; //[60] HC-wise statistics
      ClassDef(AliTRDrawStatsSector, 1);
    };


    AliTRDrawStatsSector fStatsSector[18]; //[18] sector-wise statistics
    Int_t fBytesRead;			   // number of bytes read
    ClassDef(AliTRDrawStats, 1);
  };

  AliTRDrawStats fStats; 	     // event statistics, clearing must be done by the user

  AliTRDrawStats* GetStats() { return &fStats; }
  Int_t GetEventSize(Int_t sector) { return fStats.fStatsSector[sector].fBytes; }
  Int_t GetNTracklets(Int_t sector) { return fStats.fStatsSector[sector].fNTracklets; }
  Int_t GetNMCMs(Int_t sector) { return fStats.fStatsSector[sector].fNMCMs; }
  Int_t GetNChannels(Int_t sector) { return fStats.fStatsSector[sector].fNChannels; }

  // raw data dumping
  void SetDumpMCM(Int_t det, Int_t rob, Int_t mcm, Bool_t dump = kTRUE);

  Bool_t IsDumping() { return (fNDumpMCMs > 0); }
  Bool_t DumpingMCM(Int_t det, Int_t rob, Int_t mcm);

  void DumpRaw(TString title, UInt_t *start, Int_t length); 

 protected:
  Int_t ReadSmHeader();
  Int_t ReadStackIndexHeader(Int_t stack);

  Int_t ReadLinkData();
  Int_t ReadTracklets();
  Int_t ReadHcHeader();
  Int_t ReadTPData(Int_t mode = 1);
  Int_t ReadZSData();
  Int_t ReadNonZSData();

  Int_t SeekNextLink();

  // MCM header decoding
  Int_t ROB(UInt_t mcmhdr) const { return 0x7 & mcmhdr >> 28; }
  Int_t MCM(UInt_t mcmhdr) const { return 0xf & mcmhdr >> 24; }
  Int_t Row(UInt_t mcmhdr) const { return (ROB(mcmhdr) / 2) * 4 + MCM(mcmhdr) / 4; }
  Int_t AdcColOffset(UInt_t mcmhdr) const { return (MCM(mcmhdr) % 4 + 1) * 21 + (ROB(mcmhdr) % 2) * 84 - 1; }
  Int_t PadColOffset(UInt_t mcmhdr) const { return (MCM(mcmhdr) % 4 + 1) * 18 + (ROB(mcmhdr) % 2) * 72 + 1; }
  Int_t EvNo(UInt_t mcmhdr) const { return 0xfffff & mcmhdr >> 4; }
  Int_t Check(UInt_t mcmhdr) const { return 0xf & mcmhdr; }
  Int_t CouldBeMCMhdr(UInt_t mcmhdr) const { return ((0xf & mcmhdr) == 0xc); }

  Int_t GetMCMReadoutPos(Int_t mcm) const { return (mcm > -1 && mcm < 16) ? fgkMcmOrder[mcm] : -1; }
  Int_t GetROBReadoutPos(Int_t rob) const { return (rob > -1 && rob < 4) ? fgkRobOrder[rob] : -1; }

  // ADC mask decoding
  Int_t GetActiveChannels(UInt_t adcmask) const { return 0x1fffff & adcmask >> 4; }
  inline Int_t GetNActiveChannelsFromMask(UInt_t adcmask) const; // { Int_t nch = 0; for (Int_t i = 0; i < 21; i++) if ((GetActiveChannels(adcmask) & 1 << i)) nch++; return nch; }
  Int_t GetNActiveChannels(UInt_t adcmask) const { return (0x1f & ~(adcmask >> 25)); }
  Int_t CouldBeADCmask(UInt_t adcmask) const { return ((0xf & adcmask) == 0xc && (0x3 & adcmask >> 30) == 0x1); }
      
  // error message generation
  TString EquipmentError(ErrorCode_t err = kUnknown, TString msg = "", ...); 
  TString StackError    (ErrorCode_t err = kUnknown, TString msg = "", ...);     
  TString LinkError     (ErrorCode_t err = kUnknown, TString msg = "", ...); 
  TString ROBError      (ErrorCode_t err = kUnknown, TString msg = "", ...); 
  TString MCMError      (ErrorCode_t err = kUnknown, TString msg = "", ...); 

  static const char* fgErrorMessages[kLastErrorCode];     // error messages corresponding to the error codes
  static const Int_t fgErrorDebugLevel[kLastErrorCode];   // error debug level
  static       ErrorBehav_t fgErrorBehav[kLastErrorCode]; // bevhaviour in case of error of given type

  // I/O
  AliRawReader *fRawReader;                     // pointer to the raw reader to take the data from 
  AliTRDdigitsManager *fDigitsManager;          // pointer to the digitsManager to fill the data
  AliTRDdigitsParam   *fDigitsParam;            // pointer to the parameters belonging to the digits

  TTree *fErrors;                               // tree containing the occured error codes
  struct { Int_t fSector; Int_t fStack; Int_t fLink; Int_t fError; Int_t fRob; Int_t fMcm; } 
  fLastError;                                   // last error which occured
  UInt_t fErrorFlags;                           // error flags used to steer subsequent reading
  char   fErrorBuffer[100];                     // buffer for error message

  UInt_t *fPayloadStart;                        // pointer to start of data payload
  UInt_t *fPayloadCurr;                         // pointer to current reading position in the payload
  Int_t   fPayloadSize;                         // size of the payload (in UInt_t words)

  static const Int_t fgkNlinks;                 // number of links to read
  static const Int_t fgkNstacks;                // number of stacks to read
  static const UInt_t fgkDataEndmarker;         // data endmarker 
  static const UInt_t fgkTrackletEndmarker;     // tracklet endmarker
  static const Int_t fgkMcmOrder [];            // expected readout order of the MCMs
  static const Int_t fgkRobOrder [];            // expected readout order of the ROBs

  // persistent information
  Int_t  fNtimebins;                            // number of timebins
  Int_t  fLastEvId;                             // Event ID of last event 

  // information valid at current reader position
  // all the variables fCurr... refer to the value at the current
  // reading position
  Int_t fCurrSlot;                              // current slot 
  Int_t fCurrLink;				// current link
  Int_t fCurrRobPos; 				// current ROB number
  Int_t fCurrMcmPos;				// current MCM number

  // DDL header
  UInt_t fCurrEquipmentId;			// current Equipment ID

  // SMU index header
  UInt_t fCurrSmuIndexHeaderSize;		// current size of the SMU index header
  UInt_t fCurrSmuIndexHeaderVersion;		// current version of the SMU index header
  UInt_t fCurrTrackEnable;			// current value of track enable
  UInt_t fCurrTrackletEnable; 			// current value of tracklet enable
  UInt_t fCurrStackMask;			// current mask of active stacks

  // Stack index header
  UInt_t *fCurrStackIndexWord;			// current stack index words
  UInt_t *fCurrStackHeaderSize;			// current stack index sizes
  UInt_t *fCurrStackHeaderVersion;		// current stack header versions
  UInt_t *fCurrLinkMask;			// current link masks
  UInt_t *fCurrCleanCheckout;			// current clean checkout flags
  UInt_t *fCurrBoardId;				// current board IDs
  UInt_t *fCurrHwRev;				// current hardware revision
  UInt_t *fCurrLinkMonitorFlags;		// current link monitor flags
  UInt_t *fCurrLinkDataTypeFlags;		// current link data flags
  UInt_t *fCurrLinkDebugFlags;			// current link debug flags

  // HC information
  Int_t fCurrSpecial;				// current value of the special flag
  Int_t fCurrMajor;				// current major version
  Int_t fCurrMinor;				// current minor version
  Int_t fCurrAddHcWords;			// current number of additional HC-header words
  Int_t fCurrSm;				// current sector
  Int_t fCurrStack;				// current stack
  Int_t fCurrLayer;				// current layer
  Int_t fCurrSide;				// current side
  Int_t fCurrHC;				// current HC
  Int_t fCurrCheck;				// current check bits
  Int_t fCurrNtimebins;				// current number of timebins
  Int_t fCurrBC;				// current BC
  Int_t fCurrPtrgCnt;				// current pretrigger count
  Int_t fCurrPtrgPhase;				// current pretrigger phase

  // settings for dumping
  Int_t fDumpMCM[100];		                // MCMs to dump
  Int_t fNDumpMCMs;                             // number of MCMs to dump

  // tracklet information
  TClonesArray *fTrackletArray;			// pointer to array for tracklet storage

  // output data
  AliTRDarrayADC *fAdcArray;			// pointer to ADC array
  AliTRDSignalIndex *fSignalIndex;		// pointer to the signal index
  TTree *fTrackletTree; 			// pointer to the tree for tracklet storage

  AliTRDrawStream(const AliTRDrawStream&);           // not implemented
  AliTRDrawStream& operator=(const AliTRDrawStream&); // not implemented

  ClassDef(AliTRDrawStream, 0);
};

Int_t AliTRDrawStream::GetNActiveChannelsFromMask(UInt_t adcmask) const
{
  // return number of active bits in the ADC mask

  adcmask = GetActiveChannels(adcmask);
  adcmask = adcmask - ((adcmask >> 1) & 0x55555555);
  adcmask = (adcmask & 0x33333333) + ((adcmask >> 2) & 0x33333333);
  return (((adcmask + (adcmask >> 4)) & 0xF0F0F0F) * 0x1010101) >> 24;
}

#endif
