#ifndef ALITRDRAWSTREAM_H
#define ALITRDRAWSTREAM_H

#include "TObject.h"

#include "AliTRDrawStreamBase.h"
//#include "AliRawReader.h"

class TObjArray;
class TString;
class TTree;

class AliRawReader;
class AliTRDdigitsManager;
class AliTRDdigitsParam;
class AliTRDarrayADC;
class AliTRDSignalIndex;

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

  TTree* GetErrorTree() { return fErrors; }
  static const char* GetErrorMessage(ErrorCode_t errCode);

 protected:
  Int_t ReadSmHeader();
  Int_t ReadStackIndexHeader(Int_t stack);

  Int_t ReadLinkData();
  Int_t ReadTracklets();
  Int_t ReadHcHeader();
  Int_t ReadTPData(Int_t mode = 1);
  Int_t ReadZSData();
  Int_t ReadNonZSData();

  // MCM header decoding
  Int_t ROB(UInt_t mcmhdr) { return 0x7 & mcmhdr >> 28; }
  Int_t MCM(UInt_t mcmhdr) { return 0xf & mcmhdr >> 24; }
  Int_t Row(UInt_t mcmhdr) { return (ROB(mcmhdr) / 2) * 4 + MCM(mcmhdr) / 4; }
  Int_t AdcColOffset(UInt_t mcmhdr) { return (MCM(mcmhdr) % 4 + 1) * 21 + (ROB(mcmhdr) % 2) * 84 - 1; }
  Int_t PadColOffset(UInt_t mcmhdr) { return (MCM(mcmhdr) % 4 + 1) * 18 + (ROB(mcmhdr) % 2) * 72 + 1; }
  Int_t EvNo(UInt_t mcmhdr) { return 0xfffff & mcmhdr >> 4; }
  Int_t Check(UInt_t mcmhdr) { return 0xf & mcmhdr; }
  Int_t CouldBeMCMhdr(UInt_t mcmhdr) { return ((0xf & mcmhdr) == 0xc); }

  Int_t GetMCMReadoutPos(Int_t mcm) { return (mcm > -1 && mcm < 16) ? fgkMcmOrder[mcm] : -1; }
  Int_t GetROBReadoutPos(Int_t rob) { return (rob > -1 && rob < 4) ? fgkRobOrder[rob] : -1; }

  // ADC mask decoding
  Int_t GetActiveChannels(UInt_t adcmask) { return 0x1fffff & adcmask >> 4; }
  Int_t GetNActiveChannelsFromMask(UInt_t adcmask); // { Int_t nch = 0; for (Int_t i = 0; i < 21; i++) if ((GetActiveChannels(adcmask) & 1 << i)) nch++; return nch; }
  Int_t GetNActiveChannels(UInt_t adcmask) { return (0x1f & ~(adcmask >> 25)); }
  Int_t CouldBeADCmask(UInt_t adcmask) { return ((0xf & adcmask) == 0xc && (0x3 & adcmask >> 30) == 0x1); }
      
  // error message generation
  TString EquipmentError(ErrorCode_t err = kUnknown, TString msg = ""); 
  TString StackError    (ErrorCode_t err = kUnknown, TString msg = "");     
  TString LinkError     (ErrorCode_t err = kUnknown, TString msg = ""); 
  TString ROBError      (ErrorCode_t err = kUnknown, TString msg = ""); 
  TString MCMError      (ErrorCode_t err = kUnknown, TString msg = ""); 

  static char* fgErrorMessages[kLastErrorCode];

  // I/O
  AliRawReader *fRawReader;
  AliTRDdigitsManager *fDigitsManager;
  AliTRDdigitsParam   *fDigitsParam;

  TTree *fErrors; 
  struct { Int_t fSector; Int_t fStack; Int_t fLink; Int_t fError; Int_t fRob; Int_t fMcm; } fLastError; 

  UInt_t *fPayloadStart;
  UInt_t *fPayloadCurr;
  Int_t   fPayloadSize;

  static const Int_t fgkNlinks;
  static const Int_t fgkNstacks;
  static const UInt_t fgkDataEndmarker;
  static const UInt_t fgkTrackletEndmarker;
  static const Int_t fgkMcmOrder [];
  static const Int_t fgkRobOrder [];

  // persistent information
  Int_t  fNtimebins;
  Int_t  fLastEvId;

  // information valid at current reader position
  Int_t fCurrSlot;
  Int_t fCurrLink;
  Int_t fCurrRobPos; 
  Int_t fCurrMcmPos;

  // DDL header
  UInt_t fCurrEquipmentId;

  // SMU index header
  UInt_t fCurrSmuIndexHeaderSize;
  UInt_t fCurrSmuIndexHeaderVersion;
  UInt_t fCurrTrackEnable;
  UInt_t fCurrTrackletEnable; 
  UInt_t fCurrStackMask;

  // Stack index header
  UInt_t *fCurrStackIndexWord;
  UInt_t *fCurrStackHeaderSize;
  UInt_t *fCurrStackHeaderVersion;
  UInt_t *fCurrLinkMask;
  UInt_t *fCurrCleanCheckout;
  UInt_t *fCurrBoardId;
  UInt_t *fCurrHwRev;
  UInt_t *fCurrLinkMonitorFlags;
  UInt_t *fCurrLinkDataTypeFlags;
  UInt_t *fCurrLinkDebugFlags;

  // HC information
  Int_t fCurrSpecial;
  Int_t fCurrMajor;
  Int_t fCurrMinor;
  Int_t fCurrAddHcWords;
  Int_t fCurrSm;
  Int_t fCurrStack;
  Int_t fCurrLayer;
  Int_t fCurrSide;
  Int_t fCurrHC;
  Int_t fCurrCheck;
  Int_t fCurrNtimebins;
  Int_t fCurrBC;
  Int_t fCurrPtrgCnt;
  Int_t fCurrPtrgPhase;

  // tracklet information
  TClonesArray *fTrackletArray;

  // output data
  AliTRDarrayADC *fAdcArray;
  AliTRDSignalIndex *fSignalIndex;
  TTree *fTrackletTree; 

  AliTRDrawStream(const AliTRDrawStream&);           // not implemented
  AliTRDrawStream operator=(const AliTRDrawStream&); // not implemented

  ClassDef(AliTRDrawStream, 0);
};

#endif
