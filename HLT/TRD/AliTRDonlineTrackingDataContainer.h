#ifndef ALITRDONLINETRACKINGDATACONTAINER_H
#define ALITRDONLINETRACKINGDATACONTAINER_H

#include "AliESDTrdTracklet.h"
#include "AliESDTrdTrack.h"
#include "AliHLTLogging.h"

class AliTRDonlineTrackingDataContainer : public AliHLTLogging {
 public:
  AliTRDonlineTrackingDataContainer();
  AliTRDonlineTrackingDataContainer(const AliTRDonlineTrackingDataContainer& cont);
  ~AliTRDonlineTrackingDataContainer();

  void SetLogPrefix(const char* prefix) { *fLogPrefix = prefix; };
  void SetGtuPtMultiplier(Double_t mult) { fGtuPtMultiplier = mult; }
  void SetGtuPtMultiplierFromMagField(Double_t magField);

  void Clear(const Option_t* = "");
  Int_t AddTracklet(UInt_t HCId, UInt_t trackletWord);
  Int_t AddTracklet(const AliESDTrdTracklet* tracklet);
  Int_t AddTrack(UShort_t stack, ULong64_t trackWord, ULong64_t extTrackWord, const UInt_t trackletWords[6], Int_t addInfo = -1);
  Int_t AddTrack(const AliESDTrdTrack* track, Int_t addInfo = -1);

  void SetTrackAddInfo(UShort_t stack, UInt_t trackIndex, Int_t addInfo);
  void SetSectorTrgWord(UShort_t sector, UInt_t trgWord) { fSectorTrgWords[sector] = trgWord; };
  void SetStackTrgWord(UShort_t sector, UShort_t stack, ULong64_t trgWord) { fStackTrgWords[sector][stack] = trgWord; };

  Int_t GetNumTracklets();
  Int_t GetNumTracklets(UInt_t det);
  Int_t GetNumTracks();
  Int_t GetNumTracks(UShort_t stack);

  Int_t GetTrackletHCId(UInt_t det, UInt_t trackletIndex) { return fTrackletHCId[det][trackletIndex]; };
  Int_t GetTrackletBinY(UInt_t det, UInt_t trackletIndex);
  Int_t GetTrackletBinDy(UInt_t det, UInt_t trackletIndex);
  Int_t GetTrackletBinZ(UInt_t det, UInt_t trackletIndex);
  Int_t GetTrackletPID(UInt_t det, UInt_t trackletIndex);
  Float_t GetTrackletLocalY(UInt_t det, UInt_t trackletIndex);

  AliESDTrdTracklet* GetTracklet(UInt_t det, UInt_t trackletIndex);
  AliESDTrdTrack* GetTrack(UInt_t stack, UInt_t trackIndex, Bool_t constructTracklets = kTRUE);
  Int_t GetTrackAddInfo(UShort_t stack, UInt_t trackIndex);

  inline Int_t GetTrackLayerMask(UInt_t stack, UInt_t trackIndex){
    return (fTrackWords[stack][trackIndex] >> 56) & 0x3f;
  }

  inline Int_t GetTrackA(UInt_t stack, UInt_t trackIndex){
    return (((fTrackWords[stack][trackIndex] >> 38) & 0x3ffff) ^ 0x20000) - 0x20000;
  }

  inline Int_t GetTrackPID(UInt_t stack, UInt_t trackIndex){
    return fTrackWords[stack][trackIndex] & 0xff;
  }

  Double_t GetTrackPt(UInt_t stack, UInt_t trackIndex);
  UShort_t GetTrackLayerNum(UInt_t stack, UInt_t trackIndex);
  UInt_t GetTrackTrackletWord(UInt_t stack, UInt_t trackIndex, UShort_t layer);
  Int_t GetTrackTrackletBinY(UInt_t stack, UInt_t trackIndex, UShort_t layer);
  Float_t GetTrackTrackletLocalY(UInt_t stack, UInt_t trackIndex, UShort_t layer);
  Int_t GetTrackTrackletBinZ(UInt_t stack, UInt_t trackIndex, UShort_t layer);
  UShort_t GetTrackTrackletPID(UInt_t stack, UInt_t trackIndex, UShort_t layer);

  UInt_t GetSectorTrgWord(UShort_t sector) { return fSectorTrgWords[sector]; };
  ULong64_t GetStackTrgWord(UShort_t sector, UShort_t stack) { return fStackTrgWords[sector][stack]; };
  UInt_t GetSectorTrgContribs(UShort_t sector) { return fSectorTrgWords[sector] & 0x3ff; };
  Float_t GetTrackletStartTime(UShort_t sector, UShort_t stack) { return ((fStackTrgWords[sector][stack] >> 20) & 0x3ff)*1./125.; };
  Float_t GetTrackletEndTime(UShort_t sector, UShort_t stack) { return ((fStackTrgWords[sector][stack] >> 10) & 0x3ff)*1./125.; };
  Float_t GetTMUTrackingDoneTime(UShort_t sector, UShort_t stack) { return ((fStackTrgWords[sector][stack] >> 0) & 0x3ff)*1./50.; };
  Float_t GetSMUTrackingDoneTime(UShort_t sector) { return ((fSectorTrgWords[sector] >> 12) & 0x3ff)*1./120.; };

  Bool_t Compress(void* &buffer, UInt_t &sizeInBytes);
  Bool_t Decompress(const void* buffer, UInt_t sizeInBytes, Bool_t cumulative = kFALSE, Bool_t verbose = kFALSE);

  void PrintBuffer(const void* buffer, UInt_t sizeInBytes, const char* identifier);
  void PrintSummary(const char* identifier = "none");

 protected:

  AliTRDonlineTrackingDataContainer& operator=(const AliTRDonlineTrackingDataContainer &rhs); // not implemented

  static const unsigned int fgkLeadingMagicConst = 0x1234face;    //! header magic constant
  static const unsigned int fgkHeaderTypeData = 0xe;              //! header type constant
  static const unsigned int fgkHeaderTypeStack = 0xa;             //! header type constant
  static const unsigned int fgkHeaderTypeTracklet = 0xb;          //! header type constant
  static const unsigned int fgkHeaderTypeTrack = 0xc;             //! header type constant
  static const unsigned int fgkHeaderTypeGtu = 0xd;               //! header type constant
  static const unsigned int fgkCrosscheckSeed = 0xc5f3a19b;       //! crosscheck seed constant
  static const unsigned int fgkInvalidTrackletWord = 0x10001000;  //! marking invalid tracklet word
  static const unsigned int fgkTrailingMagicConst = 0x5678beaf;   //! trailer magic constant

  static const unsigned int fgkNumLayers = 6;                     //! number of TRD layers
  static const unsigned int fgkNumStacks = 90;                    //! number of TRD stacks
  static const unsigned int fgkNumChambers = 540;                 //! number of TRD chambers
  static const unsigned int fgkNumStacksPerSector = 5;            //! number of TRD stacks per sector
  static const unsigned int fgkNumSectors = 18;                   //! number of TRD sectors
  static const unsigned int fgkMaxTrackletsPerChamber = 1024;     //! maximum number of tracklets per chamber
  static const unsigned int fgkMaxTracksPerStack = 32;            //! maximum number of tracks per stack
  static const Float_t fgkBinWidthY;                              //! geometric TRD constant

  Double_t fGtuPtMultiplier;                                         //! factor for sign correction
  Bool_t fStoreGtuInfo;                                              //! controls storage of GTU header info upon compress
  TString* fLogPrefix;                                               //! log message prefix

  UInt_t fNumTracklets[fgkNumChambers];                              //! number of tracklets by chamber
  UInt_t fTrackletWords[fgkNumChambers][fgkMaxTrackletsPerChamber];  //! tracklet words by chamber
  UInt_t fTrackletHCId[fgkNumChambers][fgkMaxTrackletsPerChamber];   //! half-chambeer to tracklets

  UInt_t fNumTracks[fgkNumStacks];                                               //! number of tracks by stack
  ULong64_t fTrackWords[fgkNumStacks][fgkMaxTracksPerStack];                     //! track words by stack
  ULong64_t fTrackExtWords[fgkNumStacks][fgkMaxTracksPerStack];                  //! extended track words by stack
  UInt_t fTrackTrackletWords[fgkNumStacks][fgkMaxTracksPerStack][fgkNumLayers];  //! tracklet words by track and stack
  Int_t fTrackAddInfo[fgkNumStacks][fgkMaxTracksPerStack];                       //! additional info to tracks by stack

  UInt_t fSectorTrgWords[fgkNumSectors];                                         //! sector trigger words
  ULong64_t fStackTrgWords[fgkNumSectors][fgkNumStacksPerSector];                //! stack trigger words (=TMU tracking done words)

  UInt_t DataWordsNeeded();
  UInt_t MakeDataHeader();
  UInt_t MakeStackHeader(UShort_t stack, Bool_t trackletsPresent,
			 Bool_t tracksPresent, UInt_t size);
  UInt_t MakeTrackletHeader(UShort_t det);
  UInt_t MakeTrackHeader(UShort_t stack);
  UInt_t MakeGtuHeader(Bool_t storeInfo);

  Int_t StackFromStackHeader(UInt_t stackHeader);
  UInt_t SizeFromStackHeader(UInt_t stackHeader);
  Int_t HCIdFromTrackletHeader(UInt_t trackletHeader);
  Bool_t StoreGtuInfoFromDataHeader(UInt_t dataHeader);
  Bool_t GtuInfoPresentFromGtuHeader(UInt_t gtuHeader);

  UShort_t ExtractTrackLayerMask(ULong64_t trackWord);
  UShort_t ExtractTrackPID(ULong64_t trackWord);

  ClassDef(AliTRDonlineTrackingDataContainer, 1);
};

#endif
