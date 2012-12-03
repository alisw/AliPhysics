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
  void SetGtuPtMultiplier(const Double_t mult) { fGtuPtMultiplier = mult; }
  void SetGtuPtMultiplierFromMagField(const Double_t magField);

  void Clear(const Option_t* = "");
  Int_t AddTracklet(const UInt_t HCId, const UInt_t trackletWord);
  Int_t AddTracklet(const AliESDTrdTracklet* tracklet);
  Int_t AddTrack(const UShort_t stack, const ULong64_t trackWord, const ULong64_t extTrackWord, const UInt_t trackletWords[6], const Int_t addInfo = -1);
  Int_t AddTrack(const AliESDTrdTrack* track, const Int_t addInfo = -1);

  void SetTrackAddInfo(const UShort_t stack, const UInt_t trackIndex, const Int_t addInfo);
  void SetSectorTrgWord(const UShort_t sector, const UInt_t trgWord) { fSectorTrgWords[sector] = trgWord; };
  void SetStackTrgWord(const UShort_t sector, const UShort_t stack, const ULong64_t trgWord) { fStackTrgWords[sector][stack] = trgWord; };

  Int_t GetNumTracklets();
  Int_t GetNumTracklets(const UInt_t det);
  Int_t GetNumTracks();
  Int_t GetNumTracks(const UShort_t stack);

  Int_t GetTrackletHCId(const UInt_t det, const UInt_t trackletIndex) { return fTrackletHCId[det][trackletIndex]; };
  Int_t GetTrackletBinY(const UInt_t det, const UInt_t trackletIndex);
  Int_t GetTrackletBinDy(const UInt_t det, const UInt_t trackletIndex);
  Int_t GetTrackletBinZ(const UInt_t det, const UInt_t trackletIndex);
  Int_t GetTrackletPID(const UInt_t det, const UInt_t trackletIndex);
  Float_t GetTrackletLocalY(const UInt_t det, const UInt_t trackletIndex);

  AliESDTrdTracklet* GetTracklet(const UInt_t det, const UInt_t trackletIndex);
  AliESDTrdTrack* GetTrack(const UInt_t stack, const UInt_t trackIndex, const Bool_t constructTracklets = kTRUE);
  Int_t GetTrackAddInfo(const UShort_t stack, const UInt_t trackIndex);

  inline Int_t GetTrackLayerMask(const UInt_t stack, const UInt_t trackIndex){
    return (fTrackWords[stack][trackIndex] >> 56) & 0x3f;
  }

  inline Int_t GetTrackA(const UInt_t stack, const UInt_t trackIndex){
    return (((fTrackWords[stack][trackIndex] >> 38) & 0x3ffff) ^ 0x20000) - 0x20000;
  }

  inline Int_t GetTrackPID(const UInt_t stack, const UInt_t trackIndex){
    return fTrackWords[stack][trackIndex] & 0xff;
  }

  Double_t GetTrackPt(const UInt_t stack, const UInt_t trackIndex);
  UShort_t GetTrackLayerNum(const UInt_t stack, const UInt_t trackIndex);
  UInt_t GetTrackTrackletWord(const UInt_t stack, const UInt_t trackIndex, const UShort_t layer);
  Int_t GetTrackTrackletBinY(const UInt_t stack, const UInt_t trackIndex, const UShort_t layer);
  Float_t GetTrackTrackletLocalY(const UInt_t stack, const UInt_t trackIndex, const UShort_t layer);
  Int_t GetTrackTrackletBinZ(const UInt_t stack, const UInt_t trackIndex, const UShort_t layer);
  UShort_t GetTrackTrackletPID(const UInt_t stack, const UInt_t trackIndex, const UShort_t layer);

  UInt_t GetSectorTrgWord(const UShort_t sector) { return fSectorTrgWords[sector]; };
  ULong64_t GetStackTrgWord(const UShort_t sector, const UShort_t stack) { return fStackTrgWords[sector][stack]; };
  UInt_t GetSectorTrgContribs(const UShort_t sector) { return fSectorTrgWords[sector] & 0x3ff; };
  Float_t GetTrackletStartTime(const UShort_t sector, const UShort_t stack) { return ((fStackTrgWords[sector][stack] >> 20) & 0x3ff)*1./125.; };
  Float_t GetTrackletEndTime(const UShort_t sector, const UShort_t stack) { return ((fStackTrgWords[sector][stack] >> 10) & 0x3ff)*1./125.; };
  Float_t GetTMUTrackingDoneTime(const UShort_t sector, const UShort_t stack) { return ((fStackTrgWords[sector][stack] >> 0) & 0x3ff)*1./50.; };
  Float_t GetSMUTrackingDoneTime(const UShort_t sector) { return ((fSectorTrgWords[sector] >> 12) & 0x3ff)*1./120.; };

  Bool_t Compress(void* &buffer, UInt_t &sizeInBytes);
  Bool_t Decompress(const void* buffer, const UInt_t sizeInBytes, const Bool_t cumulative = kFALSE, const Bool_t verbose = kFALSE);

  void PrintBuffer(const void* buffer, const UInt_t sizeInBytes, const char* identifier);
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
  UInt_t MakeStackHeader(const UShort_t stack, const Bool_t trackletsPresent,
			 const Bool_t tracksPresent, const UInt_t size);
  UInt_t MakeTrackletHeader(const UShort_t det);
  UInt_t MakeTrackHeader(const UShort_t stack);
  UInt_t MakeGtuHeader(const Bool_t storeInfo);

  Int_t StackFromStackHeader(const UInt_t stackHeader);
  UInt_t SizeFromStackHeader(const UInt_t stackHeader);
  Int_t HCIdFromTrackletHeader(const UInt_t trackletHeader);
  Bool_t StoreGtuInfoFromDataHeader(const UInt_t dataHeader);
  Bool_t GtuInfoPresentFromGtuHeader(const UInt_t gtuHeader);

  UShort_t ExtractTrackLayerMask(const ULong64_t trackWord);
  UShort_t ExtractTrackPID(const ULong64_t trackWord);

  ClassDef(AliTRDonlineTrackingDataContainer, 1);
};

#endif
