#ifndef ALITOFRAWSTREAM_H
#define ALITOFRAWSTREAM_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

///////////////////////////////////////////////////////////////
//                                                           //
//   This class provides the key-reading for TOF raw data.   //
//                                                           //
///////////////////////////////////////////////////////////////

#include "TObject.h"

//#include "AliTOFHitData.h"
#include "AliTOFHitDataBuffer.h"
#include "AliTOFDecoder.h"
//#include "AliTOFCableLengthMap.h"

class AliTOFHitData;

class AliTOFDecoderV2;

/**********************************
 * OLD DEFINITIONS 
 **********************************/ 

/******************************************
GENERAL DATA FORMAT                    
******************************************/

//filler
//#ifndef FILLER
#define FILLER 0x70000000
//#endif

//word type mask/position
#define WORD_TYPE_MASK 0xf0000000
#define WORD_TYPE_POSITION 28

//global header word required bit pattern
#define GLOBAL_HEADER 0x40000000

//global trailer word required bit pattern
#define GLOBAL_TRAILER 0x50000000

//error word required bit pattern
// already defined in AliTOFDecoder
//#define ERROR 0x30000000

//header slot ID mask/position
#define HEADER_SLOT_ID_MASK 0x0000000f
#define HEADER_SLOT_ID_POSITION 0

//word types
#define GLOBAL_HEADER_TYPE 4
#define GLOBAL_TRAILER_TYPE 5
#define ERROR_TYPE 6
#define FILLER_TYPE 7
#define TRM_CHAIN0_HEADER_TYPE 0
#define TRM_CHAIN0_TRAILER_TYPE 1
#define TRM_CHAIN1_HEADER_TYPE 2
#define TRM_CHAIN1_TRAILER_TYPE 3

//slot types
#define DRM_ID_NUMBER 1
#define LTM_ID_NUMBER 2


/******************************************
DRM DATA FORMAT                        
******************************************/

//DRM global header word required bit pattern
#define DRM_GLOBAL_HEADER 0x40000001

//DRM event words mask/position
#define DRM_EVENT_WORDS_MASK 0x001ffff0
#define DRM_EVENT_WORDS_POSITION 4

//DRM DRM ID mask/position
#define DRM_DRM_ID_MASK 0x0fe00000
#define DRM_DRM_ID_POSITION 21

//DRM status header 1 word required bit pattern
#define DRM_STATUS_HEADER_1 0x40000001

//DRM slot ID mask/position
#define DRM_SLOT_ID_MASK 0x00007ff0
#define DRM_SLOT_ID_POSITION 4

//DRM C-bit mask/position
#define DRM_C_BIT_MASK 0x00008000
#define DRM_C_BIT_POSITION 15

//DRM Vers-ID mask/position
#define DRM_VERS_ID_MASK 0x001f0000
#define DRM_VERS_ID_POSITION 16

//DRM DRM Header size mask/position
#define DRM_HEADER_SIZE_MASK 0x01e00000
#define DRM_HEADER_SIZE_POSITION 21

//DRM status header 2 word required bit pattern
#define DRM_STATUS_HEADER_2 0x40000001

//DRM enable ID mask/position
#define DRM_ENABLE_ID_MASK 0x00007ff0
#define DRM_ENABLE_ID_POSITION 4

//DRM zero in word2 mask/position
#define DRM_ZERO_WORD2_MASK 0x00008000
#define DRM_ZERO_WORD2_POSITION 15

//DRM fault ID mask/position
#define DRM_FAULT_ID_MASK 0x07ff0000
#define DRM_FAULT_ID_POSITION 16

//DRM RTO bit mask/position
#define DRM_RTO_BIT_MASK 0x08000000
#define DRM_RTO_BIT_POSITION 27

//DRM status header 3 word required bit pattern
#define DRM_STATUS_HEADER_3 0x40000001

//DRM L0 BCID mask/position
#define DRM_L0_BCID_MASK 0x0000fff0
#define DRM_L0_BCID_POSITION 4

//DRM Run Time Info mask/position
#define DRM_RUNTIME_INFO_MASK 0x0fff0000
#define DRM_RUNTIME_INFO_POSITION 16

//DRM status header 4 word required bit pattern
#define DRM_STATUS_HEADER_4 0x40000001

//DRM Temperature mask/position
#define DRM_TEMPERATURE_MASK 0x00003ff0
#define DRM_TEMPERATURE_POSITION 4

//DRM 1st zero in word4 mask/position
#define DRM_ZERO_1_WORD4_MASK 0x00004000
#define DRM_ZERO_1_WORD4_POSITION 14

//DRM ACK mask/position
#define DRM_ACK_MASK 0x00008000
#define DRM_ACK_POSITION 15

//DRM Sens AD mask/position
#define DRM_SENS_AD_MASK 0x00070000
#define DRM_SENS_AD_POSITION 16

//DRM 2nd zero in word4 mask/position
#define DRM_ZERO_2_WORD4_MASK 0x00080000
#define DRM_ZERO_2_WORD4_POSITION 19

//DRM event CRC mask/position
#define DRM_EVENT_CRC_MASK 0x000ffff0
#define DRM_EVENT_CRC_POSITION 4

//DRM global trailer word required bit pattern
#define DRM_GLOBAL_TRAILER 0x50000001

//DRM local event counter mask/position
#define DRM_LOCAL_EVENT_COUNTER_MASK 0x0000fff0
#define DRM_LOCAL_EVENT_COUNTER_POSITION 4


/******************************************
TRM DATA FORMAT                        
******************************************/

//TRM global header word required bit pattern
#define TRM_GLOBAL_HEADER 0x40000000

//TRM slot ID mask/position
#define TRM_SLOT_ID_MASK 0x0000000f
#define TRM_SLOT_ID_POSITION 0

//TRM event words mask/position
#define TRM_EVENT_WORDS_MASK 0x0001fff0
#define TRM_EVENT_WORDS_POSITION 4

//TRM ACQ-bits mask/position
#define TRM_ACQ_BITS_MASK 0x00060000
#define TRM_ACQ_BITS_POSITION 17

//TRM L-bit mask/position
#define TRM_L_BIT_MASK 0x00080000
#define TRM_L_BIT_POSITION 19

//TRM chain-0 header word required bit pattern
#define TRM_CHAIN_0_HEADER 0x00000000

//TRM chain-1 header word required bit pattern
#define TRM_CHAIN_1_HEADER 0x20000000

//TRM bunch ID mask/position
#define TRM_BUNCH_ID_MASK 0x0000fff0
#define TRM_BUNCH_ID_POSITION 4

//TRM PB24 temp mask/position
#define TRM_PB24_TEMP_MASK 0x00ff0000
#define TRM_PB24_TEMP_POSITION 16

//TRM PB24 ID mask/position
#define TRM_PB24_ID_MASK 0x07000000
#define TRM_PB24_ID_POSITION 24

//TRM TS-bit mask/position
#define TRM_TS_BIT_MASK 0x08000000
#define TRM_TS_BIT_POSITION 27

//TRM chain-0 trailer word required bit pattern
#define TRM_CHAIN_0_TRAILER 0x10000000

//TRM chain-1 trailer word required bit pattern
#define TRM_CHAIN_1_TRAILER 0x30000000

//TRM status mask/position
#define TRM_STATUS_MASK 0x0000000f
#define TRM_STATUS_POSITION 0


//TDC digit

//TRM TDC digit word required bit pattern
#define TRM_TDC_DIGIT 0x8000000

//TRM digit time mask/position
#define TRM_DIGIT_TIME_MASK 0x00001fff
#define TRM_DIGIT_TIME_POSITION 0

//TRM long digit time mask/position
#define TRM_LONG_DIGIT_TIME_MASK 0x001fffff
#define TRM_LONG_DIGIT_TIME_POSITION 0

//TRM TOT width mask/position
#define TRM_TOT_WIDTH_MASK 0x001fe000
#define TRM_TOT_WIDTH_POSITION 13

//TRM chan mask/position
#define TRM_CHAN_MASK 0x00e00000
#define TRM_CHAN_POSITION 21

//TRM TDC ID mask/position
#define TRM_TDC_ID_MASK 0x0f000000
#define TRM_TDC_ID_POSITION 24

//TRM E-bit mask/position
#define TRM_E_BIT_MASK 0x10000000
#define TRM_E_BIT_POSITION 28

//TRM PS-bits mask/position
#define TRM_PS_BITS_MASK 0x60000000
#define TRM_PS_BITS_POSITION 29

#define TRM_FIRST_SLOT_ID 3

//define hptdc time bin width
#define TIME_BIN_WIDTH 24.4e-3 //ns

//define hptdc tot bin width
// already defined in AliTOFDecoder
//#define TOT_BIN_WIDTH 48.4e-3 //ns

//TRM errors

//TRM TDC error word required bit pattern
#define TRM_TDC_ERROR 0x6000000

//TRM TDC diagnostic error word required bit pattern
#define TRM_TDC_DIAGNOSTIC_ERROR 0x6f00000

//TRM TDC error flags mask/position
#define TRM_TDC_ERROR_FLAGS_MASK 0x00007fff
#define TRM_TDC_ERROR_FLAGS_POSITION 0

//TRM TDC error TDC ID mask/position
#define TRM_TDC_ERROR_TDC_ID_MASK 0x0f00000
#define TRM_TDC_ERROR_TDC_ID_POSITION 24

//TRM TDC fault chip flag ID mask/position
#define TRM_TDC_ERROR_FAULT_CHIP_FLAG_ID_MASK 0x00007fff
#define TRM_TDC_ERROR_FAULT_CHIP_FLAG_ID_POSITION 0

//TRM TDC error C-bit mask/position
#define TRM_TDC_ERROR_C_BIT_MASK 0x00008000
#define TRM_TDC_ERROR_C_BIT_POSITION 15

//TRM TDC JTAG error code mask/position
#define TRM_TDC_ERROR_JTAG_ERROR_CODE_MASK 0x000007ff
#define TRM_TDC_ERROR_JTAG_ERROR_CODE_POSITION 0

//TRM TDC disgnostic error TDC ID mask/position
#define TRM_TDC_DIAGNOSTIC_ERROR_TDC_ID_MASK 0x00007800
#define TRM_TDC_DIAGNOSTIC_ERROR_TDC_ID_POSITION 11 

//TRM global trailer word required bit pattern
//#define TRM_GLOBAL_TRAILER 0x50000000
#define TRM_GLOBAL_TRAILER 0x5000000f

//TRM event CRC mask/position
#define TRM_EVENT_CRC_MASK 0x0000fff0
#define TRM_EVENT_CRC_POSITION 4

//TRM event counter mask/position
#define TRM_EVENT_COUNTER_MASK 0x0fff0000
#define TRM_EVENT_COUNTER_POSITION 16


/******************************************
LTM DATA FORMAT                        
******************************************/

//LTM global header word required bit pattern
#define LTM_GLOBAL_HEADER 0x40000002

//LTM event words mask/position
#define LTM_EVENT_WORDS_MASK 0x0001fff0
#define LTM_EVENT_WORDS_POSITION 4

//LTM C-bit mask/position
#define LTM_C_BIT_MASK 0x00020000
#define LTM_C_BIT_POSITION 17

//LTM fault mask/position
#define LTM_FAULT_MASK 0x00fc0000
#define LTM_FAULT_POSITION 18

//PDL data 

//PDL value 1 mask/position
#define LTM_PDL_VALUE_1_MASK 0x000000ff
#define LTM_PDL_VALUE_1_POSITION 0

//PDL value 2 mask/position
#define LTM_PDL_VALUE_2_MASK 0x0000ff00
#define LTM_PDL_VALUE_2_POSITION 8

//PDL value 3 mask/position
#define LTM_PDL_VALUE_3_MASK 0x00ff0000
#define LTM_PDL_VALUE_3_POSITION 16

//PDL value 4 mask/position
#define LTM_PDL_VALUE_4_MASK 0xff000000
#define LTM_PDL_VALUE_4_POSITION 24

//ADC data 

//ADC value 1 mask/position
#define LTM_ADC_VALUE_1_MASK 0x000003ff
#define LTM_ADC_VALUE_1_POSITION 0

//ADC value 2 mask/position
#define LTM_ADC_VALUE_2_MASK 0x000ffc00
#define LTM_ADC_VALUE_2_POSITION 10

//ADC value 3 mask/position
#define LTM_ADC_VALUE_3_MASK 0x3ff00000
#define LTM_ADC_VALUE_3_POSITION 20

//LTM global trailer word required bit pattern
#define LTM_GLOBAL_TRAILER 0x50000002

//LTM event CRC mask/position
#define LTM_EVENT_CRC_MASK 0x0000fff0
#define LTM_EVENT_CRC_POSITION 4

//LTM event number mask/position
#define LTM_EVENT_NUMBER_MASK 0x0fff0000
#define LTM_EVENT_NUMBER_POSITION 16

/******************************************
 * END OF OLD DEFINITIONS
 ******************************************/

class TClonesArray;
class AliRawReader;
class AliTOFrawData;

class AliTOFRawStream: public TObject {
 public:

  AliTOFRawStream(); // default ctr
  AliTOFRawStream(AliRawReader* rawReader); // ctr
  virtual ~AliTOFRawStream(); // default dtr

  AliTOFRawStream(const AliTOFRawStream& stream); // copy ctr
  AliTOFRawStream& operator = (const AliTOFRawStream& stream); // ass. op.

  virtual Bool_t Next();
  
  virtual void   LoadRawData(Int_t indexDDL);

  Int_t GetDDL()        const {return fDDL;};
  Int_t GetTRM()        const {return fTRM;};
  Int_t GetTDC()        const {return fTDC;};
  Int_t GetTRMchain()   const {return fTRMchain;};
  Int_t GetTDCchannel() const {return fTDCchannel;};
  
  Int_t GetSector() const {return fSector;};
  Int_t GetPlate()  const {return fPlate;};
  Int_t GetStrip()  const {return fStrip;};
  Int_t GetPadZ()   const {return fPadZ;};
  Int_t GetPadX()   const {return fPadX;};
  
  Int_t GetTofBin() const {return fTime;};
  Int_t GetToTbin() const {return fToT;};
  Float_t GetLeadingEdge() const {return fLeadingEdge;};
  Float_t GetTrailingEdge() const {return fTrailingEdge;};

  Int_t GetPSbit()  const {return fPSbit;};
  Int_t GetACQ()    const {return fACQ;};
    
  Int_t GetErrorFlag()  const {return fErrorFlag;};

  Bool_t GetDecoderVersion() const {return fNewDecoderVersion;};

  void SetDDL(Int_t nDDL)            {fDDL = nDDL;};
  void SetTRM(Int_t nTRM)            {fTRM = nTRM;};
  void SetTDC(Int_t nTDC)            {fTDC = nTDC;};
  void SetTRMchain(Int_t nChain)     {fTRMchain = nChain;};
  void SetTDCchannel(Int_t nChannel) {fTDCchannel = nChannel;};

  void SetDecoderVersion(Bool_t version) {fNewDecoderVersion = version;};

  TClonesArray *GetRawData() const {return fTOFrawData;};

  void SetSector();
  void SetPlate();
  void SetStrip();
  void SetPadZ();
  void SetPadX();

  Bool_t GetBCCorrections() const {return fgApplyBCCorrections;}; // getter for the BC application switch
  Int_t GetDDLBCshift(Int_t ddl) const {return fgkddlBCshift[ddl];}; // getter for the DDL BC shift
  Int_t GetLocalEventCounterDRM() const {return fLocalEventCounterDRM;}; // getter for the DRM event counter
  Int_t GetLocalEventCounterLTM() const {return fLocalEventCounterLTM;}; // getter for the LTM event counter
  Int_t GetLocalEventCounterTRM(Int_t trm) const {return fLocalEventCounterTRM[trm];}; // getter for the TRM event counter
  Int_t GetLocalEventCounterChain(Int_t trm, Int_t chain) const {return fLocalEventCounterChain[trm][chain];}; // getter for the chain event counter
  Int_t GetChainBunchID(Int_t trm, Int_t chain) const {return fChainBunchID[trm][chain];}; // getter for the chain BC ID


  void Raw2Digits(AliRawReader* rawReader, TClonesArray * const digitsArray);
  void Raw2SDigits(AliRawReader* rawReader, TClonesArray * const sdigitsArray);

  static void  EquipmentId2VolumeId(Int_t nDDL, Int_t nTRM, Int_t iChain,
				    Int_t iTDC, Int_t iCH, Int_t *volume);
  void  EquipmentId2VolumeId(AliTOFHitData *hitData, Int_t *volume) const;
  static Int_t Equip2VolNplate(Int_t iDDL, Int_t nTRM, Int_t nTDC);
  static Int_t Equip2VolNstrip(Int_t iDDL, Int_t nTRM, Int_t nTDC);
  static Int_t Equip2VolNpad(Int_t iDDL, Int_t iChain, Int_t nTDC, Int_t iCH);
  static Int_t Equip2VolNpadX(Int_t iDDL, Int_t iChain, Int_t nTDC, Int_t iCH);
  static Int_t Equip2VolNpadZ(Int_t iDDL, Int_t iChain, Int_t nTDC, Int_t iCH);
  static Int_t GetDDLnumberPerSector(Int_t nDDL);
  static Int_t GetSectorNumber(Int_t nDDL);

  static Int_t Geant2DDL(Int_t vol[]);
  static Int_t Geant2TRM(Int_t vol[]);
  static Int_t Geant2TDC(Int_t vol[]);
  static Int_t Geant2Chain(Int_t vol[]);
  static Int_t Geant2Channel(Int_t vol[]);
  static void  Geant2EquipmentId(Int_t vol[], Int_t eqId[]);

  Bool_t DecodeDDL(Int_t DDLMin, Int_t DDLMax, Int_t verbose);
  Bool_t Decode(Int_t verbose);
  Bool_t DecodeV2(Int_t verbose);
  AliTOFDecoder *GetDecoder() const {return fDecoder;};
  AliTOFDecoderV2 *GetDecoderV2() const {return fDecoderV2;};
  void SetV2718Patch(Bool_t V2718Patch = kTRUE) {fDecoder->SetV2718Patch(V2718Patch);};

  void SetRawReader(AliRawReader * const rawReader) {fRawReader=rawReader;};

  AliTOFHitDataBuffer * GetDataBuffer(Int_t DDL) {return &fDataBuffer[DDL];};
  AliTOFHitDataBuffer * GetPackedDataBuffer(Int_t DDL) {return &fPackedDataBuffer[DDL];};

  void ResetDataBuffer(Int_t DDL) {fDataBuffer[DDL].Reset();};
  void ResetPackedDataBuffer(Int_t DDL) {fPackedDataBuffer[DDL].Reset();};

  void ResetBuffers();

  Bool_t LoadRawDataBuffers(Int_t indexDDL, Int_t verbose = 0);
  Bool_t LoadRawDataBuffersV2(Int_t indexDDL, Int_t verbose = 0);
  static void ApplyBCCorrections(Bool_t Value = kTRUE) {fgApplyBCCorrections = Value;};
  
  Int_t GetEventID() const {return fEventID;}; // getter for the eventID1 (bunch crossing) in the common data header

  void LTM2VolumeID(Int_t iDDL,
		    Int_t iTRM,
		    Int_t iChain,
		    Int_t iTDC,
		    Int_t iChannel,
		    Int_t detind0[], Int_t detind1[]) const;
  void VolumeID2LTM(Int_t detind[],
		    Int_t iDDL = -1,
		    Int_t iTRM = -1,
		    Int_t iChain = -1,
		    Int_t iTDC = -1,
		    Int_t iChannel = -1) const;

  enum ETOFRawStreamError {
    kPadXError = 0,
    kPadAlongStripError = 1,
    kPlateError = 2,
    kStripError = 3,
    kSectorError = 4,
    kDDLMinError = 5,
    kDDLMaxError = 6,
    kDDLdataReading = 7,
    kDDLDecoder = 8
  };

 private:
  
  Int_t GetField(UInt_t word, Int_t fieldMask, Int_t fieldPosition) const;
  
  AliRawReader*  fRawReader; // object for reading the raw data

  TClonesArray *fTOFrawData; // pointer to AliTOFrawData TClonesArray

  AliTOFDecoder *fDecoder; //pointer to TOF decoder
  AliTOFDecoderV2 *fDecoderV2; //pointer to TOF decoder

  Int_t         fDDL;          // DDL file number [0;71]
  Int_t         fTRM;          // TRM number [1;12]
  Int_t         fTRMchain;     // TRM chain number [0;1]
  Int_t         fTDC;          // TDC number [0;14]
  Int_t         fTDCchannel;   // TDC channel number [0;7]
  Int_t         fTime;         // time-of-flight measurement [0;8191]
  Int_t         fToT;          // time-over-threshould measurement [0;255]
  Int_t         fLeadingEdge;  // leading edge measurement
  Int_t         fTrailingEdge; // trailing edge measurement
  Int_t         fErrorFlag;    // error flag
  
  Int_t         fSector;     // sector number [0;17]
  Int_t         fPlate;      // plate number [0;4]
  Int_t         fStrip;      // strip number [0;14/18]
  Int_t         fPadX;       // pad number along the strip [0;47]
  Int_t         fPadZ;       // pad-row number [0;1]

  Int_t fPackedDigits;       // counter for packed digits

  Int_t fWordType;           // word type
  Int_t fSlotID;             // crate slot ID number
  Int_t fACQ;                // flag to identif the aquisition kind
  Int_t fPSbit;              // flag for packing 
  Int_t fTDCerrorFlag;       // TDC error flag
  Bool_t fInsideDRM;         // inside/outside DRM
  Bool_t fInsideTRM;         // inside/outside TRM
  Bool_t fInsideLTM;         // inside/outside LTM
  Bool_t fInsideTRMchain0;   // inside/outside chain 0
  Bool_t fInsideTRMchain1;   // inside/outside chain 1

  AliTOFHitDataBuffer fDataBuffer[72]; // AliTOFHitDataBuffer
  AliTOFHitDataBuffer fPackedDataBuffer[72]; // AliTOFHitDataBuffer

  Int_t fLocalEventCounterDRM;          // event counter recorded in the DRM global trailer
  Int_t fLocalEventCounterLTM;          // event counter recorded in the LTM global trailer
  Int_t fLocalEventCounterTRM[13];      // event counter recorded in the TRMs global trailer
  Int_t fLocalEventCounterChain[13][2]; // event counter recorded in the chains trailer
  Int_t fChainBunchID[13][2];           // BC ID recorded in the chains header

  //AliTOFCableLengthMap * fCableLengthMap; // Pointer to the map of Amphenol cable length

  Int_t fEventID; // event ID1 in the common data header

  Bool_t fNewDecoderVersion;   // setting whether to use the new decoder version
                               //  -true -> new version
                               //  -false ->old version (default value!!)

  static const Int_t fgkddlBCshift[72]; // DDL BC shifts
  static Bool_t fgApplyBCCorrections; // switch to choose if apply or not the BC shift corrections

  static const Int_t fgkStrip0MapCrate0[];   // 1st strip number in crate 0
  static const Int_t fgkStrip1MapCrate0[];   // 2nd strip number in crate 0
  static const Int_t fgkStrip0MapCrate1[];   // 1st strip number in crate 1
  static const Int_t fgkStrip1MapCrate1[];   // 2nd strip number in crate 1
  static const Int_t fgkStrip0MapCrate2[];   // 1st strip number in crate 2
  static const Int_t fgkStrip1MapCrate2[];   // 2nd strip number in crate 2
  static const Int_t fgkStrip0MapCrate3[];   // 1st strip number in crate 3
  static const Int_t fgkStrip1MapCrate3[];   // 2nd strip number in crate 3

  static const Int_t fgkModule0MapCrate0[];  // 1st module number in crate 0
  static const Int_t fgkModule1MapCrate0[];  // 2nd module number in crate 0
  static const Int_t fgkModule0MapCrate1[];  // 1st module number in crate 1
  static const Int_t fgkModule1MapCrate1[];  // 2nd module number in crate 1
  static const Int_t fgkModule0MapCrate2[];  // 1st module number in crate 2
  static const Int_t fgkModule1MapCrate2[];  // 2nd module number in crate 2
  static const Int_t fgkModule0MapCrate3[];  // 1st module number in crate 3
  static const Int_t fgkModule1MapCrate3[];  // 2nd module number in crate 3

  static const Int_t fgkChannelMap0[5][19];  // mapping padX<24 <-> TDC channels
  static const Int_t fgkChannelMap24[5][19]; // mapping padX>=24 <-> TDC channels
  static const Int_t fgkChainMap0[5][19];    // mapping padX<24 <-> TRM chain
  static const Int_t fgkChainMap24[5][19];   // mapping padX>=24 <-> TRM chain


  ClassDef(AliTOFRawStream, 3)  // class for reading TOF raw digits
};

#endif
