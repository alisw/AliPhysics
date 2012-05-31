#ifndef ALITOFDECODER_H
#define ALITOFDECODER_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

///////////////////////////////////////////////////////////////
//                                                           //
//   This class provides the basic TOF raw data decoder.     //
//                                                           //
///////////////////////////////////////////////////////////////

//define decoder status and bits
#define DRM_BIT                 0x1
#define LTM_BIT                 0x2
#define TRM_BIT                 0x4
#define CHAIN_A_BIT             0x8
#define CHAIN_B_BIT             0x10

#define DRM_HEADER_STATUS       0x0
#define DRM_TRAILER_STATUS      (DRM_BIT)
#define LTM_HEADER_STATUS       (DRM_BIT)
#define LTM_TRAILER_STATUS      (DRM_BIT|LTM_BIT)
#define TRM_HEADER_STATUS       (DRM_BIT)
#define TRM_TRAILER_STATUS      (DRM_BIT|TRM_BIT)
#define CHAIN_A_HEADER_STATUS   (DRM_BIT|TRM_BIT)
#define CHAIN_A_TRAILER_STATUS  (DRM_BIT|TRM_BIT|CHAIN_A_BIT)
#define CHAIN_B_HEADER_STATUS   (DRM_BIT|TRM_BIT)
#define CHAIN_B_TRAILER_STATUS  (DRM_BIT|TRM_BIT|CHAIN_B_BIT)
#define CHAIN_A_TDC_HIT_STATUS  (DRM_BIT|TRM_BIT|CHAIN_A_BIT)
#define CHAIN_B_TDC_HIT_STATUS  (DRM_BIT|TRM_BIT|CHAIN_B_BIT)

//define DRM/LTM fixed number of words
#define DRM_DATA_WORDS          5
#define LTM_DATA_WORDS          48

//define masks
#define WORD_TYPE_MASK          0xf0000000
#define SLOT_ID_MASK            0x0000000f

//define word types
#define GLOBAL_HEADER           0x40000000
#define GLOBAL_TRAILER          0x50000000
#define CHAIN_A_HEADER          0x00000000
#define CHAIN_A_TRAILER         0x10000000
#define CHAIN_B_HEADER          0x20000000
#define CHAIN_B_TRAILER         0x30000000
#define ERROR                   0x60000000
#define FILLER                  0x70000000

//define TRM ACQ status
#define PACKING_ENABLED_ACQ     0x0
#define LEADING_ONLY_ACQ        0x1
#define TRAILING_ONLY_ACQ       0x2
#define PACKING_DISABLED_ACQ    0x3

//define TDC hit PS status
#define PACKED_HIT_PS           0x0
#define LEADING_HIT_PS          0x1
#define TRAILING_HIT_PS         0x2
#define TOT_OVF_HIT_PS          0x3

//define mandatory numbers
#define N_EQUIPMENT             72
#define N_DDL                   N_EQUIPMENT
#define N_TRM                   10
#define N_CHAIN                 2
#define N_TDC                   15
#define N_CHANNEL               8
#define TRM_FIRST_SLOT_ID       3
#define TRM_LAST_SLOT_ID        12

#define TIME_BIN_WIDTH          24.4e-3//ns
#define TOT_BIN_WIDTH           48.8e-3//ns
#define TIME_TO_TOT_BIN_WIDTH   ( TIME_BIN_WIDTH / TOT_BIN_WIDTH )
#define TOT_TO_TIME_BIN_WIDTH   ( TOT_BIN_WIDTH / TIME_BIN_WIDTH )

//max hit number in a TDC
#define MAX_TDC_HIT_NUMBER 100
//max TDC errors
#define MAX_TDC_ERROR_NUMBER 1000
//max hit number in a TRM
#define MAX_TRM_HIT_NUMBER 2400

#include "TObject.h"
#include "AliTOFRawDataFormat.h"
#include "AliTOFHitData.h"
//#include "AliTOFHitDataBuffer.h"
//#include "AliTOFGeometry.h"

class AliTOFHitDataBuffer;
struct AliRawDataHeader;

class AliTOFDecoder : public TObject
{
 public:
  AliTOFDecoder(); //default constructor
  AliTOFDecoder(AliTOFHitDataBuffer *DB, AliTOFHitDataBuffer *PDB); //constructor
  AliTOFDecoder(const AliTOFDecoder &source); //copy constructor
  AliTOFDecoder &operator = (const AliTOFDecoder &source); //operator =
  ~AliTOFDecoder(); //distructor
  
  Bool_t Decode(const UInt_t *rawData, Int_t nWords, const AliRawDataHeader *cdh); //main decode function
  void   SetVerbose(Int_t Verbose = 1) {fVerbose = Verbose;}; //set verbose level
  void   SetV2718Patch(Bool_t V2718Patch = kTRUE) {fV2718Patch = V2718Patch;}; //set V2718 patch (no DRM)
  void   SetDataBuffer(AliTOFHitDataBuffer *DB) {fDataBuffer = DB;}; //set up data buffer
  void   SetPackedDataBuffer(AliTOFHitDataBuffer *PDB) {fPackedDataBuffer = PDB;}; //set up packed data buffer

  void GetArrayDDL(Int_t* array, Int_t iDDL); // method to return array of TOFchannels corresponding to a given DDL id 
  void PrintStack(const UInt_t *rawData, Int_t nWords, const AliRawDataHeader *cdh);

enum {
  kMinPlate0=0,kMaxPlate0=2,
  kMinStrip0=0,kMaxStrip0=6,
  kMinPadz0=0,kMaxPadz0=1,
  kMinPadx0=0,kMaxPadx0=23,

  kMinPlate1=0,kMaxPlate1=2,
  kMinStrip1=0,kMaxStrip1=7,
  kMinPadz1=0,kMaxPadz1=1,
  kMinPadx1=24,kMaxPadx1=47,

  kMinPlate2=2,kMaxPlate2=4,
  kMinStrip2=8,kMaxStrip2=14,
  kMinPadz2=0,kMaxPadz2=1,
  kMinPadx2=24,kMaxPadx2=47,

  kMinPlate3=2,kMaxPlate3=4,
  kMinStrip3=7,kMaxStrip3=14,
  kMinPadz3=0,kMaxPadz3=1,
  kMinPadx3=0,kMaxPadx3=23

};


 private:
  /* SPIDER
   * - Software Packing Inside Decoding Routines -
   * developed by Roberto Preghenella (R+)
   * use at your own risk
   */     
  Bool_t InitializeSpider(); //initialize SPIDER routine
  Bool_t ResetSpider(); //reset SPIDER routine
  Bool_t Spider(AliTOFHitData &hitData); //main SPIDER routine

  Int_t                fVerbose; //verbose flag
  Bool_t               fV2718Patch; //V2718 patch flag
  AliTOFHitDataBuffer *fDataBuffer; //data buffer pointer
  AliTOFHitDataBuffer *fPackedDataBuffer; //packed data buffe pointer

  //decoding objects
  //AliTOFTRMGlobalHeader          *fTRMGlobalHeader; //TRM global header
  //AliTOFTRMGlobalTrailer         *fTRMGlobalTrailer; //TRM global trailer
  //AliTOFTRMChainHeader           *fTRMChainHeader; //TRM chain header
  //AliTOFTRMChainTrailer          *fTRMChainTrailer; //TRM chain trailer
  //AliTOFTDCPackedHit             *fTDCPackedHit; //TDC packed hit
  //AliTOFTDCUnpackedHit           *fTDCUnpackedHit; //TDC unpacked hit
  //AliTOFTRMTDCError              *fTRMTDCError; //TRM TDC error
  //AliTOFTRMDiagnosticErrorWord1  *fTRMDiagnosticErrorWord1; //TRM diagnostic error word 1
  //AliTOFTRMDiagnosticErrorWord2  *fTRMDiagnosticErrorWord2; //TRM diagnostica error word 2

  //SPIDER variables
  Int_t         fSpiderCurrentSlotID; //SPIDER current slot ID
  Int_t         fSpiderCurrentChain; //SPIDER current chain
  Int_t         fSpiderCurrentTDC; //SPIDER current TDC
  Bool_t        fSpiderLeadingFlag[8]; //SPIDER channel leading flag
  AliTOFHitData fSpiderLeadingHit[8]; //SPIDER channel leading hit

  ClassDef(AliTOFDecoder, 2)
};

#endif /* ALITOFDECODER_H */
