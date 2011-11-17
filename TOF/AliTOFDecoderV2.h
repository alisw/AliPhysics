#ifndef ALITOFDECODERV2_H
#define ALITOFDECODERV2_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: AliTOFDecoder.h,v 1.2 2007/05/08 11:55:24 arcelli Exp $ */

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
#define DRM_STATUS_HEADER_WORDS 4
#define LTM_PDL_DATA_WORDS      12
#define LTM_ADC_DATA_WORDS      20
#define LTM_OR_DATA_WORDS      16

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

//define CRC macros to convert 32-bit CRC into DRM/TRM CRC
#define COMPUTE_DRM_CRC(a) ( ((a & 0x0000ffff) >> 0)  ^\
			     ((a & 0xffff0000) >> 16) )
#define COMPUTE_TRM_CRC(a) ( ((a & 0x00000fff) >> 0)  ^\
			     ((a & 0x00fff000) >> 12) ^\
			     ((a & 0xff000000) >> 24) )
#define COMPUTE_LTM_CRC(a) ( ((a & 0x00000fff) >> 0)  ^\
			     ((a & 0x00fff000) >> 12) ^\
			     ((a & 0xff000000) >> 24) )


#include "TObject.h"
#include "AliRawReader.h"
#include "AliTOFTDCHit.h"

class AliTOFDecoderSummaryData;
class AliTOFDRMSummaryData;
class AliTOFLTMSummaryData;
class AliTOFTRMSummaryData;
class AliTOFChainSummaryData;
class AliTOFTDCHitBuffer;
class AliTOFTDCErrorBuffer;
struct AliRawDataHeader;

class AliTOFDecoderV2 : public TObject
{
 public:
  AliTOFDecoderV2(AliRawReader *reader = NULL); //default constructor
  AliTOFDecoderV2(const AliTOFDecoderV2 &source); //copy constructor
  AliTOFDecoderV2 &operator = (const AliTOFDecoderV2 &source); //operator =
  ~AliTOFDecoderV2(); //distructor
  /* setters */
  void SetRawReader(AliRawReader *value) {fRawReader = value;}; // set raw reader
  void SetVerbose(Bool_t Verbose = kTRUE) {fVerbose = Verbose;}; //set verbose level
  void SetLogErrors(Bool_t Value = kTRUE) {fLogErrors = Value;}; //set log errors
  void SetV2718Patch(Bool_t V2718Patch = kTRUE) {fV2718Patch = V2718Patch;}; //set V2718 patch (no DRM)
  void SetRecoverError(Bool_t RecoverError = kTRUE) {fRecoverError = RecoverError;}; //decoder will try to recover decoding errors
  void SetRecoverErrorThr(Int_t value) {fRecoverErrorThr = value;}; // setter
  void SetSpider(Bool_t value = kTRUE) {fSpider = value;}; //set spider
  void SetRunNumber(Int_t RunNumber) {fRunNumber = RunNumber;}; //set run number
  void SetEventNumber(UInt_t EventNumber) {fEventNumber = EventNumber;}; //set event number
  void SetEquipmentID(Int_t EquipmentID) {fEquipmentID = EquipmentID;}; //set equipment ID
  /* getters */
  AliTOFDecoderSummaryData *GetDecoderSummaryData() {return fDecoderSummaryData;}; //get decoder summary data
  /* methods */
  Bool_t Decode(UInt_t *rawData, UInt_t nWords); //main decode function
  void Spider(AliTOFTDCHit &hit);
  void ResetSpider();
  /* raw reader decoder interface */
  Bool_t NextEvent() {return fRawReader ? fRawReader->NextEvent() : kFALSE;}; // next event
  UInt_t GetEventType() {return fRawReader ? fRawReader->GetType() : 0;}; // get event type
  Bool_t DecodeNext(); // decode next
  const AliRawDataHeader *GetCDH() const {return fRawReader ? fRawReader->GetDataHeader() : NULL;}; // get CDH

 private:
  AliRawReader *fRawReader; // raw reader
  Bool_t fVerbose; //verbose flag
  Bool_t fLogErrors; //log errors flag
  Bool_t fV2718Patch; //V2718 patch flag
  Bool_t fRecoverError; //recover error flag
  Int_t fRecoverErrorThr; // recover error thr
  Bool_t fSpider; //spider flag
  Int_t  fRunNumber; //run number
  UInt_t  fEventNumber; //event number
  Int_t  fEquipmentID; //equipment ID

  //summary data pointers
  AliTOFDecoderSummaryData *fDecoderSummaryData; //decoder summary data
  AliTOFDRMSummaryData     *fDRMSummaryData; //DRM summary data
  AliTOFLTMSummaryData     *fLTMSummaryData; //LTM summary data
  AliTOFTRMSummaryData     *fTRMSummaryData; //TRM summary data
  AliTOFChainSummaryData   *fChainSummaryData; //chain summary data

  //buffer pointers
  AliTOFTDCHitBuffer       *fTDCHitBuffer; //TDC hit buffer
  AliTOFTDCHitBuffer       *fTDCPackedHitBuffer; //TDC packed hit buffer
  AliTOFTDCErrorBuffer     *fTDCErrorBuffer; //TDC error buffer

  //decoding objects
  AliTOFDRMGlobalHeader          *fDRMGlobalHeader; //DRM global header
  AliTOFDRMGlobalTrailer         *fDRMGlobalTrailer; //DRM global trailer
  AliTOFDRMStatusHeader1         *fDRMStatusHeader1; //DRM status header1
  AliTOFDRMStatusHeader2         *fDRMStatusHeader2; //DRM status header2
  AliTOFDRMStatusHeader3         *fDRMStatusHeader3; //DRM status header3
  AliTOFDRMStatusHeader4         *fDRMStatusHeader4; //DRM status header4
  AliTOFDRMEventCRC              *fDRMEventCRC; //DRM event CRC
  AliTOFLTMGlobalHeader          *fLTMGlobalHeader; //LTM global header
  AliTOFLTMGlobalTrailer         *fLTMGlobalTrailer; //LTM global trailer
  AliTOFLTMPDLData               *fLTMPDLData; //LTM PDL data
  AliTOFLTMADCData               *fLTMADCData; //LTM ADC data
  AliTOFLTMORData                *fLTMORData; //LTM OR data
  AliTOFTRMGlobalHeader          *fTRMGlobalHeader; //TRM global header
  AliTOFTRMGlobalTrailer         *fTRMGlobalTrailer; //TRM global trailer
  AliTOFTRMChainHeader           *fTRMChainHeader; //TRM chain header
  AliTOFTRMChainTrailer          *fTRMChainTrailer; //TRM chain trailer
  AliTOFTDCPackedHit             *fTDCPackedHit; //TDC packed hit
  AliTOFTDCUnpackedHit           *fTDCUnpackedHit; //TDC unpacked hit
  AliTOFTRMTDCError              *fTRMTDCError; //TRM TDC error
  AliTOFTRMDiagnosticErrorWord1  *fTRMDiagnosticErrorWord1; //TRM diagnostic error word 1
  AliTOFTRMDiagnosticErrorWord2  *fTRMDiagnosticErrorWord2; //TRM diagnostica error word 2

  /* Spider data members */
  AliTOFTDCHit fSpiderBuffer[8]; // SPIDER buffer
  Bool_t fSpiderBufferFull[8]; // SPIDER buffer full flag
  Int_t fSpiderTDCID; // SPIDER TDC ID
  AliTOFTDCHitBuffer *fSpiderTDCPackedHitBuffer; // SPIDER buffer


  /* Summary Data Functions */
  //fill DRM summary data 
  void FillDRMSummaryData(const AliTOFDRMGlobalHeader *DRMGlobalHeader); //DRM global header
  void FillDRMSummaryData(const AliTOFDRMGlobalTrailer *DRMGlobalTrailer); //DRM global trailer
  void FillDRMSummaryData(const AliTOFDRMStatusHeader1 *DRMStatusHeader1); //DRM status header 1
  void FillDRMSummaryData(const AliTOFDRMStatusHeader2 *DRMStatusHeader2); //DRM status header 2
  void FillDRMSummaryData(const AliTOFDRMStatusHeader3 *DRMStatusHeader3); //DRM status header 3
  void FillDRMSummaryData(const AliTOFDRMStatusHeader4 *DRMStatusHeader4); //DRM status header 4
  void FillDRMSummaryData(const AliTOFDRMEventCRC *DRMEventCRC); //DRM event CRC
  //fill LTM summary data
  void FillLTMSummaryData(const AliTOFLTMGlobalHeader *LTMGlobalHeader); //LTM global header
  void FillLTMSummaryData(const AliTOFLTMGlobalTrailer *LTMGlobalTrailer); //LTM global trailer
  void FillLTMSummaryData(const AliTOFLTMPDLData *LTMPDLData, Int_t PDLWord); //LTM PDL data
  void FillLTMSummaryData(const AliTOFLTMADCData *LTMADCData, Int_t ADCWord); //LTM ADC data
  void FillLTMSummaryData(const AliTOFLTMORData *LTMORData, Int_t ORWord); //LTM OR data
  //fill TRM summary data
  void FillTRMSummaryData(const AliTOFTRMGlobalHeader *TRMGlobalHeader); //TRM global header
  void FillTRMSummaryData(const AliTOFTRMGlobalTrailer *TRMGlobalTrailer); //TRM global trailer
  //fill chain summary data
  void FillChainSummaryData(const AliTOFTRMChainHeader *TRMChainHeader); //TRM chain header
  void FillChainSummaryData(const AliTOFTRMChainTrailer *TRMChainTrailer); //TRM chain trailer

  ClassDef(AliTOFDecoderV2, 1);
};

#endif /* ALITOFDECODERV2_H */
