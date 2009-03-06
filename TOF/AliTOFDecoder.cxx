/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
***************************************************************************/

/*
$Log$
Revision 1.5  2007/11/24 14:58:34  zampolli
New Method implemented (#DDL <-> TOF channels)

Revision 1.4  2007/05/18 13:08:57  decaro
Coding convention: RS1 violation -> suppression

Revision 1.3  2007/05/08 11:56:05  arcelli
improved verbosity in verbose mode (R.Preghenella)

Revision 1.2  2007/05/03 11:34:43  decaro
Coding convention: RS1 violation -> suppression

Revision 1.1  2007/04/27 11:00:32  arcelli
TOF Raw Data decoder

  author: Roberto Preghenella (R+), preghenella@bo.infn.it
*/


//////////////////////////////////////////////////////////////////////
//                                                                  //
//                                                                  //
//   Class for raw data decoding                                    //
//                                                                  //
//                                                                  //
//////////////////////////////////////////////////////////////////////
                               


#include "AliLog.h"
#include "AliTOFHitData.h"
#include "AliTOFHitDataBuffer.h"
#include "AliTOFDecoder.h"
#include "AliTOFGeometry.h"
#include "AliRawDataHeader.h"
#include "AliTOFRawDataFormat.h"

ClassImp(AliTOFDecoder)

//_________________________________________________________________

AliTOFDecoder::AliTOFDecoder() :
  TObject(),
  fVerbose(0),
  fV2718Patch(kFALSE),
  fDataBuffer(0x0),
  fPackedDataBuffer(0x0),
  //fTRMGlobalHeader(0x0),
  //fTRMGlobalTrailer(0x0),
  //fTRMChainHeader(0x0),
  //fTRMChainTrailer(0x0),
  //fTDCPackedHit(0x0),
  //fTDCUnpackedHit(0x0),
  //fTRMTDCError(0x0),
  //fTRMDiagnosticErrorWord1(0x0),
  //fTRMDiagnosticErrorWord2(0x0),
  fSpiderCurrentSlotID(-1),
  fSpiderCurrentChain(-1),
  fSpiderCurrentTDC(-1)
{
  //default constructor
}

//_________________________________________________________________

AliTOFDecoder::AliTOFDecoder(AliTOFHitDataBuffer *DataBuffer, AliTOFHitDataBuffer *PackedDataBuffer) :
  TObject(),
  fVerbose(0),
  fV2718Patch(kFALSE),
  fDataBuffer(DataBuffer),
  fPackedDataBuffer(PackedDataBuffer),
  //fTRMGlobalHeader(0x0),
  //fTRMGlobalTrailer(0x0),
  //fTRMChainHeader(0x0),
  //fTRMChainTrailer(0x0),
  //fTDCPackedHit(0x0),
  //fTDCUnpackedHit(0x0),
  //fTRMTDCError(0x0),
  //fTRMDiagnosticErrorWord1(0x0),
  //fTRMDiagnosticErrorWord2(0x0),
  fSpiderCurrentSlotID(-1),
  fSpiderCurrentChain(-1),
  fSpiderCurrentTDC(-1)
{
  //another constructor
}

//_________________________________________________________________

AliTOFDecoder::AliTOFDecoder(const AliTOFDecoder &source) : 
  TObject(source),
  fVerbose(source.fVerbose),
  fV2718Patch(source.fV2718Patch),
  fDataBuffer(source.fDataBuffer),
  fPackedDataBuffer(source.fPackedDataBuffer),
  //fTRMGlobalHeader(source.fTRMGlobalHeader),
  //fTRMGlobalTrailer(source.fTRMGlobalTrailer),
  //fTRMChainHeader(source.fTRMChainHeader),
  //fTRMChainTrailer(source.fTRMChainTrailer),
  //fTDCPackedHit(source.fTDCPackedHit),
  //fTDCUnpackedHit(source.fTDCUnpackedHit),
  //fTRMTDCError(source.fTRMTDCError),
  //fTRMDiagnosticErrorWord1(source.fTRMDiagnosticErrorWord1),
  //fTRMDiagnosticErrorWord2(source.fTRMDiagnosticErrorWord2),
  fSpiderCurrentSlotID(source.fSpiderCurrentSlotID),
  fSpiderCurrentChain(source.fSpiderCurrentChain),
  fSpiderCurrentTDC(source.fSpiderCurrentTDC)
{
  //copy constructor
  
}

//_________________________________________________________________

AliTOFDecoder &
AliTOFDecoder::operator = (const AliTOFDecoder &source)
{
  //operator =

  if (this == &source)
    return *this;

  TObject::operator=(source);
  fVerbose = source.fVerbose;
  fV2718Patch = source.fV2718Patch;
  fDataBuffer = source.fDataBuffer;
  fPackedDataBuffer = source.fPackedDataBuffer;
  //fTRMGlobalHeader = source.fTRMGlobalHeader;
  //fTRMGlobalTrailer = source.fTRMGlobalTrailer;
  //fTRMChainHeader = source.fTRMChainHeader;
  //fTRMChainTrailer = source.fTRMChainTrailer;
  //fTDCPackedHit = source.fTDCPackedHit;
  //fTDCUnpackedHit = source.fTDCUnpackedHit;
  //fTRMTDCError = source.fTRMTDCError;
  //fTRMDiagnosticErrorWord1 = source.fTRMDiagnosticErrorWord1;
  //fTRMDiagnosticErrorWord2 = source.fTRMDiagnosticErrorWord2;
  fSpiderCurrentSlotID = source.fSpiderCurrentSlotID;
  fSpiderCurrentChain = source.fSpiderCurrentChain;
  fSpiderCurrentTDC = source.fSpiderCurrentTDC;
  return *this;
}

AliTOFDecoder::~AliTOFDecoder()
{}

//_________________________________________________________________

Bool_t
AliTOFDecoder::Decode(UInt_t *rawData, Int_t nWords, const AliRawDataHeader *cdh)
{
  /* main decoding routine.
   * it loops over nWords 32-bit words 
   * starting at *rawData and decodes them.
   * it also fills some buffers in order to
   * have the decoded data available for other
   * classes.
   */

  AliTOFTRMGlobalHeader          *lTRMGlobalHeader; //TRM global header
  AliTOFTRMGlobalTrailer         *lTRMGlobalTrailer; //TRM global trailer
  AliTOFTRMChainHeader           *lTRMChainHeader; //TRM chain header
  //AliTOFTRMChainTrailer          *lTRMChainTrailer; //TRM chain trailer
  AliTOFTDCPackedHit             *lTDCPackedHit; //TDC packed hit
  AliTOFTDCUnpackedHit           *lTDCUnpackedHit; //TDC unpacked hit
  //AliTOFTRMTDCError              *lTRMTDCError; //TRM TDC error
  //AliTOFTRMDiagnosticErrorWord1  *lTRMDiagnosticErrorWord1; //TRM diagnostic error word 1
  //AliTOFTRMDiagnosticErrorWord2  *lTRMDiagnosticErrorWord2; //TRM diagnostica error word 2


  AliTOFHitData hitData;
  
  //useful variables
  Int_t   status;
  Short_t tempPS;
  Float_t tempTOT; //ns
  Int_t   tempTOTBin; //TOT_BIN_WIDTH

  //decoder variables
  UShort_t decodeStatus = 0x0;
  Short_t  currentDDL = -1;
  Short_t  currentSlotID = -1;
  Short_t  currentACQ = -1;
  Short_t  currentChain = -1;
  Short_t  currentBunchID = -1;
  Short_t  currentMiniEventID = cdh ? cdh->GetMiniEventID() : (Short_t)-1;
  Short_t  currentEventID1 = cdh ? cdh->GetEventID1() : (Short_t)-1;
  AliDebug(1, Form("EvID1 = %d, EvID2 = %d, currentMiniEventID = %d", currentEventID1, cdh->GetEventID2(), currentMiniEventID));
  if (!cdh)
    AliWarning("CDH not valid: deltaBunchID not reliable ");

  /*** V2718 patch ***/
  if (fV2718Patch){
    decodeStatus = decodeStatus | DRM_BIT;
    if (fVerbose)
      AliInfo("DRM not present: - V2718 patch decoding -");
  }
  /*** V2718 patch ***/

  if (fVerbose==2)
    AliInfo("Initialize SPIDER function");
  status = InitializeSpider();
  
  if (fVerbose)
    AliInfo("Start decoding");
  
  if (fVerbose)
    AliInfo("Loop over the data and decode");
  
  if (fVerbose)
    AliInfo("  St    Hex Word \t   Decoded Word");
  
  //loop over raw data
  for (Int_t iWord = 0; iWord < nWords; iWord++, rawData++){
    
    //switch word type
    switch (*rawData & WORD_TYPE_MASK){
      
    case GLOBAL_HEADER:
      
      //switch slot ID
      switch (*rawData & SLOT_ID_MASK){
	
	//DRM global header (slotID=1)
      case 1:
	//check decode status
	if ( decodeStatus != DRM_HEADER_STATUS ){
	  AliError(Form("  %02x - 0x%08x [ERROR] Unexpected DRM global header",decodeStatus,*rawData));
	  return kTRUE;
	}
	//decode status ok
	if (fVerbose)
	  AliInfo(Form("  %02x - 0x%08x \t  DRM global header",decodeStatus,*rawData));
	//change decode status
	decodeStatus = decodeStatus | DRM_BIT;
	
	//skip DRM data
	for (Int_t i = 0; i < DRM_DATA_WORDS; i++, iWord++, rawData++){
	  if (fVerbose)
	    AliInfo(Form("  %02x - 0x%08x \t  DRM data",decodeStatus,*rawData));
	}
	break;
	
	//LTM global header (slotID=2)
      case 2:
	//check decode status
	if ( decodeStatus != LTM_HEADER_STATUS ){
	  AliError(Form("  %02x - 0x%08x [ERROR] Unexpected LTM global header",decodeStatus,*rawData));
	  return kTRUE;
	}
	//decode status ok
	if (fVerbose)
	  AliInfo(Form("  %02x - 0x%08x \t  LTM global header",decodeStatus,*rawData));
	//change decode status
	decodeStatus = decodeStatus | LTM_BIT;
	
	//skip LTM data
	for (Int_t i = 0; i < LTM_DATA_WORDS; i++, iWord++, rawData++){
	  if (fVerbose)
	    AliInfo(Form("  %02x - 0x%08x \t  LTM data",decodeStatus,*rawData));
	}
	break;
	
	//TRM global header (slotID=3-12)
      case 3: case 4: case 5: case 6: case 7: case 8: case 9: case 10: case 11: case 12:
	//check decode status
	if ( decodeStatus != TRM_HEADER_STATUS ){
	  AliError(Form("  %02x - 0x%08x [ERROR] Unexpected TRM global header",decodeStatus,*rawData));
	  return kTRUE;
	}
	//decode status ok
	//set TRM global header
	lTRMGlobalHeader = (AliTOFTRMGlobalHeader*)rawData;	
	//set current TRM
	currentSlotID = lTRMGlobalHeader->GetSlotID();
	currentACQ = lTRMGlobalHeader->GetACQBits();
	if (fVerbose)
	  AliInfo(Form("  %02x - 0x%08x \t  TRM global header \t slotID=%02d ACQ=%01d L=%01d",decodeStatus,*rawData,lTRMGlobalHeader->GetSlotID(),lTRMGlobalHeader->GetACQBits(),lTRMGlobalHeader->GetLBit()));
	//change decode status
	decodeStatus = decodeStatus | TRM_BIT;
	break;
	
      default:
	AliError(Form("  %02x - 0x%08x [ERROR] Not valid slotID in global header",decodeStatus,*rawData));
	return kTRUE;
	break;
	
      }
      //end switch slotID
      break;
      
    case GLOBAL_TRAILER:
      
      //switch slot ID
      switch (*rawData & SLOT_ID_MASK){
	
	//DRM global trailer (slotID=1)
      case 1:
	//check decode status
	if ( decodeStatus != DRM_TRAILER_STATUS ){
	  AliError(Form("  %02x - 0x%08x [ERROR] Unexpected DRM global trailer",decodeStatus,*rawData));
	  return kTRUE;
	}
	//decode status ok
	if (fVerbose)
	  AliInfo(Form("  %02x - 0x%08x \t  DRM global trailer",decodeStatus,*rawData));
	//change decode status
	decodeStatus = decodeStatus & ~DRM_BIT;
	break;
	
	//LTM global trailer (slotID=2)
      case 2:
	//check decode status
	if ( decodeStatus != LTM_TRAILER_STATUS ){
	  AliError(Form("  %02x - 0x%08x [ERROR] Unexpected LTM global trailer",decodeStatus,*rawData));
	  return kTRUE;
	}
	//decode status ok
	if (fVerbose)
	  AliInfo(Form("  %02x - 0x%08x \t  LTM global trailer",decodeStatus,*rawData));
	//change decode status
	decodeStatus = decodeStatus & ~LTM_BIT;
	break;
	
	//TRM global trailer (slotID=15)
      case 15:
	//check decode status
	if ( decodeStatus != TRM_TRAILER_STATUS ){
	  AliError(Form("  %02x - 0x%08x [ERROR] Unexpected TRM global trailer",decodeStatus,*rawData));
	  return kTRUE;
	}
	//decode status ok
	//set TRM global trailer
	lTRMGlobalTrailer = (AliTOFTRMGlobalTrailer *)rawData;	
	if (fVerbose)
	  AliInfo(Form("  %02x - 0x%08x \t  TRM global trailer \t CRC=%04d eventCounter=%04d",decodeStatus,*rawData,lTRMGlobalTrailer->GetEventCRC(),lTRMGlobalTrailer->GetEventCounter()));
	//change decode status
	decodeStatus = decodeStatus & ~TRM_BIT;
	break; 
	
      default:
	AliError(Form("  %02x - 0x%08x [ERROR] Not valid slotID/pattern in global trailer",decodeStatus,*rawData));
	return kTRUE;
	break;
      }
      break;
      
    case CHAIN_A_HEADER:
      //check decode status
      if ( (decodeStatus != CHAIN_A_HEADER_STATUS) ){
	AliError(Form("  %02x - 0x%08x [ERROR] Unexpected TRM chain A header",decodeStatus,*rawData));
	return kTRUE;
      }
      //decode status ok
      lTRMChainHeader = (AliTOFTRMChainHeader *)rawData;
      currentChain = 0;
      currentBunchID = lTRMChainHeader->GetBunchID();
      if (fVerbose)
	AliInfo(Form("  %02x - 0x%08x \t  TRM chain A header \t chain=%01d bunchID=%04d",decodeStatus,*rawData,currentChain,currentBunchID));
      //change decode status
      decodeStatus = decodeStatus | CHAIN_A_BIT;
      break;
      
    case CHAIN_A_TRAILER:
      //check decode status
      if ( decodeStatus != CHAIN_A_TRAILER_STATUS ){
	AliError(Form("  %02x - 0x%08x [ERROR] Unexpected TRM chain A trailer",decodeStatus,*rawData));
	return kTRUE;
      }
      //decode status ok
      if (fVerbose)
	AliInfo(Form("  %02x - 0x%08x \t  TRM chain A trailer",decodeStatus,*rawData));
      //change decode status
      decodeStatus = decodeStatus & ~CHAIN_A_BIT;
      break;
      
    case CHAIN_B_HEADER:
      //check decode status
      if ( decodeStatus != CHAIN_B_HEADER_STATUS ){
	AliError(Form("  %02x - 0x%08x [ERROR] Unexpected TRM chain B header",decodeStatus,*rawData));
	return kTRUE;
      }
      //decode status ok
      lTRMChainHeader = (AliTOFTRMChainHeader *)rawData;
      currentChain = 1;
      currentBunchID = lTRMChainHeader->GetBunchID();
      if (fVerbose)
	AliInfo(Form("  %02x - 0x%08x \t  TRM chain B header \t chain=%01d bunchID=%04d",decodeStatus,*rawData,currentChain,currentBunchID));
      //change decode status
      decodeStatus = decodeStatus | CHAIN_B_BIT;
      break;
      
    case CHAIN_B_TRAILER:
      //check decode status
      if ( decodeStatus != CHAIN_B_TRAILER_STATUS ){
	AliError(Form("  %02x - 0x%08x [ERROR] Unexpected TRM chain B trailer",decodeStatus,*rawData));
	return kTRUE;
      }
      //decode status ok
      if (fVerbose)
	AliInfo(Form("  %02x - 0x%08x \t  TRM chain B trailer",decodeStatus,*rawData));
      //change decode status
      decodeStatus = decodeStatus & ~CHAIN_B_BIT;
      break;
      
    case ERROR:
      if (fVerbose)
	AliInfo(Form("  %02x - 0x%08x \t  TDC error",decodeStatus,*rawData));
      break;
      
    case FILLER:
      if (fVerbose)
	AliInfo(Form("  %02x - 0x%08x \t  Filler",decodeStatus,*rawData));
      break;
      
    default:
      //check decode status
      if ( decodeStatus != CHAIN_A_TDC_HIT_STATUS &&
	   decodeStatus != CHAIN_B_TDC_HIT_STATUS ){
	AliError(Form("  %02x - 0x%08x [ERROR] Unexpected or unknown word",decodeStatus,*rawData));
	return kTRUE;
      }
      //decode status ok
      
      //switch TRM ACQ
      switch (currentACQ){
	
      case PACKING_ENABLED_ACQ:
	//decode TDC packed hit
	lTDCPackedHit = (AliTOFTDCPackedHit *)rawData;
	lTDCUnpackedHit = (AliTOFTDCUnpackedHit *)rawData;
	//set hit in the equipment data
	hitData.SetDDLID(currentDDL);
	hitData.SetSlotID(currentSlotID);
	hitData.SetACQ(currentACQ);
	hitData.SetChain(currentChain);
	hitData.SetPS(lTDCPackedHit->GetPSBits());
	hitData.SetTDC(lTDCPackedHit->GetTDCID());
	hitData.SetChan(lTDCPackedHit->GetChan());
	hitData.SetTime((float)lTDCPackedHit->GetHitTime() * TIME_BIN_WIDTH);
	hitData.SetTimeBin(lTDCPackedHit->GetHitTime());
	hitData.SetTOT((float)lTDCPackedHit->GetTOTWidth() * TOT_BIN_WIDTH);
	hitData.SetTOTBin(lTDCPackedHit->GetTOTWidth());
	hitData.SetDeltaBunchID(currentBunchID - currentEventID1);
	//orphane leading hit
	if (hitData.GetPS()==LEADING_HIT_PS){
	  hitData.SetTime((float)lTDCUnpackedHit->GetHitTime() * TIME_BIN_WIDTH);
	  hitData.SetTimeBin(lTDCUnpackedHit->GetHitTime());
	  //set TOT to zero
	  hitData.SetTOT(0);
	  hitData.SetTOTBin(0);
	  //push hit data in packed data buffer
	  if (fPackedDataBuffer != 0x0)
	    fPackedDataBuffer->Add(hitData);
	  //set TOT to not measured
	  hitData.SetTOT(-1);
	  hitData.SetTOTBin(-1);
	  //push hit data in packed data buffer
	  if (fDataBuffer != 0x0)
	    fDataBuffer->Add(hitData);
	}
	//orphane trailing hit
	else if (hitData.GetPS()==TRAILING_HIT_PS){
	  hitData.SetTime((float)lTDCUnpackedHit->GetHitTime() * TIME_BIN_WIDTH);
	  hitData.SetTimeBin(lTDCUnpackedHit->GetHitTime());
	  //set TOT to not measured
	  hitData.SetTOT(-1);
	  hitData.SetTOTBin(-1);
	  //push hit data in data buffer
	  if (fDataBuffer != 0x0)
	    fDataBuffer->Add(hitData);
	}
	//packed hit and OVF
	else{
	  //push hit data in packed data buffer
	  if (fPackedDataBuffer != 0x0)
	    fPackedDataBuffer->Add(hitData);
	  //save PS temporary
	  tempPS = hitData.GetPS();
	  //save TOT temporary
	  tempTOT = hitData.GetTOT();
	  tempTOTBin = hitData.GetTOTBin();
	  //unpack the hit: leading hit
	  hitData.SetPS(LEADING_HIT_PS);
	  //set TOT to not measured
	  hitData.SetTOT(-1);
	  hitData.SetTOTBin(-1);
	  //push leading hit data in data buffer
	  if (fDataBuffer != 0x0)
	    fDataBuffer->Add(hitData);
	  //unpack the hit: trailing hit
	  hitData.SetPS(TRAILING_HIT_PS);
	  hitData.SetTime(hitData.GetTime() + tempTOT);
	  hitData.SetTimeBin(hitData.GetTimeBin() + (Int_t)(tempTOTBin * TOT_TO_TIME_BIN_WIDTH));
	  //push trailing hit data in data buffer
	  if (fDataBuffer != 0x0)
	    fDataBuffer->Add(hitData);
	  //restore packed hit
	  hitData.SetPS(tempPS);
	  hitData.SetTime(hitData.GetTime() - tempTOT);
	  hitData.SetTimeBin(hitData.GetTimeBin() - (Int_t)(tempTOTBin * TOT_TO_TIME_BIN_WIDTH));
	  hitData.SetTOT(tempTOT);
	  hitData.SetTOTBin(tempTOTBin);
	}
	
	if (fVerbose)
	  switch (hitData.GetPS()){
	  case PACKED_HIT_PS:
	    AliInfo(Form("  %02x - 0x%08x \t  TDC hit [packed] \t PS=%01d TDC=%01d chan=%01d TOT=%3.1fns (%d) time=%4.1fns (%d)",decodeStatus,*rawData,hitData.GetPS(),hitData.GetTDC(),hitData.GetChan(),hitData.GetTOT(),hitData.GetTOTBin(),hitData.GetTime(),hitData.GetTimeBin()));
	    break;
	  case LEADING_HIT_PS:
	    AliInfo(Form("  %02x - 0x%08x \t  TDC hit [orp.lead] \t PS=%01d TDC=%01d chan=%01d time=%4.1fns (%d)",decodeStatus,*rawData,hitData.GetPS(),hitData.GetTDC(),hitData.GetChan(),hitData.GetTime(),hitData.GetTimeBin()));
	    break;
	  case TRAILING_HIT_PS:
	    AliInfo(Form("  %02x - 0x%08x \t  TDC hit [orp.trai] \t PS=%01d TDC=%01d chan=%01d time=%4.1fns (%d)",decodeStatus,*rawData,hitData.GetPS(),hitData.GetTDC(),hitData.GetChan(),hitData.GetTime(),hitData.GetTimeBin()));
	    break;
	  case TOT_OVF_HIT_PS:
	    AliInfo(Form("  %02x - 0x%08x \t  TDC hit [TOT ovfl] \t PS=%01d TDC=%01d chan=%01d TOT=%3.1fns (%d) time=%4.1fns (%d)",decodeStatus,*rawData,hitData.GetPS(),hitData.GetTDC(),hitData.GetChan(),hitData.GetTOT(),hitData.GetTOTBin(),hitData.GetTime(),hitData.GetTimeBin()));
	    break;
	  }
	break;
	
      case LEADING_ONLY_ACQ: case TRAILING_ONLY_ACQ:
	//decode TDC unpacked hit
	lTDCUnpackedHit = (AliTOFTDCUnpackedHit *)rawData;
	//set hit in the equipment data
	hitData.SetDDLID(currentDDL);
	hitData.SetSlotID(currentSlotID);
	hitData.SetACQ(currentACQ);
	hitData.SetChain(currentChain);
	hitData.SetPS(lTDCUnpackedHit->GetPSBits());
	hitData.SetTDC(lTDCUnpackedHit->GetTDCID());
	hitData.SetChan(lTDCUnpackedHit->GetChan());
	hitData.SetTime((float)lTDCUnpackedHit->GetHitTime() * TIME_BIN_WIDTH);
	hitData.SetTimeBin(lTDCUnpackedHit->GetHitTime());
	hitData.SetTOT(-1.);
	hitData.SetTOTBin(-1);
	hitData.SetDeltaBunchID(currentBunchID - currentEventID1);
	//push hit data in data buffer
	  if (fDataBuffer != 0x0)
	    fDataBuffer->Add(hitData);
	
	if (fVerbose)
	  switch (hitData.GetPS()){
	  case PACKED_HIT_PS:
	    AliInfo(Form("  %02x - 0x%08x \t  TDC hit [packed] \t PS=%01d TDC=%01d chan=%01d TOT=%3.1fns (%d) time=%4.1fns (%d)",decodeStatus,*rawData,hitData.GetPS(),hitData.GetTDC(),hitData.GetChan(),hitData.GetTOT(),hitData.GetTOTBin(),hitData.GetTime(),hitData.GetTimeBin()));
	    break;
	  case LEADING_HIT_PS:
	    AliInfo(Form("  %02x - 0x%08x \t  TDC hit [leading] \t PS=%01d TDC=%01d chan=%01d time=%4.1fns (%d)",decodeStatus,*rawData,hitData.GetPS(),hitData.GetTDC(),hitData.GetChan(),hitData.GetTime(),hitData.GetTimeBin()));
	    break;
	  case TRAILING_HIT_PS:
	    AliInfo(Form("  %02x - 0x%08x \t  TDC hit [trailing] \t PS=%01d TDC=%01d chan=%01d time=%4.1fns (%d)",decodeStatus,*rawData,hitData.GetPS(),hitData.GetTDC(),hitData.GetChan(),hitData.GetTime(),hitData.GetTimeBin()));
	    break;
	  case TOT_OVF_HIT_PS:
	    AliInfo(Form("  %02x - 0x%08x \t  TDC hit [TOT ovfl] \t PS=%01d TDC=%01d chan=%01d TOT=%3.1fns (%d) time=%4.1fns (%d)",decodeStatus,*rawData,hitData.GetPS(),hitData.GetTDC(),hitData.GetChan(),hitData.GetTOT(),hitData.GetTOTBin(),hitData.GetTime(),hitData.GetTimeBin()));
	    break;
	  }
	break;
	
      case PACKING_DISABLED_ACQ:
	//decode TDC unpacked hit
	lTDCUnpackedHit = (AliTOFTDCUnpackedHit *)rawData;
	//set hit in the equipment data
	hitData.SetDDLID(currentDDL);
	hitData.SetSlotID(currentSlotID);
	hitData.SetACQ(currentACQ);
	hitData.SetChain(currentChain);
	hitData.SetPS(lTDCUnpackedHit->GetPSBits());
	hitData.SetTDC(lTDCUnpackedHit->GetTDCID());
	hitData.SetChan(lTDCUnpackedHit->GetChan());
	hitData.SetTime((float)lTDCUnpackedHit->GetHitTime() * TIME_BIN_WIDTH);
	hitData.SetTimeBin(lTDCUnpackedHit->GetHitTime());
	hitData.SetTOT(-1.);
	hitData.SetTOTBin(-1);
	hitData.SetDeltaBunchID(currentBunchID - currentEventID1);
	//push hit data in data buffer
	  if (fDataBuffer != 0x0)
	    fDataBuffer->Add(hitData);
	
	if (fVerbose)
	  switch (hitData.GetPS()){
	  case PACKED_HIT_PS:
	    AliInfo(Form("  %02x - 0x%08x \t  TDC hit [packed] \t PS=%01d TDC=%01d chan=%01d TOT=%3.1fns (%d) time=%4.1fns (%d)",decodeStatus,*rawData,hitData.GetPS(),hitData.GetTDC(),hitData.GetChan(),hitData.GetTOT(),hitData.GetTOTBin(),hitData.GetTime(),hitData.GetTimeBin()));
	    break;
	  case LEADING_HIT_PS:
	    AliInfo(Form("  %02x - 0x%08x \t  TDC hit [leading] \t PS=%01d TDC=%01d chan=%01d time=%4.1fns (%d)",decodeStatus,*rawData,hitData.GetPS(),hitData.GetTDC(),hitData.GetChan(),hitData.GetTime(),hitData.GetTimeBin()));
	    break;
	  case TRAILING_HIT_PS:
	    AliInfo(Form("  %02x - 0x%08x \t  TDC hit [trailing] \t PS=%01d TDC=%01d chan=%01d time=%4.1fns (%d)",decodeStatus,*rawData,hitData.GetPS(),hitData.GetTDC(),hitData.GetChan(),hitData.GetTime(),hitData.GetTimeBin()));
	    break;
	  case TOT_OVF_HIT_PS:
	    AliInfo(Form("  %02x - 0x%08x \t  TDC hit [TOT ovfl] \t PS=%01d TDC=%01d chan=%01d TOT=%3.1fns (%d) time=%4.1fns (%d)",decodeStatus,*rawData,hitData.GetPS(),hitData.GetTDC(),hitData.GetChan(),hitData.GetTOT(),hitData.GetTOTBin(),hitData.GetTime(),hitData.GetTimeBin()));
	    break;
	  }
	//call spider function
	if (fVerbose==2)
	  AliInfo("Calling SPIDER function");
	Spider(hitData);
	break;
      }
      //end switch TRM ACQ
      
      
    }
    //end switch word type
    
  }
  //end equipment data loop
  
  if (fVerbose)
    AliInfo("End of data loop");
  
  if (fVerbose==2)
    AliInfo("Reset SPIDER function");
  status = ResetSpider();
  
  if (fVerbose)
    AliInfo("Decoder is exiting succesfully.");

  return kFALSE;  
}

//_________________________________________________________________

Bool_t 
AliTOFDecoder::InitializeSpider(){

  /* SPIDER initialization routine.
     it initializes SPIDER variables in order
     to have SPIDER ready to pack tof data 
     in packed data objects
  */
  
  if (fVerbose==2)
    AliInfo("Initializing SPIDER");
  
  fSpiderCurrentSlotID=-1;
  fSpiderCurrentChain=-1;
  fSpiderCurrentTDC=-1;
  
  for (Int_t chan=0;chan<N_CHANNEL;chan++)
    fSpiderLeadingFlag[chan] = kFALSE;
  
  return kFALSE;
}

//_________________________________________________________________

Bool_t 
AliTOFDecoder::ResetSpider(){

  /* SPIDER reset routine.
     it resets SPIDER buffers and 
     variables in order to empty full 
     buffers a set up SIPDER for new
     HPTDC data
  */

  if (fVerbose==2)
    AliInfo("Resetting SPIDER buffers");

  for (Int_t chan=0;chan<N_CHANNEL;chan++){
    if (fSpiderLeadingFlag[chan]){
      if (fVerbose==2)
	AliInfo("Buffer non empty: put leading hit into buffer as orphane");
      //set TOT to zero
      fSpiderLeadingHit[chan].SetACQ(4);
      fSpiderLeadingHit[chan].SetPS(1);
      fSpiderLeadingHit[chan].SetTOT(0);
      fSpiderLeadingHit[chan].SetTOTBin(0);
      //push hit into packed buffer
      if (fPackedDataBuffer != 0x0)
	fPackedDataBuffer->Add(fSpiderLeadingHit[chan]);
      if (fVerbose==2)
	AliInfo(Form("Packed hit: slotID=%d chain=%d PS=%01d TDC=%01d chan=%01d TOT=%3.1fns (%d) time=%4.1fns (%d)",fSpiderLeadingHit[chan].GetSlotID(),fSpiderLeadingHit[chan].GetChain(),fSpiderLeadingHit[chan].GetPS(),fSpiderLeadingHit[chan].GetTDC(),fSpiderLeadingHit[chan].GetChan(),fSpiderLeadingHit[chan].GetTOT(),fSpiderLeadingHit[chan].GetTOTBin(),fSpiderLeadingHit[chan].GetTime(),fSpiderLeadingHit[chan].GetTimeBin()));
      
    }
    fSpiderLeadingFlag[chan]=kFALSE;
  }
  
  return kFALSE;
}

//_________________________________________________________________

Bool_t 
AliTOFDecoder::Spider(AliTOFHitData &hitData){

  /* main SPIDER routine.
     it receives, reads, stores and packs
     unpacked HPTDC data in packed data
     object. it also fills buffers.
  */
 
  Int_t status;

  if (fVerbose==2)
    AliInfo("Hit data received");

  //check if TDC is changed (slotID,chain,TDC triplet)
  if (fSpiderCurrentSlotID!=hitData.GetSlotID() ||
      fSpiderCurrentChain!=hitData.GetChain() ||
      fSpiderCurrentTDC!=hitData.GetTDC() ){
    if (fVerbose==2)
      AliInfo("Data coming from a new TDC: reset buffers");
    //reset spider
    status = ResetSpider();
    //set current TDC 
    fSpiderCurrentSlotID=hitData.GetSlotID();
    fSpiderCurrentChain=hitData.GetChain();
    fSpiderCurrentTDC=hitData.GetTDC();
  }
  
  //switch PS bits
  switch (hitData.GetPS()){

  case LEADING_HIT_PS:
    if (fVerbose==2)
      AliInfo(Form("Leading hit: slotID=%d chain=%d PS=%01d TDC=%01d chan=%01d TOT=%3.1fns (%d) time=%4.1fns (%d)",hitData.GetSlotID(),hitData.GetChain(),hitData.GetPS(),hitData.GetTDC(),hitData.GetChan(),hitData.GetTOT(),hitData.GetTOTBin(),hitData.GetTime(),hitData.GetTimeBin()));
    //check spider leading flag
    if (fSpiderLeadingFlag[hitData.GetChan()]){
      if (fVerbose==2)
	AliInfo("Leading hit: buffer full, put previous in buffers as orphane and keep current");
      //set TOT at zero for previous hit
      fSpiderLeadingHit[hitData.GetChan()].SetACQ(4);
      fSpiderLeadingHit[hitData.GetChan()].SetPS(1);
      fSpiderLeadingHit[hitData.GetChan()].SetTOT(0);
      fSpiderLeadingHit[hitData.GetChan()].SetTOTBin(0);
      //push previous hit into packed buffer
      if (fPackedDataBuffer != 0x0)
	fPackedDataBuffer->Add(fSpiderLeadingHit[hitData.GetChan()]);
      //set current hit
      fSpiderLeadingHit[hitData.GetChan()]=hitData;
      if (fVerbose==2)
	AliInfo(Form("Packed hit: slotID=%d chain=%d PS=%01d TDC=%01d chan=%01d TOT=%3.1fns (%d) time=%4.1fns (%d)",fSpiderLeadingHit[hitData.GetChan()].GetSlotID(),fSpiderLeadingHit[hitData.GetChan()].GetChain(),fSpiderLeadingHit[hitData.GetChan()].GetPS(),fSpiderLeadingHit[hitData.GetChan()].GetTDC(),fSpiderLeadingHit[hitData.GetChan()].GetChan(),fSpiderLeadingHit[hitData.GetChan()].GetTOT(),fSpiderLeadingHit[hitData.GetChan()].GetTOTBin(),fSpiderLeadingHit[hitData.GetChan()].GetTime(),fSpiderLeadingHit[hitData.GetChan()].GetTimeBin()));
    }
    else{
      if (fVerbose==2)
	AliInfo("Leading hit: buffer empty, keep current hit and set flag");
      fSpiderLeadingHit[hitData.GetChan()]=hitData;
      //set spider leading flag
      fSpiderLeadingFlag[hitData.GetChan()]=kTRUE;
    }
    break;

  case TRAILING_HIT_PS:
    if (fVerbose==2)
      AliInfo(Form("Trailing hit: slotID=%d chain=%d PS=%01d TDC=%01d chan=%01d TOT=%3.1fns (%d) time=%4.1fns (%d)",hitData.GetSlotID(),hitData.GetChain(),hitData.GetPS(),hitData.GetTDC(),hitData.GetChan(),hitData.GetTOT(),hitData.GetTOTBin(),hitData.GetTime(),hitData.GetTimeBin()));
    //check spider leading flag
    if (fSpiderLeadingFlag[hitData.GetChan()]){
      if (fVerbose==2)
	AliInfo("Trailing hit: buffer full, pack leading and trailing");
      hitData.SetACQ(4);
      hitData.SetPS(0);
      hitData.SetTOT(hitData.GetTime()-fSpiderLeadingHit[hitData.GetChan()].GetTime());
      hitData.SetTOTBin((Int_t)((hitData.GetTimeBin()-fSpiderLeadingHit[hitData.GetChan()].GetTimeBin())*TIME_TO_TOT_BIN_WIDTH));
      hitData.SetTime(fSpiderLeadingHit[hitData.GetChan()].GetTime());
      hitData.SetTimeBin(fSpiderLeadingHit[hitData.GetChan()].GetTimeBin());
      //check TOT and set TOT overflow if TOT < 0
      if (hitData.GetTOT() < 0){
	hitData.SetPS(3);
	hitData.SetTOT(0);
	hitData.SetTOTBin(0);
      }
      if (fPackedDataBuffer != 0x0)
	fPackedDataBuffer->Add(hitData);      
      if (fVerbose==2)
	AliInfo(Form("Packed hit: slotID=%d chain=%d PS=%01d TDC=%01d chan=%01d TOT=%3.1fns (%d) time=%4.1fns (%d)",hitData.GetSlotID(),hitData.GetChain(),hitData.GetPS(),hitData.GetTDC(),hitData.GetChan(),hitData.GetTOT(),hitData.GetTOTBin(),hitData.GetTime(),hitData.GetTimeBin()));
      //unset spider leading flag
      fSpiderLeadingFlag[hitData.GetChan()]=kFALSE;
    }
    else{
      if (fVerbose==2)
	AliInfo("Trailing hit: buffer empty, throw hit away");
    }
    break;
  }
  //end switch PS bits

  return kFALSE;
}
//_____________________________________________________________________________
void AliTOFDecoder::GetArrayDDL(Int_t* array, Int_t ddl){

  // method that fills array with the
  // TOF channels indexes corresponding
  // to DDL iDDL

  AliTOFGeometry *geom = new AliTOFGeometry();
  Int_t indexDDL = ddl%4;
  Int_t iSector = Int_t(ddl/4);
  if (fVerbose)
    AliInfo(Form(" Sector = %i, DDL within sector = %i",iSector, indexDDL));

  Int_t volume[5];
  volume[0]=iSector;
  Int_t minPlate=0, maxPlate=0, minStrip2=0, maxStrip2=0, minPadz=0, maxPadz=0, minPadx=0, maxPadx=0;

  if (indexDDL==0){
    minPlate=kMinPlate0;
    maxPlate=kMaxPlate0;
    minStrip2=kMinStrip0;
    maxStrip2=kMaxStrip0;
    minPadz=kMinPadz0;
    maxPadz=kMaxPadz0;
    minPadx=kMinPadx0;
    maxPadx=kMaxPadx0;
  }

  else if (indexDDL==1){
    minPlate=kMinPlate1;
    maxPlate=kMaxPlate1;
    minStrip2=kMinStrip1;
    maxStrip2=kMaxStrip1;
    minPadz=kMinPadz1;
    maxPadz=kMaxPadz1;
    minPadx=kMinPadx1;
    maxPadx=kMaxPadx1;
  }

  else if (indexDDL==2){
    minPlate=kMinPlate2;
    maxPlate=kMaxPlate2;
    minStrip2=kMinStrip2;
    maxStrip2=kMaxStrip2;
    minPadz=kMinPadz2;
    maxPadz=kMaxPadz2;
    minPadx=kMinPadx2;
    maxPadx=kMaxPadx2;
  }

  else if (indexDDL==3){
    minPlate=kMinPlate3;
    maxPlate=kMaxPlate3;
    minStrip2=kMinStrip3;
    maxStrip2=kMaxStrip3;
    minPadz=kMinPadz3;
    maxPadz=kMaxPadz3;
    minPadx=kMinPadx3;
    maxPadx=kMaxPadx3;
  }

  Int_t ichTOF=0;

  Int_t minStrip=0;
  Int_t maxStrip=18;  
  for (Int_t iPlate=minPlate;iPlate<=maxPlate;iPlate++){
    if (iPlate==2) {
      maxStrip = maxStrip2;
      minStrip = minStrip2;
    }
    else {
      maxStrip = 18;
      minStrip = 0;
    }
    for (Int_t iStrip=minStrip;iStrip<=maxStrip;iStrip++){
      for (Int_t iPadz=minPadz;iPadz<=maxPadz;iPadz++){
	for (Int_t iPadx=minPadx;iPadx<=maxPadx;iPadx++){
	  volume[1]=iPlate;
	  volume[2]=iStrip;
	  volume[3]=iPadz;
	  volume[4]=iPadx;
	  if (fVerbose)
	    AliInfo(Form(" volume[0] = %i, volume[1] = %i, volume[2] = %i, volume[3] = %i, volume[4] = %i",volume[0],volume[1],volume[2],volume[3],volume[4]));

	  if (indexDDL==0 || indexDDL==2){
	    array[ichTOF]=geom->GetIndex(volume);
	    if (fVerbose)
	      AliInfo(Form(" ichTOF = %i, TOFChannel = %i",ichTOF,array[ichTOF]));

	  }
	  else {
	    array[ichTOF]=geom->GetIndex(volume);
	    if (fVerbose)
	      AliInfo(Form(" ichTOF = %i, TOFChannel = %i",ichTOF,array[ichTOF]));

	  }
	  ichTOF++;
	}
      }
    }
  }
  //AliInfo(Form("ichTOF = %i",ichTOF));
  if ((indexDDL%2==0 && ichTOF!=2160) ||
      (indexDDL%2==1 && ichTOF!=2208)) {
    AliWarning(Form("Something strange occurred, number of entries in array different from expected! Please, check! ichTOF = %i",ichTOF));
  }
  return;
}
