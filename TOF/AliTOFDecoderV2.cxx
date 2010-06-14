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
$Log: AliTOFDecoder.cxx,v $
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
#include "AliTOFDecoderV2.h"
#include "AliTOFTDCHit.h"

ClassImp(AliTOFDecoderV2)

//_________________________________________________________________

AliTOFDecoderV2::AliTOFDecoderV2(AliRawReader *reader) :
  TObject(),
  fRawReader(reader),
  fVerbose(kFALSE),
  fLogErrors(kTRUE),
  fV2718Patch(kFALSE),
  fRecoverError(kTRUE),
  fRecoverErrorThr(0),
  fSpider(kTRUE),
  fRunNumber(0),
  fEventNumber(0),
  fEquipmentID(0),
  fDecoderSummaryData(0x0),
  fDRMSummaryData(0x0),
  fLTMSummaryData(0x0),
  fTRMSummaryData(0x0),
  fChainSummaryData(0x0),
  fTDCHitBuffer(0x0),
  fTDCPackedHitBuffer(0x0),
  fTDCErrorBuffer(0x0),
  fDRMGlobalHeader(0x0),
  fDRMGlobalTrailer(0x0),
  fDRMStatusHeader1(0x0),
  fDRMStatusHeader2(0x0),
  fDRMStatusHeader3(0x0),
  fDRMStatusHeader4(0x0),
  fDRMEventCRC(0x0),
  fLTMGlobalHeader(0x0),
  fLTMGlobalTrailer(0x0),
  fLTMPDLData(0x0),
  fLTMADCData(0x0),
  fLTMORData(0x0),
  fTRMGlobalHeader(0x0),
  fTRMGlobalTrailer(0x0),
  fTRMChainHeader(0x0),
  fTRMChainTrailer(0x0),
  fTDCPackedHit(0x0),
  fTDCUnpackedHit(0x0),
  fTRMTDCError(0x0),
  fTRMDiagnosticErrorWord1(0x0),
  fTRMDiagnosticErrorWord2(0x0),
  fSpiderTDCID(-1),
  fSpiderTDCPackedHitBuffer(0x0)
{
  //default constructor
  if (fRawReader)
    fRawReader->Select("TOF", 0, 71);

  fDecoderSummaryData = new AliTOFDecoderSummaryData();

  for (Int_t iChan = 0; iChan < N_CHANNEL; iChan++)
    fSpiderBufferFull[iChan] = kFALSE;
}

//_________________________________________________________________

AliTOFDecoderV2::AliTOFDecoderV2(const AliTOFDecoderV2 &source) : 
  TObject(),
  fRawReader(source.fRawReader),
  fVerbose(source.fVerbose),
  fLogErrors(source.fLogErrors),
  fV2718Patch(source.fV2718Patch),
  fRecoverError(source.fRecoverError),
  fRecoverErrorThr(source.fRecoverErrorThr),
  fSpider(source.fSpider),
  fRunNumber(source.fRunNumber),
  fEventNumber(source.fEventNumber),
  fEquipmentID(source.fEquipmentID),
  fDecoderSummaryData(0x0),
  fDRMSummaryData(0x0),
  fLTMSummaryData(0x0),
  fTRMSummaryData(0x0),
  fChainSummaryData(0x0),
  fTDCHitBuffer(0x0),
  fTDCPackedHitBuffer(0x0),
  fTDCErrorBuffer(0x0),
  fDRMGlobalHeader(0x0),
  fDRMGlobalTrailer(0x0),
  fDRMStatusHeader1(0x0),
  fDRMStatusHeader2(0x0),
  fDRMStatusHeader3(0x0),
  fDRMStatusHeader4(0x0),
  fDRMEventCRC(0x0),
  fLTMGlobalHeader(0x0),
  fLTMGlobalTrailer(0x0),
  fLTMPDLData(0x0),
  fLTMADCData(0x0),
  fLTMORData(0x0),
  fTRMGlobalHeader(0x0),
  fTRMGlobalTrailer(0x0),
  fTRMChainHeader(0x0),
  fTRMChainTrailer(0x0),
  fTDCPackedHit(0x0),
  fTDCUnpackedHit(0x0),
  fTRMTDCError(0x0),
  fTRMDiagnosticErrorWord1(0x0),
  fTRMDiagnosticErrorWord2(0x0),
  fSpiderTDCID(-1),
  fSpiderTDCPackedHitBuffer(0x0)
{
  //copy constructor
  fDecoderSummaryData = new AliTOFDecoderSummaryData(*source.fDecoderSummaryData);
  
  for (Int_t iChan = 0; iChan < N_CHANNEL; iChan++)
    fSpiderBufferFull[iChan] = kFALSE;
}

//_________________________________________________________________

AliTOFDecoderV2 &
AliTOFDecoderV2::operator = (const AliTOFDecoderV2 &source)
{
  //operator =
  fRawReader = source.fRawReader;
  fVerbose = source.fVerbose;
  fLogErrors = source.fLogErrors;
  fV2718Patch = source.fV2718Patch;
  fRecoverError = source.fRecoverError;
  fRecoverErrorThr = source.fRecoverErrorThr;
  fSpider = source.fSpider;
  fRunNumber = source.fRunNumber;
  fEventNumber = source.fEventNumber;
  fEquipmentID = source.fEquipmentID;
  *fDecoderSummaryData = *source.fDecoderSummaryData;
  return *this;
}

AliTOFDecoderV2::~AliTOFDecoderV2()
{
  delete fDecoderSummaryData;
}

//_________________________________________________________________

Bool_t
AliTOFDecoderV2::Decode(UInt_t *rawData, UInt_t nWords)
{
  /* main decoding routine.
   * it loops over nWords 32-bit words 
   * starting at *rawData and decodes them.
   * it also fills some buffers in order to
   * have the decoded data available for other
   * classes.
   */

  //reset decoder summary data
  fDecoderSummaryData->Reset();

  //fill decoder summary data
  fDecoderSummaryData->SetRunNumber(fRunNumber);
  fDecoderSummaryData->SetEventNumber(fEventNumber);
  fDecoderSummaryData->SetEquipmentID(fEquipmentID);
  fDecoderSummaryData->SetInputWords(nWords);
  fDecoderSummaryData->SetRecoverError(fRecoverError);
  fDecoderSummaryData->SetSpider(fSpider);

  AliTOFTDCHit hit;
  AliTOFTDCError error;

  //decoder status
  UShort_t decoderStatus = 0x0;

  //CRC variables
  UInt_t DRMCRC = 0x0;
  UInt_t LTMCRC = 0x0;
  UInt_t TRMCRC = 0x0;

  // error warning counter
  Int_t errorWarning = 0;
  
  if (fRecoverError && fVerbose)
    AliInfo("Recover error option enabled: potentially dangerous!");

  /*** V2718 patch ***/
  if (fV2718Patch){
    decoderStatus = decoderStatus | DRM_BIT;
    fDecoderSummaryData->SetDecoderStatus(decoderStatus);
    fDecoderSummaryData->SetCurrentDRMID(0x0);
    fDecoderSummaryData->SetV2718Patch(kTRUE);
    fDRMSummaryData = fDecoderSummaryData->GetDRMSummaryData();
    fDRMSummaryData->SetHeader(kTRUE);
    fDRMSummaryData->SetDRMID(0x0);
    if (fVerbose)
      AliInfo("DRM not present: - V2718 patch decoding -");
  }
  /*** V2718 patch ***/

  if (fVerbose)
    AliInfo("Start decoding");
  
  if (fVerbose)
    AliInfo("Loop over the data and decode");
  
  if (fVerbose)
    AliInfo("  St    Hex Word \t   Decoded Word");
  
  //loop over raw data
  for (UInt_t iWord = 0; 
       iWord < nWords; 
       iWord++, rawData++, fDecoderSummaryData->SetDecodedWords(iWord)){
    
    //try to recover error
    if (fDecoderSummaryData->GetRecoveringError() && fVerbose)
      AliInfo(Form("  %02x - 0x%08x",decoderStatus,*rawData));
    
    //compute CRC with current data
    DRMCRC ^= *rawData;
    LTMCRC ^= *rawData;
    TRMCRC ^= *rawData;

    //switch word type
    switch (*rawData & WORD_TYPE_MASK){
      
    case GLOBAL_HEADER:
      
      //switch slot ID
      switch (*rawData & SLOT_ID_MASK){
	
	//DRM global header (slotID=1)
      case 1:
	//try to recover error
	if (fDecoderSummaryData->GetRecoveringError())
	  continue;
	//check decode status
	if ( decoderStatus != DRM_HEADER_STATUS ){
	  if (fLogErrors)
	    AliError(Form("  %02x - 0x%08x [ERROR] Unexpected DRM global header (curslot=%d)",decoderStatus,*rawData,fDecoderSummaryData->GetCurrentSlotID()));
	  fDecoderSummaryData->SetErrorDetected(kTRUE);
	  fDecoderSummaryData->SetErrorSlotID(fDecoderSummaryData->GetCurrentSlotID());
	  errorWarning++;
	  //try to recover error
	  if (fRecoverError){
	    if (errorWarning > fRecoverErrorThr) {
	      if (fVerbose)
		AliInfo("Trying to recover the error: searching for the next header");
	      fDecoderSummaryData->SetRecoveringError(kTRUE);
	      continue;
	    }
	    else {
	      if (fVerbose)
		AliInfo("Do not try to recover error yet, go on with decoding process");
	      continue;
	    }
	  }
	  return(fDecoderSummaryData->GetErrorDetected());
	}
	//decode status ok
	errorWarning = 0;
	//set DRM global header
	fDRMGlobalHeader = (AliTOFDRMGlobalHeader *)rawData;
	//reset DRM CRC
	DRMCRC = 0x0;
	//fill decoder summary data
	fDecoderSummaryData->SetCurrentDRMID(fDRMGlobalHeader->GetDRMID());
	fDecoderSummaryData->SetCurrentSlotID(fDRMGlobalHeader->GetSlotID());
	//get DRM summary data
	fDRMSummaryData = fDecoderSummaryData->GetDRMSummaryData();
	//reset DRM summary data
	fDRMSummaryData->Reset();
	//fill DRM summary data
	FillDRMSummaryData(fDRMGlobalHeader);
	//print verbose
	if (fVerbose)
	  AliInfo(Form("  %02x - 0x%08x \t  DRM global header",decoderStatus,*rawData));
	//change decode status
	decoderStatus = decoderStatus | DRM_BIT;
	fDecoderSummaryData->SetDecoderStatus(decoderStatus);
	//decode DRM status headers
	for (Int_t i = 0; i < DRM_STATUS_HEADER_WORDS; i++){
	  iWord++;
	  rawData++;
	  DRMCRC ^= *rawData;

	  switch (i){
	  case 0: //DRM status header 1
	    fDRMStatusHeader1 = (AliTOFDRMStatusHeader1 *)rawData;
	    FillDRMSummaryData(fDRMStatusHeader1);
	    if (fVerbose)
	      AliInfo(Form("  %02x - 0x%08x \t  DRM status header 1",decoderStatus,*rawData));
	    break;
	  case 1: //DRM status header 2
	    fDRMStatusHeader2 = (AliTOFDRMStatusHeader2 *)rawData;
	    FillDRMSummaryData(fDRMStatusHeader2);
	    if (fVerbose)
	      AliInfo(Form("  %02x - 0x%08x \t  DRM status header 2",decoderStatus,*rawData));
	    break;
	  case 2: //DRM status header 3
	    fDRMStatusHeader3 = (AliTOFDRMStatusHeader3 *)rawData;
	    FillDRMSummaryData(fDRMStatusHeader3);
	    if (fVerbose)
	      AliInfo(Form("  %02x - 0x%08x \t  DRM status header 3",decoderStatus,*rawData));
	    break;
	  case 3: //DRM status header 4
	    fDRMStatusHeader4 = (AliTOFDRMStatusHeader4 *)rawData;
	    FillDRMSummaryData(fDRMStatusHeader4);
	    if (fVerbose)
	      AliInfo(Form("  %02x - 0x%08x \t  DRM status header 4",decoderStatus,*rawData));
	    break;
	  }
	}
	//decode DRM event CRC
	iWord++;
	rawData++;
	DRMCRC ^= *rawData;
	//remove DRM event CRC from DRM CRC
	DRMCRC ^= *rawData;
	fDRMEventCRC = (AliTOFDRMEventCRC *)rawData;
	FillDRMSummaryData(fDRMEventCRC);
	if (fVerbose)
	  AliInfo(Form("  %02x - 0x%08x \t  DRM event CRC",decoderStatus,*rawData));
	break;
	
	//LTM global header (slotID=2)
      case 2:
	//recover error
	if (fDecoderSummaryData->GetRecoveringError()){
	  //change decode status
	  decoderStatus = LTM_HEADER_STATUS;
	  fDecoderSummaryData->SetDecoderStatus(decoderStatus);
	  fDecoderSummaryData->SetRecoveringError(kFALSE);
	  if (fVerbose)
	    AliInfo("LTM global header found: error probably recovered");
	}
	//check decode status
	if ( decoderStatus != LTM_HEADER_STATUS ){
	  if (fLogErrors)
	    AliError(Form("  %02x - 0x%08x [ERROR] Unexpected LTM global header (curslot=%d)",decoderStatus,*rawData,fDecoderSummaryData->GetCurrentSlotID()));
	  fDecoderSummaryData->SetErrorDetected(kTRUE);
	  fDecoderSummaryData->SetErrorSlotID(fDecoderSummaryData->GetCurrentSlotID());
	  errorWarning++;
	  //try to recover error
	  if (fRecoverError){
	    if (errorWarning > fRecoverErrorThr) {
	      if (fVerbose)
		AliInfo("Trying to recover the error: searching for the next header");
	      fDecoderSummaryData->SetRecoveringError(kTRUE);
	      continue;
	    }
	    else {
	      if (fVerbose)
		AliInfo("Do not try to recover error yet, go on with decoding process");
	      continue;
	    }
	  }
	  return(fDecoderSummaryData->GetErrorDetected());
	}
	//decode status ok
	errorWarning = 0;
	//set LTM global header
	fLTMGlobalHeader = (AliTOFLTMGlobalHeader *)rawData;
	//reset LTM CRC
	LTMCRC = 0x0;
	//fill decoder summary data
	fDecoderSummaryData->SetCurrentSlotID(fLTMGlobalHeader->GetSlotID());
	//get LTM summary data
	fLTMSummaryData = fDRMSummaryData->GetLTMSummaryData();
	//reset LTM summary data
	fLTMSummaryData->Reset();
	//fill LTM summary data
	FillLTMSummaryData(fLTMGlobalHeader);
	//set DRM slot enable mask bit
	fDRMSummaryData->SetDecoderSlotEnableMaskBit(fLTMGlobalHeader->GetSlotID() - 2);
	//print verbose
	if (fVerbose)
	  AliInfo(Form("  %02x - 0x%08x \t  LTM global header",decoderStatus,*rawData));
	//change decode status
	decoderStatus = decoderStatus | LTM_BIT;
	fDecoderSummaryData->SetDecoderStatus(decoderStatus);
	
	//decode LTM PDL data
	for (Int_t iPDLWord = 0; iPDLWord < LTM_PDL_DATA_WORDS; iPDLWord++){
	  iWord++;
	  rawData++;
	  DRMCRC ^= *rawData;
	  LTMCRC ^= *rawData;
	  //set LTM PDL data
	  fLTMPDLData = (AliTOFLTMPDLData *)rawData;
	  //fill LTM summary data
	  FillLTMSummaryData(fLTMPDLData, iPDLWord);
	  //print verbose
	  if (fVerbose)
	    AliInfo(Form("  %02x - 0x%08x \t  LTM PDL data \t\t PDL1=%03d PDL2=%03d PDL3=%03d PDL4=%03d",decoderStatus,*rawData,fLTMPDLData->GetPDLValue1(),fLTMPDLData->GetPDLValue2(),fLTMPDLData->GetPDLValue3(),fLTMPDLData->GetPDLValue4()));
	}
	//decode LTM ADC data
	for (Int_t iADCWord = 0; iADCWord < LTM_ADC_DATA_WORDS; iADCWord++){
	  iWord++;
	  rawData++;
	  DRMCRC ^= *rawData;
	  LTMCRC ^= *rawData;
	  //set LTM ADC data
	  fLTMADCData = (AliTOFLTMADCData *)rawData;
	  //fill LTM summary data
	  FillLTMSummaryData(fLTMADCData, iADCWord);
	  //print verbose
	  if (fVerbose)
	    AliInfo(Form("  %02x - 0x%08x \t  LTM ADC data \t\t ADC1=%04d ADC2=%04d ADC3=%04d",decoderStatus,*rawData,fLTMADCData->GetADCValue1(),fLTMADCData->GetADCValue2(),fLTMADCData->GetADCValue3()));
	}
	//decode LTM OR data
	for (Int_t iORWord = 0; iORWord < LTM_OR_DATA_WORDS; iORWord++){
	  iWord++;
	  rawData++;
	  DRMCRC ^= *rawData;
	  LTMCRC ^= *rawData;
	  //set LTM OR data
	  fLTMORData = (AliTOFLTMORData *)rawData;
	  //fill LTM summary data
	  FillLTMSummaryData(fLTMORData, iORWord);
	  //print verbose
	  if (fVerbose)
	    AliInfo(Form("  %02x - 0x%08x \t  LTM OR data \t\t ADC1=%04d ADC2=%04d ADC3=%04d",decoderStatus,*rawData,fLTMORData->GetORValue1(),fLTMORData->GetORValue2(),fLTMORData->GetORValue3()));
	}
	break;
	
	//TRM global header (slotID=3-12)
      case 3: case 4: case 5: case 6: case 7: case 8: case 9: case 10: case 11: case 12:
	//recover error
	if (fDecoderSummaryData->GetRecoveringError()){
	  //change decode status
	  decoderStatus = TRM_HEADER_STATUS;
	  fDecoderSummaryData->SetDecoderStatus(decoderStatus);
	  fDecoderSummaryData->SetRecoveringError(kFALSE);
	  if (fVerbose)
	    AliInfo("TRM global header found: error probably recovered");
	}
	//check decode status
	if ( decoderStatus != TRM_HEADER_STATUS ){
	  if (fLogErrors)
	    AliError(Form("  %02x - 0x%08x [ERROR] Unexpected TRM global header (curslot=%d)",decoderStatus,*rawData,fDecoderSummaryData->GetCurrentSlotID()));
	  fDecoderSummaryData->SetErrorDetected(kTRUE);
	  fDecoderSummaryData->SetErrorSlotID(fDecoderSummaryData->GetCurrentSlotID());
	  errorWarning++;
	  //try to recover error
	  if (fRecoverError){
	    if (errorWarning > fRecoverErrorThr) {
	      if (fVerbose)
		AliInfo("Trying to recover the error: searching for the next header");
	      fDecoderSummaryData->SetRecoveringError(kTRUE);
	      continue;
	    }
	    else {
	      if (fVerbose)
		AliInfo("Do not try to recover error yet, go on with decoding process");
	      continue;
	    }
	  }
	  return(fDecoderSummaryData->GetErrorDetected());
	}
	//decode status ok
	errorWarning = 0;
	//set TRM global header
	fTRMGlobalHeader = (AliTOFTRMGlobalHeader *)rawData;	
	//reset TRM CRC
	TRMCRC = 0x0;
	//fill decoder summary data
	fDecoderSummaryData->SetCurrentSlotID(fTRMGlobalHeader->GetSlotID());
	//get TRM summary data
	fTRMSummaryData = fDRMSummaryData->GetTRMSummaryData(fTRMGlobalHeader->GetSlotID() - TRM_FIRST_SLOT_ID);
	//reset TRM summary data
	fTRMSummaryData->Reset();
	//fill TRM summary data
	FillTRMSummaryData(fTRMGlobalHeader);
	//set DRM slot enable mask bit
	fDRMSummaryData->SetDecoderSlotEnableMaskBit(fTRMGlobalHeader->GetSlotID() - 2);
	//print verbose
	if (fVerbose)
	  AliInfo(Form("  %02x - 0x%08x \t  TRM global header \t slotID=%02d ACQ=%01d L=%01d",decoderStatus,*rawData,fTRMGlobalHeader->GetSlotID(),fTRMGlobalHeader->GetACQBits(),fTRMGlobalHeader->GetLBit()));
	//change decode status
	decoderStatus = decoderStatus | TRM_BIT;
	fDecoderSummaryData->SetDecoderStatus(decoderStatus);
	break;
	
      default:
	if (fLogErrors)
	  AliError(Form("  %02x - 0x%08x [ERROR] Not valid slotID in global header (curslot=%d)",decoderStatus,*rawData,fDecoderSummaryData->GetCurrentSlotID()));
	fDecoderSummaryData->SetErrorDetected(kTRUE);
	fDecoderSummaryData->SetErrorSlotID(fDecoderSummaryData->GetCurrentSlotID());
	  //try to recover error
	  if (fRecoverError){
	    if (errorWarning > fRecoverErrorThr) {
	      if (fVerbose)
		AliInfo("Trying to recover the error: searching for the next header");
	      fDecoderSummaryData->SetRecoveringError(kTRUE);
	      continue;
	    }
	    else {
	      if (fVerbose)
		AliInfo("Do not try to recover error yet, go on with decoding process");
	      continue;
	    }
	  }
	return(fDecoderSummaryData->GetErrorDetected());
	break;
	
      }
      //end switch slotID
      break;
      
    case GLOBAL_TRAILER:
   
   //switch slot ID
   switch (*rawData & SLOT_ID_MASK){
     
     //DRM global trailer (slotID=1)
   case 1:
     //recover error
     if (fDecoderSummaryData->GetRecoveringError()){
       //change decode status
       decoderStatus = DRM_TRAILER_STATUS;
       fDecoderSummaryData->SetDecoderStatus(decoderStatus);
       fDecoderSummaryData->SetRecoveringError(kFALSE);
       if (fVerbose)
	 AliInfo("DRM global trailer found: error probably recovered");
     }
     //check decode status
     if ( decoderStatus != DRM_TRAILER_STATUS ){
       if (fLogErrors)
	 AliError(Form("  %02x - 0x%08x [ERROR] Unexpected DRM global trailer (curslot=%d)",decoderStatus,*rawData,fDecoderSummaryData->GetCurrentSlotID()));
       fDecoderSummaryData->SetErrorDetected(kTRUE);
       fDecoderSummaryData->SetErrorSlotID(fDecoderSummaryData->GetCurrentSlotID());
	  errorWarning++;
	  //try to recover error
	  if (fRecoverError){
	    if (errorWarning > fRecoverErrorThr) {
	      if (fVerbose)
		AliInfo("Trying to recover the error: searching for the next header");
	      fDecoderSummaryData->SetRecoveringError(kTRUE);
	      continue;
	    }
	    else {
	      if (fVerbose)
		AliInfo("Do not try to recover error yet, go on with decoding process");
	      continue;
	    }
	  }
	  return(fDecoderSummaryData->GetErrorDetected());
	}
	//decode status ok
	errorWarning = 0;
	//set DRM global trailer
	fDRMGlobalTrailer = (AliTOFDRMGlobalTrailer *)rawData;
	//remove global trailer from DRM CRC
	DRMCRC ^= *rawData;
	//fill DRM summary data
	FillDRMSummaryData(fDRMGlobalTrailer);
	fDRMSummaryData->SetDecoderCRC(COMPUTE_DRM_CRC(DRMCRC));
	//print verbose
	if (fVerbose)
	  AliInfo(Form("  %02x - 0x%08x \t  DRM global trailer",decoderStatus,*rawData));
	//change decode status
	decoderStatus = decoderStatus & ~DRM_BIT;
	fDecoderSummaryData->SetDecoderStatus(decoderStatus);
	break;
	
	//LTM global trailer (slotID=2)
      case 2:
	//try to recover error
	if (fDecoderSummaryData->GetRecoveringError())
	  continue;
   	//check decode status
	if ( decoderStatus != LTM_TRAILER_STATUS ){
	  if (fLogErrors)
	    AliError(Form("  %02x - 0x%08x [ERROR] Unexpected LTM global trailer (curslot=%d)",decoderStatus,*rawData,fDecoderSummaryData->GetCurrentSlotID()));
	  fDecoderSummaryData->SetErrorDetected(kTRUE);
	  fDecoderSummaryData->SetErrorSlotID(fDecoderSummaryData->GetCurrentSlotID());
	  errorWarning++;
	  //try to recover error
	  if (fRecoverError){
	    if (errorWarning > fRecoverErrorThr) {
	      if (fVerbose)
		AliInfo("Trying to recover the error: searching for the next header");
	      fDecoderSummaryData->SetRecoveringError(kTRUE);
	      continue;
	    }
	    else {
	      if (fVerbose)
		AliInfo("Do not try to recover error yet, go on with decoding process");
	      continue;
	    }
	  }
	  return(fDecoderSummaryData->GetErrorDetected());
	}
	//decode status ok
	errorWarning = 0;
	//set LTM global trailer
	fLTMGlobalTrailer = (AliTOFLTMGlobalTrailer *)rawData;
	//remove global trailer from LTM CRC
	LTMCRC ^= *rawData;
	//fill LTM summary data
	FillLTMSummaryData(fLTMGlobalTrailer);
	fLTMSummaryData->SetDecoderCRC(COMPUTE_LTM_CRC(LTMCRC));
	//print verbose
	if (fVerbose)
	  AliInfo(Form("  %02x - 0x%08x \t  LTM global trailer",decoderStatus,*rawData));
	//change decode status
	decoderStatus = decoderStatus & ~LTM_BIT;
	fDecoderSummaryData->SetDecoderStatus(decoderStatus);
	break;
	
	//TRM global trailer (slotID=15)
      case 15:
	//try to recover error
	if (fDecoderSummaryData->GetRecoveringError())
	  continue;
	//check decode status
	if ( decoderStatus != TRM_TRAILER_STATUS ){
	  if (fLogErrors)
	    AliError(Form("  %02x - 0x%08x [ERROR] Unexpected TRM global trailer (curslot=%d)",decoderStatus,*rawData,fDecoderSummaryData->GetCurrentSlotID()));
	  fDecoderSummaryData->SetErrorDetected(kTRUE);
	  fDecoderSummaryData->SetErrorSlotID(fDecoderSummaryData->GetCurrentSlotID());
	  errorWarning++;
	  //try to recover error
	  if (fRecoverError){
	    if (errorWarning > fRecoverErrorThr) {
	      if (fVerbose)
		AliInfo("Trying to recover the error: searching for the next header");
	      fDecoderSummaryData->SetRecoveringError(kTRUE);
	      continue;
	    }
	    else {
	      if (fVerbose)
		AliInfo("Do not try to recover error yet, go on with decoding process");
	      continue;
	    }
	  }
	  return(fDecoderSummaryData->GetErrorDetected());
	}
	//decode status ok
	errorWarning = 0;
	//set TRM global trailer
	fTRMGlobalTrailer = (AliTOFTRMGlobalTrailer *)rawData;	
	//remove global trailer from TRM CRC
	TRMCRC ^= *rawData;
	//fill TRM summary data
	FillTRMSummaryData(fTRMGlobalTrailer);
	fTRMSummaryData->SetDecoderCRC(COMPUTE_TRM_CRC(TRMCRC));
	//print verbose
	if (fVerbose)
	  AliInfo(Form("  %02x - 0x%08x \t  TRM global trailer \t CRC=%04d eventCounter=%04d",decoderStatus,*rawData,fTRMGlobalTrailer->GetEventCRC(),fTRMGlobalTrailer->GetEventCounter()));
	//change decode status
	decoderStatus = decoderStatus & ~TRM_BIT;
	fDecoderSummaryData->SetDecoderStatus(decoderStatus);
	break; 
	
      default:
	//try to recover error
	if (fDecoderSummaryData->GetRecoveringError())
	  continue;
	if (fLogErrors)
	  AliError(Form("  %02x - 0x%08x [ERROR] Not valid slotID/pattern in global trailer (curslot=%d)",decoderStatus,*rawData,fDecoderSummaryData->GetCurrentSlotID()));
	fDecoderSummaryData->SetErrorDetected(kTRUE);
	fDecoderSummaryData->SetErrorSlotID(fDecoderSummaryData->GetCurrentSlotID());
	  errorWarning++;
	//try to recover error
	if (fRecoverError){
	  if (fVerbose)
	    AliInfo("Trying to recover the error: searching for the next header");
	  fDecoderSummaryData->SetRecoveringError(kTRUE);
	  continue;
	}
	return(fDecoderSummaryData->GetErrorDetected());
	break;
      }
      break;
      
    case CHAIN_A_HEADER:

      //try to recover error
      if (fDecoderSummaryData->GetRecoveringError())
	continue;
      //check decode status
      if ( decoderStatus != CHAIN_A_HEADER_STATUS  && !fDecoderSummaryData->GetRecoveringError() ){
	if (fLogErrors)
	  AliError(Form("  %02x - 0x%08x [ERROR] Unexpected TRM chain A header (curslot=%d)",decoderStatus,*rawData,fDecoderSummaryData->GetCurrentSlotID()));
	fDecoderSummaryData->SetErrorDetected(kTRUE);
	fDecoderSummaryData->SetErrorSlotID(fDecoderSummaryData->GetCurrentSlotID());
	  errorWarning++;
	  //try to recover error
	  if (fRecoverError){
	    if (errorWarning > fRecoverErrorThr) {
	      if (fVerbose)
		AliInfo("Trying to recover the error: searching for the next header");
	      fDecoderSummaryData->SetRecoveringError(kTRUE);
	      continue;
	    }
	    else {
	      if (fVerbose)
		AliInfo("Do not try to recover error yet, go on with decoding process");
	      continue;
	    }
	  }
	return(fDecoderSummaryData->GetErrorDetected());
      }
      //decode status ok
	errorWarning = 0;
      //set TRM chain header
      fTRMChainHeader = (AliTOFTRMChainHeader *)rawData;
      //fill decoder summary data
      fDecoderSummaryData->SetCurrentChain(0);
      //get chain summary data
      fChainSummaryData = fTRMSummaryData->GetChainSummaryData(0);
      //reset chain summary data
      fChainSummaryData->Reset();
      //fill chain summary data
      FillChainSummaryData(fTRMChainHeader);
      //get tdc hit buffer
      fTDCHitBuffer = fChainSummaryData->GetTDCHitBuffer();
      //reset tdc hit buffer
      fTDCHitBuffer->Reset();
      //get tdc packed hit buffer
      fTDCPackedHitBuffer = fChainSummaryData->GetTDCPackedHitBuffer();
      //reset tdc packed hit buffer
      fTDCPackedHitBuffer->Reset();
      //get tdc error buffer
      fTDCErrorBuffer = fChainSummaryData->GetTDCErrorBuffer();
      //reset tdc error buffer
      fTDCErrorBuffer->Reset();
      //print verbose
      if (fVerbose)
	AliInfo(Form("  %02x - 0x%08x \t  TRM chain A header \t chain=%01d bunchID=%04d",decoderStatus,*rawData,0,fTRMChainHeader->GetBunchID()));
      //change decode status
      decoderStatus = decoderStatus | CHAIN_A_BIT;
      fDecoderSummaryData->SetDecoderStatus(decoderStatus);
      //reset spider
      if (fSpider)
	ResetSpider();
      break;
      
    case CHAIN_A_TRAILER:

      //try to recover error
      if (fDecoderSummaryData->GetRecoveringError())
	continue;
      //check decode status
      if ( decoderStatus != CHAIN_A_TRAILER_STATUS  && !fDecoderSummaryData->GetRecoveringError()){
	  if (fLogErrors)
	    AliError(Form("  %02x - 0x%08x [ERROR] Unexpected TRM chain A trailer (curslot=%d)",decoderStatus,*rawData,fDecoderSummaryData->GetCurrentSlotID()));
	fDecoderSummaryData->SetErrorDetected(kTRUE);
	fDecoderSummaryData->SetErrorSlotID(fDecoderSummaryData->GetCurrentSlotID());
	  errorWarning++;
	  //try to recover error
	  if (fRecoverError){
	    if (errorWarning > fRecoverErrorThr) {
	      if (fVerbose)
		AliInfo("Trying to recover the error: searching for the next header");
	      fDecoderSummaryData->SetRecoveringError(kTRUE);
	      continue;
	    }
	    else {
	      if (fVerbose)
		AliInfo("Do not try to recover error yet, go on with decoding process");
	      continue;
	    }
	  }
	return(fDecoderSummaryData->GetErrorDetected());
      }
      //decode status ok
	errorWarning = 0;
      //set TRM chain trailer
      fTRMChainTrailer = (AliTOFTRMChainTrailer *)rawData;
      //fill chain summary data
      FillChainSummaryData(fTRMChainTrailer);
      //print verbose
      if (fVerbose)
	AliInfo(Form("  %02x - 0x%08x \t  TRM chain A trailer",decoderStatus,*rawData));
      //change decode status
      decoderStatus = decoderStatus & ~CHAIN_A_BIT;
      fDecoderSummaryData->SetDecoderStatus(decoderStatus);
      break;
      
    case CHAIN_B_HEADER:

      //try to recover error
      if (fDecoderSummaryData->GetRecoveringError())
	continue;
      //check decode status
      if ( decoderStatus != CHAIN_B_HEADER_STATUS  && !fDecoderSummaryData->GetRecoveringError()){
	  if (fLogErrors)
	    AliError(Form("  %02x - 0x%08x [ERROR] Unexpected TRM chain B header (curslot=%d)",decoderStatus,*rawData,fDecoderSummaryData->GetCurrentSlotID()));
	fDecoderSummaryData->SetErrorDetected(kTRUE);
	fDecoderSummaryData->SetErrorSlotID(fDecoderSummaryData->GetCurrentSlotID());
	  errorWarning++;
	  //try to recover error
	  if (fRecoverError){
	    if (errorWarning > fRecoverErrorThr) {
	      if (fVerbose)
		AliInfo("Trying to recover the error: searching for the next header");
	      fDecoderSummaryData->SetRecoveringError(kTRUE);
	      continue;
	    }
	    else {
	      if (fVerbose)
		AliInfo("Do not try to recover error yet, go on with decoding process");
	      continue;
	    }
	  }
	return(fDecoderSummaryData->GetErrorDetected());
      }
      //decode status ok
	errorWarning = 0;
      //set TRM chain header
      fTRMChainHeader = (AliTOFTRMChainHeader *)rawData;
      //fill decoder summary data
      fDecoderSummaryData->SetCurrentChain(1);
      //get chain summary data
      fChainSummaryData = fTRMSummaryData->GetChainSummaryData(1);
      //reset chain summary data
      fChainSummaryData->Reset();
      //fill chain summary data
      FillChainSummaryData(fTRMChainHeader);
      //get tdc hit buffer
      fTDCHitBuffer = fChainSummaryData->GetTDCHitBuffer();
      //reset tdc hit buffer
      fTDCHitBuffer->Reset();
      //get tdc packed hit buffer
      fTDCPackedHitBuffer = fChainSummaryData->GetTDCPackedHitBuffer();
      //reset tdc packed hit buffer
      fTDCPackedHitBuffer->Reset();
      //get tdc error buffer
      fTDCErrorBuffer = fChainSummaryData->GetTDCErrorBuffer();
      //reset tdc error buffer
      fTDCErrorBuffer->Reset();
      //print verbose
      if (fVerbose)
	AliInfo(Form("  %02x - 0x%08x \t  TRM chain B header \t chain=%01d bunchID=%04d",decoderStatus,*rawData,1,fTRMChainHeader->GetBunchID()));
      //change decode status
      decoderStatus = decoderStatus | CHAIN_B_BIT;
      fDecoderSummaryData->SetDecoderStatus(decoderStatus);
      //reset spider
      if (fSpider)
	ResetSpider();
      break;
      
    case CHAIN_B_TRAILER:

      //try to recover error
      if (fDecoderSummaryData->GetRecoveringError())
	continue;
      //check decode status
      if ( decoderStatus != CHAIN_B_TRAILER_STATUS  && !fDecoderSummaryData->GetRecoveringError()){
	  if (fLogErrors)
	    AliError(Form("  %02x - 0x%08x [ERROR] Unexpected TRM chain B trailer (curslot=%d)",decoderStatus,*rawData,fDecoderSummaryData->GetCurrentSlotID()));
	fDecoderSummaryData->SetErrorDetected(kTRUE);
	fDecoderSummaryData->SetErrorSlotID(fDecoderSummaryData->GetCurrentSlotID());
	  errorWarning++;
	  //try to recover error
	  if (fRecoverError){
	    if (errorWarning > fRecoverErrorThr) {
	      if (fVerbose)
		AliInfo("Trying to recover the error: searching for the next header");
	      fDecoderSummaryData->SetRecoveringError(kTRUE);
	      continue;
	    }
	    else {
	      if (fVerbose)
		AliInfo("Do not try to recover error yet, go on with decoding process");
	      continue;
	    }
	  }
	return(fDecoderSummaryData->GetErrorDetected());
      }
      //decode status ok
	errorWarning = 0;
      //set TRM chain trailer
      fTRMChainTrailer = (AliTOFTRMChainTrailer *)rawData;
      //fill chain summary data
      FillChainSummaryData(fTRMChainTrailer);
      //print verbose
      if (fVerbose)
	AliInfo(Form("  %02x - 0x%08x \t  TRM chain B trailer",decoderStatus,*rawData));
      //change decode status
      decoderStatus = decoderStatus & ~CHAIN_B_BIT;
      fDecoderSummaryData->SetDecoderStatus(decoderStatus);
      break;
      
    case ERROR:

      //try to recover error
      if (fDecoderSummaryData->GetRecoveringError())
	continue;
      //decode TRM TDC error
      fTRMTDCError = (AliTOFTRMTDCError *)rawData;
      //check diagnostic word
      if (fTRMTDCError->GetTDCID() == 15) {
	if (fVerbose)
	  AliInfo(Form("  %02x - 0x%08x \t  Diagnostic error word",decoderStatus,*rawData));
	break;
      }
      //set error data
      error.SetErrorFlags(fTRMTDCError->GetErrorFlags());
      error.SetTDCID(fTRMTDCError->GetTDCID());
      //fill TDC error buffer
      fTDCErrorBuffer->Add(error);
      if (fVerbose)
	AliInfo(Form("  %02x - 0x%08x \t  TDC error",decoderStatus,*rawData));
      break;
      
    case FILLER:

      //try to recover error
      if (fDecoderSummaryData->GetRecoveringError())
	continue;
      if (fVerbose)
	AliInfo(Form("  %02x - 0x%08x \t  Filler",decoderStatus,*rawData));
      break;
      
    default:

      //try to recover error
      if (fDecoderSummaryData->GetRecoveringError())
	continue;
      //check decode status
      if ( decoderStatus != CHAIN_A_TDC_HIT_STATUS &&
	   decoderStatus != CHAIN_B_TDC_HIT_STATUS  && !fDecoderSummaryData->GetRecoveringError()){
	  if (fLogErrors)
	    AliError(Form("  %02x - 0x%08x [ERROR] Unexpected or unknown word (curslot=%d)",decoderStatus,*rawData,fDecoderSummaryData->GetCurrentSlotID()));
	fDecoderSummaryData->SetErrorDetected(kTRUE);
	fDecoderSummaryData->SetErrorSlotID(fDecoderSummaryData->GetCurrentSlotID());
	  errorWarning++;
	  //try to recover error
	  if (fRecoverError){
	    if (errorWarning > fRecoverErrorThr) {
	      if (fVerbose)
		AliInfo("Trying to recover the error: searching for the next header");
	      fDecoderSummaryData->SetRecoveringError(kTRUE);
	      continue;
	    }
	    else {
	      if (fVerbose)
		AliInfo("Do not try to recover error yet, go on with decoding process");
	      continue;
	    }
	  }
	return(fDecoderSummaryData->GetErrorDetected());
      }
      //decode status ok
	errorWarning = 0;
      
      //switch TRM ACQ
      switch (fTRMSummaryData->GetACQBits()){
	
      case PACKING_ENABLED_ACQ:
	//decode TDC packed/unpacked hit
	fTDCPackedHit = (AliTOFTDCPackedHit *)rawData;
	fTDCUnpackedHit = (AliTOFTDCUnpackedHit *)rawData;
	//set hit data
	hit.SetChan(fTDCUnpackedHit->GetChan());
	hit.SetTDCID(fTDCUnpackedHit->GetTDCID());
	hit.SetEBit(fTDCUnpackedHit->GetEBit());
	hit.SetPSBits(fTDCUnpackedHit->GetPSBits());
	//switch PS bits
	switch (hit.GetPSBits()){
	  //packed hit or overflow hit
	case PACKED_HIT_PS: case TOT_OVF_HIT_PS:
	  hit.SetHitTime(fTDCPackedHit->GetHitTime());
	  hit.SetTOTWidth(fTDCPackedHit->GetTOTWidth());
	  //add hit
	  fTDCHitBuffer->Add(hit);
	  fTDCPackedHitBuffer->Add(hit);
	  break; 
	  //orphane leading
	case LEADING_HIT_PS:
	  hit.SetHitTime(fTDCUnpackedHit->GetHitTime());
	  hit.SetTOTWidth(0);
	  //add hit
	  fTDCHitBuffer->Add(hit);
	  fTDCPackedHitBuffer->Add(hit);
	  break;
	  //orphane trailing
	case TRAILING_HIT_PS:
	  hit.SetHitTime(fTDCUnpackedHit->GetHitTime());
	  hit.SetTOTWidth(0);
	  //add hit
	  fTDCHitBuffer->Add(hit);
	  break;
	}
	//end switch PS bits
	//print verbose
	if (fVerbose)
	  switch (hit.GetPSBits()){
	  case PACKED_HIT_PS:
	    AliInfo(Form("  %02x - 0x%08x \t  TDC hit [packed] \t PS=%1d TDC=%1d chan=%1d TOT=%3d time=%4d (%7.1f ns)",decoderStatus,*rawData,hit.GetPSBits(),hit.GetTDCID(),hit.GetChan(),hit.GetTOTWidth(),hit.GetHitTime(),hit.GetHitTime()*TIME_BIN_WIDTH));
	    break;
	  case LEADING_HIT_PS:
	    AliInfo(Form("  %02x - 0x%08x \t  TDC hit [orp.lead] \t PS=%1d TDC=%1d chan=%1d time=%4d (%7.1f ns)",decoderStatus,*rawData,hit.GetPSBits(),hit.GetTDCID(),hit.GetChan(),hit.GetHitTime(),hit.GetHitTime()*TIME_BIN_WIDTH));
	    break;
	  case TRAILING_HIT_PS:
	    AliInfo(Form("  %02x - 0x%08x \t  TDC hit [orp.trai] \t PS=%1d TDC=%1d chan=%1d time=%4d (%7.1f ns)",decoderStatus,*rawData,hit.GetPSBits(),hit.GetTDCID(),hit.GetChan(),hit.GetHitTime(),hit.GetHitTime()*TIME_BIN_WIDTH));
	    break;
	  case TOT_OVF_HIT_PS:
	    AliInfo(Form("  %02x - 0x%08x \t  TDC hit [TOT ovfl] \t PS=%1d TDC=%1d chan=%1d TOT=%3d time=%4d (%7.1f ns)",decoderStatus,*rawData,hit.GetPSBits(),hit.GetTDCID(),hit.GetChan(),hit.GetTOTWidth(),hit.GetHitTime(),hit.GetHitTime()*TIME_BIN_WIDTH));
	    break;
	  }
	break;
	
      case LEADING_ONLY_ACQ: case TRAILING_ONLY_ACQ:
	//decode TDC unpacked hit
	fTDCUnpackedHit = (AliTOFTDCUnpackedHit *)rawData;
	//set hit data
	hit.SetChan(fTDCUnpackedHit->GetChan());
	hit.SetTDCID(fTDCUnpackedHit->GetTDCID());
	hit.SetEBit(fTDCUnpackedHit->GetEBit());
	hit.SetPSBits(fTDCUnpackedHit->GetPSBits());
	hit.SetHitTime(fTDCUnpackedHit->GetHitTime());
	hit.SetTOTWidth(0);
	//add hit
	fTDCHitBuffer->Add(hit);
	//print verbose
	if (fVerbose)
	  switch (hit.GetPSBits()){
	  case LEADING_HIT_PS:
	    AliInfo(Form("  %02x - 0x%08x \t  TDC hit [leading] \t PS=%1d TDC=%1d chan=%1d time=%4d (%7.1f ns)",decoderStatus,*rawData,hit.GetPSBits(),hit.GetTDCID(),hit.GetChan(),hit.GetHitTime(),hit.GetHitTime()*TIME_BIN_WIDTH));
	    break;
	  case TRAILING_HIT_PS:
	    AliInfo(Form("  %02x - 0x%08x \t  TDC hit [trailing] \t PS=%1d TDC=%1d chan=%1d time=%4d (%7.1f ns)",decoderStatus,*rawData,hit.GetPSBits(),hit.GetTDCID(),hit.GetChan(),hit.GetHitTime(),hit.GetHitTime()*TIME_BIN_WIDTH));
	  }
	break;
	
      case PACKING_DISABLED_ACQ:
	//decode TDC unpacked hit
	fTDCUnpackedHit = (AliTOFTDCUnpackedHit *)rawData;
	//set hit data
	hit.SetChan(fTDCUnpackedHit->GetChan());
	hit.SetTDCID(fTDCUnpackedHit->GetTDCID());
	hit.SetEBit(fTDCUnpackedHit->GetEBit());
	hit.SetPSBits(fTDCUnpackedHit->GetPSBits());
	hit.SetHitTime(fTDCUnpackedHit->GetHitTime());
	hit.SetTOTWidth(0);
	//add hit
	fTDCHitBuffer->Add(hit);
	//print verbose
	if (fVerbose)
	  switch (hit.GetPSBits()){
	  case LEADING_HIT_PS:
	    AliInfo(Form("  %02x - 0x%08x \t  TDC hit [leading] \t PS=%1d TDC=%1d chan=%1d time=%4d (%7.1f ns)",decoderStatus,*rawData,hit.GetPSBits(),hit.GetTDCID(),hit.GetChan(),hit.GetHitTime(),hit.GetHitTime()*TIME_BIN_WIDTH));
	    break;
	  case TRAILING_HIT_PS:
	    AliInfo(Form("  %02x - 0x%08x \t  TDC hit [trailing] \t PS=%1d TDC=%1d chan=%1d time=%4d (%7.1f ns)",decoderStatus,*rawData,hit.GetPSBits(),hit.GetTDCID(),hit.GetChan(),hit.GetHitTime(),hit.GetHitTime()*TIME_BIN_WIDTH));
	  } //print verbose
	//spider
	if (fSpider)
	  Spider(hit);

	break;
      }
      //end switch TRM ACQ

      break;
      
    }
    
  }
  //end equipment data loop
  
  if (fVerbose)
    AliInfo("End of data loop");

  //reset spider
  if (fSpider)
    ResetSpider();
  
  /*** V2718 patch ***/
  if (fV2718Patch){
    decoderStatus = decoderStatus & ~DRM_BIT;
    fDecoderSummaryData->SetDecoderStatus(decoderStatus);
    fDRMSummaryData->SetTrailer(kTRUE);
    fDRMSummaryData->SetSlotEnableMask(fDRMSummaryData->GetDecoderSlotEnableMask());
    fDRMSummaryData->SetCBit(1);
    fDRMSummaryData->SetLocalEventCounter(fTRMSummaryData->GetEventCounter());
    if (fVerbose)
      AliInfo("DRM was not present: - V2718 end patch decoding -");
  }
  /*** V2718 patch ***/
 
  if (fVerbose)
    AliInfo("Decoder is exiting succesfully.");

  return(fDecoderSummaryData->GetErrorDetected());  
}

//_________________________________________________________________

void
AliTOFDecoderV2::FillDRMSummaryData(AliTOFDRMGlobalHeader *DRMGlobalHeader)
{
  fDRMSummaryData->SetHeader(kTRUE);
  fDRMSummaryData->SetSlotID(DRMGlobalHeader->GetSlotID());
  fDRMSummaryData->SetEventWords(DRMGlobalHeader->GetEventWords());
  fDRMSummaryData->SetDRMID(DRMGlobalHeader->GetDRMID());
}

//_________________________________________________________________

void
AliTOFDecoderV2::FillDRMSummaryData(AliTOFDRMGlobalTrailer *DRMGlobalTrailer)
{
  fDRMSummaryData->SetTrailer(kTRUE);
  fDRMSummaryData->SetLocalEventCounter(DRMGlobalTrailer->GetLocalEventCounter());
}

//_________________________________________________________________

void
AliTOFDecoderV2::FillDRMSummaryData(AliTOFDRMStatusHeader1 *DRMStatusHeader1)
{
  fDRMSummaryData->SetPartecipatingSlotID(DRMStatusHeader1->GetPartecipatingSlotID());
  fDRMSummaryData->SetCBit(DRMStatusHeader1->GetCBit());
  fDRMSummaryData->SetVersID(DRMStatusHeader1->GetVersID());
  fDRMSummaryData->SetDRMhSize(DRMStatusHeader1->GetDRMhSize());
}

//_________________________________________________________________

void
AliTOFDecoderV2::FillDRMSummaryData(AliTOFDRMStatusHeader2 *DRMStatusHeader2)
{
  fDRMSummaryData->SetSlotEnableMask(DRMStatusHeader2->GetSlotEnableMask());
  fDRMSummaryData->SetFaultID(DRMStatusHeader2->GetFaultID());
  fDRMSummaryData->SetRTOBit(DRMStatusHeader2->GetRTOBit());
}

//_________________________________________________________________

void
AliTOFDecoderV2::FillDRMSummaryData(AliTOFDRMStatusHeader3 *DRMStatusHeader3)
{
  fDRMSummaryData->SetL0BCID(DRMStatusHeader3->GetL0BCID());
  fDRMSummaryData->SetRunTimeInfo(DRMStatusHeader3->GetRunTimeInfo());
}

//_________________________________________________________________

void
AliTOFDecoderV2::FillDRMSummaryData(AliTOFDRMStatusHeader4 *DRMStatusHeader4)
{
  fDRMSummaryData->SetTemperature(DRMStatusHeader4->GetTemperature());
  fDRMSummaryData->SetACKBit(DRMStatusHeader4->GetACKBit());
  fDRMSummaryData->SetSensAD(DRMStatusHeader4->GetSensAD());
}

//_________________________________________________________________

void
AliTOFDecoderV2::FillDRMSummaryData(AliTOFDRMEventCRC *DRMEventCRC)
{
  fDRMSummaryData->SetEventCRC(DRMEventCRC->GetEventCRC());
}

//_________________________________________________________________

void
AliTOFDecoderV2::FillLTMSummaryData(AliTOFLTMGlobalHeader *LTMGlobalHeader)
{
  fLTMSummaryData->SetHeader(kTRUE);
  fLTMSummaryData->SetSlotID(LTMGlobalHeader->GetSlotID());
  fLTMSummaryData->SetEventWords(LTMGlobalHeader->GetEventWords());
  fLTMSummaryData->SetCBit(LTMGlobalHeader->GetCBit());
  fLTMSummaryData->SetFault(LTMGlobalHeader->GetFault());
}

//_________________________________________________________________

void
AliTOFDecoderV2::FillLTMSummaryData(AliTOFLTMGlobalTrailer *LTMGlobalTrailer)
{
  fLTMSummaryData->SetTrailer(kTRUE);
  fLTMSummaryData->SetEventCRC(LTMGlobalTrailer->GetEventCRC());
  fLTMSummaryData->SetEventNumber(LTMGlobalTrailer->GetEventNumber());
}

//_________________________________________________________________

void
AliTOFDecoderV2::FillLTMSummaryData(AliTOFLTMPDLData *LTMPDLData, Int_t PDLWord)
{
  fLTMSummaryData->SetPDL(4 * PDLWord + 0, LTMPDLData->GetPDLValue1());
  fLTMSummaryData->SetPDL(4 * PDLWord + 1, LTMPDLData->GetPDLValue2());
  fLTMSummaryData->SetPDL(4 * PDLWord + 2, LTMPDLData->GetPDLValue3());
  fLTMSummaryData->SetPDL(4 * PDLWord + 3, LTMPDLData->GetPDLValue4());
}

//_________________________________________________________________

void
AliTOFDecoderV2::FillLTMSummaryData(AliTOFLTMADCData *LTMADCData, Int_t ADCWord)
{
  fLTMSummaryData->SetADC(3 * ADCWord + 0, LTMADCData->GetADCValue1());
  fLTMSummaryData->SetADC(3 * ADCWord + 1, LTMADCData->GetADCValue2());
  fLTMSummaryData->SetADC(3 * ADCWord + 2, LTMADCData->GetADCValue3());
}

//_________________________________________________________________

void
AliTOFDecoderV2::FillLTMSummaryData(AliTOFLTMORData *LTMORData, Int_t ORWord)
{
  fLTMSummaryData->SetOR(3 * ORWord + 0, LTMORData->GetORValue1());
  fLTMSummaryData->SetOR(3 * ORWord + 1, LTMORData->GetORValue2());
  fLTMSummaryData->SetOR(3 * ORWord + 2, LTMORData->GetORValue3());
}

//_________________________________________________________________

void
AliTOFDecoderV2::FillTRMSummaryData(AliTOFTRMGlobalHeader *TRMGlobalHeader)
{
  fTRMSummaryData->SetHeader(kTRUE);
  fTRMSummaryData->SetSlotID(TRMGlobalHeader->GetSlotID());
  fTRMSummaryData->SetEventWords(TRMGlobalHeader->GetEventWords());
  fTRMSummaryData->SetACQBits(TRMGlobalHeader->GetACQBits());
  fTRMSummaryData->SetLBit(TRMGlobalHeader->GetLBit());
  fTRMSummaryData->SetEBit(TRMGlobalHeader->GetEBit());
}

//_________________________________________________________________

void
AliTOFDecoderV2::FillTRMSummaryData(AliTOFTRMGlobalTrailer *TRMGlobalTrailer)
{
  fTRMSummaryData->SetTrailer(kTRUE);
  fTRMSummaryData->SetEventCRC(TRMGlobalTrailer->GetEventCRC());
  fTRMSummaryData->SetEventCounter(TRMGlobalTrailer->GetEventCounter());
}

//_________________________________________________________________

void
AliTOFDecoderV2::FillChainSummaryData(AliTOFTRMChainHeader *TRMChainHeader)
{
  fChainSummaryData->SetHeader(kTRUE);
  switch (*(UInt_t *)TRMChainHeader & WORD_TYPE_MASK){
  case CHAIN_A_HEADER:
    fChainSummaryData->SetChain(0);
    break;
  case CHAIN_B_HEADER:
    fChainSummaryData->SetChain(1);
    break;
  }
  fChainSummaryData->SetBunchID(TRMChainHeader->GetBunchID());
  fChainSummaryData->SetPB24Temp(TRMChainHeader->GetPB24Temp());
  fChainSummaryData->SetPB24ID(TRMChainHeader->GetPB24ID());
  fChainSummaryData->SetTSBit(TRMChainHeader->GetTSBit());
}

//_________________________________________________________________

void
AliTOFDecoderV2::FillChainSummaryData(AliTOFTRMChainTrailer *TRMChainTrailer)
{
  fChainSummaryData->SetTrailer(kTRUE);
  fChainSummaryData->SetStatus(TRMChainTrailer->GetStatus());
  fChainSummaryData->SetEventCounter(TRMChainTrailer->GetEventCounter());
}

//_________________________________________________________________

void
AliTOFDecoderV2::ResetSpider()
{		      
  //reset condition
  if (fVerbose)
    AliInfo("Reset signal received, empty and reset buffer");
  for (Int_t iChan = 0; iChan < N_CHANNEL; iChan++){
    if (fSpiderBufferFull[iChan]) {
      if (fVerbose)
	AliInfo(Form("Spider buffer is full for channel %d", iChan));
      fSpiderTDCPackedHitBuffer->Add(fSpiderBuffer[iChan]);
    }
    fSpiderBufferFull[iChan] = kFALSE;
  }
  fSpiderTDCID = -1;
  return;
}

void
AliTOFDecoderV2::Spider(AliTOFTDCHit &hit){
  
  if (fVerbose)
    AliInfo("Hit has been received from decode main routine");

  //check new TDC
  if (fSpiderTDCID != hit.GetTDCID()){
    if (fVerbose)
      AliInfo("Data coming from a new TDC, empty and reset buffer");
    for (Int_t iChan = 0; iChan < N_CHANNEL; iChan++){
      if (fSpiderBufferFull[iChan])
	fSpiderTDCPackedHitBuffer->Add(fSpiderBuffer[iChan]);
      fSpiderBufferFull[iChan] = kFALSE;
    }
    fSpiderTDCPackedHitBuffer = fTDCPackedHitBuffer;
    fSpiderTDCID = hit.GetTDCID();
  }	      

  //switch PS bits
  switch(hit.GetPSBits()){
    //leading hit
  case LEADING_HIT_PS:
    //check buffer status
    if (fSpiderBufferFull[hit.GetChan()]){ //buffer full
      fSpiderTDCPackedHitBuffer->Add(fSpiderBuffer[hit.GetChan()]); //buffered hit is orphane
      fSpiderBuffer[hit.GetChan()] = hit; //current hit into buffer
      if (fVerbose)
	AliInfo("Leading hit and buffer full, buffered hit is a orphane leading hit");
    } 
    else{ //buffer empty
      fSpiderBuffer[hit.GetChan()] = hit; //current hit into buffer
      fSpiderBufferFull[hit.GetChan()] = kTRUE; //set buffer full
    }
    break;
    //trailing hit
  case TRAILING_HIT_PS:
    //check buffer status
    if (fSpiderBufferFull[hit.GetChan()]){ //buffer full
      fSpiderTDCPackedHitBuffer->Add(fSpiderBuffer[hit.GetChan()] << hit); //pack hits (Leading << Trailing) and save
      fSpiderBufferFull[hit.GetChan()] = kFALSE; //unset buffer full
      if (fVerbose)
	AliInfo("Trailing hit and buffer full, pack leading and trailing hit");
    } 
    else{ //buffer empty
      ; //do nothing
      if (fVerbose)
	AliInfo("Trailing hit and buffer empty, trow trailing hit away");
    }
    break;
  } //switch PS bits

}

//_________________________________________________________________

Bool_t 
AliTOFDecoderV2::DecodeNext()
{
  /* decode next */

  if (!fRawReader || !fRawReader->ReadHeader())
    return kFALSE;

  const Int_t size = fRawReader->GetDataSize(); 
  UChar_t *data = new UChar_t[size];
  if (fRawReader->ReadNext(data, size) != 1) {
    delete [] data;
    return kFALSE;
  }
      
  /* decode equipment data */
  SetEquipmentID(fRawReader->GetEquipmentId());
  Decode((UInt_t *)data, size / 4);

  delete [] data;
  return kTRUE;
}
