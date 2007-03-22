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
 **************************************************************************/

/* $Id$ */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// This class provides access to TRD digits in raw data.                     //
//                                                                           //
// It loops over all TRD digits in the raw data given by the AliRawReader.   //
// The Next method goes to the next digit. If there are no digits left       //
// it returns kFALSE.                                                        //
// Several getters provide information about the current digit.              //
//                                                                           //
// Author: C. Lippmann (C.Lippmann@gsi.de)                                   //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "AliLog.h"
#include "AliRawReader.h"
#include "AliTRDdigitsManager.h"
#include "AliTRDdataArrayI.h"

#include "AliTRDRawStream.h"
#include "AliTRDgeometry.h"
#include "AliTRDCommonParam.h"
#include "AliTRDcalibDB.h"

ClassImp(AliTRDRawStream)

//_____________________________________________________________________________
AliTRDRawStream::AliTRDRawStream() 
  :TObject()
  ,fSig()
  ,fADC(0)
  ,fTB(0)
  ,fEv(0)
  ,fROB(0)
  ,fMCM(0)
  ,fSM(0)
  ,fLAYER(0)
  ,fSTACK(0)
  ,fROC(0)
  ,fSIDE(0)
  ,fDCS(0)
  ,fROW(0)
  ,fCOL(0)
  ,fBCctr(0)
  ,fPTctr(0)
  ,fPTphase(0)
  ,fRVmajor(0)
  ,fRVminor(0)
  ,fHCHWords(0)
  ,fTBins(0)
  ,fTCon(0)
  ,fPEDon(0)
  ,fGAINon(0)
  ,fXTon(0)
  ,fNonLinOn(0)
  ,fBypass(0)
  ,fCommonAdditive(0)
  ,fZeroSuppressed(0)
  ,fHCHctr1(0)
  ,fHCHctr2(0)
  ,fMCMHctr1(0)
  ,fMCMHctr2(0)
  ,fGTUctr1(0)
  ,fGTUctr2(0)
  ,fHCdataCtr(0)
  ,fTracklPID(0.)
  ,fTracklDefL(0.)
  ,fTracklPadPos(0.)
  ,fTracklPadRow(0)
  ,fGTUlinkMask()
  ,fRawReader(NULL)
  ,fCount(0)
  ,fDetector(-1)
  ,fPrevDetector(-1)
  ,fNPads(-1)
  ,fRow(-1)
  ,fPrevRow(-1)
  ,fColumn(-1)
  ,fPrevColumn(-1)
  ,fTime(-1)
  ,fSignal(-1)
  ,fRawVersion(2)
  ,fDataWord(0)
  ,fStatus(0)
  ,fTbSwitch(0)
  ,fTbSwitchCtr(0)
  ,fTimeWords(0)
  ,fWordCtr(0)
  ,fRowMax(0)
  ,fColMax(0)
  ,fADCmask()
  ,fChamberDone()
  ,fGeo(NULL)
  ,fDigitsManager(NULL) 
  ,fDigits(NULL) 
  ,fTrack0(NULL) 
  ,fTrack1(NULL) 
  ,fTrack2(NULL) 
{
  //
  // Default constructor
  //

  for (Int_t i = 0; i < 540; i++) {
    fChamberDone[i] = 0;
  }

}

//_____________________________________________________________________________
AliTRDRawStream::AliTRDRawStream(AliRawReader *rawReader
                               , AliTRDdigitsManager *man
			       , AliTRDdataArrayI *dig) 
  :TObject()
  ,fSig()
  ,fADC(0)
  ,fTB(0)
  ,fEv(0)
  ,fROB(0)
  ,fMCM(0)
  ,fSM(0)
  ,fLAYER(0)
  ,fSTACK(0)
  ,fROC(0)
  ,fSIDE(0)
  ,fDCS(0)
  ,fROW(0)
  ,fCOL(0)
  ,fBCctr(0)
  ,fPTctr(0)
  ,fPTphase(0)
  ,fRVmajor(0)
  ,fRVminor(0)
  ,fHCHWords(0)
  ,fTBins(0)
  ,fTCon(0)
  ,fPEDon(0)
  ,fGAINon(0)
  ,fXTon(0)
  ,fNonLinOn(0)
  ,fBypass(0)
  ,fCommonAdditive(0)
  ,fZeroSuppressed(0)
  ,fHCHctr1(0)
  ,fHCHctr2(0)
  ,fMCMHctr1(0)
  ,fMCMHctr2(0)
  ,fGTUctr1(0)
  ,fGTUctr2(0)
  ,fHCdataCtr(0)
  ,fTracklPID(0.)
  ,fTracklDefL(0.)
  ,fTracklPadPos(0.)
  ,fTracklPadRow(0)
  ,fGTUlinkMask()
  ,fRawReader(rawReader)
  ,fCount(0)
  ,fDetector(-1)
  ,fPrevDetector(-1)
  ,fNPads(-1)
  ,fRow(-1)
  ,fPrevRow(-1)
  ,fColumn(-1)
  ,fPrevColumn(-1)
  ,fTime(-1)
  ,fSignal(-1)
  ,fRawVersion(2)
  ,fDataWord(0)
  ,fStatus(0)
  ,fTbSwitch(0)
  ,fTbSwitchCtr(0)
  ,fTimeWords(0)
  ,fWordCtr(0)
  ,fRowMax(0)
  ,fColMax(0)
  ,fADCmask()
  ,fChamberDone()
  ,fGeo(NULL) 
  ,fDigitsManager(man) 
  ,fDigits(dig) 
  ,fTrack0(NULL) 
  ,fTrack1(NULL) 
  ,fTrack2(NULL) 

{
  //
  // Create an object to read TRD raw digits
  //

  fGeo = new AliTRDgeometry();

  fRawReader->Select("TRD");

  for (Int_t i = 0; i < 540; i++) {
    fChamberDone[i] = 0;
  }

}

//_____________________________________________________________________________
AliTRDRawStream::AliTRDRawStream(AliRawReader *rawReader) 
  :TObject()
  ,fSig()
  ,fADC(0)
  ,fTB(0)
  ,fEv(0)
  ,fROB(0)
  ,fMCM(0)
  ,fSM(0)
  ,fLAYER(0)
  ,fSTACK(0)
  ,fROC(0)
  ,fSIDE(0)
  ,fDCS(0)
  ,fROW(0)
  ,fCOL(0)
  ,fBCctr(0)
  ,fPTctr(0)
  ,fPTphase(0)
  ,fRVmajor(0)
  ,fRVminor(0)
  ,fHCHWords(0)
  ,fTBins(0)
  ,fTCon(0)
  ,fPEDon(0)
  ,fGAINon(0)
  ,fXTon(0)
  ,fNonLinOn(0)
  ,fBypass(0)
  ,fCommonAdditive(0)
  ,fZeroSuppressed(0)
  ,fHCHctr1(0)
  ,fHCHctr2(0)
  ,fMCMHctr1(0)
  ,fMCMHctr2(0)
  ,fGTUctr1(0)
  ,fGTUctr2(0)
  ,fHCdataCtr(0)
  ,fTracklPID(0.)
  ,fTracklDefL(0.)
  ,fTracklPadPos(0.)
  ,fTracklPadRow(0)
  ,fGTUlinkMask()
  ,fRawReader(rawReader)
  ,fCount(0)
  ,fDetector(-1)
  ,fPrevDetector(-1)
  ,fNPads(-1)
  ,fRow(-1)
  ,fPrevRow(-1)
  ,fColumn(-1)
  ,fPrevColumn(-1)
  ,fTime(-1)
  ,fSignal(-1)
  ,fRawVersion(2)
  ,fDataWord(0)
  ,fStatus(0)
  ,fTbSwitch(0)
  ,fTbSwitchCtr(0)
  ,fTimeWords(0)
  ,fWordCtr(0)
  ,fRowMax(0)
  ,fColMax(0)
  ,fADCmask()
  ,fChamberDone()
  ,fGeo(NULL) 
  ,fDigitsManager(NULL) 
  ,fDigits(NULL) 
  ,fTrack0(NULL) 
  ,fTrack1(NULL) 
  ,fTrack2(NULL) 

{
  //
  // Create an object to read TRD raw digits
  //

  fRawReader->Select("TRD");

  for (Int_t i = 0; i < 540; i++) {
    fChamberDone[i] = 0;
  }

}

//_____________________________________________________________________________
AliTRDRawStream::AliTRDRawStream(const AliTRDRawStream& stream)
  :TObject(stream)
  ,fSig()
  ,fADC(-1)
  ,fTB(-1)
  ,fEv(-1)
  ,fROB(-1)
  ,fMCM(-1)
  ,fSM(-1)
  ,fLAYER(-1)
  ,fSTACK(-1)
  ,fROC(-1)
  ,fSIDE(-1)
  ,fDCS(-1)
  ,fROW(-1)
  ,fCOL(-1)
  ,fBCctr(-1)
  ,fPTctr(-1)
  ,fPTphase(-1)
  ,fRVmajor(-1)
  ,fRVminor(-1)
  ,fHCHWords(-1)
  ,fTBins(-1)
  ,fTCon(0)
  ,fPEDon(0)
  ,fGAINon(0)
  ,fXTon(0)
  ,fNonLinOn(-1)
  ,fBypass(-1)
  ,fCommonAdditive(-1)
  ,fZeroSuppressed(0)
  ,fHCHctr1(-1)
  ,fHCHctr2(-1)
  ,fMCMHctr1(-1)
  ,fMCMHctr2(-1)
  ,fGTUctr1(-1)
  ,fGTUctr2(-1)
  ,fHCdataCtr(-1)
  ,fTracklPID(-1.)
  ,fTracklDefL(-1.)
  ,fTracklPadPos(-1.)
  ,fTracklPadRow(-1)
  ,fGTUlinkMask()
  ,fRawReader(NULL)
  ,fCount(-1)
  ,fDetector(-1)
  ,fPrevDetector(-1)
  ,fNPads(-1)
  ,fRow(-1)
  ,fPrevRow(-1)
  ,fColumn(-1)
  ,fPrevColumn(-1)
  ,fTime(-1)
  ,fSignal(-1)
  ,fRawVersion(-1)
  ,fDataWord(0)
  ,fStatus(-1)
  ,fTbSwitch(0)
  ,fTbSwitchCtr(0)
  ,fTimeWords(0)
  ,fWordCtr(0)
  ,fRowMax(-1)
  ,fColMax(-1)
  ,fADCmask()
  ,fChamberDone()
  ,fGeo(NULL) 
  ,fDigitsManager(NULL) 
  ,fDigits(NULL) 
  ,fTrack0(NULL) 
  ,fTrack1(NULL) 
  ,fTrack2(NULL) 

{
  //
  // Copy constructor
  //

  AliFatal("Copy constructor not implemented");

}

//_____________________________________________________________________________
AliTRDRawStream& AliTRDRawStream::operator = (const AliTRDRawStream& 
					      /* stream */)
{
  //
  // Assigment operator
  //

  Fatal("operator =", "assignment operator not implemented");
  return *this;

}

//_____________________________________________________________________________
AliTRDRawStream::~AliTRDRawStream()
{
  //
  // Destructor
  //

}

//_____________________________________________________________________________
Bool_t AliTRDRawStream::SetRawVersion(Int_t rv)
{
  //
  // Set the raw data version
  //

  if ( rv >= 0 && rv <= 3 ) {
    fRawVersion = rv;
    return kTRUE;
  }

  return kFALSE;

}

//_____________________________________________________________________________
Bool_t AliTRDRawStream::Next()
{
  //
  // This is Bogdans code for reading raw data (offline use only).
  // It is used for fRawVersion == 0. This funcyion read the next raw digit.
  // Returns kFALSE if there is no digit left
  //

  fPrevDetector = fDetector;
  fPrevRow      = fRow;
  fPrevColumn   = fColumn;
  UChar_t data;

  AliTRDcalibDB *calibration = AliTRDcalibDB::Instance();
  if (!calibration) return kFALSE;
  
  Int_t timeBins = calibration->GetNumberOfTimeBins();
  
  while (fCount >= 0) {

    while (fCount == 0) {  // next detector

      // read the flag
      if (!fRawReader->ReadNextChar(data)) {
        return kFALSE;
      }
      if (data != 0xBB) {
	AliError(Form("wrong flag: %x", data));
	fCount = -1;
	return kFALSE;
      }

      // read the detector number
      if (!fRawReader->ReadNextChar(data)) {
	AliError("Could not read detector number");
	fCount = -1;
	return kFALSE;
      }
      fDetector = data;
      if (!fRawReader->ReadNextChar(data)) {
	AliError("Could not read detector number");
	fCount = -1;
	return kFALSE;
      }
      fDetector += (UInt_t(data) << 8);

      // read the number of byts
      if (!fRawReader->ReadNextChar(data)) {
	AliError("Could not read number of bytes");
	fCount = -1;
	return kFALSE;
      }
      fCount = data;
      if (!fRawReader->ReadNextChar(data)) {
	AliError("Could not read number of bytes");
	fCount = -1;
	return kFALSE;
      }
      fCount += (UInt_t(data) << 8);
      if (!fRawReader->ReadNextChar(data)) {
        AliError("Could not read number of bytes");
        fCount = -1;
        return kFALSE;
      }
      fCount += (UInt_t(data) << 16);

      // read the number of active pads
      if (!fRawReader->ReadNextChar(data)) {
	AliError("Could not read number of active pads");
	fCount = -1;
	return kFALSE;
      }
      fNPads = data;
      if (!fRawReader->ReadNextChar(data)) {
	AliError("Could not read number of active pads");
	fCount = -1;
	return kFALSE;
      }
      fNPads += (UInt_t(data) << 8);

      fTime = timeBins;

    }

    // read the pad row and column number
    if ((fTime >= timeBins) && (fCount > 2)) {
      if (!fRawReader->ReadNextChar(data)) {
	AliError("Could not read row number");
	fCount = -1;
	return kFALSE;
      }
      fCount--;
      fRow = data - 1;
      if (!fRawReader->ReadNextChar(data)) {
	AliError("Could not read column number");
	fCount = -1;
	return kFALSE;
      }
      fCount--;
      fColumn = data - 1;
      fTime = 0;
    }

    // read the next data byte
    if (!fRawReader->ReadNextChar(data)) {
      AliError("Could not read data");
      fCount = -1;
      return kFALSE;
    }
    fCount--;

    if (data == 0) {  // zeros
      if (!fRawReader->ReadNextChar(data)) {
	AliError("Could not read time value");
	fCount = -1;
	return kFALSE;
      }
      fCount--;
      fTime += data + 1;

    } 
    else {          // signal
      fSignal = (UInt_t(data & 0x7F) << 8);
      if (!fRawReader->ReadNextChar(data)) {
	AliError("Could not read ADC value");
	fCount = -1;
	return kFALSE;
      }
      fCount--;
      fSignal += data;
      fTime++;
      return kTRUE;
    }
  }

  return kFALSE;

}

//____________________________________________________________________________
Int_t AliTRDRawStream::ReadAll()
{

  //
  // Read all TRD raw data word (32 bits). This is for all fRawVersion > 0.
  // Return 0 if something is not cool, 2 if data seems empty
  //
  // by C. Lippmann
  //
 
  AliTRDCommonParam *commonParam = AliTRDCommonParam::Instance();
  if (!commonParam) {
    AliError("Could not get common parameters");
    return 0;
  }

  AliTRDcalibDB     *calibration = AliTRDcalibDB::Instance();
  if (!calibration) {
    AliError("Could not get calibration object");
    return 0;
  }
  
  UInt_t timeTotal = calibration->GetNumberOfTimeBins();
  AliDebug(2, Form("Number of Timebins read from CDB: %d", timeTotal));

  // The number of data words needed for this number of time bins (there
  // are 3 time bins in one word)
  fTimeWords = (timeTotal-1)/3 + 1;

  fTbSwitch    = 3;
  fTbSwitchCtr = 0;

  Int_t  EqID     = 0;
  UInt_t datasize = 0;
  Int_t  iDET     = 0;

  Int_t retval = 0;

  Bool_t sizeOK = kFALSE;

  fHCHctr1 = fHCHctr2 =  0;
  fGTUctr1 = fGTUctr2 = -1;

  fHCdataCtr = 0;
  fWordCtr   = 0;

  AliInfo("Converting TRD raw data to digits ...");

  while ( 1 ) { // loop over all supermodules

    fWordCtr   = 0;
    fHCHctr1  = 0;
    fMCMHctr1 = 0;

    //
    // 0) Find first GTU Link Mask and test if we can read data
    //
    do {

      if ( !fRawReader->ReadNextInt( fDataWord ) ) {  // This is the standard exit point
	// Compress also the digits from the last detector
	if ( fChamberDone[iDET] == 2 ) {
	  //printf("Compressing data for det %d\n", iDET);
	  fDigits->Compress(1,0);
	  fTrack0->Compress(1,0);
	  fTrack1->Compress(1,0);
	  fTrack2->Compress(1,0);
	}
	AliInfo(Form("Finished processing TRD raw data: Found %d Half-Chambers", fHCHctr2));
	//
	/*
	fDigits = fDigitsManager->GetDigits(iDET+1);
	fTrack0 = fDigitsManager->GetDictionary(iDET+1,0);
	fTrack1 = fDigitsManager->GetDictionary(iDET+1,1);
	fTrack2 = fDigitsManager->GetDictionary(iDET+1,2);
	fDigits->Allocate(fRowMax,fColMax,timeTotal);
	fTrack0->Allocate(fRowMax,fColMax,timeTotal);
	fTrack1->Allocate(fRowMax,fColMax,timeTotal);
	fTrack2->Allocate(fRowMax,fColMax,timeTotal);
	fDigits->SetDataUnchecked(0, 0, 0, 50);
	fTrack0->SetDataUnchecked(0, 0, 0, 0);
	fTrack1->Set
	fTrack2->SetDataUnchecked(0, 0, 0, 0);	
	fDigits->Compress(1,0);
	fTrack0->Compress(1,0);
	fTrack1->Compress(1,0);
	fTrack2->Compress(1,0);
	*/
	//
	if ( sizeOK ) return 1;
	else          return 2;
      }

      fWordCtr++;

      // After reading the first word check for size of this data and get Eq. ID
      if ( fWordCtr == 1 ) {
	datasize = fRawReader->GetDataSize()/4;  // Size of this payload in 32bit words
	EqID     = fRawReader->GetEquipmentId(); // Get Equipment ID
	//if ( sizeOK = kFALSE && datasize > 0 ) { sizeOK = kTRUE; printf("Yo\n"); }
	if ( datasize > 0 ) sizeOK = kTRUE;
      }

      // GTU Link Mask?
      if ( (fDataWord & 0xfffff000) ==  0xe0000000 ) {
	fStatus = 1;      // GTU link mask found
	DecodeGTUlinkMask();
	break;
      } 
      else {
	AliError(Form("Equipment %d: First data word is not GTU Link Mask!", EqID));
	return 0;
      }

    } 
    while ( fWordCtr < datasize );

    //
    // loop over all half chambers in one supermodule
    //
    while ( fWordCtr < datasize ) {

      //
      // 1) Find end_of_tracklet_marker
      //
      while ( fWordCtr < datasize ) {

	if ( !fRawReader->ReadNextInt( fDataWord ) ) {
	  AliError("Could not read data");
	  return 0;
	}
	fWordCtr++;

	// GTU Link Mask?
	if ( (fDataWord & 0xfffff000) ==  0xe0000000 ) {
	  DecodeGTUlinkMask();
	  continue;
	}

	// endoftrackletmarker?
	if ( fDataWord == kEndoftrackletmarker ) {
	  AliDebug(3, "end-of-tracklet-marker found");
	  fStatus = 1;
	  break;
	} 
        else {
	  // Tracklets found
	  AliDebug(3, "Tracklet found");
	  DecodeTracklet();
	}

      }

      if ( fStatus == 0 ) break;
    
      //
      // 2) Look for non-end_of_tracklet_marker
      //
      fStatus = 0;
      while ( fWordCtr < datasize ) { 

	if ( !fRawReader->ReadNextInt( fDataWord ) ) {
	  AliError("Could not read data");
	  return 0;
	}
	fWordCtr++;
	//printf("Word %d: 0x%08x\n", fWordCtr, fDataWord); 

	if ( fDataWord != kEndoftrackletmarker ) {
	  fStatus = 1;
	  break;
	}

      }

      if ( fStatus == 0 ) break;
    
      //
      // 3) This Word must be Half Chamber Header
      //
      fStatus = 0;
      if ( (fDataWord & 0x00000003) == 1 ) { // HC header

	// If both half chambers of chamber corresponding to previous header
	// were already processed, we can compress these digits
	iDET = fGeo->GetDetector(fLAYER, fSTACK, fSM); // !!this is still the previous HC!!
	if ( fChamberDone[iDET] == 2 ) {
	  //printf("Compressing data for det %d\n", iDET);
	  fDigits->Compress(1,0);
	  fTrack0->Compress(1,0);
	  fTrack1->Compress(1,0);
	  fTrack2->Compress(1,0);
	}
	// Read from new HC header the chamber position (fLAYER, fSTACK, fSM)
	DecodeHCheader(timeTotal); // This is the new header!
	iDET    = fGeo->GetDetector(fLAYER, fSTACK, fSM);
	fRowMax = commonParam->GetRowMax(fLAYER,fSTACK,fSM);
	fColMax = commonParam->GetColMax(fROC);

	// The container for the digits of this detector
	fDigits = fDigitsManager->GetDigits(iDET);
	fTrack0 = fDigitsManager->GetDictionary(iDET,0);
	fTrack1 = fDigitsManager->GetDictionary(iDET,1);
	fTrack2 = fDigitsManager->GetDictionary(iDET,2);
	
	// Allocate memory if it was not already done
	if (fDigits->GetNtime() == 0) {
	  //printf("Allocating digits memory for det %d\n", iDET);
	  fDigits->Allocate(fRowMax, fColMax, fTBins);
	  fTrack0->Allocate(fRowMax, fColMax, fTBins);
	  fTrack1->Allocate(fRowMax, fColMax, fTBins);
	  fTrack2->Allocate(fRowMax, fColMax, fTBins);
	}
	if ( fZeroSuppressed ) {  // Make sure digits in this HC are 0
	  for ( Int_t colctr = fSIDE*fColMax/2; colctr < fColMax*(1+fSIDE)/2; colctr++ ) {
	    for ( Int_t rowctr = 0; rowctr < fRowMax; rowctr++ ) {
	      for ( Int_t timectr = 0; timectr < fTBins; timectr++ ) {
		fDigits->SetDataUnchecked(rowctr, colctr, timectr, 0);
		fTrack0->SetDataUnchecked(rowctr, colctr, timectr, 0);
		fTrack1->SetDataUnchecked(rowctr, colctr, timectr, 0);
		fTrack2->SetDataUnchecked(rowctr, colctr, timectr, 0);
	      }
	    }
	  }
	}
	fMCMHctr2 = 0;
	fHCdataCtr = 0;

	fChamberDone[iDET]++;

      }
    
      //
      // 4) Scan MCM data
      //
      fStatus = 0;
      while ( fWordCtr < datasize ) {

	if ( !fRawReader->ReadNextInt( fDataWord ) ) {
	  AliError("Could not read data");
	  return 0;
	}
	fWordCtr++;
	//printf("Word %d: 0x%08x\n", fWordCtr, fDataWord); 
      
	fHCdataCtr += 4;

	//if( (fDataWord & 0x0000000f) == 0xC ) { // MCM Header
	if( ((fDataWord & 0x80000000) == 0x0) && ((fDataWord & 0x0000000f) == 0xC) ) { // MCM Header
	  DecodeMCMheader();
	  if ( fMCM < 0 || fMCM > 15 || fROB < 0 || fROB > 7 ) {
	    AliError("Wrong fMCM or fROB. Skip this data");
	    break;
	  }
	  fTbSwitch    = 3;  // For first adc channel we expect: (fDataWord & 3) = 3
	  fTbSwitchCtr = 0;  // 
	  fADC = fTB   = 0;  // Reset Counter
	  fStatus      = 1;  // Now 1 means MCM header is found
	  continue;
	}
    
	if ( fDataWord == kEndofrawdatamarker ) {  // End of half-chamber data, finished
	  fGTUctr1 = -1;
	  break;
	}

	if ( fStatus == 1 ) {       // MCM header is set, ADC data is valid.
    
	  // Found some data. Decode it now:
	  retval = DecodeDataWord();
	  if ( retval ==  0 ) continue;
	  if ( retval == -1 ) break;
	  
	  if ( fADC > 1 && fADC < (Int_t)fGeo->ADCmax()-1 ) {

	    // Write Digits
	    if ( fCOL >= 0 && fCOL < fColMax && fROW >= 0 && fROW < fRowMax ) {  // A real pad
	      for ( Int_t ctr = 0; ctr < 3; ctr++ ) {
		if ( fTB+ctr < (Int_t)timeTotal ) {
		  /*
		  fDigits->SetDataUnchecked(fROW, fCOL, fTB+ctr, fSig[ctr]);
		  fTrack0->SetDataUnchecked(fROW, fCOL, fTB+ctr, 0);
		  fTrack1->SetDataUnchecked(fROW, fCOL, fTB+ctr, 0);
		  fTrack2->SetDataUnchecked(fROW, fCOL, fTB+ctr, 0);
		  */

		  //
		  // Pedestal subtraction done here (???):
		  //
		  //fDigits->SetData(fROW, fCOL, fTB+ctr, fSig[ctr] - fCommonAdditive);
		  fDigits->SetData(fROW, fCOL, fTB+ctr, fSig[ctr]);
		  fTrack0->SetData(fROW, fCOL, fTB+ctr, 0);
		  fTrack1->SetData(fROW, fCOL, fTB+ctr, 0);
		  fTrack2->SetData(fROW, fCOL, fTB+ctr, 0);
		}
	      }
	    }

	    fTB += 3;

	  }
          else {

            fCOL = -1;

	  }

	}

      }

      AliDebug(2, Form("SM%02d L%dS%d side %x: Processed %dMCMs=%dbyte", fSM, fLAYER, fSTACK, fSIDE+10,
		       fMCMHctr2, fHCdataCtr));

    } // End Half-Chamber loop

    AliDebug(1, Form("SM%02d (Eq %d): Processed %d HC (%dMCMs=%dbyte)", fSM, EqID, fHCHctr1, fMCMHctr1,
		     datasize*4));
    
  } // End Super Module loop

  // Compress also the digits from the last detector
  if ( fChamberDone[iDET] == 2 ) {
    //printf("Compressing data for det %d\n", iDET);
    fDigits->Compress(1,0);
    fTrack0->Compress(1,0);
    fTrack1->Compress(1,0);
    fTrack2->Compress(1,0);
  }

  if ( sizeOK ) return 1;
  else          return 2;

}

//============================================================================
// Decoding functions
//============================================================================

//____________________________________________________________________________
void AliTRDRawStream::DecodeHCheader(Int_t timeBins)
{
  //
  // Decode a half chamber header
  //

  if ( (fDataWord >> 31) == 0 )  { // Can only happen for fRawVersion == 1

    if ( fRawVersion != 1 ) {

      AliWarning("===============================================================================");
      AliWarning(Form("Mismatch between fRawVersion (%d) and HC header signature", fRawVersion));
      AliWarning("Setting fRawVersion to 1");
      AliWarning("===============================================================================");
      fRawVersion = 1;

    }

    DecodeHCheaderV1();
    return;

  } 
  else {

    fRVmajor = (fDataWord >> 24) & 0x7f;
    fRVminor = (fDataWord >> 17) & 0x7f;

    if ( fRawVersion != fRVmajor ) {

      AliWarning("===============================================================================");
      AliWarning(Form("Mismatch between fRawVersion (%d) and fRVmajor from HC header (%d)"
		      ,fRawVersion,fRVmajor));
      AliWarning(Form("Setting fRawVersion to %d", fRVmajor));
      AliWarning("===============================================================================");
      fRawVersion = fRVmajor;

    }
    if ( fRawVersion >= 2 && fRawVersion <= 3 ) {
      DecodeHCheaderV2V3(timeBins);
    }

    //
    // check for zero suppression
    if ( fRawVersion >= 3 || fRawVersion <= 4 ) fZeroSuppressed = kTRUE;
    else                                        fZeroSuppressed = kFALSE;

    return;

  }

  AliError(Form(" Unsupported raw version: %d", fRawVersion));
  return;

}

//____________________________________________________________________________
void AliTRDRawStream::DecodeHCheaderV1()
{

  //
  // Decode the HC header (fRawVersion == 1, SM I Commissioning 06)
  //

  if ( (fDataWord & 0x3) == 1 ) {

    fDCS   = (fDataWord >> 20);
    fSM    = (fDataWord >> 15) & 0x1f;
    fLAYER = (fDataWord >> 12) & 0x7;
    fSTACK = (fDataWord >>  9) & 0x7;
    fSIDE  = (fDataWord >>  8) & 0x1;

    fROC   = fGeo->GetDetectorSec(fLAYER, fSTACK);

    //AliDebug(3, Form("0x%08x: HC header: dcs=%d; sm=%d; roc=%d; side=%x", fDataWord, fDCS, fSM, fROC, fSIDE+10));

    if ((fSM    <  0) || 
        (fSM    > 17) || 
        (fLAYER <  0) || 
        (fLAYER >  5) || 
        (fSTACK <  0) || 
        (fSTACK >  4) || 
        (fSIDE  <  0) || 
        (fSIDE  >  1)) {
      AliWarning(Form("0x%08x: Strange HC header: dcs=%d; sm=%d; layer=%d; stack=%d.",
		      fDataWord, fDCS, fSM, fLAYER, fSTACK));
    } 
    else {
      fStatus = 1;
      fHCHctr1++;
      fHCHctr2++;
    }
    fHCHWords = 0;

  } 
  else { 

    AliError(Form("0x%08x: No HC header when it was expected.", fDataWord)); 

  }

}


//____________________________________________________________________________
void AliTRDRawStream::DecodeHCheaderV2V3(Int_t timeBins)
{
  //
  // Decode the HC header (fRawVersion == 2, 3, 4, ???)
  //

  // 1st word (h[0])
  if ( (fDataWord & 0x3) == 1 ) {

    fHCHWords = (fDataWord >> 14) & 0x7;
    fSM       = (fDataWord >>  9) & 0x1f;
    fLAYER    = (fDataWord >>  6) & 0x7;
    fSTACK    = (fDataWord >>  3) & 0x7;
    fSIDE     = (fDataWord >>  2) & 0x1;

    fROC      = fGeo->GetDetectorSec(fLAYER, fSTACK);

    AliDebug(3, Form("0x%08x: HC header: sm=%d; roc=%d; side=%x", fDataWord, fSM, fROC, fSIDE+10));

    if ((fSM    <  0) || 
        (fSM    > 17) || 
        (fLAYER <  0) || 
        (fLAYER >  5) || 
        (fSTACK <  0) || 
        (fSTACK >  4) || 
        (fSIDE  <  0) || 
        (fSIDE  >  1)) {
      AliError(Form("0x%08x: Strange HC header: dcs=%d; sm=%d; layer=%d; stack=%d.",
		    fDataWord, fDCS, fSM, fLAYER, fSTACK));
    } 
    else {
      fStatus = 1;
      fHCHctr1++;
      fHCHctr2++;
    }
  } 
  else { 
    AliError(Form("0x%08x: No HC header when it was expected.", fDataWord)); 
  }

  // 2nd word (h[1])
  if ( fHCHWords >= 1 ) {
    // read one more word
    if ( !fRawReader->ReadNextInt( fDataWord ) ) {
      AliError("Could not read data");
      return;
    }
    fWordCtr++;
    if ( (fDataWord & 0x3) == 1 ) {
      
      fBCctr   =  (fDataWord >> 16);
      fPTctr   =  (fDataWord >> 12) & 0xf;
      fPTphase =  (fDataWord >>  8) & 0xf;
      fTBins   = ((fDataWord >>  2) & 0x3f) + 1;

      AliDebug(3, Form("0x%08x: HC header 2: BCctr=%d PTctr=%d PTph=%d TB=%d"
                      , fDataWord, fBCctr, fPTctr, fPTphase, fTBins));

      if( fTBins != timeBins ) {

	AliWarning("===============================================================================");
	AliError(Form("Mismatch between nNTB from CDB (%d) and from HC header (%d)"
		      , timeBins, fTBins));
	AliWarning(Form("We will use the value from the raw data (HC header): %d", fTBins));
	AliWarning("===============================================================================");

      }

    }

  }

  // 3nd word (h[2])
  if ( fHCHWords >= 2 ) {
    // read one more word
    if ( !fRawReader->ReadNextInt( fDataWord ) ) {
      AliError("Could not read data");
      return;
    }
    fWordCtr++;
    if ( (fDataWord & 0x3) == 1 ) {
       
      fTCon     = (fDataWord >> 29) & 0x1;
      fPEDon    = (fDataWord >> 31) & 0x1;
      fGAINon   = (fDataWord >> 30) & 0x1;
      fXTon     = (fDataWord >> 28) & 0x1;
      fNonLinOn = (fDataWord >> 27) & 0x1;
      fBypass   = (fDataWord >> 26) & 0x1;

      fCommonAdditive = (fDataWord >> 20) & 0x3f;

      AliDebug(3, Form("0x%08x: HC header 3: TC=%d, PED=%d, GAIN=%d, XT=%d, NonLin=%d, Bypass=%d, Add=%d"
		      , fTCon, fPEDon, fGAINon, fXTon, fNonLinOn, fBypass, fCommonAdditive));

      if( fTBins != timeBins ) {
	AliError(Form("Mismatch between Number of Time Bins from CDB (%d) and from HC header (%d)"
                     , timeBins, fTBins));
      }

    }

  }

}  

//____________________________________________________________________________
void AliTRDRawStream::DecodeMCMheader()
{
  //
  //
  //

  if ( fRawVersion >= 1 && fRawVersion <= 3 ) {
    DecodeMCMheaderVx();
    return;
  }

  AliError(Form(" Unsupported raw version: %d", fRawVersion));
  return;

}

//____________________________________________________________________________
void AliTRDRawStream::DecodeMCMheaderVx()
{

  //
  // Decode the MCM header
  //

  fMCM  = (fDataWord & 0xff000000) >> 24;
  fEv   = (fDataWord & 0x00fffff0) >> 4;

  fROB  = fMCM / 16;
  fMCM  = fMCM % 16;

  fROW  = fGeo->GetPadRowFromMCM(fROB, fMCM);

  AliDebug(4, Form("0x%08x: SM%d L%dS%d. MCM Header: fROB=%d fMCM=%02d fEv=%02d"
		  , fDataWord, fSM, fLAYER, fSTACK, fROB, fMCM, fEv));

  if ( fROB % 2 == 0 && fSIDE == 1 ) {
    AliError(Form("SM%d L%dS%d: Mismatch between fROB (%d) and fSIDE (%d): fMCM=%02d"
                 , fSM, fLAYER, fSTACK, fROB, fSIDE, fMCM ));
  }
  if ( fROB % 2 != 0 && fSIDE == 0 ) {
    AliError(Form("SM%d L%dS%d: Mismatch between fROB (%d) and fSIDE (%d): fMCM=%02d"
                 , fSM, fLAYER, fSTACK, fROB, fSIDE, fMCM ));
  }
  if ( (fSTACK == 2 && fROW >= fGeo->RowmaxC0()) ||
       (fSTACK != 2 && fROW >= fGeo->RowmaxC1()) || fROW < 0 ) {
    AliError(Form("SM%d L%dS%d: Wrong Padrow (%d) fROB=%d, fSIDE=%d, fMCM=%02d"
                 , fSM, fLAYER, fSTACK, fROW, fROB, fSIDE, fMCM ));
  }
  
  fMCMHctr1++;
  fMCMHctr2++;

  // AdcMask for Zero supressed data
  if ( fRawVersion == 3 ) {
    // read one more word
    if ( !fRawReader->ReadNextInt( fDataWord ) ) {
      AliError("Could not read data");
      return;
    }
    fWordCtr++;
    if ( (fDataWord & 0x000007ff) == 0xC ) {     // at the moment bits 4-10 are empty
    //if ( (fDataWord & 0x0000000f) == 0xC ) {
      
      for ( Int_t ctr = 0; ctr < fGeo->ADCmax(); ctr++ ) {
	if ( (fDataWord >> (11+ctr)) == 0x1 ) fADCmask[ctr] = kTRUE;
	else                                  fADCmask[ctr] = kFALSE;
      }

      AliDebug(4, Form("0x%08x: ADC mask", fDataWord));

    }
    else {
      AliError("Expected ADC mask but did not find one!");
    }

  }

}

//____________________________________________________________________________
void AliTRDRawStream::DecodeTracklet()
{
  //
  //
  //

  if ( fRawVersion >= 1 && fRawVersion <= 3 ) {
    DecodeTrackletVx();
    return;
  }

  AliError(Form(" Unsupported raw version: %d", fRawVersion));
  return;

}

//____________________________________________________________________________
void AliTRDRawStream::DecodeTrackletVx()
{

  //
  // Decode the Tracklet
  //
  // this function is not tested yet on real tracklets
  //

  fTracklPID    = (fDataWord >> 24) & 0xff;
  fTracklPadRow = (fDataWord >> 20) & 0xf;    // 0:15
  fTracklDefL   = (fDataWord >> 13) & 0x7f;
  fTracklPadPos = (fDataWord)       & 0x1fff;

  fTracklPID    /= (Float_t)((1<<8) - 1);                      // 0:1 (steps of 0.39%)
  fTracklDefL    = (fTracklDefL  - ((1<< 7)-1)/2.) * 140.e-4;  // -0.889:0.889cm 
  fTracklPadPos  = (fTracklPadPos - ((1<<13)-1)/2.) * 160.e-4; // -65.528:65.528 cm

  AliDebug(4, Form("0x%08x: Tracklet found: SM%d L%dS%d side %x: PadRow=%d PadPos=%f DefL=%f PID=%f"
		  , fDataWord, fSM, fLAYER, fSTACK, fSIDE+10
                  , fTracklPadRow, fTracklPadPos, fTracklDefL, fTracklPID));

  if( (fSTACK == 2) && (fTracklPadRow >= (Int_t) fGeo->RowmaxC0()) ||
      (fSTACK != 2) && (fTracklPadRow >= (Int_t) fGeo->RowmaxC1()) ) {
    AliError(Form("Strange Row read from Tracklet Word: %d", fTracklPadRow));
  }

}

//____________________________________________________________________________
void AliTRDRawStream::DecodeGTUlinkMask()
{
  //
  //
  //

  if ( fRawVersion >= 1 && fRawVersion <= 3 ) {
    DecodeGTUlinkMaskVx();
    return;
  }

  AliError(Form(" Unsupported raw version: %d", fRawVersion));
  return;

}

//____________________________________________________________________________
void AliTRDRawStream::DecodeGTUlinkMaskVx()
{

  //
  // Decode the link masks sent by the GTU. These marke the active optical links
  // between GTU and Super Module. Up to now only fully active links are found
  // (0xfff = 12 active links).
  //

  if ( fGTUctr1 == -1 ) fGTUctr2++;
  fGTUctr1++;

  //printf("fGTUctr=%d",fGTUctr);
  //printf("fGTUctr1=%d, fGTUctr2=%d",fGTUctr1, fGTUctr2);

  if ( (fGTUctr1 >= 0) && (fGTUctr1 < 5) && (fGTUctr2 >= 0) && (fGTUctr2 < 18) ) {
    fGTUlinkMask[fGTUctr2][fGTUctr1] = (fDataWord & 0xfff);
  }

  for ( Int_t ctr = 0; ctr < 12; ctr++ ) {
    if ( IsGTULinkActive(fGTUctr2, ctr/2, fGTUctr1, ctr%2) ) {
      AliDebug(3, Form("SM %2d Stack %d: GTU Link %2d is active!", fGTUctr2, fGTUctr1, ctr)); 
    }
  }

}

//____________________________________________________________________________
Int_t  AliTRDRawStream::DecodeDataWord()
{

  //
  // Decode the Data
  //

  if      ( fRawVersion >= 1 && fRawVersion <= 2 ) {
    return DecodeDataWordV1V2();
  }
  else if ( fRawVersion >= 3 && fRawVersion <= 3 ) {
    return DecodeDataWordV3();
  }

  AliError(Form(" Unsupported raw version: %d", fRawVersion));
  return -1;

}

//____________________________________________________________________________
Int_t  AliTRDRawStream::DecodeDataWordV1V2()
{

  //
  // Decode the Data (full raw data. No zero suppression. 21 adc channels)
  //
  // return  0 means continue to next data word
  // return -1 means break data loop
  //

  if ( (fDataWord & 0x00000003) != 0x2 && (fDataWord & 0x00000003) != 0x3) {
    AliError(Form("Data %08x : Data Word ends neither with b11 nor b10", (Int_t)fDataWord));
    return -1;
  }

  if ( (fDataWord & 0x00000003) != fTbSwitch ) {    // Next ADC channel found
    //if ( fTB+1 != timeBins ) AliError(Form("Time bins in data (%d) != DB (%d)", fTB+1, timeBins));
    fTbSwitch = (fTbSwitch & 2) | !(fTbSwitch & 1);   // 0x3 <--> 0x2
    fTbSwitchCtr = 0;
    fADC++;
    fTB=0;
  }

  fTbSwitchCtr++; // Just read one word

  // We have only timeTotal time bins
  if ( fTbSwitchCtr > fTimeWords ) {
    AliError(Form("Data is strange. Already found %d words for this ADC channel", (Int_t)fTbSwitchCtr));
    return 0;
  }

  // We have only 21 ADC channels.
  if ( fADC > (Int_t)fGeo->ADCmax()-1 ) {
    AliError(Form("Data %08x : Data is strange. fADC is already %d", (Int_t)fDataWord, (Int_t)fADC));
    return 0;
  }

  // There are 18 pads connected to each MCM ADC channels 2...19. The other channels cross to other
  // MCMs and are good for online tracking in the MCM.
  if ( fADC > 1 && fADC < (Int_t)fGeo->ADCmax()-1 ) {

    // Get Pad column
    fCOL = fGeo->GetPadColFromADC(fROB, fMCM, fADC);

    // We have only 144 Pad Columns
    if ( fCOL > fColMax-1 || fCOL < 0 ) {
      AliError(Form("SM%d L%dS%d: Wrong Pad column (%d) fROB=%d, fSIDE=%d, fMCM=%02d", fSM,
		    fLAYER, fSTACK, fCOL, fROB, fSIDE, fMCM ));
    }

    // Decode 32 bit data words with information from 3 time bins and copy the data
    fSig[0] = (fDataWord & 0x00000ffc) >> 2;
    fSig[1] = (fDataWord & 0x003ff000) >> 12;
    fSig[2] = (fDataWord & 0xffc00000) >> 22;

    // Print data to screen:
    AliDebug(5, Form("SM%d L%dS%d: ROB%d MCM=%d ADC=%d (ROW=%d COL=%d): Data %04d %04d %04d\n",
		     fSM, fLAYER, fSTACK, fROB, fMCM, fADC, fROW, fCOL, fSig[0], fSig[1], fSig[2]));
    
  }
  else {
    
    fCOL = -1;
    
  }

  return 1;

}

//____________________________________________________________________________
Int_t  AliTRDRawStream::DecodeDataWordV3()
{

  //
  // Decode the data (Zero suppresses data. 21 adc channels)
  //
  // return  0 means continue to next data word
  // return -1 means break data loop
  //
  // NOT TESTED YET!!!!!!!!
  //

  if ( (fDataWord & 0x00000003) != 0x2 && (fDataWord & 0x00000003) != 0x3) {
    AliError(Form("Data %08x : Data Word ends neither with b11 nor b10", (Int_t)fDataWord));
    return -1;
  }

  if ( (fDataWord & 0x00000003) != fTbSwitch ) {    // Next ADC channel found
    //if ( fTB+1 != timeBins ) AliError(Form("Time bins in data (%d) != DB (%d)", fTB+1, timeBins));
    fTbSwitch = (fTbSwitch & 2) | !(fTbSwitch & 1);   // 0x3 <--> 0x2
    fTbSwitchCtr = 0;
    //
    // Jump to next ADC channel that is not masked
    do {
      fADC++;
    } while ( ((fADC < fGeo->ADCmax()) && (fADCmask[fADC] == kFALSE)) || (fADC >= fGeo->ADCmax()) );
    fTB=0;
  }

  fTbSwitchCtr++; // Just read one word

  // We have only timeTotal time bins
  if ( fTbSwitchCtr > fTimeWords ) {
    AliError(Form("Data is strange. Already found %d words for this ADC channel", (Int_t)fTbSwitchCtr));
    return 0;
  }

  // We have only 21 ADC channels.
  if ( fADC > (Int_t)fGeo->ADCmax()-1 ) {
    AliError(Form("Data %08x : Data is strange. fADC is already %d", (Int_t)fDataWord, (Int_t)fADC));
    return 0;
  }

  // There are 18 pads connected to each MCM ADC channels 2...19. The other channels cross to other
  // MCMs and are good for online tracking in the MCM.
  if ( fADC > 1 && fADC < (Int_t)fGeo->ADCmax()-1 ) {

    // Get Pad column
    fCOL = fGeo->GetPadColFromADC(fROB, fMCM, fADC);

    // We have only 144 Pad Columns
    if ( fCOL > fColMax-1 || fCOL < 0 ) {
      AliError(Form("SM%d L%dS%d: Wrong Pad column (%d) fROB=%d, fSIDE=%d, fMCM=%02d", fSM,
		    fLAYER, fSTACK, fCOL, fROB, fSIDE, fMCM ));
    }

    // Decode 32 bit data words with information from 3 time bins and copy the data
    fSig[0] = (fDataWord & 0x00000ffc) >> 2;
    fSig[1] = (fDataWord & 0x003ff000) >> 12;
    fSig[2] = (fDataWord & 0xffc00000) >> 22;

    // Print data to screen:
    AliDebug(5, Form("SM%d L%dS%d: ROB%d MCM=%d ADC=%d (ROW=%d COL=%d): Data %04d %04d %04d\n",
		     fSM, fLAYER, fSTACK, fROB, fMCM, fADC, fROW, fCOL, fSig[0], fSig[1], fSig[2]));
    
  }
  else {
    
    fCOL = -1;
    
  }

  return 1;

}
