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
  ,fFiltered(0)
  ,fHCHctr1(0)
  ,fHCHctr2(0)
  ,fMCMHctr1(0)
  ,fMCMHctr2(0)
  ,fGTUctr1(0)
  ,fGTUctr2(0)
  ,fTracklPID(0)
  ,fTracklDefL(0)
  ,fTracklPadPos(0)
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
  ,fRawVersion(1)
  ,fDataWord(0)
  ,fStatus(0)
  ,fRowMax(0)
  ,fColMax(0)
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
  ,fFiltered(0)
  ,fHCHctr1(0)
  ,fHCHctr2(0)
  ,fMCMHctr1(0)
  ,fMCMHctr2(0)
  ,fGTUctr1(0)
  ,fGTUctr2(0)
  ,fTracklPID(0)
  ,fTracklDefL(0)
  ,fTracklPadPos(0)
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
  ,fRawVersion(1)
  ,fDataWord(0)
  ,fStatus(0)
  ,fRowMax(0)
  ,fColMax(0)
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
  ,fFiltered(0)
  ,fHCHctr1(0)
  ,fHCHctr2(0)
  ,fMCMHctr1(0)
  ,fMCMHctr2(0)
  ,fGTUctr1(0)
  ,fGTUctr2(0)
  ,fTracklPID(0)
  ,fTracklDefL(0)
  ,fTracklPadPos(0)
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
  ,fRawVersion(1)
  ,fDataWord(0)
  ,fStatus(0)
  ,fRowMax(0)
  ,fColMax(0)
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
  ,fFiltered(0)
  ,fHCHctr1(0)
  ,fHCHctr2(0)
  ,fMCMHctr1(0)
  ,fMCMHctr2(0)
  ,fGTUctr1(0)
  ,fGTUctr2(0)
  ,fTracklPID(0)
  ,fTracklDefL(0)
  ,fTracklPadPos(0)
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
  ,fRawVersion(1)
  ,fDataWord(0)
  ,fStatus(0)
  ,fRowMax(0)
  ,fColMax(0)
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

  delete fGeo;
  delete fRawReader;
  //delete fDigitsManager;
  delete fDigits;
  delete fTrack0;
  delete fTrack1;
  delete fTrack2;

}

//_____________________________________________________________________________
Bool_t AliTRDRawStream::SetRawVersion(Int_t rv)
{
  //
  // Set the raw data version
  //

  if ( rv >= 0 && rv <= 2 ) {
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
Bool_t AliTRDRawStream::ReadAll()
{

  //
  // Read all TRD raw data word (32 bits). This is for all FrawVersion > 0.
  // Return kFALSE if something is not cool
  //
  // by C. Lippmann
  //
 
  AliTRDCommonParam *commonParam = AliTRDCommonParam::Instance();
  if (!commonParam) {
    AliError("Could not get common parameters");
    return kFALSE;
  }

  AliTRDcalibDB     *calibration = AliTRDcalibDB::Instance();
  if (!calibration) {
    AliError("Could not get calibration object");
    return kFALSE;
  }
  
  UInt_t timeTotal = calibration->GetNumberOfTimeBins();

  // The number of data words needed for this number of time bins (there
  // are 3 time bins in one word)
  UInt_t timeWords = ( timeTotal%3 == 0 ) ? timeTotal/3 :  timeTotal/3 + 1;

  AliDebug(2, Form("Number of Timebins read from CDB: %d", timeTotal));

  UInt_t TBswitch    = 3;
  UInt_t TBswitchCtr = 0;
  Int_t  WordCtr     = 0;
  Int_t  EqID        = 0;
  Int_t  datasize    = 0;
  Int_t  iDET        = 0;

  fHCHctr1 = fHCHctr2 =  0;
  fGTUctr1 = fGTUctr2 = -1;

  AliInfo("Converting TRD raw data to digits ...");

  while ( 1 ) { // loop over all supermodules

    WordCtr   = 0;
    fHCHctr1  = 0;
    fMCMHctr1 = 0;

    //
    // 0) Find first GTU Link Mask and test if we can read data
    //
    do {

      if ( !fRawReader->ReadNextInt( fDataWord ) ) {
	AliInfo(Form("Finished processing TRD raw data: Found %d Half-Chambers", fHCHctr2));
	return kTRUE;
      }
      WordCtr++;

      // After reading the first word check for size of this data and get Eq. ID
      if ( WordCtr == 1 ) {
	datasize = fRawReader->GetDataSize()/4;  // Size of this payload is in 32bit words
	EqID     = fRawReader->GetEquipmentId(); // Get Equipment ID
      }

      // GTU Link Mask?
      if ( (fDataWord & 0xfffff000) ==  0xe0000000 ) {
	fStatus = 1;      // GTU link mask found
	DecodeGTUlinkMask();
	break;
      } 
      else {
	AliError(Form("Equipment %d: First data word is not GTU Link Mask!", EqID));
	return kFALSE;
      }

    } 
    while ( WordCtr < datasize );

    //
    // loop over all half chambers in one supermodule
    //
    while ( WordCtr < datasize ) {

      //
      // 1) Find end_of_tracklet_marker
      //
      while ( WordCtr < datasize ) {

	if ( !fRawReader->ReadNextInt( fDataWord ) ) {
	  AliError("Could not read data");
	  return kFALSE;
	}
	WordCtr++;

	// GTU Link Mask?
	if ( (fDataWord & 0xfffff000) ==  0xe0000000 ) {
	  DecodeGTUlinkMask();
	  continue;
	}

	// end_of_tracklet_marker?
	if ( fDataWord == end_of_tracklet_marker ) {
	  AliDebug(3, "end_of_tracklet_marker found");
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
      while ( WordCtr < datasize ) { 

	if ( !fRawReader->ReadNextInt( fDataWord ) ) {
	  AliError("Could not read data");
	  return kFALSE;
	}
	WordCtr++;
	//printf("Word %d: 0x%08x\n", WordCtr, fDataWord); 

	if ( fDataWord != end_of_tracklet_marker ) {
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
	iDET = fGeo->GetDetector(fLAYER, fSTACK, fSM);
	if ( fChamberDone[iDET] == 2 ) {
	  fDigits->Compress(1,0);
	  fTrack0->Compress(1,0);
	  fTrack1->Compress(1,0);
	  fTrack2->Compress(1,0);
	}

	// Read from new HC header the chamber position (fLAYER, fSTACK, fSM)
	DecodeHCheader(timeTotal);
	WordCtr += fHCHWords;
	iDET    = fGeo->GetDetector(fLAYER, fSTACK, fSM);
	fRowMax = commonParam->GetRowMax(fLAYER,fSTACK,fSM);
	fColMax = commonParam->GetColMax(fROC);

	// Add a container for the digits of this detector
	fDigits = fDigitsManager->GetDigits(iDET);
	fTrack0 = fDigitsManager->GetDictionary(iDET,0);
	fTrack1 = fDigitsManager->GetDictionary(iDET,1);
	fTrack2 = fDigitsManager->GetDictionary(iDET,2);
	
	fChamberDone[iDET]++;
	
	// Allocate memory if it was not already done
	if (fDigits->GetNtime() == 0) {
	  fDigits->Allocate(fRowMax,fColMax,timeTotal);
	  fTrack0->Allocate(fRowMax,fColMax,timeTotal);
	  fTrack1->Allocate(fRowMax,fColMax,timeTotal);
	  fTrack2->Allocate(fRowMax,fColMax,timeTotal);
	}

	fMCMHctr2 = 0;

      }
    
      //
      // 4) Scan MCM data
      //
      fStatus = 0;
      while ( WordCtr < datasize ) {

	if ( !fRawReader->ReadNextInt( fDataWord ) ) {
	  AliError("Could not read data");
	  return kFALSE;
	}
	WordCtr++;
	//printf("Word %d: 0x%08x\n", WordCtr, fDataWord); 
      
	//if ( WordCtr == 4*datasize ) AliInfo(Form("Achtung! WordCtr=%d (%d)", WordCtr, 4*datasize));

	if( (fDataWord & 0x0000000f) == 0xC ) { // MCM Header
	  DecodeMCMheader();
	  if ( fMCM < 0 || fMCM > 15 || fROB < 0 || fROB > 7 ) {
	    AliError("Wrong fMCM or fROB. Skip this data");
	    break;
	  }
	  TBswitch    = 3;  // For first adc channel we expect: (fDataWord & 3) = 3
	  TBswitchCtr = 0;  // 
	  fADC = fTB  = 0;  // Reset Counter
	  fStatus     = 1;  // Now 1 means MCM header is found
	  continue;
	}
    
	// End of half-chamber data, finished:
	if ( fDataWord == end_of_event_marker ) {
	  fGTUctr1 = -1;
	  break;
	}

	if ( fStatus == 1 ) {       // MCM header is set, ADC data is valid.
    
	  //
	  // Found some data. Decode it now:
	  //
	  if ( (fDataWord & 0x00000003) != 0x2 && (fDataWord & 0x00000003) != 0x3) {
	    AliError(Form("Data %08x : Data Word ends neither with 11 nor 10", (Int_t)fDataWord));
	    break;
	  }

	  if ( (fDataWord & 0x00000003) != TBswitch ) {    // Next ADC channel found
	    //if ( fTB+1 != timeBins ) AliError(Form("Time bins in data (%d) != DB (%d)", fTB+1, timeBins));
	    TBswitch = (TBswitch & 2) | !(TBswitch & 1);   // 0x3 <--> 0x2
	    TBswitchCtr = 0;
	    fADC++;
	    fTB=0;
	  }

      	  TBswitchCtr++; // Just read one word
	
	  // We have only timeTotal time bins
	  if ( TBswitchCtr > timeWords ) {
	    AliError(Form("Data is strange. Already found %d words for this ADC channel", (Int_t)TBswitchCtr));
	    continue;
	  }

	  // We have only 21 ADC channels.
	  if ( fADC > 20 ) {
	    AliError(Form("Data %08x : Data is strange. fADC is already %d", (Int_t)fDataWord,
			  (Int_t)fADC));
	    continue;
	  }

	  // There are 18 pads connected to each MCM ADC channels 2...19. The
	  // other channels cross to other MCMs and are good for online tracking
	  // in the MCM.
	  if ( fADC > 1 && fADC < (Int_t)fGeo->ADCmax()-1 ) {

	    fCOL = fGeo->GetPadCol(fROB, fMCM, fADC);

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

	    // Write Digits
	    if ( fCOL >= 0 && fCOL < fColMax && fROW >= 0 && fROW < fRowMax ) {  // A real pad
	      for ( Int_t ctr = 0; ctr <3; ctr++ ) {
		if ( fTB+ctr < (Int_t)timeTotal ) {
		  fDigits->SetDataUnchecked(fROW, fCOL, fTB+ctr, fSig[ctr]);
		  fTrack0->SetDataUnchecked(fROW, fCOL, fTB+ctr, 0);
		  fTrack1->SetDataUnchecked(fROW, fCOL, fTB+ctr, 0);
		  fTrack2->SetDataUnchecked(fROW, fCOL, fTB+ctr, 0);
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

      AliDebug(2, Form("SM%d L%dS%d side %x: Processed %d MCMs.", fSM, fLAYER, fSTACK, fSIDE+10, fMCMHctr2));

    } // End Half-Chamber loop

    AliDebug(1, Form("SM%d (Eq %d): Processed %d HC (%d MCMs)", fSM, EqID, fHCHctr1, fMCMHctr1));

  } // End Super Module loop

  return kTRUE;

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

  if ( (fDataWord >> 31) != 1 )  {

    if ( fRawVersion == 1 ) {
      DecodeHCheaderV1();
    }
    else {
      AliError(Form("Mismatch between fRawVersion (%d) and HC header signature", fRawVersion));
    }
    return;

  } 
  else {

    fRVmajor = (fDataWord >> 24) & 0x7f;
    fRVminor = (fDataWord >> 17) & 0x7f;
    if ( fRawVersion != fRVmajor ) {
      AliError(Form("Mismatch between fRawVersion (%d) and fRVmajor from HC header(%d)"
                   ,fRawVersion,fRVmajor));
    }
    if (fRawVersion == 2 ) {
      DecodeHCheaderV2(timeBins);
    }
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
      AliError(Form("0x%08x: Strange HC header: dcs=%d; sm=%d; layer=%d; stack=%d.",
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
void AliTRDRawStream::DecodeHCheaderV2(Int_t timeBins)
{
  //
  // Decode the HC header (fRawVersion == 2, Full raw production)
  //

  // 1st word
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

  // 2nd word
  if ( fHCHWords >= 1 ) {
    // read one more word
    if ( !fRawReader->ReadNextInt( fDataWord ) ) {
      AliError("Could not read data");
      return;
    }
    if ( (fDataWord & 0x3) == 1 ) {
      
      fBCctr   =  (fDataWord >> 16);
      fPTctr   =  (fDataWord >> 12) & 0xf;
      fPTphase =  (fDataWord >>  8) & 0xf;
      fTBins   = ((fDataWord >>  2) & 0x3f) + 1;

      AliDebug(3, Form("0x%08x: HC header 2: BCctr=%d PTctr=%d PTph=%d TB=%d"
                      , fDataWord, fBCctr, fPTctr, fPTphase, fTBins));

      if( fTBins != timeBins ) {
	AliError(Form("Mismatch between Number of Time Bins from CDB (%d) and from HC header (%d)"
                     , timeBins, fTBins));
      }

    }

  }

  // 3rd word
  if ( fHCHWords >= 2 ) {
    // read one more word
    if ( !fRawReader->ReadNextInt( fDataWord ) ) {
      AliError("Could not read data");
      return;
    }
    if ( (fDataWord & 0x3) == 1 ) {
      /*
      Not finished. Next to come:
      fTCon 
      fPEDon
      fGAINon
      fFiltered
      .....
    */
    }

  }

}  

//____________________________________________________________________________
void AliTRDRawStream::DecodeMCMheader()
{
  //
  //
  //

  if ( fRawVersion >= 1 && fRawVersion <= 2 ) {
    DecodeMCMheaderV1();
    return;
  }

  AliError(Form(" Unsupported raw version: %d", fRawVersion));
  return;

}

//____________________________________________________________________________
void AliTRDRawStream::DecodeMCMheaderV1()
{

  //
  // Decode the MCM header
  //

  fMCM  = (fDataWord & 0xff000000) >> 24;
  fEv   = (fDataWord & 0x00fffff0) >> 4;

  fROB  = fMCM / 16;
  fMCM  = fMCM % 16;

  fROW  = fGeo->GetPadRow(fROB, fMCM);

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
  if ( (fSTACK == 2 && fROW > 11) || (fSTACK != 2 && fROW > 15) || fROW < 0 ) {
    AliError(Form("SM%d L%dS%d: Wrong Padrow (%d) fROB=%d, fSIDE=%d, fMCM=%02d"
                 , fSM, fLAYER, fSTACK, fROW, fROB, fSIDE, fMCM ));
  }
  
  fMCMHctr1++;
  fMCMHctr2++;

}

//____________________________________________________________________________
void AliTRDRawStream::DecodeTracklet()
{
  //
  //
  //

  if ( fRawVersion >= 1 && fRawVersion <= 2 ) {
    DecodeTrackletV1();
    return;
  }

  AliError(Form(" Unsupported raw version: %d", fRawVersion));
  return;

}

//____________________________________________________________________________
void AliTRDRawStream::DecodeTrackletV1()
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

  if( (fSTACK == 2) && (fTracklPadRow >= (Int_t)fGeo->RowmaxC0) ||
      (fSTACK != 2) && (fTracklPadRow >= (Int_t)fGeo->RowmaxC1) ) {
    AliError(Form("Strange Row read from Tracklet Word: %d", fTracklPadRow));
  }

}

//____________________________________________________________________________
void AliTRDRawStream::DecodeGTUlinkMask()
{
  //
  //
  //

  if ( fRawVersion >= 1 && fRawVersion <= 2 ) {
    DecodeGTUlinkMaskV1();
    return;
  }

  AliError(Form(" Unsupported raw version: %d", fRawVersion));
  return;

}

//____________________________________________________________________________
void AliTRDRawStream::DecodeGTUlinkMaskV1()
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

