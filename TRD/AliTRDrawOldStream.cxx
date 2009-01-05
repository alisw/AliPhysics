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

/* $Id: AliTRDrawOldStream.cxx 23562 2008-01-25 15:35:04Z cblume $ */

////////////////////////////////////////////////////////////////////////////
//                                                                        //
// This class provides access to TRD digits in raw data.                  //
//                                                                        //
// It loops over all TRD digits in the raw data given by the AliRawReader //
// The Next method goes to the next digit. If there are no digits left    //
// it returns kFALSE.                                                     //
// Several getters provide information about the current digit.           //
//                                                                        //
// Author:                                                                //
//   Christian Lippmann (C.Lippmann@gsi.de)                               //
//                                                                        //
////////////////////////////////////////////////////////////////////////////

#include "AliLog.h"
#include "AliRawReader.h"

#include "AliTRDrawOldStream.h"
#include "AliTRDgeometry.h"
#include "AliTRDcalibDB.h"
#include "AliTRDdigitsManager.h"
#include "AliTRDarrayDictionary.h"
#include "AliTRDarrayADC.h"
#include "AliTRDSignalIndex.h"
#include "AliTRDfeeParam.h"

ClassImp(AliTRDrawOldStream)

//_____________________________________________________________________________
AliTRDrawOldStream::AliTRDrawOldStream() 
  :AliTRDrawStreamBase()
//  :TObject()
  ,fGeo(NULL) 
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
  ,fDET(0)
  ,fLastDET(-1)
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
  ,fRawVersion(2)
  ,fNextStatus(0)
  ,fTbSwitch(0)
  ,fTbSwitchCtr(0)
  ,fTimeWords(0)
  ,fWordCtr(0)
  ,fRowMax(0)
  ,fColMax(0)
  ,fADCmask()
  ,fChamberDone()
  ,fRetVal(0)
  ,fEqID(0)
  ,fDataSize(0)
  ,fSizeOK(kFALSE)
  ,fCountBytes(0)
  ,fBufSize(0)
  ,fBufferSet(kFALSE)
  ,fPos(NULL)
  ,fDataWord(NULL)
  ,fTimeBinsCalib(0)
{
  //
  // Default constructor
  //

  for (Int_t i = 0; i < 540; i++) {
    fChamberDone[i] = 0;
  }

}

//_____________________________________________________________________________
AliTRDrawOldStream::AliTRDrawOldStream(AliRawReader *rawReader) 
  :AliTRDrawStreamBase(rawReader)
//  :TObject()
  ,fGeo(NULL) 
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
  ,fDET(0)
  ,fLastDET(-1)
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
  ,fRawVersion(2)
  ,fNextStatus(0)
  ,fTbSwitch(0)
  ,fTbSwitchCtr(0)
  ,fTimeWords(0)
  ,fWordCtr(0)
  ,fRowMax(0)
  ,fColMax(0)
  ,fADCmask()
  ,fChamberDone()
  ,fRetVal(0)
  ,fEqID(0)
  ,fDataSize(0)
  ,fSizeOK(kFALSE)
  ,fCountBytes(0)
  ,fBufSize(0)
  ,fBufferSet(kFALSE)
  ,fPos(NULL)
  ,fDataWord(NULL)
  ,fTimeBinsCalib(0)
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
AliTRDrawOldStream::AliTRDrawOldStream(const AliTRDrawOldStream& stream)
  :AliTRDrawStreamBase(stream)
//  :TObject()
  ,fGeo(NULL)
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
  ,fDET(0)
  ,fLastDET(-1)
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
  ,fRawVersion(-1)
  ,fNextStatus(0)
  ,fTbSwitch(0)
  ,fTbSwitchCtr(0)
  ,fTimeWords(0)
  ,fWordCtr(0)
  ,fRowMax(-1)
  ,fColMax(-1)
  ,fADCmask()
  ,fChamberDone()
  ,fRetVal(0)
  ,fEqID(0)
  ,fDataSize(0)
  ,fSizeOK(kFALSE)
  ,fCountBytes(0)
  ,fBufSize(0)
  ,fBufferSet(kFALSE)
  ,fPos(NULL)
  ,fDataWord(NULL)
  ,fTimeBinsCalib(0)
{
  //
  // Copy constructor
  //

  AliFatal("Copy constructor not implemented");

}

//_____________________________________________________________________________
AliTRDrawOldStream& AliTRDrawOldStream::operator = (const AliTRDrawOldStream& 
					      /* stream */)
{
  //
  // Assigment operator
  //

  AliFatal("Assignment operator not implemented");
  return *this;

}

//_____________________________________________________________________________
AliTRDrawOldStream::~AliTRDrawOldStream()
{
  //
  // Destructor
  //

  if (fGeo) {  
    delete fGeo;
  }

}

//_____________________________________________________________________________
void AliTRDrawOldStream::SetRawReader(AliRawReader *rawReader) 
{
  //
  // Sets the raw reader
  //

  if (rawReader)
    {
      fRawReader = rawReader;
    }

}

//_____________________________________________________________________________
Bool_t AliTRDrawOldStream::SetRawVersion(Int_t rv)
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

//____________________________________________________________________________
Bool_t AliTRDrawOldStream::Init()
{
  //
  // Initialization
  //

  if (!AliTRDcalibDB::Instance()) {
    AliError("Could not get calibration object");
    return kFALSE;
  }

  if (!fGeo) {
    fGeo = new AliTRDgeometry();
  }
  
  fTimeBinsCalib = AliTRDcalibDB::Instance()->GetNumberOfTimeBins();
  AliDebug(2, Form("Number of Timebins read from CDB: %d", fTimeBinsCalib));

  // The number of data words needed for this number of time bins (there
  // are 3 time bins in one word)
  fTimeWords = (fTimeBinsCalib-1)/3 + 1;

  fTbSwitch    = 3;
  fTbSwitchCtr = 0;

  fHCHctr1 = fHCHctr2 =  0;
  fGTUctr1 = fGTUctr2 = -1;

  fHCdataCtr = 0;
  fWordCtr   = 0;  

  fDET     = 0;
  fLastDET     = -1;
  fRetVal = 0;
  fEqID     = 0;
  fDataSize = 0;
  fSizeOK = kFALSE;
  
  fNextStatus = kStart;

  fCountBytes = 0;
  fBufSize = 0;
  fDataWord = NULL;
  fPos = NULL;
  fWordCtr = 0;
  fBufferSet = kFALSE;

  return kTRUE;

}

//____________________________________________________________________________
Int_t AliTRDrawOldStream::NextData()
{
  //
  // Updates the next data word pointer
  //

  if (fCountBytes + kSizeWord >= fBufSize)
    {
      fBufferSet = fRawReader->ReadNextData(fPos);
      if (fBufferSet == kTRUE)
	{
	  fBufSize = fRawReader->GetDataSize();
	  fCountBytes = 0;	  
	  fDataWord = (UInt_t*)fPos;
	  fNextStatus = kNextSM;
	  fWordCtr = 0;
	  return kNextSM;
	}
      else
	{
	  fNextStatus = kStop;
	  return kNoMoreData;
	}
    }
  else
    {

      fPos += kSizeWord;
      fCountBytes += kSizeWord;	  
      fDataWord = (UInt_t*)fPos;
      fWordCtr++;
      return kWordOK;
    }
}

//____________________________________________________________________________
Bool_t AliTRDrawOldStream::Next()
{
  //
  // Updates the next data word pointer
  //

  if (fNextStatus == kStart)
    {
      Init();
    }

  while (fNextStatus != kStop)
    { // !kStop
      NextData();
      if (fNextStatus == kNextMCM || fNextStatus == kNextData)
	{
	  fHCdataCtr += 4;
	  
	  if( ((*fDataWord & 0x80000000) == 0x0) && ((*fDataWord & 0x0000000f) == 0xC) )
	    { // MCM Header
	      DecodeMCMheader();
	      if ( fMCM < 0 || fMCM > 15 || fROB < 0 || fROB > 7 ) 
		{
		  AliWarning("Wrong fMCM or fROB. Skip this data");
		  fRawReader->AddMajorErrorLog(kWrongMCMorROB,Form("MCM=%d, ROB=%d",fMCM,fROB));
		  fNextStatus = kNextHC;
		  continue;
		}
	      fTbSwitch    = 3;  // For first adc channel we expect: (*fDataWord & 3) = 3
	      fTbSwitchCtr = 0;  // 
	      fADC = fTB   = 0;  // Reset Counter
	      fNextStatus = kNextData;
	      continue;
	    }
	  
	  if ( *fDataWord == kEndofrawdatamarker ) 
	    {  // End of half-chamber data, finished
	      fGTUctr1 = -1;
	      fNextStatus = kNextHC;
	      continue;
	    }
	  
	  if (fNextStatus == kNextData )
	    {       // MCM header is set, ADC data is valid.
	      
	      // Found some data. Decode it now:
	      fRetVal = DecodeDataWord();
	      if ( fRetVal ==  0 ) continue;
	      if ( fRetVal == -1 ) 
		{
		  fNextStatus = kNextHC;
		  continue;
		}
	      if ( fRetVal == 1)
		{
		  {  // A real pad
		    fTB += 3;		
		    return kTRUE;
		  }	       		
		}
	      // following ifs have been moved to DEcodeDatawordV1V2
	      // 	    if ( fADC > 1 && fADC < (Int_t)fGeo->ADCmax()-1 ) 
	      // 	      {	      
	      // 		// Write Digits
	      // 		if ( fCOL >= 0 && fCOL < fColMax && fROW >= 0 && fROW < fRowMax ) 
	      // 		  {  // A real pad
	      // 		    fTB += 3;		
	      // 		    return kTRUE;
	      // 		  }	       
	      // 	      }
	      // 	    else 
	      // 	      {
	      // 		fCOL = -1;	       
	      // 	      }
	    }// kNextData  
	  
	  continue;
	} //next mcm

      if ( fNextStatus == kNextHC )
	{
	  //
	  // 1) Find end_of_tracklet_marker
	  //
	  // GTU Link Mask?
	  if ( (*fDataWord & 0xfffff000) ==  0xe0000000 ) 
	    {
	      DecodeGTUlinkMask();
	      continue;
	    }
	  
	  // endoftrackletmarker?
	  if ( *fDataWord == kEndoftrackletmarker ) 
	    {
	      AliDebug(3, "end-of-tracklet-marker found");
	      fNextStatus = kSeekNonEoTracklet;
	      continue;
	    } 
	  else 
	    {
	      // Tracklets found
	      AliDebug(3, "Tracklet found");
	      DecodeTracklet();
	      continue;
	    }
	} //if next HC

      if (fNextStatus == kSeekNonEoTracklet)
	{
	  //
	  // 2) Look for non-end_of_tracklet_marker
	  //
	  //printf("Word %d: 0x%08x\n", fWordCtr, *fDataWord); 

	  if ( *fDataWord != kEndoftrackletmarker ) 
	    {
	      fNextStatus = kDecodeHC;
	      AliDebug(3, "NON end-of-tracklet-marker found");
	      //// no do not continue - this should be the hcheader
	    }
	  else
	    {
	      //just go on and find the non-end_of_tracklet_marker
	      continue;
	    }
	}

      if ( fNextStatus == kDecodeHC )
	{
	  AliDebug(3, "Decode HC");

	  //
	  // 3) This Word must be Half Chamber Header
	  //
	  if ( (*fDataWord & 0x00000003) == 1 ) 
	    { // HC header
	      DecodeHCheader(fTimeBinsCalib); // This is the new header!
	      fLastDET = fDET;
	      fDET    = fGeo->GetDetector(fLAYER, fSTACK, fSM);
	      fRowMax = fGeo->GetRowMax(fLAYER,fSTACK,fSM);
	      fColMax = fGeo->GetColMax(fROC);
	      
	      fMCMHctr2 = 0;
	      fHCdataCtr = 0;
	      
	      fChamberDone[fDET]++;
	      fNextStatus = kNextMCM;
	      AliDebug(3, "Decode HC OK");	      
	      continue;
	    } //HC header
	  else
	    {
	      AliDebug(3, "Decode HC NOT OK");	      
	      fNextStatus = kNextSM;
	      continue;
	    }
	} // if decode HC

      if (fNextStatus == kNextSM)
	{
	  
	  fDET     = 0;
	  fRetVal = 0;
	  fEqID     = 0;
	  fDataSize = 0;
	  fSizeOK = kFALSE;
	  
	  // After reading the first word check for size of this data and get Eq. ID
	  if ( fWordCtr == 1 ) 
	    {
	      fDataSize = fRawReader->GetDataSize()/4;  // Size of this payload in 32bit words
	      fEqID     = fRawReader->GetEquipmentId(); // Get Equipment ID
	      if ( fDataSize > 0 ) fSizeOK = kTRUE;
	    }
	  
	  // GTU Link Mask?
	  if ( (*fDataWord & 0xfffff000) ==  0xe0000000 ) 
	    {
	      DecodeGTUlinkMask();
	      fNextStatus = kNextHC;
	      continue;
	    } 
	  else 
	    {
	      AliWarning(Form("Equipment %d: First data word is not GTU Link Mask!", fEqID));
              fRawReader->AddMajorErrorLog(kGTULinkMaskMissing,Form("Equipment %d",fEqID));
	      fNextStatus = kStop;
	    }	    
	}// if nextSM

    } // not kStop

  AliDebug(1, Form("That's all folks! %d", fSM));
  return kFALSE;
}

//____________________________________________________________________________
Int_t AliTRDrawOldStream::NextChamber(AliTRDdigitsManager *man, UInt_t** /*trackletContainer*/)
{
  //
  // Updates the next data word pointer
  //

  AliTRDarrayADC *digits = 0;
  AliTRDarrayDictionary *track0 = 0;
  AliTRDarrayDictionary *track1 = 0;
  AliTRDarrayDictionary *track2 = 0; 
  AliTRDSignalIndex *indexes = 0;
	  
  if (fNextStatus == kStart)
    {
      Init();
    }

  while (fNextStatus != kStop)
    { // !kStop
      NextData();
      if (fNextStatus == kNextMCM || fNextStatus == kNextData)
      //while (fNextStatus == kNextMCM || fNextStatus == kNextData)
	{
	  fHCdataCtr += 4;
	  
	  if( ((*fDataWord & 0x80000000) == 0x0) && ((*fDataWord & 0x0000000f) == 0xC) )
	    { // MCM Header
	      DecodeMCMheader();
	      if ( fMCM < 0 || fMCM > 15 || fROB < 0 || fROB > 7 ) 
		{
		  AliWarning("Wrong fMCM or fROB. Skip this data");
		  fRawReader->AddMajorErrorLog(kWrongMCMorROB,Form("MCM=%d, ROB=%d",fMCM,fROB));
		  fNextStatus = kNextHC;
		  continue;
		}
	      fTbSwitch    = 3;  // For first adc channel we expect: (*fDataWord & 3) = 3
	      fTbSwitchCtr = 0;  // 
	      fADC = fTB   = 0;  // Reset Counter
	      fNextStatus  = kNextData;

// 	      NextData(); // if while loop!
	      continue; // if if
	    }
	  
	  if ( *fDataWord == kEndofrawdatamarker ) 
	    {  // End of half-chamber data, finished
	      fGTUctr1 = -1;
	      fNextStatus = kNextHC;
	      // full chamber processed ?
	      if (fChamberDone[fDET] == 2)
		{
		  return fDET;
		}
	      else
		{
// 		  break; // if while loop
		  continue; // if if
		}
	    }
	  
	  if (fNextStatus == kNextData )
	    {       // MCM header is set, ADC data is valid.
	      
	      // Found some data. Decode it now:
	      fRetVal = DecodeDataWord();
	      if ( fRetVal ==  0 ) continue;
	      if ( fRetVal == -1 ) 
		{
		  fNextStatus = kNextHC;

// 		  NextData(); // if while loop!
// 		  break; //if while loop!
		  continue;// if if
		}
	      
	    if ( fRetVal == 1)
	      {
		{  // A real pad
		  // here fill the data arrays
		  //return kTRUE;
		  for (Int_t it = 0; it < 3; it++)
		    {
		      if ( GetTimeBin() + it < GetNumberOfTimeBins() )
			{
			  if (GetSignals()[it] > 0)
			    {
			      digits->SetData(fROW, fCOL, fTB + it, fSig[it]);
			      indexes->AddIndexRC(fROW, fCOL);
			      if (man->UsesDictionaries())
				{
				  track0->SetData(fROW, fCOL, fTB + it, 0);
				  track1->SetData(fROW, fCOL, fTB + it, 0);
				  track2->SetData(fROW, fCOL, fTB + it, 0);
				}
			    }
			} // check the tbins range
		    } // for each tbin of current 3
		  fTB += 3;		
		}// real pad	       		
	      } // if fRetVal == 1
	    
	    // following ifs have been moved to DEcodeDatawordV1V2
// 	    if ( fADC > 1 && fADC < (Int_t)fGeo->ADCmax()-1 ) 
// 	      {	      
// 		// Write Digits
// 		if ( fCOL >= 0 && fCOL < fColMax && fROW >= 0 && fROW < fRowMax ) 
// 		  {  // A real pad
// 		    fTB += 3;		
// 		    return kTRUE;
// 		  }	       
// 	      }
// 	    else 
// 	      {
// 		fCOL = -1;	       
// 	      }
	    }// kNextData  
	  
// 	  NextData(); // if while loop!
	  continue; //if if
	} //next mcm

      if ( fNextStatus == kNextHC )
	{
	  //
	  // 1) Find end_of_tracklet_marker
	  //
	  // GTU Link Mask?
	  if ( (*fDataWord & 0xfffff000) ==  0xe0000000 ) 
	    {
	      DecodeGTUlinkMask();
	      continue;
	    }
	  
	  // endoftrackletmarker?
	  if ( *fDataWord == kEndoftrackletmarker ) 
	    {
	      AliDebug(3, "end-of-tracklet-marker found");
	      fNextStatus = kSeekNonEoTracklet;
	      continue;
	    } 
	  else 
	    {
	      // Tracklets found
	      AliDebug(3, "Tracklet found");
	      DecodeTracklet();
	      continue;
	    }
	} //if next HC

      if (fNextStatus == kSeekNonEoTracklet)
	{
	  //
	  // 2) Look for non-end_of_tracklet_marker
	  //
	  //printf("Word %d: 0x%08x\n", fWordCtr, *fDataWord); 

	  if ( *fDataWord != kEndoftrackletmarker ) 
	    {
	      fNextStatus = kDecodeHC;
	      AliDebug(3, "NON end-of-tracklet-marker found");
	      //// no do not continue - this should be the hcheader
	    }
	  else
	    {
	      //just go on and find the non-end_of_tracklet_marker
	      continue;
	    }
	}

      if ( fNextStatus == kDecodeHC )
	{
	  AliDebug(3, "Decode HC");

	  //
	  // 3) This Word must be Half Chamber Header
	  //
	  if ( (*fDataWord & 0x00000003) == 1 ) 
	    { // HC header
	      DecodeHCheader(fTimeBinsCalib); // This is the new header!
	      fDET    = fGeo->GetDetector(fLAYER, fSTACK, fSM);
	      fRowMax = fGeo->GetRowMax(fLAYER,fSTACK,fSM);
	      fColMax = fGeo->GetColMax(fROC);

	      if (fLastDET != fDET)
		{
		  AliDebug(4, "New DET!");	      
		  // allocate stuff for the new det
		  //man->ResetArrays();
		  digits = (AliTRDarrayADC *) man->GetDigits(fDET);
		  track0 = (AliTRDarrayDictionary *) man->GetDictionary(fDET,0);
		  track1 = (AliTRDarrayDictionary *) man->GetDictionary(fDET,1);
		  track2 = (AliTRDarrayDictionary *) man->GetDictionary(fDET,2);
		  
		  // Allocate memory space for the digits buffer
		  if (digits->GetNtime() == 0) 
		    {
		      AliDebug(4, "Allocating digits");	      
		      //AliDebug(5, Form("Alloc digits for det %d", det));
		      digits->Allocate(fRowMax, fColMax, fTBins);
		      if (man->UsesDictionaries())
			{
			  track0->Allocate(fRowMax, fColMax, fTBins);
			  track1->Allocate(fRowMax, fColMax, fTBins);
			  track2->Allocate(fRowMax, fColMax, fTBins);
			}
		    }
		  
		  indexes = man->GetIndexes(fDET);
		  indexes->SetSM(fSM);
		  indexes->SetStack(fSTACK);
		  indexes->SetLayer(fLAYER);
		  indexes->SetDetNumber(fDET);
		  
		  if (indexes->IsAllocated() == kFALSE)
		    {
		      AliDebug(4, "Allocating indexes");	      
		      indexes->Allocate(fRowMax, fColMax, fTBins);
		    }
		  fLastDET = fDET;
		}
	      
	      fMCMHctr2 = 0;
	      fHCdataCtr = 0;
	      
	      fChamberDone[fDET]++;
	      fNextStatus = kNextMCM;
	      AliDebug(3, "Decode HC OK");	      
	      continue;
	    } //HC header
	  else
	    {
	      AliDebug(3, "Decode HC NOT OK");	      
	      fNextStatus = kNextSM;
	      continue;
	    }
	} // if decode HC

      if (fNextStatus == kNextSM)
	{
	  
	  fDET     = 0;
	  fRetVal = 0;
	  fEqID     = 0;
	  fDataSize = 0;
	  fSizeOK = kFALSE;
	  
	  // After reading the first word check for size of this data and get Eq. ID
	  if ( fWordCtr == 1 ) 
	    {
	      fDataSize = fRawReader->GetDataSize()/4;  // Size of this payload in 32bit words
	      fEqID     = fRawReader->GetEquipmentId(); // Get Equipment ID
	      if ( fDataSize > 0 ) fSizeOK = kTRUE;
	    }
	  
	  // GTU Link Mask?
	  if ( (*fDataWord & 0xfffff000) ==  0xe0000000 ) 
	    {
	      DecodeGTUlinkMask();
	      fNextStatus = kNextHC;
	      continue;
	    } 
	  else 
	    {
	      AliWarning(Form("Equipment %d: First data word is not GTU Link Mask!", fEqID));
              fRawReader->AddMajorErrorLog(kGTULinkMaskMissing,Form("Equipment %d",fEqID));
	      fNextStatus = kStop;
	    }	    
	}// if nextSM

    } // not kStop

  AliDebug(1, Form("That's all folks! %d", fSM));
  //return kFALSE;
  return -1;
}

//============================================================================
// Decoding functions
//============================================================================


//____________________________________________________________________________
void AliTRDrawOldStream::DecodeHCheader(Int_t timeBins)
{
  //
  // Decode the HC header (fRawVersion == 2, 3, 4, ???)
  //

  fRVmajor = (*fDataWord >> 24) & 0x7f;
  fRVminor = (*fDataWord >> 17) & 0x7f;

  if (fRVmajor < 2 || fRVmajor > 4)
    AliError(Form(" Unsupported raw version: %d", fRawVersion))
  
  if ( fRawVersion != fRVmajor ) {
    
    AliWarning("===============================================================================");
    AliWarning(Form("Mismatch between fRawVersion (%d) and fRVmajor from HC header (%d)"
		    ,fRawVersion,fRVmajor));
    AliWarning(Form("Setting fRawVersion to %d", fRVmajor));
    AliWarning("===============================================================================");
    fRawVersion = fRVmajor;

  }

  //
  // check for zero suppression
  if ( fRawVersion >= 3 || fRawVersion <= 4 ) fZeroSuppressed = kTRUE;
  else                                        fZeroSuppressed = kFALSE;
  
  // 1st word (h[0])
  if ( (*fDataWord & 0x3) == 1 ) {

    fHCHWords = (*fDataWord >> 14) & 0x7;
    fSM       = (*fDataWord >>  9) & 0x1f;
    fLAYER    = (*fDataWord >>  6) & 0x7;
    fSTACK    = (*fDataWord >>  3) & 0x7;
    fSIDE     = (*fDataWord >>  2) & 0x1;

    fROC      = fGeo->GetDetectorSec(fLAYER, fSTACK);

    AliDebug(3, Form("0x%08x: HC header: sm=%d; roc=%d; side=%x", *fDataWord, fSM, fROC, fSIDE+10));

    if ((fSM    <  0) || 
        (fSM    > 17) || 
        (fLAYER <  0) || 
        (fLAYER >  5) || 
        (fSTACK <  0) || 
        (fSTACK >  4) || 
        (fSIDE  <  0) || 
        (fSIDE  >  1)) {
      AliWarning(Form("0x%08x: Strange HC header: dcs=%d; sm=%d; layer=%d; stack=%d.",
		     *fDataWord, fDCS, fSM, fLAYER, fSTACK));
      fRawReader->AddMajorErrorLog(kHCHeaderCorrupt,Form("0x%08x:dcs=%d; sm=%d; layer=%d; stack=%d.",
							 *fDataWord, fDCS, fSM, fLAYER, fSTACK));
    } 
    else {
      fHCHctr1++;
      fHCHctr2++;
    }
  } 
  else { 
    AliWarning(Form("0x%08x: No HC header when it was expected.", *fDataWord)); 
    fRawReader->AddMajorErrorLog(kHCHeaderMissing,Form("0x%08x", *fDataWord));
  }

  // 2nd word (h[1])
  if ( fHCHWords >= 1 ) {
    // read one more word
    if (NextData() != kWordOK)
      {
	AliWarning("Next HC word missing");
        fRawReader->AddMajorErrorLog(kHCWordMissing,"Next HC word missing"); 
	fNextStatus = kNextHC;
	return;
      }
    if ( (*fDataWord & 0x3) == 1 ) {
      
      fBCctr   =  (*fDataWord >> 16);
      fPTctr   =  (*fDataWord >> 12) & 0xf;
      fPTphase =  (*fDataWord >>  8) & 0xf;
      fTBins   = ((*fDataWord >>  2) & 0x3f) + 1;

      AliDebug(3, Form("0x%08x: HC header 2: BCctr=%d PTctr=%d PTph=%d TB=%d"
                      , *fDataWord, fBCctr, fPTctr, fPTphase, fTBins));

      if( fTBins != timeBins ) {

	AliWarning("===============================================================================");
	AliError(Form("Mismatch between nNTB from CDB (%d) and from HC header (%d)"
		      , timeBins, fTBins));
	AliWarning(Form("We will use the value from the raw data (HC header): %d", fTBins));
	AliWarning("===============================================================================");

	fTimeWords = (fTBins - 1)/3 + 1;	
      }

    }

  }

  // 3nd word (h[2])
  if ( fHCHWords >= 2 ) {
    // read one more word
    if (NextData() != kWordOK)
      {
	AliWarning("Next HC word missing");
        fRawReader->AddMajorErrorLog(kHCWordMissing,"Next HC word missing"); 
	fNextStatus = kNextHC;
	return;
      }
    if ( (*fDataWord & 0x3) == 1 ) {
       
      fTCon     = (*fDataWord >> 29) & 0x1;
      fPEDon    = (*fDataWord >> 31) & 0x1;
      fGAINon   = (*fDataWord >> 30) & 0x1;
      fXTon     = (*fDataWord >> 28) & 0x1;
      fNonLinOn = (*fDataWord >> 27) & 0x1;
      fBypass   = (*fDataWord >> 26) & 0x1;

      fCommonAdditive = (*fDataWord >> 20) & 0x3f;

      AliDebug(3, Form("0x%08x: HC header 3: TC=%d, PED=%d, GAIN=%d, XT=%d, NonLin=%d, Bypass=%d, Add=%d"
		      , fTCon, fPEDon, fGAINon, fXTon, fNonLinOn, fBypass, fCommonAdditive));
    }
  }

}  

//____________________________________________________________________________
void AliTRDrawOldStream::DecodeMCMheader()
{

  //
  // Decode the MCM header
  //

  if ( fRawVersion < 1 || fRawVersion > 3 ) 
    {
      AliError(Form(" Unsupported raw version: %d", fRawVersion));      
    }

  fMCM  = (*fDataWord & 0xff000000) >> 24;
  fEv   = (*fDataWord & 0x00fffff0) >> 4;

  fROB  = fMCM / 16;
  fMCM  = fMCM % 16;

  fROW  = AliTRDfeeParam::Instance()->GetPadRowFromMCM(fROB, fMCM);

  AliDebug(4, Form("0x%08x: SM%d L%dS%d. MCM Header: fROB=%d fMCM=%02d fEv=%02d"
		  , *fDataWord, fSM, fLAYER, fSTACK, fROB, fMCM, fEv));

  if ( fROB % 2 == 0 && fSIDE == 1 ) {
    AliWarning(Form("SM%d L%dS%d: Mismatch between fROB (%d) and fSIDE (%d): fMCM=%02d"
                   , fSM, fLAYER, fSTACK, fROB, fSIDE, fMCM ));
    fRawReader->AddMajorErrorLog(kROBSideMismatch,Form("SM%d L%dS%d: fROB (%d) fSIDE (%d): fMCM=%02d"
						       , fSM, fLAYER, fSTACK, fROB, fSIDE, fMCM ));
  }
  if ( fROB % 2 != 0 && fSIDE == 0 ) {
    AliWarning(Form("SM%d L%dS%d: Mismatch between fROB (%d) and fSIDE (%d): fMCM=%02d"
                   , fSM, fLAYER, fSTACK, fROB, fSIDE, fMCM ));
    fRawReader->AddMajorErrorLog(kROBSideMismatch,Form("SM%d L%dS%d: fROB (%d) fSIDE (%d): fMCM=%02d"
						       , fSM, fLAYER, fSTACK, fROB, fSIDE, fMCM ));
  }
  if ( (fSTACK == 2 && fROW >= fGeo->RowmaxC0()) ||
       (fSTACK != 2 && fROW >= fGeo->RowmaxC1()) || fROW < 0 ) {
    AliWarning(Form("SM%d L%dS%d: Wrong Padrow (%d) fROB=%d, fSIDE=%d, fMCM=%02d"
                   , fSM, fLAYER, fSTACK, fROW, fROB, fSIDE, fMCM ));
    fRawReader->AddMajorErrorLog(kWrongPadrow,Form("SM%d L%dS%d: Padrow (%d) fROB=%d, fSIDE=%d, fMCM=%02d"
						  , fSM, fLAYER, fSTACK, fROW, fROB, fSIDE, fMCM ));
  }
  
  fMCMHctr1++;
  fMCMHctr2++;

  // AdcMask for Zero supressed data
  if ( fRawVersion == 3 ) {
    // read one more word
    if (NextData() != kWordOK)
      {
	AliWarning("MCM ADC mask missing");
        fRawReader->AddMajorErrorLog(kMCMADCMaskMissing,"Missing"); 
	fNextStatus = kNextHC;
	return;
      }
    if ( (*fDataWord & 0x000007ff) == 0xC ) {     // at the moment bits 4-10 are empty
      
      for ( Int_t ctr = 0; ctr < fGeo->ADCmax(); ctr++ ) {
	if ( (*fDataWord >> (11+ctr)) == 0x1 ) fADCmask[ctr] = kTRUE;
	else                                  fADCmask[ctr] = kFALSE;
      }

      AliDebug(4, Form("0x%08x: ADC mask", *fDataWord));

    }
    else {
      AliWarning("Expected ADC mask but did not find one!");
      fRawReader->AddMajorErrorLog(kMCMADCMaskMissing,"Decode error"); 
    }

  }

}

//____________________________________________________________________________
void AliTRDrawOldStream::DecodeTracklet()
{

  //
  // Decode the Tracklet
  //
  // this function is not tested yet on real tracklets
  //

  if ( fRawVersion < 1 || fRawVersion > 3 ) 
    {
      AliError(Form(" Unsupported raw version: %d", fRawVersion));      
    }

  fTracklPID    = (*fDataWord >> 24) & 0xff;
  fTracklPadRow = (*fDataWord >> 20) & 0xf;    // 0:15
  fTracklDefL   = (*fDataWord >> 13) & 0x7f;
  fTracklPadPos = (*fDataWord)       & 0x1fff;

  fTracklPID    /= (Float_t)((1<<8) - 1);                      // 0:1 (steps of 0.39%)
  fTracklDefL    = (fTracklDefL  - ((1<< 7)-1)/2.) * 140.e-4;  // -0.889:0.889cm 
  fTracklPadPos  = (fTracklPadPos - ((1<<13)-1)/2.) * 160.e-4; // -65.528:65.528 cm

  //AliDebug(4, Form("0x%08x: Tracklet found: SM%d L%dS%d side %x: PadRow=%d PadPos=%f DefL=%f PID=%f"
  //		  , *fDataWord, fSM, fLAYER, fSTACK, fSIDE+10
  //                , fTracklPadRow, fTracklPadPos, fTracklDefL, fTracklPID));

  if( ((fSTACK == 2) && (fTracklPadRow >= (Int_t) fGeo->RowmaxC0())) ||
      ((fSTACK != 2) && (fTracklPadRow >= (Int_t) fGeo->RowmaxC1())) ) {
    AliWarning(Form("Strange Row read from Tracklet Word: %d", fTracklPadRow));
    fRawReader->AddMajorErrorLog(kTrackletRowMismatch,Form("Word: %d", fTracklPadRow));
  }

}

//____________________________________________________________________________
void AliTRDrawOldStream::DecodeGTUlinkMask()
{

  //
  // Decode the link masks sent by the GTU. These marke the active optical links
  // between GTU and Super Module. Up to now only fully active links are found
  // (0xfff = 12 active links).
  //

  if ( fRawVersion < 1 || fRawVersion > 3 ) 
    {
      AliError(Form(" Unsupported raw version: %d", fRawVersion));      
    }

  if ( fGTUctr1 == -1 ) fGTUctr2++;
  fGTUctr1++;

  if ( (fGTUctr1 >= 0) && (fGTUctr1 < 5) && (fGTUctr2 >= 0) && (fGTUctr2 < 18) ) {
    fGTUlinkMask[fGTUctr2][fGTUctr1] = (*fDataWord & 0xfff);
  }

}

//____________________________________________________________________________
Int_t  AliTRDrawOldStream::DecodeDataWord()
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
Int_t  AliTRDrawOldStream::DecodeDataWordV1V2()
{

  //
  // Decode the Data (full raw data. No zero suppression. 21 adc channels)
  //
  // return  0 means continue to next data word
  // return -1 means break data loop
  //

//   //  check the content first! - something wrong with that...
//   // Decode 32 bit data words with information from 3 time bins and copy the data
//   fSig[0] = (*fDataWord & 0x00000ffc) >> 2;
//   fSig[1] = (*fDataWord & 0x003ff000) >> 12;
//   fSig[2] = (*fDataWord & 0xffc00000) >> 22;  
//   if (fSig[0] <= 0 && fSig[1] <= 0 && fSig[2] <= 0)
//     return 0;

  if ( (*fDataWord & 0x00000003) != 0x2 && (*fDataWord & 0x00000003) != 0x3) {
    //AliWarning(Form("Data %08x : Data Word ends neither with b11 nor b10", (Int_t)*fDataWord));
    fRawReader->AddMinorErrorLog(kDataMaskError,Form("Data %08x", (Int_t)*fDataWord));
    return -1;
  }

  if ( (*fDataWord & 0x00000003) != fTbSwitch ) {    // Next ADC channel found
    fTbSwitch = (fTbSwitch & 2) | !(fTbSwitch & 1);   // 0x3 <--> 0x2
    fTbSwitchCtr = 0;
    fADC++;
    fTB=0;
  }

  fTbSwitchCtr++; // Just read one word

  // We have only timeTotal time bins
  if ( fTbSwitchCtr > fTimeWords ) {
    //AliWarning(Form("Data is strange. Already found %d words for this ADC channel", (Int_t)fTbSwitchCtr));
    fRawReader->AddMinorErrorLog(kADCNumberOverflow,Form("%d words", (Int_t)fTbSwitchCtr));
    return 0;
  }

  // We have only 21 ADC channels.
  if ( fADC > (Int_t)fGeo->ADCmax()-1 ) {
    //AliWarning(Form("Data %08x : Data is strange. fADC is already %d", (Int_t)*fDataWord, (Int_t)fADC));
    fRawReader->AddMinorErrorLog(kADCChannelOverflow,Form("Data %08x : fADC=%d", (Int_t)*fDataWord, (Int_t)fADC));
    return 0;
  }

  // There are 18 pads connected to each MCM ADC channels 2...19. The other channels cross to other
  // MCMs and are good for online tracking in the MCM.
  if ( fADC > 1 && fADC < (Int_t)fGeo->ADCmax()-1 ) {

    // Get Pad column
    fCOL = AliTRDfeeParam::Instance()->GetPadColFromADC(fROB, fMCM, fADC);

    // We have only 144 Pad Columns
    //if ( fCOL > fColMax-1 || fCOL < 0 ) {
    if ( fCOL >= 0 && fCOL < fColMax && fROW >= 0 && fROW < fRowMax ) 
      {
	// Decode 32 bit data words with information from 3 time bins and copy the data
	fSig[0] = (*fDataWord & 0x00000ffc) >> 2;
	fSig[1] = (*fDataWord & 0x003ff000) >> 12;
	fSig[2] = (*fDataWord & 0xffc00000) >> 22;
	
	if (fSig[0] > 0 || fSig[1] > 0 || fSig[2] > 0)
	  return 1;
	else
	  return 0;
      }
    else
      {
// 	AliWarning(Form("SM%d L%dS%d: Wrong Pad column (%d) fROB=%d, fSIDE=%d, fMCM=%02d", fSM,
// 			fLAYER, fSTACK, fCOL, fROB, fSIDE, fMCM ));
	fRawReader->AddMajorErrorLog(kWrongPadcolumn,Form("SM%d L%dS%d: column (%d) fROB=%d, fSIDE=%d, fMCM=%02d", fSM,
							  fLAYER, fSTACK, fCOL, fROB, fSIDE, fMCM ));
	return 0;
      }
    // Print data to screen:
    // Do NOT switch on for default production, it is VERY slow
    //    AliDebug(5, Form("SM%d L%dS%d: ROB%d MCM=%d ADC=%d (ROW=%d COL=%d): Data %04d %04d %04d\n",
    //		     fSM, fLAYER, fSTACK, fROB, fMCM, fADC, fROW, fCOL, fSig[0], fSig[1], fSig[2]));
    
  }
  else {
    
    fCOL = -1;
    return 0;
  }

  return 1;

}

//____________________________________________________________________________
Int_t  AliTRDrawOldStream::DecodeDataWordV3()
{

  //
  // Decode the data (Zero suppresses data. 21 adc channels)
  //
  // return  0 means continue to next data word
  // return -1 means break data loop
  //
  // NOT TESTED YET!!!!!!!!
  //

  if ( (*fDataWord & 0x00000003) != 0x2 && (*fDataWord & 0x00000003) != 0x3) {
    AliWarning(Form("Data %08x : Data Word ends neither with b11 nor b10", (Int_t)*fDataWord));
    fRawReader->AddMinorErrorLog(kDataMaskError,Form("Data %08x", (Int_t)*fDataWord));
    return -1;
  }

  if ( (*fDataWord & 0x00000003) != fTbSwitch ) {    // Next ADC channel found
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
    AliWarning(Form("Data is strange. Already found %d words for this ADC channel", (Int_t)fTbSwitchCtr));
    fRawReader->AddMinorErrorLog(kADCNumberOverflow,Form("%d words", (Int_t)fTbSwitchCtr));
    return 0;
  }

  // We have only 21 ADC channels.
  if ( fADC > (Int_t)fGeo->ADCmax()-1 ) {
    AliWarning(Form("Data %08x : Data is strange. fADC is already %d", (Int_t)*fDataWord, (Int_t)fADC));
    fRawReader->AddMinorErrorLog(kADCChannelOverflow,Form("Data %08x : fADC=%d", (Int_t)*fDataWord, (Int_t)fADC));
    return 0;
  }

  // There are 18 pads connected to each MCM ADC channels 2...19. The other channels cross to other
  // MCMs and are good for online tracking in the MCM.
  if ( fADC > 1 && fADC < (Int_t)fGeo->ADCmax()-1 ) {

    // Get Pad column
    //fCOL = fGeo->GetPadColFromADC(fROB, fMCM, fADC);
    fCOL = AliTRDfeeParam::Instance()->GetPadColFromADC(fROB, fMCM, fADC);

    // We have only 144 Pad Columns
    if ( fCOL > fColMax-1 || fCOL < 0 ) {
      AliWarning(Form("SM%d L%dS%d: Wrong Pad column (%d) fROB=%d, fSIDE=%d, fMCM=%02d", fSM,
		    fLAYER, fSTACK, fCOL, fROB, fSIDE, fMCM ));
      fRawReader->AddMajorErrorLog(kWrongPadcolumn,Form("SM%d L%dS%d: column (%d) fROB=%d, fSIDE=%d, fMCM=%02d", fSM,
		    fLAYER, fSTACK, fCOL, fROB, fSIDE, fMCM ));
    }

    // Decode 32 bit data words with information from 3 time bins and copy the data
    fSig[0] = (*fDataWord & 0x00000ffc) >> 2;
    fSig[1] = (*fDataWord & 0x003ff000) >> 12;
    fSig[2] = (*fDataWord & 0xffc00000) >> 22;

    // Print data to screen:
    AliDebug(5, Form("SM%d L%dS%d: ROB%d MCM=%d ADC=%d (ROW=%d COL=%d): Data %04d %04d %04d\n",
		     fSM, fLAYER, fSTACK, fROB, fMCM, fADC, fROW, fCOL, fSig[0], fSig[1], fSig[2]));
    
  }
  else {
    
    fCOL = -1;
    
  }

  return 1;

}
