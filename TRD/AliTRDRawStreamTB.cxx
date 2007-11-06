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
#include "AliTRDRawStreamTB.h"
#include "AliTRDgeometry.h"
#include "AliTRDcalibDB.h"
#include "AliTRDfeeParam.h"

#include "AliTRDdigitsManager.h"
#include "AliTRDdataArrayI.h"
#include "AliTRDSignalIndex.h"

#include "TMath.h"

using namespace AliTRDrawDataUtilsTB;

ClassImp(AliTRDRawStreamTB)

Bool_t AliTRDRawStreamTB::fgStackIndexBug=kFALSE;
Int_t  AliTRDRawStreamTB::fgForceRawVersion = -1;
Bool_t AliTRDRawStreamTB::fgSupressWarnings=kFALSE;
Bool_t AliTRDRawStreamTB::fgRawDataHack=kFALSE;
Bool_t AliTRDRawStreamTB::fgExtraDebug=kFALSE;

//_____________________________________________________________________________
AliTRDRawStreamTB::AliTRDRawStreamTB() 
  :TObject()
  ,fSig()
  ,fADC(0)
  ,fMaxADCgeom(-1)
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
  ,fDET(-1)
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
  ,fMCMWordCrt(0)
  ,fMCMWordsExpected(0)
  ,fRawReader(NULL)
  ,fRawVersion(0)
  ,fRawDigitThreshold(0)
  ,fNextStatus(0)
  ,fLastStatus(0)
  ,fTbSwitch(0)
  ,fTbSwitchCtr(0)
  ,fTimeWords(0)
  ,fWordCtr(0)
  ,fRowMax(0)
  ,fColMax(0)
  ,fADCmask()
  ,fLastADCmask(0)
  ,fChamberDone()
  ,fRetVal(0)
  ,fEqID(0)
  ,fDataSize(0)
  ,fSizeOK(kFALSE)
  ,fCountBytes(0)
  ,fBufSize(0)
  ,fkBufferSet(kFALSE)
  ,fPos(NULL)
  ,fDataWord(NULL)
  ,fTimeBinsCalib(0)
  ,fADClookup()
  ,fNActiveADCs(0)
  ,fEndOfDataFlag(kFALSE)
  ,fHCInfo()
  ,fMCMInfo()
  ,fiSMx(0)
  ,fSharedPadsOn(kFALSE)
  ,fIsPadShared(kFALSE)
  ,fGeo(NULL) 
{
  //
  // Default constructor
  //

  for (Int_t i = 0; i < 540; i++) {
    fChamberDone[i] = 0;
  }

}

//_____________________________________________________________________________
AliTRDRawStreamTB::AliTRDRawStreamTB(AliRawReader *rawReader) 
  :TObject()
  ,fSig()
  ,fADC(0)
  ,fMaxADCgeom(-1)
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
  ,fDET(-1)
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
  ,fMCMWordCrt(0)
  ,fMCMWordsExpected(0)
  ,fRawReader(rawReader)
  ,fRawVersion(0)
  ,fRawDigitThreshold(0)
  ,fNextStatus(0)
  ,fLastStatus(0)
  ,fTbSwitch(0)
  ,fTbSwitchCtr(0)
  ,fTimeWords(0)
  ,fWordCtr(0)
  ,fRowMax(0)
  ,fColMax(0)
  ,fADCmask()
  ,fLastADCmask(0)
  ,fChamberDone()
  ,fRetVal(0)
  ,fEqID(0)
  ,fDataSize(0)
  ,fSizeOK(kFALSE)
  ,fCountBytes(0)
  ,fBufSize(0)
  ,fkBufferSet(kFALSE)
  ,fPos(NULL)
  ,fDataWord(NULL)
  ,fTimeBinsCalib(0)
  ,fADClookup()
  ,fNActiveADCs(0)
  ,fEndOfDataFlag(kFALSE)
  ,fHCInfo()
  ,fMCMInfo()
  ,fiSMx(0)
  ,fSharedPadsOn(kFALSE)
  ,fIsPadShared(kFALSE)
  ,fGeo(NULL) 
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
AliTRDRawStreamTB::AliTRDRawStreamTB(const AliTRDRawStreamTB& stream)
  :TObject(stream)
  ,fSig()
  ,fADC(-1)
  ,fMaxADCgeom(-1)
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
  ,fDET(-1)
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
  ,fMCMWordCrt(0)
  ,fMCMWordsExpected(0)
  ,fRawReader(NULL)
  ,fRawVersion(0)
  ,fRawDigitThreshold(0)
  ,fNextStatus(0)
  ,fLastStatus(0)
  ,fTbSwitch(0)
  ,fTbSwitchCtr(0)
  ,fTimeWords(0)
  ,fWordCtr(0)
  ,fRowMax(-1)
  ,fColMax(-1)
  ,fADCmask()
  ,fLastADCmask(0)
  ,fChamberDone()
  ,fRetVal(0)
  ,fEqID(0)
  ,fDataSize(0)
  ,fSizeOK(kFALSE)
  ,fCountBytes(0)
  ,fBufSize(0)
  ,fkBufferSet(kFALSE)
  ,fPos(NULL)
  ,fDataWord(NULL)
  ,fTimeBinsCalib(0)
  ,fADClookup()
  ,fNActiveADCs(0)
  ,fEndOfDataFlag(kFALSE)
  ,fHCInfo()
  ,fMCMInfo()
  ,fiSMx(0)
  ,fSharedPadsOn(kFALSE)
  ,fIsPadShared(kFALSE)
  ,fGeo(NULL)
{
  //
  // Copy constructor
  //

  AliFatal("Copy constructor not implemented");

}

//_____________________________________________________________________________
AliTRDRawStreamTB& AliTRDRawStreamTB::operator = (const AliTRDRawStreamTB& 
					      /* stream */)
{
  //
  // Assigment operator
  //

  Fatal("operator =", "assignment operator not implemented");
  return *this;

}

//_____________________________________________________________________________
AliTRDRawStreamTB::~AliTRDRawStreamTB()
{
  //
  // Destructor
  //

  if (fGeo) 
    {  
      delete fGeo;
    }

}

//_____________________________________________________________________________
void AliTRDRawStreamTB::SetRawReader(AliRawReader *rawReader) 
{
  //
  // Set the rawreader
  //

  if (rawReader)
    {
      fRawReader = rawReader;
    }
}

//_____________________________________________________________________________
Bool_t AliTRDRawStreamTB::SetRawVersion(Int_t rv)
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
Int_t AliTRDRawStreamTB::Init()
{
  //
  // Initialization
  //

  if (!AliTRDcalibDB::Instance()) {
    AliError("Could not get calibration object");
    return 0;
  }

  if (!fGeo) {
    fGeo = new AliTRDgeometry();
  }

  fMaxADCgeom = (Int_t)fGeo->ADCmax();
  
  fTimeBinsCalib = AliTRDcalibDB::Instance()->GetNumberOfTimeBins();
  //AliDebug(2, Form("Number of Timebins read from CDB: %d", fTimeBinsCalib));

  // The number of data words needed for this number of time bins (there
  // are 3 time bins in one word)
  fTimeWords = (fTimeBinsCalib-1)/3 + 1;

  fTbSwitch    = 3;
  fTbSwitchCtr = 0;

  fHCHctr1 = fHCHctr2 =  0;
  fGTUctr1 = fGTUctr2 = -1;

  fHCdataCtr = 0;
  fWordCtr   = 0;  

  fDET     =  -1;
  fLastDET =  -1;
  fRetVal = 0;
  fEqID     = 0;
  fDataSize = 0;
  fSizeOK = kFALSE;
  
  fLastStatus = kStart;
  fNextStatus = kStart;

  fCountBytes = 0;
  fBufSize = 0;
  fDataWord = NULL;
  fPos = NULL;
  fWordCtr = 0;
  fkBufferSet = kFALSE;

  fMCMWordCrt = 0;
  fMCMWordsExpected = 0;

  fEndOfDataFlag = kFALSE;
  // set all ADC active
  fNActiveADCs = ChannelsToRead(0x1fffff); // should be 1111 1111 1111 1111 1111 1 = 21 bits active (0-20)

  fLastADCmask = 0;

  //default value overwritten by HC header word h[2]
  fCommonAdditive = 10;
  //fSharedPadsOn = kFALSE;
  //fSharedPadsOn = kTRUE;
  fIsPadShared = kFALSE;
  fZeroSuppressed = kFALSE;

  return kTRUE;
}

//____________________________________________________________________________
void AliTRDRawStreamTB::SwapOnEndian()
{
  //
  // Check the endian and swap if needed
  //

  int itemp = 1;
  char* ptemp = (char*) &itemp;
  if (ptemp[0] != 1)
    {
      // need to swap...
      // assume we are at the begining of the buffer!
      //AliDebug(5, "Swapping.");
      UInt_t *pbegin = (UInt_t*)fPos;
      UInt_t iutmp = 0;
      for (UInt_t i = 0; i < fBufSize / AliTRDrawDataUtilsTB::kSizeWord; i++)
	{
	  fDataWord = pbegin + i;
	  iutmp = (((*fDataWord & 0x000000ffU) << 24) | ((*fDataWord & 0x0000ff00U) <<  8) |
		   ((*fDataWord & 0x00ff0000U) >>  8) | ((*fDataWord & 0xff000000U) >> 24));
	  // here we override the value in the buffer!
	  *fDataWord = iutmp;
	}
      fDataWord = pbegin;
    }
}

//____________________________________________________________________________
Int_t AliTRDRawStreamTB::NextData()
{
  //
  // Updates the next data word pointer
  //

  if (fCountBytes + kSizeWord >= fBufSize)
    {
      fkBufferSet = fRawReader->ReadNextData(fPos);
      if (fkBufferSet == kTRUE)
	{
	  fBufSize = fRawReader->GetDataSize();
	  fCountBytes = 0;	  
	  fDataWord = (UInt_t*)fPos;
	  SwapOnEndian();
	  if (fgRawDataHack == kTRUE)
	    {
	      SkipWords(24);
	      fCountBytes += 24 * kSizeWord;
	      //fgRawDataHack = kFALSE;
	    }
	  ChangeStatus(kNextSM);
	  fWordCtr = 0;
	  AliDebug(3, "NextSM. Buffer is set.");
	  return kNextSM;
	}
      else
	{
	  AliDebug(3, "No more data!");
	  //ChangeStatus(kStop);
	  ChangeStatus(kNoMoreData);
	  return kNoMoreData;
	}
    }
  else
    {
      fPos += kSizeWord;
      fCountBytes += kSizeWord;	  
      fDataWord = (UInt_t*)fPos;
      fWordCtr++;
      AliDebug(10, Form("Current word %d : 0x%x at 0x%x. Count bytes %d of %d", 
			fWordCtr, *fDataWord, fPos, fCountBytes, fBufSize));
      //       if (fCountBytes > 3080000)
      // 	printf("Current word %d : 0x%x at 0x%x. Count bytes %d of %d\n", 
      // 	       fWordCtr, *fDataWord, fPos, fCountBytes, fBufSize);
      return kWordOK;
    }
}

//____________________________________________________________________________
Int_t AliTRDRawStreamTB::RewindWord()
{
  //
  // Updates the next data word pointer
  //
  
  if (fkBufferSet == kFALSE)
    return kStart;

  if (fCountBytes > 0)
    {
      fPos -= kSizeWord;
      fCountBytes -= kSizeWord;	  
      fDataWord = (UInt_t*)fPos;
      fWordCtr--;      
    }
  else
    {
      AliWarning("Already at the begining of the buffer");
    }
  
  return kWordOK;
}

//============================================================================
// Decoding functions
//============================================================================

//____________________________________________________________________________
void AliTRDRawStreamTB::DecodeHCheader(Int_t timeBins)
{
  //
  // Decode the HC header (fRawVersion == 2, 3, 4, ???)
  //
  //AliDebug(3, "Here");
    
  // 1st word (h[0])
  if ( (*fDataWord & 0x3) == 1 ) 
    {
      
      fHCInfo.DecodeH0(fDataWord);
      if (fgExtraDebug)
	fHCInfo.Dump();
      
      fHCHWords = fHCInfo.fNExtraWords;
      fSM       = fHCInfo.fSM;
      fLAYER    = fHCInfo.fLayer;
      fSTACK    = fHCInfo.fStack;
      fSIDE     = fHCInfo.fSide;

      fRVmajor  = fHCInfo.fRawVMajor;
      fRVminor  = fHCInfo.fRawVMinor;
      
      //assuming valid raw version w/o ZS
      if (fRVmajor <= 0)
	fRVmajor = 2;

      if ( fRawVersion != fRVmajor ) 
	{
	  if (fgSupressWarnings == kTRUE)
	    {
	      if (fgForceRawVersion > 0)
		{
		  fRawVersion = fgForceRawVersion;
		}
	      else
		{
		  fRawVersion = fRVmajor;	
		}      
	    }
	  else
	    {
	      AliWarning("===============================================================================");
	      AliWarning(Form("Mismatch between fRawVersion (%d) and fRVmajor from HC header (%d)"
			      ,fRawVersion,fRVmajor));
	      if (fgForceRawVersion > 0)
		{
		  AliWarning(Form("Forcing fRawVersion to %d", fgForceRawVersion));
		  fRawVersion = fgForceRawVersion;
		}
	      else
		{
		  AliWarning(Form("Setting fRawVersion to %d", fRVmajor));
		  fRawVersion = fRVmajor;	
		}
	      AliWarning("===============================================================================");
	    } //warn or not
	}//if raw version missmatch

      //AliInfo(Form("Raw Version %d", fRawVersion));
      
      //
      // check for zero suppression
      if ( fRawVersion >= 3 || fRawVersion <= 4 ) fZeroSuppressed = kTRUE;
      else                                        fZeroSuppressed = kFALSE;
    
      fROC      = fGeo->GetDetectorSec(fLAYER, fSTACK);
      
      AliDebug(3, Form("0x%08x: HC header: sm=%d; roc=%d; side=%x", *fDataWord, fSM, fROC, fSIDE+10));
      AliDebug(5, Form("0x%08x: HC header: expecting %d HC words", *fDataWord, fHCHWords));

      if ((fSM    <  0) || 
	  (fSM    > 17) || 
	  (fLAYER <  0) || 
	  (fLAYER >  5) || 
	  (fSTACK <  0) || 
	  (fSTACK >  4) || 
	  (fSIDE  <  0) || 
	  (fSIDE  >  1)) 
	{
	  AliWarning(Form("0x%08x: Strange HC header: dcs=%d; sm=%d; layer=%d; stack=%d.",
			  *fDataWord, fDCS, fSM, fLAYER, fSTACK));
	  fRawReader->AddMajorErrorLog(kHCHeaderCorrupt,Form("0x%08x:dcs=%d; sm=%d; layer=%d; stack=%d.",
							     *fDataWord, fDCS, fSM, fLAYER, fSTACK));
	} 
      else 
	{
	  fHCHctr1++;
	  fHCHctr2++;
	}
    } 
  else 
    { 
      AliWarning(Form("0x%08x: No HC header when it was expected.", *fDataWord)); 
      fRawReader->AddMajorErrorLog(kHCHeaderMissing,Form("0x%08x", *fDataWord));
    }

  // 2nd word (h[1])
  if ( fHCHWords >= 1 ) 
    {
      // read one more word
      if (NextData() != kWordOK)
	{
	  AliWarning("Next HC ( H[1] ) word missing");
	  fRawReader->AddMajorErrorLog(kHCWordMissing,"Next HC ( H[1] )word missing"); 
	  ChangeStatus(kNextHC);
	  return;
	}

      if ( (*fDataWord & 0x3) == 1 ) 
	{

	  fHCInfo.DecodeH1(fDataWord);
 	  
	  fBCctr   = fHCInfo.fBunchCrossCounter;
 	  fPTctr   = fHCInfo.fPreTriggerCounter;
 	  fPTphase = fHCInfo.fPreTriggerPhase;
 	  fTBins   = fHCInfo.fTimeBins;
	  
	  fTimeWords = (fTBins - 1)/3 + 1;	
	  
 	  AliDebug(3, Form("0x%08x: HC header [1]: BCctr=%d PTctr=%d PTph=%d TB=%d"
 			   , *fDataWord, fBCctr, fPTctr, fPTphase, fTBins));

	  if( fTBins != timeBins ) 
	    {	      
	      if (fgSupressWarnings == kTRUE)
		{
		  fTimeWords = (fTBins - 1)/3 + 1;	
		}
	      else
		{
		  AliWarning("===============================================================================");
		  AliError(Form("Mismatch between nNTB from CDB (%d) and from HC header (%d)"
				, timeBins, fTBins));
		  AliWarning(Form("We will use the value from the raw data (HC header): %d", fTBins));
		  AliWarning("===============================================================================");   		  
		  fTimeWords = (fTBins - 1)/3 + 1;	
		}
	    }
	}      
    }

  // 3nd word (h[2])
  if ( fHCHWords >= 2 ) 
    {
      // read one more word
      if (NextData() != kWordOK)
	{
	  AliWarning("Next HC ( [2] )word missing");
	  fRawReader->AddMajorErrorLog(kHCWordMissing,"Next HC ( [2] ) word missing"); 
	  ChangeStatus(kNextHC);
	  return;
	}

      if ( (*fDataWord & 0x3) == 1 ) 
	{
       
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

  if (fTBins <= 0)
    {
      fTBins = 30;
      fTimeWords = (fTBins - 1)/3 + 1;	
      AliDebug(5, Form("Forced tbins %d timewords %d", fTBins, fTimeWords));
    }
}  

//____________________________________________________________________________
Int_t AliTRDRawStreamTB::ChannelsToRead(Int_t ADCmask)
{
  //AliDebug(3, "Here");
  memset(fADClookup, -1, 32 * sizeof(Int_t));
  fADClookup[0] = 0; // count entries
  fADClookup[1] = 2; // index - data start at 2
  UInt_t mask = 0;
  for (Int_t i = 0; i < 30; i++)
    {
      mask = 1 << i;
      if ((ADCmask & mask))
	{
	  //AliDebug(9, Form("fDataWord=0x%08x mask=0x%08x i=%d", *fDataWord, mask, i));
	  fADClookup[fADClookup[1]] = i;
	  ++fADClookup[0];
	  ++fADClookup[1];
	}
    }

  // test the iteration - comment out for production
  // begin of comment out section
//   char schannels[512];
//   sprintf(schannels, "ADC Channels to read: ");
//   fADClookup[1] = 2;
//   while(fADClookup[1] - 2 < fADClookup[0])
//     {
//       AliDebug(9, Form("max=%d index=%d adc=%d", fADClookup[0], fADClookup[1], fADClookup[fADClookup[1]]));
//       strcat(schannels, Form("%d ", fADClookup[fADClookup[1]]));
//       fADClookup[1]++;
//     }
  //AliDebug(9, Form("%s", schannels));
  //AliDebug(9, Form("ADC channels = %d", fADClookup[0]));
  // end of comment out section

  fADClookup[1] = 2;
  return fADClookup[0];
}

//____________________________________________________________________________
void AliTRDRawStreamTB::DecodeTracklet()
{

  //
  // Decode the Tracklet
  //
  // this function is not tested yet on real tracklets
  //
  //AliDebug(3, "Here");
//   if ( fRawVersion < 1 || fRawVersion > 3 ) 
//     {
//       AliError(Form(" Unsupported raw version: %d", fRawVersion));      
//     }

  fTracklPID    = (*fDataWord >> 24) & 0xff;
  fTracklPadRow = (*fDataWord >> 20) & 0xf;    // 0:15
  fTracklDefL   = (*fDataWord >> 13) & 0x7f;
  fTracklPadPos = (*fDataWord)       & 0x1fff;

  fTracklPID    /= (Float_t)((1<<8) - 1);                      // 0:1 (steps of 0.39%)
  fTracklDefL    = (fTracklDefL  - ((1<< 7)-1)/2.) * 140.e-4;  // -0.889:0.889cm 
  fTracklPadPos  = (fTracklPadPos - ((1<<13)-1)/2.) * 160.e-4; // -65.528:65.528 cm

  AliDebug(4, Form("0x%08x: Tracklet found: SM%d L%dS%d side %x: PadRow=%d PadPos=%f DefL=%f PID=%f"
		   , *fDataWord, fSM, fLAYER, fSTACK, fSIDE+10
		   , fTracklPadRow, fTracklPadPos, fTracklDefL, fTracklPID));

  if( (fSTACK == 2) && (fTracklPadRow >= (Int_t) fGeo->RowmaxC0()) ||
      (fSTACK != 2) && (fTracklPadRow >= (Int_t) fGeo->RowmaxC1()) ) {
    AliWarning(Form("Strange Row read from Tracklet Word: %d", fTracklPadRow));
    fRawReader->AddMajorErrorLog(kTrackletRowMismatch,Form("Word: %d", fTracklPadRow));
  }

}

//____________________________________________________________________________
void AliTRDRawStreamTB::DecodeMCMheader()
{

  //
  // Decode the MCM header
  //
  //AliDebug(3, "Here");

  if ( fRawVersion < 0 || fRawVersion > 3 ) 
    {
      AliError(Form(" Unsupported raw version: %d", fRawVersion));      
    }

  fMCMInfo.Decode(fDataWord);
  
  fMCM = fMCMInfo.fMCM;
  fROB = fMCMInfo.fROB;
  fEv  = fMCMInfo.fEvCounter;
  if (fgExtraDebug)
    fMCMInfo.Dump();

  fROW  = AliTRDfeeParam::Instance()->GetPadRowFromMCM(fROB, fMCM);
  
  AliDebug(4, Form("0x%08x: SM%d L%dS%d. MCM Header: fROB=%d fMCM=%02d fEv=%02d"
		   , *fDataWord, fSM, fLAYER, fSTACK, fROB, fMCM, fEv));

  if ( fROB % 2 == 0 && fSIDE == 1 ) 
    {
      AliWarning(Form("SM%d L%dS%d: Mismatch between fROB (%d) and fSIDE (%d): fMCM=%02d"
		      , fSM, fLAYER, fSTACK, fROB, fSIDE, fMCM ));
      fRawReader->AddMajorErrorLog(kROBSideMismatch,Form("SM%d L%dS%d: fROB (%d) fSIDE (%d): fMCM=%02d"
							 , fSM, fLAYER, fSTACK, fROB, fSIDE, fMCM ));
    }

  if ( fROB % 2 != 0 && fSIDE == 0 ) 
    {
      AliWarning(Form("SM%d L%dS%d: Mismatch between fROB (%d) and fSIDE (%d): fMCM=%02d"
		      , fSM, fLAYER, fSTACK, fROB, fSIDE, fMCM ));
      fRawReader->AddMajorErrorLog(kROBSideMismatch,Form("SM%d L%dS%d: fROB (%d) fSIDE (%d): fMCM=%02d"
							 , fSM, fLAYER, fSTACK, fROB, fSIDE, fMCM ));
    }

  if ( (fSTACK == 2 && fROW >= fGeo->RowmaxC0()) ||
       (fSTACK != 2 && fROW >= fGeo->RowmaxC1()) || fROW < 0 ) 
    {
      AliWarning(Form("SM%d L%dS%d: Wrong Padrow (%d) fROB=%d, fSIDE=%d, fMCM=%02d"
		      , fSM, fLAYER, fSTACK, fROW, fROB, fSIDE, fMCM ));
      fRawReader->AddMajorErrorLog(kWrongPadrow,Form("SM%d L%dS%d: Padrow (%d) fROB=%d, fSIDE=%d, fMCM=%02d"
						     , fSM, fLAYER, fSTACK, fROW, fROB, fSIDE, fMCM ));
    }
  
  fMCMHctr1++;
  fMCMHctr2++;

  fMCMWordCrt = 1; // MCM header

  // AdcMask for Zero supressed data
  if ( fRawVersion == 3 ) 
    {
      // read one more word
      if (NextData() != kWordOK)
	{
	  AliWarning("MCM ADC mask missing");
	  fRawReader->AddMajorErrorLog(kMCMADCMaskMissing,"Missing"); 
	  fNextStatus = kNextHC;
	  return;
	}
      else
	{
	  ++fMCMWordCrt;

	  for ( Int_t ctr = 0; ctr < fMaxADCgeom; ctr++ ) {
	    if ( (*fDataWord >> (11+ctr)) == 0x1 ) fADCmask[ctr] = kTRUE;
	    else                                  fADCmask[ctr] = kFALSE;
	  }

	  if (*fDataWord & 0xf != 0xc)
	    {
	      AliWarning(Form("ADC mask does not end with 0xc : 0x%x", *fDataWord));
	    }

	  //AliDebug(4, Form("0x%08x: ADC mask", *fDataWord));
	  // 7 MSbits are ignored!
	  UInt_t maskWord = (*fDataWord >> 4) & 0x1fffff;
	  fNActiveADCs = ChannelsToRead(maskWord);	  
	}
    }

  if (fRawVersion <= 2)
    {
      fNActiveADCs = ChannelsToRead(0x1fffff); // should be 1111 1111 1111 1111 1111 1 = 21 bits active (0-20)
      //directly get it like that:
      fMCMWordsExpected = 1 + fNActiveADCs * fTBins / 3;      
    }

  if (fRawVersion >= 3)
    {
      // raw version 3:
      //directly get it like that:
      fMCMWordsExpected = 1 + 1 + (fTBins * fNActiveADCs) / 3;
    }
  
  AliDebug(8, Form("TBins %d fNActiveADCs %d => We expect %d MCM words. We read %d so far.", fTBins, fNActiveADCs, fMCMWordsExpected, fMCMWordCrt));
}

//____________________________________________________________________________
Bool_t AliTRDRawStreamTB::DecodeADCWord()
{
  // Decode ADC word for any pad
  // Rearange this function!!!
  //AliDebug(3, "Here");
  Bool_t kIsDataOK = kFALSE;
  
  // Get Pad column
  // fCOL = fFee->GetPadColFromADC(fROB, fMCM, fADC);
  if ( fADC >= fMaxADCgeom - 1)
    {
      // let us guess the Column
      // take the one before last ADC and shift by one column
      // later we check if we are inside the limits of the chamber
      fCOL = AliTRDfeeParam::Instance()->GetPadColFromADC(fROB, fMCM, fADC - 1);
      fCOL--;
      //AliDebug(8, Form("-x fADC %d fCOL %d", fADC, fCOL));
      if (fCOL >= fColMax || fCOL < 0)
	{
	  //AliDebug(8, Form("-xx fADC %d fCOL %d", fADC, fCOL));
	  fCOL = -1;
	  kIsDataOK = kFALSE;
	  return kIsDataOK;
	}
    }
  else
    {
      //AliDebug(8, Form("-y fADC %d fCOL %d", fADC, fCOL));
      fCOL = AliTRDfeeParam::Instance()->GetPadColFromADC(fROB, fMCM, fADC);
      if (fCOL >= fColMax || fCOL < 0)
	{
	  //AliDebug(8, Form("-xx fADC %d fCOL %d", fADC, fCOL));
	  fCOL = -1;
	  kIsDataOK = kFALSE;
	  return kIsDataOK;
	}
    }

  // We have only 144 Pad Columns
  //if ( fCOL > fColMax-1 || fCOL < 0 ) 
  if ( fCOL >= 0 && fCOL < fColMax && fROW >= 0 && fROW < fRowMax ) 
    {
      // Decode 32 bit data words with information from 3 time bins and copy the data
      fSig[0] = (*fDataWord & 0x00000ffc) >> 2;
      fSig[1] = (*fDataWord & 0x003ff000) >> 12;
      fSig[2] = (*fDataWord & 0xffc00000) >> 22;
      
      // Print data to screen:
      AliDebug(5, Form("DATA : 0x%x tTbin %d", *fDataWord, fTB));
      AliDebug(5, Form("SM%d L%dS%d: ROB%d MCM=%d ADC=%d (ROW=%d COL=%d): Data %04d %04d %04d\n",
		       fSM, fLAYER, fSTACK, fROB, fMCM, fADC, fROW, fCOL, fSig[0], fSig[1], fSig[2]));	      
      kIsDataOK = kTRUE;
    }
  else
    {
      AliWarning(Form("SM%d L%dS%d: Wrong Pad column (%d) fROB=%d, fSIDE=%d, fMCM=%02d", fSM,
		      fLAYER, fSTACK, fCOL, fROB, fSIDE, fMCM ));
      fRawReader->AddMajorErrorLog(kWrongPadcolumn,Form("SM%d L%dS%d: column (%d) fROB=%d, fSIDE=%d, fMCM=%02d", fSM,
							fLAYER, fSTACK, fCOL, fROB, fSIDE, fMCM ));
      kIsDataOK = kFALSE;
    }
  
  return kIsDataOK;
}

//____________________________________________________________________________
Bool_t AliTRDRawStreamTB::DecodeNextRawWord()
{
  //AliDebug(8, Form("-----------------------------------------"));
  //AliDebug(8, Form("DATA IS 0x%x", *fDataWord));
  
  //AliDebug(3, "Here");
  if ( *fDataWord == kEndOfRawDataMarker ) 
    {  // End of half-chamber data, finished
      AliDebug(3, "We have reached the eof data. Next should be HC header.");
      fGTUctr1 = -1;
      ChangeStatus(kNextHC);
      fEndOfDataFlag = kTRUE;
      return kFALSE;
    }

  if (fADClookup[1] - 2 > fADClookup[0])
    {
      AliDebug(8, Form("Overflow Index ADC = %d Max Index = %d Value = %d. Done with ADCs in this MCM. Is this already MCM header 0x%x?", 
 		       fADClookup[1] - 2, fADClookup[0], fADClookup[fADClookup[1]], *fDataWord));
      fTB = 0;    
      fTbSwitchCtr = 0;
      fMCMWordsExpected = 0;
      AliWarning("Trying to recover. Fall back to DecodeMCM.");
      RewindWord();
      ChangeStatus(kNextMCM);
      return kFALSE;
    }

  if ( (*fDataWord & 0x00000003) != 0x2 && (*fDataWord & 0x00000003) != 0x3) 
    {
      AliWarning(Form("Data %08x : Data Word ends neither with b11 nor b10", (Int_t)*fDataWord));
      fRawReader->AddMinorErrorLog(kDataMaskError,Form("Data %08x", (Int_t)*fDataWord));
      fMCMWordsExpected = 0;
      fTB = 0;    
      fTbSwitchCtr = 0;
      AliWarning("Trying to recover. Fall back to DecodeMCM.");
      RewindWord();
      ChangeStatus(kNextMCM);
      return kFALSE;
  }

  if ( (*fDataWord & 0x3) != fLastADCmask || fTbSwitchCtr >= fTimeWords) 
    {    
      fADC = fADClookup[fADClookup[1]];
      AliDebug(8, Form("Next fADC = %d at index = %d MCM Word Number: %d Max MCM Words is %d", 
		       fADC, fADClookup[1] - 2, fMCMWordCrt, fMCMWordsExpected));
      AliDebug(10, Form("LastMask 0x%x this data 0x%x OR fTbSwitchCtr %d  >= fTimeWords %d", 
			fLastADCmask, (*fDataWord) & 0x3, fTbSwitchCtr, fTimeWords));
      if ((*fDataWord & 0x3) != fLastADCmask && fTbSwitchCtr < fTimeWords && fTbSwitchCtr > 0)
	{
	  AliWarning(Form("Change of mask! But not all words read!"));
	  AliWarning(Form("Next fADC = %d at index = %d MCM Word Number: %d Max MCM Words is %d", 
			     fADC, fADClookup[1] - 2, fMCMWordCrt, fMCMWordsExpected));
	  AliWarning(Form("LastMask 0x%x this data 0x%x OR fTbSwitchCtr %d  >= fTimeWords %d", 
			  fLastADCmask, (*fDataWord) & 0x3, fTbSwitchCtr, fTimeWords));
	}
      ++fADClookup[1];
      fTB = 0;    
      fTbSwitchCtr = 0;
      ChangeStatus(kNextData);
      fLastADCmask = (*fDataWord) & 0x3;
    }

  ++fTbSwitchCtr;

  //decode data here
  Bool_t kIsDataOK = kFALSE;

  // We have only 21 ADC channels.
  if ( fADC > fMaxADCgeom - 1 || fADC < 0) 
    {
      AliWarning(Form("Data 0x%08x : Data is strange. fADC is %d", (Int_t)*fDataWord, (Int_t)fADC));
      AliWarning(Form("fADClookup[0] %d fADClookup[1] %d", fADClookup[0], fADClookup[1]));
      fRawReader->AddMinorErrorLog(kADCChannelOverflow,Form("Data %08x : fADC=%d", (Int_t)*fDataWord, (Int_t)fADC));
    }
  else
    {
      // There are 18 pads connected to each MCM ADC channels 2...19. The other channels cross to other
      // MCMs and are good for online tracking in the MCM.
      if (fSharedPadsOn == kTRUE)
	{
	  kIsDataOK = DecodeADCWord();
	  if (fADC <= 1 || fADC == fMaxADCgeom - 1)
	    {
	      AliDebug(9, Form("Shared Pad fADC == %d", fADC));
	      fIsPadShared = kTRUE;
	    }
	  else
	    {
	      fIsPadShared = kFALSE;	      
	    }
	}
      else
	{
	  if ( fADC > 1 && fADC < fMaxADCgeom - 1 ) 
	    {
	      kIsDataOK = DecodeADCWord();
	    }
	  else 
	    {      
	      AliDebug(9, Form("fADC not accepted - shared pad %d - DATA : 0x%x", fADC, *fDataWord));
	      fCOL = -1;
	      kIsDataOK = kFALSE;
	    }
	}
    }// if fADC is ok

  ++fMCMWordCrt;
  AliDebug(8, Form("We expect %d MCM words. We read %d so far.Current word is 0x%x", fMCMWordsExpected, fMCMWordCrt, *fDataWord));

  // all mcm data processed go to next one
  if ( fMCMWordCrt >= fMCMWordsExpected)
    {
      AliDebug(8, Form("We expect %d MCM words. We read %d so far. Going to next MCM.", fMCMWordsExpected, fMCMWordCrt));
      ChangeStatus(kNextMCM);
    }
  
  return kIsDataOK;
}

//____________________________________________________________________________
Bool_t AliTRDRawStreamTB::DecodeMCM()
{
  fTbSwitch    = 3;  // For first adc channel we expect: (*fDataWord & 3) = 3
  fTbSwitchCtr = 0;  // 
  fADC = fTB   = 0;  // Reset Counter
  fLastADCmask = 0; // Reset

  //AliDebug(3, "Here");

  //if( ((*fDataWord & 0x80000000) == 0x0) && ((*fDataWord & 0x0000000f) == 0xC) )
  //if( ((*fDataWord & 0xf0000000) == 0x80000000) && ((*fDataWord & 0x0000000f) == 0xC) )
  if( (*fDataWord & 0x0000000f) == 0xC ) 
    { // MCM Header
      // this changed the status
      DecodeMCMheader();
      if ( fMCM < 0 || fMCM > 15 || fROB < 0 || fROB > 7 ) 
	{
	  AliWarning("Wrong fMCM or fROB. Skip this data");
	  fRawReader->AddMajorErrorLog(kWrongMCMorROB,Form("MCM=%d, ROB=%d",fMCM,fROB));
	  ChangeStatus(kNextHC);
	  return kFALSE;
	}

      if (fMCMWordCrt < fMCMWordsExpected)
 	{
 	  AliDebug(5, Form("Going to read data."));
 	  ChangeStatus(kNextData);
 	}
      else
 	{
 	  //AliDebug(5, Form("! fMCMWordCrt < fMCMWordsExpected"));
 	  ChangeStatus(kNextMCM);
 	}
      //fEndOfDataFlag = kFALSE;
      return kTRUE;
    }

  if ( *fDataWord == kEndOfRawDataMarker ) 
    {  // End of half-chamber data, finished
      AliDebug(10, "We have reached the eof data. Next should be HC header.");
      fGTUctr1 = -1;

      ChangeStatus(kNextHC);
      fEndOfDataFlag = kTRUE;
      //AliDebug(5, "Expecting MCM header but got End-Of-Raw-Data Marker");
      if (fMCMWordsExpected == 0 || fMCMWordsExpected == fMCMWordCrt)
	{
	  AliDebug(5, Form("Got all mcm words. Returning true."));
	  return kTRUE;
	}
      else
	{
	  AliDebug(5, Form("MCM words missing? %d [expected=%d got=%d] ", fMCMWordsExpected - fMCMWordCrt, fMCMWordsExpected, fMCMWordCrt));	  
	  //AliWarning(Form("MCM words missing? %d [expected=%d got=%d] ", fMCMWordsExpected - fMCMWordCrt, fMCMWordsExpected, fMCMWordCrt));	  
	  return kFALSE;
	}
    }

  //AliDebug(3, Form("Expecting MCM header but got 0x%x. Going to Next MCM header.", *fDataWord));
  AliWarning(Form("Expecting MCM header but got 0x%x. Fall back: Next MCM header.", *fDataWord));
  ChangeStatus(kNextMCM);      

  return kFALSE;
}

//____________________________________________________________________________
Bool_t AliTRDRawStreamTB::DecodeHC()
{
  //AliDebug(3, "Here");

  // left overs from the last chamber?
  if (*fDataWord == kEndOfRawDataMarker)
    {
      AliDebug(3, "We have reached the eof data. Next should be HC header or end of the event!");
      ChangeStatus(kNextHC);
      return kFALSE;
    }

//   if (kFALSE == fSMinfo[fiSMx].fTrackletEnable)
//     {
//       // jump to data decoding - no traklets to expect
//       ChangeStatus(kDecodeHC);
//     }
//   else
//     {
//       ChangeStatus(kNextHC);
//     }

  if ( fNextStatus == kNextHC && fSMinfo[fiSMx].fTrackletEnable == kTRUE)
    {
      //AliDebug(5, "kNextHC");
      //
      // 1) Find end_of_tracklet_marker if tracklets present!
      //
      // GTU Link Mask?

      // endoftrackletmarker?
      if ( *fDataWord == kEndOfTrackletMarker ) 
	{
	  AliDebug(3, "End-of-tracklet-marker found");
	  AliDebug(5, Form("Data 0x%x", *fDataWord));
	  ChangeStatus(kSeekNonEoTracklet);
	  return kTRUE;
	} 
      else 
	{
	  // Tracklets found
	  //AliDebug(3, "Tracklet found");
	  AliDebug(5, Form("Tracklet data 0x%x", *fDataWord));
	  DecodeTracklet();
	  return kTRUE;
	}
    } // if next HC
  else
    {
      if (fSMinfo[fiSMx].fTrackletEnable == kFALSE)
	ChangeStatus(kDecodeHC);
    }

  if (fNextStatus == kSeekNonEoTracklet)
    {
      AliDebug(5, "kSeekNonEoTracklet");

      //
      // 2) Look for non-end_of_tracklet_marker
      //
      //printf("Word %d: 0x%08x\n", fWordCtr, *fDataWord); 
      
      if ( *fDataWord != kEndOfTrackletMarker ) 
	{
	  ChangeStatus(kDecodeHC);
	  AliDebug(3, "NON end-of-tracklet-marker found");
	  AliDebug(5, Form("Data 0x%x", *fDataWord));
	  //// no do not continue - this should be the hcheader
	}
      else
	{
	  //just go on and find the non-end_of_tracklet_marker
	  return kTRUE;
	}
    }

  if ( fNextStatus == kDecodeHC )
    {
      AliDebug(5, "kDecodeHC");
      
      //
      // 3) This Word must be Half Chamber Header
      //
      if ( (*fDataWord & 0xf0000000) == 0x80000000 && (*fDataWord & 0x00000003) == 1 ) 
	{ // HC header
	  AliDebug(5, Form("Is this the HC header? 0x%x", *fDataWord));
	  DecodeHCheader(fTimeBinsCalib); // This is the new header!
// 	  fLastDET = fDET;
	  fDET    = fGeo->GetDetector(fLAYER, fSTACK, fSM);
	  fRowMax = fGeo->GetRowMax(fLAYER,fSTACK,fSM);
	  fColMax = fGeo->GetColMax(fROC);
	  
	  fMCMHctr2 = 0;
	  fHCdataCtr = 0;
	  fChamberDone[fDET]++;
	  AliDebug(2, Form("--------------      DET %d fChamberDone[fDET]=%d", fDET, fChamberDone[fDET]));

	  ChangeStatus(kNextMCM);
	  return kTRUE;
	} //HC header
      else
	{
	  AliWarning(Form("Expecting HC header mask but got 0x%x. Fall back: Next HC.", *fDataWord));
	  ChangeStatus(kNextHC);
	  // before we went to //ChangeStatus(kNextSM);
	}
    } // if decode HC
  
  return kFALSE;
}
//____________________________________________________________________________
Bool_t AliTRDRawStreamTB::DecodeGTUlinkMask()
{
  //
  // Decode the link masks sent by the GTU. These marke the active optical links
  // between GTU and Super Module. Up to now only fully active links are found
  // (0xfff = 12 active links).
  //

  if ( (*fDataWord & 0xfffff000) ==  0xe0000000 )
    {
      if ( fRawVersion < 1 || fRawVersion > 3 ) 
	{
	  AliError(Form(" Unsupported raw version: %d", fRawVersion));      
	}
      
      if ( fGTUctr1 == -1 ) fGTUctr2++;
      fGTUctr1++;

      if ( (fGTUctr1 >= 0) && (fGTUctr1 < 5) && (fGTUctr2 >= 0) && (fGTUctr2 < 18) ) 
	{
	  fGTUlinkMask[fGTUctr2][fGTUctr1] = (*fDataWord & 0xfff);
	}
     
      //AliDebug(5, Form("GTU link mask 0x%x decoded 0x%x", *fDataWord, fGTUlinkMask[fGTUctr2][fGTUctr1]));
      return kTRUE;
    }

  return kFALSE;
}

// //____________________________________________________________________________
void AliTRDRawStreamTB::ChangeStatus(Int_t kstat)
{
  fLastStatus = fNextStatus;
  fNextStatus = kstat;  
}

//____________________________________________________________________________
Bool_t AliTRDRawStreamTB::DecodeSM()
{
  fLastDET  = -1;
  fDET      = -1;
  fRetVal   = 0;
  fEqID     = 0;
  fDataSize = 0;
  fSizeOK = kFALSE;

  Int_t status = DecodeHeadingInfo();
  if (status == kWordOK)
    {    
      ChangeStatus(kNextHC);
      return kTRUE;
    } 
  else 
    {
      AliWarning(Form("Decoding SM info failed. Fall back: None. Stop.", fEqID));
      fRawReader->AddMajorErrorLog(kGTULinkMaskMissing,Form("Equipment %d",fEqID));
      ChangeStatus(kStop);
    }	    

  return kFALSE;
}

//____________________________________________________________________________
Bool_t AliTRDRawStreamTB::Next()
{
  //
  // Updates the next data word pointer
  //

  if (fNextStatus == kStart)
    {
      Init();
    }

  while (fNextStatus != kStop && fNextStatus != kNoMoreData)
    { // !kStop
      NextData();
      
      switch (fNextStatus)
	{
	case kNextData:
	  {
	    if (DecodeNextRawWord() == kTRUE)
	      {
		fTB += 3;
		if (fTB - 3 > 29)
		  {
		    AliWarning(Form("Invalid time bin %d. sm %d det %d rob %d col %d row %d mcm=%d adc=%d ", fTB-3, fSM, fDET, fROB, fCOL, fROW, fMCM, fADC));
		    AliWarning(Form("max=%d index=%d adc=%d", fADClookup[0], fADClookup[1], fADClookup[fADClookup[1]-1]));
		  }
		if (fSig[0] > fRawDigitThreshold || fSig[1] > fRawDigitThreshold || fSig[2] > fRawDigitThreshold) 
		  return kTRUE;
	      }
	  }; break;
	case kNextMCM:
	  {
	    if (DecodeMCM() == kFALSE)
	      AliWarning(Form("Decode MCM unsuccessfull. Current Word 0x%x at pos 0x%x", *fDataWord, fPos));	  
	  }; break;
	case kNextHC:
	case kSeekNonEoTracklet:
	case kDecodeHC:
	  {
	    if (DecodeHC() == kFALSE)
	      {
		if (*fDataWord == kEndOfRawDataMarker)
		  {
		    AliDebug(2, Form("End of data at 0x%x", fPos));
		  }
		else
		  {
		    AliWarning(Form("Decode HC unsuccessfull. Current Word 0x%x at pos 0x%x", *fDataWord, fPos));	  
		  }
	      }
	  }; break;
	case kNextSM:
	  {
	    if (DecodeSM() == kFALSE)
		AliWarning(Form("Decode SM unsuccessfull. Current Word 0x%x at pos 0x%x", *fDataWord, fPos));	  
	  }; break;
	case kNoMoreData:
	case kStop:
	  ; break;
	default:
	  AliWarning(Form("Unknown state %d. Last state %d. Current Word 0x%x at pos 0x%x", fNextStatus, fLastStatus, *fDataWord, fPos));  
	  ChangeStatus(kStop);
	};

    } // not kStop

  //AliDebug(1, Form("That's all folks! %d", fSM));
  return kFALSE;
}

//____________________________________________________________________________
Int_t AliTRDRawStreamTB::NextChamber(AliTRDdigitsManager *digitsManager)
{
  //
  // Fills single chamber digit array 
  // Return value is the detector number
  //

  AliTRDdataArrayI *digits = 0;
  AliTRDdataArrayI *track0 = 0;
  AliTRDdataArrayI *track1 = 0;
  AliTRDdataArrayI *track2 = 0; 
  AliTRDSignalIndex *indexes = 0;

  // Loop through the digits
  Int_t lastdet = -1;
  Int_t det    = -1;
//   Int_t returnDet = -1;
  Int_t it = 0;
  while (Next()) 
    {      
      det    = fDET;
    
      if (det != lastdet) 
	{ 
	  // If new detector found
	  if (lastdet == -1)
	    {
	      lastdet = det;
	    }
	  else
	    {
	      return lastdet;
	    }
	  
	  // Add a container for the digits of this detector
	  digits = digitsManager->GetDigits(det);
	  track0 = digitsManager->GetDictionary(det,0);
	  track1 = digitsManager->GetDictionary(det,1);
	  track2 = digitsManager->GetDictionary(det,2);

	  // Allocate memory space for the digits buffer
	  if (digits->GetNtime() == 0) 
	    {
	      digits->Allocate(fRowMax,fColMax, fTBins);
	      track0->Allocate(fRowMax,fColMax, fTBins);
	      track1->Allocate(fRowMax,fColMax, fTBins);
	      track2->Allocate(fRowMax,fColMax, fTBins);
	    }

	  indexes = digitsManager->GetIndexes(det);
	  indexes->SetSM(GetSM());
	  indexes->SetStack(GetStack());
	  indexes->SetLayer(GetLayer());
	  indexes->SetDetNumber(det);
	  if (indexes->IsAllocated() == kFALSE)
	    indexes->Allocate(fRowMax, fColMax, fTBins);
	}
    
      // 3 timebin data are stored per word
      for (it = 0; it < 3; it++)
	{
// 	  if ( GetTimeBin() + it < fTBins )
// 	    {
	      if (fSig[it] > 0)
		{
		  digits->SetDataUnchecked(fROW, fCOL, fTB - 3 + it, fSig[it]);

		  indexes->AddIndexTBin(fROW, fCOL, fTB - 3 + it);
		  track0->SetDataUnchecked(fROW, fCOL, fTB - 3 + it, 0);
		  track1->SetDataUnchecked(fROW, fCOL, fTB - 3 + it, 0);
		  track2->SetDataUnchecked(fROW, fCOL, fTB - 3 + it, 0);
		}
// 	    }
	} // tbins
    }// while Next()

  // what happens if the last HC is turned off?
  return det;
  //return -1;
}

// new stuff
//____________________________________________________________________________
Int_t AliTRDRawStreamTB::SkipWords(UInt_t iw)
{
  Int_t status = kWordOK;
  for (UInt_t i = 0; i < iw; i++)
    {
      status = NextData();
      AliDebug(5, Form("Skipping word %d of %d [0x%x]", i+1, iw, *fDataWord));
      if (status != kWordOK)
	return status;
    }

  //status = NextData();
  return status;
}

//____________________________________________________________________________
Int_t AliTRDRawStreamTB::DecodeHeadingInfo()
{
  // SM info...
  fSMinfo[fiSMx].Decode(fDataWord);
  if (fgExtraDebug)
    fSMinfo[fiSMx].Dump();

  Int_t status = SkipWords(fSMinfo[fiSMx].fHeaderSize);
  if (status != kWordOK)
    return status;
  else
    status = NextData();

  //fSMinfo[fiSMx].fStackActive[0] = kTRUE;
  // now stack info..
  if (status == kWordOK)
    {
      for (Int_t i = 0; i < 5; i++)
	{
	  if (fSMinfo[fiSMx].fStackActive[i] == kFALSE)
	    continue;
	  AliDebug(5, Form("Decode stack info %d 0x%x", i, *fDataWord));
	  if (fgStackIndexBug == kFALSE)
	    fStackInfo[fiSMx][i].Decode(fDataWord);
	  else
	    fStackInfo[fiSMx][i].DecodeBug(fDataWord);
	  
	  if (fgExtraDebug)
	    fStackInfo[fiSMx][i].Dump();
	  
	  status = SkipWords(fStackInfo[fiSMx][i].fHeaderSize);
	  if (status != kWordOK)
	    return status;      
	  else
	    if (i < 4)
	      status = NextData();

	  if (status != kWordOK)
	    return status;      
	}
    }

  RewindWord();
  return status;
}
