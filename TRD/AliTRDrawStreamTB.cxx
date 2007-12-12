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
// Author: M. Ploskon (ploskon@ikf.uni-frankfurt.de)                         //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "TString.h"
#include "TFile.h"
#include "TTreeStream.h"

#include "AliTRDrawStreamTB.h"
#include "AliTRDgeometry.h"
#include "AliTRDfeeParam.h"
#include "AliTRDdigitsManager.h"
#include "AliTRDdataArrayI.h"
#include "AliTRDSignalIndex.h"

#include "AliLog.h"
#include "AliRawReader.h"

#define END_OF_TRACKLET_MARKEROLD 0xaaaaaaaa
#define END_OF_TRACKLET_MARKERNEW 0x10001000
#define ENDOFRAWDATAMARKER 0x00000000
#define WORD_SIZE sizeof(UInt_t)           // size of a word in bytes
#define EXTRA_LEAD_WORDS 24
#define CDH_WORDS 8

#define IS_BIT_SET(w,b) ( ((w) >> (b)) & 0x1 ) // 1 if bit b is set in word w
#define GET_VALUE_AT(w,m,s) (( (w) >> (s)) & (m) )      // get value of word w rshifted by s and mask with m

// SM index word masks:
#define SM_HEADER_SIZE(w) GET_VALUE_AT(w,0xffff,16) 
#define TRACKLETS_ENABLED(w) IS_BIT_SET(w,5)
#define STACK_MASK(w) ((w) & 0x1f)

// Stack word masks
#define STACK_HEADER_SIZE(w) GET_VALUE_AT(w,0xffff,16)
#define STACK_LINK_WORD(w) ((w) & 0xfff)

// HC word masks
//#define HC_HEADER_MASK_ERR(w) ( ((w) & (0x80000003)) == (0x80000001) ? 0 : 1) // 0 if OK!!!
#define HC_HEADER_MASK_ERR(w) ( ((w) & (0x3)) == (0x1) ? 0 : 1) // 0 if OK!!!

// HC word 0
#define HC_SPECIAL_RAW_VERSION(w) IS_BIT_SET(w,31)
#define HC_MAJOR_RAW_VERSION(w) GET_VALUE_AT(w,0x7f,24)
#define HC_MINOR_RAW_VERSION(w) GET_VALUE_AT(w,0x7f,17)
#define HC_EXTRA_WORDS(w) GET_VALUE_AT(w,0x7,14)
#define HC_DCS_BOARD(w) GET_VALUE_AT(w,0xfff<<20,20)
#define HC_SM_NUMBER(w) GET_VALUE_AT(w,0x1f,9)
#define HC_LAYER_NUMBER(w) GET_VALUE_AT(w,0x7,6)
#define HC_STACK_NUMBER(w) GET_VALUE_AT(w,0x7,3)
#define HC_SIDE_NUMBER(w) IS_BIT_SET(w,2)

// HC word 1
#define HC_NTIMEBINS(w) GET_VALUE_AT(w,0x3f,26)
#define HC_BUNCH_CROSS_COUNTER(w) GET_VALUE_AT(w,0xffff,10)
#define HC_PRETRIGGER_COUNTER(w) GET_VALUE_AT(w,0xf,6)
#define HC_PRETRIGGER_PHASE(w) GET_VALUE_AT(w,0xf,6)

//--------------------------------------------------------

// MCM word and ADC mask
#define MCM_HEADER_MASK_ERR(w) ( ((w) & (0xf)) == (0xc) ? 0 : 1) // 0 if OK!!!
#define MCM_ADCMASK_MASK_ERR(w) ( ((w) & (0xf)) == (0xc) ? 0 : 1) // 0 if OK!!!
#define MCM_MCM_NUMBER(w) GET_VALUE_AT(w,0x0f,24)
#define MCM_ROB_NUMBER(w) GET_VALUE_AT(w,0x7,28)
#define MCM_EVENT_COUNTER(w) GET_VALUE_AT(w,0x00fffff,4)
#define MCM_ADCMASK_VAL(w) GET_VALUE_AT(w,0x1fffff,4)

//this should give 0x1fffff as a mask
#define MCM_DUMMY_ADCMASK_VAL 0x001fffffc 
//--------------------------------------------------------

#define MAX_TRACKLETS_PERHC 256 // max number of tracklets per HC - large number for now

//--------------------------------------------------------
#define ADC_WORD_MASK(w) ((w) & 0x3)
//--------------------------------------------------------
ClassImp(AliTRDrawStreamTB)

Bool_t AliTRDrawStreamTB::fgExtraSkip = kFALSE;
Bool_t AliTRDrawStreamTB::fgSkipCDH = kFALSE;
Bool_t AliTRDrawStreamTB::fgWarnError = kTRUE;
Bool_t AliTRDrawStreamTB::fgCleanDataOnly = kTRUE;
Bool_t AliTRDrawStreamTB::fgDebugFlag = kTRUE;
Bool_t AliTRDrawStreamTB::fgDebugStreamFlag = kFALSE;
TTreeSRedirector *AliTRDrawStreamTB::fgDebugStreamer = 0;
UInt_t AliTRDrawStreamTB::fgStreamEventCounter = 0;
UInt_t AliTRDrawStreamTB::fgDumpHead = 0;

Int_t AliTRDrawStreamTB::fgEmptySignals[] = 
  {
    -1, -1, -1, -1, -1,  -1, -1, -1, -1, -1,  -1, -1, -1, -1, -1
    -1, -1, -1, -1, -1,  -1, -1, -1, -1, -1,  -1, -1, -1, -1, -1
  };

AliTRDrawStreamTB::AliTRDrawStreamTB()
  : TObject()
  , fSM()
  , fStack(0)
  , fHC(0)
  , fMCM(0)
  , fADC(0)
  , fpPos(0)
  , fpBegin(0)
  , fpEnd(0)
  , fWordLength(0)
  , fStackNumber(-1)
  , fStackLinkNumber(-1)
  , fhcMCMcounter(0)
  , fmcmADCcounter(0)
  , fLinkTrackletCounter(-1)
  , fEndOfTrackletCount(-1)
  , fMaskADCword(0)
  , fTbinADC(0)
  , fDecodedADCs(-1)
  , fEventCounter(0)
  , fLastEventCounter(0)
  , fCommonAdditive(0)
  , fSharedPadsOn(kFALSE)
  , fMaxADCgeom(0)
  , fGeometry(0)
  , fRawReader(0)
  , fTRDfeeParam(0)
  , fDebugStreamOwned(kFALSE)
{
  //
  // default constructor
  //

  if (Init() == kFALSE)
    {
      AliWarning("Unable to Init.");	  
    }
}

//--------------------------------------------------------
AliTRDrawStreamTB::AliTRDrawStreamTB(AliRawReader *rawReader)
  : TObject()
  , fSM()
  , fStack(0)
  , fHC(0)
  , fMCM(0)
  , fADC(0)
  , fpPos(0)
  , fpBegin(0)
  , fpEnd(0)
  , fWordLength(0)
  , fStackNumber(-1)
  , fStackLinkNumber(-1)
  , fhcMCMcounter(0)
  , fmcmADCcounter(0)
  , fLinkTrackletCounter(-1)
  , fEndOfTrackletCount(-1)
  , fMaskADCword(0)
  , fTbinADC(0)
  , fDecodedADCs(-1)
  , fEventCounter(0)
  , fLastEventCounter(0)
  , fCommonAdditive(0)
  , fSharedPadsOn(kFALSE)
  , fMaxADCgeom(0)
  , fGeometry(0)
  , fRawReader(rawReader)
  , fTRDfeeParam(0)
  , fDebugStreamOwned(kFALSE)
{
  //
  // default constructor
  //
  if (fRawReader)
    {	  
      if (Init() == kFALSE)
	{
	  AliWarning("Unable to Init. Try setting up the reader with SetReader or buffer with Init(void *, UInt_t )");	  
	}
    }
  else
    {
      AliWarning("Unable to setup reader. Use SetReader(AliRawReader*).");
    }
}

//------------------------------------------------------------

AliTRDrawStreamTB::AliTRDrawStreamTB(const AliTRDrawStreamTB& /*st*/)
  : TObject()
  , fSM()
  , fStack(0)
  , fHC(0)
  , fMCM(0)
  , fADC(0)
  , fpPos(0)
  , fpBegin(0)
  , fpEnd(0)
  , fWordLength(0)
  , fStackNumber(-1)
  , fStackLinkNumber(-1)
  , fhcMCMcounter(0)
  , fmcmADCcounter(0)
  , fLinkTrackletCounter(-1)
  , fEndOfTrackletCount(-1)
  , fMaskADCword(0)
  , fTbinADC(0)
  , fDecodedADCs(-1)
  , fEventCounter(0)
  , fLastEventCounter(0)
  , fCommonAdditive(0)
  , fSharedPadsOn(kFALSE)
  , fMaxADCgeom(0)
  , fGeometry(0)
  , fRawReader(0)
  , fTRDfeeParam(0)
  , fDebugStreamOwned(kFALSE)
{
  //
  // Copy constructor
  // 
  AliError("Not implemeneted.");
}

//------------------------------------------------------------
void AliTRDrawStreamTB::SetRawVersion(Int_t fraw)
{
  //
  // function provided for backward compatibility
  //
  AliWarning("Raw data version is read from raw data stream! No point of setting it in here.");
  fraw = 0; // avoid warnings
}

//------------------------------------------------------------
AliTRDrawStreamTB::~AliTRDrawStreamTB()
{
  //
  // destructor
  //
  delete fGeometry;

  if (fDebugStreamOwned == kTRUE)
    {
      if (fgDebugStreamer)
	{
	  delete fgDebugStreamer;
	  fgDebugStreamer = 0;
	}
    }
}

//------------------------------------------------------------

AliTRDrawStreamTB &
AliTRDrawStreamTB::operator=(const AliTRDrawStreamTB &)
{
  //
  // we are not using this functionality
  //
  AliFatal("May not use.");
  return *this;
}

//------------------------------------------------------------

void    
AliTRDrawStreamTB::EnableDebug(TTreeSRedirector *debugStream)
{
  //
  // replace the current debug streamer or create a new one if owned
  //

  if (debugStream == 0)
    {
      if (fDebugStreamOwned == kTRUE)
	{
	  if (fgDebugStreamer)
	    delete fgDebugStreamer;
	}

      fDebugStreamOwned = kTRUE;
      fgDebugStreamer    = new TTreeSRedirector("TRDrawDataDebug.root");
      fgStreamEventCounter = 0;

    }
  else
    {
      if (fDebugStreamOwned == kTRUE)
	{
	  if (fgDebugStreamer)
	    delete fgDebugStreamer;
	}

      fgDebugStreamer = debugStream;
      fDebugStreamOwned = kFALSE;
    }
}

//___________________________________________________________
void 
AliTRDrawStreamTB::SwapOnEndian()
{
  //
  // Check the endian and swap if needed
  //

  int itemp = 1;
  char* ptemp = (char*) &itemp;
  if (ptemp[0] != 1)
    {
      // need to swap...
      if (fgDebugFlag) AliDebug(8, "Swapping.");

      fpPos = fpBegin;
      UInt_t iutmp = 0;
      while (fpPos < fpEnd)
	{
	  fpPos += 1;
	  iutmp = (((*fpPos & 0x000000ffU) << 24) | ((*fpPos & 0x0000ff00U) <<  8) |
		   ((*fpPos & 0x00ff0000U) >>  8) | ((*fpPos & 0xff000000U) >> 24));
	  // here we override the value in the buffer!
	  *fpPos = iutmp;	  
	}
      fpPos = fpBegin;
    }
}

//------------------------------------------------------------
void 
AliTRDrawStreamTB::DumpErrorCount()
{
  //
  // print the error statistics into the stdout
  //

  ;
  //   AliInfo(Form("Error counts: HC0 %d HC1 %d MCM %d ADCmask %d ADC %d", 
  // 	       fHC->fH0ErrorCounter,  fHC->fH1ErrorCounter, fMCM->fErrorCounter, fMCM->fMaskErrorCounter, adc.fErrorCounter));
}

//------------------------------------------------------------
Bool_t 
AliTRDrawStreamTB::DumpWords(UInt_t *px, UInt_t iw, UInt_t marker)
{
  //
  // skip few words
  // note: can be made faster
  //
    
  TString tsreturn = Form("\n[ Dump Sequence at 0x%08x ] : ", px);
  for (UInt_t i = 0; i < iw; i++)
    {
      if (px + iw >= fpEnd)
	return kFALSE;

      if (i % 8 == 0)
	tsreturn += "\n                              ";
      if (marker != 0 && marker == px[i])
	tsreturn += Form(" *>0x%08x<* ", px[i]);
      else
	tsreturn += Form("0x%08x ", px[i]);
    }
  tsreturn += "\n";

  AliInfo(tsreturn.Data());

  return kTRUE;
}

//------------------------------------------------------------
Bool_t 
AliTRDrawStreamTB::SkipWords(UInt_t iw)
{
  //
  // skip few words
  // note: can be made faster
  //
  
  if ( fpPos + iw < fpEnd )
    {
      fpPos += iw;
      return kTRUE;
    }
  else
    {
      if (fgWarnError) AliWarning(Form("Skip %d words failed. %d available", iw, fpEnd - fpPos - 1));
      if (fgDebugStreamer)
	{
	  TTreeSRedirector &cstream = *fgDebugStreamer;
	  cstream << "SkipWords"
		  << "Event=" << fgStreamEventCounter
		  << ".nwords=" << iw
// 		  << ".read=" << (Char_t*)(fpPos - fpBegin)
// 		  << ".togo=" << (Char_t*)(fpEnd - fpPos)
// 		  << ".length=" << (Char_t*)(fpEnd - fpBegin)
		  << "\n";
	}
      return kFALSE;
    }

  return kTRUE;
}

//------------------------------------------------------------
Bool_t 
AliTRDrawStreamTB::SetReader(AliRawReader *reader)
{
  //
  //
  //

  if (reader != 0)
    {
      fRawReader = reader;
      if (fRawReader)
	{	  
	  return Init();
	}
      else
	{
	  AliWarning("Unable to setup reader.");
	  return kFALSE;
	}
    }
  else
    {
      AliWarning("AliRawReader argument is 0.");
      fRawReader = 0;
    }

  return kFALSE;
}

//------------------------------------------------------------
Int_t 
AliTRDrawStreamTB::NextBuffer()
{
  //
  // return -1 if no more buffers available
  // return 0 DecodeSM failed (clean data required for example) but still maybe more data to come
  // return 1 DecodeSM OK
  // 

  if (fRawReader != 0)
    {
      UChar_t *buffer = 0;
      UInt_t length = 0;
      Bool_t kBufferSet = fRawReader->ReadNextData(buffer);
      if (kBufferSet == kTRUE)
	{
	  if (fgDebugFlag)  AliDebug(9, "Buffer is set.");
	  length = fRawReader->GetDataSize();
	  if (fgExtraSkip == kTRUE)
	    {
	      buffer += EXTRA_LEAD_WORDS * WORD_SIZE;
	      length -= EXTRA_LEAD_WORDS * WORD_SIZE;
	    }

	  if (fgSkipCDH == kTRUE)
	    {
	      buffer += CDH_WORDS * WORD_SIZE;
	      length -= CDH_WORDS * WORD_SIZE;	      
	    }

	  if (length > 0)
	    {
	      if (fgDebugFlag)  AliDebug(9, Form("Buffer length : %d", length));
	      //return InitBuffer((void*)buffer, length);
	      if (DecodeSM((void*)buffer, length) == kTRUE)
		return 1;
	      else
		return 0;
	    }
	}
      else
	{
	  return -1;
	}
    }

  return -1;
}

//------------------------------------------------------------
void 
AliTRDrawStreamTB::ResetCounters()
{
  //
  // reset some global counters
  //

  fStackNumber = 0;
  fStackLinkNumber = 0;

  fDecodedADCs = 0;

  fSM.fActiveStacks = 0;
  fSM.fNexpectedHalfChambers = 0;

  fLinkTrackletCounter = 0;

  fLastEventCounter = 0;
  fEventCounter = 0;
}

//------------------------------------------------------------
void 
AliTRDrawStreamTB::ResetIterators()
{
  //
  // reset the data iterators used in the Next()
  //
  fStackNumber = 0;
  fStackLinkNumber = 0;
  fhcMCMcounter = 0;  
  fmcmADCcounter = 0;
}

//------------------------------------------------------------

Bool_t 
AliTRDrawStreamTB::Next()
{
  //
  // returns with true on next adc read
  // returns false on errors and end of buffer
  // 

  while (fStackNumber < 5 && fSM.fActiveStacks > 0)
    {
      if(fSM.fStackActive[fStackNumber] == kTRUE)
	{
	  fStack = &fSM.fStacks[fStackNumber];
	  while (fStackLinkNumber < 12)
	    {
	      if (fStack->fLinksActive[fStackLinkNumber] == kTRUE)
		{
		  fHC = &fStack->fHalfChambers[fStackLinkNumber];
		  if (!fHC)
		    {
		      AliError(Form("Super Strange. HC missing at stack %d link %d", fStackNumber, fStackLinkNumber));
		      return kFALSE;
		    }
		  //if (fgDebugFlag)  AliDebug(2, DumpHCinfoH0(fHC));		  
		  if (fHC->fCorrupted == 0)
		    {
		      while (fhcMCMcounter < fHC->fMCMmax)
			{
			  fMCM = &fHC->fMCMs[fhcMCMcounter];
			  //if (fgDebugFlag)  AliDebug(2, DumpMCMinfo(fMCM));		  			  
			  if (!fMCM)
			    {
			      AliError(Form("Super Strange. HC missing at stack %d link %d atMCMslot %d", 
					    fStackNumber, fStackLinkNumber, fhcMCMcounter));
			      return kFALSE;
			    }
			  while(fmcmADCcounter < fMCM->fADCmax)
			    {
			      //printf("%d %d %d %d\n", fStackNumber, fStackLinkNumber, fhcMCMcounter, fmcmADCcounter);
			      fADC = &fMCM->fADCs[fmcmADCcounter];
			      if (!fADC)
				{
				  AliError(Form("Super Strange. ADC missing at stack %d link %d MCMslot %d ADCslot %d", 
						fStackNumber, fStackLinkNumber, fhcMCMcounter, fmcmADCcounter));
				  return kFALSE;
				}
			      fmcmADCcounter++;
			      //printf("ADC count %d \n",  fmcmADCcounter);
			      if (fSharedPadsOn)
				{
				  return kTRUE;
				}
			      else
				{
				  if (fADC->fIsShared == kFALSE)
				    return kTRUE;
				}
			    } //while adc in MCM
			  //printf("MCM count %d \n", fhcMCMcounter);
			  fhcMCMcounter++;
			  // next MCM should go through all active ADCs
			  fmcmADCcounter = 0;
			} // while mcm
		    } // if HC OK
		}// if link active
	      fStackLinkNumber++;
	      // next stack link (HC) should go through all active MCMs
	      fhcMCMcounter = 0;
	    }// while links
	}// if stack active
      fStackNumber++;
      // next stack should go through all links - start from 0
      fStackLinkNumber = 0;
    }

  // in case rawreader manages the mem buffers
  // lets go for the next buffer
  if (fRawReader)
    {
      Int_t nextBuff = NextBuffer();
      while (nextBuff != -1)
	{
	  if (nextBuff > 0)
	    return Next();	   	  
	  nextBuff = NextBuffer();
	}
    }

  return kFALSE;
}

//------------------------------------------------------------

Int_t 
AliTRDrawStreamTB::NextChamber(AliTRDdigitsManager *digitsManager)
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
      det    = GetDet();
    
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

	  Int_t rowMax = GetRowMax();
	  Int_t colMax = GetColMax();
	  Int_t ntbins = GetNumberOfTimeBins();

	  // Allocate memory space for the digits buffer
	  if (digits->GetNtime() == 0) 
	    {
	      digits->Allocate(rowMax, colMax, ntbins);
	      track0->Allocate(rowMax, colMax, ntbins);
	      track1->Allocate(rowMax, colMax, ntbins);
	      track2->Allocate(rowMax, colMax, ntbins);
	    }

	  indexes = digitsManager->GetIndexes(det);
	  indexes->SetSM(GetSM());
	  indexes->SetStack(GetStack());
	  indexes->SetLayer(GetLayer());
	  indexes->SetDetNumber(det);
	  if (indexes->IsAllocated() == kFALSE)
	    indexes->Allocate(rowMax, colMax, ntbins);
	}
    
      // ntimebins data are ready to read
      for (it = 0; it < GetNumberOfTimeBins(); it++)
	{
	  if (GetSignals()[it] > 0)
	    {
	      digits->SetDataUnchecked(GetRow(), GetCol(), it, GetSignals()[it]);
	      
	      indexes->AddIndexTBin(GetRow(), GetCol(), it);
	      track0->SetDataUnchecked(GetRow(), GetCol(), it, 0);
	      track1->SetDataUnchecked(GetRow(), GetCol(), it, 0);
	      track2->SetDataUnchecked(GetRow(), GetCol(), it, 0);
	    }
	} // tbins
    }// while Next()

  // what happens if the last HC is turned off?
  return det;
  //return -1;
}

//------------------------------------------------------------
Bool_t
AliTRDrawStreamTB::Init()
{
  //
  // general init
  //

  TDirectory *saveDir = gDirectory; 
  
  if (!fGeometry) 
    {
      fGeometry = new AliTRDgeometry();
    }
  
  if (!fGeometry) 
    {
      AliError("Geometry FAILED!");
      return kFALSE;
    }

  fTRDfeeParam = AliTRDfeeParam::Instance();
  if (!fTRDfeeParam)
    {
      AliError("AliTRDfeeParam FAILED!");
      return kFALSE;
    }

  fMaxADCgeom = (Int_t)fGeometry->ADCmax();

  fCommonAdditive = 10;

  ResetCounters();

  // in case rawreader manages the mem buffers
  // lets go for the next buffer
  if (fRawReader)
    {
      Int_t nextBuff = NextBuffer();
      while (nextBuff != -1)
	{
	  if (nextBuff > 0)
	    return kTRUE;	   	  
	  nextBuff = NextBuffer();
	}
    }

  if (fgDebugStreamFlag == kTRUE)
    {
      if (!fgDebugStreamer)
	{
	  fgDebugStreamer    = new TTreeSRedirector("TRDrawDataDebug.root");
	  fgStreamEventCounter = 0;
      
	  if (fgDebugStreamer)
	    {
	      AliInfo(Form("Debug Streamer Initialized! Remember to delete! %s::DeleteDebugStream()", this->IsA()->GetName()));
	      //fDebugStreamOwned = kTRUE;
	    }
	  else
	    {
	      AliError("Unable to init debug stream");
	    }
	}
    }

  saveDir->cd();

  return kTRUE;
}

void AliTRDrawStreamTB::DeleteDebugStream()
{
  // 
  // static helper function
  //

  if (fgDebugStreamer)
    {
      delete fgDebugStreamer;
      fgDebugStreamer = 0;
    }
  
}

//------------------------------------------------------------
Bool_t 
AliTRDrawStreamTB::InitBuffer(void *buffer, UInt_t length)
{
  // 
  // init the class before reading from the buffer
  // and read SM heading info:
  // -- super module header
  // -- stack headers + link status words
  //

  // not in the pre scan mode
  fgStreamEventCounter++;
  
  ResetCounters();

  fpBegin = (UInt_t *)buffer;

  if (WORD_SIZE == 0)
    {
      AliFatal("Strange! Size of UInt_t == 0");
      return kFALSE;
    }

  fWordLength = length/WORD_SIZE;
  fpEnd = fpBegin + fWordLength;
  fpPos = fpBegin;

  if (fpBegin == 0 || length <= 0)
    {
      AliError(Form("This will not work! Pointer to the buffer is 0x%08x of size %d", fpBegin, length));
      return kFALSE;
    }

  SwapOnEndian();

  if (fgDumpHead > 0)
    {
      AliInfo("------------------------------------------------");
      if (DumpWords(fpBegin, fgDumpHead) == kFALSE)
	{
	  AliError("Dump failed. Not enough data.");	  
	}
      AliInfo("------------------------------------------------");
    }

  //decode SM
  DecodeSMInfo(fpPos, &fSM);

  if (fgDebugFlag)  AliDebug(5, DumpSMInfo(&fSM));

  fpPos++;
  if (fpPos < fpEnd)
    {	
      if (SkipWords(fSM.fHeaderSize) == kTRUE)
	{
	  //decode stacks info
	  //for (Int_t istack = 0; istack < fSM.fActiveStacks; istack++)
	  for (Int_t istack = 0; istack < 5; istack++)
	    {
	      if (fSM.fStackActive[istack] == kFALSE)
		continue;

	      fStack = &fSM.fStacks[istack];
	      DecodeStackInfo(fpPos, fStack);
	      fpPos++;

	      fSM.fNexpectedHalfChambers += fStack->fActiveLinks;
	      
	      if (fgDebugFlag)  AliDebug(5, DumpStackInfo(fStack));
	      
	      if (SkipWords(fStack->fHeaderSize) == kFALSE)
		{
		  if (fRawReader) fRawReader->AddMajorErrorLog(kDecodeStackInfo, "Stack head words missing");
		  return kFALSE;
		}
	    }
	}
      else
	{
	  return kFALSE;
	}
    }
  else
    {
      if (fgWarnError) AliWarning("No Stack info present. Strange.");
      if (fRawReader) fRawReader->AddMajorErrorLog(kDecodeStackInfo, "Stack info missing");
      return kFALSE;
    }

  if (fgDebugFlag)  AliDebug(5, Form("Expected half chambers : %d", fSM.fNexpectedHalfChambers));
  
  //fpPos++;
  if (fpPos < fpEnd)
    {
      if (fgDebugFlag)  AliDebug(5, "Init OK.");
    }
  else
    {
      if (fgWarnError) AliWarning("No data just after init. Strange.");
      if (fRawReader) fRawReader->AddMajorErrorLog(kMissingData, "Missing data");
      return kFALSE;
    }

  return kTRUE;
}

//------------------------------------------------------------
Bool_t 
AliTRDrawStreamTB::DecodeSM(void *buffer, UInt_t length)
{
  //
  // decode all sm at once
  // still here for devel and debug
  //
  
  //memset(&fSM, 0, sizeof(fSM));

  ResetIterators();

  fSM.fClean = kTRUE;
  if (InitBuffer(buffer, length) == kTRUE)
    {
      // no we are already set!
      // fpPos++;      
    }
  else
    {
      if (fgWarnError) AliError("--- INIT failed. ---------------");      
      fSM.fClean = kFALSE;
      return kFALSE;
    }

  //decode data here
  //fSM.fActiveStacks
  for (Int_t istack = 0; istack < 5; istack++)
    {
      fStackNumber = istack;
      if (fSM.fStackActive[istack] == kFALSE)
	continue;
      
      fStack = &fSM.fStacks[istack];

      //fStack[istack].fActiveLinks // max is 12
      for (Int_t ilink = 0; ilink < 12; ilink++)
	{
	  fStackLinkNumber = ilink;
	  if (fStack->fLinksActive[ilink] == kFALSE)
	    continue;

	  if (fpPos >= fpEnd)
	    {
	      if (fRawReader) fRawReader->AddMajorErrorLog(kLinkDataMissing, "Link data missing");	      
	      fSM.fClean = kFALSE;
	      break;
	    }

	  fHC = &fStack->fHalfChambers[ilink];

	  if (fSM.fTrackletEnable == kTRUE)
	    {
	      if (DecodeTracklets() == kFALSE)
		{
		  // 	      if (fgWarnError) AliError(Form("--- Decode Tracklets failed. --------------- Link %d of %d", 
		  // 					    ilink + 1, fStack->fActiveLinks));
		  if (fgDebugStreamer)
		    {
		      TTreeSRedirector &cstream = *fgDebugStreamer;
		      cstream << "TrackletDecodeError"
			      << "Event=" << fgStreamEventCounter
			      << ".Stack=" << fStackNumber
			      << ".Link=" << fStackLinkNumber
			      << ".EndOfTrackletCount=" << fEndOfTrackletCount
			      << "\n";
		    }
		  
		  fSM.fClean = kFALSE;
		  SeekEndOfData();

		  if (fgWarnError) 
		    {
		      AliError(Form("Failed stack %d link %d", fStackNumber, fStackLinkNumber));
		      AliError(Form("Debug Event Counter : %d", fgStreamEventCounter)); 
		    }
		  continue;
		  //return kFALSE;
		}
	    }

	  if (fpPos >= fpEnd)
	    {
	      if (fRawReader) fRawReader->AddMajorErrorLog(kHCdataMissing, "HC data missing");	      
	      fSM.fClean = kFALSE;
	      break;
	    }
	  
	  if (DecodeHC() == kFALSE)
	    {
// 	      if (fgWarnError) AliError(Form("--- Decode HC failed. --------------- Link %d of %d", 
// 					    ilink + 1, fStack->fActiveLinks));
	      fSM.fClean = kFALSE;
	      fHC->fCorrupted += 100;
	      SeekEndOfData();	

	      if (fgWarnError) 
		{
		  AliError(Form("Failed HC : %s", DumpHCinfoH0(fHC)));
		  AliError(Form("Failed HC : %s", DumpHCinfoH1(fHC)));
		  AliError(Form("Debug Event Counter : %d", fgStreamEventCounter)); 
		}
	      // let us assume that we have the HC data although link mask says different
	      // recovery
	      if (fStackLinkNumber == -1)
		{
		  if (fgWarnError) AliWarning("Trying to recover to the right Link Mask.");
		  ilink -= 1;
		}
	      	      
	      continue;
	      //return kFALSE;
	    }
	  else
	    {
	      // just finish off with the end of data markers
	      SeekEndOfData();
	    }

	  //if (fgDebugFlag)  AliDebug(5, Form("++++ Done with link %d of %d", ilink + 1, fStack->fActiveLinks));

	} // ilink
    } // istack

  // for next()
  ResetIterators();

  if (fSM.fClean == kTRUE)
    return kTRUE;
  
  if (fgCleanDataOnly && (fSM.fClean == kFALSE))
    {
      if (fgWarnError) AliWarning("Buffer with errors. Returning FALSE.");
      fSM.fActiveStacks = 0; // Next will not give data
      return kFALSE;
    }

  return kTRUE;
}

//------------------------------------------------------------
Int_t 
AliTRDrawStreamTB::DecodeSM()
{
  //
  // decode SM data in case AliRawReader is in use

  if (fRawReader)
    {      
      Int_t nextBuff = NextBuffer();
      while (nextBuff != -1)
	{
	  if (nextBuff > 0)
	    return nextBuff;	   	  
	  nextBuff = NextBuffer();
	}
      return -1;
    }
  else
    {
      AliWarning("AliRawReader not set.");
    }
  return kFALSE;
}

//------------------------------------------------------------
Int_t 
AliTRDrawStreamTB::DecodeSM(AliRawReader *reader)
{
  //
  // decode SM with the AliRawReader
  //
  if (reader != 0)
    {
      fRawReader = reader;
      return DecodeSM();
    }
  else
    {
      AliWarning("Argument AliRawReader is 0.");
    }
  return kFALSE;
}

//------------------------------------------------------------
Bool_t 
AliTRDrawStreamTB::SeekEndOfData()
{
  //
  // go to end of data marker
  //
  Int_t fEndOfDataCount = 0;

  while ( *fpPos != ENDOFRAWDATAMARKER && fpPos < fpEnd )
    {
      fpPos++;
    }
  
  while (*fpPos == ENDOFRAWDATAMARKER && fpPos < fpEnd )
    {
      fEndOfDataCount++;
      fpPos++;      
    }
  
//   if (fEndOfDataCount == 0)
//     {
//       if (fgWarnError) AliWarning("End of buffer reached first. No end of data marker?");
//       return kFALSE;
//     }

  return kTRUE;
}

//------------------------------------------------------------
Bool_t 
AliTRDrawStreamTB::SeekNextMCMheader()
{
  //
  // go to mcm marker
  //

  // move back to the mcm pos
  fpPos = fMCM->fPos + 1;

  while ( *fpPos != ENDOFRAWDATAMARKER && fpPos < fpEnd )
    {
      //if (CheckMCMMask(fpPos) == 0)
      if (MCM_HEADER_MASK_ERR(*fpPos) == 0)      
	{
	  if (fgWarnError) AliWarning(Form("^^^ Found : Pos 0x%08x : Val 0x%08x", fpPos, *fpPos));
	  return kTRUE;
	}
      fpPos++;
    }

  SeekEndOfData();
  return kFALSE;
}

//------------------------------------------------------------

Bool_t 
AliTRDrawStreamTB::DecodeTracklets()
{
  //
  // decode tracklets
  //

  fLinkTrackletCounter = 0;
  fEndOfTrackletCount = 0;

  if (fgDebugFlag)  AliDebug(10, Form("Decode tracklets at 0x%08x : 0x%08x", fpPos, *fpPos));
  //cout << Form("Decode tracklets..") << endl;;
  //printf("[I] Tracklet... ? \n");

  while ( *fpPos != END_OF_TRACKLET_MARKEROLD && *fpPos != END_OF_TRACKLET_MARKERNEW && fpPos < fpEnd )
    {

      if (fgDebugFlag)  AliDebug(10, Form("Tracklet found at 0x%08x : 0x%08x", fpPos, *fpPos));

      //decode tracklets here...
      fLinkTrackletCounter++;

      if (fLinkTrackletCounter > MAX_TRACKLETS_PERHC)
	{
	  if (fgWarnError) AliWarning(Form("Max number of tracklets exceeded %d > %d. Something is wrong with the input data!", 
			  fLinkTrackletCounter, MAX_TRACKLETS_PERHC));

	  if (fRawReader) fRawReader->AddMajorErrorLog(kTrackletOverflow,"Too many tracklets"); 
	  return kFALSE;
	}

      if (fgDebugStreamer)
	{
	  // jan!
	  unsigned pid;
	  unsigned row;
	  signed defl;
	  signed ypos;
	  
	  double deflm;
	  double yposm;
	  
	  unsigned long val;
	  
	  val = *fpPos;

	  pid = val >> 24;
	  row = (val >> 20) & 0xF;
	  defl = (val >> 13) & 0x7F;
	  defl = defl << 25 >> 25;
	  ypos = val & 0x1FFF;
	  ypos = ypos << 19 >> 19;
	  
	  deflm = defl * 140.0e-6; /* m */
	  yposm = ypos * 160.0e-6; /* m */
	  
	  TTreeSRedirector &cstream = *fgDebugStreamer;
	  cstream << "MCMtracklets"
		  << "Event=" << fgStreamEventCounter
		  << ".stack=" << fStackNumber
		  << ".link=" << fStackLinkNumber
		  << ".pid=" << pid
		  << ".row=" << row
		  << ".defl=" << defl
		  << ".ypos=" << ypos
		  << ".deflm=" << deflm
		  << ".yposm=" << yposm
		  << "\n";
	}
      fpPos++;
    }

  while ( ( *fpPos == END_OF_TRACKLET_MARKEROLD || *fpPos == END_OF_TRACKLET_MARKERNEW ) && fpPos < fpEnd )
    {
      //seek non end of tracklets
      //printf("[I] End of tracklets. \n");
      if (fgDebugFlag)  AliDebug(10, Form("EoTracklets found at 0x%08x : 0x%08x", fpPos, *fpPos));

      fEndOfTrackletCount++;
      fpPos++;
    }

  if ( fEndOfTrackletCount < 2 )
    {
      //if (fgWarnError) AliWarning(Form("End of tracklet marker words missing %d", 2 - fEndOfTrackletCount));
      if (fRawReader) fRawReader->AddMajorErrorLog(kEOTrackeltsMissing, "End of tracklets word missing"); 
      return kFALSE;
    }

  return kTRUE;
}

//------------------------------------------------------------
Bool_t 
AliTRDrawStreamTB::IsRowValid()
{
  //SLOW GEOM
  if ( (fHC->fStack == 2 && fMCM->fROW >= fGeometry->RowmaxC0()) ||
       (fHC->fStack != 2 && fMCM->fROW >= fGeometry->RowmaxC1()) || fMCM->fROW < 0 ) 
    {
      if (fgWarnError) AliWarning(Form("SM%d L%dS%d: Wrong Padrow (%d) fROB=%d, fSIDE=%d, fMCM=%02d"
		      , fHC->fSM, fHC->fLayer, fHC->fStack, fMCM->fROW, fMCM->fROB, fHC->fSide, fMCM->fMCM ));
      if (fRawReader) fRawReader->AddMajorErrorLog(kWrongPadrow, "Wrong Row");
      return kFALSE;
    }
  return kTRUE;
}

//------------------------------------------------------------
Bool_t 
AliTRDrawStreamTB::IsMCMheaderOK()
{
  //
  // check the mcm header
  //

  if ( fMCM->fCorrupted > 10 )
    {
      if (fRawReader) fRawReader->AddMajorErrorLog(kMCMheaderCorrupted,"ADC mask Corrupted"); 
      if (fgWarnError) AliWarning(Form("Wrong ADC Mask word 0x%08x %s. Error : %d", *fMCM->fPos, DumpMCMadcMask(fMCM), fMCM->fCorrupted));
      return kFALSE;
    }

  if ( fMCM->fCorrupted > 0 )
    {
      if (fRawReader) fRawReader->AddMajorErrorLog(kMCMheaderCorrupted,"Corrupted"); 
      if (fgWarnError) AliWarning(Form("Wrong MCM word 0x%08x %s. Error : %d", *fMCM->fPos, DumpMCMinfo(fMCM), fMCM->fCorrupted));
      return kFALSE;
    }

  if ( fMCM->fMCM < 0 || fMCM->fMCM > 15 || fMCM->fROB < 0 || fMCM->fROB > 7 ) 
    {
      if (fRawReader) fRawReader->AddMajorErrorLog(kWrongMCMorROB,"Wrong ROB or MCM"); 
      if (fgWarnError) AliWarning(Form("Wrong fMCM or fROB. %s Skip this data.", DumpMCMinfo(fMCM)));
      return kFALSE;
    }

  if (IsRowValid() == kFALSE)
    return kFALSE;

  return kTRUE;
}
//------------------------------------------------------------
Bool_t 
AliTRDrawStreamTB::DecodeMCMheader()
{
  //
  // decode the mcm header
  //

  DecodeMCMheader(fpPos, fMCM);

  fMCM->fROW = fTRDfeeParam->GetPadRowFromMCM(fMCM->fROB, fMCM->fMCM);

  if (fHC->fRawVMajor > 2)
    {
      fpPos++;
      if ( fpPos < fpEnd )
	{
	  DecodeMask(fpPos, fMCM);
	  MCMADCwordsWithTbins(fHC->fTimeBins, fMCM);
	  fMCM->fAdcDataPos = fpPos + 1;
	}
      else
	{
	  if (fgWarnError) AliWarning("Expected ADC mask word. Fail due to buffer END.");	  
	  if (fRawReader) fRawReader->AddMajorErrorLog(kMCMADCMaskMissing,"Missing"); 
	  return kFALSE;
	}
    }
  else
    {
      UInt_t dummyMask = MCM_DUMMY_ADCMASK_VAL;
      DecodeMask(&dummyMask, fMCM);
      MCMADCwordsWithTbins(fHC->fTimeBins, fMCM);
      fMCM->fAdcDataPos = fpPos + 1;
    }

  if (fgDebugFlag)  
    {
      AliDebug(6, DumpMCMinfo(fMCM));
      AliDebug(7, DumpMCMadcMask(fMCM));
    }

  if (IsMCMheaderOK() == kFALSE)
    {
      return kFALSE;
    }

  if (fEventCounter == 0)
    {
      fEventCounter = fMCM->fEvCounter;
    }
  
  if (fEventCounter < fLastEventCounter)
    {
      if (fgWarnError) AliWarning(Form("Weird. Event from the past? Current %d Last %d %s", fEventCounter, fLastEventCounter, DumpMCMinfo(fMCM)));
      if (fRawReader) fRawReader->AddMajorErrorLog(kMCMeventMissmatch, "Wrong MCM event counter ? Past-Future"); 
      return kFALSE;
    }
  
  if (fEventCounter != fMCM->fEvCounter)
    {
      if (fgWarnError) AliWarning(Form("Weird. Event missmatch? %d %s", fEventCounter, DumpMCMinfo(fMCM)));
      if (fRawReader) fRawReader->AddMajorErrorLog(kMCMeventMissmatch, "Wrong MCM event counter ?"); 
      return kFALSE;
    }

  //fLastEventCounter = fEventCounter;
  //AliInfo(Form("Current MCM %s", DumpMCMinfo(fMCM)));
    
  return kTRUE;
  //return IsMCMheaderOK();
}

//------------------------------------------------------------
Bool_t 
AliTRDrawStreamTB::IsHCheaderOK()
{
  //
  // check the half chamber header
  //

  if (fHC->fCorrupted > 0)
    {
      if (fgWarnError) AliWarning(Form("Wrong HC Header word. Word 0x%08x Error : %d", *fHC->fPos, fHC->fCorrupted));
      if (fRawReader) fRawReader->AddMajorErrorLog(kHCHeaderCorrupt, "Corrupted"); 

      return kFALSE;
    }

  if (fHC->fStack < 0 || fHC->fStack > 4)
    {
      if (fgWarnError) AliWarning(Form("Wrong Stack %d", fHC->fStack));
      if (fRawReader) fRawReader->AddMajorErrorLog(kHCHeaderWrongStack, "Wrong Stack");       
      return kFALSE;
    }

  if (fHC->fStack != fStackNumber) 
    {
      if (fgWarnError) AliWarning(Form("Missmatch: Stack in HC header %d HW-stack %d", 
				       fHC->fStack, fStackNumber));
      if (fRawReader) fRawReader->AddMajorErrorLog(kHCHeaderWrongStack, "Stack-HWstack");       
      // Try this for recovery in DecodeSM(void*,UInt_t) after DecodeHC failed
      // buffer will still will be marked as NOT clean
      fStackNumber = -1;
      return kFALSE;
    }

  if (fHC->fLayer < 0 || fHC->fLayer >= AliTRDgeometry::kNplan)
    {
      if (fgWarnError) AliWarning(Form("Wrong plane(layer) %d", fHC->fLayer));
      if (fRawReader) fRawReader->AddMajorErrorLog(kHCHeaderWrongLayer, "Wrong Plane"); 
      return kFALSE;
    }

  if ((fHC->fLayer * 2 != fStackLinkNumber) && (fHC->fLayer * 2 != fStackLinkNumber - 1)) 
    {
      if (fgWarnError) AliWarning(Form("Missmatch: plane(layer) in HCheader %d HW-Link %d | %s", 
				       fHC->fLayer, fStackLinkNumber, DumpStackInfo(fStack)));
      if (fRawReader) fRawReader->AddMajorErrorLog(kHCHeaderWrongLayer, "Plane-Link missmatch"); 
      // Try this for recovery in DecodeSM(void*,UInt_t) after DecodeHC failed
      // buffer will still will be marked as NOT clean
      fStackLinkNumber = -1;
      return kFALSE;      
    }

  if (fHC->fSide < 0 || fHC->fSide > 1)
    {
      if (fgWarnError) AliWarning(Form("Wrong Side %d", fHC->fSide));
      if (fRawReader) fRawReader->AddMajorErrorLog(kHCHeaderWrongSide, "Wrong Side");       
      return kFALSE;
    }
  
  // SLOW GEOM
  fHC->fDET = fGeometry->GetDetector(fHC->fLayer, fHC->fStack, fHC->fSM);
  if (fHC->fDET < 0 || fHC->fDET >= AliTRDgeometry::kNdet)
    {
      if (fgWarnError) AliWarning(Form("Wrong detector %d", fHC->fDET));      
      if (fRawReader) fRawReader->AddMajorErrorLog(kHCHeaderWrongDet, "Wrong Det");       
      return kFALSE;
    }

  // SLOW GEOM
  // this check fails - raw data inconsistent with geometry?
  if (fHC->fSM != fGeometry->GetSector(fHC->fDET)
      || fHC->fSM <0 || fHC->fSM >= AliTRDgeometry::kNsect)
    {
      if (fgWarnError) AliWarning(Form("Wrong SM(sector) %d (Geometry says: %d) Stack=%d Layer=%d Det=%d", 
				       fHC->fSM, fGeometry->GetSector(fHC->fDET),
				       fHC->fStack, fHC->fLayer, fHC->fDET));      
      if (fRawReader) fRawReader->AddMajorErrorLog(kHCHeaderWrongSM, "Wrong SM");       
      return kFALSE;
    }

  // SLOW GEOM
  // CPU EXPENSIVE!!!
  fHC->fROC    = fGeometry->GetDetectorSec(fHC->fLayer, fHC->fStack);
  if (fHC->fROC < 0)
    {
      if (fRawReader) fRawReader->AddMajorErrorLog(kHCHeaderWrongROC, "Wrong ROC");       
      return kFALSE;
    }

  fHC->fRowMax = fGeometry->GetRowMax(fHC->fLayer, fHC->fStack, fHC->fSM);
  if (fHC->fRowMax < 1)
    {
      if (fRawReader) fRawReader->AddMajorErrorLog(kHCHeaderWrongROC, "Wrong ROC Row Max");       
      return kFALSE;
    }

  fHC->fColMax = fGeometry->GetColMax(fHC->fROC);
  if (fHC->fColMax < 1)
    {
      if (fRawReader) fRawReader->AddMajorErrorLog(kHCHeaderWrongROC, "Wrong ROC Col Max");       
      return kFALSE;
    }
  // extra - not needed
  //   if (fHC->fSM <0 || fHC->fSM >= AliTRDgeometry::kNsect)
  //     {
  //       if (fgWarnError) AliWarning(Form("Wrong SM(sector) %d (Geometry says: %d) Stack=%d Layer=%d", 
  // 		      fHC->fSM, fGeometry->GetDetectorSec(fHC->fLayer, fHC->fStack),
  // 		      fHC->fStack, fHC->fLayer));      
  //       return kFALSE;
  //     }
  
  return kTRUE;
}

//------------------------------------------------------------
Bool_t 
AliTRDrawStreamTB::DecodeHCheader()
{  
  //
  // decode the half chamber header
  //

  DecodeHCwordH0(fpPos, fHC);
    
  if (fHC->fNExtraWords > 0)
    {
      fpPos++;
      if (fpPos < fpEnd)
	{
	  DecodeHCwordH1(fpPos, fHC);
	}
      else
	{
	  if (fgWarnError) AliWarning("Expected HC header word 1. Fail due to buffer END.");
	  if (fRawReader) fRawReader->AddMajorErrorLog(kHCWordMissing,"Next HC word 1 (count from 0) missing"); 
	  return kFALSE;
	}
    }

  if (fgDebugFlag)  AliDebug(5, DumpHCinfoH0(fHC));
  if (fgDebugFlag)  AliDebug(5, DumpHCinfoH1(fHC));

  fHC->fDET = -1;
  if (IsHCheaderOK() == kFALSE)
    {
      return kFALSE;
    }
  
  return kTRUE;
}

//------------------------------------------------------------
Bool_t 
AliTRDrawStreamTB::DecodeHC()
{
  //
  // decode the hc data
  // function for devel & debug
  //

  if (DecodeHCheader() == kFALSE)
    {
      if (fgDebugStreamer)
	{
	  TTreeSRedirector &cstream = *fgDebugStreamer;
	  cstream << "HCDecodeError"
		  << "Event=" << fgStreamEventCounter
		  << ".stack=" << fStackNumber
		  << ".link=" << fStackLinkNumber
		  << ".hcCorrupted" << fHC->fCorrupted
		  << ".hcRVmajor=" << fHC->fRawVMajor
		  << ".hcDCSboard=" << fHC->fDCSboard
		  << ".hcSM=" << fHC->fSM
		  << ".hcStack=" << fHC->fStack
		  << ".hcLayer=" << fHC->fLayer
		  << ".hcSide=" << fHC->fSide
		  << ".hcTbins=" << fHC->fTimeBins
		  << ".hcBunchCross=" << fHC->fBunchCrossCounter
		  << ".hcPreTrig=" << fHC->fPreTriggerCounter
		  << ".hcPreTrigPh=" << fHC->fPreTriggerPhase
		  << ".hcDet=" << fHC->fDET
		  << ".hcROC=" << fHC->fROC
		  << ".hc=RowMax" << fHC->fRowMax
		  << ".hc=ColMax" << fHC->fColMax
		  << "\n";
	}
      
      if (fgWarnError) AliWarning("Header decode failed.");
      return kFALSE;
    }
  else
    {
      fpPos++;
      if (fpPos >= fpEnd)
	{
	  if (fgWarnError) AliWarning("No MCM data? Not enough data in the buffer.");
	  if (fRawReader) fRawReader->AddMajorErrorLog(kMCMdataMissing, "MCM data missing"); 
	  return kFALSE;
	}
    }

  fHC->fMCMmax = 0;
  while (*fpPos != ENDOFRAWDATAMARKER && fpPos < fpEnd)
    {
      if (fHC->fMCMmax > TRD_MAX_MCM)
	{
	  if (fgWarnError) AliWarning("More mcm data than expected!");
	  if (fRawReader) fRawReader->AddMajorErrorLog(kMCMoverflow, "Too many mcms found!"); 
	  return kFALSE;
	}

      fMCM = &fHC->fMCMs[fHC->fMCMmax];

      if (DecodeMCMheader() == kFALSE)
	{
	  if (fgDebugStreamer)
	    {
	      TTreeSRedirector &cstream = *fgDebugStreamer;
	      cstream << "MCMDecodeError"
		      << "Event=" << fgStreamEventCounter
		      << ".stack=" << fStackNumber
		      << ".link=" << fStackLinkNumber
		      << ".hcDCSboard=" << fHC->fDCSboard
		      << ".hcSM=" << fHC->fSM
		      << ".hcStack=" << fHC->fStack
		      << ".hcLayer=" << fHC->fLayer
		      << ".hcSide=" << fHC->fSide
		      << ".hcBunchCross=" << fHC->fBunchCrossCounter
		      << ".hcPreTrig=" << fHC->fPreTriggerCounter
		      << ".hcPreTrigPh=" << fHC->fPreTriggerPhase
		      << ".hcDet=" << fHC->fDET

		      << ".mcmCorrupted=" << fMCM->fCorrupted
		      << ".mcmROB=" << fMCM->fROB
		      << ".mcmMCM=" << fMCM->fMCM
		      << ".mcmADCMaskWord=" << fMCM->fADCMaskWord
		      << ".mcmEventCounter=" << fMCM->fEvCounter

		      << ".fEventCounter=" << fEventCounter 
		      << ".fLastEventCounter=" << fLastEventCounter
		      << "\n";
	    }
	  
	  return kFALSE;
	}

      fHC->fMCMmax++;

      if (fMCM->fADCmax > 0)
	{
	  fpPos++;
	  if (fpPos >= fpEnd)
	    {
	      if (fgDebugFlag)  AliDebug(9, Form("Buffer short of data. ADC data expected."));	  
	      if (fRawReader) fRawReader->AddMajorErrorLog(kADCdataMissing, "ADC data missing"); 
	    }
	  
	  //	  for (Int_t iadc = 0; iadc < TRD_MAX_ADC; iadc++)
	  for (Int_t iadc = 0; iadc < fMCM->fADCmax; iadc++)
	    {
	      fADC = &fMCM->fADCs[iadc];
	      fADC->fADCnumber = fMCM->fADCchannel[iadc];

	      if (fgDebugFlag)  AliDebug(9, Form("This is ADC %d of %d. ADC number is %d.", 
						iadc+1, fMCM->fADCmax, fMCM->fADCchannel[iadc]));

	      if (fpPos + fMCM->fSingleADCwords >= fpEnd)
		{
		  if (fgDebugStreamer)
		    {
		      TTreeSRedirector &cstream = *fgDebugStreamer;
		      cstream << "ADCDecodeError"
			      << "Event=" << fgStreamEventCounter
			      << ".hcSM=" << fHC->fSM
			      << ".hcStack=" << fHC->fStack
			      << ".hcLayer=" << fHC->fLayer
			      << ".mcmROB=" << fMCM->fROB
			      << ".mcmMCM=" << fMCM->fMCM
			      << ".mcmADCMaskWord=" << fMCM->fADCMaskWord
			      << "\n";			
		    }
		  
		  if (fgWarnError) AliWarning(Form("This is ADC %d of %d. ADC number is %d.", iadc+1, fMCM->fADCmax, fMCM->fADCchannel[iadc]));
		  if (fgWarnError) AliWarning("--> ADC (10 words) expected. Not enough data in the buffer.");
		  if (fgWarnError) AliWarning(Form("--> ADC mask : %s", DumpMCMadcMask(fMCM)));

		  if (fRawReader) fRawReader->AddMajorErrorLog(kADCdataMissing, "ADC data missing - less than expected"); 
		  return kFALSE;
		}

	      //DECODE the data here
	      if (DecodeADC() == kFALSE)
		{
		  // check if we are out of the det when the pad is shared
		  if (fADC->fIsShared && fADC->fCorrupted == 111)
		    {
		      fADC->fCOL = -1;
		      fpPos = fADC->fPos + fMCM->fSingleADCwords;
		    }
		  else
		    {
		      if (fgWarnError) AliWarning(Form("ADC decode failed."));
		      //fpPos = fMCM->fPpos + fMCMADCWords;
		      fpPos = fADC->fPos + fMCM->fSingleADCwords;
		      return kFALSE;
		    }
		}
	      //decode the ADC words here
	    }	  
	} //mcm data present
      else
	{
	  fpPos++;
	}
    }//while eof data

  if (fpPos >= fpEnd)
    {
      if (fgWarnError) AliWarning("We are at the end of buffer. There should be one more word left.");
      //SeekEndOfData();  
      return kFALSE;
    }

  return kTRUE;
}

//------------------------------------------------------------
Bool_t
AliTRDrawStreamTB::DecodeADC()
{
  //
  // decode single ADC channel
  //

  fADC->fCorrupted = 0;
  fMaskADCword = ADC_WORD_MASK(*fpPos);
  fADC->fPos = fpPos;

  fTbinADC = 0;

  for (Int_t i = 0; i < TRD_MAX_TBINS; i++)
    fADC->fSignals[i] = 0;

  for (Int_t iw = 0; iw < fMCM->fSingleADCwords; iw++)
    {
      if (fMaskADCword != ADC_WORD_MASK(*fpPos))
	{
	  //this is corrupted data
	  fADC->fCorrupted += 1;
	  if (fgWarnError) AliWarning(Form("Mask Change in ADC data 0x%08x 0x%08x", *(fpPos-1), fMaskADCword, *fpPos, ADC_WORD_MASK(*fpPos)));
	  if (fRawReader) fRawReader->AddMajorErrorLog(kADCmaskMissmatch, "Mask change inside single channel"); 

	  break;
	}
      
      fADC->fSignals[fTbinADC + 0] = (*fpPos & 0x00000ffc) >> 2;
      fADC->fSignals[fTbinADC + 1] = (*fpPos & 0x003ff000) >> 12;
      fADC->fSignals[fTbinADC + 2] = (*fpPos & 0xffc00000) >> 22;

      fTbinADC += 3;
      fpPos++;
    }

  if (fADC->fADCnumber <= 1 || fADC->fADCnumber == fMaxADCgeom - 1)
    {
      fADC->fIsShared = kTRUE;
    }
  else
    {
      fADC->fIsShared = kFALSE;
    }

  if ( fADC->fADCnumber >= fMaxADCgeom - 1)
    {
      // let us guess the Column
      // take the one before last ADC and shift by one column
      // later we check if we are inside the limits of the chamber
      fADC->fCOL = AliTRDfeeParam::Instance()->GetPadColFromADC(fMCM->fROB, fMCM->fMCM, fADC->fADCnumber - 1);
      fADC->fCOL--;
    }
  else
    {
      fADC->fCOL = fTRDfeeParam->GetPadColFromADC(fMCM->fROB, fMCM->fMCM, fADC->fADCnumber);
    }

  if (fADC->fCOL >= fHC->fColMax || fADC->fCOL < 0)
    {
      if (fADC->fIsShared == kFALSE)
	{
	  fADC->fCorrupted += 100;

	  if (fgWarnError) AliWarning(Form("Wrong column! ADCnumber %d MaxIs %d Col %d MaxIs %d MCM= %s", 
					   fADC->fADCnumber, fMaxADCgeom, fADC->fCOL, fHC->fColMax, DumpMCMinfo(fMCM)));
	  if (fRawReader) fRawReader->AddMajorErrorLog(kWrongPadcolumn, "Wrong column"); 	  
	}
      else
	{
	  // flag it - we are out of the det when the pad is shared
	  if (fgDebugFlag) AliDebug(10, Form("Column out of the detector! ADCnumber %d MaxIs %d Col %d MaxIs %d MCM= %s", 
					     fADC->fADCnumber, fMaxADCgeom, fADC->fCOL, fHC->fColMax, DumpMCMinfo(fMCM)));
	  fADC->fCorrupted += 111;
	}
    }

  if (fADC->fCorrupted > 0)
    {
      return kFALSE;
    }

  fDecodedADCs++;
  return kTRUE;
}

//--------------------------------------------------------
void AliTRDrawStreamTB::DecodeSMInfo(const UInt_t *word, struct AliTRDrawSM *sm) const
{
  //
  // Decode the data *word into SM info structure
  //
    
  sm->fPos = (UInt_t*)word;

  // do it once here
  UInt_t vword = *word;
  sm->fHeaderSize = SM_HEADER_SIZE(vword);
    
  if (TRACKLETS_ENABLED(vword) > 0)
    sm->fTrackletEnable = kTRUE;
  else
    sm->fTrackletEnable = kFALSE;
    
  UInt_t stackMask = STACK_MASK(vword);
  sm->fActiveStacks = 0;
  for (Int_t i = 0; i < 5; i++)
    {
      //if ((stackMask >> i) & 0x1 > 0)
      if (IS_BIT_SET(stackMask,i) > 0)
	{
	  sm->fStackActive[i] = kTRUE;
	  sm->fActiveStacks++;
	}
      else
	{
	  sm->fStackActive[i] = kFALSE;
	}
    }
}

//--------------------------------------------------------
const char *AliTRDrawStreamTB::DumpSMInfo(const struct AliTRDrawSM *sm)
{
  //
  // Get SM structure into a const char
  //
    
  return Form("[ SM Info 0x%08x] : Hsize %d TrackletEnable %d Stacks %d %d %d %d %d",
	      *sm->fPos,
	      sm->fHeaderSize, sm->fTrackletEnable,
	      sm->fStackActive[0], sm->fStackActive[1], sm->fStackActive[2],
	      sm->fStackActive[3], sm->fStackActive[4]);      
}

//--------------------------------------------------------
void AliTRDrawStreamTB::DecodeStackInfo(const UInt_t *word, struct AliTRDrawStack *st) const
{
  //
  // Decode stack info - active links
  //

  st->fPos = (UInt_t*)word;
      
  // do it once here
  UInt_t vword = *word;
  st->fHeaderSize = STACK_HEADER_SIZE(vword);

  UInt_t linkMask = STACK_LINK_WORD(vword);
  st->fActiveLinks = 0;
  for (Int_t i = 0; i < 12; i++)
    {
      //if (linkMask & (0x1 << i) > 0)
      if (IS_BIT_SET(linkMask,i) > 0)
	{
	  st->fLinksActive[i] = kTRUE;
	  st->fTrackletDecode[i] = kTRUE;
	  st->fHCDecode[i] = kTRUE;
	  st->fActiveLinks++;
	}
      else
	{
	  st->fTrackletDecode[i] = kFALSE;
	  st->fLinksActive[i] = kFALSE;
	  st->fHCDecode[i] = kFALSE;
	}
    }
}
  
//--------------------------------------------------------
const char *AliTRDrawStreamTB::DumpStackInfo(const struct AliTRDrawStack *st)
{
  //
  // format the string with the stack info
  //
  return Form("[ Stack Info 0x%08x ] : Hsize %d Links Active %d %d %d %d %d %d %d %d %d %d %d %d",
	      *st->fPos,
	      st->fHeaderSize,
	      st->fLinksActive[0], st->fLinksActive[1], st->fLinksActive[2], st->fLinksActive[3],
	      st->fLinksActive[4], st->fLinksActive[5], st->fLinksActive[6], st->fLinksActive[7],
	      st->fLinksActive[8], st->fLinksActive[9], st->fLinksActive[10], st->fLinksActive[11]);
}

//--------------------------------------------------------
void AliTRDrawStreamTB::DecodeHCwordH0(const UInt_t *word, struct AliTRDrawHC *hc) const
{
  //
  // decode the hc header word 0
  //

  // do it once here
  UInt_t vword = *word;

  hc->fCorrupted = HC_HEADER_MASK_ERR(vword);
  if (hc->fCorrupted > 0)
    {
      hc->fH0ErrorCounter++;
    }

  hc->fSpecialRawV =  HC_SPECIAL_RAW_VERSION(vword);
  hc->fRawVMajor = HC_MAJOR_RAW_VERSION(vword);
  hc->fRawVMinor = HC_MINOR_RAW_VERSION(vword);
  hc->fNExtraWords = HC_EXTRA_WORDS(vword);
  hc->fDCSboard = HC_DCS_BOARD(vword);
  hc->fSM = HC_SM_NUMBER(vword);
  hc->fStack = HC_STACK_NUMBER(vword);
  hc->fLayer = HC_LAYER_NUMBER(vword);
  hc->fSide = HC_SIDE_NUMBER(vword);

  hc->fPos[0] = (UInt_t*)word;

}

//--------------------------------------------------------
void AliTRDrawStreamTB::DecodeHCwordH1(const UInt_t *word, struct AliTRDrawHC *hc) const
{
  //
  // 
  //

  // do it once here
  UInt_t vword = *word;

  hc->fCorrupted += 10 * HC_HEADER_MASK_ERR(vword);
  if (hc->fCorrupted > 10)
    {
      hc->fH1ErrorCounter++;
    }

  hc->fTimeBins = HC_NTIMEBINS(vword);
  hc->fBunchCrossCounter = HC_BUNCH_CROSS_COUNTER(vword);
  hc->fPreTriggerCounter = HC_PRETRIGGER_COUNTER(vword);
  hc->fPreTriggerPhase = HC_PRETRIGGER_PHASE(vword);

  hc->fPos[1] = (UInt_t*)word;
}
  
//--------------------------------------------------------
const char *AliTRDrawStreamTB::DumpHCinfoH0(const struct AliTRDrawHC *hc)
{
  //
  // 
  //
  if (!hc)
    return Form("Unable to dump. Null received as parameter!?!");
  else
    return Form("[ HC[0] at 0x%08x ] : 0x%08x Info is : RawV %d SM %d Stack %d Layer %d Side %d DCSboard %d",
		hc->fPos[0], *(hc->fPos[0]), hc->fRawVMajor, hc->fSM, hc->fStack, hc->fLayer, hc->fSide, hc->fDCSboard);
}

//--------------------------------------------------------
const char *AliTRDrawStreamTB::DumpHCinfoH1(const struct AliTRDrawHC *hc)
{
  //
  // 
  //
  if (!hc)
    return Form("Unable to dump. Null received as parameter!?!");
  else
    return Form("[ HC[1] at 0x%08x ] : 0x%08x Info is : TBins %d BCcount %d PreTrigCount %d PreTrigPhase %d",
		hc->fPos[1], *(hc->fPos[1]), hc->fTimeBins, hc->fBunchCrossCounter, hc->fPreTriggerCounter, hc->fPreTriggerPhase);
}

//--------------------------------------------------------
void AliTRDrawStreamTB::DecodeMCMheader(const UInt_t *word, struct AliTRDrawMCM *mcm) const
{
  //
  // decode the mcm header
  //

  // do it once here
  UInt_t vword = *word;

  mcm->fCorrupted = MCM_HEADER_MASK_ERR(vword);
  //printf("0x%08x %d\n", word, MCM_HEADER_MASK_ERR(word));
  if (mcm->fCorrupted)
    mcm->fErrorCounter++;
  mcm->fROB = MCM_ROB_NUMBER(vword);
  mcm->fMCM = MCM_MCM_NUMBER(vword);
  mcm->fEvCounter = MCM_EVENT_COUNTER(vword);

  mcm->fPos = (UInt_t*)word;
}

//--------------------------------------------------------
UInt_t AliTRDrawStreamTB::GetMCMadcMask(const UInt_t *word, struct AliTRDrawMCM *mcm) const
{
  //
  // get the adc mask
  //

  // do it once here
  UInt_t vword = *word;

  mcm->fADCindex = 0;
  mcm->fADCmax   = 0;
  mcm->fADCMask  = 0;
  mcm->fADCMaskWord = vword;
  //memset(mcm->fADCchannel, 0, sizeof(UInt_t) * 30);
  //if ((word & 0x0000000F) == 0x0000000C)
  if ( MCM_ADCMASK_MASK_ERR(vword) == 0 )
    {
      //mcm->fADCMask = (word >> 4) & skADCmaskBits;
      mcm->fADCMask = MCM_ADCMASK_VAL(vword);
    }
  else
    {
      mcm->fADCMask = 0xffffffff;
      mcm->fCorrupted += 10;
      mcm->fMaskErrorCounter++;
    }

  return mcm->fADCMask;
}

//--------------------------------------------------------
void AliTRDrawStreamTB::DecodeMask(const UInt_t *word, struct AliTRDrawMCM *mcm) const
{
  //
  // decode the adc mask - adcs to be read out
  //

  mcm->fMCMADCWords = 0;
  mcm->fSingleADCwords = 0;
  mcm->fADCmax = 0;
  mcm->fADCMask = GetMCMadcMask(word, mcm);
  if (mcm->fADCMask > 0)
    {
      for (Int_t i = 0; i < TRD_MAX_ADC; i++)
	{
	  mcm->fADCchannel[mcm->fADCmax] = 0;
	  if( IS_BIT_SET(mcm->fADCMask,i) )
	    {
	      mcm->fADCchannel[mcm->fADCmax] = i;
	      mcm->fADCmax++;
	    }
	}
    }
}

//--------------------------------------------------------
void AliTRDrawStreamTB::MCMADCwordsWithTbins(UInt_t fTbins, struct AliTRDrawMCM *mcm) const
{
  //
  //  count the expected mcm words for a given tbins
  //
  mcm->fMCMADCWords = ( mcm->fADCmax ) * ( fTbins / 3 );
  mcm->fSingleADCwords = 0;
  if (mcm->fADCmax > 0)
    {
      mcm->fSingleADCwords = mcm->fMCMADCWords/mcm->fADCmax;
    }
}
  
//--------------------------------------------------------
const char *AliTRDrawStreamTB::DumpMCMinfo(const struct AliTRDrawMCM *mcm)
{
  //
  // mcm info in a string
  //
  if (!mcm)
    return Form("Unable to dump. Null received as parameter!?!");
  else
    return Form("[ MCM 0x%08x ] : ROB %d MCM %d EvCounter %d", *(mcm->fPos), mcm->fROB, mcm->fMCM, mcm->fEvCounter);
}
  
//--------------------------------------------------------
//TString DumpMCMadcMask(const struct AliTRDrawMCMInfo *mcm)
const char *AliTRDrawStreamTB::DumpMCMadcMask(const struct AliTRDrawMCM *mcm)
{
  //
  // mcm adc mask in a string
  //
  if (!mcm)
    return Form("Unable to dump. Null received as parameter!?!");

  TString tsreturn = Form("[Word] : 0x%08x => [Mask] : 0x%08x : ", mcm->fADCMaskWord, mcm->fADCMask);
  for (Int_t i = 0; i < 21; i++)
    {
      tsreturn += Form("%d ", mcm->fADCchannel[i]);
    }
  tsreturn += "";
  return tsreturn.Data();
}
