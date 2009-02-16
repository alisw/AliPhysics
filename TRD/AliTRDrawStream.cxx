/**************************************************************************
* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
q*                                                                        *
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

/* $Id: AliTRDrawStream.cxx 27797 2008-08-05 14:37:22Z cblume $ */

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

#include "AliTRDrawStream.h"
#include "AliTRDgeometry.h"
#include "AliTRDfeeParam.h"
#include "AliTRDdigitsManager.h"
#include "AliTRDarrayDictionary.h"
#include "AliTRDarrayADC.h"
#include "AliTRDSignalIndex.h"
//#include "AliTRDcalibDB.h" 
//#include "Cal/AliTRDCalPadStatus.h" //is now read in clusterizer
#include "AliTRDrawTPStream.h"

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
#define LINK0_DATA_TYPE_FLAG(w) (GET_VALUE_AT(w,0x3,4) == (0x0) ? 0 : 1) // 0 if physics data
#define LINK1_DATA_TYPE_FLAG(w) (GET_VALUE_AT(w,0x3,20) == (0x0) ? 0 : 1) // 0 if physics data
#define LINK0_MONITOR_FLAG(w) (GET_VALUE_AT(w,0xf,0) == (0x0) ? 0 : 1) // 0 if OK
#define LINK1_MONITOR_FLAG(w) (GET_VALUE_AT(w,0xf,16) == (0x0) ? 0 : 1) // 0 if OK

// HC word masks
//#define HC_HEADER_MASK_ERR(w) ( ((w) & (0x80000003)) == (0x80000001) ? 0 : 1) // 0 if OK!!!
#define HC_HEADER_MASK_ERR(w) ( ((w) & (0x3)) == (0x1) ? 0 : 1) // 0 if OK!!!

// HC word 0
#define HC_SPECIAL_RAW_VERSION(w) IS_BIT_SET(w,31)
#define HC_MAJOR_RAW_VERSION(w) GET_VALUE_AT(w,0x7f,24)
#define HC_MAJOR_RAW_VERSION_OPT(w) GET_VALUE_AT(w,0x7,24)
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

// MCM word and ADC mask
#define MCM_HEADER_MASK_ERR(w) ( ((w) & (0xf)) == (0xc) ? 0 : 1) // 0 if OK!!!
#define MCM_ADCMASK_MASK_ERR(w) ( ((w) & (0xf)) == (0xc) ? 0 : 1) // 0 if OK!!!
#define MCM_MCM_NUMBER(w) GET_VALUE_AT(w,0x0f,24)
#define MCM_ROB_NUMBER(w) GET_VALUE_AT(w,0x7,28)
#define MCM_EVENT_COUNTER(w) GET_VALUE_AT(w,0x00fffff,4)
#define MCM_ADCMASK_VAL(w) GET_VALUE_AT(w,0x1fffff,4)
#define MCM_ADCMASK_NADC(w) GET_VALUE_AT(w,0x1f,25)

#define MCM_DUMMY_ADCMASK_VAL 0x015fffffc  // updated 
#define ADCDATA_VAL1 0x2  // updated 
#define ADCDATA_VAL2 0x3  // updated 

//--------------------------------------------------------
#define ADC_WORD_MASK(w) ((w) & 0x3)
//--------------------------------------------------------


ClassImp(AliTRDrawStream)

Bool_t AliTRDrawStream::fgExtraSkip = kFALSE;
Bool_t AliTRDrawStream::fgSkipCDH = kFALSE;
Bool_t AliTRDrawStream::fgWarnError = kTRUE;
Bool_t AliTRDrawStream::fgCleanDataOnly = kFALSE;
Bool_t AliTRDrawStream::fgDebugFlag = kTRUE;
Bool_t AliTRDrawStream::fgEnableMemoryReset = kTRUE;
Bool_t AliTRDrawStream::fgStackNumberChecker = kTRUE;
Bool_t AliTRDrawStream::fgStackLinkNumberChecker = kTRUE;
Bool_t AliTRDrawStream::fgSkipData = kTRUE;
Bool_t AliTRDrawStream::fgEnableDecodeConfigData = kFALSE;
Int_t AliTRDrawStream::fgDumpHead = -1;
Int_t  AliTRDrawStream::fgCommonAdditive = 0;
Int_t AliTRDrawStream::fgEmptySignals[] = 
  {
    -1, -1, -1, -1, -1,  -1, -1, -1, -1, -1,  -1, -1, -1, -1, -1
    -1, -1, -1, -1, -1,  -1, -1, -1, -1, -1,  -1, -1, -1, -1, -1
  };
Short_t AliTRDrawStream::fgMCMordering[] =
  {
    12, 13, 14, 15, 8, 9, 10, 11, 4, 5, 6, 7, 0, 1, 2, 3 
  };
Short_t AliTRDrawStream::fgROBordering[] =
  {
    0, 1, 2, 3
  };
Int_t  AliTRDrawStream::fgLastHC = -1;
Int_t  AliTRDrawStream::fgLastROB = -1;
Int_t  AliTRDrawStream::fgLastIndex = -1;
Bool_t  AliTRDrawStream::fDumpingEnable = kFALSE;
Int_t  AliTRDrawStream::fDumpingSM = -1;
Int_t  AliTRDrawStream::fDumpingStack = -1;
Int_t  AliTRDrawStream::fDumpingLayer = -1;
Int_t  AliTRDrawStream::fDumpingROB = -1;
Int_t  AliTRDrawStream::fDumpingMCM = -1;


AliTRDrawStream::AliTRDrawStream()
  : AliTRDrawStreamBase()
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
  , fNWordsCounter(-1)
  , fMaskADCword(0)
  , fTbinADC(0)
  , fDecodedADCs(-1)
  , fEventCounter(0)
  , fLastEventCounter(0)
  , fSharedPadsOn(kFALSE)
  , fMaxADCgeom(0)
  , fBufferRead(0)
  , fGeometry(0)
  , fRawReader(0)
  , fTRDfeeParam(0)
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
AliTRDrawStream::AliTRDrawStream(AliRawReader *rawReader)
  : AliTRDrawStreamBase(rawReader)
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
  , fNWordsCounter(-1)
  , fMaskADCword(0)
  , fTbinADC(0)
  , fDecodedADCs(-1)
  , fEventCounter(0)
  , fLastEventCounter(0)
  , fSharedPadsOn(kFALSE)
  , fMaxADCgeom(0)
  , fBufferRead(0)
  , fGeometry(0)
  , fRawReader(rawReader)
  , fTRDfeeParam(0)
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

AliTRDrawStream::AliTRDrawStream(const AliTRDrawStream& /*st*/)
  : AliTRDrawStreamBase()
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
  , fNWordsCounter(-1)
  , fMaskADCword(0)
  , fTbinADC(0)
  , fDecodedADCs(-1)
  , fEventCounter(0)
  , fLastEventCounter(0)
  , fSharedPadsOn(kFALSE)
  , fMaxADCgeom(0)
  , fBufferRead(0)
  , fGeometry(0)
  , fRawReader(0)
  , fTRDfeeParam(0)
{
  //
  // Copy constructor
  // 
  AliError("Not implemeneted.");
}

//------------------------------------------------------------
Bool_t AliTRDrawStream::SetRawVersion(Int_t fraw)
{
  //
  // function provided for backward compatibility
  //
  AliWarning("Raw data version is read from raw data stream! No point of setting it in here.");
  fraw = 0; // avoid warnings
  return kFALSE;
}

//------------------------------------------------------------
AliTRDrawStream::~AliTRDrawStream()
{
  //
  // destructor
  //
  delete fGeometry;
}

//------------------------------------------------------------

AliTRDrawStream &
AliTRDrawStream::operator=(const AliTRDrawStream &)
{
  //
  // we are not using this functionality
  //
  AliFatal("May not use.");
  return *this;
}

//___________________________________________________________
void 
AliTRDrawStream::SwapOnEndian()
{
  //
  // Check the endian and swap if needed
  //
  int itemp = 1;
  char* ptemp = (char*) &itemp;
  if (ptemp[0] != 1)
    {
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
Bool_t
AliTRDrawStream::DumpWords(UInt_t *px, UInt_t iw, UInt_t marker)
{
    //UInt_t iEnd = iw;
    //if ( iw ==0 ) iEnd = (fpEnd - px)/sizeof(UInt_t); // if iw is 0, dump all words

  TString tsreturn = Form("\n[ Dump Sequence at 0x%08x ] : ", px);
  for (UInt_t i = 0; i < iw; i++)
  //for (UInt_t i = 0; i < iEnd; i++)
    {
      if ( iw != 0 && px + iw > fpEnd) return kFALSE;

      if (i % 8 == 0) tsreturn += "\n                              ";
      if (marker != 0 && marker == px[i]) tsreturn += Form(" *>0x%08x<* ", px[i]);
      else tsreturn += Form("0x%08x ", px[i]);
    }
  tsreturn += "\n";

  AliInfo(tsreturn.Data());

  return kTRUE;
}

//------------------------------------------------------------
Bool_t 
AliTRDrawStream::SkipWords(UInt_t iw)
{
  //
  // Skip words corresponding to iw
  //
  if ( fpPos + iw < fpEnd )
    {
      fpPos += iw;
      return kTRUE;
    }
  else
    {
      if (fgWarnError) AliWarning(Form("Skip %d words failed. %d available", iw, fpEnd - fpPos - 1));
      return kFALSE;
    }

  return kTRUE;
}

//------------------------------------------------------------
Bool_t 
AliTRDrawStream::SetReader(AliRawReader *reader)
{

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
AliTRDrawStream::NextBuffer()
{
  //
  // return -1 if no more buffers available
  // return  0 DecodeSM failed (clean data required for example) but still maybe more data to come
  // return  1 DecodeSM OK
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
              if (fgEnableMemoryReset) ResetMemory(); //[mj]
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
AliTRDrawStream::ResetCounters()
{
  //
  // reset some global counters
  //
  fBufferRead = kFALSE; // important to read buffer

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
AliTRDrawStream::ResetIterators()
{
  //
  // reset data which should be reset every sm
  //
  fStackNumber = 0;     // reset for Next() function 
  fStackLinkNumber = 0; // reset for Next() function
  fhcMCMcounter = 0;  
  fmcmADCcounter = 0;
}

//------------------------------------------------------------
void
AliTRDrawStream::ResetPerSM()
{
  //
  // reset every SM
  //

  fSM.fHeaderSize = 0;
  fSM.fTrackletEnable = kFALSE;
  fSM.fCorrupted = 0;
  fSM.fNexpectedHalfChambers = 0;
  fSM.fNexpectedHalfChambers = 0;
  fSM.fClean = kTRUE;
  fSM.fPos = NULL;
  for (Int_t i=0; i<5; i++){
    fSM.fStackActive[i] = kFALSE;
  }
}     

//------------------------------------------------------------
void
AliTRDrawStream::ResetPerStack()
{
  //
  // reset every Stack
  //

  fStack->fHeaderSize = 0;
  fStack->fActiveLinks = 0;
  fStack->fPos = NULL;
  for (Int_t i=0; i<12; i++){
    fStack->fLinksActive[i] = kFALSE;
    fStack->fLinksDataType[i] = 0;
    fStack->fLinksMonitor[i] = 0;
    fStack->fLinkMonitorError[i] = 0;
  }
}

//------------------------------------------------------------
void 
AliTRDrawStream::ResetPerHC()
{
  //
  // reset every HC
  //
  fEventCounter = 0;
  fHC->fTrackletError = 0;
  fHC->fNTracklets = 0;
  fHC->fSpecialRawV = 0;
  fHC->fRawVMajor = 0;
  fHC->fRawVMajorOpt = 0;
  fHC->fRawVMinor = 0;
  fHC->fNExtraWords = 0;
  fHC->fDCSboard = 0;
  fHC->fSM = 0;
  fHC->fStack = 0;
  fHC->fLayer = 0;
  fHC->fSide = 0;
  fHC->fTimeBins = 0;
  fHC->fBunchCrossCounter = 0;
  fHC->fPreTriggerCounter = 0;
  fHC->fPreTriggerPhase = 0;
  fHC->fDET = 0;
  fHC->fROC = 0;
  fHC->fRowMax = 0;
  fHC->fColMax = 0;
  fHC->fMCMmax = 0;

  fHC->fH0Corrupted = 0;
  fHC->fH1Corrupted = 0;
  fHC->fCorrupted = 0;
}

//------------------------------------------------------------
void
AliTRDrawStream::ResetPerMCM()
{
  //
  // reset every MCM 
  //
  fMCM->fROB = 0;
  fMCM->fMCM = 0;
  fMCM->fROW = 0;
  fMCM->fEvCounter = 0;
  fMCM->fADCMask = 0;
  fMCM->fADCMaskWord = 0;
  fMCM->fADCmax = 0;
  fMCM->fADCcount = 0;
  fMCM->fMCMADCWords = 0;
  fMCM->fSingleADCwords = 0;
  fMCM->fMCMhdCorrupted = 0;
  fMCM->fADCmaskCorrupted = 0;
  fMCM->fCorrupted = 0;
  fMCM->fPos = NULL;
  fMCM->fAdcDataPos = NULL;
  fMCM->fADCcounter = 0;

  for (Int_t i=0; i<21; i++){
    fMCM->fADCchannel[i] = 0;
  }
}

//------------------------------------------------------------
void
AliTRDrawStream::ResetPerADC()
{
  //
  // reset every ADC 
  //
  fADC->fPos = NULL;
  fADC->fADCnumber = 0;
  fADC->fCOL = 0;
  fADC->fIsShared = kTRUE;
  fADC->fCorrupted = 0;

  for (Int_t i=0; i<30; i++){
    fADC->fSignals[i] = 0;
  }
}

//------------------------------------------------------------
void
AliTRDrawStream::ResetMemory()
{                 
  //              
  // initialize all the data members to prevent read data
  // from previous buffer
  //              
  ResetPerSM();
  for (Int_t istack=0; istack<5; istack++){
    fStack = &fSM.fStacks[istack];
    ResetPerStack();
    for (Int_t ilink=0; ilink<12; ilink++){
        fHC = &fStack->fHalfChambers[ilink];
        ResetPerHC();
        for (Int_t imcm=0; imcm<12; imcm++){
          fMCM = &fHC->fMCMs[imcm];
          ResetPerMCM();
          for (Int_t iadc=0; iadc<12; iadc++){
              fADC = &fMCM->fADCs[iadc];
              ResetPerADC();
          }
        }
    }      
  }
}         

//------------------------------------------------------------

Bool_t 
AliTRDrawStream::Next()
{
  //
  // returns with true on next adc read
  // returns false on errors and end of buffer
  // 
if (fBufferRead)
  {

    while (fStackNumber < 5 && fSM.fActiveStacks > 0)
      {
        if(fSM.fStackActive[fStackNumber] == kTRUE)
    {
      fStack = &fSM.fStacks[fStackNumber];
      while (fStackLinkNumber < 12)
        {
          if (fStack->fLinksActive[fStackLinkNumber] == kTRUE && fStack->fLinksMonitor[fStackLinkNumber] == 0)
      {
        fHC = &fStack->fHalfChambers[fStackLinkNumber];
        if (!fHC)
          {
            AliError(Form("HC missing at stack %d link %d", fStackNumber, fStackLinkNumber));
            return kFALSE;
          }
        if (fHC->fCorrupted == 0 && (fHC->fH0Corrupted == 0 && fHC->fH1Corrupted == 0)) // if HC data corrupted(in any case), we don't read data at all from this HC 
          {
            while (fhcMCMcounter < fHC->fMCMmax)
        {
          fMCM = &fHC->fMCMs[fhcMCMcounter];
          if (!fMCM)
            {
              AliError(Form("HC missing at stack %d link %d atMCMslot %d", 
                fStackNumber, fStackLinkNumber, fhcMCMcounter));
              return kFALSE;
            }
          while(fmcmADCcounter < fMCM->fADCmax)
            {
              fADC = &fMCM->fADCs[fmcmADCcounter];
              if (!fADC)
          {
            AliError(Form("ADC missing at stack %d link %d MCMslot %d ADCslot %d", 
              fStackNumber, fStackLinkNumber, fhcMCMcounter, fmcmADCcounter));
            return kFALSE;
          }
              fmcmADCcounter++;
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
  } // fBufferRead

  // in case rawreader manages the mem buffers, go for the next buffer 
  if (fRawReader)
    {
      Int_t nextBuff = NextBuffer();
      while (nextBuff != -1)
    {
      if (nextBuff > 0)
              {
                fBufferRead = kTRUE;
          return Next();	   	  
              }
      nextBuff = NextBuffer();
    }
    }

  return kFALSE;
}

//------------------------------------------------------------
Int_t 
AliTRDrawStream::NextChamber(AliTRDdigitsManager *digitsManager, UInt_t **trackletContainer) 
{
  //
  // Fills single chamber digit array 
  // Return value is the detector number
  //

  //AliTRDcalibDB *cal = AliTRDcalibDB::Instance();
  AliTRDarrayADC *digits = 0;
  AliTRDarrayDictionary *track0 = 0;
  AliTRDarrayDictionary *track1 = 0;
  AliTRDarrayDictionary *track2 = 0; 
  AliTRDSignalIndex *indexes = 0;

  // Loop through the digits
  Int_t lastdet = -1;
  Int_t det     = -1;
  Int_t lastside = -1; 
  Int_t side     = -1;
  Int_t it = 0;
  Int_t ntracklets = 0;

  if (trackletContainer){ 
    for (Int_t i = 0; i < 2; i++) 
      for (Int_t j = 0; j < MAX_TRACKLETS_PERHC; j++) 
          trackletContainer[i][j] = 0; 
  }

  while (Next()) 
    {      
      det    = GetDet();
      side   = GetSide();

      if (trackletContainer)
        {
        if ((det + side*AliTRDgeometry::kNdet) != (lastdet + lastside*AliTRDgeometry::kNdet))
        {
          if (det != lastdet)
            {
              if (lastdet != -1)
                {
                return lastdet;
                }
            }
          ntracklets = GetNTracklets();
          if(ntracklets > 0) memcpy(trackletContainer[side], GetTrackletWords(), sizeof(UInt_t) * ntracklets); //copy tracklet words to trackletContainer array
          lastside = side; 
          } 
        } 

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

    if (det < 0 || det >= AliTRDgeometry::kNdet)
      {
        if (fSM.fClean == kTRUE)
    {
      AliError(Form("Strange Det Number %d BUT event buffer seems to be clean.", det));
    }
        else
    {
      AliError(Form("Strange Det Number %d. Event buffer marked NOT clean!", det));
    }
        continue;
      }

    // Add a container for the digits of this detector
    digits = (AliTRDarrayADC *) digitsManager->GetDigits(det);

          if (digitsManager->UsesDictionaries()) 
            {
        track0 = (AliTRDarrayDictionary *) digitsManager->GetDictionary(det,0);
        track1 = (AliTRDarrayDictionary *) digitsManager->GetDictionary(det,1);
        track2 = (AliTRDarrayDictionary *) digitsManager->GetDictionary(det,2);
      }

    if (!digits)
      {
        if (fSM.fClean == kTRUE)
    {
      AliError(Form("Unable to get digits for det %d BUT event buffer seems to be clean.", det));
    }
        else
    {
      AliError(Form("Unable to get digits for det %d. Event buffer is NOT clean!", det));
    }
        return -1;
      }

    Int_t rowMax = GetRowMax();
    Int_t colMax = GetColMax();
    Int_t ntbins = GetNumberOfTimeBins();

    // Allocate memory space for the digits buffer
    if (digits->GetNtime() == 0) 
      {
        digits->Allocate(rowMax, colMax, ntbins);
              if (digitsManager->UsesDictionaries()) 
                {
            track0->Allocate(rowMax, colMax, ntbins);
            track1->Allocate(rowMax, colMax, ntbins);
            track2->Allocate(rowMax, colMax, ntbins);
    }
      }

    indexes = digitsManager->GetIndexes(det);
    indexes->SetSM(GetSM());
    indexes->SetStack(GetStack());
    indexes->SetLayer(GetLayer());
    indexes->SetDetNumber(det);
    if (indexes->IsAllocated() == kFALSE)
      indexes->Allocate(rowMax, colMax, ntbins);
  }

      //Char_t padStatus =  cal->GetPadStatus(det, GetCol(), GetRow());

      // ntimebins data are ready to read
      for (it = 0; it < GetNumberOfTimeBins(); it++)
  {
    if (GetSignals()[it] > 0)
      {
        digits->SetData(GetRow(), GetCol(), it, GetSignals()[it]);
        /*if(padStatus)
	  digits->SetPadStatus(GetRow(), GetCol(), it, padStatus);*/
                
        indexes->AddIndexRC(GetRow(), GetCol());
              if (digitsManager->UsesDictionaries()) 
                {
            track0->SetData(GetRow(), GetCol(), it, 0);
            track1->SetData(GetRow(), GetCol(), it, 0);
            track2->SetData(GetRow(), GetCol(), it, 0);
    }
      }
  } // tbins
    }// while Next()

  return det;
}

//------------------------------------------------------------
Bool_t
AliTRDrawStream::Init()
{
  //
  // Initialize geometry and fee parameters 
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

  ResetCounters(); // fBufferRead is set to kFALSE - important

  saveDir->cd();

  return kTRUE;
}

//------------------------------------------------------------
Bool_t 
AliTRDrawStream::InitBuffer(void *buffer, UInt_t length)
{
  // 
  // set initial information about the buffer
  //

  if (fgDebugFlag)  AliDebug(5, Form("Equipment ID: %d",fRawReader->GetEquipmentId()));
  if (fRawReader->GetEquipmentId()<1024 || fRawReader->GetEquipmentId()>1041) //tmp protection
    return kFALSE; 

  ResetCounters();

  fpBegin = (UInt_t *)buffer;

  if (WORD_SIZE == 0)
    {
      AliFatal("Strange word size. size of UInt_t == 0");
      return kFALSE;
    }

  fWordLength = length/WORD_SIZE;
  fpEnd = fpBegin + fWordLength;
  fpPos = fpBegin;

  if (fpBegin == 0 || length <= 0)
    {
      AliError(Form("Buffer size or pointer is strange. pointer to the buffer is 0x%08x of size %d", fpBegin, length));
      return kFALSE;
    }

  SwapOnEndian();

  if (fgDumpHead >= 0){
    if ( fgDumpHead == 0 ){ // dump all words
        AliInfo(Form("---------- Dumping all words from the beginnig of the buffer ----------"));
        if (DumpWords(fpBegin, fWordLength) == kFALSE) AliError("Dump failed. Not enough data.");
    } else {
        AliInfo(Form("---------- Dumping %u words from the beginnig of the buffer ----------",fgDumpHead));
        if (DumpWords(fpBegin, fgDumpHead) == kFALSE) AliError("Dump failed. Not enough data.");
    }
    AliInfo(Form("---------- Dumping ended ----------------------------------------------"));
  }

  return kTRUE;
}

//------------------------------------------------------------
Bool_t 
AliTRDrawStream::DecodeGTUheader()
{
  // Decode Supermodule Index Word
  DecodeSMInfo(fpPos, &fSM);

  if (fgDebugFlag)  AliDebug(5, DumpSMInfo(&fSM));

  fpPos++;
  if (fpPos < fpEnd)
    {	
      // fSM.fHeaderSize represent additional Supermodule header size which contains additional information regarding hardware design.
      // For the moment, we skip decoding these words 
      if (SkipWords(fSM.fHeaderSize) == kTRUE)
  {
    for (Int_t istack = 0; istack < 5; istack++)
      {
        if (fSM.fStackActive[istack] == kFALSE)
    continue;

        fStack = &fSM.fStacks[istack];

              // Decode Stack Index Word of given stack
        DecodeStackInfo(fpPos, fStack);
        fpPos++;

        fSM.fNexpectedHalfChambers += fStack->fActiveLinks;
        
        if (fgDebugFlag)  AliDebug(5, DumpStackInfo(fStack));
        
        if (SkipWords(fStack->fHeaderSize-6) == kFALSE) // 6 is the 6 stack header words for 12 links 
    {
      if (fRawReader) fRawReader->AddMajorErrorLog(kDecodeStackInfo, "Stack header words missing");
      return kFALSE;
    }
              for (Int_t iword=0; iword<6; iword++) // decode 6 stack header words
                {
                  // Decode Stack Header Word of given stack
            DecodeStackHeader(fpPos, fStack, iword); 
            fpPos++;
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
      if (fgWarnError) AliWarning("No additional sm headers and stack index words present.");
      if (fRawReader) fRawReader->AddMajorErrorLog(kDecodeStackInfo, "Stack info missing");
      return kFALSE;
    }

  if (fpPos < fpEnd)
    {
      if (fgDebugFlag)  AliDebug(5, "GTU headers are OK.");
    }
  else
    {
      if (fgWarnError) AliWarning("No data just after GTU headers.");
      if (fRawReader) fRawReader->AddMajorErrorLog(kMissingData, "Missing sm data");
      return kFALSE;
    }

  if (fgDebugFlag)  AliDebug(5, Form("Expected half chambers from GTU header: %d", fSM.fNexpectedHalfChambers));

  return kTRUE;
}

//------------------------------------------------------------
Bool_t 
AliTRDrawStream::DecodeSM(void *buffer, UInt_t length)
{
  //
  // decode one sm data in buffer
  //
  
  ResetIterators(); 

  fSM.fClean = kTRUE;
  if (InitBuffer(buffer, length) == kFALSE)
    {
      if (fgWarnError) AliError("InitBuffer failed.");      
      fSM.fClean = kFALSE;
      return kFALSE;
    }

  if (DecodeGTUheader()== kFALSE)
    return kFALSE;

  for (Int_t istack = 0; istack < 5; istack++)
    {
      fStackNumber = istack; 
      if (fSM.fStackActive[istack] == kFALSE)
  continue;
      
      fStack = &fSM.fStacks[istack];

      fgLastHC  = -1; // to check rob number odering 
      for (Int_t ilink = 0; ilink < 12; ilink++)
  {
    fStackLinkNumber = ilink; 
    if (fStack->fLinksActive[ilink] == kFALSE)
      continue;

          // check GTU link monitor 
          if (!(fStack->fLinksDataType[ilink] == 0 && fStack->fLinksMonitor[ilink] == 0))
            {
              fStack->fLinkMonitorError[ilink] = 1;
              SeekEndOfData(); // skip this HC data if GTU link monitor report error
              fStack->fLinkMonitorError[ilink] += fNWordsCounter; // counts words of given hc having link monitor error
              continue; 
            }

    if (fpPos >= fpEnd)
      {
        if (fRawReader) fRawReader->AddMajorErrorLog(kLinkDataMissing, "Link data missing");	      
              if (fgWarnError) AliError("Link data missing.");      
        fSM.fClean = kFALSE;
        break;
      }

    fHC = &fStack->fHalfChambers[ilink];
          ResetPerHC();

    if (fSM.fTrackletEnable == kTRUE)
      {
        if (DecodeTracklets() == kFALSE)
    {
      
      fSM.fClean = kFALSE;
      SeekEndOfData();

      if (fgWarnError) 
        {
          AliError(Form("Tracklet decoding failed stack %d link %d", fStackNumber, fStackLinkNumber));
        }
      continue;
    }
      }

    if (fpPos >= fpEnd)
      {
        if (fRawReader) fRawReader->AddMajorErrorLog(kHCdataMissing, "HC data missing");	      
              if (fgWarnError) AliError("HC data missing.");      
        fSM.fClean = kFALSE;
        break;
      }
    
          fgLastROB   = -1; // to check mcm number odering 
          fgLastIndex = -1 ; // to check mcm number odering 
    if (DecodeHC() == kFALSE)
      {
        fSM.fClean = kFALSE;
              if (fHC->fCorrupted < 16)  SeekEndOfData(); // In case that we meet END_OF_TRACKLET_MARKERNEW 
                                                          // during ADC data decoding or MCM header decoding
                                                          // we don't seek ENDOFRAWDATAMARKER

        if (fgWarnError) 
    {
      AliError(Form("Failed HC : %s", DumpHCinfoH0(fHC)));
      AliError(Form("Failed HC : %s", DumpHCinfoH1(fHC)));
    }
                
        continue;
      }
    else
      {
        SeekEndOfData(); // make sure that finish off with the end of data markers
      }

  } // ilink
    } // istack

  ResetIterators(); // need to do it again for Next() function 

  if (fSM.fClean == kTRUE)
    return kTRUE;
  
  if (fgCleanDataOnly && (fSM.fClean == kFALSE))
    {
      if (fgWarnError) 
  {
    AliWarning("Buffer with errors. Returning FALSE.");
    AliWarning(Form("--- Failed SM : %s ---", DumpSMInfo(&fSM)));
  }
      fSM.fActiveStacks = 0; // Next() will not give data
      return kFALSE;
    }

  return kTRUE;
}

//------------------------------------------------------------
Int_t 
AliTRDrawStream::DecodeSM()
{
  //
  // decode SM data in case AliRawReader is in use
  //    
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
AliTRDrawStream::DecodeSM(AliRawReader *reader)
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
AliTRDrawStream::SeekEndOfData()
{
  //
  // go to end of data marker
  //
  Int_t fEndOfDataCount = 0;
  fNWordsCounter = 0;

  while ( *fpPos != ENDOFRAWDATAMARKER && fpPos < fpEnd )
    {
      fpPos++;
      fNWordsCounter++;
    }
  while (*fpPos == ENDOFRAWDATAMARKER && fpPos < fpEnd )
    {
      fEndOfDataCount++;
      fpPos++;      
    }
  
  return kTRUE;
}

//------------------------------------------------------------
Bool_t
AliTRDrawStream::SkipMCMdata(UInt_t iw)
{

  if (fgDebugFlag) AliDebug(11,Form("Skip %d words due to MCM header corruption.",iw));
  UInt_t iwcounter = 0;  
  while ( *fpPos != ENDOFRAWDATAMARKER && iwcounter < iw)
    {
      if ( *fpPos == END_OF_TRACKLET_MARKERNEW) 
        {  
          if (fgDebugFlag) AliDebug(11,"Met END_OF_TRACKLET_MARKERNEW");
          fMCM->fCorrupted += 16;
          fHC->fCorrupted += 16;
          return kFALSE;
        } 
      fpPos++;
      iwcounter++; 
    }

  if (iwcounter == iw)
    {
      fpPos++;
      return kTRUE;
    }

  if (fgDebugFlag) AliDebug(11,"Met ENDOFRAWDATAMARKER");
  return kFALSE;
}

//------------------------------------------------------------
Bool_t 
AliTRDrawStream::SeekNextMCMheader()
{
  //
  // go to mcm marker
  //

  fpPos++;

  while ( *fpPos != ENDOFRAWDATAMARKER && fpPos < fpEnd )
    {
      if (MCM_HEADER_MASK_ERR(*fpPos) == 0 && MCM_HEADER_MASK_ERR(*(fpPos+1)) == 0)      
  {
    if (fgDebugFlag) AliDebug(11,Form("^^^ Found : Pos 0x%08x : Val 0x%08x", fpPos, *fpPos));
    return kTRUE;
  }
      if ( *fpPos == END_OF_TRACKLET_MARKERNEW) 
        {  
          fMCM->fCorrupted += 16;
          fHC->fCorrupted += 16;
          return kFALSE;
        } 
      fpPos++;
    }

  SeekEndOfData();
  return kFALSE;
}

//------------------------------------------------------------

Bool_t 
AliTRDrawStream::DecodeTracklets()
{
  //
  // decode tracklets
  //

  fLinkTrackletCounter = 0;
  fEndOfTrackletCount = 0;
  fHC->fNTracklets = 0;

  for (Int_t i = 0; i < MAX_TRACKLETS_PERHC; i++) 
  fHC->fTrackletWords[i] = 0; 

  if (fgDebugFlag)  AliDebug(10, Form("Decode tracklets at 0x%08x : 0x%08x", fpPos, *fpPos));

  while ( *fpPos != END_OF_TRACKLET_MARKEROLD && *fpPos != END_OF_TRACKLET_MARKERNEW && fpPos < fpEnd )
    {
      if (fgDebugFlag)  AliDebug(10, Form("Tracklet found at 0x%08x : 0x%08x", fpPos, *fpPos));

      fLinkTrackletCounter++;

      if (fLinkTrackletCounter > MAX_TRACKLETS_PERHC)
  {
    if (fgDebugFlag) AliDebug(11,Form("Max number of tracklets exceeded %d > %d.", 
        fLinkTrackletCounter, MAX_TRACKLETS_PERHC));
    if (fRawReader) fRawReader->AddMajorErrorLog(kTrackletOverflow,"Too many tracklets"); 
          fHC->fTrackletError = 1;
    return kFALSE;
  }

      fHC->fTrackletWords[fLinkTrackletCounter-1] = UInt_t(*fpPos); //store tracklet words into array  
      fHC->fNTracklets = fLinkTrackletCounter;
      fpPos++;
    }

  while ( ( *fpPos == END_OF_TRACKLET_MARKEROLD || *fpPos == END_OF_TRACKLET_MARKERNEW ) && fpPos < fpEnd )
    {
      if (fgDebugFlag)  AliDebug(10, Form("EoTracklets found at 0x%08x : 0x%08x", fpPos, *fpPos));

      fEndOfTrackletCount++;
      fpPos++;
    }

  if ( fEndOfTrackletCount < 2 )
    {
      if (fgDebugFlag) AliDebug(11,"End of tracklets word missing"); 
      if (fRawReader) fRawReader->AddMajorErrorLog(kEOTrackeltsMissing, "End of tracklets word missing"); 
      fHC->fTrackletError += 2;
      return kFALSE;
    }

  return kTRUE;
}

//------------------------------------------------------------
Bool_t 
AliTRDrawStream::IsRowValid()
{
  if ( (fHC->fStack == 2 && fMCM->fROW >= fGeometry->RowmaxC0()) ||
      (fHC->fStack != 2 && fMCM->fROW >= fGeometry->RowmaxC1()) || fMCM->fROW < 0 ) 
    {
      if (fgDebugFlag) AliDebug(11,Form("SM%d L%dS%d: Wrong Padrow (%d) fROB=%d, fSIDE=%d, fMCM=%02d"
          , fHC->fSM, fHC->fLayer, fHC->fStack, fMCM->fROW, fMCM->fROB, fHC->fSide, fMCM->fMCM ));
      return kFALSE;
    }
  return kTRUE;
}

//------------------------------------------------------------
Bool_t 
AliTRDrawStream::IsMCMheaderOK()
{
  //
  // check the mcm header
  //

  if (fgLastROB != fMCM->fROB) 
    {
      fgLastIndex = 0;
      if (fgLastROB== -1) fgLastROB = fMCM->fROB;
    }
  else
    {
      Int_t matchingcounter = 0; 
      for (Int_t i=fgLastIndex+1; i<16; i++)
        { 
          if ( fMCM->fMCM == fgMCMordering[i] )
            {
              fgLastIndex = i;
              matchingcounter++;
              break;
            }
        }
      if (matchingcounter == 0)   
        {
          fMCM->fMCMhdCorrupted += 2;
          AliDebug(11,Form("MCM number from last MCM is larger: MCM # from last MCM %d,  MCM # from current MCM %d \n",(fMCM-1)->fMCM, fMCM->fMCM));
        }
    }

  if ( fgLastHC == fHC->fLayer*2 + fHC->fSide )
    {
      if ( fMCM->fROB < (fMCM-1)->fROB )
        {
        fMCM->fMCMhdCorrupted += 2;
        AliDebug(11,Form("ROB number from last MCM is larger: ROB # from last MCM %d,  ROB # from current MCM %d \n",(fMCM-1)->fROB, fMCM->fROB));
        }
      else fgLastROB = fMCM->fROB; 
    }

  fgLastHC = fHC->fLayer*2 + fHC->fSide; 

  /*
  // this check will come back later again when we have "patched MCM map"
  int expectedROB = -1;
  if(!fHC->fSide) expectedROB = int(fHC->fMCMmax/16)*2;
  else expectedROB = int(fHC->fMCMmax/16)*2 + 1;
  int expectedMCM = 4*(3-int((fHC->fMCMmax%16)/4)) + fHC->fMCMmax%4;

  if ( expectedROB != fMCM->fROB || expectedMCM != fMCM->fMCM)
    {
      fMCM->fMCMhdCorrupted += 2;
      AliDebug(11,Form("ROB expected %d ROB read %d,  MCM expected %d MCM read %d\n",expectedROB, fMCM->fROB, expectedMCM, fMCM->fMCM));
    }
  */

  // below two conditions are redundant  
  /*
  if ( fMCM->fMCM < 0 || fMCM->fMCM > 15 || fMCM->fROB < 0 || fMCM->fROB > 7 ) 
    {
      fMCM->fMCMhdCorrupted += 8;  // need to assign new number
      if (fgDebugFlag) AliDebug(11,Form("ROB or MCM number is out of range. %s\n", DumpMCMinfo(fMCM)));
    }
  if (IsRowValid() == kFALSE)
    {
      fMCM->fMCMhdCorrupted += 16; // need to assign new number
    }
  */  
    
  if (fEventCounter == 0)
    {
      fEventCounter = fMCM->fEvCounter;
    }

  if (fEventCounter != fMCM->fEvCounter)
    {
      fMCM->fMCMhdCorrupted += 4;      if (fgDebugFlag) AliDebug(11,Form("Event number(%d) of current MCM is different from that(%d) of reference MCM %s.\n", fMCM->fEvCounter, fEventCounter, DumpMCMinfo(fMCM)));
    }

  if (fEventCounter < fLastEventCounter)
    {
      fMCM->fMCMhdCorrupted += 8;      if (fgDebugFlag) AliDebug(11,Form("Event from the past? Current %d Last %d %s.\n", fEventCounter, fLastEventCounter, DumpMCMinfo(fMCM)));
    }

  if ( fMCM->fADCmaskCorrupted > 0 )
      return kFALSE;

  if ( fMCM->fMCMhdCorrupted > 0 )
      return kFALSE;

  return kTRUE;
}

//------------------------------------------------------------
Bool_t 
AliTRDrawStream::DecodeMCMheader()
{
  //
  // decode the mcm header
  //

  DecodeMCMheader(fpPos, fMCM); 

  if (fDumpingEnable) 
    {
      if (fMCM->fMCM == fDumpingMCM) 
        {
          if (fMCM->fROB == fDumpingROB && fHC->fLayer == fDumpingLayer)
            {
              if (fHC->fSM == fDumpingSM && fHC->fStack == fDumpingStack)
                { 
                  if (fgDebugFlag) {
                    AliDebug(5,DumpHCinfoH0(fHC));
                    AliDebug(5,DumpMCMinfo(fMCM));
                  }
                  DumpWords(fpPos, 212);
                }  
            }
        }
    }

  if (fHC->fCorrupted >= 16)
    {
      fpPos--; 
      return kFALSE;
    }

  fMCM->fROW = fTRDfeeParam->GetPadRowFromMCM(fMCM->fROB, fMCM->fMCM); 

  if ((fHC->fRawVMajor > 2 && fHC->fRawVMajor <5) || ((fHC->fRawVMajor & 32) == 32)) //cover old and new version definition of ZS data
    {
      fpPos++;
      if ( fpPos < fpEnd )
  {
    DecodeMask(fpPos, fMCM); 
          if (fHC->fCorrupted >= 16)
            {
              fpPos--; 
              return kFALSE;
            }
    MCMADCwordsWithTbins(fHC->fTimeBins, fMCM);
    fMCM->fAdcDataPos = fpPos + 1;
  }
      else
  {
    if (fgDebugFlag) AliDebug(11,"Expected ADC mask word. Fail due to buffer END.");	  
    if (fRawReader) fRawReader->AddMajorErrorLog(kMCMADCMaskMissing,"Missing"); 
          fHC->fCorrupted += 32;
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
      return kFALSE;
    
  return kTRUE;
}

//------------------------------------------------------------
Bool_t 
AliTRDrawStream::IsHCheaderOK()
{
  //
  // check insanity of half chamber header
  //

  if (fHC->fStack < 0 || fHC->fStack > 4)
    {
      if (fgDebugFlag) AliDebug(11,Form("Wrong Stack %d", fHC->fStack));
      return kFALSE;
    }

  if (fHC->fLayer < 0 || fHC->fLayer >= AliTRDgeometry::kNlayer)
    {
      if (fgDebugFlag) AliDebug(11,Form("Wrong layer %d", fHC->fLayer));
      return kFALSE;
    }

  if (fHC->fSide < 0 || fHC->fSide > 1)
    {
      if (fgDebugFlag) AliDebug(11,Form("Wrong Side %d", fHC->fSide));
      return kFALSE;
    }

  if (fgStackNumberChecker)
    {
    if (fHC->fStack != fStackNumber) 
      {
        if (fgDebugFlag) AliDebug(11,Form("Missmatch: Stack number between HC header %d and GTU link mask %d", 
              fHC->fStack, fStackNumber));
        fStackNumber = -1;
        return kFALSE;
    }
    }

  if (fgStackLinkNumberChecker)
    {
    //if (fHC->fLayer * 2 + fHC->fSide != fStackLinkNumber) 
    // let it make flexible to consider known fiber swapping
    if ((fHC->fLayer * 2 != fStackLinkNumber) && (fHC->fLayer * 2 != fStackLinkNumber - 1))  
      {
        if (fgDebugFlag) AliDebug(11,Form("Missmatch: Layer number between HC header %d and GTU link mask %d | %s", 
                    fHC->fLayer, fStackLinkNumber, DumpStackInfo(fStack)));
        fStackLinkNumber = -1;
        return kFALSE;      
      }
    }

  // SLOW GEOM : consistancy check with geometry
  fHC->fDET = fGeometry->GetDetector(fHC->fLayer, fHC->fStack, fHC->fSM);
  if (fHC->fDET < 0 || fHC->fDET >= AliTRDgeometry::kNdet)
    {
      if (fgDebugFlag) AliDebug(11,Form("Wrong detector %d", fHC->fDET));      
      if (fRawReader) fRawReader->AddMajorErrorLog(kHCHeaderWrongDet, "Wrong Det");       
      return kFALSE;
    }

  if (fHC->fSM != fGeometry->GetSector(fHC->fDET)
      || fHC->fSM <0 || fHC->fSM >= AliTRDgeometry::kNsector)
    {
      if (fgDebugFlag) AliDebug(11,Form("Wrong SM(sector) %d (Geometry says: %d) Stack=%d Layer=%d Det=%d", 
              fHC->fSM, fGeometry->GetSector(fHC->fDET),
              fHC->fStack, fHC->fLayer, fHC->fDET));      
      if (fRawReader) fRawReader->AddMajorErrorLog(kHCHeaderWrongSM, "Wrong SM");       
      return kFALSE;
    }

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
  
  return kTRUE;
}

//------------------------------------------------------------
Bool_t 
AliTRDrawStream::DecodeHCheader()
{  
  //
  // decode the half chamber header
  //

  if (DecodeHCwordH0(fpPos, fHC) == kFALSE)
    return kFALSE;
    
  if (fHC->fNExtraWords > 0)
    {
      fpPos++;
      if (fpPos < fpEnd)
  {
    if (DecodeHCwordH1(fpPos, fHC) == kFALSE)
            return kFALSE;
  }
      else
  {
    if (fgDebugFlag) AliDebug(11,"Expected HC header word 1. Fail due to buffer END.");
    if (fRawReader) fRawReader->AddMajorErrorLog(kHCWordMissing,"Next HC word 1 (count from 0) missing"); 
    return kFALSE;
  }
    }

  if (fgDebugFlag)  AliDebug(5, DumpHCinfoH0(fHC));
  if (fgDebugFlag)  AliDebug(5, DumpHCinfoH1(fHC));

  fHC->fDET = -1;
  if (IsHCheaderOK() == kFALSE)
    {
      fHC->fH0Corrupted += 2;
      if (fgDebugFlag) AliDebug(11,Form("H0 Header Insane. Word 0x%08x", *fHC->fPos));
      return kFALSE;
    }
  
  return kTRUE;
}

//------------------------------------------------------------
Bool_t 
AliTRDrawStream::DecodeHC()
{
  //
  // decode hc header and data
  //

  if (DecodeHCheader() == kFALSE)
    {
      if (fgWarnError) AliWarning(Form("HC Header decode failed. H0 Error: %d H1 Error: %d",fHC->fH0Corrupted,fHC->fH1Corrupted));
      if (fRawReader) fRawReader->AddMajorErrorLog(kHCHeaderCorrupt, "HC header corrupted"); 
      return kFALSE;
    }
  else
    {
      fpPos++;
      if (fpPos >= fpEnd)
  {
          fHC->fCorrupted += 1;
    if (fgDebugFlag) AliDebug(11,"No MCM data? Not enough data in the buffer.");
    if (fRawReader) fRawReader->AddMajorErrorLog(kMCMdataMissing, "MCM data missing"); 
    return kFALSE;
  }
    }

  if ((fHC->fRawVMajor & 64) == 64) // test pattern data
    {
      AliTRDrawTPStream *tpStream = new AliTRDrawTPStream(fHC->fRawVMajorOpt, fpPos);
      if (tpStream->DecodeTPdata() == kFALSE)
        {
         if (fgWarnError) AliError("failed to decode test pattern data");
         return kFALSE; 
        }
      return kTRUE;
    } 

  fHC->fMCMmax = 0;
  while (*fpPos != ENDOFRAWDATAMARKER && fpPos < fpEnd)
    {
      if (fHC->fMCMmax > TRD_MAX_MCM)
  {
          fHC->fCorrupted += 2;
    if (fgDebugFlag) AliDebug(11,"More mcm data than expected!");
    if (fRawReader) fRawReader->AddMajorErrorLog(kMCMoverflow, "Too many mcms found!"); 
    return kFALSE;
  }

      fMCM = &fHC->fMCMs[fHC->fMCMmax];

      if (DecodeMCMheader() == kFALSE)
  {
          if (fHC->fCorrupted < 4) fHC->fCorrupted += 4; // benchmark hc data corruption as 4
    
    if (fgSkipData == kTRUE || fHC->fCorrupted >= 16)
              return kFALSE; // stop HC data reading
          
          fHC->fMCMmax++; // increase mcm counter to match with expected rob/mcm number

          // in case we decide to keep reading data, skip this mcm data and find next mcm header 
          if (fMCM->fADCmaskCorrupted < 2) 
            {  
              if (SkipMCMdata(fMCM->fADCcount*fMCM->fSingleADCwords) == kFALSE)
            return kFALSE;
              continue;
            }
          else 
            {
              if (SeekNextMCMheader() == kFALSE)
            return kFALSE;
              continue;
            }
  }

      fHC->fMCMmax++;

      if (fMCM->fADCmax > 0)
  {
    fpPos++;
    if (fpPos >= fpEnd)
      {
              fMCM->fCorrupted += 1;
              if (fHC->fCorrupted < 4) fHC->fCorrupted += 4; // benchmark hc data corruption as 4
        if (fgDebugFlag)  AliDebug(9, Form("Buffer short of data. ADC data expected."));	  
        return kFALSE;
      }

    for (Int_t iadc = 0; iadc < fMCM->fADCmax; iadc++)
      {
        fADC = &fMCM->fADCs[iadc];
        fADC->fADCnumber = fMCM->fADCchannel[iadc];

        if (fgDebugFlag)  AliDebug(9, Form("This is ADC %d of %d. ADC number is %d.", 
            iadc+1, fMCM->fADCmax, fMCM->fADCchannel[iadc]));

        if (fpPos + fMCM->fSingleADCwords >= fpEnd)
    {
      
                  fMCM->fCorrupted += 2;
                  if (fHC->fCorrupted < 4) fHC->fCorrupted += 4; // benchmark hc data corruption as 4
      if (fgDebugFlag) AliDebug(11,"ADC (10 words) expected. Not enough data in the buffer.");
      return kFALSE;
    }

              if (fHC->fRawVMajor < 64) // normal(real) ADC data
                {
            if (DecodeADC() == kFALSE)
        {
                      if (fMCM->fCorrupted < 4) fMCM->fCorrupted += 4; // benchmark mcm data corruption as 4
                      if (fHC->fCorrupted < 4) fHC->fCorrupted += 4;   // benchmark hc data corruption as 4
          if (fADC->fIsShared && fADC->fCorrupted == 16)   // check if we are out of the det when the pad is shared
            {
              fADC->fCOL = -1;
              fpPos = fADC->fPos + fMCM->fSingleADCwords;
            }
          else
            {
              if (fgDebugFlag) AliDebug(11,Form("ADC decode failed."));
                    if (fgSkipData == kTRUE || fHC->fCorrupted >= 16) 
                            return kFALSE; // stop HC data reading
            }
        }
                } 
              else // test pattern data
                {
                  if (fgWarnError) AliError("These are test pattern data. You need other reader"); // will be served in other class
                }
      } 
  } 
      else
  {
    fpPos++;
  }
    }//while eof data

  if (fpPos >= fpEnd)
    {
      if (fgDebugFlag) AliDebug(11,"We are at the end of buffer. There should be one more word left.");
      return kFALSE;
    }

  return kTRUE;
}
//------------------------------------------------------------

Bool_t
AliTRDrawStream::DecodeADC()
{
  //
  // decode single ADC channel
  //

  fADC->fCorrupted = 0;
  if(fADC->fADCnumber%2==1) fMaskADCword = ADC_WORD_MASK(ADCDATA_VAL1);
  if(fADC->fADCnumber%2==0) fMaskADCword = ADC_WORD_MASK(ADCDATA_VAL2);

  fADC->fPos = fpPos;
  fTbinADC = 0;

  for (Int_t i = 0; i < TRD_MAX_TBINS; i++)
    fADC->fSignals[i] = 0;

  for (Int_t iw = 0; iw < fMCM->fSingleADCwords; iw++)
    {
      if (HC_HEADER_MASK_ERR(*fpPos) == 0 || *fpPos == END_OF_TRACKLET_MARKERNEW)
        {
          if (fgWarnError) AliError(Form("There should be ADC data. We meet HC header or END_OF_TRACKLET_MARKER 0x%08x",*fpPos));
    fADC->fCorrupted += 16;
          fHC->fCorrupted += 16; 
          fpPos--;

          return kFALSE;
        }
      if (fMaskADCword != ADC_WORD_MASK(*fpPos))
  {
    fADC->fCorrupted += 1;
          if (fgDebugFlag) AliDebug(11,Form("Wrong ADC data mask! ADC channel number: %02d [Expected mask: 0x%08x  Current mask: 0x%08x] MCM= %s Error : %d",
                                          fADC->fADCnumber, fMaskADCword, ADC_WORD_MASK(*fpPos),DumpMCMinfo(fMCM),fADC->fCorrupted));
          fpPos++;
    continue;
  }

      // here we subtract the baseline ( == common additive)
      fADC->fSignals[fTbinADC + 0] = ((*fpPos & 0x00000ffc) >>  2) - fgCommonAdditive;
      fADC->fSignals[fTbinADC + 1] = ((*fpPos & 0x003ff000) >> 12) - fgCommonAdditive;
      fADC->fSignals[fTbinADC + 2] = ((*fpPos & 0xffc00000) >> 22) - fgCommonAdditive;

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
    fADC->fCorrupted += 32;
    if (fgDebugFlag) AliDebug(11,Form("Wrong column! ADCnumber %d MaxIs %d Col %d MaxIs %d MCM= %s", 
            fADC->fADCnumber, fMaxADCgeom, fADC->fCOL, fHC->fColMax, DumpMCMinfo(fMCM)));
  }
      //else
  //{
    // we are out of the det when the pad is shared
    //if (fgDebugFlag) AliDebug(11, Form("Column out of the detector! ADCnumber %d MaxIs %d Col %d MaxIs %d MCM= %s", 
    //				     fADC->fADCnumber, fMaxADCgeom, fADC->fCOL, fHC->fColMax, DumpMCMinfo(fMCM)));
    //fADC->fCorrupted += 32;
  //}
    }

  if (fADC->fCorrupted > 0)
    {
      return kFALSE;
    }

  fDecodedADCs++;
  return kTRUE;
}

//--------------------------------------------------------


void AliTRDrawStream::DecodeSMInfo(const UInt_t *word, struct AliTRDrawSM *sm) const
{
  //
  // Decode Supermodule Index Word
  // The Supermodule Index Word is a 32-Bit word wit following structure
  // ssssssss ssssssss vvvv rrrr r d t mmmm
  // s: Size of the Supermodule Header, v: Supermodule Header Version, r: Reserved for future use
  // d: Track Data Enabled Bit, t: Tracklet Data Enabled Bit, m: Stack Mask 
  //
  sm->fPos = (UInt_t*)word; 

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
const char *AliTRDrawStream::DumpSMInfo(const struct AliTRDrawSM *sm)
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
void AliTRDrawStream::DecodeStackInfo(const UInt_t *word, struct AliTRDrawStack *st) const
{
  //
  // Decode Stack #i Index Word
  // The Stack #i Index Word is a 32-Bit word wit following structure
  // ssssssss ssssssss vvvv mmmm mmmmmmmm
  // s: Size of the Stack #i Header, v: Supermodule Header Version, m: Link Mask
  //
  st->fPos = (UInt_t*)word;
      
  UInt_t vword = *word;
  st->fHeaderSize = STACK_HEADER_SIZE(vword);

  UInt_t linkMask = STACK_LINK_WORD(vword);
  st->fActiveLinks = 0;
  for (Int_t i = 0; i < 12; i++)
    {
      if (IS_BIT_SET(linkMask,i) > 0)
  {
    st->fLinksActive[i] = kTRUE;
    st->fActiveLinks++;
  }
      else
  {
    st->fLinksActive[i] = kFALSE;
  }
    }
}
  
//--------------------------------------------------------
void AliTRDrawStream::DecodeStackHeader(const UInt_t *word, struct AliTRDrawStack *st, Int_t iword) const
{
      st->fPos = (UInt_t*)word;
      
      UInt_t vword = *word;
      st->fLinksDataType[2*iword]    = LINK0_DATA_TYPE_FLAG(vword);
      st->fLinksMonitor[2*iword]     = LINK0_MONITOR_FLAG(vword);
      st->fLinksDataType[2*iword+1]  = LINK1_DATA_TYPE_FLAG(vword);
      st->fLinksMonitor[2*iword+1]   = LINK1_MONITOR_FLAG(vword);
}

//--------------------------------------------------------
const char *AliTRDrawStream::DumpStackInfo(const struct AliTRDrawStack *st)
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
Bool_t AliTRDrawStream::DecodeHCwordH0(const UInt_t *word, struct AliTRDrawHC *hc) const
{
  //
  // decode the hc header word 0
  //
  UInt_t vword = *word;

  hc->fH0Corrupted = HC_HEADER_MASK_ERR(vword);
  if (hc->fH0Corrupted > 0)
    {
    if (fgDebugFlag) AliDebug(11,Form("H0 Header Mask Error. Word 0x%08x", *fHC->fPos));
    return kFALSE;
    }

  hc->fSpecialRawV =  HC_SPECIAL_RAW_VERSION(vword);
  hc->fRawVMajor = HC_MAJOR_RAW_VERSION(vword);
  hc->fRawVMajorOpt = HC_MAJOR_RAW_VERSION_OPT(vword); 
  hc->fRawVMinor = HC_MINOR_RAW_VERSION(vword);
  hc->fNExtraWords = HC_EXTRA_WORDS(vword);
  hc->fDCSboard = HC_DCS_BOARD(vword);
  hc->fSM = HC_SM_NUMBER(vword);
  hc->fStack = HC_STACK_NUMBER(vword);
  hc->fLayer = HC_LAYER_NUMBER(vword);
  hc->fSide = HC_SIDE_NUMBER(vword);

  hc->fPos[0] = (UInt_t*)word;

  return kTRUE;
}

//--------------------------------------------------------
Bool_t AliTRDrawStream::DecodeHCwordH1(const UInt_t *word, struct AliTRDrawHC *hc) const
{
  //
  // decode the hc header word 1
  //

  UInt_t vword = *word;

  hc->fH1Corrupted = HC_HEADER_MASK_ERR(vword);
  if (hc->fH1Corrupted > 0)
    { 
    if (fgDebugFlag) AliDebug(11,Form("H1 Header Mask Error. Word 0x%08x", *fHC->fPos));
    return kFALSE;
    }

  hc->fTimeBins = HC_NTIMEBINS(vword);
  hc->fBunchCrossCounter = HC_BUNCH_CROSS_COUNTER(vword);
  hc->fPreTriggerCounter = HC_PRETRIGGER_COUNTER(vword);
  hc->fPreTriggerPhase = HC_PRETRIGGER_PHASE(vword);

  hc->fPos[1] = (UInt_t*)word;

  return kTRUE;
}
  
//--------------------------------------------------------
const char *AliTRDrawStream::DumpHCinfoH0(const struct AliTRDrawHC *hc)
{
  //
  // dump the hc header word 0
  //
  if (!hc)
    return Form("Unable to dump. Null received as parameter!?!");
  else
    return Form("[ HC[0] at 0x%08x ] : 0x%08x Info is : RawV %d SM %d Stack %d Layer %d Side %d DCSboard %d",
    hc->fPos[0], *(hc->fPos[0]), hc->fRawVMajor, hc->fSM, hc->fStack, hc->fLayer, hc->fSide, hc->fDCSboard);
}

//--------------------------------------------------------
const char *AliTRDrawStream::DumpHCinfoH1(const struct AliTRDrawHC *hc)
{
  //
  // dump the hc header word 1
  //
  if (!hc)
    return Form("Unable to dump. Null received as parameter!?!");
  else
    return Form("[ HC[1] at 0x%08x ] : 0x%08x Info is : TBins %d BCcount %d PreTrigCount %d PreTrigPhase %d",
    hc->fPos[1], *(hc->fPos[1]), hc->fTimeBins, hc->fBunchCrossCounter, hc->fPreTriggerCounter, hc->fPreTriggerPhase);
}

//--------------------------------------------------------
void AliTRDrawStream::DecodeMCMheader(const UInt_t *word, struct AliTRDrawMCM *mcm) const
{
  //
  // decode the mcm header
  //
  UInt_t vword = *word;

  if (vword == END_OF_TRACKLET_MARKERNEW) 
    {
      if (fgWarnError) AliError(Form("There should be MCM header. We meet END_OF_TRACKLET_MARKER 0x%08x",vword));
      mcm->fMCMhdCorrupted += 16;
      fHC->fCorrupted += 16; //to finish data reading of this HC
    }

  mcm->fMCMhdCorrupted = MCM_HEADER_MASK_ERR(vword); //if MCM header mask has error
  if (fgDebugFlag && mcm->fMCMhdCorrupted != 0) AliDebug(11,Form("Wrong MCM header mask 0x%08x.\n", *fpPos));

  mcm->fROB = MCM_ROB_NUMBER(vword);
  mcm->fMCM = MCM_MCM_NUMBER(vword);
  mcm->fEvCounter = MCM_EVENT_COUNTER(vword);
  mcm->fPos = (UInt_t*)word;
}

//--------------------------------------------------------
UInt_t AliTRDrawStream::GetMCMadcMask(const UInt_t *word, struct AliTRDrawMCM *mcm) const
{
  //
  // get the adc mask
  //
  UInt_t vword = *word;

  mcm->fADCmax    = 0;
  mcm->fADCMask   = 0;
  mcm->fADCcount  = 0;
  mcm->fADCMaskWord = vword;

  if (vword == END_OF_TRACKLET_MARKERNEW)
    {
      if (fgWarnError) AliError(Form("There should be MCMadcMask. We meet END_OF_TRACKLET_MARKER 0x%08x",vword));
      mcm->fADCmaskCorrupted += 16;
      fHC->fCorrupted += 16; //to finish data reading of this HC
    }

  if ( MCM_ADCMASK_MASK_ERR(vword) == 0 )
    {
      mcm->fADCMask  = MCM_ADCMASK_VAL(vword);
      mcm->fADCcount = MCM_ADCMASK_NADC(~vword);
    }
  else
    {
      mcm->fADCMask = 0xffffffff;
      mcm->fADCmaskCorrupted = 1; // mcm adc mask error
      if (fgDebugFlag) AliDebug(11,Form("Wrong ADC Mask word 0x%08x.\n", *fpPos));
    }

  return mcm->fADCMask;
}

//--------------------------------------------------------
void AliTRDrawStream::DecodeMask(const UInt_t *word, struct AliTRDrawMCM *mcm) const
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
  if (mcm->fADCcount != mcm->fADCmax && fHC->fRawVMajor >= 32) // backward compatibility
    {
      mcm->fADCmaskCorrupted += 2; 
      if (fgDebugFlag) AliDebug(11,Form("ADC counts from ADCMask are different %d %d : ADCMask word 0x%08x\n", mcm->fADCcount, mcm->fADCmax, *fMCM->fPos));
    }
}

//--------------------------------------------------------
void AliTRDrawStream::MCMADCwordsWithTbins(UInt_t fTbins, struct AliTRDrawMCM *mcm) const
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
const char *AliTRDrawStream::DumpMCMinfo(const struct AliTRDrawMCM *mcm)
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
const char *AliTRDrawStream::DumpMCMadcMask(const struct AliTRDrawMCM *mcm)
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
