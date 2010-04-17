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

/* $Id: AliTRDrawFastStream.cxx 27797 2008-08-05 14:37:22Z cblume $ */

////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// This class provides access to TRD digits in raw data in a way of streaming //
//                                                                            //
// It loops over all TRD digits in the raw data given by the AliRawReader.    //
// The Next method goes to the next digit. If there are no digits left        //
// it returns kFALSE.                                                         //
// Several getters provide information about the current digit.               //
//                                                                            //
// Author: M. Ploskon (ploskon@ikf.uni-frankfurt.de)                          //  
// Author: MinJung Kweon (minjung@physi.uni-heidelberg.de)                    //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#include "TFile.h"

#include "AliTRDrawFastStream.h"
#include "AliTRDgeometry.h"
#include "AliTRDfeeParam.h"
#include "AliTRDdigitsManager.h"
#include "AliTRDarrayDictionary.h"
#include "AliTRDarrayADC.h"
#include "AliTRDSignalIndex.h"
#include "AliTRDdigitsParam.h"

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
#define HC_PRETRIGGER_PHASE(w) GET_VALUE_AT(w,0xf,2)

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
ClassImp(AliTRDrawFastStream)

Bool_t AliTRDrawFastStream::fgExtraSkip = kFALSE;
Bool_t AliTRDrawFastStream::fgSkipCDH = kFALSE;
Bool_t AliTRDrawFastStream::fgCleanDataOnly = kFALSE;
Bool_t AliTRDrawFastStream::fgDebugFlag = kTRUE;
Bool_t AliTRDrawFastStream::fgEnableMemoryReset = kTRUE;
Bool_t AliTRDrawFastStream::fgStackNumberChecker = kTRUE;
Bool_t AliTRDrawFastStream::fgStackLinkNumberChecker = kFALSE;
Bool_t AliTRDrawFastStream::fgSkipData = kTRUE;
Bool_t AliTRDrawFastStream::fgEnableDecodeConfigData = kFALSE;
Int_t AliTRDrawFastStream::fgDumpHead = -1;
Short_t AliTRDrawFastStream::fgMCMordering[] =
  {
    12, 13, 14, 15, 8, 9, 10, 11, 4, 5, 6, 7, 0, 1, 2, 3 
  };
Short_t AliTRDrawFastStream::fgROBordering[] =
  {
    0, 1, 2, 3
  };

Int_t  AliTRDrawFastStream::fgLastHC = -1;
Int_t  AliTRDrawFastStream::fgLastROB = -1;
Int_t  AliTRDrawFastStream::fgLastIndex = -1;

//--------------------------------------------------------
AliTRDrawFastStream::AliTRDrawFastStream()
  : AliTRDrawStreamBase()
  , fSM()
  , fStack(0)
  , fHC(0)
  , fLastHC(0)
  , fMCM()
  , fpPos(0)
  , fpBegin(0)
  , fpEnd(0)
  , fWordLength(0)
  , fpPosTemp(0)
  , fGlobalNTimeBins(0)
  , fIsTimeBinSet(kFALSE)
  , fStackNumber(-1)
  , fStackLinkNumber(-1)
  , fLinkTrackletCounter(-1)
  , fEndOfTrackletCount(-1)
  , fNWordsCounter(-1)
  , fMaskADCword(0)
  , fTbinADC(0)
  , fEventCounter(0)
  , fLastEventCounter(0)
  , fSharedPadsOn(kTRUE)
  , fMaxADCgeom(0)
  , fADCnumber(0)
  , fCOL(0)
  , fExtendedCOL(0)
  , fIsShared(0)
  , fWarnError(kTRUE)
  , fWarnWarning(kFALSE)
  , fBufferRead(0)
  , fGeometry(0)
  , fRawReader(0)
  , fTRDfeeParam(0)
  , fCommonAdditive(0)
{
  //
  // default constructor
  //

  if (Init() == kFALSE) {
    AliWarning("Unable to Init.");	  
  }
}

//--------------------------------------------------------
AliTRDrawFastStream::AliTRDrawFastStream(AliRawReader *rawReader)
  : AliTRDrawStreamBase(rawReader)
  , fSM()
  , fStack(0)
  , fHC(0)
  , fLastHC(0)
  , fMCM()
  , fpPos(0)
  , fpBegin(0)
  , fpEnd(0)
  , fWordLength(0)
  , fpPosTemp(0)
  , fGlobalNTimeBins(0)
  , fIsTimeBinSet(kFALSE)
  , fStackNumber(-1)
  , fStackLinkNumber(-1)
  , fLinkTrackletCounter(-1)
  , fEndOfTrackletCount(-1)
  , fNWordsCounter(-1)
  , fMaskADCword(0)
  , fTbinADC(0)
  , fEventCounter(0)
  , fLastEventCounter(0)
  , fSharedPadsOn(kTRUE)
  , fMaxADCgeom(0)
  , fADCnumber(0)
  , fCOL(0)
  , fExtendedCOL(0)
  , fIsShared(0)
  , fWarnError(kTRUE)
  , fWarnWarning(kFALSE)
  , fBufferRead(0)
  , fGeometry(0)
  , fRawReader(rawReader)
  , fTRDfeeParam(0)
  , fCommonAdditive(0)
{
  //
  // default constructor
  //
  if (fRawReader) {
    if (Init() == kFALSE) {
      AliWarning("Unable to Init. Try setting up the reader with SetReader or buffer with Init(void *, UInt_t )");	  
    }
  }
  else {
    AliWarning("Unable to setup reader. Use SetReader(AliRawReader*).");
  }
}

//------------------------------------------------------------
AliTRDrawFastStream::AliTRDrawFastStream(const AliTRDrawFastStream& /*st*/)
  : AliTRDrawStreamBase()
  , fSM()
  , fStack(0)
  , fHC(0)
  , fLastHC(0)
  , fMCM()
  , fpPos(0)
  , fpBegin(0)
  , fpEnd(0)
  , fWordLength(0)
  , fpPosTemp(0)
  , fGlobalNTimeBins(0)
  , fIsTimeBinSet(kFALSE)
  , fStackNumber(-1)
  , fStackLinkNumber(-1)
  , fLinkTrackletCounter(-1)
  , fEndOfTrackletCount(-1)
  , fNWordsCounter(-1)
  , fMaskADCword(0)
  , fTbinADC(0)
  , fEventCounter(0)
  , fLastEventCounter(0)
  , fSharedPadsOn(kTRUE)
  , fMaxADCgeom(0)
  , fADCnumber(0)
  , fCOL(0)
  , fExtendedCOL(0)
  , fIsShared(0)
  , fWarnError(kTRUE)
  , fWarnWarning(kFALSE)
  , fBufferRead(0)
  , fGeometry(0)
  , fRawReader(0)
  , fTRDfeeParam(0)
  , fCommonAdditive(0)
{
  //
  // Copy constructor
  // 
  AliError("Not implemeneted.");
}

//------------------------------------------------------------
Bool_t AliTRDrawFastStream::SetRawVersion(Int_t fraw)
{
  //
  // function provided for backward compatibility
  //
  AliWarning("Raw data version is read from raw data stream! No point of setting it in here.");
  fraw = 0; // avoid warnings
  return kFALSE;
}

//------------------------------------------------------------
AliTRDrawFastStream::~AliTRDrawFastStream()
{
  //
  // destructor
  //
  delete fGeometry;
}

//------------------------------------------------------------

AliTRDrawFastStream &
AliTRDrawFastStream::operator=(const AliTRDrawFastStream &)
{
  //
  // we are not using this functionality
  //
  AliFatal("May not use.");
  return *this;
}

//___________________________________________________________
void AliTRDrawFastStream::SwapOnEndian()
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
Bool_t AliTRDrawFastStream::SkipWords(UInt_t iw)
{
  //
  // Skip words corresponding to iw
  //
  if ( fpPos + iw < fpEnd ) {
    fpPos += iw;
    return kTRUE;
  }
  else {
    if (fWarnWarning) AliWarning(Form("Skip %d words failed. %d available", iw, fpEnd - fpPos - 1));
    return kFALSE;
  }

  return kTRUE;
}

//------------------------------------------------------------
Bool_t AliTRDrawFastStream::SetReader(AliRawReader *reader)
{
  // 
  // set raw reader pointer
  // 
  if (reader != 0) {
    fRawReader = reader;
    if (fRawReader) {
      return Init();
    }
    else {
      AliWarning("Unable to setup reader.");
      return kFALSE;
    }
  }
  else {
    AliWarning("AliRawReader argument is 0.");
    fRawReader = 0;
  }

  return kFALSE;
}

//------------------------------------------------------------
Int_t AliTRDrawFastStream::NextBuffer()
{
  //
  // return -1 if no more buffers available
  // return  0 if SMHeader decoding failed 
  // return  1 if SMHeader dedoding is OK
  // 
  if (fRawReader != 0) {
    UChar_t *buffer = 0;
    UInt_t length = 0;
    Bool_t kBufferSet = fRawReader->ReadNextData(buffer);
    if (kBufferSet == kTRUE) {
      if (fgDebugFlag)  AliDebug(9, "Buffer is set.");
      length = fRawReader->GetDataSize();
      if (fgExtraSkip == kTRUE) {
        buffer += EXTRA_LEAD_WORDS * WORD_SIZE;
        length -= EXTRA_LEAD_WORDS * WORD_SIZE;
      }

      if (fgSkipCDH == kTRUE) {
        buffer += CDH_WORDS * WORD_SIZE;
        length -= CDH_WORDS * WORD_SIZE;	      
      }

      if (length > 0) {
        if (fgDebugFlag)  AliDebug(9, Form("Buffer length : %d", length));
        if (fgEnableMemoryReset) ResetMemory(); 
        if (DecodeSMHeader((void*)buffer, length) == kTRUE)
          return 1;
        else
          return 0;
      }
    }
    else {
      return -1;
    }
  }

  return -1;
}

//------------------------------------------------------------
void AliTRDrawFastStream::ResetCounters()
{
  //
  // reset some global counters
  //
  fBufferRead = kFALSE; // kFALSE if no buffer read

  fSM.fActiveStacks = 0;
  fSM.fNexpectedHalfChambers = 0;

  fLastEventCounter = 0;
  fEventCounter = 0;

  ResetIterators();
}

//------------------------------------------------------------
void AliTRDrawFastStream::ResetIterators()
{
  //
  // reset data which should be reset every sm
  //
  fStackNumber = 0;     
  fStackLinkNumber = 0; 
}

//------------------------------------------------------------
void AliTRDrawFastStream::ResetPerSM()
{
  //
  // reset every SM
  //
  fSM.fHeaderSize = 0;
  fSM.fTrackletEnable = kFALSE;
  fSM.fCorrupted = 0;
  fSM.fNexpectedHalfChambers = 0;
  fSM.fNexpectedHalfChambers = 0;
  fSM.fPos = NULL;
  for (Int_t i=0; i<5; i++) {
     fSM.fStackActive[i] = kFALSE;
  }
}

//------------------------------------------------------------
void AliTRDrawFastStream::ResetPerStack()
{
  //
  // reset every Stack
  //
  fStack->fHeaderSize = 0;
  fStack->fActiveLinks = 0;
  fStack->fPos = NULL;
  for (Int_t i=0; i<12; i++) {
     fStack->fLinksActive[i] = kFALSE;
     fStack->fLinksDataType[i] = 0;
     fStack->fLinksMonitor[i] = 0;
     fStack->fLinkMonitorError[i] = 0;
  }
}

//------------------------------------------------------------
void AliTRDrawFastStream::ResetPerHC()
{
  //
  // reset every HC
  //
  fEventCounter = 0;
  fHC->fNTracklets = 0;
  fHC->fSpecialRawV = 0;
  fHC->fRawVMajor = 0;
  fHC->fRawVMajorOpt = 0;
  fHC->fRawVMinor = 0;
  fHC->fNExtraWords = 0;
  fHC->fDCSboard = 0;
  fHC->fTimeBins = 0;
  fHC->fBunchCrossCounter = 0;
  fHC->fPreTriggerCounter = 0;
  fHC->fPreTriggerPhase = 0;
  fHC->fMCMmax = 0;

  fHC->fSM = 0;
  fHC->fStack = 0;
  fHC->fStackHCheader = 0;
  fHC->fLayer = 0;
  fHC->fLayerHCheader = 0;
  fHC->fSide = 0;
  fHC->fSideHCheader = 0;
  fHC->fDET = 0;
  fHC->fROC = 0;
  fHC->fRowMax = 0;
  fHC->fColMax = 0;

  fHC->fH0Corrupted = 0;
  fHC->fH1Corrupted = 0;
  fHC->fCorrupted = 0;
  fHC->fEOTECorrupted = kFALSE;
  fHC->fBufferCorrupted = kFALSE;
  fHC->fDataCorrupted = kFALSE;

  fHC->fNErrors= 0;
  memset(fHC->fErrorCodes, 0, 1411*sizeof(UShort_t)); // initialize error container
  memset(fHC->fTrackletWords, 0, MAXTRACKLETSPERHC*sizeof(UInt_t)); // initialize tracklet container
}

//------------------------------------------------------------
void AliTRDrawFastStream::ResetPerMCM()
{
  //
  // reset every MCM
  //
  fMCM.fROB = 0;
  fMCM.fMCM = 0;
  fMCM.fROW = 0;
  fMCM.fEvCounter = 0;
  fMCM.fADCMask = 0;
  fMCM.fADCMaskWord = 0;
  fMCM.fADCmax = 0;
  fMCM.fADCcount = 0;
  fMCM.fSingleADCwords = 0;
  fMCM.fMCMhdCorrupted = 0;
  fMCM.fADCmaskCorrupted = 0;
  fMCM.fDataCorrupted = kFALSE;
  fMCM.fPos = NULL;
  fMCM.fAdcDataPos = NULL;
  fMCM.fADCcounter = 0;

  memset(fMCM.fADCchannel, 0, TRDMAXADC*sizeof(UInt_t));
}

//------------------------------------------------------------
void AliTRDrawFastStream::ResetMemory()
{
  //
  // initialize all the data members to prevent read data
  // from previous buffer
  //
  ResetPerSM();
  for (Int_t istack=0; istack<5; istack++) {
     fStack = &fSM.fStacks[istack];
     ResetPerStack();
     for (Int_t ilink=0; ilink<12; ilink++) {
        fHC = &fStack->fHalfChambers[ilink];
        ResetPerHC();
     }
  }
}

//------------------------------------------------------------
Bool_t AliTRDrawFastStream::Next()
{
  //
  // returns with true on next adc read, returns false on errors and end of buffer
  // 
  if (fBufferRead) {
    while (fStackNumber < 5 && fSM.fActiveStacks > 0) {
      if (fSM.fStackActive[fStackNumber] == kTRUE) {
        fStack = &fSM.fStacks[fStackNumber];
        while (fStackLinkNumber < 12) {
          if (fStack->fLinksActive[fStackLinkNumber] == kTRUE) {
          //if (fStack->fLinksActive[fStackLinkNumber] == kTRUE && fStack->fLinksMonitor[fStackLinkNumber] == 0)
            fHC = &fStack->fHalfChambers[fStackLinkNumber];
		        //ResetPerHC(); // [mj - you don't need to do? ]
            if (!fHC) {
              AliError(Form("HC missing at stack %d link %d", fStackNumber, fStackLinkNumber));
              return kFALSE;
            }
            fStackLinkNumber++;
            return kTRUE;
          } //link active ?
          else fStackLinkNumber++;
        } //stack link number loop 
      } //stack active ?
      fStackNumber++;
      fStackLinkNumber = 0;
    } //stack number loop
  } //fBufferRead

  // go for the next buffer 
  if (fRawReader) {
    Int_t nextBuff = NextBuffer();
    while (nextBuff != -1) {
      if (nextBuff > 0) {
        fBufferRead = kTRUE;
        return Next();	   	  
      }
      nextBuff = NextBuffer();
    }
  }

  return kFALSE;
}

//------------------------------------------------------------
Int_t AliTRDrawFastStream::NextChamber(AliTRDdigitsManager *digitsManager, UInt_t **trackletContainer, UShort_t **errorCodeContainer) 
{
  
  //
  // Fills single chamber digit array 
  // Return value is the detector number
  //
  // first of all, you do the SM header decoding only at the beginning of the SM data reading
  // then, every HC, you call Next() which points next HC. then, there you decode the given HC 
  // and at the same time, put the digit into digitmanager 
  //
  AliTRDarrayADC *digits = 0;
  AliTRDarrayDictionary *track0 = 0;
  AliTRDarrayDictionary *track1 = 0;
  AliTRDarrayDictionary *track2 = 0; 
  AliTRDSignalIndex *indexes = 0;
  AliTRDdigitsParam *digitsparam = 0;

  Int_t lastdet = -1;
  Int_t det     = -1;
  Int_t side    = -1;
  //Int_t it = 0;
  Int_t ntracklets = 0;
  Int_t nErrors = 0;

  if (trackletContainer){
    for (Int_t i = 0; i < 2; i++)
       memset(trackletContainer[i], 0, MAXTRACKLETSPERHC*sizeof(UInt_t));
  }

  if (errorCodeContainer){
    for (Int_t i = 0; i < 2; i++)
       memset(errorCodeContainer[i], 0, 1411*sizeof(UShort_t));
  }

  while ( Next() ) { // loop over HC 

    // get this information from the GTU header
    det    = GetDet();
    side   = GetSide();

    if (det != lastdet) {
      // If new detector found
      if (lastdet == -1) {lastdet = det; fLastHC = fHC;}
      else {fStackLinkNumber--; fHC = fLastHC ; return lastdet;}

      if (det < 0 || det >= AliTRDgeometry::kNdet) continue;

      // Add a container for the digits of this detector
      digits = (AliTRDarrayADC *) digitsManager->GetDigits(det);

      if (digitsManager->UsesDictionaries()) {
        track0 = (AliTRDarrayDictionary *) digitsManager->GetDictionary(det,0);
        track1 = (AliTRDarrayDictionary *) digitsManager->GetDictionary(det,1);
        track2 = (AliTRDarrayDictionary *) digitsManager->GetDictionary(det,2);
      }

      if (!digits) return -1;

      //Int_t rowMax = GetRowMax();
      Int_t rowMax = fGeometry->RowmaxC1(); // we use maximum row number among all detectors to reuse memory
      Int_t colMax = GetColMax();
      Int_t ntbins = GetGlobalNTimeBins(); 

      // Set digitsparam variables
      digitsparam = (AliTRDdigitsParam *) digitsManager->GetDigitsParam();
      digitsparam->SetNTimeBins(det,ntbins);
      fCommonAdditive=10;
      digitsparam->SetADCbaseline(det,fCommonAdditive);

      // Allocate memory space for the digits buffer
      //if (digits->GetNtime() == 0) {
      if (ntbins != digits->GetNtime()) {
        digits->Allocate(rowMax, colMax, ntbins);
        if (digitsManager->UsesDictionaries()) {
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
      if (indexes->GetNtime() != ntbins)
        indexes->Allocate(rowMax, colMax, ntbins);
    }

    if (fSM.fTrackletEnable == kTRUE) { 
      if (DecodeTracklets() == kFALSE) {
        SeekEndOfData();

        if (fWarnWarning) AliWarning(Form("Tracklet decoding failed stack %d link %d", GetStack(), fStackLinkNumber));

        // copy error codes in memory into error container
        if (errorCodeContainer) {
          nErrors = GetNErrors();
          if(nErrors > 0) memcpy(errorCodeContainer[side], GetErrorCodes(), sizeof(UShort_t) * 1411); // [mj temp - optimize] 
        }
        
        continue; // if it fails to decode tracklets of this HC, it skips further decoding and goes to next HC  
      }     
    }   

    // decode hc data
    fgLastROB   = -1; // to check mcm number odering 
    fgLastIndex = -1 ; // to check mcm number odering 
    if (DecodeHC(digitsManager, digits, track0, track1, track2, indexes) == kFALSE) {
      // encode HC error code
      fHC->fErrorCodes[2] += fHC->fH0Corrupted;
      fHC->fErrorCodes[2] += (fHC->fH1Corrupted << 2);
      fHC->fErrorCodes[2] += (fHC->fCorrupted << 3);
      fHC->fErrorCodes[2] += ((fHC->fBufferCorrupted & 1) << 6);  
      fHC->fErrorCodes[2] += ((fHC->fEOTECorrupted & 1) << 7);  
      fHC->fErrorCodes[2] += ((fHC->fDataCorrupted & 1) << 8);  

	    if (fHC->fEOTECorrupted != kTRUE)  SeekEndOfData();

/*
      if (fWarnError) {
        AliError(Form("Failed HC : %s", DumpHCinfoH0(fHC)));
        AliError(Form("Failed HC : %s", DumpHCinfoH1(fHC)));
      }
*/ // [mj temp]
	  }
    else SeekEndOfData(); // make sure that finish off with the end of data markers

    // set pretrigger phase since it is only avaliable after decoding HC header 
    digitsparam->SetPretriggerPhase(det,GetPreTriggerPhase());

    // copy tracklets in memory into tracklet container
    if (trackletContainer) {
      ntracklets = GetNTracklets();
      // copy tracklet words to trackletContainer array 
      if(ntracklets > 0) memcpy(trackletContainer[side], GetTrackletWords(), sizeof(UInt_t) * ntracklets); 
    }

    // copy error codes in memory into error container
    if (errorCodeContainer) {
      nErrors = GetNErrors();
      if(nErrors > 0) memcpy(errorCodeContainer[side], GetErrorCodes(), sizeof(UShort_t) * 1411); 
    }

  }// end of while 

  return det;
}

//------------------------------------------------------------
Bool_t AliTRDrawFastStream::Init()
{
  //
  // Initialize geometry and fee parameters 
  //
  TDirectory *saveDir = gDirectory; 
  
  if (!fGeometry) {
    fGeometry = new AliTRDgeometry();
  }
  
  if (!fGeometry) {
    AliError("Geometry FAILED!");
    return kFALSE;
  }

  fTRDfeeParam = AliTRDfeeParam::Instance();
  if (!fTRDfeeParam) {
    AliError("AliTRDfeeParam FAILED!");
    return kFALSE;
  }

  fMaxADCgeom = (Int_t)fGeometry->ADCmax();

  ResetCounters(); // fBufferRead is set to kFALSE - important

  saveDir->cd();

  return kTRUE;
}

//------------------------------------------------------------
Bool_t AliTRDrawFastStream::InitBuffer(void *buffer, UInt_t length)
{
  // 
  // set initial information about the buffer
  //
  if (fgDebugFlag)  AliDebug(5, Form("Equipment ID: %d",fRawReader->GetEquipmentId()));
  if (fRawReader->GetEquipmentId()<1024 || fRawReader->GetEquipmentId()>1041) 
    return kFALSE; 

  ResetCounters();

  fpBegin = (UInt_t *)buffer;

  if (WORD_SIZE == 0) {
    AliFatal("Strange word size. size of UInt_t == 0");
    return kFALSE;
  }

  fWordLength = length/WORD_SIZE;
  fpEnd = fpBegin + fWordLength;
  fpPos = fpBegin;

  if (fpBegin == 0 || length <= 0) {
    AliError(Form("Buffer size or pointer is strange. pointer to the buffer is 0x%08x of size %d", fpBegin, length));
    return kFALSE;
  }

  SwapOnEndian();

  if (fgDumpHead >= 0) {
    if ( fgDumpHead == 0 ) { // dump all words
      AliInfo(Form("---------- Dumping all words from the beginnig of the buffer ----------"));
      if (DumpWords(fpBegin, fWordLength) == kFALSE) AliError("Dump failed. Not enough data.");
    } 
    else { 
      AliInfo(Form("---------- Dumping %u words from the beginnig of the buffer ----------",fgDumpHead));
      if (DumpWords(fpBegin, fgDumpHead) == kFALSE) AliError("Dump failed. Not enough data.");
    }
    AliInfo(Form("---------- Dumping ended ----------------------------------------------"));
  }

  return kTRUE;
}

//------------------------------------------------------------
Bool_t AliTRDrawFastStream::DecodeSMHeader(void *buffer, UInt_t length)
{
  // 
  // decode one sm data in buffer
  // 
  ResetIterators(); 

  if (InitBuffer(buffer, length) == kFALSE) {
    if (fWarnError) AliError("InitBuffer failed.");      
    return kFALSE;
  }

  if (DecodeGTUheader()== kFALSE)
    return kFALSE;

  Int_t nLinkErrors=0;
  for (Int_t istack = 0; istack < 5; istack++) {
     fStackNumber = istack; 
     if (fSM.fStackActive[istack] == kFALSE) continue;
      
     fStack = &fSM.fStacks[istack];
     
     fgLastHC  = -1; // to check rob number odering 
     for (Int_t ilink = 0; ilink < 12; ilink++) {
        fStackLinkNumber = ilink; 
        if (fStack->fLinksActive[ilink] == kFALSE) continue;
    
        // check GTU link monitor 
        if (!(fStack->fLinksDataType[ilink] == 0 && fStack->fLinksMonitor[ilink] == 0)) {
          fStack->fLinkMonitorError[ilink] = 1;
          fStack->fLinkMonitorError[ilink] += fNWordsCounter; // counts words of given hc having link monitor error
          nLinkErrors++;
          //continue;
        }

  	    if (fpPos >= fpEnd) {
          if (fRawReader) fRawReader->AddMajorErrorLog(kLinkDataMissing, "Link data missing");          
          if (fWarnWarning) AliWarning("Link data missing.");
          break;
        }

	      // set det number using SM header 
        fHC = &fStack->fHalfChambers[ilink];
	      fHC->fSM = fRawReader->GetEquipmentId() - 1024;
	      fHC->fStack = fStackNumber;
	      fHC->fLayer = Int_t(fStackLinkNumber/2.);
	      fHC->fSide = fStackLinkNumber%2;
	      fHC->fDET = fGeometry->GetDetector(fHC->fLayer, fHC->fStack, fHC->fSM);
	      fHC->fRowMax = fGeometry->GetRowMax(fHC->fLayer, fHC->fStack, fHC->fSM);
	      fHC->fROC    = fGeometry->GetDetectorSec(fHC->fLayer, fHC->fStack);
	      fHC->fColMax = fGeometry->GetColMax(fHC->fROC);
	   }
  }	

  // set number of timebin to be used in the digit container 
  if (nLinkErrors>30) fGlobalNTimeBins=30;
  else if (!fIsTimeBinSet) {
    fpPosTemp = fpPos;
    SetGlobalNTimebins();
    fIsTimeBinSet = kTRUE;
  }

  ResetIterators(); // need to do it again for Next() function

  return kTRUE;
}

//------------------------------------------------------------
Bool_t AliTRDrawFastStream::DecodeGTUheader()
{
  //
  // decode Supermodule Index Word
  //
  DecodeSMInfo(fpPos, &fSM);

  if (fgDebugFlag)  AliDebug(5, DumpSMInfo(&fSM));

  fpPos++;
  if (fpPos < fpEnd) {
    // fSM.fHeaderSize represent additional Supermodule header size which contains additional information regarding hardware design.
    // For the moment, we skip decoding these words 
    if (SkipWords(fSM.fHeaderSize) == kTRUE) {
      for (Int_t istack = 0; istack < 5; istack++) {
         if (fSM.fStackActive[istack] == kFALSE)
           continue;

         fStack = &fSM.fStacks[istack];

         // Decode Stack Index Word of given stack
         DecodeStackInfo(fpPos, fStack);
         fpPos++;

         fSM.fNexpectedHalfChambers += fStack->fActiveLinks;
        
         if (fgDebugFlag)  AliDebug(5, DumpStackInfo(fStack));
        
         if (SkipWords(fStack->fHeaderSize-6) == kFALSE) { // 6 is the 6 stack header words for 12 links 
           if (fRawReader) fRawReader->AddMajorErrorLog(kDecodeStackInfo, "Stack header words missing");
           return kFALSE;
         }
         for (Int_t iword=0; iword<6; iword++) { // decode 6 stack header words
            // Decode Stack Header Word of given stack
            DecodeStackHeader(fpPos, fStack, iword); 
            fpPos++;
         } // iword
      } // istack
    }
    else {
      return kFALSE;
    }
  }
  else {
    if (fWarnWarning) AliWarning("No additional sm headers and stack index words present.");
    if (fRawReader) fRawReader->AddMajorErrorLog(kDecodeStackInfo, "Stack info missing");
    return kFALSE;
  }

  if (fpPos < fpEnd) {
    if (fgDebugFlag)  AliDebug(5, "GTU headers are OK.");
  }
  else {
    if (fWarnError) AliError("No data just after GTU headers.");
    if (fRawReader) fRawReader->AddMajorErrorLog(kMissingData, "Missing sm data");
    return kFALSE;
  }

  if (fgDebugFlag)  AliDebug(5, Form("Expected half chambers from GTU header: %d", fSM.fNexpectedHalfChambers));

  return kTRUE;
}

//--------------------------------------------------------
void AliTRDrawFastStream::DecodeStackInfo(const UInt_t *word, struct AliTRDrawStack *st) const
{
  //
  // decode Stack #i Index Word
  // The Stack #i Index Word is a 32-Bit word with following structure
  // ssssssss ssssssss vvvv mmmm mmmmmmmm
  // s: Size of the Stack #i Header, v: Supermodule Header Version, m: Link Mask
  //
  st->fPos = (UInt_t*)word;

  UInt_t vword = *word;
  st->fHeaderSize = STACK_HEADER_SIZE(vword);

  UInt_t linkMask = STACK_LINK_WORD(vword);
  st->fActiveLinks = 0;
  for (Int_t i = 0; i < 12; i++) {
     if (IS_BIT_SET(linkMask,i) > 0) {
       st->fLinksActive[i] = kTRUE;
       st->fActiveLinks++;
     }
     else {
       st->fLinksActive[i] = kFALSE;
     }
  }
}

//--------------------------------------------------------
void AliTRDrawFastStream::DecodeStackHeader(const UInt_t *word, struct AliTRDrawStack *st, Int_t iword) const
{
  //
  // decode stack header
  //
  st->fPos = (UInt_t*)word;

  UInt_t vword = *word;
  st->fLinksDataType[2*iword]    = LINK0_DATA_TYPE_FLAG(vword);
  st->fLinksMonitor[2*iword]     = LINK0_MONITOR_FLAG(vword);
  st->fLinksDataType[2*iword+1]  = LINK1_DATA_TYPE_FLAG(vword);
  st->fLinksMonitor[2*iword+1]   = LINK1_MONITOR_FLAG(vword);
}

//------------------------------------------------------------
Bool_t AliTRDrawFastStream::DecodeTracklets()
{
  //
  // decode tracklets
  //
  fLinkTrackletCounter = 0; // tracklet counter of this link 
  fEndOfTrackletCount = 0;  // tracklet endmarker counter of this link
  fHC->fNTracklets = 0;     // number of tracklet of this link, should be less than 256

  if (fgDebugFlag)  AliDebug(10, Form("Decode tracklets at 0x%08x : 0x%08x", fpPos, *fpPos));

  while (!(*fpPos == END_OF_TRACKLET_MARKEROLD || *fpPos == END_OF_TRACKLET_MARKERNEW) && fpPos < fpEnd) {
    if (fgDebugFlag)  AliDebug(10, Form("Tracklet found at 0x%08x : 0x%08x", fpPos, *fpPos));

    fLinkTrackletCounter++;

    if (fLinkTrackletCounter > MAXTRACKLETSPERHC) {
      if (fgDebugFlag) AliDebug(11,Form("Max number of tracklets exceeded %d > %d. Tracklets are wrong either GTU header has problem",
                                        fLinkTrackletCounter, MAXTRACKLETSPERHC));
      if (fRawReader) fRawReader->AddMajorErrorLog(kTrackletOverflow,"Too many tracklets"); 
      fHC->fErrorCodes[1] = 1;          
      return kFALSE;
    }

    fHC->fTrackletWords[fLinkTrackletCounter-1] = UInt_t(*fpPos); //store tracklet words into memory 
    fpPos++;
  }

  fHC->fNTracklets = fLinkTrackletCounter;

  while ((*fpPos == END_OF_TRACKLET_MARKEROLD || *fpPos == END_OF_TRACKLET_MARKERNEW) && fpPos < fpEnd) {
    if (fgDebugFlag)  AliDebug(10, Form("EoTracklets found at 0x%08x : 0x%08x", fpPos, *fpPos));

    fEndOfTrackletCount++;
    fpPos++;
  }

  if (fEndOfTrackletCount < 2) {
    if (fgDebugFlag) AliDebug(11,"End of tracklets word missing"); 
    if (fRawReader) fRawReader->AddMajorErrorLog(kEOTrackeltsMissing, "End of tracklets word missing"); 
    fHC->fErrorCodes[1] += 2;          
    return kFALSE;
  }

  return kTRUE;
}

//--------------------------------------------------------
void AliTRDrawFastStream::DecodeSMInfo(const UInt_t *word, struct AliTRDrawSM *sm) const
{
  //
  // decode Supermodule Index Word
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

//------------------------------------------------------------
Bool_t AliTRDrawFastStream::DecodeHC(AliTRDdigitsManager *digitsManager, AliTRDarrayADC *digits, 
                             AliTRDarrayDictionary *track0, AliTRDarrayDictionary *track1, AliTRDarrayDictionary *track2, 
                             AliTRDSignalIndex *indexes)
{
  //
  // decode hc header and data
  //
  if (fpPos >= fpEnd) {
    fHC->fCorrupted += 1; 
    if (fWarnError) AliError("No data(including HC header) in the buffer");
    return kFALSE;;
  }

  if (DecodeHCheader() == kFALSE) {
    if (fWarnWarning) AliWarning(Form("HC Header decode failed. H0 Error: %d H1 Error: %d",fHC->fH0Corrupted,fHC->fH1Corrupted));
    return kFALSE;
  }
  else {
    fpPos++;
    if (fpPos >= fpEnd) {
      fHC->fCorrupted += 2; 
      if (fWarnError) AliError("No data right after HC header in the buffer");
      return kFALSE;
    }
  }

//   if ((fHC->fRawVMajor & 64) == 64) { // test pattern data
//     AliTRDrawTPStream *tpStream = new AliTRDrawTPStream(fHC->fRawVMajorOpt, fpPos);
//     if (tpStream->DecodeTPdata() == kFALSE) {
//       if (fWarnError) AliError("Failed to decode test pattern data");
//       return kFALSE;
//     }
//     return kTRUE;
//   }

  fHC->fMCMmax = 0; // count number of mcms in a hc 
  while (*fpPos != ENDOFRAWDATAMARKER && fpPos < fpEnd) {

    ResetPerMCM(); // reset for every mcm 

    if (fHC->fMCMmax > TRDMAXMCM) {
      fHC->fCorrupted += 4; 
      if (fgDebugFlag) AliDebug(11,"More mcm data than expected");
      return kFALSE;
    }

    if (DecodeMCMheader() == kFALSE) {

      // encode mcm level error codes
      fHC->fErrorCodes[fHC->fMCMmax+2] += fMCM.fMCMhdCorrupted;
      fHC->fErrorCodes[fHC->fMCMmax+2] += (fMCM.fADCmaskCorrupted << 4);
      fHC->fErrorCodes[fHC->fMCMmax+2] += ((fMCM.fDataCorrupted & 1) << 6);
      fHC->fErrorCodes[fHC->fMCMmax+2] += (fMCM.fMCM << 7);  // encode MCM number
      fHC->fErrorCodes[fHC->fMCMmax+2] += (fMCM.fROB << 11); // encode ROB number

      fHC->fMCMmax++; // increase mcm counter to match with expected rob/mcm number

      // in case we decide to keep reading data, skip this mcm data and find next mcm header
      if (fMCM.fADCmaskCorrupted < 2) {
        if (SkipMCMdata(fMCM.fADCcount*fMCM.fSingleADCwords) == kFALSE)
          return kFALSE;
        continue;
      }
      else {
        if (SeekNextMCMheader() == kFALSE)
          return kFALSE;
        continue;
      }
	  }

    fHC->fMCMmax++;

    if (fMCM.fADCmax > 0) {
      fpPos++;
      if (fpPos >= fpEnd) {
        fHC->fBufferCorrupted = kTRUE; 
        if (fgDebugFlag)  AliDebug(11, Form("Buffer ends in the middle of data"));
        return kFALSE;
      }

      fADCnumber = 0;
      for (Int_t iadc = 0; iadc < fMCM.fADCmax; iadc++) {
         fADCnumber = fMCM.fADCchannel[iadc];
           fExtendedCOL = fTRDfeeParam->GetExtendedPadColFromADC(fMCM.fROB, fMCM.fMCM, fADCnumber);
           fCOL = fTRDfeeParam->GetPadColFromADC(fMCM.fROB, fMCM.fMCM, fADCnumber);

         if (fADCnumber <= 1 || fADCnumber == fMaxADCgeom - 1)  // if adc number = 0, 1, 20
		       fIsShared = kTRUE;
	       else 
           fIsShared = kFALSE;

         if (fpPos + fMCM.fSingleADCwords >= fpEnd) {
           fHC->fBufferCorrupted = kTRUE;
           if (fgDebugFlag) AliDebug(11,"ADC (10 words) expected. Not enough data in the buffer.");
           return kFALSE;
         }

         //if (GetGlobalNTimeBins() < 31){
         if (fHC->fTimeBins < 31){
           if (DecodeADC(digitsManager, digits, track0, track1, track2, indexes) == kFALSE) {
             return kFALSE;
           }
         }
         //else if (GetGlobalNTimeBins() > 32) {
         else if (fHC->fTimeBins > 32) {
           if (DecodeADCExtended(digitsManager, digits, track0, track1, track2, indexes) == kFALSE) {
             return kFALSE;
           }
         }
         else { // nsamples = 31, 32 are not implemented in the TRAP and should never happen  
           if (fWarnError) AliError("nsamples are 31 or 32. These are not implemented in the TRAP and should never happen!");
         }

      } // iadc
    }
    else { // if there is no adc activated
      fpPos++;
    }
  } // mcm loop

  if (fpPos >= fpEnd) {
    fHC->fBufferCorrupted = kTRUE; 
    if (fgDebugFlag) AliDebug(11,"We are at the end of buffer. Not enough data in the buffer.");
    return kFALSE;
  }

  return kTRUE;
}

//------------------------------------------------------------
Bool_t AliTRDrawFastStream::DecodeHCheader()
{
  //
  // decode the half chamber header
  // if it fails to decode HC header for both H0 and H1, return kFALSE
  //
  if (DecodeHCwordH0(fpPos, fHC) == kFALSE)
    return kFALSE;

  if (fHC->fNExtraWords > 0) {
    fpPos++;
    if (fpPos < fpEnd) {
      if (DecodeHCwordH1(fpPos, fHC) == kFALSE)
        return kFALSE;
    }
    else {
      fHC->fBufferCorrupted = kTRUE;
      if (fgDebugFlag) AliDebug(11,"Expected HC header word H1. Fail due to buffer END.");
      return kFALSE;
    }
  }

  if (fgDebugFlag)  AliDebug(5, DumpHCinfoH0(fHC));
  if (fgDebugFlag)  AliDebug(5, DumpHCinfoH1(fHC));

  if (IsHCheaderOK() == kFALSE) {
    fHC->fH0Corrupted += 2;
    if (fgDebugFlag) AliDebug(11,Form("H0 Header Insane. Word 0x%08x", *fHC->fPos));
    return kFALSE;
  }

  return kTRUE;
}

//--------------------------------------------------------
Bool_t AliTRDrawFastStream::DecodeHCwordH0(const UInt_t *word, struct AliTRDrawHC *hc) const
{
  //
  // decode the hc header word 0
  //
  UInt_t vword = *word;
  hc->fPos[0] = (UInt_t*)word;

  hc->fH0Corrupted = HC_HEADER_MASK_ERR(vword);
  if (hc->fH0Corrupted > 0) {
    if (fgDebugFlag) AliDebug(11,Form("H0 Header Mask Error. Word 0x%08x",*fHC->fPos ));
    return kFALSE;
  }
  hc->fSpecialRawV =  HC_SPECIAL_RAW_VERSION(vword);
  hc->fRawVMajor = HC_MAJOR_RAW_VERSION(vword);
  hc->fRawVMajorOpt = HC_MAJOR_RAW_VERSION_OPT(vword);
  hc->fRawVMinor = HC_MINOR_RAW_VERSION(vword);
  hc->fNExtraWords = HC_EXTRA_WORDS(vword);
  hc->fDCSboard = HC_DCS_BOARD(vword);
  hc->fSMHCheader = HC_SM_NUMBER(vword);
  hc->fStackHCheader = HC_STACK_NUMBER(vword);
  hc->fLayerHCheader = HC_LAYER_NUMBER(vword);
  hc->fSideHCheader = HC_SIDE_NUMBER(vword);

  return kTRUE;
}

//--------------------------------------------------------
Bool_t AliTRDrawFastStream::DecodeHCwordH1(const UInt_t *word, struct AliTRDrawHC *hc) const
{
  //
  // decode the hc header word 1
  //
  UInt_t vword = *word;
  hc->fPos[1] = (UInt_t*)word;

  hc->fH1Corrupted = HC_HEADER_MASK_ERR(vword);
  if (hc->fH1Corrupted > 0) {
    if (fgDebugFlag) AliDebug(11,Form("H1 Header Mask Error. Word 0x%08x", *(fHC->fPos+1) ));
    return kFALSE;
  }
  hc->fTimeBins = HC_NTIMEBINS(vword);
  hc->fBunchCrossCounter = HC_BUNCH_CROSS_COUNTER(vword);
  hc->fPreTriggerCounter = HC_PRETRIGGER_COUNTER(vword);
  hc->fPreTriggerPhase = HC_PRETRIGGER_PHASE(vword);

  return kTRUE;
}

//------------------------------------------------------------
Bool_t AliTRDrawFastStream::IsHCheaderOK()
{
  //
  // check insanity of half chamber header
  //
  if (fHC->fStackHCheader < 0 || fHC->fStackHCheader > 4) {
    if (fgDebugFlag) AliDebug(11,Form("Wrong Stack %d", fHC->fStackHCheader));
      return kFALSE;
  }

  if (fHC->fLayerHCheader < 0 || fHC->fLayerHCheader >= AliTRDgeometry::kNlayer) {
    if (fgDebugFlag) AliDebug(11,Form("Wrong layer %d", fHC->fLayerHCheader));
      return kFALSE;
  }

  if (fHC->fSideHCheader < 0 || fHC->fSideHCheader > 1) {
    if (fgDebugFlag) AliDebug(11,Form("Wrong Side %d", fHC->fSideHCheader));
      return kFALSE;
  } 
    
  if (fHC->fSMHCheader != fHC->fSM) {
    if (fgDebugFlag) AliDebug(11,Form("Missmatch: SM number between HC header %d and GTU link mask %d",
                                      fHC->fSMHCheader, fHC->fSM));
    return kFALSE;
  }

  if (fgStackNumberChecker) {
    if (fHC->fStackHCheader != fHC->fStack) {
      if (fgDebugFlag) AliDebug(11,Form("Missmatch: Stack number between HC header %d and GTU link mask %d",
                                        fHC->fStackHCheader, fHC->fStack));
      return kFALSE;
    }
  }

  if (fgStackLinkNumberChecker) {
    if (fHC->fLayerHCheader * 2 + fHC->fSideHCheader != fHC->fLayer * 2 + fHC->fSide) {
      if (fgDebugFlag) AliDebug(11,Form("Missmatch: Layer number between HC header %d and GTU link mask %d | %s",
                                        fHC->fLayerHCheader, fHC->fLayer, DumpStackInfo(fStack)));
      return kFALSE;
    }
  }

  // SLOW GEOM : consistancy check with geometry
  if (fHC->fDET < 0 || fHC->fDET >= AliTRDgeometry::kNdet) {
    if (fgDebugFlag) AliDebug(11,Form("Wrong detector %d", fHC->fDET));
    if (fRawReader) fRawReader->AddMajorErrorLog(kHCHeaderWrongDet, "Wrong Det");
    return kFALSE;
  }

  if (fHC->fSM != fGeometry->GetSector(fHC->fDET) || fHC->fSM <0 || fHC->fSM >= AliTRDgeometry::kNsector) {
    if (fgDebugFlag) AliDebug(11,Form("Wrong SM(sector) %d (Geometry says: %d) Stack=%d Layer=%d Det=%d",
                                      fHC->fSM, fGeometry->GetSector(fHC->fDET), fHC->fStack, fHC->fLayer, fHC->fDET));
    if (fRawReader) fRawReader->AddMajorErrorLog(kHCHeaderWrongSM, "Wrong SM");
    return kFALSE;
  }

  if (fHC->fROC < 0) {
    if (fRawReader) fRawReader->AddMajorErrorLog(kHCHeaderWrongROC, "Wrong ROC");
    return kFALSE;
  }

  if (fHC->fRowMax < 1) {
    if (fRawReader) fRawReader->AddMajorErrorLog(kHCHeaderWrongROC, "Wrong ROC Row Max");
    return kFALSE;
  }

  if (fHC->fColMax < 1) {
    if (fRawReader) fRawReader->AddMajorErrorLog(kHCHeaderWrongROC, "Wrong ROC Col Max");
    return kFALSE;
  }

  return kTRUE;
}

//------------------------------------------------------------
Bool_t AliTRDrawFastStream::DecodeMCMheader()
{   
  //
  // decode the mcm header
  //
  DecodeMCMheader(fpPos, &fMCM);

  if (fHC->fEOTECorrupted == kTRUE) {
    fpPos--;
    return kFALSE;
  }

  fMCM.fROW = fTRDfeeParam->GetPadRowFromMCM(fMCM.fROB, fMCM.fMCM);

  if ((fHC->fRawVMajor > 2 && fHC->fRawVMajor <5) || ((fHC->fRawVMajor & 32) == 32)) { //cover old and new version definition of ZS data
    fpPos++;
    if ( fpPos < fpEnd ) {
      DecodeMask(fpPos, &fMCM);
      if (fHC->fEOTECorrupted == kTRUE) {
        fpPos--;
        return kFALSE;
      }
      MCMADCwordsWithTbins(fHC->fTimeBins, &fMCM);
      fMCM.fAdcDataPos = fpPos + 1;
    }
    else {
      if (fgDebugFlag) AliDebug(11,"Expected ADC mask word. Fail due to buffer END.");
      if (fRawReader) fRawReader->AddMajorErrorLog(kMCMADCMaskMissing,"Missing");
      fHC->fBufferCorrupted = kTRUE;
      return kFALSE;
    }
  }
  else {
    UInt_t dummyMask = MCM_DUMMY_ADCMASK_VAL;
    DecodeMask(&dummyMask, &fMCM);
    MCMADCwordsWithTbins(fHC->fTimeBins, &fMCM);
    fMCM.fAdcDataPos = fpPos + 1;
  }
  if (IsMCMheaderOK() == kFALSE)
      return kFALSE;

  return kTRUE;
}

//--------------------------------------------------------
void AliTRDrawFastStream::DecodeMCMheader(const UInt_t *word, struct AliTRDrawMCM *mcm) const
{
  //
  // decode the mcm header
  //
  UInt_t vword = *word;

  if (vword == END_OF_TRACKLET_MARKERNEW) {
    if (fWarnError) AliError(Form("There should be MCM header. We meet END_OF_TRACKLET_MARKER 0x%08x",vword));
    fHC->fEOTECorrupted = kTRUE; //to finish data reading of this HC
  }

  mcm->fMCMhdCorrupted = MCM_HEADER_MASK_ERR(vword); //if MCM header mask has error
  if (fgDebugFlag && mcm->fMCMhdCorrupted != 0) { 
    fHC->fDataCorrupted = kTRUE;
    AliDebug(11,Form("Wrong MCM header mask 0x%08x.\n", *fpPos));
  }

  mcm->fROB = MCM_ROB_NUMBER(vword);
  mcm->fMCM = MCM_MCM_NUMBER(vword);
  mcm->fEvCounter = MCM_EVENT_COUNTER(vword);
  mcm->fPos = (UInt_t*)word;
}

//--------------------------------------------------------
UInt_t AliTRDrawFastStream::GetMCMadcMask(const UInt_t *word, struct AliTRDrawMCM *mcm) const
{
  //
  // get the adc mask
  //
  UInt_t vword = *word;

  mcm->fADCMask   = 0;
  mcm->fADCcount  = 0;
  mcm->fADCMaskWord = vword;

  if (vword == END_OF_TRACKLET_MARKERNEW) {
    if (fWarnError) AliError(Form("There should be MCMadcMask. We meet END_OF_TRACKLET_MARKER 0x%08x",vword));
    fHC->fEOTECorrupted = kTRUE; //to finish data reading of this HC
  }

  if ( MCM_ADCMASK_MASK_ERR(vword) == 0 ) {
    mcm->fADCMask  = MCM_ADCMASK_VAL(vword);
    mcm->fADCcount = MCM_ADCMASK_NADC(~vword);
  }
  else {
    mcm->fADCMask = 0xffffffff;
    mcm->fADCmaskCorrupted = 1; // mcm adc mask error
    fHC->fDataCorrupted = kTRUE;
    if (fgDebugFlag) AliDebug(11,Form("Wrong ADC Mask word 0x%08x.\n", *fpPos));
  }

  return mcm->fADCMask;
}

//--------------------------------------------------------
void AliTRDrawFastStream::DecodeMask(const UInt_t *word, struct AliTRDrawMCM *mcm) const
{
  //
  // decode the adc mask - adcs to be read out
  //
  mcm->fSingleADCwords = 0;
  mcm->fADCmax = 0;
  mcm->fADCMask = GetMCMadcMask(word, mcm);

  if (mcm->fADCMask > 0) {
    for (Int_t i = 0; i < TRDMAXADC; i++) {
       mcm->fADCchannel[mcm->fADCmax] = 0;
       if ( IS_BIT_SET(mcm->fADCMask,i) ) {
         mcm->fADCchannel[mcm->fADCmax] = i;
         mcm->fADCmax++;
       }
    }
  }
  if (mcm->fADCcount != mcm->fADCmax && fHC->fRawVMajor >= 32) { // backward compatibility
    mcm->fADCmaskCorrupted += 2;
    fHC->fDataCorrupted = kTRUE;
    if (fgDebugFlag) AliDebug(11,Form("ADC counts from ADCMask are different %d %d : ADCMask word 0x%08x\n", 
                                      mcm->fADCcount, mcm->fADCmax, *fMCM.fPos));
  }
}

//--------------------------------------------------------
void AliTRDrawFastStream::MCMADCwordsWithTbins(UInt_t fTbins, struct AliTRDrawMCM *mcm) const
{
  //
  //  count the expected mcm words for a given tbins
  //
  mcm->fSingleADCwords = 0;
  mcm->fSingleADCwords = (fTbins-1)/3+1;
  if (fTbins > 32) mcm->fSingleADCwords = 10; // if the timebin is more than 30, then fix the number of adc words to 10
}

//------------------------------------------------------------
Bool_t AliTRDrawFastStream::IsMCMheaderOK()
{
  //
  // check the mcm header
  //
  if (fgLastROB != fMCM.fROB) {
    fgLastIndex = 0;
    if (fgLastROB== -1) fgLastROB = fMCM.fROB;
  }
  else {
    Int_t matchingcounter = 0; 
    for (Int_t i=fgLastIndex+1; i<16; i++) {
       if ( fMCM.fMCM == fgMCMordering[i] ) {
         fgLastIndex = i;
         matchingcounter++;
         break;
       }
    }
    if (matchingcounter == 0) {
      fMCM.fMCMhdCorrupted += 2;
      AliDebug(11,Form("MCM number from last MCM is larger: MCM # from current MCM %d \n", fMCM.fMCM));
    }
  } 
  
  if ( fgLastHC == fHC->fLayer*2 + fHC->fSide ) {
    if ( fMCM.fROB < fgLastROB ) {
      if((fMCM.fMCMhdCorrupted & 2) == 0) fMCM.fMCMhdCorrupted += 2;
      AliDebug(11,Form("ROB number from last MCM is larger: ROB # from current MCM %d \n", fMCM.fROB));
    }
    else fgLastROB = fMCM.fROB;
  }

  fgLastHC = fHC->fLayer*2 + fHC->fSide;

  if (fEventCounter == 0) {
    fEventCounter = fMCM.fEvCounter;
  }

  if (fEventCounter != fMCM.fEvCounter) {
    fMCM.fMCMhdCorrupted += 4;
    if (fgDebugFlag) AliDebug(11,Form("Event number(%d) of current MCM is different from that(%d) of reference MCM %s.\n"
                                      , fMCM.fEvCounter, fEventCounter, DumpMCMinfo(&fMCM)));
  }

  if (fEventCounter < fLastEventCounter) {
    fMCM.fMCMhdCorrupted += 8;
    if (fgDebugFlag) AliDebug(11,Form("Event from the past? Current %d Last %d %s.\n", fEventCounter, fLastEventCounter, DumpMCMinfo(&fMCM)));
  }

  if ( fMCM.fADCmaskCorrupted > 0 )
    return kFALSE;

  if ( fMCM.fMCMhdCorrupted > 0 )
    return kFALSE;

  return kTRUE;
}

//------------------------------------------------------------
Bool_t AliTRDrawFastStream::DecodeADC(AliTRDdigitsManager *digitsManager, AliTRDarrayADC *digits, 
                              AliTRDarrayDictionary *track0, AliTRDarrayDictionary *track1, AliTRDarrayDictionary *track2, 
                              AliTRDSignalIndex *indexes)
{
  //
  // decode single ADC channel
  //
  if(fADCnumber%2==1) fMaskADCword = ADC_WORD_MASK(ADCDATA_VAL1);
  if(fADCnumber%2==0) fMaskADCword = ADC_WORD_MASK(ADCDATA_VAL2);

  fTbinADC = 0;
  Bool_t isWritten = kFALSE;

  for (Int_t iw = 0; iw < fMCM.fSingleADCwords; iw++) {
     if (HC_HEADER_MASK_ERR(*fpPos) == 0 || *fpPos == END_OF_TRACKLET_MARKERNEW) {
       if (fWarnError) AliError(Form("There should be ADC data. We meet HC header or END_OF_TRACKLET_MARKER 0x%08x",*fpPos));
       fHC->fEOTECorrupted = kTRUE; 
       fpPos--;
       return kFALSE;
     }
     if (fMaskADCword != ADC_WORD_MASK(*fpPos)) {
       if (fgDebugFlag) AliDebug(11,Form("Wrong ADC data mask! [Expected mask: 0x%08x  Current mask: 0x%08x] ADC channel number: %02d MCM= %s ",
                                         fMaskADCword, ADC_WORD_MASK(*fpPos), fADCnumber, DumpMCMinfo(&fMCM)));
      // encode adc level error codes
       Int_t index = 21*(fMCM.fMCM + 16*int(fMCM.fROB/2)) + fADCnumber;
       fHC->fErrorCodes[index+66] += 1; 
       if (!isWritten) { 
         fHC->fErrorCodes[index+66] += (fADCnumber << 4);; 
         fHC->fErrorCodes[index+66] += (fMCM.fMCM << 9);; 
         fHC->fErrorCodes[index+66] += (fMCM.fROB << 13);; 
         isWritten = kTRUE; 
       }
       fMCM.fDataCorrupted = kTRUE;
       fHC->fDataCorrupted = kTRUE;
       fpPos++;
       continue;
     }
     // decode and put into the digit container
     Int_t adcSignals[3];
     adcSignals[0] = ((*fpPos & 0x00000ffc) >>  2);
     adcSignals[1] = ((*fpPos & 0x003ff000) >> 12);
     adcSignals[2] = ((*fpPos & 0xffc00000) >> 22);

     if(GetCol() < 0 || (!fSharedPadsOn & fIsShared)) {fpPos++; continue;};	

     for (Int_t i = 0; i < 3; i++) {
        digits->SetDataByAdcCol(GetRow(), GetExtendedCol(), fTbinADC + i, adcSignals[i]);
        indexes->AddIndexRC(GetRow(), GetCol());
        if (digitsManager->UsesDictionaries()) {
	  track0->SetData(GetRow(), GetCol(), fTbinADC + i, 0);
	  track1->SetData(GetRow(), GetCol(), fTbinADC + i, 0);
	  track2->SetData(GetRow(), GetCol(), fTbinADC + i, 0);
        }
     } // i
     fTbinADC += 3;
     fpPos++;
  } // iw

  return kTRUE;
}

//------------------------------------------------------------
Bool_t AliTRDrawFastStream::DecodeADCExtended(AliTRDdigitsManager *digitsManager, AliTRDarrayADC *digits, 
                              AliTRDarrayDictionary *track0, AliTRDarrayDictionary *track1, AliTRDarrayDictionary *track2, 
                              AliTRDSignalIndex *indexes)
{
  //
  // decode single ADC channel
  //
  if(fADCnumber%2==1) fMaskADCword = ADC_WORD_MASK(ADCDATA_VAL1);
  if(fADCnumber%2==0) fMaskADCword = ADC_WORD_MASK(ADCDATA_VAL2);

  Bool_t isWritten = kFALSE; //for error code recording

  fTbinADC = ((*fpPos & 0x000000fc) >>  2);
  fMCM.fSingleADCwords  = ((*fpPos & 0x00000f00) >>  8);
  
  Int_t adcFirst2Signals[2];
  adcFirst2Signals[0] = ((*fpPos & 0x003ff000) >> 12);
  adcFirst2Signals[1] = ((*fpPos & 0xffc00000) >> 22);

  for (Int_t i = 0; i < 2; i++) {
     digits->SetDataByAdcCol(GetRow(), GetExtendedCol(), fTbinADC + i, adcFirst2Signals[i]);
     indexes->AddIndexRC(GetRow(), GetCol());
     if (digitsManager->UsesDictionaries()) {
       track0->SetData(GetRow(), GetCol(), fTbinADC + i, 0);
       track1->SetData(GetRow(), GetCol(), fTbinADC + i, 0);
       track2->SetData(GetRow(), GetCol(), fTbinADC + i, 0);
     }
  } // i

  fpPos++;
  for (Int_t iw = 0; iw < fMCM.fSingleADCwords-1; iw++) {
     if (HC_HEADER_MASK_ERR(*fpPos) == 0 || *fpPos == END_OF_TRACKLET_MARKERNEW) {
       if (fWarnError) AliError(Form("There should be ADC data. We meet HC header or END_OF_TRACKLET_MARKER 0x%08x",*fpPos));
       fHC->fEOTECorrupted = kTRUE; 
       fpPos--;
       return kFALSE;
     }
     if (fMaskADCword != ADC_WORD_MASK(*fpPos)) {
       if (fgDebugFlag) AliDebug(11,Form("Wrong ADC data mask! [Expected mask: 0x%08x  Current mask: 0x%08x] ADC channel number: %02d MCM= %s ",
                                         fMaskADCword, ADC_WORD_MASK(*fpPos), fADCnumber, DumpMCMinfo(&fMCM)));
      // encode adc level error codes
       Int_t index = 21*(fMCM.fMCM + 16*int(fMCM.fROB/2)) + fADCnumber;
       fHC->fErrorCodes[index+66] += 1; 
       if (!isWritten) { 
         fHC->fErrorCodes[index+66] += (fADCnumber << 4);; 
         fHC->fErrorCodes[index+66] += (fMCM.fMCM << 9);; 
         fHC->fErrorCodes[index+66] += (fMCM.fROB << 13);; 
         isWritten = kTRUE; 
       }
       fMCM.fDataCorrupted = kTRUE;
       fHC->fDataCorrupted = kTRUE;
       fpPos++;
       continue;
     }
     // decode and put into the digit container
     Int_t adcSignals[3];
     adcSignals[0] = ((*fpPos & 0x00000ffc) >>  2);
     adcSignals[1] = ((*fpPos & 0x003ff000) >> 12);
     adcSignals[2] = ((*fpPos & 0xffc00000) >> 22);

     if(GetCol() < 0 || (!fSharedPadsOn & fIsShared)) {fpPos++; continue;};	
     for (Int_t i = 0; i < 3; i++) {
        digits->SetDataByAdcCol(GetRow(), GetExtendedCol(), fTbinADC + 2 + i, adcSignals[i]);
        indexes->AddIndexRC(GetRow(), GetCol());
        if (digitsManager->UsesDictionaries()) {
 	  track0->SetData(GetRow(), GetCol(), fTbinADC + 2 + i, 0);
	  track1->SetData(GetRow(), GetCol(), fTbinADC + 2 + i, 0);
	  track2->SetData(GetRow(), GetCol(), fTbinADC + 2 + i, 0);
        }
     } // i
     fTbinADC += 3;
     fpPos++;
  } // iw

  return kTRUE;
}

//------------------------------------------------------------
Bool_t AliTRDrawFastStream::SeekEndOfData()
{
  //
  // go to end of data marker
  //
  Int_t fEndOfDataCount = 0;
  fNWordsCounter = 0;

  while ( *fpPos != ENDOFRAWDATAMARKER && fpPos < fpEnd ) {
    fpPos++;
    fNWordsCounter++;
  }
  while (*fpPos == ENDOFRAWDATAMARKER && fpPos < fpEnd ) {
    fEndOfDataCount++;
    fpPos++;
  }

  return kTRUE;
}

//------------------------------------------------------------
Bool_t AliTRDrawFastStream::SeekNextMCMheader()
{
  //
  // go to mcm marker
  //
  fpPos++;

  while ( *fpPos != ENDOFRAWDATAMARKER && fpPos < fpEnd ) {
    if (MCM_HEADER_MASK_ERR(*fpPos) == 0 && MCM_HEADER_MASK_ERR(*(fpPos+1)) == 0) {
      if (fgDebugFlag) AliDebug(11,Form("Found : Pos 0x%08x : Val 0x%08x", fpPos, *fpPos));
      return kTRUE;
    }
    if ( *fpPos == END_OF_TRACKLET_MARKERNEW) {
      fHC->fEOTECorrupted = kTRUE;
      return kFALSE;
    }
    fpPos++;
  }

  SeekEndOfData();
  return kFALSE;
}

//------------------------------------------------------------
Bool_t AliTRDrawFastStream::SkipMCMdata(UInt_t iw)
{
  //
  // skip mcm data words due to corruption
  //
  if (fgDebugFlag) AliDebug(11,Form("Skip %d words due to MCM header corruption.",iw));
  UInt_t iwcounter = 0;
  while ( *fpPos != ENDOFRAWDATAMARKER && iwcounter < iw) {
    if ( *fpPos == END_OF_TRACKLET_MARKERNEW) {
      if (fgDebugFlag) AliDebug(11,"Met END_OF_TRACKLET_MARKERNEW");
      fHC->fEOTECorrupted = kTRUE;
      return kFALSE;
    }
    fpPos++;
    iwcounter++;
  } // while

  if (iwcounter == iw) {
    fpPos++;
    return kTRUE;
  }

  if (fgDebugFlag) AliDebug(11,"Met ENDOFRAWDATAMARKER");
  return kFALSE;
}

//------------------------------------------------------------
Bool_t AliTRDrawFastStream::SetGlobalNTimebins()
{
  // get number of time bin info from HC headers then set 
  Int_t nHCs=0;
  while (SetNTimebins()==kFALSE){
    if (fgDebugFlag) AliDebug(11,Form("Failed to get number of time bin information from the %sth HC",nHCs));
    nHCs++;
  }

  return kTRUE;
}

//------------------------------------------------------------
Bool_t AliTRDrawFastStream::SetNTimebins()
{
  // goes to the HC header position
  while (!(*fpPosTemp == END_OF_TRACKLET_MARKERNEW) && fpPosTemp < fpEnd) {
    fpPosTemp++;
  }
  while (*fpPosTemp == END_OF_TRACKLET_MARKERNEW) {
    fpPosTemp++;
  }
  // skip H0 
  fpPosTemp++;

  UInt_t vword = 0;
  if (!(vword = *fpPosTemp)) {
    fGlobalNTimeBins = 30; // default number of timebins 
    return kFALSE;
  }

  // get the number of time bins 
  if (HC_HEADER_MASK_ERR(vword) == 0) {
    fGlobalNTimeBins = HC_NTIMEBINS(vword);
    if (fGlobalNTimeBins > 64 || fGlobalNTimeBins < 10) return kFALSE; // minimal protection
  }
  else
    return kFALSE;

  return kTRUE;
}

//------------------------------------------------------------
Bool_t AliTRDrawFastStream::DumpWords(UInt_t *px, UInt_t iw, UInt_t marker)
{
  //
  // dump given number of words for debugging
  //
  TString tsreturn = Form("\n[ Dump Sequence at 0x%08x ] : ", px);
  for (UInt_t i = 0; i < iw; i++) {
     if ( iw != 0 && px + iw > fpEnd) 
       return kFALSE;
     if (i % 8 == 0) tsreturn += "\n                              ";
     if (marker != 0 && marker == px[i]) tsreturn += Form(" *>0x%08x<* ", px[i]);
     else tsreturn += Form("0x%08x ", px[i]);
  }
  tsreturn += "\n";

  AliInfo(tsreturn.Data());

  return kTRUE;
}

//--------------------------------------------------------
const char *AliTRDrawFastStream::DumpSMInfo(const struct AliTRDrawSM *sm)
{
  //
  // format the string with the sm info
  //
  return Form("[ SM Info 0x%08x] : Hsize %d TrackletEnable %d Stacks %d %d %d %d %d",
              *sm->fPos, sm->fHeaderSize, sm->fTrackletEnable,
              sm->fStackActive[0], sm->fStackActive[1], sm->fStackActive[2],
              sm->fStackActive[3], sm->fStackActive[4]);      
}

//--------------------------------------------------------
const char *AliTRDrawFastStream::DumpStackInfo(const struct AliTRDrawStack *st)
{
  //
  // format the string with the stack info
  //
  return Form("[ Stack Info 0x%08x ] : Hsize %d Links Active %d %d %d %d %d %d %d %d %d %d %d %d",
              *st->fPos, st->fHeaderSize,
              st->fLinksActive[0], st->fLinksActive[1], st->fLinksActive[2], st->fLinksActive[3],
              st->fLinksActive[4], st->fLinksActive[5], st->fLinksActive[6], st->fLinksActive[7],
              st->fLinksActive[8], st->fLinksActive[9], st->fLinksActive[10], st->fLinksActive[11]);

}
//--------------------------------------------------------
const char *AliTRDrawFastStream::DumpHCinfoH0(const struct AliTRDrawHC *hc)
{
  //
  // dump the hc header word 0 in strings
  //
  if (!hc)
    return Form("Unable to dump. Null received as parameter!?!");
  else
     return Form("[ HC[0] at 0x%08x ] : 0x%08x Info is : RawV %d SM %d Stack %d Layer %d Side %d DCSboard %d",
                hc->fPos[0], *(hc->fPos[0]), hc->fRawVMajor, fRawReader->GetEquipmentId()-1024, 
                hc->fStack, hc->fLayer, hc->fSide, hc->fDCSboard);
}

//--------------------------------------------------------
const char *AliTRDrawFastStream::DumpHCinfoH1(const struct AliTRDrawHC *hc)
{
  //
  // dump the hc header word 1 in strings
  //
  if (!hc)
    return Form("Unable to dump. Null received as parameter!?!");
  else
    return Form("[ HC[1] at 0x%08x ] : 0x%08x Info is : TBins %d BCcount %d PreTrigCount %d PreTrigPhase %d",
                hc->fPos[1], *(hc->fPos[1]), hc->fTimeBins, hc->fBunchCrossCounter, hc->fPreTriggerCounter, hc->fPreTriggerPhase);
}

//--------------------------------------------------------
const char *AliTRDrawFastStream::DumpMCMinfo(const struct AliTRDrawMCM *mcm)
{
  //
  // dump mcm info in strings
  //
  if (!mcm)
    return Form("Unable to dump. Null received as parameter!?!");
  else
    return Form("[ MCM 0x%08x ] : ROB %d MCM %d EvCounter %d", *(mcm->fPos), mcm->fROB, mcm->fMCM, mcm->fEvCounter);
}

//--------------------------------------------------------
const char *AliTRDrawFastStream::DumpMCMadcMask(const struct AliTRDrawMCM *mcm)
{
  //
  // mcm adc mask in strings
  //
  if (!mcm)
    return Form("Unable to dump. Null received as parameter!?!");

  TString tsreturn = Form("[Word] : 0x%08x => [Mask] : 0x%08x : ", mcm->fADCMaskWord, mcm->fADCMask);
  for (Int_t i = 0; i < TRDMAXADC; i++) {
     tsreturn += Form("%d ", mcm->fADCchannel[i]);
  }
  tsreturn += "";
  return tsreturn.Data();
}


