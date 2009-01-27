#ifndef ALITRDRAWSTREAMTB_H
#define ALITRDRAWSTREAMTB_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: AliTRDrawStream.h 27696 2008-07-31 09:18:53Z cblume $ */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// This class provides access to TRD digits in raw data.                     //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "TObject.h"
#include "TString.h"
#include "AliTRDrawStreamBase.h"

class AliTRDgeometry;
class AliRawReader;
class AliTRDdigitsManager;
class TTreeSRedirector;
class AliTRDfeeParam;


class AliTRDrawStream : public AliTRDrawStreamBase
{ // class def begin

 public:
    
  //--------------------------------------------------------

  // THE STRUCTURES
  
  //--------------------------------------------------------

  struct AliTRDrawADC
  {//adc struct

    UInt_t             *fPos;                    //! position of ADC 1st word in the buffer
    Short_t             fADCnumber;              // number of the ADC 0 .. 20
    Short_t             fCOL;                    // column - row from MCM
    Int_t               fSignals[TRD_MAX_TBINS]; // signals for this adc
    Bool_t              fIsShared;               // is pad chared between MCMs
    Short_t             fCorrupted;              // is adc word corrupted

    AliTRDrawADC()
      : fPos(0)
      , fADCnumber(0)
      , fCOL(0)
      , fSignals()
      , fIsShared(kFALSE)
      , fCorrupted(0)
    {
      // default constructor
    };

    AliTRDrawADC(const AliTRDrawADC& p): 
        fPos(p.fPos)
      , fADCnumber(p.fADCnumber)
      , fCOL(p.fCOL)
      , fSignals()
      , fIsShared(p.fIsShared)
      , fCorrupted(p.fCorrupted)
    {
      // copy constructor
    };

    AliTRDrawADC &operator=(const AliTRDrawADC &) 
    {
      // assignment operator
      // not implemented
      return *this;
    };

  };
  
  //--------------------------------------------------------

  struct AliTRDrawMCM
  { // mcm struct
    Int_t               fROB;                     // ROB number
    Int_t               fMCM;                     // MCM number
    Int_t               fROW;                     // row number filed during decoding
      
    UInt_t              fEvCounter;               // MCM event counter
    UInt_t              fADCMask;                 // ADC mask
    UInt_t              fADCMaskWord;             // word with ADC mask in
    UInt_t              fADCchannel[TRD_MAX_ADC]; // channels to be decoded accrording to ADC mask
      
    Int_t               fADCmax;                  // number of ADCs fired
    Int_t               fADCcount;                // number of ADCs fired from double checking bit
    Int_t               fMCMADCWords;             // mcm words to expect
    Int_t               fSingleADCwords;          // n of words per ADC
      
    Int_t               fMCMhdCorrupted;          // is mcm header corrupted
    Int_t               fADCmaskCorrupted;        // is mcm adc mask corrupted
    Int_t               fCorrupted;               // is mcm data missing
      
    UInt_t             *fPos;                     //! position of mcm header in the buffer
    UInt_t             *fAdcDataPos;              //! start of ADC data for this mcm

    Int_t               fADCcounter;              // count the adcs decoded
    AliTRDrawADC        fADCs[TRD_MAX_ADC];       // 21 adcs
      
    AliTRDrawMCM()
      : fROB(-1)
      , fMCM(-1)
      , fROW(-1)
      , fEvCounter(0)
      , fADCMask(0)
      , fADCMaskWord(0)
      , fADCchannel()      
      , fADCmax(0)
      , fADCcount(0)
      , fMCMADCWords(0)      
      , fSingleADCwords(0)
      , fMCMhdCorrupted(0)      
      , fADCmaskCorrupted(0)      
      , fCorrupted(0)      
      , fPos(0)
      , fAdcDataPos(0)
      , fADCcounter(0)
      , fADCs()
    {
      // default constructor
    };

    AliTRDrawMCM(const AliTRDrawMCM & p):
        fROB(p.fROB)
      , fMCM(p.fMCM)
      , fROW(p.fROW)
      , fEvCounter(p.fEvCounter)
      , fADCMask(p.fADCMask)
      , fADCMaskWord(p.fADCMaskWord)
      , fADCchannel()      
      , fADCmax(p.fADCmax)
      , fADCcount(p.fADCcount)
      , fMCMADCWords(p.fMCMADCWords)      
      , fSingleADCwords(p.fSingleADCwords)
      , fMCMhdCorrupted(p.fMCMhdCorrupted)      
      , fADCmaskCorrupted(p.fADCmaskCorrupted)      
      , fCorrupted(p.fCorrupted)      
      , fPos(p.fPos)
      , fAdcDataPos(p.fAdcDataPos)
      , fADCcounter(p.fADCcounter)
      , fADCs()
    {
      // copy constructor
    };

    AliTRDrawMCM &operator=(const AliTRDrawMCM &)
    {
      // assignment operator
      // not implemented
      return *this;
    };

  };

  //--------------------------------------------------------

  struct AliTRDrawHC
  { // hc struct
      
    //tacklet words of given HC
    UInt_t              fTrackletWords[MAX_TRACKLETS_PERHC]; // array to keep tracklet words [mj]
    Short_t             fTrackletError;                      // tracklet error 
    Short_t             fNTracklets;                    // number of tracklet  

    // header word 0
    Int_t               fSpecialRawV;       // Raw data version
    Int_t               fRawVMajor;         // Raw data version
    Int_t               fRawVMajorOpt;      // Raw data version
    Int_t               fRawVMinor;         // Raw data version
    Int_t               fNExtraWords;       // N extra HC header words
    Int_t               fDCSboard;          // DCS board number
    Int_t               fSM;                // Super Module number
    Int_t               fStack;             // Stack number (some people might call it a chamber)
    Int_t               fLayer;             // Layer number (some people might call it a plane)
    Int_t               fSide;              // Side of HC
      
    // header word 1
    Int_t               fTimeBins;          // N of t bins
    UInt_t              fBunchCrossCounter; // Bunch crossing counter
    UInt_t              fPreTriggerCounter; // Pre Trigger counter
    UInt_t              fPreTriggerPhase;   // Pre Trigger phase
      
    // error 
    Int_t               fH0Corrupted;       // is hc header 0 corrupted 
    Int_t               fH1Corrupted;       // is hc header 1 corrupted
    Int_t               fCorrupted;         // is hc data corrupted 

    UInt_t             *fPos[2];            //! position of the header words in buffer
      
    Int_t               fDET;               // filled while decoding
    Int_t               fROC;               // filled while decoding
    Int_t               fRowMax;            // filled while decoding
    Int_t               fColMax;            // filled while decoding

    // hc data
    Int_t               fMCMmax;            // number of mcm found
    AliTRDrawMCM        fMCMs[TRD_MAX_MCM]; // 4 ROBS 16 each 

    AliTRDrawHC()
      : fTrackletWords() //[mj]
      , fTrackletError(0)
      , fNTracklets(0)
      , fSpecialRawV(0)
      , fRawVMajor(0)
      , fRawVMajorOpt(0)
      , fRawVMinor(0)
      , fNExtraWords(0)
      , fDCSboard(-1)
      , fSM(-1)
      , fStack(-1)
      , fLayer(-1)
      , fSide(-1)
      , fTimeBins(0)
      , fBunchCrossCounter(0)
      , fPreTriggerCounter(0)
      , fPreTriggerPhase(0)
      , fH0Corrupted(0)
      , fH1Corrupted(0)
      , fCorrupted(0)
      , fPos()
      , fDET(-1)
      , fROC(-1)
      , fRowMax(-1)
      , fColMax(-1)
      , fMCMmax(0)
      , fMCMs()
    {
      // default constructor hc info 
    };

    AliTRDrawHC(const AliTRDrawHC & p):
        fTrackletWords() //[mj]
      , fTrackletError(p.fTrackletError)
      , fNTracklets(p.fNTracklets)
      , fSpecialRawV(p.fSpecialRawV)
      , fRawVMajor(p.fRawVMajor)
      , fRawVMajorOpt(p.fRawVMajorOpt)
      , fRawVMinor(p.fRawVMinor)
      , fNExtraWords(p.fNExtraWords)
      , fDCSboard(p.fDCSboard)
      , fSM(p.fSM)
      , fStack(p.fStack)
      , fLayer(p.fLayer)
      , fSide(p.fSide)
      , fTimeBins(p.fTimeBins)
      , fBunchCrossCounter(p.fBunchCrossCounter)
      , fPreTriggerCounter(p.fPreTriggerCounter)
      , fPreTriggerPhase(p.fPreTriggerPhase)
      , fH0Corrupted(p.fH0Corrupted)
      , fH1Corrupted(p.fH1Corrupted)
      , fCorrupted(p.fCorrupted)
      , fPos()
      , fDET(p.fDET)
      , fROC(p.fROC)
      , fRowMax(p.fRowMax)
      , fColMax(p.fColMax)
      , fMCMmax(p.fMCMmax)
      , fMCMs()
    {
      // copy constructor
    };

    AliTRDrawHC &operator=(const AliTRDrawHC &)
    {
      // assignment operator
      // not implemented
      return *this;
    };

  };

  //--------------------------------------------------------
    
  struct AliTRDrawStack
  {
    UInt_t           fHeaderSize;            // header size of the stack info
    Bool_t           fLinksActive[12];       // data links active - 1 per half chamber
    Short_t          fLinksDataType[12];     // 0 indicating real data for the front-end electronics 
    Short_t          fLinksMonitor[12];      // 0 indicating properly operating link 
    Short_t          fLinkMonitorError[12];  // record link monitor error
    Int_t            fActiveLinks;           // number of active links
    UInt_t          *fPos;                   //! position in the buffer
	            
    AliTRDrawHC      fHalfChambers[12];      // 12 half chambers in a stack
      
    AliTRDrawStack()
      : fHeaderSize(0)
      , fLinksActive()
      , fLinksDataType()
      , fLinksMonitor()
      , fLinkMonitorError()
      , fActiveLinks(0)
      , fPos(0)
      , fHalfChambers()
    {
      // default constructor
    };      

    AliTRDrawStack(const AliTRDrawStack & p):
        fHeaderSize(p.fHeaderSize)
      , fLinksActive()
      , fLinksDataType()
      , fLinksMonitor()
      , fLinkMonitorError()
      , fActiveLinks(p.fActiveLinks)
      , fPos(p.fPos)
      , fHalfChambers()
    {
      // copy constructor
    };

    AliTRDrawStack &operator=(const AliTRDrawStack &)
    {
      // assignment operator
      // not implemented
      return *this;
    };

  };

  //--------------------------------------------------------

  struct AliTRDrawSM
  {
    UInt_t            fHeaderSize;            // size of the header in words
    Bool_t            fTrackletEnable;        // tracklet enable bit
    Bool_t            fStackActive[5];        // map of active/expected stacks
    Int_t             fActiveStacks;          // number of active stacks
    Int_t             fCorrupted;             // is sm info corrupted
    Int_t             fNexpectedHalfChambers; // number of half chambers to be read out in this sm
    Bool_t            fClean;                 // true if everything went OK - false is some error occured
    UInt_t           *fPos;                   // location of the sm info - should be the first word (after CDH if not DDL buffer)

    AliTRDrawStack    fStacks[5];             

    AliTRDrawSM()
      : fHeaderSize(0)
      , fTrackletEnable(0)
      , fStackActive()
      , fActiveStacks(0)
      , fCorrupted(0)
      , fNexpectedHalfChambers(0)
      , fClean(kTRUE)
      , fPos(0)
      , fStacks()
      {
	// Default constructor
      };      

    AliTRDrawSM(const AliTRDrawSM & p):
        fHeaderSize(p.fHeaderSize)
      , fTrackletEnable(p.fTrackletEnable)
      , fStackActive()
      , fActiveStacks(p.fActiveStacks)
      , fCorrupted(p.fCorrupted)
      , fNexpectedHalfChambers(p.fNexpectedHalfChambers)
      , fClean(p.fClean)
      , fPos(p.fPos)
      , fStacks()
    {
      // copy constructor
    };

    AliTRDrawSM &operator=(const AliTRDrawSM &)
    {
      // assignment operator
      // not implemented
      return *this;
    };

  };
  
  //--------------------------------------------------------
     
  AliTRDrawStream();
  AliTRDrawStream(AliRawReader *rawReader);
  virtual ~AliTRDrawStream();

  //--------------------------------------------------------

  virtual Bool_t       Next();                           // read the next data in the memory
  //virtual Int_t      NextChamber(AliTRDdigitsManager *man); // read next chamber data in the momory
  virtual Int_t        NextChamber(AliTRDdigitsManager *man, UInt_t **trackletContainer); // read next chamber data in the memory
  virtual Bool_t       Init();                           // initialize some internal variables

  Int_t    NextBuffer(); // go and init next buffer if available - check the implementation file for return values

  Bool_t   SetRawVersion(Int_t fraw); // set the raw version - used for backward compat.
  
  Bool_t   IsCurrentPadShared() const {return fADC->fIsShared;} // is current pad shared between mcms
  void     SetSharedPadReadout(Bool_t fv) {fSharedPadsOn = fv;} // set the flag on if the reader should return the shared pads
  
  Bool_t   DecodeSM(void *buffer, UInt_t length); // decode a buffer
  Int_t    DecodeSM();                            // used with raw reader
  Int_t    DecodeSM(AliRawReader *reader);        // used with raw reader
	   
  Bool_t   SetReader(AliRawReader *reader); // set the raw reader to use
    	   
  // info from Supermodule Index Word
  Bool_t    IsTrackletEnableBitSet() const {return fSM.fTrackletEnable;} // get status of tracklet enable bit
  Bool_t    IsStackActive(Int_t is) const {return fSM.fStackActive[is];} // get status of stack enable bit
  Int_t     GetNofActiveStacks() const {return fSM.fActiveStacks;}       // get number of active stacks from stack mask
  UInt_t   *GetGTUheaderWords() const {return fSM.fPos;}       // get number of active stacks from stack mask

  // info from Stack Index Word
  Int_t     GetNexpectedHalfChambers() const {return fSM.fNexpectedHalfChambers;}                    // get number of expected HC in a sm
  Int_t     GetNofActiveLinksInStack(Int_t is) const {return fSM.fStacks[is].fActiveLinks;}          // get number of active links in a stack
  Bool_t    IsLinkActiveInStack(Int_t is, Int_t il) const {return fSM.fStacks[is].fLinksActive[il];} // check whether the link is active

  // info from Stack Header Word
  Short_t    GetLinkMonitorError(Int_t is, Int_t il) const {return fSM.fStacks[is].fLinkMonitorError[il];} // get link monitor error

  // info from Tracklet Data
  Int_t     GetTrackletErrorCode(Int_t is, Int_t il) const {return fSM.fStacks[is].fHalfChambers[il].fTrackletError;}
  Int_t     GetNTracklets(Int_t is, Int_t il) const {return fSM.fStacks[is].fHalfChambers[il].fNTracklets;} // get number of tracklets 

  // info from HC Header Word
  Int_t     GetSM(Int_t is, Int_t il) const {return fSM.fStacks[is].fHalfChambers[il].fSM;}
  Int_t     GetLayer(Int_t is, Int_t il) const {return fSM.fStacks[is].fHalfChambers[il].fLayer;}
  Int_t     GetStack(Int_t is, Int_t il) const {return fSM.fStacks[is].fHalfChambers[il].fStack;}
  Int_t     GetSide(Int_t is, Int_t il) const {return fSM.fStacks[is].fHalfChambers[il].fSide;}
  Int_t     GetH0ErrorCode(Int_t is, Int_t il) const {return fSM.fStacks[is].fHalfChambers[il].fH0Corrupted;}
  Int_t     GetH1ErrorCode(Int_t is, Int_t il) const {return fSM.fStacks[is].fHalfChambers[il].fH1Corrupted;}
  Int_t     GetNumberOfTimeBins(Int_t is, Int_t il) const { return fSM.fStacks[is].fHalfChambers[il].fTimeBins;}
  UInt_t   *GetTrackletWords(Int_t is, Int_t il) { return fSM.fStacks[is].fHalfChambers[il].fTrackletWords;}

  // info from HC data
  Int_t     GetHCErrorCode(Int_t is, Int_t il) const {return fSM.fStacks[is].fHalfChambers[il].fCorrupted;}
  Int_t     GetHCMCMmax(Int_t is, Int_t il) const {return fSM.fStacks[is].fHalfChambers[il].fMCMmax;}

  // from MCM Header Word
  // rob and mcm ordering
  // side 0(even link) - ROB: 0 2 4 6  MCM: 12 13 14 15 8 9 10 11 4 5 6 7 0 1 2 3  
  // side 1( odd link) - ROB: 1 3 5 7  MCM: 12 13 14 15 8 9 10 11 4 5 6 7 0 1 2 3  
  Int_t     GetMCM(Int_t stack, Int_t link, Int_t mcm) {return fSM.fStacks[stack].fHalfChambers[link].fMCMs[mcm].fMCM;}
  Int_t     GetROB(Int_t stack, Int_t link, Int_t mcm) {return fSM.fStacks[stack].fHalfChambers[link].fMCMs[mcm].fROB;}
  Int_t     GetMCMhdErrorCode(Int_t stack, Int_t link, Int_t mcm) {return fSM.fStacks[stack].fHalfChambers[link].fMCMs[mcm].fMCMhdCorrupted;}
  Int_t     GetMCMADCMaskErrorCode(Int_t stack, Int_t link, Int_t mcm) {return fSM.fStacks[stack].fHalfChambers[link].fMCMs[mcm].fADCmaskCorrupted;}
  Int_t     GetEventNumber(Int_t stack, Int_t link, Int_t mcm) {return fSM.fStacks[stack].fHalfChambers[link].fMCMs[mcm].fEvCounter;}
  Int_t     GetADCcount(Int_t stack, Int_t link, Int_t mcm) {return fSM.fStacks[stack].fHalfChambers[link].fMCMs[mcm].fADCcount;}

  // info from MCM data words
  Int_t     GetMCMErrorCode(Int_t stack, Int_t link, Int_t mcm) {return fSM.fStacks[stack].fHalfChambers[link].fMCMs[mcm].fCorrupted;} // get MCM data error code
  Int_t     GetADCErrorCode(Int_t stack, Int_t link, Int_t mcm, Int_t adc) {return fSM.fStacks[stack].fHalfChambers[link].fMCMs[mcm].fADCs[adc].fCorrupted;} // get ADC error code
  Int_t     GetADCnumber(Int_t stack, Int_t link, Int_t mcm, Int_t adc) {return fSM.fStacks[stack].fHalfChambers[link].fMCMs[mcm].fADCs[adc].fADCnumber;} // get ADC error code

  Int_t     GetRow(Int_t stack, Int_t link, Int_t mcm) {return fSM.fStacks[stack].fHalfChambers[link].fMCMs[mcm].fROW;}         // get current row number
  Int_t     GetCol(Int_t stack, Int_t link, Int_t mcm, Int_t adc) {return fSM.fStacks[stack].fHalfChambers[link].fMCMs[mcm].fADCs[adc].fCOL;}         // get current column number

  // info from ADC data words
  Int_t    *GetSignalDirect(Int_t stack, Int_t link, Int_t mcm, Int_t adc) {return fSM.fStacks[stack].fHalfChambers[link].fMCMs[mcm].fADCs[adc].fSignals;}


  // from here, only works with returning ADC channel pointer using Next() 
  UInt_t   *GetTrackletWords() const { return fHC->fTrackletWords;}                 // return tracklet words pointer per hc [mj]
  Int_t     GetTrackletErrorCode() const {return fHC ? fHC->fTrackletError : -1;}   // get tracklet error code
  Int_t     GetNTracklets() const {return fHC ? fHC->fNTracklets : -1;}             // get number of tracklets 

  Int_t     GetSpecialRawVersion() const {return fHC ? fHC->fSpecialRawV : -1;}     // return special raw version
  Int_t     GetMajorRawVersion() const {return fHC ? fHC->fRawVMajor : -1;}         // major raw version getter
  Int_t     GetRawVersion() const {return fHC ? fHC->fRawVMajor : -1;}              // compatibility see funtion above
  Int_t     GetMinorRawVersion() const {return fHC ? fHC->fRawVMinor : -1;}         // minor raw version

  Int_t     GetSM() const {return fHC ? fHC->fSM : -1;}                              //  SM Position of CURRENT half chamber in full TRD
  Int_t     GetLayer() const {return fHC ? fHC->fLayer : -1;}                        //  Layer Position of CURRENT half chamber in full TRD
  Int_t     GetStack() const {return fHC ? fHC->fStack : -1;}                        //  Stack Position of CURRENT half chamber in full TRD
  Int_t     GetSide() const {return fHC ? fHC->fSide : -1;}                          // get side
  Int_t     GetDCS() const { return fHC ? fHC->fDCSboard : -1;}                      //  DCS board number read from data (HC header)
  Int_t     GetROC() const { return fHC ? fHC->fROC : -1;}                           //  ROB Position of CURRENT half chamber in full TRD
  Int_t     GetNumberOfTimeBins() const { return fHC ? fHC->fTimeBins : 0;}          // Get Ntime bins
  UInt_t    GetBunchCrossCounter() const {return fHC ? fHC->fBunchCrossCounter : 0;} // get bunch cross counter
  UInt_t    GetPreTriggerCounter() const {return fHC ? fHC->fPreTriggerCounter : 0;} // get pre trigger info
  UInt_t    GetPreTriggerPhase() const {return fHC ? fHC->fPreTriggerPhase : 0;}     // get trigger phase

  Int_t     GetRow() const {return fMCM ? fMCM->fROW : -1;}         // get current row number
  Int_t     GetCol() const {return fADC ? fADC->fCOL : -1;}         // get current column number
  Int_t     GetRowMax() const { return fHC ? fHC->fRowMax : -1;}    // Get maximum rows in the current HC
  Int_t     GetColMax() const { return fHC ? fHC->fColMax : -1;}    // Get maximum cols in the current HC
  // compatibility
  Int_t     GetMaxRow() const { return fHC ? fHC->fRowMax : -1;}    // Get maximum rows in the current HC
  Int_t     GetMaxCol() const { return fHC ? fHC->fColMax : -1;}    // Get maximum cols in the current HC

  UInt_t    GetHCword0() const {return fHC ? *fHC->fPos[0] : 0;}    // get the HC word 0
  UInt_t    GetHCword1() const {return fHC ? *fHC->fPos[1] : 0;}    // get the HC word 1
	    
  Int_t     GetDET() const {return fHC ? fHC->fDET : -1;}           // get current det number
  Int_t     GetDet() const {return fHC ? fHC->fDET : -1;}           // get current det number
    	    
  Int_t     GetROB() const {return fMCM ? fMCM->fROB : -1;}         // get current ROB number
  Int_t     GetMCM() const {return fMCM ? fMCM->fMCM : -1;}         // get current MCM number
  Int_t     GetEventNumber() const { return fMCM->fEvCounter;}      //  MCM Event number and position of current MCM on TRD chamber

  Int_t     GetADC() const { return fADC ? fADC->fADCnumber : -1;}  //  MCM ADC channel and Time Bin of word 1
  Int_t     GetTimeBin() const { return 0;}                         //  MCM ADC channel and Time Bin of word 1
  Int_t    *GetSignals() const { return fADC ? fADC->fSignals : (Int_t *)fgEmptySignals;} // signals in the n-time bins from data word

  Int_t     GetHCErrorCode() const {return fHC ? fHC->fCorrupted : -1;}    // get HC error code
  Int_t     GetH0ErrorCode() const {return fHC ? fHC->fH0Corrupted : -1;}  // get HC header word0 error code
  Int_t     GetH1ErrorCode() const {return fHC ? fHC->fH1Corrupted : -1;}  // get HC header word1 error code
  Int_t     GetMCMErrorCode() const {return fMCM ? fMCM->fCorrupted : -1;} // get MCM data error code
  Int_t     GetADCErrorCode() const {return fADC ? fADC->fCorrupted : -1;} // get ADC data error code
  Int_t     GetMCMhdErrorCode() const {return fMCM ? fMCM->fMCMhdCorrupted: -1;} // get MCM header word error code
  Int_t     GetMCMADCMaskErrorCode() const {return fMCM ? fMCM->fADCmaskCorrupted: -1;} // get MCM adc mask error code

  UInt_t   *GetSMstreamPosition() const {return fSM.fPos;} // get position of the SM index word in the buffer

  Bool_t    IsSMbufferClean() const {return fSM.fClean;}   // is data clean

  //----------------------------------------------------------
 
  static void    SetNoDebug() {fgDebugFlag = kFALSE;} // allow debug info
  static void    EnableMemoryReset() {fgEnableMemoryReset = kTRUE;} // allow memory reset
  static void    SetNoErrorWarning() {fgWarnError = kFALSE;} // disable warning and error info
  static void    SetForceCleanDataOnly() {fgCleanDataOnly = kTRUE;} // clean data only
  static void    AllowCorruptedData() {fgCleanDataOnly = kFALSE;} // accept corrupted data

  static void    SetExtraWordsFix() {fgExtraSkip = kTRUE;} // extra skip of 24 32-bit words 
  static void    SetSkipCDH() {fgSkipCDH = kTRUE;} // skip of 8 32-bit words 
  static void    SetDumpHead(Int_t iv) {fgDumpHead = iv;}
  static void    DisableStackNumberChecker() {fgStackNumberChecker = kFALSE;}  // set false to cleanroom data 
  static void    DisableStackLinkNumberChecker() {fgStackLinkNumberChecker = kFALSE;}  
  static void    DisableSkipData() {fgSkipData = kFALSE;} // keep reading next words even previous words were corrupted - debugging purpose  
  static void    SetDumpingEnable() {fDumpingEnable = kTRUE;} 
  static void    SetDumpingMCM(Int_t sm, Int_t stack, Int_t layer, Int_t rob, Int_t mcm) {fDumpingSM = sm; fDumpingStack = stack; fDumpingLayer = layer; fDumpingROB = rob; fDumpingMCM = mcm;}

  // this is a temporary solution!
  // baseline should come with the HC header word 2 (count from 0!)
  static void    SetSubtractBaseline(Int_t baseline) {fgCommonAdditive = baseline;}
  Int_t          GetCommonAdditive() const {return fgCommonAdditive;}           // return the common additive

  //--------------------------------------------------------
  // Decoding functions
  //--------------------------------------------------------

  void DecodeSMInfo(const UInt_t *word, struct AliTRDrawSM *sm) const ;
  const char *DumpSMInfo(const struct AliTRDrawSM *sm);
  void DecodeStackInfo(const UInt_t *word, struct AliTRDrawStack *st) const;
  void DecodeStackHeader(const UInt_t *word, struct AliTRDrawStack *st, Int_t iword) const;
  const char *DumpStackInfo(const struct AliTRDrawStack *st);
  Bool_t DecodeHCwordH0(const UInt_t *word, struct AliTRDrawHC *hc) const;
  Bool_t DecodeHCwordH1(const UInt_t *word, struct AliTRDrawHC *hc) const;
  const char *DumpHCinfoH0(const struct AliTRDrawHC *hc);
  const char *DumpHCinfoH1(const struct AliTRDrawHC *hc);
  void DecodeMCMheader(const UInt_t *word, struct AliTRDrawMCM *mcm) const;
  UInt_t GetMCMadcMask(const UInt_t *word, struct AliTRDrawMCM *mcm) const;
  void DecodeMask(const UInt_t *word, struct AliTRDrawMCM *mcm) const;
  void MCMADCwordsWithTbins(UInt_t fTbins, struct AliTRDrawMCM *mcm) const;
  const char *DumpMCMinfo(const struct AliTRDrawMCM *mcm);
  const char *DumpMCMadcMask(const struct AliTRDrawMCM *mcm);


 protected:

  Bool_t InitBuffer(void *buffer, UInt_t length); // init the buffer - called by DecodeSM(void*, UInt_t)
  Bool_t DumpWords(UInt_t *px, UInt_t iw, UInt_t marker = 0); // dump some words onto the screen

  void   SwapOnEndian();         // swap if endian is BIG
  Bool_t SkipWords(UInt_t iw);   // skip number of words
  Bool_t DecodeGTUheader();      // decode data in GTU header
  Bool_t DecodeTracklets();      // decode tracklets
  Bool_t DecodeHC();             // decode data in HC

  Bool_t DecodeADC();            // decode 10 ADC words

  Bool_t DecodeHCheader();       // decode HC  header
  Bool_t SeekEndOfData();        // go to next end of raw data marker (actually 1 word after)
  Bool_t SkipMCMdata(UInt_t iw);  // skip this mcm data due to mcm header corruption
  Bool_t SeekNextMCMheader();    // go to next mcm header
  Bool_t DecodeMCMheader();      // decode mcm header

  Bool_t IsRowValid();       // check if row within the range
  Bool_t IsHCheaderOK();     // check if current hc header data make sense
  Bool_t IsMCMheaderOK();    // check if current mcm header data make sense
    
  void   ResetCounters();    // reset some counters
  void   ResetIterators();   // needed for Next()
  void   ResetPerSM();       // reset every SM 
  void   ResetPerStack();    // reset every Stack 
  void   ResetPerHC();       // reset every HC 
  void   ResetPerMCM();      // reset every MCM  
  void   ResetPerADC();      // reset every ADC  
  void   ResetMemory();      // reset all data members

  AliTRDrawStream(const AliTRDrawStream& st);
  AliTRDrawStream &operator=(const AliTRDrawStream &);

  // ----------------- DATA MEMBERS START

  struct AliTRDrawSM       fSM;    // one SM per buffer
  struct AliTRDrawStack   *fStack; //! pointer to the current stack
  struct AliTRDrawHC      *fHC;    //! current HC
  struct AliTRDrawMCM     *fMCM;   //! current MCM
  struct AliTRDrawADC     *fADC;   //! current ADC
  
  UInt_t *fpPos;   // current position in the buffer
  UInt_t *fpBegin; // begin - pointer to the buffer word 0
  UInt_t *fpEnd;   // end of the buffer

  UInt_t  fWordLength;  // length of the buffer in 32bit words

  Int_t   fStackNumber;     // current stack number
  Int_t   fStackLinkNumber; // current link in the stack

  Int_t   fhcMCMcounter;    // mcm counter inside single hc - used in Next()  
  Int_t   fmcmADCcounter;   // adc counrer inside single adc - used in Next() 

  Int_t   fLinkTrackletCounter; // count the tracklets in the current HC
  Int_t   fEndOfTrackletCount;  // count link by link (hc by hc) used for debug
  Int_t   fNWordsCounter;       // counts words of given hc having link monitor error

  UInt_t  fMaskADCword; // temp mask when decoding adcs
  UInt_t  fTbinADC;     // temp adc 
  Int_t   fDecodedADCs; // counter of decoded adcs [mj] do we need?

  UInt_t  fEventCounter;     // stores the valid/current MCM event counter
  UInt_t  fLastEventCounter; // last known event counter of MCM

  Bool_t  fSharedPadsOn; // do we want to output shared pads - default is off
  Int_t   fMaxADCgeom;   // maximum ADC channels per mcm

  Bool_t  fBufferRead;

  AliTRDgeometry *fGeometry;  //! TRD geometry
  AliRawReader   *fRawReader; //! raw reader    

  AliTRDfeeParam      *fTRDfeeParam; // pointer to the fee params

  // STATIC 

  static Bool_t fgExtraSkip; // whether we should skip the leading 24 words
  static Bool_t fgSkipCDH; // whether we should skip CDH (8 words)
  static Bool_t fgWarnError; // no errors no warnings
  static Bool_t fgCleanDataOnly; // release only clean events = no errors
  static Bool_t fgDebugFlag; // allow debugging info
  static Bool_t fgEnableMemoryReset; // allow memory reset
  static Bool_t fgStackNumberChecker; // decide if we check stack number insanity - set false to cleanroom data
  static Bool_t fgStackLinkNumberChecker; // decide if we check stack link number insanity - debuging purpose
  static Bool_t fgSkipData; // decide if we skip corrupted data of given HC
  static Int_t fgDumpHead; // number of words to dump (from the start of the buffer) on each Init
  static Int_t  fgEmptySignals[30]; // empty signals in case of ADC pointer = NULL
  static Short_t  fgMCMordering[16]; // mcm number odering for mcm header corruption check
  static Short_t  fgROBordering[16]; // mcm number odering for mcm header corruption check
  static Int_t  fgLastHC; 
  static Int_t  fgLastROB; 
  static Int_t  fgLastIndex; 

  static Bool_t fDumpingEnable; 

  static Int_t  fDumpingSM;
  static Int_t  fDumpingStack;
  static Int_t  fDumpingLayer;
  static Int_t  fDumpingROB;
  static Int_t  fDumpingMCM;

  // this is a temporary solution!
  // baseline should come with the HC header word 2 (count from 0!)
  static Int_t   fgCommonAdditive; // common additive  - should be decoded! from HC word2

  // ----------------- DATA MEMBERS STOP

  enum ETRDzRawStreamError 
    {
      kDecodeStackInfo          = 1 //
      , kMissingData	        = 2 //
      , kLinkDataMissing        = 3 //
      , kHCdataMissing	        = 4 //
      , kTrackletOverflow       = 5 //
      , kEOTrackeltsMissing     = 6 //
      , kWrongPadrow	        = 7 //
      , kMCMheaderCorrupted     = 8 //
      , kWrongMCMorROB	        = 9 //
      , kMCMeventMissmatch      = 10 //
      , kMCMADCMaskMissing      = 11 //
      , kHCHeaderCorrupt        = 12 //
      , kHCHeaderWrongStack     = 13 //
      , kHCHeaderWrongLayer     = 14 //
      , kHCHeaderWrongSide      = 15 //
      , kHCHeaderWrongSM        = 16 //
      , kHCHeaderWrongDet       = 17 //
      , kHCHeaderWrongROC       = 18 //
      , kHCWordMissing	        = 19 //
      , kMCMdataMissing	        = 20 //
      , kMCMoverflow	        = 21 //
      , kADCdataMissing	        = 22 //
      , kADCmaskMissmatch       = 23 //
      , kWrongPadcolumn	        = 24 //
    };			       

  ClassDef(AliTRDrawStream, 0)
}; //clas def end

#endif
