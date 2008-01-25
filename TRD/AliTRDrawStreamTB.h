#ifndef ALITRDRAWSTREAMTB_H
#define ALITRDRAWSTREAMTB_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

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

// definitions in AliTRDrawStreamBase.h:
/* #define TRD_MAX_TBINS 30 */
/* #define TRD_MAX_ADC   21 */
/* #define TRD_MAX_MCM   4 * 16 */

//class AliTRDrawStreamTB : public TObject
class AliTRDrawStreamTB : public AliTRDrawStreamBase
{ // class def begin

 public:
    
  //--------------------------------------------------------

  // THE STRUCTURES
  
  //--------------------------------------------------------

  struct AliTRDrawADC
  {//adc struct

    UInt_t             *fPos; //! position of ADC 1st word in the buffer
    Short_t             fADCnumber; // number of the ADC 0 .. 20
    Short_t             fCOL; // column - row from MCM
    Int_t               fSignals[TRD_MAX_TBINS]; // signals for this adc
    Bool_t              fIsShared; // is pad chared between MCMs
    Short_t             fCorrupted; // is adc word corrupted

    AliTRDrawADC()
      : fPos(0)
      , fADCnumber(0)
      , fCOL(0)
      , fSignals()
      , fIsShared(kFALSE)
      , fCorrupted(0)
    {
      // default constructor
      ;
    }

    AliTRDrawADC(const AliTRDrawADC& p): 
        fPos(p.fPos)
      , fADCnumber(p.fADCnumber)
      , fCOL(p.fCOL)
      , fSignals()
      , fIsShared(p.fIsShared)
      , fCorrupted(p.fCorrupted)
    {
      // copy constructor
      ; 
    }

    AliTRDrawADC &operator=(const AliTRDrawADC &) 
    {
      // assignment operator
      // not implemented
      return *this;
    }


  };
  
  //--------------------------------------------------------

  struct AliTRDrawMCM
  { // mcm struct
    Int_t               fROB; // ROB number
    Int_t               fMCM; // MCM number
    Int_t               fROW; //filed during decoding!
      
    UInt_t              fEvCounter; // MCM event counter
    UInt_t              fADCMask;   // ADC mask
    UInt_t              fADCMaskWord; // word with ADC mask in
    UInt_t              fADCchannel[30]; // channels to be decoded accrording to ADC mask
      
    Int_t               fADCindex; // index of current ADC (comment: 1 ADC is 1 pad)
    Int_t               fADCmax;   // number of ADCs fired
    Int_t               fMCMADCWords; // mcm words to expect
    Int_t               fSingleADCwords; // n of words per ADC
      
    Int_t               fCorrupted; // is mcm word corrupted
      
    UInt_t              fErrorCounter; // count the mcm header errors
    UInt_t              fMaskErrorCounter; // count the adc mask errors
      
    UInt_t             *fPos; //! position of mcm header in the buffer
    UInt_t             *fAdcDataPos; //! start of ADC data for this mcm

    Int_t               fADCcounter; // count the adcs decoded
    AliTRDrawADC        fADCs[TRD_MAX_ADC]; // 21 adcs
      
    AliTRDrawMCM()
      : fROB(-1)
      , fMCM(-1)
      , fROW(-1)
      , fEvCounter(0)
      , fADCMask(0)
      , fADCMaskWord(0)
      , fADCchannel()      
      , fADCindex(0)
      , fADCmax(0)
      , fMCMADCWords(0)      
      , fSingleADCwords(0)
      , fCorrupted(0)      
      , fErrorCounter(0)
      , fMaskErrorCounter(0)
      , fPos(0)
      , fAdcDataPos(0)
      , fADCcounter(0)
      , fADCs()
    {
      // default constructor
      /* the following violates coding conventions */
      //for (Int_t i = 0; i < 30; i++) fADCchannel[i] = 0;
    }

    AliTRDrawMCM(const AliTRDrawMCM & p):
        fROB(p.fROB)
      , fMCM(p.fMCM)
      , fROW(p.fROW)
      , fEvCounter(p.fEvCounter)
      , fADCMask(p.fADCMask)
      , fADCMaskWord(p.fADCMaskWord)
      , fADCchannel()      
      , fADCindex(p.fADCindex)
      , fADCmax(p.fADCmax)
      , fMCMADCWords(p.fMCMADCWords)      
      , fSingleADCwords(p.fSingleADCwords)
      , fCorrupted(p.fCorrupted)      
      , fErrorCounter(p.fErrorCounter)
      , fMaskErrorCounter(p.fMaskErrorCounter)
      , fPos(p.fPos)
      , fAdcDataPos(p.fAdcDataPos)
      , fADCcounter(p.fADCcounter)
      , fADCs()
    {
      // copy constructor
      ;
    }

    AliTRDrawMCM &operator=(const AliTRDrawMCM &)
    {
      // assignment operator
      // not implemented
      return *this;
    }

  };

  //--------------------------------------------------------

  struct AliTRDrawHC
  { // hc struct
    // global
    Int_t               fCorrupted;         // Zero if not corrupted

    UInt_t              fH0ErrorCounter;  // Count the H0 word errors
    UInt_t              fH1ErrorCounter;  // Count the H1 word errors
      
    // word 0
    Int_t               fSpecialRawV;       // Raw data version
    Int_t               fRawVMajor;         // Raw data version
    Int_t               fRawVMinor;         // Raw data version
    Int_t               fNExtraWords;       // N extra HC header words
    Int_t               fDCSboard;          // DCS board number
    Int_t               fSM;                // Super Module number
    Int_t               fStack;             // Stack number (some people might call it a chamber)
    Int_t               fLayer;             // Layer number (some people might call it a plane)
    Int_t               fSide;              // Side of HC
      
    //word 1
    Int_t               fTimeBins;          // N of t bins
    UInt_t              fBunchCrossCounter; // Bunch crossing counter
    UInt_t              fPreTriggerCounter; // Pre Trigger counter
    UInt_t              fPreTriggerPhase;   // Pre Trigger phase
      
    UInt_t             *fPos[2];           //! position of the header words in buffer
      
    Int_t               fDET; // filled while decoding
    Int_t               fROC; // filled while decoding
    Int_t               fRowMax; // filled while decoding
    Int_t               fColMax; // filled while decoding

    //AliTRDrawMCM        fMCMs[4][16]; // 4 ROBS 16 each - that is max!
    Int_t               fMCMmax; // number of mcm found
    AliTRDrawMCM        fMCMs[TRD_MAX_MCM]; // 4 ROBS 16 each - that is max!

    AliTRDrawHC()
      : fCorrupted(0)
      , fH0ErrorCounter(0)
      , fH1ErrorCounter(0)
      , fSpecialRawV(0)
      , fRawVMajor(0)
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
      , fPos()
      , fDET(-1)
      , fROC(-1)
      , fRowMax(-1)
      , fColMax(-1)
      , fMCMmax(0)
      , fMCMs()
    {
      // default constructor hc info 
      /* the following violates coding conventions */
      //fPos[0] = fPos[0] = 0;
    }

    AliTRDrawHC(const AliTRDrawHC & p):
        fCorrupted(p.fCorrupted)
      , fH0ErrorCounter(p.fH0ErrorCounter)
      , fH1ErrorCounter(p.fH1ErrorCounter)
      , fSpecialRawV(p.fSpecialRawV)
      , fRawVMajor(p.fRawVMajor)
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
      , fPos()
      , fDET(p.fDET)
      , fROC(p.fROC)
      , fRowMax(p.fRowMax)
      , fColMax(p.fColMax)
      , fMCMmax(p.fMCMmax)
      , fMCMs()
    {
      // copy constructor
      ;
    }

    AliTRDrawHC &operator=(const AliTRDrawHC &)
    {
      // assignment operator
      // not implemented
      return *this;
    }

  };

  //--------------------------------------------------------
    
  struct AliTRDrawStack
  {
    UInt_t           fHeaderSize;         // header size of the stack info
    Bool_t           fLinksActive[12];    // data links active - 1 per half chamber
    Bool_t           fTrackletDecode[12]; // book keeping while decoding - set false after decoding
    Bool_t           fHCDecode[12];       // book keeping while decoding - set false after HC header decoding
    Int_t            fActiveLinks;        // number of active links
    UInt_t          *fPos;                //! position in the buffer
	            
    AliTRDrawHC      fHalfChambers[12];   // 6 chambers in a stack
      
    AliTRDrawStack()
      : fHeaderSize(0)
      , fLinksActive()
      , fTrackletDecode() //book keeping while decoding - set false after decoding
      , fHCDecode() //book keeping while decoding - set false after HC header decoding
      , fActiveLinks(0)
      , fPos(0)
      , fHalfChambers()
    {
      // default constructor
    }      

    AliTRDrawStack(const AliTRDrawStack & p):
        fHeaderSize(p.fHeaderSize)
      , fLinksActive()
      , fTrackletDecode() //book keeping while decoding - set false after decoding
      , fHCDecode() //book keeping while decoding - set false after HC header decoding
      , fActiveLinks(p.fActiveLinks)
      , fPos(p.fPos)
      , fHalfChambers()
    {
      // copy constructor
      ;
    }

    AliTRDrawStack &operator=(const AliTRDrawStack &)
    {
      // assignment operator
      // not implemented
      return *this;
    }

  };

  /* the following violates coding conventions */
/*       for (Int_t i = 0; i < 12; i++) */
/* 	{ */
/* 	  fLinksActive[i] = fTrackletDecode[i] = fHCDecode[i] = kFALSE; */
/* 	} */

  //--------------------------------------------------------

  struct AliTRDrawSM
  {
    UInt_t            fHeaderSize;            // size of the header in words
    Bool_t            fTrackletEnable;        // tracklet enable bit
    Bool_t            fStackActive[5];        // map of active/expected stacks
    Int_t             fActiveStacks;          // number of active stacks
    Int_t             fCorrupted;             // is sm info corrupted
    Int_t             fNexpectedHalfChambers; // number of half chambers to be read out
    Bool_t            fClean;                 // true if everything went OK - false is some error occured
    UInt_t           *fPos;                   //! location of the sm info - should be the first word (after CDH if not DDL buffer)

    AliTRDrawStack    fStacks[5];             // we do have only five stacks ;)

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
	// coding rule violation in the next line
	//for (Int_t i = 0; i < 5; i++) fStackActive[i] = kFALSE;
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
      ;
    }

    AliTRDrawSM &operator=(const AliTRDrawSM &)
    {
      // assignment operator
      // not implemented
      return *this;
    }

  };
  
  //--------------------------------------------------------
     
  AliTRDrawStreamTB();
  AliTRDrawStreamTB(AliRawReader *rawReader);
  virtual ~AliTRDrawStreamTB();

  //--------------------------------------------------------

  virtual Bool_t       Next();              // Read the next data
  virtual Int_t NextChamber(AliTRDdigitsManager *man); // read next chamber data
  virtual Bool_t       Init();              // Init some internal variables

  Bool_t   SetRawVersion(Int_t fraw); // set the raw version - used for backward compat.
  
  Bool_t   IsCurrentPadShared() const {return fADC->fIsShared;}          // is current pad shared between mcms
  void     SetSharedPadReadout(Bool_t fv) {fSharedPadsOn = fv;}          //  set the flag on if the reader should return the shared pads
  Bool_t   IsDataZeroSuppressed() const {return (fHC->fRawVMajor > 2) ? kTRUE : kFALSE;} // check the version and tell if ZS is on
  
  Bool_t   DecodeSM(void *buffer, UInt_t length); //decode a buffer
  Int_t    DecodeSM(); // used with raw reader
  Int_t    DecodeSM(AliRawReader *reader); // used with raw reader
	   
  Bool_t   SetReader(AliRawReader *reader); // set the raw reader to use
    	   
  Bool_t    IsTrackletEnableBitSet() const {return fSM.fTrackletEnable;} // get status of the tracklets enable bit
  Bool_t    IsStackActive(Int_t is) const {return fSM.fStackActive[is];} // get status of a stack
  Int_t     GetNofActiveStacks() const {return fSM.fActiveStacks;} // get number of active stacks
  UInt_t   *GetSMstreamPosition() const {return fSM.fPos;} // get position of the SM index word in the buffer
	    
  Bool_t    IsSMbufferClean() const {return fSM.fClean;} // is data clean
    
  Bool_t    IsLinkActiveInStack(Int_t is, Int_t il) const {return fSM.fStacks[is].fLinksActive[il];} // check whether the link is active
  Int_t     GetActiveLinksInStack(Int_t is) const {return fSM.fStacks[is].fActiveLinks;} //get active links in a stack
    	    
  Int_t     GetSpecialRawVersion() const {return fHC ? fHC->fSpecialRawV : -1;} // return special raw version
  Int_t     GetMajorRawVersion() const {return fHC ? fHC->fRawVMajor : -1;} // major raw version getter
  Int_t     GetRawVersion() const {return fHC ? fHC->fRawVMajor : -1;} // compatibility see funtion above
  Int_t     GetMinorRawVersion() const {return fHC ? fHC->fRawVMinor : -1;} // minor raw version

  Int_t     GetSM() const {return fHC ? fHC->fSM : -1;} //  Position of CURRENT half chamber in full TRD
  Int_t     GetLayer() const {return fHC ? fHC->fLayer : -1;} //  PLANE = Position of CURRENT half chamber in full TRD
  Int_t     GetStack() const {return fHC ? fHC->fStack : -1;} //  CHAMBER = Position of CURRENT half chamber in full TRD
  Int_t     GetSide() const {return fHC ? fHC->fSide : -1;} // get side
  Int_t     GetDCS() const { return fHC ? fHC->fDCSboard : -1;} //  DCS board number read from data (HC header)
  //Int_t     GetDCSboard() {return fHC ? fHC->fDCSboard : -1;} //  DCS board number read from data (HC header)
  Int_t     GetROC() const { return fHC ? fHC->fROC : -1;} //  Position of CURRENT half chamber in full TRD
  Int_t     GetNumberOfTimeBins() const { return fHC ? fHC->fTimeBins : 0;} // Get Ntime bins
  UInt_t    GetBunchCrossCounter() const {return fHC ? fHC->fBunchCrossCounter : 0;} // get bunch cross counter
  UInt_t    GetPreTriggerCounter() const {return fHC ? fHC->fPreTriggerCounter : 0;} // get pre trigger info
  UInt_t    GetPreTriggerPhase() const {return fHC ? fHC->fPreTriggerPhase : 0;} // get trigger phase

  Int_t     GetRow() const {return fMCM ? fMCM->fROW : -1;} // get current row number
  Int_t     GetCol() const {return fADC ? fADC->fCOL : -1;} // get current column number
  Int_t     GetRowMax() const { return fHC ? fHC->fRowMax : -1;} // Get maximum rows in the current HC
  Int_t     GetColMax() const { return fHC ? fHC->fColMax : -1;} // Get maximum cols in the current HC
  // compatibility
  Int_t     GetMaxRow() const { return fHC ? fHC->fRowMax : -1;} // Get maximum rows in the current HC
  Int_t     GetMaxCol() const { return fHC ? fHC->fColMax : -1;} // Get maximum cols in the current HC

  UInt_t    GetHCword0() const {return fHC ? *fHC->fPos[0] : 0;} // get the HC word 0
  UInt_t    GetHCword1() const {return fHC ? *fHC->fPos[1] : 0;} // get the HC word 1
	    
  Int_t     GetDET() const {return fHC ? fHC->fDET : -1;} // get current det number
  Int_t     GetDet() const {return fHC ? fHC->fDET : -1;} // get current det number
    	    
  Int_t     GetROB() const {return fMCM ? fMCM->fROB : -1;} // get current ROB number
  Int_t     GetMCM() const {return fMCM ? fMCM->fMCM : -1;} // get current MCM number
  Int_t     GetEventNumber() const { return fMCM->fEvCounter;} //  MCM Event number and position of current MCM on TRD chamber

  Int_t     IsMCMcorrupted() const {return fMCM ? fMCM->fCorrupted : -1;} // is current MCM header corrupted

  Int_t    *GetSignals() const { return fADC ? fADC->fSignals : (Int_t *)fgEmptySignals;}//Signals in the n-time bins from Data Word
  Int_t     GetADC() const { return fADC ? fADC->fADCnumber : -1 ;}        //  MCM ADC channel and Time Bin of word 1
  Int_t     GetTimeBin() const { return 0;}                                //  MCM ADC channel and Time Bin of word 1
  
  //----------------------------------------------------------
 
  static void    SetNoDebug() {fgDebugFlag = kFALSE;} // allow debug info
  static void    SetNoErrorWarning() {fgWarnError = kFALSE;} // disable warning and error info
  static void    SetForceCleanDataOnly() {fgCleanDataOnly = kTRUE;} // clean data only
  static void    AllowCorruptedData() {fgCleanDataOnly = kFALSE;} // accept corrupted data

  static void    SetExtraWordsFix() {fgExtraSkip = kTRUE;} // extra skip of 24 32-bit words 
  static void    SetSkipCDH() {fgSkipCDH = kTRUE;} // skip of 8 32-bit words 
  void           EnableDebug(TTreeSRedirector *debugStream = 0); // enable the debug stream - use the paramter if non zero
  static void    EnableDebugStream() {fgDebugStreamFlag = kTRUE;} //global enable of the dbug stream
  static void    DeleteDebugStream(); // helper function to delete the debug streamer
  static void    SetDumpHead(UInt_t iv) {fgDumpHead = iv;}
  static void    DisableStackNumberChecker() {fgStackNumberChecker = kFALSE;}  // set false to cleanroom data 

  // this is a temporary solution!
  // baseline should come with the HC header word 2 (count from 0!)
  static void    SetSubtractBaseline(Int_t baseline) {fgCommonAdditive = baseline;}
  Int_t          GetCommonAdditive() const {return fgCommonAdditive;}           // return the common additive

  void    DumpErrorCount();
  void    ReSetStreamEventCounter(Int_t ival = 0) {fgStreamEventCounter = ival;} // reset the event counter for the debug streamer

  //--------------------------------------------------------
  // Decoding functions
  //--------------------------------------------------------

  void DecodeSMInfo(const UInt_t *word, struct AliTRDrawSM *sm) const ;
  const char *DumpSMInfo(const struct AliTRDrawSM *sm);
  void DecodeStackInfo(const UInt_t *word, struct AliTRDrawStack *st) const;
  const char *DumpStackInfo(const struct AliTRDrawStack *st);
  void DecodeHCwordH0(const UInt_t *word, struct AliTRDrawHC *hc) const;
  void DecodeHCwordH1(const UInt_t *word, struct AliTRDrawHC *hc) const;
  const char *DumpHCinfoH0(const struct AliTRDrawHC *hc);
  const char *DumpHCinfoH1(const struct AliTRDrawHC *hc);
  void DecodeMCMheader(const UInt_t *word, struct AliTRDrawMCM *mcm) const;
  UInt_t GetMCMadcMask(const UInt_t *word, struct AliTRDrawMCM *mcm) const;
  void DecodeMask(const UInt_t *word, struct AliTRDrawMCM *mcm) const;
  void MCMADCwordsWithTbins(UInt_t fTbins, struct AliTRDrawMCM *mcm) const;
  const char *DumpMCMinfo(const struct AliTRDrawMCM *mcm);
  const char *DumpMCMadcMask(const struct AliTRDrawMCM *mcm);

 protected:

 private:

  Bool_t InitBuffer(void *buffer, UInt_t length); // init the buffer - called by DecodeSM(void*, UInt_t)
  Bool_t DumpWords(UInt_t *px, UInt_t iw, UInt_t marker = 0); // dump some words onto the screen

  Int_t  NextBuffer(); // go and init next buffer if available - check the implementation file for return values

  void   SwapOnEndian();       // swap if endian is BIG
  Bool_t SkipWords(UInt_t iw); // skip number of words
  Bool_t DecodeTracklets();    // decode tracklets
  Bool_t DecodeHC();           // decode data in HC

  Bool_t DecodeADC();          // decode 10 ADC words

  Bool_t DecodeHCheader();       // decode HC  header
  Bool_t SeekEndOfData();      // go to next end of raw data marker (actually 1 word after)
  Bool_t SeekNextMCMheader();    // go to next mcm header
  Bool_t DecodeMCMheader();      // decode mcm header

  Bool_t IsRowValid(); // check if row within the range
  Bool_t IsHCheaderOK(); // check if current hc header data make sense
  Bool_t IsMCMheaderOK(); // check if current mcm header data make sense
    
  void   ResetCounters(); // reset some counters
  void   ResetIterators(); // needed for Next()
  AliTRDrawStreamTB(const AliTRDrawStreamTB& st);
  AliTRDrawStreamTB &operator=(const AliTRDrawStreamTB &);

  // ----------------- DATA MEMBERS START

  struct AliTRDrawSM       fSM;    // one SM per buffer
  struct AliTRDrawStack   *fStack; //! pointer to the current stack
  struct AliTRDrawHC      *fHC;    //! current HC
  struct AliTRDrawMCM     *fMCM;   //! current MCM
  struct AliTRDrawADC     *fADC;   //! current ADC
  
  UInt_t *fpPos;   // current position in the buffer
  UInt_t *fpBegin; // begin - pointer to the buffer word 0
  UInt_t *fpEnd;   // end of the buffer

  UInt_t  fWordLength; // length of the buffer in 32bit words

  Int_t   fStackNumber;     // current stack number
  Int_t   fStackLinkNumber; // current link in the stack

  Int_t   fhcMCMcounter; // mcm counter inside single hc - used in Next()
  Int_t   fmcmADCcounter; // adc counrer inside single adc - used in Next()

  Int_t   fLinkTrackletCounter; // count the tracklets in the current HC
  Int_t   fEndOfTrackletCount;  // count link by link (hc by hc) used for debug

  UInt_t  fMaskADCword; // temp mask when decoding adcs
  UInt_t  fTbinADC;     // temp adc 
  Int_t   fDecodedADCs; // counter of decoded adcs

  UInt_t  fEventCounter; // stores the valid/current MCM event counter
  UInt_t  fLastEventCounter; // last known event counter of MCM

  Bool_t  fSharedPadsOn; // do we want to output shared pads - default is off
  Int_t   fMaxADCgeom; // maximum ADC channels per mcm

  AliTRDgeometry *fGeometry; //! TRD geometry
  AliRawReader   *fRawReader; //! raw reader    

  AliTRDfeeParam      *fTRDfeeParam; // pointer to the fee params

  Bool_t               fDebugStreamOwned; // created in this case by this - allowed to delete?

  // STATIC 

  static Bool_t fgExtraSkip; // whether we should skip the leading 24 words
  static Bool_t fgSkipCDH; // whether we should skip CDH (8 words)
  static Bool_t fgWarnError; // no errors no warnings
  static Bool_t fgCleanDataOnly; // release only clean events = no errors
  static Bool_t fgDebugFlag; // allow debugging info
  static Bool_t fgDebugStreamFlag; // set on debug streamer
  static Bool_t fgStackNumberChecker; // decide if we check stack number insanity - set false to cleanroom data
  static TTreeSRedirector *fgDebugStreamer; //!Debug streamer
  static UInt_t fgStreamEventCounter; // event counter for debug streamer
  static UInt_t fgDumpHead; // number of words to dump (from the start of the buffer) on each Init
  static Int_t  fgEmptySignals[30]; // empty signals in case of ADC pointer = NULL


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

  ClassDef(AliTRDrawStreamTB, 0)
}; //clas def end

#endif
























