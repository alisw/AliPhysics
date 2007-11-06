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

#include <TObject.h>

class AliTRDgeometry;
class AliRawReader;
class AliTRDdigitsManager;

namespace AliTRDrawDataUtilsTB
{
  //--------------------------------------------------------
  // Some constants:
  //static const UInt_t kEndoftrackletmarker = 0xAAAAAAAA; /*This marks the end of tracklet data words*/
  //static const UInt_t kEndOfTrackletMarker = 0x10001000; /*This marks the end of tracklet data words*/
  static const UInt_t kEndOfTrackletMarker = 0xaaaaaaaa; /*This marks the end of tracklet data words*/
  //static const UInt_t kEndofrawdatamarker  = 0x00000000; /*This marks the end of half-chamber-data*/
  static const UInt_t kEndOfRawDataMarker  = 0x00000000; /*This marks the end of half-chamber-data*/
  static const UInt_t kSizeWord = sizeof(UInt_t); // size of a word in bytes

  //--------------------------------------------------------
  // SM index word masks
  static const UInt_t skSMheaderSizeBits = 0xffff0000;
  static const UInt_t skTrackletsEnableBits = 0x1 << 5;
  static const UInt_t skStackMasksBits = 0x1f;

  static UInt_t GetSMheaderSize(UInt_t *word)
  {
    return (*word & skSMheaderSizeBits) >> 16;
  }

  static UInt_t GetTrackletEnable(UInt_t *word)
  {
    return (*word & skTrackletsEnableBits) >> 5;    
  }
  
  static UInt_t GetStackMask(UInt_t *word)
  {
    return (*word & skStackMasksBits);    
  }

/*   static UInt_t GetStackStatus(UInt_t *word, Int_t istack) */
/*   { */
/*     return (*word & (0x1 << istack)); */
/*   } */

  //--------------------------------------------------------
  struct AliTRDrawSMinfo
  {
    UInt_t fHeaderSize;
    Bool_t fTrackletEnable;
    Bool_t fStackActive[5];
    Int_t  fActiveStacks;

    void Decode(UInt_t *word)
    {
      fHeaderSize = GetSMheaderSize(word);

      if (GetTrackletEnable(word) > 0)
	fTrackletEnable = kTRUE;
      else
	fTrackletEnable = kFALSE;

      UInt_t stackMask = GetStackMask(word);
      fActiveStacks = 0;
      for (Int_t i = 0; i < 5; i++)
	{
	  //if (stackMask & (0x1 << i) > 0)
	  if ((stackMask >> i) & 0x1 > 0)
	    {
	      fStackActive[i] = kTRUE;
	      fActiveStacks++;
	    }
	  else
	    {
	      fStackActive[i] = kFALSE;
	    }
	}
    }
    
    void Dump()
    {
      printf("[I] SM Info is : Hsize %d TrackletEnable %d Stacks %d %d %d %d %d\n",
	     fHeaderSize, fTrackletEnable, 
	     fStackActive[0], fStackActive[1], fStackActive[2], 
	     fStackActive[3], fStackActive[4]); 
      
    }
  };

  //--------------------------------------------------------
  // Stack word masks
  static const UInt_t skStackHeaderSizeBits = 0xffff0000;
  static const UInt_t skStackLinkMaskBits = 0xfff;

  static UInt_t GetStackHeaderSize(UInt_t *word)
    {
      return (*word & skStackHeaderSizeBits) >> 16;
    }

  static UInt_t GetStackLinkWord(UInt_t *word)
    {
      return (*word & skStackLinkMaskBits);
    }

  static UInt_t GetStackHeaderSizeBug(UInt_t *word)
    {
      return (*word & 0x0000ffff);
    }

  static UInt_t GetStackLinkWordBug(UInt_t *word)
    {
      return (*word & 0x0fff) >> 16;
    }

  struct AliTRDrawStackInfo
  {
    UInt_t fHeaderSize;
    Bool_t fLinksActive[12];
    Int_t  fActiveLinks;

    void Decode(UInt_t *word)
    {
      fHeaderSize = GetStackHeaderSize(word);

      UInt_t linkMask = GetStackLinkWord(word);
      fActiveLinks = 0;
      for (Int_t i = 0; i < 12; i++)
	{
	  if (linkMask & (0x1 << i) > 0)
	    {
	      fLinksActive[i] = kTRUE;
	      fActiveLinks++;
	    }
	  else
	    {
	      fLinksActive[i] = kFALSE;
	    }
	}
    }

    void DecodeBug(UInt_t *word)
    {
      fHeaderSize = GetStackHeaderSizeBug(word);

      UInt_t linkMask = GetStackLinkWordBug(word);
      for (Int_t i = 0; i < 12; i++)
	{
	  if (linkMask & (0x1 << i) > 0)
	    fLinksActive[i] = kTRUE;
	  else
	    fLinksActive[i] = kFALSE;
	}
    }
    
    void Dump()
      {
	printf("[I] Stack Info is : Hsize %d Links Active %d %d %d %d %d %d %d %d %d %d %d %d\n",
	       fHeaderSize,
	       fLinksActive[0], fLinksActive[1], fLinksActive[2], fLinksActive[3], 
	       fLinksActive[4], fLinksActive[5], fLinksActive[6], fLinksActive[7], 
	       fLinksActive[8], fLinksActive[9], fLinksActive[10], fLinksActive[11]);
      }
  };

  //--------------------------------------------------------
  // HC word masks
  // word 0
  static const UInt_t skSpecialRawVBits = 0x1 << 31;
  static const UInt_t skRawVMajorBits = 0x7f << 24;
  static const UInt_t skRawVMinorBits = 0x7f << 17;
  static const UInt_t skHCExtraWordsBits = 0x7 << 14;
  static const UInt_t skHCSMBits = 0x1f << 9;
  static const UInt_t skHCLayerBits = 0x7 << 6;
  static const UInt_t skHCStackBits = 0x7 << 3;
  static const UInt_t skHCSideBits = 0x1 << 2;

  static const UInt_t skDCSBoardBits = 0xfff << 20;

  static UInt_t GetSpecialRawV(UInt_t *word)
    {
      return (*word & skSpecialRawVBits) >> 31;
    }

  static UInt_t GetRawVMajor(UInt_t *word)
    {
      return (*word & skRawVMajorBits) >> 24;
    }

  static UInt_t GetRawVMinor(UInt_t *word)
    {
      return (*word & skRawVMinorBits) >> 31;
    }

  static UInt_t GetDCSBoard(UInt_t *word)
    {
      return (*word & skDCSBoardBits) >> 20;
    }

  static UInt_t GetSM(UInt_t *word)
    {
      return (*word & skHCSMBits) >> 9;
    }

  static UInt_t GetStack(UInt_t *word)
    {
      return (*word & skHCStackBits) >> 3;
    }

  static UInt_t GetLayer(UInt_t *word)
    {
      return (*word & skHCLayerBits) >> 6;
    }

  static UInt_t GetSide(UInt_t *word)
    {
      return (*word & skHCSideBits) >> 2;
    }

  static UInt_t GetExtraWords(UInt_t *word)
    {
      return (*word & skHCExtraWordsBits) >> 14;
      //fHCHWords = (*fDataWord >> 14) & 0x7;
    }

  // word 1
  //static const UInt_t skTimeBinsBits = 0xfc << 26;
  static const UInt_t skTimeBinsBits = 0x3f << 26;
  static const UInt_t skBunchCrossCounterBits = 0xffff << 10;
  static const UInt_t skPreTriggerCounterBits = 0xf << 6;
  static const UInt_t skPreTriggerPhase = 0xf << 2;
  
  static UInt_t GetTimeBins(UInt_t *word)
    {
      return (*word & skTimeBinsBits) >> 26;
    }

  static UInt_t GetBunchCrossCounter(UInt_t *word)
    {
      return (*word & skBunchCrossCounterBits) >> 10;
    }

  static UInt_t GetPreTriggerCounter(UInt_t *word)
    {
      return (*word & skPreTriggerCounterBits) >> 6;
    }

  static UInt_t GetPreTriggerPhase(UInt_t *word)
    {
      return (*word & skPreTriggerPhase) >> 2;
    }

  struct AliTRDrawHCInfo
  {
    UInt_t fSpecialRawV;
    UInt_t fRawVMajor;
    UInt_t fRawVMinor;
    UInt_t fNExtraWords;
    UInt_t fDCSboard;
    UInt_t fSM;
    UInt_t fStack;
    UInt_t fLayer;
    UInt_t fSide;

/*     fBCctr   =  (*fDataWord >> 16); */
/*     fPTctr   =  (*fDataWord >> 12) & 0xf; */
/*     fPTphase =  (*fDataWord >>  8) & 0xf; */
/*     fTBins   = ((*fDataWord >>  2) & 0x3f) + 1; */
/*     fTimeWords = (fTBins - 1)/3 + 1;	 */

    UInt_t fTimeBins;
    UInt_t fBunchCrossCounter;
    UInt_t fPreTriggerCounter;
    UInt_t fPreTriggerPhase;

    void DecodeH0(UInt_t *word)
    {
      fSpecialRawV = GetSpecialRawV(word);
      fRawVMajor = GetRawVMajor(word);
      fRawVMinor = GetRawVMinor(word);
      fNExtraWords = GetExtraWords(word);
      fDCSboard = GetDCSBoard(word);
      fSM = GetSM(word);
      fStack = GetStack(word);
      fLayer = GetLayer(word);
      fSide = GetSide(word);
    }

    void DecodeH1(UInt_t *word)
    {
      fTimeBins = GetTimeBins(word);
      fBunchCrossCounter = GetBunchCrossCounter(word);
      fPreTriggerCounter = GetPreTriggerCounter(word);
      fPreTriggerPhase = GetPreTriggerPhase(word);
    }

    void Dump()
    {
/*       printf("[I] HC Info is : DCSboard %d SM %d Stack %d Layer %d Side %d \n", */
/* 	     fDCSboard, fSM, fStack, fLayer, fSide);	      */
      printf("[I] HC Info is : RawV %d SM %d Stack %d Layer %d Side %d \n",
	     fRawVMajor, fSM, fStack, fLayer, fSide);	     
    }
  };

  //--------------------------------------------------------
  // MCM word masks
  static const UInt_t skROBBits = 0x70000000; // 0x7 << 28
  static const UInt_t skMCMBits = 0x0f000000; // 0xf << 24
  static const UInt_t skHCEvCounterBits = 0x00fffff0; //0x00fffff << 4

  static UInt_t GetROB(UInt_t *word)
    {
      return (*word & skROBBits) >> 28;
    }

  static UInt_t GetMCM(UInt_t *word)
    {
      return (*word & skMCMBits) >> 24;
    }

  static UInt_t GetEvCounter(UInt_t *word)
    {
      return (*word & skHCEvCounterBits) >> 4;
    }

  struct AliTRDrawMCMInfo
  {
    UInt_t fROB;
    UInt_t fMCM;
    UInt_t fEvCounter;

    void Decode(UInt_t *word)
    {
      fROB = GetROB(word);
      fMCM = GetMCM(word);
      fEvCounter = GetEvCounter(word);
    }

    void Dump()
    {
      printf("[I] MCM Info is : ROB %d MCM %d EvCounter %d\n", fROB, fMCM, fEvCounter);
    }
  };

}

class AliTRDRawStreamTB: public TObject {

  public :

    AliTRDRawStreamTB();
    AliTRDRawStreamTB(AliRawReader *rawReader);
    virtual ~AliTRDRawStreamTB();

    virtual Bool_t       Next();              // Read the next data
    virtual Int_t NextChamber(AliTRDdigitsManager *man); // read next chamber data
    virtual Int_t        Init();              // Init for the fRawVersion > 1

    enum { kDDLOffset = 0x400 };              // Offset for DDL numbers

    Bool_t               SetRawVersion(Int_t rv);
    Int_t                GetRawVersion() const                      { return fRawVersion;     };
    void                 SetRawReader(AliRawReader *rawReader);

    // Get Filter settings (does not make a lot of sense):
    Int_t                TRAPfilterTCon() const                     { return fTCon;           };
    Int_t                TRAPfilterPEDon() const                    { return fPEDon;          };
    Int_t                TRAPfilterGAINon() const                   { return fGAINon;         };
    Int_t                TRAPsendsUnfilteredData() const            { return fBypass;         };

    // Get Tracklet parameters (does not make a lot of sense):
    Float_t              GetTrackletPID() const                     { return fTracklPID;      };
    Float_t              GetTrackletDeflLength() const              { return fTracklDefL;     };
    Float_t              GetTrackletPadPos() const                  { return fTracklPadPos;   };
    Int_t                GetTrackletPadRow() const                  { return fTracklPadRow;   };

    // Check if the link has optical power (HC sends data)
    Bool_t               IsGTULinkActive(Int_t sm, Int_t la, Int_t sta, Int_t side)
      { return ( ((fGTUlinkMask[sm][sta]) >> (2*la+side)) & 0x1 ); };

    static void SetStackIndexBug(Bool_t val) {fgStackIndexBug = val;}
    static void SetForceRawVersion(Int_t val) {fgForceRawVersion = val;}
    static void SupressWarnings(Bool_t val) {fgSupressWarnings = val;}
    static void ExtraDebug(Bool_t val) {fgExtraDebug = val;}
    static void RawBufferMissAligned(Bool_t val) {fgRawDataHack = val;}

    Int_t *GetSignals() { return fSig;}                         //  Signals in the three time bins from Data Word
    Int_t GetADC() const { return fADC;}                            //  MCM ADC channel and Time Bin of word 1
    Int_t GetTimeBin() const { return fTB - 3;}                             //  MCM ADC channel and Time Bin of word 1
    Int_t GetEventNumber() const { return fEv;}                             //  MCM Event number and position of current MCM on TRD chamber
    Int_t GetROB() const { return fROB;}                            //  MCM Event number and position of current MCM on TRD chamber
    Int_t GetMCM() const { return fMCM;}                           //  MCM Event number and position of current MCM on TRD chamber
    Int_t GetSM() const { return fSM;}                             //  Position of CURRENT half chamber in full TRD
    Int_t GetLayer() const { return fLAYER;}                          //  PLANE = Position of CURRENT half chamber in full TRD
    Int_t GetStack() const { return fSTACK;}                          //  CHAMBER = Position of CURRENT half chamber in full TRD
    Int_t GetROC() const { return fROC;}                            //  Position of CURRENT half chamber in full TRD
    Int_t GetSide() const { return fSIDE;}                           //  Position of CURRENT half chamber in full TRD
    Int_t GetDCS() const { return fDCS;}                            //  DCS board number read from data (HC header)
    Int_t GetRow() const { return fROW;}
    Int_t GetCol() const { return fCOL;}                            //  Detector Pad coordinates
    Int_t GetDet() const { return fDET;} // helper
    Int_t GetLastDet() const { return fLastDET;} // helper
    Int_t GetMaxRow() const { return fRowMax;}
    Int_t GetMaxCol() const { return fColMax;}
    Int_t GetNumberOfTimeBins() const {return fTBins;}

    Int_t  GetCommonAdditive() const {return fCommonAdditive;}
    Bool_t IsCurrentPadShared() const {return fIsPadShared;}
    void   SetSharedPadReadout(Bool_t fv) {fSharedPadsOn = fv;}
    Bool_t IsDataZeroSuppressed() const {return fZeroSuppressed;}
      
  private :

    Int_t    fSig[3];                         //  Signals in the three time bins from Data Word
    Int_t    fADC;                            //  MCM ADC channel and Time Bin of word 1
    Int_t    fMaxADCgeom;                     //  Max ADC number as returned by AliTRDgeometry::ADCmax()
    Int_t    fTB;                             //  MCM ADC channel and Time Bin of word 1
    Int_t    fEv;                             //  MCM Event number and position of current MCM on TRD chamber
    Int_t    fROB;                            //  MCM Event number and position of current MCM on TRD chamber
    Int_t    fMCM;                            //  MCM Event number and position of current MCM on TRD chamber
    Int_t    fSM;                             //  Position of CURRENT half chamber in full TRD
    Int_t    fLAYER;                          //  Position of CURRENT half chamber in full TRD
    Int_t    fSTACK;                          //  Position of CURRENT half chamber in full TRD
    Int_t    fROC;                            //  Position of CURRENT half chamber in full TRD
    Int_t    fSIDE;                           //  Position of CURRENT half chamber in full TRD
    Int_t    fDCS;                            //  DCS board number read from data (HC header)
    Int_t    fROW;                            //  Detector Row coordinates
    Int_t    fCOL;                            //  Detector Pad coordinates
    Int_t    fDET;                            //  Current detector - version > 1
    Int_t    fLastDET;                        //  Previous detector - version > 1

    Int_t    fBCctr;                          //  Counters from HC header (>=V2)
    Int_t    fPTctr;                          //  Counters from HC header (>=V2)
    Int_t    fPTphase;                        //  Counters from HC header (>=V2)
    Int_t    fRVmajor;                        //  Raw version numbers and number of additional HC headerwords (>=V2)
    Int_t    fRVminor;                        //  Raw version numbers and number of additional HC headerwords (>=V2)
    Int_t    fHCHWords;                       //  Raw version numbers and number of additional HC headerwords (>=V2)
    Int_t    fTBins;                          //  Number of time bins read from HC header (>=V2)
    Bool_t   fTCon;                           //  Filter settings read from HC header (>=V2)
    Bool_t   fPEDon;                          //  Filter settings read from HC header (>=V2)
    Bool_t   fGAINon;                         //  Filter settings read from HC header (>=V2)
    Bool_t   fXTon;                           //  Filter settings read from HC header (>=V2)
    Bool_t   fNonLinOn;                       //  Filter settings read from HC header (>=V2)
    Bool_t   fBypass;                         //  Filter settings read from HC header (>=V2)
    Int_t    fCommonAdditive;                 //  Common baseline additive read from HC header (>=V2)

    Bool_t   fZeroSuppressed;                 // Data is zero suppressed

    Int_t    fHCHctr1;                        //  HC and MCM Header counter
    Int_t    fHCHctr2;                        //  HC and MCM Header counter
    Int_t    fMCMHctr1;                       //  HC and MCM Header counter
    Int_t    fMCMHctr2;                       //  HC and MCM Header counter
    Int_t    fGTUctr1;                        //  GTU LinkMask Counter
    Int_t    fGTUctr2;                        //  GTU LinkMask Counter
    Int_t    fHCdataCtr;                      //  Data Counter for half chamber
    
    Float_t  fTracklPID;                      //  Tracklet parameters
    Float_t  fTracklDefL;                     //  Tracklet parameters
    Float_t  fTracklPadPos;                   //  Tracklet parameters
    Int_t    fTracklPadRow;                   //  Tracklet parameters

    UShort_t fGTUlinkMask[18][5];             //  Mask with active links

    Int_t    fMCMWordCrt;                     //  Word Counter for a single MCM
    Int_t    fMCMWordsExpected;                //  Expected Words from MCM 

    AliTRDRawStreamTB(const AliTRDRawStreamTB &stream);
    AliTRDRawStreamTB &operator=(const AliTRDRawStreamTB &stream);
    void SwapOnEndian();

    AliRawReader *fRawReader;              //  Object for reading the raw data

    Int_t    fRawVersion;                  //  Which version of raw data decoding is used
    Int_t    fRawDigitThreshold;           //  Naive "zero"(threshold) supression. Usefull for Raw Data ver 2.

    Int_t    fNextStatus;                  //  Navigation in raw versions > 1
    Int_t    fLastStatus;                  //  Navigation in raw versions > 1
    UInt_t   fTbSwitch;                    //  Time Bin Switch for internal use
    UInt_t   fTbSwitchCtr;                 //  Counter for internal use
    UInt_t   fTimeWords;                   //  Number of Words needed to store the data of 1 ADC ch.
    UInt_t   fWordCtr;                     //  Word Counter

    Int_t    fRowMax;                      //  Maximum number of pad rows and columns
    Int_t    fColMax;                      //  Maximum number of pad rows and columns

    Bool_t   fADCmask[21];                 //  Mask of active ADCs for zero suppression
    UInt_t   fLastADCmask;                 //  Last ADC read out

    Int_t    fChamberDone[540];            //  Chamber was processed already (1=1HC, 2=full chamber)

    Int_t    fRetVal;                      //  Datadecode return
    Int_t    fEqID;                        //  Equipment id
    UInt_t   fDataSize;                    //  Size of the data available in the current buffer
    Bool_t   fSizeOK;                      //  Did we read anything
    UInt_t   fCountBytes;                  //  Bytes traversed in the buffer
    UInt_t   fBufSize;                     //  Size of the current RawReader buffer
    Bool_t   fkBufferSet;                  //  Is the RawReader buffer available
    UChar_t *fPos;                         //  Position in the buffer of the RawReader
    UInt_t  *fDataWord;                    //  The pointer to the current 32 bit data word
    UInt_t   fTimeBinsCalib;               //  N of time bins retrieved from calibration

    Int_t    fADClookup[32];               //  Lookup for version 3 (1[entries]+1[index]+30[fADC channels] = 32)
    Int_t    fNActiveADCs;                 //  Number of active ADC channels
    Bool_t   fEndOfDataFlag;               //  End of data flag

    AliTRDrawDataUtilsTB::AliTRDrawSMinfo    fSMinfo[18];            // Super module info
    AliTRDrawDataUtilsTB::AliTRDrawStackInfo fStackInfo[18][5];      // Stack info
    AliTRDrawDataUtilsTB::AliTRDrawHCInfo    fHCInfo;                // Half chamber info
    AliTRDrawDataUtilsTB::AliTRDrawMCMInfo   fMCMInfo;               // MCM info

    Int_t    fiSMx;                        //  Position of CURRENT half chamber in full TRD
    static Bool_t   fgStackIndexBug;        //  Buggy index word flag

    Bool_t   fSharedPadsOn;                // Should we return on shared pad readout
    Bool_t   fIsPadShared;                 // Set if the current pad is shared
    static Int_t   fgForceRawVersion;             // Force the raw version specified
    static Bool_t  fgSupressWarnings;       // Superss warnings (raw data version and TB missmatch)
    static Bool_t  fgExtraDebug;           // Extra info on the screen
    static Bool_t  fgRawDataHack;          // skip 23 words - a hack for the DATE raw format
    enum ETRDzRawStreamError {
       kHCWordMissing = 1                  //
      ,kMCMADCMaskMissing = 2              //
      ,kWrongMCMorROB = 3                  //
      ,kGTULinkMaskMissing = 4             //
      ,kHCHeaderCorrupt = 5                //
      ,kHCHeaderMissing = 6                //
      ,kROBSideMismatch = 7                //
      ,kWrongPadrow = 8                    //
      ,kWrongPadcolumn = 9                 //
      ,kTrackletRowMismatch = 10           //
      ,kDataMaskError = 11                 //
      ,kADCNumberOverflow = 12             //
      ,kADCChannelOverflow = 13            //
    };

 protected:

    AliTRDgeometry *fGeo;                  //  TRD geometry

    Bool_t DecodeGTUlinkMask();
    Bool_t DecodeADCWord();
    Bool_t DecodeNextRawWord();
    Bool_t DecodeMCM();
    Bool_t DecodeHC();
    Bool_t DecodeSM();
    inline void  ChangeStatus(Int_t kstat);

    void  DecodeHCheader(Int_t timeBins = 0);
    void  DecodeMCMheader();
    void  DecodeTracklet();

    void  SetRawDigitThreshold (Int_t ith) {fRawDigitThreshold = ith;} // set the naive zero suppression threshold

    Int_t NextData(); // get the next piece of memory
    Int_t RewindWord(); //go back one word

    Int_t ChannelsToRead(Int_t ADCmask); // Get the active ADC channels from the mask (V3 and 2)

    Int_t SkipWords(UInt_t iw);
    Int_t DecodeHeadingInfo();

    enum { kStart
         , kStop
         , kWordOK
         , kNoMoreData
         , kNextSM
         , kNextHC
         , kSeekNonEoTracklet
         , kDecodeHC
         , kNextMCM
         , kNextData
         , kReading };
    
    ClassDef(AliTRDRawStreamTB, 1)               // Class for reading TRD raw digits

};
#endif
