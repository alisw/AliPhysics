#ifndef ALITRDRAWSTREAMV2_H
#define ALITRDRAWSTREAMV2_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// This class provides access to TRD digits in raw data.                     //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include <TObject.h>
#include "AliTRDrawStreamBase.h"

class AliTRDgeometry;
class AliRawReader;
class AliTRDdigitsManager;

//class AliTRDRawStreamV2: public TObject {
class AliTRDRawStreamV2 : public AliTRDrawStreamBase {

  public :

    AliTRDRawStreamV2();
    AliTRDRawStreamV2(AliRawReader *rawReader);
    virtual ~AliTRDRawStreamV2();

    virtual Bool_t       Next();                                //  Read the next data
    virtual Int_t        NextChamber(AliTRDdigitsManager *man, UInt_t **trackletContainer); //  Read next chamber data
    virtual Bool_t        Init();                                //  Init for the fRawVersion > 1

    enum { kDDLOffset = 0x400 };                                //  Offset for DDL numbers

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

    Int_t *GetSignals() const    { return (Int_t*)fSig;     } //  Signals in the three time bins from Data Word
    Int_t GetADC() const         { return fADC;     } //  MCM ADC channel and Time Bin of word 1
    Int_t GetTimeBin() const     { return fTB - 3;  } //  MCM ADC channel and Time Bin of word 1
    Int_t GetEventNumber() const { return fEv;      } //  MCM Event number and position of current MCM on TRD chamber
    Int_t GetROB() const         { return fROB;     } //  MCM Event number and position of current MCM on TRD chamber
    Int_t GetMCM() const         { return fMCM;     } //  MCM Event number and position of current MCM on TRD chamber
    Int_t GetSM() const          { return fSM;      } //  Position of CURRENT half chamber in full TRD
    Int_t GetLayer() const       { return fLAYER;   } //  PLANE = Position of CURRENT half chamber in full TRD
    Int_t GetStack() const       { return fSTACK;   } //  CHAMBER = Position of CURRENT half chamber in full TRD
    Int_t GetROC() const         { return fROC;     } //  Position of CURRENT half chamber in full TRD
    Int_t GetSide() const        { return fSIDE;    } //  Position of CURRENT half chamber in full TRD
    Int_t GetDCS() const         { return fDCS;     } //  DCS board number read from data (HC header)
    Int_t GetRow() const         { return fROW;     }
    Int_t GetCol() const         { return fCOL;     } //  Detector Pad coordinates
    Int_t GetDet() const         { return fDET;     } //  Helper
    Int_t GetLastDet() const     { return fLastDET; } //  Helper
    Int_t GetMaxRow() const      { return fRowMax;  }
    Int_t GetMaxCol() const      { return fColMax;  }
    Int_t GetNumberOfTimeBins() const 
                                 { return fTBins;   }
    
    void SwapOnEndian();
	static void    SetDumpHead(UInt_t iv) {fgDumpHead = iv;} // set number of words to be dumped on the screen 

 protected:

    AliTRDgeometry *fGeo;                     //  TRD geometry

    Bool_t DecodeGTUlinkMask();
    Bool_t DecodeNextRawWord();
    Bool_t DecodeMCM();
    Bool_t DecodeHC();
    Bool_t DecodeSM();
	Bool_t DumpWords(UInt_t *px, UInt_t iw, UInt_t marker = 0); // dump some words onto the screen - debugging purpose
    inline void ChangeStatus(Int_t kstat);
/*       { */
/* 	fLastStatus = fNextStatus; */
/* 	fNextStatus = kstat;   */
/*       } */

    void  DecodeHCheader(Int_t timeBins = 0);
    void  DecodeMCMheader();
    void  DecodeTracklet();

    void  SetRawDigitThreshold (Int_t ith) 
                { fRawDigitThreshold = ith; } //  Set the naive zero suppression threshold

    Int_t NextData();                         //  Get the next piece of memory
    Int_t ChannelsToRead(Int_t ADCmask);      //  Get the active ADC channels from the mask (V3 and 2)

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

    // Some constants:
    static const UInt_t fgkEndoftrackletmarker = 0xAAAAAAAA;     //  This marks the end of tracklet data words
    static const UInt_t fgkEndofrawdatamarker  = 0x00000000;     //  This marks the end of half-chamber-data
    static const UInt_t fgkSizeWord            = sizeof(UInt_t); //  Size of a word in bytes

	static UInt_t fgDumpHead;                                    // number of words to be dumped (from the start of the buffer)

  private :

    Int_t    fSig[3];                                            //  Signals in the three time bins from Data Word
    Int_t    fADC;                                               //  MCM ADC channel and Time Bin of word 1
    Int_t    fTB;                                                //  MCM ADC channel and Time Bin of word 1
    Int_t    fEv;                                                //  MCM Event number and position of current MCM on TRD chamber
    Int_t    fROB;                                               //  MCM Event number and position of current MCM on TRD chamber
    Int_t    fMCM;                                               //  MCM Event number and position of current MCM on TRD chamber
    Int_t    fSM;                                                //  Position of CURRENT half chamber in full TRD
    Int_t    fLAYER;                                             //  Position of CURRENT half chamber in full TRD
    Int_t    fSTACK;                                             //  Position of CURRENT half chamber in full TRD
    Int_t    fROC;                                               //  Position of CURRENT half chamber in full TRD
    Int_t    fSIDE;                                              //  Position of CURRENT half chamber in full TRD
    Int_t    fDCS;                                               //  DCS board number read from data (HC header)
    Int_t    fROW;                                               //  Detector Row coordinates
    Int_t    fCOL;                                               //  Detector Pad coordinates
    Int_t    fDET;                                               //  Current detector - version > 1
    Int_t    fLastDET;                                           //  Previous detector - version > 1

    Int_t    fBCctr;                                             //  Counters from HC header (>=V2)
    Int_t    fPTctr;                                             //  Counters from HC header (>=V2)
    Int_t    fPTphase;                                           //  Counters from HC header (>=V2)
    Int_t    fRVmajor;                                           //  Raw version numbers and number of additional HC headerwords (>=V2)
    Int_t    fRVminor;                                           //  Raw version numbers and number of additional HC headerwords (>=V2)
    Int_t    fHCHWords;                                          //  Raw version numbers and number of additional HC headerwords (>=V2)
    Int_t    fTBins;                                             //  Number of time bins read from HC header (>=V2)
    Bool_t   fTCon;                                              //  Filter settings read from HC header (>=V2)
    Bool_t   fPEDon;                                             //  Filter settings read from HC header (>=V2)
    Bool_t   fGAINon;                                            //  Filter settings read from HC header (>=V2)
    Bool_t   fXTon;                                              //  Filter settings read from HC header (>=V2)
    Bool_t   fNonLinOn;                                          //  Filter settings read from HC header (>=V2)
    Bool_t   fBypass;                                            //  Filter settings read from HC header (>=V2)
    Int_t    fCommonAdditive;                                    //  Common baseline additive read from HC header (>=V2)

    Bool_t   fZeroSuppressed;                                    //  Data is zero suppressed

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
    Int_t    fMCMWordsExpected;               //  Expected Words from MCM 

    AliTRDRawStreamV2(const AliTRDRawStreamV2 &stream);
    AliTRDRawStreamV2 &operator=(const AliTRDRawStreamV2 &stream);

    AliRawReader *fRawReader;                 //  Object for reading the raw data

    Int_t    fRawVersion;                     //  Which version of raw data decoding is used
    Int_t    fRawDigitThreshold;              //  Naive "zero"(threshold) supression. Usefull for Raw Data ver 2.

    Int_t    fNextStatus;                     //  Navigation in raw versions > 1
    Int_t    fLastStatus;                     //  Navigation in raw versions > 1
    UInt_t   fTbSwitch;                       //  Time Bin Switch for internal use
    UInt_t   fTbSwitchCtr;                    //  Counter for internal use
    UInt_t   fTimeWords;                      //  Number of Words needed to store the data of 1 ADC ch.
    UInt_t   fWordCtr;                        //  Word Counter

    Int_t    fRowMax;                         //  Maximum number of pad rows and columns
    Int_t    fColMax;                         //  Maximum number of pad rows and columns

    Bool_t   fADCmask[21];                    //  Mask of active ADCs for zero suppression
    UInt_t   fLastADCmask;                    //  Last ADC read out

    UShort_t fChamberDone[540];               //  Chamber was processed already (1=1HC, 2=full chamber)

    Int_t    fRetVal;                         //  Datadecode return
    Int_t    fEqID;                           //  Equipment id
    UInt_t   fDataSize;                       //  Size of the data available in the current buffer
    Bool_t   fSizeOK;                         //  Did we read anything
    UInt_t   fCountBytes;                     //  Bytes traversed in the buffer
    UInt_t   fBufSize;                        //  Size of the current RawReader buffer
    Bool_t   fkBufferSet;                     //  Is the RawReader buffer available
    UChar_t *fPos;                            //  Position in the buffer of the RawReader
	UInt_t  *fpBegin; 						  // begin - pointer to the buffer word 0
	UInt_t  *fpEnd;   					      // end of the buffer
    UInt_t  *fDataWord;                       //  The pointer to the current 32 bit data word
    UInt_t   fTimeBinsCalib;                  //  N of time bins retrieved from calibration

    Int_t    fADClookup[32];                  //  Lookup for version 3 (1[entries]+1[index]+30[fADC channels] = 32)
    Int_t    fNActiveADCs;                    //  Number of active ADC channels
    Bool_t   fEndOfDataFlag;                  //  End of data flag

    enum ETRDzRawStreamError {
       kHCWordMissing       =  1              //
      ,kMCMADCMaskMissing   =  2              //
      ,kWrongMCMorROB       =  3              //
      ,kGTULinkMaskMissing  =  4              //
      ,kHCHeaderCorrupt     =  5              //
      ,kHCHeaderMissing     =  6              //
      ,kROBSideMismatch     =  7              //
      ,kWrongPadrow         =  8              //
      ,kWrongPadcolumn      =  9              //
      ,kTrackletRowMismatch = 10              //
      ,kDataMaskError       = 11              //
      ,kADCNumberOverflow   = 12              //
      ,kADCChannelOverflow  = 13              //
    };
    
    ClassDef(AliTRDRawStreamV2,1)                                //  Class for reading TRD raw digits

};
#endif
