#ifndef ALITRDRAWSTREAM_H
#define ALITRDRAWSTREAM_H
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
class AliTRDdataArrayI;

// Some constants:
const UInt_t kEndoftrackletmarker = 0xAAAAAAAA; /*This marks the end of tracklet data words*/
const UInt_t kEndofrawdatamarker  = 0x00000000; /*This marks the end of half-chamber-data*/

class AliTRDRawStream: public TObject {

  public :

    AliTRDRawStream();
    AliTRDRawStream(AliRawReader *rawReader);
    AliTRDRawStream(AliRawReader *rawReader, AliTRDdigitsManager *man, AliTRDdataArrayI *dig);
    virtual ~AliTRDRawStream();

    virtual Bool_t       Next();              // Next function (for fRawVersion = 0 (Bogdans first version))
    virtual Int_t        ReadAll();           // Read function (for fRawVersion > 0)

    Int_t                GetDetectorV0() const     /* only for v0 */  { return fDetector;       };
    Int_t                GetPrevDetectorV0() const /* only for v0 */  { return fPrevDetector;   };
    Bool_t               IsNewDetectorV0() const   /* only for v0 */  { return fDetector != fPrevDetector; };
    Int_t                GetNPadsV0() const        /* only for v0 */  { return fNPads;          };
    Int_t                GetRowV0() const          /* only for v0 */  { return fRow;            };
    Int_t                GetPrevRowV0() const      /* only for v0 */  { return fPrevRow;        };
    Bool_t               IsNewRowV0() const        /* only for v0 */  { return (fRow != fPrevRow) || IsNewDetectorV0();  };
    Int_t                GetColumnV0() const       /* only for v0 */  { return fColumn;         };
    Int_t                GetPrevColumnV0() const   /* only for v0 */  { return fPrevColumn;     };
    Bool_t               IsNewColumnV0() const     /* only for v0 */  { return (fColumn != fPrevColumn) || IsNewRowV0(); };
    Int_t                GetTimeV0() const         /* only for v0 */  { return fTime-1;         };
    Int_t                GetSignalV0() const       /* only for v0 */  { return fSignal;         };

    enum { kDDLOffset = 0x400 };              // Offset for DDL numbers

    void                 SetDigitsManager(AliTRDdigitsManager *man) { fDigitsManager   = man; };
    void                 SetDigits(AliTRDdataArrayI *dig)           { fDigits          = dig; };

    AliTRDdigitsManager *GetDigitsManager() const                   { return fDigitsManager;  };

    Bool_t               SetRawVersion(Int_t rv);
    Int_t                GetRawVersion() const                      { return fRawVersion;     };

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


  private :

    Int_t    fSig[3];                         //  Signals in the three time bins from Data Word
    Int_t    fADC;                            //  MCM ADC channel and Time Bin of word 1
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
    Int_t    fROW;
    Int_t    fCOL;                            //  Detector Pad coordinates

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

    AliTRDRawStream(const AliTRDRawStream &stream);
    AliTRDRawStream &operator=(const AliTRDRawStream &stream);

    AliRawReader *fRawReader;              //  Object for reading the raw data

    // The following is used for v0:
    Int_t    fCount;                       //  Counter of bytes to be read for current detector
    Int_t    fDetector;                    //  Index of current detector
    Int_t    fPrevDetector;                //  Index of previous detector
    Int_t    fNPads;                       //  Number of active pads
    Int_t    fRow;                         //  Index of current pad row
    Int_t    fPrevRow;                     //  Index of previous pad row
    Int_t    fColumn;                      //  Index of current pad column
    Int_t    fPrevColumn;                  //  Index of previous pad column
    Int_t    fTime;                        //  Index of current time bin
    Int_t    fSignal;                      //  Signal in ADC counts

    // This is again new:
    Int_t    fRawVersion;                  //  Which version of raw data decoding is used
    UInt_t   fDataWord;                    //  The current 32 bit data word
    Int_t    fStatus;                      //  Status word used by some functions
    UInt_t   fTbSwitch;                    //  Time Bin Switch for internal use
    UInt_t   fTbSwitchCtr;                 //  Counter for internal use
    UInt_t   fTimeWords;                   //  Number of Words needed to store the data of 1 ADC ch.
    UInt_t   fWordCtr;                     //  Word Counter

    Int_t    fRowMax;                      //  Maximum number of pad rows and columns
    Int_t    fColMax;                      //  Maximum number of pad rows and columns

    Bool_t   fADCmask[21];                 //  Mask of active ADCs for zero suppression
    UShort_t fChamberDone[540];            //  Chamber was processed already (1=1HC, 2=full chamber)

 protected:

    AliTRDgeometry *fGeo;                  //  TRD geometry

    AliTRDdigitsManager *fDigitsManager;   // Manager for the output digits
    AliTRDdataArrayI    *fDigits;          // Output digits
    AliTRDdataArrayI    *fTrack0;          // Track dictionary
    AliTRDdataArrayI    *fTrack1;          // Track dictionary
    AliTRDdataArrayI    *fTrack2;          // Track dictionary

    void  DecodeHCheader(Int_t timeBins);
    void  DecodeHCheaderV1();                  // Valid for fRawversion = 1
    void  DecodeHCheaderV2V3(Int_t timeBins);  // Valid for fRawversion = 2,3,?

    void  DecodeMCMheader();
    void  DecodeMCMheaderVx();                 // So far valid for all fRawversion = 1,2,3, ...

    void  DecodeTracklet();
    void  DecodeTrackletVx();                  // So far valid for all fRawversion = 1,2,3, ...

    void  DecodeGTUlinkMask();
    void  DecodeGTUlinkMaskVx();               // So far valid for all fRawversion = 1,2,3, ...

    Int_t DecodeDataWord();
    Int_t DecodeDataWordV1V2();                // Valid for fRawversion = 1, 2, ...
    Int_t DecodeDataWordV3();                  // Valid for fRawversion = 3, ...

    ClassDef(AliTRDRawStream, 3)               // Class for reading TRD raw digits

};
#endif
