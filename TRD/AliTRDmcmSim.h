#ifndef ALITRDMCMSIM_H
#define ALITRDMCMSIM_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

///////////////////////////////////////////////////////
//                                                   //
//  Multi Chip Module Simulation Class               //
//                                                   //
///////////////////////////////////////////////////////

#include <TObject.h>
#include "AliTRDCommonParam.h"
#include "AliTRDcalibDB.h"

class TClonesArray;

class AliRunLoader;
class AliTRDfeeParam;
class AliTRDSimParam;
class AliTRDtrapConfig;
class AliTRDcalibDB;
class AliTRDgeometry;
class AliTRDpadPlane;
class AliTRDarrayADC;
class AliTRDdigitsManager;

class AliTRDmcmSim : public TObject {
 public:
                    AliTRDmcmSim();
  virtual          ~AliTRDmcmSim();

          void      Init( Int_t cha, Int_t rob, Int_t mcm, Bool_t newEvent = kFALSE );   // Initialize MCM by the position parameters
	  void      Reset();                                   // clears filter registers and internal data

	  Bool_t    LoadMCM(AliRunLoader* const runloader, Int_t det, Int_t rob, Int_t mcm);
	  void      NoiseTest(Int_t nsamples, Int_t mean, Int_t sigma, Int_t inputGain = 1, Int_t inputTail = 2);

          void      SetData(Int_t iadc, Int_t *adc);           // Set ADC data with array 
          void      SetData(Int_t iadc, Int_t it, Int_t adc ); // Set ADC data
	  void      SetData(AliTRDarrayADC *adcArray, 
			    AliTRDdigitsManager *digitsManager = 0x0);         // Set ADC data from adcArray
          void      SetDataPedestal(Int_t iadc );              // Fill ADC data with pedestal values
  static  void      SetApplyCut(Bool_t applyCut) { fgApplyCut = applyCut; }

          Int_t     GetDetector() const  { return fDetector;  };     // Returns Chamber ID (0-539)
          Int_t     GetRobPos() const { return fRobPos; };     // Returns ROB position (0-7)
          Int_t     GetMcmPos() const { return fMcmPos; };     // Returns MCM position (0-17) (16,17 are mergers)
          Int_t     GetRow() const    { return fRow;    };     // Returns Row number on chamber where the MCM is sitting
	  Int_t     GetCol( Int_t iadc );                      // Get corresponding column (0-143) from for ADC channel iadc = [0:20]
	  // for the ADC/Col mapping, see: http://wiki.kip.uni-heidelberg.de/ti/TRD/index.php/Image:ROB_MCM_numbering.pdf
  static  Bool_t    GetApplyCut() { return fgApplyCut; }

	  void WriteData(AliTRDarrayADC *digits);

	  Int_t     ProduceRawStream( UInt_t *buf, Int_t bufsize, UInt_t iEv = 0 ); // Produce raw data stream - Read data format
	  Int_t     ProduceTrackletStream( UInt_t *buf, Int_t bufsize ); // produce the tracklet stream for this MCM
		
	  // different stages of processing in the TRAP
	  void      Filter();                                  // Apply digital filters for existing data (according to configuration)
	  void      ZSMapping();                               // Do ZS mapping for existing data
	  void      Tracklet();                                // Run tracklet preprocessor and perform tracklet fit

	  // apply individual filters to all channels and timebins
	  void      FilterPedestal();                   // Apply pedestal filter
	  void      FilterGain();                       // Apply gain filter
	  void      FilterTail();                       // Apply tail filter

	  // filter initialization (resets internal registers)
	  void      FilterPedestalInit();
	  void      FilterGainInit();
	  void      FilterTailInit(Int_t baseline); //??? automatic baseline??

	  // feed single sample to individual filter
	  // this changes the internal registers
	  // all filters operate on 12-bit values!
	  UShort_t  FilterPedestalNextSample(Int_t adc, Int_t timebin, UShort_t value);
	  UShort_t  FilterGainNextSample(Int_t adc, UShort_t value);
	  UShort_t  FilterTailNextSample(Int_t adc, UShort_t value);

	  // tracklet calculation
	  void      AddHitToFitreg(Int_t adc, UShort_t timebin, UShort_t qtot, Short_t ypos, Int_t label);
	  void      CalcFitreg();
	  void      TrackletSelection();
	  void      FitTracklet();

	  TClonesArray* GetTrackletArray() const { return fTrackletArray; }

	  // data display
	  void      Print(Option_t* const option="") const;   // print stored data to stdout
	  void      Draw(Option_t* const option ="");         // draw data (ADC data, hits and tracklets)
	  void      DumpData( char *f, char *target );  // Dump data stored (only for debugging)

 protected:
	  Bool_t    CheckInitialized();                 // Check whether the class is initialized
	  
	  Bool_t    fInitialized;                       // Status whether the class is initialized or not
	  Int_t     fMaxTracklets;                      // maximum number of tracklet-words submitted per mcm **new** //???
	  Int_t     fDetector;                             // Chamber ID
	  Int_t     fRobPos;                            // ROB Position on chamber 
	  Int_t     fMcmPos;                            // MCM Position on chamber 
	  Int_t     fRow;                               // Pad row number (0-11 or 0-15) of the MCM on chamber
	  Int_t     fNADC;                              // Number of ADC (usually 21) //??? static const
	  Int_t     fNTimeBin;                          // Number of Timebins (variable)  //??? why stored here? taken from TRAPconfig
	  Int_t   **fADCR;                              // Array with MCM ADC values (Raw, 12 bit)
	  Int_t   **fADCF;                              // Array with MCM ADC values (Filtered, 12 bit)
	  UInt_t   *fMCMT;                              // tracklet word for one mcm/trap-chip **new** //??? needed?
	  TClonesArray *fTrackletArray;                 // Array of AliTRDtrackletMCM which contains MC information in addition to the tracklet word
	  Int_t   **fZSM;                               // Zero suppression map
	  Int_t    *fZSM1Dim;                           // Zero suppression map (1 dimensional projection)

	  static const Int_t fgkAddDigits = 2;          // additional digits used for internal representation
	                                                // all internal data as after data control block (i.e. 12 bit), s. TRAP manual
	  static const Int_t fgkNCPU = 4;               // Number of CPUs in the TRAP
	  Int_t     fFitPtr[fgkNCPU];                   // pointer to the tracklet to be calculated by CPU i
	  static const Int_t fgkNHitsMC = 100;          // maximum number of hits for which MC information is kept

	  // Parameter classes
	  AliTRDfeeParam    *fFeeParam;                 // FEE parameters
	  AliTRDtrapConfig  *fTrapConfig;               // TRAP config
	  AliTRDSimParam    *fSimParam;                 // Simulation parameters
	  AliTRDCommonParam *fCommonParam;              // common parameters
	  AliTRDcalibDB     *fCal;                      // Calibration interface
	  AliTRDgeometry    *fGeo;                      // Geometry

	  AliTRDdigitsManager *fDigitsManager;          // pointer to digits manager used for MC label calculation

	  // internal filter registers
	  UInt_t*   fPedAcc;                            // Accumulator for pedestal filter
	  UInt_t*   fGainCounterA;                      // Counter for values above FGTA in the gain filter
	  UInt_t*   fGainCounterB;                      // Counter for values above FGTB in the gain filter
 	  UShort_t* fTailAmplLong;                      // Amplitude of the long component in the tail filter
 	  UShort_t* fTailAmplShort;                     // Amplitude of the short component in the tail filter

	  // hit detection
	  // individual hits can be stored as MC info
	  struct Hit_t {                                // Array of detected hits (only available in MC)
	    Int_t fChannel;                             // ADC channel of the hit
	    Int_t fTimebin;                             // timebin of the hit
	    Int_t fQtot;                                // total charge of the hit
	    Int_t fYpos;                                // calculated y-position
	    Int_t fLabel;                               // label (only in MC)
	  } fHits[fgkNHitsMC];
	  Int_t fNHits;                                 // Number of detected hits

	  // tracklet calculation
	  struct FitReg_t {                             // pointer to the 18 fit registers 
	    Int_t   fNhits;                             // number of hits
	    UInt_t fQ0;                                 // charge accumulated in first window
	    UInt_t fQ1;                                 // charge accumulated in second window
	    UInt_t fSumX;                               // sum x
	    Int_t fSumY;                                // sum y 
	    UInt_t fSumX2;                              // sum x**2
	    UInt_t fSumY2;                              // sum y**2
	    Int_t fSumXY;                               // sum x*y
	  } *fFitReg;

	  //??? cleaning up
	  void Sort2(UShort_t  idx1i, UShort_t  idx2i, UShort_t  val1i, UShort_t  val2i, 
		     UShort_t *idx1o, UShort_t *idx2o, UShort_t *val1o, UShort_t *val2o) const;
	  void Sort3(UShort_t  idx1i, UShort_t  idx2i, UShort_t  idx3i, 
		     UShort_t  val1i, UShort_t  val2i, UShort_t  val3i, 
		     UShort_t *idx1o, UShort_t *idx2o, UShort_t *idx3o, 
		     UShort_t *val1o, UShort_t *val2o, UShort_t *val3o);
	  void Sort6To4(UShort_t  idx1i, UShort_t  idx2i, UShort_t  idx3i, UShort_t  idx4i, UShort_t  idx5i, UShort_t  idx6i, 
			UShort_t  val1i, UShort_t  val2i, UShort_t  val3i, UShort_t  val4i, UShort_t  val5i, UShort_t  val6i, 
			UShort_t *idx1o, UShort_t *idx2o, UShort_t *idx3o, UShort_t *idx4o, 
			UShort_t *val1o, UShort_t *val2o, UShort_t *val3o, UShort_t *val4o);
	  void Sort6To2Worst(UShort_t  idx1i, UShort_t  idx2i, UShort_t  idx3i, UShort_t  idx4i, UShort_t  idx5i, UShort_t  idx6i, 
			     UShort_t  val1i, UShort_t  val2i, UShort_t  val3i, UShort_t  val4i, UShort_t  val5i, UShort_t  val6i, 
			     UShort_t *idx5o, UShort_t *idx6o);

	  UInt_t AddUintClipping(UInt_t a, UInt_t b, UInt_t nbits) const;  // Add a and b (unsigned) with clipping to the maximum value representable by nbits

 private:
	  AliTRDmcmSim(const AliTRDmcmSim &m);             // not implemented
	  AliTRDmcmSim &operator=(const AliTRDmcmSim &m);  // not implemented

	  static Bool_t fgApplyCut;

  ClassDef(AliTRDmcmSim,4)
};

#endif
