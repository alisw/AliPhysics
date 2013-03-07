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

#include <iostream>

class TObject;
class TClonesArray;
class TH2F;

class AliRunLoader;
class AliTRDfeeParam;
class AliTRDtrapConfig;
class AliTRDarrayADC;
class AliTRDarrayDictionary;
class AliTRDdigitsManager;

class AliTRDmcmSim : public TObject {
 public:
                    AliTRDmcmSim();
  virtual          ~AliTRDmcmSim();

          void      Init(Int_t det, Int_t rob, Int_t mcm, Bool_t newEvent = kFALSE);
	  // Initialize MCM by the position parameters

	  void      Reset();
	  // clears filter registers and internal data

	  Bool_t    LoadMCM(AliRunLoader* const runloader, Int_t det, Int_t rob, Int_t mcm);

	  void      NoiseTest(Int_t nsamples, Int_t mean, Int_t sigma, Int_t inputGain = 1, Int_t inputTail = 2);

	  Int_t     GetDataRaw(Int_t iadc, Int_t timebin)      const { return (fADCR[iadc][timebin] >> 2); }
	  // Get unfiltered ADC data
	  Int_t     GetDataFiltered(Int_t iadc, Int_t timebin) const { return (fADCF[iadc][timebin] >> 2); }
	  // Get filtered ADC data

          void      SetData(Int_t iadc, const Int_t* const adc);           // Set ADC data with array
          void      SetData(Int_t iadc, Int_t it, Int_t adc); // Set ADC data
	  void      SetData(AliTRDarrayADC * const adcArray,
			    AliTRDdigitsManager * const digitsManager = 0x0);         // Set ADC data from adcArray
	  void      SetDataByPad(const AliTRDarrayADC *const adcArray,
				 AliTRDdigitsManager * const digitsManager = 0x0);    // Set ADC data from adcArray
          void      SetDataPedestal(Int_t iadc);              // Fill ADC data with pedestal values

  static  Bool_t    GetApplyCut() { return fgApplyCut; }
  static  void      SetApplyCut(Bool_t applyCut) { fgApplyCut = applyCut; }

  static  Int_t     GetAddBaseline() { return fgAddBaseline; }
  static  void      SetAddBaseline(Int_t baseline) { fgAddBaseline = baseline; }
  // Additional baseline which is added for the processing
  // in the TRAP and removed when writing back the data.
  // This is needed to run with TRAP parameters set for a
  // different baseline but it will not change the baseline
  // of the output.

  static  void      SetStoreClusters(Bool_t storeClusters) { fgStoreClusters = storeClusters; }
  static  Bool_t    GetStoreClusters() { return fgStoreClusters; }

          Int_t     GetDetector() const  { return fDetector;  };     // Returns Chamber ID (0-539)
          Int_t     GetRobPos() const { return fRobPos; };           // Returns ROB position (0-7)
          Int_t     GetMcmPos() const { return fMcmPos; };           // Returns MCM position (0-17) (16,17 are mergers)
          Int_t     GetRow() const    { return fRow;    };           // Returns Row number on chamber where the MCM is sitting
	  Int_t     GetCol( Int_t iadc );                            // Get corresponding column (0-143) from for ADC channel iadc = [0:20]
	  // for the ADC/Col mapping, see: http://wiki.kip.uni-heidelberg.de/ti/TRD/index.php/Image:ROB_MCM_numbering.pdf

	  void      WriteData(AliTRDarrayADC *digits);
	  Bool_t    StoreTracklets();                          // Stores tracklets via runloader
	  TString   GetTrklBranchName() const { return fTrklBranchName; }
	  void      SetTrklBranchName(TString name) { fTrklBranchName = name; }

	  Int_t     ProduceRawStream( UInt_t *buf, Int_t bufsize, UInt_t iEv = 0 ) const; // Produce raw data stream - Real data format
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
	  void      FilterPedestalInit(Int_t baseline = 10);
	  void      FilterGainInit();
	  void      FilterTailInit(Int_t baseline = -1);

	  // feed single sample to individual filter
	  // this changes the internal registers
	  // all filters operate on (10+2)-bit values!
	  UShort_t  FilterPedestalNextSample(Int_t adc, Int_t timebin, UShort_t value);
	  UShort_t  FilterGainNextSample(Int_t adc, UShort_t value);
	  UShort_t  FilterTailNextSample(Int_t adc, UShort_t value);

	  // tracklet calculation
	  void      AddHitToFitreg(Int_t adc, UShort_t timebin, UShort_t qtot, Short_t ypos, Int_t label[]);
	  void      CalcFitreg();
	  void      TrackletSelection();
	  void      FitTracklet();

	  Int_t         GetNHits() const { return fNHits; }
	  Bool_t        GetHit(Int_t index, Int_t &channel, Int_t &timebin, Int_t &qtot, Int_t &ypos, Float_t &y, Int_t &label) const;
	  TClonesArray* GetTrackletArray() const { return fTrackletArray; }

	  // data display
	  void      Print(Option_t* const option="") const;   // print stored data to stdout
	  void      Draw(Option_t* const option ="");         // draw data (ADC data, hits and tracklets)

  friend  std::ostream& operator<<(std::ostream &os, const AliTRDmcmSim &mcm); // data output using ostream (e.g. cout << mcm;)
  static  ostream&  Cfdat(ostream &os);                       // manipulator to activate cfdat output
  static  ostream&  Raw  (ostream &os);                       // manipulator to activate raw output
  static  ostream&  Text (ostream &os);                       // manipulator to activate text output

	  // PID
	  Int_t GetPID(Int_t q0, Int_t q1);
	  void PrintPidLutHuman();

	  // I/O
	  void PrintFitRegXml(ostream& os) const;
	  void PrintTrackletsXml(ostream& os) const;
	  void PrintAdcDatHuman(ostream& os) const;
	  void PrintAdcDatXml(ostream& os) const;
	  void PrintAdcDatDatx(ostream& os, Bool_t broadcast=kFALSE, Int_t timeBinOffset = -1) const;

  static  Bool_t ReadPackedConfig(AliTRDtrapConfig *cfg, Int_t hc, UInt_t *data, Int_t size);

  // DMEM addresses
  static const Int_t fgkDmemAddrLUTcor0       = 0xC02A;
  static const Int_t fgkDmemAddrLUTcor1       = 0xC028;
  static const Int_t fgkDmemAddrLUTnbins      = 0xC029;

  static const Int_t fgkDmemAddrLUTStart      = 0xC100; // LUT start address
  static const Int_t fgkDmemAddrLUTEnd        = 0xC3FF; // maximum possible end address for the LUT table
  static const Int_t fgkDmemAddrLUTLength     = 0xC02B; // address where real size of the LUT table is stored

  static const Int_t fgkDmemAddrTrackletStart = 0xC0E0; // Storage area for tracklets, start address
  static const Int_t fgkDmemAddrTrackletEnd   = 0xC0E3; // Storage area for tracklets, end address

  static const Int_t fgkDmemAddrDeflCorr      = 0xc022; // DMEM address of deflection correction
  static const Int_t fgkDmemAddrNdrift        = 0xc025; // DMEM address of Ndrift
  static const Int_t fgkDmemAddrDeflCutStart  = 0xc030; // DMEM start address of deflection cut
  static const Int_t fgkDmemAddrDeflCutEnd    = 0xc055; // DMEM end address of deflection cut

 protected:
	  Bool_t    CheckInitialized() const;           // Check whether the class is initialized

	  void      SetNTimebins(Int_t ntimebins);      // allocate data arrays corr. to the no. of timebins

 static const Int_t fgkFormatIndex;                     // index for format settings in stream

 static const Int_t fgkMaxTracklets = 4;                // maximum number of tracklet-words submitted per MCM (one per CPU)
 static const Int_t fgkAddDigits = 2;                   // additional digits used for internal representation of ADC data
	                                                // all internal data as after data control block (i.e. 12 bit), s. TRAP manual
 static const Int_t fgkNCPU = 4;                        // Number of CPUs in the TRAP
 static const Int_t fgkNHitsMC = 100;                   // maximum number of hits for which MC information is kept

 static const UShort_t fgkFPshifts[4];                  // shifts for pedestal filter

	  Bool_t    fInitialized;                       // memory is allocated if initialized
	  Int_t     fDetector;                          // Chamber ID
	  Int_t     fRobPos;                            // ROB Position on chamber
	  Int_t     fMcmPos;                            // MCM Position on chamber
	  Int_t     fRow;                               // Pad row number (0-11 or 0-15) of the MCM on chamber
	  Int_t     fNTimeBin;                          // Number of timebins currently allocated
	  Int_t   **fADCR;                              // Array with MCM ADC values (Raw, 12 bit)
	  Int_t   **fADCF;                              // Array with MCM ADC values (Filtered, 12 bit)
	  UInt_t   *fMCMT;                              // tracklet word for one mcm/trap-chip
	  TClonesArray *fTrackletArray;                 // Array of AliTRDtrackletMCM which contains MC information in addition to the tracklet word
	  Int_t    *fZSMap;                             // Zero suppression map (1 dimensional projection)

	  Int_t     fFitPtr[fgkNCPU];                   // pointer to the tracklet to be calculated by CPU i

          TString   fTrklBranchName;   	                // name of the tracklet branch to right to

	  // Parameter classes
	  AliTRDfeeParam    *fFeeParam;                 // FEE parameters
	  AliTRDtrapConfig  *fTrapConfig;               // TRAP config

	  AliTRDdigitsManager *fDigitsManager;          // pointer to digits manager used for MC label calculation
	  AliTRDarrayDictionary* fDict[3];              // pointers to label dictionaries


	  // internal filter registers
	  UInt_t*   fPedAcc;                            // Accumulator for pedestal filter
	  UInt_t*   fGainCounterA;                      // Counter for values above FGTA in the gain filter
	  UInt_t*   fGainCounterB;                      // Counter for values above FGTB in the gain filter
 	  UShort_t* fTailAmplLong;                      // Amplitude of the long component in the tail filter
 	  UShort_t* fTailAmplShort;                     // Amplitude of the short component in the tail filter

	  // hit detection
	  // individual hits can be stored as MC info
	  struct Hit_t {                                // Array of detected hits (only available in MC)
	  Hit_t() : fChannel(0), fTimebin(0), fQtot(0), fYpos(0) { fLabel[0] = 0; fLabel[1] = 0; fLabel[2] = 0; }
	    Int_t fChannel;                             // ADC channel of the hit
	    Int_t fTimebin;                             // timebin of the hit
	    Int_t fQtot;                                // total charge of the hit
	    Int_t fYpos;                                // calculated y-position
	    Int_t fLabel[3];                            // up to 3 labels (only in MC)
	  } fHits[fgkNHitsMC];
	  Int_t fNHits;                                 // Number of detected hits

	  // tracklet calculation
	  struct FitReg_t {                             // pointer to the 18 fit registers
	    Int_t  fNhits;                              // number of hits
	    UInt_t fQ0;                                 // charge accumulated in first window
	    UInt_t fQ1;                                 // charge accumulated in second window
	    UInt_t fSumX;                               // sum x
	    Int_t  fSumY;                               // sum y
	    UInt_t fSumX2;                              // sum x**2
	    UInt_t fSumY2;                              // sum y**2
	    Int_t  fSumXY;                              // sum x*y
	  } *fFitReg;

	  // Sort functions as in TRAP
	  void Sort2(UShort_t  idx1i, UShort_t  idx2i, UShort_t  val1i, UShort_t  val2i,
		     UShort_t * const idx1o, UShort_t * const idx2o, UShort_t * const val1o, UShort_t * const val2o) const;
	  void Sort3(UShort_t  idx1i, UShort_t  idx2i, UShort_t  idx3i,
		     UShort_t  val1i, UShort_t  val2i, UShort_t  val3i,
		     UShort_t * const idx1o, UShort_t * const idx2o, UShort_t * const idx3o,
		     UShort_t * const val1o, UShort_t * const val2o, UShort_t * const val3o);
	  void Sort6To4(UShort_t  idx1i, UShort_t  idx2i, UShort_t  idx3i, UShort_t  idx4i, UShort_t  idx5i, UShort_t  idx6i,
			UShort_t  val1i, UShort_t  val2i, UShort_t  val3i, UShort_t  val4i, UShort_t  val5i, UShort_t  val6i,
			UShort_t * const idx1o, UShort_t * const idx2o, UShort_t * const idx3o, UShort_t * const idx4o,
			UShort_t * const val1o, UShort_t * const val2o, UShort_t * const val3o, UShort_t * const val4o);
	  void Sort6To2Worst(UShort_t  idx1i, UShort_t  idx2i, UShort_t  idx3i, UShort_t  idx4i, UShort_t  idx5i, UShort_t  idx6i,
			     UShort_t  val1i, UShort_t  val2i, UShort_t  val3i, UShort_t  val4i, UShort_t  val5i, UShort_t  val6i,
			     UShort_t * const idx5o, UShort_t * const idx6o);

	  UInt_t AddUintClipping(UInt_t a, UInt_t b, UInt_t nbits) const;
	  // Add a and b (unsigned) with clipping to the maximum value representable by nbits

 private:
	  AliTRDmcmSim(const AliTRDmcmSim &m);             // not implemented
	  AliTRDmcmSim &operator=(const AliTRDmcmSim &m);  // not implemented

  static Bool_t fgApplyCut;               // apply cut on deflection length

  static Int_t fgAddBaseline;             // add baseline to the ADC values

  static Bool_t fgStoreClusters;          // whether to store all clusters in the tracklets

  ClassDef(AliTRDmcmSim,7)
};

std::ostream& operator<<(std::ostream& os, const AliTRDmcmSim& mcm);

#endif
