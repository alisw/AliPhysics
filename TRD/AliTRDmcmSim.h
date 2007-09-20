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

class AliTRDfeeParam;
class AliTRDSimParam;
class AliTRDcalibDB;
class AliTRDgeometry;

class AliTRDmcmSim : public TObject {

 public:

                    AliTRDmcmSim();
                    AliTRDmcmSim(const AliTRDmcmSim &m);
  virtual          ~AliTRDmcmSim();
  AliTRDmcmSim      &operator=(const AliTRDmcmSim &m);

  virtual void      Copy(TObject &m) const;

          void      Init( Int_t cha, Int_t rob, Int_t mcm );   // Initialize MCM by the position parameters
          void      SetData(Int_t iadc, Int_t *adc);           // Set ADC data with array 
          void      SetData(Int_t iadc, Int_t it, Int_t adc ); // Set ADC data
          void      SetDataPedestal(Int_t iadc );              // Fill ADC data with pedestal values

          Int_t     GetChaId()  { return fChaId;  };           // Returns Chamber ID (0-539)
          Int_t     GetRobPos() { return fRobPos; };           // Returns ROB position (0-7)
          Int_t     GetMcmPos() { return fMcmPos; };           // Returns MCM position (0-17) (16,17 are mergers)
          Int_t     GetRow()    { return fRow; };              // Returns Row number on chamber where the MCM is sitting
	  Int_t     GetCol( Int_t iadc );                      // Get corresponding column (0-143) from for ADC channel iadc = [0:20]
	  // for the ADC/Col mapping, see: http://wiki.kip.uni-heidelberg.de/ti/TRD/index.php/Image:ROB_MCM_numbering.pdf

	  Int_t     ProduceRawStream( UInt_t *buf, Int_t bufsize ); // Produce raw data stream from this MCM
	  void      Filter();                                  // Apply digital filters for existing data
	  void      ZSMapping();                               // Do ZS mapping for existing data
	  void      DumpData( char *f, char *target );         // Dump data stored (only for debugging)

 protected:

          Bool_t    fInitialized;                       // Status whether the class is initialized or not
	  Int_t     fChaId;                             // Chamber ID
	  Int_t     fSector;                            // Sector number of the supermodule
	  Int_t     fStack;                             // Chamber stack ID
	  Int_t     fLayer;                             // Chamber layer ID
          Int_t     fRobPos;                            // ROB Position on chamber
          Int_t     fMcmPos;                            // MCM Position on chamber
	  Int_t     fNADC;                              // Number of ADC (usually 21)
	  Int_t     fNTimeBin;                          // Number of Timebins (variable)
          Int_t     fRow;                               // Pad row number (0-11 or 0-15) of the MCM on chamber
          Int_t   **fADCR;                              // Array with MCM ADC values (Raw)
          Int_t   **fADCF;                              // Array with MCM ADC values (Filtered)
          Int_t   **fZSM;                               // Zero suppression map
          Int_t    *fZSM1Dim;                           // Zero suppression map (1 dimensional projection)

	  // parameter classes
	  AliTRDfeeParam *fFeeParam;
	  AliTRDSimParam *fSimParam;
	  AliTRDcalibDB  *fCal;
	  AliTRDgeometry *fGeo;

	  Bool_t    CheckInitialized();                 // Check whether the class is initialized
	  void      FilterPedestal();                   // Apply pedestal filter
	  void      FilterGain();                       // Apply gain filter
	  void      FilterTail();                       // Apply tail filter

	  // Here are the several sub functions used only in FilterTail
	  void      FilterSimDeConvExpA(Int_t *source, Double_t *target, Int_t n, Int_t nexp);
	  void      FilterSimDeConvExpD(Int_t *source, Int_t *target, Int_t n, Int_t nexp);
	  void      FilterSimDeConvExpMI(Int_t *source, Double_t *target, Int_t n);
	  void      FilterSimTailMakerSpline(Double_t *ampin, Double_t *ampout, Double_t lambda, Int_t n);
	  void      FilterSimTailCancelationMI(Double_t *ampin, Double_t *ampout, Double_t norm, Double_t lambda, Int_t n);

  ClassDef(AliTRDmcmSim,3)

};

#endif
