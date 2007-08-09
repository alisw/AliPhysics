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
class AliTRDgeometry;

class AliTRDmcmSim : public TObject {

 public:

  AliTRDmcmSim();
  virtual         ~AliTRDmcmSim();

  //PH  virtual void      Copy(TObject &m) const;

          void      Init( Int_t cha_id, Int_t rob_pos, Int_t mcm_pos ); // Initialize MCM by the position parameters
          void      SetData(Int_t iadc, Int_t *adc);           // Set ADC data with array 
          void      SetData(Int_t iadc, Int_t it, Int_t adc ); // Set ADC data one by one
          void      SetDataPedestal(Int_t iadc );              // Fill ADC data with pedestal values

          Int_t     GetChaId()  { return fChaId;  };    // Returns Chamber ID (0-539)
          Int_t     GetRobPos() { return fRobPos; };    // Returns ROB position (0-7)
          Int_t     GetMcmPos() { return fMcmPos; };    // Returns MCM position (0-17) (16,17 are mergers)
          Int_t     GetRow()    { return fRow; };       // Returns Row number on chamber where it is sitting
	  Int_t     GetCol( Int_t iadc );               // Get corresponding column (0-143) from adc ID (= 0-20) 
                    // for the ADC/Col mapping, see: http://wiki.kip.uni-heidelberg.de/ti/TRD/index.php/Image:ROB_MCM_numbering.pdf

	  Int_t     ProduceRawStream( UInt_t *buf, Int_t bufsize );    // Produce raw data stream from this MCM
	  void      Filter();                           //  Apply digital filters
	  void      ZSMapping();                        //  Do ZS mapping

 protected:

          Bool_t    fInitialized;                       // Status

	  AliTRDfeeParam *fFeeParam;
	  AliTRDgeometry *fGeo;

	  Int_t     fChaId;                             // Chamber ID
	  Int_t     fSector;                            // Sector number of the supermodule
	  Int_t     fStack;                             // Chamber stack ID
	  Int_t     fLayer;                             // Chamber layer ID
          Int_t     fRobPos;                            // ROB Position
          Int_t     fMcmPos;                            // MCM Position
	  Int_t     fNADC;                              // Number of ADC (usually 21)
	  Int_t     fNTimeBin;                          // Number of Timebin (valiable)
          Int_t     fRow;                               // Pad row number (0-11 or 0-15)
          Int_t     fColOfADCbeg;                       // Column corresponds to ADC channel 0
          Int_t     fColOfADCend;                       // Column corresponds to ADC channel [fNADC-1] (usually 20)
          Int_t   **fADCR;                              // Array with MCM ADC values (Raw)
          Int_t   **fADCF;                              // Array with MCM ADC values (Filtered)
          Int_t   **fZSM;                               // Zero Suppression Map
          Int_t    *fZSM1Dim;                           // Zero Suppression (1 dimensional projection)

	  void      FilterPedestal();                   // Apply pedestal filter
	  void      FilterGain();                       // Apply gain filter
	  void      FilterTail();                       // Apply tail filter

	  void      FilterSimDeConvExpA(Int_t *source, Double_t *target, Int_t n, Int_t nexp);
	  void      FilterSimDeConvExpD(Int_t *source, Int_t *target, Int_t n, Int_t nexp);
	  void      FilterSimDeConvExpMI(Int_t *source, Double_t *target, Int_t n);
	  void      FilterSimTailMakerSpline(Double_t *ampin, Double_t *ampout, Double_t lambda, Int_t n);
	  void      FilterSimTailCancelationMI(Double_t *ampin, Double_t *ampout, Double_t norm, Double_t lambda, Int_t n);

  ClassDef(AliTRDmcmSim,3)

};

#endif
