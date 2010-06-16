#ifndef ALIVZERODigitizer_H
#define ALIVZERODigitizer_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
  
///_________________________________________________________________________
///
///  Class for making Digits in VZERO 
///_________________________________________________________________________   


// --- Standard library ---

// --- AliRoot header files ---

#include "AliDigitizer.h"

#include "AliVZEROConst.h"

class TClonesArray;
class TF1;
class AliRunDigitizer;
class AliCDBManager;
class AliCDBStorage;
class AliVZEROCalibData;

class AliVZERODigitizer: public AliDigitizer {

 public:

   AliVZERODigitizer() ;                       // constructor
   AliVZERODigitizer(AliRunDigitizer *manager);// constructor
   virtual ~AliVZERODigitizer() ;              // destructor
  
   virtual Bool_t Init();
   virtual void   Exec(Option_t* option=0);

   void AddDigit(Int_t PMnumber, Float_t time, Float_t width, Bool_t integrator, Short_t *chargeADC, Int_t *labels);
   void ResetDigit();
						
   AliVZEROCalibData *GetCalibData() const;

   TF1*   GetSignalShape() const { return fSignalShape; }
   TF1*   GetPMResponse() const { return fPMResponse; }
   TF1*   GetSinglePhESpectrum() const { return fSinglePhESpectrum; }
   double SignalShape(double *x, double *par);
   double PMResponse(double *x, double *par);
   double SinglePhESpectrum(double *x, double *par);

   Int_t  Cell2Pmt(Int_t cell) const;

 protected:
 
   AliVZEROCalibData *fCalibData;  //! calibration data
 
 private:
 
   AliVZERODigitizer(const AliVZERODigitizer& /*digitizer*/); 
      
   AliVZERODigitizer& operator = (const AliVZERODigitizer& /*digitizer*/); 
  
   Float_t  fPhotoCathodeEfficiency; // Photocathode efficiency

   Int_t    fNdigits;                //! Number of digits
   TClonesArray *fDigits;            //! List of digits
   
   TF1*     fSignalShape;            // function which describes the PMT signal shape
   TF1*     fPMResponse;             // function which describes the PM time response
   TF1*     fSinglePhESpectrum;      // function which describes the single ph.e. PM response

   Float_t  fAdc[64][kNClocks];      //! Container for ADC samples
   Float_t  fLeadingTime[64];        //! Leading time container
   Float_t  fTimeWidth[64];          //! Time width container
   Float_t  fAdcPedestal[64][2];     //! Pedestals, one per integrator
   Float_t  fAdcSigma[64][2];        //! Sigma of pedestals
   Float_t  fPmGain[64];             //! PMT gains
   Int_t    fNBins[64];              //! Number of bins in fTime container
   Int_t    fNBinsLT[64];            //! Number of bins in fTime container (match window only)
   Float_t  fBinSize[64];            //! Bin size in fTime container
   Float_t  fHptdcOffset[64];        //! HPTDC time offsets channel by channel

   Float_t *fTime[64];               //! Main container used in digitization
   
   ClassDef(AliVZERODigitizer,5)     // digitizer for VZERO

};

#endif // AliVZERODigitizer_H
