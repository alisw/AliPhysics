#ifndef ALIADDigitizer_H
#define ALIADDigitizer_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
  
///_________________________________________________________________________
///
///  Class for making Digits in AD 
///_________________________________________________________________________   


// --- Standard library ---

// --- AliRoot header files ---

#include "AliDigitizer.h"
#include "AliADConst.h"

class TClonesArray;
class TF1;
class AliDigitizationInput;
class AliCDBManager;
class AliCDBStorage;
class AliADCalibData;
class AliAD;

class AliADDigitizer: public AliDigitizer {

public:
   enum DigiTask_t { 
     kHits2Digits, 
     kHits2SDigits
   };
  
   AliADDigitizer() ;                       // default constructor
   AliADDigitizer(AliAD *AD, DigiTask_t task);         // constructor
   AliADDigitizer(AliDigitizationInput* digInput); // constructor
   virtual ~AliADDigitizer() ;              // destructor

   virtual Bool_t Init();
   virtual void   Digitize(Option_t* option=0);

   void DigitizeHits();
   void DigitizeSDigits();
   void WriteDigits(AliLoader *loader);
   void WriteSDigits(AliLoader *loader);
   void ReadSDigits();

   void AddDigit(Int_t pmnumber, Float_t time, Float_t width, Bool_t integrator, Short_t *chargeADC, Int_t *labels);
   void AddSDigit(Int_t pmnumber, Int_t nbins, Float_t *charges, Int_t *labels);
   TClonesArray* DigitsArray(); 
   TClonesArray* SDigitsArray(); 
   void ResetDigits();
						
   AliADCalibData *GetCalibData() const;

   TF1*   GetSignalShape() const { return fSignalShape; }
   TF1*   GetPMResponse() const { return fPMResponse; }
   TF1*   GetSinglePhESpectrum() const { return fSinglePhESpectrum; }
   double SignalShape(double *x, double *par);
   double PMResponse(double *x, double *par);
   double SinglePhESpectrum(double *x, double *par);

 protected:
 
   AliADCalibData *fCalibData;  //! calibration data
 
 private:
 
   AliADDigitizer(const AliADDigitizer& /*digitizer*/); 
      
   AliADDigitizer& operator = (const AliADDigitizer& /*digitizer*/); 

   Int_t    fNdigits;                //! Number of digits
   TClonesArray *fDigits;            //! List of digits
   
   TF1*     fSignalShape;            // function which describes the PMT signal shape
   TF1*     fPMResponse;             // function which describes the PM time response
   TF1*     fSinglePhESpectrum;      // function which describes the single ph.e. PM response

   Float_t  fAdc[16][kNClocks];      //! Container for ADC samples
   Float_t  fLeadingTime[16];        //! Leading time container
   Float_t  fTimeWidth[16];          //! Time width container
   Float_t  fAdcPedestal[16][2];     //! Pedestals, one per integrator
   Float_t  fAdcSigma[16][2];        //! Sigma of pedestals
   Float_t  fPmGain[16];             //! PMT gains
   Int_t    fNBins[16];              //! Number of bins in fTime container
   Int_t    fNBinsLT[16];            //! Number of bins in fTime container (match window only)
   Float_t  fBinSize[16];            //! Bin size in fTime container
   Float_t  fHptdcOffset[16];        //! HPTDC time offsets channel by channel
   Float_t  fClockOffset[16];        //! Clock offsets channel by channel

   Float_t *fTime[16];               //! Main container used in digitization
   Int_t    fLabels[16][3];          //! Container for MC labels
   Bool_t   fEvenOrOdd;              //! Choise of integrator in central ADC sample

   DigiTask_t fTask;                 //! The task (to be) executed by the digitizer
   AliAD  *fAD;                //! Pointer to AliDetector object
   ClassDef(AliADDigitizer,1)     // digitizer for AD

};

#endif // AliADDigitizer_H
