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

class TClonesArray;
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

   void AddDigit(Int_t PMnumber, Float_t adc, Float_t time);
   void ResetDigit();
   void GetCollisionMode();
   void GetCollisionMode(Int_t collisionMode, Float_t beamEnergy) 
                        {fCollisionMode=collisionMode; fBeamEnergy=beamEnergy;}
						
   AliVZEROCalibData *GetCalibData() const;
   Int_t GetPMNumber(Int_t cell) const;

 protected:
 
   AliVZEROCalibData *fCalibData;  //! calibration data
 
 private:
 
   AliVZERODigitizer(const AliVZERODigitizer& /*digitizer*/); 
      
   AliVZERODigitizer& operator = (const AliVZERODigitizer& /*digitizer*/); 
  
   Float_t  fPhotoCathodeEfficiency; // Photocathode efficiency
   Float_t  fPMVoltage ;             // Photomultiplier voltage
   Float_t  fPMGain;                 // Photomultiplier gain

   Int_t    fNdigits;                //! Number of digits
   TClonesArray *fDigits;            //! List of digits
   
   Int_t    fCollisionMode;          // =0->p-p, =1->A-A
   Float_t  fBeamEnergy;	     // beam energy
   
   ClassDef(AliVZERODigitizer,2)     // digitizer for VZERO

};

#endif // AliVZERODigitizer_H
