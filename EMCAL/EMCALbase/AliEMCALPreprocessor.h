#ifndef ALIEMCALPREPROCESSOR_H
#define ALIEMCALPREPROCESSOR_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */
///////////////////////////////////////////////////////////////////////////////
// Class AliEMCALPreprocessor
///////////////////////////////////////////////////////////////////////////////


#include "AliPreprocessor.h"

class AliEMCALSensorTempArray;
class TEnv;

class AliEMCALPreprocessor : public AliPreprocessor {

 public:
  
  AliEMCALPreprocessor(); //! ctor
  AliEMCALPreprocessor(AliShuttleInterface* shuttle); //! overloaded ctor
  AliEMCALPreprocessor(const AliEMCALPreprocessor &); //! copy ctor
  AliEMCALPreprocessor& operator = (const  AliEMCALPreprocessor &source); //! assignment operator
  virtual ~AliEMCALPreprocessor();//! dtor

 protected:

  virtual void Initialize(Int_t run, UInt_t startTime, UInt_t endTime);//!
  virtual UInt_t Process(TMap* dcsAliasMap);//!
  UInt_t  MapTemperature(TMap* dcsAliasMap);//!
  UInt_t  MapTriggerConfig(TMap* dcsAliasMap);//!
  UInt_t  ExtractPedestals(Int_t sourceFXS);//!
  UInt_t  ExtractSignal(Int_t sourceFXS);//!

 private:
  TEnv                   *fConfEnv;  // Preprocessor configuration map
  AliEMCALSensorTempArray  *fTemp;     // CDB class for temperature sensors
  Bool_t                 fConfigOK;  // Identify succesful reading of OCDB Config
    
  ClassDef(AliEMCALPreprocessor,1);

};

#endif
