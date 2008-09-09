//-*- Mode: C++ -*-
// $Id$

#ifndef ALIHLTPREDICTIONPROCESSORTPC_H
#define ALIHLTPREDICTIONPROCESSORTPC_H
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               *

/** @file   AliHLTPredictionProcessorTPC.h
    @author Haavard Helstrup
    @date   
    @brief  
*/

#include "AliHLTPredictionProcessorInterface.h"
#include "TTree.h"
#include "AliTPCSensorTempArray.h"

class TMap;

class TObjArray;

class AliHLTPredictionProcessorTPC : public AliHLTPredictionProcessorInterface
{
 public:
   AliHLTPredictionProcessorTPC(const char* detector, 
                                AliHLTPendolino* pendolino);
   virtual ~AliHLTPredictionProcessorTPC();
   
   UInt_t makePrediction (Bool_t doPrediction);
   void Initialize(Int_t run, UInt_t startTime, UInt_t endTime);
   UInt_t Process (TMap* dcsAliasMap);
   UInt_t ExtractTemperature(TMap* dcsAliasMap);
   TMap* produceTestData(TString aliasName);
   virtual Bool_t ProcessDCS();
   
  private: 
    Bool_t                 fConfigOK;     // Identify succesful reading of OCDB Config
    TTree*                 fConfTreeTemp; // TTree holding temperature configuration  
    AliTPCSensorTempArray* fTemp;         // Array holding temperature readings
    Bool_t                 fPredict;      // Indicates whether predictions should be made
    Int_t                  fRun;          // run number
    UInt_t                 fStartTime;    // start time
    UInt_t                 fEndTime;      // end time

    AliHLTPredictionProcessorTPC
             (const AliHLTPredictionProcessorTPC& predictPro); // Disabled Copy Constructor

    AliHLTPredictionProcessorTPC& 
             operator=(const AliHLTPredictionProcessorTPC& rhs); // Disabled Assignment Operator


  ClassDef(AliHLTPredictionProcessorTPC,2)
};

inline Bool_t AliHLTPredictionProcessorTPC::ProcessDCS() {
	return true;
}

#endif

