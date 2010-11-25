#ifndef ALI_TRD_PREPROCESSOR_H
#define ALI_TRD_PREPROCESSOR_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

////////////////////////////////////////////////////////////////////////////
//                                                                        //
// TRD preprocessor for the database SHUTTLE                              //
//                                                                        //
////////////////////////////////////////////////////////////////////////////

#include "AliPreprocessor.h"

class TMap;
class AliTRDCalDCS;
class AliTRDCalROC;
class AliTRDCalSingleChamberStatus;

class AliTRDPreprocessor : public AliPreprocessor
{

 public:

  AliTRDPreprocessor(AliShuttleInterface *shuttle);
  AliTRDPreprocessor(const AliTRDPreprocessor &org);
  virtual ~AliTRDPreprocessor();


 protected:

  virtual void    Initialize(Int_t run, UInt_t startTime, UInt_t endTime);
  virtual UInt_t  Process(TMap *dcsAliasMap);

          Bool_t  ExtractHalfChamberStatusDAQ();
          Bool_t  ExtractPedestals();
          Bool_t  ExtractDriftVelocityDAQ();
          Bool_t  ExtractHLT();
          Bool_t  ProcessDCS();
          Bool_t  ProcessDCS(TMap *dcsAliasMap);
	  AliTRDPreprocessor& operator = (const AliTRDPreprocessor& rhs);

 private:
	  
	  AliTRDCalDCS* fCalDCSObjSOR;    // 
	  AliTRDCalDCS* fCalDCSObjEOR;    // 

	  Bool_t  fVdriftHLT;             // HLT Vdrift
	  UInt_t  ProcessDCSConfigData(); // process DCS configuration

	  Bool_t AreThereDataPedestal(AliTRDCalSingleChamberStatus * const calROCStatus, Bool_t second);
	  void   SetDefaultStatus(AliTRDCalSingleChamberStatus &calROCStatus, Bool_t second);
	  void   SetStatus(AliTRDCalSingleChamberStatus &calROCStatus, AliTRDCalSingleChamberStatus *calROCStatusPrevious,Bool_t second);
	  void   SetDefaultNoise(AliTRDCalROC &calROCNoise, Bool_t second);
	  void   SetNoise(AliTRDCalROC &calROCNoise, AliTRDCalROC *calROCNoisePrevious, Bool_t second);

	  ClassDef(AliTRDPreprocessor,1)          // The SHUTTLE preprocessor for TRD

};
#endif

