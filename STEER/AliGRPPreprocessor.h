#ifndef ALIGRPPREPROCESSOR_H
#define ALIGRPPREPROCESSOR_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//-------------------------------------------------------------------------
//                          Class AliGRPPreprocessor
//                  Global Run Parameters (GRP) preprocessor
//
//    Origin: Panos Christakoglou, UOA-CERN, Panos.Christakoglou@cern.ch
//-------------------------------------------------------------------------



//////////////////////////////////////////////////////////////////////////
//                                                                      //
//                        AliGRPPreprocessor                            //
//                                                                      //
//           Implementation of the GRP preprocessor                     //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "AliPreprocessor.h"

class TList;
class TString;
class AliDCSSensorArray;

class AliGRPPreprocessor: public AliPreprocessor {
 public:
  AliGRPPreprocessor(AliShuttleInterface* shuttle);
  virtual ~AliGRPPreprocessor();
  
 protected:

  virtual void Initialize(Int_t run, UInt_t startTime, UInt_t endTime);
  
  virtual UInt_t Process(TMap* valueSet);

  TList *ProcessDaqLB();
  UInt_t ProcessDaqFxs();
  UInt_t ProcessDcsFxs();
  TList *ProcessDcsDPs(TMap* valueSet);
  AliDCSSensorArray *GetPressureMap(TMap *dcsAliasMap, AliDCSSensorArray *fPressure);

 private:
  static const char* fgkDCSDataPoints[12]; //! names of dcs dps
  AliDCSSensorArray *fPressure; //pressure array

  AliGRPPreprocessor(const AliGRPPreprocessor&); // Not implemented
  AliGRPPreprocessor& operator=(const AliGRPPreprocessor&); // Not implemented

  ClassDef(AliGRPPreprocessor, 0);
};

#endif
