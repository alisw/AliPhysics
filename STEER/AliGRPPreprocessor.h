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

class AliGRPPreprocessor: public AliPreprocessor {
 public:
  AliGRPPreprocessor();
  AliGRPPreprocessor(AliShuttleInterface* shuttle);
  virtual ~AliGRPPreprocessor();
  
 protected:

  virtual void Initialize(Int_t run, UInt_t startTime, UInt_t endTime);
  
  virtual UInt_t Process(TMap* valueSet);
  
  ClassDef(AliGRPPreprocessor, 0);
};

#endif
