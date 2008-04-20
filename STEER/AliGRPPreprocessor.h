#ifndef ALIGRPPREPROCESSOR_H
#define ALIGRPPREPROCESSOR_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//-------------------------------------------------------------------------
//                          Class AliGRPPreprocessor
//                  Global Run Parameters (GRP) preprocessor
//
//    Origin: Panos Christakoglou, UOA-CERN, Panos.Christakoglou@cern.ch
//    Modified: Ernesto.Lopez.Torres@cern.ch  CEADEN-CERN
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
  virtual            ~AliGRPPreprocessor();
  
  static      Int_t   ReceivePromptRecoParameters(
                                  UInt_t run,
                                  const char* dbHost,
                                  Int_t dbPort,
                                  const char* dbName,
                                  const char* user,
                                  const char* password,
                                  const char *cdbRoot
                                 );

 protected:

  virtual      void   Initialize(Int_t run, UInt_t startTime, UInt_t endTime);
  
  virtual     UInt_t   Process(TMap* valueSet);

                TMap*  ProcessDaqLB();
              UInt_t   ProcessDaqFxs();
              UInt_t   ProcessDcsFxs();
               Int_t   ProcessDcsDPs(TMap* valueSet, TMap* grpmap);
   AliDCSSensorArray*  GetPressureMap(TMap *dcsAliasMap, AliDCSSensorArray *fPressure);
  
 private:
 
  static const Int_t   fgknDAQLbPar;            //! number of DAQ lb parameters
  static const Int_t   fgknDCSDP;               //! number of dcs dps
  static const char*   fgkDCSDataPoints[];      //! names of dcs dps
  static const char*   fgkLHCState[];           //! names of LHC States
  
  AliDCSSensorArray*   fPressure; //pressure array

                       AliGRPPreprocessor(const AliGRPPreprocessor&); // Not implemented
                       AliGRPPreprocessor& operator=(const AliGRPPreprocessor&); // Not implemented

  ClassDef(AliGRPPreprocessor, 1);
};

#endif
