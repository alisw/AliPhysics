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
class AliGRPObject;
class AliSplineFit;

class AliGRPPreprocessor: public AliPreprocessor {
 public:

	enum DP {kLHCState = 0, kL3Polarity, kDipolePolarity, kLHCLuminosity, kBeamIntensity, 
		 kL3Current, kDipoleCurrent, 
		 kL3_BSF17_H1, kL3_BSF17_H2, kL3_BSF17_H3, kL3_BSF17_Temperature, 
		 kL3_BSF4_H1, kL3_BSF4_H2, kL3_BSF4_H3, kL3_BSF4_Temperature, 
		 kL3_BKF17_H1, kL3_BKF17_H2, kL3_BKF17_H3, kL3_BKF17_Temperature, 
		 kL3_BKF4_H1, kL3_BKF4_H2, kL3_BKF4_H3, kL3_BKF4_Temperature, 
		 kL3_BSF13_H1, kL3_BSF13_H2, kL3_BSF13_H3, kL3_BSF13_Temperature,
		 kL3_BSF8_H1, kL3_BSF8_H2, kL3_BSF8_H3, kL3_BSF8_Temperature,
		 kL3_BKF13_H1, kL3_BKF13_H2, kL3_BKF13_H3, kL3_BKF13_Temperature,
		 kL3_BKF8_H1, kL3_BKF8_H2, kL3_BKF8_H3, kL3_BKF8_Temperature,
		 kDipole_Inside_H1, kDipole_Inside_H2, kDipole_Inside_H3, kDipole_Inside_Temperature,
		 kDipole_Outside_H1, kDipole_Outside_H2, kDipole_Outside_H3, kDipole_Outside_Temperature,
                 kCavernTemperature, kCavernAtmosPressure, kSurfaceAtmosPressure};

	enum DP_HallProbes { 
		 k_HP_L3_BSF17_H1= 0 , k_HP_L3_BSF17_H2, k_HP_L3_BSF17_H3, k_HP_L3_BSF17_Temperature, 
		 k_HP_L3_BSF4_H1, k_HP_L3_BSF4_H2, k_HP_L3_BSF4_H3, k_HP_L3_BSF4_Temperature, 
		 k_HP_L3_BKF17_H1, k_HP_L3_BKF17_H2, k_HP_L3_BKF17_H3, k_HP_L3_BKF17_Temperature, 
		 k_HP_L3_BKF4_H1, k_HP_L3_BKF4_H2, k_HP_L3_BKF4_H3, k_HP_L3_BKF4_Temperature, 
		 k_HP_L3_BSF13_H1, k_HP_L3_BSF13_H2, k_HP_L3_BSF13_H3, k_HP_L3_BSF13_Temperature,
		 k_HP_L3_BSF8_H1, k_HP_L3_BSF8_H2, k_HP_L3_BSF8_H3, k_HP_L3_BSF8_Temperature,
		 k_HP_L3_BKF13_H1, k_HP_L3_BKF13_H2, k_HP_L3_BKF13_H3, k_HP_L3_BKF13_Temperature,
		 k_HP_L3_BKF8_H1, k_HP_L3_BKF8_H2, k_HP_L3_BKF8_H3, k_HP_L3_BKF8_Temperature,
		 k_HP_Dipole_Inside_H1, k_HP_Dipole_Inside_H2, k_HP_Dipole_Inside_H3, k_HP_Dipole_Inside_Temperature,
		 k_HP_Dipole_Outside_H1, k_HP_Dipole_Outside_H2, k_HP_Dipole_Outside_H3, k_HP_Dipole_Outside_Temperature};

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

                AliGRPObject*  ProcessDaqLB();
              UInt_t   ProcessDaqFxs();
              UInt_t   ProcessDcsFxs();
               Int_t   ProcessDcsDPs(TMap* valueSet, AliGRPObject* grpobj);
               Int_t   ProcessLHCDPs(TMap* valueSet, AliGRPObject* grpobj);
               Int_t   ProcessL3DPs(TMap* valueSet, AliGRPObject* grpobj);
               Int_t   ProcessDipoleDPs(TMap* valueSet, AliGRPObject* grpobj);
               Int_t   ProcessEnvDPs(TMap* valueSet, AliGRPObject* grpobj);
               Int_t   ProcessHPDPs(TMap* valueSet, AliGRPObject* grpobj);
	       //               Int_t   ProcessDcsDPs(TMap* valueSet, TMap* grpmap);
   AliDCSSensorArray*  GetPressureMap(TMap *dcsAliasMap);
   AliSplineFit* GetSplineFit(TObjArray *array, const TString& stringID);
   //AliSplineFit* GetSplineFit(TMap* mapDCS, const TString& stringID);
   TString ProcessChar(TObjArray *array);
   Char_t ProcessBool(TObjArray *array);
   Float_t ProcessInt(TObjArray *array);
   Float_t ProcessUInt(TObjArray *array);
   Float_t* ProcessFloatAll(TObjArray* array);

 private:
 
  static const Int_t   fgknDAQLbPar;            //! number of DAQ lb parameters
  static const Int_t   fgknDCSDP;               //! number of dcs dps
  static const char*   fgkDCSDataPoints[];      //! names of dcs dps
  static const char*   fgkLHCState[];           //! names of LHC States
  static const char*   fgkDCSDataPoints_HallProbes[];      //! names of dcs dps for Hall Probes
  static const Int_t   fgknDCSDP_HallProbes;           //! names of LHC States for Hall Probes

  AliDCSSensorArray*   fPressure; //pressure array

                       AliGRPPreprocessor(const AliGRPPreprocessor&); // Not implemented
                       AliGRPPreprocessor& operator=(const AliGRPPreprocessor&); // Not implemented

  ClassDef(AliGRPPreprocessor, 0);
};

#endif
