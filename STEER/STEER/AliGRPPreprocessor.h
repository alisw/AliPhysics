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
class TObjArray;
class AliDCSSensorArray;
class AliGRPObject;
class AliSplineFit;
class AliLHCClockPhase;

class AliGRPPreprocessor: public AliPreprocessor {
 public:

	enum DP {kL3Polarity = 0, kDipolePolarity,  
		 kL3Current, kDipoleCurrent, 
		 kL3bsf17H1, kL3bsf17H2, kL3bsf17H3, kL3bsf17Temperature, 
		 kL3bsf4H1, kL3bsf4H2, kL3bsf4H3, kL3bsf4Temperature, 
		 kL3bkf17H1, kL3bkf17H2, kL3bkf17H3, kL3bkf17Temperature, 
		 kL3bkf4H1, kL3bkf4H2, kL3bkf4H3, kL3bkf4Temperature, 
		 kL3bsf13H1, kL3bsf13H2, kL3bsf13H3, kL3bsf13Temperature,
		 kL3bsf8H1, kL3bsf8H2, kL3bsf8H3, kL3bsf8Temperature,
		 kL3bkf13H1, kL3bkf13H2, kL3bkf13H3, kL3bkf13Temperature,
		 kL3bkf8H1, kL3bkf8H2, kL3bkf8H3, kL3bkf8Temperature,
		 kDipoleInsideH1, kDipoleInsideH2, kDipoleInsideH3, kDipoleInsideTemperature,
		 kDipoleOutsideH1, kDipoleOutsideH2, kDipoleOutsideH3, kDipoleOutsideTemperature,
                 kCavernTemperature, kCavernAtmosPressure, kSurfaceAtmosPressure, 
		 kCavernAtmosPressure2};

	enum DPHallProbes { 
		 khpL3bsf17H1=0, khpL3bsf17H2, khpL3bsf17H3, khpL3bsf17Temperature, 
		 khpL3bsf4H1, khpL3bsf4H2, khpL3bsf4H3, khpL3bsf4Temperature, 
		 khpL3bkf17H1, khpL3bkf17H2, khpL3bkf17H3, khpL3bkf17Temperature, 
		 khpL3bkf4H1, khpL3bkf4H2, khpL3bkf4H3, khpL3bkf4Temperature, 
		 khpL3bsf13H1, khpL3bsf13H2, khpL3bsf13H3, khpL3bsf13Temperature,
		 khpL3bsf8H1, khpL3bsf8H2, khpL3bsf8H3, khpL3bsf8Temperature,
		 khpL3bkf13H1, khpL3bkf13H2, khpL3bkf13H3, khpL3bkf13Temperature,
		 khpL3bkf8H1, khpL3bkf8H2, khpL3bkf8H3, khpL3bkf8Temperature,
		 khpDipoleInsideH1, khpDipoleInsideH2, khpDipoleInsideH3, khpDipoleInsideTemperature,
		 khpDipoleOutsideH1, khpDipoleOutsideH2, khpDipoleOutsideH3, khpDipoleOutsideTemperature};

                      AliGRPPreprocessor(AliShuttleInterface* shuttle);
  virtual            ~AliGRPPreprocessor();
  
  static      Int_t   ReceivePromptRecoParameters(
                                  UInt_t run,
                                  const char* dbHost,
                                  Int_t dbPort,
                                  const char* dbName,
                                  const char* user,
                                  const char* password,
                                  const char *cdbRoot,
				  TString &gdc
                                 );

 protected:

  virtual      void   Initialize(Int_t run, UInt_t startTime, UInt_t endTime);
  
  virtual     UInt_t   Process(TMap* valueSet);

               Int_t   ProcessDaqLB(AliGRPObject* grpobj);
              UInt_t   ProcessDaqFxs();
              UInt_t   ProcessDqmFxs();
              UInt_t   ProcessSPDMeanVertex();
              UInt_t   ProcessLHCData(AliGRPObject* grpobj);
              UInt_t   ProcessDcsFxs(TString partition="", TString detector="");
               Int_t   ProcessDcsDPs(TMap* valueSet, AliGRPObject* grpobj);
               Int_t   ProcessL3DPs(const TMap* valueSet, AliGRPObject* grpobj);
               Int_t   ProcessDipoleDPs(const TMap* valueSet, AliGRPObject* grpobj);
               Int_t   ProcessEnvDPs(TMap* valueSet, AliGRPObject* grpobj);
               Int_t   ProcessHPDPs(const TMap* valueSet, AliGRPObject* grpobj);
	       //               Int_t   ProcessDcsDPs(TMap* valueSet, TMap* grpmap);
   AliDCSSensorArray*  GetPressureMap(TMap *dcsAliasMap);
   AliSplineFit* GetSplineFit(const TObjArray *array, const TString& stringID);
   //AliSplineFit* GetSplineFit(TMap* mapDCS, const TString& stringID);
   TString ProcessChar(const TObjArray *array);
   Char_t ProcessBool(const TObjArray *array, Bool_t &change);
   Float_t ProcessInt(const TObjArray *array);
   Float_t ProcessUInt(const TObjArray *array);
   Float_t* ProcessFloatAll(const TObjArray* array);
   Float_t* ProcessFloatAllMagnet(const TObjArray* array, Int_t indexDP, Bool_t &isZero);
   Bool_t GetDAQStartEndTimeOk() const {return fdaqStartEndTimeOk;}
   void SetDAQStartEndTimeOk(Bool_t daqStartEndTimeOk) {fdaqStartEndTimeOk = daqStartEndTimeOk;}

   Float_t ProcessEnergy(TObjArray* array, Double_t timeStart);
   AliLHCClockPhase* ProcessLHCClockPhase(TObjArray *beam1phase,TObjArray *beam2phase,
					  Double_t timeEnd);

 private:
 
  static const Int_t   fgknDAQLbPar;            //! number of DAQ lb parameters for PHYSICS runs
  static const Int_t   fgknDCSDP;               //! number of dcs dps
  static const char*   fgkDCSDataPoints[];      //! names of dcs dps
  static const char*   fgkDCSDataPointsHallProbes[];      //! names of dcs dps for Hall Probes
  static const Int_t   fgknDCSDPHallProbes;           //! number of Hall Probes

  AliDCSSensorArray*   fPressure; //pressure array

  AliGRPPreprocessor(const AliGRPPreprocessor&); // Not implemented
  AliGRPPreprocessor& operator=(const AliGRPPreprocessor&); // Not implemented

  Float_t fmaxFloat; // maximum float accepted
  Float_t fminFloat; // minimum float accepted
  Double_t fmaxDouble; // maximum double accepted
  Double_t fminDouble; // minimum double accepted
  Int_t fmaxInt; // maximum int accepted
  Int_t fminInt; // minimum int accepted
  UInt_t fmaxUInt; // maximum uint accepted
  UInt_t fminUInt; // minimum uint accepted
  Bool_t fdaqStartEndTimeOk; // flag to set whether the DAQ_time_start/end fields are set

  TObjArray* ffailedDPs; //TObjArray of failed DPs names

  static const Int_t   fgknLHCDP;               //! number of lhc dps
  static const char*   fgkLHCDataPoints[];      //! names of lhc dps
  static const Int_t   fgkDCSDPHallTopShift;    //! shift from the top to the Hall Probes from the list of DCS DPs names
  static const Int_t   fgkDCSDPNonWorking;    //! Number of non working DCS DPs 

  ClassDef(AliGRPPreprocessor, 0);
};

#endif
