#ifndef AliEMCALSENSORTEMPARRAY_H
#define AliEMCALSENSORTEMPARRAY_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */


///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  EMCAL calibration class for temperature sensors                            //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "TSystem.h"

#include "AliDCSSensorArray.h"
#include "AliEMCALSensorTemp.h"

class TTimeStamp;
class TMap;
class TGraph;
class TObjString;
class AliSplineFit;
class AliDCSSensor;

#include "TString.h"


class AliEMCALSensorTempArray : public AliDCSSensorArray {
 public:
  AliEMCALSensorTempArray();
  AliEMCALSensorTempArray(Int_t run);
  AliEMCALSensorTempArray(const char *fname,
			  const TString& amandaString = kAmandaString);
  AliEMCALSensorTempArray (UInt_t startTime, UInt_t endTime, TTree* confTree,
			   const TString& amandaString = kAmandaString);
  AliEMCALSensorTempArray(const AliEMCALSensorTempArray &c);
  virtual ~AliEMCALSensorTempArray();
  AliEMCALSensorTempArray &operator=(const AliEMCALSensorTempArray &c);
  void ReadSensors  (const char *dbEntry);
  AliEMCALSensorTemp* GetSensor (Int_t side, Int_t sector, Int_t num);
  AliEMCALSensorTemp* GetSensor (Int_t IdDCS);
  AliEMCALSensorTemp* GetSensor (Double_t x, Double_t y, Double_t z);
  Double_t GetTempGradientY(UInt_t timeSec, Int_t side);
  
 protected:
  
  ClassDef(AliEMCALSensorTempArray,1)       //  EMCAL calibration class for saved temperature sensor parameters
    
};

#endif
