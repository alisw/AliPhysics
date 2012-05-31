#ifndef AliTPCSensorTempArray_H
#define AliTPCSensorTempArray_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */


///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  TPC calibration class for temperature sensors                            //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "TSystem.h"

#include "AliDCSSensorArray.h"
#include "AliTPCSensorTemp.h"

class TTimeStamp;
class TMap;
class TGraph;
class TObjString;
class AliSplineFit;
class AliDCSSensor;

#include "TString.h"


class AliTPCSensorTempArray : public AliDCSSensorArray {
 public:
  AliTPCSensorTempArray();
  AliTPCSensorTempArray(Int_t run);
  AliTPCSensorTempArray(const char *fname,
                        const TString& amandaString = kAmandaString);
  AliTPCSensorTempArray (UInt_t startTime, UInt_t endTime, TTree* confTree,
                         const TString& amandaString = kAmandaString);
  AliTPCSensorTempArray(const AliTPCSensorTempArray &c);
  virtual ~AliTPCSensorTempArray();
  AliTPCSensorTempArray &operator=(const AliTPCSensorTempArray &c);
  void ReadSensors  (const char *dbEntry);
  AliTPCSensorTemp* GetSensor (Int_t type, Int_t side, Int_t sector, Int_t num);
  AliTPCSensorTemp* GetSensor (Int_t IdDCS);
  AliTPCSensorTemp* GetSensor (Double_t x, Double_t y, Double_t z);
  Double_t GetTempGradientY(UInt_t timeSec, Int_t side);

 protected:

  ClassDef(AliTPCSensorTempArray,1)       //  TPC calibration class for parameters which are saved per pad

};

#endif
