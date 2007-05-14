#ifndef AliTPCSensorPressureArray_H
#define AliTPCSensorPressureArray_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */


///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  TPC calibration class for pressure sensors                               //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "TSystem.h"

#include "AliDCSSensorArray.h"
#include "AliTPCSensorPressure.h"

class TTimeStamp;
class TMap;
class TGraph;
class TObjString;
class AliSplineFit;
class AliDCSSensor;

#include "TString.h"


class AliTPCSensorPressureArray : public AliDCSSensorArray {
 public:
  AliTPCSensorPressureArray();
  AliTPCSensorPressureArray(Int_t prevRun); 
  AliTPCSensorPressureArray(const char *fname);
  AliTPCSensorPressureArray (UInt_t startTime, UInt_t endTime, const char *filepath=".");
  AliTPCSensorPressureArray(const AliTPCSensorPressureArray &c);   
  virtual ~AliTPCSensorPressureArray();
  AliTPCSensorPressureArray &operator=(const AliTPCSensorPressureArray &c);
  virtual void Copy (TObject &c) const;
  void SetGraph     (TMap *map);
  void MakeSplineFit(TMap *map);
  void ReadSensors  (const char *fname);
  const char* GetAmandaString() { return fAmandaString.Data(); }
  void SetAmandaString(const char* string) {fAmandaString=string;}
  TMap* ExtractDCS  (TMap *dcsMap);
  AliTPCSensorPressure* GetSensor (Int_t type, Int_t side, Int_t sector, Int_t num);
  AliTPCSensorPressure* GetSensor (Int_t IdDCS);
  AliTPCSensorPressure* GetSensor (Double_t x, Double_t y, Double_t z);

 protected:

  TString fAmandaString;        //! Amanda string to identify temperature entries
  ClassDef(AliTPCSensorPressureArray,1)       //  TPC calibration class for parameters which are saved per pad

};

#endif
