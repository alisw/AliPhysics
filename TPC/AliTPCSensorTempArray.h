#ifndef AliTPCSensorTempArray_H
#define AliTPCSensorTempArray_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */


///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  TPC calibration class for temperature sensors                            //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////


#include "AliTPCSensorTemp.h"
#include "AliDCSSensorArray.h"

class TTimeStamp;
class TMap;
class TGraph;
class TObjString;
class AliSplineFit;
class AliDCSSensor;

class AliTPCSensorTempArray : public AliDCSSensorArray {
 public:
  AliTPCSensorTempArray();
  AliTPCSensorTempArray(Int_t prevRun); 
  AliTPCSensorTempArray(const char *fname);
  AliTPCSensorTempArray (UInt_t startTime, UInt_t endTime);
  AliTPCSensorTempArray(const AliTPCSensorTempArray &c);   
  virtual ~AliTPCSensorTempArray();
  AliTPCSensorTempArray &operator=(const AliTPCSensorTempArray &c);
  virtual void Copy (TObject &c) const;
  void SetGraph     (TMap *map);
  void MakeSplineFit(TMap *map);
  void ReadSensors  (const char *fname);
  TMap* ExtractDCS  (TMap *dcsMap);

 protected:

  ClassDef(AliTPCSensorTempArray,1)       //  TPC calibration class for parameters which are saved per pad

};

#endif
