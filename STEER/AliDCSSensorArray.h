#ifndef AliDCSSensorArray_H
#define AliDCSSensorArray_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */


///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  Calibration class for DCS sensors                                        //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "TNamed.h"
#include "TMap.h"
#include "TObjString.h"
#include "AliCDBManager.h"
#include "AliCDBEntry.h"
#include "AliDCSValue.h"
#include "AliDCSSensor.h"

class TGraph;
class TTimeStamp;


class AliDCSSensorArray : public TNamed {
 public:
  AliDCSSensorArray();
  AliDCSSensorArray(Int_t prevRun, const char* dbEntry);
  AliDCSSensorArray(const AliDCSSensorArray &c);   
  virtual ~AliDCSSensorArray();
  AliDCSSensorArray &operator=(const AliDCSSensorArray &c);
  virtual void Copy (TObject &c) const;
  void SetStartTime (const TTimeStamp& start) { fStartTime = start; }
  void SetEndTime   (const TTimeStamp& end) { fEndTime = end; }
  void SetGraph     (TMap *map, const char* amandaString);
  void MakeSplineFit(TMap *map, const char* amandaString);
  TMap* ExtractDCS  (TMap *dcsMap, const char* amandaString);
  TGraph* MakeGraph (TObjArray *valueSet);
  Double_t GetValue  (UInt_t timeSec, Int_t sensor);
  AliDCSSensor* GetSensor (Int_t IdDCS);
  AliDCSSensor* GetSensor (Double_t x, Double_t y, Double_t z);

 protected:
  Int_t fFirstSensor;             // DCS ID of first sensor
  Int_t fLastSensor;              // DCS ID of last sensor
  TTimeStamp  fStartTime;         // start time for measurements in this entry
  TTimeStamp  fEndTime;           // end time for measurements in this entry
  TClonesArray *fSensors;         // Array of sensors

  ClassDef(AliDCSSensorArray,1)       //  TPC calibration class for parameters which are saved per pad

};

#endif
