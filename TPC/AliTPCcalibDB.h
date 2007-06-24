#ifndef ALITPCCALIBDB_H
#define ALITPCCALIBDB_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// Class providing the calibration parameters by accessing the CDB           //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////


#include "TObject.h"
class AliTPCCalPad;
class AliTPCSensorTempArray;
class AliDCSSensorArray;
class AliCDBEntry;
class AliTPCParam;
//class AliCDBStorage;

class AliTPCcalibDB : public TObject
{
 public: 
  static AliTPCcalibDB* Instance();
  AliTPCcalibDB();
  virtual ~AliTPCcalibDB();
  static void Terminate();
  void   SetRun(Long64_t run);   
  //
  AliTPCCalPad* GetPadGainFactor() {return fPadGainFactor;}
  AliTPCCalPad* GetPadTime0() {return fPadTime0;}
  AliTPCCalPad* GetPadPRFWidth() {return fPadPRFWidth;}
  AliTPCCalPad* GetPadNoise() {return fPadNoise;}
  AliTPCCalPad* GetPedestals() {return fPedestals;}
  AliTPCSensorTempArray* GetTemperature() {return fTemperature;}
  AliDCSSensorArray* GetPressure() {return fPressure;}
  AliTPCParam*  GetParameters(){return fParam;}
  static void     CreateObjectList(const Char_t *filename, TObjArray *calibObjects);
  static void MakeTree(const char * fileName, TObjArray * array, const char * mapFileName = 0, AliTPCCalPad* outlierPad = 0, Float_t ltmFraction = 0.9);
  
  //
protected:
  void         Update();  //update entries
  AliCDBEntry* GetCDBEntry(const char* cdbPath);   
  Long64_t        fRun;         // current run number            
//  AliCDBStorage* fLocator;      // Storage locator retrieved from AliCDBManager
  //
  // calibration parameters per pad
  //
  AliTPCCalPad* fPadGainFactor;
  AliTPCCalPad* fPadTime0;
  AliTPCCalPad* fPadPRFWidth;
  AliTPCCalPad* fPadNoise;
  AliTPCCalPad* fPedestals;
  AliTPCSensorTempArray* fTemperature;
  AliDCSSensorArray *fPressure;
  //
  //
  AliTPCParam * fParam;
  //
  static AliTPCcalibDB* fgInstance;
  static Bool_t       fgTerminated;
  ClassDef(AliTPCcalibDB, 0)
 private:
   AliTPCcalibDB (const AliTPCcalibDB& org);
   AliTPCcalibDB& operator= (const AliTPCcalibDB& rhs);
};


#endif
