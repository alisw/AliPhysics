#ifndef ALITPCCALIBDB_H
#define ALITPCCALIBDB_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// Class providing the calibration parameters by accessing the CDB           //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////


class AliTPCTransform;
class AliTPCExB;
#include "TObject.h"
class AliTPCCalPad;
class AliTPCSensorTempArray;
class AliDCSSensorArray;
class AliCDBEntry;
class AliTPCParam;
class AliTPCAltroMapping;
class AliTPCClusterParam;
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
  AliTPCTransform* GetTransform() {return fTransform;}
  AliTPCExB* GetExB() {return fExB;}
  AliTPCCalPad* GetPadGainFactor() {return fPadGainFactor;}
  AliTPCCalPad* GetDedxGainFactor() {return fDedxGainFactor;}
  AliTPCCalPad* GetPadTime0() {return fPadTime0;}
  AliTPCCalPad* GetPadNoise() {return fPadNoise;}
  AliTPCCalPad* GetPedestals() {return fPedestals;}
  AliTPCSensorTempArray* GetTemperature() {return fTemperature;}
  AliTPCParam*  GetParameters(){return fParam;}
  AliTPCAltroMapping ** GetMapping(){ return fMapping;}
  AliTPCClusterParam *GetClusterParam(){ return fClusterParam;}
  //
  static void     CreateObjectList(const Char_t *filename, TObjArray *calibObjects);
  static void MakeTree(const char * fileName, TObjArray * array, const char * mapFileName = 0, AliTPCCalPad* outlierPad = 0, Float_t ltmFraction = 0.9);
  
  //
protected:
  void         Update();  //update entries
  AliCDBEntry* GetCDBEntry(const char* cdbPath);   
  Long64_t        fRun;         // current run number
  AliTPCTransform *fTransform;      // object responsible for spacial corrections
  AliTPCExB *fExB;              // ExB correction factor
//  AliCDBStorage* fLocator;      // Storage locator retrieved from AliCDBManager
  //
  // calibration parameters per pad
  //
  AliTPCCalPad* fPadGainFactor;   // Gain calibration entry
  AliTPCCalPad* fDedxGainFactor;   // Gain calibration entry - for dEdx
  AliTPCCalPad* fPadTime0;        // Time0 calibration entry
  AliTPCCalPad* fPadNoise;        // Noise calibration entry
  AliTPCCalPad* fPedestals;       // Pedestal calibration entry
  AliTPCSensorTempArray* fTemperature; // Temperature calibration entry
  AliTPCAltroMapping **fMapping;   // Altro mapping   
  //
  //
  AliTPCParam * fParam;           // TPC parameters
  AliTPCClusterParam * fClusterParam;  // TPC cluster error, shape and Q parameterization
  //
  static AliTPCcalibDB* fgInstance;  // singleton control
  static Bool_t       fgTerminated;  // termination control 
  ClassDef(AliTPCcalibDB, 0)
 private:
   AliTPCcalibDB (const AliTPCcalibDB& );
   AliTPCcalibDB& operator= (const AliTPCcalibDB& );
};


#endif
