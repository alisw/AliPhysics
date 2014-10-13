#ifndef ALITPCTEMPMAP_H
#define ALITPCTEMPMAP_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */


///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  TPC calibration class for temperature maps and tendencies                //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "TSystem.h"

class TGraph;
class TGraph2D;
class TLinearFitter;
class TString;
class AliTPCSensorTempArray;
class TTimeStamp;

class AliTPCTempMap : public TNamed  {
 public:
  AliTPCTempMap(AliTPCSensorTempArray *SensorsDCS);
  AliTPCTempMap(const AliTPCTempMap &c);   
  virtual ~AliTPCTempMap();
  AliTPCTempMap &operator=(const AliTPCTempMap &c);
  virtual void Copy (TObject &c) const;
  TLinearFitter *GetLinearFitter(Int_t type, Int_t side, UInt_t timeSec);
  TLinearFitter *GetLinearFitter(Int_t type, Int_t side, TTimeStamp& stamp);
  //
  Double_t GetTempGradientY(UInt_t timeSec, Int_t side);
  TGraph2D *GetTempMapsViaSensors(Int_t type, Int_t side, UInt_t timeSec);
  TGraph *MakeGraphGradient(Int_t axis, Int_t side, Int_t nPoints);

  Double_t GetTemperature(Double_t x, Double_t y, Double_t z, UInt_t timeSec);
  Double_t GetTemperature(Double_t x, Double_t y, Double_t z, TTimeStamp &stamp);
  Bool_t  IsOK(Float_t value);
 protected:
  
  AliTPCSensorTempArray *fTempArray;   // Array of Sensors (initialized in Constructor)
  TString fStringFEsimulation; // Placeholder for file of FiniteElement 
                               // Simulation under ideal conditions - not existing yet

 private:

  AliTPCTempMap(const char *fname);

  ClassDef(AliTPCTempMap,2)      //  

};

#endif
