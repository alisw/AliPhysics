#ifndef AliT0CalibLatency_H
#define AliT0CalibLatency_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

////////////////////////////////////////////////
//  class for T0 calibration                 //
////////////////////////////////////////////////

#include "TNamed.h"

class AliT0CalibLatency: public TNamed {

 public:
  AliT0CalibLatency();
  AliT0CalibLatency(const char* name);
  AliT0CalibLatency(const AliT0CalibLatency &calibda);
  AliT0CalibLatency& operator= (const AliT0CalibLatency &calibda);
  virtual ~AliT0CalibLatency();
  
  virtual void  Print(Option_t* option= "") const; 
  
  Float_t   GetLatencyL1() const {return fLatencyL1;}
  void      SetLatencyL1(Float_t lat) {fLatencyL1 = lat;}
  Float_t   GetLatencyL1A() const {return fLatencyL1A;}
  void      SetLatencyL1A(Float_t lat) {fLatencyL1A = lat;}
  Float_t   GetLatencyL1C() const {return fLatencyL1C;}
  void      SetLatencyL1C(Float_t lat) {fLatencyL1C = lat;}
  Float_t   GetLatencyHPTDC() const {return fLatencyHPTDC;}
  void      SetLatencyHPTDC(Float_t lat) {fLatencyHPTDC = lat;}
  
 protected:
  Float_t   fLatencyL1;         //Latency L1
  Float_t   fLatencyL1A;        //Latency L1 for OrA
  Float_t   fLatencyL1C;        //Latency L1 for orC
  Float_t   fLatencyHPTDC;      //Latency HPTDC


  ClassDef(AliT0CalibLatency,1)    // T0 Sensor Calibration data
};

typedef AliT0CalibLatency AliSTARTCalibLatency; // for backward compatibility

#endif

