#ifndef ALIPHOSCALIBRATIONDB_H
#define ALIPHOSCALIBRATIONDB_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//_________________________________________________________________________    
//                  
//*-- Author: D.Peressounko (RRC KI & SUBATECH)


// --- ROOT system ---
#include "TNamed.h"
#include "TString.h" 

// --- Standard library ---

// --- AliRoot header files ---
#include "AliPHOSCalibrationData.h" 

class AliPHOSCalibrationDB:public TNamed {

public:
  AliPHOSCalibrationDB() ;          // ctor
  AliPHOSCalibrationDB(const char * database) ;
  virtual ~AliPHOSCalibrationDB() ; // dtor

  //Main method: calibrates if gains are known, otherwise - returns 0
  Float_t Calibrate(Int_t amp, Int_t absId)const ;

  //Get calibration parameters using AliPHOSCalibrManager
  void GetParameters(void) ; 

  AliPHOSCalibrationDB & operator = (const AliPHOSCalibrationDB & ) ;
private:

  AliPHOSCalibrationData fPedestals ;
  AliPHOSCalibrationData fGains ;

  ClassDef(AliPHOSCalibrationDB,1)  // description 

};

#endif // AliPHOSCALIBRATIONDB_H
