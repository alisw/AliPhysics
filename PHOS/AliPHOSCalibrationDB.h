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
#include "TArrayF.h" 
#include "TString.h" 

// --- Standard library ---

// --- AliRoot header files ---
class AliPHOSConTableDB ;

class AliPHOSCalibrationDB:public TNamed {

public:
  AliPHOSCalibrationDB() ;          // ctor
  AliPHOSCalibrationDB(const char * database) ;
  virtual ~AliPHOSCalibrationDB() ; // dtor

  //Main method: calibrates if gains are known, otherwise - returns 0
  Float_t Calibrate(Int_t amp, Int_t absId)const ;

  //Read gains of pedestals from ascii file
  void ReadCalibrationParameters(const char * filename = "gains.dat",Option_t* opt = "gains") ;

  //Sets the same parameters for all channels  
  void SetAll(Float_t pedestal = 0, Float_t slope = 0.01) ; 

  //To know correspondance when reads list of gains from ascii file 
  void SetConTableDB(AliPHOSConTableDB * ctdb){fctdb = ctdb; }

  //Set parameters for particlular channel
  void SetParameters(Int_t AbsId, Float_t pedestal = 0, Float_t slope = 0.01)
    {if(fPedestals){fPedestals->AddAt(pedestal,AbsId-1) ; fSlopes->AddAt(slope,AbsId-1) ;} }

  //To be replaced in real DB when updating will really be necessary
  void Update(Int_t /*event*/,Int_t /*run*/){} 

  AliPHOSCalibrationDB & operator = (const AliPHOSCalibrationDB & ) ;
private:
  Int_t     fCurentRun ;       //! 
  Int_t     fNChannels ;
  TString   fFileName ;
  TArrayF * fPedestals ;
  TArrayF * fSlopes ;
  AliPHOSConTableDB * fctdb ;  //!

  ClassDef(AliPHOSCalibrationDB,1)  // description 

};

#endif // AliPHOSCALIBRATIONDB_H
