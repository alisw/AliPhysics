#ifndef ALIPHOSCALIBRMANAGER_H
#define ALIPHOSCALIBRMANAGER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */



//_________________________________________________________________________    
//                  
//*-- Author: D.Peressounko (RRC KI & SUBATECH)


// --- ROOT system ---
#include "TNamed.h"
#include "TString.h" 
class TArrayF ; 
// --- Standard library ---

// --- AliRoot header files ---
class AliPHOSConTableDB ;
class AliPHOSCalibrationData ;

class AliPHOSCalibrManager:public TNamed {

public:
  AliPHOSCalibrManager() ;          // ctor
  AliPHOSCalibrManager(const AliPHOSCalibrManager & manager);
 
  virtual ~AliPHOSCalibrManager() ; // dtor
  static AliPHOSCalibrManager * GetInstance() ;
  static AliPHOSCalibrManager * GetInstance(const char * filename,const char * kind = "root" ) ; 

  void GetParameters(AliPHOSCalibrationData &data) ; 
  
  void SetConTable(AliPHOSConTableDB * ct){fctdb = ct ;}

  void WriteData(AliPHOSCalibrationData &data) ;

  AliPHOSCalibrManager & operator = (const AliPHOSCalibrManager & right) ;

private:
  //Read gains of pedestals from ascii file
  void ReadFromASCII(AliPHOSCalibrationData & data) ;

  void ReadFromRoot(AliPHOSCalibrationData &data) ;

  AliPHOSCalibrManager(const char* filename,const char * kind = "root") ;          

private:
  TString   fFileName ;        //Name of file with calibration data
  Int_t     fInputKind;        //Kind of input to read/write data
  AliPHOSConTableDB * fctdb ;  //Connection table used to read from ASCII file
  
  static AliPHOSCalibrManager * fgCaMa ; // pointer to the unique instance of singleton
	    
 
  ClassDef(AliPHOSCalibrManager,1)  // description 

};

#endif // AliPHOSCALIBRMANAGER_H
