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
  AliPHOSCalibrManager(const AliPHOSCalibrManager & manager) : TNamed(manager) {
    // cpy ctor: no need
    // requested by the Coding Convention
    Fatal("cpy ctor", "not implemented") ;
  }
 
  virtual ~AliPHOSCalibrManager() ; // dtor
  static AliPHOSCalibrManager * GetInstance() ;
  static AliPHOSCalibrManager * GetInstance(const char * dbfilename ) ; 

  //To know correspondance when reads list of gains from ascii file 
  void SetConTableDB(AliPHOSConTableDB * ctdb){fctdb = ctdb; }

  //Read gains of pedestals from ascii file
  void ReadFromASCII(AliPHOSCalibrationData & data,const char * filename = "gains.dat") ;

  void ReadFromRoot(AliPHOSCalibrationData &data,Int_t run) ;

  void WriteData(AliPHOSCalibrationData *data) ;

  AliPHOSCalibrManager & operator = (const AliPHOSCalibrManager & right) ;

private:
  AliPHOSCalibrManager(const char* filename) ;          

private:
  TString   fFileName ;        //Name of file with calibration data
  AliPHOSConTableDB * fctdb ;  //! Connection table of PHOS
  static AliPHOSCalibrManager * fgCaMa ; // pointer to the unique instance of singleton
	    
 
  ClassDef(AliPHOSCalibrManager,1)  // description 

};

#endif // AliPHOSCALIBRMANAGER_H
