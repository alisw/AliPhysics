#ifndef ALIPHOSCONTABLEDB_H
#define ALIPHOSCONTABLEDB_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//_________________________________________________________________________
//  Base Class for PHOS     
//                  
//*-- Author: (SUBATECH)


// --- ROOT system ---
#include "TNamed.h"
class TArrayS ;

// --- Standard library ---

// --- AliRoot header files ---
class AliPHOSGeometry ;

class AliPHOSConTableDB: public TNamed {

public:
  AliPHOSConTableDB() ;          // ctor
  AliPHOSConTableDB(const char * title) ;          // ctor
  virtual ~AliPHOSConTableDB() ; // dtor

  //Calculate table from known numbe of raws/columns 
  //assuming that prototype is situated in the center of 3 PHOS mod.
  void BuildDB(void) ;
  Int_t GetNchanels(){return fNcrInProto ; } 
  Int_t GetNRaws(){return fProtoRaws ;} 
  Int_t GetNColumns(){return fProtoColumns ;}

  //set the number of columns in prototype
  void SetNCols(Int_t ncolumns){fProtoColumns = ncolumns ;}
  //Set the number of raw in prototype
  void SetNRaws(Int_t nraws){fProtoRaws = nraws ;}

  //Plot correspondance between Prototype Id and PHOS
  //Options are "Zoom" - only proto region is plotted
  //            "PHOS" - presents both PHOS and Proto ids
  void PlotProtoMap(Option_t * opt="Zoom") ; 

  //Transforms channel number in prototype into AbsId number in PHOS
  Int_t Raw2AbsId(Int_t raw) ;

  //Transforms AbsId number in PHOS into channel number in prototype 
  Int_t AbsId2Raw(Int_t AbsId) ;

  virtual void Print(Option_t * option="") const ;

private:
  AliPHOSGeometry * fGeom ;   //!

  Int_t     fProtoRaws ;        //Parameters
  Int_t     fProtoColumns ;     //used to calculate
  Int_t     fRawOffset ;        //correspondance
  Int_t     fColOffset ;        //map
  Int_t     fNcrInProto ;    //Number of channels in prototype
  Int_t     fMinAbsId ;
  Int_t     fMaxAbsId ;
  TArrayS * fAbsIdMap ;      //Map of correspondance between Raw and PHOS ID
  TArrayS * fRawIdMap ;      //Map of correspondance between AbsId and Raw

  ClassDef(AliPHOSConTableDB,1)  // description 


};

#endif // AliPHOSCONTABLEDB_H
