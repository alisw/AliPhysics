#ifndef ALIPHOSCONTABLEDB_H
#define ALIPHOSCONTABLEDB_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

/* History of cvs commits:
 *
 * $Log$
 */

//_________________________________________________________________________
//  Class for PHOS connection table, i.e. correspondence between 
//  "cable number" and PHOS AbsId number.
//                  
//*-- Author: D.Peressounko "RRC Kurchatov Institute"


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

  AliPHOSConTableDB(const AliPHOSConTableDB & cdb) ;   // cpy ctor
  virtual ~AliPHOSConTableDB() ; // dtor

  //Calculate table from known number of raws/columns 
  //assuming that prototype is situated in the center of 3 PHOS mod.
  void BuildDB(void) ;
  Int_t GetNchanels()const {return fNcrInProto ; } 
  Int_t GetNRaws()const{return fProtoRaws ;} 
  Int_t GetNColumns()const{return fProtoColumns ;}

  //set the number of columns in prototype
  void SetNCols(Int_t ncolumns){fProtoColumns = ncolumns ;}
  //Set the number of raw in prototype
  void SetNRaws(Int_t nraws){fProtoRaws = nraws ;}

  //Plot correspondance between Prototype Id and PHOS
  //Options are "Zoom" - only proto region is plotted
  //            "PHOS" - presents both PHOS and Proto ids
  void PlotProtoMap(Option_t * opt="Zoom") ; 

  //Transforms channel number in prototype into AbsId number in PHOS
  Int_t Raw2AbsId(Int_t raw)const ;

  //Transforms AbsId number in PHOS into channel number in prototype 
  Int_t AbsId2Raw(Int_t AbsId)const ;

  virtual void Print(const Option_t * = "") const ;


  AliPHOSConTableDB& operator=(const AliPHOSConTableDB& cdb) ;

private:
  AliPHOSGeometry * fGeom ;   //! PHOS geometry class

  Int_t     fProtoRaws ;        //  Parameters
  Int_t     fProtoColumns ;     //  used to calculate
  Int_t     fRawOffset ;        //  correspondance
  Int_t     fColOffset ;        //  map
  Int_t     fNcrInProto ;       //Number of channels in prototype
  Int_t     fMinAbsId ;         //Minimal AbsId, corresponding to some prototype cristall.
  Int_t     fMaxAbsId ;         //Maximal AbsId, corresponding to some prototype cristall
  TArrayS * fAbsIdMap ;         //Map of correspondance between Raw and PHOS ID
  TArrayS * fRawIdMap ;         //Map of correspondance between AbsId and Raw

  ClassDef(AliPHOSConTableDB,1)  // description 


};

#endif // AliPHOSCONTABLEDB_H
