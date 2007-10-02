#ifndef ALITOFALIGNMENT_H
#define ALITOFALIGNMENT_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//////////////////////////////////////////////////////////////////
//  class for TOF Alignment::                                   //
//////////////////////////////////////////////////////////////////
#include <TObject.h>
#include <TString.h>
#include "TTask.h"
#include "TGeoManager.h"
#include "TGeoVolume.h"
#include "TGeoBBox.h"
#include "TGeoTrd1.h"
#include "TGeoPhysicalNode.h"
#include "TGeoNode.h"


class TObjArray;

class AliTOFAlignment :public TTask{

  enum {kMaxAlignObj=2000}; //maximal number of the TOF Alignable Objects

 public:

  AliTOFAlignment(); 
  AliTOFAlignment(const AliTOFAlignment &t); //Copy Ctor 
  AliTOFAlignment& operator=(const AliTOFAlignment &source); // Assignment Operator
  virtual ~AliTOFAlignment();
  virtual void WriteParOnCDB(Char_t *sel, Int_t minrun, Int_t maxrun);
  virtual void ReadParFromCDB(Char_t *sel, Int_t nrun);
  virtual void WriteSimParOnCDB(Char_t *sel, Int_t minrun, Int_t maxrun);
  virtual void ReadSimParFromCDB(Char_t *sel, Int_t nrun);
  virtual void Smear(Float_t *tr=0, Float_t *rot=0); // create a set of AlignObj for TOF
  virtual void Align(Float_t *tr=0, Float_t *rot=0); // create a set of AlignObj for TOF
  TObjArray * GetTOFAlignArray() const {return fTOFAlignObjArray;}

  //methods for survey
  virtual void WriteOnCDBforDC();
  virtual void ReadFromCDBforDC();
  virtual void BuildGeomForSurvey();      //Build ideal geometry (FTOA in BTOF)
  virtual void InsertMisAlignment( Float_t *mis); //To test align. from Survey
  virtual void MakeDefData(const Int_t nf,TString namefiles[]); //Combines survey data from different files
  virtual void WriteCombData(const Char_t *nomefile, Int_t option); //Write combined data on a file 
  virtual void ReadSurveyDataAndAlign(); //Read survey data and call the right Alignement procedure  
  virtual void WriteSimSurveyData(const Char_t *nomefile); // Write sim data in standard format

private:
  
  static const Double_t fgkRorigTOF;  //Radius of the TOF ext. volume, cm
  static const Double_t fgkX1BTOF;    //x1 size of BTOF
  static const Double_t fgkX2BTOF;    //x2 size of BTOF
  static const Double_t fgkYBTOF;     //y size of BTOF
  static const Double_t fgkZBTOF;     //z size of BTOF
  
  // Four fiducial marks on SM, expressed in local coordinates (origin=center of TOF SM)
  // They are positioned at x=+/- 38 cm, y=+/- 457.3 cm, z=11.2 cm

  static const Double_t fgkXFM; //x pos of FM in BTOF, cm 
  static const Double_t fgkYFM; //y pos of FM in BTOF, cm
  static const Double_t fgkZFM; //z pos of FM in BTOF, cm
  
  Int_t fNTOFAlignObj;           //Number of Alignable Objects
  TGeoManager *fTOFmgr;          //Pointer to TGeoManager
  TObjArray *fTOFAlignObjArray;  //Pointer to the TOF alignable objects
  TGeoHMatrix* fTOFMatrixId[18]; //Ideal Matrices of TOF Volumes in the GRS
  Float_t fCombFMData[72][6];    //Combined survey data
  Int_t fNFMforSM[18][5];        //Number of FM for each SM

  void AlignFromSurveyABC(Int_t iSM);  //From Survey data of FM ABC, derive the needed transformations to get the Alignment Objects. 
  void AlignFromSurveyABD(Int_t iSM);  //From Survey data of FM ABD, derive the needed transformations to get the Alignment Objects. 
  void AlignFromSurveyACD(Int_t iSM);  //From Survey data of FM ACD, derive the needed transformations to get the Alignment Objects. 
  void AlignFromSurveyBCD(Int_t iSM);  //From Survey data of FM BCD, derive the needed transformations to get the Alignment Objects. 

  ClassDef(AliTOFAlignment,2)   // TOF Alignment 
};

#endif
