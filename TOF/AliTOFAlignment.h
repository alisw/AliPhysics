#ifndef ALITOFALIGNMENT_H
#define ALITOFALIGNMENT_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//////////////////////////////////////////////////////////////////
//  class for TOF Alignment::                                   //
//////////////////////////////////////////////////////////////////

#include "TTask.h"
#include "TGeoManager.h"
#include "TGeoVolume.h"
#include "TGeoBBox.h"
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
  virtual void WriteOnCDBforDC();
  virtual void ReadFromCDBforDC();
  virtual void BuildGeomForSurvey(); //Generate ideal volumes for Survey
  virtual void AlignFromSurvey(); //Prototype to get Align Obj from Survey
  virtual void InsertMisAlignment( Float_t *mis); //To test align. from Survey
  TObjArray * GetTOFAlignArray() const {return fTOFAlignObjArray;}

private:


  // This is the extended box (the TOF master sensitive volume+services)
  static const Double_t fgkXsizeTOF;  // x size of the TOF ext. volume, cm
  static const Double_t fgkYsizeTOF;  // y size of the TOF ext. volume, cm
  static const Double_t fgkZsizeTOF;  // x size of the TOF ext. volume, cm
  static const Double_t fgkRorigTOF;  // Radius of the TOF ext. volume, cm

  // Four fiducial marks on SM, expressed in local coordinates (origin=center of TOF SM)
  // They are positioned at x=+/- 38 cm, y=11.2, z=+/- 456.94 cm

  static const Double_t fgkXFM; //x FM pos, cm 
  static const Double_t fgkYFM; //y FM pos, cm
  static const Double_t fgkZFM; //z FM pos, cm

  // This is the z size of the TOF master sensitive volume (it is shorter)
  static const Double_t fgkZsizeTOFSens; //cm

  Int_t fNTOFAlignObj; // Number of Alignable Objects
  TGeoManager *fTOFmgr; //pointer to TGeoManager
  TObjArray *fTOFAlignObjArray; // Pointer to the TOF alignable objects
  TGeoHMatrix* fTOFMatrixId[18]; // Ideal Matrices of TOF Volumes in the GRS
  Float_t fTOFSurveyFM[18][4][3]; // Temporary container for Survey data xyz (4/TOF SuperModule)
  ClassDef(AliTOFAlignment,2)   // TOF Alignment 
};

#endif
