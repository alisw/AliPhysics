#ifndef ALIZDCSATURATIONCALIB_H
#define ALIZDCSATURATIONCALIB_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//////////////////////////////////////////////////////////
//  class for ZDC calibration -> p-A high rate run      //
//////////////////////////////////////////////////////////

#include "TNamed.h"
#include "AliCDBEntry.h"

class AliZDC;

class AliZDCSaturationCalib: public TNamed {

 public:
  AliZDCSaturationCalib();
  AliZDCSaturationCalib(const char* name);
  AliZDCSaturationCalib(const AliZDCSaturationCalib &calibda);
  AliZDCSaturationCalib& operator= (const AliZDCSaturationCalib &calibda);
  virtual ~AliZDCSaturationCalib();
  void Reset();
  virtual void  Print(Option_t *) const; 
  //
  Float_t* GetZNASatCalib()   		const {return (float*)fZNASatCalibration;}
  Float_t GetZNASatCalib(int i)   	const {return fZNASatCalibration[i];}
  void 	SetZNASatCalib(Float_t* EnCalib);
  
  Float_t* GetZNCSatCalib()   		const {return (float*)fZNCSatCalibration;}
  Float_t GetZNCSatCalib(int i)   	const {return fZNCSatCalibration[i];}
  void 	SetZNCSatCalib(Float_t* EnCalib);
  
 protected:
  // 
  Float_t  fZNASatCalibration[4];	 // Coeff. for ZNA calibration
  Float_t  fZNCSatCalibration[4];	 // Coeff. for ZNC calibration
  //
  ClassDef(AliZDCSaturationCalib,1)    // ZDC calibration calibration data
};

#endif
