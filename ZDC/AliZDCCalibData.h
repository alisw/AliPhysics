#ifndef ALIZDCCALIBDATA_H
#define ALIZDCCALIBDATA_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

////////////////////////////////////////////////
//  class for ZDC calibration                 //
////////////////////////////////////////////////

#include "TNamed.h"
#include "TH1.h"
#include "AliZDC.h"
#include "AliCDBEntry.h"

class AliZDCCalibData: public TNamed {

 public:
  AliZDCCalibData();
  AliZDCCalibData(const char* name);
  AliZDCCalibData(const AliZDCCalibData &calibda);
  AliZDCCalibData& operator= (const AliZDCCalibData &calibda);
  virtual ~AliZDCCalibData();
  void Reset();
  virtual void  Print(Option_t *) const; 
  //
  Float_t  GetMeanPed(Int_t channel)   	const {return fMeanPedestal[channel];}
  Float_t* GetMeanPed()   const {return (float*)fMeanPedestal;}
  Float_t  GetEnCalib(Int_t channel)	const {return fEnCalibration[channel];}
  Float_t* GetEnCalib()   const {return (float*)fEnCalibration;}
  //
  void     SetMeanPed(Float_t val, Int_t channel) {fMeanPedestal[channel]=val;}
  void     SetMeanPed(Float_t* MeanPed);
  void 	   SetEnCalib(Float_t val, Int_t channel) {fEnCalibration[channel]=val;}
  void 	   SetEnCalib(Float_t* EnCalib);
  //
  void     PrepHistos();
  TH1F*    GetHistMeanPed() const {return fHistMeanPed;}
  void     CleanHistos();

 protected:
  Float_t  fMeanPedestal[47];	// Mean pedestal values
  Float_t  fEnCalibration[4];	// Coeff. for energy calibration (4 different ZDC's?)
  TH1F*    fHistMeanPed;        //! histos for drawing mean pedestals
  //
  ClassDef(AliZDCCalibData,1)    // ZDC Sensor Calibration data
};

#endif
