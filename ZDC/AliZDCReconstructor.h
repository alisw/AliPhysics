#ifndef ALIZDCRECONSTRUCTOR_H
#define ALIZDCRECONSTRUCTOR_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

/////////////////////////////////////////////
//                                   	   //
//       class for ZDC reconstruction      //
//                                   	   //
/////////////////////////////////////////////

#include "AliReconstructor.h"
#include "AliCDBManager.h"
#include "AliCDBStorage.h"
#include "AliZDCPedestals.h"
#include "AliZDCEnCalib.h"
#include "AliZDCTowerCalib.h"
#include "AliZDCRecoParam.h"
#include "AliZDCRecoParampp.h"
#include "AliZDCRecoParamPbPb.h"
#include "AliLog.h"

class AliLoader;

class AliZDCReconstructor: public AliReconstructor {
public:
  AliZDCReconstructor();
  virtual ~AliZDCReconstructor();

  virtual void   Init();
  virtual Bool_t HasDigitConversion() const {return kFALSE;};
  
  virtual void Reconstruct(TTree* digitsTree, TTree* clustersTree) const; 
  virtual void Reconstruct(AliRawReader* rawReader, TTree* clustersTree) const;

  virtual void FillESD(TTree* /*digitsTree*/, TTree* clustersTree, AliESDEvent* esd) const 
  	        {FillZDCintoESD(clustersTree, esd);}
  virtual void FillESD(AliRawReader* /*rawReader*/, TTree* clustersTree, AliESDEvent* esd) const 
  	        {FillZDCintoESD(clustersTree, esd);}
  
  // parameter settings for reconstruction
  void SetRecoMode(Int_t recoMode, Float_t beamEnergy) 
  		  {fRecoMode=recoMode; fBeamEnergy=beamEnergy;}
  static void SetRecoParam(AliZDCRecoParam * param){fRecoParam = param;}
  
  Int_t   GetRecoMode() {return fRecoMode;}
  Float_t GetBeamEnergy() {return fBeamEnergy;}
  
  static const AliZDCRecoParam* GetRecoParam() {return fRecoParam;}

  void  SetPedSubMode(Int_t pedsubMode) {fPedSubMode=pedsubMode;}
  Int_t GetPedSubMode() {return fPedSubMode;}
  
  void    SetSignalThreshold(Float_t val) {fSignalThreshold=val;}
  Float_t GetSignalThreshold() {return fSignalThreshold;}
  
  // OCDB objects for reconstruction
  AliCDBStorage       *SetStorage(const char* uri);
  AliZDCPedestals     *GetPedestalData() const; 
  AliZDCEnCalib       *GetEnergyCalibData() const; 
  AliZDCTowerCalib    *GetTowerCalibData() const; 
  AliZDCRecoParampp   *GetppRecoParamFromOCDB() const;  
  AliZDCRecoParamPbPb *GetPbPbRecoParamFromOCDB() const;  
  
  void WritePbPbRecoParamInOCDB() const;

private:
  AliZDCReconstructor(const AliZDCReconstructor&);
  AliZDCReconstructor& operator =(const AliZDCReconstructor&);

  void   ReconstructEventpp(TTree *clustersTree, 
  	    Float_t* ZN1ADCCorr, Float_t* ZP1ADCCorr, Float_t* ZN2ADCCorr, Float_t* ZP2ADCCorr,
	    Float_t* ZEM1ADCCorr, Float_t* ZEM2ADCCorr, Float_t* PMRef1, Float_t* PMRef2,
	    Bool_t channelsOff, Bool_t chUnderflow, Bool_t chOverflow) const;
  void   ReconstructEventPbPb(TTree *clustersTree, 
  	    Float_t* ZN1ADCCorr, Float_t* ZP1ADCCorr, Float_t* ZN2ADCCorr, Float_t* ZP2ADCCorr,
	    Float_t* ZEM1ADCCorr, Float_t* ZEM2ADCCorr, Float_t* PMRef1, Float_t* PMRef2,
	    Bool_t channelsOff, Bool_t chUnderflow, Bool_t chOverflow) const;
  void   BuildRecoParam(Float_t ZDCC, Float_t ZDCA, Float_t ZEM) const;
  
  void   FillZDCintoESD(TTree *clustersTree, AliESDEvent*esd) const;


  static AliZDCRecoParam *fRecoParam; // reconstruction parameters

  AliZDCPedestals  *fPedData; 	    //! pedestal calibration data
  AliZDCEnCalib    *fEnCalibData;   //! energy calibration data
  AliZDCTowerCalib *fTowCalibData;  //! equalization calibration data
  
  Int_t   fRecoMode;	    // =1->p-p, =2->A-A
  Float_t fBeamEnergy;	    // beam energy
  Int_t	  fNRun;	    // Run Number (from raw data)
  Bool_t  fIsCalibrationMB; // true if run type = "CALIBRATION_MB"
  Int_t   fPedSubMode;	    // =0->mean values, =1->from correlations
  Float_t fSignalThreshold; // Threshold value for "triggering" in p-p

  ClassDef(AliZDCReconstructor, 9)   // class for the ZDC reconstruction
};

#endif
