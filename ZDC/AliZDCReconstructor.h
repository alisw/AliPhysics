#ifndef ALIZDCRECONSTRUCTOR_H
#define ALIZDCRECONSTRUCTOR_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// class for ZDC reconstruction                                              //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "AliReconstructor.h"
#include "AliCDBManager.h"
#include "AliCDBStorage.h"
#include "AliZDCPedestals.h"
#include "AliZDCCalib.h"
#include "AliZDCRecoParam.h"
#include "AliZDCRecoParampp.h"
#include "AliZDCRecoParamPbPb.h"
#include "AliLog.h"

class AliLoader;

class AliZDCReconstructor: public AliReconstructor {
public:
  AliZDCReconstructor();
  virtual ~AliZDCReconstructor();

  virtual Bool_t HasDigitConversion() const {return kFALSE;};

  virtual void Reconstruct(TTree* digitsTree, TTree* clustersTree) const; 
  virtual void Reconstruct(AliRawReader* rawReader, TTree* clustersTree) const;

  virtual void FillESD(TTree* /*digitsTree*/, TTree* clustersTree, AliESDEvent* esd) const 
  	        {FillZDCintoESD(clustersTree,esd);}
  virtual void FillESD(AliRawReader* /*rawReader*/, TTree* clustersTree, AliESDEvent* esd) const 
  	        {FillZDCintoESD(clustersTree,esd);}

  // parameter settings for reconstruction
  void SetRecoMode();
  static void SetRecoParam(AliZDCRecoParam * param){fRecoParam = param;}
  
  Int_t   GetRecoMode() {return fRecoMode;}
  static const AliZDCRecoParam* GetRecoParam(){return fRecoParam;}
  
  // OCDB objects for reconstruction
  AliCDBStorage   *SetStorage(const char* uri);
  AliZDCPedestals *GetPedData() const; 
  AliZDCCalib     *GetECalibData() const; 
  
private:
  AliZDCReconstructor(const AliZDCReconstructor&);
  AliZDCReconstructor& operator =(const AliZDCReconstructor&);

  void   ReconstructEventpp(TTree *clustersTree, 
  	    Float_t* ZN1ADCCorr, Float_t* ZP1ADCCorr, Float_t* ZN2ADCCorr, Float_t* ZP2ADCCorr,
	    Float_t* ZEM1ADCCorr, Float_t* ZEM2ADCCorr, Float_t* PMRef1, Float_t* PMRef2) const;
  void   ReconstructEventPbPb(TTree *clustersTree, 
  	    Float_t* ZN1ADCCorr, Float_t* ZP1ADCCorr, Float_t* ZN2ADCCorr, Float_t* ZP2ADCCorr,
	    Float_t* ZEM1ADCCorr, Float_t* ZEM2ADCCorr, Float_t* PMRef1, Float_t* PMRef2) const;
  void   FillZDCintoESD(TTree *clustersTree, AliESDEvent*esd) const;

  static AliZDCRecoParam *fRecoParam; // reconstruction parameters

  AliZDCPedestals *fPedData; 	//! pedestal calibration data
  AliZDCCalib     *fECalibData; //! energy and equalization calibration data
  Int_t           fRecoMode;	// =0->p-p, =1->A-A
  Float_t         fBeamEnergy;	// beam energy

  ClassDef(AliZDCReconstructor, 5)   // class for the ZDC reconstruction
};

#endif
