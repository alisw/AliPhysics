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
#include "AliZDCRecoParam.h"
#include "AliESDZDC.h"

class AliCDBManager;
class AliCDBStorage;
class AliZDCPedestals;
class AliZDCEnCalib;
class AliZDCTowerCalib;
class AliZDCMBCalib;
class AliZDCTDCCalib;
class AliZDCRecoParampp;
class AliZDCRecoParamPbPb;
class AliLog;
class AliLoader;

class AliZDCReconstructor: public AliReconstructor {
public:
  AliZDCReconstructor();
  virtual ~AliZDCReconstructor();

  virtual void   Init();
  virtual void   Init(TString beamType, Float_t beamEnergy);
  virtual Bool_t HasDigitConversion() const {return kFALSE;};
  
  virtual void Reconstruct(TTree*digitsTree, TTree* clustersTree) const; 
  virtual void Reconstruct(AliRawReader* rawReader, TTree* clustersTree) const;

  virtual void FillESD(TTree* /*digitsTree*/, TTree* clustersTree, 
  		       AliESDEvent* esd) const {FillZDCintoESD(clustersTree, esd);}
  virtual void FillESD(AliRawReader* /*rawReader*/, TTree* clustersTree, 
		       AliESDEvent* esd) const {FillZDCintoESD(clustersTree, esd);}
  
  void   FillZDCintoESD(TTree *clustersTree, AliESDEvent *esd) const;
  
  // parameter settings for reconstruction
  void SetRecoMode(Int_t recoMode, Float_t beamEnergy) 
       {fRecoMode=recoMode; fBeamEnergy=beamEnergy;}
  static void SetRecoParam(AliZDCRecoParam * const param)
  	      {fgRecoParam = param;}
  
  Int_t   GetRecoMode() const {return fRecoMode;}
  Float_t GetBeamEnergy() const {return fBeamEnergy;}
  
  AliESDZDC* GetZDCESDData() const {return fESDZDC;}
  
  static const AliZDCRecoParam* GetRecoParam() 
  {return dynamic_cast<const AliZDCRecoParam*>(AliReconstructor::GetRecoParam(9));}
    
  void  SetPedSubMode(Int_t pedsubMode) {fPedSubMode=pedsubMode;}
  Int_t GetPedSubMode() const {return fPedSubMode;}
  
  void    SetSignalThreshold(Float_t val) {fSignalThreshold=val;}
  Float_t GetSignalThreshold() const {return fSignalThreshold;}
  
  // OCDB objects for reconstruction
  AliCDBStorage       *SetStorage(const char* uri);
  AliZDCPedestals     *GetPedestalData() const; 
  AliZDCEnCalib       *GetEnergyCalibData() const; 
  AliZDCTowerCalib    *GetTowerCalibData() const; 
  AliZDCMBCalib       *GetMBCalibData() const; 
  AliZDCTDCCalib      *GetTDCCalibData() const; 
  
private:
  AliZDCReconstructor(const AliZDCReconstructor&); //Not implemented
  AliZDCReconstructor& operator =(const AliZDCReconstructor&); //Not implemented

  void   ReconstructEventpp(TTree *clustersTree, 
	 const Float_t* const corrADCZN1, const Float_t* const corrADCZP1, 
	 const Float_t* const corrADCZN2, const Float_t* const corrADCZP2,
	 const Float_t* const corrADCZEM1, const Float_t* const corrADCZEM2,
	 Float_t* sPMRef1, Float_t* sPMRef2, Bool_t isScalerOn, UInt_t* scaler, 
	 Int_t tdcData[32][4], const Int_t* const evQualityBlock, 
	 const Int_t* const triggerBlock, const Int_t* const chBlock, UInt_t puBits) const;
  void   ReconstructEventPbPb(TTree *clustersTree, 
	 const Float_t* const corrADCZN1, const Float_t* const corrADCZP1, 
	 const Float_t* const corrADCZN2, const Float_t* const corrADCZP2,
	 const Float_t* const corrADCZEM1, const Float_t* const corrADCZEM2,
	 Float_t* sPMRef1, Float_t* sPMRef2, Bool_t isScalerOn, UInt_t* scaler, 
	 Int_t tdcData[32][4], const Int_t* const evQualityBlock, 
	 const Int_t* const triggerBlock, const Int_t* const chBlock, UInt_t puBits) const;

  static AliZDCRecoParam *fgRecoParam; // reconstruction parameters

  static AliZDCMBCalib   *fgMBCalibData;   //! mb calibration data
  AliZDCPedestals  *fPedData; 	    	  //! pedestal calibration data
  AliZDCEnCalib    *fEnCalibData;   	  //! energy calibration data
  AliZDCTowerCalib *fTowCalibData;  	  //! equalization calibration data
  AliZDCTDCCalib   *fTDCCalibData;  	  //! TDC offset data
  
  Int_t    fRecoMode;	    // =1->p-p, =2->A-A
  Float_t  fBeamEnergy;	    // beam energy
  Int_t	   fNRun;	    // Run Number (from raw data)
  Bool_t   fIsCalibrationMB; // true if run type = "CALIBRATION_MB"
  Int_t    fPedSubMode;	    // =0->mean values, =1->from correlations
  Float_t  fSignalThreshold; // Threshold value for "triggering" in p-p
  Double_t fMeanPhase;      // LHC clock phase
  
  AliESDZDC* fESDZDC;       // ESD output object  

  ClassDef(AliZDCReconstructor, 13)   // class for the ZDC reconstruction
};

#endif
