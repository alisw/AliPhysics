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
class AliZDCSaturationCalib;
class AliZDCTowerCalib;
class AliZDCMBCalib;
class AliZDCTDCCalib;
class AliZDCChMap;
class AliZDCRecoParampp;
class AliZDCRecoParamPbPb;
class AliLog;
class AliLoader;

class AliZDCReconstructor: public AliReconstructor {
public:
  //  MAPPING !!!! NB !!!!! MUST BE THE SAME AS FOR AliZDCRawStream.h !!!!
  // Module type codes
  enum ZDCModules{kV965=1, kV830=2, kTRG=3, kTRGI=4, kPU=5, KV1290=6, kV775N=7}; 
  
  // Module type codes
  enum ZDCGeoAddr{kFirstADCGeo=0, kLastADCGeo=3, kADDADCGeo=5,
       kTDCFakeGeo=8, kZDCTDCGeo=4, kADDTDCGeo=6,
       kScalerGeo=16, kPUGeo=29, kTrigScales=30, kTrigHistory=31};
  
  // Signal codes for ZDC 
  // Same codes used in DAQ configuration file
  // To be changed ONLY IF this file is changed!!! 
  // **** DO NOT CHANGE THE FOLLOWING LINES!!! ****
  enum ZDCSignal{
       kNotConnected=0, kVoid=1,
       kZNAC=2, kZNA1=3, kZNA2=4, kZNA3=5, kZNA4=6,
       kZPAC=7, kZPA1=8, kZPA2=9, kZPA3=10, kZPA4=11,
       kZNCC=12, kZNC1=13, kZNC2=14, kZNC3=15, kZNC4=16,
       kZPCC=17, kZPC1=18, kZPC2=19, kZPC3=20, kZPC4=21,
       kZEM1=22, kZEM2=23,
       kZDCAMon=24, kZDCCMon=25,
       kZNACoot=26, kZNA1oot=27, kZNA2oot=28, kZNA3oot=29, kZNA4oot=30,
       kZPACoot=31, kZPA1oot=32, kZPA2oot=33, kZPA3oot=34, kZPA4oot=35,
       kZNCCoot=36, kZNC1oot=37, kZNC2oot=38, kZNC3oot=39, kZNC4oot=40,
       kZPCCoot=41, kZPC1oot=42, kZPC2oot=43, kZPC3oot=44, kZPC4oot=45,
       kZEM1oot=46, kZEM2oot=47,
       kZDCAMonoot=48, kZDCCMonoot=49,
       kL1MBI=50, kL1CNI=51, kL1SCI=52, kL1EMDI=53, kL0I=54, 
       kL1MBO=55, kL1CNO=56, kL1SCO=57, kL1EMDO=58, 
       kHMBCN=59, kHSCEMD=60,
       kZNACD=61, kZNA1D=62, kZNA2D=63, kZNA3D=64, kZNA4D=65,
       kZPACD=66, kZPA1D=67, kZPA2D=68, kZPA3D=69, kZPA4D=70,
       kZNCCD=71, kZNC1D=72, kZNC2D=73, kZNC3D=74, kZNC4D=75,
       kZPCCD=76, kZPC1D=77, kZPC2D=78, kZPC3D=79, kZPC4D=80,
       kZEM1D=81, kZEM2D=82,
       kZDCAMonD=83, kZDCCMonD=84,
       kZNAD=85, kZPAD=86, kZNCD=87, kZPCD=88, kZEMD=89,
       kZNA0D=90, kZPA0D=91, kZNC0D=92, kZPC0D=93, k1kHzD=94, kGate=95, kAD=96, kCD=97, 
       kAorCD=98, kAandCD=99, kZEMORD=100, kAorCorZEMORD=101, kAorCorZEMD=102, kAD0=103, kAD1=104, kAD2=105, 
       kAD3=106, kAD4=107, kAD5=108, kAD6=109, kAD7=110, kAD8=111, kAD9=112, kAD10=113, 
       kAD11=114, kAD12=115, kAD13=116, kAD14=117, kAD15=118, kAD0D=119, kAD1D=120, kAD2D=121,
       kAD3D=122, kAD4D=123, kAD5D=124, kAD6D=125, kAD7D=126, kAD8D=127, kAD9D=128, kAD10D=129,
       kAD11D=130, kAD12D=131, kAD13D=132, kAD14D=133, kAD15D=134, kL0=135,
       k1ZAC=136, k1ZED=137, k1ZMD=138, k1ZMB=139
       };
  //  MAPPING !!!! NB !!!!! MUST BE THE SAME AS FOR AliZDCRawStream.h !!!!

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
  
  int    GetChannelSignal(int det, int quad, Bool_t intime=kTRUE) const;
  
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
    
  void    SetSignalThreshold(Float_t val) {fSignalThreshold=val;}
  Float_t GetSignalThreshold() const {return fSignalThreshold;}
  
  // OCDB objects for reconstruction
  AliCDBStorage       *SetStorage(const char* uri);
  AliZDCPedestals     *GetPedestalData() const; 
  AliZDCEnCalib       *GetEnergyCalibData() const; 
  AliZDCSaturationCalib *GetSaturationCalibData() const; 
  AliZDCTowerCalib    *GetTowerCalibData() const; 
  AliZDCMBCalib       *GetMBCalibData() const; 
  AliZDCTDCCalib      *GetTDCCalibData() const; 
  AliZDCChMap         *GetMapping() const; 
  
private:
  AliZDCReconstructor(const AliZDCReconstructor&); //Not implemented
  AliZDCReconstructor& operator =(const AliZDCReconstructor&); //Not implemented

  void   ReconstructEventpp(TTree *clustersTree, 
	 Float_t corrADC[24][2], Int_t signalCodeADC[24], Int_t signalCodeTDC[7], Bool_t isScalerOn, UInt_t* scaler, 
	 Int_t tdcData[32][4], const Int_t* const evQualityBlock, 
	 const Int_t* const triggerBlock, const Int_t* const chBlock, UInt_t puBits) const;
	 
  void   ReconstructEventPbPb(TTree *clustersTree, 
	 Float_t corrADC[24][2], Int_t signalCodeADC[24], Int_t signalCodeTDC[7], Bool_t isScalerOn, UInt_t* scaler, 
	 Int_t tdcData[32][4], const Int_t* const evQualityBlock, 
	 const Int_t* const triggerBlock, const Int_t* const chBlock, UInt_t puBits) const;

  static AliZDCRecoParam *fgRecoParam; // reconstruction parameters

  static AliZDCMBCalib   *fgMBCalibData;  //! mb calibration data
  AliZDCPedestals  *fPedData; 	    	  //! pedestal calibration data
  AliZDCEnCalib    *fEnCalibData;   	  //! energy calibration data
  AliZDCSaturationCalib  *fSatCalibData;  //! energy calibration data
  AliZDCTowerCalib *fTowCalibData;  	  //! equalization calibration data
  AliZDCTDCCalib   *fTDCCalibData;  	  //! TDC offset data
  AliZDCChMap      *fMapping;  	  	  //! Mapping from OCDB
  
  Int_t    fRecoMode;	     // =1->p-p, p-A (A-p), =2->A-A
  Float_t  fBeamEnergy;	     // beam energy
  Int_t	   fNRun;	     // Run Number (from raw data)
  Bool_t   fIsCalibrationMB; // true if run type = "CALIBRATION_MB"
  Float_t  fSignalThreshold; // Threshold value for "triggering" in p-p
  Double_t fMeanPhase;       // LHC clock phase
  
  AliESDZDC* fESDZDC;        // ESD output object  

  ClassDef(AliZDCReconstructor, 15)   // class for the ZDC reconstruction
};

#endif
