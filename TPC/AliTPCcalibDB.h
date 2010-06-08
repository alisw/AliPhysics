#ifndef ALITPCCALIBDB_H
#define ALITPCCALIBDB_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// Class providing the calibration parameters by accessing the CDB           //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////


class AliTPCTransform;
class AliTPCExB;
#include "TObject.h"
#include "TObjArray.h"
#include "TArrayI.h"
#include "TVectorD.h"
#include "TGraph.h"
#include "AliTPCCalPad.h"
#include "TString.h"
#include "AliSplineFit.h"

class AliTPCSensorTempArray;
class AliDCSSensorArray;
class AliCDBEntry;
class AliTPCParam;
class AliTPCAltroMapping;
class AliTPCClusterParam;
class AliDCSSensor;
class AliDCSSensorArray;
class AliTPCCalibVdrift;
class AliGRPObject;
class AliTPCCalibRaw;
class AliTPCdataQA;
class TMap;
class AliMagF;
class AliTPCcalibDButil;
class AliCTPTimeParams;
class AliTPCCorrection;
//class AliCDBStorage;

class AliTPCcalibDB : public TObject
{
 public: 
  static AliTPCcalibDB* Instance();
  AliTPCcalibDB();
  virtual ~AliTPCcalibDB();
  static void Terminate();
  void   SetRun(Long64_t run);   
  void   Update();  //update entries
  void   UpdateRunInformations(Int_t run, Bool_t force=kFALSE);
  void   UpdateNonRec();
  //
  Long64_t GetRun() const {return fRun;}
  //
  //
  //
  AliTPCTransform* GetTransform() const {return fTransform;}
  AliTPCExB*    GetExB() const {return fExB;}
  void          SetExBField(Float_t bz);
  void          SetExBField( const AliMagF*   bmap);
  static AliTPCExB*    GetExB(Float_t bz,Bool_t bdelete);
  AliTPCCalPad* GetPadGainFactor() const {return fPadGainFactor;}
  AliTPCCalPad* GetDedxGainFactor() const {return fDedxGainFactor;}
  AliTPCCalPad* GetPadTime0() const {return fPadTime0;}
  AliTPCCalPad* GetDistortionMap(Int_t i) const;
  AliTPCCorrection * GetTPCComposedCorrection() const { return fComposedCorrection;}
  void          SetTPCComposedCorrection(AliTPCCorrection *compCorr) { fComposedCorrection=compCorr;}
  AliTPCCorrection * GetTPCComposedCorrection(Float_t field) const;

  AliTPCCalPad* GetPadNoise() const {return fPadNoise;}
  AliTPCCalPad* GetPedestals() const{return fPedestals;}
  //ALTRO config data
  TObjArray* GetAltroConfigData()  const {return fALTROConfigData;}
  AliTPCCalPad* GetALTROAcqStart() const {return fALTROConfigData?static_cast<AliTPCCalPad*>(fALTROConfigData->FindObject("AcqStart")):0;}
  AliTPCCalPad* GetALTROZsThr()    const {return fALTROConfigData?static_cast<AliTPCCalPad*>(fALTROConfigData->FindObject("ZsThr")):0;}
  AliTPCCalPad* GetALTROFPED()     const {return fALTROConfigData?static_cast<AliTPCCalPad*>(fALTROConfigData->FindObject("FPED")):0;}
  AliTPCCalPad* GetALTROAcqStop()  const {return fALTROConfigData?static_cast<AliTPCCalPad*>(fALTROConfigData->FindObject("AcqStop")):0;}
  AliTPCCalPad* GetALTROMasked()   const {return fALTROConfigData?static_cast<AliTPCCalPad*>(fALTROConfigData->FindObject("Masked")):0;}
  TMap* GetRCUconfig() const {return fALTROConfigData?(TMap*)(fALTROConfigData->FindObject("RCUconfig")):0;}
  Int_t GetRCUTriggerConfig() const;
  Bool_t IsTrgL0();
  Bool_t IsTrgL1();
    
    
  //Pulser data
  TObjArray*    GetPulserData()  const {return fPulserData;}
  AliTPCCalPad* GetPulserTmean() const {return fPulserData?static_cast<AliTPCCalPad*>(fPulserData->FindObject("PulserTmean")):0;}
  AliTPCCalPad* GetPulserTrms()  const {return fPulserData?static_cast<AliTPCCalPad*>(fPulserData->FindObject("PulserTrms")):0;}
  AliTPCCalPad* GetPulserQmean() const {return fPulserData?static_cast<AliTPCCalPad*>(fPulserData->FindObject("PulserQmean")):0;}
  //CE data
  TObjArray*    GetCEData()     const {return fCEData;}
  AliTPCCalPad* GetCETmean()    const {return fCEData?static_cast<AliTPCCalPad*>(fCEData->FindObject("CETmean")):0;}
  AliTPCCalPad* GetCETrms()     const {return fCEData?static_cast<AliTPCCalPad*>(fCEData->FindObject("CETrms")):0;}
  AliTPCCalPad* GetCEQmean()    const {return fCEData?static_cast<AliTPCCalPad*>(fCEData->FindObject("CEQmean")):0;}
  TObjArray*    GetCErocTtime() const {return fCEData?static_cast<TObjArray*>(fCEData->FindObject("rocTtime")):0;}
  TObjArray*    GetCErocQtime() const {return fCEData?static_cast<TObjArray*>(fCEData->FindObject("rocQtime")):0;}
  TGraph*       GetCErocTgraph(const Int_t roc)const {return GetCErocTtime()?static_cast<TGraph*>(GetCErocTtime()->At(roc)):0;}
  TGraph*       GetCErocQgraph(const Int_t roc)const {return GetCErocQtime()?static_cast<TGraph*>(GetCErocQtime()->At(roc)):0;}
  static Float_t GetCEdriftTime(Int_t run, Int_t sector, Double_t timeStamp=-1., Int_t *entries=0);
  static Float_t GetCEchargeTime(Int_t run, Int_t sector, Double_t timeStamp=-1., Int_t *entries=0);
  //Raw calibration
  AliTPCCalibRaw* GetCalibRaw() const {return fCalibRaw;}
  //QA object
  AliTPCdataQA*   GetDataQA() const {return fDataQA;}
  //
  AliTPCSensorTempArray* GetTemperature() const {return fTemperature;}
  AliTPCParam*  GetParameters() const {return fParam;}
  AliTPCAltroMapping ** GetMapping() const{ return fMapping;}
  AliTPCClusterParam *GetClusterParam() const { return fClusterParam;}
  TObjArray *GetTimeGainSplines() const { return fTimeGainSplines;}  
  //
  //GRP information
  static AliGRPObject * GetGRP(Int_t run);
  static TMap *  GetGRPMap(Int_t run);
  static Float_t GetPressure(Int_t timeStamp, Int_t run, Int_t type=0);
  static Float_t GetL3Current(Int_t run, Int_t statType=0);
  static Float_t GetBz(Int_t run);
  static Char_t  GetL3Polarity(Int_t run);
  static TString GetRunType(Int_t run);
  //
  static Float_t GetDCSSensorValue(AliDCSSensorArray *arr, Int_t timeStamp, const char * sensorName, Int_t sigDigits=-1);
  static Float_t GetDCSSensorMeanValue(AliDCSSensorArray *arr, const char * sensorName, Int_t sigDigits=-1);
  //Voltage information
  static Float_t GetChamberHighVoltage(Int_t run, Int_t sector, Int_t timeStamp=-1, Int_t sigDigits=0);
  static Float_t GetSkirtVoltage(Int_t run, Int_t sector, Int_t timeStamp=-1, Int_t sigDigits=0);
  static Float_t GetCoverVoltage(Int_t run, Int_t sector, Int_t timeStamp=-1, Int_t sigDigits=0);
  static Float_t GetGGoffsetVoltage(Int_t run, Int_t sector, Int_t timeStamp=-1, Int_t sigDigits=0);
  static Float_t GetGGnegVoltage(Int_t run, Int_t sector, Int_t timeStamp=-1, Int_t sigDigits=0);
  static Float_t GetGGposVoltage(Int_t run, Int_t sector, Int_t timeStamp=-1, Int_t sigDigits=0);
  //Goofie Values
  static Float_t GetValueGoofie(Int_t timeStamp, Int_t run, Int_t type);
  //
  static Bool_t  GetTemperatureFit(Int_t timeStamp, Int_t run, Int_t side,TVectorD& fit);
  static Float_t GetTemperature(Int_t timeStamp, Int_t run, Int_t side);
  static Double_t GetPTRelative(UInt_t timeSec, Int_t run,  Int_t side);
  AliDCSSensor * GetPressureSensor(Int_t run, Int_t type=0);
  //AliDCSSensor * GetVoltageSensor(Int_t run, Int_t type=0);
  AliTPCSensorTempArray * GetTemperatureSensor(Int_t run);
  AliDCSSensorArray *     GetGoofieSensors(Int_t run);
  AliDCSSensorArray *     GetVoltageSensors(Int_t run);
  AliTPCCalibVdrift *     GetVdrift(Int_t run);
  TObjArray *							GetTimeGainSplinesRun(Int_t run);
  TObjArray*              GetTimeVdriftSplineRun(Int_t run);
  static Float_t GetGain(Int_t sector, Int_t row, Int_t pad);
  //
  Double_t      GetVDriftCorrectionTime(Int_t timeStamp, Int_t run, Int_t side, Int_t mode);
  Double_t      GetTime0CorrectionTime(Int_t timeStamp, Int_t run, Int_t side, Int_t mode);
  Double_t      GetVDriftCorrectionGy(Int_t timeStamp, Int_t run, Int_t side, Int_t mode);
  //
  AliSplineFit* GetVdriftSplineFit(const char* name, Int_t run);
  AliSplineFit* CreateVdriftSplineFit(const char* graphName, Int_t run);
  //
  static void     CreateObjectList(const Char_t *filename, TObjArray *calibObjects);
  static void MakeTree(const char * fileName, TObjArray * array, const char * mapFileName = 0, AliTPCCalPad* outlierPad = 0, Float_t ltmFraction = 0.9);
  static void RegisterExB(Int_t index, Float_t bz, Bool_t bdelete);
  //
  //
  
  AliTPCCalPad* MakeDeadMap(Double_t notInMap=1, const char *nameMappingFile="$ALICE_ROOT/TPC/Calib/tpcMapping.root" );
  AliGRPObject * MakeGRPObjectFromMap(TMap *map);
  AliCTPTimeParams* GetCTPTimeParams() const {return fCTPTimeParams;}
  //Create a tree suited for diplaying with the AliTPCCalibViewerGUI
  Bool_t CreateGUITree(const char* filename="");
  static Bool_t CreateGUITree(Int_t run, const char* filename="");
  static Bool_t CreateRefFile(Int_t run, const char* filename="");
  //
protected:
  
  AliCDBEntry* GetCDBEntry(const char* cdbPath);   
  Long64_t        fRun;         // current run number
  AliTPCTransform *fTransform;      // object responsible for spacial corrections
  AliTPCExB *fExB;              // ExB correction factor
//  AliCDBStorage* fLocator;      // Storage locator retrieved from AliCDBManager
  //
  // calibration parameters per pad
  //
  AliTPCCalPad* fPadGainFactor;   // Gain calibration entry
  AliTPCCalPad* fDedxGainFactor;   // Gain calibration entry - for dEdx
  AliTPCCalPad* fPadTime0;        // Time0 calibration entry
  TObjArray   *fDistortionMap;    // distortion map
  AliTPCCorrection *fComposedCorrection;  // general space point corrections
  TObjArray *      fComposedCorrectionArray; //space point corrections for different field setting
  AliTPCCalPad* fPadNoise;        // Noise calibration entry
  AliTPCCalPad* fPedestals;       // Pedestal calibration entry
  AliTPCCalibRaw *fCalibRaw;      // raw data calibration entry
  AliTPCdataQA  *fDataQA;         // qa object
  TObjArray *fALTROConfigData;    // ALTRO configuration data
  TObjArray *fPulserData;         // Calibration Pulser data
  TObjArray *fCEData;             // CE data
  //
  //
  //
  //
  AliTPCSensorTempArray* fTemperature; // Temperature calibration entry
  AliTPCAltroMapping **fMapping;   // Altro mapping   
  //
  //
  AliTPCParam * fParam;           // TPC parameters
  AliTPCClusterParam * fClusterParam;  // TPC cluster error, shape and Q parameterization
  TObjArray * fTimeGainSplines; // Array of AliSplineFits: at 0 MIP position in time ; at 1 Fermi Plateau from cosmics
  //
  // Get the corssrun information
  //
  TObjArray      fTimeGainSplinesArray; //! array Array of AliSplineFits: at 0 MIP position in time ; at 1 Fermi Plateau from cosmics
  TObjArray      fGRPArray;							//! array of GRPs  -  per run
  TObjArray      fGRPMaps;							//! array of GRPs maps  -  per run - old data  
  TObjArray      fGoofieArray;					//! array of GOOFIE values -per run
  TObjArray      fVoltageArray;					//! array of Chamber HV values -per run
  TObjArray      fTemperatureArray;			//! array of temperature sensors - per run
  TObjArray      fVdriftArray;					//! array of v drift interfaces
  TObjArray      fDriftCorrectionArray;                //! array of drift correction

  TArrayI        fRunList;							//! run list - indicates try to get the run param
  //
  static AliTPCcalibDB* fgInstance;  // singleton control
  static Bool_t       fgTerminated;  // termination control 
  static TObjArray    fgExBArray;    // array of ExB corrections
  AliTPCcalibDButil   *fDButil;       // utility class
  //ctp info
  AliCTPTimeParams *fCTPTimeParams;   //CTP timing parameters
  
  ClassDef(AliTPCcalibDB, 0)
 private:
   AliTPCcalibDB (const AliTPCcalibDB& );
   AliTPCcalibDB& operator= (const AliTPCcalibDB& );
};


#endif
