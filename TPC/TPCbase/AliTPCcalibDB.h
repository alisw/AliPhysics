#ifndef ALITPCCALIBDB_H
#define ALITPCCALIBDB_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/// \class AliTPCcalibDB
/// \brief Class providing the calibration parameters by accessing the CDB


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
#include "AliRecoParam.h"
#include "TMap.h"

class TGraphErrors;
class AliTPCSensorTempArray;
class AliDCSSensorArray;
class AliCDBEntry;
class AliTPCParam;
class AliTPCAltroMapping;
class AliTPCClusterParam;
class AliTPCRecoParam;
class AliDCSSensor;
class AliDCSSensorArray;
class AliTPCCalibVdrift;
class AliGRPObject;
class AliTPCCalibRaw;
class AliTPCdataQA;
class AliMagF;
class AliTPCcalibDButil;
class AliCTPTimeParams;
class AliTPCCorrection;
class AliTPCChebCorr;
//class AliCDBStorage;

class AliTPCcalibDB : public TObject
{
 public:
  enum EDcsGasSensor { kNeon=0, kArgon, kCO2, kN2, kH2O, kO2, kNGasSensor };

  static AliTPCcalibDB* Instance();
  AliTPCcalibDB();
  virtual ~AliTPCcalibDB();
  static void Terminate();
  void   SetRun(Long64_t run);   
  void   Update();  //update entries
  void   UpdateRunInformations(Int_t run, Bool_t force=kFALSE);
  void   UpdateNonRec();
  Bool_t   GetTailcancelationGraphs(Int_t sector, TGraphErrors ** graphRes, Float_t * indexAmpGraphs);
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
  AliTPCCalPad* GetPadGainFactorOld() const {return fPadGainFactor;}
  AliTPCCalPad* GetPadGainFactor() const {return fActiveChannelMap;}
  AliTPCCalPad* GetActiveChannelMap() const { return fActiveChannelMap; }
  AliTPCCalPad* GetDedxGainFactor() const {return fDedxGainFactor;}
  AliTPCCalPad* GetPadTime0() const {return fPadTime0;}
  AliTPCCalPad* GetDistortionMap(Int_t i) const;
  AliTPCCorrection * GetTPCComposedCorrection() const { return fComposedCorrection;}
  TObjArray * GetTPCComposedCorrectionArray() const { return fComposedCorrectionArray;}
  TObjArray*  GetCorrectionMaps()             const {return fCorrectionMaps;}
  void          SetTPCComposedCorrection(AliTPCCorrection *compCorr) { fComposedCorrection=compCorr;}
  AliTPCCorrection * GetTPCComposedCorrection(Float_t field) const;
  AliTPCCorrection * GetTPCComposedCorrectionDelta() const;
  Bool_t      HasAlignmentOCDB() const { return fBHasAlignmentOCDB;}

  AliTPCCalPad* GetPadNoise() const {return fPadNoise;}
  AliTPCCalPad* GetPedestals() const{return fPedestals;}

  void LoadCorrectionMaps();

  // ===| ALTRO config data |===================================================
  TObjArray* GetAltroConfigData()  const {return fALTROConfigData;}
  AliTPCCalPad* GetALTROAcqStart() const {return fALTROConfigData?static_cast<AliTPCCalPad*>(fALTROConfigData->FindObject("AcqStart")):0;}
  AliTPCCalPad* GetALTROZsThr()    const {return fALTROConfigData?static_cast<AliTPCCalPad*>(fALTROConfigData->FindObject("ZsThr")):0;}
  AliTPCCalPad* GetALTROFPED()     const {return fALTROConfigData?static_cast<AliTPCCalPad*>(fALTROConfigData->FindObject("FPED")):0;}
  AliTPCCalPad* GetALTROAcqStop()  const {return fALTROConfigData?static_cast<AliTPCCalPad*>(fALTROConfigData->FindObject("AcqStop")):0;}
  AliTPCCalPad* GetALTROMasked()   const {return fALTROConfigData?static_cast<AliTPCCalPad*>(fALTROConfigData->FindObject("Masked")):0;}
  TMap* GetDDLMap()   const {return fALTROConfigData?static_cast<TMap*>(fALTROConfigData->FindObject("DDLArray")):0;}
  TMap* GetRCUconfig() const {return fALTROConfigData?(TMap*)(fALTROConfigData->FindObject("RCUconfig")):0;}
  Int_t GetRCUTriggerConfig() const;
  Bool_t IsTrgL0();
  Bool_t IsTrgL1();
  Int_t GetMaxTimeBinAllPads() const { return fMaxTimeBinAllPads; }

  //
  TObjArray*    GetIonTailArray()  const {return fIonTailArray;}
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
  TObjArray*    GetCEfitsDrift()const {return fCEData?static_cast<TObjArray*>(fCEData->FindObject("ceFitsDrift")):0;}
  TGraph*       GetCErocTgraph(const Int_t roc)const {return GetCErocTtime()?static_cast<TGraph*>(GetCErocTtime()->At(roc)):0;}
  TGraph*       GetCErocQgraph(const Int_t roc)const {return GetCErocQtime()?static_cast<TGraph*>(GetCErocQtime()->At(roc)):0;}
  static Float_t GetCEdriftTime(Int_t run, Int_t sector, Double_t timeStamp=-1., Int_t *entries=0);
  static Float_t GetCEchargeTime(Int_t run, Int_t sector, Double_t timeStamp=-1., Int_t *entries=0);
  //Raw calibration
  AliTPCCalibRaw* GetCalibRaw() const {return fCalibRaw;}

  //QA object
  AliTPCdataQA*   GetDataQA() const {return fDataQA;}
  //
  // Gas sensor data
  //
  AliDCSSensorArray* GetGasSensors() const { return fGasSensorArray; }
  Float_t GetGasSensorValue(EDcsGasSensor type, Int_t timeStamp=-1, Int_t sigDigits=-1);

  //
  AliTPCSensorTempArray* GetTemperature() const {return fTemperature;}
  AliTPCParam*  GetParameters() const {return fParam;}
  AliTPCAltroMapping ** GetMapping() const{ return fMapping;}
  AliTPCClusterParam *GetClusterParam() const { return fClusterParam;}
  TObjArray *GetTimeGainSplines() const { return fTimeGainSplines;}  
  //
  //Reco param getter
  AliTPCRecoParam *GetRecoParam(Int_t i) const;
  AliTPCRecoParam *GetRecoParamFromGRP() const;
  AliRecoParam::EventSpecie_t GetEventSpecieFromGRP() const {return fRunEventSpecie;}
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
  static Float_t GetChamberHighVoltage(Int_t run, Int_t sector, Int_t timeStamp=-1, Int_t sigDigits=0, Bool_t current=kFALSE);
  static Float_t GetSkirtVoltage(Int_t run, Int_t sector, Int_t timeStamp=-1, Int_t sigDigits=0);
  static Float_t GetCoverVoltage(Int_t run, Int_t sector, Int_t timeStamp=-1, Int_t sigDigits=0);
  static Float_t GetGGoffsetVoltage(Int_t run, Int_t sector, Int_t timeStamp=-1, Int_t sigDigits=0);
  static Float_t GetGGnegVoltage(Int_t run, Int_t sector, Int_t timeStamp=-1, Int_t sigDigits=0);
  static Float_t GetGGposVoltage(Int_t run, Int_t sector, Int_t timeStamp=-1, Int_t sigDigits=0);
  //
  Bool_t  GetChamberHVStatus(UInt_t roc)                  const { return (roc<72)?fChamberHVStatus[roc]:kFALSE;   }
  Float_t GetChamberHighVoltageMedian(UInt_t roc)         const { return (roc<72)?fChamberHVmedian[roc]:0.;       }
  Float_t GetChamberCurrentNominalHighVoltage(UInt_t roc) const { return (roc<72)?fCurrentNominalVoltage[roc]:0.; }
  Float_t GetChamberGoodHighVoltageFraction(UInt_t roc)   const { return (roc<72)?fChamberHVgoodFraction[roc]:0.; }
  AliDCSSensor* GetChamberHVSensor(UInt_t roc)            const { return (roc<72)?fHVsensors[roc]:0x0;            }
  Double_t GetGainCorrectionHVandPT(Int_t timeStamp, Int_t run, Int_t sector, Int_t deltaCache, Int_t mode);
  Bool_t  IsDataTakingActive(time_t timeStamp);
  //
  //Goofie Values
  //
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
  TObjArray *             GetTimeGainSplinesRun(Int_t run);
  TObjArray*              GetTimeVdriftSplineRun(Int_t run);
  static Float_t GetGain(Int_t sector, Int_t row, Int_t pad);
  //
  // Drift velocity information
  //
  Double_t      GetVDriftCorrectionTime(Int_t timeStamp, Int_t run, Int_t side, Int_t mode);
  Double_t      GetTime0CorrectionTime(Int_t timeStamp, Int_t run, Int_t side, Int_t mode);
  Double_t      GetVDriftCorrectionGy(Int_t timeStamp, Int_t run, Int_t side, Int_t mode);
  Double_t      GetVDriftCorrectionDeltaZ(Int_t timeStamp, Int_t run, Int_t side, Int_t mode);
  //
  AliSplineFit* GetVdriftSplineFit(const char* name, Int_t run);
  AliSplineFit* CreateVdriftSplineFit(const char* graphName, Int_t run);
  //
  static void     CreateObjectList(const Char_t *filename, TObjArray *calibObjects);
  static void MakeTree(const char * fileName, TObjArray * array, const char * mapFileName = 0, AliTPCCalPad* outlierPad = 0, Float_t ltmFraction = 0.9);
  static void RegisterExB(Int_t index, Float_t bz, Bool_t bdelete);
  //
  // Dead channel map functions
  //
  Int_t GetMaskedChannelsFromCorrectionMaps(TBits maskedPads[72]);
  //
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
  void         UpdateChamberHighVoltageData();
  Int_t        InitDeadMap();
  void         InitAltroData();

  Int_t        fRun;         ///< current run number
  AliTPCTransform *fTransform;      ///< object responsible for spacial corrections
  AliTPCExB *fExB;              ///< ExB correction factor
//  AliCDBStorage* fLocator;      // Storage locator retrieved from AliCDBManager
  //
  // calibration parameters per pad
  //
  AliTPCCalPad* fPadGainFactor;   ///< Gain calibration entry
  AliTPCCalPad* fActiveChannelMap; ///< Map of active channels calculated on the fly
  AliTPCCalPad* fDedxGainFactor;   ///< Gain calibration entry - for dEdx
  AliTPCCalPad* fPadTime0;        ///< Time0 calibration entry
  TObjArray   *fDistortionMap;    ///< distortion map
  AliTPCCorrection *fComposedCorrection;  ///< general space point corrections
  TObjArray *      fComposedCorrectionArray; ///< space point corrections for different field setting
  TObjArray*       fCorrectionMaps;          ///< RS: new fast Chebyshev parameterization maps
  AliTPCCalPad* fPadNoise;        ///< Noise calibration entry
  AliTPCCalPad* fPedestals;       ///< Pedestal calibration entry
  AliTPCCalibRaw *fCalibRaw;      ///< raw data calibration entry
  AliTPCdataQA  *fDataQA;         ///< qa object
  TObjArray *fALTROConfigData;    ///< ALTRO configuration data
  TObjArray * fIonTailArray;      ///< array of graphs with the ion tail
  TObjArray *fPulserData;         ///< Calibration Pulser data
  TObjArray *fCEData;             ///< CE data
  //
  // Gas sensor data
  //
  static const char* fgkGasSensorNames[kNGasSensor]; ///< DCS gas sensor names
  AliDCSSensorArray* fGasSensorArray;                ///< Gas sensors

  //
  // Defived ALTRO information
  //
  Int_t fMaxTimeBinAllPads;       ///< Maximum Time bin in whole TPC extracted from AltroConfig
  //
  // Chamber HV info
  //
  Bool_t  fChamberHVStatus[72];       ///< Status of the Chamber, HV wise (on/off)
  Float_t fChamberHVmedian[72];       ///< median chamber high voltage
  Float_t fCurrentNominalVoltage[72]; ///< current nominal voltages
  Float_t fChamberHVgoodFraction[72]; ///< fraction of time the chamber has a good HV (wrt. robust median)
  AliDCSSensor *fHVsensors[72];       ///< HV sensors
  TGraph *fGrRunState;                ///< store information if run is active or paused
  //
  //
  //
  AliTPCSensorTempArray* fTemperature; ///< Temperature calibration entry
  AliTPCAltroMapping **fMapping;   ///< Altro mapping
  //
  //
  AliRecoParam::EventSpecie_t fRunEventSpecie;    ///< Event specie suggested for the run according to GRP
  AliTPCParam * fParam;                ///< TPC parameters
  AliTPCClusterParam * fClusterParam;  ///< TPC cluster error, shape and Q parameterization
  TObjArray * fRecoParamList;          ///< List of TPC reco param objects
  TObjArray * fTimeGainSplines;        ///< Array of AliSplineFits: at 0 MIP position in time ; at 1 Fermi Plateau from cosmics
  //
  // Get the corssrun information
  //
  TMap      fTimeGainSplinesArray; //!<! array Array of AliSplineFits: at 0 MIP position in time ; at 1 Fermi Plateau from cosmics
  TMap      fGRPArray;							//!<! array of GRPs  -  per run
  TMap      fGRPMaps;							//!<! array of GRPs maps  -  per run - old data
  TMap      fGoofieArray;					//!<! array of GOOFIE values -per run
  TMap      fVoltageArray;					//!<! array of Chamber HV values -per run
  TMap      fTemperatureArray;			//!<! array of temperature sensors - per run
  TMap      fVdriftArray;					//!<! array of v drift interfaces
  TMap      fDriftCorrectionArray;                //!<! array of drift correction

  TArrayI        fRunList;							//!<! run list - indicates try to get the run param
  Bool_t         fBHasAlignmentOCDB;                ///< Flag - alignment from the Transformation class
  //
  static AliTPCcalibDB* fgInstance;  ///< singleton control
  static Bool_t       fgTerminated;  ///< termination control
  static TObjArray    fgExBArray;    ///< array of ExB corrections
  AliTPCcalibDButil   *fDButil;       ///< utility class
  //ctp info
  AliCTPTimeParams *fCTPTimeParams;   ///< CTP timing parameters
  Int_t            fMode;             ///< RCU trigger config mode

 private:
   AliTPCcalibDB (const AliTPCcalibDB& );
   AliTPCcalibDB& operator= (const AliTPCcalibDB& );
  
   /// \cond CLASSIMP
   ClassDef(AliTPCcalibDB, 3)
   /// \endcond
};

#endif
