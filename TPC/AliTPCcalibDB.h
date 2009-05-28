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

class AliTPCCalPad;
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
class TMap;
//class AliCDBStorage;

class AliTPCcalibDB : public TObject
{
 public: 
  static AliTPCcalibDB* Instance();
  AliTPCcalibDB();
  virtual ~AliTPCcalibDB();
  static void Terminate();
  void   SetRun(Long64_t run);   
  //
  AliTPCTransform* GetTransform() {return fTransform;}
  AliTPCExB*    GetExB() {return fExB;}
  void          SetExBField(Float_t bz);
  static AliTPCExB*    GetExB(Float_t bz,Bool_t bdelete);
  AliTPCCalPad* GetPadGainFactor() {return fPadGainFactor;}
  AliTPCCalPad* GetDedxGainFactor() {return fDedxGainFactor;}
  AliTPCCalPad* GetPadTime0() {return fPadTime0;}
  AliTPCCalPad* GetPadNoise() {return fPadNoise;}
  AliTPCCalPad* GetPedestals() {return fPedestals;}
  AliTPCSensorTempArray* GetTemperature() {return fTemperature;}
  AliTPCParam*  GetParameters(){return fParam;}
  AliTPCAltroMapping ** GetMapping(){ return fMapping;}
  AliTPCClusterParam *GetClusterParam(){ return fClusterParam;}
  TObjArray *GetTimeGainSplines(){ return fTimeGainSplines;}  
  //
  //
  static AliGRPObject * GetGRP(Int_t run);
  static TMap *  GetGRPMap(Int_t run);
  static Float_t GetPressure(Int_t timeStamp, Int_t run, Int_t type=0);
  static Float_t GetChamberHighVoltage(Int_t timeStamp, Int_t run, Int_t sector);
  static Float_t GetValueGoofie(Int_t timeStamp, Int_t run, Int_t type);
  static Bool_t  GetTemperatureFit(Int_t timeStamp, Int_t run, Int_t side,TVectorD& fit);
  static Float_t GetTemperature(Int_t timeStamp, Int_t run, Int_t side);
  static Double_t GetPTRelative(UInt_t timeSec, Int_t run,  Int_t side);
  AliDCSSensor * GetPressureSensor(Int_t run, Int_t type=0);
  //AliDCSSensor * GetVoltageSensor(Int_t run, Int_t type=0);
  AliTPCSensorTempArray * GetTemperatureSensor(Int_t run);
  AliDCSSensorArray *     GetGoofieSensors(Int_t run);
  AliDCSSensorArray *     GetVoltageSensors(Int_t run);
  AliTPCCalibVdrift *     GetVdrift(Int_t run);
  static Float_t GetGain(Int_t sector, Int_t row, Int_t pad);
  //
  static void     CreateObjectList(const Char_t *filename, TObjArray *calibObjects);
  static void MakeTree(const char * fileName, TObjArray * array, const char * mapFileName = 0, AliTPCCalPad* outlierPad = 0, Float_t ltmFraction = 0.9);
  static void RegisterExB(Int_t index, Float_t bz, Bool_t bdelete);
  //
  //
  static  void ProcessGoofie( AliDCSSensorArray* goofieArray, TVectorD & vecEntries, TVectorD & vecMedian, TVectorD &vecMean, TVectorD &vecRMS);
  static void ProcessEnv(const char * runList);

  AliGRPObject * MakeGRPObjectFromMap(TMap *map);
  //
protected:
  
  void         Update();  //update entries
  AliCDBEntry* GetCDBEntry(const char* cdbPath);   
  void GetRunInformations(Int_t run); // JUST FOR CALIBRATION STUDIES
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
  AliTPCCalPad* fPadNoise;        // Noise calibration entry
  AliTPCCalPad* fPedestals;       // Pedestal calibration entry
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
  TObjArray      fGRPArray;           //! array of GRPs  -  per run
  TObjArray      fGRPMaps;            //! array of GRPs maps  -  per run - old data  
  TObjArray      fGoofieArray;        //! array of GOOFIE values -per run
  TObjArray      fVoltageArray;       //! array of Chamber HV values -per run
  TObjArray      fTemperatureArray;   //! array of temperature sensors - per run
  TObjArray      fVdriftArray;        //! array of v drift interfaces
  TArrayI        fRunList;            //! run list - indicates try to get the run param
  //
  static AliTPCcalibDB* fgInstance;  // singleton control
  static Bool_t       fgTerminated;  // termination control 
  static TObjArray    fgExBArray;    // array of ExB corrections
  ClassDef(AliTPCcalibDB, 0)
 private:
   AliTPCcalibDB (const AliTPCcalibDB& );
   AliTPCcalibDB& operator= (const AliTPCcalibDB& );
};


#endif
