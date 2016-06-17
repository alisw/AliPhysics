/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/// \class AliTPCcalibDB
/// \brief Class providing the calibration parameters by accessing the CDB
///
/// Request an instance with AliTPCcalibDB::Instance()
/// If a new event is processed set the event number with SetRun
/// Then request the calibration data
///
/// Calibration data:
/// 0.)  Altro mapping
///          Simulation      - not yet
///          Reconstruction  - AliTPCclusterer::Digits2Clusters(AliRawReader* rawReader)
///
/// 1.)  pad by pad calibration -  AliTPCCalPad
///
///      a.) fPadGainFactor
///          Simulation: AliTPCDigitizer::ExecFast - Multiply by gain
///          Reconstruction : AliTPCclusterer::Digits2Clusters - Divide by gain
///
///      b.) fPadNoise -
///          Simulation:        AliTPCDigitizer::ExecFast
///          Reconstruction:    AliTPCclusterer::FindClusters(AliTPCCalROC * noiseROC)
///                             Noise depending cut on clusters charge (n sigma)
///      c.) fPedestal:
///          Simulation:     Not used yet - To be impleneted - Rounding to the nearest integer
///          Reconstruction: Used in AliTPCclusterer::Digits2Clusters(AliRawReader* rawReader)
///                          if data taken without zero suppression
///                          Currently switch in  fRecoParam->GetCalcPedestal();
///
///      d.) fPadTime0
///          Simulation:      applied in the AliTPC::MakeSector - adding offset
///          Reconstruction:  AliTPCTransform::Transform() - remove offset
///                           AliTPCTransform::Transform() - to be called
///                           in AliTPCtracker::Transform()
///
/// 2.)  Space points transformation:
///
///      a.) General coordinate tranformation - AliTPCtransform (see $ALICE_ROOT/TPC/AliTPCtransform.cxx)
///          Created on fly - use the other calibration components
///                 Unisochronity  - (substract time0 - pad by pad)
///                 Drift velocity - Currently common drift velocity - functionality of AliTPCParam
///                 ExB effect
///          Simulation     - Not used directly (the effects are applied one by one (see AliTPC::MakeSector)
///          Reconstruction -
///                           AliTPCclusterer::AddCluster
///                           AliTPCtracker::Transform
///      b.) ExB effect calibration -
///             classes (base class AliTPCExB, implementation- AliTPCExBExact.h  AliTPCExBFirst.h)
///             a.a) Simulation:   applied in the AliTPC::MakeSector -
///                                calib->GetExB()->CorrectInverse(dxyz0,dxyz1);
///             a.b) Reconstruction -
///
///                  in AliTPCtransform::Correct() - called calib->GetExB()->Correct(dxyz0,dxyz1)
///
///  3.)   cluster error, shape and Q parameterization

#include <iostream>
#include <fstream>


#include <AliCDBManager.h>
#include <AliCDBEntry.h>
#include <AliCDBId.h>
#include <AliLog.h>
#include <AliMagF.h>
#include <AliSplineFit.h>
#include <AliCTPTimeParams.h>

#include "TGraphErrors.h"
#include "AliTPCcalibDB.h"
#include "AliTPCdataQA.h"
#include "AliTPCcalibDButil.h"
#include "AliTPCAltroMapping.h"
#include "AliTPCExB.h"

#include "AliTPCCalROC.h"
#include "AliTPCCalPad.h"
#include "AliTPCSensorTempArray.h"
#include "AliGRPObject.h"
#include "AliTPCTransform.h"
#include "AliTPCmapper.h"
#include "AliTPCclusterMI.h"

class AliCDBStorage;
class AliTPCCalDet;
//
//

#include "TFile.h"
#include "TKey.h"
#include "TGraphErrors.h"
#include "TGeoGlobalMagField.h"

#include "TObjArray.h"
#include "TObjString.h"
#include "TString.h"
#include "TDirectory.h"
#include "TArrayI.h"
#include "AliTPCCalPad.h"
#include "AliTPCCalibPulser.h"
#include "AliTPCCalibPedestal.h"
#include "AliTPCCalibCE.h"
#include "AliTPCExBFirst.h"
#include "AliTPCTempMap.h"
#include "AliTPCCalibVdrift.h"
#include "AliTPCCalibRaw.h"
#include "AliTPCParam.h"
#include "AliTPCCorrection.h"
#include "AliTPCComposedCorrection.h"
#include "AliTPCPreprocessorOnline.h"
#include "AliTimeStamp.h"
#include "AliTriggerRunScalers.h"
#include "AliTriggerScalers.h"
#include "AliTriggerScalersRecord.h"
#include "AliDAQ.h"
#include "AliTPCRecoParam.h"
/// \cond CLASSIMP
ClassImp(AliTPCcalibDB)
/// \endcond

AliTPCcalibDB* AliTPCcalibDB::fgInstance = 0;
Bool_t AliTPCcalibDB::fgTerminated = kFALSE;
TObjArray    AliTPCcalibDB::fgExBArray;    ///< array of ExB corrections


//_ singleton implementation __________________________________________________
AliTPCcalibDB* AliTPCcalibDB::Instance()
{
  /// Singleton implementation
  /// Returns an instance of this class, it is created if necessary

  if (fgTerminated != kFALSE)
    return 0;

  if (fgInstance == 0)
    fgInstance = new AliTPCcalibDB();

  return fgInstance;
}

void AliTPCcalibDB::Terminate()
{
  /// Singleton implementation
  /// Deletes the instance of this class and sets the terminated flag, instances cannot be requested anymore
  /// This function can be called several times.

  fgTerminated = kTRUE;

  if (fgInstance != 0)
  {
    delete fgInstance;
    fgInstance = 0;
  }
}

//_____________________________________________________________________________
AliTPCcalibDB::AliTPCcalibDB():
  TObject(),
  fRun(-1),
  fTransform(0),
  fExB(0),
  fPadGainFactor(0),
  fActiveChannelMap(0),
  fDedxGainFactor(0),
  fPadTime0(0),
  fDistortionMap(0),
  fComposedCorrection(0),
  fComposedCorrectionArray(0),
  fCorrectionMaps(0),
  fPadNoise(0),
  fPedestals(0),
  fCalibRaw(0),
  fDataQA(0),
  fALTROConfigData(0),
  fIonTailArray(0),
  fPulserData(0),
  fCEData(0),
  fMaxTimeBinAllPads(-1),
  fHVsensors(),
  fGrRunState(0x0),
  fTemperature(0),
  fMapping(0),
  fParam(0),
  fClusterParam(0),
  fRecoParamList(0),
  fTimeGainSplines(0),
  fTimeGainSplinesArray(1),
  fGRPArray(1),            //! array of GRPs  -  per run  - JUST for calibration studies
  fGRPMaps(1),            //! array of GRPs  -  per run  - JUST for calibration studies
  fGoofieArray(1),         //! array of GOOFIE values -per run - Just for calibration studies
  fVoltageArray(1),
  fTemperatureArray(1),    //! array of temperature sensors - per run - Just for calibration studies
  fVdriftArray(1),                 //! array of v drift interfaces
  fDriftCorrectionArray(1),  //! array of drift correction
  fRunList(1),              //! run list - indicates try to get the run param
  fBHasAlignmentOCDB(kFALSE),    // Flag  - has the alignment on the composed correction ?
  fDButil(0),
  fCTPTimeParams(0),
  fMode(-1)
{
  /// constructor

  fgInstance=this;
  for (Int_t i=0;i<72;++i){
    fChamberHVStatus[i]=kTRUE;
    fChamberHVmedian[i]=-1;
    fCurrentNominalVoltage[i]=0.;
    fChamberHVgoodFraction[i]=0.;
  }
  Update();    // temporary
  fTimeGainSplinesArray.SetOwner(); //own the keys
  fGRPArray.SetOwner(); //own the keys
  fGRPMaps.SetOwner(); //own the keys
  fGoofieArray.SetOwner(); //own the keys
  fVoltageArray.SetOwner(); //own the keys
  fTemperatureArray.SetOwner(); //own the keys
  fVdriftArray.SetOwner(); //own the keys
  fDriftCorrectionArray.SetOwner(); //own the keys
}

AliTPCcalibDB::AliTPCcalibDB(const AliTPCcalibDB& ):
  TObject(),
  fRun(-1),
  fTransform(0),
  fExB(0),
  fPadGainFactor(0),
  fActiveChannelMap(0),
  fDedxGainFactor(0),
  fPadTime0(0),
  fDistortionMap(0),
  fComposedCorrection(0),
  fComposedCorrectionArray(0),
  fCorrectionMaps(0),
  fPadNoise(0),
  fPedestals(0),
  fCalibRaw(0),
  fDataQA(0),
  fALTROConfigData(0),
  fIonTailArray(0),
  fPulserData(0),
  fCEData(0),
  fMaxTimeBinAllPads(-1),
  fHVsensors(),
  fGrRunState(0x0),
  fTemperature(0),
  fMapping(0),
  fParam(0),
  fClusterParam(0),
  fRecoParamList(0),
  fTimeGainSplines(0),
  fTimeGainSplinesArray(1),
  fGRPArray(0),          //! array of GRPs  -  per run  - JUST for calibration studies
  fGRPMaps(0),          //! array of GRPs  -  per run  - JUST for calibration studies
  fGoofieArray(0),        //! array of GOOFIE values -per run - Just for calibration studies
  fVoltageArray(0),
  fTemperatureArray(0),   //! array of temperature sensors - per run - Just for calibration studies
  fVdriftArray(0),         //! array of v drift interfaces
  fDriftCorrectionArray(0), //! array of v drift corrections
  fRunList(0),              //! run list - indicates try to get the run param
  fBHasAlignmentOCDB(kFALSE),    // Flag  - has the alignment on the composed correction ?
  fDButil(0),
  fCTPTimeParams(0),
  fMode(-1)
{
  /// Copy constructor invalid -- singleton implementation

  Error("copy constructor","invalid -- singleton implementation");
  for (Int_t i=0;i<72;++i){
    fChamberHVStatus[i]=kTRUE;
    fChamberHVmedian[i]=-1;
    fCurrentNominalVoltage[i]=0.;
    fChamberHVgoodFraction[i]=0.;
  }
  fTimeGainSplinesArray.SetOwner(); //own the keys
  fGRPArray.SetOwner(); //own the keys
  fGRPMaps.SetOwner(); //own the keys
  fGoofieArray.SetOwner(); //own the keys
  fVoltageArray.SetOwner(); //own the keys
  fTemperatureArray.SetOwner(); //own the keys
  fVdriftArray.SetOwner(); //own the keys
  fDriftCorrectionArray.SetOwner(); //own the keys
}

AliTPCcalibDB& AliTPCcalibDB::operator= (const AliTPCcalibDB& )
{
/// Singleton implementation - no assignment operator

  Error("operator =", "assignment operator not implemented");
  return *this;
}



//_____________________________________________________________________________
AliTPCcalibDB::~AliTPCcalibDB()
{
  /// destructor
  ///
  /// delete fIonTailArray;

  delete fActiveChannelMap;
  delete fGrRunState;
  fgInstance = 0;
}

AliTPCCalPad* AliTPCcalibDB::GetDistortionMap(Int_t i) const {
  /// get distortion map - due E field distortions

  return (fDistortionMap) ? (AliTPCCalPad*)fDistortionMap->At(i):0;
}

AliTPCRecoParam* AliTPCcalibDB::GetRecoParam(Int_t i) const {
  return (fRecoParamList) ? (AliTPCRecoParam*)fRecoParamList->At(i):0;
}

//_____________________________________________________________________________
AliCDBEntry* AliTPCcalibDB::GetCDBEntry(const char* cdbPath)
{
  /// Retrieves an entry with path <cdbPath> from the CDB.

  char chinfo[1000];

  AliCDBEntry* entry = AliCDBManager::Instance()->Get(cdbPath, fRun);
  if (!entry)
  {
    snprintf(chinfo,1000,"AliTPCcalibDB: Failed to get entry:\t%s ", cdbPath);
    AliError(chinfo);
    return 0;
  }
  return entry;
}


//_____________________________________________________________________________
void AliTPCcalibDB::SetRun(Long64_t run)
{
  /// Sets current run number. Calibration data is read from the corresponding file.

  if (fRun == run)
    return;
	fRun = run;
  Update();
}



void AliTPCcalibDB::Update(){
  /// cache the OCDB entries for simulation, reconstruction, calibration


  AliCDBEntry * entry=0x0;
  Bool_t cdbCache = AliCDBManager::Instance()->GetCacheFlag(); // save cache status
  AliCDBManager::Instance()->SetCacheFlag(kTRUE); // activate CDB cache
  fDButil = new AliTPCcalibDButil;
  //
  fRun = AliCDBManager::Instance()->GetRun();

  //RS: new TPC correction map (array of AliTPCChebCorr) will be loaded on demand
  fCorrectionMaps = 0; // assuming that the object is managed by the OCDB
  //
  
  entry          = GetCDBEntry("TPC/Calib/Parameters");
  if (entry){
    //if (fPadNoise) delete fPadNoise;
    entry->SetOwner(kTRUE);
    fParam = (AliTPCParam*)(entry->GetObject());
  }else{
    AliFatal("TPC - Missing calibration entry TPC/Calib/Parameters");
  }

  //
  // check the presence of the detectors
  try {
    entry = AliCDBManager::Instance()->Get("GRP/GRP/Data");
  } catch(...) {
    AliInfo("No GRP entry found");
    entry = 0x0;
  }
  if (entry) {
    AliGRPObject* grpData = dynamic_cast<AliGRPObject*>(entry->GetObject());
    if (!grpData) {printf("Failed to get GRP data for run %d\n",fRun); return;}
    Int_t activeDetectors = grpData->GetDetectorMask();
    TString detStr = AliDAQ::ListOfTriggeredDetectors(activeDetectors);
    //printf("Detectors in the data:\n%s\n",detStr.Data());
    if ( detStr.Contains("TPC")==0){
      AliInfo("TPC not present in the run");
      return;
    }
  }

 


  entry          = GetCDBEntry("TPC/Calib/PadGainFactor");
  if (entry){
    //if (fPadGainFactor) delete fPadGainFactor;
    entry->SetOwner(kTRUE);
    fPadGainFactor = (AliTPCCalPad*)entry->GetObject();
  }else{
    AliFatal("TPC - Missing calibration entry TPC/Calib/PadGainFactor");
  }
  //
  entry          = GetCDBEntry("TPC/Calib/TimeGain");
  if (entry){
    //if (fTimeGainSplines) delete fTimeGainSplines;
    entry->SetOwner(kTRUE);
    fTimeGainSplines = (TObjArray*)entry->GetObject();
  }else{
    AliFatal("TPC - Missing calibration entry TPC/Calib/Timegain");
  }
  //
  entry          = GetCDBEntry("TPC/Calib/GainFactorDedx");
  if (entry){
    entry->SetOwner(kTRUE);
    fDedxGainFactor = (AliTPCCalPad*)entry->GetObject();
  }else{
    AliFatal("TPC - Missing calibration entry TPC/Calib/gainFactordEdx");
  }
  //
  entry          = GetCDBEntry("TPC/Calib/PadTime0");
  if (entry){
    //if (fPadTime0) delete fPadTime0;
    entry->SetOwner(kTRUE);
    fPadTime0 = (AliTPCCalPad*)entry->GetObject();
  }else{
    AliFatal("TPC - Missing calibration entry");
  }

  entry          = GetCDBEntry("TPC/Calib/Distortion");
  if (entry){
    //if (fPadTime0) delete fPadTime0;
    entry->SetOwner(kTRUE);
    fDistortionMap =dynamic_cast<TObjArray*>(entry->GetObject());
  }else{
    //AliFatal("TPC - Missing calibration entry")
  }


  //
  //
  entry          = GetCDBEntry("TPC/Calib/PadNoise");
  if (entry){
    //if (fPadNoise) delete fPadNoise;
    entry->SetOwner(kTRUE);
    fPadNoise = (AliTPCCalPad*)entry->GetObject();
  }else{
    AliFatal("TPC - Missing calibration entry");
  }

  entry          = GetCDBEntry("TPC/Calib/Pedestals");
  if (entry){
    //if (fPedestals) delete fPedestals;
    entry->SetOwner(kTRUE);
    fPedestals = (AliTPCCalPad*)entry->GetObject();
  }

  entry          = GetCDBEntry("TPC/Calib/Temperature");
  if (entry){
    //if (fTemperature) delete fTemperature;
    entry->SetOwner(kTRUE);
    fTemperature = (AliTPCSensorTempArray*)entry->GetObject();
  }

 
  entry          = GetCDBEntry("TPC/Calib/ClusterParam");
  if (entry){
    entry->SetOwner(kTRUE);
    fClusterParam = (AliTPCClusterParam*)(entry->GetObject());
  }else{
    AliFatal("TPC - Missing calibration entry");
  }

  entry          = GetCDBEntry("TPC/Calib/RecoParam");
  if (entry){
    //PH    entry->SetOwner(kTRUE);
    fRecoParamList = dynamic_cast<TObjArray*>(entry->GetObject());

  }else{
    AliFatal("TPC - Missing calibration entry TPC/Calib/RecoParam");
  }


  //ALTRO configuration data
  entry          = GetCDBEntry("TPC/Calib/AltroConfig");
  if (entry){
    entry->SetOwner(kTRUE);
    fALTROConfigData=(TObjArray*)(entry->GetObject());
  }else{
    AliFatal("TPC - Missing calibration entry");
  }

  //Calibration Pulser data
  entry          = GetCDBEntry("TPC/Calib/Pulser");
  if (entry){
    entry->SetOwner(kTRUE);
    fPulserData=(TObjArray*)(entry->GetObject());
  }

  //Calibration ION tail data
  entry          = GetCDBEntry("TPC/Calib/IonTail");
  if (entry){
    //delete fIonTailArray; fIonTailArray=NULL;
    entry->SetOwner(kTRUE);
     fIonTailArray=(TObjArray*)(entry->GetObject());
     fIonTailArray->SetOwner(); //own the keys
  }

  //CE data
  entry          = GetCDBEntry("TPC/Calib/CE");
  if (entry){
    entry->SetOwner(kTRUE);
    fCEData=(TObjArray*)(entry->GetObject());
  }
  //RAW calibration data
 //  entry          = GetCDBEntry("TPC/Calib/Raw");

  entry          = GetCDBEntry("TPC/Calib/Mapping");
  if (entry){
    //if (fPadNoise) delete fPadNoise;
    entry->SetOwner(kTRUE);
    TObjArray * array = dynamic_cast<TObjArray*>(entry->GetObject());
    if (array && array->GetEntriesFast()==6){
      fMapping = new AliTPCAltroMapping*[6];
      for (Int_t i=0; i<6; i++){
        fMapping[i] =  dynamic_cast<AliTPCAltroMapping*>(array->At(i));
      }
    }
  }

  //CTP calibration data
  entry          = GetCDBEntry("GRP/CTP/CTPtiming");
  if (entry){
    //entry->SetOwner(kTRUE);
    fCTPTimeParams=dynamic_cast<AliCTPTimeParams*>(entry->GetObject());
  }else{
    AliError("TPC - Missing calibration entry");
  }
  //TPC space point correction data
  entry          = GetCDBEntry("TPC/Calib/Correction");
  if (entry){
    //entry->SetOwner(kTRUE);
    fComposedCorrection=dynamic_cast<AliTPCCorrection*>(entry->GetObject());
    if (fComposedCorrection) fComposedCorrection->Init();
    fComposedCorrectionArray=dynamic_cast<TObjArray*>(entry->GetObject());
    if (fComposedCorrectionArray){
      for (Int_t i=0; i<fComposedCorrectionArray->GetEntries(); i++){
	AliTPCComposedCorrection* composedCorrection= dynamic_cast<AliTPCComposedCorrection*>(fComposedCorrectionArray->At(i));
	if (composedCorrection) {
	  composedCorrection->Init();
	  if (composedCorrection->GetCorrections()){
	    if (composedCorrection->GetCorrections()->FindObject("FitAlignTPC")){
	      fBHasAlignmentOCDB=kTRUE;
	    }
	  }
	}
      }
    }
  }else{
    AliError("TPC - Missing calibration entry-  TPC/Calib/Correction");
  }

  //RCU trigger config mode
  fMode=GetRCUTriggerConfig();
  //
  if (!fTransform) {
    fTransform=new AliTPCTransform();
    fTransform->SetCurrentRun(AliCDBManager::Instance()->GetRun());
  }

  // Chamber HV data
  // needs to be called before InitDeadMap
  UpdateChamberHighVoltageData();

  // Create Dead Channel Map
  InitDeadMap();

  // Calculate derived ALTRO information
  InitAltroData();

  //
  AliCDBManager::Instance()->SetCacheFlag(cdbCache); // reset original CDB cache
}


void AliTPCcalibDB::LoadCorrectionMaps()
{
  // TPC fast Chebyshev correction map, loaded on 1st demand
  AliCDBEntry* entry = GetCDBEntry("TPC/Calib/CorrectionMaps");
  if (entry) {
    fCorrectionMaps = dynamic_cast<TObjArray*>(entry->GetObject());
  }
  else{
    AliError("TPC - Missing calibration entry-  TPC/Calib/CorrectionMaps");
  }
}

void AliTPCcalibDB::UpdateNonRec(){
  /// Update/Load the parameters which are important for QA studies
  /// and not used yet for the reconstruction
  ///
  /// RAW calibration data

  AliCDBEntry * entry=0;
  entry          = GetCDBEntry("TPC/Calib/Raw");
  if (entry){
    entry->SetOwner(kTRUE);
    TObjArray *arr=dynamic_cast<TObjArray*>(entry->GetObject());
    if (arr) fCalibRaw=(AliTPCCalibRaw*)arr->At(0);
    else fCalibRaw = (AliTPCCalibRaw*)(entry->GetObject());
  }
  //QA calibration data
  entry          = GetCDBEntry("TPC/Calib/QA");
  if (entry){
    entry->SetOwner(kTRUE);
    fDataQA=dynamic_cast<AliTPCdataQA*>(entry->GetObject());
  }
  // High voltage
  if (fRun>=0 && !fVoltageArray.GetValue(Form("%i",fRun))){
    entry = AliCDBManager::Instance()->Get("TPC/Calib/HighVoltage",fRun);
    if (entry)  {
      fVoltageArray.Add(new TObjString(Form("%i",fRun)),entry->GetObject());
    }
  }

}

Bool_t AliTPCcalibDB::GetTailcancelationGraphs(Int_t sector, TGraphErrors ** graphRes, Float_t * indexAmpGraphs){

///   Read OCDB entry object of Iontail (TObjArray of TGraphErrors of TRFs)
///   Naming of the TRF objects is: "gr_<chamber_type>_<voltage>_<laser_track_angle>_<distance_to_COG>" --> "gr_iroc_1240_1_1"

  //Int_t run = fTransform->GetCurrentRunNumber();
  //SetRun(run);
  //Float_t rocVoltage = GetChamberHighVoltage(run,sector, -1);      // Get the voltage from OCDB with a getter (old function)
//   Float_t rocVoltage=GetChamberHighVoltageMedian(sector);                      // Get the voltage from OCDB, new function from Jens

  Int_t nominalVoltage = (sector<36) ? 1240 : 1470 ;     // nominal voltage of 2012 when the TRF functions were produced

  Float_t rocVoltage = nominalVoltage;

  if ( rocVoltage < nominalVoltage/2. || rocVoltage > nominalVoltage*2. )
  {
    AliInfo(Form("rocVoltage out of range: roc: %.2f, nominal: %i", rocVoltage, nominalVoltage));
    return kFALSE;
  }

  Int_t tempVoltage = 0;
  Int_t trackAngle  = 4;                                 // (1=first, 2=second, 3=third, 4=first+second, 5=all tracks) note: 3rd is distorted by low freq
  TString rocType   = (sector<36) ? "iroc" : "oroc";
  const Int_t ngraph=fIonTailArray->GetLast();

  // create array of voltages in order to select the proper TRF with closest voltage
  Int_t voltages[ngraph];     // array of voltages
  for (Int_t i=0; i<ngraph; i++){
    voltages[i]=0;
  }

  // loop over response functions in the TObjarray
  Int_t nvoltages=0;
  for (Int_t i=0;i<=ngraph;i++){

    // read the TRF object name in order to select proper TRF for the given sector
    TString objname(fIonTailArray->At(i)->GetName());
    if (!objname.Contains(rocType)) continue;

    TObjArray *objArr = objname.Tokenize("_");

    // select the roc type (IROC or OROC) and the trackAngle
    if ( atoi(static_cast<TObjString*>(objArr->At(3))->GetName())==trackAngle )
    {
      // Create the voltage array for proper voltage value selection
      voltages[nvoltages]=atoi(static_cast<TObjString*>(objArr->At(2))->GetName());
      nvoltages++;
    }
    delete objArr;
  }

  // find closest voltage value to ROC voltage (among the TRF' voltage array --> to select proper t.r.f.)
  Int_t ampIndex     = 0;
  Int_t diffVoltage  = TMath::Abs(rocVoltage - voltages[0]);
  for (Int_t k=0;k<ngraph;k++) {
    if (diffVoltage >= TMath::Abs(rocVoltage-voltages[k]) && voltages[k]!=0)
      {
        diffVoltage    = TMath::Abs(rocVoltage-voltages[k]);
        ampIndex   = k;
      }
  }
  tempVoltage = voltages[ampIndex];    // use closest voltage to current voltage
  //if (run<140000) tempVoltage = nominalVoltage;    // for 2010 data

  // assign TGraphErrors
  Int_t igraph=0;
  for (Int_t i=0; i<=ngraph; i++){

    // read TRFs for TObjArray and select the roc type (IROC or OROC) and the trackAngle
    TGraphErrors * trfObj = static_cast<TGraphErrors*>(fIonTailArray->At(i));
    TString objname(trfObj->GetName());
    if (!objname.Contains(rocType)) continue; //choose ROC type

    TObjArray *objArr1 = objname.Tokenize("_");

    // TRF eleminations
    TObjString* angleString = static_cast<TObjString*>(objArr1->At(3));
    TObjString* voltageString = static_cast<TObjString*>(objArr1->At(2));
    //choose angle and voltage
    if ((atoi(angleString->GetName())==trackAngle) && (atoi(voltageString->GetName())==tempVoltage) )
    {
      // Apply Voltage scaling
      Int_t voltage       = atoi(voltageString->GetName());
      Double_t voltageScaled = 1;
      if (rocVoltage>0)  voltageScaled = Double_t(voltage)/Double_t(rocVoltage); // for jens how it can happen that we have clusters at 0 HV ?
      const Int_t nScaled          = TMath::Nint(voltageScaled*trfObj->GetN())-1;
      Double_t x;
      Double_t y;

      delete graphRes[igraph];
      graphRes[igraph]       = new TGraphErrors(nScaled);

      for (Int_t j=0; j<nScaled; j++){
        x = TMath::Nint(j*(voltageScaled));
        y = (j<trfObj->GetN()) ? (1./voltageScaled)*trfObj->GetY()[j] : 0.;
        graphRes[igraph]->SetPoint(j,x,y);
      }

      // fill arrays for proper position and amplitude selections
      TObjString* distanceToCenterOfGravity = static_cast<TObjString*>(objArr1->At(4));
      indexAmpGraphs[igraph] = (distanceToCenterOfGravity->GetString().Atof())/10.;
      // smooth voltage scaled graph
      for (Int_t m=1; m<nScaled;m++){
        if (graphRes[igraph]->GetY()[m]==0) graphRes[igraph]->GetY()[m] = graphRes[igraph]->GetY()[m-1];
      }
      igraph++;
    }
  delete objArr1;
  }
  return kTRUE;
}

void AliTPCcalibDB::CreateObjectList(const Char_t *filename, TObjArray *calibObjects)
{
/// Create calibration objects and read contents from OCDB

   if ( calibObjects == 0x0 ) return;
   ifstream in;
   in.open(filename);
   if ( !in.is_open() ){
      fprintf(stderr,"Error: cannot open list file '%s'", filename);
      return;
   }

   AliTPCCalPad *calPad=0x0;

   TString sFile;
   sFile.ReadFile(in);
   in.close();

   TObjArray *arrFileLine = sFile.Tokenize("\n");

   TIter nextLine(arrFileLine);

   TObjString *sObjLine=0x0;
   while ( (sObjLine = (TObjString*)nextLine()) ){
      TString sLine(sObjLine->GetString());

      TObjArray *arrNextCol = sLine.Tokenize("\t");

      TObjString *sObjType     = (TObjString*)(arrNextCol->At(0));
      TObjString *sObjFileName = (TObjString*)(arrNextCol->At(1));
      delete arrNextCol;

      if ( !sObjType || ! sObjFileName ) continue;
      TString sType(sObjType->GetString());
      TString sFileName(sObjFileName->GetString());
//       printf("%s\t%s\n",sType.Data(),sFileName.Data());

      TFile *fIn = TFile::Open(sFileName);
      if ( !fIn ){
         fprintf(stderr,"File not found: '%s'", sFileName.Data());
         continue;
      }

      if ( sType == "CE" ){
         AliTPCCalibCE *ce = (AliTPCCalibCE*)fIn->Get("AliTPCCalibCE");

         calPad = new AliTPCCalPad((TObjArray*)ce->GetCalPadT0());
         calPad->SetNameTitle("CETmean","CETmean");
         calibObjects->Add(calPad);

         calPad = new AliTPCCalPad((TObjArray*)ce->GetCalPadQ());
         calPad->SetNameTitle("CEQmean","CEQmean");
         calibObjects->Add(calPad);

         calPad = new AliTPCCalPad((TObjArray*)ce->GetCalPadRMS());
         calPad->SetNameTitle("CETrms","CETrms");
         calibObjects->Add(calPad);

      } else if ( sType == "Pulser") {
         AliTPCCalibPulser *sig = (AliTPCCalibPulser*)fIn->Get("AliTPCCalibPulser");

         calPad = new AliTPCCalPad((TObjArray*)sig->GetCalPadT0());
         calPad->SetNameTitle("PulserTmean","PulserTmean");
         calibObjects->Add(calPad);

         calPad = new AliTPCCalPad((TObjArray*)sig->GetCalPadQ());
         calPad->SetNameTitle("PulserQmean","PulserQmean");
         calibObjects->Add(calPad);

         calPad = new AliTPCCalPad((TObjArray*)sig->GetCalPadRMS());
         calPad->SetNameTitle("PulserTrms","PulserTrms");
         calibObjects->Add(calPad);

      } else if ( sType == "Pedestals") {
         AliTPCCalibPedestal *ped = (AliTPCCalibPedestal*)fIn->Get("AliTPCCalibPedestal");

         calPad = new AliTPCCalPad((TObjArray*)ped->GetCalPadPedestal());
         calPad->SetNameTitle("Pedestals","Pedestals");
         calibObjects->Add(calPad);

         calPad = new AliTPCCalPad((TObjArray*)ped->GetCalPadRMS());
         calPad->SetNameTitle("Noise","Noise");
         calibObjects->Add(calPad);

      } else {
         fprintf(stderr,"Undefined Type: '%s'",sType.Data());

      }
      delete fIn;
   }
   delete arrFileLine;
}

Int_t AliTPCcalibDB::InitDeadMap()
{
  /// Initialize DeadChannel Map
  /// Source of information:
  /// -  HV (see UpdateChamberHighVoltageData())
  /// -  Altro disabled channels. Noisy channels.
  /// -  DDL list
  ///
  /// List of required OCDB Entries (See also UpdateChamberHighVoltageData())
  /// - TPC/Calib/AltroConfig
  /// - TPC/Calib/HighVoltage

  // check necessary information
  const Int_t run=GetRun();
  if (run<0){
    AliError("run not set in CDB manager. Cannot create active channel map");
    return 0;
  }
  AliDCSSensorArray* voltageArray = GetVoltageSensors(run);
  AliTPCCalPad*          altroMap = GetALTROMasked();
  TMap*                    mapddl = GetDDLMap();

  if (!voltageArray && !altroMap && !mapddl) {
    AliError("All necessary information to create the activate channel are map missing.");
    AliError(" -> Check existance of the OCDB entries: 'TPC/Calib/AltroConfig', 'TPC/Calib/HighVoltage'");
    return 0;
  }

  // mapping handler
  AliTPCmapper map(gSystem->ExpandPathName("$ALICE_ROOT/TPC/mapping/"));

  //=============================================================
  // Setup DDL map

  Bool_t ddlMap[216]={0};
  for (Int_t iddl=0; iddl<216; ++iddl) ddlMap[iddl]=1;
  if (mapddl){
    TObjString *s = (TObjString*)mapddl->GetValue("DDLArray");
    if (s){
      for (Int_t iddl=0; iddl<216; ++iddl) {
        ddlMap[iddl]=TString(s->GetString()(iddl))!="0";
        if (!ddlMap[iddl]) {
          Int_t roc = map.GetRocFromEquipmentID(iddl+768);
          AliWarning(Form("Inactive DDL (#%d, ROC %2d) detected based on the 'DDLArray' in 'TPC/Calib/AltroConfig'. This will deactivate many channels.", iddl, roc));
        }
      }
    }
  } else {
    AliError("DDL map missing. ActiveChannelMap can only be created with parts of the information.");
    AliError(" -> Check existance of 'DDLArray' in the OCDB entry: 'TPC/Calib/AltroConfig'");
  }
  // Setup DDL map done
  // ============================================================

  // ============================================================
  // Set up channel masking from correction maps

  TBits maskedPads[72];
  GetMaskedChannelsFromCorrectionMaps(maskedPads);

  // channel masking done
  // ============================================================

  //=============================================================
  // Setup active chnnel map
  //

  if (!fActiveChannelMap) fActiveChannelMap=new AliTPCCalPad("ActiveChannelMap","ActiveChannelMap");

  if (!altroMap) {
    AliError("ALTRO dead channel map missing. ActiveChannelMap can only be created with parts of the information.");
    AliError(" -> Check existance of 'Masked' in the OCDB entry: 'TPC/Calib/AltroConfig'");
  }

  AliTPCROC* tpcROC = AliTPCROC::Instance();

  for (Int_t iROC=0;iROC<AliTPCCalPad::kNsec;++iROC){
    AliTPCCalROC *roc=fActiveChannelMap->GetCalROC(iROC);
    if (!roc){
      AliError(Form("No ROC %d in active channel map",iROC));
      continue;
    }

    // check for bad voltage
    // see UpdateChamberHighVoltageData()
    if (!fChamberHVStatus[iROC]){
      roc->Multiply(0.);
      AliWarning(Form("Turning off all channels of ROC %2d due to a bad HV status", iROC));
      AliWarning(" -> Check messages in UpdateChamberHighVoltageData()");
      continue;
    }

    AliTPCCalROC *masked=0x0;
    if (altroMap) masked=altroMap->GetCalROC(iROC);

    Int_t numberOfDeactivatedChannels=0;
    for (UInt_t irow=0; irow<roc->GetNrows(); ++irow){
      // first channel in row
      const Int_t channel0 = tpcROC->GetRowIndexes(iROC)[irow];

      for (UInt_t ipad=0; ipad<roc->GetNPads(irow); ++ipad){
        //per default the channel is on
        roc->SetValue(irow,ipad,1);

        // apply altro dead channel mask (inverse logik, it is not active, but inactive channles)
        if (masked && masked->GetValue(irow, ipad)) roc->SetValue(irow, ipad ,0);

        // mask channels if a DDL is inactive
        Int_t ddlId=map.GetEquipmentID(iROC, irow, ipad)-768;
        if (ddlId>=0 && !ddlMap[ddlId]) roc->SetValue(irow, ipad ,0);

        // mask channels if error on space point coorection is too large
        if (maskedPads && maskedPads[iROC].TestBitNumber(channel0+ipad)) roc->SetValue(irow, ipad ,0);

        if (roc->GetValue(irow, ipad)<0.0001) {
          ++numberOfDeactivatedChannels;
        }
      }
    }

    if (numberOfDeactivatedChannels>0) {
      AliInfo(Form("Deactivated %4d channels in ROC %2d due to altro and DDL map states",
                   numberOfDeactivatedChannels, iROC));
    }
  }

  return 1;
}

void AliTPCcalibDB::InitAltroData()
{
  /// Initialise derived ALTRO data
  ///
  /// List of required OCDB Entries
  /// - TPC/Calib/AltroConfig
  /// - TPC/Calib/Parameters

  // ===| Maximum time bin |====================================================
  //
  // Calculate the maximum time using the 'AcqStart' cal pad object from
  // TPC/Calib/AltroConfig
  // if this object is not available, the value will return the max time bin
  // stored in the AliTPCParam object from TPC/Calib/Parameters

  fMaxTimeBinAllPads=-1;

  const AliTPCCalPad *calPadAcqStop = GetALTROAcqStop();

  if (calPadAcqStop) {
    //find max elememt
    // TODO: change this once the GetMaxElement function is implemented in AliTPCCalPad
    Float_t maxBin=-1;
    for (Int_t iroc=0; iroc<AliTPCCalPad::kNsec; ++iroc) {
      const AliTPCCalROC *roc = calPadAcqStop->GetCalROC(iroc);
      if (!roc) continue;
      for (Int_t ichannel=0; ichannel<roc->GetNchannels(); ++ichannel) {
        const Float_t val=roc->GetValue(ichannel);
        if (val>maxBin) maxBin=val;
      }
    }
    fMaxTimeBinAllPads = TMath::Nint(maxBin);
  }

  if (fMaxTimeBinAllPads<0) {
    if (fParam) {
      AliWarning("Could not access 'AcqStop' map from AltroConfig or invalid max time bine. fMaxTimeBinAllPads will be set from AliTPCParam.");
      fMaxTimeBinAllPads = fParam->GetMaxTBin();
    } else {
      // last fallback
      AliWarning("Could neither access 'Parameters' nor 'AcqStop' map from AltroConfig. fMaxTimeBinAllPads will be set to 1000.");
      fMaxTimeBinAllPads=1000;
    }
  }

  AliInfo(TString::Format("fMaxTimeBinAllPads set to %d", fMaxTimeBinAllPads).Data());
}

void AliTPCcalibDB::MakeTree(const char * fileName, TObjArray * array, const char * mapFileName, AliTPCCalPad* outlierPad, Float_t ltmFraction) {
  /// Write a tree with all available information
  /// if mapFileName is specified, the Map information are also written to the tree
  /// pads specified in outlierPad are not used for calculating statistics
  ///  - the same function as AliTPCCalPad::MakeTree -

   AliTPCROC* tpcROCinstance = AliTPCROC::Instance();

   TObjArray* mapIROCs = 0;
   TObjArray* mapOROCs = 0;
   TVectorF *mapIROCArray = 0;
   TVectorF *mapOROCArray = 0;
   Int_t mapEntries = 0;
   TString* mapNames = 0;

   if (mapFileName) {
      TFile mapFile(mapFileName, "read");

      TList* listOfROCs = mapFile.GetListOfKeys();
      mapEntries = listOfROCs->GetEntries()/2;
      mapIROCs = new TObjArray(mapEntries*2);
      mapOROCs = new TObjArray(mapEntries*2);
      mapIROCArray = new TVectorF[mapEntries];
      mapOROCArray = new TVectorF[mapEntries];

      mapNames = new TString[mapEntries];
      for (Int_t ivalue = 0; ivalue < mapEntries; ivalue++) {
	TString nameROC(((TKey*)(listOfROCs->At(ivalue*2)))->GetName());
         nameROC.Remove(nameROC.Length()-4, 4);
         mapIROCs->AddAt((AliTPCCalROC*)mapFile.Get((nameROC + "IROC").Data()), ivalue);
         mapOROCs->AddAt((AliTPCCalROC*)mapFile.Get((nameROC + "OROC").Data()), ivalue);
         mapNames[ivalue].Append(nameROC);
      }

      for (Int_t ivalue = 0; ivalue < mapEntries; ivalue++) {
         mapIROCArray[ivalue].ResizeTo(tpcROCinstance->GetNChannels(0));
         mapOROCArray[ivalue].ResizeTo(tpcROCinstance->GetNChannels(36));

         for (UInt_t ichannel = 0; ichannel < tpcROCinstance->GetNChannels(0); ichannel++)
            (mapIROCArray[ivalue])[ichannel] = ((AliTPCCalROC*)(mapIROCs->At(ivalue)))->GetValue(ichannel);
         for (UInt_t ichannel = 0; ichannel < tpcROCinstance->GetNChannels(36); ichannel++)
            (mapOROCArray[ivalue])[ichannel] = ((AliTPCCalROC*)(mapOROCs->At(ivalue)))->GetValue(ichannel);
      }

   } //  if (mapFileName)

   TTreeSRedirector cstream(fileName);
   Int_t arrayEntries = array->GetEntries();

   TString* names = new TString[arrayEntries];
   for (Int_t ivalue = 0; ivalue < arrayEntries; ivalue++)
      names[ivalue].Append(((AliTPCCalPad*)array->At(ivalue))->GetName());

   for (UInt_t isector = 0; isector < tpcROCinstance->GetNSectors(); isector++) {
      //
      // get statistic for given sector
      //
      TVectorF median(arrayEntries);
      TVectorF mean(arrayEntries);
      TVectorF rms(arrayEntries);
      TVectorF ltm(arrayEntries);
      TVectorF ltmrms(arrayEntries);
      TVectorF medianWithOut(arrayEntries);
      TVectorF meanWithOut(arrayEntries);
      TVectorF rmsWithOut(arrayEntries);
      TVectorF ltmWithOut(arrayEntries);
      TVectorF ltmrmsWithOut(arrayEntries);

      TVectorF *vectorArray = new TVectorF[arrayEntries];
      for (Int_t ivalue = 0; ivalue < arrayEntries; ivalue++)
         vectorArray[ivalue].ResizeTo(tpcROCinstance->GetNChannels(isector));

      for (Int_t ivalue = 0; ivalue < arrayEntries; ivalue++) {
         AliTPCCalPad* calPad = (AliTPCCalPad*) array->At(ivalue);
         AliTPCCalROC* calROC = calPad->GetCalROC(isector);
         AliTPCCalROC* outlierROC = 0;
         if (outlierPad) outlierROC = outlierPad->GetCalROC(isector);
         if (calROC) {
            median[ivalue] = calROC->GetMedian();
            mean[ivalue] = calROC->GetMean();
            rms[ivalue] = calROC->GetRMS();
            Double_t ltmrmsValue = 0;
            ltm[ivalue] = calROC->GetLTM(&ltmrmsValue, ltmFraction);
            ltmrms[ivalue] = ltmrmsValue;
            if (outlierROC) {
               medianWithOut[ivalue] = calROC->GetMedian(outlierROC);
               meanWithOut[ivalue] = calROC->GetMean(outlierROC);
               rmsWithOut[ivalue] = calROC->GetRMS(outlierROC);
               ltmrmsValue = 0;
               ltmWithOut[ivalue] = calROC->GetLTM(&ltmrmsValue, ltmFraction, outlierROC);
               ltmrmsWithOut[ivalue] = ltmrmsValue;
            }
         }
         else {
            median[ivalue] = 0.;
            mean[ivalue] = 0.;
            rms[ivalue] = 0.;
            ltm[ivalue] = 0.;
            ltmrms[ivalue] = 0.;
            medianWithOut[ivalue] = 0.;
            meanWithOut[ivalue] = 0.;
            rmsWithOut[ivalue] = 0.;
            ltmWithOut[ivalue] = 0.;
            ltmrmsWithOut[ivalue] = 0.;
         }
      }

      //
      // fill vectors of variable per pad
      //
      TVectorF *posArray = new TVectorF[8];
      for (Int_t ivalue = 0; ivalue < 8; ivalue++)
         posArray[ivalue].ResizeTo(tpcROCinstance->GetNChannels(isector));

      Float_t posG[3] = {0};
      Float_t posL[3] = {0};
      Int_t ichannel = 0;
      for (UInt_t irow = 0; irow < tpcROCinstance->GetNRows(isector); irow++) {
         for (UInt_t ipad = 0; ipad < tpcROCinstance->GetNPads(isector, irow); ipad++) {
            tpcROCinstance->GetPositionLocal(isector, irow, ipad, posL);
            tpcROCinstance->GetPositionGlobal(isector, irow, ipad, posG);
            posArray[0][ichannel] = irow;
            posArray[1][ichannel] = ipad;
            posArray[2][ichannel] = posL[0];
            posArray[3][ichannel] = posL[1];
            posArray[4][ichannel] = posG[0];
            posArray[5][ichannel] = posG[1];
            posArray[6][ichannel] = (Int_t)(ipad - (Double_t)(tpcROCinstance->GetNPads(isector, irow))/2);
            posArray[7][ichannel] = ichannel;

            // loop over array containing AliTPCCalPads
            for (Int_t ivalue = 0; ivalue < arrayEntries; ivalue++) {
               AliTPCCalPad* calPad = (AliTPCCalPad*) array->At(ivalue);
               AliTPCCalROC* calROC = calPad->GetCalROC(isector);
               if (calROC)
                  (vectorArray[ivalue])[ichannel] = calROC->GetValue(irow, ipad);
               else
                  (vectorArray[ivalue])[ichannel] = 0;
            }
            ichannel++;
         }
      }

      cstream << "calPads" <<
         "sector=" << isector;

      for  (Int_t ivalue = 0; ivalue < arrayEntries; ivalue++) {
         cstream << "calPads" <<
            (Char_t*)((names[ivalue] + "_Median=").Data()) << median[ivalue] <<
            (Char_t*)((names[ivalue] + "_Mean=").Data()) << mean[ivalue] <<
            (Char_t*)((names[ivalue] + "_RMS=").Data()) << rms[ivalue] <<
            (Char_t*)((names[ivalue] + "_LTM=").Data()) << ltm[ivalue] <<
            (Char_t*)((names[ivalue] + "_RMS_LTM=").Data()) << ltmrms[ivalue];
         if (outlierPad) {
            cstream << "calPads" <<
               (Char_t*)((names[ivalue] + "_Median_OutlierCutted=").Data()) << medianWithOut[ivalue] <<
               (Char_t*)((names[ivalue] + "_Mean_OutlierCutted=").Data()) << meanWithOut[ivalue] <<
               (Char_t*)((names[ivalue] + "_RMS_OutlierCutted=").Data()) << rmsWithOut[ivalue] <<
               (Char_t*)((names[ivalue] + "_LTM_OutlierCutted=").Data()) << ltmWithOut[ivalue] <<
               (Char_t*)((names[ivalue] + "_RMS_LTM_OutlierCutted=").Data()) << ltmrmsWithOut[ivalue];
         }
      }

      for  (Int_t ivalue = 0; ivalue < arrayEntries; ivalue++) {
         cstream << "calPads" <<
            (Char_t*)((names[ivalue] + ".=").Data()) << &vectorArray[ivalue];
      }

      if (mapFileName) {
         for  (Int_t ivalue = 0; ivalue < mapEntries; ivalue++) {
            if (isector < 36)
               cstream << "calPads" <<
                  (Char_t*)((mapNames[ivalue] + ".=").Data()) << &mapIROCArray[ivalue];
            else
               cstream << "calPads" <<
                  (Char_t*)((mapNames[ivalue] + ".=").Data()) << &mapOROCArray[ivalue];
         }
      }

      cstream << "calPads" <<
         "row.=" << &posArray[0] <<
         "pad.=" << &posArray[1] <<
         "lx.=" << &posArray[2] <<
         "ly.=" << &posArray[3] <<
         "gx.=" << &posArray[4] <<
         "gy.=" << &posArray[5] <<
         "rpad.=" << &posArray[6] <<
         "channel.=" << &posArray[7];

      cstream << "calPads" <<
         "\n";

      delete[] posArray;
      delete[] vectorArray;
   }


   delete[] names;
   if (mapFileName) {
      delete mapIROCs;
      delete mapOROCs;
      delete[] mapIROCArray;
      delete[] mapOROCArray;
      delete[] mapNames;
   }
}

Int_t AliTPCcalibDB::GetRCUTriggerConfig() const
{
  /// return the RCU trigger configuration register

  TMap *map=GetRCUconfig();
  if (!map) return -1;
  TVectorF *v=(TVectorF*)map->GetValue("TRGCONF_TRG_MODE");
  Float_t mode=-1;
  for (Int_t i=0; i<v->GetNrows(); ++i){
    Float_t newmode=v->GetMatrixArray()[i];
    if (newmode>-1){
      if (mode>-1&&newmode!=mode) AliWarning("Found different RCU trigger configurations!!!");
      mode=newmode;
    }
  }
  return (Int_t)mode;
}

Bool_t AliTPCcalibDB::IsTrgL0()
{
  /// return if the FEE readout was triggered on L0

  if (fMode<0) return kFALSE;
  return (fMode==1);
}

Bool_t AliTPCcalibDB::IsTrgL1()
{
  /// return if the FEE readout was triggered on L1

  if (fMode<0) return kFALSE;
  return (fMode==0);
}

void AliTPCcalibDB::RegisterExB(Int_t index, Float_t bz, Bool_t bdelete){
  /// Register static ExB correction map
  /// index - registration index - used for visualization
  /// bz    - bz field in kGaus

  //  Float_t factor =  bz/(-5.);  // default b filed in Cheb with minus sign
  Float_t factor =  bz/(5.);  // default b filed in Cheb with minus sign
                              // was chenged in the Revision ???? (Ruben can you add here number)

  AliMagF*   bmap = new AliMagF("MapsExB","MapsExB", factor,TMath::Sign(1.f,factor),AliMagF::k5kG);

  AliTPCExBFirst *exb  = new  AliTPCExBFirst(bmap,0.88*2.6400e+04,50,50,50);
  AliTPCExB::SetInstance(exb);

  if (bdelete){
    delete bmap;
  }else{
    AliTPCExB::RegisterField(index,bmap);
  }
  if (index>=fgExBArray.GetEntries()) fgExBArray.Expand((index+1)*2+11);
  fgExBArray.AddAt(exb,index);
}


AliTPCExB*    AliTPCcalibDB::GetExB(Float_t bz, Bool_t deleteB) {
  /// bz filed in KGaus not in tesla
  /// Get ExB correction map
  /// if doesn't exist - create it

  Int_t index = TMath::Nint(5+bz);
  if (index>fgExBArray.GetEntries()) fgExBArray.Expand((index+1)*2+11);
  if (!fgExBArray.At(index)) AliTPCcalibDB::RegisterExB(index,bz,deleteB);
  return (AliTPCExB*)fgExBArray.At(index);
}


void  AliTPCcalibDB::SetExBField(Float_t bz){
  /// Set magnetic filed for ExB correction

  fExB = GetExB(bz,kFALSE);
}

void  AliTPCcalibDB::SetExBField(const AliMagF*   bmap){
  /// Set magnetic field for ExB correction

  AliTPCExBFirst *exb  = new  AliTPCExBFirst(bmap,0.88*2.6400e+04,50,50,50);
  AliTPCExB::SetInstance(exb);
  fExB=exb;
}



void AliTPCcalibDB::UpdateRunInformations( Int_t run, Bool_t force){
  /// - > Don't use it for reconstruction - Only for Calibration studies

  if (run<=0) return;
  TObjString runstr(Form("%i",run));
  fRun=run;
  AliCDBEntry * entry = 0;
  if (run>= fRunList.fN){
    fRunList.Set(run*2+1);
    //
    //
    fALTROConfigData->Expand(run*2+1);    // ALTRO configuration data
    fPulserData->Expand(run*2+1);         // Calibration Pulser data
    fCEData->Expand(run*2+1);             // CE data
    if (!fTimeGainSplines) fTimeGainSplines = new TObjArray(run*2+1);
    fTimeGainSplines->Expand(run*2+1); // Array of AliSplineFits: at 0 MIP position in
  }
  if (fRunList[run]>0 &&force==kFALSE) return;

  fRunList[run]=1;  // sign as used

  //
  entry = AliCDBManager::Instance()->Get("GRP/GRP/Data",run);
  if (entry)  {
    AliGRPObject * grpRun = dynamic_cast<AliGRPObject*>(entry->GetObject());
    if (!grpRun){
      TMap* map = dynamic_cast<TMap*>(entry->GetObject());
      if (map){
	//grpRun = new AliGRPObject;
	//grpRun->ReadValuesFromMap(map);
	grpRun =  MakeGRPObjectFromMap(map);

	fGRPMaps.Add(new TObjString(runstr),map);
      }
    }
    fGRPArray.Add(new TObjString(runstr),(AliGRPObject*)grpRun->Clone());
  }
  entry = AliCDBManager::Instance()->Get("TPC/Calib/Goofie",run);
  if (entry){
    fGoofieArray.Add(new TObjString(runstr),entry->GetObject());
  }
  //

  //
  entry = AliCDBManager::Instance()->Get("TPC/Calib/TimeGain",run);
  if (entry)  {
    fTimeGainSplinesArray.Add(new TObjString(runstr),entry->GetObject());
  }else{
    AliFatal("TPC - Missing calibration entry TimeGain");
  }
  //
  entry = AliCDBManager::Instance()->Get("TPC/Calib/TimeDrift",run);
  if (entry)  {
    TObjArray * timeArray = (TObjArray*)entry->GetObject();
    fDriftCorrectionArray.Add(new TObjString(runstr),entry->GetObject());
    AliTPCCorrection * correctionTime = (AliTPCCorrection *)timeArray->FindObject("FitCorrectionTime");
    if (correctionTime && fComposedCorrectionArray){
      correctionTime->Init();
      if (fComposedCorrectionArray->GetEntriesFast()<4) fComposedCorrectionArray->Expand(40);
      fComposedCorrectionArray->AddAt(correctionTime,4); //add time dependent correction to the list of available corrections
    }
  }else{
    AliFatal("TPC - Missing calibration entry TimeDrift");
  }
  //
  entry = AliCDBManager::Instance()->Get("TPC/Calib/Temperature",run);
  if (entry)  {
    fTemperatureArray.Add(new TObjString(runstr),entry->GetObject());
  }

  // High voltage
  entry = AliCDBManager::Instance()->Get("TPC/Calib/HighVoltage",run);
  if (!fVoltageArray.GetValue(runstr.GetName()) && entry)  {
    fVoltageArray.Add(new TObjString(runstr),entry->GetObject());
  }

  //apply fDButil filters

  fDButil->UpdateFromCalibDB();
  if (fTemperature) fDButil->FilterTemperature(fTemperature);

  AliDCSSensor * press = GetPressureSensor(run,0);
  AliTPCSensorTempArray * temp = GetTemperatureSensor(run);
  Bool_t accept=kTRUE;
  if (temp) {
    accept = fDButil->FilterTemperature(temp)>0.1;
  }
  if (press) {
    const Double_t kMinP=900.;
    const Double_t kMaxP=1050.;
    const Double_t kMaxdP=10.;
    const Double_t kSigmaCut=4.;
    fDButil->FilterSensor(press,kMinP,kMaxP,kMaxdP,kSigmaCut);
    if (press->GetFit()==0) accept=kFALSE;
  }

  if (press && temp &&accept){
    AliTPCCalibVdrift * vdrift = new AliTPCCalibVdrift(temp, press,0);
    fVdriftArray.Add(new TObjString(runstr),vdrift);
  }

  fDButil->FilterCE(120., 3., 4.,0);
  fDButil->FilterTracks(run, 10.,0);

}


Float_t AliTPCcalibDB::GetGain(Int_t sector, Int_t row, Int_t pad){
  /// Get Gain factor for given pad

  AliTPCCalPad *calPad = Instance()->fDedxGainFactor;;
  if (!calPad) return 0;
  return calPad->GetCalROC(sector)->GetValue(row,pad);
}

AliSplineFit* AliTPCcalibDB::GetVdriftSplineFit(const char* name, Int_t run){
  /// GetDrift velocity spline fit

  TObjArray *arr=GetTimeVdriftSplineRun(run);
  if (!arr) return 0;
  return dynamic_cast<AliSplineFit*>(arr->FindObject(name));
}

AliSplineFit* AliTPCcalibDB::CreateVdriftSplineFit(const char* graphName, Int_t run){
  /// create spline fit from the drift time graph in TimeDrift

  TObjArray *arr=GetTimeVdriftSplineRun(run);
  if (!arr) return 0;
  TGraph *graph=dynamic_cast<TGraph*>(arr->FindObject(graphName));
  if (!graph) return 0;
  AliSplineFit *fit = new AliSplineFit();
  fit->SetGraph(graph);
  fit->SetMinPoints(graph->GetN()+1);
  fit->InitKnots(graph,2,0,0.001);
  fit->SplineFit(0);
  return fit;
}

AliGRPObject *AliTPCcalibDB::GetGRP(Int_t run){
  /// Get GRP object for given run

  AliGRPObject * grpRun = dynamic_cast<AliGRPObject *>((Instance()->fGRPArray).GetValue(Form("%i",run)));
  if (!grpRun) {
    Instance()->UpdateRunInformations(run);
    grpRun = dynamic_cast<AliGRPObject *>(Instance()->fGRPArray.GetValue(Form("%i",run)));
    if (!grpRun) return 0;
  }
  return grpRun;
}

TMap *  AliTPCcalibDB::GetGRPMap(Int_t run){
  /// Get GRP map for given run

  TMap * grpRun = dynamic_cast<TMap *>((Instance()->fGRPMaps).GetValue(Form("%i",run)));
  if (!grpRun) {
    Instance()->UpdateRunInformations(run);
    grpRun = dynamic_cast<TMap *>(Instance()->fGRPMaps.GetValue(Form("%i",run)));
    if (!grpRun) return 0;
  }
  return grpRun;
}


AliDCSSensor * AliTPCcalibDB::GetPressureSensor(Int_t run, Int_t type){
  /// Get Pressure sensor
  /// run  = run number
  /// type = 0 - Cavern pressure
  ///        1 - Suface pressure
  /// First try to get if trom map - if existing  (Old format of data storing)


  TMap *map = GetGRPMap(run);
  if (map){
    AliDCSSensor * sensor = 0;
    TObject *osensor=0;
    if (type==0) osensor = ((*map)("fCavernPressure"));
    if (type==1) osensor = ((*map)("fP2Pressure"));
    sensor =dynamic_cast<AliDCSSensor *>(osensor);
    if (sensor) return sensor;
  }
  //
  // If not map try to get it from the GRPObject
  //
  AliGRPObject * grpRun = dynamic_cast<AliGRPObject *>(fGRPArray.GetValue(Form("%i",run)));
  if (!grpRun) {
    UpdateRunInformations(run);
    grpRun = dynamic_cast<AliGRPObject *>(fGRPArray.GetValue(Form("%i",run)));
    if (!grpRun) return 0;
  }
  AliDCSSensor * sensor = grpRun->GetCavernAtmosPressure();
  if (type==1) sensor = grpRun->GetSurfaceAtmosPressure();
  return sensor;
}

AliTPCSensorTempArray * AliTPCcalibDB::GetTemperatureSensor(Int_t run){
  /// Get temperature sensor array

  AliTPCSensorTempArray * tempArray = (AliTPCSensorTempArray *)fTemperatureArray.GetValue(Form("%i",run));
  if (!tempArray) {
    UpdateRunInformations(run);
    tempArray = (AliTPCSensorTempArray *)fTemperatureArray.GetValue(Form("%i",run));
  }
  return tempArray;
}


TObjArray * AliTPCcalibDB::GetTimeGainSplinesRun(Int_t run){
  /// Get temperature sensor array

  TObjArray * gainSplines = (TObjArray *)fTimeGainSplinesArray.GetValue(Form("%i",run));
  if (!gainSplines) {
    UpdateRunInformations(run);
    gainSplines = (TObjArray *)fTimeGainSplinesArray.GetValue(Form("%i",run));
  }
  return gainSplines;
}

TObjArray * AliTPCcalibDB::GetTimeVdriftSplineRun(Int_t run){
  /// Get drift spline array

  TObjArray * driftSplines = (TObjArray *)fDriftCorrectionArray.GetValue(Form("%i",run));
  if (!driftSplines) {
    UpdateRunInformations(run);
    driftSplines = (TObjArray *)fDriftCorrectionArray.GetValue(Form("%i",run));
  }
  return driftSplines;
}

AliDCSSensorArray * AliTPCcalibDB::GetVoltageSensors(Int_t run){
  /// Get temperature sensor array

  AliDCSSensorArray * voltageArray = (AliDCSSensorArray *)fVoltageArray.GetValue(Form("%i",run));
  if (!voltageArray) {
    UpdateRunInformations(run);
    voltageArray = (AliDCSSensorArray *)fVoltageArray.GetValue(Form("%i",run));
  }
  return voltageArray;
}

AliDCSSensorArray * AliTPCcalibDB::GetGoofieSensors(Int_t run){
  /// Get temperature sensor array

  AliDCSSensorArray * goofieArray = (AliDCSSensorArray *)fGoofieArray.GetValue(Form("%i",run));
  if (!goofieArray) {
    UpdateRunInformations(run);
    goofieArray = (AliDCSSensorArray *)fGoofieArray.GetValue(Form("%i",run));
  }
  return goofieArray;
}



AliTPCCalibVdrift *     AliTPCcalibDB::GetVdrift(Int_t run){
  /// Get the interface to the the vdrift

  AliTPCCalibVdrift  * vdrift = (AliTPCCalibVdrift*)fVdriftArray.GetValue(Form("%i",run));
  if (!vdrift) {
    UpdateRunInformations(run);
    vdrift= (AliTPCCalibVdrift*)fVdriftArray.GetValue(Form("%i",run));
  }
  return vdrift;
}

Float_t AliTPCcalibDB::GetCEdriftTime(Int_t run, Int_t sector, Double_t timeStamp, Int_t *entries)
{
  /// GetCE drift time information for 'sector'
  /// sector 72 is the mean drift time of the A-Side
  /// sector 73 is the mean drift time of the C-Side
  /// it timestamp==-1 return mean value

  AliTPCcalibDB::Instance()->SetRun(run);
  TGraph *gr=AliTPCcalibDB::Instance()->GetCErocTgraph(sector);
  if (!gr||sector<0||sector>73) {
    if (entries) *entries=0;
    return 0.;
  }
  Float_t val=0.;
  if (timeStamp==-1.){
    val=gr->GetMean(2);
  }else{
    for (Int_t ipoint=0;ipoint<gr->GetN();++ipoint){
      Double_t x,y;
      gr->GetPoint(ipoint,x,y);
      if (x<timeStamp) continue;
      val=y;
      break;
    }
  }
  return val;
}

Float_t AliTPCcalibDB::GetCEchargeTime(Int_t run, Int_t sector, Double_t timeStamp, Int_t *entries)
{
  /// GetCE mean charge for 'sector'
  /// it timestamp==-1 return mean value

  AliTPCcalibDB::Instance()->SetRun(run);
  TGraph *gr=AliTPCcalibDB::Instance()->GetCErocQgraph(sector);
  if (!gr||sector<0||sector>71) {
    if (entries) *entries=0;
    return 0.;
  }
  Float_t val=0.;
  if (timeStamp==-1.){
    val=gr->GetMean(2);
  }else{
    for (Int_t ipoint=0;ipoint<gr->GetN();++ipoint){
      Double_t x,y;
      gr->GetPoint(ipoint,x,y);
      if (x<timeStamp) continue;
      val=y;
      break;
    }
  }
  return val;
}

Float_t AliTPCcalibDB::GetDCSSensorValue(AliDCSSensorArray *arr, Int_t timeStamp, const char * sensorName, Int_t sigDigits)
{
  /// Get Value for a DCS sensor 'sensorName', run 'run' at time 'timeStamp'

  Float_t val=0;
  const TString sensorNameString(sensorName);
  AliDCSSensor *sensor = arr->GetSensor(sensorNameString);
  const Int_t startTime = Int_t(sensor->GetStartTime());
  const Int_t endTime   = Int_t(sensor->GetEndTime());

  if (!sensor) return val;
  //use the dcs graph if possible
  TGraph *gr=sensor->GetGraph();
  if (gr){
    for (Int_t ipoint=0;ipoint<gr->GetN();++ipoint){
      Double_t x,y;
      gr->GetPoint(ipoint,x,y);
      const Int_t time=TMath::Nint(startTime+x*3600.); //time in graph is hours
      if (time<=timeStamp && timeStamp<=endTime) {
        val=y;
        continue;
      }
      break;
    }
    //if val is still 0, test if if the requested time is within 5min of the first/last
    //data point or start/end time of the sensor. If this is the case return the firs/last entry
    //the timestamps might not be syncronised for all calibration types, sometimes a 'pre'
    //and 'pos' period is requested. Especially to the HV this is not the case!
    //first point
    if (val==0 ){
      Double_t x,y;
      gr->GetPoint(0,x,y);
      const Int_t time=TMath::Min(TMath::Nint(startTime+x*3600.), startTime); //time in graph is hours
      const Int_t dtime=time-timeStamp;
      if ( (dtime>=0) && (dtime<5*60) ) val=y;
    }
    //last point
    if (val==0 ){
      Double_t x,y;
      gr->GetPoint(gr->GetN()-1,x,y);
      const Int_t time=TMath::Max(TMath::Nint(startTime+x*3600.), endTime); //time in graph is hours
      const Int_t dtime=timeStamp-time;
      if ( (dtime>=0) && (dtime<5*60) ) val=y;
    }
  } else {
    val=sensor->GetValue(timeStamp);
  }
  if (sigDigits>=0){
    val=(Float_t)TMath::Floor(val * TMath::Power(10., sigDigits) + .5) / TMath::Power(10., sigDigits);
  }
  return val;
}

Float_t AliTPCcalibDB::GetDCSSensorMeanValue(AliDCSSensorArray *arr, const char * sensorName, Int_t sigDigits)
{
  /// Get mean Value for a DCS sensor 'sensorName' during run 'run'

  Float_t val=0;
  const TString sensorNameString(sensorName);
  AliDCSSensor *sensor = arr->GetSensor(sensorNameString);
  if (!sensor) return val;

  //use dcs graph if it exists
  TGraph *gr=sensor->GetGraph();
  if (gr){
    val=gr->GetMean(2);
  } else {
    //if we don't have the dcs graph, try to get some meaningful information
    if (!sensor->GetFit()) return val;
    Int_t nKnots=sensor->GetFit()->GetKnots();
    Double_t tMid=(sensor->GetEndTime()-sensor->GetStartTime())/2.;
    for (Int_t iKnot=0;iKnot<nKnots;++iKnot){
      if (sensor->GetFit()->GetX()[iKnot]>tMid/3600.) break;
      val=(Float_t)sensor->GetFit()->GetY0()[iKnot];
    }
  }
  if (sigDigits>=0){
    // val/=10;
    val=(Float_t)TMath::Floor(val * TMath::Power(10., sigDigits) + .5) / TMath::Power(10., sigDigits);
    //    val*=10;
  }
  return val;
}

Int_t AliTPCcalibDB::GetMaskedChannelsFromCorrectionMaps(TBits maskedPads[72])
{
  // set bits for masked pads
  int nrowsMasked=0,npadsMasked=0,npadsTotMasked=0;
  if (AliTPCRecoParam::GetPrimaryDCACut()) {
    AliInfo("Special reconstruction for residuals extraction, masking is disabled");
    return 0;
  }
  AliCDBManager* man = AliCDBManager::Instance();
  AliMagF* fld = (AliMagF*)TGeoGlobalMagField::Instance()->GetField();
  if (!fld || !man->IsDefaultStorageSet() || man->GetRun()<0) {
    AliFatal("OCDB of B-field is not initialized");
  }
  const int run = GetRun();
  // pick the recoparam matching to field (as low or high flux)
  AliRecoParam::EventSpecie_t spec = fld->GetBeamType()==AliMagF::kBeamTypeAA ? AliRecoParam::kHighMult : AliRecoParam::kLowMult;
  AliTPCRecoParam* param = 0;
  int parID=0;
  while( (param=GetRecoParam(parID++)) ) { if (param->GetEventSpecie()&spec) break;}
  if (!param) AliFatal("Failed to extract recoparam");
  //
  if (!param->GetUseCorrectionMap()) {
    AliInfo("Residual correction maps are not used, no masking needed");
    return 0;
  }

  AliTPCTransform* transform = GetTransform();
  if (!transform) AliFatal("Failed to extract Transform");
  transform->SetCurrentRecoParam(param);

  // Load maps covering the run and build masks for disable full rows
  //
  // get reference map
  AliTPCChebCorr *mapRef = AliTPCTransform::LoadFieldDependendStaticCorrectionMap(kTRUE);
  //
  // get correction maps
  TObjArray* arrMaps = AliTPCTransform::LoadCorrectionMaps(kFALSE);
  arrMaps->SetOwner(kTRUE);
  if (!((AliTPCChebCorr*)arrMaps->UncheckedAt(0))->GetTimeDependent()) {
    // static maps are field-dependent
    AliTPCChebCorr* mapKeep = AliTPCTransform::LoadFieldDependendStaticCorrectionMap(kFALSE,arrMaps);
    arrMaps->Remove((TObject*)mapKeep);
    arrMaps->Delete();
    arrMaps->Add(mapKeep);  // leave only maps needed for this run if the object was default one
  }
  int nMaps = arrMaps->GetEntriesFast();

  // mask full rows disabled in at least one of the maps
  TBits maskedRows[72];
  for (int isec=72;isec--;) {  // flag masked full rows
    mapRef->GetNMaskedRows(isec,&maskedRows[isec]);
    for (int imap=nMaps;imap--;) ((AliTPCChebCorr*)arrMaps->UncheckedAt(imap))->GetNMaskedRows(isec,&maskedRows[isec]);
    nrowsMasked += maskedRows[isec].CountBits();
    //       printf("ROC%d masked %d rows\n",isec,maskedRows[isec].CountBits());
  }
  //
  // Now we need to mask individual pads where the assigned errors or distortions are too large
  // Since >1 map may cover the run, we query all off them in their central time stamps
  // 1) Make sure the maps are unique for this run
  AliGRPObject* grp = GetGRP(run);
  time_t tGRPmin = grp->GetTimeStart();
  time_t tGRPmax = grp->GetTimeEnd();
  time_t tCentGRP = (tGRPmin+tGRPmax)/2;  // center of the run according to grp
  //
  time_t mapsT[nMaps];
  for (int i=0;i<nMaps;i++) {
    mapsT[i] = ((AliTPCChebCorr*)arrMaps->At(i))->GetTimeStampCenter();
    if (mapsT[i]<tGRPmin || mapsT[i]>tCentGRP) mapsT[i] = tCentGRP;
  }
  //
  delete arrMaps; // we don't need anymore these maps, Transform will load according to time stamp
  delete mapRef;
  //
  AliTPCROC* roc = AliTPCROC::Instance();
  double testv[3]={0.};

  // determine ranges for time-bin to be assigned to test cluster
  transform->SetCurrentTimeStamp(mapsT[0]);
  const double tbinMin=0.,tbinMax=800.0;
  testv[0] = 62.; testv[1] = 0.; testv[2] = tbinMin; // at highest pad of IROC, min drift-time
  transform->Local2RotatedGlobal(0,testv);
  const double ztMin = testv[2];
  testv[0] = 62.; testv[1] = 0.; testv[2] = tbinMax; // at highest pad of IROC, max drift-time
  transform->Local2RotatedGlobal(0,testv);
  const double ztMax = testv[2];
  const double tzSlope = (tbinMax-tbinMin)/(ztMax-ztMin);
  //
  // for test only >>>>>>>>>>
  //  double *mxdist = (double*)param->GetBadPadMaxDistXYZ();
  //  double *mxerr  = (double*)param->GetBadPadMaxErrYZ();
  //  mxdist[0] = 7.; mxdist[1] = 7.; mxdist[2] = 7.;
  //  mxerr[0] = 2.; mxerr[1] = 2.;
  // for test only <<<<<<<<<<
  const double *padMaxDist = param->GetBadPadMaxDistXYZD();
  const double *padMaxErr  = param->GetBadClusMaxErrYZ();
  const double maxErrY2    = padMaxErr[0]>0 ? padMaxErr[0]*padMaxErr[0] : -1.;
  const double maxErrZ2    = padMaxErr[1]>0 ? padMaxErr[1]*padMaxErr[1] : -1.;
  //
  const float tgLabTest = 0.2f; // check distortions/errors for this pad for high pt track at tgLab=0.2
  AliTPCclusterMI clProbe;
  for (int imap=0;imap<nMaps;imap++) {
    AliInfoF("Querying maps at time %ld",mapsT[imap]);
    transform->SetCurrentTimeStamp(mapsT[imap]);
    //
    for (int isec=72;isec--;) {
      TBits& padStatus = maskedPads[isec];
      for (int irow=roc->GetNRows(isec);irow--;) {
        const int npads = roc->GetNPads(isec, irow), channel0 = roc->GetRowIndexes(isec)[irow];
        if (maskedRows[isec].TestBitNumber(irow)) { // full row is masked
          for (int ipad=npads;ipad--;) padStatus.SetBitNumber(channel0+ipad); // mask all pads of the row
          //printf("ROC%d masked %d pads of full row %d\n",isec,npads,irow);
          npadsTotMasked += npads;
          continue;
        }
        const double xRow = roc->GetPadRowRadii(isec,irow), zTest = xRow*tgLabTest;
        const double tbinTest = tbinMin+(zTest-ztMin)*tzSlope; // define tbin at which distortions/errors will be evaluated
        for (int ipad=npads;ipad--;) {
          testv[0] = irow;
          testv[1] = 0.5+ipad; // pad center
          testv[2] = tbinTest;
          transform->Transform(testv,&isec,0,0);
          const float* clCorr    = transform->GetLastMapCorrection();
          const float* clCorrRef = transform->GetLastMapCorrectionRef();
          Bool_t bad = kFALSE;

          // does the distortion exceed max.allowed?
          for (int dir=4;dir--;) if (padMaxDist[dir]>0. && TMath::Abs(clCorr[dir])>padMaxDist[dir]) {
            bad=kTRUE;
            // printf("ROC%2d row%2d pad%2d bad: distortion[%d]=%.2f exceeds threshold %.2f\n",
            //       isec,npads,irow,dir,clCorr[dir],padMaxDist[dir]);
            break;
          }

          if (!bad) { // check errors assigned to such cluster
            clProbe.SetDetector(isec);
            clProbe.SetX(testv[0]); // store coordinates
            clProbe.SetY(testv[1]);
            clProbe.SetZ(testv[2]);
            clProbe.SetDistortions(clCorr[0]-clCorrRef[0],clCorr[1]-clCorrRef[1],clCorr[2]-clCorrRef[2]);
            clProbe.SetDistortionDispersion(clCorr[3]); // store distortions wrt ref and dispersion (ref already subtracted)
            const double errY2 = transform->ErrY2Syst(&clProbe,testv[1]/testv[0]);
            const double errZ2 = transform->ErrZ2Syst(&clProbe,testv[2]/testv[0]);
            if ( (maxErrY2>0 && errY2>maxErrY2) || (maxErrZ2>0 && errZ2>maxErrZ2) ) {
              bad=kTRUE;
              //printf("ROC%2d row%2d pad%2d bad: Errors^2 for Y:%.3e or Z:%.3e exceeds threshold %.3e|%.3e\n",
              //     isec,npads,irow,errY2,errZ2,maxErrY2,maxErrZ2);
            }
          }
          //
          if (bad) {
            if (!padStatus.TestBitNumber(channel0+ipad)) npadsMasked++;
            padStatus.SetBitNumber(channel0+ipad);
          }
        } // loop over all pads of the row
      } // loop over all rows of the ROC
    } // loop over all ROCs
  } // loop over all maps
  //
  npadsTotMasked += npadsMasked;
  AliInfoF("masked %d pads (%d full padrows, %d individual pads)",
           npadsTotMasked,nrowsMasked,npadsMasked);
  return npadsTotMasked;
  //
}

Bool_t AliTPCcalibDB::IsDataTakingActive(time_t timeStamp)
{
  //
  // Check if the data taking is active.
  // This information ist based on the trigger scalers and calculated in UpdateChamberHighVoltageData() below.
  // in case there is no GRP object or no trigger scalers fGrRunState should be a NULL pointer
  //   if this is the case we assume by default that the data taking is active
  // NOTE: The logik changed. Before v5-06-03-79-gc804e5a we assumed by default the data taking is inactive
  //
  if (!fGrRunState) return kTRUE;
  Double_t time=Double_t(timeStamp);
  Int_t currentPoint=0;
  Bool_t currentVal=fGrRunState->GetY()[currentPoint]>0.5;
  Bool_t retVal=currentVal;
  Double_t currentTime=fGrRunState->GetX()[currentPoint];

  while (time>currentTime){
    retVal=currentVal;
    if (currentPoint==fGrRunState->GetN()) break;
    currentVal=fGrRunState->GetY()[currentPoint]>0.5;
    currentTime=fGrRunState->GetX()[currentPoint];
    ++currentPoint;
  }

  return retVal;
}

void AliTPCcalibDB::UpdateChamberHighVoltageData()
{
  /// set chamber high voltage data
  /// 1. Robust median (sampling the hv graphs over time)
  /// 2. Current nominal voltages (nominal voltage corrected for common HV offset)
  /// 3. Fraction of good HV values over time (deviation from robust median)
  /// 4. HV status, based on the above
  ///
  /// List of required OCDB Entries
  /// - GRP/GRP/Data
  /// - GRP/CTP/Scalers
  /// - TPC/Calib/HighVoltage
  /// - TPC/Calib/Parameters

  // reset active run state graph
  delete fGrRunState;
  fGrRunState=0x0;

  // start and end time of the run
  const Int_t run=GetRun();
  if (run<0) return;

  // if no valid run information - return
  AliGRPObject* grp = GetGRP(run);
  if (!grp) return;

  const Int_t startTimeGRP = grp->GetTimeStart();
  const Int_t stopTimeGRP  = grp->GetTimeEnd();

  //
  // In case we use a generated GRP we cannot make use of the start time and end time information
  // therefore we cannot calculate proper HV information and will skip this
  //
  if (startTimeGRP==0 && stopTimeGRP==0) {
    AliWarning("Using a generated GRP with 'GetTimeStart()' and 'GetTimeEnd()' == 0. Cannot calculate HV information.");
    return;
  }

  //
  // check active state by analysing the scalers
  //
  // initialise graph with active running
  const char* hltMode = NULL;
  hltMode = gSystem->Getenv("HLT_ONLINE_MODE");

  AliCDBEntry *entry = NULL;
  if (!hltMode) entry = GetCDBEntry("GRP/CTP/Scalers");
  if (entry) {
    // entry->SetOwner(kTRUE);
    AliTriggerRunScalers *sca = (AliTriggerRunScalers*)entry->GetObject();
    Int_t nchannels = sca->GetNumClasses(); // number of scaler channels (i.e. trigger classes)
    Int_t npoints = sca->GetScalersRecords()->GetEntries(); // number of samples

    // require at least two points from the scalers.
    if (npoints>1) {
      fGrRunState=new TGraph;
      fGrRunState->SetPoint(fGrRunState->GetN(),Double_t(startTimeGRP)-.001,0);
      fGrRunState->SetPoint(fGrRunState->GetN(),Double_t(startTimeGRP),1);
      ULong64_t lastSum=0;
      Double_t timeLast=Double_t(startTimeGRP);
      Bool_t active=kTRUE;
      for (int i=0; i<npoints; i++) {
        AliTriggerScalersRecord *rec = (AliTriggerScalersRecord *) sca->GetScalersRecord(i);
        Double_t time = ((AliTimeStamp*) rec->GetTimeStamp())->GetSeconds();
        // check if time is inside the grp times. For dummy scaler entries the time might be compatible with 0
        if ( time<startTimeGRP || time>stopTimeGRP ){
          AliWarning(Form("Time of scaler record %d: %.0f is outside the GRP times (%d, %d). Skipping this record.", i, time, startTimeGRP, stopTimeGRP));
          continue;
        }
        ULong64_t sum=0;
        for (int j=0; j<nchannels; j++) sum += ((AliTriggerScalers*) rec->GetTriggerScalers()->At(j))->GetL2CA();
        if (TMath::Abs(time-timeLast)<.001 && sum==lastSum ) continue;
        if (active && sum==lastSum){
          fGrRunState->SetPoint(fGrRunState->GetN(),timeLast-.01,1);
          fGrRunState->SetPoint(fGrRunState->GetN(),timeLast,0);
          active=kFALSE;
        } else if (!active && sum>lastSum ){
          fGrRunState->SetPoint(fGrRunState->GetN(),timeLast-.01,0);
          fGrRunState->SetPoint(fGrRunState->GetN(),timeLast,1);
          active=kTRUE;
        }
        lastSum=sum;
        timeLast=time;
      }
      fGrRunState->SetPoint(fGrRunState->GetN(),Double_t(stopTimeGRP),active);
      fGrRunState->SetPoint(fGrRunState->GetN(),Double_t(stopTimeGRP)+.001,0);
    } else {
      AliWarning("Only one entry found in the trigger scalers. Most probably this is a dummy entry. Scaler information will not be used!");
    }
  }


  // reset all values
  for (Int_t iROC=0;iROC<72;++iROC) {
    fChamberHVmedian[iROC]       = -1;
    fChamberHVgoodFraction[iROC] = 0.;
    fCurrentNominalVoltage[iROC] = -999.;
    fChamberHVStatus[iROC]       = kFALSE;
  }

  AliDCSSensorArray* voltageArray = GetVoltageSensors(run);
  if (!voltageArray) {
    AliError("Voltage Array missing. Cannot calculate HV information!");
    AliError(" -> Check OCDB entry: 'TPC/Calib/HighVoltage'");
    return;
  }

  // max HV diffs before a chamber is masked
  const Float_t maxVdiff      = fParam->GetMaxVoltageDeviation();
  const Float_t maxDipVoltage = fParam->GetMaxDipVoltage();
  const Float_t maxFracHVbad  = fParam->GetMaxFractionHVbad();

  const Int_t samplingPeriod=1;

  // array with sampled voltages
  const Int_t maxSamples=(stopTimeGRP-startTimeGRP)/samplingPeriod + 10*samplingPeriod;
  Float_t *vSampled = new Float_t[maxSamples];

  // deviation of the median from the nominal voltage
  Double_t chamberMedianDeviation[72]={0.};

  for (Int_t iROC=0; iROC<72; ++iROC){
    chamberMedianDeviation[iROC]=0.;
    TString sensorName="";
    Char_t sideName='A';
    if ((iROC/18)%2==1) sideName='C';
    if (iROC<36) sensorName=Form("TPC_ANODE_I_%c%02d_VMEAS",sideName,iROC%18);
    else        sensorName=Form("TPC_ANODE_O_%c%02d_0_VMEAS",sideName,iROC%18);

    AliDCSSensor *sensor = voltageArray->GetSensor(sensorName);

    fHVsensors[iROC]=sensor;
    if (!sensor) continue;

    Int_t nPointsSampled=0;

    TGraph *gr=sensor->GetGraph();
    if ( gr && gr->GetN()>0 ){
      //1. sample voltage over time
      //   get a robust median
      //   buffer sampled voltages

      // current sampling time
      Int_t time=startTimeGRP;

      // input graph sampling point
      const Int_t nGraph=gr->GetN();
      Int_t pointGraph=0;

      //initialise graph information
      Int_t timeGraph=stopTimeGRP;
      if (gr->GetN()>1) timeGraph=TMath::Nint(gr->GetX()[pointGraph+1]*3600+sensor->GetStartTime());
      Double_t sampledHV=gr->GetY()[pointGraph++];

      while (time<stopTimeGRP){
        while (timeGraph<=time && pointGraph+1<nGraph){
          timeGraph=TMath::Nint(gr->GetX()[pointGraph+1]*3600+sensor->GetStartTime());
          sampledHV=gr->GetY()[pointGraph++];
        }
        time+=samplingPeriod;
        if (!IsDataTakingActive(time-samplingPeriod)) continue;
        vSampled[nPointsSampled++]=sampledHV;
      }

      if (nPointsSampled<1) continue;

      fChamberHVmedian[iROC]=TMath::Median(nPointsSampled,vSampled);
      chamberMedianDeviation[iROC]=fChamberHVmedian[iROC]-fParam->GetNominalVoltage(iROC);

      //2. calculate good HV fraction
      Int_t ngood=0;
      for (Int_t ipoint=0; ipoint<nPointsSampled; ++ipoint) {
        if (TMath::Abs(vSampled[ipoint]-fChamberHVmedian[iROC])<maxDipVoltage) ++ngood;
      }

      fChamberHVgoodFraction[iROC]=Float_t(ngood)/Float_t(nPointsSampled);
    } else if (!gr && !sensor->GetFit() ){
      // This is an exception handling.
      // It was observed that for some rund in the 2010 data taking no HV info is available
      //    for some sectors. However they were active. So take care about this
      fChamberHVmedian[iROC]       = fParam->GetNominalVoltage(iROC);
      fChamberHVgoodFraction[iROC] = 1.;
      AliWarning(Form("ROC %d detected without HV Splines and HV graph. Will set median HV to nominal voltage",iROC));
    } else {
      AliError(Form("No Graph or graph without points found for HV sensor of ROC %d",iROC));
    }
  }

  delete [] vSampled;
  vSampled=0x0;

  // get median deviation from all chambers (detect e.g. -50V)
  const Double_t medianIROC=TMath::Median( 36, chamberMedianDeviation );
  const Double_t medianOROC=TMath::Median( 36, chamberMedianDeviation+36 );

  // Define current default voltages
  for (Int_t iROC=0;iROC<72/*AliTPCCalPad::kNsec*/;++iROC){
    const Float_t averageDeviation=(iROC<36)?medianIROC:medianOROC;
    fCurrentNominalVoltage[iROC]=fParam->GetNominalVoltage(iROC)+averageDeviation;
  }

  //
  // Check HV status
  //
  Int_t nbad=0;
  for (Int_t iROC=0;iROC<72/*AliTPCCalPad::kNsec*/;++iROC){
    fChamberHVStatus[iROC]=kTRUE;
    const Float_t averageDeviation=(iROC<36)?medianIROC:medianOROC;

    //a. Deviation of median from current nominal voltage
    //   allow larger than nominal voltages
    if (fCurrentNominalVoltage[iROC]-fChamberHVmedian[iROC] >  maxVdiff) {
      AliWarning(Form("Low voltage detected for ROC %2d",iROC));
      AliWarning(Form(" -> Based on nominal voltage: %.2f + averge deviation: %.2f = current nominal voltage: %.2f - meadian HV: %.2f > max diff: %.2f",
                      fParam->GetNominalVoltage(iROC), averageDeviation, fCurrentNominalVoltage[iROC], fChamberHVmedian[iROC],  maxVdiff));
      AliWarning(" -> Check consistent usage of HV information in the OCDB entry: 'TPC/Calib/HighVoltage'");
      AliWarning(" ->                 and the nominal voltages in the OCDB entry: 'TPC/Calib/Parameters'");
      fChamberHVStatus[iROC]=kFALSE;
    }

    //b. Fraction of bad hv values
    if ( 1.-fChamberHVgoodFraction[iROC] > maxFracHVbad ) {
      AliWarning(Form("Large fraction of low HV readings detected in ROC %2d: %.2f > %.2f",
                   iROC, 1.-fChamberHVgoodFraction[iROC], maxFracHVbad));
      AliWarning(Form(" -> Based on HV information from OCDB entry: 'TPC/Calib/HighVoltage' with median voltage: %.2f", fChamberHVmedian[iROC]));
      AliWarning(     " -> Check with experts if this chamber had HV problems in this run");
      fChamberHVStatus[iROC]=kFALSE;
    }

    if (!fChamberHVStatus[iROC]) {
      ++nbad;
    }
  }

  // check if all chamber are off
  if (nbad==72) {
    AliFatal("Something went wrong in the chamber HV status calculation. Check warning messages above. All chambers would be deactivated!");
  }
}

Float_t AliTPCcalibDB::GetChamberHighVoltage(Int_t run, Int_t sector, Int_t timeStamp, Int_t sigDigits, Bool_t current) {
  /// return the chamber HV for given run and time: 0-35 IROC, 36-72 OROC
  /// if timeStamp==-1 return mean value

  Float_t val=0;
  TString sensorName="";
  TTimeStamp stamp(timeStamp);
  AliDCSSensorArray* voltageArray = AliTPCcalibDB::Instance()->GetVoltageSensors(run);
  if (!voltageArray || (sector<0) || (sector>71)) return val;
  Char_t sideName='A';
  if ((sector/18)%2==1) sideName='C';
  if (sector<36){
    //IROC
    sensorName=Form("TPC_ANODE_I_%c%02d_VMEAS",sideName,sector%18);
  }else{
    //OROC
    sensorName=Form("TPC_ANODE_O_%c%02d_0_VMEAS",sideName,sector%18);
  }
  if (current){
    if (sector<36){
      //IROC
      sensorName=Form("TPC_ANODE_I_%c%02d_IMEAS",sideName,sector%18);
    }else{
      //OROC
      sensorName=Form("TPC_ANODE_O_%c%02d_0_IMEAS",sideName,sector%18);
    }

  }
  if (timeStamp==-1){
    val=AliTPCcalibDB::GetDCSSensorMeanValue(voltageArray, sensorName.Data(),sigDigits);
  } else {
    val=AliTPCcalibDB::GetDCSSensorValue(voltageArray, timeStamp, sensorName.Data(),sigDigits);
  }
  return val;
}
Float_t AliTPCcalibDB::GetSkirtVoltage(Int_t run, Int_t sector, Int_t timeStamp, Int_t sigDigits)
{
  /// Get the skirt voltage for 'run' at 'timeStamp' and 'sector': 0-35 IROC, 36-72 OROC
  /// type corresponds to the following: 0 - IROC A-Side; 1 - IROC C-Side; 2 - OROC A-Side; 3 - OROC C-Side
  /// if timeStamp==-1 return the mean value for the run

  Float_t val=0;
  TString sensorName="";
  TTimeStamp stamp(timeStamp);
  AliDCSSensorArray* voltageArray = AliTPCcalibDB::Instance()->GetVoltageSensors(run);
  if (!voltageArray || (sector<0) || (sector>71)) return val;
  Char_t sideName='A';
  if ((sector/18)%2==1) sideName='C';
  sensorName=Form("TPC_SKIRT_%c_VMEAS",sideName);
  if (timeStamp==-1){
    val=AliTPCcalibDB::GetDCSSensorMeanValue(voltageArray, sensorName.Data(),sigDigits);
  } else {
    val=AliTPCcalibDB::GetDCSSensorValue(voltageArray, timeStamp, sensorName.Data(),sigDigits);
  }
  return val;
}

Float_t AliTPCcalibDB::GetCoverVoltage(Int_t run, Int_t sector, Int_t timeStamp, Int_t sigDigits)
{
  /// Get the cover voltage for run 'run' at time 'timeStamp'
  /// type corresponds to the following: 0 - IROC A-Side; 1 - IROC C-Side; 2 - OROC A-Side; 3 - OROC C-Side
  /// if timeStamp==-1 return the mean value for the run

  Float_t val=0;
  TString sensorName="";
  TTimeStamp stamp(timeStamp);
  AliDCSSensorArray* voltageArray = AliTPCcalibDB::Instance()->GetVoltageSensors(run);
  if (!voltageArray || (sector<0) || (sector>71)) return val;
  Char_t sideName='A';
  if ((sector/18)%2==1) sideName='C';
  if (sector<36){
    //IROC
    sensorName=Form("TPC_COVER_I_%c_VMEAS",sideName);
  }else{
    //OROC
    sensorName=Form("TPC_COVER_O_%c_VMEAS",sideName);
  }
  if (timeStamp==-1){
    val=AliTPCcalibDB::GetDCSSensorMeanValue(voltageArray, sensorName.Data(),sigDigits);
  } else {
    val=AliTPCcalibDB::GetDCSSensorValue(voltageArray, timeStamp, sensorName.Data(),sigDigits);
  }
  return val;
}

Float_t AliTPCcalibDB::GetGGoffsetVoltage(Int_t run, Int_t sector, Int_t timeStamp, Int_t sigDigits)
{
  /// Get the GG offset voltage for run 'run' at time 'timeStamp'
  /// type corresponds to the following: 0 - IROC A-Side; 1 - IROC C-Side; 2 - OROC A-Side; 3 - OROC C-Side
  /// if timeStamp==-1 return the mean value for the run

  Float_t val=0;
  TString sensorName="";
  TTimeStamp stamp(timeStamp);
  AliDCSSensorArray* voltageArray = AliTPCcalibDB::Instance()->GetVoltageSensors(run);
  if (!voltageArray || (sector<0) || (sector>71)) return val;
  Char_t sideName='A';
  if ((sector/18)%2==1) sideName='C';
  if (sector<36){
    //IROC
    sensorName=Form("TPC_GATE_I_%c_OFF_VMEAS",sideName);
  }else{
    //OROC
    sensorName=Form("TPC_GATE_O_%c_OFF_VMEAS",sideName);
  }
  if (timeStamp==-1){
    val=AliTPCcalibDB::GetDCSSensorMeanValue(voltageArray, sensorName.Data(),sigDigits);
  } else {
    val=AliTPCcalibDB::GetDCSSensorValue(voltageArray, timeStamp, sensorName.Data(),sigDigits);
  }
  return val;
}

Float_t AliTPCcalibDB::GetGGnegVoltage(Int_t run, Int_t sector, Int_t timeStamp, Int_t sigDigits)
{
  /// Get the GG offset voltage for run 'run' at time 'timeStamp'
  /// type corresponds to the following: 0 - IROC A-Side; 1 - IROC C-Side; 2 - OROC A-Side; 3 - OROC C-Side
  /// if timeStamp==-1 return the mean value for the run

  Float_t val=0;
  TString sensorName="";
  TTimeStamp stamp(timeStamp);
  AliDCSSensorArray* voltageArray = AliTPCcalibDB::Instance()->GetVoltageSensors(run);
  if (!voltageArray || (sector<0) || (sector>71)) return val;
  Char_t sideName='A';
  if ((sector/18)%2==1) sideName='C';
  if (sector<36){
    //IROC
    sensorName=Form("TPC_GATE_I_%c_NEG_VMEAS",sideName);
  }else{
    //OROC
    sensorName=Form("TPC_GATE_O_%c_NEG_VMEAS",sideName);
  }
  if (timeStamp==-1){
    val=AliTPCcalibDB::GetDCSSensorMeanValue(voltageArray, sensorName.Data(),sigDigits);
  } else {
    val=AliTPCcalibDB::GetDCSSensorValue(voltageArray, timeStamp, sensorName.Data(),sigDigits);
  }
  return val;
}

Float_t AliTPCcalibDB::GetGGposVoltage(Int_t run, Int_t sector, Int_t timeStamp, Int_t sigDigits)
{
  /// Get the GG offset voltage for run 'run' at time 'timeStamp'
  /// type corresponds to the following: 0 - IROC A-Side; 1 - IROC C-Side; 2 - OROC A-Side; 3 - OROC C-Side
  /// if timeStamp==-1 return the mean value for the run

  Float_t val=0;
  TString sensorName="";
  TTimeStamp stamp(timeStamp);
  AliDCSSensorArray* voltageArray = AliTPCcalibDB::Instance()->GetVoltageSensors(run);
  if (!voltageArray || (sector<0) || (sector>71)) return val;
  Char_t sideName='A';
  if ((sector/18)%2==1) sideName='C';
  if (sector<36){
    //IROC
    sensorName=Form("TPC_GATE_I_%c_POS_VMEAS",sideName);
  }else{
    //OROC
    sensorName=Form("TPC_GATE_O_%c_POS_VMEAS",sideName);
  }
  if (timeStamp==-1){
    val=AliTPCcalibDB::GetDCSSensorMeanValue(voltageArray, sensorName.Data(),sigDigits);
  } else {
    val=AliTPCcalibDB::GetDCSSensorValue(voltageArray, timeStamp, sensorName.Data(),sigDigits);
  }
  return val;
}

Float_t AliTPCcalibDB::GetPressure(Int_t timeStamp, Int_t run, Int_t type){
  /// GetPressure for given time stamp and runt

  TTimeStamp stamp(timeStamp);
  AliDCSSensor * sensor = Instance()->GetPressureSensor(run,type);
  if (!sensor) return 0;
  return sensor->GetValue(stamp);
}

Float_t AliTPCcalibDB::GetL3Current(Int_t run, Int_t statType){
  /// return L3 current
  /// stat type is: AliGRPObject::Stats: kMean = 0, kTruncMean = 1, kMedian = 2, kSDMean = 3, kSDMedian = 4

  Float_t current=-1;
  AliGRPObject *grp=AliTPCcalibDB::GetGRP(run);
  if (grp) current=grp->GetL3Current((AliGRPObject::Stats)statType);
  return current;
}

Float_t AliTPCcalibDB::GetBz(Int_t run){
  /// calculate BZ in T from L3 current

  Float_t bz=-1;
  Float_t current=AliTPCcalibDB::GetL3Current(run);
  if (current>-1) bz=5*current/30000.*.1;
  return bz;
}

Char_t  AliTPCcalibDB::GetL3Polarity(Int_t run) {
  /// get l3 polarity from GRP

  Char_t pol=-100;
  AliGRPObject *grp=AliTPCcalibDB::GetGRP(run);
  if (grp) pol=grp->GetL3Polarity();
  return pol;
}

TString AliTPCcalibDB::GetRunType(Int_t run){
  /// return run type from grp

//   TString type("UNKNOWN");
  AliGRPObject *grp=AliTPCcalibDB::GetGRP(run);
  if (grp) return grp->GetRunType();
  return "UNKNOWN";
}

Float_t AliTPCcalibDB::GetValueGoofie(Int_t timeStamp, Int_t run, Int_t type){
  /// GetPressure for given time stamp and runt

  TTimeStamp stamp(timeStamp);
  AliDCSSensorArray* goofieArray = AliTPCcalibDB::Instance()->GetGoofieSensors(run);
  if (!goofieArray) return 0;
  AliDCSSensor *sensor = goofieArray->GetSensor(type);
  return sensor->GetValue(stamp);
}






Bool_t  AliTPCcalibDB::GetTemperatureFit(Int_t timeStamp, Int_t run, Int_t side,TVectorD& fit){
  /// GetTmeparature fit at parameter for given time stamp

  TTimeStamp tstamp(timeStamp);
  AliTPCSensorTempArray* tempArray  = Instance()->GetTemperatureSensor(run);
  if (! tempArray) return kFALSE;
  AliTPCTempMap * tempMap = new AliTPCTempMap(tempArray);
  TLinearFitter * fitter = tempMap->GetLinearFitter(3,side,tstamp);
  if (fitter){
    fitter->Eval();
    fitter->GetParameters(fit);
  }
  delete fitter;
  delete tempMap;
  if (!fitter) return kFALSE;
  return kTRUE;
}

Float_t AliTPCcalibDB::GetTemperature(Int_t timeStamp, Int_t run, Int_t side){
  /// Get mean temperature

  TVectorD vec(5);
  if (side==0) {
    GetTemperatureFit(timeStamp,run,0,vec);
    return vec[0];
  }
  if (side==1){
    GetTemperatureFit(timeStamp,run,0,vec);
    return vec[0];
  }
  return 0;
}


Double_t AliTPCcalibDB::GetPTRelative(UInt_t timeSec, Int_t run, Int_t side){
  /// Get relative P/T
  /// time - absolute time
  /// run  - run number
  /// side - 0 - A side   1-C side

  AliTPCCalibVdrift * vdrift =  Instance()->GetVdrift(run);
  if (!vdrift) return 0;
  return vdrift->GetPTRelative(timeSec,side);
}

AliGRPObject * AliTPCcalibDB::MakeGRPObjectFromMap(TMap *map){
  /// Function to covert old GRP run information from TMap to GRPObject
  ///
  ///  TMap * map = AliTPCcalibDB::GetGRPMap(52406);

  if (!map) return 0;
  AliDCSSensor * sensor = 0;
  TObject *osensor=0;
  osensor = ((*map)("fP2Pressure"));
  sensor  =dynamic_cast<AliDCSSensor *>(osensor);
  //
  if (!sensor) return 0;
  //
  AliDCSSensor * sensor2 = new AliDCSSensor(*sensor);
  osensor = ((*map)("fCavernPressure"));
  TGraph * gr = new TGraph(2);
  gr->GetX()[0]= -100000.;
  gr->GetX()[1]= 1000000.;
  gr->GetY()[0]= atof(osensor->GetName());
  gr->GetY()[1]= atof(osensor->GetName());
  sensor2->SetGraph(gr);
  sensor2->SetFit(0);


  AliGRPObject *grpRun = new AliGRPObject;
  grpRun->ReadValuesFromMap(map);
  grpRun->SetCavernAtmosPressure(sensor2);
  grpRun->SetCavernAtmosPressure(sensor2);
  grpRun->SetSurfaceAtmosPressure(sensor);
  return grpRun;
}

Bool_t AliTPCcalibDB::CreateGUITree(Int_t run, const char* filename)
{
  /// Create a gui tree for run number 'run'

  if (!AliCDBManager::Instance()->GetDefaultStorage()){
    AliLog::Message(AliLog::kError, "Default Storage not set. Cannot create Calibration Tree!",
                    MODULENAME(), "AliTPCcalibDB", FUNCTIONNAME(), __FILE__, __LINE__);
    return kFALSE;
  }
  //db instance
  AliTPCcalibDB *db=AliTPCcalibDB::Instance();
  // retrieve cal pad objects
  db->SetRun(run);
  db->CreateGUITree(filename);
  return kTRUE;
}

Bool_t AliTPCcalibDB::CreateGUITree(const char* filename){
  ///

  if (!AliCDBManager::Instance()->GetDefaultStorage()){
    AliError("Default Storage not set. Cannot create calibration Tree!");
    return kFALSE;
  }
  UpdateNonRec();  // load all infromation now

  AliTPCPreprocessorOnline prep;
  if (GetActiveChannelMap()) prep.AddComponent(new AliTPCCalPad(*GetActiveChannelMap()));

  // gain map
  if (GetDedxGainFactor()) prep.AddComponent(new AliTPCCalPad(*GetDedxGainFactor()));
  //noise and pedestals
  if (GetPedestals()) prep.AddComponent(new AliTPCCalPad(*(GetPedestals())));
  if (GetPadNoise() ) prep.AddComponent(new AliTPCCalPad(*(GetPadNoise())));
  //pulser data
  if (GetPulserTmean()) prep.AddComponent(new AliTPCCalPad(*(GetPulserTmean())));
  if (GetPulserTrms() ) prep.AddComponent(new AliTPCCalPad(*(GetPulserTrms())));
  if (GetPulserQmean()) prep.AddComponent(new AliTPCCalPad(*(GetPulserQmean())));
  //CE data
  if (GetCETmean()) prep.AddComponent(new AliTPCCalPad(*(GetCETmean())));
  if (GetCETrms() ) prep.AddComponent(new AliTPCCalPad(*(GetCETrms())));
  if (GetCEQmean()) prep.AddComponent(new AliTPCCalPad(*(GetCEQmean())));
  //Altro data
  if (GetALTROAcqStart() ) prep.AddComponent(new AliTPCCalPad(*(GetALTROAcqStart() )));
  if (GetALTROZsThr()    ) prep.AddComponent(new AliTPCCalPad(*(GetALTROZsThr()    )));
  if (GetALTROFPED()     ) prep.AddComponent(new AliTPCCalPad(*(GetALTROFPED()     )));
  if (GetALTROAcqStop()  ) prep.AddComponent(new AliTPCCalPad(*(GetALTROAcqStop()  )));
  if (GetALTROMasked()   ) prep.AddComponent(new AliTPCCalPad(*(GetALTROMasked()   )));
  //QA
  AliTPCdataQA *dataQA=GetDataQA();
  if (dataQA) {
    if (dataQA->GetNLocalMaxima())
      prep.AddComponent(new AliTPCCalPad(*(dataQA->GetNLocalMaxima())));
    if (dataQA->GetMaxCharge())
      prep.AddComponent(new AliTPCCalPad(*(dataQA->GetMaxCharge())));
    if (dataQA->GetMeanCharge())
      prep.AddComponent(new AliTPCCalPad(*(dataQA->GetMeanCharge())));
    if (dataQA->GetNoThreshold())
      prep.AddComponent(new AliTPCCalPad(*(dataQA->GetNoThreshold())));
    if (dataQA->GetNTimeBins())
      prep.AddComponent(new AliTPCCalPad(*(dataQA->GetNTimeBins())));
    if (dataQA->GetNPads())
      prep.AddComponent(new AliTPCCalPad(*(dataQA->GetNPads())));
    if (dataQA->GetTimePosition())
      prep.AddComponent(new AliTPCCalPad(*(dataQA->GetTimePosition())));
  }

  //
  TString file(filename);
  if (file.IsNull()) file=Form("guiTreeRun_%i.root",fRun);
  prep.DumpToFile(file.Data());
  return kTRUE;
}

Bool_t AliTPCcalibDB::CreateRefFile(Int_t run, const char* filename)
{
  /// Create a gui tree for run number 'run'

  if (!AliCDBManager::Instance()->GetDefaultStorage()){
    AliLog::Message(AliLog::kError, "Default Storage not set. Cannot create Calibration Tree!",
                    MODULENAME(), "AliTPCcalibDB", FUNCTIONNAME(), __FILE__, __LINE__);
    return kFALSE;
  }
  TString file(filename);
  if (file.IsNull()) file=Form("RefCalPads_%d.root",run);
  TDirectory *currDir=gDirectory;
  //db instance
  AliTPCcalibDB *db=AliTPCcalibDB::Instance();
  // retrieve cal pad objects
  db->SetRun(run);
  //open file
  TFile f(file.Data(),"recreate");
  //noise and pedestals
  db->GetPedestals()->Write("Pedestals");
  db->GetPadNoise()->Write("PadNoise");
  //pulser data
  db->GetPulserTmean()->Write("PulserTmean");
  db->GetPulserTrms()->Write("PulserTrms");
  db->GetPulserQmean()->Write("PulserQmean");
  //CE data
  db->GetCETmean()->Write("CETmean");
  db->GetCETrms()->Write("CETrms");
  db->GetCEQmean()->Write("CEQmean");
  //Altro data
  db->GetALTROAcqStart() ->Write("ALTROAcqStart");
  db->GetALTROZsThr()    ->Write("ALTROZsThr");
  db->GetALTROFPED()     ->Write("ALTROFPED");
  db->GetALTROAcqStop()  ->Write("ALTROAcqStop");
  db->GetALTROMasked()   ->Write("ALTROMasked");
  //
  f.Close();
  currDir->cd();
  return kTRUE;
}



Double_t AliTPCcalibDB::GetVDriftCorrectionTime(Int_t timeStamp, Int_t run, Int_t /*side*/, Int_t mode){
  /// Get time dependent drift velocity correction
  /// multiplication factor        vd = vdnom *(1+vdriftcorr)
  /// Arguments:
  /// mode determines the algorith how to combine the Laser Track, LaserCE and physics tracks
  /// timestamp - timestamp
  /// run       - run number
  /// side      - the drift velocity per side (possible for laser and CE)
  ///
  /// Notice - Extrapolation outside of calibration range  - using constant function

  Double_t result=0;
  // mode 1  automatic mode - according to the distance to the valid calibration
  //                        -
  Double_t deltaP=0,  driftP=0,      wP  = 0.;
  Double_t deltaITS=0,driftITS=0,    wITS= 0.;
  Double_t deltaLT=0, driftLT=0,     wLT = 0.;
  Double_t deltaCE=0, driftCE=0,     wCE = 0.;
  driftP  = fDButil->GetVDriftTPC(deltaP,run,timeStamp);
  driftITS= fDButil->GetVDriftTPCITS(deltaITS,run,timeStamp);
  driftCE = fDButil->GetVDriftTPCCE(deltaCE, run,timeStamp,36000,2);
  driftLT = fDButil->GetVDriftTPCLaserTracks(deltaLT,run,timeStamp,36000,2);
  deltaITS = TMath::Abs(deltaITS);
  deltaP   = TMath::Abs(deltaP);
  deltaLT  = TMath::Abs(deltaLT);
  deltaCE  = TMath::Abs(deltaCE);
  if (mode==1) {
    const Double_t kEpsilon=0.00000000001;
    const Double_t kdeltaT=360.; // 10 minutes
    if(TMath::Abs(deltaITS) < 12*kdeltaT) {
      result = driftITS;
    } else {
    wITS  = 64.*kdeltaT/(deltaITS +kdeltaT);
    wLT   = 16.*kdeltaT/(deltaLT  +kdeltaT);
    wP    = 0. *kdeltaT/(deltaP   +kdeltaT);
    wCE   = 1. *kdeltaT/(deltaCE  +kdeltaT);
    //
    //
    if (TMath::Abs(driftP)<kEpsilon)  wP=0;  // invalid calibration
    if (TMath::Abs(driftITS)<kEpsilon)wITS=0;  // invalid calibration
    if (TMath::Abs(driftLT)<kEpsilon) wLT=0;  // invalid calibration
    if (TMath::Abs(driftCE)<kEpsilon) wCE=0;  // invalid calibration
    if (wP+wITS+wLT+wCE<kEpsilon) return 0;
    result = (driftP*wP+driftITS*wITS+driftLT*wLT+driftCE*wCE)/(wP+wITS+wLT+wCE);
   }


  }

  return result;
}

Double_t AliTPCcalibDB::GetTime0CorrectionTime(Int_t timeStamp, Int_t run, Int_t /*side*/, Int_t mode){
  /// Get time dependent time 0 (trigger delay in cm) correction
  /// additive correction        time0 = time0+ GetTime0CorrectionTime
  /// Value etracted combining the vdrift correction using laser tracks and CE and the physics track matchin
  /// Arguments:
  /// mode determines the algorith how to combine the Laser Track and physics tracks
  /// timestamp - timestamp
  /// run       - run number
  /// side      - the drift velocity per side (possible for laser and CE)
  ///
  /// Notice - Extrapolation outside of calibration range  - using constant function

  Double_t result=0;
  if (mode==2) {
    // TPC-TPC mode
    result=fDButil->GetTriggerOffsetTPC(run,timeStamp);
    result  *=fParam->GetZLength();
  }
  if (mode==1){
    // TPC-ITS mode
    Double_t dist=0;
    result= -fDButil->GetTime0TPCITS(dist, run, timeStamp)*fParam->GetDriftV()/1000000.;
  }
  return result;

}




Double_t AliTPCcalibDB::GetVDriftCorrectionGy(Int_t timeStamp, Int_t run, Int_t side, Int_t /*mode*/){
  /// Get global y correction drift velocity correction factor
  /// additive factor        vd = vdnom*(1+GetVDriftCorrectionGy *gy)
  /// Value etracted combining the vdrift correction using laser tracks and CE or TPC-ITS
  /// Arguments:
  /// mode determines the algorith how to combine the Laser Track, LaserCE or TPC-ITS
  /// timestamp - timestamp
  /// run       - run number
  /// side      - the drift velocity gy correction per side (CE and Laser tracks)
  ///
  /// Notice - Extrapolation outside of calibration range  - using constant function

  if (run<=0 && fTransform) run = fTransform->GetCurrentRunNumber();
  UpdateRunInformations(run,kFALSE);
  TObjArray *array =AliTPCcalibDB::Instance()->GetTimeVdriftSplineRun(run);
  if (!array) return 0;
  Double_t result=0;

  // use TPC-ITS if present
  TGraphErrors *gr= (TGraphErrors*)array->FindObject("ALIGN_ITSB_TPC_VDGY");
  if (!gr) gr = (TGraphErrors*)array->FindObject("ALIGN_TOFB_TPC_VDGY");
  if(gr) {
    result = AliTPCcalibDButil::EvalGraphConst(gr,timeStamp);

    // transform from [(cm/mus)/ m] to [1/cm]
    result /= (fParam->GetDriftV()/1000000.);
    result /= 100.;

    //printf("result %e \n", result);
    return result;
  }

  // use laser if ITS-TPC not present
  TGraphErrors *laserA= (TGraphErrors*)array->FindObject("GRAPH_MEAN_GLOBALYGRADIENT_LASER_ALL_A");
  TGraphErrors *laserC= (TGraphErrors*)array->FindObject("GRAPH_MEAN_GLOBALYGRADIENT_LASER_ALL_C");

  if (laserA && laserC){
   result= (laserA->Eval(timeStamp)+laserC->Eval(timeStamp))*0.5;
  }
  if (laserA && side==0){
    result = (laserA->Eval(timeStamp));
  }
  if (laserC &&side==1){
    result = (laserC->Eval(timeStamp));
  }
  //printf("laser result %e \n", -result/250.);

  return -result/250.; //normalized before
}


Double_t AliTPCcalibDB::GetVDriftCorrectionDeltaZ(Int_t /*timeStamp*/, Int_t run, Int_t /*side*/, Int_t /*mode*/){
  /// Get deltaZ run/by/run  correction - as fitted together with drift velocity
  /// Value extracted  form the TPC-ITS, mean value is used

  // Arguments:
  // mode determines the algorith how to combine the Laser Track, LaserCE or TPC-ITS
  // timestamp - not used
  // run       - run number
  // side      - common for boith sides
  //
  if (run<=0 && fTransform) run = fTransform->GetCurrentRunNumber();
  UpdateRunInformations(run,kFALSE);
  TObjArray *array =AliTPCcalibDB::Instance()->GetTimeVdriftSplineRun(run);
  if (!array) return 0;
  Double_t result=0;

  // use TPC-ITS if present
  TGraphErrors *gr= (TGraphErrors*)array->FindObject("ALIGN_ITSB_TPC_DELTAZ");
  if(gr) {
    result = TMath::Mean(gr->GetN(), gr->GetY());
  }
  return result;
}




AliTPCCalPad* AliTPCcalibDB::MakeDeadMap(Double_t notInMap, const char* nameMappingFile) {
///   Read list of active DDLs from OCDB entry
///   Generate and return AliTPCCalPad containing 1 for all pads in active DDLs,
///   0 for all pads in non-active DDLs.
///   For DDLs with missing status information (no DCS input point to Shuttle),
///     the value of the AliTPCCalPad entry is determined by the parameter
///     notInMap (default value 1)

  char chinfo[1000];

  TFile *fileMapping = new TFile(nameMappingFile, "read");
  AliTPCmapper *mapping = (AliTPCmapper*) fileMapping->Get("tpcMapping");
  if (!mapping) {
    snprintf(chinfo,1000,"Failed to get mapping object from %s.  ...\n", nameMappingFile);
    AliError (chinfo);
    return 0;
  }

  AliTPCCalPad *deadMap = new AliTPCCalPad("deadMap","deadMap");
  if (!deadMap) {
     AliError("Failed to allocate dead map AliTPCCalPad");
     return 0;
  }

  /// get list of active DDLs from OCDB entry
  Int_t idDDL=0;
  if (!fALTROConfigData ) {
     AliError("No ALTRO config OCDB entry available");
     return 0;
  }
  TMap *activeDDL = (TMap*)fALTROConfigData->FindObject("DDLArray");
  TObjString *ddlArray=0;
  if (activeDDL) {
    ddlArray = (TObjString*)activeDDL->GetValue("DDLArray");
    if (!ddlArray) {
      AliError("Empty list of active DDLs in OCDB entry");
      return 0;
    }
  } else {
    AliError("List of active DDLs not available in OCDB entry");
    return 0;
  }
  TString arrDDL=ddlArray->GetString();
  Int_t offset = mapping->GetTpcDdlOffset();
  Double_t active;
  for (Int_t i=0; i<mapping->GetNumDdl(); i++) {
    idDDL= i+offset;
    if (idDDL<0) continue;
    Int_t patch = mapping->GetPatchFromEquipmentID(idDDL);
    if (patch<0) continue;
    Int_t roc=mapping->GetRocFromEquipmentID(idDDL);
    if (roc<0) continue;
    AliTPCCalROC *calRoc=deadMap->GetCalROC(roc);
    if (calRoc) {
     for ( Int_t branch = 0; branch < 2; branch++ ) {
      for ( Int_t fec = 0; fec < mapping->GetNfec(patch, branch); fec++ ) {
        for ( Int_t altro = 0; altro < 8; altro++ ) {
         for ( Int_t channel = 0; channel < 16; channel++ ) {
           Int_t hwadd     = mapping->CodeHWAddress(branch, fec, altro, channel);
           Int_t row       = mapping->GetPadRow(patch, hwadd);        // row in a ROC (IROC or OROC)
//              Int_t globalrow = mapping.GetGlobalPadRow(patch, hwadd);  // row in full sector (IROC plus OROC)
           Int_t pad       = mapping->GetPad(patch, hwadd);
           if (!TString(arrDDL[i]).IsDigit()) {
	      active = notInMap;
           } else {
              active=TString(arrDDL[i]).Atof();
	   }
           calRoc->SetValue(row,pad,active);
         } // end channel for loop
        } // end altro for loop
      } // end fec for loop
     } // end branch for loop
    } // valid calROC
   } // end loop on active DDLs
   return deadMap;
}



AliTPCCorrection * AliTPCcalibDB::GetTPCComposedCorrection(Float_t field) const{
  /// GetComposed correction for given field setting
  /// If not specific correction for field used return correction for all field
  ///        - Complication needed to gaurantee OCDB back compatibility
  ///        - Not neeeded for the new space point correction

  if (!fComposedCorrectionArray) return 0;
  if (field>0.1 && fComposedCorrectionArray->At(1)) {
    return (AliTPCCorrection *)fComposedCorrectionArray->At(1);
  }
  if (field<-0.1 &&fComposedCorrectionArray->At(2)) {
    return (AliTPCCorrection *)fComposedCorrectionArray->At(2);
  }
  return (AliTPCCorrection *)fComposedCorrectionArray->At(0);

}


AliTPCCorrection * AliTPCcalibDB::GetTPCComposedCorrectionDelta() const{
  /// GetComposedCorrection delta
  /// Delta is time dependent - taken form the CalibTime OCDB entry

  if (!fComposedCorrectionArray) return 0;
  if (fRun<0) return 0;
  if (fDriftCorrectionArray.GetValue(Form("%i",fRun))==0) return 0;
  if (fComposedCorrectionArray->GetEntriesFast()<=4) {
    fComposedCorrectionArray->Expand(5);
    TObjArray * timeArray =(TObjArray*)(fDriftCorrectionArray.GetValue(Form("%i",fRun)));
     AliTPCCorrection * correctionTime = (AliTPCCorrection *)timeArray->FindObject("FitCorrectionTime");
     if (correctionTime){
       correctionTime->Init();
       fComposedCorrectionArray->AddAt(correctionTime,4); //add time dependent c
     }
  }
  return (AliTPCCorrection *)fComposedCorrectionArray->At(4);  //
}

Double_t AliTPCcalibDB::GetGainCorrectionHVandPT(Int_t timeStamp, Int_t run, Int_t sector, Int_t deltaCache, Int_t mode){
  /// Correction for  changes of gain caused by change of the HV and by relative change of the gas density
  /// Function is slow some kind of caching needed
  /// Cache implemented using the static TVectorD
  ///
  /// Input paremeters:
  ///  deltaCache - maximal time differnce above which the cache is recaclulated
  ///  mode       - mode==0 by default return combined correction
  ///                       actual HV and Pt correction has to be present in the run calibration otherwise it is ignored.
  ///                       (retrun value differnt than 1 only in case calibration present in the OCDB entry CalibTimeGain
  ///               mode==1 return combined correction ( important for calibration pass)
  ///                       (in case thereis  no calibration in  CalibTimeGain, default value from the AliTPCParam (Parameters) is used
  ///                       this mode is used in the CPass0
  ///               mode==2 return HV correction
  ///               mode==3 return P/T correction
  ///  Usage in the simulation/reconstruction
  ///  MC:     Qcorr  = Qorig*GetGainCorrectionHVandPT   ( in AliTPC.cxx )
  ///  Rec:    dEdx   = dEdx/GetGainCorrectionHVandPT    ( in aliTPCseed.cxx )

  static Float_t gGainCorrection[72];
  static Float_t gGainCorrectionPT[72];
  static Float_t gGainCorrectionHV[72];
  static Int_t    gTimeStamp=-99999999;
  static Bool_t   hasTimeDependent=kFALSE;
  if ( TMath::Abs(timeStamp-gTimeStamp)> deltaCache){
    //
    TGraphErrors * graphGHV = 0;
    TGraphErrors * graphGPT = 0;
    TObjArray *timeGainSplines = GetTimeGainSplinesRun(run);
    if (timeGainSplines){
      graphGHV  = (TGraphErrors*) timeGainSplines->FindObject("GainSlopesHV");
      graphGPT  = (TGraphErrors*) timeGainSplines->FindObject("GainSlopesPT");
      if (graphGHV) hasTimeDependent=kTRUE;
    }
    if (!graphGHV) graphGHV = fParam->GetGainSlopesHV();
    if (!graphGPT) graphGPT = fParam->GetGainSlopesPT();
    //
    for (Int_t isec=0; isec<72; isec++){
      Double_t HV= GetChamberHighVoltage(run,isec, timeStamp);
      if (HV<=0){ // check if the HV was available
        HV=GetChamberCurrentNominalHighVoltage(isec);
        AliWarningF("Could not get proper HV for run,sec,time (%d, %2d, %d), using current nominal voltage: %.2f", run, isec, timeStamp, HV);
//         AliDCSSensor* sensor = GetChamberHVSensor(isec);
//         if (sensor && sensor->GetGraph()==NULL &&  sensor->GetFit()==NULL){
//           HV=fParam->GetNominalVoltage(isec);
//         }
      }
      Double_t deltaHV= HV - fParam->GetNominalVoltage(isec);
      Double_t deltaGHV=0;
      Double_t deltaGPT=0;
      if (graphGHV) deltaGHV = graphGHV->GetY()[isec]*deltaHV;
      if (graphGPT) deltaGPT = graphGPT->GetY()[isec]*GetPTRelative(timeStamp,run,0);
      gGainCorrection[isec]=(1.+deltaGHV)*(1.+deltaGPT);
      gGainCorrectionPT[isec]=1+deltaGPT;
      gGainCorrectionHV[isec]=1+deltaGHV;
    }
    gTimeStamp=timeStamp;
  }
  if (mode==0){
    if (hasTimeDependent) return gGainCorrection[sector];
    if (!hasTimeDependent) return 1;
  }
  if (mode==1) return gGainCorrection[sector];
  if (mode==2) return gGainCorrectionPT[sector];
  if (mode==3) return gGainCorrectionHV[sector];
  return 1;
}
