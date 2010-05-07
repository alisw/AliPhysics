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


///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// Class providing the calibration parameters by accessing the CDB           //
//                                                                           //
// Request an instance with AliTPCcalibDB::Instance()                        //
// If a new event is processed set the event number with SetRun              //
// Then request the calibration data                                         ////
//
//
// Calibration data:
// 0.)  Altro mapping
//          Simulation      - not yet 
//          Reconstruction  - AliTPCclustererMI::Digits2Clusters(AliRawReader* rawReader)
//
// 1.)  pad by pad calibration -  AliTPCCalPad
//      
//      a.) fPadGainFactor
//          Simulation: AliTPCDigitizer::ExecFast - Multiply by gain
//          Reconstruction : AliTPCclustererMI::Digits2Clusters - Divide by gain  
//
//      b.) fPadNoise -
//          Simulation:        AliTPCDigitizer::ExecFast
//          Reconstruction:    AliTPCclustererMI::FindClusters(AliTPCCalROC * noiseROC)
//                             Noise depending cut on clusters charge (n sigma)
//      c.) fPedestal:
//          Simulation:     Not used yet - To be impleneted - Rounding to the nearest integer
//          Reconstruction: Used in AliTPCclustererMI::Digits2Clusters(AliRawReader* rawReader) 
//                          if data taken without zero suppression  
//                          Currently switch in  fRecoParam->GetCalcPedestal();
//      
//      d.) fPadTime0
//          Simulation:      applied in the AliTPC::MakeSector - adding offset
//          Reconstruction:  AliTPCTransform::Transform() - remove offset
//                           AliTPCTransform::Transform() - to be called
//                           in AliTPCtracker::Transform()      
//
// 
// 2.)  Space points transformation:
//
//      a.) General coordinate tranformation - AliTPCtransform (see $ALICE_ROOT/TPC/AliTPCtransform.cxx)
//          Created on fly - use the other calibration components
//                 Unisochronity  - (substract time0 - pad by pad)
//                 Drift velocity - Currently common drift velocity - functionality of AliTPCParam
//                 ExB effect    
//          Simulation     - Not used directly (the effects are applied one by one (see AliTPC::MakeSector)
//          Reconstruction - 
//                           AliTPCclustererMI::AddCluster
//                           AliTPCtrackerMI::Transform
//      b.) ExB effect calibration - 
//             classes (base class AliTPCExB, implementation- AliTPCExBExact.h  AliTPCExBFirst.h)
//             a.a) Simulation:   applied in the AliTPC::MakeSector - 
//                                calib->GetExB()->CorrectInverse(dxyz0,dxyz1);
//             a.b) Reconstruction -  
//                  
//                  in AliTPCtransform::Correct() - called calib->GetExB()->Correct(dxyz0,dxyz1)
//
//  3.)   cluster error, shape and Q parameterization
//
//
//
///////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <fstream>


#include <AliCDBManager.h>
#include <AliCDBEntry.h>
#include <AliCDBId.h>
#include <AliLog.h>
#include <AliMagF.h>
#include <AliSplineFit.h>
#include <AliCTPTimeParams.h>

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

class AliCDBStorage;
class AliTPCCalDet;
//
//

#include "TFile.h"
#include "TKey.h"
#include "TGraphErrors.h"

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
#include "AliTPCPreprocessorOnline.h"


ClassImp(AliTPCcalibDB)

AliTPCcalibDB* AliTPCcalibDB::fgInstance = 0;
Bool_t AliTPCcalibDB::fgTerminated = kFALSE;
TObjArray    AliTPCcalibDB::fgExBArray;    // array of ExB corrections


//_ singleton implementation __________________________________________________
AliTPCcalibDB* AliTPCcalibDB::Instance()
{
  //
  // Singleton implementation
  // Returns an instance of this class, it is created if neccessary
  //
  
  if (fgTerminated != kFALSE)
    return 0;

  if (fgInstance == 0)
    fgInstance = new AliTPCcalibDB();
  
  return fgInstance;
}

void AliTPCcalibDB::Terminate()
{
  //
  // Singleton implementation
  // Deletes the instance of this class and sets the terminated flag, instances cannot be requested anymore
  // This function can be called several times.
  //
  
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
  fDedxGainFactor(0),
  fPadTime0(0),
  fDistortionMap(0),
  fComposedCorrection(0),
  fPadNoise(0),
  fPedestals(0),
  fCalibRaw(0),
  fDataQA(0),
  fALTROConfigData(0),
  fPulserData(0),
  fCEData(0),
  fTemperature(0),
  fMapping(0),
  fParam(0),
  fClusterParam(0),
  fTimeGainSplines(0),
  fTimeGainSplinesArray(100000),
  fGRPArray(100000),            //! array of GRPs  -  per run  - JUST for calibration studies
  fGRPMaps(100000),            //! array of GRPs  -  per run  - JUST for calibration studies
  fGoofieArray(100000),         //! array of GOOFIE values -per run - Just for calibration studies
  fVoltageArray(100000),
  fTemperatureArray(100000),    //! array of temperature sensors - per run - Just for calibration studies
  fVdriftArray(100000),                 //! array of v drift interfaces
  fDriftCorrectionArray(100000),  //! array of drift correction
  fRunList(100000),              //! run list - indicates try to get the run param 
  fDButil(0),
  fCTPTimeParams(0)
{
  //
  // constructor
  //  
  //
  fgInstance=this;
  Update();    // temporary
}

AliTPCcalibDB::AliTPCcalibDB(const AliTPCcalibDB& ):
  TObject(),
  fRun(-1),
  fTransform(0),
  fExB(0),
  fPadGainFactor(0),
  fDedxGainFactor(0),
  fPadTime0(0),
  fDistortionMap(0),
  fComposedCorrection(0),
  fPadNoise(0),
  fPedestals(0),
  fCalibRaw(0),
  fDataQA(0),
  fALTROConfigData(0),
  fPulserData(0),
  fCEData(0),
  fTemperature(0),
  fMapping(0),
  fParam(0),
  fClusterParam(0),
  fTimeGainSplines(0),
  fTimeGainSplinesArray(100000),
  fGRPArray(0),          //! array of GRPs  -  per run  - JUST for calibration studies
  fGRPMaps(0),          //! array of GRPs  -  per run  - JUST for calibration studies
  fGoofieArray(0),        //! array of GOOFIE values -per run - Just for calibration studies
  fVoltageArray(0),
  fTemperatureArray(0),   //! array of temperature sensors - per run - Just for calibration studies
  fVdriftArray(0),         //! array of v drift interfaces
  fDriftCorrectionArray(0),         //! array of v drift interfaces
  fRunList(0),              //! run list - indicates try to get the run param 
  fDButil(0),
  fCTPTimeParams(0)
{
  //
  // Copy constructor invalid -- singleton implementation
  //
   Error("copy constructor","invalid -- singleton implementation");
}

AliTPCcalibDB& AliTPCcalibDB::operator= (const AliTPCcalibDB& )
{
//
// Singleton implementation - no assignment operator
//
  Error("operator =", "assignment operator not implemented");
  return *this;
}



//_____________________________________________________________________________
AliTPCcalibDB::~AliTPCcalibDB() 
{
  //
  // destructor
  //
  
}
AliTPCCalPad* AliTPCcalibDB::GetDistortionMap(Int_t i) const {
  //
  // get distortion map - due E field distortions
  //
  return (fDistortionMap) ? (AliTPCCalPad*)fDistortionMap->At(i):0;
}

//_____________________________________________________________________________
AliCDBEntry* AliTPCcalibDB::GetCDBEntry(const char* cdbPath)
{
  // 
  // Retrieves an entry with path <cdbPath> from the CDB.
  //
  char chinfo[1000];
    
  AliCDBEntry* entry = AliCDBManager::Instance()->Get(cdbPath, fRun); 
  if (!entry) 
  { 
    sprintf(chinfo,"AliTPCcalibDB: Failed to get entry:\t%s ", cdbPath);
    AliError(chinfo); 
    return 0; 
  }
  return entry;
}


//_____________________________________________________________________________
void AliTPCcalibDB::SetRun(Long64_t run)
{
  //
  // Sets current run number. Calibration data is read from the corresponding file. 
  //  
  if (fRun == run)
    return;  
	fRun = run;
  Update();
}
  


void AliTPCcalibDB::Update(){
  //
  // cache the OCDB entries for simulation, reconstruction, calibration
  //  
  //
  AliCDBEntry * entry=0;
  Bool_t cdbCache = AliCDBManager::Instance()->GetCacheFlag(); // save cache status
  AliCDBManager::Instance()->SetCacheFlag(kTRUE); // activate CDB cache
  fDButil = new AliTPCcalibDButil;   
  //
  entry          = GetCDBEntry("TPC/Calib/PadGainFactor");
  if (entry){
    //if (fPadGainFactor) delete fPadGainFactor;
    entry->SetOwner(kTRUE);
    fPadGainFactor = (AliTPCCalPad*)entry->GetObject();
  }else{
    AliFatal("TPC - Missing calibration entry TPC/Calib/PadGainFactor")
  }
  //
  entry          = GetCDBEntry("TPC/Calib/TimeGain");
  if (entry){
    //if (fTimeGainSplines) delete fTimeGainSplines;
    entry->SetOwner(kTRUE);
    fTimeGainSplines = (TObjArray*)entry->GetObject();
  }else{
    AliFatal("TPC - Missing calibration entry TPC/Calib/Timegain")
  }
  //
  entry          = GetCDBEntry("TPC/Calib/GainFactorDedx");
  if (entry){
    entry->SetOwner(kTRUE);
    fDedxGainFactor = (AliTPCCalPad*)entry->GetObject();
  }else{
    AliFatal("TPC - Missing calibration entry TPC/Calib/gainFactordEdx")
  }
  //
  entry          = GetCDBEntry("TPC/Calib/PadTime0");
  if (entry){
    //if (fPadTime0) delete fPadTime0;
    entry->SetOwner(kTRUE);
    fPadTime0 = (AliTPCCalPad*)entry->GetObject();
  }else{
    AliFatal("TPC - Missing calibration entry")
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
    AliFatal("TPC - Missing calibration entry")
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

  entry          = GetCDBEntry("TPC/Calib/Parameters");
  if (entry){
    //if (fPadNoise) delete fPadNoise;
    entry->SetOwner(kTRUE);
    fParam = (AliTPCParam*)(entry->GetObject()->Clone());
  }else{
    AliFatal("TPC - Missing calibration entry TPC/Calib/Parameters")
  }

  entry          = GetCDBEntry("TPC/Calib/ClusterParam");
  if (entry){
    entry->SetOwner(kTRUE);
    fClusterParam = (AliTPCClusterParam*)(entry->GetObject()->Clone());
  }else{
    AliFatal("TPC - Missing calibration entry")
  }

  //ALTRO configuration data
  entry          = GetCDBEntry("TPC/Calib/AltroConfig");
  if (entry){
    entry->SetOwner(kTRUE);
    fALTROConfigData=(TObjArray*)(entry->GetObject());
  }else{
    AliFatal("TPC - Missing calibration entry")
  }
  
  //Calibration Pulser data
  entry          = GetCDBEntry("TPC/Calib/Pulser");
  if (entry){
    entry->SetOwner(kTRUE);
    fPulserData=(TObjArray*)(entry->GetObject());
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
    AliError("TPC - Missing calibration entry")
  }  
  //TPC space point correction data
  entry          = GetCDBEntry("TPC/Calib/Correction");
  if (entry){
    //entry->SetOwner(kTRUE);
    fComposedCorrection=dynamic_cast<AliTPCCorrection*>(entry->GetObject());
    fComposedCorrection->Init();
  }else{
    AliError("TPC - Missing calibration entry-  TPC/Calib/Correction")
  }  

  //
  if (!fTransform) {
    fTransform=new AliTPCTransform(); 
    fTransform->SetCurrentRun(AliCDBManager::Instance()->GetRun());
  }

  //
  AliCDBManager::Instance()->SetCacheFlag(cdbCache); // reset original CDB cache
}

void AliTPCcalibDB::UpdateNonRec(){
  //
  // Update/Load the parameters which are important for QA studies
  // and not used yet for the reconstruction
  //
   //RAW calibration data
  AliCDBEntry * entry=0;
  entry          = GetCDBEntry("TPC/Calib/Raw");
  if (entry){
    entry->SetOwner(kTRUE);
    TObjArray *arr=(TObjArray*)(entry->GetObject());
    if (arr) fCalibRaw=(AliTPCCalibRaw*)arr->At(0);
  }
  //QA calibration data
  entry          = GetCDBEntry("TPC/Calib/QA");
  if (entry){
    entry->SetOwner(kTRUE);
    fDataQA=dynamic_cast<AliTPCdataQA*>(entry->GetObject());
  }
  // High voltage
  if (fRun>=0){
    entry = AliCDBManager::Instance()->Get("TPC/Calib/HighVoltage",fRun);
    if (entry)  {
      fVoltageArray.AddAt(entry->GetObject(),fRun);
    }
  }

}



void AliTPCcalibDB::CreateObjectList(const Char_t *filename, TObjArray *calibObjects)
{
//
// Create calibration objects and read contents from OCDB
//
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
      
      if ( !sObjType || ! sObjFileName ) continue;
      TString sType(sObjType->GetString());
      TString sFileName(sObjFileName->GetString());
      printf("%s\t%s\n",sType.Data(),sFileName.Data());
      
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
}



void AliTPCcalibDB::MakeTree(const char * fileName, TObjArray * array, const char * mapFileName, AliTPCCalPad* outlierPad, Float_t ltmFraction) {
  //
  // Write a tree with all available information
  // if mapFileName is specified, the Map information are also written to the tree
  // pads specified in outlierPad are not used for calculating statistics
  //  - the same function as AliTPCCalPad::MakeTree - 
  //
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
  //
  // return the RCU trigger configuration register
  //
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
  //
  // return if the FEE readout was triggered on L0
  //
  Int_t mode=GetRCUTriggerConfig();
  if (mode<0) return kFALSE;
  return (mode==1);
}

Bool_t AliTPCcalibDB::IsTrgL1()
{
  //
  // return if the FEE readout was triggered on L1
  //
  Int_t mode=GetRCUTriggerConfig();
  if (mode<0) return kFALSE;
  return (mode==0);
}

void AliTPCcalibDB::RegisterExB(Int_t index, Float_t bz, Bool_t bdelete){
  //
  // Register static ExB correction map
  // index - registration index - used for visualization
  // bz    - bz field in kGaus

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
  //
  // bz filed in KGaus not in tesla
  // Get ExB correction map
  // if doesn't exist - create it
  //
  Int_t index = TMath::Nint(5+bz);
  if (index>fgExBArray.GetEntries()) fgExBArray.Expand((index+1)*2+11);
  if (!fgExBArray.At(index)) AliTPCcalibDB::RegisterExB(index,bz,deleteB);
  return (AliTPCExB*)fgExBArray.At(index);
}


void  AliTPCcalibDB::SetExBField(Float_t bz){
  //
  // Set magnetic filed for ExB correction
  //
  fExB = GetExB(bz,kFALSE);
}

void  AliTPCcalibDB::SetExBField(const AliMagF*   bmap){
  //
  // Set magnetic field for ExB correction
  //
  AliTPCExBFirst *exb  = new  AliTPCExBFirst(bmap,0.88*2.6400e+04,50,50,50);
  AliTPCExB::SetInstance(exb);
  fExB=exb;
}





void AliTPCcalibDB::UpdateRunInformations( Int_t run, Bool_t force){
  //
  // - > Don't use it for reconstruction - Only for Calibration studies
  //
  if (run<=0) return;
  fRun=run;
  AliCDBEntry * entry = 0;
  if (run>= fRunList.fN){
    fRunList.Set(run*2+1);
    fGRPArray.Expand(run*2+1);
    fGRPMaps.Expand(run*2+1);
    fGoofieArray.Expand(run*2+1);
    fVoltageArray.Expand(run*2+1); 
    fTemperatureArray.Expand(run*2+1);
    fVdriftArray.Expand(run*2+1);
    fDriftCorrectionArray.Expand(run*2+1);
    fTimeGainSplinesArray.Expand(run*2+1);
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

	fGRPMaps.AddAt(map,run);
      }
    }
    fGRPArray.AddAt(grpRun,run);
  }
  entry = AliCDBManager::Instance()->Get("TPC/Calib/Goofie",run);
  if (entry){
    fGoofieArray.AddAt(entry->GetObject(),run);
  }
  //
  
  //
  entry = AliCDBManager::Instance()->Get("TPC/Calib/TimeGain",run);
  if (entry)  {
    fTimeGainSplinesArray.AddAt(entry->GetObject(),run);
  }else{
    AliFatal("TPC - Missing calibration entry TimeGain")
  }
  //
  entry = AliCDBManager::Instance()->Get("TPC/Calib/TimeDrift",run);
  if (entry)  {
    fDriftCorrectionArray.AddAt(entry->GetObject(),run);
  }else{
    AliFatal("TPC - Missing calibration entry TimeDrift")
  }
  //
  entry = AliCDBManager::Instance()->Get("TPC/Calib/Temperature",run);
  if (entry)  {
    fTemperatureArray.AddAt(entry->GetObject(),run);
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
    const Double_t kMinP=950.;
    const Double_t kMaxP=1050.;
    const Double_t kMaxdP=10.;
    const Double_t kSigmaCut=4.;
    fDButil->FilterSensor(press,kMinP,kMaxP,kMaxdP,kSigmaCut);
    if (press->GetFit()==0) accept=kFALSE;
  }
  if (press && temp &&accept){
    AliTPCCalibVdrift * vdrift = new AliTPCCalibVdrift(temp, press,0);
    fVdriftArray.AddAt(vdrift,run);
  }
  fDButil->FilterCE(120., 3., 4.,0);
  fDButil->FilterTracks(run, 10.,0);
}


Float_t AliTPCcalibDB::GetGain(Int_t sector, Int_t row, Int_t pad){
  //
  // Get Gain factor for given pad
  //
  AliTPCCalPad *calPad = Instance()->fDedxGainFactor;;
  if (!calPad) return 0;
  return calPad->GetCalROC(sector)->GetValue(row,pad);
}

AliSplineFit* AliTPCcalibDB::GetVdriftSplineFit(const char* name, Int_t run){
  //
  // GetDrift velocity spline fit
  //
  TObjArray *arr=GetTimeVdriftSplineRun(run);
  if (!arr) return 0;
  return dynamic_cast<AliSplineFit*>(arr->FindObject(name));
}

AliSplineFit* AliTPCcalibDB::CreateVdriftSplineFit(const char* graphName, Int_t run){
  //
  // create spline fit from the drift time graph in TimeDrift
  //
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
  //
  // Get GRP object for given run 
  //
  if (run>= ((Instance()->fGRPArray)).GetEntriesFast()){
    Instance()->UpdateRunInformations(run);    
  }
  AliGRPObject * grpRun = dynamic_cast<AliGRPObject *>((Instance()->fGRPArray).At(run));
  if (!grpRun) {
    Instance()->UpdateRunInformations(run);
    grpRun = dynamic_cast<AliGRPObject *>(Instance()->fGRPArray.At(run));
    if (!grpRun) return 0; 
  }
  return grpRun;
}

TMap *  AliTPCcalibDB::GetGRPMap(Int_t run){
  //
  // Get GRP map for given run
  //
  TMap * grpRun = dynamic_cast<TMap *>((Instance()->fGRPMaps).At(run));
  if (!grpRun) {
    Instance()->UpdateRunInformations(run);
    grpRun = dynamic_cast<TMap *>(Instance()->fGRPMaps.At(run));
    if (!grpRun) return 0; 
  }
  return grpRun;
}


AliDCSSensor * AliTPCcalibDB::GetPressureSensor(Int_t run, Int_t type){
  //
  // Get Pressure sensor
  // run  = run number
  // type = 0 - Cavern pressure
  //        1 - Suface pressure
  // First try to get if trom map - if existing  (Old format of data storing)
  //


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
  AliGRPObject * grpRun = dynamic_cast<AliGRPObject *>(fGRPArray.At(run)); 
  if (!grpRun) {
    UpdateRunInformations(run);
    grpRun = dynamic_cast<AliGRPObject *>(fGRPArray.At(run));
    if (!grpRun) return 0; 
  }
  AliDCSSensor * sensor = grpRun->GetCavernAtmosPressure();
  if (type==1) sensor = grpRun->GetSurfaceAtmosPressure();
  return sensor; 
}

AliTPCSensorTempArray * AliTPCcalibDB::GetTemperatureSensor(Int_t run){
  //
  // Get temperature sensor array
  //
  AliTPCSensorTempArray * tempArray = (AliTPCSensorTempArray *)fTemperatureArray.At(run);
  if (!tempArray) {
    UpdateRunInformations(run);
    tempArray = (AliTPCSensorTempArray *)fTemperatureArray.At(run);
  }
  return tempArray;
}


TObjArray * AliTPCcalibDB::GetTimeGainSplinesRun(Int_t run){
  //
  // Get temperature sensor array
  //
  TObjArray * gainSplines = (TObjArray *)fTimeGainSplinesArray.At(run);
  if (!gainSplines) {
    UpdateRunInformations(run);
    gainSplines = (TObjArray *)fTimeGainSplinesArray.At(run);
  }
  return gainSplines;
}

TObjArray * AliTPCcalibDB::GetTimeVdriftSplineRun(Int_t run){
  //
  // Get drift spline array
  //
  TObjArray * driftSplines = (TObjArray *)fDriftCorrectionArray.At(run);
  if (!driftSplines) {
    UpdateRunInformations(run);
    driftSplines = (TObjArray *)fDriftCorrectionArray.At(run);
  }
  return driftSplines;
}

AliDCSSensorArray * AliTPCcalibDB::GetVoltageSensors(Int_t run){
  //
  // Get temperature sensor array
  //
  AliDCSSensorArray * voltageArray = (AliDCSSensorArray *)fVoltageArray.At(run);
  if (!voltageArray) {
    UpdateRunInformations(run);
    voltageArray = (AliDCSSensorArray *)fVoltageArray.At(run);
  }
  return voltageArray;
}

AliDCSSensorArray * AliTPCcalibDB::GetGoofieSensors(Int_t run){
  //
  // Get temperature sensor array
  //
  AliDCSSensorArray * goofieArray = (AliDCSSensorArray *)fGoofieArray.At(run);
  if (!goofieArray) {
    UpdateRunInformations(run);
    goofieArray = (AliDCSSensorArray *)fGoofieArray.At(run);
  }
  return goofieArray;
}



AliTPCCalibVdrift *     AliTPCcalibDB::GetVdrift(Int_t run){
  //
  // Get the interface to the the vdrift 
  //
  AliTPCCalibVdrift  * vdrift = (AliTPCCalibVdrift*)fVdriftArray.At(run);
  if (!vdrift) {
    UpdateRunInformations(run);
    vdrift= (AliTPCCalibVdrift*)fVdriftArray.At(run);
  }
  return vdrift;
}

Float_t AliTPCcalibDB::GetCEdriftTime(Int_t run, Int_t sector, Double_t timeStamp, Int_t *entries)
{
  //
  // GetCE drift time information for 'sector'
  // sector 72 is the mean drift time of the A-Side
  // sector 73 is the mean drift time of the C-Side
  // it timestamp==-1 return mean value
  //
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
  //
  // GetCE mean charge for 'sector'
  // it timestamp==-1 return mean value
  //
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
  //
  // Get Value for a DCS sensor 'sensorName', run 'run' at time 'timeStamp'
  //
  Float_t val=0;
  const TString sensorNameString(sensorName);
  AliDCSSensor *sensor = arr->GetSensor(sensorNameString);
  if (!sensor) return val;
  //use the dcs graph if possible
  TGraph *gr=sensor->GetGraph();
  if (gr){
    for (Int_t ipoint=0;ipoint<gr->GetN();++ipoint){
      Double_t x,y;
      gr->GetPoint(ipoint,x,y);
      Int_t time=TMath::Nint(sensor->GetStartTime()+x*3600); //time in graph is hours
      if (time<timeStamp) continue;
      val=y;
      break;
    }
    //if val is still 0, test if if the requested time if within 5min of the first/last
    //data point. If this is the case return the firs/last entry
    //the timestamps might not be syncronised for all calibration types, sometimes a 'pre'
    //and 'pos' period is requested. Especially to the HV this is not the case!
    //first point
    if (val==0 ){
      Double_t x,y;
      gr->GetPoint(0,x,y);
      Int_t time=TMath::Nint(sensor->GetStartTime()+x*3600); //time in graph is hours
      if ((time-timeStamp)<5*60) val=y;
    }
    //last point
    if (val==0 ){
      Double_t x,y;
      gr->GetPoint(gr->GetN()-1,x,y);
      Int_t time=TMath::Nint(sensor->GetStartTime()+x*3600); //time in graph is hours
      if ((timeStamp-time)<5*60) val=y;
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
  //
  // Get mean Value for a DCS sensor 'sensorName' during run 'run'
  //
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
    val/=10;
    val=(Float_t)TMath::Floor(val * TMath::Power(10., sigDigits) + .5) / TMath::Power(10., sigDigits);
    val*=10;
  }
  return val;
}

Float_t AliTPCcalibDB::GetChamberHighVoltage(Int_t run, Int_t sector, Int_t timeStamp, Int_t sigDigits) {
  //
  // return the chamber HV for given run and time: 0-35 IROC, 36-72 OROC
  // if timeStamp==-1 return mean value
  //
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
  if (timeStamp==-1){
    val=AliTPCcalibDB::GetDCSSensorMeanValue(voltageArray, sensorName.Data(),sigDigits);
  } else {
    val=AliTPCcalibDB::GetDCSSensorValue(voltageArray, timeStamp, sensorName.Data(),sigDigits);
  }
  return val;
}
Float_t AliTPCcalibDB::GetSkirtVoltage(Int_t run, Int_t sector, Int_t timeStamp, Int_t sigDigits)
{
  //
  // Get the skirt voltage for 'run' at 'timeStamp' and 'sector': 0-35 IROC, 36-72 OROC
  // type corresponds to the following: 0 - IROC A-Side; 1 - IROC C-Side; 2 - OROC A-Side; 3 - OROC C-Side
  // if timeStamp==-1 return the mean value for the run
  //
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
  //
  // Get the cover voltage for run 'run' at time 'timeStamp'
  // type corresponds to the following: 0 - IROC A-Side; 1 - IROC C-Side; 2 - OROC A-Side; 3 - OROC C-Side
  // if timeStamp==-1 return the mean value for the run
  //
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
  //
  // Get the GG offset voltage for run 'run' at time 'timeStamp'
  // type corresponds to the following: 0 - IROC A-Side; 1 - IROC C-Side; 2 - OROC A-Side; 3 - OROC C-Side
  // if timeStamp==-1 return the mean value for the run
  //
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
  //
  // Get the GG offset voltage for run 'run' at time 'timeStamp'
  // type corresponds to the following: 0 - IROC A-Side; 1 - IROC C-Side; 2 - OROC A-Side; 3 - OROC C-Side
  // if timeStamp==-1 return the mean value for the run
  //
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
  //
  // Get the GG offset voltage for run 'run' at time 'timeStamp'
  // type corresponds to the following: 0 - IROC A-Side; 1 - IROC C-Side; 2 - OROC A-Side; 3 - OROC C-Side
  // if timeStamp==-1 return the mean value for the run
  //
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
  //
  // GetPressure for given time stamp and runt
  //
  TTimeStamp stamp(timeStamp);
  AliDCSSensor * sensor = Instance()->GetPressureSensor(run,type);
  if (!sensor) return 0;
  return sensor->GetValue(stamp);
}

Float_t AliTPCcalibDB::GetL3Current(Int_t run, Int_t statType){
  //
  // return L3 current
  // stat type is: AliGRPObject::Stats: kMean = 0, kTruncMean = 1, kMedian = 2, kSDMean = 3, kSDMedian = 4
  //
  Float_t current=-1;
  AliGRPObject *grp=AliTPCcalibDB::GetGRP(run);
  if (grp) current=grp->GetL3Current((AliGRPObject::Stats)statType);
  return current;
}

Float_t AliTPCcalibDB::GetBz(Int_t run){
  //
  // calculate BZ in T from L3 current
  //
  Float_t bz=-1;
  Float_t current=AliTPCcalibDB::GetL3Current(run);
  if (current>-1) bz=5*current/30000.*.1;
  return bz;
}

Char_t  AliTPCcalibDB::GetL3Polarity(Int_t run) {
  //
  // get l3 polarity from GRP
  //
  Char_t pol=-100;
  AliGRPObject *grp=AliTPCcalibDB::GetGRP(run);
  if (grp) pol=grp->GetL3Polarity();
  return pol;
}

TString AliTPCcalibDB::GetRunType(Int_t run){
  //
  // return run type from grp
  //

//   TString type("UNKNOWN");
  AliGRPObject *grp=AliTPCcalibDB::GetGRP(run);
  if (grp) return grp->GetRunType();
  return "UNKNOWN";
}

Float_t AliTPCcalibDB::GetValueGoofie(Int_t timeStamp, Int_t run, Int_t type){
  //
  // GetPressure for given time stamp and runt
  //
  TTimeStamp stamp(timeStamp);
  AliDCSSensorArray* goofieArray = AliTPCcalibDB::Instance()->GetGoofieSensors(run);
  if (!goofieArray) return 0;
  AliDCSSensor *sensor = goofieArray->GetSensor(type);
  return sensor->GetValue(stamp);
}






Bool_t  AliTPCcalibDB::GetTemperatureFit(Int_t timeStamp, Int_t run, Int_t side,TVectorD& fit){
  //
  // GetTmeparature fit at parameter for given time stamp
  //
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
  //
  // Get mean temperature
  // 
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
  //
  // Get relative P/T 
  // time - absolute time
  // run  - run number
  // side - 0 - A side   1-C side
  AliTPCCalibVdrift * vdrift =  Instance()->GetVdrift(run);
  if (!vdrift) return 0;
  return vdrift->GetPTRelative(timeSec,side);
}

AliGRPObject * AliTPCcalibDB::MakeGRPObjectFromMap(TMap *map){
  //
  // Function to covert old GRP run information from TMap to GRPObject
  //
  //  TMap * map = AliTPCcalibDB::GetGRPMap(52406);
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
  grpRun->SetSurfaceAtmosPressure(sensor);
  return grpRun;
}

Bool_t AliTPCcalibDB::CreateGUITree(Int_t run, const char* filename)
{
  //
  // Create a gui tree for run number 'run'
  //

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
  //
  //
  //
  if (!AliCDBManager::Instance()->GetDefaultStorage()){
    AliError("Default Storage not set. Cannot create calibration Tree!");
    return kFALSE;
  }
  UpdateNonRec();  // load all infromation now

  AliTPCPreprocessorOnline prep;
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
  if (file.IsNull()) file=Form("guiTreeRun_%d.root",fRun);
  prep.DumpToFile(file.Data());
  return kTRUE;
}

Bool_t AliTPCcalibDB::CreateRefFile(Int_t run, const char* filename)
{
  //
  // Create a gui tree for run number 'run'
  //
  
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
  //
  // Get time dependent drift velocity correction
  // multiplication factor        vd = vdnom *(1+vdriftcorr)
  // Arguments:
  // mode determines the algorith how to combine the Laser Track, LaserCE and physics tracks
  // timestamp - timestamp
  // run       - run number
  // side      - the drift velocity per side (possible for laser and CE)
  //
  // Notice - Extrapolation outside of calibration range  - using constant function
  //
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

  return result;
}

Double_t AliTPCcalibDB::GetTime0CorrectionTime(Int_t timeStamp, Int_t run, Int_t /*side*/, Int_t mode){
  //
  // Get time dependent time 0 (trigger delay in cm) correction
  // additive correction        time0 = time0+ GetTime0CorrectionTime
  // Value etracted combining the vdrift correction using laser tracks and CE and the physics track matchin
  // Arguments:
  // mode determines the algorith how to combine the Laser Track and physics tracks
  // timestamp - timestamp
  // run       - run number
  // side      - the drift velocity per side (possible for laser and CE)
  //
  // Notice - Extrapolation outside of calibration range  - using constant function
  //
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
  //
  // Get global y correction drift velocity correction factor
  // additive factor        vd = vdnom*(1+GetVDriftCorrectionGy *gy)
  // Value etracted combining the vdrift correction using laser tracks and CE
  // Arguments:
  // mode determines the algorith how to combine the Laser Track, LaserCE
  // timestamp - timestamp
  // run       - run number
  // side      - the drift velocity gy correction per side (CE and Laser tracks)
  //
  // Notice - Extrapolation outside of calibration range  - using constant function
  // 
  if (run<=0 && fTransform) run = fTransform->GetCurrentRunNumber();
  UpdateRunInformations(run,kFALSE);
  TObjArray *array =AliTPCcalibDB::Instance()->GetTimeVdriftSplineRun(run);
  if (!array) return 0;
  TGraphErrors *laserA= (TGraphErrors*)array->FindObject("GRAPH_MEAN_GLOBALYGRADIENT_LASER_ALL_A");
  TGraphErrors *laserC= (TGraphErrors*)array->FindObject("GRAPH_MEAN_GLOBALYGRADIENT_LASER_ALL_C");
  
  Double_t result=0;
  if (laserA && laserC){
   result= (laserA->Eval(timeStamp)+laserC->Eval(timeStamp))*0.5;
  }
  if (laserA && side==0){
    result = (laserA->Eval(timeStamp));
  }
  if (laserC &&side==1){
    result = (laserC->Eval(timeStamp));
  }
  return -result/250.; //normalized before
}

AliTPCCalPad* AliTPCcalibDB::MakeDeadMap(Double_t notInMap, const char* nameMappingFile) {
//
//   Read list of active DDLs from OCDB entry
//   Generate and return AliTPCCalPad containing 1 for all pads in active DDLs,
//   0 for all pads in non-active DDLs. 
//   For DDLs with missing status information (no DCS input point to Shuttle),
//     the value of the AliTPCCalPad entry is determined by the parameter
//     notInMap (default value 1)
//
  char chinfo[1000];
   
  TFile *fileMapping = new TFile(nameMappingFile, "read");
  AliTPCmapper *mapping = (AliTPCmapper*) fileMapping->Get("tpcMapping");
  if (!mapping) {
    sprintf(chinfo,"Failed to get mapping object from %s.  ...\n", nameMappingFile);
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
    Int_t patch = mapping->GetPatchFromEquipmentID(idDDL);   
    Int_t roc=mapping->GetRocFromEquipmentID(idDDL);
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


