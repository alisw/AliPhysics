// --- ROOT system
#include <TFile.h>
#include <TClonesArray.h>
#include <TList.h>
#include <TObjString.h>
#include <TTimeStamp.h>

#include "AliZDCPreprocessor.h"
#include "AliCDBManager.h"
#include "AliCDBEntry.h"
#include "AliCDBMetaData.h"
#include "AliDCSValue.h"
#include "AliAlignObj.h"
#include "AliAlignObjParams.h"
#include "AliLog.h"
#include "AliZDCDataDCS.h"
#include "AliZDCChMap.h"
#include "AliZDCPedestals.h"
#include "AliZDCLaserCalib.h"
#include "AliZDCEnCalib.h"
#include "AliZDCTowerCalib.h"
#include "AliZDCMBCalib.h"
#include "AliZDCTDCCalib.h"

/////////////////////////////////////////////////////////////////////
//								   //
// Class implementing Shuttle ZDC pre-processor.	           //
// It takes data from DCS and DAQ and writes calibration objects   //
// in the OCDB and reference values/histos in the ReferenceData.   //
//								   //
/////////////////////////////////////////////////////////////////////

// ******************************************************************
//    RETURN CODES:
// return 0 : everything OK
// return 1 : no DCS input data Map
// return 2 : error storing DCS data in RefData 
// return 3 : error storing alignment object in OCDB
// return 4 : error in ZDCMapping.dat file retrieved from DAQ FXS (not existing|empty|corrupted)
// return 5 : error storing mapping obj. in OCDB
// return 6 : error storing energy calibration obj. in OCDB
// return 7 : error storing tower inter-calibration obj. in OCDB
// return 8 : error in ZDCEnergyCalib.dat file retrieved from DAQ FXS 
// return 9 : error in ZDCTowerCalib.dat file retrieved from DAQ FXS 
// return 10: error in ZDCPedestal.dat file retrieved from DAQ FXS 
// return 11: error storing pedestal calibration obj. in OCDB
// return 12: error in ZDCPedHisto.root file retrieved from DAQ FXS 
// return 13: error storing pedestal histos in RefData 
// return 14: error in ZDCLaserCalib.dat file retrieved from DAQ FXS 
// return 15: error storing laser calibration obj. in OCDB
// return 16: error in ZDCLaserHisto.root file retrieved from DAQ FXS 
// return 17: error storing laser histos in RefData 
// return 18: error in ZDCMBCalib.root file retrieved from DAQ FXS 
// return 19: error storing MB calibration obj. in OCDB
// return 20: error in ZDCTDCCalib.root file retrieved from DAQ FXS
// return 21: error in storing TDC calibration obj. in OCDB
// return 22: error in ZDCTDCHisto.root file retrieved from DAQ FXS
// Return 23: error storing TDC reference histos in RefData
// ******************************************************************

ClassImp(AliZDCPreprocessor)

//______________________________________________________________________________________________
AliZDCPreprocessor::AliZDCPreprocessor(AliShuttleInterface* shuttle) :
  AliPreprocessor("ZDC", shuttle),
  fData(0)
{
  // constructor
  // May 2009 - run types updated according to
  // http://alice-ecs.web.cern.ch/alice-ecs/runtypes_3.36.html
  AddRunType("STANDALONE_PEDESTAL");
  AddRunType("STANDALONE_LASER");
  AddRunType("STANDALONE_COSMIC");
  AddRunType("CALIBRATION_EMD");
  AddRunType("CALIBRATION_MB");
  AddRunType("CALIBRATION_CENTRAL");
  AddRunType("CALIBRATION_SEMICENTRAL");
  AddRunType("CALIBRATION_BC");
  AddRunType("PHYSICS");
}


//______________________________________________________________________________________________
AliZDCPreprocessor::~AliZDCPreprocessor()
{
  // destructor
}


//______________________________________________________________________________________________
void AliZDCPreprocessor::Initialize(Int_t run, UInt_t startTime, UInt_t endTime)
{
  // Creates AliZDCDataDCS object

  AliPreprocessor::Initialize(run, startTime, endTime);

  AliDebug(2,Form("\n\tRun %d \n\tStartTime %s \n\tEndTime %s \n\tStartTime DCS Query %s \n\tEndTime DCS Query %s", run,
		TTimeStamp(startTime).AsString(),
		TTimeStamp(endTime).AsString(), ((TTimeStamp)GetStartTimeDCSQuery()).AsString(), ((TTimeStamp)GetEndTimeDCSQuery()).AsString()));

  fRun = run;
  fStartTime = startTime;
  fEndTime = endTime;

  fData = new AliZDCDataDCS(fRun, fStartTime, fEndTime, GetStartTimeDCSQuery(), GetEndTimeDCSQuery());
}

//_____________________________________________________________________________
Bool_t AliZDCPreprocessor::ProcessDCS(){

  // tells whether DCS should be processed or not

  TString runType = GetRunType();
  Log(Form("RunType %s",runType.Data()));

  if (runType=="STANDALONE_COSMIC" || runType=="STANDALONE_PEDESTAL"){
    return kFALSE;
  }

  return kTRUE;
}

//______________________________________________________________________________________________
UInt_t AliZDCPreprocessor::ProcessDCSData(TMap* dcsAliasMap)
{
  
  // Fills data into a AliZDCDataDCS object
  if(!dcsAliasMap){
    Log(" No DCS map found: ZDC exiting from Shuttle");
    if(fData){
      delete fData;
      fData = 0;
    }
    return 1;
  }

  Log(Form("Processing data from DCS"));
   
  // The processing of the DCS input data is forwarded to AliZDCDataDCS
  //dcsAliasMap->Print(""); 
  Bool_t resDCSProcess = fData->ProcessData(*dcsAliasMap);
  if(resDCSProcess==kFALSE){
    Log(" Problems in processing DCS DP");
    return 1;
  }  
  
  // ------------------------------------------------------
  // Change introduced 26/9/09 in order NOT to process the
  // HV DP since some of them are never found in amanda DB
  // ------------------------------------------------------
  // Store DCS data as reference
  AliCDBMetaData metadata;
  metadata.SetResponsible("Chiara Oppedisano");
  metadata.SetComment("DCS DP TMap for ZDC");
  Bool_t resDCSRef = kTRUE;
  resDCSRef = StoreReferenceData("DCS","Data", dcsAliasMap, &metadata);
  
  if(resDCSRef==kFALSE) return 2;

  // --- Writing ZDC table positions into alignment object
  TClonesArray *array = new TClonesArray("AliAlignObjParams",10);
  TClonesArray &alobj = *array;
  AliAlignObjParams a;
  Double_t dx=0., dz=0., dpsi=0., dtheta=0., dphi=0.;
  // Vertical table position in mm from DCS
  Double_t dyZN1 = (Double_t) (fData->GetAlignData(0)/10.);
  Double_t dyZP1 = (Double_t) (fData->GetAlignData(1)/10.);
  Double_t dyZN2 = (Double_t) (fData->GetAlignData(2)/10.);
  Double_t dyZP2 = (Double_t) (fData->GetAlignData(3)/10.);
  //
  const char *n1ZDC="ZDC/NeutronZDC_C";  
  const char *p1ZDC="ZDC/ProtonZDC_C";
  const char *n2ZDC="ZDC/NeutronZDC_A";
  const char *p2ZDC="ZDC/ProtonZDC_A";
  //
  UShort_t iIndex=0;
  AliGeomManager::ELayerID iLayer = AliGeomManager::kInvalidLayer;
  UShort_t volid = AliGeomManager::LayerToVolUID(iLayer,iIndex);
  //
  new(alobj[0]) AliAlignObjParams(n1ZDC, volid, dx, dyZN1, dz, dpsi, dtheta, dphi, kTRUE);
  new(alobj[1]) AliAlignObjParams(p1ZDC, volid, dx, dyZP1, dz, dpsi, dtheta, dphi, kTRUE);
  new(alobj[2]) AliAlignObjParams(n2ZDC, volid, dx, dyZN2, dz, dpsi, dtheta, dphi, kTRUE);
  new(alobj[3]) AliAlignObjParams(p2ZDC, volid, dx, dyZP2, dz, dpsi, dtheta, dphi, kTRUE);
 
  // save in CDB storage
  AliCDBMetaData mdDCS;
  mdDCS.SetResponsible("Chiara Oppedisano");
  mdDCS.SetComment("Alignment object for ZDC");
  Bool_t resultAl = Store("Align","Data", array, &mdDCS, 0, kFALSE);
  if(resultAl==kFALSE)  return 3;
  
  return 0;
}

//______________________________________________________________________________________________
UInt_t AliZDCPreprocessor::ProcessChMap()
{ 
  const int kNModules=10, kNch=48, kNScch=32, kNtdcch=32;
  
  // Reading the file for mapping from FXS
  TList* daqSource = GetFileSources(kDAQ, "MAPPING");
  if(!daqSource){
    AliError(Form("No sources for file ZDCChMapping.dat in run %d ", fRun));
    return 4;
  }
  if(daqSource->GetEntries()==0) return 4;
  Log("\t List of DAQ sources for MAPPING id: "); daqSource->Print();
  //
  TIter iter(daqSource);
  TObjString* source = 0;
  Int_t isou = 0;
  Int_t modMap[kNModules][3], adcMap[kNch][6], scMap[kNScch][6], tdcMap[kNtdcch][4]; 
  //
  while((source = dynamic_cast<TObjString*> (iter.Next()))){
     TString fileName = GetFile(kDAQ, "MAPPING", source->GetName());
     Log(Form("\t Getting file #%d: ZDCChMapping.dat from %s\n",++isou, source->GetName()));

     if(fileName.Length() <= 0){
       Log(Form("No file from source %s!", source->GetName()));
       return 4;
     }
     // --- Reading file with calibration data
     //const char* fname = fileName.Data();
     if(fileName){
       FILE *file;
       if((file = fopen(fileName,"r")) == NULL){
	 printf("Cannot open file %s \n",fileName.Data());
         return 4;
       }
       Log(Form("File %s connected to process data for ADC mapping", fileName.Data()));
       //
       for(Int_t j=0; j<kNch; j++){	  
           for(Int_t k=0; k<6; k++){
             int read = fscanf(file,"%d",&adcMap[j][k]);
	     if(read == 0) AliDebug(3," Failing in reading data from mapping file");
           }
       }
       for(Int_t j=kNch; j<kNch+kNScch; j++){	  
           for(Int_t k=0; k<6; k++){
             int read = fscanf(file,"%d",&scMap[j-kNch][k]);
	     if(read == 0) AliDebug(3," Failing in reading data from mapping file");
           }
       }
       for(Int_t j=kNch+kNScch; j<kNch+kNScch+kNtdcch; j++){	  
           for(Int_t k=0; k<4; k++){
             int read = fscanf(file,"%d",&tdcMap[j-kNch-kNScch][k]);
	     if(read == 0) AliDebug(3," Failing in reading data from mapping file");
           }
       }
       for(Int_t j=kNch+kNScch+kNtdcch; j<kNch+kNScch+kNtdcch+kNModules; j++){	  
           for(Int_t k=0; k<3; k++){
             int read = fscanf(file,"%d",&modMap[j-kNch-kNScch-kNtdcch][k]);
	     if(read == 0) AliDebug(3," Failing in reading data from mapping file");
           }
       }
       fclose(file);
     }
     else{
       Log(Form("File %s not found", fileName.Data()));
       return 4;
     }
  }
  
  // Store the currently read map ONLY IF it is different
  // from the entry in the OCDB
  Bool_t adcMapUpdated=kFALSE, scMapUpdated=kFALSE, tdcMapUpdated=kFALSE;
  Bool_t updateOCDB = kFALSE;
  
  AliCDBEntry *cdbEntry = GetFromOCDB("Calib","ChMap");
  if(!cdbEntry){
    Log(" No existing CDB entry for ADC mapping");
    updateOCDB = kTRUE;
  }
  else{
    AliZDCChMap *chMap = (AliZDCChMap*) cdbEntry->GetObject();
    for(Int_t i=0; i<kNch; i++){
      if(  (adcMap[i][1] != chMap->GetADCModule(i)) 
        || (adcMap[i][2] != chMap->GetADCChannel(i)) 
	|| (adcMap[i][3] != chMap->GetADCSignalCode(i)) 
	|| (adcMap[i][4] != chMap->GetDetector(i)) 
	|| (adcMap[i][5] != chMap->GetSector(i))) 
	 adcMapUpdated = kTRUE;
    }
    for(Int_t i=0; i<kNScch; i++){
      if(  (scMap[i][2] != chMap->GetScChannel(i)) 
	|| (scMap[i][3] != chMap->GetScSignalCode(i)) )
	 scMapUpdated = kTRUE;
    }
    for(Int_t i=0; i<kNtdcch; i++){
      if(  (tdcMap[i][2] != chMap->GetTDCChannel(i)) 
	|| (tdcMap[i][3] != chMap->GetTDCSignalCode(i)))
	 tdcMapUpdated = kTRUE;
    }
  }
  if(adcMapUpdated || scMapUpdated || tdcMapUpdated) updateOCDB = kTRUE;
  //
  Bool_t resChMapStore = kTRUE;
  if(updateOCDB==kTRUE){
    Log(" A new entry ZDC/Calib/ChMap will be created");
    //
    // --- Initializing mapping calibration object
    AliZDCChMap *mapCalib = new AliZDCChMap("ZDC");
    // Writing channel map in the OCDB
    for(Int_t k=0; k<kNModules; k++){
      mapCalib->SetModuleMap(k, modMap[k][0], modMap[k][1], modMap[k][2]);
    }
    for(Int_t k=0; k<kNch; k++){
      mapCalib->SetADCModule(k,adcMap[k][1]);
      mapCalib->SetADCChannel(k,adcMap[k][2]);
      mapCalib->SetADCSignalCode(k,adcMap[k][3]);
      mapCalib->SetDetector(k,adcMap[k][4]);
      mapCalib->SetSector(k,adcMap[k][5]);
    }
    for(Int_t k=0; k<kNScch; k++){
       mapCalib->SetScChannel(k, scMap[k][2]);
       mapCalib->SetScSignalCode(k, scMap[k][3]);
       mapCalib->SetScDetector(k, scMap[k][4]);
       mapCalib->SetScSector(k, scMap[k][5]);
    }
    for(Int_t k=0; k<kNtdcch; k++){
       mapCalib->SetTDCChannel(k, tdcMap[k][2]);
       mapCalib->SetTDCSignalCode(k, tdcMap[k][3]);
    }
    //
    //mapCalib->Print("");
    // 
    AliCDBMetaData metaData;
    metaData.SetBeamPeriod(0);
    metaData.SetResponsible("Chiara Oppedisano");
    metaData.SetComment("AliZDCChMap object created by ZDC preprocessor");  
    //
    resChMapStore = Store("Calib","ChMap",mapCalib, &metaData, 0, kTRUE);
    printf("  Mapping object stored in OCDB\n");
  }
  else{
    Log(" ZDC/Calib/ChMap entry in OCDB is valid and won't be updated");
    resChMapStore = kTRUE;
  }
  delete daqSource; daqSource=0;
  
  TString runType = GetRunType();
  if(runType.CompareTo("PHYSICS")==0){
    Log(Form("RunType %s -> producing TDC calibration data",runType.Data()));
    
    // Reading the file for mapping from FXS
    TList* daqSourcetdc = GetFileSources(kDAQ, "TDCDATA");
    if(!daqSourcetdc){
      AliError(Form("No sources for file ZDCChMappingTDCCalib.dat in run %d ", fRun));
      return 20;
    }
    if(daqSourcetdc->GetEntries()==0) return 20;
    Log("\t List of DAQ sources for TDCDATA id: "); daqSourcetdc->Print();
    //
    Bool_t resTDCcal = kTRUE;
    TIter itertdc(daqSourcetdc);
    TObjString* sourcetdc = 0;
    Int_t isoutdc = 0;
    //
    while((sourcetdc = dynamic_cast<TObjString*> (itertdc.Next()))){
     TString fileNametdc = GetFile(kDAQ, "TDCDATA", sourcetdc->GetName());
     Log(Form("\t Getting file #%d: ZDCTDCdata.dat from %s\n",++isoutdc, sourcetdc->GetName()));

     if(fileNametdc.Length() <= 0){
       Log(Form("No file from source %s!", sourcetdc->GetName()));
       return 20;
     }
     // --- Initializing TDC calibration object
     AliZDCTDCCalib *tdcCalib = new AliZDCTDCCalib("ZDC");
     // --- Reading file with calibration data
     //const char* fname = fileName.Data();
     if(fileNametdc){
       FILE *filetdc;
       if((filetdc = fopen(fileNametdc,"r")) == NULL){
	 printf("Cannot open file %s \n",fileNametdc.Data());
         return 20;
       }
       Log(Form("File %s connected to process TDC data", fileNametdc.Data()));
       //
       Float_t tdcMean[6][2];
       for(Int_t it=0; it<6; it++){
         for(Int_t iu=0; iu<2; iu++) tdcMean[it][iu]=0.;
       }
       for(Int_t k=0; k<6; k++){
        for(Int_t j=0; j<2; j++){
           int leggi = fscanf(filetdc,"%f",&tdcMean[k][j]);
	   if(leggi==0) AliDebug(3," Failing reading data from tdc file");
	   tdcCalib->SetMeanTDC(k, tdcMean[k][0]);
	   tdcCalib->SetWidthTDC(k, tdcMean[k][1]);
	}
       }
       fclose(filetdc);
     }
     else{
       Log(Form("File %s not found", fileNametdc.Data()));
       return 20;
     }
     //
     AliCDBMetaData metaData;
     metaData.SetBeamPeriod(0);
     metaData.SetResponsible("Chiara Oppedisano");
     metaData.SetComment("Filling AliZDCTDCCalib object");  
     //
     resTDCcal = Store("Calib","TDCCalib",tdcCalib, &metaData, 0, kTRUE);
     if(resTDCcal==kFALSE) return 21;
    }
    delete daqSourcetdc; daqSourcetdc = 0;

    Bool_t restdcHist = kTRUE;
    TList* daqSourceH = GetFileSources(kDAQ, "TDCHISTOS");
    if(!daqSourceH){
      Log(Form("No source for TDCHISTOS id run %d !", fRun));
      return 22;
    }
    Log("\t List of DAQ sources for TDCHISTOS id: "); daqSourceH->Print();
    //
    TIter iterH(daqSourceH);
    TObjString* sourceH = 0;
    Int_t iH=0;
    while((sourceH = dynamic_cast<TObjString*> (iterH.Next()))){
     TString stringTDCFileName = GetFile(kDAQ, "TDCHISTOS", sourceH->GetName());
     if(stringTDCFileName.Length() <= 0){
     	Log(Form("No TDCHISTOS file from source %s!", sourceH->GetName()));
        return 22;
     }
     const char* tdcFileName = stringTDCFileName.Data();
     Log(Form("\t Getting file #%d: %s from %s\n",++iH, tdcFileName, sourceH->GetName()));
     restdcHist = StoreReferenceFile(tdcFileName, "tdcReference.root");
     if(restdcHist==kFALSE) return 23;
    }
    delete daqSourceH; daqSourceH=0;
  }
  
    if(resChMapStore==kFALSE) return 5;
    else return 0;

}

//______________________________________________________________________________________________
UInt_t AliZDCPreprocessor::ProcessppData()
{
   Bool_t resEnCal=kTRUE, resTowCal=kTRUE;
  
   // *********** Energy calibration
   // --- Cheking if there is already the entry in the OCDB
   AliCDBEntry *cdbEnEntry = GetFromOCDB("Calib", "EnergyCalib");
   if(!cdbEnEntry){   
     Log(Form(" ZDC/Calib/EnergyCalib entry will be created"));
     // --- Initializing calibration object
     AliCDBMetaData metaData;
     metaData.SetBeamPeriod(0);
     metaData.SetResponsible("Chiara Oppedisano");
     //
     AliZDCEnCalib *eCalib = new AliZDCEnCalib("ZDC");
     for(Int_t j=0; j<6; j++) eCalib->SetEnCalib(j,1.);
     metaData.SetComment("AliZDCEnCalib object");  
     //eCalib->Print("");
     resEnCal = Store("Calib", "EnergyCalib", eCalib, &metaData, 0, kTRUE);
   }
   else{ 
     // if entry exists it is still valid (=1 for all runs!)
     Log(Form(" Valid ZDC/Calib/EnergyCalib object already existing in OCDB!!!"));
     resEnCal = kTRUE;
   }
   
   if(resEnCal==kFALSE) return 6;

   //
   // *********** Tower inter-calibration
   // --- Cheking if there is already the entry in the OCDB
   AliCDBEntry *cdbTowEntry = GetFromOCDB("Calib", "TowerCalib");
   if(!cdbTowEntry){   
     AliZDCTowerCalib *towCalib = new AliZDCTowerCalib("ZDC");
     for(Int_t j=0; j<5; j++){  
        towCalib->SetZN1EqualCoeff(j, 1.);
        towCalib->SetZP1EqualCoeff(j, 1.);
        towCalib->SetZN2EqualCoeff(j, 1.);
        towCalib->SetZP2EqualCoeff(j, 1.);  
     }
     //towCalib->Print("");
     // 
     AliCDBMetaData metaData;
     metaData.SetBeamPeriod(0);
     metaData.SetResponsible("Chiara Oppedisano");
     metaData.SetComment("AliZDCTowerCalib object");  
     //
     resTowCal = Store("Calib", "TowerCalib", towCalib, &metaData, 0, kTRUE);
   }
   else{ 
     // if entry exists it is still valid (=1 for all runs!)
     Log(Form(" Valid ZDC/Calib/TowerCalib object already existing in OCDB!!!"));
     resTowCal = kTRUE;
   }
   
   if(resTowCal==kFALSE) return 7;
   
   return 0;
}

//______________________________________________________________________________________________
UInt_t AliZDCPreprocessor::ProcessCalibData(Float_t beamEnergy)
{
  TList* daqSources = GetFileSources(kDAQ, "EMDENERGYCALIB");
  if(!daqSources){
    AliError(Form("No sources for CALIBRATION_EMD run %d !", fRun));
    return 8;
  }
  Log("\t List of DAQ sources for EMDENERGYCALIB id: "); daqSources->Print();
  //
  TIter iter2(daqSources);
  TObjString* source = 0;
  Int_t i=0;
  Bool_t resEnCal=kTRUE, resTowCal=kTRUE;
  
  while((source = dynamic_cast<TObjString*> (iter2.Next()))){
    TString stringEMDFileName = GetFile(kDAQ, "EMDENERGYCALIB", source->GetName());
    if(stringEMDFileName.Length() <= 0){
      Log(Form("No file from source %s!", source->GetName()));
      return 8;
    }
    const char* emdFileName = stringEMDFileName.Data();
    Log(Form("\t Getting file #%d: %s from %s\n",++i,emdFileName,source->GetName()));
    //
    // --- Initializing energy calibration object
    AliZDCEnCalib *eCalib = new AliZDCEnCalib("ZDC");
    // --- Reading file with calibration data
    if(emdFileName){
      FILE *file;
      if((file = fopen(emdFileName,"r")) == NULL){
        printf("Cannot open file %s \n",emdFileName);
        return 8;
      }
      Log(Form("File %s connected to process data from EM dissociation events", emdFileName));
      //
      Float_t fitValEMD[6];
      for(Int_t j=0; j<6; j++){        
        if(j<6){
          int iread = fscanf(file,"%f",&fitValEMD[j]);
          if(iread==0) AliDebug(3," Failing reading data from EMD calibration data file");
          if(fitValEMD[j]!=1.) eCalib->SetEnCalib(j, beamEnergy/fitValEMD[j]);
	  else eCalib->SetEnCalib(j, fitValEMD[j]);
        }
      }
      //
      fclose(file);
    }
    else{
      Log(Form("File %s not found", emdFileName));
      return 8;
    }
    //eCalib->Print("");
    // 
    AliCDBMetaData metaData;
    metaData.SetBeamPeriod(0);
    metaData.SetResponsible("Chiara Oppedisano");
    metaData.SetComment("Filling AliZDCEnCalib object");  
    //
    resEnCal = Store("Calib","EnergyCalib",eCalib, &metaData, 0, kTRUE);
    if(resEnCal==kFALSE) return 6;
  }
  delete daqSources; daqSources = 0;
  
  TList* daqSourcesH = GetFileSources(kDAQ, "EMDTOWERCALIB");
  if(!daqSourcesH){
    AliError(Form("No sources for CALIBRATION_EMD run %d !", fRun));
    return 9;
  }
  Log("\t List of DAQ sources for EMDTOWERCALIB id: "); daqSourcesH->Print();
  //
  TIter iter2H(daqSourcesH);
  TObjString* sourceH = 0;
  Int_t iH=0;
  while((sourceH = dynamic_cast<TObjString*> (iter2H.Next()))){
    TString stringtowEMDFileName = GetFile(kDAQ, "EMDTOWERCALIB", sourceH->GetName());
    if(stringtowEMDFileName.Length() <= 0){
      Log(Form("No file from source %s!", sourceH->GetName()));
      return 9;
    }
    const char * towEMDFileName = stringtowEMDFileName.Data();
    Log(Form("\t Getting file #%d: %s from source %s\n",++iH,towEMDFileName,sourceH->GetName()));
    // --- Initializing energy calibration object
    AliZDCTowerCalib *towCalib = new AliZDCTowerCalib("ZDC");
    // --- Reading file with calibration data
    if(towEMDFileName){
      FILE *file;
      if((file = fopen(towEMDFileName,"r")) == NULL){
        printf("Cannot open file %s \n",towEMDFileName);
        return 9;
      }
      //
      Float_t equalCoeff[4][5];
      for(Int_t j=0; j<4; j++){        
         for(Int_t k=0; k<5; k++){
           int leggi = fscanf(file,"%f",&equalCoeff[j][k]);
           if(leggi==0) AliDebug(3," Failing reading data from EMD calibration file");
           if(j==0)	 towCalib->SetZN1EqualCoeff(k, equalCoeff[j][k]);
           else if(j==1) towCalib->SetZP1EqualCoeff(k, equalCoeff[j][k]);
           else if(j==2) towCalib->SetZN2EqualCoeff(k, equalCoeff[j][k]);
           else if(j==3) towCalib->SetZP2EqualCoeff(k, equalCoeff[j][k]);  
         }
      }
      //
      fclose(file);
    }
    else{
      Log(Form("File %s not found", towEMDFileName));
      return 9;
    }
    //towCalib->Print("");
    // 
    AliCDBMetaData metaData;
    metaData.SetBeamPeriod(0);
    metaData.SetResponsible("Chiara Oppedisano");
    metaData.SetComment("Filling AliZDCTowerCalib object");  
    //
    resTowCal = Store("Calib","TowerCalib",towCalib, &metaData, 0, kTRUE);
    if(resTowCal==kFALSE) return 7;
  }
  delete daqSourcesH; daqSourcesH = 0;
  
     
  return 0;
}

//______________________________________________________________________________________________
UInt_t AliZDCPreprocessor::ProcessPedestalData()
{
  TList* daqSources = GetFileSources(kDAQ, "PEDESTALDATA");
  if(!daqSources){
    Log(Form("No source for STANDALONE_PEDESTAL run %d !", fRun));
    return 10;
  }
  if(daqSources->GetEntries()==0) return 10;
  Log("\t List of DAQ sources for PEDESTALDATA id: "); daqSources->Print();
  //
  TIter iter(daqSources);
  TObjString* source;
  Int_t i=0;
  Bool_t resPedCal=kTRUE, resPedHist=kTRUE;
  
  while((source = dynamic_cast<TObjString*> (iter.Next()))){
     TString stringPedFileName = GetFile(kDAQ, "PEDESTALDATA", source->GetName());
     if(stringPedFileName.Length() <= 0){
     	Log(Form("No PEDESTALDATA file from source %s!", source->GetName()));
        return 10;
     }
     const char* pedFileName = stringPedFileName.Data();
     Log(Form("\t Getting file #%d: %s from %s\n",++i,pedFileName,source->GetName()));
     //
     // --- Initializing pedestal calibration object
     AliZDCPedestals *pedCalib = new AliZDCPedestals("ZDC");
     // --- Reading file with pedestal calibration data
     // no. ADCch = (22 signal ch. + 2 reference PMs) * 2 gain chain = 48
     const Int_t knZDCch = 48;
     FILE *file;
     if((file = fopen(pedFileName,"r")) == NULL){
       printf("Cannot open file %s \n",pedFileName);
       return 10;
     }
     Log(Form("File %s connected to process pedestal data", pedFileName));
     Float_t pedVal[(3*knZDCch)][2];
     for(Int_t k=0; k<(3*knZDCch); k++){
        for(Int_t j=0; j<2; j++){
           int aleggi = fscanf(file,"%f",&pedVal[k][j]);
           if(aleggi==0) AliDebug(3," Failing reading data from pedestal file");
           //if(j==1) printf("pedVal[%d] -> %f, %f \n",k,pedVal[k][0],pedVal[k][1]);
        }
        if(k<knZDCch){
          pedCalib->SetMeanPed(k,pedVal[k][0]);
          pedCalib->SetMeanPedWidth(k,pedVal[k][1]);
        }
        else if(k>=knZDCch && k<(2*knZDCch)){
          pedCalib->SetOOTPed(k-knZDCch,pedVal[k][0]);
          pedCalib->SetOOTPedWidth(k-knZDCch,pedVal[k][1]);
        }
        else if(k>=(2*knZDCch) && k<(3*knZDCch)){
          pedCalib->SetPedCorrCoeff(k-(2*knZDCch),pedVal[k][0],pedVal[k][1]);
        }
     }
     fclose(file);
     //pedCalib->Print("");
     // 
     AliCDBMetaData metaData;
     metaData.SetBeamPeriod(0);
     metaData.SetResponsible("Chiara Oppedisano");
     metaData.SetComment("Filling AliZDCPedestals object");  
     //
     resPedCal = Store("Calib","Pedestals",pedCalib, &metaData, 0, kTRUE);
     if(resPedCal==kFALSE) return 11;
  }
  delete daqSources; daqSources = 0;
  
  TList* daqSourceH = GetFileSources(kDAQ, "PEDESTALHISTOS");
  if(!daqSourceH){
    Log(Form("No source for PEDESTALHISTOS id run %d !", fRun));
    return 12;
  }
  Log("\t List of DAQ sources for PEDESTALHISTOS id: "); daqSourceH->Print();
  //
  TIter iterH(daqSourceH);
  TObjString* sourceH = 0;
  Int_t iH=0;
  while((sourceH = dynamic_cast<TObjString*> (iterH.Next()))){
     TString stringPedFileName = GetFile(kDAQ, "PEDESTALHISTOS", sourceH->GetName());
     if(stringPedFileName.Length() <= 0){
     	Log(Form("No PEDESTALHISTOS file from source %s!", sourceH->GetName()));
        return 12;
     }
     const char* pedFileName = stringPedFileName.Data();
     Log(Form("\t Getting file #%d: %s from %s\n",++iH, pedFileName, sourceH->GetName()));
     resPedHist = StoreReferenceFile(pedFileName, "pedestalReference.root");
     if(resPedHist==kFALSE) return 13;
  }
  delete daqSourceH; daqSourceH=0;
  
  return 0;
}

//______________________________________________________________________________________________
UInt_t AliZDCPreprocessor::ProcessLaserData()
{
  TList* daqSources = GetFileSources(kDAQ, "LASERDATA");
  if(!daqSources){
    AliError(Form("No sources for STANDALONE_LASER run %d !", fRun));
    return 14;
  }
  if(daqSources->GetEntries()==0) return 14;
  Log("\t List of DAQ sources for LASERDATA id: "); daqSources->Print();
  //
  TIter iter2(daqSources);
  TObjString* source = 0;
  Int_t i=0;
  Bool_t resLaserCal=kTRUE, resLaserHist=kTRUE;
  
  while((source = dynamic_cast<TObjString*> (iter2.Next()))){
     TString stringLaserFileName = GetFile(kDAQ, "LASERDATA", source->GetName());
     if(stringLaserFileName.Length() <= 0){
       Log(Form("No LASER file from source %s!", source->GetName()));
       return 14;
     }
     const char* laserFileName = stringLaserFileName.Data();
     Log(Form("\t Getting file #%d: %s from %s\n",++i,laserFileName,source->GetName()));
     //
     // --- Initializing pedestal calibration object
     AliZDCLaserCalib *lCalib = new AliZDCLaserCalib("ZDC");
     // --- Reading file with pedestal calibration data
     if(laserFileName){
       FILE *file;
       if((file = fopen(laserFileName,"r")) == NULL){
         printf("Cannot open file %s \n",laserFileName);
         return 14;
       }
       Log(Form("File %s connected to process data from LASER events", laserFileName));
       //
       Float_t ivalRead[22][4]; 
       for(Int_t j=0; j<22; j++){
          for(Int_t k=0; k<4; k++){
            int aleggi = fscanf(file,"%f",&ivalRead[j][k]);
            if(aleggi==0) AliDebug(3," Failng reading data from laser file");
            //printf(" %d %1.0f  ",k, ivalRead[j][k]);
          }
          lCalib->SetDetector(j, (Int_t) ivalRead[j][0]);
          lCalib->SetSector(j, (Int_t) ivalRead[j][1]);
          lCalib->SetfPMValue(j, ivalRead[j][2]);
          lCalib->SetfPMWidth(j, ivalRead[j][3]);
       }
       fclose(file);
     }
     else{
       Log(Form("File %s not found", laserFileName));
       return 14;
     }
     //lCalib->Print("");
     // 
     AliCDBMetaData metaData;
     metaData.SetBeamPeriod(0);
     metaData.SetResponsible("Chiara Oppedisano");
     metaData.SetComment("Filling AliZDCLaserCalib object");  
     //
     resLaserCal = Store("Calib","LaserCalib",lCalib, &metaData, 0, kTRUE);
     if(resLaserCal==kFALSE) return 15;
  }
  delete daqSources; daqSources = 0;

  TList* daqSourceH = GetFileSources(kDAQ, "LASERHISTOS");
  if(!daqSourceH){
    AliError(Form("No sources for STANDALONE_LASER run %d !", fRun));
    return 16;
  }
  Log("\t List of DAQ sources for LASERHISTOS id: "); daqSourceH->Print();
  //
  TIter iter2H(daqSourceH);
  TObjString* sourceH = 0;
  Int_t iH=0;
  while((sourceH = dynamic_cast<TObjString*> (iter2H.Next()))){
     Log(Form("\t Getting file #%d\n",++iH));
     TString stringLaserFileName = GetFile(kDAQ, "LASERHISTOS", sourceH->GetName());
     if(stringLaserFileName.Length() <= 0){
       Log(Form("No LASER file from source %s!", sourceH->GetName()));
       return 16;
     }
     resLaserHist = StoreReferenceFile(stringLaserFileName.Data(), "laserReference.root");
     //
     if(resLaserHist==kFALSE) return 17;
  }
  delete daqSourceH; daqSourceH = 0;
  
  return 0;
}


//______________________________________________________________________________________________
UInt_t AliZDCPreprocessor::ProcessMBCalibData()
{
  TList* daqSources = GetFileSources(kDAQ, "MBCALIB");
  if(!daqSources){
    AliError(Form("No sources for CALIBRATION_MB run %d !", fRun));
    return 18;
  }
  if(daqSources->GetEntries()==0) return 18;
  Log("\t List of DAQ sources for MBCALIB id: "); daqSources->Print();
  //
  TIter iter2(daqSources);
  TObjString* source = 0;
  Int_t i=0;
  Bool_t resMBCal=kTRUE;
  
  while((source = dynamic_cast<TObjString*> (iter2.Next()))){
     TString stringMBFileName = GetFile(kDAQ, "MBCALIB", source->GetName());
     if(stringMBFileName.Length() <= 0){
       Log(Form("No MBCALIB file from source %s!", source->GetName()));
       return 18;
     }
     const char* mbFileName = stringMBFileName.Data();
     Log(Form("\t Getting file #%d: %s from %s\n",++i,mbFileName,source->GetName()));
     //
     // --- Initializing calibration object
     AliZDCMBCalib *mbCalib = new AliZDCMBCalib("ZDC");
     // --- Reading file with calibration data
     if(mbFileName){
       TFile * fileHistos = TFile::Open(mbFileName);
       Log(Form("File %s connected to process data from CALIBRATION_MB events", mbFileName));
       //
       fileHistos->cd();
       TH2F *hZDCvsZEM = (TH2F*)  fileHistos->Get("hZDCvsZEM");
       TH2F *hZDCCvsZEM = (TH2F*) fileHistos->Get("hZDCCvsZEM");
       TH2F *hZDCAvsZEM = (TH2F*) fileHistos->Get("hZDCAvsZEM");
       //
       mbCalib->SetZDCvsZEM(hZDCvsZEM);
       mbCalib->SetZDCCvsZEM(hZDCCvsZEM);
       mbCalib->SetZDCAvsZEM(hZDCAvsZEM);
       //
       //fileHistos->Close();
     }
     else{
       Log(Form("File %s not found", mbFileName));
       return 14;
     }
     // 
     AliCDBMetaData metaData;
     metaData.SetBeamPeriod(0);
     metaData.SetResponsible("Chiara Oppedisano");
     metaData.SetComment("Filling AliZDCMBCalib object");  
     //
     //mbCalib->Dump();
     //
     resMBCal = Store("Calib","MBCalib",mbCalib, &metaData, 0, kTRUE);
       printf(" here 1000\n");
     if(resMBCal==kFALSE) return 19;
  }
  delete daqSources; daqSources = 0;
  
  return 0;
}

//______________________________________________________________________________________________
UInt_t AliZDCPreprocessor::Process(TMap* dcsAliasMap)
{
 UInt_t resDCS = 0;
 UInt_t resChMap=0;
 UInt_t resEnergyCalib=0, resPedestalCalib=0, resLaserCalib=0, resMBCalib=0;

 // ************************* Process DCS data ****************************
 if(ProcessDCS()) resDCS = ProcessDCSData(dcsAliasMap);
  
 // ********************************* From DAQ ************************************

 const char* beamType = GetRunParameter("beamType");
 TString runType = GetRunType();
 Float_t beamEnergy = (Float_t)(((TString)GetRunParameter("beamEnergy")).Atof()); 
 printf("\t **** AliZDCPreprocessor -> runType %s, beamType %s,  beamEnergy %1.0f ****\n",
 	runType.Data(),beamType,beamEnergy);

 // ******************************************
 // ADC channel mapping
 // ******************************************
 resChMap = ProcessChMap();
 
 // ******************************************
 // Calibration param. for p-p data (all = 1)
 // ******************************************
 // NO ENERGY CALIBRATION -> coefficients set to 1.
 // Temp -> also inter-calibration coefficients are set to 1.
 if((strcmp(beamType,"p-p")==0) || (strcmp(beamType,"P-P")==0)) resEnergyCalib = ProcessppData();
 
 // *****************************************************
 // CALIBRATION_EMD -> Energy calibration and equalization
 // *****************************************************
 if((strcmp(beamType,"A-A")==0) && (runType.CompareTo("CALIBRATION_EMD")==0)) 
   resEnergyCalib =  ProcessCalibData(beamEnergy);
 
 // *****************************************************
 // STANDALONE_PEDESTALS -> Pedestal subtraction 
 // *****************************************************
 if(runType.CompareTo("STANDALONE_PEDESTAL")==0) resPedestalCalib = ProcessPedestalData();
 
 // *****************************************************
 // STANDALONE_LASER -> Signal stability and ageing 
 // *****************************************************
 else if(runType.CompareTo("STANDALONE_LASER")==0) resLaserCalib = ProcessLaserData();

 // *****************************************************
 // CALIBRATION_MB -> Signal stability and ageing 
 // *****************************************************
 else if(runType.CompareTo("CALIBRATION_MB")==0) resMBCalib = ProcessMBCalibData();
 
 if(resDCS!=0)  	      return resDCS;
 else if(resChMap!=0)	      return resChMap;
 else if(resEnergyCalib!=0)   return resEnergyCalib;
 else if(resPedestalCalib!=0) return resPedestalCalib;
 else if(resLaserCalib!=0)    return resLaserCalib;
 else if(resMBCalib!=0)       return resMBCalib;
 
 return 0;
  
}
