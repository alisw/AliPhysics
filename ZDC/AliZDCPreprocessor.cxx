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

/////////////////////////////////////////////////////////////////////
//								   //
// Class implementing Shuttle ZDC pre-processor.	           //
// It takes data from DCS and DAQ and writes calibration objects   //
// in the OCDB and reference values/histos in the ReferenceData.   //
//								   //
/////////////////////////////////////////////////////////////////////

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
void AliZDCPreprocessor::Initialize(Int_t run, UInt_t startTime,
	UInt_t endTime)
{
  // Creates AliZDCDataDCS object

  AliPreprocessor::Initialize(run, startTime, endTime);

  Log(Form("\n\tRun %d \n\tStartTime %s \n\tEndTime %s", run,
  	  TTimeStamp(startTime).AsString(),
  	  TTimeStamp(endTime).AsString()));

  fRun = run;
  fStartTime = startTime;
  fEndTime = endTime;

  fData = new AliZDCDataDCS(fRun, fStartTime, fEndTime);
}

//______________________________________________________________________________________________
UInt_t AliZDCPreprocessor::ProcessChMap(TString runType)
{ 
  // Reading the file for mapping from FXS
  TList* daqSource = GetFileSources(kDAQ, runType.Data());
  if(!daqSource){
    AliError(Form("No sources run %d for run type %s!", fRun, runType.Data()));
    return 1;
  }
  Log("\t List of DAQ sources for mapping "); daqSource->Print();
  //
  TIter iter(daqSource);
  TObjString* source = 0;
  Int_t isou = 0;
  Int_t res = 999;
  Int_t readMap[48][6]; 
  //
  while((source = dynamic_cast<TObjString*> (iter.Next()))){
     Log(Form("\n\t Getting file #%d\n",++isou));
     TString fileName = "ZDCChMapping.dat";

     if(fileName.Length() <= 0){
       Log(Form("No file from source %s!", source->GetName()));
       return 1;
     }
     // --- Reading file with calibration data
     //const char* fname = fileName.Data();
     if(fileName){
       FILE *file;
       if((file = fopen(fileName,"r")) == NULL){
	 printf("Cannot open file %s \n",fileName.Data());
         return 1;
       }
       Log(Form("File %s connected to process data for ADC mapping", fileName.Data()));
       //
       for(Int_t j=0; j<48; j++){	  
           for(Int_t k=0; k<6; k++){
             int read = fscanf(file,"%d",&readMap[j][k]);
	     if(read == 0){
	       printf("\t WARNING! Failing in reading datum [%d][%d]\n",j,k);
	     }
           }
       }
       fclose(file);
     }
     else{
       Log(Form("File %s not found", fileName.Data()));
       return 1;
     }
  }
  delete daqSource; daqSource=0;
  
  // Store the currently read map ONLY IF it is different
  // from the entry in the OCDB
  Bool_t updateOCDB = kFALSE;
  
  AliCDBEntry *cdbEntry = GetFromOCDB("Calib","ChMap");
  if(!cdbEntry){
    Log(" No existing CDB entry for ADC mapping");
    updateOCDB = kTRUE;
  }
  else{
    AliZDCChMap *chMap = (AliZDCChMap*) cdbEntry->GetObject();
    for(Int_t i=0; i<48; i++){
      if(  (readMap[i][1] == chMap->GetADCModule(i)) 
        && (readMap[i][2] == chMap->GetADCChannel(i)) 
	&& (readMap[i][4] == chMap->GetDetector(i)) 
	&& (readMap[i][5] == chMap->GetSector(i))){
	 updateOCDB = kFALSE;
      }
      else updateOCDB = kTRUE;
    }
  }
  //
  if(updateOCDB==kTRUE){
    Log(" A new entry ZDC/Calib/ChMap will be created");
    //
    // --- Initializing mapping calibration object
    AliZDCChMap *mapCalib = new AliZDCChMap("ZDC");
    // Writing channel map in the OCDB
    for(Int_t k=0; k<48; k++){
      mapCalib->SetADCModule(k,readMap[k][1]);
      mapCalib->SetADCChannel(k,readMap[k][2]);
      mapCalib->SetDetector(k,readMap[k][4]);
      mapCalib->SetSector(k,readMap[k][5]);
    }
    //mapCalib->Print("");
    // 
    AliCDBMetaData metaData;
    metaData.SetBeamPeriod(0);
    metaData.SetResponsible("Chiara Oppedisano");
    metaData.SetComment("Filling AliZDCChMap object");  
    //
    res = Store("Calib","ChMap",mapCalib, &metaData, 0, 1);
  }
  else{
    Log(" ZDC/Calib/ChMap entry in OCDB is valid and won't be updated");
    res = kTRUE;
  }

  
  return res;

}

//______________________________________________________________________________________________
UInt_t AliZDCPreprocessor::Process(TMap* dcsAliasMap)
{
 // ********************************* From DCS ************************************
 // Fills data into a AliZDCDataDCS object
 if(!dcsAliasMap) return 1;
 Log(Form("Processing data from DCS"));
 
 // The processing of the DCS input data is forwarded to AliZDCDataDCS
 Float_t dcsValues[28]; // DCSAliases=28
 fData->SetCalibData(dcsValues);
 fData->ProcessData(*dcsAliasMap);
 // Store DCS data for reference
 AliCDBMetaData metadata;
 metadata.SetResponsible("Chiara Oppedisano");
 metadata.SetComment("DCS data for ZDC");
 Bool_t resDCSRef = kTRUE;
 resDCSRef = StoreReferenceData("DCS","Data",fData,&metadata);
 //dcsAliasMap->Print("");

 // --- Writing ZDC table positions into alignment object
 TClonesArray *array = new TClonesArray("AliAlignObjParams",10);
 TClonesArray &alobj = *array;
 AliAlignObjParams a;
 Double_t dx=0., dz=0., dpsi=0., dtheta=0., dphi=0.;
 // Vertical table position in mm from DCS
 Double_t dyZN1 = (Double_t) (dcsValues[0]/10.);
 Double_t dyZP1 = (Double_t) (dcsValues[1]/10.);
 Double_t dyZN2 = (Double_t) (dcsValues[2]/10.);
 Double_t dyZP2 = (Double_t) (dcsValues[3]/10.);
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
 AliCDBMetaData md;
 md.SetResponsible("Chiara Oppedisano");
 md.SetComment("Alignment object for ZDC");
 Bool_t resultAl = kTRUE;
 resultAl = Store("Align","Data", array, &md, 0, 0);
  
 // ********************************* From DAQ ************************************
 Bool_t resChMap=kTRUE, resPedCal=kTRUE, resPedHist=kTRUE;
 Bool_t resLaserCal=kTRUE, resLaserHist=kTRUE;
 Bool_t resEnCal=kTRUE, resTowCal=kTRUE;
 // 
 const char* beamType = GetRunParameter("beamType");
 TString runType = GetRunType();
 printf("\n\t AliZDCPreprocessor -> beamType %s\n",beamType);
 printf("\t AliZDCPreprocessor -> runType  %s\n\n",runType.Data());
 
 // -------------- p-p data -------------
 // NO ENERGY CALIBRATION -> coefficients set to 1.
 // Temp -> also inter-calibration coefficients are set to 1.
 if(strcmp(beamType,"p-p")==0){
 
   // *********** Energy calibration
   // --- Cheking if there is already the entry in the OCDB
   AliCDBEntry *cdbEnEntry = GetFromOCDB("Calib", "EnergyCalib");
   if(!cdbEnEntry){   
     Log(Form(" ZDC/Calib/EnergyCalib entry will be created"));
     // --- Initializing calibration object
     // 
     AliCDBMetaData metaData;
     metaData.SetBeamPeriod(0);
     metaData.SetResponsible("Chiara Oppedisano");
     //
     AliZDCEnCalib *eCalib = new AliZDCEnCalib("ZDC");
     for(Int_t j=0; j<6; j++) eCalib->SetEnCalib(j,1.);
     metaData.SetComment("AliZDCEnCalib object");  
     //eCalib->Print("");
     resEnCal = Store("Calib", "EnergyCalib", eCalib, &metaData, 0, 1);
   }
   else{ 
     // if entry exists it is still valid (=1 for all runs!)
     Log(Form(" Valid ZDC/Calib/EnergyCalib object already existing in OCDB!!!"));
     resEnCal = kTRUE;
   }
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
     resTowCal = Store("Calib", "TowerCalib", towCalib, &metaData, 0, 1);
   }
   else{ 
     // if entry exists it is still valid (=1 for all runs!)
     Log(Form(" Valid ZDC/Calib/TowerCalib object already existing in OCDB!!!"));
     resTowCal = kTRUE;
   }
   
   // *********** Tower relative calibration
   // --- Cheking if there is already the entry in the OCDB
   AliCDBEntry *cdbSEntry = GetFromOCDB("Calib", "TowerCalib");
   if(!cdbSEntry){   
     Log(Form(" ZDC/Calib/TowerCalib entry will be created"));
     // --- Initializing calibration object
     AliZDCTowerCalib *towCalib = new AliZDCTowerCalib("ZDC");
     //
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
     resTowCal = Store("Calib", "TowerCalib", towCalib, &metaData, 0, 1);
   }
   else{ 
     // if entry exists it is still valid (=1 for all runs!)
     Log(Form(" Valid ZDC/Calib/TowerCalib object already existing in OCDB!!!"));
     resEnCal = kTRUE;
   }
 }
 
 // ******************************************
 //   ZDC ADC channel mapping
 // ******************************************
 resChMap = ProcessChMap(runType);
 // 
 // *****************************************************
 // [a] PEDESTALS -> Pedestal subtraction
 // *****************************************************
 // 
 if(runType=="STANDALONE_PEDESTAL"){
  TList* daqSources = GetFileSources(kDAQ, "PEDESTALDATA");
  if(!daqSources){
    Log(Form("No source for STANDALONE_PEDESTAL run %d !", fRun));
    return 1;
  }
  Log("\t List of DAQ sources for STANDALONE_PEDESTAL run");
  daqSources->Print();
  //
  TIter iter(daqSources);
  TObjString* source = 0;
  Int_t i=0;
  while((source = dynamic_cast<TObjString*> (iter.Next()))){
     Log(Form("\n\t Getting file #%d\n",++i));
     TString stringPedFileName = GetFile(kDAQ, "PEDESTALDATA", source->GetName());
     if(stringPedFileName.Length() <= 0){
     	Log(Form("No PEDESTAL file from source %s!", source->GetName()));
        return 1;
     }
     // calibration data
       // --- Initializing pedestal calibration object
       AliZDCPedestals *pedCalib = new AliZDCPedestals("ZDC");
       // --- Reading file with pedestal calibration data
       const char* pedFileName = stringPedFileName.Data();
       // no. ADCch = (22 signal ch. + 2 reference PMs) * 2 gain chain = 48
       const Int_t knZDCch = 48;

       FILE *file;
       if((file = fopen(pedFileName,"r")) == NULL){
         printf("Cannot open file %s \n",pedFileName);
         return 1;
       }
       Log(Form("File %s connected to process pedestal data", pedFileName));
       Float_t pedVal[(2*knZDCch)][2];
       for(Int_t k=0; k<(2*knZDCch); k++){
          for(Int_t j=0; j<2; j++){
             fscanf(file,"%f",&pedVal[k][j]);
             //if(j==1) printf("pedVal[%d] -> %f, %f \n",k,pedVal[k][0],pedVal[k][1]);
          }
          if(k<knZDCch){
            pedCalib->SetMeanPed(k,pedVal[k][0]);
            pedCalib->SetMeanPedWidth(i,pedVal[k][1]);
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
       resPedCal = Store("Calib","Pedestals",pedCalib, &metaData, 0, 1);
  }
  delete daqSources; daqSources = 0;
  //
  TList* daqSourceH = GetFileSources(kDAQ, "PEDHISTOS");
  if(!daqSourceH){
    Log(Form("No source for STANDALONE_PEDESTAL run %d !", fRun));
    return 1;
  }
  Log("\t List of DAQ sources for STANDALONE_PEDESTAL run");
  daqSourceH->Print();
  //
  TIter iterH(daqSourceH);
  TObjString* sourceH = 0;
  Int_t iH=0;
  while((sourceH = dynamic_cast<TObjString*> (iterH.Next()))){
     Log(Form("\n\t Getting file #%d\n",++iH));
     TString stringPedFileName = GetFile(kDAQ, "PEDHISTOS", sourceH->GetName());
     if(stringPedFileName.Length() <= 0){
     	Log(Form("No PEDESTAL file from source %s!", sourceH->GetName()));
        return 1;
     }
     TFile *histoFile = TFile::Open(stringPedFileName.Data());
     //
     AliCDBMetaData metadata1;
     metadata1.SetResponsible("Chiara Oppedisano");
     metadata1.SetComment("Pedestal histos");
     resPedHist = StoreReferenceData("Histos","Pedestal", histoFile, &metadata1);
  }
  delete daqSourceH; daqSourceH = 0;
 }
 // *****************************************************
 // [b] STANDALONE_LASER EVENTS -> Signal stability
 // *****************************************************
 else if(runType=="STANDALONE_LASER"){
  TList* daqSources = GetFileSources(kDAQ, "LASERDATA");
  if(!daqSources){
    AliError(Form("No sources for STANDALONE_LASER run %d !", fRun));
    return 1;
  }
  Log("\t List of DAQ sources for STANDALONE_LASER run");
  daqSources->Print();
  //
  TIter iter2(daqSources);
  TObjString* source = 0;
  Int_t i=0;
  while((source = dynamic_cast<TObjString*> (iter2.Next()))){
     Log(Form("\n\t Getting file #%d\n",++i));
     TString stringLASERFileName = GetFile(kDAQ, "LASERDATA", source->GetName());
     if(stringLASERFileName.Length() <= 0){
       Log(Form("No LASER file from source %s!", source->GetName()));
       return 1;
     }
     // --- Initializing pedestal calibration object
     AliZDCLaserCalib *lCalib = new AliZDCLaserCalib("ZDC");
     // --- Reading file with pedestal calibration data
     const char* laserFileName = stringLASERFileName.Data();
     if(laserFileName){
       FILE *file;
       if((file = fopen(laserFileName,"r")) == NULL){
         printf("Cannot open file %s \n",laserFileName);
         return 1;
       }
       Log(Form("File %s connected to process data from LASER events", laserFileName));
       //
       Float_t ivalRead[22][4]; 
       for(Int_t j=0; j<22; j++){
     	  for(Int_t k=0; k<4; k++){
     	    fscanf(file,"%f",&ivalRead[j][k]);
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
       return 1;
     }
     //lCalib->Print("");
     // 
     AliCDBMetaData metaData;
     metaData.SetBeamPeriod(0);
     metaData.SetResponsible("Chiara Oppedisano");
     metaData.SetComment("Filling AliZDCLaserCalib object");  
     //
     resLaserCal = Store("Calib","LaserCalib",lCalib, &metaData, 0, 1);
  }
  delete daqSources; daqSources = 0;
  //
  TList* daqSourceH = GetFileSources(kDAQ, "LASERHISTOS");
  if(!daqSourceH){
    AliError(Form("No sources for STANDALONE_LASER run %d !", fRun));
    return 1;
  }
  Log("\t List of DAQ sources for STANDALONE_LASER run");
  daqSourceH->Print();
  //
  TIter iter2H(daqSourceH);
  TObjString* sourceH = 0;
  Int_t iH=0;
  while((sourceH = dynamic_cast<TObjString*> (iter2H.Next()))){
     Log(Form("\n\t Getting file #%d\n",++iH));
     TString stringLASERFileName = GetFile(kDAQ, "LASERHISTOS", sourceH->GetName());
     if(stringLASERFileName.Length() <= 0){
       Log(Form("No LASER file from source %s!", sourceH->GetName()));
       return 1;
     }
     TFile *histoFile = TFile::Open(stringLASERFileName.Data());
     //
     AliCDBMetaData metadata2;
     metadata2.SetResponsible("Chiara Oppedisano");
     metadata2.SetComment("Laser run histos");
     resPedHist = StoreReferenceData("Histos", "Laser", histoFile, &metadata2);
  }
  delete daqSourceH; daqSourceH = 0;
 }
 // *****************************************************
 // [c] EMD EVENTS -> Energy calibration and equalization
 // *****************************************************
 else if(runType=="CALIBRATION_EMD" && strcmp(beamType,"Pb-Pb")==0){
  TList* daqSources = GetFileSources(kDAQ, "ENERGYCALIB");
  if(!daqSources){
    AliError(Form("No sources for CALIBRATION_EMD run %d !", fRun));
    return 1;
  }
  Log("\t List of DAQ sources for CALIBRATION_EMD run");
  daqSources->Print();
  //
  TIter iter2(daqSources);
  TObjString* source = 0;
  Int_t i=0;
  while((source = dynamic_cast<TObjString*> (iter2.Next()))){
    Log(Form("\n\t Getting file #%d\n",++i));
    TString stringEMDFileName = GetFile(kDAQ, "ENERGYCALIB", source->GetName());
    if(stringEMDFileName.Length() <= 0){
      Log(Form("No file from source %s!", source->GetName()));
      return 1;
    }
    // --- Initializing energy calibration object
    AliZDCEnCalib *eCalib = new AliZDCEnCalib("ZDC");
    // --- Reading file with pedestal calibration data
    const char* emdFileName = stringEMDFileName.Data();
    if(emdFileName){
      FILE *file;
      if((file = fopen(emdFileName,"r")) == NULL){
    	printf("Cannot open file %s \n",emdFileName);
    	return 1;
      }
      Log(Form("File %s connected to process data from EM dissociation events", emdFileName));
      //
      Float_t fitValEMD[6];
      for(Int_t j=0; j<6; j++){        
    	if(j<6){
    	  fscanf(file,"%f",&fitValEMD[j]);
    	  eCalib->SetEnCalib(j,fitValEMD[j]);
    	}
      }
      //
      fclose(file);
    }
    else{
      Log(Form("File %s not found", emdFileName));
      return 1;
    }
    //eCalib->Print("");
    // 
    AliCDBMetaData metaData;
    metaData.SetBeamPeriod(0);
    metaData.SetResponsible("Chiara Oppedisano");
    metaData.SetComment("Filling AliZDCEnCalib object");  
    //
    resEnCal = Store("Calib","EnergyCalib",eCalib, &metaData, 0, 1);
  }
  delete daqSources; daqSources = 0;
  //
  TList* daqSourcesH = GetFileSources(kDAQ, "TOWERCALIB");
  if(!daqSourcesH){
    AliError(Form("No sources for CALIBRATION_EMD run %d !", fRun));
    return 1;
  }
  Log("\t List of DAQ sources for CALIBRATION_EMD run");
  daqSourcesH->Print();
  //
  TIter iter2H(daqSourcesH);
  TObjString* sourceH = 0;
  Int_t iH=0;
  while((sourceH = dynamic_cast<TObjString*> (iter2H.Next()))){
    Log(Form("\n\t Getting file #%d\n",++iH));
    TString stringEMDFileName = GetFile(kDAQ, "TOWERCALIB", sourceH->GetName());
    if(stringEMDFileName.Length() <= 0){
      Log(Form("No file from source %s!", sourceH->GetName()));
      return 1;
    }
    // --- Initializing energy calibration object
    AliZDCTowerCalib *towCalib = new AliZDCTowerCalib("ZDC");
    // --- Reading file with pedestal calibration data
    const char* emdFileName = stringEMDFileName.Data();
    if(emdFileName){
      FILE *file;
      if((file = fopen(emdFileName,"r")) == NULL){
    	printf("Cannot open file %s \n",emdFileName);
    	return 1;
      }
      Log(Form("File %s connected to process data from EM dissociation events", emdFileName));
      //
      Float_t equalCoeff[4][5];
      for(Int_t j=0; j<4; j++){        
    	 for(Int_t k=0; k<5; k++){
    	    fscanf(file,"%f",&equalCoeff[j][k]);
    	    if(j==0)	  towCalib->SetZN1EqualCoeff(k, equalCoeff[j][k]);
    	    else if(j==1) towCalib->SetZP1EqualCoeff(k, equalCoeff[j][k]);
    	    else if(j==2) towCalib->SetZN2EqualCoeff(k, equalCoeff[j][k]);
    	    else if(j==3) towCalib->SetZP2EqualCoeff(k, equalCoeff[j][k]);  
    	 }
      }
      //
      fclose(file);
    }
    else{
      Log(Form("File %s not found", emdFileName));
      return 1;
    }
    //towCalib->Print("");
    // 
    AliCDBMetaData metaData;
    metaData.SetBeamPeriod(0);
    metaData.SetResponsible("Chiara Oppedisano");
    metaData.SetComment("Filling AliZDCTowerCalib object");  
    //
    resTowCal = Store("Calib","TowerCalib",towCalib, &metaData, 0, 1);
  }
  delete daqSources; daqSources = 0;
 }//Pb-Pb if


  // note that the parameters are returned as character strings!
  const char* nEvents = GetRunParameter("totalEvents");
  if(nEvents) Log(Form("Number of events for run %d: %s",fRun, nEvents));
  else Log(Form("Number of events not put in logbook!"));
 
  UInt_t result = 0;
  if(resDCSRef==kFALSE || resultAl==kFALSE || resPedCal==kFALSE || resPedHist==kFALSE ||
     resLaserCal==kFALSE || resLaserHist==kFALSE || resEnCal==kFALSE || resTowCal || 
     resChMap==kFALSE){
      if(resDCSRef == kFALSE)         result = 1;
      else if(resultAl == kFALSE)     result = 2;
      else if(resChMap == kFALSE)     result = 3;
      else if(resPedCal == kFALSE)    result = 4;
      else if(resPedHist == kFALSE)   result = 5;
      else if(resLaserCal == kFALSE)  result = 6;
      else if(resLaserHist == kFALSE) result = 7;
      else if(resEnCal == kFALSE)     result = 8;
      else if(resTowCal == kFALSE)    result = 9;
  }
  
  return result;
  
}
