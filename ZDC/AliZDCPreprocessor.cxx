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
#include "AliZDCCalib.h"

/////////////////////////////////////////////////////////////////////
//								   //
// Class implementing ZDC pre-processor.			   //
// It takes data from DCS and passes it to the class AliZDCDataDCS //
// The class is then written to the CDB.			   //
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
  Log("\t List of sources "); daqSource->Print();
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
             fscanf(file,"%d",&readMap[j][k]);
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
    Log("\t AliZDCPreprocessor -> WARNING! No CDB entry for ch. mapping\n");
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
    Log("\t AliZDCPreprocessor -> A new entry ZDC/Calib/ChMap will be created");
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
    Log("\t AliZDCPreprocessor -> ZDC/Calib/ChMap entry in OCDB is valid and won't be updated\n");
    res = kTRUE;
  }

  
  return res;

}

//______________________________________________________________________________________________
UInt_t AliZDCPreprocessor::Process(TMap* dcsAliasMap)
{
  // *************** From DCS ******************
  // Fills data into a AliZDCDataDCS object
  if(!dcsAliasMap) return 1;
  printf("Processing data from DCS\n");
  
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
  dcsAliasMap->Print("");

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
  
 // *************** From DAQ ******************
 Bool_t resChMap=kTRUE, resPedCal=kTRUE, resLaserCal=kTRUE, resECal=kTRUE;
 // 
 const char* beamType = GetRunParameter("beamType");
 TString runType = GetRunType();
 printf("\n\t AliZDCPreprocessor -> beamType %s\n",beamType);
 printf("\t AliZDCPreprocessor -> runType  %s\n\n",runType.Data());
 
 // -------------- p-p data -------------
 if(strcmp(beamType,"p-p")==0){
   
   // --- Cheking if there is already the entry in the OCDB
   AliCDBEntry *cdbEntry = GetFromOCDB("Calib", "EMDCalib");
   if(!cdbEntry){   
     printf("\t AliZDCPreprocessor -> ZDC/Calib/EMDCalib entry will be created\n");
     // --- Initializing calibration object
     AliZDCCalib *eCalib = new AliZDCCalib("ZDC");
     //
     for(Int_t j=0; j<6; j++) eCalib->SetEnCalib(j,1.);
     for(Int_t j=0; j<5; j++){  
        eCalib->SetZN1EqualCoeff(j, 1.);
        eCalib->SetZP1EqualCoeff(j, 1.);
        eCalib->SetZN2EqualCoeff(j, 1.);
        eCalib->SetZP2EqualCoeff(j, 1.);  
     }
     //eCalib->Print("");
     // 
     AliCDBMetaData metaData;
     metaData.SetBeamPeriod(0);
     metaData.SetResponsible("Chiara Oppedisano");
     metaData.SetComment("AliZDCCalib object");  
     //
     resECal = Store("Calib","EMDCalib",eCalib, &metaData, 0, 1);
   }
   else{
     printf("\t AliZDCPreprocessor -> ZDC/Calib/EMDCalib object already existing in OCDB!!!\n");
     resECal = kTRUE;
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
  TList* daqSources = GetFileSources(kDAQ, "PEDESTALS");
  if(!daqSources){
    Log(Form("No source for STANDALONE_PEDESTAL run %d !", fRun));
    return 1;
  }
  Log("\t List of sources for STANDALONE_PEDESTAL");
  daqSources->Print();
  //
  TIter iter(daqSources);
  TObjString* source = 0;
  Int_t i=0;
  while((source = dynamic_cast<TObjString*> (iter.Next()))){
       Log(Form("\n\t Getting file #%d\n",++i));
       TString stringPedFileName = GetFile(kDAQ, "PEDESTALS", source->GetName());
       if(stringPedFileName.Length() <= 0){
          Log(Form("No PEDESTAL file from source %s!", source->GetName()));
	  return 1;
       }
       // --- Initializing pedestal calibration object
       AliZDCPedestals *pedCalib = new AliZDCPedestals("ZDC");
       // --- Reading file with pedestal calibration data
       const char* pedFileName = stringPedFileName.Data();
       // no. ADCch = (22 signal ch. + 2 reference PMs) * 2 gain chain = 48
       const Int_t knZDCch = 48;
       if(pedFileName){
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
       }
       else{
          Log(Form("File %s not found", pedFileName));
          return 1;
       }
       //
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
 }
 // *****************************************************
 // [b] STANDALONE_LASER EVENTS -> Signal stability
 // *****************************************************
 else if(runType=="STANDALONE_LASER"){
  TList* daqSources = GetFileSources(kDAQ, "LASER");
  if(!daqSources){
    AliError(Form("No sources for STANDALONE_LASER run %d !", fRun));
    return 1;
  }
  Log("\t List of sources for STANDALONE_LASER");
  daqSources->Print();
  //
  TIter iter2(daqSources);
  TObjString* source = 0;
  Int_t i=0;
  while((source = dynamic_cast<TObjString*> (iter2.Next()))){
       Log(Form("\n\t Getting file #%d\n",++i));
       TString stringLASERFileName = GetFile(kDAQ, "LASER", source->GetName());
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
  
 }
 // *****************************************************
 // [c] EMD EVENTS -> Energy calibration and equalization
 // *****************************************************
 else if(runType=="CALIBRATION_EMD" && strcmp(beamType,"Pb-Pb")==0){
  TList* daqSources = GetFileSources(kDAQ, "EMDCALIB");
  if(!daqSources){
    AliError(Form("No sources for CALIBRATION_EMD run %d !", fRun));
    return 1;
  }
  Log("\t List of sources for CALIBRATION_EMD");
  daqSources->Print();
  //
  TIter iter2(daqSources);
  TObjString* source = 0;
  Int_t i=0;
  while((source = dynamic_cast<TObjString*> (iter2.Next()))){
      Log(Form("\n\t Getting file #%d\n",++i));
      TString stringEMDFileName = GetFile(kDAQ, "EMDCALIB", source->GetName());
      if(stringEMDFileName.Length() <= 0){
        Log(Form("No EMDCALIB file from source %s!", source->GetName()));
        return 1;
      }
      // --- Initializing pedestal calibration object
      AliZDCCalib *eCalib = new AliZDCCalib("ZDC");
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
        Float_t fitValEMD[6]; Float_t equalCoeff[4][5];
        Float_t calibVal[4];
        for(Int_t j=0; j<10; j++){	   
          if(j<6){
            fscanf(file,"%f",&fitValEMD[j]);
            if(j<4){
              calibVal[j] = 2760./fitValEMD[j];
              eCalib->SetEnCalib(j,calibVal[j]);
            }
            else eCalib->SetEnCalib(j,fitValEMD[j]);
          }
          else{
            for(Int_t k=0; k<5; k++){
               fscanf(file,"%f",&equalCoeff[j][k]);
               if(j==6)      eCalib->SetZN1EqualCoeff(k, equalCoeff[j][k]);
               else if(j==7) eCalib->SetZP1EqualCoeff(k, equalCoeff[j][k]);
               else if(j==8) eCalib->SetZN2EqualCoeff(k, equalCoeff[j][k]);
               else if(j==9) eCalib->SetZP2EqualCoeff(k, equalCoeff[j][k]);  
            }
          }
        }
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
      metaData.SetComment("Filling AliZDCCalib object");  
      //
      resECal = Store("Calib","EMDCalib",eCalib, &metaData, 0, 1);
  }
 }


  // note that the parameters are returned as character strings!
  const char* nEvents = GetRunParameter("totalEvents");
  if(nEvents) Log(Form("Number of events for run %d: %s",fRun, nEvents));
  else Log(Form("Number of events not put in logbook!"));
 
  UInt_t result = 0;
  if(resDCSRef==kFALSE || resultAl==kFALSE || resPedCal==kFALSE ||
     resLaserCal==kFALSE || resECal==kFALSE || resChMap==kFALSE){
    if(resDCSRef == kFALSE)        result = 1;
    else if(resultAl == kFALSE)    result = 2;
    else if(resChMap == kFALSE)    result = 3;
    else if(resPedCal == kFALSE)   result = 4;
    else if(resLaserCal == kFALSE) result = 5;
    else if(resECal == kFALSE)     result = 6;
  }
  
  return result;
  
}
