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
  AddRunType("STANDALONE_PEDESTAL");
  AddRunType("STANDALONE_LASER");
  AddRunType("STANDALONE_EMD");
  AddRunType("STANDALONE_COSMIC");
  AddRunType("STANDALONE_BC");
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
  // Writing channel map in the OCDB
  TList* daqSource=0;
  if(runType.CompareTo("STANDALONE_PEDESTAL"))
    daqSource = GetFileSources(kDAQ, "PEDESTALS");
  else if(runType.CompareTo("STANDALONE_LASER"))
    daqSource = GetFileSources(kDAQ, "LASER");
  else if(runType.CompareTo("STANDALONE_EMD"))
    daqSource = GetFileSources(kDAQ, "EMDCALIB");
  else if(runType.CompareTo("STANDALONE_COSMIC"))
    daqSource = GetFileSources(kDAQ, "COSMICS");
  else if(runType.CompareTo("STANDALONE_BC"))
    daqSource = GetFileSources(kDAQ, "BC");
  else if(runType.CompareTo("PHYSICS"))
    daqSource = GetFileSources(kDAQ, "PHYSICS");
  
  if(!daqSource){
    AliError(Form("No sources run %d for run type %s!", fRun, runType));
    return 1;
  }
  Log("\t List of sources "); daqSource->Print();
  //
  TIter iter(daqSource);
  TObjString* source = 0;
  Int_t isou=0;
  Int_t res=999;
  while((source = dynamic_cast<TObjString*> (iter.Next()))){
     Log(Form("\n\t Getting file #%d\n",++isou));
     TString fileName;
     if(runType.CompareTo("STANDALONE_PEDESTAL"))
      fileName = GetFile(kDAQ, "PEDESTALS", source->GetName());
     else if(runType.CompareTo("STANDALONE_LASER"))
      fileName = GetFile(kDAQ, "LASER", source->GetName());
     else if(runType.CompareTo("STANDALONE_EMD"))
      fileName = GetFile(kDAQ, "EMDCALIB", source->GetName());
     else if(runType.CompareTo("STANDALONE_COSMIC"))
      fileName = GetFile(kDAQ, "COSMICS", source->GetName());
     else if(runType.CompareTo("STANDALONE_BC"))
      fileName = GetFile(kDAQ, "BC", source->GetName());
     else if(runType.CompareTo("PHYSICS"))
      fileName = GetFile(kDAQ, "PHYSICS", source->GetName());

     if(fileName.Length() <= 0){
       Log(Form("No file from source %s!", source->GetName()));
       return 1;
     }
     // --- Initializing pedestal calibration object
     AliZDCChMap *mapCalib = new AliZDCChMap("ZDC");
     // --- Reading file with pedestal calibration data
     const char* fname = fileName.Data();
     if(fname){
       FILE *file;
       if((file = fopen(fname,"r")) == NULL){
	 printf("Cannot open file %s \n",fname);
         return 1;
       }
       Log(Form("File %s connected to process data for ADC mapping", fname));
       //
       Int_t chMap[48][6]; 
       for(Int_t j=0; j<48; j++){	  
           for(Int_t k=0; k<6; k++){
             fscanf(file,"%d",&chMap[j][k]);
           }
	   mapCalib->SetADCModule(j,chMap[j][1]);
	   mapCalib->SetADCChannel(j,chMap[j][2]);
	   mapCalib->SetDetector(j,chMap[j][4]);
	   mapCalib->SetSector(j,chMap[j][5]);
       }
       fclose(file);
     }
     else{
       Log(Form("File %s not found", fname));
       return 1;
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
  delete daqSource; daqSource=0;
  
  return res;

}

//______________________________________________________________________________________________
UInt_t AliZDCPreprocessor::Process(TMap* dcsAliasMap)
{
  // *************** From DCS ******************
  // Fills data into a AliZDCDataDCS object
  if(!dcsAliasMap) return 1;

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
  //
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
  const char *n1ZDC="ZDC/NeutronZDC1";
  const char *p1ZDC="ZDC/ProtonZDC1";
  const char *n2ZDC="ZDC/NeutronZDC2";
  const char *p2ZDC="ZDC/ProtonZDC2";
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

if(strcmp(beamType,"p-p")==0){

   // --- Initializing pedestal calibration object
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
         Float_t pedVal[(3*knZDCch)][3];
         for(Int_t k=0; k<(3*knZDCch); k++){
            for(Int_t j=0; j<3; j++){
               fscanf(file,"%f",&pedVal[k][j]);
	       //if(j==1) printf("pedVal[%d] -> %f, %f \n",k,pedVal[k][0],pedVal[k][1]);
            }
            if(k<knZDCch){
              pedCalib->SetMeanPed(k,pedVal[k][1]);
              pedCalib->SetMeanPedWidth(i,pedVal[k][2]);
            }
            else if(k>=knZDCch && k<(2*knZDCch)){
              pedCalib->SetOOTPed(k-knZDCch,pedVal[k][1]);
              pedCalib->SetOOTPedWidth(k-knZDCch,pedVal[k][2]);
            }
            else if(k>=(2*knZDCch) && k<(3*knZDCch)){
              pedCalib->SetPedCorrCoeff(k-(2*knZDCch),pedVal[k][1],pedVal[k][2]);
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
	 Float_t ivalRead[4][2]; 
    	 for(Int_t j=0; j<4; j++){
            for(Int_t k=0; k<2; k++){
              fscanf(file,"%f",&ivalRead[j][k]);
    	      //printf(" %d %1.0f  ",k, ivalRead[j][k]);
    	    }
	    if(j==0 || j==1){
	      lCalib->SetGain(j, 0);
	      if(j==0) lCalib->SetSector(j, 1);
	      else lCalib->SetSector(j, 4);
	    }
	    else if(j==2 || j==3){
	      lCalib->SetGain(j, 1);
	      if(j==2) lCalib->SetSector(j, 1);
	      else lCalib->SetSector(j, 4);
	    }
	    lCalib->SetfPMRefValue(j, ivalRead[j][0]);
	    lCalib->SetfPMRefWidth(j, ivalRead[j][1]);
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
else if(runType=="STANDALONE_EMD" && strcmp(beamType,"Pb-Pb")==0){
  TList* daqSources = GetFileSources(kDAQ, "EMDCALIB");
  if(!daqSources){
    AliError(Form("No sources for STANDALONE_EMD run %d !", fRun));
    return 1;
  }
  Log("\t List of sources for STANDALONE_EMD");
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
        Float_t fitValEMD[6]; Float_t equalCoeff[5][4];
        Float_t calibVal[4];
        for(Int_t j=0; j<10; j++){	   
          if(j<6){
            fscanf(file,"%f",&fitValEMD[j]);
            if(j<4){
              calibVal[j] = 2.76/fitValEMD[j];
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
