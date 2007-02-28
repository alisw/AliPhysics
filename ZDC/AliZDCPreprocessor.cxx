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
#include "AliAlignObjAngles.h"
#include "AliLog.h"
#include "AliZDCDataDCS.h"
#include "AliZDCCalibData.h"

//
// Class implementing ZDC pre-processor.
// It takes data from DCS and passes it to the class AliZDCDataDCS.
// The class is then written to the CDB.
//

ClassImp(AliZDCPreprocessor)

//______________________________________________________________________________________________
AliZDCPreprocessor::AliZDCPreprocessor(const char* detector, AliShuttleInterface* shuttle) :
  AliPreprocessor(detector, shuttle),
  fData(0)
{
  // constructor
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
UInt_t AliZDCPreprocessor::Process(TMap* dcsAliasMap)
{
  // *************** From DCS ******************
  // Fills data into a AliZDCDataDCS object
  if(!dcsAliasMap) return 0;

  // The processing of the DCS input data is forwarded to AliZDCDataDCS
  Float_t DCSValues[26];
  fData->ProcessData(*dcsAliasMap, DCSValues);
  //dcsAliasMap->Print("");
  //
  // --- Writing ZDC table positions into alignment object
  TClonesArray *array = new TClonesArray("AliAlignObjAngles",10);
  TClonesArray &alobj = *array;
  AliAlignObjAngles a;
  Double_t dx=0., dz=0., dpsi=0., dtheta=0., dphi=0.;
  Double_t dyZN1 = (Double_t) DCSValues[0];
  Double_t dyZP1 = (Double_t) DCSValues[1];
  Double_t dyZN2 = (Double_t) DCSValues[2];
  Double_t dyZP2 = (Double_t) DCSValues[3];
  const char *ZDCn1="ZDC/NeutronZDC1";
  const char *ZDCp1="ZDC/ProtonZDC1";
  const char *ZDCn2="ZDC/NeutronZDC2";
  const char *ZDCp2="ZDC/ProtonZDC2";
  UShort_t iIndex=0;
  AliAlignObj::ELayerID iLayer = AliAlignObj::kInvalidLayer;
  UShort_t volid = AliAlignObj::LayerToVolUID(iLayer,iIndex);
  //
  new(alobj[0]) AliAlignObjAngles(ZDCn1, volid, dx, dyZN1, dz, dpsi, dtheta, dphi, kTRUE);
  new(alobj[1]) AliAlignObjAngles(ZDCp1, volid, dx, dyZP1, dz, dpsi, dtheta, dphi, kTRUE);
  new(alobj[2]) AliAlignObjAngles(ZDCn2, volid, dx, dyZN2, dz, dpsi, dtheta, dphi, kTRUE);
  new(alobj[3]) AliAlignObjAngles(ZDCp2, volid, dx, dyZP2, dz, dpsi, dtheta, dphi, kTRUE);
  // save in CDB storage
  AliCDBMetaData md;
  md.SetResponsible("Chiara Oppedisano");
  md.SetComment("Alignment object for ZDC");
  Store("Align","Data", array, &md, 0, 0);
  
  // --- Writing ZDC PTMs HV values into calibration object
  AliZDCCalibData *calibdata = new AliZDCCalibData("ZDC");
  for(Int_t j=0; j<22; j++)
     calibdata->SetPMTHVVal(j,DCSValues[j+4]);

  
  // *************** From DAQ ******************
  // [a] PEDESTALS
  TList* Pedfilesources = GetFileSources(kDAQ, "PEDESTALS");
  if(!Pedfilesources){
    AliError(Form("No PEDESTALS file source for run %d !", fRun));
    return 0;
  }
  Log("\t List of sources for PEDESTALS");
  Pedfilesources->Print();
  //
  TIter iter(Pedfilesources);
  TObjString* source;
  Int_t i=0;
  UInt_t result = 0;
  while((source=dynamic_cast<TObjString*> (iter.Next()))){
       Log(Form("\n\t Getting file #%d\n",++i));
       TString stringp = GetFile(kDAQ, "PEDESTALS", source->GetName());
       const char * PedFileName = stringp.Data();
       const Int_t NZDCch = 44;
       if(PedFileName){
         FILE *file;
         if((file = fopen(PedFileName,"r")) == NULL){
           printf("Cannot open file %s \n",PedFileName);
	   return 0;
         }
         Log(Form("File %s connected to analyze pedestal events", PedFileName));
         Float_t PedVal[(3*NZDCch)][2];
         for(Int_t i=0; i<(3*NZDCch); i++){
            for(Int_t j=0; j<2; j++){
               fscanf(file,"%f",&PedVal[i][j]);
	       //if(j==1) printf("PedVal[%d] -> %f, %f \n",i,PedVal[i][0],PedVal[i][1]);
            }
            if(i<NZDCch){
              calibdata->SetMeanPed(i,PedVal[i][0]);
              calibdata->SetMeanPedWidth(i,PedVal[i][1]);
            }
            else if(i>=NZDCch && i<(2*NZDCch)){
              calibdata->SetOOTPed(i-NZDCch,PedVal[i][0]);
              calibdata->SetOOTPedWidth(i-NZDCch,PedVal[i][1]);
            }
            else if(i>=(2*NZDCch) && i<(3*NZDCch)){
              calibdata->SetPedCorrCoeff(i-(2*NZDCch),PedVal[i][0],PedVal[i][1]);
            }
         }
       }
       else{
          Log(Form("File %s not found", PedFileName));
          return 0;
       }
       //
      //calibdata->Print("");
  }

  // [a] EMD EVENTS
  TList* EMDfilesources = GetFileSources(kDAQ, "EMDCALIB");
  if(!EMDfilesources){
    AliError(Form("No EMDCALIB file source for run %d !", fRun));
    return 0;
  }
  Log("\t List of sources for EMDCALIB");
  EMDfilesources->Print();
  //
  TIter iter2(EMDfilesources);
  TObjString* source2;
  Int_t j=0;
  while((source2=dynamic_cast<TObjString*> (iter2.Next()))){
       Log(Form("\n\t Getting file #%d\n",++j));
       TString stringe = GetFile(kDAQ, "EMDCALIB", source2->GetName());
       const char* EMDFileName = stringe.Data();
       if(EMDFileName){
    	 FILE *file;
    	 if((file = fopen(EMDFileName,"r")) == NULL){
    	   printf("Cannot open file %s \n",EMDFileName);
	   return 0;
    	 }
    	 Log(Form("File %s connected to analyze EM dissociation events", EMDFileName));
    	 Float_t EMDFitVal[2];
    	 for(Int_t j=0; j<2; j++){	    
    	   fscanf(file,"%f",&EMDFitVal[j]);
         }
         calibdata->SetEnCalib(EMDFitVal);
       }
       else{
         Log(Form("File %s not found", EMDFileName));
         return 0;
       }
       //calibdata->Print("");
  }
  
  // note that the parameters are returned as character strings!
  const char* nEvents = GetRunParameter("totalEvents");
  if(nEvents) Log(Form("Number of events for run %d: %s",fRun, nEvents));
  else Log(Form("Number of events not put in logbook!"));

  // Storing the final CDB file
  AliCDBMetaData metaData;
  metaData.SetBeamPeriod(0);
  metaData.SetResponsible("Chiara");
  metaData.SetComment("Filling AliZDCCalibData object");

  result = Store("Calib","Data",calibdata, &metaData, 0, 0);
 
  return result;
}

