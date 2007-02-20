// --- ROOT system
#include <TFile.h>
#include <TTimeStamp.h>

#include "AliZDCPreprocessor.h"
#include "AliCDBManager.h"
#include "AliCDBEntry.h"
#include "AliCDBMetaData.h"
#include "AliDCSValue.h"
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

	AliInfo(Form("\n\tRun %d \n\tStartTime %s \n\tEndTime %s", run,
		TTimeStamp(startTime).AsString(),
		TTimeStamp(endTime).AsString()));

	fData = new AliZDCDataDCS(fRun, fStartTime, fEndTime);
}

//______________________________________________________________________________________________
UInt_t AliZDCPreprocessor::Process(TMap* dcsAliasMap)
{
  // Fills data into a AliZDCDataDCS object
  if(!dcsAliasMap) return 0;

  // The processing of the DCS input data is forwarded to AliZDCDataDCS
  Float_t DCSValues[26];
  fData->ProcessData(*dcsAliasMap, DCSValues);
  //dcsAliasMap->Print("");
  //
  AliZDCCalibData *calibdata = new AliZDCCalibData("ZDC");
  calibdata->SetDCSCalibData(DCSValues);

  const char* PedFileName = GetFile(kDAQ, "PEDESTALS", "LDC0");
  const Int_t NZDCch = 44;
  if(PedFileName){
    AliInfo(Form("File %s connected to analyze pedestal events", PedFileName));
    Float_t PedVal[(3*NZDCch)][2];
    for(Int_t i=0; i<(3*NZDCch); i++){
       for(Int_t j=0; j<2; j++){
          fscanf(file,"%f",&PedVal[i][j]);
	  printf("PedVal[%d][%d] -> %f \n",i,j,PedVal[i][j]);
       }
       if(i<NZDCch){
         calibdata->SetMeanPed(i,PedVal[i][0]);
         calibdata->SetMeanPedWidth(i,PedVal[i][1]);
       }
       else if(i>=NZDCch && i<(2*NZDCch)){
         calibdata->SetOOTPed(i,PedVal[i][0]);
         calibdata->SetOOTPedWidth(i,PedVal[i][1]);
       }
       else if(i>=(2*NZDCch) && i<(3*NZDCch)){
         calibdata->SetPedCorrCoeff(i,PedVal[i][0],PedVal[i][1]);
       }
    }
  }
  else AliInfo(Form("File %s not found", PedFileName));

  const char* EMDFileName = GetFile(kDAQ, "EMDCALIB", "LDC0");
  if(EMDFileName){
    AliInfo(Form("File %s connected to analyze EM dissociation events", EMDFileName));
    Float_t EMDFitVal[2];
    for(Int_t j=0; j<2; j++){          
      fscanf(file,"%f",&EMDFitVal[j]);
    }
    calibdata->SetEnCalib(EMDFitVal);
  }
  else AliInfo(Form("File %s not found", EMDFileName));
  //
  calibdata->Print("");
  
  // note that the parameters are returned as character strings!
  const char* nEvents = GetRunParameter("totalEvents");
  if (nEvents) {
  	Log(Form("Number of events for run %d: %s",fRun, nEvents));
  } else {
	Log(Form("Number of events not put in logbook!"));
  }

  //Now we have to store the final CDB file
  AliCDBMetaData metaData;
  metaData.SetBeamPeriod(0);
  metaData.SetResponsible("Chiara");
  metaData.SetComment("This preprocessor fills an AliZDCDataDCS object.");

  UInt_t result = Store("SHUTTLE","Data",fData, &metaData, 0, 0);
  delete fData;
  fData = 0;

  return result;
}

