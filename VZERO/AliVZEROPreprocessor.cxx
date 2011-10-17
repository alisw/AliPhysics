#include "AliVZEROPreprocessor.h"
#include "AliVZEROCalibData.h"
#include "AliVZEROTriggerData.h"
#include "AliCDBMetaData.h"
#include "AliCDBEntry.h"
#include "AliDCSValue.h"
#include "AliLog.h"
#include "AliShuttleInterface.h"
#include "AliVZERODataFEE.h"
#include "AliVZERODataDCS.h"

#include <TFile.h>
#include <TTimeStamp.h>
#include <TObjString.h>
#include <TSystem.h>
#include <TH1F.h>


class Tlist;

//
//  This class is  a simple preprocessor for VZERO detector.
//
//  It gets High Voltage values for a given run from DCS and Pedestal values from DAQ 
//  and writes them as Calibration MetaData into OCDB/VZERO/Calib/Data
//  It also retrieves FEE parameters from DCS archive   
//  and writes them as Trigger MetaData into OCDB/VZERO/Trigger/Data 
//  (to be used for trigger simulation)
//

ClassImp(AliVZEROPreprocessor)

//______________________________________________________________________________________________
AliVZEROPreprocessor::AliVZEROPreprocessor(AliShuttleInterface* shuttle) :
	AliPreprocessor("V00", shuttle),
	fData(0),
	fFEEData(0)
 
{
  // constructor  
  
  AddRunType("STANDALONE_PULSER");
  AddRunType("STANDALONE_BC");
  AddRunType("PHYSICS");
    
}

//______________________________________________________________________________________________
AliVZEROPreprocessor::~AliVZEROPreprocessor()
{
  // destructor
	delete fFEEData;
	delete fData;
	
}

//______________________________________________________________________________________________
void AliVZEROPreprocessor::Initialize(Int_t run, UInt_t startTime,
	UInt_t endTime)
{
  // Creates AliVZERODataDCS object

   AliPreprocessor::Initialize(run, startTime, endTime);
  
   Log(Form("\n\tRun %d \n\tStartTime %s \n\tEndTime %s", run,
		TTimeStamp(startTime).AsString(),
		TTimeStamp(endTime).AsString()));

   fRun       = run;
   // fStartTime = startTime;
   // fEndTime   = endTime;
   fStartTime = GetStartTimeDCSQuery ();
   fEndTime   = GetEndTimeDCSQuery ();
   time_t daqStart = (time_t) (((TString)GetRunParameter("DAQ_time_start")).Atoi());
   time_t daqEnd   = (time_t) (((TString)GetRunParameter("DAQ_time_end")).Atoi());
   
	fData      = new AliVZERODataDCS(fRun, fStartTime, fEndTime,(UInt_t)daqStart, (UInt_t)daqEnd);
	fFEEData   = new AliVZERODataFEE(fRun, fStartTime, fEndTime);		
   
}

//______________________________________________________________________________________________
UInt_t AliVZEROPreprocessor::Process(TMap* dcsAliasMap)
{
  // Fills data retrieved from DCS and DAQ into a AliVZEROCalibData object and 
  // stores it into CalibrationDB


  // *** GET RUN TYPE ***
  TString runType = GetRunType();


  // *** REFERENCE DATA *** 
  
  TString fileName; 
  AliVZEROCalibData *calibData = new AliVZEROCalibData();
  
  // *************** HV From DCS ******************
  // Fills data into a AliVZERODataDCS object
  if(!dcsAliasMap) return 1;

 	// The Processing of the DCS input data is forwarded to AliVZERODataDCS
  if (!fData->ProcessData(*dcsAliasMap)) return 1;

	// Writes VZERO PMs HV values into VZERO calibration object and Timing resolution parameters
  	calibData->FillDCSData(fData);
	   
   // *************** From DAQ ******************
   
	TString sourcesId = "V00da_results";

	TList* sourceList = GetFileSources(kDAQ, sourcesId.Data());
  	if (!sourceList)  {
		Log(Form("No sources found for id %s", sourcesId.Data()));      		
      		return 1; }
	Log(Form("The following sources produced files with the id %s",sourcesId.Data()));
	sourceList->Print();    

  	TIter iter(sourceList);
  	TObjString *source;
		
	while((source=dynamic_cast<TObjString*> (iter.Next()))){
  		fileName = GetFile(kDAQ, sourcesId.Data(), source->GetName());
  		if (fileName.Length() > 0)
    		Log(Form("Got the file %s, now we can extract some values.", fileName.Data()));
		FILE *file;
		if((file = fopen(fileName.Data(),"r")) == NULL){
            	                   Log(Form("Cannot open file %s",fileName.Data()));
	    	  	           return 1;}
		Float_t pedMean[128], pedSigma[128], adcMean[128], adcSigma[128] ;
		for(Int_t j=0; j<128; j++)fscanf(file,"%f %f %f %f",
			                  &pedMean[j], &pedSigma[j], &adcMean[j], &adcSigma[j]);
		fclose(file);
	    	
		calibData->SetPedestal(pedMean);
		calibData->SetSigma(pedSigma);			
		calibData->SetADCsigma(adcSigma);
		}				

	delete source;      
  
  // Check that everything was properly transmitted

//   for(Int_t j=0; j<128; j++){printf("Pedestal[%d] -> %f \n",j,calibData->GetPedestal(j));}
//   for(Int_t j=0; j<128; j++){printf("pedSigma[%d] -> %f \n",j,calibData->GetSigma(j));}
//   for(Int_t j=0; j<128; j++){printf("Gain[%d] -> %f \n",j,calibData->GetGain(j));}
//   for(Int_t j=0; j<128; j++){printf("adcSigma[%d] -> %f \n",j,calibData->GetADCsigma(j));}
//   for(Int_t j=0; j<64; j++){printf("MeanHV[%d] -> %f \n",j,calibData->GetMeanHV(j));}
//   for(Int_t j=0; j<64; j++){printf("WidthHV[%d] -> %f \n",j,calibData->GetWidthHV(j));}
  
  // Now we store the VZERO Calibration Object into CalibrationDB

  Bool_t resECal=kTRUE;
  
  Bool_t result = 0;
//  if(sourceList && sourceList->GetEntries()>0)
//  {
  AliCDBMetaData metaData;
  metaData.SetBeamPeriod(0);
  metaData.SetResponsible("Brigitte Cheynis");
  metaData.SetComment("This preprocessor fills an AliVZEROCalibData object");

  resECal = Store("Calib", "Data", calibData, &metaData, 0, kTRUE);
//  }
  if(resECal==kFALSE ) result = 1;
  

  delete calibData;
  delete sourceList; 

 // -----------------------------------------------------------------------
 // Retrieves Front End Electronics Parameters from  DCS archive
 // -----------------------------------------------------------------------
	AliVZEROTriggerData *triggerData = new AliVZEROTriggerData();

 	// The processing of the DCS input data is forwarded to AliVZERODataFEE
	fFEEData->ProcessData(*dcsAliasMap);

	// Writes VZERO FEE parameters values into VZERO  Trigger parametrization object
	triggerData->FillData(fFEEData);

	// Stores the VZERO Trigger Object into TriggerDB
	
	resECal=kTRUE;
	
	metaData.SetBeamPeriod(0);
	metaData.SetResponsible("Brigitte Cheynis");
	metaData.SetComment("This preprocessor fills an AliVZEROTriggerData object");
	
	resECal = Store("Trigger", "Data", triggerData, &metaData, 0, kTRUE);
	if(resECal==kFALSE ) result = 1;
	

   // *************** From DAQ DA - Equalization factors ******************
  
  TH1F *eqFactors = new TH1F("VZEROEqualizationFactors","VZERO Equalization Factors for Pb-Pb",64,-0.5,63.5);
  sourcesId = "V00DAEqualFactors";

  TList* sourceList2 = GetFileSources(kDAQ, sourcesId.Data());
  if (!sourceList2)  {
    Log(Form("No sources found for id %s", sourcesId.Data()));      		
    return 1; }
  Log(Form("The following sources produced files with the id %s",sourcesId.Data()));
  sourceList2->Print();

  TIter iter2(sourceList2);
  TObjString *source2;
		
  while((source2=dynamic_cast<TObjString*> (iter2.Next()))){
    fileName = GetFile(kDAQ, sourcesId.Data(), source2->GetName());
    if (fileName.Length() > 0)
      Log(Form("Got the file %s, now we can extract some values.", fileName.Data()));
    FILE *file2;
    if((file2 = fopen(fileName.Data(),"r")) == NULL){
      Log(Form("Cannot open file %s",fileName.Data()));
      return 1;}

    Double_t alpha[66];
    alpha[0] = alpha[65] = 0;
    Int_t tempCh;
    Float_t tempAlpha;
    for(Int_t j=0; j<64; ++j) {
      fscanf(file2,"%d %f", &tempCh, &tempAlpha);
      alpha[tempCh+1] = tempAlpha;
    }
    fclose(file2);

    // Check that everything was properly transmitted
    printf("Equalization factors (0->64): ");
    for(Int_t j=0; j<64; ++j) printf("%.3f ",alpha[j+1]);
    printf("\n");

    eqFactors->SetContent(alpha);
  }

  delete source2;      
  
  // Now we store the VZERO Equalization Factors Object into OCDB

  resECal=kTRUE;
  
  AliCDBMetaData metaData2;
  metaData2.SetBeamPeriod(0);
  metaData2.SetResponsible("Brigitte Cheynis");
  metaData2.SetComment("VZERO Equalization Factors object filled by VZERO preprocessor");

  resECal = Store("Calib", "EqualizationFactors", eqFactors, &metaData2, 0, kTRUE);

  if(resECal==kFALSE ) result = 1;

  delete eqFactors;
  delete sourceList2; 

	
  return result;
}

