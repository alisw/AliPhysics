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

/* $Id$ */

////////////////////////////////////////////////////////////////////////////
//                                                                        //
// This class is a first implementation for the TRD.                      //
// It takes data from HLT and computes the parameters                     //
// and stores both reference data and online calibration                  //
// parameters in the CDB                                                  //
// It alsotakes DCS data, does spline fits                                //
// and stores both reference data and spline fits results                 //
// in the CDB                                                             //
//                                                                        //
// Authors:                                                               //
//   R. Bailhache (R.Bailhache@gsi.de)                                    //
//   J. Book (jbook@ikf.uni-frankfurt.de)                                 //
//   W. Monange   (w.monange@gsi.de)                                      //
//   F. Kramer    (kramer@ikf.uni-frankfurt.de)                           //
// Maintainers:                                                           //
//   Hans.Beck@cern.ch                                                    //
//   R. Bailhache (R.Bailhache@gsi.de)                                    //
//                                                                        //
//                                                                        //
////////////////////////////////////////////////////////////////////////////

#include <fstream>

#include <TFile.h>
#include <TProfile2D.h>
#include <TObjString.h>
#include <TString.h>
#include <TList.h>
#include <TSAXParser.h>

#include "AliCDBMetaData.h"
#include "AliLog.h"

#include "AliTRDPreprocessor.h"
#include "AliTRDSensorArray.h"
#include "AliTRDCalibraFit.h"
#include "AliTRDCalibraMode.h"
#include "AliTRDCalibPadStatus.h"
#include "AliTRDSaxHandler.h"
#include "AliTRDgeometry.h"
#include "AliTRDCalibChamberStatus.h"
#include "AliTRDCalPad.h"
#include "AliTRDCalPadStatus.h"
#include "AliTRDCalDCSv2.h"
#include "AliTRDCalSingleChamberStatus.h"
#include "AliTRDCalChamberStatus.h"
#include "AliTRDCalROC.h"

ClassImp(AliTRDPreprocessor)

//______________________________________________________________________________________________
AliTRDPreprocessor::AliTRDPreprocessor(AliShuttleInterface *shuttle)
  :AliPreprocessor("TRD", shuttle)
  ,fCalDCSObjSOR(0)
  ,fCalDCSObjEOR(0)
  ,fVdriftHLT(0)
{
  //
  // Constructor
  //


  AddRunType("PHYSICS");
  AddRunType("STANDALONE");
  AddRunType("PEDESTAL");
  AddRunType("DAQ");
  
  
}
//______________________________________________________________________________________________
 AliTRDPreprocessor::AliTRDPreprocessor(const AliTRDPreprocessor&  ) :
   AliPreprocessor("TRD",0),
   fCalDCSObjSOR(0),
   fCalDCSObjEOR(0),
   fVdriftHLT(0)
{
  
  Fatal("AliTRDPreprocessor", "copy constructor not implemented");
  
}
//______________________________________________________________________________________________
AliTRDPreprocessor::~AliTRDPreprocessor()
{
  //
  // Destructor
  //

  if (fCalDCSObjSOR ) delete fCalDCSObjSOR;
  if (fCalDCSObjEOR ) delete fCalDCSObjEOR;

}
//______________________________________________________________________________________________
AliTRDPreprocessor& AliTRDPreprocessor::operator = (const AliTRDPreprocessor& )
{
  Fatal("operator =", "assignment operator not implemented");
  return *this;
}
//______________________________________________________________________________________________
void AliTRDPreprocessor::Initialize(Int_t run, UInt_t startTime, UInt_t endTime)
{
  //
  // Initialization routine for the TRD preprocessor
  //

  AliPreprocessor::Initialize(run,startTime,endTime);
  
}

//______________________________________________________________________________________________
UInt_t AliTRDPreprocessor::Process(TMap* dcsAliasMap)
{
  //
  // Process DCS and calibration part for HLT
  //

  const TString runType = GetRunType();
  Log(Form("runtype %s\n",runType.Data()));

  // We initialize the return values of all 
  // subfunctions to zero, ie 'good'
  Int_t DCSConfigReturn = 0;
  Int_t pedastalReturn  = 0;
  Int_t DCSReturn       = 0;
  Int_t HLTReturn       = 0;
  Int_t DAQReturn       = 0;
  
  // Always process the configuration data
  DCSConfigReturn = ProcessDCSConfigData();

  // Extract pedastals in pedastal runs
  if (runType=="PEDESTAL"){
    pedastalReturn = ExtractPedestals();
  }

  // Process the DCS sensors
  if ((runType=="PHYSICS") || (runType=="STANDALONE") || (runType=="DAQ")){
    // DCS
    DCSReturn = ProcessDCS(dcsAliasMap);

    /*
    if(runType=="PHYSICS"){
      // HLT if On
      //TString runPar = GetRunParameter("HLTStatus");
      //if(runPar=="1") {
      if(GetHLTStatus()) {
        HLTReturn = ExtractHLT();
      } 
      // DAQ if HLT failed
      if(!fVdriftHLT) {
        DAQReturn = ExtractDriftVelocityDAQ();
      }
    }
    */
   
  }

  // Check for any errors
  if(DCSConfigReturn || pedastalReturn || DCSReturn
     || HLTReturn || DAQReturn){
    // Error
    Log(Form("ERROR - Listing return values for subfunctions (0=good)"));
    Log(Form("DCSConfig %d, pedastal %d, DCS %d, HLT %d, DAQ %d"
	     ,DCSConfigReturn,pedastalReturn,DCSReturn
	     ,HLTReturn,DAQReturn));
    return 1;
  }
  else{
    // Good
    Log(Form("SUCCESS All subfuntions returned without error"));
    return 0;
  }
}
//______________________________________________________________________________
Bool_t AliTRDPreprocessor::ProcessDCS()
{
  //
  // Default process DCS method
  //

  TString runType = GetRunType();
  if ((runType == "PHYSICS") || (runType == "STANDALONE")) {
    return kTRUE;
  }
  return kFALSE;

}

//______________________________________________________________________________
Bool_t AliTRDPreprocessor::ProcessDCS(TMap *dcsAliasMap)
{
  //
  // Process DCS method
  //

  Bool_t error = kFALSE;

  AliCDBMetaData metaData;
  metaData.SetBeamPeriod(0);
  metaData.SetResponsible("Wilfried Monange/Raphaelle Bailhache");
  metaData.SetComment("TRD calib test");

  Log("****** DCS ******\n");
	
  TObjArray * list=AliTRDSensorArray::GetList();
	
  if (list == 0x0) {
    Log ("Error during AliTRDSensorArray::GetList");
    Log ("DCS will not be processing");
    return kTRUE;
  }

  const Int_t nEntries = list->GetEntries ();
  Log (Form ("%d alias loaded", nEntries));
		
  Bool_t results[nEntries];
  Int_t  nGraph[nEntries];
		
  for (Int_t iAlias = 0; iAlias < nEntries; iAlias++) {
			
    AliTRDSensorArray * oneTRDDCS = (AliTRDSensorArray *)list->At (iAlias);
			
    oneTRDDCS->SetStartTime (TTimeStamp (fStartTime));
    oneTRDDCS->SetEndTime (TTimeStamp (fEndTime));
			
    Log(Form("Processing DCS : \"%s\"", oneTRDDCS->GetStoreName ().Data ()));
			
    TMap *map=oneTRDDCS->ExtractDCS (dcsAliasMap);
		
    nGraph [iAlias] = map->GetEntries ();
		
    if (nGraph [iAlias] == 0) {
      Log("No TGraph for this dcsDatapointAlias : not stored");
      results [iAlias] = kFALSE;
      //error  = kTRUE;
      delete map;
      map=0;
      continue;
    }
		
    oneTRDDCS->SetGraph(map);

    // Remove duplicate entries where the value did not change
    // more than a relative precision. We take a std precision 
    // of 1e-4 and apply a higher precision for special cases
    Double_t precision=1e-4;
    // precision 1e-5 for HV
    if( oneTRDDCS->GetStoreName().Contains("trd_hvAnodeUmon") ||
	oneTRDDCS->GetStoreName().Contains("trd_hvDriftUmon") ){
      precision=1e-5;
    }
    oneTRDDCS->RemoveGraphDuplicates(precision);
    // Remove HV anode current points below a threshold (in uA)
    // when both neighbors are also below threshold
    if(oneTRDDCS->GetStoreName().Contains("trd_hvAnodeImon")){
      oneTRDDCS->RemoveAbsBelowThreshold(0.010);
    }
    if(oneTRDDCS->GetStoreName().Contains("trd_hvDriftImon")){
      oneTRDDCS->RemoveAbsBelowThreshold(0.10);
    }

    // Store the data point
    results[iAlias]=Store("Calib", oneTRDDCS->GetStoreName().Data(), oneTRDDCS, &metaData, 0, kFALSE); 
    delete map;
    map=0;
    oneTRDDCS->ClearGraph();
    //results [iAlias] = StoreReferenceData("Calib", oneTRDDCS->GetStoreName ().Data (), oneTRDDCS, &metaData); 

    if (!results[iAlias]) {
      AliError("Problem during StoreRef DCS");
      //error=kTRUE;
    }

    //BEGIN TEST (should not be removed ...)
    /*
    oneTRDDCS->ClearGraph();
    oneTRDDCS->ClearFit();
    oneTRDDCS->SetDiffCut2 (0.1);
    map=oneTRDDCS->ExtractDCS (dcsAliasMap);
    oneTRDDCS->SetGraph (map);
    Store("Calib", ("cut_"+oneTRDDCS->GetStoreName()).Data(), oneTRDDCS, &metaData, 0, kTRUE); 
    delete map;


    if(iAlias==1 || iAlias==19) continue;
    
    oneTRDDCS->ClearGraph();
    oneTRDDCS->ClearFit();
    oneTRDDCS->SetDiffCut2(0);
    map=oneTRDDCS->ExtractDCS(dcsAliasMap);
    oneTRDDCS->MakeSplineFit(map);
    Store("Calib", ("fit_"+oneTRDDCS->GetStoreName()).Data() , oneTRDDCS, &metaData, 0, kTRUE); 
    delete map;

     
    oneTRDDCS->ClearGraph(); 
    oneTRDDCS->ClearFit();
    oneTRDDCS->SetDiffCut2 (0.1);
    map=oneTRDDCS->ExtractDCS (dcsAliasMap);
    oneTRDDCS->MakeSplineFit(map);
    Store("Calib", ("cutfit_"+oneTRDDCS->GetStoreName()).Data() , oneTRDDCS, &metaData, 0, kTRUE); 
    delete map;
    */    
    //END TEST

  }

  // Check errors
  Int_t nbCount = 0;
  for (Int_t iAlias = 0; iAlias < nEntries; iAlias++) {
    if (results[iAlias]) {
      nbCount++;
    }
  }
  if(nbCount == 0) error = kTRUE;

		
  Log ("         Summury of DCS :\n");
  Log (Form("%30s %10s %10s", "dcsDatapointAlias", "Stored ?", "# graph"));
  for (Int_t iAlias = 0; iAlias < nEntries; iAlias++) {
    AliTRDSensorArray * oneTRDDCS = (AliTRDSensorArray *)list->At (iAlias);
    Log (Form ("%30s %10s %4d", 
	       oneTRDDCS->GetStoreName ().Data (),
	       results[iAlias] ? "ok" : "X",
	       nGraph [iAlias]));
  }
  Log ("*********** End of DCS **********");
  
  // Clean up
  delete list;
  list=0;
  
  return error;

}
//______________________________________________________________________________________________
Bool_t AliTRDPreprocessor::ExtractHalfChamberStatusDAQ()
{
  //
  // Half chamber status algorithm running on DAQ
  //

  Bool_t error = kFALSE;

  AliCDBMetaData metaData;
  metaData.SetBeamPeriod(0);
  metaData.SetResponsible("Raphaelle Bailhache");
  metaData.SetComment("TRD calib test");
  
  // Take the output of the DA on DAQ
  TList * listpad = GetFileSources(kDAQ,"HALFCHAMBERSTATUS");
  if (!listpad) {
    Log("No list found for the HalfChamberStatus");
    return kTRUE;
  }

  AliTRDCalibChamberStatus *calPed = 0x0;
    
  // loop through all files (only one normally)
  UInt_t index = 0;
  while (listpad->At(index)!=NULL) {
    TObjString* fileNameEntry = (TObjString*) listpad->At(index);
    if (fileNameEntry != NULL)
      {
        TString fileName = GetFile(kDAQ, "HALFCHAMBERSTATUS",
				   fileNameEntry->GetString().Data());
	if(fileName.Length() ==0){
	  Log(Form("Error by retrieving the file %d for the halfchamberstatus",(Int_t)index));
	  delete listpad;
	  return kTRUE;
	}
	
	TFile *f = TFile::Open(fileName);
	f->GetObject("calibchamberstatus",calPed);
	
	if(calPed) {
	  
	  // store as reference data
	  TString name("HalfChamberStatus");
	  if(!StoreReferenceData("DAQData",(const char *)name,(TObject *) calPed,&metaData)){
	    Log(Form("Error storing AliTRDCalibPadStatus object %d as reference data",(Int_t)index));
	    error = kTRUE;
	  }
	} // calPed
	else Log(Form("Error getting AliTRDCalibChamberStatus onject from file"));
		 
      } // fileNameEntry
    ++index;
  }// while (list)

  Log(Form("%d elements found in the list for the halfchamberstatus",(Int_t)index));
  if(index!=1){
    delete listpad;
    return kTRUE;
  }

  //
  // Produce the AliTRDCalChamberStatus  name calHalfChamberStatus
  //
  AliTRDCalChamberStatus *calHalfChamberStatus = 0x0;
  if(calPed) {
    //calPed->AnalyseHisto();   // check number of events, create calHalfChamberStatus (done on DAQ)
    if(fCalDCSObjEOR) {
      calPed->CheckEORStatus((AliTRDCalDCSv2 *)fCalDCSObjEOR);
    }
    calHalfChamberStatus=(AliTRDCalChamberStatus *)calPed->GetCalChamberStatus();
  }

  //
  // Store  
  //  

  AliCDBMetaData md3; 
  md3.SetObjectClassName("AliTRDCalChamberStatus");
  md3.SetResponsible("Raphaelle Bailhache");
  md3.SetBeamPeriod(1);
  md3.SetComment("TRD calib test");
  if(!Store("Calib","ChamberStatus"    ,(TObject *)calHalfChamberStatus, &md3, 0, kTRUE)){
    Log("Error storing the pedestal");
    delete listpad;
    return kTRUE;
  }
  
  delete listpad;
  return error; 
  
}
//______________________________________________________________________________________________
Bool_t AliTRDPreprocessor::ExtractPedestals()
{
  //
  // Pedestal running on LDCs at the DAQ
  //

  //
  // The reference data are stored in:
  // PadStatus1 for sm-00-01-02-09-10-11
  // PadStatus2 for sm-03-04-05-12-13-14
  // PadStatus3 for sm-06-07-08-15-16-17
  // PadStatus0 if nothing found..means problems
  //

  Bool_t error = kFALSE;

  // Init a AliTRDCalibPadStatus
  AliTRDCalibPadStatus calPedSum = AliTRDCalibPadStatus();

  AliCDBMetaData metaData;
  metaData.SetBeamPeriod(0);
  metaData.SetResponsible("Raphaelle Bailhache");
  metaData.SetComment("TRD calib test");
  
  // Sum the contributions of the LDCs
  TList * listpad = GetFileSources(kDAQ,"PADSTATUS");
  if (!listpad) {
    Log("No list found for the PEDESTRAL Run");
    return kTRUE;
  }
  
  // loop through all files from LDCs  
  UInt_t index = 0;
  while (listpad->At(index)!=NULL) {
    TObjString* fileNameEntry = (TObjString*) listpad->At(index);
    if (fileNameEntry != NULL)
      {
        TString fileName = GetFile(kDAQ, "PADSTATUS",
				   fileNameEntry->GetString().Data());
	if(fileName.Length() ==0){
	  Log(Form("Error by retrieving the file %d for the pedestal",(Int_t)index));
	  delete listpad;
	  return kTRUE;
	}
	
	TFile *f = TFile::Open(fileName);
	AliTRDCalibPadStatus *calPed;
	f->GetObject("calibpadstatus",calPed);
	
	if(calPed){
	  
	  // analyse
	  //calPed->AnalyseHisto();
	    	  
	  // Add to the calPedSum
	  for (Int_t idet=0; idet<540; idet++) {
	    AliTRDCalROC *rocMean  = calPed->GetCalRocMean(idet, kFALSE);
	    if ( rocMean )  {
	      calPedSum.SetCalRocMean(rocMean,idet);
	    }
	    AliTRDCalROC *rocRMS = calPed->GetCalRocRMS(idet, kFALSE);
	    if ( rocRMS )  {
	      calPedSum.SetCalRocRMS(rocRMS,idet);
	    }
	    AliTRDCalROC *rocMeand  = calPed->GetCalRocMeand(idet, kFALSE);
	    if ( rocMeand )  {
	      calPedSum.SetCalRocMeand(rocMeand,idet);
	    }
	    AliTRDCalROC *rocRMSd = calPed->GetCalRocRMSd(idet, kFALSE);
	    if ( rocRMSd )  {
	      calPedSum.SetCalRocRMSd(rocRMSd,idet);
	    }
	  }// det loop

	  // store as reference data
	  TString name("PadStatus");
	  name += index;
	  if(!StoreReferenceData("DAQData",(const char *)name,(TObject *) calPed,&metaData)){
	    Log(Form("Error storing AliTRDCalibPadStatus object %d as reference data",(Int_t)index));
	    error = kTRUE;
	  }

	} // calPed
      } // fileNameEntry
    ++index;
  }// while (list)

  Log(Form("%d elements found in the list for the pedestal",(Int_t)index));
  if(index==0){
    delete listpad;
    return kTRUE;
  }

  //
  // Create pedestal 
  //
    
  // Create Pad Status
  AliTRDCalPadStatus *calPadStatus = calPedSum.CreateCalPadStatus();
  // Create Noise 
  //Make the AliTRDCalPad
  AliTRDCalPad *calPad2 = calPedSum.CreateCalPad();
  //Make the AliTRDCalDet correspondant
  AliTRDCalDet *calDet = calPedSum.CreateCalDet();
 
  //
  // Take the noise and Pad status from the previous OCDB
  //

  AliTRDCalPad *calPadPrevious=0;
  AliCDBEntry* entry = GetFromOCDB("Calib", "PadNoise");
  if (entry) calPadPrevious = (AliTRDCalPad*)entry->GetObject();
  if ( calPadPrevious==NULL ) {
     Log("AliTRDPreprocsessor: No previous TRD pad noise entry available.\n");
     calPadPrevious = new AliTRDCalPad("PadNoise", "PadNoise");
  }

  AliTRDCalPadStatus *calPadStatusPrevious=0;
  entry = GetFromOCDB("Calib", "PadStatus");
  if (entry) calPadStatusPrevious = (AliTRDCalPadStatus*)entry->GetObject();
  if ( calPadStatusPrevious==NULL ) {
    Log("AliTRDPreprocsessor: No previous TRD pad status entry available.\n");
    calPadStatusPrevious = new AliTRDCalPadStatus("padstatus", "padstatus");
    for (Int_t idet=0; idet<540; ++idet)
      {
	AliTRDCalSingleChamberStatus *calROC = calPadStatusPrevious->GetCalROC(idet);
	for(Int_t k = 0; k < calROC->GetNchannels(); k++){
	  calROC->SetStatus(k,AliTRDCalPadStatus::kMasked);
	}
      }
  }

  
  // Loop over detectors for check
  for (Int_t det=0; det<AliTRDgeometry::kNdet; ++det)  {
    
    // noise
    AliTRDCalROC *calROCPreviousNoise = calPadPrevious->GetCalROC(det);
    AliTRDCalROC *calROCNoise         = calPad2->GetCalROC(det);

    // padstatus
    AliTRDCalSingleChamberStatus *calROCPreviousStatus = calPadStatusPrevious->GetCalROC(det);
    AliTRDCalSingleChamberStatus *calROCStatus         = calPadStatus->GetCalROC(det);
    
    
    // loop over first half and second half chamber
    for(Int_t half = 0; half < 2; half++){

      Bool_t data         = AreThereDataPedestal(calROCStatus,(Bool_t)half);
      //printf("There are data for the detector %d the half %d: %d\n",det,half,data);
      if(!data){
	// look if data in the OCDB
	Bool_t dataPrevious = AreThereDataPedestal(calROCPreviousStatus,(Bool_t)half);
	// if no data at all, set to default value
	if(!dataPrevious){
	  SetDefaultStatus(*calROCStatus,(Bool_t)half);
	  SetDefaultNoise(*calROCNoise,(Bool_t)half);
	}
	else{
	  // if data, set to previous value
	  SetStatus(*calROCStatus,calROCPreviousStatus,(Bool_t)half);
	  SetNoise(*calROCNoise,calROCPreviousNoise,(Bool_t)half);
	}
      }
    }
  }
  
  //
  // Store  
  //  

  AliCDBMetaData md3; 
  md3.SetObjectClassName("AliTRDCalPadStatus");
  md3.SetResponsible("Raphaelle Bailhache");
  md3.SetBeamPeriod(1);
  md3.SetComment("TRD calib test");
  if(!Store("Calib","PadStatus"    ,(TObject *)calPadStatus, &md3, 0, kTRUE)){
    Log("Error storing the pedestal");
    delete listpad;
    return kTRUE;
  }

  AliCDBMetaData md4; 
  md4.SetObjectClassName("AliTRDCalPad");
  md4.SetResponsible("Raphaelle Bailhache");
  md4.SetBeamPeriod(1);
  md4.SetComment("TRD calib test");
  if(!Store("Calib","PadNoise"    ,(TObject *)calPad2, &md4, 0, kTRUE)){
    Log("Error storing the pedestal");
    delete listpad;
    return kTRUE;
  }
  
  AliCDBMetaData md5; 
  md5.SetObjectClassName("AliTRDCalDet");
  md5.SetResponsible("Raphaelle Bailhache");
  md5.SetBeamPeriod(1);
  md5.SetComment("TRD calib test");
  if(!Store("Calib","DetNoise"    ,(TObject *)calDet, &md5, 0, kTRUE)){
    Log("Error storing the pedestal");
    delete listpad;
    return kTRUE;
  }  

  delete listpad;
  return error; 
  
}

//__________________________________________________________________
Bool_t AliTRDPreprocessor::AreThereDataPedestal(const AliTRDCalSingleChamberStatus * const calROCStatus
                                              , Bool_t second)
{

  //
  // Data for this half chamber
  //

  Bool_t data         = kFALSE;
  Int_t nCols         = calROCStatus->GetNcols();
  Int_t nCol0         = 0;
  Int_t nColE         = (Int_t) nCols/2 - 2;
  if(second) {
    nCol0 = nColE + 4;
    nColE = nCols;
  }

  Int_t totalnumberofpads = 0;
  Int_t totalnumberofdata = 0; 

  for(Int_t col = nCol0; col < nColE; col++){
    for(Int_t row = 0; row < calROCStatus->GetNrows(); row++){
      totalnumberofpads++;
      //printf("ismasked %d\n",(Int_t)calROCStatus->IsMasked(col,row));
      if(!calROCStatus->GetStatus(col,row)) {
	data = kTRUE;
	totalnumberofdata++;
      }
    }
  }
  if(totalnumberofdata < (Int_t)(totalnumberofpads/2)) data = kFALSE;

  return data;
  
}
//__________________________________________________________________
void AliTRDPreprocessor::SetDefaultStatus(AliTRDCalSingleChamberStatus &calROCStatus, Bool_t second){

  //
  // default status for this half chamber
  //

  Int_t nCols         = calROCStatus.GetNcols();
  Int_t nCol0         = 0;
  Int_t nColE         = (Int_t) nCols/2;
  if(second) {
    nCol0 = nColE;
    nColE = nCols;
  }
  for(Int_t col = nCol0; col < nColE; col++){
    for(Int_t row = 0; row < calROCStatus.GetNrows(); row++){
      calROCStatus.SetStatus(col,row,0);
    }
  }
}
//__________________________________________________________________
void AliTRDPreprocessor::SetStatus(AliTRDCalSingleChamberStatus &calROCStatus, AliTRDCalSingleChamberStatus *calROCStatusPrevious,Bool_t second){

  //
  // previous status for this half chamber
  //

  Int_t nCols         = calROCStatus.GetNcols();
  Int_t nCol0         = 0;
  Int_t nColE         = (Int_t) nCols/2;
  if(second) {
    nCol0 = nColE;
    nColE = nCols;
  }
  for(Int_t col = nCol0; col < nColE; col++){
    for(Int_t row = 0; row < calROCStatus.GetNrows(); row++){
      calROCStatus.SetStatus(col,row,calROCStatusPrevious->GetStatus(col,row));
    }
  }
}
//__________________________________________________________________
void AliTRDPreprocessor::SetDefaultNoise(AliTRDCalROC &calROCNoise, Bool_t second){

  //
  // default noise for this half chamber
  //

  Int_t nCols         = calROCNoise.GetNcols();
  Int_t nCol0         = 0;
  Int_t nColE         = (Int_t) nCols/2;
  if(second) {
    nCol0 = nColE;
    nColE = nCols;
  }
  for(Int_t col = nCol0; col < nColE; col++){
    for(Int_t row = 0; row < calROCNoise.GetNrows(); row++){
      calROCNoise.SetValue(col,row,0.12);
    }
  }
}
//__________________________________________________________________
void AliTRDPreprocessor::SetNoise(AliTRDCalROC &calROCNoise, AliTRDCalROC *calROCNoisePrevious, Bool_t second){

  //
  // previous noise for this half chamber
  //

  Int_t nCols         = calROCNoise.GetNcols();
  Int_t nCol0         = 0;
  Int_t nColE         = (Int_t) nCols/2;
  if(second) {
    nCol0 = nColE;
    nColE = nCols;
  }
  for(Int_t col = nCol0; col < nColE; col++){
    for(Int_t row = 0; row < calROCNoise.GetNrows(); row++){
      calROCNoise.SetValue(col,row,calROCNoisePrevious->GetValue(col,row));
    }
  }
}
//______________________________________________________________________________________________
Bool_t AliTRDPreprocessor::ExtractDriftVelocityDAQ()
{
  //
  // Drift velocity DA running on monitoring servers at the DAQ
  //

  Bool_t error = kFALSE; 

  // Objects for HLT and DAQ zusammen
  AliTRDCalibraFit *calibra = AliTRDCalibraFit::Instance();
  AliCDBMetaData metaData;
  metaData.SetBeamPeriod(0);
  metaData.SetResponsible("Raphaelle Bailhache");
  metaData.SetComment("TRD calib test");
  // Store the infos for the detector
  AliCDBMetaData md1; 
  md1.SetObjectClassName("AliTRDCalDet");
  md1.SetResponsible("Raphaelle Bailhache");
  md1.SetBeamPeriod(0);
  md1.SetComment("TRD calib test");
  // Store the infos for the pads
  AliCDBMetaData md2; 
  md2.SetObjectClassName("AliTRDCalPad");
  md2.SetResponsible("Raphaelle Bailhache");
  md2.SetBeamPeriod(0);
  md2.SetComment("TRD calib test");

  // Take the file from the DAQ file exchange server
  TList *listdaq = GetFileSources(kDAQ,"VDRIFT");
  if (!listdaq) {
    Log("No list found for vdrift (DAQ)");
    return kTRUE;
  }
  
  if(listdaq->GetSize() !=1){
    Log(Form("Problem on the size of the list: %d (DAQ)",listdaq->GetSize()));
    delete listdaq;
    return kTRUE;
  }
  
  TObjString* fileNameEntry = (TObjString*) listdaq->At(0);
  if(fileNameEntry != NULL){
    TString fileName = GetFile(kDAQ, "VDRIFT",
			       fileNameEntry->GetString().Data());
    if(fileName.Length() ==0){
      Log("Error retrieving the file vdrift (DAQ)");
      delete listdaq;
      return kTRUE;
    }
    TFile *filedaq = TFile::Open(fileName);
    TProfile2D *histodriftvelocity = (TProfile2D *) filedaq->Get("PH2d");
    if (histodriftvelocity) {
      
      // store as reference data
      if(!StoreReferenceData("DAQData","VdriftT0",(TObject *) histodriftvelocity,&metaData)){
	Log("Error storing 2D Profile for vdrift from the DAQ");
	error = kTRUE;
      }
      
      // analyse
      
      Log("Take the PH reference data. Now we will try to fit\n");
      calibra->SetMinEntries(2000); // If there is less than 2000
      calibra->AnalysePH(histodriftvelocity);
      
      Int_t nbtg = 6*4*18*((Int_t) ((AliTRDCalibraMode *)calibra->GetCalibraMode())->GetDetChamb0(1))
	+ 6*  18*((Int_t) ((AliTRDCalibraMode *)calibra->GetCalibraMode())->GetDetChamb2(1));
      Int_t nbfit        = calibra->GetNumberFit();
      Int_t nbfitSuccess = calibra->GetNumberFitSuccess();
      Int_t nbE        = calibra->GetNumberEnt();
      
      // if enough statistics store the results
      if ((nbtg >                  0) && 
	  (nbfit        >= 0.5*nbE) && (nbE > 30) && (nbfitSuccess > 30)) {
	// create the cal objects
	calibra->RemoveOutliers(1,kTRUE);
	calibra->PutMeanValueOtherVectorFit(1,kTRUE);
	calibra->RemoveOutliers2(kTRUE);
	calibra->PutMeanValueOtherVectorFit2(1,kTRUE);
	TObjArray object      = calibra->GetVectorFit();
	AliTRDCalDet *objdriftvelocitydet = calibra->CreateDetObjectVdrift(&object,kTRUE);
	//	TObject *objdriftvelocitypad = calibra->CreatePadObjectVdrift();
	object              = calibra->GetVectorFit2();
	AliTRDCalDet *objtime0det         = calibra->CreateDetObjectT0(&object,kTRUE);
	//	TObject *objtime0pad         = calibra->CreatePadObjectT0();
	calibra->ResetVectorFit();
	// store
	 if(!Store("Calib","ChamberVdrift"    ,(TObject *) objdriftvelocitydet,&md1,0,kTRUE)){
	 	  Log("Error storing the calibration object for the chamber vdrift (DAQ)");
	 	  error = kTRUE;
	 	}
	 	if(!Store("Calib","ChamberT0"        ,(TObject *) objtime0det        ,&md1,0,kTRUE)){
	 	  Log("Error storing the calibration object for the chamber t0 (DAQ)");
	 	  error = kTRUE;
	 	}
	// 	if(!Store("Calib","LocalVdrift"      ,(TObject *) objdriftvelocitypad,&md2,0,kTRUE)){
	// 	  Log("Error storing the calibration object for the local drift velocity (DAQ)");
	// 	  error = kTRUE;
	// 	}
	// 	if(!Store("Calib","LocalT0"          ,(TObject *) objtime0pad        ,&md2,0,kTRUE)){
	// 	  Log("Error storing the calibration object for the local time0 (DAQ)");
	// 	  error = kTRUE;
	// 	}
      }
      else{
	Log("Not enough statistics for the average pulse height (DAQ)");
      }
    } // histo here
  }// if fileNameEntry
  
  delete listdaq; 
  return error; 
  
}

//______________________________________________________________________________________________
Bool_t AliTRDPreprocessor::ExtractHLT()
{
  //
  // Gain, vdrift and PRF calibration running on HLT
  // return kTRUE if NULL pointer to the list
  //

  Bool_t error = kFALSE;

  // Objects for HLT and DAQ zusammen
  AliTRDCalibraFit *calibra = AliTRDCalibraFit::Instance();
  AliCDBMetaData metaData;
  metaData.SetBeamPeriod(0);
  metaData.SetResponsible("Raphaelle Bailhache");
  metaData.SetComment("TRD calib test");
  // Store the infos for the detector
  AliCDBMetaData md1; 
  md1.SetObjectClassName("AliTRDCalDet");
  md1.SetResponsible("Raphaelle Bailhache");
  md1.SetBeamPeriod(0);
  md1.SetComment("TRD calib test");
  // Store the infos for the pads
  AliCDBMetaData md2; 
  md2.SetObjectClassName("AliTRDCalPad");
  md2.SetResponsible("Raphaelle Bailhache");
  md2.SetBeamPeriod(0);
  md2.SetComment("TRD calib test");
  
  // Take the file from the HLT file exchange server
  TList *listhlt = GetFileSources(kHLT,"GAINDRIFTPRF");
  if (!listhlt) {
    Log("No list found for the HLT");
    return kTRUE;
  }

  if(listhlt->GetSize() != 1) {
    Log(Form("Problem on the size of the list: %d (HLT)",listhlt->GetSize()));
    delete listhlt;
    return kTRUE;
  }
  
  TObjString* fileNameEntry = (TObjString*) listhlt->At(0);
  if(fileNameEntry != NULL){
    TString fileName = GetFile(kHLT, "GAINDRIFTPRF",
			       fileNameEntry->GetString().Data());
    if(fileName.Length() ==0){
      Log("Error retrieving the file (HLT)");
      delete listhlt;
      return kTRUE;
    }
    // Take the file
    TFile *filehlt = TFile::Open(fileName);
    
    // gain
    TH2I *histogain = (TH2I *) filehlt->Get("CH2d");
    if (histogain) {
      histogain->SetDirectory(0);
      // store the reference data
      if(!StoreReferenceData("HLTData","Gain",(TObject *) histogain,&metaData)){
	Log("Error storing 2D histos for gain");
	error = kTRUE;
      }
      // analyse
      Log("Take the CH reference data. Now we will try to fit\n");
      calibra->SetMinEntries(800); // If there is less than 1000 entries in the histo: no fit
      calibra->AnalyseCH(histogain);
      Int_t nbtg = 6*4*18*((Int_t) ((AliTRDCalibraMode *)calibra->GetCalibraMode())->GetDetChamb0(0))
	+ 6*  18*((Int_t) ((AliTRDCalibraMode *)calibra->GetCalibraMode())->GetDetChamb2(0));
      Int_t nbfit       = calibra->GetNumberFit();
      Int_t nbE         = calibra->GetNumberEnt();
      // enough statistics
      if ((nbtg >                  0) && 
	  (nbfit        >= 0.5*nbE) && (nbE > 30)) {
	// create the cal objects
	calibra->PutMeanValueOtherVectorFit(1,kTRUE);
	TObjArray object           = calibra->GetVectorFit();
	AliTRDCalDet *objgaindet   = calibra->CreateDetObjectGain(&object);
	//	TObject *objgainpad        = calibra->CreatePadObjectGain();
	// store them
	if(!Store("Calib","ChamberGainFactor",(TObject *) objgaindet         ,&md1,0,kTRUE)){
	  Log("Error storing the calibration object for the chamber gain");
	  error = kTRUE;
	}
	// 	if(!Store("Calib","LocalGainFactor"  ,(TObject *) objgainpad         ,&md2,0,kTRUE)){
	// 	  Log("Error storing the calibration object for the local gain factor");
	// 	  error = kTRUE;
	// 	}
      }
      calibra->ResetVectorFit();
    }// if histogain
    
    // vdrift
    fVdriftHLT = kFALSE;
    TProfile2D *histodriftvelocity = (TProfile2D *) filehlt->Get("PH2d");
    if (histodriftvelocity) {
      histodriftvelocity->SetDirectory(0);
      // store the reference data
      if(!StoreReferenceData("HLTData","VdriftT0",(TObject *) histodriftvelocity,&metaData)){
	Log("Error storing 2D Profile for average pulse height (HLT)");
	error = kTRUE;
      }
      // analyse
      Log("Take the PH reference data. Now we will try to fit\n");
      calibra->SetMinEntries(800*20); // If there is less than 20000
      calibra->AnalysePH(histodriftvelocity);
      Int_t nbtg = 6*4*18*((Int_t) ((AliTRDCalibraMode *)calibra->GetCalibraMode())->GetDetChamb0(1))
	+ 6*  18*((Int_t) ((AliTRDCalibraMode *)calibra->GetCalibraMode())->GetDetChamb2(1));
      Int_t nbfit        = calibra->GetNumberFit();
      Int_t nbfitSuccess = calibra->GetNumberFitSuccess();
      Int_t nbE          = calibra->GetNumberEnt();
      // enough statistics
      if ((nbtg >                  0) && 
	  (nbfit        >= 0.5*nbE) && (nbE > 30) && (nbfitSuccess > 30)) {
	// create the cal objects
	calibra->RemoveOutliers(1,kTRUE);
	calibra->PutMeanValueOtherVectorFit(1,kTRUE);
	calibra->RemoveOutliers2(kTRUE);
	calibra->PutMeanValueOtherVectorFit2(1,kTRUE);
	TObjArray object  = calibra->GetVectorFit();
	AliTRDCalDet *objdriftvelocitydet = calibra->CreateDetObjectVdrift(&object,kTRUE);
	//	TObject *objdriftvelocitypad      = calibra->CreatePadObjectVdrift();
	object              = calibra->GetVectorFit2();
	AliTRDCalDet *objtime0det  = calibra->CreateDetObjectT0(&object,kTRUE);
	//TObject *objtime0pad       = calibra->CreatePadObjectT0();
	// store them
	if(!Store("Calib","ChamberVdrift"    ,(TObject *) objdriftvelocitydet,&md1,0,kTRUE)){
	  Log("Error storing the calibration object for the chamber vdrift (HLT)");
	  error = kTRUE;		
	}
	if(!Store("Calib","ChamberT0"        ,(TObject *) objtime0det        ,&md1,0,kTRUE)){
	  Log("Error storing the calibration object for the chamber t0 (HLT)");
	  error = kTRUE;
	}
	// if(!Store("Calib","LocalVdrift"      ,(TObject *) objdriftvelocitypad,&md2,0,kTRUE)){
	// 	  Log("Error storing the calibration object for the local drift velocity (HLT)");
	// 	  error = kTRUE;
	// 	}
	// 	if(!Store("Calib","LocalT0"          ,(TObject *) objtime0pad        ,&md2,0,kTRUE)){
	// 	  Log("Error storing the calibration object for the local time0 (HLT)");
	// 	  error = kTRUE;
	// 	}
	fVdriftHLT = kTRUE;
      }
      calibra->ResetVectorFit();
    }// if TProfile2D
    
    // prf
    /*
    TProfile2D *histoprf = (TProfile2D *) filehlt->Get("PRF2d");
    if (histoprf) {
      histoprf->SetDirectory(0);    
      // store reference data
      if(!StoreReferenceData("HLTData","PRF",(TObject *) histoprf,&metaData)){
	Log("Error storing the 2D Profile for Pad Response Function");
	error = kTRUE;
      }
      // analyse
      Log("Take the PRF reference data. Now we will try to fit\n");
      calibra->SetMinEntries(600); // If there is less than 20000
      calibra->AnalysePRFMarianFit(histoprf);
      Int_t nbtg = 6*4*18*((Int_t) ((AliTRDCalibraMode *)calibra->GetCalibraMode())->GetDetChamb0(2))
	+ 6*  18*((Int_t) ((AliTRDCalibraMode *)calibra->GetCalibraMode())->GetDetChamb2(2));
      Int_t nbfit        = calibra->GetNumberFit();
      Int_t nbE          = calibra->GetNumberEnt();
      // enough statistics
      if ((nbtg >                  0) && 
	  (nbfit        >= 0.95*nbE) && (nbE > 30)) {
	// create cal pad objects 
	TObjArray object            = calibra->GetVectorFit();
	TObject *objPRFpad          = calibra->CreatePadObjectPRF(&object);
	// store them
	if(!Store("Calib","PRFWidth"         ,(TObject *) objPRFpad          ,&md2,0,kTRUE)){
	  Log("Error storing the calibration object for the Pad Response Function");
	  error = kTRUE;
	}
      }
      calibra->ResetVectorFit();
    }// if PRF
    */
  }// if fileNameEntry
  
  delete listhlt;
  return error;
  
}

//_____________________________________________________________________________
UInt_t AliTRDPreprocessor::ProcessDCSConfigData()
{
  // process the configuration of FEE, PTR and GTU
  // reteive XML filei(s) from the DCS FXS
  // parse it/them and store TObjArrays in the CDB

  Log("Processing the DCS config summary files.");

  if(fCalDCSObjSOR) delete fCalDCSObjSOR; fCalDCSObjSOR=0;
  if(fCalDCSObjEOR) delete fCalDCSObjEOR; fCalDCSObjEOR=0;

  TString xmlFile[2];
  TString esor[2] = {"SOR", "EOR"};
  // get the XML files
  xmlFile[0] = GetFile(kDCS,"CONFIGSUMMARYSOR","");
  xmlFile[1] = GetFile(kDCS,"CONFIGSUMMARYEOR","");

  // check both files
  for (Int_t iFile=0; iFile<2; iFile++) {

    // check if the file are there
    if (xmlFile[iFile].IsNull()) {
      Log(Form("Warning: %s file not found!", esor[iFile].Data()));
      continue;
    }
    Log(Form("%s file found: %s", esor[iFile].Data(), xmlFile[iFile].Data()));

    // test the file
    std::ifstream fileTest;
    fileTest.open(xmlFile[iFile].Data(), std::ios_base::binary | std::ios_base::in);
    if (!fileTest.good() || fileTest.eof() || !fileTest.is_open()) {
      Log(Form("Warning: %s file not valid!", esor[iFile].Data()));
      continue;
    } 
    // check if the file is empty
    fileTest.seekg(0, std::ios_base::end);
    if (static_cast<int>(fileTest.tellg()) < 2) {
      Log(Form("Warning: %s file empty!", esor[iFile].Data()));
      continue;
    }
    Log(Form("%s file is valid.", esor[iFile].Data()));
    // Close the file descriptor
    fileTest.close();
    
    // make a robust XML validation
    TSAXParser testParser;
    if (testParser.ParseFile(xmlFile[iFile].Data()) < 0 ) {
      Log(Form("Warning: %s XML content is not well-formed!", esor[iFile].Data()));
      continue;
    }
    Log(Form("%s XML content is well-formed.", esor[iFile].Data()));

    // create parser and parse
    TSAXParser saxParser;
    AliTRDSaxHandler saxHandler;
    saxParser.ConnectToHandler("AliTRDSaxHandler", &saxHandler);
    saxParser.ParseFile(xmlFile[iFile].Data());

    // report errors of the parser if present
    if (saxParser.GetParseCode() != 0) {
      Log(Form("Warning: %s XML file validation failed! Parse code: %d", esor[iFile].Data(), saxParser.GetParseCode()));
      continue;
    }
    Log(Form("%s XML validation OK.", esor[iFile].Data()));

    // report errors of the handler if present
    if (saxHandler.GetHandlerStatus() != 0) {
      Log(Form("Warning: Creating %s calibration object failed! Error code: %d", esor[iFile].Data(), saxHandler.GetHandlerStatus()));
      continue;
    }
    Log(Form("%s SAX handler reports no errors.", esor[iFile].Data()));

    // get the calibration object storing the data from the handler
    AliTRDCalDCSv2 *calDCSObj = (AliTRDCalDCSv2*)saxHandler.GetCalDCSObj()->Clone();
    calDCSObj->EvaluateGlobalParameters();
    calDCSObj->SetRunType(GetRunType());
    calDCSObj->SetStartTime(GetStartTimeDCSQuery());
    calDCSObj->SetEndTime(GetEndTimeDCSQuery());
    if (iFile == 0) fCalDCSObjSOR = calDCSObj;
    if (iFile == 1) fCalDCSObjEOR = calDCSObj;
  }

  // If none of the two objects exists we don't even store a CDB entry
  if (!fCalDCSObjSOR && !fCalDCSObjEOR) {
    Log("ERROR: Failed reading both files!");
    return 1;
  }

  // put both objects in one TObjArray to store them
  TObjArray calObjArray(2);
  calObjArray.SetOwner();

  if (fCalDCSObjSOR) {
    calObjArray.AddAt(fCalDCSObjSOR,0);
    Log("TRDCalDCS object for SOR created.");
  }
  if (fCalDCSObjEOR) {
    calObjArray.AddAt(fCalDCSObjEOR,1);
    Log("TRDCalDCS object for EOR created.");
  }

  // store the DCS calib data in the CDB
  AliCDBMetaData metaData1;
  metaData1.SetBeamPeriod(0);
  metaData1.SetResponsible("Frederick Kramer");
  metaData1.SetComment("DCS configuration data in two AliTRDCalDCSv2 objects in one TObjArray (0:SOR, 1:EOR).");
  if (!Store("Calib", "DCS", &calObjArray, &metaData1, 0, kTRUE)) { Log("ERROR: Storing DCS config data object failed!"); return 1; }

  // ALICE policy is to always have an SOR object
  if(!fCalDCSObjSOR){
    Log("ERROR: Could not build an SOR object.");  
    return 1;
  }
  
  Log("SUCCESS: Processing of the DCS config summary file DONE.");  
  return 0;
}

