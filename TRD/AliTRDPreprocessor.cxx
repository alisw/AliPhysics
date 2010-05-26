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
//   W. Monange   (w.monange@gsi.de)                                      //
//   F. Kramer    (kramer@ikf.uni-frankfurt.de)                           //
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
#include "Cal/AliTRDCalPad.h"
#include "Cal/AliTRDCalPadStatus.h"
#include "Cal/AliTRDCalDCS.h"
#include "Cal/AliTRDCalSingleChamberStatus.h"
#include "Cal/AliTRDCalROC.h"

ClassImp(AliTRDPreprocessor)

//______________________________________________________________________________________________
AliTRDPreprocessor::AliTRDPreprocessor(AliShuttleInterface *shuttle)
  :AliPreprocessor("TRD", shuttle)
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
AliTRDPreprocessor::~AliTRDPreprocessor()
{
  //
  // Destructor
  //

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

  TString runType = GetRunType();
  Log(Form("runtype %s\n",runType.Data()));
  
  // always process the configuration data
  Int_t DCSConfigReturn = ProcessDCSConfigData();
  if(DCSConfigReturn) return DCSConfigReturn; 
  
  if (runType=="PEDESTAL"){
    if(ExtractPedestals()) return 1;
    return 0;
  } 

  if ((runType=="PHYSICS") || (runType=="STANDALONE") || (runType=="DAQ")){
    // DCS
    if(ProcessDCS(dcsAliasMap)) return 1; 
    if(runType=="PHYSICS"){
      // HLT if On
      //TString runPar = GetRunParameter("HLTStatus");
      //if(runPar=="1") {
      if(GetHLTStatus()) {
	//if(ExtractHLT()) return 1; // for testing!
	ExtractHLT();
      } 
      // DAQ if HLT failed
      if(!fVdriftHLT) {
	//if(ExtractDriftVelocityDAQ()) return 1; 
	ExtractDriftVelocityDAQ(); // for testing!
      }
    }
  }
  
  return 0;  
  
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
	
  TObjArray * list=AliTRDSensorArray::GetList ();
	
  if (list == 0x0) {
    Log ("Error during AliTRDSensorArray::GetList");
    Log ("DCS will not be processing");
    return kTRUE;
  }

  Int_t nEntries = list->GetEntries ();
  Log (Form ("%d alias loaded", nEntries));
		
  Bool_t * results=new Bool_t [nEntries];
  Int_t  * nGraph=new Int_t [nEntries];
		
  for (Int_t iAlias = 0; iAlias < nEntries; iAlias++) {
			
    AliTRDSensorArray * oneTRDDCS = (AliTRDSensorArray *)list->At (iAlias);
			
    oneTRDDCS->SetStartTime (TTimeStamp (fStartTime));
    oneTRDDCS->SetEndTime (TTimeStamp (fEndTime));
			
    Log(Form("Processing DCS : \"%s\"", oneTRDDCS->GetStoreName ().Data ()));
			
    TMap * map;

    map=oneTRDDCS->ExtractDCS (dcsAliasMap);
		
    nGraph [iAlias] = map->GetEntries ();
		
    if (nGraph [iAlias] == 0) {
      Log("No TGraph for this dcsDatapointAlias : not stored");
      results [iAlias] = kFALSE;
      //error  = kTRUE;
      continue;
    }
		
    oneTRDDCS->SetGraph(map);
    results[iAlias]=Store("Calib", oneTRDDCS->GetStoreName().Data(), oneTRDDCS, &metaData, 0, kFALSE); 
    delete map;		

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
  
  delete results;
  delete nGraph;

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
	  
	  Int_t ldc = 0; 

	  // analyse
	  //calPed->AnalyseHisto();
	    	  
	  // Add to the calPedSum
	  for (Int_t idet=0; idet<540; idet++) {
	    AliTRDCalROC *rocMean  = calPed->GetCalRocMean(idet, kFALSE);
	    if ( rocMean )  {
	      calPedSum.SetCalRocMean(rocMean,idet);
	      ldc = (Int_t) (idet / 30);
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

	  if((ldc==0) || (ldc==1) || (ldc==2)) ldc = 1;
	  if((ldc==3) || (ldc==4) || (ldc==5)) ldc = 2;
	  if((ldc==6) || (ldc==7) || (ldc==8)) ldc = 3;
	  if((ldc==9) || (ldc==10) || (ldc==11)) ldc = 4;
	  if((ldc==12) || (ldc==13) || (ldc==14)) ldc = 5;
	  if((ldc==15) || (ldc==16) || (ldc==17)) ldc = 6;
	
	  // store as reference data
	  TString name("PadStatus");
	  name += ldc;
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
Bool_t AliTRDPreprocessor::AreThereDataPedestal(AliTRDCalSingleChamberStatus * const calROCStatus
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
	//TObject *objdriftvelocitypad = calibra->CreatePadObjectVdrift();
	object              = calibra->GetVectorFit2();
	AliTRDCalDet *objtime0det         = calibra->CreateDetObjectT0(&object,kTRUE);
	//TObject *objtime0pad         = calibra->CreatePadObjectT0();
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
	//if(!Store("Calib","LocalVdrift"      ,(TObject *) objdriftvelocitypad,&md2,0,kTRUE)){
	//  Log("Error storing the calibration object for the local drift velocity (DAQ)");
	//  error = kTRUE;
	//}
	//if(!Store("Calib","LocalT0"          ,(TObject *) objtime0pad        ,&md2,0,kTRUE)){
	//  Log("Error storing the calibration object for the local time0 (DAQ)");
	//  error = kTRUE;
	//}
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
	//TObject *objgainpad        = calibra->CreatePadObjectGain();
	// store them
	if(!Store("Calib","ChamberGainFactor",(TObject *) objgaindet         ,&md1,0,kTRUE)){
	  Log("Error storing the calibration object for the chamber gain");
	  error = kTRUE;
	}
	//if(!Store("Calib","LocalGainFactor"  ,(TObject *) objgainpad         ,&md2,0,kTRUE)){
	//  Log("Error storing the calibration object for the local gain factor");
	//  error = kTRUE;
	//}
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
	//TObject *objdriftvelocitypad      = calibra->CreatePadObjectVdrift();
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
	//if(!Store("Calib","LocalVdrift"      ,(TObject *) objdriftvelocitypad,&md2,0,kTRUE)){
	//  Log("Error storing the calibration object for the local drift velocity (HLT)");
	//  error = kTRUE;
	//}
	//if(!Store("Calib","LocalT0"          ,(TObject *) objtime0pad        ,&md2,0,kTRUE)){
	//  Log("Error storing the calibration object for the local time0 (HLT)");
	//  error = kTRUE;
	//}
	fVdriftHLT = kTRUE;
      }
      calibra->ResetVectorFit();
    }// if TProfile2D
    
    // prf
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
  }// if fileNameEntry
  
  delete listhlt;
  return error;
  
}

//_____________________________________________________________________________
UInt_t AliTRDPreprocessor::ProcessDCSConfigData()
{
  // 
  // process the configuration of FEE, PTR and GTU
  // reteive XML filei(s) from the DCS FXS
  // parse it/them and store TObjArrays in the CDB
  //
  // return 0 for success, otherwise:
  //  5 : Could not get the (SOR and EOR) or SOR file from the FXS
  //  8 : Both files are not valid
  // 10 : SOR XML is not well-formed
  // 11 : EOR XML is not well-formed
  // 12 : ERROR in SOR XML SAX validation: something wrong with the content
  // 12 : ERROR in EOR XML SAX validation: something wrong with the content
  // 14 : ERROR while creating SOR calibration objects in the handler
  // 15 : ERROR while creating EOR calibration objects in the handler
  // 16 : ERROR while storing data in the CDB
  //

  Log("Processing the DCS config summary files.");

  // get the XML files
  Log("Requesting the 2 summaryfiles from the FXS..");
  TString xmlFileS = GetFile(kDCS,"CONFIGSUMMARYSOR","");
  TString xmlFileE = GetFile(kDCS,"CONFIGSUMMARYEOR","");
  // Request EOR and SOR files from the fxs, if both are not found exit

  Bool_t fileExistE = kTRUE, fileExistS = kTRUE;

  // check if the files are there
  if (xmlFileS.IsNull()) {
	  Log(Form("Warning: SOR file %s not found!", xmlFileS.Data()));
    fileExistS = kFALSE;
  } else Log(Form("SOR file found: %s", xmlFileS.Data()));

  if (xmlFileE.IsNull()) {
    Log(Form("Warning: EOR file %s not found!", xmlFileE.Data()));
    fileExistE = kFALSE;
  } else Log(Form("EOR file found: %s", xmlFileE.Data()));

  if (fileExistS==0 && fileExistE==0) {
    Log(Form("ERROR: SOR and EOR files not found: %s and %s", xmlFileS.Data(), xmlFileE.Data()));
    return 5;
  }

  // test the files
  if (fileExistS) {
    Log("Checking if SOR file is valid.");
    std::ifstream fileTestS;
    fileTestS.open(xmlFileS.Data(), std::ios_base::binary | std::ios_base::in);
    if (!fileTestS.good() || fileTestS.eof() || !fileTestS.is_open()) {
      Log(Form("Warning: File %s not valid!",xmlFileS.Data()));
      fileExistS = kFALSE;
      return 5;
    }
    Log("Checking if SOR file is not empty.");
    fileTestS.seekg(0, std::ios_base::end);
    if (static_cast<int>(fileTestS.tellg()) < 2) {
      Log(Form("Warning: File %s is empty!",xmlFileS.Data()));
      fileExistS = kFALSE;
      return 5;
    }
  }


  if (fileExistE) {
    Log("Checking if EOR file is valid.");
    std::ifstream fileTestE;
    fileTestE.open(xmlFileE.Data(), std::ios_base::binary | std::ios_base::in);
    if (!fileTestE.good() || fileTestE.eof() || !fileTestE.is_open()) {
      Log(Form("Warning: File %s not valid!",xmlFileE.Data()));
      fileExistE = kFALSE;
    }
    Log("Checking if EOR file is not empty.");
    fileTestE.seekg(0, std::ios_base::end);
    if (static_cast<int>(fileTestE.tellg()) < 2) {
      Log(Form("Warning: File %s is empty!",xmlFileE.Data()));
      fileExistE = kFALSE;
    }
  }

  if (fileExistS==0 && fileExistE==0) {
    Log("ERROR: Both files (SOR/EOR) are not valid!");
    return 8;
  }

  Log("One or both of the tested files are valid.");

  // make a robust XML validation
  TSAXParser testParser;
  if (fileExistS && testParser.ParseFile(xmlFileS.Data()) < 0 ) {
    Log("ERROR: XML content (SOR) is not well-formed.");
    return 10;
  } else if (fileExistE && testParser.ParseFile(xmlFileE.Data()) < 0 ) {
    Log("ERROR: XML content (EOR) is not well-formed.");
    return 11;
  }
  Log("XML contents are well-formed.");

  // create parser and parse
  TSAXParser saxParserS, saxParserE;
  AliTRDSaxHandler saxHandlerS, saxHandlerE;
  if (fileExistS) {
    saxParserS.ConnectToHandler("AliTRDSaxHandler", &saxHandlerS);
    saxParserS.ParseFile(xmlFileS.Data());
  }
  if (fileExistE) {
    saxParserE.ConnectToHandler("AliTRDSaxHandler", &saxHandlerE);
    saxParserE.ParseFile(xmlFileE.Data());
  }
  // report errors of the parser if present
  if (fileExistS && saxParserS.GetParseCode() != 0) {
    Log(Form("ERROR in XML file validation. SOR Parse Code: %s", saxParserS.GetParseCode()));
    return 12;
  }
  if (fileExistE && saxParserE.GetParseCode() != 0) {
    Log(Form("ERROR in XML file validation. EOR Parse Code: %s", saxParserE.GetParseCode()));
    return 13;
  }
  Log("XML file validation OK.");
  // report errors of the handler if present
  if (fileExistS && saxHandlerS.GetHandlerStatus() != 0) {
    Log(Form("ERROR while creating calibration objects. SOR Error code: %s", saxHandlerS.GetHandlerStatus()));
    return 14;  
  }
  if (fileExistE && saxHandlerE.GetHandlerStatus() != 0) {
    Log(Form("ERROR while creating calibration objects. EOR Error code: %s", saxHandlerE.GetHandlerStatus()));
    return 15;
  }
  Log("SAX handler reports no errors.");

  // put both objects in one TObjArray to store them
  TObjArray* fCalObjArray = new TObjArray(2);
  fCalObjArray->SetOwner();

  // get the calibration object storing the data from the handler
  if (fileExistS) {
    AliTRDCalDCS* fCalDCSObjSOR = saxHandlerS.GetCalDCSObj();
    fCalDCSObjSOR->EvaluateGlobalParameters();
    fCalDCSObjSOR->SetRunType(GetRunType());
    fCalDCSObjSOR->SetStartTime(GetStartTimeDCSQuery());
    fCalDCSObjSOR->SetEndTime(GetEndTimeDCSQuery());
    fCalObjArray->AddAt(fCalDCSObjSOR,0);
    Log("TRDCalDCS object for SOR created.");
  }

  if (fileExistE) {
    AliTRDCalDCS* fCalDCSObjEOR = saxHandlerE.GetCalDCSObj();
    fCalDCSObjEOR->EvaluateGlobalParameters();
    fCalDCSObjEOR->SetRunType(GetRunType());
    fCalDCSObjEOR->SetStartTime(GetStartTimeDCSQuery());
    fCalDCSObjEOR->SetEndTime(GetEndTimeDCSQuery());
    fCalObjArray->AddAt(fCalDCSObjEOR,1);
    Log("TRDCalDCS object for EOR created.");
  }

  // store the DCS calib data in the CDB
  AliCDBMetaData metaData1;
  metaData1.SetBeamPeriod(0);
  metaData1.SetResponsible("Frederick Kramer");
  metaData1.SetComment("DCS configuration data in two AliTRDCalDCS objects in one TObjArray (0:SOR, 1:EOR).");
  if (!Store("Calib", "DCS", fCalObjArray, &metaData1, 0, kTRUE)) {
    Log("problems while storing DCS config data object");
    return 16;
  } else {
    Log("DCS config data object stored.");
  }

  //delete fCalObjArray;

  Log("SUCCESS: Processing of the DCS config summary file DONE.");  
  return 0;
}



