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
// Author:                                                                //
//   R. Bailhache (R.Bailhache@gsi.de)                                    //
//   W. Monange   (w.monange@gsi.de)                                      //
//                                                                        //
////////////////////////////////////////////////////////////////////////////

#include "AliTRDPreprocessor.h"

#include <TFile.h>
#include <TProfile2D.h>
#include <TStopwatch.h>
#include <TObjString.h>
#include <TString.h>
#include <TList.h>
#include <TCollection.h>

#include "AliCDBMetaData.h"
#include "AliLog.h"

#include "AliTRDSensorArray.h"
#include "AliTRDCalibraFit.h"
#include "AliTRDCalibraMode.h"
#include "AliTRDCalibPadStatus.h"
#include "Cal/AliTRDCalDet.h"
#include "Cal/AliTRDCalPadStatus.h"

ClassImp(AliTRDPreprocessor)

//______________________________________________________________________________________________
AliTRDPreprocessor::AliTRDPreprocessor(AliShuttleInterface *shuttle)
  :AliPreprocessor("TRD", shuttle),
  fResult(0),
  fVdriftHLT(0)
{
  //
  // Constructor
  //

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

  fResult = 0;

  //
  // DCS
  //
  
  ProcessDCS(dcsAliasMap);

  //
  // Process the calibration data for the HLT and DAQ part
  //

  TString runType = GetRunType();
  printf("runtype %s\n",(const char*)runType);


  if (strcmp(runType, "PEDESTAL") == 0){
    ExtractPedestals();
  }
  
  if (strcmp(runType, "PHYSICS") == 0){
    printf("test\n");
    //if(GetRunParameter("HLTStatus")==1) ExtractHLT();
    ExtractHLT();
    ExtractDriftVelocityDAQ();
  }
  printf("result is %d\n",fResult);  
  return fResult;  
  
}


//______________________________________________________________________________
void AliTRDPreprocessor::ProcessDCS(TMap * dcsAliasMap)
{

  AliCDBMetaData metaData;
  metaData.SetBeamPeriod(0);
  metaData.SetResponsible("Wilfried Monange/Raphaelle Bailhache");
  metaData.SetComment("TRD calib test");
	
	
  Log ("****** DCS ******\n");
	
  TObjArray * list = AliTRDSensorArray::GetList ();
	
  if (list == 0x0) {
    Log ("Error during AliTRDSensorArray::GetList");
    Log ("DCS will not be processing");
    return;
  }

  Int_t nEntries = list->GetEntries ();
  Log (Form ("%d alias loaded", nEntries));
		
  Bool_t * results=new Bool_t [nEntries];
  Int_t  * nGraph=new Int_t [nEntries];
		
  for (Int_t iAlias = 0; iAlias < nEntries; iAlias++) {
			
    AliTRDSensorArray * oneTRDDCS = (AliTRDSensorArray *)list->At (iAlias);
			
    oneTRDDCS->SetStartTime (TTimeStamp (fStartTime));
    oneTRDDCS->SetEndTime (TTimeStamp (fEndTime));
			
    Log (Form("Processing DCS : \"%s\"", oneTRDDCS->GetStoreName ().Data ()));
			
    TMap * map;

    map=oneTRDDCS->ExtractDCS (dcsAliasMap);
		
    nGraph [iAlias] = map->GetEntries ();
		
    if (nGraph [iAlias] == 0) {
      Log("No TGraph for this dcsDatapointAlias : not stored");
      results [iAlias] = kFALSE;
      continue;
    }
		
    oneTRDDCS->SetGraph(map);
    results[iAlias]=Store("Calib", oneTRDDCS->GetStoreName().Data(), oneTRDDCS, &metaData, 0, kTRUE); 
    delete map;		

    //results [iAlias] = StoreReferenceData("Calib", oneTRDDCS->GetStoreName ().Data (), oneTRDDCS, &metaData); 


    if (!results[iAlias]) {
      AliError("Problem during StoreRef DCS"); 
      fResult|=kEStore;
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


}

//______________________________________________________________________________________________
void AliTRDPreprocessor::ExtractPedestals()
{
  //
  // Pedestal running on LDCs at the DAQ
  //

  // Init a AliTRDCalibPadStatus
  AliTRDCalibPadStatus calPedSum = AliTRDCalibPadStatus();

  AliCDBMetaData metaData;
  metaData.SetBeamPeriod(0);
  metaData.SetResponsible("Raphaelle Bailhache");
  metaData.SetComment("TRD calib test");
  
  // Sum the contributions of the LDCs
  TList * listpad = GetFileSources(kDAQ,"PADSTATUS");
  if (listpad) {
    
    // loop through all files from LDCs

    UInt_t index = 0;
    while (listpad->At(index)!=NULL) {
      TObjString* fileNameEntry = (TObjString*) listpad->At(index);
      if (fileNameEntry != NULL) {
        TString fileName = GetFile(kDAQ, "PADSTATUS",
				   fileNameEntry->GetString().Data());
	if(fileName.Length() !=0){
	  TFile *f = TFile::Open(fileName);
	  AliTRDCalibPadStatus *calPed;
	  f->GetObject("calibpadstatus",calPed);
	  
	  if(calPed){

	    // analyse
	    calPed->AnalyseHisto();
	    
	    // store as reference data
	    TString name("PadStatus");
	    name += (Int_t)index;
	    if(!StoreReferenceData("DAQData",(const char *)name,(TObject *) calPed,&metaData)){
	      Log(Form("Error storing AliTRDCalibPadStatus object %d as reference data",(Int_t)index));
	      fResult |= kEStore;
	    }
	    
	       
	    // Add to the calPedSum
	    for (Int_t idet=0; idet<540; idet++) {
	      AliTRDCalROC *rocMean  = calPed->GetCalRocMean(idet, kFALSE);
	      if ( rocMean )  calPedSum.SetCalRocMean(rocMean,idet);
	      AliTRDCalROC *rocRMS = calPed->GetCalRocRMS(idet, kFALSE);
	      if ( rocRMS )  calPedSum.SetCalRocRMS(rocRMS,idet);
	    }// det loop
	    
	  } // calPed
	}
	else{
	  Log(Form("Error by retrieving the file %d for the pedestal",(Int_t)index));
	  fResult |= kEGetFileDAQ;  
	}
      }// fileNameEntry length
      ++index;
    }// while (list)
    Log(Form("%d elements found in the list for the pedestal",(Int_t)index));
    if(index==0){
      fResult |= kEEmptyListDAQ;
    }
    //
    // Store pedestal entry to OCDB
    //


    // Create Pad Status
    AliTRDCalPadStatus *calPadStatus = calPedSum.CreateCalPadStatus();
    AliCDBMetaData md3; 
    md3.SetObjectClassName("AliTRDCalPadStatus");
    md3.SetResponsible("Raphaelle Bailhache");
    md3.SetBeamPeriod(1);
    md3.SetComment("TRD calib test");
    if(!Store("Calib","PadStatus"    ,(TObject *)calPadStatus, &md3, 0, kTRUE)){
      Log("Error storing the pedestal");
      fResult |= kEStore;    
    }
    delete listpad;
  }
  else {
    Log("No list found for the PEDESTRAL Run");
    fResult |= kEGetFileDAQ;    
  }
  

}
//______________________________________________________________________________________________
void AliTRDPreprocessor::ExtractDriftVelocityDAQ()
{
  //
  // Drift velocity DA running on monitoring servers at the DAQ
  //


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
  if (listdaq) {
    if(listdaq->GetSize() ==1){
      TObjString* fileNameEntry = (TObjString*) listdaq->At(0);
      if(fileNameEntry != NULL){
	TString fileName = GetFile(kDAQ, "VDRIFT",
				   fileNameEntry->GetString().Data());
	if(fileName.Length() !=0){
	  TFile *filedaq = TFile::Open(fileName);
	  TProfile2D *histodriftvelocity = (TProfile2D *) filedaq->Get("PH2d");
	  if (histodriftvelocity) {
	    
	    // store as reference data
	    if(!StoreReferenceData("DAQData","VdriftT0",(TObject *) histodriftvelocity,&metaData)){
	      Log("Error storing 2D Profile for vdrift from the DAQ");
	      fResult |= kEStore;
	    }
	    
	    // analyse
	    if(!fVdriftHLT){
	      Log("Take the PH reference data. Now we will try to fit\n");
	      calibra->SetMinEntries(2000); // If there is less than 2000
	      calibra->AnalysePH(histodriftvelocity);
	      
	      Int_t nbtg = 6*4*18*((Int_t) ((AliTRDCalibraMode *)calibra->GetCalibraMode())->GetDetChamb0(1))
		+ 6*  18*((Int_t) ((AliTRDCalibraMode *)calibra->GetCalibraMode())->GetDetChamb2(1));
	      Int_t nbfit        = calibra->GetNumberFit();
	      Int_t nbE        = calibra->GetNumberEnt();
	      
	      // if enough statistics store the results
	      if ((nbtg >                  0) && 
		  (nbfit        >= 0.95*nbE)) {
		// create the cal objects
		TObjArray object      = calibra->GetVectorFit();
		AliTRDCalDet *objdriftvelocitydet = calibra->CreateDetObjectVdrift(&object,kTRUE);
		TObject *objdriftvelocitypad = calibra->CreatePadObjectVdrift();
		object              = calibra->GetVectorFit2();
		AliTRDCalDet *objtime0det         = calibra->CreateDetObjectT0(&object,kTRUE);
		TObject *objtime0pad         = calibra->CreatePadObjectT0();
		calibra->ResetVectorFit();
		// store
		if(!Store("Calib","ChamberVdrift"    ,(TObject *) objdriftvelocitydet,&md1,0,kTRUE)){
		  Log("Error storing the calibration object for the chamber vdrift (DAQ)");
		  fResult |= kEStore;    
		}
		if(!Store("Calib","ChamberT0"        ,(TObject *) objtime0det        ,&md1,0,kTRUE)){
		  Log("Error storing the calibration object for the chamber t0 (DAQ)");
		  fResult |= kEStore;    
		}
		if(!Store("Calib","LocalVdrift"      ,(TObject *) objdriftvelocitypad,&md2,0,kTRUE)){
		  Log("Error storing the calibration object for the local drift velocity (DAQ)");
		  fResult |= kEStore;
		}
		if(!Store("Calib","LocalT0"          ,(TObject *) objtime0pad        ,&md2,0,kTRUE)){
		  Log("Error storing the calibration object for the local time0 (DAQ)");
		  fResult |= kEStore;
		}
	      }
	      else{
		Log("Not enough statistics for the average pulse height (DAQ)");
	      }
	    } // analyse
	  } // histo here
	}
	else{
	  Log("Error retrieving the file vdrift (DAQ)");
	  if(!fVdriftHLT) fResult |= kEGetFileDAQ;
	}
      }// if fileNameEntry
    } // size DAQ list
    else{
      Log(Form("Problem on the size of the list: %d (DAQ)",listdaq->GetSize()));
      if(!fVdriftHLT) fResult |= kEEmptyListDAQ;
    }
    delete listdaq; 
  } 
  else{
    Log("No list found for vdrift (DAQ)");
    if(!fVdriftHLT){
      fResult |= kEGetFileDAQ;
    }
  }
   
}
//______________________________________________________________________________________________
void AliTRDPreprocessor::ExtractHLT()
{
  //
  // Gain, vdrift and PRF calibration running on HLT
  // return kFALSE if NULL pointer to the list
  //

 
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
  if (listhlt) {
    if (listhlt->GetSize() == 1) {
      TObjString* fileNameEntry = (TObjString*) listhlt->At(0);
      if(fileNameEntry != NULL){
	TString fileName = GetFile(kHLT, "GAINDRIFTPRF",
				   fileNameEntry->GetString().Data());
	if(fileName.Length() !=0){
	  // Take the file
	  TFile *filehlt = TFile::Open(fileName);
	  
	  
	  // gain
	  TH2I *histogain = (TH2I *) filehlt->Get("CH2d");
	  histogain->SetDirectory(0);
	  if (histogain) {
	    // store the reference data
	    if(!StoreReferenceData("HLTData","Gain",(TObject *) histogain,&metaData)){
	      Log("Error storing 2D histos for gain");
	      fResult |= kEStore;
	    }
	    // analyse
	    Log("Take the CH reference data. Now we will try to fit\n");
	    calibra->SetMinEntries(1000); // If there is less than 1000 entries in the histo: no fit
	    calibra->AnalyseCH(histogain);
	    Int_t nbtg = 6*4*18*((Int_t) ((AliTRDCalibraMode *)calibra->GetCalibraMode())->GetDetChamb0(0))
	      + 6*  18*((Int_t) ((AliTRDCalibraMode *)calibra->GetCalibraMode())->GetDetChamb2(0));
	    Int_t nbfit       = calibra->GetNumberFit();
	    Int_t nbE         = calibra->GetNumberEnt();
	    // enough statistics
	    if ((nbtg >                  0) && 
		(nbfit        >= 0.95*nbE)) {
	      // create the cal objects
	      TObjArray object           = calibra->GetVectorFit();
	      AliTRDCalDet *objgaindet   = calibra->CreateDetObjectGain(&object,calibra->GetScaleFitFactor(),kTRUE);
	      TObject *objgainpad        = calibra->CreatePadObjectGain();
	      // store them
	      if(!Store("Calib","ChamberGainFactor",(TObject *) objgaindet         ,&md1,0,kTRUE)){
		Log("Error storing the calibration object for the chamber gain");
		fResult |= kEStore;
	      }
	      if(!Store("Calib","LocalGainFactor"  ,(TObject *) objgainpad         ,&md2,0,kTRUE)){
		Log("Error storing the calibration object for the local gain factor");
		fResult |= kEStore;
	      }
	    }
	    calibra->ResetVectorFit();
	  }// if histogain
	  
	  
	  
	  // vdrift
	  fVdriftHLT = kFALSE;
	  TProfile2D *histodriftvelocity = (TProfile2D *) filehlt->Get("PH2d");
	  histodriftvelocity->SetDirectory(0);
	  if (histodriftvelocity) {
	    // store the reference data
	    if(!StoreReferenceData("HLTData","VdriftT0",(TObject *) histodriftvelocity,&metaData)){
	      Log("Error storing 2D Profile for average pulse height (HLT)");
	      fResult |= kEStore;
	    }
	    // analyse
	    Log("Take the PH reference data. Now we will try to fit\n");
	    calibra->SetMinEntries(1000*20); // If there is less than 20000
	    calibra->AnalysePH(histodriftvelocity);
	    Int_t nbtg = 6*4*18*((Int_t) ((AliTRDCalibraMode *)calibra->GetCalibraMode())->GetDetChamb0(1))
	      + 6*  18*((Int_t) ((AliTRDCalibraMode *)calibra->GetCalibraMode())->GetDetChamb2(1));
	    Int_t nbfit        = calibra->GetNumberFit();
	    Int_t nbE          = calibra->GetNumberEnt();
	    // enough statistics
	    if ((nbtg >                  0) && 
		(nbfit        >= 0.95*nbE)) {
	      // create the cal objects
	      TObjArray object  = calibra->GetVectorFit();
	      AliTRDCalDet *objdriftvelocitydet = calibra->CreateDetObjectVdrift(&object,kTRUE);
	      TObject *objdriftvelocitypad      = calibra->CreatePadObjectVdrift();
	      object              = calibra->GetVectorFit2();
	      AliTRDCalDet *objtime0det  = calibra->CreateDetObjectT0(&object,kTRUE);
	      TObject *objtime0pad       = calibra->CreatePadObjectT0();
	      // store them
	      if(!Store("Calib","ChamberVdrift"    ,(TObject *) objdriftvelocitydet,&md1,0,kTRUE)){
		Log("Error storing the calibration object for the chamber vdrift (HLT)");
		fResult |= kEStore;    
	      }
	      if(!Store("Calib","ChamberT0"        ,(TObject *) objtime0det        ,&md1,0,kTRUE)){
		Log("Error storing the calibration object for the chamber t0 (HLT)");
		fResult |= kEStore;    
	      }
	      if(!Store("Calib","LocalVdrift"      ,(TObject *) objdriftvelocitypad,&md2,0,kTRUE)){
		Log("Error storing the calibration object for the local drift velocity (HLT)");
		fResult |= kEStore;
	      }
	      if(!Store("Calib","LocalT0"          ,(TObject *) objtime0pad        ,&md2,0,kTRUE)){
		Log("Error storing the calibration object for the local time0 (HLT)");
		fResult |= kEStore;
	      }
	      fVdriftHLT = kTRUE;
	    }
	    calibra->ResetVectorFit();
	  }// if TProfile2D
	  
	  
	  // prf
	  TProfile2D *histoprf = (TProfile2D *) filehlt->Get("PRF2d");
	  histoprf->SetDirectory(0);
	  if (histoprf) {
	    // store reference data
	    if(!StoreReferenceData("HLTData","PRF",(TObject *) histoprf,&metaData)){
	      Log("Error storing the 2D Profile for Pad Response Function");
	      fResult |= kEStore;
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
		(nbfit        >= 0.95*nbE)) {
	      // create cal pad objects 
	      TObjArray object            = calibra->GetVectorFit();
	      TObject *objPRFpad          = calibra->CreatePadObjectPRF(&object);
	      // store them
	      if(!Store("Calib","PRFWidth"         ,(TObject *) objPRFpad          ,&md2,0,kTRUE)){
		Log("Error storing the calibration object for the Pad Response Function");
		fResult |= kEStore;
	      }
	    }
	    calibra->ResetVectorFit();
	  }// if PRF
	}
	else{
	  Log("Error retrieving the file (HLT)");
	  fResult |= kEGetFileHLT;
	} 
      } // fileNameEntry
    } // if HLT size of tlist
    else{
      Log(Form("Problem on the size of the list: %d (HLT)",listhlt->GetSize()));
      fResult |= kEEmptyListHLT;
    }
    delete listhlt;
  } //if HLT tlist
  else{
    Log("No list found for the HLT");
    fResult |= kEGetFileHLT;
  }
  
}
