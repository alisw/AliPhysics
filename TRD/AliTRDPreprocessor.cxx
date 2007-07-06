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
                   :AliPreprocessor("TRD", shuttle)
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

  UInt_t result = 0;

  //
  // DCS
  //
  
  AliCDBMetaData metaData;
  metaData.SetBeamPeriod(0);
  metaData.SetResponsible("Wilfried Monange/Raphaelle Bailhache");
  metaData.SetComment("TRD calib test");
	
	
  Log ("****** DCS ******\n");
	
  TObjArray * list = AliTRDSensorArray::GetList ();
	
  if (list == 0x0) {
	  Log ("Error during AliTRDSensorArray::GetList");
	  Log ("DCS will not be processing");
  }else { 
	
	Int_t nEntries = list->GetEntries ();
	Log (Form ("%d alias loaded", nEntries));
		
	Bool_t * results = new Bool_t [nEntries];
	Int_t  * nGraph = new Int_t [nEntries];
		
	for (Int_t iAlias = 0; iAlias < nEntries; iAlias++) {
			
		AliTRDSensorArray * oneTRDDCS = (AliTRDSensorArray *)list->At (iAlias);
			
		oneTRDDCS->SetStartTime (TTimeStamp (fStartTime));
		oneTRDDCS->SetEndTime (TTimeStamp (fEndTime));
			
		Log (Form("Processing DCS : \"%s\"", 
			oneTRDDCS->GetStoreName ().Data ()));
			
		TMap * map = oneTRDDCS->ExtractDCS (dcsAliasMap);
		
		nGraph [iAlias] = map->GetEntries ();
		
		if (nGraph [iAlias] == 0) {
			Log("No TGraph for this dcsDatapointAlias : not stored");
			results [iAlias] = kFALSE;
			result |= kEStoreRefDCS;
			continue;
		}
		
		oneTRDDCS->SetGraph (map);

		results [iAlias] = Store("Calib", 
								 oneTRDDCS->GetStoreName ().Data (), 
								 oneTRDDCS,  
								 &metaData,
								 0, 
								 kTRUE); 
		
		/*	
		results [iAlias] = StoreReferenceData("Calib", 
												oneTRDDCS->GetStoreName ().Data (), 
												oneTRDDCS, 
												&metaData); 
		*/
		if (!results [iAlias]) {
			AliError("Problem during StoreRef DCS"); 
			result |= kEStoreRefDCS;
		}
			
		delete map;
			
		/*
			//BEGIN TEST
		oneTRDDCS->SetDiffCut2 (0.1);
		map = oneTRDDCS->ExtractDCS (dcsAliasMap);
		oneTRDDCS->SetGraph (map);
			
		StoreReferenceData("Calib", 
							(oneTRDDCS->GetStoreName ()+"Cut").Data(), 
							oneTRDDCS, &metaData); 
		delete map;
			//END TEST
		*/
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

  //
  // Process the calibration data for the HLT part
  //

 // Objects for HLT and DAQ zusammen
  AliTRDCalibraFit *calibra = AliTRDCalibraFit::Instance();
  Bool_t hltVdrifthisto = kFALSE;
  // Store the infos for the detector
  AliCDBMetaData *md1= new AliCDBMetaData(); 
  md1->SetObjectClassName("AliTRDCalDet");
  md1->SetResponsible("Raphaelle Bailhache");
  md1->SetBeamPeriod(0);
  md1->SetComment("TRD calib test");
  // Store the infos for the pads
  AliCDBMetaData *md2= new AliCDBMetaData(); 
  md2->SetObjectClassName("AliTRDCalPad");
  md2->SetResponsible("Raphaelle Bailhache");
  md2->SetBeamPeriod(0);
  md2->SetComment("TRD calib test");
 
  //
  // Process the calibration data for the HLT part
  //

  // How long does it take for the HLT part?
  TStopwatch timer;
  timer.Start();

  //Run type
  TString runType = GetRunType();
  Log(Form("Run type for run %d: %s", fRun, runType.Data()));
  if (strcmp(runType, "PHYSICS") != 0){
    Log("Nothing to do for HLT!");
  }
  else{
    // note that the parameters are returned as character strings!
    const char* nEvents = GetRunParameter("totalEvents");
    if (nEvents) {
      Log(Form("Number of events for run %d: %s",fRun, nEvents));
    } else {
      Log(Form("Number of events not put in logbook!"));
    }

    // Take the file from the HLT file exchange server
    TList *filesourceshlt = GetFileSources(kHLT,"GAINDRIFTPRF");
    if (!filesourceshlt) {
      Log(Form("No sources found for GAINDRIFTPRF for run %d !",fRun));
      result |= kEListFileHLT;
    }
    else{
      if (filesourceshlt->GetSize() != 1) {
	Log(Form("More than one source found for GAINDRIFTPRF for run %d!",fRun));
	filesourceshlt->Print();
	result |= kEListFileHLT;
      }
      else {

	//Debug mode
	//calibra->SetDebugLevel(2);       //Debug

	// Loop over the files taken from the HLT
	TIter iter(filesourceshlt);
	TObjString *sourcehlt;
	while ((sourcehlt = dynamic_cast<TObjString *> (iter.Next()))) {
    
	  TString filenamehlt = GetFile(kHLT,"GAINDRIFTPRF",sourcehlt->GetName());
	  if (filenamehlt.Length() == 0) {
	    Log(Form("Error retrieving file from source %d failed!", sourcehlt->GetName()));
	    result |= kEOpenFileHLT;
	  }
	  else{

	    // Init some things
	    TH2I         *histogain           = 0x0; // Histogram taken from HLT for gain factor
	    TProfile2D   *histodriftvelocity  = 0x0; // Profile taken from HLT for drift velocity and time0
	    TProfile2D   *histoprf            = 0x0; // Profile taken from HLT for prf


	    // Take the histos
	    TFile *filehlt = TFile::Open(filenamehlt);
	    histogain = (TH2I *) filehlt->Get("CH2d");
	    histogain->SetDirectory(0);
	    if (!histogain) {
	      Log("Error retrieving 2D histos for gain failed!");
	      result |= kETakeHistoHLT;
	    }
	    histodriftvelocity = (TProfile2D *) filehlt->Get("PH2d");
	    histodriftvelocity->SetDirectory(0);
	    if (!histodriftvelocity) {
	      Log("Error retrieving 2D Profile for average pulse height failed!");
	      result |= kETakeHistoHLT;
	    }
	    histoprf = (TProfile2D *) filehlt->Get("PRF2d");
	    histoprf->SetDirectory(0);
	    if (!histoprf) {
	      Log("Error retrieving 2D Profile for Pad Response Function failed!");
	      result |= kETakeHistoHLT;
	    }
	    filehlt->Close();
	    
	    // try to fit them and store
	    if (histogain) {
	      if(!StoreReferenceData("HLTData","Gain",(TObject *) histogain,&metaData)){
		Log("Error storing 2D histos for gain as reference data");
		result |= kEStoreHistoHLT;
	      }
	      Log("Take the CH reference data. Now we will try to fit\n");
	      calibra->SetMinEntries(100); // If there is less than 100 entries in the histo: no fit
	      calibra->AnalyseCH(histogain);
	      Int_t nbtg = 6*4*18*((Int_t) ((AliTRDCalibraMode *)calibra->GetCalibraMode())->GetDetChamb0(0))
		+ 6*  18*((Int_t) ((AliTRDCalibraMode *)calibra->GetCalibraMode())->GetDetChamb2(0));
	      Int_t nbfit       = calibra->GetNumberFit();
	      Int_t nbE         = calibra->GetNumberEnt();
	      if ((nbtg >                  0) && 
		  (nbfit        >= 0.95*nbE)) {
		TObjArray object             = calibra->GetVectorFit();
		AliTRDCalDet *objgaindet   = calibra->CreateDetObjectGain(&object,calibra->GetScaleFitFactor(),kTRUE);
		TObject *objgainpad        = calibra->CreatePadObjectGain();
		if(!Store("Calib","ChamberGainFactor",(TObject *) objgaindet         ,md1,0,kTRUE)){
		  Log("Error storing the calibration object for the chamber gain");
		  result |= kEStoreCalHLT;
		}
		if(!Store("Calib","LocalGainFactor"  ,(TObject *) objgainpad         ,md2,0,kTRUE)){
		  Log("Error storing the calibration object for the local gain factor");
		  result |= kEStoreCalHLT;
		}
	      }
	      else{
		Log("Not enough statistics for the gain");
		result |= kEFitHistoHLT;
	      }
	      calibra->ResetVectorFit();
	    }
    
	    if (histodriftvelocity) {
	      if(!StoreReferenceData("HLTData","VdriftT0",(TObject *) histodriftvelocity,&metaData)){
		Log("Error storing 2D Profile for average pulse height as reference data");
		result |= kEStoreHistoHLT;
	      }
	      Log("Take the PH reference data. Now we will try to fit\n");
	      calibra->SetMinEntries(100*20); // If there is less than 2000
	      calibra->AnalysePH(histodriftvelocity);
	      Int_t nbtg = 6*4*18*((Int_t) ((AliTRDCalibraMode *)calibra->GetCalibraMode())->GetDetChamb0(1))
		+ 6*  18*((Int_t) ((AliTRDCalibraMode *)calibra->GetCalibraMode())->GetDetChamb2(1));
	      Int_t nbfit        = calibra->GetNumberFit();
	      Int_t nbE          = calibra->GetNumberEnt();
	      if ((nbtg >                  0) && 
		  (nbfit        >= 0.95*nbE)) {
		TObjArray object  = calibra->GetVectorFit();
		AliTRDCalDet *objdriftvelocitydet = calibra->CreateDetObjectVdrift(&object,kTRUE);
		TObject *objdriftvelocitypad      = calibra->CreatePadObjectVdrift();
		object              = calibra->GetVectorFit2();
		AliTRDCalDet *objtime0det  = calibra->CreateDetObjectT0(&object,kTRUE);
		TObject *objtime0pad       = calibra->CreatePadObjectT0();
		if(!Store("Calib","ChamberVdrift"    ,(TObject *) objdriftvelocitydet,md1,0,kTRUE)){
		  Log("Error storing the calibration object for the chamber vdrift");
		  result |= kEStoreCalHLT;    
		}
		if(!Store("Calib","ChamberT0"        ,(TObject *) objtime0det        ,md1,0,kTRUE)){
		  Log("Error storing the calibration object for the chamber t0");
		  result |= kEStoreCalHLT;    
		}
		if(!Store("Calib","LocalVdrift"      ,(TObject *) objdriftvelocitypad,md2,0,kTRUE)){
		  Log("Error storing the calibration object for the local drift velocity");
		  result |= kEStoreCalHLT;
		}
		if(!Store("Calib","LocalT0"          ,(TObject *) objtime0pad        ,md2,0,kTRUE)){
		  Log("Error storing the calibration object for the local time0");
		  result |= kEStoreCalHLT;
		}
		hltVdrifthisto = kTRUE;
	      }
	      else{
		Log("Not enough statistics for the average pulse height");
		result |= kEFitHistoHLT;
	      }	     
	      calibra->ResetVectorFit();
	    }
	    
	    if (histoprf) {
	      if(!StoreReferenceData("HLTData","PRF",(TObject *) histoprf,&metaData)){
		Log("Error storing the 2D Profile for Pad Response Function as reference data");
		result |= kEStoreHistoHLT;
	      }
	      Log("Take the PRF reference data. Now we will try to fit\n");
	      calibra->SetMinEntries(100*20); // If there is less than 2000
	      calibra->AnalysePRFMarianFit(histoprf);
	      Int_t nbtg = 6*4*18*((Int_t) ((AliTRDCalibraMode *)calibra->GetCalibraMode())->GetDetChamb0(2))
		+ 6*  18*((Int_t) ((AliTRDCalibraMode *)calibra->GetCalibraMode())->GetDetChamb2(2));
	      Int_t nbfit        = calibra->GetNumberFit();
	      Int_t nbE          = calibra->GetNumberEnt();
	      if ((nbtg >                  0) && 
		  (nbfit        >= 0.95*nbE)) {
		TObjArray object              = calibra->GetVectorFit();
		TObject *objPRFpad          = calibra->CreatePadObjectPRF(&object);
		if(!Store("Calib","PRFWidth"         ,(TObject *) objPRFpad          ,md2,0,kTRUE)){
		  Log("Error storing the calibration object for the Pad Response Function");
		  result |= kEStoreCalHLT;
		}
	      }
	      else{
		Log("Not enough statistics for the Pad Response Function");
		result |= kEFitHistoHLT;
	      }
	      calibra->ResetVectorFit();
	    }
	  } // if HLT openfile
	} // while iter list
      } // if HLT size of tlist
    } //if HLT tlist
    delete filesourceshlt;
  } //if run type physics
  // time
  timer.Stop();
  timer.Print();

  //
  // Process the calibration data for the DAQ part
  //

  // How long does it take for the DAQ part?
  timer.Reset();
  timer.Start();
  
  // Take the file from the DAQ file exchange server
  TList *filesourcesdaq = GetFileSources(kDAQ,"PADSTATUSVDRIFT");
  if (!filesourcesdaq) {
    Log(Form("No sources found for GAINDRIFTPRF for run %d !",fRun));
    result |= kEListFileDAQ;
  }
  else{
    if (filesourcesdaq->GetSize() != 1) {
      Log(Form("More than one source found for PADSTATUSVDRIFT for run %d!",fRun));
      filesourcesdaq->Print();
      result |= kEListFileDAQ;
    }
    else{
      // Loop over the files taken from the DAQ
      TIter iter(filesourcesdaq);
      TObjString *sourcedaq;
      while ((sourcedaq = dynamic_cast<TObjString *> (iter.Next()))) {
    
	TString filenamedaq = GetFile(kDAQ,"PADSTATUSVDRIFT",sourcedaq->GetName());
	if (filenamedaq.Length() == 0) {
	  Log(Form("Error retrieving file from source %d failed!", sourcedaq->GetName()));
	  result |= kEOpenFileDAQ;
	}
	else {

	  // Init some things
	  TProfile2D   *histodriftvelocity  = 0x0; // Profile taken from DAQ for drift velocity and time0
	  AliTRDCalibPadStatus *calibpadstatus = 0x0; // AliTRDCalibPadStatus from DAQ for pad status
	  
	  // Take the histos
	  Bool_t something = kFALSE;
	  TFile *filedaq = TFile::Open(filenamedaq);
	  calibpadstatus = (AliTRDCalibPadStatus *) filedaq->Get("calibpadstatus");
	  if (!calibpadstatus) {
	    Log("No pedetral run!");
	  }
	  else something = kTRUE;
	  histodriftvelocity = (TProfile2D *) filedaq->Get("PH2d");
	  if (!histodriftvelocity) {
	    Log("No Vdrift TProfile2D!");
	  }
	  else{
	    histodriftvelocity->SetDirectory(0);
	  }
	  if(histodriftvelocity) something = kTRUE;
	  if(!something){
	    Log("Error DAQ, nothing in the file!");
	    result |= kETakeObjectDAQ;
	  }
	 	  
	  // try to fit and store reference data
	  if(calibpadstatus){
	    calibpadstatus->AnalyseHisto();
	    if(!StoreReferenceData("DAQData","PadStatus",(TObject *) calibpadstatus,&metaData)){
	      Log("Error storing AliTRDCalibPadStatus object as reference data");
	      result |= kEStoreRefDAQ;
	    }
	    AliTRDCalPadStatus *calPadStatus = calibpadstatus->CreateCalPadStatus();
	    AliCDBMetaData *md3= new AliCDBMetaData(); 
	    md3->SetObjectClassName("AliTRDCalPadStatus");
	    md3->SetResponsible("Raphaelle Bailhache");
	    md3->SetBeamPeriod(1);
       	    md3->SetComment("TRD calib test");
	    if(!Store("Calib","PadStatus"    ,(TObject *)calPadStatus, md3, 0, kTRUE)){
	      Log("Error storing the calibration object for the chamber vdrift");
	      result |= kEStoreCalDAQ;    
	    }
	  }
	  if (histodriftvelocity) {
	    if(!StoreReferenceData("DAQData","VdriftT0",(TObject *) histodriftvelocity,&metaData)){
	      Log("Error storing 2D Profile for average pulse height as reference data");
	      result |= kEStoreRefDAQ;
	    }
	    if(!hltVdrifthisto){
	      Log("Take the PH reference data. Now we will try to fit\n");
	      calibra->SetMinEntries(100*20); // If there is less than 2000
	      calibra->AnalysePH(histodriftvelocity);
	      Int_t nbtg = 6*4*18*((Int_t) ((AliTRDCalibraMode *)calibra->GetCalibraMode())->GetDetChamb0(1))
		+ 6*  18*((Int_t) ((AliTRDCalibraMode *)calibra->GetCalibraMode())->GetDetChamb2(1));
	      Int_t nbfit        = calibra->GetNumberFit();
	      Int_t nbE        = calibra->GetNumberEnt();
	      if ((nbtg >                  0) && 
		  (nbfit        >= 0.95*nbE)) {
		TObjArray object      = calibra->GetVectorFit();
		AliTRDCalDet *objdriftvelocitydet = calibra->CreateDetObjectVdrift(&object,kTRUE);
		TObject *objdriftvelocitypad = calibra->CreatePadObjectVdrift();
		object              = calibra->GetVectorFit2();
		AliTRDCalDet *objtime0det         = calibra->CreateDetObjectT0(&object,kTRUE);
		TObject *objtime0pad         = calibra->CreatePadObjectT0();
		calibra->ResetVectorFit();
		if(!Store("Calib","ChamberVdrift"    ,(TObject *) objdriftvelocitydet,md1,0,kTRUE)){
		  Log("Error storing the calibration object for the chamber vdrift");
		  result |= kEStoreCalDAQ;    
		}
		if(!Store("Calib","ChamberT0"        ,(TObject *) objtime0det        ,md1,0,kTRUE)){
		  Log("Error storing the calibration object for the chamber t0");
		  result |= kEStoreCalDAQ;    
		}
		if(!Store("Calib","LocalVdrift"      ,(TObject *) objdriftvelocitypad,md2,0,kTRUE)){
		  Log("Error storing the calibration object for the local drift velocity");
		  result |= kEStoreCalDAQ;
		}
		if(!Store("Calib","LocalT0"          ,(TObject *) objtime0pad        ,md2,0,kTRUE)){
		  Log("Error storing the calibration object for the local time0");
		  result |= kEStoreCalDAQ;
		}
	      }
	      else{
		Log("Not enough statistics for the average pulse height");
		result |= kEFitObjectDAQ;
	      }
	    }
	  }
	  filedaq->Close();
	}// if DAQ open file
      } // while iter DAQ list
    } // size DAQ list
  } // tlist
  delete filesourcesdaq;
  // time
  timer.Stop();
  timer.Print();
  
  return result;  

}
