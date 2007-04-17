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
//                                                                        //
// Author:                                                                //
//   R. Bailhache (R.Bailhache@gsi.de)                                    //
//                                                                        //
////////////////////////////////////////////////////////////////////////////

#include "AliTRDPreprocessor.h"

#include <TTimeStamp.h>
#include <TFile.h>
#include <TProfile2D.h>
#include <TH2I.h>
#include <TStopwatch.h>
#include <TObjString.h>
#include <TString.h>
#include <TList.h>
#include <TCollection.h>

#include "AliCDBMetaData.h"
#include "AliDCSValue.h"
#include "AliLog.h"

#include "AliTRDCalibraFit.h"
#include "AliTRDCalibraMode.h"
#include "Cal/AliTRDCalDet.h"

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
UInt_t AliTRDPreprocessor::Process(TMap* /*dcsAliasMap*/)
{
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
    Log("Nothing to do!");
    return 0;
  }


  // note that the parameters are returned as character strings!
  const char* nEvents = GetRunParameter("totalEvents");
  if (nEvents) {
    Log(Form("Number of events for run %d: %s",fRun, nEvents));
  } else {
    Log(Form("Number of events not put in logbook!"));
  }


  // Metadata for the reference data
  AliCDBMetaData metaData;
  metaData.SetBeamPeriod(1);
  metaData.SetResponsible("Raphaelle Bailhache");
  metaData.SetComment("This preprocessor fills reference data.");

  // Take the file from the HLT file exchange server
  TList *filesources = GetFileSources(kHLT,"GAINDRIFTPRF");
  if (!filesources) {
    Log(Form("No sources found for GAINDRIFTPRF for run %d !",fRun));
    return 1;
  }
  if (filesources->GetSize() != 1) {
    Log(Form("More than one source found for GAINDRIFTPRF for run %d!",fRun));
    filesources->Print();
    delete filesources;
    return 1;
  }

  // Call a AliTRDCalibra instance for fit
  AliTRDCalibraFit *calibra = AliTRDCalibraFit::Instance();

  //Choose the fit methods
  calibra->SetFitChargeNDB(4); //for the relative gain
  calibra->SetFitMeanWOn();   //weighted mean
  calibra->SetFitPHNDB(3);    //for the average pulse height
  calibra->SetFitLagrPolOn(); //LagrPol
  calibra->SetFitPRFNDB(0);   //for the PRF
  calibra->SetFitPRFOn();     //gaussian fit

  //Debug mode
  //calibra->SetDebug(1);       //Debug

  // Init some things
  AliTRDCalDet *objgaindet          = 0x0; // Object for det average gain factor
  AliTRDCalDet *objdriftvelocitydet = 0x0; // Object for det average drift velocity
  AliTRDCalDet *objtime0det         = 0x0; // Object for det average time0 
  TObject      *objgainpad          = 0x0; // Object for pad (relative to the det) gain factor
  TObject      *objdriftvelocitypad = 0x0; // Object for pad (relative to the det) drift velocity
  TObject      *objtime0pad         = 0x0; // Object for pad (relative to the det) time0
  TObject      *objPRFpad           = 0x0; // Object for pad prf width
  TH2I         *histogain           = 0x0; // Histogram taken from HLT for gain factor
  TProfile2D   *histodriftvelocity  = 0x0; // Profile taken from HLT for drift velocity and time0
  TProfile2D   *histoprf            = 0x0; // Profile taken from HLT for prf

  Int_t    numberfit[3]        = { 0,   0,   0   }; // Number of histos fitted for gain, drift velocity and prf
  Int_t    numberEnt[3]        = { 0,   0,   0   }; // Number of histos with entries
  Double_t statisticmean[3]    = { 0.0, 0.0, 0.0 }; // Mean values of the number of entries in these histos
  Int_t    numbertotalgroup[3] = { 0,   0,   0   }; // Total number of groups

  // Loop over the files taken from the HLT
  TIter iter(filesources);
  TObjString *source;
  while ((source = dynamic_cast<TObjString *> (iter.Next()))) {
    
    TString filename = GetFile(kHLT,"GAINDRIFTPRF",source->GetName());
    if (filename.Length() == 0) {
      Log(Form("Error retrieving file from source %d failed!", source->GetName()));
      delete filesources;
      return 2;
    }

    // Take the histos
    TFile *file = TFile::Open(filename);
    histogain = (TH2I *) file->Get("CH2d");
    histogain->SetDirectory(0);
    if (!histogain) {
      Log("Error retrieving 2D histos for gain failed!");
      delete filesources;
      return 2;
    }
    histodriftvelocity = (TProfile2D *) file->Get("PH2d");
    histodriftvelocity->SetDirectory(0);
    if (!histodriftvelocity) {
      Log("Error retrieving 2D Profile for average pulse height failed!");
      delete filesources;
      return 2;
    }
    histoprf = (TProfile2D *) file->Get("PRF2d");
    histoprf->SetDirectory(0);
    if (!histoprf) {
      Log("Error retrieving 2D Profile for Pad Response Function failed!");
      delete filesources;
      return 2;
    }
    file->Close();

    // Set the mode of calibration from the TObject, store the reference data and try to fit them
    if (histogain) {
      calibra->SetModeCalibrationFromTObject((TObject *) histogain,0);
      if(!StoreReferenceData("HLTData","Gain",(TObject *) histogain,&metaData)){
	Log("Error storing 2D histos for gain as reference data");
	delete filesources;
	return 3;
      }
      Log("Take the CH reference data. Now we will try to fit\n");
      calibra->SetMinEntries(100); // If there is less than 100 entries in the histo: no fit
      calibra->FitCHOnline(histogain);
      numbertotalgroup[0] = 6*4*18*((Int_t) ((AliTRDCalibraMode *)calibra->GetCalibraMode())->GetDetChamb0(0))
      	                  + 6*  18*((Int_t) ((AliTRDCalibraMode *)calibra->GetCalibraMode())->GetDetChamb2(0));
      numberfit[0]        = calibra->GetNumberFit();
      statisticmean[0]    = calibra->GetStatisticMean(); 
      numberEnt[0]        = calibra->GetNumberEnt();
      objgaindet          = calibra->CreateDetObjectTree(calibra->GetGain(),0);
      objgainpad          = calibra->CreatePadObjectTree(calibra->GetGain(),0,objgaindet);
    }
    
    if (histodriftvelocity) {
      calibra->SetModeCalibrationFromTObject((TObject *) histodriftvelocity,1);
      if(!StoreReferenceData("HLTData","VdriftT0",(TObject *) histodriftvelocity,&metaData)){
	Log("Error storing 2D Profile for average pulse height as reference data");
	delete filesources;
	return 3;
      }
      Log("Take the PH reference data. Now we will try to fit\n");
      calibra->SetMinEntries(100*20); // If there is less than 2000
      calibra->FitPHOnline(histodriftvelocity);
      numbertotalgroup[1] = 6*4*18*((Int_t) ((AliTRDCalibraMode *)calibra->GetCalibraMode())->GetDetChamb0(1))
	                  + 6*  18*((Int_t) ((AliTRDCalibraMode *)calibra->GetCalibraMode())->GetDetChamb2(1));
      numberfit[1]        = calibra->GetNumberFit();
      statisticmean[1]    = calibra->GetStatisticMean(); 
      numberEnt[1]        = calibra->GetNumberEnt();
      objdriftvelocitydet = calibra->CreateDetObjectTree(calibra->GetVdrift(),1);
      objdriftvelocitypad = calibra->CreatePadObjectTree(calibra->GetVdrift(),1,objdriftvelocitydet);
      objtime0det         = calibra->CreateDetObjectTree(calibra->GetT0(),3);
      objtime0pad         = calibra->CreatePadObjectTree(calibra->GetT0(),3,objtime0det);
    }
    
    if (histoprf) {
      calibra->SetModeCalibrationFromTObject((TObject *) histoprf,2);
      if(!StoreReferenceData("HLTData","PRF",(TObject *) histoprf,&metaData)){
	Log("Error storing the 2D Profile for Pad Response Function as reference data");
	delete filesources;
	return 3;
      }
      Log("Take the PRF reference data. Now we will try to fit\n");
      calibra->SetMinEntries(100*20); // If there is less than 2000
      calibra->SetRangeFitPRF(0.5);
      calibra->FitPRFOnline(histoprf);
      numbertotalgroup[2] = 6*4*18*((Int_t) ((AliTRDCalibraMode *)calibra->GetCalibraMode())->GetDetChamb0(2))
                           + 6*  18*((Int_t) ((AliTRDCalibraMode *)calibra->GetCalibraMode())->GetDetChamb2(2));
      numberfit[2]        = calibra->GetNumberFit();
      statisticmean[2]    = calibra->GetStatisticMean(); 
      numberEnt[2]        = calibra->GetNumberEnt();
      objPRFpad           = calibra->CreatePadObjectTree(calibra->GetPRF());
    }

  }
  
  // Bilan of the fit statistic
  Log(Form("The mean number of entries required for a fit is: %d"
              ,(Int_t) calibra->GetMinEntries()));
  Log(Form("FOR THE CH: There is a mean statistic of: %f, with %d fits for %d groups and %d histos with entries"
              ,statisticmean[0],numberfit[0],numbertotalgroup[0],numberEnt[0]));
  Log(Form("FOR THE PH: There is a mean statistic of: %f, with %d fits for %d groups and %d histos with entries"
              ,statisticmean[1],numberfit[1],numbertotalgroup[1],numberEnt[1]));
  Log(Form("FOR THE PRF: There is a mean statistic of: %f, with %d fits for %d groups and %d histos with entries"
              ,statisticmean[2],numberfit[2],numbertotalgroup[2],numberEnt[2]));
  
  
  //
  // Store the coefficients in the grid OCDB if enough statistics
  //
  
  // Store the infos for the detector
  AliCDBMetaData *md1= new AliCDBMetaData(); 
  md1->SetObjectClassName("AliTRDCalDet");
  md1->SetResponsible("Raphaelle Bailhache");
  md1->SetBeamPeriod(1);
  md1->SetAliRootVersion("01-10-06"); // root version
  md1->SetComment("The dummy values in this calibration file are for testing only");
  // Gain
  if ((numbertotalgroup[0] >                  0) && 
      (numberfit[0]        >= 0.95*numberEnt[0])) {
    if(!Store("Calib","ChamberGainFactor",(TObject *) objgaindet         ,md1,0,kTRUE)){
      Log("Error storing the calibration object for the chamber gain");
      delete filesources;
      return 4;
    }
  }
  else{
    Log("Not enough statistics for the gain");
  }
  // Vdrift and time0
  if ((numbertotalgroup[1] >                  0) && 
      (numberfit[1]        >= 0.95*numberEnt[1])) {
    if(!Store("Calib","ChamberVdrift"    ,(TObject *) objdriftvelocitydet,md1,0,kTRUE)){
      Log("Error storing the calibration object for the chamber vdrift");
      delete filesources;
      return 4;
    }
    if(!Store("Calib","ChamberT0"        ,(TObject *) objtime0det        ,md1,0,kTRUE)){
      Log("Error storing the calibration object for the chamber t0");
      delete filesources;
      return 4;
    }
  }
  else{
    Log("Not enough statistics for the average pulse height");
  }
  
  // Store the infos for the pads
  AliCDBMetaData *md2= new AliCDBMetaData(); 
  md2->SetObjectClassName("AliTRDCalPad");
  md2->SetResponsible("Raphaelle Bailhache");
  md2->SetBeamPeriod(1);
  md2->SetAliRootVersion("01-10-06"); //root version
  md2->SetComment("The dummy values in this calibration file are for testing only");
  // Gain
  if ((numbertotalgroup[0] >                  0) && 
      (numberfit[0]        >= 0.95*numberEnt[0])) {
    if(!Store("Calib","LocalGainFactor"  ,(TObject *) objgainpad         ,md2,0,kTRUE)){
      Log("Error storing the calibration object for the local gain factor");
      delete filesources;
      return 4;
    }
  }
  // Vdrift and time0
  if ((numbertotalgroup[1] >                  0) && 
      (numberfit[1]        >= 0.95*numberEnt[1])) {
    if(!Store("Calib","LocalVdrift"      ,(TObject *) objdriftvelocitypad,md2,0,kTRUE)){
      Log("Error storing the calibration object for the local drift velocity");
      delete filesources;
      return 4;
    }
    if(!Store("Calib","LocalT0"          ,(TObject *) objtime0pad        ,md2,0,kTRUE)){
      Log("Error storing the calibration object for the local time0");
      delete filesources;
      return 4;
    }
  }
  // Pad Response Width
  if ((numbertotalgroup[2] >                  0) && 
      (numberfit[2]        >= 0.95*numberEnt[2])) {
    if(!Store("Calib","PRFWidth"         ,(TObject *) objPRFpad          ,md2,0,kTRUE)){
      Log("Error storing the calibration object for the Pad Response Function");
      delete filesources;
      return 4;
    }
  }
  else{
    Log("Not enough statistics for the Pad Response Function");
  }
  
  // End
  delete filesources;
  timer.Stop();
  timer.Print();
  return 0;  

}
