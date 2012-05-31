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

/*
$Log: AliT0Preprocessor.cxx,v $
Revision 1.8  2007/12/07 15:22:51  alla
bug fixed by Alberto
 
Revision 1.7  2007/12/06 16:35:24  alla
new bugs fixed by Tomek

Revision 1.5  2007/11/23 19:28:52  alla
bug fixed

Version 2.1  2007/11/21 
Preprocessor storing data to OCDB (T.Malkiewicz)

Version 1.1  2006/10   
Preliminary test version (T.Malkiewicz)
*/   
// T0 preprocessor:
// 1) takes data from DCS and passes it to the class AliTOFDataDCS 
// for processing and writes the result to the Reference DB.
// 2) takes data form DAQ (both from Laser Calibration and Physics runs), 
// processes it, and stores either to OCDB or to Reference DB.


#include "AliT0Preprocessor.h"
#include "AliT0DataDCS.h"
#include "AliT0CalibWalk.h"
#include "AliT0CalibTimeEq.h"

#include "AliCDBMetaData.h"
#include "AliDCSValue.h"
#include "AliCDBEntry.h"
#include "AliLog.h"

#include <TTimeStamp.h>
#include <TFile.h>
#include <TObjString.h>
#include <TNamed.h>
#include "AliT0Dqclass.h"
#include "TClass.h"


ClassImp(AliT0Preprocessor)

//____________________________________________________
AliT0Preprocessor::AliT0Preprocessor(AliShuttleInterface* shuttle) : 
  AliPreprocessor("T00", shuttle), 
  fData(0)
{
  //constructor
  AddRunType("PHYSICS");
  AddRunType("STANDALONE");
  AddRunType("AMPLITUDE_CALIBRATION");
}
//____________________________________________________

AliT0Preprocessor::~AliT0Preprocessor()
{
  //destructor
  delete fData;
  fData = 0;
}
//____________________________________________________

void AliT0Preprocessor::Initialize(Int_t run, UInt_t startTime, UInt_t endTime)
{
  // Creates AliT0DataDCS object
  AliPreprocessor::Initialize(run, startTime, endTime);
  AliInfo(Form("\n\tRun %d \n\tStartTime %s \n\tEndTime %s", run, TTimeStamp(startTime).AsString(), TTimeStamp(endTime).AsString()));
  fData = new AliT0DataDCS(fRun, fStartTime, fEndTime, GetStartTimeDCSQuery(), GetEndTimeDCSQuery());
}
//____________________________________________________

Bool_t AliT0Preprocessor::ProcessDCS(){
	// Check whether DCS should be processed or not...
	TString runType = GetRunType();
	Log(Form("ProcessDCS - RunType: %s",runType.Data()));

	if(runType == "PHYSICS" )
	  return kTRUE;
	else
	  return kFALSE;
	
}
//____________________________________________________

UInt_t AliT0Preprocessor::ProcessDCSDataPoints(TMap* dcsAliasMap){
	// Fills data into AliT0DataDCS object
	Log("Processing DCS DP");
	Bool_t resultDCSMap=kFALSE;
	Bool_t resultDCSStore=kFALSE;

        if(!dcsAliasMap)
        {
          Log("No DCS input data");
          return 1;
        }
        else
        {
          resultDCSMap=fData->ProcessData(*dcsAliasMap);
          if(!resultDCSMap)
          {
            Log("Error when processing DCS data");
            return 2;// return error Code for processed DCS data not stored
          }
          else
          {
            AliCDBMetaData metaDataDCS;
            metaDataDCS.SetBeamPeriod(0);
            metaDataDCS.SetResponsible("Tomasz Malkiewicz");
            metaDataDCS.SetComment("This preprocessor fills an AliTODataDCS object.");
            AliInfo("Storing DCS Data");
            resultDCSStore = StoreReferenceData("Calib","DCSData",fData, &metaDataDCS);
            if (!resultDCSStore)
            {
              Log("Some problems occurred while storing DCS data results in ReferenceDB");
              return 2;// return error Code for processed DCS data not stored
            }
          }
        }
	return 0;
}
//____________________________________________________

UInt_t AliT0Preprocessor::ProcessLaser()
{
  // Processing data from DAQ Standalone run
  Log("Processing Laser calibration - Walk Correction");
  Bool_t resultLaser = kFALSE;
  Bool_t writeok = kFALSE;
  //processing DAQ
  TList* list = GetFileSources(kDAQ, "AMPLITUDE_CALIBRATION");
  AliT0CalibWalk *laser = new AliT0CalibWalk();
  TObjString *source;
  if (list)
    {
      TIter iter(list);
      while ((source = dynamic_cast<TObjString *> (iter.Next())))
	{
	  const char *laserFile = GetFile(kDAQ, "AMPLITUDE_CALIBRATION", source->GetName());
	  if (laserFile)
	    {
	      Log(Form("File with Id AMPLITUDE_CALIBRAION found in source %s!", source->GetName()));
	     writeok = laser->MakeWalkCorrGraph(laserFile);
	      
	    }
	}
      
      AliCDBMetaData metaData;
      metaData.SetBeamPeriod(0);
      metaData.SetResponsible("Alla");
      metaData.SetComment("Walk correction from laser runs.");
      if (writeok) resultLaser=Store("Calib","Slewing_Walk", laser, &metaData, 0, 1);
      else {
 	
	Log(Form("writeok = %d no peaks in CFD spectra",writeok));
	return 0;
      }		  
      Log(Form("resultLaser = %d",resultLaser));
      if (!resultLaser)
	{
	  Log("No Laser Data stored");
	  return 3;//return error code for failure in storing Laser Data
	}
    }
  else
    {
      Log(Form("Could not find file with Id  AMPLITUDE_CALIBRAION "));
      return 1;
    }
  
   return 0;
}

//____________________________________________________

UInt_t AliT0Preprocessor::ProcessPhysics(){
	//Processing data from DAQ Physics run
	Log("Processing Physics");
	
	Bool_t resultOnline=kFALSE; 
	//processing DAQ
	TList* listPhys = GetFileSources(kDAQ, "PHYSICS");
        if (listPhys)
          {
            TIter iter(listPhys);
            TObjString *sourcePhys;
            while ((sourcePhys = dynamic_cast<TObjString *> (iter.Next())))
            {
              const char *filePhys = GetFile(kDAQ, "PHYSICS", sourcePhys->GetName());
              if (filePhys)
              {
                AliT0CalibTimeEq *online = new AliT0CalibTimeEq();
                online->Reset();
                Bool_t writeok = online->ComputeOnlineParams(filePhys);
                AliCDBMetaData metaData;
                metaData.SetBeamPeriod(0);
                metaData.SetResponsible("Alla Maevskaya");
                metaData.SetComment("Time equalizing result.");

                if (writeok) resultOnline = Store("Calib","TimeDelay", online, &metaData, 0, 1);
		else {
		  
		  Log(Form("writeok = %d not enough data for equalizing",resultOnline));
		  return 0;
		}		  
		Log(Form("resultOnline = %d",resultOnline));
                delete online;
              }
		else
              {
                Log(Form("Could not find file with Id PHYSICS in source %s!", sourcePhys->GetName()));
                return 1;
              }
              
            }
            if (!resultOnline)
            {
              Log("No Data stored");
              return 4;//return error code for failure in storing OCDB Data
            }
          } else {
	  	Log("No sources found for id PHYSICS!");
		return 1;
	  }
	return 0;
}
//____________________________________________________

UInt_t AliT0Preprocessor::Process(TMap* dcsAliasMap )
{
  // T0 preprocessor return codes:
  // return=0 : all ok
  // return=1 : no DCS input data 
  // return=2 : failed to store DCS data
  // return=3 : no Laser data (Walk correction)
  // return=4 : failed to store OCDB time equalized data
  // return=5 : no DAQ input for OCDB
  // return=6 : failed to retrieve DAQ data from OCDB
  // return=7 : failed to store T0 OCDB data
  // return=8 : not enough data for equalizing
  Bool_t dcsDP = ProcessDCS();
  Log(Form("dcsDP = %d",dcsDP));	
  TString runType = GetRunType();
  Log(Form("RunType: %s",runType.Data()));
  //processing
  if(runType == "STANDALONE"){
    if(dcsDP==1){
      Int_t iresultDCS = ProcessDCSDataPoints(dcsAliasMap);
      return iresultDCS;
    }
  }
  
  if(runType == "AMPLITUDE_CALIBRATION"){
    Int_t iresultLaser = ProcessLaser();
    if(dcsDP==1){
      Int_t iresultDCS = ProcessDCSDataPoints(dcsAliasMap);
      return iresultDCS;
    }
    
    Log(Form("iresultLaser = %d",iresultLaser));
    return iresultLaser;
  }
  
  else if(runType == "PHYSICS"){
    Int_t iresultPhysics = ProcessPhysics();
    if(dcsDP==1){
      Int_t iresultDCS = ProcessDCSDataPoints(dcsAliasMap);
      return iresultDCS;
    }
    Log(Form("iresultPhysics = %d",iresultPhysics));
	  return iresultPhysics; 
      }
  
	
	
	return 0;
}
