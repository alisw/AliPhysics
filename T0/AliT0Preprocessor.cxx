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
$Log$
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
// 1) takes data from DCS and passes it to the class AliTOFDataDCS for processing and writes the result to the Reference DB.
// 2) takes data form DAQ (both from Laser Calibration and Physics runs), processes it, and stores either to OCDB or to Reference DB.

#include "AliT0Preprocessor.h"
#include "AliT0DataDCS.h"
#include "AliT0CalibWalk.h"
#include "AliT0CalibTimeEq.h"

#include "AliCDBMetaData.h"
#include "AliDCSValue.h"
#include "AliLog.h"

#include <TTimeStamp.h>
#include <TFile.h>
#include <TObjString.h>
#include <TNamed.h>
#include "AliT0Dqclass.h"


ClassImp(AliT0Preprocessor)

//____________________________________________________
AliT0Preprocessor::AliT0Preprocessor(AliShuttleInterface* shuttle) : AliPreprocessor("T00", shuttle), fData(0)
{
  //constructor
}
//____________________________________________________

AliT0Preprocessor::~AliT0Preprocessor()
{
  delete fData;
  fData = 0;
}
//____________________________________________________

void AliT0Preprocessor::Initialize(Int_t run, UInt_t startTime, UInt_t endTime)
{
  AliPreprocessor::Initialize(run, startTime, endTime);
  AliInfo(Form("\n\tRun %d \n\tStartTime %s \n\tEndTime %s", run, TTimeStamp(startTime).AsString(), TTimeStamp(endTime).AsString()));
  fData = new AliT0DataDCS(fRun, fStartTime, fEndTime);
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

	Bool_t resultDCSMap=kFALSE;
	Bool_t resultDCSStore=kFALSE;
	Bool_t resultLaser=kFALSE;
	Bool_t resultOnline=kFALSE;  
     
        if(!dcsAliasMap)
        {
          Log("No DCS input data");
          return 1;
        }
        else
        {
	  /*
          resultDCSMap=fData->ProcessData(*dcsAliasMap);
          if(!resultDCSMap)
          {
            Log("Error when processing DCS data");
            return 2;// return error Code for processed DCS data not stored
          }
          else
          {
	    Float_t *meanScaler[24] = fData->GetScalerMean();
	     	
            AliCDBMetaData metaDataDCS;
            metaDataDCS.SetBeamPeriod(0);
            metaDataDCS.SetResponsible("Tomasz Malkiewicz");
            metaDataDCS.SetComment("This preprocessor fills an AliTODataDCS object.");
            AliInfo("Storing DCS Data");
            resultDCSStore = Store("Calib","DCSData",meanScaler, &metaDataDCS);
            if (!resultDCSStore)
            {
              Log("Some problems occurred while storing DCS data results in ReferenceDB");
              return 2;// return error Code for processed DCS data not stored
            }

          }
	  */
        }

        // processing DAQ

        TString runType = GetRunType();

        if(runType == "T0_STANDALONE_LASER")
        {
          TList* list = GetFileSources(kDAQ, "LASER");
          if (list)
          {
            TIter iter(list);
            TObjString *source;
            while ((source = dynamic_cast<TObjString *> (iter.Next())))
            {
              const char *laserFile = GetFile(kDAQ, "LASER", source->GetName());
              if (laserFile)
              {
                Log(Form("File with Id LASER found in source %s!", source->GetName()));
                AliT0CalibWalk *laser = new AliT0CalibWalk();
                laser->MakeWalkCorrGraph(laserFile);
                AliCDBMetaData metaData;
                metaData.SetBeamPeriod(0);
                metaData.SetResponsible("Tomek&Michal");
                metaData.SetComment("Walk correction from laser runs.");
		//TObjArray* arrLaser = laser->GetfWalk();
		resultLaser=Store("Calib","Slewing_Walk", laser, &metaData, 0, 1);
                delete laser;
              }
              else
              {
                Log(Form("Could not find file with Id LASER in source %s!", source->GetName()));
                return 1;
              }
            }
            if (!resultLaser)
            {
              Log("No Laser Data stored");
              return 3;//return error code for failure in storing Laser Data
            }
          } else {
	  	Log("No sources found for id LASER!");
		return 1;
	  }
        }
        else if(runType == "PHYSICS")
        {
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
                online->ComputeOnlineParams("CFD", 20, 4., filePhys);
                AliCDBMetaData metaData;
                metaData.SetBeamPeriod(0);
                metaData.SetResponsible("Tomek&Michal");
                metaData.SetComment("Time equalizing result.");
                resultOnline = Store("Calib","TimeDelay", online, &metaData, 0, 1);
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
              Log("No Laser Data stored");
              return 4;//return error code for failure in storing OCDB Data
            }
          } else {
	  	Log("No sources found for id PHYSICS!");
		return 1;
	  }
        }

  return 0;
}

