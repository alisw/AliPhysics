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



#include "AliACORDEPreprocessor.h"
#include "TRandom.h"
#include "TFile.h"
#include "AliCDBMetaData.h"
#include "AliCDBEntry.h"
#include "AliLog.h"
#include "AliACORDECalibData.h"

#include <TTimeStamp.h>
#include <TObjString.h>
#include <TList.h>
#include <TH1F.h>

//
// This is the first version of ACORDE Preprocessor
// It takes data from DAQ and passes it to the class AliACORDECalibModule and 
// stores reference data.
//
// Authors
// Pedro Gonzalez pedro.gonzalez@fcfm.buap.mx
// Irais Bautista irais@fcfm.buap.mx
// Arturo Fernandez Tellez afernan@cern.ch

ClassImp(AliACORDEPreprocessor)

//______________________________________________________________________________________________
AliACORDEPreprocessor::AliACORDEPreprocessor(AliShuttleInterface* shuttle) :
  AliPreprocessor("ACO", shuttle),
  fCalData(0)
{
  // constructor
}

//______________________________________________________________________________________________
AliACORDEPreprocessor::~AliACORDEPreprocessor()
{
  // destructor
}

//______________________________________________________________________________________________
void AliACORDEPreprocessor::Initialize(Int_t run, UInt_t startTime,
	UInt_t endTime)
{
  // Creates AliACORDECalibModule object

  AliPreprocessor::Initialize(run, startTime, endTime);

	Log(Form("\n\tRun %d \n\tStartTime %s \n\tEndTime %s", run,
		TTimeStamp(startTime).AsString(),
		TTimeStamp(endTime).AsString()));

        fCalData = new AliACORDECalibData();
}

//______________________________________________________________________________________________
UInt_t AliACORDEPreprocessor::Process(TMap* /*dcsAliasMap*/)
{
  

 
 
  TH1D *fH1,*fH2,*fH3,*fH4; //Histogram of the rates per module
   TFile *daqFile=0x0;


   // retrieve the run type from the Shuttle,

 
   TString runType = GetRunType();


   //acorde	STANDALONE_BC
   //acorde	STANDALONE_PULSER
    
   if(runType !="STANDALONE_PULSER")
   {

   Log("RunType is not STANDALONE_PULSER, nothing to do");
   return 1;

   }


   Log(Form("Run type for run %d: %s", fRun, runType.Data()));
   TString SourcesId = "CALIB";

  

   //retrieve the list of sources that produced the file with id RATES
 
    TList* sourceList = GetFileSources(kDAQ,SourcesId.Data());

   if (!sourceList)
   {
  	Log(Form("Error: No sources found for id %s", SourcesId.Data()));
	return 2;
   }
  
   // TODO We have the list of sources that produced the files with Id RATES 
   // Now we will loop on the list and we'll query the files one by one. 



   Log(Form("The following sources produced files with the id %s",SourcesId.Data()));
   sourceList->Print();
     
   TIter iter(sourceList);
   TObjString *source = 0;
 
 

   while((source=dynamic_cast<TObjString*> (iter.Next())))
   {
  	
        TString fileName = GetFile(kDAQ,SourcesId.Data(), source->GetName());

  	if (fileName.Length() > 0)
   		Log(Form("Got the file %s, now we can extract some values.", fileName.Data()));

                daqFile = new TFile(fileName.Data(),"READ");
              
              if(!daqFile)
              {
                            
              Log(Form("There are not histos     1"));
	      return 3;

              }

           
                
            fH1 = (TH1D*)daqFile->Get("fHist1");
            fH2 = (TH1D*)daqFile->Get("fHist2");
            fH3 = (TH1D*)daqFile->Get("fHist3");
            fH4 = (TH1D*)daqFile->Get("fHist4");
          

             
             if(fH1!=NULL&&fH2!=NULL&&fH3!=NULL&&fH4!=NULL)  
             {  
             fCalData->AddHHits(fH1);
             fCalData->AddHTHits(fH2);
             fCalData->AddHMultiHits(fH3);
             fCalData->AddHTMultiHits(fH4);
             }
            
             else
            {
             Log(Form("There are not histos     2"));
             return 4;
            }
    

   }                   
                                                                          
 
 
  delete sourceList;
  

        //Now we have to store

        AliCDBMetaData metaData;
	metaData.SetBeamPeriod(0);
	metaData.SetResponsible("Pedro and Irais");
	metaData.SetComment("This preprocessor fills an AliACORDECalibModule object.");

          	Bool_t result = StoreReferenceData("Calib", "Data",fCalData, &metaData);
      
	delete fCalData;
	fCalData = 0;
  

  if (!result)
  return 4;
  
  return 0;
}

