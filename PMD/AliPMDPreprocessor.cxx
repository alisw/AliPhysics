/***************************************************************************
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
//-----------------------------------------------------//
//                                                     //
//  Source File : AliPMDPreprocessor.cxx               //
//                                                     //
//                                                     //
//-----------------------------------------------------//

// --- ROOT system
#include <TFile.h>
#include <TTimeStamp.h>
#include <TObjString.h>
#include <TTree.h>

#include "AliLog.h"
#include "AliShuttleInterface.h"
#include "AliCDBMetaData.h"
#include "AliPMDCalibData.h"
#include "AliPMDPedestal.h"
#include "AliPMDPreprocessor.h"


ClassImp(AliPMDPreprocessor)
  
//______________________________________________________________________________________________
AliPMDPreprocessor::AliPMDPreprocessor(AliShuttleInterface* shuttle) :
  AliPreprocessor("PMD", shuttle)
{
  // constructor
}

//______________________________________________________________________________________________
AliPMDPreprocessor::~AliPMDPreprocessor()
{
  // destructor
}

//______________________________________________________________________________________________
void AliPMDPreprocessor::Initialize(Int_t run, UInt_t startTime,
	UInt_t endTime)
{
  // Creates AliPMDDataDAQ object

    AliPreprocessor::Initialize(run, startTime, endTime);

    AliInfo(Form("\n\tRun %d \n\tStartTime %s \n\tEndTime %s", run,
		 TTimeStamp(startTime).AsString(),
		 TTimeStamp(endTime).AsString()));

    fRun = run;
    fStartTime = startTime;
    fEndTime = endTime;

}

//______________________________________________________________________________________________
UInt_t AliPMDPreprocessor::Process(TMap* pdaqAliasMap)
{
    
    if(!pdaqAliasMap) return 1;
    TString runType = GetRunType();
    if(runType == "PEDESTAL_RUN"){
	AliPMDPedestal *pedestal = new AliPMDPedestal();
	
        TList* filesources = GetFileSources(kDAQ, "PMD_PED");
	
        if(!filesources) {
	    Log(Form("No sources found for PMD_PED!"));
	    return 1;
	}
	
        AliInfo("Here's the list of sources for PMD_PED");
        filesources->Print();
	
        TIter iter(filesources);
        TObjString* source;
        UInt_t result = 0;
	TString filename;
        while((source=dynamic_cast<TObjString*> (iter.Next()))){
	    filename = GetFile(kDAQ, "PMD_PED", source->GetName());
	    if(filename.Length() == 0) {
		Log(Form("Error retrieving file from source %s failed!", source->GetName()));
		delete filesources;
		return 1;
	    }
	    
	    Log(Form("File with id PMD_PED got from %s", source->GetName()));

	    Int_t det, sm, row, col;
	    Float_t mean, rms;

	    TFile *f= new TFile(filename.Data());
	    if(!f || !f->IsOpen()) 
	    {
		Log(Form("Error opening file with Id PMD_PED from source %s!", source->GetName()));
		return 1;
	    } 
	    TTree *tree = dynamic_cast<TTree *> (f->Get("ped"));
	    if (!tree) 
	    {
		Log("Could not find object \"ped\" in PED file!");
		return 1;
	    }
	    
	    tree->SetBranchAddress("det",  &det);
	    tree->SetBranchAddress("sm",   &sm);
	    tree->SetBranchAddress("row",  &row);
	    tree->SetBranchAddress("col",  &col);
	    tree->SetBranchAddress("mean", &mean);
	    tree->SetBranchAddress("rms",  &rms);


	    Int_t nEntries = (Int_t) tree->GetEntries();
	    for(Int_t i = 0; i < nEntries; i++)
	    {
		tree->GetEntry(i);
		pedestal->SetPedMeanRms(det,sm,row,col,mean,rms);
	    }
	    f->Close();
	    delete f;
	}
	AliCDBMetaData metaData;
	metaData.SetBeamPeriod(0);
	metaData.SetComment("test PMD preprocessor");
	
	result = Store("Calib","Ped", pedestal, &metaData);
	delete pedestal;
	if(result==0)
	{
	    Log("Error storing");                        
	    return 1;
	}
	else
	{
	    return 0;
	}
	
    }else if (runType == "PHYSICS"){
	
	AliPMDCalibData *calibda = new AliPMDCalibData();
	
        TList* filesources = GetFileSources(kDAQ, "PMDGAINS");
	
        if(!filesources) {
	    Log(Form("No sources found for PMDGAINS!"));
	    return 1;
	}
	
        AliInfo("Here's the list of sources for PMDGAINS");
        filesources->Print();
	
        TIter iter(filesources);
        TObjString* source;
        UInt_t result = 0;
	TString filename;
        while((source=dynamic_cast<TObjString*> (iter.Next()))){
	    filename = GetFile(kDAQ, "PMDGAINS", source->GetName());
	    if(filename.Length() == 0) {
		Log(Form("Error retrieving file from source %s failed!", source->GetName()));
		delete filesources;
		return 1;
	    }
	    
	    Log(Form("File with id PMDGAINS got from %s", source->GetName()));

	    Int_t det, sm, row, col;
	    Float_t gain;

	    TFile *f1= new TFile(filename.Data());
	    if(!f1 || !f1->IsOpen()) 
	    {
		Log(Form("Error opening file with Id PMDGAINS from source %s!", source->GetName()));
		return 1;
	    } 
	    TTree *tree = dynamic_cast<TTree *> (f1->Get("ic"));
	    if (!tree) 
	    {
		Log("Could not find object \"ic\" in DAQ file!");
		return 1;
	    }
	    
	    tree->SetBranchAddress("det",  &det);
	    tree->SetBranchAddress("sm",   &sm);
	    tree->SetBranchAddress("row",  &row);
	    tree->SetBranchAddress("col",  &col);
	    tree->SetBranchAddress("gain", &gain);

	    Int_t nEntries = (Int_t) tree->GetEntries();
	    for(Int_t i = 0; i < nEntries; i++)
	    {
		tree->GetEntry(i);

		//if(DET>1 || SM>23 || ROW>95 || COL>95) {
		//    printf("Error! gain[%d,%d,%d,%d] = %f\n",
		//   DET,SM,ROW,COL,GAIN);
		//  continue;
		//    		 	}

		calibda->SetGainFact(det,sm,row,col,gain);
	    }
	    f1->Close();
	    delete f1;
	}

	AliCDBMetaData metaData;
	metaData.SetBeamPeriod(0);
	metaData.SetComment("test PMD preprocessor");
	result = Store("Calib","Gain", calibda, &metaData);
	delete calibda;
	if(result==0)
	{
	    Log("Error storing");                        
	    return 1;
	}
	else
	{
	    return 0;
	}
	
  // Store DCS data for reference
  AliCDBMetaData metadata;
  metadata.SetComment("DCS data for PMD");
  Bool_t resStore = kFALSE;
  resStore = StoreReferenceData("DCS","Data",pdaqAliasMap,&metadata);
	if(resStore==0)
	{
	    Log("Error storing");                        
	    return 1;
	}
	else
	{
	    return 0;
	}

    }
    
    return 2;
}
