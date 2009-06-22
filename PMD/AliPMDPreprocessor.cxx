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
#include "AliPMDHotData.h"
#include "AliPMDMeanSm.h"
#include "AliPMDPedestal.h"
#include "AliPMDPreprocessor.h"


ClassImp(AliPMDPreprocessor)
  
//______________________________________________________________________________________________
AliPMDPreprocessor::AliPMDPreprocessor(AliShuttleInterface* shuttle) :
  AliPreprocessor("PMD", shuttle)
{
  // constructor
  AddRunType("PHYSICS");
  AddRunType("PEDESTAL");
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

//-----------------------------------------
Bool_t AliPMDPreprocessor::ProcessDAQ()
{
  TString RunType = GetRunType();
  Log(Form("RunType %s",RunType.Data()));
  if (RunType !="PHYSICS" || RunType != "PEDESTAL") {
    return kFALSE;
  }

  return kTRUE;
}


//_____________________________________________________________________
UInt_t AliPMDPreprocessor::Process(TMap* pdaqAliasMap)
{
  
  if(!pdaqAliasMap) return 1;
  TString runType = GetRunType();
  if(runType == "PEDESTAL"){
    AliPMDPedestal *pedestal = new AliPMDPedestal();
    
    TList* filesources = GetFileSources(kDAQ, "PMD_PED.root");
    
    if(!filesources) {
      Log(Form("No sources found for PMD_PED.root!"));
      return 1;
    }
    
    AliInfo("Here's the list of sources for PMD_PED.root");
    filesources->Print();
    
    TIter iter(filesources);
    TObjString* source;
    UInt_t result = 0;
    TString filename;
    while((source=dynamic_cast<TObjString*> (iter.Next()))){
      filename = GetFile(kDAQ, "PMD_PED.root", source->GetName());
      if(filename.Length() == 0) {
	Log(Form("Error retrieving file from source %s failed!", source->GetName()));
	delete filesources;
	return 1;
      }
      
      Log(Form("File with id PMD_PED.root got from %s", source->GetName()));
      
      Int_t det, sm, row, col;
      Float_t mean, rms;
      
      TFile *f= new TFile(filename.Data());
      if(!f || !f->IsOpen()) 
	{
	  Log(Form("Error opening file with Id PMD_PED.root from source %s!", source->GetName()));
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
    
    result = Store("Calib","Ped", pedestal, &metaData,0,kTRUE);
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
    
    TList* filesources = GetFileSources(kDAQ, "PMDGAINS.root");
    
    if(!filesources) {
      Log(Form("No sources found for PMDGAINS.root!"));
      return 1;
    }
    
    AliInfo("Here's the list of sources for PMDGAINS.root");
    filesources->Print();
    
    TIter iter(filesources);
    TObjString* source;
    UInt_t result = 0;
    TString filename;
    while((source=dynamic_cast<TObjString*> (iter.Next()))){
      filename = GetFile(kDAQ, "PMDGAINS.root", source->GetName());
      if(filename.Length() == 0) {
	Log(Form("Error retrieving file from source %s failed!", source->GetName()));
	delete filesources;
	return 1;
      }
      
      Log(Form("File with id PMDGAINS.root got from %s", source->GetName()));
      
      Int_t det, sm, row, col;
      Float_t gain;

      TFile *f1= new TFile(filename.Data());
      if(!f1 || !f1->IsOpen()) 
	{
	  Log(Form("Error opening file with Id PMDGAINS.root from source %s!", source->GetName()));
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
	  //  printf("Error! gain[%d,%d,%d,%d] = %f\n",
	  //   DET,SM,ROW,COL,GAIN);
	  //  continue;
	  //}
	  
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
    //------------------For Storing HOT Data-------------------------//
    AliPMDHotData *hotda = new AliPMDHotData();
    TList* filesource = GetFileSources(kDAQ, "PMD_HOT.root");
	
    if(!filesource) {
      Log(Form("No sources found for PMD_HOT.root!"));
      return 1;
    }
    
    AliInfo("Here's the list of sources for PMD_HOT.root");
    filesource->Print();
    
    TIter iter2(filesource);
    TObjString* sources;
    UInt_t hotresult = 0;
    TString filenames;
    while((sources=dynamic_cast<TObjString*> (iter2.Next()))){
      filenames = GetFile(kDAQ, "PMD_HOT.root", sources->GetName());
      if(filenames.Length() == 0) {
	Log(Form("Error retrieving file from source %s failed!", sources->GetName()));
	delete filesource;
	return 1;
      }
      
      Log(Form("File with id PMD_HOT.root got from %s", sources->GetName()));
      
      Int_t det, sm, row, col;
      Float_t flag;
      
      TFile *f2= new TFile(filenames.Data());
      if(!f2 || !f2->IsOpen()) 
	{
	  Log(Form("Error opening file with Id PMD_HOT.root from source %s!", sources->GetName()));
	  return 1;
	} 
      TTree *tree1 = dynamic_cast<TTree *> (f2->Get("hot"));
      if (!tree1) 
	{
	  Log("Could not find object \"hot\" in DAQ file!");
	  return 1;
	}
      
      tree1->SetBranchAddress("det",  &det);
      tree1->SetBranchAddress("sm",   &sm);
      tree1->SetBranchAddress("row",  &row);
      tree1->SetBranchAddress("col",  &col);
      tree1->SetBranchAddress("flag", &flag);
      
      
      
      Int_t nEntries = (Int_t) tree1->GetEntries();
      for(Int_t j = 0; j < nEntries; j++)
	{
	  tree1->GetEntry(j);
	  
	  //if(det>1 || sm>23 || row>95 || col>95) {
	  // printf("Error! gain[%d,%d,%d,%d] = %f\n",
	  // det,sm,row,col,flag);
	  //continue;
	  //}
	  
	  hotda->SetHotChannel(det,sm,row,col,flag);
	}
      f2->Close();
      delete f2;
    }
    hotresult = Store("Calib","Hot", hotda, &metaData);
    delete hotda;
    if(hotresult==0)
      {
	Log("Error storing");                        
	return 1;
      }
    //-------------------------------------------------------------------  
//-----------------------------------for storing SM MEAN--------------//
	
    //------------------For Storing HOT Data-------------------------//
    AliPMDMeanSm *smmeanda = new AliPMDMeanSm();
    TList* filesourc = GetFileSources(kDAQ, "PMD_MEAN_SM.root");
	
    if(!filesourc) {
      Log(Form("No sources found for PMD_MEAN_SM.root!"));
      return 1;
    }
    
    AliInfo("Here's the list of sources for PMD_MEAN_SM.root");
    filesourc->Print();
    
    TIter iter3(filesourc);
    TObjString* sourc;
    UInt_t meanresult = 0;
    TString filenam;
    while((sourc=dynamic_cast<TObjString*> (iter3.Next()))){
      filenam = GetFile(kDAQ, "PMD_MEAN_SM.root", sourc->GetName());
      if(filenam.Length() == 0) {
	Log(Form("Error retrieving file from source %s failed!", sourc->GetName()));
	delete filesourc;
	return 1;
      }
      
      Log(Form("File with id PMD_MEAN_SM.root got from %s", sourc->GetName()));
      
      Int_t det, sm ;
      Float_t smmean;
      
      TFile *f3= new TFile(filenam.Data());
      if(!f3 || !f3->IsOpen()) 
	{
	  Log(Form("Error opening file with Id PMD_MEAN_SM.root from source %s!", sourc->GetName()));
	  return 1;
	} 
      TTree *tree2 = dynamic_cast<TTree *> (f3->Get("mean"));
      if (!tree2) 
	{
	  Log("Could not find object \"hot\" in DAQ file!");
	  return 1;
	}
      
      tree2->SetBranchAddress("det",  &det);
      tree2->SetBranchAddress("sm",   &sm);
      tree2->SetBranchAddress("smmean", &smmean);
      
      Int_t nEntries = (Int_t) tree2->GetEntries();
      for(Int_t j = 0; j < nEntries; j++)
	{
	  tree2->GetEntry(j);
	  smmeanda->SetMeanSm(det,sm,smmean);
	}
      f3->Close();
      delete f3;
    }
    meanresult = Store("Calib","SMMEAN", smmeanda, &metaData);
    delete smmeanda;
    if(meanresult==0)
      {
	Log("Error storing");                        
	return 1;
      }
    // -------------Store DCS data for reference------------//
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


