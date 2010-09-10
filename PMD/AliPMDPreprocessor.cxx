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

/*********************************************
 
 *   PMD Preproccessor Source Code      *
    
    
     0 --> Run Successesfully
     1 --> No pmd Alias is available
 
**********************************************/
 

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

AliPMDPreprocessor::AliPMDPreprocessor(AliShuttleInterface* shuttle) :
  AliPreprocessor("PMD", shuttle)
{
  AddRunType("PHYSICS");
  AddRunType("PEDESTAL");
}

//________________________________________________________ 
AliPMDPreprocessor::~AliPMDPreprocessor()
{
  
}

//________________________________________________________ 
void AliPMDPreprocessor::Initialize(Int_t run, UInt_t startTime, UInt_t endTime)
{
  AliPreprocessor::Initialize(run, startTime, endTime);
  
  AliInfo(Form("\nRun       : %d \nStart Time: %s \nEnd Time  : %s", run,
	       TTimeStamp(startTime).AsString(),
	       TTimeStamp(endTime).AsString()));
  
  fRun = run;
  fStartTime = startTime;
  fEndTime   = endTime;

}

//________________________________________________________ 
Bool_t AliPMDPreprocessor::ProcessDAQ()
{
  TString sRunType = GetRunType();
  Log(Form("RunType %s",sRunType.Data()));
  if (sRunType !="PHYSICS" || sRunType != "PEDESTAL") 
    {
      return kFALSE;
    }
  return kTRUE;
}

//________________________________________________________ 
Bool_t AliPMDPreprocessor::StorePmdPED()
{
  
  AliPMDPedestal *pedestal = new AliPMDPedestal();
  TList* gfsPmdPed = GetFileSources(kDAQ, "PMD_PED.root");
  
  if(!gfsPmdPed) 
    {
      Log(Form("PMDPED: No Shuttle List for PMD PED "));
      return kFALSE;
    }
  else
    {
      AliInfo("PMDPED: Here's the list of sources for PMD_PED.root");
      gfsPmdPed->Print();
      
      TIter iter(gfsPmdPed);
      TObjString* srcPed;
      TString nameOfFile;
      UInt_t resultPmdPed = 0;
      
      while((srcPed = dynamic_cast<TObjString*> (iter.Next())))
	{
	  nameOfFile = GetFile(kDAQ, "PMD_PED.root", srcPed->GetName());
	  if(nameOfFile.Length() == 0) 
	    {
	      Log(Form("PMDPED: Error retrieving file from source %s failed!", srcPed->GetName()));
	      delete gfsPmdPed;
	      return kFALSE;
	    }
	  
	  Log(Form("PMDPED: File with id PMD_PED.root got from %s", srcPed->GetName()));
	  
	  Int_t det, sm, row, col;
	  Float_t mean, rms;
	  
	  TFile *opnFile= new TFile(nameOfFile.Data());
	  if(!opnFile || !opnFile->IsOpen()) 
	    {
	      Log(Form("PMDPED: Error opening file with Id PMD_PED.root from source %s!", srcPed->GetName()));
	      return kFALSE;
	    } 
	  
	  TTree *tree = dynamic_cast<TTree *> (opnFile->Get("ped"));
	  if (!tree) 
	    {
	      Log("PMDPED: Could not find object \"ped\" in PED file!");
	       return kFALSE;
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
	  opnFile->Close();
	  delete opnFile;
	}
      
      AliCDBMetaData mdPED;
      mdPED.SetBeamPeriod(0);
      mdPED.SetResponsible("Satyajit Jena");
      mdPED.SetComment("PMDPED: PMD preprocessor");
      
      resultPmdPed = Store("Calib","Ped", pedestal, &mdPED,0,kTRUE);
      delete pedestal;
      if(resultPmdPed==0)
	{
	  Log("PMDPED: Error storing");                        
	  return kFALSE;
	}
      else
	{
	  return kTRUE;
	}

    } 
 
}

//________________________________________________________
Bool_t AliPMDPreprocessor::StorePmdGAIN()
{
TList* gfsPmdGain = GetFileSources(kDAQ, "PMDGAINS.root");

if(!gfsPmdGain)
    {
      Log(Form("PMDGAIN: No sources found for PMDGAINS.root!"));
      return kFALSE;
    }    
  else
    {
      AliInfo("PMDGAIN: Here's the list of sources for PMDGAINS.root");
      gfsPmdGain->Print();
      
      TIter iter(gfsPmdGain);
      TObjString* source;
      TString nameOfFile;
      UInt_t result = 0;
 
      while((source = dynamic_cast<TObjString*> (iter.Next())))
	{
	  nameOfFile = GetFile(kDAQ, "PMDGAINS.root", source->GetName());
	  if(nameOfFile.Length() == 0) 
	    {
	      Log(Form("PMDGAIN: Error retrieving file from source %s failed!", source->GetName()));
	      delete gfsPmdGain;
	      return kFALSE;
	    }
	  
	  Log(Form("PMDGAIN: File with id PMDGAINS.root got from %s", source->GetName()));
	  
	  Int_t det, sm, row, col;
	  Float_t gain;
	  
	  TFile *opnFile = new TFile(nameOfFile.Data());
	  if (!opnFile || !opnFile->IsOpen()) 
	    {
	      Log(Form("PMDGAIN: Error opening file with Id PMDGAINS.root from source %s!", source->GetName()));
	      return kFALSE;
	    }
	  
	  TTree *tree = dynamic_cast<TTree *> (opnFile->Get("ic"));
	  if (!tree) 
	    {
	      Log("PMDGAIN: Could not find object \"ic\" in DAQ file!");
	      return kFALSE;
	    }
	  
	  tree->SetBranchAddress("det",  &det);
	  tree->SetBranchAddress("sm",   &sm);
	  tree->SetBranchAddress("row",  &row);
	  tree->SetBranchAddress("col",  &col);
	  tree->SetBranchAddress("gain", &gain);
	  
	  Int_t nEntries = (Int_t) tree->GetEntries();
	  AliPMDCalibData *calibda = new AliPMDCalibData();
	  
	  for(Int_t i = 0; i < nEntries; i++)
	    {
	      tree->GetEntry(i);
	      calibda->SetGainFact(det,sm,row,col,gain);
	    }
	  opnFile->Close();
	  delete opnFile;
	  
	  AliCDBMetaData mdGAIN;
	  mdGAIN.SetBeamPeriod(0);
	  mdGAIN.SetResponsible("Satyajit Jena");
	  mdGAIN.SetComment("PMDGAIN: PMD GAIN Data");
	  result = Store("Calib","Gain", calibda, &mdGAIN);
	  delete calibda;
	}
      
      if (result==0)
	{
	  Log("PMDGAIN: Error storing");                        
	  return kFALSE;
	}
      else
	{
	  return kTRUE;
	}
      
    }
}


//___________________________________________________
Bool_t AliPMDPreprocessor::StorePmdHOT()
{
  
  AliPMDHotData *hotda = new AliPMDHotData();
  TList* fsPmdHot = GetFileSources(kDAQ, "PMD_HOT.root");
  
  if(!fsPmdHot) 
    {
      Log(Form("PMDHOT: No sources found for PMD_HOT.root!"));
      return kFALSE;
    }
  else
    {
      AliInfo("PMDHOT: Here's the list of sources for PMD_HOT.root");
      fsPmdHot->Print();
      
      TIter iter(fsPmdHot);
      TObjString* source;
      UInt_t hotresult = 0;
      TString nameOfFile;

      while((source = dynamic_cast<TObjString*> (iter.Next())))
	{
	  nameOfFile = GetFile(kDAQ, "PMD_HOT.root", source->GetName());
	  if(nameOfFile.Length() == 0) 
	    {
	      Log(Form("PMDHOT: Error retrieving file from source %s failed!", source->GetName()));
	      delete fsPmdHot;
	      return kFALSE;
	    }
	  
	  Log(Form("PMDHOT: File with id PMD_HOT.root got from %s", source->GetName()));
	  
	  Int_t det, sm, row, col;
	  Float_t flag;
	  
	  TFile *opnFile = new TFile(nameOfFile.Data());
	  if(!opnFile || !opnFile->IsOpen()) 
	    {
	      Log(Form("PMDHOT: Error opening file with Id PMD_HOT.root from source %s!", source->GetName()));
	      return kFALSE;
	    } 

	  TTree *tree = dynamic_cast<TTree *> (opnFile->Get("hot"));
	  if (!tree) 
	    {
	      Log("PMDHOT: Could not find object \"hot\" in DAQ file!");
	      return kFALSE;
	    }
	  
	  tree->SetBranchAddress("det",  &det);
	  tree->SetBranchAddress("sm",   &sm);
	  tree->SetBranchAddress("row",  &row);
	  tree->SetBranchAddress("col",  &col);
	  tree->SetBranchAddress("flag", &flag);
	  
	  Int_t nEntries = (Int_t) tree->GetEntries();
	  for(Int_t j = 0; j < nEntries; j++)
	    {
	      tree->GetEntry(j);
	      hotda->SetHotChannel(det,sm,row,col,flag);
	    }
	  opnFile->Close();
	  delete opnFile;
	}
      
      AliCDBMetaData metaData;
      metaData.SetBeamPeriod(0);
      metaData.SetResponsible("Satyajit Jena");
      metaData.SetComment("PMDHOT: PMD preprocessor");
      hotresult = Store("Calib","Hot", hotda, &metaData);
      delete hotda;
      if(hotresult==0)
	{
	  Log("PMDHOT: Error storing");                        
	  return kFALSE;
	}
      else
	{
	  return kTRUE;
	}
      
    }
}

//________________________________________________________
Bool_t AliPMDPreprocessor::StorePmdMEAN()
{
  AliPMDMeanSm *smmeanda = new AliPMDMeanSm();
  TList* gfsPmdMean = GetFileSources(kDAQ, "PMD_MEAN_SM.root");
  
  if(!gfsPmdMean) 
    {
      Log(Form("PMDMEAN: No sources found for PMD_MEAN_SM.root!"));
      return kFALSE;
    }
  else
    {
      AliInfo("PMDMEAN: Here's the list of sources for PMD_MEAN_SM.root");
      gfsPmdMean->Print();
      
      TIter iter(gfsPmdMean);
      TObjString* sourc;
      UInt_t storeMeanData = 0;
      TString filenam;
      
      while((sourc=dynamic_cast<TObjString*> (iter.Next())))
	{
	  filenam = GetFile(kDAQ, "PMD_MEAN_SM.root", sourc->GetName());
	  if(filenam.Length() == 0) 
	    {
	      Log(Form("PMDMEAN: Error retrieving file from source %s failed!", sourc->GetName()));
	      delete gfsPmdMean;
	      return kFALSE;
	    }
	  
	  Log(Form("PMDMEAN: File with id PMD_MEAN_SM.root got from %s", sourc->GetName()));
	  
	  Int_t det = 0, sm = 0;
	  Float_t smmean = 0.;
	  
	  TFile *opnFile = new TFile(filenam.Data());
	  if(!opnFile || !opnFile->IsOpen()) 
	    {
	      Log(Form("PMDMEAN: Error opening file with Id PMD_MEAN_SM.root from source %s!", sourc->GetName()));
	      return kFALSE;
	    }
	  
	  TTree *tree = dynamic_cast<TTree *> (opnFile->Get("mean"));
	  if (!tree) 
	    {
	      Log("PMDMEAN: Could not find object \"hot\" in DAQ file!");
	      return kFALSE;
	    }
	  
	  tree->SetBranchAddress("det",  &det);
	  tree->SetBranchAddress("sm",   &sm);
	  tree->SetBranchAddress("smmean", &smmean);
	  
	  Int_t nEntries = (Int_t) tree->GetEntries();
	  for(Int_t j = 0; j < nEntries; j++)
	    {
	      tree->GetEntry(j);
	      smmeanda->SetMeanSm(det,sm,smmean);
	    }
	  opnFile->Close();
	  delete opnFile;
	}
      
      AliCDBMetaData mdMEAN;
      mdMEAN.SetBeamPeriod(0);
      mdMEAN.SetResponsible("Satyajit Jena");
      mdMEAN.SetComment("PMDMEAN: PMD preprocessor");
      
      storeMeanData = Store("Calib","SMMEAN", smmeanda, &mdMEAN);
      delete smmeanda;

      if(storeMeanData==0)
	{
	  Log("PMDMEAN: Error storing");                        
	  return kFALSE;
	}
      else
	{
	  return kTRUE;
	}
    }
  
}

//_____________________________________________________________
Bool_t AliPMDPreprocessor::StorePmdDCS(TMap *sDaqAM)
{
	
  AliCDBMetaData mdDCS;
  mdDCS.SetResponsible("Satyajit Jena");
  mdDCS.SetComment("DCS data for PMD");

  Bool_t resStore = kFALSE;
  resStore = StoreReferenceData("DCS","Data",sDaqAM,&mdDCS);

  if(resStore == 0)
    {
      Log("PMDDP: Error storing");                        
      return kFALSE;
    }
  else
    {
      return kTRUE;
    }

}

//_____________________________________________________________________

UInt_t AliPMDPreprocessor::Process(TMap* pmdDaqAliasMap)
{
  
  if(!pmdDaqAliasMap)
    {
      return 1;
    }
  
  TString runType = GetRunType();

  if(runType == "PEDESTAL")
    {
 
      Log(Form("------------------ PMD Pedestal --------------"));
      Bool_t pmdPed = StorePmdPED();
      if(!pmdPed)
	{
	  Log(Form("ERROR:  Couldn't write PMD pedestal file"));
	  return 1;
	}
      else
	{
	  Log(Form("Storing of PMD Pedestal File is Successful"));
	  return 0;
	}
    } 
  
  else if (runType == "PHYSICS")
      {
	Log(Form("------------------- PMD GAIN----------------"));	
	Bool_t pmdGAIN = StorePmdGAIN();
	if (!pmdGAIN)
	  {
	    Log(Form("ERROR:  Couldn't write PMD GAIN file"));
	   
	  }
	else
	  {
	    Log(Form("Storing of PMD GAIN File is Successful"));
	   
	  }
	
	Log(Form("------------------- PMD HOT ----------------"));
	Bool_t pmdHOT = StorePmdHOT();
	if (!pmdHOT)
	  {
	    Log(Form("ERROR:  Couldn't write PMD HOT file"));
	   
	  }
	else
	  {
	    Log(Form("Storing of PMD HOT File is Successful"));
	   
	  }
	
	Log(Form("------------------- SM MEAN ----------------"));
	Bool_t pmdMEAN = StorePmdMEAN();
	if (!pmdMEAN)
	  {
	    Log(Form("ERROR:  Couldn't write PMD SM MEAN file"));
	   
	  }
	else
	  {
	    Log(Form("Storing of PMD SM MEAN File is Successful"));
	   
	  }
	
	Log(Form("------------------- DCS DP -----------------"));
	Bool_t pmdDCS = StorePmdDCS(pmdDaqAliasMap);
	if (!pmdDCS)
	  {
	    Log(Form("ERROR:  Couldn't write PMD DCS DP file"));
	   
	  }
	else
	  {
	    Log(Form("Storing of PMD DCS dp is Successul"));
	   
	  }
	
	if (pmdGAIN || pmdHOT || pmdMEAN || pmdDCS) 
	  return 0;
	else 
	  return 1;
	
      }
  return 2;
}


