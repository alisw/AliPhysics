/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE HLT Project.                                    *
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
// Author: Mikolaj Krzewicki, mikolaj.krzewicki@cern.ch, 2014-12-15

#include "TFile.h"
#include "TDirectory.h"
#include "TList.h"
#include "TCollection.h"
#include "TH1.h"
#include "TTree.h"
#include "AliAnalysisTask.h"
#include "AliAnalysisDataContainer.h"
#include "AliHLTAnalysisManager.h"

//_____________________________________________________________________________
AliHLTAnalysisManager::AliHLTAnalysisManager(const char* name, const char* title)
  :
  AliAnalysisManager(name,title)
{
  //ctor
  //on HLT by default we have an external event loop
  SetExternalLoop(kTRUE);
}

//_____________________________________________________________________________
AliHLTAnalysisManager::~AliHLTAnalysisManager()
{
  //dtor
}

//_____________________________________________________________________________
Bool_t AliHLTAnalysisManager::InitAnalysis()
{
  //initialize everything
  if (AliAnalysisManager::InitAnalysis() == kFALSE)
  {
    return kFALSE;
  }
  Bool_t dirStatus = TH1::AddDirectoryStatus();  
  TIter nextT(GetTasks());
  AliAnalysisTask* task=NULL;
  //inititialize the tasks
  while ((task=(AliAnalysisTask*)nextT())) 
  {
    TH1::AddDirectory(kFALSE);
    task->CreateOutputObjects();
  }
  TH1::AddDirectory(dirStatus);
  return kTRUE;
}

//_____________________________________________________________________________
Int_t AliHLTAnalysisManager::WriteAnalysisToFile()
{
  //write the data to file
  TDirectory *opwd = gDirectory;  
  TIter nextOutput(GetOutputs());
  TList listOfOpenFiles;
  while (AliAnalysisDataContainer* output=(AliAnalysisDataContainer*)nextOutput())
  {
    const char* filename   = output->GetFileName();
    const char* openoption = "RECREATE";
    TFile* file=NULL;
    if (!(file=(TFile*)listOfOpenFiles.FindObject(filename)))
      { file = new TFile(filename,openoption); }
    output->SetFile(file);
    listOfOpenFiles.Add(file);
    file->cd();
    TString dir = output->GetFolderName();
    if (!dir.IsNull())
    {
      if (!file->GetDirectory(dir)) file->mkdir(dir);
      file->cd(dir);
    }
    TObject* outputData=output->GetData();
    if (!outputData) 
    {
      continue;
    }
    outputData->Print();
    if (outputData->InheritsFrom(TCollection::Class())) 
    {
      // If data is a collection, we set the name of the collection 
      // as the one of the container and we save as a single key.
      TCollection *coll = (TCollection*)output->GetData();
      coll->SetName(output->GetName());
      coll->Write(output->GetName(), TObject::kSingleKey);
    } 
    else 
    {
      if (outputData->InheritsFrom(TTree::Class())) 
      {
        TTree *tree = (TTree*)output->GetData();
        tree->SetDirectory(gDirectory);
        tree->AutoSave();
      } 
      else 
      {
        output->GetData()->Write();
      }   
    }
    opwd->cd();
  }
  return 0;
}

//_____________________________________________________________________________
Int_t AliHLTAnalysisManager::ResetOutputData()
{
  //call tasks ResetOutputData() methods to clear the stats in the output data
  //after sending them out
  TIter nextOutput(GetTasks());
  Int_t n=0;
  while (AliAnalysisTask* output=(AliAnalysisTask*)nextOutput())
  {
    output->ResetOutputData();
    n++;
  }
  return n;
}
