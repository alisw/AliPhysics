/**********************************************************************************
* Copyright (C) 2016, Copyright Holders of the ALICE Collaboration                *
* All rights reserved.                                                            *
*                                                                                 *
* Redistribution and use in source and binary forms, with or without              *
* modification, are permitted provided that the following conditions are met:     *
*   * Redistributions of source code must retain the above copyright              *
*     notice, this list of conditions and the following disclaimer.               *
*   * Redistributions in binary form must reproduce the above copyright           *
*     notice, this list of conditions and the following disclaimer in the         *
*     documentation and/or other materials provided with the distribution.        *
*   * Neither the name of the <organization> nor the                              *
*     names of its contributors may be used to endorse or promote products        *
*     derived from this software without specific prior written permission.       *
*                                                                                 *
* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND *
* ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED   *
* WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE          *
* DISCLAIMED. IN NO EVENT SHALL ALICE COLLABORATION BE LIABLE FOR ANY             *
* DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES      *
* (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;    *
* LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND     *
* ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT      *
* (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS   *
* SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.                    *
* *********************************************************************************/

#include "AliAnalysisTaskPythiaBranchEA.h"

#include <TClonesArray.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TF1.h>
#include <TGrid.h>

#include <AliAnalysisManager.h>
#include <AliVEventHandler.h>
#include <AliTLorentzVector.h>
#include <AliBasicParticle.h>

/// \cond CLASSIMP
ClassImp(PWGJE::EMCALJetTasks::AliAnalysisTaskPythiaBranchEA);

using namespace PWGJE::EMCALJetTasks;
using namespace std;

/// \endcond

/**
 * Default constructor. Needed by ROOT I/O
 */
AliAnalysisTaskPythiaBranchEA::AliAnalysisTaskPythiaBranchEA() :
  AliAnalysisTaskSE(),
  fEventInitialized(false),
  fEvent(nullptr),
  fRandom(0),
  fOutputCollectionName(""),
  fThermalParticlesArray(),
  fMinEta(-0.9),
  fMaxEta(0.9),
  fPyFilePath(""),
  fPyFileMask(""),
  fPyFile(""),
  fNumber(1),
  fNfiles(1000),
  fOutput(0),
  fHistManager(),
  fInput(0x0)
{}

/**
 * Standard constructor. Should be used by the user.
 * @param[in] name Name of the task
 */
AliAnalysisTaskPythiaBranchEA::AliAnalysisTaskPythiaBranchEA(const char* name) :
  AliAnalysisTaskSE(name),
  fEventInitialized(false),
  fEvent(nullptr),
  fRandom(0),
  fOutputCollectionName(""),
  fThermalParticlesArray(),
  fMinEta(-0.9),
  fMaxEta(0.9),
  fPyFilePath(""),
  fPyFileMask(""),
  fPyFile(""),
  fNumber(1),
  fNfiles(1000),
  fOutput(0),
  fHistManager(name),
  fInput(0x0)
{
  DefineOutput(1, TList::Class());
}

/**
 * Destructor
 */
AliAnalysisTaskPythiaBranchEA::~AliAnalysisTaskPythiaBranchEA() {}

/**
 * Performing run-independent initialization. Here the histograms should be instantiated.
 */
void AliAnalysisTaskPythiaBranchEA::UserCreateOutputObjects()
{
  //
  TString title, histname; 
  /////////////////////////////////////
  // Allocate histograms
  
  // Eta-phi of particles
  histname = "hPyEtaPhi";
  title = histname + ";#eta;#phi";
  fHistManager.CreateTH2(histname.Data(), title.Data(), 100, -1., 1., 100, 0, 2*TMath::Pi());
  
  // Particle pT spectrum
  histname = "hPyPt";
  title = histname + ";#it{p}_{T}";
  fHistManager.CreateTH1(histname.Data(), title.Data(), 200, 0, 20);
  
  
  
  /////////////////////////////////////
  
  // Allow for output files
  OpenFile(1);
  fOutput = new TList();
  fOutput->SetOwner();
  
  TIter next(fHistManager.GetListOfHistograms());
  TObject* obj = 0;
  while ((obj = next())) {
    fOutput->Add(obj);
  }
  
  PostData(1, fOutput); // Post data for ALL output slots > 0 here.
}

/**
 * This function is executed automatically for the first event. Some extra initialization can be performed here.
 */
void AliAnalysisTaskPythiaBranchEA::ExecOnce()
{
  fEvent = InputEvent();
  if (!fEvent) {
    AliError("Could not retrieve event! Returning!");
    return;
  }

  fRandom.SetSeed(0);
  gRandom = &fRandom; // It seems this is necessary in order to use the TF1::GetRandom() function...


  if (fPyFile.BeginsWith("alien://")) {
    TGrid::Connect("alien://");
  }

  fNumber = fRandom.Integer(fNfiles); //0-999
  fPyFile = fPyFileMask;
  fPyFile.ReplaceAll("XXX",Form("%d",fNumber)); 
  fPyFile = Form("%s/%s", fPyFilePath.Data(), fPyFile.Data());
  fInput = fopen(fPyFile.Data(),"r");
 
  if (!fInput) {
    AliError(Form("Could not retrieve PYTHIA file %s!",fPyFile.Data()));
    return;
  }
 
  // Create a new collection during the first event
  CreateNewObjectBranch();
  
  fEventInitialized = true;
}

/**
 * Steers creation of a new collection in the event. Adapted from AliEmcalCopyCollection.
 */
void AliAnalysisTaskPythiaBranchEA::CreateNewObjectBranch()
{
  // Check to ensure that we are not trying to create a new branch on top of the old branch
  TObject* existingObject = fEvent->FindListObject(fOutputCollectionName.c_str());
  if (existingObject) {
    AliFatal(TString::Format("Attempted to create a new branch, \"%s\", with the same name as an existing branch!", fOutputCollectionName.c_str()));
  }
  
  // Create new branch and add it to the event. It will then be there for every event, and cleared for each event.
  fThermalParticlesArray = new TClonesArray("AliBasicParticle", 1000);
  fThermalParticlesArray->SetName(fOutputCollectionName.c_str());
  fEvent->AddObject(fThermalParticlesArray);
}

/**
 * Steers each event. It enforces that the event is initialized before executing the main analysis of the event.
 */
void AliAnalysisTaskPythiaBranchEA::UserExec(Option_t *option)
{
  // Initialize the event if not initialized
  if (!fEventInitialized) {
    ExecOnce();
  }
  
  // Only continue if we are initialized successfully
  if (!fEventInitialized) {
    return;
  }
  
  // Generate the thermal particles, and write them to the event
  Run();
  
  // Plot some properties of the generated thermal particles
  FillHistograms();
  
}

/**
 * Run analysis code here.
 */
void AliAnalysisTaskPythiaBranchEA::Run()
{
  // Get the TClonesArray
  fThermalParticlesArray = dynamic_cast<TClonesArray*>(fEvent->FindListObject(fOutputCollectionName.c_str()));
  if (!fThermalParticlesArray) {
    AliError(Form("%s: Could not retrieve array with name %s!", GetName(), fOutputCollectionName.c_str()));
    return;
  }

  //read event from text file 
  Bool_t evtLoop = 1;
  Float_t pt, eta, phi;
  Int_t index=0;

  while(evtLoop){

     if(!fInput){
        if (fPyFilePath.BeginsWith("alien://")) {
          TGrid::Connect("alien://");
        }

        fNumber = fRandom.Integer(fNfiles); //0-999
        fPyFile = fPyFileMask;
        fPyFile.ReplaceAll("XXX",Form("%d",fNumber)); 
        fPyFile = Form("%s/%s", fPyFilePath.Data(), fPyFile.Data());

        fInput = fopen(fPyFile.Data(),"r");

        if (!fInput) {
          AliError(Form("Could not retrieve PYTHIA file %s!",fPyFile.Data()));
          return;
        }
 
     }

     if(fscanf(fInput,"%f %f %f", &pt, &eta, &phi) == EOF){

        fclose(fInput);
        fInput=NULL;

        continue;
     }

     if(phi > 900){
        evtLoop=0;
        break;
     }else{
        if(fMinEta < eta && eta < fMaxEta){
           new((*fThermalParticlesArray)[index]) AliBasicParticle(eta, phi, pt, 3);   //use dummy value of charge
           index++;
        }
     }
  }

  fThermalParticlesArray->Compress();
}

/**
 * Loop over particles to fill histograms.
 */
void AliAnalysisTaskPythiaBranchEA::FillHistograms()
{
   // Loop over particles.
   Int_t N = fThermalParticlesArray->GetEntries();
   for(Int_t i=0; i<N; i++) {
       AliBasicParticle* thermalParticle = dynamic_cast<AliBasicParticle*>(fThermalParticlesArray->At(i));
       Double_t pT = thermalParticle->Pt();
       Double_t eta = thermalParticle->Eta();
       Double_t phi = thermalParticle->Phi();
       
       TString histname = "hPyPt";
       fHistManager.FillTH1(histname, pT);
       
       histname = "hPyEtaPhi";
       fHistManager.FillTH2(histname, eta, phi);
   }
}



//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
/**
 * AddTask.
 */
AliAnalysisTaskPythiaBranchEA* AliAnalysisTaskPythiaBranchEA::AddTaskPythiaBranchEA(
  const char *outputCollectionName,
  const char *pyfilepath,
  const char *pyfilemask,
  const char *suffix)
{
  // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr)
  {
    ::Error("AddTaskPythiaBranchEA", "No analysis manager to connect to.");
    return 0;
  }
  
  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  AliVEventHandler* handler = mgr->GetInputEventHandler();
  if (!handler)
  {
    ::Error("AddTaskPythiaBranchEA", "This task requires an input event handler");
    return 0;
  }
  
  enum EDataType_t {
    kUnknown,
    kESD,
    kAOD
  };
  
  EDataType_t dataType = kUnknown;
  
  if (handler->InheritsFrom("AliESDInputHandler")) {
    dataType = kESD;
  }
  else if (handler->InheritsFrom("AliAODInputHandler")) {
    dataType = kAOD;
  }
  
  //-------------------------------------------------------
  // Init the task and do settings
  //-------------------------------------------------------
  TString name("AliAnalysisTaskPythiaBranchEA");
  if (strcmp(outputCollectionName,"") != 0) {
    name += "_";
    name += outputCollectionName;
  }
  if (strcmp(suffix,"") != 0) {
    name += "_";
    name += suffix;
  }

  /////////////////////////////////////////////////////////////
  // Configure the task
  AliAnalysisTaskPythiaBranchEA* task = new AliAnalysisTaskPythiaBranchEA(name.Data());

  task->SetOutputCollectionName(outputCollectionName);
  task->SetPythiaFilePath(pyfilepath);
  task->SetPythiaFileMask(pyfilemask);
  
  //-------------------------------------------------------
  // Final settings, pass to manager and set the containers
  //-------------------------------------------------------
  
  mgr->AddTask(task);
  
  // Create containers for input/output
  AliAnalysisDataContainer *cinput1  = mgr->GetCommonInputContainer()  ;
  TString contname(name);
  contname += "_histos";
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(contname.Data(),
                                                            TList::Class(),AliAnalysisManager::kOutputContainer,
                                                            Form("%s", AliAnalysisManager::GetCommonFileName()));
  mgr->ConnectInput  (task, 0,  cinput1 );
  mgr->ConnectOutput (task, 1, coutput1 );
  
  return task;
}
