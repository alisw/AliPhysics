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
#include <TFile.h>
#include <TSystem.h>

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
  fPyFileName(""),
  fPyFile(""),
  fNumber(1),
  fNfiles(1000),
  fEffMomSmear(kFALSE),
  fMomSmearFilePath(""),
  fMomSmearFileName(""),
  fMomSmearingFile(0x0),
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
  fPyFileName(""),
  fPyFile(""),
  fNumber(1),
  fNfiles(1000),
  fEffMomSmear(kFALSE),
  fMomSmearFilePath(""),
  fMomSmearFileName(""),
  fMomSmearingFile(0x0),
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
  fHistManager.CreateTH2(histname.Data(), title.Data(), 100, -1., 1., 100, -TMath::Pi(), TMath::Pi());
  
  // Particle pT spectrum (after eff and mom smearing)
  histname = "hPyPtSmeared";
  title = histname + ";#it{p}_{T}";
  fHistManager.CreateTH1(histname.Data(), title.Data(), 120, 0, 120);
 
  // Particle pT spectrum origianal
  histname = "hPyPtOrg";
  title = histname + ";#it{p}_{T}";
  fHistManager.CreateTH1(histname.Data(), title.Data(), 120, 0, 120);
  
  // Particle pT spectrum origianal
  histname = "hPyPtEff";
  title = histname + ";#it{p}_{T}";
  fHistManager.CreateTH1(histname.Data(), title.Data(), 120, 0, 120);
 
  // ptsmeared-ptorg of particles
  histname = "hPyPtSmearedPtOrg";
  title = histname + ";#it{p}_{T,smeared};#it{p}_{T,org}";
  fHistManager.CreateTH2(histname.Data(), title.Data(), 50, 0, 50, 50, 0, 50);
 
  
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
 * This function will open new pythia text file
 */
void AliAnalysisTaskPythiaBranchEA::GetNewPythiaFile(){

  if (fPyFilePath.BeginsWith("alien://")) {
     TGrid::Connect("alien://");
 
  }

  for(Int_t k=0;k<100;k++){ //in case that some of the files are missing try other files
     fNumber = fRandom.Integer(fNfiles); //0-999
     fPyFileName = fPyFileMask;
     fPyFileName.ReplaceAll("XXX",Form("%d",fNumber)); 
     fPyFile = Form("%s/%s", fPyFilePath.Data(), fPyFileName.Data());
     TFile::Cp(fPyFile.Data(), fPyFileName.Data());
     fInput = fopen(fPyFileName.Data(),"r");
     
     if(fInput) break;
  } 
 
  if (!fInput) {
    AliError(Form("Could not retrieve PYTHIA file %s!",fPyFileName.Data()));
    return;
  }
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

  GetNewPythiaFile();

  if(fEffMomSmear) ReadMomentumSmearingFile();  
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
 * Read parameters for momentum smearing 
 */
void AliAnalysisTaskPythiaBranchEA::ReadMomentumSmearingFile(){

   if(fMomSmearFilePath.BeginsWith("alien://")) {
      TGrid::Connect("alien://");
   }
   
   TString name = Form("%s/%s", fMomSmearFilePath.Data(), fMomSmearFileName.Data());
   TFile::Cp(name.Data(), fMomSmearFileName.Data());
   fMomSmearingFile = new TFile(fMomSmearFileName.Data(),"READ");
   if(!fMomSmearingFile){
      AliFatal(TString::Format("Momentum smearing file %s not found!", name.Data()));
   }
   
   for(Int_t i=0; i<40; i++){
      name = Form("hPtResolution%d_%d",i,i+1);
      fhSigmaPt[i] = (TH1D*) fMomSmearingFile->Get(name.Data());
   
      if(!fhSigmaPt[i]){
         AliFatal(TString::Format("Momentum smearing histogram %s not found!", name.Data()));
      }
      fhSigmaPt[i]->SetDirectory(0);
   }

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
  Float_t pt, eta, phi, newpt;
  Int_t index=0;

  while(evtLoop){

     if(!fInput){
        GetNewPythiaFile();
     }

     if(fscanf(fInput,"%f %f %f", &pt, &eta, &phi) == EOF){

        fclose(fInput);
        fInput=NULL;

        gSystem->Exec(Form("rm %s", fPyFileName.Data()));        

        continue;
     }

     if(phi > 900){
        evtLoop=0;
        break;
     }else{
        if(fMinEta < eta && eta < fMaxEta){
           newpt = pt;

           fHistManager.FillTH1("hPyPtOrg", pt);

           if(fEffMomSmear){
              if(!PassedDetectorEfficiency(pt)) continue;

              fHistManager.FillTH1("hPyPtEff", pt);
              newpt = SmearPt(pt);

              fHistManager.FillTH2("hPyPtSmearedPtOrg", pt, newpt);
           } 


           new((*fThermalParticlesArray)[index]) AliBasicParticle(eta, phi, newpt, 3);   //use dummy value of charge
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
       
       TString histname = "hPyPtSmeared";
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
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
/**
 *  Track reconstruction efficiency 
 */

Bool_t AliAnalysisTaskPythiaBranchEA::PassedDetectorEfficiency(Float_t pt){

   Double_t eff = 0.;

   if(pt<1.){
      eff = 0.774694;
   }else if(pt < 2.){
      eff = 0.8466;
   }else if(pt < 3.){
      eff = 0.833747;
   }else if(pt < 4.){
      eff = 0.81656;
   }else if(pt < 5.){
      eff = 0.816461; 
   }else if(pt < 7.){
      eff = 0.821438;
   }else if(pt < 10.){
      eff = 0.828123;
   }else if(pt < 100.){
      eff = 0.838103 + (-0.000292742)*pt + (-3.3135e-05)*pt*pt +  (1.41748e-07)*pt*pt*pt; 
   }else{
      eff = 0.61922724; 
   } 

   if(fRandom.Uniform(0.,1.) < eff) return kTRUE; 
   else return kFALSE; 

}
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
/**
 *  Smear Pt 
 */

Float_t AliAnalysisTaskPythiaBranchEA::SmearPt(Float_t pt){

   Double_t sigmaFrac; //sigma_1/pT /  1/pT 
   Int_t index;

   if(pt < 30.){ //read sigma from histogram
      index = TMath::Nint(TMath::Floor(pt));
      if(fhSigmaPt[index]){
         sigmaFrac = fhSigmaPt[index]->GetRandom();
      }else{
         ::Error("AddTaskPythiaBranchEA", "Histogram with momentum smearing is missing.");
         return 0.;
      }
   }else{     // take parameterized value
      sigmaFrac = -0.0291678 + 0.00995427*sqrt( pt + 5.90289); 
   }

   if(pt<1e-5) return 0.;
    
   Double_t oneoverpT    = (Double_t) (1./pt); 
   Double_t sigma        = sigmaFrac*oneoverpT; //sigma_1/pT 
   Double_t oneoverpTnew = fRandom.Gaus(oneoverpT, sigma); 
   Float_t  newpt        = (Float_t) (1./oneoverpTnew);
   if(newpt < 0){
      for(Int_t i=0; i<100; i++){
         oneoverpTnew = fRandom.Gaus(oneoverpT, sigma); 
         newpt        = (Float_t) (1./oneoverpTnew);
         if(newpt>0) break;
      }
   } 
 
   return newpt;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
/**
 *  Terminate clean up
 */
void AliAnalysisTaskPythiaBranchEA::Terminate(Option_t *){

  if(fInput){
        fclose(fInput);
        fInput=NULL;
   }
   gSystem->Exec(Form("rm %s", fPyFileName.Data()));        

   if(fMomSmearingFile){
      fMomSmearingFile->Close();
      delete fMomSmearingFile;
   }
} 

