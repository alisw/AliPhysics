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

#include "AliAnalysisTaskGenerateThermalBackgroundMC.h"

#include <TClonesArray.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TF1.h>

#include <AliAnalysisManager.h>
#include <AliVEventHandler.h>
#include <AliTLorentzVector.h>
#include <AliBasicParticle.h>

/// \cond CLASSIMP
ClassImp(AliAnalysisTaskGenerateThermalBackgroundMC);
/// \endcond

/**
 * Default constructor. Needed by ROOT I/O
 */
AliAnalysisTaskGenerateThermalBackgroundMC::AliAnalysisTaskGenerateThermalBackgroundMC() :
  AliAnalysisTaskSE(),
  fEventInitialized(false),
  fEvent(nullptr),
  fOutputCollectionName(""),
  fThermalParticlesArray(),
  fChargedParticleFraction(0.67),
  fRandom(0),
  fPtDistribution(nullptr),
  fAlpha(2),
  fBeta(0.3),
  fMinChargedPt(0.15),
  fMinNeutralPt(0.3),
  fMaxPt(100.),
  fUseGaussianForN(true),
  fNGaussianMean(3500),
  fNGaussianSigma(500),
  fMinN(3000),
  fMaxN(4500),
  fMinEta(-0.9),
  fMaxEta(0.9),
  fOutput(0),
  fHistManager()
{}

/**
 * Standard constructor. Should be used by the user.
 * @param[in] name Name of the task
 */
AliAnalysisTaskGenerateThermalBackgroundMC::AliAnalysisTaskGenerateThermalBackgroundMC(const char* name) :
  AliAnalysisTaskSE(name),
  fEventInitialized(false),
  fEvent(nullptr),
  fOutputCollectionName(""),
  fThermalParticlesArray(),
  fChargedParticleFraction(0.67),
  fRandom(0),
  fPtDistribution(nullptr),
  fAlpha(2),
  fBeta(0.3),
  fMinChargedPt(0.15),
  fMinNeutralPt(0.3),
  fMaxPt(100.),
  fUseGaussianForN(true),
  fNGaussianMean(3500),
  fNGaussianSigma(500),
  fMinN(3000),
  fMaxN(4500),
  fMinEta(-0.9),
  fMaxEta(0.9),
  fOutput(0),
  fHistManager(name)
{
  DefineOutput(1, TList::Class());
}

/**
 * Destructor
 */
AliAnalysisTaskGenerateThermalBackgroundMC::~AliAnalysisTaskGenerateThermalBackgroundMC() {}

/**
 * Performing run-independent initialization. Here the histograms should be instantiated.
 */
void AliAnalysisTaskGenerateThermalBackgroundMC::UserCreateOutputObjects()
{
  // Define gamma distribution f_gamma to sample single particle pT.
  // We use the convention in TMath::GammaDist(x, alpha, mu, beta). We fix mu = 0.
  // Then f_gamma(x;a,b) = x^(a-1) * (1/b)^a * e^(-x/b) ~  x^(a-1) * e^(-x/b).
  fPtDistribution = new TF1("fGamma", "TMath::GammaDist(x, [0], 0, [1])", TMath::Min(fMinChargedPt, fMinNeutralPt), fMaxPt);
  fPtDistribution->SetParameter(0, fAlpha);
  fPtDistribution->SetParameter(1, fBeta);
  
  /////////////////////////////////////
  // Allocate histograms
  
  // Number of particles per event vs. mean pT per event
  TString histname = "hNvsMeanPt";
  TString title = histname + ";#it{N}_{particles};<#it{p}_{T}>";
  if (fUseGaussianForN) {
    fHistManager.CreateTH2(histname.Data(), title.Data(), 200, fNGaussianMean-5*fNGaussianSigma, fNGaussianMean+5*fNGaussianSigma, 200, 0., 4*fBeta);
  }
  else {
    fHistManager.CreateTH2(histname.Data(), title.Data(), 200, fMinN*0.9, fMaxN*1.1, 200, 0., 4*fBeta);
  }
    
  // Eta-phi of particles
  histname = "hEtaPhi";
  title = histname + ";#eta;#phi";
  fHistManager.CreateTH2(histname.Data(), title.Data(), 100, -1., 1., 100, 0, 2*TMath::Pi());
  
  // Particle pT spectrum
  histname = "hPt";
  title = histname + ";#it{p}_{T}";
  fHistManager.CreateTH1(histname.Data(), title.Data(), 200, 0, 20);
  
  // Plot highest-pT particle in the event
  histname = "hLeadingPt";
  title = histname + ";#it{p}_{T,leading}";
  fHistManager.CreateTH1(histname.Data(), title.Data(), 200, 0, 20);
  
  // Background density
  histname = "hBackgroundDensity";
  title = histname + ";#rho (GeV/#it{c});counts";
  fHistManager.CreateTH1(histname.Data(), title.Data(), 200, 0, 500);
  
  // DeltaPt
  histname = "hDeltaPt02";
  title = histname + ";#delta#it{p}_{T} (GeV/#it{c});counts";
  fHistManager.CreateTH1(histname.Data(), title.Data(), 400, -100, 100);
  
  histname = "hDeltaPt04";
  title = histname + ";#delta#it{p}_{T} (GeV/#it{c});counts";
  fHistManager.CreateTH1(histname.Data(), title.Data(), 400, -100, 100);
  
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
void AliAnalysisTaskGenerateThermalBackgroundMC::ExecOnce()
{
  fEvent = InputEvent();
  if (!fEvent) {
    AliError("Could not retrieve event! Returning!");
    return;
  }
  
  // Initialize random number generator to random seed
  fRandom.SetSeed(0);
  gRandom = &fRandom; // It seems this is necessary in order to use the TF1::GetRandom() function...
  
  // Create a new collection during the first event
  CreateNewObjectBranch();
  
  fEventInitialized = true;
}

/**
 * Steers creation of a new collection in the event. Adapted from AliEmcalCopyCollection.
 */
void AliAnalysisTaskGenerateThermalBackgroundMC::CreateNewObjectBranch()
{
  // Check to ensure that we are not trying to create a new branch on top of the old branch
  TObject* existingObject = fEvent->FindListObject(fOutputCollectionName.c_str());
  if (existingObject) {
    AliFatal(TString::Format("Attempted to create a new branch, \"%s\", with the same name as an existing branch!", fOutputCollectionName.c_str()));
  }
  
  // Create new branch and add it to the event. It will then be there for every event, and cleared for each event.
  fThermalParticlesArray = new TClonesArray("AliBasicParticle", fMaxN+1);
  fThermalParticlesArray->SetName(fOutputCollectionName.c_str());
  fEvent->AddObject(fThermalParticlesArray);
}

/**
 * Steers each event. It enforces that the event is initialized before executing the main analysis of the event.
 */
void AliAnalysisTaskGenerateThermalBackgroundMC::UserExec(Option_t *option)
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
void AliAnalysisTaskGenerateThermalBackgroundMC::Run()
{
  // Get the TClonesArray
  fThermalParticlesArray = dynamic_cast<TClonesArray*>(fEvent->FindListObject(fOutputCollectionName.c_str()));
  if (!fThermalParticlesArray) {
    AliError(Form("%s: Could not retrieve array with name %s!", GetName(), fOutputCollectionName.c_str()));
    return;
  }
  
  // Set number of particles in event
  Int_t N = 0;
  if (fUseGaussianForN) {
    N = fRandom.Gaus(fNGaussianMean, fNGaussianSigma);
  }
  else {
    N = fRandom.Uniform(fMinN, fMaxN);
  }
  
  // Loop through particles to:
  // (1) Sample their pT from the thermal distribution, and assign a random eta/phi
  // (2) Create an AliBasicParticle, and write it to the TClonesArray
  for (Int_t i=0; i<N; i++) {
    
    Double_t eta = fRandom.Uniform(fMinEta, fMaxEta);
    Double_t phi = fRandom.Uniform(0., 2*TMath::Pi());
    Double_t pT = fPtDistribution->GetRandom();
    
    Int_t charge = 1;
    if (fRandom.Rndm() < fChargedParticleFraction) {
      if (pT < fMinChargedPt) {
        continue;
      }
    }
    else {
      charge = 0;
      if (pT < fMinNeutralPt) {
        continue;
      }
    }
    
    new((*fThermalParticlesArray)[i]) AliBasicParticle(eta, phi, pT, charge);
    
  }
  fThermalParticlesArray->Compress();
}

/**
 * Loop over particles to fill histograms.
 */
void AliAnalysisTaskGenerateThermalBackgroundMC::FillHistograms()
{
  // Initialize the various sums to 0
  Double_t ptSum = 0;
  Double_t ptSumRC02 = 0;
  Double_t ptSumRC04 = 0;
  Double_t maxPt = 0;
  
  // Generate random cone eta-phi
  Double_t jetR02 = 0.2;
  Double_t jetR04 = 0.4;
  Double_t etaRC02 = fRandom.Uniform(fMinEta + jetR02, fMaxEta - jetR02);
  Double_t phiRC02 = fRandom.Uniform(jetR02, 2*TMath::Pi() - jetR02);
  Double_t etaRC04 = fRandom.Uniform(fMinEta + jetR04, fMaxEta - jetR04);
  Double_t phiRC04 = fRandom.Uniform(jetR04, 2*TMath::Pi() - jetR04);
    
  // Loop over particles.
  Int_t N = fThermalParticlesArray->GetEntries();
  for (Int_t i=0; i<N; i++) {
    
    AliBasicParticle* thermalParticle = dynamic_cast<AliBasicParticle*>(fThermalParticlesArray->At(i));
    Double_t pT = thermalParticle->Pt();
    Double_t eta = thermalParticle->Eta();
    Double_t phi = thermalParticle->Phi();
  
    TString histname = "hPt";
    fHistManager.FillTH1(histname, pT);
    
    histname = "hEtaPhi";
    fHistManager.FillTH2(histname, eta, phi);
    
    ptSum += pT;
    
    if (pT > maxPt) {
      maxPt = pT;
    }
    
    if (GetDeltaR(eta, phi, etaRC02, phiRC02) < jetR02) {
      ptSumRC02 += pT;
    }
    if (GetDeltaR(eta, phi, etaRC04, phiRC04) < jetR04) {
      ptSumRC04 += pT;
    }
  }
  
  // Compute mean pT
  Double_t meanPt = ptSum / N;
  TString histname = "hNvsMeanPt";
  fHistManager.FillTH2(histname, N, meanPt);
  
  // Compute leading particle pT in the event
  histname = "hLeadingPt";
  fHistManager.FillTH1(histname, maxPt);
  
  // Compute the background density
  Double_t area = (fMaxEta - fMinEta) * 2 * TMath::Pi();
  Double_t backgroundDensity = ptSum / area;
  histname = "hBackgroundDensity";
  fHistManager.FillTH1(histname, backgroundDensity);
    
  // Compute delta pT
  Double_t deltaPt02 = ptSumRC02 - backgroundDensity * TMath::Pi() * jetR02 * jetR02;
  histname = "hDeltaPt02";
  fHistManager.FillTH1(histname, deltaPt02);
  
  Double_t deltaPt04 = ptSumRC04 - backgroundDensity * TMath::Pi() * jetR04 * jetR04;
  histname = "hDeltaPt04";
  fHistManager.FillTH1(histname, deltaPt04);
  
}

/**
 * Get deltaR between two points in eta/phi.
 */
Double_t AliAnalysisTaskGenerateThermalBackgroundMC::GetDeltaR(Double_t eta, Double_t phi, Double_t etaRef, Double_t phiRef)
{
  Double_t deltaPhi = TMath::Abs(phi - phiRef);
  Double_t deltaEta = TMath::Abs(eta - etaRef);
  Double_t deltaR = TMath::Sqrt( deltaPhi*deltaPhi + deltaEta*deltaEta );
  return deltaR;
}

/**
 * AddTask.
 */
AliAnalysisTaskGenerateThermalBackgroundMC* AliAnalysisTaskGenerateThermalBackgroundMC::AddTaskGenerateThermalBackgroundMC(
  const char *outputCollectionName,
  const Double_t beta,
  const char *suffix)
{
  // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr)
  {
    ::Error("AddTaskGenerateThermalBackgroundMC", "No analysis manager to connect to.");
    return 0;
  }
  
  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  AliVEventHandler* handler = mgr->GetInputEventHandler();
  if (!handler)
  {
    ::Error("AddTaskGenerateThermalBackgroundMC", "This task requires an input event handler");
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
  TString name("AliAnalysisTaskGenerateThermalBackgroundMC");
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
  AliAnalysisTaskGenerateThermalBackgroundMC* task = new AliAnalysisTaskGenerateThermalBackgroundMC(name.Data());

  task->SetOutputCollectionName(outputCollectionName);
  task->SetBeta(beta);
  
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
