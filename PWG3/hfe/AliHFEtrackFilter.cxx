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

/* $Id$ */

//
// Track filter class 
// Apply cut steps to all tracks in one event and returns a list of
// filtered tracks
//
// Author:
//   Markus Fasel <M.Fasel@gsi.de>
//
#include <TMath.h>
#include <TString.h>

#include "AliAnalysisCuts.h"
#include "AliCFContainer.h"
#include "AliCFParticleGenCuts.h"
#include "AliCFTrackKineCuts.h"
#include "AliCFTrackIsPrimaryCuts.h"
#include "AliCFTrackQualityCuts.h"
#include "AliESDEvent.h"
#include "AliESDtrack.h"
#include "AliMCEvent.h"
#include "AliMCParticle.h"

#include "AliHFEcontainer.h"
#include "AliHFEcutStep.h"
#include "AliHFEextraCuts.h"
#include "AliHFEtrackFilter.h"
#include "AliHFEtools.h"

ClassImp(AliHFEtrackFilter)

//__________________________________________________________________
AliHFEtrackFilter::AliHFEtrackFilter(const Char_t *name) :
  TNamed(name, ""),
  fFilteredTracks(NULL),
  fCutSteps(NULL),
  fEfficiencyContainers(NULL),
  fMC(NULL),
  fMCsignal(NULL),
  fPtBins(0),
  fEtaBins(0),
  fPhiBins(0),
  fPtBinning(NULL),
  fEtaBinning(NULL),
  fPhiBinning(NULL)
{
  //
  // Constructor
  //
  fCutSteps = new TObjArray;
  fFilteredTracks = new TObjArray;
  fEfficiencyContainers = new TObjArray(4);
}

//__________________________________________________________________
AliHFEtrackFilter::AliHFEtrackFilter(const AliHFEtrackFilter &ref):
  TNamed(ref),
  fFilteredTracks(NULL),
  fCutSteps(NULL),
  fEfficiencyContainers(NULL),
  fMC(NULL),
  fMCsignal(NULL),
  fPtBins(0),
  fEtaBins(0),
  fPhiBins(0),
  fPtBinning(NULL),
  fEtaBinning(NULL),
  fPhiBinning(NULL)
{
  //
  // Copy constructor
  //
  ref.Copy(*this);
}

//__________________________________________________________________
AliHFEtrackFilter &AliHFEtrackFilter::operator=(const AliHFEtrackFilter &o){
  //
  // Assignment operator
  //
  if(this != &o){
    o.Copy(*this);
  }
  return *this;
}

//__________________________________________________________________
AliHFEtrackFilter::~AliHFEtrackFilter(){
  //
  // Destructor
  //

  // remove correction framework content
  if(fPtBinning) delete fPtBinning;
  if(fEtaBinning) delete fEtaBinning;
  if(fPhiBinning) delete fPhiBinning;
  if(TestBit(kOwnCFContainers)) 
    fEfficiencyContainers->Delete();
  delete fEfficiencyContainers;
  // remove default content
  delete fFilteredTracks;
  delete fCutSteps;
}

//__________________________________________________________________
void AliHFEtrackFilter::Copy(TObject &o) const{
  //
  // Copy content into the object o
  //
  TNamed::Copy(o);
  AliHFEtrackFilter &target = dynamic_cast<AliHFEtrackFilter &>(o);

  // Make copy
  target.fFilteredTracks = dynamic_cast<TObjArray *>(fFilteredTracks->Clone());
  target.fCutSteps = dynamic_cast<TObjArray *>(fCutSteps->Clone());
  target.fEfficiencyContainers = dynamic_cast<TObjArray *>(fEfficiencyContainers->Clone());
  target.fMC = fMC;
  target.fMCsignal = dynamic_cast<AliHFEcutStep *>(fMCsignal->Clone());
  
  target.fPtBins = fPtBins;
  target.fEtaBins = fEtaBins;
  target.fPhiBins = fPhiBins;
  target.fPtBinning = new Double_t[fPtBins];
  memcpy(target.fPtBinning, fPtBinning, sizeof(Double_t) * fPtBins);
  target.fEtaBinning = new Double_t[fEtaBins];
  memcpy(target.fEtaBinning, fEtaBinning, sizeof(Double_t) * fEtaBins);
  target.fPhiBinning = new Double_t[fPhiBins];
  memcpy(target.fPhiBinning, fPhiBinning, sizeof(Double_t) * fPhiBins);
}

//__________________________________________________________________
void AliHFEtrackFilter::InitCF(){
  //
  // Initialize correction framework container
  // No HFE container given
  // Only for testing purpose
  //
   const Char_t *cnames[4] = {"_container", "_container_signal", "_container_MC", "_container_signalMC"};
  const Char_t *ctitlesAppends[4] = {"(All Rec)", "(All Signals)", "(All Rec (MC))", "(All Signals(MC))"};

  // Create the binning if not done from outside
  if(!fPtBins){
    fPtBins = 40;
    fPtBinning = AliHFEtools::MakeLogarithmicBinning(fPtBins, 0.1, 10); 
  }
  if(!fEtaBins){
    fEtaBins = 8;
    fEtaBinning = AliHFEtools::MakeLinearBinning(fEtaBins, -0.9, 0.9); 
  }
  if(!fPhiBins){
    fPhiBins = 18;
    fPhiBinning = AliHFEtools::MakeLinearBinning(fPhiBins, 0, 2*TMath::Pi());
  }
  Double_t chargeBins[3] = {-1.1, 0., 1.1};
  Int_t nStep = fCutSteps->GetEntriesFast()+1;
  Int_t nBins[4] = {fPtBins, fEtaBins, fPhiBins + 2};
  AliCFContainer *ctmp = NULL;
  for(Int_t icont = 0; icont < 4; icont++){
    TString containername =  GetName() + TString(cnames[icont]);
    TString containertitle = TString("Container for filter ") + GetName() + TString(ctitlesAppends[icont]);
    ctmp = new AliCFContainer(containername.Data(), containertitle.Data(), nStep, 4, nBins);
    SetBit(kOwnCFContainers, kTRUE);

    // Set the binning
    ctmp->SetBinLimits(0, fPtBinning);
    ctmp->SetBinLimits(1, fEtaBinning);
    ctmp->SetBinLimits(2, fPhiBinning);
    ctmp->SetBinLimits(3, chargeBins);

    // Label variable names
    ctmp->SetVarTitle(0, "pt");
    ctmp->SetVarTitle(1, "eta");
    ctmp->SetVarTitle(2, "phi");
    ctmp->SetVarTitle(3, "charge");

    // Label step name
    ctmp->SetStepTitle(0, "No Cuts");
    AliHFEcutStep *cutStep = NULL;
    for(Int_t istep = 0; istep < fCutSteps->GetEntriesFast(); istep++){
      cutStep = dynamic_cast<AliHFEcutStep *>(fCutSteps->UncheckedAt(istep)); 
      if(cutStep) ctmp->SetStepTitle(istep + 1, cutStep->GetName());
    }
    fEfficiencyContainers->AddAt(ctmp, 0);
  }
  OwnContainers();
}

//__________________________________________________________________
void AliHFEtrackFilter::InitCF(AliHFEcontainer *cont){
  //
  // Initialize Correction Framework container 
  // Use the given HFE container
  // HFE standard way to use the correction Framework
  //
  
  const Char_t *cnames[4] = {"_container", "_container_signal", "_container_MC", "_container_signalMC"};
  const Char_t *ctitlesAppends[4] = {"(All Rec)", "(All Signals)", "(All Rec (MC))", "(All Signals(MC))"};
  Int_t nStep = fCutSteps->GetEntriesFast()+1;
  AliCFContainer *ctmp = NULL;
  for(Int_t icont = 0; icont < 4; icont++){
    TString contname = GetName() + TString(cnames[icont]); 
    TString contTitle = TString("Container for filter ") + GetName() + TString(ctitlesAppends[icont]);
    //printf("Adding container %s: %s\n", contname.Data(), contTitle.Data());
    cont->CreateContainer(contname.Data(), contTitle.Data(), nStep);
    fEfficiencyContainers->AddAt((ctmp = cont->GetCFContainer(contname.Data())), icont);
    
    // Label step name
    ctmp->SetStepTitle(0, "No Cuts");
    AliHFEcutStep *cutStep = NULL;
    for(Int_t istep = 0; istep < fCutSteps->GetEntriesFast(); istep++){
      cutStep = dynamic_cast<AliHFEcutStep *>(fCutSteps->UncheckedAt(istep)); 
      if(cutStep) ctmp->SetStepTitle(istep + 1, cutStep->GetName());
    }
  }
  ReleaseContainers();
}

//__________________________________________________________________
void AliHFEtrackFilter::FilterTracks(AliESDEvent * const event){
  //
  // Perform track filtering
  // Check each cut step one by one and select tracks which pass
  // all cuts. If the correction framework is initialized, a correction#
  // framework container will be filled for each cut step.
  //
  AliHFEcutStep *cutStep = NULL;
  AliESDtrack *track = NULL;
  AliMCParticle *mctrack = NULL;
  Double_t cont[4] = {0., 0., 0., 0.}, contMC[4] = {0., 0., 0., 0.};
  Bool_t goodTrack = kTRUE, signal = kFALSE;
  Int_t nStep = fCutSteps->GetEntriesFast();
  AliCFContainer *call = dynamic_cast<AliCFContainer *>(fEfficiencyContainers->At(0));
  AliCFContainer *cSignal = dynamic_cast<AliCFContainer *>(fEfficiencyContainers->At(1));
  AliCFContainer *callMC = dynamic_cast<AliCFContainer *>(fEfficiencyContainers->At(2));
  AliCFContainer *cSignalMC = dynamic_cast<AliCFContainer *>(fEfficiencyContainers->At(3));
  for(Int_t itrack = 0; itrack < event->GetNumberOfTracks(); itrack++){
    signal = kFALSE;
    track = event->GetTrack(itrack);
    // check Monte-Carlo Information if available 
    if(fMC){
      mctrack = dynamic_cast<AliMCParticle *>(fMC->GetTrack(TMath::Abs(track->GetLabel())));
      if(mctrack){
        //AliMCParticle *mother = dynamic_cast<AliMCParticle *>(fMC->GetTrack(mctrack->Particle()->GetFirstMother()));
        //AliInfo(Form("Label %d, Mother %d", track->GetLabel(), mother->Particle()->GetPdgCode()));
        contMC[0] = mctrack->Pt();
        contMC[1] = mctrack->Eta();
        contMC[2] = mctrack->Phi();
        contMC[3] = mctrack->Charge();
        if(fMCsignal->IsSelected(mctrack)) signal = kTRUE;
        //if(TMath::Abs(mother->Particle()->GetPdgCode()) != 443) signal = kFALSE;
        //AliInfo(Form("Signal? %s", signal ? "Yes": "No"));
      }
    }

    // Fill Array without cut
    cont[0] = track->Pt();
    cont[1] = track->Eta();
    cont[2] = track->Phi();
    cont[3] = track->Charge();
    if(call) call->Fill(cont, 0);
    if(callMC) callMC->Fill(contMC, 0);
    if(signal){
      if(cSignal) cSignal->Fill(cont, 0);
      if(cSignalMC) cSignalMC->Fill(contMC, 0);
    }
    // cut the track
    goodTrack = kTRUE;
    for(Int_t icut = 0; icut < nStep; icut++){
      cutStep = dynamic_cast<AliHFEcutStep *>(fCutSteps->UncheckedAt(icut));
      if(cutStep && (!cutStep->IsSelected(track))){
        // track cut away
        goodTrack = kFALSE;
        break;
      }

      // Track survived cut step: Fill container
      if(call) call->Fill(cont, icut + 1);
      if(callMC) callMC->Fill(contMC, icut + 1);
      if(signal){
        if(cSignal) cSignal->Fill(cont, icut + 1);
        if(cSignalMC) cSignalMC->Fill(contMC, icut + 1);
      }
    }

    // Append track to the list of filtered tracks
    if(goodTrack) fFilteredTracks->Add(track);
  }
}

//__________________________________________________________________
void AliHFEtrackFilter::GenerateCutSteps(){
  //
  // Make the default cuts
  //
  MakeCutStepRecKineITSTPC();
  MakeCutStepPrimary();
  MakeCutStepHFEITS();
  MakeCutStepHFETRD();

  MakeMCSignalCuts();
}

//__________________________________________________________________
void AliHFEtrackFilter::AddCutStep(AliHFEcutStep *step){
  //
  // Add new cut step to the filter
  //
  if(!fCutSteps->FindObject(step->GetName())) fCutSteps->Add(step);
}

//__________________________________________________________________
AliHFEcutStep *AliHFEtrackFilter::GetCutStep(Int_t istep){
  //
  // Getter for single cut step
  //
  return dynamic_cast<AliHFEcutStep *>(fCutSteps->At(istep));
}

//__________________________________________________________________
AliHFEcutStep *AliHFEtrackFilter::GetCutStep(const Char_t *name){
  //
  // Getter for single cut step (by step name)
  //
  return dynamic_cast<AliHFEcutStep *>(fCutSteps->FindObject(name));
}

//__________________________________________________________________
void AliHFEtrackFilter::Flush(){
  //
  // Empty track container
  //
  fFilteredTracks->Clear();
}

//__________________________________________________________________
void AliHFEtrackFilter::SetPtBins(Int_t nBins, Double_t *binning){
  //
  // User defined pt binning
  //
  fPtBins = nBins;
  fPtBinning = new Double_t[fPtBins + 1];
  memcpy(fPtBinning, binning, sizeof(Double_t) * nBins);
}

//__________________________________________________________________
void AliHFEtrackFilter::SetEtaBins(Int_t nBins, Double_t *binning){
  //
  // User defined eta binning
  //
  fEtaBins = nBins;
  fEtaBinning = new Double_t[fEtaBins + 1];
  memcpy(fEtaBinning, binning, sizeof(Double_t) * nBins);
}

//__________________________________________________________________
void AliHFEtrackFilter::SetPhiBins(Int_t nBins, Double_t *binning){
  //
  // User defined phi binning
  //
  fPhiBins = nBins;
  fPhiBinning = new Double_t[fPhiBins + 1];
  memcpy(fPhiBinning, binning, sizeof(Double_t) * nBins);
}

//__________________________________________________________________
AliHFEcutStep *AliHFEtrackFilter::MakeCutStepRecKineITSTPC(){
  //
  // Make the cut step for Rec Kine
  // Cut step is already included in the filter
  //
  AliHFEcutStep *fCutStep = new AliHFEcutStep("RecKineITSTPC");

  AliCFTrackQualityCuts *trackQuality = new AliCFTrackQualityCuts((Char_t *)"QualityRec", (Char_t *)"REC Track Quality Cuts");
  trackQuality->SetMinNClusterTPC(80);
  trackQuality->SetMaxChi2PerClusterTPC(3.5);
  trackQuality->SetStatus(AliESDtrack::kTPCrefit | AliESDtrack::kITSrefit);
  trackQuality->SetMaxCovDiagonalElements(2., 2., 0.5, 0.5, 2); 
  fCutStep->AddCut(trackQuality);

  AliHFEextraCuts *hfecuts = new AliHFEextraCuts("HFETPC","Extra cuts from the HFE group");
  hfecuts->SetClusterRatioTPC(0.6, AliHFEextraCuts::kFoundOverCR);
  fCutStep->AddCut(hfecuts);
  
  AliCFTrackKineCuts *kineCuts = new AliCFTrackKineCuts((Char_t *)"RecKine", (Char_t *)"REC Kine Cuts");
  kineCuts->SetPtRange(0.1, 20);
  kineCuts->SetEtaRange(-0.8, 0.8);
  fCutStep->AddCut(kineCuts); 

  AddCutStep(fCutStep);
  return fCutStep;
}

//__________________________________________________________________
AliHFEcutStep *AliHFEtrackFilter::MakeCutStepPrimary(){
  // 
  // Make cut on primaries
  // Cut step is already included in the filter
  //
  AliHFEcutStep *fCutStep = new AliHFEcutStep("Primary");

  AliCFTrackIsPrimaryCuts *primaryCut = new AliCFTrackIsPrimaryCuts((Char_t *)"PrimaryCuts", (Char_t *)"REC Primary Cuts");
  primaryCut->SetMaxDCAToVertexXY(0.5);
  primaryCut->SetMaxDCAToVertexZ(5.);
  primaryCut->SetAcceptKinkDaughters(kFALSE);
  fCutStep->AddCut(primaryCut);

  AddCutStep(fCutStep);
  return fCutStep;
}

//__________________________________________________________________
AliHFEcutStep *AliHFEtrackFilter::MakeCutStepHFEITS(){
  //
  // Add special ITS cuts
  // Cut step is already included in the filter
  //
  AliHFEcutStep *fCutStep = new AliHFEcutStep("HFEITS");

  AliHFEextraCuts *hfecuts = new AliHFEextraCuts((Char_t *)"HFEPixelsCuts",(Char_t *)"Extra cuts from the HFE group");
  hfecuts->SetRequireITSpixel(AliHFEextraCuts::kFirst);
  //hfecuts->SetCheckITSstatus(kTRUE);
  fCutStep->AddCut(hfecuts);

  AddCutStep(fCutStep);
  return fCutStep;
}

//__________________________________________________________________
AliHFEcutStep *AliHFEtrackFilter::MakeCutStepHFETRD(){
  //
  // Add special TRD cut
  // Cut step is already included in the filter
  //
  AliHFEcutStep *fCutStep = new AliHFEcutStep("HFETRD");

  AliHFEextraCuts *hfecuts = new AliHFEextraCuts("HFETRDCuts","Extra cuts from the HFE group");
  hfecuts->SetMinTrackletsTRD(0);
  fCutStep->AddCut(hfecuts);

  AddCutStep(fCutStep);
  return fCutStep;
}

//__________________________________________________________________
AliHFEcutStep *AliHFEtrackFilter::MakeMCSignalCuts(){
  //
  // Define MC Signal
  // Cut step is already included in the filter
  //
  fMCsignal = new AliHFEcutStep("MCSignal");
  AliCFParticleGenCuts *genCuts = new AliCFParticleGenCuts((Char_t *)"fCutsGenMC", (Char_t *)"Particle Generation Cuts");
  genCuts->SetRequireIsCharged();
  genCuts->SetRequireIsPrimary();
  genCuts->SetProdVtxRange2D();
  genCuts->SetProdVtxRangeX(0, 1);
  genCuts->SetProdVtxRangeY(0, 1);
  genCuts->SetRequirePdgCode(11, kTRUE);
  fMCsignal->AddCut(genCuts);
  
  AliCFTrackKineCuts *kineMCcuts = new AliCFTrackKineCuts((Char_t *)"fCutsKineMC",(Char_t *)"MC Kine Cuts");
  kineMCcuts->SetPtRange(0.1, 20.);
  kineMCcuts->SetEtaRange(-0.8, 0.8);
  fMCsignal->AddCut(kineMCcuts);

  return fMCsignal;
}

//__________________________________________________________________
void AliHFEtrackFilter::SetMC(AliMCEvent * const mc){
  //
  // Publish MC event to the single cut steps
  //
  fMC = mc;
  fMCsignal->SetMC(fMC);
  AliHFEcutStep *cs = NULL;
  for(Int_t icut = 0; icut < fCutSteps->GetEntriesFast(); icut++)
    if((cs = dynamic_cast<AliHFEcutStep *>(fCutSteps->UncheckedAt(icut)))) cs->SetMC(fMC);
}

//__________________________________________________________________
void AliHFEtrackFilter::SetRecEvent(AliVEvent *rec){
  //
  // Publish MC event to the single cut steps
  //
  AliHFEcutStep *cs = NULL;
  for(Int_t icut = 0; icut < fCutSteps->GetEntriesFast(); icut++)
    if((cs = dynamic_cast<AliHFEcutStep *>(fCutSteps->UncheckedAt(icut)))) cs->SetRecEvent(rec);
}


//__________________________________________________________________
AliCFContainer *AliHFEtrackFilter::GetEfficiencyContainer(Int_t icont){
  //
  // return EfficiencyContainer
  //
  if(icont >= 4 || icont < 0) return NULL;
  return dynamic_cast<AliCFContainer *>(fEfficiencyContainers->At(icont));
} 
