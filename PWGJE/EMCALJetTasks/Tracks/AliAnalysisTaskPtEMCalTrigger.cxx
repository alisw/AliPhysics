/**************************************************************************
 * Copyright(c) 1998-2014, ALICE Experiment at CERN, All rights reserved. *
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

#include <map>
#include <cstring>
#include <iostream>
#include <memory>
#include <vector>
#include <string>
#include <sstream>

#include <TClonesArray.h>
#include <TDirectory.h>
#include <TH1.h>
#include <THashList.h>
#include <THistManager.h>
#include <TList.h>
#include <TKey.h>
#include <TMath.h>
#include <TObjArray.h>
#include <TObjString.h>
#include <TString.h>
#include <TVector2.h>

#include "AliAODEvent.h"
#include "AliESDEvent.h"
#include "AliInputEventHandler.h"
#include "AliMCEvent.h"
#include "AliParticleContainer.h"
#include "AliVCluster.h"
#include "AliVParticle.h"
#include "AliVTrack.h"
#include "AliVVertex.h"

#include "AliClusterContainer.h"
#include "AliEmcalJet.h"
#include "AliEMCALTriggerPatchInfo.h"
#include "AliEmcalPhysicsSelection.h"
#include "AliAnalysisTaskPtEMCalTrigger.h"
#include "AliEmcalTrackSelection.h"
#include "AliEmcalTrackSelectionESD.h"
#include "AliEmcalTrackSelectionAOD.h"
#include "AliEmcalCutBase.h"
#include "AliEmcalVCutsWrapper.h"
#include "AliEmcalESDtrackCutsWrapper.h"
#include "AliJetContainer.h"
#include "AliParticleContainer.h"
#include "AliPicoTrack.h"

ClassImp(PWGJE::EMCALJetTasks::AliAnalysisTaskPtEMCalTrigger)

using namespace PWGJE::EMCALJetTasks;

  /*
   * constants
   */
  const Int_t AliAnalysisTaskPtEMCalTrigger::kNJetRadii = 4;
  const Double_t jetRadVals[AliAnalysisTaskPtEMCalTrigger::kNJetRadii] = {0.2, 0.3, 0.4, 0.5};
  const Double_t *AliAnalysisTaskPtEMCalTrigger::kJetRadii = jetRadVals;

  /**
   * Dummy constructor, initialising the values with default (NULL) values
   */
  AliAnalysisTaskPtEMCalTrigger::AliAnalysisTaskPtEMCalTrigger():
                		                AliAnalysisTaskEmcalJet(),
                		                fHistos(NULL),
                		                fListTrackCuts(NULL),
                		                fEtaRange(),
                		                fPtRange(),
                		                fEnergyRange(),
                		                fVertexRange(),
                		                fJetContainersMC(),
                		                fJetContainersData(),
                		                fSelectAllTracks(kFALSE),
                		                fSwapEta(kFALSE),
                		                fUseTriggersFromTriggerMaker(kFALSE)
  {
  }

  /**
   * Main constructor, setting default values for eta and zvertex cut
   * @param name Name of the task
   */
  AliAnalysisTaskPtEMCalTrigger::AliAnalysisTaskPtEMCalTrigger(const char *name):
                		                AliAnalysisTaskEmcalJet(name, kTRUE),
                		                fHistos(NULL),
                		                fListTrackCuts(NULL),
                		                fEtaRange(),
                		                fPtRange(),
                		                fEnergyRange(),
                		                fVertexRange(),
                		                fJetContainersMC(),
                		                fJetContainersData(),
                		                fSelectAllTracks(kFALSE),
                		                fSwapEta(kFALSE),
                		                fUseTriggersFromTriggerMaker(kFALSE)
  {

    fListTrackCuts = new TList;
    fListTrackCuts->SetOwner(false);

    // Set default cuts
    fEtaRange.SetLimits(-0.8, 0.8);
    fPtRange.SetLimits(0.15, 100.);
    fEnergyRange.SetLimits(0., 1000.);
    fVertexRange.SetLimits(-10., 10.);
    SetMakeGeneralHistograms(kTRUE);
  }

  /**
   * Destructor, deleting output
   */
  AliAnalysisTaskPtEMCalTrigger::~AliAnalysisTaskPtEMCalTrigger(){
    //if(fTrackSelection) delete fTrackSelection;
    if(fHistos) delete fHistos;
    if(fListTrackCuts) delete fListTrackCuts;
  }

  /**
   * Create the list of output objects and define the histograms.
   * Also adding the track cuts to the list of histograms.
   */
  void AliAnalysisTaskPtEMCalTrigger::UserCreateOutputObjects(){
    AliAnalysisTaskEmcal::UserCreateOutputObjects();
    SetCaloTriggerPatchInfoName("EmcalTriggers");
    fHistos = new THistManager("PtEMCalTriggerHistograms");
    fHistos->ReleaseOwner();

    if(fJetContainersMC.GetEntries()){
      AliDebug(1,"Jet containers for MC truth:");
      TObjString *contname(NULL);
      TIter contMCIter(&fJetContainersMC);
      while((contname = dynamic_cast<TObjString *>(contMCIter.Next())))
        AliDebug(1, Form("Next container: %s", contname->String().Data()));
    }
    if(fJetContainersData.GetEntries()){
      AliDebug(1, "Jet containers for Data:");
      TObjString *contname(NULL);
      TIter contDataIter(&fJetContainersData);
      while((contname = dynamic_cast<TObjString *>(contDataIter.Next())))
        AliDebug(1, Form("Next container: %s", contname->String().Data()));
    }
    if(fJetCollArray.GetEntries()){
      AliDebug(1, "Jet containers attached to this task:");
      AliJetContainer *cont(NULL);
      TIter contIter(&fJetCollArray);
      while((cont = dynamic_cast<AliJetContainer *>(contIter.Next())))
        AliDebug(1, Form("Container: %s", cont->GetName()));
    }

    std::map<std::string, std::string> triggerCombinations;
    const char *triggernames[12] = {"MinBias", "EMCJHigh", "EMCJLow", "EMCGHigh", "EMCGLow", "NoEMCal", "EMCHighBoth", "EMCHighGammaOnly", "EMCHighJetOnly", "EMCLowBoth", "EMCLowGammaOnly", "EMCLowJetOnly"},
        *bitnames[4] = {"CINT7", "EMC7", "kEMCEGA", "kEMCEJE"};
    // Define axes for the trigger correlation histogram
    const TAxis *triggeraxis[5]; memset(triggeraxis, 0, sizeof(const TAxis *) * 5);
    const TAxis *bitaxes[4]; memset(bitaxes, 0, sizeof(TAxis *) * 4);
    const char *binlabels[2] = {"OFF", "ON"};
    TAxis mytrgaxis[5], mybitaxis[4];
    for(int itrg = 0; itrg < 5; ++itrg){
      DefineAxis(mytrgaxis[itrg], triggernames[itrg], triggernames[itrg], 2, -0.5, 1.5, binlabels);
      triggeraxis[itrg] = mytrgaxis+itrg;
      if(itrg < 4){
        DefineAxis(mybitaxis[itrg], bitnames[itrg], bitnames[itrg], 2, -0.5, 1.5, binlabels);
        bitaxes[itrg] = mybitaxis+itrg;
      }
    }
    // Define names and titles for different triggers in the histogram container
    triggerCombinations.insert(std::pair<std::string,std::string>(triggernames[0], "min. bias events"));
    triggerCombinations.insert(std::pair<std::string,std::string>(triggernames[1], "jet-triggered events (high threshold)"));
    triggerCombinations.insert(std::pair<std::string,std::string>(triggernames[2], "jet-triggered events (low threshold)"));
    triggerCombinations.insert(std::pair<std::string,std::string>(triggernames[3], "gamma-triggered events (high threshold)"));
    triggerCombinations.insert(std::pair<std::string,std::string>(triggernames[4], "gamma-triggered events (low threshold)"));
    triggerCombinations.insert(std::pair<std::string,std::string>(triggernames[5], "non-EMCal-triggered events"));
    triggerCombinations.insert(std::pair<std::string,std::string>(triggernames[6], "jet and gamma triggered events (high threshold)"));
    triggerCombinations.insert(std::pair<std::string,std::string>(triggernames[7], "exclusively gamma-triggered events (high threshold)"));
    triggerCombinations.insert(std::pair<std::string,std::string>(triggernames[8], "exclusively jet-triggered events (high threshold)"));
    triggerCombinations.insert(std::pair<std::string,std::string>(triggernames[9], "jet and gamma triggered events (low threshold)"));
    triggerCombinations.insert(std::pair<std::string,std::string>(triggernames[10], "exclusively gamma-triggered events (low threshold)"));
    triggerCombinations.insert(std::pair<std::string,std::string>(triggernames[11], "exclusively-triggered events (low threshold)"));
    // Define axes for the pt histogram
    // Dimensions:
    // 1. pt
    // 2. eta
    // 3. phi
    // 4. vertex
    // 5. pileup (0 = all events, 1 = after pileup rejection)
    // 6. track cuts (0 = no cuts; 1 = after std cuts)
    TArrayD ptbinning, zvertexBinning, etabinning, pileupaxis(3);
    pileupaxis[0] = -0.5; pileupaxis[1] = 0.5; pileupaxis[2] = 1.5;
    CreateDefaultPtBinning(ptbinning);
    CreateDefaultZVertexBinning(zvertexBinning);
    CreateDefaultEtaBinning(etabinning);
    TAxis htrackaxes[7];
    DefineAxis(htrackaxes[0], "pt", "p_{t} (GeV/c)", ptbinning);
    DefineAxis(htrackaxes[1], "eta", "#eta", etabinning);
    DefineAxis(htrackaxes[2], "phi", "#phi", 20, 0, 2 * TMath::Pi());
    DefineAxis(htrackaxes[3], "zvertex", "z_{V} (cm)", zvertexBinning);
    DefineAxis(htrackaxes[4], "pileup", "Pileup rejection", 2, -0.5, 1.5);
    DefineAxis(htrackaxes[5], "trackcuts", "Track Cuts", (fListTrackCuts ? fListTrackCuts->GetEntries() : 0) + 1, -0.5, (fListTrackCuts ? fListTrackCuts->GetEntries() : 0) + 0.5);
    DefineAxis(htrackaxes[6], "mbtrigger", "Has MB trigger", 2, -0.5, 1.5);
    const TAxis *trackaxes[7];
    for(int iaxis = 0; iaxis < 7; ++iaxis) trackaxes[iaxis] = htrackaxes + iaxis;
    TAxis hclusteraxes[4];
    DefineAxis(hclusteraxes[0], "energy", "E (GeV)", ptbinning);
    DefineAxis(hclusteraxes[1], "zvertex", "z_{V} (cm)", zvertexBinning);
    DefineAxis(hclusteraxes[2], "pileup", "Pileup rejection", 2,1 -0.5, 1.5);
    DefineAxis(hclusteraxes[3], "mbtrigger", "Has MB trigger", 2, -0.5, 1.5);
    const TAxis *clusteraxes[4];
    for(int iaxis = 0; iaxis < 4; ++iaxis) clusteraxes[iaxis] = hclusteraxes + iaxis;
    TAxis hpatchenergyaxes[5];
    DefineAxis(hpatchenergyaxes[0], "energy", "Patch energy (GeV)", 100, 0., 100);
    DefineAxis(hpatchenergyaxes[1], "eta", "#eta", etabinning);
    DefineAxis(hpatchenergyaxes[2], "phi", "#phi",  20, 0, 2 * TMath::Pi());
    DefineAxis(hpatchenergyaxes[3], "isMain", "Main trigger", 2, -0.5, 1.5);
    DefineAxis(hpatchenergyaxes[4],  "emcalgood", "EMCAL good event", 2, -0.5, 1.5);
    const TAxis *patchenergyaxes[5];
    for(int iaxis = 0; iaxis < 5; ++iaxis) patchenergyaxes[iaxis] = hpatchenergyaxes + iaxis;
    TAxis hpatchampaxes[5];
    DefineAxis(hpatchampaxes[0], "amplitude", "Patch energy (GeV)", 10000, 0., 10000.);
    DefineAxis(hpatchampaxes[1], "eta", "#eta", etabinning);
    DefineAxis(hpatchampaxes[2], "phi", "#phi",  20, 0, 2 * TMath::Pi());
    DefineAxis(hpatchampaxes[3], "isMain", "Main trigger", 2, -0.5, 1.5);
    DefineAxis(hpatchampaxes[4],  "emcalgood", "EMCAL good event", 2, -0.5, 1.5);
    const TAxis *patchampaxes[5];
    for(int iaxis = 0; iaxis < 5; ++iaxis) patchampaxes[iaxis] = hpatchampaxes + iaxis;
    std::string patchnames[] = {"Level0", "JetHigh", "JetLow", "GammaHigh", "GammaLow"};
    for(std::string * triggerpatch = patchnames; triggerpatch < patchnames + sizeof(patchnames)/sizeof(std::string); ++triggerpatch){
      fHistos->CreateTHnSparse(Form("Energy%s", triggerpatch->c_str()), Form("Patch energy for %s trigger patches", triggerpatch->c_str()), 5, patchenergyaxes, "s");
      fHistos->CreateTHnSparse(Form("EnergyRough%s", triggerpatch->c_str()), Form("Rough patch energy for %s trigger patches", triggerpatch->c_str()), 5, patchenergyaxes, "s");
      fHistos->CreateTHnSparse(Form("Amplitude%s", triggerpatch->c_str()), Form("Patch amplitude for %s trigger patches", triggerpatch->c_str()), 5, patchampaxes, "s");
    }

    // Create histogram for MC-truth
    fHistos->CreateTHnSparse("hMCtrueParticles", "Particle-based histogram for MC-true particles", 5, trackaxes, "s");
    if(fJetCollArray.GetEntries()){
      for(int irad = 0; irad < kNJetRadii; irad++){
        fHistos->CreateTHnSparse(Form("hMCtrueParticlesRad%02d", int(kJetRadii[irad]*10)),
            Form("Particle-based histogram for MC-true particles in Jets with radius %.1f", kJetRadii[irad]*10), 5, trackaxes, "s");
      }
      // histogram for isolated particles
      fHistos->CreateTHnSparse("hMCtrueParticlesIsolated", "Particle-based histogram for isolated MC-true particles", 5, trackaxes, "s");
    }
    for(std::map<std::string,std::string>::iterator it = triggerCombinations.begin(); it != triggerCombinations.end(); ++it){
      const std::string name = it->first, &title = it->second;
      // Create event-based histogram
      fHistos->CreateTH2(Form("hEventHist%s", name.c_str()), Form("Event-based data for %s events; pileup rejection; z_{V} (cm)", title.c_str()), pileupaxis, zvertexBinning);
      // Create track-based histogram
      fHistos->CreateTHnSparse(Form("hTrackHist%s", name.c_str()), Form("Track-based data for %s events", title.c_str()), 7, trackaxes, "s");
      fHistos->CreateTHnSparse(Form("hTrackInAcceptanceHist%s", name.c_str()), Form("Track-based data for %s events  for tracks matched to EMCal clusters", title.c_str()), 7, trackaxes, "s");
      fHistos->CreateTHnSparse(Form("hMCTrackHist%s", name.c_str()), Form("Track-based data for %s events with MC kinematics", title.c_str()), 7, trackaxes, "s");
      fHistos->CreateTHnSparse(Form("hMCTrackInAcceptanceHist%s", name.c_str()), Form("Track-based data for %s events with MC kinematics for tracks matched to EMCal clusters", title.c_str()), 7, trackaxes, "s");
      // Create cluster-based histogram (Uncalibrated and calibrated clusters)
      fHistos->CreateTHnSparse(Form("hClusterCalibHist%s", name.c_str()), Form("Calib. cluster-based histogram for %s events", title.c_str()), 4, clusteraxes, "s");
      fHistos->CreateTHnSparse(Form("hClusterUncalibHist%s", name.c_str()), Form("Uncalib. cluster-based histogram for %s events", title.c_str()), 4, clusteraxes, "s");
      if(fJetCollArray.GetEntries()){
        // Create histograms for particles in different jetfinder radii, only track-based
        for(int irad = 0; irad < kNJetRadii; irad++){
          fHistos->CreateTHnSparse(Form("hTrackHist%sRad%02d", name.c_str(), int(kJetRadii[irad]*10)),
              Form("Track-based data for %s events for tracks in jets with Radius %.1f", title.c_str(), kJetRadii[irad]), 7, trackaxes, "s");
          fHistos->CreateTHnSparse(Form("hTrackInAcceptanceHist%sRad%02d", name.c_str(), int(kJetRadii[irad]*10)),
              Form("Track-based data for %s events for tracks matched to EMCal clusters in jets with Radius %.1f", title.c_str(), kJetRadii[irad]),
              7, trackaxes, "s");
          fHistos->CreateTHnSparse(Form("hMCTrackHist%sRad%02d", name.c_str(), int(kJetRadii[irad]*10)),
              Form("Track-based data for %s events with MC kinematics for tracks in jets with Radius %.1f", title.c_str(), kJetRadii[irad]),
              7, trackaxes, "s");
          fHistos->CreateTHnSparse(Form("hMCTrackInAcceptanceHist%sRad%02d", name.c_str(), int(kJetRadii[irad]*10)),
              Form("Track-based data for %s events with MC kinematics for tracks matched to EMCal clusters in jets with Radius %.1f", title.c_str(), kJetRadii[irad]),
              7, trackaxes, "s");
          fHistos->CreateTHnSparse(Form("hClusterCalibHist%sRad%02d", name.c_str(), int(kJetRadii[irad]*10)),
              Form("Calib. cluster-based histogram for %s events for clusters in Jets with radius %.1f", title.c_str(), kJetRadii[irad]), 4, clusteraxes, "s");
        }
      }
      // Create also histograms for isolated particles
      fHistos->CreateTHnSparse(Form("hTrackHist%sIsolated", name.c_str()), Form("Track-based data for %s events for isolated tracks", title.c_str()), 7, trackaxes, "s");
      fHistos->CreateTHnSparse(Form("hTrackInAcceptanceHist%sIsolated", name.c_str()), Form("Track-based data for %s events for isolated tracks matched to EMCal clusters", title.c_str()), 7, trackaxes, "s");
      fHistos->CreateTHnSparse(Form("hMCTrackHist%sIsolated", name.c_str()), Form("Track-based data for %s events with MC kinematics for isolated tracks", title.c_str()), 7, trackaxes, "s");
      fHistos->CreateTHnSparse(Form("hMCTrackInAcceptanceHist%sIsolated", name.c_str()), Form("Track-based data for %s events with MC kinematics for isolated tracks matched to EMCal clusters", title.c_str()), 7, trackaxes, "s");
    }
    fHistos->CreateTHnSparse("hEventTriggers", "Trigger type per event", 5, triggeraxis);
    fHistos->CreateTHnSparse("hEventsTriggerbit", "Trigger bits for the different events", 4, bitaxes);
    fOutput->Add(fHistos->GetListOfHistograms());
    if(fListTrackCuts && fListTrackCuts->GetEntries()){
      TIter cutIter(fListTrackCuts);
      AliEmcalTrackSelection *cutObject(NULL);
      while((cutObject = dynamic_cast<AliEmcalTrackSelection *>(cutIter()))){
        PWG::EMCAL::AliEmcalVCutsWrapper *trackcuts = dynamic_cast<PWG::EMCAL::AliEmcalVCutsWrapper *>(cutObject->GetTrackCuts(0));
        if(trackcuts){
          PWG::EMCAL::AliEmcalESDtrackCutsWrapper *esdcuts = dynamic_cast<PWG::EMCAL::AliEmcalESDtrackCutsWrapper *>(trackcuts->GetCutObject());
          if(esdcuts) {
            AliESDtrackCuts *mycuts = esdcuts->GetTrackCuts();
            mycuts->DefineHistograms();
            fOutput->Add(mycuts);
          }
        }
      }
    }
    PostData(1, fOutput);
  }

  /**
   * Runs the event loop
   *
   * @param option Additional options
   */
  Bool_t AliAnalysisTaskPtEMCalTrigger::Run(){
    // Common checks: Have SPD vertex and primary vertex from tracks, and both need to have at least one contributor
    AliDebug(1,Form("Number of calibrated clusters: %d", fCaloClusters->GetEntries()));
    AliDebug(1,Form("Number of matched tracks: %d", fTracks->GetEntries()));
    if(fMCEvent){
      // Build always trigger strig from the trigger maker in case of MC
      fUseTriggersFromTriggerMaker = kTRUE;
    }

    Bool_t emcalGood = fInputHandler->IsEventSelected() & AliEmcalPhysicsSelection::kEmcalOk;

    // Loop over trigger patches, fill patch energy
    AliEMCALTriggerPatchInfo *triggerpatch(NULL);
    TIter patchIter(this->fTriggerPatchInfo);
    while((triggerpatch = dynamic_cast<AliEMCALTriggerPatchInfo *>(patchIter()))){
      double triggerpatchinfo[5] = {triggerpatch->GetPatchE(), triggerpatch->GetEtaGeo(), triggerpatch->GetPhiGeo(), triggerpatch->IsMainTrigger() ? 1. : 0., emcalGood ? 1. : 0.};
      double triggerpatchinfoamp[5] = {static_cast<double>(triggerpatch->GetADCAmp()), triggerpatch->GetEtaGeo(), triggerpatch->GetPhiGeo(), triggerpatch->IsMainTrigger() ? 1. : 0., emcalGood ? 1. : 0.};
      double triggerpatchinfoer[5] = {triggerpatch->GetADCAmpGeVRough(), triggerpatch->GetEtaGeo(), triggerpatch->GetPhiGeo(), triggerpatch->IsMainTrigger() ? 1. : 0., emcalGood ? 1. : 0.};
      if(triggerpatch->IsJetHigh()){
        fHistos->FillTHnSparse("EnergyJetHigh", triggerpatchinfo);
        fHistos->FillTHnSparse("AmplitudeJetHigh", triggerpatchinfoamp);
        fHistos->FillTHnSparse("EnergyRoughJetHigh", triggerpatchinfoer);
      }
      if(triggerpatch->IsJetLow()){
        fHistos->FillTHnSparse("EnergyJetLow", triggerpatchinfo);
        fHistos->FillTHnSparse("AmplitudeJetLow", triggerpatchinfoamp);
        fHistos->FillTHnSparse("EnergyRoughJetLow", triggerpatchinfoer);
      }
      if(triggerpatch->IsGammaHigh()){
        fHistos->FillTHnSparse("EnergyGammaHigh", triggerpatchinfo);
        fHistos->FillTHnSparse("AmplitudeGammaHigh", triggerpatchinfoamp);
        fHistos->FillTHnSparse("EnergyRoughGammaHigh", triggerpatchinfoer);
      }
      if(triggerpatch->IsGammaLow()){
        fHistos->FillTHnSparse("EnergyGammaLow", triggerpatchinfo);
        fHistos->FillTHnSparse("AmplitudeGammaLow", triggerpatchinfoamp);
        fHistos->FillTHnSparse("EnergyRoughGammaLow", triggerpatchinfoer);
      }
      if(triggerpatch->IsLevel0()){
        fHistos->FillTHnSparse("EnergyLevel0", triggerpatchinfo);
        fHistos->FillTHnSparse("AmplitudeLevel0", triggerpatchinfoamp);
        fHistos->FillTHnSparse("EnergyRoughLevel0", triggerpatchinfoer);
      }
    }

    const AliVVertex *vtxTracks = fInputEvent->GetPrimaryVertex(),
        *vtxSPD = GetSPDVertex();
    if(!(vtxTracks && vtxSPD)) return false;
    if(!fVertexRange.IsInRange(vtxTracks->GetZ())) return false;
    if(vtxTracks->GetNContributors() < 1 || vtxSPD->GetNContributors() < 1) return false;

    double triggers[5]; memset(triggers, 0, sizeof(double) * 5);
    double triggerbits[4]; memset(triggerbits, 0, sizeof(double) * 4);
    if(fInputHandler->IsEventSelected() & AliVEvent::kINT7){
      triggers[0] = 1.;
      triggerbits[0] = 1.;
    }

    // check triggerbits
    if(fInputHandler->IsEventSelected() & AliVEvent::kEMC7){
      triggerbits[1] = 1.;
    }
    if(fInputHandler->IsEventSelected() & AliVEvent::kEMCEGA){
      triggerbits[2] = 1.;
    }
    if(fInputHandler->IsEventSelected() & AliVEvent::kEMCEJE){
      triggerbits[3] = 1.;
    }
    fHistos->FillTHnSparse("hEventsTriggerbit", triggerbits);

    std::vector<std::string> triggerstrings;
    // EMCal-triggered event, distinguish types
    TString trgstr(fUseTriggersFromTriggerMaker ? BuildTriggerString() : fInputEvent->GetFiredTriggerClasses());
    AliDebug(1, Form("Triggerstring: %s\n", trgstr.Data()));
    if(trgstr.Contains("EJ1")){
      triggerstrings.push_back("EMCJHigh");
      triggers[1] = 1;
      if(trgstr.Contains("EG1"))
        triggerstrings.push_back("EMCHighBoth");
      else
        triggerstrings.push_back("EMCHighJetOnly");
    }
    if(trgstr.Contains("EJ2")){
      triggerstrings.push_back("EMCJLow");
      triggers[2] = 1;
      if(trgstr.Contains("EG2"))
        triggerstrings.push_back("EMCLowBoth");
      else
        triggerstrings.push_back("EMCLowJetOnly");
    }
    if(trgstr.Contains("EG1")){
      triggerstrings.push_back("EMCGHigh");
      triggers[3] = 1;
      if(!trgstr.Contains("EJ1"))
        triggerstrings.push_back("EMCHighGammaOnly");
    }
    if(trgstr.Contains("EG2")){
      triggerstrings.push_back("EMCGLow");
      triggers[4] = 1;
      if(!trgstr.Contains("EJ2"))
        triggerstrings.push_back("EMCLowGammaOnly");
    }

    fHistos->FillTHnSparse("hEventTriggers", triggers);

    // apply event selection: Combine the Pileup cut from SPD with the other pA Vertex selection cuts.
    bool isPileupEvent = fInputEvent->IsPileupFromSPD(3, 0.8, 3., 2., 5.);
    isPileupEvent = isPileupEvent || (TMath::Abs(vtxTracks->GetZ() - vtxSPD->GetZ()) > 0.5);
    double covSPD[6]; vtxSPD->GetCovarianceMatrix(covSPD);
    isPileupEvent = isPileupEvent || (vtxSPD->IsFromVertexerZ() && TMath::Sqrt(covSPD[5]) > 0.25); // selection effectively inactive in old versions of the code

    // Fill event-based histogram
    const double &zv = vtxTracks->GetZ();
    if(triggers[0]) FillEventHist("MinBias", zv, isPileupEvent);
    if(!triggerstrings.size()) // Non-EMCal-triggered
      FillEventHist("NoEMCal", zv, isPileupEvent);
    else{
      // EMCal-triggered events
      for(std::vector<std::string>::iterator it = triggerstrings.begin(); it != triggerstrings.end(); ++it)
        FillEventHist(it->c_str(), zv, isPileupEvent);
    }

    const AliEmcalJet *foundJet(NULL);
    char histname[1024];
    std::vector<double> coneradii;
    AliJetContainer *jetContMC(NULL), *jetContData(NULL);

    // Fill MC truth
    AliVParticle *part(NULL);
    if(fMCEvent){
      if(fJetCollArray.GetEntries() &&
          (jetContMC = dynamic_cast<AliJetContainer *>(fJetCollArray.FindObject((static_cast<TObjString *>(fJetContainersMC.At(0)))->String().Data())))){
        // In case we have a jet array we loop over all MC selected particles in the particle container of the jet array
        TIter particles(jetContMC->GetParticleContainer()->GetArray());
        while((part = dynamic_cast<AliVParticle *>(particles()))){
          if(part->Charge() == 0) continue;
          if(!fEtaRange.IsInRange(part->Eta())) continue;
          if(!fPtRange.IsInRange(part->Pt())) continue;
          FillMCParticleHist("hMCtrueParticles", part, zv, isPileupEvent);

          /*
           * Jet part: Find track in jet container,
           * check according to number of particles in jet, and
           * check for different cone radii
           */
          foundJet = FoundTrackInJet(part, jetContMC);
          if(foundJet && foundJet->GetNumberOfConstituents() > 1){
            for(int irad = 0; irad < kNJetRadii; irad++){
              if(IsInRadius(part, foundJet, kJetRadii[irad])){
                sprintf(histname, "hMCtrueParticlesRad%02d", int(kJetRadii[irad]*10));
                FillMCParticleHist(histname,  part, zv, isPileupEvent);
              }
            }
          } else {
            // isolated track
            FillMCParticleHist("hMCtrueParticlesIsolated", part, zv, isPileupEvent);
          }
        }
      } else {
        // Use MC Event
        for(int ipart = 0; ipart < fMCEvent->GetNumberOfTracks(); ipart++){
          // Select only physical primary particles
          part = fMCEvent->GetTrack(ipart);
          if(part->Charge() == 0) continue;
          if(!fEtaRange.IsInRange(part->Eta())) continue;
          if(!fPtRange.IsInRange(part->Pt())) continue;
          if(!fMCEvent->IsPhysicalPrimary(ipart)) continue;
          FillMCParticleHist("hMCtrueParticles", part, zv, isPileupEvent);
        }
      }
    }

    AliVTrack *track(NULL);
    AliPicoTrack *picoTrack(NULL);
    TObject *containerObject(NULL);
    // Loop over all tracks (No cuts applied)
    if(fSelectAllTracks){
      // loop over all tracks only if requested
      TIter allTrackIter(fTracks);
      while((containerObject = dynamic_cast<TObject *>(allTrackIter()))){
        if((picoTrack = dynamic_cast<AliPicoTrack *>(containerObject))){
          track = picoTrack->GetTrack();
        } else
          track = dynamic_cast<AliVTrack *>(containerObject);
        if(!IsTrueTrack(track)) continue;
        if(!fEtaRange.IsInRange(track->Eta())) continue;
        if(!fPtRange.IsInRange(track->Pt())) continue;
        if(triggers[0]) FillTrackHist("MinBias", track, zv, isPileupEvent, 0, triggers[0]);
        if(!triggerstrings.size()) // Non-EMCal-triggered
          FillTrackHist("NoEMCal", track, zv, isPileupEvent, 0, triggers[0]);
        else {
          // EMCal-triggered events
          for(std::vector<std::string>::iterator it = triggerstrings.begin(); it != triggerstrings.end(); ++it)
            FillTrackHist(it->c_str(), track, zv, isPileupEvent, 0, triggers[0]);
        }
      }
    }

    // Now apply track selection cuts
    // allow for several track selections to be tested at the same time
    // each track selection gets a different cut ID starting from 1
    // cut ID 0 is reserved for the case of no cuts
    // In case we have a jet container attached we check whether the track is
    // found in a jet, and check for different cone radii around a jet
    if(fListTrackCuts && fListTrackCuts->GetEntries()){
      for(int icut = 0; icut < fListTrackCuts->GetEntries(); icut++){
        AliEmcalTrackSelection *trackSelection = static_cast<AliEmcalTrackSelection *>(fListTrackCuts->At(icut));
        TIter trackIter(trackSelection->GetAcceptedTracks(fTracks));
        while((track = dynamic_cast<AliVTrack *>(trackIter()))){
          if(fMCEvent && !IsTrueTrack(track)) continue;   // Reject fake tracks in case of MC
          if(!fEtaRange.IsInRange(track->Eta())) continue;
          if(!fPtRange.IsInRange(track->Pt())) continue;
          coneradii.clear();
          if(this->fJetCollArray.GetEntries() &&
              (jetContData = dynamic_cast<AliJetContainer *>(this->fJetCollArray.FindObject((static_cast<TObjString *>(this->fJetContainersData.At(0)))->String().Data())))){
            foundJet = FoundTrackInJet(track, jetContData);
            if(foundJet){
              for(int irad = 0; irad < kNJetRadii; irad++){
                if(IsInRadius(track, foundJet, kJetRadii[irad])) coneradii.push_back(kJetRadii[irad]);
              }
            }
          }
          if(triggers[0]){
            FillTrackHist("MinBias", track, zv, isPileupEvent, icut + 1, triggers[0]);
            if(coneradii.size()){
              for(std::vector<double>::iterator radIter = coneradii.begin(); radIter != coneradii.end(); radIter++)
                FillTrackHist("MinBias", track, zv, isPileupEvent, icut + 1, triggers[0], *radIter);
            }
          }
          if(!triggerstrings.size()){ // Non-EMCal-triggered
            FillTrackHist("NoEMCal", track, zv, isPileupEvent, icut + 1, triggers[0]);
            for(std::vector<double>::iterator radIter = coneradii.begin(); radIter != coneradii.end(); radIter++)
              FillTrackHist("NoEMCal", track, zv, isPileupEvent, icut + 1, triggers[0], *radIter);
          } else {
            // EMCal-triggered events
            for(std::vector<std::string>::iterator it = triggerstrings.begin(); it != triggerstrings.end(); ++it){
              FillTrackHist(it->c_str(), track, zv, isPileupEvent, icut + 1, triggers[0]);
              for(std::vector<double>::iterator radIter = coneradii.begin(); radIter != coneradii.end(); radIter++)
                FillTrackHist(it->c_str(), track, zv, isPileupEvent, icut + 1, triggers[0], *radIter);
            }
          }
        }
      }
    }

    // Next step we loop over the (uncalibrated) emcal clusters and fill histograms with the cluster energy
    const AliVCluster *clust(NULL);
    for(int icl = 0; icl < fInputEvent->GetNumberOfCaloClusters(); icl++){
      clust = fInputEvent->GetCaloCluster(icl);
      if(!clust->IsEMCAL()) continue;
      if(!fEnergyRange.IsInRange(clust->E())) continue;
      if(triggers[0]) FillClusterHist("hClusterUncalibHistMinBias", clust, zv, isPileupEvent, triggers[0]);
      if(!triggerstrings.size()){	// Non-EMCal-triggered
        FillClusterHist("hClusterUncalibHistNoEMCal", clust, zv, isPileupEvent, triggers[0]);
      } else{
        for(std::vector<std::string>::iterator it = triggerstrings.begin(); it != triggerstrings.end(); ++it){
          sprintf(histname, "hClusterUncalibHist%s", it->c_str());
          FillClusterHist(histname, clust, zv, isPileupEvent, triggers[0]);
        }
      }
    }

    if(fCaloClusters){
      TIter clustIter(fCaloClusters);
      while((clust = dynamic_cast<const AliVCluster *>(clustIter()))){
        if(!clust->IsEMCAL()) continue;
        if(!fEnergyRange.IsInRange(clust->E())) continue;
        if(triggers[0]) FillClusterHist("hClusterCalibHistMinBias", clust, zv, isPileupEvent, triggers[0]);
        if(!triggerstrings.size())	// Non-EMCal-triggered
          FillClusterHist("hClusterCalibHistNoEMCal", clust, zv, isPileupEvent, triggers[0]);
        else{
          for(std::vector<std::string>::iterator it = triggerstrings.begin(); it != triggerstrings.end(); ++it){
            sprintf(histname, "hClusterCalibHist%s", it->c_str());
            FillClusterHist(histname, clust, zv, isPileupEvent, triggers[0]);
          }
        }
      }
    }

    PostData(1, fOutput);
    return true;
  }

  /**
   * Creating the default \f$ p_{t} \f$ binning.
   *
   * @param binning Array where to store the results.
   */
  void AliAnalysisTaskPtEMCalTrigger::CreateDefaultPtBinning(TArrayD &binning) const{
    std::vector<double> mybinning;
    std::map<double,double> definitions;
    definitions.insert(std::pair<double,double>(2.5, 0.1));
    definitions.insert(std::pair<double,double>(7., 0.25));
    definitions.insert(std::pair<double,double>(15., 0.5));
    definitions.insert(std::pair<double,double>(25., 1.));
    definitions.insert(std::pair<double,double>(40., 2.5));
    definitions.insert(std::pair<double,double>(50., 5.));
    definitions.insert(std::pair<double,double>(100., 10.));
    double currentval = 0;
    for(std::map<double,double>::iterator id = definitions.begin(); id != definitions.end(); ++id){
      double limit = id->first, binwidth = id->second;
      while(currentval < limit){
        currentval += binwidth;
        mybinning.push_back(currentval);
      }
    }
    binning.Set(mybinning.size());
    int ib = 0;
    for(std::vector<double>::iterator it = mybinning.begin(); it != mybinning.end(); ++it)
      binning[ib++] = *it;
  }

  /**
   * Creating default z-Vertex binning.
   *
   * @param binning Array where to store the results.
   */
  void AliAnalysisTaskPtEMCalTrigger::CreateDefaultZVertexBinning(TArrayD &binning) const {
    std::vector<double> mybinning;
    double currentval = -10;
    mybinning.push_back(currentval);
    while(currentval <= 10.){
      currentval += 5.;
      mybinning.push_back(currentval);
    }
    binning.Set(mybinning.size());
    int ib = 0;
    for(std::vector<double>::iterator it = mybinning.begin(); it != mybinning.end(); ++it)
      binning[ib++] = *it;
  }

  /**
   * Creating default z-Vertex binning.
   *
   * @param binning Array where to store the results.
   */
  void AliAnalysisTaskPtEMCalTrigger::CreateDefaultEtaBinning(TArrayD& binning) const {
    std::vector<double> mybinning;
    double currentval = -0.8;
    mybinning.push_back(currentval);
    while(currentval <= 0.8){
      currentval += 0.1;
      mybinning.push_back(currentval);
    }
    binning.Set(mybinning.size());
    int ib = 0;
    for(std::vector<double>::iterator it = mybinning.begin(); it != mybinning.end(); ++it)
      binning[ib++] = *it;
  }

  /**
   * Define an axis with a given binning
   *
   * @param axis Axis to be defined
   * @param name Name of the axis
   * @param title Title of the axis
   * @param binning axis binning
   * @param labels array of bin labels
   */
  void AliAnalysisTaskPtEMCalTrigger::DefineAxis(TAxis& axis, const char* name,
      const char* title, const TArrayD& binning, const char** labels) {
    axis.Set(binning.GetSize()-1, binning.GetArray());
    axis.SetName(name);
    axis.SetTitle(title);
    if(labels){
      for(int ib = 1; ib <= axis.GetNbins(); ++ib)
        axis.SetBinLabel(ib, labels[ib-1]);
    }
  }

  /**
   * Define an axis with number of bins from min to max
   *
   * @param axis Axis to be defined
   * @param name Name of the axis
   * @param title Title of the axis
   * @param nbins Number of bins
   * @param min lower limit of the axis
   * @param max upper limit of the axis
   * @param labels array of bin labels
   */
  void AliAnalysisTaskPtEMCalTrigger::DefineAxis(TAxis& axis, const char* name,
      const char* title, int nbins, double min, double max,
      const char** labels) {
    axis.Set(nbins, min, max);
    axis.SetName(name);
    axis.SetTitle(title);
    if(labels){
      for(int ib = 1; ib <= axis.GetNbins(); ++ib)
        axis.SetBinLabel(ib, labels[ib-1]);
    }
  }

  /**
   * Fill event-based histogram
   *
   * @param trigger name of the trigger configuration to be processed
   * @param vz z-position of the vertex
   * @param isPileup signalises if the event is flagged as pileup event
   */
  void AliAnalysisTaskPtEMCalTrigger::FillEventHist(const char* trigger,
      double vz, bool isPileup) {
    char histname[1024];
    sprintf(histname, "hEventHist%s", trigger);
    fHistos->FillTH2(histname, 0., vz);
    if(!isPileup){
      fHistos->FillTH2(histname, 1., vz);
    }
  }

  /**
   * Fill track-based histogram with corresponding information
   *
   * @param trigger name of the trigger
   * @param track ESD track selected
   * @param vz z-position of the vertex
   * @param isPileup flag event as pileup event
   * @param cut id of the cut (0 = no cut)
   */
  void AliAnalysisTaskPtEMCalTrigger::FillTrackHist(const char* trigger,
      const AliVTrack* track, double vz, bool isPileup, int cut, bool isMinBias, double jetradius) {
    double etasign = fSwapEta ? -1. : 1.;
    double data[7] = {TMath::Abs(track->Pt()), etasign * track->Eta(), track->Phi(), vz, 0, static_cast<double>(cut), isMinBias ? 1. : 0.};
    double dataMC[7] = {0., 0., 0., vz, 0, static_cast<double>(cut), isMinBias ? 1. : 0.};
    AliVParticle *assocMC(NULL);
    if(fMCEvent && (assocMC = fMCEvent->GetTrack(TMath::Abs(track->GetLabel())))){
      // Select track onl
      dataMC[0] = TMath::Abs(assocMC->Pt());
      dataMC[1] = etasign * assocMC->Eta();
      dataMC[2] = assocMC->Phi();
    }
    char histname[1024], histnameAcc[1024], histnameMC[1024], histnameMCAcc[1024];
    sprintf(histname, "hTrackHist%s", trigger);
    sprintf(histnameAcc, "hTrackInAcceptanceHist%s", trigger);
    sprintf(histnameMC, "hMCTrackHist%s", trigger);
    sprintf(histnameMCAcc, "hMCTrackInAcceptanceHist%s", trigger);
    if(jetradius > 0.){
      char *hnames[] = {histname, histnameAcc, histnameMC, histnameMCAcc};
      for(unsigned int iname = 0; iname < sizeof(hnames)/sizeof(char *); iname++){
        char *myhname = hnames[iname];
        sprintf(myhname, "%sRad%02d", myhname, int(jetradius * 10.));
      }
    }
    Bool_t isEMCAL = kFALSE;
    if(track->IsEMCAL()){
      // Check if the cluster is matched to only one track
      AliVCluster *emcclust(NULL);
      AliDebug(2, Form("cluster id: %d\n", track->GetEMCALcluster()));
      if(fCaloClusters) {
        AliDebug(2, "Using calibrated clusters");
        emcclust = dynamic_cast<AliVCluster *>(fCaloClusters->At(track->GetEMCALcluster()));
      } else {
        AliDebug(2, "Using uncalibrated clusters");
        emcclust = fInputEvent->GetCaloCluster(track->GetEMCALcluster());
      }
      if(!emcclust) AliError("Null pointer to EMCal cluster");
      if(emcclust && emcclust->GetNTracksMatched() <= 1){
        isEMCAL = kTRUE;
      }
    }
    fHistos->FillTHnSparse(histname, data);
    if(fMCEvent) fHistos->FillTHnSparse(histnameMC, dataMC);
    if(isEMCAL){
      fHistos->FillTHnSparse(histnameAcc, data);
      if(fMCEvent) fHistos->FillTHnSparse(histnameMCAcc, dataMC);
    }
    if(!isPileup){
      data[4] = 1;
      dataMC[4] = 1;
      fHistos->FillTHnSparse(histname, data);
      if(fMCEvent) fHistos->FillTHnSparse(histnameMC, dataMC);
      if(isEMCAL){
        fHistos->FillTHnSparse(histnameAcc, data);
        if(fMCEvent) fHistos->FillTHnSparse(histnameMCAcc, dataMC);
      }
    }
  }

  /**
   * Fill cluster-based histogram with corresponding information
   *
   * @param trigger name of the trigger
   * @param cluster the EMCal cluster information
   * @param vz z-position of the vertex
   * @param isPileup flag event as pileup event
   */
  void AliAnalysisTaskPtEMCalTrigger::FillClusterHist(const char* histname,
      const AliVCluster* clust, double vz, bool isPileup, bool isMinBias) {
    double data[4] =  {clust->E(), vz, 0, isMinBias ? 1. : 0.};
    fHistos->FillTHnSparse(histname, data);
    if(!isPileup){
      data[2] = 1.;
      fHistos->FillTHnSparse(histname, data);
    }
  }

  /**
   * Fill histogram for MC-true particles with the information pt, eta and phi
   *
   * @param track the Monte-Carlo track
   */
  void AliAnalysisTaskPtEMCalTrigger::FillMCParticleHist(const char *histname, const AliVParticle * const track, double vz, bool isPileup){
    double data[5] = {TMath::Abs(track->Pt()), track->Eta(), track->Phi(), vz, 0.};
    fHistos->FillTHnSparse(histname, data);
    if(!isPileup){
      data[4] = 1.;
      fHistos->FillTHnSparse(histname, data);
    }
  }

  /**
   * Check if the track has an associated MC particle, and that the particle is a physical primary
   * In case of data we do not do the selection at that step (always return true)
   *
   * @param track Track to check
   * @result true primary track (true or false)
   */
  bool AliAnalysisTaskPtEMCalTrigger::IsTrueTrack(const AliVTrack *const track) const{
    if(!fMCEvent) return true;
    AliVParticle *mcassociate = fMCEvent->GetTrack(TMath::Abs(track->GetLabel()));
    if(!mcassociate) return false;
    return fMCEvent->IsPhysicalPrimary(TMath::Abs(track->GetLabel()));
  }

  /**
   * Add new track cuts to the task
   *
   * @param trackCuts Object of type AliESDtrackCuts
   */
  void AliAnalysisTaskPtEMCalTrigger::AddESDTrackCuts(AliESDtrackCuts* trackCuts) {
    fListTrackCuts->AddLast(new AliEmcalTrackSelectionESD(trackCuts));
  }

  /**
   * Add new track cuts to the task
   *
   * @param trackCuts Object of type AliESDtrackCuts
   */
  void AliAnalysisTaskPtEMCalTrigger::AddCutsForAOD(AliESDtrackCuts* trackCuts, UInt_t filterbits) {
    fListTrackCuts->AddLast(new AliEmcalTrackSelectionAOD(trackCuts, filterbits));
  }


  /**
   * Build trigger string from the trigger maker
   *
   * @return blank-separated string of fired trigger classes
   */
  TString AliAnalysisTaskPtEMCalTrigger::BuildTriggerString() {
    AliDebug(1, "trigger checking");
    TString result = "";
    if(HasTriggerType(kJ1)) result += "EJ1 ";
    if(HasTriggerType(kJ2)) result += "EJ2 ";
    if(HasTriggerType(kG1)) result += "EG1 ";
    if(HasTriggerType(kG2)) result += "EG2 ";
    return result;
  }

  /**
   * Accessor for the SPD vertex, creating transparency for ESDs and AODs
   *
   * @return the spd vertex
   */
  const AliVVertex* AliAnalysisTaskPtEMCalTrigger::GetSPDVertex() const {
    AliESDEvent *esd = dynamic_cast<AliESDEvent *>(fInputEvent);
    if(esd){
      return esd->GetPrimaryVertexSPD();
    } else {
      AliAODEvent *aod = dynamic_cast<AliAODEvent *>(fInputEvent);
      if(aod){
        return aod->GetPrimaryVertexSPD();
      }
    }
    return NULL;
  }

  /**
   * Correlate track to reconstructed jet
   *
   * @param track particle to be checked
   * @param jets container of recontructed jets
   * @return The matched jet (NULL if not found)
   */
  const AliEmcalJet * AliAnalysisTaskPtEMCalTrigger::FoundTrackInJet(
      const AliVParticle * const track, AliJetContainer *const jets) const
  {

    const AliEmcalJet *result = NULL;
    jets->ResetCurrentID();
    const AliEmcalJet *tmp = jets->GetNextAcceptJet();
    while(tmp){
      if(TrackInJet(track, tmp, jets->GetParticleContainer())){
        result = tmp;
        break;
      }
      tmp = jets->GetNextAcceptJet();
    }
    return result;
  }

  /**
   * Check if track is in radius around a given jet
   *
   * @param track Track to check
   * @param reconstructed jet jet to probe
   * @param radius cone radius
   * @return result of the test (true if track is inside the cone radius, false otherwise)
   */
  bool AliAnalysisTaskPtEMCalTrigger::IsInRadius(const AliVParticle *const track, const AliEmcalJet *reconstructedJet, Double_t radius) const {
    return reconstructedJet->DeltaR(track) < radius;
  }

  /**
   * Check if track is in radius around a given jet
   *
   * @param track Track to check
   * @param reconstructed jet jet to probe
   * @param radius cone radius
   * @return result of the test (true if track is inside the cone radius, false otherwise)
   */
  bool AliAnalysisTaskPtEMCalTrigger::IsInRadius(const AliVCluster *const clust, const AliEmcalJet *reconstructedJet, Double_t radius) const {
    double vertexPos[3];
    fInputEvent->GetPrimaryVertex()->GetXYZ(vertexPos);
    TLorentzVector clustVect;
    clust->GetMomentum(clustVect, vertexPos);

    Double_t dPhi = reconstructedJet->Phi() - clustVect.Phi();
    Double_t dEta = reconstructedJet->Eta() - clustVect.Eta();
    dPhi = TVector2::Phi_mpi_pi ( dPhi );

    double deltaR = TMath::Sqrt ( dPhi * dPhi + dEta * dEta );
    return deltaR < radius;
  }

  /**
   * Check whether track is among the jet constituents
   *
   * @param track track to check
   * @param reconstructedJet reconstructed jet to check
   * @param tracks container with tracks used for jetfinding
   * @return true if found, false otherwise
   */
  bool AliAnalysisTaskPtEMCalTrigger::TrackInJet(const AliVParticle* const track,
      const AliEmcalJet* reconstructedJet, const AliParticleContainer * const particles) const {
    bool found = false;
    const AliPicoTrack *picotmp(NULL);
    const AliVParticle *tmp(NULL);
    for(int ipart = 0; ipart < reconstructedJet->GetNumberOfTracks(); ipart++){
      tmp = dynamic_cast<const AliVParticle *>(reconstructedJet->TrackAt(ipart, particles->GetArray()));
      if((picotmp = dynamic_cast<const AliPicoTrack *>(tmp)))   // handle pico tracks
        tmp = picotmp->GetTrack();
      if(!tmp->Compare(track)){
        found = true;
        break;
      }
    }
    return found;
  }

  /**
   * Check whether a cluster is in a radius around a given jet
   *
   * @param clust The cluster to check
   * @param reconstructedJet reconstructed jet to check
   * @return the jet containing the cluster (null otherwise)
   */
  const AliEmcalJet* AliAnalysisTaskPtEMCalTrigger::FoundClusterInJet(const AliVCluster* const clust, AliJetContainer* const jets) const {

    const AliEmcalJet *result = NULL;
    jets->ResetCurrentID();
    const AliEmcalJet *tmp = jets->GetNextAcceptJet();
    while(tmp){
      if(ClusterInJet(clust, tmp, jets->GetClusterContainer())){
        result = tmp;
        break;
      }
      tmp =jets->GetNextAcceptJet();
    }
    return result;
  }

  /**
   * Check whether cluster is among the jet constituents
   *
   * @param track track to check
   * @param reconstructedJet reconstructed jet to check
   * @param clusters the cluster container
   * @return true if found, false otherwise
   */
  bool AliAnalysisTaskPtEMCalTrigger::ClusterInJet(const AliVCluster* const clust,
      const AliEmcalJet* reconstructedJet, const AliClusterContainer* const clusters) const {
    bool found = false;
    const AliVCluster *tmp(NULL);
    for(int ipart = 0; ipart < reconstructedJet->GetNumberOfTracks(); ipart++){
      tmp = dynamic_cast<const AliVCluster *>(reconstructedJet->ClusterAt(ipart, clusters->GetArray()));
      if(!tmp->Compare(clust)){
        found = true;
        break;
      }
    }
    return found;
  }

  /**
   * Add new Jet input container to the analysis task
   *
   * @param contname Name of the container
   * @param isMC Defines whether the container is for MC truth or not
   */
  void AliAnalysisTaskPtEMCalTrigger::AddJetContainerName(const Char_t* contname, Bool_t isMC) {
    TList &mycontainer = isMC ? fJetContainersMC : fJetContainersData;
    mycontainer.Add(new TObjString(contname));
  }
