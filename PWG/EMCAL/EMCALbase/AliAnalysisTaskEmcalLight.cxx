/**************************************************************************
 * Copyright(c) 1998-2016, ALICE Experiment at CERN, All rights reserved. *
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
#include <sstream>
#include <array>
#include <memory>

#include <RVersion.h>
#include <TClonesArray.h>
#include <TList.h>
#include <TObject.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TProfile.h>
#include <TSystem.h>
#include <TFile.h>
#include <TChain.h>
#include <TKey.h>

#include "AliGenCocktailEventHeader.h"
#include "AliStack.h"
#include "AliAODEvent.h"
#include "AliAnalysisManager.h"
#include "AliCentrality.h"
#include "AliEmcalList.h"
#include "AliEMCALGeometry.h"
#include "AliESDEvent.h"
#include "AliEmcalParticle.h"
#include "AliEventplane.h"
#include "AliInputEventHandler.h"
#include "AliLog.h"
#include "AliMCParticle.h"
#include "AliVCluster.h"
#include "AliVEventHandler.h"
#include "AliVParticle.h"
#include "AliAODTrack.h"
#include "AliVCaloTrigger.h"
#include "AliGenPythiaEventHeader.h"
#include "AliGenEventHeader.h"
#include "AliAODMCHeader.h"
#include "AliMCEvent.h"
#include "AliEMCALTriggerPatchInfo.h"

#include "AliMultSelection.h"

#include "AliAnalysisTaskEmcalLight.h"

Double_t AliAnalysisTaskEmcalLight::fgkEMCalDCalPhiDivide = 4.;

ClassImp(AliAnalysisTaskEmcalLight)

AliAnalysisTaskEmcalLight::AliAnalysisTaskEmcalLight() :
  AliAnalysisTaskSE(),
  fForceBeamType(kNA),
  fGeneralHistograms(kFALSE),
  fCreateHisto(kTRUE),
  fNeedEmcalGeom(kTRUE),
  fUseBuiltinEventSelection(kFALSE),
  fCentBins(),
  fCentralityEstimation(kNewCentrality),
  fIsPythia(kFALSE),
  fIsMonteCarlo(kFALSE),
  fMCEventHeaderName(),
  fCaloCellsName(),
  fCaloTriggersName(),
  fCaloTriggerPatchInfoName(),
  fCentEst("V0M"),
  fParticleCollArray(),
  fClusterCollArray(),
  fTriggerSelectionBitMap(0),
  fMinCent(-1),
  fMaxCent(-1),
  fMinVz(-999),
  fMaxVz(999),
  fMaxVzDiff(-1),
  fMinNVertCont(0),
  fMinPtHard(-1),
  fMaxPtHard(-1),
  fMaxMinimumBiasPtHard(-1),
  fAcceptedTriggerClasses(),
  fRejectedTriggerClasses(),
  fMCRejectFilter(kFALSE),
  fPtHardAndJetPtFactor(0.),
  fPtHardAndClusterPtFactor(0.),
  fPtHardAndTrackPtFactor(0.),
  fSwitchOffLHC15oFaultyBranches(kFALSE),
  fEventSelectionAfterRun(kFALSE),
  fUseAliEmcalList(kFALSE),
  fUsePtHardBinScaling(kFALSE),
  fSelectGeneratorName(),
  fMinimumEventWeight(1e-6),
  fMaximumEventWeight(1e6),
  fInhibit(kFALSE),
  fLocalInitialized(kFALSE),
  fWarnMissingCentrality(kTRUE),
  fDataType(kAOD),
  fGeom(0),
  fCaloCells(0),
  fCaloTriggers(0),
  fTriggerPatchInfo(0),
  fCent(-1),
  fCentBin(-1),
  fEPV0(-1.0),
  fEPV0A(-1.0),
  fEPV0C(-1.0),
  fNVertCont(0),
  fNVertSPDCont(0),
  fFiredTriggerBitMap(0),
  fFiredTriggerClasses(),
  fBeamType(kNA),
  fMCHeader(0),
  fPythiaHeader(0),
  fUseXsecFromHeader(false),
  fPtHardBin(0),
  fPtHard(0),
  fNTrials(0),
  fXsection(0),
  fEventWeight(1),
  fGeneratorName(),
  fOutput(0),
  fHistograms()
{
  fVertex[0] = 0;
  fVertex[1] = 0;
  fVertex[2] = 0;
  fVertexSPD[0] = 0;
  fVertexSPD[1] = 0;
  fVertexSPD[2] = 0;
}

AliAnalysisTaskEmcalLight::AliAnalysisTaskEmcalLight(const char *name, Bool_t histo) :
  AliAnalysisTaskSE(name),
  fForceBeamType(kNA),
  fGeneralHistograms(kFALSE),
  fCreateHisto(kTRUE),
  fNeedEmcalGeom(kTRUE),
  fUseBuiltinEventSelection(kFALSE),
  fCentBins(6),
  fCentralityEstimation(kNewCentrality),
  fIsPythia(kFALSE),
  fIsMonteCarlo(kFALSE),
  fMCEventHeaderName(),
  fCaloCellsName(),
  fCaloTriggersName(),
  fCaloTriggerPatchInfoName(),
  fCentEst("V0M"),
  fParticleCollArray(),
  fClusterCollArray(),
  fTriggerSelectionBitMap(0),
  fMinCent(-1),
  fMaxCent(-1),
  fMinVz(-999),
  fMaxVz(999),
  fMaxVzDiff(-1),
  fMinNVertCont(0),
  fMinPtHard(-1),
  fMaxPtHard(-1),
  fMaxMinimumBiasPtHard(-1),
  fAcceptedTriggerClasses(),
  fRejectedTriggerClasses(),
  fMCRejectFilter(kFALSE),
  fPtHardAndJetPtFactor(0.),
  fPtHardAndClusterPtFactor(0.),
  fPtHardAndTrackPtFactor(0.),
  fSwitchOffLHC15oFaultyBranches(kFALSE),
  fEventSelectionAfterRun(kFALSE),
  fUseAliEmcalList(kFALSE),
  fUsePtHardBinScaling(kFALSE),
  fSelectGeneratorName(),
  fMinimumEventWeight(1e-6),
  fMaximumEventWeight(1e6),
  fInhibit(kFALSE),
  fLocalInitialized(kFALSE),
  fWarnMissingCentrality(kTRUE),
  fDataType(kAOD),
  fGeom(0),
  fCaloCells(0),
  fCaloTriggers(0),
  fTriggerPatchInfo(0),
  fCent(0),
  fCentBin(-1),
  fEPV0(-1.0),
  fEPV0A(-1.0),
  fEPV0C(-1.0),
  fNVertCont(0),
  fNVertSPDCont(0),
  fFiredTriggerBitMap(0),
  fFiredTriggerClasses(),
  fBeamType(kNA),
  fMCHeader(0),
  fPythiaHeader(0),
  fUseXsecFromHeader(false),
  fPtHardBin(0),
  fPtHard(0),
  fNTrials(0),
  fXsection(0),
  fEventWeight(1),
  fGeneratorName(),
  fOutput(0),
  fHistograms()
{
  fVertex[0] = 0;
  fVertex[1] = 0;
  fVertex[2] = 0;
  fVertexSPD[0] = 0;
  fVertexSPD[1] = 0;
  fVertexSPD[2] = 0;

  fCentBins[0] = 0;
  fCentBins[1] = 10;
  fCentBins[2] = 30;
  fCentBins[3] = 50;
  fCentBins[4] = 90;
  fCentBins[5] = 100;

  fAliEventCuts.OverrideAutomaticTriggerSelection(AliVEvent::kAny, true);

  if (fCreateHisto) DefineOutput(1, TList::Class());
}

AliAnalysisTaskEmcalLight::~AliAnalysisTaskEmcalLight()
{
  for (auto cont_it : fParticleCollArray) delete cont_it.second;
  for (auto cont_it : fClusterCollArray) delete cont_it.second;
}

void AliAnalysisTaskEmcalLight::UserCreateOutputObjects()
{
  if (fInhibit) {
    AliWarningStream() << "The execution of this task is inhibited. Returning." << std::endl;
    return;
  }
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (mgr) {
    AliVEventHandler *evhand = mgr->GetInputEventHandler();
    if (evhand) {
      if (evhand->InheritsFrom("AliESDInputHandler") || evhand->InheritsFrom("AliDummyHandler")) {
        fDataType = kESD;
      }
      else {
        fDataType = kAOD;
      }
    }
    else {
      AliError("Event handler not found!");
    }
  }
  else {
    AliError("Analysis manager not found!");
  }

  if (!fCreateHisto)
    return;

  OpenFile(1);
  if(fUseAliEmcalList) {
    auto emclist = new AliEmcalList;
    if(fUsePtHardBinScaling) emclist->SetUseScaling(true);
    emclist->SetNameXsec("fHistXsectionExternalFile");
    emclist->SetNameTrials("fHistTrialsExternalFile");
    fOutput = emclist;
  } else {
    fOutput = new TList();
  }
  fOutput->SetOwner(); // @suppress("Ambiguous problem")

  if (fCentralityEstimation == kNoCentrality) fCentBins.clear();

  if (!fGeneralHistograms) return;

  TH1* h = nullptr;

  if (fIsMonteCarlo) {
    auto weight_bins = GenerateLogFixedBinArray(1000, fMinimumEventWeight, fMaximumEventWeight, true);

    h = new TH1F("fHistEventsVsPtHard", "fHistEventsVsPtHard", 1000, 0, 1000);
    h->GetXaxis()->SetTitle("#it{p}_{T,hard} (GeV/#it{c})");
    h->GetYaxis()->SetTitle("events");
    fOutput->Add(h);
    fHistograms["fHistEventsVsPtHard"] = h;

    h = new TH1F("fHistTrialsVsPtHard", "fHistTrialsVsPtHard", 1000, 0, 1000);
    h->GetXaxis()->SetTitle("#it{p}_{T,hard} (GeV/#it{c})");
    h->GetYaxis()->SetTitle("trials");
    fOutput->Add(h);
    fHistograms["fHistTrialsVsPtHard"] = h;

    h = new TProfile("fHistXsection", "fHistXsection", 50, 0, 50);
    h->GetXaxis()->SetTitle("#it{p}_{T,hard} bin");
    h->GetYaxis()->SetTitle("total integrated cross section (mb)");
    fOutput->Add(h);
    fHistograms["fHistXsection"] = h;

    h = new TH1F("fHistXsectionDistribution", "fHistXsectionDistribution", 1000, &weight_bins[0]);
    h->GetXaxis()->SetTitle("total integrated cross section (mb)");
    h->GetYaxis()->SetTitle("events");
    fOutput->Add(h);
    fHistograms["fHistXsectionDistribution"] = h;

    h = new TH1F("fHistEventWeights", "fHistEventWeights", 1000, &weight_bins[0]);
    h->GetXaxis()->SetTitle("weight");
    h->GetYaxis()->SetTitle("events");
    fOutput->Add(h);
    fHistograms["fHistEventWeights"] = h;

    h = new TH2F("fHistEventWeightsVsPtHard", "fHistEventWeightsVsPtHard", 1000, 0, 1000, 1000, &weight_bins[0]);
    h->GetXaxis()->SetTitle("#it{p}_{T,hard} (GeV/#it{c})");
    h->GetYaxis()->SetTitle("event weight");
    fOutput->Add(h);
    fHistograms["fHistEventWeightsVsPtHard"] = h;

    h = new TH1F("fHistEventsVsPtHardNoSel", "fHistEventsVsPtHardNoSel", 1000, 0, 1000);
    h->GetXaxis()->SetTitle("#it{p}_{T,hard} (GeV/#it{c})");
    h->GetYaxis()->SetTitle("events");
    fOutput->Add(h);
    fHistograms["fHistEventsVsPtHardNoSel"] = h;

    h = new TH1F("fHistTrialsVsPtHardNoSel", "fHistTrialsVsPtHardNoSel", 1000, 0, 1000);
    h->GetXaxis()->SetTitle("#it{p}_{T,hard} (GeV/#it{c})");
    h->GetYaxis()->SetTitle("trials");
    fOutput->Add(h);
    fHistograms["fHistTrialsVsPtHardNoSel"] = h;

    h = new TProfile("fHistXsectionNoSel", "fHistXsectionNoSel", 50, 0, 50);
    h->GetXaxis()->SetTitle("#it{p}_{T,hard} bin");
    h->GetYaxis()->SetTitle("total integrated cross section (mb)");
    fOutput->Add(h);
    fHistograms["fHistXsectionNoSel"] = h;

    h = new TH1F("fHistXsectionDistributionNoSel", "fHistXsectionDistributionNoSel", 1000, &weight_bins[0]);
    h->GetXaxis()->SetTitle("total integrated cross section (mb)");
    h->GetYaxis()->SetTitle("events");
    fOutput->Add(h);
    fHistograms["fHistXsectionDistributionNoSel"] = h;

    h = new TH1F("fHistEventWeightsNoSel", "fHistEventWeightsNoSel", 1000, &weight_bins[0]);
    h->GetXaxis()->SetTitle("weight");
    h->GetYaxis()->SetTitle("events");
    fOutput->Add(h);
    fHistograms["fHistEventWeightsNoSel"] = h;

    h = new TH2F("fHistEventWeightsVsPtHardNoSel", "fHistEventWeightsVsPtHardNoSel", 1000, 0, 1000, 1000, &weight_bins[0]);
    h->GetXaxis()->SetTitle("#it{p}_{T,hard} (GeV/#it{c})");
    h->GetYaxis()->SetTitle("event weight");
    fOutput->Add(h);
    fHistograms["fHistEventWeightsVsPtHardNoSel"] = h;

    h = new TH1F("fHistTrialsExternalFile", "fHistTrialsExternalFile", 50, 0, 50);
    h->GetXaxis()->SetTitle("#it{p}_{T,hard} bin");
    h->GetYaxis()->SetTitle("trials");
    fOutput->Add(h);
    fHistograms["fHistTrialsExternalFile"] = h;

    h = new TH1F("fHistEventsExternalFile", "fHistEventsExternalFile", 50, 0, 50);
    h->GetXaxis()->SetTitle("#it{p}_{T,hard} bin");
    h->GetYaxis()->SetTitle("total events");
    fOutput->Add(h);
    fHistograms["fHistEventsExternalFile"] = h;

    h = new TProfile("fHistXsectionExternalFile", "fHistXsectionExternalFile", 50, 0, 50);
    h->GetXaxis()->SetTitle("#it{p}_{T,hard} bin");
    h->GetYaxis()->SetTitle("total integrated cross section (mb)");
    fOutput->Add(h);
    fHistograms["fHistXsectionExternalFile"] = h;
  }

  h = new TH1F("fHistZVertex","Z vertex position", 60, -30, 30);
  h->GetXaxis()->SetTitle("V_{#it{z}}");
  h->GetYaxis()->SetTitle("counts");
  fOutput->Add(h);
  fHistograms["fHistZVertex"] = h;

  h = new TH1F("fHistZVertexNoSel","Z vertex position (no event selection)", 60, -30, 30);
  h->GetXaxis()->SetTitle("V_{#it{z}}");
  h->GetYaxis()->SetTitle("counts");
  fOutput->Add(h);
  fHistograms["fHistZVertexNoSel"] = h;

  if (fCentralityEstimation != kNoCentrality) {
    h = new TH1F("fHistCentrality","Event centrality distribution", 100, 0, 100);
    h->GetXaxis()->SetTitle("Centrality (%)");
    h->GetYaxis()->SetTitle("counts");
    fOutput->Add(h);
    fHistograms["fHistCentrality"] = h;

    h = new TH1F("fHistCentralityNoSel","Event centrality distribution (no event selection)", 100, 0, 100);
    h->GetXaxis()->SetTitle("Centrality (%)");
    h->GetYaxis()->SetTitle("counts");
    fOutput->Add(h);
    fHistograms["fHistCentralityNoSel"] = h;
  }

  if (fForceBeamType != kpp) {
    h = new TH1F("fHistEventPlane","Event plane", 120, -TMath::Pi(), TMath::Pi());
    h->GetXaxis()->SetTitle("event plane");
    h->GetYaxis()->SetTitle("counts");
    fOutput->Add(h);
    fHistograms["fHistEventPlane"] = h;

    h = new TH1F("fHistEventPlaneNoSel","Event plane (no event selection)", 120, -TMath::Pi(), TMath::Pi());
    h->GetXaxis()->SetTitle("event plane");
    h->GetYaxis()->SetTitle("counts");
    fOutput->Add(h);
    fHistograms["fHistEventPlaneNoSel"] = h;
  }

  h = new TH1F("fHistEventRejection","Reasons to reject event",30,0,30);
#if ROOT_VERSION_CODE < ROOT_VERSION(6,4,2)
  h->SetBit(TH1::kCanRebin);
#else
  h->SetCanExtend(TH1::kAllAxes);
#endif
  std::array<std::string, 10> labels = {"PhysSel", "Evt Gen Name", "Trg class (acc)", "Trg class (rej)", "Cent", "vertex contr.", "Vz", "VzSPD", "SelPtHardBin", "MCOutlier"};
  int i = 1;
  for (auto label : labels) {
    h->GetXaxis()->SetBinLabel(i, label.c_str());
    i++;
  }
  h->GetYaxis()->SetTitle("counts");
  fOutput->Add(h);
  fHistograms["fHistEventRejection"] = h;

  h = new TH1F("fHistTriggerClasses","fHistTriggerClasses",3,0,3);
#if ROOT_VERSION_CODE < ROOT_VERSION(6,4,2)
  h->SetBit(TH1::kCanRebin);
#else
  h->SetCanExtend(TH1::kAllAxes);
#endif
  fOutput->Add(h);
  fHistograms["fHistTriggerClasses"] = h;

  h = new TH1F("fHistTriggerClassesNoSel","fHistTriggerClassesNoSel",3,0,3);
#if ROOT_VERSION_CODE < ROOT_VERSION(6,4,2)
  h->SetBit(TH1::kCanRebin);
#else
  h->SetCanExtend(TH1::kAllAxes);
#endif
  fOutput->Add(h);
  fHistograms["fHistTriggerClassesNoSel"] = h;

  h = new TH1F("fHistEventCount","fHistEventCount",2,0,2);
  h->GetXaxis()->SetBinLabel(1,"Accepted");
  h->GetXaxis()->SetBinLabel(2,"Rejected");
  h->GetYaxis()->SetTitle("counts");
  fOutput->Add(h);
  fHistograms["fHistEventCount"] = h;

    // Finish setting up AliEventCuts
  if (!fUseBuiltinEventSelection) {
    fAliEventCuts.AddQAplotsToList(fOutput);
  }

  PostData(1, fOutput);
}

Bool_t AliAnalysisTaskEmcalLight::FillGeneralHistograms(Bool_t eventSelected)
{
  if (eventSelected) {
    if (fIsMonteCarlo) {
      GetGeneralTH1("fHistEventsVsPtHard", true)->Fill(fPtHard);
      GetGeneralTH1("fHistTrialsVsPtHard", true)->Fill(fPtHard, fNTrials);
      GetGeneralTH1("fHistEventWeights", true)->Fill(fEventWeight);
      GetGeneralTH2("fHistEventWeightsVsPtHard", true)->Fill(fPtHard, fEventWeight);
      GetGeneralTH1("fHistXsectionDistribution", true)->Fill(fXsection);
      GetGeneralTProfile("fHistXsection", true)->Fill(fPtHardBin, fXsection);
    }

    GetGeneralTH1("fHistZVertex")->Fill(fVertex[2]);

    TH1* hCent = GetGeneralTH1("fHistCentrality");
    if (hCent) hCent->Fill(fCent);

    TH1* hEventPlane = GetGeneralTH1("fHistEventPlane");
    if (hEventPlane) hEventPlane->Fill(fEPV0);

    TH1* hTriggerClasses = GetGeneralTH1("fHistTriggerClasses");
    for (auto fired_trg : fFiredTriggerClasses) hTriggerClasses->Fill(fired_trg.c_str(), 1);
  }
  else {
    if (fIsMonteCarlo) {
      GetGeneralTH1("fHistEventsVsPtHardNoSel", true)->Fill(fPtHard);
      GetGeneralTH1("fHistTrialsVsPtHardNoSel", true)->Fill(fPtHard, fNTrials);
      GetGeneralTH1("fHistEventWeightsNoSel", true)->Fill(fEventWeight);
      GetGeneralTH2("fHistEventWeightsVsPtHardNoSel", true)->Fill(fPtHard, fEventWeight);
      GetGeneralTH1("fHistXsectionDistributionNoSel", true)->Fill(fXsection);
      GetGeneralTProfile("fHistXsectionNoSel", true)->Fill(fPtHardBin, fXsection);
    }

    GetGeneralTH1("fHistZVertexNoSel", true)->Fill(fVertex[2]);

    TH1* hCent = GetGeneralTH1("fHistCentralityNoSel");
    if (hCent) hCent->Fill(fCent);

    TH1* hEventPlane = GetGeneralTH1("fHistEventPlaneNoSel");
    if (hEventPlane) hEventPlane->Fill(fEPV0);

    TH1* hTriggerClasses = GetGeneralTH1("fHistTriggerClassesNoSel", true);
    for (auto fired_trg : fFiredTriggerClasses) hTriggerClasses->Fill(fired_trg.c_str(), 1);
  }

  return kTRUE;
}

void AliAnalysisTaskEmcalLight::UserExec(Option_t *option)
{
  if (fInhibit) {
    AliWarningStream() << "The execution of this task is inhibited. Returning." << std::endl;
    return;
  }

  if (!fLocalInitialized) ExecOnce();

  if (!fLocalInitialized) return;

  if (!RetrieveEventObjects()) return;

  Bool_t eventSelected = IsEventSelected();

  if (fGeneralHistograms && fCreateHisto) {
    if (eventSelected) {
      GetGeneralTH1("fHistEventCount", true)->Fill("Accepted",1);
    }
    else {
      GetGeneralTH1("fHistEventCount", true)->Fill("Rejected",1);
    }

    FillGeneralHistograms(kFALSE);
    if (eventSelected) FillGeneralHistograms(kTRUE);
  }

  Bool_t runOk = kFALSE;
  if (eventSelected || fEventSelectionAfterRun) runOk = Run();

  if (fCreateHisto && eventSelected && runOk) FillHistograms();

  if (fCreateHisto && fOutput) {
    // information for this iteration of the UserExec in the container
    PostData(1, fOutput);
  }
}

Bool_t AliAnalysisTaskEmcalLight::PythiaInfoFromFile(const char* currFile, Float_t &xsec, Float_t &trials, Int_t &pthard, Bool_t &useXsecFromHeader)
{

  TString file(currFile);
  xsec = 0;
  trials = 1;

  // Determine archive type
  TString archivetype;
  std::unique_ptr<TObjArray> walk(file.Tokenize("/"));
  for(auto t : *walk){
    TString &tok = static_cast<TObjString *>(t)->String();
    if(tok.Contains(".zip")){
      archivetype = tok;
      Int_t pos = archivetype.Index(".zip");
      archivetype.Replace(pos, archivetype.Length() - pos, "");
    }
  }
  if(archivetype.Length()){
    AliDebugStream(1) << "Auto-detected archive type " << archivetype << std::endl;
    Ssiz_t pos1 = file.Index(archivetype,archivetype.Length(),0,TString::kExact);
    Ssiz_t pos = file.Index("#",1,pos1,TString::kExact);
    Ssiz_t pos2 = file.Index(".root",5,TString::kExact);
    file.Replace(pos+1,pos2-pos1,"");
  } else {
    // not an archive take the basename....
    file.ReplaceAll(gSystem->BaseName(file.Data()),"");
  }
  AliDebugStream(1) << "File name: " << file << std::endl;

  // Build virtual file name
  // Support for train tests
  TString virtualFileName;
  if(file.Contains("__alice")){
    TString tmp(file);
    Int_t pos = tmp.Index("__alice");
    tmp.Replace(0, pos, "");
    tmp.ReplaceAll("__", "/");
    // cut out tag for archive and root file
    // this needs a determin
    std::unique_ptr<TObjArray> toks(tmp.Tokenize("/"));
    TString tag = "_" + archivetype;
    for(auto t : *toks){
      TString &path = static_cast<TObjString *>(t)->String();
      if(path.Contains(tag)){
        Int_t posTag = path.Index(tag);
        path.Replace(posTag, path.Length() - posTag, "");
      }
      virtualFileName += "/" + path;
    }
  } else {
    virtualFileName = file;
  }

  AliDebugStream(1) << "Physical file name " << file << ", virtual file name " << virtualFileName << std::endl;

  // Get the pt hard bin
  TString strPthard(virtualFileName);

  /*
  // Dead code - to be removed after testing phase
  // Procedure will fail for everything else than the expected path name
  strPthard.Remove(strPthard.Last('/'));
  strPthard.Remove(strPthard.Last('/'));
  if (strPthard.Contains("AOD")) strPthard.Remove(strPthard.Last('/'));    
  strPthard.Remove(0,strPthard.Last('/')+1);
  if (strPthard.IsDec()) pthard = strPthard.Atoi();
  else 
    AliWarningStream() << "Could not extract file number from path " << strPthard << std::endl;
  */

  // New implementation : pattern matching
  // Reason: Implementation valid only for old productions (new productions swap run number and pt-hard bin)
  // Idea: Don't use the position in the string but the match different informations
  // + Year clearly 2000+
  // + Run number can be match to the one in the event
  // + If we know it is not year or run number, it must be the pt-hard bin if we start from the beginning
  // The procedure is only valid for the current implementations and unable to detect non-pt-hard bins
  // It will also fail in case of arbitrary file names

  bool binfound = false;
  std::unique_ptr<TObjArray> tokens(strPthard.Tokenize("/"));
  for(auto t : *tokens) {
    TString &tok = static_cast<TObjString *>(t)->String();
    if(tok.IsDec()){
      Int_t number = tok.Atoi();
      if(number > 2000 && number < 3000){
        // Year
        continue;
      } else if(number == fInputHandler->GetEvent()->GetRunNumber()){
        // Run number
        continue;
      } else {
        if(!binfound){
          // the first number that is not one of the two must be the pt-hard bin
          binfound = true;
          pthard = number;
          break;
        }
      }
    }
  }
  if(!binfound) {
    AliErrorStream() << "Could not extract file number from path " << strPthard << std::endl;
  } else {
    AliInfoStream() << "Auto-detecting pt-hard bin " << pthard << std::endl;
  }

  AliInfoStream() << "File: " << file << std::endl;

  // problem that we cannot really test the existance of a file in a archive so we have to live with open error message from root
  std::unique_ptr<TFile> fxsec(TFile::Open(Form("%s%s",file.Data(),"pyxsec.root")));

  if (!fxsec) {
    // next trial fetch the histgram file
    fxsec = std::unique_ptr<TFile>(TFile::Open(Form("%s%s",file.Data(),"pyxsec_hists.root")));
    if (!fxsec){
      AliErrorStream() << "Failed reading cross section from file " << file << std::endl;
      useXsecFromHeader = true;
      return kFALSE; // not a severe condition but inciate that we have no information
    }
    else {
      // find the tlist we want to be independtent of the name so use the Tkey
      TKey* key = (TKey*)fxsec->GetListOfKeys()->At(0); 
      if (!key) return kFALSE;
      TList *list = dynamic_cast<TList*>(key->ReadObj());
      if (!list) return kFALSE;
      TProfile *xSecHist = static_cast<TProfile*>(list->FindObject("h1Xsec"));
      // check for failure
      if(!xSecHist->GetEntries()) {
        // No cross seciton information available - fall back to raw
        AliErrorStream() << "No cross section information available in file " << fxsec->GetName() <<" - fall back to cross section in PYTHIA header" << std::endl;
        useXsecFromHeader = true;
      } else {
        // Cross section histogram filled - take it from there
        xsec = xSecHist->GetBinContent(1);
        if(!xsec) AliErrorStream() << GetName() << ": Cross section 0 for file " << file << std::endl;
        useXsecFromHeader = false;
      }
      trials  = ((TH1F*)list->FindObject("h1Trials"))->GetBinContent(1);
    }
  } else { // no tree pyxsec.root
    TTree *xtree = (TTree*)fxsec->Get("Xsection");
    if (!xtree) return kFALSE;
    UInt_t   ntrials  = 0;
    Double_t  xsection  = 0;
    xtree->SetBranchAddress("xsection",&xsection);
    xtree->SetBranchAddress("ntrials",&ntrials);
    xtree->GetEntry(0);
    trials = ntrials;
    xsec = xsection;
  }
  return kTRUE;
}

Bool_t AliAnalysisTaskEmcalLight::UserNotify()
{
  if (!fIsPythia) return kTRUE;

  TTree *tree = AliAnalysisManager::GetAnalysisManager()->GetTree();
  if (!tree) {
    AliError(Form("%s - UserNotify: No current tree!",GetName()));
    return kFALSE;
  }

  Float_t xsection    = 0;
  Float_t trials      = 0;
  Int_t   pthardbin   = 0;
  Bool_t  useXsecFromHeader = false;

  TFile *curfile = tree->GetCurrentFile();
  if (!curfile) {
    AliError(Form("%s - UserNotify: No current file!",GetName()));
    return kFALSE;
  }

  TChain *chain = dynamic_cast<TChain*>(tree);
  if (chain) tree = chain->GetTree();

  Int_t nevents = tree->GetEntriesFast();

  Bool_t res = PythiaInfoFromFile(curfile->GetName(), xsection, trials, pthardbin, useXsecFromHeader);

  fPtHardBin = pthardbin >= 0 ? pthardbin : 0;
  fUseXsecFromHeader = useXsecFromHeader;

  if (!res) return kTRUE;

  if (fGeneralHistograms  && fCreateHisto) {
    GetGeneralTH1("fHistTrialsExternalFile", true)->Fill(fPtHardBin, trials);
    if(!useXsecFromHeader) GetGeneralTProfile("fHistXsectionExternalFile", true)->Fill(fPtHardBin, xsection);
    GetGeneralTH1("fHistEventsExternalFile", true)->Fill(fPtHardBin, nevents);
  }
  
  return kTRUE;
}

void AliAnalysisTaskEmcalLight::ExecOnce()
{
  if (!InputEvent()) {
    AliError(Form("%s: Could not retrieve event! Returning!", GetName()));
    return;
  }

  if (fNeedEmcalGeom && !fGeom) {
    fGeom = AliEMCALGeometry::GetInstanceFromRunNumber(InputEvent()->GetRunNumber());
    if (!fGeom) {
      AliFatal(Form("%s: Can not get EMCal geometry instance. If you do not need the EMCal geometry, disable it by setting task->SetNeedEmcalGeometry(kFALSE).", GetName()));
      return;
    }
  }

  if (fSwitchOffLHC15oFaultyBranches && dynamic_cast<AliAODEvent*>(InputEvent())) {
    TTree *aodTree = AliAnalysisManager::GetAnalysisManager()->GetTree();
    aodTree->SetBranchStatus("D0toKpi.fPx", 0);
    aodTree->SetBranchStatus("D0toKpi.fPy", 0);
    aodTree->SetBranchStatus("D0toKpi.fPz", 0);
    aodTree->SetBranchStatus("D0toKpi.fd0", 0);
    aodTree->SetBranchStatus("Charm3Prong.fPx", 0);
    aodTree->SetBranchStatus("Charm3Prong.fPy", 0);
    aodTree->SetBranchStatus("Charm3Prong.fPz", 0);
    aodTree->SetBranchStatus("Charm3Prong.fd0", 0);
    aodTree->SetBranchStatus("Dstar.fPx", 0);
    aodTree->SetBranchStatus("Dstar.fPy", 0);
    aodTree->SetBranchStatus("Dstar.fPz", 0);
    aodTree->SetBranchStatus("Dstar.fd0", 0);
  }

  //Load all requested track branches - each container knows name already
  for (auto cont_it : fParticleCollArray) cont_it.second->SetArray(InputEvent());

  //Load all requested cluster branches - each container knows name already
  for (auto cont_it : fClusterCollArray) cont_it.second->SetArray(InputEvent());

  if (!fCaloCellsName.IsNull() && !fCaloCells) {
    fCaloCells =  dynamic_cast<AliVCaloCells*>(InputEvent()->FindListObject(fCaloCellsName));
    if (!fCaloCells) {
      AliError(Form("%s: Could not retrieve cells %s!", GetName(), fCaloCellsName.Data()));
      return;
    }
  }

  if (!fCaloTriggersName.IsNull() && !fCaloTriggers) {
    fCaloTriggers =  dynamic_cast<AliVCaloTrigger*>(InputEvent()->FindListObject(fCaloTriggersName));
    if (!fCaloTriggers) {
      AliError(Form("%s: Could not retrieve calo triggers %s!", GetName(), fCaloTriggersName.Data()));
      return;
    }
  }

  if (!fCaloTriggerPatchInfoName.IsNull() && !fTriggerPatchInfo) {
    fTriggerPatchInfo = GetArrayFromEvent(fCaloTriggerPatchInfoName.Data(),"AliEMCALTriggerPatchInfo");
    if (!fTriggerPatchInfo) {
      AliError(Form("%s: Could not retrieve calo trigger patch info %s!", GetName(), fCaloTriggerPatchInfoName.Data()));
      return;
    }

  }

  fLocalInitialized = kTRUE;
}

AliAnalysisTaskEmcalLight::EBeamType_t AliAnalysisTaskEmcalLight::GetBeamType()
{
  if (fForceBeamType != kNA)
    return fForceBeamType;

  AliESDEvent *esd = dynamic_cast<AliESDEvent*>(InputEvent());
  if (esd) {
    const AliESDRun *run = esd->GetESDRun();
    TString beamType = run->GetBeamType();
    if (beamType == "p-p")
      return kpp;
    else if (beamType == "A-A")
      return kAA;
    else if (beamType == "p-A")
      return kpA;
    else
      return kNA;
  } else {
    Int_t runNumber = InputEvent()->GetRunNumber();
    // All run number ranges taken from the RCT
    if ((runNumber >= 136833 && runNumber <= 139517) ||  // LHC10h
        (runNumber >= 167693 && runNumber <= 170593) || // LHC11h
        (runNumber >= 244824 && runNumber <= 246994)) { // LHC15o
      return kAA;
    } else if ((runNumber >= 188356 && runNumber <= 188366) ||   // LHC12g
               (runNumber >= 195164 && runNumber <= 197388) ||  // LHC13b-f
               (runNumber >= 265015 && runNumber <= 267166)) {  // LHC16q-t
      return kpA;
    } else {
      return kpp;
    }
  }
}

Bool_t AliAnalysisTaskEmcalLight::IsEventSelected(){
  if(!IsTriggerSelected()) return false;
  if(fUseBuiltinEventSelection) return IsEventSelectedInternal();
  if(!CheckMCOutliers()) return false;
  return fAliEventCuts.AcceptEvent(fInputEvent);
}

Bool_t AliAnalysisTaskEmcalLight::IsEventSelectedInternal()
{
  TH1* hEventRejection = GetGeneralTH1("fHistEventRejection", true);

  if (fTriggerSelectionBitMap != 0 && (fFiredTriggerBitMap & fTriggerSelectionBitMap) == 0) {
    if (fGeneralHistograms) hEventRejection->Fill("PhysSel",1);
    return kFALSE;
  }

  if (!fSelectGeneratorName.IsNull() && !fGeneratorName.IsNull()) {
    if (!fGeneratorName.Contains(fSelectGeneratorName)) {
      if (fGeneralHistograms) hEventRejection->Fill("Evt Gen Name",1);
      return kFALSE;
    }
  }

  if (fMinCent < fMaxCent && fMaxCent > 0) {
    if (fCent < fMinCent || fCent > fMaxCent) {
      if (fGeneralHistograms) hEventRejection->Fill("Cent",1);
      return kFALSE;
    }
  }

  if (fNVertCont < fMinNVertCont) {
    if (fGeneralHistograms) hEventRejection->Fill("vertex contr.",1);
    return kFALSE;
  }

  if (fMinVz < fMaxVz) {
    if (fVertex[2] < fMinVz || fVertex[2] > fMaxVz) {
      if (fGeneralHistograms) hEventRejection->Fill("Vz",1);
      return kFALSE;
    }
  }

  if (fMaxVzDiff >= 0) {
    if (fNVertSPDCont > 0) {
      Double_t vzSPD = fVertexSPD[2];
      Double_t dvertex = TMath::Abs(fVertex[2] - vzSPD);
      //if difference larger than fZvertexDiff
      if (dvertex > fMaxVzDiff) {
        if (fGeneralHistograms) hEventRejection->Fill("VzSPD",1);
        return kFALSE;
      }
    }
  }

  if (fMinPtHard >= 0 && fPtHard < fMinPtHard)  {
    if (fGeneralHistograms) hEventRejection->Fill("SelPtHardBin",1);
    return kFALSE;
  }

  if (fMaxPtHard >= 0 && fPtHard >= fMaxPtHard)  {
    if (fGeneralHistograms) hEventRejection->Fill("SelPtHardBin",1);
    return kFALSE;
  }

  if (fPtHardBin == 0 && fMaxMinimumBiasPtHard >= 0 && fPtHard > fMaxMinimumBiasPtHard) {
    if (fGeneralHistograms) hEventRejection->Fill("SelPtHardBin",1);
    return kFALSE;
  }

  // Reject filter for MC data
  if (!CheckMCOutliers()) {
    if (fGeneralHistograms) hEventRejection->Fill("MCOutlier",1);
    return kFALSE;
  }

  return kTRUE;
}

Bool_t AliAnalysisTaskEmcalLight::IsTriggerSelected(){
  TH1* hEventRejection = GetGeneralTH1("fHistEventRejection", true);
  Bool_t acceptedTrgClassFound = kFALSE;
  if (fAcceptedTriggerClasses.size() > 0) {
    for (const auto &acc_trg : fAcceptedTriggerClasses) {
      std::string teststring(acc_trg);
      bool fullmatch(false);
      auto posexact = acc_trg.find("EXACT");
      if(posexact != std::string::npos) {
        fullmatch = true;
        teststring.erase(posexact, 5);
      }
      for (const auto &fired_trg : fFiredTriggerClasses) {
        bool classmatch = fullmatch ? teststring == fired_trg : fired_trg.find(teststring) != std::string::npos;
        if (classmatch) {
          acceptedTrgClassFound = kTRUE;
          break;
        }
      }
      if (acceptedTrgClassFound) break;
    }

    if (!acceptedTrgClassFound) {
      if (fGeneralHistograms) hEventRejection->Fill("Trg class (acc)",1);
      return kFALSE;
    }
  }

  if (fRejectedTriggerClasses.size() > 0) {
    for (const auto &rej_trg : fRejectedTriggerClasses) {
      std::string teststring(rej_trg);
      bool fullmatch(false);
      auto posexact = rej_trg.find("EXACT");
      if(posexact != std::string::npos) {
        fullmatch = true;
        teststring.erase(posexact, 5);
      }
      for (const auto &fired_trg : fFiredTriggerClasses) {
        bool classmatch = fullmatch ? teststring == fired_trg : fired_trg.find(teststring) != std::string::npos;
        if (classmatch) {
          if (fGeneralHistograms) hEventRejection->Fill("Trg class (rej)",1);
          return kFALSE;
        }
      }
    }
  }
  return kTRUE;
}

TClonesArray *AliAnalysisTaskEmcalLight::GetArrayFromEvent(const char *name, const char *clname)
{
  TClonesArray *arr = 0;
  TString sname(name);
  if (!sname.IsNull()) {
    arr = dynamic_cast<TClonesArray*>(InputEvent()->FindListObject(sname));
    if (!arr) {
      AliWarning(Form("%s: Could not retrieve array with name %s!", GetName(), name));
      return 0;
    }
  } else {
    return 0;
  }

  if (!clname)
    return arr;

  TString objname(arr->GetClass()->GetName());
  TClass cls(objname);
  if (!cls.InheritsFrom(clname)) {
    AliWarning(Form("%s: Objects of type %s in %s are not inherited from %s!",
        GetName(), cls.GetName(), name, clname));
    return 0;
  }
  return arr;
}

Bool_t AliAnalysisTaskEmcalLight::RetrieveEventObjects()
{
  fVertex[0] = 0;
  fVertex[1] = 0;
  fVertex[2] = 0;
  fNVertCont = 0;

  fVertexSPD[0] = 0;
  fVertexSPD[1] = 0;
  fVertexSPD[2] = 0;
  fNVertSPDCont = 0;

  fFiredTriggerClasses.clear();
  std::stringstream firedClasses(InputEvent()->GetFiredTriggerClasses().Data());
  while (firedClasses.good()) {
    std::string trgClass;
    firedClasses >> trgClass;
    if (!trgClass.empty()) fFiredTriggerClasses.push_back(trgClass);
  }

  if (fDataType == kESD) {
    fFiredTriggerBitMap = static_cast<AliInputEventHandler*>(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler())->IsEventSelected();
  }
  else {
    fFiredTriggerBitMap = static_cast<AliVAODHeader*>(InputEvent()->GetHeader())->GetOfflineTrigger();
  }

  const AliVVertex *vert = InputEvent()->GetPrimaryVertex();
  if (vert) {
    vert->GetXYZ(fVertex);
    fNVertCont = vert->GetNContributors();
  }

  const AliVVertex *vertSPD = InputEvent()->GetPrimaryVertexSPD();
  if (vertSPD) {
    vertSPD->GetXYZ(fVertexSPD);
    fNVertSPDCont = vertSPD->GetNContributors();
  }

  fBeamType = GetBeamType();

  fCent    = 99;
  fCentBin = -1;
  fEPV0    = -999;
  fEPV0A   = -999;
  fEPV0C   = -999;

  if (fCentralityEstimation == kNewCentrality) {
    // New centrality estimation (AliMultSelection)
    // See https://twiki.cern.ch/twiki/bin/viewauth/ALICE/AliMultSelectionCalibStatus for calibration status period-by-period)
    AliMultSelection *MultSelection = static_cast<AliMultSelection*>(InputEvent()->FindListObject("MultSelection"));
    if (MultSelection) {
      fCent = MultSelection->GetMultiplicityPercentile(fCentEst.Data());
    }
    else {
      if(fWarnMissingCentrality) AliWarning(Form("%s: Could not retrieve centrality information! Assuming 99", GetName()));
    }
  }
  else if (fCentralityEstimation == kOldCentrality) {
    // Old centrality estimation (AliCentrality, works only on Run-1 PbPb and pPb)
    AliCentrality *aliCent = InputEvent()->GetCentrality();
    if (aliCent) {
      fCent = aliCent->GetCentralityPercentile(fCentEst.Data());
    }
    else {
      if(fWarnMissingCentrality) AliWarning(Form("%s: Could not retrieve centrality information! Assuming 99", GetName()));
    }
  }
  if (!fCentBins.empty() && fCentralityEstimation != kNoCentrality) {
    for (auto cent_it = fCentBins.begin(); cent_it != fCentBins.end() - 1; cent_it++) {
      if (fCent >= *cent_it && fCent < *(cent_it+1)) fCentBin = cent_it - fCentBins.begin();
    }
  }
  else {
    fCentBin = 0;
  }

  if (fBeamType == kAA || fBeamType == kpA ) {
    AliEventplane *aliEP = InputEvent()->GetEventplane();
    if (aliEP) {
      fEPV0  = aliEP->GetEventplane("V0" ,InputEvent());
      fEPV0A = aliEP->GetEventplane("V0A",InputEvent());
      fEPV0C = aliEP->GetEventplane("V0C",InputEvent());
    } else {
      AliWarning(Form("%s: Could not retrieve event plane information!", GetName()));
    }
  }

  if (fIsMonteCarlo && MCEvent()) {
    AliGenEventHeader* header = MCEvent()->GenEventHeader();
    if (fMCEventHeaderName.IsNull()) {
      fMCHeader = header;
    }
    else {
      if (header->InheritsFrom(fMCEventHeaderName)) {
        fMCHeader = header;
      }
      else if (header->InheritsFrom("AliGenCocktailEventHeader")) {
        AliGenCocktailEventHeader* cocktailHeader = static_cast<AliGenCocktailEventHeader*>(header);
        TList* headers = cocktailHeader->GetHeaders();
        for (auto obj : *headers) { // @suppress("Symbol is not resolved")
          if (obj->InheritsFrom(fMCEventHeaderName)){
            fMCHeader = static_cast<AliGenEventHeader*>(obj);
            break;
          }
        }
      }
    }
    if (fMCHeader) {
      fEventWeight = fMCHeader->EventWeight();
      if (fIsPythia) {
        fPythiaHeader = static_cast<AliGenPythiaEventHeader*>(fMCHeader);
        fPtHard = fPythiaHeader->GetPtHard();
        fXsection = fPythiaHeader->GetXsection();
        fNTrials = fPythiaHeader->Trials();
        if(fUseXsecFromHeader) GetGeneralTProfile("fHistXsectionExternalFile", true)->Fill(fPtHardBin, fXsection);
      }
    }
  }

  for (auto cont_it : fParticleCollArray) cont_it.second->NextEvent(InputEvent());
  for (auto cont_it : fClusterCollArray) cont_it.second->NextEvent(InputEvent());

  return kTRUE;
}

AliParticleContainer* AliAnalysisTaskEmcalLight::AddParticleContainer(EMCAL_STRINGVIEW branchName, EMCAL_STRINGVIEW contName)
{
  if (branchName.size() == 0) return 0;

  AliParticleContainer* cont = 0;

#if ROOT_VERSION_CODE > ROOT_VERSION(6,10,0) 
#define EMCAL_STRINGVIEW_NONCONST std::string_view
#else 
#define EMCAL_STRINGVIEW_NONCONST std::string
#endif

  if (branchName == EMCAL_STRINGVIEW_NONCONST("tracks") || branchName == EMCAL_STRINGVIEW_NONCONST("Tracks")) cont = new AliTrackContainer(branchName.data());
  else if (branchName == EMCAL_STRINGVIEW_NONCONST("mcparticles")) cont = new AliMCParticleContainer(branchName.data());
  else cont = new AliParticleContainer(branchName.data());

  if (contName.size() > 0) cont->SetName(contName.data());

  AdoptParticleContainer(cont);

  return cont;
}

AliClusterContainer* AliAnalysisTaskEmcalLight::AddClusterContainer(EMCAL_STRINGVIEW branchName, EMCAL_STRINGVIEW contName)
{
  if (branchName.size() == 0) return 0;

  AliClusterContainer* cont = new AliClusterContainer(branchName.data());

  if (contName.size() > 0) cont->SetName(contName.data());

  AdoptClusterContainer(cont);

  return cont;
}

AliParticleContainer* AliAnalysisTaskEmcalLight::GetParticleContainer(EMCAL_STRINGVIEW name) const
{
  std::map<std::string, AliParticleContainer*>::const_iterator cont_it = fParticleCollArray.find(std::string(name));
  if (cont_it != fParticleCollArray.end()) return cont_it->second;
  else return nullptr;
}

AliClusterContainer* AliAnalysisTaskEmcalLight::GetClusterContainer(EMCAL_STRINGVIEW name) const
{
  std::map<std::string, AliClusterContainer*>::const_iterator cont_it = fClusterCollArray.find(std::string(name));
  if (cont_it != fClusterCollArray.end()) return cont_it->second;
  else return nullptr;
}

void AliAnalysisTaskEmcalLight::AddObjectToEvent(TObject *obj, Bool_t attempt)
{
  if (!(InputEvent()->FindListObject(obj->GetName()))) {
    InputEvent()->AddObject(obj);
  }
  else {
    if (!attempt) {
      AliFatal(Form("%s: Container with name %s already present. Aborting", GetName(), obj->GetName()));
    }
  }
}

Bool_t AliAnalysisTaskEmcalLight::IsTrackInEmcalAcceptance(AliVParticle* part, Double_t edges) const
{

  if (!fGeom) {
    AliWarning(Form("%s - AliAnalysisTaskEmcalBase::IsTrackInEmcalAcceptance - Geometry is not available!", GetName()));
    return kFALSE;
  }

  Double_t minPhi = fGeom->GetArm1PhiMin() - edges;
  Double_t maxPhi = fGeom->GetArm1PhiMax() + edges;

  if (part->Phi() > minPhi && part->Phi() < maxPhi) {
    return kTRUE;
  }
  else {
    return kFALSE;
  }
}

void AliAnalysisTaskEmcalLight::SetRejectionReasonLabels(TAxis* axis)
{
  axis->SetBinLabel(1,  "NullObject");
  axis->SetBinLabel(2,  "Pt");
  axis->SetBinLabel(3,  "Acceptance");
  axis->SetBinLabel(4,  "MCLabel");
  axis->SetBinLabel(5,  "BitMap");
  axis->SetBinLabel(6,  "HF cut");
  axis->SetBinLabel(7,  "Bit6");
  axis->SetBinLabel(8,  "NotHybridTrack");
  axis->SetBinLabel(9,  "MCFlag");
  axis->SetBinLabel(10, "MCGenerator");
  axis->SetBinLabel(11, "ChargeCut");
  axis->SetBinLabel(12, "MinDistanceTPCSectorEdge");
  axis->SetBinLabel(13, "Bit12");
  axis->SetBinLabel(14, "IsEMCal");
  axis->SetBinLabel(15, "Time");
  axis->SetBinLabel(16, "Energy");
  axis->SetBinLabel(17, "ExoticCut");
  axis->SetBinLabel(18, "Bit17");
  axis->SetBinLabel(19, "Area");
  axis->SetBinLabel(20, "AreaEmc");
  axis->SetBinLabel(21, "ZLeadingCh");
  axis->SetBinLabel(22, "ZLeadingEmc");
  axis->SetBinLabel(23, "NEF");
  axis->SetBinLabel(24, "MinLeadPt");
  axis->SetBinLabel(25, "MaxTrackPt");
  axis->SetBinLabel(26, "MaxClusterPt");
  axis->SetBinLabel(27, "Flavour");
  axis->SetBinLabel(28, "TagStatus");
  axis->SetBinLabel(29, "MinNConstituents");
  axis->SetBinLabel(30, "Bit29");
  axis->SetBinLabel(31, "Bit30");
  axis->SetBinLabel(32, "Bit31");
}

Double_t AliAnalysisTaskEmcalLight::GetParallelFraction(AliVParticle* part1, AliVParticle* part2)
{
  TVector3 vect1(part1->Px(), part1->Py(), part1->Pz());
  TVector3 vect2(part2->Px(), part2->Py(), part2->Pz());
  Double_t z = (vect1 * vect2) / (vect2 * vect2);
  return z;
}

Double_t AliAnalysisTaskEmcalLight::GetParallelFraction(const TVector3& vect1, AliVParticle* part2)
{
  TVector3 vect2(part2->Px(), part2->Py(), part2->Pz());
  Double_t z = (vect1 * vect2) / (vect2 * vect2);
  return z;
}

void AliAnalysisTaskEmcalLight::GetEtaPhiDiff(const AliVTrack *t, const AliVCluster *v, Double_t &phidiff, Double_t &etadiff)
{
  phidiff = 999;
  etadiff = 999;

  if (!t||!v) return;

  Double_t veta = t->GetTrackEtaOnEMCal();
  Double_t vphi = t->GetTrackPhiOnEMCal();

  Float_t pos[3] = {0};
  v->GetPosition(pos);
  TVector3 cpos(pos);
  Double_t ceta     = cpos.Eta();
  Double_t cphi     = cpos.Phi();
  etadiff=veta-ceta;
  phidiff=TVector2::Phi_mpi_pi(vphi-cphi);
}

Byte_t AliAnalysisTaskEmcalLight::GetTrackType(const AliVTrack *t)
{
  Byte_t ret = 0;
  if (t->TestBit(BIT(22)) && !t->TestBit(BIT(23)))
    ret = 1;
  else if (!t->TestBit(BIT(22)) && t->TestBit(BIT(23)))
    ret = 2;
  else if (t->TestBit(BIT(22)) && t->TestBit(BIT(23)))
    ret = 3;
  return ret;
}

Byte_t AliAnalysisTaskEmcalLight::GetTrackType(const AliAODTrack *aodTrack, UInt_t filterBit1, UInt_t filterBit2)
{

  Int_t res = 0;

  if (aodTrack->TestFilterBit(filterBit1)) {
    res = 0;
  }
  else if (aodTrack->TestFilterBit(filterBit2)) {
    if ((aodTrack->GetStatus()&AliVTrack::kITSrefit)!=0) {
      res = 1;
    }
    else {
      res = 2;
    }
  }
  else {
    res = 3;
  }

  return res;
}

AliAnalysisTaskEmcalLight::EBeamType_t AliAnalysisTaskEmcalLight::BeamTypeFromRunNumber(Int_t runnumber)
{
  EBeamType_t b = kpp;
  if ((runnumber >= 136833 && runnumber <= 139517) || // LHC10h Run-1 (Pb-Pb)
      (runnumber >= 167693 && runnumber <= 170593) || // LHC11h Run-1 (Pb-Pb)
      (runnumber >= 244824 && runnumber <= 246994) //|| // LHC15o Run-2 (Pb-Pb)
      //(runnumber >= 295581 && runnumber <= 297624)    // LHC18q+r Run-2 (Pb-Pb)
      ) 
  {     
    b = kAA;
  }
  else if ((runnumber > 188356 && runnumber <= 188503) ||  // LHC12g Run-1 (p-Pb pilot)
      (runnumber >= 195164 && runnumber <= 197388) ||      // LHC13b,c,d,e,f Run-1 (p-Pb)
      (runnumber >= 265077 && runnumber <= 267166)) {      // LHC16 Run-2 (p-Pb)
    b = kpA;
  }
  return b;
}

Bool_t AliAnalysisTaskEmcalLight::CheckMCOutliers()
{
  if (!fPythiaHeader || !fMCRejectFilter) return kTRUE;

  // Condition 1: Pythia jet / pT-hard > factor
  if (fPtHardAndJetPtFactor > 0.) {
    AliTLorentzVector jet;

    Int_t nTriggerJets =  fPythiaHeader->NTriggerJets();

    AliDebug(1,Form("Njets: %d, pT Hard %f",nTriggerJets, fPtHard));

    Float_t tmpjet[]={0,0,0,0};
    for (Int_t ijet = 0; ijet< nTriggerJets; ijet++) {
      fPythiaHeader->TriggerJet(ijet, tmpjet);

      jet.SetPxPyPzE(tmpjet[0],tmpjet[1],tmpjet[2],tmpjet[3]);

      AliDebug(1,Form("jet %d; pycell jet pT %f",ijet, jet.Pt()));

      //Compare jet pT and pt Hard
      if (jet.Pt() > fPtHardAndJetPtFactor * fPtHard) {
        AliInfo(Form("Reject jet event with : pT Hard %2.2f, pycell jet pT %2.2f, rejection factor %1.1f\n", fPtHard, jet.Pt(), fPtHardAndJetPtFactor));
        return kFALSE;
      }
    }
  }
  // end condition 1

  // Condition 2 : Reconstructed EMCal cluster pT / pT-hard > factor
  if (fPtHardAndClusterPtFactor > 0.) {
    AliClusterContainer* mccluscont = fClusterCollArray.begin()->second;
    if ((Bool_t)mccluscont) {
      for (auto cluster : mccluscont->all()) {// Not cuts applied ; use accept for cuts
        Float_t ecluster = cluster->E();

        if (ecluster > (fPtHardAndClusterPtFactor * fPtHard)) {
          AliInfo(Form("Reject : ecluster %2.2f, calo %d, factor %2.2f, ptHard %f",ecluster,cluster->GetType(),fPtHardAndClusterPtFactor,fPtHard));
          return kFALSE;
        }
      }
    }
  }
  // end condition 2

  // condition 3 : Reconstructed track pT / pT-hard >factor
  std::vector<AliMCParticleContainer *> mcpcont;
  for(auto cont : fParticleCollArray) {
    AliMCParticleContainer *mccont = dynamic_cast<AliMCParticleContainer *>(cont.second);
    if(mccont) mcpcont.push_back(mccont);
  }
  if (fPtHardAndTrackPtFactor > 0.) {
    AliMCParticleContainer* mcpartcont = *mcpcont.begin();
    if ((Bool_t)mcpartcont) {
      for (auto mctrack : mcpartcont->all()) {// Not cuts applied ; use accept for cuts
        Float_t trackpt = mctrack->Pt();
        if (trackpt > (fPtHardAndTrackPtFactor * fPtHard) ) {
          AliInfo(Form("Reject : track %2.2f, factor %2.2f, ptHard %f", trackpt, fPtHardAndTrackPtFactor, fPtHard));
          return kFALSE;
        }
      }
    }
  }
  // end condition 3

  return kTRUE;
}

Double_t AliAnalysisTaskEmcalLight::DeltaPhi(Double_t phia, Double_t phib, Double_t rangeMin, Double_t rangeMax)
{
  Double_t dphi = -999;
  const Double_t tpi = TMath::TwoPi();

  if (phia < 0)         phia += tpi;
  else if (phia > tpi) phia -= tpi;
  if (phib < 0)         phib += tpi;
  else if (phib > tpi) phib -= tpi;
  dphi = phib - phia;
  if (dphi < rangeMin)      dphi += tpi;
  else if (dphi > rangeMax) dphi -= tpi;

  return dphi;
}

void AliAnalysisTaskEmcalLight::GenerateFixedBinArray(int n, double min, double max, std::vector<double>& array, bool last)
{
  double binWidth = (max - min) / n;
  double v = min;
  if (last) n++;
  for (int i = 0; i < n; i++) {
    array.push_back(v);
    v += binWidth;
  }
}

std::vector<double> AliAnalysisTaskEmcalLight::GenerateFixedBinArray(int n, double min, double max, bool last)
{
  std::vector<double> array;
  GenerateFixedBinArray(n, min, max, array, last);
  return array;
}

void AliAnalysisTaskEmcalLight::GenerateLogFixedBinArray(int n, double min, double max, std::vector<double>& array, bool last)
{
  if (min <= 0 || max < min) {
    AliErrorClassStream() << "Cannot generate a log scale fixed-bin array with limits " << min << ", " << max << std::endl;
    return;
  }
  double binWidth = std::pow(max / min, 1.0 / n);
  double v = min;
  if (last) n++;
  for (int i = 0; i < n; i++) {
    array.push_back(v);
    v *= binWidth;
  }
}

std::vector<double> AliAnalysisTaskEmcalLight::GenerateLogFixedBinArray(int n, double min, double max, bool last)
{
  std::vector<double> array;
  GenerateLogFixedBinArray(n, min, max, array, last);
  return array;
}


TH1* AliAnalysisTaskEmcalLight::GetGeneralTH1(const char* name, bool warn)
{
  auto search = fHistograms.find(name);
  if (search != fHistograms.end()) {
    return search->second;
  }
  else {
    if (warn) AliErrorStream() << "Could not find histogram '" << name << "'" << std::endl;
    return nullptr;
  }
}

TH2* AliAnalysisTaskEmcalLight::GetGeneralTH2(const char* name, bool warn)
{
  return static_cast<TH2*>(GetGeneralTH1(name, warn));
}

TProfile* AliAnalysisTaskEmcalLight::GetGeneralTProfile(const char* name, bool warn)
{
  return static_cast<TProfile*>(GetGeneralTH1(name, warn));
}

void AliAnalysisTaskEmcalLight::SetIsPythia(Bool_t i)
{ 
  fIsPythia = i;
  if (fIsPythia) { 
    fIsMonteCarlo = kTRUE; 
    fMCEventHeaderName = "AliGenPythiaEventHeader"; 
  }
  else {
    if (fMCEventHeaderName == "AliGenPythiaEventHeader") {
      fMCEventHeaderName = "";
    }
  }
}

void AliAnalysisTaskEmcalLight::SetMCEventHeaderName(const char* name)
{ 
  TClass gen_header_class(name);
  if (gen_header_class.InheritsFrom("AliGenEventHeader")) {
    fMCEventHeaderName = name;
  }
  else {
    AliWarningStream() << "Class name '" << name << "' does not inherit from 'AliGenEventHeader'. Not setting it." << std::endl;
  }
}