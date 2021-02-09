/************************************************************************************
 * Copyright (C) 2017, Copyright Holders of the ALICE Collaboration                 *
 * All rights reserved.                                                             *
 *                                                                                  *
 * Redistribution and use in source and binary forms, with or without               *
 * modification, are permitted provided that the following conditions are met:      *
 *     * Redistributions of source code must retain the above copyright             *
 *       notice, this list of conditions and the following disclaimer.              *
 *     * Redistributions in binary form must reproduce the above copyright          *
 *       notice, this list of conditions and the following disclaimer in the        *
 *       documentation and/or other materials provided with the distribution.       *
 *     * Neither the name of the <organization> nor the                             *
 *       names of its contributors may be used to endorse or promote products       *
 *       derived from this software without specific prior written permission.      *
 *                                                                                  *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND  *
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED    *
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE           *
 * DISCLAIMED. IN NO EVENT SHALL ALICE COLLABORATION BE LIABLE FOR ANY              *
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES       *
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;     *
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND      *
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT       *
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS    *
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.                     *
 ************************************************************************************/
#include <RVersion.h>
#include <iostream>
#include <memory>

#include <TClonesArray.h>
#include <AliEmcalList.h>
#include <TObject.h>
#include <TObjString.h>
#include <TH1F.h>
#include <TProfile.h>
#include <TSystem.h>
#include <TFile.h>
#include <TChain.h>
#include <TKey.h>

#include "AliAnalysisTaskEmcal.h"
#include "AliAnalysisUtils.h"
#include "AliAODEvent.h"
#include "AliAODMCHeader.h"
#include "AliAODTrack.h"
#include "AliAnalysisManager.h"
#include "AliCentrality.h"
#include "AliEmcalDownscaleFactorsOCDB.h"
#include "AliEMCALGeometry.h"
#include "AliEmcalMCPartonInfo.h"
#include "AliEmcalPythiaInfo.h"
#include "AliEMCALTriggerPatchInfo.h"
#include "AliESDEvent.h"
#include "AliAODInputHandler.h"
#include "AliESDInputHandler.h"
#include "AliEventplane.h"
#include "AliGenPythiaEventHeader.h"
#include "AliGenHerwigEventHeader.h"
#include "AliGenHepMCEventHeader.h"
#include "AliInputEventHandler.h"
#include "AliLog.h"
#include "AliMCEvent.h"
#include "AliMCParticle.h"
#include "AliMultiInputEventHandler.h"
#include "AliMultSelection.h"
#include "AliStack.h"
#include "AliVCaloTrigger.h"
#include "AliVCluster.h"
#include "AliVEventHandler.h"
#include "AliVParticle.h"
#include "AliNanoAODHeader.h"
#include "AliAnalysisTaskEmcalEmbeddingHelper.h"

Double_t AliAnalysisTaskEmcal::fgkEMCalDCalPhiDivide = 4.;

/// \cond CLASSIMP
ClassImp(AliAnalysisTaskEmcal);
/// \endcond

AliAnalysisTaskEmcal::AliAnalysisTaskEmcal() : 
  AliAnalysisTaskSE("AliAnalysisTaskEmcal"),
  fPythiaInfoName(""),
  fNameMCPartonInfo(""),
  fForceBeamType(kNA),
  fGeneralHistograms(kFALSE),
  fLocalInitialized(kFALSE),
  fFileChanged(kTRUE),
  fCreateHisto(kTRUE),
  fCaloCellsName(),
  fCaloTriggersName(),
  fCaloTriggerPatchInfoName(),
  fMinCent(-999),
  fMaxCent(-999),
  fMinVz(-999),
  fMaxVz(999),
  fMinVertexContrib(1),
  fTrackPtCut(0),
  fMinNTrack(0),
  fZvertexDiff(0.5),
  fUseAliAnaUtils(kFALSE),
  fRejectPileup(kFALSE),
  fTklVsClusSPDCut(kFALSE),
  fOffTrigger(AliVEvent::kAny),
  fTrigClass(),
  fMinBiasRefTrigger("CINT7-B-NOPF-ALLNOTRD"),
  fTriggerTypeSel(kND),
  fNbins(250),
  fMinBinPt(0),
  fMaxBinPt(250),
  fMinPtTrackInEmcal(0),
  fEventPlaneVsEmcal(-1),
  fMinEventPlane(-1e6),
  fMaxEventPlane(1e6),
  fMinPtHard(-1e10),
  fMaxPtHard(1e10),
  fCentEst("V0M"),
  fIsEmbedded(kFALSE),
  fIsPythia(kFALSE),
  fIsHerwig(kFALSE),
  fIsHepMC(kFALSE),
  fGetPtHardBinFromName(kTRUE),
  fSelectPtHardBin(-999),
  fMinMCLabel(0),
  fMCLabelShift(0),
  fNcentBins(4),
  fNeedEmcalGeom(kTRUE),
  fParticleCollArray(),
  fClusterCollArray(),
  fTriggers(0),
  fEMCalTriggerMode(kOverlapWithLowThreshold),
  fUseNewCentralityEstimation(kFALSE),
  fGeneratePythiaInfoObject(kFALSE),
  fUsePtHardBinScaling(kFALSE),
  fUseXsecFromHeader(kFALSE),
  fMCRejectFilter(kFALSE),
  fCountDownscaleCorrectedEvents(kFALSE),
  fUseBuiltinEventSelection(kFALSE),
  fPtHardAndJetPtFactor(0.),
  fPtHardAndClusterPtFactor(0.),
  fPtHardAndTrackPtFactor(0.),
  fRunNumber(-1),
  fAliEventCuts(kFALSE),
  fAliAnalysisUtils(nullptr),
  fIsEsd(kFALSE),
  fGeom(nullptr),
  fTracks(nullptr),
  fCaloClusters(nullptr),
  fCaloCells(nullptr),
  fCaloTriggers(nullptr),
  fTriggerPatchInfo(nullptr),
  fCent(0),
  fCentBin(-1),
  fEPV0(-1.0),
  fEPV0A(-1.0),
  fEPV0C(-1.0),
  fNVertCont(0),
  fNVertSPDCont(0),
  fBeamType(kNA),
  fPythiaHeader(nullptr),
  fHerwigHeader(nullptr),
  fHepMCHeader(nullptr),
  fPtHard(0),
  fPtHardBin(0),
  fPtHardBinGlobal(-1),
  fPtHardInitialized(false),
  fDoCheckPtHardBin(true),
  fNPtHardBins(11),
  fPtHardBinning(),
  fNTrials(0),
  fXsection(0),
  fPythiaInfo(nullptr),
  fMCPartonInfo(nullptr),
  fOutput(nullptr),
  fHistEventCount(nullptr),
  fHistTrialsAfterSel(nullptr),
  fHistEventsAfterSel(nullptr),
  fHistXsectionAfterSel(nullptr),
  fHistTrials(nullptr),
  fHistEvents(nullptr),
  fHistXsection(nullptr),
  fHistPtHard(nullptr),
  fHistPtHardCorr(nullptr),
  fHistPtHardCorrGlobal(nullptr),
  fHistPtHardBinCorr(nullptr),
  fHistCentrality(nullptr),
  fHistZVertex(nullptr),
  fHistEventPlane(nullptr),
  fHistEventRejection(nullptr),
  fHistTriggerClasses(nullptr),
  fHistTriggerClassesCorr(nullptr)
{
  fVertex[0] = 0;
  fVertex[1] = 0;
  fVertex[2] = 0;
  fVertexSPD[0] = 0;
  fVertexSPD[1] = 0;
  fVertexSPD[2] = 0;

  fParticleCollArray.SetOwner(kTRUE);
  fClusterCollArray.SetOwner(kTRUE);
}

AliAnalysisTaskEmcal::AliAnalysisTaskEmcal(const char *name, Bool_t histo) : 
  AliAnalysisTaskSE(name),
  fPythiaInfoName(""),
  fNameMCPartonInfo(""),
  fForceBeamType(kNA),
  fGeneralHistograms(kFALSE),
  fLocalInitialized(kFALSE),
  fFileChanged(kFALSE),
  fCreateHisto(histo),
  fCaloCellsName(),
  fCaloTriggersName(),
  fCaloTriggerPatchInfoName(),
  fMinCent(-999),
  fMaxCent(-999),
  fMinVz(-999),
  fMaxVz(999),
  fMinVertexContrib(1),
  fTrackPtCut(0),
  fMinNTrack(0),
  fZvertexDiff(0.5),
  fUseAliAnaUtils(kFALSE),
  fRejectPileup(kFALSE),
  fTklVsClusSPDCut(kFALSE),
  fOffTrigger(AliVEvent::kAny),
  fTrigClass(),
  fMinBiasRefTrigger("CINT7-B-NOPF-ALLNOTRD"),
  fTriggerTypeSel(kND),
  fNbins(250),
  fMinBinPt(0),
  fMaxBinPt(250),
  fMinPtTrackInEmcal(0),
  fEventPlaneVsEmcal(-1),
  fMinEventPlane(-1e6),
  fMaxEventPlane(1e6),
  fMinPtHard(-1e10),
  fMaxPtHard(1e10),
  fCentEst("V0M"),
  fIsEmbedded(kFALSE),
  fIsPythia(kFALSE),
  fIsHerwig(kFALSE),
  fIsHepMC(kFALSE),
  fGetPtHardBinFromName(kTRUE),
  fSelectPtHardBin(-999),
  fMinMCLabel(0),
  fMCLabelShift(0),
  fNcentBins(4),
  fNeedEmcalGeom(kTRUE),
  fParticleCollArray(),
  fClusterCollArray(),
  fTriggers(0),
  fEMCalTriggerMode(kOverlapWithLowThreshold),
  fUseNewCentralityEstimation(kFALSE),
  fGeneratePythiaInfoObject(kFALSE),
  fUsePtHardBinScaling(kFALSE),
  fUseXsecFromHeader(kFALSE),
  fMCRejectFilter(kFALSE),
  fCountDownscaleCorrectedEvents(kFALSE),
  fUseBuiltinEventSelection(kFALSE),
  fPtHardAndJetPtFactor(0.),
  fPtHardAndClusterPtFactor(0.),
  fPtHardAndTrackPtFactor(0.),
  fRunNumber(-1),
  fAliEventCuts(kFALSE),
  fAliAnalysisUtils(nullptr),
  fIsEsd(kFALSE),
  fGeom(nullptr),
  fTracks(nullptr),
  fCaloClusters(nullptr),
  fCaloCells(nullptr),
  fCaloTriggers(nullptr),
  fTriggerPatchInfo(nullptr),
  fCent(0),
  fCentBin(-1),
  fEPV0(-1.0),
  fEPV0A(-1.0),
  fEPV0C(-1.0),
  fNVertCont(0),
  fNVertSPDCont(0),
  fBeamType(kNA),
  fPythiaHeader(nullptr),
  fHerwigHeader(nullptr),
  fHepMCHeader(nullptr),
  fPtHard(0),
  fPtHardBin(0),
  fPtHardBinGlobal(-1),
  fPtHardInitialized(false),
  fDoCheckPtHardBin(true),
  fNPtHardBins(11),
  fPtHardBinning(),
  fNTrials(0),
  fXsection(0),
  fPythiaInfo(nullptr),
  fMCPartonInfo(nullptr),
  fOutput(nullptr),
  fHistEventCount(nullptr),
  fHistTrialsAfterSel(nullptr),
  fHistEventsAfterSel(nullptr),
  fHistXsectionAfterSel(nullptr),
  fHistTrials(nullptr),
  fHistEvents(nullptr),
  fHistXsection(nullptr),
  fHistPtHard(nullptr),
  fHistPtHardCorr(nullptr),
  fHistPtHardCorrGlobal(nullptr),
  fHistPtHardBinCorr(nullptr),
  fHistCentrality(nullptr),
  fHistZVertex(nullptr),
  fHistEventPlane(nullptr),
  fHistEventRejection(nullptr),
  fHistTriggerClasses(nullptr),
  fHistTriggerClassesCorr(nullptr)
{
  fVertex[0] = 0;
  fVertex[1] = 0;
  fVertex[2] = 0;
  fVertexSPD[0] = 0;
  fVertexSPD[1] = 0;
  fVertexSPD[2] = 0;
  fParticleCollArray.SetOwner(kTRUE);
  fClusterCollArray.SetOwner(kTRUE);
  // Do not perform trigger selection in the AliEvent cuts but let the task do this before
  fAliEventCuts.OverrideAutomaticTriggerSelection(AliVEvent::kAny, true);

  if (fCreateHisto) {
    DefineOutput(1, AliEmcalList::Class());
  }
}

AliAnalysisTaskEmcal::~AliAnalysisTaskEmcal()
{
  if(fOutput) delete fOutput;
  if(fAliAnalysisUtils) delete fAliAnalysisUtils;
  if(fPythiaInfo) delete fPythiaInfo;
}

void AliAnalysisTaskEmcal::SetClusPtCut(Double_t cut, Int_t c)
{
  AliClusterContainer *cont = GetClusterContainer(c);
  if (cont) cont->SetClusPtCut(cut);
  else AliError(Form("%s in SetClusPtCut(...): container %d not found",GetName(),c));
}

void AliAnalysisTaskEmcal::SetClusTimeCut(Double_t min, Double_t max, Int_t c)
{
  AliClusterContainer *cont = GetClusterContainer(c);
  if (cont) cont->SetClusTimeCut(min,max);
  else AliError(Form("%s in SetClusTimeCut(...): container %d not found",GetName(),c));
}

void AliAnalysisTaskEmcal::SetTrackPtCut(Double_t cut, Int_t c)
{
  AliParticleContainer *cont = GetParticleContainer(c);
  if (cont) cont->SetParticlePtCut(cut);
  else AliError(Form("%s in SetTrackPtCut(...): container %d not found",GetName(),c));

  fTrackPtCut = cut;
}

void AliAnalysisTaskEmcal::SetTrackEtaLimits(Double_t min, Double_t max, Int_t c)
{
  AliParticleContainer *cont = GetParticleContainer(c);
  if (cont) cont->SetParticleEtaLimits(min,max);
  else AliError(Form("%s in SetTrackPtCut(...): container %d not found",GetName(),c));
}

void AliAnalysisTaskEmcal::SetTrackPhiLimits(Double_t min, Double_t max, Int_t c)
{
  AliParticleContainer *cont = GetParticleContainer(c);
  if (cont) cont->SetParticlePhiLimits(min,max);
  else AliError(Form("%s in SetTrackPhiLimits(...): container %d not found",GetName(),c));
}

void AliAnalysisTaskEmcal::UserCreateOutputObjects()
{
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (mgr) {
    AliVEventHandler *evhand = mgr->GetInputEventHandler();
    if (evhand) {
      if (evhand->InheritsFrom("AliESDInputHandler")) {
        fIsEsd = kTRUE;
      }
      else {
        fIsEsd = kFALSE;        
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
  fOutput = new AliEmcalList();
  fOutput->SetUseScaling(fUsePtHardBinScaling);
  fOutput->SetOwner();

  if (fForceBeamType == kpp)
    fNcentBins = 1;

  if (!fGeneralHistograms)
    return;

  if (fIsPythia || fIsHerwig || fIsHepMC) {
    fHistTrialsAfterSel = new TH1F("fHistTrialsAfterSel", "fHistTrialsAfterSel", fNPtHardBins, 0, fNPtHardBins);
    fHistTrialsAfterSel->GetXaxis()->SetTitle("p_{T} hard bin");
    fHistTrialsAfterSel->GetYaxis()->SetTitle("trials");
    fOutput->Add(fHistTrialsAfterSel);

    fHistEventsAfterSel = new TH1F("fHistEventsAfterSel", "fHistEventsAfterSel", fNPtHardBins, 0, fNPtHardBins);
    fHistEventsAfterSel->GetXaxis()->SetTitle("p_{T} hard bin");
    fHistEventsAfterSel->GetYaxis()->SetTitle("total events");
    fOutput->Add(fHistEventsAfterSel);

    fHistXsectionAfterSel = new TProfile("fHistXsectionAfterSel", "fHistXsectionAfterSel", fNPtHardBins, 0, fNPtHardBins);
    fHistXsectionAfterSel->GetXaxis()->SetTitle("p_{T} hard bin");
    fHistXsectionAfterSel->GetYaxis()->SetTitle("xsection");
    fOutput->Add(fHistXsectionAfterSel);

    fHistTrials = new TH1F("fHistTrials", "fHistTrials", fNPtHardBins, 0, fNPtHardBins);
    fHistTrials->GetXaxis()->SetTitle("p_{T} hard bin");
    fHistTrials->GetYaxis()->SetTitle("trials");
    fOutput->Add(fHistTrials);

    fHistEvents = new TH1F("fHistEvents", "fHistEvents", fNPtHardBins, 0, fNPtHardBins);
    fHistEvents->GetXaxis()->SetTitle("p_{T} hard bin");
    fHistEvents->GetYaxis()->SetTitle("total events");
    fOutput->Add(fHistEvents);

    fHistXsection = new TProfile("fHistXsection", "fHistXsection", fNPtHardBins, 0, fNPtHardBins);
    fHistXsection->GetXaxis()->SetTitle("p_{T} hard bin");
    fHistXsection->GetYaxis()->SetTitle("xsection");
    fOutput->Add(fHistXsection);

    // Set the bin labels
    Bool_t binningAvailable = false;
    if(fPtHardBinning.GetSize() > 0) {
      AliInfoStream() << "Using custom pt-hard binning" << std::endl;
      if(fPtHardBinning.GetSize() == fNPtHardBins + 1) binningAvailable = true;
      else AliErrorStream() << "Pt-hard binning (" << fPtHardBinning.GetSize() -1 << ") does not match the amount of bins (" << fNPtHardBins << ")" << std::endl;
    } else {
      // Check if we fall back to the default binning
      if(fNPtHardBins == 11) {
        AliInfoStream() << "11 pt-hard bins - fall back to default binning for bin labels" << std::endl;
        const Int_t kDefaultPtHardBinning[12] = {0,5,11,21,36,57, 84,117,152,191,234,1000000};
        fPtHardBinning.Set(12);
        for(Int_t ib = 0; ib < 12; ib++) fPtHardBinning[ib] = kDefaultPtHardBinning[ib];
        binningAvailable = true;
      } else {
        AliErrorStream() << "No default binning available for " << fNPtHardBins << " pt-hard bins - bin labels will not be set." << std::endl;
      }
    }

    if(binningAvailable){
      for (Int_t i = 0; i < fNPtHardBins; i++) {
        fHistTrialsAfterSel->GetXaxis()->SetBinLabel(i+1, Form("%d-%d",fPtHardBinning[i],fPtHardBinning[i+1]));
        fHistEventsAfterSel->GetXaxis()->SetBinLabel(i+1, Form("%d-%d",fPtHardBinning[i],fPtHardBinning[i+1]));
        fHistXsectionAfterSel->GetXaxis()->SetBinLabel(i+1, Form("%d-%d",fPtHardBinning[i],fPtHardBinning[i+1]));

        fHistTrials->GetXaxis()->SetBinLabel(i+1, Form("%d-%d",fPtHardBinning[i],fPtHardBinning[i+1]));
        fHistXsection->GetXaxis()->SetBinLabel(i+1, Form("%d-%d",fPtHardBinning[i],fPtHardBinning[i+1]));
        fHistEvents->GetXaxis()->SetBinLabel(i+1, Form("%d-%d",fPtHardBinning[i],fPtHardBinning[i+1]));
      }
    } else {
      AliErrorStream() << "No suitable binning available - skipping bin labels" << std::endl;
    }


    fHistPtHard = new TH1F("fHistPtHard", "fHistPtHard", fNbins*2, fMinBinPt, fMaxBinPt*4);
    fHistPtHard->GetXaxis()->SetTitle("p_{T,hard} (GeV/c)");
    fHistPtHard->GetYaxis()->SetTitle("counts");
    fOutput->Add(fHistPtHard);

    // Pt-hard correlation histograms (debugging)
    fHistPtHardCorr = new TH2F("fHistPtHardCorr", "Correlation between pt-hard value and binning", fNPtHardBins, -0.5, fNPtHardBins - 0.5, 1000, 0., 500.);
    fHistPtHardCorr->SetXTitle("p_{T,hard bin}");
    fHistPtHardCorr->SetYTitle("p_{T,hard value}");
    fOutput->Add(fHistPtHardCorr);

    fHistPtHardCorrGlobal = new TH2F("fHistPtHardCorrGlobal", "Correlation between global pt-hard value and binning", fNPtHardBins, -0.5, fNPtHardBins - 0.5, 1000, 0., 500.);
    fHistPtHardCorrGlobal->SetXTitle("p_{T,hard} bin_{global}");
    fHistPtHardCorrGlobal->SetYTitle("p_{T,hard} value");
    fOutput->Add(fHistPtHardCorrGlobal);

    fHistPtHardBinCorr = new TH2F("fHistPtHardBinCorr", "Correlation between global and local pt-hard bin", fNPtHardBins, -0.5, fNPtHardBins - 0.5, fNPtHardBins, -0.5, fNPtHardBins);
    fHistPtHardBinCorr->SetXTitle("p_{T,hard} bin_{local}");
    fHistPtHardBinCorr->SetYTitle("p_{T,hard} bin_{global}");
    fOutput->Add(fHistPtHardBinCorr);
  }

  fHistZVertex = new TH1F("fHistZVertex","Z vertex position", 60, -30, 30);
  fHistZVertex->GetXaxis()->SetTitle("z");
  fHistZVertex->GetYaxis()->SetTitle("counts");
  fOutput->Add(fHistZVertex);

  if (fForceBeamType != kpp) {
    fHistCentrality = new TH1F("fHistCentrality","Event centrality distribution", 200, 0, 100);
    fHistCentrality->GetXaxis()->SetTitle("Centrality (%)");
    fHistCentrality->GetYaxis()->SetTitle("counts");
    fOutput->Add(fHistCentrality);

    fHistEventPlane = new TH1F("fHistEventPlane","Event plane", 120, -TMath::Pi(), TMath::Pi());
    fHistEventPlane->GetXaxis()->SetTitle("event plane");
    fHistEventPlane->GetYaxis()->SetTitle("counts");
    fOutput->Add(fHistEventPlane);
  }

  if(fUseBuiltinEventSelection){
    fHistEventRejection = new TH1F("fHistEventRejection","Reasons to reject event",20,0,20);
#if ROOT_VERSION_CODE < ROOT_VERSION(6,4,2)
    fHistEventRejection->SetBit(TH1::kCanRebin);
#else
    fHistEventRejection->SetCanExtend(TH1::kAllAxes);
#endif
    fHistEventRejection->GetXaxis()->SetBinLabel(1,"PhysSel");
    fHistEventRejection->GetXaxis()->SetBinLabel(2,"trigger");
    fHistEventRejection->GetXaxis()->SetBinLabel(3,"trigTypeSel");
    fHistEventRejection->GetXaxis()->SetBinLabel(4,"Cent");
    fHistEventRejection->GetXaxis()->SetBinLabel(5,"vertex contr.");
    fHistEventRejection->GetXaxis()->SetBinLabel(6,"Vz");
    fHistEventRejection->GetXaxis()->SetBinLabel(7,"VzSPD");
    fHistEventRejection->GetXaxis()->SetBinLabel(8,"trackInEmcal");
    fHistEventRejection->GetXaxis()->SetBinLabel(9,"minNTrack");
    fHistEventRejection->GetXaxis()->SetBinLabel(10,"VtxSel2013pA");
    fHistEventRejection->GetXaxis()->SetBinLabel(11,"PileUp");
    fHistEventRejection->GetXaxis()->SetBinLabel(12,"EvtPlane");
    fHistEventRejection->GetXaxis()->SetBinLabel(13,"SelPtHardBin");
    fHistEventRejection->GetXaxis()->SetBinLabel(14,"Bkg evt");
    fHistEventRejection->GetYaxis()->SetTitle("counts");
    fOutput->Add(fHistEventRejection);
  }
  else {
    fAliEventCuts.AddQAplotsToList(fOutput);
  }

  fHistTriggerClasses = new TH1F("fHistTriggerClasses","fHistTriggerClasses",3,0,3);
#if ROOT_VERSION_CODE < ROOT_VERSION(6,4,2)
  fHistTriggerClasses->SetBit(TH1::kCanRebin);
#else
  fHistTriggerClasses->SetCanExtend(TH1::kAllAxes);
#endif
  fOutput->Add(fHistTriggerClasses);

  if(fCountDownscaleCorrectedEvents){
    fHistTriggerClassesCorr = new TH1F("fHistTriggerClassesCorr","fHistTriggerClassesCorr",3,0,3);
#if ROOT_VERSION_CODE < ROOT_VERSION(6,4,2)
    fHistTriggerClassesCorr->SetBit(TH1::kCanRebin);
#else
    fHistTriggerClassesCorr->SetCanExtend(TH1::kAllAxes);
#endif
    fOutput->Add(fHistTriggerClassesCorr);
  }

  fHistEventCount = new TH1F("fHistEventCount","fHistEventCount",2,0,2);
  fHistEventCount->GetXaxis()->SetBinLabel(1,"Accepted");
  fHistEventCount->GetXaxis()->SetBinLabel(2,"Rejected");
  fHistEventCount->GetYaxis()->SetTitle("counts");
  fOutput->Add(fHistEventCount);

  PostData(1, fOutput);
}

Bool_t AliAnalysisTaskEmcal::FillGeneralHistograms()
{
  if (fIsPythia || fIsHerwig || fIsHepMC) {
    // Protection: In case the pt-hard bin handling is not initialized we fall back to the
    // global pt-hard bin (usually 0) in order to aviod mismatch between histograms before
    // and after selection
    fHistEventsAfterSel->Fill(fPtHardInitialized ? fPtHardBinGlobal : fPtHardBin, 1);
    fHistTrialsAfterSel->Fill(fPtHardInitialized ? fPtHardBinGlobal : fPtHardBin, fNTrials);
    fHistXsectionAfterSel->Fill(fPtHardInitialized ? fPtHardBinGlobal : fPtHardBin, fXsection);
    fHistPtHard->Fill(fPtHard);
    if(fPtHardInitialized){
    	fHistPtHardCorr->Fill(fPtHardBin, fPtHard);
    	fHistPtHardCorrGlobal->Fill(fPtHardBinGlobal, fPtHard);
    	fHistPtHardBinCorr->Fill(fPtHardBin, fPtHardBinGlobal);
    }
  }


  fHistZVertex->Fill(fVertex[2]);

  if (fForceBeamType != kpp) {
    fHistCentrality->Fill(fCent);
    fHistEventPlane->Fill(fEPV0);
  }

  std::unique_ptr<TObjArray> triggerClasses(InputEvent()->GetFiredTriggerClasses().Tokenize(" "));
  TObjString* triggerClass(nullptr);
  for(auto trg : *triggerClasses){
    triggerClass = static_cast<TObjString*>(trg);
    fHistTriggerClasses->Fill(triggerClass->GetString(), 1);
  }

  if(fCountDownscaleCorrectedEvents){
    // downscale-corrected number of events are calculated based on the min. bias reference
    // Formula: N_corr = N_MB * d_Trg/d_{Min_Bias}
    if(InputEvent()->GetFiredTriggerClasses().Contains(fMinBiasRefTrigger)){
      PWG::EMCAL::AliEmcalDownscaleFactorsOCDB *downscalefactors = PWG::EMCAL::AliEmcalDownscaleFactorsOCDB::Instance();
      Double_t downscaleref = downscalefactors->GetDownscaleFactorForTriggerClass(fMinBiasRefTrigger);
      for(auto t : downscalefactors->GetTriggerClasses()){
        Double_t downscaletrg = downscalefactors->GetDownscaleFactorForTriggerClass(t);
        fHistTriggerClassesCorr->Fill(t, downscaletrg/downscaleref);
      }
    }
  }

  return kTRUE;
}

void AliAnalysisTaskEmcal::UserExec(Option_t *option)
{
  if (!fLocalInitialized){
    ExecOnce();
    UserExecOnce();
  }

  if (!fLocalInitialized)
    return;

  if(fFileChanged){
    FileChanged();
    fFileChanged = kFALSE;
  }

  if (!RetrieveEventObjects())
    return;

  if(InputEvent()->GetRunNumber() != fRunNumber){
    fRunNumber = InputEvent()->GetRunNumber();
    RunChanged(fRunNumber);
    if(fCountDownscaleCorrectedEvents) PWG::EMCAL::AliEmcalDownscaleFactorsOCDB::Instance()->SetRun(fRunNumber);
  }
 
  // Cut on the min. and max. pt-hard:
  // This is of relevance in particular when combining
  // min. bias events with pt-hard events by replacing low
  // pt-hard bins with min. bias events, where in min.
  // bias events with large pt-hard need to be removed in
  // order to not double count them when merging with the
  // pt-hard production in its bias-free region.
  // This should not be part of the normal event selection
  // as event counting for the cross section normalization
  // depends on it.
  if(fPtHard < fMinPtHard || fPtHard > fMaxPtHard) return;

  // Apply fallback for pythia cross section if needed
  if(fIsPythia && fUseXsecFromHeader && fPythiaHeader){
    AliDebugStream(1) << "Fallback to cross section from pythia header required" << std::endl;
    /*
    // Get the pthard bin
    Float_t pthard = fPythiaHeader->GetPtHard();
    int pthardbin = 0;
    if(fPtHardBinning.GetSize()){
      for(int ib = 0; ib < fNPtHardBins; ib++){
        if(pthard >= static_cast<Float_t>(fPtHardBinning[ib]) && pthard < static_cast<Float_t>(fPtHardBinning[ib+1])) {
          pthardbin = ib;
          break;
        }
      }
    }
    */
    fHistXsection->Fill(fPtHardBinGlobal, fPythiaHeader->GetXsection());
    fHistTrials->Fill(fPtHardBin);
    fHistEvents->Fill(fPtHardBin);
  }

  if(fIsHepMC && fHepMCHeader) {
    fHistXsection->Fill(fPtHardBinGlobal, fHepMCHeader->sigma_gen());
    fHistTrials->Fill(fPtHardBinGlobal, fHepMCHeader->ntrials());
    fHistEvents->Fill(fPtHardBinGlobal);
  }

  if (IsEventSelected()) {
    if (fGeneralHistograms) fHistEventCount->Fill("Accepted",1);
  }
  else {
    if (fGeneralHistograms) fHistEventCount->Fill("Rejected",1);
    return;
  }

  if (fGeneralHistograms && fCreateHisto) {
    if (!FillGeneralHistograms())
      return;
  }

  if (!Run())
    return;

  if (fCreateHisto) {
    if (!FillHistograms())
      return;
  }

  if (fCreateHisto && fOutput) {
    // information for this iteration of the UserExec in the container
    PostData(1, fOutput);
  }
}

Bool_t AliAnalysisTaskEmcal::AcceptCluster(AliVCluster *clus, Int_t c) const
{
  AliWarning("AliAnalysisTaskEmcal::AcceptCluster method is deprecated. Please use GetCusterContainer(c)->AcceptCluster(clus).");

  if (!clus) return kFALSE;

  AliClusterContainer *cont = GetClusterContainer(c);
  if (!cont) {
    AliError(Form("%s:Container %d not found",GetName(),c));
    return 0;
  }
  UInt_t rejectionReason = 0;
  return cont->AcceptCluster(clus, rejectionReason);
}

Bool_t AliAnalysisTaskEmcal::AcceptTrack(AliVParticle *track, Int_t c) const
{
  AliWarning("AliAnalysisTaskEmcal::AcceptTrack method is deprecated. Please use GetParticleContainer(c)->AcceptParticle(clus).");

  if (!track) return kFALSE;

  AliParticleContainer *cont = GetParticleContainer(c);
  if (!cont) {
    AliError(Form("%s:Container %d not found",GetName(),c));
    return 0;
  }

  UInt_t rejectionReason = 0;
  return cont->AcceptParticle(track, rejectionReason);
}

Int_t AliAnalysisTaskEmcal::ParsePtHardBinFromPath(const char *currentfile) {
  TString file(currentfile);
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

  Int_t pthard = -1;
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
  return pthard;
}


Bool_t AliAnalysisTaskEmcal::PythiaInfoFromFile(const char* currFile, Float_t &fXsec, Float_t &fTrials)
{
  fXsec = 0;
  fTrials = 1;

  TString file(currFile);
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

  // problem that we cannot really test the existance of a file in a archive so we have to live with open error message from root
  std::unique_ptr<TFile> fxsec(TFile::Open(Form("%s%s",file.Data(),"pyxsec.root")));

  if (!fxsec) {
    // next trial fetch the histgram file
    fxsec = std::unique_ptr<TFile>(TFile::Open(Form("%s%s",file.Data(),"pyxsec_hists.root")));
    if (!fxsec){
      AliErrorStream() << "Failed reading cross section from file " << file << std::endl;
      fUseXsecFromHeader = true;
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
        fUseXsecFromHeader = true;
      } else {
        // Cross section histogram filled - take it from there
        fXsec = xSecHist->GetBinContent(1);
        if(!fXsec) AliErrorStream() << GetName() << ": Cross section 0 for file " << file << std::endl;
        fUseXsecFromHeader = false;
      }
      fTrials  = ((TH1F*)list->FindObject("h1Trials"))->GetBinContent(1);
    }
  } else { // no tree pyxsec.root
    TTree *xtree = (TTree*)fxsec->Get("Xsection");
    if (!xtree) return kFALSE;
    UInt_t   ntrials  = 0;
    Double_t  xsection  = 0;
    xtree->SetBranchAddress("xsection",&xsection);
    xtree->SetBranchAddress("ntrials",&ntrials);
    xtree->GetEntry(0);
    fTrials = ntrials;
    fXsec = xsection;
  }
  return kTRUE;
}

Bool_t AliAnalysisTaskEmcal::UserNotify(){
  fPtHardInitialized = kFALSE;
  fFileChanged = kTRUE;
  return kTRUE;
}

Bool_t AliAnalysisTaskEmcal::FileChanged(){
  if (!(fIsPythia || fIsHepMC) || !fGeneralHistograms || !fCreateHisto)
    return kTRUE;

  // Handling of the pt-hard path common for pythia and HepMC pt-hard productions
  
  // Debugging:
  AliInfoStream() << "FileChanged called for run " << InputEvent()->GetRunNumber() << std::endl;

  TTree *tree = AliAnalysisManager::GetAnalysisManager()->GetTree();
  if (!tree) {
    AliErrorStream() << GetName() << " - FileChanged: No current tree!" << std::endl;
    return kFALSE;
  }
  TChain *chain = dynamic_cast<TChain*>(tree);
  if (chain) tree = chain->GetTree();

  TFile *curfile = tree->GetCurrentFile();
  if (!curfile) {
    AliErrorStream() << GetName() << " - FileChanged: No current file!" << std::endl;
    return kFALSE;
  }

  if(fGetPtHardBinFromName) {
    // Use the bin obtained from the path name
    fPtHardBinGlobal = ParsePtHardBinFromPath(curfile->GetName());
    fPtHardInitialized = kTRUE;
  } else {
    // Put everything in the first bin
    fPtHardBinGlobal = 0;
  }

  if ((fPtHardBinGlobal < 0) || (fPtHardBinGlobal > fNPtHardBins-1)){
    AliErrorStream() << GetName() << ": Invalid global pt-hard bin " << fPtHardBinGlobal << " detected" << std::endl;
    fPtHardBinGlobal = 0;
  }

  if(!fIsPythia) return kTRUE;

  Float_t xsection    = 0;
  Float_t trials      = 0;
  Int_t nevents = tree->GetEntriesFast();

  fUseXsecFromHeader = false;
  PythiaInfoFromFile(curfile->GetName(), xsection, trials);

  if(!fUseXsecFromHeader){
    AliDebugStream(1) << "Using cross section from file pyxsec.root" << std::endl;
    fHistXsection->Fill(fPtHardBinGlobal, xsection);
    fHistTrials->Fill(fPtHardBinGlobal, trials);
    fHistEvents->Fill(fPtHardBinGlobal, nevents);
  }

  return kTRUE;
}

void AliAnalysisTaskEmcal::LoadPythiaInfo(AliVEvent *event)
{
  if (!fPythiaInfoName.IsNull() && !fPythiaInfo) {
    fPythiaInfo = dynamic_cast<AliEmcalPythiaInfo*>(event->FindListObject(fPythiaInfoName));
    if (!fPythiaInfo) {
      AliError(Form("%s: Could not retrieve parton infos! %s!", GetName(), fPythiaInfoName.Data()));
      return;
    }
  }
}

void AliAnalysisTaskEmcal::LoadMCPartonInfo(AliVEvent *event)
{
  if (!fNameMCPartonInfo.IsNull() && !fMCPartonInfo) {
    fMCPartonInfo = dynamic_cast<PWG::EMCAL::AliEmcalMCPartonInfo*>(event->FindListObject(fNameMCPartonInfo));
    if (!fMCPartonInfo) {
      AliError(Form("%s: Could not retrieve parton infos! %s!", GetName(), fNameMCPartonInfo.Data()));
      return;
    }
  }
} 

void AliAnalysisTaskEmcal::ExecOnce()
{
  if (!InputEvent()) {
    AliError(Form("%s: Could not retrieve event! Returning!", GetName()));
    return;
  }

  LoadPythiaInfo(InputEvent());
  LoadMCPartonInfo(InputEvent());

  if (fNeedEmcalGeom) {
    fGeom = AliEMCALGeometry::GetInstanceFromRunNumber(InputEvent()->GetRunNumber());
    if (!fGeom) {
      AliFatal(Form("%s: Can not get EMCal geometry instance. If you do not need the EMCal geometry, disable it by setting task->SetNeedEmcalGeometry(kFALSE).", GetName()));
      return;
    }
  }

  if (fEventPlaneVsEmcal >= 0) {
    if (fGeom) {
      Double_t ep = (fGeom->GetArm1PhiMax() + fGeom->GetArm1PhiMin()) / 2 * TMath::DegToRad() + fEventPlaneVsEmcal - TMath::Pi();
      fMinEventPlane = ep - TMath::Pi() / 4;
      fMaxEventPlane = ep + TMath::Pi() / 4;
    }
    else {
      AliWarning("Could not set event plane limits because EMCal geometry was not loaded!");
    }
  }

  //Load all requested track branches - each container knows name already
  for (Int_t i =0; i<fParticleCollArray.GetEntriesFast(); i++) {
    AliParticleContainer *cont = static_cast<AliParticleContainer*>(fParticleCollArray.At(i));
    cont->SetArray(InputEvent());
  }

  if (fParticleCollArray.GetEntriesFast()>0) {
    fTracks = GetParticleArray(0);
    if (!fTracks) {
      AliError(Form("%s: Could not retrieve first track branch!", GetName()));
      return;
    }
  }

  //Load all requested cluster branches - each container knows name already
  for (Int_t i =0; i<fClusterCollArray.GetEntriesFast(); i++) {
    AliClusterContainer *cont = static_cast<AliClusterContainer*>(fClusterCollArray.At(i));
    cont->SetArray(InputEvent());
  }

  if (fClusterCollArray.GetEntriesFast()>0) {
    fCaloClusters = GetClusterArray(0);
    if (!fCaloClusters) {
      AliError(Form("%s: Could not retrieve first cluster branch!", GetName()));
      return;
    }
  }

  if (!fCaloCellsName.IsNull() && !fCaloCells) {
    TString objectname = fCaloCellsName;
    if(fCaloCellsName == "usedefault") {
      TString datatype; 
      if(fInputHandler->IsA() == AliAODInputHandler::Class()) {
        objectname = "emcalCells";
        datatype = "AOD"; 
      } else {
        objectname = "EMCALCells";
        datatype = "ESD";
      }
      AliInfoStream() << GetName() << ": [Cell container] usedefault: Using container " << objectname << " for data type " << datatype << std::endl;
    } else {
      AliInfoStream() << GetName() << ": [Cell container] user-defined: Using container " << objectname << std::endl;
    }
    fCaloCells =  dynamic_cast<AliVCaloCells*>(InputEvent()->FindListObject(objectname));
    if (!fCaloCells) {
      AliErrorStream() << GetName() << ": Could not retrieve cells " << objectname << "!" << std::endl; 
      return;
    }
  }

  if (!fCaloTriggersName.IsNull() && !fCaloTriggers) {
    TString objectname = fCaloTriggersName;
    if(fCaloTriggersName == "usedefault") {
      TString datatype; 
      if(fInputHandler->IsA() == AliAODInputHandler::Class()) {
        objectname = "emcalTrigger";
        datatype = "AOD"; 
      } else {
        objectname = "EMCALTrigger";
        datatype = "ESD";
      }
      AliInfoStream() << GetName() << ": [Trigger container] usedefault: Using container " << objectname << " for data type " << datatype << std::endl;
    } else {
      AliInfoStream() << GetName() << ": [Trigger container] user-defined: Using container " << objectname << std::endl;
    }
    fCaloTriggers =  dynamic_cast<AliVCaloTrigger*>(InputEvent()->FindListObject(objectname));
    if (!fCaloTriggers) {
      AliErrorStream() << GetName() <<": Could not retrieve calo triggers " << objectname << "!" << std::endl;
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

AliAnalysisTaskEmcal::BeamType AliAnalysisTaskEmcal::GetBeamType() const
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
    if ((runNumber >= 136833 && runNumber <= 139517) ||   // LHC10h
        (runNumber >= 167693 && runNumber <= 170593) ||   // LHC11h
        (runNumber >= 244824 && runNumber <= 246994) ||   // LHC15o
        (runNumber >= 295581 && runNumber <= 297624)) {   // LHC18p-q
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

ULong_t AliAnalysisTaskEmcal::GetTriggerList()
{
  if (!fTriggerPatchInfo)
    return 0;

  //number of patches in event
  Int_t nPatch = fTriggerPatchInfo->GetEntries();

  //loop over patches to define trigger type of event
  Int_t nG1 = 0;
  Int_t nG2 = 0;
  Int_t nJ1 = 0;
  Int_t nJ2 = 0;
  Int_t nL0 = 0;
  AliEMCALTriggerPatchInfo *patch;
  for (Int_t iPatch = 0; iPatch < nPatch; iPatch++) {
    patch = (AliEMCALTriggerPatchInfo*)fTriggerPatchInfo->At( iPatch );
    if (patch->IsGammaHigh()) nG1++;
    if (patch->IsGammaLow())  nG2++;
    if (patch->IsJetHigh()) nJ1++;
    if (patch->IsJetLow())  nJ2++;
    if (patch->IsLevel0())  nL0++;
  }

  AliDebug(2, "Patch summary: ");
  AliDebug(2, Form("Number of patches: %d", nPatch));
  AliDebug(2, Form("Jet:   low[%d], high[%d]" ,nJ2, nJ1));
  AliDebug(2, Form("Gamma: low[%d], high[%d]" ,nG2, nG1));

  ULong_t triggers(0);
  if (nL0>0) SETBIT(triggers, kL0);
  if (nG1>0) SETBIT(triggers, kG1);
  if (nG2>0) SETBIT(triggers, kG2);
  if (nJ1>0) SETBIT(triggers, kJ1);
  if (nJ2>0) SETBIT(triggers, kJ2);
  return triggers;
}

Bool_t AliAnalysisTaskEmcal::HasTriggerType(TriggerType trigger)
{
  //
  if(trigger==kND) {
    AliWarning(Form("%s: Requesting undefined trigger type!", GetName())); 
    return kFALSE;
  }
  //MV: removing this logic which as far as I can see doesn't make any sense
  // if(trigger & kND){
  //   return fTriggers == 0;
  // }
  return TESTBIT(fTriggers, trigger);
}

Bool_t AliAnalysisTaskEmcal::IsEventSelected(){
  if(fUseBuiltinEventSelection) return IsEventSelectedInternal();
  if(!IsTriggerSelected()) return false;
  if(!CheckMCOutliers()) return false;
  return fAliEventCuts.AcceptEvent(fInputEvent);
}

Bool_t AliAnalysisTaskEmcal::IsEventSelectedInternal()
{
  AliDebugStream(3) << "Using default event selection" << std::endl;
  if (fOffTrigger != AliVEvent::kAny) {
    UInt_t res = 0;
    const AliESDEvent *eev = dynamic_cast<const AliESDEvent*>(InputEvent());
    if (eev) {
      res = ((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected();
    } else {
      const AliAODEvent *aev = dynamic_cast<const AliAODEvent*>(InputEvent());
      if (aev) {
	      res = ((AliVAODHeader*)aev->GetHeader())->GetOfflineTrigger();
      }
    }
    if ((res & fOffTrigger) == 0) {
      if (fGeneralHistograms) fHistEventRejection->Fill("PhysSel",1);
      return kFALSE;
    }
  }

  if(!IsTriggerSelected()) {
    if (fGeneralHistograms) fHistEventRejection->Fill("trigger",1);
    return kFALSE;
  }

  if (fTriggerTypeSel != kND) {
    if (!HasTriggerType(fTriggerTypeSel)) {
      if (fGeneralHistograms) fHistEventRejection->Fill("trigTypeSel",1);
      return kFALSE;
    }
  }

  if ((fMinCent != -999) && (fMaxCent != -999)) {
    if (fCent<fMinCent || fCent>fMaxCent) {
      if (fGeneralHistograms) fHistEventRejection->Fill("Cent",1);
      return kFALSE;
    }
  }

  if (fUseAliAnaUtils) {
    if (!fAliAnalysisUtils)
      fAliAnalysisUtils = new AliAnalysisUtils();
    fAliAnalysisUtils->SetMinVtxContr(fMinVertexContrib);
    fAliAnalysisUtils->SetMaxVtxZ(999);
    if(fMinVz<-998.) fMinVz = -10.;
    if(fMaxVz>998.)  fMaxVz = 10.;

    if (!fAliAnalysisUtils->IsVertexSelected2013pA(InputEvent())) {
      if (fGeneralHistograms) fHistEventRejection->Fill("VtxSel2013pA",1);
      return kFALSE;
    }

    if (fRejectPileup && fAliAnalysisUtils->IsPileUpEvent(InputEvent())) {
      if (fGeneralHistograms) fHistEventRejection->Fill("PileUp",1);
      return kFALSE;
    }

    if(fTklVsClusSPDCut && fAliAnalysisUtils->IsSPDClusterVsTrackletBG(InputEvent())) {
      if (fGeneralHistograms) fHistEventRejection->Fill("Bkg evt",1);
      return kFALSE;
    }
  }

  if ((fMinVz > -998.) && (fMaxVz < 998.)) {
    if (fNVertCont == 0 ) {
      if (fGeneralHistograms) fHistEventRejection->Fill("vertex contr.",1);
      return kFALSE;
    }
    Double_t vz = fVertex[2];
    if (vz < fMinVz || vz > fMaxVz) {
      if (fGeneralHistograms) fHistEventRejection->Fill("Vz",1);
      return kFALSE;
    }

    if (fNVertSPDCont > 0 && fZvertexDiff < 999) {
      Double_t vzSPD = fVertexSPD[2];
      Double_t dvertex = TMath::Abs(vz-vzSPD);
      //if difference larger than fZvertexDiff
      if (dvertex > fZvertexDiff) {
        if (fGeneralHistograms) fHistEventRejection->Fill("VzSPD",1);
        return kFALSE;
      }
    }
  }

  if (fMinPtTrackInEmcal > 0 && fGeom) {
    Bool_t trackInEmcalOk = kFALSE;
    Int_t ntracks = GetNParticles(0);
    for (Int_t i = 0; i < ntracks; i++) {
      AliVParticle *track = GetAcceptParticleFromArray(i,0);
      if (!track)
        continue;

      Double_t phiMin = fGeom->GetArm1PhiMin() * TMath::DegToRad();
      Double_t phiMax = fGeom->GetArm1PhiMax() * TMath::DegToRad();
      Int_t runNumber = InputEvent()->GetRunNumber();
      if (runNumber>=177295 && runNumber<=197470) { //small SM masked in 2012 and 2013
        phiMin = 1.4;
        phiMax = TMath::Pi();
      }

      if (track->Eta() < fGeom->GetArm1EtaMin() || track->Eta() > fGeom->GetArm1EtaMax() || track->Phi() < phiMin || track->Phi() > phiMax)
        continue;
      if (track->Pt() > fMinPtTrackInEmcal) {
        trackInEmcalOk = kTRUE;
        break;
      }
    }
    if (!trackInEmcalOk) {
      if (fGeneralHistograms) fHistEventRejection->Fill("trackInEmcal",1);
      return kFALSE;
    }
  }

  if (fMinNTrack > 0) {
    Int_t nTracksAcc = 0;
    Int_t ntracks = GetNParticles(0);
    for (Int_t i = 0; i < ntracks; i++) {
      AliVParticle *track = GetAcceptParticleFromArray(i,0);
      if (!track)
        continue;
      if (track->Pt() > fTrackPtCut) {
        nTracksAcc++;
        if (nTracksAcc>=fMinNTrack)
          break;
      }
    }
    if (nTracksAcc<fMinNTrack) {
      if (fGeneralHistograms) fHistEventRejection->Fill("minNTrack",1);
      return kFALSE;
    }
  }

  if (!(fEPV0 > fMinEventPlane && fEPV0 <= fMaxEventPlane) &&
      !(fEPV0 + TMath::Pi() > fMinEventPlane && fEPV0 + TMath::Pi() <= fMaxEventPlane) &&
      !(fEPV0 - TMath::Pi() > fMinEventPlane && fEPV0 - TMath::Pi() <= fMaxEventPlane)) 
  {
    if (fGeneralHistograms) fHistEventRejection->Fill("EvtPlane",1);
    return kFALSE;
  }

  if (fSelectPtHardBin != -999 && fSelectPtHardBin != fPtHardBin)  {
    if (fGeneralHistograms) fHistEventRejection->Fill("SelPtHardBin",1);
    return kFALSE;
  }

  // Reject filter for MC data
  if (!CheckMCOutliers()) return kFALSE;

  return kTRUE;
}

Bool_t AliAnalysisTaskEmcal::IsTriggerSelected(){
  // Default implementation of trigger selection
  // same as previously (code moved from IsEventSelected
  // to trigger selection). Users should re-implement
  // this function in case they have certain needs, in
  // particular for EMCAL triggers
  AliDebugStream(3) << "Using default trigger selection" << std::endl;
  if (!fTrigClass.IsNull()) {
    TString fired = InputEvent()->GetFiredTriggerClasses();
    if (!fired.Contains("-B-")) return kFALSE;

    std::unique_ptr<TObjArray> arr(fTrigClass.Tokenize("|"));
    if (!arr) return kFALSE;
    Bool_t match = false;
    for (Int_t i=0;i<arr->GetEntriesFast();++i) {
      TObject *obj = arr->At(i);
      if (!obj) continue;

      //Check if requested trigger was fired
      TString objStr = obj->GetName();
      if(fEMCalTriggerMode == kOverlapWithLowThreshold &&
          (objStr.Contains("J1") || objStr.Contains("J2") || objStr.Contains("G1") || objStr.Contains("G2"))) {
        // This is relevant for EMCal triggers with 2 thresholds
        // If the kOverlapWithLowThreshold was requested than the overlap between the two triggers goes with the lower threshold trigger
        TString trigType1 = "J1";
        TString trigType2 = "J2";
        if(objStr.Contains("G")) {
          trigType1 = "G1";
          trigType2 = "G2";
        }
        if(objStr.Contains(trigType2) && fired.Contains(trigType2.Data())) { //requesting low threshold + overlap
          match = 1;
          break;
        } else if(objStr.Contains(trigType1) && fired.Contains(trigType1.Data()) && !fired.Contains(trigType2.Data())) { //high threshold only
          match = 1;
          break;
        }
      }
      else {
        // If this is not an EMCal trigger, or no particular treatment of EMCal triggers was requested,
        // simply check that the trigger was fired
        if (fired.Contains(obj->GetName())) {
          match = 1;
          break;
        }
      }
    }
    if (!match) return kFALSE;
  }
  return kTRUE;
}

Bool_t AliAnalysisTaskEmcal::CheckMCOutliers()
{
  if (!fPythiaHeader || !fMCRejectFilter) return kTRUE;
  AliDebugStream(2) << "Using custom outlier rejection" << std::endl;

  // Condition 1: Pythia jet / pT-hard > factor
  if (fPtHardAndJetPtFactor > 0.) {
    AliTLorentzVector jet;

    Int_t nTriggerJets =  fPythiaHeader->NTriggerJets();

    AliDebug(2,Form("Njets: %d, pT Hard %f",nTriggerJets, fPtHard));

    Float_t tmpjet[]={0,0,0,0};
    for (Int_t ijet = 0; ijet< nTriggerJets; ijet++) {
      fPythiaHeader->TriggerJet(ijet, tmpjet);

      jet.SetPxPyPzE(tmpjet[0],tmpjet[1],tmpjet[2],tmpjet[3]);

      AliDebug(2,Form("jet %d; pycell jet pT %f",ijet, jet.Pt()));

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
    AliClusterContainer* mccluscont = GetClusterContainer(0);
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
  if (fPtHardAndTrackPtFactor > 0.) {
    AliMCParticleContainer* mcpartcont = dynamic_cast<AliMCParticleContainer*>(GetParticleContainer(0));
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

TClonesArray *AliAnalysisTaskEmcal::GetArrayFromEvent(const char *name, const char *clname)
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

Bool_t AliAnalysisTaskEmcal::RetrieveEventObjects()
{
  fVertex[0] = 0;
  fVertex[1] = 0;
  fVertex[2] = 0;
  fNVertCont = 0;

  fVertexSPD[0] = 0;
  fVertexSPD[1] = 0;
  fVertexSPD[2] = 0;
  fNVertSPDCont = 0;

  if (fGeneratePythiaInfoObject && MCEvent()) {
    GeneratePythiaInfoObject(MCEvent());
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
  TObject * header = InputEvent()->GetHeader();
  if (fBeamType == kAA || fBeamType == kpA ) {
    if (fUseNewCentralityEstimation) {
      if (header->InheritsFrom("AliNanoAODStorage")){
        AliNanoAODHeader *nanoHead = (AliNanoAODHeader*)header;
        fCent=nanoHead->GetCentr(fCentEst.Data());
      }else{
        AliMultSelection *MultSelection = static_cast<AliMultSelection*>(InputEvent()->FindListObject("MultSelection"));
        if (MultSelection) {
          fCent = MultSelection->GetMultiplicityPercentile(fCentEst.Data());
        }
        else {
          AliWarning(Form("%s: Could not retrieve centrality information! Assuming 99", GetName()));
        }
      }
    }
    else { // old centrality estimation < 2015
      if (header->InheritsFrom("AliNanoAODStorage")){
        AliNanoAODHeader *nanoHead = (AliNanoAODHeader*)header;
        fCent=nanoHead->GetCentr(fCentEst.Data());
      }else{
        AliCentrality *aliCent = InputEvent()->GetCentrality();
        if (aliCent) {
          fCent = aliCent->GetCentralityPercentile(fCentEst.Data());
        }
        else {
          AliWarning(Form("%s: Could not retrieve centrality information! Assuming 99", GetName()));
        }
      }
    }

    if (fNcentBins==4) {
      if      (fCent >=  0 && fCent <   10) fCentBin = 0;
      else if (fCent >= 10 && fCent <   30) fCentBin = 1;
      else if (fCent >= 30 && fCent <   50) fCentBin = 2;
      else if (fCent >= 50 && fCent <= 100) fCentBin = 3;
      else {
        AliWarning(Form("%s: Negative centrality: %f. Assuming 99", GetName(), fCent));
        fCentBin = fNcentBins-1;
      }
    }
    else if (fNcentBins==5) {  // for PbPb 2015
      if      (fCent >=  0 && fCent <   10) fCentBin = 0;
      else if (fCent >= 10 && fCent <   30) fCentBin = 1;
      else if (fCent >= 30 && fCent <   50) fCentBin = 2;
      else if (fCent >= 50 && fCent <= 90) fCentBin = 3;
      else if (fCent > 90) {
        fCent = 99;
        fCentBin = 4;
      }
      else {
        AliWarning(Form("%s: Negative centrality: %f. Assuming 99", GetName(), fCent));
        fCentBin = fNcentBins-1;
      }
    }
    else {
      Double_t centWidth = (fMaxCent-fMinCent)/(Double_t)fNcentBins;
      if(centWidth>0.) {
        fCentBin = TMath::FloorNint(fCent/centWidth);
      }
      else {
        fCentBin = 0;
      }
      if (fCentBin>=fNcentBins) {
        AliWarning(Form("%s: fCentBin too large: cent = %f fCentBin = %d. Assuming 99", GetName(),fCent,fCentBin));
        fCentBin = fNcentBins-1;
      }
    }
    if (header->InheritsFrom("AliNanoAODStorage")){
        AliNanoAODHeader *nanoHead = (AliNanoAODHeader*)header;
        fEPV0=nanoHead->GetVar(nanoHead->GetVarIndex("cstEvPlaneV0"));
        fEPV0A=nanoHead->GetVar(nanoHead->GetVarIndex("cstEvPlaneV0A"));
        fEPV0C=nanoHead->GetVar(nanoHead->GetVarIndex("cstEvPlaneV0C"));
    }else{
    AliEventplane *aliEP = InputEvent()->GetEventplane();
    if (aliEP) {
      fEPV0  = aliEP->GetEventplane("V0" ,InputEvent());
      fEPV0A = aliEP->GetEventplane("V0A",InputEvent());
      fEPV0C = aliEP->GetEventplane("V0C",InputEvent());
    } else {
      AliWarning(Form("%s: Could not retrieve event plane information!", GetName()));
    }
    }
  }
  else {
    fCent = 99;
    fCentBin = 0;
  }

  if (fIsPythia) {
    if (MCEvent()) {
      fPythiaHeader = dynamic_cast<AliGenPythiaEventHeader*>(MCEvent()->GenEventHeader());
      if (!fPythiaHeader) {
        // Check if AOD
        AliAODMCHeader* aodMCH = dynamic_cast<AliAODMCHeader*>(InputEvent()->FindListObject(AliAODMCHeader::StdBranchName()));

        if (aodMCH) {
          for (UInt_t i = 0;i<aodMCH->GetNCocktailHeaders();i++) {
            fPythiaHeader = dynamic_cast<AliGenPythiaEventHeader*>(aodMCH->GetCocktailHeader(i));
            if (fPythiaHeader) break;
          }
        }
      }
    }
  }

  if (fPythiaHeader) {
    fPtHard = fPythiaHeader->GetPtHard();

    if(fPtHardBinning.GetSize()){
      // pt-hard binning defined for the corresponding dataset - automatically determine the bin
      for (fPtHardBin = 0; fPtHardBin < fNPtHardBins; fPtHardBin++) {
        if (fPtHard >= static_cast<float>(fPtHardBinning[fPtHardBin]) && fPtHard < static_cast<float>(fPtHardBinning[fPtHardBin+1]))
          break;
      }
    } else {
      // No pt-hard binning defined for the dataset - leaving the bin to 0
      fPtHardBin = 0;
    }

    if(fPtHardInitialized && fDoCheckPtHardBin){
      // do check only in case the global pt-hard bin is initialized
      if(fPtHardBin != fPtHardBinGlobal){
        AliErrorStream() << GetName() << ": Mismatch in pt-hard bin determination. Local: " << fPtHardBin << ", Global: " << fPtHardBinGlobal << std::endl;
      }
    }

    fXsection = fPythiaHeader->GetXsection();
    fNTrials = fPythiaHeader->Trials();
  }

  if (fIsHerwig) {
    if (MCEvent()) {
      fHerwigHeader = dynamic_cast<AliGenHerwigEventHeader*>(MCEvent()->GenEventHeader());
     
      if (!fHerwigHeader) {
        // Check if AOD
        AliAODMCHeader* aodMCH = dynamic_cast<AliAODMCHeader*>(InputEvent()->FindListObject(AliAODMCHeader::StdBranchName()));

        if (aodMCH) {
          for (UInt_t i = 0;i<aodMCH->GetNCocktailHeaders();i++) {
            fHerwigHeader = dynamic_cast<AliGenHerwigEventHeader*>(aodMCH->GetCocktailHeader(i));
            if (fHerwigHeader) break;
          }
        }
      }
    }
  }

  if (fHerwigHeader) {
    fPtHard = fHerwigHeader->GetPtHard();

    if(fPtHardBinning.GetSize()){
      // pt-hard binning defined for the corresponding dataset - automatically determine the bin
      for (fPtHardBin = 0; fPtHardBin < fNPtHardBins; fPtHardBin++) {
        if (fPtHard >= fPtHardBinning[fPtHardBin] && fPtHard < fPtHardBinning[fPtHardBin+1])
          break;
      }
    } else {
      // No pt-hard binning defined for the dataset - leaving the bin to 0
      fPtHardBin = 0;
    }
    if(fPtHardInitialized && fDoCheckPtHardBin){
      // do check only in case the global pt-hard bin is initialized
      if(fPtHardBin != fPtHardBinGlobal){
        AliErrorStream() << GetName() << ": Mismatch in pt-hard bin determination. Local: " << fPtHardBin << ", Global: " << fPtHardBinGlobal << std::endl;
      }
    }
    fXsection = fHerwigHeader->Weight();
    fNTrials = fHerwigHeader->Trials();
  }

  if (fIsHepMC) {
    if (MCEvent()) {
      fHepMCHeader = dynamic_cast<AliGenHepMCEventHeader*>(MCEvent()->GenEventHeader());
     
      if (!fHepMCHeader) {
        // Check if AOD
        AliAODMCHeader* aodMCH = dynamic_cast<AliAODMCHeader*>(InputEvent()->FindListObject(AliAODMCHeader::StdBranchName()));

        if (aodMCH) {
          for (UInt_t i = 0;i<aodMCH->GetNCocktailHeaders();i++) {
            fHepMCHeader = dynamic_cast<AliGenHepMCEventHeader*>(aodMCH->GetCocktailHeader(i));
            if (fHepMCHeader) break;
          }
        }
      }
    }
  }
  

  if (fHepMCHeader) {
    fPtHard = fHepMCHeader->pthard();

    if(fPtHardBinning.GetSize()){
      // pt-hard binning defined for the corresponding dataset - automatically determine the bin
      for (fPtHardBin = 0; fPtHardBin < fNPtHardBins; fPtHardBin++) {
        if (fPtHard >= fPtHardBinning[fPtHardBin] && fPtHard < fPtHardBinning[fPtHardBin+1])
          break;
      }
    } else {
      // No pt-hard binning defined for the dataset - leaving the bin to 0
      fPtHardBin = 0;
    }

    if(fPtHardInitialized && fDoCheckPtHardBin){
      // do check only in case the global pt-hard bin is initialized
      if(fPtHardBin != fPtHardBinGlobal){
        AliErrorStream() << GetName() << ": Mismatch in pt-hard bin determination. Local: " << fPtHardBin << ", Global: " << fPtHardBinGlobal << std::endl;
      }
    }
    fXsection = fHepMCHeader->sigma_gen();
    fNTrials = fHepMCHeader->ntrials();
  }


  fTriggers = GetTriggerList();

  AliEmcalContainer* cont = 0;

  TIter nextPartColl(&fParticleCollArray);
  while ((cont = static_cast<AliEmcalContainer*>(nextPartColl()))){
    cont->NextEvent(InputEvent());
  }

  TIter nextClusColl(&fClusterCollArray);
  while ((cont = static_cast<AliParticleContainer*>(nextClusColl()))){
    cont->NextEvent(InputEvent());
  }

  UserRetrieveEventObjects();

  return kTRUE;
}

AliMCParticleContainer* AliAnalysisTaskEmcal::AddMCParticleContainer(const char *n)
{
  if (TString(n).IsNull()) return 0;

  AliMCParticleContainer* cont = new AliMCParticleContainer(n);

  fParticleCollArray.Add(cont);

  return cont;
}

AliTrackContainer* AliAnalysisTaskEmcal::AddTrackContainer(const char *n)
{
  if (TString(n).IsNull()) return 0;

  AliTrackContainer* cont = new AliTrackContainer(n);

  fParticleCollArray.Add(cont);

  return cont;
}

AliParticleContainer* AliAnalysisTaskEmcal::AddParticleContainer(const char *n) 
{
  if (TString(n).IsNull()) return 0;

  AliParticleContainer* cont = new AliParticleContainer(n);

  fParticleCollArray.Add(cont);

  return cont;
}

AliClusterContainer* AliAnalysisTaskEmcal::AddClusterContainer(const char *n) 
{
  if (TString(n).IsNull()) return 0;

  AliClusterContainer* cont = new AliClusterContainer(n);

  fClusterCollArray.Add(cont);

  return cont;
}

AliParticleContainer* AliAnalysisTaskEmcal::GetParticleContainer(Int_t i) const 
{
  if (i<0 || i>fParticleCollArray.GetEntriesFast()) return 0;
  AliParticleContainer *cont = static_cast<AliParticleContainer*>(fParticleCollArray.At(i));
  return cont;
}

AliClusterContainer* AliAnalysisTaskEmcal::GetClusterContainer(Int_t i) const 
{
  if (i<0 || i>fClusterCollArray.GetEntriesFast()) return 0;
  AliClusterContainer *cont = static_cast<AliClusterContainer*>(fClusterCollArray.At(i));
  return cont;
}

AliParticleContainer* AliAnalysisTaskEmcal::GetParticleContainer(const char *name) const 
{
  AliParticleContainer *cont = static_cast<AliParticleContainer*>(fParticleCollArray.FindObject(name));
  return cont;
}

AliClusterContainer* AliAnalysisTaskEmcal::GetClusterContainer(const char *name) const 
{
  AliClusterContainer *cont = static_cast<AliClusterContainer*>(fClusterCollArray.FindObject(name));
  return cont;
}

TClonesArray* AliAnalysisTaskEmcal::GetParticleArray(Int_t i) const 
{
  AliParticleContainer *cont = GetParticleContainer(i);
  if (!cont) {
    AliError(Form("%s: Particle container %d not found",GetName(),i));
    return 0;
  }
  TString contName = cont->GetArrayName();
  return cont->GetArray();
}

TClonesArray* AliAnalysisTaskEmcal::GetClusterArray(Int_t i) const 
{
  AliClusterContainer *cont = GetClusterContainer(i);
  if (!cont) {
    AliError(Form("%s:Cluster container %d not found",GetName(),i));
    return 0;
  }
  return cont->GetArray();
}

AliVParticle* AliAnalysisTaskEmcal::GetAcceptParticleFromArray(Int_t p, Int_t c) const 
{

  AliParticleContainer *cont = GetParticleContainer(c);
  if (!cont) {
    AliError(Form("%s: Particle container %d not found",GetName(),c));
    return 0;
  }
  AliVParticle *vp = cont->GetAcceptParticle(p);

  return vp;
}

AliVCluster* AliAnalysisTaskEmcal::GetAcceptClusterFromArray(Int_t cl, Int_t c) const 
{
  AliClusterContainer *cont = GetClusterContainer(c);
  if (!cont) {
    AliError(Form("%s: Cluster container %d not found",GetName(),c));
    return 0;
  }
  AliVCluster *vc = cont->GetAcceptCluster(cl);

  return vc;
}

Int_t AliAnalysisTaskEmcal::GetNParticles(Int_t i) const 
{
  AliParticleContainer *cont = GetParticleContainer(i);
  if (!cont) {
    AliError(Form("%s: Particle container %d not found",GetName(),i));
    return 0;
  }
  return cont->GetNEntries();
}

Int_t AliAnalysisTaskEmcal::GetNClusters(Int_t i) const 
{
  AliClusterContainer *cont = GetClusterContainer(i);
  if (!cont) {
    AliError(Form("%s: Cluster container %d not found",GetName(),i));
    return 0;
  }
  return cont->GetNEntries();
}

AliEMCALTriggerPatchInfo* AliAnalysisTaskEmcal::GetMainTriggerPatch(TriggerCategory trigger, Bool_t doSimpleOffline)
{

  if (!fTriggerPatchInfo) {
    AliError(Form("%s: fTriggerPatchInfo not available",GetName()));
    return 0;
  }

  //number of patches in event
  Int_t nPatch = fTriggerPatchInfo->GetEntries();

  //extract main trigger patch(es)
  AliEMCALTriggerPatchInfo *patch(NULL), *selected(NULL);
  for (Int_t iPatch = 0; iPatch < nPatch; iPatch++) {

    patch = (AliEMCALTriggerPatchInfo*)fTriggerPatchInfo->At( iPatch );
    if (patch->IsMainTrigger()) {
      if(doSimpleOffline){
        if(patch->IsOfflineSimple()){
          switch(trigger){
          case kTriggerLevel0:
            // option not yet implemented in the trigger maker
            if(patch->IsLevel0()) selected = patch;
            break;
          case kTriggerLevel1Jet: 
            if(patch->IsJetHighSimple() || patch->IsJetLowSimple()){
              if(!selected) selected = patch;
              else if(patch->GetADCOfflineAmp() > selected->GetADCOfflineAmp()) selected = patch;
            }
            break;
          case kTriggerLevel1Gamma:
            if(patch->IsGammaHighSimple() || patch->IsGammaLowSimple()){
              if(!selected) selected = patch;
              else if(patch->GetADCOfflineAmp() > selected->GetADCOfflineAmp()) selected = patch;
            }
            break;
          default:   // Silence compiler warnings
            AliError("Untreated case: Main Patch is recalculated; should be in 'else' branch");
          };
        }
      } else {  // Not OfflineSimple
        switch(trigger){
        case kTriggerLevel0:
          if(patch->IsLevel0()) selected = patch;
          break;
        case kTriggerLevel1Jet:
          if(patch->IsJetHigh() || patch->IsJetLow()){
            if(!selected) selected = patch;
            else if (patch->GetADCAmp() > selected->GetADCAmp())
              selected = patch;
          }
          break;
        case kTriggerLevel1Gamma:
          if(patch->IsGammaHigh() || patch->IsGammaLow()){
            if(!selected) selected = patch;
            else if (patch->GetADCAmp() > selected->GetADCAmp())
              selected = patch;
          }
          break;
        default:
          AliError("Untreated case: Main Patch is recalculated; should be in 'else' branch");
        };
      }
    }
    else if ((trigger == kTriggerRecalcJet &&  patch->IsRecalcJet()) || 
        (trigger == kTriggerRecalcGamma && patch->IsRecalcGamma())) {  // recalculated patches
      if (doSimpleOffline && patch->IsOfflineSimple()) {
        if(!selected) selected = patch;
        else if (patch->GetADCOfflineAmp() > selected->GetADCOfflineAmp())  // this in fact should not be needed, but we have it in teh other branches as well, so keeping it for compleness
          selected = patch;
      }
      else if (!doSimpleOffline && !patch->IsOfflineSimple()) {
        if(!selected) selected = patch;
        else if (patch->GetADCAmp() > selected->GetADCAmp()) 
          selected = patch;
      }
    }
  }
  return selected;
}

void AliAnalysisTaskEmcal::AddObjectToEvent(TObject *obj, Bool_t attempt)
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

Bool_t AliAnalysisTaskEmcal::IsTrackInEmcalAcceptance(AliVParticle* part, Double_t edges) const
{

  if (!fGeom) {
    AliWarning(Form("%s - AliAnalysisTaskEmcal::IsTrackInEmcalAcceptance - Geometry is not available!", GetName()));
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

void AliAnalysisTaskEmcal::SetMCProductionType(MCProductionType_t prodtype) {
  if(prodtype == MCProductionType_t::kNoMC) return;
  switch(prodtype) {
    case MCProductionType_t::kMCPythiaMB: SetIsPythia(true); break;
    case MCProductionType_t::kMCPythiaPtHard: SetIsPythia(true); break;
    case MCProductionType_t::kMCHerwig6: SetIsHerwig(true); break;
    case MCProductionType_t::kMCHepMCMB: SetIsHepMC(true); break;
    case MCProductionType_t::kMCHepMCPtHard: SetIsHepMC(true); break;
    case MCProductionType_t::kNoMC: break;
  };
  // In case of min. bias production reduce to 1 pt-hard bin
  if(prodtype == MCProductionType_t::kMCPythiaMB || prodtype == MCProductionType_t::kMCHepMCMB) {
    SetNumberOfPtHardBins(1);
    TArrayI mbbinning(2);
    mbbinning[0] = 0;
    mbbinning[1] = 10000;
    SetUserPtHardBinning(mbbinning);
    SetGetPtHardBinFromPath(false);
  }
}


void AliAnalysisTaskEmcal::SetRejectionReasonLabels(TAxis* axis)
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

AliAnalysisTaskEmcal::MCProductionType_t AliAnalysisTaskEmcal::ConfigureMCDataset(const char *dataset) {  
  TString namedataset(dataset);
  namedataset.ToLower();
  PtHardBinning_t binningtype = PtHardBinning_t::kBinningUnknown;
  MCProductionType_t prodtype = MCProductionType_t::kNoMC;
  std::vector<TString> datasetsPthard20Pythia = {"lhc16c2", "lhc16h3", "lhc18b8", "lhc18f5", "lhc18g2", "lhc19a1", "lhc19d3", "lhc19f4", "lhc20g4"};
  std::vector<TString> datasetsPthard20HepMC = {"lhc20j3", "lhc20k1"};
  std::vector<TString> datasetsPthard13Pythia = {"lhc18i4a", "lhc18i4b2", "lhc19k3a", "lhc19k3b", "lhc19k3c"};
  std::vector<TString> datasetsPthard10Pythia = {"lhc12a15a", "lhc13b4"};
  std::vector<TString> datasetsPthard06Pythia = {"lhc17h6e2", "lhc17h6f2"};
  std::vector<TString> datasetsMBPythia = {
    "lhc17c3a", "lhc17h8a","lhc18l4a", "lhc18l4b",                                                                                // D-Mesons pp 13 TeV, 2016-18
    "lhc15h1", "lhc15h2",                                                                                                         // MB pp 8 TeV, 2012
    "lhc17f6", "lhc17f9", "lhc17d17", "lhc17f5", "lhc17d3", "lhc17e5", "lhc18f1", "lhc18d8", "lhc17d16", "lhc17d18",              // MB pp 13 TeV, 2016
    "lhc18d3", "lhc17h1", "lhc18c12", "lhc17k4", "lhc17h11", "lhc18c13", "lhc18a8", "lhc17l5", "lhc18a9", "lhc18a1",               // MB pp 13 TeV, 2017
    "lhc18g4", "lhc18g5", "lhc18g6", "lhc18g2", "lhc18h2", "lhc18h4", "lhc18j1", "lhc18j4", "lhc18k1", "lhc18k2", "lhc18k3",      // MB pp 13 TeV, 2018 
    "lhc18h1"                                                                                                                     // MB pp 13 TeV, 2018, low-B

  };

  bool foundDataset = false;
  for(auto dset : datasetsPthard20Pythia) {
    if(namedataset.Contains(dset)) {
      binningtype = PtHardBinning_t::kBinning20;
      prodtype = MCProductionType_t::kMCPythiaPtHard;
      foundDataset = true;
      break;
    }
  }

  if(!foundDataset) {
    for(auto dset : datasetsPthard13Pythia) {
      if(namedataset.Contains(dset)) {
        binningtype = PtHardBinning_t::kBinning13;
        prodtype = MCProductionType_t::kMCPythiaPtHard;
        foundDataset = true;
        break;
      }
    }
  }

  if(!foundDataset) {
    for(auto dset : datasetsPthard10Pythia) {
      if(namedataset.Contains(dset)) {
        binningtype = PtHardBinning_t::kBinning10;
        prodtype = MCProductionType_t::kMCPythiaPtHard;
        foundDataset = true;
        break;
      }
    }
  }

  if(!foundDataset) {
    for(auto dset : datasetsPthard06Pythia) {
      if(namedataset.Contains(dset)) {
        binningtype = PtHardBinning_t::kBinning06;
        prodtype = MCProductionType_t::kMCPythiaPtHard;
        foundDataset = true;
        break;
      }
    }
  }

  if(!foundDataset) {
    for(auto dset : datasetsPthard20HepMC) {
      if(namedataset.Contains(dset)) {
        binningtype = PtHardBinning_t::kBinning20;
        prodtype = MCProductionType_t::kMCHepMCPtHard;
        foundDataset = true;
        break;
      }
    }
  }

  if(!foundDataset) {
    for(auto dset : datasetsMBPythia){
      if(namedataset.Contains(dset)) {
        prodtype = MCProductionType_t::kMCPythiaMB;
        foundDataset = true;
        break;
      }
    }
  }

  SetMCProductionType(prodtype);
  if(binningtype != PtHardBinning_t::kBinningUnknown) {
    SetUsePtHardBinScaling(true);
    SetUserPtHardBinning(GetPtHardBinningForProd(binningtype));
  }
  return prodtype;
}

TArrayI AliAnalysisTaskEmcal::GetPtHardBinningForProd(PtHardBinning_t binningtype) {
  TArrayI binning;
  switch(binningtype) {
    case PtHardBinning_t::kBinning06: {
      const Int_t kNBinLimits = 8;
      binning.Set(kNBinLimits);
      const Int_t binlimits[] = {0, 5, 11, 21, 36, 57, 84, 1000000};
      memcpy(binning.GetArray(), binlimits, sizeof(int) * kNBinLimits);
      break; 
    }
    case PtHardBinning_t::kBinning10: {
      const Int_t kNBinLimits = 12;
      binning.Set(kNBinLimits);
      const Int_t binlimits[] = {0, 5, 11, 21, 36, 57, 84, 117, 152, 191, 234, 1000000};
      memcpy(binning.GetArray(), binlimits, sizeof(int) * kNBinLimits);
      break; 
    }
    case PtHardBinning_t::kBinning13: {
      const Int_t kNBinLimits = 15;
      binning.Set(kNBinLimits);
      const Int_t binlimits[] = {0, 5, 7, 9, 12, 16, 21, 28, 36, 45, 57, 70, 85, 100, 1000000};
      memcpy(binning.GetArray(), binlimits, sizeof(int) * kNBinLimits);
      break; 
    }
    case PtHardBinning_t::kBinning20: {
      const Int_t kNBinLimits = 22;
      binning.Set(kNBinLimits);
      const Int_t binlimits[] = {0, 5, 7, 9, 12, 16, 21, 28, 36, 45, 57, 70, 85, 99, 115, 132, 150, 169, 190, 212, 235, 1000000};
      memcpy(binning.GetArray(), binlimits, sizeof(int) * kNBinLimits);
      break;
    }
    default:
      AliErrorGeneralStream("AliAnalysisTaskEmcal::GetPtHardBinningForProd") << "Requested binning type not implemented" << std::endl;
  };

  return binning;
}

Double_t AliAnalysisTaskEmcal::GetParallelFraction(AliVParticle* part1, AliVParticle* part2)
{
  TVector3 vect1(part1->Px(), part1->Py(), part1->Pz());
  TVector3 vect2(part2->Px(), part2->Py(), part2->Pz());
  Double_t z = (vect1 * vect2) / (vect2 * vect2);
  return z;
}

Double_t AliAnalysisTaskEmcal::GetParallelFraction(const TVector3& vect1, AliVParticle* part2)
{
  TVector3 vect2(part2->Px(), part2->Py(), part2->Pz());
  Double_t z = (vect1 * vect2) / (vect2 * vect2);
  return z;
}

void AliAnalysisTaskEmcal::GetEtaPhiDiff(const AliVTrack *t, const AliVCluster *v, Double_t &phidiff, Double_t &etadiff)
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

Byte_t AliAnalysisTaskEmcal::GetTrackType(const AliVTrack *t)
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

Byte_t AliAnalysisTaskEmcal::GetTrackType(const AliAODTrack *aodTrack, UInt_t filterBit1, UInt_t filterBit2)
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

void AliAnalysisTaskEmcal::GeneratePythiaInfoObject(AliMCEvent* mcEvent)
{
  if (!fPythiaInfo) {
    fPythiaInfo = new AliEmcalPythiaInfo();
  }

  AliStack* stack = mcEvent->Stack();

  const Int_t nprim = stack->GetNprimary();
  // reject if partons are missing from stack for some reason
  if (nprim < 8) return;

  TParticle *part6 = stack->Particle(6);
  TParticle *part7 = stack->Particle(7);

  fPythiaInfo->SetPartonFlag6(TMath::Abs(part6->GetPdgCode()));
  fPythiaInfo->SetParton6(part6->Pt(), part6->Eta(), part6->Phi(), part6->GetMass());

  fPythiaInfo->SetPartonFlag7(TMath::Abs(part7->GetPdgCode()));
  fPythiaInfo->SetParton7(part7->Pt(), part7->Eta(), part7->Phi(), part7->GetMass());

  AliGenPythiaEventHeader *pythiaGenHeader = dynamic_cast<AliGenPythiaEventHeader*>(mcEvent->GenEventHeader());
  if(pythiaGenHeader){ 
    Float_t ptWeight=pythiaGenHeader->EventWeight(); 
    fPythiaInfo->SetPythiaEventWeight(ptWeight);}
}

AliAODInputHandler* AliAnalysisTaskEmcal::AddAODHandler()
{
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddAODHandler", "No analysis manager to connect to.");
    return NULL;
  }

  AliAODInputHandler* aodHandler = new AliAODInputHandler();

  AliVEventHandler *inputHandler=mgr->GetInputEventHandler();
  if (inputHandler && (inputHandler->IsA() == AliMultiInputEventHandler::Class())) {
    AliMultiInputEventHandler *multiInputHandler=(AliMultiInputEventHandler*)inputHandler;
    multiInputHandler->AddInputEventHandler(aodHandler);
  }
  else {
    if (!inputHandler) {
      mgr->SetInputEventHandler(aodHandler);
    }
    else {
      ::Error("AddAODHandler", "inputHandler is NOT null. AOD handler was NOT added !!!");
      return NULL;
    }
  }

  return aodHandler;
}

AliESDInputHandler* AliAnalysisTaskEmcal::AddESDHandler()
{
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddESDHandler", "No analysis manager to connect to.");
    return NULL;
  }

  AliESDInputHandler *esdHandler = new AliESDInputHandler();

  AliVEventHandler *inputHandler=mgr->GetInputEventHandler();
  if (inputHandler && (inputHandler->IsA() == AliMultiInputEventHandler::Class())) {
    AliMultiInputEventHandler *multiInputHandler=(AliMultiInputEventHandler*)inputHandler;
    multiInputHandler->AddInputEventHandler(esdHandler);
  }
  else {
    if (!inputHandler) {
      mgr->SetInputEventHandler(esdHandler);
    }
    else {
      ::Error("AddESDHandler", "inputHandler is NOT null. ESD handler was NOT added !!!");
      return NULL;
    }
  }

  return esdHandler;
}

