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

#include <RVersion.h>
#include <TClonesArray.h>
#include <TList.h>
#include <TObject.h>
#include <TH1F.h>
#include <TProfile.h>
#include <TSystem.h>
#include <TFile.h>
#include <TChain.h>
#include <TKey.h>

#include "AliStack.h"
#include "AliAODEvent.h"
#include "AliAnalysisManager.h"
#include "AliCentrality.h"
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
#include "AliAODMCHeader.h"
#include "AliMCEvent.h"
#include "AliEMCALTriggerPatchInfo.h"

#include "AliMultSelection.h"

#include "AliAnalysisTaskEmcalLight.h"

Double_t AliAnalysisTaskEmcalLight::fgkEMCalDCalPhiDivide = 4.;

/// \cond CLASSIMP
ClassImp(AliAnalysisTaskEmcalLight)
/// \endcond

/**
 * Default constructor.
 */
AliAnalysisTaskEmcalLight::AliAnalysisTaskEmcalLight() :
  AliAnalysisTaskSE(),
  fForceBeamType(kNA),
  fGeneralHistograms(kFALSE),
  fCreateHisto(kTRUE),
  fNeedEmcalGeom(kTRUE),
  fCentBins(),
  fCentralityEstimation(kNewCentrality),
  fIsPythia(kFALSE),
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
  fLocalInitialized(kFALSE),
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
  fPythiaHeader(0),
  fPtHardBin(-1),
  fPtHard(0),
  fNTrials(0),
  fXsection(0),
  fOutput(0),
  fHistTrialsVsPtHardNoSel(0),
  fHistEventsVsPtHardNoSel(0),
  fHistXsectionVsPtHardNoSel(0),
  fHistTriggerClassesNoSel(0),
  fHistZVertexNoSel(0),
  fHistCentralityNoSel(0),
  fHistEventPlaneNoSel(0),
  fHistTrialsVsPtHard(0),
  fHistEventsVsPtHard(0),
  fHistXsectionVsPtHard(0),
  fHistTriggerClasses(0),
  fHistZVertex(0),
  fHistCentrality(0),
  fHistEventPlane(0),
  fHistEventCount(0),
  fHistEventRejection(0),
  fHistTrials(0),
  fHistEvents(0),
  fHistXsection(0)
{
  fVertex[0] = 0;
  fVertex[1] = 0;
  fVertex[2] = 0;
  fVertexSPD[0] = 0;
  fVertexSPD[1] = 0;
  fVertexSPD[2] = 0;
}

/**
 * Standard constructor. Should be used by the user.
 *
 * Note: This constructor also handles the general histograms. In
 * case the second parameter is true, then general histograms (see
 * UserCreateOutputObjects and FillHistograms) are created and filled
 * by the task, and a container is provided handling the user histograms.
 * @param[in] name Name of the task
 * @param[in] histo If true then general histograms are filled by the task
 */
AliAnalysisTaskEmcalLight::AliAnalysisTaskEmcalLight(const char *name, Bool_t histo) :
  AliAnalysisTaskSE(name),
  fForceBeamType(kNA),
  fGeneralHistograms(kFALSE),
  fCreateHisto(kTRUE),
  fNeedEmcalGeom(kTRUE),
  fCentBins(6),
  fCentralityEstimation(kNewCentrality),
  fIsPythia(kFALSE),
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
  fLocalInitialized(kFALSE),
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
  fPythiaHeader(0),
  fPtHardBin(-1),
  fPtHard(0),
  fNTrials(0),
  fXsection(0),
  fOutput(0),
  fHistTrialsVsPtHardNoSel(0),
  fHistEventsVsPtHardNoSel(0),
  fHistXsectionVsPtHardNoSel(0),
  fHistTriggerClassesNoSel(0),
  fHistZVertexNoSel(0),
  fHistCentralityNoSel(0),
  fHistEventPlaneNoSel(0),
  fHistTrialsVsPtHard(0),
  fHistEventsVsPtHard(0),
  fHistXsectionVsPtHard(0),
  fHistTriggerClasses(0),
  fHistZVertex(0),
  fHistCentrality(0),
  fHistEventPlane(0),
  fHistEventCount(0),
  fHistEventRejection(0),
  fHistTrials(0),
  fHistEvents(0),
  fHistXsection(0)
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

  if (fCreateHisto) DefineOutput(1, TList::Class());
}

/**
 * Destructor
 */
AliAnalysisTaskEmcalLight::~AliAnalysisTaskEmcalLight()
{
  for (auto cont_it : fParticleCollArray) delete cont_it.second;
  for (auto cont_it : fClusterCollArray) delete cont_it.second;
}

/**
 * Performing run-independent initialization. This consists of
 * - Determining data type (ESD/AOD)
 * - Creating general QA histograms
 *
 * Attention: Histograms are only created in case the task is
 * configured for this (second argument in the named constructor).
 * In this case the container fOuput is created which can be used
 * by the users to handle and store their histograms. In this case
 * the users must overwrite this function in their tasks and call
 * this function right at the beginning of their function.
 *
 * The general QA histograms monitor event related observables like
 * the z-position of the primary vertex before and after event selection,
 * the trigger classes selecting the event and the event rejection
 * reason, but also Monte-Carlo related observables like the cross
 * section, the number of trials and the \f$ p_{t} \f$-hard bin in
 * case of a corresponding production.
 */
void AliAnalysisTaskEmcalLight::UserCreateOutputObjects()
{
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (mgr) {
    AliVEventHandler *evhand = mgr->GetInputEventHandler();
    if (evhand) {
      if (evhand->InheritsFrom("AliESDInputHandler")) {
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
  fOutput = new TList();
  fOutput->SetOwner();

  if (fCentralityEstimation == kNoCentrality) fCentBins.clear();

  if (!fGeneralHistograms) return;

  if (fIsPythia) {
    fHistEventsVsPtHard = new TH1F("fHistEventsVsPtHard", "fHistEventsVsPtHard", 1000, 0, 1000);
    fHistEventsVsPtHard->GetXaxis()->SetTitle("#it{p}_{T,hard} (GeV/#it{c})");
    fHistEventsVsPtHard->GetYaxis()->SetTitle("events");
    fOutput->Add(fHistEventsVsPtHard);

    fHistTrialsVsPtHard = new TH1F("fHistTrialsVsPtHard", "fHistTrialsVsPtHard", 1000, 0, 1000);
    fHistTrialsVsPtHard->GetXaxis()->SetTitle("#it{p}_{T,hard} (GeV/#it{c})");
    fHistTrialsVsPtHard->GetYaxis()->SetTitle("trials");
    fOutput->Add(fHistTrialsVsPtHard);

    fHistXsectionVsPtHard = new TProfile("fHistXsectionVsPtHard", "fHistXsectionVsPtHard", 1000, 0, 1000);
    fHistXsectionVsPtHard->GetXaxis()->SetTitle("#it{p}_{T,hard} (GeV/#it{c})");
    fHistXsectionVsPtHard->GetYaxis()->SetTitle("xsection");
    fOutput->Add(fHistXsectionVsPtHard);

    fHistEventsVsPtHardNoSel = new TH1F("fHistEventsVsPtHardNoSel", "fHistEventsVsPtHardNoSel", 1000, 0, 1000);
    fHistEventsVsPtHardNoSel->GetXaxis()->SetTitle("#it{p}_{T,hard} (GeV/#it{c})");
    fHistEventsVsPtHardNoSel->GetYaxis()->SetTitle("events");
    fOutput->Add(fHistEventsVsPtHardNoSel);

    fHistTrialsVsPtHardNoSel = new TH1F("fHistTrialsVsPtHardNoSel", "fHistTrialsVsPtHardNoSel", 1000, 0, 1000);
    fHistTrialsVsPtHardNoSel->GetXaxis()->SetTitle("#it{p}_{T,hard} (GeV/#it{c})");
    fHistTrialsVsPtHardNoSel->GetYaxis()->SetTitle("trials");
    fOutput->Add(fHistTrialsVsPtHardNoSel);

    fHistXsectionVsPtHardNoSel = new TProfile("fHistXsectionVsPtHardNoSel", "fHistXsectionVsPtHardNoSel", 1000, 0, 1000);
    fHistXsectionVsPtHardNoSel->GetXaxis()->SetTitle("#it{p}_{T,hard} (GeV/#it{c})");
    fHistXsectionVsPtHardNoSel->GetYaxis()->SetTitle("xsection");
    fOutput->Add(fHistXsectionVsPtHardNoSel);

    fHistTrials = new TH1F("fHistTrials", "fHistTrials", 50, 0, 50);
    fHistTrials->GetXaxis()->SetTitle("#it{p}_{T,hard} bin");
    fHistTrials->GetYaxis()->SetTitle("trials");
    fOutput->Add(fHistTrials);

    fHistEvents = new TH1F("fHistEvents", "fHistEvents", 50, 0, 50);
    fHistEvents->GetXaxis()->SetTitle("#it{p}_{T,hard} bin");
    fHistEvents->GetYaxis()->SetTitle("total events");
    fOutput->Add(fHistEvents);

    fHistXsection = new TProfile("fHistXsection", "fHistXsection", 50, 0, 50);
    fHistXsection->GetXaxis()->SetTitle("#it{p}_{T,hard} bin");
    fHistXsection->GetYaxis()->SetTitle("xsection");
    fOutput->Add(fHistXsection);
  }

  fHistZVertex = new TH1F("fHistZVertex","Z vertex position", 60, -30, 30);
  fHistZVertex->GetXaxis()->SetTitle("V_{#it{z}}");
  fHistZVertex->GetYaxis()->SetTitle("counts");
  fOutput->Add(fHistZVertex);

  fHistZVertexNoSel = new TH1F("fHistZVertexNoSel","Z vertex position (no event selection)", 60, -30, 30);
  fHistZVertexNoSel->GetXaxis()->SetTitle("V_{#it{z}}");
  fHistZVertexNoSel->GetYaxis()->SetTitle("counts");
  fOutput->Add(fHistZVertexNoSel);

  if (fCentralityEstimation != kNoCentrality) {
    fHistCentrality = new TH1F("fHistCentrality","Event centrality distribution", 100, 0, 100);
    fHistCentrality->GetXaxis()->SetTitle("Centrality (%)");
    fHistCentrality->GetYaxis()->SetTitle("counts");
    fOutput->Add(fHistCentrality);

    fHistCentralityNoSel = new TH1F("fHistCentralityNoSel","Event centrality distribution (no event selection)", 100, 0, 100);
    fHistCentralityNoSel->GetXaxis()->SetTitle("Centrality (%)");
    fHistCentralityNoSel->GetYaxis()->SetTitle("counts");
    fOutput->Add(fHistCentralityNoSel);
  }

  if (fForceBeamType != kpp) {
    fHistEventPlane = new TH1F("fHistEventPlane","Event plane", 120, -TMath::Pi(), TMath::Pi());
    fHistEventPlane->GetXaxis()->SetTitle("event plane");
    fHistEventPlane->GetYaxis()->SetTitle("counts");
    fOutput->Add(fHistEventPlane);

    fHistEventPlaneNoSel = new TH1F("fHistEventPlaneNoSel","Event plane (no event selection)", 120, -TMath::Pi(), TMath::Pi());
    fHistEventPlaneNoSel->GetXaxis()->SetTitle("event plane");
    fHistEventPlaneNoSel->GetYaxis()->SetTitle("counts");
    fOutput->Add(fHistEventPlaneNoSel);
  }

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
  fHistEventRejection->GetXaxis()->SetBinLabel(14,"MCOutlier");
  fHistEventRejection->GetYaxis()->SetTitle("counts");
  fOutput->Add(fHistEventRejection);

  fHistTriggerClasses = new TH1F("fHistTriggerClasses","fHistTriggerClasses",3,0,3);
#if ROOT_VERSION_CODE < ROOT_VERSION(6,4,2)
  fHistTriggerClasses->SetBit(TH1::kCanRebin);
#else
  fHistTriggerClasses->SetCanExtend(TH1::kAllAxes);
#endif
  fOutput->Add(fHistTriggerClasses);

  fHistTriggerClassesNoSel = new TH1F("fHistTriggerClassesNoSel","fHistTriggerClassesNoSel",3,0,3);
#if ROOT_VERSION_CODE < ROOT_VERSION(6,4,2)
  fHistTriggerClassesNoSel->SetBit(TH1::kCanRebin);
#else
  fHistTriggerClassesNoSel->SetCanExtend(TH1::kAllAxes);
#endif
  fOutput->Add(fHistTriggerClassesNoSel);

  fHistEventCount = new TH1F("fHistEventCount","fHistEventCount",2,0,2);
  fHistEventCount->GetXaxis()->SetBinLabel(1,"Accepted");
  fHistEventCount->GetXaxis()->SetBinLabel(2,"Rejected");
  fHistEventCount->GetYaxis()->SetTitle("counts");
  fOutput->Add(fHistEventCount);

  PostData(1, fOutput);
}

/**
 * Filling general histrograms. Among the general histograms
 * that are filled only in case of running over MC productions
 * are
 * - \f$ p_{t} \f$-hard bin
 * - Cross section after event selection
 * - Number of trials after event selection
 * - Number of events after event selection
 * In any case the vertex distribution is filled as general
 * histograms. For heavy ion collisions also the centrality
 * distribution and the event plane distribution are filled.
 * @param[in] eventSelected flag that tells the method whether event selection has been performed already (a different set of histograms is filled)
 * @return Always true
 */
Bool_t AliAnalysisTaskEmcalLight::FillGeneralHistograms(Bool_t eventSelected)
{
  if (eventSelected) {
    if (fIsPythia) {
      fHistEventsVsPtHard->Fill(fPtHard, 1);
      fHistTrialsVsPtHard->Fill(fPtHard, fNTrials);
      fHistXsectionVsPtHard->Fill(fPtHard, fXsection);
    }

    fHistZVertex->Fill(fVertex[2]);

    if (fHistCentrality) fHistCentrality->Fill(fCent);
    if (fHistEventPlane) fHistEventPlane->Fill(fEPV0);


    for (auto fired_trg : fFiredTriggerClasses) fHistTriggerClasses->Fill(fired_trg.c_str(), 1);
  }
  else {
    if (fIsPythia) {
      fHistEventsVsPtHardNoSel->Fill(fPtHard, 1);
      fHistTrialsVsPtHardNoSel->Fill(fPtHard, fNTrials);
      fHistXsectionVsPtHardNoSel->Fill(fPtHard, fXsection);
    }

    fHistZVertexNoSel->Fill(fVertex[2]);

    if (fHistCentralityNoSel) fHistCentralityNoSel->Fill(fCent);
    if (fHistEventPlaneNoSel) fHistEventPlaneNoSel->Fill(fEPV0);

    for (auto fired_trg : fFiredTriggerClasses) fHistTriggerClassesNoSel->Fill(fired_trg.c_str(), 1);
  }

  return kTRUE;
}

/**
 * Event loop, called for each event. The function consists of three
 * steps:
 * -# Event selection
 * -# Running the user code
 * -# Filling general (QA) histograms
 * The event selection steps are documented in the function IsEventSelected.
 *
 * Users must not overwrite this function. Instead the virtual function Run
 * should be user and implemented by the user. The return value of the Run
 * function decides on whether general histograms are filled.
 *
 * In case the task is not yet initialized, which is the case for the first
 * event, the UserExec performs several basic initilization steps, documented
 * in the functions ExecOnce. Note that this is only done for the first event
 * and only for properties which need the presence of an input event.
 *
 * @param[in] option Not used
 */
void AliAnalysisTaskEmcalLight::UserExec(Option_t *option)
{
  if (!fLocalInitialized) ExecOnce();

  if (!fLocalInitialized) return;

  if (!RetrieveEventObjects()) return;

  Bool_t eventSelected = IsEventSelected();

  if (fGeneralHistograms && fCreateHisto) {
    if (eventSelected) {
      fHistEventCount->Fill("Accepted",1);
    }
    else {
      fHistEventCount->Fill("Rejected",1);
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

/**
 * Get the cross section and the trails either from pyxsec.root or from pysec_hists.root
 * Get the pt hard bin from the file path
 * This is to called in Notify and should provide the path to the AOD/ESD file
 * (Partially copied from AliAnalysisHelperJetTasks)
 * @param[in] currFile Name of the current ESD/AOD file
 * @param[out] fXsec Cross section calculated by PYTHIA
 * @param[out] fTrials Number of trials needed by PYTHIA
 * @param[out] pthard \f$ p_{t} \f$-hard bin, extracted from path name
 * @return True if parameters were obtained successfully, false otherwise
 */
Bool_t AliAnalysisTaskEmcalLight::PythiaInfoFromFile(const char* currFile, Float_t &fXsec, Float_t &fTrials, Int_t &pthard)
{

  TString file(currFile);
  fXsec = 0;
  fTrials = 1;

  if (file.Contains(".zip#")) {
    Ssiz_t pos1 = file.Index("root_archive",12,0,TString::kExact);
    Ssiz_t pos = file.Index("#",1,pos1,TString::kExact);
    Ssiz_t pos2 = file.Index(".root",5,TString::kExact);
    file.Replace(pos+1,pos2-pos1,"");
  } else {
    // not an archive take the basename....
    file.ReplaceAll(gSystem->BaseName(file.Data()),"");
  }
  AliDebug(1,Form("File name: %s",file.Data()));

  // Get the pt hard bin
  TString strPthard(file);

  strPthard.Remove(strPthard.Last('/'));
  strPthard.Remove(strPthard.Last('/'));
  if (strPthard.Contains("AOD")) strPthard.Remove(strPthard.Last('/'));
  strPthard.Remove(0,strPthard.Last('/')+1);
  if (strPthard.IsDec()) {
    pthard = strPthard.Atoi();
  }
  else {
    AliWarning(Form("Could not extract file number from path %s", strPthard.Data()));
    pthard = -1;
  }

  // problem that we cannot really test the existance of a file in a archive so we have to live with open error message from root
  TFile *fxsec = TFile::Open(Form("%s%s",file.Data(),"pyxsec.root"));

  if (!fxsec) {
    // next trial fetch the histgram file
    fxsec = TFile::Open(Form("%s%s",file.Data(),"pyxsec_hists.root"));
    if (!fxsec) {
      // not a severe condition but inciate that we have no information
      return kFALSE;
    } else {
      // find the tlist we want to be independtent of the name so use the Tkey
      TKey* key = static_cast<TKey*>(fxsec->GetListOfKeys()->At(0));
      if (!key) {
        fxsec->Close();
        return kFALSE;
      }
      TList *list = dynamic_cast<TList*>(key->ReadObj());
      if (!list) {
        fxsec->Close();
        return kFALSE;
      }
      fXsec = static_cast<TProfile*>(list->FindObject("h1Xsec"))->GetBinContent(1);
      fTrials = static_cast<TH1F*>(list->FindObject("h1Trials"))->GetBinContent(1);
      fxsec->Close();
    }
  } else { // no tree pyxsec.root
    TTree *xtree = static_cast<TTree*>(fxsec->Get("Xsection"));
    if (!xtree) {
      fxsec->Close();
      return kFALSE;
    }
    UInt_t   ntrials  = 0;
    Double_t  xsection  = 0;
    xtree->SetBranchAddress("xsection",&xsection);
    xtree->SetBranchAddress("ntrials",&ntrials);
    xtree->GetEntry(0);
    fTrials = ntrials;
    fXsec = xsection;
    fxsec->Close();
  }
  return kTRUE;
}

/**
 * Notifying the user that the input data file has
 * changed and performing steps needed to be done.
 *
 * This function is of relevance for analysis of
 * Monte-Carlo productions in \f$ p_{t} \f$-hard
 * bins as it reads the pythia cross section and
 * the number of trials from the file pyxsec.root
 * and fills the relevant distributions with
 * the values obtained.
 * @return False if the data tree or the data file
 * doesn't exist, true otherwise
 */
Bool_t AliAnalysisTaskEmcalLight::UserNotify()
{
  if (!fIsPythia || !fGeneralHistograms || !fCreateHisto)
    return kTRUE;

  TTree *tree = AliAnalysisManager::GetAnalysisManager()->GetTree();
  if (!tree) {
    AliError(Form("%s - UserNotify: No current tree!",GetName()));
    return kFALSE;
  }

  Float_t xsection    = 0;
  Float_t trials      = 0;
  Int_t   pthardbin   = 0;

  TFile *curfile = tree->GetCurrentFile();
  if (!curfile) {
    AliError(Form("%s - UserNotify: No current file!",GetName()));
    return kFALSE;
  }

  TChain *chain = dynamic_cast<TChain*>(tree);
  if (chain) tree = chain->GetTree();

  Int_t nevents = tree->GetEntriesFast();

  Bool_t res = PythiaInfoFromFile(curfile->GetName(), xsection, trials, pthardbin);

  fPtHardBin = pthardbin;

  if (!res) return kTRUE;

  fHistTrials->Fill(fPtHardBin, trials);
  fHistXsection->Fill(fPtHardBin, xsection);
  fHistEvents->Fill(fPtHardBin, nevents);

  return kTRUE;
}

/**
 * Perform steps needed to initialize the analysis.
 * This function relies on the presence of an input
 * event (ESD or AOD event). Consequently it is called
 * internally by UserExec for the first event.
 *
 * This function connects all containers attached to
 * this task to the corresponding arrays in the
 * input event. Furthermore it initializes the geometry.
 */
void AliAnalysisTaskEmcalLight::ExecOnce()
{
  if (!InputEvent()) {
    AliError(Form("%s: Could not retrieve event! Returning!", GetName()));
    return;
  }

  if (fNeedEmcalGeom) {
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

/**
 * Get beam type : pp-AA-pA
 * ESDs have it directly, AODs get it from hardcoded run number ranges
 * @return Beam type of the run.
 */
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

/**
 * Performing event selection. This contains
 * - Selection of the trigger class
 * - Selection according to the centrality class
 * - Selection of event with good vertex quality
 * - Selection of the event plane orientation
 * - Selection of the multiplicity (including
 *   above minimum \f$ p_{t} \f$ and tracks in the
 *   EMCAL acceptance
 *
 * Note that for the vertex selection both the usage
 * of the analysis util and the range of the z-position
 * of the primary vertex need to be specified.
 *
 * In case the event is rejected, a histogram
 * monitoring the rejeciton reason is filled with
 * the bin corresponding to the source of the rejection
 * of the current event.
 *
 * @return True if the event is selected.
 */
Bool_t AliAnalysisTaskEmcalLight::IsEventSelected()
{
  if (fTriggerSelectionBitMap != 0 && (fFiredTriggerBitMap & fTriggerSelectionBitMap) == 0) {
    if (fGeneralHistograms) fHistEventRejection->Fill("PhysSel",1);
    return kFALSE;
  }

  Bool_t acceptedTrgClassFound = kFALSE;
  if (fAcceptedTriggerClasses.size() > 0) {
    for (auto acc_trg : fAcceptedTriggerClasses) {
      for (auto fired_trg : fFiredTriggerClasses) {
        if (fired_trg.find(acc_trg) != std::string::npos) {
          acceptedTrgClassFound = kTRUE;
          break;
        }
      }
      if (acceptedTrgClassFound) break;
    }

    if (!acceptedTrgClassFound) {
      if (fGeneralHistograms) fHistEventRejection->Fill("Trg class (acc)",1);
      return kFALSE;
    }
  }

  if (fRejectedTriggerClasses.size() > 0) {
    for (auto rej_trg : fRejectedTriggerClasses) {
      for (auto fired_trg : fFiredTriggerClasses) {
        if (fired_trg.find(rej_trg) != std::string::npos) {
          if (fGeneralHistograms) fHistEventRejection->Fill("Trg class (rej)",1);
          return kFALSE;
        }
      }
    }
  }

  if (fMinCent < fMaxCent && fMaxCent > 0) {
    if (fCent < fMinCent || fCent > fMaxCent) {
      if (fGeneralHistograms) fHistEventRejection->Fill("Cent",1);
      return kFALSE;
    }
  }

  if (fNVertCont < fMinNVertCont) {
    if (fGeneralHistograms) fHistEventRejection->Fill("vertex contr.",1);
    return kFALSE;
  }

  if (fMinVz < fMaxVz) {
    if (fVertex[2] < fMinVz || fVertex[2] > fMaxVz) {
      if (fGeneralHistograms) fHistEventRejection->Fill("Vz",1);
      return kFALSE;
    }
  }

  if (fMaxVzDiff >= 0) {
    if (fNVertSPDCont > 0) {
      Double_t vzSPD = fVertexSPD[2];
      Double_t dvertex = TMath::Abs(fVertex[2] - vzSPD);
      //if difference larger than fZvertexDiff
      if (dvertex > fMaxVzDiff) {
        if (fGeneralHistograms) fHistEventRejection->Fill("VzSPD",1);
        return kFALSE;
      }
    }
  }

  if (fMinPtHard >= 0 && fPtHard < fMinPtHard)  {
    if (fGeneralHistograms) fHistEventRejection->Fill("SelPtHardBin",1);
    return kFALSE;
  }

  if (fMaxPtHard >= 0 && fPtHard >= fMaxPtHard)  {
    if (fGeneralHistograms) fHistEventRejection->Fill("SelPtHardBin",1);
    return kFALSE;
  }

  if (fPtHardBin == 0 && fMaxMinimumBiasPtHard >= 0 && fPtHard > fMaxMinimumBiasPtHard) {
    if (fGeneralHistograms) fHistEventRejection->Fill("SelPtHardBin",1);
    return kFALSE;
  }

  // Reject filter for MC data
  if (!CheckMCOutliers()) {
    if (fGeneralHistograms) fHistEventRejection->Fill("MCOutlier",1);
    return kFALSE;
  }

  return kTRUE;
}

/**
 * Read a TClonesArray from event. Attention: Both the
 * name of the array and the name of the object stored inside
 * must match.
 * @param[in] name Name of the array to be read in
 * @param[in] clname Name of the type of the objects stored in the array
 * @return Pointer to the TClonesArray (NULL if not found)
 */
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

/**
 * Retrieve objects from event.
 * @return
 */
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
      AliWarning(Form("%s: Could not retrieve centrality information! Assuming 99", GetName()));
    }
  }
  else if (fCentralityEstimation == kOldCentrality) {
    // Old centrality estimation (AliCentrality, works only on Run-1 PbPb and pPb)
    AliCentrality *aliCent = InputEvent()->GetCentrality();
    if (aliCent) {
      fCent = aliCent->GetCentralityPercentile(fCentEst.Data());
    }
    else {
      AliWarning(Form("%s: Could not retrieve centrality information! Assuming 99", GetName()));
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

    if (fPythiaHeader) {
      fPtHard = fPythiaHeader->GetPtHard();
      fXsection = fPythiaHeader->GetXsection();
      fNTrials = fPythiaHeader->Trials();
    }
  }

  for (auto cont_it : fParticleCollArray) cont_it.second->NextEvent();
  for (auto cont_it : fClusterCollArray) cont_it.second->NextEvent();

  return kTRUE;
}

/**
 * Create new particle container and attach it to the task. The name
 * provided to this function must match the name of the array attached
 * to the new container inside the input event.
 * @param[in] branchName Name of the array the container points to
 * @param[in] contName Name of the container points to (optional)
 * @return Pointer to the new particle container
 */
AliParticleContainer* AliAnalysisTaskEmcalLight::AddParticleContainer(std::string branchName, std::string contName)
{
  if (branchName.size() == 0) return 0;

  AliParticleContainer* cont = 0;

  if (branchName == "tracks" || branchName == "Tracks") cont = new AliTrackContainer(branchName.c_str());
  else if (branchName == "mcparticles") cont = new AliMCParticleContainer(branchName.c_str());
  else cont = new AliParticleContainer(branchName.c_str());

  if (contName.size() > 0) cont->SetName(contName.c_str());

  AdoptParticleContainer(cont);

  return cont;
}

/**
 * Create new cluster container and attach it to the task. The name
 * provided to this function must match the name of the array attached
 * to the new container inside the input event.
 * @param[in] branchName Name of the array the container points to
 * @param[in] contName Name of the container points to (optional)
 * @return Pointer to the new cluster container
 */
AliClusterContainer* AliAnalysisTaskEmcalLight::AddClusterContainer(std::string branchName, std::string contName)
{
  if (branchName.size() == 0) return 0;

  AliClusterContainer* cont = new AliClusterContainer(branchName.c_str());

  if (contName.size() > 0) cont->SetName(contName.c_str());

  AdoptClusterContainer(cont);

  return cont;
}

/**
 * Find particle container attached to this task according to its name
 * @param[in] name Name of the particle container
 * @return Particle container found under the given name
 */
AliParticleContainer* AliAnalysisTaskEmcalLight::GetParticleContainer(std::string name) const
{
  std::map<std::string, AliParticleContainer*>::const_iterator cont_it = fParticleCollArray.find(name);
  if (cont_it != fParticleCollArray.end()) return cont_it->second;
  else return nullptr;
}

/**
 * Find cluster container attached to this task according to its name
 * @param[in] name Name of the cluster container
 * @return Cluster container found under the given name
 */
AliClusterContainer* AliAnalysisTaskEmcalLight::GetClusterContainer(std::string name) const
{
  std::map<std::string, AliClusterContainer*>::const_iterator cont_it = fClusterCollArray.find(name);
  if (cont_it != fClusterCollArray.end()) return cont_it->second;
  else return nullptr;
}

/**
 * Add object to event
 * @param[in] obj Object to be added
 * @param[in] attempt If true don't handle error
 */
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

/**
 * Determines if a track is inside the EMCal acceptance, using \f$\eta\f$/\f$\phi\f$ at the vertex (no propagation).
 * Includes +/- edges. Useful to determine whether track propagation should be attempted.
 * @param[in] part Particle to check
 * @param[in] edges Size of the edges in \f$\phi\f$ excluded from the EMCAL acceptance
 * @return True if a particle is inside the EMCAL acceptance, false otherwise
 */
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

/**
 * Calculates the fraction of momentum z of part 1 w.r.t. part 2 in the direction of part 2.
 * @param[in] part1 Momentum vector for which the relative fraction is calculated
 * @param[in] part2 Reference momentum vector for the calculation
 * @return Relative fraction of momentum of particle 1 with respect to particle 2
 */
Double_t AliAnalysisTaskEmcalLight::GetParallelFraction(AliVParticle* part1, AliVParticle* part2)
{
  TVector3 vect1(part1->Px(), part1->Py(), part1->Pz());
  TVector3 vect2(part2->Px(), part2->Py(), part2->Pz());
  Double_t z = (vect1 * vect2) / (vect2 * vect2);
  return z;
}

/**
 * Calculates the fraction of momentum z of vect 1 w.r.t. part 2 in the direction of part 2.
 * @param[in] vect1 Momentum vector for which the relative fraction is calculated
 * @param[in] part2 Reference momentum vector for the calculation
 * @return Relative fraction of momentum of particle 1 with respect to particle 2
 */
Double_t AliAnalysisTaskEmcalLight::GetParallelFraction(const TVector3& vect1, AliVParticle* part2)
{
  TVector3 vect2(part2->Px(), part2->Py(), part2->Pz());
  Double_t z = (vect1 * vect2) / (vect2 * vect2);
  return z;
}

/**
 * Calculate \f$\phi\f$ and \f$\eta\f$ difference between a track (t) and a cluster (c). The
 * position of the track is obtained on the EMCAL surface
 * @param[in] t Track to check
 * @param[in] v Cluster to check
 * @param[out] phidiff Distance in \f$\phi\f$ between cluster and track
 * @param[out] etadiff Distance in \f$\eta\f$ between cluster and track
 */
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

/**
 * Get track type encoded from bits 20 and 21.
 * @param[in] t Track to check
 * @return
 */
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

/**
 * Return track type: 0 = filterBit1, 1 = filterBit2 && ITS, 2 = filterBit2 && !ITS.
 * Returns 3 if filterBit1 and filterBit2 do not test.
 * WARNING: only works with AOD tracks and AOD filter bits must be provided. Otherwise will always return 0.
 * @param aodTrack
 * @param filterBit1
 * @param filterBit2
 * @return
 */
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

/**
 * Determine the beam type based on hard-coded run ranges
 * \param runnumber run number
 * \return enumeration value corresponding to the beam type
 */
AliAnalysisTaskEmcalLight::EBeamType_t AliAnalysisTaskEmcalLight::BeamTypeFromRunNumber(Int_t runnumber)
{
  EBeamType_t b = kpp;
  if ((runnumber >= 136833 && runnumber <= 139517) || // LHC10h Run-1 (Pb-Pb)
      (runnumber >= 167693 && runnumber <= 170593) || // LHC11h Run-1 (Pb-Pb)
      (runnumber >= 244824 && runnumber <= 246994)) { // LHC15o Run-2 (Pb-Pb)
    b = kAA;
  }
  else if ((runnumber > 188356 && runnumber <= 188503) ||  // LHC12g Run-1 (p-Pb pilot)
      (runnumber >= 195164 && runnumber <= 197388) ||      // LHC13b,c,d,e,f Run-1 (p-Pb)
      (runnumber >= 265077 && runnumber <= 267166)) {      // LHC16 Run-2 (p-Pb)
    b = kpA;
  }
  return b;
}

/**
 * Filter the mc tails in pt-hard distributions
 * See https://twiki.cern.ch/twiki/bin/view/ALICE/JetMCProductionsCrossSections#How_to_reject_tails_in_the_pT_ha
 * @return kTRUE if it is not a MC outlier
 */
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

