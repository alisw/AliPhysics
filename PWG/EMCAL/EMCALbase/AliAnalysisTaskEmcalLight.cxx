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
  fNcentBins(4),
  fUseNewCentralityEstimation(kFALSE),
  fIsPythia(kFALSE),
  fCaloCellsName(),
  fCaloTriggersName(),
  fCaloTriggerPatchInfoName(),
  fCentEst("V0M"),
  fParticleCollArray(),
  fClusterCollArray(),
  fTriggerSelectionBitMap(0),
  fMinCent(-999),
  fMaxCent(-999),
  fMinVz(-999),
  fMaxVz(999),
  fZvertexDiff(0.5),
  fMinPtTrack(0),
  fMinNTrack(0),
  fMinPtTrackInEmcal(0),
  fSelectPtHardBin(-999),
  fAcceptedTriggerClasses(),
  fRejectedTriggerClasses(),
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
  fPtHard(0),
  fPtHardBin(0),
  fNTrials(0),
  fXsection(0),
  fOutput(0),
  fHistEventCount(0),
  fHistTrialsAfterSel(0),
  fHistEventsAfterSel(0),
  fHistXsectionAfterSel(0),
  fHistTrials(0),
  fHistEvents(0),
  fHistXsection(0),
  fHistPtHard(0),
  fHistCentrality(0),
  fHistZVertex(0),
  fHistEventPlane(0),
  fHistEventRejection(0),
  fHistTriggerClasses(0)
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
  fNcentBins(4),
  fUseNewCentralityEstimation(kFALSE),
  fIsPythia(kFALSE),
  fCaloCellsName(),
  fCaloTriggersName(),
  fCaloTriggerPatchInfoName(),
  fCentEst("V0M"),
  fParticleCollArray(),
  fClusterCollArray(),
  fTriggerSelectionBitMap(0),
  fMinCent(-999),
  fMaxCent(-999),
  fMinVz(-999),
  fMaxVz(999),
  fZvertexDiff(0.5),
  fMinPtTrack(0),
  fMinNTrack(0),
  fMinPtTrackInEmcal(0),
  fSelectPtHardBin(-999),
  fAcceptedTriggerClasses(),
  fRejectedTriggerClasses(),
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
  fPtHard(0),
  fPtHardBin(0),
  fNTrials(0),
  fXsection(0),
  fOutput(0),
  fHistEventCount(0),
  fHistTrialsAfterSel(0),
  fHistEventsAfterSel(0),
  fHistXsectionAfterSel(0),
  fHistTrials(0),
  fHistEvents(0),
  fHistXsection(0),
  fHistPtHard(0),
  fHistCentrality(0),
  fHistZVertex(0),
  fHistEventPlane(0),
  fHistEventRejection(0),
  fHistTriggerClasses(0)
{
  fVertex[0] = 0;
  fVertex[1] = 0;
  fVertex[2] = 0;
  fVertexSPD[0] = 0;
  fVertexSPD[1] = 0;
  fVertexSPD[2] = 0;
  fParticleCollArray.SetOwner(kTRUE);
  fClusterCollArray.SetOwner(kTRUE);

  if (fCreateHisto) {
    DefineOutput(1, TList::Class()); 
  }
}

/**
 * Destructor
 */
AliAnalysisTaskEmcalLight::~AliAnalysisTaskEmcalLight()
{
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

  if (fForceBeamType == kpp)
    fNcentBins = 1;

  if (!fGeneralHistograms)
    return;

  if (fIsPythia) {
    fHistTrialsAfterSel = new TH1F("fHistTrialsAfterSel", "fHistTrialsAfterSel", 11, 0, 11);
    fHistTrialsAfterSel->GetXaxis()->SetTitle("p_{T} hard bin");
    fHistTrialsAfterSel->GetYaxis()->SetTitle("trials");
    fOutput->Add(fHistTrialsAfterSel);

    fHistEventsAfterSel = new TH1F("fHistEventsAfterSel", "fHistEventsAfterSel", 11, 0, 11);
    fHistEventsAfterSel->GetXaxis()->SetTitle("p_{T} hard bin");
    fHistEventsAfterSel->GetYaxis()->SetTitle("total events");
    fOutput->Add(fHistEventsAfterSel);

    fHistXsectionAfterSel = new TProfile("fHistXsectionAfterSel", "fHistXsectionAfterSel", 11, 0, 11);
    fHistXsectionAfterSel->GetXaxis()->SetTitle("p_{T} hard bin");
    fHistXsectionAfterSel->GetYaxis()->SetTitle("xsection");
    fOutput->Add(fHistXsectionAfterSel);

    fHistTrials = new TH1F("fHistTrials", "fHistTrials", 11, 0, 11);
    fHistTrials->GetXaxis()->SetTitle("p_{T} hard bin");
    fHistTrials->GetYaxis()->SetTitle("trials");
    fOutput->Add(fHistTrials);

    fHistEvents = new TH1F("fHistEvents", "fHistEvents", 11, 0, 11);
    fHistEvents->GetXaxis()->SetTitle("p_{T} hard bin");
    fHistEvents->GetYaxis()->SetTitle("total events");
    fOutput->Add(fHistEvents);

    fHistXsection = new TProfile("fHistXsection", "fHistXsection", 11, 0, 11);
    fHistXsection->GetXaxis()->SetTitle("p_{T} hard bin");
    fHistXsection->GetYaxis()->SetTitle("xsection");
    fOutput->Add(fHistXsection);

    const Int_t ptHardLo[11] = { 0, 5,11,21,36,57, 84,117,152,191,234};
    const Int_t ptHardHi[11] = { 5,11,21,36,57,84,117,152,191,234,1000000};

    for (Int_t i = 1; i < 12; i++) {
      fHistTrialsAfterSel->GetXaxis()->SetBinLabel(i, Form("%d-%d",ptHardLo[i-1],ptHardHi[i-1]));
      fHistEventsAfterSel->GetXaxis()->SetBinLabel(i, Form("%d-%d",ptHardLo[i-1],ptHardHi[i-1]));

      fHistTrials->GetXaxis()->SetBinLabel(i, Form("%d-%d",ptHardLo[i-1],ptHardHi[i-1]));
      fHistXsection->GetXaxis()->SetBinLabel(i, Form("%d-%d",ptHardLo[i-1],ptHardHi[i-1]));
      fHistEvents->GetXaxis()->SetBinLabel(i, Form("%d-%d",ptHardLo[i-1],ptHardHi[i-1]));
    }

    fHistPtHard = new TH1F("fHistPtHard", "fHistPtHard", 250, 0, 1000);
    fHistPtHard->GetXaxis()->SetTitle("p_{T,hard} (GeV/c)");
    fHistPtHard->GetYaxis()->SetTitle("counts");
    fOutput->Add(fHistPtHard);
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

  fHistTriggerClasses = new TH1F("fHistTriggerClasses","fHistTriggerClasses",3,0,3);
#if ROOT_VERSION_CODE < ROOT_VERSION(6,4,2)
  fHistTriggerClasses->SetBit(TH1::kCanRebin);
#else
  fHistTriggerClasses->SetCanExtend(TH1::kAllAxes);
#endif
  fOutput->Add(fHistTriggerClasses);

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
 * @return Always true
 */
Bool_t AliAnalysisTaskEmcalLight::FillGeneralHistograms()
{
  if (fIsPythia) {
    fHistEventsAfterSel->Fill(fPtHardBin, 1);
    fHistTrialsAfterSel->Fill(fPtHardBin, fNTrials);
    fHistXsectionAfterSel->Fill(fPtHardBin, fXsection);
    fHistPtHard->Fill(fPtHard);
  }

  fHistZVertex->Fill(fVertex[2]);

  if (fForceBeamType != kpp) {
    fHistCentrality->Fill(fCent);
    fHistEventPlane->Fill(fEPV0);
  }

  TObjArray* triggerClasses = InputEvent()->GetFiredTriggerClasses().Tokenize(" ");
  TIter next(triggerClasses);
  TObjString* triggerClass = 0;
  while ((triggerClass = static_cast<TObjString*>(next()))) {
    fHistTriggerClasses->Fill(triggerClass->GetString(), 1);
  }
  delete triggerClasses;
  triggerClasses = 0;

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
  if (!fLocalInitialized)
    ExecOnce();

  if (!fLocalInitialized)
    return;

  if (!RetrieveEventObjects())
    return;

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
  if (strPthard.IsDec()) 
    pthard = strPthard.Atoi();
  else 
    AliWarning(Form("Could not extract file number from path %s", strPthard.Data()));

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
      TKey* key = (TKey*)fxsec->GetListOfKeys()->At(0); 
      if (!key) {
        fxsec->Close();
        return kFALSE;
      }
      TList *list = dynamic_cast<TList*>(key->ReadObj());
      if (!list) {
        fxsec->Close();
        return kFALSE;
      }
      fXsec = ((TProfile*)list->FindObject("h1Xsec"))->GetBinContent(1);
      fTrials  = ((TH1F*)list->FindObject("h1Trials"))->GetBinContent(1);
      fxsec->Close();
    }
  } else { // no tree pyxsec.root
    TTree *xtree = (TTree*)fxsec->Get("Xsection");
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

  if (!res) return kTRUE;

  fHistTrials->Fill(pthardbin, trials);
  fHistXsection->Fill(pthardbin, xsection);
  fHistEvents->Fill(pthardbin, nevents);

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

  //Load all requested track branches - each container knows name already
  for (Int_t i =0; i<fParticleCollArray.GetEntriesFast(); i++) {
    AliParticleContainer *cont = static_cast<AliParticleContainer*>(fParticleCollArray.At(i));
    cont->SetArray(InputEvent());
  }

  //Load all requested cluster branches - each container knows name already
  for (Int_t i =0; i<fClusterCollArray.GetEntriesFast(); i++) {
    AliClusterContainer *cont = static_cast<AliClusterContainer*>(fClusterCollArray.At(i));
    cont->SetArray(InputEvent());
  }

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
    if ((runNumber >= 136851 && runNumber <= 139517) ||  // LHC10h
        (runNumber >= 166529 && runNumber <= 170593)) {  // LHC11h
      return kAA;
    } else if ((runNumber>=188365 && runNumber <= 188366) ||   // LHC12g
        (runNumber >= 195344 && runNumber <= 196608)) { // LHC13b-f
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
  if (fAcceptedTriggerClasses.GetEntriesFast() > 0) {
    TIter acceptedTrigger(&fAcceptedTriggerClasses);
    TObject* obj = 0;
    while ((obj = acceptedTrigger())) {
      if (fFiredTriggerClasses.Contains(obj->GetName())) {
        acceptedTrgClassFound = kTRUE;
        break;
      }
    }

    if (!acceptedTrgClassFound) {
      if (fGeneralHistograms) fHistEventRejection->Fill("Trg class (acc)",1);
      return kFALSE;
    }
  }

  if (fRejectedTriggerClasses.GetEntriesFast() > 0) {
    TIter rejectedTrigger(&fRejectedTriggerClasses);
    TObject* obj = 0;
    while ((obj = rejectedTrigger())) {
      if (fFiredTriggerClasses.Contains(obj->GetName())) {
        if (fGeneralHistograms) fHistEventRejection->Fill("Trg class (rej)",1);
        return kFALSE;
      }
    }
  }

  if ((fMinCent != -999) && (fMaxCent != -999)) {
    if (fCent < fMinCent || fCent > fMaxCent) {
      if (fGeneralHistograms) fHistEventRejection->Fill("Cent",1);
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

  if (fMinPtTrackInEmcal > 0 && fGeom && GetParticleContainer(0)) {
    Bool_t trackInEmcalOk = kFALSE;
    Int_t ntracks = GetParticleContainer(0)->GetNParticles();
    for (Int_t i = 0; i < ntracks; i++) {
      AliVParticle *track = GetParticleContainer(0)->GetAcceptParticle(i);
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
    Int_t ntracks = GetParticleContainer(0)->GetNParticles();
    for (Int_t i = 0; i < ntracks; i++) {
      AliVParticle *track = GetParticleContainer(0)->GetAcceptParticle(i);
      if (!track)
        continue;
      if (track->Pt() > fMinPtTrack) {
        nTracksAcc++;
        if (nTracksAcc >= fMinNTrack)
          break;
      }
    }
    if (nTracksAcc<fMinNTrack) {
      if (fGeneralHistograms) fHistEventRejection->Fill("minNTrack",1);
      return kFALSE;
    }
  }

  if (fSelectPtHardBin != -999 && fSelectPtHardBin != fPtHardBin)  {
    if (fGeneralHistograms) fHistEventRejection->Fill("SelPtHardBin",1);
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

  fFiredTriggerClasses = InputEvent()->GetFiredTriggerClasses();

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

  if (fBeamType == kAA || fBeamType == kpA ) {
    if (fUseNewCentralityEstimation) {
      AliMultSelection *MultSelection = static_cast<AliMultSelection*>(InputEvent()->FindListObject("MultSelection"));
      if (MultSelection) {
        fCent = MultSelection->GetMultiplicityPercentile(fCentEst.Data());
      }
      else {
        AliWarning(Form("%s: Could not retrieve centrality information! Assuming 99", GetName()));
      }
    }
    else { // old centrality estimation < 2015
      AliCentrality *aliCent = InputEvent()->GetCentrality();
      if (aliCent) {
        fCent = aliCent->GetCentralityPercentile(fCentEst.Data());
      }
      else {
        AliWarning(Form("%s: Could not retrieve centrality information! Assuming 99", GetName()));
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

    AliEventplane *aliEP = InputEvent()->GetEventplane();
    if (aliEP) {
      fEPV0  = aliEP->GetEventplane("V0" ,InputEvent());
      fEPV0A = aliEP->GetEventplane("V0A",InputEvent());
      fEPV0C = aliEP->GetEventplane("V0C",InputEvent());
    } else {
      AliWarning(Form("%s: Could not retrieve event plane information!", GetName()));
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

    if (fPythiaHeader) {
      fPtHard = fPythiaHeader->GetPtHard();

      const Int_t ptHardLo[11] = { 0, 5,11,21,36,57, 84,117,152,191,234};
      const Int_t ptHardHi[11] = { 5,11,21,36,57,84,117,152,191,234,1000000};
      for (fPtHardBin = 0; fPtHardBin < 11; fPtHardBin++) {
        if (fPtHard >= ptHardLo[fPtHardBin] && fPtHard < ptHardHi[fPtHardBin])
          break;
      }

      fXsection = fPythiaHeader->GetXsection();
      fNTrials = fPythiaHeader->Trials();
    }
  }

  AliEmcalContainer* cont = 0;

  TIter nextPartColl(&fParticleCollArray);
  while ((cont = static_cast<AliEmcalContainer*>(nextPartColl()))) cont->NextEvent();

  TIter nextClusColl(&fClusterCollArray);
  while ((cont = static_cast<AliParticleContainer*>(nextClusColl()))) cont->NextEvent();

  return kTRUE;
}

/**
 * Create new container for MC particles and attach it to the task. The name
 * provided to this function must match the name of the array attached
 * to the new container inside the input event.
 * @param[in] n Name of the container and the array the container points to
 * @return Pointer to the new container for MC particles
 */
AliMCParticleContainer* AliAnalysisTaskEmcalLight::AddMCParticleContainer(const char *n)
{
  if (TString(n).IsNull()) return 0;

  AliMCParticleContainer* cont = new AliMCParticleContainer(n);

  fParticleCollArray.Add(cont);

  return cont;
}

/**
 * Create new track container and attach it to the task. The name
 * provided to this function must match the name of the array attached
 * to the new container inside the input event.
 * @param[in] n Name of the container and the array the container points to
 * @return Pointer to the new track container
 */
AliTrackContainer* AliAnalysisTaskEmcalLight::AddTrackContainer(const char *n)
{
  if (TString(n).IsNull()) return 0;

  AliTrackContainer* cont = new AliTrackContainer(n);

  fParticleCollArray.Add(cont);

  return cont;
}

/**
 * Create new particle container and attach it to the task. The name
 * provided to this function must match the name of the array attached
 * to the new container inside the input event.
 * @param[in] n Name of the container and the array the container points to
 * @return Pointer to the new particle container
 */
AliParticleContainer* AliAnalysisTaskEmcalLight::AddParticleContainer(const char *n)
{
  if (TString(n).IsNull()) return 0;

  AliParticleContainer* cont = new AliParticleContainer(n);

  fParticleCollArray.Add(cont);

  return cont;
}

/**
 * Create new cluster container and attach it to the task. The name
 * provided to this function must match the name of the array attached
 * to the new container inside the input event.
 * @param[in] n Name of the container and the array the container points to
 * @return Pointer to the new cluster container
 */
AliClusterContainer* AliAnalysisTaskEmcalLight::AddClusterContainer(const char *n)
{
  if (TString(n).IsNull()) return 0;

  AliClusterContainer* cont = new AliClusterContainer(n);

  fClusterCollArray.Add(cont);

  return cont;
}

/**
 * Get \f$ i^{th} \f$ particle container attached to this task
 * @param[in] i Index of the particle container
 * @return Particle container found for the given index (NULL if no particle container exists for that index)
 */
AliParticleContainer* AliAnalysisTaskEmcalLight::GetParticleContainer(Int_t i) const
{
  if (i<0 || i>fParticleCollArray.GetEntriesFast()) return 0;
  AliParticleContainer *cont = static_cast<AliParticleContainer*>(fParticleCollArray.At(i));
  return cont;
}

/**
 * Get \f$ i^{th} \f$ cluster container attached to this task
 * @param[in] i Index of the cluster container
 * @return Cluster container found for the given index (NULL if no cluster container exists for that index)
 */
AliClusterContainer* AliAnalysisTaskEmcalLight::GetClusterContainer(Int_t i) const
{
  if (i<0 || i>fClusterCollArray.GetEntriesFast()) return 0;
  AliClusterContainer *cont = static_cast<AliClusterContainer*>(fClusterCollArray.At(i));
  return cont;
}

/**
 * Find particle container attached to this task according to its name
 * @param[in] name Name of the particle container
 * @return Particle container found under the given name
 */
AliParticleContainer* AliAnalysisTaskEmcalLight::GetParticleContainer(const char *name) const
{
  AliParticleContainer *cont = static_cast<AliParticleContainer*>(fParticleCollArray.FindObject(name));
  return cont;
}

/**
 * Find cluster container attached to this task according to its name
 * @param[in] name Name of the cluster container
 * @return Cluster container found under the given name
 */
AliClusterContainer* AliAnalysisTaskEmcalLight::GetClusterContainer(const char *name) const
{
  AliClusterContainer *cont = static_cast<AliClusterContainer*>(fClusterCollArray.FindObject(name));
  return cont;
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
      (runnumber >= 265077 && runnumber <= 999999)) {      // LHC16 Run-2 (p-Pb)
    b = kpA;
  }
  return b;
}
