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

/* AliAnalysisTaskAO2Dconverter
 *
 * Convert Run 2 ESDs to Run 3 prototype AODs (AliAO2D.root).
 */

#include "TChain.h"
#include "TTree.h"
#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliAnalysisTaskAO2Dconverter.h"
#include "AliVHeader.h"

#include "AliESDCaloCells.h"
#include "AliESDHeader.h"
#include "AliESDtrack.h"

#include "AliMCEvent.h"
#include "AliMCEventHandler.h"
#include "AliPIDResponse.h"

#include "AliGenCocktailEventHeader.h"
#include "AliGenDPMjetEventHeader.h"
#include "AliGenEpos3EventHeader.h"
#include "AliGenEposEventHeader.h"
#include "AliGenEventHeader.h"
#include "AliGenEventHeaderTunedPbPb.h"
#include "AliGenGeVSimEventHeader.h"
#include "AliGenHepMCEventHeader.h"
#include "AliGenHerwigEventHeader.h"
#include "AliGenHijingEventHeader.h"
#include "AliGenPythiaEventHeader.h"
#include "AliGenToyEventHeader.h"

ClassImp(AliAnalysisTaskAO2Dconverter);

namespace
{

ULong64_t GetEventIdAsLong(AliVHeader *header)
{
  return ((ULong64_t)header->GetBunchCrossNumber() +
          (ULong64_t)header->GetOrbitNumber() * 3564 +
          (ULong64_t)header->GetPeriodNumber() * 16777215 * 3564);
}

} // namespace

AliAnalysisTaskAO2Dconverter::AliAnalysisTaskAO2Dconverter(const char* name)
    : AliAnalysisTaskSE(name)
    , fTrackFilter(Form("AO2Dconverter%s", name), Form("fTrackFilter%s", name))
    , fEventCuts{}
{
  DefineInput(0, TChain::Class());
  for (Int_t i = 0; i < kTrees; i++) {
    fTreeStatus[i] = kTRUE;
    DefineOutput(1 + i, TTree::Class());
  }
}

AliAnalysisTaskAO2Dconverter::~AliAnalysisTaskAO2Dconverter()
{
  for (Int_t i = 0; i < kTrees; i++)
    if (fTree[i])
      delete fTree[i];
}

const TString AliAnalysisTaskAO2Dconverter::TreeName[kTrees] = { "O2events", "O2tracks", "O2calo", "O2tof", "O2kine" };

const TString AliAnalysisTaskAO2Dconverter::TreeTitle[kTrees] = { "Event tree", "Barrel tracks", "Calorimeter cells", "TOF hits", "Kinematics" };

const TClass* AliAnalysisTaskAO2Dconverter::Generator[kGenerators] = { AliGenEventHeader::Class(), AliGenCocktailEventHeader::Class(), AliGenDPMjetEventHeader::Class(), AliGenEpos3EventHeader::Class(), AliGenEposEventHeader::Class(), AliGenEventHeaderTunedPbPb::Class(), AliGenGeVSimEventHeader::Class(), AliGenHepMCEventHeader::Class(), AliGenHerwigEventHeader::Class(), AliGenHijingEventHeader::Class(), AliGenPythiaEventHeader::Class(), AliGenToyEventHeader::Class() };

TTree* AliAnalysisTaskAO2Dconverter::CreateTree(TreeIndex t)
{
  fTree[t] = new TTree(TreeName[t], TreeTitle[t]);
  if (fTreeStatus[t])
    fTree[t]->Branch("fEventId", &fEventId, "fEventId/l"); // Branch common to all trees
  return fTree[t];
}

void AliAnalysisTaskAO2Dconverter::PostTree(TreeIndex t)
{
  if (!fTreeStatus[t])
    return;
  PostData(t + 1, fTree[t]);
}

void AliAnalysisTaskAO2Dconverter::FillTree(TreeIndex t)
{
  if (!fTreeStatus[t])
    return;
  fTree[t]->Fill();
}

void AliAnalysisTaskAO2Dconverter::UserCreateOutputObjects()
{
  switch (fTaskMode) { // Setting active/inactive containers based on the TaskMode
  case kStandard:
    DisableTree(kKinematics);
    break;
  default:
    break;
  }

  // create output objects
  OpenFile(1); // Necessary for large outputs

  // Associate branches for fEventTree
  TTree* Events = CreateTree(kEvents);
  Events->SetAutoFlush(fNumberOfEventsPerCluster);
  if (fTreeStatus[kEvents]) {
    Events->Branch("fVtxX", &fVtxX, "fVtxX/F");
    Events->Branch("fVtxY", &fVtxY, "fVtxY/F");
    Events->Branch("fVtxZ", &fVtxZ, "fVtxZ/F");
    Events->Branch("fCentFwd", &fCentFwd, "fCentFwd/F");
    Events->Branch("fCentBarrel", &fCentBarrel, "fCentBarrel/F");
    Events->Branch("fEventTime", &fEventTime, "fEventTime[10]/F");
    Events->Branch("fEventTimeRes", &fEventTimeRes, "fEventTimeRes[10]/F");
    Events->Branch("fEventTimeMask", &fEventTimeMask, "fEventTimeMask[10]/b");
    if (fTaskMode == kMC) {
      Events->Branch("fGeneratorID", &fGeneratorID, "fGeneratorID/S");
      Events->Branch("fMCVtxX", &fMCVtxX, "fMCVtxX/F");
      Events->Branch("fMCVtxY", &fMCVtxY, "fMCVtxY/F");
      Events->Branch("fMCVtxZ", &fMCVtxZ, "fMCVtxZ/F");
    }
  }
  PostTree(kEvents);

  // Associate branches for fTrackTree
  TTree* Tracks = CreateTree(kTracks);
  Tracks->SetAutoFlush(fNumberOfEventsPerCluster);
  if (fTreeStatus[kTracks]) {
    Tracks->Branch("fX", &fX, "fX/F");
    Tracks->Branch("fAlpha", &fAlpha, "fAlpha/F");
    Tracks->Branch("fY", &fY, "fY/F");
    Tracks->Branch("fZ", &fZ, "fZ/F");
    Tracks->Branch("fSnp", &fSnp, "fSnp/F");
    Tracks->Branch("fTgl", &fTgl, "fTgl/F");
    Tracks->Branch("fSigned1Pt", &fSigned1Pt, "fSigned1Pt/F");
    Tracks->Branch("fCYY", &fCYY, "fCYY/F");
    Tracks->Branch("fCZY", &fCZY, "fCZY/F");
    Tracks->Branch("fCZZ", &fCZZ, "fCZZ/F");
    Tracks->Branch("fCSnpY", &fCSnpY, "fCSnpY/F");
    Tracks->Branch("fCSnpZ", &fCSnpZ, "fCSnpZ/F");
    Tracks->Branch("fCSnpSnp", &fCSnpSnp, "fCSnpSnp/F");
    Tracks->Branch("fCTglY", &fCTglY, "fCTglY/F");
    Tracks->Branch("fCTglZ", &fCTglZ, "fCTglZ/F");
    Tracks->Branch("fCTglSnp", &fCTglSnp, "fCTglSnp/F");
    Tracks->Branch("fCTglTgl", &fCTglTgl, "fCTglTgl/F");
    Tracks->Branch("fC1PtY", &fC1PtY, "fC1PtY/F");
    Tracks->Branch("fC1PtZ", &fC1PtZ, "fC1PtZ/F");
    Tracks->Branch("fC1PtSnp", &fC1PtSnp, "fC1PtSnp/F");
    Tracks->Branch("fC1PtTgl", &fC1PtTgl, "fC1PtTgl/F");
    Tracks->Branch("fC1Pt21Pt2", &fC1Pt21Pt2, "fC1Pt21Pt2/F");
    Tracks->Branch("fTPCinnerP", &fTPCinnerP, "fTPCinnerP/F");
    Tracks->Branch("fFlags", &fFlags, "fFlags/l");
    Tracks->Branch("fITSClusterMap", &fITSClusterMap, "fITSClusterMap/b");
    Tracks->Branch("fTPCncls", &fTPCncls, "fTPCncls/s");
    Tracks->Branch("fTRDntracklets", &fTRDntracklets, "fTRDntracklets/b");
    Tracks->Branch("fITSchi2Ncl", &fITSchi2Ncl, "fITSchi2Ncl/F");
    Tracks->Branch("fTPCchi2Ncl", &fTPCchi2Ncl, "fTPCchi2Ncl/F");
    Tracks->Branch("fTRDchi2", &fTRDchi2, "fTRDchi2/F");
    Tracks->Branch("fTOFchi2", &fTOFchi2, "fTOFchi2/F");
    Tracks->Branch("fTPCsignal", &fTPCsignal, "fTPCsignal/F");
    Tracks->Branch("fTRDsignal", &fTRDsignal, "fTRDsignal/F");
    Tracks->Branch("fTOFsignal", &fTOFsignal, "fTOFsignal/F");
    Tracks->Branch("fLength", &fLength, "fLength/F");
    Tracks->Branch("fLabel", &fLabel, "fLabel/I");
    Tracks->Branch("fTOFLabel", &fTOFLabel, "fTOFLabel[3]/I");
  }
  PostTree(kTracks);

  // Associate branches for Calo
  TTree* Calo = CreateTree(kCalo);
  Calo->SetAutoFlush(fNumberOfEventsPerCluster);
  if (fTreeStatus[kCalo]) {
    Calo->Branch("fCellNumber", &fCellNumber, "fCellNumber/S");
    Calo->Branch("fAmplitude", &fAmplitude, "fAmplitude/F");
    Calo->Branch("fTime", &fTime, "fTime/F");
    Calo->Branch("fType", &fType, "fType/B");
  }
  PostTree(kCalo);

  // Associate branches for TOF
  TTree* TOF = CreateTree(kTOF);
  TOF->SetAutoFlush(fNumberOfEventsPerCluster);
  if (fTreeStatus[kTOF]) {
    TOF->Branch("fTOFChannel", &fTOFChannel, "fTOFChannel/I");
    TOF->Branch("fTOFncls", &fTOFncls, "fTOFncls/S");
    TOF->Branch("fDx", &fDx, "fDx/F");
    TOF->Branch("fDz", &fDz, "fDz/F");
    TOF->Branch("fToT", &fToT, "fToT/F");
  }
  PostTree(kTOF);

  // Associate branches for Kinematics
  TTree* Kinematics = CreateTree(kKinematics);
  Kinematics->SetAutoFlush(fNumberOfEventsPerCluster);
  if (fTreeStatus[kMC]) {
    Kinematics->Branch("fPdgCode", &fPdgCode, "fPdgCode/I");
    Kinematics->Branch("fMother", &fMother, "fMother[2]/I");
    Kinematics->Branch("fDaughter", &fDaughter, "fDaughter[2]/I");

    Kinematics->Branch("fPx", &fPx, "fPx/F");
    Kinematics->Branch("fPy", &fPy, "fPy/F");
    Kinematics->Branch("fPz", &fPz, "fPz/F");

    Kinematics->Branch("fVx", &fVx, "fVx/F");
    Kinematics->Branch("fVy", &fVy, "fVy/F");
    Kinematics->Branch("fVz", &fVz, "fVz/F");
    Kinematics->Branch("fVt", &fVt, "fVt/F");
  }
  PostTree(kKinematics);

  Prune(); //Removing all unwanted branches (if any)
}

void AliAnalysisTaskAO2Dconverter::Prune()
{
  if (fPruneList.IsNull() || fPruneList.IsWhitespace())
    return;
  TObjArray* arr = fPruneList.Tokenize(" ");
  for (Int_t i = 0; i < arr->GetEntries(); i++) {
    Bool_t found = kFALSE;
    for (Int_t j = 0; j < kTrees; j++) {
      TObjArray* branches = fTree[j]->GetListOfBranches();
      for (Int_t k = 0; k < branches->GetEntries(); k++) {
        TString bname = branches->At(k)->GetName();
        if (!bname.EqualTo(arr->At(i)->GetName()))
          continue;
        fTree[j]->SetBranchStatus(bname, 0);
        found = kTRUE;
      }
    }
    if (!found)
      AliFatal(Form("Did not find Branch %s", arr->At(i)->GetName()));
  }
  fPruneList = "";
}

void AliAnalysisTaskAO2Dconverter::UserExec(Option_t *)
{
  fESD = dynamic_cast<AliESDEvent *>(InputEvent());
  if (!fESD) {
    ::Fatal("AliAnalysisTaskAO2Dconverter::UserExec", "Something is wrong with the event handler");
  }

  // We avoid cases where we have zero reconstructed tracks
  if (!fEventCuts.AcceptEvent(fESD))
  {
    return;
  }

  // Configuration of the PID response
  AliPIDResponse* PIDResponse = (AliPIDResponse*)((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->GetPIDResponse();
  PIDResponse->SetTOFResponse(fESD, AliPIDResponse::kBest_T0);
  AliTOFPIDResponse TOFResponse = PIDResponse->GetTOFResponse();

  // Configuration of the MC event (if needed)
  AliMCEvent* MCEvt = nullptr;
  if (fTaskMode == kMC) {
    AliMCEventHandler* eventHandler = dynamic_cast<AliMCEventHandler*>(AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler()); //Get the MC handler

    if (!eventHandler) //Check on the MC handler
      AliFatal("Could not retrieve MC event handler");
    MCEvt = eventHandler->MCEvent(); //Get the MC Event

    if (!MCEvt) // Check on the MC Event
      AliFatal("Could not retrieve MC event");
    PIDResponse->SetCurrentMCEvent(MCEvt); //Set The PID response on the current MC event
  }
  const AliVVertex *vtx = fEventCuts.GetPrimaryVertex();
  if (!vtx) {
    ::Fatal("AliAnalysisTaskAO2Dconverter::UserExec", "Vertex not defined");
  }
  fEventId = GetEventIdAsLong(fESD->GetHeader());
  fVtxX = vtx->GetX();
  fVtxY = vtx->GetY();
  fVtxZ = vtx->GetZ();
  fCentFwd = fEventCuts.GetCentrality(0);
  fCentBarrel = fEventCuts.GetCentrality(1);
  for (Int_t i = 0; i < TOFResponse.GetNmomBins(); i++) {
    if (i >= 10)
      AliFatal("Index is too high!");
    Float_t mom = (TOFResponse.GetMinMom(i) + TOFResponse.GetMaxMom(i)) / 2.f;
    fEventTime[i] = TOFResponse.GetStartTime(mom);
    fEventTimeRes[i] = TOFResponse.GetStartTimeRes(mom);
    if (TOFResponse.GetStartTimeMask(mom) & 0x1)
      SETBIT(fEventTimeMask[i], 0);
    else
      CLRBIT(fEventTimeMask[i], 0);
    //
    if (TOFResponse.GetStartTimeMask(mom) & 0x2)
      SETBIT(fEventTimeMask[i], 1);
    else
      CLRBIT(fEventTimeMask[i], 1);
    //
    if (TOFResponse.GetStartTimeMask(mom) & 0x3)
      SETBIT(fEventTimeMask[i], 2);
    else
      CLRBIT(fEventTimeMask[i], 2);
  }
  if (MCEvt) {
    const AliVVertex* MCvtx = MCEvt->GetPrimaryVertex();
    if (!MCvtx) //Check on the MC vertex
      AliFatal("Could not retrieve MC vertex");
    fMCVtxX = MCvtx->GetX();
    fMCVtxY = MCvtx->GetY();
    fMCVtxZ = MCvtx->GetZ();
    AliGenEventHeader* mcGenH = MCEvt->GenEventHeader();
    for (Int_t gen = 0; gen < kGenerators; gen++) {
      if (mcGenH->InheritsFrom(Generator[gen]))
        SETBIT(fGeneratorID, gen);
      else
        CLRBIT(fGeneratorID, gen);
    }
    if (mcGenH->InheritsFrom(Generator[kAliGenCocktailEventHeader])) {
      TList* headers = ((AliGenCocktailEventHeader*)mcGenH)->GetHeaders();
      for (Int_t cocktail = 0; cocktail < headers->GetEntries(); headers++) {
        for (Int_t gen = 0; gen < kGenerators; gen++) {
          if (mcGenH->InheritsFrom(Generator[gen]))
            SETBIT(fGeneratorID, gen);
        }
      }
    }
  }
  FillTree(kEvents);

  // Fill fTrackTree
  Int_t ntrk = fESD->GetNumberOfTracks();
  for (Int_t itrk = 0; itrk < ntrk; itrk++)
  {
    AliESDtrack *track = fESD->GetTrack(itrk);
    if (!fTrackFilter.IsSelected(track))
      continue;

    fX = track->GetX();
    fAlpha = track->GetAlpha();

    fY = track->GetY();
    fZ = track->GetZ();
    fSnp = track->GetSnp();
    fTgl = track->GetTgl();
    fSigned1Pt = track->GetSigned1Pt();

    fCYY = track->GetSigmaY2();
    fCZY = track->GetSigmaZY();
    fCZZ = track->GetSigmaZ2();
    fCSnpY = track->GetSigmaSnpY();
    fCSnpZ = track->GetSigmaSnpZ();
    fCSnpSnp = track->GetSigmaSnp2();
    fCTglY = track->GetSigmaTglY();
    fCTglZ = track->GetSigmaTglZ();
    fCTglSnp = track->GetSigmaTglSnp();
    fCTglTgl = track->GetSigmaTgl2();
    fC1PtY = track->GetSigma1PtY();
    fC1PtZ = track->GetSigma1PtZ();
    fC1PtSnp = track->GetSigma1PtSnp();
    fC1PtTgl = track->GetSigma1PtTgl();
    fC1Pt21Pt2 = track->GetSigma1Pt2();

    const AliExternalTrackParam *intp = track->GetTPCInnerParam();
    fTPCinnerP = (intp ? intp->GetP() : 0); // Set the momentum to 0 if the track did not reach TPC

    fFlags = track->GetStatus();

    fITSClusterMap = track->GetITSClusterMap();
    fTPCncls = track->GetTPCNcls();
    fTRDntracklets = track->GetTRDntracklets();

    fITSchi2Ncl = (track->GetITSNcls() ? track->GetITSchi2() / track->GetITSNcls() : 0);
    fTPCchi2Ncl = (track->GetTPCNcls() ? track->GetTPCchi2() / track->GetTPCNcls() : 0);
    fTRDchi2 = track->GetTRDchi2();
    fTOFchi2 = track->GetTOFchi2();

    fTPCsignal = track->GetTPCsignal();
    fTRDsignal = track->GetTRDsignal();
    fTOFsignal = track->GetTOFsignal();
    fLength = track->GetIntegratedLength();

    fTOFncls = track->GetNTOFclusters();

    if (fTOFncls > 0) {
      Int_t* TOFclsIndex = track->GetTOFclusterArray(); //Index of the matchable cluster (there are fNTOFClusters of them)
      for (Int_t icls = 0; icls < fTOFncls; icls++) {
        AliESDTOFCluster* TOFcls = (AliESDTOFCluster*)fESD->GetESDTOFClusters()->At(TOFclsIndex[icls]);
        fToT = TOFcls->GetTOFsignalToT(0);
        fTOFChannel = TOFcls->GetTOFchannel();
        for (Int_t mtchbl = 0; mtchbl < TOFcls->GetNMatchableTracks(); mtchbl++) {
          if (TOFcls->GetTrackIndex(mtchbl) != track->GetID())
            continue;
          fDx = TOFcls->GetDx(mtchbl);
          fDz = TOFcls->GetDz(mtchbl);
          fLengthRatio = fLength > 0 ? TOFcls->GetLength(mtchbl) / fLength : -1;
          break;
        }
        FillTree(kTOF);
      }
    }

    fLabel = track->GetLabel();
    track->GetTOFLabel(fTOFLabel);

    FillTree(kTracks);

  } // end loop on tracks

  // Fill fCaloTree
  AliESDCaloCells *cells = fESD->GetEMCALCells();
  Short_t nCells = cells->GetNumberOfCells();
  for (Short_t ice = 0; ice < nCells; ++ice)
  {
    Short_t cellNumber;
    Double_t amplitude;
    Double_t time;
    Int_t mclabel;
    Double_t efrac;

    cells->GetCell(ice, cellNumber, amplitude, time, mclabel, efrac);
    fCellNumber = cellNumber;
    fAmplitude = amplitude;
    fTime = time;
    fType = cells->GetType(); // common for all cells

    FillTree(kCalo);
  } // end loop on calo cells

  cells = fESD->GetPHOSCells();
  nCells = cells->GetNumberOfCells();
  for (Short_t icp = 0; icp < nCells; ++icp)
  {
    Short_t cellNumber;
    Double_t amplitude;
    Double_t time;
    Int_t mclabel;
    Double_t efrac;

    cells->GetCell(icp, cellNumber, amplitude, time, mclabel, efrac);
    fCellNumber = cellNumber;
    fAmplitude = amplitude;
    fTime = time;
    fType = cells->GetType(); // common for all cells

    FillTree(kCalo);
  } // end loop on PHOS cells

  if (MCEvt) {
    TParticle* particle = nullptr;
    for (Int_t i = 0; i < MCEvt->GetNumberOfTracks(); i++) { //loop on primary MC tracks Before Event Selection
      particle = MCEvt->Particle(i);

      //Get the kinematic values of the particles
      fPdgCode = particle->GetPdgCode();
      fMother[0] = particle->GetFirstMother();
      fMother[1] = particle->GetSecondMother();
      fDaughter[0] = particle->GetFirstDaughter();
      fDaughter[1] = particle->GetLastDaughter();

      fPx = particle->Px();
      fPy = particle->Py();
      fPz = particle->Pz();

      fVx = particle->Vx();
      fVy = particle->Vy();
      fVz = particle->Vz();
      fVt = particle->T();

      FillTree(kKinematics);
    }
  }
  //Posting data
  for (Int_t i = 0; i < kTrees; i++)
    PostTree((TreeIndex)i);
}

void AliAnalysisTaskAO2Dconverter::Terminate(Option_t *)
{
  // terminate
  // called at the END of the analysis (when all events are processed)
}

AliAnalysisTaskAO2Dconverter *AliAnalysisTaskAO2Dconverter::AddTask(TString suffix)
{
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr)
  {
    return nullptr;
  }
  // get the input event handler, again via a static method.
  // this handler is part of the managing system and feeds events
  // to your task
  if (!mgr->GetInputEventHandler())
  {
    return nullptr;
  }
  // by default, a file is open for writing. here, we get the filename
  TString fileName = "AO2D.root";
  if (!suffix.IsNull())
    fileName += ":" + suffix; // create a subfolder in the file
  // now we create an instance of your task
  AliAnalysisTaskAO2Dconverter *task = new AliAnalysisTaskAO2Dconverter((TString("AO2D") + suffix).Data());
  if (!task)
    return nullptr;
  // add your task to the manager
  mgr->AddTask(task);
  // your task needs input: here we connect the manager to your task
  mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
  // same for the output
  for (Int_t i = 0; i < kTrees; i++)
    mgr->ConnectOutput(task, 1 + i, mgr->CreateContainer(TreeName[i], TTree::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));
  // in the end, this macro returns a pointer to your task. this will be convenient later on
  // when you will run your analysis in an analysis train on grid
  return task;
}
