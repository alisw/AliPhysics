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

#include "AliPIDResponse.h"

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

const TString AliAnalysisTaskAO2Dconverter::TreeName[kTrees] = { "O2events", "O2tracks", "O2calo", "O2tof" };

const TString AliAnalysisTaskAO2Dconverter::TreeTitle[kTrees] = { "Event tree", "Barrel tracks", "Calorimeter cells", "TOF hits" };

TTree* AliAnalysisTaskAO2Dconverter::CreateTree(TreeIndex t)
{
  fTree[t] = new TTree(TreeName[t], TreeTitle[t]);
  return fTree[t];
}

void AliAnalysisTaskAO2Dconverter::FillTree(TreeIndex t)
{
  if (!fTreeStatus[t])
    return;
  fTree[t]->Fill();
}

void AliAnalysisTaskAO2Dconverter::UserCreateOutputObjects()
{
  // create output objects
  OpenFile(1); // Necessary for large outputs

  // Associate branches for fEventTree
  TTree* Events = CreateTree(kEvents);
  Events->SetAutoFlush(fNumberOfEventsPerCluster);
  if (fTreeStatus[kEvents]) {
    Events->Branch("fEventId", &fEventId, "fEventId/l");
    Events->Branch("fVtxX", &fVtxX, "fVtxX/F");
    Events->Branch("fVtxY", &fVtxY, "fVtxY/F");
    Events->Branch("fVtxZ", &fVtxZ, "fVtxZ/F");
    Events->Branch("fCentFwd", &fCentFwd, "fCentFwd/F");
    Events->Branch("fCentBarrel", &fCentBarrel, "fCentBarrel/F");
    Events->Branch("fEventTime", &fEventTime, "fEventTime[10]/F");
    Events->Branch("fEventTimeRes", &fEventTimeRes, "fEventTimeRes[10]/F");
    Events->Branch("fEventTimeMask", &fEventTimeMask, "fEventTimeMask[10]/b");
  }
  PostData(1, Events);

  // Associate branches for fTrackTree
  TTree* Tracks = CreateTree(kTracks);
  Tracks->SetAutoFlush(fNumberOfEventsPerCluster);
  if (fTreeStatus[kTracks]) {
    Tracks->Branch("fEventId", &fEventId, "fEventId/l"); // same
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
  }
  PostData(2, Tracks);

  // Associate branches for Calo
  TTree* Calo = CreateTree(kCalo);
  Calo->SetAutoFlush(fNumberOfEventsPerCluster);
  if (fTreeStatus[kCalo]) {
    Calo->Branch("fEventId", &fEventId, "fEventId/l"); // same
    Calo->Branch("fCellNumber", &fCellNumber, "fCellNumber/S");
    Calo->Branch("fAmplitude", &fAmplitude, "fAmplitude/F");
    Calo->Branch("fTime", &fTime, "fTime/F");
    Calo->Branch("fType", &fType, "fType/B");
  }
  PostData(3, Calo);

  // Associate branches for TOF
  TTree* TOF = CreateTree(kTOF);
  TOF->SetAutoFlush(fNumberOfEventsPerCluster);
  if (fTreeStatus[kTOF]) {
    TOF->Branch("fEventId", &fEventId, "fEventId/l"); // same
    TOF->Branch("fTOFChannel", &fTOFChannel, "fTOFChannel/I");
    TOF->Branch("fTOFncls", &fTOFncls, "fTOFncls/S");
    TOF->Branch("fDx", &fDx, "fDx/F");
    TOF->Branch("fDz", &fDz, "fDz/F");
    TOF->Branch("fToT", &fToT, "fToT/F");
  }
  PostData(4, TOF);

  Prune(); //Removing all unwanted branches (if any)
}

void AliAnalysisTaskAO2Dconverter::Prune()
{
  if (fPruneList.IsNull() || fPruneList.IsWhitespace())
    return;
  TObjArray* arr = fPruneList.Tokenize(" ");
  for (Int_t i = 0; i < arr->GetEntries(); i++)
    for (Int_t j = 0; j < kTrees; j++) {
      TObjArray* branches = fTree[j]->GetListOfBranches();
      for (Int_t k = 0; k < branches->GetEntries(); k++) {
        TString bname = branches->At(k)->GetName();
        if (!bname.EqualTo(arr->At(i)->GetName()))
          continue;
        fTree[j]->SetBranchStatus(bname, 0);
      }
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
  AliPIDResponse* PIDResponse = (AliPIDResponse*)((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->GetPIDResponse();
  PIDResponse->SetTOFResponse(fESD, AliPIDResponse::kBest_T0);
  AliTOFPIDResponse TOFResponse = PIDResponse->GetTOFResponse();
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
  //Posting data
  for (Int_t i = 0; i < kTrees; i++)
    PostData(1 + i, fTree[i]);
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
