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

AliAnalysisTaskAO2Dconverter::AliAnalysisTaskAO2Dconverter(const char *name) : AliAnalysisTaskSE(name),
                                                                               fEventCuts{}
{
  DefineInput(0, TChain::Class());
  DefineOutput(1, TTree::Class());
  DefineOutput(2, TTree::Class());
  DefineOutput(3, TTree::Class());
}

AliAnalysisTaskAO2Dconverter::~AliAnalysisTaskAO2Dconverter()
{
  if (fEventTree)
    delete fEventTree;
  if (fTrackTree)
    delete fTrackTree;
  if (fCaloTree)
    delete fCaloTree;
}
void AliAnalysisTaskAO2Dconverter::UserCreateOutputObjects()
{
  // create output objects
  OpenFile(1); // Necessary for large outputs

  // Associate branches for fEventTree
  fEventTree = new TTree("O2events", "Event tree");
  fEventTree->SetAutoFlush(fNumberOfEventsPerCluster);
  fEventTree->Branch("fEventId", &fEventId, "fEventId/l");
  fEventTree->Branch("fVtxX", &fVtxX, "fVtxX/F");
  fEventTree->Branch("fVtxY", &fVtxY, "fVtxY/F");
  fEventTree->Branch("fVtxZ", &fVtxZ, "fVtxZ/F");
  fEventTree->Branch("fCentFwd", &fCentFwd, "fCentFwd/F");
  fEventTree->Branch("fCentBarrel", &fCentBarrel, "fCentBarrel/F");
  PostData(1, fEventTree);

  // Associate branches for fTrackTree
  fTrackTree = new TTree("O2tracks", "Barrel tracks");
  fTrackTree->SetAutoFlush(fNumberOfEventsPerCluster);
  fTrackTree->Branch("fEventId", &fEventId, "fEventId/l"); // same
  fTrackTree->Branch("fX", &fX, "fX/F");
  fTrackTree->Branch("fAlpha", &fAlpha, "fAlpha/F");
  fTrackTree->Branch("fY", &fY, "fY/F");
  fTrackTree->Branch("fZ", &fZ, "fZ/F");
  fTrackTree->Branch("fSnp", &fSnp, "fSnp/F");
  fTrackTree->Branch("fTgl", &fTgl, "fTgl/F");
  fTrackTree->Branch("fSigned1Pt", &fSigned1Pt, "fSigned1Pt/F");
  fTrackTree->Branch("fCYY", &fCYY, "fCYY/F");
  fTrackTree->Branch("fCZY", &fCZY, "fCZY/F");
  fTrackTree->Branch("fCZZ", &fCZZ, "fCZZ/F");
  fTrackTree->Branch("fCSnpY", &fCSnpY, "fCSnpY/F");
  fTrackTree->Branch("fCSnpZ", &fCSnpZ, "fCSnpZ/F");
  fTrackTree->Branch("fCSnpSnp", &fCSnpSnp, "fCSnpSnp/F");
  fTrackTree->Branch("fCTglY", &fCTglY, "fCTglY/F");
  fTrackTree->Branch("fCTglZ", &fCTglZ, "fCTglZ/F");
  fTrackTree->Branch("fCTglSnp", &fCTglSnp, "fCTglSnp/F");
  fTrackTree->Branch("fCTglTgl", &fCTglTgl, "fCTglTgl/F");
  fTrackTree->Branch("fC1PtY", &fC1PtY, "fC1PtY/F");
  fTrackTree->Branch("fC1PtZ", &fC1PtZ, "fC1PtZ/F");
  fTrackTree->Branch("fC1PtSnp", &fC1PtSnp, "fC1PtSnp/F");
  fTrackTree->Branch("fC1PtTgl", &fC1PtTgl, "fC1PtTgl/F");
  fTrackTree->Branch("fC1Pt21Pt2", &fC1Pt21Pt2, "fC1Pt21Pt2/F");
  fTrackTree->Branch("fTPCinnerP", &fTPCinnerP, "fTPCinnerP/F");
  fTrackTree->Branch("fFlags", &fFlags, "fFlags/l");
  fTrackTree->Branch("fITSClusterMap", &fITSClusterMap, "fITSClusterMap/b");
  fTrackTree->Branch("fTPCncls", &fTPCncls, "fTPCncls/s");
  fTrackTree->Branch("fTRDntracklets", &fTRDntracklets, "fTRDntracklets/b");
  fTrackTree->Branch("fITSchi2Ncl", &fITSchi2Ncl, "fITSchi2Ncl/F");
  fTrackTree->Branch("fTPCchi2Ncl", &fTPCchi2Ncl, "fTPCchi2Ncl/F");
  fTrackTree->Branch("fTRDchi2", &fTRDchi2, "fTRDchi2/F");
  fTrackTree->Branch("fTOFchi2", &fTOFchi2, "fTOFchi2/F");
  fTrackTree->Branch("fTPCsignal", &fTPCsignal, "fTPCsignal/F");
  fTrackTree->Branch("fTRDsignal", &fTRDsignal, "fTRDsignal/F");
  fTrackTree->Branch("fTOFsignal", &fTOFsignal, "fTOFsignal/F");
  fTrackTree->Branch("fLength", &fLength, "fLength/F");
  PostData(2, fTrackTree);

  fCaloTree = new TTree("O2calo", "Calorimeter cells");
  fCaloTree->Branch("fEventId", &fEventId, "fEventId/l"); // same
  fCaloTree->Branch("fCellNumber", &fCellNumber, "fCellNumber/S");
  fCaloTree->Branch("fAmplitude", &fAmplitude, "fAmplitude/F");
  fCaloTree->Branch("fTime", &fTime, "fTime/F");
  fCaloTree->Branch("fType", &fType, "fType/B");
  PostData(3, fCaloTree);
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
  fEventTree->Fill();

  // Fill fTrackTree
  Int_t ntrk = fESD->GetNumberOfTracks();
  for (Int_t itrk = 0; itrk < ntrk; itrk++)
  {
    AliESDtrack *track = fESD->GetTrack(itrk);

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

    fTrackTree->Fill();

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

    fCaloTree->Fill();
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

    fCaloTree->Fill();
  } // end loop on PHOS cells

  PostData(1, fEventTree);
  PostData(2, fTrackTree);
  PostData(3, fCaloTree);
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
  mgr->ConnectOutput(task, 1, mgr->CreateContainer("O2event", TTree::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));
  mgr->ConnectOutput(task, 2, mgr->CreateContainer("O2tracks", TTree::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));
  mgr->ConnectOutput(task, 3, mgr->CreateContainer("O2calo", TTree::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));
  // in the end, this macro returns a pointer to your task. this will be convenient later on
  // when you will run your analysis in an analysis train on grid
  return task;
}
