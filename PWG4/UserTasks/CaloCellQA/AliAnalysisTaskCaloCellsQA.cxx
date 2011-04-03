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
// Container class for bad channels & bad runs identification
//
// See AddTaskCaloCellsQA.C for usage example.
//
//----
//  Author: Olga Driga (SUBATECH)

// --- ROOT system ---
#include <TFile.h>

// --- AliRoot header files ---
#include <AliAnalysisTaskCaloCellsQA.h>
#include <AliVEvent.h>
#include <AliVCaloCells.h>
#include <AliVCluster.h>
#include <AliVVertex.h>

ClassImp(AliAnalysisTaskCaloCellsQA)

//________________________________________________________________
AliAnalysisTaskCaloCellsQA::AliAnalysisTaskCaloCellsQA(const char *name) : AliAnalysisTaskSE(name),
  fClusArray(0),
  fCellsQA(0),
  fOutfile(0)
{
  fClusArray = new TObjArray;
  fCellsQA = NULL;
  fOutfile = NULL;
}

//________________________________________________________________
AliAnalysisTaskCaloCellsQA::~AliAnalysisTaskCaloCellsQA()
{
  delete fClusArray;
  if (fCellsQA) delete fCellsQA;
  if (fOutfile) delete fOutfile;
}

//________________________________________________________________
void AliAnalysisTaskCaloCellsQA::InitCaloCellsQA(char* fname, Int_t nmods, Int_t det)
{
  // Must be called at the very beginning.
  //
  // fname -- output file name;
  // nmods -- number of supermodules + 1;
  // det -- detector;
  // initGeom -- if true, initialize geometry for you.

  if (det == kEMCAL)
    fCellsQA = new AliCaloCellsQA(nmods, AliCaloCellsQA::kEMCAL);
  else if (det == kPHOS)
    fCellsQA = new AliCaloCellsQA(nmods, AliCaloCellsQA::kPHOS);
  else
    Fatal("AliAnalysisTaskCaloCellsQA::InitCellsQA", "Wrong detector provided");

  fOutfile = new TString(fname);
}

//________________________________________________________________
void AliAnalysisTaskCaloCellsQA::UserCreateOutputObjects()
{
  // Per run histograms cannot be initialized here

  fCellsQA->InitSummaryHistograms();
}

//________________________________________________________________
void AliAnalysisTaskCaloCellsQA::UserExec(Option_t *)
{
  // Does the job for one event

  // event
  AliVEvent *event = InputEvent();
  if (!event) {
    Warning("AliAnalysisTaskCaloCellsQA::UserExec", "Could not get event");
    return;
  }

  // cells
  AliVCaloCells *cells;
  if (fCellsQA->GetDetector() == AliCaloCellsQA::kEMCAL)
    cells = event->GetEMCALCells();
  else
    cells = event->GetPHOSCells();

  if (!cells) {
    Warning("AliAnalysisTaskCaloCellsQA::UserExec", "Could not get cells");
    return;
  }

  // primary vertex
  AliVVertex *vertex = (AliVVertex*) event->GetPrimaryVertex();
  if (!vertex) {
    Warning("AliAnalysisTaskCaloCellsQA::UserExec", "Could not get primary vertex");
    return;
  }

  // collect clusters
  fClusArray->Clear();
  for (Int_t i = 0; i < event->GetNumberOfCaloClusters(); i++) {
    AliVCluster *clus = event->GetCaloCluster(i);
    if (!clus) {
      Warning("AliAnalysisTaskCaloCellsQA::UserExec", "Could not get cluster");
      return;
    }

    fClusArray->Add(clus);
  }

  Double_t vertexXYZ[3];
  vertex->GetXYZ(vertexXYZ);
  fCellsQA->Fill(event->GetRunNumber(), fClusArray, cells, vertexXYZ);
}

//________________________________________________________________
void AliAnalysisTaskCaloCellsQA::Terminate(Option_t*)
{
  // The AliCaloCellsQA analysis output should not be saved through
  // AliAnalysisManager, because it will write with TObject::kSingleKey
  // option. Such an output cannot be merged later: we have histograms
  // which are run-dependent, i.e. not the same between every analysis
  // instance.

  TFile f(fOutfile->Data(), "RECREATE");
  fCellsQA->GetListOfHistos()->Write();
}
