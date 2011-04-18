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
// By default, pileup events are not processed, use SetAvoidPileup()
// to change this.
//
// Clusters containing a bad cell can be rejected, use SetBadCells().
//
// See AddTaskCaloCellsQA.C for usage example.
//
//----
//  Author: Olga Driga (SUBATECH)

// --- ROOT system ---
#include <TFile.h>
#include <TObjArray.h>

// --- AliRoot header files ---
#include <AliAnalysisTaskCaloCellsQA.h>
#include <AliVEvent.h>
#include <AliVCaloCells.h>
#include <AliVCluster.h>
#include <AliVVertex.h>

ClassImp(AliAnalysisTaskCaloCellsQA)

//________________________________________________________________
AliAnalysisTaskCaloCellsQA::AliAnalysisTaskCaloCellsQA() : AliAnalysisTaskSE(),
  fkAvoidPileup(kTRUE),
  fCellsQA(0),
  fOutfile(new TString),
  fBadCells(0),
  fNBad(0)
{
}

//________________________________________________________________
AliAnalysisTaskCaloCellsQA::AliAnalysisTaskCaloCellsQA(const char *name) : AliAnalysisTaskSE(name),
  fkAvoidPileup(kTRUE),
  fCellsQA(0),
  fOutfile(new TString),
  fBadCells(0),
  fNBad(0)
{
}

//________________________________________________________________
AliAnalysisTaskCaloCellsQA::~AliAnalysisTaskCaloCellsQA()
{
  if (fCellsQA) delete fCellsQA;
  delete fOutfile;
  if (fBadCells) delete [] fBadCells;
}

//________________________________________________________________
void AliAnalysisTaskCaloCellsQA::InitCaloCellsQA(char* fname, Int_t nmods, Int_t det)
{
  // Must be called at the very beginning.
  //
  // fname -- output file name;
  // nmods -- number of supermodules + 1;
  // det -- detector;

  if (det == kEMCAL)
    fCellsQA = new AliCaloCellsQA(nmods, AliCaloCellsQA::kEMCAL);
  else if (det == kPHOS)
    fCellsQA = new AliCaloCellsQA(nmods, AliCaloCellsQA::kPHOS);
  else
    Fatal("AliAnalysisTaskCaloCellsQA::InitCellsQA", "Wrong detector provided");

  *fOutfile = fname;
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

  // pileup;  FIXME: why AliVEvent does not allow a version without arguments?
  if (fkAvoidPileup && event->IsPileupFromSPD(3,0.8,3.,2.,5.))
    return;

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
  TObjArray clusArray;
  for (Int_t i = 0; i < event->GetNumberOfCaloClusters(); i++) {
    AliVCluster *clus = event->GetCaloCluster(i);
    if (!clus) {
      Warning("AliAnalysisTaskCaloCellsQA::UserExec", "Could not get cluster");
      return;
    }

    // basic filtering
    if (fCellsQA->GetDetector() == AliCaloCellsQA::kEMCAL && !clus->IsEMCAL()) continue;
    if (fCellsQA->GetDetector() == AliCaloCellsQA::kPHOS && !clus->IsPHOS()) continue;
    if (IsClusterBad(clus)) continue;

    clusArray.Add(clus);
  }

  Double_t vertexXYZ[3];
  vertex->GetXYZ(vertexXYZ);
  fCellsQA->Fill(event->GetRunNumber(), &clusArray, cells, vertexXYZ);
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

//____________________________________________________________
void AliAnalysisTaskCaloCellsQA::SetBadCells(Int_t badcells[], Int_t nbad)
{
  // Set absId numbers for bad cells;
  // clusters which contain a bad cell will be rejected.

  if (fBadCells) delete [] fBadCells;

  // switch off bad cells, if asked
  if (nbad <= 0) {
    fNBad = 0;
    return;
  }

  fNBad = nbad;
  fBadCells = new Int_t[nbad];

  for (Int_t i = 0; i < nbad; i++)
    fBadCells[i] = badcells[i];
}

//________________________________________________________________
Bool_t AliAnalysisTaskCaloCellsQA::IsClusterBad(AliVCluster *clus)
{
  // Returns true if cluster contains a bad cell

  for (Int_t b = 0; b < fNBad; b++)
    for (Int_t c = 0; c < clus->GetNCells(); c++)
      if (clus->GetCellAbsId(c) == fBadCells[b])
        return kTRUE;

  return kFALSE;
}
