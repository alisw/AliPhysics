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
// The output can be directly saved into a file (see the class second
// constructor). Saving through AliAnalysisManager is always happen with
// TObject::kSingleKey option. In such case an output from different runs
// cannot be merged later, because here we have histograms which are
// run-dependent, i.e. not the same between every analysis instance.
//
// If you forget to initialize EMCAL/PHOS geometry, the class will do it
// for you.
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
#include "AliAnalysisManager.h"
#include <AliVEvent.h>
#include <AliVCaloCells.h>
#include <AliVCluster.h>
#include <AliVVertex.h>
#include <AliEMCALGeometry.h>
#include <AliPHOSGeometry.h>
#include <AliLog.h>

ClassImp(AliAnalysisTaskCaloCellsQA)

//________________________________________________________________
AliAnalysisTaskCaloCellsQA::AliAnalysisTaskCaloCellsQA() : AliAnalysisTaskSE(),
  fkAvoidPileup(kTRUE),
  fCellsQA(0),
  fOutfile(""),
  fNBad(0),
  fBadCells(0)
{
  // Constructor for root I/O, do not use it
}

//________________________________________________________________
AliAnalysisTaskCaloCellsQA::AliAnalysisTaskCaloCellsQA(const char *name, Int_t nmods, Int_t det, char *outfile) :
  AliAnalysisTaskSE(name),
  fkAvoidPileup(kTRUE),
  fCellsQA(0),
  fOutfile(""),
  fNBad(0),
  fBadCells(0)
{
  // Main constructor.
  //
  // nmods -- number of supermodules + 1;
  // det -- detector;
  // outfile -- file name to write to; if NULL, write into output container;
  //   allows to avoid limitation of AliAnalysisManager, which always
  //   writes with TObject::kSingleKey option.


  if (det == kEMCAL)
    fCellsQA = new AliCaloCellsQA(nmods, AliCaloCellsQA::kEMCAL);
  else if (det == kPHOS)
    fCellsQA = new AliCaloCellsQA(nmods, AliCaloCellsQA::kPHOS);
  else
    AliFatal("Wrong detector provided");

  if (outfile) fOutfile = outfile;
  else
    DefineOutput(1, TObjArray::Class());
}

//________________________________________________________________
AliAnalysisTaskCaloCellsQA::~AliAnalysisTaskCaloCellsQA()
{
  if (!AliAnalysisManager::GetAnalysisManager()->IsProofMode()) delete fCellsQA;
  if (fBadCells) delete [] fBadCells;
}

//________________________________________________________________
void AliAnalysisTaskCaloCellsQA::UserCreateOutputObjects()
{
  // Per run histograms cannot be initialized here

  fCellsQA->InitSummaryHistograms();

  if (fOutfile.Length() == 0)
    PostData(1, fCellsQA->GetListOfHistos());
}

//________________________________________________________________
void AliAnalysisTaskCaloCellsQA::UserExec(Option_t *)
{
  // Does the job for one event

  // event
  AliVEvent *event = InputEvent();
  if (!event) {
    AliWarning("Could not get event");
    return;
  }

  fCellsQA->InitTransientFindCurrentRun(event->GetRunNumber());

  // check geometry
  if (fCellsQA->GetDetector() == AliCaloCellsQA::kEMCAL) {
    if (!AliEMCALGeometry::GetInstance()) {
      AliInfo("EMCAL geometry not initialized, initializing it for you");
      AliEMCALGeometry::GetInstance("EMCAL_COMPLETEV1");
    }
  } else {
    if (!AliPHOSGeometry::GetInstance()) {
      AliInfo("PHOS geometry not initialized, initializing it for you");
      AliPHOSGeometry::GetInstance("IHEP");
    }
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
    AliWarning("Could not get cells");
    return;
  }

  // primary vertex
  AliVVertex *vertex = (AliVVertex*) event->GetPrimaryVertex();
  if (!vertex) {
    AliWarning("Could not get primary vertex");
    return;
  }

  // collect clusters
  TObjArray clusArray;
  for (Int_t i = 0; i < event->GetNumberOfCaloClusters(); i++) {
    AliVCluster *clus = event->GetCaloCluster(i);
    if (!clus) {
      AliWarning("Could not get cluster");
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

  if (fOutfile.Length() == 0)
    PostData(1, fCellsQA->GetListOfHistos());
}

//________________________________________________________________
void AliAnalysisTaskCaloCellsQA::Terminate(Option_t*)
{
  // Handle direct saving of the histograms into a file

  if (fOutfile.Length() > 0) {
    TFile f(fOutfile.Data(), "RECREATE");
    fCellsQA->GetListOfHistos()->Write();
  }
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
