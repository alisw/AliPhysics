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

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// class for PMD reconstruction                                              //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "Riostream.h"
#include "AliPMDReconstructor.h"
#include "AliRunLoader.h"
#include "AliRun.h"
#include "AliPMDClusterFinder.h"
#include "AliPMDtracker.h"
#include "AliRawReader.h"
#include "AliESDPmdTrack.h"
#include "AliESDEvent.h"
#include "AliLog.h"

ClassImp(AliPMDReconstructor)

//_____________________________________________________________________________
void AliPMDReconstructor::Reconstruct(AliRunLoader* runLoader) const
{
// reconstruct clusters from digits file

  AliPMDClusterFinder *pmdClus = new AliPMDClusterFinder(runLoader);
  pmdClus->Load();
  for (Int_t iEvent = 0; iEvent < runLoader->GetNumberOfEvents(); iEvent++)
    {
      pmdClus->Digits2RecPoints(iEvent);
    }
  pmdClus->UnLoad();
  delete pmdClus;

}
// ------------------------------------------------------------------------ //

void AliPMDReconstructor::Reconstruct(AliRunLoader* runLoader,
				      AliRawReader *rawReader) const
{
// reconstruct clusters from Raw Data

  AliPMDClusterFinder pmdClus(runLoader);
  pmdClus.LoadClusters();

  Int_t iEvent = 0;
  while (rawReader->NextEvent()) {
    pmdClus.Digits2RecPoints(iEvent,rawReader);
    iEvent++;
  }
  pmdClus.UnLoadClusters();
  
}


// ------------------------------------------------------------------------ //

void AliPMDReconstructor::Reconstruct(AliRawReader *rawReader,
				      TTree *clustersTree) const
{
// reconstruct clusters from Raw Data

  AliPMDClusterFinder pmdClus;
  pmdClus.Digits2RecPoints(rawReader, clustersTree);

}

// ------------------------------------------------------------------------ //
void AliPMDReconstructor::Reconstruct(TTree *digitsTree,
				      TTree *clustersTree) const
{
// reconstruct clusters from Raw Data

  AliPMDClusterFinder pmdClus;
  pmdClus.Digits2RecPoints(digitsTree, clustersTree);

}

// ------------------------------------------------------------------------ //

//void AliPMDReconstructor::FillESD(AliRunLoader* runLoader,AliESDEvent* esd) const
//{
//  AliLoader* loader = runLoader->GetLoader("PMDLoader");
//  if (!loader) {
//    AliError("PMD loader not found");
//    return;
//  }
//  loader->LoadRecPoints("READ");
//  TTree *treeR = loader->TreeR();
//  AliPMDtracker pmdtracker;
//  pmdtracker.LoadClusters(treeR);
//  pmdtracker.Clusters2Tracks(esd);
//  loader->UnloadRecPoints();
//}

// ------------------------------------------------------------------------ //
void AliPMDReconstructor::FillESD(AliRawReader* /*rawReader*/,
				  TTree* clustersTree, AliESDEvent* esd) const
{
  AliPMDtracker pmdtracker;
  pmdtracker.LoadClusters(clustersTree);
  pmdtracker.Clusters2Tracks(esd);
}
// ------------------------------------------------------------------------ //
void AliPMDReconstructor::FillESD(TTree * /*digitsTree*/,
				  TTree* clustersTree, AliESDEvent* esd) const
{
  AliPMDtracker pmdtracker;
  pmdtracker.LoadClusters(clustersTree);
  pmdtracker.Clusters2Tracks(esd);
}


