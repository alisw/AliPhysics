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

/* $Id$ */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  This class creates and fills the monitor histograms for the TPC          //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////


#include "AliMonitorTPC.h"
#include "AliMonitorDataTPC.h"
#include "AliMonitorHisto.h"
#include "AliMonitorTrend.h"
#include "AliTPCParam.h"
#include "AliTPCRawStream.h"
#include "AliTPCClustersRow.h"
#include "AliTPCclusterMI.h"
#include "AliTPCtrack.h"
#include "AliRunLoader.h"
#include <TFolder.h>
#include <TTree.h>


ClassImp(AliMonitorTPC) 


//_____________________________________________________________________________
AliMonitorTPC::AliMonitorTPC(AliTPCParam* param)
{
// create a TPC monitor object with the given parameters

  fParam = param;
  fData = new AliMonitorDataTPC(10000);
}

//_____________________________________________________________________________
AliMonitorTPC::AliMonitorTPC(const AliMonitorTPC& monitor) :
  AliMonitor(monitor)
{
  Fatal("AliMonitorTPC", "copy constructor not implemented");
}

//_____________________________________________________________________________
AliMonitorTPC& AliMonitorTPC::operator = (const AliMonitorTPC& /*monitor*/)
{
  Fatal("operator =", "assignment operator not implemented");
  return *this;
}

//_____________________________________________________________________________
AliMonitorTPC::~AliMonitorTPC()
{
  delete fData;
}


//_____________________________________________________________________________
void AliMonitorTPC::CreateHistos(TFolder* folder)
{
// create the TPC monitor histograms

  fFolder = folder->AddFolder("TPC", "TPC");

  fPadsCharge = CreateHisto1("PadsCharge", "charge distribution of pads", 
			     100, 0, 200, "charge", "#Delta N/N",
			     AliMonitorHisto::kNormEvents);

  fClustersCharge = CreateHisto1("ClustersCharge", 
				 "charge distribution of clusters", 
				 100, 0, 1000, "charge", "#Delta N/N",
				 AliMonitorHisto::kNormEvents);

  Int_t nRows = fParam->GetNRowLow() + fParam->GetNRowUp();
  fNClustersVsRow = CreateHisto1("NClustersVsRow", 
				 "mean number of clusters per pad row", 
				 nRows, -0.5, nRows-0.5,
				 "pad row", "<N_{clusters}>",
				 AliMonitorHisto::kNormEvents);

  Int_t nSector = fParam->GetNInnerSector();
  fNClustersVsSector = CreateHisto1("NClustersVsSector", 
				    "mean number of clusters per sector", 
				    nSector, -0.5, nSector-0.5, 
				    "sector", "<N_{clusters}>",
				    AliMonitorHisto::kNormEvents);

  fNTracks = CreateTrend("NTracks", "number of tracks per event", 
			 "N_{tracks}");

  fTrackPt = CreateHisto1("TrackPt", "pt distribution of tracks", 
			  90, 0, 3, "p_{t} [GeV/c]", "#Delta N/N",
			  AliMonitorHisto::kNormNone);

  fTrackEta = CreateHisto1("TrackEta", "eta distribution of tracks", 
			   100, -2, 2, "#eta", "#Delta N/N",
			   AliMonitorHisto::kNormEntries);

  fTrackPhi = CreateHisto1("TrackPhi", "phi distribution of tracks", 
			   120, 0, 360, "#phi [#circ]", "#Delta N/N",
			   AliMonitorHisto::kNormEntries);
  fTrackPhi->SetDescription("The phi distribution should be flat on average.\nIf it is not flat check for dead TPC sectors.");

  fTrackNCl = CreateHisto1("TrackNCl", "Number of clusters per track", 
			   200, 0, 200, "N_{clusters}", "#Delta N/N",
			   AliMonitorHisto::kNormNone);

  fTrackDEdxVsP = CreateHisto2("TrackDEdxVsP", "dE/dx of tracks", 
			       100, 0, 3, 100, 0, 200, 
			       "p [GeV/c]", "dE/dx", "#Delta N/N",
			       AliMonitorHisto::kNormEntries);

  fTrackDEdx = CreateHisto1("TrackDEdx", "dE/dx of tracks with 0.4<p<1.0 GeV/c", 
			       50, 0, 100, 
			       "dE/dx", "#Delta N/N",
			       AliMonitorHisto::kNormEntries);

  fTrackEtaVsPhi = CreateHisto2("TrackEtaVsPhi", "#phi vs #eta", 
			       20, -1, 1, 25, 0, 360, 
			       "#eta", "#phi", "#Delta N/N",
			       AliMonitorHisto::kNormNone);

  fPtEtaVsPhi = CreateHisto2("PtEtaVsPhi", "#phi vs #eta", 
			       20, -1, 1, 25, 0, 360, 
			       "#eta", "#phi", "#Delta N/N",
			       AliMonitorHisto::kNormNone);

}


//_____________________________________________________________________________
void AliMonitorTPC::CreateBranches(TTree* tree)
{
// create a branch with TPC variables

  tree->Branch("TPC", "AliMonitorDataTPC", &fData);
}


//_____________________________________________________________________________
void AliMonitorTPC::FillHistos(AliRunLoader* runLoader, 
			       AliRawReader* rawReader)
{
// fill the TPC monitor histogrms

  rawReader->Reset();
  AliTPCRawStream input(rawReader);
  while (input.Next()) {
    fPadsCharge->Fill(input.GetSignal());
  }

  AliLoader* tpcLoader = runLoader->GetLoader("TPCLoader");
  if (!tpcLoader) return;


  tpcLoader->LoadRecPoints();
  TTree* clusters = tpcLoader->TreeR();
  if (!clusters) return;
  AliTPCClustersRow* clustersRow = new AliTPCClustersRow;
  clustersRow->SetClass("AliTPCclusterMI");
  clusters->SetBranchAddress("Segment", &clustersRow);

  for (Int_t i = 0; i < clusters->GetEntries(); i++) {
    clusters->GetEntry(i);
    Int_t iSector, iRow;
    fParam->AdjustSectorRow(clustersRow->GetID(), iSector, iRow);
    if (iSector >= fParam->GetNInnerSector()) {
      iRow += fParam->GetNRowLow();
      iSector -= fParam->GetNInnerSector();
    }

    TObjArray* array = clustersRow->GetArray();
    for (Int_t j = 0; j < array->GetEntriesFast(); j++) {
      AliTPCclusterMI* cluster = (AliTPCclusterMI*) array->At(j);
      fClustersCharge->Fill(cluster->GetQ());
      fNClustersVsRow->Fill(iRow);
      fNClustersVsSector->Fill(iSector);
    }
  }
  fNClustersVsSector->ScaleErrorBy(10.);

  delete clustersRow;
  tpcLoader->UnloadRecPoints();


  tpcLoader->LoadTracks();
  TTree* tracks = tpcLoader->TreeT();
  if (!tracks) return;
  AliTPCtrack* track = new AliTPCtrack;
  tracks->SetBranchAddress("tracks", &track);

  Int_t nTracks = (Int_t) tracks->GetEntries();
  fNTracks->Fill(nTracks);
  fData->SetNTracks(nTracks);
  for (Int_t i = 0; i < nTracks; i++) {
    tracks->GetEntry(i);
    fTrackPt->Fill(track->Pt());
    fTrackEta->Fill(track->Eta());
    fTrackPhi->Fill(track->Phi() * TMath::RadToDeg());
    if(track->Pt()>3.) {
      fTrackEtaVsPhi->Fill(track->Eta(),track->Phi() * TMath::RadToDeg());
      fPtEtaVsPhi->Fill(track->Eta(),track->Phi() * TMath::RadToDeg(),track->Pt());
    }
    fTrackNCl->Fill(track->GetNumberOfClusters());
    fTrackDEdxVsP->Fill(track->P(), track->GetdEdx());
    if(track->P()>0.4 && track->P()<1.0)
      fTrackDEdx->Fill(track->GetdEdx());

    fData->SetData(i, track->Pt(), track->Eta(), 
		   track->Phi() * TMath::RadToDeg());
  }

  delete track;
  tpcLoader->UnloadTracks();
}
