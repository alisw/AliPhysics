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
#include "AliTPCRawStream.h"
#include "AliTPCClustersRow.h"
#include "AliTPCclusterMI.h"
#include "AliTPCtrack.h"


ClassImp(AliMonitorDataTPC) 


//_____________________________________________________________________________
AliMonitorDataTPC::AliMonitorDataTPC()
{
  fPt = fEta = fPhi = NULL;
  fSize = 0;
}

//_____________________________________________________________________________
AliMonitorDataTPC::AliMonitorDataTPC(Int_t size)
{
  fPt = new Float_t[size];
  fEta = new Float_t[size];
  fPhi = new Float_t[size];
  fSize = size;
}

//_____________________________________________________________________________
AliMonitorDataTPC::~AliMonitorDataTPC()
{
  delete[] fPt;
  delete[] fEta;
  delete[] fPhi;
}

//_____________________________________________________________________________
void AliMonitorDataTPC::SetSize(Int_t size)
{
  if (size > fSize) {
    delete[] fPt;
    delete[] fEta;
    delete[] fPhi;
    fPt = new Float_t[size];
    fEta = new Float_t[size];
    fPhi = new Float_t[size];
    fSize = size;
  }
}



ClassImp(AliMonitorTPC) 


//_____________________________________________________________________________
AliMonitorTPC::AliMonitorTPC(AliTPCParam* param)
{
// create a TPC monitor object with the given parameters

  fParam = param;
  fData = new AliMonitorDataTPC(10000);
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
			  AliMonitorHisto::kNormEntries);

  fTrackEta = CreateHisto1("TrackEta", "eta distribution of tracks", 
			   100, -2, 2, "#eta", "#Delta N/N",
			   AliMonitorHisto::kNormEntries);

  fTrackPhi = CreateHisto1("TrackPhi", "phi distribution of tracks", 
			   120, 0, 360, "#phi [#circ]", "#Delta N/N",
			   AliMonitorHisto::kNormEntries);
  fTrackPhi->SetDescription("The phi distribution should be flat on average.\nIf it is not flat check for dead TPC sectors.");

  fTrackNCl = CreateHisto1("TrackNCl", "Number of clusters per track", 
			   200, 0, 200, "N_{clusters}", "#Delta N/N",
			   AliMonitorHisto::kNormEntries);

  fTrackDEdxVsP = CreateHisto2("TrackDEdxVsP", "dE/dx of tracks", 
			       100, 0, 3, 100, 0, 200, 
			       "p [GeV/c]", "dE/dx", "#Delta N/N",
			       AliMonitorHisto::kNormEntries);
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

  fNTracks->Fill(tracks->GetEntries());
  fData->fNTracks = (Int_t) tracks->GetEntries();
  fData->SetSize(fData->fNTracks);
  for (Int_t i = 0; i < tracks->GetEntries(); i++) {
    tracks->GetEntry(i);
    fTrackPt->Fill(track->Pt());
    fTrackEta->Fill(track->Eta());
    fTrackPhi->Fill(track->Phi() * TMath::RadToDeg());
    fTrackNCl->Fill(track->GetNumberOfClusters());
    fTrackDEdxVsP->Fill(track->P(), track->GetdEdx());

    fData->fPt[i] = track->Pt();
    fData->fEta[i] = track->Eta();
    fData->fPhi[i] = track->Phi() * TMath::RadToDeg();
  }

  delete track;
  tpcLoader->UnloadTracks();
}
