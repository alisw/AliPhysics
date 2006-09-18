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

#include <TFolder.h>
#include <TTree.h>
#include <TVector3.h>

#include "AliLog.h"
#include "AliESD.h"
#include "AliMonitorDataTPC.h"
#include "AliMonitorHisto.h"
#include "AliMonitorTPC.h"
#include "AliMonitorTrend.h"
#include "AliRawReader.h"
#include "AliRunLoader.h"
#include "AliTPCClustersRow.h"
#include "AliTPCParam.h"
#include "AliTPCRawStream.h"
#include "AliTPCclusterMI.h"

ClassImp(AliMonitorTPC) 


//_____________________________________________________________________________
AliMonitorTPC::AliMonitorTPC(AliTPCParam* param):
  AliMonitor(),
  fParam(param),
  fPadsCharge(NULL),
  fClustersCharge(NULL),
  fNClustersVsRow(NULL),
  fNClustersVsSector(NULL),
  fNTracks(NULL),
  fTrackPt(NULL),
  fTrackEta(NULL),
  fTrackPhi(NULL),
  fTrackNCl(NULL),
  fTrackDEdxVsP(NULL),
  fTrackDEdx(NULL),
  fTrackEtaVsPhi(NULL),
  fPtEtaVsPhi(NULL),
  fData(new AliMonitorDataTPC(10000))
{
// create a TPC monitor object with the given parameters

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
			   120, -180, 180, "#phi [#circ]", "#Delta N/N",
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
			       AliRawReader* rawReader, AliESD* esd)
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


  Int_t nTracks = 0;
  for (Int_t i = 0; i < esd->GetNumberOfTracks(); i++) {
    AliESDtrack* track = esd->GetTrack(i);
    if (!track || ((track->GetStatus() | AliESDtrack::kTPCin) == 0)) continue;
    nTracks++;

    Double_t pxyz[3];
    track->GetInnerPxPyPz(pxyz);
    TVector3 pTrack(pxyz);
    Double_t p = pTrack.Mag();
    Double_t pt = pTrack.Pt();
    Double_t eta = pTrack.Eta();
    Double_t phi = pTrack.Phi() * TMath::RadToDeg();

    fTrackPt->Fill(pt);
    fTrackEta->Fill(eta);
    fTrackPhi->Fill(phi);
    if (pt > 3.) {
      fTrackEtaVsPhi->Fill(eta, phi);
      fPtEtaVsPhi->Fill(eta, phi, pTrack.Pt());
    }
    fTrackNCl->Fill(track->GetTPCclusters(NULL));
    fTrackDEdxVsP->Fill(p, track->GetTPCsignal());
    if(p>0.4 && p<1.0)
      fTrackDEdx->Fill(track->GetTPCsignal());

    fData->SetData(i, pt, eta, phi); 
  }
  fNTracks->Fill(nTracks);
  fData->SetNTracks(nTracks);
}
