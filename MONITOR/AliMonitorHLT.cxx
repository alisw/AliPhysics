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
//  This class creates and fills the monitor histograms for the HLT          //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////


#include "AliMonitorHLT.h"
#ifdef ALI_HLT
#include <stdlib.h>
#include "AliL3MemHandler.h"
#include "AliL3SpacePointData.h"
#include "AliL3TrackArray.h"
#include "AliL3Track.h"
#endif

//_____________________________________________________________________________
AliMonitorHLT::AliMonitorHLT(AliTPCParam* param)
{
// create a HLT monitor object with the given parameters

  fParam = param;
}

//_____________________________________________________________________________
AliMonitorHLT::~AliMonitorHLT()
{
}


//_____________________________________________________________________________
void AliMonitorHLT::CreateHistos(TFolder* folder)
{
// create the HLT monitor histograms

  fFolder = folder->AddFolder("HLT", "HLT");

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
}


#include <TCanvas.h>
//_____________________________________________________________________________
void AliMonitorHLT::FillHistos(AliRunLoader* /*runLoader*/, 
			       AliRawReader* /*rawReader*/)
{
// fill the HLT monitor histogrms

#ifndef ALI_HLT
  Warning("FillHistos", "the code was compiled without HLT support");

#else
  AliL3MemHandler memHandler;
  for (Int_t iSector = 0; iSector < fParam->GetNInnerSector(); iSector++) {
    char fileName[256];
    sprintf(fileName, "hlt/points_%d_-1.raw", iSector);
    if (!memHandler.SetBinaryInput(fileName)) {
      Warning("FillHistos", "could not open file %s", fileName);
      continue;
    }
    AliL3SpacePointData* clusters = 
      (AliL3SpacePointData*) memHandler.Allocate();
    UInt_t nClusters = 0;
    memHandler.Binary2Memory(nClusters, clusters);

    for (UInt_t iCluster = 0; iCluster < nClusters; iCluster++) {
      AliL3SpacePointData& cluster = clusters[iCluster];
      fClustersCharge->Fill(cluster.fCharge);
      fNClustersVsRow->Fill(cluster.fPadRow);
      fNClustersVsSector->Fill(iSector);
    }

    memHandler.Free();
    memHandler.CloseBinaryInput();
  }

  fNClustersVsSector->ScaleErrorBy(10.);

  if (!memHandler.SetBinaryInput("hlt/tracks.raw")) {
    Warning("FillHistos", "could not open file hlt/tracks.raw");
    return;
  }
  AliL3TrackArray* tracks = new AliL3TrackArray;
  memHandler.Binary2TrackArray(tracks);

  fNTracks->Fill(tracks->GetNTracks());
  for (Int_t iTrack = 0; iTrack < tracks->GetNTracks(); iTrack++) {
    AliL3Track* track = tracks->GetTrack(iTrack);
    fTrackPt->Fill(track->GetPt());
    fTrackEta->Fill(track->GetPseudoRapidity());
    fTrackPhi->Fill(track->GetPsi() * TMath::RadToDeg());
  }

  delete tracks;
  memHandler.CloseBinaryInput();
#endif
}
