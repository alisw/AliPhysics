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
#include "AliMonitorTrend.h"
#include "AliTPCParam.h"
#include <TFolder.h>
#ifdef ALI_HLT
#include <stdlib.h>
#include "AliL3MemHandler.h"
#include "AliL3SpacePointData.h"
#include "AliL3TrackArray.h"
#include "AliL3Track.h"
#include "AliL3Transform.h"
#endif

//_____________________________________________________________________________
AliMonitorHLT::AliMonitorHLT(AliTPCParam* param)
{
// create a HLT monitor object with the given parameters

  fParam = param;
}

//_____________________________________________________________________________
AliMonitorHLT::AliMonitorHLT(const AliMonitorHLT& monitor) :
  AliMonitor(monitor)
{
  Fatal("AliMonitorHLT", "copy constructor not implemented");
}

//_____________________________________________________________________________
AliMonitorHLT& AliMonitorHLT::operator = (const AliMonitorHLT& /*monitor*/)
{
  Fatal("operator =", "assignment operator not implemented");
  return *this;
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

  fTrackNHits = CreateHisto1("TrackNHits", "Number of hits per track", 
			   200, 0, 200, "N_{hits}", "#Delta N/N",
			   AliMonitorHisto::kNormEntries);

  fTrackDEdxVsP = CreateHisto2("TrackDEdxVsP", "dE/dx of tracks", 
			       100, 0, 3, 100, 0, 1000, 
			       "p [GeV/c]", "dE/dx", "#Delta N/N",
			       AliMonitorHisto::kNormEntries);

  fTrackDz0 = CreateHisto1("TrackDz0", "Dz0 of tracks", 
			   100, -100, 100, "#Delta z0 [cm]", "#Delta N/N",
			   AliMonitorHisto::kNormEntries);

  fTrackDr0 = CreateHisto1("TrackDr0", "Dr0 of tracks", 
			   130, 80, 250, "#Delta r0 [cm]", "#Delta N/N",
			   AliMonitorHisto::kNormEntries);

  fTrackAngle = CreateHisto1("TrackAngle", "azimutal distribution of tracks", 
			   100, -4, 4, " ", "#Delta N/N",
			   AliMonitorHisto::kNormEntries);

}


//_____________________________________________________________________________
void AliMonitorHLT::FillHistos(AliRunLoader* /*runLoader*/, 
			       AliRawReader* /*rawReader*/)
{
// fill the HLT monitor histogrms

#ifndef ALI_HLT
  Warning("FillHistos", "the code was compiled without HLT support");

#else
  AliL3MemHandler clusterHandler[36];
  AliL3SpacePointData *clusters[36];
  for (Int_t iSector = 0; iSector < fParam->GetNInnerSector(); iSector++) {
    char fileName[256];
    sprintf(fileName, "hlt/points_%d_-1.raw", iSector);
    if (!clusterHandler[iSector].SetBinaryInput(fileName)) {
      Warning("FillHistos", "could not open file %s", fileName);
      continue;
    }
    clusters[iSector] = (AliL3SpacePointData*) clusterHandler[iSector].Allocate();
    UInt_t nClusters = 0;
    clusterHandler[iSector].Binary2Memory(nClusters, clusters[iSector]);

    for (UInt_t iCluster = 0; iCluster < nClusters; iCluster++) {
      AliL3SpacePointData& cluster = clusters[iSector][iCluster];
      fClustersCharge->Fill(cluster.fCharge);
      fNClustersVsRow->Fill(cluster.fPadRow);
      fNClustersVsSector->Fill(iSector);
    }

    clusterHandler[iSector].CloseBinaryInput();
  }

  fNClustersVsSector->ScaleErrorBy(10.);

  AliL3MemHandler memHandler;
  if (!memHandler.SetBinaryInput("hlt/tracks.raw")) {
    Warning("FillHistos", "could not open file hlt/tracks.raw");
    return;
  }
  AliL3TrackArray* tracks = new AliL3TrackArray;
  memHandler.Binary2TrackArray(tracks);

  fNTracks->Fill(tracks->GetNTracks());
  for (Int_t iTrack = 0; iTrack < tracks->GetNTracks(); iTrack++) {
    AliL3Track* track = tracks->GetCheckedTrack(iTrack);
    if(!track) continue;
    track->CalculateHelix();
    fTrackPt->Fill(track->GetPt());
    fTrackEta->Fill(track->GetPseudoRapidity());
    fTrackPhi->Fill(track->GetPsi() * TMath::RadToDeg());
    fTrackDz0->Fill(track->GetZ0());
    fTrackDr0->Fill(track->GetR0());
    fTrackNHits->Fill(track->GetNHits());

    // Track dEdx
    Int_t nc=track->GetNHits();
    UInt_t *hits = track->GetHitNumbers();
    Float_t sampleDEdx[159];
    for (Int_t iHit = 0; iHit < nc; iHit++) {
      UInt_t hitID = hits[iHit];
      Int_t iSector = (hitID>>25) & 0x7f;
      UInt_t position = hitID&0x3fffff;
      UChar_t padrow = clusters[iSector][position].fPadRow;
      Float_t pWidth = AliL3Transform::GetPadPitchWidthLow();
      if (padrow>63)
	pWidth = AliL3Transform::GetPadPitchWidthUp(); 
      Float_t corr=1.; if (padrow>63) corr=0.67;
      sampleDEdx[iHit] = clusters[iSector][position].fCharge/pWidth*corr;
      Double_t crossingangle = track->GetCrossingAngle(padrow,iSector);
      fTrackAngle->Fill(crossingangle);
      Double_t s = sin(crossingangle);
      Double_t t = track->GetTgl();
      sampleDEdx[iHit] *= TMath::Sqrt((1-s*s)/(1+t*t));
    }

    /* Cook dEdx */
    Int_t i;
    Int_t swap;//stupid sorting
    do {
      swap=0;
      for (i=0; i<nc-1; i++) {
	if (sampleDEdx[i]<=sampleDEdx[i+1]) continue;
	Float_t tmp=sampleDEdx[i];
	sampleDEdx[i]=sampleDEdx[i+1]; sampleDEdx[i+1]=tmp;
	swap++;
      }
    } while (swap);

    Double_t low=0.05; Double_t up=0.7;
    Int_t nl=Int_t(low*nc), nu=Int_t(up*nc);
    Float_t trackDEdx=0;
    for (i=nl; i<=nu; i++) trackDEdx += sampleDEdx[i];
    trackDEdx /= (nu-nl+1);

    fTrackDEdxVsP->Fill(track->GetP(),trackDEdx);
  }

  delete tracks;
  memHandler.CloseBinaryInput();
#endif
}
