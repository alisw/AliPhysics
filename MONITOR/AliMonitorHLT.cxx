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
#include "AliLog.h"
#include <TFolder.h>
#include <stdlib.h>
#include "AliHLTMemHandler.h"
#include "AliHLTSpacePointData.h"
#include "AliHLTTrackArray.h"
#include "AliHLTTrack.h"
#include "AliHLTTransform.h"
#include "AliHLTVertex.h"

//_____________________________________________________________________________
AliMonitorHLT::AliMonitorHLT(AliTPCParam* param):
  AliMonitor(),
  fParam(param),
  fClustersCharge(NULL),
  fNClustersVsRow(NULL),
  fNClustersVsSector(NULL),
  fNTracks(NULL),
  fTrackPt(NULL),
  fTrackEta(NULL),
  fTrackPhi(NULL),
  fTrackNHits(NULL),
  fTrackDEdxVsP(NULL),
  fTrackDEdx(NULL),
  fTrackDz0(NULL),
  fTrackDr0(NULL),
  fTrackEtaVsPhi(NULL),
  fPtEtaVsPhi(NULL),
  fTrackZvsNHits(NULL),
  fTrackXYvsNHits(NULL)
{
// create a HLT monitor object with the given parameters

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
			  AliMonitorHisto::kNormNone);

  fTrackEta = CreateHisto1("TrackEta", "eta distribution of tracks", 
			   100, -2, 2, "#eta", "#Delta N/N",
			   AliMonitorHisto::kNormEntries);

  fTrackPhi = CreateHisto1("TrackPhi", "phi distribution of tracks", 
			   120, 0, 360, "#phi [#circ]", "#Delta N/N",
			   AliMonitorHisto::kNormEntries);

  fTrackNHits = CreateHisto1("TrackNHits", "Number of hits per track", 
			   200, 0, 200, "N_{hits}", "#Delta N/N",
			   AliMonitorHisto::kNormNone);

  fTrackXYvsNHits = CreateHisto2("TrackXYvsNHits", "XY vs Number of hits per track", 
				 50, 0, 200,50,0,10,
				 "N_{hits}","Radius [cm]","#Delta N/N",
				 AliMonitorHisto::kNormNone);

  fTrackZvsNHits = CreateHisto2("TrackZvsNHits", "Z vs Number of hits per track", 
				 50, 0, 200,50,-20,20,
				 "N_{hits}","Z [cm]","#Delta N/N",
				 AliMonitorHisto::kNormNone);

  fTrackDEdxVsP = CreateHisto2("TrackDEdxVsP", "dE/dx of tracks", 
			       100, 0, 3, 100, 0, 1000, 
			       "p [GeV/c]", "dE/dx", "#Delta N/N",
			       AliMonitorHisto::kNormEntries);

  fTrackDEdx = CreateHisto1("TrackDEdx", "dE/dx for tracks with 0.4<p<1.0 GeV/c", 
			       50, 0, 300, 
			       "dE/dx", "#Delta N/N",
			       AliMonitorHisto::kNormEntries);

  fTrackDz0 = CreateHisto1("TrackDz0", "Dz0 of tracks", 
			   100, -100, 100, "#Delta z0 [cm]", "#Delta N/N",
			   AliMonitorHisto::kNormEntries);

  fTrackDr0 = CreateHisto1("TrackDr0", "Dr0 of tracks", 
			   130, 80, 250, "#Delta r0 [cm]", "#Delta N/N",
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
void AliMonitorHLT::FillHistos(AliRunLoader* /*runLoader*/, 
			       AliRawReader* /*rawReader*/, AliESD* /*esd*/)
{
// fill the HLT monitor histogrms

  AliHLTMemHandler clusterHandler[36];
  AliHLTSpacePointData *clusters[36];
  for (Int_t iSector = 0; iSector < fParam->GetNInnerSector(); iSector++) {
    char fileName[256];
    sprintf(fileName, "hlt/points_%d_-1.raw", iSector);
    if (!clusterHandler[iSector].SetBinaryInput(fileName)) {
      AliWarning(Form("could not open file %s", fileName));
      continue;
    }
    clusters[iSector] = (AliHLTSpacePointData*) clusterHandler[iSector].Allocate();
    UInt_t nClusters = 0;
    clusterHandler[iSector].Binary2Memory(nClusters, clusters[iSector]);

    for (UInt_t iCluster = 0; iCluster < nClusters; iCluster++) {
      AliHLTSpacePointData& cluster = clusters[iSector][iCluster];
      fClustersCharge->Fill(cluster.fCharge);
      fNClustersVsRow->Fill(cluster.fPadRow);
      fNClustersVsSector->Fill(iSector);
    }

    clusterHandler[iSector].CloseBinaryInput();
  }

  fNClustersVsSector->ScaleErrorBy(10.);

  AliHLTMemHandler memHandler;
  if (!memHandler.SetBinaryInput("hlt/tracks.raw")) {
    AliWarning("could not open file hlt/tracks.raw");
    return;
  }
  AliHLTTrackArray* tracks = new AliHLTTrackArray;
  memHandler.Binary2TrackArray(tracks);
  Double_t xc,yc,zc;
  AliHLTVertex vertex;

  fNTracks->Fill(tracks->GetNTracks());
  for (Int_t iTrack = 0; iTrack < tracks->GetNTracks(); iTrack++) {
    AliHLTTrack* track = tracks->GetCheckedTrack(iTrack);
    if(!track) continue;
    track->CalculateHelix();
    track->GetClosestPoint(&vertex,xc,yc,zc);
    if(TMath::Abs(zc)>10.) continue;
    fTrackPt->Fill(track->GetPt());
    fTrackEta->Fill(track->GetPseudoRapidity());
    fTrackPhi->Fill(track->GetPsi() * TMath::RadToDeg());
    if(track->GetPt()>3.) {
      fTrackEtaVsPhi->Fill(track->GetPseudoRapidity(),track->GetPsi() * TMath::RadToDeg());
      fPtEtaVsPhi->Fill(track->GetPseudoRapidity(),track->GetPsi() * TMath::RadToDeg(),track->GetPt());
    }
    fTrackDz0->Fill(track->GetZ0());
    fTrackDr0->Fill(track->GetR0());
    fTrackNHits->Fill(track->GetNHits());
    fTrackXYvsNHits->Fill(track->GetNHits(),TMath::Sqrt(xc*xc+yc*yc));
    fTrackZvsNHits->Fill(track->GetNHits(),zc);

    // Track dEdx
    Int_t nc=track->GetNHits();
    UInt_t *hits = track->GetHitNumbers();
    Float_t sampleDEdx[159];
    for (Int_t iHit = 0; iHit < nc; iHit++) {
      UInt_t hitID = hits[iHit];
      Int_t iSector = (hitID>>25) & 0x7f;
      UInt_t position = hitID&0x3fffff;
      UChar_t padrow = clusters[iSector][position].fPadRow;
      Float_t pWidth = AliHLTTransform::GetPadPitchWidthLow();
      if (padrow>63)
	pWidth = AliHLTTransform::GetPadPitchWidthUp(); 
      Float_t corr=1.; if (padrow>63) corr=0.67;
      sampleDEdx[iHit] = clusters[iSector][position].fCharge/pWidth*corr;
      Double_t crossingangle = track->GetCrossingAngle(padrow,iSector);
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
    if(track->GetP()>0.4 && track->GetP()<1.0)
      fTrackDEdx->Fill(trackDEdx);
  }

  delete tracks;
  memHandler.CloseBinaryInput();
}
