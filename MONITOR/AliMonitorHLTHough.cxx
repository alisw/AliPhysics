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
// This class creates and fills monitor histograms for HLT Hough transform   //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////


#include "AliMonitorHLTHough.h"
#include "AliMonitorHisto.h"
#include "AliMonitorTrend.h"
#include "AliTPCParam.h"
#include <TFolder.h>
#ifdef ALI_HLT
#include <stdlib.h>
#include <AliL3MemHandler.h>
#include <AliL3TrackArray.h>
#include <AliL3HoughTrack.h>
#endif

//_____________________________________________________________________________
AliMonitorHLTHough::AliMonitorHLTHough(AliTPCParam* param)
{
// create a HLT monitor object with the given parameters

  fParam = param;
}

//_____________________________________________________________________________
AliMonitorHLTHough::AliMonitorHLTHough(const AliMonitorHLTHough& monitor) :
  AliMonitor(monitor)
{
  Fatal("AliMonitorHLTHough", "copy constructor not implemented");
}

//_____________________________________________________________________________
AliMonitorHLTHough& AliMonitorHLTHough::operator = (const AliMonitorHLTHough& 
						    /*monitor*/)
{
  Fatal("operator =", "assignment operator not implemented");
  return *this;
}

//_____________________________________________________________________________
AliMonitorHLTHough::~AliMonitorHLTHough()
{
}


//_____________________________________________________________________________
void AliMonitorHLTHough::CreateHistos(TFolder* folder)
{
// create the HLT Hough transform monitor histograms

  fFolder = folder->AddFolder("HLTHOUGH", "HLTHOUGH");

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
void AliMonitorHLTHough::FillHistos(AliRunLoader* /*runLoader*/, 
				    AliRawReader* /*rawReader*/)
{
// fill the HLT Hough transform monitor histograms

#ifndef ALI_HLT
  Warning("FillHistos", "the code was compiled without HLT support");

#else
  AliL3MemHandler memHandler[36];
  Int_t nHoughTracks = 0;
  for (Int_t iSector = 0; iSector < fParam->GetNInnerSector(); iSector++) {
    char fileName[256];
    sprintf(fileName, "hlt/tracks_ho_%d.raw", iSector);
    if (!memHandler[iSector].SetBinaryInput(fileName)) {
      Warning("FillHistos", "could not open file hlt/tracks.raw");
      return;
    }
    AliL3TrackArray* tracks = new AliL3TrackArray;
    memHandler[iSector].Binary2TrackArray(tracks);

    nHoughTracks += tracks->GetNTracks();
    for (Int_t iTrack = 0; iTrack < tracks->GetNTracks(); iTrack++) {
      AliL3HoughTrack *track = (AliL3HoughTrack*)tracks->GetCheckedTrack(iTrack);
      if(!track) continue;

      track->CalculateHelix();
      track->Rotate(iSector);

      fTrackPt->Fill(track->GetPt());
      fTrackEta->Fill(track->GetPseudoRapidity());
      fTrackPhi->Fill(track->GetPsi() * TMath::RadToDeg());
      if(track->GetPt()>3.) {
	fTrackEtaVsPhi->Fill(track->GetPseudoRapidity(),track->GetPsi() * TMath::RadToDeg());
	fPtEtaVsPhi->Fill(track->GetPseudoRapidity(),track->GetPsi() * TMath::RadToDeg(),track->GetPt());
      }
    }
    fNTracks->Fill(nHoughTracks);

    delete tracks;
    memHandler[iSector].CloseBinaryInput();
  }
#endif
}
