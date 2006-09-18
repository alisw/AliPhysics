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
//  This class creates and fills the monitor histograms for the ITS          //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include <TFolder.h>
#include <TTree.h>
#include <TVector3.h>

#include "AliLog.h"
#include "AliESD.h"
#include "AliITSRawStreamSDD.h"
#include "AliITSRawStreamSSD.h"
#include "AliITSclusterV2.h"
#include "AliITSgeom.h"
#include "AliMonitorHisto.h"
#include "AliMonitorITS.h"
#include "AliRawReader.h"
#include "AliRunLoader.h"

ClassImp(AliMonitorITS) 


//_____________________________________________________________________________
AliMonitorITS::AliMonitorITS(AliITSgeom* geom):
  AliMonitor(),
  fGeom(geom),
  fSDDDigitsCharge(NULL),
  fSSDDigitsCharge(NULL),
  fSDDClustersCharge(NULL),
  fSSDClustersCharge(NULL),
  fSPDNClustersVsModule(NULL),
  fSDDNClustersVsModule(NULL),
  fSSDNClustersVsModule(NULL),
  fNClustersVsLayer(NULL),
  fNTracks(NULL),
  fNTracksITSTPC(NULL),
  fTrackPt(NULL),
  fTrackEta(NULL),
  fTrackPhi(NULL),
  fTrackDEdxVsP(NULL)

{
// create a ITS monitor object with the given geometry

}


//_____________________________________________________________________________
void AliMonitorITS::CreateHistos(TFolder* folder)
{
// create the ITS monitor histograms

  fFolder = folder->AddFolder("ITS", "ITS");

  fSDDDigitsCharge = CreateHisto1("SDDDigitsCharge", 
				  "charge distribution of SDD digits", 
				  128, 0, 256, "charge", "#Delta N/N",
				  AliMonitorHisto::kNormEvents);

  fSSDDigitsCharge = CreateHisto1("SSDDigitsCharge", 
				  "charge distribution of SSD digits", 
				  100, 0, 200, "charge", "#Delta N/N",
				  AliMonitorHisto::kNormEvents);

  fSDDClustersCharge = CreateHisto1("SDDClustersCharge", 
				    "charge distribution of SDD clusters", 
				    100, 0, 200, "charge", "#Delta N/N",
				    AliMonitorHisto::kNormEvents);

  fSSDClustersCharge = CreateHisto1("SSDClustersCharge", 
				    "charge distribution of SSD clusters", 
				    100, 0, 400, "charge", "#Delta N/N",
				    AliMonitorHisto::kNormEvents);

  Int_t nModulesSPD = fGeom->GetLastSPD() - fGeom->GetStartSPD() + 1;
  fSPDNClustersVsModule = CreateHisto1("SPDNClustersVsModule", 
				       "mean number of clusters per SPD module", 
				       nModulesSPD, -0.5, nModulesSPD-0.5,
				       "module", "<N_{clusters}>",
				       AliMonitorHisto::kNormEvents);

  Int_t nModulesSDD = fGeom->GetLastSDD() - fGeom->GetStartSDD() + 1;
  fSDDNClustersVsModule = CreateHisto1("SDDNClustersVsModule", 
				       "mean number of clusters per SDD module", 
				       nModulesSDD, -0.5, nModulesSDD-0.5,
				       "module", "<N_{clusters}>",
				       AliMonitorHisto::kNormEvents);

  Int_t nModulesSSD = fGeom->GetLastSSD() - fGeom->GetStartSSD() + 1;
  fSSDNClustersVsModule = CreateHisto1("SSDNClustersVsModule", 
				       "mean number of clusters per SSD module", 
				       nModulesSSD, -0.5, nModulesSSD-0.5,
				       "module", "<N_{clusters}>",
				       AliMonitorHisto::kNormEvents);

  Int_t nLayer = fGeom->GetNlayers();
  fNClustersVsLayer = CreateHisto1("NClustersVsLayer", 
				   "mean number of clusters per layer", 
				   nLayer, 0.5, nLayer+0.5, 
				   "layer", "<N_{clusters}>",
				   AliMonitorHisto::kNormEvents);

  fNTracks = CreateHisto1("NTracks", "number of tracks per event", 
			  300, 0, 30000, "N_{tracks}", "#Delta N/N",
			  AliMonitorHisto::kNormEvents);

  fNTracksITSTPC = CreateHisto2("NTracksTPCITS", 
				"correlation of number of ITS and TPC tracks per event", 
				30, 0, 30000, 30, 0, 30000, 
				"N_{tracks}(ITS)", "N_{tracks}(TPC)", 
				"#Delta N/N", AliMonitorHisto::kNormEvents);

  fTrackPt = CreateHisto1("TrackPt", "pt distribution of tracks", 
			  90, 0, 3, "p_{t} [GeV/c]", "#Delta N/N",
			  AliMonitorHisto::kNormEntries);

  fTrackEta = CreateHisto1("TrackEta", "eta distribution of tracks", 
			   100, -2, 2, "#eta", "#Delta N/N",
			   AliMonitorHisto::kNormEntries);

  fTrackPhi = CreateHisto1("TrackPhi", "phi distribution of tracks", 
			   120, -180, 180, "#phi [#circ]", "#Delta N/N",
			   AliMonitorHisto::kNormEntries);

  fTrackDEdxVsP = CreateHisto2("TrackDEdxVsP", "dE/dx of tracks", 
			       100, 0, 3, 100, 0, 150, 
			       "p [GeV/c]", "dE/dx", "#Delta N/N",
			       AliMonitorHisto::kNormEntries);
}


//_____________________________________________________________________________
void AliMonitorITS::FillHistos(AliRunLoader* runLoader, 
			       AliRawReader* rawReader, AliESD* esd)
{
// fill the ITS monitor histogrms

  rawReader->Reset();
  AliITSRawStreamSDD inputSDD(rawReader);
  while (inputSDD.Next()) {
    fSDDDigitsCharge->Fill(inputSDD.GetSignal());
  }

  rawReader->Reset();
  AliITSRawStreamSSD inputSSD(rawReader);
  while (inputSSD.Next()) {
    fSSDDigitsCharge->Fill(inputSSD.GetSignal());
  }

  AliLoader* itsLoader = runLoader->GetLoader("ITSLoader");
  if (!itsLoader) return;


  itsLoader->LoadRecPoints();
  TTree* clustersTree = itsLoader->TreeR();
  if (!clustersTree) return;
  TClonesArray* clusters = new TClonesArray("AliITSclusterV2");
  clustersTree->SetBranchAddress("Clusters", &clusters);

  for (Int_t iModule = 0; iModule < clustersTree->GetEntries(); iModule++) {
    clustersTree->GetEntry(iModule);
    Int_t iLayer, iLadder, iDetector;
    fGeom->GetModuleId(iModule, iLayer, iLadder, iDetector);

    for (Int_t j = 0; j < clusters->GetEntriesFast(); j++) {
      AliITSclusterV2* cluster = (AliITSclusterV2*) clusters->At(j);
      switch (iLayer) {
      case 1: case 2: {
	fSPDNClustersVsModule->Fill(iModule - fGeom->GetStartSPD());
	break;
      }
      case 3: case 4: {
	fSDDClustersCharge->Fill(cluster->GetQ());
	fSDDNClustersVsModule->Fill(iModule - fGeom->GetStartSDD());
	break;
      }
      case 5: case 6: {
	fSSDClustersCharge->Fill(cluster->GetQ());
	fSSDNClustersVsModule->Fill(iModule - fGeom->GetStartSSD());
	break;
      }
      }
      fNClustersVsLayer->Fill(iLayer);
    }
  }

  delete clusters;
  itsLoader->UnloadRecPoints();


  Int_t nTracks = 0;
  Int_t nTPCTracks = 0;
  for (Int_t i = 0; i < esd->GetNumberOfTracks(); i++) {
    AliESDtrack* track = esd->GetTrack(i);
    if (!track) continue;
    if ((track->GetStatus() | AliESDtrack::kTPCin) != 0) nTPCTracks++;
    if ((track->GetStatus() | AliESDtrack::kITSin) == 0) continue;
    nTracks++;

    Double_t pxyz[3];
    track->GetPxPyPz(pxyz);
    TVector3 pTrack(pxyz);
    fTrackPt->Fill(pTrack.Pt());
    fTrackEta->Fill(pTrack.Eta());
    fTrackPhi->Fill(pTrack.Phi() * TMath::RadToDeg());
    fTrackDEdxVsP->Fill(pTrack.Mag(), track->GetITSsignal());
  }
  fNTracks->Fill(nTracks);
  fNTracksITSTPC->Fill(nTracks, nTPCTracks);
}
