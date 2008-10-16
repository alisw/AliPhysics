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

/* $Id: AliTRDtrackingEfficiencyCombined.cxx 27496 2008-07-22 08:35:45Z cblume $ */

////////////////////////////////////////////////////////////////////////////
//                                                                        //
//  Reconstruction QA                                                     //
//                                                                        //
//  Authors:                                                              //
//    Markus Fasel <M.Fasel@gsi.de>                                       //
//                                                                        //
////////////////////////////////////////////////////////////////////////////

#include <TObjArray.h>
#include <TProfile.h>
#include <TMath.h>
#include <TCanvas.h>
#include "TTreeStream.h"

#include "AliMagFMaps.h"
#include "AliTracker.h"
#include "AliTrackReference.h"
#include "AliAnalysisManager.h"

#include "AliTRDseedV1.h"
#include "AliTRDtrackV1.h"
#include "AliTRDtrackInfo/AliTRDtrackInfo.h"
#include "AliTRDtrackingEfficiencyCombined.h"

ClassImp(AliTRDtrackingEfficiencyCombined)

//_____________________________________________________________________________
AliTRDtrackingEfficiencyCombined::AliTRDtrackingEfficiencyCombined()
  :AliTRDrecoTask("TrackingEffMC", "Combined Tracking Efficiency")
{
  //
  // Default constructor
  //
}


//_____________________________________________________________________________
void AliTRDtrackingEfficiencyCombined::CreateOutputObjects(){
  //
  // Create output objects
  //

  OpenFile(0, "RECREATE");

  const Int_t nbins = 11;
  Float_t xbins[nbins+1] = {.5, .7, .9, 1.3, 1.7, 2.4, 3.5, 4.5, 5.5, 7., 9., 11.};
  
  fContainer = new TObjArray();
  fContainer->Add(new TProfile("trEffComb", "Combined Tracking Efficiency", nbins, xbins));
  fContainer->Add(new TProfile("trContComb", "Combined Tracking Contamination", nbins, xbins));
}

//_____________________________________________________________________________
void AliTRDtrackingEfficiencyCombined::Exec(Option_t *){
  //
  // Do it
  //

	const Float_t kAlpha = 0.349065850;
	Int_t naccepted = 0, nrejected = 0, ndoublecounted = 0;
	Int_t labelsacc[10000];
	Int_t labelsrej[10000];
	Float_t momacc[10000];
	Float_t momrej[10000];
	TProfile *efficiency = (TProfile *)fContainer->At(0);
	TProfile *contamination = (TProfile *)fContainer->At(1);
	
	Int_t nTrackInfos = fTracks->GetEntriesFast();
	Double_t mom = 0;
	AliTRDtrackV1 *TRDtrack = 0x0;
	AliTRDtrackInfo *trkInf = 0x0;
	AliTrackReference *trackRef = 0x0;
	for(Int_t itinf = 0; itinf < nTrackInfos; itinf++){
		mom = 0.;
		trkInf = dynamic_cast<AliTRDtrackInfo *>(fTracks->UncheckedAt(itinf));
		if(!trkInf) continue;
		if((TRDtrack = trkInf->GetTrack()) || trkInf->GetNumberOfClustersRefit()){
			// check if allready found by the tracker
			Bool_t found = kFALSE;
			for(Int_t il = 0; il < naccepted; il++){
				if(labelsacc[il] == trkInf->GetLabel()) found = kTRUE;
			}
			if(found){
				mom =  trackRef ? trackRef->P() : TRDtrack->P();
				contamination->Fill(mom, 1);
				ndoublecounted++;
				continue;
			}
			if(trkInf->GetNTrackRefs()){
				Int_t iref = 0;
				while(!(trackRef = trkInf->GetTrackRef(iref++)));
			}
			if(!trackRef) printf("Error: Track Reference missing for Track %d\n", trkInf->GetLabel());
			mom =  trackRef ? trackRef->P() : trkInf->GetOuterParam()->P();

//           Accept track
			if(fDebugLevel > 3)printf("Accept track\n");
			momacc[naccepted] = mom; 
			labelsacc[naccepted++] = trkInf->GetLabel();
/*			printf("Reconstructed: event %3d Tracks: MC[%d] ESD[%d] NRefs[%d]\n", (Int_t)AliAnalysisManager::GetAnalysisManager()->GetCurrentEntry(), trkInf->GetLabel(),  trkInf->GetTrackId(), trkInf->GetNTrackRefs());*/
    } else{
			if(fDebugLevel>10) printf("Analysing Track References\n");
      // Check if track is findable
			Float_t xmin = 10000.0, xmax = 0.0; 
			Float_t ymin = 0.0, ymax = 0.0;
			Float_t zmin = 0.0, zmax = 0.0;
			Float_t lastx = 0.0, x = 0.0;
			Int_t nLayers = 0;
/*			trackRef = trkInf->GetTrackRef(0);*/
/*			xmin = trackRef->LocalX(); xmax = trackRef->LocalX();
      ymin = trackRef->LocalY(); ymax = trackRef->LocalY();
      mom = trackRef->P();*/
      Int_t *sector = new Int_t[trkInf->GetNTrackRefs()];
      for(Int_t itr = 0; itr < trkInf->GetNTrackRefs(); itr++){
        trackRef = trkInf->GetTrackRef(itr);
				if(fDebugLevel>10) printf("%d. x[%f], y[%f], z[%f]\n", itr, trackRef->LocalX(), trackRef->LocalY(), trackRef->Z());
				x = trackRef->LocalX(); 
				if(x < 250. || x > 370.) continue;	// Be Sure that we are inside TRD
        sector[itr] = Int_t(trackRef->Alpha()/kAlpha);
        if(x < xmin){
          xmin = trackRef->LocalX();
          ymin = trackRef->LocalY();
          zmin = trackRef->Z();
          mom = trackRef->P();
        }
        if(x > xmax){
          xmax = trackRef->LocalX();
          ymax = trackRef->LocalY();
          zmax = trackRef->Z();
        }
				if(itr > 0){
					Float_t dist = TMath::Abs(x - lastx);
					if(fDebugLevel>10) printf("x = %f, lastx = %f, dist = %f\n", x, lastx, dist);
					if(TMath::Abs(dist - 3.7) < 0.1) nLayers++; 	// ref(i+1) has to be larger than ref(i)
				}
				lastx = x;
      }
      // Apply cuts
      Bool_t findable = kTRUE;
      if(trkInf->GetNTrackRefs() > 2 && xmax > xmin){
        if(mom < 0.55) findable = kFALSE;									// momentum cut at 0.6
 				Double_t yangle = (ymax -ymin)/(xmax - xmin);
				Double_t zangle = (zmax -zmin)/(xmax - xmin);
				if(fDebugLevel>10) printf("track: y-Angle = %f, z-Angle = %f\n", yangle, zangle);
				if(fDebugLevel>10) printf("nLayers = %d\n", nLayers);
				if(TMath::ATan(TMath::Abs(yangle)) > 45.) findable = kFALSE;
				if(TMath::ATan(TMath::Abs(zangle)) > 45.) findable = kFALSE;
				if(nLayers < 4) findable = kFALSE;
        if(!trkInf->IsPrimary()) findable = kFALSE;
        Bool_t samesec = kTRUE;
        for(Int_t iref = 1; iref < trkInf->GetNTrackRefs(); iref++)
          if(sector[iref] != sector[0]) samesec = kFALSE;
        if(!samesec) findable = kFALSE;		// Discard all tracks which are passing more than one sector
        if(fDebugLevel){
          Double_t trackAngle = TMath::ATan(yangle);
          Bool_t primary = trkInf->IsPrimary();
          (*fDebugStream) << "NotFoundTrack"
            << "Momentum=" 	<< mom
            << "trackAngle="<< trackAngle
            << "NLayers="	<< nLayers
            << "Primary="	<< primary
            << "\n";
        }
      }
      else
        findable = kFALSE;
      delete[] sector;
      if(findable){
				momrej[nrejected] = mom;
				labelsrej[nrejected++] = trkInf->GetLabel();
/*			  printf("Not Reconstructed: event %3d Tracks: MC[%d] ESD[%d] NRefs[%d]\n", (Int_t)AliAnalysisManager::GetAnalysisManager()->GetCurrentEntry(), trkInf->GetLabel(),  trkInf->GetTrackId(), trkInf->GetNTrackRefs());*/
      }
    }
  }
	for(Int_t itk = 0; itk < naccepted; itk++){
		if(fDebugLevel > 2)printf("Accepted MC track: %d\n", labelsacc[itk]);
		efficiency->Fill(momacc[itk], 1);
		contamination->Fill(momacc[itk], 0);
	}
	Int_t nall = naccepted;
	for(Int_t imis = 0; imis < nrejected; imis++){
		Bool_t found = kFALSE;
		for(Int_t ifound = 0; ifound < naccepted; ifound++){
			if(labelsacc[ifound] == labelsrej[imis]){
				found = kTRUE;
				break;
			}
		}
		if(!found){
			efficiency->Fill(momrej[imis], 0);
			contamination->Fill(momrej[imis], 0);
			if(fDebugLevel > 2)printf("Rejected MC track: %d\n", labelsrej[imis]);
			nall++;
		}
	}
  //if(fDebugLevel>=1)
  printf("%3d Tracks: MC[%3d] TRD[%3d | %5.2f%%] \n", (Int_t)AliAnalysisManager::GetAnalysisManager()->GetCurrentEntry(), nall, naccepted, 1.E2*Float_t(naccepted)/Float_t(nall));
  printf("%3d Tracks: ALL[%3d] DoubleCounted[%3d | %5.2f%%] \n", (Int_t)AliAnalysisManager::GetAnalysisManager()->GetCurrentEntry(), nall + ndoublecounted, ndoublecounted, 1.E2*Float_t(ndoublecounted)/Float_t(nall + ndoublecounted));

  PostData(0, fContainer);
}

//_____________________________________________________________________________
void AliTRDtrackingEfficiencyCombined::Terminate(Option_t *)
{
  //
  // Termination
  //

  if(fDebugStream){ 
    delete fDebugStream;
    fDebugStream = 0x0;
    fDebugLevel = 0;
  }

  fContainer = dynamic_cast<TObjArray*>(GetOutputData(0));
  if (!fContainer) {
    Printf("ERROR: list not available");
    return;
  }

  /*TProfile *hEff = (TProfile*)fContainer->At(0);
  TProfile *hEffCont = (TProfile*)fContainer->At(1);
  Printf("Eff[%p] EffCont[%p]\n", (void*)hEff, (void*)hEffCont);


  TCanvas *c2 = new TCanvas("c2","",800,400);
  c2->Divide(2,1);

  c2->cd(1);
  hEff->DrawCopy("e1");
  c2->cd(2);
  hEffCont->DrawCopy("e1");*/
}
