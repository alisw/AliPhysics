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

/* $Id: AliTRDtrackingEfficiency.cxx 27496 2008-07-22 08:35:45Z cblume $ */

////////////////////////////////////////////////////////////////////////////
//                                                                        //
//  Reconstruction QA                                                     //
//                                                                        //
//  Authors:                                                              //
//    Markus Fasel <M.Fasel@gsi.de>                                       //
//                                                                        //
////////////////////////////////////////////////////////////////////////////

#include <TClonesArray.h>
#include <TObjArray.h>
#include <TProfile.h>
#include "TTreeStream.h"

#include "AliPID.h"
#include "AliESDtrack.h"
#include "AliTrackReference.h"
#include "AliExternalTrackParam.h"
#include "AliTracker.h"
#include "AliMagF.h"
#include "AliAnalysisManager.h"

#include "Cal/AliTRDCalPID.h"
#include "AliTRDtrackingEfficiency.h"
#include "AliTRDtrackInfo/AliTRDtrackInfo.h"

ClassImp(AliTRDtrackingEfficiency)

//____________________________________________________________________
AliTRDtrackingEfficiency::AliTRDtrackingEfficiency()
  :AliTRDrecoTask("TrackingEff", "Barrel Tracking Efficiency")
  ,fMissed(0x0)
{
  //
  // Default constructor
  //
}

//____________________________________________________________________
AliTRDtrackingEfficiency::~AliTRDtrackingEfficiency()
{
  if(fMissed){
    fMissed->Delete();
    delete fMissed;
  }
}

//____________________________________________________________________
void  AliTRDtrackingEfficiency::CreateOutputObjects()
{
  //
  // Create output objects
  //

  OpenFile(0, "RECREATE");
  const Int_t nbins = AliTRDCalPID::kNMom;
  Float_t xbins[nbins+1] = {.5, .7, .9, 1.3, 1.7, 2.4, 3.5, 4.5, 5.5, 7., 9., 11.};

  TH1 *h = 0x0;
  fContainer = new TObjArray();
  for(Int_t is=0; is<AliPID::kSPECIES; is++){
    fContainer->Add(h = new TProfile(Form("h%s", AliTRDCalPID::GetPartSymb(is)), "", nbins, xbins));
    h->SetLineColor(AliTRDCalPID::GetPartColor(is));
    h->SetMarkerColor(AliTRDCalPID::GetPartColor(is));
    h->SetMarkerStyle(7);
  }
  fContainer->Add(h = new TProfile("h", "", nbins, xbins));
  h->SetMarkerStyle(7);
} 

//____________________________________________________________________
void AliTRDtrackingEfficiency::Exec(Option_t *)
{
  //
  // Do it
  //

  Int_t labelsacc[10000]; 
  memset(labelsacc, 0, sizeof(Int_t) * 10000);
	
  if(!fMissed){ 
    fMissed = new TClonesArray("AliTRDtrackInfo", 10);
    fMissed->SetOwner();
  }

  Float_t mom;
  Int_t selection[10000], nselect = 0;
  ULong_t status; Int_t pidx;
  Int_t nTRD = 0, nTPC = 0, nMiss = 0;
  AliTRDtrackInfo     *track = 0x0;
  AliTrackReference     *ref = 0x0;
  AliExternalTrackParam *esd = 0x0;
  for(Int_t itrk=0; itrk<fTracks->GetEntriesFast(); itrk++){
    track = (AliTRDtrackInfo*)fTracks->UncheckedAt(itrk);

		if(!track->HasESDtrack()) continue;
    status = track->GetStatus();

    // missing TPC propagation - interesting for SA
    if(!(status&AliESDtrack::kTPCout)) continue;

    // missing MC info.
    if(HasMCdata() && track->GetNTrackRefs() <= 1) continue;
   
    nTPC++;
    selection[nselect++]=itrk;
    ref  = track->GetTrackRef(0);
    esd  = track->GetESDinfo()->GetOuterParam();
    mom  = ref ? ref->P(): esd->P();
    pidx = AliTRDCalPID::GetPartIndex(track->GetPDG());
    pidx = TMath::Max(pidx, 0);

    //Int_t n = track->GetNumberOfClusters(); 
    // where are this tracklets ???
    //if(ncls0 > ncls1) printf("%3d ESD[%3d] TRD[%3d|%3d]\n", itrk, ncls0, ncls1, n);
    if(track->GetNumberOfClustersRefit()){ 
      ((TProfile*)fContainer->At(pidx))->Fill(mom, 1.);
			labelsacc[nTRD] = track->GetLabel();
      nTRD++;
      continue;
    }



    Float_t xmed, xleng;
    Int_t iref = 1; Bool_t found = kFALSE;
    while((ref = track->GetTrackRef(iref))){
      xmed = .5*(ref->LocalX() + track->GetTrackRef(iref-1)->LocalX());
      xleng= (ref->LocalX() - track->GetTrackRef(iref-1)->LocalX());
      if(TMath::Abs(xmed - 298.5) < .5 &&
        TMath::Abs(xleng - 3.7) < .1){ 
        found = kTRUE;
        break;
      }
      iref++;
    }
    if(!found){ 
      nTPC--;
      // track missing first layer. Maybe interesting for SA.
      continue;
    }
    nselect--;
    new ((*fMissed)[nMiss]) AliTRDtrackInfo(*track);
    nMiss++;
  }
  if(fDebugLevel>=2) printf("%3d Tracks: ESD[%3d] TPC[%3d] TRD[%3d | %5.2f%%] Off[%d]\n", (Int_t)AliAnalysisManager::GetAnalysisManager()->GetCurrentEntry(), fTracks->GetEntriesFast(), nTPC, nTRD, nTPC ? 1.E2*nTRD/float(nTPC) : 0., fMissed->GetEntriesFast());


  // Find double tracks
  Float_t threshold = 10.;
  AliTrackReference *refMiss = 0x0;
  AliExternalTrackParam *op = 0x0;
  AliTRDtrackInfo       *tt = 0x0;
  for(Int_t imiss=0; imiss<nMiss; imiss++){
    //printf("Searching missing %d ...\n", imiss);

    // get outer param of missed
    tt = (AliTRDtrackInfo*)fMissed->UncheckedAt(imiss);
    op = tt->GetESDinfo()->GetOuterParam();
    Double_t alpha = op->GetAlpha(), cosa = TMath::Cos(alpha), sina = TMath::Sin(alpha);

    Double_t xyz[3], x0, y0, z0, x, y, z, dx, dy, dz, d;

    Bool_t FOUND = kFALSE;
    for(Int_t iselect=0; iselect<nselect; iselect++){
      track = (AliTRDtrackInfo*)fTracks->UncheckedAt(selection[iselect]);

      // check first MC ... if available
      d = 0;
      for(Int_t iref=0; iref<track->GetNTrackRefs(); iref++){
        if(!(ref = track->GetTrackRef(iref))) continue;
        if((refMiss = tt->GetTrackRef(iref))){
          dy = ref->LocalY() - refMiss->LocalY();
          dz = ref->Z() - refMiss->Z();
        } else {
          // compare missOP with refTrackRef in LTC
          x0 = ref->LocalX();
          op->GetYAt(x0, AliTracker::GetBz(), y0);
          op->GetZAt(x0, AliTracker::GetBz(), z0);
          dy = y0 - ref->LocalY();
          dz = z0 - ref->Z();
        }
        d += (dy*dy + dz*dz);
      }
      //printf("\td[%d] = %f N[%d]\n", selection[iselect], d, track->GetNTrackRefs());
      if((track->GetNTrackRefs())){ 
        d /= track->GetNTrackRefs();
        if(d < threshold){
          //printf("\t\tFound %2d in ref[%3d] : d[%f]\n", imiss, selection[iselect], d/track->GetNTrackRefs());
          FOUND = kTRUE; break;
        }
      }

      // process outer param ... always available
      // compare missOP with OP in GTC
      esd = track->GetESDinfo()->GetOuterParam();
      esd->GetXYZ(xyz);
      x0 = esd->GetX();
      op->GetYAt(x0, AliTracker::GetBz(), y0);
      op->GetZAt(x0, AliTracker::GetBz(), z0);
      x = x0*cosa - y0*sina;
      y = x0*sina + y0*cosa;
      z = z0;
      dx=xyz[0]-x;
      dy=xyz[1]-y;
      dz=xyz[2]-z;
      d = dx*dx+dy*dy+dz*dz;
      //printf("\td[%d] = %f op\n", selection[iselect], d);
      if(d < threshold){
        //printf("\t\tFound %2d in op[%3d]  : d[%f] dx[%5.2f] dy[%5.2f] dz[%5.2f]\n", imiss, selection[iselect], d, dx, dy, dz);
        FOUND = kTRUE; break;
      }
    }
    if(FOUND) nTPC--;
    else{ 
      ref = tt->GetTrackRef(0);
      mom = ref ? ref->P(): op->P();
      pidx = AliTRDCalPID::GetPartIndex(tt->GetPDG());
      pidx = TMath::Max(pidx, 0);
      ((TProfile*)fContainer->At(pidx))->Fill(mom, 0.);
      if(fDebugLevel>=2) printf("\tNOT FOUND Id[%d] Mom[%f]\n", tt->GetTrackId(), mom);
    }
  }

  if(fDebugLevel>=2) printf("%3d Tracks: ESD[%3d] TPC[%3d] TRD[%3d | %5.2f%%] Off[%d]\n", (Int_t)AliAnalysisManager::GetAnalysisManager()->GetCurrentEntry(), fTracks->GetEntriesFast(), nTPC, nTRD, nTPC ? 1.E2*nTRD/float(nTPC) : 0., fMissed->GetEntriesFast());

  //fMissed->Delete();
	// check for double countings
	Int_t indices[10000]; memset(indices, 0, sizeof(Int_t) * 10000);
	TMath::Sort(nTRD, labelsacc, indices);
	if(fDebugLevel > 2){
	for(Int_t itk = 0; itk < nTRD - 1; itk++)
		if(labelsacc[indices[itk]] ==labelsacc[indices[itk + 1]]) printf("Double counted MC track: %d\n", labelsacc[indices[itk]]);
	}
  PostData(0, fContainer);
}

//____________________________________________________________________
void AliTRDtrackingEfficiency::Terminate(Option_t *)
{
  //
  // Terminate
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

}



//____________________________________________________________________
Bool_t AliTRDtrackingEfficiency::GetRefFigure(Int_t ifig)
{
  Bool_t FIRST = kTRUE;
  TProfile *h = 0x0;
  switch(ifig){
  case 0:
    h = (TProfile*)fContainer->At(AliPID::kSPECIES);
    for(Int_t is=0; is<AliPID::kSPECIES; is++){
      h->Add((TProfile*)fContainer->At(is));
    }
    h->SetMarkerColor(kBlack);
    h->SetLineColor(kBlack);
    h->GetXaxis()->SetTitle("p [GeV/c]");
    h->Draw("e1");
    break;
  case 1:
    FIRST = kTRUE;
    for(Int_t is=0; is<AliPID::kSPECIES; is++){
      if(!(h = (TProfile*)fContainer->At(is))) continue;
      if(FIRST){
        h->Draw("e1");
        h->GetXaxis()->SetTitle("p [GeV/c]");
      } else h->Draw("same e1");
      FIRST = kFALSE;
    }
    break;
  }
  return kTRUE;
}


//____________________________________________________________________
Bool_t AliTRDtrackingEfficiency::PostProcess()
{
  fNRefFigures = HasMCdata() ? 2 : 1; 
  return kTRUE;
}
