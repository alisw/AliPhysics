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


////////////////////////////////////////////////////////////////////////////
//                                                                        //
//  TRD multiplicity                                                      //
//                                                                        //
// Task to select true single tracks
// 
// Program flow:
// Part A - task AliTRDmultiplicity (within the framework qaRec):
//   For TRD standalone or TRD&TOF tracks I make a FitRiemanTilt
//   The resulting parameterization is dumped into a DebugStreamer written to a file
// 
// Part B – $ALICE_ROOT/TRD/qaRec/macros/TrackletsinTRD.C[h] (analysis macro):
//   The clusters are read in
//   The above produced file is read in
//   I make a fit through the parameterization of the  FitRiemanTilt -> in order to get a straight line (representing the track)
//   the distance of each cluster with respect to the straight line is calculated
//   If there is no cluster at a distance of 0.6  to 2 with respect to the track, the track is defined as a good track, i.e. clean single track
//                                                                        //
//  Authors: Yvonne C. Pachmayer <pachmay@physi.uni-heidelberg.de>        //
//                                                                        //
////////////////////////////////////////////////////////////////////////////

#include <TClonesArray.h>
#include <TObjArray.h>
#include <TProfile.h>
#include <TH1.h>
#include <TH1F.h>
#include <TFile.h>
#include "AliTRDtracker.h"
#include "TTreeStream.h"

#include "AliPID.h"
#include "AliESDtrack.h"
#include "AliTrackReference.h"
#include "AliExternalTrackParam.h"
#include "AliTracker.h"
#include "AliTRDseedV1.h"
#include "AliTRDtrackV1.h"
#include "AliTRDtrackerV1.h"
#include "AliTrackPointArray.h"
//#include "AliMagFMaps.h"
#include "AliTRDcluster.h"
#include "AliAnalysisManager.h"

#include "Cal/AliTRDCalPID.h"
#include "AliTRDgeometry.h"
#include "AliTRDmultiplicity.h"
#include "AliTRDtrackInfo/AliTRDtrackInfo.h"

ClassImp(AliTRDmultiplicity)

//____________________________________________________________________
AliTRDmultiplicity::AliTRDmultiplicity()
    :AliTRDrecoTask("Multiplicity", "Barrel Tracking Multiplicity")
  ,fEventCounter(0)
{
  //
  // Default constructor
    //
}

//____________________________________________________________________
AliTRDmultiplicity::~AliTRDmultiplicity()
{
}

//____________________________________________________________________
void  AliTRDmultiplicity::CreateOutputObjects()
{
  //
  // Create output objects
  //

  OpenFile(0, "RECREATE");

  TH1 *h = 0x0;
  fContainer = new TObjArray();
  for(Int_t is=0; is<AliTRDgeometry::kNsector; is++){
  fContainer->Add(h = new TH1F(Form("h_sector_%i", is), " ", 100,-10,10));
  }
} 

//____________________________________________________________________
void AliTRDmultiplicity::Exec(Option_t *)
{
  //
  // Do it
  //



  ULong_t status;
  AliTRDtrackInfo     *track = 0x0;
  AliExternalTrackParam *esd = 0x0;
  AliTRDtrackV1 *TRDtrack = 0x0;
  Double_t x_anode[AliTRDgeometry::kNlayer] = {300.2, 312.8, 325.4, 338.0, 350.6, 363.2}; // Take the default X0
  Int_t stand_alone=1;
  fEventCounter++;


  for(Int_t itrk=0; itrk<fTracks->GetEntriesFast(); itrk++){
    track = (AliTRDtrackInfo*)fTracks->UncheckedAt(itrk);
    if(!track->HasESDtrack()) continue;
    status = track->GetStatus();

    AliTrackPoint points[6];
    Float_t xyz[3];
    memset(xyz, 0, sizeof(Float_t) * 3);
    Float_t tracklet_x[6];
    Float_t tracklet_y[6];
    Float_t tracklet_z[6];
    for(Int_t a=0;a<6;a++)
    {
        xyz[0] = x_anode[a];
        points[a].SetXYZ(xyz);
        tracklet_x[a]=-1000;
        tracklet_y[a]=-1000;
        tracklet_z[a]=-1000;
    }
    Int_t det_tracklet=600;



    // TRD standalone
    if(((status&AliESDtrack::kTRDout)>0) && !((status&AliESDtrack::kTRDin)>0))
    {
        TRDtrack = track->GetTrack();

        if(TRDtrack)
        {
            AliTRDseedV1 *tracklet = 0x0;
            Int_t counter_tracklet=0;
            for(Int_t itl = 0; itl < AliTRDgeometry::kNlayer; itl++)
            {
                tracklet = TRDtrack->GetTracklet(itl);
                if(tracklet)
                {
                    if(tracklet->IsOK())
                    {
                        counter_tracklet++;
                        det_tracklet=tracklet->GetDetector();
                        //printf("%i %f %f \n",itl,tracklet->GetYfit(0),tracklet->GetZfit(0));
                        tracklet_x[itl]=tracklet->GetX0();
                        tracklet_y[itl]=tracklet->GetYfit(0);
                        tracklet_z[itl]=tracklet->GetZfit(0);
                    }
                }
            }
            // this cut is needed otherwise the order of tracklets in the fit function is not correct
            if(counter_tracklet==AliTRDgeometry::kNlayer) AliTRDtrackerV1::FitRiemanTilt(const_cast<AliTRDtrackV1 *>(TRDtrack), 0x0, kTRUE, counter_tracklet, points);
            else continue;

            if(fDebugLevel>=1){
              for(Int_t il = 0; il < AliTRDgeometry::kNlayer; il++){
                //printf("---------------------------------------- %i %i %f %f %f \n",counter_tracklet,il,points[il].GetX(),points[il].GetY(),points[il].GetZ());
                Double_t point_x=points[il].GetX();
                Double_t point_y=points[il].GetY();
                Double_t point_z=points[il].GetZ();


                (*fDebugStream)   << "TrackletsinTRD"
                    << "standalone=" << stand_alone
                    << "eventcounter=" << fEventCounter
                    << "layer="     << il
                    << "dettracklet=" << det_tracklet
                    << "xtracklet=" << tracklet_x[il]
                    << "xtrack="    << point_x
                    << "ytracklet=" << tracklet_y[il]
                    << "ytrack="    << point_y
                    << "ztracklet=" << tracklet_z[il]
                    << "ztrack="    << point_z
                    << "num_tracklets=" << counter_tracklet
                    << "\n";

              }
            }
        }
    } // end TRD standalone



    // TPC cluster selection for cosmic data
    if((track->GetTPCncls())<40) continue;
    // missing TPC propagation
    if(!(status&AliESDtrack::kTPCout)) continue;
    stand_alone=0;

    esd = track->GetESDinfo()->GetOuterParam();
    // calculate sector - does not matter if track ends in TPC or TRD - sector number independent
    Double_t alpha;
    alpha=-100;
    Int_t fSector=100;
    if(esd){
        alpha=esd->GetAlpha();
        fSector = (Int_t)((alpha+ TMath::Pi())/(20*TMath::Pi()/180));
        if((fSector>-1) && (fSector<18))
        {
        }
        else
        {
            continue;
        }
    }



    // only these sectors have a TRD detector at the moment
   if(fSector==0||fSector==8||fSector==9||fSector==17)
    {
        TRDtrack = track->GetTrack();
        if(!TRDtrack) continue;
        AliTRDseedV1 *tracklet = 0x0;
        Int_t counter_tracklet=0;


        for(Int_t itl = 0; itl < AliTRDgeometry::kNlayer; itl++){
            tracklet = TRDtrack->GetTracklet(itl);
            if(!tracklet || !tracklet->IsOK()) continue;
            counter_tracklet++;
            det_tracklet=tracklet->GetDetector();
            tracklet_x[itl]=tracklet->GetX0();
            tracklet_y[itl]=tracklet->GetYfit(0);
            tracklet_z[itl]=tracklet->GetZfit(0);

            /*
            AliTRDcluster *c = 0x0;
            for(Int_t itime = 0; itime < AliTRDseedV1::kNtb; itime++){
                c = tracklet->GetClusters(itime);
                if(!c) continue;
//                Float_t cluster_x = c->GetX();
//                Bool_t isprop;
//                isprop=kFALSE;
//                isprop=AliTracker::PropagateTrackTo(esd,(Double_t)cluster_x,0.105,5,kFALSE);
//                if(isprop)
//                {
//                    Int_t detector = c->GetDetector();
//                    ((TH1F*)fContainer->At(fSector))->Fill((c->GetY())-(esd->GetY()));
//                    if(c->GetY()-esd->GetY()>);
//                    printf("diff: %f\n",c->GetY()-(esd->GetY()));
//                }
            }
            */

        } // loop over tracklets ends

        if(counter_tracklet==AliTRDgeometry::kNlayer)
        {
            AliTRDtrackerV1::FitRiemanTilt(const_cast<AliTRDtrackV1 *>(TRDtrack), 0x0, kTRUE, counter_tracklet, points);
        }
        else continue;

       

        if(fDebugLevel>=1){
          for(Int_t il = 0; il < AliTRDgeometry::kNlayer; il++){
            //printf("---------------------------------------- %i %i %f %f %f \n",counter_tracklet,il,points[il].GetX(),points[il].GetY(),points[il].GetZ());
            Double_t point_x=points[il].GetX();
            Double_t point_y=points[il].GetY();
            Double_t point_z=points[il].GetZ();
            
            (*fDebugStream)   << "TrackletsinTRD"
                << "standalone=" << stand_alone
                << "eventcounter=" << fEventCounter
                << "layer="     << il
                << "dettracklet=" << det_tracklet
                << "xtracklet=" << tracklet_x[il]
                << "xtrack="    << point_x
                << "ytracklet=" << tracklet_y[il]
                << "ytrack="    << point_y
                << "ztracklet=" << tracklet_z[il]
                << "ztrack="    << point_z
                << "num_tracklets=" << counter_tracklet
                << "\n";
          }
        }

    }
  }

  
  PostData(0, fContainer);
}

//____________________________________________________________________
void AliTRDmultiplicity::Terminate(Option_t *)
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



