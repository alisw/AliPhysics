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

#include <TObjArray.h>
#include <TH1F.h>
#include <TTreeStream.h>

#include "AliESDtrack.h"
#include "AliExternalTrackParam.h"
#include "AliTRDseedV1.h"
#include "AliTRDtrackV1.h"
#include "AliTRDtrackerV1.h"
#include "AliTrackPointArray.h"

#include "AliTRDmultiplicity.h"
#include "AliTRDgeometry.h"
#include "info/AliTRDtrackInfo.h"

ClassImp(AliTRDmultiplicity)

//____________________________________________________________________
AliTRDmultiplicity::AliTRDmultiplicity()
    :AliTRDrecoTask()
  ,fEventCounter(0)
{
  //
  // Default constructor
    //
}

AliTRDmultiplicity::AliTRDmultiplicity(char* name)
    :AliTRDrecoTask(name, "Barrel Tracking Multiplicity")
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
void  AliTRDmultiplicity::UserCreateOutputObjects()
{
  //
  // Create output objects
  //

  TH1 *h = 0x0;
  fContainer = new TObjArray();
  for(Int_t is=0; is<AliTRDgeometry::kNsector; is++){
  fContainer->Add(h = new TH1F(Form("h_sector_%i", is), " ", 100,-10,10));
  }
} 

//____________________________________________________________________
void AliTRDmultiplicity::UserExec(Option_t *)
{
  //
  // Do it
  //



  ULong_t status;
  AliTRDtrackInfo     *track = 0x0;
  AliExternalTrackParam *esd = 0x0;
  AliTRDtrackV1 *trackTRD = 0x0;
  Double_t x_anode[AliTRDgeometry::kNlayer] = {300.2, 312.8, 325.4, 338.0, 350.6, 363.2}; // Take the default X0
  Int_t standAlone=1;
  fEventCounter++;


  for(Int_t itrk=0; itrk<fTracks->GetEntriesFast(); itrk++){
    track = (AliTRDtrackInfo*)fTracks->UncheckedAt(itrk);
    if(!track->HasESDtrack()) continue;
    status = track->GetStatus();

    AliTrackPoint points[6];
    Float_t xyz[3];
    memset(xyz, 0, sizeof(Float_t) * 3);
    Float_t trackletX[6];
    Float_t trackletY[6];
    Float_t trackletZ[6];
    for(Int_t a=0;a<6;a++)
    {
        xyz[0] = x_anode[a];
        points[a].SetXYZ(xyz);
        trackletX[a]=-1000;
        trackletY[a]=-1000;
        trackletZ[a]=-1000;
    }
    Int_t detTracklet=600;



    // TRD standalone
    if(((status&AliESDtrack::kTRDout)>0) && !((status&AliESDtrack::kTRDin)>0))
    {
        trackTRD = track->GetTrack();

        if(trackTRD)
        {
            AliTRDseedV1 *tracklet = 0x0;
            Int_t counterTracklet=0;
            for(Int_t itl = 0; itl < AliTRDgeometry::kNlayer; itl++)
            {
                tracklet = trackTRD->GetTracklet(itl);
                if(tracklet)
                {
                    if(tracklet->IsOK())
                    {
                        counterTracklet++;
                        detTracklet=tracklet->GetDetector();
                        //printf("%i %f %f \n",itl,tracklet->GetYfit(0),tracklet->GetZfit(0));
                        trackletX[itl]=tracklet->GetX0();
                        trackletY[itl]=tracklet->GetYfit(0);
                        trackletZ[itl]=tracklet->GetZfit(0);
                    }
                }
            }
            // this cut is needed otherwise the order of tracklets in the fit function is not correct
            if(counterTracklet==AliTRDgeometry::kNlayer) AliTRDtrackerV1::FitRiemanTilt(const_cast<AliTRDtrackV1 *>(trackTRD), 0x0, kTRUE, counterTracklet, points);
            else continue;

            if(DebugLevel()>=1){
              for(Int_t il = 0; il < AliTRDgeometry::kNlayer; il++){
                //printf("---------------------------------------- %i %i %f %f %f \n",counterTracklet,il,points[il].GetX(),points[il].GetY(),points[il].GetZ());
                Double_t pointX=points[il].GetX();
                Double_t pointY=points[il].GetY();
                Double_t pointZ=points[il].GetZ();


                (*DebugStream())   << "TrackletsinTRD"
                    << "standalone=" << standAlone
                    << "eventcounter=" << fEventCounter
                    << "layer="     << il
                    << "dettracklet=" << detTracklet
                    << "xtracklet=" << trackletX[il]
                    << "xtrack="    << pointX
                    << "ytracklet=" << trackletY[il]
                    << "ytrack="    << pointY
                    << "ztracklet=" << trackletZ[il]
                    << "ztrack="    << pointZ
                    << "num_tracklets=" << counterTracklet
                    << "\n";

              }
            }
        }
    } // end TRD standalone



    // TPC cluster selection for cosmic data
    if((track->GetTPCncls())<40) continue;
    // missing TPC propagation
    if(!(status&AliESDtrack::kTPCout)) continue;
    standAlone=0;

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
        trackTRD = track->GetTrack();
        if(!trackTRD) continue;
        AliTRDseedV1 *tracklet = 0x0;
        Int_t counterTracklet=0;


        for(Int_t itl = 0; itl < AliTRDgeometry::kNlayer; itl++){
            tracklet = trackTRD->GetTracklet(itl);
            if(!tracklet || !tracklet->IsOK()) continue;
            counterTracklet++;
            detTracklet=tracklet->GetDetector();
            trackletX[itl]=tracklet->GetX0();
            trackletY[itl]=tracklet->GetYfit(0);
            trackletZ[itl]=tracklet->GetZfit(0);

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

        if(counterTracklet==AliTRDgeometry::kNlayer)
        {
            AliTRDtrackerV1::FitRiemanTilt(const_cast<AliTRDtrackV1 *>(trackTRD), 0x0, kTRUE, counterTracklet, points);
        }
        else continue;

       

        if(DebugLevel()>=1){
          for(Int_t il = 0; il < AliTRDgeometry::kNlayer; il++){
            //printf("---------------------------------------- %i %i %f %f %f \n",counterTracklet,il,points[il].GetX(),points[il].GetY(),points[il].GetZ());
            Double_t pointX=points[il].GetX();
            Double_t pointY=points[il].GetY();
            Double_t pointZ=points[il].GetZ();
            
            (*DebugStream())   << "TrackletsinTRD"
                << "standalone=" << standAlone
                << "eventcounter=" << fEventCounter
                << "layer="     << il
                << "dettracklet=" << detTracklet
                << "xtracklet=" << trackletX[il]
                << "xtrack="    << pointX
                << "ytracklet=" << trackletY[il]
                << "ytrack="    << pointY
                << "ztracklet=" << trackletZ[il]
                << "ztrack="    << pointZ
                << "num_tracklets=" << counterTracklet
                << "\n";
          }
        }

    }
  }

  
  PostData(1, fContainer);
}


