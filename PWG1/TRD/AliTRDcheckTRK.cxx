/**************************************************************************
* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
*                                                                        *
* Author: The ALICE Off-line Project.                                    *
* Contributors are mentioned in the code where appropriate.              *
*                                                                        *
* Permission to use, copy, modify and distribute this software and its   *
* documentation strictly for non-commercialf purposes is hereby granted   *
* without fee, provided that the above copyright notice appears in all   *
* copies and that both the copyright notice and this permission notice   *
* appear in the supporting documentation. The authors make no claims     *
* about the suitability of this software for any purpose. It is          *
* provided "as is" without express or implied warranty.                  *
**************************************************************************/


////////////////////////////////////////////////////////////////////////////
//                                                                        //
//  TRD tracker systematic                                                //
//
//
//  Authors:                                                              //
//    Alexandru Bercuci <A.Bercuci@gsi.de>                                //
//                                                                        //
////////////////////////////////////////////////////////////////////////////

#include "TROOT.h"
#include "TAxis.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TObjArray.h"
#include "THnSparse.h"
#include <TVectorT.h>

#include "AliLog.h"

#include "AliTRDcluster.h"
#include "AliTRDseedV1.h"
#include "AliTRDtrackV1.h"
#include "AliTRDtrackerV1.h"
#include "AliTRDtransform.h"

#include "AliTRDcheckTRK.h"

ClassImp(AliTRDcheckTRK)

Bool_t  AliTRDcheckTRK::fgKalmanUpdate = kTRUE;
Bool_t  AliTRDcheckTRK::fgClRecalibrate = kFALSE;
Float_t AliTRDcheckTRK::fgKalmanStep = 2.;
Float_t AliTRDcheckTRK::fgCalib[540][2];
//__________________________________________________________________________
AliTRDcheckTRK::AliTRDcheckTRK()
  : AliTRDresolution()
{
// Default constructor
  SetNameTitle("TRDtrackerSys", "TRD Tracker Systematic");
  memset(AliTRDcheckTRK::fgCalib, 0, 540*2*sizeof(Float_t));
}

//__________________________________________________________________________
AliTRDcheckTRK::AliTRDcheckTRK(char* name)
  : AliTRDresolution(name, kFALSE)
{
// User constructor
  SetTitle("TRD Tracker Systematic");
  memset(AliTRDcheckTRK::fgCalib, 0, 540*2*sizeof(Float_t));
  InitFunctorList();
}

//__________________________________________________________________________
AliTRDcheckTRK::~AliTRDcheckTRK()
{
// Destructor
}

//__________________________________________________________________________
Int_t AliTRDcheckTRK::GetSpeciesByMass(Float_t m)
{
// Find particle index by mass
// 0 electron
// 1 muon
// 2 pion
// 3 kaon
// 4 proton

  for(Int_t is(0); is<AliPID::kSPECIES; is++) if(TMath::Abs(m-AliPID::ParticleMass(is))<1.e-4) return is;
  return -1;
}



//__________________________________________________________________________
TH1* AliTRDcheckTRK::PlotTrack(const AliTRDtrackV1 *track)
{
// comment needed

  if(track) fkTrack = track;
  if(!fkTrack){
    AliDebug(4, "No Track defined.");
    return NULL;
  }
  // make a local copy of current track
  AliTRDtrackV1 lt(*fkTrack);
  if(!PropagateKalman(lt)) return NULL;
  PlotCluster(&lt);
  PlotTracklet(&lt);
  PlotTrackIn(&lt);
  return NULL;
}


//___________________________________________________
Bool_t AliTRDcheckTRK::PropagateKalman(AliTRDtrackV1 &t)
{
// Propagate Back Kalman from the TPC input parameter to the last tracklet attached to track.
// On the propagation recalibration of clusters, tracklet refit and material budget are recalculated (on demand)
// On output the track is updated with the new info
//
// A.Bercuci@gsi.de

//  printf("PropagateKalman()\n");
  Int_t ntracklets(t.GetNumberOfTracklets());
  if(!ntracklets){
    printf("E - AliTRDcheckTRK::PropagateKalman :: No tracklets attached to track.\n");
    return kFALSE;
  }

  AliTRDseedV1 *tr(NULL);
  AliExternalTrackParam *ref(NULL);
  if(!(ref = t.GetTrackIn())){
    printf("E - AliTRDcheckTRK::PropagateKalman :: Track did not entered TRD fiducial volume.\n");
    return kFALSE;
  }
  if(ref->Pt()<1.e-3){printf("small refpt\n"); return kFALSE;}


  // Initialize TRD track to the reference
  AliTRDtrackV1 tt;
  tt.Set(ref->GetX(), ref->GetAlpha(), ref->GetParameter(), ref->GetCovariance());
  tt.SetMass(t.GetMass());
  tt.SetTrackIn();tt.SetTrackOut(t.GetTrackOut());

  for(Int_t ily(0); ily<AliTRDgeometry::kNlayer; ily++){
    if(!(tr = t.GetTracklet(ily))) continue;
    Int_t det(tr->GetDetector());
    Float_t *calib = GetCalib(det);
    if(fgClRecalibrate && calib[0]>0.){
      AliTRDtransform trans(det);
      AliTRDcluster *c(NULL);
      Float_t exb, vd, t0, s2, dl, dt; tr->GetCalibParam(exb, vd, t0, s2, dl, dt);
      tr->ResetClusterIter(kFALSE);
      while((c = tr->PrevCluster())){
        if(!trans.Transform(c/*, GetCalib(det)*/)){
          printf("W - AliTRDcheckTRK::PropagateKalman :: Transform() failed for Det[%03d]\n", det);
          break;
        }
      }
      if(!tr->FitRobust()) printf("W - AliTRDcheckTRK::PropagateKalman :: FitRobust() failed for Det[%03d]\n", det);
    }
    if(!AliTRDtrackerV1::PropagateToX(tt, tr->GetX0(), fgKalmanStep)) continue;
    tr->Update(&tt);
    if(fgKalmanUpdate){
      Double_t x(tr->GetX0()),
               p[2] = { tr->GetYfit(0), tr->GetZfit(0)},
               covTrklt[3];
      tr->GetCovAt(x, covTrklt);
      if(!((AliExternalTrackParam&)tt).Update(p, covTrklt)) continue;
      //tr->Update(&tt);
      tt.SetTracklet(tr, 0);
      tt.SetNumberOfClusters();
      tt.UpdateChi2(((AliExternalTrackParam)tt).GetPredictedChi2(p, covTrklt));
    }
  }
  //tt.Print("a");
  t.~AliTRDtrackV1();
  new(&t) AliTRDtrackV1(tt);
  return kTRUE;
}


