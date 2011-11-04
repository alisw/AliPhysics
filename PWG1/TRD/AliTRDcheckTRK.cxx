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

#include "AliLog.h"

#include "AliTRDcluster.h"
#include "AliTRDseedV1.h"
#include "AliTRDtrackletOflHelper.h"
#include "AliTRDtrackV1.h"
#include "AliTRDtrackerV1.h"
#include "AliTRDtransform.h"

#include "AliTRDcheckTRK.h"

ClassImp(AliTRDcheckTRK)

Bool_t  AliTRDcheckTRK::fgKalmanUpdate = kTRUE;
Bool_t  AliTRDcheckTRK::fgClRecalibrate = kFALSE;
Float_t AliTRDcheckTRK::fgKalmanStep = 2.;
//__________________________________________________________________________
AliTRDcheckTRK::AliTRDcheckTRK()
  : AliTRDresolution()
{
// Default constructor
  SetNameTitle("TRDtrackerSys", "TRD Tracker Systematic");
}

//__________________________________________________________________________
AliTRDcheckTRK::AliTRDcheckTRK(char* name)
  : AliTRDresolution(name, kFALSE)
{
// User constructor
  SetTitle("TRD Tracker Systematic");
  MakePtCalib();
  InitFunctorList();
}

//__________________________________________________________________________
AliTRDcheckTRK::~AliTRDcheckTRK()
{
// Destructor
}


//__________________________________________________________________________
TObjArray* AliTRDcheckTRK::Histos()
{
// Build extra calibration plots
  if(!(fContainer = AliTRDresolution::Histos())) return NULL;
  //fContainer->Expand(AliTRDresolution::kNclasses+1);

  THnSparse *H(NULL);
  if(!(H = (THnSparseI*)gROOT->FindObject("Roads"))){
    const Char_t *title[kNdim] = {"layer", "charge", fgkTitle[kPt], fgkTitle[kYrez], fgkTitle[kPrez], "#sigma^{*}/<#sigma_{y}> [a.u.]", "n_{cl}"};
    const Int_t nbins[kNdim]   = {AliTRDgeometry::kNlayer, 2, kNptBins, fgkNbins[kYrez], fgkNbins[kPrez], kNSigmaBins, kNclusters};
    const Double_t min[kNdim]  = {-0.5, -0.5, -0.5, -1., -5., 0., 8.5},
                   max[kNdim]  = {AliTRDgeometry::kNlayer-0.5, 1.5, kNptBins-0.5, 1., 5., 5., min[6]+kNclusters};
    TString st("Tracking Roads Calib;");
    // define minimum info to be saved in non debug mode
    for(Int_t idim(0); idim<kNdim; idim++){ st += title[idim]; st+=";";}
    H = new THnSparseI("Roads", st.Data(), kNdim, nbins, min, max);
  } else H->Reset();
  fContainer->AddAt(H, fContainer->GetEntries()/*AliTRDresolution::kNclasses*/);
  return fContainer;
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
  if(!PropagateKalman(lt, fkESD->GetTPCoutParam())) return NULL;
  PlotCluster(&lt);
  PlotTracklet(&lt);
  PlotTrackIn(&lt);
  DoRoads(&lt);
  return NULL;
}

//________________________________________________________
void AliTRDcheckTRK::MakePtCalib(Float_t pt0, Float_t dpt)
{
// Build pt segments
  for(Int_t j(0); j<=kNptBins; j++){
    pt0+=(TMath::Exp(j*j*dpt)-1.);
    fPtBinCalib[j]=pt0;
  }
}

//__________________________________________________________________________
Int_t AliTRDcheckTRK::GetPtBinCalib(Float_t pt)
{
// Find pt bin according to local pt segmentation
  Int_t ipt(-1);
  while(ipt<kNptBins){
    if(pt<fPtBinCalib[ipt+1]) break;
    ipt++;
  }
  return ipt;
}


//__________________________________________________________________________
TH1* AliTRDcheckTRK::DoRoads(const AliTRDtrackV1 *track)
{
// comment needed
  if(track) fkTrack = track;
  if(!fkTrack){
    AliDebug(4, "No Track defined.");
    return NULL;
  }
  if(TMath::Abs(fkESD->GetTOFbc())>1){
    AliDebug(4, Form("Track with BC_index[%d] not used.", fkESD->GetTOFbc()));
    return NULL;
  }
  THnSparse *H(NULL);
  if(!fContainer || !(H = (THnSparse*)fContainer->At(3))){
    AliWarning("No output container defined.");
    return NULL;
  }
//  return NULL;
  Double_t val[kNdim];
  AliTRDseedV1 *fTracklet(NULL);
  AliTRDtrackletOflHelper helper; TObjArray cl(AliTRDseedV1::kNclusters); cl.SetOwner(kFALSE);
  for(Int_t il(0); il<AliTRDgeometry::kNlayer; il++){
    if(!(fTracklet = fkTrack->GetTracklet(il))) continue;
    if(!fTracklet->IsOK() || !fTracklet->IsChmbGood()) continue;
    Int_t ipt(GetPtBinCalib(fTracklet->GetPt()));
    if(ipt<0) continue;
    val[0] = il;
    val[1] = fkTrack->Charge()<0?0:1;
    val[2] = ipt;
    Double_t dyt(fTracklet->GetYfit(0) - fTracklet->GetYref(0)),
             dzt(fTracklet->GetZfit(0) - fTracklet->GetZref(0)),
             dydx(fTracklet->GetYfit(1)),
             tilt(fTracklet->GetTilt());
    // correct for tilt rotation
    val[3] = dyt - dzt*tilt;
    dydx+= tilt*fTracklet->GetZref(1);
    val[4] = TMath::ATan((dydx - fTracklet->GetYref(1))/(1.+ fTracklet->GetYref(1)*dydx)) * TMath::RadToDeg();
    fTracklet->ResetClusterIter(kTRUE); AliTRDcluster *c(NULL);
    while((c = fTracklet->NextCluster())) cl.AddLast(c);
    helper.Init(AliTRDtransform::Geometry().GetPadPlane(fTracklet->GetDetector()), &cl);
    Double_t r, y, s; helper.GetRMS(r, y, s, fTracklet->GetX0());
    val[5] = s/helper.GetSyMean();
    val[6] = fTracklet->GetN();
    H->Fill(val);
  }
  return NULL;
}

//___________________________________________________
Bool_t AliTRDcheckTRK::PropagateKalman(AliTRDtrackV1 &t, AliExternalTrackParam *ref)
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
  AliExternalTrackParam *trdin(NULL);
  if(!(trdin = t.GetTrackIn())){
    printf("E - AliTRDcheckTRK::PropagateKalman :: Track did not entered TRD fiducial volume.\n");
    return kFALSE;
  }
  if(!ref){
    printf("E - AliTRDcheckTRK::PropagateKalman :: Missing TPC out param.\n");
    return kFALSE;
  }
  if(ref->Pt()<1.e-3) return kFALSE;


  // Initialize TRD track to the reference
  AliTRDtrackV1 tt;
  tt.Set(ref->GetX(), ref->GetAlpha(), ref->GetParameter(), ref->GetCovariance());
  tt.SetMass(t.GetMass());
  tt.SetTrackOut(t.GetTrackOut());

  for(Int_t ily(0); ily<AliTRDgeometry::kNlayer; ily++){
    if(!(tr = t.GetTracklet(ily))) continue;
    Int_t det(tr->GetDetector());
    //Float_t *calib = GetCalib(det);
    if(fgClRecalibrate/* && calib[0]>0.*/){
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
    if(!tt.GetTrackIn()) tt.SetTrackIn();
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


