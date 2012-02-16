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

/* $Id: AliTRDresolution.cxx 27496 2008-07-22 08:35:45Z cblume $ */

////////////////////////////////////////////////////////////////////////////
//                                                                        //
//  TRD tracking resolution                                               //
//
// The class performs resolution and residual studies 
// of the TRD tracks for the following quantities :
//   - spatial position (y, [z])
//   - angular (phi) tracklet
//   - momentum at the track level
// 
// The class has to be used for regular detector performance checks using the official macros:
//   - $ALICE_ROOT/TRD/qaRec/run.C
//   - $ALICE_ROOT/TRD/qaRec/makeResults.C
// 
// For stand alone usage please refer to the following example: 
// {  
//   gSystem->Load("libANALYSIS.so");
//   gSystem->Load("libTRDqaRec.so");
//   AliTRDresolution *res = new AliTRDresolution();
//   //res->SetMCdata();
//   //res->SetVerbose();
//   //res->SetVisual();
//   res->Load();
//   if(!res->PostProcess()) return;
//   res->GetRefFigure(0);
// }  
//
//  Authors:                                                              //
//    Alexandru Bercuci <A.Bercuci@gsi.de>                                //
//    Markus Fasel <M.Fasel@gsi.de>                                       //
//                                                                        //
////////////////////////////////////////////////////////////////////////////

#include <TSystem.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TObjArray.h>
#include <TH3.h>
#include <TH2.h>
#include <TH1.h>
#include <THnSparse.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TGaxis.h>
#include <TBox.h>
#include <TLegend.h>
#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>
#include <TLinearFitter.h>
#include <TMath.h>
#include <TMatrixT.h>
#include <TVectorT.h>
#include <TTreeStream.h>
#include <TGeoManager.h>
#include <TDatabasePDG.h>

#include "AliPID.h"
#include "AliLog.h"
#include "AliESDtrack.h"
#include "AliMathBase.h"
#include "AliTrackPointArray.h"

#include "AliTRDresolution.h"
#include "AliTRDgeometry.h"
#include "AliTRDtransform.h"
#include "AliTRDpadPlane.h"
#include "AliTRDcluster.h"
#include "AliTRDseedV1.h"
#include "AliTRDtrackV1.h"
#include "AliTRDReconstructor.h"
#include "AliTRDrecoParam.h"
#include "AliTRDpidUtil.h"
#include "AliTRDinfoGen.h"
#include "AliAnalysisManager.h"
#include "info/AliTRDclusterInfo.h"
#include "info/AliTRDeventInfo.h"

ClassImp(AliTRDresolution)
ClassImp(AliTRDresolution::AliTRDresolutionProjection)

Int_t const   AliTRDresolution::fgkNbins[kNdim] = {
  Int_t(kNbunchCross)/*bc*/,
  180/*phi*/,
  50/*eta*/,
  50/*dy*/,
  40/*dphi*/,
  50/*dz*/,
  3/*chg*species*/,
  kNpt/*pt*/
};  //! no of bins/projection
Double_t const AliTRDresolution::fgkMin[kNdim] = {
  -1.5,
  -TMath::Pi(),
  -1.,
  -1.5,
  -10.,
  -2.5,
  -1.5,
  -0.5
};    //! low limits for projections
Double_t const AliTRDresolution::fgkMax[kNdim] = {
  1.5,
  TMath::Pi(),
  1.,
  1.5,
  10.,
  2.5,
  1.5,
  kNpt-0.5
};    //! high limits for projections
Char_t const *AliTRDresolution::fgkTitle[kNdim] = {
  "bunch cross",
  "#phi [rad]",
  "#eta",
  "#Deltay [cm]",
  "#Delta#phi [deg]",
  "#Deltaz [cm]",
  "chg*spec*rc",
  "bin_p_{t}"
};  //! title of projection

Char_t const * AliTRDresolution::fgPerformanceName[kNclasses] = {
    "Det2Cluster"
    ,"Cluster2Track"
    ,"Tracklet2Track"
    ,"Tracklet2TRDin"
    ,"Cluster2MC"
    ,"Tracklet2MC"
    ,"TRDin2MC"
    ,"TRD2MC"
//    ,"Tracklet2TRDout"
//    ,"TRDout2MC"
};
Float_t AliTRDresolution::fgPtBin[25];

//________________________________________________________
AliTRDresolution::AliTRDresolution()
  :AliTRDrecoTask()
  ,fSteer(0)
  ,fIdxFrame(0)
  ,fPtThreshold(.3)
  ,fDyRange(0.75)
  ,fBCbinTOF(0)
  ,fBCbinFill(0)
  ,fBsign(kFALSE)
  ,fProj(NULL)
  ,fDBPDG(NULL)
  ,fCl(NULL)
  ,fMCcl(NULL)
{
  //
  // Default constructor
  //
  SetNameTitle("TRDresolution", "TRD spatial and momentum resolution");
  MakePtSegmentation();
  SetProcesses(kTRUE, kTRUE, kTRUE, kTRUE);
}

//________________________________________________________
AliTRDresolution::AliTRDresolution(char* name, Bool_t xchange)
  :AliTRDrecoTask(name, "TRD spatial and momentum resolution")
  ,fSteer(0)
  ,fIdxFrame(0)
  ,fPtThreshold(.3)
  ,fDyRange(0.75)
  ,fBCbinTOF(0)
  ,fBCbinFill(0)
  ,fBsign(kFALSE)
  ,fProj(NULL)
  ,fDBPDG(NULL)
  ,fCl(NULL)
  ,fMCcl(NULL)
{
  //
  // Default constructor
  //

  InitFunctorList();
  MakePtSegmentation();
  SetProcesses(kTRUE, kTRUE, kTRUE, kTRUE);
  if(xchange){
    SetUseExchangeContainers();
    DefineOutput(kClToTrk, TObjArray::Class()); // cluster2track
    DefineOutput(kClToMC, TObjArray::Class()); // cluster2mc
  }
}

//________________________________________________________
AliTRDresolution::~AliTRDresolution()
{
  //
  // Destructor
  //
  if (AliAnalysisManager::GetAnalysisManager() && AliAnalysisManager::GetAnalysisManager()->IsProofMode()) return;
  if(fProj){fProj->Delete(); delete fProj;}
  if(fCl){fCl->Delete(); delete fCl;}
  if(fMCcl){fMCcl->Delete(); delete fMCcl;}
}


//________________________________________________________
void AliTRDresolution::UserCreateOutputObjects()
{
  // spatial resolution

  AliTRDrecoTask::UserCreateOutputObjects();
  if(UseExchangeContainers()) InitExchangeContainers();
}

//________________________________________________________
void AliTRDresolution::InitExchangeContainers()
{
// Init containers for subsequent tasks (AliTRDclusterResolution)

  fCl = new TObjArray(200); fCl->SetOwner(kTRUE);
  fMCcl = new TObjArray(); fMCcl->SetOwner(kTRUE);
  PostData(kClToTrk, fCl);
  PostData(kClToMC, fMCcl);
}

//________________________________________________________
void AliTRDresolution::UserExec(Option_t *opt)
{
  //
  // Execution part
  //

  if(fCl) fCl->Delete();
  if(fMCcl) fMCcl->Delete();
  AliTRDrecoTask::UserExec(opt);
}

//________________________________________________________
Bool_t AliTRDresolution::Pulls(Double_t* /*dyz[2]*/, Double_t* /*cov[3]*/, Double_t /*tilt*/) const
{
// Helper function to calculate pulls in the yz plane
// using proper tilt rotation
// Uses functionality defined by AliTRDseedV1.

  return kTRUE;
/*
  Double_t t2(tilt*tilt);
  // exit door until a bug fix is found for AliTRDseedV1::GetCovSqrt

  // rotate along pad
  Double_t cc[3];
  cc[0] = cov[0] - 2.*tilt*cov[1] + t2*cov[2];
  cc[1] = cov[1]*(1.-t2) + tilt*(cov[0] - cov[2]);
  cc[2] = t2*cov[0] + 2.*tilt*cov[1] + cov[2];
  // do sqrt
  Double_t sqr[3]={0., 0., 0.};
  if(AliTRDseedV1::GetCovSqrt(cc, sqr)) return kFALSE;
  Double_t invsqr[3]={0., 0., 0.};
  if(AliTRDseedV1::GetCovInv(sqr, invsqr)<1.e-40) return kFALSE;
  Double_t tmp(dyz[0]);
  dyz[0] = invsqr[0]*tmp + invsqr[1]*dyz[1];
  dyz[1] = invsqr[1]*tmp + invsqr[2]*dyz[1];
  return kTRUE;
*/
}

//________________________________________________________
TH1* AliTRDresolution::DetCluster(const TObjArray *cls)
{
  //
  // Plot the cluster distributions
  //

  if(cls) fkClusters = cls;
  if(!fkClusters){
    AliDebug(4, "No Clusters defined.");
    return NULL;
  }
  Int_t ncl(0);
  if(!(ncl = fkClusters->GetEntriesFast())){
    AliDebug(1, "No RecPoints defined.");
    return NULL;
  }
  Int_t det(-1);
  AliTRDcluster *cl(NULL);
  for(Int_t icl(0); icl<ncl; icl++){
    if(!(cl = (AliTRDcluster*)(*fkClusters)[icl])) continue;
    det = cl->GetDetector(); break;
  }
  if(det<0){
    AliDebug(1, "No useful clusters defined.");
    return NULL;
  }
  THnSparse *H(NULL);
  if(!fContainer || !(H = ((THnSparse*)fContainer->At(kDetector)))){
    AliWarning("No output container defined.");
    return NULL;
  }
  Int_t ly(AliTRDgeometry::GetLayer(det)),
        stk(AliTRDgeometry::GetStack(det));
  Double_t val[kNdim],
           alpha((0.5+AliTRDgeometry::GetSector(det))*AliTRDgeometry::GetAlpha()),
           cs(TMath::Cos(alpha)),
           sn(TMath::Sin(alpha));
  for(Int_t icl(0); icl<ncl; icl++){
    if(!(cl = (AliTRDcluster*)(*fkClusters)[icl])) continue;
    val[kBC]  = ly;
    val[kPhi] = TMath::ATan2(cl->GetX()*sn + cl->GetY()*cs, cl->GetX()*cs - cl->GetY()*sn);
    val[kEta] = (5-stk)*16-cl->GetPadRow()-1-(stk<3?4:0);
    val[kYrez]= fEvent->GetMultiplicity();
    val[4]    = TMath::Abs(cl->GetQ());
    val[5]    = cl->IsFivePad()?1:cl->GetNPads();
    H->Fill(val);
  }
  return NULL;
}

//________________________________________________________
TH1* AliTRDresolution::PlotCluster(const AliTRDtrackV1 *track)
{
  //
  // Plot the cluster distributions
  //

  if(track) fkTrack = track;
  if(!fkTrack){
    AliDebug(4, "No Track defined.");
    return NULL;
  }
  if(fkESD && TMath::Abs(fkESD->GetTOFbc()) > 1){
    AliDebug(4, Form("Track with BC_index[%d] not used.", fkESD->GetTOFbc()));
    return NULL;
  }
  if(fkESD && fPt<fPtThreshold){
    AliDebug(4, Form("Track with pt[%6.4f] under threshold.", fPt));
    return NULL;
  }
  THnSparse *H(NULL);
  if(!fContainer || !(H = ((THnSparse*)fContainer->At(kCluster)))){
    AliWarning("No output container defined.");
    return NULL;
  }

  AliTRDgeometry *geo(AliTRDinfoGen::Geometry());
  Double_t val[kNdim+2],
           alpha(0.), cs(-2.), sn(0.); //Float_t exb, vd, t0, s2, dl, dt;
  TObjArray     *clInfoArr(NULL);
  AliTRDseedV1  *fTracklet(NULL);
  AliTRDcluster *c(NULL), *cc(NULL);
  for(Int_t ily=0; ily<AliTRDgeometry::kNlayer; ily++){
    if(!(fTracklet = fkTrack->GetTracklet(ily))) continue;
    if(!fTracklet->IsOK()) continue;
    //fTracklet->GetCalibParam(exb, vd, t0, s2, dl, dt);
    val[kBC]  = ily;
    if(cs<-1.){
      alpha = (0.5+AliTRDgeometry::GetSector(fTracklet->GetDetector()))*AliTRDgeometry::GetAlpha();
      cs    = TMath::Cos(alpha);
      sn    = TMath::Sin(alpha);
    }
    val[kPhi] = TMath::ATan2(fTracklet->GetX()*sn + fTracklet->GetY()*cs, fTracklet->GetX()*cs - fTracklet->GetY()*sn);
    Float_t tgl = fTracklet->GetZ()/fTracklet->GetX()/TMath::Sqrt(1.+fTracklet->GetY()*fTracklet->GetY()/fTracklet->GetX()/fTracklet->GetX());
    val[kEta] = -TMath::Log(TMath::Tan(0.5 *  (0.5*TMath::Pi() - TMath::ATan(tgl))));
    val[kPt]  = TMath::ATan(fTracklet->GetYref(1))*TMath::RadToDeg();
    Float_t corr = 1./TMath::Sqrt(1.+fTracklet->GetYref(1)*fTracklet->GetYref(1)+fTracklet->GetZref(1)*fTracklet->GetZref(1));
    Int_t row0(-1);
    Float_t padCorr(fTracklet->GetTilt()*fTracklet->GetPadLength());
    fTracklet->ResetClusterIter(kTRUE);
    while((c = fTracklet->NextCluster())){
      Float_t xc(c->GetX()),
              q(TMath::Abs(c->GetQ()));
      if(row0<0) row0 = c->GetPadRow();

      val[kYrez] = c->GetY() + padCorr*(c->GetPadRow() - row0) -fTracklet->GetYat(xc);
      val[kPrez] = fTracklet->GetX0()-xc;
      val[kZrez] = 0.; Int_t ic(0), tb(c->GetLocalTimeBin());;
      if((cc = fTracklet->GetClusters(tb-1))) {val[kZrez] += TMath::Abs(cc->GetQ()); ic++;}
      if((cc = fTracklet->GetClusters(tb-2))) {val[kZrez] += TMath::Abs(cc->GetQ()); ic++;}
      if(ic) val[kZrez] /= (ic*q);
      val[kSpeciesChgRC]= fTracklet->IsRowCross()?0.:(TMath::Max(q*corr, Float_t(3.)));
      val[kNdim]   = fEvent->GetMultiplicity();
      val[kNdim+1] = c->IsFivePad()?1:c->GetNPads();
      H->Fill(val);
/*      // tilt rotation of covariance for clusters
      Double_t sy2(c->GetSigmaY2()), sz2(c->GetSigmaZ2());
      cov[0] = (sy2+t2*sz2)*corr;
      cov[1] = tilt*(sz2 - sy2)*corr;
      cov[2] = (t2*sy2+sz2)*corr;
      // sum with track covariance
      cov[0]+=covR[0]; cov[1]+=covR[1]; cov[2]+=covR[2];
      Double_t dyz[2]= {dy[1], dz[1]};
      Pulls(dyz, cov, tilt);*/

      // Get z-position with respect to anode wire
      Float_t yt(fTracklet->GetYref(0)-val[kZrez]*fTracklet->GetYref(1)),
              zt(fTracklet->GetZref(0)-val[kZrez]*fTracklet->GetZref(1));
      Int_t istk = geo->GetStack(c->GetDetector());
      AliTRDpadPlane *pp = geo->GetPadPlane(ily, istk);
      Float_t rowZ = pp->GetRow0();
      Float_t d  = rowZ - zt + pp->GetAnodeWireOffset();
      d -= ((Int_t)(2 * d)) / 2.0;
      if (d > 0.25) d  = 0.5 - d;

      AliTRDclusterInfo *clInfo(NULL);
      clInfo = new AliTRDclusterInfo;
      clInfo->SetCluster(c);
      //Float_t covF[] = {cov[0], cov[1], cov[2]};
      clInfo->SetGlobalPosition(yt, zt, fTracklet->GetYref(1), fTracklet->GetZref(1)/*, covF*/);
      clInfo->SetResolution(val[kYrez]);
      clInfo->SetAnisochronity(d);
      clInfo->SetDriftLength(val[kZrez]);
      clInfo->SetTilt(fTracklet->GetTilt());
      if(fCl) fCl->Add(clInfo);
      //else AliDebug(1, "Cl exchange container missing. Activate by calling \"InitExchangeContainers()\"");

      if(DebugLevel()>=2){
        if(!clInfoArr){
          clInfoArr=new TObjArray(AliTRDseedV1::kNclusters);
          clInfoArr->SetOwner(kFALSE);
        }
        clInfoArr->Add(clInfo);
      }
    }
    if(DebugLevel()>=2 && clInfoArr){
      ULong_t status = fkESD->GetStatus();
      (*DebugStream()) << "cluster"
        <<"status="  << status
        <<"clInfo.=" << clInfoArr
        << "\n";
      clInfoArr->Clear();
    }
  }
  if(clInfoArr) delete clInfoArr;

  if(!track) return NULL;
  // special care for EVE usage
  TH1 *h(NULL);
  if((h = (TH1*)gDirectory->Get(Form("%s_proj_%d", H->GetName(), kYrez)))) delete h;
  return H->Projection(kYrez);
}


//________________________________________________________
TH1* AliTRDresolution::PlotTracklet(const AliTRDtrackV1 *track)
{
// Plot normalized residuals for tracklets to track.
//
// We start from the result that if X=N(|m|, |Cov|)
// BEGIN_LATEX
// (Cov^{-1})^{1/2}X = N((Cov^{-1})^{1/2}*|m|, |1|)
// END_LATEX
// in our case X=(y_trklt - y_trk z_trklt - z_trk) and |Cov| = |Cov_trklt| + |Cov_trk| at the radial
// reference position.
  if(track) fkTrack = track;
  if(!fkTrack){
    AliDebug(4, "No Track defined.");
    return NULL;
  }
  if(fkESD && TMath::Abs(fkESD->GetTOFbc())>1){
    AliDebug(4, Form("Track with BC_index[%d] not used.", fkESD->GetTOFbc()));
    //return NULL;
  }
  THnSparse *H(NULL);
  if(!fContainer || !(H = (THnSparse*)fContainer->At(kTracklet))){
    AliWarning("No output container defined.");
    return NULL;
  }

  const Int_t ndim(kNdim+8);
  Double_t val[ndim],
           alpha(0.), cs(-2.), sn(0.);
  Float_t sz[AliTRDseedV1::kNtb], pos[AliTRDseedV1::kNtb];
  Int_t ntbGap[AliTRDseedV1::kNtb];
  AliTRDseedV1 *fTracklet(NULL);
  for(Int_t il(0); il<AliTRDgeometry::kNlayer; il++){
    if(!(fTracklet = fkTrack->GetTracklet(il))) continue;
    if(!fTracklet->IsOK() || !fTracklet->IsChmbGood()) continue;
    val [kBC] = il;
    if(cs<-1.){
      alpha = (0.5+AliTRDgeometry::GetSector(fTracklet->GetDetector()))*AliTRDgeometry::GetAlpha();
      cs    = TMath::Cos(alpha);
      sn    = TMath::Sin(alpha);
    }
    val[kPhi] = TMath::ATan2(fTracklet->GetX()*sn + fTracklet->GetY()*cs, fTracklet->GetX()*cs - fTracklet->GetY()*sn);
    Float_t tgl = fTracklet->GetZ()/fTracklet->GetX()/TMath::Sqrt(1.+fTracklet->GetY()*fTracklet->GetY()/fTracklet->GetX()/fTracklet->GetX());//fTracklet->GetTgl();
    val[kEta] = -TMath::Log(TMath::Tan(0.5 *  (0.5*TMath::Pi() - TMath::ATan(tgl))));

    val[kSpeciesChgRC]= fTracklet->IsRowCross()?0:fkTrack->Charge();// fSpecies;
    val[kPt]  = fPt<0.8?0:(fPt<1.5?1:2);//GetPtBin(fTracklet->GetMomentum());
    Double_t dyt(fTracklet->GetYfit(0) - fTracklet->GetYref(0)),
             dzt(fTracklet->GetZfit(0) - fTracklet->GetZref(0)),
             dydx(fTracklet->GetYfit(1)),
             tilt(fTracklet->GetTilt());
    // correct for tilt rotation
    val[kYrez] = dyt - dzt*tilt;
    val[kZrez] = fTracklet->IsRowCross()?(dzt + dyt*tilt):(fTracklet->GetdQdl()*3.e-4-1.5);
    dydx+= tilt*fTracklet->GetZref(1);
    val[kPrez] = TMath::ATan((dydx - fTracklet->GetYref(1))/(1.+ fTracklet->GetYref(1)*dydx)) * TMath::RadToDeg();
    val[kNdim] = fEvent?fEvent->GetMultiplicity():0;
    val[kNdim+1] = 1.e2*fTracklet->GetTBoccupancy()/AliTRDseedV1::kNtb;
    Int_t n = fTracklet->GetChargeGaps(sz, pos, ntbGap);
    val[kNdim+2] = 0.; for(Int_t igap(0); igap<n; igap++) val[kNdim+2] += sz[igap];
    for(Int_t ifill(0); ifill<3; ifill++){ val[kNdim+3+ifill]=0.;val[kNdim+4+ifill]=0.;}
    for(Int_t igap(0), ifill(0); igap<n&&ifill<2; igap++){
      if(ntbGap[igap]<2) continue;
      val[kNdim+3+ifill] = sz[igap];
      val[kNdim+4+ifill] = pos[igap];
      ifill++;
    }
    H->Fill(val);

//     // compute covariance matrix
//     fTracklet->GetCovAt(x, cov);
//     fTracklet->GetCovRef(covR);
//     cov[0] += covR[0]; cov[1] += covR[1]; cov[2] += covR[2];
//     Double_t dyz[2]= {dy[1], dz[1]};
//     Pulls(dyz, cov, tilt);
//     ((TH3S*)arr->At(1))->Fill(sgm[fSegmentLevel], dyz[0], dyz[1]);
//     ((TH3S*)arr->At(3))->Fill(tht, dyz[1], rc);

    if(DebugLevel()>=3){
      Bool_t rc(fTracklet->IsRowCross());
      UChar_t err(fTracklet->GetErrorMsg());
      Double_t x(fTracklet->GetX()),
               pt(fTracklet->GetPt()),
               yt(fTracklet->GetYref(0)),
               zt(fTracklet->GetZref(0)),
               phi(fTracklet->GetYref(1)),
               tht(fTracklet->GetZref(1));
      Int_t ncl(fTracklet->GetN()), det(fTracklet->GetDetector());
      (*DebugStream()) << "tracklet"
        <<"pt="  << pt
        <<"x="   << x
        <<"yt="  << yt
        <<"zt="  << zt
        <<"phi=" << phi
        <<"tht=" << tht
        <<"det=" << det
        <<"n="   << ncl
        <<"dy0=" << dyt
        <<"dz0=" << dzt
        <<"dy="  << val[kYrez]
        <<"dz="  << val[kZrez]
        <<"dphi="<< val[kPrez]
        <<"dQ  ="<< val[kNdim]
        <<"rc="  << rc
        <<"err=" << err
        << "\n";
    }
  }
  if(!track) return NULL;
  // special care for EVE usage
  TH1 *h(NULL);
  if((h = (TH1*)gDirectory->Get(Form("%s_proj_%d", H->GetName(), kYrez)))) delete h;
  return H->Projection(kYrez);
}


//________________________________________________________
TH1* AliTRDresolution::PlotTrackIn(const AliTRDtrackV1 *track)
{
// Store resolution/pulls of Kalman before updating with the TRD information
// at the radial position of the first tracklet. The following points are used
// for comparison
//  - the (y,z,snp) of the first TRD tracklet
//  - the (y, z, snp, tgl, pt) of the MC track reference
//
// Additionally the momentum resolution/pulls are calculated for usage in the
// PID calculation.
  //printf("AliTRDresolution::PlotTrackIn() :: track[%p]\n", (void*)track);

  if(track) fkTrack = track;
  if(!fkTrack){
    AliDebug(4, "No Track defined.");
    return NULL;
  }
  //fkTrack->Print();
  TH1 *h(NULL); // EVE projection
  // check container
  THnSparseI *H=(THnSparseI*)fContainer->At(kTrackIn);
  if(!H){
    AliError(Form("Missing container @ %d", Int_t(kTrackIn)));
    return NULL;
  }
  // check input track status
  AliExternalTrackParam *tin(NULL);
  if(!(tin = fkTrack->GetTrackIn())){
    AliError("Track did not entered TRD fiducial volume.");
    return NULL;
  }
  // check first tracklet
  AliTRDseedV1 *fTracklet(fkTrack->GetTracklet(0));
  if(!fTracklet){
    AliDebug(3, "No Tracklet in ly[0]. Skip track.");
    if(!track) return NULL;
    // special care for EVE usage
    if((h = (TH1*)gDirectory->Get(Form("%s_proj_%d", H->GetName(), kYrez)))) delete h;
    return H->Projection(kYrez);
  }
  if(!fTracklet->IsOK() || !fTracklet->IsChmbGood()){
    AliDebug(3, "Tracklet or Chamber not OK. Skip track.");
    if(!track) return NULL;
    // special care for EVE usage
    if((h = (TH1*)gDirectory->Get(Form("%s_proj_%d", H->GetName(), kYrez)))) delete h;
    return H->Projection(kYrez);
  }
  // check radial position
  Double_t x = tin->GetX();
  if(TMath::Abs(x-fTracklet->GetX())>1.e-3){
    AliDebug(1, Form("Tracklet did not match Track. dx[cm]=%+4.1f", x-fTracklet->GetX()));
    if(!track) return NULL;
    // special care for EVE usage
    if((h = (TH1*)gDirectory->Get(Form("%s_proj_%d", H->GetName(), kYrez)))) delete h;
    return H->Projection(kYrez);
  }
//  printf("USE y[%+f] dydx[%+f]\n", fTracklet->GetYfit(0), fTracklet->GetYfit(1));

  Int_t bc(fkESD?fkESD->GetTOFbc()/2:0);
  const Double_t *parR(tin->GetParameter());
  Double_t dyt(fTracklet->GetYfit(0)-parR[0]), dzt(fTracklet->GetZfit(0)-parR[1]),
            phit(fTracklet->GetYfit(1)),
            tilt(fTracklet->GetTilt()),
            norm(1./TMath::Sqrt((1.-parR[2])*(1.+parR[2])));

  // correct for tilt rotation
  Double_t dy  = dyt - dzt*tilt,
           dz  = dzt + dyt*tilt,
           dx  = dy/(parR[2]*norm-parR[3]*norm*tilt);
  phit       += tilt*parR[3];
  Double_t dphi = TMath::ATan(phit) - TMath::ASin(parR[2]);

  Double_t val[kNdim+3];
  val[kBC]          = bc==0?0:(bc<0?-1.:1.);
  Double_t alpha = (0.5+AliTRDgeometry::GetSector(fTracklet->GetDetector()))*AliTRDgeometry::GetAlpha(),
           cs    = TMath::Cos(alpha),
           sn    = TMath::Sin(alpha);
  val[kPhi] = TMath::ATan2(fTracklet->GetX()*sn + fTracklet->GetY()*cs, fTracklet->GetX()*cs - fTracklet->GetY()*sn);
  Float_t tgl = fTracklet->GetZ()/fTracklet->GetX()/TMath::Sqrt(1.+fTracklet->GetY()*fTracklet->GetY()/fTracklet->GetX()/fTracklet->GetX());
  val[kEta] = -TMath::Log(TMath::Tan(0.5 *  (0.5*TMath::Pi() - TMath::ATan(tgl))));
  val[kYrez]        = dy;
  val[kZrez]        = fTracklet->IsRowCross()?dz:(fTracklet->GetdQdl()*5.e-4 - 2.5);
  val[kPrez]        = dphi*TMath::RadToDeg();
  val[kSpeciesChgRC]= fTracklet->IsRowCross()?0:fSpecies; //printf("val[%d] = %4.2f\n", kSpeciesChgRC, val[kSpeciesChgRC]);
  val[kPt]          = fPt<0.8?0:(fPt<1.5?1:2);//GetPtBin(fPt);
  val[kNdim]        = GetPtBin(fTracklet->GetMomentum());
  val[kNdim+1]      = dx;
  val[kNdim+2]      = fEvent?fEvent->GetBunchFill():0;
  H->Fill(val);
  if(DebugLevel()>=3){
    (*DebugStream()) << "trackIn"
      <<"tracklet.="  << fTracklet
      <<"trackIn.="   << tin
      << "\n";
  }

  if(!track) return NULL;
  // special care for EVE usage
    if((h = (TH1*)gDirectory->Get(Form("%s_proj_%d", H->GetName(), kYrez)))) delete h;
    return H->Projection(kYrez);
}

/*
//________________________________________________________
TH1* AliTRDresolution::PlotTrackOut(const AliTRDtrackV1 *track)
{
// Store resolution/pulls of Kalman after last update with the TRD information
// at the radial position of the first tracklet. The following points are used
// for comparison
//  - the (y,z,snp) of the first TRD tracklet
//  - the (y, z, snp, tgl, pt) of the MC track reference
//
// Additionally the momentum resolution/pulls are calculated for usage in the
// PID calculation.

  if(track) fkTrack = track;
  return NULL;
}
*/
//________________________________________________________
TH1* AliTRDresolution::PlotMC(const AliTRDtrackV1 *track)
{
  //
  // Plot MC distributions
  //

  if(!HasMCdata()) return NULL;
  if(track) fkTrack = track;
  if(!fkTrack){
    AliDebug(4, "No Track defined.");
    return NULL;
  }
  Int_t bc(TMath::Abs(fkESD->GetTOFbc()));

  THnSparse *H(NULL);
  if(!fContainer){
    AliWarning("No output container defined.");
    return NULL;
  }
  // retriev track characteristics
  Int_t pdg = fkMC->GetPDG(),
        sIdx(AliTRDpidUtil::Pdg2Pid(TMath::Abs(pdg))+1), // species index
        sign(0),
//        sgm[3],
        label(fkMC->GetLabel());
//        fSegmentLevel(0);
  if(!fDBPDG) fDBPDG=TDatabasePDG::Instance();
  TParticlePDG *ppdg(fDBPDG->GetParticle(pdg));
  if(ppdg) sign = ppdg->Charge() > 0. ? 1 : -1;

  TH1 *h(NULL);
  AliTRDgeometry *geo(AliTRDinfoGen::Geometry());
  AliTRDseedV1 *fTracklet(NULL); TObjArray *clInfoArr(NULL);
  UChar_t s;
  Double_t x, y, z, pt, dydx, dzdx, dzdl;
  Float_t pt0, x0, y0, z0, dx, dy, dz, dydx0, dzdx0;
  Double_t covR[7]/*, cov[3]*/;

/*  if(DebugLevel()>=3){
    // get first detector
    Int_t det = -1;
    for(Int_t ily=0; ily<AliTRDgeometry::kNlayer; ily++){
      if(!(fTracklet = fkTrack->GetTracklet(ily))) continue;
      det = fTracklet->GetDetector();
      break;
    }
    if(det>=0){
      TVectorD X(12), Y(12), Z(12), dX(12), dY(12), dZ(12), vPt(12), dPt(12), budget(12), cCOV(12*15);
      Double_t m(-1.);
      m = fkTrack->GetMass();
      if(fkMC->PropagateKalman(&X, &Y, &Z, &dX, &dY, &dZ, &vPt, &dPt, &budget, &cCOV, m)){
        (*DebugStream()) << "MCkalman"
          << "pdg=" << pdg
          << "det=" << det
          << "x="   << &X
          << "y="   << &Y
          << "z="   << &Z
          << "dx="  << &dX
          << "dy="  << &dY
          << "dz="  << &dZ
          << "pt="  << &vPt
          << "dpt=" << &dPt
          << "bgt=" << &budget
          << "cov=" << &cCOV
          << "\n";
      }
    }
  }*/
  AliTRDcluster *c(NULL);
  Double_t val[kNdim+1];
  for(Int_t ily=0; ily<AliTRDgeometry::kNlayer; ily++){
    if(!(fTracklet = fkTrack->GetTracklet(ily)))/* ||
       !fTracklet->IsOK())*/ continue;

    x= x0 = fTracklet->GetX();
    Bool_t rc(fTracklet->IsRowCross()); Float_t eta, phi;
    if(!fkMC->GetDirections(x0, y0, z0, dydx0, dzdx0, pt0, eta, phi, s)) continue;

    // MC track position at reference radial position
    dx  = x0 - x;
    Float_t ymc = y0 - dx*dydx0;
    Float_t zmc = z0 - dx*dzdx0;
    //phi -= TMath::Pi();

    val[kBC]  = ily;
    val[kPhi] = phi;
    val[kEta] = eta;
    val[kSpeciesChgRC]= rc?0.:sign*sIdx;
    val[kPt]  = pt0<0.8?0:(pt0<1.5?1:2);//GetPtBin(pt0);
    Double_t tilt(fTracklet->GetTilt());
//             ,t2(tilt*tilt)
//             ,corr(1./(1. + t2))
//             ,cost(TMath::Sqrt(corr));

    AliExternalTrackParam *tin(fkTrack->GetTrackIn());
    if(ily==0 && tin){ // trackIn residuals
      // check radial position
      if(TMath::Abs(tin->GetX()-x)>1.e-3) AliDebug(1, Form("TrackIn radial mismatch. dx[cm]=%+4.1f", tin->GetX()-x));
      else{
        val[kBC]          = (bc>=kNbunchCross)?(kNbunchCross-1):bc;
        val[kYrez]        = tin->GetY()-ymc;
        val[kZrez]        = rc?(tin->GetZ()-zmc):(fTracklet->GetdQdl()*1.8e-4 - 0.9);
        val[kPrez]        = (TMath::ASin(tin->GetSnp())-TMath::ATan(dydx0))*TMath::RadToDeg();
        if((H = (THnSparseI*)fContainer->At(kMCtrackIn))) H->Fill(val);
      }
    }
    //if(bc>1) break; // do nothing for the rest of TRD objects if satellite bunch

    // track residuals
    dydx = fTracklet->GetYref(1);
    dzdx = fTracklet->GetZref(1);
    dzdl = fTracklet->GetTgl();
    y  = fTracklet->GetYref(0);
    dy = y - ymc;
    z  = fTracklet->GetZref(0);
    dz = z - zmc;
    pt = TMath::Abs(fTracklet->GetPt());
    fTracklet->GetCovRef(covR);

    val[kYrez] = dy;
    val[kPrez] = TMath::ATan((dydx - dydx0)/(1.+ dydx*dydx0))*TMath::RadToDeg();
    val[kZrez] = dz;
    val[kNdim] = 1.e2*(pt/pt0-1.);
    if((H = (THnSparse*)fContainer->At(kMCtrack))) H->Fill(val);
/*      // theta resolution/ tgl pulls
      Double_t dzdl0 = dzdx0/TMath::Sqrt(1.+dydx0*dydx0),
                dtgl = (dzdl - dzdl0)/(1.- dzdl*dzdl0);
      ((TH2I*)arr->At(6))->Fill(dzdl0, TMath::ATan(dtgl));
      ((TH2I*)arr->At(7))->Fill(dzdl0, (dzdl - dzdl0)/TMath::Sqrt(covR[4]));
      // pt resolution  \\ 1/pt pulls \\ p resolution for PID
      Double_t p0 = TMath::Sqrt(1.+ dzdl0*dzdl0)*pt0,
              p  = TMath::Sqrt(1.+ dzdl*dzdl)*pt;
      ((TH3S*)((TObjArray*)arr->At(8)))->Fill(pt0, pt/pt0-1., sign*sIdx);
      ((TH3S*)((TObjArray*)arr->At(9)))->Fill(1./pt0, (1./pt-1./pt0)/TMath::Sqrt(covR[6]), sign*sIdx);
      ((TH3S*)((TObjArray*)arr->At(10)))->Fill(p0, p/p0-1., sign*sIdx);*/

    // Fill Debug stream for MC track
    if(DebugLevel()>=4){
      Int_t det(fTracklet->GetDetector());
      (*DebugStream()) << "MC"
        << "det="     << det
        << "pdg="     << pdg
        << "sgn="     << sign
        << "pt="      << pt0
        << "x="       << x0
        << "y="       << y0
        << "z="       << z0
        << "dydx="    << dydx0
        << "dzdx="    << dzdx0
        << "\n";

      // Fill Debug stream for Kalman track
      (*DebugStream()) << "MCtrack"
        << "pt="      << pt
        << "x="       << x
        << "y="       << y
        << "z="       << z
        << "dydx="    << dydx
        << "dzdx="    << dzdx
        << "s2y="     << covR[0]
        << "s2z="     << covR[2]
        << "\n";
    }

    // tracklet residuals
    dydx = fTracklet->GetYfit(1) + tilt*dzdx0;
    dzdx = fTracklet->GetZfit(1);
    y  = fTracklet->GetYfit(0);
    dy = y - ymc;
    z  = fTracklet->GetZfit(0);
    dz = z - zmc;
    val[kYrez] = dy - dz*tilt;
    val[kPrez] = TMath::ATan((dydx - dydx0)/(1.+ dydx*dydx0))*TMath::RadToDeg();
    val[kZrez] = rc?(dz + dy*tilt):(fTracklet->GetdQdl()*3.e-4 - 1.5);
//      val[kNdim] = pt/pt0-1.;
    if((H = (THnSparse*)fContainer->At(kMCtracklet))) H->Fill(val);


    // Fill Debug stream for tracklet
    if(DebugLevel()>=4){
      Float_t s2y = fTracklet->GetS2Y();
      Float_t s2z = fTracklet->GetS2Z();
      (*DebugStream()) << "MCtracklet"
        << "rc="    << rc
        << "x="     << x
        << "y="     << y
        << "z="     << z
        << "dydx="  << dydx
        << "s2y="   << s2y
        << "s2z="   << s2z
        << "\n";
    }

    AliTRDpadPlane *pp = geo->GetPadPlane(ily, AliTRDgeometry::GetStack(fTracklet->GetDetector()));
    Float_t zr0 = pp->GetRow0() + pp->GetAnodeWireOffset();
    //Double_t exb = AliTRDCommonParam::Instance()->GetOmegaTau(1.5);

    H = (THnSparse*)fContainer->At(kMCcluster);
    val[kPt]  = TMath::ATan(dydx0)*TMath::RadToDeg();
    //Float_t corr = 1./TMath::Sqrt(1.+dydx0*dydx0+dzdx0*dzdx0);
    Int_t row0(-1);
    Float_t padCorr(tilt*fTracklet->GetPadLength());
    fTracklet->ResetClusterIter(kTRUE);
    while((c = fTracklet->NextCluster())){
      if(row0<0) row0 = c->GetPadRow();
      x = c->GetX();//+fXcorr[c->GetDetector()][c->GetLocalTimeBin()];
      y = c->GetY()  + padCorr*(c->GetPadRow() - row0);
      z = c->GetZ();
      dx = x0 - x;
      ymc= y0 - dx*dydx0;
      zmc= z0 - dx*dzdx0;
      dy = y - ymc;
      dz = z - zmc;
      val[kYrez] = dy - dz*tilt;
      val[kPrez] = dx;
      val[kZrez] = 0.; AliTRDcluster *cc(NULL); Int_t ic(0), tb(c->GetLocalTimeBin()); Float_t  q(TMath::Abs(c->GetQ()));
      if((cc = fTracklet->GetClusters(tb-1))) {val[kZrez] += TMath::Abs(cc->GetQ()); ic++;}
      if((cc = fTracklet->GetClusters(tb-2))) {val[kZrez] += TMath::Abs(cc->GetQ()); ic++;}
      if(ic) val[kZrez] /= (ic*q);
      if(H) H->Fill(val);


      // Fill calibration container
      Float_t d = zr0 - zmc;
      d -= ((Int_t)(2 * d)) / 2.0;
      if (d > 0.25) d  = 0.5 - d;
      AliTRDclusterInfo *clInfo = new AliTRDclusterInfo;
      clInfo->SetCluster(c);
      clInfo->SetMC(pdg, label);
      clInfo->SetGlobalPosition(ymc, zmc, dydx0, dzdx0);
      clInfo->SetResolution(dy);
      clInfo->SetAnisochronity(d);
      clInfo->SetDriftLength(dx);
      clInfo->SetTilt(tilt);
      if(fMCcl) fMCcl->Add(clInfo);
      else AliDebug(1, "MCcl exchange container missing. Activate by calling \"InitExchangeContainers()\"");
      if(DebugLevel()>=5){
        if(!clInfoArr){
          clInfoArr=new TObjArray(AliTRDseedV1::kNclusters);
          clInfoArr->SetOwner(kFALSE);
        }
        clInfoArr->Add(clInfo);
      }
    }
    // Fill Debug Tree
    if(DebugLevel()>=5 && clInfoArr){
      (*DebugStream()) << "MCcluster"
        <<"clInfo.=" << clInfoArr
        << "\n";
      clInfoArr->Clear();
    }
  }
  if(clInfoArr) delete clInfoArr;
  return h;
}


//__________________________________________________________________________
Int_t AliTRDresolution::GetPtBin(Float_t pt)
{
// Find pt bin according to local pt segmentation
  Int_t ipt(-1);
  while(ipt<24){
    if(pt<fgPtBin[ipt+1]) break;
    ipt++;
  }
  return ipt;
}

//________________________________________________________
Float_t AliTRDresolution::GetMeanStat(TH1 *h, Float_t cut, Option_t *opt)
{
// return mean number of entries/bin of histogram "h"
// if option "opt" is given the following values are accepted:
//   "<" : consider only entries less than "cut"
//   ">" : consider only entries greater than "cut"

  //Int_t dim(h->GetDimension());
  Int_t nbx(h->GetNbinsX()), nby(h->GetNbinsY()), nbz(h->GetNbinsZ());
  Double_t sum(0.); Int_t n(0);
  for(Int_t ix(1); ix<=nbx; ix++)
    for(Int_t iy(1); iy<=nby; iy++)
      for(Int_t iz(1); iz<=nbz; iz++){
        if(strcmp(opt, "")==0){sum += h->GetBinContent(ix, iy, iz); n++;}
        else{
          if(strcmp(opt, "<")==0) {
            if(h->GetBinContent(ix, iy, iz)<cut) {sum += h->GetBinContent(ix, iy, iz); n++;}
          } else if(strcmp(opt, ">")==0){
            if(h->GetBinContent(ix, iy, iz)>cut) {sum += h->GetBinContent(ix, iy, iz); n++;}
          } else {sum += h->GetBinContent(ix, iy, iz); n++;}
        }
      }
  return n>0?sum/n:0.;
}

//________________________________________________________
void AliTRDresolution::GetRangeZ(TH2 *h2, Float_t &min, Float_t &max)
{
// Get range on Z axis such to avoid outliers

  Double_t cnt[20000], c, m, s;
  Int_t nx(h2->GetXaxis()->GetNbins()), ny(h2->GetYaxis()->GetNbins()), nc(0);
  for(Int_t ix(1); ix<=nx; ix++){
    for(Int_t iy(1); iy<=ny; iy++){
      if((c = h2->GetBinContent(ix, iy))<10) continue;
      cnt[nc++] = c;
      if(nc==20000) break;
    }
    if(nc==20000) break;
  }
  AliMathBase::EvaluateUni(nc, cnt, m, s, 0);
  min = m-s; max = m+2.*s;
}

//________________________________________________________
Bool_t AliTRDresolution::GetRefFigure(Int_t ifig)
{
  //
  // Get the reference figures
  //

  if(!gPad){
    AliWarning("Please provide a canvas to draw results.");
    return kFALSE;
  }
/*  Int_t selection[100], n(0), selStart(0); //
  Int_t ly0(0), dly(5);
  TList *l(NULL); TVirtualPad *pad(NULL); */
  switch(ifig){
  case 0:
    break;
  }
  AliWarning(Form("Reference plot [%d] missing result", ifig));
  return kFALSE;
}


//________________________________________________________
void AliTRDresolution::MakePtSegmentation(Float_t pt0, Float_t dpt)
{
// Build pt segments
  for(Int_t j(0); j<=24; j++){
    pt0+=(TMath::Exp(j*j*dpt)-1.);
    fgPtBin[j]=pt0;
  }
}

//________________________________________________________
void AliTRDresolution::MakeSummary()
{
// Build summary plots

  if(!fProj){
    AliError("Missing results");
    return;
  }
  TVirtualPad *p(NULL); TCanvas *cOut(NULL);
  TObjArray *arr(NULL); TH2 *h2(NULL);

  // cluster resolution
  // define palette
  gStyle->SetPalette(1);
  const Int_t nClViews(9);
  const Char_t *vClName[nClViews] = {"ClY", "ClYn", "ClYp", "ClQn", "ClQp", "ClYXTCp", "ClYXTCn", "ClYXPh", "ClYXPh"};
  const UChar_t vClOpt[nClViews] = {1, 1, 1, 0, 0, 0, 0, 0, 1};
  const Float_t vClSignMin[2] = {2.6e2, 4.4e2},
                vClSignMax[2] = {4.4e2, 6.2e2},
                vClMin[nClViews] = {3.2e2, vClSignMin[Int_t(fBsign)], vClSignMin[Int_t(!fBsign)], 0., 0., 0., 0., 0., 3.2e2},
                vClMax[nClViews] = {5.e2, vClSignMax[Int_t(fBsign)], vClSignMax[Int_t(!fBsign)], 0., 0., 0., 0., 0., 5.e2};
  const Int_t nTrkltViews(20);
  const Char_t *vTrkltName[nTrkltViews] = {
    "TrkltY", "TrkltYn", "TrkltYp", "TrkltY", "TrkltYn", "TrkltYp",
    "TrkltRCZ",
    "TrkltPh", "TrkltPhn", "TrkltPhp",
    "TrkltQ", "TrkltQn", "TrkltQp",
    "TrkltQS", "TrkltQSn", "TrkltQSp",
    "TrkltTBn", "TrkltTBp", "TrkltTBn", "TrkltTBp"
//    "TrkltYnl0", "TrkltYpl0", "TrkltPhnl0", "TrkltPhpl0", "TrkltQnl0", "TrkltQpl0",  // electrons low pt
/*    "TrkltYnl1", "TrkltYpl1", "TrkltPhnl1", "TrkltPhpl1", "TrkltQnl1", "TrkltQpl1",  // muons low pt
    "TrkltYnl2", "TrkltYpl2", "TrkltPhnl2", "TrkltPhpl2", "TrkltQnl2", "TrkltQpl2"  // pions low pt*/
  };
  const UChar_t vTrkltOpt[nTrkltViews] = {
    0, 0, 0, 1, 1, 1,
    0,
    0, 0, 0,
    0, 0, 0,
    0, 0, 0,
    0, 0, 2, 2
  };
  const Int_t nTrkInViews(5);
  const Char_t *vTrkInName[nTrkInViews][6] = {
    {"TrkInY", "TrkInYn", "TrkInYp", "TrkInRCZ", "TrkInPhn", "TrkInPhp"},
    {"TrkInY", "TrkInYn", "TrkInYp", "TrkInRCZ", "TrkInPhn", "TrkInPhp"},
    {"TrkInYnl", "TrkInYni", "TrkInYnh", "TrkInYpl", "TrkInYpi", "TrkInYph"},
    {"TrkInXnl", "TrkInXpl", "TrkInXl", "TrkInYnh", "TrkInYph", "TrkInYh"},
    {"TrkInPhnl", "TrkInPhni", "TrkInPhnh", "TrkInPhpl", "TrkInPhpi", "TrkInPhph"},
    //{"TrkInRCX", "TrkInRCY", "TrkInRCPh", "TrkInRCZl", "TrkInRCZi", "TrkInRCZh"}
  };
  const UChar_t vTrkInOpt[nTrkInViews] = {0, 1, 0, 0, 0};
  const Float_t min[6] = {0.15, 0.15, 0.15, 0.15, 0.5, 0.5};
  const Float_t max[6] = {0.6, 0.6, 0.6, 0.6, 2.3, 2.3};
  const Char_t *ttt[6] = {"#sigma(#Deltay) [cm]", "#sigma(#Deltay) [cm]", "#sigma(#Deltay) [cm]", "#sigma(#Deltaz) [cm]", "#sigma(#Delta#phi) [deg]", "#sigma(#Delta#phi) [deg]"};

  const Int_t nTrkViews(27);
  const Char_t *vTrkName[nTrkViews] = {
    "TrkY", "TrkYn", "TrkYp",
    "TrkPh", "TrkPhn", "TrkPhp",
    "TrkDPt", "TrkDPtn", "TrkDPtp",
    "TrkYnl0", "TrkYpl0", "TrkPhnl0", "TrkPhpl0", "TrkDPtnl0", "TrkDPtpl0",  // electrons low pt
    "TrkYnl1", "TrkYpl1", "TrkPhnl1", "TrkPhpl1", "TrkDPtnl1", "TrkDPtpl1",  // muons low pt
    "TrkYnl2", "TrkYpl2", "TrkPhnl2", "TrkPhpl2", "TrkDPtnl2", "TrkDPtpl2"  // pions low pt
  };
  const Char_t *typName[] = {"", "MC"};
  const Int_t nx(2048), ny(1536);

  if((arr = (TObjArray*)fProj->At(kDetector))){
    cOut = new TCanvas(Form("%s_DetOccupancy", GetName()), "Detector performance", 2*nx, 2*ny);
    cOut->Divide(AliTRDgeometry::kNlayer,AliTRDeventInfo::kCentralityClasses, 1.e-5, 1.e-5);
    for(Int_t icen(0); icen<AliTRDeventInfo::kCentralityClasses; icen++){
      for(Int_t ily(0); ily<AliTRDgeometry::kNlayer; ily++){
        p=cOut->cd(icen*AliTRDgeometry::kNlayer+ily+1); p->SetRightMargin(0.1572581);p->SetTopMargin(0.08262712);
        if(!(h2 = (TH2*)arr->FindObject(Form("HDet%d%dEn", ily, icen)))) continue;
        SetNormZ(h2, 1, 11, 1, -1, 10.);
        SetRangeZ(h2, 0., 150, -25.);
        h2->GetZaxis()->SetTitle("Rel. Det. Occup. [%]");
        h2->GetZaxis()->CenterTitle();
        h2->Draw("colz");
        MakeDetectorPlot(ily, "pad");
      }
    }
    cOut->SaveAs(Form("%s.gif", cOut->GetName()));
    cOut = new TCanvas(Form("%s_DetCharge", GetName()), "Detector performance", nx, ny);
    cOut->Divide(3,2, 1.e-5, 1.e-5);
    for(Int_t ily(0); ily<AliTRDgeometry::kNlayer; ily++){
      p=cOut->cd(ily+1); p->SetRightMargin(0.1572581);p->SetTopMargin(0.08262712);
      if(!(h2 = (TH2*)arr->FindObject(Form("HDet%d_2D", ily)))) continue;
      SetNormZ(h2, 1, 11, 1, -1, 10.);
      SetRangeZ(h2, 0., 45., -15.);
      h2->GetZaxis()->SetTitle("Rel. Mean(q) [%]");
      h2->GetZaxis()->CenterTitle();
      h2->Draw("colz");
      MakeDetectorPlot(ily, "pad");
    }
    cOut->SaveAs(Form("%s.gif", cOut->GetName()));
  }
  for(Int_t ityp(0); ityp<(HasMCdata()?2:1); ityp++){
    if((arr = (TObjArray*)fProj->At(ityp?kMCcluster:kCluster))){
      for(Int_t iview(0); iview<nClViews; iview++){
        cOut = new TCanvas(Form("%s_%s%s_%d", GetName(), typName[ityp], vClName[iview], vClOpt[iview]), "Cluster Resolution", nx, ny);
        cOut->Divide(3,2, 1.e-5, 1.e-5);
        Int_t nplot(0);
        for(Int_t ily(0); ily<AliTRDgeometry::kNlayer; ily++){
          p=cOut->cd(ily+1);    p->SetRightMargin(0.1572581);p->SetTopMargin(0.08262712);
          if(!(h2 = (TH2*)arr->FindObject(Form("H%s%s%d_2D", typName[ityp], vClName[iview], ily)))) continue;
          nplot++;
          if(vClOpt[iview]==0) h2->Draw("colz");
          else if(vClOpt[iview]==1) DrawSigma(h2, "#sigma(#Deltay) [#mum]", vClMin[iview], vClMax[iview], 1.e4);
          if(iview<5) MakeDetectorPlot(ily);
        }
        if(nplot==AliTRDgeometry::kNlayer) cOut->SaveAs(Form("%s.gif", cOut->GetName()));
        else delete cOut;
      }
    }
    // tracklet systematic
    if((arr = (TObjArray*)fProj->At(ityp?kMCtracklet:kTracklet))){
      for(Int_t iview(0); iview<nTrkltViews; iview++){
        cOut = new TCanvas(Form("%s_%s%s_%d", GetName(), typName[ityp], vTrkltName[iview], vTrkltOpt[iview]), "Tracklet Resolution", nx, ny);
        cOut->Divide(3,2, 1.e-5, 1.e-5);
        Int_t nplot(0);
        for(Int_t iplot(0); iplot<6; iplot++){
          p=cOut->cd(iplot+1); p->SetRightMargin(0.1572581); p->SetTopMargin(0.08262712);
          if(!(h2 = (TH2*)arr->FindObject(Form("H%s%s%d_2D", typName[ityp], vTrkltName[iview], iplot)))) continue;
          nplot++;
          if(vTrkltOpt[iview]==0) h2->Draw("colz");
          else if (vTrkltOpt[iview]==1) DrawSigma(h2, "#sigma(#Deltay) [cm]", .15, .4);
          else if (vTrkltOpt[iview]==2) DrawSigma(h2, "#sigma(occupancy) [%]", 10.5, 15.);
          MakeDetectorPlot(iplot);
        }
        if(nplot==6){
          cOut->Modified();cOut->Update();
          cOut->SaveAs(Form("%s.gif", cOut->GetName()));
        }
        delete cOut;
      }
    }
    // trackIn systematic
    if((arr = (TObjArray*)fProj->At(ityp?kMCtrackIn:kTrackIn))){
      for(Int_t iview(0); iview<nTrkInViews; iview++){
        cOut = new TCanvas(Form("%s_%s%s_%d", GetName(), typName[ityp], vTrkInName[iview][0], vTrkInOpt[iview]), "Track IN Resolution", nx, ny);
        cOut->Divide(3,2, 1.e-5, 1.e-5);
        Int_t nplot(0);
        for(Int_t iplot(0); iplot<6; iplot++){
          p=cOut->cd(iplot+1); p->SetRightMargin(0.1572581); p->SetTopMargin(0.08262712);
          if(!(h2 = (TH2*)arr->FindObject(Form("H%s%s_2D", typName[ityp], vTrkInName[iview][iplot])))){
            AliInfo(Form("Missing H%s%s_2D", typName[ityp], vTrkInName[iview][iplot]));
            continue;
          }
          nplot++;
          if(vTrkInOpt[iview]==0) h2->Draw("colz");
          else DrawSigma(h2, ttt[iplot], min[iplot], max[iplot]);
          MakeDetectorPlot(0);
        }
        if(nplot==6) cOut->SaveAs(Form("%s.gif", cOut->GetName()));
        else delete cOut;
      }
      // species
      const Float_t zmin[] = {1., 20., 25., 10., 5.},
                    zmax[] = {10., 38., 43., 28., 14.};
      cOut = new TCanvas(Form("%s_%sTrkInSpc", GetName(), typName[ityp]), "Track IN PID", Int_t(1.5*nx), Int_t(1.5*ny));
      cOut->Divide(5,3, 1.e-5, 1.e-5);
      Int_t nplot(0); const Char_t *chName[] = {"p", "n", ""};
      for(Int_t ich(0), ipad(1); ich<3; ich++){
        TH2 *h2s(NULL);
        if(!(h2s = (TH2*)arr->FindObject(Form("H%sTrkInY%sEn", typName[ityp], chName[ich])))) {
          AliInfo(Form("Missing H%sTrkIn%sEn", typName[ityp], chName[ich]));
          continue;
        }
        for(Int_t ispec(0); ispec<AliPID::kSPECIES; ispec++){
          p=cOut->cd(ipad++); p->SetRightMargin(0.1572581); p->SetTopMargin(0.08262712);
          if(!(h2 = (TH2*)arr->FindObject(Form("H%sTrkInY%s%dEn", typName[ityp], chName[ich], ispec)))) {
            AliInfo(Form("Missing H%sTrkIn%s%dEn", typName[ityp], chName[ich], ispec));
            continue;
          }
          nplot++;
          h2->Divide(h2, h2s, 1.e2);
          h2->SetContour(9);
          h2->GetZaxis()->SetRangeUser(zmin[ispec], zmax[ispec]);
          h2->GetZaxis()->SetTitle("Rel. Abundancy [%]");h2->GetZaxis()->CenterTitle();
          TString tit(h2->GetTitle()); h2->SetTitle(Form("%s :: Relative Abundancy", ((*(tit.Tokenize("::")))[0])->GetName()));
          h2->Draw("colz");
          MakeDetectorPlot(0);
        }
      }
      if(nplot==15) cOut->SaveAs(Form("%s.gif", cOut->GetName()));
      else delete cOut;
    }
  }
  // track MC systematic
  if((arr = (TObjArray*)fProj->At(kMCtrack))) {
    for(Int_t iview(0); iview<nTrkViews; iview++){
      cOut = new TCanvas(Form("%s_MC%s", GetName(), vTrkName[iview]), "Track Resolution", nx, ny);
      cOut->Divide(3,2, 1.e-5, 1.e-5);
      Int_t nplot(0);
      for(Int_t iplot(0); iplot<6; iplot++){
        p=cOut->cd(iplot+1); p->SetRightMargin(0.1572581); p->SetTopMargin(0.08262712);
        if(!(h2 = (TH2*)arr->FindObject(Form("HMC%s%d_2D", vTrkName[iview], iplot)))) continue;
        h2->Draw("colz"); nplot++;
      }
      if(nplot==6) cOut->SaveAs(Form("%s.gif", cOut->GetName()));
      else delete cOut;
    }
  }


  gStyle->SetPalette(1);
}

//________________________________________________________
void AliTRDresolution::DrawSigma(TH2 *h2, const Char_t *title, Float_t m, Float_t M, Float_t scale)
{
  // Draw error bars scaled with "scale" instead of content values
  //use range [m,M] if limits are specified

  if(!h2) return;
  TAxis *ax(h2->GetXaxis()), *ay(h2->GetYaxis());
  TH2F *h2e = new TH2F(Form("%s_E", h2->GetName()),
                Form("%s;%s;%s;%s", h2->GetTitle(), ax->GetTitle(), ay->GetTitle(), title),
                ax->GetNbins(), ax->GetXmin(), ax->GetXmax(),
                ay->GetNbins(), ay->GetXmin(), ay->GetXmax());
  h2e->SetContour(9);
  TAxis *az(h2e->GetZaxis());
  if(M>m) az->SetRangeUser(m, M);
  az->CenterTitle();
  az->SetTitleOffset(1.5);
  for(Int_t ix(1); ix<=h2->GetNbinsX(); ix++){
    for(Int_t iy(1); iy<=h2->GetNbinsY(); iy++){
      if(h2->GetBinContent(ix, iy)<-100.) continue;
      Float_t v(scale*h2->GetBinError(ix, iy));
      if(M>m && v<m) v=m+TMath::Abs((M-m)*1.e-3);
      h2e->SetBinContent(ix, iy, v);
    }
  }
  h2e->Draw("colz");
}

//________________________________________________________
void AliTRDresolution::GetRange(TH2 *h2, Char_t mod, Float_t *range)
{
// Returns the range of the bulk of data in histogram h2. Removes outliers.
// The "range" vector should be initialized with 2 elements
// Option "mod" can be any of
//   - 0 : gaussian like distribution
//   - 1 : tailed distribution

  Int_t nx(h2->GetNbinsX())
      , ny(h2->GetNbinsY())
      , n(nx*ny);
  Double_t *data=new Double_t[n];
  for(Int_t ix(1), in(0); ix<=nx; ix++){
    for(Int_t iy(1); iy<=ny; iy++)
      data[in++] = h2->GetBinContent(ix, iy);
  }
  Double_t mean, sigm;
  AliMathBase::EvaluateUni(n, data, mean, sigm, Int_t(n*.8));

  range[0]=mean-3.*sigm; range[1]=mean+3.*sigm;
  if(mod==1) range[0]=TMath::Max(Float_t(1.e-3), range[0]);
  AliDebug(2, Form("h[%s] range0[%f %f]", h2->GetName(), range[0], range[1]));
  TH1S h1("h1SF0", "", 100, range[0], range[1]);
  h1.FillN(n,data,0);
  delete [] data;

  switch(mod){
  case 0:// gaussian distribution
  {
    TF1 fg("fg", "gaus", mean-3.*sigm, mean+3.*sigm);
    h1.Fit(&fg, "QN");
    mean=fg.GetParameter(1); sigm=fg.GetParameter(2);
    range[0] = mean-2.5*sigm;range[1] = mean+2.5*sigm;
    AliDebug(2, Form("     rangeG[%f %f]", range[0], range[1]));
    break;
  }
  case 1:// tailed distribution
  {
    Int_t bmax(h1.GetMaximumBin());
    Int_t jBinMin(1), jBinMax(100);
    for(Int_t ibin(bmax); ibin--;){
      if(h1.GetBinContent(ibin)<1.){
        jBinMin=ibin; break;
      }
    }
    for(Int_t ibin(bmax); ibin++;){
      if(h1.GetBinContent(ibin)<1.){
        jBinMax=ibin; break;
      }
    }
    range[0]=h1.GetBinCenter(jBinMin); range[1]=h1.GetBinCenter(jBinMax);
    AliDebug(2, Form("     rangeT[%f %f]", range[0], range[1]));
    break;
  }
  }

  return;
}

//________________________________________________________
Bool_t AliTRDresolution::MakeProjectionDetector()
{
// Analyse cluster
  const Int_t kNcontours(9);
  const Int_t kNstat(100);
  if(fProj && fProj->At(kDetector)) return kTRUE;
  if(!fContainer){
    AliError("Missing data container.");
    return kFALSE;
  }
  THnSparse *H(NULL);
  if(!(H = (THnSparse*)fContainer->FindObject("hDet2Cluster"))){
    AliInfo(Form("Missing/Wrong data @ hDet2Cluster."));
    return kTRUE;
  }
  Int_t ndim(H->GetNdimensions());
  Int_t coord[kNdim]; memset(coord, 0, sizeof(Int_t) * kNdim); Double_t v = 0.;
  TAxis *aa[kNdim], *an(NULL); memset(aa, 0, sizeof(TAxis*) * kNdim);
  for(Int_t id(0); id<ndim; id++) aa[id] = H->GetAxis(id);
  if(ndim < 5) aa[4] = new TAxis(1, -0.5, 0.5);
  Int_t nPad(1);
  if(ndim > 5){
    an = H->GetAxis(5);
    nPad = kNpads+1;
  }
  // build list of projections
  const Int_t nsel(8*AliTRDgeometry::kNlayer*AliTRDeventInfo::kCentralityClasses);
  // define rebinning strategy
  const Int_t nEtaPhi(4); Int_t rebinEtaPhiX[nEtaPhi] = {1, 2, 2, 1}, rebinEtaPhiY[nEtaPhi] = {2, 1, 1, 5};
  //const Char_t *cenName[AliTRDeventInfo::kCentralityClasses] = {"0-10%", "10-20%", "20-50%", "50-80%", "80-100%"};
  const Char_t *cenName[AliTRDeventInfo::kCentralityClasses] = {"2800-inf", "2100-2799", "1400-2099", "700-1399", "0-699"};
  AliTRDresolutionProjection hp[kDetNproj];  TObjArray php(kDetNproj);
  Int_t ih(0), isel(-1), np[nsel]; memset(np, 0, nsel*sizeof(Int_t));
  for(Int_t ipad(0); ipad<nPad; ipad++){
    for(Int_t icen(0); icen<AliTRDeventInfo::kCentralityClasses; icen++){
      for(Int_t ily(0); ily<AliTRDgeometry::kNlayer; ily++){
        isel++; // new selection
        hp[ih].Build(Form("HDet%d%d%d", ily, icen, ipad),
                    Form("Detectors :: Ly[%d] Cen[%s] Pads[%s]", ily, cenName[icen], ipad?(ipad<kNpads?Form("%d", ipad+1):Form(">%d", kNpads)):"deconv"),
                    kEta, kPhi, 4, aa);
        hp[ih].SetRebinStrategy(nEtaPhi, rebinEtaPhiX, rebinEtaPhiY);
        hp[ih].SetShowRange(10., 55.);
        php.AddLast(&hp[ih++]); np[isel]++;
      }
    }
  }
  AliInfo(Form("Build %3d 3D projections.", ih));

  Int_t ly(0), cen(0), npad(0);
  for (Long64_t ib(0); ib < H->GetNbins(); ib++) {
    v = H->GetBinContent(ib, coord); if(v<1.) continue;
    ly = coord[kBC]-1;    // layer selection
    cen = coord[kYrez]-1; // centrality selection
    npad = 0;             // no. of pads selection
    if(an) npad = TMath::Min(coord[5]-1, Int_t(kNpads));
    isel = npad*AliTRDeventInfo::kCentralityClasses*AliTRDgeometry::kNlayer+cen*AliTRDgeometry::kNlayer+ly;
    ((AliTRDresolutionProjection*)php.At(isel))->Increment(coord, v);
    //Int_t ioff=isel;for(Int_t jh(0); jh<np[isel]; jh++) ((AliTRDresolutionProjection*)php.At(ioff+jh))->Increment(coord, v);
  }
  TObjArray *arr(NULL);
  fProj->AddAt(arr = new TObjArray(kDetNproj), kDetector);

  TH2 *h2(NULL);  Int_t jh(0);
  for(; ih--; ){
    if(!hp[ih].fH) continue;
    if(!(h2 = hp[ih].Projection2D(kNstat, kNcontours, 0, kFALSE))) continue;
    arr->AddAt(h2, jh++);
    if(!(h2 = (TH2*)gDirectory->Get(Form("%sEn", hp[ih].fH->GetName())))) continue;
    arr->AddAt(h2, jh++);
  }
  AliTRDresolutionProjection *pr0(NULL), *pr1(NULL);
  for(Int_t ily(0); ily<AliTRDgeometry::kNlayer; ily++){
    for(Int_t icen(0); icen<AliTRDeventInfo::kCentralityClasses; icen++){
      if((pr0 = (AliTRDresolutionProjection*)php.FindObject(Form("HDet%d%d%d", ily, icen, 0)))){
        for(Int_t ipad(1); ipad<nPad; ipad++){
          if((pr1 = (AliTRDresolutionProjection*)php.FindObject(Form("HDet%d%d%d", ily, icen, ipad)))){
            (*pr0)+=(*pr1);
          }
        }
        pr0->fH->SetNameTitle(Form("HDet%d%d", ily, icen), Form("Detectors :: Ly[%d] Cen[%s]", ily, cenName[icen]));
        if((h2 = pr0->Projection2D(kNstat, kNcontours, 0, kFALSE))) arr->AddAt(h2, jh++);
        if((h2 = (TH2*)gDirectory->Get(Form("%sEn", pr0->fH->GetName())))) arr->AddAt(h2, jh++);
        if(icen && (pr1 = (AliTRDresolutionProjection*)php.FindObject(Form("HDet%d%d%d", ily, 0, 0)))) (*pr1)+=(*pr0);
      }
    }
    if((pr0 = (AliTRDresolutionProjection*)php.FindObject(Form("HDet%d%d%d", ily, 0, 0)))){
      pr0->fH->SetNameTitle(Form("HDet%d", ily), Form("Detectors :: Ly[%d]", ily));
      if((h2 = pr0->Projection2D(kNstat, kNcontours))) arr->AddAt(h2, jh++);
    }
  }
  AliInfo(Form("Done %3d 2D projections.", jh));
  return kTRUE;
}

//________________________________________________________
Bool_t AliTRDresolution::MakeProjectionCluster(Bool_t mc)
{
// Analyse cluster
  const Int_t kNcontours(9);
  const Int_t kNstat(100);
  Int_t cidx=mc?kMCcluster:kCluster;
  if(fProj && fProj->At(cidx)) return kTRUE;
  if(!fContainer){
    AliError("Missing data container.");
    return kFALSE;
  }
  const Char_t *projName[] = {"hCluster2Track", "hCluster2MC"};
  THnSparse *H(NULL);
  if(!(H = (THnSparse*)fContainer->FindObject(projName[Int_t(mc)]))){
    AliError(Form("Missing/Wrong data @ %s.", projName[Int_t(mc)]));
    return kFALSE;
  }
  Int_t ndim(H->GetNdimensions()); Bool_t debug(ndim>Int_t(kNdimCl));
  Int_t coord[kNdim]; memset(coord, 0, sizeof(Int_t) * kNdim); Double_t v = 0.;
  TAxis *aa[kNdim], *as(NULL), *apt(NULL), *acen(NULL), *anp(NULL); memset(aa, 0, sizeof(TAxis*) * kNdim);
  for(Int_t id(0); id<ndim; id++) aa[id] = H->GetAxis(id);
  if(ndim > Int_t(kPt)) apt = H->GetAxis(kPt);
  if(ndim > Int_t(kSpeciesChgRC)) as  = H->GetAxis(kSpeciesChgRC);
  if(ndim > Int_t(kNdim)) acen  = H->GetAxis(kNdim);
  if(ndim > Int_t(kNdim)+1) anp  = H->GetAxis(kNdim+1);
  // calculate size depending on debug level
  const Int_t nCh(apt?2:1);
  const Int_t nCen(acen?Int_t(AliTRDeventInfo::kCentralityClasses):1);
  const Int_t nNpad(anp?(Int_t(kNpads)+1):1);
  
  // build list of projections
  const Int_t nsel(AliTRDgeometry::kNlayer*kNcharge*(kNpads+1)*AliTRDeventInfo::kCentralityClasses);
  // define rebinning strategy
  const Int_t nEtaPhi(4); Int_t rebinEtaPhiX[nEtaPhi] = {1, 2, 5, 1}, rebinEtaPhiY[nEtaPhi] = {2, 1, 1, 5};
  AliTRDresolutionProjection hp[kClNproj];  TObjArray php(kClNproj);
  const Char_t chName[kNcharge] = {'n', 'p'};const Char_t chSgn[kNcharge] = {'-', '+'};
  const Char_t *cenName[AliTRDeventInfo::kCentralityClasses] = {"2800-inf", "2100-2799", "1400-2099", "700-1399", "0-699"};
  Int_t ih(0), isel(-1), np[nsel]; memset(np, 0, nsel*sizeof(Int_t));
  for(Int_t ily(0); ily<AliTRDgeometry::kNlayer; ily++){
    for(Int_t ich(0); ich<nCh; ich++){
      for(Int_t icen(0); icen<nCen; icen++){
        for(Int_t ipad(0); ipad<nNpad; ipad++){
          isel++; // new selection
          hp[ih].Build(Form("H%sClY%c%d%d%d", mc?"MC":"", chName[ich], ily, icen, ipad),
                      Form("Clusters[%c] :: #Deltay Ly[%d] Cen[%s] Pads[%s]", chSgn[ich], ily, cenName[icen], ipad?(ipad<kNpads?Form("%d", ipad+1):Form(">%d", kNpads)):"deconv"),
                      kEta, kPhi, kYrez, aa);
          hp[ih].SetRebinStrategy(nEtaPhi, rebinEtaPhiX, rebinEtaPhiY);
          php.AddLast(&hp[ih++]); np[isel]++;
          if(!debug) break;
          hp[ih].Build(Form("H%sClQ%c%d%d%d", mc?"MC":"", chName[ich], ily, icen, ipad),
                      Form("Clusters[%c] :: Q Ly[%d] Cen[%s] Pads[%s]", chSgn[ich], ily, cenName[icen], ipad?(ipad<kNpads?Form("%d", ipad+1):Form(">%d", kNpads)):"deconv"),
                      kEta, kPhi, kSpeciesChgRC, aa);
          hp[ih].SetRebinStrategy(nEtaPhi, rebinEtaPhiX, rebinEtaPhiY);
          hp[ih].SetShowRange(24., 33.);
          php.AddLast(&hp[ih++]); np[isel]++;
          hp[ih].Build(Form("H%sClYXTC%c%d%d%d", mc?"MC":"", chName[ich], ily, icen, ipad),
                      Form("Clusters[%c] :: #Deltay(x,TC) Ly[%d] Cen[%s] Pads[%s]", chSgn[ich], ily, cenName[icen], ipad?(ipad<kNpads?Form("%d", ipad+1):Form(">%d", kNpads)):"deconv"),
                      kPrez, kZrez, kYrez, aa);
          php.AddLast(&hp[ih++]); np[isel]++;
          hp[ih].Build(Form("H%sClYXPh%c%d%d%d", mc?"MC":"", chName[ich], ily, icen, ipad),
                        Form("Clusters[%c] :: #Deltay(x,#Phi) Ly[%d] Cen[%s] Pads[%s]", chSgn[ich], ily, cenName[icen], ipad?(ipad<kNpads?Form("%d", ipad+1):Form(">%d", kNpads)):"deconv"),
                        kPrez, kPt, kYrez, aa);
          php.AddLast(&hp[ih++]); np[isel]++;
        }
      }
    }
  }
  AliInfo(Form("Build %3d 3D projections.", ih));

  AliTRDresolutionProjection *pr0(NULL), *pr1(NULL);
  Int_t ly(0), ch(0), rcBin(as?as->FindBin(0.):-1), chBin(apt?apt->FindBin(0.):-1), ioff(0), cen(0), npad(0);
  for (Long64_t ib(0); ib < H->GetNbins(); ib++) {
    v = H->GetBinContent(ib, coord); if(v<1.) continue;
    ly = coord[kBC]-1;
    // RC selection
    if(rcBin>0 && coord[kSpeciesChgRC] == rcBin) continue;

    // charge selection
    ch = 0; // [-] track
    if(chBin>0 && coord[kPt] > chBin) ch = 1;  // [+] track
    cen = 0; // centrality selection
    if(acen) cen = coord[kNdim]-1;
    npad = 0;             // no. of pads selection
    if(anp) npad = TMath::Min(coord[kNdim+1]-1, Int_t(kNpads));

    if(debug){
      isel = ly*nCh*nCen*nNpad
            +ch*nCen*nNpad
            +cen*nNpad
            +npad;
      ioff=isel*4;
    } else {
      isel = ly; ioff = isel;
    }
    if(ioff>=ih){
      AliError(Form("Wrong selection %d [%3d]", ioff, ih));
      return kFALSE;
    }
    if(!(pr0=(AliTRDresolutionProjection*)php.At(ioff))) {
      AliError(Form("Missing projection %d", ioff));
      return kFALSE;
    }
    if(strcmp(pr0->fH->GetName(), Form("H%sClY%c%d%d%d", mc?"MC":"", chName[ch], ly, cen, npad))!=0){
      AliError(Form("Projection mismatch :: request[H%sClY%c%d%d%d] found[%s]", mc?"MC":"", chName[ch], ly, cen, npad, pr0->fH->GetName()));
      return kFALSE;
    }
    for(Int_t jh(0); jh<np[isel]; jh++) ((AliTRDresolutionProjection*)php.At(ioff+jh))->Increment(coord, v);
  }
  TObjArray *arr(NULL);
  fProj->AddAt(arr = new TObjArray(kClNproj), cidx);

  TH2 *h2(NULL);  Int_t jh(0);
  for(; ih--; ){
    if(!hp[ih].fH) continue;
    if(strchr(hp[ih].fH->GetName(), 'Q')){
      if(!(h2 = hp[ih].Projection2D(kNstat, kNcontours, 0, kFALSE))) continue;
      arr->AddAt(h2, jh++);
      if(!(h2 = (TH2*)gDirectory->Get(Form("%sEn", hp[ih].fH->GetName())))) continue;
      arr->AddAt(h2, jh++);
    } else {
      if((h2 = hp[ih].Projection2D(kNstat, kNcontours, 1))) arr->AddAt(h2, jh++);
    }
  }
  for(Int_t ily(0); ily<AliTRDgeometry::kNlayer; ily++){
    for(Int_t ich(0); ich<nCh; ich++){
      for(Int_t icen(0); icen<nCen; icen++){
        /*!dy*/
        if((pr0 = (AliTRDresolutionProjection*)php.FindObject(Form("H%sClY%c%d%d%d", mc?"MC":"", chName[ich], ily, icen, 0)))){
          for(Int_t ipad(1); ipad<nNpad; ipad++){
            if(!(pr1 = (AliTRDresolutionProjection*)php.FindObject(Form("H%sClY%c%d%d%d", mc?"MC":"", chName[1], ily, icen, ipad)))) continue;
            (*pr0)+=(*pr1);
          }
          pr0->fH->SetNameTitle(Form("H%sClY%c%d%d", mc?"MC":"", chName[ich], ily, icen), Form("Clusters[%c] :: #Deltay Ly[%d] Cen[%s]", chSgn[ich], ily, cenName[icen]));
          if((h2 = pr0->Projection2D(kNstat, kNcontours, 1))) arr->AddAt(h2, jh++);
          if(icen && (pr1 = (AliTRDresolutionProjection*)php.FindObject(Form("H%sClY%c%d%d%d", mc?"MC":"", chName[ich], ily, 0, 0)))) (*pr1)+=(*pr0);
        }
        /*!Q*/
        if((pr0 = (AliTRDresolutionProjection*)php.FindObject(Form("H%sClQ%c%d%d%d", mc?"MC":"", chName[ich], ily, icen, 0)))){
          for(Int_t ipad(1); ipad<nNpad; ipad++){
            if(!(pr1 = (AliTRDresolutionProjection*)php.FindObject(Form("H%sClQ%c%d%d%d", mc?"MC":"", chName[1], ily, icen, ipad)))) continue;
            (*pr0)+=(*pr1);
          }
          pr0->fH->SetNameTitle(Form("H%sClQ%c%d%d", mc?"MC":"", chName[ich], ily, icen), Form("Clusters[%c] :: Q Ly[%d] Cen[%s]", chSgn[ich], ily, cenName[icen]));
          if((h2 = pr0->Projection2D(kNstat, kNcontours, 2, kFALSE))) arr->AddAt(h2, jh++);
          if((h2 = (TH2*)gDirectory->Get(Form("%sEn", pr0->fH->GetName())))) arr->AddAt(h2, jh++);
          pr0->fH->SetName(Form("H%sClQS%c%d%d", mc?"MC":"", chName[ich], ily, icen));
          if((h2 = pr0->Projection2D(kNstat, kNcontours, 0))) arr->AddAt(h2, jh++);
          if(icen && (pr1 = (AliTRDresolutionProjection*)php.FindObject(Form("H%sClQ%c%d%d%d", mc?"MC":"", chName[ich], ily, 0, 0)))) (*pr1)+=(*pr0);
        }
        /*!YXTC*/
        if((pr0 = (AliTRDresolutionProjection*)php.FindObject(Form("H%sClYXTC%c%d%d%d", mc?"MC":"", chName[ich], ily, icen, 0)))){
          for(Int_t ipad(1); ipad<nNpad; ipad++){
            if(!(pr1 = (AliTRDresolutionProjection*)php.FindObject(Form("H%sClYXTC%c%d%d%d", mc?"MC":"", chName[1], ily, icen, ipad)))) continue;
            (*pr0)+=(*pr1);
          }
          pr0->fH->SetNameTitle(Form("H%sClYXTC%c%d%d", mc?"MC":"", chName[ich], ily, icen), Form("Clusters[%c] :: #Deltay(x,TC) Ly[%d] Cen[%s]", chSgn[ich], ily, cenName[icen]));
          if((h2 = pr0->Projection2D(kNstat, kNcontours, 1))) arr->AddAt(h2, jh++);
          if(icen && (pr1 = (AliTRDresolutionProjection*)php.FindObject(Form("H%sClYXTC%c%d%d%d", mc?"MC":"", chName[ich], ily, 0, 0)))) (*pr1)+=(*pr0);
        }
        /*!YXPh*/
        if((pr0 = (AliTRDresolutionProjection*)php.FindObject(Form("H%sClYXPh%c%d%d%d", mc?"MC":"", chName[ich], ily, icen, 0)))){
          for(Int_t ipad(1); ipad<nNpad; ipad++){
            if(!(pr1 = (AliTRDresolutionProjection*)php.FindObject(Form("H%sClYXPh%c%d%d%d", mc?"MC":"", chName[1], ily, icen, ipad)))) continue;
            (*pr0)+=(*pr1);
          }
          pr0->fH->SetNameTitle(Form("H%sClYXPh%c%d%d", mc?"MC":"", chName[ich], ily, icen), Form("Clusters[%c] :: #Deltay(x,#Phi) Ly[%d] Cen[%s]", chSgn[ich], ily, cenName[icen]));
          if((h2 = pr0->Projection2D(kNstat, kNcontours, 1))) arr->AddAt(h2, jh++);
          if(icen && (pr1 = (AliTRDresolutionProjection*)php.FindObject(Form("H%sClYXPh%c%d%d%d", mc?"MC":"", chName[ich], ily, 0, 0)))) (*pr1)+=(*pr0);
        }
      } // end centrality integration
      /*!dy*/
      if((pr0 = (AliTRDresolutionProjection*)php.FindObject(Form("H%sClY%c%d%d%d", mc?"MC":"", chName[ich], ily, 0, 0)))){
        pr0->fH->SetNameTitle(Form("H%sClY%c%d", mc?"MC":"", chName[ich], ily),
                              Form("Clusters[%c]:: #Deltay Ly[%d]", chSgn[ich], ily));
        if((h2 = pr0->Projection2D(kNstat, kNcontours, 1))) arr->AddAt(h2, jh++);
        if(ich && (pr1 = (AliTRDresolutionProjection*)php.FindObject(Form("H%sClY%c%d%d%d", mc?"MC":"", chName[0], ily, 0, 0)))) (*pr1)+=(*pr0);
      }
      /*!Q*/
      if((pr0 = (AliTRDresolutionProjection*)php.FindObject(Form("H%sClQ%c%d%d%d", mc?"MC":"", chName[ich], ily, 0, 0)))){
        pr0->fH->SetNameTitle(Form("H%sClQ%c%d", mc?"MC":"", chName[ich], ily),
                              Form("Clusters[%c]:: Q Ly[%d]", chSgn[ich], ily));
        if((h2 = pr0->Projection2D(kNstat, kNcontours, 2))) arr->AddAt(h2, jh++);
        pr0->fH->SetName(Form("H%sClQS%c%d", mc?"MC":"", chName[ich], ily));
        if((h2 = pr0->Projection2D(kNstat, kNcontours, 0))) arr->AddAt(h2, jh++);
        if(ich && (pr1 = (AliTRDresolutionProjection*)php.FindObject(Form("H%sClQ%c%d%d%d", mc?"MC":"", chName[0], ily, 0, 0)))) (*pr1)+=(*pr0);
      }
      /*!YXTC*/
      if((pr0 = (AliTRDresolutionProjection*)php.FindObject(Form("H%sClYXTC%c%d%d%d", mc?"MC":"", chName[ich], ily, 0, 0)))){
        pr0->fH->SetNameTitle(Form("H%sClYXTC%c%d", mc?"MC":"", chName[ich], ily),
                              Form("Clusters[%c]:: #Deltay(x,TC) Ly[%d]", chSgn[ich], ily));
        if((h2 = pr0->Projection2D(kNstat, kNcontours, 1))) arr->AddAt(h2, jh++);
        if(ich && (pr1 = (AliTRDresolutionProjection*)php.FindObject(Form("H%sClYXTC%c%d%d%d", mc?"MC":"", chName[0], ily, 0, 0)))) (*pr1)+=(*pr0);
      }
      /*!YXPh*/
      if((pr0 = (AliTRDresolutionProjection*)php.FindObject(Form("H%sClYXPh%c%d%d%d", mc?"MC":"", chName[ich], ily, 0, 0)))){
        pr0->fH->SetNameTitle(Form("H%sClYXPh%c%d", mc?"MC":"", chName[ich], ily),
                              Form("Clusters[%c]:: #Deltay(x,#Phi) Ly[%d]", chSgn[ich], ily));
        if((h2 = pr0->Projection2D(kNstat, kNcontours, 1))) arr->AddAt(h2, jh++);
        if(ich && (pr1 = (AliTRDresolutionProjection*)php.FindObject(Form("H%sClYXPh%c%d%d%d", mc?"MC":"", chName[0], ily, 0, 0)))) (*pr1)+=(*pr0);
      }
    } // end charge integration
    /*!dy*/
    if((pr0 = (AliTRDresolutionProjection*)php.FindObject(Form("H%sClY%c%d%d%d", mc?"MC":"", chName[0], ily, 0, 0)))){
      pr0->fH->SetNameTitle(Form("H%sClY%d", mc?"MC":"", ily), Form("Clusters :: #Deltay Ly[%d]", ily));
      if((h2 = pr0->Projection2D(kNstat, kNcontours, 1))) arr->AddAt(h2, jh++);
    }
    /*!YXPh*/
    if((pr0 = (AliTRDresolutionProjection*)php.FindObject(Form("H%sClYXPh%c%d%d%d", mc?"MC":"", chName[0], ily, 0, 0)))){
      pr0->fH->SetNameTitle(Form("H%sClYXPh%d", mc?"MC":"", ily), Form("Clusters :: #Deltay Ly[%d]", ily));
      if((h2 = pr0->Projection2D(kNstat, kNcontours, 1))) arr->AddAt(h2, jh++);
    }

  }
  AliInfo(Form("Done %3d 2D projections.", jh));
  return kTRUE;
}

//________________________________________________________
Bool_t AliTRDresolution::MakeProjectionTracklet(Bool_t mc)
{
// Analyse tracklet
  const Int_t kNcontours(9);
  const Int_t kNstat(30);
  const Int_t kNstatQ(30);
  Int_t cidx=mc?kMCtracklet:kTracklet;
  if(fProj && fProj->At(cidx)) return kTRUE;
  if(!fContainer){
    AliError("Missing data container.");
    return kFALSE;
  }
  const Char_t *projName[] = {"hTracklet2Track", "hTracklet2MC"};
  THnSparse *H(NULL);
  if(!(H = (THnSparse*)fContainer->FindObject(projName[Int_t(mc)]))){
    AliError(Form("Missing/Wrong data @ %s.", projName[Int_t(mc)]));
    return kFALSE;
  }
  const Int_t mdim(kNdim+8);
  Int_t ndim(H->GetNdimensions()); Bool_t debug(ndim>Int_t(kNdimTrklt));
  Int_t coord[mdim]; memset(coord, 0, sizeof(Int_t) * mdim); Double_t v = 0.;
  TAxis *aa[mdim], *as(NULL), *ap(NULL), *ac(NULL); memset(aa, 0, sizeof(TAxis*) * mdim);
  for(Int_t id(0); id<ndim; id++) aa[id] = H->GetAxis(id);
  if(ndim > Int_t(kSpeciesChgRC)) as = H->GetAxis(kSpeciesChgRC); // init species/charge selection
  if(ndim > Int_t(kPt))           ap = H->GetAxis(kPt);           // init pt selection
  if(ndim > Int_t(kNdim))         ac = H->GetAxis(kNdim);         // init centrality selection
  // calculate size depending on debug level
  const Int_t nCen(debug?Int_t(AliTRDeventInfo::kCentralityClasses):1);
  const Int_t nPt(debug?Int_t(kNpt):1);
  const Int_t nSpc(1);//ndim>kNdimTrklt?fgkNbins[kSpeciesChgRC]:1);
  const Int_t nCh(debug?Int_t(kNcharge):1);

  // build list of projections
  const Int_t nsel(AliTRDeventInfo::kCentralityClasses*AliTRDgeometry::kNlayer*kNpt*(AliPID::kSPECIES*kNcharge + 1));
  // define rebinning strategy
  const Int_t nEtaPhi(4); Int_t rebinEtaPhiX[nEtaPhi] = {1, 2, 5, 1}, rebinEtaPhiY[nEtaPhi] = {2, 1, 1, 5};
  AliTRDresolutionProjection hp[kTrkltNproj]; TObjArray php(kTrkltNproj);
  Int_t ih(0), isel(-1), np[nsel]; memset(np, 0, nsel*sizeof(Int_t));
  const Char_t chName[kNcharge] = {'n', 'p'};const Char_t chSgn[kNcharge] = {'-', '+'};
  const Char_t ptName[kNpt] = {'l', 'i', 'h'};
  const Char_t *ptCut[kNpt] = {"p_{t}[GeV/c]<0.8", "0.8<=p_{t}[GeV/c]<1.5", "p_{t}[GeV/c]>=1.5"};
//  const Char_t *cenName[AliTRDeventInfo::kCentralityClasses] = {"0-10%", "10-20%", "20-50%", "50-80%", "80-100%"};
  const Char_t *cenName[AliTRDeventInfo::kCentralityClasses] = {"2800-inf", "2100-2799", "1400-2099", "700-1399", "0-699"};
  for(Int_t icen(0); icen<nCen; icen++){
    for(Int_t ily(0); ily<AliTRDgeometry::kNlayer; ily++){
      for(Int_t ipt(0); ipt<nPt; ipt++){
        for(Int_t isp(0); isp<nSpc; isp++){
          for(Int_t ich(0); ich<nCh; ich++){
            isel++; // new selection
            hp[ih].Build(Form("H%sTrkltY%c%c%d%d%d", mc?"MC":"", chName[ich], ptName[ipt], isp, ily, icen),
                        Form("Tracklets[%s%c]:: #Deltay{%s} Ly[%d] Cen[%s]", AliPID::ParticleLatexName(isp), chSgn[ich], ptCut[ipt], ily, cenName[icen]),
                        kEta, kPhi, kYrez, aa);
            //hp[ih].SetShowRange(-0.1,0.1);
            hp[ih].SetRebinStrategy(nEtaPhi, rebinEtaPhiX, rebinEtaPhiY);
            php.AddLast(&hp[ih++]); np[isel]++;
            hp[ih].Build(Form("H%sTrkltPh%c%c%d%d%d", mc?"MC":"", chName[ich], ptName[ipt], isp, ily, icen),
                        Form("Tracklets[%s%c]:: #Delta#phi{%s} Ly[%d] Cen[%s]", AliPID::ParticleLatexName(isp), chSgn[ich], ptCut[ipt], ily, cenName[icen]),
                        kEta, kPhi, kPrez, aa);
            //hp[ih].SetShowRange(-0.5,0.5);
            hp[ih].SetRebinStrategy(nEtaPhi, rebinEtaPhiX, rebinEtaPhiY);
            php.AddLast(&hp[ih++]); np[isel]++;
            hp[ih].Build(Form("H%sTrkltQ%c%c%d%d%d", mc?"MC":"", chName[ich], ptName[ipt], isp, ily, icen),
                        Form("Tracklets[%s%c]:: dQdl{%s} Ly[%d] Cen[%s]", AliPID::ParticleLatexName(isp), chSgn[ich], ptCut[ipt], ily, cenName[icen]),
                        kEta, kPhi, kZrez, aa);
            hp[ih].SetShowRange(1.,2.3);
            hp[ih].SetRebinStrategy(nEtaPhi, rebinEtaPhiX, rebinEtaPhiY);
            php.AddLast(&hp[ih++]); np[isel]++;
            hp[ih].Build(Form("H%sTrkltTB%c%c%d%d%d", mc?"MC":"", chName[ich], ptName[ipt], isp, ily, icen),
                        Form("Tracklets[%s%c]:: OccupancyTB{%s} Ly[%d] Cen[%s]", AliPID::ParticleLatexName(isp), chSgn[ich], ptCut[ipt], ily, cenName[icen]),
                        kEta, kPhi, kNdim+1, aa);
            hp[ih].SetShowRange(30., 70.);
            hp[ih].SetRebinStrategy(nEtaPhi, rebinEtaPhiX, rebinEtaPhiY);
            php.AddLast(&hp[ih++]); np[isel]++;
          }
        }
        if(ndim==kNdimTrklt) continue;

        isel++; // new selection
        hp[ih].Build(Form("H%sTrkltRCZ%c%d%d", mc?"MC":"", ptName[ipt], ily, icen),
                    Form("Tracklets[RC]:: #Deltaz{%s} Ly[%d] Cen[%s]", ptCut[ipt], ily, cenName[icen]),
                    kEta, kPhi, kZrez, aa);
  //      hp[ih].SetShowRange(-0.1,0.1);
        hp[ih].SetRebinStrategy(nEtaPhi, rebinEtaPhiX, rebinEtaPhiY);
        php.AddLast(&hp[ih++]); np[isel]++;
/*        hp[ih].Build(Form("H%sTrkltRCY%c%d%d", mc?"MC":"", ptName[ipt], ily, icen),
                    Form("Tracklets[RC]:: #Deltay{%s} Ly[%d] Cen[%s]", ptCut[ipt], ily, cenName[icen]),
                    kEta, kPhi, kYrez, aa);
        //hp[ih].SetShowRange(-0.1,0.1);
        hp[ih].SetRebinStrategy(nEtaPhi, rebinEtaPhiX, rebinEtaPhiY);
        php.AddLast(&hp[ih++]); np[isel]++;
        hp[ih].Build(Form("H%sTrkltRCPh%c%d%d", mc?"MC":"", ptName[ipt], ily, icen),
                    Form("Tracklets[RC]:: #Delta#phi{%s} Ly[%d] Cen[%s]", ptCut[ipt], ily, cenName[icen]),
                    kEta, kPhi, kPrez, aa);
        //hp[ih].SetShowRange(-0.1,0.1);
        hp[ih].SetRebinStrategy(nEtaPhi, rebinEtaPhiX, rebinEtaPhiY);
        php.AddLast(&hp[ih++]); np[isel]++;
        hp[ih].Build(Form("H%sTrkltRCQ%c%d%d", mc?"MC":"", ptName[ipt], ily, icen),
                    Form("Tracklets[RC]:: dQdl{%s} Ly[%d] Cen[%s]", ptCut[ipt], ily, cenName[icen]),
                    kEta, kPhi, kNdim, aa);
        //hp[ih].SetShowRange(-0.1,0.1);
        hp[ih].SetRebinStrategy(nEtaPhi, rebinEtaPhiX, rebinEtaPhiY);
        php.AddLast(&hp[ih++]); np[isel]++;*/
      }
    }
  }
  AliInfo(Form("Build %3d 3D projections.", ih));

  AliTRDresolutionProjection *pr0(NULL), *pr1(NULL);
  Int_t ly(0), ch(0), sp(2), rcBin(as?as->FindBin(0.):-1), pt(0), cen(0), ioff(0), jspc(nSpc*nCh+1), kspc(nSpc*nCh*4+1);
  for (Long64_t ib(0); ib < H->GetNbins(); ib++) {
    v = H->GetBinContent(ib, coord);
    if(v<1.) continue;
    ly = coord[kBC]-1; // layer selection
    // charge selection
    ch = 0; sp=0;// [e-] track [dafault]
    if(rcBin>0){ // debug mode in which species/charge are also saved
      sp = Int_t(TMath::Abs(as->GetBinCenter(coord[kSpeciesChgRC])))-1;
      if(coord[kSpeciesChgRC] > rcBin) ch = 1;  // [+] track
      else if(coord[kSpeciesChgRC] == rcBin) ch = 2;  // [RC] track
    }
    // pt selection
    pt = 0; // low pt
    if(ap) pt = TMath::Min(coord[kPt]-1, Int_t(kNpt)-1);
    // centrality selection
    cen = 0; // default
    if(ac) cen = coord[kNdim]-1;
    // global selection
    if(ndim==kNdimTrklt){
      ioff = ly*4;
      isel = ly;
    } else {
      isel = cen*AliTRDgeometry::kNlayer*nPt*jspc+ly*nPt*jspc+pt*jspc; isel+=sp<0?(nSpc*nCh):ch;
      ioff = cen*AliTRDgeometry::kNlayer*nPt*kspc+ly*nPt*kspc+pt*kspc; ioff+=sp<0?((nSpc*nCh)*4):(ch*4);
    }
    if(ioff>=ih){
      AliError(Form("Wrong selection %d [%3d]", ioff, ih));
      return kFALSE;
    }
    if(!(pr0=(AliTRDresolutionProjection*)php.At(ioff))) {
      AliError(Form("Missing projection %d", ioff));
      return kFALSE;
    }
    if(sp>=0){
      if(strcmp(pr0->fH->GetName(), Form("H%sTrkltY%c%c%d%d%d", mc?"MC":"", chName[ch], ptName[pt], sp, ly, cen))!=0){
        AliError(Form("Projection mismatch :: request[H%sTrkltY%c%c%d%d%d] found[%s]", mc?"MC":"", chName[ch], ptName[pt], sp, ly, cen, pr0->fH->GetName()));
        return kFALSE;
      }
    } else {
      if(strcmp(pr0->fH->GetName(), Form("H%sTrkltRCZ%c%d%d", mc?"MC":"", ptName[pt], ly, cen))!=0){
        AliError(Form("Projection mismatch :: request[H%sTrkltRCZ%c%d%d] found[%s]", mc?"MC":"", ptName[pt], ly, cen, pr0->fH->GetName()));
        return kFALSE;
      }
    }
    for(Int_t jh(0); jh<np[isel]; jh++) ((AliTRDresolutionProjection*)php.At(ioff+jh))->Increment(coord, v);
  }
  TObjArray *arr(NULL);
  fProj->AddAt(arr = new TObjArray(kTrkltNproj), cidx);

  TH2 *h2(NULL); Int_t jh(0);
  for(; ih--; ){
    if(!hp[ih].fH) continue;
    Int_t mid(0), nstat(kNstat);
    if(strchr(hp[ih].fH->GetName(), 'Q')){ mid=2; nstat=kNstatQ;}
    if(!(h2 = hp[ih].Projection2D(nstat, kNcontours, mid))) continue;
    arr->AddAt(h2, jh++);
  }
  // build combined performance plots
  for(Int_t ily(0); ily<AliTRDgeometry::kNlayer; ily++){
    for(Int_t ich(0); ich<nCh; ich++){
      for(Int_t icen(0); icen<nCen; icen++){
        for(Int_t ipt(0); ipt<nPt; ipt++){
          /*!dy*/
          if((pr0 = (AliTRDresolutionProjection*)php.FindObject(Form("H%sTrkltY%c%c%d%d%d", mc?"MC":"", chName[ich], ptName[ipt], 0, ily, icen)))){
            for(Int_t isp(1); isp<nSpc; isp++){
              if(!(pr1 = (AliTRDresolutionProjection*)php.FindObject(Form("H%sTrkltY%c%c%d%d%d", mc?"MC":"", chName[ich], ptName[ipt], isp, ily, icen)))) continue;
              (*pr0)+=(*pr1);
            }
            pr0->fH->SetNameTitle(Form("H%sTrkltY%c%c%d%d", mc?"MC":"", chName[ich], ptName[ipt], ily, icen),
                                      Form("Tracklets[%c]:: #Deltay{%s} Ly[%d] Cen[%s]", chSgn[ich], ptCut[ipt], ily, cenName[icen]));
            if((h2 = pr0->Projection2D(kNstat, kNcontours, 1))) arr->AddAt(h2, jh++);
            if(ipt && (pr1 = (AliTRDresolutionProjection*)php.FindObject(Form("H%sTrkltY%c%c%d%d%d", mc?"MC":"", chName[ich], ptName[0], 0, ily, icen)))) (*pr1)+=(*pr0);
          }
          /*!dphi*/
          if((pr0 = (AliTRDresolutionProjection*)php.FindObject(Form("H%sTrkltPh%c%c%d%d%d", mc?"MC":"", chName[ich], ptName[ipt], 0, ily, icen)))){
            for(Int_t isp(1); isp<nSpc; isp++){
              if(!(pr1 = (AliTRDresolutionProjection*)php.FindObject(Form("H%sTrkltPh%c%c%d%d%d", mc?"MC":"", chName[ich], ptName[ipt], isp, ily, icen)))) continue;
              (*pr0)+=(*pr1);
            }
            pr0->fH->SetNameTitle(Form("H%sTrkltPh%c%c%d%d", mc?"MC":"", chName[ich], ptName[ipt], ily, icen),
                                      Form("Tracklets[%c]:: #Delta#phi{%s} Ly[%d] Cen[%s]", chSgn[ich], ptCut[ipt], ily, cenName[icen]));
            if((h2 = pr0->Projection2D(kNstat, kNcontours, 1))) arr->AddAt(h2, jh++);
            if(ipt && (pr1 = (AliTRDresolutionProjection*)php.FindObject(Form("H%sTrkltPh%c%c%d%d%d", mc?"MC":"", chName[ich], ptName[0], 0, ily, icen)))) (*pr1)+=(*pr0);
          }
          /*!dQ/dl*/
          if((pr0 = (AliTRDresolutionProjection*)php.FindObject(Form("H%sTrkltQ%c%c%d%d%d", mc?"MC":"", chName[ich], ptName[ipt], 0, ily, icen)))){
            for(Int_t isp(1); isp<nSpc; isp++){
              if(!(pr1 = (AliTRDresolutionProjection*)php.FindObject(Form("H%sTrkltQ%c%c%d%d%d", mc?"MC":"", chName[ich], ptName[ipt], isp, ily, icen)))) continue;
              (*pr0)+=(*pr1);
            }
            pr0->fH->SetNameTitle(Form("H%sTrkltQ%c%c%d%d", mc?"MC":"", chName[ich], ptName[ipt], ily, icen),
                                      Form("Tracklets[%c]:: dQdl{%s} Ly[%d] Cen[%s]", chSgn[ich], ptCut[ipt], ily, cenName[icen]));
            if((h2 = pr0->Projection2D(kNstatQ, kNcontours, 2))) arr->AddAt(h2, jh++);
            pr0->fH->SetNameTitle(Form("H%sTrkltQS%c%c%d%d", mc?"MC":"", chName[ich], ptName[ipt], ily, icen),
                                      Form("Tracklets[%c]:: dQdl{%s} Ly[%d] Cen[%s]", chSgn[ich], ptCut[ipt], ily, cenName[icen]));
            pr0->SetShowRange(2.4, 5.1);
            if((h2 = pr0->Projection2D(kNstat, kNcontours, 0))) arr->AddAt(h2, jh++);
            if(ipt && (pr1 = (AliTRDresolutionProjection*)php.FindObject(Form("H%sTrkltQ%c%c%d%d%d", mc?"MC":"", chName[ich], ptName[0], 0, ily, icen)))) (*pr1)+=(*pr0);
          }
          /*!TB occupancy*/
          if((pr0 = (AliTRDresolutionProjection*)php.FindObject(Form("H%sTrkltTB%c%c%d%d%d", mc?"MC":"", chName[ich], ptName[ipt], 0, ily, icen)))){
            for(Int_t isp(1); isp<nSpc; isp++){
              if(!(pr1 = (AliTRDresolutionProjection*)php.FindObject(Form("H%sTrkltTB%c%c%d%d%d", mc?"MC":"", chName[ich], ptName[ipt], isp, ily, icen)))) continue;
              (*pr0)+=(*pr1);
            }
            pr0->fH->SetNameTitle(Form("H%sTrkltTB%c%c%d%d", mc?"MC":"", chName[ich], ptName[ipt], ily, icen),
                                      Form("Tracklets[%c]:: OccupancyTB{%s} Ly[%d] Cen[%s]", chSgn[ich], ptCut[ipt], ily, cenName[icen]));
            if((h2 = pr0->Projection2D(kNstat, kNcontours))) arr->AddAt(h2, jh++);
            if(ipt && (pr1 = (AliTRDresolutionProjection*)php.FindObject(Form("H%sTrkltTB%c%c%d%d%d", mc?"MC":"", chName[ich], ptName[0], 0, ily, icen)))) (*pr1)+=(*pr0);
          }
        } // end pt integration
        /*!dy*/
        if((pr0 = (AliTRDresolutionProjection*)php.FindObject(Form("H%sTrkltY%c%c%d%d%d", mc?"MC":"", chName[ich], ptName[0], 0, ily, icen)))){
          pr0->fH->SetNameTitle(Form("H%sTrkltY%c%d%d", mc?"MC":"", chName[ich], ily, icen),
                                Form("Tracklets[%c]:: #Deltay Ly[%d] Cen[%s]", chSgn[ich], ily, cenName[icen]));
          if((h2 = pr0->Projection2D(kNstat, kNcontours, 1))) arr->AddAt(h2, jh++);
          if(icen && (pr1 = (AliTRDresolutionProjection*)php.FindObject(Form("H%sTrkltY%c%c%d%d%d", mc?"MC":"", chName[ich], ptName[0], 0, ily, 0)))) (*pr1)+=(*pr0);
        }
        /*!dphi*/
        if((pr0 = (AliTRDresolutionProjection*)php.FindObject(Form("H%sTrkltPh%c%c%d%d%d", mc?"MC":"", chName[ich], ptName[0], 0, ily, icen)))){
          pr0->fH->SetNameTitle(Form("H%sTrkltPh%c%d%d", mc?"MC":"", chName[ich], ily, icen),
                                Form("Tracklets[%c]:: #Delta#phi Ly[%d] Cen[%s]", chSgn[ich], ily, cenName[icen]));
          if((h2 = pr0->Projection2D(kNstat, kNcontours, 1))) arr->AddAt(h2, jh++);
          if(icen && (pr1 = (AliTRDresolutionProjection*)php.FindObject(Form("H%sTrkltPh%c%c%d%d%d", mc?"MC":"", chName[ich], ptName[0], 0, ily, 0)))) (*pr1)+=(*pr0);
        }
        /*!dQ/dl*/
        if((pr0 = (AliTRDresolutionProjection*)php.FindObject(Form("H%sTrkltQ%c%c%d%d%d", mc?"MC":"", chName[ich], ptName[0], 0, ily, icen)))){
          pr0->fH->SetNameTitle(Form("H%sTrkltQ%c%d%d", mc?"MC":"", chName[ich], ily, icen),
                                Form("Tracklets[%c]:: dQdl Ly[%d] Cen[%s]", chSgn[ich], ily, cenName[icen]));
          pr0->SetShowRange(1.,2.3);
          if((h2 = pr0->Projection2D(kNstatQ, kNcontours, 2))) arr->AddAt(h2, jh++);
          pr0->fH->SetNameTitle(Form("H%sTrkltQS%c%d%d", mc?"MC":"", chName[ich], ily, icen),
                                Form("Tracklets[%c]:: dQdl Ly[%d] Cen[%s]", chSgn[ich], ily, cenName[icen]));
          pr0->SetShowRange(2.4,5.1);
          if((h2 = pr0->Projection2D(kNstat, kNcontours, 0))) arr->AddAt(h2, jh++);
          if(icen && (pr1 = (AliTRDresolutionProjection*)php.FindObject(Form("H%sTrkltQ%c%c%d%d%d", mc?"MC":"", chName[ich], ptName[0], 0, ily, 0)))) (*pr1)+=(*pr0);
        }
        /*!TB occupancy*/
        if((pr0 = (AliTRDresolutionProjection*)php.FindObject(Form("H%sTrkltTB%c%c%d%d%d", mc?"MC":"", chName[ich], ptName[0], 0, ily, icen)))){
          pr0->fH->SetNameTitle(Form("H%sTrkltTB%c%d%d", mc?"MC":"", chName[ich], ily, icen),
                                Form("Tracklets[%c]:: OccupancyTB Ly[%d] Cen[%s]", chSgn[ich], ily, cenName[icen]));
          if((h2 = pr0->Projection2D(kNstat, kNcontours))) arr->AddAt(h2, jh++);
          if(icen && (pr1 = (AliTRDresolutionProjection*)php.FindObject(Form("H%sTrkltTB%c%c%d%d%d", mc?"MC":"", chName[ich], ptName[0], 0, ily, 0)))) (*pr1)+=(*pr0);
        }
      } // end centrality integration
      /*!dy*/
      if((pr0 = (AliTRDresolutionProjection*)php.FindObject(Form("H%sTrkltY%c%c%d%d%d", mc?"MC":"", chName[ich], ptName[0], 0, ily, 0)))){
        pr0->fH->SetNameTitle(Form("H%sTrkltY%c%d", mc?"MC":"", chName[ich], ily), Form("Tracklets[%c] :: #Deltay Ly[%d]", chSgn[ich], ily));
        if((h2 = pr0->Projection2D(kNstat, kNcontours, 1))) arr->AddAt(h2, jh++);
        if(ich && (pr1 = (AliTRDresolutionProjection*)php.FindObject(Form("H%sTrkltY%c%c%d%d%d", mc?"MC":"", chName[0], ptName[0], 0, ily, 0)))) (*pr1)+=(*pr0);
      }
      /*!dphi*/
      if((pr0 = (AliTRDresolutionProjection*)php.FindObject(Form("H%sTrkltPh%c%c%d%d%d", mc?"MC":"", chName[ich], ptName[0], 0, ily, 0)))){
        pr0->fH->SetNameTitle(Form("H%sTrkltPh%c%d", mc?"MC":"", chName[ich], ily), Form("Tracklets[%c] :: #Delta#phi Ly[%d]", chSgn[ich], ily));
        pr0->SetShowRange(-.9,.9);
        if((h2 = pr0->Projection2D(kNstat, kNcontours, 1))) arr->AddAt(h2, jh++);
        if(ich && (pr1 = (AliTRDresolutionProjection*)php.FindObject(Form("H%sTrkltPh%c%c%d%d%d", mc?"MC":"", chName[0], ptName[0], 0, ily, 0)))) (*pr1)+=(*pr0);
      }
      /*!dQ/dl*/
      if((pr0 = (AliTRDresolutionProjection*)php.FindObject(Form("H%sTrkltQ%c%c%d%d%d", mc?"MC":"", chName[ich], ptName[0], 0, ily, 0)))){
        pr0->fH->SetNameTitle(Form("H%sTrkltQ%c%d", mc?"MC":"", chName[ich], ily), Form("Tracklets[%c] :: dQdl Ly[%d]", chSgn[ich], ily));
        pr0->SetShowRange(1.,2.3);
        if((h2 = pr0->Projection2D(kNstatQ, kNcontours, 2))) arr->AddAt(h2, jh++);
        pr0->fH->SetNameTitle(Form("H%sTrkltQS%c%d", mc?"MC":"", chName[ich], ily), Form("Tracklets[%c] :: dQdl Ly[%d]", chSgn[ich], ily));
        pr0->SetShowRange(2.4,5.1);
        if((h2 = pr0->Projection2D(kNstat, kNcontours, 0))) arr->AddAt(h2, jh++);
        if(ich && (pr1 = (AliTRDresolutionProjection*)php.FindObject(Form("H%sTrkltQ%c%c%d%d%d", mc?"MC":"", chName[0], ptName[0], 0, ily, 0)))) (*pr1)+=(*pr0);
      }
      /*!TB occupancy*/
      if((pr0 = (AliTRDresolutionProjection*)php.FindObject(Form("H%sTrkltTB%c%c%d%d%d", mc?"MC":"", chName[ich], ptName[0], 0, ily, 0)))){
        pr0->fH->SetNameTitle(Form("H%sTrkltTB%c%d", mc?"MC":"", chName[ich], ily), Form("Tracklets[%c] :: OccupancyTB Ly[%d]", chSgn[ich], ily));
        if((h2 = pr0->Projection2D(kNstat, kNcontours))) arr->AddAt(h2, jh++);
        if(ich && (pr1 = (AliTRDresolutionProjection*)php.FindObject(Form("H%sTrkltTB%c%c%d%d%d", mc?"MC":"", chName[0], ptName[0], 0, ily, 0)))) (*pr1)+=(*pr0);
      }
    } // end charge integration
    /*!dy*/
    if((pr0 = (AliTRDresolutionProjection*)php.FindObject(Form("H%sTrkltY%c%c%d%d%d", mc?"MC":"", chName[0], ptName[0], 0, ily, 0)))){
      pr0->fH->SetNameTitle(Form("H%sTrkltY%d", mc?"MC":"", ily), Form("Tracklets :: #Deltay Ly[%d]", ily));
      if((h2 = pr0->Projection2D(kNstat, kNcontours, 1))) arr->AddAt(h2, jh++);
    }
    /*!dphi*/
    if((pr0 = (AliTRDresolutionProjection*)php.FindObject(Form("H%sTrkltPh%c%c%d%d%d", mc?"MC":"", chName[0], ptName[0], 0, ily, 0)))){
      pr0->fH->SetNameTitle(Form("H%sTrkltPh%d", mc?"MC":"", ily), Form("Tracklets :: #Delta#phi Ly[%d]", ily));
      pr0->SetShowRange(-.45,.45);
      if((h2 = pr0->Projection2D(kNstat, kNcontours, 1))) arr->AddAt(h2, jh++);
    }
    /*!dQdl*/
    if((pr0 = (AliTRDresolutionProjection*)php.FindObject(Form("H%sTrkltQ%c%c%d%d%d", mc?"MC":"", chName[0], ptName[0], 0, ily, 0)))){
      pr0->fH->SetNameTitle(Form("H%sTrkltQ%d", mc?"MC":"", ily), Form("Tracklets :: dQdl Ly[%d]", ily));
      pr0->SetShowRange(1.,2.3);
      if((h2 = pr0->Projection2D(kNstat, kNcontours, 2))) arr->AddAt(h2, jh++);
      pr0->fH->SetName(Form("H%sTrkltQS%d", mc?"MC":"", ily));
      pr0->SetShowRange(2.4,5.1);
      if((h2 = pr0->Projection2D(kNstat, kNcontours, 0))) arr->AddAt(h2, jh++);
    }
    /*!TB occupancy*/
    if((pr0 = (AliTRDresolutionProjection*)php.FindObject(Form("H%sTrkltTB%c%c%d%d%d", mc?"MC":"", chName[0], ptName[0], 0, ily, 0)))){
      pr0->fH->SetNameTitle(Form("H%sTrkltTB%d", mc?"MC":"", ily), Form("Tracklets :: OccupancyTB Ly[%d]", ily));
      if((h2 = pr0->Projection2D(kNstat, kNcontours))) arr->AddAt(h2, jh++);
    }

    /*! Row Cross processing*/
    for(Int_t icen(0); icen<nCen; icen++){
      /*!RC dz*/
      if((pr0 = (AliTRDresolutionProjection*)php.FindObject(Form("H%sTrkltRCZ%c%d%d", mc?"MC":"", ptName[0], ily, icen)))){
        for(Int_t ipt(0); ipt<kNpt; ipt++){
          if(!(pr1 = (AliTRDresolutionProjection*)php.FindObject(Form("H%sTrkltRCZ%c%d%d", mc?"MC":"", ptName[ipt], ily, icen)))) continue;
          (*pr0)+=(*pr1);
        }
        pr0->fH->SetNameTitle(Form("H%sTrkltRCZ%d%d", mc?"MC":"", ily, icen), Form("Tracklets[RC]:: #Deltaz Ly[%d] Cen[%s]", ily, cenName[icen]));
        if((h2 = pr0->Projection2D(kNstat, kNcontours, 1))) arr->AddAt(h2, jh++);
        if(icen && (pr1 = (AliTRDresolutionProjection*)php.FindObject(Form("H%sTrkltRCZ%c%d%d", mc?"MC":"", ptName[0], ily, 0)))) (*pr1)+=(*pr0);
      }
    } // end centrality integration for row cross
    /*!RC dz*/
    if((pr0 = (AliTRDresolutionProjection*)php.FindObject(Form("H%sTrkltRCZ%c%d%d", mc?"MC":"", ptName[0], ily, 0)))){
      pr0->fH->SetNameTitle(Form("H%sTrkltRCZ%d", mc?"MC":"", ily), Form("Tracklets[RC] :: #Deltaz Ly[%d]", ily));
      if((h2 = pr0->Projection2D(kNstat, kNcontours, 1))) arr->AddAt(h2, jh++);
    }
  } // end layer loop
  AliInfo(Form("Done %3d 2D projections.", jh));
  return kTRUE;
}

//________________________________________________________
Bool_t AliTRDresolution::MakeProjectionTrackIn(Bool_t mc)
{
// Analyse track in
  const Int_t kNcontours(9);
  const Int_t kNstat(30);
  Int_t cidx=mc?kMCtrackIn:kTrackIn;
  if(fProj && fProj->At(cidx)) return kTRUE;
  if(!fContainer){
    AliError("Missing data container.");
    return kFALSE;
  }
  const Char_t *projName[] = {"hTracklet2TRDin", "hTRDin2MC"};
  THnSparse *H(NULL);
  if(!(H = (THnSparse*)fContainer->FindObject(projName[Int_t(mc)]))){
    AliError(Form("Missing/Wrong data @ %s.", projName[Int_t(mc)]));
    return kFALSE;
  }

  const Int_t mdim(kNdim+3);
  Int_t coord[mdim]; memset(coord, 0, sizeof(Int_t) * mdim); Double_t v = 0.;
  Int_t ndim(H->GetNdimensions());
  TAxis *aa[mdim], *as(NULL), *ap(NULL), *ax(NULL), *abf(NULL); memset(aa, 0, sizeof(TAxis*) * (mdim));
  for(Int_t id(0); id<ndim; id++) aa[id] = H->GetAxis(id);
  if(ndim > Int_t(kSpeciesChgRC)) as = H->GetAxis(kSpeciesChgRC);
  if(ndim > Int_t(kPt))           ap = H->GetAxis(kPt);
  if(ndim > Int_t(kNdim)+1)       ax = H->GetAxis(kNdim+1);
  if(ndim > Int_t(kNdim)+2)      abf = H->GetAxis(kNdim+2);
  //AliInfo(Form("Using : Species[%c] Pt[%c] BunchFill[%c]", as?'y':'n', ap?'y':'n', abf?'y':'n'));
  const Int_t nPt(ndim>Int_t(kNdimTrkIn)?Int_t(kNpt):1);

  // build list of projections
  const Int_t nsel(kNpt*(AliPID::kSPECIES*kNcharge + 1));
  const Char_t chName[kNcharge] = {'n', 'p'};const Char_t chSgn[kNcharge] = {'-', '+'};
  const Char_t ptName[kNpt] = {'l', 'i', 'h'};
  const Char_t *ptCut[kNpt] = {"p_{t}[GeV/c]<0.8", "0.8<=p_{t}[GeV/c]<1.5", "p_{t}[GeV/c]>=1.5"};
  // define rebinning strategy
  const Int_t nEtaPhi(4); Int_t rebinEtaPhiX[nEtaPhi] = {1, 2, 5, 1}, rebinEtaPhiY[nEtaPhi] = {2, 1, 1, 5};
  AliTRDresolutionProjection hp[kMCTrkInNproj]; TObjArray php(kMCTrkInNproj+2);
  Int_t ih(0), isel(-1), np[nsel]; memset(np, 0, nsel*sizeof(Int_t));
  // define list of projections
  for(Int_t ipt(0); ipt<nPt; ipt++){
    for(Int_t isp(0); isp<AliPID::kSPECIES; isp++){
      for(Int_t ich(0); ich<kNcharge; ich++){
        isel++; // new selection
        hp[ih].Build(Form("H%sTrkInY%c%c%d", mc?"MC":"", chName[ich], ptName[ipt], isp),
                     Form("TrackIn[%s%c]:: #Deltay{%s}", AliPID::ParticleLatexName(isp), chSgn[ich], ptCut[ipt]),
                     kEta, kPhi, kYrez, aa);
        hp[ih].SetRebinStrategy(nEtaPhi, rebinEtaPhiX, rebinEtaPhiY);
        php.AddLast(&hp[ih++]); np[isel]++;
        hp[ih].Build(Form("H%sTrkInPh%c%c%d", mc?"MC":"", chName[ich], ptName[ipt], isp),
                     Form("TrackIn[%s%c]:: #Delta#phi{%s}", AliPID::ParticleLatexName(isp), chSgn[ich], ptCut[ipt]),
                     kEta, kPhi, kPrez, aa);
        hp[ih].SetRebinStrategy(nEtaPhi, rebinEtaPhiX, rebinEtaPhiY);
        php.AddLast(&hp[ih++]); np[isel]++;
        hp[ih].Build(Form("H%sTrkInQ%c%c%d", mc?"MC":"", chName[ich], ptName[ipt], isp),
                     Form("TrackIn[%s%c]:: dQdl {%s}", AliPID::ParticleLatexName(isp), chSgn[ich], ptCut[ipt]),
                     kEta, kPhi, kZrez, aa);
        hp[ih].SetRebinStrategy(nEtaPhi, rebinEtaPhiX, rebinEtaPhiY);
        php.AddLast(&hp[ih++]); np[isel]++;
        if(!ax) continue;
        hp[ih].Build(Form("H%sTrkInX%c%c%d", mc?"MC":"", chName[ich], ptName[ipt], isp),
                     Form("TrackIn[%s%c]:: #Deltax{%s}", AliPID::ParticleLatexName(isp), chSgn[ich], ptCut[ipt]),
                     kEta, kPhi, kNdim+1, aa);
        hp[ih].SetRebinStrategy(nEtaPhi, rebinEtaPhiX, rebinEtaPhiY);
        php.AddLast(&hp[ih++]); np[isel]++;
      }
    }
    isel++; // RC tracks
    hp[ih].Build(Form("H%sTrkInRCZ%c", mc?"MC":"", ptName[ipt]),
                  Form("TrackIn[RC]:: #Deltaz{%s}", ptCut[ipt]),
                  kEta, kPhi, kZrez, aa);
    hp[ih].SetRebinStrategy(nEtaPhi, rebinEtaPhiX, rebinEtaPhiY);
    php.AddLast(&hp[ih++]); np[isel]++;
    hp[ih].Build(Form("H%sTrkInRCY%c", mc?"MC":"", ptName[ipt]),
                  Form("TrackIn[RC]:: #Deltay{%s}", ptCut[ipt]),
                  kEta, kPhi, kYrez, aa);
    hp[ih].SetRebinStrategy(nEtaPhi, rebinEtaPhiX, rebinEtaPhiY);
    php.AddLast(&hp[ih++]); np[isel]++;
    hp[ih].Build(Form("H%sTrkInRCPh%c", mc?"MC":"", ptName[ipt]),
                  Form("TrackIn[RC]:: #Delta#phi{%s}", ptCut[ipt]),
                  kEta, kPhi, kPrez, aa);
    hp[ih].SetRebinStrategy(nEtaPhi, rebinEtaPhiX, rebinEtaPhiY);
    php.AddLast(&hp[ih++]); np[isel]++;
    if(!ax) continue;
    hp[ih].Build(Form("H%sTrkInRCX%c", mc?"MC":"", ptName[ipt]),
                  Form("TrackIn[RC]:: #Deltax{%s}", ptCut[ipt]),
                  kEta, kPhi, kNdim+1, aa);
    hp[ih].SetRebinStrategy(nEtaPhi, rebinEtaPhiX, rebinEtaPhiY);
    php.AddLast(&hp[ih++]); np[isel]++;
  }
  AliInfo(Form("Build %3d 3D projections.", ih));

  // fill projections
  Int_t ch(0), pt(0), sp(1), rcBin(as?as->FindBin(0.):-1), ioff(0);
  AliTRDresolutionProjection *pr0(NULL), *pr1(NULL);
  for (Long64_t ib(0); ib < H->GetNbins(); ib++) {
    v = H->GetBinContent(ib, coord);
    if(v<1.) continue;
    if(fBCbinTOF>0 && coord[kBC]!=fBCbinTOF) continue; // TOF bunch cross cut
    if(fBCbinFill>0 && abf && coord[kNdim+2]!=fBCbinTOF) continue; // Fill bunch cut
    // charge selection
    ch = 0; sp=1;// [-] track
    if(rcBin>0){ // debug mode in which species are also saved
      sp = Int_t(TMath::Abs(as->GetBinCenter(coord[kSpeciesChgRC])))-1;
      if(sp>=AliPID::kSPECIES){
        AliDebug(2, Form("Wrong SpeciesIndex[%d]. Rescale", sp));
        sp = AliPID::kSPECIES-1;
      }
      if(coord[kSpeciesChgRC] > rcBin) ch = 1;  // [+] track
      else if(coord[kSpeciesChgRC] == rcBin) ch = 2;  // [RC] track
    }
    // pt selection
    pt = 0; // low pt
    if(ap) pt = TMath::Min(coord[kPt]-1, Int_t(kNpt)-1);
    // global selection
    isel = pt*(AliPID::kSPECIES*kNcharge+1); isel+=sp<0?(AliPID::kSPECIES*kNcharge):(sp*kNcharge+ch);
    ioff = isel*(ax?4:3);
    if(ioff>=ih){
      AliError(Form("Wrong selection %d [%3d]", ioff, ih));
      return kFALSE;
    }
    if(!(pr0=(AliTRDresolutionProjection*)php.At(ioff))) {
      AliError(Form("Missing projection %d", ioff));
      return kFALSE;
    }
    if(sp>=0){
      if(strcmp(pr0->fH->GetName(), Form("H%sTrkInY%c%c%d", mc?"MC":"", chName[ch], ptName[pt], sp))!=0){
        AliError(Form("Projection mismatch :: request[H%sTrkInY%c%c%d] found[%s]", mc?"MC":"", chName[ch], ptName[pt], sp, pr0->fH->GetName()));
        return kFALSE;
      }
    } else {
      if(strcmp(pr0->fH->GetName(), Form("H%sTrkInRCZ%c", mc?"MC":"", ptName[pt]))!=0){
        AliError(Form("Projection mismatch :: request[H%sTrkltRCZ%c] found[%s]", mc?"MC":"", ptName[pt], pr0->fH->GetName()));
        return kFALSE;
      }
    }
    for(Int_t jh(0); jh<np[isel]; jh++) ((AliTRDresolutionProjection*)php.At(ioff+jh))->Increment(coord, v);
  }
  TObjArray *arr(NULL);
  fProj->AddAt(arr = new TObjArray(mc?kMCTrkInNproj:kTrkInNproj), cidx);

  TH2 *h2(NULL); Int_t jh(0);
  for(; ih--; ){
    if(!hp[ih].fH) continue;
    if(!(h2 = hp[ih].Projection2D(kNstat, kNcontours))) continue;
    arr->AddAt(h2, jh++);
  }
  // build combined performance plots
  // combine up the tree of projections
  AliTRDresolutionProjection xlow[2], specY[kNcharge*AliPID::kSPECIES], specPh[kNcharge*AliPID::kSPECIES], specQ[kNcharge*AliPID::kSPECIES];
  for(Int_t ich(0); ich<kNcharge; ich++){
    // PID dependency - summation over pt
    for(Int_t isp(0); isp<AliPID::kSPECIES; isp++){
      /*!dy*/
      Int_t idx(ich*AliPID::kSPECIES+isp);
      if((pr0 = (AliTRDresolutionProjection*)php.FindObject(Form("H%sTrkInY%c%c%d", mc?"MC":"", chName[ich], ptName[0], isp)))){
        specY[idx] = (*pr0);
        specY[idx].SetNameTitle(Form("H%sTrkInY%c%d", mc?"MC":"", chName[ich], isp), "Sum over pt");
        specY[idx].fH->SetNameTitle(Form("H%sTrkInY%c%d", mc?"MC":"", chName[ich], isp),
                              Form("TrackIn[%s%c]:: #Deltay", AliPID::ParticleLatexName(isp), chSgn[ich]));
        for(Int_t ipt(1); ipt<nPt; ipt++){
          if(!(pr1 = (AliTRDresolutionProjection*)php.FindObject(Form("H%sTrkInY%c%c%d", mc?"MC":"", chName[ich], ptName[ipt], isp)))) continue;
          specY[idx]+=(*pr1);
        }
        php.AddLast(&specY[idx]);
        if((h2 = specY[idx].Projection2D(kNstat, kNcontours, 1, kFALSE))) arr->AddAt(h2, jh++);
        if((h2 = (TH2*)gDirectory->Get(Form("%sEn", specY[idx].fH->GetName())))) arr->AddAt(h2, jh++);
        if(ich && (pr1 = (AliTRDresolutionProjection*)php.FindObject(Form("H%sTrkInY%c%d", mc?"MC":"", chName[0], isp)))) (*pr1)+=specY[idx];
      }
      /*!dphi*/
      if((pr0 = (AliTRDresolutionProjection*)php.FindObject(Form("H%sTrkInPh%c%c%d", mc?"MC":"", chName[ich], ptName[0], isp)))){
        specPh[idx] = (*pr0);
        specPh[idx].SetNameTitle(Form("H%sTrkInPh%c%d", mc?"MC":"", chName[ich], isp), "Sum over pt");
        specPh[idx].fH->SetNameTitle(Form("H%sTrkInPh%c%d", mc?"MC":"", chName[ich], isp),
                              Form("TrackIn[%s%c]:: #Delta#phi", AliPID::ParticleLatexName(isp), chSgn[ich]));
        specPh[idx].SetShowRange(-1.5, 1.5);
        for(Int_t ipt(1); ipt<nPt; ipt++){
          if(!(pr1 = (AliTRDresolutionProjection*)php.FindObject(Form("H%sTrkInPh%c%c%d", mc?"MC":"", chName[ich], ptName[ipt], isp)))) continue;
          specPh[idx]+=(*pr1);
        }
        php.AddLast(&specPh[idx]);
        if((h2 = specPh[idx].Projection2D(kNstat, kNcontours, 1))) arr->AddAt(h2, jh++);
        if(ich && (pr1 = (AliTRDresolutionProjection*)php.FindObject(Form("H%sTrkInPh%c%d", mc?"MC":"", chName[0], isp)))) (*pr1)+=specPh[idx];
      }
      /*!dQdl*/
      if((pr0 = (AliTRDresolutionProjection*)php.FindObject(Form("H%sTrkInQ%c%c%d", mc?"MC":"", chName[ich], ptName[0], isp)))){
        specQ[idx] = (*pr0);
        specQ[idx].SetNameTitle(Form("H%sTrkInQ%c%d", mc?"MC":"", chName[ich], isp), "Sum over pt");
        specQ[idx].fH->SetNameTitle(Form("H%sTrkInQ%c%d", mc?"MC":"", chName[ich], isp),
                              Form("TrackIn[%s%c]:: dQdl", AliPID::ParticleLatexName(isp), chSgn[ich]));
        specQ[idx].SetShowRange(-2.2, -1.75);
        specQ[idx].fH->GetZaxis()->SetTitle("dQdl [a.u.]");
        for(Int_t ipt(1); ipt<nPt; ipt++){
          if(!(pr1 = (AliTRDresolutionProjection*)php.FindObject(Form("H%sTrkInQ%c%c%d", mc?"MC":"", chName[ich], ptName[ipt], isp)))) continue;
          specQ[idx]+=(*pr1);
        }
        php.AddLast(&specQ[idx]);
        if((h2 = specQ[idx].Projection2D(kNstat, kNcontours, 2))) arr->AddAt(h2, jh++);
        specQ[idx].fH->SetName(Form("H%sTrkInQS%c%d", mc?"MC":"", chName[ich], isp));
        specQ[idx].SetShowRange(-1.75, -1.25);
        if((h2 = specQ[idx].Projection2D(kNstat, kNcontours, 0))) arr->AddAt(h2, jh++);
        if(ich && (pr1 = (AliTRDresolutionProjection*)php.FindObject(Form("H%sTrkInQ%c%d", mc?"MC":"", chName[0], isp)))) (*pr1)+=specQ[idx];
      }
    } // end PID loop for pt integration

    // pt dependency - summation over PID
    for(Int_t ipt(0); ipt<nPt; ipt++){
      /*!dy*/
      if((pr0 = (AliTRDresolutionProjection*)php.FindObject(Form("H%sTrkInY%c%c%d", mc?"MC":"", chName[ich], ptName[ipt], 0)))){
        for(Int_t isp(1); isp<AliPID::kSPECIES; isp++){
          if(!(pr1 = (AliTRDresolutionProjection*)php.FindObject(Form("H%sTrkInY%c%c%d", mc?"MC":"", chName[ich], ptName[ipt], isp)))) continue;
          (*pr0)+=(*pr1);
        }
        pr0->fH->SetNameTitle(Form("H%sTrkInY%c%c", mc?"MC":"", chName[ich], ptName[ipt]),
                                  Form("TrackIn[%c]:: #Deltay{%s}", chSgn[ich], ptCut[ipt]));
        pr0->SetShowRange(-0.3, 0.3);
        if((h2 = pr0->Projection2D(kNstat, kNcontours, 1))) arr->AddAt(h2, jh++);
        if(ipt && (pr1 = (AliTRDresolutionProjection*)php.FindObject(Form("H%sTrkInY%c%c%d", mc?"MC":"", chName[ich], ptName[0], 0)))) (*pr1)+=(*pr0);
      }
      /*!dphi*/
      if((pr0 = (AliTRDresolutionProjection*)php.FindObject(Form("H%sTrkInPh%c%c%d", mc?"MC":"", chName[ich], ptName[ipt], 0)))){
        for(Int_t isp(1); isp<AliPID::kSPECIES; isp++){
          if(!(pr1 = (AliTRDresolutionProjection*)php.FindObject(Form("H%sTrkInPh%c%c%d", mc?"MC":"", chName[ich], ptName[ipt], isp)))) continue;
          (*pr0)+=(*pr1);
        }
        pr0->fH->SetNameTitle(Form("H%sTrkInPh%c%c", mc?"MC":"", chName[ich], ptName[ipt]),
                                  Form("TrackIn[%c]:: #Delta#phi{%s}", chSgn[ich], ptCut[ipt]));
        if((h2 = pr0->Projection2D(kNstat, kNcontours, 1))) arr->AddAt(h2, jh++);
        if(ipt && (pr1 = (AliTRDresolutionProjection*)php.FindObject(Form("H%sTrkInPh%c%c%d", mc?"MC":"", chName[ich], ptName[0], 0)))) (*pr1)+=(*pr0);
      }
      /*!dx*/
      if((pr0 = (AliTRDresolutionProjection*)php.FindObject(Form("H%sTrkInX%c%c%d", mc?"MC":"", chName[ich], ptName[ipt], 0)))){
        for(Int_t isp(1); isp<AliPID::kSPECIES; isp++){
          if(!(pr1 = (AliTRDresolutionProjection*)php.FindObject(Form("H%sTrkInX%c%c%d", mc?"MC":"", chName[ich], ptName[ipt], isp)))) continue;
          (*pr0)+=(*pr1);
        }
        pr0->fH->SetNameTitle(Form("H%sTrkInX%c%c", mc?"MC":"", chName[ich], ptName[ipt]),
                                  Form("TrackIn[%c]:: #Deltax{%s}", chSgn[ich], ptCut[ipt]));
        if((h2 = pr0->Projection2D(kNstat, kNcontours, 1))) arr->AddAt(h2, jh++);
        if(!ipt){
          xlow[ich] = (*pr0);
          xlow[ich].SetNameTitle(Form("H%sTrkInX%c%c%d", mc?"MC":"", chName[ich], ptName[0], 5),
                                 Form("TrackIn[%c]:: #Deltax{%s}", chSgn[ich], ptCut[0]));
          php.AddLast(&xlow[ich]);
        }
        if(ipt && (pr1 = (AliTRDresolutionProjection*)php.FindObject(Form("H%sTrkInX%c%c%d", mc?"MC":"", chName[ich], ptName[0], 0)))) (*pr1)+=(*pr0);
      }
    } // end pt loop for PID integration

    /*!dy*/
    if((pr0 = (AliTRDresolutionProjection*)php.FindObject(Form("H%sTrkInY%c%c%d", mc?"MC":"", chName[ich], ptName[0], 0)))){
      pr0->fH->SetNameTitle(Form("H%sTrkInY%c", mc?"MC":"", chName[ich]),
                            Form("TrackIn[%c]:: #Deltay", chSgn[ich]));
      pr0->SetShowRange(-0.3, 0.3);
      if((h2 = pr0->Projection2D(kNstat, kNcontours, 1, kFALSE))) arr->AddAt(h2, jh++);
      if((h2 = (TH2*)gDirectory->Get(Form("%sEn", pr0->fH->GetName())))) arr->AddAt(h2, jh++);
      if(ich && (pr1 = (AliTRDresolutionProjection*)php.FindObject(Form("H%sTrkInY%c%c%d", mc?"MC":"", chName[0], ptName[0], 0)))) (*pr1)+=(*pr0);
    }
    /*!dy high pt*/
    if(ich && (pr0 = (AliTRDresolutionProjection*)php.FindObject(Form("H%sTrkInY%c%c%d", mc?"MC":"", chName[0], ptName[2], 0)))){
      if((pr1 = (AliTRDresolutionProjection*)php.FindObject(Form("H%sTrkInY%c%c%d", mc?"MC":"", chName[ich], ptName[2], 0)))){
        (*pr0)+=(*pr1);
        pr0->fH->SetNameTitle(Form("H%sTrkInY%c", mc?"MC":"", ptName[2]), Form("TrackIn :: #Deltay{%s}", ptCut[2]));
        pr0->SetShowRange(-0.3, 0.3);
        if((h2 = pr0->Projection2D(kNstat, kNcontours, 1))) arr->AddAt(h2, jh++);
      }
    }
    /*!dphi*/
    if((pr0 = (AliTRDresolutionProjection*)php.FindObject(Form("H%sTrkInPh%c%c%d", mc?"MC":"", chName[ich], ptName[0], 0)))){
      pr0->fH->SetNameTitle(Form("H%sTrkInPh%c", mc?"MC":"", chName[ich]),
                            Form("TrackIn[%c]:: #Delta#phi", chSgn[ich]));
      pr0->SetShowRange(-1., 1.);
      if((h2 = pr0->Projection2D(kNstat, kNcontours, 1))) arr->AddAt(h2, jh++);
      if(ich==1 && (pr1 = (AliTRDresolutionProjection*)php.FindObject(Form("H%sTrkInPh%c%c%d", mc?"MC":"", chName[0], ptName[0], 0)))) (*pr1)+=(*pr0);
    }
    /*!dx*/
    if((pr0 = (AliTRDresolutionProjection*)php.FindObject(Form("H%sTrkInX%c%c%d", mc?"MC":"", chName[ich], ptName[0], 0)))){
      pr0->fH->SetNameTitle(Form("H%sTrkInX%c", mc?"MC":"", chName[ich]),
                            Form("TrackIn[%c]:: #Deltax", chSgn[ich]));
      if((h2 = pr0->Projection2D(kNstat, kNcontours, 1))) arr->AddAt(h2, jh++);
      if(ich==1 && (pr1 = (AliTRDresolutionProjection*)php.FindObject(Form("H%sTrkInX%c%c%d", mc?"MC":"", chName[0], ptName[0], 0)))) (*pr1)+=(*pr0);
    }
    /*!dx low pt*/
    if(ich && (pr0 = (AliTRDresolutionProjection*)php.FindObject(Form("H%sTrkInX%c%c%d", mc?"MC":"", chName[0], ptName[0], 5)))){
      if((pr1 = (AliTRDresolutionProjection*)php.FindObject(Form("H%sTrkInX%c%c%d", mc?"MC":"", chName[ich], ptName[0], 5)))){
        (*pr0)+=(*pr1);
        pr0->fH->SetNameTitle(Form("H%sTrkInX%c", mc?"MC":"", ptName[0]), Form("TrackIn :: #Deltax{%s}", ptCut[0]));
        if((h2 = pr0->Projection2D(kNstat, kNcontours, 1))) arr->AddAt(h2, jh++);
      }
    }
  } // end charge loop

  for(Int_t isp(0); isp<AliPID::kSPECIES; isp++){
    /*!dy*/
    if((pr0 = (AliTRDresolutionProjection*)php.FindObject(Form("H%sTrkInY%c%d", mc?"MC":"", chName[0], isp)))){
      pr0->fH->SetNameTitle(Form("H%sTrkInY%d", mc?"MC":"", isp), Form("TrackIn[%s] :: #Deltay", AliPID::ParticleLatexName(isp)));
      pr0->SetShowRange(-0.3, 0.3);
      if((h2 = pr0->Projection2D(kNstat, kNcontours, 1, kFALSE))) arr->AddAt(h2, jh++);
      if((h2 = (TH2*)gDirectory->Get(Form("%sEn", pr0->fH->GetName())))) arr->AddAt(h2, jh++);
    }
    /*!dphi*/
    if((pr0 = (AliTRDresolutionProjection*)php.FindObject(Form("H%sTrkInPh%c%d", mc?"MC":"", chName[0], isp)))){
      pr0->fH->SetNameTitle(Form("H%sTrkInPh%d", mc?"MC":"", isp), Form("TrackIn[%s] :: #Delta#phi", AliPID::ParticleLatexName(isp)));
      pr0->SetShowRange(-1., 1.);
      if((h2 = pr0->Projection2D(kNstat, kNcontours, 1))) arr->AddAt(h2, jh++);
    }
    /*!dQdl*/
    if((pr0 = (AliTRDresolutionProjection*)php.FindObject(Form("H%sTrkInQ%c%d", mc?"MC":"", chName[0], isp)))){
      pr0->fH->SetNameTitle(Form("H%sTrkInQ%d", mc?"MC":"", isp), Form("TrackIn[%s] :: dQdl", AliPID::ParticleLatexName(isp)));
      pr0->SetShowRange(-2.2, -1.75);
      if((h2 = pr0->Projection2D(kNstat, kNcontours, 2))) arr->AddAt(h2, jh++);
      pr0->fH->SetName(Form("H%sTrkInQS%d", mc?"MC":"", isp));
      pr0->SetShowRange(-1.7, -1.25);
      if((h2 = pr0->Projection2D(kNstat, kNcontours, 0))) arr->AddAt(h2, jh++);
    }
  } // end PID processing

  /*!dy*/
  if((pr0 = (AliTRDresolutionProjection*)php.FindObject(Form("H%sTrkInY%c%c%d", mc?"MC":"", chName[0], ptName[0], 0)))){
    pr0->fH->SetNameTitle(Form("H%sTrkInY", mc?"MC":""), "TrackIn :: #Deltay");
    pr0->SetShowRange(-0.3, 0.3);
    if((h2 = pr0->Projection2D(kNstat, kNcontours, 1, kFALSE))) arr->AddAt(h2, jh++);
    if((h2 = (TH2*)gDirectory->Get(Form("%sEn", pr0->fH->GetName())))) arr->AddAt(h2, jh++);
  }
  /*!dphi*/
  if((pr0 = (AliTRDresolutionProjection*)php.FindObject(Form("H%sTrkInPh%c%c%d", mc?"MC":"", chName[0], ptName[0], 0)))){
    pr0->fH->SetNameTitle(Form("H%sTrkInPh", mc?"MC":""), "TrackIn :: #Delta#phi");
    if((h2 = pr0->Projection2D(kNstat, kNcontours, 1))) arr->AddAt(h2, jh++);
  }
  /*!dx*/
  if((pr0 = (AliTRDresolutionProjection*)php.FindObject(Form("H%sTrkInX%c%c%d", mc?"MC":"", chName[0], ptName[0], 0)))){
    pr0->fH->SetNameTitle(Form("H%sTrkInX", mc?"MC":""), "TrackIn :: #Deltax");
    if((h2 = pr0->Projection2D(kNstat, kNcontours, 1))) arr->AddAt(h2, jh++);
  }

  // Row Cross processing
  /*!RC dz*/
  if((pr0 = (AliTRDresolutionProjection*)php.FindObject(Form("H%sTrkInRCZ%c", mc?"MC":"", ptName[0])))){
    for(Int_t ipt(0); ipt<kNpt; ipt++){
      if(!(pr1 = (AliTRDresolutionProjection*)php.FindObject(Form("H%sTrkInRCZ%c", mc?"MC":"", ptName[ipt])))) continue;
      (*pr0)+=(*pr1);
    }
    pr0->fH->SetNameTitle(Form("H%sTrkInRCZ", mc?"MC":""), "TrackIn[RC]:: #Deltaz");
    if((h2 = pr0->Projection2D(kNstat, kNcontours, 1))) arr->AddAt(h2, jh++);
  }
  /*!RC dy*/
  if((pr0 = (AliTRDresolutionProjection*)php.FindObject(Form("H%sTrkInRCY%c", mc?"MC":"", ptName[0])))){
    for(Int_t ipt(0); ipt<kNpt; ipt++){
      if(!(pr1 = (AliTRDresolutionProjection*)php.FindObject(Form("H%sTrkInRCY%c", mc?"MC":"", ptName[ipt])))) continue;
      (*pr0)+=(*pr1);
    }
    pr0->fH->SetNameTitle(Form("H%sTrkInRCY", mc?"MC":""), "TrackIn[RC]:: #Deltay");
    if((h2 = pr0->Projection2D(kNstat, kNcontours, 1))) arr->AddAt(h2, jh++);
  }
  /*!RC dphi*/
  if((pr0 = (AliTRDresolutionProjection*)php.FindObject(Form("H%sTrkInRCPh%c", mc?"MC":"", ptName[0])))){
    for(Int_t ipt(0); ipt<kNpt; ipt++){
      if(!(pr1 = (AliTRDresolutionProjection*)php.FindObject(Form("H%sTrkInRCPh%c", mc?"MC":"", ptName[ipt])))) continue;
      (*pr0)+=(*pr1);
    }
    pr0->fH->SetNameTitle(Form("H%sTrkInRCPh", mc?"MC":""), "TrackIn[RC]:: #Delta#phi");
    if((h2 = pr0->Projection2D(kNstat, kNcontours, 1))) arr->AddAt(h2, jh++);
  }
  /*!RC dx*/
  if((pr0 = (AliTRDresolutionProjection*)php.FindObject(Form("H%sTrkInRCX%c", mc?"MC":"", ptName[0])))){
    for(Int_t ipt(0); ipt<kNpt; ipt++){
      if(!(pr1 = (AliTRDresolutionProjection*)php.FindObject(Form("H%sTrkInRCX%c", mc?"MC":"", ptName[ipt])))) continue;
      (*pr0)+=(*pr1);
    }
    pr0->fH->SetNameTitle(Form("H%sTrkInRCX", mc?"MC":""), "TrackIn[RC]:: #Deltax");
    if((h2 = pr0->Projection2D(kNstat, kNcontours, 1))) arr->AddAt(h2, jh++);
  }
  AliInfo(Form("Done %3d 2D projections.", jh));
  return kTRUE;
}


//________________________________________________________
Bool_t AliTRDresolution::MakeProjectionTrack()
{
// Analyse tracklet
  const Int_t kNcontours(9);
  const Int_t kNstat(30);
  Int_t cidx(kMCtrack);
  if(fProj && fProj->At(cidx)) return kTRUE;
  if(!fContainer){
    AliError("Missing data container.");
    return kFALSE;
  }
  THnSparse *H(NULL);
  if(!(H = (THnSparse*)fContainer->FindObject("hTRD2MC"))){
    AliError("Missing/Wrong data @ hTRD2MC.");
    return kFALSE;
  }
  Int_t ndim(H->GetNdimensions());
  Int_t coord[kNdim+1]; memset(coord, 0, sizeof(Int_t) * (kNdim+1)); Double_t v = 0.;
  TAxis *aa[kNdim+1], *as(NULL), *ap(NULL); memset(aa, 0, sizeof(TAxis*) * (kNdim+1));
  for(Int_t id(0); id<ndim; id++) aa[id] = H->GetAxis(id);
  if(ndim > kSpeciesChgRC) as = H->GetAxis(kSpeciesChgRC);
  if(ndim > kPt) ap = H->GetAxis(kPt);

  // build list of projections
  const Int_t nsel(AliTRDgeometry::kNlayer*kNpt*AliPID::kSPECIES*7);//, npsel(3);
  const Char_t chName[kNcharge] = {'n', 'p'};const Char_t chSgn[kNcharge] = {'-', '+'};
  const Char_t ptName[kNpt] = {'l', 'i', 'h'};
  const Char_t *ptCut[kNpt] = {"p_{t}[GeV/c]<0.8", "0.8<=p_{t}[GeV/c]<1.5", "p_{t}[GeV/c]>=1.5"};
  // define rebinning strategy
  const Int_t nEtaPhi(4); Int_t rebinEtaPhiX[nEtaPhi] = {1, 2, 5, 1}, rebinEtaPhiY[nEtaPhi] = {2, 1, 1, 5};
  AliTRDresolutionProjection hp[kTrkNproj]; TObjArray php(kTrkNproj);
  Int_t ih(0), isel(-1), np[nsel]; memset(np, 0, nsel*sizeof(Int_t));
  for(Int_t ily(0); ily<AliTRDgeometry::kNlayer; ily++){
    for(Int_t ipt(0); ipt<kNpt; ipt++){
      for(Int_t isp(0); isp<AliPID::kSPECIES; isp++){
        for(Int_t ich(0); ich<kNcharge; ich++){
          isel++; // new selection
          hp[ih].Build(Form("HMCTrkY%c%c%d%d", chName[ich], ptName[ipt], isp, ily),
                       Form("Tracks[%s%c]:: #Deltay{%s} Ly[%d]", AliPID::ParticleLatexName(isp), chSgn[ich], ptCut[ipt], ily),
                       kEta, kPhi, kYrez, aa);
          hp[ih].SetRebinStrategy(nEtaPhi, rebinEtaPhiX, rebinEtaPhiY);
          php.AddLast(&hp[ih++]); np[isel]++;
          hp[ih].Build(Form("HMCTrkPh%c%c%d%d", chName[ich], ptName[ipt], isp, ily),
                       Form("Tracks[%s%c]:: #Delta#phi{%s} Ly[%d]", AliPID::ParticleLatexName(isp), chSgn[ich], ptCut[ipt], ily),
                       kEta, kPhi, kPrez, aa);
          hp[ih].SetRebinStrategy(nEtaPhi, rebinEtaPhiX, rebinEtaPhiY);
          php.AddLast(&hp[ih++]); np[isel]++;
          hp[ih].Build(Form("HMCTrkDPt%c%c%d%d", chName[ich], ptName[ipt], isp, ily),
                       Form("Tracks[%s%c]:: #Deltap_{t}/p_{t}{%s} Ly[%d]", AliPID::ParticleLatexName(isp), chSgn[ich], ptCut[ipt], ily),
                       kEta, kPhi, kNdim, aa);
          hp[ih].SetShowRange(0.,10.);
          hp[ih].SetRebinStrategy(nEtaPhi, rebinEtaPhiX, rebinEtaPhiY);
          php.AddLast(&hp[ih++]); np[isel]++;
        }
      }
      isel++; // new selection
      hp[ih].Build(Form("HMCTrkZ%c%d", ptName[ipt], ily),
                    Form("Tracks[RC]:: #Deltaz{%s} Ly[%d]", ptCut[ipt], ily),
                    kEta, kPhi, kZrez, aa);
      hp[ih].SetRebinStrategy(nEtaPhi, rebinEtaPhiX, rebinEtaPhiY);
      php.AddLast(&hp[ih++]); np[isel]++;
    }
  }

  Int_t ly(0), ch(0), pt(0), sp(2), rcBin(as?as->FindBin(0.):-1);
  for (Long64_t ib(0); ib < H->GetNbins(); ib++) {
    v = H->GetBinContent(ib, coord);
    if(v<1.) continue;
    ly = coord[kBC]-1; // layer selection
    // charge selection
    ch=0; sp=2;// [pi-] track [dafault]
    if(rcBin>0){ // debug mode in which species are also saved
      sp = Int_t(TMath::Abs(as->GetBinCenter(coord[kSpeciesChgRC])))-1;
      if(coord[kSpeciesChgRC] > rcBin) ch = 1;        // [+] track
      else if(coord[kSpeciesChgRC] == rcBin) ch = 2;  // [RC] track
    }
    // pt selection
    pt = 0; // low pt [default]
    if(ap) pt = coord[kPt]-1;
    // global selection
    Int_t ioff = ly*kNpt*31+pt*31; ioff+=3*(sp<0?10:(sp*kNcharge+ch));
    isel = ly*kNpt*11+pt*11; isel+=sp<0?10:(sp*kNcharge+ch);
    AliDebug(4, Form("SELECTION[%d] :: ch[%c] pt[%c] sp[%d] ly[%d]\n", np[isel], ch==2?'Z':chName[ch], ptName[pt], sp, ly));
    for(Int_t jh(0); jh<np[isel]; jh++) ((AliTRDresolutionProjection*)php.At(ioff+jh))->Increment(coord, v);
  }
  TObjArray *arr(NULL);
  fProj->AddAt(arr = new TObjArray(kTrkNproj), cidx);

  TH2 *h2(NULL); Int_t jh(0);
  for(; ih--; ){
    if(!hp[ih].fH) continue;
    if(!(h2 = hp[ih].Projection2D(kNstat, kNcontours))) continue;
    arr->AddAt(h2, jh++);
  }

  // combine up the tree of projections
  AliTRDresolutionProjection *pr0(NULL), *pr1(NULL);
  //Int_t iproj(0), jproj(0);
  for(Int_t ily(0); ily<AliTRDgeometry::kNlayer; ily++){
    for(Int_t ich(0); ich<kNcharge; ich++){
      for(Int_t ipt(0); ipt<kNpt; ipt++){
        /*!dy*/
        if((pr0 = (AliTRDresolutionProjection*)php.FindObject(Form("HMCTrkY%c%c%d%d", chName[ich], ptName[ipt], 0, ily)))){
          for(Int_t isp(1); isp<AliPID::kSPECIES; isp++){
            if(!(pr1 = (AliTRDresolutionProjection*)php.FindObject(Form("HMCTrkY%c%c%d%d", chName[ich], ptName[ipt], isp, ily)))) continue;
            (*pr0)+=(*pr1);
          }
          AliDebug(2, Form("Rename %s to HMCTrkY%c%c%d", pr0->fH->GetName(), chName[ich], ptName[ipt], ily));
          pr0->fH->SetNameTitle(Form("HMCTrkY%c%c%d", chName[ich], ptName[ipt], ily),
                                    Form("Tracks[%c]:: #Deltay{%s} Ly[%d]", chSgn[ich], ptCut[ipt], ily));
          if((h2 = pr0->Projection2D(kNstat, kNcontours, 1))) arr->AddAt(h2, jh++);
          if(ipt && (pr1 = (AliTRDresolutionProjection*)php.FindObject(Form("HMCTrkY%c%c%d%d", chName[ich], ptName[0], 0, ily)))) (*pr1)+=(*pr0);
        }
        /*!dphi*/
        if((pr0 = (AliTRDresolutionProjection*)php.FindObject(Form("HMCTrkPh%c%c%d%d", chName[ich], ptName[ipt], 0, ily)))){
          for(Int_t isp(1); isp<AliPID::kSPECIES; isp++){
            if(!(pr1 = (AliTRDresolutionProjection*)php.FindObject(Form("HMCTrkPh%c%c%d%d", chName[ich], ptName[ipt], isp, ily)))) continue;
            (*pr0)+=(*pr1);
          }
          AliDebug(2, Form("Rename %s to HMCTrkPh%c%c%d", pr0->fH->GetName(), chName[ich], ptName[ipt], ily));
          pr0->fH->SetNameTitle(Form("HMCTrkPh%c%c%d", chName[ich], ptName[ipt], ily),
                                    Form("Tracks[%c]:: #Delta#phi{%s} Ly[%d]", chSgn[ich], ptCut[ipt], ily));
          if((h2 = pr0->Projection2D(kNstat, kNcontours, 1))) arr->AddAt(h2, jh++);
          if(ipt && (pr1 = (AliTRDresolutionProjection*)php.FindObject(Form("HMCTrkPh%c%c%d%d", chName[ich], ptName[0], 0, ily)))) (*pr1)+=(*pr0);
        }

        /*!dpt/pt*/
        if((pr0 = (AliTRDresolutionProjection*)php.FindObject(Form("HMCTrkDPt%c%c%d%d", chName[ich], ptName[ipt], 0, ily)))){
          for(Int_t isp(1); isp<AliPID::kSPECIES; isp++){
            if(!(pr1 = (AliTRDresolutionProjection*)php.FindObject(Form("HMCTrkDPt%c%c%d%d", chName[ich], ptName[ipt], isp, ily)))) continue;
            (*pr0)+=(*pr1);
          }
          AliDebug(2, Form("Rename %s to HMCTrkDPt%c%c%d", pr0->fH->GetName(), chName[ich], ptName[ipt], ily));
          pr0->fH->SetNameTitle(Form("HMCTrkDPt%c%c%d", chName[ich], ptName[ipt], ily),
                                    Form("Tracks[%c]:: #Deltap_{t}/p_{t}{%s} Ly[%d]", chSgn[ich], ptCut[ipt], ily));
          if((h2 = pr0->Projection2D(kNstat, kNcontours, 1))) arr->AddAt(h2, jh++);
          if(ipt && (pr1 = (AliTRDresolutionProjection*)php.FindObject(Form("HMCTrkDPt%c%c%d%d", chName[ich], ptName[0], 0, ily)))) (*pr1)+=(*pr0);
        }
      }
      /*!dy*/
      if((pr0 = (AliTRDresolutionProjection*)php.FindObject(Form("HMCTrkY%c%c%d%d", chName[ich], ptName[0], 0, ily)))){
        pr0->fH->SetNameTitle(Form("HMCTrkY%c%d", chName[ich], ily),
                              Form("Tracks[%c]:: #Deltay Ly[%d]", chSgn[ich], ily));
        if((h2 = pr0->Projection2D(kNstat, kNcontours, 1))) arr->AddAt(h2, jh++);
        if(ich==1 && (pr1 = (AliTRDresolutionProjection*)php.FindObject(Form("HMCTrkY%c%c%d%d", chName[0], ptName[0], 0, ily)))) (*pr1)+=(*pr0);
      }
      /*!dphi*/
      if((pr0 = (AliTRDresolutionProjection*)php.FindObject(Form("HMCTrkPh%c%c%d%d", chName[ich], ptName[0], 0, ily)))){
        pr0->fH->SetNameTitle(Form("HMCTrkPh%c%d", chName[ich], ily),
                              Form("Tracks[%c]:: #Delta#phi Ly[%d]", chSgn[ich], ily));
        if((h2 = pr0->Projection2D(kNstat, kNcontours, 1))) arr->AddAt(h2, jh++);
        if(ich==1 && (pr1 = (AliTRDresolutionProjection*)php.FindObject(Form("HMCTrkPh%c%c%d%d", chName[0], ptName[0], 0, ily)))) (*pr1)+=(*pr0);
      }
      /*!dpt/pt*/
      if((pr0 = (AliTRDresolutionProjection*)php.FindObject(Form("HMCTrkDPt%c%c%d%d", chName[ich], ptName[0], 0, ily)))){
        pr0->fH->SetNameTitle(Form("HMCTrkDPt%c%d", chName[ich], ily),
                              Form("Tracks[%c]:: #Deltap_{t}/p_{t} Ly[%d]", chSgn[ich], ily));
        if((h2 = pr0->Projection2D(kNstat, kNcontours, 1))) arr->AddAt(h2, jh++);
        if(ich==1 && (pr1 = (AliTRDresolutionProjection*)php.FindObject(Form("HMCTrkDPt%c%c%d%d", chName[0], ptName[0], 0, ily)))) (*pr1)+=(*pr0);
      }
    }
    /*!dy*/
    if((pr0 = (AliTRDresolutionProjection*)php.FindObject(Form("HMCTrkY%c%c%d%d", chName[0], ptName[0], 0, ily)))){
      pr0->fH->SetNameTitle(Form("HMCTrkY%d", ily), Form("Tracks :: #Deltay Ly[%d]", ily));
      if((h2 = pr0->Projection2D(kNstat, kNcontours, 1))) arr->AddAt(h2, jh++);
    }
    /*!dphi*/
    if((pr0 = (AliTRDresolutionProjection*)php.FindObject(Form("HMCTrkPh%c%c%d%d", chName[0], ptName[0], 0, ily)))){
      pr0->fH->SetNameTitle(Form("HMCTrkPh%d", ily), Form("Tracks :: #Delta#phi Ly[%d]", ily));
      if((h2 = pr0->Projection2D(kNstat, kNcontours, 1))) arr->AddAt(h2, jh++);
    }
    /*!dpt/pt*/
    if((pr0 = (AliTRDresolutionProjection*)php.FindObject(Form("HMCTrkDPt%c%c%d%d", chName[0], ptName[0], 0, ily)))){
      pr0->fH->SetNameTitle(Form("HMCTrkDPt%d", ily), Form("Tracks :: #Deltap_{t}/p_{t} Ly[%d]", ily));
      if((h2 = pr0->Projection2D(kNstat, kNcontours, 1))) arr->AddAt(h2, jh++);
    }
  }
  AliInfo(Form("Done %3d 2D projections.", jh));
  return kTRUE;
}

//________________________________________________________
Bool_t AliTRDresolution::PostProcess()
{
// Fit, Project, Combine, Extract values from the containers filled during execution

  if (!fContainer) {
    AliError("ERROR: list not available");
    return kFALSE;
  }
  if(!fProj){
    AliInfo("Building array of projections ...");
    fProj = new TObjArray(kNclasses); fProj->SetOwner(kTRUE);
  }

  //PROCESS EXPERIMENTAL DISTRIBUTIONS
  // Clusters detector
  if(HasProcessDetector() && !MakeProjectionDetector()) return kFALSE;
  // Clusters residuals
  if(HasProcessCluster() && !MakeProjectionCluster()) return kFALSE;
  fNRefFigures = 3;
  // Tracklet residual/pulls
  if(HasProcessTrklt() && !MakeProjectionTracklet()) return kFALSE;
  fNRefFigures = 7;
  // TRDin residual/pulls
  if(HasProcessTrkIn() && !MakeProjectionTrackIn()) return kFALSE;
  fNRefFigures = 11;

  if(!HasMCdata()) return kTRUE;
  //PROCESS MC RESIDUAL DISTRIBUTIONS

  // CLUSTER Y RESOLUTION/PULLS
  if(!MakeProjectionCluster(kTRUE)) return kFALSE;
  fNRefFigures = 17;

  // TRACKLET RESOLUTION/PULLS
  if(!MakeProjectionTracklet(kTRUE)) return kFALSE;
  fNRefFigures = 21;

  // TRACK RESOLUTION/PULLS
  if(!MakeProjectionTrack()) return kFALSE;
  fNRefFigures+=16;

  // TRACK TRDin RESOLUTION/PULLS
  if(!MakeProjectionTrackIn(kTRUE)) return kFALSE;
  fNRefFigures+=8;

  return kTRUE;
}


//________________________________________________________
void AliTRDresolution::Terminate(Option_t *opt)
{
  AliTRDrecoTask::Terminate(opt);
  if(HasPostProcess()) PostProcess();
}

//________________________________________________________
void AliTRDresolution::AdjustF1(TH1 *h, TF1 *f)
{
// Helper function to avoid duplication of code
// Make first guesses on the fit parameters

  // find the intial parameters of the fit !! (thanks George)
  Int_t nbinsy = Int_t(.5*h->GetNbinsX());
  Double_t sum = 0.;
  for(Int_t jbin=nbinsy-4; jbin<=nbinsy+4; jbin++) sum+=h->GetBinContent(jbin); sum/=9.;
  f->SetParLimits(0, 0., 3.*sum);
  f->SetParameter(0, .9*sum);
  Double_t rms = h->GetRMS();
  f->SetParLimits(1, -rms, rms);
  f->SetParameter(1, h->GetMean());

  f->SetParLimits(2, 0., 2.*rms);
  f->SetParameter(2, rms);
  if(f->GetNpar() <= 4) return;

  f->SetParLimits(3, 0., sum);
  f->SetParameter(3, .1*sum);

  f->SetParLimits(4, -.3, .3);
  f->SetParameter(4, 0.);

  f->SetParLimits(5, 0., 1.e2);
  f->SetParameter(5, 2.e-1);
}

//________________________________________________________
TObjArray* AliTRDresolution::BuildMonitorContainerCluster(const char* name, Bool_t expand, Float_t range)
{
// Build performance histograms for AliTRDcluster.vs TRD track or MC
//  - y reziduals/pulls

  TObjArray *arr = new TObjArray(2);
  arr->SetName(name); arr->SetOwner();
  TH1 *h(NULL); char hname[100], htitle[300];

  // tracklet resolution/pull in y direction
  snprintf(hname, 100, "%s_%s_Y", GetNameId(), name);
  snprintf(htitle, 300, "Y res for \"%s\" @ %s;tg(#phi);#Delta y [cm];%s", GetNameId(), name, "Detector");
  Float_t rr = range<0.?fDyRange:range;
  if(!(h = (TH3S*)gROOT->FindObject(hname))){
    Int_t nybins=50;
    if(expand) nybins*=2;
    h = new TH3S(hname, htitle,
                 48, -.48, .48,            // phi
                 60, -rr, rr,              // dy
                 nybins, -0.5, nybins-0.5);// segment
  } else h->Reset();
  arr->AddAt(h, 0);
  snprintf(hname, 100, "%s_%s_YZpull", GetNameId(), name);
  snprintf(htitle, 300, "YZ pull for \"%s\" @ %s;%s;#Delta y  / #sigma_{y};#Delta z  / #sigma_{z}", GetNameId(), name, "Detector");
  if(!(h = (TH3S*)gROOT->FindObject(hname))){
    h = new TH3S(hname, htitle, 540, -0.5, 540-0.5, 100, -4.5, 4.5, 100, -4.5, 4.5);
  } else h->Reset();
  arr->AddAt(h, 1);

  return arr;
}

//________________________________________________________
TObjArray* AliTRDresolution::BuildMonitorContainerTracklet(const char* name, Bool_t expand)
{
// Build performance histograms for AliExternalTrackParam.vs TRD tracklet
//  - y reziduals/pulls
//  - z reziduals/pulls
//  - phi reziduals
  TObjArray *arr = BuildMonitorContainerCluster(name, expand, 0.05);
  arr->Expand(5);
  TH1 *h(NULL); char hname[100], htitle[300];

  // tracklet resolution/pull in z direction
  snprintf(hname, 100, "%s_%s_Z", GetNameId(), name);
  snprintf(htitle, 300, "Z res for \"%s\" @ %s;tg(#theta);#Delta z [cm]", GetNameId(), name);
  if(!(h = (TH2S*)gROOT->FindObject(hname))){
    h = new TH2S(hname, htitle, 50, -1., 1., 100, -.05, .05);
  } else h->Reset();
  arr->AddAt(h, 2);
  snprintf(hname, 100, "%s_%s_Zpull", GetNameId(), name);
  snprintf(htitle, 300, "Z pull for \"%s\" @ %s;tg(#theta);#Delta z  / #sigma_{z};row cross", GetNameId(), name);
  if(!(h = (TH3S*)gROOT->FindObject(hname))){
    h = new TH3S(hname, htitle, 50, -1., 1., 100, -5.5, 5.5, 2, -0.5, 1.5);
    h->GetZaxis()->SetBinLabel(1, "no RC");
    h->GetZaxis()->SetBinLabel(2, "RC");
  } else h->Reset();
  arr->AddAt(h, 3);

  // tracklet to track phi resolution
  snprintf(hname, 100, "%s_%s_PHI", GetNameId(), name);
  snprintf(htitle, 300, "#Phi res for \"%s\" @ %s;tg(#phi);#Delta #phi [rad];%s", GetNameId(), name, "Detector");
  Int_t nsgms=540;
  if(!(h = (TH3S*)gROOT->FindObject(hname))){
    h = new TH3S(hname, htitle, 48, -.48, .48, 100, -.5, .5, nsgms, -0.5, nsgms-0.5);
  } else h->Reset();
  arr->AddAt(h, 4);

  return arr;
}

//________________________________________________________
TObjArray* AliTRDresolution::BuildMonitorContainerTrack(const char* name)
{
// Build performance histograms for AliExternalTrackParam.vs MC
//  - y resolution/pulls
//  - z resolution/pulls
//  - phi resolution, snp pulls
//  - theta resolution, tgl pulls
//  - pt resolution, 1/pt pulls, p resolution

  TObjArray *arr = BuildMonitorContainerTracklet(name);
  arr->Expand(11);
  TH1 *h(NULL); char hname[100], htitle[300];
  //TAxis *ax(NULL);

  // snp pulls
  snprintf(hname, 100, "%s_%s_SNPpull", GetNameId(), name);
  snprintf(htitle, 300, "SNP pull for \"%s\" @ %s;tg(#phi);#Delta snp  / #sigma_{snp};entries", GetNameId(), name);
  if(!(h = (TH2I*)gROOT->FindObject(hname))){
    h = new TH2I(hname, htitle, 60, -.3, .3, 100, -4.5, 4.5);
  } else h->Reset();
  arr->AddAt(h, 5);

  // theta resolution
  snprintf(hname, 100, "%s_%s_THT", GetNameId(), name);
  snprintf(htitle, 300, "#Theta res for \"%s\" @ %s;tg(#theta);#Delta #theta [rad];entries", GetNameId(), name);
  if(!(h = (TH2I*)gROOT->FindObject(hname))){
    h = new TH2I(hname, htitle, 100, -1., 1., 100, -5e-3, 5e-3);
  } else h->Reset();
  arr->AddAt(h, 6);
  // tgl pulls
  snprintf(hname, 100, "%s_%s_TGLpull", GetNameId(), name);
  snprintf(htitle, 300, "TGL pull for \"%s\" @ %s;tg(#theta);#Delta tgl  / #sigma_{tgl};entries", GetNameId(), name);
  if(!(h = (TH2I*)gROOT->FindObject(hname))){
    h = new TH2I(hname, htitle, 100, -1., 1., 100, -4.5, 4.5);
  } else h->Reset();
  arr->AddAt(h, 7);

  const Int_t kNdpt(150);
  const Int_t kNspc = 2*AliPID::kSPECIES+1;
  Float_t lPt=0.1, lDPt=-.1, lSpc=-5.5;
  Float_t binsPt[kNpt+1], binsSpc[kNspc+1], binsDPt[kNdpt+1];
  for(Int_t i=0;i<kNpt+1; i++,lPt=TMath::Exp(i*.15)-1.) binsPt[i]=lPt;
  for(Int_t i=0; i<kNspc+1; i++,lSpc+=1.) binsSpc[i]=lSpc;
  for(Int_t i=0; i<kNdpt+1; i++,lDPt+=2.e-3) binsDPt[i]=lDPt;

  // Pt resolution
  snprintf(hname, 100, "%s_%s_Pt", GetNameId(), name);
  snprintf(htitle, 300, "#splitline{P_{t} res for}{\"%s\" @ %s};p_{t} [GeV/c];#Delta p_{t}/p_{t}^{MC};SPECIES", GetNameId(), name);
  if(!(h = (TH3S*)gROOT->FindObject(hname))){
    h = new TH3S(hname, htitle,
                 kNpt, binsPt, kNdpt, binsDPt, kNspc, binsSpc);
    //ax = h->GetZaxis();
    //for(Int_t ib(1); ib<=ax->GetNbins(); ib++) ax->SetBinLabel(ib, fgParticle[ib-1]);
  } else h->Reset();
  arr->AddAt(h, 8);
  // 1/Pt pulls
  snprintf(hname, 100, "%s_%s_1Pt", GetNameId(), name);
  snprintf(htitle, 300, "#splitline{1/P_{t} pull for}{\"%s\" @ %s};1/p_{t}^{MC} [c/GeV];#Delta(1/p_{t})/#sigma(1/p_{t});SPECIES", GetNameId(), name);
  if(!(h = (TH3S*)gROOT->FindObject(hname))){
    h = new TH3S(hname, htitle,
                 kNpt, 0., 2., 100, -4., 4., kNspc, -5.5, 5.5);
    //ax = h->GetZaxis();
    //for(Int_t ib(1); ib<=ax->GetNbins(); ib++) ax->SetBinLabel(ib, fgParticle[ib-1]);
  } else h->Reset();
  arr->AddAt(h, 9);
  // P resolution
  snprintf(hname, 100, "%s_%s_P", GetNameId(), name);
  snprintf(htitle, 300, "P res for \"%s\" @ %s;p [GeV/c];#Delta p/p^{MC};SPECIES", GetNameId(), name);
  if(!(h = (TH3S*)gROOT->FindObject(hname))){
    h = new TH3S(hname, htitle,
                 kNpt, binsPt, kNdpt, binsDPt, kNspc, binsSpc);
    //ax = h->GetZaxis();
    //for(Int_t ib(1); ib<=ax->GetNbins(); ib++) ax->SetBinLabel(ib, fgParticle[ib-1]);
  } else h->Reset();
  arr->AddAt(h, 10);

  return arr;
}


//________________________________________________________
TObjArray* AliTRDresolution::Histos()
{
  //
  // Define histograms
  //

  if(fContainer) return fContainer;

  fContainer  = new TObjArray(kNclasses); fContainer->SetOwner(kTRUE);
  THnSparse *H(NULL);
  const Int_t nhn(100); Char_t hn[nhn]; TString st;

  //++++++++++++++++++++++
  // cluster to detector
  snprintf(hn, nhn, "h%s", fgPerformanceName[kDetector]);
  if(!(H = (THnSparseI*)gROOT->FindObject(hn))){
    const Int_t mdim(6);
    const Char_t *cldTitle[mdim] = {"layer", fgkTitle[kPhi], "pad row", "centrality", "q [a.u.]", "no. of pads"/*, "tb [10 Hz]"*/};
    const Int_t cldNbins[mdim]   = {AliTRDgeometry::kNlayer, fgkNbins[kPhi], 76, AliTRDeventInfo::kCentralityClasses, 50, kNpads};
    const Double_t cldMin[mdim]  = {-0.5, fgkMin[kPhi], -0.5, -0.5, 0., 0.5},
                   cldMax[mdim]  = {AliTRDgeometry::kNlayer-0.5, fgkMax[kPhi], 75.5, AliTRDeventInfo::kCentralityClasses - 0.5, 1200., kNpads+.5};
    st = "cluster proprieties;";
    // define minimum info to be saved in non debug mode
    Int_t ndim=DebugLevel()>=1?mdim:Int_t(kNdimDet);
    for(Int_t idim(0); idim<ndim; idim++){ st += cldTitle[idim]; st+=";";}
    H = new THnSparseI(hn, st.Data(), ndim, cldNbins, cldMin, cldMax);
  } else H->Reset();
  fContainer->AddAt(H, kDetector);
  //++++++++++++++++++++++
  // cluster to tracklet residuals/pulls
  snprintf(hn, nhn, "h%s", fgPerformanceName[kCluster]);
  if(!(H = (THnSparseI*)gROOT->FindObject(hn))){
    const Int_t mdim(10);
    const Char_t *clTitle[mdim] = {"layer", fgkTitle[kPhi], fgkTitle[kEta], fgkTitle[kYrez], "#Deltax [cm]", "Q</Q", "Q/angle", "#Phi [deg]", "centrality", "no. of pads"};
    const Int_t clNbins[mdim]   = {AliTRDgeometry::kNlayer, fgkNbins[kPhi], fgkNbins[kEta], fgkNbins[kYrez], 45, 30, 30, 15, AliTRDeventInfo::kCentralityClasses, kNpads};
    const Double_t clMin[mdim]  = {-0.5, fgkMin[kPhi], fgkMin[kEta], fgkMin[kYrez]/10., -.5, 0.1, -2., -45, -0.5, 0.5},
                   clMax[mdim]  = {AliTRDgeometry::kNlayer-0.5, fgkMax[kPhi], fgkMax[kEta], fgkMax[kYrez]/10., 4., 2.1, 118., 45, AliTRDeventInfo::kCentralityClasses - 0.5, kNpads+.5};
    st = "cluster spatial&charge resolution;";
    // define minimum info to be saved in non debug mode
    Int_t ndim=DebugLevel()>=1?mdim:Int_t(kNdimCl);
    for(Int_t idim(0); idim<ndim; idim++){ st += clTitle[idim]; st+=";";}
    H = new THnSparseI(hn, st.Data(), ndim, clNbins, clMin, clMax);
  } else H->Reset();
  fContainer->AddAt(H, kCluster);
  //++++++++++++++++++++++
  // tracklet to TRD track
  snprintf(hn, nhn, "h%s", fgPerformanceName[kTracklet]);
  if(!(H = (THnSparseI*)gROOT->FindObject(hn))){
    const Int_t mdim(kNdim+8);
    Char_t *trTitle[mdim]; memcpy(trTitle, fgkTitle, kNdim*sizeof(Char_t*));
    Int_t trNbins[mdim]; memcpy(trNbins, fgkNbins, kNdim*sizeof(Int_t));
    Double_t trMin[mdim]; memcpy(trMin, fgkMin, kNdim*sizeof(Double_t));
    Double_t trMax[mdim]; memcpy(trMax, fgkMax, kNdim*sizeof(Double_t));
    // set specific fields
    trMin[kYrez] = -0.45; trMax[kYrez] = 0.45;
    trMin[kPrez] = -4.5; trMax[kPrez] = 4.5;
    trMin[kZrez] = -1.5; trMax[kZrez] = 1.5;
    trTitle[kBC]=StrDup("layer"); trNbins[kBC] = AliTRDgeometry::kNlayer; trMin[kBC] = -0.5; trMax[kBC] = AliTRDgeometry::kNlayer-0.5;
    Int_t jdim(kNdim);
    trTitle[jdim]=StrDup("centrality"); trNbins[jdim] = AliTRDeventInfo::kCentralityClasses; trMin[jdim] = -.5; trMax[jdim] = AliTRDeventInfo::kCentralityClasses - 0.5;
    jdim++; trTitle[jdim]=StrDup("occupancy [%]"); trNbins[jdim] = 12; trMin[jdim] = 25.; trMax[jdim] = 85.;
//    jdim++; trTitle[jdim]=StrDup("gap [cm]"); trNbins[jdim] = 80; trMin[jdim] = 0.1; trMax[jdim] = 2.1;
//     jdim++; trTitle[jdim]=StrDup("size_{0} [cm]"); trNbins[jdim] = 16; trMin[jdim] = 0.15; trMax[jdim] = 1.75;
//     jdim++; trTitle[jdim]=StrDup("pos_{0} [cm]"); trNbins[jdim] = 10; trMin[jdim] = 0.; trMax[jdim] = 3.5;
//     jdim++; trTitle[jdim]=StrDup("size_{1} [cm]"); trNbins[jdim] = 16; trMin[jdim] = 0.15; trMax[jdim] = 1.75;
//     jdim++; trTitle[jdim]=StrDup("pos_{1} [cm]"); trNbins[jdim] = 10; trMin[jdim] = 0.; trMax[jdim] = 3.5;

    st = "tracklet spatial&charge resolution;";
    // define minimum info to be saved in non debug mode
    Int_t ndim=DebugLevel()>=1?(jdim+1):kNdimTrklt;
    for(Int_t idim(0); idim<ndim; idim++){ st += trTitle[idim]; st+=";";}
    H = new THnSparseI(hn, st.Data(), ndim, trNbins, trMin, trMax);
  } else H->Reset();
  fContainer->AddAt(H, kTracklet);
  //++++++++++++++++++++++
  // tracklet to TRDin
  snprintf(hn, nhn, "h%s", fgPerformanceName[kTrackIn]);
  if(!(H = (THnSparseI*)gROOT->FindObject(hn))){
    // set specific fields
    const Int_t mdim(kNdim+3);
    Char_t *trinTitle[mdim]; memcpy(trinTitle, fgkTitle, kNdim*sizeof(Char_t*));
    Int_t trinNbins[mdim];   memcpy(trinNbins, fgkNbins, kNdim*sizeof(Int_t));
    Double_t trinMin[mdim];  memcpy(trinMin, fgkMin, kNdim*sizeof(Double_t));
    Double_t trinMax[mdim];  memcpy(trinMax, fgkMax, kNdim*sizeof(Double_t));
    trinNbins[kSpeciesChgRC] = Int_t(kNcharge)*AliPID::kSPECIES+1; trinMin[kSpeciesChgRC] = -AliPID::kSPECIES-0.5; trinMax[kSpeciesChgRC] = AliPID::kSPECIES+0.5;
    trinTitle[kNdim]=StrDup("p [GeV/c]"); trinNbins[kNdim] = 24; trinMin[kNdim] = -0.5; trinMax[kNdim] = 23.5;
    trinTitle[kNdim+1]=StrDup("dx [cm]"); trinNbins[kNdim+1]=48; trinMin[kNdim+1]=-2.4; trinMax[kNdim+1]=2.4;
    trinTitle[kNdim+2]=StrDup("Fill Bunch"); trinNbins[kNdim+2]=3500; trinMin[kNdim+2]=-0.5; trinMax[kNdim+2]=3499.5;
    st = "r-#phi/z/angular residuals @ TRD entry;";
    // define minimum info to be saved in non debug mode
    Int_t ndim=DebugLevel()>=1?mdim:kNdim;//kNdimTrkIn;
    for(Int_t idim(0); idim<ndim; idim++){st+=trinTitle[idim]; st+=";";}
    H = new THnSparseI(hn, st.Data(), ndim, trinNbins, trinMin, trinMax);
  } else H->Reset();
  fContainer->AddAt(H, kTrackIn);
  // tracklet to TRDout
//  fContainer->AddAt(BuildMonitorContainerTracklet("TrkOUT"), kTrackOut);


  // Resolution histos
  if(!HasMCdata()) return fContainer;

  //++++++++++++++++++++++
  // cluster to TrackRef residuals/pulls
  snprintf(hn, nhn, "h%s", fgPerformanceName[kMCcluster]);
  if(!(H = (THnSparseI*)gROOT->FindObject(hn))){
    const Char_t *clTitle[kNdim] = {"layer", fgkTitle[kPhi], fgkTitle[kEta], fgkTitle[kYrez], "#Deltax [cm]", "Q</Q", fgkTitle[kSpeciesChgRC], "#Phi [deg]"};
    const Int_t clNbins[kNdim]   = {AliTRDgeometry::kNlayer, fgkNbins[kPhi], fgkNbins[kEta], fgkNbins[kYrez], 20, 10, Int_t(kNcharge)*AliPID::kSPECIES+1, 15};
    const Double_t clMin[kNdim]  = {-0.5, fgkMin[kPhi], fgkMin[kEta], fgkMin[kYrez]/10., 0., 0.1, -AliPID::kSPECIES-0.5, -45},
                   clMax[kNdim]  = {AliTRDgeometry::kNlayer-0.5, fgkMax[kPhi], fgkMax[kEta], fgkMax[kYrez]/10., 4., 2.1, AliPID::kSPECIES+0.5, 45};
    st = "MC cluster spatial resolution;";
    // define minimum info to be saved in non debug mode
    Int_t ndim=DebugLevel()>=1?kNdim:4;
    for(Int_t idim(0); idim<ndim; idim++){ st += clTitle[idim]; st+=";";}
    H = new THnSparseI(hn, st.Data(), ndim, clNbins, clMin, clMax);
  } else H->Reset();
  fContainer->AddAt(H, kMCcluster);
  //++++++++++++++++++++++
  // tracklet to TrackRef
  snprintf(hn, nhn, "h%s", fgPerformanceName[kMCtracklet]);
  if(!(H = (THnSparseI*)gROOT->FindObject(hn))){
    Char_t *trTitle[kNdim]; memcpy(trTitle, fgkTitle, kNdim*sizeof(Char_t*));
    Int_t trNbins[kNdim]; memcpy(trNbins, fgkNbins, kNdim*sizeof(Int_t));
    Double_t trMin[kNdim]; memcpy(trMin, fgkMin, kNdim*sizeof(Double_t));
    Double_t trMax[kNdim]; memcpy(trMax, fgkMax, kNdim*sizeof(Double_t));
    // set specific fields
    trTitle[kBC]=StrDup("layer"); trNbins[kBC] = AliTRDgeometry::kNlayer; trMin[kBC] = -0.5; trMax[kBC] = AliTRDgeometry::kNlayer-0.5;
    trMin[kYrez] = -0.54; trMax[kYrez] = -trMin[kYrez];
    trMin[kPrez] = -4.5; trMax[kPrez] = -trMin[kPrez];
    trMin[kZrez] = -1.5; trMax[kZrez] = -trMin[kZrez];
    trNbins[kSpeciesChgRC] = Int_t(kNcharge)*AliPID::kSPECIES+1;trMin[kSpeciesChgRC] = -AliPID::kSPECIES-0.5; trMax[kSpeciesChgRC] = AliPID::kSPECIES+0.5;

    st = "MC tracklet spatial resolution;";
    // define minimum info to be saved in non debug mode
    Int_t ndim=DebugLevel()>=1?kNdim:4;
    for(Int_t idim(0); idim<ndim; idim++){ st += trTitle[idim]; st+=";";}
    H = new THnSparseI(hn, st.Data(), ndim, trNbins, trMin, trMax);
  } else H->Reset();
  fContainer->AddAt(H, kMCtracklet);
  //++++++++++++++++++++++
  // TRDin to TrackRef
  snprintf(hn, nhn, "h%s", fgPerformanceName[kMCtrackIn]);
  if(!(H = (THnSparseI*)gROOT->FindObject(hn))){
    st = "MC r-#phi/z/angular residuals @ TRD entry;";
    // set specific fields
    Int_t trNbins[kNdim]; memcpy(trNbins, fgkNbins, kNdim*sizeof(Int_t));
    Double_t trMin[kNdim]; memcpy(trMin, fgkMin, kNdim*sizeof(Double_t));
    Double_t trMax[kNdim]; memcpy(trMax, fgkMax, kNdim*sizeof(Double_t));
    trMin[kYrez] = -0.54; trMax[kYrez] = -trMin[kYrez];
    trMin[kPrez] = -2.4; trMax[kPrez] = -trMin[kPrez];
    trMin[kZrez] = -0.9; trMax[kZrez] = -trMin[kZrez];
    trNbins[kSpeciesChgRC] = Int_t(kNcharge)*AliPID::kSPECIES+1;trMin[kSpeciesChgRC] = -AliPID::kSPECIES-0.5; trMax[kSpeciesChgRC] = AliPID::kSPECIES+0.5;
    // define minimum info to be saved in non debug mode
    Int_t ndim=DebugLevel()>=1?kNdim:7;
    for(Int_t idim(0); idim<ndim; idim++){ st += fgkTitle[idim]; st+=";";}
    H = new THnSparseI(hn, st.Data(), ndim, trNbins, trMin, trMax);
  } else H->Reset();
  fContainer->AddAt(H, kMCtrackIn);
  //++++++++++++++++++++++
  // track to TrackRef
  snprintf(hn, nhn, "h%s", fgPerformanceName[kMCtrack]);
  if(!(H = (THnSparseI*)gROOT->FindObject(hn))){
    Char_t *trTitle[kNdim+1]; memcpy(trTitle, fgkTitle, kNdim*sizeof(Char_t*));
    Int_t trNbins[kNdim+1]; memcpy(trNbins, fgkNbins, kNdim*sizeof(Int_t));
    Double_t trMin[kNdim+1]; memcpy(trMin, fgkMin, kNdim*sizeof(Double_t));
    Double_t trMax[kNdim+1]; memcpy(trMax, fgkMax, kNdim*sizeof(Double_t));
    // set specific fields
    trTitle[kBC]=StrDup("layer"); trNbins[kBC] = AliTRDgeometry::kNlayer; trMin[kBC] = -0.5; trMax[kBC] = AliTRDgeometry::kNlayer-0.5;
    trMin[kYrez] = -0.9; trMax[kYrez] = -trMin[kYrez];
    trMin[kPrez] = -1.5; trMax[kPrez] = -trMin[kPrez];
    trMin[kZrez] = -0.9; trMax[kZrez] = -trMin[kZrez];
    trNbins[kSpeciesChgRC] = Int_t(kNcharge)*AliPID::kSPECIES+1;trMin[kSpeciesChgRC] = -AliPID::kSPECIES-0.5; trMax[kSpeciesChgRC] = AliPID::kSPECIES+0.5;
    trTitle[kNdim]=StrDup("#Deltap_{t}/p_{t} [%]"); trNbins[kNdim] = 25; trMin[kNdim] = -4.5; trMax[kNdim] = 20.5;

    st = "MC track spatial&p_{t} resolution;";
    // define minimum info to be saved in non debug mode
    Int_t ndim=DebugLevel()>=1?(kNdim+1):4;
    for(Int_t idim(0); idim<ndim; idim++){ st += trTitle[idim]; st+=";";}
    H = new THnSparseI(hn, st.Data(), ndim, trNbins, trMin, trMax);
  } else H->Reset();
  fContainer->AddAt(H, kMCtrack);

  return fContainer;
}

//____________________________________________________________________
Bool_t AliTRDresolution::FitTrack(const Int_t np, AliTrackPoint *points, Float_t param[10])
{
//
// Fit track with a staight line using the "np" clusters stored in the array "points".
// The following particularities are stored in the clusters from points:
//   1. pad tilt as cluster charge
//   2. pad row cross or vertex constrain as fake cluster with cluster type 1
// The parameters of the straight line fit are stored in the array "param" in the following order :
//     param[0] - x0 reference radial position
//     param[1] - y0 reference r-phi position @ x0
//     param[2] - z0 reference z position @ x0
//     param[3] - slope dy/dx
//     param[4] - slope dz/dx
//
// Attention :
// Function should be used to refit tracks for B=0T
//

  if(np<40){
    if(AliLog::GetDebugLevel("PWGPP", "AliTRDresolution")>1) printf("D-AliTRDresolution::FitTrack: Not enough clusters to fit a track [%d].\n", np);
    return kFALSE;
  }
  TLinearFitter yfitter(2, "pol1"), zfitter(2, "pol1");

  Double_t x0(0.);
  for(Int_t ip(0); ip<np; ip++) x0+=points[ip].GetX();
  x0/=Float_t(np);

  Double_t x, y, z, dx, tilt(0.);
  for(Int_t ip(0); ip<np; ip++){
    x = points[ip].GetX(); z = points[ip].GetZ();
    dx = x - x0;
    zfitter.AddPoint(&dx, z, points[ip].GetClusterType()?1.e-3:1.);
  }
  if(zfitter.Eval() != 0) return kFALSE;

  Double_t z0    = zfitter.GetParameter(0);
  Double_t dzdx  = zfitter.GetParameter(1);
  for(Int_t ip(0); ip<np; ip++){
    if(points[ip].GetClusterType()) continue;
    x    = points[ip].GetX();
    dx   = x - x0;
    y    = points[ip].GetY();
    z    = points[ip].GetZ();
    tilt = points[ip].GetCharge();
    y -= tilt*(-dzdx*dx + z - z0);
    Float_t xyz[3] = {x, y, z}; points[ip].SetXYZ(xyz);
    yfitter.AddPoint(&dx, y, 1.);
  }
  if(yfitter.Eval() != 0) return kFALSE;
  Double_t y0   = yfitter.GetParameter(0);
  Double_t dydx = yfitter.GetParameter(1);

  param[0] = x0; param[1] = y0; param[2] = z0; param[3] = dydx; param[4] = dzdx;
  if(AliLog::GetDebugLevel("PWGPP", "AliTRDresolution")>3) printf("D-AliTRDresolution::FitTrack: x0[%f] y0[%f] z0[%f] dydx[%f] dzdx[%f].\n", x0, y0, z0, dydx, dzdx);
  return kTRUE;
}

//____________________________________________________________________
Bool_t AliTRDresolution::FitTracklet(const Int_t ly, const Int_t np, const AliTrackPoint *points, const Float_t param[10], Float_t par[3])
{
//
// Fit tracklet with a staight line using the coresponding subset of clusters out of the total "np" clusters stored in the array "points".
// See function FitTrack for the data stored in the "clusters" array

// The parameters of the straight line fit are stored in the array "param" in the following order :
//     par[0] - x0 reference radial position
//     par[1] - y0 reference r-phi position @ x0
//     par[2] - slope dy/dx
//
// Attention :
// Function should be used to refit tracks for B=0T
//

  TLinearFitter yfitter(2, "pol1");

  // grep data for tracklet
  Double_t x0(0.), x[60], y[60], dy[60];
  Int_t nly(0);
  for(Int_t ip(0); ip<np; ip++){
    if(points[ip].GetClusterType()) continue;
    if(points[ip].GetVolumeID() != ly) continue;
    Float_t xt(points[ip].GetX())
           ,yt(param[1] + param[3] * (xt - param[0]));
    x[nly] = xt;
    y[nly] = points[ip].GetY();
    dy[nly]= y[nly]-yt;
    x0    += xt;
    nly++;
  }
  if(nly<10){
    if(AliLog::GetDebugLevel("PWGPP", "AliTRDresolution")>1) printf("D-AliTRDresolution::FitTracklet: Not enough clusters to fit a tracklet [%d].\n", nly);
    return kFALSE;
  }
  // set radial reference for fit
  x0 /= Float_t(nly);

  // find tracklet core
  Double_t mean(0.), sig(1.e3);
  AliMathBase::EvaluateUni(nly, dy, mean, sig, 0);

  // simple cluster error parameterization
  Float_t kSigCut = TMath::Sqrt(5.e-4 + param[3]*param[3]*0.018);

  // fit tracklet core
  for(Int_t jly(0); jly<nly; jly++){
    if(TMath::Abs(dy[jly]-mean)>kSigCut) continue;
    Double_t dx(x[jly]-x0);
    yfitter.AddPoint(&dx, y[jly], 1.);
  }
  if(yfitter.Eval() != 0) return kFALSE;
  par[0] = x0;
  par[1] = yfitter.GetParameter(0);
  par[2] = yfitter.GetParameter(1);
  return kTRUE;
}

//____________________________________________________________________
Bool_t AliTRDresolution::UseTrack(const Int_t np, const AliTrackPoint *points, Float_t param[10])
{
//
// Global selection mechanism of tracksbased on cluster to fit residuals
// The parameters are the same as used ni function FitTrack().

  const Float_t kS(0.6), kM(0.2);
  TH1S h("h1", "", 100, -5.*kS, 5.*kS);
  Float_t dy, dz, s, m;
  for(Int_t ip(0); ip<np; ip++){
    if(points[ip].GetClusterType()) continue;
    Float_t x0(points[ip].GetX())
          ,y0(param[1] + param[3] * (x0 - param[0]))
          ,z0(param[2] + param[4] * (x0 - param[0]));
    dy=points[ip].GetY() - y0; h.Fill(dy);
    dz=points[ip].GetZ() - z0;
  }
  TF1 fg("fg", "gaus", -5.*kS, 5.*kS);
  fg.SetParameter(1, 0.);
  fg.SetParameter(2, 2.e-2);
  h.Fit(&fg, "QN");
  m=fg.GetParameter(1); s=fg.GetParameter(2);
  if(s>kS || TMath::Abs(m)>kM) return kFALSE;
  return kTRUE;
}

//________________________________________________________
void AliTRDresolution::GetLandauMpvFwhm(TF1 * const f, Float_t &mpv, Float_t &xm, Float_t &xM)
{
  //
  // Get the most probable value and the full width half mean
  // of a Landau distribution
  //

  const Float_t dx = 1.;
  mpv = f->GetParameter(1);
  Float_t fx, max = f->Eval(mpv);

  xm = mpv - dx;
  while((fx = f->Eval(xm))>.5*max){
    if(fx>max){
      max = fx;
      mpv = xm;
    }
    xm -= dx;
  }

  xM += 2*(mpv - xm);
  while((fx = f->Eval(xM))>.5*max) xM += dx;
}


// #include "TFile.h"
// //________________________________________________________
// Bool_t AliTRDresolution::LoadCorrection(const Char_t *file)
// {
//   if(!file){
//     AliWarning("Use cluster position as in reconstruction.");
//     SetLoadCorrection();
//     return kTRUE;
//   }
//   TDirectory *cwd(gDirectory);
//   TString fileList;
//   FILE *filePtr = fopen(file, "rt");
//   if(!filePtr){
//     AliWarning(Form("Couldn't open correction list \"%s\". Use cluster position as in reconstruction.", file));
//     SetLoadCorrection();
//     return kFALSE;
//   }
//   TH2 *h2 = new TH2F("h2", ";time [time bins];detector;dx [#mum]", 30, -0.5, 29.5, AliTRDgeometry::kNdet, -0.5, AliTRDgeometry::kNdet-0.5);
//   while(fileList.Gets(filePtr)){
//     if(!TFile::Open(fileList.Data())) {
//       AliWarning(Form("Couldn't open \"%s\"", fileList.Data()));
//       continue;
//     } else AliInfo(Form("\"%s\"", fileList.Data()));
//
//     TTree *tSys = (TTree*)gFile->Get("tSys");
//     h2->SetDirectory(gDirectory); h2->Reset("ICE");
//     tSys->Draw("det:t>>h2", "dx", "goff");
//     for(Int_t idet(0); idet<AliTRDgeometry::kNdet; idet++){
//       for(Int_t it(0); it<30; it++) fXcorr[idet][it]+=(1.e-4*h2->GetBinContent(it+1, idet+1));
//     }
//     h2->SetDirectory(cwd);
//     gFile->Close();
//   }
//   cwd->cd();
//
//   if(AliLog::GetDebugLevel("PWGPP", "AliTRDresolution")>=2){
//     for(Int_t il(0); il<184; il++) printf("-"); printf("\n");
//     printf("DET|");for(Int_t it(0); it<30; it++) printf(" tb%02d|", it); printf("\n");
//     for(Int_t il(0); il<184; il++) printf("-"); printf("\n");
//     FILE *fout = fopen("TRD.ClusterCorrection.txt", "at");
//     fprintf(fout, "  static const Double_t dx[AliTRDgeometry::kNdet][30]={\n");
//     for(Int_t idet(0); idet<AliTRDgeometry::kNdet; idet++){
//       printf("%03d|", idet);
//       fprintf(fout, "    {");
//       for(Int_t it(0); it<30; it++){
//         printf("%+5.0f|", 1.e4*fXcorr[idet][it]);
//         fprintf(fout, "%+6.4f%c", fXcorr[idet][it], it==29?' ':',');
//       }
//       printf("\n");
//       fprintf(fout, "}%c\n", idet==AliTRDgeometry::kNdet-1?' ':',');
//     }
//     fprintf(fout, "  };\n");
//   }
//   SetLoadCorrection();
//   return kTRUE;
// }

//________________________________________________________
AliTRDresolution::AliTRDresolutionProjection::AliTRDresolutionProjection()
  :TNamed()
  ,fH(NULL)
  ,fNrebin(0)
  ,fRebinX(NULL)
  ,fRebinY(NULL)
{
  // constructor
  memset(fAx, 0, 3*sizeof(Int_t));
  memset(fRange, 0, 4*sizeof(Float_t));
}

//________________________________________________________
AliTRDresolution::AliTRDresolutionProjection::~AliTRDresolutionProjection()
{
  // destructor
  if(fH) delete fH;
}

//________________________________________________________
void AliTRDresolution::AliTRDresolutionProjection::Build(const Char_t *n, const Char_t *t, Int_t ix, Int_t iy, Int_t iz, TAxis *aa[])
{
// check and build (if neccessary) projection determined by axis "ix", "iy" and "iz"
  if(!aa[ix] || !aa[iy] || !aa[iz]) return;
  TAxis *ax(aa[ix]), *ay(aa[iy]), *az(aa[iz]);
  // check ax definiton to protect against older versions of the data
  if(ax->GetNbins()<=0 || (ax->GetXmax()-ax->GetXmin())<=0.){
    AliWarning(Form("Wrong definition of axis[%d] \"%s\"[%d](%f %f).", ix, ax->GetTitle(), ax->GetNbins(), ax->GetXmin(), ax->GetXmax()));
    return;
  }
  if(ay->GetNbins()<=0 || (ay->GetXmax()-ay->GetXmin())<=0.){
    AliWarning(Form("Wrong definition of axis[%d] \"%s\"[%d](%f %f).", ix, ay->GetTitle(), ay->GetNbins(), ay->GetXmin(), ay->GetXmax()));
    return;
  }
  if(az->GetNbins()<=0 || (az->GetXmax()-az->GetXmin())<=0.){
    AliWarning(Form("Wrong definition of axis[%d] \"%s\"[%d](%f %f).", ix, az->GetTitle(), az->GetNbins(), az->GetXmin(), az->GetXmax()));
    return;
  }
  SetNameTitle(n,t);
  fH = new TH3I(n, Form("%s;%s;%s;%s", t, ax->GetTitle(), ay->GetTitle(), az->GetTitle()),
    ax->GetNbins(), ax->GetXmin(), ax->GetXmax(),
    ay->GetNbins(), ay->GetXmin(), ay->GetXmax(),
    az->GetNbins(), az->GetXmin(), az->GetXmax());
  fAx[0] = ix; fAx[1] = iy; fAx[2] = iz;
  fRange[0] = az->GetXmin()/3.; fRange[1] = az->GetXmax()/3.;
  AliDebug(2, Form("H3(%s, %s) :: %s[%3d %4.2f %4.2f]%s[%3d %4.2f %4.2f]%s[%3d %4.2f %4.2f]", n, t,
    ax->GetTitle(), ax->GetNbins(), ax->GetXmin(), ax->GetXmax(),
    ay->GetTitle(), ay->GetNbins(), ay->GetXmin(), ay->GetXmax(),
    az->GetTitle(), az->GetNbins(), az->GetXmin(), az->GetXmax()));
}

//________________________________________________________
AliTRDresolution::AliTRDresolutionProjection& AliTRDresolution::AliTRDresolutionProjection::operator=(const AliTRDresolutionProjection& rhs)
{
// copy projections
  if(this == &rhs) return *this;

  TNamed::operator=(rhs);
  if(fNrebin){fNrebin=0; delete [] fRebinX; delete [] fRebinY;}
  if(rhs.fNrebin) SetRebinStrategy(rhs.fNrebin, rhs.fRebinX, rhs.fRebinY);
  memcpy(fAx, rhs.fAx, 3*sizeof(Int_t));
  memcpy(fRange, rhs.fRange, 4*sizeof(Float_t));
  if(fH) delete fH;
  if(rhs.fH) fH=(TH3I*)rhs.fH->Clone(Form("%s_CLONE", rhs.fH->GetName()));
  return *this;
}

//________________________________________________________
AliTRDresolution::AliTRDresolutionProjection& AliTRDresolution::AliTRDresolutionProjection::operator+=(const AliTRDresolutionProjection& other)
{
// increment projections
  if(!fH || !other.fH) return *this;
  AliDebug(2, Form("%s+=%s [%s+=%s]", GetName(), other.GetName(), fH->GetName(), (other.fH)->GetName()));
  fH->Add(other.fH);
  return *this;
}

//________________________________________________________
void AliTRDresolution::AliTRDresolutionProjection::Increment(Int_t bin[], Double_t v)
{
// increment bin with value "v" pointed by general coord in "bin"
  if(!fH) return;
  AliDebug(4, Form("  %s[%2d]", fH->GetName(), Int_t(v)));
  fH->AddBinContent(fH->GetBin(bin[fAx[0]],bin[fAx[1]],bin[fAx[2]]), Int_t(v));
}

//________________________________________________________
TH2* AliTRDresolution::AliTRDresolutionProjection::Projection2D(const Int_t nstat, const Int_t ncol, const Int_t mid, Bool_t del)
{
// build the 2D projection and adjust binning

  const Char_t *title[] = {"Mean", "#mu", "MPV"};
  if(!fH) return NULL;
  TAxis *ax(fH->GetXaxis()), *ay(fH->GetYaxis()), *az(fH->GetZaxis());
  TH2D *h2s(NULL), *hyx(NULL);
  if(!(h2s = (TH2D*)fH->Project3D("yx"))) return NULL;
  // save a copy of the original distribution
  if(!del) hyx = (TH2D*)h2s->Clone(Form("%sEn", fH->GetName()));
  Int_t irebin(0), dxBin(1), dyBin(1);
  while(irebin<fNrebin && (AliTRDresolution::GetMeanStat(h2s, .5, ">")<nstat)){
    h2s->Rebin2D(fRebinX[irebin], fRebinY[irebin]);
    dxBin*=fRebinX[irebin];dyBin*=fRebinY[irebin];
    irebin++;
  }
  Int_t nx(h2s->GetNbinsX()), ny(h2s->GetNbinsY());
  delete h2s;

  // start projection
  TH1 *h(NULL); Int_t n(0);
  Float_t dz=(fRange[1]-fRange[1])/ncol;
  TString titlez(az->GetTitle()); TObjArray *tokenTitle(titlez.Tokenize(" "));
  Int_t nt(tokenTitle->GetEntriesFast());
  TH2 *h2 = new TH2F(Form("%s_2D", fH->GetName()),
            Form("%s;%s;%s;%s(%s) %s", fH->GetTitle(), ax->GetTitle(), ay->GetTitle(), title[mid], nt>0?(*tokenTitle)[0]->GetName():"", nt>1?(*tokenTitle)[1]->GetName():""),
            nx, ax->GetXmin(), ax->GetXmax(), ny, ay->GetXmin(), ay->GetXmax());
  h2->SetContour(ncol);
  h2->GetZaxis()->CenterTitle();
  h2->GetZaxis()->SetTitleOffset(1.4);
  h2->GetZaxis()->SetRangeUser(fRange[0], fRange[1]);
  AliDebug(2, Form("%s[%s] nx[%d] ny[%d]", h2->GetName(), h2->GetTitle(), nx, ny));
  for(Int_t iy(0); iy<ny; iy++){
    for(Int_t ix(0); ix<nx; ix++){
      h = fH->ProjectionZ(Form("%s_z", h2->GetName()), ix*dxBin+1, (ix+1)*dxBin+1, iy*dyBin+1, (iy+1)*dyBin+1);
      Int_t ne((Int_t)h->Integral());
      if(ne<nstat/2){
        h2->SetBinContent(ix+1, iy+1, -999);
        h2->SetBinError(ix+1, iy+1, 1.);
        n++;
      }else{
        Float_t v(h->GetMean()), ve(h->GetRMS());
        if(mid==1){
          TF1 fg("fg", "gaus", az->GetXmin(), az->GetXmax());
          fg.SetParameter(0, Float_t(ne)); fg.SetParameter(1, v); fg.SetParameter(2, ve);
          h->Fit(&fg, "WQ");
          v = fg.GetParameter(1); ve = fg.GetParameter(2);
        } else if (mid==2) {
          TF1 fl("fl", "landau", az->GetXmin(), az->GetXmax());
          fl.SetParameter(0, Float_t(ne)); fl.SetParameter(1, v); fl.SetParameter(2, ve);
          h->Fit(&fl, "WQ");
          v = fl.GetMaximumX(); ve = fl.GetParameter(2);
/*          TF1 fgle("gle", "[0]*TMath::Landau(x, [1], [2], 1)*TMath::Exp(-[3]*x/[1])", az->GetXmin(), az->GetXmax());
          fgle.SetParameter(0, fl.GetParameter(0));
          fgle.SetParameter(1, fl.GetParameter(1));
          fgle.SetParameter(2, fl.GetParameter(2));
          fgle.SetParameter(3, 1.);fgle.SetParLimits(3, 0., 5.);
          h->Fit(&fgle, "WQ");
          v = fgle.GetMaximumX(); ve = fgle.GetParameter(2);*/
        }
        if(v<fRange[0]) h2->SetBinContent(ix+1, iy+1, fRange[0]+0.1*dz);
        else h2->SetBinContent(ix+1, iy+1, v);
        h2->SetBinError(ix+1, iy+1, ve);
      }
    }
  }
  if(h) delete h;
  if(n==nx*ny){delete h2; h2=NULL;} // clean empty projections
  return h2;
}

//________________________________________________________
void AliTRDresolution::SetNormZ(TH2 *h2, Int_t bxmin, Int_t bxmax, Int_t bymin, Int_t bymax, Float_t thr)
{
// Normalize histo content to the mean value in the range specified by bin ranges
// [bxmin, bxmax] on the x axis and [bymin, bymax] on the y axis.
// Optionally a threshold "thr" can be specified to disregard entries with no meaning 

  Float_t s = 0., c=0.; Int_t is(0);
  for(Int_t ix(bxmin); ix<=(bxmax>0?bxmax:(h2->GetXaxis()->GetNbins())); ix++){
    for(Int_t iy(bymin); iy<=(bymax>0?bymax:(h2->GetYaxis()->GetNbins())); iy++){
      if((c = h2->GetBinContent(ix, iy))<thr) continue;
      s += c; is++;
    }
  }
  s/= is?is:1;
  for(Int_t ix(1); ix<=h2->GetXaxis()->GetNbins(); ix++){
    for(Int_t iy(1); iy<=h2->GetYaxis()->GetNbins(); iy++){
      if((c = h2->GetBinContent(ix, iy))<thr) h2->SetBinContent(ix, iy, thr-100);
      else h2->SetBinContent(ix, iy, 100.*(c/s-1.));
    }
  }
}

//________________________________________________________
void AliTRDresolution::SetProcesses(Bool_t det, Bool_t cl, Bool_t trklt, Bool_t trkin)
{
// steer processes
  if(det) SETBIT(fSteer, kDetector); else CLRBIT(fSteer, kDetector);
  if(cl)  SETBIT(fSteer, kCluster); else CLRBIT(fSteer, kCluster);
  if(trklt) SETBIT(fSteer, kTracklet); else CLRBIT(fSteer, kTracklet);
  if(trkin) SETBIT(fSteer, kTrackIn); else CLRBIT(fSteer, kTrackIn);
}

//________________________________________________________
void AliTRDresolution::SetRangeZ(TH2 *h2, Float_t min, Float_t max, Float_t thr)
{
// Set range on Z axis such to avoid outliers

  Float_t c(0.), dz(1.e-3*(max-min));
  for(Int_t ix(1); ix<=h2->GetXaxis()->GetNbins(); ix++){
    for(Int_t iy(1); iy<=h2->GetYaxis()->GetNbins(); iy++){
      if((c = h2->GetBinContent(ix, iy))<thr) continue;
      if(c<=min) h2->SetBinContent(ix, iy, min+dz);
    }
  }
  h2->GetZaxis()->SetRangeUser(min, max);
}

//________________________________________________________
void AliTRDresolution::AliTRDresolutionProjection::SetRebinStrategy(Int_t n, Int_t rebx[], Int_t reby[])
{
// define rebinning strategy for this projection
  fNrebin = n;
  fRebinX = new Int_t[n]; memcpy(fRebinX, rebx, n*sizeof(Int_t));
  fRebinY = new Int_t[n]; memcpy(fRebinY, reby, n*sizeof(Int_t));
}


