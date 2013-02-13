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

#include <TFile.h>
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
#include "AliESDHeader.h"
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
ClassImp(AliTRDresolution::AliTRDrecoProjection)

Int_t const   AliTRDresolution::fgkNbins[kNdim] = {
  Int_t(kNbunchCross)/*bc*/,
  144/*phi*/,
  45/*eta*/,
  50/*dy*/,
  40/*dphi*/,
  50/*dz*/,
  5/*chg*species*/,
  3/*pt*/
};  //! no of bins/projection
Double_t const AliTRDresolution::fgkMin[kNdim] = {
  -1.5,        /*bc*/
  -TMath::Pi(),/*phi*/
  -0.9,        /*eta*/
  -1.5,         /*dy*/
  -5.,         /*dphi*/
  -2.5,        /*dz*/
  -2.5,        /*chg*species*/
  -0.5         /*pt*/
};    //! low limits for projections
Double_t const AliTRDresolution::fgkMax[kNdim] = {
  1.5,        /*bc*/
  TMath::Pi(),/*phi*/
  0.9,        /*eta*/
  1.5,         /*dy*/
  5.,         /*dphi*/
  2.5,        /*dz*/
  2.5,        /*chg*species*/
  2.5         /*pt*/
};    //! high limits for projections
Char_t const *AliTRDresolution::fgkTitle[kNdim] = {
  "bunch cross",
  "#phi [rad]",
  "#eta",
  "#Deltay [cm]",
  "#Delta#phi [deg]",
  "#Deltaz [cm]",
  "chg*spec*rc",
  "p_{t} [bin]"
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

//________________________________________________________
AliTRDresolution::AliTRDresolution()
  :AliTRDrecoTask()
  ,fSteer(0)
  ,fPtThreshold(.3)
  ,fBCbinTOF(0)
  ,fBCbinFill(0)
  ,fLYselect(-1)
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
  SetProcesses(kTRUE, kTRUE, kTRUE, kTRUE);
}

//________________________________________________________
AliTRDresolution::AliTRDresolution(char* name, Bool_t xchange)
  :AliTRDrecoTask(name, "TRD spatial and momentum resolution")
  ,fSteer(0)
  ,fPtThreshold(.3)
  ,fBCbinTOF(0)
  ,fBCbinFill(0)
  ,fLYselect(-1)
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
  if(fProj){
    fProj->Delete();
    delete fProj;
  }
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
    if(!(cl = (AliTRDcluster*)(*fkClusters)[icl])){
      AliDebug(2, Form("Missing clusters @ %d", icl));
      continue;
    }
    if(cl->IsMasked()) printf("Mask %03d[%02d_%d_%d] %d[%d]", det, AliTRDgeometry::GetSector(det), stk, ly, cl->GetPadMaskedPosition(), cl->GetPadMaskedStatus());
    val[kBC]  = ly;
    val[kPhi] = TMath::ATan2(cl->GetX()*sn + cl->GetY()*cs, cl->GetX()*cs - cl->GetY()*sn);
    val[kEta] = (5-stk)*16-cl->GetPadRow()-1-(stk<3?4:0);
    val[kYrez]= fEvent->GetMultiplicity();
    val[4]    = TMath::Abs(cl->GetQ());
    val[5]    = fTriggerList?fTriggerSlot:(cl->IsFivePad()?1:cl->GetNPads());
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
      val[kNdim]   = fEvent?fEvent->GetMultiplicity():0.;
      val[kNdim+1] = fTriggerList?fTriggerSlot:(c->IsFivePad()?1:c->GetNPads());
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

  // down scale PID resolution
  Int_t spc(fSpecies); if(spc==3) spc=2; if(spc==-3) spc=-2; if(spc>3) spc=3; if(spc<-3) spc=-3;
  const Int_t ndim(kNdim+8);
  Double_t val[ndim],
           alpha(0.), cs(-2.), sn(0.);
//  Float_t sz[AliTRDseedV1::kNtb], pos[AliTRDseedV1::kNtb]; Int_t ntbGap[AliTRDseedV1::kNtb];
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

    val[kSpeciesChgRC]= fTracklet->IsRowCross()?0:spc;// fSpecies;
    val[kPt]  = GetPtBin(fPt); //fPt<0.8?0:(fPt<1.5?1:2);//GetPtBin(fTracklet->GetMomentum());
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
    val[kNdim+1] = fTriggerList?fTriggerSlot:(1.e2*fTracklet->GetTBoccupancy()/AliTRDseedV1::kNtb);
/*    Int_t n = fTracklet->GetChargeGaps(sz, pos, ntbGap);
    val[kNdim+2] = 0.; for(Int_t igap(0); igap<n; igap++) val[kNdim+2] += sz[igap];
    for(Int_t ifill(0); ifill<3; ifill++){ val[kNdim+3+ifill]=0.;val[kNdim+4+ifill]=0.;}
    for(Int_t igap(0), ifill(0); igap<n&&ifill<2; igap++){
      if(ntbGap[igap]<2) continue;
      val[kNdim+3+ifill] = sz[igap];
      val[kNdim+4+ifill] = pos[igap];
      ifill++;
    }*/
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

  // down scale PID resolution
  Int_t spc(fSpecies); if(spc==3) spc=2; if(spc==-3) spc=-2; if(spc>3) spc=3; if(spc<-3) spc=-3;
  // read V0 PID if available
  Int_t v0pid(-2);
  if(fkESD->IsElectron()) v0pid = -1;
  else if(fkESD->IsPion()) v0pid = 0;
  else if(fkESD->IsProton()) v0pid = 1;
  if(DebugLevel()>=1/* && v0pid>-2*/){
    Float_t tpc(fkESD->GetTPCdedx());
    Float_t tof(fkESD->GetTOFbeta());
    Int_t ev(fEvent->GetEventHeader()->GetEventNumberInFile());
    AliTRDtrackV1 t(*fkTrack); t.SetOwner();
    (*DebugStream()) << "trackIn"
      <<"ev="       << ev
      <<"pid="      << v0pid
      <<"tpc="      << tpc
      <<"tof="      << tof
      <<"track.="   << &t
      <<"trackIn.=" << tin
      << "\n";
  }

  //Int_t bc(fkESD?fkESD->GetTOFbc()/2:0);
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

  Double_t val[kNdim+2];
  val[kBC]          = fTriggerSlot?fTriggerSlot:v0pid;//bc==0?0:(bc<0?-1.:1.);
  Double_t alpha = (0.5+AliTRDgeometry::GetSector(fTracklet->GetDetector()))*AliTRDgeometry::GetAlpha(),
           cs    = TMath::Cos(alpha),
           sn    = TMath::Sin(alpha);
  val[kPhi]         = TMath::ATan2(fTracklet->GetX()*sn + fTracklet->GetY()*cs, fTracklet->GetX()*cs - fTracklet->GetY()*sn);
  Float_t tgl       = fTracklet->GetZ()/fTracklet->GetX()/TMath::Sqrt(1.+fTracklet->GetY()*fTracklet->GetY()/fTracklet->GetX()/fTracklet->GetX());
  val[kEta]         = -TMath::Log(TMath::Tan(0.5 *  (0.5*TMath::Pi() - TMath::ATan(tgl))));
  val[kYrez]        = dy;
  val[kZrez]        = fTracklet->IsRowCross()?dz:(fTracklet->GetdQdl()*5.e-4 - 2.5);
  val[kPrez]        = dphi*TMath::RadToDeg();
  val[kSpeciesChgRC]= fTracklet->IsRowCross()?0:spc;
  Float_t pz(0.); if(fTracklet->GetMomentum()-fPt>1.e-5) pz = TMath::Sqrt((fTracklet->GetMomentum()-fPt)*(fTracklet->GetMomentum()+fPt));
  val[kPt]          = fTracklet->IsRowCross()?GetPtBin(pz):GetPtBin(fPt);
  val[kNdim]        = GetPtBin(fTracklet->GetMomentum());
  val[kNdim+1]      = dx;
  //val[kNdim+2]      = fEvent?fEvent->GetBunchFill():0;
  H->Fill(val);

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
//  Int_t bc(TMath::Abs(fkESD->GetTOFbc()));

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
  if(sIdx==3) sIdx=2; if(sIdx>3) sIdx=3; // downscale PID

  TH1 *h(NULL);
  AliTRDgeometry *geo(AliTRDinfoGen::Geometry());
  AliTRDseedV1 *fTracklet(NULL); TObjArray *clInfoArr(NULL);
  UChar_t s;
  Double_t x, y, z, pt, dydx, dzdx/*, dzdl*/;
  Float_t pt0, p0, x0, y0, z0, dx, dy, dz, dydx0, dzdx0;
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
  Double_t val[kNdim+2];
  for(Int_t ily=0; ily<AliTRDgeometry::kNlayer; ily++){
    if(!(fTracklet = fkTrack->GetTracklet(ily)))/* ||
       !fTracklet->IsOK())*/ continue;

    x= x0 = fTracklet->GetX();
    Bool_t rc(fTracklet->IsRowCross()); Float_t eta, phi;
    if(!fkMC->GetDirections(x0, y0, z0, dydx0, dzdx0, pt0, p0, eta, phi, s)) continue;

    // MC track position at reference radial position
    dx  = x0 - x;
    Float_t ymc = y0 - dx*dydx0;
    Float_t zmc = z0 - dx*dzdx0;
    //phi -= TMath::Pi();

    val[kBC]  = ily;
    val[kPhi] = phi;
    val[kEta] = eta;
    val[kSpeciesChgRC]= rc?0.:sign*sIdx;
    val[kPt]  = GetPtBin(pt0);
    val[kNdim]= GetPtBin(p0);
    Double_t tilt(fTracklet->GetTilt());
//             ,t2(tilt*tilt)
//             ,corr(1./(1. + t2))
//             ,cost(TMath::Sqrt(corr));

    AliExternalTrackParam *tin(fkTrack->GetTrackIn());
    if(ily==0 && tin){ // trackIn residuals
      // check radial position
      if(TMath::Abs(tin->GetX()-x)>1.e-3) AliDebug(1, Form("TrackIn radial mismatch. dx[cm]=%+4.1f", tin->GetX()-x));
      else{
        Int_t exactPID = -2;
        switch(TMath::Abs(pdg)){
        case 11: exactPID = -1;break;
        case 211: exactPID = 0;break;
        case 2212: exactPID = 1;break;
        }
        val[kBC]          = exactPID;
        val[kYrez]        = tin->GetY()-ymc;
        val[kZrez]        = rc?(tin->GetZ()-zmc):(fTracklet->GetdQdl()*5e-4 - 2.5);
        val[kPrez]        = (TMath::ASin(tin->GetSnp())-TMath::ATan(dydx0))*TMath::RadToDeg();
        val[kNdim+1]      = 0.;//dx;
        if((H = (THnSparseI*)fContainer->At(kMCtrackIn))) H->Fill(val);
        val[kBC]          = ily; // reset for subsequent components

        if(DebugLevel()>=1 && exactPID>-2){
          Float_t tpc(fkESD->GetTPCdedx());
          Float_t tof(fkESD->GetTOFbeta());
          AliTRDtrackV1 t(*fkTrack); t.SetOwner();
          (*DebugStream()) << "MCtrackIn"
            <<"pid="      << exactPID
            <<"tpc="      << tpc
            <<"tof="      << tof
            <<"track.="   << &t
            <<"trackIn.=" << tin
            << "\n";
        }
      }
    }
    //if(bc>1) break; // do nothing for the rest of TRD objects if satellite bunch

    // track residuals
    dydx = fTracklet->GetYref(1);
    dzdx = fTracklet->GetZref(1);
    //dzdl = fTracklet->GetTgl();
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
    val[kNdim] = fEvent?fEvent->GetMultiplicity():0;
    val[kNdim+1] = 1.e2*fTracklet->GetTBoccupancy()/AliTRDseedV1::kNtb;
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
    Float_t padCorr(tilt*fTracklet->GetPadLength()),
            corr(1./TMath::Sqrt(1.+dydx0*dydx0+dzdx0*dzdx0));
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
      val[kSpeciesChgRC]= rc?0.:(TMath::Max(q*corr, Float_t(3.)));
      val[kNdim+1] = c->IsFivePad()?1:c->GetNPads();
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
  if(!track) return NULL;
  // special care for EVE usage
  if(H && (h = (TH1*)gDirectory->Get(Form("%s_proj_%d", H->GetName(), kYrez)))) delete h;
  return H?H->Projection(kYrez):NULL;
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
void AliTRDresolution::MakeSummary()
{
// Build summary plots

  if(!fProj){
    AliError("Missing results");
    return;
  }
  TVirtualPad *p(NULL); TCanvas *cOut(NULL);
  TObjArray *arr(NULL); TH2 *h2(NULL);
  TH2 *h2e[100] = {NULL}; Int_t ih2e(0); // save sigma histos for later deletion

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

  if((arr = (TObjArray*)fProj->FindObject("hDet2Cluster"))){
    cOut = new TCanvas(Form("%s_DetOccupancy", GetName()), "Detector performance", 2*nx, 2*ny);
    cOut->Divide(AliTRDgeometry::kNlayer,AliTRDeventInfo::kCentralityClasses, 1.e-5, 1.e-5);
    Int_t n=0;
    for(Int_t icen(0); icen<AliTRDeventInfo::kCentralityClasses; icen++){
      for(Int_t ily(0); ily<AliTRDgeometry::kNlayer; ily++){
        p=cOut->cd(icen*AliTRDgeometry::kNlayer+ily+1); p->SetRightMargin(0.1572581);p->SetTopMargin(0.08262712);
        if(!(h2 = (TH2*)arr->FindObject(Form("HDet%d%dEn", ily, icen)))) continue;
        if(SetNormZ(h2, 1, -1, 1, -1, 10.) < 1.e3) continue;  // cut all bins with < 10 entries
        SetRangeZ(h2, -90., 90, -200.);
        //h2->GetXaxis()->SetNdivisions(1010);
        h2->GetZaxis()->SetTitle("Rel. Det. Occup. [%]");
        h2->GetZaxis()->CenterTitle();
        h2->SetContour(9); h2->Draw("colz"); n++;
        MakeDetectorPlot(ily, "p");
      }
    }
    if(n>=AliTRDgeometry::kNlayer) cOut->SaveAs(Form("%s.gif", cOut->GetName()));

    cOut = new TCanvas(Form("%s_DetCharge", GetName()), "Detector performance", nx, ny);
    cOut->Divide(3,2, 1.e-5, 1.e-5);
    for(Int_t ily(0); ily<AliTRDgeometry::kNlayer; ily++){
      p=cOut->cd(ily+1); p->SetRightMargin(0.1572581);p->SetTopMargin(0.08262712);
      if(!(h2 = (TH2*)arr->FindObject(Form("HDet%d_2D", UseLYselectTrklt()?fLYselect:ily)))) continue;
      SetNormZ(h2, 1, -1, 1, -1, 10.); // cut all <q> < 10
      SetRangeZ(h2, -30., 30., -200.);
      //h2->GetXaxis()->SetNdivisions(1010);
      h2->GetZaxis()->SetTitle("Rel. Mean(q) [%]");
      h2->GetZaxis()->CenterTitle();
      h2->Draw("colz");
      MakeDetectorPlot(ily, "p");
    }
    cOut->SaveAs(Form("%s.gif", cOut->GetName()));
  }
  for(Int_t ityp(0); ityp<(HasMCdata()?2:1); ityp++){
    if((arr = (TObjArray*)fProj->FindObject(ityp?"hCluster2MC":"hCluster2Track"))){
      for(Int_t iview(0); iview<nClViews; iview++){
        cOut = new TCanvas(Form("%s_%s%s_%d", GetName(), typName[ityp], vClName[iview], vClOpt[iview]), "Cluster Resolution", nx, ny);
        cOut->Divide(3,2, 1.e-5, 1.e-5);
        Int_t nplot(0);
        for(Int_t ily(0); ily<AliTRDgeometry::kNlayer; ily++){
          p=cOut->cd(ily+1);    p->SetRightMargin(0.1572581);p->SetTopMargin(0.08262712);
          if(!(h2 = (TH2*)arr->FindObject(Form("H%s%s%d_2D", typName[ityp], vClName[iview], UseLYselectTrklt()?fLYselect:ily)))) continue;
          nplot++;
          if(vClOpt[iview]==0) h2->Draw("colz");
          else if(vClOpt[iview]==1) h2e[ih2e++] = DrawSigma(h2, "#sigma(#Deltay) [#mum]", vClMin[iview], vClMax[iview], 1.e4);
          if(iview<5) MakeDetectorPlot(ily);
        }
        if(nplot==AliTRDgeometry::kNlayer) cOut->SaveAs(Form("%s.gif", cOut->GetName()));
        else delete cOut;
      }
    }
    // tracklet systematic
    if((arr = (TObjArray*)fProj->FindObject(ityp?"hTracklet2MC":"hTracklet2Track"))){
      for(Int_t iview(0); iview<nTrkltViews; iview++){
        cOut = new TCanvas(Form("%s_%s%s_%d", GetName(), typName[ityp], vTrkltName[iview], vTrkltOpt[iview]), "Tracklet Resolution", nx, ny);
        cOut->Divide(3,2, 1.e-5, 1.e-5);
        Int_t nplot(0);
        for(Int_t iplot(0); iplot<6; iplot++){
          p=cOut->cd(iplot+1); p->SetRightMargin(0.1572581); p->SetTopMargin(0.08262712);
          if(!(h2 = (TH2*)arr->FindObject(Form("H%s%s%d_2D", typName[ityp], vTrkltName[iview], iplot)))) continue;
          nplot++;
          if(vTrkltOpt[iview]==0) h2->Draw("colz");
          else if (vTrkltOpt[iview]==1) h2e[ih2e++] = DrawSigma(h2, "#sigma(#Deltay) [cm]", .15, .4);
          else if (vTrkltOpt[iview]==2) h2e[ih2e++] = DrawSigma(h2, "#sigma(occupancy) [%]", 10.5, 15.);
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
    for(Int_t iv0(0); iv0<2; iv0++){
      const Char_t *prefix = ityp?"MC":(iv0?"V0":"");
      if(ityp && !iv0) arr = (TObjArray*)fProj->FindObject("hTRDin2MC");
      else if(!ityp && !iv0) arr = (TObjArray*)fProj->FindObject("hTracklet2TRDin");
      else if(!ityp && iv0) arr = (TObjArray*)fProj->FindObject("hTracklet2TRDinV0");
      else continue;
      if(arr){
        for(Int_t iview(0); iview<nTrkInViews; iview++){
          cOut = new TCanvas(Form("%s_%s%s_%d", GetName(), prefix, vTrkInName[iview][0], vTrkInOpt[iview]), "Track IN Resolution", nx, ny);
          cOut->Divide(3,2, 1.e-5, 1.e-5);
          Int_t nplot(0);
          for(Int_t iplot(0); iplot<6; iplot++){
            p=cOut->cd(iplot+1); p->SetRightMargin(0.1572581); p->SetTopMargin(0.08262712);
            if(!(h2 = (TH2*)arr->FindObject(Form("H%s%s_2D", prefix, vTrkInName[iview][iplot])))){
              AliInfo(Form("Missing H%s%s_2D", prefix, vTrkInName[iview][iplot]));
              continue;
            }
            nplot++;
            if(vTrkInOpt[iview]==0) h2->Draw("colz");
            else h2e[ih2e++] = DrawSigma(h2, ttt[iplot], min[iplot], max[iplot]);
            MakeDetectorPlot(0);
          }
          if(nplot==6) cOut->SaveAs(Form("%s.gif", cOut->GetName()));
          else delete cOut;
        }
        // species
        const Float_t zmin[] = {1., 61., 15.},
                      zmax[] = {10., 79., 33.};
        cOut = new TCanvas(Form("%s_%sTrkInSpc", GetName(), prefix), "Track IN PID", Int_t(1.5*ny), Int_t(1.5*ny));
        cOut->Divide(3,3, 1.e-5, 1.e-5);
        Int_t nplot(0); const Char_t *chName[] = {"p", "n", ""};
        for(Int_t ich(0), ipad(1); ich<3; ich++){
          TH2 *h2s(NULL);
          if(!(h2s = (TH2*)arr->FindObject(Form("H%sTrkInY%sEn", prefix, chName[ich])))) {
            AliInfo(Form("Missing H%sTrkIn%sEn", prefix, chName[ich]));
            continue;
          }
          Int_t irebin(0), dxBin(1), dyBin(1);
          const Int_t nrebin(5); Int_t rebinX[nrebin] = {1, 2, 1, 2, 1}, rebinY[nrebin] = {2, 1, 2, 1, 2};
          if(h2s->GetNbinsY()==180){ // old binning
            rebinX[0] = 1; rebinY[0] = 2;
            rebinX[1] = 2; rebinY[1] = 1;
            rebinX[2] = 2; rebinY[2] = 1;
            rebinX[3] = 1; rebinY[3] = 5;
            rebinX[4] = 1; rebinY[4] = 1; // dummy
          }
          while(irebin<nrebin && (AliTRDrecoTask::GetMeanStat(h2s, .5, 1)<200)){
            h2s->Rebin2D(rebinX[irebin], rebinY[irebin]);
            dxBin*=rebinX[irebin];dyBin*=rebinY[irebin];irebin++;
          }
          AliDebug(2, Form("Rebin level[%d] @ chg[%d]. Binning[%2dx%2d]", irebin, ich, h2s->GetNbinsX(), h2s->GetNbinsY()));
          for(Int_t ispec(0); ispec<kNspc; ispec++){
            p=cOut->cd(ipad++); p->SetRightMargin(0.1572581); p->SetTopMargin(0.08262712);
            if(!(h2 = (TH2*)arr->FindObject(Form("H%sTrkInY%s%dEn", prefix, chName[ich], ispec)))) {
              AliInfo(Form("Missing H%sTrkIn%s%dEn", prefix, chName[ich], ispec));
              continue;
            }
            nplot++;
            h2->Rebin2D(dxBin,dyBin);
            h2->Divide(h2, h2s, 1.e2);
            h2->SetContour(9);
            h2->GetZaxis()->SetRangeUser(zmin[ispec], zmax[ispec]);
            h2->GetZaxis()->SetTitle("Rel. Abundancy [%]");h2->GetZaxis()->CenterTitle();
            TString tit(h2->GetTitle()); TObjArray *atit = tit.Tokenize("::");
            h2->SetTitle(Form("%s :: Relative Abundancy", ((*atit)[0])->GetName()));
            atit->Delete(); delete atit;
            h2->Draw("colz");
            MakeDetectorPlot(0);
          }
        }
        if(nplot==9) cOut->SaveAs(Form("%s.gif", cOut->GetName()));
        else delete cOut;
        // pt resolution plot
        cOut = new TCanvas(Form("%s_%sTrkInPt", GetName(), prefix), "TrackIn Pt", Int_t(1.5*ny), Int_t(1.5*ny));
        cOut->Divide(3,2, 1.e-5, 1.e-5);
        for(Int_t ich(0), ipad(1); ich<2; ich++){
          for(Int_t ispc(0); ispc<kNspc; ispc++){
            p=cOut->cd(ipad++); p->SetRightMargin(0.1572581);p->SetTopMargin(0.08262712);
            if(!(h2 = (TH2*)arr->FindObject(Form("H%sTrkInPt%s%d_2D", prefix, chName[ich], ispc)))) continue;
            h2->GetZaxis()->SetRangeUser(0.5, 1.2); h2->SetContour(7); h2->Draw("colz");
            MakeDetectorPlot(0);
          }
        }
        cOut->SaveAs(Form("%s.gif", cOut->GetName()));
        // MPV(Q) & <Q plots>
        const char *chQ[] = {"Q", "QS"};
        for(Int_t iq(0); iq<2; iq++){
          cOut = new TCanvas(Form("%s_%sTrkIn%s", GetName(), prefix, chQ[iq]), "Track IN PID", Int_t(1.5*ny), Int_t(1.5*ny));
          cOut->Divide(3,3, 1.e-5, 1.e-5);
          nplot=0;
          for(Int_t ich(0), ipad(1); ich<3; ich++){
            for(Int_t ispec(0); ispec<kNspc; ispec++){
              p=cOut->cd(ipad++); p->SetRightMargin(0.1572581); p->SetTopMargin(0.08262712);
              if(!(h2 = (TH2*)arr->FindObject(Form("H%sTrkIn%s%s%d_2D", prefix, chQ[iq], chName[ich], ispec)))) {
                AliInfo(Form("Missing H%sTrkIn%s%s%d_2D", prefix, chQ[iq], chName[ich], ispec));
                continue;
              }
              nplot++;
              h2->Draw("colz");
              MakeDetectorPlot(0);
            }
          }
          if(nplot==9) cOut->SaveAs(Form("%s.gif", cOut->GetName()));
          else delete cOut;
        }
      }
    }
  }
  // track MC systematic
  if((arr = (TObjArray*)fProj->FindObject("hTRD2MC"))) {
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

  // clean histos
  for(Int_t ih(ih2e); ih--;) delete h2e[ih];
}

//________________________________________________________
TH2* AliTRDresolution::DrawSigma(TH2 *h2, const Char_t *title, Float_t m, Float_t M, Float_t scale)
{
  // Draw error bars scaled with "scale" instead of content values
  //use range [m,M] if limits are specified

  if(!h2) return NULL;
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
  return h2e;
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
    AliInfo("Missing/Wrong data @ hDet2Cluster.");
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
  const Int_t nEtaPhi(5); Int_t rebinEtaPhiX[nEtaPhi] = {1, 2, 1, 2, 1}, rebinEtaPhiY[nEtaPhi] = {2, 1, 2, 1, 2};
  if(aa[1]->GetNbins()==180){ // old binning
    rebinEtaPhiX[0] = 1; rebinEtaPhiY[0] = 2;
    rebinEtaPhiX[1] = 2; rebinEtaPhiY[1] = 1;
    rebinEtaPhiX[2] = 2; rebinEtaPhiY[2] = 1;
    rebinEtaPhiX[3] = 1; rebinEtaPhiY[3] = 5;
    rebinEtaPhiX[4] = 1; rebinEtaPhiY[4] = 1; // dummy
  }
  //const Char_t *cenName[AliTRDeventInfo::kCentralityClasses] = {"0-10%", "10-20%", "20-50%", "50-80%", "80-100%"};
  const Char_t *cenName[AliTRDeventInfo::kCentralityClasses] = {"2800-inf", "2100-2799", "1400-2099", "700-1399", "0-699"};
  AliTRDrecoProjection hp[kDetNproj];  TObjArray php(kDetNproj);
  Int_t ih(0), isel(-1), np[nsel]={0};
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
    ((AliTRDrecoProjection*)php.At(isel))->Increment(coord, v);
    //Int_t ioff=isel;for(Int_t jh(0); jh<np[isel]; jh++) ((AliTRDrecoProjection*)php.At(ioff+jh))->Increment(coord, v);
  }
  TObjArray *arr(NULL);
  fProj->AddAt(arr = new TObjArray(kDetNproj), kDetector);
  arr->SetName("hDet2Cluster");  arr->SetOwner();

  TH2 *h2(NULL);  Int_t jh(0);
  for(; ih--; ){
    if(!hp[ih].H()) continue;
    if((h2 = hp[ih].Projection2D(kNstat, kNcontours, 0, kFALSE))) arr->AddAt(h2, jh++);
    if((h2 = (TH2*)gDirectory->Get(Form("%sEn", hp[ih].H()->GetName())))) arr->AddAt(h2, jh++);
  }
  AliTRDrecoProjection *pr0(NULL), *pr1(NULL);
  for(Int_t ily(0); ily<AliTRDgeometry::kNlayer; ily++){
    for(Int_t icen(0); icen<AliTRDeventInfo::kCentralityClasses; icen++){
      if((pr0 = (AliTRDrecoProjection*)php.FindObject(Form("HDet%d%d%d", ily, icen, 0)))){
        for(Int_t ipad(1); ipad<nPad; ipad++){
          if((pr1 = (AliTRDrecoProjection*)php.FindObject(Form("HDet%d%d%d", ily, icen, ipad)))){
            (*pr0)+=(*pr1);
          }
        }
        pr0->H()->SetNameTitle(Form("HDet%d%d", ily, icen), Form("Detectors :: Ly[%d] Cen[%s]", ily, cenName[icen]));
        if((h2 = pr0->Projection2D(kNstat, kNcontours, 0, kFALSE))) arr->AddAt(h2, jh++);
        if((h2 = (TH2*)gDirectory->Get(Form("%sEn", pr0->H()->GetName())))) arr->AddAt(h2, jh++);
        if(icen && (pr1 = (AliTRDrecoProjection*)php.FindObject(Form("HDet%d%d%d", ily, 0, 0)))) (*pr1)+=(*pr0);
      }
    }
    if((pr0 = (AliTRDrecoProjection*)php.FindObject(Form("HDet%d%d%d", ily, 0, 0)))){
      pr0->H()->SetNameTitle(Form("HDet%d", UseLYselectTrklt()?fLYselect:ily), Form("Detectors :: Ly[%d]", UseLYselectTrklt()?fLYselect:ily));
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
    return kTRUE;
  }
  Int_t ndim(H->GetNdimensions()); Bool_t debug(ndim>Int_t(kNdimCl));
  Int_t coord[10]; memset(coord, 0, sizeof(Int_t) * 10); Double_t v = 0.;
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
  const Int_t nEtaPhi(5); Int_t rebinEtaPhiX[nEtaPhi] = {1, 3, 1, 3, 1}, rebinEtaPhiY[nEtaPhi] = {2, 1, 2, 1, 2};
  if(aa[1]->GetNbins()==180){ // old binning
    rebinEtaPhiX[0] = 1; rebinEtaPhiY[0] = 2;
    rebinEtaPhiX[1] = 2; rebinEtaPhiY[1] = 1;
    rebinEtaPhiX[2] = 5; rebinEtaPhiY[2] = 1;
    rebinEtaPhiX[3] = 1; rebinEtaPhiY[3] = 5;
    rebinEtaPhiX[4] = 1; rebinEtaPhiY[4] = 1; // dummy
  }
  AliTRDrecoProjection hp[kClNproj];  TObjArray php(kClNproj);
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
  AliInfo(Form("Build %3d 3D %s projections.", ih, mc?"MC":""));

  AliTRDrecoProjection *pr0(NULL), *pr1(NULL);
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
      if(mc && ly<0) ly=0; // fix for bug in PlotMC
      isel = ly; ioff = isel;
    }
    if(ioff>=ih){
      AliError(Form("Wrong selection %d [%3d]", ioff, ih));
      return kFALSE;
    }
    if(!(pr0=(AliTRDrecoProjection*)php.At(ioff))) {
      AliError(Form("Missing projection %d", ioff));
      return kFALSE;
    }
    if(strcmp(pr0->H()->GetName(), Form("H%sClY%c%d%d%d", mc?"MC":"", chName[ch], ly, cen, npad))!=0){
      AliError(Form("Projection mismatch :: request[H%sClY%c%d%d%d] found[%s]", mc?"MC":"", chName[ch], ly, cen, npad, pr0->H()->GetName()));
      return kFALSE;
    }
    for(Int_t jh(0); jh<np[isel]; jh++) ((AliTRDrecoProjection*)php.At(ioff+jh))->Increment(coord, v);
  }
  if(HasDump3DFor(kCluster)){
    TDirectory *cwd = gDirectory;
    TFile::Open(Form("DumpRes_%s.root", H->GetName()), "RECREATE");
    for(Int_t ip(0); ip<php.GetEntriesFast(); ip++){
      if(!(pr0 = (AliTRDrecoProjection*)php.At(ip))) continue;
      if(!pr0->H()) continue;
      TH3 *h3=(TH3*)pr0->H()->Clone();
      h3->Write();
    }
    gFile->Close();
    cwd->cd();
  }

  TObjArray *arr(NULL);
  fProj->AddAt(arr = new TObjArray(kClNproj), cidx);
  arr->SetName(projName[Int_t(mc)]);  arr->SetOwner();

  TH2 *h2(NULL);  Int_t jh(0);
  for(; ih--; ){
    if(!hp[ih].H()) continue;
    if(strchr(hp[ih].H()->GetName(), 'Q')){
      if((h2 = hp[ih].Projection2D(kNstat, kNcontours, 0, kFALSE))) arr->AddAt(h2, jh++);
      if((h2 = (TH2*)gDirectory->Get(Form("%sEn", hp[ih].H()->GetName())))) arr->AddAt(h2, jh++);
    } else {
      if((h2 = hp[ih].Projection2D(kNstat, kNcontours, 1))) arr->AddAt(h2, jh++);
    }
  }
  Double_t m(0.), s(0.), se(0.), trend(0.);
  for(Int_t ily(0); ily<AliTRDgeometry::kNlayer; ily++){
    for(Int_t ich(0); ich<nCh; ich++){
      for(Int_t icen(0); icen<nCen; icen++){
        /*!dy*/
        if((pr0 = (AliTRDrecoProjection*)php.FindObject(Form("H%sClY%c%d%d%d", mc?"MC":"", chName[ich], ily, icen, 0)))){
          for(Int_t ipad(1); ipad<nNpad; ipad++){
            if(!(pr1 = (AliTRDrecoProjection*)php.FindObject(Form("H%sClY%c%d%d%d", mc?"MC":"", chName[ich], ily, icen, ipad)))) continue;
            (*pr0)+=(*pr1);
          }
          pr0->H()->SetNameTitle(Form("H%sClY%c%d%d", mc?"MC":"", chName[ich], ily, icen), Form("Clusters[%c] :: #Deltay Ly[%d] Cen[%s]", chSgn[ich], ily, cenName[icen]));
          if((h2 = pr0->Projection2D(kNstat, kNcontours, 1))) arr->AddAt(h2, jh++);
          if(icen && (pr1 = (AliTRDrecoProjection*)php.FindObject(Form("H%sClY%c%d%d%d", mc?"MC":"", chName[ich], ily, 0, 0)))) (*pr1)+=(*pr0);
        }
        /*!Q*/
        if((pr0 = (AliTRDrecoProjection*)php.FindObject(Form("H%sClQ%c%d%d%d", mc?"MC":"", chName[ich], ily, icen, 0)))){
          for(Int_t ipad(1); ipad<nNpad; ipad++){
            if(!(pr1 = (AliTRDrecoProjection*)php.FindObject(Form("H%sClQ%c%d%d%d", mc?"MC":"", chName[ich], ily, icen, ipad)))) continue;
            (*pr0)+=(*pr1);
          }
          pr0->H()->SetNameTitle(Form("H%sClQ%c%d%d", mc?"MC":"", chName[ich], ily, icen), Form("Clusters[%c] :: Q Ly[%d] Cen[%s]", chSgn[ich], ily, cenName[icen]));
          if((h2 = pr0->Projection2D(kNstat, kNcontours, 2, kFALSE))) arr->AddAt(h2, jh++);
          if((h2 = (TH2*)gDirectory->Get(Form("%sEn", pr0->H()->GetName())))) arr->AddAt(h2, jh++);
          pr0->H()->SetName(Form("H%sClQS%c%d%d", mc?"MC":"", chName[ich], ily, icen));
          if((h2 = pr0->Projection2D(kNstat, kNcontours, 0))) arr->AddAt(h2, jh++);
          if(icen && (pr1 = (AliTRDrecoProjection*)php.FindObject(Form("H%sClQ%c%d%d%d", mc?"MC":"", chName[ich], ily, 0, 0)))) (*pr1)+=(*pr0);
        }
        /*!YXTC*/
        if((pr0 = (AliTRDrecoProjection*)php.FindObject(Form("H%sClYXTC%c%d%d%d", mc?"MC":"", chName[ich], ily, icen, 0)))){
          for(Int_t ipad(1); ipad<nNpad; ipad++){
            if(!(pr1 = (AliTRDrecoProjection*)php.FindObject(Form("H%sClYXTC%c%d%d%d", mc?"MC":"", chName[ich], ily, icen, ipad)))) continue;
            (*pr0)+=(*pr1);
          }
          pr0->H()->SetNameTitle(Form("H%sClYXTC%c%d%d", mc?"MC":"", chName[ich], ily, icen), Form("Clusters[%c] :: #Deltay(x,TC) Ly[%d] Cen[%s]", chSgn[ich], ily, cenName[icen]));
          if((h2 = pr0->Projection2D(kNstat, kNcontours, 1))) arr->AddAt(h2, jh++);
          if(icen && (pr1 = (AliTRDrecoProjection*)php.FindObject(Form("H%sClYXTC%c%d%d%d", mc?"MC":"", chName[ich], ily, 0, 0)))) (*pr1)+=(*pr0);
        }
        /*!YXPh*/
        if((pr0 = (AliTRDrecoProjection*)php.FindObject(Form("H%sClYXPh%c%d%d%d", mc?"MC":"", chName[ich], ily, icen, 0)))){
          for(Int_t ipad(1); ipad<nNpad; ipad++){
            if(!(pr1 = (AliTRDrecoProjection*)php.FindObject(Form("H%sClYXPh%c%d%d%d", mc?"MC":"", chName[ich], ily, icen, ipad)))) continue;
            (*pr0)+=(*pr1);
          }
          pr0->H()->SetNameTitle(Form("H%sClYXPh%c%d%d", mc?"MC":"", chName[ich], ily, icen), Form("Clusters[%c] :: #Deltay(x,#Phi) Ly[%d] Cen[%s]", chSgn[ich], ily, cenName[icen]));
          if((h2 = pr0->Projection2D(kNstat, kNcontours, 1))) arr->AddAt(h2, jh++);
          if(icen && (pr1 = (AliTRDrecoProjection*)php.FindObject(Form("H%sClYXPh%c%d%d%d", mc?"MC":"", chName[ich], ily, 0, 0)))) (*pr1)+=(*pr0);
        }
      } // end centrality integration
      /*!dy*/
      if((pr0 = (AliTRDrecoProjection*)php.FindObject(Form("H%sClY%c%d%d%d", mc?"MC":"", chName[ich], ily, 0, 0)))){
        pr0->H()->SetNameTitle(Form("H%sClY%c%d", mc?"MC":"", chName[ich], UseLYselectTrklt()?fLYselect:ily),
                              Form("Clusters[%c]:: #Deltay Ly[%d]", chSgn[ich], UseLYselectTrklt()?fLYselect:ily));
        if((h2 = pr0->Projection2D(kNstat, kNcontours, 1))) arr->AddAt(h2, jh++);
        if(ich && (pr1 = (AliTRDrecoProjection*)php.FindObject(Form("H%sClY%c%d%d%d", mc?"MC":"", chName[0], ily, 0, 0)))) (*pr1)+=(*pr0);
      }
      /*!Q*/
      if((pr0 = (AliTRDrecoProjection*)php.FindObject(Form("H%sClQ%c%d%d%d", mc?"MC":"", chName[ich], ily, 0, 0)))){
        pr0->H()->SetNameTitle(Form("H%sClQ%c%d", mc?"MC":"", chName[ich], UseLYselectTrklt()?fLYselect:ily),
                              Form("Clusters[%c]:: Q Ly[%d]", chSgn[ich], UseLYselectTrklt()?fLYselect:ily));
        if((h2 = pr0->Projection2D(kNstat, kNcontours, 2))) arr->AddAt(h2, jh++);
        pr0->H()->SetName(Form("H%sClQS%c%d", mc?"MC":"", chName[ich], UseLYselectTrklt()?fLYselect:ily));
        if((h2 = pr0->Projection2D(kNstat, kNcontours, 0))) arr->AddAt(h2, jh++);
        if(ich && (pr1 = (AliTRDrecoProjection*)php.FindObject(Form("H%sClQ%c%d%d%d", mc?"MC":"", chName[0], ily, 0, 0)))) (*pr1)+=(*pr0);
      }
      /*!YXTC*/
      if((pr0 = (AliTRDrecoProjection*)php.FindObject(Form("H%sClYXTC%c%d%d%d", mc?"MC":"", chName[ich], ily, 0, 0)))){
        pr0->H()->SetNameTitle(Form("H%sClYXTC%c%d", mc?"MC":"", chName[ich], UseLYselectTrklt()?fLYselect:ily),
                              Form("Clusters[%c]:: #Deltay(x,TC) Ly[%d]", chSgn[ich], UseLYselectTrklt()?fLYselect:ily));
        if((h2 = pr0->Projection2D(kNstat, kNcontours, 1))) arr->AddAt(h2, jh++);
        if(ich && (pr1 = (AliTRDrecoProjection*)php.FindObject(Form("H%sClYXTC%c%d%d%d", mc?"MC":"", chName[0], ily, 0, 0)))) (*pr1)+=(*pr0);
      }
      /*!YXPh*/
      if((pr0 = (AliTRDrecoProjection*)php.FindObject(Form("H%sClYXPh%c%d%d%d", mc?"MC":"", chName[ich], ily, 0, 0)))){
        pr0->H()->SetNameTitle(Form("H%sClYXPh%c%d", mc?"MC":"", chName[ich], UseLYselectTrklt()?fLYselect:ily),
                              Form("Clusters[%c]:: #Deltay(x,#Phi) Ly[%d]", chSgn[ich], UseLYselectTrklt()?fLYselect:ily));
        if((h2 = pr0->Projection2D(kNstat, kNcontours, 1))) arr->AddAt(h2, jh++);
        if(ich && (pr1 = (AliTRDrecoProjection*)php.FindObject(Form("H%sClYXPh%c%d%d%d", mc?"MC":"", chName[0], ily, 0, 0)))) (*pr1)+=(*pr0);
      }
    } // end charge integration
    /*!dy*/
    if((pr0 = (AliTRDrecoProjection*)php.FindObject(Form("H%sClY%c%d%d%d", mc?"MC":"", chName[0], ily, 0, 0)))){
      pr0->H()->SetNameTitle(Form("H%sClY%d", mc?"MC":"", UseLYselectTrklt()?fLYselect:ily), Form("Clusters :: #Deltay Ly[%d]", UseLYselectTrklt()?fLYselect:ily));
      if((h2 = pr0->Projection2D(kNstat, kNcontours, 1))) arr->AddAt(h2, jh++);
      if((trend=pr0->GetTrendValue(1,&m,&s,&se))>-100.) PutTrendValue(Form("%sClS%d", mc?"MC":"", UseLYselectTrklt()?fLYselect:ily), s, se);
    }
    /*!YXPh*/
    if((pr0 = (AliTRDrecoProjection*)php.FindObject(Form("H%sClYXPh%c%d%d%d", mc?"MC":"", chName[0], ily, 0, 0)))){
      pr0->H()->SetNameTitle(Form("H%sClYXPh%d", mc?"MC":"", UseLYselectTrklt()?fLYselect:ily), Form("Clusters :: #Deltay Ly[%d]", UseLYselectTrklt()?fLYselect:ily));
      if((h2 = pr0->Projection2D(kNstat, kNcontours, 1))) arr->AddAt(h2, jh++);
    }

  }
  AliInfo(Form("Done %3d 2D %s projections.", jh, mc?"MC":""));
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
  if(!fProj){
    AliError("Missing results container.");
    return kFALSE;
  }
  if(!fContainer){
    AliError("Missing data container.");
    return kFALSE;
  }
  const Char_t *projName[] = {"hTracklet2Track", "hTracklet2MC"};
  THnSparse *H(NULL);
  if(!(H = (THnSparse*)fContainer->FindObject(projName[Int_t(mc)]))){
    AliError(Form("Missing/Wrong data @ %s.", projName[Int_t(mc)]));
    return kTRUE;
  }
  const Int_t mdim(kNdim+8);
  Int_t ndim(H->GetNdimensions()); //Bool_t debug(ndim>Int_t(kNdimTrklt));
  Int_t coord[mdim]; memset(coord, 0, sizeof(Int_t) * mdim); Double_t v = 0.;
  TAxis *aa[mdim], *as(NULL), *ap(NULL), *ac(NULL); memset(aa, 0, sizeof(TAxis*) * mdim);
  for(Int_t id(0); id<ndim; id++) aa[id] = H->GetAxis(id);
  if(ndim > Int_t(kSpeciesChgRC)) as = H->GetAxis(kSpeciesChgRC); // init species/charge selection
  if(ndim > Int_t(kPt))           ap = H->GetAxis(kPt);           // init pt selection
//  if(ndim > Int_t(kNdim))         ac = H->GetAxis(kNdim);         // init centrality selection
  // calculate size depending on debug level
  const Int_t nCen(ac?Int_t(AliTRDeventInfo::kCentralityClasses):1);
  const Int_t nPt(ap?Int_t(fNpt+2):1);
  const Int_t nSpc(as?Int_t(kNspc):1);
  const Int_t nCh(as?Int_t(kNcharge):1);
  const Int_t nLy(UseLYselectTrklt()?1:AliTRDgeometry::kNlayer);

  // build list of projections
  const Int_t nsel(AliTRDeventInfo::kCentralityClasses*AliTRDgeometry::kNlayer*(fgNPt+2)*(AliPID::kSPECIES*kNcharge + 1));
  // define rebinning strategy
  const Int_t nEtaPhi(5); Int_t rebinEtaPhiX[nEtaPhi] = {1, 3, 1, 3, 1}, rebinEtaPhiY[nEtaPhi] = {2, 1, 2, 1, 2};
  if(aa[1]->GetNbins()==180){ // old binning
    rebinEtaPhiX[0] = 1; rebinEtaPhiY[0] = 2;
    rebinEtaPhiX[1] = 2; rebinEtaPhiY[1] = 1;
    rebinEtaPhiX[2] = 5; rebinEtaPhiY[2] = 1;
    rebinEtaPhiX[3] = 1; rebinEtaPhiY[3] = 5;
    rebinEtaPhiX[4] = 1; rebinEtaPhiY[4] = 1; // dummy
  }
  AliTRDrecoProjection hp[kTrkltNproj]; TObjArray php(kTrkltNproj);
  Int_t ih(0), isel(-1), np[nsel]; memset(np, 0, nsel*sizeof(Int_t));
  const Char_t chName[kNcharge] = {'n', 'p'};const Char_t chSgn[kNcharge] = {'-', '+'};
  const Char_t ptShortName[5] = {'L', 'l', 'i', 'h', 'H'};
  Char_t ptName[fgNPt+2] = {0};
  Char_t *ptCut[fgNPt+2] = {NULL};

//  const Char_t *cenName[AliTRDeventInfo::kCentralityClasses] = {"0-10%", "10-20%", "20-50%", "50-80%", "80-100%"};
  const Char_t *cenName[AliTRDeventInfo::kCentralityClasses] = {"2800-inf", "2100-2799", "1400-2099", "700-1399", "0-699"};
  for(Int_t icen(0); icen<nCen; icen++){
    for(Int_t ily(0); ily<nLy; ily++){
      for(Int_t ipt(0); ipt<nPt; ipt++){
        if(!ptName[ipt]){
          ptName[ipt]= nPt>5?(64+ipt):ptShortName[ipt];
          ptCut[ipt] = StrDup(Form("#it{%4.2f<=p_{t}^{%s}[GeV/c]<%4.2f}",ipt?fgPt[ipt-1]:0., mc?"MC":"", ipt>fgNPt?99.99:fgPt[ipt]));
        }
        for(Int_t isp(0); isp<nSpc; isp++){
          for(Int_t ich(0); ich<nCh; ich++){
            isel++; // new selection
            AliDebug(3, Form("Building sel[%3d|%4d]", isel, ih));
            hp[ih].Build(Form("H%sTrkltY%c%c%d%d%d", mc?"MC":"", chName[ich], ptName[ipt], isp, UseLYselectTrklt()?fLYselect:ily, icen),
                        Form("Tracklets[%s%c]:: #Deltay{%s} Ly[%d] Cen[%s]", AliPID::ParticleLatexName(isp), chSgn[ich], ptCut[ipt], UseLYselectTrklt()?fLYselect:ily, cenName[icen]),
                        kEta, kPhi, kYrez, aa);
            //hp[ih].SetShowRange(-0.1,0.1);
            hp[ih].SetRebinStrategy(nEtaPhi, rebinEtaPhiX, rebinEtaPhiY);
            php.AddLast(&hp[ih++]); np[isel]++;
            hp[ih].Build(Form("H%sTrkltPh%c%c%d%d%d", mc?"MC":"", chName[ich], ptName[ipt], isp, UseLYselectTrklt()?fLYselect:ily, icen),
                        Form("Tracklets[%s%c]:: #Delta#phi{%s} Ly[%d] Cen[%s]", AliPID::ParticleLatexName(isp), chSgn[ich], ptCut[ipt], UseLYselectTrklt()?fLYselect:ily, cenName[icen]),
                        kEta, kPhi, kPrez, aa);
            //hp[ih].SetShowRange(-0.5,0.5);
            hp[ih].SetRebinStrategy(nEtaPhi, rebinEtaPhiX, rebinEtaPhiY);
            php.AddLast(&hp[ih++]); np[isel]++;
            hp[ih].Build(Form("H%sTrkltQ%c%c%d%d%d", mc?"MC":"", chName[ich], ptName[ipt], isp, UseLYselectTrklt()?fLYselect:ily, icen),
                        Form("Tracklets[%s%c]:: dQdl{%s} Ly[%d] Cen[%s]", AliPID::ParticleLatexName(isp), chSgn[ich], ptCut[ipt], UseLYselectTrklt()?fLYselect:ily, cenName[icen]),
                        kEta, kPhi, kZrez, aa);
            hp[ih].SetShowRange(1.,2.3);
            hp[ih].SetRebinStrategy(nEtaPhi, rebinEtaPhiX, rebinEtaPhiY);
            php.AddLast(&hp[ih++]); np[isel]++;
//             hp[ih].Build(Form("H%sTrkltTB%c%c%d%d%d", mc?"MC":"", chName[ich], ptName[ipt], isp, ily, icen),
//                         Form("Tracklets[%s%c]:: OccupancyTB{%s} Ly[%d] Cen[%s]", AliPID::ParticleLatexName(isp), chSgn[ich], ptCut[ipt], ily, cenName[icen]),
//                         kEta, kPhi, kNdim+1, aa);
//             hp[ih].SetShowRange(30., 70.);
//             hp[ih].SetRebinStrategy(nEtaPhi, rebinEtaPhiX, rebinEtaPhiY);
//             php.AddLast(&hp[ih++]); np[isel]++;
          }
        }
        if(ndim==kNdimTrklt) continue;

        isel++; // new selection
        AliDebug(3, Form("Building selRC[%3d|%4d]", isel, ih));
        hp[ih].Build(Form("H%sTrkltRCZ%c%d%d", mc?"MC":"", ptName[ipt], UseLYselectTrklt()?fLYselect:ily, icen),
                    Form("Tracklets[RC]:: #Deltaz{%s} Ly[%d] Cen[%s]", ptCut[ipt], UseLYselectTrklt()?fLYselect:ily, cenName[icen]),
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
  AliInfo(Form("Build %3d 3D %s projections%s", ih, mc?"MC":"", UseLYselectTrklt()?Form(" for Ly[%d].", fLYselect):"."));

  AliTRDrecoProjection *pr0(NULL), *pr1(NULL);
  Int_t ly(0), ch(0), sp(2), rcBin(as?as->FindBin(0.):-1), pt(0), cen(0), ioff(0), jspc(nSpc*nCh+1), kspc(nSpc*nCh*3/*4*/+1);
  for (Long64_t ib(0); ib < H->GetNbins(); ib++) {
    v = H->GetBinContent(ib, coord);
    if(v<1.) continue;
    ly = coord[kBC]-1; // layer selection
    if(UseLYselectTrklt()&& (ly!=fLYselect)) continue;
    // charge selection
    ch = 0; sp=0;// [e-] track [dafault]
    if(rcBin>0){ // debug mode in which species/charge are also saved
      sp = Int_t(TMath::Abs(as->GetBinCenter(coord[kSpeciesChgRC])))-1;
      if(coord[kSpeciesChgRC] > rcBin) ch = 1;  // [+] track
      else if(coord[kSpeciesChgRC] == rcBin) ch = 2;  // [RC] track
    }
    // pt selection
    pt = 0; // low pt
    if(ap) pt = coord[kPt];//TMath::Min(coord[kPt], Int_t(kNpt)+1);
    // centrality selection
    cen = 0; // default
    if(ac) cen = coord[kNdim]-1;
    // global selection
    if(ndim==kNdimTrklt){
      if(mc && ly<0) ly=0; // fix for bug in PlotMC
      ioff = ly*3/*4*/;
      isel = ly;
    } else {
      isel = cen*nLy*nPt*jspc+(!UseLYselectTrklt())*ly*nPt*jspc+pt*jspc; isel+=sp<0?(nSpc*nCh):(sp*nCh+ch);
      ioff = cen*nLy*nPt*kspc+(!UseLYselectTrklt())*ly*nPt*kspc+pt*kspc; ioff+=sp<0?((nSpc*nCh)*3/*4*/):(sp*nCh*3+ch*3/*4*/);
    }
    if(ioff>=ih){
      AliError(Form("Wrong selection %d [%3d]", ioff, ih));
      return kFALSE;
    }
    if(!(pr0=(AliTRDrecoProjection*)php.At(ioff))) {
      AliError(Form("Missing projection %d", ioff));
      return kFALSE;
    }
    if(sp>=0){
      if(strcmp(pr0->H()->GetName(), Form("H%sTrkltY%c%c%d%d%d", mc?"MC":"", chName[ch], ptName[pt], sp, UseLYselectTrklt()?fLYselect:ly, cen))!=0){
        AliError(Form("Projection mismatch :: request[H%sTrkltY%c%c%d%d%d] found[%s]", mc?"MC":"", chName[ch], ptName[pt], sp, UseLYselectTrklt()?fLYselect:ly, cen, pr0->H()->GetName()));
        return kFALSE;
      }
    } else {
      if(strcmp(pr0->H()->GetName(), Form("H%sTrkltRCZ%c%d%d", mc?"MC":"", ptName[pt], UseLYselectTrklt()?fLYselect:ly, cen))!=0){
        AliError(Form("Projection mismatch :: request[H%sTrkltRCZ%c%d%d] found[%s]", mc?"MC":"", ptName[pt], UseLYselectTrklt()?fLYselect:ly, cen, pr0->H()->GetName()));
        return kFALSE;
      }
    }
    for(Int_t jh(0); jh<np[isel]; jh++) ((AliTRDrecoProjection*)php.At(ioff+jh))->Increment(coord, v);
  }
  if(HasDump3DFor(kTracklet)){
    TDirectory *cwd = gDirectory;
    TFile::Open(Form("DumpRes_%s.root", H->GetName()), "RECREATE");
    for(Int_t ip(0); ip<php.GetEntriesFast(); ip++){
      if(!(pr0 = (AliTRDrecoProjection*)php.At(ip))) continue;
      if(!pr0->H()) continue;
      TH3 *h3=(TH3*)pr0->H()->Clone();
      h3->Write();
    }
    gFile->Close();
    cwd->cd();
  }

  TH2 *h2(NULL); Int_t jh(0);
  TObjArray *arr(NULL);
  if(!(arr = (TObjArray*)fProj->At(cidx))){
    fProj->AddAt(arr = new TObjArray(kTrkltNproj), cidx);
    arr->SetName(projName[Int_t(mc)]); arr->SetOwner();
  } else jh = arr->GetEntriesFast();

  for(; ih--; ){
    if(!hp[ih].H()) continue;
    Int_t mid(0), nstat(kNstat);
    if(strchr(hp[ih].H()->GetName(), 'Q')){ mid=2; nstat=kNstatQ;}
    if(!(h2 = hp[ih].Projection2D(nstat, kNcontours, mid))) continue;
    arr->AddAt(h2, jh++);
  }
  // build combined performance plots
  Double_t m(0.), s(0.), se(0.), trend(0.);
  for(Int_t ily(0); ily<nLy; ily++){
    for(Int_t ich(0); ich<nCh; ich++){
      for(Int_t ipt(0); ipt<nPt; ipt++){
        for(Int_t icen(0); icen<nCen; icen++){
          /*!dy*/
          if((pr0 = (AliTRDrecoProjection*)php.FindObject(Form("H%sTrkltY%c%c%d%d%d", mc?"MC":"", chName[ich], ptName[ipt], 0, UseLYselectTrklt()?fLYselect:ily, icen)))){
            for(Int_t isp(1); isp<nSpc; isp++){
              if(!(pr1 = (AliTRDrecoProjection*)php.FindObject(Form("H%sTrkltY%c%c%d%d%d", mc?"MC":"", chName[ich], ptName[ipt], isp, UseLYselectTrklt()?fLYselect:ily, icen)))) continue;
              (*pr0)+=(*pr1);
            }
            pr0->H()->SetNameTitle(Form("H%sTrkltY%c%c%d%d", mc?"MC":"", chName[ich], ptName[ipt], UseLYselectTrklt()?fLYselect:ily, icen),
                                      Form("Tracklets[%c]:: #Deltay{%s} Ly[%d] Cen[%s]", chSgn[ich], ptCut[ipt], UseLYselectTrklt()?fLYselect:ily, cenName[icen]));
            if((h2 = pr0->Projection2D(kNstat, kNcontours, 1))) arr->AddAt(h2, jh++);
            if(icen && (pr1 = (AliTRDrecoProjection*)php.FindObject(Form("H%sTrkltY%c%c%d%d%d", mc?"MC":"", chName[ich], ptName[ipt], 0, UseLYselectTrklt()?fLYselect:ily, 0)))) (*pr1)+=(*pr0);
          }
          /*!dphi*/
          if((pr0 = (AliTRDrecoProjection*)php.FindObject(Form("H%sTrkltPh%c%c%d%d%d", mc?"MC":"", chName[ich], ptName[ipt], 0, UseLYselectTrklt()?fLYselect:ily, icen)))){
            for(Int_t isp(1); isp<nSpc; isp++){
              if(!(pr1 = (AliTRDrecoProjection*)php.FindObject(Form("H%sTrkltPh%c%c%d%d%d", mc?"MC":"", chName[ich], ptName[ipt], isp, UseLYselectTrklt()?fLYselect:ily, icen)))) continue;
              (*pr0)+=(*pr1);
            }
            pr0->H()->SetNameTitle(Form("H%sTrkltPh%c%c%d%d", mc?"MC":"", chName[ich], ptName[ipt], UseLYselectTrklt()?fLYselect:ily, icen),
                                      Form("Tracklets[%c]:: #Delta#phi{%s} Ly[%d] Cen[%s]", chSgn[ich], ptCut[ipt], UseLYselectTrklt()?fLYselect:ily, cenName[icen]));
            if((h2 = pr0->Projection2D(kNstat, kNcontours, 1))) arr->AddAt(h2, jh++);
            if(icen && (pr1 = (AliTRDrecoProjection*)php.FindObject(Form("H%sTrkltPh%c%c%d%d%d", mc?"MC":"", chName[ich], ptName[ipt], 0, UseLYselectTrklt()?fLYselect:ily, 0)))) (*pr1)+=(*pr0);
          }
          /*!dQ/dl*/
          if((pr0 = (AliTRDrecoProjection*)php.FindObject(Form("H%sTrkltQ%c%c%d%d%d", mc?"MC":"", chName[ich], ptName[ipt], 0, UseLYselectTrklt()?fLYselect:ily, icen)))){
            for(Int_t isp(1); isp<nSpc; isp++){
              if(!(pr1 = (AliTRDrecoProjection*)php.FindObject(Form("H%sTrkltQ%c%c%d%d%d", mc?"MC":"", chName[ich], ptName[ipt], isp, UseLYselectTrklt()?fLYselect:ily, icen)))) continue;
              (*pr0)+=(*pr1);
            }
            pr0->H()->SetNameTitle(Form("H%sTrkltQ%c%c%d%d", mc?"MC":"", chName[ich], ptName[ipt], UseLYselectTrklt()?fLYselect:ily, icen),
                                      Form("Tracklets[%c]:: dQdl{%s} Ly[%d] Cen[%s]", chSgn[ich], ptCut[ipt], UseLYselectTrklt()?fLYselect:ily, cenName[icen]));
            if((h2 = pr0->Projection2D(kNstatQ, kNcontours, 2))) arr->AddAt(h2, jh++);
            pr0->H()->SetNameTitle(Form("H%sTrkltQS%c%c%d%d", mc?"MC":"", chName[ich], ptName[ipt], UseLYselectTrklt()?fLYselect:ily, icen),
                                      Form("Tracklets[%c]:: dQdl{%s} Ly[%d] Cen[%s]", chSgn[ich], ptCut[ipt], UseLYselectTrklt()?fLYselect:ily, cenName[icen]));
            pr0->SetShowRange(2.4, 5.1);
            if((h2 = pr0->Projection2D(kNstat, kNcontours, 0))) arr->AddAt(h2, jh++);
            if(icen && (pr1 = (AliTRDrecoProjection*)php.FindObject(Form("H%sTrkltQ%c%c%d%d%d", mc?"MC":"", chName[ich], ptName[ipt], 0, UseLYselectTrklt()?fLYselect:ily, 0)))) (*pr1)+=(*pr0);
          }
          /*!TB occupancy*/
//           if((pr0 = (AliTRDrecoProjection*)php.FindObject(Form("H%sTrkltTB%c%c%d%d%d", mc?"MC":"", chName[ich], ptName[ipt], 0, UseLYselectTrklt()?fLYselect:ily, icen)))){
//             for(Int_t isp(1); isp<nSpc; isp++){
//               if(!(pr1 = (AliTRDrecoProjection*)php.FindObject(Form("H%sTrkltTB%c%c%d%d%d", mc?"MC":"", chName[ich], ptName[ipt], isp, UseLYselectTrklt()?fLYselect:ily, icen)))) continue;
//               (*pr0)+=(*pr1);
//             }
//             pr0->H()->SetNameTitle(Form("H%sTrkltTB%c%c%d%d", mc?"MC":"", chName[ich], ptName[ipt], UseLYselectTrklt()?fLYselect:ily, icen),
//                                       Form("Tracklets[%c]:: OccupancyTB{%s} Ly[%d] Cen[%s]", chSgn[ich], ptCut[ipt], UseLYselectTrklt()?fLYselect:ily, cenName[icen]));
//             if((h2 = pr0->Projection2D(kNstat, kNcontours))) arr->AddAt(h2, jh++);
//             if(icen && (pr1 = (AliTRDrecoProjection*)php.FindObject(Form("H%sTrkltTB%c%c%d%d%d", mc?"MC":"", chName[ich], ptName[ipt], 0, UseLYselectTrklt()?fLYselect:ily, 0)))) (*pr1)+=(*pr0);
//           }
        } // end centrality integration
        /*!dy*/
        if((pr0 = (AliTRDrecoProjection*)php.FindObject(Form("H%sTrkltY%c%c%d%d%d", mc?"MC":"", chName[ich], ptName[ipt], 0, UseLYselectTrklt()?fLYselect:ily, 0)))){
          pr0->H()->SetNameTitle(Form("H%sTrkltY%c%c%d", mc?"MC":"", chName[ich], ptName[ipt], UseLYselectTrklt()?fLYselect:ily),
                                Form("Tracklets[%c]:: #Deltay Ly[%d] %s", chSgn[ich], UseLYselectTrklt()?fLYselect:ily, ptCut[ipt]));
          if((h2 = pr0->Projection2D(kNstat, kNcontours, 1))) arr->AddAt(h2, jh++);
          if((trend=pr0->GetTrendValue(1,&m,&s,&se))>-100.){
            PutTrendValue(Form("%sTrkltY%c%c%d", mc?"MC":"", chName[ich], ptName[ipt], UseLYselectTrklt()?fLYselect:ily), trend, m);
            PutTrendValue(Form("%sTrkltYS%c%c%d", mc?"MC":"", chName[ich], ptName[ipt], UseLYselectTrklt()?fLYselect:ily), s, se);
          }
          if(ipt && (pr1 = (AliTRDrecoProjection*)php.FindObject(Form("H%sTrkltY%c%c%d%d%d", mc?"MC":"", chName[ich], ptName[0], 0, UseLYselectTrklt()?fLYselect:ily, 0)))) (*pr1)+=(*pr0);
        }
        /*!dphi*/
        if((pr0 = (AliTRDrecoProjection*)php.FindObject(Form("H%sTrkltPh%c%c%d%d%d", mc?"MC":"", chName[ich], ptName[ipt], 0, UseLYselectTrklt()?fLYselect:ily, 0)))){
          pr0->H()->SetNameTitle(Form("H%sTrkltPh%c%c%d", mc?"MC":"", chName[ich], ptName[ipt], UseLYselectTrklt()?fLYselect:ily),
                                Form("Tracklets[%c]:: #Delta#phi Ly[%d] %s", chSgn[ich], UseLYselectTrklt()?fLYselect:ily, ptCut[ipt]));
          if((h2 = pr0->Projection2D(kNstat, kNcontours, 1))) arr->AddAt(h2, jh++);
          if((trend=pr0->GetTrendValue(1,&m,&s,&se))>-100.){
            PutTrendValue(Form("%sTrkltPh%c%c%d", mc?"MC":"", chName[ich], ptName[ipt], UseLYselectTrklt()?fLYselect:ily), trend, m);
            PutTrendValue(Form("%sTrkltPhS%c%c%d", mc?"MC":"", chName[ich], ptName[ipt], UseLYselectTrklt()?fLYselect:ily), s, se);
          }
          if(ipt && (pr1 = (AliTRDrecoProjection*)php.FindObject(Form("H%sTrkltPh%c%c%d%d%d", mc?"MC":"", chName[ich], ptName[0], 0, UseLYselectTrklt()?fLYselect:ily, 0)))) (*pr1)+=(*pr0);
        }
        /*!dQ/dl*/
        if((pr0 = (AliTRDrecoProjection*)php.FindObject(Form("H%sTrkltQ%c%c%d%d%d", mc?"MC":"", chName[ich], ptName[ipt], 0, UseLYselectTrklt()?fLYselect:ily, 0)))){
          pr0->H()->SetNameTitle(Form("H%sTrkltQ%c%c%d", mc?"MC":"", chName[ich], ptName[ipt], UseLYselectTrklt()?fLYselect:ily),
                                Form("Tracklets[%c]:: dQdl Ly[%d] %s", chSgn[ich], UseLYselectTrklt()?fLYselect:ily, ptCut[ipt]));
          pr0->SetShowRange(1.,2.3);
          if((h2 = pr0->Projection2D(kNstatQ, kNcontours, 2))) arr->AddAt(h2, jh++);
          pr0->H()->SetNameTitle(Form("H%sTrkltQS%c%c%d", mc?"MC":"", chName[ich], ptName[ipt], UseLYselectTrklt()?fLYselect:ily),
                                Form("Tracklets[%c]:: dQdl Ly[%d] %s", chSgn[ich], UseLYselectTrklt()?fLYselect:ily, ptCut[ipt]));
          pr0->SetShowRange(2.4,5.1);
          if((h2 = pr0->Projection2D(kNstat, kNcontours, 0))) arr->AddAt(h2, jh++);
          if((trend=pr0->GetTrendValue(2,&m,&s,&se))>-100.){
            PutTrendValue(Form("%sTrkltQ%c%c%d", mc?"MC":"", chName[ich], ptName[ipt], UseLYselectTrklt()?fLYselect:ily), trend, m);
            PutTrendValue(Form("%sTrkltQS%c%c%d", mc?"MC":"", chName[ich], ptName[ipt], UseLYselectTrklt()?fLYselect:ily), s, se);
          }
          if(ipt && (pr1 = (AliTRDrecoProjection*)php.FindObject(Form("H%sTrkltQ%c%c%d%d%d", mc?"MC":"", chName[ich], ptName[0], 0, UseLYselectTrklt()?fLYselect:ily, 0)))) (*pr1)+=(*pr0);
        }
        /*!TB occupancy*/
//         if((pr0 = (AliTRDrecoProjection*)php.FindObject(Form("H%sTrkltTB%c%c%d%d%d", mc?"MC":"", chName[ich], ptName[ipt], 0, UseLYselectTrklt()?fLYselect:ily, 0)))){
//           pr0->H()->SetNameTitle(Form("H%sTrkltTB%c%c%d", mc?"MC":"", chName[ich], ptName[ipt], UseLYselectTrklt()?fLYselect:ily),
//                                 Form("Tracklets[%c]:: OccupancyTB Ly[%d] %s", chSgn[ich], UseLYselectTrklt()?fLYselect:ily, ptCut[ipt]));
//           if((h2 = pr0->Projection2D(kNstat, kNcontours))) arr->AddAt(h2, jh++);
//           if(ipt && (pr1 = (AliTRDrecoProjection*)php.FindObject(Form("H%sTrkltTB%c%c%d%d%d", mc?"MC":"", chName[ich], ptName[0], 0, UseLYselectTrklt()?fLYselect:ily, 0)))) (*pr1)+=(*pr0);
//         }
      } // end pt integration
      /*!dy*/
      if((pr0 = (AliTRDrecoProjection*)php.FindObject(Form("H%sTrkltY%c%c%d%d%d", mc?"MC":"", chName[ich], ptName[0], 0, UseLYselectTrklt()?fLYselect:ily, 0)))){
        pr0->H()->SetNameTitle(Form("H%sTrkltY%c%d", mc?"MC":"", chName[ich], UseLYselectTrklt()?fLYselect:ily), Form("Tracklets[%c] :: #Deltay Ly[%d]", chSgn[ich], UseLYselectTrklt()?fLYselect:ily));
        if((h2 = pr0->Projection2D(kNstat, kNcontours, 1))) arr->AddAt(h2, jh++);
        if((trend=pr0->GetTrendValue(1,&m,&s,&se))>-100.){
          PutTrendValue(Form("%sTrkltY%c%d", mc?"MC":"", chName[ich], UseLYselectTrklt()?fLYselect:ily), trend, m);
          PutTrendValue(Form("%sTrkltYS%c%d", mc?"MC":"", chName[ich], UseLYselectTrklt()?fLYselect:ily), s, se);
        }
        if(ich && (pr1 = (AliTRDrecoProjection*)php.FindObject(Form("H%sTrkltY%c%c%d%d%d", mc?"MC":"", chName[0], ptName[0], 0, UseLYselectTrklt()?fLYselect:ily, 0)))) (*pr1)+=(*pr0);
      }
      /*!dphi*/
      if((pr0 = (AliTRDrecoProjection*)php.FindObject(Form("H%sTrkltPh%c%c%d%d%d", mc?"MC":"", chName[ich], ptName[0], 0, UseLYselectTrklt()?fLYselect:ily, 0)))){
        pr0->H()->SetNameTitle(Form("H%sTrkltPh%c%d", mc?"MC":"", chName[ich], UseLYselectTrklt()?fLYselect:ily), Form("Tracklets[%c] :: #Delta#phi Ly[%d]", chSgn[ich], UseLYselectTrklt()?fLYselect:ily));
        pr0->SetShowRange(-.9,.9);
        if((h2 = pr0->Projection2D(kNstat, kNcontours, 1))) arr->AddAt(h2, jh++);
        if((trend=pr0->GetTrendValue(1,&m,&s,&se))>-100.){
          PutTrendValue(Form("%sTrkltPh%c%d", mc?"MC":"", chName[ich], UseLYselectTrklt()?fLYselect:ily), trend, m);
          PutTrendValue(Form("%sTrkltPhS%c%d", mc?"MC":"", chName[ich], UseLYselectTrklt()?fLYselect:ily), s, se);
        }
        if(ich && (pr1 = (AliTRDrecoProjection*)php.FindObject(Form("H%sTrkltPh%c%c%d%d%d", mc?"MC":"", chName[0], ptName[0], 0, UseLYselectTrklt()?fLYselect:ily, 0)))) (*pr1)+=(*pr0);
      }
      /*!dQ/dl*/
      if((pr0 = (AliTRDrecoProjection*)php.FindObject(Form("H%sTrkltQ%c%c%d%d%d", mc?"MC":"", chName[ich], ptName[0], 0, UseLYselectTrklt()?fLYselect:ily, 0)))){
        pr0->H()->SetNameTitle(Form("H%sTrkltQ%c%d", mc?"MC":"", chName[ich], UseLYselectTrklt()?fLYselect:ily), Form("Tracklets[%c] :: dQdl Ly[%d]", chSgn[ich], UseLYselectTrklt()?fLYselect:ily));
        pr0->SetShowRange(1.,2.3);
        if((h2 = pr0->Projection2D(kNstatQ, kNcontours, 2))) arr->AddAt(h2, jh++);
        pr0->H()->SetNameTitle(Form("H%sTrkltQS%c%d", mc?"MC":"", chName[ich], UseLYselectTrklt()?fLYselect:ily), Form("Tracklets[%c] :: dQdl Ly[%d]", chSgn[ich], UseLYselectTrklt()?fLYselect:ily));
        pr0->SetShowRange(2.4,5.1);
        if((h2 = pr0->Projection2D(kNstat, kNcontours, 0))) arr->AddAt(h2, jh++);
        if((trend=pr0->GetTrendValue(2,&m,&s,&se))>-100.){
          PutTrendValue(Form("%sTrkltQ%c%d", mc?"MC":"", chName[ich], UseLYselectTrklt()?fLYselect:ily), trend, m);
          PutTrendValue(Form("%sTrkltQS%c%d", mc?"MC":"", chName[ich], UseLYselectTrklt()?fLYselect:ily), s, se);
        }
        if(ich && (pr1 = (AliTRDrecoProjection*)php.FindObject(Form("H%sTrkltQ%c%c%d%d%d", mc?"MC":"", chName[0], ptName[0], 0, UseLYselectTrklt()?fLYselect:ily, 0)))) (*pr1)+=(*pr0);
      }
      /*!TB occupancy*/
//       if((pr0 = (AliTRDrecoProjection*)php.FindObject(Form("H%sTrkltTB%c%c%d%d%d", mc?"MC":"", chName[ich], ptName[0], 0, UseLYselectTrklt()?fLYselect:ily, 0)))){
//         pr0->H()->SetNameTitle(Form("H%sTrkltTB%c%d", mc?"MC":"", chName[ich], UseLYselectTrklt()?fLYselect:ily), Form("Tracklets[%c] :: OccupancyTB Ly[%d]", chSgn[ich], UseLYselectTrklt()?fLYselect:ily));
//         if((h2 = pr0->Projection2D(kNstat, kNcontours))) arr->AddAt(h2, jh++);
//         if(ich && (pr1 = (AliTRDrecoProjection*)php.FindObject(Form("H%sTrkltTB%c%c%d%d%d", mc?"MC":"", chName[0], ptName[0], 0, UseLYselectTrklt()?fLYselect:ily, 0)))) (*pr1)+=(*pr0);
//       }
    } // end charge integration
    /*!dy*/
    if((pr0 = (AliTRDrecoProjection*)php.FindObject(Form("H%sTrkltY%c%c%d%d%d", mc?"MC":"", chName[0], ptName[0], 0, UseLYselectTrklt()?fLYselect:ily, 0)))){
      pr0->H()->SetNameTitle(Form("H%sTrkltY%d", mc?"MC":"", UseLYselectTrklt()?fLYselect:ily), Form("Tracklets :: #Deltay Ly[%d]", UseLYselectTrklt()?fLYselect:ily));
      if((h2 = pr0->Projection2D(kNstat, kNcontours, 1))) arr->AddAt(h2, jh++);
      if((trend=pr0->GetTrendValue(1,&m,&s,&se))>-100.){
        PutTrendValue(Form("%sTrkltY%d", mc?"MC":"", UseLYselectTrklt()?fLYselect:ily), trend, m);
        PutTrendValue(Form("%sTrkltYS%d", mc?"MC":"", UseLYselectTrklt()?fLYselect:ily), s, se);
      }
    }
    /*!dphi*/
    if((pr0 = (AliTRDrecoProjection*)php.FindObject(Form("H%sTrkltPh%c%c%d%d%d", mc?"MC":"", chName[0], ptName[0], 0, UseLYselectTrklt()?fLYselect:ily, 0)))){
      pr0->H()->SetNameTitle(Form("H%sTrkltPh%d", mc?"MC":"", UseLYselectTrklt()?fLYselect:ily), Form("Tracklets :: #Delta#phi Ly[%d]", UseLYselectTrklt()?fLYselect:ily));
      pr0->SetShowRange(-.45,.45);
      if((h2 = pr0->Projection2D(kNstat, kNcontours, 1))) arr->AddAt(h2, jh++);
    }
    /*!dQdl*/
    if((pr0 = (AliTRDrecoProjection*)php.FindObject(Form("H%sTrkltQ%c%c%d%d%d", mc?"MC":"", chName[0], ptName[0], 0, UseLYselectTrklt()?fLYselect:ily, 0)))){
      pr0->H()->SetNameTitle(Form("H%sTrkltQ%d", mc?"MC":"", UseLYselectTrklt()?fLYselect:ily), Form("Tracklets :: dQdl Ly[%d]", UseLYselectTrklt()?fLYselect:ily));
      pr0->SetShowRange(1.,2.3);
      if((h2 = pr0->Projection2D(kNstat, kNcontours, 2))) arr->AddAt(h2, jh++);
      pr0->H()->SetName(Form("H%sTrkltQS%d", mc?"MC":"", UseLYselectTrklt()?fLYselect:ily));
      pr0->SetShowRange(2.4,5.1);
      if((h2 = pr0->Projection2D(kNstat, kNcontours, 0))) arr->AddAt(h2, jh++);
    }
    /*!TB occupancy*/
//     if((pr0 = (AliTRDrecoProjection*)php.FindObject(Form("H%sTrkltTB%c%c%d%d%d", mc?"MC":"", chName[0], ptName[0], 0, UseLYselectTrklt()?fLYselect:ily, 0)))){
//       pr0->H()->SetNameTitle(Form("H%sTrkltTB%d", mc?"MC":"", UseLYselectTrklt()?fLYselect:ily), Form("Tracklets :: OccupancyTB Ly[%d]", UseLYselectTrklt()?fLYselect:ily));
//       if((h2 = pr0->Projection2D(kNstat, kNcontours))) arr->AddAt(h2, jh++);
//     }

    /*! Row Cross processing*/
    for(Int_t icen(0); icen<nCen; icen++){
      /*!RC dz*/
      if((pr0 = (AliTRDrecoProjection*)php.FindObject(Form("H%sTrkltRCZ%c%d%d", mc?"MC":"", ptName[0], UseLYselectTrklt()?fLYselect:ily, icen)))){
        for(Int_t ipt(0); ipt<nPt; ipt++){
          if(!(pr1 = (AliTRDrecoProjection*)php.FindObject(Form("H%sTrkltRCZ%c%d%d", mc?"MC":"", ptName[ipt], UseLYselectTrklt()?fLYselect:ily, icen)))) continue;
          (*pr0)+=(*pr1);
        }
        pr0->H()->SetNameTitle(Form("H%sTrkltRCZ%d%d", mc?"MC":"", UseLYselectTrklt()?fLYselect:ily, icen), Form("Tracklets[RC]:: #Deltaz Ly[%d] Cen[%s]", UseLYselectTrklt()?fLYselect:ily, cenName[icen]));
        if((h2 = pr0->Projection2D(kNstat, kNcontours, 1))) arr->AddAt(h2, jh++);
        if(icen && (pr1 = (AliTRDrecoProjection*)php.FindObject(Form("H%sTrkltRCZ%c%d%d", mc?"MC":"", ptName[0], UseLYselectTrklt()?fLYselect:ily, 0)))) (*pr1)+=(*pr0);
      }
    } // end centrality integration for row cross
    /*!RC dz*/
    if((pr0 = (AliTRDrecoProjection*)php.FindObject(Form("H%sTrkltRCZ%c%d%d", mc?"MC":"", ptName[0], UseLYselectTrklt()?fLYselect:ily, 0)))){
      pr0->H()->SetNameTitle(Form("H%sTrkltRCZ%d", mc?"MC":"", UseLYselectTrklt()?fLYselect:ily), Form("Tracklets[RC] :: #Deltaz Ly[%d]", UseLYselectTrklt()?fLYselect:ily));
      if((h2 = pr0->Projection2D(kNstat, kNcontours, 1))) arr->AddAt(h2, jh++);
    }
  } // end layer loop
  AliInfo(Form("Done %3d 2D %s projections.", jh, mc?"MC":""));

// clean local memory allocation
  for(Int_t i(fgNPt+2);i--;) if(ptCut[i]) delete [] ptCut[i];
  return kTRUE;
}

//________________________________________________________
Bool_t AliTRDresolution::MakeProjectionTrackIn(Bool_t mc, Bool_t v0)
{
// Analyse track in
  const Int_t kNcontours(9);
  const Int_t kNstat(30);
  Int_t cidx=mc?kMCtrackIn:(v0?(kV0TrackIn-1):kTrackIn);
  if(fProj && fProj->At(cidx)){
    AliInfo(Form("Nothing to do for container %s.", ((TObjArray*)fProj->At(cidx))->GetName() ));
    return kTRUE;
  }
  if(!fContainer){
    AliError("Missing data container.");
    return kFALSE;
  }
  const Char_t *projName[] = {"hTracklet2TRDin", "hTRDin2MC", "hTracklet2TRDinV0"};
  THnSparse *H(NULL);
  if(!(H = (THnSparse*)fContainer->FindObject(projName[Int_t(mc)]))){
    AliError(Form("Missing/Wrong data @ %s.", projName[Int_t(mc)]));
    return kTRUE;
  }
  const Char_t *prefix = mc?"MC":(v0?"V0":"");

  const Int_t mdim(kNdim+3);
  Int_t coord[mdim]; memset(coord, 0, sizeof(Int_t) * mdim); Double_t v = 0.;
  Int_t ndim(H->GetNdimensions());
  TAxis *aa[mdim], *as(NULL), *ap(NULL), *apt(NULL), *ax(NULL), *abf(NULL); memset(aa, 0, sizeof(TAxis*) * (mdim));
  for(Int_t id(0); id<ndim; id++) aa[id] = H->GetAxis(id);
  if(ndim > Int_t(kSpeciesChgRC)) as = H->GetAxis(kSpeciesChgRC);
  if(ndim > Int_t(kPt))          apt = H->GetAxis(kPt);
  if(ndim > Int_t(kNdim))         ap = H->GetAxis(kNdim);
  if(ndim > Int_t(kNdim)+1)       ax = H->GetAxis(kNdim+1);
  if(ndim > Int_t(kNdim)+2)      abf = H->GetAxis(kNdim+2);
  //AliInfo(Form("Using : Species[%c] Pt[%c] BunchFill[%c]", as?'y':'n', ap?'y':'n', abf?'y':'n'));
  const Int_t nPt(apt?(apt->GetNbins()+2):1);

  // build list of projections
  const Int_t nsel((fgNPt+2)*(kNspc*kNcharge + 1) + kNspc*kNcharge);
  const Int_t nprj((fgNPt+2)*(kNspc*kNcharge + 1)*4 + kNspc*kNcharge);
  const Char_t chName[kNcharge] = {'n', 'p'};const Char_t chSgn[kNcharge] = {'-', '+'};
  const Char_t *spcName[2][kNspc] = {{"e", "#mu#pi", "Kp"},
                                     {"e", "#pi", "p"}};
  const Char_t ptShortName[5] = {'L', 'l', 'i', 'h', 'H'};
  Char_t ptName[fgNPt+2] = {0};
  Char_t *ptCut[fgNPt+2] = {NULL};
  Char_t  *pCut[fgNPt+2] = {NULL};
  // define rebinning strategy
  const Int_t nEtaPhi(5); Int_t rebinEtaPhiX[nEtaPhi] = {1, 3, 1, 3, 1}, rebinEtaPhiY[nEtaPhi] = {2, 1, 2, 1, 2};
  if(aa[1]->GetNbins()==180){ // old binning
    rebinEtaPhiX[0] = 1; rebinEtaPhiY[0] = 2;
    rebinEtaPhiX[1] = 2; rebinEtaPhiY[1] = 1;
    rebinEtaPhiX[2] = 5; rebinEtaPhiY[2] = 1;
    rebinEtaPhiX[3] = 1; rebinEtaPhiY[3] = 5;
    rebinEtaPhiX[4] = 1; rebinEtaPhiY[4] = 1; // dummy
  }
  AliTRDrecoProjection hp[nprj]; TObjArray php(nprj+2);
  Int_t ih(0), isel(-1), np[nsel]; memset(np, 0, nsel*sizeof(Int_t));
  // define list of projections
  for(Int_t ipt(0); ipt<nPt; ipt++){
    ptName[ipt]= nPt>5?(64+ipt):ptShortName[ipt];
    ptCut[ipt] = StrDup(Form("#it{%4.2f<=p_{t}^{%s}[GeV/c]<%4.2f}",ipt?fgPt[ipt-1]:0., mc?"MC":"", ipt>fgNPt?99.99:fgPt[ipt]));
    pCut[ipt]  = StrDup(Form("#it{%4.2f<=p^{%s}[GeV/c]<%4.2f}",ipt?fgPt[ipt-1]:0., mc?"MC":"", ipt>fgNPt?99.99:fgPt[ipt]));
    for(Int_t isp(0); isp<kNspc; isp++){
      for(Int_t ich(0); ich<kNcharge; ich++){
        isel++; // new selection
        AliDebug(3, Form("Building sel[%3d|%4d] spc[%s%c] pt[%s]", isel, ih, spcName[v0][isp], chSgn[ich], ptCut[ipt]));
        hp[ih].Build(Form("H%sTrkInY%c%c%d", prefix, chName[ich], ptName[ipt], isp),
                     Form("TrackIn[%s%c]:: #Deltay{%s}", spcName[v0][isp], chSgn[ich], ptCut[ipt]),
                     kEta, kPhi, kYrez, aa);
        hp[ih].SetRebinStrategy(nEtaPhi, rebinEtaPhiX, rebinEtaPhiY);
        php.AddLast(&hp[ih++]); np[isel]++;
        hp[ih].Build(Form("H%sTrkInPh%c%c%d", prefix, chName[ich], ptName[ipt], isp),
                     Form("TrackIn[%s%c]:: #Delta#phi{%s}", spcName[v0][isp], chSgn[ich], ptCut[ipt]),
                     kEta, kPhi, kPrez, aa);
        hp[ih].SetRebinStrategy(nEtaPhi, rebinEtaPhiX, rebinEtaPhiY);
        php.AddLast(&hp[ih++]); np[isel]++;
        hp[ih].Build(Form("H%sTrkInQ%c%c%d", prefix, chName[ich], ptName[ipt], isp),
                     Form("TrackIn[%s%c]:: dQdl {%s}", spcName[v0][isp], chSgn[ich], pCut[ipt]),
                     kEta, kPhi, kZrez, aa);
        hp[ih].SetRebinStrategy(nEtaPhi, rebinEtaPhiX, rebinEtaPhiY);
        php.AddLast(&hp[ih++]); np[isel]++;
        if(!ax) continue;
        hp[ih].Build(Form("H%sTrkInX%c%c%d", prefix, chName[ich], ptName[ipt], isp),
                     Form("TrackIn[%s%c]:: #Deltax{%s}", spcName[v0][isp], chSgn[ich], ptCut[ipt]),
                     kEta, kPhi, kNdim+1, aa);
        hp[ih].SetRebinStrategy(nEtaPhi, rebinEtaPhiX, rebinEtaPhiY);
        php.AddLast(&hp[ih++]); np[isel]++;
      }
    }
    isel++; // RC projections
    AliDebug(3, Form("Building RCsel[%3d] pt[%s]", isel, ptCut[ipt]));
    hp[ih].Build(Form("H%sTrkInRCZ%c", prefix, ptName[ipt]),
                  Form("TrackIn[RC]:: #Deltaz{%s}", ptCut[ipt]),
                  kEta, kPhi, kZrez, aa);
    hp[ih].SetRebinStrategy(nEtaPhi, rebinEtaPhiX, rebinEtaPhiY);
    php.AddLast(&hp[ih++]); np[isel]++;
    hp[ih].Build(Form("H%sTrkInRCY%c", prefix, ptName[ipt]),
                  Form("TrackIn[RC]:: #Deltay{%s}", ptCut[ipt]),
                  kEta, kPhi, kYrez, aa);
    hp[ih].SetRebinStrategy(nEtaPhi, rebinEtaPhiX, rebinEtaPhiY);
    php.AddLast(&hp[ih++]); np[isel]++;
    hp[ih].Build(Form("H%sTrkInRCPh%c", prefix, ptName[ipt]),
                  Form("TrackIn[RC]:: #Delta#phi{%s}", ptCut[ipt]),
                  kEta, kPhi, kPrez, aa);
    hp[ih].SetRebinStrategy(nEtaPhi, rebinEtaPhiX, rebinEtaPhiY);
    php.AddLast(&hp[ih++]); np[isel]++;
    if(!ax) continue;
    hp[ih].Build(Form("H%sTrkInRCX%c", prefix, ptName[ipt]),
                  Form("TrackIn[RC]:: #Deltax{%s}", ptCut[ipt]),
                  kEta, kPhi, kNdim+1, aa);
    hp[ih].SetRebinStrategy(nEtaPhi, rebinEtaPhiX, rebinEtaPhiY);
    php.AddLast(&hp[ih++]); np[isel]++;
  }
  // pt projections
  for(Int_t isp(0); isp<kNspc; isp++){
    for(Int_t ich(0); ich<kNcharge; ich++){
      isel++; // new selection
      AliDebug(3, Form("Building PTsel[%3d] spc[%s%c]", isel, spcName[v0][isp], chSgn[ich]));
      hp[ih].Build(Form("H%sTrkInPt%c%d", prefix, chName[ich], isp),
                    Form("TrackIn[%s%c]:: P_{t}[GeV/c]", spcName[v0][isp], chSgn[ich]),
                    kEta, kPhi, kPt, aa);
      hp[ih].SetRebinStrategy(nEtaPhi, rebinEtaPhiX, rebinEtaPhiY);
      php.AddLast(&hp[ih++]); //np[isel]++;
    }
  }
  AliInfo(Form("Build %3d 3D %s projections.", ih, prefix));

  // fill projections
  Int_t ch(0), pt(0), p(0), sp(1), rcBin(as?as->FindBin(0.):-1), ioff(0), joff(0), jsel(0), ksel(0);
  AliTRDrecoProjection *pr0(NULL), *pr1(NULL);
  for (Long64_t ib(0); ib < H->GetNbins(); ib++) {
    v = H->GetBinContent(ib, coord);
    if(v<1.) continue;
    if(fBCbinTOF>0 && coord[kBC]!=fBCbinTOF) continue; // TOF bunch cross cut
    if(fBCbinFill>0 && abf && coord[kNdim+2]!=fBCbinTOF) continue; // Fill bunch cut
    if(coord[kBC]==3) continue;
    // charge selection
    ch = 0; sp=1;// [pi-] track
    if(rcBin>0){
      if(v0){ // V0 PID for dQdl
        if(!coord[kBC]) continue;
        sp   = coord[kBC]-1;
      } else {// TPC PID for statistics
        sp = Int_t(TMath::Abs(as->GetBinCenter(coord[kSpeciesChgRC])))-1;
      }

      // take care of old data format (2*AliPID::kSPECIES+1)
      if(as->GetNbins() == kNcharge*AliPID::kSPECIES+1){
        if(sp>2) sp=2;
        else if(sp==2) sp=1;
      }
      if(coord[kSpeciesChgRC] > rcBin) ch = 1;  // [+] track
      else if(coord[kSpeciesChgRC] == rcBin) ch = 2;  // [RC] track
    }
    // pt selection
    pt = 0; p = 0; // low pt
    if(apt) pt = coord[kPt];//TMath::Min(coord[kPt], Int_t(kNpt)+1);
    if(ap ) p  = coord[kNdim];//TMath::Min(coord[kNdim], Int_t(kNpt)+1);
    // global selection
    isel = pt*(kNspc*kNcharge+1); isel+=(ch==2?(kNspc*kNcharge):(sp*kNcharge+ch)); ioff = isel*(ax?4:3);
    jsel = p *(kNspc*kNcharge+1); jsel+=(ch==2?(kNspc*kNcharge):(sp*kNcharge+ch)); joff = jsel*(ax?4:3);
    ksel = ch==2?0:nPt*(kNspc*kNcharge+1); ksel *= (ax?4:3); ksel+=sp*kNcharge+ch;
    if(ioff>=ih){
      AliError(Form("Wrong selection %d [%3d]", ioff, ih));
      return kFALSE;
    }
    if(!(pr0=(AliTRDrecoProjection*)php.At(ioff))) {
      AliError(Form("Missing projection @ %d", ioff));
      return kFALSE;
    }
    if(ch<2){
      if(strcmp(pr0->H()->GetName(), Form("H%sTrkInY%c%c%d", prefix, chName[ch], ptName[pt], sp))!=0){
        AliError(Form("Projection mismatch :: request[H%sTrkInY%c%c%d] found[%s]", prefix, chName[ch], ptName[pt], sp, pr0->H()->GetName()));
        return kFALSE;
      }
    } else {
      if(strcmp(pr0->H()->GetName(), Form("H%sTrkInRCZ%c", prefix, ptName[pt]))!=0){
        AliError(Form("Projection mismatch :: request[H%sTrkltRCZ%c] found[%s]", prefix, ptName[pt], pr0->H()->GetName()));
        return kFALSE;
      }
    }
    AliDebug(2, Form("Found %s for selection sp[%d] ch[%d] pt[%d]", pr0->H()->GetName(), sp, ch, pt));
    for(Int_t jh(0); jh<np[isel]; jh++){
      if(ch<2 && jh==2) ((AliTRDrecoProjection*)php.At(joff+jh))->Increment(coord, v); // special care for dQdl=f(p)
      else ((AliTRDrecoProjection*)php.At(ioff+jh))->Increment(coord, v); // all = f(p_t)
    }
    if(ksel){
      if(!(pr0 = (AliTRDrecoProjection*)php.At(ksel))){
        AliError(Form("Missing P_t projection @ %d", ksel));
        return kFALSE;
      }
      AliDebug(2, Form("Found %s for selection[%d] sp[%d] ch[%d]", pr0->H()->GetName(), ksel, sp, ch));
      pr0->Increment(coord, v); // p_t spectra
    }
  }
  if(HasDump3DFor(kTrackIn)){
    TDirectory *cwd = gDirectory;
    TFile::Open(Form("Dump%s_%s.root", GetName(), H->GetName()), "UPDATE");
    for(Int_t ip(0); ip<php.GetEntriesFast(); ip++){
      if(!(pr0 = (AliTRDrecoProjection*)php.At(ip))) continue;
      if(!pr0->H()) continue;
      TH3 *h3=(TH3*)pr0->H()->Clone();
      h3->Write();
    }
    gFile->Close();
    cwd->cd();
  }

  TObjArray *arr(NULL);
  fProj->AddAt(arr = new TObjArray(mc?kMCTrkInNproj:kTrkInNproj), cidx);
  arr->SetName(mc?projName[1]:(v0?projName[2]:projName[0]));  arr->SetOwner();

  TH2 *h2(NULL); Int_t jh(0);
  for(; ih--; ){
    if(!hp[ih].H()) continue;
    if(strstr(hp[ih].H()->GetName(), "TrkInPt")){
      for(Int_t ipt(0); ipt<nPt; ipt++) arr->AddAt(hp[ih].Projection2Dbin(ipt), jh++);
    } else {
      if((h2 = hp[ih].Projection2D(kNstat, kNcontours))) arr->AddAt(h2, jh++);
    }
  }
  // build combined performance plots
  // combine up the tree of projections
  Double_t m(0.), s(0.), se(0.), trend(0.);
  AliTRDrecoProjection xlow[2], specY[kNcharge*kNspc], specPh[kNcharge*kNspc], specQ[kNcharge*kNspc];
  for(Int_t ich(0); ich<kNcharge; ich++){
    // PID dependency - summation over pt
    for(Int_t isp(0); isp<kNspc; isp++){
      /*!dy*/
      Int_t idx(ich*kNspc+isp);
      if((pr0 = (AliTRDrecoProjection*)php.FindObject(Form("H%sTrkInY%c%c%d", prefix, chName[ich], ptName[0], isp)))){
        specY[idx] = (*pr0);
        specY[idx].SetNameTitle(Form("H%sTrkInY%c%d", prefix, chName[ich], isp), "Sum over pt");
        specY[idx].H()->SetNameTitle(Form("H%sTrkInY%c%d", prefix, chName[ich], isp),
                              Form("TrackIn[%s%c]:: #Deltay", spcName[v0][isp], chSgn[ich]));
        if((trend=pr0->GetTrendValue(1,&m,&s,&se))>-100.){
          PutTrendValue(Form("%sTrkInY%c%c%d", prefix, chName[ich], ptName[0], isp), trend, m);
          PutTrendValue(Form("%sTrkInYS%c%c%d", prefix, chName[ich], ptName[0], isp), s, se);
        }
        for(Int_t ipt(1); ipt<nPt; ipt++){
          if(!(pr1 = (AliTRDrecoProjection*)php.FindObject(Form("H%sTrkInY%c%c%d", prefix, chName[ich], ptName[ipt], isp)))) continue;
          if((trend=pr1->GetTrendValue(1,&m,&s, &se))>-100.){
            PutTrendValue(Form("%sTrkInY%c%c%d", prefix, chName[ich], ptName[ipt], isp), trend, m);
            PutTrendValue(Form("%sTrkInYS%c%c%d", prefix, chName[ich], ptName[ipt], isp), s, se);
          }
          specY[idx]+=(*pr1);
        }
        php.AddLast(&specY[idx]);
        if((h2 = specY[idx].Projection2D(kNstat, kNcontours, 1, kFALSE))) arr->AddAt(h2, jh++);
        if((trend=specY[idx].GetTrendValue(1,&m,&s,&se))>-100.){
          PutTrendValue(Form("%sTrkInY%c%d", prefix, chName[ich], isp), trend, m);
          PutTrendValue(Form("%sTrkInYS%c%d", prefix, chName[ich], isp), s,se);
        }
        if((h2 = (TH2*)gDirectory->Get(Form("%sEn", specY[idx].H()->GetName())))) arr->AddAt(h2, jh++);
        if(ich && (pr1 = (AliTRDrecoProjection*)php.FindObject(Form("H%sTrkInY%c%d", prefix, chName[0], isp)))) (*pr1)+=specY[idx];
      }
      /*!dphi*/
      if((pr0 = (AliTRDrecoProjection*)php.FindObject(Form("H%sTrkInPh%c%c%d", prefix, chName[ich], ptName[0], isp)))){
        specPh[idx] = (*pr0);
        specPh[idx].SetNameTitle(Form("H%sTrkInPh%c%d", prefix, chName[ich], isp), "Sum over pt");
        specPh[idx].H()->SetNameTitle(Form("H%sTrkInPh%c%d", prefix, chName[ich], isp),
                              Form("TrackIn[%s%c]:: #Delta#phi", spcName[v0][isp], chSgn[ich]));
        specPh[idx].SetShowRange(-1.5, 1.5);
        if((trend=pr0->GetTrendValue(1,&m,&s,&se))>-100.){
          PutTrendValue(Form("%sTrkInPh%c%c%d", prefix, chName[ich], ptName[0], isp), trend, m);
          PutTrendValue(Form("%sTrkInPhS%c%c%d", prefix, chName[ich], ptName[0], isp), s, se);
        }
        for(Int_t ipt(1); ipt<nPt; ipt++){
          if(!(pr1 = (AliTRDrecoProjection*)php.FindObject(Form("H%sTrkInPh%c%c%d", prefix, chName[ich], ptName[ipt], isp)))) continue;
          if((trend=pr1->GetTrendValue(1,&m,&s,&se))>-100.){
            PutTrendValue(Form("%sTrkInPh%c%c%d", prefix, chName[ich], ptName[ipt], isp), trend, m);
            PutTrendValue(Form("%sTrkInPhS%c%c%d", prefix, chName[ich], ptName[ipt], isp), s, se);
          }
          specPh[idx]+=(*pr1);
        }
        php.AddLast(&specPh[idx]);
        if((h2 = specPh[idx].Projection2D(kNstat, kNcontours, 1))) arr->AddAt(h2, jh++);
        if((trend=specPh[idx].GetTrendValue(1,&m,&s,&se))>-100.){
          PutTrendValue(Form("%sTrkInPh%c%d", prefix, chName[ich], isp), trend, m);
          PutTrendValue(Form("%sTrkInPhS%c%d", prefix, chName[ich], isp), s,se);
        }
        if(ich && (pr1 = (AliTRDrecoProjection*)php.FindObject(Form("H%sTrkInPh%c%d", prefix, chName[0], isp)))) (*pr1)+=specPh[idx];
      }
      /*!dQdl*/
      if((pr0 = (AliTRDrecoProjection*)php.FindObject(Form("H%sTrkInQ%c%c%d", prefix, chName[ich], ptName[0], isp)))){
        specQ[idx] = (*pr0);
        specQ[idx].SetNameTitle(Form("H%sTrkInQ%c%d", prefix, chName[ich], isp), "Sum over p");
        specQ[idx].H()->SetNameTitle(Form("H%sTrkInQ%c%d", prefix, chName[ich], isp),
                              Form("TrackIn[%s%c]:: dQdl", spcName[v0][isp], chSgn[ich]));
        specQ[idx].SetShowRange(-2.2, -1.75);
        specQ[idx].H()->GetZaxis()->SetTitle("dQdl [a.u.]");
        if((trend = pr0->GetTrendValue(2, &m))>-100.) PutTrendValue(Form("%sTrkInQ%c%c%d", prefix, chName[ich], ptName[0], isp), trend, m);
        if((trend = pr0->GetTrendValue(0, &m))>-100.) PutTrendValue(Form("%sTrkInQS%c%c%d", prefix, chName[ich], ptName[0], isp), trend, m);
        for(Int_t ipt(1); ipt<nPt; ipt++){
          if(!(pr1 = (AliTRDrecoProjection*)php.FindObject(Form("H%sTrkInQ%c%c%d", prefix, chName[ich], ptName[ipt], isp)))) continue;
          if((trend=pr1->GetTrendValue(2, &m))>-100.) PutTrendValue(Form("%sTrkInQ%c%c%d", prefix, chName[ich], ptName[ipt], isp), trend, m);
          if((trend=pr1->GetTrendValue(0, &m))>-100.) PutTrendValue(Form("%sTrkInQS%c%c%d", prefix, chName[ich], ptName[ipt], isp), trend, m);
          specQ[idx]+=(*pr1);
        }
        php.AddLast(&specQ[idx]);
        if((h2 = specQ[idx].Projection2D(kNstat, kNcontours, 2))) arr->AddAt(h2, jh++);
        specQ[idx].H()->SetName(Form("H%sTrkInQS%c%d", prefix, chName[ich], isp));
        specQ[idx].SetShowRange(-1.85, -1.4);
        if((h2 = specQ[idx].Projection2D(kNstat, kNcontours, 0))) arr->AddAt(h2, jh++);
        if((trend=specQ[idx].GetTrendValue(2, &m))>-100.) PutTrendValue(Form("%sTrkInQ%c%d", prefix, chName[ich], isp), trend, m);
        if((trend=specQ[idx].GetTrendValue(0, &m))>-100.) PutTrendValue(Form("%sTrkInQS%c%d", prefix, chName[ich], isp), trend, m);
        if(ich && (pr1 = (AliTRDrecoProjection*)php.FindObject(Form("H%sTrkInQ%c%d", prefix, chName[0], isp)))) (*pr1)+=specQ[idx];
      }
    } // end PID loop for pt integration

    // pt dependency - summation over PID
    for(Int_t ipt(0); ipt<nPt; ipt++){
      /*!dy*/
      if((pr0 = (AliTRDrecoProjection*)php.FindObject(Form("H%sTrkInY%c%c%d", prefix, chName[ich], ptName[ipt], 0)))){
        for(Int_t isp(1); isp<kNspc; isp++){
          if(!(pr1 = (AliTRDrecoProjection*)php.FindObject(Form("H%sTrkInY%c%c%d", prefix, chName[ich], ptName[ipt], isp)))) continue;
          (*pr0)+=(*pr1);
        }
        pr0->H()->SetNameTitle(Form("H%sTrkInY%c%c", prefix, chName[ich], ptName[ipt]),
                                  Form("TrackIn[%c]:: #Deltay{%s}", chSgn[ich], ptCut[ipt]));
        pr0->SetShowRange(-0.3, 0.3);
        if((h2 = pr0->Projection2D(kNstat, kNcontours, 1))) arr->AddAt(h2, jh++);
        if((trend=pr0->GetTrendValue(1,&m,&s,&se))>-100.){
          PutTrendValue(Form("%sTrkInY%c%c", prefix, chName[ich], ptName[ipt]), trend, m);
          PutTrendValue(Form("%sTrkInYS%c%c", prefix, chName[ich], ptName[ipt]), s, se);
        }
        if(ipt && (pr1 = (AliTRDrecoProjection*)php.FindObject(Form("H%sTrkInY%c%c%d", prefix, chName[ich], ptName[0], 0)))) (*pr1)+=(*pr0);
      }
      /*!dphi*/
      if((pr0 = (AliTRDrecoProjection*)php.FindObject(Form("H%sTrkInPh%c%c%d", prefix, chName[ich], ptName[ipt], 0)))){
        for(Int_t isp(1); isp<kNspc; isp++){
          if(!(pr1 = (AliTRDrecoProjection*)php.FindObject(Form("H%sTrkInPh%c%c%d", prefix, chName[ich], ptName[ipt], isp)))) continue;
          (*pr0)+=(*pr1);
        }
        pr0->H()->SetNameTitle(Form("H%sTrkInPh%c%c", prefix, chName[ich], ptName[ipt]),
                                  Form("TrackIn[%c]:: #Delta#phi{%s}", chSgn[ich], ptCut[ipt]));
        if((h2 = pr0->Projection2D(kNstat, kNcontours, 1))) arr->AddAt(h2, jh++);
        if((trend=pr0->GetTrendValue(1,&m,&s,&se))>-100.){
          PutTrendValue(Form("%sTrkInPh%c%c", prefix, chName[ich], ptName[ipt]), trend, m);
          PutTrendValue(Form("%sTrkInPhS%c%c", prefix, chName[ich], ptName[ipt]), s, se);
        }
        if(ipt && (pr1 = (AliTRDrecoProjection*)php.FindObject(Form("H%sTrkInPh%c%c%d", prefix, chName[ich], ptName[0], 0)))) (*pr1)+=(*pr0);
      }
      /*!dx*/
      if((pr0 = (AliTRDrecoProjection*)php.FindObject(Form("H%sTrkInX%c%c%d", prefix, chName[ich], ptName[ipt], 0)))){
        for(Int_t isp(1); isp<kNspc; isp++){
          if(!(pr1 = (AliTRDrecoProjection*)php.FindObject(Form("H%sTrkInX%c%c%d", prefix, chName[ich], ptName[ipt], isp)))) continue;
          (*pr0)+=(*pr1);
        }
        pr0->H()->SetNameTitle(Form("H%sTrkInX%c%c", prefix, chName[ich], ptName[ipt]),
                                  Form("TrackIn[%c]:: #Deltax{%s}", chSgn[ich], ptCut[ipt]));
        if((h2 = pr0->Projection2D(kNstat, kNcontours, 1))) arr->AddAt(h2, jh++);
        PutTrendValue(Form("%sTrkInX%c%c", prefix, chName[ich], ptName[ipt]), pr0->GetTrendValue(1));
        if(!ipt){
          xlow[ich] = (*pr0);
          xlow[ich].SetNameTitle(Form("H%sTrkInX%c%c%d", prefix, chName[ich], ptName[0], 5),
                                 Form("TrackIn[%c]:: #Deltax{%s}", chSgn[ich], ptCut[0]));
          php.AddLast(&xlow[ich]);
        }
        if(ipt && (pr1 = (AliTRDrecoProjection*)php.FindObject(Form("H%sTrkInX%c%c%d", prefix, chName[ich], ptName[0], 0)))) (*pr1)+=(*pr0);
      }
    } // end pt loop for PID integration

    /*!dy*/
    if((pr0 = (AliTRDrecoProjection*)php.FindObject(Form("H%sTrkInY%c%c%d", prefix, chName[ich], ptName[0], 0)))){
      pr0->H()->SetNameTitle(Form("H%sTrkInY%c", prefix, chName[ich]),
                            Form("TrackIn[%c]:: #Deltay", chSgn[ich]));
      pr0->SetShowRange(-0.3, 0.3);
      if((h2 = pr0->Projection2D(kNstat, kNcontours, 1, kFALSE))) arr->AddAt(h2, jh++);
      if((h2 = (TH2*)gDirectory->Get(Form("%sEn", pr0->H()->GetName())))) arr->AddAt(h2, jh++);
      if((trend=pr0->GetTrendValue(1, &m, &s, &se))>-100.){
        PutTrendValue(Form("%sTrkInY%c", prefix, chName[ich]), trend, m);
        PutTrendValue(Form("%sTrkInYS%c", prefix, chName[ich]), s, se);
      }
      if(ich && (pr1 = (AliTRDrecoProjection*)php.FindObject(Form("H%sTrkInY%c%c%d", prefix, chName[0], ptName[0], 0)))) (*pr1)+=(*pr0);
    }
    /*!dy high pt*/
    if(ich && (pr0 = (AliTRDrecoProjection*)php.FindObject(Form("H%sTrkInY%c%c%d", prefix, chName[0], ptName[3], 0)))){
      if((pr1 = (AliTRDrecoProjection*)php.FindObject(Form("H%sTrkInY%c%c%d", prefix, chName[ich], ptName[3], 0)))){
        (*pr0)+=(*pr1);
        pr0->H()->SetNameTitle(Form("H%sTrkInY%c", prefix, ptName[3]), Form("TrackIn :: #Deltay{%s}", ptCut[3]));
        pr0->SetShowRange(-0.3, 0.3);
        if((h2 = pr0->Projection2D(kNstat, kNcontours, 1))) arr->AddAt(h2, jh++);
      }
    }
    /*!dphi*/
    if((pr0 = (AliTRDrecoProjection*)php.FindObject(Form("H%sTrkInPh%c%c%d", prefix, chName[ich], ptName[0], 0)))){
      pr0->H()->SetNameTitle(Form("H%sTrkInPh%c", prefix, chName[ich]),
                            Form("TrackIn[%c]:: #Delta#phi", chSgn[ich]));
      pr0->SetShowRange(-1., 1.);
      if((h2 = pr0->Projection2D(kNstat, kNcontours, 1))) arr->AddAt(h2, jh++);
      if((trend=pr0->GetTrendValue(1, &m, &s, &se))>-100.){
        PutTrendValue(Form("%sTrkInPh%c", prefix, chName[ich]), trend, m);
        PutTrendValue(Form("%sTrkInPhS%c", prefix, chName[ich]), s, se);
      }
      if(ich==1 && (pr1 = (AliTRDrecoProjection*)php.FindObject(Form("H%sTrkInPh%c%c%d", prefix, chName[0], ptName[0], 0)))) (*pr1)+=(*pr0);
    }
    /*!dx*/
    if((pr0 = (AliTRDrecoProjection*)php.FindObject(Form("H%sTrkInX%c%c%d", prefix, chName[ich], ptName[0], 0)))){
      pr0->H()->SetNameTitle(Form("H%sTrkInX%c", prefix, chName[ich]),
                            Form("TrackIn[%c]:: #Deltax", chSgn[ich]));
      if((h2 = pr0->Projection2D(kNstat, kNcontours, 1))) arr->AddAt(h2, jh++);
      PutTrendValue(Form("%sTrkInX%c", prefix, chName[ich]), pr0->GetTrendValue(1));
      if(ich==1 && (pr1 = (AliTRDrecoProjection*)php.FindObject(Form("H%sTrkInX%c%c%d", prefix, chName[0], ptName[0], 0)))) (*pr1)+=(*pr0);
    }
    /*!dx low pt*/
    if(ich && (pr0 = (AliTRDrecoProjection*)php.FindObject(Form("H%sTrkInX%c%c%d", prefix, chName[0], ptName[1], 5)))){
      if((pr1 = (AliTRDrecoProjection*)php.FindObject(Form("H%sTrkInX%c%c%d", prefix, chName[ich], ptName[1], 5)))){
        (*pr0)+=(*pr1);
        pr0->H()->SetNameTitle(Form("H%sTrkInX%c", prefix, ptName[1]), Form("TrackIn :: #Deltax{%s}", ptCut[1]));
        if((h2 = pr0->Projection2D(kNstat, kNcontours, 1))) arr->AddAt(h2, jh++);
      }
    }
  } // end charge loop

  for(Int_t isp(0); isp<kNspc; isp++){
    /*!dy*/
    if((pr0 = (AliTRDrecoProjection*)php.FindObject(Form("H%sTrkInY%c%d", prefix, chName[0], isp)))){
      pr0->H()->SetNameTitle(Form("H%sTrkInY%d", prefix, isp), Form("TrackIn[%s] :: #Deltay", spcName[v0][isp]));
      pr0->SetShowRange(-0.3, 0.3);
      if((h2 = pr0->Projection2D(kNstat, kNcontours, 1, kFALSE))) arr->AddAt(h2, jh++);
      PutTrendValue(Form("%sTrkInY%d", prefix, isp), pr0->GetTrendValue(1));
      if((h2 = (TH2*)gDirectory->Get(Form("%sEn", pr0->H()->GetName())))) arr->AddAt(h2, jh++);
    }
    /*!dphi*/
    if((pr0 = (AliTRDrecoProjection*)php.FindObject(Form("H%sTrkInPh%c%d", prefix, chName[0], isp)))){
      pr0->H()->SetNameTitle(Form("H%sTrkInPh%d", prefix, isp), Form("TrackIn[%s] :: #Delta#phi", spcName[v0][isp]));
      pr0->SetShowRange(-1., 1.);
      if((h2 = pr0->Projection2D(kNstat, kNcontours, 1))) arr->AddAt(h2, jh++);
      PutTrendValue(Form("%sTrkInPh%d", prefix, isp), pr0->GetTrendValue(1));
    }
    /*!dQdl*/
    if((pr0 = (AliTRDrecoProjection*)php.FindObject(Form("H%sTrkInQ%c%d", prefix, chName[0], isp)))){
      pr0->H()->SetNameTitle(Form("H%sTrkInQ%d", prefix, isp), Form("TrackIn[%s] :: dQdl", spcName[v0][isp]));
      pr0->SetShowRange(-2.2, -1.75);
      if((h2 = pr0->Projection2D(kNstat, kNcontours, 2))) arr->AddAt(h2, jh++);
      pr0->H()->SetName(Form("H%sTrkInQS%d", prefix, isp));
      pr0->SetShowRange(-1.85, -1.4);
      if((h2 = pr0->Projection2D(kNstat, kNcontours, 0))) arr->AddAt(h2, jh++);
      if((trend=pr0->GetTrendValue(2, &m))>-100.) PutTrendValue(Form("%sTrkInQ%d", prefix, isp), trend, m);
      if((trend=pr0->GetTrendValue(0, &m))>-100.) PutTrendValue(Form("%sTrkInQS%d", prefix, isp), trend, m);
    }
  } // end PID processing

  /*!dy*/
  if((pr0 = (AliTRDrecoProjection*)php.FindObject(Form("H%sTrkInY%c%c%d", prefix, chName[0], ptName[0], 0)))){
    pr0->H()->SetNameTitle(Form("H%sTrkInY", prefix), "TrackIn :: #Deltay");
    pr0->SetShowRange(-0.3, 0.3);
    if((h2 = pr0->Projection2D(kNstat, kNcontours, 1, kFALSE))) arr->AddAt(h2, jh++);
    if((h2 = (TH2*)gDirectory->Get(Form("%sEn", pr0->H()->GetName())))) arr->AddAt(h2, jh++);
  }
  /*!dphi*/
  if((pr0 = (AliTRDrecoProjection*)php.FindObject(Form("H%sTrkInPh%c%c%d", prefix, chName[0], ptName[0], 0)))){
    pr0->H()->SetNameTitle(Form("H%sTrkInPh", prefix), "TrackIn :: #Delta#phi");
    if((h2 = pr0->Projection2D(kNstat, kNcontours, 1))) arr->AddAt(h2, jh++);
  }
  /*!dx*/
  if((pr0 = (AliTRDrecoProjection*)php.FindObject(Form("H%sTrkInX%c%c%d", prefix, chName[0], ptName[0], 0)))){
    pr0->H()->SetNameTitle(Form("H%sTrkInX", prefix), "TrackIn :: #Deltax");
    if((h2 = pr0->Projection2D(kNstat, kNcontours, 1))) arr->AddAt(h2, jh++);
  }

  // Row Cross processing
  /*!RC dz*/
  if((pr0 = (AliTRDrecoProjection*)php.FindObject(Form("H%sTrkInRCZ%c", prefix, ptName[0])))){
    for(Int_t ipt(0); ipt<nPt; ipt++){
      if(!(pr1 = (AliTRDrecoProjection*)php.FindObject(Form("H%sTrkInRCZ%c", prefix, ptName[ipt])))) continue;
      (*pr0)+=(*pr1);
    }
    pr0->H()->SetNameTitle(Form("H%sTrkInRCZ", prefix), "TrackIn[RC]:: #Deltaz");
    if((h2 = pr0->Projection2D(kNstat, kNcontours, 1))) arr->AddAt(h2, jh++);
    if((trend=pr0->GetTrendValue(1,&m,&s,&se))>-100.){
      PutTrendValue(Form("%sTrkInRCZ", prefix), trend, m);
      PutTrendValue(Form("%sTrkInRCZS", prefix), s, se);
    }
  }
  /*!RC dy*/
  if((pr0 = (AliTRDrecoProjection*)php.FindObject(Form("H%sTrkInRCY%c", prefix, ptName[0])))){
    for(Int_t ipt(0); ipt<nPt; ipt++){
      if(!(pr1 = (AliTRDrecoProjection*)php.FindObject(Form("H%sTrkInRCY%c", prefix, ptName[ipt])))) continue;
      (*pr0)+=(*pr1);
    }
    pr0->H()->SetNameTitle(Form("H%sTrkInRCY", prefix), "TrackIn[RC]:: #Deltay");
    if((h2 = pr0->Projection2D(kNstat, kNcontours, 1))) arr->AddAt(h2, jh++);
  }
  /*!RC dphi*/
  if((pr0 = (AliTRDrecoProjection*)php.FindObject(Form("H%sTrkInRCPh%c", prefix, ptName[0])))){
    for(Int_t ipt(0); ipt<nPt; ipt++){
      if(!(pr1 = (AliTRDrecoProjection*)php.FindObject(Form("H%sTrkInRCPh%c", prefix, ptName[ipt])))) continue;
      (*pr0)+=(*pr1);
    }
    pr0->H()->SetNameTitle(Form("H%sTrkInRCPh", prefix), "TrackIn[RC]:: #Delta#phi");
    if((h2 = pr0->Projection2D(kNstat, kNcontours, 1))) arr->AddAt(h2, jh++);
  }
  /*!RC dx*/
  if((pr0 = (AliTRDrecoProjection*)php.FindObject(Form("H%sTrkInRCX%c", prefix, ptName[0])))){
    for(Int_t ipt(0); ipt<nPt; ipt++){
      if(!(pr1 = (AliTRDrecoProjection*)php.FindObject(Form("H%sTrkInRCX%c", prefix, ptName[ipt])))) continue;
      (*pr0)+=(*pr1);
    }
    pr0->H()->SetNameTitle(Form("H%sTrkInRCX", prefix), "TrackIn[RC]:: #Deltax");
    if((h2 = pr0->Projection2D(kNstat, kNcontours, 1))) arr->AddAt(h2, jh++);
  }

  // P_t processing
  for(Int_t isp(0); isp<kNspc; isp++){
    for(Int_t ich(0); ich<kNcharge; ich++){
      TH2 *h2pt[26]={NULL}, *h2s(NULL);
      for(Int_t ipt(0); ipt<nPt; ipt++){
        if(!(h2pt[ipt] = (TH2*)arr->FindObject(Form("H%sTrkInPt%c%d%d_2D", prefix, chName[ich], isp, ipt)))){
          AliWarning(Form("Missing \"H%sTrkInPt%c%d%d_2D\"", prefix, chName[ich], isp, ipt));
          continue;
        }
        if(!h2s) h2s = (TH2F*)h2pt[ipt]->Clone("h2s");
        else h2s->Add(h2pt[ipt]);
      }
      if(!h2s) continue;
      Int_t irebin(0), dxBin(1), dyBin(1);
      while(irebin<nEtaPhi && (AliTRDrecoTask::GetMeanStat(h2s, .5, 1)<200)){
        h2s->Rebin2D(rebinEtaPhiX[irebin], rebinEtaPhiY[irebin]);
        dxBin*=rebinEtaPhiX[irebin];dyBin*=rebinEtaPhiY[irebin];irebin++;
      }
      AliDebug(2, Form("Rebin level[%d] @ chg[%d] spc[%d]. Binning[%2dx%2d]", irebin, ich, isp, h2s->GetNbinsX(), h2s->GetNbinsY()));
      Int_t nx(h2s->GetNbinsX()), ny(h2s->GetNbinsY());
      for(Int_t ipt(0); ipt<nPt; ipt++) h2pt[ipt]->Rebin2D(dxBin, dyBin);
      delete h2s;

      h2 = new TH2F(Form("H%sTrkInPt%c%d_2D", prefix, chName[ich], isp),
                    Form("TrackIn[%s%c]:: <P_{t}>[GeV/c];%s;%s", spcName[v0][isp], chSgn[ich], aa[2]->GetTitle(), aa[1]->GetTitle()),
                    nx, aa[2]->GetXmin(), aa[2]->GetXmax(),
                    ny, aa[1]->GetXmin(), aa[1]->GetXmax());
      arr->AddAt(h2, jh++);
      for(Int_t ix(1); ix<=nx; ix++){ // eta
        for(Int_t iy(1); iy<=ny; iy++){ // phi
          Float_t w[26]={0.}, sw(0.);
          for(Int_t ipt(0); ipt<nPt; ipt++){
            if(!h2pt[ipt]) continue;
            w[ipt] = h2pt[ipt]->GetBinContent(ix, iy);
            sw    += w[ipt];
          }
          if(sw<=0.) h2->SetBinContent(ix, iy, 0.);
          else{
            Float_t ptm(0.);
            for(Int_t ipt(0); ipt<nPt; ipt++){
              w[ipt]/=sw;
              ptm   += w[ipt]*fgPt[ipt?ipt-1:0];
            }
            h2->SetBinContent(ix, iy, ptm);
          }
        }
      }
      PutTrendValue(Form("%sTrkInPt%c%d", prefix, chName[ich], isp), GetMeanStat(h2, 0.01, 1));
    }
  }
  AliInfo(Form("Done %3d 2D %s projections.", jh, prefix));

// clean local memory allocation
  for(Int_t i(fgNPt+2);i--;){
    if(ptCut[i]) delete [] ptCut[i];
    if(pCut[i]) delete [] pCut[i];
  }
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
    return kTRUE;
  }
  Int_t ndim(H->GetNdimensions());
  Int_t coord[kNdim+1]; memset(coord, 0, sizeof(Int_t) * (kNdim+1)); Double_t v = 0.;
  TAxis *aa[kNdim+1], *as(NULL), *ap(NULL); memset(aa, 0, sizeof(TAxis*) * (kNdim+1));
  for(Int_t id(0); id<ndim; id++) aa[id] = H->GetAxis(id);
  if(ndim > kSpeciesChgRC) as = H->GetAxis(kSpeciesChgRC);
  if(ndim > kPt) ap = H->GetAxis(kPt);
  const Int_t nPt(ap?(ap->GetNbins()+2):1);

  // build list of projections
  const Int_t nsel(AliTRDgeometry::kNlayer*fgNPt*AliPID::kSPECIES*7);//, npsel(3);
  const Char_t chName[kNcharge] = {'n', 'p'};const Char_t chSgn[kNcharge] = {'-', '+'};
  const Char_t ptShortName[5] = {'L', 'l', 'i', 'h', 'H'};
  Char_t ptName[fgNPt+2] = {0};
  Char_t *ptCut[fgNPt+2] = {NULL};
  // define rebinning strategy
  const Int_t nEtaPhi(5); Int_t rebinEtaPhiX[nEtaPhi] = {1, 3, 1, 3, 1}, rebinEtaPhiY[nEtaPhi] = {2, 1, 2, 1, 2};
  if(aa[1]->GetNbins()==180){ // old binning
    rebinEtaPhiX[0] = 1; rebinEtaPhiY[0] = 2;
    rebinEtaPhiX[1] = 2; rebinEtaPhiY[1] = 1;
    rebinEtaPhiX[2] = 5; rebinEtaPhiY[2] = 1;
    rebinEtaPhiX[3] = 1; rebinEtaPhiY[3] = 5;
    rebinEtaPhiX[4] = 1; rebinEtaPhiY[4] = 1; // dummy
  }
  AliTRDrecoProjection hp[kTrkNproj]; TObjArray php(kTrkNproj);
  Int_t ih(0), isel(-1), np[nsel]; memset(np, 0, nsel*sizeof(Int_t));
  for(Int_t ily(0); ily<AliTRDgeometry::kNlayer; ily++){
    for(Int_t ipt(0); ipt<nPt; ipt++){
      if(!ptName[ipt]){
        ptName[ipt]= nPt>5?(64+ipt):ptShortName[ipt];
        ptCut[ipt] = StrDup(Form("#it{%4.2f<=p_{t}^{MC}[GeV/c]<%4.2f}",ipt?fgPt[ipt-1]:0., ipt>fgNPt?99.99:fgPt[ipt]));
      }
      for(Int_t isp(0); isp<AliPID::kSPECIES; isp++){
        for(Int_t ich(0); ich<kNcharge; ich++){
          isel++; // new selection
          hp[ih].Build(Form("HMCTrkY%c%c%d%d", chName[ich], ptName[ipt], isp, UseLYselectTrklt()?fLYselect:ily),
                       Form("Tracks[%s%c]:: #Deltay{%s} Ly[%d]", AliPID::ParticleLatexName(isp), chSgn[ich], ptCut[ipt], UseLYselectTrklt()?fLYselect:ily),
                       kEta, kPhi, kYrez, aa);
          hp[ih].SetRebinStrategy(nEtaPhi, rebinEtaPhiX, rebinEtaPhiY);
          php.AddLast(&hp[ih++]); np[isel]++;
          hp[ih].Build(Form("HMCTrkPh%c%c%d%d", chName[ich], ptName[ipt], isp, UseLYselectTrklt()?fLYselect:ily),
                       Form("Tracks[%s%c]:: #Delta#phi{%s} Ly[%d]", AliPID::ParticleLatexName(isp), chSgn[ich], ptCut[ipt], UseLYselectTrklt()?fLYselect:ily),
                       kEta, kPhi, kPrez, aa);
          hp[ih].SetRebinStrategy(nEtaPhi, rebinEtaPhiX, rebinEtaPhiY);
          php.AddLast(&hp[ih++]); np[isel]++;
          hp[ih].Build(Form("HMCTrkDPt%c%c%d%d", chName[ich], ptName[ipt], isp, UseLYselectTrklt()?fLYselect:ily),
                       Form("Tracks[%s%c]:: #Deltap_{t}/p_{t}{%s} Ly[%d]", AliPID::ParticleLatexName(isp), chSgn[ich], ptCut[ipt], UseLYselectTrklt()?fLYselect:ily),
                       kEta, kPhi, kNdim, aa);
          hp[ih].SetShowRange(0.,10.);
          hp[ih].SetRebinStrategy(nEtaPhi, rebinEtaPhiX, rebinEtaPhiY);
          php.AddLast(&hp[ih++]); np[isel]++;
        }
      }
      isel++; // new selection
      hp[ih].Build(Form("HMCTrkZ%c%d", ptName[ipt], UseLYselectTrklt()?fLYselect:ily),
                    Form("Tracks[RC]:: #Deltaz{%s} Ly[%d]", ptCut[ipt], UseLYselectTrklt()?fLYselect:ily),
                    kEta, kPhi, kZrez, aa);
      hp[ih].SetRebinStrategy(nEtaPhi, rebinEtaPhiX, rebinEtaPhiY);
      php.AddLast(&hp[ih++]); np[isel]++;
    }
  }
  AliInfo(Form("Build %3d 3D projections.", ih));

  Int_t ly(0), ch(0), pt(0), sp(2), rcBin(as?as->FindBin(0.):-1);
  for (Long64_t ib(0); ib < H->GetNbins(); ib++) {
    v = H->GetBinContent(ib, coord);
    if(v<1.) continue;
    ly = coord[kBC]-1; // layer selection
    if(ly<0) ly=0; // fix for bug in PlotMC
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
    Int_t ioff = ly*nPt*31+pt*31; ioff+=3*(sp<0?10:(sp*kNcharge+ch));
    isel = ly*nPt*11+pt*11; isel+=sp<0?10:(sp*kNcharge+ch);
    AliDebug(4, Form("SELECTION[%d] :: ch[%c] pt[%c] sp[%d] ly[%d]\n", np[isel], ch==2?'Z':chName[ch], ptName[pt], sp, ly));
    for(Int_t jh(0); jh<np[isel]; jh++) ((AliTRDrecoProjection*)php.At(ioff+jh))->Increment(coord, v);
  }
  TObjArray *arr(NULL);
  fProj->AddAt(arr = new TObjArray(kTrkNproj), cidx);
  arr->SetName("hTRD2MC");  arr->SetOwner();

  TH2 *h2(NULL); Int_t jh(0);
  for(; ih--; ){
    if(!hp[ih].H()) continue;
    if(!(h2 = hp[ih].Projection2D(kNstat, kNcontours))) continue;
    arr->AddAt(h2, jh++);
  }

  // combine up the tree of projections
  AliTRDrecoProjection *pr0(NULL), *pr1(NULL);
  //Int_t iproj(0), jproj(0);
  for(Int_t ily(0); ily<AliTRDgeometry::kNlayer; ily++){
    for(Int_t ich(0); ich<kNcharge; ich++){
      for(Int_t ipt(0); ipt<nPt; ipt++){
        /*!dy*/
        if((pr0 = (AliTRDrecoProjection*)php.FindObject(Form("HMCTrkY%c%c%d%d", chName[ich], ptName[ipt], 0, UseLYselectTrklt()?fLYselect:ily)))){
          for(Int_t isp(1); isp<AliPID::kSPECIES; isp++){
            if(!(pr1 = (AliTRDrecoProjection*)php.FindObject(Form("HMCTrkY%c%c%d%d", chName[ich], ptName[ipt], isp, UseLYselectTrklt()?fLYselect:ily)))) continue;
            (*pr0)+=(*pr1);
          }
          AliDebug(2, Form("Rename %s to HMCTrkY%c%c%d", pr0->H()->GetName(), chName[ich], ptName[ipt], UseLYselectTrklt()?fLYselect:ily));
          pr0->H()->SetNameTitle(Form("HMCTrkY%c%c%d", chName[ich], ptName[ipt], UseLYselectTrklt()?fLYselect:ily),
                                    Form("Tracks[%c]:: #Deltay{%s} Ly[%d]", chSgn[ich], ptCut[ipt], UseLYselectTrklt()?fLYselect:ily));
          if((h2 = pr0->Projection2D(kNstat, kNcontours, 1))) arr->AddAt(h2, jh++);
          if(ipt && (pr1 = (AliTRDrecoProjection*)php.FindObject(Form("HMCTrkY%c%c%d%d", chName[ich], ptName[0], 0, UseLYselectTrklt()?fLYselect:ily)))) (*pr1)+=(*pr0);
        }
        /*!dphi*/
        if((pr0 = (AliTRDrecoProjection*)php.FindObject(Form("HMCTrkPh%c%c%d%d", chName[ich], ptName[ipt], 0, UseLYselectTrklt()?fLYselect:ily)))){
          for(Int_t isp(1); isp<AliPID::kSPECIES; isp++){
            if(!(pr1 = (AliTRDrecoProjection*)php.FindObject(Form("HMCTrkPh%c%c%d%d", chName[ich], ptName[ipt], isp, UseLYselectTrklt()?fLYselect:ily)))) continue;
            (*pr0)+=(*pr1);
          }
          AliDebug(2, Form("Rename %s to HMCTrkPh%c%c%d", pr0->H()->GetName(), chName[ich], ptName[ipt], UseLYselectTrklt()?fLYselect:ily));
          pr0->H()->SetNameTitle(Form("HMCTrkPh%c%c%d", chName[ich], ptName[ipt], UseLYselectTrklt()?fLYselect:ily),
                                    Form("Tracks[%c]:: #Delta#phi{%s} Ly[%d]", chSgn[ich], ptCut[ipt], UseLYselectTrklt()?fLYselect:ily));
          if((h2 = pr0->Projection2D(kNstat, kNcontours, 1))) arr->AddAt(h2, jh++);
          if(ipt && (pr1 = (AliTRDrecoProjection*)php.FindObject(Form("HMCTrkPh%c%c%d%d", chName[ich], ptName[0], 0, UseLYselectTrklt()?fLYselect:ily)))) (*pr1)+=(*pr0);
        }

        /*!dpt/pt*/
        if((pr0 = (AliTRDrecoProjection*)php.FindObject(Form("HMCTrkDPt%c%c%d%d", chName[ich], ptName[ipt], 0, UseLYselectTrklt()?fLYselect:ily)))){
          for(Int_t isp(1); isp<AliPID::kSPECIES; isp++){
            if(!(pr1 = (AliTRDrecoProjection*)php.FindObject(Form("HMCTrkDPt%c%c%d%d", chName[ich], ptName[ipt], isp, UseLYselectTrklt()?fLYselect:ily)))) continue;
            (*pr0)+=(*pr1);
          }
          AliDebug(2, Form("Rename %s to HMCTrkDPt%c%c%d", pr0->H()->GetName(), chName[ich], ptName[ipt], UseLYselectTrklt()?fLYselect:ily));
          pr0->H()->SetNameTitle(Form("HMCTrkDPt%c%c%d", chName[ich], ptName[ipt], UseLYselectTrklt()?fLYselect:ily),
                                    Form("Tracks[%c]:: #Deltap_{t}/p_{t}{%s} Ly[%d]", chSgn[ich], ptCut[ipt], UseLYselectTrklt()?fLYselect:ily));
          if((h2 = pr0->Projection2D(kNstat, kNcontours, 1))) arr->AddAt(h2, jh++);
          if(ipt && (pr1 = (AliTRDrecoProjection*)php.FindObject(Form("HMCTrkDPt%c%c%d%d", chName[ich], ptName[0], 0, UseLYselectTrklt()?fLYselect:ily)))) (*pr1)+=(*pr0);
        }
      }
      /*!dy*/
      if((pr0 = (AliTRDrecoProjection*)php.FindObject(Form("HMCTrkY%c%c%d%d", chName[ich], ptName[0], 0, UseLYselectTrklt()?fLYselect:ily)))){
        pr0->H()->SetNameTitle(Form("HMCTrkY%c%d", chName[ich], UseLYselectTrklt()?fLYselect:ily),
                              Form("Tracks[%c]:: #Deltay Ly[%d]", chSgn[ich], UseLYselectTrklt()?fLYselect:ily));
        if((h2 = pr0->Projection2D(kNstat, kNcontours, 1))) arr->AddAt(h2, jh++);
        if(ich==1 && (pr1 = (AliTRDrecoProjection*)php.FindObject(Form("HMCTrkY%c%c%d%d", chName[0], ptName[0], 0, UseLYselectTrklt()?fLYselect:ily)))) (*pr1)+=(*pr0);
      }
      /*!dphi*/
      if((pr0 = (AliTRDrecoProjection*)php.FindObject(Form("HMCTrkPh%c%c%d%d", chName[ich], ptName[0], 0, UseLYselectTrklt()?fLYselect:ily)))){
        pr0->H()->SetNameTitle(Form("HMCTrkPh%c%d", chName[ich], UseLYselectTrklt()?fLYselect:ily),
                              Form("Tracks[%c]:: #Delta#phi Ly[%d]", chSgn[ich], UseLYselectTrklt()?fLYselect:ily));
        if((h2 = pr0->Projection2D(kNstat, kNcontours, 1))) arr->AddAt(h2, jh++);
        if(ich==1 && (pr1 = (AliTRDrecoProjection*)php.FindObject(Form("HMCTrkPh%c%c%d%d", chName[0], ptName[0], 0, UseLYselectTrklt()?fLYselect:ily)))) (*pr1)+=(*pr0);
      }
      /*!dpt/pt*/
      if((pr0 = (AliTRDrecoProjection*)php.FindObject(Form("HMCTrkDPt%c%c%d%d", chName[ich], ptName[0], 0, UseLYselectTrklt()?fLYselect:ily)))){
        pr0->H()->SetNameTitle(Form("HMCTrkDPt%c%d", chName[ich], UseLYselectTrklt()?fLYselect:ily),
                              Form("Tracks[%c]:: #Deltap_{t}/p_{t} Ly[%d]", chSgn[ich], UseLYselectTrklt()?fLYselect:ily));
        if((h2 = pr0->Projection2D(kNstat, kNcontours, 1))) arr->AddAt(h2, jh++);
        if(ich==1 && (pr1 = (AliTRDrecoProjection*)php.FindObject(Form("HMCTrkDPt%c%c%d%d", chName[0], ptName[0], 0, UseLYselectTrklt()?fLYselect:ily)))) (*pr1)+=(*pr0);
      }
    }
    /*!dy*/
    if((pr0 = (AliTRDrecoProjection*)php.FindObject(Form("HMCTrkY%c%c%d%d", chName[0], ptName[0], 0, UseLYselectTrklt()?fLYselect:ily)))){
      pr0->H()->SetNameTitle(Form("HMCTrkY%d", UseLYselectTrklt()?fLYselect:ily), Form("Tracks :: #Deltay Ly[%d]", UseLYselectTrklt()?fLYselect:ily));
      if((h2 = pr0->Projection2D(kNstat, kNcontours, 1))) arr->AddAt(h2, jh++);
    }
    /*!dphi*/
    if((pr0 = (AliTRDrecoProjection*)php.FindObject(Form("HMCTrkPh%c%c%d%d", chName[0], ptName[0], 0, UseLYselectTrklt()?fLYselect:ily)))){
      pr0->H()->SetNameTitle(Form("HMCTrkPh%d", UseLYselectTrklt()?fLYselect:ily), Form("Tracks :: #Delta#phi Ly[%d]", UseLYselectTrklt()?fLYselect:ily));
      if((h2 = pr0->Projection2D(kNstat, kNcontours, 1))) arr->AddAt(h2, jh++);
    }
    /*!dpt/pt*/
    if((pr0 = (AliTRDrecoProjection*)php.FindObject(Form("HMCTrkDPt%c%c%d%d", chName[0], ptName[0], 0, UseLYselectTrklt()?fLYselect:ily)))){
      pr0->H()->SetNameTitle(Form("HMCTrkDPt%d", UseLYselectTrklt()?fLYselect:ily), Form("Tracks :: #Deltap_{t}/p_{t} Ly[%d]", UseLYselectTrklt()?fLYselect:ily));
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
    fProj = new TObjArray(kNclasses+1); fProj->SetOwner(kTRUE);
  }
  // set pt/p segmentation. guess from data
  THnSparse *H(NULL);
  if(!(H = (THnSparse*)fContainer->FindObject("hTracklet2TRDin"))){
    AliWarning("Missing/Wrong data @ hTracklet2TRDin needed for infering pt/p segmentation.");
    return kFALSE;
  }
  // protect for backward compatibility
  fNpt=H->GetAxis(kPt)?(H->GetAxis(kPt)->GetNbins()+1):1;
  if(!MakeMomSegmentation()) return kFALSE;

  //PROCESS EXPERIMENTAL DISTRIBUTIONS
  // Clusters detector
  if(HasProcess(kDetector)) MakeProjectionDetector();
  // Clusters residuals
  if(HasProcess(kCluster)) MakeProjectionCluster();
  fNRefFigures = 3;
  // Tracklet residual/pulls
  if(HasProcess(kTracklet)) MakeProjectionTracklet();
  fNRefFigures = 7;
  // TRDin residual/pulls
  if(HasProcess(kTrackIn)){
    MakeProjectionTrackIn();
    MakeProjectionTrackIn(kFALSE, kTRUE);
  }
  fNRefFigures = 11;

  if(!HasMCdata()) return kTRUE;
  //PROCESS MC RESIDUAL DISTRIBUTIONS

  // CLUSTER Y RESOLUTION/PULLS
  if(HasProcess(kCluster)) MakeProjectionCluster(kTRUE);
  fNRefFigures = 17;

  // TRACKLET RESOLUTION/PULLS
  if(HasProcess(kTracklet)) if(!MakeProjectionTracklet(kTRUE)) return kFALSE;
  fNRefFigures = 21;

  // TRACK RESOLUTION/PULLS
//  if(HasProcess(kTracklet)) if(!MakeProjectionTrack()) return kFALSE;
//  fNRefFigures+=16;

  // TRACK TRDin RESOLUTION/PULLS
  if(HasProcess(kTrackIn)) if(!MakeProjectionTrackIn(kTRUE)) return kFALSE;
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
TObjArray* AliTRDresolution::Histos()
{
  //
  // Define histograms
  //

  if(fContainer) return fContainer;

  fContainer  = new TObjArray(kNclasses); fContainer->SetOwner(kTRUE);
  THnSparse *H(NULL);
  const Int_t nhn(100); Char_t hn[nhn]; TString st;
  Int_t triggerBins(fTriggerList?(TMath::Power(2.,fTriggerList->GetEntriesFast())-1):0);

  //++++++++++++++++++++++
  // cluster to detector
  snprintf(hn, nhn, "h%s", fgPerformanceName[kDetector]);
  if(!(H = (THnSparseI*)gROOT->FindObject(hn))){
    const Int_t mdim(6);
    const Char_t *cldTitle[mdim] = {"layer", fgkTitle[kPhi], "pad row", "centrality", "q [a.u.]", triggerBins?"trigger":"no. of pads"/*, "tb [10 Hz]"*/};
    const Int_t cldNbins[mdim]   = {AliTRDgeometry::kNlayer, fgkNbins[kPhi], 76, AliTRDeventInfo::kCentralityClasses, 50, triggerBins?triggerBins:kNpads};
    const Double_t cldMin[mdim]  = {-0.5, fgkMin[kPhi], -0.5, -0.5, 0., 0.5},
                   cldMax[mdim]  = {AliTRDgeometry::kNlayer-0.5, fgkMax[kPhi], 75.5, AliTRDeventInfo::kCentralityClasses - 0.5, 1200., triggerBins?(triggerBins+.5):(kNpads+.5)};
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
    const Char_t *clTitle[mdim] = {"layer", fgkTitle[kPhi], fgkTitle[kEta], fgkTitle[kYrez], "#Deltax [cm]", "Q</Q", "Q/angle", "#Phi [deg]", "centrality", triggerBins?"trigger":"no. of pads"};
    const Int_t clNbins[mdim]   = {AliTRDgeometry::kNlayer, fgkNbins[kPhi], fgkNbins[kEta], fgkNbins[kYrez], 45, 30, 30, 15, AliTRDeventInfo::kCentralityClasses, triggerBins?triggerBins:kNpads};
    const Double_t clMin[mdim]  = {-0.5, fgkMin[kPhi], fgkMin[kEta], fgkMin[kYrez]/10., -.5, 0.1, -2., -45, -0.5, 0.5},
                   clMax[mdim]  = {AliTRDgeometry::kNlayer-0.5, fgkMax[kPhi], fgkMax[kEta], fgkMax[kYrez]/10., 4., 2.1, 118., 45, AliTRDeventInfo::kCentralityClasses - 0.5, triggerBins?(triggerBins+.5):(kNpads+.5)};
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
    trNbins[kPt]=fNpt-1; trMax[kPt] = trNbins[kPt]-.5;
    Int_t jdim(kNdim);
    trTitle[jdim]=StrDup("centrality"); trNbins[jdim] = AliTRDeventInfo::kCentralityClasses; trMin[jdim] = -.5; trMax[jdim] = AliTRDeventInfo::kCentralityClasses - 0.5;
    jdim++; 
    if(triggerBins){
      trTitle[jdim]=StrDup("trigger"); trNbins[jdim] = triggerBins;
      trMin[jdim] = 0.5; trMax[jdim] = triggerBins+.5;
    } else {
      trTitle[jdim]=StrDup("occupancy [%]"); trNbins[jdim] = 12;
      trMin[jdim] = 25.; trMax[jdim] = 85.;
    }

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
    const Int_t mdim(kNdim+2);
    Char_t *trinTitle[mdim]; memcpy(trinTitle, fgkTitle, kNdim*sizeof(Char_t*));
    Int_t trinNbins[mdim];   memcpy(trinNbins, fgkNbins, kNdim*sizeof(Int_t));
    Double_t trinMin[mdim];  memcpy(trinMin, fgkMin, kNdim*sizeof(Double_t));
    Double_t trinMax[mdim];  memcpy(trinMax, fgkMax, kNdim*sizeof(Double_t));
    if(triggerBins){
      trinTitle[kBC]=StrDup("trigger"); trinNbins[kBC] = triggerBins;
      trinMin[kBC] = 0.5; trinMax[kBC] = triggerBins+.5;
    }
//    trinNbins[kSpeciesChgRC] = Int_t(kNcharge)*(kNspc-1)+1; trinMin[kSpeciesChgRC] = -kNspc+0.5; trinMax[kSpeciesChgRC] = kNspc-0.5;
    trinNbins[kPt]=fNpt-1; trinMax[kPt] = trinNbins[kPt]-.5;
    trinTitle[kNdim]=StrDup("p [bin]"); trinNbins[kNdim] = fNpt-1; trinMin[kNdim] = -0.5; trinMax[kNdim] = trinNbins[kNdim]-.5;
    trinTitle[kNdim+1]=StrDup("dx [cm]"); trinNbins[kNdim+1]=48; trinMin[kNdim+1]=-2.4; trinMax[kNdim+1]=2.4;
    //trinTitle[kNdim+2]=StrDup("Fill Bunch"); trinNbins[kNdim+2]=3500; trinMin[kNdim+2]=-0.5; trinMax[kNdim+2]=3499.5;
    st = "r-#phi/z/angular residuals @ TRD entry;";
    // define minimum info to be saved in non debug mode
    Int_t ndim=(DebugLevel()>=1?mdim:(kNdim+1));//kNdimTrkIn;
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
    H = (THnSparseI*)((THnSparseI*)fContainer->At(kCluster))->Clone(hn);
    H->SetTitle("MC cluster spatial resolution");
  } else H->Reset();
  fContainer->AddAt(H, kMCcluster);
  //++++++++++++++++++++++
  // tracklet to TrackRef
  snprintf(hn, nhn, "h%s", fgPerformanceName[kMCtracklet]);
  if(!(H = (THnSparseI*)gROOT->FindObject(hn))){
    H = (THnSparseI*)((THnSparseI*)fContainer->At(kTracklet))->Clone(hn);
    H->SetTitle("MC tracklet spatial resolution");
  } else H->Reset();
  fContainer->AddAt(H, kMCtracklet);
  //++++++++++++++++++++++
  // TRDin to TrackRef
  snprintf(hn, nhn, "h%s", fgPerformanceName[kMCtrackIn]);
  if(!(H = (THnSparseI*)gROOT->FindObject(hn))){
    H = (THnSparseI*)((THnSparseI*)fContainer->At(kTrackIn))->Clone(hn);
    H->SetTitle("MC r-#phi/z/angular residuals @ TRD entry");
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
    trNbins[kPt] =fNpt-1;trMax[kPt]   = trNbins[kPt]-.5;
    //trNbins[kSpeciesChgRC] = Int_t(kNcharge)*AliPID::kSPECIES+1;trMin[kSpeciesChgRC] = -AliPID::kSPECIES-0.5; trMax[kSpeciesChgRC] = AliPID::kSPECIES+0.5;
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
    Float_t xyz[3] = {(Float_t)x, (Float_t)y, (Float_t)z}; points[ip].SetXYZ(xyz);
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
  Float_t dy/*, dz*/, s, m;
  for(Int_t ip(0); ip<np; ip++){
    if(points[ip].GetClusterType()) continue;
    Float_t x0(points[ip].GetX())
          ,y0(param[1] + param[3] * (x0 - param[0]))
          /*,z0(param[2] + param[4] * (x0 - param[0]))*/;
    dy=points[ip].GetY() - y0; h.Fill(dy);
    //dz=points[ip].GetZ() - z0;
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
void AliTRDresolution::SetProcesses(Bool_t det, Bool_t cl, Bool_t trklt, Bool_t trkin)
{
// steer processes
  if(det) SETBIT(fSteer, kDetector); else CLRBIT(fSteer, kDetector);
  if(cl)  SETBIT(fSteer, kCluster); else CLRBIT(fSteer, kCluster);
  if(trklt) SETBIT(fSteer, kTracklet); else CLRBIT(fSteer, kTracklet);
  if(trkin) SETBIT(fSteer, kTrackIn); else CLRBIT(fSteer, kTrackIn);
}


//________________________________________________________
void AliTRDresolution::SetDump3D(Bool_t det, Bool_t cl, Bool_t trklt, Bool_t trkin)
{
// steer processes
  if(det) SETBIT(fSteer, 4+kDetector); else CLRBIT(fSteer, 4+kDetector);
  if(cl)  SETBIT(fSteer, 4+kCluster); else CLRBIT(fSteer, 4+kCluster);
  if(trklt) SETBIT(fSteer, 4+kTracklet); else CLRBIT(fSteer, 4+kTracklet);
  if(trkin) SETBIT(fSteer, 4+kTrackIn); else CLRBIT(fSteer, 4+kTrackIn);
}

