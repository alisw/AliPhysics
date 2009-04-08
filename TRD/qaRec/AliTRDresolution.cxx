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
//   res->Load("TRD.TaskResolution.root");
//   if(!res->PostProcess()) return;
//   res->GetRefFigure(0);
// }  
//
//  Authors:                                                              //
//    Alexandru Bercuci <A.Bercuci@gsi.de>                                //
//    Markus Fasel <M.Fasel@gsi.de>                                       //
//                                                                        //
////////////////////////////////////////////////////////////////////////////

#include <cstring>

#include <TROOT.h>
#include <TSystem.h>
#include <TPDGCode.h>
#include <TObjArray.h>
#include <TH3.h>
#include <TH2.h>
#include <TH1.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TBox.h>
#include <TProfile.h>
#include <TGraphErrors.h>
#include <TMath.h>
#include "TTreeStream.h"
#include "TGeoManager.h"

#include "AliAnalysisManager.h"
#include "AliTrackReference.h"
#include "AliTrackPointArray.h"
#include "AliCDBManager.h"

#include "AliTRDcalibDB.h"
#include "AliTRDCommonParam.h"
#include "AliTRDSimParam.h"
#include "AliTRDgeometry.h"
#include "AliTRDpadPlane.h"
#include "AliTRDcluster.h"
#include "AliTRDseedV1.h"
#include "AliTRDtrackV1.h"
#include "AliTRDtrackerV1.h"
#include "AliTRDReconstructor.h"
#include "AliTRDrecoParam.h"

#include "AliTRDtrackInfo/AliTRDclusterInfo.h"
#include "AliTRDtrackInfo/AliTRDtrackInfo.h"
#include "AliTRDresolution.h"

ClassImp(AliTRDresolution)

//________________________________________________________
AliTRDresolution::AliTRDresolution()
  :AliTRDrecoTask("Resolution", "Spatial and Momentum Resolution")
  ,fStatus(0)
  ,fReconstructor(0x0)
  ,fGeo(0x0)
  ,fGraphS(0x0)
  ,fGraphM(0x0)
  ,fCl(0x0)
  ,fTrklt(0x0)
  ,fMCcl(0x0)
  ,fMCtrklt(0x0)
{
  fReconstructor = new AliTRDReconstructor();
  fReconstructor->SetRecoParam(AliTRDrecoParam::GetLowFluxParam());
  fGeo = new AliTRDgeometry();

  InitFunctorList();

  DefineOutput(1, TObjArray::Class()); // cluster2track
  DefineOutput(2, TObjArray::Class()); // tracklet2track
  DefineOutput(3, TObjArray::Class()); // cluster2mc
  DefineOutput(4, TObjArray::Class()); // tracklet2mc
}

//________________________________________________________
AliTRDresolution::~AliTRDresolution()
{
  if(fGraphS){fGraphS->Delete(); delete fGraphS;}
  if(fGraphM){fGraphM->Delete(); delete fGraphM;}
  delete fGeo;
  delete fReconstructor;
  if(gGeoManager) delete gGeoManager;
  if(fCl){fCl->Delete(); delete fCl;}
  if(fTrklt){fTrklt->Delete(); delete fTrklt;}
  if(fMCcl){fMCcl->Delete(); delete fMCcl;}
  if(fMCtrklt){fMCtrklt->Delete(); delete fMCtrklt;}
}


//________________________________________________________
void AliTRDresolution::CreateOutputObjects()
{
  // spatial resolution
  OpenFile(0, "RECREATE");

  fContainer = Histos();

  fCl = new TObjArray();
  fCl->SetOwner(kTRUE);
  fTrklt = new TObjArray();
  fTrklt->SetOwner(kTRUE);
  fMCcl = new TObjArray();
  fMCcl->SetOwner(kTRUE);
  fMCtrklt = new TObjArray();
  fMCtrklt->SetOwner(kTRUE);
}

//________________________________________________________
void AliTRDresolution::Exec(Option_t *opt)
{
  fCl->Delete();
  fTrklt->Delete();
  fMCcl->Delete();
  fMCtrklt->Delete();

  AliTRDrecoTask::Exec(opt);

  PostData(1, fCl);
  PostData(2, fTrklt);
  PostData(3, fMCcl);
  PostData(4, fMCtrklt);
}

//________________________________________________________
TH1* AliTRDresolution::PlotCluster(const AliTRDtrackV1 *track)
{
  if(track) fTrack = track;
  if(!fTrack){
    AliWarning("No Track defined.");
    return 0x0;
  }
  TH1 *h = 0x0;
  if(!(h = ((TH2I*)fContainer->At(kCluster)))){
    AliWarning("No output histogram defined.");
    return 0x0;
  }

  Double_t cov[3];
  Float_t x0, y0, z0, dy, dydx, dzdx;
  AliTRDseedV1 *fTracklet = 0x0;  
  for(Int_t ily=0; ily<AliTRDgeometry::kNlayer; ily++){
    if(!(fTracklet = fTrack->GetTracklet(ily))) continue;
    if(!fTracklet->IsOK()) continue;
    x0 = fTracklet->GetX0();

    // retrive the track angle with the chamber
    y0   = fTracklet->GetYref(0);
    z0   = fTracklet->GetZref(0);
    dydx = fTracklet->GetYref(1);
    dzdx = fTracklet->GetZref(1);
    fTracklet->GetCovRef(cov);
    Float_t tilt = fTracklet->GetTilt();
    AliTRDcluster *c = 0x0;
    fTracklet->ResetClusterIter(kFALSE);
    while((c = fTracklet->PrevCluster())){
      Float_t xc = c->GetX();
      Float_t yc = c->GetY();
      Float_t zc = c->GetZ();
      Float_t dx = x0 - xc; 
      Float_t yt = y0 - dx*dydx;
      Float_t zt = z0 - dx*dzdx; 
      yc -= tilt*(zc-zt); // tilt correction
      dy = yt - yc;

      h->Fill(dydx, dy/TMath::Sqrt(cov[0] + c->GetSigmaY2()));
  
      if(fDebugLevel>=1){
        // Get z-position with respect to anode wire
        //AliTRDSimParam    *simParam    = AliTRDSimParam::Instance();
        Int_t istk = fGeo->GetStack(c->GetDetector());
        AliTRDpadPlane *pp = fGeo->GetPadPlane(ily, istk);
        Float_t row0 = pp->GetRow0();
        Float_t d  =  row0 - zt + pp->GetAnodeWireOffset();
        d -= ((Int_t)(2 * d)) / 2.0;
        if (d > 0.25) d  = 0.5 - d;

/*        AliTRDclusterInfo *clInfo = new AliTRDclusterInfo;
        fCl->Add(clInfo);
        clInfo->SetCluster(c);
        clInfo->SetGlobalPosition(yt, zt, dydx, dzdx);
        clInfo->SetResolution(dy);
        clInfo->SetAnisochronity(d);
        clInfo->SetDriftLength(dx);
        (*fDebugStream) << "ClusterResiduals"
          <<"clInfo.=" << clInfo
          << "\n";*/
      }
    }
  }
  return h;
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
  if(track) fTrack = track;
  if(!fTrack){
    AliWarning("No Track defined.");
    return 0x0;
  }
  TH1 *h = 0x0;
  if(!(h = ((TH2I*)fContainer->At(kTracklet)))){
    AliWarning("No output histogram defined.");
    return 0x0;
  }

  Double_t cov[3], covR[3], sqr[3], inv[3];
  Float_t x, dx, dy, dz;
  AliTRDseedV1 *fTracklet = 0x0;  
  for(Int_t il=AliTRDgeometry::kNlayer; il--;){
    if(!(fTracklet = fTrack->GetTracklet(il))) continue;
    if(!fTracklet->IsOK()) continue;
    x    = fTracklet->GetX();
    dx   = fTracklet->GetX0() - x;
    // compute dy^2 and dz^2
    dy   = fTracklet->GetYref(0)-dx*fTracklet->GetYref(1) - fTracklet->GetY();
    dz   = fTracklet->GetZref(0)-dx*fTracklet->GetZref(1) - fTracklet->GetZ();
    // compute covariance matrix
    fTracklet->GetCovAt(x, cov);
    fTracklet->GetCovRef(covR);
    cov[0] += covR[0]; cov[1] += covR[1]; cov[2] += covR[2]; 
    // compute square root matrix
    if(AliTRDseedV1::GetCovInv(cov, inv)==0.) continue;
    if(AliTRDseedV1::GetCovSqrt(inv, sqr)<0.) continue;
    
    Double_t y = sqr[0]*dy+sqr[1]*dz;
    Double_t z = sqr[1]*dy+sqr[2]*dz;
    ((TH3*)h)->Fill(y, z, fTracklet->GetYref(1));
  }
  return h;
}

//________________________________________________________
TH1* AliTRDresolution::PlotTrackletPhi(const AliTRDtrackV1 *track)
{
  if(track) fTrack = track;
  if(!fTrack){
    AliWarning("No Track defined.");
    return 0x0;
  }
  TH1 *h = 0x0;
  if(!(h = ((TH2I*)fContainer->At(kTrackletPhi)))){
    AliWarning("No output histogram defined.");
    return 0x0;
  }

  AliTRDseedV1 *tracklet = 0x0;
  for(Int_t il=0; il<AliTRDgeometry::kNlayer; il++){
    if(!(tracklet = fTrack->GetTracklet(il))) continue;
    if(!tracklet->IsOK()) continue;
    h->Fill(tracklet->GetYref(1), TMath::ATan(tracklet->GetYref(1))-TMath::ATan(tracklet->GetYfit(1)));
  }
  return h;
}

//________________________________________________________
TH1* AliTRDresolution::PlotTrackIn(const AliTRDtrackV1 *track)
{
  if(track) fTrack = track;
  if(!fTrack){
    AliWarning("No Track defined.");
    return 0x0;
  }
  AliExternalTrackParam *tin = 0x0;
  if(!(tin = fTrack->GetTrackLow())){
    AliWarning("Track did not entered TRD fiducial volume.");
    return 0x0;
  }
  TH1 *h = 0x0;
  
  Double_t x = tin->GetX();
  AliTRDseedV1 *tracklet = 0x0;  
  for(Int_t ily=0; ily<AliTRDgeometry::kNlayer; ily++){
    if(!(tracklet = fTrack->GetTracklet(ily))) continue;
    break;
  }
  if(!tracklet || TMath::Abs(x-tracklet->GetX())>1.e-3){
    AliWarning("Tracklet did not match TRD entrance.");
    return 0x0;
  }
  const Int_t kNPAR(5);
  Double_t parR[kNPAR]; memcpy(parR, tin->GetParameter(), kNPAR*sizeof(Double_t));
  Double_t covR[3*kNPAR]; memcpy(covR, tin->GetCovariance(), 3*kNPAR*sizeof(Double_t));
  Double_t cov[3]; tracklet->GetCovAt(x, cov);

  // define sum covariances
  TMatrixDSym COV(kNPAR);
  Double_t *pc = &covR[0];
  for(Int_t ir=0, idx=0; ir<kNPAR; ir++){
    for(Int_t ic = 0; ic<=ir; ic++,pc++,idx++){ 
      COV(ir,ic) = (*pc);
      COV(ic,ir) = (*pc);
      //printf("c%2d[%f] ", idx, *pc);
    }
    //printf("\n");
  }
  //COV.Print();

  COV(0,0) += cov[0]; 
  COV(1,0) += cov[1]; COV(0, 1) += cov[1];
  COV(1,1) += cov[2];
  Double_t dy = parR[0] - tracklet->GetY(); 
  Double_t dz = parR[1] - tracklet->GetZ(); 

  if(fDebugLevel>=1){
    Double_t dydx = 1./TMath::Sqrt(1.-parR[2]*parR[2]),
             pt   = 1./parR[4];
    Double_t c0 = covR[0]+cov[0],
             c1 = covR[1]+cov[1],
             c2 = covR[2]+cov[2];
    (*fDebugStream) << "trackIn"
      << "x="       << x
      << "dy="      << dy
      << "dz="      << dz
      << "dydx="    << dydx
      << "dzdx="    << parR[3]
      << "pt="      << pt
      << "s2y="     << c0
      << "cyz="     << c1
      << "s2z="     << c2
      << "\n";
  }


  if(!HasMCdata()) return h;
  UChar_t s;
  Float_t dx, pt0, x0=tracklet->GetX0(), y0, z0, dydx0, dzdx0;
  if(!fMC->GetDirections(x0, y0, z0, dydx0, dzdx0, pt0, s)) return h;
  // translate to reference radial position
  dx = x0 - x; y0 -= dx*dydx0; z0 -= dx*dzdx0;
  // re-set covariance matrix
  pc = &covR[0];
  for(Int_t ir = 0; ir<kNPAR; ir++){
    for(Int_t ic = 0; ic<=ir; ic++, pc++){ 
      COV(ir,ic) = (*pc); COV(ic,ir) = (*pc);
    }
  }
//   TMatrixDSymEigen eigen(COV);
//   TVectorD evals = eigen.GetEigenValues();
//   TMatrixDSym evalsm(kNPAR);
//   for(Int_t ir=0; ir<kNPAR; ir++) for(Int_t ic=0; ic<kNPAR; ic++) evalsm(ir,ic) = (ir==ic ? evals(ir): 0.);
//   TMatrixD evecs = eigen.GetEigenVectors();
//   TMatrixD sqrcov(evecs, TMatrixD::kMult, TMatrixD(evalsm, TMatrixD::kMult, evecs.T()));

  // calculate delta
  parR[0]-=y0;parR[1]-=z0;
  parR[2]-=1./TMath::Sqrt(1.+dydx0*dydx0); parR[3]-=dzdx0;
  parR[4]-=1./pt0;
  // 
  if(fDebugLevel>=1){
    (*fDebugStream) << "trackInMC"
      << "x="     << x
      << "dp0="   << parR[0]
      << "dp1="   << parR[1]
      << "dp2="   << parR[2]
      << "dp3="   << parR[3]
      << "dp4="   << parR[4]
      << "cov="   << &COV
      << "\n";
  }
  return h;
}

//________________________________________________________
TH1* AliTRDresolution::PlotMC(const AliTRDtrackV1 *track)
{
  if(!HasMCdata()){ 
    AliWarning("No MC defined. Results will not be available.");
    return 0x0;
  }
  if(track) fTrack = track;
  if(!fTrack){
    AliWarning("No Track defined.");
    return 0x0;
  }
  TH1 *h = 0x0;
  UChar_t s;
  Int_t pdg = fMC->GetPDG(), det=-1;
  Int_t label = fMC->GetLabel();
  Double_t xAnode, x, y, z, pt, dydx, dzdx;
  Float_t pt0, x0, y0, z0, dx, dy, dz, dydx0, dzdx0;
  Double_t covR[3]/*, cov[3]*/;

  if(fDebugLevel>=1){
    Double_t DX[12], DY[12], DZ[12], DPt[12], COV[12][15];
    fMC->PropagateKalman(DX, DY, DZ, DPt, COV);
    (*fDebugStream) << "MCkalman"
      << "pdg="  << pdg
      << "dx0="  << DX[0]
      << "dx1="  << DX[1]
      << "dx2="  << DX[2]
      << "dy0="  << DY[0]
      << "dy1="  << DY[1]
      << "dy2="  << DY[2]
      << "dz0="  << DZ[0]
      << "dz1="  << DZ[1]
      << "dz2="  << DZ[2]
      << "dpt0=" << DPt[0]
      << "dpt1=" << DPt[1]
      << "dpt2=" << DPt[2]
      << "\n";
  }

  AliTRDseedV1 *fTracklet = 0x0;  
  for(Int_t ily=0; ily<AliTRDgeometry::kNlayer; ily++){
    if(!(fTracklet = fTrack->GetTracklet(ily)))/* ||
       !fTracklet->IsOK())*/ continue;

    det = fTracklet->GetDetector();
    x0  = fTracklet->GetX0();
    //radial shift with respect to the MC reference (radial position of the pad plane)
    x= fTracklet->GetX();
    if(!fMC->GetDirections(x0, y0, z0, dydx0, dzdx0, pt0, s)) continue;
    xAnode  = fTracklet->GetX0();

    // MC track position at reference radial position
    dx  = x0 - x;
    if(fDebugLevel>=1){
      (*fDebugStream) << "MC"
        << "det="     << det
        << "pdg="     << pdg
        << "pt="      << pt0
        << "x="       << x0
        << "y="       << y0
        << "z="       << z0
        << "dydx="    << dydx0
        << "dzdx="    << dzdx0
        << "\n";
    }
    Float_t yt = y0 - dx*dydx0;
    Float_t zt = z0 - dx*dzdx0;
    //p = pt0*TMath::Sqrt(1.+dzdx0*dzdx0); // pt -> p

    // Kalman position at reference radial position
    dx = xAnode - x;
    y  = fTracklet->GetYref(0) - dx*fTracklet->GetYref(1);
    dy = yt - y;
    z  = fTracklet->GetZref(0) - dx*fTracklet->GetZref(1);
    dz = zt - z;
    dzdx = fTracklet->GetTgl();
    pt = fTracklet->GetPt();
    fTracklet->GetCovRef(covR);

    ((TH2I*)fContainer->At(kMCtrackY))->Fill(dydx0, dy);
    ((TH2I*)fContainer->At(kMCtrackYPull))->Fill(dydx0, dy/TMath::Sqrt(covR[0]));
    if(ily==0){
      ((TH2I*)fContainer->At(kMCtrackZIn))->Fill(dzdx0, dz);
      ((TH2I*)fContainer->At(kMCtrackZInPull))->Fill(dzdx0, dz/TMath::Sqrt(covR[2]));
    } else if(ily==AliTRDgeometry::kNlayer-1) {
      ((TH2I*)fContainer->At(kMCtrackZOut))->Fill(dzdx0, dz);
      ((TH2I*)fContainer->At(kMCtrackZOutPull))->Fill(dzdx0, dz/TMath::Sqrt(covR[2]));
    }
    if(pdg!=kElectron && pdg!=kPositron) ((TH2I*)fContainer->At(kMCtrackPt))->Fill(1./pt0, pt-pt0);
    // Fill Debug stream for Kalman track
    if(fDebugLevel>=1){
      dydx = fTracklet->GetYref(1);
      (*fDebugStream) << "MCtrack"
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

    // recalculate tracklet based on the MC info
    AliTRDseedV1 tt(*fTracklet);
    tt.SetZref(0, z0 - (x0-xAnode)*dzdx0);
    tt.SetZref(1, dzdx0); 
    tt.Fit(kTRUE, kTRUE);
    x= tt.GetX();y= tt.GetY();z= tt.GetZ();
    dydx = tt.GetYfit(1);
    dx = x0 - x;
    yt = y0 - dx*dydx0;
    zt = z0 - dx*dzdx0;
    Bool_t rc = tt.IsRowCross(); 
    
    // add tracklet residuals for y and dydx
    if(!rc){
      dy    = yt-y;

      Float_t dphi  = (dydx - dydx0);
      dphi /= 1.- dydx*dydx0;

      ((TH2I*)fContainer->At(kMCtrackletY))->Fill(dydx0, dy);
      if(tt.GetS2Y()>0.) ((TH2I*)fContainer->At(kMCtrackletYPull))->Fill(dydx0, dy/TMath::Sqrt(tt.GetS2Y()));
      ((TH2I*)fContainer->At(kMCtrackletPhi))->Fill(dydx0, dphi*TMath::RadToDeg());
    } else {
      // add tracklet residuals for z
      dz = zt-z;
      ((TH2I*)fContainer->At(kMCtrackletZ))->Fill(dzdx0, dz);
      if(tt.GetS2Z()>0.) ((TH2I*)fContainer->At(kMCtrackletZPull))->Fill(dzdx0, dz/TMath::Sqrt(tt.GetS2Z()));
    }
  
    // Fill Debug stream for tracklet
    if(fDebugLevel>=1){
      Float_t s2y = tt.GetS2Y();
      Float_t s2z = tt.GetS2Z();
      (*fDebugStream) << "MCtracklet"
        << "rc="    << rc
        << "x="     << x
        << "y="     << y
        << "z="     << z
        << "dydx="  << dydx
        << "s2y="   << s2y
        << "s2z="   << s2z
        << "\n";
    }

    Int_t istk = AliTRDgeometry::GetStack(det); 
    AliTRDpadPlane *pp = fGeo->GetPadPlane(ily, istk);
    Float_t zr0 = pp->GetRow0() + pp->GetAnodeWireOffset();
    Float_t tilt = fTracklet->GetTilt();
    Double_t exb = AliTRDCommonParam::Instance()->GetOmegaTau(1.5);

    AliTRDcluster *c = 0x0;
    fTracklet->ResetClusterIter(kFALSE);
    while((c = fTracklet->PrevCluster())){
      Float_t  q = TMath::Abs(c->GetQ());
      x = c->GetX(); y = c->GetY();
//       Int_t col = c->GetPadCol();
//       Int_t row = c->GetPadRow();
//       Double_t pw = pp->GetColSize(col);
//       Double_t y0 = (pp->GetColPos(col) + 0.5) * pw;
//       Double_t s2 = AliTRDcalibDB::Instance()->GetPRFWidth(det, col, row); s2 *= s2; s2 -= - 1.5e-1;
//       y = c->GetYloc(y0, s2, pw); y-=(xAnode-x)*exb;

      z = c->GetZ();
      dx = x0 - x; 
      yt = y0 - dx*dydx0;
      zt = z0 - dx*dzdx0;
      dy = yt - (y - tilt*(z-zt));

      // Fill Histograms
      if(q>20. && q<250.) ((TH2I*)fContainer->At(kMCcluster))->Fill(dydx0, dy);

      // Fill calibration container
      Float_t d = zr0 - zt;
      d -= ((Int_t)(2 * d)) / 2.0;
      if (d > 0.25) d  = 0.5 - d;
      AliTRDclusterInfo *clInfo = new AliTRDclusterInfo;
      fMCcl->Add(clInfo);
      clInfo->SetCluster(c);
      clInfo->SetMC(pdg, label);
      clInfo->SetGlobalPosition(yt, zt, dydx0, dzdx0);
      clInfo->SetResolution(dy);
      clInfo->SetAnisochronity(d);
      clInfo->SetDriftLength(((c->GetPadTime()+0.5)*.1)*1.5);
      //dx-.5*AliTRDgeometry::CamHght());
      clInfo->SetTilt(tilt);

      // Fill Debug Tree
      if(fDebugLevel>=2){
        //clInfo->Print();
        (*fDebugStream) << "MCcluster"
          <<"clInfo.=" << clInfo
          << "\n";
      }
    }
  }
  return h;
}


//________________________________________________________
Bool_t AliTRDresolution::GetRefFigure(Int_t ifig)
{
  Float_t y[2] = {0., 0.};
  TBox *b = 0x0;
  TAxis *ax = 0x0;
  TGraphErrors *g = 0x0;
  switch(ifig){
  case kCluster:
    if(!(g = (TGraphErrors*)fGraphS->At(ifig))) break;
    g->Draw("apl");
    ax = g->GetHistogram()->GetYaxis();
    y[0] = -0.5; y[1] = 2.5;
    ax->SetRangeUser(y[0], y[1]);
    ax->SetTitle("Cluster-Track Pulls #sigma/#mu [mm]");
    ax = g->GetHistogram()->GetXaxis();
    ax->SetTitle("tg(#phi)");
    if(!(g = (TGraphErrors*)fGraphM->At(ifig))) break;
    g->Draw("pl");
    b = new TBox(-.15, y[0], .15, y[1]);
    b->SetFillStyle(3002);b->SetFillColor(kGreen);
    b->SetLineColor(0); b->Draw();
    return kTRUE;
  case kTracklet:
    if(!(g = (TGraphErrors*)fGraphS->At(ifig))) break;
    g->Draw("apl");
    ax = g->GetHistogram()->GetYaxis();
    ax->SetRangeUser(-.5, 3.);
    ax->SetTitle("Tracklet-Track YZ-Pulls #sigma/#mu [mm]");
    ax = g->GetHistogram()->GetXaxis();
    ax->SetTitle("tg(#phi)");
    if(!(g = (TGraphErrors*)fGraphM->At(ifig))) break;
    g->Draw("pl");
    b = new TBox(-.15, y[0], .15, y[1]);
    b->SetFillStyle(3002);b->SetFillColor(kGreen);
    b->SetLineColor(0); b->Draw();
    return kTRUE;
  case kTrackletPhi:
    if(!(g = (TGraphErrors*)fGraphS->At(ifig))) break;
    g->Draw("apl");
    ax = g->GetHistogram()->GetYaxis();
    ax->SetRangeUser(-.5, 2.);
    ax->SetTitle("Tracklet-Track Phi-Residuals #sigma/#mu [rad]");
    ax = g->GetHistogram()->GetXaxis();
    ax->SetTitle("tg(#phi)");
    if(!(g = (TGraphErrors*)fGraphM->At(ifig))) break;
    g->Draw("pl");
    b = new TBox(-.15, y[0], .15, y[1]);
    b->SetFillStyle(3002);b->SetFillColor(kGreen);
    b->SetLineColor(0); b->Draw();
    return kTRUE;
  case kMCcluster:
    if(!(g = (TGraphErrors*)fGraphS->At(ifig))) break;
    ax = g->GetHistogram()->GetYaxis();
    y[0] = -50.; y[1] = 600.;
    ax->SetRangeUser(y[0], y[1]);
    ax->SetTitle("Y_{cluster} #sigma/#mu [#mum]");
    ax = g->GetHistogram()->GetXaxis();
    ax->SetRangeUser(-.3, .3);
    ax->SetTitle("tg(#phi)");
    g->Draw("apl");
    if(!(g = (TGraphErrors*)fGraphM->At(ifig))) break;
    g->Draw("pl");
    b = new TBox(-.15, y[0], .15, y[1]);
    b->SetFillStyle(3002);b->SetFillColor(kBlue);
    b->SetLineColor(0); b->Draw();
    return kTRUE;
  case kMCtrackletY:
    if(!(g = (TGraphErrors*)fGraphS->At(ifig))) break;
    ax = g->GetHistogram()->GetYaxis();
    y[0] = -50.; y[1] = 250.;
    ax->SetRangeUser(y[0], y[1]);
    ax->SetTitle("Y_{tracklet} #sigma/#mu [#mum]");
    ax = g->GetHistogram()->GetXaxis();
    ax->SetRangeUser(-.2, .2);
    ax->SetTitle("tg(#phi)");
    g->Draw("apl");
    if(!(g = (TGraphErrors*)fGraphM->At(ifig))) break;
    g->Draw("pl");
    b = new TBox(-.15, y[0], .15, y[1]);
    b->SetFillStyle(3002);b->SetFillColor(kBlue);
    b->SetLineColor(0); b->Draw();
    return kTRUE;
  case kMCtrackletZ:
    if(!(g = (TGraphErrors*)fGraphS->At(ifig))) break;
    ax = g->GetHistogram()->GetYaxis();
    ax->SetRangeUser(-50., 700.);
    ax->SetTitle("Z_{tracklet} #sigma/#mu [#mum]");
    ax = g->GetHistogram()->GetXaxis();
    ax->SetTitle("tg(#theta)");
    g->Draw("apl");
    if(!(g = (TGraphErrors*)fGraphM->At(ifig))) break;
    g->Draw("pl");
    return kTRUE;
  case kMCtrackletPhi:
    if(!(g = (TGraphErrors*)fGraphS->At(ifig))) break;
    ax = g->GetHistogram()->GetYaxis();
    y[0] = -.05; y[1] = .2;
    ax->SetRangeUser(y[0], y[1]);
    ax->SetTitle("#Phi_{tracklet} #sigma/#mu [deg]");
    ax = g->GetHistogram()->GetXaxis();
    ax->SetTitle("tg(#phi)");
    g->Draw("apl");
    if(!(g = (TGraphErrors*)fGraphM->At(ifig))) break;
    g->Draw("pl");
    return kTRUE;
  case kMCtrackY:
    if(!(g = (TGraphErrors*)fGraphS->At(ifig))) break;
    ax = g->GetHistogram()->GetYaxis();
    y[0] = -50.; y[1] = 200.;
    ax->SetRangeUser(y[0], y[1]);
    ax->SetTitle("Y_{track} #sigma/#mu [#mum]");
    ax = g->GetHistogram()->GetXaxis();
    ax->SetRangeUser(-.2, .2);
    ax->SetTitle("tg(#phi)");
    g->Draw("apl");
    if(!(g = (TGraphErrors*)fGraphM->At(ifig))) break;
    g->Draw("pl");
    b = new TBox(-.15, y[0], .15, y[1]);
    b->SetFillStyle(3002);b->SetFillColor(kBlue);
    b->SetLineColor(0); b->Draw();
    return kTRUE;
  case kMCtrackZIn:
  case kMCtrackZOut:
    if(!(g = (TGraphErrors*)fGraphS->At(ifig))) break;
    ax = g->GetHistogram()->GetYaxis();
    ax->SetRangeUser(-500., 2000.);
    ax->SetTitle(Form("Z_{track}^{%s} #sigma/#mu [#mum]", ifig==kMCtrackZIn ? "in" : "out"));
    ax = g->GetHistogram()->GetXaxis();
    ax->SetTitle("tg(#theta)");
    g->Draw("apl");
    if(!(g = (TGraphErrors*)fGraphM->At(ifig))) break;
    g->Draw("pl");
    return kTRUE;
  case kMCtrackPt:
    if(!(g = (TGraphErrors*)fGraphS->At(ifig))) break;
    ax = g->GetHistogram()->GetYaxis();
    ax->SetRangeUser(-.5, 2.);
    ax->SetTitle("#epsilon_{P_{t}}^{track} / #mu [%]");
    ax = g->GetHistogram()->GetXaxis();
    ax->SetTitle("1/p_{t}");
    g->Draw("apl");
    if(!(g = (TGraphErrors*)fGraphM->At(ifig))) break;
    g->Draw("pl");
    return kTRUE;
  case kMCtrackletYPull:
  case kMCtrackletZPull:
  case kMCtrackYPull:
  case kMCtrackZInPull:
  case kMCtrackZOutPull:
    if(!(g = (TGraphErrors*)fGraphS->At(ifig))) break;
    ax = g->GetHistogram()->GetYaxis();
    ax->SetRangeUser(-.5, 2.);
    ax->SetTitle("MC Pulls");
    g->Draw("apl");
    if(!(g = (TGraphErrors*)fGraphM->At(ifig))) break;
    g->Draw("pl");
    return kTRUE;
  }
  AliInfo(Form("Reference plot [%d] missing result", ifig));
  return kFALSE;
}


//________________________________________________________
Bool_t AliTRDresolution::PostProcess()
{
  //fContainer = dynamic_cast<TObjArray*>(GetOutputData(0));
  if (!fContainer) {
    Printf("ERROR: list not available");
    return kFALSE;
  }
  fNRefFigures = fContainer->GetEntriesFast();
  TGraphErrors *gm= 0x0, *gs= 0x0;
  if(!fGraphS){ 
    fGraphS = new TObjArray(fNRefFigures);
    fGraphS->SetOwner();
    for(Int_t ig=0; ig<fNRefFigures; ig++){
      gs = new TGraphErrors();
      gs->SetLineColor(kRed);
      gs->SetMarkerStyle(23);
      gs->SetMarkerColor(kRed);
      gs->SetNameTitle(Form("s_%d", ig), "");
      fGraphS->AddAt(gs, ig);
    }
  }
  if(!fGraphM){ 
    fGraphM = new TObjArray(fNRefFigures);
    fGraphM->SetOwner();
    for(Int_t ig=0; ig<fNRefFigures; ig++){
      gm = new TGraphErrors();
      gm->SetLineColor(kBlue);
      gm->SetMarkerStyle(7);
      gm->SetMarkerColor(kBlue);
      gm->SetNameTitle(Form("m_%d", ig), "");
      fGraphM->AddAt(gm, ig);
    }
  }

  // DEFINE MODELS
  // simple gauss
  TF1 f("f1", "gaus", -.5, .5);  
  // gauss on a constant background
  TF1 fb("fb", "[0]*exp(-0.5*((x-[1])/[2])**2)+[3]", -.5, .5);
  // gauss on a gauss background
  TF1 fc("fc", "[0]*exp(-0.5*((x-[1])/[2])**2)+[3]*exp(-0.5*((x-[4])/[5])**2)", -.5, .5);


  //PROCESS EXPERIMENTAL DISTRIBUTIONS

  // Clusters residuals
  Process(kCluster, &f);

  // Tracklet y residuals
  Process3D(kTracklet, &f);

  // Tracklet phi residuals
  Process(kTrackletPhi, &f);

  if(!HasMCdata()) return kTRUE;


  //PROCESS MC RESIDUAL DISTRIBUTIONS

  // cluster y resolution
  Process(kMCcluster, &f, 1.e4);

  // tracklet resolution
  Process(kMCtrackletY, &f, 1.e4); // y
  Process(kMCtrackletZ, &f, 1.e4); // z
  Process(kMCtrackletPhi, &f); // phi

  // tracklet pulls
  Process(kMCtrackletYPull, &f); // y
  Process(kMCtrackletZPull, &f); // z

  // track resolution
  Process(kMCtrackY, &f, 1.e4);    // y
  Process(kMCtrackZIn, &f, 1.e4);  // z towards TPC
  Process(kMCtrackZOut, &f, 1.e4); // z towards TOF

  // track pulls
  Process(kMCtrackYPull, &f);    // y
  Process(kMCtrackZInPull, &f);  // z towards TPC
  Process(kMCtrackZOutPull, &f); // z towards TOF

  // track Pt resolution
  TH2I *h2 = (TH2I*)fContainer->At(kMCtrackPt);
  TAxis *ax = h2->GetXaxis();
  gm = (TGraphErrors*)fGraphM->At(kMCtrackPt);
  gs = (TGraphErrors*)fGraphS->At(kMCtrackPt);
  TF1 fg("fg", "gaus", -1.5, 1.5);
  TF1 fl("fl", "landau", -4., 15.);
  TF1 fgl("fgl", "gaus(0)+landau(3)", -5., 20.);
  for(Int_t ip=1; ip<=ax->GetNbins(); ip++){
    TH1D *h = h2->ProjectionY("ppt", ip, ip);
    if(h->GetEntries()<70) continue;

    h->Fit(&fg, "QN", "", -1.5, 1.5);
    fgl.SetParameter(0, fg.GetParameter(0));
    fgl.SetParameter(1, fg.GetParameter(1));
    fgl.SetParameter(2, fg.GetParameter(2));
    h->Fit(&fl, "QN", "", -4., 15.);
    fgl.SetParameter(3, fl.GetParameter(0));
    fgl.SetParameter(4, fl.GetParameter(1));
    fgl.SetParameter(5, fl.GetParameter(2));

    h->Fit(&fgl, "NQ", "", -5., 20.);

    Float_t invpt = ax->GetBinCenter(ip);
    Int_t jp = gm->GetN();
    gm->SetPoint(jp, invpt, fgl.GetParameter(1));
    gm->SetPointError(jp, 0., fgl.GetParError(1));
    gs->SetPoint(jp, invpt, fgl.GetParameter(2)*invpt);
    gs->SetPointError(jp, 0., fgl.GetParError(2));
    // fgl.GetParameter(4) // Landau MPV
    // fgl.GetParameter(5) // Landau Sigma
  }


  return kTRUE;
}


//________________________________________________________
void AliTRDresolution::Terminate(Option_t *)
{
  if(fDebugStream){ 
    delete fDebugStream;
    fDebugStream = 0x0;
    fDebugLevel = 0;
  }
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

  f->SetParLimits(1, -.2, .2);
  f->SetParameter(1, -0.1);

  f->SetParLimits(2, 0., 4.e-1);
  f->SetParameter(2, 2.e-2);
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
  if(fContainer) return fContainer;

  fContainer  = new TObjArray(kNhistos);
  //fContainer->SetOwner(kTRUE);

  TH1 *h = 0x0;
  // cluster to tracklet residuals [2]
  if(!(h = (TH2I*)gROOT->FindObject("hCl"))){
    h = new TH2I("hCl", "Clusters-Track Residuals", 21, -.33, .33, 100, -.5, .5);
    h->GetXaxis()->SetTitle("tg(#phi)");
    h->GetYaxis()->SetTitle("#Delta y [cm]");
    h->GetZaxis()->SetTitle("entries");
  } else h->Reset();
  fContainer->AddAt(h, kCluster);

  // tracklet to track residuals [2]
  if(!(h = (TH3I*)gROOT->FindObject("hTrklt"))){
    h = new TH3I("hTrklt", "Tracklets-Track Residuals", 100, -.5, .5, 100, -.5, .5, 21, -.33, .33);
    h->GetXaxis()->SetTitle("#Delta y [cm]");
    h->GetYaxis()->SetTitle("#Delta z [cm]");
    h->GetZaxis()->SetTitle("tg(#phi)");
  } else h->Reset();
  fContainer->AddAt(h, kTracklet);

  // tracklet to track residuals angular [2]
  if(!(h = (TH2I*)gROOT->FindObject("hTrkltPhi"))){
    h = new TH2I("hTrkltPhi", "Tracklets-Track Residuals (#Phi)", 21, -.33, .33, 100, -.5, .5);
    h->GetXaxis()->SetTitle("tg(#phi)");
    h->GetYaxis()->SetTitle("#Delta phi [#circ]");
    h->GetZaxis()->SetTitle("entries");
  } else h->Reset();
  fContainer->AddAt(h, kTrackletPhi);


  // Resolution histos
  if(!HasMCdata()) return fContainer;

  // cluster y resolution [0]
  if(!(h = (TH2I*)gROOT->FindObject("hMCcl"))){
    h = new TH2I("hMCcl", "Cluster Resolution", 31, -.48, .48, 100, -.5, .5);
    h->GetXaxis()->SetTitle("tg(#phi)");
    h->GetYaxis()->SetTitle("#Delta y [cm]");
    h->GetZaxis()->SetTitle("entries");
  } else h->Reset();
  fContainer->AddAt(h, kMCcluster);

  // tracklet y resolution [0]
  if(!(h = (TH2I*)gROOT->FindObject("hMCtrkltY"))){
    h = new TH2I("hMCtrkltY", "Tracklet Resolution (Y)", 31, -.48, .48, 100, -.2, .2);
    h->GetXaxis()->SetTitle("tg(#phi)");
    h->GetYaxis()->SetTitle("#Delta y [cm]");
    h->GetZaxis()->SetTitle("entries");
  } else h->Reset();
  fContainer->AddAt(h, kMCtrackletY);

  // tracklet y resolution [0]
  if(!(h = (TH2I*)gROOT->FindObject("hMCtrkltYPull"))){
    h = new TH2I("hMCtrkltYPull", "Tracklet Pulls (Y)", 31, -.48, .48, 100, -3.5, 3.5);
    h->GetXaxis()->SetTitle("tg(#phi)");
    h->GetYaxis()->SetTitle("#Delta y / #sigma_{y}");
    h->GetZaxis()->SetTitle("entries");
  } else h->Reset();
  fContainer->AddAt(h, kMCtrackletYPull);

  // tracklet y resolution [0]
  if(!(h = (TH2I*)gROOT->FindObject("hMCtrkltZ"))){
    h = new TH2I("hMCtrkltZ", "Tracklet Resolution (Z)", 50, -1., 1., 100, -.5, .5);
    h->GetXaxis()->SetTitle("tg(#theta)");
    h->GetYaxis()->SetTitle("#Delta z [cm]");
    h->GetZaxis()->SetTitle("entries");
  } else h->Reset();
  fContainer->AddAt(h, kMCtrackletZ);

  // tracklet y resolution [0]
  if(!(h = (TH2I*)gROOT->FindObject("hMCtrkltZPull"))){
    h = new TH2I("hMCtrkltZPull", "Tracklet Pulls (Z)", 31, -1., 1., 100, -3.5, 3.5);
    h->GetXaxis()->SetTitle("tg(#theta)");
    h->GetYaxis()->SetTitle("#Delta z / #sigma_{z}");
    h->GetZaxis()->SetTitle("entries");
  } else h->Reset();
  fContainer->AddAt(h, kMCtrackletZPull);

  // tracklet angular resolution [1]
  if(!(h = (TH2I*)gROOT->FindObject("hMCtrkltPhi"))){
    h = new TH2I("hMCtrkltPhi", "Tracklet Resolution (#Phi)", 31, -.48, .48, 100, -10., 10.);
    h->GetXaxis()->SetTitle("tg(#phi)");
    h->GetYaxis()->SetTitle("#Delta #phi [deg]");
    h->GetZaxis()->SetTitle("entries");
  } else h->Reset();
  fContainer->AddAt(h, kMCtrackletPhi);

  // Kalman track y resolution
  if(!(h = (TH2I*)gROOT->FindObject("hMCtrkY"))){
    h = new TH2I("hMCtrkY", "Kalman Track Resolution (Y)", 31, -.48, .48, 100, -.2, .2);
    h->GetXaxis()->SetTitle("tg(#phi)");
    h->GetYaxis()->SetTitle("#Delta y [cm]");
    h->GetZaxis()->SetTitle("entries");
  } else h->Reset();
  fContainer->AddAt(h, kMCtrackY);

  // Kalman track y pulls
  if(!(h = (TH2I*)gROOT->FindObject("hMCtrkYPull"))){
    h = new TH2I("hMCtrkYPull", "Kalman Track Pulls (Y)", 31, -.48, .48, 100, -3.5, 3.5);
    h->GetXaxis()->SetTitle("tg(#phi)");
    h->GetYaxis()->SetTitle("#Delta y / #sigma_{y}");
    h->GetZaxis()->SetTitle("entries");
  } else h->Reset();
  fContainer->AddAt(h, kMCtrackYPull);

  // Kalman track Z resolution [IN]
  if(!(h = (TH2I*)gROOT->FindObject("hMCtrkZIn"))){
    h = new TH2I("hMCtrkZIn", "Kalman Track Resolution (Zin)", 20, -1., 1., 100, -1., 1.);
    h->GetXaxis()->SetTitle("tg(#theta)");
    h->GetYaxis()->SetTitle("#Delta z [cm]");
    h->GetZaxis()->SetTitle("entries");
  } else h->Reset();
  fContainer->AddAt(h, kMCtrackZIn);

  // Kalman track Z resolution [OUT]
  if(!(h = (TH2I*)gROOT->FindObject("hMCtrkZOut"))){
    h = new TH2I("hMCtrkZOut", "Kalman Track Resolution (Zout)", 20, -1., 1., 100, -1., 1.);
    h->GetXaxis()->SetTitle("tg(#theta)");
    h->GetYaxis()->SetTitle("#Delta z [cm]");
    h->GetZaxis()->SetTitle("entries");
  } else h->Reset();
  fContainer->AddAt(h, kMCtrackZOut);

  // Kalman track Z pulls
  if(!(h = (TH2I*)gROOT->FindObject("hMCtrkZInPull"))){
    h = new TH2I("hMCtrkZInPull", "Kalman Track Pulls (Zin)", 20, -1., 1., 100, -4.5, 4.5);
    h->GetXaxis()->SetTitle("tg(#theta)");
    h->GetYaxis()->SetTitle("#Delta z / #sigma_{z}");
    h->GetZaxis()->SetTitle("entries");
  } else h->Reset();
  fContainer->AddAt(h, kMCtrackZInPull);

  // Kalman track Z pulls
  if(!(h = (TH2I*)gROOT->FindObject("hMCtrkZOutPull"))){
    h = new TH2I("hMCtrkZOutPull", "Kalman Track Pulls (Zout)", 20, -1., 1., 100, -4.5, 4.5);
    h->GetXaxis()->SetTitle("tg(#theta)");
    h->GetYaxis()->SetTitle("#Delta z / #sigma_{z}");
    h->GetZaxis()->SetTitle("entries");
  } else h->Reset();
  fContainer->AddAt(h, kMCtrackZOutPull);

  // Kalman track Pt resolution
  if(!(h = (TH2I*)gROOT->FindObject("hMCtrkPt"))){
    h = new TH2I("hMCtrkPt", "Kalman Track Resolution (Pt)", 100, 0., 2., 150, -1., 10.);
    h->GetXaxis()->SetTitle("1/p_{t}");
    h->GetYaxis()->SetTitle("#Delta p_{t} [GeV/c]");
    h->GetZaxis()->SetTitle("entries");
  } else h->Reset();
  fContainer->AddAt(h, kMCtrackPt);

  return fContainer;
}


//________________________________________________________
Bool_t AliTRDresolution::Process(ETRDresolutionPlot plot, TF1 *f, Float_t k)
{
  if(!fContainer || !fGraphS || !fGraphM) return kFALSE;
  Bool_t kBUILD = kFALSE;
  if(!f){ 
    f = new TF1("f1", "gaus", -.5, .5);
    kBUILD = kTRUE;
  }

  TH2I *h2 = 0x0;
  if(!(h2 = (TH2I *)(fContainer->At(plot)))) return kFALSE;
  TGraphErrors *gm = 0x0, *gs = 0x0;
  if(!(gm=(TGraphErrors*)fGraphM->At(plot))) return kFALSE;
  if(gm->GetN()) for(Int_t ip=gm->GetN(); ip--;) gm->RemovePoint(ip);
  if(!(gs=(TGraphErrors*)fGraphS->At(plot))) return kFALSE;
  if(gs->GetN()) for(Int_t ip=gs->GetN(); ip--;) gs->RemovePoint(ip);
  Char_t pn[10]; sprintf(pn, "p%02d", plot);
  for(Int_t ibin = 1; ibin <= h2->GetNbinsX(); ibin++){
    Double_t x = h2->GetXaxis()->GetBinCenter(ibin);
    TH1D *h = h2->ProjectionY(pn, ibin, ibin);
    if(h->GetEntries()<100) continue;
    AdjustF1(h, f);

    h->Fit(f, "QN");
    
    Int_t ip = gm->GetN();
    gm->SetPoint(ip, x, k*f->GetParameter(1));
    gm->SetPointError(ip, 0., k*f->GetParError(1));
    gs->SetPoint(ip, x, k*f->GetParameter(2));
    gs->SetPointError(ip, 0., k*f->GetParError(2));
  }

  if(kBUILD) delete f;
  return kTRUE;
}

//________________________________________________________
Bool_t AliTRDresolution::Process3D(ETRDresolutionPlot plot, TF1 *f, Float_t k)
{
  if(!fContainer || !fGraphS || !fGraphM) return kFALSE;
  Bool_t kBUILD = kFALSE;
  if(!f){ 
    f = new TF1("f1", "gaus", -.5, .5);
    kBUILD = kTRUE;
  }

  TH3I *h3 = 0x0;
  if(!(h3 = (TH3I*)(fContainer->At(plot)))) return kFALSE;
  TGraphErrors *gm = 0x0, *gs = 0x0;
  if(!(gm=(TGraphErrors*)fGraphM->At(plot))) return kFALSE;
  if(gm->GetN()) for(Int_t ip=gm->GetN(); ip--;) gm->RemovePoint(ip);
  if(!(gs=(TGraphErrors*)fGraphS->At(plot))) return kFALSE;
  if(gs->GetN()) for(Int_t ip=gs->GetN(); ip--;) gs->RemovePoint(ip);
  Char_t pn[10]; sprintf(pn, "p%02d", plot);
  for(Int_t ibin = 1; ibin <= h3->GetNbinsX(); ibin++){
    Double_t x = h3->GetXaxis()->GetBinCenter(ibin);
    TH1D *h = h3->ProjectionZ(pn, ibin, ibin);
    if(h->GetEntries()<100) continue;
    AdjustF1(h, f);

    h->Fit(f, "QN");
    
    Int_t ip = gm->GetN();
    gm->SetPoint(ip, x, k*f->GetParameter(1));
    gm->SetPointError(ip, 0., k*f->GetParError(1));
    gs->SetPoint(ip, x, k*f->GetParameter(2));
    gs->SetPointError(ip, 0., k*f->GetParError(2));
  }

  if(kBUILD) delete f;
  return kTRUE;
}

//________________________________________________________
void AliTRDresolution::SetRecoParam(AliTRDrecoParam *r)
{

  fReconstructor->SetRecoParam(r);
}
