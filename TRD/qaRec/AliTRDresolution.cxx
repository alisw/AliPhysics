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
#include <TGaxis.h>
#include <TBox.h>
#include <TProfile.h>
#include <TGraphErrors.h>
#include <TMath.h>
#include <TMatrixT.h>
#include <TVectorT.h>
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

#include "info/AliTRDclusterInfo.h"
#include "info/AliTRDtrackInfo.h"
#include "AliTRDresolution.h"

ClassImp(AliTRDresolution)
UChar_t AliTRDresolution::fNElements[kNhistos] = {
  1, 1, 1,
  2, 5, 10, 2, 10
};
Char_t *AliTRDresolution::fAxTitle[32][3] = {
  // ESD
  {"tg(#phi)", "PULL: #mu_{y}^{cl}", "PULL: #sigma_{y}^{cl}"}
 ,{"tg(#phi)", "PULL: #mu_{y}^{trklt}", "PULL: #sigma_{y}^{trklt}"}
 ,{"tg(#phi)", "PULL: #mu_{#phi}^{trklt}", "PULL: #sigma_{#phi}^{trklt}"}
  // MC cluster
 ,{"tg(#phi)", "#mu_{y}^{cl} [#mum]", "#sigma_{y}^{cl} [#mum]"}
 ,{"tg(#phi)", "PULL: #mu_{y}^{cl}", "PULL: #sigma_{y}^{cl}"}
  // MC tracklet
 ,{"tg(#phi)", "#mu_{y}^{trklt} [#mum]", "#sigma_{y}^{trklt} [#mum]"}
 ,{"tg(#phi)", "PULL: #mu_{y}^{trklt}", "PULL: #sigma_{y}^{trklt}"}
 ,{"tg(#theta)", "#mu_{z}^{trklt} [#mum]", "#sigma_{z}^{trklt} [#mum]"}
 ,{"tg(#theta)", "PULL: #mu_{z}^{trklt}", "PULL: #sigma_{z}^{trklt}"}
 ,{"tg(#phi)", "#mu_{#phi}^{trklt} [mrad]", "#sigma_{#phi}^{trklt} [mrad]"}
  // MC track TPC
 ,{"tg(#phi)", "#mu_{y}^{TPC trk} [#mum]", "#sigma_{y}^{TPC trk} [#mum]"}
 ,{"tg(#phi)", "PULL: #mu_{y}^{TPC trk}", "PULL: #sigma_{y}^{TPC trk}"}
 ,{"tg(#theta)", "#mu_{z}^{TPC trk} [#mum]", "#sigma_{z}^{TPC trk} [#mum]"}
 ,{"tg(#theta)", "PULL: #mu_{z}^{TPC trk}", "PULL: #sigma_{z}^{TPC trk}"}
 ,{"tg(#phi)", "#mu_{#phi}^{TPC trk} [mrad]", "#sigma_{#phi}^{TPC trk} [mrad]"}
 ,{"tg(#phi)", "PULL: #mu_{snp}^{TPC trk}", "PULL: #sigma_{snp}^{TPC trk}"}
 ,{"tg(#theta)", "#mu_{#theta}^{TPC trk} [mrad]", "#sigma_{#theta}^{TPC trk} [mrad]"}
 ,{"tg(#theta)", "PULL: #mu_{tgl}^{TPC trk}", "PULL: #sigma_{tgl}^{TPC trk}"}
 ,{"p_{t}^{MC} [GeV/c]", "#mu^{TPC trk}(#Deltap_{t}/p_{t}^{MC}) [%]", "#sigma^{TPC trk}(#Deltap_{t}/p_{t}^{MC}) [%]"}
 ,{"1/p_{t}^{MC} [c/GeV]", "PULL: #mu_{1/p_{t}}^{TPC trk}", "PULL: #sigma_{1/p_{t}}^{TPC trk}"}
  // MC track HMPID
 ,{"tg(#theta)", "#mu_{z}^{trk} [#mum]", "#sigma_{z}^{trk} [#mum]"}
 ,{"tg(#theta)", "PULL: #mu_{z}^{trk}", "PULL: #sigma_{z}^{trk}"}
  // MC track in TRD
 ,{"tg(#phi)", "#mu_{y}^{trk} [#mum]", "#sigma_{y}^{trk} [#mum]"}
 ,{"tg(#phi)", "PULL: #mu_{y}^{trk}", "PULL: #sigma_{y}^{trk}"}
 ,{"tg(#theta)", "#mu_{z}^{Trk} [#mum]", "#sigma_{z}^{Trk} [#mum]"}
 ,{"tg(#theta)", "PULL: #mu_{z}^{Trk}", "PULL: #sigma_{z}^{Trk}"}
 ,{"tg(#phi)", "#mu_{#phi}^{Trk} [mrad]", "#sigma_{#phi}^{Trk} [mrad]"}
 ,{"tg(#phi)", "PULL: #mu_{snp}^{Trk}", "PULL: #sigma_{snp}^{Trk}"}
 ,{"tg(#theta)", "#mu_{#theta}^{Trk} [mrad]", "#sigma_{#theta}^{Trk} [mrad]"}
 ,{"tg(#theta)", "PULL: #mu_{tgl}^{Trk}", "PULL: #sigma_{tgl}^{Trk}"}
 ,{"p_{t}^{MC} [GeV/c]", "#mu^{Trk}(#Deltap_{t}/p_{t}^{MC}) [%]", "#sigma^{Trk}(#Deltap_{t}/p_{t}^{MC}) [%]"}
 ,{"1/p_{t}^{MC} [c/GeV]", "PULL: #mu_{1/p_{t}}^{Trk}", "PULL: #sigma_{1/p_{t}}^{Trk}"}
};

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
    Double_t dydx0 = tracklet->GetYref(1),
             dydx  = tracklet->GetYfit(1);
    Double_t dtgPhi = (dydx - dydx0)/(1.- dydx*dydx0);    h->Fill(dydx0, TMath::ATan(dtgPhi));
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
  TMatrixDSym COV(kNPAR); TVectorD PAR(kNPAR);
  Double_t *pc = &covR[0], *pp = &parR[0];
  for(Int_t ir=0; ir<kNPAR; ir++, pp++){
    PAR(ir) = (*pp);
    for(Int_t ic = 0; ic<=ir; ic++,pc++){ 
      COV(ir,ic) = (*pc); COV(ic,ir) = (*pc);
    }
  }
  PAR[4] = TMath::Abs(PAR[4]); // remove sign of pt !!

  //COV.Print(); PAR.Print();

//   COV(0,0) += cov[0]; 
//   COV(1,0) += cov[1]; COV(0, 1) += cov[1];
//   COV(1,1) += cov[2];
//   Double_t dy = parR[0] - tracklet->GetY(); 
//   Double_t dz = parR[1] - tracklet->GetZ(); 

  if(fDebugLevel>=1){
/*    Double_t dydx = 1./TMath::Sqrt(1.-parR[2]*parR[2]),
             pt   = 1./parR[4];
    Double_t c0 = covR[0]+cov[0],
             c1 = covR[1]+cov[1],
             c2 = covR[2]+cov[2];*/
    (*fDebugStream) << "trackIn"
      << "x="       << x
      << "P="       << &PAR
      << "C="       << &COV
      << "\n";

    Double_t y = tracklet->GetY(); 
    Double_t z = tracklet->GetZ(); 
    (*fDebugStream) << "trackletIn"
      << "y="       << y
      << "z="       << z
      << "Vy="      << cov[0]
      << "Cyz="     << cov[1]
      << "Vz="      << cov[2]
      << "\n";
  }


  if(!HasMCdata()) return h;
  UChar_t s;
  Float_t dx, pt0, x0=tracklet->GetX0(), y0, z0, dydx0, dzdx0;
  if(!fMC->GetDirections(x0, y0, z0, dydx0, dzdx0, pt0, s)) return h;
  // translate to reference radial position
  dx = x0 - x; y0 -= dx*dydx0; z0 -= dx*dzdx0;
  //Fill MC info
  TVectorD PARMC(kNPAR);
  PARMC[0]=y0; PARMC[1]=z0;
  PARMC[2]=dydx0/TMath::Sqrt(1.+dydx0*dydx0); PARMC[3]=dzdx0;
  PARMC[4]=1./pt0;

//   TMatrixDSymEigen eigen(COV);
//   TVectorD evals = eigen.GetEigenValues();
//   TMatrixDSym evalsm(kNPAR);
//   for(Int_t ir=0; ir<kNPAR; ir++) for(Int_t ic=0; ic<kNPAR; ic++) evalsm(ir,ic) = (ir==ic ? evals(ir): 0.);
//   TMatrixD evecs = eigen.GetEigenVectors();
//   TMatrixD sqrcov(evecs, TMatrixD::kMult, TMatrixD(evalsm, TMatrixD::kMult, evecs.T()));
  
  // fill histos
  TObjArray *arr = (TObjArray*)fContainer->At(kMCtrackTPC);
  // y resolution/pulls
  ((TH2I*)arr->At(0))->Fill(dydx0, PARMC[0]-PAR[0]);
  ((TH2I*)arr->At(1))->Fill(dydx0, (PARMC[0]-PAR[0])/TMath::Sqrt(COV(0,0)));
  // z resolution/pulls
  ((TH2I*)arr->At(2))->Fill(dzdx0, PARMC[1]-PAR[1]);
  ((TH2I*)arr->At(3))->Fill(dzdx0, (PARMC[1]-PAR[1])/TMath::Sqrt(COV(1,1)));
  // phi resolution/snp pulls
  ((TH2I*)arr->At(4))->Fill(dydx0, TMath::ASin(PARMC[2])-TMath::ASin(PAR[2]));
  ((TH2I*)arr->At(5))->Fill(dydx0, (PARMC[2]-PAR[2])/TMath::Sqrt(COV(2,2)));
  // theta resolution/tgl pulls
  ((TH2I*)arr->At(6))->Fill(dzdx0, TMath::ATan((PARMC[3]-PAR[3])/(1-PARMC[3]*PAR[3])));
  ((TH2I*)arr->At(7))->Fill(dzdx0, (PARMC[3]-PAR[3])/TMath::Sqrt(COV(3,3)));
  // pt resolution/1/pt pulls
  ((TH2I*)arr->At(8))->Fill(pt0, 1.-PARMC[4]/PAR[4]);
  ((TH2I*)arr->At(9))->Fill(PARMC[4], (PARMC[4]-PAR[4])/TMath::Sqrt(COV(4,4)));

  // fill debug for MC 
  if(fDebugLevel>=1){
    (*fDebugStream) << "trackInMC"
      << "P="   << &PARMC
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
  TObjArray *arr = 0x0;
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
    dydx = fTracklet->GetYref(1);
    dzdx = fTracklet->GetTgl();
    pt = TMath::Abs(fTracklet->GetPt());
    fTracklet->GetCovRef(covR);

    arr = (TObjArray*)fContainer->At(kMCtrack);
    // y resolution/pulls
    ((TH2I*)arr->At(0))->Fill(dydx0, dy);
    ((TH2I*)arr->At(1))->Fill(dydx0, dy/TMath::Sqrt(covR[0]));
    // z resolution/pulls
    ((TH2I*)arr->At(2))->Fill(dzdx0, dz);
    ((TH2I*)arr->At(3))->Fill(dzdx0, dz/TMath::Sqrt(covR[2]));
    // phi resolution/ snp pulls
    Double_t dtgp = (dydx - dydx0)/(1.- dydx*dydx0);
    ((TH2I*)arr->At(4))->Fill(dydx0, TMath::ATan(dtgp));
    //TODO ((TH2I*)arr->At(5))->Fill(dydx0, );
    // theta resolution/ tgl pulls
    Double_t dtgl = (dzdx - dzdx0)/(1.- dzdx*dzdx0);
    ((TH2I*)arr->At(6))->Fill(dzdx0, 
    TMath::ATan(dtgl));
    //TODO ((TH2I*)arr->At(7))->Fill(dydx0, );
    // pt resolution/ 1/pt pulls
    if(pdg!=kElectron && pdg!=kPositron){ 
      ((TH2I*)arr->At(8))->Fill(pt0, 1.-pt/pt0);
      //TODO ((TH2I*)arr->At(9))->Fill(1./pt0, (pt0/pt-1.)/TMath::Sqrt(covR[4]));
    }
    // Fill Debug stream for Kalman track
    if(fDebugLevel>=1){
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
    arr = (TObjArray*)fContainer->At(kMCtracklet);
    if(!rc){
      dy    = yt-y;

      Float_t dphi  = (dydx - dydx0);
      dphi /= 1.- dydx*dydx0;

      ((TH2I*)arr->At(0))->Fill(dydx0, dy);
      if(tt.GetS2Y()>0.) ((TH2I*)arr->At(1))->Fill(dydx0, dy/TMath::Sqrt(tt.GetS2Y()));
      ((TH2I*)arr->At(4))->Fill(dydx0, TMath::ATan(dphi));
    } else {
      // add tracklet residuals for z
      dz = zt-z;
      ((TH2I*)arr->At(2))->Fill(dzdx0, dz);
      if(tt.GetS2Z()>0.) ((TH2I*)arr->At(3))->Fill(dzdx0, dz/TMath::Sqrt(tt.GetS2Z()));
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
    //Double_t exb = AliTRDCommonParam::Instance()->GetOmegaTau(1.5);

    arr = (TObjArray*)fContainer->At(kMCcluster);
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
      if(q>20. && q<250.){ 
        ((TH2I*)arr->At(0))->Fill(dydx0, dy);
        ((TH2I*)arr->At(1))->Fill(dydx0, dy/TMath::Sqrt(c->GetSigmaY2()));
      }

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
  Float_t xy[4] = {0., 0., 0., 0.};
  if(!gPad){
    AliWarning("Please provide a canvas to draw results.");
    return kFALSE;
  }
  TList *l = 0x0;
  switch(ifig){
  case kCluster:
    xy[0] = -.3; xy[1] = -0.5; xy[2] = .3; xy[3] = 2.5;
    if(GetGraphPlot(&xy[0], ETRDresolutionPlot(ifig))) return kTRUE;
    break;
  case kTracklet:
    xy[0] = -.3; xy[1] = -0.5; xy[2] = .3; xy[3] = 2.5;
    if(GetGraphPlot(&xy[0], ETRDresolutionPlot(ifig))) return kTRUE;
    break;
  case kTrackletPhi:
    xy[0] = -.3; xy[1] = -0.5; xy[2] = .3; xy[3] = 2.5;
    if(GetGraphPlot(&xy[0], ETRDresolutionPlot(ifig))) return kTRUE;
    break;
  case kMCcluster:
    gPad->Divide(2, 1); l=gPad->GetListOfPrimitives(); 
    xy[0]=-.3; xy[1]=-50.; xy[2]=.3; xy[3]=650.;
    ((TVirtualPad*)l->At(0))->cd();
    if(!GetGraphPlot(&xy[0], ETRDresolutionPlot(ifig), 0)) break;
    xy[0] = -.3; xy[1] = -0.5; xy[2] = .3; xy[3] = 2.5;
    ((TVirtualPad*)l->At(1))->cd();
    if(!GetGraphPlot(&xy[0], ETRDresolutionPlot(ifig), 1)) break;
    return kTRUE;
  case 4: //kMCtracklet [y]
    gPad->Divide(2, 1); l=gPad->GetListOfPrimitives(); 
    xy[0]=-.2; xy[1]=-50.; xy[2]=.2; xy[3] =250.;
    ((TVirtualPad*)l->At(0))->cd();
    if(!GetGraphPlot(&xy[0], kMCtracklet, 0)) break;
    xy[0] = -.2; xy[1] = -0.5; xy[2] = .2; xy[3] = 2.5;
    ((TVirtualPad*)l->At(1))->cd();
    if(!GetGraphPlot(&xy[0], kMCtracklet, 1)) break;
    return kTRUE;
  case 5: //kMCtracklet [z]
    gPad->Divide(2, 1); l=gPad->GetListOfPrimitives(); 
    xy[0]=-1.; xy[1]=-100.; xy[2]=1.; xy[3] =2500.;
    ((TVirtualPad*)l->At(0))->cd();
    if(!GetGraphPlot(&xy[0], kMCtracklet, 2)) break;
    xy[0] = -1.; xy[1] = -0.5; xy[2] = 1.; xy[3] = 2.5;
    ((TVirtualPad*)l->At(1))->cd();
    if(!GetGraphPlot(&xy[0], kMCtracklet, 3)) break;
    return kTRUE;
  case 6: //kMCtracklet [phi]
    xy[0]=-.3; xy[1]=-3.; xy[2]=.3; xy[3] =25.;
    if(!GetGraphPlot(&xy[0], kMCtracklet, 4)) break;
    return kTRUE;
  case 7: //kMCtrack [y]
    gPad->Divide(2, 1); l=gPad->GetListOfPrimitives(); 
    xy[0]=-.2; xy[1]=-50.; xy[2]=.2; xy[3] =200.;
    ((TVirtualPad*)l->At(0))->cd();
    if(!GetGraphPlot(&xy[0], kMCtrack, 0)) break;
    xy[0] = -.2; xy[1] = -0.5; xy[2] = .2; xy[3] = 2.5;
    ((TVirtualPad*)l->At(1))->cd();
    if(!GetGraphPlot(&xy[0], kMCtrack, 1)) break;
    return kTRUE;
  case 8: //kMCtrack [pt]
    xy[0] = 0.; xy[1] = -0.5; xy[2] = 15.; xy[3] = 5.5;
    if(!GetGraphPlot(&xy[0], kMCtrack, 2)) break;
    return kTRUE;
  case 9: // kMCtrackTPC [y]
    gPad->Divide(2, 1); l=gPad->GetListOfPrimitives(); 
    xy[0]=-.25; xy[1]=-50.; xy[2]=.25; xy[3] =800.;
    ((TVirtualPad*)l->At(0))->cd();
    if(!GetGraphPlot(&xy[0], kMCtrackTPC, 0)) break;
    xy[0] = -.25; xy[1] = -0.5; xy[2] = .25; xy[3] = 2.5;
    ((TVirtualPad*)l->At(1))->cd();
    if(!GetGraphPlot(&xy[0], kMCtrackTPC, 1)) break;
    return kTRUE;
  case 10: // kMCtrackTPC [z]
    gPad->Divide(2, 1); l=gPad->GetListOfPrimitives(); 
    xy[0]=-1.; xy[1]=-500.; xy[2]=1.; xy[3] =800.;
    ((TVirtualPad*)l->At(0))->cd();
    if(!GetGraphPlot(&xy[0], kMCtrackTPC, 2)) break;
    xy[0] = -1.; xy[1] = -0.5; xy[2] = 1.; xy[3] = 2.5;
    ((TVirtualPad*)l->At(1))->cd();
    if(!GetGraphPlot(&xy[0], kMCtrackTPC, 3)) break;
    return kTRUE;
  case 11: // kMCtrackTPC [phi|snp]
    gPad->Divide(2, 1); l=gPad->GetListOfPrimitives(); 
    xy[0]=-.25; xy[1]=-0.5; xy[2]=.25; xy[3] =2.5;
    ((TVirtualPad*)l->At(0))->cd();
    if(!GetGraphPlot(&xy[0], kMCtrackTPC, 4)) break;
    //xy[0] = -.25; xy[1] = -0.5; xy[2] = .25; xy[3] = 2.5;
    ((TVirtualPad*)l->At(1))->cd();
    if(!GetGraphPlot(&xy[0], kMCtrackTPC, 5)) break;
    return kTRUE;
  case 12: // kMCtrackTPC [theta|tgl]
    gPad->Divide(2, 1); l=gPad->GetListOfPrimitives(); 
    xy[0]=-1.; xy[1]=-10.5; xy[2]=1.; xy[3] =20.5;
    ((TVirtualPad*)l->At(0))->cd();
    if(!GetGraphPlot(&xy[0], kMCtrackTPC, 6)) break;
    xy[0] = -1.; xy[1] = -0.5; xy[2] = 1.; xy[3] = 2.5;
    ((TVirtualPad*)l->At(1))->cd();
    if(!GetGraphPlot(&xy[0], kMCtrackTPC, 7)) break;
    return kTRUE;
  case 13: // kMCtrackTPC [pt]
    gPad->Divide(2, 1); l=gPad->GetListOfPrimitives(); 
    xy[0] = 0.; xy[1] = -.5; xy[2] = 15.; xy[3] = 2.5;
    ((TVirtualPad*)l->At(0))->cd();
    if(!GetGraphPlot(xy, AliTRDresolution::kMCtrackTPC, 8)) break;
    xy[0]=0.; xy[1]=-0.5; xy[2]=2.; xy[3] =2.5;
    ((TVirtualPad*)l->At(1))->cd();
    if(!GetGraphPlot(xy, AliTRDresolution::kMCtrackTPC, 9)) break;
    return kTRUE;
  case 14:  // kMCtrackHMPID [z]
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
  TGraphErrors *gm= 0x0, *gs= 0x0;
  if(!fGraphS && !fGraphM){ 
    TObjArray *aM(0x0), *aS(0x0);
    Int_t n = fContainer->GetEntriesFast();
    fGraphS = new TObjArray(n); fGraphS->SetOwner();
    fGraphM = new TObjArray(n); fGraphM->SetOwner();
    for(Int_t ig=0; ig<n; ig++){
      if(fNElements[ig]>1){
        fGraphM->AddAt(aM = new TObjArray(fNElements[ig]), ig);
        fGraphS->AddAt(aS = new TObjArray(fNElements[ig]), ig);
      } else {
        aM = fGraphM;aS = fGraphS;
      }
      for(Int_t ic=0; ic<fNElements[ig]; ic++){
        aS->AddAt(gs = new TGraphErrors(), fNElements[ig]>1?ic:ig);
        gs->SetMarkerStyle(23);
        gs->SetMarkerColor(kRed);
        gs->SetLineColor(kRed);
        gs->SetNameTitle(Form("s_%d%02d", ig, ic), "");

        aM->AddAt(gm = new TGraphErrors(), fNElements[ig]>1?ic:ig);
        gm->SetLineColor(kBlue);
        gm->SetMarkerStyle(7);
        gm->SetMarkerColor(kBlue);
        gm->SetNameTitle(Form("m_%d%02d", ig, ic), "");
      }
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
  Process(kCluster, -1, &f); fNRefFigures = 1;
  // Tracklet y residuals
  Process3D(kTracklet, &f); fNRefFigures = 2;
  // Tracklet phi residuals
  Process(kTrackletPhi, -1, &f, 1.e3); fNRefFigures = 3;

  if(!HasMCdata()) return kTRUE;


  //PROCESS MC RESIDUAL DISTRIBUTIONS

  // CLUSTER Y RESOLUTION/PULLS
  Process(kMCcluster, 0, &f, 1.e4);
  Process(kMCcluster, 1, &f);
  fNRefFigures = 4;

  // TRACKLET RESOLUTION/PULLS
  Process(kMCtracklet, 0, &f, 1.e4); // y
  Process(kMCtracklet, 1, &f);       // y pulls
  Process(kMCtracklet, 2, &f, 1.e4); // z
  Process(kMCtracklet, 3, &f);       // z pulls
  Process(kMCtracklet, 4, &f, 1.e3); // phi
  fNRefFigures = 7;

  // TRACK RESOLUTION/PULLS
  Process(kMCtrack, 0, &f, 1.e4);   // y
  Process(kMCtrack, 1, &f);         // y PULL
  Process(kMCtrack, 2, &f, 1.e4);   // z
  Process(kMCtrack, 3, &f);         // z PULL
  Process(kMCtrack, 4, &f, 1.e3);   // phi
  //Process(kMCtrack, 5, &f);         // snp PULL
  Process(kMCtrack, 6, &f, 1.e3);   // theta
  //Process(kMCtrack, 7, &f);         // tgl PULL
  Process(kMCtrack, 8, &f, 1.e2);   // pt resolution
  //Process(kMCtrack, 9, &f);         // 1/pt pulls
  fNRefFigures = 12;


  // TRACK TPC RESOLUTION/PULLS
  Process(kMCtrackTPC, 0, &f, 1.e4);// y resolution
  Process(kMCtrackTPC, 1, &f);      // y pulls
  Process(kMCtrackTPC, 2, &f, 1.e4);// z resolution
  Process(kMCtrackTPC, 3, &f);      // z pulls
  Process(kMCtrackTPC, 4, &f, 1.e3);// phi resolution
  Process(kMCtrackTPC, 5, &f);      // snp pulls
  Process(kMCtrackTPC, 6, &f, 1.e3);// theta resolution
  Process(kMCtrackTPC, 7, &f);      // tgl pulls
  Process(kMCtrackTPC, 8, &f, 1.e2);// pt resolution
  Process(kMCtrackTPC, 9, &f);      // 1/pt pulls
  fNRefFigures = 17;
  return kTRUE;

  // TRACK HMPID RESOLUTION/PULLS
  Process(kMCtrackHMPID, 0, &f, 1.e4); // z towards TOF
  Process(kMCtrackHMPID, 1, &f);       // z towards TOF
  fNRefFigures = 19;

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
    h = new TH2I("hCl", "Cluster Residuals", 21, -.33, .33, 100, -.5, .5);
    h->GetXaxis()->SetTitle("tg(#phi)");
    h->GetYaxis()->SetTitle("#Delta y [cm]");
    h->GetZaxis()->SetTitle("entries");
  } else h->Reset();
  fContainer->AddAt(h, kCluster);

  // tracklet to track residuals [2]
  if(!(h = (TH3I*)gROOT->FindObject("hTrklt"))){
    h = new TH3I("hTrklt", "Tracklet YZ Residuals", 100, -.5, .5, 100, -.5, .5, 21, -.33, .33);
    h->GetXaxis()->SetTitle("#Delta y [cm]");
    h->GetYaxis()->SetTitle("#Delta z [cm]");
    h->GetZaxis()->SetTitle("tg(#phi)");
  } else h->Reset();
  fContainer->AddAt(h, kTracklet);

  // tracklet to track residuals angular [2]
  if(!(h = (TH2I*)gROOT->FindObject("hTrkltPhi"))){
    h = new TH2I("hTrkltPhi", "Tracklet #phi Residuals (#Phi)", 21, -.33, .33, 100, -.5, .5);
    h->GetXaxis()->SetTitle("tg(#phi)");
    h->GetYaxis()->SetTitle("#Delta phi [rad]");
    h->GetZaxis()->SetTitle("entries");
  } else h->Reset();
  fContainer->AddAt(h, kTrackletPhi);


  // Resolution histos
  if(!HasMCdata()) return fContainer;

  // cluster y resolution [0]
  TObjArray *arr = 0x0;
  fContainer->AddAt(arr = new TObjArray(fNElements[kMCcluster]), kMCcluster);
  arr->SetName("Cl");
  if(!(h = (TH2I*)gROOT->FindObject("hMCcl"))){
    h = new TH2I("hMCcl", "Cluster Resolution", 31, -.48, .48, 100, -.3, .3);
    h->GetXaxis()->SetTitle("tg(#phi)");
    h->GetYaxis()->SetTitle("#Delta y [cm]");
    h->GetZaxis()->SetTitle("entries");
  } else h->Reset();
  arr->AddAt(h, 0);
  if(!(h = (TH2I*)gROOT->FindObject("hMCclPull"))){
    h = new TH2I("hMCclPull", "Cluster Pulls", 31, -.48, .48, 100, -4.5, 4.5);
    h->GetXaxis()->SetTitle("tg(#phi)");
    h->GetYaxis()->SetTitle("#Deltay/#sigma_{y}");
    h->GetZaxis()->SetTitle("entries");
  } else h->Reset();
  arr->AddAt(h, 1);


  // TRACKLET RESOLUTION
  fContainer->AddAt(arr = new TObjArray(fNElements[kMCtracklet]), kMCtracklet);
  arr->SetName("Trklt");
  // tracklet y resolution
  if(!(h = (TH2I*)gROOT->FindObject("hMCtrkltY"))){
    h = new TH2I("hMCtrkltY", "Tracklet Resolution (Y)", 31, -.48, .48, 100, -.2, .2);
    h->GetXaxis()->SetTitle("tg(#phi)");
    h->GetYaxis()->SetTitle("#Delta y [cm]");
    h->GetZaxis()->SetTitle("entries");
  } else h->Reset();
  arr->AddAt(h, 0);
  // tracklet y pulls
  if(!(h = (TH2I*)gROOT->FindObject("hMCtrkltYPull"))){
    h = new TH2I("hMCtrkltYPull", "Tracklet Pulls (Y)", 31, -.48, .48, 100, -4.5, 4.5);
    h->GetXaxis()->SetTitle("tg(#phi)");
    h->GetYaxis()->SetTitle("#Delta y / #sigma_{y}");
    h->GetZaxis()->SetTitle("entries");
  } else h->Reset();
  arr->AddAt(h, 1);
  // tracklet z resolution
  if(!(h = (TH2I*)gROOT->FindObject("hMCtrkltZ"))){
    h = new TH2I("hMCtrkltZ", "Tracklet Resolution (Z)", 50, -1., 1., 100, -1., 1.);
    h->GetXaxis()->SetTitle("tg(#theta)");
    h->GetYaxis()->SetTitle("#Delta z [cm]");
    h->GetZaxis()->SetTitle("entries");
  } else h->Reset();
  arr->AddAt(h, 2);
  // tracklet z pulls
  if(!(h = (TH2I*)gROOT->FindObject("hMCtrkltZPull"))){
    h = new TH2I("hMCtrkltZPull", "Tracklet Pulls (Z)", 31, -1., 1., 100, -3.5, 3.5);
    h->GetXaxis()->SetTitle("tg(#theta)");
    h->GetYaxis()->SetTitle("#Delta z / #sigma_{z}");
    h->GetZaxis()->SetTitle("entries");
  } else h->Reset();
  arr->AddAt(h, 3);
  // tracklet phi resolution
  if(!(h = (TH2I*)gROOT->FindObject("hMCtrkltPhi"))){
    h = new TH2I("hMCtrkltPhi", "Tracklet Resolution (#Phi)", 31, -.48, .48, 100, -.15, .15);
    h->GetXaxis()->SetTitle("tg(#phi)");
    h->GetYaxis()->SetTitle("#Delta #phi [rad]");
    h->GetZaxis()->SetTitle("entries");
  } else h->Reset();
  arr->AddAt(h, 4);


  // KALMAN TRACK RESOLUTION
  fContainer->AddAt(arr = new TObjArray(fNElements[kMCtrack]), kMCtrack);
  arr->SetName("Track");
  // Kalman track y resolution
  if(!(h = (TH2I*)gROOT->FindObject("hMCtrkY"))){
    h = new TH2I("hMCtrkY", "Track Y Resolution", 31, -.48, .48, 100, -.2, .2);
    h->GetXaxis()->SetTitle("tg(#phi)");
    h->GetYaxis()->SetTitle("#Delta y [cm]");
    h->GetZaxis()->SetTitle("entries");
  } else h->Reset();
  arr->AddAt(h, 0);
  // Kalman track y pulls
  if(!(h = (TH2I*)gROOT->FindObject("hMCtrkYPull"))){
    h = new TH2I("hMCtrkYPull", "Track Y Pulls", 31, -.48, .48, 100, -4., 4.);
    h->GetXaxis()->SetTitle("tg(#phi)");
    h->GetYaxis()->SetTitle("#Delta y / #sigma_{y}");
    h->GetZaxis()->SetTitle("entries");
  } else h->Reset();
  arr->AddAt(h, 1);
  // Kalman track Z
  if(!(h = (TH2I*)gROOT->FindObject("hMCtrkZ"))){
    h = new TH2I("hMCtrkZ", "Track Z Resolution", 20, -1., 1., 100, -1., 1.);
    h->GetXaxis()->SetTitle("tg(#theta)");
    h->GetYaxis()->SetTitle("#Delta z [cm]");
    h->GetZaxis()->SetTitle("entries");
  } else h->Reset();
  arr->AddAt(h, 2);
  // Kalman track Z pulls
  if(!(h = (TH2I*)gROOT->FindObject("hMCtrkZPull"))){
    h = new TH2I("hMCtrkZPull", "Track Z Pulls", 20, -1., 1., 100, -4.5, 4.5);
    h->GetXaxis()->SetTitle("tg(#theta)");
    h->GetYaxis()->SetTitle("#Delta z / #sigma_{z}");
    h->GetZaxis()->SetTitle("entries");
  } else h->Reset();
  arr->AddAt(h, 3);
  // Kalman track SNP
  if(!(h = (TH2I*)gROOT->FindObject("hMCtrkSNP"))){
    h = new TH2I("hMCtrkSNP", "Track Phi Resolution", 20, -.3, .3, 100, -.02, .02);
    h->GetXaxis()->SetTitle("tg(#phi)");
    h->GetYaxis()->SetTitle("#Delta #phi [rad]");
    h->GetZaxis()->SetTitle("entries");
  } else h->Reset();
  arr->AddAt(h, 4);
  // Kalman track SNP pulls
  if(!(h = (TH2I*)gROOT->FindObject("hMCtrkSNPPull"))){
    h = new TH2I("hMCtrkSNPPull", "Track SNP Pulls", 20, -.3, .3, 100, -4.5, 4.5);
    h->GetXaxis()->SetTitle("tg(#phi)");
    h->GetYaxis()->SetTitle("#Delta(sin(#phi)) / #sigma_{sin(#phi)}");
    h->GetZaxis()->SetTitle("entries");
  } else h->Reset();
  arr->AddAt(h, 5);
  // Kalman track TGL
  if(!(h = (TH2I*)gROOT->FindObject("hMCtrkTGL"))){
    h = new TH2I("hMCtrkTGL", "Track Theta Resolution", 20, -1., 1., 100, -.1, .1);
    h->GetXaxis()->SetTitle("tg(#theta)");
    h->GetYaxis()->SetTitle("#Delta#theta [rad]");
    h->GetZaxis()->SetTitle("entries");
  } else h->Reset();
  arr->AddAt(h, 6);
  // Kalman track TGL pulls
  if(!(h = (TH2I*)gROOT->FindObject("hMCtrkTGLPull"))){
    h = new TH2I("hMCtrkTGLPull", "Track TGL  Pulls", 20, -1., 1., 100, -4.5, 4.5);
    h->GetXaxis()->SetTitle("tg(#theta)");
    h->GetYaxis()->SetTitle("#Delta(tg(#theta)) / #sigma_{tg(#theta)}");
    h->GetZaxis()->SetTitle("entries");
  } else h->Reset();
  arr->AddAt(h, 7);
  // Kalman track Pt resolution
  if(!(h = (TH2I*)gROOT->FindObject("hMCtrkPt"))){
    h = new TH2I("hMCtrkPt", "Track Pt Resolution", 40, 0., 20., 150, -.15, .15);
    h->GetXaxis()->SetTitle("p_{t} [GeV/c]");
    h->GetYaxis()->SetTitle("#Delta p_{t}/p_{t}^{MC}");
    h->GetZaxis()->SetTitle("entries");
  } else h->Reset();
  arr->AddAt(h, 8);
  // Kalman track Pt pulls
  if(!(h = (TH2I*)gROOT->FindObject("hMCtrkPtPulls"))){
    h = new TH2I("hMCtrkPtPulls", "Track 1/Pt Pulls", 40, 0., 2., 100, -4., 4.);
    h->GetXaxis()->SetTitle("1/p_{t}^{MC} [c/GeV]");
    h->GetYaxis()->SetTitle("#Delta(1/p_{t})/#sigma(1/p_{t}) ");
    h->GetZaxis()->SetTitle("entries");
  } else h->Reset();
  arr->AddAt(h, 9);


  // TPC TRACK RESOLUTION
  fContainer->AddAt(arr = new TObjArray(fNElements[kMCtrackTPC]), kMCtrackTPC);
  arr->SetName("TrackTPC");
  // Kalman track Y
  if(!(h = (TH2I*)gROOT->FindObject("hMCtrkYIn"))){
    h = new TH2I("hMCtrkYIn", "Track[TPC] Y Resolution", 20, -.3, .3, 100, -.5, .5);
    h->GetXaxis()->SetTitle("tg(#phi)");
    h->GetYaxis()->SetTitle("#Delta y [cm]");
    h->GetZaxis()->SetTitle("entries");
  } else h->Reset();
  arr->AddAt(h, 0);
  // Kalman track Y pulls
  if(!(h = (TH2I*)gROOT->FindObject("hMCtrkYInPull"))){
    h = new TH2I("hMCtrkYInPull", "Track[TPC] Y Pulls", 20, -.3, .3, 100, -4.5, 4.5);
    h->GetXaxis()->SetTitle("tg(#phi)");
    h->GetYaxis()->SetTitle("#Delta y / #sigma_{y}");
    h->GetZaxis()->SetTitle("entries");
  } else h->Reset();
  arr->AddAt(h, 1);
  // Kalman track Z
  if(!(h = (TH2I*)gROOT->FindObject("hMCtrkZIn"))){
    h = new TH2I("hMCtrkZIn", "Track[TPC] Z Resolution", 20, -1., 1., 100, -1., 1.);
    h->GetXaxis()->SetTitle("tg(#theta)");
    h->GetYaxis()->SetTitle("#Delta z [cm]");
    h->GetZaxis()->SetTitle("entries");
  } else h->Reset();
  arr->AddAt(h, 2);
  // Kalman track Z pulls
  if(!(h = (TH2I*)gROOT->FindObject("hMCtrkZInPull"))){
    h = new TH2I("hMCtrkZInPull", "Track[TPC] Z Pulls", 20, -1., 1., 100, -4.5, 4.5);
    h->GetXaxis()->SetTitle("tg(#theta)");
    h->GetYaxis()->SetTitle("#Delta z / #sigma_{z}");
    h->GetZaxis()->SetTitle("entries");
  } else h->Reset();
  arr->AddAt(h, 3);
  // Kalman track SNP
  if(!(h = (TH2I*)gROOT->FindObject("hMCtrkSNPIn"))){
    h = new TH2I("hMCtrkSNPIn", "Track[TPC] Phi Resolution", 60, -.3, .3, 100, -.02, .02);
    h->GetXaxis()->SetTitle("tg(#phi)");
    h->GetYaxis()->SetTitle("#Delta #phi [rad]");
    h->GetZaxis()->SetTitle("entries");
  } else h->Reset();
  arr->AddAt(h, 4);
  // Kalman track SNP pulls
  if(!(h = (TH2I*)gROOT->FindObject("hMCtrkSNPInPull"))){
    h = new TH2I("hMCtrkSNPInPull", "Track[TPC] SNP Pulls", 60, -.3, .3, 100, -4.5, 4.5);
    h->GetXaxis()->SetTitle("tg(#phi)");
    h->GetYaxis()->SetTitle("#Delta(sin(#phi)) / #sigma_{sin(#phi)}");
    h->GetZaxis()->SetTitle("entries");
  } else h->Reset();
  arr->AddAt(h, 5);
  // Kalman track TGL
  if(!(h = (TH2I*)gROOT->FindObject("hMCtrkTGLIn"))){
    h = new TH2I("hMCtrkTGLIn", "Track[TPC] Theta Resolution", 60, -1., 1., 100, -.1, .1);
    h->GetXaxis()->SetTitle("tg(#theta)");
    h->GetYaxis()->SetTitle("#Delta#theta [rad]");
    h->GetZaxis()->SetTitle("entries");
  } else h->Reset();
  arr->AddAt(h, 6);
  // Kalman track TGL pulls
  if(!(h = (TH2I*)gROOT->FindObject("hMCtrkTGLInPull"))){
    h = new TH2I("hMCtrkTGLInPull", "Track[TPC] TGL  Pulls", 60, -1., 1., 100, -4.5, 4.5);
    h->GetXaxis()->SetTitle("tg(#theta)");
    h->GetYaxis()->SetTitle("#Delta(tg(#theta)) / #sigma_{tg(#theta)}");
    h->GetZaxis()->SetTitle("entries");
  } else h->Reset();
  arr->AddAt(h, 7);
  // Kalman track Pt resolution
  if(!(h = (TH2I*)gROOT->FindObject("hMCtrkPtIn"))){
    h = new TH2I("hMCtrkPtIn", "Track[TPC] Pt Resolution", 80, 0., 20., 150, -.15, .15);
    h->GetXaxis()->SetTitle("p_{t} [GeV/c]");
    h->GetYaxis()->SetTitle("#Delta p_{t}/p_{t}^{MC}");
    h->GetZaxis()->SetTitle("entries");
  } else h->Reset();
  arr->AddAt(h, 8);
  // Kalman track Pt pulls
  if(!(h = (TH2I*)gROOT->FindObject("hMCtrkPtInPulls"))){
    h = new TH2I("hMCtrkPtInPulls", "Track[TPC] 1/Pt Pulls", 80, 0., 2., 100, -4., 4.);
    h->GetXaxis()->SetTitle("1/p_{t}^{MC} [c/GeV]");
    h->GetYaxis()->SetTitle("#Delta(1/p_{t})/#sigma(1/p_{t}) ");
    h->GetZaxis()->SetTitle("entries");
  } else h->Reset();
  arr->AddAt(h, 9);



  // Kalman track Z resolution [OUT]
  fContainer->AddAt(arr = new TObjArray(fNElements[kMCtrackHMPID]), kMCtrackHMPID);
  arr->SetName("TrackHMPID");
  if(!(h = (TH2I*)gROOT->FindObject("hMCtrkZOut"))){
    h = new TH2I("hMCtrkZOut", "Track[TOF] Z Resolution", 20, -1., 1., 100, -1., 1.);
    h->GetXaxis()->SetTitle("tg(#theta)");
    h->GetYaxis()->SetTitle("#Delta z [cm]");
    h->GetZaxis()->SetTitle("entries");
  } else h->Reset();
  arr->AddAt(h, 0);
  // Kalman track Z pulls
  if(!(h = (TH2I*)gROOT->FindObject("hMCtrkZOutPull"))){
    h = new TH2I("hMCtrkZOutPull", "Track[TOF] Z Pulls", 20, -1., 1., 100, -4.5, 4.5);
    h->GetXaxis()->SetTitle("tg(#theta)");
    h->GetYaxis()->SetTitle("#Delta z / #sigma_{z}");
    h->GetZaxis()->SetTitle("entries");
  } else h->Reset();
  arr->AddAt(h, 1);

  return fContainer;
}


//________________________________________________________
Bool_t AliTRDresolution::Process(ETRDresolutionPlot plot, Int_t idx, TF1 *f, Float_t k)
{
  if(!fContainer || !fGraphS || !fGraphM) return kFALSE;
  Bool_t kBUILD = kFALSE;
  if(!f){ 
    f = new TF1("f1", "gaus", -.5, .5);
    kBUILD = kTRUE;
  }

  // retrive containers
  TH2I *h2 = idx<0 ? (TH2I*)(fContainer->At(plot)) : (TH2I*)((TObjArray*)(fContainer->At(plot)))->At(idx);
  if(!h2) return kFALSE;
  TGraphErrors *gm = idx<0 ? (TGraphErrors*)fGraphM->At(plot) : (TGraphErrors*)((TObjArray*)(fGraphM->At(plot)))->At(idx);
  if(!gm) return kFALSE;
  if(gm->GetN()) for(Int_t ip=gm->GetN(); ip--;) gm->RemovePoint(ip);

  TGraphErrors *gs = idx<0 ? (TGraphErrors*)fGraphS->At(plot) : (TGraphErrors*)((TObjArray*)(fGraphS->At(plot)))->At(idx);
  if(!gs) return kFALSE;
  if(gs->GetN()) for(Int_t ip=gs->GetN(); ip--;) gs->RemovePoint(ip);

  Char_t pn[10]; sprintf(pn, "p%02d%02d", plot, idx);
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
Bool_t AliTRDresolution::GetGraphPlot(Float_t *bb, ETRDresolutionPlot ip, Int_t idx)
{
  if(!fGraphS || !fGraphM) return kFALSE;
  TGraphErrors *gm = idx<0 ? (TGraphErrors*)fGraphM->At(ip) : (TGraphErrors*)((TObjArray*)(fGraphM->At(ip)))->At(idx);
  if(!gm) return kFALSE;
  TGraphErrors *gs = idx<0 ? (TGraphErrors*)fGraphS->At(ip) : (TGraphErrors*)((TObjArray*)(fGraphS->At(ip)))->At(idx);
  if(!gs) return kFALSE;
  gs->Draw("apl"); gm->Draw("pl");

  //printf("bb[%f %f %f %f]\n", bb[0], bb[1], bb[2], bb[3]);

  // axis range
  TAxis *ax = 0x0;
  ax = gs->GetHistogram()->GetYaxis();
  ax->SetRangeUser(bb[1], bb[3]);
  ax = gs->GetHistogram()->GetXaxis();
  ax->SetRangeUser(bb[0], bb[2]);

  // axis titles
  Int_t nref = 0;
  for(Int_t jp=0; jp<(Int_t)ip; jp++) nref+=fNElements[jp];
  UChar_t jdx = idx<0?0:idx;
  for(Int_t jc=0; jc<TMath::Min(jdx,fNElements[ip]-1); jc++) nref++;
  Char_t **at = fAxTitle[nref];
  ax->SetTitle(*at);ax->CenterTitle();
  TGaxis *gax = 0x0;
  gax = new TGaxis(bb[2], bb[1], bb[2], bb[3], bb[1], bb[3], 510, "+L");
  gax->SetLineColor(kBlue);gax->SetLineWidth(2);gax->SetTextColor(kBlue);
  //gax->SetVertical();
  gax->CenterTitle();
  gax->SetTitle(*(++at)); gax->Draw();

  gax = new TGaxis(bb[0], bb[1], bb[0], bb[3], bb[1], bb[3], 510, "");
  gax->SetLineColor(kRed);gax->SetLineWidth(2);gax->SetTextColor(kRed);
  /*gax->SetVertical();*/gax->CenterTitle();
  gax->SetTitle(*(++at)); gax->Draw();
  //return kTRUE;

  // bounding box
  TBox *b = new TBox(-.15, bb[1], .15, bb[3]);
  b->SetFillStyle(3002);b->SetLineColor(0);
  b->SetFillColor(ip<=Int_t(kMCcluster)?kGreen:kBlue);
  b->Draw();
  return kTRUE;
}

// gSystem->Load("libANALYSIS.so")
// gSystem->Load("libTRDqaRec.so")
// AliTRDresolution res
void AliTRDresolution::DumpAxTitle(Int_t ip, Int_t idx)
{
  Int_t nref = 0;
  for(Int_t jp=0; jp<(Int_t)ip; jp++) nref+=fNElements[jp];
  UChar_t jdx = idx<0?0:idx;
  printf("ref entry* %d jdx[%d]\n", nref, jdx);
  for(Int_t jc=0; jc<TMath::Min(jdx,fNElements[ip]-1); jc++) nref++;
  printf("ref entry %d\n", nref);
  Char_t **at = fAxTitle[nref];

  printf("%s || ", (*at));
  at++; printf("%s || ", (*at));
  at++; printf("%s\n", (*at));
}


//________________________________________________________
void AliTRDresolution::SetRecoParam(AliTRDrecoParam *r)
{

  fReconstructor->SetRecoParam(r);
}
