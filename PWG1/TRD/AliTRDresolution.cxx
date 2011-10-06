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

#include "info/AliTRDclusterInfo.h"

ClassImp(AliTRDresolution)
//ClassImp(AliTRDresolution::AliTRDresolutionProjection)

Int_t const   AliTRDresolution::fgkNbins[kNdim] = {
  Int_t(kNbunchCross)/*bc*/,
  180/*phi*/,
  50/*eta*/,
  50/*dy*/,
  40/*dphi*/,
  50/*dz*/,
  Int_t(kNcharge)*AliPID::kSPECIES+1/*chg*species*/,
  kNpt/*pt*/
};  //! no of bins/projection
Double_t const AliTRDresolution::fgkMin[kNdim] = {
  -0.5,
  -TMath::Pi(),
  -1.,
  -1.5,
  -10.,
  -2.5,
  -AliPID::kSPECIES-0.5,
  -0.5
};    //! low limits for projections
Double_t const AliTRDresolution::fgkMax[kNdim] = {
  Int_t(kNbunchCross)-0.5,
  TMath::Pi(),
  1.,
  1.5,
  10.,
  2.5,
  AliPID::kSPECIES+0.5,
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

Int_t const AliTRDresolution::fgkNproj[kNclasses] = {
  48, 72, 8, 5,
  2, 5, 11, 11, 11
};
Char_t const * AliTRDresolution::fgPerformanceName[kNclasses] = {
    "Cluster2Track"
    ,"Tracklet2Track"
    ,"Tracklet2TRDin" 
    ,"Tracklet2TRDout" 
    ,"Cluster2MC"
    ,"Tracklet2MC"
    ,"TRDin2MC"
    ,"TRDout2MC"
    ,"TRD2MC"
};
Float_t AliTRDresolution::fgPtBin[kNpt+1];

//________________________________________________________
AliTRDresolution::AliTRDresolution()
  :AliTRDrecoTask()
  ,fIdxPlot(0)
  ,fIdxFrame(0)
  ,fPtThreshold(1.)
  ,fDyRange(0.75)
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
}

//________________________________________________________
AliTRDresolution::AliTRDresolution(char* name, Bool_t xchange)
  :AliTRDrecoTask(name, "TRD spatial and momentum resolution")
  ,fIdxPlot(0)
  ,fIdxFrame(0)
  ,fPtThreshold(1.)
  ,fDyRange(0.75)
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
  if(TMath::Abs(fkESD->GetTOFbc()) > 1){
    AliDebug(4, Form("Track with BC_index[%d] not used.", fkESD->GetTOFbc()));
    return NULL;
  }
  if(fPt<fPtThreshold){
    AliDebug(4, Form("Track with pt[%6.4f] under threshold.", fPt));
    return NULL;
  }
  THnSparse *H(NULL);
  if(!fContainer || !(H = ((THnSparse*)fContainer->At(kCluster)))){
    AliWarning("No output container defined.");
    return NULL;
  }

  AliTRDgeometry *geo(AliTRDinfoGen::Geometry());
  Double_t val[kNdim]; //Float_t exb, vd, t0, s2, dl, dt;
  TObjArray     *clInfoArr(NULL);
  AliTRDseedV1  *fTracklet(NULL);
  AliTRDcluster *c(NULL), *cc(NULL);
  for(Int_t ily=0; ily<AliTRDgeometry::kNlayer; ily++){
    if(!(fTracklet = fkTrack->GetTracklet(ily))) continue;
    if(!fTracklet->IsOK()) continue;
    //fTracklet->GetCalibParam(exb, vd, t0, s2, dl, dt);
    val[kBC]  = ily;
    val[kPhi] = fPhi;
    val[kEta] = fEta;
    val[kPt]  = TMath::ATan(fTracklet->GetYref(1))*TMath::RadToDeg();
    Float_t corr = 1./TMath::Sqrt(1.+fTracklet->GetYref(1)*fTracklet->GetYref(1)+fTracklet->GetZref(1)*fTracklet->GetZref(1));
    Int_t row0(-1);
    Float_t padCorr(fTracklet->GetTilt()*fTracklet->GetPadLength());
    fTracklet->ResetClusterIter(kTRUE);
    while((c = fTracklet->NextCluster())){
      Float_t xc(c->GetX()),
              q(c->GetQ());
      Int_t tb(c->GetLocalTimeBin());
      if(row0<0) row0 = c->GetPadRow();

      val[kYrez] = c->GetY() + padCorr*(c->GetPadRow() - row0) -fTracklet->GetYat(xc);
      val[kPrez] = fTracklet->GetX0()-xc;
      val[kZrez] = 0.; Int_t ic(0);
      if((cc = fTracklet->GetClusters(tb-1))) {val[kZrez] += cc->GetQ(); ic++;}
      if((cc = fTracklet->GetClusters(tb-2))) {val[kZrez] += cc->GetQ(); ic++;}
      if(ic) val[kZrez] /= (ic*q);
      val[kSpeciesChgRC]= fTracklet->IsRowCross()?0.:(TMath::Max(q*corr, Float_t(3.)));
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

  return NULL;//H->Projection(kEta, kPhi);
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
  if(TMath::Abs(fkESD->GetTOFbc())>1){
    AliDebug(4, Form("Track with BC_index[%d] not used.", fkESD->GetTOFbc()));
    return NULL;
  }
  THnSparse *H(NULL);
  if(!fContainer || !(H = (THnSparse*)fContainer->At(kTracklet))){
    AliWarning("No output container defined.");
    return NULL;
  }
//  return NULL;
  Double_t val[kNdim+1];
  AliTRDseedV1 *fTracklet(NULL);
  for(Int_t il(0); il<AliTRDgeometry::kNlayer; il++){
    if(!(fTracklet = fkTrack->GetTracklet(il))) continue;
    if(!fTracklet->IsOK()) continue;
    val [kBC] = il; 
    val[kPhi] = fPhi;
    val[kEta] = fEta;
    val[kSpeciesChgRC]= fSpecies;
    val[kPt]  = GetPtBin(fTracklet->GetMomentum());
    Double_t dyt(fTracklet->GetYref(0) - fTracklet->GetYfit(0)),
             dzt(fTracklet->GetZref(0) - fTracklet->GetZfit(0)),
             dydx(fTracklet->GetYfit(1)),
             tilt(fTracklet->GetTilt());
    // correct for tilt rotation
    val[kYrez] = dyt - dzt*tilt;
    val[kZrez] = dzt + dyt*tilt;
    dydx+= tilt*fTracklet->GetZref(1);
    val[kPrez] = TMath::ATan((fTracklet->GetYref(1) - dydx)/(1.+ fTracklet->GetYref(1)*dydx)) * TMath::RadToDeg();
    if(fTracklet->IsRowCross()){
      val[kSpeciesChgRC]= 0.;
//      val[kPrez] = fkTrack->Charge(); // may be better defined
    }/* else {
      Float_t exb, vd, t0, s2, dl, dt;
      fTracklet->GetCalibParam(exb, vd, t0, s2, dl, dt);
      val[kZrez] = TMath::ATan((fTracklet->GetYref(1) - exb)/(1+fTracklet->GetYref(1)*exb));
    }*/
    val[kNdim] = fTracklet->GetdQdl();
    if(DebugLevel()>=1) H->Fill(val);

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
      Int_t ncl(fTracklet->GetN()),
            det(fTracklet->GetDetector());
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
  return NULL;//H->Projection(kEta, kPhi);
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
    return NULL;
  }
  // check radial position
  Double_t x = tin->GetX();
  if(TMath::Abs(x-fTracklet->GetX())>1.e-3){
    AliDebug(1, Form("Tracklet did not match Track. dx[cm]=%+4.1f", x-fTracklet->GetX()));
    return NULL;
  }
  //printf("USE y[%+f] dydx[%+f]\n", fTracklet->GetYfit(0), fTracklet->GetYfit(1));

  Int_t bc(TMath::Abs(fkESD->GetTOFbc())%2);
  const Double_t *parR(tin->GetParameter());
  Double_t dyt(parR[0] - fTracklet->GetYfit(0)), dzt(parR[1] - fTracklet->GetZfit(0)),
            phit(fTracklet->GetYfit(1)),
            tilt(fTracklet->GetTilt());

  // correct for tilt rotation
  Double_t dy  = dyt - dzt*tilt,
           dz  = dzt + dyt*tilt;
  phit       += tilt*parR[3];
  Double_t dphi = TMath::ASin(parR[2])-TMath::ATan(phit);

  Double_t val[kNdim];
  val[kBC]          = bc;
  val[kPhi]         = fPhi;
  val[kEta]         = fEta;
  val[kSpeciesChgRC]= fTracklet->IsRowCross()?0:fSpecies;
  val[kPt]          = GetPtBin(fPt);
  val[kYrez]        = dy;
  val[kZrez]        = dz;
  val[kPrez]        = dphi*TMath::RadToDeg();
  H->Fill(val);
  if(DebugLevel()>=3){
    (*DebugStream()) << "trackIn"
      <<"tracklet.="  << fTracklet
      <<"trackIn.="   << tin
      << "\n";
  }

  if(!HasMCdata()) return NULL; // H->Projection(kEta, kPhi);
  if(!(H = (THnSparseI*)fContainer->At(kMCtrackIn))) {
    AliError(Form("Missing container @ %d", Int_t(kMCtrackIn)));
    return NULL;
  }

  // get MC info
  UChar_t s;
  Float_t pt0, eta, x0=fTracklet->GetX0(), y0, z0, dydx0, dzdx0;
  if(!fkMC->GetDirections(x0, y0, z0, dydx0, dzdx0, pt0, eta, s)) return NULL;
  dyt = y0 - fTracklet->GetYfit(0);
  dzt = z0 - fTracklet->GetZfit(0);
  phit= fTracklet->GetYfit(1) + tilt*dzdx0;
  Float_t phi = TMath::ATan2(y0, x0);
  dy  = dyt - dzt*tilt;
  dz  = dzt + dyt*tilt;
  dphi= TMath::ASin(dydx0)-TMath::ATan(phit);

  Int_t pdg = fkMC->GetPDG(),
        sIdx(AliTRDpidUtil::Pdg2Pid(TMath::Abs(pdg))+1), // species index
        sign(0);
  if(!fDBPDG) fDBPDG=TDatabasePDG::Instance();
  TParticlePDG *ppdg(fDBPDG->GetParticle(pdg));
  if(ppdg) sign = ppdg->Charge() > 0. ? 1 : -1;


  val[kBC]          = (bc>=kNbunchCross)?(kNbunchCross-1):bc;
  val[kPhi]         = phi;
  val[kEta]         = eta;
  val[kSpeciesChgRC]= fTracklet->IsRowCross()?0:sign*(sIdx+1);
  val[kPt]          = GetPtBin(pt0);
  val[kYrez]        = dy;
  val[kZrez]        = dz;
  val[kPrez]        = dphi*TMath::RadToDeg();
  H->Fill(val);

  return NULL; //H->Projection(kEta, kPhi);
}

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

//________________________________________________________
TH1* AliTRDresolution::PlotMC(const AliTRDtrackV1 *track)
{
  //
  // Plot MC distributions
  //

  if(!HasMCdata()){ 
    AliDebug(2, "No MC defined. Results will not be available.");
    return NULL;
  }
  if(track) fkTrack = track;
  if(!fkTrack){
    AliDebug(4, "No Track defined.");
    return NULL;
  }
  if(!fContainer){
    AliWarning("No output container defined.");
    return NULL;
  }
  // retriev track characteristics
  Int_t pdg = fkMC->GetPDG(),
        sIdx(AliTRDpidUtil::Pdg2Pid(TMath::Abs(pdg))+1), // species index
        sign(0),
        sgm[3],
        label(fkMC->GetLabel()),
        fSegmentLevel(0);
  if(!fDBPDG) fDBPDG=TDatabasePDG::Instance();
  TParticlePDG *ppdg(fDBPDG->GetParticle(pdg));
  if(ppdg) sign = ppdg->Charge() > 0. ? 1 : -1;

  TObjArray *arr(NULL);TH1 *h(NULL);
  AliTRDgeometry *geo(AliTRDinfoGen::Geometry());
  AliTRDseedV1 *fTracklet(NULL); TObjArray *clInfoArr(NULL);
  UChar_t s;
  Double_t xAnode, x, y, z, pt, dydx, dzdx, dzdl;
  Float_t pt0, x0, y0, z0, dx, dy, dz, dydx0, dzdx0;
  Double_t covR[7]/*, cov[3]*/;
  
  if(DebugLevel()>=3){
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
  }
  for(Int_t ily=0; ily<AliTRDgeometry::kNlayer; ily++){
    if(!(fTracklet = fkTrack->GetTracklet(ily)))/* ||
       !fTracklet->IsOK())*/ continue;

    sgm[2] = fTracklet->GetDetector();
    sgm[0] = AliTRDgeometry::GetSector(sgm[2]);
    sgm[1] = sgm[0] * AliTRDgeometry::kNstack + AliTRDgeometry::GetStack(sgm[2]);
    Double_t tilt(fTracklet->GetTilt())
            ,t2(tilt*tilt)
            ,corr(1./(1. + t2))
            ,cost(TMath::Sqrt(corr));
    x0  = fTracklet->GetX0();
    //radial shift with respect to the MC reference (radial position of the pad plane)
    x= fTracklet->GetX();
    Bool_t rc(fTracklet->IsRowCross()); Float_t eta;
    if(!fkMC->GetDirections(x0, y0, z0, dydx0, dzdx0, pt0, eta, s)) continue;
    xAnode  = fTracklet->GetX0();

    // MC track position at reference radial position
    dx  = x0 - x;
    if(DebugLevel()>=4){
      (*DebugStream()) << "MC"
        << "det="     << sgm[2]
        << "pdg="     << pdg
        << "sgn="     << sign
        << "pt="      << pt0
        << "x="       << x0
        << "y="       << y0
        << "z="       << z0
        << "dydx="    << dydx0
        << "dzdx="    << dzdx0
        << "\n";
    }
    Float_t ymc = y0 - dx*dydx0;
    Float_t zmc = z0 - dx*dzdx0;
    //p = pt0*TMath::Sqrt(1.+dzdx0*dzdx0); // pt -> p

    // Kalman position at reference radial position
    dx = xAnode - x;
    dydx = fTracklet->GetYref(1);
    dzdx = fTracklet->GetZref(1);
    dzdl = fTracklet->GetTgl();
    y  = fTracklet->GetYref(0) - dx*dydx;
    dy = y - ymc;
    z  = fTracklet->GetZref(0) - dx*dzdx;
    dz = z - zmc;
    pt = TMath::Abs(fTracklet->GetPt());
    fTracklet->GetCovRef(covR);

    arr = (TObjArray*)((TObjArray*)fContainer->At(kMCtrack))->At(ily);
    // y resolution/pulls
    if(pt0>fPtThreshold) ((TH3S*)arr->At(0))->Fill(dydx0, dy, sgm[fSegmentLevel]);
    ((TH3S*)arr->At(1))->Fill(sgm[fSegmentLevel], dy/TMath::Sqrt(covR[0]), dz/TMath::Sqrt(covR[2]));
    // z resolution/pulls
    ((TH2S*)arr->At(2))->Fill(dzdx0, dz);
    ((TH3S*)arr->At(3))->Fill(dzdx0, dz/TMath::Sqrt(covR[2]), 0);
    // phi resolution/ snp pulls
    Double_t dtgp = (dydx - dydx0)/(1.- dydx*dydx0);
    ((TH2I*)arr->At(4))->Fill(dydx0, TMath::ATan(dtgp));
    Double_t dsnp = dydx/TMath::Sqrt(1.+dydx*dydx) - dydx0/TMath::Sqrt(1.+dydx0*dydx0);
    ((TH2I*)arr->At(5))->Fill(dydx0, dsnp/TMath::Sqrt(covR[3]));
    // theta resolution/ tgl pulls
    Double_t dzdl0 = dzdx0/TMath::Sqrt(1.+dydx0*dydx0),
              dtgl = (dzdl - dzdl0)/(1.- dzdl*dzdl0);
    ((TH2I*)arr->At(6))->Fill(dzdl0, 
    TMath::ATan(dtgl));
    ((TH2I*)arr->At(7))->Fill(dzdl0, (dzdl - dzdl0)/TMath::Sqrt(covR[4]));
    // pt resolution  \\ 1/pt pulls \\ p resolution for PID
    Double_t p0 = TMath::Sqrt(1.+ dzdl0*dzdl0)*pt0,
             p  = TMath::Sqrt(1.+ dzdl*dzdl)*pt;
    ((TH3S*)((TObjArray*)arr->At(8)))->Fill(pt0, pt/pt0-1., sign*sIdx);
    ((TH3S*)((TObjArray*)arr->At(9)))->Fill(1./pt0, (1./pt-1./pt0)/TMath::Sqrt(covR[6]), sign*sIdx);
    ((TH3S*)((TObjArray*)arr->At(10)))->Fill(p0, p/p0-1., sign*sIdx);

    // Fill Debug stream for Kalman track
    if(DebugLevel()>=4){
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

    // recalculate tracklet based on the MC info
    AliTRDseedV1 tt(*fTracklet);
    tt.SetZref(0, z0 - (x0-xAnode)*dzdx0);
    tt.SetZref(1, dzdx0); 
    tt.SetReconstructor(AliTRDinfoGen::Reconstructor());
    tt.Fit(1);
    x= tt.GetX();y= tt.GetY();z= tt.GetZ();
    dydx = tt.GetYfit(1);
    dx = x0 - x;
    ymc = y0 - dx*dydx0;
    zmc = z0 - dx*dzdx0;
    dy = y-ymc;
    dz = z-zmc;
    Float_t dphi  = (dydx - dydx0);
    dphi /= (1.- dydx*dydx0);

    // add tracklet residuals for y and dydx
    arr = (TObjArray*)fContainer->At(kMCtracklet);

    if(pt0>fPtThreshold) ((TH3S*)arr->At(0))->Fill(dydx0, dy, sgm[fSegmentLevel]);
    if(tt.GetS2Y()>0. && tt.GetS2Z()>0.) ((TH3S*)arr->At(1))->Fill(sgm[fSegmentLevel], dy/TMath::Sqrt(tt.GetS2Y()), dz/TMath::Sqrt(tt.GetS2Z()));
    ((TH3S*)arr->At(2))->Fill(dzdl0, dz, rc);
    if(tt.GetS2Z()>0.) ((TH3S*)arr->At(3))->Fill(dzdl0, dz/TMath::Sqrt(tt.GetS2Z()), rc);
    ((TH2I*)arr->At(4))->Fill(dydx0, TMath::ATan(dphi));
  
    // Fill Debug stream for tracklet
    if(DebugLevel()>=4){
      Float_t s2y = tt.GetS2Y();
      Float_t s2z = tt.GetS2Z();
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

    AliTRDpadPlane *pp = geo->GetPadPlane(ily, AliTRDgeometry::GetStack(sgm[2]));
    Float_t zr0 = pp->GetRow0() + pp->GetAnodeWireOffset();
    //Double_t exb = AliTRDCommonParam::Instance()->GetOmegaTau(1.5);

    arr = (TObjArray*)fContainer->At(kMCcluster);
    AliTRDcluster *c = NULL;
    tt.ResetClusterIter(kFALSE);
    while((c = tt.PrevCluster())){
      Float_t  q = TMath::Abs(c->GetQ());
      x = c->GetX();//+fXcorr[c->GetDetector()][c->GetLocalTimeBin()]; y = c->GetY();z = c->GetZ();
      dx = x0 - x; 
      ymc= y0 - dx*dydx0;
      zmc= z0 - dx*dzdx0;
      dy = cost*(y - ymc - tilt*(z-zmc));
      dz = cost*(z - zmc + tilt*(y-ymc));
      
      // Fill Histograms
      if(q>20. && q<250. && pt0>fPtThreshold && c->IsInChamber()){ 
        ((TH3S*)arr->At(0))->Fill(dydx0, dy, sgm[fSegmentLevel]);
        ((TH3S*)arr->At(1))->Fill(sgm[fSegmentLevel], dy/TMath::Sqrt(c->GetSigmaY2()), dz/TMath::Sqrt(c->GetSigmaZ2()));
      }

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
  while(ipt<AliTRDresolution::kNpt){
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
  for(Int_t j(0); j<=kNpt; j++){
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
  const Int_t nClViews(8);
  const Char_t *vClName[nClViews] = {"HClY", "HClYn", "HClYp", "HClQn", "HClQp", "HClYXTCp", "HClYXTCn", "HClYXPh"};
  const UChar_t vClOpt[nClViews] = {1, 1, 1, 0, 0, 0, 0, 0};
  if((arr = (TObjArray*)fProj->At(kCluster))){
    for(Int_t iview(0); iview<nClViews; iview++){
      cOut = new TCanvas(Form("TRDsummary%s_Cl%02d", GetName(), iview), "Cluster Resolution", 1024, 768);
      cOut->Divide(3,2, 1.e-5, 1.e-5);
      Int_t nplot(0);
      for(Int_t iplot(0); iplot<6; iplot++){
        p=cOut->cd(iplot+1);    p->SetRightMargin(0.1572581);p->SetTopMargin(0.08262712);
        if(!(h2 = (TH2*)arr->FindObject(Form("%s%d_2D", vClName[iview], iplot)))) continue;
        nplot++;
        if(vClOpt[iview]==0) h2->Draw("colz");
        else if(vClOpt[iview]==1) DrawSigma(h2, 1.e4, 2.e2, 6.e2, "#sigma(#Deltay) [#mum]");
      }
      if(nplot) cOut->SaveAs(Form("%s.gif", cOut->GetName()));
      else delete cOut;
    }
  }

  // tracklet systematic
  const Int_t nTrkltViews(10);
  const Char_t *vTrkltName[nTrkltViews] = {"HTrkltY", "HTrkltYn", "HTrkltYp", "HTrkltPhn", "HTrkltPhp", "HTrkltZ", "HTrkltQn", "HTrkltQp", "HTrkltPn", "HTrkltPp"};
  if((arr = (TObjArray*)fProj->At(kTracklet))){
    for(Int_t iview(0); iview<nTrkltViews; iview++){
      cOut = new TCanvas(Form("TRDsummary%s_Trklt%02d", GetName(), iview), "Tracklet Resolution", 1024, 768);
      cOut->Divide(3,2, 1.e-5, 1.e-5);
      Int_t nplot(0);
      for(Int_t iplot(0); iplot<6; iplot++){
        p=cOut->cd(iplot+1); p->SetRightMargin(0.1572581); p->SetTopMargin(0.08262712);
        if(!(h2 = (TH2*)arr->FindObject(Form("%s%d_2D", vTrkltName[iview], iplot)))) continue;
        h2->Draw("colz"); nplot++;
      }
      if(nplot) cOut->SaveAs(Form("%s.gif", cOut->GetName()));
      else delete cOut;
    }
  }
  // trackIn systematic
  const Char_t *hname[] = {"HTrkInY", "HTrkInYn", "HTrkInYp", "HTrkInZ", "HTrkInPhn", "HTrkInPhp"};
  if((arr = (TObjArray*)fProj->At(kTrackIn))){
    cOut = new TCanvas(Form("TRDsummary%s_TrkIn", GetName()), "Track IN Resolution", 1024, 768);
    cOut->Divide(3,2, 1.e-5, 1.e-5);
    Int_t nplot(0);
    for(Int_t iplot(0); iplot<6; iplot++){
      p=cOut->cd(iplot+1);    p->SetRightMargin(0.1572581);p->SetTopMargin(0.08262712);
        if(!(h2 = (TH2*)arr->FindObject(Form("%s_2D", hname[iplot])))) continue;
        h2->Draw("colz"); nplot++;
    }
    if(nplot) cOut->SaveAs(Form("%s.gif", cOut->GetName()));
    else delete cOut;
  }
  gStyle->SetPalette(1);
}

//________________________________________________________
void AliTRDresolution::DrawSigma(TH2 *h2, Float_t scale, Float_t m, Float_t M, const Char_t *title)
{
  // Draw error bars scaled with "scale" instead of content values
  //use range [m,M] if limits are specified

  if(!h2) return;
  TH2 *h2e = (TH2F*)h2->Clone(Form("%s_E", h2->GetName()));
  h2e->SetContour(10);
  if(M>m) h2e->GetZaxis()->SetRangeUser(m, M);
  if(title) h2e->GetZaxis()->SetTitle(title);
  
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
Bool_t AliTRDresolution::MakeProjectionCluster()
{
// Analyse cluster
  const Int_t kNcontours(9);
  const Int_t kNstat(300);
  Int_t cidx = kCluster;
  if(fProj && fProj->At(cidx)) return kTRUE;
  if(!fContainer){
    AliError("Missing data container.");
    return kFALSE;
  }
  THnSparse *H(NULL);
  if(!(H = (THnSparse*)fContainer->At(cidx))){
    AliError(Form("Missing/Wrong data @ %d.", cidx));
    return kFALSE;
  }
  Int_t ndim(H->GetNdimensions());
  Int_t coord[kNdim]; memset(coord, 0, sizeof(Int_t) * kNdim); Double_t v = 0.;
  TAxis *aa[kNdim], *as(NULL), *apt(NULL); memset(aa, 0, sizeof(TAxis*) * kNdim);
  for(Int_t id(0); id<ndim; id++) aa[id] = H->GetAxis(id);
  if(ndim > kPt) apt = H->GetAxis(kPt);
  if(ndim > kSpeciesChgRC) as  = H->GetAxis(kSpeciesChgRC);
  // build list of projections
  const Int_t nsel(12), npsel(5);
  // define rebinning strategy
  const Int_t nEtaPhi(4); Int_t rebinEtaPhiX[nEtaPhi] = {1, 2, 5, 1}, rebinEtaPhiY[nEtaPhi] = {2, 1, 1, 5};
  AliTRDresolutionProjection hp[fgkNproj[cidx]], *php[nsel][npsel]; memset(php, 0, nsel*npsel*sizeof(AliTRDresolutionProjection*));
  Int_t ih(0), isel(-1), np[nsel]; memset(np, 0, nsel*sizeof(Int_t));
  for(Int_t ily(0); ily<AliTRDgeometry::kNlayer; ily++){
    isel++; // new selection
    hp[ih].Build(Form("HClY%d", ily), Form("Clusters :: r-#phi residuals ly%d", ily), kEta, kPhi, kYrez, aa);
    hp[ih].SetRebinStrategy(nEtaPhi, rebinEtaPhiX, rebinEtaPhiY);
      php[isel][np[isel]++] = &hp[ih++];
    hp[ih].Build(Form("HClYn%d", ily), Form("Clusters[-]:: r-#phi residuals ly%d", ily), kEta, kPhi, kYrez, aa);
    hp[ih].SetRebinStrategy(nEtaPhi, rebinEtaPhiX, rebinEtaPhiY);
      php[isel][np[isel]++] = &hp[ih++];
    hp[ih].Build(Form("HClQn%d", ily), Form("Clusters[-]:: Charge distribution ly%d", ily), kEta, kPhi, kSpeciesChgRC, aa);
    hp[ih].SetRebinStrategy(nEtaPhi, rebinEtaPhiX, rebinEtaPhiY);
    hp[ih].SetShowRange(20., 40.);
      php[isel][np[isel]++] = &hp[ih++];
    hp[ih].Build(Form("HClYXTCn%d", ily), Form("Clusters[-]:: r-#phi(x,TC) residuals ly%d", ily), kPrez, kZrez, kYrez, aa);
//    hp[ih].SetRebinStrategy(nEtaPhi, rebinEtaPhiX, rebinEtaPhiY);
      php[isel][np[isel]++] = &hp[ih++];
    hp[ih].Build(Form("HClYXPh%d", ily), Form("Clusters :: r-#phi(x,#Phi) residuals ly%d", ily), kPrez, kPt, kYrez, aa);
//    hp[ih].SetRebinStrategy(nEtaPhi, rebinEtaPhiX, rebinEtaPhiY);
      php[isel][np[isel]++] = &hp[ih++];
    isel++; // new selection
      php[isel][np[isel]++] = &hp[ih-5]; // relink HClY
    hp[ih].Build(Form("HClYp%d", ily), Form("Clusters[+]:: r-#phi residuals ly%d", ily), kEta, kPhi, kYrez, aa);
    hp[ih].SetRebinStrategy(nEtaPhi, rebinEtaPhiX, rebinEtaPhiY);
      php[isel][np[isel]++] = &hp[ih++];
    hp[ih].Build(Form("HClQp%d", ily), Form("Clusters[+]:: Charge distribution ly%d", ily), kEta, kPhi, kSpeciesChgRC, aa);
    hp[ih].SetShowRange(20., 40.);
    hp[ih].SetRebinStrategy(nEtaPhi, rebinEtaPhiX, rebinEtaPhiY);
      php[isel][np[isel]++] = &hp[ih++];
    hp[ih].Build(Form("HClYXTCp%d", ily), Form("Clusters[+]:: r-#phi(x,TC) residuals ly%d", ily), kPrez, kZrez, kYrez, aa);
//    hp[ih].SetRebinStrategy(nEtaPhi, rebinEtaPhiX, rebinEtaPhiY);
      php[isel][np[isel]++] = &hp[ih++];
      php[isel][np[isel]++] = &hp[ih-4]; // relink HClYXPh
  }

  Int_t ly(0), ch(0), rcBin(as?as->FindBin(0.):-1), chBin(apt?apt->FindBin(0.):-1);
  for (Long64_t ib(0); ib < H->GetNbins(); ib++) {
    v = H->GetBinContent(ib, coord); if(v<1.) continue;
    ly = coord[kBC]-1;
    // RC selection
    if(rcBin>0 && coord[kSpeciesChgRC] == rcBin) continue;

    // charge selection
    ch = 0; // [-] track
    if(chBin>0 && coord[kPt] > chBin) ch = 1;  // [+] track

    isel = ly*2+ch;
    for(Int_t jh(0); jh<np[isel]; jh++) php[isel][jh]->Increment(coord, v);
  }

  if(!fProj){
    AliInfo("Building array of projections ...");
    fProj = new TObjArray(kNclasses); fProj->SetOwner(kTRUE);
  }
  TObjArray *arr(NULL);
  fProj->AddAt(arr = new TObjArray(fgkNproj[cidx]), cidx);

  TH2 *h2(NULL);
  for(; ih--; ){
    if(!hp[ih].fH) continue;
    Int_t mid(1), nstat(kNstat);
    if(strchr(hp[ih].fH->GetName(), 'Q')){ mid=2; nstat=300;}
    if(!(h2 = hp[ih].Projection2D(nstat, kNcontours, mid))) continue;
    arr->AddAt(h2, ih);
  }

  return kTRUE;
}

//________________________________________________________
Bool_t AliTRDresolution::MakeProjectionTracklet()
{
// Analyse tracklet
  const Int_t kNcontours(9);
  const Int_t kNstat(100);
  Int_t cidx = kTracklet;
  if(fProj && fProj->At(cidx)) return kTRUE;
  if(!fContainer){
    AliError("Missing data container.");
    return kFALSE;
  }
  THnSparse *H(NULL);
  if(!(H = (THnSparse*)fContainer->At(cidx))){
    AliError(Form("Missing/Wrong data @ %d.", cidx));
    return kFALSE;
  }
  Int_t ndim(H->GetNdimensions());
  Int_t coord[kNdim+1]; memset(coord, 0, sizeof(Int_t) * (kNdim+1)); Double_t v = 0.;
  TAxis *aa[kNdim+1], *as(NULL); memset(aa, 0, sizeof(TAxis*) * (kNdim+1));
  for(Int_t id(0); id<ndim; id++) aa[id] = H->GetAxis(id);
  if(ndim > kSpeciesChgRC) as = H->GetAxis(kSpeciesChgRC);
  // build list of projections
  const Int_t nsel(18), npsel(6);
  // define rebinning strategy
  const Int_t nEtaPhi(4); Int_t rebinEtaPhiX[nEtaPhi] = {1, 2, 5, 1}, rebinEtaPhiY[nEtaPhi] = {2, 1, 1, 5};
  AliTRDresolutionProjection hp[fgkNproj[cidx]], *php[nsel][npsel]; memset(php, 0, nsel*npsel*sizeof(AliTRDresolutionProjection*));
  Int_t ih(0), isel(-1), np[nsel]; memset(np, 0, nsel*sizeof(Int_t));
  for(Int_t ily(0); ily<AliTRDgeometry::kNlayer; ily++){
    isel++; // new selection
    hp[ih].Build(Form("HTrkltY%d", ily), Form("Tracklets   :: r-#phi residuals ly%d", ily), kEta, kPhi, kYrez, aa);
    hp[ih].SetRebinStrategy(nEtaPhi, rebinEtaPhiX, rebinEtaPhiY);
      php[isel][np[isel]++] = &hp[ih++];
    hp[ih].Build(Form("HTrkltYn%d", ily), Form("Tracklets[-]:: r-#phi residuals ly%d", ily), kEta, kPhi, kYrez, aa);
    hp[ih].SetRebinStrategy(nEtaPhi, rebinEtaPhiX, rebinEtaPhiY);
      php[isel][np[isel]++] = &hp[ih++];
    hp[ih].Build(Form("HTrkltPhn%d", ily), Form("Tracklets[-]:: #Delta#phi residuals ly%d", ily), kEta, kPhi, kPrez, aa);
    hp[ih].SetRebinStrategy(nEtaPhi, rebinEtaPhiX, rebinEtaPhiY);
      php[isel][np[isel]++] = &hp[ih++];
    hp[ih].Build(Form("HTrkltPn%d", ily), Form("Tracklets[-]:: Momentum distribution ly%d", ily), kEta, kPhi, kPt, aa);
    hp[ih].SetShowRange(6.,12.);
    hp[ih].SetRebinStrategy(nEtaPhi, rebinEtaPhiX, rebinEtaPhiY);
      php[isel][np[isel]++] = &hp[ih++];
    hp[ih].Build(Form("HTrkltYPn%d", ily), Form("Tracklets[-]:: r-#phi/p_{t} residuals ly%d", ily), kPt, kPhi, kYrez, aa);
      php[isel][np[isel]++] = &hp[ih++];
    hp[ih].Build(Form("HTrkltQn%d", ily), Form("Tracklets[-]:: dQdl ly%d", ily), kEta, kPhi, kNdim, aa);
    hp[ih].SetShowRange(700.,1100.);
    hp[ih].SetRebinStrategy(nEtaPhi, rebinEtaPhiX, rebinEtaPhiY);
      php[isel][np[isel]++] = &hp[ih++];
    isel++; // new selection
    php[isel][np[isel]++] = &hp[ih-6]; // relink first histo
    hp[ih].Build(Form("HTrkltYp%d", ily), Form("Tracklets[+]:: r-#phi residuals ly%d", ily), kEta, kPhi, kYrez, aa);
    hp[ih].SetRebinStrategy(nEtaPhi, rebinEtaPhiX, rebinEtaPhiY);
      php[isel][np[isel]++] = &hp[ih++];
    hp[ih].Build(Form("HTrkltPhp%d", ily), Form("Tracklets[+]:: #Delta#phi residuals ly%d", ily), kEta, kPhi, kPrez, aa);
    hp[ih].SetRebinStrategy(nEtaPhi, rebinEtaPhiX, rebinEtaPhiY);
      php[isel][np[isel]++] = &hp[ih++];
    hp[ih].Build(Form("HTrkltPp%d", ily), Form("Tracklets[+]:: Momentum distribution ly%d", ily), kEta, kPhi, kPt, aa);
    hp[ih].SetShowRange(6.,12.);
    hp[ih].SetRebinStrategy(nEtaPhi, rebinEtaPhiX, rebinEtaPhiY);
      php[isel][np[isel]++] = &hp[ih++];
    hp[ih].Build(Form("HTrkltYPp%d", ily), Form("Tracklets[+]:: r-#phi/p_{t} residuals ly%d", ily), kPt, kPhi, kYrez, aa);
      php[isel][np[isel]++] = &hp[ih++];
    hp[ih].Build(Form("HTrkltQp%d", ily), Form("Tracklets[+]:: dQdl ly%d", ily), kEta, kPhi, kNdim, aa);
    hp[ih].SetShowRange(700.,1100.);
    hp[ih].SetRebinStrategy(nEtaPhi, rebinEtaPhiX, rebinEtaPhiY);
      php[isel][np[isel]++] = &hp[ih++];
    isel++; // new selection
    hp[ih].Build(Form("HTrkltZ%d", ily), Form("Tracklets[RC]:: z residuals ly%d", ily), kEta, kPhi, kZrez, aa);
    hp[ih].SetRebinStrategy(nEtaPhi, rebinEtaPhiX, rebinEtaPhiY);
      php[isel][np[isel]++] = &hp[ih++];
  }

  Int_t ly(0), ch(0), rcBin(as?as->FindBin(0.):-1);
  for (Long64_t ib(0); ib < H->GetNbins(); ib++) {
    v = H->GetBinContent(ib, coord);
    if(v<1.) continue;
    ly = coord[kBC]-1; // layer selection
    // charge selection
    ch = 0; // [-] track
    if(rcBin>0){ // debug mode in which species are also saved
      if(coord[kSpeciesChgRC] > rcBin) ch = 1;  // [+] track
      else if(coord[kSpeciesChgRC] == rcBin) ch = 2;  // [RC] track
    }
    isel = ly*3+ch;
    for(Int_t jh(0); jh<np[isel]; jh++) php[isel][jh]->Increment(coord, v);
  }

  if(!fProj){
    AliInfo("Building array of projections ...");
    fProj = new TObjArray(kNclasses); fProj->SetOwner(kTRUE);
  }
  TObjArray *arr(NULL);
  fProj->AddAt(arr = new TObjArray(fgkNproj[cidx]), cidx);

  TH2 *h2(NULL);
  for(; ih--; ){
    if(!hp[ih].fH) continue;
    Int_t mid(0), nstat(kNstat);
    if(strchr(hp[ih].fH->GetName(), 'Q')){ mid=2; nstat=300;}
    if(!(h2 = hp[ih].Projection2D(nstat, kNcontours, mid))) continue;
    arr->AddAt(h2, ih);
  }
  return kTRUE;
}

//________________________________________________________
Bool_t AliTRDresolution::MakeProjectionTrackIn()
{
// Analyse track in

  const Int_t kNcontours(9);
  const Int_t kNstat(30);

  Int_t cidx = kTrackIn;
  if(fProj && fProj->At(cidx)) return kTRUE;
  if(!fContainer){
    AliError("Missing data container.");
    return kFALSE;
  }
  THnSparse *H(NULL);
  if(!(H = (THnSparse*)fContainer->At(cidx))){
    AliError(Form("Missing/Wrong data @ %d.", Int_t(cidx)));
    return kFALSE;
  }

  Int_t coord[kNdim]; memset(coord, 0, sizeof(Int_t) * kNdim); Double_t v = 0.;
  Int_t ndim(H->GetNdimensions());
  TAxis *aa[kNdim+1], *as(NULL); memset(aa, 0, sizeof(TAxis*) * (kNdim+1));
  for(Int_t id(0); id<ndim; id++) aa[id] = H->GetAxis(id);
  if(ndim > kSpeciesChgRC) as = H->GetAxis(kSpeciesChgRC);
  // build list of projections
  const Int_t nsel(3), npsel(4);
  // define rebinning strategy
  const Int_t nEtaPhi(4); Int_t rebinEtaPhiX[nEtaPhi] = {1, 2, 5, 1}, rebinEtaPhiY[nEtaPhi] = {2, 1, 1, 5};
  const Int_t nPtPhi(2); Int_t rebinPtPhiX[nEtaPhi] = {1, 1}, rebinPtPhiY[nEtaPhi] = {2, 5};
  AliTRDresolutionProjection hp[fgkNproj[cidx]], *php[nsel][npsel]; memset(php, 0, nsel*npsel*sizeof(AliTRDresolutionProjection*));
  Int_t ih(0), isel(-1), np[nsel]; memset(np, 0, nsel*sizeof(Int_t));
  // define list of projections
  isel++;  // negative tracks
  hp[ih].Build("HTrkInY", "TrackIn :: r-#phi residuals", kEta, kPhi, kYrez, aa);
  hp[ih].SetRebinStrategy(nEtaPhi, rebinEtaPhiX, rebinEtaPhiY);
    php[isel][np[isel]++] = &hp[ih++];
  hp[ih].Build("HTrkInYn", "TrackIn[-]:: r-#phi residuals", kEta, kPhi, kYrez, aa);
  hp[ih].SetRebinStrategy(nEtaPhi, rebinEtaPhiX, rebinEtaPhiY);
    php[isel][np[isel]++] = &hp[ih++];
  hp[ih].Build("HTrkInPhn", "TrackIn[-]:: #Delta#phi residuals", kEta, kPhi, kPrez, aa);
  hp[ih].SetRebinStrategy(nEtaPhi, rebinEtaPhiX, rebinEtaPhiY);
    php[isel][np[isel]++] = &hp[ih++];
  hp[ih].Build("HTrkInYPn", "TrackIn[-]:: r-#phi/p_{t} residuals", kPt, kPhi, kYrez, aa);
  hp[ih].SetRebinStrategy(nPtPhi, rebinPtPhiX, rebinPtPhiY);
    php[isel][np[isel]++] = &hp[ih++];
  isel++; // positive tracks
  php[isel][np[isel]++] = &hp[ih-4]; // relink first histo
  hp[ih].Build("HTrkInYp", "TrackIn[+]:: r-#phi residuals", kEta, kPhi, kYrez, aa);
  hp[ih].SetRebinStrategy(nEtaPhi, rebinEtaPhiX, rebinEtaPhiY);
    php[isel][np[isel]++] = &hp[ih++];
  hp[ih].Build("HTrkInPhp", "TrackIn[+]:: #Delta#phi residuals", kEta, kPhi, kPrez, aa);
  hp[ih].SetRebinStrategy(nEtaPhi, rebinEtaPhiX, rebinEtaPhiY);
    php[isel][np[isel]++] = &hp[ih++];
  hp[ih].Build("HTrkInYPp", "TrackIn[+]:: r-#phi/p_{t} residuals", kPt, kPhi, kYrez, aa);
  hp[ih].SetRebinStrategy(nPtPhi, rebinPtPhiX, rebinPtPhiY);
    php[isel][np[isel]++] = &hp[ih++];
  isel++; // RC tracks
  hp[ih].Build("HTrkInZ", "TrackIn[RC]:: z residuals", kEta, kPhi, kZrez, aa);
  hp[ih].SetRebinStrategy(nEtaPhi, rebinEtaPhiX, rebinEtaPhiY);
    php[isel][np[isel]++] = &hp[ih++];

  // fill projections
  Int_t ch(0), rcBin(as?as->FindBin(0.):-1);
  for (Long64_t ib(0); ib < H->GetNbins(); ib++) {
    v = H->GetBinContent(ib, coord);
    if(v<1.) continue;
    if(coord[kBC]>1) continue; // bunch cross cut
    // charge selection
    ch = 0; // [-] track
    if(rcBin>0){ // debug mode in which species are also saved
      if(coord[kSpeciesChgRC] > rcBin) ch = 1;  // [+] track
      else if(coord[kSpeciesChgRC] == rcBin) ch = 2;  // [RC] track
    }
    for(Int_t jh(0); jh<np[ch]; jh++) php[ch][jh]->Increment(coord, v);
  }
  if(!fProj){
    AliInfo("Building array of projections ...");
    fProj = new TObjArray(kNclasses); fProj->SetOwner(kTRUE);
  }
  TObjArray *arr(NULL);
  fProj->AddAt(arr = new TObjArray(fgkNproj[cidx]), cidx);

  TH2 *h2(NULL);
  for(; ih--; ){
    if(!hp[ih].fH) continue;
    if(!(h2 = hp[ih].Projection2D(kNstat, kNcontours/*, mid*/))) continue;
    arr->AddAt(h2, ih);
  }
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

  // DEFINE MODELS
  // simple gauss
  TF1 fg("fGauss", "gaus", -.5, .5);  
  // Landau for charge resolution
  TF1 fch("fClCh", "landau", 0., 1000.);  
  // Landau for e+- pt resolution
  TF1 fpt("fPt", "landau", -0.1, 0.2);  

  //PROCESS EXPERIMENTAL DISTRIBUTIONS
  // Clusters residuals
  if(!MakeProjectionCluster()) return kFALSE;
  fNRefFigures = 3;
  // Tracklet residual/pulls
  if(!MakeProjectionTracklet()) return kFALSE;
  fNRefFigures = 7;
  // TRDin residual/pulls
  if(!MakeProjectionTrackIn()) return kFALSE;
  fNRefFigures = 11;
  // TRDout residual/pulls
//  if(!MakeProjectionTrackOut()) return kFALSE;
  fNRefFigures = 15;

  if(!HasMCdata()) return kTRUE;


  //PROCESS MC RESIDUAL DISTRIBUTIONS

  // CLUSTER Y RESOLUTION/PULLS
//  if(!MakeProjectionClusterMC()) return kFALSE;
  fNRefFigures = 17;

  // TRACKLET RESOLUTION/PULLS
//  if(!MakeProjectionTrackletMC()) return kFALSE;
  fNRefFigures = 21;

  // TRACK RESOLUTION/PULLS
/*  Process3Darray(kMCtrack, 0, &fg, 1.e4);   // y
  Process3DlinkedArray(kMCtrack, 1, &fg);   // y PULL
  Process3Darray(kMCtrack, 2, &fg, 1.e4);   // z
  Process3Darray(kMCtrack, 3, &fg);         // z PULL
  Process2Darray(kMCtrack, 4, &fg, 1.e3);   // phi
  Process2Darray(kMCtrack, 5, &fg);         // snp PULL
  Process2Darray(kMCtrack, 6, &fg, 1.e3);   // theta
  Process2Darray(kMCtrack, 7, &fg);         // tgl PULL
  Process3Darray(kMCtrack, 8, &fg, 1.e2);   // pt resolution
  Process3Darray(kMCtrack, 9, &fg);         // 1/pt pulls
  Process3Darray(kMCtrack, 10, &fg, 1.e2);  // p resolution*/
//  if(!MakeProjectionTrackMC(kMCtrack)) return kFALSE;
  fNRefFigures+=16;

  // TRACK TRDin RESOLUTION/PULLS
//  if(!MakeProjectionTrackMC(kMCtrackIn)) return kFALSE;
  fNRefFigures+=8;

  // TRACK TRDout RESOLUTION/PULLS
//  if(!MakeProjectionTrackMC(kMCtrackOut)) return kFALSE;
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
  // cluster to tracklet residuals/pulls
  snprintf(hn, nhn, "h%s", fgPerformanceName[kCluster]);
  if(!(H = (THnSparseI*)gROOT->FindObject(hn))){
    const Char_t *clTitle[kNdim] = {"layer", fgkTitle[kPhi], fgkTitle[kEta], fgkTitle[kYrez], "#Deltax [cm]", "Q</Q", "Q/angle", "#Phi [deg]"};
    const Int_t clNbins[kNdim]   = {AliTRDgeometry::kNlayer, fgkNbins[kPhi], fgkNbins[kEta], fgkNbins[kYrez], 33, 10, 30, 15};
    const Double_t clMin[kNdim]  = {-0.5, fgkMin[kPhi], fgkMin[kEta], fgkMin[kYrez]/10., 0., 0.1, -2., -45},
                   clMax[kNdim]  = {AliTRDgeometry::kNlayer-0.5, fgkMax[kPhi], fgkMax[kEta], fgkMax[kYrez]/10., 3.3, 2.1, 118., 45};
    st = "cluster spatial&charge resolution;";
    // define minimum info to be saved in non debug mode
    Int_t ndim=DebugLevel()>=1?kNdim:4;
    for(Int_t idim(0); idim<ndim; idim++){ st += clTitle[idim]; st+=";";}
    H = new THnSparseI(hn, st.Data(), ndim, clNbins, clMin, clMax);
  } else H->Reset();
  fContainer->AddAt(H, kCluster);
  //++++++++++++++++++++++
  // tracklet to TRD track
  snprintf(hn, nhn, "h%s", fgPerformanceName[kTracklet]);
  if(!(H = (THnSparseI*)gROOT->FindObject(hn))){
    Char_t *trTitle[kNdim+1]; memcpy(trTitle, fgkTitle, kNdim*sizeof(Char_t*));
    Int_t trNbins[kNdim+1]; memcpy(trNbins, fgkNbins, kNdim*sizeof(Int_t));
    Double_t trMin[kNdim+1]; memcpy(trMin, fgkMin, kNdim*sizeof(Double_t));
    Double_t trMax[kNdim+1]; memcpy(trMax, fgkMax, kNdim*sizeof(Double_t));
    // set specific fields
    trTitle[kBC]=StrDup("layer"); trNbins[kBC] = AliTRDgeometry::kNlayer; trMin[kBC] = -0.5; trMax[kBC] = AliTRDgeometry::kNlayer-0.5;
    trTitle[kNdim]=StrDup("dq/dl [a.u.]"); trNbins[kNdim] = 30; trMin[kNdim] = 100.; trMax[kNdim] = 3100;

    st = "tracklet spatial&charge resolution;";
    // define minimum info to be saved in non debug mode
    Int_t ndim=DebugLevel()>=1?(kNdim+1):4;
    for(Int_t idim(0); idim<ndim; idim++){ st += trTitle[idim]; st+=";";}
    H = new THnSparseI(hn, st.Data(), ndim, trNbins, trMin, trMax);
  } else H->Reset();
  fContainer->AddAt(H, kTracklet);
  //++++++++++++++++++++++
  // tracklet to TRDin
  snprintf(hn, nhn, "h%s", fgPerformanceName[kTrackIn]);
  if(!(H = (THnSparseI*)gROOT->FindObject(hn))){
    st = "r-#phi/z/angular residuals @ TRD entry;";
    // define minimum info to be saved in non debug mode
    Int_t ndim=DebugLevel()>=1?kNdim:7;
    for(Int_t idim(0); idim<ndim; idim++){ st += fgkTitle[idim]; st+=";";}
    H = new THnSparseI(hn, st.Data(), ndim, fgkNbins, fgkMin, fgkMax);
  } else H->Reset();
  fContainer->AddAt(H, kTrackIn);
  // tracklet to TRDout
//  fContainer->AddAt(BuildMonitorContainerTracklet("TrkOUT"), kTrackOut);


  // Resolution histos
  if(!HasMCdata()) return fContainer;

  // cluster resolution 
  fContainer->AddAt(BuildMonitorContainerCluster("MCcl"),  kMCcluster);

  // tracklet resolution
  fContainer->AddAt(BuildMonitorContainerTracklet("MCtracklet"), kMCtracklet);

  // track resolution
  TObjArray *arr(NULL);
  fContainer->AddAt(arr = new TObjArray(AliTRDgeometry::kNlayer), kMCtrack);
  arr->SetName("MCtrk");
  for(Int_t il(0); il<AliTRDgeometry::kNlayer; il++) arr->AddAt(BuildMonitorContainerTrack(Form("MCtrk_Ly%d", il)), il);

  // TRDin TRACK RESOLUTION
  fContainer->AddAt(H, kMCtrackIn);

  // TRDout TRACK RESOLUTION
  fContainer->AddAt(BuildMonitorContainerTrack("MCtrkOUT"), kMCtrackOut);

  return fContainer;
}

//________________________________________________________
Bool_t AliTRDresolution::Process(TH2* const h2, TGraphErrors **g, Int_t stat)
{
// Robust function to process sigma/mean for 2D plot dy(x)
// For each x bin a gauss fit is performed on the y projection and the range
// with the minimum chi2/ndf is choosen

  if(!h2) {
    if(AliLog::GetDebugLevel("PWG1", "AliTRDresolution")>0) printf("D-AliTRDresolution::Process() : NULL pointer input container.\n");
    return kFALSE;
  }
  if(!Int_t(h2->GetEntries())){
    if(AliLog::GetDebugLevel("PWG1", "AliTRDresolution")>0) printf("D-AliTRDresolution::Process() : Empty h[%s - %s].\n", h2->GetName(), h2->GetTitle());
    return kFALSE;
  }
  if(!g || !g[0]|| !g[1]) {
    if(AliLog::GetDebugLevel("PWG1", "AliTRDresolution")>0) printf("D-AliTRDresolution::Process() : NULL pointer output container.\n");
    return kFALSE;
  }
  // prepare
  TAxis *ax(h2->GetXaxis()), *ay(h2->GetYaxis());
  Float_t ymin(ay->GetXmin()), ymax(ay->GetXmax()), dy(ay->GetBinWidth(1)), y0(0.), y1(0.);
  TF1 f("f", "gaus", ymin, ymax);
  Int_t n(0);
  if((n=g[0]->GetN())) for(;n--;) g[0]->RemovePoint(n);
  if((n=g[1]->GetN())) for(;n--;) g[1]->RemovePoint(n);
  TH1D *h(NULL);
  if((h=(TH1D*)gROOT->FindObject("py"))) delete h;
  Double_t x(0.), y(0.), ex(0.), ey(0.), sy(0.), esy(0.);
  

  // do actual loop
  Float_t chi2OverNdf(0.);
  for(Int_t ix = 1, np=0; ix<=ax->GetNbins(); ix++){
    x = ax->GetBinCenter(ix); ex= ax->GetBinWidth(ix)*0.288; // w/sqrt(12)
    ymin = ay->GetXmin(); ymax = ay->GetXmax();

    h = h2->ProjectionY("py", ix, ix);
    if((n=(Int_t)h->GetEntries())<stat){
      if(AliLog::GetDebugLevel("PWG1", "AliTRDresolution")>1) printf("I-AliTRDresolution::Process() : Low statistics @ x[%f] stat[%d]=%d [%d].\n", x, ix, n, stat);
      continue;
    }
    // looking for a first order mean value
    f.SetParameter(1, 0.); f.SetParameter(2, 3.e-2);
    h->Fit(&f, "QNW");
    chi2OverNdf = f.GetChisquare()/f.GetNDF();
    printf("x[%f] range[%f %f] chi2[%f] ndf[%d] chi2/ndf[%f]\n", x, ymin, ymax, f.GetChisquare(),f.GetNDF(),chi2OverNdf);
    y = f.GetParameter(1); y0 = y-4*dy; y1 = y+4*dy;
    ey  = f.GetParError(1);
    sy = f.GetParameter(2); esy = f.GetParError(2);
//     // looking for the best chi2/ndf value
//     while(ymin<y0 && ymax>y1){
//       f.SetParameter(1, y);
//       f.SetParameter(2, sy);
//       h->Fit(&f, "QNW", "", y0, y1);
//       printf("   range[%f %f] chi2[%f] ndf[%d] chi2/ndf[%f]\n", y0, y1, f.GetChisquare(),f.GetNDF(),f.GetChisquare()/f.GetNDF());
//       if(f.GetChisquare()/f.GetNDF() < Chi2OverNdf){
//         chi2OverNdf = f.GetChisquare()/f.GetNDF();
//         y  = f.GetParameter(1); ey  = f.GetParError(1);
//         sy = f.GetParameter(2); esy = f.GetParError(2);
//         printf("    set y[%f] sy[%f] chi2/ndf[%f]\n", y, sy, chi2OverNdf);
//       }
//       y0-=dy; y1+=dy;
//     }
    g[0]->SetPoint(np, x, y);
    g[0]->SetPointError(np, ex, ey);
    g[1]->SetPoint(np, x, sy);
    g[1]->SetPointError(np, ex, esy);
    np++;
  }
  return kTRUE;
}


//________________________________________________________
Bool_t AliTRDresolution::Process(TH2 * const h2, TF1 *f, Float_t k, TGraphErrors **g)
{
  //
  // Do the processing
  //

  Char_t pn[10]; snprintf(pn, 10, "p%03d", fIdxPlot);
  Int_t n = 0;
  if((n=g[0]->GetN())) for(;n--;) g[0]->RemovePoint(n);
  if((n=g[1]->GetN())) for(;n--;) g[1]->RemovePoint(n);
  if(Int_t(h2->GetEntries())){ 
    AliDebug(4, Form("%s: g[%s %s]", pn, g[0]->GetName(), g[0]->GetTitle()));
  } else {
    AliDebug(2, Form("%s: g[%s %s]: Missing entries.", pn, g[0]->GetName(), g[0]->GetTitle()));
    fIdxPlot++;
    return kTRUE;
  }

  const Int_t kINTEGRAL=1;
  for(Int_t ibin = 0; ibin < Int_t(h2->GetNbinsX()/kINTEGRAL); ibin++){
    Int_t abin(ibin*kINTEGRAL+1),bbin(abin+kINTEGRAL-1),mbin(abin+Int_t(kINTEGRAL/2));
    Double_t x = h2->GetXaxis()->GetBinCenter(mbin);
    TH1D *h = h2->ProjectionY(pn, abin, bbin);
    if((n=(Int_t)h->GetEntries())<150){ 
      AliDebug(4, Form("  x[%f] range[%d %d] stat[%d] low statistics !", x, abin, bbin, n));
      continue;
    }
    h->Fit(f, "QN");
    Int_t ip = g[0]->GetN();
    AliDebug(4, Form("  x_%d[%f] range[%d %d] stat[%d] M[%f] Sgm[%f]", ip, x, abin, bbin, n, f->GetParameter(1), f->GetParameter(2)));
    g[0]->SetPoint(ip, x, k*f->GetParameter(1));
    g[0]->SetPointError(ip, 0., k*f->GetParError(1));
    g[1]->SetPoint(ip, x, k*f->GetParameter(2));
    g[1]->SetPointError(ip, 0., k*f->GetParError(2));
/*  
    g[0]->SetPoint(ip, x, k*h->GetMean());
    g[0]->SetPointError(ip, 0., k*h->GetMeanError());
    g[1]->SetPoint(ip, x, k*h->GetRMS());
    g[1]->SetPointError(ip, 0., k*h->GetRMSError());*/
  }
  fIdxPlot++;
  return kTRUE;
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
    if(AliLog::GetDebugLevel("PWG1", "AliTRDresolution")>1) printf("D-AliTRDresolution::FitTrack: Not enough clusters to fit a track [%d].\n", np);
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
  if(AliLog::GetDebugLevel("PWG1", "AliTRDresolution")>3) printf("D-AliTRDresolution::FitTrack: x0[%f] y0[%f] z0[%f] dydx[%f] dzdx[%f].\n", x0, y0, z0, dydx, dzdx);
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
    if(AliLog::GetDebugLevel("PWG1", "AliTRDresolution")>1) printf("D-AliTRDresolution::FitTracklet: Not enough clusters to fit a tracklet [%d].\n", nly);
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
//   if(AliLog::GetDebugLevel("PWG1", "AliTRDresolution")>=2){
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
  :fH(NULL)
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
  fH = new TH3I(n, Form("%s;%s;%s;%s", t, ax->GetTitle(), ay->GetTitle(), az->GetTitle()),
    ax->GetNbins(), ax->GetXmin(), ax->GetXmax(),
    ay->GetNbins(), ay->GetXmin(), ay->GetXmax(),
    az->GetNbins(), az->GetXmin(), az->GetXmax());
  fAx[0] = ix; fAx[1] = iy; fAx[2] = iz;
  fRange[0] = az->GetXmin()/3.; fRange[1] = az->GetXmax()/3.;
}

//________________________________________________________
void AliTRDresolution::AliTRDresolutionProjection::Increment(Int_t bin[], Double_t v)
{
// increment bin with value "v" pointed by general coord in "bin"
  if(!fH) return;
  fH->AddBinContent(
        fH->GetBin(bin[fAx[0]],bin[fAx[1]],bin[fAx[2]]), v);
}

//________________________________________________________
TH2* AliTRDresolution::AliTRDresolutionProjection::Projection2D(const Int_t nstat, const Int_t ncol, const Int_t mid)
{
// build the 2D projection and adjust binning

  const Char_t *title[] = {"Mean", "#mu", "MPV"};
  if(!fH) return NULL;
  TAxis *ax(fH->GetXaxis()), *ay(fH->GetYaxis()), *az(fH->GetZaxis());
  TH2 *h2s = (TH2*)fH->Project3D("yx");
  Int_t irebin(0), dxBin(1), dyBin(1);
  while(irebin<fNrebin && (AliTRDresolution::GetMeanStat(h2s, .5, ">")<nstat)){
    h2s->Rebin2D(fRebinX[irebin], fRebinY[irebin]);
    dxBin*=fRebinX[irebin];dyBin*=fRebinY[irebin];
    irebin++;
  }
  Int_t nx(h2s->GetNbinsX()), ny(h2s->GetNbinsY());
  if(h2s) delete h2s;

  // start projection
  TH1 *h(NULL);
  Float_t dz=(fRange[1]-fRange[1])/ncol;
  TString titlez(az->GetTitle()); TObjArray *tokenTitle(titlez.Tokenize(" "));
  TH2 *h2 = new TH2F(Form("%s_2D", fH->GetName()),
            Form("%s;%s;%s;%s(%s) %s", fH->GetTitle(), ax->GetTitle(), ay->GetTitle(), title[mid], (*tokenTitle)[0]->GetName(), tokenTitle->GetEntriesFast()>1?(*tokenTitle)[1]->GetName():""),
            nx, ax->GetXmin(), ax->GetXmax(), ny, ay->GetXmin(), ay->GetXmax());
  h2->SetContour(ncol);
  h2->GetZaxis()->CenterTitle();
  h2->GetZaxis()->SetRangeUser(fRange[0], fRange[1]);
  //printf("%s[%s] nx[%d] ny[%d]\n", h2->GetName(), h2->GetTitle(), nx, ny);
  for(Int_t iy(0); iy<ny; iy++){
    for(Int_t ix(0); ix<nx; ix++){
      h = fH->ProjectionZ(Form("%s_z", h2->GetName()), ix*dxBin+1, (ix+1)*dxBin+1, iy*dyBin+1, (iy+1)*dyBin+1);
      Int_t ne((Int_t)h->Integral());
      if(ne<nstat/2){
        h2->SetBinContent(ix+1, iy+1, -999);
        h2->SetBinError(ix+1, iy+1, 1.);
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
  return h2;
}

void AliTRDresolution::AliTRDresolutionProjection::SetRebinStrategy(Int_t n, Int_t rebx[], Int_t reby[])
{
// define rebinning strategy for this projection
  fNrebin = n;
  fRebinX = new Int_t[n]; memcpy(fRebinX, rebx, n*sizeof(Int_t));
  fRebinY = new Int_t[n]; memcpy(fRebinY, reby, n*sizeof(Int_t));
}


