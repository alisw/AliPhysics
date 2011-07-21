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
//__________________________________________________________________________
AliTRDcheckTRK::AliTRDcheckTRK()
  : AliTRDrecoTask()
{
// Default constructor
  SetNameTitle("TRDtrackerSys", "TRD Tracker Systematic");
  Float_t pt0(0.2);
  for(Int_t j(0); j<=kNpt; j++){
    pt0+=(TMath::Exp(j*j*.002)-1.);
    fPtBins[j]=pt0;
  }
  memset(fProj, 0, 10*sizeof(TH1*));
}

//__________________________________________________________________________
AliTRDcheckTRK::AliTRDcheckTRK(char* name)
  : AliTRDrecoTask(name, "TRD Tracker Systematic")
{
// User constructor
  InitFunctorList();
  Float_t pt0(0.2);
  for(Int_t j(0); j<=kNpt; j++){
    pt0+=(TMath::Exp(j*j*.002)-1.);
    fPtBins[j]=pt0;
  }
  memset(fProj, 0, 10*sizeof(TH1*));
}

//__________________________________________________________________________
AliTRDcheckTRK::~AliTRDcheckTRK()
{
// Destructor
}

//__________________________________________________________________________
Int_t AliTRDcheckTRK::GetPtBin(Float_t pt)
{
// Find pt bin according to local pt segmentation
  Int_t ipt(0);
  while(ipt<kNpt){
    if(pt<fPtBins[ipt]) break;
    ipt++;
  }
  return TMath::Max(0,ipt);
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
Bool_t AliTRDcheckTRK::GetRefFigure(Int_t ifig)
{
  switch(ifig){
  case 0:
    if(!MakeProjectionEtaPhi()) return kFALSE;
    break;
  }
  return kTRUE;
}

//__________________________________________________________________________
TObjArray* AliTRDcheckTRK::Histos()
{
  //
  // Create QA histogram tree
  //

  if(fContainer) return fContainer;

  THnSparseI *h(NULL);
  fContainer = new TObjArray(kNclasses);
  fContainer->SetOwner(kTRUE);

  const Char_t *ccn[kNclasses] = {"Entry", "Propag"};
  const Char_t *ctt[kNclasses] = {"r-#phi/z/angular residuals @ entry", "r-#phi/z/angular residuals each ly"};

  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  Int_t   bins[kNdim+1] = {Int_t(kNbunchCross)/*bc*/, 180/*phi*/, 50/*eta*/, Int_t(kNcharge)*AliPID::kSPECIES+1/*chg*species*/, kNpt/*pt*/, 50/*dy*/, 50/*dz*/, 40/*dphi*/, AliTRDgeometry::kNlayer};
  Double_t min[kNdim+1] = {-0.5, -TMath::Pi(), -1., -AliPID::kSPECIES-0.5, -0.5, -1.5, -2.5, -10., -0.5},
           max[kNdim+1] = {Int_t(kNbunchCross)-0.5, TMath::Pi(), 1., AliPID::kSPECIES+0.5, kNpt-0.5, 1.5, 2.5, 10., AliTRDgeometry::kNlayer-0.5};
  Char_t hn[100], ht[700];
  snprintf(hn, 100, "h%s", ccn[kEntry]);
  if(!(h = (THnSparseI*)gROOT->FindObject(hn))){
    snprintf(ht, 700, "%s;bunch cross;#phi [rad];#eta;chg*spec*rc;bin_p_{t};#Deltay [cm];#Deltaz [cm];#Delta#phi [deg];", ctt[kEntry]);
    h = new THnSparseI(hn, ht, kNdim, bins, min, max);
  } else h->Reset();
  fContainer->AddAt(h, kEntry);

  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  min[5] = -0.15; max[5] = 0.15;
  snprintf(hn, 100, "h%s", ccn[kPropagation]);
  if(!(h = (THnSparseI*)gROOT->FindObject(hn))){
    snprintf(ht, 700, "%s;bunch cross;#phi [rad];#eta;chg*spec*rc;bin_p_{t};#Deltay [cm];#Deltaz [cm];#Delta#phi [deg];layer;", ctt[kPropagation]);
    h = new THnSparseI(hn, ht, kNdim+1, bins, min, max);
  } else h->Reset();
  fContainer->AddAt(h, kPropagation);

  return fContainer;
}

//__________________________________________________________________________
TH1* AliTRDcheckTRK::PlotEntry(const AliTRDtrackV1 *track)
{
// comment needed
  if(track) fkTrack = track;
  if(!fkTrack){
    AliDebug(4, "No Track defined.");
    return NULL;
  }
  // check container
  THnSparseI *h=(THnSparseI*)fContainer->At(kEntry);
  if(!h){
    AliError(Form("Missing container @ %d", Int_t(kEntry)));
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

  Double_t xyz[3];
  if(!tin->GetXYZ(xyz)){
    AliDebug(1, "Failed getting global track postion");
    xyz[0]=x; xyz[1]=0.; xyz[2]=0.;
  }
  Float_t mass(fkTrack->GetMass()),
          eta(tin->Eta()),
          phi(TMath::ATan2(xyz[1], xyz[0]));
  Int_t charge(fkTrack->Charge()),
        species(GetSpeciesByMass(mass)),
        bc(fkESD->GetTOFbc());
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
  val[kBC]          = (bc>=kNbunchCross)?(kNbunchCross-1):bc;
  val[kPhi]         = phi;
  val[kEta]         = eta;
  val[kSpeciesChgRC]= fTracklet->IsRowCross()?0:charge*(species+1);
  val[kPt]          = GetPtBin(tin->Pt());
  val[kYrez]        = dy;
  val[kZrez]        = dz;
  val[kPrez]        = dphi*TMath::RadToDeg();
  h->Fill(val);
  return NULL;
}

//__________________________________________________________________________
TH1* AliTRDcheckTRK::PlotPropagation(const AliTRDtrackV1 *track)
{
// comment needed

  if(track) fkTrack = track;
  if(!fkTrack){
    AliDebug(4, "No Track defined.");
    return NULL;
  }
  // check container
  THnSparseI *h=(THnSparseI*)fContainer->At(kPropagation);
  if(!h){
    AliError(Form("Missing container @ %d", Int_t(kPropagation)));
    return NULL;
  }
  // check input track status
  AliExternalTrackParam *tin(NULL);
  if(!(tin = fkTrack->GetTrackIn())){
    AliError("Track did not entered TRD fiducial volume.");
    return NULL;
  }
  
  Float_t mass(fkTrack->GetMass());
  Int_t charge(fkTrack->Charge()),
        species(GetSpeciesByMass(mass)),
        bc(fkESD->GetTOFbc()+1);
  TVectorD dX(AliTRDgeometry::kNlayer),
           dY(AliTRDgeometry::kNlayer),
           dZ(AliTRDgeometry::kNlayer),
           dPhi(AliTRDgeometry::kNlayer),
           vPt(AliTRDgeometry::kNlayer),
           vPhi(AliTRDgeometry::kNlayer),
           vEta(AliTRDgeometry::kNlayer),
           budget(AliTRDgeometry::kNlayer),
           cCOV(AliTRDgeometry::kNlayer*15);
  const Int_t nopt(1000); Char_t opt[nopt];

  // propagation
  //snprintf(opt, nopt, "");

  AliTRDseedV1 *tr(NULL);
  Double_t val[kNdim+1];
  if(PropagateKalman(fkTrack, tin, &dX, &dY, &dZ, &dPhi, &vPt, &vPhi, &vEta, &budget, &cCOV, opt)){
    val[kBC]=(bc>=kNbunchCross)?(kNbunchCross-1):bc;
    for(Int_t ily(0); ily<AliTRDgeometry::kNlayer; ily++){
      if(dX[ily]<0.) continue;
      if(!(tr = fkTrack->GetTracklet(ily))) continue;
      val[kPhi]         = vPhi[ily];
      val[kEta]         = vEta[ily];
      val[kSpeciesChgRC]= tr->IsRowCross()?0:charge*(species+1);
      val[kPt]          = GetPtBin(vPt[ily]);
      val[kYrez]        = dY[ily];
      val[kZrez]        = dZ[ily];
      val[kPrez]        = dPhi[ily]*TMath::RadToDeg();
      val[kNdim]        = ily;
      h->Fill(val);
    }
  }

  return NULL;
}


//___________________________________________________
Bool_t AliTRDcheckTRK::PropagateKalman(const AliTRDtrackV1 *t, AliExternalTrackParam *ref, TVectorD *dx, TVectorD *dy, TVectorD *dz, TVectorD *dphi, TVectorD *pt, TVectorD *eta, TVectorD *phi, TVectorD *budget, TVectorD *cov, Option_t */*opt*/)
{
// Propagate Kalman from the first TRD tracklet to
// last one and save residuals in the y, z and pt.
// The track is intialized from the reference "ref".
//
// This is to calibrate the dEdx and MS corrections

  Int_t ntracklets(t->GetNumberOfTracklets());
  if(!ntracklets) return kFALSE;
  if(ref->Pt()<1.e-3) return kFALSE;

  AliTRDseedV1 *tr(NULL);
  for(Int_t itr(AliTRDgeometry::kNlayer); itr--;){
    (*dx)[itr] = -1.; (*dy)[itr] = 100.; (*dz)[itr] = 100.; (*dphi)[itr] = 100.;
  }

  // Initialize TRD track to the reference
  AliTRDtrackV1 tt;
  tt.Set(ref->GetX(), ref->GetAlpha(), ref->GetParameter(), ref->GetCovariance());
  tt.SetMass(t->GetMass());

  Double_t x0(ref->GetX()), xyz[3] = {0., 0., 0.};
  for(Int_t ily(0); ily<AliTRDgeometry::kNlayer; ily++){
    if(!(tr = t->GetTracklet(ily))) continue;
    if(fgClRecalibrate){
      AliTRDtransform trans(tr->GetDetector());
      AliTRDcluster *c(NULL);
      tr->ResetClusterIter(kFALSE);
      while((c = tr->PrevCluster())) trans.Recalibrate(c, kFALSE);
      tr->FitRobust();
    }
    if(!AliTRDtrackerV1::PropagateToX(tt, tr->GetX0(), fgKalmanStep)) continue;
    if(fgKalmanUpdate){
      Double_t x(tr->GetX0()),
               p[2] = { tr->GetYfit(0), tr->GetZfit(0)},
               covTrklt[3];
      tr->GetCovAt(x, covTrklt);
      if(!((AliExternalTrackParam&)tt).Update(p, covTrklt)) continue;
    }
    if(!tt.GetXYZ(xyz)) continue;
    (*phi)[ily] = TMath::ATan2(xyz[1], xyz[0]);
    (*eta)[ily] = tt.Eta();
    (*dx)[ily]  = tt.GetX() - x0;
    (*pt)[ily]  = tt.Pt();
    if(budget) (*budget)[ily] = tt.GetBudget(0);
    if(cov){
      const Double_t *cc(tt.GetCovariance());
      for(Int_t ic(0), jp(ily*15); ic<15; ic++, jp++) (*cov)[jp]=cc[ic];
    }
    const Double_t *parR(tt.GetParameter());
    Double_t dyt(parR[0] - tr->GetYfit(0)), dzt(parR[1] - tr->GetZfit(0)),
             dydx(tr->GetYfit(1)),
             tilt(tr->GetTilt());
    // correct for tilt rotation
    (*dy)[ily]  = dyt - dzt*tilt;
    (*dz)[ily]  = dzt + dyt*tilt;
    dydx       += tilt*parR[3];
    (*dphi)[ily]= TMath::ASin(parR[2])-TMath::ATan(dydx);
  }
  return kTRUE;
}


//__________________________________________________________________________
Bool_t AliTRDcheckTRK::MakeProjectionEtaPhi()
{
// Get dy residual projection as a function of eta and phi

  if(fProj[0] && fProj[1]) return kTRUE;
  if(!fContainer || fContainer->GetEntries() != kNclasses){
    AliError("Missing/Wrong data container.");
    return kFALSE;
  }
  THnSparse *H(NULL);
  if(!(H = (THnSparse*)fContainer->At(kEntry))){
    AliError(Form("Missing/Wrong data @ %d.", Int_t(kEntry)));
    return kFALSE;
  }

  Int_t  coord[kNdim]; memset(coord, 0, sizeof(Int_t) * kNdim); Double_t v = 0.;
  TAxis //*abc(H->GetAxis(kBC)),
        *aphi(H->GetAxis(kPhi)),
        *aeta(H->GetAxis(kEta)),
        //*as(H->GetAxis(kSpeciesChgRC)),
        //*apt(H->GetAxis(kPt)),
        *ay(H->GetAxis(kYrez)),
        *az(H->GetAxis(kZrez));
        //*ap(H->GetAxis(kPrez));
  Int_t neta(aeta->GetNbins()), nphi(aphi->GetNbins());
  TH3I *h3[3];
  h3[0] = new TH3I("h30", Form("r-#phi residuals for neg tracks;%s;%s;%s", aeta->GetTitle(), aphi->GetTitle(), ay->GetTitle()),
            neta, aeta->GetXmin(), aeta->GetXmax(),
            nphi, aphi->GetXmin(), aphi->GetXmax(),
            ay->GetNbins(), ay->GetXmin(), ay->GetXmax());
  h3[1] = new TH3I("h31", Form("z residuals for row cross;%s;%s;%s", aeta->GetTitle(), aphi->GetTitle(), az->GetTitle()),
            neta, aeta->GetXmin(), aeta->GetXmax(),
            nphi, aphi->GetXmin(), aphi->GetXmax(),
            az->GetNbins(), az->GetXmin(), az->GetXmax());
  h3[2] = new TH3I("h32", Form("r-#phi residuals for pos tracks;%s;%s;%s", aeta->GetTitle(), aphi->GetTitle(), ay->GetTitle()),
            neta, aeta->GetXmin(), aeta->GetXmax(),
            nphi, aphi->GetXmin(), aphi->GetXmax(),
            ay->GetNbins(), ay->GetXmin(), ay->GetXmax());
  for (Long64_t ib(0); ib < H->GetNbins(); ib++) {
    v = H->GetBinContent(ib, coord);
    if(coord[kBC]>1) continue; // bunch cross cut
    // species selection
    if(coord[kSpeciesChgRC]<6){
      h3[0]->AddBinContent(
        h3[0]->GetBin(coord[kEta], coord[kPhi], coord[kYrez]), v);
    } else if(coord[kSpeciesChgRC]==6) {
      h3[1]->AddBinContent(
        h3[1]->GetBin(coord[kEta], coord[kPhi], coord[kZrez]), v);
    } else if(coord[kSpeciesChgRC]>6) {
      h3[2]->AddBinContent(
        h3[2]->GetBin(coord[kEta], coord[kPhi], coord[kYrez]), v);
    }
  }
  printf("start Z projection\n");
  TH2F *h2[3];
  for(Int_t iproj(0); iproj<3; iproj++){
    h2[iproj] = new TH2F(Form("h2%d", iproj),
            Form("%s;%s;%s;%s", h3[iproj]->GetTitle(), aeta->GetTitle(), aphi->GetTitle(), h3[iproj]->GetZaxis()->GetTitle()),
            neta, aeta->GetXmin(), aeta->GetXmax(),
            nphi, aphi->GetXmin(), aphi->GetXmax());
    for(Int_t iphi(1); iphi<=nphi; iphi++){
      for(Int_t ieta(1); ieta<=neta; ieta++){
        TH1 *h = h3[iproj]->ProjectionZ(Form("hy%d", iproj), ieta, ieta, iphi, iphi);
        h2[iproj]->SetBinContent(ieta, iphi, h->GetEntries()>20?h->GetMean():-999.);
      }
    }
  }

  return kTRUE;
}
