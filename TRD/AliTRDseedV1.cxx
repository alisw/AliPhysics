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

/* $Id$ */

////////////////////////////////////////////////////////////////////////////
////
//  The TRD offline tracklet
//
// The running horse of the TRD reconstruction. The following tasks are preformed:
//   1. Clusters attachment to tracks based on prior information stored at tracklet level (see AttachClusters)
//   2. Clusters position recalculation based on track information (see GetClusterXY and Fit)
//   3. Cluster error parametrization recalculation (see Fit)
//   4. Linear track approximation (Fit)
//   5. Optimal position (including z estimate for pad row cross tracklets) and covariance matrix of the track fit inside one TRD chamber (Fit)
//   6. Tilt pad correction and systematic effects (GetCovAt)
//   7. dEdx calculation (CookdEdx)
//   8. PID probabilities estimation (CookPID)
//
//  Authors:                                                              //
//    Alex Bercuci <A.Bercuci@gsi.de>                                     //
//    Markus Fasel <M.Fasel@gsi.de>                                       //
//                                                                        //
////////////////////////////////////////////////////////////////////////////

#include "TMath.h"
#include <TTreeStream.h>

#include "AliLog.h"
#include "AliMathBase.h"
#include "AliCDBManager.h"
#include "AliTracker.h"

#include "AliTRDpadPlane.h"
#include "AliTRDcluster.h"
#include "AliTRDseedV1.h"
#include "AliTRDtrackV1.h"
#include "AliTRDcalibDB.h"
#include "AliTRDchamberTimeBin.h"
#include "AliTRDtrackingChamber.h"
#include "AliTRDtrackerV1.h"
#include "AliTRDrecoParam.h"
#include "AliTRDCommonParam.h"

#include "Cal/AliTRDCalPID.h"
#include "Cal/AliTRDCalROC.h"
#include "Cal/AliTRDCalDet.h"

ClassImp(AliTRDseedV1)

//____________________________________________________________________
AliTRDseedV1::AliTRDseedV1(Int_t det) 
  :AliTRDtrackletBase()
  ,fkReconstructor(NULL)
  ,fClusterIter(NULL)
  ,fExB(0.)
  ,fVD(0.)
  ,fT0(0.)
  ,fS2PRF(0.)
  ,fDiffL(0.)
  ,fDiffT(0.)
  ,fClusterIdx(0)
  ,fErrorMsg(0)
  ,fN(0)
  ,fDet(det)
  ,fPt(0.)
  ,fdX(0.)
  ,fX0(0.)
  ,fX(0.)
  ,fY(0.)
  ,fZ(0.)
  ,fS2Y(0.)
  ,fS2Z(0.)
  ,fChi2(0.)
{
  //
  // Constructor
  //
  memset(fIndexes,0xFF,kNclusters*sizeof(fIndexes[0]));
  memset(fClusters, 0, kNclusters*sizeof(AliTRDcluster*));
  memset(fPad, 0, 4*sizeof(Float_t));
  fYref[0] = 0.; fYref[1] = 0.; 
  fZref[0] = 0.; fZref[1] = 0.; 
  fYfit[0] = 0.; fYfit[1] = 0.; 
  fZfit[0] = 0.; fZfit[1] = 0.; 
  memset(fdEdx, 0, kNslices*sizeof(Float_t)); 
  for(int ispec=0; ispec<AliPID::kSPECIES; ispec++) fProb[ispec]  = -1.;
  fLabels[0]=-1; fLabels[1]=-1; // most freq MC labels
  fLabels[2]=0;  // number of different labels for tracklet
  memset(fRefCov, 0, 7*sizeof(Double_t));
  // stand alone curvature
  fC[0] = 0.; fC[1] = 0.; 
  // covariance matrix [diagonal]
  // default sy = 200um and sz = 2.3 cm 
  fCov[0] = 4.e-4; fCov[1] = 0.; fCov[2] = 5.3; 
  SetStandAlone(kFALSE);
}

//____________________________________________________________________
AliTRDseedV1::AliTRDseedV1(const AliTRDseedV1 &ref)
  :AliTRDtrackletBase((AliTRDtrackletBase&)ref)
  ,fkReconstructor(NULL)
  ,fClusterIter(NULL)
  ,fExB(0.)
  ,fVD(0.)
  ,fT0(0.)
  ,fS2PRF(0.)
  ,fDiffL(0.)
  ,fDiffT(0.)
  ,fClusterIdx(0)
  ,fErrorMsg(0)
  ,fN(0)
  ,fDet(-1)
  ,fPt(0.)
  ,fdX(0.)
  ,fX0(0.)
  ,fX(0.)
  ,fY(0.)
  ,fZ(0.)
  ,fS2Y(0.)
  ,fS2Z(0.)
  ,fChi2(0.)
{
  //
  // Copy Constructor performing a deep copy
  //
  if(this != &ref){
    ref.Copy(*this);
  }
  SetBit(kOwner, kFALSE);
  SetStandAlone(ref.IsStandAlone());
}


//____________________________________________________________________
AliTRDseedV1& AliTRDseedV1::operator=(const AliTRDseedV1 &ref)
{
  //
  // Assignment Operator using the copy function
  //

  if(this != &ref){
    ref.Copy(*this);
  }
  SetBit(kOwner, kFALSE);

  return *this;
}

//____________________________________________________________________
AliTRDseedV1::~AliTRDseedV1()
{
  //
  // Destructor. The RecoParam object belongs to the underlying tracker.
  //

  //printf("I-AliTRDseedV1::~AliTRDseedV1() : Owner[%s]\n", IsOwner()?"YES":"NO");

  if(IsOwner()) {
    for(int itb=0; itb<kNclusters; itb++){
      if(!fClusters[itb]) continue; 
      //AliInfo(Form("deleting c %p @ %d", fClusters[itb], itb));
      delete fClusters[itb];
      fClusters[itb] = NULL;
    }
  }
}

//____________________________________________________________________
void AliTRDseedV1::Copy(TObject &ref) const
{
  //
  // Copy function
  //

  //AliInfo("");
  AliTRDseedV1 &target = (AliTRDseedV1 &)ref; 

  target.fkReconstructor = fkReconstructor;
  target.fClusterIter   = NULL;
  target.fExB           = fExB;
  target.fVD            = fVD;
  target.fT0            = fT0;
  target.fS2PRF         = fS2PRF;
  target.fDiffL         = fDiffL;
  target.fDiffT         = fDiffT;
  target.fClusterIdx    = 0;
  target.fErrorMsg      = fErrorMsg;
  target.fN             = fN;
  target.fDet           = fDet;
  target.fPt            = fPt;
  target.fdX            = fdX;
  target.fX0            = fX0;
  target.fX             = fX;
  target.fY             = fY;
  target.fZ             = fZ;
  target.fS2Y           = fS2Y;
  target.fS2Z           = fS2Z;
  target.fChi2          = fChi2;
  
  memcpy(target.fIndexes, fIndexes, kNclusters*sizeof(Int_t));
  memcpy(target.fClusters, fClusters, kNclusters*sizeof(AliTRDcluster*));
  memcpy(target.fPad, fPad, 4*sizeof(Float_t));
  target.fYref[0] = fYref[0]; target.fYref[1] = fYref[1]; 
  target.fZref[0] = fZref[0]; target.fZref[1] = fZref[1]; 
  target.fYfit[0] = fYfit[0]; target.fYfit[1] = fYfit[1]; 
  target.fZfit[0] = fZfit[0]; target.fZfit[1] = fZfit[1]; 
  memcpy(target.fdEdx, fdEdx, kNslices*sizeof(Float_t)); 
  memcpy(target.fProb, fProb, AliPID::kSPECIES*sizeof(Float_t)); 
  memcpy(target.fLabels, fLabels, 3*sizeof(Int_t)); 
  memcpy(target.fRefCov, fRefCov, 7*sizeof(Double_t)); 
  target.fC[0] = fC[0]; target.fC[1] = fC[1];
  memcpy(target.fCov, fCov, 3*sizeof(Double_t)); 
  
  TObject::Copy(ref);
}


//____________________________________________________________
Bool_t AliTRDseedV1::Init(AliTRDtrackV1 *track)
{
// Initialize this tracklet using the track information
//
// Parameters:
//   track - the TRD track used to initialize the tracklet
// 
// Detailed description
// The function sets the starting point and direction of the
// tracklet according to the information from the TRD track.
// 
// Caution
// The TRD track has to be propagated to the beginning of the
// chamber where the tracklet will be constructed
//

  Double_t y, z; 
  if(!track->GetProlongation(fX0, y, z)) return kFALSE;
  Update(track);
  return kTRUE;
}


//_____________________________________________________________________________
void AliTRDseedV1::Reset(Option_t *opt)
{
//
// Reset seed. If option opt="c" is given only cluster arrays are cleared.
//
  for(Int_t ic=kNclusters; ic--;) fIndexes[ic] = -1;
  memset(fClusters, 0, kNclusters*sizeof(AliTRDcluster*));
  fN=0; SetBit(kRowCross, kFALSE);
  if(strcmp(opt, "c")==0) return;

  fExB=0.;fVD=0.;fT0=0.;fS2PRF=0.;
  fDiffL=0.;fDiffT=0.;
  fClusterIdx=0;
  fErrorMsg = 0;
  fDet=-1;
  fPt=0.;
  fdX=0.;fX0=0.; fX=0.; fY=0.; fZ=0.;
  fS2Y=0.; fS2Z=0.;
  fC[0]=0.; fC[1]=0.; 
  fChi2 = 0.;

  memset(fPad, 0, 4*sizeof(Float_t));
  fYref[0] = 0.; fYref[1] = 0.; 
  fZref[0] = 0.; fZref[1] = 0.; 
  fYfit[0] = 0.; fYfit[1] = 0.; 
  fZfit[0] = 0.; fZfit[1] = 0.; 
  memset(fdEdx, 0, kNslices*sizeof(Float_t)); 
  for(int ispec=0; ispec<AliPID::kSPECIES; ispec++) fProb[ispec]  = -1.;
  fLabels[0]=-1; fLabels[1]=-1; // most freq MC labels
  fLabels[2]=0;  // number of different labels for tracklet
  memset(fRefCov, 0, 7*sizeof(Double_t));
  // covariance matrix [diagonal]
  // default sy = 200um and sz = 2.3 cm 
  fCov[0] = 4.e-4; fCov[1] = 0.; fCov[2] = 5.3; 
}

//____________________________________________________________________
void AliTRDseedV1::Update(const AliTRDtrackV1 *trk)
{ 
  // update tracklet reference position from the TRD track

  Double_t fSnp = trk->GetSnp();
  Double_t fTgl = trk->GetTgl();
  fPt = trk->Pt();
  Double_t norm =1./TMath::Sqrt((1.-fSnp)*(1.+fSnp)); 
  fYref[1] = fSnp*norm;
  fZref[1] = fTgl*norm;
  SetCovRef(trk->GetCovariance());

  Double_t dx = trk->GetX() - fX0;
  fYref[0] = trk->GetY() - dx*fYref[1];
  fZref[0] = trk->GetZ() - dx*fZref[1];
}

//_____________________________________________________________________________
void AliTRDseedV1::UpdateUsed()
{
  //
  // Calculate number of used clusers in the tracklet
  //

  Int_t nused = 0, nshared = 0;
  for (Int_t i = kNclusters; i--; ) {
    if (!fClusters[i]) continue;
    if(fClusters[i]->IsUsed()){ 
      nused++;
    } else if(fClusters[i]->IsShared()){
      if(IsStandAlone()) nused++;
      else nshared++;
    }
  }
  SetNUsed(nused);
  SetNShared(nshared);
}

//_____________________________________________________________________________
void AliTRDseedV1::UseClusters()
{
  //
  // Use clusters
  //
  // In stand alone mode:
  // Clusters which are marked as used or shared from another track are
  // removed from the tracklet
  //
  // In barrel mode:
  // - Clusters which are used by another track become shared
  // - Clusters which are attached to a kink track become shared
  //
  AliTRDcluster **c = &fClusters[0];
  for (Int_t ic=kNclusters; ic--; c++) {
    if(!(*c)) continue;
    if(IsStandAlone()){
      if((*c)->IsShared() || (*c)->IsUsed()){ 
        if((*c)->IsShared()) SetNShared(GetNShared()-1);
        else SetNUsed(GetNUsed()-1);
        (*c) = NULL;
        fIndexes[ic] = -1;
        SetN(GetN()-1);
        continue;
      }
    } else {
      if((*c)->IsUsed() || IsKink()){
        (*c)->SetShared();
        continue;
      }
    }
    (*c)->Use();
  }
}



//____________________________________________________________________
void AliTRDseedV1::CookdEdx(Int_t nslices)
{
// Calculates average dE/dx for all slices and store them in the internal array fdEdx. 
//
// Parameters:
//  nslices : number of slices for which dE/dx should be calculated
// Output:
//  store results in the internal array fdEdx. This can be accessed with the method
//  AliTRDseedV1::GetdEdx()
//
// Detailed description
// Calculates average dE/dx for all slices. Depending on the PID methode 
// the number of slices can be 3 (LQ) or 8(NN). 
// The calculation of dQ/dl are done using the tracklet fit results (see AliTRDseedV1::GetdQdl(Int_t))
//
// The following effects are included in the calculation:
// 1. calibration values for t0 and vdrift (using x coordinate to calculate slice)
// 2. cluster sharing (optional see AliTRDrecoParam::SetClusterSharing())
// 3. cluster size
//

  memset(fdEdx, 0, kNslices*sizeof(Float_t));
  const Double_t kDriftLength = (.5 * AliTRDgeometry::AmThick() + AliTRDgeometry::DrThick());

  AliTRDcluster *c(NULL);
  for(int ic=0; ic<AliTRDtrackerV1::GetNTimeBins(); ic++){
    if(!(c = fClusters[ic]) && !(c = fClusters[ic+kNtb])) continue;
    Float_t dx = TMath::Abs(fX0 - c->GetX());
    
    // Filter clusters for dE/dx calculation
    
    // 1.consider calibration effects for slice determination
    Int_t slice;
    if(dx<kDriftLength){ // TODO should be replaced by c->IsInChamber() 
      slice = Int_t(dx * nslices / kDriftLength);
    } else slice = c->GetX() < fX0 ? nslices-1 : 0;


    // 2. take sharing into account
    Float_t w = /*c->IsShared() ? .5 :*/ 1.;
    
    // 3. take into account large clusters TODO
    //w *= c->GetNPads() > 3 ? .8 : 1.;
    
    //CHECK !!!
    fdEdx[slice]   += w * GetdQdl(ic); //fdQdl[ic];
  } // End of loop over clusters
}

//_____________________________________________________________________________
void AliTRDseedV1::CookLabels()
{
  //
  // Cook 2 labels for seed
  //

  Int_t labels[200];
  Int_t out[200];
  Int_t nlab = 0;
  for (Int_t i = 0; i < kNclusters; i++) {
    if (!fClusters[i]) continue;
    for (Int_t ilab = 0; ilab < 3; ilab++) {
      if (fClusters[i]->GetLabel(ilab) >= 0) {
        labels[nlab] = fClusters[i]->GetLabel(ilab);
        nlab++;
      }
    }
  }

  fLabels[2] = AliMathBase::Freq(nlab,labels,out,kTRUE);
  fLabels[0] = out[0];
  if ((fLabels[2]  > 1) && (out[3] > 1)) fLabels[1] = out[2];
}

//____________________________________________________________
Float_t AliTRDseedV1::GetAnodeWireOffset(Float_t zt)
{
  Float_t d = fPad[3] - zt;
  if(d<0.){
    AliError(Form("Fail AnodeWireOffset calculation z0[%+7.2f] zt[%+7.2f] d[%+7.2f].", fPad[3], zt, d));
    return 0.125;
  } 
  d -= ((Int_t)(2 * d)) / 2.0;
  if(d > 0.25) d = 0.5 - d;
  return d;
}


//____________________________________________________________________
Float_t AliTRDseedV1::GetdQdl(Int_t ic, Float_t *dl) const
{
// Using the linear approximation of the track inside one TRD chamber (TRD tracklet) 
// the charge per unit length can be written as:
// BEGIN_LATEX
// #frac{dq}{dl} = #frac{q_{c}}{dx * #sqrt{1 + #(){#frac{dy}{dx}}^{2}_{fit} + #(){#frac{dz}{dx}}^{2}_{ref}}}
// END_LATEX
// where qc is the total charge collected in the current time bin and dx is the length 
// of the time bin. 
// The following correction are applied :
//   - charge : pad row cross corrections
//              [diffusion and TRF assymetry] TODO
//   - dx     : anisochronity, track inclination - see Fit and AliTRDcluster::GetXloc() 
//              and AliTRDcluster::GetYloc() for the effects taken into account
// 
//Begin_Html
//<img src="TRD/trackletDQDT.gif">
//End_Html
// In the picture the energy loss measured on the tracklet as a function of drift time [left] and respectively 
// drift length [right] for different particle species is displayed.
// Author : Alex Bercuci <A.Bercuci@gsi.de>
//
  Float_t dq = 0.;
  // check whether both clusters are inside the chamber
  Bool_t hasClusterInChamber = kFALSE;
  if(fClusters[ic] && fClusters[ic]->IsInChamber()){
    hasClusterInChamber = kTRUE;
    dq += TMath::Abs(fClusters[ic]->GetQ());
  }else if(fClusters[ic+kNtb] && fClusters[ic+kNtb]->IsInChamber()){
    hasClusterInChamber = kTRUE;
    dq += TMath::Abs(fClusters[ic+kNtb]->GetQ());
  }
  if(!hasClusterInChamber) return 0.;
  if(dq<1.e-3) return 0.;

  Double_t dx = fdX;
  if(ic-1>=0 && ic+1<kNtb){
    Float_t x2(0.), x1(0.);
    // try to estimate upper radial position (find the cluster which is inside the chamber)
    if(fClusters[ic-1] && fClusters[ic-1]->IsInChamber()) x2 = fClusters[ic-1]->GetX(); 
    else if(fClusters[ic-1+kNtb] && fClusters[ic-1+kNtb]->IsInChamber()) x2 = fClusters[ic-1+kNtb]->GetX(); 
    else if(fClusters[ic] && fClusters[ic]->IsInChamber()) x2 = fClusters[ic]->GetX()+fdX;
    else x2 = fClusters[ic+kNtb]->GetX()+fdX;
    // try to estimate lower radial position (find the cluster which is inside the chamber)
    if(fClusters[ic+1] && fClusters[ic+1]->IsInChamber()) x1 = fClusters[ic+1]->GetX();
    else if(fClusters[ic+1+kNtb] && fClusters[ic+1+kNtb]->IsInChamber()) x1 = fClusters[ic+1+kNtb]->GetX();
    else if(fClusters[ic] && fClusters[ic]->IsInChamber()) x1 = fClusters[ic]->GetX()-fdX;
    else x1 = fClusters[ic+kNtb]->GetX()-fdX;

    dx = .5*(x2 - x1);
  }
  dx *= TMath::Sqrt(1. + fYfit[1]*fYfit[1] + fZref[1]*fZref[1]);
  if(dl) (*dl) = dx;
  if(dx>1.e-9) return dq/dx;
  else return 0.;
}

//____________________________________________________________
Float_t AliTRDseedV1::GetMomentum(Float_t *err) const
{ 
// Returns momentum of the track after update with the current tracklet as:
// BEGIN_LATEX
// p=#frac{1}{1/p_{t}} #sqrt{1+tgl^{2}}
// END_LATEX
// and optionally the momentum error (if err is not null). 
// The estimated variance of the momentum is given by:
// BEGIN_LATEX
// #sigma_{p}^{2} = (#frac{dp}{dp_{t}})^{2} #sigma_{p_{t}}^{2}+(#frac{dp}{dtgl})^{2} #sigma_{tgl}^{2}+2#frac{dp}{dp_{t}}#frac{dp}{dtgl} cov(tgl,1/p_{t})
// END_LATEX
// which can be simplified to
// BEGIN_LATEX
// #sigma_{p}^{2} = p^{2}p_{t}^{4}tgl^{2}#sigma_{tgl}^{2}-2p^{2}p_{t}^{3}tgl cov(tgl,1/p_{t})+p^{2}p_{t}^{2}#sigma_{1/p_{t}}^{2}
// END_LATEX
//

  Double_t p = fPt*TMath::Sqrt(1.+fZref[1]*fZref[1]);
  Double_t p2 = p*p;
  Double_t tgl2 = fZref[1]*fZref[1];
  Double_t pt2 = fPt*fPt;
  if(err){
    Double_t s2 = 
      p2*tgl2*pt2*pt2*fRefCov[4]
     -2.*p2*fZref[1]*fPt*pt2*fRefCov[5]
     +p2*pt2*fRefCov[6];
    (*err) = TMath::Sqrt(s2);
  }
  return p;
}

//____________________________________________________________________
Float_t AliTRDseedV1::GetOccupancyTB() const
{
// Returns procentage of TB occupied by clusters

  Int_t n(0);
  AliTRDcluster *c(NULL);
  for(int ic=0; ic<AliTRDtrackerV1::GetNTimeBins(); ic++){
    if(!(c = fClusters[ic]) && !(c = fClusters[ic+kNtb])) continue;
    n++;
  }

  return Float_t(n)/AliTRDtrackerV1::GetNTimeBins();
}

//____________________________________________________________________
Float_t* AliTRDseedV1::GetProbability(Bool_t force)
{	
  if(!force) return &fProb[0];
  if(!CookPID()) return NULL;
  return &fProb[0];
}

//____________________________________________________________
Bool_t AliTRDseedV1::CookPID()
{
// Fill probability array for tracklet from the DB.
//
// Parameters
//
// Output
//   returns pointer to the probability array and NULL if missing DB access 
//
// Retrieve PID probabilities for e+-, mu+-, K+-, pi+- and p+- from the DB according to tracklet information:
// - estimated momentum at tracklet reference point 
// - dE/dx measurements
// - tracklet length
// - TRD layer
// According to the steering settings specified in the reconstruction one of the following methods are used
// - Neural Network [default] - option "nn"  
// - 2D Likelihood - option "!nn"  

  AliTRDcalibDB *calibration = AliTRDcalibDB::Instance();
  if (!calibration) {
    AliError("No access to calibration data");
    return kFALSE;
  }

  if (!fkReconstructor) {
    AliError("Reconstructor not set.");
    return kFALSE;
  }

  // Retrieve the CDB container class with the parametric detector response
  const AliTRDCalPID *pd = calibration->GetPIDObject(fkReconstructor->GetPIDMethod());
  if (!pd) {
    AliError("No access to AliTRDCalPID object");
    return kFALSE;
  }

  // calculate tracklet length TO DO
  Float_t length = (AliTRDgeometry::AmThick() + AliTRDgeometry::DrThick())/ TMath::Sqrt((1.0 - GetSnp()*GetSnp()) / (1.0 + GetTgl()*GetTgl()));
  
  //calculate dE/dx
  CookdEdx(AliTRDCalPID::kNSlicesNN);
  AliDebug(4, Form("p=%6.4f[GeV/c] dEdx{%7.2f %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f} l=%4.2f[cm]", GetMomentum(), fdEdx[0], fdEdx[1], fdEdx[2], fdEdx[3], fdEdx[4], fdEdx[5], fdEdx[6], fdEdx[7], length));

  // Sets the a priori probabilities
  Bool_t kPIDNN(fkReconstructor->GetPIDMethod()==AliTRDpidUtil::kNN);
  for(int ispec=0; ispec<AliPID::kSPECIES; ispec++)
    fProb[ispec] = pd->GetProbability(ispec, GetMomentum(), &fdEdx[0], length, kPIDNN?GetPlane():fkReconstructor->GetRecoParam()->GetPIDLQslices());
  
  return kTRUE;
}

//____________________________________________________________________
Float_t AliTRDseedV1::GetQuality(Bool_t kZcorr) const
{
  //
  // Returns a quality measurement of the current seed
  //

  Float_t zcorr = kZcorr ? GetTilt() * (fZfit[0] - fZref[0]) : 0.;
  return 
      .5 * TMath::Abs(18.0 - GetN())
    + 10.* TMath::Abs(fYfit[1] - fYref[1])
    + 5. * TMath::Abs(fYfit[0] - fYref[0] + zcorr)
    + 2. * TMath::Abs(fZfit[0] - fZref[0]) / GetPadLength();
}

//____________________________________________________________________
void AliTRDseedV1::GetCovAt(Double_t x, Double_t *cov) const
{
// Computes covariance in the y-z plane at radial point x (in tracking coordinates) 
// and returns the results in the preallocated array cov[3] as :
//   cov[0] = Var(y)
//   cov[1] = Cov(yz)
//   cov[2] = Var(z)
//
// Details
//
// For the linear transformation
// BEGIN_LATEX
// Y = T_{x} X^{T}
// END_LATEX
//   The error propagation has the general form
// BEGIN_LATEX
// C_{Y} = T_{x} C_{X} T_{x}^{T} 
// END_LATEX
//  We apply this formula 2 times. First to calculate the covariance of the tracklet 
// at point x we consider: 
// BEGIN_LATEX
// T_{x} = (1 x); X=(y0 dy/dx); C_{X}=#(){#splitline{Var(y0) Cov(y0, dy/dx)}{Cov(y0, dy/dx) Var(dy/dx)}} 
// END_LATEX
// and secondly to take into account the tilt angle
// BEGIN_LATEX
// T_{#alpha} = #(){#splitline{cos(#alpha) __ sin(#alpha)}{-sin(#alpha) __ cos(#alpha)}}; X=(y z); C_{X}=#(){#splitline{Var(y)    0}{0   Var(z)}} 
// END_LATEX
//
// using simple trigonometrics one can write for this last case
// BEGIN_LATEX
// C_{Y}=#frac{1}{1+tg^{2}#alpha} #(){#splitline{(#sigma_{y}^{2}+tg^{2}#alpha#sigma_{z}^{2}) __ tg#alpha(#sigma_{z}^{2}-#sigma_{y}^{2})}{tg#alpha(#sigma_{z}^{2}-#sigma_{y}^{2}) __ (#sigma_{z}^{2}+tg^{2}#alpha#sigma_{y}^{2})}} 
// END_LATEX
// which can be aproximated for small alphas (2 deg) with
// BEGIN_LATEX
// C_{Y}=#(){#splitline{#sigma_{y}^{2} __ (#sigma_{z}^{2}-#sigma_{y}^{2})tg#alpha}{((#sigma_{z}^{2}-#sigma_{y}^{2})tg#alpha __ #sigma_{z}^{2}}} 
// END_LATEX
//
// before applying the tilt rotation we also apply systematic uncertainties to the tracklet 
// position which can be tunned from outside via the AliTRDrecoParam::SetSysCovMatrix(). They might 
// account for extra misalignment/miscalibration uncertainties. 
//
// Author :
// Alex Bercuci <A.Bercuci@gsi.de> 
// Date : Jan 8th 2009
//


  Double_t xr     = fX0-x; 
  Double_t sy2    = fCov[0] +2.*xr*fCov[1] + xr*xr*fCov[2];
  Double_t sz2    = fS2Z;
  //GetPadLength()*GetPadLength()/12.;

  // insert systematic uncertainties
  if(fkReconstructor){
    Double_t sys[15]; memset(sys, 0, 15*sizeof(Double_t));
    fkReconstructor->GetRecoParam()->GetSysCovMatrix(sys);
    sy2 += sys[0];
    sz2 += sys[1];
  }

  // rotate covariance matrix if no RC
  if(!IsRowCross()){
    Double_t t2 = GetTilt()*GetTilt();
    Double_t correction = 1./(1. + t2);
    cov[0] = (sy2+t2*sz2)*correction;
    cov[1] = GetTilt()*(sz2 - sy2)*correction;
    cov[2] = (t2*sy2+sz2)*correction;
  } else {
    cov[0] = sy2; cov[1] = 0.; cov[2] = sy2;
  }

  AliDebug(4, Form("C(%6.1f %+6.3f %6.1f)  RC[%c]", 1.e4*TMath::Sqrt(cov[0]), cov[1], 1.e4*TMath::Sqrt(cov[2]), IsRowCross()?'y':'n'));
}

//____________________________________________________________
Int_t AliTRDseedV1::GetCovSqrt(const Double_t * const c, Double_t *d)
{
// Helper function to calculate the square root of the covariance matrix. 
// The input matrix is stored in the vector c and the result in the vector d. 
// Both arrays have to be initialized by the user with at least 3 elements. Return negative in case of failure.
// 
// For calculating the square root of the symmetric matrix c
// the following relation is used:
// BEGIN_LATEX
// C^{1/2} = VD^{1/2}V^{-1}
// END_LATEX
// with V being the matrix with the n eigenvectors as columns. 
// In case C is symmetric the followings are true:
//   - matrix D is diagonal with the diagonal given by the eigenvalues of C
//   - V = V^{-1}
//
// Author A.Bercuci <A.Bercuci@gsi.de>
// Date   Mar 19 2009

  const Double_t kZero(1.e-20);
  Double_t l[2], // eigenvalues
           v[3]; // eigenvectors
  // the secular equation and its solution :
  // (c[0]-L)(c[2]-L)-c[1]^2 = 0
  // L^2 - L*Tr(c)+DET(c) = 0
  // L12 = [Tr(c) +- sqrt(Tr(c)^2-4*DET(c))]/2
  Double_t tr = c[0]+c[2],           // trace
          det = c[0]*c[2]-c[1]*c[1]; // determinant
  if(TMath::Abs(det)<kZero) return 1;
  Double_t dd = TMath::Sqrt(tr*tr - 4*det);
  l[0] = .5*(tr + dd*(c[0]>c[2]?-1.:1.));
  l[1] = .5*(tr + dd*(c[0]>c[2]?1.:-1.));
  if(l[0]<kZero || l[1]<kZero) return 2;
  // the sym V matrix
  // | v00   v10|
  // | v10   v11|
  Double_t den = (l[0]-c[0])*(l[0]-c[0])+c[1]*c[1];
  if(den<kZero){ // almost diagonal
    v[0] = TMath::Sign(0., c[1]);
    v[1] = TMath::Sign(1., (l[0]-c[0]));
    v[2] = TMath::Sign(0., c[1]*(l[0]-c[0])*(l[1]-c[2]));
  } else {
    Double_t tmp = 1./TMath::Sqrt(den);
    v[0] = c[1]* tmp;
    v[1] = (l[0]-c[0])*tmp;
    if(TMath::Abs(l[1]-c[2])<kZero) v[2] = TMath::Sign(v[0]*(l[0]-c[0])/kZero, (l[1]-c[2]));
    else v[2] = v[0]*(l[0]-c[0])/(l[1]-c[2]);
  }
  // the VD^{1/2}V is: 
  l[0] = TMath::Sqrt(l[0]); l[1] = TMath::Sqrt(l[1]);
  d[0] = v[0]*v[0]*l[0]+v[1]*v[1]*l[1];
  d[1] = v[0]*v[1]*l[0]+v[1]*v[2]*l[1];
  d[2] = v[1]*v[1]*l[0]+v[2]*v[2]*l[1];

  return 0;
}

//____________________________________________________________
Double_t AliTRDseedV1::GetCovInv(const Double_t * const c, Double_t *d)
{
// Helper function to calculate the inverse of the covariance matrix.
// The input matrix is stored in the vector c and the result in the vector d. 
// Both arrays have to be initialized by the user with at least 3 elements
// The return value is the determinant or 0 in case of singularity.
//
// Author A.Bercuci <A.Bercuci@gsi.de>
// Date   Mar 19 2009

  Double_t det = c[0]*c[2] - c[1]*c[1];
  if(TMath::Abs(det)<1.e-20) return 0.;
  Double_t invDet = 1./det;
  d[0] = c[2]*invDet;
  d[1] =-c[1]*invDet;
  d[2] = c[0]*invDet;
  return det;
}

//____________________________________________________________________
UShort_t AliTRDseedV1::GetVolumeId() const
{
  for(Int_t ic(0);ic<kNclusters; ic++){
    if(fClusters[ic]) return fClusters[ic]->GetVolumeId();
  }
  return 0;
}


//____________________________________________________________________
void AliTRDseedV1::Calibrate()
{
// Retrieve calibration and position parameters from OCDB. 
// The following information are used
//  - detector index
//  - column and row position of first attached cluster. If no clusters are attached 
// to the tracklet a random central chamber position (c=70, r=7) will be used.
//
// The following information is cached in the tracklet
//   t0 (trigger delay)
//   drift velocity
//   PRF width
//   omega*tau = tg(a_L)
//   diffusion coefficients (longitudinal and transversal)
//
// Author :
// Alex Bercuci <A.Bercuci@gsi.de> 
// Date : Jan 8th 2009
//

  AliCDBManager *cdb = AliCDBManager::Instance();
  if(cdb->GetRun() < 0){
    AliError("OCDB manager not properly initialized");
    return;
  }

  AliTRDcalibDB *calib = AliTRDcalibDB::Instance();
  AliTRDCalROC  *vdROC = calib->GetVdriftROC(fDet),
                *t0ROC = calib->GetT0ROC(fDet);;
  const AliTRDCalDet *vdDet = calib->GetVdriftDet();
  const AliTRDCalDet *t0Det = calib->GetT0Det();

  Int_t col = 70, row = 7;
  AliTRDcluster **c = &fClusters[0];
  if(GetN()){ 
    Int_t ic = 0;
    while (ic<kNclusters && !(*c)){ic++; c++;} 
    if(*c){
      col = (*c)->GetPadCol();
      row = (*c)->GetPadRow();
    }
  }

  fT0    = (t0Det->GetValue(fDet) + t0ROC->GetValue(col,row)) / AliTRDCommonParam::Instance()->GetSamplingFrequency();
  fVD    = vdDet->GetValue(fDet) * vdROC->GetValue(col, row);
  fS2PRF = calib->GetPRFWidth(fDet, col, row); fS2PRF *= fS2PRF;
  fExB   = AliTRDCommonParam::Instance()->GetOmegaTau(fVD);
  AliTRDCommonParam::Instance()->GetDiffCoeff(fDiffL,
  fDiffT, fVD);
  AliDebug(4, Form("Calibration params for Det[%3d] Col[%3d] Row[%2d]\n  t0[%f]  vd[%f]  s2PRF[%f]  ExB[%f]  Dl[%f]  Dt[%f]", fDet, col, row, fT0, fVD, fS2PRF, fExB, fDiffL, fDiffT));


  SetBit(kCalib, kTRUE);
}

//____________________________________________________________________
void AliTRDseedV1::SetOwner()
{
  //AliInfo(Form("own [%s] fOwner[%s]", own?"YES":"NO", fOwner?"YES":"NO"));
  
  if(TestBit(kOwner)) return;
  for(int ic=0; ic<kNclusters; ic++){
    if(!fClusters[ic]) continue;
    fClusters[ic] = new AliTRDcluster(*fClusters[ic]);
  }
  SetBit(kOwner);
}

//____________________________________________________________
void AliTRDseedV1::SetPadPlane(AliTRDpadPlane *p)
{
// Shortcut method to initialize pad geometry.
  fPad[0] = p->GetLengthIPad();
  fPad[1] = p->GetWidthIPad();
  fPad[2] = TMath::Tan(TMath::DegToRad()*p->GetTiltingAngle());
  fPad[3] = p->GetRow0() + p->GetAnodeWireOffset();
}


//____________________________________________________________________
Bool_t	AliTRDseedV1::AttachClusters(AliTRDtrackingChamber *const chamber, Bool_t tilt)
{
//
// Projective algorithm to attach clusters to seeding tracklets. The following steps are performed :
// 1. Collapse x coordinate for the full detector plane
// 2. truncated mean on y (r-phi) direction
// 3. purge clusters
// 4. truncated mean on z direction
// 5. purge clusters
//
// Parameters
//  - chamber : pointer to tracking chamber container used to search the tracklet
//  - tilt    : switch for tilt correction during road building [default true]
// Output
//  - true    : if tracklet found successfully. Failure can happend because of the following:
//      -
// Detailed description
//	
// We start up by defining the track direction in the xy plane and roads. The roads are calculated based
// on tracking information (variance in the r-phi direction) and estimated variance of the standard 
// clusters (see AliTRDcluster::SetSigmaY2()) corrected for tilt (see GetCovAt()). From this the road is
// BEGIN_LATEX
// r_{y} = 3*#sqrt{12*(#sigma^{2}_{Trk}(y) + #frac{#sigma^{2}_{cl}(y) + tg^{2}(#alpha_{L})#sigma^{2}_{cl}(z)}{1+tg^{2}(#alpha_{L})})}
// r_{z} = 1.5*L_{pad}
// END_LATEX
// 
// Author : Alexandru Bercuci <A.Bercuci@gsi.de>
// Debug  : level >3

  const AliTRDrecoParam* const recoParam = fkReconstructor->GetRecoParam(); //the dynamic cast in GetRecoParam is slow, so caching the pointer to it

  if(!recoParam){
    AliError("Tracklets can not be used without a valid RecoParam.");
    return kFALSE;
  }
  // Initialize reco params for this tracklet
  // 1. first time bin in the drift region
  Int_t t0 = 14;
  Int_t kClmin = Int_t(recoParam->GetFindableClusters()*AliTRDtrackerV1::GetNTimeBins());

  Double_t sysCov[5]; recoParam->GetSysCovMatrix(sysCov);	
  Double_t s2yTrk= fRefCov[0], 
           s2yCl = 0., 
           s2zCl = GetPadLength()*GetPadLength()/12., 
           syRef = TMath::Sqrt(s2yTrk),
           t2    = GetTilt()*GetTilt();
  //define roads
  Double_t kroady = 1., //recoParam->GetRoad1y();
           kroadz = GetPadLength() * recoParam->GetRoadzMultiplicator() + 1.;
  // define probing cluster (the perfect cluster) and default calibration
  Short_t sig[] = {0, 0, 10, 30, 10, 0,0};
  AliTRDcluster cp(fDet, 6, 75, 0, sig, 0);
  if(fkReconstructor->IsHLT()) cp.SetRPhiMethod(AliTRDcluster::kCOG);
  if(!IsCalibrated()) Calibrate();

  AliDebug(4, "");
  AliDebug(4, Form("syKalman[%f] rY[%f] rZ[%f]", syRef, kroady, kroadz));

  // working variables
  const Int_t kNrows = 16;
  const Int_t kNcls  = 3*kNclusters; // buffer size
  AliTRDcluster *clst[kNrows][kNcls];
  Bool_t blst[kNrows][kNcls];
  Double_t cond[4], dx, dy, yt, zt, yres[kNrows][kNcls];
  Int_t idxs[kNrows][kNcls], ncl[kNrows], ncls = 0;
  memset(ncl, 0, kNrows*sizeof(Int_t));
  memset(yres, 0, kNrows*kNcls*sizeof(Double_t));
  memset(blst, 0, kNrows*kNcls*sizeof(Bool_t));   //this is 8 times faster to memset than "memset(clst, 0, kNrows*kNcls*sizeof(AliTRDcluster*))"

  // Do cluster projection
  AliTRDcluster *c = NULL;
  AliTRDchamberTimeBin *layer = NULL;
  Bool_t kBUFFER = kFALSE;
  for (Int_t it = 0; it < kNtb; it++) {
    if(!(layer = chamber->GetTB(it))) continue;
    if(!Int_t(*layer)) continue;
    // get track projection at layers position
    dx   = fX0 - layer->GetX();
    yt = fYref[0] - fYref[1] * dx;
    zt = fZref[0] - fZref[1] * dx;
    // get standard cluster error corrected for tilt
    cp.SetLocalTimeBin(it);
    cp.SetSigmaY2(0.02, fDiffT, fExB, dx, -1./*zt*/, fYref[1]);
    s2yCl = (cp.GetSigmaY2() + sysCov[0] + t2*s2zCl)/(1.+t2);
    // get estimated road
    kroady = 3.*TMath::Sqrt(12.*(s2yTrk + s2yCl));

    AliDebug(5, Form("  %2d x[%f] yt[%f] zt[%f]", it, dx, yt, zt));

    AliDebug(5, Form("  syTrk[um]=%6.2f syCl[um]=%6.2f syClTlt[um]=%6.2f Ry[mm]=%f", 1.e4*TMath::Sqrt(s2yTrk), 1.e4*TMath::Sqrt(cp.GetSigmaY2()), 1.e4*TMath::Sqrt(s2yCl), 1.e1*kroady));

    // select clusters
    cond[0] = yt; cond[2] = kroady;
    cond[1] = zt; cond[3] = kroadz;
    Int_t n=0, idx[6];
    layer->GetClusters(cond, idx, n, 6);
    for(Int_t ic = n; ic--;){
      c  = (*layer)[idx[ic]];
      dy = yt - c->GetY();
      dy += tilt ? GetTilt() * (c->GetZ() - zt) : 0.;
      // select clusters on a 3 sigmaKalman level
/*      if(tilt && TMath::Abs(dy) > 3.*syRef){ 
        printf("too large !!!\n");
        continue;
      }*/
      Int_t r = c->GetPadRow();
      AliDebug(5, Form("   -> dy[%f] yc[%f] r[%d]", TMath::Abs(dy), c->GetY(), r));
      clst[r][ncl[r]] = c;
      blst[r][ncl[r]] = kTRUE;
      idxs[r][ncl[r]] = idx[ic];
      yres[r][ncl[r]] = dy;
      ncl[r]++; ncls++;

      if(ncl[r] >= kNcls) {
        AliWarning(Form("Cluster candidates row[%d] reached buffer limit[%d]. Some may be lost.", r, kNcls));
        kBUFFER = kTRUE;
        break;
      }
    }
    if(kBUFFER) break;
  }
  AliDebug(4, Form("Found %d clusters. Processing ...", ncls));
  if(ncls<kClmin){ 
    AliDebug(1, Form("CLUSTERS FOUND %d LESS THAN THRESHOLD %d.", ncls, kClmin));
    SetErrorMsg(kAttachClFound);
    return kFALSE;
  }

  // analyze each row individualy
  Bool_t kRowSelection(kFALSE);
  Double_t mean[]={1.e3, 1.e3, 1.3}, syDis[]={1.e3, 1.e3, 1.3};
  Int_t nrow[] = {0, 0, 0}, rowId[] = {-1, -1, -1}, nr = 0, lr=-1;
  TVectorD vdy[3];
  for(Int_t ir=0; ir<kNrows; ir++){
    if(!(ncl[ir])) continue;
    if(lr>0 && ir-lr != 1){ 
      AliDebug(2, "Rows attached not continuous. Turn on selection."); 
      kRowSelection=kTRUE;
    }

    AliDebug(5, Form("  r[%d] n[%d]", ir, ncl[ir]));
    // Evaluate truncated mean on the y direction
    if(ncl[ir] < 4) continue;
    AliMathBase::EvaluateUni(ncl[ir], yres[ir], mean[nr], syDis[nr], Int_t(ncl[ir]*.8));

    // TODO check mean and sigma agains cluster resolution !!
    AliDebug(4, Form("  m_%d[%+5.3f (%5.3fs)] s[%f]", nr, mean[nr], TMath::Abs(mean[nr]/syDis[nr]), syDis[nr]));
    // remove outliers based on a 3 sigmaDistr level
    Bool_t kFOUND = kFALSE;
    for(Int_t ic = ncl[ir]; ic--;){
      if(yres[ir][ic] - mean[nr] > 3. * syDis[nr]){ 
        blst[ir][ic] = kFALSE; continue;
      }
      nrow[nr]++; rowId[nr]=ir; kFOUND = kTRUE;
    }
    if(kFOUND){ 
      vdy[nr].Use(nrow[nr], yres[ir]);
      nr++; 
    }
    lr = ir; if(nr>=3) break;
  }
  if(recoParam->GetStreamLevel(AliTRDrecoParam::kTracker) > 3 && fkReconstructor->IsDebugStreaming()){
    TTreeSRedirector &cstreamer = *fkReconstructor->GetDebugStream(AliTRDrecoParam::kTracker);
    UChar_t stat(0);
    if(IsKink()) SETBIT(stat, 1);
    if(IsStandAlone()) SETBIT(stat, 2);
    cstreamer << "AttachClusters"
        << "stat="   << stat
        << "det="    << fDet
        << "pt="     << fPt
        << "s2y="    << s2yTrk
        << "r0="     << rowId[0]
        << "dy0="    << &vdy[0]
        << "m0="     << mean[0]
        << "s0="     << syDis[0]
        << "r1="     << rowId[1]
        << "dy1="    << &vdy[1]
        << "m1="     << mean[1]
        << "s1="     << syDis[1]
        << "r2="     << rowId[2]
        << "dy2="    << &vdy[2]
        << "m2="     << mean[2]
        << "s2="     << syDis[2]
        << "\n";
  }


  // analyze gap in rows attached 
  if(kRowSelection){
    SetErrorMsg(kAttachRowGap);
    Int_t rowRemove(-1); 
    if(nr==2){ // select based on minimum distance to track projection
      if(TMath::Abs(mean[0])<TMath::Abs(mean[1])){ 
        if(nrow[1]>nrow[0]) AliDebug(2, Form("Conflicting mean[%f < %f] but ncl[%d < %d].", TMath::Abs(mean[0]), TMath::Abs(mean[1]), nrow[0], nrow[1]));
      }else{
        if(nrow[1]<nrow[0]) AliDebug(2, Form("Conflicting mean[%f > %f] but ncl[%d > %d].", TMath::Abs(mean[0]), TMath::Abs(mean[1]), nrow[0], nrow[1]));
        Swap(nrow[0],nrow[1]); Swap(rowId[0],rowId[1]);
        Swap(mean[0],mean[1]); Swap(syDis[0],syDis[1]);
      }
      rowRemove=1; nr=1; 
    } else if(nr==3){ // select based on 2 consecutive rows
      if(rowId[1]==rowId[0]+1 && rowId[1]!=rowId[2]-1){ 
        nr=2;rowRemove=2;
      } else if(rowId[1]!=rowId[0]+1 && rowId[1]==rowId[2]-1){ 
        Swap(nrow[0],nrow[2]); Swap(rowId[0],rowId[2]);
        Swap(mean[0],mean[2]); Swap(syDis[0],syDis[2]);
        nr=2; rowRemove=2;
      }
    }
    if(rowRemove>0){nrow[rowRemove]=0; rowId[rowRemove]=-1;}
  }
  AliDebug(4, Form("  Ncl[%d[%d] + %d[%d] + %d[%d]]", nrow[0], rowId[0],  nrow[1], rowId[1], nrow[2], rowId[2]));

  if(nr==3){
    SetBit(kRowCross, kTRUE); // mark pad row crossing
    SetErrorMsg(kAttachRow);
    const Float_t am[]={TMath::Abs(mean[0]), TMath::Abs(mean[1]), TMath::Abs(mean[2])};
    AliDebug(4, Form("complex row configuration\n"
      "  r[%d] n[%d] m[%6.3f] s[%6.3f]\n"
      "  r[%d] n[%d] m[%6.3f] s[%6.3f]\n"
      "  r[%d] n[%d] m[%6.3f] s[%6.3f]\n"
      , rowId[0], nrow[0], am[0], syDis[0]
      , rowId[1], nrow[1], am[1], syDis[1]
      , rowId[2], nrow[2], am[2], syDis[2]));
    Int_t id[]={0,1,2}; TMath::Sort(3, am, id, kFALSE);
    // backup
    Int_t rnn[3]; memcpy(rnn, nrow, 3*sizeof(Int_t));
    Int_t rid[3]; memcpy(rid, rowId, 3*sizeof(Int_t));
    Double_t rm[3]; memcpy(rm, mean, 3*sizeof(Double_t));
    Double_t rs[3]; memcpy(rs, syDis, 3*sizeof(Double_t));
    nrow[0]=rnn[id[0]]; rowId[0]=rid[id[0]]; mean[0]=rm[id[0]]; syDis[0]=rs[id[0]];
    nrow[1]=rnn[id[1]]; rowId[1]=rid[id[1]]; mean[1]=rm[id[1]]; syDis[1]=rs[id[1]];
    nrow[2]=0;          rowId[2]=-1; mean[2] = 1.e3; syDis[2] = 1.e3;
    AliDebug(4, Form("solved configuration\n"
      "  r[%d] n[%d] m[%+6.3f] s[%6.3f]\n"
      "  r[%d] n[%d] m[%+6.3f] s[%6.3f]\n"
      "  r[%d] n[%d] m[%+6.3f] s[%6.3f]\n"
      , rowId[0], nrow[0], mean[0], syDis[0]
      , rowId[1], nrow[1], mean[1], syDis[1]
      , rowId[2], nrow[2], mean[2], syDis[2]));
    nr=2;
  } else if(nr==2) {
    SetBit(kRowCross, kTRUE); // mark pad row crossing
    if(nrow[1] > nrow[0]){ // swap row order
      Swap(nrow[0],nrow[1]); Swap(rowId[0],rowId[1]);
      Swap(mean[0],mean[1]); Swap(syDis[0],syDis[1]);
    }
  }

  // Select and store clusters 
  // We should consider here :
  //  1. How far is the chamber boundary
  //  2. How big is the mean
  Int_t n(0); Float_t dyc[kNclusters]; memset(dyc,0,kNclusters*sizeof(Float_t));
  for (Int_t ir = 0; ir < nr; ir++) {
    Int_t jr(rowId[ir]);
    AliDebug(4, Form("  Attaching Ncl[%d]=%d ...", jr, ncl[jr]));
    for (Int_t ic = 0; ic < ncl[jr]; ic++) {
      if(!blst[jr][ic])continue;
      c = clst[jr][ic];
      Int_t it(c->GetPadTime());
      Int_t idx(it+kNtb*ir);
      if(fClusters[idx]){
        AliDebug(4, Form("Many cluster candidates on row[%2d] tb[%2d].", jr, it));
        // TODO should save also the information on where the multiplicity happened and its size
        SetErrorMsg(kAttachMultipleCl);
        // TODO should also compare with mean and sigma for this row
        if(yres[jr][ic] > dyc[idx]) continue;
      }

      // TODO proper indexing of clusters !!
      fIndexes[idx]  = chamber->GetTB(it)->GetGlobalIndex(idxs[jr][ic]);
      fClusters[idx] = c;
      dyc[idx]        = yres[jr][ic];
      n++;
    }
  }
  SetN(n);

  // number of minimum numbers of clusters expected for the tracklet
  if (GetN() < kClmin){
    AliDebug(1, Form("NOT ENOUGH CLUSTERS %d ATTACHED TO THE TRACKLET [min %d] FROM FOUND %d.", GetN(), kClmin, n));
    SetErrorMsg(kAttachClAttach);
    return kFALSE;
  }

  // Load calibration parameters for this tracklet  
  Calibrate();

  // calculate dx for time bins in the drift region (calibration aware)
  Float_t x[2] = {0.,0.}; Int_t tb[2]={0,0};
  for (Int_t it = t0, irp=0; irp<2 && it < AliTRDtrackerV1::GetNTimeBins(); it++) {
    if(!fClusters[it]) continue;
    x[irp]  = fClusters[it]->GetX();
    tb[irp] = fClusters[it]->GetLocalTimeBin();
    irp++;
  }  
  Int_t dtb = tb[1] - tb[0];
  fdX = dtb ? (x[0] - x[1]) / dtb : 0.15;
  return kTRUE;
}

//____________________________________________________________
void AliTRDseedV1::Bootstrap(const AliTRDReconstructor *rec)
{
//   Fill in all derived information. It has to be called after recovery from file or HLT.
//   The primitive data are
//   - list of clusters
//   - detector (as the detector will be removed from clusters)
//   - position of anode wire (fX0) - temporary
//   - track reference position and direction
//   - momentum of the track
//   - time bin length [cm]
// 
//   A.Bercuci <A.Bercuci@gsi.de> Oct 30th 2008
//
  fkReconstructor = rec;
  AliTRDgeometry g;
  SetPadPlane(g.GetPadPlane(fDet));

  //fSnp = fYref[1]/TMath::Sqrt(1+fYref[1]*fYref[1]);
  //fTgl = fZref[1];
  Int_t n = 0, nshare = 0, nused = 0;
  AliTRDcluster **cit = &fClusters[0];
  for(Int_t ic = kNclusters; ic--; cit++){
    if(!(*cit)) return;
    n++;
    if((*cit)->IsShared()) nshare++;
    if((*cit)->IsUsed()) nused++;
  }
  SetN(n); SetNUsed(nused); SetNShared(nshare);
  Fit();
  CookLabels();
  GetProbability();
}


//____________________________________________________________________
Bool_t AliTRDseedV1::Fit(UChar_t opt)
{
//
// Linear fit of the clusters attached to the tracklet
//
// Parameters :
//   - opt : switch for tilt pad correction of cluster y position. Options are
//           0 no correction [default]
//           1 full tilt correction [dz/dx and z0]
//           2 pseudo tilt correction [dz/dx from pad-chamber geometry]
//
// Output :
//  True if successful
//
// Detailed description
//
//            Fit in the xy plane
// 
// The fit is performed to estimate the y position of the tracklet and the track 
// angle in the bending plane. The clusters are represented in the chamber coordinate 
// system (with respect to the anode wire - see AliTRDtrackerV1::FollowBackProlongation() 
// on how this is set). The x and y position of the cluster and also their variances 
// are known from clusterizer level (see AliTRDcluster::GetXloc(), AliTRDcluster::GetYloc(), 
// AliTRDcluster::GetSX() and AliTRDcluster::GetSY()). 
// If gaussian approximation is used to calculate y coordinate of the cluster the position 
// is recalculated taking into account the track angle. The general formula to calculate the 
// error of cluster position in the gaussian approximation taking into account diffusion and track
// inclination is given for TRD by:
// BEGIN_LATEX
// #sigma^{2}_{y} = #sigma^{2}_{PRF} + #frac{x#delta_{t}^{2}}{(1+tg(#alpha_{L}))^{2}} + #frac{x^{2}tg^{2}(#phi-#alpha_{L})tg^{2}(#alpha_{L})}{12}
// END_LATEX
//
// Since errors are calculated only in the y directions, radial errors (x direction) are mapped to y
// by projection i.e.
// BEGIN_LATEX
// #sigma_{x|y} = tg(#phi) #sigma_{x}
// END_LATEX
// and also by the lorentz angle correction
//
//            Fit in the xz plane
//
// The "fit" is performed to estimate the radial position (x direction) where pad row cross happens. 
// If no pad row crossing the z position is taken from geometry and radial position is taken from the xy 
// fit (see below).
// 
// There are two methods to estimate the radial position of the pad row cross:
//   1. leading cluster radial position : Here the lower part of the tracklet is considered and the last 
// cluster registered (at radial x0) on this segment is chosen to mark the pad row crossing. The error 
// of the z estimate is given by :
// BEGIN_LATEX
// #sigma_{z} = tg(#theta) #Delta x_{x_{0}}/12
// END_LATEX
// The systematic errors for this estimation are generated by the following sources:
//   - no charge sharing between pad rows is considered (sharp cross)
//   - missing cluster at row cross (noise peak-up, under-threshold signal etc.).
// 
//   2. charge fit over the crossing point : Here the full energy deposit along the tracklet is considered 
// to estimate the position of the crossing by a fit in the qx plane. The errors in the q directions are 
// parameterized as s_q = q^2. The systematic errors for this estimation are generated by the following sources:
//   - no general model for the qx dependence
//   - physical fluctuations of the charge deposit 
//   - gain calibration dependence
//
//            Estimation of the radial position of the tracklet
//
// For pad row cross the radial position is taken from the xz fit (see above). Otherwise it is taken as the 
// interpolation point of the tracklet i.e. the point where the error in y of the fit is minimum. The error
// in the y direction of the tracklet is (see AliTRDseedV1::GetCovAt()):
// BEGIN_LATEX
// #sigma_{y} = #sigma^{2}_{y_{0}} + 2xcov(y_{0}, dy/dx) + #sigma^{2}_{dy/dx}
// END_LATEX
// and thus the radial position is:
// BEGIN_LATEX
// x = - cov(y_{0}, dy/dx)/#sigma^{2}_{dy/dx}
// END_LATEX
//
//            Estimation of tracklet position error 
//
// The error in y direction is the error of the linear fit at the radial position of the tracklet while in the z 
// direction is given by the cluster error or pad row cross error. In case of no pad row cross this is given by:
// BEGIN_LATEX
// #sigma_{y} = #sigma^{2}_{y_{0}} - 2cov^{2}(y_{0}, dy/dx)/#sigma^{2}_{dy/dx} + #sigma^{2}_{dy/dx}
// #sigma_{z} = Pad_{length}/12
// END_LATEX
// For pad row cross the full error is calculated at the radial position of the crossing (see above) and the error 
// in z by the width of the crossing region - being a matter of parameterization. 
// BEGIN_LATEX
// #sigma_{z} = tg(#theta) #Delta x_{x_{0}}/12
// END_LATEX
// In case of no tilt correction (default in the barrel tracking) the tilt is taken into account by the rotation of
// the covariance matrix. See AliTRDseedV1::GetCovAt() for details.
//
// Author 
// A.Bercuci <A.Bercuci@gsi.de>

  if(!fkReconstructor){
    AliError("The tracklet needs the reconstruction setup. Please initialize by SetReconstructor().");
    return kFALSE;
  }
  if(!IsCalibrated()) Calibrate();
  if(opt>2){
    AliWarning(Form("Option [%d] outside range [0, 2]. Using default"));
    opt=0;
  }

  const Int_t kClmin = 8;
  const Float_t kScalePulls = 10.; // factor to scale y pulls - NOT UNDERSTOOD
  // get track direction
  Double_t y0   = fYref[0];
  Double_t dydx = fYref[1]; 
  Double_t z0   = fZref[0];
  Double_t dzdx = fZref[1];

  AliTRDtrackerV1::AliTRDLeastSquare fitterY;
  AliTRDtrackerV1::AliTRDLeastSquare fitterZ;

  // book cluster information
  Double_t qc[kNclusters], xc[kNclusters], yc[kNclusters], zc[kNclusters], sy[kNclusters];

  Bool_t tilt(opt==1)       // full tilt correction
        ,pseudo(opt==2)     // pseudo tilt correction
        ,rc(IsRowCross())   // row cross candidate 
        ,kDZDX(IsPrimary());// switch dzdx calculation for barrel primary tracks
  Int_t n(0);   // clusters used in fit 
  AliTRDcluster *c(NULL), *cc(NULL), **jc = &fClusters[0];
  const AliTRDrecoParam* const recoParam = fkReconstructor->GetRecoParam(); //the dynamic cast in GetRecoParam is slow, so caching the pointer to it

  const Char_t *tcName[]={"NONE", "FULL", "HALF"};
  AliDebug(2, Form("Options : TC[%s] dzdx[%c]", tcName[opt], kDZDX?'Y':'N'));

  for (Int_t ic=0; ic<kNclusters; ic++, ++jc) {
    xc[ic]  = -1.; yc[ic]  = 999.; zc[ic]  = 999.; sy[ic]  = 0.;
    if(!(c = (*jc))) continue;
    if(!c->IsInChamber()) continue;
    // compute pseudo tilt correction
    if(kDZDX){ 
      fZfit[0] = c->GetZ();
      if(rc){
        for(Int_t kc=AliTRDseedV1::kNtb; kc<AliTRDseedV1::kNclusters; kc++){
          if(!(cc=fClusters[kc])) continue;
          if(!cc->IsInChamber()) continue;
          fZfit[0] += cc->GetZ(); fZfit[0] *= 0.5;
          break;
        }
      }
      fZfit[1] = fZfit[0]/fX0;
      if(rc){
        fZfit[0] += fZfit[1]*0.5*AliTRDgeometry::CdrHght();
        fZfit[1] = fZfit[0]/fX0;
      }
      kDZDX=kFALSE;
    }

    Float_t w = 1.;
    if(c->GetNPads()>4) w = .5;
    if(c->GetNPads()>5) w = .2;

    // cluster charge
    qc[n]   = TMath::Abs(c->GetQ());
    // pad row of leading 

    xc[n]   = fX0 - c->GetX();

    // Recalculate cluster error based on tracking information
    c->SetSigmaY2(fS2PRF, fDiffT, fExB, xc[n], -1./*zcorr?zt:-1.*/, dydx);
    c->SetSigmaZ2(fPad[0]*fPad[0]/12.); // for HLT
    sy[n]  = TMath::Sqrt(c->GetSigmaY2());

    yc[n]  = recoParam->UseGAUS() ? 
      c->GetYloc(y0, sy[n], GetPadWidth()): c->GetY();
    zc[n]   = c->GetZ();

    //optional r-phi correction
    //printf("   n[%2d] yc[%7.5f] ", n, yc[n]);
    Float_t correction(0.);
    if(tilt) correction = fPad[2]*(xc[n]*dzdx + zc[n] - z0);
    else if(pseudo) correction = fPad[2]*(xc[n]*fZfit[1] + zc[n]-fZfit[0]);
    yc[n]-=correction;
    //printf("corr(%s%s)[%7.5f] yc1[%7.5f]\n", (tilt?"TC":""), (zcorr?"PC":""), correction, yc[n]);

    AliDebug(5, Form("  tb[%2d] dx[%6.3f] y[%6.2f+-%6.3f]", c->GetLocalTimeBin(), xc[n], yc[n], sy[n]));
    fitterY.AddPoint(&xc[n], yc[n], sy[n]);
    if(rc) fitterZ.AddPoint(&xc[n], qc[n]*(ic<kNtb?1.:-1.), 1.);
    n++;
  }

  // to few clusters
  if (n < kClmin){ 
    SetErrorMsg(kFitCl);
    return kFALSE; 
  }

  // fit XY
  if(!fitterY.Eval()){
    SetErrorMsg(kFitFailedY);
    return kFALSE;
  }
  fYfit[0] = fitterY.GetFunctionParameter(0);
  fYfit[1] = -fitterY.GetFunctionParameter(1);
  // store covariance
  Double_t p[3];
  fitterY.GetCovarianceMatrix(p);
  fCov[0] = kScalePulls*p[1]; // variance of y0
  fCov[1] = kScalePulls*p[2]; // covariance of y0, dydx
  fCov[2] = kScalePulls*p[0]; // variance of dydx
  // the ref radial position is set at the minimum of 
  // the y variance of the tracklet
  fX   = -fCov[1]/fCov[2];
  fS2Y = fCov[0] +2.*fX*fCov[1] + fX*fX*fCov[2];

  Float_t xs=fX+.5*AliTRDgeometry::CamHght();
  if(xs < 0. || xs > AliTRDgeometry::CamHght()+AliTRDgeometry::CdrHght()){
    AliDebug(1, Form("Ref radial position ouside chamber x[%5.2f].", fX));
    SetErrorMsg(kFitFailedY);
    return kFALSE;
  }

/*    // THE LEADING CLUSTER METHOD for z fit
    Float_t xMin = fX0;
    Int_t ic=n=kNclusters-1; jc = &fClusters[ic];
    AliTRDcluster *c0 =0x0, **kc = &fClusters[kNtb-1];
    for(; ic>kNtb; ic--, --jc, --kc){
      if((c0 = (*kc)) && c0->IsInChamber() && (xMin>c0->GetX())) xMin = c0->GetX();
      if(!(c = (*jc))) continue;
      if(!c->IsInChamber()) continue;
      zc[kNclusters-1] = c->GetZ(); 
      fX = fX0 - c->GetX();
    }
    fZfit[0] = .5*(zc[0]+zc[kNclusters-1]); fZfit[1] = 0.;
    // Error parameterization
    fS2Z     = fdX*fZref[1];
    fS2Z    *= fS2Z; fS2Z    *= 0.2887; //  1/sqrt(12)*/

  // fit QZ
  if(opt!=1 && IsRowCross()){
    if(!fitterZ.Eval()) SetErrorMsg(kFitFailedZ);
    if(!HasError(kFitFailedZ) && fitterZ.GetFunctionParameter(1)!=0.){ 
      // TODO - one has to recalculate xy fit based on
      // better knowledge of z position
//       Double_t x = -fitterZ.GetFunctionParameter(0)/fitterZ.GetFunctionParameter(1);
//       Double_t z0 = .5*(zc[0]+zc[n-1]);
//       fZfit[0] = z0 + fZfit[1]*x; 
//       fZfit[1] = fZfit[0]/fX0; 
//       redo fit on xy plane
    }
    // temporary external error parameterization
    fS2Z     = 0.05+0.4*TMath::Abs(fZref[1]); fS2Z *= fS2Z;
    // TODO correct formula
    //fS2Z     = sigma_x*TMath::Abs(fZref[1]);
  } else {
    //fZfit[0] = zc[0] + dzdx*0.5*AliTRDgeometry::CdrHght();
    fS2Z     = GetPadLength()*GetPadLength()/12.;
  }
  return kTRUE;
}


/*
//_____________________________________________________________________________
void AliTRDseedV1::FitMI()
{
//
// Fit the seed.
// Marian Ivanov's version 
//
// linear fit on the y direction with respect to the reference direction. 
// The residuals for each x (x = xc - x0) are deduced from:
// dy = y - yt             (1)
// the tilting correction is written :
// y = yc + h*(zc-zt)      (2)
// yt = y0+dy/dx*x         (3)
// zt = z0+dz/dx*x         (4)
// from (1),(2),(3) and (4)
// dy = yc - y0 - (dy/dx + h*dz/dx)*x + h*(zc-z0)
// the last term introduces the correction on y direction due to tilting pads. There are 2 ways to account for this:
// 1. use tilting correction for calculating the y
// 2. neglect tilting correction here and account for it in the error parametrization of the tracklet.
  const Float_t kRatio  = 0.8;
  const Int_t   kClmin  = 5;
  const Float_t kmaxtan = 2;

  if (TMath::Abs(fYref[1]) > kmaxtan){
		//printf("Exit: Abs(fYref[1]) = %3.3f, kmaxtan = %3.3f\n", TMath::Abs(fYref[1]), kmaxtan);
		return;              // Track inclined too much
	}

  Float_t  sigmaexp  = 0.05 + TMath::Abs(fYref[1] * 0.25); // Expected r.m.s in y direction
  Float_t  ycrosscor = GetPadLength() * GetTilt() * 0.5;           // Y correction for crossing 
  Int_t fNChange = 0;

  Double_t sumw;
  Double_t sumwx;
  Double_t sumwx2;
  Double_t sumwy;
  Double_t sumwxy;
  Double_t sumwz;
  Double_t sumwxz;

	// Buffering: Leave it constant fot Performance issues
  Int_t    zints[kNtb];            // Histograming of the z coordinate 
                                         // Get 1 and second max probable coodinates in z
  Int_t    zouts[2*kNtb];       
  Float_t  allowedz[kNtb];         // Allowed z for given time bin
  Float_t  yres[kNtb];             // Residuals from reference
  //Float_t  anglecor = GetTilt() * fZref[1];  // Correction to the angle
  
  Float_t pos[3*kNtb]; memset(pos, 0, 3*kNtb*sizeof(Float_t));
  Float_t *fX = &pos[0], *fY = &pos[kNtb], *fZ = &pos[2*kNtb];
  
  Int_t fN  = 0; AliTRDcluster *c = 0x0; 
  fN2 = 0;
  for (Int_t i = 0; i < AliTRDtrackerV1::GetNTimeBins(); i++) {
    yres[i] = 10000.0;
    if (!(c = fClusters[i])) continue;
    if(!c->IsInChamber()) continue;
    // Residual y
    //yres[i] = fY[i] - fYref[0] - (fYref[1] + anglecor) * fX[i] + GetTilt()*(fZ[i] - fZref[0]);
    fX[i] = fX0 - c->GetX();
    fY[i] = c->GetY();
    fZ[i] = c->GetZ();
    yres[i] = fY[i] - GetTilt()*(fZ[i] - (fZref[0] - fX[i]*fZref[1]));
    zints[fN] = Int_t(fZ[i]);
    fN++;
  }

  if (fN < kClmin){
    //printf("Exit fN < kClmin: fN = %d\n", fN);
    return; 
  }
  Int_t nz = AliTRDtrackerV1::Freq(fN, zints, zouts, kFALSE);
  Float_t fZProb   = zouts[0];
  if (nz <= 1) zouts[3] = 0;
  if (zouts[1] + zouts[3] < kClmin) {
    //printf("Exit zouts[1] = %d, zouts[3] = %d\n",zouts[1],zouts[3]);
    return;
  }
  
  // Z distance bigger than pad - length
  if (TMath::Abs(zouts[0]-zouts[2]) > 12.0) zouts[3] = 0;
  
  Int_t  breaktime = -1;
  Bool_t mbefore   = kFALSE;
  Int_t  cumul[kNtb][2];
  Int_t  counts[2] = { 0, 0 };
  
  if (zouts[3] >= 3) {

    //
    // Find the break time allowing one chage on pad-rows
    // with maximal number of accepted clusters
    //
    fNChange = 1;
    for (Int_t i = 0; i < AliTRDtrackerV1::GetNTimeBins(); i++) {
      cumul[i][0] = counts[0];
      cumul[i][1] = counts[1];
      if (TMath::Abs(fZ[i]-zouts[0]) < 2) counts[0]++;
      if (TMath::Abs(fZ[i]-zouts[2]) < 2) counts[1]++;
    }
    Int_t  maxcount = 0;
    for (Int_t i = 0; i < AliTRDtrackerV1::GetNTimeBins(); i++) {
      Int_t after  = cumul[AliTRDtrackerV1::GetNTimeBins()][0] - cumul[i][0];
      Int_t before = cumul[i][1];
      if (after + before > maxcount) { 
        maxcount  = after + before; 
        breaktime = i;
        mbefore   = kFALSE;
      }
      after  = cumul[AliTRDtrackerV1::GetNTimeBins()-1][1] - cumul[i][1];
      before = cumul[i][0];
      if (after + before > maxcount) { 
        maxcount  = after + before; 
        breaktime = i;
        mbefore   = kTRUE;
      }
    }
    breaktime -= 1;
  }

  for (Int_t i = 0; i < AliTRDtrackerV1::GetNTimeBins()+1; i++) {
    if (i >  breaktime) allowedz[i] =   mbefore  ? zouts[2] : zouts[0];
    if (i <= breaktime) allowedz[i] = (!mbefore) ? zouts[2] : zouts[0];
  }  

  if (((allowedz[0] > allowedz[AliTRDtrackerV1::GetNTimeBins()]) && (fZref[1] < 0)) ||
      ((allowedz[0] < allowedz[AliTRDtrackerV1::GetNTimeBins()]) && (fZref[1] > 0))) {
    //
    // Tracklet z-direction not in correspondance with track z direction 
    //
    fNChange = 0;
    for (Int_t i = 0; i < AliTRDtrackerV1::GetNTimeBins()+1; i++) {
      allowedz[i] = zouts[0];  // Only longest taken
    } 
  }
  
  if (fNChange > 0) {
    //
    // Cross pad -row tracklet  - take the step change into account
    //
    for (Int_t i = 0; i < AliTRDtrackerV1::GetNTimeBins()+1; i++) {
      if (!fClusters[i]) continue; 
      if(!fClusters[i]->IsInChamber()) continue;
      if (TMath::Abs(fZ[i] - allowedz[i]) > 2) continue;
      // Residual y
      //yres[i] = fY[i] - fYref[0] - (fYref[1] + anglecor) * fX[i] + GetTilt()*(fZ[i] - fZref[0]);   
      yres[i] = fY[i] - GetTilt()*(fZ[i] - (fZref[0] - fX[i]*fZref[1]));
//       if (TMath::Abs(fZ[i] - fZProb) > 2) {
//         if (fZ[i] > fZProb) yres[i] += GetTilt() * GetPadLength();
//         if (fZ[i] < fZProb) yres[i] -= GetTilt() * GetPadLength();
      }
    }
  }
  
  Double_t yres2[kNtb];
  Double_t mean;
  Double_t sigma;
  for (Int_t i = 0; i < AliTRDtrackerV1::GetNTimeBins()+1; i++) {
    if (!fClusters[i]) continue;
    if(!fClusters[i]->IsInChamber()) continue;
    if (TMath::Abs(fZ[i] - allowedz[i]) > 2) continue;
    yres2[fN2] = yres[i];
    fN2++;
  }
  if (fN2 < kClmin) {
		//printf("Exit fN2 < kClmin: fN2 = %d\n", fN2);
    fN2 = 0;
    return;
  }
  AliMathBase::EvaluateUni(fN2,yres2,mean,sigma, Int_t(fN2*kRatio-2.));
  if (sigma < sigmaexp * 0.8) {
    sigma = sigmaexp;
  }
  //Float_t fSigmaY = sigma;

  // Reset sums
  sumw   = 0; 
  sumwx  = 0; 
  sumwx2 = 0;
  sumwy  = 0; 
  sumwxy = 0; 
  sumwz  = 0;
  sumwxz = 0;

  fN2    = 0;
  Float_t fMeanz = 0;
  Float_t fMPads = 0;
  fUsable = 0;
  for (Int_t i = 0; i < AliTRDtrackerV1::GetNTimeBins()+1; i++) {
    if (!fClusters[i]) continue;
    if (!fClusters[i]->IsInChamber()) continue;
    if (TMath::Abs(fZ[i] - allowedz[i]) > 2){fClusters[i] = 0x0; continue;}
    if (TMath::Abs(yres[i] - mean) > 4.0 * sigma){fClusters[i] = 0x0;  continue;}
    SETBIT(fUsable,i);
    fN2++;
    fMPads += fClusters[i]->GetNPads();
    Float_t weight = 1.0;
    if (fClusters[i]->GetNPads() > 4) weight = 0.5;
    if (fClusters[i]->GetNPads() > 5) weight = 0.2;
   
   	
    Double_t x = fX[i];
    //printf("x = %7.3f dy = %7.3f fit %7.3f\n", x, yres[i], fY[i]-yres[i]);
    
    sumw   += weight; 
    sumwx  += x * weight; 
    sumwx2 += x*x * weight;
    sumwy  += weight * yres[i];  
    sumwxy += weight * (yres[i]) * x;
    sumwz  += weight * fZ[i];    
    sumwxz += weight * fZ[i] * x;

  }

  if (fN2 < kClmin){
		//printf("Exit fN2 < kClmin(2): fN2 = %d\n",fN2);
    fN2 = 0;
    return;
  }
  fMeanz = sumwz / sumw;
  Float_t correction = 0;
  if (fNChange > 0) {
    // Tracklet on boundary
    if (fMeanz < fZProb) correction =  ycrosscor;
    if (fMeanz > fZProb) correction = -ycrosscor;
  }

  Double_t det = sumw * sumwx2 - sumwx * sumwx;
  fYfit[0]    = (sumwx2 * sumwy  - sumwx * sumwxy) / det;
  fYfit[1]    = (sumw   * sumwxy - sumwx * sumwy)  / det;
  
  fS2Y = 0;
  for (Int_t i = 0; i < AliTRDtrackerV1::GetNTimeBins()+1; i++) {
    if (!TESTBIT(fUsable,i)) continue;
    Float_t delta = yres[i] - fYfit[0] - fYfit[1] * fX[i];
    fS2Y += delta*delta;
  }
  fS2Y = TMath::Sqrt(fS2Y / Float_t(fN2-2));
	// TEMPORARY UNTIL covariance properly calculated
	fS2Y = TMath::Max(fS2Y, Float_t(.1));
  
  fZfit[0]   = (sumwx2 * sumwz  - sumwx * sumwxz) / det;
  fZfit[1]   = (sumw   * sumwxz - sumwx * sumwz)  / det;
//   fYfitR[0] += fYref[0] + correction;
//   fYfitR[1] += fYref[1];
//  fYfit[0]   = fYfitR[0];
  fYfit[1]   = -fYfit[1];

  UpdateUsed();
}*/

//___________________________________________________________________
void AliTRDseedV1::Print(Option_t *o) const
{
  //
  // Printing the seedstatus
  //

  AliInfo(Form("Det[%3d] X0[%7.2f] Pad{L[%5.2f] W[%5.2f] Tilt[%+6.2f]}", fDet, fX0, GetPadLength(), GetPadWidth(), GetTilt()));
  AliInfo(Form("N[%2d] Nused[%2d] Nshared[%2d] [%d]", GetN(), GetNUsed(), GetNShared(), fN));
  AliInfo(Form("FLAGS : RC[%c] Kink[%c] SA[%c]", IsRowCross()?'y':'n', IsKink()?'y':'n', IsStandAlone()?'y':'n'));
  AliInfo(Form("CALIB PARAMS :  T0[%5.2f]  Vd[%5.2f]  s2PRF[%5.2f]  ExB[%5.2f]  Dl[%5.2f]  Dt[%5.2f]", fT0, fVD, fS2PRF, fExB, fDiffL, fDiffT));

  Double_t cov[3], x=GetX();
  GetCovAt(x, cov);
  AliInfo("    |  x[cm]  |      y[cm]       |      z[cm]      |  dydx |  dzdx |");
  AliInfo(Form("Fit | %7.2f | %7.2f+-%7.2f | %7.2f+-%7.2f| %5.2f | ----- |", x, GetY(), TMath::Sqrt(cov[0]), GetZ(), TMath::Sqrt(cov[2]), fYfit[1]));
  AliInfo(Form("Ref | %7.2f | %7.2f+-%7.2f | %7.2f+-%7.2f| %5.2f | %5.2f |", x, fYref[0]-fX*fYref[1], TMath::Sqrt(fRefCov[0]), fZref[0]-fX*fYref[1], TMath::Sqrt(fRefCov[2]), fYref[1], fZref[1]))
  AliInfo(Form("P / Pt [GeV/c] = %f / %f", GetMomentum(), fPt));
  if(IsStandAlone()) AliInfo(Form("C Rieman / Vertex [1/cm] = %f / %f", fC[0], fC[1]));
  AliInfo(Form("dEdx [a.u.]    = %f / %f / %f / %f / %f/ %f / %f / %f", fdEdx[0], fdEdx[1], fdEdx[2], fdEdx[3], fdEdx[4], fdEdx[5], fdEdx[6], fdEdx[7]));
  AliInfo(Form("PID            = %5.3f / %5.3f / %5.3f / %5.3f / %5.3f", fProb[0], fProb[1], fProb[2], fProb[3], fProb[4]));

  if(strcmp(o, "a")!=0) return;

  AliTRDcluster* const* jc = &fClusters[0];
  for(int ic=0; ic<kNclusters; ic++, jc++) {
    if(!(*jc)) continue;
    (*jc)->Print(o);
  }
}


//___________________________________________________________________
Bool_t AliTRDseedV1::IsEqual(const TObject *o) const
{
  // Checks if current instance of the class has the same essential members
  // as the given one

  if(!o) return kFALSE;
  const AliTRDseedV1 *inTracklet = dynamic_cast<const AliTRDseedV1*>(o);
  if(!inTracklet) return kFALSE;

  for (Int_t i = 0; i < 2; i++){
    if ( fYref[i] != inTracklet->fYref[i] ) return kFALSE;
    if ( fZref[i] != inTracklet->fZref[i] ) return kFALSE;
  }
  
  if ( fS2Y != inTracklet->fS2Y ) return kFALSE;
  if ( GetTilt() != inTracklet->GetTilt() ) return kFALSE;
  if ( GetPadLength() != inTracklet->GetPadLength() ) return kFALSE;
  
  for (Int_t i = 0; i < kNclusters; i++){
//     if ( fX[i] != inTracklet->GetX(i) ) return kFALSE;
//     if ( fY[i] != inTracklet->GetY(i) ) return kFALSE;
//     if ( fZ[i] != inTracklet->GetZ(i) ) return kFALSE;
    if ( fIndexes[i] != inTracklet->fIndexes[i] ) return kFALSE;
  }
//   if ( fUsable != inTracklet->fUsable ) return kFALSE;

  for (Int_t i=0; i < 2; i++){
    if ( fYfit[i] != inTracklet->fYfit[i] ) return kFALSE;
    if ( fZfit[i] != inTracklet->fZfit[i] ) return kFALSE;
    if ( fLabels[i] != inTracklet->fLabels[i] ) return kFALSE;
  }
  
/*  if ( fMeanz != inTracklet->GetMeanz() ) return kFALSE;
  if ( fZProb != inTracklet->GetZProb() ) return kFALSE;*/
  if ( fN != inTracklet->fN ) return kFALSE;
  //if ( fNUsed != inTracklet->fNUsed ) return kFALSE;
  //if ( fFreq != inTracklet->GetFreq() ) return kFALSE;
  //if ( fNChange != inTracklet->GetNChange() ) return kFALSE;
   
  if ( fC != inTracklet->fC ) return kFALSE;
  //if ( fCC != inTracklet->GetCC() ) return kFALSE;
  if ( fChi2 != inTracklet->fChi2 ) return kFALSE;
  //  if ( fChi2Z != inTracklet->GetChi2Z() ) return kFALSE;

  if ( fDet != inTracklet->fDet ) return kFALSE;
  if ( fPt != inTracklet->fPt ) return kFALSE;
  if ( fdX != inTracklet->fdX ) return kFALSE;
  
  for (Int_t iCluster = 0; iCluster < kNclusters; iCluster++){
    AliTRDcluster *curCluster = fClusters[iCluster];
    AliTRDcluster *inCluster = inTracklet->fClusters[iCluster];
    if (curCluster && inCluster){
      if (! curCluster->IsEqual(inCluster) ) {
        curCluster->Print();
        inCluster->Print();
        return kFALSE;
      }
    } else {
      // if one cluster exists, and corresponding 
      // in other tracklet doesn't - return kFALSE
      if(curCluster || inCluster) return kFALSE;
    }
  }
  return kTRUE;
}

