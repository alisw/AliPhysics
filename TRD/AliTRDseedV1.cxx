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
//                                                                        //
//  The TRD track seed                                                    //
//                                                                        //
//  Authors:                                                              //
//    Alex Bercuci <A.Bercuci@gsi.de>                                     //
//    Markus Fasel <M.Fasel@gsi.de>                                       //
//                                                                        //
////////////////////////////////////////////////////////////////////////////

#include "TMath.h"
#include "TLinearFitter.h"
#include "TClonesArray.h" // tmp
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
#include "AliTRDReconstructor.h"
#include "AliTRDrecoParam.h"
#include "AliTRDCommonParam.h"

#include "Cal/AliTRDCalPID.h"
#include "Cal/AliTRDCalROC.h"
#include "Cal/AliTRDCalDet.h"

ClassImp(AliTRDseedV1)

//____________________________________________________________________
AliTRDseedV1::AliTRDseedV1(Int_t det) 
  :AliTRDtrackletBase()
  ,fReconstructor(0x0)
  ,fClusterIter(0x0)
  ,fExB(0.)
  ,fVD(0.)
  ,fT0(0.)
  ,fS2PRF(0.)
  ,fDiffL(0.)
  ,fDiffT(0.)
  ,fClusterIdx(0)
  ,fN(0)
  ,fDet(det)
  ,fTilt(0.)
  ,fPadLength(0.)
  ,fMom(0.)
  ,fdX(0.)
  ,fX0(0.)
  ,fX(0.)
  ,fY(0.)
  ,fZ(0.)
  ,fS2Y(0.)
  ,fS2Z(0.)
  ,fC(0.)
  ,fChi2(0.)
{
  //
  // Constructor
  //
  for(Int_t ic=kNTimeBins; ic--;) fIndexes[ic] = -1;
  memset(fClusters, 0, kNTimeBins*sizeof(AliTRDcluster*));
  fYref[0] = 0.; fYref[1] = 0.; 
  fZref[0] = 0.; fZref[1] = 0.; 
  fYfit[0] = 0.; fYfit[1] = 0.; 
  fZfit[0] = 0.; fZfit[1] = 0.; 
  memset(fdEdx, 0, kNSlices*sizeof(Float_t)); 
  for(int ispec=0; ispec<AliPID::kSPECIES; ispec++) fProb[ispec]  = -1.;
  fLabels[0]=-1; fLabels[1]=-1; // most freq MC labels
  fLabels[2]=0;  // number of different labels for tracklet
  fRefCov[0] = 1.; fRefCov[1] = 0.; fRefCov[2] = 1.;
  // covariance matrix [diagonal]
  // default sy = 200um and sz = 2.3 cm 
  fCov[0] = 4.e-4; fCov[1] = 0.; fCov[2] = 5.3; 
  SetStandAlone(kFALSE);
}

//____________________________________________________________________
AliTRDseedV1::AliTRDseedV1(const AliTRDseedV1 &ref)
  :AliTRDtrackletBase((AliTRDtrackletBase&)ref)
  ,fReconstructor(0x0)
  ,fClusterIter(0x0)
  ,fExB(0.)
  ,fVD(0.)
  ,fT0(0.)
  ,fS2PRF(0.)
  ,fDiffL(0.)
  ,fDiffT(0.)
  ,fClusterIdx(0)
  ,fN(0)
  ,fDet(-1)
  ,fTilt(0.)
  ,fPadLength(0.)
  ,fMom(0.)
  ,fdX(0.)
  ,fX0(0.)
  ,fX(0.)
  ,fY(0.)
  ,fZ(0.)
  ,fS2Y(0.)
  ,fS2Z(0.)
  ,fC(0.)
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
    for(int itb=0; itb<kNTimeBins; itb++){
      if(!fClusters[itb]) continue; 
      //AliInfo(Form("deleting c %p @ %d", fClusters[itb], itb));
      delete fClusters[itb];
      fClusters[itb] = 0x0;
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

  target.fReconstructor = fReconstructor;
  target.fClusterIter   = 0x0;
  target.fExB           = fExB;
  target.fVD            = fVD;
  target.fT0            = fT0;
  target.fS2PRF         = fS2PRF;
  target.fDiffL         = fDiffL;
  target.fDiffT         = fDiffT;
  target.fClusterIdx    = 0;
  target.fN             = fN;
  target.fDet           = fDet;
  target.fTilt          = fTilt;
  target.fPadLength     = fPadLength;
  target.fMom           = fMom;
  target.fdX            = fdX;
  target.fX0            = fX0;
  target.fX             = fX;
  target.fY             = fY;
  target.fZ             = fZ;
  target.fS2Y           = fS2Y;
  target.fS2Z           = fS2Z;
  target.fC             = fC;
  target.fChi2          = fChi2;
  
  memcpy(target.fIndexes, fIndexes, kNTimeBins*sizeof(Int_t));
  memcpy(target.fClusters, fClusters, kNTimeBins*sizeof(AliTRDcluster*));
  target.fYref[0] = fYref[0]; target.fYref[1] = fYref[1]; 
  target.fZref[0] = fZref[0]; target.fZref[1] = fZref[1]; 
  target.fYfit[0] = fYfit[0]; target.fYfit[1] = fYfit[1]; 
  target.fZfit[0] = fZfit[0]; target.fZfit[1] = fZfit[1]; 
  memcpy(target.fdEdx, fdEdx, kNSlices*sizeof(Float_t)); 
  memcpy(target.fProb, fProb, AliPID::kSPECIES*sizeof(Float_t)); 
  memcpy(target.fLabels, fLabels, 3*sizeof(Int_t)); 
  memcpy(target.fRefCov, fRefCov, 3*sizeof(Double_t)); 
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
  UpDate(track);
  return kTRUE;
}


//_____________________________________________________________________________
void AliTRDseedV1::Reset()
{
  //
  // Reset seed
  //
  fExB=0.;fVD=0.;fT0=0.;fS2PRF=0.;
  fDiffL=0.;fDiffT=0.;
  fClusterIdx=0;
  fN=0;
  fDet=-1;fTilt=0.;fPadLength=0.;
  fMom=0.;
  fdX=0.;fX0=0.; fX=0.; fY=0.; fZ=0.;
  fS2Y=0.; fS2Z=0.;
  fC=0.; fChi2 = 0.;

  for(Int_t ic=kNTimeBins; ic--;) fIndexes[ic] = -1;
  memset(fClusters, 0, kNTimeBins*sizeof(AliTRDcluster*));
  fYref[0] = 0.; fYref[1] = 0.; 
  fZref[0] = 0.; fZref[1] = 0.; 
  fYfit[0] = 0.; fYfit[1] = 0.; 
  fZfit[0] = 0.; fZfit[1] = 0.; 
  memset(fdEdx, 0, kNSlices*sizeof(Float_t)); 
  for(int ispec=0; ispec<AliPID::kSPECIES; ispec++) fProb[ispec]  = -1.;
  fLabels[0]=-1; fLabels[1]=-1; // most freq MC labels
  fLabels[2]=0;  // number of different labels for tracklet
  fRefCov[0] = 1.; fRefCov[1] = 0.; fRefCov[2] = 1.;
  // covariance matrix [diagonal]
  // default sy = 200um and sz = 2.3 cm 
  fCov[0] = 4.e-4; fCov[1] = 0.; fCov[2] = 5.3; 
}

//____________________________________________________________________
void AliTRDseedV1::UpDate(const AliTRDtrackV1 *trk)
{ 
  // update tracklet reference position from the TRD track
  // Funny name to avoid the clash with the function AliTRDseed::Update() (has to be made obselete)

  Double_t fSnp = trk->GetSnp();
  Double_t fTgl = trk->GetTgl();
  fMom = trk->GetP();
  fYref[1] = fSnp/(1. - fSnp*fSnp);
  fZref[1] = fTgl;
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
  for (Int_t i = kNTimeBins; i--; ) {
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
  for (Int_t ic=kNTimeBins; ic--; c++) {
    if(!(*c)) continue;
    if(IsStandAlone()){
      if((*c)->IsShared() || (*c)->IsUsed()){ 
        (*c) = 0x0;
        fIndexes[ic] = -1;
        SetN(GetN()-1);
        if((*c)->IsShared()) SetNShared(GetNShared()-1);
        else SetNUsed(GetNUsed()-1);
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

  Int_t nclusters[kNSlices]; 
  memset(nclusters, 0, kNSlices*sizeof(Int_t));
  memset(fdEdx, 0, kNSlices*sizeof(Float_t));

  const Double_t kDriftLength = (.5 * AliTRDgeometry::AmThick() + AliTRDgeometry::DrThick());

  AliTRDcluster *c = 0x0;
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
    nclusters[slice]++;
  } // End of loop over clusters

  //if(fReconstructor->GetPIDMethod() == AliTRDReconstructor::kLQPID){
  if(nslices == AliTRDpidUtil::kLQslices){
  // calculate mean charge per slice (only LQ PID)
    for(int is=0; is<nslices; is++){ 
      if(nclusters[is]) fdEdx[is] /= nclusters[is];
    }
  }
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
  for (Int_t i = 0; i < kNTimeBins; i++) {
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


//____________________________________________________________________
void AliTRDseedV1::GetClusterXY(const AliTRDcluster *c, Double_t &x, Double_t &y)
{
// Return corrected position of the cluster taking into 
// account variation of the drift velocity with drift length.


  // drift velocity correction TODO to be moved to the clusterizer
  const Float_t cx[] = {
    -9.6280e-02, 1.3091e-01,-1.7415e-02,-9.9221e-02,-1.2040e-01,-9.5493e-02,
    -5.0041e-02,-1.6726e-02, 3.5756e-03, 1.8611e-02, 2.6378e-02, 3.3823e-02,
     3.4811e-02, 3.5282e-02, 3.5386e-02, 3.6047e-02, 3.5201e-02, 3.4384e-02,
     3.2864e-02, 3.1932e-02, 3.2051e-02, 2.2539e-02,-2.5154e-02,-1.2050e-01,
    -1.2050e-01
  };

  // PRF correction TODO to be replaced by the gaussian 
  // approximation with full error parametrization and // moved to the clusterizer
  const Float_t cy[AliTRDgeometry::kNlayer][3] = {
    { 4.014e-04, 8.605e-03, -6.880e+00},
    {-3.061e-04, 9.663e-03, -6.789e+00},
    { 1.124e-03, 1.105e-02, -6.825e+00},
    {-1.527e-03, 1.231e-02, -6.777e+00},
    { 2.150e-03, 1.387e-02, -6.783e+00},
    {-1.296e-03, 1.486e-02, -6.825e+00}
  }; 

  Int_t ily = AliTRDgeometry::GetLayer(c->GetDetector());
  x = c->GetX() - cx[c->GetLocalTimeBin()];
  y = c->GetY() + cy[ily][0] + cy[ily][1] * TMath::Sin(cy[ily][2] * c->GetCenter());
  return;
}

//____________________________________________________________________
Float_t AliTRDseedV1::GetdQdl(Int_t ic) const
{
// Using the linear approximation of the track inside one TRD chamber (TRD tracklet) 
// the charge per unit length can be written as:
// BEGIN_LATEX
// #frac{dq}{dl} = #frac{q_{c}}{dx * #sqrt{1 + #(){#frac{dy}{dx}}^{2}_{fit} + #(){#frac{dy}{dx}}^{2}_{ref}}}
// END_LATEX
// where qc is the total charge collected in the current time bin and dx is the length 
// of the time bin. For the moment (Jan 20 2009) only pad row cross corrections are 
// considered for the charge but none are applied for drift velocity variations along 
// the drift region or assymetry of the TRF
// 
// Author : Alex Bercuci <A.Bercuci@gsi.de>
//
  Float_t dq = 0.;
  if(fClusters[ic]) dq += TMath::Abs(fClusters[ic]->GetQ());
  if(fClusters[ic+kNtb]) dq += TMath::Abs(fClusters[ic+kNtb]->GetQ());
  if(dq<1.e-3 || fdX < 1.e-3) return 0.;

  return dq/fdX/TMath::Sqrt(1. + fYfit[1]*fYfit[1] + fZref[1]*fZref[1]);
}

//____________________________________________________________________
Float_t* AliTRDseedV1::GetProbability(Bool_t force)
{	
  if(!force) return &fProb[0];
  if(!CookPID()) return 0x0;
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
//   returns pointer to the probability array and 0x0 if missing DB access 
//
// Detailed description

  
  // retrive calibration db
  AliTRDcalibDB *calibration = AliTRDcalibDB::Instance();
  if (!calibration) {
    AliError("No access to calibration data");
    return kFALSE;
  }

  if (!fReconstructor) {
    AliError("Reconstructor not set.");
    return kFALSE;
  }

  // Retrieve the CDB container class with the parametric detector response
  const AliTRDCalPID *pd = calibration->GetPIDObject(fReconstructor->GetPIDMethod());
  if (!pd) {
    AliError("No access to AliTRDCalPID object");
    return kFALSE;
  }
  //AliInfo(Form("Method[%d] : %s", fReconstructor->GetRecoParam() ->GetPIDMethod(), pd->IsA()->GetName()));

  // calculate tracklet length TO DO
  Float_t length = (AliTRDgeometry::AmThick() + AliTRDgeometry::DrThick());
  /// TMath::Sqrt((1.0 - fSnp[iPlane]*fSnp[iPlane]) / (1.0 + fTgl[iPlane]*fTgl[iPlane]));
  
  //calculate dE/dx
  CookdEdx(fReconstructor->GetNdEdxSlices());
  
  // Sets the a priori probabilities
  for(int ispec=0; ispec<AliPID::kSPECIES; ispec++) {
    fProb[ispec] = pd->GetProbability(ispec, fMom, &fdEdx[0], length, GetPlane());	
  }

  return kTRUE;
}

//____________________________________________________________________
Float_t AliTRDseedV1::GetQuality(Bool_t kZcorr) const
{
  //
  // Returns a quality measurement of the current seed
  //

  Float_t zcorr = kZcorr ? fTilt * (fZfit[0] - fZref[0]) : 0.;
  return 
      .5 * TMath::Abs(18.0 - GetN())
    + 10.* TMath::Abs(fYfit[1] - fYref[1])
    + 5. * TMath::Abs(fYfit[0] - fYref[0] + zcorr)
    + 2. * TMath::Abs(fZfit[0] - fZref[0]) / fPadLength;
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
  Double_t sz2    = fPadLength*fPadLength/12.;

  // insert systematic uncertainties
  Double_t sys[15];
  fReconstructor->GetRecoParam()->GetSysCovMatrix(sys);
  sy2 += sys[0];
  sz2 += sys[1];

  // rotate covariance matrix
  Double_t t2 = fTilt*fTilt;
  Double_t correction = 1./(1. + t2);
  cov[0] = (sy2+t2*sz2)*correction;
  cov[1] = fTilt*(sz2 - sy2)*correction;
  cov[2] = (t2*sy2+sz2)*correction;
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
    while (ic<kNTimeBins && !(*c)){ic++; c++;} 
    if(*c){
      col = (*c)->GetPadCol();
      row = (*c)->GetPadRow();
    }
  }

  fT0    = t0Det->GetValue(fDet) + t0ROC->GetValue(col,row);
  fVD    = vdDet->GetValue(fDet) * vdROC->GetValue(col, row);
  fS2PRF = calib->GetPRFWidth(fDet, col, row); fS2PRF *= fS2PRF;
  fExB   = AliTRDCommonParam::Instance()->GetOmegaTau(fVD);
  AliTRDCommonParam::Instance()->GetDiffCoeff(fDiffL,
  fDiffT, fVD);
  SetBit(kCalib, kTRUE);
}

//____________________________________________________________________
void AliTRDseedV1::SetOwner()
{
  //AliInfo(Form("own [%s] fOwner[%s]", own?"YES":"NO", fOwner?"YES":"NO"));
  
  if(TestBit(kOwner)) return;
  for(int ic=0; ic<kNTimeBins; ic++){
    if(!fClusters[ic]) continue;
    fClusters[ic] = new AliTRDcluster(*fClusters[ic]);
  }
  SetBit(kOwner);
}

// //____________________________________________________________________
// Bool_t	AliTRDseedV1::AttachClustersIter(AliTRDtrackingChamber *chamber, Float_t quality, Bool_t kZcorr, AliTRDcluster *c)
// {
//   //
//   // Iterative process to register clusters to the seed.
//   // In iteration 0 we try only one pad-row and if quality not
//   // sufficient we try 2 pad-rows (about 5% of tracks cross 2 pad-rows)
//   //
//   // debug level 7
//   //
//   
//   if(!fReconstructor->GetRecoParam() ){
//     AliError("Seed can not be used without a valid RecoParam.");
//     return kFALSE;
//   }
// 
//   AliTRDchamberTimeBin *layer = 0x0;
//   if(fReconstructor->GetStreamLevel(AliTRDReconstructor::kTracker)>=7){
//     AliTRDtrackingChamber ch(*chamber);
//     ch.SetOwner(); 
//     TTreeSRedirector &cstreamer = *fReconstructor->GetDebugStream(AliTRDReconstructor::kTracker);
//     cstreamer << "AttachClustersIter"
//       << "chamber.="   << &ch
//       << "tracklet.="  << this
//       << "\n";	
//   }
// 
//   Float_t  tquality;
//   Double_t kroady = fReconstructor->GetRecoParam() ->GetRoad1y();
//   Double_t kroadz = fPadLength * .5 + 1.;
//   
//   // initialize configuration parameters
//   Float_t zcorr = kZcorr ? fTilt * (fZfit[0] - fZref[0]) : 0.;
//   Int_t   niter = kZcorr ? 1 : 2;
//   
//   Double_t yexp, zexp;
//   Int_t ncl = 0;
//   // start seed update
//   for (Int_t iter = 0; iter < niter; iter++) {
//     ncl = 0;
//     for (Int_t iTime = 0; iTime < AliTRDtrackerV1::GetNTimeBins(); iTime++) {
//       if(!(layer = chamber->GetTB(iTime))) continue;
//       if(!Int_t(*layer)) continue;
//       
//       // define searching configuration
//       Double_t dxlayer = layer->GetX() - fX0;
//       if(c){
//         zexp = c->GetZ();
//         //Try 2 pad-rows in second iteration
//         if (iter > 0) {
//           zexp = fZref[0] + fZref[1] * dxlayer - zcorr;
//           if (zexp > c->GetZ()) zexp = c->GetZ() + fPadLength*0.5;
//           if (zexp < c->GetZ()) zexp = c->GetZ() - fPadLength*0.5;
//         }
//       } else zexp = fZref[0] + (kZcorr ? fZref[1] * dxlayer : 0.);
//       yexp  = fYref[0] + fYref[1] * dxlayer - zcorr;
//       
//       // Get and register cluster
//       Int_t    index = layer->SearchNearestCluster(yexp, zexp, kroady, kroadz);
//       if (index < 0) continue;
//       AliTRDcluster *cl = (*layer)[index];
//       
//       fIndexes[iTime]  = layer->GetGlobalIndex(index);
//       fClusters[iTime] = cl;
// //       fY[iTime]        = cl->GetY();
// //       fZ[iTime]        = cl->GetZ();
//       ncl++;
//     }
//     if(fReconstructor->GetStreamLevel(AliTRDReconstructor::kTracker)>=7) AliInfo(Form("iter = %d ncl [%d] = %d", iter, fDet, ncl));
//     
//     if(ncl>1){	
//       // calculate length of the time bin (calibration aware)
//       Int_t irp = 0; Float_t x[2]={0., 0.}; Int_t tb[2] = {0,0};
//       for (Int_t iTime = 0; iTime < AliTRDtrackerV1::GetNTimeBins(); iTime++) {
//         if(!fClusters[iTime]) continue;
//         x[irp]  = fClusters[iTime]->GetX();
//         tb[irp] = iTime;
//         irp++;
//         if(irp==2) break;
//       } 
//       Int_t dtb = tb[1] - tb[0];
//       fdX = dtb ? (x[0] - x[1]) / dtb : 0.15;
// 
//       // update X0 from the clusters (calibration/alignment aware)
//       for (Int_t iTime = 0; iTime < AliTRDtrackerV1::GetNTimeBins(); iTime++) {
//         if(!(layer = chamber->GetTB(iTime))) continue;
//         if(!layer->IsT0()) continue;
//         if(fClusters[iTime]){ 
//           fX0 = fClusters[iTime]->GetX();
//           break;
//         } else { // we have to infere the position of the anode wire from the other clusters
//           for (Int_t jTime = iTime+1; jTime < AliTRDtrackerV1::GetNTimeBins(); jTime++) {
//             if(!fClusters[jTime]) continue;
//             fX0 = fClusters[jTime]->GetX() + fdX * (jTime - iTime);
//             break;
//           }
//         }
//       }	
//       
//       // update YZ reference point
//       // TODO
//       
//       // update x reference positions (calibration/alignment aware)
// //       for (Int_t iTime = 0; iTime < AliTRDtrackerV1::GetNTimeBins(); iTime++) {
// //         if(!fClusters[iTime]) continue;
// //         fX[iTime] = fX0 - fClusters[iTime]->GetX();
// //       } 
//       
//       FitMI();
//     }
//     if(fReconstructor->GetStreamLevel(AliTRDReconstructor::kTracker)>=7) AliInfo(Form("iter = %d nclFit [%d] = %d", iter, fDet, fN2));
//     
//     if(IsOK()){
//       tquality = GetQuality(kZcorr);
//       if(tquality < quality) break;
//       else quality = tquality;
//     }
//     kroadz *= 2.;
//   } // Loop: iter
//   if (!IsOK()) return kFALSE;
// 
//   if(fReconstructor->GetStreamLevel(AliTRDReconstructor::kTracker)>=1) CookLabels();
// 
//   // load calibration params
//   Calibrate();
//   UpdateUsed();
//   return kTRUE;	
// }

//____________________________________________________________________
Bool_t	AliTRDseedV1::AttachClusters(AliTRDtrackingChamber *chamber, Bool_t tilt)
{
  //
  // Projective algorithm to attach clusters to seeding tracklets
  //
  // Parameters
  //
  // Output
  //
  // Detailed description
  // 1. Collapse x coordinate for the full detector plane
  // 2. truncated mean on y (r-phi) direction
  // 3. purge clusters
  // 4. truncated mean on z direction
  // 5. purge clusters
  // 6. fit tracklet
  //	
  Bool_t kPRINT = kFALSE;
  if(!fReconstructor->GetRecoParam() ){
    AliError("Seed can not be used without a valid RecoParam.");
    return kFALSE;
  }
  // Initialize reco params for this tracklet
  // 1. first time bin in the drift region
  Int_t t0 = 4;
  Int_t kClmin = Int_t(fReconstructor->GetRecoParam() ->GetFindableClusters()*AliTRDtrackerV1::GetNTimeBins());

  Double_t syRef  = TMath::Sqrt(fRefCov[0]);
  //define roads
  Double_t kroady = 1.; 
  //fReconstructor->GetRecoParam() ->GetRoad1y();
  Double_t kroadz = fPadLength * 1.5 + 1.;
  if(kPRINT) printf("AttachClusters() sy[%f] road[%f]\n", syRef, kroady);

  // working variables
  const Int_t kNrows = 16;
  AliTRDcluster *clst[kNrows][kNTimeBins];
  Double_t cond[4], dx, dy, yt, zt,
    yres[kNrows][kNTimeBins];
  Int_t idxs[kNrows][kNTimeBins], ncl[kNrows], ncls = 0;
  memset(ncl, 0, kNrows*sizeof(Int_t));
  memset(clst, 0, kNrows*kNTimeBins*sizeof(AliTRDcluster*));

  // Do cluster projection
  AliTRDcluster *c = 0x0;
  AliTRDchamberTimeBin *layer = 0x0;
  Bool_t kBUFFER = kFALSE;
  for (Int_t it = 0; it < AliTRDtrackerV1::GetNTimeBins(); it++) {
    if(!(layer = chamber->GetTB(it))) continue;
    if(!Int_t(*layer)) continue;
    
    dx   = fX0 - layer->GetX();
    yt = fYref[0] - fYref[1] * dx;
    zt = fZref[0] - fZref[1] * dx;
    if(kPRINT) printf("\t%2d dx[%f] yt[%f] zt[%f]\n", it, dx, yt, zt);

    // select clusters on a 5 sigmaKalman level
    cond[0] = yt; cond[2] = kroady;
    cond[1] = zt; cond[3] = kroadz;
    Int_t n=0, idx[6];
    layer->GetClusters(cond, idx, n, 6);
    for(Int_t ic = n; ic--;){
      c  = (*layer)[idx[ic]];
      dy = yt - c->GetY();
      dy += tilt ? fTilt * (c->GetZ() - zt) : 0.;
      // select clusters on a 3 sigmaKalman level
/*      if(tilt && TMath::Abs(dy) > 3.*syRef){ 
        printf("too large !!!\n");
        continue;
      }*/
      Int_t r = c->GetPadRow();
      if(kPRINT) printf("\t\t%d dy[%f] yc[%f] r[%d]\n", ic, TMath::Abs(dy), c->GetY(), r);
      clst[r][ncl[r]] = c;
      idxs[r][ncl[r]] = idx[ic];
      yres[r][ncl[r]] = dy;
      ncl[r]++; ncls++;

      if(ncl[r] >= kNTimeBins) {
        AliWarning(Form("Cluster candidates reached limit %d. Some may be lost.", kNTimeBins));
        kBUFFER = kTRUE;
        break;
      }
    }
    if(kBUFFER) break;
  }
  if(kPRINT) printf("Found %d clusters\n", ncls);
  if(ncls<kClmin) return kFALSE;
 
  // analyze each row individualy
  Double_t mean, syDis;
  Int_t nrow[] = {0, 0, 0}, nr = 0, lr=-1;
  for(Int_t ir=kNrows; ir--;){
    if(!(ncl[ir])) continue;
    if(lr>0 && lr-ir != 1){
      if(kPRINT) printf("W - gap in rows attached !!\n"); 
    }
    if(kPRINT) printf("\tir[%d] lr[%d] n[%d]\n", ir, lr, ncl[ir]);
    // Evaluate truncated mean on the y direction
    if(ncl[ir] > 3) AliMathBase::EvaluateUni(ncl[ir], yres[ir], mean, syDis, Int_t(ncl[ir]*.8));
    else {
      mean = 0.; syDis = 0.;
    } 

    // TODO check mean and sigma agains cluster resolution !!
    if(kPRINT) printf("\tr[%2d] m[%f %5.3fsigma] s[%f]\n", ir, mean, TMath::Abs(mean/syRef), syDis);
    // select clusters on a 3 sigmaDistr level
    Bool_t kFOUND = kFALSE;
    for(Int_t ic = ncl[ir]; ic--;){
      if(yres[ir][ic] - mean > 3. * syDis){ 
        clst[ir][ic] = 0x0; continue;
      }
      nrow[nr]++; kFOUND = kTRUE;
    }
    // exit loop
    if(kFOUND) nr++; 
    lr = ir; if(nr>=3) break;
  }
  if(kPRINT) printf("lr[%d] nr[%d] nrow[0]=%d nrow[1]=%d nrow[2]=%d\n", lr, nr, nrow[0], nrow[1], nrow[2]);

  // classify cluster rows
  Int_t row = -1;
  switch(nr){
  case 1:
    row = lr;
    break;
  case 2:
    SetBit(kRowCross, kTRUE); // mark pad row crossing
    if(nrow[0] > nrow[1]){ row = lr+1; lr = -1;}
    else{ 
      row = lr; lr = 1;
      nrow[2] = nrow[1];
      nrow[1] = nrow[0];
      nrow[0] = nrow[2];
    }
    break;
  case 3:
    SetBit(kRowCross, kTRUE); // mark pad row crossing
    break;
  }
  if(kPRINT) printf("\trow[%d] n[%d]\n\n", row, nrow[0]);
  if(row<0) return kFALSE;

  // Select and store clusters 
  // We should consider here :
  //  1. How far is the chamber boundary
  //  2. How big is the mean
  Int_t n = 0;
  for (Int_t ir = 0; ir < nr; ir++) {
    Int_t jr = row + ir*lr; 
    if(kPRINT) printf("\tattach %d clusters for row %d\n", ncl[jr], jr);
    for (Int_t ic = 0; ic < ncl[jr]; ic++) {
      if(!(c = clst[jr][ic])) continue;
      Int_t it = c->GetPadTime();
      // TODO proper indexing of clusters !!
      fIndexes[it+kNtb*ir]  = chamber->GetTB(it)->GetGlobalIndex(idxs[jr][ic]);
      fClusters[it+kNtb*ir] = c;
  
      //printf("\tid[%2d] it[%d] idx[%d]\n", ic, it, fIndexes[it]);
  
      n++;
    }
  }  

  // number of minimum numbers of clusters expected for the tracklet
  if (n < kClmin){
    //AliWarning(Form("Not enough clusters to fit the tracklet %d [%d].", n, kClmin));
    return kFALSE;
  }
  SetN(n);

  // Load calibration parameters for this tracklet  
  Calibrate();

  // calculate dx for time bins in the drift region (calibration aware)
  Int_t irp = 0; Float_t x[2] = {0.,0.}; Int_t tb[2]={0,0};
  for (Int_t it = t0; it < AliTRDtrackerV1::GetNTimeBins(); it++) {
    if(!fClusters[it]) continue;
    x[irp]  = fClusters[it]->GetX();
    tb[irp] = it;
    irp++;
    if(irp==2) break;
  }  
  Int_t dtb = tb[1] - tb[0];
  fdX = dtb ? (x[0] - x[1]) / dtb : 0.15;

  // update X0 from the clusters (calibration/alignment aware) TODO remove dependence on x0 !!
  for (Int_t it = 0; it < AliTRDtrackerV1::GetNTimeBins(); it++) {
    if(!(layer = chamber->GetTB(it))) continue;
    if(!layer->IsT0()) continue;
    if(fClusters[it]){ 
      fX0 = fClusters[it]->GetX();
      break;
    } else { // we have to infere the position of the anode wire from the other clusters
      for (Int_t jt = it+1; jt < AliTRDtrackerV1::GetNTimeBins(); jt++) {
        if(!fClusters[jt]) continue;
        fX0 = fClusters[jt]->GetX() + fdX * (jt - it);
        break;
      }
    }
  }	

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
  fReconstructor = rec;
  AliTRDgeometry g;
  AliTRDpadPlane *pp = g.GetPadPlane(fDet);
  fTilt      = TMath::Tan(TMath::DegToRad()*pp->GetTiltingAngle());
  fPadLength = pp->GetLengthIPad();
  //fSnp = fYref[1]/TMath::Sqrt(1+fYref[1]*fYref[1]);
  //fTgl = fZref[1];
  Int_t n = 0, nshare = 0, nused = 0;
  AliTRDcluster **cit = &fClusters[0];
  for(Int_t ic = kNTimeBins; ic--; cit++){
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
Bool_t AliTRDseedV1::Fit(Bool_t tilt, Int_t errors)
{
  //
  // Linear fit of the tracklet
  //
  // Parameters :
  //
  // Output :
  //  True if successful
  //
  // Detailed description
  // 2. Check if tracklet crosses pad row boundary
  // 1. Calculate residuals in the y (r-phi) direction
  // 3. Do a Least Square Fit to the data
  //

  if(!IsCalibrated()){
    AliWarning("Tracklet fit failed. Call Calibrate().");
    return kFALSE;
  }

  const Int_t kClmin = 8;


  // cluster error parametrization parameters 
  // 1. sy total charge
  const Float_t sq0inv = 0.019962; // [1/q0]
  const Float_t sqb    = 1.0281564;    //[cm]
  // 2. sy for the PRF
  const Float_t scy[AliTRDgeometry::kNlayer][4] = {
    {2.827e-02, 9.600e-04, 4.296e-01, 2.271e-02},
    {2.952e-02,-2.198e-04, 4.146e-01, 2.339e-02},
    {3.090e-02, 1.514e-03, 4.020e-01, 2.402e-02},
    {3.260e-02,-2.037e-03, 3.946e-01, 2.509e-02},
    {3.439e-02,-3.601e-04, 3.883e-01, 2.623e-02},
    {3.510e-02, 2.066e-03, 3.651e-01, 2.588e-02},
  };
  // 3. sy parallel to the track
  const Float_t sy0 =  2.649e-02; // [cm]
  const Float_t sya = -8.864e-04; // [cm]
  const Float_t syb = -2.435e-01; // [cm]

  // 4. sx parallel to the track
  const Float_t sxgc = 5.427e-02;
  const Float_t sxgm = 7.783e-01;
  const Float_t sxgs = 2.743e-01;
  const Float_t sxe0 =-2.065e+00;
  const Float_t sxe1 =-2.978e-02;

  // 5. sx perpendicular to the track
//   const Float_t sxd0 = 1.881e-02;
//   const Float_t sxd1 =-4.101e-01;
//   const Float_t sxd2 = 1.572e+00;

  // get track direction
  Double_t y0   = fYref[0];
  Double_t dydx = fYref[1]; 
  Double_t z0   = fZref[0];
  Double_t dzdx = fZref[1];
  Double_t yt, zt;

  const Int_t kNtb = AliTRDtrackerV1::GetNTimeBins();
  // calculation of tg^2(phi - a_L) and tg^2(a_L)
  Double_t tgg = (dydx-fExB)/(1.+dydx*fExB); tgg *= tgg;
  //Double_t exb2= fExB*fExB;

  //AliTRDtrackerV1::AliTRDLeastSquare fitterZ;
  TLinearFitter  fitterY(1, "pol1");
  // convertion factor from square to gauss distribution for sigma
  //Double_t convert = 1./TMath::Sqrt(12.);
  
  // book cluster information
  Double_t qc[kNTimeBins], xc[kNTimeBins], yc[kNTimeBins], zc[kNTimeBins], sy[kNTimeBins];

  Int_t ily = AliTRDgeometry::GetLayer(fDet);
  Int_t fN = 0;
  AliTRDcluster *c=0x0, **jc = &fClusters[0];
  for (Int_t ic=0; ic<kNtb; ic++, ++jc) {
    //zRow[ic] = -1;
    xc[ic]  = -1.;
    yc[ic]  = 999.;
    zc[ic]  = 999.;
    sy[ic]  = 0.;
    //sz[ic]  = 0.;
    if(!(c = (*jc))) continue;
    if(!c->IsInChamber()) continue;

    Float_t w = 1.;
    if(c->GetNPads()>4) w = .5;
    if(c->GetNPads()>5) w = .2;

    //zRow[fN] = c->GetPadRow();
    qc[fN]   = TMath::Abs(c->GetQ());
    // correct cluster position for PRF and v drift
    //Int_t jc = TMath::Max(fN-3, 0);
    //xc[fN]   = c->GetXloc(fT0, fVD, &qc[jc], &xc[jc]/*, z0 - c->GetX()*dzdx*/);
    //Double_t s2 = fS2PRF + fDiffL*fDiffL*xc[fN]/(1.+2.*exb2)+tgg*xc[fN]*xc[fN]*exb2/12.;
    //yc[fN]   = c->GetYloc(s2, fPadLength, xc[fN], fExB);
    
    // uncalibrated cluster correction 
    // TODO remove
    Double_t x, y; GetClusterXY(c, x, y);
    xc[fN]   = fX0 - x;
    yc[fN]   = y;
    zc[fN]   = c->GetZ();

    // extrapolated y value for the track
    yt = y0 - xc[fN]*dydx; 
    // extrapolated z value for the track
    zt = z0 - xc[fN]*dzdx; 
    // tilt correction
    if(tilt) yc[fN] -= fTilt*(zc[fN] - zt); 

    // ELABORATE CLUSTER ERROR
    // TODO to be moved to AliTRDcluster
    // basic y error (|| to track).
    sy[fN]  = xc[fN] < AliTRDgeometry::CamHght() ? 2. : sy0 + sya*TMath::Exp(1./(xc[fN]+syb));
    //printf("cluster[%d]\n\tsy[0] = %5.3e [um]\n", fN,  sy[fN]*1.e4);
    // y error due to total charge
    sy[fN] += sqb*(1./qc[fN] - sq0inv);
    //printf("\tsy[1] = %5.3e [um]\n", sy[fN]*1.e4);
    // y error due to PRF
    sy[fN] += scy[ily][0]*TMath::Gaus(c->GetCenter(), scy[ily][1], scy[ily][2]) - scy[ily][3];
    //printf("\tsy[2] = %5.3e [um]\n", sy[fN]*1.e4);

    sy[fN] *= sy[fN];

    // ADD ERROR ON x
    // error of drift length parallel to the track
    Double_t sx = sxgc*TMath::Gaus(xc[fN], sxgm, sxgs) + TMath::Exp(sxe0+sxe1*xc[fN]); // [cm]
    //printf("\tsx[0] = %5.3e [um]\n", sx*1.e4);
    // error of drift length perpendicular to the track
    //sx += sxd0 + sxd1*d + sxd2*d*d;
    sx *= sx; // square sx

    // add error from ExB 
    if(errors>0) sy[fN] += fExB*fExB*sx;
    //printf("\tsy[3] = %5.3e [um^2]\n", sy[fN]*1.e8);

    // global radial error due to misalignment/miscalibration
    Double_t sx0  = 0.; sx0 *= sx0;
    // add sx contribution to sy due to track angle
    if(errors>1) sy[fN] += tgg*(sx+sx0);
    // TODO we should add tilt pad correction here
    //printf("\tsy[4] = %5.3e [um^2]\n", sy[fN]*1.e8);
    c->SetSigmaY2(sy[fN]);

    sy[fN]  = TMath::Sqrt(sy[fN]);
    fitterY.AddPoint(&xc[fN], yc[fN], sy[fN]);
    fN++;
  }
  // to few clusters
  if (fN < kClmin) return kFALSE; 

  // fit XY
  fitterY.Eval();
  fYfit[0] = fitterY.GetParameter(0);
  fYfit[1] = -fitterY.GetParameter(1);
  // store covariance
  Double_t *p = fitterY.GetCovarianceMatrix();
  fCov[0] = p[0]; // variance of y0
  fCov[1] = p[1]; // covariance of y0, dydx
  fCov[2] = p[3]; // variance of dydx
  // the ref radial position is set at the minimum of 
  // the y variance of the tracklet
  fX   = -fCov[1]/fCov[2]; //fXref = fX0 - fXref;
  fS2Y = fCov[0] +2.*fX*fCov[1] + fX*fX*fCov[2];

  // fit XZ
  if(IsRowCross()){ 
    // TODO pad row cross position estimation !!!
    //AliInfo(Form("Padrow cross in detector %d", fDet));
    fZfit[0] = .5*(zc[0]+zc[fN-1]); fZfit[1] = 0.;
    fS2Z     = 0.02+1.55*fZref[1]; fS2Z *= fS2Z;
  } else {
    fZfit[0] = zc[0]; fZfit[1] = 0.;
    fS2Z     = fPadLength*fPadLength/12.;
  }


//   // determine z offset of the fit
//   Float_t zslope = 0.;
//   Int_t nchanges = 0, nCross = 0;
//   if(nz==2){ // tracklet is crossing pad row
//     // Find the break time allowing one chage on pad-rows
//     // with maximal number of accepted clusters
//     Int_t padRef = zRow[0];
//     for (Int_t ic=1; ic<fN; ic++) {
//       if(zRow[ic] == padRef) continue;
//       
//       // debug
//       if(zRow[ic-1] == zRow[ic]){
//         printf("ERROR in pad row change!!!\n");
//       }
//     
//       // evaluate parameters of the crossing point
//       Float_t sx = (xc[ic-1] - xc[ic])*convert;
//       fCross[0] = .5 * (xc[ic-1] + xc[ic]);
//       fCross[2] = .5 * (zc[ic-1] + zc[ic]);
//       fCross[3] = TMath::Max(dzdx * sx, .01);
//       zslope    = zc[ic-1] > zc[ic] ? 1. : -1.;
//       padRef    = zRow[ic];
//       nCross    = ic;
//       nchanges++;
//     }
//   }
// 
//   // condition on nCross and reset nchanges TODO
// 
//   if(nchanges==1){
//     if(dzdx * zslope < 0.){
//       AliInfo("Tracklet-Track mismatch in dzdx. TODO.");
//     }
// 
// 
//     //zc[nc] = fitterZ.GetFunctionParameter(0); 
//     fCross[1] = fYfit[0] - fCross[0] * fYfit[1];
//     fCross[0] = fX0 - fCross[0];
//   }

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
  Float_t  ycrosscor = fPadLength * fTilt * 0.5;           // Y correction for crossing 
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
  //Float_t  anglecor = fTilt * fZref[1];  // Correction to the angle
  
  Float_t pos[3*kNtb]; memset(pos, 0, 3*kNtb*sizeof(Float_t));
  Float_t *fX = &pos[0], *fY = &pos[kNtb], *fZ = &pos[2*kNtb];
  
  Int_t fN  = 0; AliTRDcluster *c = 0x0; 
  fN2 = 0;
  for (Int_t i = 0; i < AliTRDtrackerV1::GetNTimeBins(); i++) {
    yres[i] = 10000.0;
    if (!(c = fClusters[i])) continue;
    if(!c->IsInChamber()) continue;
    // Residual y
    //yres[i] = fY[i] - fYref[0] - (fYref[1] + anglecor) * fX[i] + fTilt*(fZ[i] - fZref[0]);
    fX[i] = fX0 - c->GetX();
    fY[i] = c->GetY();
    fZ[i] = c->GetZ();
    yres[i] = fY[i] - fTilt*(fZ[i] - (fZref[0] - fX[i]*fZref[1]));
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
      //yres[i] = fY[i] - fYref[0] - (fYref[1] + anglecor) * fX[i] + fTilt*(fZ[i] - fZref[0]);   
      yres[i] = fY[i] - fTilt*(fZ[i] - (fZref[0] - fX[i]*fZref[1]));
//       if (TMath::Abs(fZ[i] - fZProb) > 2) {
//         if (fZ[i] > fZProb) yres[i] += fTilt * fPadLength;
//         if (fZ[i] < fZProb) yres[i] -= fTilt * fPadLength;
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

  AliInfo(Form("Det[%3d] Tilt[%+6.2f] Pad[%5.2f]", fDet, fTilt, fPadLength));
  AliInfo(Form("N[%2d] Nused[%2d] Nshared[%2d]", GetN(), GetNUsed(), GetNShared()));
  AliInfo(Form("x[%7.2f] y[%7.2f] z[%7.2f] dydx[%5.2f] dzdx[%5.2f]", fX0, fYfit[0], fZfit[0], fYfit[1], fZfit[1]));
  AliInfo(Form("Ref        y[%7.2f] z[%7.2f] dydx[%5.2f] dzdx[%5.2f]", fYref[0], fZref[0], fYref[1], fZref[1]))


  if(strcmp(o, "a")!=0) return;

  AliTRDcluster* const* jc = &fClusters[0];
  for(int ic=0; ic<kNTimeBins; ic++, jc++) {
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
  if ( fTilt != inTracklet->fTilt ) return kFALSE;
  if ( fPadLength != inTracklet->fPadLength ) return kFALSE;
  
  for (Int_t i = 0; i < kNTimeBins; i++){
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
  if ( fMom != inTracklet->fMom ) return kFALSE;
  if ( fdX != inTracklet->fdX ) return kFALSE;
  
  for (Int_t iCluster = 0; iCluster < kNTimeBins; iCluster++){
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
