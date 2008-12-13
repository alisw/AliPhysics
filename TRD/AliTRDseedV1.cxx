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
#include "Cal/AliTRDCalPID.h"

ClassImp(AliTRDseedV1)

//____________________________________________________________________
AliTRDseedV1::AliTRDseedV1(Int_t det) 
  :AliTRDseed()
  ,fReconstructor(0x0)
  ,fClusterIter(0x0)
  ,fClusterIdx(0)
  ,fDet(det)
  ,fMom(0.)
  ,fSnp(0.)
  ,fTgl(0.)
  ,fdX(0.)
  ,fXref(0.)
{
  //
  // Constructor
  //
  //printf("AliTRDseedV1::AliTRDseedV1()\n");

  for(int islice=0; islice < knSlices; islice++) fdEdx[islice] = 0.;
  for(int ispec=0; ispec<AliPID::kSPECIES; ispec++) fProb[ispec]  = -1.;
  fRefCov[0] = 1.; fRefCov[1] = 0.; fRefCov[2] = 1.;
}

//____________________________________________________________________
AliTRDseedV1::AliTRDseedV1(const AliTRDseedV1 &ref)
  :AliTRDseed((AliTRDseed&)ref)
  ,fReconstructor(ref.fReconstructor)
  ,fClusterIter(0x0)
  ,fClusterIdx(0)
  ,fDet(ref.fDet)
  ,fMom(ref.fMom)
  ,fSnp(ref.fSnp)
  ,fTgl(ref.fTgl)
  ,fdX(ref.fdX)
  ,fXref(ref.fXref)
{
  //
  // Copy Constructor performing a deep copy
  //

  //printf("AliTRDseedV1::AliTRDseedV1(const AliTRDseedV1 &)\n");
  SetBit(kOwner, kFALSE);
  for(int islice=0; islice < knSlices; islice++) fdEdx[islice] = ref.fdEdx[islice];
  for(int ispec=0; ispec<AliPID::kSPECIES; ispec++) fProb[ispec] = ref.fProb[ispec];
  memcpy(fRefCov, ref.fRefCov, 3*sizeof(Double_t));
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

  if(IsOwner()) 
    for(int itb=0; itb<knTimebins; itb++){
      if(!fClusters[itb]) continue; 
      //AliInfo(Form("deleting c %p @ %d", fClusters[itb], itb));
      delete fClusters[itb];
      fClusters[itb] = 0x0;
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

  target.fClusterIter   = 0x0;
  target.fClusterIdx    = 0;
  target.fDet           = fDet;
  target.fMom           = fMom;
  target.fSnp           = fSnp;
  target.fTgl           = fTgl;
  target.fdX            = fdX;
  target.fXref          = fXref;
  target.fReconstructor = fReconstructor;
  
  for(int islice=0; islice < knSlices; islice++) target.fdEdx[islice] = fdEdx[islice];
  for(int ispec=0; ispec<AliPID::kSPECIES; ispec++) target.fProb[ispec] = fProb[ispec];
  memcpy(target.fRefCov, fRefCov, 3*sizeof(Double_t));
  
  AliTRDseed::Copy(target);
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
  fYref[0] = y;
  fYref[1] = track->GetSnp()/(1. - track->GetSnp()*track->GetSnp());
  fZref[0] = z;
  fZref[1] = track->GetTgl();
  
  const Double_t *cov = track->GetCovariance();
  fRefCov[0] = cov[0]; // Var(y)
  fRefCov[1] = cov[1]; // Cov(yz)
  fRefCov[2] = cov[5]; // Var(z)

  //printf("Tracklet ref x[%7.3f] y[%7.3f] z[%7.3f], snp[%f] tgl[%f]\n", fX0, fYref[0], fZref[0], track->GetSnp(), track->GetTgl());
  return kTRUE;
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
// The calculation of dQ/dl are done using the tracklet fit results (see AliTRDseedV1::GetdQdl(Int_t)) i.e.
//
// dQ/dl = qc/(dx * sqrt(1 + dy/dx^2 + dz/dx^2))
//
// The following effects are included in the calculation:
// 1. calibration values for t0 and vdrift (using x coordinate to calculate slice)
// 2. cluster sharing (optional see AliTRDrecoParam::SetClusterSharing())
// 3. cluster size
//

  Int_t nclusters[knSlices];
  for(int i=0; i<knSlices; i++){ 
    fdEdx[i]     = 0.;
    nclusters[i] = 0;
  }
  Float_t clength = (/*.5 * */AliTRDgeometry::AmThick() + AliTRDgeometry::DrThick());

  AliTRDcluster *cluster = 0x0;
  for(int ic=0; ic<AliTRDtrackerV1::GetNTimeBins(); ic++){
    if(!(cluster = fClusters[ic])) continue;
    Float_t x = cluster->GetX();
    
    // Filter clusters for dE/dx calculation
    
    // 1.consider calibration effects for slice determination
    Int_t slice; 
    if(cluster->IsInChamber()) slice = Int_t(TMath::Abs(fX0 - x) * nslices / clength);
    else slice = x < fX0 ? 0 : nslices-1;
    
    // 2. take sharing into account
    Float_t w = cluster->IsShared() ? .5 : 1.;
    
    // 3. take into account large clusters TODO
    //w *= c->GetNPads() > 3 ? .8 : 1.;
    
    //CHECK !!!
    fdEdx[slice]   += w * GetdQdl(ic); //fdQdl[ic];
    nclusters[slice]++;
  } // End of loop over clusters

  //if(fReconstructor->GetPIDMethod() == AliTRDReconstructor::kLQPID){
  if(nslices == AliTRDReconstructor::kLQslices){
  // calculate mean charge per slice (only LQ PID)
    for(int is=0; is<nslices; is++){ 
      if(nclusters[is]) fdEdx[is] /= nclusters[is];
    }
  }
}


//____________________________________________________________________
Float_t AliTRDseedV1::GetdQdl(Int_t ic) const
{
  return fClusters[ic] ? TMath::Abs(fClusters[ic]->GetQ()) /fdX / TMath::Sqrt(1. + fYfit[1]*fYfit[1] + fZref[1]*fZref[1]) : 0.;
}

//____________________________________________________________________
Double_t* AliTRDseedV1::GetProbability()
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
    return 0x0;
  }

  if (!fReconstructor) {
    AliError("Reconstructor not set.");
    return 0x0;
  }

  // Retrieve the CDB container class with the parametric detector response
  const AliTRDCalPID *pd = calibration->GetPIDObject(fReconstructor->GetPIDMethod());
  if (!pd) {
    AliError("No access to AliTRDCalPID object");
    return 0x0;
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

  return &fProb[0];
}

//____________________________________________________________________
Float_t AliTRDseedV1::GetQuality(Bool_t kZcorr) const
{
  //
  // Returns a quality measurement of the current seed
  //

  Float_t zcorr = kZcorr ? fTilt * (fZProb - fZref[0]) : 0.;
  return 
      .5 * TMath::Abs(18.0 - fN2)
    + 10.* TMath::Abs(fYfit[1] - fYref[1])
    + 5. * TMath::Abs(fYfit[0] - fYref[0] + zcorr)
    + 2. * TMath::Abs(fMeanz - fZref[0]) / fPadLength;
}

//____________________________________________________________________
void AliTRDseedV1::GetCovAt(Double_t /*x*/, Double_t *cov) const
{
// Computes covariance in the y-z plane at radial point x

  Int_t ic = 0; while (!fClusters[ic]) ic++; 
  AliTRDcalibDB *fCalib = AliTRDcalibDB::Instance();
  Double_t exB         = fCalib->GetOmegaTau(fCalib->GetVdriftAverage(fClusters[ic]->GetDetector()), -AliTracker::GetBz()*0.1);

  Double_t sy2    = fSigmaY2*fSigmaY2 + .2*(fYfit[1]-exB)*(fYfit[1]-exB);
  Double_t sz2    = fPadLength/12.;


  //printf("Yfit[1] %f sy20 %f SigmaY2 %f\n", fYfit[1], sy20, fSigmaY2);

  cov[0] = sy2;
  cov[1] = fTilt*(sy2-sz2);
  cov[2] = sz2;

  // insert systematic uncertainties calibration and misalignment
  Double_t sys[15];
  fReconstructor->GetRecoParam()->GetSysCovMatrix(sys);
  cov[0] += (sys[0]*sys[0]);
  cov[2] += (sys[1]*sys[1]);
}


//____________________________________________________________________
void AliTRDseedV1::SetOwner()
{
  //AliInfo(Form("own [%s] fOwner[%s]", own?"YES":"NO", fOwner?"YES":"NO"));
  
  if(TestBit(kOwner)) return;
  for(int ic=0; ic<knTimebins; ic++){
    if(!fClusters[ic]) continue;
    fClusters[ic] = new AliTRDcluster(*fClusters[ic]);
  }
  SetBit(kOwner);
}

//____________________________________________________________________
Bool_t	AliTRDseedV1::AttachClustersIter(AliTRDtrackingChamber *chamber, Float_t quality, Bool_t kZcorr, AliTRDcluster *c)
{
  //
  // Iterative process to register clusters to the seed.
  // In iteration 0 we try only one pad-row and if quality not
  // sufficient we try 2 pad-rows (about 5% of tracks cross 2 pad-rows)
  //
  // debug level 7
  //
  
  if(!fReconstructor->GetRecoParam() ){
    AliError("Seed can not be used without a valid RecoParam.");
    return kFALSE;
  }

  AliTRDchamberTimeBin *layer = 0x0;
  if(fReconstructor->GetStreamLevel(AliTRDReconstructor::kTracker)>=7){
    AliTRDtrackingChamber ch(*chamber);
    ch.SetOwner(); 
    TTreeSRedirector &cstreamer = *fReconstructor->GetDebugStream(AliTRDReconstructor::kTracker);
    cstreamer << "AttachClustersIter"
      << "chamber.="   << &ch
      << "tracklet.="  << this
      << "\n";	
  }

  Float_t  tquality;
  Double_t kroady = fReconstructor->GetRecoParam() ->GetRoad1y();
  Double_t kroadz = fPadLength * .5 + 1.;
  
  // initialize configuration parameters
  Float_t zcorr = kZcorr ? fTilt * (fZProb - fZref[0]) : 0.;
  Int_t   niter = kZcorr ? 1 : 2;
  
  Double_t yexp, zexp;
  Int_t ncl = 0;
  // start seed update
  for (Int_t iter = 0; iter < niter; iter++) {
    ncl = 0;
    for (Int_t iTime = 0; iTime < AliTRDtrackerV1::GetNTimeBins(); iTime++) {
      if(!(layer = chamber->GetTB(iTime))) continue;
      if(!Int_t(*layer)) continue;
      
      // define searching configuration
      Double_t dxlayer = layer->GetX() - fX0;
      if(c){
        zexp = c->GetZ();
        //Try 2 pad-rows in second iteration
        if (iter > 0) {
          zexp = fZref[0] + fZref[1] * dxlayer - zcorr;
          if (zexp > c->GetZ()) zexp = c->GetZ() + fPadLength*0.5;
          if (zexp < c->GetZ()) zexp = c->GetZ() - fPadLength*0.5;
        }
      } else zexp = fZref[0] + (kZcorr ? fZref[1] * dxlayer : 0.);
      yexp  = fYref[0] + fYref[1] * dxlayer - zcorr;
      
      // Get and register cluster
      Int_t    index = layer->SearchNearestCluster(yexp, zexp, kroady, kroadz);
      if (index < 0) continue;
      AliTRDcluster *cl = (*layer)[index];
      
      fIndexes[iTime]  = layer->GetGlobalIndex(index);
      fClusters[iTime] = cl;
      fY[iTime]        = cl->GetY();
      fZ[iTime]        = cl->GetZ();
      ncl++;
    }
    if(fReconstructor->GetStreamLevel(AliTRDReconstructor::kTracker)>=7) AliInfo(Form("iter = %d ncl [%d] = %d", iter, fDet, ncl));
    
    if(ncl>1){	
      // calculate length of the time bin (calibration aware)
      Int_t irp = 0; Float_t x[2]; Int_t tb[2];
      for (Int_t iTime = 0; iTime < AliTRDtrackerV1::GetNTimeBins(); iTime++) {
        if(!fClusters[iTime]) continue;
        x[irp]  = fClusters[iTime]->GetX();
        tb[irp] = iTime;
        irp++;
        if(irp==2) break;
      } 
      fdX = (x[1] - x[0]) / (tb[0] - tb[1]);
  
      // update X0 from the clusters (calibration/alignment aware)
      for (Int_t iTime = 0; iTime < AliTRDtrackerV1::GetNTimeBins(); iTime++) {
        if(!(layer = chamber->GetTB(iTime))) continue;
        if(!layer->IsT0()) continue;
        if(fClusters[iTime]){ 
          fX0 = fClusters[iTime]->GetX();
          break;
        } else { // we have to infere the position of the anode wire from the other clusters
          for (Int_t jTime = iTime+1; jTime < AliTRDtrackerV1::GetNTimeBins(); jTime++) {
            if(!fClusters[jTime]) continue;
            fX0 = fClusters[jTime]->GetX() + fdX * (jTime - iTime);
            break;
          }
        }
      }	
      
      // update YZ reference point
      // TODO
      
      // update x reference positions (calibration/alignment aware)
      for (Int_t iTime = 0; iTime < AliTRDtrackerV1::GetNTimeBins(); iTime++) {
        if(!fClusters[iTime]) continue;
        fX[iTime] = fClusters[iTime]->GetX() - fX0;
      } 
      
      AliTRDseed::Update();
    }
    if(fReconstructor->GetStreamLevel(AliTRDReconstructor::kTracker)>=7) AliInfo(Form("iter = %d nclFit [%d] = %d", iter, fDet, fN2));
    
    if(IsOK()){
      tquality = GetQuality(kZcorr);
      if(tquality < quality) break;
      else quality = tquality;
    }
    kroadz *= 2.;
  } // Loop: iter
  if (!IsOK()) return kFALSE;

  if(fReconstructor->GetStreamLevel(AliTRDReconstructor::kTracker)>=1) CookLabels();
  UpdateUsed();
  return kTRUE;	
}

//____________________________________________________________________
Bool_t	AliTRDseedV1::AttachClusters(AliTRDtrackingChamber *chamber
                                      ,Bool_t kZcorr)
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

  if(!fReconstructor->GetRecoParam() ){
    AliError("Seed can not be used without a valid RecoParam.");
    return kFALSE;
  }

  const Int_t kClusterCandidates = 2 * knTimebins;
  
  //define roads
  Double_t kroady = fReconstructor->GetRecoParam() ->GetRoad1y();
  Double_t kroadz = fPadLength * 1.5 + 1.;
  // correction to y for the tilting angle
  Float_t zcorr = kZcorr ? fTilt * (fZProb - fZref[0]) : 0.;

  // working variables
  AliTRDcluster *clusters[kClusterCandidates];
  Double_t cond[4], yexp[knTimebins], zexp[knTimebins],
    yres[kClusterCandidates], zres[kClusterCandidates];
  Int_t ncl, *index = 0x0, tboundary[knTimebins];
  
  // Do cluster projection
  AliTRDchamberTimeBin *layer = 0x0;
  Int_t nYclusters = 0; Bool_t kEXIT = kFALSE;
  for (Int_t iTime = 0; iTime < AliTRDtrackerV1::GetNTimeBins(); iTime++) {
    if(!(layer = chamber->GetTB(iTime))) continue;
    if(!Int_t(*layer)) continue;
    
    fX[iTime] = layer->GetX() - fX0;
    zexp[iTime] = fZref[0] + fZref[1] * fX[iTime];
    yexp[iTime] = fYref[0] + fYref[1] * fX[iTime] - zcorr;
    
    // build condition and process clusters
    cond[0] = yexp[iTime] - kroady; cond[1] = yexp[iTime] + kroady;
    cond[2] = zexp[iTime] - kroadz; cond[3] = zexp[iTime] + kroadz;
    layer->GetClusters(cond, index, ncl);
    for(Int_t ic = 0; ic<ncl; ic++){
      AliTRDcluster *c = layer->GetCluster(index[ic]);
      clusters[nYclusters] = c;
      yres[nYclusters++] = c->GetY() - yexp[iTime];
      if(nYclusters >= kClusterCandidates) {
        AliWarning(Form("Cluster candidates reached limit %d. Some may be lost.", kClusterCandidates));
        kEXIT = kTRUE;
        break;
      }
    }
    tboundary[iTime] = nYclusters;
    if(kEXIT) break;
  }
  
  // Evaluate truncated mean on the y direction
  Double_t mean, sigma;
  AliMathBase::EvaluateUni(nYclusters, yres, mean, sigma, Int_t(nYclusters*.8)-2);
  // purge cluster candidates
  Int_t nZclusters = 0;
  for(Int_t ic = 0; ic<nYclusters; ic++){
    if(yres[ic] - mean > 4. * sigma){
      clusters[ic] = 0x0;
      continue;
    }
    zres[nZclusters++] = clusters[ic]->GetZ() - zexp[clusters[ic]->GetLocalTimeBin()];
  }
  
  // Evaluate truncated mean on the z direction
  AliMathBase::EvaluateUni(nZclusters, zres, mean, sigma, Int_t(nZclusters*.8)-2);
  // purge cluster candidates
  for(Int_t ic = 0; ic<nZclusters; ic++){
    if(zres[ic] - mean > 4. * sigma){
      clusters[ic] = 0x0;
      continue;
    }
  }

  
  // Select only one cluster/TimeBin
  Int_t lastCluster = 0;
  fN2 = 0;
  for (Int_t iTime = 0; iTime < AliTRDtrackerV1::GetNTimeBins(); iTime++) {
    ncl = tboundary[iTime] - lastCluster;
    if(!ncl) continue;
    Int_t iptr = lastCluster;
    if(ncl > 1){
      Float_t dold = 9999.;
      for(int ic=lastCluster; ic<tboundary[iTime]; ic++){
        if(!clusters[ic]) continue;
        Float_t y = yexp[iTime] - clusters[ic]->GetY();
        Float_t z = zexp[iTime] - clusters[ic]->GetZ();
        Float_t d = y * y + z * z;
        if(d > dold) continue;
        dold = d;
        iptr = ic;
      }
    }
    fIndexes[iTime]  = chamber->GetTB(iTime)->GetGlobalIndex(iptr);
    fClusters[iTime] = clusters[iptr];
    fY[iTime]        = clusters[iptr]->GetY();
    fZ[iTime]        = clusters[iptr]->GetZ();
    lastCluster      = tboundary[iTime];
    fN2++;
  }
  
  // number of minimum numbers of clusters expected for the tracklet
  Int_t kClmin = Int_t(fReconstructor->GetRecoParam() ->GetFindableClusters()*AliTRDtrackerV1::GetNTimeBins());
  if (fN2 < kClmin){
    AliWarning(Form("Not enough clusters to fit the tracklet %d [%d].", fN2, kClmin));
    fN2 = 0;
    return kFALSE;
  }

  // update used clusters
  fNUsed = 0;
  for (Int_t iTime = 0; iTime < AliTRDtrackerV1::GetNTimeBins(); iTime++) {
    if(!fClusters[iTime]) continue;
    if((fClusters[iTime]->IsUsed())) fNUsed++;
  }

  if (fN2-fNUsed < kClmin){
    AliWarning(Form("Too many clusters already in use %d (from %d).", fNUsed, fN2));
    fN2 = 0;
    return kFALSE;
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
  fSnp = fYref[1]/TMath::Sqrt(1+fYref[1]*fYref[1]);
  fTgl = fZref[1];
  fN = 0; fN2 = 0; fMPads = 0.;
  AliTRDcluster **cit = &fClusters[0];
  for(Int_t ic = knTimebins; ic--; cit++){
    if(!(*cit)) return;
    fN++; fN2++;
    fX[ic] = (*cit)->GetX() - fX0;
    fY[ic] = (*cit)->GetY();
    fZ[ic] = (*cit)->GetZ();
  }
  Update(); // Fit();
  CookLabels();
  GetProbability();
}


//____________________________________________________________________
Bool_t AliTRDseedV1::Fit(Bool_t tilt)
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

  const Int_t kClmin = 8;
//   const Float_t q0 = 100.;
//   const Float_t clSigma0 = 2.E-2;    //[cm]
//   const Float_t clSlopeQ = -1.19E-2; //[1/cm]

  // get track direction
  Double_t y0   = fYref[0];
  Double_t dydx = fYref[1]; 
  Double_t z0   = fZref[0];
  Double_t dzdx = fZref[1];
  Double_t yt, zt;

  const Int_t kNtb = AliTRDtrackerV1::GetNTimeBins();
  AliTRDtrackerV1::AliTRDLeastSquare fitterZ;
  TLinearFitter  fitterY(1, "pol1");
  // convertion factor from square to gauss distribution for sigma
  Double_t convert = 1./TMath::Sqrt(12.);
  
  // book cluster information
  Double_t xc[knTimebins], yc[knTimebins], zc[knTimebins], sy[knTimebins], sz[knTimebins];
  Int_t zRow[knTimebins];


  fN = 0;
  AliTRDcluster *c=0x0, **jc = &fClusters[0];
  for (Int_t ic=0; ic<kNtb; ic++, ++jc) {
    zRow[ic] = -1;
    xc[ic]  = -1.;
    yc[ic]  = 999.;
    zc[ic]  = 999.;
    sy[ic]  = 0.;
    sz[ic]  = 0.;
    if(!(c = (*jc))) continue;
    if(!c->IsInChamber()) continue;
    Float_t w = 1.;
    if(c->GetNPads()>4) w = .5;
    if(c->GetNPads()>5) w = .2;
    zRow[fN] = c->GetPadRow();
    xc[fN]   = fX0 - c->GetX();
    yc[fN]   = c->GetY();
    zc[fN]   = c->GetZ();

    // extrapolated y value for the track
    yt = y0 - xc[fN]*dydx; 
    // extrapolated z value for the track
    zt = z0 - xc[fN]*dzdx; 
    // tilt correction
    if(tilt) yc[fN] -= fTilt*(zc[fN] - zt); 

    // elaborate cluster error
    //Float_t qr = c->GetQ() - q0;
    sy[fN]   = 1.;//qr < 0. ? clSigma0*TMath::Exp(clSlopeQ*qr) : clSigma0;
    fitterY.AddPoint(&xc[fN], yc[fN]/*-yt*/, sy[fN]);

    sz[fN]   = fPadLength*convert;
    fitterZ.AddPoint(&xc[fN], zc[fN], sz[fN]);
    fN++;
  }
  // to few clusters
  if (fN < kClmin) return kFALSE; 

  // fit XY plane
  fitterY.Eval();
  fYfit[0] = /*y0+*/fitterY.GetParameter(0);
  fYfit[1] = /*dydx-*/-fitterY.GetParameter(1);

  // check par row crossing
  Int_t zN[2*AliTRDseed::knTimebins];
  Int_t nz = AliTRDtrackerV1::Freq(fN, zRow, zN, kFALSE);
  // more than one pad row crossing
  if(nz>2) return kFALSE; 


  // determine z offset of the fit
  Float_t zslope = 0.;
  Int_t nchanges = 0, nCross = 0;
  if(nz==2){ // tracklet is crossing pad row
    // Find the break time allowing one chage on pad-rows
    // with maximal number of accepted clusters
    Int_t padRef = zRow[0];
    for (Int_t ic=1; ic<fN; ic++) {
      if(zRow[ic] == padRef) continue;
      
      // debug
      if(zRow[ic-1] == zRow[ic]){
        printf("ERROR in pad row change!!!\n");
      }
    
      // evaluate parameters of the crossing point
      Float_t sx = (xc[ic-1] - xc[ic])*convert;
      fCross[0] = .5 * (xc[ic-1] + xc[ic]);
      fCross[2] = .5 * (zc[ic-1] + zc[ic]);
      fCross[3] = TMath::Max(dzdx * sx, .01);
      zslope    = zc[ic-1] > zc[ic] ? 1. : -1.;
      padRef    = zRow[ic];
      nCross    = ic;
      nchanges++;
    }
  }

  // condition on nCross and reset nchanges TODO

  if(nchanges==1){
    if(dzdx * zslope < 0.){
      AliInfo("tracklet direction does not correspond to the track direction. TODO.");
    }
    SetBit(kRowCross, kTRUE); // mark pad row crossing
    fitterZ.AddPoint(&fCross[0], fCross[2], fCross[3]);
    fitterZ.Eval();
    //zc[nc] = fitterZ.GetFunctionParameter(0); 
    fCross[1] = fYfit[0] - fCross[0] * fYfit[1];
    fCross[0] = fX0 - fCross[0];
  } else if(nchanges > 1){ // debug
    AliError("N pad row crossing > 1.");
    return kFALSE;
  }

  UpdateUsed();

  return kTRUE;
}


//___________________________________________________________________
void AliTRDseedV1::Print(Option_t *o) const
{
  //
  // Printing the seedstatus
  //

  AliInfo(Form("Det[%3d] Tilt[%+6.2f] Pad[%5.2f]", fDet, fTilt, fPadLength));
  AliInfo(Form("Nattach[%2d] Nfit[%2d] Nuse[%2d] pads[%f]", fN, fN2, fNUsed, fMPads));
  AliInfo(Form("x[%7.2f] y[%7.2f] z[%7.2f] dydx[%5.2f] dzdx[%5.2f]", fX0, fYfit[0], fZfit[0], fYfit[1], fZfit[1]));
  AliInfo(Form("Ref        y[%7.2f] z[%7.2f] dydx[%5.2f] dzdx[%5.2f]", fYref[0], fZref[0], fYref[1], fZref[1]))


  if(strcmp(o, "a")!=0) return;

  AliTRDcluster* const* jc = &fClusters[0];
  for(int ic=0; ic<AliTRDtrackerV1::GetNTimeBins(); ic++, jc++) {
    if(!(*jc)) continue;
    (*jc)->Print(o);
  }

/*  printf("  fSigmaY =%f\n", fSigmaY);
  printf("  fSigmaY2=%f\n", fSigmaY2);            
  printf("  fMeanz  =%f\n", fMeanz);
  printf("  fZProb  =%f\n", fZProb);
  printf("  fLabels[0]=%d fLabels[1]=%d\n", fLabels[0], fLabels[1]);*/
  
/*  printf("  fC      =%f\n", fC);        
  printf("  fCC     =%f\n",fCC);      
  printf("  fChi2   =%f\n", fChi2);  
  printf("  fChi2Z  =%f\n", fChi2Z);*/
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
    if ( fYref[i] != inTracklet->GetYref(i) ) return kFALSE;
    if ( fZref[i] != inTracklet->GetZref(i) ) return kFALSE;
  }
  
  if ( fSigmaY != inTracklet->GetSigmaY() ) return kFALSE;
  if ( fSigmaY2 != inTracklet->GetSigmaY2() ) return kFALSE;
  if ( fTilt != inTracklet->GetTilt() ) return kFALSE;
  if ( fPadLength != inTracklet->GetPadLength() ) return kFALSE;
  
  for (Int_t i = 0; i < knTimebins; i++){
    if ( fX[i] != inTracklet->GetX(i) ) return kFALSE;
    if ( fY[i] != inTracklet->GetY(i) ) return kFALSE;
    if ( fZ[i] != inTracklet->GetZ(i) ) return kFALSE;
    if ( fIndexes[i] != inTracklet->GetIndexes(i) ) return kFALSE;
    if ( fUsable[i] != inTracklet->IsUsable(i) ) return kFALSE;
  }

  for (Int_t i=0; i < 2; i++){
    if ( fYfit[i] != inTracklet->GetYfit(i) ) return kFALSE;
    if ( fZfit[i] != inTracklet->GetZfit(i) ) return kFALSE;
    if ( fYfitR[i] != inTracklet->GetYfitR(i) ) return kFALSE;
    if ( fZfitR[i] != inTracklet->GetZfitR(i) ) return kFALSE;
    if ( fLabels[i] != inTracklet->GetLabels(i) ) return kFALSE;
  }
  
  if ( fMeanz != inTracklet->GetMeanz() ) return kFALSE;
  if ( fZProb != inTracklet->GetZProb() ) return kFALSE;
  if ( fN2 != inTracklet->GetN2() ) return kFALSE;
  if ( fNUsed != inTracklet->GetNUsed() ) return kFALSE;
  if ( fFreq != inTracklet->GetFreq() ) return kFALSE;
  if ( fNChange != inTracklet->GetNChange() ) return kFALSE;
  if ( fNChange != inTracklet->GetNChange() ) return kFALSE;
   
  if ( fC != inTracklet->GetC() ) return kFALSE;
  if ( fCC != inTracklet->GetCC() ) return kFALSE;
  if ( fChi2 != inTracklet->GetChi2() ) return kFALSE;
  //  if ( fChi2Z != inTracklet->GetChi2Z() ) return kFALSE;

  if ( fDet != inTracklet->GetDetector() ) return kFALSE;
  if ( fMom != inTracklet->GetMomentum() ) return kFALSE;
  if ( fdX != inTracklet->GetdX() ) return kFALSE;
  
  for (Int_t iCluster = 0; iCluster < knTimebins; iCluster++){
    AliTRDcluster *curCluster = fClusters[iCluster];
    AliTRDcluster *inCluster = inTracklet->GetClusters(iCluster);
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
