
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

#include "AliLog.h"
#include "AliESDtrack.h"
#include "AliTracker.h"

#include "AliTRDtrackV1.h"
#include "AliTRDcluster.h"
#include "AliTRDcalibDB.h"
#include "AliTRDReconstructor.h"
#include "AliTRDPIDResponse.h"
#include "AliTRDrecoParam.h"

ClassImp(AliTRDtrackV1)

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  Represents a reconstructed TRD track                                     //
//  Local TRD Kalman track                                                   //
//                                                                           //
//  Authors:                                                                 //
//    Alex Bercuci <A.Bercuci@gsi.de>                                        //
//    Markus Fasel <M.Fasel@gsi.de>                                          //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

//_______________________________________________________________
AliTRDtrackV1::AliTRDtrackV1() : AliKalmanTrack()
  ,fStatus(0)
  ,fESDid(0)
  ,fDE(0.)
  ,fkReconstructor(NULL)
  ,fBackupTrack(NULL)
  ,fTrackLow(NULL)
  ,fTrackHigh(NULL)
{
  //
  // Default constructor
  //
  //printf("AliTRDtrackV1::AliTRDtrackV1()\n");

  for(int i =0; i<3; i++) fBudget[i] = 0.;
  
  Float_t pid = 1./AliPID::kSPECIES;
  for(int is =0; is<AliPID::kSPECIES; is++) fPID[is] = pid;

  for(int ip=0; ip<kNplane; ip++){
    fTrackletIndex[ip] = -1;
    fTracklet[ip]      = NULL;
  }
}

//_______________________________________________________________
AliTRDtrackV1::AliTRDtrackV1(const AliTRDtrackV1 &ref) : AliKalmanTrack(ref)
  ,fStatus(ref.fStatus)
  ,fESDid(ref.fESDid)
  ,fDE(ref.fDE)
  ,fkReconstructor(ref.fkReconstructor)
  ,fBackupTrack(NULL)
  ,fTrackLow(NULL)
  ,fTrackHigh(NULL)
{
  //
  // Copy constructor
  //

  //printf("AliTRDtrackV1::AliTRDtrackV1(const AliTRDtrackV1 &)\n");
  SetBit(kOwner, kFALSE);
  for(int ip=0; ip<kNplane; ip++){ 
    fTrackletIndex[ip] = ref.fTrackletIndex[ip];
    fTracklet[ip]      = ref.fTracklet[ip];
  }
  if(ref.fTrackLow) fTrackLow = new AliExternalTrackParam(*ref.fTrackLow);
  if(ref.fTrackHigh) fTrackHigh = new AliExternalTrackParam(*ref.fTrackHigh);
 
  for (Int_t i = 0; i < 3;i++) fBudget[i]      = ref.fBudget[i];

	for(Int_t is = 0; is<AliPID::kSPECIES; is++) fPID[is] = ref.fPID[is];
  
  AliKalmanTrack::SetNumberOfClusters(ref.GetNumberOfClusters());
}

//_______________________________________________________________
AliTRDtrackV1::AliTRDtrackV1(const AliESDtrack &t) : AliKalmanTrack()
  ,fStatus(0)
  ,fESDid(0)
  ,fDE(0.)
  ,fkReconstructor(NULL)
  ,fBackupTrack(NULL)
  ,fTrackLow(NULL)
  ,fTrackHigh(NULL)
{
  //
  // Constructor from AliESDtrack
  //

  SetESDid(t.GetID());
  SetLabel(t.GetLabel());
  SetChi2(0.0);

  SetMass(t.GetMass()/*0.000510*/);
  AliKalmanTrack::SetNumberOfClusters(t.GetTRDncls()); 
  Int_t ti[]={-1, -1, -1, -1, -1, -1}; t.GetTRDtracklets(&ti[0]);
  for(int ip=0; ip<kNplane; ip++){ 
    fTrackletIndex[ip] = ti[ip];
    fTracklet[ip]      = NULL;
  }
  for(int i =0; i<3; i++) fBudget[i] = 0.;
  
  Float_t pid = 1./AliPID::kSPECIES;
  for(int is =0; is<AliPID::kSPECIES; is++) fPID[is] = pid;

  const AliExternalTrackParam *par = &t;
  if (t.GetStatus() & AliESDtrack::kTRDbackup) { 
    par = t.GetOuterParam();
    if (!par) {
      AliError("No backup info!"); 
      par = &t;
    }
  }
  Set(par->GetX() 
     ,par->GetAlpha()
     ,par->GetParameter()
     ,par->GetCovariance());

  if(t.GetStatus() & AliESDtrack::kTIME) {
    StartTimeIntegral();
    Double_t times[10]; 
    t.GetIntegratedTimes(times); 
    SetIntegratedTimes(times);
    SetIntegratedLength(t.GetIntegratedLength());
  }
}

//_______________________________________________________________
AliTRDtrackV1::AliTRDtrackV1(AliTRDseedV1 * const trklts, const Double_t p[5], const Double_t cov[15]
             , Double_t x, Double_t alpha) : AliKalmanTrack()
  ,fStatus(0)
  ,fESDid(0)
  ,fDE(0.)
  ,fkReconstructor(NULL)
  ,fBackupTrack(NULL)
  ,fTrackLow(NULL)
  ,fTrackHigh(NULL)
{
  //
  // The stand alone tracking constructor
  // TEMPORARY !!!!!!!!!!!
  // to check :
  // 1. covariance matrix
  // 2. dQdl calculation
  //

  Double_t b(GetBz());
  Double_t cnv = (TMath::Abs(b) < 1.e-5) ? 1.e5 : 1./GetBz()/kB2C;
  
  Double_t pp[5] = { p[0]    
                    , p[1]
                    , p[2]
                    , p[3]
                    , p[4]*cnv      };
  
  Double_t c22 = x*x*cov[14] - 2*x*cov[12] + cov[ 5];
  Double_t c32 =   x*cov[13] -     cov[ 8];
  Double_t c20 =   x*cov[10] -     cov[ 3];
  Double_t c21 =   x*cov[11] -     cov[ 4];
  Double_t c42 =   x*cov[14] -     cov[12];
  
  Double_t cc[15] = { cov[ 0]
                    , cov[ 1],     cov[ 2]
                    , c20,         c21,         c22
                    , cov[ 6],     cov[ 7],     c32,     cov[ 9]
                    , cov[10]*cnv, cov[11]*cnv, c42*cnv, cov[13]*cnv, cov[14]*cnv*cnv };
  
  Double_t mostProbablePt=AliExternalTrackParam::GetMostProbablePt();
  Double_t p0=TMath::Sign(1/mostProbablePt,pp[4]);
  Double_t w0=cc[14]/(cc[14] + p0*p0), w1=p0*p0/(cc[14] + p0*p0);
  AliDebug(3, Form("Pt mixing : w0[%4.2f] pt0[%5.3f] w1[%4.2f] pt[%5.3f]", w0, 1./p0, w1, 1./pp[4]));
  
  pp[4] = w0*p0 + w1*pp[4];
  cc[10]*=w1; cc[11]*=w1; cc[12]*=w1; cc[13]*=w1; cc[14]*=w1;
  Set(x,alpha,pp,cc);
  AliDebug(2, Form("Init @ x[%6.2f] pt[%5.3f]", x, 1./pp[4]));
  Int_t ncls = 0;
  for(int iplane=0; iplane<kNplane; iplane++){
    fTrackletIndex[iplane] = -1;
	  if(!trklts[iplane].IsOK()) fTracklet[iplane] = NULL;
    else{ 
      fTracklet[iplane] = &trklts[iplane];
      ncls += fTracklet[iplane]->GetN();
    }
  }
  AliKalmanTrack::SetNumberOfClusters(ncls);		
  for(int i =0; i<3; i++) fBudget[i] = 0.;
  
  Float_t pid = 1./AliPID::kSPECIES;
  for(int is =0; is<AliPID::kSPECIES; is++) fPID[is] = pid;
}

//_______________________________________________________________
AliTRDtrackV1::~AliTRDtrackV1()
{
  // Clean up all objects allocated by the track during its lifetime.
  AliDebug(2, Form("Deleting track[%d]\n   fBackupTrack[%p] fTrackLow[%p] fTrackHigh[%p] Owner[%c].", fESDid, (void*)fBackupTrack, (void*)fTrackLow, (void*)fTrackHigh, TestBit(kOwner)?'y':'n'));

  if(fBackupTrack) delete fBackupTrack; fBackupTrack = NULL;

  if(fTrackLow) delete fTrackLow; fTrackLow = NULL;
  if(fTrackHigh) delete fTrackHigh; fTrackHigh = NULL;

  for(Int_t ip=0; ip<kNplane; ip++){
    if(TestBit(kOwner) && fTracklet[ip]) delete fTracklet[ip];
    fTracklet[ip] = NULL;
    fTrackletIndex[ip] = -1;
  }
}
	
//_______________________________________________________________
AliTRDtrackV1 &AliTRDtrackV1::operator=(const AliTRDtrackV1 &t)
{
  //
  // Assignment operator
  //

  if (this != &t) {
    AliKalmanTrack::operator=(t);
    ((AliTRDtrackV1 &) t).Copy(*this);
  }

  return *this;

}

//_____________________________________________________________________________
void AliTRDtrackV1::Copy(TObject &t) const
{
  //
  // Copy function
  //

  ((AliTRDtrackV1 &) t).fStatus         = fStatus;
  ((AliTRDtrackV1 &) t).fESDid          = fESDid;
  ((AliTRDtrackV1 &) t).fDE             = fDE;
  ((AliTRDtrackV1 &) t).fkReconstructor = fkReconstructor;
  ((AliTRDtrackV1 &) t).fBackupTrack    = 0x0;
  ((AliTRDtrackV1 &) t).fTrackLow       = 0x0;
  ((AliTRDtrackV1 &) t).fTrackHigh      = 0x0;

  for(Int_t ip = 0; ip < kNplane; ip++) { 
    ((AliTRDtrackV1 &) t).fTrackletIndex[ip] = fTrackletIndex[ip];
    ((AliTRDtrackV1 &) t).fTracklet[ip]      = fTracklet[ip];
  }
  if (fTrackLow) {
    ((AliTRDtrackV1 &) t).fTrackLow  = new AliExternalTrackParam(*fTrackLow);
  }
  if (fTrackHigh){
    ((AliTRDtrackV1 &) t).fTrackHigh = new AliExternalTrackParam(*fTrackHigh);
  }
 
  for (Int_t i = 0; i < 3; i++) {
    ((AliTRDtrackV1 &) t).fBudget[i] = fBudget[i];
  }
  for (Int_t is = 0; is < AliPID::kSPECIES; is++) {
    ((AliTRDtrackV1 &) t).fPID[is] = fPID[is];
  }  

}

//_______________________________________________________________
Bool_t AliTRDtrackV1::CookLabel(Float_t wrong)
{
  // set MC label for this track
  if(!GetNumberOfClusters()) return kFALSE;

  Int_t s[kMAXCLUSTERSPERTRACK][2];
  for (Int_t i = 0; i < kMAXCLUSTERSPERTRACK; i++) {
    s[i][0] = -1;
    s[i][1] =  0;
  }
  
  Bool_t labelAdded;
  Int_t label;
  AliTRDcluster *c    = NULL;
  for (Int_t ip = 0; ip < kNplane; ip++) {
    if(fTrackletIndex[ip]<0 || !fTracklet[ip]) continue;
    for (Int_t ic = 0; ic < AliTRDseedV1::kNclusters; ic++) {
      if(!(c = fTracklet[ip]->GetClusters(ic))) continue;
      for (Int_t k = 0; k < 3; k++) { 
        label      = c->GetLabel(k);
        labelAdded = kFALSE; 
        Int_t j = 0;
        if (label >= 0) {
          while ((!labelAdded) && (j < kMAXCLUSTERSPERTRACK)) {
            if ((s[j][0] == label) || 
                (s[j][1] ==     0)) {
              s[j][0] = label; 
              s[j][1]++; 
              labelAdded = kTRUE;
            }
            j++;
          }
        }
      }
    }
  }
  
  Int_t max = 0;
  label = -123456789;
  for (Int_t i = 0; i < kMAXCLUSTERSPERTRACK; i++) {
    if (s[i][1] <= max) continue;
    max   = s[i][1]; 
    label = s[i][0];
  }
  if ((1. - Float_t(max)/GetNumberOfClusters()) > wrong) label = -label;
  
  SetLabel(label); 
  
  return kTRUE;
}

//_______________________________________________________________
Bool_t AliTRDtrackV1::CookPID()
{
//
// Cook the PID information for the track by delegating the omonim function of the tracklets. 
// Computes the number of tracklets used. The tracklet information are considered independent. 
// For the moment no global track measurement of PID is performed as for example to estimate 
// bremsstrahlung probability based on global chi2 of the track.
//
// The status bit AliESDtrack::kTRDpid is set during the call of AliTRDtrackV1::UpdateESDtrack().The PID performance of the 
//TRD for tracks with 6 tacklets is displayed below.
//Begin_Html
//<img src="TRD/trackPID.gif">
//End_Html
//
  const AliTRDPIDResponse *pidResponse = AliTRDcalibDB::Instance()->GetPIDResponse(fkReconstructor->GetRecoParam()->GetPIDmethod());
  if(!pidResponse){
    AliError("PID Response not available");
    return kFALSE;
  }
  Int_t nslices = pidResponse->GetNumberOfSlices();
  Double_t dEdx[kNplane * (Int_t)AliTRDPIDResponse::kNslicesNN];
  Float_t trackletP[kNplane];
  memset(dEdx, 0, sizeof(Double_t) * kNplane * (Int_t)AliTRDPIDResponse::kNslicesNN);
  memset(trackletP, 0, sizeof(Float_t)*kNplane);
  for(Int_t iseed = 0; iseed < kNplane; iseed++){
    if(!fTracklet[iseed]) continue;
    trackletP[iseed] = fTracklet[iseed]->GetMomentum();
    fTracklet[iseed]->SetPID();
    if(pidResponse->GetPIDmethod() == AliTRDPIDResponse::kLQ1D){
      dEdx[iseed] = fTracklet[iseed]->GetdQdl();
    } else {
      fTracklet[iseed]->CookdEdx(nslices);
      const Float_t *trackletdEdx = fTracklet[iseed]->GetdEdx();
      for(Int_t islice = 0; islice < nslices; islice++){
        dEdx[iseed*nslices + islice] = trackletdEdx[islice];
      }
    }
  }
  pidResponse->GetResponse(nslices, dEdx, trackletP, fPID);
  return kTRUE;
}

//___________________________________________________________
UChar_t AliTRDtrackV1::GetNumberOfTrackletsPID() const
{
// Retrieve number of tracklets used for PID calculation. 

  UChar_t nPID = 0;
  for(int ip=0; ip<kNplane; ip++){
    if(fTrackletIndex[ip]<0 || !fTracklet[ip]) continue;
    if(!fTracklet[ip]->IsOK()) continue;
    
    nPID++;
  }
  return nPID;
}

//_______________________________________________________________
AliTRDcluster* AliTRDtrackV1::GetCluster(Int_t id)
{
  // Get the cluster at a certain position in the track
  Int_t n = 0;
  for(Int_t ip=0; ip<kNplane; ip++){
    if(!fTracklet[ip]) continue;
    if(n+fTracklet[ip]->GetN() <= id){ 
      n+=fTracklet[ip]->GetN();
      continue;
    }
    AliTRDcluster *c = NULL;
    for(Int_t ic=AliTRDseedV1::kNclusters; ic--;){
      if(!(c = fTracklet[ip]->GetClusters(ic))) continue;

      if(n<id){n++; continue;}
      return c;
    }
  }
  return NULL;
}

//_______________________________________________________________
Int_t  AliTRDtrackV1::GetClusterIndex(Int_t id) const
{
  // Get the cluster index at a certain position in the track
  Int_t n = 0;
  for(Int_t ip=0; ip<kNplane; ip++){
    if(!fTracklet[ip]) continue;
    if(n+fTracklet[ip]->GetN() <= id){ 
      n+=fTracklet[ip]->GetN();
      continue;
    }
    for(Int_t ic=AliTRDseedV1::kNclusters; ic--;){
      if(!(fTracklet[ip]->GetClusters(ic))) continue;
      if(n<id){n++; continue;}
      return fTracklet[ip]->GetIndexes(ic);
    }
  }
  return -1;
}

//_______________________________________________________________
Double_t AliTRDtrackV1::GetPredictedChi2(const AliTRDseedV1 *trklt, Double_t *cov) const
{
// Compute chi2 between tracklet and track. The value is calculated at the radial position of the track
// equal to the reference radial position of the tracklet (see AliTRDseedV1)
// 
// The chi2 estimator is computed according to the following formula
// BEGIN_LATEX
// #chi^{2}=(X_{trklt}-X_{track})(C_{trklt}+C_{track})^{-1}(X_{trklt}-X_{track})^{T}
// END_LATEX
// where X=(y z), the position of the track/tracklet in the yz plane
// 

  Double_t p[2]   = { trklt->GetY(), trklt->GetZ()};
  trklt->GetCovAt(trklt->GetX(), cov);
  return AliExternalTrackParam::GetPredictedChi2(p, cov);
}

//_______________________________________________________________
Int_t AliTRDtrackV1::GetSector() const
{
  return Int_t(GetAlpha()/AliTRDgeometry::GetAlpha() + (GetAlpha()>0. ? 0 : AliTRDgeometry::kNsector));
}

//_______________________________________________________________
Bool_t AliTRDtrackV1::IsEqual(const TObject *o) const
{
  // Checks whether two tracks are equal
  if (!o) return kFALSE;
  const AliTRDtrackV1 *inTrack = dynamic_cast<const AliTRDtrackV1*>(o);
  if (!inTrack) return kFALSE;
  
  //if ( fPIDquality != inTrack->GetPIDquality() ) return kFALSE;
  
  if(memcmp(fPID, inTrack->fPID, AliPID::kSPECIES*sizeof(Double32_t))) return kFALSE;
  if(memcmp(fBudget, inTrack->fBudget, 3*sizeof(Double32_t))) return kFALSE;
  if(memcmp(&fDE, &inTrack->fDE, sizeof(Double32_t))) return kFALSE;
  if(memcmp(&fFakeRatio, &inTrack->fFakeRatio, sizeof(Double32_t))) return kFALSE;
  if(memcmp(&fChi2, &inTrack->fChi2, sizeof(Double32_t))) return kFALSE;
  if(memcmp(&fMass, &inTrack->fMass, sizeof(Double32_t))) return kFALSE;
  if( fLab != inTrack->fLab ) return kFALSE;
  if( fN != inTrack->fN ) return kFALSE;
  Double32_t l(0.), in(0.);
  l = GetIntegratedLength(); in = inTrack->GetIntegratedLength();
  if(memcmp(&l, &in, sizeof(Double32_t))) return kFALSE;
  l=GetX(); in=inTrack->GetX();
  if(memcmp(&l, &in, sizeof(Double32_t))) return kFALSE;
  l = GetAlpha(); in = inTrack->GetAlpha();
  if(memcmp(&l, &in, sizeof(Double32_t))) return kFALSE;
  if(memcmp(GetParameter(), inTrack->GetParameter(), 5*sizeof(Double32_t))) return kFALSE;
  if(memcmp(GetCovariance(), inTrack->GetCovariance(), 15*sizeof(Double32_t))) return kFALSE;
  
  for (Int_t iTracklet = 0; iTracklet < kNplane; iTracklet++){
    AliTRDseedV1 *curTracklet = fTracklet[iTracklet];
    AliTRDseedV1 *inTracklet = inTrack->GetTracklet(iTracklet);
    if (curTracklet && inTracklet){
      if (! curTracklet->IsEqual(inTracklet) ) {
        curTracklet->Print();
        inTracklet->Print();
        return kFALSE;
      }
    } else {
      // if one tracklet exists, and corresponding 
      // in other track doesn't - return kFALSE
      if(inTracklet || curTracklet) return kFALSE;
    }
  }

  return kTRUE;
}

//_______________________________________________________________
Bool_t AliTRDtrackV1::IsElectron() const
{
  if(GetPID(0) > fkReconstructor->GetRecoParam()->GetPIDThreshold(GetP())) return kTRUE;
  return kFALSE;
}

	
//_____________________________________________________________________________
Int_t AliTRDtrackV1::MakeBackupTrack()
{
//
// Creates a backup track
// TO DO update quality check of the track.
//

  Float_t occupancy(0.); Int_t n(0), ncls(0);
  for(Int_t il(AliTRDgeometry::kNlayer); il--;){ 
    if(!fTracklet[il]) continue;
    n++; 
    occupancy+=fTracklet[il]->GetTBoccupancy()/AliTRDseedV1::kNtb;
    ncls += fTracklet[il]->GetN();
  }
  if(!n) return -1;
  occupancy/=n;

  //Float_t ratio1 = Float_t(t.GetNumberOfClusters()+1) / Float_t(t.GetNExpected()+1);  
  
  Int_t failedCutId(0);
  if(GetChi2()/n > 5.0) failedCutId=1; 
  if(occupancy < 0.7) failedCutId=2;
  //if(ratio1 >   0.6) && 
  //if(ratio0+ratio1  >   1.5) && 
  if(GetNCross() != 0)  failedCutId=3;
  if(TMath::Abs(GetSnp()) > 0.85) failedCutId=4;
  if(ncls < 20) failedCutId=5;

  if(failedCutId){ 
    AliDebug(2, Form("\n"
      "chi2/tracklet < 5.0   [%c] %5.2f\n"
      "occupancy > 0.7       [%c] %4.2f\n"
      "NCross == 0           [%c] %d\n"
      "Abs(snp) < 0.85       [%c] %4.2f\n"
      "NClusters > 20        [%c] %d"
      ,(GetChi2()/n<5.0)?'y':'n', GetChi2()/n
      ,(occupancy>0.7)?'y':'n', occupancy
      ,(GetNCross()==0)?'y':'n', GetNCross()
      ,(TMath::Abs(GetSnp())<0.85)?'y':'n', TMath::Abs(GetSnp())
      ,(ncls>20)?'y':'n', ncls
    ));
    return failedCutId;
  }

  if(fBackupTrack) {
    fBackupTrack->~AliTRDtrackV1();
    new(fBackupTrack) AliTRDtrackV1((AliTRDtrackV1&)(*this));
    return 0;
  }
  fBackupTrack = new AliTRDtrackV1((AliTRDtrackV1&)(*this));
  return 0;
}

//_____________________________________________________________________________
Int_t AliTRDtrackV1::GetProlongation(Double_t xk, Double_t &y, Double_t &z) const
{
  //
  // Find a prolongation at given x
  // Return -1 if it does not exist
  //  

  Double_t bz = GetBz();
  if (!AliExternalTrackParam::GetYAt(xk,bz,y)) return -1;
  if (!AliExternalTrackParam::GetZAt(xk,bz,z)) return -1;

  return 1;  

}

//_____________________________________________________________________________
Bool_t AliTRDtrackV1::PropagateTo(Double_t xk, Double_t xx0, Double_t xrho)
{
  //
  // Propagates this track to a reference plane defined by "xk" [cm] 
  // correcting for the mean crossed material.
  //
  // "xx0"  - thickness/rad.length [units of the radiation length] 
  // "xrho" - thickness*density    [g/cm^2] 
  // 

  if (TMath::Abs(xk - GetX())<AliTRDReconstructor::GetEpsilon()*0.1) return kTRUE; // 10% of the tracker precision

  Double_t xyz0[3] = {GetX(), GetY(), GetZ()}, // track position BEFORE propagation 
           b[3];    // magnetic field 
  GetBxByBz(b);
  if(!AliExternalTrackParam::PropagateToBxByBz(xk,b)) return kFALSE;
 
  // local track position AFTER propagation 
  Double_t xyz1[3] = {GetX(), GetY(), GetZ()};
//  printf("x0[%6.2f] -> x1[%6.2f] dx[%6.2f] rho[%f]\n", xyz0[0], xyz1[0], xyz0[0]-xk, xrho/TMath::Abs(xyz0[0]-xk));
  if(xyz0[0] < xk) {
    xrho = -xrho;
    if (IsStartedTimeIntegral()) {
      Double_t l2  = TMath::Sqrt((xyz1[0]-xyz0[0])*(xyz1[0]-xyz0[0]) 
                               + (xyz1[1]-xyz0[1])*(xyz1[1]-xyz0[1]) 
                               + (xyz1[2]-xyz0[2])*(xyz1[2]-xyz0[2]));
      Double_t crv = AliExternalTrackParam::GetC(b[2]);
      if (TMath::Abs(l2*crv) > 0.0001) {
        // Make correction for curvature if neccesary
        l2 = 0.5 * TMath::Sqrt((xyz1[0]-xyz0[0])*(xyz1[0]-xyz0[0]) 
                             + (xyz1[1]-xyz0[1])*(xyz1[1]-xyz0[1]));
        l2 = 2.0 * TMath::ASin(l2 * crv) / crv;
        l2 = TMath::Sqrt(l2*l2 + (xyz1[2]-xyz0[2])*(xyz1[2]-xyz0[2]));
      }
      AddTimeStep(l2);
    }
  }
  if (!AliExternalTrackParam::CorrectForMeanMaterial(xx0, xrho, GetMass())) return kFALSE;


  {

    // Energy losses
    Double_t p2    = (1.0 + GetTgl()*GetTgl()) / (GetSigned1Pt()*GetSigned1Pt());
    Double_t beta2 = p2 / (p2 + GetMass()*GetMass());
    if ((beta2 < 1.0e-10) || 
        ((5940.0 * beta2/(1.0 - beta2 + 1.0e-10) - beta2) < 0.0)) {
      return kFALSE;
    }

    Double_t dE    = 0.153e-3 / beta2 
                   * (TMath::Log(5940.0 * beta2/(1.0 - beta2 + 1.0e-10)) - beta2)
                   * xrho;
    fBudget[0] += xrho;

    /*
    // Suspicious part - think about it ?
    Double_t kinE =  TMath::Sqrt(p2);
    if (dE > 0.8*kinE) dE = 0.8 * kinE;  //      
    if (dE < 0)        dE = 0.0;         // Not valid region for Bethe bloch 
    */
 
    fDE += dE;

    /*
    // Suspicious ! I.B.
    Double_t sigmade = 0.07 * TMath::Sqrt(TMath::Abs(dE));   // Energy loss fluctuation 
    Double_t sigmac2 = sigmade*sigmade*fC*fC*(p2+GetMass()*GetMass())/(p2*p2);
    fCcc += sigmac2;
    fCee += fX*fX * sigmac2;  
    */

  }

  return kTRUE;
}

//_____________________________________________________________________________
Int_t   AliTRDtrackV1::PropagateToR(Double_t r,Double_t step)
{
  //
  // Propagate track to the radial position
  // Rotation always connected to the last track position
  //

  Double_t xyz0[3];
  Double_t xyz1[3];
  Double_t y;
  Double_t z; 

  Double_t radius = TMath::Sqrt(GetX()*GetX() + GetY()*GetY());
  // Direction +-
  Double_t dir    = (radius > r) ? -1.0 : 1.0;   

  for (Double_t x = radius+dir*step; dir*x < dir*r; x += dir*step) {

    GetXYZ(xyz0);	
    Double_t alpha = TMath::ATan2(xyz0[1],xyz0[0]);
    Rotate(alpha,kTRUE);
    GetXYZ(xyz0);	
    if(GetProlongation(x,y,z)<0) return -1;
    xyz1[0] = x * TMath::Cos(alpha) + y * TMath::Sin(alpha); 
    xyz1[1] = x * TMath::Sin(alpha) - y * TMath::Cos(alpha);
    xyz1[2] = z;
    Double_t param[7];
    if(AliTracker::MeanMaterialBudget(xyz0,xyz1,param)<=0.) return -1;
    if (param[1] <= 0) {
      param[1] = 100000000;
    }
    PropagateTo(x,param[1],param[0]*param[4]);

  } 

  GetXYZ(xyz0);	
  Double_t alpha = TMath::ATan2(xyz0[1],xyz0[0]);
  Rotate(alpha,kTRUE);
  GetXYZ(xyz0);	
  if(GetProlongation(r,y,z)<0) return -1;
  xyz1[0] = r * TMath::Cos(alpha) + y * TMath::Sin(alpha); 
  xyz1[1] = r * TMath::Sin(alpha) - y * TMath::Cos(alpha);
  xyz1[2] = z;
  Double_t param[7];
  if(AliTracker::MeanMaterialBudget(xyz0,xyz1,param) <= 0.) return -1;

  if (param[1] <= 0) {
    param[1] = 100000000;
  }
  PropagateTo(r,param[1],param[0]*param[4]);

  return 0;

}

//_____________________________________________________________________________
void AliTRDtrackV1::Print(Option_t *o) const
{
  // Print track status
  AliInfo(Form("PID [%4.1f %4.1f %4.1f %4.1f %4.1f]", 1.E2*fPID[0], 1.E2*fPID[1], 1.E2*fPID[2], 1.E2*fPID[3], 1.E2*fPID[4]));
  AliInfo(Form("Material[%5.2f %5.2f %5.2f]", fBudget[0], fBudget[1], fBudget[2]));

  AliInfo(Form("x[%7.2f] t[%7.4f] alpha[%f] mass[%f]", GetX(), GetIntegratedLength(), GetAlpha(), fMass));
  AliInfo(Form("Ntr[%1d] NtrPID[%1d] Ncl[%3d] lab[%3d]", GetNumberOfTracklets(), GetNumberOfTrackletsPID(), fN, fLab));

  printf("|X| = (");
  const Double_t *curP = GetParameter();
  for (Int_t i = 0; i < 5; i++) printf("%7.2f ", curP[i]);
  printf(")\n");

  printf("|V| = \n");
  const Double_t *curC = GetCovariance();
  for (Int_t i = 0, j=4, k=0; i<15; i++, k++){
    printf("%7.2f ", curC[i]);
    if(k==j){ 
      printf("\n");
      k=-1; j--;
    }
  }
  if(strcmp(o, "a")!=0) return;

  for(Int_t ip=0; ip<kNplane; ip++){
    if(!fTracklet[ip]) continue;
    fTracklet[ip]->Print(o);
  }
}


//_____________________________________________________________________________
Bool_t AliTRDtrackV1::Rotate(Double_t alpha, Bool_t absolute)
{
  //
  // Rotates track parameters in R*phi plane
  // if absolute rotation alpha is in global system
  // otherwise alpha rotation is relative to the current rotation angle
  //  

  if (absolute) alpha -= GetAlpha();
  //else fNRotate++;

  return AliExternalTrackParam::Rotate(GetAlpha()+alpha);
}

//___________________________________________________________
void AliTRDtrackV1::SetNumberOfClusters() 
{
// Calculate the number of clusters attached to this track
	
  Int_t ncls = 0;
  for(int ip=0; ip<kNplane; ip++){
    if(fTracklet[ip] && fTrackletIndex[ip] >= 0) ncls += fTracklet[ip]->GetN();
  }
  AliKalmanTrack::SetNumberOfClusters(ncls);	
}

	
//_______________________________________________________________
void AliTRDtrackV1::SetOwner()
{
  //
  // Toggle ownership of tracklets
  //

  if(TestBit(kOwner)) return;
  for (Int_t ip = 0; ip < kNplane; ip++) {
    if(fTrackletIndex[ip]<0 || !fTracklet[ip]) continue;
    fTracklet[ip] = new AliTRDseedV1(*fTracklet[ip]);
    fTracklet[ip]->SetOwner();
  }
  SetBit(kOwner);
}

//_______________________________________________________________
void AliTRDtrackV1::SetTracklet(AliTRDseedV1 *const trklt, Int_t index)
{
  //
  // Set the tracklets
  //
  Int_t plane = trklt->GetPlane();

  fTracklet[plane]      = trklt;
  fTrackletIndex[plane] = index;
}

//_______________________________________________________________
void AliTRDtrackV1::SetTrackIn()
{
//  Save location of birth for the TRD track
//  If the pointer is not valid allocate memory
//
  const AliExternalTrackParam *op = dynamic_cast<const AliExternalTrackParam*>(this);

  //printf("SetTrackIn() : fTrackLow[%p]\n", (void*)fTrackLow);
  if(fTrackLow){
    fTrackLow->~AliExternalTrackParam();
    new(fTrackLow) AliExternalTrackParam(*op);
  } else fTrackLow = new AliExternalTrackParam(*op);
}

//_______________________________________________________________
void AliTRDtrackV1::SetTrackOut(const AliExternalTrackParam *op)
{
//  Save location of death for the TRD track
//  If the pointer is not valid allocate memory
//
  if(!op) op = dynamic_cast<const AliExternalTrackParam*>(this);
  if(fTrackHigh){
    fTrackHigh->~AliExternalTrackParam();
    new(fTrackHigh) AliExternalTrackParam(*op);
  } else fTrackHigh = new AliExternalTrackParam(*op);
}

//_______________________________________________________________
void AliTRDtrackV1::UnsetTracklet(Int_t plane)
{
  if(plane<0) return;
  fTrackletIndex[plane] = -1;
  fTracklet[plane] = NULL;
}


//_______________________________________________________________
void AliTRDtrackV1::UpdateChi2(Float_t chi2)
{
// Update chi2/track with one tracklet contribution
  SetChi2(GetChi2() + chi2);
}

//_______________________________________________________________
void AliTRDtrackV1::UpdateESDtrack(AliESDtrack *track)
{
  //
  // Update the TRD PID information in the ESD track
  //

//   Int_t nslices = AliTRDcalibDB::Instance()->GetPIDResponse(fkReconstructor->GetRecoParam()->GetPIDmethod())->GetNumberOfSlices();
  // number of tracklets used for PID calculation
  UChar_t nPID = GetNumberOfTrackletsPID();
  // number of tracklets attached to the track
  UChar_t nTrk = GetNumberOfTracklets();
  // pack the two numbers together and store them in the ESD
  track->SetTRDntracklets(nPID | (nTrk<<3));
  // allocate space to store raw PID signals dEdx & momentum
  track->SetNumberOfTRDslices((AliTRDPIDResponse::kNslicesNN+3)*AliTRDgeometry::kNlayer);
  // store raw signals
  Float_t p, sp; Double_t spd;
  for (Int_t ip = 0; ip < kNplane; ip++) {
    if(fTrackletIndex[ip]<0 || !fTracklet[ip]) continue;
    if(!fTracklet[ip]->HasPID()) continue;
    fTracklet[ip]->CookdEdx(AliTRDPIDResponse::kNslicesNN);
    const Float_t *dedx = fTracklet[ip]->GetdEdx();
    for (Int_t js = 0; js < AliTRDPIDResponse::kNslicesNN; js++, dedx++){
      track->SetTRDslice(*dedx, ip, js+1);
    }
    p = fTracklet[ip]->GetMomentum(&sp); 
    // store global quality per tracklet instead of momentum error
    // 26.09.11 A.Bercuci
    // first implementation store no. of time bins filled in tracklet (5bits  see "y" bits) and
    // no. of double clusters in case of pad row cross (4bits see "x" bits)
    // bit map for tracklet quality xxxxyyyyy
    Int_t nCross(fTracklet[ip]->IsRowCross()?fTracklet[ip]->GetTBcross():0);
    spd = Double_t(fTracklet[ip]->GetTBoccupancy() | (nCross<<5));
    track->SetTRDmomentum(p, ip, &spd);
    track->SetTRDslice(fTracklet[ip]->GetdQdl(), ip, 0); // Set Summed dEdx into the first slice
  }
  // store PID probabilities
  track->SetTRDpid(fPID);
}
