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

////
// This class stores a tracklet (a track that lives only in a single TPC
// sector). Its objects can be constructed out of TPCseeds, that are
// holding the necessary cluster information.
////
////
//// 

#include "AliTPCTracklet.h"
#include "TObjArray.h"
#include "TLinearFitter.h"
#include "AliTPCseed.h"
#include "AliESDVertex.h"
#include "AliTracker.h"
#include "TTreeStream.h"
#include "TRandom3.h"
#include "TDecompChol.h"

#include <iostream>
using namespace std;

ClassImp(AliTPCTracklet)

const Double_t AliTPCTracklet::kB2C=0.299792458e-3;
Float_t  AliTPCTracklet::fgEdgeCutY=3;
Float_t  AliTPCTracklet::fgEdgeCutX=0;

AliTPCTracklet::AliTPCTracklet() 
  : fNClusters(0),fNStoredClusters(0),fClusters(0),fSector(-1),fOuter(0),
    fInner(0),fPrimary(0) {
  ////
  // The default constructor. It is intended to be used for I/O only.
  ////
}

AliTPCTracklet::AliTPCTracklet(const AliTPCseed *track,Int_t sector,
			       TrackType type,Bool_t storeClusters)
  : fNClusters(0),fNStoredClusters(0),fClusters(0),fSector(sector),fOuter(0),
    fInner(0),fPrimary(0) {
  ////
  // Contructor for a tracklet out of a track. Only clusters within a given 
  // sector are used.
  ///

  //TODO: only kalman works
  
  for (Int_t i=0;i<160;++i) {
    AliTPCclusterMI *c=track->GetClusterPointer(i);
    if (c && RejectCluster(c)) continue;
    if (c&&c->GetDetector()==sector)
      ++fNClusters;
  }

  if (storeClusters) {
    fClusters=new AliTPCclusterMI[fNClusters];
    for (Int_t i=0;i<160;++i) {
      AliTPCclusterMI *c=track->GetClusterPointer(i);
      if (c && RejectCluster(c)) continue;
      if (c&&c->GetDetector()==sector)
	fClusters[fNStoredClusters]=*c;
      ++fNStoredClusters;
    }
  }

  switch (type) {
  case kKalman:
    FitKalman(track,sector);
    break;
  case kLinear:
  case kQuadratic:
    FitLinear(track,sector,type);
    break;
  case kRiemann:
    FitRiemann(track,sector);
    break;
  }

}

AliTPCTracklet::AliTPCTracklet(const TObjArray &/*clusters*/,Int_t sector,
			       TrackType /*type*/,Bool_t /*storeClusters*/)
  : fNClusters(0),fNStoredClusters(0),fClusters(0),fSector(sector),fOuter(0),
    fInner(0),fPrimary(0) {
  //TODO: write it!
}

AliTPCTracklet::AliTPCTracklet(const AliTPCTracklet &t)
  : TObject(t),fNClusters(t.fNClusters),fNStoredClusters(t.fNStoredClusters),fClusters(0),
    fSector(t.fSector),fOuter(0),fInner(0),
    fPrimary(0) {
  ////
  // The copy constructor. You can copy tracklets! 
  ////

  if (t.fClusters) {
    fClusters=new AliTPCclusterMI[t.fNStoredClusters];
    for (int i=0;i<t.fNStoredClusters;++i)
      fClusters[i]=t.fClusters[i];
  }
  if (t.fOuter)
    fOuter=new AliExternalTrackParam(*t.fOuter);
  if (t.fInner)
    fInner=new AliExternalTrackParam(*t.fInner);
  if (t.fPrimary)
    fPrimary=new AliExternalTrackParam(*t.fPrimary);
}

AliTPCTracklet& AliTPCTracklet::operator=(const AliTPCTracklet &t) {
  ////
  // The assignment constructor. You can assign tracklets!
  ////
  if (this!=&t) {
    fNClusters=t.fNClusters;
    fNStoredClusters=fNStoredClusters;
    delete fClusters;
    if (t.fClusters) {
      fClusters=new AliTPCclusterMI[t.fNStoredClusters];
      for (int i=0;i<t.fNStoredClusters;++i)
	fClusters[i]=t.fClusters[i];
    }
    else
      fClusters=0;
    fSector=t.fSector;
    if (t.fOuter) {
      if (fOuter)
	*fOuter=*t.fOuter;
      else
	fOuter=new AliExternalTrackParam(*t.fOuter);
    }
    else {
      delete fOuter;
      fOuter=0;
    }

    if (t.fInner) {
      if (fInner)
	*fInner=*t.fInner;
      else
	fInner=new AliExternalTrackParam(*t.fInner);
    }
    else {
      delete fInner;
      fInner=0;
    }

    if (t.fPrimary) {
      if (fPrimary)
	*fPrimary=*t.fPrimary;
      else
	fPrimary=new AliExternalTrackParam(*t.fPrimary);
    }
    else {
      delete fPrimary;
      fPrimary=0;
    }
  }
  return *this;
}

AliTPCTracklet::~AliTPCTracklet() {
  //
  // The destructor. Yes, you can even destruct tracklets.
  //
  delete fClusters;
  delete fOuter;
  delete fInner;
  delete fPrimary;
}





void AliTPCTracklet::FitKalman(const AliTPCseed *seed,Int_t sector) {
  //
  // Fit using Kalman filter
  //
  AliTPCseed *track=new AliTPCseed(*seed);
  if (!track->Rotate(TMath::DegToRad()*(sector%18*20.+10.)-track->GetAlpha())) {
    delete track;
    return;
  }
  // fit from inner to outer row
  Double_t covar[15];
  for (Int_t i=0;i<15;i++) covar[i]=0;
  covar[0]=1.*1.;
  covar[2]=1.*1.;
  covar[5]=1.*1./(64.*64.);
  covar[9]=1.*1./(64.*64.);
  covar[14]=0;  // keep pt
  Float_t xmin=1000, xmax=-10000;
  Int_t imin=158, imax=0;
  for (Int_t i=0;i<160;i++) {
    AliTPCclusterMI *c=track->GetClusterPointer(i);
    if (!c) continue;
    if (c->GetDetector()!=sector)  continue;
    if (c->GetX()<xmin) xmin=c->GetX();
    if (c->GetX()>xmax) xmax=c->GetX();
    if (i<imin) imin=i;
    if (i>imax) imax=i;
  }
  if(imax-imin<10) {
    delete track;
    return;
  }

  for (Float_t x=track->GetX(); x<xmin; x++) track->PropagateTo(x);
  track->AddCovariance(covar);
  //
  AliExternalTrackParam paramIn;
  AliExternalTrackParam paramOut;
  Bool_t isOK=kTRUE;
  //
  //
  //
  for (Int_t i=imin; i<=imax; i++){
    AliTPCclusterMI *c=track->GetClusterPointer(i);
    if (!c) continue;
    Double_t r[3]={c->GetX(),c->GetY(),c->GetZ()};
    Double_t cov[3]={0.01,0.,0.01}; //TODO: correct error parametrisation
    AliTPCseed::GetError(c, track,cov[0],cov[2]);
    cov[0]*=cov[0];
    cov[2]*=cov[2];
    if (!track->PropagateTo(r[0])) {
      isOK=kFALSE;
      break;
    }
    if (RejectCluster(c,track)) continue;
    if ( !((static_cast<AliExternalTrackParam*>(track)->Update(&r[1],cov)))) isOK=kFALSE;
  }
  if (!isOK) { delete track; return;}
  track->AddCovariance(covar);
  //
  //
  //
  for (Int_t i=imax; i>=imin; i--){
    AliTPCclusterMI *c=track->GetClusterPointer(i);
    if (!c) continue;
    Double_t r[3]={c->GetX(),c->GetY(),c->GetZ()};
    Double_t cov[3]={0.01,0.,0.01}; 
    AliTPCseed::GetError(c, track,cov[0],cov[2]);
    cov[0]*=cov[0];
    cov[2]*=cov[2];
    if (!track->PropagateTo(r[0])) {
      isOK=kFALSE;
      break;
    }
    if (RejectCluster(c,track)) continue;
    if ( !((static_cast<AliExternalTrackParam*>(track)->Update(&r[1],cov)))) isOK=kFALSE;
  }
  if (!isOK) { delete track; return;}
  paramIn = *track;
  track->AddCovariance(covar);
  //
  //
  for (Int_t i=imin; i<=imax; i++){
    AliTPCclusterMI *c=track->GetClusterPointer(i);
    if (!c) continue;
    Double_t r[3]={c->GetX(),c->GetY(),c->GetZ()};
    Double_t cov[3]={0.01,0.,0.01}; 
    AliTPCseed::GetError(c, track,cov[0],cov[2]);
    cov[0]*=cov[0];
    cov[2]*=cov[2];
    if (!track->PropagateTo(r[0])) {
      isOK=kFALSE;
      break;
    }
    if (RejectCluster(c,track)) continue;
    if ( !((static_cast<AliExternalTrackParam*>(track)->Update(&r[1],cov)))) isOK=kFALSE;
  }
  if (!isOK) { delete track; return;}
  paramOut=*track;
  //
  //
  //
  fOuter=new AliExternalTrackParam(paramOut);
  fInner=new AliExternalTrackParam(paramIn);
  //
  delete track;
}




void AliTPCTracklet::FitLinear(const AliTPCseed *track,Int_t sector,
			       TrackType type) {
  TLinearFitter fy(1);
  TLinearFitter fz(1);
  fy.StoreData(kFALSE);
  fz.StoreData(kFALSE);
  switch (type) {
  case kLinear:
    fy.SetFormula("1 ++ x");
    fz.SetFormula("1 ++ x");
    break;
  case kQuadratic:
    fy.SetFormula("1 ++ x ++ x*x");
    fz.SetFormula("1 ++ x");
    break;
  case kKalman:
  case kRiemann:
    break;
  }
  Double_t xmax=-1.;
  Double_t xmin=1000.;
  for (Int_t i=0;i<160;++i) {
    AliTPCclusterMI *c=track->GetClusterPointer(i);
    if (c && RejectCluster(c)) continue;
    if (c&&c->GetDetector()==sector) {
      Double_t x=c->GetX();
      fy.AddPoint(&x,c->GetY());
      fz.AddPoint(&x,c->GetZ());
      xmax=TMath::Max(xmax,x);
      xmin=TMath::Min(xmin,x);
    }
  }
  fy.Eval();
  fz.Eval();
  Double_t a[3]={fy.GetParameter(0),
		 fy.GetParameter(1),
		 type==kQuadratic?fy.GetParameter(2):0.};
  Double_t ca[6]={fy.GetCovarianceMatrixElement(0,0),
		  fy.GetCovarianceMatrixElement(1,0),
		  fy.GetCovarianceMatrixElement(1,1),
		  type==kQuadratic?fy.GetCovarianceMatrixElement(2,0):0.,
		  type==kQuadratic?fy.GetCovarianceMatrixElement(2,1):0.,
		  type==kQuadratic?fy.GetCovarianceMatrixElement(2,2):0.};
  for (int i=0;i<6;++i) ca[i]*=fy.GetChisquare()/fNClusters;
  Double_t b[2]={fz.GetParameter(0),
		 fz.GetParameter(1)};
  Double_t cb[3]={fz.GetCovarianceMatrixElement(0,0),
		  fz.GetCovarianceMatrixElement(1,0),
		  fz.GetCovarianceMatrixElement(1,1)};
  for (int i=0;i<3;++i) cb[i]*=fz.GetChisquare()/fNClusters;
  Double_t p[5];
  Double_t c[15];
  Double_t alpha=track->GetAlpha();
  Quadratic2Helix(a,ca,b,cb,0.,p,c);
  fPrimary=new AliExternalTrackParam(0.,alpha,p,c);
  Quadratic2Helix(a,ca,b,cb,xmin,p,c);
  fInner=new AliExternalTrackParam(xmin,alpha,p,c);
  Quadratic2Helix(a,ca,b,cb,xmax,p,c);
  fOuter=new AliExternalTrackParam(xmax,alpha,p,c);
}
  
void AliTPCTracklet::Quadratic2Helix(Double_t *a,Double_t *ca,
				     Double_t *b,Double_t *cb,
				     Double_t x0,
				     Double_t *p,Double_t *c) {
  // y(x)=a[0]+a[1]*x+a[2]*x^2
  // z(x)=b[0]+b[1]*x
  // parametrises the corosponding helix at x0

  // get the polynoms at x0
  Double_t a0=x0*x0*a[2] + x0*a[1] + a[0];
  Double_t a1=2.*x0*a[2] +     a[1];
  Double_t a2=      a[2];
  Double_t ca00=ca[0]+x0*(2.*ca[1]+x0*(ca[2]+2.*ca[3]+x0*(2.*ca[4]+x0*ca[5])));
  Double_t ca10=ca[1]+x0*(ca[2]+2.*ca[3]+x0*(3.*ca[4]+x0*2.*ca[5]));
  Double_t ca11=ca[2]+x0*4.*(ca[4]+x0*ca[5]);
  Double_t ca20=ca[3]+x0*(ca[4]+x0*ca[5]);
  Double_t ca21=ca[3]+x0*2.*ca[5];
  Double_t ca22=ca[5];

  Double_t b0=x0*b[1] + b[0];
  Double_t b1=   b[1];
  Double_t cb00=cb[0]+x0*(2.*cb[1]+x0*cb[2]);
  Double_t cb10=cb[1]+x0*cb[2];
  Double_t cb11=cb[2];

  // transform to helix parameters
  Double_t f   =1.+a1*a1;
  Double_t f2  =f*f;
  Double_t fi  =1./f; 
  Double_t fi12=TMath::Sqrt(fi);
  Double_t fi32=fi*fi12;
  Double_t fi2 =fi*fi;
  Double_t fi52=fi2*fi12;
  Double_t fi3 =fi2*fi;
  Double_t fi5 =fi2*fi3;
  
  Double_t xyz[3]={0.}; // TODO...
  Double_t fc=1./(GetBz(xyz)*kB2C);

  p[0]=a0;            // y0
  p[1]=b0;            // z0
  p[2]=a1*fi12;       // snp
  p[3]=b1;            // tgl
  p[4]=2.*a2*fi32*fc; // 1/pt

  c[0] =ca00;      //  y0-y0
  c[1] =0.;        //  z0-y0
  c[2] =cb00;      //  z0-z0
  c[3] =ca10*fi32; // snp-y0
  c[4] =0.;        // snp-z0
  c[5] =ca11*fi3;  // snp-snp
  c[6] =0.;        // tgl-y0
  c[7] =cb10;      // tgl-z0
  c[8] =0.;        // tgl-snp
  c[9] =cb11;      // tgl-tgl
  c[10]=2.*(-3.*a1*a2*ca10+f*ca20)*fi3*fc;  // 1/pt-y0
  c[11]=0.;                                 // 1/pt-z0
  c[12]=2.*(-3.*a1*a2*ca11+f*ca21)*fi52*fc; // 1/pt-snp
  c[13]=0.;                                 // 1/pt-tgl
  c[14]=(-12.*a1*a2*(-3.*a1*a2*ca11+2.*f*ca21)+4.*f2*ca22)*fi5
    *fc*fc;        // 1/pt-1/pt
}


void AliTPCTracklet::FitRiemann(const AliTPCseed *track,Int_t sector) {
  TLinearFitter fy(2);
  fy.StoreData(kFALSE);
  fy.SetFormula("hyp2");
  Double_t xmax=-1.;
  Double_t xmin=1000.;
  for (Int_t i=0;i<160;++i) {
    AliTPCclusterMI *c=track->GetClusterPointer(i);
    if (c && RejectCluster(c)) continue;
    if (c&&c->GetDetector()==sector) {
      Double_t x=c->GetX();
      Double_t y=c->GetY();
      Double_t xy[2]={x,y};
      Double_t r=x*x+y*y;
      Double_t errx=1.,erry=1.;//TODO!
      Double_t err=TMath::Sqrt(4.*x*x*errx+4.*y*y*erry);
      err=1.;
      fy.AddPoint(xy,r,err);
      xmax=TMath::Max(xmax,x);
      xmin=TMath::Min(xmin,x);
    }
  }
  fy.Eval();
  Double_t a[3]={fy.GetParameter(0),
		 fy.GetParameter(1),
		 fy.GetParameter(2)};
  Double_t ca[6]={fy.GetCovarianceMatrixElement(0,0),
		  fy.GetCovarianceMatrixElement(1,0),
		  fy.GetCovarianceMatrixElement(1,1),
		  fy.GetCovarianceMatrixElement(2,0),
		  fy.GetCovarianceMatrixElement(2,1),
		  fy.GetCovarianceMatrixElement(2,2)};

  TLinearFitter fz(1);
  fz.StoreData(kFALSE);
  fz.SetFormula("hyp1");
  Double_t R=.5*TMath::Sqrt(4.*a[0]+a[1]*a[1]+a[2]*a[2]);
  Double_t oldx=0.;
  Double_t oldy=R;
  Double_t phi=0.;
  for (Int_t i=0;i<160;++i) {
    AliTPCclusterMI *c=track->GetClusterPointer(i);
    if (c && RejectCluster(c)) continue;
    if (c&&c->GetDetector()==sector) {
      Double_t x=c->GetX();
      Double_t y=c->GetY();
      Double_t dx=x-oldx;
      Double_t dy=y-oldy;
      phi+=2.*TMath::Abs(TMath::ATan2(.5*TMath::Sqrt(dx*dx+dy*dy),R));
      Double_t err=1.;
      fz.AddPoint(&phi,c->GetZ(),err);
      oldx=x;
      oldy=y;
    }
  }
  fz.Eval();
  Double_t b[2]={fz.GetParameter(0),
		 fz.GetParameter(1)};
  Double_t cb[3]={fz.GetCovarianceMatrixElement(0,0),
		  fz.GetCovarianceMatrixElement(1,0),
		  fz.GetCovarianceMatrixElement(1,1)};

  Double_t p[5];
  Double_t c[15];
  Double_t alpha=track->GetAlpha();
  if (Riemann2Helix(a,ca,b,cb,0.,p,c))
    fPrimary=new AliExternalTrackParam(0.,alpha,p,c);
  if (Riemann2Helix(a,ca,b,cb,xmin,p,c))
    fInner=new AliExternalTrackParam(xmin,alpha,p,c);
  if (Riemann2Helix(a,ca,b,cb,xmax,p,c))
    fOuter=new AliExternalTrackParam(xmax,alpha,p,c);
}

Bool_t AliTPCTracklet::Riemann2Helix(Double_t *a,Double_t */*ca*/,
				     Double_t *b,Double_t */*cb*/,
				     Double_t x0,
				     Double_t *p,Double_t *c) {
  //TODO: signs!

  Double_t xr0=.5*a[1];
  Double_t yr0=.5*a[2];
  Double_t R=.5*TMath::Sqrt(4.*a[0]+a[1]*a[1]+a[2]*a[2]);
  Double_t dx=x0-xr0;
  if (dx*dx>=R*R) return kFALSE;
  Double_t dy=TMath::Sqrt((R-dx)*(R+dx)); //sign!!
  if (TMath::Abs(yr0+dy)>TMath::Abs(yr0-dy))
    dy=-dy;
  Double_t y0=yr0+dy; 
  Double_t tgp=-dx/dy; //TODO: dy!=0
  Double_t z0=b[0]+TMath::ATan(tgp)*b[1];
  Double_t xyz[3]={x0,y0,z0};
  Double_t fc=1./(GetBz(xyz)*kB2C);
  fc=1;
  p[0]=y0;  // y0
  p[1]=z0; // z0
  p[2]=tgp/TMath::Sqrt(1.+tgp*tgp); // snp
  p[3]=b[1];       // tgl
  p[4]=1./R*fc;    // 1/pt

  c[0] =0.;      //  y0-y0
  c[1] =0.;        //  z0-y0
  c[2] =0.;      //  z0-z0
  c[3] =0.; // snp-y0
  c[4] =0.;        // snp-z0
  c[5] =0.;  // snp-snp
  c[6] =0.;        // tgl-y0
  c[7] =0.;      // tgl-z0
  c[8] =0.;        // tgl-snp
  c[9] =0.;      // tgl-tgl
  c[10]=0.;  // 1/pt-y0
  c[11]=0.;                                 // 1/pt-z0
  c[12]=0.; // 1/pt-snp
  c[13]=0.;                                 // 1/pt-tgl
  c[14]=0.;        // 1/pt-1/pt

  return kTRUE;
}

TObjArray AliTPCTracklet::CreateTracklets(const AliTPCseed *track,
					  TrackType type,
					  Bool_t storeClusters,
					  Int_t minClusters,
					  Int_t maxTracklets) {
// The tracklet factory: It creates several tracklets out of a track. They
// are created for sectors that fullfill the constraint of having enough
// clusters inside. Futhermore you can specify the maximum amount of
// tracklets that are to be created.
// The tracklets appear in a sorted fashion, beginning with those having the
// most clusters.

  Int_t sectors[72]={0};
  for (Int_t i=0;i<160;++i) {
    AliTPCclusterMI *c=track->GetClusterPointer(i);
    if (c && RejectCluster(c)) continue;
    if (c)
      ++sectors[c->GetDetector()];
  }
  Int_t indices[72];
  TMath::Sort(72,sectors,indices);
  TObjArray tracklets;
  if (maxTracklets>72) maxTracklets=72; // just to protect against "users".
  for (Int_t i=0;i<maxTracklets&&sectors[indices[i]]>=minClusters;++i) {
    tracklets.Add(new AliTPCTracklet(track,indices[i],type,storeClusters));
  }
  return tracklets;
}

TObjArray AliTPCTracklet::CreateTracklets(const TObjArray &/*clusters*/,
					  TrackType /*type*/,
					  Bool_t /*storeClusters*/,
					  Int_t /*minClusters*/,
					  Int_t /*maxTracklets*/) {
  // TODO!

  TObjArray tracklets;
  return tracklets;
}

Bool_t AliTPCTracklet::PropagateToMeanX(const AliTPCTracklet &t1,
					const AliTPCTracklet &t2,
					AliExternalTrackParam *&t1m,
					AliExternalTrackParam *&t2m) {
  // This function propagates two Tracklets to a common x-coordinate. This
  // x is dermined as the one that is in the middle of the two tracklets (they
  // are assumed to live on two distinct x-intervalls).
  // The inner parametrisation of the outer Tracklet and the outer 
  // parametrisation of the inner Tracklet are used and propagated to this
  // common x. This result is saved not inside the Tracklets but two new
  // ExternalTrackParams are created (that means you might want to delete
  // them).
  // In the case that the alpha angles of the Tracklets differ both angles
  // are tried out for this propagation.
  // In case of any failure kFALSE is returned, no AliExternalTrackParam
  // is created und the pointers are set to 0.

  if (t1.GetInner() && t1.GetOuter() && 
      t2.GetInner() && t2.GetOuter()) {
    if (t1.GetOuter()->GetX()<t2.GetInner()->GetX()) {
      t1m=new AliExternalTrackParam(*t1.GetOuter());
      t2m=new AliExternalTrackParam(*t2.GetInner());
    }
    else {
      t1m=new AliExternalTrackParam(*t1.GetInner());
      t2m=new AliExternalTrackParam(*t2.GetOuter());
    }
    Double_t mx=.5*(t1m->GetX()+t2m->GetX());
    //Double_t b1,b2;
    Double_t xyz[3];
    t1m->GetXYZ(xyz);
    //b1=GetBz(xyz);
    Double_t b1[3]; AliTracker::GetBxByBz(xyz,b1);
    t2m->GetXYZ(xyz);
    //b2=GetBz(xyz);
    Double_t b2[3]; AliTracker::GetBxByBz(xyz,b2);
    if (t1m->Rotate(t2m->GetAlpha()) 
	//&& t1m->PropagateTo(mx,b1) 
	//&& t2m->PropagateTo(mx,b2));
	&& t1m->PropagateToBxByBz(mx,b1) 
	&& t2m->PropagateToBxByBz(mx,b2));
    else
      if (t2m->Rotate(t1m->GetAlpha())
	  //&& t1m->PropagateTo(mx,b1) 
	  //&& t2m->PropagateTo(mx,b2));
	  && t1m->PropagateToBxByBz(mx,b1) 
	  && t2m->PropagateToBxByBz(mx,b2));
      else {
	delete t1m;
	delete t2m;
	t1m=t2m=0;
      }
  }
  else {
    t1m=t2m=0;
  }
  return t1m&&t2m;
}

double AliTPCTracklet::GetBz(Double_t *xyz) 
{
  return AliTracker::GetBz(xyz);
}

void AliTPCTracklet::RandomND(Int_t ndim,const Double_t *p,const Double_t *c,
			      Double_t *x) {
  // This function generates a n-dimensional random variable x with mean
  // p and covariance c.
  // That is done using the cholesky decomposition of the covariance matrix,
  // Begin_Latex C=U^{t} U End_Latex, with Begin_Latex U End_Latex being an
  // upper triangular matrix. Given a vector v of iid gausian random variables
  // with variance 1 one obtains the asked result as: Begin_Latex x=U^t v 
  // End_Latex.
  // c is expected to be in a lower triangular format:
  // c[0]
  // c[1] c[2]
  // c[3] c[4] c[5]
  // etc.
  static TRandom3 random;
  Double_t *c2= new Double_t[ndim*ndim];
  Int_t k=0;
  for (Int_t i=0;i<ndim;++i)
    for (Int_t j=0;j<=i;++j)
      c2[i*ndim+j]=c2[j*ndim+i]=c[k++];
  TMatrixDSym cm(ndim,c2);
  delete[] c2;
  TDecompChol chol(cm);
  chol.Decompose();
  const TVectorD pv(ndim);
  const_cast<TVectorD*>(&pv)->Use(ndim,const_cast<Double_t*>(p));
  TVectorD xv(ndim);
  xv.Use(ndim,x);
  for (Int_t i=0;i<ndim;++i)
    xv[i]=random.Gaus();
  TMatrixD L=chol.GetU();
  L.T();
  xv=L*xv+pv;
}

TEllipse AliTPCTracklet::ErrorEllipse(Double_t x,Double_t y,
				      Double_t sx,Double_t sy,Double_t sxy) {
  /* Begin_Latex
     r_{1,2}=1/2 (a+c#pm#sqrt{(a-c)^{2}+(2b)^{2}})
  End_Latex */
  Double_t det1=1./(sx*sy-sxy*sxy);
  Double_t a=sy*det1;
  Double_t b=-sxy*det1;
  Double_t c=sx*det1;
  Double_t d=c-a;
  Double_t s=TMath::Sqrt(d*d+4.*b*b);
  Double_t r1=TMath::Sqrt(.5*(a+c-s));
  Double_t r2=TMath::Sqrt(.5*(a+c+s));
  Double_t alpha=.5*TMath::ATan2(2.*b,d);
  return TEllipse(x,y,r1,r2,0.,360.,alpha*TMath::RadToDeg());
}

void AliTPCTracklet::Test(const char* filename) {
  /*
    aliroot
    AliTPCTracklet::Test("");
    TFile f("AliTPCTrackletDebug.root");
    TTree *t=f.Get("AliTPCTrackletDebug");
    t->Draw("p0:p4");
    TEllipse e=AliTPCTracklet::ErrorEllipse(0.,0.,4.,1.,1.8);
    e.Draw();
 */
  TTreeSRedirector ds(filename);
  Double_t p[5]={0.};
  Double_t c[15]={4.,
		  0.,4.,
		  0.,0.,9.,
		  0.,0.,0.,16.,
		  1.8,0.,0.,0.,1.};
  for (Int_t i=0;i<10000;++i) {
    Double_t x[5];
    RandomND(5,p,c,x);
    ds<<"AliTPCTrackletDebug"
      <<"p0="<<x[0]
      <<"p1="<<x[1]
      <<"p2="<<x[2]
      <<"p3="<<x[3]
      <<"p4="<<x[4]
      <<"\n";
  }

  /*
  Double_t b;
  Double_t x=0.;
  Double_t alpha=0.;
  Double_t param[5]={0.};
  Double_t covar[15]={1.,
		      0.,1.,
		      0.,0.,1.,
		      0.,0.,0.,1.,
		      0.,0.,0.,0.,1.};
  AliExternalTrackParam track(x,alpha,param,covar);

  

  for (Int_t i=0;i<points.GetNPoints();++i) {
    Double_t x=0.;
    Double_t alpha=0.;
    Double_t param[5]={0.};
    Double_t covar[15]={1.,
			0.,1.,
			0.,0.,1.,
			0.,0.,0.,1.,
			0.,0.,0.,0.,1.};
    AliExternalTrackParam track(x,alpha,param,covar);
    for (x=90.;x<250.;x+=1.) {
      track.PropagateTo(x,b);
      AliTPCclusterMI c();
    }
  }
  */
}


Bool_t AliTPCTracklet::RejectCluster(AliTPCclusterMI* cl, AliExternalTrackParam * param){
  //
  // check the acceptance of cluster
  // Cut on edge effects
  //
  Bool_t isReject = kFALSE;
  Float_t edgeY = cl->GetX()*TMath::Tan(TMath::Pi()/18);
  Float_t dist  = edgeY - TMath::Abs(cl->GetY());
  if (param)  dist  = edgeY - TMath::Abs(param->GetY());
  if (dist<fgEdgeCutY) isReject=kTRUE;
  if (cl->GetType()<0) isReject=kTRUE;
  return isReject;
}
