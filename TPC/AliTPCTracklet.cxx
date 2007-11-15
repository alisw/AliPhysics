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
#include "AliTPCseed.h"
#include "AliESDVertex.h"

ClassImp(AliTPCTracklet)

AliTPCTracklet::AliTPCTracklet() 
  : fNClusters(0),fSector(-1),fOuter(0),fInner(0),fPrimary(0) {
  ////
  // The default constructor. It is intended to be used for I/O only.
  ////
}

AliTPCTracklet::AliTPCTracklet(const AliTPCseed *track,Int_t sector)
  : fNClusters(0),fSector(sector),fOuter(0),fInner(0),fPrimary(0) {
  ////
  // Contructor for a tracklet out of a track. Only clusters within a given 
  // sector are used.
  ///
  
  AliTPCseed *t=new AliTPCseed(*track);

  if (!t->Rotate(TMath::DegToRad()*(sector%18*20.+10.)-t->GetAlpha())) {
    delete t;
    return;
  }

  // fit from inner to outer row
  AliTPCseed *outerSeed=new AliTPCseed(*t);
  Int_t n=0;
  for (Int_t i=0;i<160;++i) {
    AliTPCclusterMI *c=t->GetClusterPointer(i);
    if (c&&c->GetDetector()==sector) {
      if (n==1)	{
	outerSeed->ResetCovariance(100.);
      }
      ++n;
      Double_t r[3]={c->GetX(),c->GetY(),c->GetZ()};
      Double_t cov[3]={0.1,0.,0.1}; //TODO: correct error parametrisation
      if (!outerSeed->PropagateTo(r[0])
	  || !static_cast<AliExternalTrackParam*>(outerSeed)->Update(&r[1],cov)) {
	delete outerSeed;
	outerSeed=0;
	break;
      }
    }
  }
  fNClusters=n;
  if (outerSeed)
    fOuter=new AliExternalTrackParam(*outerSeed);
  delete outerSeed;
  // fit from outer to inner rows
  AliTPCseed *innerSeed=new AliTPCseed(*t);
  n=0;
  for (Int_t i=159;i>=0;--i) {
    AliTPCclusterMI *c=t->GetClusterPointer(i);
    if (c&&c->GetDetector()==sector) {
      if (n==1)	{
	innerSeed->ResetCovariance(100.);
      }
      ++n;
      Double_t r[3]={c->GetX(),c->GetY(),c->GetZ()};
      Double_t cov[3]={0.1,0.,0.1};
      if (!innerSeed->PropagateTo(r[0])
	  || !static_cast<AliExternalTrackParam*>(innerSeed)->Update(&r[1],cov)) {
	delete innerSeed;
	innerSeed=0;
	break;
      }
    }
  }
  fNClusters=TMath::Max(fNClusters,n);
  if (innerSeed)
    fInner=new AliExternalTrackParam(*outerSeed);
  // propagate to the primary vertex
  if (innerSeed) {
    AliTPCseed *primarySeed=new AliTPCseed(*innerSeed);
    Double_t pos[]={0.,0.,0.};
    Double_t sigma[]={.1,.1,.1}; //TODO: is this correct?
    AliESDVertex vertex(pos,sigma);
    if (primarySeed->PropagateToVertex(&vertex))
      fPrimary=new AliExternalTrackParam(*primarySeed);
    delete primarySeed;
  }
  delete innerSeed;

  if (!fOuter&&!fInner)
    fNClusters=0;

  delete t;
}

AliTPCTracklet::AliTPCTracklet(const AliTPCTracklet &t)
  : fNClusters(t.fNClusters),fSector(t.fSector),fOuter(0),fInner(0),
    fPrimary(0) {
  ////
  // The copy constructor. You can copy tracklets! 
  ////

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
  fNClusters=t.fNClusters;
  fSector=t.fSector;
  if (this!=&t) {
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
  delete fOuter;
  delete fInner;
  delete fPrimary;
}

TObjArray AliTPCTracklet::CreateTracklets(const AliTPCseed *s,
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
    AliTPCclusterMI *c=s->GetClusterPointer(i);
    if (c)
      ++sectors[c->GetDetector()];
  }
  Int_t indices[72];
  TMath::Sort(72,sectors,indices);
  TObjArray tracklets;
  if (maxTracklets>72) maxTracklets=72; // just to protect against "users".
  for (Int_t i=0;i<maxTracklets&&sectors[indices[i]]>=minClusters;++i) {
    tracklets.Add(new AliTPCTracklet(s,indices[i]));
  }
  return tracklets;
}
