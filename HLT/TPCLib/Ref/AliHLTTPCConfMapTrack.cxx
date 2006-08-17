// @(#) $Id$

// Author: Anders Vestbo <mailto:vestbo@fi.uib.no>, Uli Frankenfeld <mailto:franken@fi.uib.no>
//*-- Copyright &copy ALICE HLT Group 

#include "AliHLTTPCStandardIncludes.h"

#include "AliHLTTPCRootTypes.h"
#include "AliHLTTPCLogging.h"
#include "AliHLTTPCVertex.h"
#include "AliHLTTPCConfMapPoint.h"
#include "AliHLTTPCConfMapFit.h"
#include "AliHLTTPCConfMapTrack.h"
#include "AliHLTTPCTransform.h"
#include "AliHLTTPC.h"

/** \class AliHLTTPCConfMapTrack
<pre>
//_____________________________________________________________
// AliHLTTPCConfMapTrack
//
// Track class for conformal mapper
</pre>
*/

ClassImp(AliHLTTPCConfMapTrack)


AliHLTTPCConfMapTrack::AliHLTTPCConfMapTrack()
{
  //Constructor

  fChiSq[0] = 0.;
  fChiSq[1] = 0.;
}

AliHLTTPCConfMapTrack::~AliHLTTPCConfMapTrack()
{

}

void AliHLTTPCConfMapTrack::DeleteCandidate()
{
  //Deletes this track by resetting all its parameters. Does not delete
  //the object itself.

  AliHLTTPCConfMapPoint *curHit = (AliHLTTPCConfMapPoint*)firstHit;
  AliHLTTPCConfMapPoint *nextHit;
  
  while(curHit != 0)
    {
      nextHit = (AliHLTTPCConfMapPoint*)curHit->nextTrackHit;
      curHit->nextTrackHit = 0;
      curHit = nextHit;
    }
  
  UInt_t *hit_numbers = GetHitNumbers();
  for(Int_t i=0; i<GetNHits(); i++)
    {
      //fHitNumbers[i] = 0;
      hit_numbers[i]=0;
    }
    
  SetRadius(0.);
  SetCenterX(0.);
  SetCenterY(0.);
  
  ComesFromMainVertex(false);

  SetNHits(0);
  SetCharge(0);
  fChiSq[0] = 0.;
  fChiSq[1] = 0.;
}


void AliHLTTPCConfMapTrack::SetProperties(Bool_t usage)
{
  //Set the hits to this track to 'usage'
  
  for(StartLoop(); LoopDone(); GetNextHit())
    {
      AliHLTTPCConfMapPoint *p = (AliHLTTPCConfMapPoint*)currentHit;
      p->SetUsage(usage);
    }
  return;
}

void AliHLTTPCConfMapTrack::Reset()
{
  //Resets the fit parameters of this track.

  //xy-plane
  s11Xy   = 0;
  s12Xy   = 0;
  s22Xy   = 0;
  g1Xy    = 0;
  g2Xy    = 0;
  fChiSq[0]  = 0.;
    
  //sz-plane
  s11Sz = 0;
  s12Sz = 0;
  s22Sz = 0;
  g1Sz  = 0;
  g2Sz  = 0;
  fChiSq[1] = 0; 
  SetLength(0);
  SetNHits(0);
}

void AliHLTTPCConfMapTrack::UpdateParam(AliHLTTPCConfMapPoint *thisHit)
{
  //Function to update fit parameters of track
  //Also, it updates the hit pointers.
  

  //Increment the number of hits assigned to this track:

  //fNHits++;
  Int_t nhits = GetNHits();
  nhits++;
  SetNHits(nhits); //SetNHits(nhits++);

  //Set the hit pointers:
  //if(fNHits == 1)
  if(GetNHits()==1)  
    firstHit = thisHit;
  else
    ((AliHLTTPCConfMapPoint*)lastHit)->nextTrackHit = thisHit;
  lastHit = thisHit;

  
  s11Xy = s11Xy + thisHit->GetXYWeight() ;
  s12Xy = s12Xy + thisHit->GetXYWeight() * thisHit->GetXprime() ;
  s22Xy = s22Xy + thisHit->GetXYWeight() * pow((thisHit->GetXprime()),2) ;
  g1Xy  = g1Xy  + thisHit->GetXYWeight() * thisHit->GetYprime() ;
  g2Xy  = g2Xy  + thisHit->GetXYWeight() * thisHit->GetXprime() * thisHit->GetYprime() ;
    
  ddXy  = s11Xy * s22Xy - pow((s12Xy),2) ;
  if ( ddXy != 0 ) 
    {
      a1Xy  = ( g1Xy * s22Xy - g2Xy * s12Xy ) / ddXy ;
      a2Xy  = ( g2Xy * s11Xy - g1Xy * s12Xy ) / ddXy ;
    }

  //     Now in the sz plane
  s11Sz = s11Sz + thisHit->GetZWeight() ;
  s12Sz = s12Sz + thisHit->GetZWeight() * thisHit->GetS() ;
  s22Sz = s22Sz + thisHit->GetZWeight() * thisHit->GetS() * thisHit->GetS() ;
  g1Sz  = g1Sz  + thisHit->GetZWeight() * thisHit->GetZ() ;
  g2Sz  = g2Sz  + thisHit->GetZWeight() * thisHit->GetS() * thisHit->GetZ() ;
  
        
  ddSz  = s11Sz * s22Sz -  s12Sz * s12Sz ;
  if ( ddSz != 0 ) {
    a1Sz  = ( g1Sz * s22Sz - g2Sz * s12Sz ) / ddSz ;
    a2Sz  = ( g2Sz * s11Sz - g1Sz * s12Sz ) / ddSz ;
  }
}


void AliHLTTPCConfMapTrack::Fill(AliHLTTPCVertex *vertex,Double_t max_Dca)
{
  //Fill track variables with or without fit.
  LOG(AliHLTTPCLog::kError,"AliHLTTPCConfMapTrack::Fill","TEST")<< "=========="<< ENDLOG;
  //fRadius = sqrt(a2Xy*a2Xy+1)/(2*fabs(a1Xy));
  Double_t radius = sqrt(a2Xy*a2Xy+1)/(2*fabs(a1Xy));
  LOG(AliHLTTPCLog::kError,"AliHLTTPCConfMapTrack::Fill","TEST")<< "====A====== r=" << GetRadius() << ENDLOG;
  SetRadius(radius);
  LOG(AliHLTTPCLog::kError,"AliHLTTPCConfMapTrack::Fill","TEST")<< "====B====== r=" << GetRadius() << ENDLOG;

  //fPt = (Double_t)(BFACT * AliHLTTPCTransform::GetBField() * fRadius);
  Double_t pt = (Double_t)(BFACT * AliHLTTPCTransform::GetBField() * GetRadius());
  SetPt(pt);

  if(GetPt() > max_Dca) //go for fit of helix in real space
    {
      AliHLTTPCConfMapFit *fit = new AliHLTTPCConfMapFit(this,vertex);
      fit->FitHelix();
      
      AliHLTTPCConfMapPoint *lHit = (AliHLTTPCConfMapPoint*)lastHit;
      AliHLTTPCConfMapPoint *fHit = (AliHLTTPCConfMapPoint*)firstHit;
      
      if(AliHLTTPC::IsTracksAtFirstPoint())
	{
	  //Z0 should not be overwritten here, since it is used as a starting
	  //point for the track swim in UpdateToFirstPoint. fFirstPointX and Y
	  //will be overwritten in UpdateToFirstPoint with the track fit 
	  //crossing point.
	  SetFirstPoint(lHit->GetX(),lHit->GetY(),GetZ0());
	  UpdateToFirstPoint();
	}
      else
	{
	  SetFirstPoint(vertex->GetX(),vertex->GetY(),vertex->GetZ());
	}
      
      SetLastPoint(fHit->GetX(),fHit->GetY(),fHit->GetZ());
      delete fit;
    }
  else if(GetPt() == 0)
    LOG(AliHLTTPCLog::kError,"AliHLTTPCConfMapTrack::Fill","Tracks")<<AliHLTTPCLog::kDec<<
      "Found track with Pt=0!!!"<<ENDLOG;
  else
    {
      LOG(AliHLTTPCLog::kError,"AliHLTTPCConfMapTrack::Fill","Tracks")<<AliHLTTPCLog::kDec<<
	"Track with pt<max_Dca :"<<GetPt()<<ENDLOG;
    }
}
/*
void AliHLTTPCConfMapTrack::UpdateToFirstPoint()
{
  //Update track parameters to the innermost point on the track.
  //Basically it justs calculates the intersection of the track, and a cylinder
  //with radius = r(innermost point). Then the parameters are updated to this point.
  //Should be called after the helixfit (in FillTracks).
  
  //AliHLTTPCConfMapPoint *lHit = (AliHLTTPCConfMapPoint*)fPoints->Last();
  AliHLTTPCConfMapPoint *lHit = (AliHLTTPCConfMapPoint*)lastHit;
  Double_t radius = sqrt(lHit->GetX()*lHit->GetX()+lHit->GetY()*lHit->GetY());
  
  //Get the track parameters
  
  Double_t tPhi0 = GetPsi() + GetCharge() * 0.5 * pi / abs(GetCharge()) ;
  Double_t x0    = GetR0() * cos(GetPhi0()) ;
  Double_t y0    = GetR0() * sin(GetPhi0()) ;
  Double_t rc    = fabs(GetPt()) / ( BFACT * AliHLTTPCTransform::GetBField() )  ;
  Double_t xc    = x0 - rc * cos(tPhi0) ;
  Double_t yc    = y0 - rc * sin(tPhi0) ;
  
  //Check helix and cylinder intersect
  
  Double_t fac1 = xc*xc + yc*yc ;
  Double_t sfac = sqrt( fac1 ) ;
    
  if ( fabs(sfac-rc) > radius || fabs(sfac+rc) < radius ) {
    LOG(AliHLTTPCLog::kError,"AliHLTTPCConfMapTrack::UpdateToLastPoint","Tracks")<<AliHLTTPCLog::kDec<<
      "Track does not intersect"<<ENDLOG;
    return;
  }
  
  //Find intersection
  
  Double_t fac2   = ( radius*radius + fac1 - rc*rc) / (2.00 * radius * sfac ) ;
  Double_t phi    = atan2(yc,xc) + GetCharge()*acos(fac2) ;
  Double_t td     = atan2(radius*sin(phi) - yc,radius*cos(phi) - xc) ;
  
  //Intersection in z
  
  if ( td < 0 ) td = td + 2. * pi ;
  Double_t deltat = fmod((-GetCharge()*td + GetCharge()*tPhi0),2*pi) ;
  if ( deltat < 0.      ) deltat += 2. * pi ;
  if ( deltat > 2.*pi ) deltat -= 2. * pi ;
  Double_t z = GetZ0() + rc * GetTgl() * deltat ;
 
  
  Double_t xExtra = radius * cos(phi) ;
  Double_t yExtra = radius * sin(phi) ;
  
  Double_t tPhi = atan2(yExtra-yc,xExtra-xc);
  
  //if ( tPhi < 0 ) tPhi += 2. * M_PI ;
  
  Double_t tPsi = tPhi - GetCharge() * 0.5 * pi / abs(GetCharge()) ;
  if ( tPsi > 2. * pi ) tPsi -= 2. * pi ;
  if ( tPsi < 0.        ) tPsi += 2. * pi ;
  
  //And finally, update the track parameters
  
  SetCenterX(xc);
  SetCenterY(yc);
  SetR0(radius);
  SetPhi0(phi);
  SetZ0(z);
  SetPsi(tPsi);
}
*/
Int_t AliHLTTPCConfMapTrack::GetMCLabel()
{
  //For evaluation study.
  //Returns the MCtrackID of the belonging clusters.
  //If MCLabel < 0, means that track is fake.

  return 0;
  /*
  Int_t num_of_clusters = GetNumberOfPoints();
  S *s=new S[num_of_clusters];
  Int_t i;
  for (i=0; i<num_of_clusters; i++) s[i].lab=s[i].max=0;
  
  Int_t lab=123456789;
  for (i=0; i<num_of_clusters; i++) {
    AliHLTTPCConfMapPoint *c=(AliHLTTPCConfMapPoint*)fPoints->UncheckedAt(i);
    lab=fabs(c->fMCTrackID[0]);
    Int_t j;
    for (j=0; j<num_of_clusters; j++)
      if (s[j].lab==lab || s[j].max==0) break;
    s[j].lab=lab;
    s[j].max++;
  }
  
  Int_t max=0;
  for (i=0; i<num_of_clusters; i++) 
    if (s[i].max>max) {max=s[i].max; lab=s[i].lab;}
    
  delete[] s;
  
  for (i=0; i<num_of_clusters; i++) {
    AliHLTTPCConfMapPoint *c=(AliHLTTPCConfMapPoint*)fPoints->UncheckedAt(i);
    if (fabs(c->fMCTrackID[1]) == lab ||
	fabs(c->fMCTrackID[2]) == lab ) max++;
  }
  
  //check if more than 10% of the clusters are incorrectly assigned (fake track):
  
  if (1.-Float_t(max)/num_of_clusters > 0.10) 
    {
      return -lab;
    }
  Int_t tail=Int_t(0.08*174);
  if (num_of_clusters < tail) return lab;
  
  max=0;
  for (i=1; i<=tail; i++) {
    AliHLTTPCConfMapPoint *c = (AliHLTTPCConfMapPoint*)fPoints->UncheckedAt(num_of_clusters-i);
    if (lab == fabs(c->fMCTrackID[0]) ||
	lab == fabs(c->fMCTrackID[1]) ||
	lab == fabs(c->fMCTrackID[2])) max++;
  }
  if (max < Int_t(0.5*tail)) 
    {
      //printf("Wrong innermost clusters\n");
      return -lab;
    }
  return lab;
  */
}

