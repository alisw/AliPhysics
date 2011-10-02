// $Id$

//**************************************************************************
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//*                                                                        *
//* Primary Authors: Matthias Richter <Matthias.Richter@ift.uib.no>        *
//*                  for The ALICE HLT Project.                            *
//*                                                                        *
//* Permission to use, copy, modify and distribute this software and its   *
//* documentation strictly for non-commercial purposes is hereby granted   *
//* without fee, provided that the above copyright notice appears in all   *
//* copies and that both the copyright notice and this permission notice   *
//* appear in the supporting documentation. The authors make no claims     *
//* about the suitability of this software for any purpose. It is          *
//* provided "as is" without express or implied warranty.                  *
//**************************************************************************

/** @file   AliHLTGlobalBarrelTrack.cxx
    @author Matthias Richter
    @date   2009-06-24
    @brief  An AliKalmanTrack implementation for global HLT barrel tracks.
*/

// see header file for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

#include <cassert>
#include <memory>
#include <iostream>
#include "AliHLTGlobalBarrelTrack.h"
#include "AliHLTSpacePointContainer.h"
#include "AliHLTTrackGeometry.h"
#include "AliHLTMisc.h"
#include "TClonesArray.h"
#include "TFile.h"
#include "TArrayC.h"
#include "TMath.h"
#include "TMarker.h"
#include "TArc.h"

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTGlobalBarrelTrack)

AliHLTGlobalBarrelTrack::AliHLTGlobalBarrelTrack()
: AliKalmanTrack()
  , fPoints()
  , fLastX(0.0)
  , fLastY(0.0)
  , fTrackID(-1)
  , fHelixRadius(0.0)
  , fHelixCenterX(0.0)
  , fHelixCenterY(0.0)
  , fSpacePoints(NULL)
  , fTrackPoints(NULL)
{
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt
}

AliHLTGlobalBarrelTrack::AliHLTGlobalBarrelTrack(const AliHLTGlobalBarrelTrack& t)
  : AliKalmanTrack(t)
  , fPoints()
  , fLastX(t.GetLastPointX())
  , fLastY(t.GetLastPointY())
  , fTrackID(t.TrackID())
  , fHelixRadius(t.fHelixRadius)
  , fHelixCenterX(t.fHelixCenterX)
  , fHelixCenterY(t.fHelixCenterY)
  , fSpacePoints(NULL)
  , fTrackPoints(NULL)
{
  // see header file for class documentation
  fPoints.assign(t.fPoints.begin(), t.fPoints.end());
}

AliHLTGlobalBarrelTrack::AliHLTGlobalBarrelTrack(const AliHLTExternalTrackParam& p )
  : AliKalmanTrack()
  , fPoints()
  , fLastX(p.fLastX)
  , fLastY(p.fLastY)
  , fTrackID(p.fTrackID)
  , fHelixRadius(0.0)
  , fHelixCenterX(0.0)
  , fHelixCenterY(0.0)
  , fSpacePoints(NULL)
  , fTrackPoints(NULL)
{
  // see header file for class documentation

  // the 5 track parameters are named in the AliHLTExternalTrackParam
  // while AliExternalTrackParam just uses an array[5]
  // the members have the same order, fY is the first one
  Set(p.fX, p.fAlpha, &p.fY, p.fC);
  SetPoints(p.fPointIDs, p.fNPoints);
  SetNumberOfClusters(p.fNPoints);
  //SetIntegratedLength(GetPathLengthTo( GetLastPointX(), b);
}

AliHLTGlobalBarrelTrack::AliHLTGlobalBarrelTrack(const AliExternalTrackParam& p )
  : AliKalmanTrack()
  , fPoints()
  , fLastX(0)
  , fLastY(0)
  , fTrackID(0)
  , fHelixRadius(0.0)
  , fHelixCenterX(0.0)
  , fHelixCenterY(0.0)
  , fSpacePoints(NULL)
  , fTrackPoints(NULL)
{
  // see header file for class documentation
  *(dynamic_cast<AliExternalTrackParam*>(this))=p;
}

AliHLTGlobalBarrelTrack::~AliHLTGlobalBarrelTrack()
{
  // see header file for class documentation
}


Double_t AliHLTGlobalBarrelTrack::GetPathLengthTo( Double_t x, Double_t b ) const
{
  // calculate the trajectory length for dx step
    
  Double_t dx = x - GetX();
  Double_t ey = GetSnp();
  if( TMath::Abs( ey )>=kAlmost1 ) return 0;

  Double_t ex = TMath::Sqrt(1-ey*ey);
  Double_t k  = GetC(b);

  Double_t ey1 = k * dx + ey;

  // check for intersection with X=x

  if ( TMath::Abs( ey1 ) >= kAlmost1  ) return 0;

  Double_t ex1 = TMath::Sqrt(1-ey1*ey1);

  Double_t ss = ey + ey1;
  Double_t cc = ex + ex1;

  if ( TMath::Abs( cc ) < 1.e-4  ) return 0;

  Double_t tg = ss / cc; 
  Double_t dl = dx * TMath::Sqrt( 1 + tg * tg );
  Double_t dSin = dl * k / 2;
  if ( dSin > 1 ) dSin = 1;
  if ( dSin < -1 ) dSin = -1;
  Double_t dS = ( TMath::Abs( k ) > 1.e-4 )  ? ( 2 * TMath::ASin( dSin ) / k ) : dl;

  return dS*TMath::Sqrt(1 + GetTgl()*GetTgl() );
}




int AliHLTGlobalBarrelTrack::ConvertTrackDataArray(const AliHLTTracksData* pTracks, unsigned sizeInByte, vector<AliHLTGlobalBarrelTrack> &tgtArray )
{
  // see header file for class documentation
  int iResult=0;
  if (!pTracks || sizeInByte<sizeof(AliHLTTracksData) || pTracks->fCount==0) return 0;

  const AliHLTUInt8_t* pEnd=reinterpret_cast<const AliHLTUInt8_t*>(pTracks);
  pEnd+=sizeInByte;

  tgtArray.resize(pTracks->fCount + tgtArray.size());
  const AliHLTUInt8_t* pCurrent=reinterpret_cast<const AliHLTUInt8_t*>(pTracks->fTracklets);
  for (unsigned i=0; i<pTracks->fCount; i++) {
    if (pCurrent+sizeof(AliHLTExternalTrackParam)>pEnd) {
      iResult=-EINVAL; break;
    }
    const AliHLTExternalTrackParam* track=reinterpret_cast<const AliHLTExternalTrackParam*>(pCurrent);
    if (pCurrent+sizeof(AliHLTExternalTrackParam)+track->fNPoints*sizeof(UInt_t)>pEnd) {
      iResult=-EINVAL; break;
    }
    tgtArray[i]=*track;
    pCurrent+=sizeof(AliHLTExternalTrackParam)+track->fNPoints*sizeof(UInt_t);
  }
  if (iResult<0) tgtArray.clear();
  else iResult=tgtArray.size();
  return iResult;
}

UInt_t AliHLTGlobalBarrelTrack::GetNumberOfPoints() const
{
  // see header file for class documentation
  return fPoints.size();
}

const UInt_t* AliHLTGlobalBarrelTrack::GetPoints() const
{
  // see header file for class documentation
  if (fPoints.size()==0) return NULL;
  return &fPoints[0];
}

int AliHLTGlobalBarrelTrack::SetPoints(const UInt_t* pArray, UInt_t arraySize)
{
  // see header file for class documentation
  if (!pArray || arraySize==0) return 0;
  fPoints.resize(arraySize);
  for (unsigned i=0; i<arraySize; i++) fPoints[i]=pArray[i];
  return fPoints.size();
}

int AliHLTGlobalBarrelTrack::CalculateHelixParams()
{
  // calculate radius and center of the helix
  // using the global magnetic field
  return CalculateHelixParams(AliHLTMisc::Instance().GetBz());
}

int AliHLTGlobalBarrelTrack::CalculateHelixParams(float bfield)
{
  // calculate radius and center of the helix
  if (TMath::Abs(bfield)<kAlmost0) {
    // no magnetic field -> straight lines
    fHelixRadius=kVeryBig;
    fHelixCenterX=kVeryBig;
    fHelixCenterY=kVeryBig;
  } else {
    fHelixRadius = GetSignedPt()/(-kB2C*bfield);
    Double_t trackPhi = Phi()-GetAlpha();

    //cout << "Helix: phi=" << trackPhi << " x=" << GetX() << " y=" << GetY() << endl;
    fHelixCenterX = GetX() + fHelixRadius *  sin(trackPhi);
    fHelixCenterY = GetY() - fHelixRadius *  cos(trackPhi);
    //cout << "Helix: center" << " x=" << fHelixCenterX << " y=" << fHelixCenterY << endl;
  }
  return 0;
}

int AliHLTGlobalBarrelTrack::CalculateCrossingPoint(float xPlane, float alphaPlane, float& y, float& z)
{
  // calculate crossing point of helix with a plane in yz
  // in the local coordinates of the plane
  int iResult=0;
  if (TMath::Abs(fHelixRadius)<kAlmost0 &&
      (iResult=CalculateHelixParams())<0) {
    return iResult;
  }

  if (TMath::Abs(fHelixRadius)>=kVeryBig) {
    // no magnetic field -> straight lines
  } else {
    // rotate helix center to local coordinates of the plane reference frame
    float cosa=TMath::Cos(alphaPlane-GetAlpha());
    float sina=TMath::Sin(alphaPlane-GetAlpha());
    float cx= fHelixCenterX * cosa + fHelixCenterY * sina;
    float cy=-fHelixCenterX * sina + fHelixCenterY * cosa;

    // crossing point of helix with plane
    Double_t aa = (cx - xPlane)*(cx - xPlane);
    Double_t r2 = fHelixRadius*fHelixRadius;
    if(aa > r2) // no crossing
      return 0;

    Double_t aa2 = sqrt(r2 - aa);
    Double_t y1 = cy + aa2;
    Double_t y2 = cy - aa2;
    y = y1;
    if(TMath::Abs(y2) < TMath::Abs(y1)) y = y2;
  
    // calculate the arc length between (x,y) and (x0,y0) with radius (cx,cy)
    // reference point is (x0,y0) rotated by the diffence of the plane angle and
    // the track reference frame angle alpha
    // 1) angle of (x,y)
    Double_t angle1 = atan2((y - cy),(xPlane - cx));
    if(angle1 < 0) angle1 += TMath::TwoPi();

    // 2) angle of (x0,y0)
    float x0= GetX() * cosa + GetY() * sina;
    float y0=-GetX() * sina + GetY() * cosa;
    Double_t angle2 = atan2((y0 - cy),(x0 - cx));
    if(angle2 < 0) angle2 += TMath::TwoPi();

    // 3) angle between (x,y) and (x0,y0)
    Double_t diffangle = angle1 - angle2;
    diffangle = fmod(diffangle,TMath::TwoPi());

    // 4) arc length
    Double_t arclength = TMath::Abs(diffangle*fHelixRadius);

    // 5) direction depending on whether going outwards or inwards
    int direction=GetX()>xPlane?-1:1;
    z = GetZ() + direction*arclength*GetTgl();

    //cout << "x=" << xPlane << " y=" << y << " cx=" << cx << " cy=" << cy << " a1=" << angle1 << " a2=" << angle2 << " diffa=" << diffangle << " s=" << arclength << " z=" << z << endl;
  }
  return 1;
}

void AliHLTGlobalBarrelTrack::Print(Option_t* option) const
{
  // see header file for class documentation
  cout << "********* Track Id: " << fTrackID << " *******************" << endl;
  AliExternalTrackParam::Print(option);
//   cout << "  Alpha "     << GetAlpha();
//   cout << "  X "         << GetX();
//   cout << "  Y "         << GetY();
//   cout << "  Z "         << GetZ() << endl;
//   cout << "  Snp "       << GetSnp();
//   cout << "  Tgl "       << GetTgl();
//   cout << "  Signed1Pt " << GetSigned1Pt() << endl;
}

void AliHLTGlobalBarrelTrack::Draw(Option_t *option)
{
  /// Inherited from TObject, draw the track
  float scale=250;
  float center[2]={0.5,0.5};

  if (TMath::Abs(fHelixRadius)<kAlmost0 &&
      (CalculateHelixParams())<0) {
    return;
  }

  TString strOption(option);
  if (strOption.IsNull()) strOption="spacepoints trackarc";
  std::auto_ptr<TObjArray> tokens(strOption.Tokenize(" "));
  if (!tokens.get()) return;
  for (int i=0; i<tokens->GetEntriesFast(); i++) {
    if (!tokens->At(i)) continue;
    const char* key="";
    TString arg=tokens->At(i)->GetName();

    key="scale=";
    if (arg.BeginsWith(key)) {
      arg.ReplaceAll(key, "");
      scale=arg.Atof();
      continue;
    }
    key="centerx=";
    if (arg.BeginsWith(key)) {
      arg.ReplaceAll(key, "");
      center[0]=arg.Atof();
      continue;
    }
    key="centery=";
    if (arg.BeginsWith(key)) {
      arg.ReplaceAll(key, "");
      center[1]=arg.Atof();
      continue;
    }
    key="spacepoints";
    if (arg.CompareTo(key)==0) {
      if (fSpacePoints) DrawProjXYSpacePoints(option, fSpacePoints, scale, center);
      continue;
    }
    key="trackarc";
    if (arg.CompareTo(key)==0) {
      DrawProjXYTrack(option, scale, center);
      continue;
    }
  }
}

int AliHLTGlobalBarrelTrack::DrawProjXYSpacePoints(Option_t */*option*/, const AliHLTSpacePointContainer* spacePoints, const float scale, float center[2])
{
  /// draw space points
  int markerColor=3;

  if (!spacePoints) return -EINVAL;

  const UInt_t* pointids=GetPoints();
  for (unsigned i=0; i<GetNumberOfPoints() && pointids; i++) {
    float clusterphi=spacePoints->GetPhi(pointids[i]);
    float cosphi=TMath::Cos(clusterphi);
    float sinphi=TMath::Sin(clusterphi);
    float clusterx=spacePoints->GetX(pointids[i]);
    float clustery=spacePoints->GetY(pointids[i]);
    // rotate
    float pointx= clusterx*sinphi + clustery*cosphi;
    float pointy=-clusterx*cosphi + clustery*sinphi;

    // FIXME: cleanup of marker objects
    TMarker* m=new TMarker(pointx/(2*scale)+center[0], pointy/(2*scale)+center[1], 3);
    m->SetMarkerColor(markerColor);
    m->Draw("same");
  }
  return 0;
}

// FIXME: make this a general geometry definition
// through an abstract class interface
const Double_t gkTPCX[159] = {
  85.195,
  85.945,
  86.695,
  87.445,
  88.195,
  88.945,
  89.695,
  90.445,
  91.195,
  91.945,
  92.695,
  93.445,
  94.195,
  94.945,
  95.695,
  96.445,
  97.195,
  97.945,
  98.695,
  99.445,
  100.195,
  100.945,
  101.695,
  102.445,
  103.195,
  103.945,
  104.695,
  105.445,
  106.195,
  106.945,
  107.695,
  108.445,
  109.195,
  109.945,
  110.695,
  111.445,
  112.195,
  112.945,
  113.695,
  114.445,
  115.195,
  115.945,
  116.695,
  117.445,
  118.195,
  118.945,
  119.695,
  120.445,
  121.195,
  121.945,
  122.695,
  123.445,
  124.195,
  124.945,
  125.695,
  126.445,
  127.195,
  127.945,
  128.695,
  129.445,
  130.195,
  130.945,
  131.695,
  135.180,
  136.180,
  137.180,
  138.180,
  139.180,
  140.180,
  141.180,
  142.180,
  143.180,
  144.180,
  145.180,
  146.180,
  147.180,
  148.180,
  149.180,
  150.180,
  151.180,
  152.180,
  153.180,
  154.180,
  155.180,
  156.180,
  157.180,
  158.180,
  159.180,
  160.180,
  161.180,
  162.180,
  163.180,
  164.180,
  165.180,
  166.180,
  167.180,
  168.180,
  169.180,
  170.180,
  171.180,
  172.180,
  173.180,
  174.180,
  175.180,
  176.180,
  177.180,
  178.180,
  179.180,
  180.180,
  181.180,
  182.180,
  183.180,
  184.180,
  185.180,
  186.180,
  187.180,
  188.180,
  189.180,
  190.180,
  191.180,
  192.180,
  193.180,
  194.180,
  195.180,
  196.180,
  197.180,
  198.180,
  199.430,
  200.930,
  202.430,
  203.930,
  205.430,
  206.930,
  208.430,
  209.930,
  211.430,
  212.930,
  214.430,
  215.930,
  217.430,
  218.930,
  220.430,
  221.930,
  223.430,
  224.930,
  226.430,
  227.930,
  229.430,
  230.930,
  232.430,
  233.930,
  235.430,
  236.930,
  238.430,
  239.930,
  241.430,
  242.930,
  244.430,
  245.930
};

int AliHLTGlobalBarrelTrack::DrawProjXYTrack(Option_t *option, const float scale, float center[2])
{
  /// draw track
  bool bDrawArc=false; // draw TArc
  bool bNoTrackPoints=false; // don't draw track points
  TString strOption(option);
  if (strOption.IsNull()) strOption="spacepoints trackarc";
  std::auto_ptr<TObjArray> tokens(strOption.Tokenize(" "));
  if (!tokens.get()) return 0;
  for (int i=0; i<tokens->GetEntriesFast(); i++) {
    if (!tokens->At(i)) continue;
    const char* key="";
    TString arg=tokens->At(i)->GetName();

    key="drawarc";
    if (arg.BeginsWith(key)) {
      bDrawArc=true;
      continue;
    }
    key="notrackpoints";
    if (arg.BeginsWith(key)) {
      bNoTrackPoints=true;
      continue;
    }
  }

  float cosa=TMath::Cos(GetAlpha());
  float sina=TMath::Sin(GetAlpha());

  // first point
  float firstpoint[2];
  firstpoint[0]= GetX()*sina + GetY()*cosa;
  firstpoint[1]=-GetX()*cosa + GetY()*sina;
  {
    //cout << " first point alpha=" << GetAlpha() << " x: " << firstpoint[0] << " y: " << firstpoint[1] << endl;
    // FIXME: cleanup of marker objects
    TMarker* m=new TMarker(firstpoint[0]/(2*scale)+center[0], firstpoint[1]/(2*scale)+center[1], 29);
    m->SetMarkerSize(2);
    m->SetMarkerColor(2);
    m->Draw("same");
  }

  // draw points in step width and remember the last point
  float lastpoint[2]={0.0, 0.0};
  int firstpadrow=0;
  for (; firstpadrow<159 && gkTPCX[firstpadrow]<GetX(); firstpadrow++);
  DrawProjXYTrackPoints(option, scale, center, firstpadrow, -1, firstpoint);
  DrawProjXYTrackPoints(option, scale, center, firstpadrow, 1, lastpoint);

  if (bDrawArc) {
    if (TMath::Abs(fHelixRadius)>=kVeryBig) {
      // no magnetic field -> straight lines
    } else {
      // rotate helix center to local coordinates of the plane reference frame
      float cx= fHelixCenterX * sina + fHelixCenterY * cosa;
      float cy=-fHelixCenterX * cosa + fHelixCenterY * sina;

      float diffx=cx-firstpoint[0];
      float diffy=cy-firstpoint[1];
      float phimin=0.0;
      float phimax=2*TMath::Pi();
      if (TMath::Abs(diffx)<kAlmost0) {
	phimin=TMath::Pi()/2;
      } else {
	phimin=TMath::ATan(diffy/diffx);
      }
      if (diffx>0) phimin+=TMath::Pi();

      diffx=cx-lastpoint[0];
      diffy=cy-lastpoint[1];
      if (TMath::Abs(diffx)<kAlmost0) {
	phimax=TMath::Pi()/2;
      } else {
	phimax=TMath::ATan(diffy/diffx);
      }
      //cout << "diffx=" << diffx << " diffy=" << diffy << " phimin=" << phimin << " phimax=" << phimax << endl;
      if (diffx>0) phimax+=TMath::Pi();
      if (phimax<0 && phimin>=0 && 
	  phimax+TMath::TwoPi()-phimin<TMath::Pi()) phimax+=TMath::TwoPi();
      //if (phimax<0 && TMath::Abs(phimax-phimin)<TMath::Pi()) phimax+=TMath::TwoPi();
      if (phimax-phimin>TMath::Pi()) phimax-=TMath::TwoPi();

      if (0/*phimin>phimax*/) {
	float tmp=phimin;
	phimin=phimax;
	phimax=tmp;
      }
      phimin*=360.0/(2*TMath::Pi());
      phimax*=360.0/(2*TMath::Pi());
      //cout << " cx=" << cx << " cy=" << cy << " r=" << fHelixRadius << " phimin=" << phimin << " phimax=" << phimax << endl;
      // FIXME: cleanup of graphics objects
      TArc* tarc=new TArc(cx/(2*scale)+center[0], cy/(2*scale)+center[1], TMath::Abs(fHelixRadius)/(2*scale), phimin, phimax);
      tarc->SetNoEdges();
      tarc->SetFillStyle(0);
      tarc->Draw("same");
    }
  }
  return 0;
}

int AliHLTGlobalBarrelTrack::DrawProjXYTrackPoints(Option_t */*option*/, const float scale, const float center[2], int firstpadrow, int step, float point[2])
{
  // draw points in step width and return the last point
  float offsetAlpha=0.0;
  float cosa=TMath::Cos(GetAlpha());
  float sina=TMath::Sin(GetAlpha());
  int markerColor=1;
  for (int padrow=firstpadrow; padrow>=0 && padrow<159; padrow+=step) {
    float x=gkTPCX[padrow];
    float y=0.0;
    float z=0.0;

    int maxshift=9;
    int shift=0;
    int result=0;
    do {
      if ((result=CalculateCrossingPoint(x, GetAlpha()-offsetAlpha, y, z))<1) break;
      float pointAlpha=TMath::ATan(y/x);
      if (TMath::Abs(pointAlpha)>TMath::Pi()/18) {
	offsetAlpha+=(pointAlpha>0?-1:1)*TMath::Pi()/9;
      	result=0;
	markerColor++;
	cosa=TMath::Cos(GetAlpha()-offsetAlpha);
	sina=TMath::Sin(GetAlpha()-offsetAlpha);
      }
    } while (result==0 && shift++<maxshift);
    if (result<1) continue;
    point[0]= x*sina + y*cosa;
    point[1]=-x*cosa + y*sina;

    //cout << x << " : x=" << x << " y=" << y << endl;
    // FIXME: cleanup of TMarker objects?
    TMarker* m=new TMarker(point[0]/(2*scale)+center[0], point[1]/(2*scale)+center[1], z>=0?2:5);
    m->SetMarkerColor(markerColor);
    m->Draw("same");
  }
  return 0;
}

void AliHLTGlobalBarrelTrack::SetTrackGeometry(AliHLTTrackGeometry* points)
{
  /// set the instance to the track points container
  fTrackPoints=points; 
  if (!fTrackPoints) return; 
  fTrackPoints->SetTrackId(GetID());
}

int AliHLTGlobalBarrelTrack::AssociateSpacePoints(AliHLTTrackGeometry* trackpoints, AliHLTSpacePointContainer& spacepoints) const
{
  /// associate the track space points to the calculated track points
  AliHLTTrackGeometry* instance=trackpoints;
  if (!instance) instance=fTrackPoints;
  if (!instance) return 0;

  UInt_t nofIds=GetNumberOfPoints();
  const UInt_t* ids=GetPoints();
  int result=instance->AssociateSpacePoints(ids, nofIds, spacepoints);
  return result;
}

int AliHLTGlobalBarrelTrack::ReadTracks(const char* filename, TClonesArray& tgt, AliHLTComponentDataType /*dt*/, unsigned /*specification*/)
{
  // open block from file and add to collection
  if (!filename) return -EINVAL;
  
  TString input=filename;
  input+="?filetype=raw";
  std::auto_ptr<TFile> pFile(new TFile(input));
  if (!pFile.get()) return -ENOMEM;
  if (pFile->IsZombie()) return -ENOENT;

  int iResult=0;
  pFile->Seek(0);
  std::auto_ptr<TArrayC> buffer(new TArrayC);
  if (!buffer.get()) return -ENOMEM;

  buffer->Set(pFile->GetSize());
  if (pFile->ReadBuffer(buffer->GetArray(), buffer->GetSize())==0) {
    const AliHLTTracksData* pTracks=reinterpret_cast<const AliHLTTracksData*>(buffer->GetArray());
    vector<AliHLTGlobalBarrelTrack> tracks;
    iResult=ConvertTrackDataArray(pTracks, buffer->GetSize(), tracks);
    if (iResult>=0) {
      int offset=tgt.GetEntriesFast();
      tgt.ExpandCreate(offset+tracks.size());
      for (unsigned i=0; i<tracks.size(); i++) {
	new (tgt[offset+i]) AliHLTGlobalBarrelTrack(tracks[i]);
      }
      iResult=tracks.size();
    } else {
      //HLTError("failed to convert tracks from file %s size %d byte(s) ", filename, pFile->GetSize());
    }
  } else {
    //HLTError("failed reading %d byte(s) from file %s", pFile->GetSize(), filename);
    iResult=-ENODATA;
  }

  return iResult;
}

int AliHLTGlobalBarrelTrack::ReadTrackList(const char* listfile, TClonesArray& tgt, AliHLTComponentDataType dt, unsigned specification)
{
  // open blank separated list of files and read tracks
  ifstream list(listfile);
  if (!list.good()) return -ENOENT;

  int count=0;
  TString file;
  while (file.ReadLine(list)) {
    //HLTInfo("adding tracks from file %s", file.Data());
    int iResult=ReadTracks(file.Data(), tgt, dt, specification);
    if (iResult<0) {
      //HLTInfo("failed to add data from file %s: error %d", file.Data(), iResult);
      return iResult;
    }
    count+=iResult;
  }

  return count;
}
