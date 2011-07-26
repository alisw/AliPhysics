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

/// @file   AliHLTTrackGeometry.cxx
/// @author Matthias Richter
/// @date   2011-05-20
/// @brief  Desciption of a track by a sequence of track points
///

#include "AliHLTTrackGeometry.h"
#include "AliHLTSpacePointContainer.h"
#include "TObjArray.h"
#include "TMarker.h"
#include "TMath.h"
#include <memory>
#include <iostream>
#include <algorithm>

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTTrackGeometry)

AliHLTTrackGeometry::AliHLTTrackGeometry()
  : TObject(), AliHLTLogging()
  , fTrackPoints()
  , fTrackId(-1)
{
  /// standard constructor
}

AliHLTTrackGeometry::AliHLTTrackGeometry(const AliHLTTrackGeometry& src)
  : TObject(src), AliHLTLogging()
  , fTrackPoints(src.fTrackPoints)
  , fTrackId(src.fTrackId)
{
  /// copy constructor
}

AliHLTTrackGeometry& AliHLTTrackGeometry::operator=(const AliHLTTrackGeometry& src)
{
  /// assignment operator
  if (this!=&src) {
    fTrackPoints.assign(src.fTrackPoints.begin(), src.fTrackPoints.end());
    fTrackId=src.fTrackId;
  }
  return *this;
}

AliHLTTrackGeometry::~AliHLTTrackGeometry()
{
  /// destructor
}

int AliHLTTrackGeometry::AddTrackPoint(const AliHLTTrackPoint& point)
{
  /// add a track point to the list
  vector<AliHLTTrackPoint>::const_iterator element = find(fTrackPoints.begin(), fTrackPoints.end(), point);
  if (element==fTrackPoints.end()) {
    fTrackPoints.push_back(point);
  } else {
    HLTError("track point of id %08x already existing", point.GetId());
    return -EEXIST;
  }
  return 0;
}

void AliHLTTrackGeometry::Clear(Option_t * /*option*/)
{
  // internal cleanup
}

void AliHLTTrackGeometry::Print(Option_t *option) const
{
  // print info
  Print(cout, option);
}

void AliHLTTrackGeometry::Print(ostream& out, Option_t */*option*/) const
{
  // print to stream
  out << "AliHLTTrackGeometry::Print" << endl;
}

void AliHLTTrackGeometry::Draw(Option_t *option)
{
  /// Inherited from TObject, draw the track
  float scale=250;
  float center[2]={0.5,0.5};
  int markerColor=1;
  int markerSize=1;
  int verbosity=0;

  TString strOption(option);
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

    key="markercolor=";
    if (arg.BeginsWith(key)) {
      arg.ReplaceAll(key, "");
      markerColor=arg.Atoi();
      continue;
    }

    key="markersize=";
    if (arg.BeginsWith(key)) {
      arg.ReplaceAll(key, "");
      markerSize=arg.Atoi();
      continue;
    }

    key="verbosity=";
    if (arg.BeginsWith(key)) {
      arg.ReplaceAll(key, "");
      verbosity=arg.Atoi();
      continue;
    }
  }

  bool bFirstPoint=true;
  float firstalpha=0.0;
  for (vector<AliHLTTrackPoint>::const_iterator point=fTrackPoints.begin();
       point!=fTrackPoints.end(); 
       point++) {
    float alpha=GetPlaneAlpha(point->GetId());
    float r=GetPlaneR(point->GetId());
    float cosa=TMath::Cos(alpha);
    float sina=TMath::Sin(alpha);
    float x = r*sina + point->GetU()*cosa;
    float y =-r*cosa + point->GetU()*sina;
    if (verbosity>0) {
      HLTInfo("ID 0x%08x: x=% .4f y=% .4f alpha=% .4f", point->GetId(), r, point->GetU(), alpha);
    }
    int color=markerColor;
    if (bFirstPoint) {
      bFirstPoint=false;
      TMarker* m=new TMarker(x/(2*scale)+center[0], y/(2*scale)+center[1], 29);
      m->SetMarkerSize(2);
      m->SetMarkerColor(2);
      m->Draw("same");
      firstalpha=alpha;
    } else {
      color+=int(9*TMath::Abs(alpha-firstalpha)/TMath::Pi());
    }
    TMarker* m=new TMarker(x/(2*scale)+center[0], y/(2*scale)+center[1], point->GetV()>0?2:5);
    m->SetMarkerColor(color);
    m->SetMarkerSize(markerSize);
    m->Draw("same");
  }
}

int AliHLTTrackGeometry::SetAssociatedSpacePoint(UInt_t planeId, UInt_t spacepointId, int status)
{
  /// set the spacepoint associated with a track point
  vector<AliHLTTrackPoint>::iterator element = find(fTrackPoints.begin(), fTrackPoints.end(), planeId);
  if (element==fTrackPoints.end()) return -ENOENT;
  element->SetAssociatedSpacePoint(spacepointId, status);
  return 0;
}

int AliHLTTrackGeometry::GetAssociatedSpacePoint(UInt_t planeId, UInt_t& spacepointId) const
{
  /// get the spacepoint associated with a track point
  /// return status flag if found, -ENOENT if no associated spacepoint found
  vector<AliHLTTrackPoint>::const_iterator element = find(fTrackPoints.begin(), fTrackPoints.end(), planeId);
  if (element==fTrackPoints.end()) return -ENOENT;
  if (!element->HaveAssociatedSpacePoint()) return -ENODATA;
  return element->GetAssociatedSpacePoint(spacepointId);
}

const AliHLTTrackGeometry::AliHLTTrackPoint* AliHLTTrackGeometry::GetTrackPoint(AliHLTUInt32_t id) const
{
  /// get const pointer to track point
  vector<AliHLTTrackPoint>::const_iterator element = find(fTrackPoints.begin(), fTrackPoints.end(), id);
  if (element==fTrackPoints.end()) return NULL;
  return &(*element);
}

AliHLTTrackGeometry::AliHLTTrackPoint* AliHLTTrackGeometry::GetTrackPoint(AliHLTUInt32_t id)
{
  /// get const pointer to track point
  vector<AliHLTTrackPoint>::iterator element = find(fTrackPoints.begin(), fTrackPoints.end(), id);
  if (element==fTrackPoints.end()) return NULL;
  return &(*element);
}

AliHLTSpacePointContainer* AliHLTTrackGeometry::ConvertToSpacePoints() const
{
  /// create a collection of all points
  HLTError("implementation of child method missing");
  return NULL;
}

int AliHLTTrackGeometry::AssociateSpacePoints(AliHLTSpacePointContainer& points)
{
  /// associate the track space points to the calculated track points
  vector<AliHLTUInt32_t> ids;
  points.GetClusterIDs(ids);
  if (ids.size()>0) return 0;
  int result=AssociateSpacePoints(&ids[0], ids.size(), points);
  if (result>0) {
    HLTInfo("associated %d space point(s) to track points", result);
  }
  return result;
}

int AliHLTTrackGeometry::AssociateSpacePoints(const AliHLTUInt32_t* trackpoints, AliHLTUInt32_t nofPoints, AliHLTSpacePointContainer& points)
{
  /// associate the track space points to the calculated track points
  if (nofPoints==0) return 0;
  if (trackpoints==NULL) return -EINVAL;
  int count=0;
  for (int i=nofPoints-1; i>=0; i--) {
    if (!points.Check(trackpoints[i])) {
      HLTWarning("can not find point id %08x", trackpoints[i]);
      continue;
    }
    float xyz[3]={points.GetX(trackpoints[i]), points.GetY(trackpoints[i]), points.GetZ(trackpoints[i])};
    AliHLTUInt32_t planeId=0;
    int result=FindMatchingTrackPoint(trackpoints[i], xyz, planeId);
    if (result<0) {
      HLTWarning("no associated track point found for space point id %08x x=%f y=%f z=%f", trackpoints[i], xyz[0], xyz[1], xyz[2]);
      continue;
    } else if (result==0) {
      HLTWarning("associated track point for space pointid %08x x=%f y=%f z=%f occupied", trackpoints[i], xyz[0], xyz[1], xyz[2]);
      continue;
    }
    SetAssociatedSpacePoint(planeId, trackpoints[i], 1);
    if (points.GetTrackID(trackpoints[i])<0 && GetTrackId()>=0) {
      points.SetTrackID(GetTrackId(), trackpoints[i]);
      HLTDebug("associating unused cluster %08x with track %d", trackpoints[i], GetTrackId());
    }
    count++;
  }
  return count;
}

int AliHLTTrackGeometry::AssociateUnusedSpacePoints(AliHLTSpacePointContainer& points)
{
  /// associate the track space points to the calculated track points
  int count=0;
  vector<AliHLTUInt32_t> ids;
  points.GetClusterIDs(ids);
  for (vector<AliHLTUInt32_t>::iterator id=ids.begin();
       id!=ids.end(); id++) {
    float xyz[3]={points.GetX(*id), points.GetY(*id), points.GetZ(*id)};
    AliHLTUInt32_t planeId=0;
    int result=FindMatchingTrackPoint(*id, xyz, planeId);
    if (result<0) {
      //HLTWarning("no associated track point found for space point id %08x x=%f y=%f z=%f", *id, xyz[0], xyz[1], xyz[2]);
      continue;
    } else if (result==0) {
      //HLTWarning("associated track point for space pointid %08x x=%f y=%f z=%f occupied", *id, xyz[0], xyz[1], xyz[2]);
      continue;
    }
    SetAssociatedSpacePoint(planeId, *id, 1);
    if (points.GetTrackID(*id)<0 && GetTrackId()>=0) {
      points.SetTrackID(GetTrackId(), *id);
      HLTDebug("associating unused cluster %08x with track %d", *id, GetTrackId());
    }
    count++;
  }
  return count;
}

ostream& operator<<(ostream &out, const AliHLTTrackGeometry& p)
{
  p.Print(out);
  return out;
}
