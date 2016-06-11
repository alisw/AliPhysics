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

/// \class AliTPCclusterMI
/// \brief Implementation of the TPC cluser
///
/// AliTPC parallel tracker -
/// Description of this class together with its intended usage
/// will follow shortly
///
/// \author Marian Ivanov   Marian.Ivanov@cern.ch

/* $Id$ */

#include <TGeoManager.h>
#include <TGeoMatrix.h>
#include "AliTPCclusterMI.h"
//#include "AliTPCclusterInfo.h"
#include "AliTrackPointArray.h"
#include "AliGeomManager.h"
#include "AliLog.h"

/// \cond CLASSIMP
ClassImp(AliTPCclusterMI)
/// \endcond


AliTPCclusterMI::AliTPCclusterMI():
  AliCluster(),
//  fInfo(0),
  fTimeBin(0),  //time bin coordinate
  fPad(0),  //pad coordinate
  fQ(0),       //Q of cluster (in ADC counts)
  fMax(0),      //maximal amplitude in cluster
  fType(0),     //type of the cluster 0 means golden
  fUsed(0),     //counter of usage
  fDisp(0),     /// dispersion of applied correction
  fDetector(0), //detector  number
  fRow(0)      //row number number
{
  //
  // default constructor
  //
}

AliTPCclusterMI::AliTPCclusterMI(const AliTPCclusterMI & cluster):
  AliCluster(cluster),
  //  fInfo(0),
  fTimeBin(cluster.fTimeBin),
  fPad(cluster.fPad),
  fQ(cluster.fQ),
  fMax(cluster.fMax),
  fType(cluster.fType),
  fUsed(cluster.fUsed),
  fDisp(cluster.fDisp),
  fDetector(cluster.fDetector),
  fRow(cluster.fRow)
{
  /// copy constructor

  // AliInfo("Copy constructor\n");

  //  if (cluster.fInfo) fInfo = new AliTPCclusterInfo(*(cluster.fInfo));
}

AliTPCclusterMI & AliTPCclusterMI::operator = (const AliTPCclusterMI & cluster)
{
  /// assignment operator

  // AliInfo("Asignment operator\n");

  if (this == &cluster) return (*this);

  (AliCluster&)(*this) = (AliCluster&)cluster;
  fQ    = cluster.fQ;
  fType = cluster.fType;
  fMax  = cluster.fMax;
  fUsed = cluster.fUsed;
  fDisp = cluster.fDisp;
  fDetector = cluster.fDetector;
  fRow  = cluster.fRow;
  fTimeBin = cluster.fTimeBin;
  fPad     = cluster.fPad;
  //  delete fInfo;
  //  fInfo = 0;
  //  if (cluster.fInfo) fInfo = new AliTPCclusterInfo(*(cluster.fInfo));
  return *this;
}




AliTPCclusterMI::AliTPCclusterMI(Int_t *lab, Float_t *hit) :
  AliCluster(0,hit,0.,0.,lab),
  //  fInfo(0),
  fTimeBin(0),  //time bin coordinate
  fPad(0),  //pad coordinate
  fQ(0),       //Q of cluster (in ADC counts)
  fMax(0),      //maximal amplitude in cluster
  fType(0),     //type of the cluster 0 means golden
  fUsed(0),     //counter of usage
  fDisp(0),     // distortion dispersion
  fDetector(0), //detector  number
  fRow(0)      //row number number
{
  /// constructor

  fQ = (UShort_t)hit[4];
  //  fInfo = 0;
}

AliTPCclusterMI::~AliTPCclusterMI() {
  /// destructor

  //  if (fInfo) delete fInfo;
  //  fInfo = 0;
}



Bool_t AliTPCclusterMI::IsSortable() const
{
  ///

  return kTRUE;

}

Int_t AliTPCclusterMI::Compare(const TObject* obj) const
{
  /// compare according y

  AliTPCclusterMI * o2 = (AliTPCclusterMI*)obj;
  return (o2->GetY()>GetY())? -1:1;
}


void AliTPCclusterMI::SetDetector(Int_t detector){
  /// set volume ID

  fDetector = (UChar_t)(detector%72);
  AliGeomManager::ELayerID id = (fDetector<36) ?
    AliGeomManager::kTPC1 :AliGeomManager::kTPC2 ;
  Int_t modId = (fDetector<36)?fDetector: fDetector-36;
  SetVolumeId(AliGeomManager::LayerToVolUID(id,modId));
}

/*
void AliTPCclusterMI::SetInfo(AliTPCclusterInfo * info) {
  ///

  if (fInfo) delete fInfo;
  fInfo = info;
}
*/

AliTPCclusterMI* AliTPCclusterMI::MakeCluster(AliTrackPoint* /*point*/) {
  /// make AliTPCclusterMI out of AliTrackPoint
  /// (not yet implemented)

  return NULL;
}


AliTrackPoint* AliTPCclusterMI::MakePoint() {
  /// make AliTrackPoint out of AliTPCclusterMI

  AliTrackPoint* point = new AliTrackPoint();
  Float_t xyz[3]={0.};
  Float_t cov[6]={0.};
  GetGlobalXYZ(xyz);
  GetGlobalCov(cov);
  // voluem ID to add later ....
  point->SetXYZ(xyz);
  point->SetCov(cov);

  return point;
}

//______________________________________________________________________________
void AliTPCclusterMI::SetGlobalTrackPoint( const AliCluster &cl, AliTrackPoint &point )
{
  /// Set global AliTrackPoint

  Float_t xyz[3]={0.};
  Float_t cov[6]={0.};
  cl.GetGlobalXYZ(xyz);
  cl.GetGlobalCov(cov);
  // voluem ID to add later ....
  point.SetXYZ(xyz);
  point.SetCov(cov);
}

//_____________________________________________________
void AliTPCclusterMI::SetDistortionDispersion(float d)
{
  // set distortion dispersion
  if (d<0) d = 0;
  UInt_t di = d*kScaleDisp;
  if (di>kMaxDisp) di = kMaxDisp;
  fDisp = di;
}

//_____________________________________________________
Float_t AliTPCclusterMI::GetDistortionDispersion() const
{
  // get distortion dispersion
  return float(fDisp)/kScaleDisp;
}

//_____________________________________________________
void AliTPCclusterMI::SetDistortions(float dx, float dy, float dz)
{
  // store distortions rounded: to 0.2 for X (9 bits) and to 0.1mm (11 bits) for y and z
  int pack = 0;
  int dxi = dx*kScaleDX;
  if (dxi<0) {
    if (dxi<-kMaxDX) dxi = -kMaxDX;
    pack |= ((-dxi)&kMaxDX)|(0x1<<(kNBitsDX-1)); // highest bit flags negative value
  }
  else {
    if (dxi>kMaxDX) dxi = kMaxDX;
    pack |= dxi&kMaxDX;
  }
  int dyi = dy*kScaleDY;
  if (dyi<0) {
    if (dyi<-kMaxDY) dyi = -kMaxDY;
    pack |= (((-dyi)&kMaxDY)|(0x1<<(kNBitsDY-1)))<<kNBitsDX; // highest bit flags negative value
  }
  else {
    if (dyi>kMaxDY) dyi = kMaxDY;
    pack |= (dyi&kMaxDY)<<kNBitsDX;
  }
  int dzi = dz*kScaleDZ;
  if (dzi<0) {
    if (dzi<-kMaxDZ) dzi = -kMaxDZ;
    pack |= (((-dzi)&kMaxDZ)|(0x1<<(kNBitsDZ-1)))<<(kNBitsDX+kNBitsDY); // highest bit flags negative value
  }
  else {
    if (dzi>kMaxDZ) dzi = kMaxDZ;
    pack |= (dzi&kMaxDZ)<<(kNBitsDX+kNBitsDY);
  }
  SetSigmaYZ(*(float*)&pack); // interpret as float
  //
}

//_____________________________________________________
void AliTPCclusterMI::GetDistortions(float &dx, float &dy, float &dz) const
{
  // Extract rounded distortions
  float v = GetSigmaYZ();
  int pack = *(int*)&v;
  int dxi = pack&kMaskDX;
  dx = dxi>kMaxDX ? -(dxi&kMaxDX) : dxi;
  dx *= 1.0f/kScaleDX; //1./50.0f;
  //
  int dyi = (pack>>kNBitsDX)&kMaskDY;
  dy = dyi>kMaxDY ? -(dyi&kMaxDY) : dyi;
  dy *= 1.0f/kScaleDY; //1./100.0f;
  //
  int dzi = (pack>>(kNBitsDX+kNBitsDY))&kMaskDZ;
  dz = dzi>kMaxDZ ? -(dzi&kMaxDZ) : dzi;
  dz *= 1.0f/kScaleDZ; //1./100.0f;
  //
}

//_____________________________________________________
Float_t AliTPCclusterMI::GetDistortionX() const
{
  // Extract rounded distortions
  float v = GetSigmaYZ();
  int pack = *(int*)&v;
  int dxi = pack&kMaskDX;
  float dx = dxi>kMaxDX ? -(dxi&kMaxDX) : dxi;
  dx *= 1.0f/kScaleDX; //1./50.0f;
  return dx;
}

//_____________________________________________________
Float_t AliTPCclusterMI::GetDistortionY() const
{
  // Extract rounded distortions
  float v = GetSigmaYZ();
  int pack = *(int*)&v;
  //
  int dyi = (pack>>kNBitsDX)&kMaskDY;
  float dy = dyi>kMaxDY ? -(dyi&kMaxDY) : dyi;
  dy *= 1.0f/kScaleDY; //1./100.0f;
  return dy;
  //
}

//_____________________________________________________
Float_t AliTPCclusterMI::GetDistortionZ() const
{
  // Extract rounded distortions
  float v = GetSigmaYZ();
  int pack = *(int*)&v;
  //
  int dzi = (pack>>(kNBitsDX+kNBitsDY))&kMaskDZ;
  float dz = dzi>kMaxDZ ? -(dzi&kMaxDZ) : dzi;
  dz *= 1.0f/kScaleDZ; //1./100.0f;
  return dz;
  //
}

//______________________________________________________________________________
Bool_t AliTPCclusterMI::GetGlobalCov(Float_t cov[6]) const
{
  // clone of the AliCluster::GetGlobalCov avoiding using kSigmaYZ (used to store
  // distortions) as a real error
    for (Int_t i = 0; i < 6; i++) cov[i] = 0;

  if (!gGeoManager || !gGeoManager->IsClosed()) {
    AliError("Can't get the global coordinates! gGeoManager doesn't exist or it is still opened!");
    return kFALSE;
  }

  const TGeoHMatrix *mt = GetTracking2LocalMatrix();
  if (!mt) return kFALSE;

  TGeoHMatrix *ml = GetMatrix();
  if (!ml) return kFALSE;

  TGeoHMatrix m;
  Double_t tcov[9] = { 0, 0, 0, 0, GetSigmaY2(), 0, 0, 0, GetSigmaZ2() };
  m.SetRotation(tcov);
  m.Multiply(&mt->Inverse());
  m.Multiply(&ml->Inverse());
  m.MultiplyLeft(mt);
  m.MultiplyLeft(ml);
  Double_t *ncov = m.GetRotationMatrix();
  cov[0] = ncov[0]; cov[1] = ncov[1]; cov[2] = ncov[2];
  cov[3] = ncov[4]; cov[4] = ncov[5];
  cov[5] = ncov[8];

  return kTRUE;
}
