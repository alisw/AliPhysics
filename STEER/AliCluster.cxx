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

//-------------------------------------------------------------------------
//                         Class AliCluster
// This is the future base for managing the clusters in barrel detectors.
// It is fully interfaced with the ROOT geometrical modeller TGeo.
// Each cluster contains XYZ coordinates in the local tracking c.s. and
// the unique ID of the sensitive detector element which continas the
// cluster. The coordinates in global c.s. are computed using the interface
// to TGeo and will be not overwritten by the derived sub-detector cluster
// classes.
//
// cvetan.cheshkov@cern.ch & jouri.belikov@cern.ch    5/3/2007
//-------------------------------------------------------------------------

#include <TGeoManager.h>
#include <TGeoMatrix.h>
#include <TGeoPhysicalNode.h>

#include "AliCluster.h"
#include "AliLog.h"
#include "AliAlignObj.h"

ClassImp(AliCluster)

//______________________________________________________________________________
AliCluster::AliCluster():
  TObject(),
  fX(0),
  fY(0),
  fZ(0),
  fSigmaY2(0),
  fSigmaZ2(0),
  fSigmaYZ(0),
  fVolumeId(0),
  fIsMisaligned(kFALSE)
{
  // Default constructor
  fTracks[0]=fTracks[1]=fTracks[2]=-3141593; 
}

//______________________________________________________________________________
AliCluster::AliCluster(UShort_t volId,
			       const Float_t *hit,
			       Float_t x,
			       Float_t sigyz,
			       const Int_t *lab):
  TObject(),
  fX(x),
  fY(hit[0]),
  fZ(hit[1]),
  fSigmaY2(hit[2]),
  fSigmaZ2(hit[3]),
  fSigmaYZ(sigyz),
  fVolumeId(volId),
  fIsMisaligned(kFALSE)
{
  // Constructor
  if (lab) {
    fTracks[0] = lab[0];
    fTracks[1] = lab[1];
    fTracks[2] = lab[2];
  }
  else
    fTracks[0]=fTracks[1]=fTracks[2]=-3141593; 
}

//______________________________________________________________________________
AliCluster::AliCluster(UShort_t volId,
			       Float_t x, Float_t y, Float_t z,
			       Float_t sy2, Float_t sz2, Float_t syz,
			       const Int_t *lab):
  TObject(),
  fX(x),
  fY(y),
  fZ(z),
  fSigmaY2(sy2),
  fSigmaZ2(sz2),
  fSigmaYZ(syz),
  fVolumeId(volId),
  fIsMisaligned(kFALSE)
{
  // Constructor
  if (lab) {
    fTracks[0] = lab[0];
    fTracks[1] = lab[1];
    fTracks[2] = lab[2];
  }
  else
    fTracks[0]=fTracks[1]=fTracks[2]=-3141593; 
}

//______________________________________________________________________________
AliCluster::AliCluster(const AliCluster& cluster):
  TObject(cluster),
  fX(cluster.fX),
  fY(cluster.fY),
  fZ(cluster.fZ),
  fSigmaY2(cluster.fSigmaY2),
  fSigmaZ2(cluster.fSigmaZ2),
  fSigmaYZ(cluster.fSigmaYZ),
  fVolumeId(cluster.fVolumeId),
  fIsMisaligned(cluster.fIsMisaligned)
{
  // Copy constructor
  fTracks[0] = cluster.fTracks[0];
  fTracks[1] = cluster.fTracks[1];
  fTracks[2] = cluster.fTracks[2];
}

//______________________________________________________________________________
AliCluster & AliCluster::operator=(const AliCluster& cluster)
{
  // Assignment operator

  if(&cluster == this) return *this;

  fX = cluster.fX;
  fY = cluster.fY;
  fZ = cluster.fZ;
  fSigmaY2 = cluster.fSigmaY2;
  fSigmaZ2 = cluster.fSigmaZ2;
  fSigmaYZ = cluster.fSigmaYZ;
  fVolumeId = cluster.fVolumeId;
  fIsMisaligned = cluster.fIsMisaligned;

  fTracks[0] = cluster.fTracks[0];
  fTracks[1] = cluster.fTracks[1];
  fTracks[2] = cluster.fTracks[2];

  return *this;
}

//______________________________________________________________________________
Bool_t AliCluster::GetGlobalXYZ(Float_t xyz[3]) const
{
  // Get the global coordinates of the cluster
  // All the needed information is taken only
  // from TGeo.

  xyz[0] = xyz[1] = xyz[2] = 0;

  if (!gGeoManager || !gGeoManager->IsClosed()) {
    AliError("Can't get the global coordinates! gGeoManager doesn't exist or it is still opened!");
    return kFALSE;
  }

  const TGeoHMatrix *mt = GetTracking2LocalMatrix();
  if (!mt) return kFALSE;
  Double_t txyz[3] = {fX, fY, fZ};
  Double_t lxyz[3] = {0, 0, 0};
  mt->LocalToMaster(txyz,lxyz);

  TGeoHMatrix *ml = GetMatrix();
  if (!ml) return kFALSE;
  Double_t gxyz[3] = {0, 0, 0};
  ml->LocalToMaster(lxyz,gxyz);
  xyz[0] = gxyz[0]; xyz[1] = gxyz[1]; xyz[2] = gxyz[2];
  return kTRUE;
}

//______________________________________________________________________________
Bool_t AliCluster::GetGlobalCov(Float_t cov[6]) const
{
  // Get the global covariance matrix of the cluster coordinates
  // All the needed information is taken only
  // from TGeo.
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
  Double_t tcov[9] = { 0, 0, 0, 0, fSigmaY2, fSigmaYZ, 0, fSigmaYZ, fSigmaZ2 };
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

//______________________________________________________________________________
Bool_t AliCluster::GetXRefPlane(Float_t &xref) const
{
  // Get the distance between the origin and the ref.plane.
  // All the needed information is taken only
  // from TGeo.
  xref = 0;

  const TGeoHMatrix *mt = GetTracking2LocalMatrix();
  if (!mt) return kFALSE;

  TGeoHMatrix *ml = GetMatrix();
  if (!ml) return kFALSE;

  TGeoHMatrix m = *mt;
  m.MultiplyLeft(ml);

  xref = (m.Inverse()).GetTranslation()[0];
  return kTRUE;
}

//______________________________________________________________________________
Bool_t AliCluster::Misalign()
{
  // ...
  // All the needed information is taken only
  // from TGeo.
  if (!gGeoManager || !gGeoManager->IsClosed()) {
    AliError("Can't get the PN entry! gGeoManager doesn't exist or it is still opened!");
    return kFALSE;
  }

  if (fIsMisaligned) {
    AliError("The cluster was already misaligned!");
    return kFALSE;
  }

  const TGeoHMatrix *mt = GetTracking2LocalMatrix();
  if (!mt) return kFALSE;

  TGeoHMatrix *ml = GetMatrix();
  if (!ml) return kFALSE;

  TGeoHMatrix *mlorig = GetMatrix(kTRUE);
  if (!mlorig) return kFALSE;

  TGeoHMatrix delta = *mt;
  delta.MultiplyLeft(ml);
  delta.MultiplyLeft(&(mlorig->Inverse()));
  delta.MultiplyLeft(&(mt->Inverse()));

  Double_t xyzorig[3] = {fX, fY, fZ};
  Double_t xyz[3] = {0, 0, 0};
  delta.LocalToMaster(xyzorig,xyz);
  fX = xyz[0]; fY = xyz[1]; fZ = xyz[2];
  fIsMisaligned = kTRUE;
  return kTRUE;
}

//______________________________________________________________________________
TGeoHMatrix* AliCluster::GetMatrix(Bool_t original) const
{
  // Get the matrix which transforms from the
  // local TGeo alignable volume c.s. to the global one.
  // In case the cluster was already misaligned, get the
  // ideal matrix from TGeo. The option 'original'
  // can be used to force the calculation of the ideal
  // matrix.
  if (!fIsMisaligned && (original == kFALSE)) {
    TGeoPNEntry *pne = GetPNEntry();
    if (!pne) return NULL;

    TGeoPhysicalNode *pnode = pne->GetPhysicalNode();
    if (pnode) return pnode->GetMatrix();

    const char* path = pne->GetTitle();
    if (!gGeoManager->cd(path)) {
      AliError(Form("Volume path %s not valid!",path));
      return NULL;
    }
    return gGeoManager->GetCurrentMatrix();
  }
  else {
    const char* symname = AliAlignObj::SymName(fVolumeId);
    if (!symname) return NULL;

    static TGeoHMatrix m;
    if (AliAlignObj::GetOrigGlobalMatrix(symname,m))
      return &m;
    else
      return NULL;
  }
}

//______________________________________________________________________________
const TGeoHMatrix* AliCluster::GetTracking2LocalMatrix() const
{
  // Get the matrix which is stored with the PN entries in TGeo.
  // The matrix makes the transformation from the tracking c.s. to
  // global one.
  TGeoPNEntry *pne = GetPNEntry();
  if (!pne) return NULL;

  const TGeoHMatrix *m = pne->GetMatrix();
  if (!m)
    AliError(Form("TGeoPNEntry (%s) contains no matrix !",pne->GetName()));

  return m;
}

//______________________________________________________________________________
TGeoPNEntry* AliCluster::GetPNEntry() const
{
  // Get a pointer to the physical node entry
  // corresponding to the alignable volume to
  // which the cluster belongs

  if (!gGeoManager || !gGeoManager->IsClosed()) {
    AliError("Can't get the PN entry! gGeoManager doesn't exist or it is still opened!");
    return NULL;
  }

  const char* symname = AliAlignObj::SymName(fVolumeId);
  if (!symname) return NULL;

  TGeoPNEntry* pne = gGeoManager->GetAlignableEntry(symname);
  if (!pne)
    AliError(Form("The symbolic volume name %s does not correspond to a physical entry!",
		  symname));
  return pne;
}
