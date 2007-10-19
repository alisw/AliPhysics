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

//-----------------------------------------------------------------
//  Implementation of the alignment object class, holding the alignment
//  constants for a single volume, through the abstract class AliAlignObj.
//  From it two derived concrete representation of alignment object class
//  (AliAlignObjParams, AliAlignObjMatrix) are derived in separate files.
//-----------------------------------------------------------------

#include <TGeoManager.h>
#include <TGeoPhysicalNode.h>
#include <TMath.h>
#include <TMatrixDSym.h>

#include "AliAlignObj.h"
#include "AliTrackPointArray.h"
#include "AliLog.h"
 
ClassImp(AliAlignObj)

//_____________________________________________________________________________
AliAlignObj::AliAlignObj():
  fVolPath(),
  fVolUID(0)
{
  // default constructor
  for(Int_t i=0; i<6; i++) fDiag[i]=-999.;
  for(Int_t i=0; i<21; i++) fODia[i]=-999.;
}

//_____________________________________________________________________________
AliAlignObj::AliAlignObj(const char* symname, UShort_t voluid) :
  TObject(),
  fVolPath(symname),
  fVolUID(voluid)
{
  // standard constructor
  //
  for(Int_t i=0; i<6; i++) fDiag[i]=-999.;
  for(Int_t i=0; i<21; i++) fODia[i]=-999.;
}

//_____________________________________________________________________________
AliAlignObj::AliAlignObj(const char* symname, UShort_t voluid, Double_t* cmat) :
  TObject(),
  fVolPath(symname),
  fVolUID(voluid)
{
  // standard constructor
  //
  SetCorrMatrix(cmat);
}

//_____________________________________________________________________________
AliAlignObj::AliAlignObj(const AliAlignObj& theAlignObj) :
  TObject(theAlignObj),
  fVolPath(theAlignObj.GetSymName()),
  fVolUID(theAlignObj.GetVolUID())
{
  //copy constructor
  for(Int_t i=0; i<6; i++) fDiag[i]=theAlignObj.fDiag[i];
  for(Int_t i=0; i<21; i++) fODia[i]=theAlignObj.fODia[i];
}

//_____________________________________________________________________________
AliAlignObj &AliAlignObj::operator =(const AliAlignObj& theAlignObj)
{
  // assignment operator
  if(this==&theAlignObj) return *this;
  fVolPath = theAlignObj.GetSymName();
  fVolUID = theAlignObj.GetVolUID();
  for(Int_t i=0; i<6; i++) fDiag[i]=theAlignObj.fDiag[i];
  for(Int_t i=0; i<21; i++) fODia[i]=theAlignObj.fODia[i];
  return *this;
}

//_____________________________________________________________________________
AliAlignObj &AliAlignObj::operator*=(const AliAlignObj& theAlignObj)
{
  // multiplication operator
  // The operator can be used to 'combine'
  // two alignment objects
  TGeoHMatrix m1;
  GetMatrix(m1);
  TGeoHMatrix m2;
  theAlignObj.GetMatrix(m2);
  m1.MultiplyLeft(&m2);
  SetMatrix(m1);
  // temporary solution: consider parameters indipendent 
  for(Int_t i=0; i<6; i++)  fDiag[i] = TMath::Sqrt((fDiag[i]*fDiag[i])+(theAlignObj.fDiag[i]*theAlignObj.fDiag[i]));
  return *this;
}

//_____________________________________________________________________________
AliAlignObj::~AliAlignObj()
{
  // dummy destructor
}

//_____________________________________________________________________________
void AliAlignObj::SetVolUID(AliGeomManager::ELayerID detId, Int_t modId)
{
  // From detector name and module number (according to detector numbering)
  // build fVolUID, unique numerical identity of that volume inside ALICE
  // fVolUID is 16 bits, first 5 reserved for detID (32 possible values),
  // remaining 11 for module ID inside det (2048 possible values).
  //
  fVolUID = AliGeomManager::LayerToVolUID(detId,modId);
}

//_____________________________________________________________________________
void AliAlignObj::GetVolUID(AliGeomManager::ELayerID &layerId, Int_t &modId) const
{
  // From the fVolUID, unique numerical identity of that volume inside ALICE,
  // (voluid is 16 bits, first 5 reserved for layerID (32 possible values),
  // remaining 11 for module ID inside det (2048 possible values)), sets
  // the argument layerId to the identity of the layer to which that volume
  // belongs and sets the argument modId to the identity of that volume
  // internally to the layer.
  //
  layerId = AliGeomManager::VolUIDToLayer(fVolUID,modId);
}

//_____________________________________________________________________________
Bool_t AliAlignObj::GetPars(Double_t tr[], Double_t angles[]) const
{
  GetTranslation(tr);
  return GetAngles(angles);
}

//_____________________________________________________________________________
Int_t AliAlignObj::GetLevel() const
{
  // Return the geometry level of the alignable volume to which
  // the alignment object is associated; this is the number of
  // slashes in the corresponding volume path
  //
  if(!gGeoManager){
    AliWarning("gGeoManager doesn't exist or it is still opened: unable to return meaningful level value.");
    return (-1);
  }
  const char* symname = GetSymName();
  const char* path;
  TGeoPNEntry* pne = gGeoManager->GetAlignableEntry(symname);
  if(pne){
    path = pne->GetTitle();
  }else{
    path = symname;
  }

  TString pathStr = path;
  if(pathStr[0]!='/') pathStr.Prepend('/');
  return pathStr.CountChar('/');
}

//_____________________________________________________________________________
Int_t AliAlignObj::Compare(const TObject *obj) const
{
  // Compare the levels of two
  // alignment objects
  // Used in the sorting during
  // the application of alignment
  // objects to the geometry
  //
  Int_t level = GetLevel();
  Int_t level2 = ((AliAlignObj *)obj)->GetLevel();
  if (level == level2)
    return 0;
  else
    return ((level > level2) ? 1 : -1);
}

//______________________________________________________________________________
void AliAlignObj::GetCovMatrix(Double_t *cmat) const
{
  // Fills the cmat argument with the coefficients of the external cov matrix (21 elements)
  // calculating them from the correlation matrix data member
  //

  for(Int_t i=0; i<6; ++i) {
    // Off diagonal elements
    for(Int_t j=0; j<i; ++j) {
      cmat[i*(i+1)/2+j] = (fDiag[j] >= 0. && fDiag[i] >= 0.) ? fODia[(i-1)*i/2+j]*fDiag[j]*fDiag[i]: -999.;
    }

    // Diagonal elements
    cmat[i*(i+1)/2+i] = (fDiag[i] >= 0.) ? fDiag[i]*fDiag[i] : -999.;
  }

  return;
}

//______________________________________________________________________________
void AliAlignObj::GetCovMatrix(TMatrixDSym& mcov) const
{
  // Fills the matrix m passed as argument as the covariance matrix calculated
  // from the coefficients of the reduced covariance matrix data members
  //

  for(Int_t i=0; i<6; ++i) {
    // Off diagonal elements
    for(Int_t j=0; j<i; ++j) {
      mcov(j,i) = mcov(i,j) = (fDiag[j] >= 0. && fDiag[i] >= 0.) ? fODia[(i-1)*i/2+j]*fDiag[j]*fDiag[i]: -999.;
    }

    // Diagonal elements
    mcov(i,i) = (fDiag[i] >= 0.) ? fDiag[i]*fDiag[i] : -999.;
  }

}

//______________________________________________________________________________
void AliAlignObj::SetCorrMatrix(Double_t *cmat)
{
  // Sets the correlation matrix data member from the coefficients of the external covariance
  // matrix (21 elements passed as argument). 
  //
  if(cmat) {

    // Diagonal elements first
    for(Int_t i=0; i<6; ++i) {
      fDiag[i] = (cmat[i*(i+1)/2+i] >= 0.) ? TMath::Sqrt(cmat[i*(i+1)/2+i]) : -999.;
    }

    // ... then the ones off diagonal
    for(Int_t i=0; i<6; ++i)
      // Off diagonal elements
      for(Int_t j=0; j<i; ++j) {
	fODia[(i-1)*i/2+j] = (fDiag[i] > 0. && fDiag[j] > 0.) ? cmat[i*(i+1)/2+j]/(fDiag[j]*fDiag[i]) : 0.;       // check for division by zero (due to diagonal element of 0) and for fDiag != -999. (due to negative input diagonal element).
	if (fODia[(i-1)*i/2+j]>1.)  fODia[(i-1)*i/2+j] =  1.; // check upper boundary
	if (fODia[(i-1)*i/2+j]<-1.) fODia[(i-1)*i/2+j] = -1.; // check lower boundary
      }
  } else {
    for(Int_t i=0; i< 6; ++i) fDiag[i]=-999.;
    for(Int_t i=0; i< 6*(6-1)/2; ++i) fODia[i]=0.;
  }

  return;
}

//______________________________________________________________________________
void AliAlignObj::SetCorrMatrix(TMatrixDSym& mcov)
{
  // Sets the correlation matrix data member from the covariance matrix mcov passed
  // passed as argument. 
  //
  if(mcov.IsValid()) {

    // Diagonal elements first
    for(Int_t i=0; i<6; ++i) {
      fDiag[i] = (mcov(i,i) >= 0.) ? TMath::Sqrt(mcov(i,i)) : -999.;
    }

    // ... then the ones off diagonal
    for(Int_t i=0; i<6; ++i)
      // Off diagonal elements
      for(Int_t j=0; j<i; ++j) {
	fODia[(i-1)*i/2+j] = (fDiag[i] > 0. && fDiag[j] > 0.) ? mcov(i,j)/(fDiag[j]*fDiag[i]) : 0.;       // check for division by zero (due to diagonal element of 0) and for fDiag != -999. (due to negative input diagonal element).
	if (fODia[(i-1)*i/2+j]>1.)  fODia[(i-1)*i/2+j] =  1.; // check upper boundary
	if (fODia[(i-1)*i/2+j]<-1.) fODia[(i-1)*i/2+j] = -1.; // check lower boundary
      }
  } else {
    for(Int_t i=0; i< 6; ++i) fDiag[i]=-999.;
    for(Int_t i=0; i< 6*(6-1)/2; ++i) fODia[i]=0.;
  }

  return;
}

//_____________________________________________________________________________
void AliAlignObj::AnglesToMatrix(const Double_t *angles, Double_t *rot) const
{
  // Calculates the rotation matrix using the 
  // Euler angles in "x y z" notation
  //
  Double_t degrad = TMath::DegToRad();
  Double_t sinpsi = TMath::Sin(degrad*angles[0]);
  Double_t cospsi = TMath::Cos(degrad*angles[0]);
  Double_t sinthe = TMath::Sin(degrad*angles[1]);
  Double_t costhe = TMath::Cos(degrad*angles[1]);
  Double_t sinphi = TMath::Sin(degrad*angles[2]);
  Double_t cosphi = TMath::Cos(degrad*angles[2]);

  rot[0] =  costhe*cosphi;
  rot[1] = -costhe*sinphi;
  rot[2] =  sinthe;
  rot[3] =  sinpsi*sinthe*cosphi + cospsi*sinphi;
  rot[4] = -sinpsi*sinthe*sinphi + cospsi*cosphi;
  rot[5] = -costhe*sinpsi;
  rot[6] = -cospsi*sinthe*cosphi + sinpsi*sinphi;
  rot[7] =  cospsi*sinthe*sinphi + sinpsi*cosphi;
  rot[8] =  costhe*cospsi;
}

//_____________________________________________________________________________
Bool_t AliAlignObj::MatrixToAngles(const Double_t *rot, Double_t *angles) const
{
  // Calculates the Euler angles in "x y z" notation
  // using the rotation matrix
  // Returns false in case the rotation angles can not be
  // extracted from the matrix
  //
  if(TMath::Abs(rot[0])<1e-7 || TMath::Abs(rot[8])<1e-7) {
    AliError("Failed to extract roll-pitch-yall angles!");
    return kFALSE;
  }
  Double_t raddeg = TMath::RadToDeg();
  angles[0]=raddeg*TMath::ATan2(-rot[5],rot[8]);
  angles[1]=raddeg*TMath::ASin(rot[2]);
  angles[2]=raddeg*TMath::ATan2(-rot[1],rot[0]);
  return kTRUE;
}

//______________________________________________________________________________
void AliAlignObj::Transform(AliTrackPoint &p) const
{
  // The method transforms the space-point coordinates using the
  // transformation matrix provided by the AliAlignObj
  // The covariance matrix is not affected since we assume
  // that the transformations are sufficiently small
  //
  if (fVolUID != p.GetVolumeID())
    AliWarning(Form("Alignment object ID is not equal to the space-point ID (%d != %d)",fVolUID,p.GetVolumeID())); 

  TGeoHMatrix m;
  GetMatrix(m);
  Double_t *rot = m.GetRotationMatrix();
  Double_t *tr  = m.GetTranslation();

  Float_t xyzin[3],xyzout[3];
  p.GetXYZ(xyzin);
  for (Int_t i = 0; i < 3; i++)
    xyzout[i] = tr[i]+
                xyzin[0]*rot[3*i]+
                xyzin[1]*rot[3*i+1]+
                xyzin[2]*rot[3*i+2];
  p.SetXYZ(xyzout);
  
}

//_____________________________________________________________________________
void AliAlignObj::Transform(AliTrackPointArray &array) const
{
  // This method is used to transform all the track points
  // from the input AliTrackPointArray
  // 
  AliTrackPoint p;
  for (Int_t i = 0; i < array.GetNPoints(); i++) {
    array.GetPoint(p,i);
    Transform(p);
    array.AddPoint(i,&p);
  }
}

//_____________________________________________________________________________
void AliAlignObj::Print(Option_t *) const
{
  // Print the contents of the
  // alignment object in angles and
  // matrix representations
  //
  Double_t tr[3];
  GetTranslation(tr);
  Double_t angles[3];
  GetAngles(angles);
  TGeoHMatrix m;
  GetMatrix(m);
  const Double_t *rot = m.GetRotationMatrix();

  printf("Volume=%s\n",GetSymName());
  if (GetVolUID() != 0) {
    AliGeomManager::ELayerID layerId;
    Int_t modId;
    GetVolUID(layerId,modId);
    printf("VolumeID=%d LayerID=%d ( %s ) ModuleID=%d\n", GetVolUID(),layerId,AliGeomManager::LayerName(layerId),modId);
  }
  printf("%12.8f%12.8f%12.8f    Tx = %12.8f    Psi   = %12.8f\n", rot[0], rot[1], rot[2], tr[0], angles[0]);
  printf("%12.8f%12.8f%12.8f    Ty = %12.8f    Theta = %12.8f\n", rot[3], rot[4], rot[5], tr[1], angles[1]);
  printf("%12.8f%12.8f%12.8f    Tz = %12.8f    Phi   = %12.8f\n", rot[6], rot[7], rot[8], tr[2], angles[2]);

}

//_____________________________________________________________________________
void AliAlignObj::SetPars(Double_t x, Double_t y, Double_t z,
			  Double_t psi, Double_t theta, Double_t phi)
{
  // Set the global delta transformation by passing 3 angles (expressed in
  // degrees) and 3 shifts (in centimeters)
  // 
  SetTranslation(x,y,z);
  SetRotation(psi,theta,phi);
}

//_____________________________________________________________________________
Bool_t AliAlignObj::SetLocalPars(Double_t x, Double_t y, Double_t z,
				 Double_t psi, Double_t theta, Double_t phi)
{
  // Set the global delta transformation by passing the parameters
  // for the local delta transformation (3 shifts and 3 angles).
  // In case that the TGeo was not initialized or not closed,
  // returns false and the object parameters are not set.
  //
  TGeoHMatrix m;
  Double_t tr[3] = {x, y, z};
  m.SetTranslation(tr);
  Double_t angles[3] = {psi, theta, phi};
  Double_t rot[9];
  AnglesToMatrix(angles,rot);
  m.SetRotation(rot);

  return SetLocalMatrix(m);

}

//_____________________________________________________________________________
Bool_t AliAlignObj::SetLocalTranslation(Double_t x, Double_t y, Double_t z)
{
  // Set the global delta transformation by passing the three shifts giving
  // the translation in the local reference system of the alignable
  // volume (known by TGeo geometry).
  // In case that the TGeo was not initialized or not closed,
  // returns false and the object parameters are not set.
  //
  TGeoHMatrix m;
  Double_t tr[3] = {x, y, z};
  m.SetTranslation(tr);

  return SetLocalMatrix(m);

}

//_____________________________________________________________________________
Bool_t AliAlignObj::SetLocalTranslation(const TGeoMatrix& m)
{
  // Set the global delta transformation by passing the matrix of
  // the local delta transformation and taking its translational part
  // In case that the TGeo was not initialized or not closed,
  // returns false and the object parameters are not set.
  //
  const Double_t* tr = m.GetTranslation();
  TGeoHMatrix mtr;
  mtr.SetTranslation(tr);

  return SetLocalMatrix(mtr);

}

//_____________________________________________________________________________
Bool_t AliAlignObj::SetLocalRotation(Double_t psi, Double_t theta, Double_t phi)
{
  // Set the global delta transformation by passing the three angles giving
  // the rotation in the local reference system of the alignable
  // volume (known by TGeo geometry).
  // In case that the TGeo was not initialized or not closed,
  // returns false and the object parameters are not set.
  //
  TGeoHMatrix m;
  Double_t angles[3] = {psi, theta, phi};
  Double_t rot[9];
  AnglesToMatrix(angles,rot);
  m.SetRotation(rot);

  return SetLocalMatrix(m);

}

//_____________________________________________________________________________
Bool_t AliAlignObj::SetLocalRotation(const TGeoMatrix& m)
{
  // Set the global delta transformation by passing the matrix of
  // the local delta transformation and taking its rotational part
  // In case that the TGeo was not initialized or not closed,
  // returns false and the object parameters are not set.
  //
  TGeoHMatrix rotm;
  const Double_t* rot = m.GetRotationMatrix();
  rotm.SetRotation(rot);

  return SetLocalMatrix(rotm);

}

//_____________________________________________________________________________
Bool_t AliAlignObj::SetLocalMatrix(const TGeoMatrix& m)
{
  // Set the global delta transformation by passing the TGeo matrix
  // for the local delta transformation.
  // In case that the TGeo was not initialized or not closed,
  // returns false and the object parameters are not set.
  //
  if (!gGeoManager || !gGeoManager->IsClosed()) {
    AliError("Can't set the alignment object parameters! gGeoManager doesn't exist or it is still opened!");
    return kFALSE;
  }

  const char* symname = GetSymName();
  TGeoPhysicalNode* node;
  TGeoPNEntry* pne = gGeoManager->GetAlignableEntry(symname);
  if(pne){
    node = gGeoManager->MakeAlignablePN(pne);
  }else{
    AliWarning(Form("The symbolic volume name %s does not correspond to a physical entry. Using it as volume path!",symname));
    node = (TGeoPhysicalNode*) gGeoManager->MakePhysicalNode(symname);
  }

  if (!node) {
    AliError(Form("Volume name or path %s not valid!",symname));
    return kFALSE;
  }
  if (node->IsAligned())
    AliWarning(Form("Volume %s has been already misaligned!",symname));

  TGeoHMatrix m1;
  const Double_t *tr = m.GetTranslation();
  m1.SetTranslation(tr);
  const Double_t* rot = m.GetRotationMatrix();
  m1.SetRotation(rot);

  TGeoHMatrix align,gprime,gprimeinv;
  gprime = *node->GetMatrix();
  gprimeinv = gprime.Inverse();
  m1.Multiply(&gprimeinv);
  m1.MultiplyLeft(&gprime);

  return SetMatrix(m1);
}

//_____________________________________________________________________________
Bool_t AliAlignObj::SetMatrix(const TGeoMatrix& m)
{
  // Set the global delta transformation by passing the TGeoMatrix
  // for it
  //
  SetTranslation(m);
  return SetRotation(m);
}

//_____________________________________________________________________________
Bool_t AliAlignObj::GetLocalPars(Double_t transl[], Double_t angles[]) const
{
  // Get the translations and angles (in degrees) expressing the
  // local delta transformation.
  // In case that the TGeo was not initialized or not closed,
  // returns false and the object parameters are not set.
  //
  if(!GetLocalTranslation(transl)) return kFALSE;
  return GetLocalAngles(angles);
}

//_____________________________________________________________________________
Bool_t AliAlignObj::GetLocalTranslation(Double_t* tr) const
{
  // Get the 3 shifts giving the translational part of the local
  // delta transformation.
  // In case that the TGeo was not initialized or not closed,
  // returns false and the object parameters are not set.
  //
  TGeoHMatrix ml;
  if(!GetLocalMatrix(ml)) return kFALSE;
  const Double_t* transl;
  transl = ml.GetTranslation();
  tr[0]=transl[0];
  tr[1]=transl[1];
  tr[2]=transl[2];
  return kTRUE;
}

//_____________________________________________________________________________
Bool_t AliAlignObj::GetLocalAngles(Double_t* angles) const
{
  // Get the 3 angles giving the rotational part of the local
  // delta transformation.
  // In case that the TGeo was not initialized or not closed,
  // returns false and the object parameters are not set.
  //
  TGeoHMatrix ml;
  if(!GetLocalMatrix(ml)) return kFALSE;
  const Double_t *rot = ml.GetRotationMatrix();
  return MatrixToAngles(rot,angles);
}

//_____________________________________________________________________________
Bool_t AliAlignObj::GetLocalMatrix(TGeoHMatrix& m) const
{
  // Get the matrix for the local delta transformation.
  // In case that the TGeo was not initialized or not closed,
  // returns false and the object parameters are not set.
  //
  if (!gGeoManager || !gGeoManager->IsClosed()) {
    AliError("Can't set the alignment object parameters! gGeoManager doesn't exist or it is still opened!");
    return kFALSE;
  }

  const char* symname = GetSymName();
  TGeoPhysicalNode* node;
  TGeoPNEntry* pne = gGeoManager->GetAlignableEntry(symname);
  if(pne){
    node = gGeoManager->MakeAlignablePN(pne);
  }else{
    AliWarning(Form("The symbolic volume name %s does not correspond to a physical entry. Using it as volume path!",symname));
    node = (TGeoPhysicalNode*) gGeoManager->MakePhysicalNode(symname);
  }

  if (!node) {
    AliError(Form("Volume name or path %s not valid!",symname));
    return kFALSE;
  }
  if (node->IsAligned())
    AliWarning(Form("Volume %s has been already misaligned!",symname));

  GetMatrix(m);
  TGeoHMatrix gprime,gprimeinv;
  gprime = *node->GetMatrix();
  gprimeinv = gprime.Inverse();
  m.Multiply(&gprime);
  m.MultiplyLeft(&gprimeinv);

  return kTRUE;
}

//_____________________________________________________________________________
Bool_t AliAlignObj::ApplyToGeometry()
{
  // Apply the current alignment object to the TGeo geometry
  // This method returns FALSE if the symname of the object was not
  // valid neither to get a TGeoPEntry nor as a volume path
  //
  if (!gGeoManager || !gGeoManager->IsClosed()) {
    AliError("Can't apply the alignment object! gGeoManager doesn't exist or it is still opened!");
    return kFALSE;
  }
  
  const char* symname = GetSymName();
  const char* path;
  TGeoPhysicalNode* node;
  TGeoPNEntry* pne = gGeoManager->GetAlignableEntry(symname);
  if(pne){
    path = pne->GetTitle();
    node = gGeoManager->MakeAlignablePN(pne);
  }else{
    AliDebug(1,Form("The symbolic volume name %s does not correspond to a physical entry. Using it as a volume path!",symname));
    path=symname;
    if (!gGeoManager->CheckPath(path)) {
      AliDebug(1,Form("Volume path %s not valid!",path));
      return kFALSE;
    }
    if (gGeoManager->GetListOfPhysicalNodes()->FindObject(path)) {
      AliError(Form("Volume %s has already been misaligned!",path));
      return kFALSE;
    }
    node = (TGeoPhysicalNode*) gGeoManager->MakePhysicalNode(path);
  }

  if (!node) {
    AliError(Form("Volume path %s not valid!",path));
    return kFALSE;
  }

  TGeoHMatrix align,gprime;
  gprime = *node->GetMatrix();
  GetMatrix(align);
  gprime.MultiplyLeft(&align);
  TGeoHMatrix *ginv = new TGeoHMatrix;
  TGeoHMatrix *g = node->GetMatrix(node->GetLevel()-1);
  *ginv = g->Inverse();
  *ginv *= gprime;
  AliGeomManager::ELayerID layerId; // unique identity for layer in the alobj
  Int_t modId; // unique identity for volume inside layer in the alobj
  GetVolUID(layerId, modId);
  AliDebug(2,Form("Aligning volume %s of detector layer %d with local ID %d",symname,layerId,modId));
  node->Align(ginv);

  return kTRUE;
}

