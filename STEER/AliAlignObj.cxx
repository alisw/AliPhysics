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
//  (AliAlignObjAngles, AliAlignObjMatrix) are derived in separate files.
//-----------------------------------------------------------------
#include <TGeoManager.h>
#include <TGeoPhysicalNode.h>
#include <TMath.h>

#include "AliAlignObj.h"
#include "AliTrackPointArray.h"
#include "AliLog.h"
#include "AliAlignObjAngles.h"
 
ClassImp(AliAlignObj)

Int_t AliAlignObj::fgLayerSize[kLastLayer - kFirstLayer] = {
  80, 160,  // ITS SPD first and second layer
  84, 176,  // ITS SDD first and second layer
  748, 950, // ITS SSD first and second layer
  36, 36,   // TPC inner and outer chambers
  90, 90, 90, 90, 90, 90,  // 6 TRD chambers' layers
  1638,     // TOF
  1, 1,     // PHOS ??
  7,        // HMPID ??
  1         // MUON ??
};

const char* AliAlignObj::fgLayerName[kLastLayer - kFirstLayer] = {
  "ITS inner pixels layer", "ITS outer pixels layer",
  "ITS inner drifts layer", "ITS outer drifts layer",
  "ITS inner strips layer", "ITS outer strips layer",
  "TPC inner chambers layer", "TPC outer chambers layer",
  "TRD chambers layer 1", "TRD chambers layer 2", "TRD chambers layer 3",
  "TRD chambers layer 4", "TRD chambers layer 5", "TRD chambers layer 6",
  "TOF layer",
  "?","?",
  "HMPID layer",
  "?"
};

TString* AliAlignObj::fgVolPath[kLastLayer - kFirstLayer] = {
  0x0,0x0,
  0x0,0x0,
  0x0,0x0,
  0x0,0x0,
  0x0,0x0,0x0,
  0x0,0x0,0x0,
  0x0,
  0x0,0x0,
  0x0,
  0x0
};

AliAlignObj** AliAlignObj::fgAlignObjs[kLastLayer - kFirstLayer] = {
  0x0,0x0,
  0x0,0x0,
  0x0,0x0,
  0x0,0x0,
  0x0,0x0,0x0,
  0x0,0x0,0x0,
  0x0,
  0x0,0x0,
  0x0,
  0x0
};

//_____________________________________________________________________________
AliAlignObj::AliAlignObj():
  fVolPath(),
  fVolUID(0)
{
  // default constructor
  InitSymNames();
}

//_____________________________________________________________________________
AliAlignObj::AliAlignObj(const char* symname, UShort_t voluid) :
  TObject(),
  fVolPath(symname),
  fVolUID(voluid)
{
  // standard constructor
  //
}

//_____________________________________________________________________________
AliAlignObj::AliAlignObj(const AliAlignObj& theAlignObj) :
  TObject(theAlignObj),
  fVolPath(theAlignObj.GetSymName()),
  fVolUID(theAlignObj.GetVolUID())
{
  //copy constructor
}

//_____________________________________________________________________________
AliAlignObj &AliAlignObj::operator =(const AliAlignObj& theAlignObj)
{
  // assignment operator
  if(this==&theAlignObj) return *this;
  fVolPath = theAlignObj.GetSymName();
  fVolUID = theAlignObj.GetVolUID();
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
  return *this;
}

//_____________________________________________________________________________
AliAlignObj::~AliAlignObj()
{
  // dummy destructor
}

//_____________________________________________________________________________
void AliAlignObj::SetVolUID(ELayerID detId, Int_t modId)
{
  // From detector name and module number (according to detector numbering)
  // build fVolUID, unique numerical identity of that volume inside ALICE
  // fVolUID is 16 bits, first 5 reserved for detID (32 possible values),
  // remaining 11 for module ID inside det (2048 possible values).
  //
  fVolUID = LayerToVolUID(detId,modId);
}

//_____________________________________________________________________________
void AliAlignObj::GetVolUID(ELayerID &layerId, Int_t &modId) const
{
  // From the fVolUID, unique numerical identity of that volume inside ALICE,
  // (voluid is 16 bits, first 5 reserved for layerID (32 possible values),
  // remaining 11 for module ID inside det (2048 possible values)), sets
  // the argument layerId to the identity of the layer to which that volume
  // belongs and sets the argument modId to the identity of that volume
  // internally to the layer.
  //
  layerId = VolUIDToLayer(fVolUID,modId);
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

  TString path_str = path;
  if(path_str[0]!='/') path_str.Prepend('/');
  return path_str.CountChar('/');
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
    ELayerID layerId;
    Int_t modId;
    GetVolUID(layerId,modId);
    printf("VolumeID=%d LayerID=%d ( %s ) ModuleID=%d\n", GetVolUID(),layerId,LayerName(layerId),modId);
  }
  printf("%12.8f%12.8f%12.8f    Tx = %12.8f    Psi   = %12.8f\n", rot[0], rot[1], rot[2], tr[0], angles[0]);
  printf("%12.8f%12.8f%12.8f    Ty = %12.8f    Theta = %12.8f\n", rot[3], rot[4], rot[5], tr[1], angles[1]);
  printf("%12.8f%12.8f%12.8f    Tz = %12.8f    Phi   = %12.8f\n", rot[6], rot[7], rot[8], tr[2], angles[2]);

}

//_____________________________________________________________________________
Int_t AliAlignObj::LayerSize(Int_t layerId)
{
  // Get the layer size for layer corresponding to layerId.
  // Implemented only for ITS,TPC,TRD,TOF and HMPID
  //
  if (layerId < kFirstLayer || layerId >= kLastLayer) {
    AliErrorClass(Form("Invalid layer index %d ! Layer range is (%d -> %d) !",layerId,kFirstLayer,kLastLayer));
    return 0;
  }
  else {
    return fgLayerSize[layerId - kFirstLayer];
 }
}

//_____________________________________________________________________________
const char* AliAlignObj::LayerName(Int_t layerId)
{
  // Get the layer name corresponding to layerId.
  // Implemented only for ITS,TPC,TRD,TOF and HMPID
  //
  if (layerId < kFirstLayer || layerId >= kLastLayer) {
    AliErrorClass(Form("Invalid layer index %d ! Layer range is (%d -> %d) !",layerId,kFirstLayer,kLastLayer));
    return "Invalid Layer!";
  }
  else {
    return fgLayerName[layerId - kFirstLayer];
 }
}

//_____________________________________________________________________________
UShort_t AliAlignObj::LayerToVolUID(ELayerID layerId, Int_t modId)
{
  // From detector (layer) name and module number (according to detector
  // internal numbering) build the unique numerical identity of that volume
  // inside ALICE
  // fVolUID is 16 bits, first 5 reserved for layerID (32 possible values),
  // remaining 11 for module ID inside det (2048 possible values).
  //
  return ((UShort_t(layerId) << 11) | UShort_t(modId));
}

//_____________________________________________________________________________
UShort_t AliAlignObj::LayerToVolUID(Int_t   layerId, Int_t modId)
{
  // From detector (layer) name and module number (according to detector
  // internal numbering) build the unique numerical identity of that volume
  // inside ALICE
  // fVolUID is 16 bits, first 5 reserved for layerID (32 possible values),
  // remaining 11 for module ID inside det (2048 possible values).
  //
  return ((UShort_t(layerId) << 11) | UShort_t(modId));
}

//_____________________________________________________________________________
AliAlignObj::ELayerID AliAlignObj::VolUIDToLayer(UShort_t voluid, Int_t &modId)
{
  // From voluid, unique numerical identity of that volume inside ALICE,
  // (voluid is 16 bits, first 5 reserved for layerID (32 possible values),
  // remaining 11 for module ID inside det (2048 possible values)), return
  // the identity of the layer to which that volume belongs and sets the
  // argument modId to the identity of that volume internally to the layer.
  //
  modId = voluid & 0x7ff;

  return VolUIDToLayer(voluid);
}

//_____________________________________________________________________________
AliAlignObj::ELayerID AliAlignObj::VolUIDToLayer(UShort_t voluid)
{
  // From voluid, unique numerical identity of that volume inside ALICE,
  // (voluid is 16 bits, first 5 reserved for layerID (32 possible values),
  // remaining 11 for module ID inside det (2048 possible values)), return
  // the identity of the layer to which that volume belongs
  //
  return ELayerID((voluid >> 11) & 0x1f);
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
    AliWarning(Form("The symbolic volume name %s does not correspond to a physical entry. Using it as a volume path!",symname));
    path=symname;
    if (!gGeoManager->CheckPath(path)) {
      AliError(Form("Volume path %s not valid!",path));
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
  AliAlignObj::ELayerID layerId; // unique identity for layer in the alobj
  Int_t modId; // unique identity for volume inside layer in the alobj
  GetVolUID(layerId, modId);
  AliDebug(2,Form("Aligning volume %s of detector layer %d with local ID %d",symname,layerId,modId));
  node->Align(ginv);

  return kTRUE;
}

//_____________________________________________________________________________
Bool_t AliAlignObj::GetFromGeometry(const char *symname, AliAlignObj &alobj)
{
  // Get the alignment object which corresponds to the symbolic volume name
  // symname (in case equal to the TGeo volume path)
  // The method is extremely slow due to the searching by string.
  // Therefore it should be used with great care!!
  // This method returns FALSE if the symname of the object was not
  // valid neither to get a TGeoPEntry nor as a volume path, or if the path
  // associated to the TGeoPNEntry was not valid.
  //

  // Reset the alignment object
  alobj.SetPars(0,0,0,0,0,0);
  alobj.SetSymName(symname);

  if (!gGeoManager || !gGeoManager->IsClosed()) {
    AliErrorClass("Can't get the alignment object! gGeoManager doesn't exist or it is still opened!");
    return kFALSE;
  }

  if (!gGeoManager->GetListOfPhysicalNodes()) {
    AliErrorClass("Can't get the alignment object! gGeoManager doesn't contain any aligned nodes!");
    return kFALSE;
  }

  const char *path;
  TGeoPNEntry* pne = gGeoManager->GetAlignableEntry(symname);
  if(pne){
    path = pne->GetTitle();
  }else{
    AliWarningClass(Form("The symbolic volume name %s does not correspond to a physical entry. Using it as a volume path!",symname));
    path = symname;
  }
  TObjArray* nodesArr = gGeoManager->GetListOfPhysicalNodes();
  TGeoPhysicalNode* node = NULL;
  for (Int_t iNode = 0; iNode < nodesArr->GetEntriesFast(); iNode++) {
    TGeoPhysicalNode* tempNode = (TGeoPhysicalNode*) nodesArr->UncheckedAt(iNode);
    const char *nodePath = tempNode->GetName();
    if (strcmp(symname,nodePath) == 0) {
      node = tempNode;
      break;
    }
  }

  if (!node) {
    if (!gGeoManager->cd(symname)) {
      AliErrorClass(Form("%s not valid neither as symbolic volume name nor as volume path!",symname));
      return kFALSE;
    }
    else {
      AliWarningClass(Form("Volume (%s) has not been misaligned!",symname));
      return kTRUE;
    }
  }

  TGeoHMatrix align,gprime,g,ginv,l;
  gprime = *node->GetMatrix();
  l = *node->GetOriginalMatrix();
  g = *node->GetMatrix(node->GetLevel()-1);
  g *= l;
  ginv = g.Inverse();
  align = gprime * ginv;

  return alobj.SetMatrix(align);
}

//_____________________________________________________________________________
void  AliAlignObj::InitAlignObjFromGeometry()
{
  // Loop over all alignable volumes and extract
  // the corresponding alignment objects from
  // the TGeo geometry

  if(fgAlignObjs[0]) return;
  
  InitSymNames();

  for (Int_t iLayer = kFirstLayer; iLayer < AliAlignObj::kLastLayer; iLayer++) {
    fgAlignObjs[iLayer-kFirstLayer] = new AliAlignObj*[AliAlignObj::LayerSize(iLayer)];
    for (Int_t iModule = 0; iModule < AliAlignObj::LayerSize(iLayer); iModule++) {
      UShort_t volid = AliAlignObj::LayerToVolUID(iLayer,iModule);
      fgAlignObjs[iLayer-kFirstLayer][iModule] = new AliAlignObjAngles("",volid,0,0,0,0,0,0,kTRUE);
      const char *symname = SymName(volid);
      if (!GetFromGeometry(symname, *fgAlignObjs[iLayer-kFirstLayer][iModule]))
	AliErrorClass(Form("Failed to extract the alignment object for the volume (ID=%d and path=%s) !",volid,symname));
    }
  }
  
}

//_____________________________________________________________________________
AliAlignObj* AliAlignObj::GetAlignObj(UShort_t voluid) {
  // Returns the alignment object for given volume ID
  //
  Int_t modId;
  ELayerID layerId = VolUIDToLayer(voluid,modId);
  return GetAlignObj(layerId,modId);
}

//_____________________________________________________________________________
AliAlignObj* AliAlignObj::GetAlignObj(ELayerID layerId, Int_t modId)
{
  // Returns pointer to alignment object given its layer and module ID
  //
  if(modId<0 || modId>=fgLayerSize[layerId-kFirstLayer]){
    AliWarningClass(Form("Module number %d not in the valid range (0->%d) !",modId,fgLayerSize[layerId-kFirstLayer]-1));
    return NULL;
  }
  InitAlignObjFromGeometry();

  return fgAlignObjs[layerId-kFirstLayer][modId];
}

//_____________________________________________________________________________
const char* AliAlignObj::SymName(UShort_t voluid) {
  // Returns the symbolic volume name for given volume ID
  //
  Int_t modId;
  ELayerID layerId = VolUIDToLayer(voluid,modId);
  return SymName(layerId,modId);
}

//_____________________________________________________________________________
const char* AliAlignObj::SymName(ELayerID layerId, Int_t modId)
{
  // Returns the symbolic volume name given for a given layer
  // and module ID
  //
  if(modId<0 || modId>=fgLayerSize[layerId-kFirstLayer]){
    AliWarningClass(Form("Module number %d not in the valid range (0->%d) !",modId,fgLayerSize[layerId-kFirstLayer]-1));
    return NULL;
  }
  InitSymNames();

  return fgVolPath[layerId-kFirstLayer][modId].Data();
}

//_____________________________________________________________________________
void AliAlignObj::InitSymNames()
{
  // Initialize the LUTs which associate the symbolic volume names
  // for each alignable volume with their unique numerical identity.
  // The LUTs are static, so they are created during the instantiation
  // of the first intance of AliAlignObj
  //
  if (fgVolPath[0]) return;

  for (Int_t iLayer = 0; iLayer < (kLastLayer - kFirstLayer); iLayer++)
    fgVolPath[iLayer] = new TString[fgLayerSize[iLayer]];

  TString symname;
  Int_t modnum; // in the following, set it to 0 at the start of each layer

  /*********************       ITS layers  ***********************/
  TString strSPD = "ITS/SPD";
  TString strSDD = "ITS/SDD";
  TString strSSD = "ITS/SSD";
  TString strStave = "/Stave";
  TString strLadder = "/Ladder";
  TString strSector = "/Sector";
  TString strSensor = "/Sensor";
  TString strEntryName1;
  TString strEntryName2;

  /*********************       SPD layer1  ***********************/
  {
    modnum = 0;

    for(Int_t c1 = 1; c1<=10; c1++){
      strEntryName1 = strSPD;
      strEntryName1 += 0;
      strEntryName1 += strSector;
      strEntryName1 += (c1-1);
      for(Int_t c2 =1; c2<=2; c2++){
	strEntryName2 = strEntryName1;
	strEntryName2 += strStave;
	strEntryName2 += (c2-1);
	for(Int_t c3 =1; c3<=4; c3++){
	  symname = strEntryName2;
	  symname += strLadder;
	  symname += (c3-1);
	  fgVolPath[kSPD1-kFirstLayer][modnum] = symname.Data();
	  modnum++;
	}
      }
    }
  }
  
  /*********************       SPD layer2  ***********************/
  {
    modnum = 0;

    for(Int_t c1 = 1; c1<=10; c1++){
      strEntryName1 = strSPD;
      strEntryName1 += 1;
      strEntryName1 += strSector;
      strEntryName1 += (c1-1);
      for(Int_t c2 =1; c2<=4; c2++){
	strEntryName2 = strEntryName1;
	strEntryName2 += strStave;
	strEntryName2 += (c2-1);
	for(Int_t c3 =1; c3<=4; c3++){
	  symname = strEntryName2;
	  symname += strLadder;
	  symname += (c3-1);
	  fgVolPath[kSPD2-kFirstLayer][modnum] = symname.Data();
	  modnum++;
	}
      }
    }
  }

  /*********************       SDD layer1  ***********************/
  {
    modnum=0;

    for(Int_t c1 = 1; c1<=14; c1++){
      strEntryName1 = strSDD;
      strEntryName1 += 2;
      strEntryName1 +=strLadder;
      strEntryName1 += (c1-1);
      for(Int_t c2 =1; c2<=6; c2++){
	symname = strEntryName1;
	symname += strSensor;
	symname += (c2-1);
	fgVolPath[kSDD1-kFirstLayer][modnum] = symname.Data();
	modnum++;
      }
    }
  }

  /*********************       SDD layer2  ***********************/
  {
    modnum=0;

    for(Int_t c1 = 1; c1<=22; c1++){
      strEntryName1 = strSDD;
      strEntryName1 += 3;
      strEntryName1 +=strLadder;
      strEntryName1 += (c1-1);
      for(Int_t c2 = 1; c2<=8; c2++){
	symname = strEntryName1;
	symname += strSensor;
	symname += (c2-1);
	fgVolPath[kSDD2-kFirstLayer][modnum] = symname.Data();
	modnum++;
      }
    }
  }

  /*********************       SSD layer1  ***********************/
  {
    modnum=0;

    for(Int_t c1 = 1; c1<=34; c1++){
      strEntryName1 = strSSD;
      strEntryName1 += 4;
      strEntryName1 +=strLadder;
      strEntryName1 += (c1-1);
      for(Int_t c2 = 1; c2<=22; c2++){
	symname = strEntryName1;
	symname += strSensor;
	symname += (c2-1);
	fgVolPath[kSSD1-kFirstLayer][modnum] = symname.Data();
	modnum++;
      }
    }
  }

  /*********************       SSD layer2  ***********************/
  {
    modnum=0;

    for(Int_t c1 = 1; c1<=38; c1++){
      strEntryName1 = strSSD;
      strEntryName1 += 5;
      strEntryName1 +=strLadder;
      strEntryName1 += (c1-1);
      for(Int_t c2 = 1; c2<=25; c2++){
	symname = strEntryName1;
	symname += strSensor;
	symname += (c2-1);
	fgVolPath[kSSD2-kFirstLayer][modnum] = symname.Data();
	modnum++;
      }
    }
  }


  /***************    TPC inner and outer layers    ****************/
  TString sAsector="TPC/EndcapA/Sector";
  TString sCsector="TPC/EndcapC/Sector";
  TString sInner="/InnerChamber";
  TString sOuter="/OuterChamber";
  
  /***************    TPC inner chambers' layer    ****************/
  {
    modnum = 0;
    
    for(Int_t cnt=1; cnt<=18; cnt++){
      symname = sAsector;
      symname += cnt;
      symname += sInner;
      fgVolPath[kTPC1-kFirstLayer][modnum] = symname.Data();
      modnum++;
    }
    for(Int_t cnt=1; cnt<=18; cnt++){
      symname = sCsector;
      symname += cnt;
      symname += sInner;
      fgVolPath[kTPC1-kFirstLayer][modnum] = symname.Data();
      modnum++;
    }
  }

  /***************    TPC outer chambers' layer    ****************/
  {
    modnum = 0;
    
    for(Int_t cnt=1; cnt<=18; cnt++){
      symname = sAsector;
      symname += cnt;
      symname += sOuter;
      fgVolPath[kTPC2-kFirstLayer][modnum] = symname.Data();
      modnum++;
    }
    for(Int_t cnt=1; cnt<=18; cnt++){
      symname = sCsector;
      symname += cnt;
      symname += sOuter;
      fgVolPath[kTPC2-kFirstLayer][modnum] = symname.Data();
      modnum++;
    }
  }    

  /*********************       TOF layer   ***********************/
  {
    modnum=0;
    
    Int_t nstrA=15;
    Int_t nstrB=19;
    Int_t nstrC=19;
    Int_t nSectors=18;
    Int_t nStrips=nstrA+2*nstrB+2*nstrC;
    
    TString snSM  = "TOF/sm";
    TString snSTRIP = "/strip";
    
    for (Int_t isect = 0; isect < nSectors; isect++) {
      for (Int_t istr = 1; istr <= nStrips; istr++) {	
	symname  = snSM;
	symname += Form("%02d",isect);
	symname += snSTRIP;
	symname += Form("%02d",istr);
	fgVolPath[kTOF-kFirstLayer][modnum] = symname.Data();	
	modnum++;
      }
    }
  } 

  /*********************      HMPID layer   ***********************/
  {
    TString str = "/HMPID/Chamber";
    TString symname;

    for (modnum=0; modnum < 7; modnum++) {
      symname = str;
      symname += modnum;
      fgVolPath[kHMPID-kFirstLayer][modnum] = symname.Data();
    }
  }

  /*********************      TRD layers 1-6   *******************/
  //!! 6 layers with index increasing in outwards direction
  {
    Int_t arTRDlayId[6] = {kTRD1, kTRD2, kTRD3, kTRD4, kTRD5, kTRD6};

    TString snStr  = "TRD/sm";
    TString snApp1 = "/st";
    TString snApp2 = "/pl";
    
    for(Int_t layer=0; layer<6; layer++){
      modnum=0;
      for (Int_t isect = 0; isect < 18; isect++) {
	for (Int_t icham = 0; icham < 5; icham++) {
	  symname  = snStr;
	  symname += Form("%02d",isect);
	  symname += snApp1;
	  symname += icham;
	  symname += snApp2;
	  symname += layer;
	  fgVolPath[arTRDlayId[layer]-kFirstLayer][modnum] = symname.Data();
	  modnum++;
	}
      }
    }
  }
}

