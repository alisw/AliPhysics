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

//-----------------------------------------------------------------
//   Implementation of the alignment object class through
//   the derived concrete representation of alignment object class:
//   AliAlignObjMatrix derived from the base class AliAlignObj
//-----------------------------------------------------------------

#include "AliAlignObj.h"
#include "AliAlignObjMatrix.h"
//#include "AliLog.h"

ClassImp(AliAlignObjMatrix)

//_____________________________________________________________________________
AliAlignObjMatrix::AliAlignObjMatrix() : 
  AliAlignObj(),
  fMatrix()
{
  // Default constructor
  //
}

//_____________________________________________________________________________
AliAlignObjMatrix::AliAlignObjMatrix(const char* volpath, UShort_t voluid, Double_t x, Double_t y, Double_t z, Double_t psi, Double_t theta, Double_t phi) :
  AliAlignObj(volpath, voluid),
  fMatrix()
{
  // standard constructor with 3 translation + 3 rotation parameters
  //
  SetTranslation(x, y, z);
  SetRotation(psi, theta, phi);
}

//_____________________________________________________________________________
AliAlignObjMatrix::AliAlignObjMatrix(const char* volpath, ELayerID layerId, Int_t volId, Double_t x, Double_t y, Double_t z, Double_t psi, Double_t theta, Double_t phi, Bool_t global) throw (const Char_t *) :
  AliAlignObj(volpath,layerId,volId),
  fMatrix()
{
  // standard constructor with 3 translation + 3 rotation parameters
  // If the user explicitly sets the global variable to kFALSE then the
  // parameters are interpreted as giving the local transformation.
  // This requires to have a gGeoMenager active instance, otherwise the
  // constructor will fail (no object created)
  // 
  if(global){
    SetTranslation(x, y, z);
    SetRotation(psi, theta, phi);
  }else{
    if(!SetLocalPars(x,y,z,psi,theta,phi)) throw "Alignment object creation failed (TGeo instance needed)!\n";
  }
}


//_____________________________________________________________________________
AliAlignObjMatrix::AliAlignObjMatrix(const char* volpath, UShort_t volUID, TGeoMatrix& m) :
  AliAlignObj(volpath,volUID),
  fMatrix()
{
  // standard constructor with TGeoMatrix
  //
  SetTranslation(m);
  SetRotation(m);
}

//_____________________________________________________________________________
AliAlignObjMatrix::AliAlignObjMatrix(const char* volpath, ELayerID layerId, Int_t volId, TGeoMatrix& m) :
  AliAlignObj(volpath,layerId,volId),
  fMatrix()
{
  // standard constructor with TGeoMatrix
  //
  SetTranslation(m);
  SetRotation(m);
}

//_____________________________________________________________________________
AliAlignObjMatrix::AliAlignObjMatrix(const AliAlignObj& theAlignObj) :
  AliAlignObj(theAlignObj),
  fMatrix()
{
  //copy constructor
  //
  Double_t tr[3];
  theAlignObj.GetTranslation(tr);
  SetTranslation(tr[0],tr[1],tr[2]);
  Double_t rot[3];
  theAlignObj.GetAngles(rot);
  SetRotation(rot[0],rot[1],rot[2]);
}

//_____________________________________________________________________________
AliAlignObjMatrix &AliAlignObjMatrix::operator =(const AliAlignObj& theAlignObj)
{  
  // assignment operator
  //
  if(this==&theAlignObj) return *this;
  ((AliAlignObj *)this)->operator=(theAlignObj);
  Double_t tr[3];
  theAlignObj.GetTranslation(tr);
  SetTranslation(tr[0],tr[1],tr[2]);
  Double_t rot[3];
  theAlignObj.GetAngles(rot);
  SetRotation(rot[0],rot[1],rot[2]);
  return *this;
}

//_____________________________________________________________________________
AliAlignObjMatrix::~AliAlignObjMatrix()
{
  // Destructor
  //
}

//_____________________________________________________________________________
void AliAlignObjMatrix::SetTranslation(Double_t x, Double_t y, Double_t z)
{
  Double_t tr[3];
  tr[0]=x; tr[1]=y; tr[2]=z;
  fMatrix.SetTranslation(tr);
}

//_____________________________________________________________________________
void AliAlignObjMatrix::SetTranslation(const TGeoMatrix& m)
{
  const Double_t *tr = m.GetTranslation();
  fMatrix.SetTranslation(tr);
}

//_____________________________________________________________________________
void AliAlignObjMatrix::SetRotation(Double_t psi, Double_t theta, Double_t phi)
{
  Double_t angles[3] = {psi, theta, phi};
  Double_t rot[9];
  AnglesToMatrix(angles,rot);
  fMatrix.SetRotation(rot);
}

//_____________________________________________________________________________
Bool_t AliAlignObjMatrix::SetRotation(const TGeoMatrix& m)
{
  const Double_t* rot = m.GetRotationMatrix();
  fMatrix.SetRotation(rot);
  return kTRUE;
}

//_____________________________________________________________________________
void AliAlignObjMatrix::SetMatrix(const TGeoMatrix& m)
{
  // Set rotation matrix and translation
  // using TGeoMatrix
  SetTranslation(m);
  SetRotation(m);
}

//_____________________________________________________________________________
void AliAlignObjMatrix::SetPars(Double_t x, Double_t y, Double_t z,
		       Double_t psi, Double_t theta, Double_t phi)
{
  // Set rotation matrix and translation
  // using 3 angles and 3 translations
  SetTranslation(x,y,z);
  SetRotation(psi,theta,phi);
}

//_____________________________________________________________________________
void AliAlignObjMatrix::GetTranslation(Double_t *tr) const
{
  // Get Translation from TGeoMatrix
  const Double_t* translation = fMatrix.GetTranslation();
  tr[0] = translation[0];
  tr[1] = translation[1];
  tr[2] = translation[2];
}

//_____________________________________________________________________________
Bool_t AliAlignObjMatrix::GetAngles(Double_t *angles) const
{
  // Get rotation angles from the TGeoHMatrix
  const Double_t* rot = fMatrix.GetRotationMatrix();
  return MatrixToAngles(rot,angles);
}

//_____________________________________________________________________________
void AliAlignObjMatrix::GetPars(Double_t tr[], Double_t angles[]) const
{
  GetTranslation(tr);
  GetAngles(angles);
}

//_____________________________________________________________________________
void AliAlignObjMatrix::GetMatrix(TGeoHMatrix& m) const
{
  // Get TGeoHMatrix
  //
  const Double_t *tr = fMatrix.GetTranslation();
  m.SetTranslation(tr);
  const Double_t *rot = fMatrix.GetRotationMatrix();
  m.SetRotation(rot);
}

//_____________________________________________________________________________
AliAlignObj& AliAlignObjMatrix::Inverse() const
{
  // Return a temporary inverse of the alignment
  // object. This means 'mis
   static AliAlignObjMatrix a;
   a = *this;

   TGeoHMatrix m;
   GetMatrix(m);
   a.SetMatrix(m.Inverse());

   return a;
}
