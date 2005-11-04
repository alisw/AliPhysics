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
//   1) the abstract class AliAlignObj
//   2) two derived concrete representation of alignment object class:
//      - AliAlignObjAngles
//      - AliAlignObjMatrix
//-----------------------------------------------------------------

#include "AliAlignObj.h"
//#include "AliLog.h"

ClassImp(AliAlignObj)

//_____________________________________________________________________________
AliAlignObj::AliAlignObj():
  fVolUID(0)
{
  // dummy constructor
}

//_____________________________________________________________________________
AliAlignObj::AliAlignObj(const AliAlignObj& theAlignObj) :
  TObject(theAlignObj)
{
  //copy constructor
  fVolPath = theAlignObj.GetVolPath();
  fVolUID = theAlignObj.GetVolUID();
}

//_____________________________________________________________________________
AliAlignObj &AliAlignObj::operator =(const AliAlignObj& theAlignObj)
{
  // assignment operator
  if(this==&theAlignObj) return *this;
  fVolPath = theAlignObj.GetVolPath();
  fVolUID = theAlignObj.GetVolUID();
  return *this;
}

//_____________________________________________________________________________
AliAlignObj::~AliAlignObj()
{
  // dummy destructor
}

//_____________________________________________________________________________
void AliAlignObj::AnglesToMatrix(const Double_t *angles, Double_t *rot) const
{
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
  if(rot[0]<1e-7 || rot[8]<1e-7) return kFALSE;
  Double_t raddeg = TMath::RadToDeg();
  angles[0]=raddeg*TMath::ATan2(-rot[5],rot[8]);
  angles[1]=raddeg*TMath::ASin(rot[2]);
  angles[2]=raddeg*TMath::ATan2(-rot[1],rot[0]);
  return kTRUE;
}

//_____________________________________________________________________________
void AliAlignObj::Print(Option_t *) const
{
  // Print the contents of the
  // alignment object in angles and
  // matrix representations
  Double_t tr[3];
  GetTranslation(tr);
  Double_t angles[3];
  GetAngles(angles);
  TGeoHMatrix m;
  GetMatrix(m);
  const Double_t *rot = m.GetRotationMatrix();
  printf("Volume=%s ID=%u\n", GetVolPath(),GetVolUID());
  printf("%12.6f%12.6f%12.6f    Tx = %12.6f    Psi   = %12.6f\n", rot[0], rot[1], rot[2], tr[0], angles[0]);
  printf("%12.6f%12.6f%12.6f    Ty = %12.6f    Theta = %12.6f\n", rot[3], rot[4], rot[5], tr[1], angles[1]);
  printf("%12.6f%12.6f%12.6f    Tz = %12.6f    Phi   = %12.6f\n", rot[6], rot[7], rot[8], tr[2], angles[2]);

}


//=============================================================================

ClassImp(AliAlignObjAngles)

//_____________________________________________________________________________
AliAlignObjAngles::AliAlignObjAngles() //: AliAlignObj()
{
  // default constructor
  fTranslation[0]=fTranslation[1]=fTranslation[2]=0.;
  fRotation[0]=fRotation[1]=fRotation[2]=0.;
}

//_____________________________________________________________________________
AliAlignObjAngles::AliAlignObjAngles(const AliAlignObjAngles& theAlignObj) :
  AliAlignObj(theAlignObj)
{
  // copy constructor
  Double_t tr[3];
  theAlignObj.GetTranslation(tr);
  SetTranslation(tr[0],tr[1],tr[2]);
  Double_t rot[3];
  theAlignObj.GetAngles(rot);
  SetRotation(rot[0],rot[1],rot[2]);
}

//_____________________________________________________________________________
AliAlignObjAngles &AliAlignObjAngles::operator =(const AliAlignObjAngles& theAlignObj)
{
  // assignment operator
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
AliAlignObjAngles::~AliAlignObjAngles()
{
  // default destructor
}

//_____________________________________________________________________________
void AliAlignObjAngles::SetTranslation(const TGeoMatrix& m)
{
  if(m.IsTranslation()){
    const Double_t* tr = m.GetTranslation();
    fTranslation[0]=tr[0];  fTranslation[1]=tr[1]; fTranslation[2]=tr[2];
  }else{
//     AliWarning("Argument matrix is not a translation! Setting zero-translation.");
    fTranslation[0] = fTranslation[1] = fTranslation[2] = 0.;
  }
}

//_____________________________________________________________________________
Bool_t AliAlignObjAngles::SetRotation(const TGeoMatrix& m)
{
  if(m.IsRotation()){
    const Double_t* rot = m.GetRotationMatrix();
    return MatrixToAngles(rot,fRotation);
  }else{
//     AliWarning("Argument matrix is not a rotation! Setting yaw-pitch-roll to zero.");
    fRotation[0] = fRotation[1] = fRotation[2] = 0.;
    return kTRUE;
  }
}

//_____________________________________________________________________________
void AliAlignObjAngles::SetMatrix(const TGeoMatrix& m)
{
  SetTranslation(m);
  SetRotation(m);
}

//_____________________________________________________________________________
void AliAlignObjAngles::GetPars(Double_t tr[], Double_t angles[]) const
{
  GetTranslation(tr);
  GetAngles(angles);
}

//_____________________________________________________________________________
void AliAlignObjAngles::GetMatrix(TGeoHMatrix& m) const
{
  m.SetTranslation(&fTranslation[0]);
  Double_t rot[9];
  AnglesToMatrix(fRotation,rot);
  m.SetRotation(rot);
}

//=============================================================================

ClassImp(AliAlignObjMatrix)

//_____________________________________________________________________________
AliAlignObjMatrix::AliAlignObjMatrix() : AliAlignObj()
{
  // Default constructor
}

AliAlignObjMatrix::AliAlignObjMatrix(const AliAlignObjMatrix& theAlignObj) :
  AliAlignObj(theAlignObj)
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

AliAlignObjMatrix &AliAlignObjMatrix::operator =(const AliAlignObjMatrix& theAlignObj)
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

