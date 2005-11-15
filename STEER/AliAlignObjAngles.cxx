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
//   the concrete representation of alignment object class
//   AliAlignObjAngles derived from the base class AliAlignObj
//-----------------------------------------------------------------

#include "AliAlignObj.h"
#include "AliAlignObjAngles.h"
//#include "AliLog.h"

ClassImp(AliAlignObjAngles)

//_____________________________________________________________________________
AliAlignObjAngles::AliAlignObjAngles() //: AliAlignObj()
{
  // default constructor
  //
  fTranslation[0]=fTranslation[1]=fTranslation[2]=0.;
  fRotation[0]=fRotation[1]=fRotation[2]=0.;
}

//_____________________________________________________________________________
AliAlignObjAngles::AliAlignObjAngles(const char* volpath, UShort_t voluid, Double_t x, Double_t y, Double_t z, Double_t psi, Double_t theta, Double_t phi)
{
  // standard constructor with 3 translation + 3 rotation parameters
  //
  fVolPath=volpath;
  fVolUID=voluid;
  fTranslation[0]=x; fTranslation[1]=y; fTranslation[2]=z;
  fRotation[0]=psi; fRotation[1]=theta; fRotation[2]=phi;
}

//_____________________________________________________________________________
AliAlignObjAngles::AliAlignObjAngles(const char* volpath, ELayerID detId, Int_t volId, Double_t x, Double_t y, Double_t z, Double_t psi, Double_t theta, Double_t phi)
{
  // standard constructor with 3 translation + 3 rotation parameters
  //
  fVolPath=volpath;
  SetVolUID(detId,volId);
  fTranslation[0]=x; fTranslation[1]=y; fTranslation[2]=z;
  fRotation[0]=psi; fRotation[1]=theta; fRotation[2]=phi;
}

//_____________________________________________________________________________
AliAlignObjAngles::AliAlignObjAngles(const char* volpath, UShort_t voluid, TGeoMatrix& m)
{
  // standard constructor with TGeoMatrix
  //
  fVolPath=volpath;
  fVolUID=voluid;
  SetTranslation(m);
  SetRotation(m);
}

//_____________________________________________________________________________
AliAlignObjAngles::AliAlignObjAngles(const AliAlignObjAngles& theAlignObj) :
  AliAlignObj(theAlignObj)
{
  // copy constructor
  //
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
AliAlignObjAngles::~AliAlignObjAngles()
{
  // default destructor
  //
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

