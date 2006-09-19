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
AliAlignObjAngles::AliAlignObjAngles() : AliAlignObj()
{
  // default constructor
  //
  fTranslation[0]=fTranslation[1]=fTranslation[2]=0.;
  fRotation[0]=fRotation[1]=fRotation[2]=0.;
}

//_____________________________________________________________________________
AliAlignObjAngles::AliAlignObjAngles(const char* volpath, UShort_t volUId, Double_t x, Double_t y, Double_t z, Double_t psi, Double_t theta, Double_t phi, Bool_t global) throw (const Char_t *) : AliAlignObj(volpath,volUId)
{
  // standard constructor with 3 translation + 3 rotation parameters
  // If the user explicitly sets the global variable to kFALSE then the
  // parameters are interpreted as giving the local transformation.
  // This requires to have a gGeoMenager active instance, otherwise the
  // constructor will fail (no object created)
  // 
  if(global){
    SetPars(x, y, z, psi, theta, phi);
  }else{
    if(!SetLocalPars(x,y,z,psi,theta,phi)) throw "Alignment object creation failed (TGeo instance needed)!\n";
  }
}

//_____________________________________________________________________________
AliAlignObjAngles::AliAlignObjAngles(const char* volpath, UShort_t volUId, TGeoMatrix& m, Bool_t global) throw (const Char_t *) : AliAlignObj(volpath,volUId)
{
  // standard constructor with TGeoMatrix
  // If the user explicitly sets the global variable to kFALSE then the
  // parameters are interpreted as giving the local transformation.
  // This requires to have a gGeoMenager active instance, otherwise the
  // constructor will fail (no object created)
  //

  if (!SetMatrix(m)) throw "Alignment object creation failed (can't extract roll-pitch-yall angles from the matrix)!\n";

  if (!global) {
    if (!SetLocalPars(fTranslation[0],fTranslation[1],fTranslation[2],fRotation[0],fRotation[1],fRotation[2])) throw "Alignment object creation failed (TGeo instance needed)!\n";
  }
}

//_____________________________________________________________________________
AliAlignObjAngles::AliAlignObjAngles(const AliAlignObj& theAlignObj) :
  AliAlignObj(theAlignObj)
{
  // copy constructor
  //
  Double_t tr[3];
  theAlignObj.GetTranslation(tr);
  SetTranslation(tr[0],tr[1],tr[2]);
  Double_t rot[3];
  if (theAlignObj.GetAngles(rot))
    SetRotation(rot[0],rot[1],rot[2]);
}

//_____________________________________________________________________________
AliAlignObjAngles &AliAlignObjAngles::operator =(const AliAlignObj& theAlignObj)
{
  // assignment operator
  //
  if(this==&theAlignObj) return *this;
  ((AliAlignObj *)this)->operator=(theAlignObj);

  Double_t tr[3];
  theAlignObj.GetTranslation(tr);
  SetTranslation(tr[0],tr[1],tr[2]);
  Double_t rot[3];
  if (theAlignObj.GetAngles(rot))
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
    fRotation[0] = fRotation[1] = fRotation[2] = 0.;
    return kTRUE;
  }
}

//_____________________________________________________________________________
void AliAlignObjAngles::GetMatrix(TGeoHMatrix& m) const
{
  m.SetTranslation(&fTranslation[0]);
  Double_t rot[9];
  AnglesToMatrix(fRotation,rot);
  m.SetRotation(rot);
}

//_____________________________________________________________________________
AliAlignObj& AliAlignObjAngles::Inverse() const
{
  // Return a temporary inverse of the alignment
  // object. This means 'mis
   static AliAlignObjAngles a;
   a = *this;

   TGeoHMatrix m;
   GetMatrix(m);
   a.SetMatrix(m.Inverse());

   return a;
}
