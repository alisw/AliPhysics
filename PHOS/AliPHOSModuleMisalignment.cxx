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

/* $Id: $ */

#include <stdexcept>
#include <iostream>

#include <TGeoMatrix.h>
#include <TObjArray.h>
#include <TString.h>
#include <TError.h>
#include <TMath.h>

#include "AliSurveyPoint.h"

#include "AliPHOSModuleMisalignment.h"
#include "AliPHOSGeometry.h"

ClassImp(AliPHOSModuleMisalignment)

namespace {

  /*
    Tiny vector/matrix utility stuff. Operates on doubles directly.
    Instead of points and vectors I use arrays of doubles with size 3.
    To make this explicit - I use references to arrays.
  */
  
  //___________________________________________________________________
  void Vector(const Double_t (&p1)[3], const Double_t (&p2)[3], Double_t (&v)[3])
  {
    for(UInt_t i = 0; i < 3; ++i)
      v[i] = p2[i] - p1[i];
  }
  
#if 0
  //___________________________________________________________________
  void MultVector(Double_t (&v)[3], Double_t m)
  {
    v[0] *= m;
    v[1] *= m;
    v[2] *= m;
  }
#endif
  
  /*
    Using points name0, name1, name2 find two orthogonal vectors.
  */
  //___________________________________________________________________
  void FindVectors(const Double_t (&pts)[3][3], Double_t (&v)[3][3])
  {
    Vector(pts[0], pts[2], v[0]);
    //v[1] will be cross-product (v[2] x v[0]).
    Vector(pts[0], pts[1], v[2]);
  }
  
  //___________________________________________________________________
  Double_t Length(const Double_t (&v)[3])
  {
    return TMath::Sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
  }
  
#if 0
  //___________________________________________________________________
  Double_t Distance(const Double_t (&p1)[3], const Double_t (&p2)[3])
  {
    return TMath::Sqrt((p2[0] - p1[0]) * (p2[0] - p1[0]) +
		       (p2[1] - p1[1]) * (p2[1] - p1[1]) +
		       (p2[2] - p1[2]) * (p2[2] - p1[2]));
  }
  
  //______________________________________________________________________________
#endif
  void CrossProduct(const Double_t (&v1)[3], const Double_t (&v2)[3], Double_t (&v3)[3])
  {
    v3[0] = v1[1] * v2[2] - v2[1] * v1[2];
    v3[1] = v1[2] * v2[0] - v2[2] * v1[0];
    v3[2] = v1[0] * v2[1] - v2[0] * v1[1];
  }
  
  //___________________________________________________________________
  void Normalize(Double_t (&v)[3])
  {
    const Double_t len = Length(v);
    if(len < 1E-10)//Threshold?
      throw std::runtime_error("Zero len vector");
    v[0] /= len;
    v[1] /= len;
    v[2] /= len;
  }
  
  //______________________________________________________________________________
  void Normalize(Double_t (&v)[3][3])
  {
    for(UInt_t i = 0; i < 3; ++i)
      Normalize(v[i]);
  }
  
  
  //___________________________________________________________________
  void FindRotation(const Double_t (&u)[3][3], const Double_t (&v)[3][3], Double_t (&r)[9])
  {
    //I have orthogonal vectors and very nice rotation matrix.
    //V = R * U, R = V * U ^ t
    r[0] = v[0][0] * u[0][0] + v[1][0] * u[1][0] + v[2][0] * u[2][0];
    r[1] = v[0][0] * u[0][1] + v[1][0] * u[1][1] + v[2][0] * u[2][1];
    r[2] = v[0][0] * u[0][2] + v[1][0] * u[1][2] + v[2][0] * u[2][2];
    
    r[3] = v[0][1] * u[0][0] + v[1][1] * u[1][0] + v[2][1] * u[2][0];
    r[4] = v[0][1] * u[0][1] + v[1][1] * u[1][1] + v[2][1] * u[2][1];
    r[5] = v[0][1] * u[0][2] + v[1][1] * u[1][2] + v[2][1] * u[2][2];
    
    r[6] = v[0][2] * u[0][0] + v[1][2] * u[1][0] + v[2][2] * u[2][0];
    r[7] = v[0][2] * u[0][1] + v[1][2] * u[1][1] + v[2][2] * u[2][1];
    r[8] = v[0][2] * u[0][2] + v[1][2] * u[1][2] + v[2][2] * u[2][2];
  }
  
  //___________________________________________________________________
  void Rotate(const Double_t (&r)[9], const Double_t (&u)[3], Double_t (&v)[3])
  {
    v[0] = r[0] * u[0] + r[1] * u[1] + r[2] * u[2];
    v[1] = r[3] * u[0] + r[4] * u[1] + r[5] * u[2];
    v[2] = r[6] * u[0] + r[7] * u[1] + r[8] * u[2];
  }
  
  //___________________________________________________________________
  void Rotate(const Double_t (&r)[9], const Double_t (&u)[3][3], Double_t (&v)[3][3])
  {
    for(UInt_t i = 0; i < 3; ++i)
      Rotate(r, u[i], v[i]);
  }
  
  /*
    PrintVector, PrintMatrix, Translate are used in "debug" mode only.
  */
  //___________________________________________________________________
  void PrintVector(const Double_t (&v)[3])
  {
    std::cout<<v[0]<<' '<<v[1]<<' '<<v[2]<<std::endl;
  }
  
  //___________________________________________________________________
  void PrintMatrix(const Double_t (&u)[3][3])
  {
    for(UInt_t i = 0; i < 3; ++i)
      PrintVector(u[i]);
  }
  
  //___________________________________________________________________
  void Translate(const Double_t (&t)[3], const Double_t (&u)[3], Double_t (&v)[3])
  {
    for(UInt_t i = 0; i < 3; ++i)
      v[i] = u[i] + t[i];
  }
  
  //___________________________________________________________________
  void Translate(const Double_t (&t)[3], const Double_t (&u)[3][3], Double_t (&v)[3][3])
  {
    for(UInt_t i = 0; i < 3; ++i)
      Translate(t, u[i], v[i]);
  }
  
}

//______________________________________________________________________________
AliPHOSModuleMisalignment::
AliPHOSModuleMisalignment(const AliPHOSGeometry & geom, Bool_t debug)
                            : fDebug(debug),
                              fAngles(),
                              fCenters(),
                              fModule(),
                              fU(),
                              fV()
{
  //Ctor.
  //Extract ideal module transformations from AliPHOSGeometry.

  //Angles.
  for (UInt_t module = 0; module < kModules; ++module)
    for (UInt_t axis = 0; axis < 3; ++axis)
      for (UInt_t angle = 0; angle < 2; ++angle)
        fAngles[module][axis][angle] = geom.GetModuleAngle(module, axis, angle);
  //Translations.
  for (UInt_t module = 0; module < kModules; ++module)
    for (UInt_t axis = 0; axis < 3; ++axis)
      fCenters[module][axis] = geom.GetModuleCenter(module, axis);
  //Points, will be rotated/translated using module angle/center.
  fModule[0][0] = -geom.GetNPhi() / 2. * geom.GetCellStep() + geom.GetCellStep() / 2.;
  fModule[0][1] = -geom.GetNZ()   / 2. * geom.GetCellStep() + geom.GetCellStep() / 2.;
  fModule[0][2] = -22.61;//This number is hardcoded, AliPHOSGeometry does not have it,
                         //only 460. but this is result of transformation applied already.
  fModule[1][0] = fModule[0][0];
  fModule[1][1] = -fModule[0][1] - geom.GetCellStep();
  fModule[1][2] = -22.61;

  fModule[2][0] = -fModule[0][0] - 7 * geom.GetCellStep();
  fModule[2][1] = fModule[0][1];
  fModule[2][2] = -22.61;
}

//______________________________________________________________________________
AliPHOSModuleMisalignment::~AliPHOSModuleMisalignment()
{
}

//______________________________________________________________________________
void AliPHOSModuleMisalignment::
DeltaTransformation(UInt_t module, const TObjArray * points, 
                    const TString & name0, const TString & name1, 
                    const TString & name2, TGeoHMatrix * delta)
{
  //Find delta transformation to misalign module. Global transformation.
  const AliSurveyPoint * pt0 = static_cast<AliSurveyPoint *>(points->FindObject(name0));
  const AliSurveyPoint * pt1 = static_cast<AliSurveyPoint *>(points->FindObject(name1));
  const AliSurveyPoint * pt2 = static_cast<AliSurveyPoint *>(points->FindObject(name2));

  if (!pt0 || !pt1 || !pt2) {
    Warning("AliPHOSModuleData::DeltaTransformation",
            "One of points not found in TObjArray");
    return;
  }

  //Transform fModule using angle and translation for module number "module".
  //fU.
  FindIdealModule(module);
  //Extract coordinates from survey.
  //fV.
  FindRealModule(pt0, pt1, pt2);
  //Find delta, using ideal module (fU) and survey data (fV).
  FindDelta(delta);
}

//______________________________________________________________________________
void AliPHOSModuleMisalignment::FindIdealModule(UInt_t module)
{
  //Ideal module described by fU.
  TGeoHMatrix matrix;
  //Rotation.
  const TGeoRotation r("", 
                       fAngles[module][0][0], fAngles[module][0][1], 
                       fAngles[module][1][0], fAngles[module][1][1], 
                       fAngles[module][2][0], fAngles[module][2][1]);
  matrix.SetRotation(r.GetRotationMatrix());
  //Translation.
  matrix.SetDx(fCenters[module][0]);
  matrix.SetDy(fCenters[module][1]);
  matrix.SetDz(fCenters[module][2]);
  //Find ideal module's points.
  matrix.LocalToMaster(fModule[0], fU[0]);
  matrix.LocalToMaster(fModule[1], fU[1]);
  matrix.LocalToMaster(fModule[2], fU[2]);
}

//______________________________________________________________________________
void AliPHOSModuleMisalignment::FindRealModule(const AliSurveyPoint * pt0, 
                                               const AliSurveyPoint * pt1,
                                               const AliSurveyPoint * pt2)
{
  //Real module - fV.
  //Survey is in millimeters.
  //AliPHOSGeometry is in centimeters.
  const Double_t scale = 0.1;

  fV[0][0] = pt0->GetX() * scale;
  fV[0][1] = pt0->GetY() * scale;
  fV[0][2] = pt0->GetZ() * scale;

  fV[1][0] = pt1->GetX() * scale;
  fV[1][1] = pt1->GetY() * scale;
  fV[1][2] = pt1->GetZ() * scale;

  fV[2][0] = pt2->GetX() * scale;
  fV[2][1] = pt2->GetY() * scale;
  fV[2][2] = pt2->GetZ() * scale;
}

//______________________________________________________________________________
void AliPHOSModuleMisalignment::FindDelta(TGeoHMatrix * delta)const
{
  //Find rotation and translation wich can
  //convert fU into fV (ideal module points into points from survey).
  Double_t u[3][3] = {};
  FindVectors(fU, u);		      
  //Find cross product u2 x u0 and save it in u[2].
  CrossProduct(u[2], u[0], u[1]);
  /*
  const Double_t lenXideal = Length(u[0]);
  const Double_t lenZideal = Length(u[2]);
  */
  Double_t v[3][3] = {};
  FindVectors(fV, v);		      
  //Find cross product (v2 x v0) and save it in v[2].
  CrossProduct(v[2], v[0], v[1]);
  /*
  const Double_t lenXreal = Length(v[0]);
  const Double_t lenZreal = Length(v[2]);
  */
  //Now, find rotation matrix.
  //1. Normalize vectors in u and v.
  try {
    Normalize(u);
    Normalize(v);
  } catch (const std::exception & e) {
    //One of lengths is zero (in principle, impossible, just to be neat).
    Error("AliPHOSModuleMisalignment::FindDelta", 
	  "\tone of vectors from ideal or real\n\tpoints have zero size\n"
	  "\tzero misalignment will be created");
    return;
  }

  //2. Rotation matrix.
  Double_t r[9] = {};
  FindRotation(u, v, r);
  delta->SetRotation(r);
  
#if 1

  //3. Now, rotate fU and look, what translation I have to add.
  Double_t idealRotated[3] = {};
  Rotate(r, fU[0], idealRotated);

  delta->SetDx(fV[0][0] - idealRotated[0]);
  delta->SetDy(fV[0][1] - idealRotated[1]);
  delta->SetDz(fV[0][2] - idealRotated[2]);

  if (fDebug) {
    const Double_t shifts[3] = 
    {fV[0][0] - idealRotated[0], 
     fV[0][1] - idealRotated[1], 
     fV[0][2] - idealRotated[2]};

    Double_t test1[3][3] = {};
    Rotate(r, fU, test1);
    Double_t test2[3][3] = {};
    Translate(shifts, test1, test2);
    std::cout<<"ideal:\n";
    PrintMatrix(fU);
    std::cout<<"misaligned:\n";
    PrintMatrix(test2);
    std::cout<<"real:\n";
    PrintMatrix(fV);
  }

#endif

#if 0
  //3. Now, rotate fU and look, what translation I have to add.
  Double_t idealRotated[3][3] = {};
  Rotate(r, fU, idealRotated);
  //Because of measurement errors, distances 
  //between points has errors. I can try to split
  //this difference (and make "final errors" smaller).
  Double_t zShift[3] = {};
  Vector(fV[0], fV[1], zShift);
  Normalize(zShift);
  
  Double_t xShift[3] = {};
  Vector(fV[0], fV[2], xShift);
  Normalize(xShift);

  MultVector(zShift, 0.5 * (lenZreal - lenZideal));
  MultVector(xShift, 0.5 * (lenXreal - lenXideal));

  Double_t pt1[3] = {};
  Translate(zShift, fV[0], pt1);
  Double_t pt2[3] = {};
  Translate(xShift, pt1, pt2);

  Double_t shifts[] = {pt2[0] - idealRotated[0][0],
                       pt2[1] - idealRotated[0][1],
                       pt2[2] - idealRotated[0][2]};

  delta->SetDx(shifts[0]);
  delta->SetDy(shifts[1]);
  delta->SetDz(shifts[2]);

  if (fDebug) {
    Double_t idealTr[3][3] = {};
    Translate(shifts, idealRotated, idealTr);

    std::cout<<"misaligned:\n";
    PrintMatrix(idealTr);
    std::cout<<"ideal1 "<<Distance(idealTr[0], idealTr[1])<<std::endl;
    std::cout<<"ideal2 "<<Distance(idealTr[0], idealTr[2])<<std::endl;
    std::cout<<"real:\n";
    PrintMatrix(fV);
    std::cout<<"real1 "<<Distance(fV[0], fV[1])<<std::endl;
    std::cout<<"real2 "<<Distance(fV[0], fV[2])<<std::endl;
  }
#endif
}
