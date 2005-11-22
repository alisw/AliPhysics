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

// $Id$
//
// Class AliMUONGeometryDetElement
// --------------------------------------
// The class defines the detection element.
//
// Author: Ivana Hrivnacova, IPN Orsay

#include <TGeoMatrix.h>
#include <Riostream.h>

#include "AliLog.h"

#include "AliMUONGeometryDetElement.h"

ClassImp(AliMUONGeometryDetElement)

//______________________________________________________________________________
AliMUONGeometryDetElement::AliMUONGeometryDetElement(
                                        Int_t detElemId,
                                        const TString& alignedVolume, 
				        const TGeoCombiTrans& relTransform)
 : TObject(),
   fAlignedVolume(alignedVolume),
   fLocalTransformation(0),
   fGlobalTransformation(0)
{ 
/// Standard constructor

  SetUniqueID(detElemId);
  fLocalTransformation = new TGeoCombiTrans(relTransform);
}

//______________________________________________________________________________
AliMUONGeometryDetElement::AliMUONGeometryDetElement()
 : TObject(),
   fAlignedVolume(),
   fLocalTransformation(0),
   fGlobalTransformation(0)
{
/// Default constructor
}

//______________________________________________________________________________
AliMUONGeometryDetElement::AliMUONGeometryDetElement(
                                   const AliMUONGeometryDetElement& rhs)
  : TObject(rhs)
{
/// Protected copy constructor

  AliFatal("Copy constructor is not implemented.");
}

//______________________________________________________________________________
AliMUONGeometryDetElement::~AliMUONGeometryDetElement() 
{
/// Destructor

  delete fLocalTransformation;
  delete fGlobalTransformation;
}

//______________________________________________________________________________
AliMUONGeometryDetElement& 
AliMUONGeometryDetElement::operator = (const AliMUONGeometryDetElement& rhs) 
{
/// Protected assignement operator

  // check assignement to self
  if (this == &rhs) return *this;

  AliFatal("Assignment operator is not implemented.");
    
  return *this;  
}

//
// private methods
//

//______________________________________________________________________________
void  AliMUONGeometryDetElement::PrintTransform(
                                            const TGeoCombiTrans* transform) const
{
/// Prints the detection element transformation

  cout << "DetElemId: " << GetUniqueID();
  cout << "  name: " << fAlignedVolume << endl;

  const double* translation = transform->GetTranslation();
  cout << "   translation: "
#if defined (__DECCXX)
       << translation[0] << ", " 
       << translation[1] << ", "
       << translation[2] << endl;
#else
       << std::fixed
       << std::setw(7) << std::setprecision(4) << translation[0] << ", " 
       << std::setw(7) << std::setprecision(4) << translation[1] << ", "
       << std::setw(7) << std::setprecision(4) << translation[2] << endl;
#endif
	 
  const double* rotation = transform->GetRotationMatrix();
  cout << "   rotation matrix:  "
#if defined (__DECCXX)
       << rotation[0] << ", " << rotation[1] << ", " << rotation[2] << endl
       << "                     "	    
       << rotation[3] << ", " << rotation[4] << ", " << rotation[5] << endl	    
       << "                     "	    
       << rotation[6] << ", " << rotation[7] << ", " << rotation[8] << endl;
#else
       << std::fixed
       << std::setw(7) << std::setprecision(4) 
       << rotation[0] << ", " << rotation[1] << ", " << rotation[2] << endl
       << "                     "	    
       << rotation[3] << ", " << rotation[4] << ", " << rotation[5] << endl	    
       << "                     "	    
       << rotation[6] << ", " << rotation[7] << ", " << rotation[8] << endl;
#endif
}     

//
// public methods
//

//______________________________________________________________________________
void  AliMUONGeometryDetElement::Global2Local(
                                  Float_t xg, Float_t yg, Float_t zg, 
                                  Float_t& xl, Float_t& yl, Float_t& zl) const
{
/// Transforms point from the global reference frame (ALIC)
/// to the local reference frame of the detection element specified
/// by detElemId.

  Double_t dxg = xg;
  Double_t dyg = yg;
  Double_t dzg = zg;

  Double_t dxl, dyl, dzl;
  Global2Local(dxg, dyg, dzg, dxl, dyl, dzl);  
  
  xl = dxl;
  yl = dyl;
  zl = dzl;
}
				  
//______________________________________________________________________________
void  AliMUONGeometryDetElement::Global2Local(
                                  Double_t xg, Double_t yg, Double_t zg, 
                                  Double_t& xl, Double_t& yl, Double_t& zl) const
{
/// Transforms point from the global reference frame (ALIC)
/// to the local reference frame of the detection element specified
/// by detElemId.

   // Check transformation
   if (!fGlobalTransformation) {
     AliError(Form("Global transformation for detection element %d not defined.",
                    GetUniqueID()));
     return;
   }  
   
   // Transform point 
   Double_t pg[3] = { xg, yg, zg };
   Double_t pl[3] = { 0., 0., 0. };
   fGlobalTransformation->MasterToLocal(pg, pl);
   
   // Return transformed point
   xl = pl[0];
   yl = pl[1];     
   zl = pl[2];  
}
				  
//______________________________________________________________________________
void  AliMUONGeometryDetElement::Local2Global(
                 Float_t xl, Float_t yl, Float_t zl, 
                 Float_t& xg, Float_t& yg, Float_t& zg) const
{
/// Transforms point from the local reference frame of the detection element 
/// specified by detElemId to the global reference frame (ALIC).

  Double_t dxl = xl;
  Double_t dyl = yl;
  Double_t dzl = zl;

  Double_t dxg, dyg, dzg;
  Local2Global(dxl, dyl, dzl, dxg, dyg, dzg);  
  
  xg = dxg;
  yg = dyg;
  zg = dzg;
}

//______________________________________________________________________________
void  AliMUONGeometryDetElement::Local2Global(
                 Double_t xl, Double_t yl, Double_t zl, 
                 Double_t& xg, Double_t& yg, Double_t& zg) const
{
/// Transforms point from the local reference frame of the detection element 
/// specified by detElemId to the global reference frame (ALIC).

   // Check transformation
   if (!fGlobalTransformation) {
     AliError(Form("Global transformation for detection element %d not defined.",
                    GetUniqueID()));
     return;
   }  
   
   // Transform point 
   Double_t pl[3] = { xl, yl, zl };
   Double_t pg[3] = { 0., 0., 0. };
   fGlobalTransformation->LocalToMaster(pl, pg);
   
   // Return transformed point
   xg = pg[0];
   yg = pg[1];     
   zg = pg[2];  
}

//______________________________________________________________________________
void AliMUONGeometryDetElement::SetGlobalTransformation(
                                                const TGeoCombiTrans& transform)
{ 
/// Sets global transformation;
/// gives warning if the global transformation is already defined.
 
  if (fGlobalTransformation) {
    delete fGlobalTransformation;
    AliWarning("Global transformation already defined was deleted.");
  }  

  fGlobalTransformation = new TGeoCombiTrans(transform);
}  
					      
//______________________________________________________________________________
void AliMUONGeometryDetElement::PrintLocalTransform() const
{
/// Prints detection element relative transformation
/// (the transformation wrt module frame)

  PrintTransform(fLocalTransformation);
}  

//______________________________________________________________________________
void AliMUONGeometryDetElement::PrintGlobalTransform() const
{
/// Prints detection element global transformation
/// (the transformation wrt global frame)

  PrintTransform(fGlobalTransformation);
}  
