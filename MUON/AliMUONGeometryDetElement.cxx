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
// --------------------------------------
// Class AliMUONGeometryDetElement
// --------------------------------------
// The class defines the detection element.
// Author: Ivana Hrivnacova, IPN Orsay

#include "AliMUONGeometryDetElement.h"

#include "AliLog.h"

#include <TGeoMatrix.h>
#include <Riostream.h>

#include <sstream>

/// \cond CLASSIMP
ClassImp(AliMUONGeometryDetElement)
/// \endcond

const TString AliMUONGeometryDetElement::fgkDENamePrefix = "DE";

//______________________________________________________________________________
AliMUONGeometryDetElement::AliMUONGeometryDetElement(
                                        Int_t detElemId,
                                        const TString& volumePath)
 : TObject(),
   fDEName(),
   fVolumePath(volumePath),
   fLocalTransformation(0),
   fGlobalTransformation(0)
{ 
/// Standard constructor

  SetUniqueID(detElemId);
  
  fDEName = fgkDENamePrefix;
  fDEName += detElemId;
}

//______________________________________________________________________________
AliMUONGeometryDetElement::AliMUONGeometryDetElement()
 : TObject(),
   fVolumePath(),
   fLocalTransformation(0),
   fGlobalTransformation(0)
{
/// Default constructor
}

//______________________________________________________________________________
AliMUONGeometryDetElement::~AliMUONGeometryDetElement() 
{
/// Destructor

  delete fLocalTransformation;
  delete fGlobalTransformation;
}

//
// private methods
//

//______________________________________________________________________________
void  AliMUONGeometryDetElement::PrintTransform(
                                            const TGeoHMatrix* transform) const
{
/// Print the detection element transformation

  cout << "DetElemId: " << GetUniqueID();
  cout << "  name: " << fVolumePath << endl;

  if ( !transform ) {
    cout << "    Transformation not defined." << endl;
    return;
  }  

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
/// Transform point from the global reference frame (ALIC)
/// to the local reference frame of this detection element.

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
/// Transform point from the global reference frame (ALIC)
/// to the local reference frame of this detection element

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
/// Transform point from the local reference frame of this detection element 
/// to the global reference frame (ALIC).

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
/// Transform point from the local reference frame of this detection element 
/// to the global reference frame (ALIC).

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
void AliMUONGeometryDetElement::SetLocalTransformation(
                                                const TGeoHMatrix& transform)
{ 
/// Set local transformation;
/// give warning if the global transformation is already defined.
 
  if (fLocalTransformation) {
    delete fLocalTransformation;
    AliWarning("Local transformation already defined was deleted.");
  }  

  fLocalTransformation = new TGeoHMatrix(transform);
}  
					      
//______________________________________________________________________________
void AliMUONGeometryDetElement::SetGlobalTransformation(
                                                const TGeoHMatrix& transform)
{ 
/// Set global transformation;
/// give warning if the global transformation is already defined.
 
  if (fGlobalTransformation) {
    delete fGlobalTransformation;
    AliWarning("Global transformation already defined was deleted.");
  }  

  fGlobalTransformation = new TGeoHMatrix(transform);
}  
					      
//______________________________________________________________________________
void AliMUONGeometryDetElement::PrintLocalTransform() const
{
/// Print detection element relative transformation
/// (the transformation wrt module frame)

  PrintTransform(fLocalTransformation);
}  

//______________________________________________________________________________
void AliMUONGeometryDetElement::PrintGlobalTransform() const
{
/// Print detection element global transformation
/// (the transformation wrt global frame)

  PrintTransform(fGlobalTransformation);
}  

//______________________________________________________________________________
TString AliMUONGeometryDetElement::GetVolumeName() const
{ 
/// Extract volume name from the path
  
  std::string volPath = fVolumePath.Data();
  std::string::size_type first = volPath.rfind('/')+1;
  std::string::size_type last = volPath.rfind('_');
  
  return volPath.substr(first, last-first );
}

//______________________________________________________________________________
Int_t AliMUONGeometryDetElement::GetVolumeCopyNo() const
{ 
/// Extract volume copyNo from the path
  
  string volPath = fVolumePath.Data();
  std::string::size_type first = volPath.rfind('_');
  std::string copyNoStr = volPath.substr(first+1, volPath.length());
  std::istringstream in(copyNoStr);
  Int_t copyNo;
  in >> copyNo;
  
  return copyNo;
}

