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
// -------------------------------------
// Class AliMUONGeometryModuleTransformer
// -------------------------------------
// Class for definition of the detector module transformations
// Author: Ivana Hrivnacova, IPN Orsay

#include "AliMUONGeometryModuleTransformer.h"
#include "AliMUONGeometryDetElement.h"	
#include "AliMUONGeometryStore.h"	

#include "AliLog.h"	

#include <TVirtualMC.h>
#include <TGeoMatrix.h>
#include <TObjArray.h>
#include <TArrayI.h>
#include <Riostream.h>

/// \cond CLASSIMP
ClassImp(AliMUONGeometryModuleTransformer)
/// \endcond

const TString AliMUONGeometryModuleTransformer::fgkModuleNamePrefix = "GM";

//______________________________________________________________________________
AliMUONGeometryModuleTransformer::AliMUONGeometryModuleTransformer(Int_t moduleId)
 : TObject(),
   fModuleId(moduleId),
   fModuleName(),
   fVolumePath(),
   fTransformation(0),
   fDetElements(0)
{
/// Standard constructor

  // Chamber transformation
  fTransformation = new TGeoHMatrix("");

  // Det elements transformation stores
  fDetElements = new AliMUONGeometryStore(true);
  
  // Compose module name
  fModuleName = fgkModuleNamePrefix;
  fModuleName += moduleId;
}


//______________________________________________________________________________
AliMUONGeometryModuleTransformer::AliMUONGeometryModuleTransformer()
 : TObject(),
   fModuleId(0),
   fModuleName(),
   fVolumePath(),
   fTransformation(0),
   fDetElements(0)
{
/// Default constructor
}


//______________________________________________________________________________
AliMUONGeometryModuleTransformer::~AliMUONGeometryModuleTransformer() 
{
/// Destructor

  delete fTransformation;
  delete fDetElements;
}

//
// public methods
//

//______________________________________________________________________________
void  AliMUONGeometryModuleTransformer::Global2Local(Int_t detElemId,
                  Float_t xg, Float_t yg, Float_t zg, 
                  Float_t& xl, Float_t& yl, Float_t& zl) const
{
/// Transform point from the global reference frame (ALIC)
/// to the local reference frame of the detection element specified
/// by detElemId.

  // Get detection element
  AliMUONGeometryDetElement* detElement = GetDetElement(detElemId);
  if (!detElement) return;
   
  // Transform point
  detElement->Global2Local(xg, yg, zg, xl, yl, zl);
}
				  
//______________________________________________________________________________
void  AliMUONGeometryModuleTransformer::Global2Local(Int_t detElemId,
                  Double_t xg, Double_t yg, Double_t zg, 
                  Double_t& xl, Double_t& yl, Double_t& zl) const
{
/// Transform point from the global reference frame (ALIC)
/// to the local reference frame of the detection element specified
/// by detElemId.

   // Get detection element
   AliMUONGeometryDetElement* detElement = GetDetElement(detElemId);
   if (!detElement) return;
   
   // Transform point
   detElement->Global2Local(xg, yg, zg, xl, yl, zl);
}
				  
//______________________________________________________________________________
void  AliMUONGeometryModuleTransformer::Local2Global(Int_t detElemId,
                 Float_t xl, Float_t yl, Float_t zl, 
                 Float_t& xg, Float_t& yg, Float_t& zg) const
{
/// Transform point from the local reference frame of the detection element 
/// specified by detElemId to the global reference frame (ALIC).

  // Get detection element
  AliMUONGeometryDetElement* detElement = GetDetElement(detElemId);
  if (!detElement) return;
   
   // Transform point
  detElement->Local2Global(xl, yl, zl, xg, yg, zg);  
}

//______________________________________________________________________________
void  AliMUONGeometryModuleTransformer::Local2Global(Int_t detElemId,
                 Double_t xl, Double_t yl, Double_t zl, 
                 Double_t& xg, Double_t& yg, Double_t& zg) const
{
/// Transform point from the local reference frame of the detection element 
/// specified by detElemId to the global reference frame (ALIC).

   // Get detection element
   AliMUONGeometryDetElement* detElement = GetDetElement(detElemId);
   if (!detElement) return;
   
   // Transform point
   detElement->Local2Global(xl, yl, zl, xg, yg, zg); 
}

//______________________________________________________________________________
void  AliMUONGeometryModuleTransformer::SetTransformation(
                                           const TGeoHMatrix& transform)
{
/// Set the module position wrt world.

  *fTransformation = transform;
}  

//______________________________________________________________________________
TString AliMUONGeometryModuleTransformer::GetVolumeName() const
{ 
/// Extract volume name from the path
  
  std::string volPath = fVolumePath.Data();
  std::string::size_type first = volPath.rfind('/')+1;
  std::string::size_type last = volPath.rfind('_');
  
  return volPath.substr(first, last-first );
}

//______________________________________________________________________________
TString AliMUONGeometryModuleTransformer::GetMotherVolumeName() const
{ 
/// Extract mother volume name from the path
  
  std::string volPath = fVolumePath.Data();
  std::string::size_type first = volPath.rfind('/');
  volPath = volPath.substr(0, first);

  std::string::size_type next = volPath.rfind('/')+1;
  std::string::size_type last = volPath.rfind('_');
  
  return volPath.substr(next, last-next );
}

//______________________________________________________________________________
AliMUONGeometryDetElement*
AliMUONGeometryModuleTransformer::GetDetElement(Int_t detElemId, Bool_t warn) const
{
/// Return the detection element specified by detElemId.
/// Give error if detection element is not defined and warn is true.

   // Get detection element
   AliMUONGeometryDetElement* detElement
     = (AliMUONGeometryDetElement*) fDetElements->Get(detElemId, warn);

   if (!detElement) {
     if (warn)
       AliErrorStream() 
         << "Detection element " << detElemId
         << " not found in module " << fModuleId << endl;
     return 0;
   }  

   return detElement;
}
