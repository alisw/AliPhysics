// $Id$
//
// Class AliMUONGeometryTransformStore
// ------------------------------------ 
// The class provides the array of transformations of the detection elements
// from a defined reference frame to the detection element frame.
// (See more in the header file.) 
//
// Author: Ivana Hrivnacova, IPN Orsay

#include <Riostream.h>
#include <TGeoMatrix.h>
#include <TObjString.h>

#include "AliMUONGeometryTransformStore.h"
#include "AliLog.h"

ClassImp(AliMUONGeometryTransformStore)

const Int_t AliMUONGeometryTransformStore::fgkHemisphere = 50; 

//______________________________________________________________________________
AliMUONGeometryTransformStore::AliMUONGeometryTransformStore(
                                  Int_t firstDetElemId, Int_t  nofDetElems,
				  AliMUONGeometrySVMap* svMap)
 : TObject(),
   fFirstDetElemId(firstDetElemId),
   fNofDetElems(nofDetElems),
   fDETransforms(100),
   fSVMap(svMap)
{ 
// Standard constructor
  
  fDETransforms.SetOwner(true);
  for (Int_t i=0; i<nofDetElems; i++) fDETransforms[i] = 0;
}

//______________________________________________________________________________
AliMUONGeometryTransformStore::AliMUONGeometryTransformStore()
 : TObject(),
   fFirstDetElemId(0),
   fNofDetElems(0),
   fDETransforms(),
   fSVMap(0)
{
// Default constructor
}

//______________________________________________________________________________
AliMUONGeometryTransformStore::AliMUONGeometryTransformStore(
                                   const AliMUONGeometryTransformStore& rhs)
  : TObject(rhs)
{
  AliFatal("Copy constructor is not implemented.");
}

//______________________________________________________________________________
AliMUONGeometryTransformStore::~AliMUONGeometryTransformStore() {
//
  fDETransforms.Delete();
}

//______________________________________________________________________________
AliMUONGeometryTransformStore& 
AliMUONGeometryTransformStore::operator = (const AliMUONGeometryTransformStore& rhs) 
{
  // check assignement to self
  if (this == &rhs) return *this;

  AliFatal("Assignment operator is not implemented.");
    
  return *this;  
}

//
// private methods
//

//______________________________________________________________________________
Int_t AliMUONGeometryTransformStore::GetDetElementIndex(Int_t detElemId) const
{
// Returns the index of detector element specified by detElemId
// ---

  Int_t index = detElemId - fFirstDetElemId;
  if (index >= fgkHemisphere) 
    index -= fgkHemisphere - fNofDetElems/2;

  return index;
}  

//______________________________________________________________________________
Int_t AliMUONGeometryTransformStore::GetDetElementId(Int_t detElemIndex) const
{
// Returns the ID of detector element specified by index
// ---

  Int_t detElemId = detElemIndex;
  
  if (detElemIndex >=  fNofDetElems/2) 
    detElemId += fgkHemisphere - fNofDetElems/2; 

  detElemId += fFirstDetElemId;

  return detElemId;
}  


//
// public methods
//

//______________________________________________________________________________
void AliMUONGeometryTransformStore::Add(Int_t detElemId,
                                        const TString& alignedVolume, 
	 	                        const TGeoCombiTrans& transform)
{
// Add transformation specified by alignedVolume in the array
// if this alignedVolume is not yet present. 
// ---
 
  TGeoCombiTrans* newTransform = new TGeoCombiTrans(transform);
  newTransform->SetName(alignedVolume);
  
  //cout << ".. adding " << detElemId 
  //     << " index: " << GetDetElementIndex(detElemId) << endl;

  // Add to the map  
  if ( !Get(detElemId)) {
  
    newTransform->SetUniqueID(detElemId);
      // Set detector element id as unique id
 
    fDETransforms.AddAt(newTransform, GetDetElementIndex(detElemId));
  } 
  else 
    AliWarning(Form("The aligned volume %s is already present", 
            alignedVolume.Data()));  
}		      
    
//______________________________________________________________________________
void  AliMUONGeometryTransformStore::Print(Option_t* /*option*/) const
{
// Prints the detector elements transformations
// ---

  for (Int_t i=0; i<fDETransforms.GetEntriesFast(); i++) {
    
    cout << "DetElemId: " << GetDetElementId(i);

    TGeoCombiTrans* matrix = (TGeoCombiTrans*)fDETransforms.At(i);
    cout << "  name: " << matrix->GetName() << endl;

    const double* translation = matrix->GetTranslation();
    cout << "   translation: "
         << std::fixed
         << std::setw(7) << std::setprecision(4) << translation[0] << ", " 
         << std::setw(7) << std::setprecision(4) << translation[1] << ", "
         << std::setw(7) << std::setprecision(4) << translation[2] << endl;
	 
    const double* rotation = matrix->GetRotationMatrix();
    cout << "   rotation matrix:  "
         << std::fixed
         << std::setw(7) << std::setprecision(4) 
         << rotation[0] << ", " << rotation[1] << ", " << rotation[2] << endl
	 << "                     "	    
         << rotation[3] << ", " << rotation[4] << ", " << rotation[5] << endl	    
	 << "                     "	    
         << rotation[6] << ", " << rotation[7] << ", " << rotation[8] << endl;

  }
}     

//______________________________________________________________________________
const TGeoCombiTrans* 
AliMUONGeometryTransformStore::Get(Int_t detElemId) const
{
// Returns transformation for the specified detector element Id
// ---

  Int_t index = GetDetElementIndex(detElemId);
  
  if ( index >= 0 && index < fNofDetElems )
    return (const TGeoCombiTrans*)fDETransforms.At(index);
  else {
    AliWarning(Form("Index %d out of limits", index));
    return 0;  
  }  
}  

//______________________________________________________________________________
const TGeoCombiTrans* 
AliMUONGeometryTransformStore::FindBySensitiveVolume(const TString& sensVolume) const
{
// Finds TGeoCombiTrans for the detector element Id specified by aligned volume 
// ---

  Int_t detElemId = fSVMap->GetDetElemId(sensVolume);

  if (!detElemId) return 0; 
        // The specified sensitive volume is not in the map   
  
  return Get(detElemId);
}  

