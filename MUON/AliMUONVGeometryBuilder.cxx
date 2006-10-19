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
// Class AliMUONVGeometryBuilder
// -----------------------------
// Abstract base class for geometry construction per geometry module(s).
// Author: Ivana Hrivnacova, IPN Orsay
// 23/01/2004

#include "AliMUONVGeometryBuilder.h"
#include "AliMUONGeometryModule.h"
#include "AliMUONGeometryDetElement.h"
#include "AliMUONGeometryEnvelopeStore.h"
#include "AliMUONGeometryEnvelope.h"
#include "AliMUONGeometryConstituent.h"
#include "AliMUONGeometryBuilder.h"
#include "AliMUONStringIntMap.h"

#include "AliMpDEManager.h"
#include "AliMpExMap.h"

#include "AliLog.h"

#include <Riostream.h>
#include <TObjArray.h>
#include <TSystem.h>
#include <TGeoMatrix.h>
#include <TVirtualMC.h>

/// \cond CLASSIMP
ClassImp(AliMUONVGeometryBuilder)
/// \endcond

//______________________________________________________________________________
AliMUONVGeometryBuilder::AliMUONVGeometryBuilder(
                            Int_t firstModuleId, 
                            Int_t nofModules)
 : TObject(),
   fGeometryModules(0),
   fReferenceFrame()
 {
/// Standard constructor

  // Create the module geometries array
  fGeometryModules = new TObjArray();

  for (Int_t i=0; i<nofModules; i++ )
    fGeometryModules->Add(new AliMUONGeometryModule(firstModuleId++));
}

//______________________________________________________________________________
AliMUONVGeometryBuilder::AliMUONVGeometryBuilder()
 : TObject(),
   fGeometryModules(0),
   fReferenceFrame()
{
/// Default constructor
}

//______________________________________________________________________________
AliMUONVGeometryBuilder::~AliMUONVGeometryBuilder() 
{
/// Destructor

  if (fGeometryModules) {
    fGeometryModules->Clear(); // Sets pointers to 0 since it is not the owner
    delete fGeometryModules;
  }
}

//
// private methods
//

//______________________________________________________________________________
TGeoHMatrix 
AliMUONVGeometryBuilder::ConvertTransform(const TGeoHMatrix& transform) const
{
/// Convert transformation into the reference frame

  if ( fReferenceFrame.IsIdentity() )
    return transform;
  else  {
    return AliMUONGeometryBuilder::Multiply( fReferenceFrame,
  				  	     transform,
    					     fReferenceFrame.Inverse() );  
  }			    
}

//______________________________________________________________________________
TGeoHMatrix 
AliMUONVGeometryBuilder::ConvertDETransform(const TGeoHMatrix& transform) const
{
/// Convert DE transformation into the reference frame

  if ( fReferenceFrame.IsIdentity() )
    return transform;
  else  {
    return AliMUONGeometryBuilder::Multiply( fReferenceFrame,
  				  	     transform );  
  }			    
}

//______________________________________________________________________________
TString  AliMUONVGeometryBuilder::ComposePath(const TString& volName,
                                              Int_t copyNo) const
{
/// Compose path from given volName and copyNo

  TString path = "/";
  path += volName;
  path += '_';
  path += copyNo;
  
  return path;
}  

//______________________________________________________________________________
void AliMUONVGeometryBuilder::MapSV(const TString& path0, 
                                    const TString& volName, Int_t detElemId) const
{
/// Update the path with all daughters volumes recursively
/// and map it to the detection element Id if it is a sensitive volume

  // Get module sensitive volumes map
  Int_t moduleId = AliMpDEManager::GetGeomModuleId(detElemId);
  AliMUONStringIntMap* svMap = GetSVMap(moduleId);     

  Int_t nofDaughters = gMC->NofVolDaughters(volName);
  if (nofDaughters == 0) {

    // Get the name of the last volume in the path
    Ssiz_t npos1 = path0.Last('/')+1; 
    Ssiz_t npos2 = path0.Last('_');
    TString volName(path0(npos1, npos2-npos1));  
    
    // Check if it is sensitive volume
    Int_t moduleId = AliMpDEManager::GetGeomModuleId(detElemId);
    AliMUONGeometryModule* geometry = GetGeometry(moduleId);
    if (  geometry->IsSensitiveVolume(volName) &&
        ! svMap->Get(path0) ) {
      //cout << ".. adding to the map  " 
      //     <<  path0 << "  "  << detElemId << endl;
    
      // Map the sensitive volume to detection element
      svMap->Add(path0, detElemId); 
    }  
    return; 
  }  

  for (Int_t i=0; i<nofDaughters; i++) {
    Int_t copyNo = gMC->VolDaughterCopyNo(volName, i);
    TString newName =  gMC->VolDaughterName(volName, i);
            
    TString path = path0;
    path += ComposePath(newName, copyNo);

    MapSV(path, newName, detElemId);
  }
}     

//
// protected methods
//

//______________________________________________________________________________
AliMUONGeometryModule*  
AliMUONVGeometryBuilder::GetGeometry(Int_t moduleId) const
{
/// Return the module geometry specified by moduleId

  for (Int_t i=0; i<fGeometryModules->GetEntriesFast(); i++) {

    AliMUONGeometryModule* geometry 
      = (AliMUONGeometryModule*)fGeometryModules->At(i);

    if ( geometry->GetModuleId() == moduleId) return geometry;
  }   
  
  return 0;
}  

//______________________________________________________________________________
AliMUONGeometryEnvelopeStore*  
AliMUONVGeometryBuilder::GetEnvelopes(Int_t moduleId) const
{
/// Return the envelope store of the module geometry specified by moduleId

  AliMUONGeometryModule* geometry = GetGeometry(moduleId);
  
  if (!geometry) {
    AliFatal(Form("Module geometry %d is not defined", moduleId)); 
    return 0;
  }
  
  return geometry->GetEnvelopeStore();
}  

//______________________________________________________________________________
AliMUONStringIntMap*  
AliMUONVGeometryBuilder::GetSVMap(Int_t moduleId) const
{
/// Return the transformation store of the module geometry specified by moduleId

  AliMUONGeometryModule* geometry = GetGeometry(moduleId);
  
  if (!geometry) {
    AliFatal(Form("Geometry %d is not defined", moduleId)); 
    return 0;
  }
  
  return geometry->GetSVMap();
}  

//______________________________________________________________________________
Int_t                          
AliMUONVGeometryBuilder::GetModuleId(const TString& envName) const
{
/// Return module Id which has the envelope with given name

  for (Int_t i=0; i<fGeometryModules->GetEntriesFast(); i++) {

    AliMUONGeometryModule* geometry 
      = (AliMUONGeometryModule*)fGeometryModules->At(i);
      
    if ( geometry->GetEnvelopeStore()->FindEnvelope(envName) ) 
      return geometry->GetModuleId();
  }   
  
  return -1;
}  


//______________________________________________________________________________
void AliMUONVGeometryBuilder::SetTranslation(Int_t moduleId, 
                                  const TGeoTranslation& translation)
{
/// Set the translation to the geometry module given by moduleId,
/// apply reference frame transformation 

  AliMUONGeometryModule* geometry = GetGeometry(moduleId);
  
  if (!geometry) {
    AliFatal(Form("Geometry %d is not defined", moduleId)); 
    return;
  }
  
  // Apply frame transform
  TGeoHMatrix newTransform = ConvertTransform(translation);

  // Set new transformation
  geometry->SetTransformation(newTransform);
}  


//______________________________________________________________________________
void AliMUONVGeometryBuilder::SetTransformation(Int_t moduleId, 
                                  const TGeoTranslation& translation,
				  const TGeoRotation& rotation)
{
/// Set the transformation to the geometry module given by moduleId,
/// apply reference frame transformation 

  AliMUONGeometryModule* geometry = GetGeometry(moduleId);
  
  if (!geometry) {
    AliFatal(Form("Geometry %d is not defined", moduleId)); 
    return;
  }
  
  TGeoCombiTrans transformation 
    = TGeoCombiTrans(translation, rotation);

  // Apply frame transform
  TGeoHMatrix newTransform = ConvertTransform(transformation);

  // Set new transformation
  geometry->SetTransformation(newTransform);
}  

//______________________________________________________________________________
void AliMUONVGeometryBuilder::SetVolume(Int_t moduleId, 
                                 const TString& volumeName, 
				 Bool_t isVirtual)
{
/// Set volume name, virtuality

  TString path = GetGeometry(moduleId)->GetVolumePath();
  // cout << "in AliMUONVGeometryBuilder::SetVolume " << path.Data() << endl;
  
  if ( path == "" ) path = "/ALIC_1";
  path += ComposePath(volumeName, 1);

  GetGeometry(moduleId)->SetVolumePath(path);
  GetGeometry(moduleId)->SetIsVirtual(isVirtual);
  // cout << "... set " << path.Data() << endl;
}  				 

//______________________________________________________________________________
void AliMUONVGeometryBuilder::SetMotherVolume(Int_t moduleId, 
                                 const TString& volumeName)
{
/// Set mother volume name

  TString motherVolumeName = ComposePath(volumeName, 1);

  TString path = GetGeometry(moduleId)->GetVolumePath();
  if ( path == "" ) path = "/ALIC_1";
  path.Insert(7, motherVolumeName);  
  
  GetGeometry(moduleId)->SetVolumePath(path);
}  				 

//
// public functions
//

//______________________________________________________________________________
void  AliMUONVGeometryBuilder::SetReferenceFrame(
                                  const TGeoCombiTrans& referenceFrame)
{ 
/// Set reference frame to builder and to all associated geometry 
/// modules

  fReferenceFrame = referenceFrame; 

  for (Int_t i=0; i<fGeometryModules->GetEntriesFast(); i++) {
    AliMUONGeometryModule* geometry 
      = (AliMUONGeometryModule*)fGeometryModules->At(i);
    AliMUONGeometryEnvelopeStore* envelopeStore 
      = geometry->GetEnvelopeStore();
      
    envelopeStore->SetReferenceFrame(referenceFrame);
  }          
}


//______________________________________________________________________________
void  AliMUONVGeometryBuilder::CreateDetElements() const
{
/// Create detection elements and fill their global and
/// local transformations from geometry.

  for (Int_t i=0; i<fGeometryModules->GetEntriesFast(); i++) {
    AliMUONGeometryModule* geometry 
      = (AliMUONGeometryModule*)fGeometryModules->At(i);
      
    const TObjArray* envelopes 
      = geometry->GetEnvelopeStore()->GetEnvelopes();    
    
    AliMpExMap* detElements 
      = geometry->GetTransformer()->GetDetElementStore(); 
      
    for (Int_t j=0; j<envelopes->GetEntriesFast(); j++) {
      AliMUONGeometryEnvelope* envelope
        = (AliMUONGeometryEnvelope*)envelopes->At(j);

      // skip envelope not corresponding to detection element
      if ( envelope->GetUniqueID() == 0) continue;
       
      // Get envelope data 
      Int_t detElemId = envelope->GetUniqueID();	

      // Compose full volume path
      TString volPath = geometry->GetVolumePath();
      volPath += ComposePath(envelope->GetName(), envelope->GetCopyNo());

      // Create detection element 
      AliMUONGeometryDetElement* detElement
        = new AliMUONGeometryDetElement(detElemId, volPath);
      detElements->Add(detElemId, detElement);
      
      // Compose  local transformation
      const TGeoCombiTrans* transform = envelope->GetTransformation(); 
      // Apply frame transform
      TGeoHMatrix localTransform = ConvertDETransform(*transform);
      detElement->SetLocalTransformation(localTransform);

      // Compose global transformation
      TGeoHMatrix globalTransform 
	= AliMUONGeometryBuilder::Multiply( 
	            (*geometry->GetTransformer()->GetTransformation()),
	             localTransform );
		    ;
      // Set the global transformation to detection element
      detElement->SetGlobalTransformation(globalTransform);
      
    }  
  }
}
//_____ _________________________________________________________________________
void  AliMUONVGeometryBuilder::RebuildSVMaps(Bool_t withEnvelopes) const
{
/// Clear the SV maps in memory and fill them from defined geometry.

  for (Int_t i=0; i<fGeometryModules->GetEntriesFast(); i++) {
    AliMUONGeometryModule* geometry 
      = (AliMUONGeometryModule*)fGeometryModules->At(i);
    
    // Clear the map   
    geometry->GetSVMap()->Clear();
     
    // Fill the map from geometry
    const TObjArray* envelopes 
      = geometry->GetEnvelopeStore()->GetEnvelopes();    

    for (Int_t j=0; j<envelopes->GetEntriesFast(); j++) {
      AliMUONGeometryEnvelope* envelope
        = (AliMUONGeometryEnvelope*)envelopes->At(j);

      // skip envelope not corresponding to detection element
      if ( envelope->GetUniqueID() == 0 ) continue;
      
      // Get volume path of detection element
      AliMUONGeometryDetElement* detElement
        = geometry->GetTransformer()->GetDetElement(envelope->GetUniqueID());
      std::string path0 = detElement->GetVolumePath().Data();	
	
      if ( ! withEnvelopes && geometry->IsVirtual() ) {
         std::string vName = geometry->GetTransformer()->GetVolumeName().Data();
	 std::string vPath = ComposePath(vName, 1).Data();
	 path0.erase(path0.find(vPath), vPath.size());
      }  
       
      if ( ! withEnvelopes && envelope->IsVirtual()) {
         std::string eName = envelope->GetName();
	 std::string ePath = ComposePath(eName, envelope->GetCopyNo()).Data();
	 path0.erase(path0.find(ePath), ePath.size());
      }

      if ( ! envelope->IsVirtual() )
        MapSV(path0, envelope->GetName(), envelope->GetUniqueID());
      else { 	
        for  (Int_t k=0; k<envelope->GetConstituents()->GetEntriesFast(); k++) {
          AliMUONGeometryConstituent* constituent
            = (AliMUONGeometryConstituent*)envelope->GetConstituents()->At(k);
         TString path = path0;
	 path += ComposePath(constituent->GetName(), constituent->GetCopyNo());
	 MapSV(path, constituent->GetName(), envelope->GetUniqueID());
        }
      }
    }  
  } 	             
}

