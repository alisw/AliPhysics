/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *      SigmaEffect_thetadegrees                                                                  *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpeateose. It is      *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

// $Id$
//
// Class AliMUONGeometryBuilder
// ----------------------------
// Manager class for geometry construction via geometry builders.
//
// Author: Ivana Hrivnacova, IPN Orsay

#include <TObjArray.h>
#include <TVirtualMC.h>

#include "AliMUONGeometryBuilder.h"
#include "AliMUONVGeometryBuilder.h"	
#include "AliMUONGeometryModule.h"	
#include "AliMUONGeometryEnvelope.h"	
#include "AliMUONGeometryEnvelopeStore.h"
#include "AliMUONGeometryDetElement.h"
#include "AliMUONGeometryStore.h"
#include "AliMUONGeometryConstituent.h"	
#include "AliModule.h"
#include "AliLog.h"

ClassImp(AliMUONGeometryBuilder)
 
//______________________________________________________________________________
AliMUONGeometryBuilder::AliMUONGeometryBuilder(AliModule* module)
  : TObject(),
    fModule(module),
    fAlign(false),
    fGlobalTransformation(), 
    fGeometryBuilders(0)
{
// Standard constructor

  fGeometryBuilders = new TObjArray(100);
}

//______________________________________________________________________________
AliMUONGeometryBuilder::AliMUONGeometryBuilder() 
  : TObject(),
    fModule(0),
    fAlign(false),
    fGlobalTransformation(),
    fGeometryBuilders(0)
{
// Default constructor
} 

//______________________________________________________________________________
AliMUONGeometryBuilder::AliMUONGeometryBuilder(const AliMUONGeometryBuilder& right) 
  : TObject(right) 
{  
  // copy constructor (not implemented)

  AliFatal("Copy constructor not provided.");
}

//______________________________________________________________________________
AliMUONGeometryBuilder::~AliMUONGeometryBuilder()
{
// Destructor
  if (fGeometryBuilders){
    fGeometryBuilders->Delete();
    delete fGeometryBuilders;
  }
}

//______________________________________________________________________________
AliMUONGeometryBuilder& 
AliMUONGeometryBuilder::operator=(const AliMUONGeometryBuilder& right)
{
  // assignement operator (not implemented)

  // check assignement to self
  if (this == &right) return *this;

  AliFatal("Assignement operator not provided.");
    
  return *this;  
}    

//
// private functions
//

//______________________________________________________________________________
void AliMUONGeometryBuilder::PlaceVolume(const TString& name, const TString& mName, 
                            Int_t copyNo, const TGeoHMatrix& matrix, 
			    Int_t npar, Double_t* param, const char* only) const
{
// Place the volume specified by name with the given transformation matrix
// ---

  // Do not apply global transformation 
  // if mother volume != ALIC
  // (as it is applied on the mother volume)
  TGeoHMatrix transform(matrix);
  if (mName != TString("ALIC"))
    transform = fGlobalTransformation.Inverse() * transform;
     
  // Decompose transformation
  const Double_t* xyz = transform.GetTranslation();
  const Double_t* rm = transform.GetRotationMatrix();
	
  //cout << "Got translation: "
  //     << xyz[0] << " " << xyz[1] << " " << xyz[2] << endl;
	
  //cout << "Got rotation: "
  //     << rm[0] << " " << rm[1] << " " << rm[2] << endl
  //     << rm[3] << " " << rm[4] << " " << rm[5] << endl
  //     << rm[6] << " " << rm[7] << " " << rm[8] << endl;

  // Check for presence of rotation
  // (will be nice to be available in TGeo)
  const Double_t kTolerance = 1e-04;
  Bool_t isRotation = true; 
  if (TMath::Abs(rm[0] - 1.) < kTolerance &&
      TMath::Abs(rm[1] - 0.) < kTolerance &&
      TMath::Abs(rm[2] - 0.) < kTolerance &&
      TMath::Abs(rm[3] - 0.) < kTolerance &&
      TMath::Abs(rm[4] - 1.) < kTolerance &&
      TMath::Abs(rm[5] - 0.) < kTolerance &&
      TMath::Abs(rm[6] - 0.) < kTolerance &&
      TMath::Abs(rm[7] - 0.) < kTolerance &&
      TMath::Abs(rm[8] - 1.) < kTolerance) isRotation = false; 

  Int_t krot = 0;
  if (isRotation) {
    TGeoRotation rot;
    rot.SetMatrix(const_cast<Double_t*>(transform.GetRotationMatrix()));
    Double_t theta1, phi1, theta2, phi2, theta3, phi3;
    rot.GetAngles(theta1, phi1, theta2, phi2, theta3, phi3);
	
    //cout << "angles: " 
    //     << theta1 << " " << phi1 << " "
    //     << theta2 << " " << phi2 << " "
    //     << theta3 << " " << phi3 << endl;
	
    fModule->AliMatrix(krot, theta1, phi1, theta2, phi2, theta3, phi3);
  }	
	
  // Place the volume in ALIC
  if (npar == 0)
    gMC->Gspos(name, copyNo, mName, xyz[0], xyz[1], xyz[2] , krot, only);
  else 
    gMC->Gsposp(name, copyNo, mName, xyz[0], xyz[1], xyz[2] , krot, only,
                param, npar);

} 

//______________________________________________________________________________
void AliMUONGeometryBuilder::FillGlobalTransformations(
                                 AliMUONVGeometryBuilder* builder)
{
// Compute and set global transformations to detection elements
// for each chamber geometry
// ---

  for (Int_t j=0; j<builder->NofGeometries(); j++) {

    AliMUONGeometryModule* geometry = builder->Geometry(j);
    AliMUONGeometryStore* detElements = geometry->GetDetElementStore();

    for (Int_t k=0; k<detElements->GetNofEntries(); k++) {
     
      AliMUONGeometryDetElement* detElement 
	= (AliMUONGeometryDetElement*)detElements->GetEntry(k);
	  
      if (!detElement) AliFatal("Detection element not found.") 
	  
      const TGeoCombiTrans* localTransform 
        = detElement->GetLocalTransformation();

      // Compose global transformation
      TGeoHMatrix total 
	= fGlobalTransformation *
	  (*geometry->GetTransformation()) * 
          fGlobalTransformation.Inverse() *
	  (*localTransform);
	          // !! The local detection element frame is 
		  // defined wrt the new ALICE coordinate system:
		  // TGL = Tglobal * Tchamber * Tde
		  //     =  Tglobal * Tchamber * Tglobal.inv  *  Tglobal * Tde
		  //     = (Tglobal * Tchamber * Tglobal.inv) * (Tglobal * Tde)
		  //     = Ttotal * Tde'
	  
      // Convert TGeoHMatrix to TGeoCombiTrans
      TGeoCombiTrans globalTransform(localTransform->GetName());
      globalTransform.SetTranslation(total.GetTranslation());  
      TGeoRotation rotation;
      rotation.SetMatrix(total.GetRotationMatrix());
      globalTransform.SetRotation(rotation);  
 
      // Set the global transformation to detection element
      detElement->SetGlobalTransformation(globalTransform);
    }
  }  			    
}	     

//
// public functions
//

//_____________________________________________________________________________
void AliMUONGeometryBuilder::AddBuilder(AliMUONVGeometryBuilder* geomBuilder)
{
// Adds the geometry builder to the list
// ---

  fGeometryBuilders->Add(geomBuilder);
}

//______________________________________________________________________________
void AliMUONGeometryBuilder::CreateGeometry()
{
//
// Construct geometry using geometry builders.
//

  for (Int_t i=0; i<fGeometryBuilders->GetEntriesFast(); i++) {

    // Get the builder
    AliMUONVGeometryBuilder* builder
      = (AliMUONVGeometryBuilder*)fGeometryBuilders->At(i);

    // Create geometry + envelopes
    //
    if (fAlign) {
      builder->ReadTransformations();
      builder->CreateGeometry();
    }
    else {  
      builder->CreateGeometry();
      builder->SetTransformations();
    } 

    // Place envelopes
    //
    for (Int_t j=0; j<builder->NofGeometries(); j++) {

      AliMUONGeometryModule* geometry = builder->Geometry(j);
  
      // Loop over envelopes
      const TObjArray* kEnvelopes = geometry->GetEnvelopeStore()->GetEnvelopes();
      for (Int_t k=0; k<kEnvelopes->GetEntriesFast(); k++) {

        // Get envelope
        AliMUONGeometryEnvelope* env = (AliMUONGeometryEnvelope*)kEnvelopes->At(k);
        const TGeoCombiTrans* kEnvTrans = env->GetTransformation();
        const char* only = "ONLY";
        if (env->IsMANY()) only = "MANY";

        if (env->IsVirtual() && env->GetConstituents()->GetEntriesFast() == 0 ) {
          // virtual envelope + nof constituents = 0 
          //         => not allowed;
          //            empty virtual envelope has no sense 
          AliFatal("Virtual envelope must have constituents.");
          return;
        }

        if (!env->IsVirtual() && env->GetConstituents()->GetEntriesFast() > 0 ) {
          // non virtual envelope + nof constituents > 0 
          //        => not allowed;
          //           use VMC to place constituents
          AliFatal("Non virtual envelope cannot have constituents.");
          return;
        }

        if (!env->IsVirtual() && env->GetConstituents()->GetEntriesFast() == 0 ) {
          // non virtual envelope + nof constituents = 0 
          //        => place envelope in ALICE by composed transformation:
          //           Tglobal * Tch * Tenv

          // Compound chamber transformation with the envelope one
          TGeoHMatrix total 
	    = fGlobalTransformation * 
	      (*geometry->GetTransformation()) * 
	      (*kEnvTrans);
          PlaceVolume(env->GetName(), geometry->GetMotherVolume(),
	              env->GetCopyNo(), total, 0, 0, only);
        }

        if (env->IsVirtual() && env->GetConstituents()->GetEntriesFast() > 0 ) {
          // virtual envelope + nof constituents > 0 
          //         => do not place envelope and place constituents
          //            in ALICE by composed transformation:
          //            Tglobal * Tch * Tenv * Tconst   

          for  (Int_t l=0; l<env->GetConstituents()->GetEntriesFast(); l++) {
            AliMUONGeometryConstituent* constituent
              = (AliMUONGeometryConstituent*)env->GetConstituents()->At(l);
 
            // Compound chamber transformation with the envelope one + the constituent one
            TGeoHMatrix total 
	      = fGlobalTransformation *
	        (*geometry->GetTransformation()) * 
	        (*kEnvTrans) * 
	        (*constituent->GetTransformation());

            PlaceVolume(constituent->GetName(), geometry->GetMotherVolume(),
	                constituent->GetCopyNo(), total,
                        constituent->GetNpar(), constituent->GetParam(), only);
          }
        }
      } // end of loop over envelopes
    } // end of loop over builder geometries
  } // end of loop over builders
}

//_____________________________________________________________________________
void AliMUONGeometryBuilder::CreateMaterials()
{
//
// Construct materials specific to modules via builders
//
  
  for (Int_t i=0; i<fGeometryBuilders->GetEntriesFast(); i++) {

    // Get the builder
    AliMUONVGeometryBuilder* builder
      = (AliMUONVGeometryBuilder*)fGeometryBuilders->At(i);

    // Create materials with each builder
    if (builder) builder->CreateMaterials();
  }
}

//______________________________________________________________________________
void AliMUONGeometryBuilder::InitGeometry()
{
 // Initialize geometry
 // ---

  // Set the chamber (sensitive region) GEANT identifier
  //
  for (Int_t i=0; i<fGeometryBuilders->GetEntriesFast(); i++) {

    // Get the builder
    AliMUONVGeometryBuilder* builder
      = (AliMUONVGeometryBuilder*)fGeometryBuilders->At(i);

    // Set sesitive volumes with each builder
    builder->SetSensitiveVolumes();
    
    // Read sensitive volume map from a file
    builder->ReadSVMap();
    if (!fAlign)  builder->FillTransformations();

    // Compute global transformations of detection elements
    FillGlobalTransformations(builder);
  }  
}

//______________________________________________________________________________
void AliMUONGeometryBuilder::WriteTransformations()
{
 // Writes transformations into files per builder
 // ---

  for (Int_t i=0; i<fGeometryBuilders->GetEntriesFast(); i++) {

    // Get the builder
    AliMUONVGeometryBuilder* builder
      = (AliMUONVGeometryBuilder*)fGeometryBuilders->At(i);

    // Write transformations
    builder->WriteTransformations();
  }
}

//______________________________________________________________________________
void AliMUONGeometryBuilder::WriteSVMaps(Bool_t rebuild)
{
 // Writes sensitive volume maps into files per builder
 // ---

  for (Int_t i=0; i<fGeometryBuilders->GetEntriesFast(); i++) {

    // Get the builder
    AliMUONVGeometryBuilder* builder
      = (AliMUONVGeometryBuilder*)fGeometryBuilders->At(i);

    // Write transformations      
    builder->WriteSVMap(rebuild);
  }
}

//_____________________________________________________________________________
void  AliMUONGeometryBuilder::SetGlobalTransformation(
                                       const TGeoCombiTrans& transform)
{
// Sets the global transformation
// ---

  fGlobalTransformation = transform;
}				       

//_____________________________________________________________________________
void AliMUONGeometryBuilder::SetAlign(Bool_t align)
{ 
// Sets the option for alignement
// ---

  fAlign = align; 

  for (Int_t i=0; i<fGeometryBuilders->GetEntriesFast(); i++) {

    // Get the builder
    AliMUONVGeometryBuilder* builder
      = (AliMUONVGeometryBuilder*)fGeometryBuilders->At(i);

    for (Int_t j=0; j<builder->NofGeometries(); j++) {

      AliMUONGeometryModule* geometry = builder->Geometry(j);
  
      geometry->SetAlign(align);
    }  	  
  }	  
}
