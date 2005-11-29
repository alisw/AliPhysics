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
#include "AliMUONGeometry.h"
#include "AliMUONGeometryTransformer.h"
#include "AliMUONGeometryModule.h"	
#include "AliMUONGeometryModuleTransformer.h"	
#include "AliMUONGeometryEnvelope.h"	
#include "AliMUONGeometryEnvelopeStore.h"
#include "AliMUONGeometryDetElement.h"
#include "AliMUONGeometryStore.h"
#include "AliMUONGeometryConstituent.h"
#include "AliModule.h"
#include "AliLog.h"
#include "AliRun.h"


ClassImp(AliMUONGeometryBuilder)
 
// static functions

//______________________________________________________________________________
TGeoHMatrix AliMUONGeometryBuilder::Multiply(const TGeoMatrix& m1, 
                                             const TGeoMatrix& m2)
{
/// Temporary fix for problem with matrix multiplication in Root 5.02/00

  if (m1.IsIdentity() && m2.IsIdentity()) return TGeoHMatrix();
  
  if (m1.IsIdentity()) return m2;
  
  if (m2.IsIdentity()) return m1;
  
  return m1 * m2;
}

//______________________________________________________________________________
TGeoHMatrix AliMUONGeometryBuilder::Multiply(const TGeoMatrix& m1, 
                                             const TGeoMatrix& m2,
                                             const TGeoMatrix& m3)
{					     
/// Temporary fix for problem with matrix multiplication in Root 5.02/00

  if (m1.IsIdentity() && m2.IsIdentity() & m3.IsIdentity())  
    return TGeoHMatrix();
  
  if (m1.IsIdentity()) return Multiply(m2, m3);
  
  if (m2.IsIdentity()) return Multiply(m1, m3);
  
  if (m3.IsIdentity()) return Multiply(m1, m2);
  
  return m1 * m2 * m3;
}

//______________________________________________________________________________
TGeoHMatrix AliMUONGeometryBuilder::Multiply(const TGeoMatrix& m1, 
                                             const TGeoMatrix& m2,
                                             const TGeoMatrix& m3,
                                             const TGeoMatrix& m4)
{					     
/// Temporary fix for problem with matrix multiplication in Root 5.02/00

  if (m1.IsIdentity() && m2.IsIdentity() & m3.IsIdentity() & m4.IsIdentity())  
    return TGeoHMatrix();
  
  if (m1.IsIdentity()) return Multiply(m2, m3, m4);
  
  if (m2.IsIdentity()) return Multiply(m1, m3, m4);
  
  if (m3.IsIdentity()) return Multiply(m1, m2, m4);
  
  if (m4.IsIdentity()) return Multiply(m1, m2, m3);
  
  return m1 * m2 * m3 * m4;
}

//______________________________________________________________________________
AliMUONGeometryBuilder::AliMUONGeometryBuilder(AliModule* module)
  : TObject(),
    fModule(module),
    fAlign(false),
    fGlobalTransformation(), 
    fGeometryBuilders(0),
    fGeometry(0)
{
/// Standard constructor

  fGeometryBuilders = new TObjArray();
  fGeometryBuilders->SetOwner(true);
  
  fGeometry = new AliMUONGeometry(true);

  // Define the global transformation:
  // Transformation from the old ALICE coordinate system to a new one:
  // x->-x, z->-z 
  TGeoRotation* rotGlobal 
    = new TGeoRotation("rotGlobal", 90., 180., 90., 90., 180., 0.);
  fGlobalTransformation = TGeoCombiTrans(0., 0., 0., rotGlobal);
}

//______________________________________________________________________________
AliMUONGeometryBuilder::AliMUONGeometryBuilder() 
  : TObject(),
    fModule(0),
    fAlign(false),
    fGlobalTransformation(),
    fGeometryBuilders(0),
    fGeometry(0)
{
/// Default constructor
} 

//______________________________________________________________________________
AliMUONGeometryBuilder::AliMUONGeometryBuilder(const AliMUONGeometryBuilder& right) 
  : TObject(right) 
{  
/// Copy constructor (not implemented)

  AliFatal("Copy constructor not provided.");
}

//______________________________________________________________________________
AliMUONGeometryBuilder::~AliMUONGeometryBuilder()
{
/// Destructor
  
  delete fGeometryBuilders;
  delete fGeometry;
}

//______________________________________________________________________________
AliMUONGeometryBuilder& 
AliMUONGeometryBuilder::operator=(const AliMUONGeometryBuilder& right)
{
/// Assignement operator (not implemented)

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
/// Place the volume specified by name with the given transformation matrix

  TGeoHMatrix transform(matrix);
  // Do not apply global transformation 
  // if mother volume was already placed in 
  // the new system of coordinates (that is MUON in negative Z)
  // (as it is applied on the mother volume)
  if (mName == TString("DDIP"))
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

//_____________________________________________________________________________
void AliMUONGeometryBuilder::SetAlign(AliMUONVGeometryBuilder* builder)
{
/// Set align option to all geometry modules associated with the builder

  for (Int_t j=0; j<builder->NofGeometries(); j++) {

    AliMUONGeometryModule* geometry = builder->Geometry(j);
  
    geometry->SetAlign(fAlign);
  }  	  
}  	     

//
// public functions
//

//_____________________________________________________________________________
void AliMUONGeometryBuilder::AddBuilder(AliMUONVGeometryBuilder* geomBuilder)
{
/// Add the geometry builder to the list

  fGeometryBuilders->Add(geomBuilder);
  
  // Pass geometry modules created in the to the geometry parametrisation
  for (Int_t i=0; i<geomBuilder->NofGeometries(); i++) {
    fGeometry->AddModule(geomBuilder->Geometry(i));
  }  
  
  if (geomBuilder->ApplyGlobalTransformation())
    geomBuilder->SetReferenceFrame(fGlobalTransformation);
  
  SetAlign(geomBuilder);
}

//______________________________________________________________________________
void AliMUONGeometryBuilder::CreateGeometry()
{
/// Construct geometry using geometry builders.

  if (fAlign) ReadTransformations();

  for (Int_t i=0; i<fGeometryBuilders->GetEntriesFast(); i++) {

    // Get the builder
    AliMUONVGeometryBuilder* builder
      = (AliMUONVGeometryBuilder*)fGeometryBuilders->At(i);

    // Create geometry + envelopes
    //
    builder->CreateGeometry();
    if (!fAlign) builder->SetTransformations();
    
    // Place module volumes and envelopes
    //
    for (Int_t j=0; j<builder->NofGeometries(); j++) {

      AliMUONGeometryModule* geometry = builder->Geometry(j);
      const TGeoCombiTrans* kModuleTransform 
       = geometry->GetTransformer()->GetTransformation();
      
      // Place the module volume
      if ( !geometry->IsVirtual() ) {
          PlaceVolume(geometry->GetVolume(), geometry->GetMotherVolume(), 
	              1, *kModuleTransform, 0, 0, "ONLY");
      }		      
  
      TGeoCombiTrans appliedGlobalTransform;
      if (builder->ApplyGlobalTransformation())
        appliedGlobalTransform = fGlobalTransformation;

      // Loop over envelopes
      const TObjArray* kEnvelopes 
        = geometry->GetEnvelopeStore()->GetEnvelopes();
      for (Int_t k=0; k<kEnvelopes->GetEntriesFast(); k++) {

        // Get envelope
        AliMUONGeometryEnvelope* env 
	  = (AliMUONGeometryEnvelope*)kEnvelopes->At(k);
	  
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
          //           Tch * [Tglobal] * Tenv

          // Compound chamber transformation with the envelope one
          if (geometry->IsVirtual()) {
             TGeoHMatrix total 
	       = Multiply( (*kModuleTransform), 
	                    appliedGlobalTransform, 
	                   (*kEnvTrans) );
             PlaceVolume(env->GetName(), geometry->GetMotherVolume(),
	                 env->GetCopyNo(), total, 0, 0, only);
          }
	  else {
             TGeoHMatrix total 
	       = Multiply( appliedGlobalTransform, 
	                   (*kEnvTrans) );
             PlaceVolume(env->GetName(), geometry->GetVolume(),
	                 env->GetCopyNo(), total, 0, 0, only);
          }			 
        }

        if (env->IsVirtual() && env->GetConstituents()->GetEntriesFast() > 0 ) {
          // virtual envelope + nof constituents > 0 
          //         => do not place envelope and place constituents
          //            in ALICE by composed transformation:
          //            Tch * [Tglobal] * Tenv * Tconst   

          for  (Int_t l=0; l<env->GetConstituents()->GetEntriesFast(); l++) {
            AliMUONGeometryConstituent* constituent
              = (AliMUONGeometryConstituent*)env->GetConstituents()->At(l);
 
            // Compound chamber transformation with the envelope one + the constituent one
            if (geometry->IsVirtual()) {
              TGeoHMatrix total 
	        = Multiply ( (*kModuleTransform),
 	                     appliedGlobalTransform, 
	                     (*kEnvTrans), 
	                     (*constituent->GetTransformation()) );

              PlaceVolume(constituent->GetName(), geometry->GetMotherVolume(),
	                  constituent->GetCopyNo(), total,
                          constituent->GetNpar(), constituent->GetParam(), only);
            }
	    else { 			  
              TGeoHMatrix total 
	        = Multiply ( appliedGlobalTransform, 
	                     (*kEnvTrans),
	                     (*constituent->GetTransformation()) );

              PlaceVolume(constituent->GetName(), geometry->GetVolume(),
	                  constituent->GetCopyNo(), total,
                          constituent->GetNpar(), constituent->GetParam(), only);
            }			  
          }
        }
      } // end of loop over envelopes
    } // end of loop over builder geometries
  } // end of loop over builders
}

//_____________________________________________________________________________
void AliMUONGeometryBuilder::CreateMaterials()
{
/// Construct materials specific to modules via builders
  
  for (Int_t i=0; i<fGeometryBuilders->GetEntriesFast(); i++) {

    // Get the builder
    AliMUONVGeometryBuilder* builder
      = (AliMUONVGeometryBuilder*)fGeometryBuilders->At(i);

    // Create materials with each builder
    if (builder) builder->CreateMaterials();
  }
}

//______________________________________________________________________________
void AliMUONGeometryBuilder::InitGeometry(const TString& svmapFileName)
{
/// Initialize geometry

  // Read alignement data if geometry is read from Root file
  if ( gAlice->IsRootGeometry() ) {
    fAlign = true;
    ReadTransformations();
  }

  // Read sensitive volume map from a file
  fGeometry->ReadSVMap(svmapFileName);
      
  // Set the chamber (sensitive region) GEANT identifier
  //
  for (Int_t i=0; i<fGeometryBuilders->GetEntriesFast(); i++) {

    // Get the builder
    AliMUONVGeometryBuilder* builder
      = (AliMUONVGeometryBuilder*)fGeometryBuilders->At(i);

    // Set sensitive volumes with each builder
    builder->SetSensitiveVolumes();

    if (!fAlign)  {
      // Fill local transformations from built geometry
      builder->FillTransformations();
    }  
  }  
}

//______________________________________________________________________________
void AliMUONGeometryBuilder::ReadTransformations(const TString& fileName)
{
/// Read transformations from ASCII files 
/// and store them in the geometry parametrisation

  // Read transformations
  //
  AliMUONGeometryTransformer* geomTransformer = fGeometry->GetTransformer();
  geomTransformer->ReadTransformations(fileName);
}

//______________________________________________________________________________
void AliMUONGeometryBuilder::WriteTransformations(const TString& fileName)
{
/// Write transformations into files per builder

  AliMUONGeometryTransformer* geomTransformer = fGeometry->GetTransformer();
  geomTransformer->WriteTransformations(fileName);
}

//______________________________________________________________________________
void AliMUONGeometryBuilder::WriteSVMaps(Bool_t rebuild, 
                                         const TString& fileName)
{
/// Write sensitive volume maps into files per builder

  // Rebuild sv maps
  //
  if (rebuild) 
    for (Int_t i=0; i<fGeometryBuilders->GetEntriesFast(); i++) {

      AliMUONVGeometryBuilder* builder
        = (AliMUONVGeometryBuilder*)fGeometryBuilders->At(i);

      builder->RebuildSVMaps();
    }  
    
  // Write maps in file
  fGeometry->WriteSVMap(fileName);
}

//_____________________________________________________________________________
void AliMUONGeometryBuilder::SetAlign(Bool_t align)
{ 
/// Set the option for alignement

  fAlign = align; 

  for (Int_t i=0; i<fGeometryBuilders->GetEntriesFast(); i++) {

    AliMUONVGeometryBuilder* builder
      = (AliMUONVGeometryBuilder*)fGeometryBuilders->At(i);
    
    SetAlign(builder); 
  }   
}
