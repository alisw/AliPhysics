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
// ----------------------------
// Class AliMUONGeometryBuilder
// ----------------------------
// Manager class for geometry construction via geometry builders.
// Author: Ivana Hrivnacova, IPN Orsay

#include "AliMUONGeometryBuilder.h"
#include "AliMUONVGeometryBuilder.h"	
#include "AliMUONGeometry.h"
#include "AliMUONGeometryTransformer.h"
#include "AliMUONGeometryModule.h"	
#include "AliMUONGeometryModuleTransformer.h"	
#include "AliMUONGeometryEnvelope.h"	
#include "AliMUONGeometryEnvelopeStore.h"
#include "AliMUONGeometryDetElement.h"
#include "AliMUONGeometryConstituent.h"

#include "AliMpDEManager.h"

#include "AliModule.h"
#include "AliLog.h"
#include "AliRun.h"

#include <TObjArray.h>
#include <TVirtualMC.h>
#include <TGeoManager.h>

// static data members
 
const TString  AliMUONGeometryBuilder::fgkDefaultVolPathsFileName = "volpath.dat";   
const TString  AliMUONGeometryBuilder::fgkDefaultTransformFileName = "transform.dat";   
const TString  AliMUONGeometryBuilder::fgkDefaultSVMapFileName = "svmap.dat";    
const TString  AliMUONGeometryBuilder::fgkOutFileNameExtension = ".out";    

/// \cond CLASSIMP
ClassImp(AliMUONGeometryBuilder)
/// \endcond

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
    fTransformFileName(fgkDefaultTransformFileName),
    fSVMapFileName(fgkDefaultSVMapFileName),
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
    fTransformFileName(),
    fSVMapFileName(),
    fGlobalTransformation(),
    fGeometryBuilders(0),
    fGeometry(0)
{
/// Default constructor
} 

//______________________________________________________________________________
AliMUONGeometryBuilder::~AliMUONGeometryBuilder()
{
/// Destructor
  
  delete fGeometryBuilders;
  delete fGeometry;
}

//
// private functions
//

//______________________________________________________________________________
void AliMUONGeometryBuilder::PlaceVolume(const TString& name, const TString& mName, 
                            Int_t copyNo, const TGeoHMatrix& matrix, 
			    Int_t npar, Double_t* param, const char* only,
			    Bool_t makeAssembly) const
{
/// Place the volume specified by name with the given transformation matrix

  if (makeAssembly)
    gGeoManager->MakeVolumeAssembly(name.Data());

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
	
  // Place the volume
  if (npar == 0)
    gMC->Gspos(name, copyNo, mName, xyz[0], xyz[1], xyz[2] , krot, only);
  else 
    gMC->Gsposp(name, copyNo, mName, xyz[0], xyz[1], xyz[2] , krot, only,
                param, npar);
} 

//______________________________________________________________________________
void AliMUONGeometryBuilder::CreateGeometryWithTGeo()
{
/// Construct geometry using geometry builders.
/// Virtual modules/envelopes are placed as TGeoVolume assembly

  if (fAlign) {
    // Read transformations from ASCII data file  
    fGeometry->GetTransformer()
      ->ReadGeometryData(fgkDefaultVolPathsFileName, fTransformFileName);
  }    
 
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
      AliMUONGeometryModuleTransformer* transformer= geometry->GetTransformer();
      const TGeoHMatrix* kModuleTransform = transformer->GetTransformation();
      TString volName       = transformer->GetVolumeName(); 
      TString motherVolName = transformer->GetMotherVolumeName(); 
      
      // Place the module volume
      PlaceVolume(volName, motherVolName, 
                  1, *kModuleTransform, 0, 0, "ONLY", geometry->IsVirtual());
  
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
	  
	// Check consistency of detElemId and module Id
	if ( env->GetUniqueID() > 0 && 
	     AliMpDEManager::GetGeomModuleId(env->GetUniqueID()) 
	     != geometry->GetModuleId() ) {
	     
	  AliErrorStream() 
	    << "Detection element " << env->GetUniqueID() 
	    << " is being placed in geometry module " << geometry->GetModuleId()
	    << " but should go in " 
	    << AliMpDEManager::GetGeomModuleId(env->GetUniqueID())
	    <<  endl;
	  AliFatal("Inconsistent IDs");
	}          
	  
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

        // Place envelope in geometry module by composed transformation:
        // [Tglobal] * Tenv
        TGeoHMatrix total 
	  = Multiply( appliedGlobalTransform, 
                     (*kEnvTrans) );
        PlaceVolume(env->GetName(), volName,
	            env->GetCopyNo(), total, 0, 0, only, env->IsVirtual());
	
        if ( env->IsVirtual() )  {
          //  Place constituents in the envelope
          for  (Int_t l=0; l<env->GetConstituents()->GetEntriesFast(); l++) {
            AliMUONGeometryConstituent* constituent
              = (AliMUONGeometryConstituent*)env->GetConstituents()->At(l);
 
            PlaceVolume(constituent->GetName(), env->GetName(),
	                constituent->GetCopyNo(),
		        *constituent->GetTransformation() ,
                        constituent->GetNpar(), constituent->GetParam(), only);
          }
        }
      } // end of loop over envelopes
    } // end of loop over builder geometries
  } // end of loop over builders
}

//______________________________________________________________________________
void AliMUONGeometryBuilder::CreateGeometryWithoutTGeo()
{
/// Construct geometry using geometry builders.
/// Virtual modules/envelopes are not placed

  if (fAlign) {
    // Read transformations from ASCII data file  
    fGeometry->GetTransformer()
      ->ReadGeometryData(fgkDefaultVolPathsFileName, fTransformFileName);
  }     

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
      AliMUONGeometryModuleTransformer* transformer= geometry->GetTransformer();
      const TGeoHMatrix* kModuleTransform = transformer->GetTransformation();
      TString volName       = transformer->GetVolumeName(); 
      TString motherVolName = transformer->GetMotherVolumeName(); 
      
      // Place the module volume
      if ( !geometry->IsVirtual() ) {
          PlaceVolume(volName, motherVolName, 
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
	  
	// Check consistency of detElemId and module Id
	if ( env->GetUniqueID() > 0 && 
	     AliMpDEManager::GetGeomModuleId(env->GetUniqueID()) 
	     != geometry->GetModuleId() ) {
	     
	  AliErrorStream() 
	    << "Detection element " << env->GetUniqueID() 
	    << " is being placed in geometry module " << geometry->GetModuleId()
	    << " but should go in " 
	    << AliMpDEManager::GetGeomModuleId(env->GetUniqueID())
	    <<  endl;
	  AliFatal("Inconsistent IDs");
	}          
	  
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
          //        => place envelope by composed transformation:
          //           Tch * [Tglobal] * Tenv

          // Compound chamber transformation with the envelope one
          if (geometry->IsVirtual()) {
             TGeoHMatrix total 
	       = Multiply( (*kModuleTransform), 
	                    appliedGlobalTransform, 
	                   (*kEnvTrans) );
             PlaceVolume(env->GetName(), motherVolName,
	                 env->GetCopyNo(), total, 0, 0, only);
          }
	  else {
             TGeoHMatrix total 
	       = Multiply( appliedGlobalTransform, 
	                   (*kEnvTrans) );
             PlaceVolume(env->GetName(), volName,
	                 env->GetCopyNo(), total, 0, 0, only);
          }			 
        }

        if (env->IsVirtual() && env->GetConstituents()->GetEntriesFast() > 0 ) {
          // virtual envelope + nof constituents > 0 
          //         => do not place envelope and place constituents
          //            by composed transformation:
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

              PlaceVolume(constituent->GetName(), motherVolName,
	                  constituent->GetCopyNo(), total,
                          constituent->GetNpar(), constituent->GetParam(), only);
            }
	    else { 			  
              TGeoHMatrix total 
	        = Multiply ( appliedGlobalTransform, 
	                     (*kEnvTrans),
	                     (*constituent->GetTransformation()) );

              PlaceVolume(constituent->GetName(), volName,
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

  if ( gMC->IsRootGeometrySupported() && 
       TString(gMC->ClassName()) != "TGeant4" ) {
       
   CreateGeometryWithTGeo();
  } 
  else
   CreateGeometryWithoutTGeo();
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

  // Load alignement data from geometry if geometry is read from Root file
  if ( gAlice->IsRootGeometry() ) {
    fAlign = true;

    fGeometry->GetTransformer()
      ->ReadGeometryData(fgkDefaultVolPathsFileName, gGeoManager);
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
      // Create detection elements from built geometry
      builder->CreateDetElements();
    }  
  }  
}

//______________________________________________________________________________
void AliMUONGeometryBuilder::WriteSVMaps(const TString& fileName, 
                                         Bool_t rebuild)
{
/// Write sensitive volume maps into files per builder

  // Rebuild sv maps
  //
  if (rebuild) 
    for (Int_t i=0; i<fGeometryBuilders->GetEntriesFast(); i++) {

      AliMUONVGeometryBuilder* builder
        = (AliMUONVGeometryBuilder*)fGeometryBuilders->At(i);

      Bool_t writeEnvelopes = false;
      if ( gMC->IsRootGeometrySupported() &&
           TString(gMC->ClassName()) != "TGeant4") writeEnvelopes = true;

      builder->RebuildSVMaps(writeEnvelopes);
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

//_____________________________________________________________________________
void AliMUONGeometryBuilder::SetAlign(const TString& fileName, Bool_t align)
{ 
/// Set the option for alignement and the transformations file name

  fTransformFileName = fileName;
  fAlign = align; 

  for (Int_t i=0; i<fGeometryBuilders->GetEntriesFast(); i++) {

    AliMUONVGeometryBuilder* builder
      = (AliMUONVGeometryBuilder*)fGeometryBuilders->At(i);
    
    SetAlign(builder); 
  }   
}
