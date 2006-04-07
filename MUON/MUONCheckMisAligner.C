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

/* $Id$ */

/*
  MUONCheckMisAligner: 
  
  This macro performs the misalignment on an existing muon arm geometry
  based on the standard definition of the detector elements in 
  $ALICE_ROOT/MUON/data

  It uses AliMUONGeometryAligner : 
  --> creates a new AliMUONGeometryTransformer and AliMUONGeometryAligner
  --> reads the transformations in from the transform.dat file (make sure that
  this file is the _standard_ one by comparing it to the one in CVS)
  --> creates a second AliMUONGeometryTransformer by misaligning the existing 
  one using AliMUONAligner::MisAlign
  --> User has to specify the magnitude of the alignments, in the Cartesian 
  co-ordiantes (which are used to apply translation misalignments) and in the
  spherical co-ordinates (which are used to apply angular displacements)
  --> User can also set misalignment ranges by hand using the methods : 
  SetMaxCartMisAlig, SetMaxAngMisAlig, SetXYAngMisAligFactor
  (last method takes account of the fact that the misalingment is greatest in 
  the XY plane, since the detection elements are fixed to a support structure
  in this plane. Misalignments in the XZ and YZ plane will be very small 
  compared to those in the XY plane, which are small already - of the order 
  of microns)
  --> Default behavior generates a "residual" misalignment using gaussian
  distributions. Uniform distributions can still be used, see 
  AliMUONGeometryAligner.
  Note : If the detection elements are allowed to be misaligned in all
  directions, this has consequences for the alignment algorithm, which 
  needs to know the number of free parameters. Eric only allowed 3 : 
  x,y,theta_xy, but in principle z and the other two angles are alignable
  as well.  

// Author:Bruce Becker

*/

void MUONCheckMisAligner(Double_t xcartmisaligm = 0.0, Double_t xcartmisaligw = 0.004, 
			 Double_t ycartmisaligm = 0.0, Double_t ycartmisaligw = 0.003, 
			 Double_t angmisaligm = 0.0, Double_t angmisaligw = 0.0023)
{
  
  AliMUONGeometryTransformer *transform = new AliMUONGeometryTransformer(true);
  transform->ReadGeometryData("volpath.dat", "transform.dat");

  AliMUONGeometryMisAligner misAligner(xcartmisaligm,xcartmisaligw,
                                       ycartmisaligm,ycartmisaligw,
				       angmisaligm,angmisaligw);

  // Generate mis alignment data
  AliMUONGeometryTransformer *newTransform = misAligner.MisAlign(transform,true); 
  newTransform->WriteTransformations("transform2.dat");
  newTransform->WriteMisAlignmentData("misalign.root");

  // Apply misAlignment via AliRoot framework
  TGeoManager::Import("geometry.root");
  AliRun::ApplyAlignObjsToGeom(
     const_cast<TClonesArray*>(newTransform->GetMisAlignmentData()));
  // Save new geometry file
  gGeoManager->Export("geometry2.root");
  
  // Extract new transformations
  AliMUONGeometryTransformer* transform3 = new AliMUONGeometryTransformer(true);
  transform3->ReadGeometryData("volpath.dat", "geometry2.root");
  transform3->WriteTransformations("transform3.dat");
               // Check that transform3.dat is equal to transform2.dat

  // Generate misaligned data in local cdb
  TClonesArray* array = newTransform->GetMisAlignmentData();
   
  // CDB manager
  AliCDBManager* cdbManager = AliCDBManager::Instance();
  cdbManager->SetDefaultStorage("local://ResMisAlignCDB");
  
  AliCDBMetaData* cdbData = new AliCDBMetaData();
  cdbData->SetResponsible("Dimuon Offline project");
  cdbData->SetComment("MUON alignment objects with residual misalignment");
  AliCDBId id("MUON/Align/Data", 0, 0); 
  cdbManager->Put(array, id, cdbData);

  // To run simulation with misaligned geometry, you have to set
  // the Align option in Config.C:
  // MUON->SetAlign("transform2.dat");
}




