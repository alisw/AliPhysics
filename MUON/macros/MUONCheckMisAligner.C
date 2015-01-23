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

/// \ingroup macros
/// \file MUONCheckMisAligner.C
/// \brief This macro performs the misalignment on an existing muon arm geometry
///  
/// This macro performs the misalignment on an existing muon arm geometry
/// based on the standard definition of the detector elements in 
/// the AliMUONGeometryTransformer class.
///
/// It uses AliMUONGeometryMisAligner : 
/// - Creates a new AliMUONGeometryTransformer and AliMUONGeometryMisAligner
/// - Loads the geometry from the specified geometry file (default is geometry.root)
/// - Creates a second AliMUONGeometryTransformer by misaligning the existing 
///   one using AliMUONGeometryMisAligner::MisAlign
/// - User has to specify the magnitude of the alignments, in the Cartesian 
///   co-ordiantes (which are used to apply translation misalignments) and in the
///   spherical co-ordinates (which are used to apply angular displacements)
/// - User can also set misalignment ranges by hand using the methods : 
///   SetMaxCartMisAlig, SetMaxAngMisAlig, SetXYAngMisAligFactor
///   (last method takes account of the fact that the misalingment is greatest in 
///   the XY plane, since the detection elements are fixed to a support structure
///   in this plane. Misalignments in the XZ and YZ plane will be very small 
///   compared to those in the XY plane, which are small already - of the order 
///   of microns)
/// - Default behavior generates a "residual" misalignment using gaussian
///   distributions. Uniform distributions can still be used, see 
///   AliMUONGeometryMisAligner.
/// - User can also generate module misalignments using SetModuleCartMisAlig
///   and SetModuleAngMisAlig
///
/// Note: If the detection elements are allowed to be misaligned in all
/// directions, this has consequences for the alignment algorithm, which 
/// needs to know the number of free parameters. Eric only allowed 3 : 
/// x,y,theta_xy, but in principle z and the other two angles are alignable
/// as well.  
///
/// \author:Bruce Becker

#if !defined(__CINT__) || defined(__MAKECINT__)

#include "AliMUONGeometryTransformer.h"
#include "AliMUONGeometryMisAligner.h"

#include "AliGeomManager.h"
#include "AliCDBManager.h"
#include "AliCDBMetaData.h"
#include "AliCDBId.h"

#include <TGeoManager.h>
#include <TClonesArray.h>

#endif

void MUONCheckMisAligner(Double_t xcartmisaligm = 0.0, Double_t xcartmisaligw = 0.004, 
			 Double_t ycartmisaligm = 0.0, Double_t ycartmisaligw = 0.003, 
			 Double_t angmisaligm = 0.0, Double_t angmisaligw = 0.0023,
			 TString nameCDB = "ResMisAlignCDB", 
                         const TString& geomFileName = "geometry.root")
{
  
  AliMUONGeometryTransformer *transform = new AliMUONGeometryTransformer();
  transform->LoadGeometryData(geomFileName.Data());

  AliMUONGeometryMisAligner misAligner(xcartmisaligm,xcartmisaligw,
                                       ycartmisaligm,ycartmisaligw,
				       angmisaligm,angmisaligw);

  // Generate mis alignment data

  // Uncomment lines below if you would like to generate module misalignments
  // misAligner.SetModuleCartMisAlig(0.0,0.1,0.0,0.1,0.0,0.1); // Full
  // misAligner.SetModuleAngMisAlig(0.0,0.02,0.0,0.04,0.0,0.02); // Full
  // misAligner.SetModuleCartMisAlig(0.0,0.003,0.0,0.003,0.0,0.003); // Res
  // misAligner.SetModuleAngMisAlig(0.0,0.0006,0.0,0.001,0.0,0.0005); // Res

  AliMUONGeometryTransformer *newTransform = misAligner.MisAlign(transform,true); 
  newTransform->WriteTransformations("transform2.dat");
  newTransform->WriteMisAlignmentData("misalign.root");

  // Apply misAlignment via AliRoot framework
  AliGeomManager::ApplyAlignObjsToGeom(
     *const_cast<TClonesArray*>(newTransform->GetMisAlignmentData()));
  // Save new geometry file
  gGeoManager->Export("geometry2.root");

  // Extract new transformations
  AliMUONGeometryTransformer* transform3 = new AliMUONGeometryTransformer();
  gGeoManager->UnlockGeometry();
  transform3->LoadGeometryData("geometry2.root");
  transform3->WriteTransformations("transform3.dat");
               // Check that transform3.dat is equal to transform2.dat

  // Generate misaligned data in local cdb
  const TClonesArray* array = newTransform->GetMisAlignmentData();
  
  // 100 mum residual resolution for chamber misalignments?
  misAligner.SetAlignmentResolution(array,-1,0.01,0.01,xcartmisaligw,ycartmisaligw);
  
  TString sLocCDB("local://");
  sLocCDB += nameCDB;
  // CDB manager
  AliCDBManager* cdbManager = AliCDBManager::Instance();
  cdbManager->SetSpecificStorage("MUON/Align/Data",sLocCDB.Data());

  AliCDBMetaData* cdbData = new AliCDBMetaData();
  cdbData->SetResponsible("Dimuon Offline project");
  cdbData->SetComment("MUON alignment objects with residual misalignment");
  AliCDBId id("MUON/Align/Data", 0, AliCDBRunRange::Infinity());
  cdbManager->Put(const_cast<TClonesArray*>(array), id, cdbData);

  // To run simulation with misaligned geometry, you have to set
  // the Align option in Config.C:
  // MUON->SetAlign("transform2.dat");
}




