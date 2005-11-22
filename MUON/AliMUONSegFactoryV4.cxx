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

////////////////////////////////////////////////////////////
//  Factory for muon chambers, segmentations and response //
////////////////////////////////////////////////////////////

/* $Id$ */

#include <Riostream.h>
#include <TObjString.h>

#include "AliRun.h"
#include "AliLog.h"

#include "AliMpPlaneType.h"
#include "AliMpExMap.h"

#include "AliMUONSegFactoryV4.h"
#include "AliMUONConstants.h"
#include "AliMUONGeometryTransformer.h"
#include "AliMUONGeometryModule.h"
#include "AliMUONSegmentation.h"
#include "AliMUONGeometrySegmentation.h"
#include "AliMUONSegmentationManager.h"
#include "AliMUONSt12QuadrantSegmentation.h"
#include "AliMUONSt345SlatSegmentationV2.h"
#include "AliMUONTriggerSegmentationV2.h"

ClassImp(AliMUONSegFactoryV4)

//__________________________________________________________________________
AliMUONSegFactoryV4::AliMUONSegFactoryV4(const char* name)
    : TNamed(name, ""),
      fSegmentation(0),
      fkGeomTransformer(0)
{  
/// Standard constructor

  fSegmentation = new AliMUONSegmentation(AliMUONConstants::NCh());
}

//__________________________________________________________________________
  AliMUONSegFactoryV4::AliMUONSegFactoryV4()
    : TNamed(),
      fSegmentation(0),
      fkGeomTransformer(0)
{
/// Default constructor
}

//__________________________________________________________________________
AliMUONSegFactoryV4::AliMUONSegFactoryV4(const AliMUONSegFactoryV4& rhs)
 : TNamed(rhs)
{
/// Protected copy constructor

  AliFatal("Not implemented.");
}

//__________________________________________________________________________

AliMUONSegFactoryV4::~AliMUONSegFactoryV4()
{
/// Destructor


  //delete fSegmentation;
        // The segmentation is supposed to be deleted in the client code
}

//__________________________________________________________________________
AliMUONSegFactoryV4&  AliMUONSegFactoryV4::operator=(const AliMUONSegFactoryV4& rhs)
{
  // Protected assignement operator

  if (this == &rhs) return *this;

  AliFatal("Not implemented.");
    
  return *this;  
}    
          
//
// Private methods
//

//__________________________________________________________________________
Bool_t AliMUONSegFactoryV4::IsGeometryDefined(Int_t ichamber)
{
// Return true, if det elements for the chamber with the given ichamber Id
// are defined in geometry (the geometry builder for this chamber was activated)

  if ( ! fkGeomTransformer ||
       ! fkGeomTransformer->GetModuleTransformer(ichamber, false) )
       
    return kFALSE;
  
  return kTRUE;
}  

//_____________________________________________________________________________
Bool_t
AliMUONSegFactoryV4::ReadDENames(const TString& fileName, AliMpExMap& map)
{ 
/// Read det element names from the file specified by name
/// and fill the map 

  // Open file
  TString filePath(gSystem->ExpandPathName("${ALICE_ROOT}/MUON/data/"));
  filePath += fileName;
  std::ifstream in(filePath);
  if (!in.good()) {
    AliErrorClass(Form("Cannot read file %s", filePath.Data()));
    return false;
  }
  
  // Read file and fill the map
  char line[80];
  while ( in.getline(line,80) )
  {    
    if ( !isdigit(line[0]) ) continue;
    TString sline(line);
    
    Ssiz_t pos = sline.First(' ');
    Int_t detelemid = TString(sline(0,pos)).Atoi();
    TObject* o = map.GetValue(detelemid);
    if (!o)
    {
      map.Add(detelemid, new TObjString(sline(pos+1,sline.Length()-pos).Data()));
    }
  }
  
  // Close file
  in.close();
  return true;
}

//_____________________________________________________________________________
void
AliMUONSegFactoryV4::BuildChamber345(Int_t firstDetElemId, Int_t lastDetElemId,
                                     const AliMpExMap& deNamesMap)
{
  // Build a single chamber for stations 345.
  // The first and lastDetElemId must correspond to the same chamber.
	
  Int_t ichamber = firstDetElemId/100 - 1;
  Int_t test = lastDetElemId/100-1;
  
  if ( test != ichamber )
	{
		AliFatal(Form("DetElemIds %d and %d not part of the same chamber !",
									firstDetElemId,lastDetElemId));
	}
	
  const Int_t kNPLANES = 2;
  const AliMpPlaneType kptypes[kNPLANES] = { kBendingPlane, kNonBendingPlane };
  
  const AliMUONGeometryModuleTransformer* kModuleTransformer 
    = fkGeomTransformer->GetModuleTransformer(ichamber);
	
  for ( Int_t iplane = 0; iplane < kNPLANES; ++iplane )
	{
		AliMUONGeometrySegmentation* segmentation = 
		new AliMUONGeometrySegmentation(kModuleTransformer);
		
		for ( Int_t d = firstDetElemId; d <= lastDetElemId; ++d ) 
		{
			if ( !deNamesMap.GetValue(d) )
	    {
	      AliWarning(Form("You are requesting an invalid detElemId = %d, I am skipping it",d));
	      continue;
	    }
			
			AliMUONVGeometryDESegmentation* slatSeg = 
	    new AliMUONSt345SlatSegmentationV2(d,kptypes[iplane]);
	    
			fSegmentation->AddDESegmentation(slatSeg);

			TString deName = ((TObjString*)deNamesMap.GetValue(d))->GetString();
			segmentation->Add(d, deName, slatSeg);
		}
		
		fSegmentation->AddModuleSegmentation(ichamber, iplane, segmentation);
	}
}

//_____________________________________________________________________________
void
AliMUONSegFactoryV4::BuildChamberTrigger(Int_t firstDetElemId, Int_t lastDetElemId,
                                         const AliMpExMap& deNamesMap)
{
  // Build a single chamber for trigger stations.
  // The first and lastDetElemId must correspond to the same chamber.
	
  Int_t ichamber = firstDetElemId/100 - 1;
  Int_t test = lastDetElemId/100-1;
  
  if ( test != ichamber )
	{
		AliFatal(Form("DetElemIds %d and %d not part of the same chamber !",
									firstDetElemId,lastDetElemId));
	}
	
  const Int_t kNPLANES = 2;
  const AliMpPlaneType kptypes[kNPLANES] = { kBendingPlane, kNonBendingPlane };
  
  const AliMUONGeometryModuleTransformer* kModuleTransformer 
    = fkGeomTransformer->GetModuleTransformer(ichamber);
	
  for ( Int_t iplane = 0; iplane < kNPLANES; ++iplane )
	{
		AliMUONGeometrySegmentation* segmentation = 
		new AliMUONGeometrySegmentation(kModuleTransformer);
		
		for ( Int_t d = firstDetElemId; d <= lastDetElemId; ++d ) 
		{
			if ( !deNamesMap.GetValue(d) )
	    {
	      AliWarning(Form("You are requesting an invalid detElemId = %d, I am skipping it",d));
	      continue;
	    }
			
			AliMUONVGeometryDESegmentation* slatSeg = 
	    new AliMUONTriggerSegmentationV2(d,kptypes[iplane]);
			      
			fSegmentation->AddDESegmentation(slatSeg);
			
			TString deName = ((TObjString*)deNamesMap.GetValue(d))->GetString();
			segmentation->Add(d, deName, slatSeg);
		}
		
		fSegmentation->AddModuleSegmentation(ichamber, iplane, segmentation);
	}
}

//__________________________________________________________________________
void AliMUONSegFactoryV4::BuildStation1() 
{
/// Station 1 

  // Quadrant segmentations:
  AliMUONSt12QuadrantSegmentation* bendSt1
    = new AliMUONSt12QuadrantSegmentation(kStation1, kBendingPlane);
  AliMUONSt12QuadrantSegmentation* nonbendSt1
    = new AliMUONSt12QuadrantSegmentation(kStation1, kNonBendingPlane);
  
  // Add in the array (for safe deleting)  
  fSegmentation->AddDESegmentation(bendSt1);  
  fSegmentation->AddDESegmentation(nonbendSt1);  

  AliMUONGeometrySegmentation* segmentation[2];

  for (Int_t chamber = 0; chamber < 2; chamber++) {

    const AliMUONGeometryModuleTransformer* kModuleTransformer 
      = fkGeomTransformer->GetModuleTransformer(chamber);
      
    segmentation[0] = new AliMUONGeometrySegmentation(kModuleTransformer);
    segmentation[1] = new AliMUONGeometrySegmentation(kModuleTransformer);
        
    // id detection elt for chamber 1
    Int_t id0 = (chamber+1)*100;

    // cathode 0
    segmentation[0]->Add(id0,      "St1_Quadrant_I",   bendSt1);
    segmentation[0]->Add(id0 +  1, "St1_Quadrant_II",  nonbendSt1); 
    segmentation[0]->Add(id0 +  2, "St1_Quadrant_III", bendSt1);
    segmentation[0]->Add(id0 +  3, "St1_Quadrant_IV",  nonbendSt1);
    fSegmentation->AddModuleSegmentation(chamber, 0, segmentation[0]);   

    // cathode 1
    segmentation[1]->Add(id0,      "St1_Quadrant_I",   nonbendSt1);
    segmentation[1]->Add(id0 +  1, "St1_Quadrant_II",  bendSt1);
    segmentation[1]->Add(id0 +  2, "St1_Quadrant_III", nonbendSt1);
    segmentation[1]->Add(id0 +  3, "St1_Quadrant_IV",  bendSt1);
    fSegmentation->AddModuleSegmentation(chamber, 1, segmentation[1]);
  }
}

//__________________________________________________________________________
void AliMUONSegFactoryV4::BuildStation2() 
{
  //
  //--------------------------------------------------------
  // Configuration for Chamber TC3/4 (Station 2) -----------
  ///^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

  // Quadrant segmentations:
  AliMUONSt12QuadrantSegmentation* bendSt2
    = new AliMUONSt12QuadrantSegmentation(kStation2, kBendingPlane);
  AliMUONSt12QuadrantSegmentation* nonbendSt2
    = new AliMUONSt12QuadrantSegmentation(kStation2, kNonBendingPlane);

  // Add in the array (for safe deleting)  
  fSegmentation->AddDESegmentation(bendSt2);  
  fSegmentation->AddDESegmentation(nonbendSt2);  

  AliMUONGeometrySegmentation* segmentation[2];

  for (Int_t chamber = 2; chamber < 4; chamber++) {

    const AliMUONGeometryModuleTransformer* kModuleTransformer 
      = fkGeomTransformer->GetModuleTransformer(chamber);
      
    segmentation[0] = new AliMUONGeometrySegmentation(kModuleTransformer);
    segmentation[1] = new AliMUONGeometrySegmentation(kModuleTransformer);
        
    // id detection elt for chamber 1
    Int_t id0 = (chamber+1)*100;

    //--------------------------------------------------------
    // Configuration for Chamber TC3/4  (Station 2) ----------           


    // cathode 0
    segmentation[0]->Add(id0,      "St2_Quadrant_I",   bendSt2);
    segmentation[0]->Add(id0 +  1, "St2_Quadrant_II",  nonbendSt2); 
    segmentation[0]->Add(id0 +  2, "St2_Quadrant_III", bendSt2);
    segmentation[0]->Add(id0 +  3, "St2_Quadrant_IV",  nonbendSt2);
    fSegmentation->AddModuleSegmentation(chamber, 0, segmentation[0]);   

    // cathode 1
    segmentation[1]->Add(id0,      "St2_Quadrant_I",   nonbendSt2);
    segmentation[1]->Add(id0 +  1, "St2_Quadrant_II",  bendSt2);
    segmentation[1]->Add(id0 +  2, "St2_Quadrant_III", nonbendSt2);
    segmentation[1]->Add(id0 +  3, "St2_Quadrant_IV",  bendSt2);
    fSegmentation->AddModuleSegmentation(chamber, 1, segmentation[1]);
  }
}       
        
//__________________________________________________________________________
void AliMUONSegFactoryV4::BuildStation3(const AliMpExMap& deNamesMap) 
{
  BuildChamber345(500,517,deNamesMap);
  BuildChamber345(600,617,deNamesMap);
}

//__________________________________________________________________________
void AliMUONSegFactoryV4::BuildStation4(const AliMpExMap& deNamesMap) 
{
  BuildChamber345(700,725,deNamesMap);
  BuildChamber345(800,825,deNamesMap);
}

//__________________________________________________________________________
void AliMUONSegFactoryV4::BuildStation5(const AliMpExMap& deNamesMap) 
{       
  BuildChamber345(900,925,deNamesMap);
  BuildChamber345(1000,1025,deNamesMap);
}

//__________________________________________________________________________
void AliMUONSegFactoryV4::BuildStation6(const AliMpExMap& deNamesMap) 
{ 
  BuildChamberTrigger(1100,1117,deNamesMap);
  BuildChamberTrigger(1200,1217,deNamesMap);
  BuildChamberTrigger(1300,1317,deNamesMap);
  BuildChamberTrigger(1400,1417,deNamesMap);
}

//__________________________________________________________________________
void AliMUONSegFactoryV4::Build(const AliMUONGeometryTransformer* geometry) 
{
/// Construct segmentation for all MUON stations
//

  fkGeomTransformer = geometry;

  AliMpExMap map1(kTRUE);
  ReadDENames("denames_slat.dat", map1);

  AliMpExMap map2(kTRUE);
  ReadDENames("denames_trigger.dat", map2);

  // Build all stations
  if (IsGeometryDefined(0))  BuildStation1();
  if (IsGeometryDefined(2))  BuildStation2();
  if (IsGeometryDefined(4))  BuildStation3(map1);
  if (IsGeometryDefined(6))  BuildStation4(map1);
  if (IsGeometryDefined(8))  BuildStation5(map1);
  if (IsGeometryDefined(10)) BuildStation6(map2);
}

//__________________________________________________________________________
void AliMUONSegFactoryV4::BuildStation(
                               const AliMUONGeometryTransformer* geometry, 
                               Int_t stationNumber) 
{
/// Construct segmentations for the given MUON station

  fkGeomTransformer = geometry;

  AliMpExMap map1(kTRUE);
  ReadDENames("denames_slat.dat", map1);

  AliMpExMap map2(kTRUE);
  ReadDENames("denames_trigger.dat", map2);

  switch (stationNumber) {    
    case 1:  BuildStation1(); break;
    case 2:  BuildStation2(); break;
    case 3:  BuildStation3(map1); break;
    case 4:  BuildStation4(map1); break;
    case 5:  BuildStation5(map1); break;
    case 6:  BuildStation6(map2); break;
    
    default: AliFatal("Wrong station number");
  }  
}         
