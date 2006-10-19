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
// Class AliMUONGeometry
// ----------------------------
// Manager class for geometry construction via geometry builders.
// Author: Ivana Hrivnacova, IPN Orsay

#include "AliMUONGeometry.h"
#include "AliMUONGeometryTransformer.h"
#include "AliMUONGeometryModule.h"
#include "AliMUONStringIntMap.h"

#include "AliMpDEManager.h"

#include "AliLog.h"

#include <TObjArray.h>
#include <Riostream.h>
#include <TSystem.h>

#include <iostream>

/// \cond CLASSIMP
ClassImp(AliMUONGeometry)
/// \endcond
 
//______________________________________________________________________________
AliMUONGeometry::AliMUONGeometry(Bool_t isOwner)
  : TObject(),
    fModules(0),
    fTransformer(0)
    
{
/// Standard constructor

  // Create array for geometry modules
  fModules = new TObjArray();
  fModules->SetOwner(isOwner);
  
  // Geometry parametrisation
  fTransformer = new AliMUONGeometryTransformer(false);  
}

//______________________________________________________________________________
AliMUONGeometry::AliMUONGeometry() 
  : TObject(),
    fModules(0),
    fTransformer(0)
{
/// Default constructor
} 

//______________________________________________________________________________
AliMUONGeometry::~AliMUONGeometry()
{
/// Destructor

  delete fModules;
  delete fTransformer;
}

//
// private methods
//

//______________________________________________________________________________
TString  AliMUONGeometry::ComposePath(const TString& volName, 
                                       Int_t copyNo) const
{
/// Compose path from given volName and copyNo

  TString path(volName);
  path += ".";
  path += copyNo;
  
  return path;
}  

//______________________________________________________________________________
void AliMUONGeometry::FillData3(const TString& sensVolumePath, 
                                Int_t detElemId)
{
/// Fill the mapping of the sensitive volume path to the detection element.

  // Module Id
  Int_t moduleId = AliMpDEManager::GetGeomModuleId(detElemId);
    
  // Get module
  AliMUONGeometryModule* module 
    = (AliMUONGeometryModule*)fModules->At(moduleId);
    
  if ( !module ) {
    AliWarningStream()
      << "Geometry module for det element " << detElemId << " not defined."
      << endl;
    return;
  }    
    
  // Get module sensitive volumes map
  AliMUONStringIntMap* svMap = module->GetSVMap();     

  // Map the sensitive volume to detection element
  svMap->Add(sensVolumePath, detElemId); 
}		   
  
//______________________________________________________________________________
TString  AliMUONGeometry::ReadData3(ifstream& in)
{
/// Read SV maps from a file.
/// Return true, if reading finished correctly.

  TString key("SV");
  while ( key == TString("SV") ) {

    // Input data
    TString   volumePath;
    Int_t     detElemId;
  
    in >> volumePath;
    in >> detElemId;

    //cout << "volumePath=" << volumePath << "  "
    //	 << "detElemId=" << detElemId 	 
    //     << endl;   

    // Fill data
    FillData3(volumePath, detElemId); 
     
    // Go to next line
    in >> key;
  } 
  
  return key;
}

//______________________________________________________________________________
void AliMUONGeometry::WriteData3(ofstream& out) const
{
/// Write association of sensitive volumes and detection elements
/// from the sensitive volume map

  for (Int_t i=0; i<fModules->GetEntriesFast(); i++) {
    AliMUONGeometryModule* geometry 
      = (AliMUONGeometryModule*)fModules->At(i);
    AliMUONStringIntMap* svMap
      = geometry->GetSVMap();

    svMap->Print("SV", out);
    out << endl;  
  }    
}

//
// public functions
//

//_____________________________________________________________________________
void AliMUONGeometry::AddModule(AliMUONGeometryModule* module)
{
/// Add the geometry module to the array

  fModules->Add(module);

  if (module)
    fTransformer->AddModuleTransformer(module->GetTransformer());
}

//______________________________________________________________________________
Bool_t  
AliMUONGeometry::ReadSVMap(const TString& fileName)
{
/// Read the sensitive volume maps from a file.
/// Return true, if reading finished correctly.

  // No reading
  // if builder is not associated with any geometry module
  if (fModules->GetEntriesFast() == 0) return false;

  // File path
  TString filePath = gSystem->Getenv("ALICE_ROOT");
  filePath += "/MUON/data/";
  filePath += fileName;
  
  // Open input file
  ifstream in(filePath, ios::in);
  if (!in) {
    cerr << filePath << endl;	
    AliFatal("File not found.");
    return false;
  }

  TString key;
  in >> key;
  while ( !in.eof() ) {
    if (key == TString("SV")) 
      key = ReadData3(in);
    else {
      AliFatal(Form("%s key not recognized",  key.Data()));
      return false;
    }
  }     

  return true;
}

//______________________________________________________________________________
Bool_t  
AliMUONGeometry::WriteSVMap(const TString& fileName) const
{
/// Write sensitive volume map into a file.
/// Return true, if writing finished correctly.

  // No writing
  // if builder is not associated with any geometry module
  if (fModules->GetEntriesFast() == 0) return false;

  // File path
  TString filePath = gSystem->Getenv("ALICE_ROOT");
  filePath += "/MUON/data/";
  filePath += fileName;
  
  // Open input file
  ofstream out(filePath, ios::out);
  if (!out) {
    cerr << filePath << endl;	
    AliError("File not found.");
    return false;
  }
#if !defined (__DECCXX)
  out.setf(std::ios::fixed);
#endif  
  WriteData3(out);
  
  return true;
}  

//_____________________________________________________________________________
const AliMUONGeometryModule* 
AliMUONGeometry::GetModule(Int_t index, Bool_t warn) const
{
/// Return the geometry module specified by index

  if (index < 0 || index >= fModules->GetEntriesFast()) {
    if (warn) {
      AliWarningStream() 
        << "Index: " << index << " outside limits" << std::endl;
    }			 
    return 0;  
  }  

  return (const AliMUONGeometryModule*) fModules->At(index);
}    

//_____________________________________________________________________________
const AliMUONGeometryModule* 
AliMUONGeometry::GetModuleByDEId(Int_t detElemId, Bool_t warn) const
{
/// Return the geometry module specified by detElemId

  // Get module index
  Int_t index = AliMpDEManager::GetGeomModuleId(detElemId);

  return GetModule(index, warn);
}    
