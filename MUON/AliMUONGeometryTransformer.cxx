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
// Class AliMUONGeometryTransformer
// ----------------------------
// Top container class for geometry transformations
// Author: Ivana Hrivnacova, IPN Orsay

#include "AliMUONGeometryTransformer.h"
#include "AliMUONGeometryModuleTransformer.h"
#include "AliMUONGeometryDetElement.h"
#include "AliMUONGeometryBuilder.h"

#include "AliMpDEManager.h"
#include "AliMpExMap.h"

#include "AliLog.h"
#include "AliAlignObjMatrix.h"
#include "AliAlignObj.h"

#include <Riostream.h>
#include <TSystem.h>
#include <TClonesArray.h>
#include <TGeoManager.h>
#include <TGeoPhysicalNode.h>
#include <TFile.h>

#include <sstream>

/// \cond CLASSIMP
ClassImp(AliMUONGeometryTransformer)
/// \endcond
 
//______________________________________________________________________________
AliMUONGeometryTransformer::AliMUONGeometryTransformer(Bool_t isOwner,
                                          const TString& detectorName)

  : TObject(),
    fDetectorName(detectorName),
    fModuleTransformers(0),
    fMisAlignArray(0)
{
/// Standard constructor

  // Create array for geometry modules
  fModuleTransformers = new TObjArray(100);
  fModuleTransformers->SetOwner(isOwner);
}

//______________________________________________________________________________
AliMUONGeometryTransformer::AliMUONGeometryTransformer() 
  : TObject(),
    fDetectorName(),
    fModuleTransformers(0),
    fMisAlignArray(0)
{
/// Default constructor
} 

//______________________________________________________________________________
AliMUONGeometryTransformer::~AliMUONGeometryTransformer()
{
/// Destructor

  delete fModuleTransformers;
  delete fMisAlignArray;
}

//
// private methods
//

//_____________________________________________________________________________
AliMUONGeometryModuleTransformer* 
AliMUONGeometryTransformer::GetModuleTransformerNonConst(
                                          Int_t index, Bool_t warn) const
{
/// Return the geometry module specified by index

  if (index < 0 || index >= fModuleTransformers->GetEntriesFast()) {
    if (warn) {
      AliWarningStream() 
        << "Index: " << index << " outside limits" << std::endl;
    }			 
    return 0;  
  }  

  return (AliMUONGeometryModuleTransformer*) fModuleTransformers->At(index);
}    

//______________________________________________________________________________
TGeoHMatrix AliMUONGeometryTransformer::GetTransform(
                  Double_t x, Double_t y, Double_t z,
		  Double_t a1, Double_t a2, Double_t a3, 
 		  Double_t a4, Double_t a5, Double_t a6) const
{		  
/// Build the transformation from the given parameters

  // Compose transform
  return TGeoCombiTrans(TGeoTranslation(x, y, z), 
                        TGeoRotation("rot", a1, a2, a3, a4, a5, a6));
}


//______________________________________________________________________________
void AliMUONGeometryTransformer::FillModuleVolPath(Int_t moduleId,
                                           const TString& volPath) 
{
/// Create module with the given moduleId and volPath

  // Get/Create geometry module transformer
  AliMUONGeometryModuleTransformer* moduleTransformer
    = GetModuleTransformerNonConst(moduleId, false);

  if ( !moduleTransformer ) {
    moduleTransformer = new AliMUONGeometryModuleTransformer(moduleId);
    AddModuleTransformer(moduleTransformer);
  }  
  moduleTransformer->SetVolumePath(volPath);
}		   
  
//______________________________________________________________________________
void AliMUONGeometryTransformer::FillDetElemVolPath(Int_t detElemId, 
                                           const TString& volPath) 
{
/// Create detection element with the given detElemId and volPath

  // Module Id
  Int_t moduleId = AliMpDEManager::GetGeomModuleId(detElemId);

  // Get detection element store
  AliMpExMap* detElements = 
    GetModuleTransformer(moduleId)->GetDetElementStore();     

  // Add detection element
  AliMUONGeometryDetElement* detElement
    = new AliMUONGeometryDetElement(detElemId, volPath);
  detElements->Add(detElemId, detElement);
}		   
  

//______________________________________________________________________________
void AliMUONGeometryTransformer::FillModuleTransform(Int_t moduleId,
                  Double_t x, Double_t y, Double_t z,
		  Double_t a1, Double_t a2, Double_t a3,
 		  Double_t a4, Double_t a5, Double_t a6) 
{
/// Fill the transformation of the module.

  AliMUONGeometryModuleTransformer* moduleTransformer
    = GetModuleTransformerNonConst(moduleId, false);

  if ( !moduleTransformer) {
    AliErrorStream() 
      << "Module " << moduleId << " has not volume path defined." << endl;
  }  
      
  // Build the transformation from the parameters
  TGeoHMatrix transform 
    = GetTransform(x, y, z, a1, a2, a3, a4, a5, a6);
      
  moduleTransformer->SetTransformation(transform);
}		   
  
//______________________________________________________________________________
void AliMUONGeometryTransformer::FillDetElemTransform(
                  Int_t detElemId, 
                  Double_t x, Double_t y, Double_t z,
		  Double_t a1, Double_t a2, Double_t a3,
 		  Double_t a4, Double_t a5, Double_t a6) 
{
/// Fill the transformation of the detection element.

  // Module Id
  Int_t moduleId = AliMpDEManager::GetGeomModuleId(detElemId);

  // Get module transformer
  const AliMUONGeometryModuleTransformer* kModuleTransformer
    = GetModuleTransformer(moduleId);

  if ( ! kModuleTransformer ) {
    AliFatal(Form("Module transformer not defined, detElemId: %d", detElemId));
    return;  
  }  

  // Get detection element
  AliMUONGeometryDetElement* detElement 
    = kModuleTransformer->GetDetElement(detElemId);     

  if ( ! detElement ) {
    AliFatal(Form("Det element %d has not volume path defined", detElemId));
    return;  
  }  
      
  // Build the transformation from the parameters
  TGeoHMatrix localTransform 
    = GetTransform(x, y, z, a1, a2, a3, a4, a5, a6);
  detElement->SetLocalTransformation(localTransform); 
   
  // Compute global transformation
  TGeoHMatrix globalTransform 
    = AliMUONGeometryBuilder::Multiply( 
                                  *kModuleTransformer->GetTransformation(),
				  localTransform );
  detElement->SetGlobalTransformation(globalTransform);
}		   

//______________________________________________________________________________
Bool_t  
AliMUONGeometryTransformer::ReadVolPaths(ifstream& in)
{
/// Read modules and detection element volume paths from stream

  Int_t id;
  TString key, volumePath;
  in >> key;
  
  while ( !in.eof() ) {

    in >> id >> volumePath;

    // cout << "id="     << id << "  "
    // 	 << "volPath= " << volumePath
    //	 << endl;   

    if ( key == AliMUONGeometryModuleTransformer::GetModuleNamePrefix() ) 
      FillModuleVolPath(id, volumePath);
  
    else if ( key == AliMUONGeometryDetElement::GetDENamePrefix() )
      FillDetElemVolPath(id, volumePath);
  
    else {
      AliFatal(Form("%s key not recognized",  key.Data()));
      return false;
    }
    in >> key;
  }     

  return true;
}

//______________________________________________________________________________
TString  AliMUONGeometryTransformer::ReadModuleTransforms(ifstream& in)
{
/// Read and fill modules transformations from the stream.
/// Return true, if reading finished correctly.

  TString key(AliMUONGeometryModuleTransformer::GetModuleNamePrefix());
  while ( key == AliMUONGeometryModuleTransformer::GetModuleNamePrefix() ) {
    Int_t id;
    Double_t  x, y, z;
    Double_t  a1, a2, a3, a4, a5, a6;
    TString dummy;
  
    in >> id;
    in >> dummy;
    in >> x;
    in >> y;
    in >> z;
    in >> dummy;
    in >> a1; 
    in >> a2; 
    in >> a3; 
    in >> a4; 
    in >> a5; 
    in >> a6; 

    //cout << "moduleId="     << id << "  "
    // 	 << "position= " << x << ", " << y << ", " << z << "  "
    //	 << "rotation= " << a1 << ", " << a2 << ", " << a3  << ", "
    //	                 << a4 << ", " << a5 << ", " << a6 
    //	 << endl;   
	 
    // Fill data
    FillModuleTransform(id, x, y, z, a1, a2, a3, a4, a5, a6);
    
    // Go to next line
    in >> key;
  }
  
  return key;   	 
}

//______________________________________________________________________________
TString  AliMUONGeometryTransformer::ReadDetElemTransforms(ifstream& in)
{
/// Read detection elements transformations from the stream.
/// Return true, if reading finished correctly.

  TString key(AliMUONGeometryDetElement::GetDENamePrefix());
  while ( key == AliMUONGeometryDetElement::GetDENamePrefix() ) {

    // Input data
    Int_t detElemId;
    Double_t  x, y, z;
    Double_t  a1, a2, a3, a4, a5, a6;
    TString dummy;
  
    in >> detElemId;
    in >> dummy;
    in >> x;
    in >> y;
    in >> z;
    in >> dummy;
    in >> a1; 
    in >> a2; 
    in >> a3; 
    in >> a4; 
    in >> a5; 
    in >> a6; 

    //cout << "detElemId=" << detElemId << "  "
    //     << "position= " << x << ", " << y << ", " << z << "  "
    //     << "rotation= " << a1 << ", " << a2 << ", " << a3  << ", "
    //	                   << a4 << ", " << a5 << ", " << a6 
    //     << endl;   
	 
    // Fill data
    FillDetElemTransform(detElemId, x, y, z, a1, a2, a3, a4, a5, a6); 	 
    
    // Go to next line
    in >> key;
  } 
  
  return key;
}

//______________________________________________________________________________
Bool_t  
AliMUONGeometryTransformer::LoadTransforms(TGeoManager* tgeoManager)
{
/// Load transformations for defined modules and detection elements
/// from the root file

  if ( !tgeoManager) {
    AliFatal("No TGeoManager defined.");
    return false;
  }   

  for (Int_t i=0; i<fModuleTransformers->GetEntriesFast(); i++) {
    AliMUONGeometryModuleTransformer* moduleTransformer 
      = (AliMUONGeometryModuleTransformer*)fModuleTransformers->At(i);

    // Module path
    TString path = moduleTransformer->GetVolumePath();
    
    // Make physical node
    TGeoPhysicalNode* moduleNode = tgeoManager->MakePhysicalNode(path);
    if ( ! moduleNode ) {
      AliErrorStream() 
        << "Module id: " << moduleTransformer->GetModuleId()
	<< " volume path: " << path << " not found in geometry." << endl;
	return false;
    }	 
    
    // Set matrix from physical node
    TGeoHMatrix matrix = *moduleNode->GetMatrix();
    moduleTransformer->SetTransformation(matrix);
    
    // Loop over detection elements
    AliMpExMap* detElements = moduleTransformer->GetDetElementStore();    
   
    for (Int_t j=0; j<detElements->GetSize(); j++) {
      AliMUONGeometryDetElement* detElement
        = (AliMUONGeometryDetElement*)detElements->GetObject(j);

      // Det element path
      TString dePath = detElement->GetVolumePath();

      // Make physical node
      TGeoPhysicalNode* deNode = tgeoManager->MakePhysicalNode(dePath);
      if ( ! deNode ) {
        AliErrorStream() 
          << "Det element id: " << detElement->GetId()
	  << " volume path: " << path << " not found in geometry." << endl;
	  return false;
      }	
	 
      // Set global matrix from physical node
      TGeoHMatrix globalMatrix = *deNode->GetMatrix();
      detElement->SetGlobalTransformation(globalMatrix);

      // Set local matrix
      TGeoHMatrix localMatrix = 
        AliMUONGeometryBuilder::Multiply(
	   matrix.Inverse(), globalMatrix );
      detElement->SetLocalTransformation(localMatrix);
    }  
  } 
  return true;    
}  

//______________________________________________________________________________
Bool_t  
AliMUONGeometryTransformer::ReadVolPaths(const TString& fileName)
{
/// Read detection element volume paths from a file.
/// Return true, if reading finished correctly.

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

  ReadVolPaths(in);
  return true;
}

//______________________________________________________________________________
Bool_t  
AliMUONGeometryTransformer::ReadTransformations(const TString& fileName)
{
/// Read transformations from a file.
/// Return true, if reading finished correctly.

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
    if ( key == AliMUONGeometryModuleTransformer::GetModuleNamePrefix() ) 
      key = ReadModuleTransforms(in);
    else if ( key == AliMUONGeometryDetElement::GetDENamePrefix() )
      key = ReadDetElemTransforms(in);
    else {
      AliFatal(Form("%s key not recognized",  key.Data()));
      return false;
    }
  }     

  return true;
}

//______________________________________________________________________________
Bool_t  
AliMUONGeometryTransformer::ReadTransformations2(const TString& fileName)
{
/// Read transformations from root geometry file.
/// Return true, if reading finished correctly.

  // File path
  TString filePath = gSystem->Getenv("ALICE_ROOT");
  filePath += "/MUON/data/";
  filePath += fileName;
  
  // Load root geometry
  TGeoManager* tgeoManager = TGeoManager::Import(fileName);

  // Retrieve matrices
  LoadTransforms(tgeoManager);     

  return true;
}

//______________________________________________________________________________
void AliMUONGeometryTransformer::WriteTransform(ofstream& out,
                                   const TGeoMatrix* transform) const
{
/// Write given transformation 

  out << "   pos: ";
  const Double_t* xyz = transform->GetTranslation();
  out << setw(10) << setprecision(4) << xyz[0] << "  " 
      << setw(10) << setprecision(4) << xyz[1] << "  " 
      << setw(10) << setprecision(4) << xyz[2];

  out << "   rot: ";
  const Double_t* rm = transform->GetRotationMatrix();
  TGeoRotation rotation;
  rotation.SetMatrix(const_cast<Double_t*>(rm));
  Double_t a1, a2, a3, a4, a5, a6;
  rotation.GetAngles(a1, a2, a3, a4, a5, a6);
      
  out << setw(8) << setprecision(4) << a1 << "  " 
      << setw(8) << setprecision(4) << a2 << "  " 
      << setw(8) << setprecision(4) << a3 << "  " 
      << setw(8) << setprecision(4) << a4 << "  " 
      << setw(8) << setprecision(4) << a5 << "  " 
      << setw(8) << setprecision(4) << a6 << "  " << endl; 
}

//______________________________________________________________________________
void AliMUONGeometryTransformer::WriteModuleVolPaths(ofstream& out) const
{
/// Write module volume paths for all module transformers

  for (Int_t i=0; i<fModuleTransformers->GetEntriesFast(); i++) {
    AliMUONGeometryModuleTransformer* moduleTransformer 
      = (AliMUONGeometryModuleTransformer*)fModuleTransformers->At(i);

    // Write data on out
    out << AliMUONGeometryModuleTransformer::GetModuleNamePrefix() << " "
        << setw(4) << moduleTransformer->GetModuleId() << "    " 
        << moduleTransformer->GetVolumePath() << endl;
  }     
  out << endl;	  	   	
}

//______________________________________________________________________________
void AliMUONGeometryTransformer::WriteDetElemVolPaths(ofstream& out) const
{
/// Write detection element volume paths for all detection elements in all 
/// module transformers

  for (Int_t i=0; i<fModuleTransformers->GetEntriesFast(); i++) {
    AliMUONGeometryModuleTransformer* moduleTransformer 
      = (AliMUONGeometryModuleTransformer*)fModuleTransformers->At(i);
    AliMpExMap* detElements = moduleTransformer->GetDetElementStore();    

    for (Int_t j=0; j<detElements->GetSize(); j++) {
      AliMUONGeometryDetElement* detElement
        = (AliMUONGeometryDetElement*)detElements->GetObject(j);
	
      // Write data on out
      out << AliMUONGeometryDetElement::GetDENamePrefix() << " " 
          << setw(4) << detElement->GetId() << "    " 
          << detElement->GetVolumePath() << endl;
    }
    out << endl;	  	   	
  }     
}

//______________________________________________________________________________
void AliMUONGeometryTransformer::WriteModuleTransforms(ofstream& out) const
{
/// Write module transformations for all module transformers

  for (Int_t i=0; i<fModuleTransformers->GetEntriesFast(); i++) {
    AliMUONGeometryModuleTransformer* moduleTransformer 
      = (AliMUONGeometryModuleTransformer*)fModuleTransformers->At(i);
    const TGeoMatrix* transform 
      = moduleTransformer->GetTransformation();    

    // Write data on out
    out << AliMUONGeometryModuleTransformer::GetModuleNamePrefix() << " " 
        << setw(4) << moduleTransformer->GetModuleId();
    
    WriteTransform(out, transform);
  }
  out << endl;
}

//______________________________________________________________________________
void AliMUONGeometryTransformer::WriteDetElemTransforms(ofstream& out) const
{
/// Write detection element transformations for all detection elements in all 
/// module transformers

  for (Int_t i=0; i<fModuleTransformers->GetEntriesFast(); i++) {
    AliMUONGeometryModuleTransformer* moduleTransformer 
      = (AliMUONGeometryModuleTransformer*)fModuleTransformers->At(i);
    AliMpExMap* detElements = moduleTransformer->GetDetElementStore();    

    for (Int_t j=0; j<detElements->GetSize(); j++) {
      AliMUONGeometryDetElement* detElement
        = (AliMUONGeometryDetElement*)detElements->GetObject(j);
      const TGeoMatrix* transform 
        = detElement->GetLocalTransformation(); 
	
      // Write data on out
      out << AliMUONGeometryDetElement::GetDENamePrefix() << " " 
          << setw(4) << detElement->GetId();
     
      WriteTransform(out, transform);
    }
    out << endl;	  	   	
  }     
}

//______________________________________________________________________________
TString AliMUONGeometryTransformer::GetModuleSymName(Int_t moduleId) const
{
/// Return the module symbolic name (use for alignment)

  const AliMUONGeometryModuleTransformer* kTransformer 
    = GetModuleTransformer(moduleId);
  if ( ! kTransformer ) {
    AliErrorStream() << "Module " << moduleId << " not found." << endl; 
    return "";
  }   
  
  return "/" + fDetectorName + "/" + kTransformer->GetModuleName();
}  

//______________________________________________________________________________
TString AliMUONGeometryTransformer::GetDESymName(Int_t detElemId) const
{
/// Return the detection element symbolic name (used for alignment)

  const AliMUONGeometryDetElement* kDetElement 
    = GetDetElement(detElemId);
  if ( ! kDetElement ) {
    AliErrorStream() << "Det element " << detElemId << " not found." << endl; 
    return "";
  }   
  
  // Module Id
  Int_t moduleId = AliMpDEManager::GetGeomModuleId(detElemId);

  return GetModuleSymName(moduleId) + "/" + kDetElement->GetDEName();
}  

//
// public functions
//

//______________________________________________________________________________
Bool_t  
AliMUONGeometryTransformer::ReadGeometryData(
                                const TString& volPathFileName,
                                const TString& transformFileName)
{
/// Read geometry data from given files;
/// if transformFileName has ".root" extension, the transformations
/// are loaded from root geometry file, otherwise ASCII file
/// format is supposed

  Bool_t result1 = ReadVolPaths(volPathFileName);

  // Get file extension
  std::string fileName = transformFileName.Data();
  std::string rootExt = fileName.substr(fileName.size()-5, fileName.size());
  Bool_t result2;
  if ( rootExt != ".root" ) 
    result2 = ReadTransformations(transformFileName);
  else   
    result2 = ReadTransformations2(transformFileName);
  
  return result1 && result2;
}  

//______________________________________________________________________________
Bool_t  
AliMUONGeometryTransformer::ReadGeometryData(
                                const TString& volPathFileName,
                                TGeoManager* tgeoManager)
{
/// Load geometry data from root geometry using defined
/// volume paths from file

  Bool_t result1 = ReadVolPaths(volPathFileName);

  Bool_t result2 = LoadTransforms(tgeoManager);
  
  return result1 && result2;
}  

//______________________________________________________________________________
Bool_t  
AliMUONGeometryTransformer::WriteGeometryData(
                                 const TString& volPathFileName,
                                 const TString& transformFileName,
				 const TString& misalignFileName) const
{
/// Write geometry data into given files

  Bool_t result1 = WriteVolumePaths(volPathFileName);
  Bool_t result2 = WriteTransformations(transformFileName);
  
  Bool_t result3 = true;
  if ( misalignFileName != "" )
    result3 = WriteMisAlignmentData(misalignFileName);
  
  return result1 && result2 && result3;
}
				 
//______________________________________________________________________________
Bool_t  
AliMUONGeometryTransformer::WriteVolumePaths(const TString& fileName) const
{
/// Write volume paths for modules and detection element volumes into a file.
/// Return true, if writing finished correctly.

  // No writing
  // if builder is not associated with any geometry module
  if (fModuleTransformers->GetEntriesFast() == 0) return false;

  // File path
  TString filePath = gSystem->Getenv("ALICE_ROOT");
  filePath += "/MUON/data/";
  filePath += fileName;
  
  // Open output file
  ofstream out(filePath, ios::out);
  if (!out) {
    cerr << filePath << endl;	
    AliError("File not found.");
    return false;
  }
#if !defined (__DECCXX)
  out.setf(std::ios::fixed);
#endif
  WriteModuleVolPaths(out);
  WriteDetElemVolPaths(out);
  
  return true;
}  

//______________________________________________________________________________
Bool_t  
AliMUONGeometryTransformer::WriteTransformations(const TString& fileName) const
{
/// Write transformations into a file.
/// Return true, if writing finished correctly.

  // No writing
  // if builder is not associated with any geometry module
  if (fModuleTransformers->GetEntriesFast() == 0) return false;

  // File path
  TString filePath = gSystem->Getenv("ALICE_ROOT");
  filePath += "/MUON/data/";
  filePath += fileName;
  
  // Open output file
  ofstream out(filePath, ios::out);
  if (!out) {
    cerr << filePath << endl;	
    AliError("File not found.");
    return false;
  }
#if !defined (__DECCXX)
  out.setf(std::ios::fixed);
#endif
  WriteModuleTransforms(out);
  WriteDetElemTransforms(out);
  
  return true;
}  

//______________________________________________________________________________
Bool_t  
AliMUONGeometryTransformer::WriteMisAlignmentData(const TString& fileName) const
{
/// Write misalignment data into a file
/// Return true, if writing finished correctly.

  // No writing
  // if builder is not associated with any geometry module
  if ( fModuleTransformers->GetEntriesFast() == 0 ) {
    AliWarningStream() << "No geometry modules defined." << endl;
    return false;
  }  
  
  // No writing
  // if builder has no mis-alignment data
  if ( ! fMisAlignArray ) {
    AliWarningStream() << "No mis-alignment data defined." << endl;
    return false;
  }  

  // File path
  TString filePath = gSystem->Getenv("ALICE_ROOT");
  filePath += "/MUON/data/";
  filePath += fileName;
  
  // Write mis-alignment data in the root file
  TFile file(fileName.Data(), "RECREATE");
  fMisAlignArray->Write();
  file.Close();
  
  return true;
}  

//_____________________________________________________________________________
void AliMUONGeometryTransformer::AddModuleTransformer(
                          AliMUONGeometryModuleTransformer* moduleTransformer)
{
/// Add the module transformer to the array

  // Expand the size if not sufficient
  Int_t moduleId = moduleTransformer->GetModuleId();
  if (  moduleId >= fModuleTransformers->GetSize() )
    fModuleTransformers->Expand(moduleId+1);

  fModuleTransformers->AddAt(moduleTransformer, moduleId);
}

//_____________________________________________________________________________
void  AliMUONGeometryTransformer::AddMisAlignModule(Int_t moduleId, 
                                              const TGeoHMatrix& matrix)
{
/// Build AliAlignObjMatrix with module ID, its volumePath
/// and the given delta transformation matrix					      

  if ( ! fMisAlignArray )
    fMisAlignArray = new TClonesArray("AliAlignObjMatrix", 200);
    
  const AliMUONGeometryModuleTransformer* kTransformer 
    = GetModuleTransformer(moduleId);
  if ( ! kTransformer ) {
    AliErrorStream() << "Module " << moduleId << " not found." << endl; 
    return;
  }   
  
  // Get unique align object ID
  Int_t volId = AliAlignObj::LayerToVolUID(AliAlignObj::kMUON, moduleId); 

  // Create mis align matrix
  TClonesArray& refArray =*fMisAlignArray;
  Int_t pos = fMisAlignArray->GetEntriesFast();
  new (refArray[pos]) AliAlignObjMatrix(GetModuleSymName(moduleId), volId, 
					const_cast<TGeoHMatrix&>(matrix),kTRUE);
}

//_____________________________________________________________________________
void  AliMUONGeometryTransformer::AddMisAlignDetElement(Int_t detElemId, 
                                              const TGeoHMatrix& matrix)
{
/// Build AliAlignObjMatrix with detection element ID, its volumePath
/// and the given delta transformation matrix					      

  if ( ! fMisAlignArray )
    fMisAlignArray = new TClonesArray("AliAlignObjMatrix", 200);

  const AliMUONGeometryDetElement* kDetElement 
    = GetDetElement(detElemId);

  if ( ! kDetElement ) {
    AliErrorStream() << "Det element " << detElemId << " not found." << endl; 
    return;
  }   
  
  // Get unique align object ID
  Int_t volId = AliAlignObj::LayerToVolUID(AliAlignObj::kMUON, detElemId); 

  // Create mis align matrix
  TClonesArray& refArray =*fMisAlignArray;
  Int_t pos = fMisAlignArray->GetEntriesFast();
  new(refArray[pos]) AliAlignObjMatrix(GetDESymName(detElemId), volId, 
				       const_cast<TGeoHMatrix&>(matrix),kTRUE);
}

//_____________________________________________________________________________
void AliMUONGeometryTransformer::AddAlignableVolumes() const
{
/// Set symbolic names to alignable objects to TGeo

  if ( ! gGeoManager ) {
    AliWarning("TGeoManager not defined.");
    return;
  }  

  // Modules 
  for (Int_t i=0; i<fModuleTransformers->GetEntriesFast(); i++) {
    AliMUONGeometryModuleTransformer* module 
      = (AliMUONGeometryModuleTransformer*)fModuleTransformers->At(i);

    // Set module symbolic name
    gGeoManager->SetAlignableEntry(GetModuleSymName(module->GetModuleId()), 
                                   module->GetVolumePath());
    //cout << "Module sym name: " << GetModuleSymName(module->GetModuleId()) 
    //     << "  volPath: " << module->GetVolumePath() << endl;

    // Detection elements
    AliMpExMap* detElements = module->GetDetElementStore();    

    for (Int_t j=0; j<detElements->GetSize(); j++) {
      AliMUONGeometryDetElement* detElement
        = (AliMUONGeometryDetElement*)detElements->GetObject(j);
	
      // Set detection element symbolic name
      gGeoManager->SetAlignableEntry(GetDESymName(detElement->GetId()), 
                                     detElement->GetVolumePath());
      //cout << "DE name: " << GetDESymName(detElement->GetId()) 
      //     << "  volPath: " << detElement->GetVolumePath() << endl;
    }  
  }     
}  	     
    
//_____________________________________________________________________________
TClonesArray* AliMUONGeometryTransformer::CreateZeroAlignmentData() const
{
/// Create array with zero alignment data
			       
  // Create array for zero-alignment objects
  TClonesArray* array = new TClonesArray("AliAlignObjMatrix", 200);
  TClonesArray& refArray =*array;
  array->SetOwner(true);

  // Identity matrix
  TGeoHMatrix matrix;

  // Modules 
  for (Int_t i=0; i<fModuleTransformers->GetEntriesFast(); i++) {
    AliMUONGeometryModuleTransformer* module 
      = (AliMUONGeometryModuleTransformer*)fModuleTransformers->At(i);

    Int_t moduleId = module->GetModuleId();
  
    // Align object ID
    Int_t volId = AliAlignObj::LayerToVolUID(AliAlignObj::kMUON, moduleId); 

    // Create mis align matrix
    Int_t pos = array->GetEntriesFast();
    new (refArray[pos]) AliAlignObjMatrix(GetModuleSymName(moduleId), volId, matrix, kTRUE);
  }     

  // Detection elements
  for (Int_t i=0; i<fModuleTransformers->GetEntriesFast(); i++) {
    AliMUONGeometryModuleTransformer* moduleTransformer 
      = (AliMUONGeometryModuleTransformer*)fModuleTransformers->At(i);
    AliMpExMap* detElements = moduleTransformer->GetDetElementStore();    

    for (Int_t j=0; j<detElements->GetSize(); j++) {
      AliMUONGeometryDetElement* detElement
        = (AliMUONGeometryDetElement*)detElements->GetObject(j);
	
      Int_t detElemId = detElement->GetId();
  
      // Align object ID
      Int_t volId = AliAlignObj::LayerToVolUID(AliAlignObj::kMUON, detElemId); 

      // Create mis align matrix
      Int_t pos = array->GetEntriesFast();
      new (refArray[pos]) AliAlignObjMatrix(GetDESymName(detElemId), volId, matrix, kTRUE);
    }
  }
  
  return array;
}       

//_____________________________________________________________________________
void AliMUONGeometryTransformer::Global2Local(Int_t detElemId,
                 Float_t xg, Float_t yg, Float_t zg, 
                 Float_t& xl, Float_t& yl, Float_t& zl) const
{
/// Transform point from the global reference frame (ALIC)
/// to the local reference frame of the detection element specified
/// by detElemId.

  const AliMUONGeometryModuleTransformer* kTransformer 
    = GetModuleTransformerByDEId(detElemId);
  
  if (kTransformer) 
    kTransformer->Global2Local(detElemId, xg, yg, zg, xl, yl, zl);
}   
		 
//_____________________________________________________________________________
void AliMUONGeometryTransformer::Global2Local(Int_t detElemId,
                 Double_t xg, Double_t yg, Double_t zg, 
                 Double_t& xl, Double_t& yl, Double_t& zl) const
{
/// Transform point from the global reference frame (ALIC)
/// to the local reference frame of the detection element specified
/// by detElemId.

  const AliMUONGeometryModuleTransformer* kTransformer 
    = GetModuleTransformerByDEId(detElemId);
  
  if (kTransformer) 
    kTransformer->Global2Local(detElemId, xg, yg, zg, xl, yl, zl);
}   

//_____________________________________________________________________________
void AliMUONGeometryTransformer::Local2Global(Int_t detElemId,
                 Float_t xl, Float_t yl, Float_t zl, 
                 Float_t& xg, Float_t& yg, Float_t& zg) const
{		 
/// Transform point from the local reference frame of the detection element 
/// specified by detElemId to the global reference frame (ALIC).

  const AliMUONGeometryModuleTransformer* kTransformer 
    = GetModuleTransformerByDEId(detElemId);
    
  if (kTransformer) 
    kTransformer->Local2Global(detElemId, xl, yl, zl, xg, yg, zg);
}   

//_____________________________________________________________________________
void AliMUONGeometryTransformer::Local2Global(Int_t detElemId,
                 Double_t xl, Double_t yl, Double_t zl, 
                 Double_t& xg, Double_t& yg, Double_t& zg) const
{		 
/// Transform point from the local reference frame of the detection element 
/// specified by detElemId to the global reference frame (ALIC).

  const AliMUONGeometryModuleTransformer* kTransformer 
    = GetModuleTransformerByDEId(detElemId);
    
  if (kTransformer) 
    kTransformer->Local2Global(detElemId, xl, yl, zl, xg, yg, zg);
}   

//_____________________________________________________________________________
const AliMUONGeometryModuleTransformer* 
AliMUONGeometryTransformer::GetModuleTransformer(Int_t index, Bool_t warn) const
{
/// Return the geometry module transformer specified by index

  return GetModuleTransformerNonConst(index, warn);
}    

//_____________________________________________________________________________
const AliMUONGeometryModuleTransformer* 
AliMUONGeometryTransformer::GetModuleTransformerByDEId(Int_t detElemId, 
                                                       Bool_t warn) const
{
/// Return the geometry module transformer specified by detection element ID

  // Get module index
  Int_t index = AliMpDEManager::GetGeomModuleId(detElemId);

  return GetModuleTransformer(index, warn);
}    

//_____________________________________________________________________________
const AliMUONGeometryDetElement* 
AliMUONGeometryTransformer::GetDetElement(Int_t detElemId, Bool_t warn) const
{
/// Return detection element with given detElemId			       

  const AliMUONGeometryModuleTransformer* kTransformer 
    = GetModuleTransformerByDEId(detElemId, warn);
    
  if (!kTransformer) return 0;
    
  return kTransformer->GetDetElement(detElemId, warn); 
}

//_____________________________________________________________________________
Bool_t  AliMUONGeometryTransformer::HasDE(Int_t detElemId) const
{
/// Return true if detection element with given detElemId is defined

  const AliMUONGeometryModuleTransformer* kTransformer 
    = GetModuleTransformerByDEId(detElemId, false);
    
  if (!kTransformer) return false;
    
  return ( kTransformer->GetDetElement(detElemId, false) != 0 );
}  
    

