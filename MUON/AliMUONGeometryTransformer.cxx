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
// Class AliMUONGeometryTransformer
// ----------------------------
// Top container class for geometry transformations
//
// Author: Ivana Hrivnacova, IPN Orsay

#include <sstream>

#include <Riostream.h>
#include <TSystem.h>

#include "AliLog.h"

#include "AliMUONGeometryTransformer.h"
#include "AliMUONGeometryModuleTransformer.h"
#include "AliMUONGeometryDetElement.h"
#include "AliMUONGeometryStore.h"
#include "AliMUONGeometryBuilder.h"


ClassImp(AliMUONGeometryTransformer)
 
//______________________________________________________________________________
AliMUONGeometryTransformer::AliMUONGeometryTransformer(Bool_t isOwner)
  : TObject(),
    fModuleTransformers(0)
{
/// Standard constructor

  // Create array for geometry modules
  fModuleTransformers = new TObjArray();
  fModuleTransformers->SetOwner(isOwner);
}

//______________________________________________________________________________
AliMUONGeometryTransformer::AliMUONGeometryTransformer() 
  : TObject(),
    fModuleTransformers(0)
{
/// Default constructor
} 

//______________________________________________________________________________
AliMUONGeometryTransformer::AliMUONGeometryTransformer(
                                   const AliMUONGeometryTransformer& right) 
  : TObject(right) 
{  
/// Copy constructor (not implemented)

  AliFatal("Copy constructor not provided.");
}

//______________________________________________________________________________
AliMUONGeometryTransformer::~AliMUONGeometryTransformer()
{
/// Destructor

  delete fModuleTransformers;
}

//______________________________________________________________________________
AliMUONGeometryTransformer& 
AliMUONGeometryTransformer::operator=(const AliMUONGeometryTransformer& right)
{
/// Assignement operator (not implemented)

  // check assignement to self
  if (this == &right) return *this;

  AliFatal("Assignement operator not provided.");
    
  return *this;  
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
TString  AliMUONGeometryTransformer::ComposePath(const TString& volName, 
                                                 Int_t copyNo) const
{
// Compose path from given volName and copyNo
// ---

  TString path(volName);
  path += ".";
  path += copyNo;
  
  return path;
}  

//______________________________________________________________________________
TGeoHMatrix AliMUONGeometryTransformer::GetTransform(
                  Double_t x, Double_t y, Double_t z,
		  Double_t a1, Double_t a2, Double_t a3, 
 		  Double_t a4, Double_t a5, Double_t a6) const
{		  
// Builds the transformation from the given parameters
// ---

  // Compose transform
  return TGeoCombiTrans(TGeoTranslation(x, y, z), 
                        TGeoRotation("rot", a1, a2, a3, a4, a5, a6));
}


//______________________________________________________________________________
void AliMUONGeometryTransformer::FillData(Int_t moduleId,
                  Double_t x, Double_t y, Double_t z,
		  Double_t a1, Double_t a2, Double_t a3,
 		  Double_t a4, Double_t a5, Double_t a6) 
{
// Fill the transformation of the module.
// ---

  // Get/Create geometry module parametrisation
  moduleId--;
      // Modules numbers in the file are starting from 1

  AliMUONGeometryModuleTransformer* moduleTransformer
    = GetModuleTransformerNonConst(moduleId, false);

  if ( !moduleTransformer) {
    moduleTransformer = new AliMUONGeometryModuleTransformer(moduleId);
    AddModuleTransformer(moduleTransformer);
  }  
      
  // Build the transformation from the parameters
  TGeoHMatrix transform 
    = GetTransform(x, y, z, a1, a2, a3, a4, a5, a6);
      
  moduleTransformer->SetTransformation(transform);

}		   
  
//______________________________________________________________________________
void AliMUONGeometryTransformer::FillData(
                  Int_t detElemId, const TString& volName, Int_t copyNo,
                  Double_t x, Double_t y, Double_t z,
		  Double_t a1, Double_t a2, Double_t a3,
 		  Double_t a4, Double_t a5, Double_t a6) 
{
// Fill the transformation of the detection element.
// ---

  // Module Id
  Int_t moduleId = AliMUONGeometryStore::GetModuleId(detElemId);

  // Compose path
  TString path = ComposePath(volName, copyNo);
  
  // Build the transformation from the parameters
  TGeoHMatrix localTransform 
    = GetTransform(x, y, z, a1, a2, a3, a4, a5, a6);
   
  // Get detection element store
  AliMUONGeometryStore* detElements = 
    GetModuleTransformer(moduleId)->GetDetElementStore();     

  // Add detection element
  AliMUONGeometryDetElement* detElement
    = new AliMUONGeometryDetElement(detElemId, path, localTransform);
  detElements->Add(detElemId, detElement);
  
  // Compute global transformation
  const AliMUONGeometryModuleTransformer* kModuleTransformer
    = GetModuleTransformer(moduleId);
  if ( ! kModuleTransformer ) {
    AliFatal(Form("Module transformation not defined, detElemId %d",
                  detElemId));
  }  

  TGeoHMatrix globalTransform 
    = AliMUONGeometryBuilder::Multiply( 
                                  *kModuleTransformer->GetTransformation(),
				  localTransform );
  detElement->SetGlobalTransformation(globalTransform);
}		   
  
//______________________________________________________________________________
TString  AliMUONGeometryTransformer::ReadData1(ifstream& in)
{
// Reads and fills modules transformations from a file
// Returns true, if reading finished correctly.
// ---

  TString key("CH");
  while ( key == TString("CH") ) {
    Int_t id, n;
    Double_t  x, y, z;
    Double_t  a1, a2, a3, a4, a5, a6;
    TString dummy;
  
    in >> id;
    in >> n;
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

    //cout << "id="     << id << "  "
    // 	 << "position= " << x << ", " << y << ", " << z << "  "
    //	 << "rotation= " << a1 << ", " << a2 << ", " << a3  << ", "
    //	                 << a4 << ", " << a5 << ", " << a6 
    //	 << endl;   
	 
    // Fill data
    FillData(id, x, y, z, a1, a2, a3, a4, a5, a6);
    
    // Go to next line
    in >> key;
  }
  
  return key;   	 
}

//______________________________________________________________________________
TString  AliMUONGeometryTransformer::ReadData2(ifstream& in)
{
// Reads detection elements transformations from a file
// Returns true, if reading finished correctly.
// ---

  TString key("DE");
  while ( key == TString("DE") ) {

    // Input data
    Int_t detElemId;
    TString   volumeName;
    Int_t     copyNo;
    Double_t  x, y, z;
    Double_t  a1, a2, a3, a4, a5, a6;
    TString dummy;
  
    in >> detElemId;
    in >> volumeName;
    in >> copyNo;
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
    //     << "volume=" << volumeName << "  "
    //     << "copyNo=" << copyNo << "  "
    //     << "position= " << x << ", " << y << ", " << z << "  "
    //     << "rotation= " << a1 << ", " << a2 << ", " << a3  << ", "
    //	                 << a4 << ", " << a5 << ", " << a6 
    //     << endl;   
	 
    // Fill data
    FillData(detElemId, volumeName, copyNo, x, y, z, a1, a2, a3, a4, a5, a6); 	 
    
    // Go to next line
    in >> key;
  } 
  
  return key;
}

//______________________________________________________________________________
void AliMUONGeometryTransformer::WriteTransform(ofstream& out,
                                   const TGeoCombiTrans* transform) const
{
// Writes the transformations 
// ---

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
void AliMUONGeometryTransformer::WriteData1(ofstream& out) const
{
// Writes modules transformations
// ---

  for (Int_t i=0; i<fModuleTransformers->GetEntriesFast(); i++) {
    AliMUONGeometryModuleTransformer* moduleTransformer 
      = (AliMUONGeometryModuleTransformer*)fModuleTransformers->At(i);
    const TGeoCombiTrans* transform 
      = moduleTransformer->GetTransformation();    

    out << "CH " 
        << setw(4) << moduleTransformer->GetModuleId() + 1 << "  "
        << setw(4) << moduleTransformer->GetDetElementStore()->GetNofEntries() << "  ";
    
    WriteTransform(out, transform);
  }
  out << endl;
}

//______________________________________________________________________________
void AliMUONGeometryTransformer::WriteData2(ofstream& out) const
{
// Writes detection elements (envelopes) transformations
// ---


  for (Int_t i=0; i<fModuleTransformers->GetEntriesFast(); i++) {
    AliMUONGeometryModuleTransformer* moduleTransformer 
      = (AliMUONGeometryModuleTransformer*)fModuleTransformers->At(i);
    AliMUONGeometryStore* detElements 
      = moduleTransformer->GetDetElementStore();    

    for (Int_t j=0; j<detElements->GetNofEntries(); j++) {
      AliMUONGeometryDetElement* detElement
        = (AliMUONGeometryDetElement*)detElements->GetEntry(j);
      const TGeoCombiTrans* transform 
        = detElement->GetLocalTransformation(); 
	
      // Get volume name & copy number from aligned volume
      string volCopyNo = detElement->GetAlignedVolume().Data();
      std::string::size_type first = volCopyNo.find('.');
      std::string volName = volCopyNo.substr(0, first);
      std::string copyNoStr = volCopyNo.substr(first+1, volCopyNo.length());
      std::istringstream in(copyNoStr);
      Int_t copyNo;
      in >> copyNo;
      
      // Write data on out
      out << "DE " 
          << setw(4) << detElement->GetId() << "    " 
          << volName << " "
	  << setw(4) << copyNo;
     
      WriteTransform(out, transform);
    }
    out << endl;	  	   	
  }     
}

//
// public functions
//

//______________________________________________________________________________
Bool_t  
AliMUONGeometryTransformer::ReadTransformations(const TString& fileName)
{
// Reads transformations from a file
// Returns true, if reading finished correctly.
// ---

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
    if (key == TString("CH")) 
      key = ReadData1(in);
    else if (key == TString("DE"))
      key = ReadData2(in);
    else {
      AliFatal(Form("%s key not recognized",  key.Data()));
      return false;
    }
  }     

  return true;
}

//______________________________________________________________________________
Bool_t  
AliMUONGeometryTransformer::WriteTransformations(const TString& fileName) const
{
// Writes transformations into a file
// Returns true, if writing finished correctly.
// ---

  // No writing
  // if builder is not associated with any geometry module
  if (fModuleTransformers->GetEntriesFast() == 0) return false;

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
  WriteData1(out);
  WriteData2(out);
  
  return true;
}  

//_____________________________________________________________________________
void AliMUONGeometryTransformer::AddModuleTransformer(
                          AliMUONGeometryModuleTransformer* moduleTransformer)
{
/// Add the geometrymodule to the array

  fModuleTransformers->AddAt(moduleTransformer, 
                             moduleTransformer->GetModuleId());
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
/// Return the geometry module specified by index

  return GetModuleTransformerNonConst(index, warn);
}    

//_____________________________________________________________________________
const AliMUONGeometryModuleTransformer* 
AliMUONGeometryTransformer::GetModuleTransformerByDEId(Int_t detElemId, 
                                                       Bool_t warn) const
{
/// Return the geometry module specified by index

  // Get module index
  Int_t index = AliMUONGeometryStore::GetModuleId(detElemId);

  return GetModuleTransformer(index, warn);
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
    

