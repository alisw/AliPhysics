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
// Abstract base class for geometry construction per chamber(s).
// Author: Ivana Hrivnacova, IPN Orsay
// 23/01/2004

#include <Riostream.h>
#include <TObjArray.h>
#include <TSystem.h>
#include <TGeoMatrix.h>
#include <TVirtualMC.h>

#include "AliMUONVGeometryBuilder.h"
#include "AliMUONChamber.h"
#include "AliMUONChamberGeometry.h"
#include "AliMUONGeometryTransformStore.h"
#include "AliMUONGeometrySVMap.h"
#include "AliMUONGeometryEnvelopeStore.h"
#include "AliMUONGeometryEnvelope.h"
#include "AliMUONGeometryConstituent.h"
#include "AliMUONConstants.h"

ClassImp(AliMUONVGeometryBuilder)

const TString AliMUONVGeometryBuilder::fgkTransformFileNamePrefix = "transform_";
const TString AliMUONVGeometryBuilder::fgkSVMapFileNamePrefix = "svmap_";
const TString AliMUONVGeometryBuilder::fgkOutFileNameSuffix = ".out";

//______________________________________________________________________________
AliMUONVGeometryBuilder::AliMUONVGeometryBuilder(const TString& fileName,
                                AliMUONChamber* ch1, AliMUONChamber* ch2,
                                AliMUONChamber* ch3, AliMUONChamber* ch4,
                                AliMUONChamber* ch5, AliMUONChamber* ch6)
 : TObject(),
   fTransformFileName(fgkTransformFileNamePrefix+fileName),
   fSVMapFileName(fgkSVMapFileNamePrefix+fileName),
   fChambers(0)
 {
// Standard constructor

  // Create the chambers array
  fChambers = new TObjArray();
  
  if (ch1) fChambers->Add(ch1);
  if (ch2) fChambers->Add(ch2);
  if (ch3) fChambers->Add(ch3);
  if (ch4) fChambers->Add(ch4);
  if (ch5) fChambers->Add(ch5);
  if (ch6) fChambers->Add(ch6);
}


//______________________________________________________________________________
AliMUONVGeometryBuilder::AliMUONVGeometryBuilder()
 : TObject(),
   fTransformFileName(),
   fSVMapFileName(),
   fChambers(0)
{
// Default constructor
}


//______________________________________________________________________________
AliMUONVGeometryBuilder::AliMUONVGeometryBuilder(const AliMUONVGeometryBuilder& rhs)
  : TObject(rhs)
{
// Protected copy constructor

  Fatal("Copy constructor", 
        "Copy constructor is not implemented.");
}

//______________________________________________________________________________
AliMUONVGeometryBuilder::~AliMUONVGeometryBuilder() {
//
  if (fChambers) {
    fChambers->Clear(); // Sets pointers to 0 sinche it is not the owner
    delete fChambers;
  }
}

//______________________________________________________________________________
AliMUONVGeometryBuilder& 
AliMUONVGeometryBuilder::operator = (const AliMUONVGeometryBuilder& rhs) 
{
// Protected assignement operator

  // check assignement to self
  if (this == &rhs) return *this;

  Fatal("operator=", 
        "Assignment operator is not implemented.");
    
  return *this;  
}

//
// private methods
//

//______________________________________________________________________________
 TString  AliMUONVGeometryBuilder::ComposePath(const TString& volName, 
                                               Int_t copyNo) const
{
// Compose path from given volName and copyNo
// ---

  TString path(volName);
  path += ".";
  path += copyNo;
  
  return path;
}  

/*
//______________________________________________________________________________
void AliMUONVGeometryBuilder::MapSV(const TString& path0, 
                                    const TString& volName, Int_t detElemId) const
{
// Update the path with all daughters volumes recursively
// and map it to the detection element Id if it is a sensitive volume
// ---

  Int_t nofDaughters = gMC->NofVolDaughters(volName);
  if (nofDaughters == 0) {

    // Get the name of the last volume in the path
    Ssiz_t npos1 = path0.Last('/')+1; 
    Ssiz_t npos2 = path0.Last('.');
    TString volName(path0(npos1, npos2-npos1));  
    
    // Check if it is sensitive volume
    Int_t chamberId = AliMUONConstants::GetChamberId(detElemId);
    AliMUONChamberGeometry* geometry = GetChamber(chamberId)->GetGeometry();
    if (geometry->IsSensitiveVolume(volName)) {
      //cout << ".. adding to the map  " 
      //     <<  path0 << "  "  << detElemId << endl;
      FillData(path0, detElemId); 
    }  
    return; 
  }  

  for (Int_t i=0; i<nofDaughters; i++) {
    Int_t copyNo = gMC->VolDaughterCopyNo(volName, i);
    TString newName =  gMC->VolDaughterName(volName, i);
            
    TString path = path0;
    path += "/";
    path += ComposePath(newName, copyNo);

    MapSV(path, newName, detElemId);
  }
}     
*/

//______________________________________________________________________________
void AliMUONVGeometryBuilder::MapSV(const TString& /*path0*/, 
                                    const TString& /*volName*/, 
				    Int_t /*detElemId*/) const
{
// Update the path with all daughters volumes recursively
// and map it to the detection element Id if it is a sensitive volume
// ---

  Warning("MapSV", "Not yet available");
}     

//______________________________________________________________________________
void AliMUONVGeometryBuilder::FillData(Int_t chamberId,
                  Double_t x, Double_t y, Double_t z,
		  Double_t a1, Double_t a2, Double_t a3,
 		  Double_t a4, Double_t a5, Double_t a6) const 
{
// Fill the transformation of the chamber.
// ---

  chamberId--;
      // Chambers numbers in the file are starting from 1

  GetChamber(chamberId)->GetGeometry()
    ->SetTranslation(TGeoTranslation(x, y, z));
  GetChamber(chamberId)->GetGeometry()
    ->SetRotation(TGeoRotation("rot", a1, a2, a3, a4, a5, a6));
}		   
  
//______________________________________________________________________________
void AliMUONVGeometryBuilder::FillData(
                  Int_t detElemId, const TString& volName, Int_t copyNo,
                  Double_t x, Double_t y, Double_t z,
		  Double_t a1, Double_t a2, Double_t a3,
 		  Double_t a4, Double_t a5, Double_t a6) const 
{
// Fill the transformation of the detection element.
// ---

  // Chamber Id
  Int_t chamberId = AliMUONConstants::GetChamberId(detElemId);

  // Get chamber transformations
  AliMUONGeometryTransformStore* transforms 
    = GetTransforms(chamberId);     

  // Compose path
  TString path = ComposePath(volName, copyNo);
  
  // Compose matrix
  TGeoCombiTrans transform(path, x, y, z, 
                           new TGeoRotation(path, a1, a2, a3, a4, a5, a6));
    
  // Add detection element transformation 
  transforms->Add(detElemId, path, transform); 
}		   
  
//______________________________________________________________________________
void AliMUONVGeometryBuilder::FillData(
                   const TString& sensVolumePath, Int_t detElemId) const
{
// Fill the mapping of the sensitive volume path to the detection element.
// ---

  // Chamber Id
  Int_t chamberId = AliMUONConstants::GetChamberId(detElemId);

  // Get chamber transformations
  AliMUONGeometrySVMap* svMap = GetSVMap(chamberId);     

  // Map the sensitive volume to detection element
  svMap->Add(sensVolumePath, detElemId); 
}		   
  
//______________________________________________________________________________
TString  AliMUONVGeometryBuilder::ReadData1(ifstream& in) const
{
// Reads and fills chambers transformations from a file
// Returns true, if reading finished correctly.
// ---

  TString key("CH");
  while ( key == TString("CH") ) {
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
TString  AliMUONVGeometryBuilder::ReadData2(ifstream& in) const
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
TString  AliMUONVGeometryBuilder::ReadData3(ifstream& in) const
{
// Reads detection elements transformations from a file
// Returns true, if reading finished correctly.
// ---

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
    FillData(volumePath, detElemId); 
     
    // Go to next line
    in >> key;
  } 
  
  return key;
}

//______________________________________________________________________________
void AliMUONVGeometryBuilder::WriteTransform(ofstream& out, 
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
  Double_t a1, a2, a3, a4, a5, a6;
  TGeoRotation* rotation = transform->GetRotation();
  rotation->GetAngles(a1, a2, a3, a4, a5, a6);
  out << setw(8) << setprecision(4) << a1 << "  " 
      << setw(8) << setprecision(4) << a2 << "  " 
      << setw(8) << setprecision(4) << a3 << "  " 
      << setw(8) << setprecision(4) << a4 << "  " 
      << setw(8) << setprecision(4) << a5 << "  " 
      << setw(8) << setprecision(4) << a6 << "  " << endl; 
}

//______________________________________________________________________________
void AliMUONVGeometryBuilder::WriteData1(ofstream& out) const
{
// Writes chamber transformations
// ---

  for (Int_t i=0; i<fChambers->GetEntriesFast(); i++) {
    AliMUONChamber* chamber 
      = (AliMUONChamber*)fChambers->At(i);
    AliMUONChamberGeometry* chamberGeometry
       = chamber->GetGeometry();
    const TGeoCombiTrans* transform 
      = chamberGeometry->GetTransformation();    

    out << "CH " 
        << setw(4) << chamber->GetId() + 1 << "  ";
    
    WriteTransform(out, transform);
  }
  out << endl;
}

//______________________________________________________________________________
void AliMUONVGeometryBuilder::WriteData2(ofstream& out) const
{
// Writes detection elements (envelopes) transformations
// ---


  for (Int_t i=0; i<fChambers->GetEntriesFast(); i++) {
    AliMUONChamber* chamber 
      = (AliMUONChamber*)fChambers->At(i);
    AliMUONChamberGeometry* chamberGeometry
       = chamber->GetGeometry();
    const TObjArray* envelopes 
      = chamberGeometry->GetEnvelopeStore()->GetEnvelopes();    

    for (Int_t j=0; j<envelopes->GetEntriesFast(); j++) {
      AliMUONGeometryEnvelope* envelope
        = (AliMUONGeometryEnvelope*)envelopes->At(j);
      const TGeoCombiTrans* transform 
        = envelope->GetTransformation(); 
      
      // skip envelope not corresponding to detection element
      if(envelope->GetUniqueID() == 0) continue;
       
      out << "DE " 
          << setw(4) << envelope->GetUniqueID() << "    " 
          << envelope->GetName() << " " 
	  << setw(4)<< envelope->GetCopyNo();
     
      WriteTransform(out, transform);
    }
    out << endl;	  	   	
  }     
}

//______________________________________________________________________________
void AliMUONVGeometryBuilder::WriteData3(ofstream& out) const
{
// Writes association of sensitive volumes and detection elements
// from the sensitive volume map
// ---

  for (Int_t i=0; i<fChambers->GetEntriesFast(); i++) {
    AliMUONChamber* chamber 
      = (AliMUONChamber*)fChambers->At(i);
    AliMUONChamberGeometry* chamberGeometry
       = chamber->GetGeometry();
    AliMUONGeometrySVMap* svMap
      = chamberGeometry->GetSVMap();

    svMap->WriteMap(out);
    out << endl;  
  }    
}

//
// protected methods
//

//______________________________________________________________________________
AliMUONChamber*  AliMUONVGeometryBuilder::GetChamber(Int_t chamberId) const
{
// Returns the chamber specified by chamberId
// ---

  for (Int_t i=0; i<fChambers->GetEntriesFast(); i++) {
    AliMUONChamber* chamber = (AliMUONChamber*)fChambers->At(i);
    if ( chamber->GetId() == chamberId) return chamber;
  }   
  
  return 0;
}  

//______________________________________________________________________________
AliMUONGeometryEnvelopeStore*  
AliMUONVGeometryBuilder::GetEnvelopes(Int_t chamberId) const
{
// Returns the envelope store of the chamber specified by chamberId
// ---

  AliMUONChamber* chamber = GetChamber(chamberId);
  
  if (!chamber) {
    Fatal("GetEnvelopes", "Chamber %d is not defined", chamberId); 
    return 0;
  }
  
  return chamber->GetGeometry()->GetEnvelopeStore();
}  

//______________________________________________________________________________
AliMUONGeometryTransformStore*  
AliMUONVGeometryBuilder::GetTransforms(Int_t chamberId) const
{
// Returns the transformation store of the chamber specified by chamberId
// ---

  AliMUONChamber* chamber = GetChamber(chamberId);
  
  if (!chamber) {
    Fatal("GetTransforms", "Chamber %d is not defined", chamberId); 
    return 0;
  }
  
  return chamber->GetGeometry()->GetTransformStore();
}  

//______________________________________________________________________________
AliMUONGeometrySVMap*  
AliMUONVGeometryBuilder::GetSVMap(Int_t chamberId) const
{
// Returns the transformation store of the chamber specified by chamberId
// ---

  AliMUONChamber* chamber = GetChamber(chamberId);
  
  if (!chamber) {
    Fatal("GetSVMap", "Chamber %d is not defined", chamberId); 
    return 0;
  }
  
  return chamber->GetGeometry()->GetSVMap();
}  

//
// public functions
//

//______________________________________________________________________________
void  AliMUONVGeometryBuilder:: FillTransformations() const
{
// Fills transformations store from defined geometry.
// ---

  for (Int_t i=0; i<fChambers->GetEntriesFast(); i++) {
    AliMUONChamber* chamber 
      = (AliMUONChamber*)fChambers->At(i);
    AliMUONChamberGeometry* chamberGeometry
       = chamber->GetGeometry();
    const TObjArray* envelopes 
      = chamberGeometry->GetEnvelopeStore()->GetEnvelopes();    
    AliMUONGeometryTransformStore* transforms 
      = chamberGeometry->GetTransformStore();     

    for (Int_t j=0; j<envelopes->GetEntriesFast(); j++) {
      AliMUONGeometryEnvelope* envelope
        = (AliMUONGeometryEnvelope*)envelopes->At(j);

      // skip envelope not corresponding to detection element
      if(envelope->GetUniqueID() == 0) continue;
       
      // Get envelope data 
      Int_t detElemId = envelope->GetUniqueID();	
      TString path = ComposePath(envelope->GetName(), 
                                 envelope->GetCopyNo());
      const TGeoCombiTrans* transform = envelope->GetTransformation(); 

      // Add detection element transformation 
      transforms->Add(detElemId, path, *transform); 
    }  
  }
}

//_____ _________________________________________________________________________
void  AliMUONVGeometryBuilder::RebuildSVMaps() const
{
// Clear the SV maps in memory and fill them from defined geometry.
// ---

  for (Int_t i=0; i<fChambers->GetEntriesFast(); i++) {
    AliMUONChamber* chamber 
      = (AliMUONChamber*)fChambers->At(i);
    AliMUONChamberGeometry* chamberGeometry
       = chamber->GetGeometry();
    
    // Clear the map   
    chamberGeometry->GetSVMap()->Clear();
     
    // Fill the map from geometry
    const TObjArray* envelopes 
      = chamberGeometry->GetEnvelopeStore()->GetEnvelopes();    

    for (Int_t j=0; j<envelopes->GetEntriesFast(); j++) {
      AliMUONGeometryEnvelope* envelope
        = (AliMUONGeometryEnvelope*)envelopes->At(j);

      // skip envelope not corresponding to detection element
      if(envelope->GetUniqueID() == 0) continue;
       
      TString path0("/ALIC.1");
      if (chamberGeometry->GetMotherVolume() != "ALIC") {
        path0 += "/";
	path0 += ComposePath(chamberGeometry->GetMotherVolume(), 1);
      }  
       
      if (!envelope->IsVirtual()) {
         TString path = path0;
         path += "/";
	 path += ComposePath(envelope->GetName(), envelope->GetCopyNo());
	 MapSV(path, envelope->GetName(), envelope->GetUniqueID());
      }
      else {	 
        for  (Int_t k=0; k<envelope->GetConstituents()->GetEntriesFast(); k++) {
          AliMUONGeometryConstituent* constituent
            = (AliMUONGeometryConstituent*)envelope->GetConstituents()->At(k);
         TString path = path0;
         path += "/";
	 path += ComposePath(constituent->GetName(), constituent->GetCopyNo());
	 MapSV(path, constituent->GetName(), envelope->GetUniqueID());
        }
      }
    }  
  } 	             
}

//______________________________________________________________________________
Bool_t  AliMUONVGeometryBuilder::ReadTransformations() const
{
// Reads transformations from a file
// Returns true, if reading finished correctly.
// ---

  // File path
  TString filePath = gSystem->Getenv("ALICE_ROOT");
  filePath += "/MUON/data/";
  filePath += fTransformFileName;
  
  // Open input file
  ifstream in(filePath, ios::in);
  if (!in) {
    cerr << filePath << endl;	
    Fatal("ReadTransformations", "File not found.");
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
      Fatal("ReadTransformations", "%s key not recognized",  key.Data());
      return false;
    }
  }     

  return true;
}

//______________________________________________________________________________
Bool_t  AliMUONVGeometryBuilder::ReadSVMap() const
{
// Reads the sensitive volume from a file
// Returns true, if reading finished correctly.
// ---

  // File path
  TString filePath = gSystem->Getenv("ALICE_ROOT");
  filePath += "/MUON/data/";
  filePath += fSVMapFileName;
  
  // Open input file
  ifstream in(filePath, ios::in);
  if (!in) {
    cerr << filePath << endl;	
    Fatal("ReadSVMap", "File not found.");
    return false;
  }

  TString key;
  in >> key;
  while ( !in.eof() ) {
    if (key == TString("SV")) 
      key = ReadData3(in);
    else {
      Fatal("ReadSVMap", "%s key not recognized",  key.Data());
      return false;
    }
  }     

  return true;
}

//______________________________________________________________________________
Bool_t  AliMUONVGeometryBuilder::WriteTransformations() const
{
// Writes transformations into a file
// Returns true, if writing finished correctly.
// ---

  // File path
  TString filePath = gSystem->Getenv("ALICE_ROOT");
  filePath += "/MUON/data/";
  filePath += fTransformFileName;
  filePath += fgkOutFileNameSuffix;
  
  // Open input file
  ofstream out(filePath, ios::out);
  if (!out) {
    cerr << filePath << endl;	
    Error("WriteTransformations", "File not found.");
    return false;
  }
#if !defined (__DECCXX)
  out.setf(std::ios::fixed);
#endif
  WriteData1(out);
  WriteData2(out);
  
  return true;
}  

//______________________________________________________________________________
Bool_t  AliMUONVGeometryBuilder::WriteSVMap(Bool_t rebuild) const
{
// Writes sensitive volume map into a file
// Returns true, if writing finished correctly.
// ---

  // File path
  TString filePath = gSystem->Getenv("ALICE_ROOT");
  filePath += "/MUON/data/";
  filePath += fSVMapFileName;
  filePath += fgkOutFileNameSuffix;
  
  // Open input file
  ofstream out(filePath, ios::out);
  if (!out) {
    cerr << filePath << endl;	
    Error("WriteTransformations", "File not found.");
    return false;
  }
#if !defined (__DECCXX)
  out.setf(std::ios::fixed);
#endif  
  if (rebuild)  RebuildSVMaps();

  WriteData3(out);
  
  return true;
}  

