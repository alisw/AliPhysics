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

//-----------------------------------------------------------------------------
// Class AliMUONGeometryTransformer
// ----------------------------
// Top container class for geometry transformations
// Author: Ivana Hrivnacova, IPN Orsay
//-----------------------------------------------------------------------------

#include "AliMUONGeometryTransformer.h"
#include "AliMUONGeometryModuleTransformer.h"
#include "AliMUONGeometryDetElement.h"
#include "AliMUONGeometryBuilder.h"

#include "AliMpDEManager.h"
#include "AliMpConstants.h"
#include "AliMpExMap.h"
#include "AliMpCDB.h"
#include "AliMpArea.h"
#include <float.h>
#include "AliMpVPadIterator.h"
#include "AliMpPad.h"
#include "AliMpDEIterator.h"
#include <TVector2.h>
#include "AliMpVSegmentation.h"
#include "AliMpSegmentation.h"
#include "AliMpExMapIterator.h"
#include "AliLog.h"
#include "AliAlignObjMatrix.h"
#include "AliAlignObj.h"

#include <Riostream.h>
#include <TSystem.h>
#include <TClonesArray.h>
#include <TGeoManager.h>
#include <TGeoPhysicalNode.h>
#include <TFile.h>
#include <TString.h>

#include <sstream>

/// \cond CLASSIMP
ClassImp(AliMUONGeometryTransformer)
/// \endcond

const TString  AliMUONGeometryTransformer::fgkDefaultDetectorName = "MUON";
 
//______________________________________________________________________________
AliMUONGeometryTransformer::AliMUONGeometryTransformer()

  : TObject(),
    fDetectorName(fgkDefaultDetectorName),
    fModuleTransformers(0),
    fMisAlignArray(0),
    fDEAreas(0x0)
{
/// Standard constructor

  // Create array for geometry modules
  fModuleTransformers = new TObjArray(100);
  fModuleTransformers->SetOwner(true);
}

//______________________________________________________________________________
AliMUONGeometryTransformer::AliMUONGeometryTransformer(TRootIOCtor* /*ioCtor*/) 
  : TObject(),
    fDetectorName(),
    fModuleTransformers(0),
    fMisAlignArray(0),
    fDEAreas(0x0)
{
/// Default constructor
} 

//______________________________________________________________________________
AliMUONGeometryTransformer::~AliMUONGeometryTransformer()
{
/// Destructor

  delete fModuleTransformers;
  delete fMisAlignArray;
  delete fDEAreas;
}

//
// private methods
//


//_____________________________________________________________________________
AliMpArea*
AliMUONGeometryTransformer::GetDEArea(Int_t detElemId) const
{
  /// Get area (in global coordinates) covered by a given detection element
  if (!fDEAreas)
  {
    CreateDEAreas();
  }
  return static_cast<AliMpArea*>(fDEAreas->GetValue(detElemId));
}

//_____________________________________________________________________________
void
AliMUONGeometryTransformer::CreateDEAreas() const
{
  /// Create DE areas
  
  fDEAreas = new AliMpExMap;
  
  AliMpDEIterator it;

  it.First();

  /// Generate the DE areas in global coordinates

  while ( !it.IsDone() )
  {
    Int_t detElemId = it.CurrentDEId();
    
    if ( !HasDE(detElemId) ) continue;
    
    const AliMpVSegmentation* seg = AliMpSegmentation::Instance()->GetMpSegmentation(detElemId,AliMp::kCath0);
    
    Double_t xg,yg,zg;
    
    AliMp::StationType stationType = AliMpDEManager::GetStationType(detElemId);
    
    Double_t xl(0.0), yl(0.0), zl(0.0);
    Double_t dx(seg->GetDimensionX());
    Double_t dy(seg->GetDimensionY());
    
    if ( stationType == AliMp::kStation12 ) 
    {
      Double_t xmin(FLT_MAX);
      Double_t xmax(-FLT_MAX);
      Double_t ymin(FLT_MAX);
      Double_t ymax(-FLT_MAX);
      
      for ( Int_t icathode = 0; icathode < 2; ++icathode ) 
      {
        const AliMpVSegmentation* cathode 
        = AliMpSegmentation::Instance()->GetMpSegmentation(detElemId,AliMp::GetCathodType(icathode));
        
        AliMpVPadIterator* itp = cathode->CreateIterator();
        
        itp->First();
        
        while ( !itp->IsDone() ) 
        {
          AliMpPad pad = itp->CurrentItem();
          AliMpArea a(pad.GetPositionX(),pad.GetPositionY(),
                      pad.GetDimensionX(), pad.GetDimensionY());
          xmin = TMath::Min(xmin,a.LeftBorder());
          xmax = TMath::Max(xmax,a.RightBorder());
          ymin = TMath::Min(ymin,a.DownBorder());
          ymax = TMath::Max(ymax,a.UpBorder());
          itp->Next();
        }
        
        delete itp;
      }
      
      xl = (xmin+xmax)/2.0;
      yl = (ymin+ymax)/2.0;
      dx = (xmax-xmin)/2.0;
      dy = (ymax-ymin)/2.0;
      
      Local2Global(detElemId,xl,yl,zl,xg,yg,zg);
    }
    else
    {
      Local2Global(detElemId,xl,yl,zl,xg,yg,zg);
    }
    
    fDEAreas->Add(detElemId,new AliMpArea(xg,yg,dx,dy));
    
    it.Next();
  }
}

//_____________________________________________________________________________
Bool_t AliMUONGeometryTransformer::LoadMapping() const
{
/// Load mapping from CDB

  if ( ! AliMpCDB::LoadMpSegmentation() ) 
  {
    AliFatal("Could not access mapping from OCDB !");
    return false;
  }
  
  return true;
}  

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
    TIter next(detElements->CreateIterator());
    AliMUONGeometryDetElement* detElement;
    while ( ( detElement = static_cast<AliMUONGeometryDetElement*>(next()) ) )
    {
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

  return "/" + fDetectorName + "/" 
             + AliMUONGeometryModuleTransformer::GetModuleName(moduleId);
}  

//______________________________________________________________________________
TString AliMUONGeometryTransformer::GetDESymName(Int_t detElemId) const
{
/// Return the detection element symbolic name (used for alignment)

  // Module Id
  Int_t moduleId = AliMpDEManager::GetGeomModuleId(detElemId);

  return GetModuleSymName(moduleId) + "/" 
         + AliMUONGeometryDetElement::GetDEName(detElemId);
}  

//
// public functions
//

//______________________________________________________________________________
Bool_t  
AliMUONGeometryTransformer::LoadTransformations()
{
/// Load transformations for defined modules and detection elements
/// using AliGeomManager

  if ( ! AliGeomManager::GetGeometry() ) {
    AliFatal("Geometry has to be laoded in AliGeomManager first.");
    return false;
  }   

  for (Int_t i=0; i<fModuleTransformers->GetEntriesFast(); i++) {
    AliMUONGeometryModuleTransformer* moduleTransformer 
      = (AliMUONGeometryModuleTransformer*)fModuleTransformers->At(i);

    // Module symbolic name
    TString symname = GetModuleSymName(moduleTransformer->GetModuleId());
    
    // Set matrix from physical node
    TGeoHMatrix* matrix = AliGeomManager::GetMatrix(symname);
    if ( ! matrix ) {
      AliErrorStream() << "Geometry module matrix not found." << endl;
      return false;
    }  
    moduleTransformer->SetTransformation(*matrix);
    
    // Loop over detection elements
    AliMpExMap* detElements = moduleTransformer->GetDetElementStore();    
    TIter next(detElements->CreateIterator());    
    AliMUONGeometryDetElement* detElement;
    
    while ( ( detElement = static_cast<AliMUONGeometryDetElement*>(next()) ) )
    {
      // Det element  symbolic name
      TString symnameDE = GetDESymName(detElement->GetId());
    
      // Set global matrix from physical node
      TGeoHMatrix* globalMatrix = AliGeomManager::GetMatrix(symnameDE);
      if ( ! globalMatrix ) {
        AliErrorStream() << "Detection element matrix not found." << endl;
        return false;
      }  
      detElement->SetGlobalTransformation(*globalMatrix, false);

      // Set local matrix
      TGeoHMatrix localMatrix = 
        AliMUONGeometryBuilder::Multiply(
	   (*matrix).Inverse(), (*globalMatrix) );
      detElement->SetLocalTransformation(localMatrix, false);
    }  
  } 
  return true;    
}  

//______________________________________________________________________________
Bool_t  
AliMUONGeometryTransformer::LoadGeometryData(const TString& fileName)
{
/// Read geometry data either from ASCII file with transformations or
/// from root geometry file (if fileName has ".root" extension)

  CreateModules();

  // Get file extension
  std::string fileName2 = fileName.Data();
  std::string rootExt = fileName2.substr(fileName2.size()-5, fileName2.size());
  
  if ( rootExt != ".root" ) 
    return ReadTransformations(fileName);
  else  { 
    // Load root geometry
    AliGeomManager::LoadGeometry(fileName.Data());
    return LoadTransformations();
  }  
}  

//______________________________________________________________________________
Bool_t  
AliMUONGeometryTransformer::LoadGeometryData()
{
/// Load geometry data from already loaded Root geometry using AliGeomManager

  if ( ! AliGeomManager::GetGeometry() ) {
    AliErrorStream() << "Geometry has not been loaded in AliGeomManager" << endl;
    return false;
  }    

  CreateModules();

  return LoadTransformations();
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
						    const TGeoHMatrix& matrix, Bool_t bGlobal)
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
  Int_t volId = AliGeomManager::LayerToVolUID(AliGeomManager::kMUON, moduleId); 

  // Create mis align matrix
  TClonesArray& refArray =*fMisAlignArray;
  Int_t pos = fMisAlignArray->GetEntriesFast();
  new (refArray[pos]) AliAlignObjMatrix(GetModuleSymName(moduleId), volId, 
					const_cast<TGeoHMatrix&>(matrix),bGlobal);
}

//_____________________________________________________________________________
void  AliMUONGeometryTransformer::AddMisAlignDetElement(Int_t detElemId, 
							const TGeoHMatrix& matrix, Bool_t bGlobal)
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
  Int_t volId = AliGeomManager::LayerToVolUID(AliGeomManager::kMUON, detElemId); 

  // Create mis align matrix
  TClonesArray& refArray =*fMisAlignArray;
  Int_t pos = fMisAlignArray->GetEntriesFast();
  new(refArray[pos]) AliAlignObjMatrix(GetDESymName(detElemId), volId, 
				       const_cast<TGeoHMatrix&>(matrix),bGlobal);
}

//______________________________________________________________________________
void AliMUONGeometryTransformer::CreateModules()
{
/// Create modules and their detection elements using info from mapping;
/// but do not fill matrices

  // Load mapping as its info is used to define modules & DEs
  LoadMapping();

  if ( fModuleTransformers->GetEntriesFast() == 0 ) {
    // Create modules only if they do not yet exist

    // Loop over geometry module
    for (Int_t moduleId = 0; moduleId < AliMpConstants::NofGeomModules(); ++moduleId ) {
    
      // Create geometry module transformer
      AliMUONGeometryModuleTransformer* moduleTransformer
        = new AliMUONGeometryModuleTransformer(moduleId);
      AddModuleTransformer(moduleTransformer);
    }
  }     
    
  // Loop over detection elements
  AliMpDEIterator it;
  for ( it.First(); ! it.IsDone(); it.Next() ) {
    
    Int_t detElemId = it.CurrentDEId();
    Int_t moduleId = AliMpDEManager::GetGeomModuleId(detElemId);

    // Get detection element store
    AliMpExMap* detElements = 
      GetModuleTransformer(moduleId)->GetDetElementStore();     

    // Add detection element
    AliMUONGeometryDetElement* detElement 
      = new AliMUONGeometryDetElement(detElemId);
    detElements->Add(detElemId, detElement);
  }   
}

//_____________________________________________________________________________
void AliMUONGeometryTransformer::AddAlignableVolumes() const
{
/// Set symbolic names and matrices to alignable objects to TGeo

  if ( ! gGeoManager ) {
    AliWarning("TGeoManager not defined.");
    return;
  }  

  // Modules 
  for (Int_t i=0; i<fModuleTransformers->GetEntriesFast(); i++) {
    AliMUONGeometryModuleTransformer* module 
      = (AliMUONGeometryModuleTransformer*)fModuleTransformers->At(i);

    // Set module symbolic name
    TGeoPNEntry* pnEntry
      = gGeoManager->SetAlignableEntry(GetModuleSymName(module->GetModuleId()),
                                       module->GetVolumePath());
    if ( ! pnEntry ) {
      AliErrorStream() 
        << "Volume path " << module->GetVolumePath().Data()
        << " for geometry module " << module->GetModuleId() << " " << module
        << " not found in geometry." << endl;
    }
    else {
      // Set module matrix
      pnEntry->SetMatrix(new TGeoHMatrix(*module->GetTransformation()));  
       // the matrix will be deleted via TGeoManager  
    }                                     

    // Detection elements
    AliMpExMap* detElements = module->GetDetElementStore();    
    TIter next(detElements->CreateIterator());    
    AliMUONGeometryDetElement* detElement;
    
    while ( ( detElement = static_cast<AliMUONGeometryDetElement*>(next()) ) )
    {
      // Set detection element symbolic name
      TGeoPNEntry* pnEntryDE
        = gGeoManager->SetAlignableEntry(GetDESymName(detElement->GetId()), 
                                         detElement->GetVolumePath());
      if ( ! pnEntryDE ) {
        AliErrorStream() 
          << "Volume path " 
          << detElement->GetVolumePath().Data() 
          << " for detection element " << detElement->GetId()
          << " not found in geometry." << endl;
      }
      else {
        // Set detection element matrix
        pnEntryDE->SetMatrix(new TGeoHMatrix(*detElement->GetGlobalTransformation()));                                      
         // the matrix will be deleted via TGeoManager 
      }                                      
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
    Int_t volId = AliGeomManager::LayerToVolUID(AliGeomManager::kMUON, moduleId); 

    // Create mis align matrix
    Int_t pos = array->GetEntriesFast();
    new (refArray[pos]) AliAlignObjMatrix(GetModuleSymName(moduleId), volId, matrix, kTRUE);
  }     

  // Detection elements
  for (Int_t i=0; i<fModuleTransformers->GetEntriesFast(); i++) {
    AliMUONGeometryModuleTransformer* moduleTransformer 
      = (AliMUONGeometryModuleTransformer*)fModuleTransformers->At(i);

    AliMpExMap* detElements = moduleTransformer->GetDetElementStore();    
    TIter next(detElements->CreateIterator());    
    AliMUONGeometryDetElement* detElement;
    
    while ( ( detElement = static_cast<AliMUONGeometryDetElement*>(next()) ) )
    {
      Int_t detElemId = detElement->GetId();
  
      // Align object ID
      Int_t volId = AliGeomManager::LayerToVolUID(AliGeomManager::kMUON, detElemId); 

      // Create mis align matrix
      Int_t pos = array->GetEntriesFast();
      new (refArray[pos]) AliAlignObjMatrix(GetDESymName(detElemId), volId, matrix, kTRUE);
    }
  }
  
  return array;
}       

//_____________________________________________________________________________
void AliMUONGeometryTransformer::ClearMisAlignmentData()
{
/// Clear the array of misalignment data

  if ( ! fMisAlignArray ) return;
  
  fMisAlignArray->Delete();
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
    

