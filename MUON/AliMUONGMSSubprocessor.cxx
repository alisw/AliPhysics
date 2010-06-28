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

//-----------------------------------------------------------------------------
// Class AliMUONGMSSubprocessor
// -----------------------------
// The shuttle subprocessor for GMS data
// Author: Ivana Hrivnacova, IPN Orsay
// 16/09/2006
//-----------------------------------------------------------------------------

#include "AliMUONGMSSubprocessor.h"
#include "AliMUONPreprocessor.h"
#include "AliMpConstants.h"

#include "AliAlignObjMatrix.h"
#include "AliGeomManager.h"
#include "AliCDBMetaData.h"
#include "AliCDBEntry.h"

#include <TTimeStamp.h>
#include <TFile.h>
#include <TArrayI.h>
#include <TClonesArray.h>
#include <TGeoManager.h>
#include <TObjString.h>
#include <Riostream.h>

/// \cond CLASSIMP
ClassImp(AliMUONGMSSubprocessor)
/// \endcond

const Int_t    AliMUONGMSSubprocessor::fgkSystem = AliPreprocessor::kDCS;
const TString  AliMUONGMSSubprocessor::fgkDataId = "GMS";
const TString  AliMUONGMSSubprocessor::fgkMatrixArrayName = "GMSarray";

//______________________________________________________________________________
AliMUONGMSSubprocessor::AliMUONGMSSubprocessor(AliMUONPreprocessor* master) 
  : AliMUONVSubprocessor(master, "GMS", "Upload GMS matrices to OCDB"),
    fTransformer(0)
{
/// Constructor
}

//______________________________________________________________________________
AliMUONGMSSubprocessor::~AliMUONGMSSubprocessor()
{
/// Destructor

  delete fTransformer;
}


//
// private methods
//


//______________________________________________________________________________
Bool_t  AliMUONGMSSubprocessor::Initialize(Int_t /*run*/, 
                                         UInt_t /*startTime*/, UInt_t /*endTime*/)
{
/// Instantiate geometry transformer

  if ( ! fTransformer ) {
    fTransformer = new AliMUONGeometryTransformer();
    fTransformer->CreateModules();
  }  
  return kTRUE;
}                                           

//______________________________________________________________________________
UInt_t AliMUONGMSSubprocessor::ProcessFile(const TString& fileName)
{
/// Convert TGeoHMatrix to AliAlignObjMatrix and fill them into AliTestDataDCS object

  Master()->Log(Form("Processing GMS file %s", fileName.Data()));
  
  // Open root file
  TFile f(fileName.Data());
  if ( ! f.IsOpen() ) {
    Master()->Log(Form("Cannot open file %s",fileName.Data()));
    return 1;
  }  
  
  // Get array with matrices
  TClonesArray* array = (TClonesArray*)f.Get(fgkMatrixArrayName);
  if ( ! array ) {
    Master()->Log(Form("TClonesArray not found in file %s",fileName.Data()));
    return 2;
  }    
  
  // Array to store correspondance between the moduleId 
  // and its corresponding entry in the GMS array.
  TArrayI moduleIdToGMSIndex;
  moduleIdToGMSIndex.Set(AliMpConstants::NofGeomModules());
  for (Int_t i=0; i<AliMpConstants::NofGeomModules(); i++){
    moduleIdToGMSIndex[i]=-1;
  }
  
  Bool_t useGlobalDelta = kTRUE;
  // Get geometry from OCDB
  AliCDBEntry* geoEntry = Master()->GetGeometryFromOCDB();
  TGeoManager *lGeometry = 0x0;
  if (geoEntry) {
    lGeometry = (TGeoManager*) geoEntry->GetObject();
    if (lGeometry) {
      // Set Geometry
      AliGeomManager::SetGeometry(lGeometry);
      useGlobalDelta = kFALSE;
    }
  }

  // Second transformer to convert GMS matrices local delta-transformation into 
  // ALICE alignment objects in the global delta convention
  AliMUONGeometryTransformer *lTransformer = new AliMUONGeometryTransformer();

  // Get mis alignment from reference run for GMS  
  AliCDBEntry* cdbEntry = Master()->GetFromOCDB("Align", "Baseline");
  TClonesArray* refArray = 0x0;
  if (cdbEntry) {
    refArray = (TClonesArray*)cdbEntry->GetObject();
    if (lGeometry && refArray) {      
      AliGeomManager::ApplyAlignObjsToGeom(*refArray);
      lTransformer->LoadGeometryData();
    }
  }
  
  // Convert matrices into Alice alignment objects
  for (Int_t i=0; i<array->GetEntriesFast(); i++ ) {
    TGeoHMatrix* matrix = (TGeoHMatrix*)array->At(i);
    printf("GMS local %i at %i \n",matrix->GetUniqueID(),i);
    matrix->Print();
    fTransformer->AddMisAlignModule(matrix->GetUniqueID(), *matrix, kTRUE);
    lTransformer->AddMisAlignModule(matrix->GetUniqueID(), *matrix, useGlobalDelta);
    moduleIdToGMSIndex[matrix->GetUniqueID()]=i;
  }
  // Store the GMS local delta-transformation 
  TClonesArray* data = const_cast< TClonesArray*>(fTransformer->GetMisAlignmentData());
  
  //Now we have to store the final CDB file
  Master()->Log("Storing GMS");
  AliCDBMetaData metaData;
  metaData.SetBeamPeriod(0);
  metaData.SetResponsible("");
  metaData.SetComment("This preprocessor fills GMS alignment objects.");
  
  Bool_t result = Master()->Store("Align", "GMS", (TObject*)data, &metaData, 0, 0);
  
  // This section apply the GMS misalignments on top of the misalignments 
  // of the GMS reference run and stores the new misalignment array in the OCDB

  // Reset the geoemtry.
  if (geoEntry) {
    lGeometry = (TGeoManager*) geoEntry->GetObject();
  }

  if (lGeometry) {
    TGeoManager::UnlockGeometry();
    // Set Geometry
    AliGeomManager::SetGeometry(lGeometry);
    useGlobalDelta = kFALSE;
  }

  AliAlignObjMatrix* refAOMat = 0x0;
  // First we need to get a copy of the local delta-transformations
  TClonesArray* matLocalArray = new TClonesArray("TGeoHMatrix",200);
  if (lGeometry && refArray) {      
    refArray->Sort(); // All Modules will be first then the DE
    // Create new misalignment array
    for (Int_t i=0; i<refArray->GetEntriesFast(); i++) {
      refAOMat = (AliAlignObjMatrix*)refArray->At(i);      
      refAOMat->Print("");
      TGeoHMatrix* reflMat = new((*matLocalArray)[i]) TGeoHMatrix(); 
      refAOMat->GetLocalMatrix(*reflMat);
      printf("ref misalignment local \n");
      reflMat->Print();
    }
  }
  
  AliAlignObjMatrix* newAOMat = 0x0;
  // Get the GMS global delta-transformation 
  data = const_cast< TClonesArray*>(lTransformer->GetMisAlignmentData());
  if (geoEntry && cdbEntry) {
    if (lGeometry && refArray) {      
      // Create new misalignment array
      TClonesArray* newArray = new TClonesArray("AliAlignObjMatrix", 200);      
      for (Int_t i=0; i<refArray->GetEntriesFast(); i++) {
        refAOMat = (AliAlignObjMatrix*)refArray->At(i);      
        TGeoHMatrix refMat;
        refAOMat->GetMatrix(refMat);
	newAOMat = new((*newArray)[i]) AliAlignObjMatrix(*refAOMat); // Copy the reference misalignment object to the new array ...
        // Need the module containing this module or detection element
        TString sName = refAOMat->GetSymName(); //Format "/MUON/GMx" or "/MUON/GMx/DEy"
        Int_t iGM = sName.Index("GM");
        Int_t iLS = sName.Last('/');
        if (iLS<iGM) { // This is a module
	  sName.Remove(0,iGM+2);
	  Int_t iMod = sName.Atoi();
	  if (moduleIdToGMSIndex[iMod]>=0) {
	    AliAlignObjMatrix* gmsAOMat = (AliAlignObjMatrix*)data->At(moduleIdToGMSIndex[iMod]);
	    TGeoHMatrix gmsMat;
	    gmsAOMat->GetMatrix(gmsMat);
	    printf("GMS global delta %i %i\n",iMod,moduleIdToGMSIndex[iMod]);
	    gmsMat.Print();
	    printf("ref misalignment \n");
	    refMat.Print();
	    refMat.MultiplyLeft(&gmsMat);
	    printf("new misalignment gms*ref\n");
	    refMat.Print();
	  }
	  else {	    
	    Master()->Log(Form("Missing GMS entry for module %d",iMod));
	  }
	  newAOMat->SetMatrix(refMat); // ... and set its matrix 	
	}
	else { // This is a module
	  newAOMat->SetLocalMatrix(*(TGeoHMatrix*)matLocalArray->At(i)); // ... and set its matrix
	}
      }
      
      //Now we have also to store this CDB file
      TString sMDComment(cdbEntry->GetMetaData()->GetComment());
      sMDComment += " GMS";
      Master()->Log("Storing MisAlignment");
      metaData.SetComment(sMDComment);
      
      result = result && Master()->Store("Align", "Data", newArray, &metaData, 0, 1);
    }
    else {
      if (!refArray)
	Master()->Log("Empty entry?");
      if (!lGeometry)
	Master()->Log("Couldn't find TGeoManager in the specified CDB entry?");
    }
  }
  else {
    if (!geoEntry)
      Master()->Log("Could not get Geometry from OCDB! Will not add a new MUON/Align/Data entry!");    
    if (!cdbEntry)
      Master()->Log("Could not get GMS reference misalignment from OCDB! Will not add a new MUON/Align/Data entry!");    
  }
  // Done with applying GMS misalignments on top of misalignment of reference run
  
  // Clear MisAlignArray in transformer
  fTransformer->ClearMisAlignmentData(); 
  
  return (result!=kTRUE);
}  

//
// public methods
//


//______________________________________________________________________________
UInt_t AliMUONGMSSubprocessor::Process(TMap* /*dcsAliasMap*/)
{
/// Process GMS alignment files.
/// Return failure (0) in case procession of some file has failed

  UInt_t result = 1;
  TList* sources = Master()->GetFileSources(fgkSystem, fgkDataId);
  TIter next(sources);
  TObjString* o(0x0);
  while ( ( o = static_cast<TObjString*>(next()) ) ) {
    TString fileName(Master()->GetFile(fgkSystem, fgkDataId, o->GetName()));
    result *= ProcessFile(fileName);
  }
  delete sources;

  return result;
}

