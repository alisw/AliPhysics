#if !defined(__CINT__) || defined(__MAKECINT__)
#include "TSystem.h"
#include "TROOT.h"
#include "TGeoManager.h"
#include "TObjString.h"
#include "TClonesArray.h"
#include "TError.h"
#include "AliGeomManager.h"
#include "AliCDBManager.h"
#include "AliCDBStorage.h"
#include "AliCDBPath.h"
#include "AliCDBEntry.h"
#include "AliCDBId.h"
#include "AliCDBMetaData.h"
#include "AliMisAligner.h"
#include "AliHMPIDMisAligner.h"
#include <TString.h>
#endif

void MakeHMPIDAlignmentObjs(Bool_t toOCDB = kFALSE, const char* misalType="residual")
{
  // Make alignment objects for HMPID detector
  // for the misalignment scenario passed as argument "misalType".
  //Input Args:  toOCDB = kFALSE -> the results are written in a local file called HMPIDMisalignObject.root
  //                      kTRUE  -> the results are written in local://$ALICE_ROOT/OCD 
  //             misalType = ideal, residual, full (see class AliHMPIDMisAligner)
  
  const char* macroname = "MakeHMPIDAlignmentObjs.C";
  
  // Load geometry from OCDB; update geometry before loading it if we are going to load
  // the alignment objects to the OCDB
  AliCDBManager* cdb = AliCDBManager::Instance();
  if(!cdb->IsDefaultStorageSet()) cdb->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
  cdb->SetRun(0);
 // AliCDBStorage* storage = 0;
      
  AliGeomManager::LoadGeometry(); //load geom from default OCDB storage

  TClonesArray* objsArray = 0;

  AliHMPIDMisAligner* misAlignerHMPID = new AliHMPIDMisAligner();
  misAlignerHMPID->SetMisalType(misalType);
  objsArray = misAlignerHMPID->MakeAlObjsArray();

    if(toOCDB)
    {
      AliCDBId id("HMPID/Align/Data",0,AliCDBRunRange::Infinity());
      AliCDBMetaData *md = misAlignerHMPID->GetCDBMetaData();
      md->SetResponsible("Domenico Di Bari");
      md->SetAliRootVersion(gSystem->Getenv("ARVERSION"));
      cdb->Put(objsArray,id,md);
    } else {
      // save on file
      TFile file("HMPIDMisalignObject.root","RECREATE");
      if(!file) {
	Error(macroName,"cannot open file for output\n");
	return;
      }
      file.cd();
      file.WriteObject(objsArray,"HMPIDAlignObjs","kSingleKey");
      file.Close();
    }
    
  return;
}
