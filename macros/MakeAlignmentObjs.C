#include "ARVersion.h"

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
#include "AliZDCMisAligner.h"
#include <TString.h>
#endif

void MakeAlignmentObjs(const char* detList="ALL", const char* CDBstorage = "local://$HOME/ResidualMisAlignment", const char* outDir="", const char* misalType="residual", Bool_t partialGeom=kFALSE){
  // Make residual misalignment objects for all detectors
  // Pass different "CDBstorage" argument if needed (e.g. to fill
  // conditions' data base on alien) or set it to null string to have
  // the objects saved locally on file 
  // This macro defines the default name and place for the detector-macros
  // in charge of producing the residual misalignment objects as 
  // $ALICE_ROOT/DET/MakeDETResidualMisAlignment.C
  //

  const char* macroName = "MakeAlignmentObjs";
  TString cdbStorage(CDBstorage);
  TString oDir(outDir);
  if(cdbStorage.IsNull() && oDir.IsNull())
  {
    Error(macroName, "Output undefined! Set either the CDB storage or the output directory!");
    return;
  }
  TString fileName;

  TMap misAligners;
  TString modList(detList);
  if(modList=="ALL") modList="ACORDE EMCAL FMD HMPID ITS MUON PMD PHOS T0 TRD TPC TOF VZERO ZDC";
  Info(macroName, "Processing detectors: %s \n", modList.Data());
  Printf("Creating %s misalignment for detectors: %s \n", misalType, modList.Data());
  if(modList.Contains("ZDC")){
    AliZDCMisAligner* misAlignerZDC = new AliZDCMisAligner();
    misAligners.Add(new TObjString("ZDC"), misAlignerZDC);
  }

  // Load geometry from CDB; update geometry before loading it if we are going to load
  // the alignment objects to the CDB
  AliCDBManager* cdb = AliCDBManager::Instance();
  if(!cdb->IsDefaultStorageSet()) cdb->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
  cdb->SetRun(0);
  AliCDBStorage* storage = 0;
  
  if(cdbStorage.IsNull()){ //if we produce the objects into a file
    AliGeomManager::LoadGeometry(); //load geom from default CDB storage
  }else{ // if we produce the objects in a CDB storage
    // update geometry in it
    Info(macroName, "Updating geometry in CDB storage %s",cdbStorage.Data());
    gROOT->ProcessLine(".L $ALICE_ROOT/GRP/UpdateCDBIdealGeom.C");
    if(partialGeom){
      UpdateCDBIdealGeom(cdbStorage.Data(),"$ALICE_ROOT/test/fpprod/Config.C");
    }else{
      UpdateCDBIdealGeom(cdbStorage.Data(),"$ALICE_ROOT/macros/Config.C");
    }
    // load the same geometry from given CDB storage
    AliCDBPath path("GRP","Geometry","Data");
    storage = cdb->GetStorage(cdbStorage.Data());
    AliCDBEntry *entry = storage->Get(path.GetPath(),cdb->GetRun());
    if(!entry) Fatal(macroName,"Couldn't load geometry data from CDB!");
    entry->SetOwner(0);
    TGeoManager* geom = (TGeoManager*) entry->GetObject();
    if (!geom) Fatal(macroName,"Couldn't find TGeoManager in the specified CDB entry!");
    AliGeomManager::SetGeometry(geom);
  }
  
  // run macro for non-sensitive modules
  // (presently generates only FRAME alignment objects)
  // gSystem->Exec("aliroot -b -q $ALICE_ROOT/GRP/MakeSTRUCTResMisAlignment.C"); !!!!!!!!!!!!!!!!!!!!!!!!!

  // run macros for sensitive modules
  TObjString *ostr;
  TString strId;
  TClonesArray* objsArray = 0;

  TObjArray *detArray = modList.Tokenize(' ');
  TIter iter(detArray);

  while((ostr = (TObjString*) iter.Next())){
    TString str(ostr->String());
    if(!oDir.IsNull())
    {
      fileName = oDir;
      fileName += str.Data();
      fileName += misalType;
      fileName += "MisAlignment.root";
    }
    TString arName(str.Data());
    arName += "AlignObjs";
    
    AliMisAligner* misAligner = dynamic_cast<AliMisAligner*> (misAligners.GetValue(str));
    misAligner->SetMisalType(misalType);
    objsArray = misAligner->MakeAlObjsArray();
    //Printf("objsArray has %d entries",objsArray->GetEntriesFast());
    //objsArray->Print();
    if(!cdbStorage.IsNull())
    {
      strId=str;
      strId+="/Align/Data";
      AliCDBId id(strId.Data(),0,AliCDBRunRange::Infinity());
      AliCDBMetaData *md = misAligner->GetCDBMetaData();
      md->SetAliRootVersion(ALIROOT_SVN_BRANCH);
      storage->Put(objsArray, id, md);
    }else{
      // save on file
      TFile file(fileName.Data(),"RECREATE");
      if(!file){
	Error(macroName,"cannot open file for output\n");
	return;
      }
      Info(macroName,"Saving alignment objects to the file %s", fileName);
      file.cd();
      file.WriteObject(objsArray,arName.Data(),"kSingleKey");
      file.Close();
    }
  }

  return;
}
