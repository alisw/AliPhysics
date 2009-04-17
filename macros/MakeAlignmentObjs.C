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
#include "AliHMPIDMisAligner.h"
#include "AliITSMisAligner.h"
#include "AliPMDMisAligner.h"
#include "AliT0MisAligner.h"
#include "AliTPCMisAligner.h"
#include "AliVZEROMisAligner.h"
#include "AliZDCMisAligner.h"
#include <TString.h>
#endif

void MakeAlignmentObjs(const char* detList="ALL", const char* ocdbOrDir = "local://$HOME/ResidualMisAlignment", const char* misalType="residual", Bool_t partialGeom=kFALSE){
  // Make alignment objects for all detectors listed in "detList"
  // for the misalignment scenario passed as argument "misalType".
  // "ocdbUriDirPath" argument is used as URI for an OCDB if it contains
  // either the string "local;//" or the string "alien://folder=",
  // otherwise it is used as the path of the directory where to
  // put the files containing the alignment objects.
  // The geometry used is the one produced with $ALICE_ROOT/macros/Config.C
  // unless "partialGeom" is set to true (=> $ALICE_ROOT/test/fpprod/Config.C).
  //

  const char* macroName = "MakeAlignmentObjs";
  Bool_t toOCDB = kFALSE;
  TString fileName("");
  TString ocdbUriDirPath(ocdbOrDir);
  if(ocdbUriDirPath.IsNull() || ocdbUriDirPath.IsWhitespace())
  {
    Error(macroName, "Output undefined! Set it either to a valid OCDB storage or to the output directory!");
    return;
  }else if(ocdbUriDirPath.Contains("local://") || ocdbUriDirPath.Contains("alien://folder=")){
      // else ocdbUriDirPath is to be interpreted as an OCDB URI
      toOCDB=kTRUE;
      Printf("Objects will be saved in the OCDB %s",ocdbUriDirPath.Data());  
  }else{ // else ocdbUriDirPath is to be interpreted as a directory path
      gSystem->ExpandPathName(ocdbUriDirPath);
      if(gSystem->AccessPathName(ocdbUriDirPath.Data()))
      {
	  Printf("Directory \"%s\" where to save files does not yet exist! ... exiting!",ocdbUriDirPath.Data());
	  return;
      }else{
	  Printf("Files with alignment objects will be saved in the directory %s",ocdbUriDirPath.Data());  
      }
  }
      
  TMap misAligners;
  TString modList(detList);
  if(modList=="ALL") modList="ACORDE EMCAL FMD HMPID ITS MUON PMD PHOS T0 TRD TPC TOF VZERO ZDC";
  Info(macroName, "Processing detectors: %s \n", modList.Data());
  Printf("Creating %s misalignment for detectors: %s \n", misalType, modList.Data());
  if(modList.Contains("EMCAL")){
    AliEMCALMisAligner* misAlignerEMCAL = new AliEMCALMisAligner();
    misAligners.Add(new TObjString("EMCAL"), misAlignerEMCAL);
  }
  if(modList.Contains("HMPID")){
    AliHMPIDMisAligner* misAlignerHMPID = new AliHMPIDMisAligner();
    misAligners.Add(new TObjString("HMPID"), misAlignerHMPID);
  }
  if(modList.Contains("ITS")){
    AliITSMisAligner* misAlignerITS = new AliITSMisAligner();
    misAligners.Add(new TObjString("ITS"), misAlignerITS);
  }
  if(modList.Contains("PMD")){
    AliPMDMisAligner* misAlignerPMD = new AliPMDMisAligner();
    misAligners.Add(new TObjString("PMD"), misAlignerPMD);
  }
  if(modList.Contains("T0")){
    AliT0MisAligner* misAlignerT0 = new AliT0MisAligner();
    misAligners.Add(new TObjString("T0"), misAlignerT0);
  }
  if(modList.Contains("TPC")){
    AliTPCMisAligner* misAlignerTPC = new AliTPCMisAligner();
    misAligners.Add(new TObjString("TPC"), misAlignerTPC);
  }
  if(modList.Contains("VZERO")){
    AliVZEROMisAligner* misAlignerVZERO = new AliVZEROMisAligner();
    misAligners.Add(new TObjString("VZERO"), misAlignerVZERO);
  }
  if(modList.Contains("ZDC")){
    AliZDCMisAligner* misAlignerZDC = new AliZDCMisAligner();
    misAligners.Add(new TObjString("ZDC"), misAlignerZDC);
  }

  // Load geometry from OCDB; update geometry before loading it if we are going to load
  // the alignment objects to the OCDB
  AliCDBManager* cdb = AliCDBManager::Instance();
  if(!cdb->IsDefaultStorageSet()) cdb->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
  cdb->SetRun(0);
  AliCDBStorage* storage = 0;
  
  if(!toOCDB){ //if we produce the objects into a file
    AliGeomManager::LoadGeometry(); //load geom from default OCDB storage
  }else{ // if we produce the objects in a OCDB storage
    // update geometry in it
    Info(macroName, "Updating geometry in OCDB storage %s",ocdbUriDirPath.Data());
    gROOT->ProcessLine(".L $ALICE_ROOT/GRP/UpdateCDBIdealGeom.C");
    if(partialGeom){
      UpdateCDBIdealGeom(ocdbUriDirPath.Data(),"$ALICE_ROOT/test/fpprod/Config.C");
    }else{
      UpdateCDBIdealGeom(ocdbUriDirPath.Data(),"$ALICE_ROOT/macros/Config.C");
    }
    // load the same geometry from given OCDB storage
    AliCDBPath path("GRP","Geometry","Data");
    storage = cdb->GetStorage(ocdbUriDirPath.Data());
    AliCDBEntry *entry = storage->Get(path.GetPath(),cdb->GetRun());
    if(!entry) Fatal(macroName,"Couldn't load geometry data from OCDB!");
    entry->SetOwner(0);
    TGeoManager* geom = (TGeoManager*) entry->GetObject();
    if (!geom) Fatal(macroName,"Couldn't find TGeoManager in the specified OCDB entry!");
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
    TString str(ostr->String()); // DET
    TString arName(str.Data());  // name of the array in case saved into the file
    arName += "AlignObjs";
    
    AliMisAligner* misAligner = dynamic_cast<AliMisAligner*> (misAligners.GetValue(str));
    misAligner->SetMisalType(misalType);
    objsArray = misAligner->MakeAlObjsArray();

    if(toOCDB)
    {
      strId=str;
      strId+="/Align/Data";
      AliCDBId id(strId.Data(),0,AliCDBRunRange::Infinity());
      AliCDBMetaData *md = misAligner->GetCDBMetaData();
      md->SetAliRootVersion(ALIROOT_SVN_BRANCH);
      storage->Put(objsArray, id, md);
    }else{
      // save on file
      fileName = ocdbUriDirPath;
      fileName += "/";
      fileName += str.Data();
      fileName += misalType;
      fileName += "MisAlignment.root";
      TFile file(fileName.Data(),"RECREATE");
      if(!file){
	Error(macroName,"cannot open file for output\n");
	return;
      }
      Info(macroName,"Saving alignment objects to the file %s", fileName.Data());
      file.cd();
      file.WriteObject(objsArray,arName.Data(),"kSingleKey");
      file.Close();
    }
  }

  return;
}
