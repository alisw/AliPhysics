void MakePHOSFullMisAlignment()
{
  // Create misalignment object for PHOS module 2 
  // from the survey measurements on P2 in May 2008.
  // To store alignment in OCDB, define the evnironment variables:
  // TOCDB=kTRUE
  // STORAGE="local://$ALICE_ROOT/OCDB"
  // Author: Timur Pocheptsov, 19.06.2008

  const char * macroName = "MakePHOS2Alignment";
  
  const AliPHOSGeometry *phosGeom = AliPHOSGeometry::GetInstance("IHEP", "IHEP");
  if (!phosGeom) {
    Error(macroName, "Cannot obtain AliPHOSGeometry singleton.\n");
    return;
  }

  //Activate CDB storage and load geometry from CDB
  //I need valid gGeoManager to create local transformations.
  
  //[Part of code, taken from ITS version of MakexxxFullMisalignment.
  //This code was required only for local delta transformations.
  //Now, transformations are global, so no gGeoManager and geometry are required.
  /*
  AliCDBManager * cdb = AliCDBManager::Instance();
  if (!cdb->IsDefaultStorageSet())
    cdb->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
  cdb->SetRun(0);

  if (TString(gSystem->Getenv("TOCDB")) == TString("kTRUE")) {
    const TString storageENV(gSystem->Getenv("STORAGE"));
    if (!storageENV.BeginsWith("local://") && !storageENV.BeginsWith("alien://")) {
      Error(macroName, "STORAGE variable set to %s is not valid.\n", storageENV.Data());
      return;
    }
    
    AliCDBStorage * storage = cdb->GetStorage(storageENV.Data());
    if (!storage) {
      Error(macroName, "Unable to open storage %s.\n", storageENV.Data());
      return;
    }
    
    AliCDBPath path("GRP","Geometry","Data");
    AliCDBEntry * entry = storage->Get(path.GetPath(), cdb->GetRun());
    if (!entry) {
      Error(macroName,"Could not get the specified CDB entry!");
      return;  
    }
    
    entry->SetOwner(0);
    AliGeomManager::SetGeometry((TGeoManager*) entry->GetObject());
  }else{
    AliGeomManager::LoadGeometry("geometry.root"); //load geom from default CDB storage
  }  
  */  
  //end of code taken from ITS version of MakexxxFullMisalignment]

  AliPHOSEMCAGeometry * emca = phosGeom->GetEMCAGeometry();
  TClonesArray alobj("AliAlignObjParams", 16);// + phosGeom->GetNModules() * emca->GetNStripX() *
                                              //   emca->GetNStripZ());
  
  const Double_t dpsi = 0., dtheta = 0., dphi = 0.;
  const Double_t displacement = 10.;
  Int_t iIndex = 0; //let all modules have index=0 in a layer with no LUT
  const AliGeomManager::ELayerID iLayer = AliGeomManager::kInvalidLayer;
  UShort_t volid = AliGeomManager::LayerToVolUID(iLayer, iIndex);
  Int_t i = 0;

  // Alignment for 5 PHOS modules

  TString surveyFileName;
  
  const Char_t * szEnv = gSystem->Getenv("ALICE_ROOT");
  if (szEnv && szEnv[0]) {
    surveyFileName += szEnv;
    if (surveyFileName[surveyFileName.Length() - 1] != '/')
      surveyFileName += '/';
  } else {
    Warning(macroName, "cannot find ALICE_ROOT environment variable\n"
		       "probably, I wan't be able to find survey file");
  }
    
  surveyFileName += "PHOS/data/phos_mod3_survey_2008.txt";

  new(alobj[i++]) AliAlignObjParams("PHOS/Module1", volid, 0., 0., 0., 0., 0., 0., kTRUE);
  new(alobj[i++]) AliAlignObjParams("PHOS/Module2", volid, 2., 0., 0., 0., 0., 0., kTRUE);
  //
  AliSurveyObj survey;
  survey.FillFromLocalFile(surveyFileName.Data());
  TGeoHMatrix module3Delta;
  AliPHOSModuleMisalignment delta(*phosGeom);
  delta.DeltaTransformation(2, survey.GetData(), "T2_10000", "T2_10027", "T2_24000", 
        		    &module3Delta);
  new(alobj[i++]) AliAlignObjParams("PHOS/Module3", volid, module3Delta, kTRUE);
  //
  new(alobj[i++]) AliAlignObjParams("PHOS/Module4", volid, 0.,  0., 0., 0., 0., 0., kTRUE);
  new(alobj[i++]) AliAlignObjParams("PHOS/Module5", volid, 0., 0., 0., 0., 0., 0., kTRUE);

  const Double_t dx = 0., dy = 0., dz = 0. ;
  // Alignment of CPV modules
  new(alobj[i++]) AliAlignObjParams("PHOS/Module1/CPV",
        volid, dx, dy, dz, dpsi, dtheta, dphi, kTRUE);
  new(alobj[i++]) AliAlignObjParams("PHOS/Module2/CPV",
        volid, dx, dy, dz, dpsi, dtheta, dphi, kTRUE);
  new(alobj[i++]) AliAlignObjParams("PHOS/Module3/CPV",
        volid, dx, dy, dz, dpsi, dtheta, dphi, kTRUE);
  new(alobj[i++]) AliAlignObjParams("PHOS/Module4/CPV",
        volid, dx, dy, dz, dpsi, dtheta, dphi, kTRUE);
  new(alobj[i++]) AliAlignObjParams("PHOS/Module5/CPV",
        volid, dx, dy, dz, dpsi, dtheta, dphi, kTRUE);
 
  // Alignment for PHOS cradle
  new(alobj[i++]) AliAlignObjParams("PHOS/Cradle0",
	  volid, 0., 0., -displacement, dpsi, dtheta, dphi, kTRUE);
  new(alobj[i++]) AliAlignObjParams("PHOS/Cradle1",
	  volid, 0., 0., +displacement, dpsi, dtheta, dphi, kTRUE);

  // Alignment for cradle wheels
  new(alobj[i++]) AliAlignObjParams("PHOS/Wheel0",
	  volid, 0., 0., -displacement, dpsi, dtheta, dphi, kTRUE);
  new(alobj[i++]) AliAlignObjParams("PHOS/Wheel1",
	  volid, 0., 0., -displacement, dpsi, dtheta, dphi, kTRUE);
  new(alobj[i++]) AliAlignObjParams("PHOS/Wheel2",
	  volid, 0., 0., +displacement, dpsi, dtheta, dphi, kTRUE);
  new(alobj[i++]) AliAlignObjParams("PHOS/Wheel3",
	  volid, 0., 0., +displacement, dpsi, dtheta, dphi, kTRUE);

  // *************************    2nd step    ***************
  if ( TString(gSystem->Getenv("TOCDB")) != TString("kTRUE") ) {
    // save on file
    const char * fileName = "PHOSfullMisalignment.root";
    TFile f(fileName,"RECREATE");
    if (!f) {
      Error(macroName, "cannot open file for output\n");
      return;
    }
    
    Info(macroName,"Saving alignment objects to the file %s", fileName);
    f.cd();
    f.WriteObject(&alobj,"PHOSAlignObjs","kSingleKey");
    f.Close();
  }else{
    // save in CDB storage
    TString storageENV = gSystem->Getenv("STORAGE");
    if(!storageENV.BeginsWith("local://") && !storageENV.BeginsWith("alien://")) {
      Error(macroName,"STORAGE variable set to %s is not valid. Exiting\n", storageENV.Data());
      return;
    }
    
    Info(macroName,"Saving alignment objects in CDB storage %s", storageENV.Data());
    AliCDBManager * cdb = AliCDBManager::Instance();
    AliCDBStorage * storage = cdb->GetStorage(storageENV.Data());
    if (!storage) {
      Error(macroName, "Unable to open storage %s\n", storageENV.Data());
      return;
    }
    
    AliCDBMetaData md;
    md.SetResponsible("Timur Pocheptsov");
    md.SetComment("Alignment objects for PHOS module 2; survey in May 2008");
    md.SetAliRootVersion(gSystem->Getenv("ARVERSION"));
    AliCDBId id("PHOS/Align/Data",0,AliCDBRunRange::Infinity());
    storage->Put(&alobj, id, &md);
  }

  alobj.Delete();
}
