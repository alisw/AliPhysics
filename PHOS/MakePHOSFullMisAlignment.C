void MakePHOSFullMisAlignment()
{
  // Create misalignment object for PHOS module 2,3,3
  // from the survey measurements on P2 in August 2009.
  // To store alignment in OCDB, define the evnironment variables:
  // TOCDB=kTRUE
  // STORAGE="local://$ALICE_ROOT/OCDB"
  // Author: Timur Pocheptsov, 19.06.2008
  // Modified: Yuri Kharlov, 11.03.2010

  const char * macroName = "MakePHOS2Alignment";
  
  const AliPHOSGeometry *phosGeom = AliPHOSGeometry::GetInstance("IHEP", "IHEP");
  if (!phosGeom) {
    Error(macroName, "Cannot obtain AliPHOSGeometry singleton.\n");
    return;
  }

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
    
  surveyFileName += "PHOS/data/Survey_1053236_PHOS.txt";

  AliSurveyObj survey;
  survey.FillFromLocalFile(surveyFileName.Data());
  TGeoHMatrix module3Delta, module2Delta, module1Delta;
  AliPHOSModuleMisalignment delta(*phosGeom);

  delta.DeltaTransformation(0, survey.GetData(), "410000", "410027", "424000", 
        		    &module1Delta);
  new(alobj[i++]) AliAlignObjParams("PHOS/Module1", volid, module1Delta, kTRUE);

  delta.DeltaTransformation(1, survey.GetData(), "310000", "310027", "324000", 
        		    &module2Delta);
  new(alobj[i++]) AliAlignObjParams("PHOS/Module2", volid, module2Delta, kTRUE);

  delta.DeltaTransformation(2, survey.GetData(), "210000", "210027", "224000", 
        		    &module3Delta);
  new(alobj[i++]) AliAlignObjParams("PHOS/Module3", volid, module3Delta, kTRUE);

  new(alobj[i++]) AliAlignObjParams("PHOS/Module4", volid, 0., 0., 0., 0., 0., 0., kTRUE);

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
    md.SetResponsible("Yuri Kharlov");
    md.SetComment("Alignment objects for PHOS modules 2,3,4; survey in August 2009");
    md.SetAliRootVersion(gSystem->Getenv("ARVERSION"));
    AliCDBId id("PHOS/Align/Data",0,AliCDBRunRange::Infinity());
    storage->Put(&alobj, id, &md);
  }

  alobj.Delete();
}
