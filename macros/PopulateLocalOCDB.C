void PopulateLocalOCDB(const char* storageUri="local://$ALICE_ROOT/OCDB", const char* dets="ALL", Int_t runMin=0, Int_t runMax=0)
{
  // dets: comma-separated list of detectors for which to make the CDB objects
  //
  /*
  TString storageString(storageUri);
  if(!storageUri || storageString.IsNull()) {
    storageUri = gSystem->ExpandPathName("local://$ALICE_ROOT/OCDB");
  }*/
  AliCDBManager *cdb = AliCDBManager::Instance();
  cdb->SetDefaultStorage(storageUri);
  if (runMax==0)
    runMax=AliCDBRunRange::Infinity();
  const Int_t nDets=20;
  const char* detectorNames[nDets] = {"ACORDE", "AD", "EMCAL", "FIT", "FMD", "GRP", "HLT", "HMPID", "ITS", "MFT", "MUON", "PHOS", "PMD", "T0", "TOF", "TPC", "TRD", "TRIGGER", "VZERO", "ZDC"};
  TString detectors(dets);
  if (dets == "ALL") {
    detectors = "";
    for (Int_t i=0; i<nDets; ++i) {
      detectors += detectorNames[i];
      detectors += ",";
    }
    detectors.Chop();
  }

  // For the 'MakeZeroMisalignment' macros
  gSystem->Setenv("TOCDB", "kTRUE");
  gSystem->Setenv("STORAGE", storageUri);

  TObjArray *detsArray = detectors.Tokenize(',');
  TString macroPath = gROOT->GetMacroPath();
  TIter iter(detsArray);
  TObjString *oStringDet = 0;
  TString sDetMacro("");
  while ((oStringDet = (TObjString*) iter.Next()))
  {
    TString stringDet(oStringDet->String());
    stringDet.ReplaceAll(" ","");
    Printf("\n   ----  Making conditions'objects for \"%s\"  ----", stringDet.Data());
    sDetMacro = "Make";
    sDetMacro += stringDet;
    sDetMacro += "CDBObjects.C";
    if (! macroPath.EndsWith(":")) {
      macroPath += ":";
    }
    macroPath += gSystem->ExpandPathName("${ALICE_ROOT}/");
    macroPath += stringDet;
    macroPath += "/macros";
    gROOT->SetMacroPath(macroPath);
    gROOT->LoadMacro(sDetMacro.Data());
    sDetMacro = (sDetMacro(0, sDetMacro.Length()-2));
    sDetMacro += "()";
    gROOT->ProcessLine(sDetMacro.Data());
  }
}
