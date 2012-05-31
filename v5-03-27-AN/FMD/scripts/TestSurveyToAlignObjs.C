void
TestSurveyToAlignObjs(Bool_t cdbStore=false)
{
  AliCDBManager* cdb = AliCDBManager::Instance();
  cdb->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
  cdb->SetRun(0);
  AliGeomManager::LoadGeometry("geometry.root");

  const char* files[] = { 
    "Survey_943928_FMD.txt", 
    "Survey_976326_FMD.txt", 
    0 
  };

  AliFMDSurveyToAlignObjs convert;
  convert.Run(files);
  convert.GetAlignObjArray()->Print();

  TClonesArray* a = convert.GetAlignObjArray();
  AliAlignObjParams* p = 0;
  for (Int_t i = 0; i < a->GetEntries(); i++) { 
    p = static_cast<AliAlignObjParams*>(a->At(i));
    Info("TestSurveyToAlignObjs", "%30s", p->GetSymName());
  }

  if (!cdbStore) 
    convert.StoreAlignObjToFile("FMD_Survey.root", "FMD");
  else 
    convert.StoreAlignObjToCDB("FMD/Align/Data", "FMD");
}

void
ShowExisting()
{
  TFile*             f = TFile::Open("$ALICE_ROOT/OCDB/FMD/Align/Data/Run0_999999999_v0_s0.root", "READ");
  AliCDBEntry*       e = static_cast<AliCDBEntry*>(f->Get("AliCDBEntry"));
  TClonesArray*      a = static_cast<TClonesArray*>(e->GetObject());
  AliAlignObjParams* p = 0;
  for (Int_t i = 0; i < a->GetEntries(); i++) { 
    p = (AliAlignObjParams*)a->At(i); 
    Info("ShowExisting", "%s %d", p->GetSymName(), p->GetVolUID()); 
  } 
}
