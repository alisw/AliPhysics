void
WriteMedArrays()
{
  TFile* file = TFile::Open("medid.root", "RECREATE");
  if (!file) {
    Warning("WriteMedArrays", "failed to open medid.root");
    return;
  }
  TObjArray* modules = gAlice->Modules();
  if (!modules) {
    Warning("WriteMedArrays", "failed to get modules");
    return;
  }
  TIter next(modules);
  AliModule* module = 0;
  while ((module = static_cast<AliModule*>(next()))) {
    Info("WriteMedArrays", "Getting medium id's for %s", module->GetName());
    TArrayI* mediumIds = module->GetIdtmed();
    if (!mediumIds) {
      Warning("WriteMedArrays", "No medium id's for %s", module->GetName());
      continue;
    }
    file->WriteObject(mediumIds,module->GetName());
  }
  file->Write();
  file->Close();
}

  
  
