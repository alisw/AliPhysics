TChain*
MakeESDChain(const char* esddir)
{
  // --- Our data chain ----------------------------------------------
  TChain* chain = new TChain("esdTree");

  // --- Get list of ESDs --------------------------------------------
  // Open source directory, and make sure we go back to were we were 
  TString oldDir(gSystem->WorkingDirectory());
  TSystemDirectory d(esddir, esddir);
  TList* files = d.GetListOfFiles();
  gSystem->ChangeDirectory(oldDir);

  // Sort list of files and check if we should add it 
  files->Sort();
  TIter next(files);
  TSystemFile* file = 0;
  while ((file = static_cast<TSystemFile*>(next()))) {
    if (file->IsDirectory()) continue;
    TString name(file->GetName());
    if (!name.EndsWith(".root")) continue;
    if (!name.Contains("AliESDs")) continue;
    TString esd(Form("%s/%s", file->GetTitle(), name.Data()));
    Info("RunManager", "Adding %s to chain", esd.Data());
    chain->Add(esd);
  }  
  return chain;
}
