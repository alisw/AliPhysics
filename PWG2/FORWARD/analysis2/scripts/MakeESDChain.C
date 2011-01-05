void
ScanDirectory(TSystemDirectory* dir, TChain* chain, bool recursive)
{
  // gROOT->IndentLevel();
  // Printf("Scanning %s ...", dir->GetName());
  // gROOT->IncreaseDirLevel();
  
  
  // Get list of files, and go back to old working directory
  TString oldDir(gSystem->WorkingDirectory());
  TList* files = dir->GetListOfFiles();
  gSystem->ChangeDirectory(oldDir);

  // Sort list of files and check if we should add it 
  files->Sort();
  TIter next(files);
  TSystemFile* file = 0;
  while ((file = static_cast<TSystemFile*>(next()))) {
    TString name(file->GetName());
    
    // Ignore special links 
    if (name == "." || name == "..") continue;

    // Check if this is a directory 
    if (file->IsDirectory()) { 
      if (recursive) 
	ScanDirectory(static_cast<TSystemDirectory*>(file),chain,recursive);
      continue;
    }
    
    // If this is not a root file, ignore 
    if (!name.EndsWith(".root")) continue;

    // If this file does not contain AliESDs, ignore 
    if (!name.Contains("AliESDs")) continue;
    
    // Get the path 
    TString esd(Form("%s/%s", file->GetTitle(), name.Data()));

    // Print and add 
    // gROOT->IndentLevel();
    // Printf("adding %s", esd.Data());
    chain->Add(esd);

  }
  // gROOT->DecreaseDirLevel();
}

TChain*
MakeESDChain(const char* esddir, bool recursive=false)
{
  // --- Our data chain ----------------------------------------------
  TChain* chain = new TChain("esdTree");

  // --- Get list of ESDs --------------------------------------------
  // Open source directory, and make sure we go back to were we were 
  TString oldDir(gSystem->WorkingDirectory());
  TSystemDirectory d(esddir, esddir);
  ScanDirectory(&d, chain, recursive);

  // chain->GetListOfFiles()->ls();

  return chain;
}
