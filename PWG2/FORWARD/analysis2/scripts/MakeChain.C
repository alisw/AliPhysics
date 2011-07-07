/** 
 * Scan a directory (optionally recursive) for data files to add to
 * the chain.  Only ROOT files, and files which name contain the
 * passed pattern are considered.
 * 
 * @param dir        Directory to scan
 * @param chain      Chain to add data to 
 * @param pattern    Pattern that the file name must contain
 * @param recursive  Whether to scan recursively 
 *
 * @ingroup pwg2_forward_scripts
 */
void
ScanDirectory(TSystemDirectory* dir, TChain* chain, 
	      const char* pattern, bool recursive)
{
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
	ScanDirectory(static_cast<TSystemDirectory*>(file),chain,
		      pattern,recursive);
      continue;
    }
    
    // If this is not a root file, ignore 
    if (!name.EndsWith(".root")) continue;

    // If this file does not contain the pattern, ignore 
    if (!name.Contains(pattern)) continue;
    if (name.Contains("friends")) continue;
    
    // Get the path 
    TString data(Form("%s/%s", file->GetTitle(), name.Data()));

    TFile* test = TFile::Open(data.Data(), "READ");
    if (!test || test->IsZombie()) { 
      Warning("ScanDirectory", "Failed to open file %s", data.Data());
      continue;
    }
    test->Close();
    chain->Add(data);

  }
}

/** 
 * Make a chain of specified data 
 * 
 * @param what       What data to chain.  Possible values are 
 *                   - ESD Event summary data (AliESD)
 *                   - AOD Analysis object data (AliAOD)
 *                   - MC  Simulation data (galice)
 * @param datadir    Data directory to scan 
 * @param recursive  Whether to recurse into sub-directories 
 * 
 * @return Pointer to newly create chain, or null
 *
 * @ingroup pwg2_forward_scripts
 */
TChain*
MakeChain(const char* what, const char* datadir, bool recursive=false)
{
  TString w(what);
  w.ToUpper();
  const char* treeName = 0;
  const char* pattern  = 0;
  if      (w.Contains("ESD")) { treeName = "esdTree"; pattern = "AliESD"; }
  else if (w.Contains("AOD")) { treeName = "aodTree"; pattern = "AliAOD"; }
  else if (w.Contains("MC"))  { treeName = "TE";      pattern = "galice"; }
  else {
    Error("MakeChain", "Unknown mode '%s' (not one of ESD,AOD, or MC)", what);
    return 0;
  }
    
  // --- Our data chain ----------------------------------------------
  TChain* chain = new TChain(treeName);

  // --- Get list of ESDs --------------------------------------------
  // Open source directory, and make sure we go back to were we were 
  TString oldDir(gSystem->WorkingDirectory());
  TSystemDirectory d(datadir, datadir);
  ScanDirectory(&d, chain, pattern, recursive);

  return chain;
}
//
// EOF
//
