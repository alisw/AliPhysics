/**
 * @file   MakeChain.C
 * @author Christian Holm Christensen <cholm@dalsgaard.hehi.nbi.dk>
 * @date   Tue Jul 12 10:20:07 2011
 * 
 * @brief  Script to generate a chain of files 
 * 
 * @ingroup pwg2_forward_scripts
 */
/** 
 * Check if a path points to a file 
 * 
 * @param path Path
 * 
 * @return True if the path points to a regular file 
 *
 * @ingroup pwg2_forward_scripts
 */
Bool_t
IsFile(const char* path)
{
  Long_t id;
  Long_t size;
  Long_t flags;
  Long_t modtime;
  gSystem->GetPathInfo(path, &id, &size, &flags, &modtime);
  return !((flags & 0x2) == 0x2);
}

/** 
 * Test if we can open a file 
 * 
 * @param name    Name of file 
 * @param pattern Pattern to check against 
 * 
 * @return True on success
 */
Bool_t
TestFile(const TString& name, const char* pattern=0)
{
  // If this is not a root file, ignore 
  if (!name.EndsWith(".root")) return false;

  // If this file does not contain the pattern, ignore 
  if (pattern && pattern[0] != '\0' && !name.Contains(pattern)) return false;
  if (name.Contains("friends")) return false;
    
  Bool_t ret  = true;
  TFile* test = TFile::Open(name.Data(), "READ");
  if (!test || test->IsZombie()) { 
    Warning("TestFile", "Failed to open file %s", name.Data());
    ret = false;
  }
  else 
    test->Close();
  return ret;
}

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
    
    // Get the path 
    TString data(Form("%s/%s", file->GetTitle(), name.Data()));

    // Check the fuile 
    if (!TestFile(data, pattern)) continue;
    chain->Add(data);
  }
}
/** 
 * Scan an input list of files 
 * 
 * @param chain Chain to add to 
 * @param path  file with list of files to add 
 * 
 * @return true on success 
 */
Bool_t 
ScanInputList(TChain* chain, const TString& path, const char* treeName=0)
{
  std::ifstream in(path.Data()); 
  if (!in) { 
    Error("ScanInputList", "Failed to open input list %s", path.Data());
    return false;
  }
  TString line;
  while (in.good()) { 
    line.ReadLine(in); // Skip white-space
    if (line.IsNull()) break; // Nothing -> EOF
    if (line[0] == '#') continue; // Ignore comment lines 
    if (!TestFile(line, 0)) continue; 
    chain->Add(line);
  }
  in.close();
  return true;
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

  // --- Get list of files --------------------------------------------
  // Open source directory, and make sure we go back to were we were 
  TString oldDir(gSystem->WorkingDirectory());
  TString path(gSystem->ExpandPathName(datadir));
  if (!IsFile(path)) {
    TSystemDirectory d(datadir, datadir);
    ScanDirectory(&d, chain, pattern, recursive);
  }
  else if (path.EndsWith(".root")) { 
    if (TestFile(path, pattern)) chain->Add(path);
  }
  else { 
    // Input seems to be a list - parse it 
    ScanInputList(chain, path);
  }
  
  // Make sure we do not make an empty chain 
  if (chain->GetListOfFiles()->GetEntries() <= 0) { 
    Warning("MakeChain", "Chain %s is empty for input %s", 
	    treeName, datadir);
    delete chain;
    chain = 0;
  }
  return chain;
}
//
// EOF
//
