#ifndef __CINT__
# include <TSystemDirectory.h>
# include <TObjArray.h>
# include <TString.h>
# include <TList.h>
# include <TSystem.h>
#else 
class TSystemDirectory;
class TList;
class TString;
#endif

//__________________________________________________________________
/** 
 * Scan directory @a dir (possibly recursive) for tree files to add
 * to the chain.    This does not follow sym-links
 * 
 * @param dir        Directory to scan
 * @param chain      Chain to add to
 * @param type       Type of tree (ESD or AOD)
 * @param recursive  Whether to scan recursively 
 * @param mc         Look also for MC files if true 
 *
 * @return true if any files where added 
 */
Bool_t ScanDirectory(const TString& pattern, TSystemDirectory* dir, TList& list,
		     Bool_t recursive=false)
{
  // Assume failure 
  Bool_t ret = false;

  // Get list of files, and go back to old working directory
  TString oldDir(gSystem->WorkingDirectory());
  TList*  files = dir->GetListOfFiles();
  if (!gSystem->ChangeDirectory(oldDir)) { 
    Error("ScanDirectory", "Failed to go back to %s", oldDir.Data());
    return false;
  }
  if (!files) return false;

  // Sort list of files and check if we should add it 
  files->Sort();
  TIter next(files);
  TSystemFile* file = 0;
  while ((file = static_cast<TSystemFile*>(next()))) {
    TString name(file->GetName());
    TString title(file->GetTitle());
    TString full(gSystem->ConcatFileName(file->GetTitle(), name.Data()));
    Bool_t  isFile(!file->IsDirectory());
    if (isFile) full = title;
    // Ignore special links 
    if (name == "." || name == "..") { 
      // Info("ScanDirectory", "Ignoring %s", name.Data());
      continue;
    }

    FileStat_t fs;
    if (gSystem->GetPathInfo(full.Data(), fs)) {
      Warning("ScanDirectory", "Cannot stat %s (%s)", full.Data(),
	      gSystem->WorkingDirectory());
      continue;
    }
    // Check if this is a directory 
    if (!isFile) { 
      if (recursive) {
	// if (title[0] == '/') 
	TSystemDirectory* d = new TSystemDirectory(file->GetName(),
						   full.Data());
	if (ScanDirectory(pattern, d, list, recursive))
	  ret = true;
	delete d;
      }
      continue;
    }
    
    // If this is not a root file, ignore 
    if (!name.EndsWith(".root")) continue;

    // If this file does not contain AliESDs, ignore 
    if (!name.Contains(pattern)) { 
      // Info("ScanDirectory", "%s does not match pattern %s", 
      //      name.Data(), fnPattern.Data());
      continue;
    }
    
    // Add 
    // Info("ScanDirectory", "Adding %s", full.Data());
    list.Add(new TObjString(full));
  }

  if (list.GetEntries() > 0) ret = true;

  gSystem->ChangeDirectory(oldDir);
  return ret;
}

//__________________________________________________________________
struct EventInfo
{
  int fRunNo;
} event;


//__________________________________________________________________
void MakeTrendTree(const char* dir=".", const char* pattern="trending")
{
  TList files;
  
  TSystemDirectory* top = new TSystemDirectory(dir, dir);
  if (!ScanDirectory(pattern, top, files)) { 
    return;
  }
  
  TFile* out = TFile::Open("trending_tree.root", "RECREATE");
  if (!out) { 
    Error("MakeTrendTree", "Failed to open output file");
    return;
  }

  TTree* tree = new TTree("T", "T");
  tree->Branch("event", "runno/I", &event);
}

//__________________________________________________________________
//
// EOF
// 
