/**
 * @file   RunFinalQA.C
 * @author Christian Holm Christensen <cholm@master.hehi.nbi.dk>
 * @date   Fri Jan  6 11:46:30 2012
 * 
 * @brief  
 * 
 * @ingroup pwglf_forward_qa_scripts
 */
//____________________________________________________________________
/** 
 * Scan directory (and possibly sub-directories) for trending files
 * 
 * @param dir        Start directory
 * @param list       List to add file names to
 * @param recursive  Whether to scan recursively
 * @param pattern    Pattern filenames must match
 * 
 * @return true on success
 * 
 * @ingroup pwglf_forward_qa_scripts
 */
Bool_t ScanDirectory(TSystemDirectory* dir, 
		     TList* list, 
		     const char* pattern, 
		     bool recursive=false)
{
  TString fnPattern = pattern;

  // Assume failure 
  Bool_t ret = false;

  // Get list of files, and go back to old working directory
  TString oldDir(gSystem->WorkingDirectory());
  TList*  files = dir->GetListOfFiles();
  if (!gSystem->ChangeDirectory(oldDir)) { 
    Error("ScanDirectory", "Failed to go back to %s", oldDir.Data());
    return false;
  }
  if (!files) {
    Warning("ScanDirectory", "No files found in %s", dir->GetName());
    return false;
  }
  // files->ls();

  TList toAdd;
    
  // Sort list of files and check if we should add it 
  // files->Sort();
  TIter next(files);
  TSystemFile* file = 0;
  while ((file = static_cast<TSystemFile*>(next()))) {
    TString name(file->GetName());
    TString title(file->GetTitle());
    TString full(gSystem->ConcatFileName(file->GetTitle(), name.Data()));
    if (file->IsA() == TSystemDirectory::Class()) full = title;
    if (name == "." || name == "..") { 
      continue;
    }

    FileStat_t fs;
    if (gSystem->GetPathInfo(full.Data(), fs)) {
      Warning("ScanDirectory", "Cannot stat %s (%s)", full.Data(),
	      gSystem->WorkingDirectory());
      continue;
    }
    // Check if this is a directory 
    if (file->IsDirectory(full)) { 
      if (recursive) {
	TSystemDirectory* d = new TSystemDirectory(file->GetName(),
						   full.Data());
	if (ScanDirectory(d,list,pattern,recursive))
	  ret = true;
	delete d;
      }
      continue;
    }
    
    // If this is not a root file, ignore 
    if (!name.EndsWith(".root")) {
      continue;
    }

    // If this file does not contain AliESDs, ignore 
    if (!name.Contains(fnPattern)) { 
      continue;
    }
    
    // Add 
    toAdd.Add(new TObjString(full));
  }

  TIter nextAdd(&toAdd);
  TObjString* s = 0;
  while ((s = static_cast<TObjString*>(nextAdd()))) {
    list->Add(s);
  }
  if (toAdd.GetEntries() > 0) ret = true;

  gSystem->ChangeDirectory(oldDir);
  return ret;
}
//____________________________________________________________________
/** 
 * Get the list of trending files
 * 
 * @param input Start directory 
 * 
 * @return List of files 
 * 
 * @ingroup pwglf_forward_qa_scripts
 */
TList*
GetListOfFiles(const char* input=".")
{
  TList* ret = new TList;
  TString dir(input);
  // if (dir == ".") dir = "";

  TString savdir(gSystem->WorkingDirectory());
  TSystemDirectory d(gSystem->BaseName(dir.Data()), dir.Data());
  if (!ScanDirectory(&d, ret, "tree_", false)) { 
    delete ret;
    ret = 0;
  }
  gSystem->ChangeDirectory(savdir);  
  if (ret) ret->Sort();
  return ret;
}

//____________________________________________________________________
/** 
 * 
 * 
 * @param dir Input directory
 * @param prodYear Production year 
 * @param prodLetter Production letter
 * 
 * @ingroup pwglf_forward_qa_scripts
 */
void 
RunFinalQA(const char* dir, Int_t prodYear=0, const char* prodLetter="")
{
   int ret = 0;
   gROOT->SetMacroPath(Form(".:%s",gROOT->GetMacroPath()));
   gSystem->Load("libGpad");
   gSystem->Load("libTree");

   gROOT->LoadMacro("QABase.h+g");
   gROOT->LoadMacro("QAPlotter.C+g");

   QAPlotter p(prodYear, prodLetter[0]);
  
   TList* l = GetListOfFiles(dir);
   TIter next(l);
   TObject* o = 0;
   while ((o = next())) {
     p.AddFile(o->GetName());
   }

   p.Run();
}
//
// EOF
//
