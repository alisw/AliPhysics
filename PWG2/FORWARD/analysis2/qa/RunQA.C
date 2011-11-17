/**
 * @file   RunQA.C
 * @author Christian Holm Christensen <cholm@nbi.dk>
 * @date   Thu Nov 17 11:35:08 2011
 * 
 * @brief  Script to run the QATrender and QAPlotter in one go
 * 
 * @ingroup pwg2_forward_qa_scripts
 */
/** 
 * Scan directory (and possibly sub-directories) for trending files
 * 
 * @param dir        Start directory
 * @param list       List to add file names to
 * @param recursive  Whether to scan recursively
 * 
 * @return true on success
 * 
 * @ingroup pwg2_forward_qa_scripts
 */
Bool_t ScanDirectory(TSystemDirectory* dir, TList* list, 
		     bool recursive=false)
{
  TString fnPattern = "trending_";

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
    // Ignore special links 
    // Info("ScanDirectory", "name=%s title=%s full=%s", 
    //      name.Data(), title.Data(), full.Data());
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
    if (file->IsDirectory(full)) { 
      if (recursive) {
	TSystemDirectory* d = new TSystemDirectory(file->GetName(),
						   full.Data());
	if (ScanDirectory(d,chain,type,recursive,mc))
	  ret = true;
	delete d;
      }
      continue;
    }
    
    // If this is not a root file, ignore 
    if (!name.EndsWith(".root")) {
      // Info("ScanDirectory", "Ignoring non-ROOT file %s", name.Data());
      continue;
    }

    // If this file does not contain AliESDs, ignore 
    if (!name.Contains(fnPattern)) { 
      // Info("ScanDirectory", "%s does not match pattern %s", 
      //      name.Data(), fnPattern.Data());
      continue;
    }
    
    // Add 
    // Info("ScanDirectory", "Adding %s", full.Data());
    toAdd.Add(new TObjString(full));
  }

  TIter nextAdd(&toAdd);
  TObjString* s = 0;
  while ((s = static_cast<TObjString*>(nextAdd()))) {
    // Info("ScanDirectory", "Adding %s", s->GetString().Data());
    list->Add(s);
  }
  if (toAdd.GetEntries() > 0) ret = true;

  gSystem->ChangeDirectory(oldDir);
  return ret;
}
/** 
 * Get the list of trending files
 * 
 * @param input Start directory 
 * 
 * @return List of files 
 * 
 * @ingroup pwg2_forward_qa_scripts
 */
TList*
GetListOfFiles(const char* input=".")
{
  TList* ret = new TList;
  TString dir(input);
  // if (dir == ".") dir = "";

  TString savdir(gSystem->WorkingDirectory());
  TSystemDirectory d(gSystem->BaseName(dir.Data()), dir.Data());
  if (!ScanDirectory(&d, ret, false)) { 
    delete ret;
    ret = 0;
  }
  gSystem->ChangeDirectory(savdir);  
  if (ret) ret->Sort();
  return ret;
}

/** 
 * Run the QATrender and QAPlotter. 
 * 
 * The QATrender is run over the list of files (runs) to produce the
 * file <tt>forward_trending.root</tt> which contains a TTree of QA
 * information - one entry per run.  
 * 
 * The QATrender will also produce two files per run: 
 * 
 * - <tt>qa_<i>runNo</i>.root</tt> which contains TCanvas objects of
 *   the finished plots.
 *
 * - <tt>qa_<i>runNo</i>.pdf</tt> which is a PDF of the TCanvases
 *   mentioned above.
 *
 * The QAPlotter is then run over the  <tt>forward_trending.root</tt>
 * file and produces two files 
 * 
 * - <tt>qa_<i>first-run</i>-<i>last-run</i>.root</tt> which contains
 *   TCanvas objects of the finished plots.  It also contains the
 *   TMultiGraph objects painted in the canvases.
 *
 * - <tt>qa_<i>first-run</i>-<i>last-run</i>.pdf</tt> which is a PDF
 *   of the TCanvases mentioned above.
 * 
 * The QAPlotter will also produce PNGs of each canvas. 
 *
 * if @a runNo is larger than zero, given, then only the that run will
 * be processed and only by QATrender.  In addition, PNGs of each
 * canvas is produced. 
 * 
 * @param runNo (optional) Run number. If greater than 0, only this
 *              run will be processed 
 * @param input Input directory 
 * @param what  (expert) Flag of what to do 
 * @param max   (expert) Maximum number of files to process. 
 * 
 * @ingroup pwg2_forward_qa_scripts
 */
void
RunQA(const char* input=".", Bool_t keep=true, Int_t runNo=-1,
      UShort_t what=0x3)
{
  gROOT->SetMacroPath(Form(".:$(ALICE_ROOT)/PWG2/FORWARD/analysis2/qa:"
			   "$(ALICE_ROOT)/PWG2/FORWARD/analysis2/corrs:%s",
			   gROOT->GetMacroPath()));
  gSystem->AddIncludePath("-I${ALICE_ROOT}/PWG2/FORWARD/analysis2/qa");
  gSystem->Load("libGpad");
  gSystem->Load("libTree");

  if (runNo > 0) what = 0x1; // Only do first pass

  if (what & 0x1) {
    gROOT->LoadMacro("QATrender.C++g");
    QATrender t(keep, runNo > 0);

    TList* l = 0;
    if (runNo < 0) l = GetListOfFiles(input);
    else {
      TObjString* s = new TObjString("");
      s->String() = gSystem->ConcatFileName(input, Form("trending_%09d.root", 
							runNo));
      l = new TList;
      l->Add(s);
    }
    if (!l) {
      Error("RunQA", "No trending files found");
      return;
    }
    // l->ls();
    TIter next(l);
    TObjString* s = 0;
    while ((s = static_cast<TObjString*>(next()))) {
      Info("Run", "Adding  file %s", s->GetName());
      t.AddFile(s->GetName());
    }
    t.Run();
  }
  if (what & 0x2) {
    gROOT->LoadMacro("QAPlotter.C++g");
    QAPlotter p;
    p.Run();
  }
}
//
// EOF
//
