#ifndef TRAINWORKER_C
#define TRAINWORKER_C
#include "TrainSetup.C"

//====================================================================
/**
 * Base class for workers 
 * 
 */
class TrainWorker
{
public:
  TrainWorker(AliAnalysisGrid* gridHandler) 
    : fHandler(gridHandler)
  {}
  virtual Bool_t Connect(const TUrl& server,
			 const TString& name) = 0;
  virtual void SetROOTVersion(const char* v) {}
  virtual void SetAliROOTVersion(const char* v) {}
  virtual void SetAliEnVersion(const char* v) {}
  virtual void SetUsePar(Bool_t usePar) {}
  virtual void SetInput(const char* d) {}
protected:
  AliAnalysisGrid* fHandler;
};

//====================================================================
/** 
 * Utility class to set-up a chain - either from a base directory, or
 * from an XML file
 */
class TrainChainWorker
{
  /** 
   * Make a chain from a base directory, pattern, and treename -
   * possibly recursively
   * 
   * @param dir         Base directory 
   * @param pattern     Pattern that file names should match
   * @param treeName    Name of tree 
   * @param recursive   Wether to scan recursively 
   * 
   * @return The created chain, or null
   */
  TChain* CreateChain(const TString& dir, 
		      const TString& pattern,
		      const TString& treeName,
		      Bool_t         mc,
		      Bool_t         recursive=true) {
    TChain* chain = new TChain(treeName.Data());
    
    if (dir == ".") dir = "";
    if (!dir.BeginsWith("/")) dir = Form("../%s", dir.Data());
    FileStat_t stat; 
    gSystem->GetPathInfo(dir, stat);
    if (!R_ISDIR(stat.fMode)) { // A file, check it 
      if (!CheckFile(dir, chain)) { 
	delete chain;
	chain = 0;
      }
      break;
    }

    // Save current directory 
    TString savdir(gSystem->WorkingDirectory());
    TSystemDirectory d(gSystem->BaseName(dir.Data()), dir.Data());
    if (!ScanDirectory(&d, chain, pattern, mc, recursive)) { 
      delete chain;
      chain = 0;
    }
    // Go back to the saved directory 
    gSystem->ChangeDirectory(savdir);
  }
  /** 
   * Check if we can add a file to the chain 
   * 
   * @param path   Full path to file 
   * @param chain  Chain 
   * 
   * @return true on success, false otherwise
   */
  Bool_t CheckFile(const TString& path, TChain* chain)
  {
    TFile* test = TFile::Open(path, "READ");
    if (!test) { 
      Warning("CheckFile", "Failed to open %s", path.Data());
      return false;
    }

    Bool_t ok = false; // Assume failure
    TObject* o = test->Get(chain->GetName());
    if (!o) 
      Warning("CheckFile", "The file %s does not contain the object %s", 
	      path.Data(), chain->GetName());
    else if (!dynamic_cast<TTree*>(o)) 
      Warning("CheckFile", "Object %s found in %s is not a TTree", 
	      o->GetName(), path.Data());
    else 
      ok = true;
    test->Close();
    if (ok) chain->AddFile(path);

    return ok;
  }
  /** 
   * Scan directory @a dir (possibly recursive) for tree files to add
   * to the chain.    This does not follow sym-links
   * 
   * @param dir        Directory to scan
   * @param chain      Chain to add to
   * @param fnPattern  File name pattern 
   * @param recursive  Whether to scan recursively 
   *
   * @return true if any files where added 
   */
  Bool_t ScanDirectory(TSystemDirectory* dir, 
		       TChain* chain, 
		       const TString& fnPattern, 
		       Bool_t mc, 
		       Bool_t recursive)
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

    TList toAdd;
    toAdd.SetOwner();
    Bool_t hasGAlice = (!(mc) ? true : false);
    Bool_t hasKine   = (!(mc) ? true : false);
    Bool_t hasTrRef  = (!(mc) ? true : false);
    
    // Sort list of files and check if we should add it 
    files->Sort();
    TIter next(files);
    TSystemFile* file = 0;
    while ((file = static_cast<TSystemFile*>(next()))) {
      TString name(file->GetName());
      TString title(file->GetTitle());
      TString full(gSystem->ConcatFileName(file->GetTitle(), name.Data()));
      if (dynamic_cast<TSystemDirectory*>(file)) full = title;
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
      if (file->IsDirectory(full)) { 
	if (recursive) {
	  // if (title[0] == '/') 
	  TSystemDirectory* d = new TSystemDirectory(file->GetName(),
						     full.Data());
	  if (ScanDirectory(&d, chain, pattern, mc, recursive)) 
	    ret = true;
	  delete d;
	}
        continue;
      }
    
      // If this is not a root file, ignore 
      if (!name.EndsWith(".root")) continue;

      // If this file does not contain AliESDs, ignore 
      if (!name.Contains(fnPattern)) { 
	// Info("ScanDirectory", "%s does not match pattern %s", 
	//      name.Data(), fnPattern.Data());
	if (mc) {
	  if (name.CompareTo("galice.root") == 0)     hasGAlice = true;
	  if (name.CompareTo("Kinematics.root") == 0) hasKine   = true;
	  if (name.CompareTo("TrackRefs.root")  == 0) hasTrRef = true;
	}
	continue;
      }
    
      // Add 
      // Info("ScanDirectory", "Adding %s", full.Data());
      toAdd.Add(new TObjString(full));
    }

    if (mc && toAdd.GetEntries() > 0 && 
	(!hasGAlice || !hasKine || !hasTrRef)) { 
      Warning("ScanDirectory", 
	      "one or more of {galice,Kinematics,TrackRefs}.root missing from "
	      "%s, not adding anything from this directory", 
	      dir->GetTitle());
      toAdd.Delete();
    }

    TIter nextAdd(&toAdd);
    TObjString* s = 0;
    Int_t added = 0;
    while ((s = static_cast<TObjString*>(nextAdd()))) {
      // Info("ScanDirectory", "Adding %s", s->GetString().Data());
      TString fn = s->GetString();
      if (!CheckFile(fn, chain)) continue;

      added++;
    }
    if (added > 0) ret = true;

    gSystem->ChangeDirectory(oldDir);
    return ret;
  }
  //------------------------------------------------------------------
  /** 
   * Create a chain from an XML containing an collection
   * 
   * @param treeName Name of tree's 
   * @param xmlFile  XML collection
   * 
   * @return Newly allocated chain or null
   */
  TChain* CreateChainFromXML(const char* treeName, 
			     const char* xmlFile) 
  {
    Long_t ret = gROOT->ProcessLine(Form("TAlienCollection(\"%s\")", 
					 xmlFile));
    if (!ret) { 
      Error("CreateChainFromXML", "Cannot create AliEn collection from "
	    "XML file %s", xmlFile);
      return 0;
    }
    
    TGridCollection* collection = reinterpret_cast<TGridCollection*>(ret);
    if (!collection) { 
      Error("CreateChainFromXML", "Cannot create AliEn collection from "
	    "XML file %s", xmlFile);
      return 0;
    }

    TChain* chain = new TChain(treeName);
    collection->Reset();
    while (collection->Next()) chain->Add(collection->GetTURL(""));
    
    return chain;
  }

};


//====================================================================
class TrainRemoteWorker : public TrainWorker
{
public:
  TrainRemoteWorker(AliAnalysisGrid* gridHandler) 
    : TrainWorker(gridHandler)
  {}
  virtual void SetROOTVersion(const char* v)
  {
    if (!v || v[0] == '\0') return;
    fHandler->SetROOTVersion(v);
  }
  virtual void SetAliROOTVersion(const char* v)
  {
    if (!v || v[0] == '\0') return;
    fHandler->SetAliROOTVersion(v);
  }
  virtual void SetUsePar(Bool_t usePar) { fUsePar = usePar; }
protected:
  Bool_t fUsePar;
};

//====================================================================
class TrainGridWorker : public TrainRemoteWorker
{
public:
  TrainGridWorker(AliAnalysisGrid* gridHandler) 
    : TrainRemoteWorker(gridHandler)
  {}
  virtual void SetAliEnVersion(const char* v) 
  {
    if (!v || v[0] == '\0') return;
    fHandler->SetAPIVersion(v);
  }
  /** 
   * Specify the input data of the form 
   *
   * @verbatim
   *   [alien://]DIR[?PATTERN]
   * @endverbatim
   * 
   * @param d 
   */
  virtual void SetInput(const char* d)
  {
    if (!d || d[0] == '\0') return;
    TString dir(d);
    
    if (dir.BeginsWith("alien://")) 
      dir.ReplaceAll("alien://", "");

    Int_t question = dir.Index("?");
    if (question != TString::kNPOS) { 
      TString pat = dir(question+1, dir.Length()-question-1);
      fHandler->SetDataPattern(pat);
      dir.Remove(question, dir.Length()-question);
    }
    fHandler->SetGridDataDir(dir);
  }
  virtual Bool_t Connect(const TUrl&, 
			 const TString& name)
  {
    TGrid::Connect("alien://");
    if (!gGrid || !gGrid->IsConnected()) { 
      // This is only fatal in grid mode 
      Error("Connect", "Failed to connect to AliEN");
      return false; 
    }

    // --- Set and make output directory -----------------------------
    // TString name = EscapedName();
    TString homeDir(gGrid->GetHomeDirectory());
    TString workDir(homeDir);
    workDir.Append("/");
    workDir.Append(name);
    
    // Make working directory 
    if (!gGrid->Cd(workDir)) { 
      gGrid->Cd(homeDir);
      if (gGrid->Mkdir(workDir)) {
	gGrid->Cd(name);
	Info("Connect", "Directory %s created", workDir.Data());
      }
    }
    // Make output directory 
    // gGrid->Mkdir("proof_output");
    // gGrid->Cd("proof_output");

    return true;
  }
};

//====================================================================
class TrainProofWorker : public TrainRemoteWorker
{
public:
  TrainProofWorker(AliAnalysisGrid* handler) 
    : TrainRemoteWorker(handler),
      fIsLite(false), 
      fInputDir("")
  {}
  virtual void SetROOTVersion(const char* v)
  {
    if (!v || v[0] == '\0') return;
    TrainRemoteWorker::SetROOTVersion(v);
    fHandler->SetRootVersionForProof(Form("VO_ALICE@ROOT::%s",v));
  }
  virtual void SetInput(const char* d)
  {
    if (!d || d[0] == '\0') return;
    TString inp(d);

    // Check if we got a dataset 
    if (inp.BeginsWith("dataset:", TString::kIgnoreCase)) { 
      Int_t idx colon = inp.Index(":");
      TString src(inp(colon+1, inp.Length()-colon-1));
      fHandler->SetProofDataSet(src);
      return;
    }

    // Check if we got a file list 
    
    // Otherwise, scan for files 
  }
  virtual Bool_t Connect(const TUrl& cluster, 
			 const TString& name)
  {
    TUrl u(cluster);
    if (!u.GetUser() || u.GetUser()[0] == '\0') { 
      UserGroup_t* ugid = gSystem->GetUserInfo();
      u.SetUser(ugid->fUser);
    }
    if (u.GetPort() == 80) u.SetPort(1093);
    TString prot(u.GetProtocol());
    if (prot.Contains("http", TString::kIgnoreCase)) {
      prot.ReplaceAll("http","proof");
      u.SetProtocol(prot);
    }
    TString workers(u.GetValueFromOptions("workers"));
    if (!workers.IsNull()) { 
      Int_t n = workers.Atoi();
      if (workers.EndsWith("x", TString::kIgnoreCase)) 
	fHandler->SetNproofWorkersPerSlave(n);
      else 
	fHandler->SetNproofWorkers(n);
    }
    u.SetOptions("");
    
    fHandler->SetProofCluster(cluster.GetUrl());
  }
private:
  Bool_t fIsLite;
  TString fInputDir;
};
#endif
//
// EOF
//
