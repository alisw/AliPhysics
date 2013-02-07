/**
 * @file   ChainBuilder.C
 * @author Christian Holm Christensen <cholm@master.hehi.nbi.dk>
 * @date   Tue Oct 16 17:54:26 2012
 * 
 * @brief  Build a chain
 * 
 * @ingroup pwglf_forward_trains_util
 */

#ifndef CHAINBUILDER_C
#define CHAINBUILDER_C
#ifndef __CINT__
# include <TString.h>
# include <TChain.h>
# include <TSystemDirectory.h>
# include <TSystem.h>
# include <TFile.h>
# include <TList.h>
# include <TError.h>
# include <TROOT.h>
# include <TGridCollection.h>
# include <TFileCollection.h>
# include <THashList.h>
# include <TKey.h>
# include <fstream>
#else 
class TString;
class TChain;
class TSystemDirectory;
#endif

// ===================================================================
/**
 * Build a chain 
 *
 * @ingroup pwglf_forward_trains_util
 */
struct ChainBuilder 
{
  enum { 
    kInvalid,
    kDirectory, 
    kXML, 
    kAscii, 
    kROOT
  };
  //------------------------------------------------------------------
  static UShort_t CheckSource(TString& src)
  {
    // Local copy 
    TString tmp(src);

    // --- Normalize the path ----------------------------------------
    if (tmp == ".") tmp = "";
    if (!tmp.BeginsWith("/")) tmp.Prepend("../");
    if (gSystem->ExpandPathName(tmp)) { 
      Error("ChainBuilder::CheckSource", 
	    "Failed to expand source %s", src.Data());
      return kInvalid;
    }

    // --- Stat the file ---------------------------------------------
    FileStat_t stat; 
    if (gSystem->GetPathInfo(tmp, stat)) return kInvalid;
    src = tmp;

    // --- Check if directory or file --------------------------------
    if (R_ISDIR(stat.fMode)) return kDirectory;

    // --- check file type -------------------------------------------
    TString type(gSystem->GetFromPipe(Form("file -b %s", src.Data())));
    if      (type.Contains("ROOT"))  return kROOT;
    else if (type.Contains("XML"))   return kXML;
    else if (type.Contains("ASCII")) return kAscii;

    Error("ChainBuilder::CheckSource", 
	  "Do not now how to process %s of type %s", 
	  src.Data(), type.Data());
    return kInvalid;
  }
  //------------------------------------------------------------------
  /** 
   * Create the chain.  User is owner. 
   * 
   * @return Null in case of problems, chain otherwise 
   */
  static TChain* Create(const TString& src, 
			const TString& treeName, 
			const TString& pattern, 
			Bool_t         mc, 
			Bool_t         recursive)
  {
    TString tmp(src);
    UShort_t type = CheckSource(tmp);

    return Create(type, tmp, treeName, pattern, mc, recursive);
  }
  //------------------------------------------------------------------
  /** 
   * Create the chain.  User is owner. 
   * 
   * @return Null in case of problems, chain otherwise 
   */
  static TChain* Create(UShort_t       type, 
			const TString& src, 
			const TString& treeName, 
			const TString& pattern, 
			Bool_t         mc, 
			Bool_t         recursive)
  {
    // --- check input -----------------------------------------------
    if (type == kInvalid) {
      Error("ChainBuilder::Create", "Source %s isn't a file or directory",
	    src.Data());
      return 0;
    }
    TString tN(treeName);
    if (tN.IsNull()) 
      Warning("ChainBuilder::Create", "No tree name specified, assuming T");

    TString pat(pattern);
    if (pat.IsNull()) {
      if      (tN.EqualTo("esdTree")) pat = "AliESD";
      else if (tN.EqualTo("aodTree")) pat = "AliAOD";
    }
      
    // --- Create output ---------------------------------------------
    TChain* chain = new TChain(tN);

    // --- execute based on type 
    Bool_t ret = true;
    switch (type) { 
    case kROOT:      ret = CreateFromFile(chain, src); break;
    case kXML:       ret = CreateFromXML(chain,  src); break;
    case kAscii:     ret = CreateFromList(chain, src); break;
    case kDirectory: ret = CreateFromDirectory(chain, src, 
					       pat, mc, 
					       recursive); break;
    default:         ret = false;
    }

    // --- Clean-up --------------------------------------------------
    if (chain->GetListOfFiles()->GetEntries() <= 0) ret = false;
    if (!ret) { 
      delete chain;
      chain = 0;
    }
    return chain;
  }
  //------------------------------------------------------------------
  /** 
   * Create a chain consiting of a single file 
   * 
   * @param chain The chain
   * @param src File name. 
   * 
   * @return Chain or null
   */
  static Bool_t CreateFromFile(TChain* chain, const TString& src)
  {
    // Info("CreateFromFile", "Making from single file %s", src.Data());
    if (!CheckFile(src, chain)) return false;
    return true;
  }
  //------------------------------------------------------------------
  /** 
   * Create a chain from an XML containing an collection
   * 
   * @return Newly allocated chain or null
   */
  static Bool_t CreateFromXML(TChain* chain, const TString& src) 
  {
    Long_t ret = gROOT->ProcessLine(Form("TAlienCollection(\"%s\")", 
					 src.Data()));
    if (!ret) { 
      Error("ChainBuilder::CreateFromXML", 
	    "Cannot create AliEn collection from XML file %s", src.Data());
      return false;
    }
    
    TGridCollection* collection = reinterpret_cast<TGridCollection*>(ret);
#if 0
    if (!collection) { 
      Error("ChainBuilder::CreateFromXML", 
	    "Cannot create AliEn collection from XML file %s", src.Data());
      return false;
    }
#endif

    collection->Reset();
    while (collection->Next()) chain->Add(collection->GetTURL(""));
    
    return true;
  }
  //------------------------------------------------------------------
  /** 
   * Create a chain from a file containing a list of files
   * 
   * @return Newly allocated chain or null
   */
  static Bool_t CreateFromList(TChain* chain, const TString& src) 
  {
    std::ifstream in(src.Data());
    if (!in) { 
      Error("ChainBuilder::CreateFromList", 
	    "Failed to open list %s", src.Data());
      return false;
    }
    
    while (in.good()) { 
      TString line;
      line.ReadToDelim(in);
      TString l(line.Strip(TString::kBoth));
      if (l.IsWhitespace() || l.BeginsWith("#")) continue;
      
      if (!CheckFile(l, chain))
	Warning("ChainBuilder::CreateFromList", 
		"Failed to add %s to chain", l.Data());
    }
    return true;
  }
  //------------------------------------------------------------------
  /** 
   * Make a chain from a base directory, pattern, and treename -
   * possibly recursively
   * 
   * @return true on success 
   */
  static Bool_t CreateFromDirectory(TChain* chain, 
				    const TString& src, 
				    const TString& pattern, 
				    Bool_t         mc, 
				    Bool_t         recursive) 
  {
    // Info("", "Scanning src=%s, pattern=%s, mc=%d recursive=%d", 
    //       src.Data(), pattern.Data(), mc, recursive);
    // Save current directory 
    TString savdir(gSystem->WorkingDirectory());
    TSystemDirectory d(gSystem->BaseName(src.Data()), src.Data());
    // Info("", "Will scan %s", d.GetTitle());
    if (!ScanDirectory(chain, &d, pattern, mc, recursive)) return false;
    // Go back to the saved directory 
    gSystem->ChangeDirectory(savdir);
    
    return true;
  }
  //------------------------------------------------------------------
  /** 
   * Check if we can add a file to the chain 
   * 
   * @param path   Full path to file 
   * @param chain  Chain 
   * 
   * @return true on success, false otherwise
   */
  static Bool_t CheckFile(const TString& path, TChain* chain)
  {
    // Info("", "Checking %s", path.Data());
    gSystem->RedirectOutput("/dev/null", "w");
    TFile* test = TFile::Open(path, "READ");
    gSystem->RedirectOutput(0);
    if (!test) { 
      Warning("ChainBuilder::CheckFile", "Failed to open %s", path.Data());
      return false;
    }

    TObject* o = test->Get(chain->GetName());
    if (!o) {
      // Let's try to find a TFileCollection 
      TList* l = test->GetListOfKeys();
      TIter next(l);
      TKey* k = 0;
      Bool_t ok = false;
      while ((k = static_cast<TKey*>(next()))) {
	TString cl(k->GetClassName());
	if (!cl.EqualTo("TFileCollection")) continue;
	TFileCollection* fc = dynamic_cast<TFileCollection*>(k->ReadObj());
	if (!fc) { 
	  Warning("", "Returned collection invalid");
	  continue;
	}
	Info("", "Adding file collection");
	chain->AddFileInfoList(fc->GetList());
	ok = true;
      }
      if (ok) { 
	test->Close();
	return true;
      }
    }
    else if (dynamic_cast<TTree*>(o)) {
      test->Close();
      chain->Add(path);
      return true;
    }
    
    Warning("ChainBuilder::CheckFile", 
	    "The file %s does not contain the tree %s or a file collection", 
	    path.Data(), chain->GetName());
    
    return false;
  }
  //------------------------------------------------------------------
  /** 
   * Scan directory @a dir (possibly recursive) for tree files to add
   * to the chain.    This does not follow sym-links
   * 
   * @param dir        Directory to scan
   * @param chain      Chain to add to
   * @param pattern    File name pattern 
   * @param mc         Simulation input 
   * @param recursive  Scan recursive 
   *
   * @return true if any files where added 
   */
  static Bool_t ScanDirectory(TChain*           chain, 
			      TSystemDirectory* dir,
			      const TString&    pattern,
			      Bool_t            mc,
			      Bool_t            recursive)
  {
    // Assume failure 
    Bool_t ret = false;

    // Get list of files, and go back to old working directory
    TString oldDir(gSystem->WorkingDirectory());
    TList*  files = dir->GetListOfFiles();
    if (!gSystem->ChangeDirectory(oldDir)) { 
      Error("ChainBuilder::ScanDirectory", "Failed to go back to %s", 
	    oldDir.Data());
      return false;
    }
    if (!files) {
      Warning("ChainBuilder::ScanDirectory", "No files");
      return false;
    }

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
      // Info("", "Got file %s", full.Data());
      if (file->IsA()->InheritsFrom(TSystemDirectory::Class())) full = title;
      // Ignore special links 
      if (name == "." || name == "..") { 
	// Info("ChainBuilder::ScanDirectory", "Ignoring %s", name.Data());
	continue;
      }
      // Info("", "Got file %s", full.Data());

      FileStat_t fs;
      if (gSystem->GetPathInfo(full.Data(), fs)) {
	Warning("ChainBuilder::ScanDirectory", "Cannot stat %s (%s)", 
		full.Data(), gSystem->WorkingDirectory());
	continue;
      }
      // Check if this is a directory 
      if (file->IsDirectory(full)) { 
	// Info("", "Recursive scan of %s", full.Data());
	if (recursive) {
	  // if (title[0] == '/') 
	  TSystemDirectory* d = new TSystemDirectory(file->GetName(),
						     full.Data());
	  if (ScanDirectory(chain, d, pattern, mc, recursive)) 
	    ret = true;
	  delete d;
	}
        continue;
      }
    
      // If this is not a root file, ignore 
      if (!name.EndsWith(".root")) {
	// Info("ScanDirectory", "File %s does not end in .root", name.Data());
	continue;
      }

      // If this file does not contain AliESDs, ignore 
      if (!name.Contains(pattern)) { 
	// Info("ChainBuilder::ScanDirectory", "%s does not match pattern %s", 
	//      name.Data(), pattern.Data());
	if (mc) {
	  if (name.CompareTo("galice.root") == 0)     hasGAlice = true;
	  if (name.CompareTo("Kinematics.root") == 0) hasKine   = true;
	  if (name.CompareTo("TrackRefs.root")  == 0) hasTrRef = true;
	}
	continue;
      }
    
      // Add 
      // Info("ChainBuilder::ScanDirectory", "Adding %s", full.Data());
      toAdd.Add(new TObjString(full));
    }

    if (mc && toAdd.GetEntries() > 0 && 
	(!hasGAlice || !hasKine || !hasTrRef)) { 
      Warning("ChainBuilder::ScanDirectory", 
	      "one or more of {galice,Kinematics,TrackRefs}.root missing from "
	      "%s, not adding anything from this directory", 
	      dir->GetTitle());
      toAdd.Delete();
    }

    TIter nextAdd(&toAdd);
    TObjString* s = 0;
    Int_t added = 0;
    while ((s = static_cast<TObjString*>(nextAdd()))) {
      // Info("ChainBuilder::ScanDirectory", 
      //      "Adding %s", s->GetString().Data());
      TString fn = s->GetString();
      if (!CheckFile(fn, chain)) continue;

      added++;
    }
    if (added > 0) ret = true;

    gSystem->ChangeDirectory(oldDir);
    return ret;
  }
};
#endif
