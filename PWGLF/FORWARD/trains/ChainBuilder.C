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
# include <TUrl.h>
# include <TString.h>
# include <TChain.h>
# include <TChainElement.h>
# include <TSystemDirectory.h>
# include <TSystem.h>
# include <TFile.h>
# include <TFileInfo.h>
# include <TList.h>
# include <TError.h>
# include <TROOT.h>
# include <TGridCollection.h>
# include <TFileCollection.h>
# include <THashList.h>
# include <TKey.h>
# include <TRegexp.h>
# include <fstream>
#else 
class TString;
class TChain;
class TSystemDirectory;
class TUrl;
class TFileCollection;
#endif

//====================================================================  
struct SuppressGuard
{
  Int_t save = 0;
  SuppressGuard(Int_t lvl=2000)
  {
    save = gErrorIgnoreLevel;
    gErrorIgnoreLevel = lvl;
  }
  ~SuppressGuard()
  {
    gErrorIgnoreLevel = save;
  }
};

// ===================================================================
/**
 * Build a chain 
 *
 * @ingroup pwglf_forward_trains_util
 */
struct ChainBuilder 
{
  enum {
    kVerbose   =  0x1,
    kRecursive =  0x2,
    kMC        =  0x4,
    kCheck     =  0x8,
    kClean     = 0x10,
    kScan      = 0x20, 
    kTrRef     = 0x40,
    kRemote    = 0x80
  };
  enum { 
    kInvalid,
    kDirectory, 
    kXML, 
    kAscii, 
    kROOT,
    kZip
  };
  //------------------------------------------------------------------
  static UShort_t CheckSource(TString& src, UShort_t flags)
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
    if (R_ISDIR(stat.fMode)) {
      if (flags & kVerbose) 
	Info("ChainBuilder::CheckSource", "%s is a directory", tmp.Data());
      return kDirectory;
    }

    // --- check file type -------------------------------------------
    TString type(gSystem->GetFromPipe(Form("file -b -L %s", src.Data())));
    if ((flags & kVerbose))
      Info("ChainBuilder::CheckSource", "file -b %s -> %s", 
	   tmp.Data(), type.Data());
    UShort_t ret = kInvalid;
    if      (type.Contains("ROOT"))  ret = kROOT;
    else if (type.Contains("XML"))   ret = kXML;
    else if (type.Contains("ASCII")) ret = kAscii;
    else if (type.Contains("Zip"))   ret = kZip;
    
    if (ret == kInvalid) {
      Error("ChainBuilder::CheckSource", 
	    "Do not now how to process %s of type %s", 
	    src.Data(), type.Data());
    }
    return ret;
  }
  //------------------------------------------------------------------
  /** 
   * Create a TChain from a URL specification.
   *
   * The URL should have the format 
   * 
   * @verbatim
   *   PROTOCOL://PATH?OPTIONS#TREENAME
   * @endverbatim
   *
   * where 
   * - @c PROTOCOL is any protocol
   * - @c PATH is a file, top-directory, file containing a
   *   TFileCollection, ASCII file with list of files, etc.
   * - @c OPTIONS is a list of options, separate by @c &
   * - @c TREENAME is the tree name 
   *
   * @c OPTIONS can be one or more of 
   * - @c mc  Also check for auxiliary MC files
   * - @c recursive When scanning directories, do so recursively
   * - @c verbose Be verbose
   * - @c check Check files by trying to open them 
   * - @c clean Remove invalid files 
   * - @c trackref For MC input, insist on TrackRefs.root presence
   * - @c pattern=PATTERN Search pattern when scanning directories 
   * 
   * @param url The input url 
   * 
   * @return Pointer to newly allocated TChain object or null
   */
  static TChain* Create(const TUrl& url)
  {
    TString     source   = url.GetFile();
    TString     treeName = url.GetAnchor();
    TString     pattern  = "";
    UShort_t    flags    = 0;
    TString     options  = url.GetOptions();
    TObjArray*  tokens   = options.Tokenize("&");
    TObjString* token    = 0; 
    TIter       next(tokens);
    while ((token = static_cast<TObjString*>(next()))) {
      const TString& str = token->String();
      TString lstr(str); lstr.ToLower();
      if      (lstr.EqualTo("mc")) 	 flags |= kMC;
      else if (lstr.EqualTo("recursive"))flags |= kRecursive;
      else if (lstr.EqualTo("verbose"))  flags |= kVerbose;
      else if (lstr.EqualTo("check"))    flags |= kCheck;
      else if (lstr.EqualTo("trackref")) flags |= kTrRef;
      else if (lstr.EqualTo("clean"))    flags |= kClean; 
      else if (lstr.EqualTo("scan"))     flags |= kScan;
      else if (lstr.EqualTo("remote"))   flags |= kRemote;
      else if (lstr.BeginsWith("pattern=")) { 
	Int_t eq = str.Index("=");
	pattern  = str(eq+1, str.Length()-eq-1);
	pattern.ReplaceAll("@", "#");
	pattern.ReplaceAll(":", "#");
      }
      else 
	Warning("", "Option %s unknown", str.Data());
    }
    delete tokens;

    TString tmp(source);
    UShort_t type = CheckSource(tmp, flags);
    
    return Create(type, tmp, treeName, pattern, flags);
  }
  //------------------------------------------------------------------
  /** 
   * Create a chain
   * 
   * @param src          Source 
   * @param treeName     Tree name 
   * @param pattern      Pattern for scans
   * @param mc           If true, check for MC files
   * @param recursive    If true, scan recursively
   * @param verbose      If true, be verbose
   * @param checkFiles   If true, check that files can be opened
   * @param removeFiles  If true, remove bad files 
   * @param trackRefs    If true, look for track references too 
   * @param remote       For remote access 
   * 
   * @return Pointer to newly allocated TChain or null
   */
  static TChain* Create(const TString& src, 
			const TString& treeName, 
			const TString& pattern, 
			Bool_t         mc, 
			Bool_t         recursive,
			Bool_t         verbose=false,
			Bool_t         checkFiles=false, 
			Bool_t         removeFiles=false,
			Bool_t         trackRefs=false,
			Bool_t         remote=false)
  {
    UShort_t flags = 0;
    if (verbose)     flags |= kVerbose;
    if (recursive)   flags |= kRecursive;
    if (mc)          flags |= kMC;
    if (checkFiles)  flags |= kCheck;
    if (removeFiles) flags |= kClean;
    if (trackRefs)   flags |= kTrRef;
    if (remote)      flags |= kRemote;

    TString tmp(src);
    UShort_t type = CheckSource(tmp, flags);

    return Create(type, tmp, treeName, pattern, flags);
  }
  //------------------------------------------------------------------
  /** 
   * Create a chain
   * 
   * @param type         Type of input
   * @param src          Source 
   * @param treeName     Tree name 
   * @param pattern      Pattern for scans
   * @param mc           If true, check for MC files
   * @param recursive    If true, scan recursively
   * @param verbose      If true, be verbose
   * @param checkFiles   If true, check that files can be opened
   * @param removeFiles  If true, remove bad files 
   * @param trackRefs    If true, look for track references too 
   * @param remote       For remote access 
   * 
   * @return Pointer to newly allocated TChain or null
   */
  static TChain* Create(UShort_t       type,
			const TString& src, 
			const TString& treeName, 
			const TString& pattern, 
			Bool_t         mc, 
			Bool_t         recursive,
			Bool_t         verbose=false,
			Bool_t         checkFiles=false, 
			Bool_t         removeFiles=false,
			Bool_t         trackRefs=false,
			Bool_t         remote=false)
  {
    // Info("ChainBuilder::Create", 
    // "src=%s treeName=%s pattern=%s mc=%s recursive=%s",
    // src.Data(), treeName.Data(), pattern.Data(), 
    // (mc ? "true" : "false"), (recursive ? "true" : "false"));
    // --- Create flags 
    UShort_t flags = 0;
    if (verbose)     flags |= kVerbose;
    if (recursive)   flags |= kRecursive;
    if (mc)          flags |= kMC;
    if (checkFiles)  flags |= kCheck;
    if (removeFiles) flags |= kClean;
    if (trackRefs)   flags |= kTrRef;
    if (remote)      flags |= kRemote;
    

    return Create(type, src, treeName, pattern, flags);
  }
  //------------------------------------------------------------------
  /** 
   * Create a chain from the inputs
   * 
   * @param type        Type of input
   * @param src         Source 
   * @param treeName    Tree name
   * @param pattern     Pattern for scans 
   * @param flags       Flags 
   * 
   * @return Pointer to newly allocated TChain object or null
   */
  static TChain* Create(UShort_t       type, 
			const TString& src, 
			const TString& treeName, 
			const TString& pattern,
			UShort_t       flags)
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
      if      (tN.EqualTo("esdTree")) pat = "AliESD*";
      else if (tN.EqualTo("aodTree")) pat = "AliAOD*";
      else if (tN.EqualTo("TE"))      pat = "galice*";
      if ((flags & kVerbose)) Info("", "Pattern set to %s", pat.Data());
    }
    if ((flags & kVerbose))
      Info("ChainBuilder::Create", "Type=%s, tree=%s, pattern=%s", 
      (type == kDirectory ? "directory" : 
      type == kXML       ? "XML" : 
      type == kAscii     ? "ASCII" : 
      type == kROOT      ? "ROOT" : "unknown"),
      tN.Data(), pat.Data());
    
    // --- Create output ---------------------------------------------
    TChain* chain = new TChain(tN);

    // --- ZIP archives ----------------------------------------------
    TString anchor;
    TString tmp(pat);
    ExtractAnchor(pat, anchor);
    if ((flags & kVerbose)) 
      Info("", "Full pattern: '%s' filename pattern: '%s' anchor: '%s'",
	   tmp.Data(), pat.Data(), anchor.Data());

    // --- execute based on type 
    Bool_t ret = true;
    switch (type) { 
    case kROOT:      ret = CreateFromFile(chain, src, anchor, flags); break;
    case kXML:       ret = CreateFromXML(chain,  src); break;
    case kAscii:     ret = CreateFromList(chain, src); break;
    case kDirectory: ret = CreateFromDirectory(chain, src, 
					       pat, anchor, 
					       flags); break;
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
   * Create a collection  
   * 
   * @param output Output file 
   * @param url    Input url 
   * @param remote For remote access 
   */
  static void CreateCollection(const TString& output, 
			       const TUrl&    url,
			       const char*    remote=0)
  {
    TChain* chain = Create(url);
    if (!chain) return;

    CreateCollection(output, chain, remote);
  }
  //------------------------------------------------------------------
  /** 
   * Create a collection 
   * 
   * @param output Input url 
   * @param chain  Chain to make collection from 
   * @param remote For remote access 
   */
  static void CreateCollection(const TString& output, 
			       const TChain*  chain,
			       const char*    remote=0)
  {
    if (!chain) return;
    TDirectory* savDir = gDirectory;
    TFile* out = TFile::Open(output, "RECREATE");
    if (!out) { 
      Error("", "Failed to open %s for output", output.Data());
      return;
    }
    TFileCollection* collection = new TFileCollection(chain->GetName());
    TObjArray*       files      = chain->GetListOfFiles();
    TChainElement*   element    = 0;
    TIter            next(files);


    collection->SetDefaultTreeName(chain->GetName());
    Long64_t nEntries = 0;
    while ((element = static_cast<TChainElement*>(next()))) {
      Info("", "Element: '%s' - '%s' %lld", 
	   element->GetName(), element->GetTitle(), 
	   element->GetEntries());
      TFileInfo*     info = new TFileInfo(element->GetTitle());
      TFileInfoMeta* meta = new TFileInfoMeta(Form("/%s",element->GetName()),
					      "TTree", element->GetEntries());
      info->AddMetaData(meta);
      info->SetBit(TFileInfo::kStaged);
      // info->AddUrl(Form("file://%s", element->GetTitle()));
      collection->Add(info);
      
      Long64_t n = element->GetEntries();
      if (n >= 0) nEntries += n;
      
    }
    Remotify(collection, remote);
    collection->Update();
    TFileInfoMeta* cMeta = new TFileInfoMeta(chain->GetName(), 
					     "TTree", nEntries);
    collection->AddMetaData(cMeta);
    out->cd();
    collection->Write();
    Printf("A total of %lld entries", nEntries);
    // collection->Print("MFL");
    out->Close();
    savDir->cd();
  }
  //------------------------------------------------------------------
  /** 
   * Exrtact the anchor 
   * 
   * @param src    Source url 
   * @param anchor On return, contains the anchor
   */
  static void ExtractAnchor(TString& src, TString& anchor)
  {
    anchor         = "";
    Int_t idxHash  = src.Index("#");
    
    if (idxHash == kNPOS) return;

    TString tmp = src(0,idxHash);
    anchor      = src(idxHash+1, src.Length()-idxHash-1);
    src = tmp;
  }

  //------------------------------------------------------------------
  /** 
   * Create a chain consiting of a single file 
   * 
   * @param chain The chain
   * @param anchor Anchor (tree name)
   * @param src File name. 
   * @param flags       Flags 
   * 
   * @return Chain or null
   */
  static Bool_t CreateFromFile(TChain*        chain, 
			       const TString& src, 
			       const TString& anchor,
			       UShort_t       flags=0)
  {
    // Info("CreateFromFile", "Making from single file %s", src.Data());
    
    if (!CheckFile(src, anchor, chain, flags)) return false;
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
    // Info("ChainBuilder::CreateFromXML", "Create from XML");
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
  static Bool_t CreateFromList(TChain*        chain, 
			       const TString& src,
			       UShort_t       flags=0)
  {
    // Info("ChainBuilder::CreateFromList", "Creating from list");
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
      
      TString anchor;
      ExtractAnchor(l, anchor);
      if (!CheckFile(l, anchor, chain, flags))
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
				    const TString& anchor,
				    UShort_t       flags)
  {
    // Info("", "Scanning src=%s, pattern=%s, mc=%d recursive=%d", 
    //      src.Data(), pattern.Data(), mc, recursive);
    // Save current directory 
    TString savdir(gSystem->WorkingDirectory());
    TSystemDirectory d(gSystem->BaseName(src.Data()), src.Data());
    if (flags & kVerbose) Info("", "Will scan %s", d.GetTitle());
    if (!ScanDirectory(chain, &d, pattern, anchor, flags))
      return false;
    // Go back to the saved directory 
    gSystem->ChangeDirectory(savdir);
    
    return true;
  }
  static void RemoveFile(const TString& path)
  {
    Info("", "Removing bad file %s", path.Data());
    SuppressGuard g;
    // gSystem->Unlink(path);
    gSystem->Rename(path, Form("%s.bad", path.Data()));
  }
  //------------------------------------------------------------------
  static TFileCollection* Remotify(TFileCollection* fc, const char* remote=0)
  {
    if (!remote) return fc;
    
    TList*           files      = fc->GetList();
    TFileInfo*       element    = 0;
    TIter            next(files);
    while ((element = static_cast<TFileInfo*>(next()))) {
      element->ResetUrl();
      TUrl* url  = 0;
      TUrl* furl = 0;
      while ((url = element->NextUrl())) {
	if (TString(url->GetProtocol()).BeginsWith("root")) {
	  furl = 0;
	  break;
	}
	if (TString(url->GetProtocol()).BeginsWith("file") && furl == 0)
	  furl = url;	  
      }
      if (furl) {
	TUrl* nurl = static_cast<TUrl*>(furl->Clone());
	nurl->SetProtocol("rootd");
	nurl->SetHost(remote);
	element->AddUrl(nurl->GetUrl(), true);
      }
    }
    return fc;
  }
  
  //------------------------------------------------------------------
  /** 
   * Check if we can add a file to the chain 
   * 
   * @param path   Full path to file 
   * @param anchor Anchor (tree name)
   * @param chain  Chain 
   * @param flags  Some flags
   * 
   * @return true on success, false otherwise
   */
  static Bool_t CheckFile(const TString& path, 
			  const TString& anchor, 
			  TChain*        chain,
			  UShort_t       flags=0)
  {
    if (flags & kVerbose) Info("", "Checking %s", path.Data());
    TString fn   = path;
    if (!anchor.IsNull()) fn.Append(TString::Format("#%s", anchor.Data()));

    TFile* test = 0;
    {
      SuppressGuard g((flags & kVerbose) ? 0 : 2000);
      test = TFile::Open(fn, "READ");
    }
    if (!test) { 
      Warning("ChainBuilder::CheckFile", "Failed to open %s", fn.Data());
      if (flags & kClean) RemoveFile(path);
      return false;
    }
    
    Bool_t           ok = false;
    TObject*         o  = test->Get(chain->GetName());
    TTree*           t  = dynamic_cast<TTree*>(o);
    TFileCollection* c  = dynamic_cast<TFileCollection*>(o);
    if (t) {
      test->Close();
      ok = true;
      if (flags & kMC) { 
	const char*  auxs[] = { "galice", "Kinematics",
				(flags & kTrRef ? "TrackRefs" : 0),
				0 };
	const char** aux    = auxs;
	while ((*aux)) { 
	  TString t1;
	  if (anchor.IsNull()) 
	    t1 = gSystem->ConcatFileName(gSystem->DirName(path.Data()),
					 Form("%s.root", *aux));
	  else 
	    t1 = TString::Format("%s#%s.root", path.Data(), *aux);
	  TFile* t2 = 0;
	  {
	    SuppressGuard g2;
	    t2 = TFile::Open(t1, "READ");
	  }
	  if (!t2) { 
	    Error("", "Needed MC file %s not found", t1.Data());
	    ok = false;
	    break;
	  }
	  t2->Close();
	  aux++;
	}
      }
      // if (flags & kRemote)
      //   fn.Prepend(Form("root://%s/", gSystem->HostName()));      
      if (ok) chain->Add(fn, kScan ? -1 : TChain::kBigNumber);
    } else if (c) {
      // chain->AddFileInfoList(Remotify(c, flags & kRemote)->GetList());
      chain->AddFileInfoList(c->GetList());
      ok = true;
    } else {
      // Let's try to find a TFileCollection 
      TList* l = test->GetListOfKeys();
      TIter next(l);
      TKey* k = 0;
      while ((k = static_cast<TKey*>(next()))) {
	TString cl(k->GetClassName());
	if (!cl.EqualTo("TFileCollection")) continue;
	c = dynamic_cast<TFileCollection*>(k->ReadObj());
	if (!c) { 
	  Warning("", "Returned collection invalid");
	  continue;
	}
	// Info("", "Adding file collection");
	// chain->AddFileInfoList(Remotify(c, flags&kRemote)->GetList());
	chain->AddFileInfoList(c->GetList());
	ok = true;
      }
      test->Close();
    }

    if (!ok) {
      Warning("ChainBuilder::CheckFile", 
	      "The file %s does not contain the tree %s or a file collection", 
	      path.Data(), chain->GetName());
      if (flags & kClean) RemoveFile(path);
    }
    return ok;
  }
  //------------------------------------------------------------------
  /** 
   * Scan directory @a dir (possibly recursive) for tree files to add
   * to the chain.    This does not follow sym-links
   * 
   * @param dir        Directory to scan
   * @param chain      Chain to add to
   * @param pattern    File name pattern 
   * @param anchor     Anchor (tree name)
   * @param flags      Flags
   *
   * @return true if any files where added 
   */
  static Bool_t ScanDirectory(TChain*           chain, 
			      TSystemDirectory* dir,
			      const TString&    pattern,
			      const TString&    anchor,
			      UShort_t          flags)
  {
    // Assume failure 
    Bool_t ret = false;
    TRegexp wild(pattern, true);

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
    
    // Sort list of files and check if we should add it 
    files->Sort();
    TIter next(files);
    TSystemFile* file = 0;
    while ((file = static_cast<TSystemFile*>(next()))) {
      TString name(file->GetName());
      TString title(file->GetTitle());
      TString full(gSystem->ConcatFileName(file->GetTitle(), name.Data()));
      if (file->IsA()->InheritsFrom(TSystemDirectory::Class())) full = title;
      // Ignore special links 
      if (name == "." || name == "..") { 
	// Info("ChainBuilder::ScanDirectory", "Ignoring %s", name.Data());
	continue;
      }
      if ((flags & kVerbose)) Info("", "Got file %s", full.Data());

      FileStat_t fs;
      if (gSystem->GetPathInfo(full.Data(), fs)) {
	Warning("ChainBuilder::ScanDirectory", "Cannot stat %s (%s)", 
		full.Data(), gSystem->WorkingDirectory());
	continue;
      }
      // Check if this is a directory 
      if (file->IsDirectory(full)) { 
	if ((flags & kVerbose)) Info("", "Recursive scan of %s", full.Data());
	if ((flags & kRecursive)) {
	  // if (title[0] == '/') 
	  TSystemDirectory* d = new TSystemDirectory(file->GetName(),
						     full.Data());
	  if (ScanDirectory(chain, d, pattern, anchor, flags))
	    ret = true;
	  delete d;
	}
        continue;
      }
    
      // If this is not a root file, ignore 
      if (!name.EndsWith(".root") && !name.EndsWith(".zip")) {
	if ((flags & kVerbose))
	  Info("ScanDirectory", "File %s does not end in .root/.zip", 
	       name.Data());
	continue;
      }

      // If this file does not contain AliESDs, ignore 
      if (!name.Contains(wild)) { 
	if ((flags & kVerbose))
	  Info("ChainBuilder::ScanDirectory", 
	       "%s does not match pattern %s", 
	       name.Data(), pattern.Data());
	continue;
      }
    
      // Add 
      // Info("ChainBuilder::ScanDirectory", "Adding %s", full.Data());
      toAdd.Add(new TObjString(full));
    }

    TIter nextAdd(&toAdd);
    TObjString* s = 0;
    Int_t added = 0;
    while ((s = static_cast<TObjString*>(nextAdd()))) {
      // Info("ChainBuilder::ScanDirectory", 
      //      "Adding %s", s->GetString().Data());
      TString fn = s->GetString();
      if (!CheckFile(fn, anchor, chain, flags)) continue;

      added++;
    }
    if (added > 0) ret = true;

    gSystem->ChangeDirectory(oldDir);
    return ret;
  }
};
#endif
//
// EOF
//
