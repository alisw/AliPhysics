#ifndef __CINT__
#include <TSystem.h>
#include <TSystemDirectory.h>
#include <TError.h>
#include <TString.h>
#include <TList.h>
#include <TFile.h>
#include <TROOT.h>
#include <TH1.h>
#include <TParameter.h>
#include "AliOADBForward.h"
#include <AliForwardCorrectionManager.h>
#include <AliCentralCorrectionManager.h>
#include <AliCorrectionManagerBase.h>
#else
class AliOADBForward;
class TSystemDirectory;
class AliCorrectionManagerBase;
class AliForwardCorrectionManager;
class AliCentralCorrectionManager;
class TFile;
#endif

struct Scanner
{
  /** 
   * Check if a path points to a file 
   * 
   * @param path Path
   * 
   * @return True if the path points to a regular file 
   *
   * @ingroup pwglf_forward_scripts
   */
  static Bool_t IsFile(const char* path)
  {
    FileStat_t stat; 
    gSystem->GetPathInfo(path, stat);
    if (stat.fIsLink || 
	R_ISDIR(stat.fMode) || 
	!R_ISREG(stat.fMode)) return false;
    return true;
    
#if 0
    Long_t id;
    Long_t size;
    Long_t flags;
    Long_t modtime;
    gSystem->GetPathInfo(path, &id, &size, &flags, &modtime);
    return !((flags & 0x2) == 0x2);
#endif
  }
  /** 
   * Test if we can open a file 
   * 
   * @param name    Name of file 
   * @param pattern Pattern to check against 
   * 
   * @return True on success
   */
  static Bool_t TestFile(const TString& name, const char* pattern=0)
  {
    // If this is not a root file, ignore 
    if (!name.EndsWith(".root")) return false;
    
    // If this file does not contain the pattern, ignore 
    if (pattern && pattern[0] != '\0' && !name.Contains(pattern)) return false;
    
    Bool_t ret  = true;
    TFile* test = TFile::Open(name.Data(), "READ");
    if (!test || test->IsZombie()) { 
      ::Warning("TestFile", "Failed to open file %s", name.Data());
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
   * @param list       List to add data to 
   * @param pattern    Pattern that the file name must contain
   * @param recursive  Whether to scan recursively 
   *
   * @ingroup pwglf_forward_scripts
   */
  void
  ScanDirectory(TSystemDirectory* dir, TList* list, 
		const char* pattern, bool recursive)
  {
    // Get list of files, and go back to old working directory
    TString oldDir(gSystem->WorkingDirectory());
    TList* files = dir->GetListOfFiles();
    if (!files) return;
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
	  ScanDirectory(static_cast<TSystemDirectory*>(file),list,
			pattern,recursive);
	continue;
      }
    
      // Get the path 
      TString data(Form("%s/%s", file->GetTitle(), name.Data()));

      // Check the fuile 
      if (!TestFile(data, pattern)) continue;
      list->Add(new TObjString(data));
    }
  }
  /** 
   * Scan data directory for files matching pattern 
   * 
   * @param datadir   Path to data directory 
   * @param pattern   File name match pattern 
   * @param recursive Recurse flag
   * 
   * @return List of file names
   */
  TList* Scan(const char* datadir, const char* pattern, bool recursive=false) 
  {
    // --- Get list of files --------------------------------------------
    // Open source directory, and make sure we go back to were we were 
    TString oldDir(gSystem->WorkingDirectory());
    TString path(gSystem->ExpandPathName(datadir));
    if (IsFile(path)) {
      Error("Scan", "%s is a file", datadir);
      return 0;
    }
    Info("Scan", "Scanning %s", path.Data());

    TList* ret = new TList;
    ret->SetOwner();
    TSystemDirectory d(datadir, path.Data());
    ScanDirectory(&d, ret, pattern, recursive);

    // Make sure we do not make an empty chain 
    if (ret->GetEntries() <= 0) { 
      Warning("Scane", "list is empty for input %s, %s", 
	      datadir, pattern);
      delete ret;
      ret = 0;
    }
    return ret;
  }
};


struct Extractor
{
  /** 
   * Constructor 
   * 
   * @param dirName   Directory name 
   * @param corrName  Correction name 
   * @param methName  Run number mode to set as default  
   * @param outFile   Output file 
   * @param cm        Correction manager to use 
   */
  Extractor(const char* dirName, 
	    const char* corrName,
	    const char* methName,
	    const char* outFile,
	    AliCorrectionManagerBase* cm)
    : fDirName(dirName), fCorrName(corrName), fMethName(methName), 
      fFile(outFile), fCM(cm)
  {
  }
  virtual ~Extractor() {}
  /** 
   * Extract files 
   * 
   * @return number of converted objects 
   */
  virtual Int_t Extract(const char* prefix="$ALICE_PHYSICS/PWGLF/FORWARD/corrections")
  {
    Scanner s;
    TString dir = TString::Format("%s/%s",
				  prefix, fDirName.Data());
    TList* l = s.Scan(dir, fCorrName);
    if (!l) {
      Warning("Extract", "No files matching %s found in %s",  
	      fCorrName.Data(), dir.Data()); 
      return 0;
    }

    TIter next(l);
    TObjString* os = 0;
    Int_t ret = 0;
    while ((os = static_cast<TObjString*>(next()))) { 
      TString& fn  = os->String();
      if (ExtractFile(fn)) ret++;
    }
    return ret;
  }
  /** 
   * Extract from a file 
   * 
   * @param fn File name 
   * 
   * @return true on success 
   */      
  virtual Bool_t ExtractFile(const TString& fn) 
  {
    UShort_t sys = 0;
    UShort_t sNN = 0;
    Short_t  fld = 0;
    Bool_t   mc  = false;
    
    ExtractFields(fn, sys, sNN, fld, mc);
    if (sNN == 2750) sNN = 2760;

    ULong_t runNo = ExtractRunNo(sys, sNN);
    if (runNo == 0xFFFFFFFF || runNo <= 0) return false;

    TObject* obj = ExtractObject(fn.Data());
    if (!obj) return false;

    // Run, sys, sNN, fld, mc, sat, obj, full, meth
    Info("", "File %s to be stored: run=%d sys=%d sNN=%d fld=%d mc=%d", 
	 fn.Data(), runNo, sys, sNN, fld, mc);
    return fCM->Store(obj, runNo, sys, sNN, fld, mc, 
		      false, fFile, fMethName);
  }
  /** 
   * Extract fields from file name 
   * 
   * @param s    File name 
   * @param sys  System 
   * @param sNN  Energy 
   * @param fld  Field 
   * @param mc   MC flag
   */    
  virtual void ExtractFields(const TString& s, 
			     UShort_t& sys, 
			     UShort_t& sNN, 
			     Short_t&  fld, 
			     Bool_t&   mc) 
  {
    TString    str    = gSystem->BaseName(s.Data());
    str.ReplaceAll(".root", "");
    TObjArray* tokens = str.Tokenize("_");
    // tokens->ls();

    TString&   sSys   = ((TObjString*)(tokens->At(1)))->String();
    TString&   sSNN   = ((TObjString*)(tokens->At(2)))->String();

    if      (sSys.EqualTo("pbpb", TString::kIgnoreCase))     sys = 2;
    else if (sSys.EqualTo("ppb",  TString::kIgnoreCase))     sys = 3;
    else if (sSys.EqualTo("pp",   TString::kIgnoreCase))     sys = 1;
    
    sSNN.ReplaceAll("GeV", "");
    Info("", "sSNN=%s -> ", sSNN.Data());
    while (sSNN[0] == '0' && sSNN.Length() > 1) sSNN.Remove(0, 1);
    sNN = sSNN.Atoi();
    Info("", "sSNN=%s sNN=%d", sSNN.Data(), sNN);

    if (tokens->GetEntries() > 3) {
      TString&   sFld   = ((TObjString*)(tokens->At(3)))->String();
      sFld.ReplaceAll("kG", "");
      while (sFld[0] == '0' && sFld.Length() > 1) sFld.Remove(0, 1);
      sFld.ReplaceAll("p", "+");
      sFld.ReplaceAll("m", "-");
      fld = sFld.Atoi();
    }

    if (tokens->GetEntries() > 4) { 
      TString& sMC = ((TObjString*)(tokens->At(4)))->String();
      mc = sMC.EqualTo("mc", TString::kIgnoreCase);
    }
    tokens->Delete();
    // Info("Extract", "%s -> %d %d %d %d", str.Data(), sys, sNN, fld, mc);
  }
  /** 
   * Get run number corresponding to arguments 
   * 
   * @param sys  System
   * @param sNN  Energy
   * 
   * @return run number 
   */
  virtual ULong_t ExtractRunNo(UShort_t sys, UShort_t sNN)
  {
    ULong_t run = 0;
    switch (sys) { 
    case 1: // pp 
      switch (sNN) { 
      case 900:   run = 118502; break;
      case 2760:  run = 146686; break;
      case 7000:  run = 114747; break;
      case 14000: run = 0xFFFFFFFF; break;
      }
      break;
    case 2: // PbPb 
      switch (sNN) { 
      case 2760: run = 137123; break;
      }
      break;
    case 3: // pPb 
      switch (sNN) {
      case 5023: run = 188246; break;
      }
      break;
    }
    if (run == 0) 
      Warning("ExtractRunNo", 
	      "Unknown energy %d for collision system %d", sNN, sys);
    return run;
  }
  /** 
   * Extract a single object from the file 
   * 
   * @param fn  File name 
   * 
   * @return Object or null
   */
  virtual TObject* ExtractObject(const TString& fn)
  {
    TFile* file = TFile::Open(fn.Data(), "READ");
    if (!file) { 
      Error("ExtractObject", "Failed to open %s", fn.Data());
      return 0;
    }
    // file->ls();
    
    TObject* obj =  file->Get(fCorrName);
    if (!obj) { 
      Error("ExtractObject", "Failed to get %s from %s", 
	    fCorrName.Data(), fn.Data());
      return 0;
    }
    file->Close();
    return obj;
  }
  TString fDirName;
  TString fCorrName;
  TString fMethName;
  TString fFile;
  AliCorrectionManagerBase* fCM;
};

//====================================================================
struct NormExtractor : public Extractor
{
  enum { 
    kINEL, 
    kNSD, 
    kINELGT0
  };
  TString fFileName;
  /** 
   * Constructor 
   * 
   * @param dirName 
   * @param corrName 
   * @param methName 
   * 
   * @return 
   */
  NormExtractor(const char* dirName, 
		const char* corrName,
		const char* methName)
    : Extractor(dirName,corrName,methName,"",0),
      fFileName("")
  {
  }
  /** 
   * Extract 
   * 
   * 
   * @return Number of converted oject
   */
  virtual Int_t Extract()
  {
    Fatal("Extract", "Cannot use this");
    return -1;
  }
  /** 
   * Extract 
   * 
   * @param s Source 
   *
   * @return Number of converted oject
   */
  virtual Bool_t ExtractFile(const TString& s)
  {
    Fatal("ExtractFile", "Cannot use this (%s)", s.Data());
    return -1;
  }
  /** 
   * Extract files 
   * 
   * @param db        Database manager 
   * @param fileName  File to store in 
   * 
   * @return number of converted objects 
   */
  virtual Int_t ExtractNorm(AliOADBForward& db, const char* fileName)
  {
    Scanner s;
    TString dir = TString::Format("$ALICE_PHYSICS/PWGLF/FORWARD/corrections/%s",
				  fDirName.Data());
    TList* l = s.Scan(dir, fCorrName);
    if (!l) {
      Warning("ExtractNorm", "No files matching %s found in %s",  
	      fCorrName.Data(), dir.Data()); 
      return 0;
    }

    fFileName = fileName;
    // if (!Open(db, fileName)) return 0;
    TIter next(l);
    TObjString* os = 0;
    Int_t ret = 0;
    while ((os = static_cast<TObjString*>(next()))) { 
      TString& fn  = os->String();
      if (ExtractNormFile(fn, db)) ret++;
    }
    return ret;
  }
  /** 
   * Overload to store file name 
   * 
   * @return true
   */
#if 0
  virtual Bool_t Open(AliOADBForward& db, 
		      const Char_t* fileName) 
  { 
    fFileName = fileName; 
    Info("Open", "file name set to %s", fFileName.Data());
    return true;
  }
#endif
  /** 
   * Store object in DB
   * 
   * @param db     Database manager
   * @param tab    Table name 
   * @param o      Object to stire 
   * @param runNo  Run number
   * @param sys    System
   * @param sNN    Energy
   * 
   * @return true on success
   */
  virtual Bool_t Store(AliOADBForward& db, const TString& tab, TObject* o, 
		       ULong_t runNo, UShort_t sys, UShort_t sNN)
  {
    Info("Store", "file name to store in %s", fFileName.Data());
    if (!db.Open(fFileName, Form("%s/%s", tab.Data(), fMethName.Data()), 
		 true, true)) { 
      Warning("Store", "Failed to open for %s/%s", tab.Data(), 
	      fMethName.Data());
      return false;
    }
    return db.Insert(tab, o, runNo, sys, sNN, 0, false, false);
  }
  /** 
   * Extract a histogram 
   * 
   * @param what  Name part
   * @param runNo Run number
   * @param sys   System
   * @param sNN   Energy 
   * @param f     File to read from 
   * @param db    Database manager 
   * 
   * @return true on success
   */    
  virtual Bool_t ExtractHist(Int_t what, ULong_t runNo, 
			     UShort_t sys, UShort_t sNN,
			     TFile& f, AliOADBForward& db)
  {
    TString oName;
    switch (what) { 
    case kINEL:    oName = "hInelNormalization";    break;
    case kNSD:     oName = "hNSDNormalization";     break;
    case kINELGT0: oName = "hINELGT0Normalization"; break;
    }                       
    TObject* obj = f.Get(oName);
    if (!obj) {
      Warning("ExtractHist", "Object %s not found", oName.Data());
      return false;
    }
    
    TString tName;
    TString ttName;
    switch (what) { 
    case kINEL:     tName = "normalizationINEL";    ttName = "INEL";   break;
    case kNSD:      tName = "normalizationNSD";     ttName = "NSD";    break;
    case kINELGT0:  tName = "normalizationINELGT0"; ttName = "INEL>0"; break;
    }
    TH1* hist = static_cast<TH1*>(obj->Clone(tName));
    hist->SetDirectory(0);
    hist->SetTitle(Form("Normalization for %s", ttName.Data()));
    // obj->SetName(tName.Data());

    return Store(db, tName, hist, runNo, sys, sNN);
  }
  /** 
   * Extract a number 
   * 
   * @param what  Name part
   * @param runNo Run number
   * @param sys   System
   * @param sNN   Energy 
   * @param f     File to read from 
   * @param db    Database manager 
   * 
   * @return true on success
   */    
  virtual Bool_t ExtractNum(Int_t what, ULong_t runNo, 
			    UShort_t sys, UShort_t sNN,
			    TFile& f, AliOADBForward& db)
  {
    TString oName;
    switch (what) { 
    case kINEL:    oName = "inelTriggerEff";    break;
    case kNSD:     oName = "nsdTriggerEff";     break;
    case kINELGT0: oName = "inelgt0TriggerEff"; break;
    }
    TObject* obj = f.Get(oName);
    if (!obj) {
      Warning("ExtractHist", "Object %s not found", oName.Data());
      return false;
    }
    
    TString tName;
    switch (what) { 
    case kINEL:     tName = "triggerEffINEL";    break;
    case kNSD:      tName = "triggerEffNSD";     break;
    case kINELGT0:  tName = "triggerEffINELGT0"; break;
    }
    TParameter<float>* p = static_cast<TParameter<float>*>(obj->Clone(tName));
    
    return Store(db, tName, p, runNo, sys, sNN); 
  }
  /** 
   * Extract from a file 
   * 
   * @param fn File name 
   * @param db database manager 
   * 
   * @return true on success 
   */      
  virtual Bool_t ExtractNormFile(const TString& fn, AliOADBForward& db) 
  {
    UShort_t sys = 0;
    UShort_t sNN = 0;
    Short_t  fld = 0;
    Bool_t   mc  = false;
    
    ExtractFields(fn, sys, sNN, fld, mc);
    if (sNN == 2750) sNN = 2760;

    ULong_t runNo = ExtractRunNo(sys, sNN);
    if (runNo == 0xFFFFFFFF || runNo <= 0) return false;


    TFile* f = TFile::Open(fn, "READ");
    if (!f) { 
      Error("ExtractFile", "Failed to open %s", fn.Data());
      return false;
    }
    ExtractHist(kINEL,   runNo, sys, sNN, *f, db);
    ExtractHist(kNSD,    runNo, sys, sNN, *f, db);
    ExtractHist(kINELGT0,runNo, sys, sNN, *f, db);
    ExtractNum(kINEL,    runNo, sys, sNN, *f, db);
    ExtractNum(kNSD,     runNo, sys, sNN, *f, db);
    ExtractNum(kINELGT0, runNo, sys, sNN, *f, db);
    return true;
  }
};

//====================================================================
Extractor*
MakeFMDExtractor(const char* dir, const char* name)
{
  return new Extractor(dir, name, "NEAR", "fmd_corrections.root", 
		       &(AliForwardCorrectionManager::Instance()));
}
Extractor*
MakeSPDExtractor(const char* dir, const char* name)
{
  return new Extractor(dir, name, "NEAR", "spd_corrections.root", 
		       &(AliCentralCorrectionManager::Instance()));
}

void
MigrateOADB(Int_t what=0x3)
{
  gROOT->Macro("$ALICE_PHYSICS/PWGLF/FORWARD/analysis2/scripts/LoadLibs.C");

  if (what & 0x1) {
    Extractor*  ee[] = {
      MakeFMDExtractor("Acceptance",        "acceptance"),
      MakeFMDExtractor("CentralAcceptance", "centralacceptance"),
      MakeFMDExtractor("CentralSecMap",     "centralsecmap"),
      MakeFMDExtractor("DoubleHit",         "doublehit"),
      MakeFMDExtractor("ELossFits",         "elossfits"),
      MakeFMDExtractor("MergingEfficiency", "merging"),
      MakeFMDExtractor("SecondaryMap",      "secondary"),
      MakeFMDExtractor("VertexBias",        "vertexbias"),
      MakeSPDExtractor("CentralSecMap",    "centralsecmap"),
      MakeSPDExtractor("CentralAcceptance","centralacceptance"),
      0 };

    gSystem->Unlink("fmd_corrections.root");
    gSystem->Unlink("spd_corrections.root");
    
    Extractor** ep   = ee;
    while (*ep) { 
      (*ep)->Extract();
      ep++;
    }
  }

  if (what & 0x2) {
    gSystem->Unlink("normalization.root");
    
    NormExtractor e7("Normalization",
		     "normalizationHists", 
		     "NEAR");
    AliOADBForward ndb;
    TString ntables;
    e7.ExtractNorm(ndb,"normalization.root");
  }
}

