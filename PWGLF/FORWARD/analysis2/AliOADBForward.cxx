#include "AliOADBForward.h"
#include <TBrowser.h>
#include <TROOT.h>
#include <TKey.h>
#include <TList.h>
#include <TDatime.h>
#include <TTree.h>
#include <TFile.h>
#include <TError.h>
#include <TSystem.h>

#ifndef ALIROOT_SVN_REVISION
# define ALIROOT_SVN_REVISION 0
#endif


ClassImp(AliOADBForward)
#if 0
; // Do not remove - for Emacs
#endif

//====================================================================
const char* 
AliOADBForward::Mode2String(ERunSelectMode mode)
{
  switch (mode) { 
  case kDefault:   return "default"; 
  case kExact:     return "exact"; 
  case kNewest:    return "newest";
  case kNear:      return "near";
  case kOlder:     return "older";
  case kNewer:     return "newer";
  }
  return "?"; // Should never get here 
}
AliOADBForward::ERunSelectMode
AliOADBForward::String2Mode(const TString& str)
{
  if      (str.EqualTo("default", TString::kIgnoreCase)) return kDefault;
  else if (str.EqualTo("exact",   TString::kIgnoreCase)) return kExact;
  else if (str.EqualTo("newest",  TString::kIgnoreCase)) return kNewest;
  else if (str.EqualTo("near",    TString::kIgnoreCase)) return kNear;
  else if (str.EqualTo("older",   TString::kIgnoreCase)) return kOlder;
  else if (str.EqualTo("newer",   TString::kIgnoreCase)) return kNewer;
  return kDefault;
}
AliOADBForward::ERunSelectMode
AliOADBForward::Int2Mode(Int_t mode)
{
  switch (mode) { 
  case kDefault:   return kDefault; 
  case kExact:     return kExact; 
  case kNewest:    return kNewest;
  case kNear:      return kNear;
  case kOlder:     return kOlder;
  case kNewer:     return kNewer;
  }
  return kDefault; // Should never get here 
}


//====================================================================
AliOADBForward::Entry::Entry(ULong_t  runNo, 
			     UShort_t sys, 
			     UShort_t sNN, 
			     Short_t  field, 
			     Bool_t   mc, 
			     Bool_t   sat,
			     TObject* o)
  : fRunNo(runNo), 
    fSys(sys), 
    fSNN(sNN), 
    fField(field),
    fMC(mc), 
    fSatellite(sat),
    fData(o),
    fTimestamp(0), 
    fAliROOTRevision(0),
    fAuthor("unknown")
{
  // 
  // Constructor 
  // 
}

//____________________________________________________________________
AliOADBForward::Entry::Entry(const Entry& o)
  : TObject(o), 
    fRunNo(o.fRunNo), 
    fSys(o.fSys), 
    fSNN(o.fSNN), 
    fField(o.fField),
    fMC(o.fMC), 
    fSatellite(o.fSatellite), 
    fData(o.fData),
    fTimestamp(0), 
    fAliROOTRevision(0),
    fAuthor(o.fAuthor)
{
  // 
  // Copy constructor 
  //
}
//____________________________________________________________________
AliOADBForward::Entry&
AliOADBForward::Entry::operator=(const Entry& o)
{
  // 
  // Assignment operator 
  //
  if (this == &o) return *this;
  fRunNo	      = o.fRunNo; 
  fSys	 	      = o.fSys; 
  fSNN	     	      = o.fSNN; 
  fField	      = o.fField;
  fMC	              = o.fMC; 
  fSatellite          = o.fSatellite;
  fData	              = o.fData;
  fTimestamp          = o.fTimestamp;
  fAliROOTRevision    = o.fAliROOTRevision;
  fAuthor             = o.fAuthor;

  return *this;
}
//____________________________________________________________________
const char*
AliOADBForward::Entry::GetTitle() const
{
  TDatime d(fTimestamp);
  return Form("%09ld, %4s, %4huGeV, %+4hd, %4s, %3s, %19s: %p (%s) by %s", 
	      (fRunNo == 0xFFFFFFFF ? -1 : fRunNo), 
	      (fSys == 1 ? "pp"   : 
	       fSys == 2 ? "PbPb" : 
	       fSys == 3 ? "pPb"  :
	       fSys == 4 ? "Pbp"  :
	       fSys == 5 ? "XeXe" : "?"), 
	      fSNN, fField, (fMC ? "mc" : "real"),
	      (fSatellite ? "sat" : "nom"), d.AsSQLString(), fData,
	      (fData ? fData->GetName() : "?"), fAuthor.Data());
  
}
//____________________________________________________________________
void
AliOADBForward::Entry::Print(Option_t* /*option*/) const 
{
  Printf("%s", GetTitle());
}
//====================================================================
AliOADBForward::Table::Table(TTree* tree, Bool_t isNew, ERunSelectMode mode)
  : fTree(tree), fEntry(0), fVerbose(false), fMode(mode), fFallBack(false)
{
  if (!tree) return;

#if 0
  Info("Table", "Making table %s (%s) with mode %s (%s)", 
       tree->GetName(), (isNew ? "new" : "old"), tree->GetTitle(), 
       Mode2String(fMode));
#endif
  if (isNew) {
    fTree->Branch("e", "AliOADBForward::Entry", &fEntry);
    fMode = String2Mode(fTree->GetTitle());
  }
  else {
    if (fMode <= kDefault || fMode > kNewer) {
      fMode = String2Mode(fTree->GetTitle());
      if (fMode == kDefault) fMode = kNear;
    }
    fTree->SetBranchAddress("e", &fEntry);
  }
#if 0
  Info("", "Mode set to %d (%s)", fMode, Mode2String(fMode));
#endif
}
//____________________________________________________________________
AliOADBForward::Table::Table(const Table& o)
  : TObject(o), 
    fTree(o.fTree), 
    fEntry(o.fEntry), 
    fVerbose(o.fVerbose),
    fMode(o.fMode), 
    fFallBack(o.fFallBack)
{
  //
  // Copy constructor 
  if (!fTree) return;
  fTree->SetBranchAddress("e", &fEntry);
}

//____________________________________________________________________
AliOADBForward::Table::~Table()
{
  // 
  // Close this table 
  //
  Close();
}
//____________________________________________________________________
AliOADBForward::Table&
AliOADBForward::Table::operator=(const Table& o)
{
  // 
  // Assignment operator 
  // 
  if (this == &o) return *this;
  fTree    = o.fTree;
  fEntry   = o.fEntry;
  fVerbose = o.fVerbose;
  fMode    = o.fMode;
  if (fTree) fTree->SetBranchAddress("e", &fEntry);

  return *this;
}

//____________________________________________________________________
const Char_t*
AliOADBForward::Table::GetTableName() const
{
  // 
  // Get the table name or null
  if (!fTree) return 0;
  return fTree->GetName();
}
//____________________________________________________________________
const Char_t*
AliOADBForward::Table::GetName() const
{
  // 
  // Overload of TObject::GetName
  //
  if (!fTree) return TObject::GetName();
  return GetTableName();
}
//____________________________________________________________________
Bool_t
AliOADBForward::Table::Update()
{
  // 
  // Flush to disk 
  //
  if (!IsOpen()) { 
    Error("Update", "No tree associated");
    return false;
  }
      
  TFile* file = fTree->GetCurrentFile();
  if (!file) { 
    Error("Update", "No file associated with tree");
    return false;
  }
  if (!file->IsWritable()) { 
    Error("Update", "File %s not writeable", file->GetName());
    return false;
  }
      
  Int_t nBytes = file->Write();
      
  return (nBytes >= 0);
}
//____________________________________________________________________
Bool_t
AliOADBForward::Table::Close()
{
  // 
  // Close the connection 
  //
  if (!IsOpen()) { 
    Error("Close", "No tree associated");
    return false;
  }
  
  // if (fTree)  delete fTree; 
  // if (fEntry) delete fEntry;
  fTree  = 0;
  fEntry = 0;
  return true;
} 
//____________________________________________________________________
Int_t
AliOADBForward::Table::Query(ULong_t        runNo,
			     ERunSelectMode mode,
			     UShort_t       sys,
			     UShort_t       sNN, 
			     Short_t        fld,
			     Bool_t         mc,
			     Bool_t         sat) const
{
  // 
  // Query the tree 
  //
  return Query(runNo, mode, Conditions(sys, sNN, fld, mc, sat));
}

//____________________________________________________________________
Int_t
AliOADBForward::Table::Query(ULong_t        runNo,
			     ERunSelectMode mode,
			     const TString& q) const
{
  // 
  // Run a query against the table 
  // 
  
  // Check the tree 
  if (!IsOpen()) { 
    Error("Close", "No tree associated");
    return -1;
  }
      
  TString query = q;
  const char* smode = "latest";
  if (runNo > 0) {
    if (mode <= kDefault || mode > kNewer) mode = fMode;
    smode = Mode2String(mode);
    switch (mode) { 
    case kExact:  
      AppendToQuery(query, Form("fRunNo == %lu", runNo)); 
      break;
    case kNewest: 
      break;
    case kNear:   
      AppendToQuery(query, Form("abs(fRunNo-%lu)<=%d",
				runNo,kMaxNearDistance)); 
      break;
    case kOlder: 
      AppendToQuery(query, Form("fRunNo <= %lu", runNo));
      break;
    case kNewer: 
      AppendToQuery(query, Form("fRunNo >= %lu", runNo));
      break;
    case kDefault: 
      Fatal("Query", "Mode should never be 'default'");
      break;
    }
  }
      
  if (query.IsNull()) {
    Warning("Query", "Empty query!");
    return -1;
  }

  if (fVerbose) 
    Printf("%s: Query is '%s'", GetName(), query.Data());
  fTree->Draw("Entry$:fRunNo:fTimestamp", query, "goff");
  Int_t nRows = fTree->GetSelectedRows();
  if (nRows <= 0) return -1;
      
  if (fVerbose) 
    Printf("Query: %s (%s)\n"
	   " Entry  |    Run    | Timestamp \n"
	   "--------+-----------+------------------------", 
	   query.Data(), smode);
      
  ULong_t  oldRun  = (mode == kNewer ? 0xFFFFFFFF : 0);
  ULong_t  oldTim  = 0;
  ULong_t  oldDist = 0xFFFFFFFF;
  Int_t    entry  = -1;  
      
  for (Int_t row = 0; row < nRows; row++) {
    Int_t    ent  = Int_t(fTree->GetV1()[row]);
    ULong_t  run  = ULong_t(fTree->GetV2()[row]);
    ULong_t  tim  = ULong_t(fTree->GetV3()[row]);
    ULong_t  dist = (run > runNo ? run - runNo : runNo - run);
	
    if (fVerbose) {
      TDatime t(tim);
      Printf(" %6d | %9ld | %19s   [dist: %9ld vs %9ld, %5s]",
	     ent, run > 0x7FFFFFFF ? -1 : run, t.AsSQLString(),
	     dist, oldDist, (tim<oldTim?"older":"newer"));
    }

    switch (mode) { 
    case kExact: break; // Done in the draw `query' 
    case kNewest: // Fall-through 
    case kOlder: 
      if (run < oldRun) continue;
      break;
    case kNewer: 
      if (run > oldRun) continue;
      break;
    case kNear: 
      if (runNo > 0 && dist > oldDist) continue;
      break;
    case kDefault:
      break;
    }
    // If we get here, then we have the best run according to the 
    // criteria 
	    
    // Now update last values and set current best entry 
    oldTim  = tim;
    oldDist = dist;
    oldRun  = run;
    entry   = ent;

    // Finally, check the timestamp 
    if (tim < oldTim) continue;
	
    // Now update last values and set current best entry 
    oldTim  = tim;
    oldDist = dist;
    oldRun  = run;
    entry   = ent;
  }

  if (fVerbose) {
    Printf("Returning entry # %d", entry);
  }
  return entry;
}

//____________________________________________________________________
Bool_t
AliOADBForward::Table::Insert(TObject* o, 
			      ULong_t  runNo, 
			      UShort_t sys, 
			      UShort_t sNN, 
			      Short_t  field, 
			      Bool_t   mc, 
			      Bool_t   sat,
			      ULong_t  aliRev,
			      const TString& author) 
{
  // 
  // Insert a new row in the table 
  //

  // Check if the file is read/write 
  if (fVerbose) 
    Info("Insert", "Inserting object %p for run=%lu, sys=%hu, sNN=%4hu, "
	 "field=%+2hd, mc=%d, sat=%d", o,runNo, sys, sNN, field, mc, sat);

  if (!IsOpen(true)) {
    Warning("Insert", "No tree, or not write-able");
    return false;
  }
      
  // If the entry doesn't exists 
  if (!fEntry) fEntry = new Entry;

  // Make author 
  TString auth(author);
  if (auth.IsNull()) { 
    UserGroup_t* u = gSystem->GetUserInfo();
    TInetAddress i = gSystem->GetHostByName(gSystem->HostName());
    auth = TString::Format("%s <%s@%s>", u->fRealName.Data(), 
			   u->fUser.Data(), i.GetHostName());
  }
    
  // Set fields 
  fEntry->fData            = o;
  fEntry->fRunNo           = runNo; // (runNo <= 0 ? 0xFFFFFFFF : runNo);
  fEntry->fSys             = sys;
  fEntry->fSNN             = sNN;
  fEntry->fField           = field;
  fEntry->fMC              = mc;
  fEntry->fSatellite       = sat;
  fEntry->fAliROOTRevision = (aliRev != 0 ? aliRev : ALIROOT_SVN_REVISION);
  fEntry->fAuthor          = auth;

  TDatime now;
  fEntry->fTimestamp       = now.Convert(true);

  // Fill into tree 
  Int_t nBytes = fTree->Fill();
  if (nBytes <= 0) {
    Warning("Insert", "Failed to insert new entry");
    return false;
  }
    
  // do an Auto-save and flush-baskets now 
  fTree->AutoSave("FlushBaskets SaveSelf");

  return true;
}

//____________________________________________________________________
Int_t
AliOADBForward::Table::GetEntry(ULong_t        run,
				ERunSelectMode mode,
				UShort_t       sys,
				UShort_t       sNN, 
				Short_t        fld,
				Bool_t         mc,
				Bool_t         sat) const
{
  // Query the tree for an object.  The strategy is as follows. 
  // 
  //  - First query with all fields 
  //    - If this returns a single entry, return that. 
  //    - If not, then ignore the run number (if given)
  //      - If this returns a single entry, return that 
  //      - If not, and fall-back is enabled, then 
  //        - Ignore the collision energy (if given) 
  //          - If this returns a single entry, return that.  
  //          - If not, ignore all passed values 
  //            - If this returns a single entry, return that.
  //            - Otherwise, give up and return -1
  //      - Otherwise, give up and return -1
  //
  // This allow us to specify default objects for a period, and for
  // collision system, energy, and field setting.
  //

  if (fVerbose)
    Printf("run=%lu mode=%s sys=%hu sNN=%hu fld=%hd mc=%d sat=%d (fall=%d)",
	   run, Mode2String(mode), sys, sNN, fld, mc, sat, fFallBack);
  Int_t entry = Query(run, mode, sys, sNN, fld, mc, sat);
  if (entry < 0 && run > 0) 
    entry = Query(0, mode, sys, sNN, fld, mc, sat);
  if (entry < 0 && fFallBack && fld != kInvalidField) {
    if (fVerbose) 
      Printf("Fall-back enabled, trying without field=%d", fld);
    entry = Query(run, mode, sys, sNN, fld, mc, sat);
  }
  if (entry < 0 && fFallBack && sNN > 0) {
    if (fVerbose)
      Printf("Fall-back enabled, will try without sNN=%d", sNN);
    entry = Query(run, mode, sys, 0, fld, mc, sat);
  }
  if (entry < 0 && fFallBack) {
    if (fVerbose)
      Printf("Fall-back enabled, will try without any fields");
    entry = Query(0, mode, 0, 0, kInvalidField, false, false);	
  }
  if (entry < 0) {
    Warning("Get", "No valid object could be found");
    return -1;
  }
  return entry;
}

//____________________________________________________________________
AliOADBForward::Entry*
AliOADBForward::Table::Get(ULong_t        run,
			   ERunSelectMode mode,
			   UShort_t       sys,
			   UShort_t       sNN, 
			   Short_t        fld,
			   Bool_t         mc,
			   Bool_t         sat) const
{
  // Query the tree for an object.  The strategy is as follows. 
  // 
  //  - First query with all fields 
  //    - If this returns a single entry, return that. 
  //    - If not, then ignore the run number (if given)
  //      - If this returns a single entry, return that 
  //      - If not, and fall-back is enabled, then 
  //        - Ignore the collision energy (if given) 
  //          - If this returns a single entry, return that.  
  //          - If not, ignore all passed values 
  //            - If this returns a single entry, return that.
  //            - Otherwise, give up and return null
  //      - Otherwise, give up and return null
  //
  // This allow us to specify default objects for a period, and for
  // collision system, energy, and field setting.
  //
  Int_t entry  = GetEntry(run, mode, sys, sNN, fld, mc, sat);
  if (entry < 0) return 0;

  Int_t nBytes = fTree->GetEntry(entry);
  if (nBytes <= 0) { 
    Warning("Get", "Failed to get entry # %d\n", entry);
    return 0;
  }
  if (fVerbose) fEntry->Print();
  return fEntry;
}
//____________________________________________________________________
TObject*
AliOADBForward::Table::GetData(ULong_t        run,
			       ERunSelectMode mode,
			       UShort_t       sys,
			       UShort_t       sNN, 
			       Short_t        fld,
			       Bool_t         mc,
			       Bool_t         sat) const
{
  // 
  // Get data associated with a row or null. 
  // See also AliOADBForward::Get
  // 
  Entry* e = Get(run, mode, sys, sNN, fld, mc, sat);
  if (!e) return 0;
  return e->fData;
}
//____________________________________________________________________
void
AliOADBForward::Table::Print(Option_t* option) const
{
  // 
  // Print the full table 
  //
  if (!IsOpen()) return;

  Printf("Table %s (default mode: %s)", GetName(), Mode2String(fMode));
  Int_t n = fTree->GetEntries();
  for (Int_t i = 0; i < n; i++) { 
    fTree->GetEntry(i);
    printf("%4d/%4d: ", i, n);
    fEntry->Print(option);
  }
}
//____________________________________________________________________
void
AliOADBForward::Table::Browse(TBrowser* b) 
{
  // Browse this table 
  if (fTree) b->Add(fTree);
}
//____________________________________________________________________
Bool_t
AliOADBForward::Table::IsOpen(Bool_t rw) const 
{ 
  if (!fTree) return false;
  if (!rw)    return true;
  
  return fTree->GetCurrentFile()->IsWritable();
}
//====================================================================
AliOADBForward::AliOADBForward() 
  : TObject(),
    fTables()
{
  //
  // Constructor 
  //
}
#if 0
//____________________________________________________________________
AliOADBForward::AliOADBForward(const AliOADBForward& other)
  : TObject(other), 
    fTables(other.fTables)
{
  // 
  // Copy constructor 
  // 
}
#endif
//____________________________________________________________________
AliOADBForward::~AliOADBForward()
{
  // 
  // Destructor 
  // 
  Close();
}
#if 0
//____________________________________________________________________
AliOADBForward&
AliOADBForward::operator=(const AliOADBForward& other)
{
  // 
  // Copy constructor 
  // 
  if (&other == this) return *this;

  fTables = other.fTables;

  return *this;
}
#endif 

//____________________________________________________________________
Bool_t
AliOADBForward::Open(const  TString& fileName, 
		     const  TString& tables, 
		     Bool_t          rw, 
		     Bool_t          verb,
		     Bool_t          fallback)
{
  TString  absPath(gSystem->ExpandPathName(fileName));
  if (absPath.IsNull()) { 
    Error("Open", "Empty path for tables %s", tables.Data());
    return false;
  }
  TObject* previous = gROOT->GetListOfFiles()->FindObject(absPath);
  TFile*   file     = 0;
  if (previous) {
    file = static_cast<TFile*>(previous);
  }
  else { 
    file = TFile::Open(fileName, (rw ? "UPDATE" : "READ"));
  }
  if (!file)  { 
    Error("Open", "Failed to open %s", GetName());
    return false;
  }
  return Open(file, tables, rw, verb, fallback);
}

//____________________________________________________________________
Bool_t
AliOADBForward::Open(TFile*         file,
		     const TString& tables,
		     Bool_t         rw, 
		     Bool_t         verb,
		     Bool_t         fallback) 
{
  // 
  // Open database file and find or create listed tables 
  // 
  if (!file) return false;
  if (rw && !file->IsWritable()) {
    Warning("Open", "Read+write access requested, but %s opened read-only",
	    file->GetName());
    if (file->ReOpen("UPDATE") < 0) { 
      Error("Open", "Failed to reopen file in read+write access mode");
      return false;
    }
  }

  if (tables.EqualTo("*")) {
    if (rw) { 
      Error("Open", "Cannot open with '*' in read/write mode");
      // return false;
    }
    TList* l = file->GetListOfKeys();
    TIter  next(l);
    TKey*  key = 0;
    while ((key = static_cast<TKey*>(next()))) { 
      TClass* cl = gROOT->GetClass(key->GetClassName());
      if (!cl) continue; 
      if (!cl->InheritsFrom(TTree::Class())) continue; 
	
      OpenTable(file, rw, key->GetName(), "DEFAULT", verb, fallback);
    }
    // file->Close();
    return true;
  }
  TObjArray*  tokens = tables.Tokenize(":,");
  TObjString* token  = 0;
  TIter       nextToken(tokens);
  while ((token = static_cast<TObjString*>(nextToken()))) {
    TString& tn = token->String();
    if (tn.IsNull()) continue;
     
    TObjArray*  parts = tn.Tokenize("/");
    TObjString* onam  = static_cast<TObjString*>(parts->At(0));
    TString&    name  = onam->String();
    TString     mode  = "DEFAULT";
    if (parts->GetEntries() > 1) 
      mode = static_cast<TObjString*>(parts->At(1))->String();
    mode.ToUpper();

    OpenTable(file, rw, name, mode, verb, fallback);

    delete parts;
  }
  delete tokens;

  return true;
}

//____________________________________________________________________
Bool_t
AliOADBForward::Close()
{
  // 
  // Flush all tables and close all files 
  // 
  TList  files;
  Int_t nFiles = GetFiles(files);
  if (nFiles <= 0) { 
    // Nothing to close 
    return true;
  }

  TIter nextFile(&files);
  TFile* file = 0;
  while ((file = static_cast<TFile*>(nextFile()))) {
    // if (file->IsWritable()) file->Write();

    file->Close();
  }
    
  fTables.DeleteAll();

  return true;
}
//____________________________________________________________________
Bool_t
AliOADBForward::Update()
{
  // 
  // Flush all tables 
  // 
  TList  files;
  Int_t nFiles = GetFiles(files);
  if (nFiles <= 0) { 
    // Nothing to close 
    return true;
  }

  TIter nextFile(&files);
  TFile* file = 0;
  Int_t  nBytes = 0;
  while ((file = static_cast<TFile*>(nextFile()))) {
    if (!file->IsWritable()) { 
      Error("Update", "File %s not writeable", file->GetName());
      continue;
    }

    nBytes += file->Write();
  }
  return (nBytes >= 0);
}
//____________________________________________________________________
AliOADBForward::Entry*
AliOADBForward::Get(const TString& table, 
		    ULong_t        run,
		    ERunSelectMode mode, 
		    UShort_t       sys,
		    UShort_t       sNN, 
		    Short_t        fld,
		    Bool_t         mc,
		    Bool_t         sat) const
{
  // 
  // Get a row from selected table 
  // 
  Table* t = FindTable(table);
  if (!t) return 0;
    
  return t->Get(run, mode, sys, sNN, fld, mc, sat);
}
//____________________________________________________________________
TObject*
AliOADBForward::GetData(const TString& table, 
			ULong_t        run,
			ERunSelectMode mode, 
			UShort_t       sys,
			UShort_t       sNN, 
			Short_t        fld,
			Bool_t         mc,
			Bool_t         sat) const
{
  Table* t = FindTable(table);
  if (!t) return 0;

  return t->GetData(run, mode, sys, sNN, fld, mc, sat);
}
//____________________________________________________________________
Bool_t
AliOADBForward::Insert(const TString& table, 
		       TObject* o, 
		       ULong_t  runNo, 
		       UShort_t sys, 
		       UShort_t sNN, 
		       Short_t  field, 
		       Bool_t   mc, 
		       Bool_t   sat,
		       ULong_t  aliRev,
		       const TString& author) 
{
  // 
  // Insert a row into the selected table 
  //
  Table* t = FindTable(table);
  if (!t) return false;

  return t->Insert(o, runNo, sys, sNN, field, mc, sat, aliRev, author);
}	    
//____________________________________________________________________
Bool_t
AliOADBForward::CopyEntry(const TString& table, 
			  ULong_t        oldRunNo, 
			  UShort_t       oldSys, 
			  UShort_t       oldSNN, 
			  Short_t        oldField, 
			  ULong_t        newRunNo, 
			  UShort_t       newSys, 
			  UShort_t       newSNN, 
			  Short_t        newField, 
			  Bool_t         mc, 
			  Bool_t         sat)
{
  Table* t = FindTable(table);
  if (!t) return false;

  Printf("Copy entry: table mode: %s", Mode2String(t->fMode));
  Entry* e = t->Get(oldRunNo, t->fMode, oldSys, oldSNN, oldField, mc, sat);
  if (!e) return false;

  return t->Insert(e->fData, newRunNo, newSys, newSNN, newField, mc, sat, 
		e->fAliROOTRevision, e->fAuthor);
}

//____________________________________________________________________
void
AliOADBForward::Print(const Option_t* option) const
{
  // 
  // Print everything 
  //
  TIter       nextTable(&fTables);
  TObjString* key   = 0;
  Table*      table  = 0;
  while ((key = static_cast<TObjString*>(nextTable()))) {
    Printf("Table: %p", key->GetName());
    table = static_cast<Table*>(fTables.GetValue(key));
    if (!table) continue;
    table->Print(option);
  }
}

//____________________________________________________________________
void
AliOADBForward::Browse(TBrowser* b) 
{
  // 
  // Browse this object
  // 
  TIter       nextTable(&fTables);
  TObjString* key   = 0;
  Table*      table  = 0;
  while ((key = static_cast<TObjString*>(nextTable()))) {
    table = static_cast<Table*>(fTables.GetValue(key));
    if (!table) continue;
    b->Add(table, key->GetName());
  }
}
//____________________________________________________________________
AliOADBForward::Table*
AliOADBForward::FindTable(const TString& name, Bool_t quite) const
{
  //
  // Find a table by name 
  // 
  TPair* p = static_cast<TPair*>(fTables.FindObject(name));
  if (!p) {
    if (!quite)
      Warning("FindTable", "Table %s not registered", name.Data());
    return 0; 
  }
  return static_cast<Table*>(p->Value());
}
//____________________________________________________________________
Int_t
AliOADBForward::GetFiles(TList& files) const
{
  // 
  // Get all associated files 
  // 
  Int_t  ret    = 0;
  TIter       nextTable(&fTables);
  TObjString* key   = 0;
  Table*      table  = 0;
  while ((key = static_cast<TObjString*>(nextTable()))) {
    table = static_cast<Table*>(fTables.GetValue(key));
    if (!table->fTree) continue; 

    TFile* f = table->fTree->GetCurrentFile();
    if (!f) continue;

    if (files.FindObject(f)) continue;
    files.Add(f);
    ret++;
  }
  return ret;
}
//____________________________________________________________________
AliOADBForward::Table*
AliOADBForward::GetTableFromFile(TFile*         file, 
				 Bool_t         rw, 
				 const TString& name,
				 const TString& mode) const
{
  // 
  // Get a table from a file, or make a new table 
  // 
  if (!file) return 0;
  if (FindTable(name, true)) return 0;

  TObject* o = file->Get(name);
  TTree*   t = 0;
  Bool_t   n = false;
  if (!o) { 
    if (!rw) { 
      // We only fail if in read-only mode 
      Error("Open", "No such object: %s in %s", name.Data(),
	    file->GetName());
      return 0;
    }
    TDirectory* savDir = gDirectory;
    file->cd();
    // Create the tree in the file 
    t = new TTree(name, mode);
    t->SetDirectory(file);
    if (savDir) savDir->cd();
    n = true;
  }
  else {
    // Get the tree, and set the branch
    t = static_cast<TTree*>(o);
  }
  Table* ret = new Table(t, n, String2Mode(mode));
  return ret;
}

//____________________________________________________________________
void
AliOADBForward::OpenTable(TFile*         file, 
			  Bool_t         rw, 
			  const TString& name,
			  const TString& mode,
			  Bool_t         verb, 
			  Bool_t         fallback)
{
  if (!file) return;

  Table* t = GetTableFromFile(file, rw, name, mode);
  if (!t) return;

  fTables.Add(new TObjString(name), t);
  t->SetVerbose(verb);
  t->SetEnableFallBack(fallback);
  if (verb) 
    Printf("Found table %s. Opened with verbose=%d, fallback=%d", 
	   name.Data(), verb, fallback);
}

//____________________________________________________________________
void
AliOADBForward::AppendToQuery(TString& q, const TString& s, Bool_t andNotOr)
{
  // 
  // Helper function 
  // 
  if (!q.IsNull()) q.Append(andNotOr ? " && " : " || ");
  q.Append(s);
}
//____________________________________________________________________
TString
AliOADBForward::Conditions(UShort_t       sys,
			   UShort_t       sNN, 
			   Short_t        fld,
			   Bool_t         mc,
			   Bool_t         sat)
{
  // Build query string 
  TString q;
  if (sys   > 0)    AppendToQuery(q, Form("fSys == %hu",          sys));
  if (sNN   > 0)    AppendToQuery(q, Form("abs(fSNN - %hu) < 11", sNN));
  if (TMath::Abs(fld) < 10) AppendToQuery(q, Form("fField == %hd",fld));
  // Boolean fields always queried.  In cases where these do not matter, 
  // we always write down the false value, so we get the correct query 
  // anyways. 
  AppendToQuery(q, Form("%sfMC",        (mc  ? " " : "!")));
  AppendToQuery(q, Form("%sfSatellite", (sat ? " " : "!")));

  return q;
}

//____________________________________________________________________
void
AliOADBForward::TestGet(AliOADBForward& t, 
			const TString& table,
			ULong_t        runNo,
			ERunSelectMode mode,
			UShort_t       sys,
			UShort_t       sNN, 
			Short_t        fld,
			Bool_t         mc,
			Bool_t         sat)
{

  Printf("=== Test query: t=%s r=%ld s=%d t=%d f=%d m=%d v=%d",
	 table.Data(), runNo, sys, sNN, fld, int(mc), int(sat));
  AliOADBForward::Entry* e = t.Get(table, runNo, mode, sys, sNN, 
				   fld, mc, sat);
  if (!e) return;
  e->Print();
}
//____________________________________________________________________
void
AliOADBForward::TestInsert(AliOADBForward& t, 
			   const TString&  table,
			   ULong_t         runNo,
			   UShort_t        sys,
			   UShort_t        sNN, 
			   Short_t         fld,
			   Bool_t          mc,
			   Bool_t          sat)
{
  static Int_t cnt = 0;
  TString what = TString::Format("%s-%03d", table.Data(), cnt++);
  Printf("=== Insert: t=%s r=%ld s=%d t=%d f=%d m=%d v=%d w=%s",
	 table.Data(), runNo, sys, sNN, fld, int(mc), int(sat), what.Data());
  t.Insert(table, new TObjString(what), runNo, sys, sNN, fld, mc, sat);
  gSystem->Sleep(500);
}

//____________________________________________________________________
void
AliOADBForward::Test()
{
  AliOADBForward* tt = new AliOADBForward();
  if (!tt->Open("db.root", "A,B", true, true)) { 
    ::Error("Test", "Failed to open DB");
    return;
  }
  AliOADBForward& t  = *tt;
  TestInsert(t, "A", 137161);
  TestInsert(t, "A", 137161);
  TestInsert(t, "A", 0     );
  TestInsert(t, "A", 999999);
  TestInsert(t, "A", 137166);


  TestInsert(t, "B", 137161);
  TestInsert(t, "B", 0     );
  t.Print();
  t.Close();

  if (!t.Open("db.root", "A,B",false,true)) {
    ::Error("Test", "Failed to open DB");
    return;
  }

  TestGet(t, "A", 137161);
  TestGet(t, "A", 137160);
  TestGet(t, "A", 0     );
  TestGet(t, "A", 137160, kNewest);
  TestGet(t, "A", 137160, kNewer);
  TestGet(t, "A", 137168, kOlder);
  TestGet(t, "A", 137161, kExact);

  new TBrowser("b", tt);
}

//
// EOF
//
