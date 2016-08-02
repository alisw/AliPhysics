#include "AliCorrectionManagerBase.h"
#include "AliOADBForward.h"
#include "AliForwardUtil.h"
#include <AliLog.h>
#include <TMath.h>
#include <iostream>
#include <TROOT.h>
#include <TSystem.h>
#include <TBrowser.h>
#include <TParameter.h>
#include <TFileMerger.h>
#include <TBits.h>
#include <TTree.h>

//____________________________________________________________________
AliCorrectionManagerBase::AliCorrectionManagerBase()
  : fCorrections(),
    fIsInit(false),
    fRun(kIgnoreValue), 
    fSys(kIgnoreValue), 
    fSNN(kIgnoreValue), 
    fField(kIgnoreField), 
    fMC(false), 
    fSatellite(false), 
    fDB(0),
    fDebug(false),
    fFallBack(false)
{
}

//____________________________________________________________________
AliCorrectionManagerBase::AliCorrectionManagerBase(Bool_t)
  : fCorrections(16),
    fIsInit(false),
    fRun(kIgnoreValue), 
    fSys(kIgnoreValue), 
    fSNN(kIgnoreValue), 
    fField(kIgnoreField), 
    fMC(false), 
    fSatellite(false), 
    fDB(0),
    fDebug(false),
    fFallBack(false)
{
  fCorrections.SetOwner(false);
  fCorrections.SetName("corrections");
}
//____________________________________________________________________
AliCorrectionManagerBase::AliCorrectionManagerBase(const 
						   AliCorrectionManagerBase& o)
  : TObject(o),
    fCorrections(),
    fIsInit(o.fIsInit),
    fRun(o.fRun), 
    fSys(o.fSys), 
    fSNN(o.fSNN), 
    fField(o.fField), 
    fMC(o.fMC), 
    fSatellite(o.fSatellite), 
    fDB(o.fDB),
    fDebug(o.fDebug),
    fFallBack(o.fFallBack)
{
  fCorrections.SetOwner(false);
  Int_t n = o.fCorrections.GetEntriesFast();
  for (Int_t i = 0; i < n; i++) { 
    fCorrections.AddAt(o.fCorrections.At(i), i);
  }
}
//____________________________________________________________________
AliCorrectionManagerBase&
AliCorrectionManagerBase::operator=(const AliCorrectionManagerBase& o)
{
  if (&o == this) return *this;

  fIsInit	= o.fIsInit;
  fRun		= o.fRun; 
  fSys		= o.fSys; 
  fSNN		= o.fSNN; 
  fField	= o.fField; 
  fMC		= o.fMC; 
  fSatellite 	= o.fSatellite;
  fDB		= o.fDB;
  fDebug        = o.fDebug;
  fFallBack     = o.fFallBack;

  fCorrections.Clear();
  Int_t n = o.fCorrections.GetEntriesFast();
  for (Int_t i = 0; i < n; i++) { 
    fCorrections.AddAt(o.fCorrections.At(i), i);
  }
  return *this;
}

//____________________________________________________________________
AliCorrectionManagerBase::~AliCorrectionManagerBase()
{
  // fCorrections.Delete();
}

#define PF(N,V,...)					\
  AliForwardUtil::PrintField(N,V, ## __VA_ARGS__)
#define PFB(N,FLAG)				\
  do {									\
    AliForwardUtil::PrintName(N);					\
    std::cout << std::boolalpha << (FLAG) << std::noboolalpha << std::endl; \
  } while(false)
#define PFV(N,VALUE)					\
  do {							\
    AliForwardUtil::PrintName(N);			\
    std::cout << (VALUE) << std::endl; } while(false)

//____________________________________________________________________
void
AliCorrectionManagerBase::Print(Option_t* option) const
{
  AliForwardUtil::PrintTask(*this);
  gROOT->IncreaseDirLevel();
  PFB("Initialized", fIsInit);
  if (fIsInit) {
    PFV("Run number", fRun);
    PFV("Collision system", AliForwardUtil::CollisionSystemString(fSys));
    PFV("Sqrt(s_NN)",AliForwardUtil::CenterOfMassEnergyString(fSNN));
    PFV("Magnetic field", AliForwardUtil::MagneticFieldString(fField));
    PFB("For simulations", fMC);
    PFB("For satellites", fSatellite);
  }
  TString opt(option);
  opt.ToUpper();
  if (opt.Contains("R")) {
    // gROOT->IncreaseDirLevel();
    Int_t n = fCorrections.GetEntriesFast();
    for (Int_t id = 0; id < n; id++) { 
      const Correction* c = GetCorrection(id);
      c->Print(option);
    }
    // gROOT->DecreaseDirLevel();  
  }
  gROOT->DecreaseDirLevel();  
}

//____________________________________________________________________
void
AliCorrectionManagerBase::Browse(TBrowser* b)
{
  b->Add(&fCorrections);
}

//____________________________________________________________________
void
AliCorrectionManagerBase::SetPrefix(const TString& prefix)
{
  Int_t n = fCorrections.GetEntriesFast();
  for (Int_t id = 0; id < n; id++) { 
    Correction* c = GetCorrection(id);
    const char* old = c->GetTitle();
    TString     oldf(gSystem->BaseName(old));
    c->SetFile(gSystem->ConcatFileName(prefix, oldf));
  }
}

//____________________________________________________________________
Bool_t
AliCorrectionManagerBase::Store(TObject*     o,
				ULong_t     runNo,
				UShort_t    sys, 
				UShort_t    sNN, 
				Short_t     field, 
				Bool_t      mc,
				Bool_t      sat, 
				const char* file,
				const char* meth) const
{
  Bool_t ret = false;
  Int_t n = fCorrections.GetEntriesFast();
  for (Int_t id = 0; id < n; id++) { 
    const Correction* c = GetCorrection(id);
    if (!o->IsA()->InheritsFrom(c->fCls)) continue;

    ret = c->StoreIt(fDB, o, runNo, sys, sNN, field, mc, sat, file, meth);
    break;
  }
  return ret;
}
    
//____________________________________________________________________
Bool_t
AliCorrectionManagerBase::Append(const TString& addition, 
				 const TString& destination) const
{
  if (addition.IsNull()) {
    AliWarning("No addition specified");
    return false;
  }
  if (destination.IsNull()) { 
    AliWarning("No destination storage specified");
    return false;
  }
  TFileMerger merger;
  merger.SetPrintLevel(1);
  merger.OutputFile(destination, "UPDATE");
  merger.AddFile(addition);
  if (!merger.PartialMerge()) {
    AliInfoF("Failed to merge %s with %s", 
	     addition.Data(), destination.Data());
    return false;
  }
  if (destination.BeginsWith("$OADB_PATH") ||
      destination.BeginsWith("$ALICE_PHYSICS"))
    AliInfoF("Now commit %s to git", destination.Data());
  return true;
}
//____________________________________________________________________
Bool_t
AliCorrectionManagerBase::CleanUp(const TString& destination,
				  Bool_t         verb,
				  Bool_t         all) const
{
  AliOADBForward* db = fDB;
  if (!db) db = new AliOADBForward;

  for (Int_t i = 0; i < 32; i++) {
    const Correction* c = GetCorrection(i);
    if (!c) continue;

    if (!c->fEnabled) {
      Info("CleanUp", "Table %s not enabled", c->GetName());
      continue;
    }
    Info("CleanUp", "Will clean up table %s", c->GetName());
    c->CleanIt(db, destination, verb, all);
  }
  return true;
}
 
  
//____________________________________________________________________
void
AliCorrectionManagerBase::EnableCorrections(UInt_t what)
{
  for (Int_t i = 0; i < 32; i++) {
    if (!GetCorrection(i)) continue;
    EnableCorrection(i, (1 << i) & what);
  }
}
//____________________________________________________________________
Bool_t
AliCorrectionManagerBase::CheckCorrections(UInt_t what, Bool_t verbose) const
{
  Bool_t ret = true;
  for (Int_t i = 0; i < 32; i++) {
    Bool_t enabled = (1 << i) & what;
    if (!enabled) continue;
    const Correction* c = GetCorrection(i);
    if (!c) {
      ret = false;
      if (verbose) 
	AliWarningF("No correction registered for bit %d", i);
      continue;
    }
    const TObject* o = c->Get();
    if (!o) { 
      ret = false;
      if (verbose) 
	AliWarningF("Correction %s (bit %d) not initialized!", c->GetName(), i);
      continue;
    }
  }
  return ret;
}

//____________________________________________________________________
const char*
AliCorrectionManagerBase::GetObjectName(Int_t what) const
{
  const Correction* c = GetCorrection(what);
  if (!c) return 0;
  const TClass* cls = c->TheClass();
  if (!cls) return 0;
  return cls->GetName();
}

//____________________________________________________________________
void
AliCorrectionManagerBase::RegisterCorrection(Int_t id, Correction* corr)
{
  fCorrections.AddAtAndExpand(corr, id);
}

//____________________________________________________________________
void
AliCorrectionManagerBase::RegisterCorrection(Int_t id, 
					     const TString& tableName, 
					     const TString& fileName, 
					     TClass*        cls, 
					     UShort_t       fields,
					     Bool_t         enabled)
{
  RegisterCorrection(id,new Correction(tableName,fileName,cls,fields,enabled));
}

//____________________________________________________________________
AliCorrectionManagerBase::Correction*
AliCorrectionManagerBase::GetCorrection(Int_t id)
{
  if (id < 0 || id > fCorrections.GetEntriesFast()) return 0;
  return static_cast<Correction*>(fCorrections.At(id));
}

//____________________________________________________________________
const AliCorrectionManagerBase::Correction*
AliCorrectionManagerBase::GetCorrection(Int_t id) const
{
  if (id < 0 || id > fCorrections.GetEntriesFast()) return 0;
  return static_cast<Correction*>(fCorrections.At(id));
}

//____________________________________________________________________
void 
AliCorrectionManagerBase::SetCorrectionFile(Int_t id, const TString& fileName)
{
  Correction* c = GetCorrection(id);
  if (!c) return;
  c->SetFile(fileName);
}

//____________________________________________________________________
Int_t
AliCorrectionManagerBase::GetId(const TString& what) const
{
  Int_t n = fCorrections.GetEntriesFast();
  for (Int_t id = 0; id < n; id++) { 
    const Correction* c = GetCorrection(id);
    if (what.EqualTo(c->GetName(), TString::kIgnoreCase)) return id;
  }
  return -1;
}

//____________________________________________________________________
void
AliCorrectionManagerBase::EnableCorrection(Int_t id, Bool_t enable)
{
  Correction* c = GetCorrection(id);
  if (!c) { 
    AliWarningF("Cannot enable non-existing correction at %d", id);
    return;
  }
  c->fEnabled = enable;
}

//____________________________________________________________________
Int_t
AliCorrectionManagerBase::GetId(const TObject* obj) const
{
  Int_t   n   = fCorrections.GetEntriesFast();
  TClass* ocl = obj->IsA();
  for (Int_t id = 0; id < n; id++) { 
    const Correction* c = GetCorrection(id);
    if (ocl->InheritsFrom(c->fCls)) return id;
  }
  return -1;
}
//____________________________________________________________________
TObject*
AliCorrectionManagerBase::Get(Int_t id)
{
  Correction* c = GetCorrection(id);
  if (!c) {
    AliWarningF("Cannot find correction with id %d", id);
    return 0;
  }
  return c->Get();
}
//____________________________________________________________________
const TObject*
AliCorrectionManagerBase::Get(Int_t id) const
{
  const Correction* c = GetCorrection(id);
  if (!c) {
    AliWarningF("Cannot find correction with id %d", id);
    return 0;
  }
  return c->Get();
}

//____________________________________________________________________
Bool_t
AliCorrectionManagerBase::InitCorrections(ULong_t    run, 
					  UShort_t   sys, 
					  UShort_t   sNN, 
					  Short_t    fld, 
					  Bool_t     mc, 
					  Bool_t     sat,
					  Bool_t     force)
{
  if (force) fIsInit = false;
  if (!CheckConditions(run, sys, sNN, fld, mc, sat)) return false;
  if (!ReadCorrections(run, sys, sNN, fld, mc, sat)) return false;
  fIsInit = true;

  if (fDB) {
    delete fDB;
    fDB = 0;
  }

  return true;
}

//____________________________________________________________________
Bool_t
AliCorrectionManagerBase::CheckConditions(ULong_t    run, 
					  UShort_t   sys, 
					  UShort_t   sNN, 
					  Short_t    fld, 
					  Bool_t     mc, 
					  Bool_t     sat)
{
  if (!fIsInit) return true;

  AliInfo("We are already initialised - checking settings...");
  Bool_t same = true;
  if (fRun != run) {
    same = false;
  }
  if (fSys != sys) { 
    AliWarningF("Initialised collision system %s (%d) and "
		"passed same %s (%d) does not match", 
		AliForwardUtil::CollisionSystemString(fSys), fSys,
		AliForwardUtil::CollisionSystemString(sys), sys);
    same = false;
  }
  if (TMath::Abs(fSNN - sNN) >= 10) {
    AliWarningF("Initialised center of mass energy per nuclean "
		"%s (%d) and passed same %s (%d) does not match",
		AliForwardUtil::CenterOfMassEnergyString(fSNN), fSNN,
		AliForwardUtil::CenterOfMassEnergyString(sNN), sNN);
    same = false;
  }
  if (fField != fld) {
      AliWarningF("Initialied L3 magnetic field %s (%d) and passed "
		  "same %s (%d) does not match", 
		  AliForwardUtil::MagneticFieldString(fField), fField,
		  AliForwardUtil::MagneticFieldString(fld), fld);
      same = false;
  }
  if (fMC != mc) {
    AliWarningF("Initialied data type (%s) and passed "
		"same (%s) does not match", 
		(fMC ? "MC" : "real"), (mc ? "MC" : "real"));
    same = false;
  }
  if (fSatellite != sat) {
    AliWarningF("Initialied collision ip type (%s) and passed "
		"same (%s) does not match", 
		(fSatellite ? "satellite" : "nominal"), 
		(sat ? "satellite" : "nominal"));
    same = false;
  }
  if (!same) {
    AliWarning("Intialised parameters and these are not the same " 
	       "- PROCEED WITH CAUTION!");
  }
  else
    AliInfo("Initialized values consistent with data");
  
  return true;

}

//____________________________________________________________________
Bool_t
AliCorrectionManagerBase::ReadCorrection(Int_t      id,
					 ULong_t    run, 
					 UShort_t   sys, 
					 UShort_t   sNN, 
					 Short_t    fld, 
					 Bool_t     mc, 
					 Bool_t     sat)
{
  if (!fDB) {
    // We should always open the database, since we're not
    // streamingthat object to disk.
    fDB = new AliOADBForward;
  }

  Correction* c = GetCorrection(id);
  if (!c->fEnabled) return true;
  return c->ReadIt(fDB, run, sys, sNN, fld, mc, sat, fDebug, fFallBack);
}

//____________________________________________________________________
Bool_t
AliCorrectionManagerBase::ReadCorrections(ULong_t    run, 
					  UShort_t   sys, 
					  UShort_t   sNN, 
					  Short_t    fld, 
					  Bool_t     mc, 
					  Bool_t     sat)
{
  if (fIsInit) return true;
  if (fRun       == run && 
      fSys       == sys && 
      fField     == fld && 
      fMC        == mc  && 
      fSatellite == sat &&
      TMath::Abs(fSNN - sNN) < 11) { 
    // Already initialized for this - return
    fIsInit = true;
    return true;
  }
  if (!fDB) {
    // We should always open the database, since we're not
    // streamingthat object to disk.
    fDB = new AliOADBForward;
  }

  fRun       = run;
  fSys       = sys; 
  fSNN       = sNN;
  fField     = fld;
  fMC        = mc;
  fSatellite = sat;
  Int_t  n   = fCorrections.GetEntriesFast();
  Bool_t ret = true;
  for (Int_t id = 0; id < n; id++) 
    if (!ReadCorrection(id, run, sys, sNN, fld, mc, sat)) ret = false;
  return ret;
}

//====================================================================
AliCorrectionManagerBase::Correction::Correction() 
  : TNamed(), 
    fCls(0), 
    fClientCls(""),
    fQueryFields(0), 
    fEnabled(false), 
    fLastEntry(),
    fObject(0)
{}

//____________________________________________________________________
AliCorrectionManagerBase::Correction::Correction(const TString& tableName, 
						 const TString& fileName, 
						 TClass*        cls,
						 UShort_t       fields,
						 Bool_t         enabled) 
  : TNamed(tableName, fileName), 
    fCls(cls), 
    fClientCls(cls->GetName()),
    fQueryFields(fields), 
    fEnabled(enabled), 
    fLastEntry(""),
    fObject(0)
{}

//____________________________________________________________________
AliCorrectionManagerBase::Correction::Correction(const Correction& o)
  : TNamed(o), 
    fCls(o.fCls), 
    fClientCls(o.fClientCls),
    fQueryFields(o.fQueryFields),
    fEnabled(o.fEnabled), 
    fLastEntry(o.fLastEntry),
    fObject(o.fObject)
{}

//____________________________________________________________________
AliCorrectionManagerBase::Correction&
AliCorrectionManagerBase::Correction::operator=(const Correction& o)
{
  if (&o == this) return *this;
  SetName(o.GetName());
  SetTitle(o.GetTitle());
  fCls       	   = o.fCls;
  //fClientCls       = o.fClientCls;
  fQueryFields     = o.fQueryFields;
  fEnabled   	   = o.fEnabled;
  fLastEntry 	   = o.fLastEntry;
  fObject    	   = o.fObject;
  return *this;
}

//____________________________________________________________________
void
AliCorrectionManagerBase::Correction::MassageFields(ULong_t&  run,
						    UShort_t& sys,
						    UShort_t& sNN, 
						    Short_t & fld, 
						    Bool_t&   mc, 
						    Bool_t&   sat) const
{
  // Massage fields according to settings 
  if (!(fQueryFields & kRun))       run = kIgnoreValue;
  if (!(fQueryFields & kSys))       sys = kIgnoreValue;
  if (!(fQueryFields & kSNN))       sNN = kIgnoreValue;
  if (!(fQueryFields & kField))     fld = AliOADBForward::kInvalidField; 
  if (!(fQueryFields & kMC))        mc  = false;
  if (!(fQueryFields & kSatellite)) sat = false;
}
//____________________________________________________________________
void
AliCorrectionManagerBase::Correction::CorrectFields(ULong_t&  run,
						    UShort_t& sys,
						    UShort_t& sNN, 
						    Short_t & fld, 
						    Bool_t&   mc, 
						    Bool_t&   sat) const
{
  // Massage fields according to settings
  if (run <= 197388 && run >= 196433) {
    // These are Pb-p runs
    sys = 4; // AliAODForwardMult::kPbp;
  }
  if (sNN == 27615) sNN = 2760;
}

//____________________________________________________________________
Bool_t
AliCorrectionManagerBase::Correction::OpenIt(AliOADBForward* db, 
					     Bool_t vrb, 
					     Bool_t fallback) const
{
  if (!db) {
    Warning("OpenIt", "No DB passed");
    return false;
  }

  // Check if table is open
  if (db->FindTable(fName, true)) return true;

  //  if not try to open it 
  if (db->Open(fTitle, fName, false, vrb, fallback)) return true;

  // Give warning 
  AliWarningF("Failed to open table %s from %s", GetName(), GetTitle());
  AliWarningF("content of %s for %s:", 
	      gSystem->WorkingDirectory(), GetName());
  gSystem->Exec("pwd; ls -l");
  return false;
}
 
//____________________________________________________________________
Bool_t
AliCorrectionManagerBase::Correction::ReadIt(AliOADBForward* db, 
					     ULong_t         run, 
					     UShort_t        sys, 
					     UShort_t        sNN, 
					     Short_t         fld, 
					     Bool_t          mc, 
					     Bool_t          sat,
					     Bool_t          vrb,
					     Bool_t          fallback)
{
  if (!fEnabled) {
    AliWarningF("Correction %s not enabled", GetName());
    return 0;
  }

  // Assume failure 
  fObject = 0;

  // Massage fields according to settings 
  MassageFields(run, sys, sNN, fld, mc, sat);

  if (!OpenIt(db, vrb, fallback)) return false;

  // Query database 
  AliOADBForward::Entry* e = db->Get(fName, run, AliOADBForward::kDefault, 
				     sys, sNN, fld, mc, sat);
  // Check return value 
  if (!e || !e->fData) {
    AliWarningF("Failed to get %s from database in %s with "
		"run=%lu sys=%hu sNN=%hu fld=%hd %s %s", 
		GetName(), GetTitle(), run, sys, sNN, fld, 
		(mc ? "MC" : "real"), (sat ? "satellite" : "nominal"));
    return false;
  }

  // Ge the returned data
  TObject* o = e->fData;
  if (!o) {
    AliWarningF("Entry %p \"%s\" has no data! (%p)", e, e->GetTitle(), o);
    return false;
  }
  const TClass* cl = TheClass();
  // Check return class 
  if (!o->IsA()->InheritsFrom(cl)) { 
    AliWarningF("%p is not pointer to a %s object but a %s", 
		o, fCls->GetName(), o->ClassName());
    return false;
  }

  // Success 
  fObject    = o;
  fLastEntry = e->GetTitle();

  return true;
}

//____________________________________________________________________
Bool_t
AliCorrectionManagerBase::Correction::StoreIt(AliOADBForward* db, 
					      TObject*        obj,
					      ULong_t         run, 
					      UShort_t        sys, 
					      UShort_t        sNN, 
					      Short_t         fld, 
					      Bool_t          mc, 
					      Bool_t          sat,
					      const char*     file, 
					      const char*     meth) const
{
  // Info("StoreIt", "Storing run=%lu sys=%hy sNN=%d fld=%d mc=%d sat=%d", 
  //       run, sys, sNN, fld, mc, sat);
  const TClass* cl = TheClass();

  // Check value class 
  if (!obj->IsA()->InheritsFrom(cl)) { 
    AliWarningF("%p is not pointer to a %s object but a %s", 
		obj, cl->GetName(), obj->ClassName());
    return false;
  }

  Bool_t          local    = file || !db;
  TString         fileName = (local ? file : fTitle.Data());
  AliOADBForward* tdb      = (local ? new AliOADBForward : db);
  
  // Try to open the table read/write 
  if (!tdb->Open(fileName, Form("%s/%s", GetName(), meth), true, false)) {
    AliWarningF("Failed to open table %s in %s", GetName(), fileName.Data());
    return false;
  }

  // Massage fields according to settings 
  MassageFields(run, sys, sNN, fld, mc, sat);
  
  // Try to insert the object 
  if (!tdb->Insert(fName, obj, run, sys, sNN, fld, mc, sat)) { 
    AliWarningF("Failed to insert into %s off database in %s with "
		"run=%lu sys=%hu sNN=%hu fld=%hd %s %s", 
		GetName(), GetTitle(), run, sys, sNN, fld, 
		(mc ? "MC" : "real"), (sat ? "satellite" : "nominal"));
    return false;
  }

  if (local) { 
    tdb->Close();
    delete tdb;

    AliInfoF("Correction object %s written to DB in %s - merge this with "
	     "%s to store for good", obj->GetName(), fileName.Data(), 
	     GetTitle());
  }

  // Success 
  return true;
}
//____________________________________________________________________
TObject*
AliCorrectionManagerBase::Correction::Get()   
{
  if (!fEnabled) {
    AliWarningF("Correction %s not enabled", GetName());
    return 0;
  }
  return fObject;
}
//____________________________________________________________________
const TObject*
AliCorrectionManagerBase::Correction::Get() const
{
  if (!fEnabled) {
    AliWarningF("Correction %s not enabled", GetName());
    return 0;
  }
  return fObject;
}

//____________________________________________________________________
const TClass*
AliCorrectionManagerBase::Correction::TheClass() const
{
  if (fCls) return fCls;
  if (fClientCls.IsNull()) { 
    AliErrorF("No class name set for correction %s", GetName());
    return 0;
  }
  fCls = gROOT->GetClass(fClientCls);
  if (!fCls) { 
    AliErrorF("Couldn't get class %s for correction %s", 
	      fClientCls.Data(), GetName());
    return 0;
  }
  return fCls;
}

//____________________________________________________________________
void
AliCorrectionManagerBase::Correction::Print(Option_t* option) const
{
  AliForwardUtil::PrintTask(*this);
  gROOT->IncreaseDirLevel();
  PFB("Enabled", fEnabled);
  if (!fEnabled) {
    gROOT->DecreaseDirLevel();
    return;
  }

  TString flds;
  if (fQueryFields & kRun)       flds.Append("run");
  if (fQueryFields & kSys)       flds.Append("|sys");
  if (fQueryFields & kSNN)       flds.Append("|sNN");
  if (fQueryFields & kField)     flds.Append("|field");
  if (fQueryFields & kMC)        flds.Append("|MC");
  if (fQueryFields & kSatellite) flds.Append("|Satellite");
  if (flds.BeginsWith("|")) flds.Remove(0,1);

  const TClass* cl = TheClass();

  PFV("Path",GetTitle());
  PFV("Data class", cl->GetName());
  PFV("Query fields", flds);
  
  if (fObject && !fLastEntry.IsNull())  PF("Entry", fLastEntry);
  
  TString opt(option);
  opt.ToUpper();
  if (opt.Contains("D") && fObject) {
    gROOT->IncreaseDirLevel();
    fObject->Print();
    gROOT->DecreaseDirLevel();
  }
  gROOT->DecreaseDirLevel();
}

//____________________________________________________________________
void
AliCorrectionManagerBase::Correction::Browse(TBrowser* b)
{
  b->Add(const_cast<TClass*>(fCls), "Class");
  TString flds;
  if (fQueryFields & kRun)       flds.Append("run");
  if (fQueryFields & kSys)       flds.Append("|sys");
  if (fQueryFields & kSNN)       flds.Append("|sNN");
  if (fQueryFields & kField)     flds.Append("|field");
  if (fQueryFields & kMC)        flds.Append("|MC");
  if (fQueryFields & kSatellite) flds.Append("|Satellite");
  if (flds.BeginsWith("|")) flds.Remove(0,1);

  b->Add(new TObjString(flds), "Query fields");
  b->Add(new TParameter<bool>("Enabled", fEnabled));
  b->Add(new TObjString(fLastEntry), "Entry");
  if (fObject) b->Add(fObject);
}
//____________________________________________________________________
Bool_t
AliCorrectionManagerBase::Correction::CleanIt(AliOADBForward* db,
					      const TString& dest,
					      Bool_t verb,
					      Bool_t all) const
{
  // Open the table for this correction 
  if (!OpenIt(db, verb , false)) {
    Warning("CleanIt", "Failed to open table for %s", GetName());
    return false;
  }

  // Get our table - should be here if the above succeeded
  AliOADBForward::Table* t = db->FindTable(fName, !verb);
  if (!t) {
    Warning("CleanIt", "Failed to get table for %s", GetName());
    return false;
  }

  // Get some pointers and make a bit mask of entries to copy
  TTree* tree = t->fTree;
  Int_t  nEnt = tree->GetEntries();
  Int_t  nCpy = 0;
  TBits  copy(nEnt);
  copy.ResetAllBits();

  // Loop over all entries 
  Info("CleanIt", "Looping over %d entries in tree", nEnt);
  for (Int_t i = 0; i < nEnt; i++) {
    if (all) {
      copy.SetBitNumber(i, true);
      continue;
    }
    
    // Read in next entry 
    tree->GetEntry(i);

    // Check if we got an object 
    AliOADBForward::Entry* e = t->fEntry;    
    if (!e) continue;
    
    // Let's see it 
    // if (verb) e->Print();

    // Now query the DB with this entry's fields 
    ULong_t  run = e->fRunNo;
    UShort_t sys = e->fSys;
    UShort_t sNN = e->fSNN;
    Short_t  fld = e->fField;
    Bool_t   mc  = e->fMC;
    Bool_t   sat = e->fSatellite;
    TString  txt = e->GetTitle();
    MassageFields(run, sys, sNN, fld, mc, sat);
    CorrectFields(run, sys, sNN, fld, mc, sat);
    
    Int_t r = t->GetEntry(run, AliOADBForward::kDefault, 
			  sys, sNN, fld, mc, sat);
    if (r < 0) { 
      Warning("CleanIt","Uh! didn't get an entry for %s (%d)",
	      txt.Data(), i);
      r = i;
      // continue;
    }
#if 0
    // Here, we hard-code some fixes to remove bad entries 
    if (fld != AliOADBForward::kInvalidField) {
      switch (sys) {
      case 1:
	if      (sNN ==  900 &&  fld <= 0) r = -1;
	else if (sNN == 7000 &&  fld >= 0) r = -1;
	else if (sNN == 2760 &&  fld >= 0) r = -1;
	break;
      case 2:
	if (fld >= 0) r = -1;
	break;
      }
    }
#endif

    Printf("%-10s (%3d %3d): %s %s", 
	   (i==r ? "copied" : "ignored"), i, r, txt.Data(), GetName());

    if (r != i) continue;

    // If the entry found by the query and this entry is the same, 
    // then we should keep it 
    copy.SetBitNumber(i,true);
    // runs.Append(Form("%7lu", run));
  }
  
  // Now loop over the entries of the tree again, this time 
  // checking if we should copy or not. 
  // Loop over all entries 
  for (Int_t i = 0; i < nEnt; i++) { 
    // If not marked for copy, continue to the next. 
    if (!copy.TestBitNumber(i)) continue; 

    // Read in next entry 
    tree->GetEntry(i);

    // Check if we got an object 
    AliOADBForward::Entry* e = t->fEntry;    
    if (!e) continue;
  
    // Let's get the entry data and store that in a new correction
    // table
    ULong_t  run = e->fRunNo;
    UShort_t sys = e->fSys;
    UShort_t sNN = e->fSNN;
    Short_t  fld = e->fField;
    Bool_t   mc  = e->fMC;
    Bool_t   sat = e->fSatellite;
    TObject* obj = e->fData;
    TString  txt = e->GetTitle();
    CorrectFields(run, sys, sNN, fld, mc, sat);
    Printf("Storing %3d: %s -> %9lu/%d/%05d/%2d/%s/%s %s",
	   i, txt.Data(), run, sys, sNN, fld, (mc?"simu":"real"),
	   (sat?"sat":"nom"), GetName());
    if (!StoreIt(0, obj, run, sys, sNN, fld, mc, sat, dest.Data(),
		 AliOADBForward::Mode2String(t->fMode))) {
      Warning("CleanIt", "Failed to write new entry to %s", dest.Data());
      continue;
    }
    nCpy++;
  }
  // Printf("Runs: %s", runs.Data());
  Printf("Copied %6d entries of %6d for %s", nCpy, nEnt, fName.Data());
  return true;
}

//
// EOF
//
