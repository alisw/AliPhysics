//
// Manager (singleton) of corrections 
// 
#include "AliForwardCorrectionManager.h"
#include "AliFMDCorrDoubleHit.h"
#include "AliFMDCorrELossFit.h"
#include "AliFMDCorrVertexBias.h"
#include "AliFMDCorrMergingEfficiency.h"
#include "AliFMDCorrAcceptance.h"
#include "AliForwardUtil.h"
#include <TString.h>
#include <AliLog.h>
#include <TFile.h>
#include <TSystem.h>
#include <TBrowser.h>
#include <TROOT.h>
#include <iostream>
#include <iomanip>
    
//____________________________________________________________________
AliForwardCorrectionManager* AliForwardCorrectionManager::fgInstance = 0;
const char* AliForwardCorrectionManager::fgkSecondaryMapSkel = "secondary";
const char* AliForwardCorrectionManager::fgkDoubleHitSkel    = "doublehit";
const char* AliForwardCorrectionManager::fgkELossFitsSkel    = "elossfits";
const char* AliForwardCorrectionManager::fgkVertexBiasSkel   = "vertexbias";
const char* AliForwardCorrectionManager::fgkMergingEffSkel   = "merging";
const char* AliForwardCorrectionManager::fgkAcceptanceSkel   = "acceptance";

#define PREFIX "$(ALICE_ROOT)/PWG2/FORWARD/corrections/"
#define ELOSSFIT_DIR   "ELossFits"
#define MERGING_DIR    "MergingEfficiency"
#define SECONDARY_DIR  "SecondaryMap"
#define DOUBLE_DIR     "DoubleHit"
#define VERTEX_DIR     "VertexBias"
#define ACCEPTANCE_DIR "Acceptance"

//____________________________________________________________________
AliForwardCorrectionManager& AliForwardCorrectionManager::Instance()
{
  // 
  // Access to the singleton object 
  // 
  // Return:
  //    Reference to the singleton object 
  //
  if (!fgInstance) fgInstance= new AliForwardCorrectionManager;
  return *fgInstance;
}

//____________________________________________________________________
AliForwardCorrectionManager::AliForwardCorrectionManager()
  : TObject(), 
    fInit(kFALSE),
    fSys(0),
    fSNN(0),
    fField(999),
    fELossFitsPath(PREFIX ELOSSFIT_DIR),
    fMergingEffPath(PREFIX MERGING_DIR), 
    fSecondaryMapPath(PREFIX SECONDARY_DIR),
    fDoubleHitPath(PREFIX DOUBLE_DIR),
    fVertexBiasPath(PREFIX VERTEX_DIR),
    fAcceptancePath(PREFIX ACCEPTANCE_DIR),
    fELossFit(0),
    fSecondaryMap(0),
    fDoubleHit(0),
    fVertexBias(0),
    fMergingEfficiency(0),
    fAcceptance(0)
{
  // 
  // Default constructor 
  //
}
//____________________________________________________________________
AliForwardCorrectionManager::AliForwardCorrectionManager(const AliForwardCorrectionManager& o)
  : TObject(o),
    fInit(o.fInit),
    fSys(o.fSys),
    fSNN(o.fSNN),
    fField(o.fField),
    fELossFitsPath(o.fELossFitsPath),
    fMergingEffPath(o.fMergingEffPath), 
    fSecondaryMapPath(o.fSecondaryMapPath),
    fDoubleHitPath(o.fDoubleHitPath),
    fVertexBiasPath(o.fVertexBiasPath),
    fAcceptancePath(o.fAcceptancePath),
    fELossFit(o.fELossFit),
    fSecondaryMap(o.fSecondaryMap),
    fDoubleHit(o.fDoubleHit),
    fVertexBias(o.fVertexBias),
    fMergingEfficiency(o.fMergingEfficiency),
    fAcceptance(o.fAcceptance)

{
  // 
  // Copy constructor 
  // 
  // Parameters:
  //    o Object to copy from 
  //
}
//____________________________________________________________________
AliForwardCorrectionManager&
AliForwardCorrectionManager::operator=(const AliForwardCorrectionManager& o)
{
  // 
  // Assignment operator 
  // 
  // Parameters:
  //    o Object to assign from 
  // 
  // Return:
  //    Reference to this object 
  //
  fInit             = o.fInit;
  fSys              = o.fSys;
  fSNN              = o.fSNN;
  fField            = o.fField;
  fELossFitsPath    = o.fELossFitsPath;
  fMergingEffPath   = o.fMergingEffPath;
  fSecondaryMapPath = o.fSecondaryMapPath;
  fDoubleHitPath    = o.fDoubleHitPath;
  fVertexBiasPath   = o.fVertexBiasPath;
  fAcceptancePath   = o.fAcceptancePath;
  fELossFit         = o.fELossFit;
  fSecondaryMap     = o.fSecondaryMap;
  fDoubleHit        = o.fDoubleHit;
  fVertexBias       = o.fVertexBias;
  fMergingEfficiency= o.fMergingEfficiency;
  fAcceptance       = o.fAcceptance;
  return *this;
}

//____________________________________________________________________
void
AliForwardCorrectionManager::SetPrefix(const char* prefix)
{
  /** 
   *
   * @param prefix Prefix to correction objects. 
   */
  fELossFitsPath    = Form("%s/%s", prefix, ELOSSFIT_DIR);
  fMergingEffPath   = Form("%s/%s", prefix, MERGING_DIR); 
  fSecondaryMapPath = Form("%s/%s", prefix, SECONDARY_DIR);
  fDoubleHitPath    = Form("%s/%s", prefix, DOUBLE_DIR);
  fVertexBiasPath   = Form("%s/%s", prefix, VERTEX_DIR);
  fAcceptancePath   = Form("%s/%s", prefix, ACCEPTANCE_DIR);
  
}
//____________________________________________________________________
void
AliForwardCorrectionManager::SetFileDir(ECorrection what, const char* dir)
{
  /** 
   * Set the file directory for a type 
   * 
   * @param what     Type 
   * @param dirname  Directory name 
   */
  TString *path = 0;
  if      (what & kSecondaryMap)        path = &fSecondaryMapPath;
  else if (what & kDoubleHit)           path = &fDoubleHitPath;
  else if (what & kELossFits)           path = &fELossFitsPath;
  else if (what & kVertexBias)          path = &fVertexBiasPath;
  else if (what & kMergingEfficiency)   path = &fMergingEffPath;
  else if (what & kAcceptance)          path = &fAcceptancePath;
  else { 
    AliWarning(Form("No such path defined for 0x%02x", what));
    return;
  }
  if (!path) {
    AliWarning(Form("Couldn't find string for path 0x%02x", what));
    return;
  }
  *path = dir;
}

//____________________________________________________________________
Bool_t
AliForwardCorrectionManager::Init(const char* cms, 
				  Float_t     sNN, 
				  Float_t     field,
				  Bool_t      mc,
				  UInt_t      what,
				  Bool_t      force)
{
  // 
  // Read in correction based on passed parameters
  // 
  // Parameters:
  //    collisionSystem Collision system string 
  //    cmsNN           Center of mass energy per nucleon pair [GeV]
  //    field           Magnetic field [kG]
  //    mc              Monte-carlo switch
  //    what            What to read in 
  //    force           Force (re-)reading of specified things
  // 
  // Return:
  //    true on success
  //
  UShort_t col = AliForwardUtil::ParseCollisionSystem(cms);
  // AliInfo(Form("Initialising with cms='%s', sNN=%fGeV field=%fkG", 
  //	       cms, sNN, field));
  return Init(col, 
	      AliForwardUtil::ParseCenterOfMassEnergy(col, sNN),
	      AliForwardUtil::ParseMagneticField(field), 
	      mc, what, force);
}

//____________________________________________________________________
Bool_t
AliForwardCorrectionManager::Init(UShort_t cms, 
				  UShort_t sNN, 
				  Short_t  field,
				  Bool_t   mc,
				  UInt_t   what,
				  Bool_t   force)
{
  // 
  // Read in corrections based on the parameters given 
  // 
  // Parameters:
  //    collisionSystem Collision system
  //    cmsNN           Center of mass energy per nuclean pair [GeV]
  //    field           Magnetic field setting [kG]
  //    mc              Monte-carlo switch
  //    what            What to read in. 
  //    force           Force (re-)reading of specified things
  // 
  // Return:
  //    
  //
  if (force) fInit = kFALSE;
  if (fInit) {
    // Check that the initialisation and the passed parameters 
    // match - if not give an error but continue - this allows 
    // users to override particular settings. 
    
    AliInfo("We are already initialised - checking settings...");
    Bool_t same = true;
    if (fSys != cms) { 
      AliWarning(Form("Initialised collision system %s (%d) and "
		      "passed same %s (%d) does not match", 
		      AliForwardUtil::CollisionSystemString(fSys), fSys,
		      AliForwardUtil::CollisionSystemString(cms), cms));
      same = false;
    }
    if (TMath::Abs(fSNN - sNN) >= 10) {
      AliWarning(Form("Initialised center of mass energy per nuclean "
		      "%s (%d) and passed same %s (%d) does not match",
		      AliForwardUtil::CenterOfMassEnergyString(fSNN), fSNN,
		      AliForwardUtil::CenterOfMassEnergyString(sNN), sNN));
      same = false;
    }
    if (fField != field) {
      AliWarning(Form("Initialied L3 magnetic field %s (%d) and passed "
		      "same %s (%d) does not match", 
		      AliForwardUtil::MagneticFieldString(fField), fField,
		      AliForwardUtil::MagneticFieldString(field), field));
      same = false;
    }
    if (!same) {
      AliWarning("Intialised parameters and these are not the same " 
		 "- PROCEED WITH CAUTION!");
    }
    else
      AliInfo("Initialized values consistent with data");
    
    return kTRUE;
  }

  Bool_t ret = kTRUE;
  if (fSys == cms && TMath::Abs(fSNN - sNN) < 10 && fField == field) { 
    // We're already initialised for these settings - do nothing and return
    fInit = kTRUE;
    return ret;
  }
  // Set cached parameters 
  fSys   = cms;
  fSNN   = sNN;
  fField = field;

  // AliInfo(Form("Initialising with cms=%d, sNN=%dGeV field=%dkG", 
  //  	       cms, sNN, field));
  // Read secondary map if requested 
  if (what & kSecondaryMap) {
    if (!ReadSecondaryMap(cms, sNN, field)) {
      AliWarning(Form("Failed to read in secondary map for "
		      "cms=%d, sNN=%dGeV, field=%dkG", cms, sNN, field));
      ret = kFALSE;
    }
  }
  // Read double hit if requested 
  if (what & kDoubleHit) {
    if (!ReadDoubleHit(cms, sNN, field)) {
      AliWarning(Form("Failed to read in double hit correction for "
		      "cms=%d, sNN=%dGeV, field=%dkG", cms, sNN, field));
      ret = kFALSE;
    }
  }
  // Read energy loss fits if requested 
  if (what & kELossFits) {
    if (!ReadELossFits(cms, sNN, field, mc)) {
      AliWarning(Form("Failed to read in energy loss fits for "
		      "cms=%d, sNN=%dGeV, field=%dkG, %s", 
		      cms, sNN, field, mc ? "MC" : "real"));
      ret = kFALSE;
    }
  }
  // Read acceptance correction if requested 
  if (what & kAcceptance) {
    if (!ReadAcceptance(cms, sNN, 0)) {
      AliWarning(Form("Failed to read in acceptance for "
		      "cms=%d, sNN=%dGeV, field=%dkG", cms, sNN, 0));
      ret = kFALSE;
    }
  }
  // Read event selection efficiencies if requested 
  if (what & kVertexBias) {
    if (!ReadVertexBias(cms, sNN, field)) {
      AliWarning(Form("Failed to read in vertex bias correction for "
		      "cms=%d, sNN=%dGeV, field=%dkG", cms, sNN, field));
      ret = kFALSE;
    }
  }
  // Read merging efficiencies if requested 
  if (what & kMergingEfficiency) {
    if (!ReadMergingEfficiency(cms, sNN, field)) {
      AliWarning(Form("Failed to read in hit merging efficiency for "
		      "cms=%d, sNN=%dGeV, field=%dkG", 
		      cms, sNN, field));
      ret = kFALSE;
    }
  }
  fInit = kTRUE;
  return ret;
}

//____________________________________________________________________
TString 
AliForwardCorrectionManager::GetFileName(ECorrection what, 
					 UShort_t    sys, 
					 UShort_t    sNN,
					 Short_t     field,
					 Bool_t      mc) const
{
  // 
  // Get the path to the specified object 
  // 
  // Parameters:
  //    what  Which stuff to get the path for 
  //    sys   Collision system
  //    sNN   Center of mass energy [GeV]
  //    field Magnetic field in the L3 magnet [kG]
  //    mc    Whether the correction objects should be valid for MC
  // 
  // Return:
  //    The full path or null 
  //
  TString fname = "";
  fname = GetObjectName(what);
  fname.Append(Form("_%s_%04dGeV_%c%1dkG_%s.root", 
		    AliForwardUtil::CollisionSystemString(sys), 
		    sNN, (field < 0 ? 'm' : 'p'), TMath::Abs(field), 
		    (mc ? "MC" : "real")));
  return fname;
}
//____________________________________________________________________
TString
AliForwardCorrectionManager::GetFileName(ECorrection what) const
{
  // 
  // Get the file name of the specified object
  // 
  // Parameters:
  //    what Which stuff to get the path for 
  // 
  // Return:
  //    The full path or null
  //
  if (!fInit) { 
    AliWarning("Corrections manager initialised, do a forced Init(...)");
    return "";
  }
  return GetFileName(what, fSys, fSNN, fField, false);
}

//____________________________________________________________________
const Char_t*
AliForwardCorrectionManager::GetFileDir(ECorrection what) const
{
  // 
  // Get the path to the specified object 
  // 
  // Parameters:
  //    what  Which stuff to get the path for 
  // 
  // Return:
  //    The full path or null 
  //
  if      (what & kSecondaryMap)        return fSecondaryMapPath;
  else if (what & kDoubleHit)           return fDoubleHitPath;
  else if (what & kELossFits)           return fELossFitsPath;
  else if (what & kVertexBias)          return fVertexBiasPath;
  else if (what & kMergingEfficiency)   return fMergingEffPath;
  else if (what & kAcceptance)          return fAcceptancePath;

  AliWarning(Form("Unknown correction: %d", what));
  return 0;
}

//____________________________________________________________________
TString 
AliForwardCorrectionManager::GetFilePath(ECorrection what, 
					 UShort_t    sys, 
					 UShort_t    sNN,
					 Short_t     field,
					 Bool_t      mc) const
{
  // 
  // Get the path to the specified object 
  // 
  // Parameters:
  //    what  Which stuff to get the path for 
  //    sys   Collision system
  //    sNN   Center of mass energy [GeV]
  //    field Magnetic field in the L3 magnet [kG]
  //    mc    Whether the correction objects should be valid for MC
  // 
  // Return:
  //    The full path or null 
  //
  TString path = "";
  const Char_t* dir = GetFileDir(what);
  if (!dir) return path;
  
  TString fname(GetFileName(what, sys, sNN, field, mc));
  if (fname.IsNull()) return path;

  path = gSystem->ConcatFileName(gSystem->ExpandPathName(dir), fname);
  
  return path;
}
//____________________________________________________________________
TString
AliForwardCorrectionManager::GetFilePath(ECorrection what) const
{
  // 
  // Get the full path to the object.  Note, the manager must be
  // initialised for this to work
  // 
  // Parameters:
  //    what Which stuff to get the path for 
  // 
  // Return:
  //    The full path or null
  //
  if (!fInit) { 
    AliWarning("Corrections manager initialised, do a forced Init(...)");
    return "";
  }
  return GetFilePath(what, fSys, fSNN, fField, false);
}

//____________________________________________________________________
TFile*
AliForwardCorrectionManager::GetFile(ECorrection what, 
				     UShort_t    sys, 
				     UShort_t    sNN, 
				     Short_t     field, 
				     Bool_t      mc, 
				     Bool_t      rw, 
				     Bool_t      newfile) const
{
  // 
  // Open the file that contains the correction object specified 
  // 
  // Parameters:
  //    what  Which stuff to get the path for 
  //    sys   Collision system
  //    sNN   Center of mass energy [GeV]
  //    field Magnetic field in the L3 magnet [kG]
  //    mc    Whether the correction objects should be valid for MC
  //    rw    Whether to open the file in read/write
  //    newfile Wheter to make the file if it doesn't exist
  // 
  // Return:
  //    The file that contains the correction object or null 
  //
  TString path = GetFilePath(what, sys, sNN, field, mc);
  if (path.IsNull()) return 0;
  
  TString opt;
  if (newfile) opt="RECREATE";
  else {
    if (gSystem->AccessPathName(path.Data(), 
				(rw ? kWritePermission : kReadPermission))) {
      AliWarning(Form("file %s cannot be found or insufficient permissions", 
		      path.Data()));
      return 0;
    }
    opt=(rw ? "UPDATE" : "READ");
  }
  TFile* file = TFile::Open(path.Data(), opt.Data());
  if (!file) { 
    AliWarning(Form("file %s cannot be opened in mode %s", 
		    path.Data(), opt.Data()));
    return 0;
  }
  return file;
}
//____________________________________________________________________
TFile*
AliForwardCorrectionManager::GetFile(ECorrection what) const
{
  // 
  // Get the file that contains the object specifed.  Note, the manager
  // must be initialised for this to work. 
  // 
  // Parameters:
  //    what Which stuff to get the path for 
  // 
  // Return:
  //    The file that contains the correction object or null
  //
  if (!fInit) { 
    AliWarning("Corrections manager initialised, do a forced Init(...)");
    return 0;
  }
  return GetFile(what, fSys, fSNN, fField, false);
}

//____________________________________________________________________
const Char_t*
AliForwardCorrectionManager::GetObjectName(ECorrection what) const
{
  // 
  // Get the object name corresponding to correction type 
  // 
  // Parameters:
  //    what Correction 
  // 
  // Return:
  //    Object name or null
  //
  if      (what & kSecondaryMap)       return fgkSecondaryMapSkel;
  else if (what & kDoubleHit)          return fgkDoubleHitSkel;
  else if (what & kELossFits)          return fgkELossFitsSkel;
  else if (what & kVertexBias)         return fgkVertexBiasSkel;
  else if (what & kMergingEfficiency)  return fgkMergingEffSkel;
  else if (what & kAcceptance)         return fgkAcceptanceSkel;
  return 0;
}

//____________________________________________________________________
TObject*
AliForwardCorrectionManager::CheckObject(TFile* file, ECorrection what) const
{
  // 
  // Check if the specified objet exists in the file, and return it
  // 
  // Parameters:
  //    file File to query 
  //    what Correction type 
  // 
  // Return:
  //    Object found, or null
  //
  TObject* o = file->Get(GetObjectName(what));
  if (!o) { 
    AliWarning(Form("Object %s not found in %s", 
		    GetObjectName(what), file->GetName()));
    file->Close();
    return 0;
  }
  return o;
}
  
//____________________________________________________________________
TObject*
AliForwardCorrectionManager::GetObject(ECorrection what, 
				       UShort_t    sys, 
				       UShort_t    sNN, 
				       Short_t     field,
				       Bool_t      mc) const
{
  // 
  // Get the path to the specified object 
  // 
  // Parameters:
  //    what  Which stuff to get the path for 
  //    sys   Collision system
  //    sNN   Center of mass energy [GeV]
  //    field Magnetic field in the L3 magnet [kG]
  //    mc    Whether the correction objects should be valid for MC
  // 
  // Return:
  //    The full path or null 
  //
  TFile* file = GetFile(what, sys, sNN, field, mc, false, false);
  if (!file) return 0;
  
  return CheckObject(file, what);
}
//____________________________________________________________________
TObject*
AliForwardCorrectionManager::GetObject(ECorrection what) const
{
  // 
  // Get the object that contaisn the specified correction
  // 
  // Parameters:
  //    what Which object to get
  // 
  // Return:
  //    The object or null
  //
  if (!fInit) { 
    AliWarning("Corrections manager initialised, do a forced Init(...)");
    return 0;
  }
  return GetObject(what, fSys, fSNN, fField, false);
}


//____________________________________________________________________
Bool_t 
AliForwardCorrectionManager::ReadSecondaryMap(UShort_t sys, UShort_t sNN, 
					      Short_t field)
{
  // 
  // Read in the secondary map 
  // 
  // Parameters:
  //    sys   Collision system
  //    sNN   Center of mass energy [GeV]
  //    field Magnetic field in the L3 magnet [kG]
  // 
  // Return:
  //    True on success, false otherwise 
  //
  if (fInit) { 
    AliWarning("Corrections manager initialised, do a forced Init(...)");
    return kFALSE;
  }

  TObject* o = GetObject(kSecondaryMap, sys, sNN, field, false);
  if (!o) return kFALSE;

  fSecondaryMap = dynamic_cast<AliFMDCorrSecondaryMap*>(o);
  if (!fSecondaryMap) {
    AliWarning(Form("Object %s (%p) is not an AliFMDCorrSecondaryMap object, "
		    "but %s", fgkSecondaryMapSkel, o, o->ClassName())); 
    return kFALSE;
  }

  // file->Close();
  return kTRUE;
}
//____________________________________________________________________
Bool_t 
AliForwardCorrectionManager::ReadDoubleHit(UShort_t sys, UShort_t sNN, 
					   Short_t field)
{
  // 
  // Read in the double hit correction
  // 
  // Parameters:
  //    sys   Collision system
  //    sNN   Center of mass energy [GeV]
  //    field Magnetic field in the L3 magnet [kG]
  // 
  // Return:
  //    True on success, false otherwise 
  //
  if (fInit) { 
    AliWarning("Corrections manager initialised, do a forced Init(...)");
    return kFALSE;
  }

  TObject* o = GetObject(kDoubleHit, sys, sNN, field, false);
  if (!o) return kFALSE;

  fDoubleHit = dynamic_cast<AliFMDCorrDoubleHit*>(o);
  if (!fDoubleHit) {
    AliWarning(Form("Object %s (%p) is not an AliFMDCorrDoubleHit object, "
		    "but %s", fgkDoubleHitSkel, o, o->ClassName())); 
    return kFALSE;
  }

  // file->Close();
  return kTRUE;
}

//____________________________________________________________________
Bool_t 
AliForwardCorrectionManager::ReadELossFits(UShort_t sys, UShort_t sNN, 
					   Short_t field, Bool_t mc)
{
  // 
  // Read in the energy loss fits 
  // 
  // Parameters:
  //    sys   Collision system
  //    sNN   Center of mass energy [GeV]
  //    field Magnetic field in the L3 magnet [kG]
  //    mc    Whether the correction objects should be valid for MC
  // 
  // Return:
  //    True on success, false otherwise 
  //
  if (fInit) { 
    AliWarning("Corrections manager initialised, do a forced Init(...)");
    return kFALSE;
  }

  TObject* o = GetObject(kELossFits, sys, sNN, field, mc);
  if (!o) return kFALSE;

  fELossFit = dynamic_cast<AliFMDCorrELossFit*>(o);
  if (!fELossFit) {
    AliWarning(Form("Object %s (%p) is not an AliFMDCorrELossFit object, "
		    "but %s", fgkELossFitsSkel, o, o->ClassName()));
    return kFALSE;
  }

  // file->Close();
  return kTRUE;
}

//____________________________________________________________________
Bool_t 
AliForwardCorrectionManager::ReadVertexBias(UShort_t sys, 
					    UShort_t sNN, 
					    Short_t field)
{
  // 
  // Read in the event selection efficiency 
  // 
  // Parameters:
  //    sys   Collision system
  //    sNN   Center of mass energy [GeV]
  //    field Magnetic field in the L3 magnet [kG]
  // 
  // Return:
  //    True on success, false otherwise 
  //
  if (fInit) { 
    AliWarning("Corrections manager initialised, do a forced Init(...)");
    return kFALSE;
  }

  TObject* o = GetObject(kVertexBias, sys, sNN, field, false);
  if (!o) return kFALSE;

  fVertexBias = dynamic_cast<AliFMDCorrVertexBias*>(o);
  if (!fVertexBias) {
    AliWarning(Form("Object %s (%p) is not an AliFMDCorrVertexBias object, "
		    "but %s", fgkVertexBiasSkel, o, o->ClassName()));
    return kFALSE;
  }

  // file->Close();
  return kTRUE;
}

//____________________________________________________________________
Bool_t 
AliForwardCorrectionManager::ReadMergingEfficiency(UShort_t sys, 
						   UShort_t sNN, 
						   Short_t field)
{
  // 
  // Read in the merging efficiency 
  // 
  // Parameters:
  //    sys   Collision system
  //    sNN   Center of mass energy [GeV]
  //    field Magnetic field in the L3 magnet [kG]
  // 
  // Return:
  //    True on success, false otherwise 
  //
  if (fInit) { 
    AliWarning("Corrections manager initialised, do a forced Init(...)");
    return kFALSE;
  }

  TObject* o = GetObject(kMergingEfficiency, sys, sNN, field, false);
  if (!o) return kFALSE;

  fMergingEfficiency = dynamic_cast<AliFMDCorrMergingEfficiency*>(o);
  if (!fMergingEfficiency) {
    AliWarning(Form("Object %s (%p) is not an AliFMDCorrMergingEfficiency "
		    "object, but %s", fgkMergingEffSkel, o, o->ClassName()));
    return kFALSE;
  }

  // file->Close();
  return kTRUE;
}

//____________________________________________________________________
Bool_t 
AliForwardCorrectionManager::ReadAcceptance(UShort_t sys, 
					    UShort_t sNN, 
					    Short_t field)
{
  // 
  // Read in the event selection efficiency 
  // 
  // Parameters:
  //    sys   Collision system
  //    sNN   Center of mass energy [GeV]
  //    field Magnetic field in the L3 magnet [kG]
  // 
  // Return:
  //    True on success, false otherwise 
  //
  if (fInit) { 
    AliWarning("Corrections manager initialised, do a forced Init(...)");
    return kFALSE;
  }

  TObject* o = GetObject(kAcceptance, sys, sNN, field, false);
  if (!o) return kFALSE;

  fAcceptance = dynamic_cast<AliFMDCorrAcceptance*>(o);
  if (!fAcceptance) {
    AliWarning(Form("Object %s (%p) is not an AliFMDCorrAcceptance object, "
		    "but %s", fgkAcceptanceSkel, o, o->ClassName()));
    return kFALSE;
  }

  // file->Close();
  return kTRUE;
}
//____________________________________________________________________
void
AliForwardCorrectionManager::Print(Option_t* option) const
{
  // 
  // Print stuff 
  //
  char ind[gROOT->GetDirLevel()+1];
  for (Int_t i = 0; i < gROOT->GetDirLevel(); i++) ind[i] = ' ';
  ind[gROOT->GetDirLevel()] = '\0';

  std::cout << ind << "AliForwardCorrectionManager:\n"
	    << ind << "  Initialised:      " 
	    << (fInit ? "yes" : "no") << std::endl;
  if (fInit) 
    std::cout << ind << "  Collision system: " 
	      << AliForwardUtil::CollisionSystemString(fSys) << "\n"
	      << ind << "  Sqrt(s_NN):       "
	      << AliForwardUtil::CenterOfMassEnergyString(fSNN) << "\n"
	      << ind << "  Magnetic field:   " 
	      << AliForwardUtil::MagneticFieldString(fField) << std::endl;
  std::cout << ind << "  Paths:\n" 
	    << ind << "    ELoss Fits:     " << fELossFitsPath << "\n"
	    << ind << "    Merging eff.:   " << fMergingEffPath << "\n"
	    << ind << "    Secondary maps: " << fSecondaryMapPath << "\n"
	    << ind << "    2-hit corr.:    " << fSecondaryMapPath << "\n"
	    << ind << "    Vertex bias:    " << fVertexBiasPath << "\n"
	    << ind << "    Acceptance:     " << fAcceptancePath << std::endl;
  TString opt(option);
  opt.ToUpper();
  if (!opt.Contains("R")) return;
  
  gROOT->IncreaseDirLevel();
  if (fELossFit)	  fELossFit->Print(option);
  else 
    std::cout << ind << "  Energy loss fits  not initialised" << std::endl;
  
  if (fSecondaryMap)	  fSecondaryMap->Print(option);
  else 
    std::cout << ind << "  Secondary particle correction not initialised" 
	      << std::endl;

  if (fDoubleHit)	  fDoubleHit->Print(option);
  else 
    std::cout << ind << "  Double hit corr. not initialised" << std::endl;

  if (fVertexBias)	  fVertexBias->Print(option);
  else 
    std::cout << ind << "  Vertex bias correction not initialised" << std::endl;
  if (fMergingEfficiency) fMergingEfficiency->Print(option);
  else 
    std::cout << ind << "  Merging eff.  not initialised" << std::endl;

  if (fAcceptance)	  fAcceptance->Print(option);
  else 
    std::cout << ind << "  Acceptance corr.  not initialised" << std::endl;
  gROOT->DecreaseDirLevel();  
}

//____________________________________________________________________
void
AliForwardCorrectionManager::Browse(TBrowser* b)
{
  // 
  // Browse thos
  // 
  if (fELossFit)	  b->Add(fELossFit,          "Energy loss fits");
  if (fSecondaryMap)	  b->Add(fSecondaryMap,      "Secondary particle corr");
  if (fDoubleHit)	  b->Add(fDoubleHit,         "Double hit corr");
  if (fVertexBias)	  b->Add(fVertexBias,        "Vertex bias corr");
  if (fMergingEfficiency) b->Add(fMergingEfficiency, "Merging eff");
  if (fAcceptance)	  b->Add(fAcceptance,        "Acceptance corr");
}

#if 1
//______________________________________________________________________________
void AliForwardCorrectionManager::Streamer(TBuffer &R__b)
{
  //
  // Stream an object of class AliForwardCorrectionManager.
  //
  if (R__b.IsReading()) {
     R__b.ReadClassBuffer(AliForwardCorrectionManager::Class(),this);
     if (fgInstance) 
       AliWarning(Form("Singleton instance already set (%p) when reading "
		       "singleton object (%p).  Read object will be new "
		       "singleton object", fgInstance, this));
     fgInstance = this;
  } else {
    R__b.WriteClassBuffer(AliForwardCorrectionManager::Class(),this);
  }
}
#endif
#if 0
//______________________________________________________________________________
void AliForwardCorrectionManager::Streamer(TBuffer &R__b)
{
   // Stream an object of class AliForwardCorrectionManager.

   UInt_t R__s, R__c;
   if (R__b.IsReading()) {
      Version_t R__v = R__b.ReadVersion(&R__s, &R__c); if (R__v) { }
      TObject::Streamer(R__b);
      R__b >> fInit;
      R__b >> fSys;
      R__b >> fSNN;
      R__b >> fField;
      fELossFitsPath.Streamer(R__b);
      fMergingEffPath.Streamer(R__b);
      fSecondaryMapPath.Streamer(R__b);
      fDoubleHitPath.Streamer(R__b);
      fVertexBiasPath.Streamer(R__b);
      fAcceptancePath.Streamer(R__b);
      R__b >> fELossFit;
      R__b >> fSecondaryMap;
      R__b >> fDoubleHit;
      R__b >> fVertexBias;
      R__b >> fMergingEfficiency;
      R__b >> fAcceptance;
      R__b.CheckByteCount(R__s, R__c, AliForwardCorrectionManager::IsA());
   } else {
      R__c = R__b.WriteVersion(AliForwardCorrectionManager::IsA(), kTRUE);
      TObject::Streamer(R__b);
      R__b << fInit;
      R__b << fSys;
      R__b << fSNN;
      R__b << fField;
      fELossFitsPath.Streamer(R__b);
      fMergingEffPath.Streamer(R__b);
      fSecondaryMapPath.Streamer(R__b);
      fDoubleHitPath.Streamer(R__b);
      fVertexBiasPath.Streamer(R__b);
      fAcceptancePath.Streamer(R__b);
      R__b << fELossFit;
      R__b << fSecondaryMap;
      R__b << fDoubleHit;
      R__b << fVertexBias;
      R__b << fMergingEfficiency;
      R__b << fAcceptance;
      R__b.SetByteCount(R__c, kTRUE);
   }
}
#endif

//____________________________________________________________________
//
// EOF
//
