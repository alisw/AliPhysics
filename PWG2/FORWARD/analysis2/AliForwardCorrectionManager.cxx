//
// Manager (singleton) of corrections 
// 
#include "AliForwardCorrectionManager.h"
#include "AliForwardUtil.h"
#include <TString.h>
#include <AliLog.h>
#include <TFile.h>
#include <TSystem.h>

    
//____________________________________________________________________
AliForwardCorrectionManager* AliForwardCorrectionManager::fgInstance = 0;
const char* AliForwardCorrectionManager::fgkSecondaryMapSkel = "secondary";
const char* AliForwardCorrectionManager::fgkDoubleHitSkel    = "doublehit";
const char* AliForwardCorrectionManager::fgkELossFitsSkel    = "elossfits";
const char* AliForwardCorrectionManager::fgkVertexBiasSkel   = "vertexbias";
const char* AliForwardCorrectionManager::fgkMergingEffSkel   = "merging";

#define PREFIX "$(ALICE_ROOT)/PWG2/FORWARD/corrections/"

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
    fELossFitsPath(PREFIX "ELossFits"),
    fMergingEffPath(PREFIX "MergingEfficiency"), 
    fSecondaryMapPath(PREFIX "SecondaryMap"),
    fDoubleHitPath(PREFIX "DoubleHit"),
    fVertexBiasPath(PREFIX "VertexBias"),
    fELossFit(0),
    fSecondaryMap(0),
    fDoubleHit(0),
    fVertexBias(0),
    fMergingEfficiency(0)
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
    fELossFit(o.fELossFit),
    fSecondaryMap(o.fSecondaryMap),
    fDoubleHit(o.fDoubleHit),
    fVertexBias(o.fVertexBias),
    fMergingEfficiency(o.fMergingEfficiency)

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
  fELossFit         = o.fELossFit;
  fSecondaryMap     = o.fSecondaryMap;
  fDoubleHit        = o.fDoubleHit;
  fVertexBias       = o.fVertexBias;
  fMergingEfficiency= o.fMergingEfficiency;
  return *this;
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
  if (fInit) return kTRUE;

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
  // Read event selection efficiencies if requested 
  if (what & kVertexBias) {
    if (!ReadVertexBias(cms, sNN, field)) {
      AliWarning(Form("Failed to read in event selection efficiency for "
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
//
// EOF
//
