//
// Manager (singleton) of corrections 
// 
#include "AliForwardCorrectionManagerOADB.h"
#include "AliFMDCorrSecondaryMap.h"
#include "AliFMDCorrDoubleHit.h"
#include "AliFMDCorrELossFit.h"
#include "AliFMDCorrVertexBias.h"
#include "AliFMDCorrMergingEfficiency.h"
#include "AliFMDCorrAcceptance.h"
#include "AliForwardUtil.h"
#include "AliOADBForward.h"
#include <TString.h>
#include <AliLog.h>
#include <TFile.h>
#include <TSystem.h>
#include <TBrowser.h>
#include <TROOT.h>
#include <iostream>
#include <iomanip>
    
//____________________________________________________________________
AliForwardCorrectionManagerOADB* AliForwardCorrectionManagerOADB::fgInstance= 0;
const char* AliForwardCorrectionManagerOADB::fgkSecondaryMapSkel = "secondary";
const char* AliForwardCorrectionManagerOADB::fgkDoubleHitSkel    = "doublehit";
const char* AliForwardCorrectionManagerOADB::fgkELossFitsSkel    = "elossfits";
const char* AliForwardCorrectionManagerOADB::fgkVertexBiasSkel   = "vertexbias";
const char* AliForwardCorrectionManagerOADB::fgkMergingEffSkel   = "merging";
const char* AliForwardCorrectionManagerOADB::fgkAcceptanceSkel   = "acceptance";

#define PREFIX  "$(ALICE_ROOT)/OADB/PWGLF/FORWARD/CORRECTIONS/data/"
#define DB_NAME "fmd_corrections.root"

//____________________________________________________________________
AliForwardCorrectionManagerOADB& AliForwardCorrectionManagerOADB::Instance()
{
  // 
  // Access to the singleton object 
  // 
  // Return:
  //    Reference to the singleton object 
  //
  if (!fgInstance) fgInstance= new AliForwardCorrectionManagerOADB;
  return *fgInstance;
}

//____________________________________________________________________
AliForwardCorrectionManagerOADB::AliForwardCorrectionManagerOADB()
  : TObject(), 
    fInit(kFALSE),
    fRunNo(0),
    fSys(0),
    fSNN(0),
    fField(999),
    fMC(false), 
    fSat(false),
    fELossFitsPath(PREFIX DB_NAME),
    fMergingEffPath(PREFIX DB_NAME), 
    fSecondaryMapPath(PREFIX DB_NAME),
    fDoubleHitPath(PREFIX DB_NAME),
    fVertexBiasPath(PREFIX DB_NAME),
    fAcceptancePath(PREFIX DB_NAME),
    fELossFit(0),
    fSecondaryMap(0),
    fDoubleHit(0),
    fVertexBias(0),
    fMergingEfficiency(0),
    fAcceptance(0),
    fDB(0)
{
  // 
  // Default constructor 
  //
}
//____________________________________________________________________
AliForwardCorrectionManagerOADB::AliForwardCorrectionManagerOADB(const AliForwardCorrectionManagerOADB& o)
  : TObject(o),
    fInit(o.fInit),
    fRunNo(o.fRunNo),
    fSys(o.fSys),
    fSNN(o.fSNN),
    fField(o.fField),
    fMC(o.fMC),
    fSat(o.fSat),
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
  fAcceptance(o.fAcceptance),
  fDB(0)
{
  // 
  // Copy constructor 
  // 
  // Parameters:
  //    o Object to copy from 
  //
}
//____________________________________________________________________
AliForwardCorrectionManagerOADB&
AliForwardCorrectionManagerOADB::operator=(const AliForwardCorrectionManagerOADB& o)
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
  fRunNo            = o.fRunNo;
  fSys              = o.fSys;
  fSNN              = o.fSNN;
  fField            = o.fField;
  fMC               = o.fMC;
  fSat              = o.fSat;
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
  fDB               = o.fDB;
  return *this;
}

//____________________________________________________________________
void
AliForwardCorrectionManagerOADB::SetPrefix(const char* prefix)
{
  /** 
   *
   * @param prefix Prefix to correction objects. 
   */
  fELossFitsPath    = Form("%s/%s", prefix, DB_NAME);
  fMergingEffPath   = Form("%s/%s", prefix, DB_NAME); 
  fSecondaryMapPath = Form("%s/%s", prefix, DB_NAME);
  fDoubleHitPath    = Form("%s/%s", prefix, DB_NAME);
  fVertexBiasPath   = Form("%s/%s", prefix, DB_NAME);
  fAcceptancePath   = Form("%s/%s", prefix, DB_NAME);
  
}
//____________________________________________________________________
void
AliForwardCorrectionManagerOADB::SetFile(ECorrection what, const char* filename)
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
  *path = filename;
}

//____________________________________________________________________
Bool_t
AliForwardCorrectionManagerOADB::Init(ULong_t     runNo, 
				      const char* sys, 
				      Float_t     sNN, 
				      Float_t     field,
				      Bool_t      mc,
				      Bool_t      sat,
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
  UShort_t col = AliForwardUtil::ParseCollisionSystem(sys);
  // AliInfo(Form("Initialising with cms='%s', sNN=%fGeV field=%fkG", 
  //	       cms, sNN, field));
  return Init(runNo, col, 
	      AliForwardUtil::ParseCenterOfMassEnergy(col, sNN),
	      AliForwardUtil::ParseMagneticField(field), 
	      mc, sat, what, force);
}

//____________________________________________________________________
Bool_t
AliForwardCorrectionManagerOADB::Init(ULong_t  runNo, 
				      UShort_t sys, 
				      UShort_t sNN, 
				      Short_t  field,
				      Bool_t   mc,
				      Bool_t   sat,
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
  if (!fDB) { 
    // We should always open the database, since we're not
    // streamingthat object to disk.
    fDB = new AliOADBForward;
  }
  if (fInit) {
    // Check that the initialisation and the passed parameters 
    // match - if not give an error but continue - this allows 
    // users to override particular settings. 
    
    AliInfo("We are already initialised - checking settings...");
    Bool_t same = true;
    if (fRunNo != runNo) {
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
    if (fField != field) {
      AliWarningF("Initialied L3 magnetic field %s (%d) and passed "
		  "same %s (%d) does not match", 
		  AliForwardUtil::MagneticFieldString(fField), fField,
		  AliForwardUtil::MagneticFieldString(field), field);
      same = false;
    }
    if (fMC != mc) {
      AliWarningF("Initialied data type (%s) and passed "
		  "same (%s) does not match", 
		  (fMC ? "MC" : "real"), (mc ? "MC" : "real"));
      same = false;
    }
    if (fSat != sat) {
      AliWarningF("Initialied collision ip type (%s) and passed "
		  "same (%s) does not match", 
		  (fSat ? "satelitte" : "nominal"), 
		  (sat ? "satellite" : "nominal"));
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
  if (fRunNo == runNo && 
      fSys   == sys   && 
      fField == field &&
      fMC    == mc    && 
      fSat   == sat   &&
      TMath::Abs(fSNN - sNN) < 10) {
    // We're already initialised for these settings - do nothing and return
    fInit = kTRUE;
    return ret;
  }
  // Set cached parameters 
  fRunNo = runNo;
  fSys   = sys;
  fSNN   = sNN;
  fField = field;
  fMC    = mc;
  fSat   = sat;

  // AliInfo(Form("Initialising with cms=%d, sNN=%dGeV field=%dkG", 
  //  	       cms, sNN, field));
  // Read secondary map if requested 
  if (what & kSecondaryMap) {
    if (!ReadSecondaryMap(runNo, sys, sNN, field, sat)) {
      AliWarningF("Failed to read in secondary map for "
		  "run=%lu, sys=%hu, sNN=%huGeV, field=%hdkG, %s", 
		  runNo, sys, sNN, field, (sat ? "satellite" : "nominal"));
      ret = kFALSE;
    }
  }
  // Read double hit if requested 
  if (what & kDoubleHit) {
    if (!ReadDoubleHit(runNo, sys, sNN, field, sat)) {
      AliWarningF("Failed to read in double hit correction for "
		  "run=%lu, sys=%hu, sNN=%huGeV, field=%hdkG, %s", 
		  runNo, sys, sNN, field, (sat ? "satellite" : "nominal"));
      ret = kFALSE;
    }
  }
  // Read energy loss fits if requested 
  if (what & kELossFits) {
    if (!ReadELossFits(runNo, sys, sNN, field, mc, sat)) {
      AliWarningF("Failed to read in energy loss fits for "
		  "run=%lu, sys=%hu, sNN=%huGeV, field=%hdkG, %s, %s", 
		  runNo, sys, sNN, field, mc ? "MC" : "real", 
		  (sat ? "satellite" : "nominal"));
      ret = kFALSE;
    }
  }
  // Read acceptance correction if requested 
  if (what & kAcceptance) {
    if (!ReadAcceptance(runNo, sys, sNN, sat)) {
      AliWarningF("Failed to read in acceptance for "
		  "run=%lu, sys=%hu, sNN=%huGeV, %s", 
		  runNo, sys, sNN, (sat ? "satellite" : "nominal"));
      ret = kFALSE;
    }
  }
  // Read event selection efficiencies if requested 
  if (what & kVertexBias) {
    if (!ReadVertexBias(runNo, sys, sNN, field, sat)) {
      AliWarningF("Failed to read in vertex bias correction for "
		  "run=%lu, sys=%hu, sNN=%huGeV, field=%hdkG, %s", 
		  runNo, sys, sNN, field, (sat ? "satellite" : "nominal"));
      ret = kFALSE;
    }
  }
  // Read merging efficiencies if requested 
  if (what & kMergingEfficiency) {
    if (!ReadMergingEfficiency(runNo, sys, sNN, field, sat)) {
      AliWarningF("Failed to read in hit merging efficiency for "
		  "run=%lu, sys=%hu, sNN=%huGeV, field=%hdkG, %s", 
		  runNo, sys, sNN, field, (sat ? "satellite" : "nominal"));
      ret = kFALSE;
    }
  }
  fInit = kTRUE;
  return ret;
}
//____________________________________________________________________
const TAxis* 
AliForwardCorrectionManagerOADB::GetEtaAxis() const
{
  if (!fSecondaryMap) return 0;
  return &(fSecondaryMap->GetEtaAxis());
}
//____________________________________________________________________
const TAxis* 
AliForwardCorrectionManagerOADB::GetVertexAxis() const
{
  if (!fSecondaryMap) return 0;
  return &(fSecondaryMap->GetVertexAxis());
}

//____________________________________________________________________
const Char_t* 
AliForwardCorrectionManagerOADB::GetFileName(ECorrection what) const
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
  return gSystem->BaseName(GetFilePath(what));
}
//____________________________________________________________________
const Char_t*
AliForwardCorrectionManagerOADB::GetFileDir(ECorrection what) const
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
  return gSystem->DirName(GetFilePath(what));
}

//____________________________________________________________________
const TString&
AliForwardCorrectionManagerOADB::GetFilePath(ECorrection what) const
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
  switch (what) {    
  case kSecondaryMap:        return fSecondaryMapPath;
  case kDoubleHit:           return fDoubleHitPath;
  case kELossFits:           return fELossFitsPath;
  case kVertexBias:          return fVertexBiasPath;
  case kMergingEfficiency:   return fMergingEffPath;
  case kAcceptance:          return fAcceptancePath;
  default: break;
  }
  static TString null;
  return null;
}


//____________________________________________________________________
const Char_t*
AliForwardCorrectionManagerOADB::GetObjectName(ECorrection what) const
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
  switch (what) {
  case kSecondaryMap:       return fgkSecondaryMapSkel;
  case kDoubleHit:          return fgkDoubleHitSkel;
  case kELossFits:          return fgkELossFitsSkel;
  case kVertexBias:         return fgkVertexBiasSkel;
  case kMergingEfficiency:  return fgkMergingEffSkel;
  case kAcceptance:         return fgkAcceptanceSkel;
  default: break;
  }
  return 0;
}

//____________________________________________________________________
const TClass*
AliForwardCorrectionManagerOADB::GetObjectClass(ECorrection what) const
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
  switch (what) {
  case kSecondaryMap:       return AliFMDCorrSecondaryMap::Class();
  case kDoubleHit:          return AliFMDCorrDoubleHit::Class();
  case kELossFits:          return AliFMDCorrELossFit::Class();
  case kVertexBias:         return AliFMDCorrVertexBias::Class();
  case kMergingEfficiency:  return AliFMDCorrMergingEfficiency::Class();
  case kAcceptance:         return AliFMDCorrAcceptance::Class();
  default: break;
  }
  return 0;
}
//____________________________________________________________________
AliForwardCorrectionManagerOADB::ECorrection
AliForwardCorrectionManagerOADB::GetObjectType(const TString& what) const
{
  TString w(what);
  w.ToLower();
  
  if      (w.EqualTo(GetObjectName(kSecondaryMap)))       
    return kSecondaryMap;
  else if (w.EqualTo(GetObjectName(kDoubleHit)))          
    return kDoubleHit;
  else if (w.EqualTo(GetObjectName(kELossFits)))          
    return kELossFits;
  else if (w.EqualTo(GetObjectName(kVertexBias)))         
    return kVertexBias;
  else if (w.EqualTo(GetObjectName(kMergingEfficiency)))  
    return kMergingEfficiency;
  else if (w.EqualTo(GetObjectName(kAcceptance)))         
    return kAcceptance;
  
  return kAll;
}
//____________________________________________________________________
AliForwardCorrectionManagerOADB::ECorrection
AliForwardCorrectionManagerOADB::GetObjectType(const TObject* obj) const
{
  TClass* cl = obj->IsA();
  
  if      (cl->InheritsFrom(GetObjectClass(kSecondaryMap)))       
    return kSecondaryMap;
  else if (cl->InheritsFrom(GetObjectClass(kDoubleHit)))          
    return kDoubleHit;
  else if (cl->InheritsFrom(GetObjectClass(kELossFits)))          
    return kELossFits;
  else if (cl->InheritsFrom(GetObjectClass(kVertexBias)))         
    return kVertexBias;
  else if (cl->InheritsFrom(GetObjectClass(kMergingEfficiency)))  
    return kMergingEfficiency;
  else if (cl->InheritsFrom(GetObjectClass(kAcceptance)))
    return kAcceptance;
  
  return kAll;
}

//____________________________________________________________________
TObject*
AliForwardCorrectionManagerOADB::GetObject(ECorrection what, 
					   ULong_t     runNo,
					   UShort_t    sys, 
					   UShort_t    sNN, 
					   Short_t     field,
					   Bool_t      mc,
					   Bool_t      sat) const
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
  if (!fDB) {
    AliWarning("Database not opened");
    return 0;
  }
  TString tableName = GetObjectName(what);
  if (!fDB->FindTable(tableName, true)) {
    if (!fDB->Open(GetFilePath(what), tableName, false, true)) { 
      AliWarningF("Failed to open table %s from %s", 
		  tableName.Data(), GetFilePath(what).Data());
      return 0;
    }
  }

  TObject* o = fDB->GetData(GetObjectName(what), runNo,
			    AliOADBForward::kDefault, 
			    sys, sNN, field, mc, sat);
  
  return o;
}
//____________________________________________________________________
TObject*
AliForwardCorrectionManagerOADB::GetObject(ECorrection what) const
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
  
  return GetObject(what, fRunNo, fSys, fSNN, fField, fMC, fSat);
}

#define CHECK_TYPE(O,RET,CL) do {	\
  RET = dynamic_cast<CL*>(O);		\
  if (!O) AliWarningF("%p is not a pointer to a %s object, but a %s", \
		      O, #CL, O->ClassName()); } while (false)

//____________________________________________________________________
Bool_t 
AliForwardCorrectionManagerOADB::ReadSecondaryMap(ULong_t  runNo, 
						  UShort_t sys, 
						  UShort_t sNN, 
						  Short_t  field,
						  Bool_t   sat)
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

  TObject* o = GetObject(kSecondaryMap, runNo, sys, sNN, field, false, sat);
  if (!o) return kFALSE;

  CHECK_TYPE(o, fSecondaryMap, AliFMDCorrSecondaryMap);
  if (!fSecondaryMap) return kFALSE;

  return kTRUE;
}
//____________________________________________________________________
Bool_t 
AliForwardCorrectionManagerOADB::ReadDoubleHit(ULong_t  runNo,
					       UShort_t sys, 
					       UShort_t sNN, 
					       Short_t  field, 
					       Bool_t   sat)
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

  TObject* o = GetObject(kDoubleHit, runNo, sys, sNN, field, false, sat);
  if (!o) return kFALSE;

  CHECK_TYPE(o, fDoubleHit, AliFMDCorrDoubleHit);
  if (!fDoubleHit) return kFALSE;

  return kTRUE;
}

//____________________________________________________________________
Bool_t 
AliForwardCorrectionManagerOADB::ReadELossFits(ULong_t  runNo,
					       UShort_t sys, 
					       UShort_t sNN, 
					       Short_t  field, 
					       Bool_t   mc,
					       Bool_t   sat)
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

  TObject* o = GetObject(kELossFits, runNo, sys, sNN, field, mc, sat);
  if (!o) return kFALSE;

  CHECK_TYPE(o, fELossFit, AliFMDCorrELossFit);
  if (!fELossFit) return kFALSE;

  return kTRUE;
}

//____________________________________________________________________
Bool_t 
AliForwardCorrectionManagerOADB::ReadVertexBias(ULong_t  runNo, 
						UShort_t sys, 
						UShort_t sNN, 
						Short_t  field, 
						Bool_t   sat)
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

  TObject* o = GetObject(kVertexBias, runNo, sys, sNN, field, false, sat);
  if (!o) return kFALSE;

  CHECK_TYPE(o, fVertexBias, AliFMDCorrVertexBias);
  if (!fVertexBias) return kFALSE;

  return kTRUE;
}

//____________________________________________________________________
Bool_t 
AliForwardCorrectionManagerOADB::ReadMergingEfficiency(ULong_t runNo, 
						       UShort_t sys, 
						       UShort_t sNN, 
						       Short_t  field, 
						       Bool_t   sat)
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

  TObject* o = GetObject(kMergingEfficiency, runNo, sys, sNN, field, false,sat);
  if (!o) return kFALSE;

  CHECK_TYPE(o, fMergingEfficiency, AliFMDCorrMergingEfficiency);
  if (!fMergingEfficiency) return kFALSE;

  return kTRUE;
}

//____________________________________________________________________
Bool_t 
AliForwardCorrectionManagerOADB::ReadAcceptance(ULong_t  runNo, 
						UShort_t sys, 
						UShort_t sNN, 
						Bool_t   sat)
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

  TObject* o = GetObject(kAcceptance, runNo, sys, sNN, 0, false, sat);
  if (!o) return kFALSE;

  CHECK_TYPE(o, fAcceptance, AliFMDCorrAcceptance);
  if (!fAcceptance) return kFALSE;

  return kTRUE;
}
//____________________________________________________________________
void
AliForwardCorrectionManagerOADB::Print(Option_t* option) const
{
  // 
  // Print stuff 
  //
  char ind[gROOT->GetDirLevel()+1];
  for (Int_t i = 0; i < gROOT->GetDirLevel(); i++) ind[i] = ' ';
  ind[gROOT->GetDirLevel()] = '\0';

  std::cout << ind << "AliForwardCorrectionManagerOADB:\n"
	    << ind << "  Initialised:      " 
	    << (fInit ? "yes" : "no") << std::endl;
  if (fInit) 
    std::cout << ind << "  Run number:       " << fRunNo << "\n"
	      << ind << "  Collision system: " 
	      << AliForwardUtil::CollisionSystemString(fSys) << "\n"
	      << ind << "  Sqrt(s_NN):       "
	      << AliForwardUtil::CenterOfMassEnergyString(fSNN) << "\n"
	      << ind << "  Magnetic field:   " 
	      << AliForwardUtil::MagneticFieldString(fField) << "\n"
	      << ind << "  For simulations:  " << (fMC ? "yes" : "no") << "\n"
	      << ind << "  For satellites:   " << (fSat ? "yes" : "no") << "\n"
	      << std::endl;
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
AliForwardCorrectionManagerOADB::Browse(TBrowser* b)
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

//____________________________________________________________________
Bool_t
AliForwardCorrectionManagerOADB::StoreObject(ULong_t        runNo,
					     UShort_t       sys, 
					     UShort_t       sNN, 
					     Short_t        fld, 
					     Bool_t         mc,
					     Bool_t         sat,
					     TObject*       obj, 
					     Bool_t         full,
					     const char*    meth) const
{
  ECorrection what = GetObjectType(obj);
  if (what == kAll) { 
    AliErrorF("Cannot deduce the correction type from object of class %s", 
	      obj->ClassName());
    return false;
  }
  return StoreObject(what, runNo, sys, sNN, fld, mc, sat, obj, full, meth);
}
//____________________________________________________________________
Bool_t
AliForwardCorrectionManagerOADB::StoreObject(const TString& what, 
					     ULong_t        runNo,
					     UShort_t       sys, 
					     UShort_t       sNN, 
					     Short_t        fld, 
					     Bool_t         mc,
					     Bool_t         sat,
					     TObject*       obj, 
					     Bool_t         full,
					     const char*    meth) const
{
  return StoreObject(GetObjectType(what), runNo, sys, 
		     sNN, fld, mc, sat, obj, full, meth);
}

//____________________________________________________________________
Bool_t
AliForwardCorrectionManagerOADB::StoreObject(ECorrection what, 
					     ULong_t     runNo,
					     UShort_t    sys, 
					     UShort_t    sNN, 
					     Short_t     fld, 
					     Bool_t      mc,
					     Bool_t      sat,
					     TObject*    obj, 
					     Bool_t      full,
					     const char* meth) const
{
  // 
  // Write correction output to (a temporary) file 
  // 
  // Parameters: 
  //   What     What to write 
  //   sys      Collision system (1: pp, 2: PbPb)
  //   sNN      Center of mass energy per nucleon (GeV)
  //   fld      Field (kG)
  //   mc       MC-only flag 
  //   obj      Object to write 
  //   full     if true, write to full path, otherwise locally
  // 
  // Return: 
  //   true on success. 
  TString         tableName = GetObjectName(what);
  Bool_t          local     = !(full && fDB);
  TString         fileName  = (local ? DB_NAME : GetFilePath(what));
  AliOADBForward* db        = (local ? new AliOADBForward : fDB);
  ECorrection     otype     = GetObjectType(obj);
  if (otype != what) { 
    AliErrorF("Correction type 0x%02x (%s) and type of object 0x%02x (%s) "
	      "does not match", what, tableName.Data(), otype, 
	      obj->ClassName());
    return false;
  }

  if (!db->Open(fileName, Form("%s/%s", tableName.Data(), meth), true, true)) {
    Error("StoreObject", 
	  "Failed to open table %s/%s in %s for read+write (%s)", 
	  tableName.Data(), meth, fileName.Data(), 
	  (local ? "local" : "global"));
    return false;
  }

  // Filter out fields we do not need
  switch (what) {
  case kELossFits:                      break; // Need all fields
  case kDoubleHit:         sat = false; break; // Depreacted
  case kAcceptance:        fld = 0;            // Fall through
  case kSecondaryMap:                          // Fall through
  case kVertexBias:                            // Fall through
  case kMergingEfficiency: mc  = false; break; // No specific for MC
  default:                              break; // Never here 
  }
    

  if (!db->Insert(tableName, obj, runNo, sys, sNN, fld, mc, sat)) { 
    AliErrorF("Failed to write %s to database for "
	      "table=%s run=%lu, sys=%hu, sNN=%hu, field=%hd, mc=%d, sat=%d", 
	      obj->GetName(), tableName.Data(), runNo, sys, sNN, 
	      fld, mc, sat);
    return false;
  }
    
  if (local) {
    db->Close();
    delete db;
    
    AliInfoF("Correction object %s written to DB in %s - merge this with "
	     "%s to store for good", obj->GetName(), DB_NAME, 
	     GetFilePath(what).Data());
  }
  return true;
}

#ifndef DOXY_INPUT
//______________________________________________________________________________
void AliForwardCorrectionManagerOADB::Streamer(TBuffer &R__b)
{
  //
  // Stream an object of class AliForwardCorrectionManagerOADB.
  //
  if (R__b.IsReading()) {
     R__b.ReadClassBuffer(AliForwardCorrectionManagerOADB::Class(),this);
     if (fgInstance) 
       AliWarning(Form("Singleton instance already set (%p) when reading "
		       "singleton object (%p).  Read object will be new "
		       "singleton object", fgInstance, this));
     fgInstance = this;
  } else {
    R__b.WriteClassBuffer(AliForwardCorrectionManagerOADB::Class(),this);
  }
}
#endif
#if 0
//______________________________________________________________________________
void AliForwardCorrectionManagerOADB::Streamer(TBuffer &R__b)
{
   // Stream an object of class AliForwardCorrectionManagerOADB.

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
      R__b.CheckByteCount(R__s, R__c, AliForwardCorrectionManagerOADB::IsA());
   } else {
      R__c = R__b.WriteVersion(AliForwardCorrectionManagerOADB::IsA(), kTRUE);
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
