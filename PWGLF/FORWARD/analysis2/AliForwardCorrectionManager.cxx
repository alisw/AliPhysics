//
// Manager (singleton) of corrections 
// 
#include "AliForwardCorrectionManager.h"
#include "AliFMDCorrSecondaryMap.h"
#include "AliFMDCorrDoubleHit.h"
#include "AliFMDCorrELossFit.h"
#include "AliFMDCorrVertexBias.h"
#include "AliFMDCorrMergingEfficiency.h"
#include "AliFMDCorrAcceptance.h"
#include "AliFMDCorrNoiseGain.h"
#include "AliForwardUtil.h"
#include "AliOADBForward.h"
#include <TString.h>
#include <AliLog.h>
#include <TFile.h>
#include <TSystem.h>
#include <TBrowser.h>
#include <TROOT.h>
#include <TClass.h>
#include <iostream>
#include <iomanip>
    
//____________________________________________________________________
AliForwardCorrectionManager* AliForwardCorrectionManager::fgInstance= 0;
const char* AliForwardCorrectionManager::fgkSecondaryMapSkel = "secondary";
const char* AliForwardCorrectionManager::fgkDoubleHitSkel    = "doublehit";
const char* AliForwardCorrectionManager::fgkELossFitsSkel    = "elossfits";
const char* AliForwardCorrectionManager::fgkVertexBiasSkel   = "vertexbias";
const char* AliForwardCorrectionManager::fgkMergingEffSkel   = "merging";
const char* AliForwardCorrectionManager::fgkAcceptanceSkel   = "acceptance";
const char* AliForwardCorrectionManager::fgkNoiseGainSkel    = "noisegain";

#define DB_NAME "fmd_corrections.root"

//____________________________________________________________________
AliForwardCorrectionManager& AliForwardCorrectionManager::Instance()
{
  // 
  // Access to the singleton object 
  // 
  // Return:
  //    Reference to the singleton object 
  //
  if (!fgInstance) fgInstance= new AliForwardCorrectionManager(false);
  return *fgInstance;
}

//____________________________________________________________________
AliForwardCorrectionManager::AliForwardCorrectionManager()
{
  // 
  // Default constructor 
  //
}
//____________________________________________________________________
AliForwardCorrectionManager::AliForwardCorrectionManager(Bool_t d)
  : AliCorrectionManagerBase(d)
{
  // 
  // Non-default constructor
  // 
  // Parameters:
  //    Not used
  //
  RegisterCorrection(kIdSecondaryMap,
		     fgkSecondaryMapSkel,
		     fgkSecondaryMapSkel,
		     AliFMDCorrSecondaryMap::Class(), 
		     kStandard|kSatellite);
  RegisterCorrection(kIdELossFits,
		     fgkELossFitsSkel,
		     fgkELossFitsSkel,
		     AliFMDCorrELossFit::Class(), 
		     kRun|kSys|kSNN|kSatellite|kMC /*kFull*/);
  RegisterCorrection(kIdVertexBias,
		     fgkVertexBiasSkel,
		     fgkVertexBiasSkel,
		     AliFMDCorrVertexBias::Class(), 
		     kStandard|kSatellite);
  RegisterCorrection(kIdMergingEfficiency,
		     fgkMergingEffSkel,
		     fgkMergingEffSkel,
		     AliFMDCorrMergingEfficiency::Class(), 
		     kStandard|kSatellite);
  RegisterCorrection(kIdDoubleHit,
		     fgkDoubleHitSkel,
		     fgkDoubleHitSkel,
		     AliFMDCorrDoubleHit::Class(),
		     kStandard|kMC);
  RegisterCorrection(kIdAcceptance,
		     fgkAcceptanceSkel,
		     fgkAcceptanceSkel,
		     AliFMDCorrAcceptance::Class(),
		     kRun|kSys|kSNN|kSatellite);
  RegisterCorrection(kIdNoiseGain,
		     fgkNoiseGainSkel,
		     fgkNoiseGainSkel,
		     AliFMDCorrNoiseGain::Class(),
		     kRun);
}
//____________________________________________________________________
Bool_t
AliForwardCorrectionManager::Init(ULong_t     runNo, 
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
AliForwardCorrectionManager::Init(ULong_t  runNo, 
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
  EnableCorrection(kIdSecondaryMap,	what & kSecondaryMap);
  EnableCorrection(kIdDoubleHit,	what & kDoubleHit);
  EnableCorrection(kIdELossFits,	what & kELossFits);
  EnableCorrection(kIdAcceptance,	what & kAcceptance);
  EnableCorrection(kIdVertexBias,	what & kVertexBias);
  EnableCorrection(kIdMergingEfficiency,what & kMergingEfficiency);
  EnableCorrection(kIdNoiseGain,	what & kNoiseGain);
  
  return InitCorrections(runNo, sys, sNN, field, mc, sat, force);
}

//____________________________________________________________________
UInt_t
AliForwardCorrectionManager::ParseFields(const TString& fields)
{
  UInt_t      ret    = 0;
  TObjArray*  tokens = fields.Tokenize(" \t,|+:;-&");
  TIter       next(tokens);
  TObjString* ostr = 0;
  while ((ostr = static_cast<TObjString*>(next()))) {
    const TString& str = ostr->String();
    
    if (str.Contains("all", TString::kIgnoreCase)) 
      ret |= kAll;
    else if (str.Contains("default", TString::kIgnoreCase)) 
      ret |= kDefault;
    else if (str.Contains(fgkSecondaryMapSkel, TString::kIgnoreCase))
      ret |= kSecondaryMap;
    else if (str.Contains(fgkDoubleHitSkel, TString::kIgnoreCase))
      ret |= kDoubleHit;
    else if (str.Contains(fgkELossFitsSkel, TString::kIgnoreCase))
      ret |= kELossFits;
    else if (str.Contains(fgkVertexBiasSkel, TString::kIgnoreCase))
      ret |= kVertexBias;
    else if (str.Contains(fgkMergingEffSkel, TString::kIgnoreCase))
      ret |= kMergingEfficiency;
    else if (str.Contains(fgkAcceptanceSkel, TString::kIgnoreCase))
      ret |= kAcceptance;
    else if (str.Contains(fgkNoiseGainSkel, TString::kIgnoreCase))
      ret |= kNoiseGain;
    else 
      AliWarningClassF("Unknown correction: %s", str.Data());
  }
  delete tokens;
  return ret;
}


//____________________________________________________________________
const AliFMDCorrELossFit*
AliForwardCorrectionManager::GetELossFit() const 
{
  /** 
   * Get the energy loss fit correction object. 
   * 
   * @return Get the energy loss fits corrections object or null pointer
   */
  return static_cast<const AliFMDCorrELossFit*>(Get(kIdELossFits)); 
}
//____________________________________________________________________
const AliFMDCorrSecondaryMap*
AliForwardCorrectionManager::GetSecondaryMap() const 
{
  /** 
   * Get the secondary correction map
   * 
   * @return Get the secondary correction map object or null
   */
  return static_cast<const AliFMDCorrSecondaryMap*>(Get(kIdSecondaryMap)); 
}
//____________________________________________________________________
const AliFMDCorrDoubleHit*
AliForwardCorrectionManager::GetDoubleHit() const 
{
  /** 
   * Get the double hit correction object
   * 
   * @return Get the double hit correction object or null 
   */
  return static_cast<const AliFMDCorrDoubleHit*>(Get(kIdDoubleHit)); 
}
//____________________________________________________________________
const AliFMDCorrVertexBias*
AliForwardCorrectionManager::GetVertexBias() const 
{
  /** 
   * Get the vertex bias correction object
   * 
   * @return Get the vertex bias correction object or null 
   */
  return static_cast<const AliFMDCorrVertexBias*>(Get(kIdVertexBias)); 
}
//____________________________________________________________________
const AliFMDCorrMergingEfficiency*
AliForwardCorrectionManager::GetMergingEfficiency() const 
{
  /** 
   * Get the merging efficiency 
   * 
   * 
   * @return Get the vertex efficiency correction 
   */
  return 
    static_cast<const AliFMDCorrMergingEfficiency*>(Get(kIdMergingEfficiency)); 
}
//____________________________________________________________________
const AliFMDCorrAcceptance*
AliForwardCorrectionManager::GetAcceptance() const 
{
  /** 
   * Get the acceptance correction due to dead channels 
   * 
   * 
   * @return Acceptance correction due to dead channels 
   */
  return static_cast<const AliFMDCorrAcceptance*>(Get(kIdAcceptance)); 
}
//____________________________________________________________________
const AliFMDCorrNoiseGain*
AliForwardCorrectionManager::GetNoiseGain() const 
{
  /** 
   * Get the noisegain calibration
   * 
   * @return NoiseGain calibration
   */
  return static_cast<const AliFMDCorrNoiseGain*>(Get(kIdNoiseGain)); 
}

//____________________________________________________________________
const TAxis* 
AliForwardCorrectionManager::GetEtaAxis() const
{
  const AliFMDCorrSecondaryMap* map = GetSecondaryMap();
  if (!map) return 0;
  return &(map->GetEtaAxis());
}
//____________________________________________________________________
const TAxis* 
AliForwardCorrectionManager::GetVertexAxis() const
{
  const AliFMDCorrSecondaryMap* map = GetSecondaryMap();
  if (!map) return 0;
  return &(map->GetVertexAxis());
}


#ifndef DOXY_INPUT
//______________________________________________________________________________
// @cond
void AliForwardCorrectionManager::Streamer(TBuffer &R__b)
{
  //
  // Stream an object of class AliForwardCorrectionManager.
  //
  if (R__b.IsReading()) {
     R__b.ReadClassBuffer(AliForwardCorrectionManager::Class(),this);
     if (fgInstance) {
       AliWarning(Form("Singleton instance already set (%p) when reading "
		       "singleton object (%p).  Read object will be new "
		       "singleton object", fgInstance, this));
       // delete fgInstance;
     }
     fgInstance = this;
     // fgInstance->fCorrections.ls();
  } else {
    R__b.WriteClassBuffer(AliForwardCorrectionManager::Class(),this);
  }
}
// @endcond
#endif

//____________________________________________________________________
//
// EOF
//
