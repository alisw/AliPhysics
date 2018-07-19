//
// Manager (singleton) of corrections 
// 
#include "AliCentralCorrectionManager.h"
#include "AliCentralCorrSecondaryMap.h"
#include "AliCentralCorrAcceptance.h"
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
AliCentralCorrectionManager* AliCentralCorrectionManager::fgInstance= 0;
const char* AliCentralCorrectionManager::fgkSecondaryMapSkel = "secondary";
const char* AliCentralCorrectionManager::fgkAcceptanceSkel   = "acceptance";

#define DB_NAME "spd_corrections.root"

//____________________________________________________________________
AliCentralCorrectionManager& AliCentralCorrectionManager::Instance()
{
  // 
  // Access to the singleton object 
  // 
  // Return:
  //    Reference to the singleton object 
  //
  if (!fgInstance) fgInstance= new AliCentralCorrectionManager(false);
  return *fgInstance;
}

//____________________________________________________________________
AliCentralCorrectionManager::AliCentralCorrectionManager()
{
  // 
  // Default constructor 
  //
}
//____________________________________________________________________
AliCentralCorrectionManager::AliCentralCorrectionManager(Bool_t d)
  : AliCorrectionManagerBase(d)
{
  // 
  // Non-default constructor
  // 
  // Parameters:
  //    Not used
  //
  RegisterCorrection(kIdSecondaryMap, fgkSecondaryMapSkel, 
		     DB_NAME, AliCentralCorrSecondaryMap::Class(), 
		     kStandard|kSatellite);
  RegisterCorrection(kIdAcceptance, fgkAcceptanceSkel, 
		     DB_NAME, AliCentralCorrAcceptance::Class(),
		     kStandard|kSatellite);
}
//____________________________________________________________________
Bool_t
AliCentralCorrectionManager::Init(ULong_t     runNo, 
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
AliCentralCorrectionManager::Init(ULong_t  runNo, 
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
  EnableCorrection(kIdAcceptance,	what & kAcceptance);
  
  return InitCorrections(runNo, sys, sNN, field, mc, sat, force);
}

//____________________________________________________________________
const AliCentralCorrSecondaryMap*
AliCentralCorrectionManager::GetSecondaryMap() const 
{
  /** 
   * Get the secondary correction map
   * 
   * @return Get the secondary correction map object or null
   */
  return static_cast<const AliCentralCorrSecondaryMap*>(Get(kIdSecondaryMap)); 
}
//____________________________________________________________________
const AliCentralCorrAcceptance*
AliCentralCorrectionManager::GetAcceptance() const 
{
  /** 
   * Get the acceptance correction due to dead channels 
   * 
   * 
   * @return Acceptance correction due to dead channels 
   */
  return static_cast<const AliCentralCorrAcceptance*>(Get(kIdAcceptance)); 
}

//____________________________________________________________________
const TAxis* 
AliCentralCorrectionManager::GetVertexAxis() const
{
  const AliCentralCorrSecondaryMap* map = GetSecondaryMap();
  if (!map) return 0;
  return &(map->GetVertexAxis());
}


#ifndef DOXY_INPUT
//______________________________________________________________________________
// @cond 
void AliCentralCorrectionManager::Streamer(TBuffer &R__b)
{
  //
  // Stream an object of class AliCentralCorrectionManager.
  //
  if (R__b.IsReading()) {
     R__b.ReadClassBuffer(AliCentralCorrectionManager::Class(),this);
     if (fgInstance) {
       AliWarning(Form("Singleton instance already set (%p) when reading "
		       "singleton object (%p).  Read object will be new "
		       "singleton object", fgInstance, this));
       // delete fgInstance;
     }
     fgInstance = this;
  } else {
    R__b.WriteClassBuffer(AliCentralCorrectionManager::Class(),this);
  }
}
// @endcond
#endif

//____________________________________________________________________
//
// EOF
//
