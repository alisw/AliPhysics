//
// Manager (singleton) of corrections 
// 
#ifndef ALICENTRALCORRECTIONMANAGER_H
#define ALICENTRALCORRECTIONMANAGER_H
/**
 * @file   AliCentralCorrectionManager.h
 * @author Christian Holm Christensen <cholm@dalsgaard.hehi.nbi.dk>
 * @date   Wed Mar 23 14:04:27 2011
 * 
 * @brief  
 * 
 * 
 * @ingroup pwglf_forward_aod
 */
#include "AliCorrectionManagerBase.h"
#include <TString.h>
class TBrowser;
class AliCentralCorrAcceptance;
class AliCentralCorrSecondaryMap;
class TAxis;

/**
 * Manager (singleton) of corrections 
 *
 * Note, that this class has a custom streamer.  That is to ensure
 * that the singleton pointer is correctly set on reading in an object
 * of this type.
 * 
 * @ingroup pwglf_forward_corr 
 */
class AliCentralCorrectionManager : public AliCorrectionManagerBase
{
private:
  /**
   * Enumeration of things that can be read in 
   */
  enum EId { 
    kIdSecondaryMap            = 0, 
    kIdAcceptance
  };
public:
  /**
   * Enumeration of things that can be read in 
   */
  enum ECorrection { 
    kSecondaryMap              = 0x01, 
    kAcceptance                = 0x02,
    kDefault                   = (kSecondaryMap|kAcceptance),
    kAll                       = (kSecondaryMap|kAcceptance) 
  };
  /** 
   * Default constructor.  This is public for the sake of the ROOT I/O
   * system, but should never be used outside of that system - that
   * is, do not use this constructor
   */
  AliCentralCorrectionManager();
  /** 
   * Access to the singleton object 
   * 
   * @return Reference to the singleton object 
   */
  static AliCentralCorrectionManager& Instance();
  /** 
   * @return name of the object 
   */
  const Char_t* GetName() const { return "centralCorrections"; }
  /** 
   * Set path to corrections 
   * 
   * @param d Path
   */
  void SetSecondaryMapPath(const char* d) 
  {
    SetCorrectionFile(kIdSecondaryMap, d);
  }
  /** 
   * Set path to corrections 
   * 
   * @param d Path
   */
  void SetAcceptancePath(const char* d)   
  {
    SetCorrectionFile(kIdAcceptance, d);
  }
  /** 
   * Read in corrections based on the parameters given 
   * 
   * @param runNumber       Run number
   * @param collisionSystem Collision system
   * @param cmsNN           Center of mass energy per nuclean pair [GeV]
   * @param field           Magnetic field setting [kG]
   * @param mc              Monte-carlo switch
   * @param what            What to read in. 
   * @param force           Force (re-)reading of specified things
   * @param satelliteCollisions For satellite collisions
   * 
   * @return 
   */
  Bool_t Init(ULong_t  runNumber,
	      UShort_t collisionSystem, 
	      UShort_t cmsNN, 
	      Short_t  field, 
	      Bool_t   mc=false,
	      Bool_t   satelliteCollisions=false,
	      UInt_t   what=kDefault,
	      Bool_t   force=false);
  /** 
   * Read in correction based on passed parameters
   * 
   * @param runNumber       Run number
   * @param collisionSystem Collision system string 
   * @param cmsNN           Center of mass energy per nucleon pair [GeV]
   * @param field           Magnetic field [kG]
   * @param mc              Monte-carlo switch
   * @param what            What to read in 
   * @param force           Force (re-)reading of specified things
   * @param satelliteCollisions For satellite collisions
   * 
   * @return true on success
   */
  Bool_t Init(ULong_t     runNumber, 
	      const char* collisionSystem, 
	      Float_t     cmsNN, 
	      Float_t     field, 
	      Bool_t      mc=false,
	      Bool_t      satelliteCollisions=false,
	      UInt_t      what=kStandard,
	      Bool_t      force=false);
  /** 
   * Get the vertex axis 
   * 
   * @return The vertex axis or null
   */
  const TAxis* GetVertexAxis() const;
  /** 
   * Get the @f$\eta@f$ axis
   * 
   * @return The @f$\eta@f$ axis or null
   */
  const TAxis* GetEtaAxis() const { return 0; }
  /** 
   * Get the secondary correction map
   * 
   * @return Get the secondary correction map object or null
   */
  const AliCentralCorrSecondaryMap* GetSecondaryMap() const;
  /** 
   * Get the acceptance correction due to dead channels 
   * 
   * 
   * @return Acceptance correction due to dead channels 
   */
  const AliCentralCorrAcceptance* GetAcceptance() const;
private:
  /** 
   * Non-default constructor - initializes corrections - used by
   * singleton access member function Instance
   * 
   * @param notUsed Ignored
   */
  AliCentralCorrectionManager(Bool_t notUsed);
  
  /** Static singleton instance */
  static AliCentralCorrectionManager* fgInstance; // Skeleton

  /** 
   * @{ 
   * @name Object name 
   */
  static const Char_t* fgkSecondaryMapSkel;  // Name of correction object 
  static const Char_t* fgkAcceptanceSkel;    // Name of correction object 
  /* 
   * @} 
   */
  ClassDef(AliCentralCorrectionManager,2) // Manager of corrections 
};

#endif
// Local Variables:
//   mode: C++ 
// End: 

