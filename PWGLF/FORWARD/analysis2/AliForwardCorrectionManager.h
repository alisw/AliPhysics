//
// Manager (singleton) of corrections 
// 
#ifndef ALIFORWARDCORRECTIONMANAGER_H
#define ALIFORWARDCORRECTIONMANAGER_H
/**
 * @file   AliForwardCorrectionManager.h
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
class AliFMDCorrELossFit;
class AliFMDCorrDoubleHit;
class AliFMDCorrVertexBias;
class AliFMDCorrMergingEfficiency;
class AliFMDCorrAcceptance;
class AliFMDCorrSecondaryMap;
class AliFMDCorrNoiseGain;
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
class AliForwardCorrectionManager : public AliCorrectionManagerBase
{
private:
  /**
   * Enumeration of things that can be read in 
   */
  enum EId { 
    kIdSecondaryMap            = 0, 
    kIdELossFits,
    kIdVertexBias,
    kIdMergingEfficiency,
    kIdDoubleHit,
    kIdAcceptance,
    kIdNoiseGain
  };
public:
  /**
   * Enumeration of things that can be read in 
   */
  enum ECorrection { 
    kSecondaryMap              = 0x01, 
    kELossFits                 = 0x02, 
    kVertexBias                = 0x04, 
    kMergingEfficiency         = 0x08,
    kDoubleHit                 = 0x10,
    kAcceptance                = 0x20,
    kNoiseGain                 = 0x40,
    kDefault                   = (kSecondaryMap|kELossFits|kAcceptance),
    kAll                       = (kSecondaryMap| 
				  kELossFits|
				  kVertexBias|
				  kMergingEfficiency|
				  kDoubleHit|
				  kAcceptance|
				  kNoiseGain)
  };
  /** 
   * Default constructor.  This is public for the sake of the ROOT I/O
   * system, but should never be used outside of that system - that
   * is, do not use this constructor
   */
  AliForwardCorrectionManager();
  /** 
   * Access to the singleton object 
   * 
   * @return Reference to the singleton object 
   */
  static AliForwardCorrectionManager& Instance();
  /** 
   * @return name of object 
   */
  const Char_t* GetName() const { return "forwardCorrections"; }
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
  void SetDoubleHitPath(const char* d)    
  {
    SetCorrectionFile(kIdDoubleHit, d);
  }
  /** 
   * Set path to corrections 
   * 
   * @param d Path
   */
  void SetELossFitsPath(const char* d)    
  {
    SetCorrectionFile(kIdELossFits, d);
  }
  /** 
   * Set path to corrections 
   * 
   * @param d Path
   */
  void SetVertexBiasPath(const char* d)   
  {
    SetCorrectionFile(kIdVertexBias, d);
  }
  /** 
   * Set path to corrections 
   * 
   * @param d Path
   */
  void SetMergingEffPath(const char* d)   
  {
    SetCorrectionFile(kIdMergingEfficiency, d);
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
	      UInt_t      what=kDefault,
	      Bool_t      force=false);
  /** 
   * Parse string with fields in it, and return the corresponding bit mask
   * 
   * @param what The string 
   * 
   * @return The corresponding bit mask
   */
  static UInt_t ParseFields(const TString& what);
  /** 
   * Get the eta axis 
   * 
   * @return Eta axis or null
   */
  const TAxis* GetEtaAxis() const;
  /** 
   * Get the vertex axis 
   * 
   * @return The vertex axis or null
   */
  const TAxis* GetVertexAxis() const;
  /** 
   * Get the energy loss fit correction object. 
   * 
   * @return Get the energy loss fits corrections object or null pointer
   */
  const AliFMDCorrELossFit* GetELossFit() const;
  /** 
   * Alias for GetELossFit
   * 
   * @return Get the energy loss fits corrections object or null pointer
   */
  const AliFMDCorrELossFit* GetELossFits() const { return GetELossFit(); }
  /** 
   * Get the secondary correction map
   * 
   * @return Get the secondary correction map object or null
   */
  const AliFMDCorrSecondaryMap* GetSecondaryMap() const;
  /** 
   * Get the double hit correction object
   * 
   * @return Get the double hit correction object or null 
   */
  const AliFMDCorrDoubleHit* GetDoubleHit() const;
  /** 
   * Get the vertex bias correction object
   * 
   * @return Get the vertex bias correction object or null 
   */
  const AliFMDCorrVertexBias* GetVertexBias() const;
  /** 
   * Get the merging efficiency 
   * 
   * @return Get the vertex efficiency correction 
   */
  const AliFMDCorrMergingEfficiency* GetMergingEfficiency() const;
  /** 
   * Get the acceptance correction due to dead channels 
   * 
   * @return Acceptance correction due to dead channels 
   */
  const AliFMDCorrAcceptance* GetAcceptance() const;
  /** 
   * Get the noise calibration.  That is, the ratio 
   *
   * @f[ 
   *   \frac{\sigma_{i}}{g_{i}k}
   * @f] 
   *
   * where @f$ k@f$ is a constant determined by the electronics of
   * units DAC/MIP, and @f$ \sigma_i, g_i@f$ are the noise and gain of
   * the @f$ i @f$ strip respectively
   * 
   * @return Noise-gain calibration
   */
  const AliFMDCorrNoiseGain* GetNoiseGain() const;
private:
  /** 
   * Non-default constructor - initializes corrections - used by
   * singleton access member function Instance
   * 
   * @param notUsed Ignored
   */
  AliForwardCorrectionManager(Bool_t notUsed);
  
  /** Static singleton instance */
  static AliForwardCorrectionManager* fgInstance; // Skeleton

  /** 
   * @{ 
   * @name Object name 
   */
  static const Char_t* fgkSecondaryMapSkel;  // Name of correction object 
  static const Char_t* fgkDoubleHitSkel;     // Name of correction object 
  static const Char_t* fgkELossFitsSkel;     // Name of correction object 
  static const Char_t* fgkVertexBiasSkel;    // Name of correction object 
  static const Char_t* fgkMergingEffSkel;    // Name of correction object 
  static const Char_t* fgkAcceptanceSkel;    // Name of correction object 
  static const Char_t* fgkNoiseGainSkel;     // Name of correction object 
  /* 
   * @} 
   */
  ClassDef(AliForwardCorrectionManager,5) // Manager of corrections 
};

#endif
// Local Variables:
//   mode: C++ 
// End: 

