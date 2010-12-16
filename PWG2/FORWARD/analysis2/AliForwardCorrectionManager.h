#ifndef ALIROOT_PWG2_FORWARD_ALIFORWARDCORRECTIONMANAGER_H
#define ALIROOT_PWG2_FORWARD_ALIFORWARDCORRECTIONMANAGER_H
#include <TObject.h>
#include "AliFMDCorrELossFit.h"
#include "AliFMDCorrSecondaryMap.h"
#include "AliFMDCorrDoubleHit.h"
#include "AliFMDCorrVertexBias.h"
#include "AliFMDCorrMergingEfficiency.h"
#include <TString.h>
class TFile;

/**
 * Manager (singleton) of corrections 
 * 
 * @ingroup pwg2_forward_corr 
 */
class AliForwardCorrectionManager : public TObject
{
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
    kAll                       = (kSecondaryMap| 
				  kELossFits|
				  kVertexBias|
				  kMergingEfficiency|
				  kDoubleHit)
  };
  /** 
   * Access to the singleton object 
   * 
   * @return Reference to the singleton object 
   */
  static AliForwardCorrectionManager& Instance();
  /** 
   * Read in corrections based on the parameters given 
   * 
   * @param collisionSystem Collision system
   * @param cmsNN           Center of mass energy per nuclean pair [GeV]
   * @param field           Magnetic field setting [kG]
   * @param mc              Monte-carlo switch
   * @param what            What to read in. 
   * @param force           Force (re-)reading of specified things
   * 
   * @return 
   */
  Bool_t Init(UShort_t collisionSystem, 
	      UShort_t cmsNN, 
	      Short_t  field, 
	      Bool_t   mc=false,
	      UInt_t   what=kAll,
	      Bool_t   force=false);
  Bool_t Init(const char* collisionSystem, 
	      Float_t     cmsNN, 
	      Float_t     field, 
	      Bool_t      mc=false,
	      UInt_t      what=kAll,
	      Bool_t      force=false);
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
  AliFMDCorrELossFit* GetELossFit() const { return fELossFit; }
  /** 
   * Get the secondary correction map
   * 
   * @return Get the secondary correction map object or null
   */
  AliFMDCorrSecondaryMap* GetSecondaryMap() const { return fSecondaryMap; }
  /** 
   * Get the double hit correction object
   * 
   * @return Get the double hit correction object or null 
   */
  AliFMDCorrDoubleHit* GetDoubleHit() const { return fDoubleHit; }
  /** 
   * Get the vertex bias correction object
   * 
   * @return Get the vertex bias correction object or null 
   */
  AliFMDCorrVertexBias* GetVertexBias() const { return fVertexBias; }
  AliFMDCorrMergingEfficiency* GetMergingEfficiency() const 
  {
    return fMergingEfficiency;
  }
  /** 
   * @{ 
   * @name Path, file, and object access utilities 
   */
  /** 
   * Get the path to the specified object 
   *
   * @param what  Which stuff to get the path for 
   * @param sys   Collision system
   * @param sNN   Center of mass energy [GeV]
   * @param field Magnetic field in the L3 magnet [kG]
   * @param mc    Whether the correction objects should be valid for MC
   * 
   * @return The full path or null 
   */
  TString GetFileName(ECorrection what, 
		      UShort_t    sys, 
		      UShort_t    sNN, 
		      Short_t     field,
		      Bool_t      mc) const;
  TString GetFileName(ECorrection what) const;
  /** 
   * Get the path to the specified object 
   *
   * @param what  Which stuff to get the path for 
   * 
   * @return The full path or null 
   */
  const Char_t* GetFileDir(ECorrection what) const;
  /** 
   * Get the path to the specified object 
   *
   * @param what  Which stuff to get the path for 
   * @param sys   Collision system
   * @param sNN   Center of mass energy [GeV]
   * @param field Magnetic field in the L3 magnet [kG]
   * @param mc    Whether the correction objects should be valid for MC
   * 
   * @return The full path or null 
   */
  TString GetFilePath(ECorrection what, 
		      UShort_t    sys, 
		      UShort_t    sNN, 
		      Short_t     field,
		      Bool_t      mc) const;
  TString GetFilePath(ECorrection what) const;
  /** 
   * Open the file that contains the correction object specified 
   * 
   * @param what  Which stuff to get the path for 
   * @param sys   Collision system
   * @param sNN   Center of mass energy [GeV]
   * @param field Magnetic field in the L3 magnet [kG]
   * @param mc    Whether the correction objects should be valid for MC
   * @param rw    Whether to open the file in read/write
   * @param newfile Wheter to make the file if it doesn't exist
   * 
   * @return The file that contains the correction object or null 
   */
  TFile* GetFile(ECorrection what, 
		 UShort_t    sys, 
		 UShort_t    sNN, 
		 Short_t     field,
		 Bool_t      mc=false, 
		 Bool_t      rw=false, 
		 Bool_t      newfile=false) const;
  TFile* GetFile(ECorrection what) const;
  /** 
   * Get the object name corresponding to correction type 
   * 
   * @param what Correction 
   * 
   * @return Object name or null
   */
  const Char_t* GetObjectName(ECorrection what) const;
  /** 
   * Check if the specified objet exists in the file, and return it
   * 
   * @param file File to query 
   * @param what Correction type 
   * 
   * @return Object found, or null
   */
  TObject* CheckObject(TFile* file,  ECorrection what) const;
  /** 
   * Get the path to the specified object 
   *
   * @param what  Which stuff to get the path for 
   * @param sys   Collision system
   * @param sNN   Center of mass energy [GeV]
   * @param field Magnetic field in the L3 magnet [kG]
   * @param mc    Whether the correction objects should be valid for MC
   * 
   * @return The full path or null 
   */
  TObject* GetObject(ECorrection what, 
		     UShort_t    sys, 
		     UShort_t    sNN, 
		     Short_t     field,
		     Bool_t      mc) const;
  TObject* GetObject(ECorrection what) const;
  /* 
   * @} 
   */
private:
  /** 
   * Default constructor 
   */
  AliForwardCorrectionManager();
  /** 
   * Copy constructor 
   * 
   * @param o Object to copy from 
   */
  AliForwardCorrectionManager(const AliForwardCorrectionManager& o);
  /** 
   * Assignment operator 
   * 
   * @param o Object to assign from 
   * 
   * @return Reference to this object 
   */
  AliForwardCorrectionManager& operator=(const AliForwardCorrectionManager& o);
  /** 
   * @{ 
   * @name Read in corrections 
   */
  /** 
   * Read in the secondary map 
   * 
   * @param sys   Collision system
   * @param sNN   Center of mass energy [GeV]
   * @param field Magnetic field in the L3 magnet [kG]
   * 
   * @return True on success, false otherwise 
   */
  Bool_t ReadSecondaryMap(UShort_t sys, UShort_t sNN, Short_t field);
  /** 
   * Read in the double hit correction
   * 
   * @param sys   Collision system
   * @param sNN   Center of mass energy [GeV]
   * @param field Magnetic field in the L3 magnet [kG]
   * 
   * @return True on success, false otherwise 
   */
  Bool_t ReadDoubleHit(UShort_t sys, UShort_t sNN, Short_t field);
  /** 
   * Read in the energy loss fits 
   * 
   * @param sys   Collision system
   * @param sNN   Center of mass energy [GeV]
   * @param field Magnetic field in the L3 magnet [kG]
   * @param mc    Whether the correction objects should be valid for MC
   * 
   * @return True on success, false otherwise 
   */
  Bool_t ReadELossFits(UShort_t sys, UShort_t sNN, Short_t field, Bool_t mc);
  /** 
   * Read in the event selection efficiency 
   * 
   * @param sys   Collision system
   * @param sNN   Center of mass energy [GeV]
   * @param field Magnetic field in the L3 magnet [kG]
   * 
   * @return True on success, false otherwise 
   */
  Bool_t ReadVertexBias(UShort_t sys, UShort_t sNN, Short_t field);
  /** 
   * Read in the merging efficiency 
   * 
   * @param sys   Collision system
   * @param sNN   Center of mass energy [GeV]
   * @param field Magnetic field in the L3 magnet [kG]
   * 
   * @return True on success, false otherwise 
   */
  Bool_t ReadMergingEfficiency(UShort_t sys, UShort_t sNN, Short_t field);
  /* 
   * @} 
   */
  
  /** Static singleton instance */
  static AliForwardCorrectionManager* fgInstance;
  Bool_t    fInit;  // whether we have been initialised 
  UShort_t  fSys;   // Collision System
  UShort_t  fSNN;   // Collision energy per nucleon (GeV)
  Short_t   fField; // L3 magnetic field (kG)
  
  /** 
   * @{
   * @name Paths 
   */ 
  TString fELossFitsPath;    // Path to energy loss fit correction 
  TString fMergingEffPath;   // Path to sharing efficiency correction 
  TString fSecondaryMapPath; // Path to secondary efficiency correction
  TString fDoubleHitPath;    // Path to double hit correction
  TString fVertexBiasPath;   // Path to event selection efficiency correction
  /* 
   * @}
   */
  /** 
   * @{ 
   * @name Object name 
   */
  static const Char_t* fgkSecondaryMapSkel;  // Name of correction object 
  static const Char_t* fgkDoubleHitSkel;     // Name of correction object 
  static const Char_t* fgkELossFitsSkel;     // Name of correction object 
  static const Char_t* fgkVertexBiasSkel;    // Name of correction object 
  static const Char_t* fgkMergingEffSkel;    // Name of correction object 
  /* 
   * @} 
   */
  /** 
   * @{ 
   * @name Correction objects 
   */
  AliFMDCorrELossFit*     fELossFit;     // Energy loss fits 
  AliFMDCorrSecondaryMap* fSecondaryMap; // Secondary particle correction
  AliFMDCorrDoubleHit*    fDoubleHit;    // Double hit correction (low flux)
  AliFMDCorrVertexBias*   fVertexBias;   // Vertex bias correction
  AliFMDCorrMergingEfficiency* fMergingEfficiency;
  /* 
   * @}
   */

  ClassDef(AliForwardCorrectionManager,1) // Manager of corrections 
};
//____________________________________________________________________
inline const TAxis* 
AliForwardCorrectionManager::GetEtaAxis() const
{
  if (!fSecondaryMap) return 0;
  return &(fSecondaryMap->GetEtaAxis());
}
//____________________________________________________________________
inline const TAxis* 
AliForwardCorrectionManager::GetVertexAxis() const
{
  if (!fSecondaryMap) return 0;
  return &(fSecondaryMap->GetVertexAxis());
}

#endif
// Local Variables:
//   mode: C++ 
// End: 

