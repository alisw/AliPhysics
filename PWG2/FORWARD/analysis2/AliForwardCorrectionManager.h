//
// Manager (singleton) of corrections 
// 
#ifndef ALIFORWARDCORRECTIONMANAGER_H
#define ALIFORWARDCORRECTIONMANAGER_H
#include <TObject.h>
// #include "AliFMDCorrELossFit.h"
#include "AliFMDCorrSecondaryMap.h"
// #include "AliFMDCorrDoubleHit.h"
// #include "AliFMDCorrVertexBias.h"
// #include "AliFMDCorrMergingEfficiency.h"
// #include "AliFMDCorrAcceptance.h"
#include <TString.h>
class TFile;
class TBrowser;
class AliFMDCorrELossFit;
// class AliFMDCorrSecondaryMap;
class AliFMDCorrDoubleHit;
class AliFMDCorrVertexBias;
class AliFMDCorrMergingEfficiency;
class AliFMDCorrAcceptance;

/**
 * Manager (singleton) of corrections 
 *
 * Note, that this class has a custom streamer.  That is to ensure
 * that the singleton pointer is correctly set on reading in an object
 * of this type.
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
    kAcceptance                = 0x20,
    kAll                       = (kSecondaryMap| 
				  kELossFits|
				  kVertexBias|
				  kMergingEfficiency|
				  kDoubleHit|
				  kAcceptance)
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
   *
   * @param prefix Prefix to correction objects. 
   */
  void SetPrefix(const char* prefix);
  /** 
   * Set the file directory for a type 
   * 
   * @param what     Type 
   * @param dirname  Directory name 
   */
  void SetFileDir(ECorrection what, const char* dirname);
  void SetSecondaryMapPath(const char* d) { SetFileDir(kSecondaryMap, d); }
  void SetDoubleHitPath(const char* d)    { SetFileDir(kDoubleHit, d); }
  void SetELossFitsPath(const char* d)    { SetFileDir(kELossFits, d); }
  void SetVertexBiasPath(const char* d)   { SetFileDir(kVertexBias, d); }
  void SetMergingEffPath(const char* d)   { SetFileDir(kMergingEfficiency, d); }
  void SetAcceptancePath(const char* d)   { SetFileDir(kAcceptance, d); }
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
  /** 
   * Read in correction based on passed parameters
   * 
   * @param collisionSystem Collision system string 
   * @param cmsNN           Center of mass energy per nucleon pair [GeV]
   * @param field           Magnetic field [kG]
   * @param mc              Monte-carlo switch
   * @param what            What to read in 
   * @param force           Force (re-)reading of specified things
   * 
   * @return true on success
   */
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
  /** 
   * Get the merging efficiency 
   * 
   * 
   * @return Get the vertex efficiency correction 
   */
  AliFMDCorrMergingEfficiency* GetMergingEfficiency() const 
  {
    return fMergingEfficiency;
  }
  /** 
   * Get the acceptance correction due to dead channels 
   * 
   * 
   * @return Acceptance correction due to dead channels 
   */
  AliFMDCorrAcceptance* GetAcceptance() const { return fAcceptance; }
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
   * @return The file name (sans directory) or null 
   */
  TString GetFileName(ECorrection what, 
		      UShort_t    sys, 
		      UShort_t    sNN, 
		      Short_t     field,
		      Bool_t      mc) const;
  /** 
   * Get the file name of the specified object
   * 
   * @param what Which stuff to get the path for 
   * 
   * @return The file name (sans directory) or null
   */
  TString GetFileName(ECorrection what) const;
  /** 
   * Get the path to the specified object 
   *
   * @param what  Which stuff to get the path for 
   * 
   * @return The files directory full path or null 
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
  /** 
   * Get the full path to the object.  Note, the manager must be
   * initialised for this to work
   * 
   * @param what Which stuff to get the path for 
   * 
   * @return The full path or null
   */
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
  /** 
   * Get the file that contains the object specifed.  Note, the manager
   * must be initialised for this to work. 
   * 
   * @param what Which stuff to get the path for 
   * 
   * @return The file that contains the correction object or null
   */
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
  /** 
   * Get the object that contaisn the specified correction
   * 
   * @param what Which object to get
   * 
   * @return The object or null
   */
  TObject* GetObject(ECorrection what) const;
  /* 
   * @} 
   */
  /**
   * @{ 
   * @name Misc 
   */
  /** 
   * Print this object 
   * 
   * @param option Passed verbatim to correction objects
   */
  void Print(Option_t* option="") const;
  /** 
   * Browse this object 
   * 
   * @param b Browser to use 
   */
  void Browse(TBrowser* b);
  /* 
   * @}
   */
private:
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
  /** 
   * Read in the acceptance correction due to dead-channels 
   * 
   * @param sys   Collision system		      
   * @param sNN   Center of mass energy [GeV]	      
   * @param field Magnetic field in the L3 magnet [kG]
   * 
   * @return True on success, false otherwise 
   */
  Bool_t ReadAcceptance(UShort_t sys, UShort_t sNN, Short_t field);
  /* 
   * @} 
   */
  
  /** Static singleton instance */
  static AliForwardCorrectionManager* fgInstance; // Skeleton
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
  TString fAcceptancePath;   // Path to acceptance correction from dead areas
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
  static const Char_t* fgkAcceptanceSkel;    // Name of correction object 
  /* 
   * @} 
   */
  /** 
   * @{ 
   * @name Correction objects 
   */
  AliFMDCorrELossFit*          fELossFit;     // Energy loss fits 
  AliFMDCorrSecondaryMap*      fSecondaryMap; // Secondary particle correction
  AliFMDCorrDoubleHit*         fDoubleHit;    // Double hit corr. (low flux)
  AliFMDCorrVertexBias*        fVertexBias;   // Vertex bias correction
  AliFMDCorrMergingEfficiency* fMergingEfficiency; // Merging eff. 
  AliFMDCorrAcceptance*        fAcceptance;   // Acceptance corr. 
  /* 
   * @}
   */
  ClassDef(AliForwardCorrectionManager,2) // Manager of corrections 
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

