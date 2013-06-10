//
// Manager (singleton) of corrections 
// 
#ifndef ALIFORWARDCORRECTIONMANAGEROADB_H
#define ALIFORWARDCORRECTIONMANAGEROADB_H
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
#include <TObject.h>
// #include "AliFMDCorrSecondaryMap.h"
#include <TString.h>
class TFile;
class TBrowser;
class TClass;
class AliFMDCorrELossFit;
class AliFMDCorrDoubleHit;
class AliFMDCorrVertexBias;
class AliFMDCorrMergingEfficiency;
class AliFMDCorrAcceptance;
class AliFMDCorrSecondaryMap;
class AliOADBForward;
class TAxis;

/**
 * Manager (singleton) of corrections 
 *
 * Note, that this class has a custom streamer.  That is to ensure
 * that the singleton pointer is correctly set on reading in an object
 * of this type.
 * 
 * @ingroup pwglf_forward_corr 
 * @ingroup pwglf_forward_aod
 */
class AliForwardCorrectionManagerOADB : public TObject
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
  AliForwardCorrectionManagerOADB();
  /** 
   * Access to the singleton object 
   * 
   * @return Reference to the singleton object 
   */
  static AliForwardCorrectionManagerOADB& Instance();
  /** 
   *
   * @param prefix Prefix to correction objects. 
   */
  void SetPrefix(const char* prefix);
  /** 
   * Set the file directory for a type 
   * 
   * @param what     Type 
   * @param filename Name of file that contains tree of corrections
   */
  void SetFile(ECorrection what, const char* filename);
  /** 
   * Set path to corrections 
   * 
   * @param d Path
   */
  void SetSecondaryMapPath(const char* d) { SetFile(kSecondaryMap, d); }
  /** 
   * Set path to corrections 
   * 
   * @param d Path
   */
  void SetDoubleHitPath(const char* d)    { SetFile(kDoubleHit, d); }
  /** 
   * Set path to corrections 
   * 
   * @param d Path
   */
  void SetELossFitsPath(const char* d)    { SetFile(kELossFits, d); }
  /** 
   * Set path to corrections 
   * 
   * @param d Path
   */
  void SetVertexBiasPath(const char* d)   { SetFile(kVertexBias, d); }
  /** 
   * Set path to corrections 
   * 
   * @param d Path
   */
  void SetMergingEffPath(const char* d)   { SetFile(kMergingEfficiency, d); }
  /** 
   * Set path to corrections 
   * 
   * @param d Path
   */
  void SetAcceptancePath(const char* d)   { SetFile(kAcceptance, d); }
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
  Bool_t Init(ULong_t  runNumber,
	      UShort_t collisionSystem, 
	      UShort_t cmsNN, 
	      Short_t  field, 
	      Bool_t   mc,
	      Bool_t   satelliteCollisions,
	      UInt_t   what,
	      Bool_t   force);
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
  Bool_t Init(ULong_t     runNumber, 
	      const char* collisionSystem, 
	      Float_t     cmsNN, 
	      Float_t     field, 
	      Bool_t      mc,
	      Bool_t      satelliteCollisions,
	      UInt_t      what,
	      Bool_t      force);
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
  /** 
   * Write a correction object to (a temporary) file.  
   *    
   * @param what   What kind of correction
   * @param sys    Collision system
   * @param sNN    Center of mass energy
   * @param field  Field 
   * @param mc     Whether this is for MC only
   * @param o      Object to write
   * @param full   If true, write to full path
   * 
   * @return True on success 
   */
  Bool_t StoreObject(ULong_t     runNo,
		     UShort_t    sys, 
		     UShort_t    sNN, 
		     Short_t     field, 
		     Bool_t      mc,
		     Bool_t      sat, 
		     TObject*    o, 
		     Bool_t      full,
		     const char* meth="NEAR") const;
  /** 
   * Write a correction object to (a temporary) file.  
   *    
   * @param what   What kind of correction
   * @param sys    Collision system
   * @param sNN    Center of mass energy
   * @param field  Field 
   * @param mc     Whether this is for MC only
   * @param o      Object to write
   * @param full   If true, write to full path
   * 
   * @return True on success 
   */
  Bool_t StoreObject(ECorrection what, 
		     ULong_t     runNo,
		     UShort_t    sys, 
		     UShort_t    sNN, 
		     Short_t     field, 
		     Bool_t      mc,
		     Bool_t      sat, 
		     TObject*    o, 
		     Bool_t      full,
		     const char* meth="NEAR") const;
  /** 
   * Write a correction object to (a temporary) file.  
   *    
   * @param what   What kind of correction
   * @param sys    Collision system
   * @param sNN    Center of mass energy
   * @param field  Field 
   * @param mc     Whether this is for MC only
   * @param o      Object to write
   * @param full   If true, write to full path
   * 
   * @return True on success 
   */
  Bool_t StoreObject(const TString& what,
		     ULong_t     runNo,
		     UShort_t    sys, 
		     UShort_t    sNN, 
		     Short_t     field, 
		     Bool_t      mc,
		     Bool_t      sat, 
		     TObject*    o, 
		     Bool_t      full,
		     const char* meth="NEAR") const;
  /** 
   * Get the object name corresponding to correction type 
   * 
   * @param what Correction 
   * 
   * @return Object name or null
   */
  const Char_t* GetObjectName(ECorrection what) const;
  /** 
   * Get the class associated with a correction
   * 
   * @param what Correction type 
   * 
   * @return Pointer to TClass object or null
   */
  const TClass* GetObjectClass(ECorrection what) const;
  /** 
   * Get the correction type from the table name
   * 
   * @param what Table name
   * 
   * @return Correction type or kAll
   */
  ECorrection   GetObjectType(const TString& what) const;
  /** 
   * Get the correction type from an object
   * 
   * @param what Object
   * 
   * @return Correction type or kAll
   */
  ECorrection   GetObjectType(const TObject* obj) const;
  /** 
   * Get the path to the specified object 
   *
   * @param what  Which stuff to get the path for 
   * 
   * @return The full path or null 
   */
  const TString& GetFilePath(ECorrection what) const;
private:
  /** 
   * Copy constructor 
   * 
   * @param o Object to copy from 
   */
  AliForwardCorrectionManagerOADB(const AliForwardCorrectionManagerOADB& o);
  /** 
   * Assignment operator 
   * 
   * @param o Object to assign from 
   * 
   * @return Reference to this object 
   */
  AliForwardCorrectionManagerOADB& operator=(const AliForwardCorrectionManagerOADB& o);
  /** 
   * @{ 
   * @name Read in corrections 
   */
  /** 
   * Read in the secondary map 
   * 
   * @param runNo Run number 
   * @param sys   Collision system
   * @param sNN   Center of mass energy [GeV]
   * @param field Magnetic field in the L3 magnet [kG]
   * @param sat   If true, get satellite collision corrections 
   * 
   * @return True on success, false otherwise 
   */
  Bool_t ReadSecondaryMap(ULong_t  runNo, 
			  UShort_t sys,
			  UShort_t sNN, 
			  Short_t  field,
			  Bool_t   sat);
  /** 
   * Read in the double hit correction
   * 
   * @param sys   Collision system
   * @param sNN   Center of mass energy [GeV]
   * @param field Magnetic field in the L3 magnet [kG]
   * 
   * @return True on success, false otherwise 
   */
  Bool_t ReadDoubleHit(ULong_t  runNo, 
		       UShort_t sys, 
		       UShort_t sNN, 
		       Short_t  field,
		       Bool_t   mc);
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
  Bool_t ReadELossFits(ULong_t  runNo, 
		       UShort_t sys, 
		       UShort_t sNN, 
		       Short_t  field, 
		       Bool_t   mc, 
		       Bool_t   sat);
  /** 
   * Read in the event selection efficiency 
   * 
   * @param sys   Collision system
   * @param sNN   Center of mass energy [GeV]
   * @param field Magnetic field in the L3 magnet [kG]
   * 
   * @return True on success, false otherwise 
   */
  Bool_t ReadVertexBias(ULong_t  runNo, 
			UShort_t sys, 
			UShort_t sNN, 
			Short_t  field, 
			Bool_t   sat);
  /** 
   * Read in the merging efficiency 
   * 
   * @param sys   Collision system
   * @param sNN   Center of mass energy [GeV]
   * @param field Magnetic field in the L3 magnet [kG]
   * 
   * @return True on success, false otherwise 
   */
  Bool_t ReadMergingEfficiency(ULong_t  runNo, 
			       UShort_t sys, 
			       UShort_t sNN, 
			       Short_t  field,
			       Bool_t   sat);
  /** 
   * Read in the acceptance correction due to dead-channels 
   * 
   * @param sys   Collision system		      
   * @param sNN   Center of mass energy [GeV]	      
   * 
   * @return True on success, false otherwise 
   */
  Bool_t ReadAcceptance(ULong_t  runNo, 
			UShort_t sys, 
			UShort_t sNN, 
			Bool_t   sat);
  /* 
   * @} 
   */
  /** 
   * @{ 
   * @name Path, file, and object access utilities 
   */
  /** 
   * Get the file name of the specified object
   * 
   * @param what Which stuff to get the path for 
   * 
   * @return The file name (sans directory) or null
   */
  const Char_t* GetFileName(ECorrection what) const;
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
  TObject* GetObject(ECorrection what, 
		     ULong_t     runNo,
		     UShort_t    sys, 
		     UShort_t    sNN, 
		     Short_t     field,
		     Bool_t      mc,
		     Bool_t      sat) const;
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
  
  /** Static singleton instance */
  static AliForwardCorrectionManagerOADB* fgInstance; // Skeleton

  /** 
   * @{
   * @name Cache variables 
   */
  Bool_t    fInit;  // whether we have been initialised 
  ULong_t   fRunNo; // Run number to use 
  UShort_t  fSys;   // Collision System
  UShort_t  fSNN;   // Collision energy per nucleon (GeV)
  Short_t   fField; // L3 magnetic field (kG)
  Bool_t    fMC;    // Are we doing MC?
  Bool_t    fSat;   // Are we doing satellite collisions 
  /* @} */

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

  /** 
   * @{ 
   * @name Database interface 
   */
  AliOADBForward* fDB; //! DB interface - do not store!
  /** @} */
  ClassDef(AliForwardCorrectionManagerOADB,1) // Manager of corrections 
};

#endif
// Local Variables:
//   mode: C++ 
// End: 

