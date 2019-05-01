// -*- mode: C++ -*-
/**
 * @file   AliCorrectionManagerBase.h
 * @author Christian Holm Christensen <cholm@master.hehi.nbi.dk>
 * @date   Sun May 19 21:56:13 2013
 * 
 * @brief  Base class for correction managers
 * 
 * 
 */
#ifndef ALICORRECTIONMANAGERBASE_H
#define ALICORRECTIONMANAGERBASE_H
#include <TString.h>
#include <TNamed.h>
#include <TObjArray.h>
class AliOADBForward;
class TBrowser;
class TAxis;

/**
 * Base class for correction managers. 
 *
 * A correction is added to the manager by calling RegisterCorrection 
 *
 * @code 
 * class MyManager : public AliCorrectionManager
 * {
 * public: 
 *   enum {
 *     kA, 
 *     kB 
 *   };
 *   MyManager() 
 *     : AliCorrectionManager(Bool_t setup=false) 
 *   {
 *     if (setup) {
 *       RegisterCorrection(kA, "A", "/some/path/to/file", 
 *                          TH2D::Class(), kStandard);
 *       RegisterCorrection(kB, "B", "/some/path/to/file", 
 *                          TParameter<float>::Class(), kStandard);
 *     }
 *   }
 *   void Init(Bool_t useA, Bool_t useB, ULong_t run, UShort_t sys,
 *             UShort_t sNN, Short_t fld, Bool_t mc, Bool_t sat, 
 *             Bool_t force=false)
 *   {
 *     if (useA) GetCorrection(kA)->fEnabled = true;
 *     if (useB) GetCorrection(kB)->fEnabled = true;
 * 
 *     return InitCorrections(run, sys, sNN, fld, mc, sat, force);
 *   }
 *   TH2D* GetA() const { return static_cast<TH2D*>(Get(kA)); }
 *   TParameter<float>* GetB() const { return static_cast<TParameter<float>*>(Get(kB)); }
 *  };
 * @endcode
 *
 * In case the derived object is to be a singleton, one has to take a
 * little extra care about when the constructor is called - especially
 * if the singleton is to be streamed to disk:
 *
 * @code 
 * class MyManager : public AliCorrectionManager
 * {
 * public: 
 *   // As above, except the constructor must be private and 
 *   MyManager& Instance() 
 *   {
 *     if (!fgInstance) fgInstance = MyManager(true);
 *     return fgInstance; 
 *   }
 *   static MyManager* fgInstance; 
 * };
 * @endcode 
 * 
 * It is important - for I/O that default construction does not
 * register the corrections.  This should only be done in-code on
 * first construction.
 *
 */
class AliCorrectionManagerBase : public TObject
{
public:
  enum EConstants { 
    kIgnoreValue = 0,
    kIgnoreField = 999 // Must be synced with AliOADBForward::kInvalidField
  };
  enum EFields {
    kRun       = 0x01, 
    kSys       = 0x02, 
    kSNN       = 0x04, 
    kField     = 0x08, 
    kMC        = 0x10, 
    kSatellite = 0x20,
    kStandard  = kRun|kSys|kSNN|kField,
    kFull      = kStandard|kMC|kSatellite
  };
  /** 
   * Destructor 
   */
  virtual ~AliCorrectionManagerBase();
  /** 
   * Set the prefix to use when looking for the input files 
   * 
   * @param prefix Prefix to use for all corrections 
   */
  virtual void SetPrefix(const TString& prefix);
  /** 
   * Set whether to enable fall-back queries 
   * 
   * @param use If true, enable fall-back queries 
   */
  virtual void SetEnableFallBack(Bool_t use=true) { fFallBack = use; }

  /** 
   * @{ 
   * @name Storing corrections 
   */
  /** 
   * Store a correction 
   * 
   * @param o      Object to store
   * @param runNo  Run number of run this was created from
   * @param sys    Collision system (1:pp, 2:PbPb, 3:pPb)
   * @param sNN    Center of mass energy in GeV
   * @param field  L3 magnetic field in kG
   * @param mc     If true, this is for simulations
   * @param sat    If true, contain info for satellite interactions
   * @param file   (optional) file name to store in 
   * @param meth   (optional) method for run look to use. 
   * 
   * @return true on success
   */
  virtual Bool_t Store(TObject*     o,
		       ULong_t     runNo,
		       UShort_t    sys, 
		       UShort_t    sNN, 
		       Short_t     field, 
		       Bool_t      mc,
		       Bool_t      sat, 
		       const char* file,
		       const char* meth="NEAR") const;
  /** 
   * Append the content of the file @a addition to the @a destination
   * file for this manager.  This used TFileMerger::PartialMerge 
   * 
   * @param destination Filename of destination storage (in OADB_PATH)
   * @param addition    Filename of addition. 
   * 
   * @return true on success 
   */
  virtual Bool_t Append(const TString& addition,
			const TString& destination="") const;
  /** 
   * Write a new database file with tables that only has one entry for
   * each query tuple.
   * 
   * @param destination Where to write the new file 
   * @param verb        Whether to be verbose 
   * @param all         If true, copy all entries
   * 
   * @return true on success, false otherwise 
   */
  virtual Bool_t CleanUp(const TString& destination,
			 Bool_t         verb=false,
			 Bool_t         all=false) const;
  /* @} */

  /** 
   * @{ 
   * @name Getting the corrections 
   */
  /** 
   * Convinience function to enable corrections on-mass.  User class
   * should overload this to properly enable corrections based on the
   * bit identifiers.
   *
   * @param what Bit mask of correction identifiers. 
   */
  virtual void EnableCorrections(UInt_t what);
  /** 
   * Check if the specified corrrections have been initialized
   * 
   * @param what     Corrections to check 
   * @param verbose  If true, be verbose 
   * 
   * @return true if all specified the corrections have been 
   */
  virtual Bool_t CheckCorrections(UInt_t what, Bool_t verbose=true) const;
  /** 
   * Read in all corrections 
   * 
   * @param run    Run number 
   * @param sys    System 
   * @param sNN    Center of mass energy 
   * @param fld    L3 magnetic field
   * @param mc     For simulations 
   * @param sat    For satellite interactions 
   * @param force  Force-reread 
   * 
   * @return true on success 
   */
  Bool_t InitCorrections(ULong_t    run, 
			 UShort_t   sys, 
			 UShort_t   sNN, 
			 Short_t    fld, 
			 Bool_t     mc, 
			 Bool_t     sat,
			 Bool_t     force=false);
  /* @} */

  /** 
   * @{ 
   * @name Get axis objects. 
   */
  /** 
   * Get the vertex axis 
   * 
   * @return The vertex axis or null
   */
  virtual const TAxis* GetVertexAxis() const { return 0; }
  /** 
   * Get the @f$\eta@f$ axis
   * 
   * @return The @f$\eta@f$ axis or null
   */
  virtual const TAxis* GetEtaAxis() const { return 0; }
  /* @} */

  /** 
   * @{ 
   * @name Misc
   */
  /** 
   * Check if the manager is initialized 
   * 
   * @return True if initialized
   */
  virtual Bool_t IsInit() const { return fIsInit; }
  /** 
   * Print information
   * 
   * @param option Options:
   * 
   *   - R  Recursive list each correction 
   *   - D  Also give details for each correction
   */
  virtual void Print(Option_t* option="") const;
  /** 
   * Browse this object
   * 
   * @param b Browser to use 
   */
  virtual void Browse(TBrowser* b);
  /** 
   * Flag that this is a folder 
   * 
   * @return Always true
   */
  virtual Bool_t IsFolder() const { return true; }
  /** 
   * Set whehter to enable debug information 
   * 
   * @param debug if true, do verbose queries 
   */
  virtual void SetDebug(Bool_t debug) { fDebug = debug; }
  /* @} */
  /** 
   * Correction registration.  These objects are streamed, so we need
   * to make them public.  That does _not_ mean that these objects
   * should be used outside of this (or derived) classe.
   */
  struct Correction : public TNamed
  {
    /** 
     * Constructor - for I/O
     */
    Correction();
    /** 
     * Constructor 
     * 
     * @param tableName   Table name
     * @param fileName    File name 
     * @param cls         Class
     * @param queryFields Enabled fields
     * @param enabled     Enabled or not 
     */
    Correction(const TString& tableName, 
	       const TString& fileName, 
	       TClass*        cls,
	       UShort_t       queryFields=kStandard,
	       Bool_t         enabled=false);
    /** 
     * Copy constructor
     * 
     * @param o Object to copy from 
     */
    Correction(const Correction& o);
    /** 
     * Assignement operator
     * 
     * @param o Object to assign from 
     * 
     * @return reference to this obejct
     */    
    Correction& operator=(const Correction& o);
    /** 
     * Destructor
     */
    ~Correction() { delete fObject; }
    /** 
     * Read the correction
     * 
     * @param db   Database interface 
     * @param run  Run number
     * @param sys  Collision system
     * @param sNN  Center of mass energy per nucleon
     * @param fld  L3 magnetic field
     * @param mc   If true, for simulated data, else real
     * @param sat  If true, for satellite interactions
     * @param vrb  If true, do verbose query
     * @param fbk  Force read
     *
     * @return true on success
     */
    Bool_t ReadIt(AliOADBForward* db,
		  ULong_t    run, 
		  UShort_t   sys, 
		  UShort_t   sNN, 
		  Short_t    fld, 
		  Bool_t     mc, 
		  Bool_t     sat,
		  Bool_t     vrb=false,
		  Bool_t     fbk=false);
    /** 
     * Store a correction
     * 
     * @param db    Possible database interface
     * @param o     Object to store 
     * @param run  Run number
     * @param sys  Collision system
     * @param sNN  Center of mass energy per nucleon
     * @param fld  L3 magnetic field
     * @param mc   If true, for simulated data, else real
     * @param sat  If true, for satellite interactions
     * @param file File to store in
     * @param meth Default run method
     * 
     * @return true on success
     */
    Bool_t StoreIt(AliOADBForward* db,
		   TObject*    o, 
		   ULong_t    run, 
		   UShort_t   sys, 
		   UShort_t   sNN, 
		   Short_t    fld, 
		   Bool_t     mc, 
		   Bool_t     sat,
		   const char* file=0,
		   const char* meth="NEAR") const;
    /** 
     * Enable this correction 
     * 
     * @param enabled If true, correction is enabled
     */
    void Enable(Bool_t enabled=true) { fEnabled = enabled; }
    /** 
     * Get the data of the correction
     * 
     * @return Data of the correction
     */
    TObject* Get();
    /** 
     * Get the data of the correction
     * 
     * @return Data of the correction
     */
    const TObject* Get() const;
    /** 
     * Set the file the table is stored in
     * 
     * @param fileName file name of file table is stored in 
     */
    void SetFile(const TString& fileName) { fTitle = fileName; }
    /** 
     * Get pointer to class meta information.  Sets the internal
     * member if not done already.
     * 
     * @return Pointer to class meta information
     */
    const TClass* TheClass() const;
    /** 
     * Print information to the screen 
     * 
     * @param option Options
     */
    void Print(Option_t* option="") const;
    /** 
     * Browse this object
     * 
     * @param b Browser to use 
     */
    void Browse(TBrowser* b);
    /** 
     * Flag that this is a folder 
     * 
     * @return Always true
     */
    Bool_t IsFolder() const { return true; }
    /**
     * 
     */
    /** 
     * Massage fields for query according to settings 
     * 
     * @param run  On return, the run number 
     * @param sys  On return, the system id 
     * @param sNN  On return, the collision energy
     * @param fld  On return, the field strength 
     * @param mc   On return, the simulation flag 
     * @param sat  On return, the satellite flag 
     */
    void MassageFields(ULong_t&  run,
		       UShort_t& sys,
		       UShort_t& sNN, 
		       Short_t & fld, 
		       Bool_t&   mc, 
		       Bool_t&   sat) const;
    /** 
     * Massage fields updating according to settings 
     * 
     * @param run  On return, the run number 
     * @param sys  On return, the system id 
     * @param sNN  On return, the collision energy
     * @param fld  On return, the field strength 
     * @param mc   On return, the simulation flag 
     * @param sat  On return, the satellite flag 
     */
    void CorrectFields(ULong_t&  run,
		       UShort_t& sys,
		       UShort_t& sNN, 
		       Short_t & fld, 
		       Bool_t&   mc, 
		       Bool_t&   sat) const;
    /** 
     * Clean the table of this correction 
     * 
     * @param db    Database 
     * @param dest  Destination 
     * @param verb  Whether to be verbose or not 
     * @param all   If true, copy all entries
     * 
     * @return true on success 
     */
    Bool_t CleanIt(AliOADBForward* db, 
		   const TString& dest, 
		   Bool_t verb=false,
		   Bool_t all=false) const;
    /** 
     * Open our data base table 
     * 
     * @param db         Database 
     * @param vrb        Whether to be verbose 
     * @param fallback   Whether to enable fall-back 
     * 
     * @return true on success 
     */
    Bool_t OpenIt(AliOADBForward* db, 
		  Bool_t vrb=false, 
		  Bool_t fallback=false) const;

    mutable TClass*  fCls;      //! Class of correction objects
    const TString fClientCls; // Class name
    UShort_t fQueryFields;    // Enabled query fields
    Bool_t   fEnabled;   // Whether we're in use 
    TString  fLastEntry; // Text representation of last entry
    TObject* fObject;    // The data 
    ClassDef(Correction,1) // Correction meta object
  };
  const char* GetObjectName(Int_t what) const;

protected:
  /**
   * Constructor 
   */
  AliCorrectionManagerBase();
  /**
   * Constructor 
   */
  AliCorrectionManagerBase(Bool_t notUsed);
  /**
   * Copy Constructor 
   */
  AliCorrectionManagerBase(const AliCorrectionManagerBase& o);
  /** 
   * Assignement operator
   * 
   * @param o Object to assign from
   * 
   * @return Reference to this 
   */
  AliCorrectionManagerBase& operator=(const AliCorrectionManagerBase& o);
  /** 
   * Register a correction 
   * 
   * @param id   Identifier 
   * @param corr Correction 
   */  
  void RegisterCorrection(Int_t id, Correction* corr);
  /** 
   * Register a new correction. 
   * 
   * @param id         Identifier 
   * @param tableName  Table name 
   * @param fileName   File name 
   * @param cls        Class 
   * @param fields     Fields 
   * @param enabled    Enabled or not 
   */
  void RegisterCorrection(Int_t id, 
			  const TString& tableName, 
			  const TString& fileName, 
			  TClass*        cls,
			  UShort_t       fields=kStandard,
			  Bool_t         enabled=false);  
  /** 
   * Enable the correction at @a id
   * 
   * @param id Identifier 
   * @param enable Whether to enable (true) or disable (false)
   */
  void EnableCorrection(Int_t id, Bool_t enable=true);
  /** 
   * Get the correction at @a id
   * 
   * @param id Identifier 
   * 
   * @return Correction or null
   */
  Correction* GetCorrection(Int_t id);
  /** 
   * Get the correction at @a id
   * 
   * @param id Identifier 
   * 
   * @return Correction or null
   */
  const Correction* GetCorrection(Int_t id) const;
  /** 
   * Set the correction file 
   * 
   * @param id        Identifier 
   * @param fileName  file 
   */
  void SetCorrectionFile(Int_t id, const TString& fileName) const;
  /** 
   * Get the id of the correction with a given name
   * 
   * @param what Name of correction to look for
   * 
   * @return Correction identifier 
   */
  Int_t GetId(const TString& what) const;
  /** 
   * Get the id of the correction to store an object
   * 
   * @param obj Correction object 
   * 
   * @return Correction identifier 
   */
  Int_t GetId(const TObject* obj) const;
  /** 
   * Get the object corresponding to ID
   * 
   * @param id Correction identifier 
   * 
   * @return Object of correction, or null if correction not found or in-active
   */
  TObject* Get(Int_t id);
  /** 
   * Get the object corresponding to ID
   * 
   * @param id Correction identifier 
   * 
   * @return Object of correction, or null if correction not found or in-active
   */
  const TObject* Get(Int_t id) const;
  /** 
   * Read in all corrections 
   * 
   * @param run    Run number 
   * @param sys    System 
   * @param sNN    Center of mass energy 
   * @param fld    L3 magnetic field
   * @param mc     For simulations 
   * @param sat    For satellite interactions 
   * 
   * @return true on success 
   */
  Bool_t CheckConditions(ULong_t    run, 
			 UShort_t   sys, 
			 UShort_t   sNN, 
			 Short_t    fld, 
			 Bool_t     mc, 
			 Bool_t     sat);
  /** 
   * Read in all corrections 
   * 
   * @param run    Run number 
   * @param sys    System 
   * @param sNN    Center of mass energy 
   * @param fld    L3 magnetic field
   * @param mc     For simulations 
   * @param sat    For satellite interactions 
   * 
   * @return true on success 
   */
  Bool_t ReadCorrections(ULong_t    run, 
			 UShort_t   sys, 
			 UShort_t   sNN, 
			 Short_t    fld, 
			 Bool_t     mc, 
			 Bool_t     sat);
  /** 
   * Read in a correction
   * 
   * @param id     Correction identifier 
   * @param run    Run number 
   * @param sys    System 
   * @param sNN    Center of mass energy 
   * @param fld    L3 magnetic field
   * @param mc     For simulations 
   * @param sat    For satellite interactions 
   * 
   * @return true on success 
   */
  Bool_t ReadCorrection(Int_t      id,
			ULong_t    run, 
			UShort_t   sys, 
			UShort_t   sNN, 
			Short_t    fld, 
			Bool_t     mc, 
			Bool_t     sat);
  /** 
   * Set the correction file name for correction at @a i
   * 
   * @param i    Identifier 
   * @param file Filename 
   */
  void SetCorrectionFile(Int_t i, const TString& file);
  
  TObjArray       fCorrections; // List of corrections
  Bool_t          fIsInit;      // Whether we're intialized or not 
  ULong_t         fRun;         // Cached run number
  UShort_t        fSys;         // Cached system (1:pp 2:PbPb 3:pPb 4:Pbp 5:XeXe)
  UShort_t        fSNN;         // Cached center of mass energy [GeV]
  Short_t         fField;       // Cached L3 magnetic field [kG]
  Bool_t          fMC;          // Cached Simulation flag
  Bool_t          fSatellite;   // Cached satellite interaction flat
  AliOADBForward* fDB;          //! do not store 
  Bool_t          fDebug;       // If true, do verbose queries 
  Bool_t          fFallBack;    // If true, enable fall-back queries 
  ClassDef(AliCorrectionManagerBase,2);
};

#endif

			  
		  
    
