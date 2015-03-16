// -*- mode: C++ -*-
#ifndef ALIOADBFORWARD_H
#define ALIOADBFORWARD_H
#include <TNamed.h>
#include <TString.h>
#include <TMap.h>
class TFile;
class TTree;
class TBrowser;
class TList;

/**
 * Container/handler of Forward calibration objects. 
 */
class AliOADBForward : public TObject
{
public:
  //=== Options ======================================================
  /**
   * Options for selecting the entries according to the run number
   * given in the query. 
   */
  enum ERunSelectMode { 
    kDefault = 0,
    /** 
     * Select only entries from exactly the run given. 
     */
    kExact = 1,
    /** 
     * Select the entry with the largest run number 
     */
    kNewest = 2,
    /**
     * Select entries from the run nearest (possibly constrained to be
     * older or newer) to the run given in the query 
     */ 
    kNear = 3, 
    /** 
     * Select entries from runs that are older than the run given in
     * the query. The oldest entry with the largest run number
     * smaller than the given run is selected.  
     */
    kOlder = 4, 
    /**
     * Select entries from runs that are newer than the run given.
     * Select the entries with the smallest run number which is larger
     * than the given run.
     */
    kNewer = 5
  };
  /** 
   * Return the mode as a string 
   * 
   * @param mode Mode 
   * 
   * @return Stringified mode 
   */
  static const char* Mode2String(ERunSelectMode mode);
  /** 
   * Parse a string to get the mode 
   * 
   * @param str Input string 
   * 
   * @return Mode 
   */
  static ERunSelectMode String2Mode(const TString& str);
  /** 
   * Return mode as an integer 
   * 
   * @param mode Mode 
   * 
   * @return Mode identifier 
   */
  static ERunSelectMode Int2Mode(Int_t mode);
  /**
   * Various flags
   * 
   */
  enum { 
    /** 
     * The maximum distance to the given run when selecting in kNear
     * mode.  Currently this is set to a large number to allow
     * selection of any entry.  This should change 
     */
    kMaxNearDistance = 1000000,
    kInvalidField    = 999
  };
  //=== Entry ========================================================
  /**
   * An entry in the FORWARD OADB database tables
   * 
   */
  class Entry : public TObject
  {
  public:
    /** 
     * Constructor 
     * 
     * @param runNo   Run number 
     * @param sys     Collision system (1:pp, 2:pbpb, 3:ppb) 
     * @param sNN     Center of mass energy (GeV)
     * @param field   L3 magnetic field (kG)
     * @param mc      True if for MC only
     * @param sat     For satellite events
     * @param o       Correction object
     */
    Entry(ULong_t  runNo  = 0, 
	  UShort_t sys    = 0, 
	  UShort_t sNN    = 0, 
	  Short_t  field  = 0, 
	  Bool_t   mc     = false, 
	  Bool_t   sat    = false,
	  TObject* o      = 0);
    /** 
     * Copy constructor
     * 
     * @param o Object to copy from
     */
    Entry(const Entry& o);
    /** 
     * Destructor 
     */
    virtual ~Entry() {}
    /** 
     * Assignment operator 
     * 
     * @param o Object to assign from 
     * 
     * @return Reference to this object
     */
    Entry& operator=(const Entry& o);
    /** 
     * Get the title of the entry 
     * 
     * @return stringified fields
     */
    const char* GetTitle() const;
    /** 
     * Print this entry 
     * 
     * @param option Not used 
     */
    void Print(Option_t* option="") const; //*MENU*
    ULong_t  fRunNo;           // Run number 
    UShort_t fSys;             // Collision system (1: pp, 2: pbpb, 3: ppb)
    UShort_t fSNN;             // Center of mass energy 
    Short_t  fField;           // L3 magnetic field
    Bool_t   fMC;              // True if only for MC 
    Bool_t   fSatellite;       // Satelitte events
    TObject* fData;            // Correction object
    UInt_t   fTimestamp;       // When the object was stored 
    ULong_t  fAliROOTRevision; // Revision of AliROOT used 
    TString  fAuthor;          // Author of calibration/correction

    ClassDef(Entry,2); // Entry in PWGLF/Forward OADB
  };

  //=== Table ========================================================
  /** 
   * A table on the Forward OADB - the underlying storage is a TTree
   * containing Entry objects.
   * 
   */
  class Table : public TObject
  {
  public:
    /** 
     * Constructor 
     * 
     * @param tree   Tree 
     * @param isNew  Whether to make the branch 
     * @param mode   How to select on the run number 
     * 
     * @return 
     */    
    Table(TTree* tree, Bool_t isNew, ERunSelectMode mode=kNear);
    /** 
     * Copy constructor 
     * 
     * @param other Object to copy from 
     */
    Table(const Table& other);
    /** 
     * Destructor. 
     *
     * Closes the corresponding file 
     */
    ~Table();
    /** 
     * Assignemt operator
     * 
     * @param other Object to assign form 
     * 
     * @return Reference to this object 
     */
    Table& operator=(const Table& other);
    /** 
     * Set the verbosity
     * 
     * @param verb 
     */
    void SetVerbose(Bool_t verb=true) { fVerbose = verb; }
    /** 
     * Set wheter to enable fall-back queries  
     * 
     * @param use If true, enable fall-back queries
     */
    void SetEnableFallBack(Bool_t use=true) { fFallBack = use; }
    // -----------------------------------------------------------------
    /** 
     * Get the name of the tree 
     * 
     * @return Table name 
     */
    const Char_t* GetTableName() const;
    /** 
     * Overload TObject::GetName 
     * 
     * @return Name 
     */
    const Char_t* GetName() const;

    // -----------------------------------------------------------------
    /** 
     * @{ 
     * @name Open/close/flush 
     */
    /** 
     * Flush to disk 
     * 
     * @return true on success
     */
    Bool_t Update();
    /** 
     * Close connection
     * 
     * @return true on success 
     */
    Bool_t Close();
    /* @} */
    // -----------------------------------------------------------------
    /** 
     * @{ 
     * @name Queries 
     */
    /** 
     * Query the tree 
     * 
     * @param runNo  Run number 
     * @param mode   Run selection mode 
     * @param sys    Collision system (1: pp, 2: PbPb, 3: pPb)
     * @param sNN    Center of mass energy (GeV)
     * @param fld    L3 magnetic field (kG)
     * @param mc     For MC only 
     * @param sat    For satellite events
     * @param mode   How to select on the run number 
     * 
     * @return Found entry number or negative number in case of problems
     */
    Int_t Query(ULong_t        runNo  = 0,
		ERunSelectMode mode   = kNear,
		UShort_t       sys    = 0,
		UShort_t       sNN    = 0, 
		Short_t        fld    = kInvalidField,
		Bool_t         mc     = false,
		Bool_t         sat    = false) const;
    /** 
     * Run a query with pre-build conditions
     * 
     * @param q     query string 
     * @param runNo The given run number 
     * @param mode  Run selection mode 
     * 
     * @return Entry number of selected entry 
     */
    Int_t Query(ULong_t        runNo,
		ERunSelectMode mode,
		const TString& q) const;
    /** 
     * Insert a new entry into the tree 
     * 
     * @param o         Object to write 
     * @param runNo     Run number 
     * @param sys       Collision system (1: pp, 2:PbPb, 3:pPb)
     * @param sNN       Center of mass energy (GeV)
     * @param field     L3 magnetic field (kG)
     * @param mc        If true, only for MC 
     * @param sat       For satellite interactions 
     * @param aliRev    AliROOT revision
     * @param author    Creater of this correction
     * 
     * @return true on success 
     */
    Bool_t Insert(TObject* o, 
		  ULong_t  runNo, 
		  UShort_t sys, 
		  UShort_t sNN, 
		  Short_t  field, 
		  Bool_t   mc=false, 
		  Bool_t   sat=false,
		  ULong_t  aliRev=0,
		  const TString& author="");
    /** 
     * Query the tree for an object.  The strategy is as follows. 
     * 
     *  - First query with all fields 
     *    - If this returns a single entry, return that. 
     *    - If not, then ignore the run number (if given)
     *      - If this returns a single entry, return that 
     *      - If not, and fall-back is enabled, then 
     *        - Ignore the collision energy (if given) 
     *          - If this returns a single entry, return that.  
     *          - If not, ignore all passed values 
     *            - If this returns a single entry, return that.
     *            - Otherwise, give up and return -1
     *      - Otherwise, give up and return -1
     *
     * This allow us to specify default objects for a period, and for
     * collision system, energy, and field setting.
     *
     * @param run    Run number 
     * @param mode   Run selection mode 
     * @param sys    Collision system (1: pp, 2: PbPb, 3: pPb)
     * @param sNN    Center of mass energy (GeV)
     * @param fld    L3 magnetic field (kG)
     * @param mc     For MC only 
     * @param sat    For satellite interactions 
     * 
     * @return Found entry number, or -1
     */
    Int_t GetEntry(ULong_t        run      = 0,
		   ERunSelectMode mode     = kNear,
		   UShort_t       sys      = 0,
		   UShort_t       sNN      = 0, 
		   Short_t        fld      = 0,
		   Bool_t         mc       = false,
		   Bool_t         sat      = false) const;
    /** 
     * Query the tree for an object.  The strategy is as follows. 
     * 
     *  - First query with all fields 
     *    - If this returns a single entry, return that. 
     *    - If not, then ignore the run number (if given)
     *      - If this returns a single entry, return that 
     *      - If not, and fall-back is enabled, then 
     *        - Ignore the collision energy (if given) 
     *          - If this returns a single entry, return that.  
     *          - If not, ignore all passed values 
     *            - If this returns a single entry, return that.
     *            - Otherwise, give up and return null
     *      - Otherwise, give up and return null
     *
     * This allow us to specify default objects for a period, and for
     * collision system, energy, and field setting.
     *
     * @param run    Run number 
     * @param mode   Run selection mode 
     * @param sys    Collision system (1: pp, 2: PbPb, 3: pPb)
     * @param sNN    Center of mass energy (GeV)
     * @param fld    L3 magnetic field (kG)
     * @param mc     For MC only 
     * @param sat    For satellite interactions 
     * 
     * @return Found entry, or null
     */
    Entry* Get(ULong_t        run      = 0,
	       ERunSelectMode mode     = kNear,
	       UShort_t       sys      = 0,
	       UShort_t       sNN      = 0, 
	       Short_t        fld      = 0,
	       Bool_t         mc       = false,
	       Bool_t         sat      = false) const;
    /** 
     * Query the tree for an object.  The strategy is as follows. 
     * 
     *  - First query with all fields 
     *    - If this returns a single entry, return that. 
     *    - If not, then ignore the run number (if given)
     *      - If this returns a single entry, return that 
     *      - Otherwise, give up and return null
     *
     * This allow us to specify default objects for a period, and for
     * collision system, energy, and field setting.
     *
     * @param run    Run number 
     * @param mode   Run selection mode 
     * @param sys    Collision system (1: pp, 2: PbPb, 3: pPb)
     * @param sNN    Center of mass energy (GeV)
     * @param fld    L3 magnetic field (kG)
     * @param mc     For MC only 
     * @param sat    For satellite interactions 
     * 
     * @return Found data, or null
     */
    TObject* GetData(ULong_t        run      = 0,
		     ERunSelectMode mode     = kNear,
		     UShort_t       sys      = 0,
		     UShort_t       sNN      = 0, 
		     Short_t        fld      = 0,
		     Bool_t         mc       = false,
		     Bool_t         sat      = false) const;
    /* @} */
    /** 
     * Print the contents of the tree 
     * 
     * @param option Passed to TTree::Scan
     */
    void Print(Option_t* option="") const; //*MENU*
    /** 
     * Browse this table 
     */
    void Browse(TBrowser* b); 
    /** 
     * Check if the tree was opened. 
     * 
     * @param rw If true, also check if the file is read/write 
     *
     * @return true if everything is dandy
     */
    Bool_t IsOpen(Bool_t rw=false) const; 

    TTree*         fTree;     // Our tree
    Entry*         fEntry;    // Entry cache 
    Bool_t         fVerbose;  // To be verbose or not 
    ERunSelectMode fMode;     // Run query mode 
    Bool_t         fFallBack; // Enable fall-back

    ClassDef(Table,1); 
  };
  // === Interface ===================================================
  /** 
   * Constructor 
   */
  AliOADBForward();
  /**
   * Destructor
   * 
   */
  ~AliOADBForward();
  // --- Open/close/update -------------------------------------------
  /** 
   * @{ 
   * @name Input/output 
   */
  /** 
   * Open a file containing tables.  Note, this member function can be
   * called multiple times to open tables in different files.
   * 
   * If a table is already associated with this handler, it will not
   * be re-associated.
   *
   * @param fileName  Path to file to get/write tables from/in
   * @param rw        if true, open read+write, otherwise read-only
   * @param tables    Tables to open 
   * @param verb      Verbosity flag 
   * @param fallback  If true allow for fall-backs
   * 
   * @return true on success 
   */
  Bool_t Open(const TString& fileName, 
	      const TString& tables  = "*", 
	      Bool_t         rw      = false, 
	      Bool_t         verb    = false,
	      Bool_t         fallback= false);
  /** 
   * Open a file containing tables.  Note, this member function can be
   * called multiple times to open tables in different files.
   *
   * If a table is already associated with this handler, it will not
   * be re-associated.
   * 
   * @param file    File to get/write tables from/in
   * @param rw      if true, open read+write, otherwise read-only
   * @param tables  Tables to open 
   * @param verb    Verbosity flag 
   * @param fallback  If true allow for fall-backs
   * 
   * @return true on success 
   */
  Bool_t Open(TFile*         file,
	      const TString& tables, 
	      Bool_t         rw      = false, 
	      Bool_t         verb    = false,
	      Bool_t         fallback= false);
  /** 
   * Close this database 
   * 
   * 
   * @return true on success
   */
  Bool_t Close();
  /** 
   * Flush to disk 
   * 
   * @return true on success
   */
  Bool_t Update();
  /* @} */
  // --- Queries -----------------------------------------------------
  /** 
   * @{ 
   * @name Queries 
   */
  /** 
   * Query the table for an object.  The strategy is as follows. 
   * 
   *  - First query with all fields 
   *    - If this returns a single entry, return that. 
   *    - If not, then ignore the run number (if given)
   *      - If this returns a single entry, return that 
   *      - Otherwise, give up and return null
   *
   * This allow us to specify default objects for a period, and for
   * collision system, energy, and field setting.
   *
   * @param table  Table name 
   * @param run    Run number 
   * @param mode   Run selection mode 
   * @param sys    Collision system (1: pp, 2: PbPb, 3: pPb)
   * @param sNN    Center of mass energy (GeV)
   * @param fld    L3 magnetic field (kG)
   * @param mc     For MC only
   * @param sat    For satellite interactions 
   * 
   * @return Found entry, or null
   */
  Entry* Get(const TString& table, 
	     ULong_t        run        = 0,
	     ERunSelectMode mode       = kNear, 
	     UShort_t       sys        = 0,
	     UShort_t       sNN        = 0, 
	     Short_t        fld        = 0,
	     Bool_t         mc         = false,
	     Bool_t         sat        = false) const;
  /** 
   * Query the table for an object.  The strategy is as follows. 
   * 
   *  - First query with all fields 
   *    - If this returns a single entry, return that. 
   *    - If not, then ignore the run number (if given)
   *      - If this returns a single entry, return that 
   *      - Otherwise, give up and return null
   *
   * This allow us to specify default objects for a period, and for
   * collision system, energy, and field setting.
   *
   * @param table  Table name 
   * @param run    Run number 
   * @param mode   Run selection mode 
   * @param sys    Collision system (1: pp, 2: PbPb, 3: pPb)
   * @param sNN    Center of mass energy (GeV)
   * @param fld    L3 magnetic field (kG)
   * @param mc     For MC only 
   * @param sat    For satellite interactions 
   * 
   * @return Found data, or null
   */
  TObject* GetData(const TString& table, 
		   ULong_t        run      = 0,
		   ERunSelectMode mode     = kNear, 
		   UShort_t       sys      = 0,
		   UShort_t       sNN      = 0, 
		   Short_t        fld      = 0,
		   Bool_t         mc       = false,
		   Bool_t         sat      = false) const;
  // --- Insert ------------------------------------------------------
  /** 
   * Insert a new entry into the table 
   * 
   * @param table     Table name 
   * @param o         Object to write 
   * @param runNo     Run number 
   * @param sys       Collision system (1: pp, 2:PbPb, 3:pPb)
   * @param sNN       Center of mass energy (GeV)
   * @param field     L3 magnetic field (kG)
   * @param mc        If true, only for MC 
   * @param sat    For satellite interactions 
   * @param aliRev    AliROOT revision
   * @param author    Creater of this correction
   * 
   * @return true on success 
   */
  Bool_t Insert(const TString& table, 
		TObject*       o, 
		ULong_t        runNo, 
		UShort_t       sys, 
		UShort_t       sNN, 
		Short_t        field, 
		Bool_t         mc=false, 
		Bool_t         sat=false,
		ULong_t        aliRev=0,
		const TString& author="");
  /** 
   * Copy one entry to another entry 
   * 
   * @param table      Table name 
   * @param oldRunNo   Old run number 
   * @param oldSys     Old collision system
   * @param oldSNN     Old center of mass energy 
   * @param oldField   Old L3 magnetic field strength
   * @param newRunNo   New run number 		     
   * @param newSys     New collision system	     
   * @param newSNN     New center of mass energy     
   * @param newField   New L3 magnetic field strength
   * @param mc         True for MC only queries
   * @param sat        True for including satellite queries
   * 
   * @return true on success 
   */
  Bool_t CopyEntry(const TString& table, 
		   ULong_t        oldRunNo, 
		   UShort_t       oldSys,
		   UShort_t       oldSNN, 
		   Short_t        oldField, 
		   ULong_t        newRunNo,
		   UShort_t       newSys,
		   UShort_t       newSNN, 
		   Short_t        newField, 
		   Bool_t         mc, 
		   Bool_t         sat);
  /* @} */
  /** 
   * Print the content of all tables
   * 
   * @param option Passed on to tables 
   */  
  void Print(const Option_t* option="") const; //*MENU*
  /** 
   * Browse this database
   * 
   * @param b Browser 
   */
  void Browse(TBrowser* b); 
  /** 
   * Find a table by name 
   * 
   * @param name  Name of table 
   * @param quite Do not print warning if not found 
   * @return Table or null 
   */
  Table* FindTable(const TString& name, Bool_t quite=false) const;
  /** 
   * Get all tables 
   * 
   * 
   * @return Map of all tables
   */
  const TMap& GetTables() const { return fTables; }
protected:
  /** 
   * Get a list of associated files 
   * 
   * @param files On return, contains list of files 
   * 
   * @return Number of associated files found 
   */
  Int_t GetFiles(TList& files) const;
  /** 
   * Get a table (TTree) from a file 
   * 
   * @param file   File to look in 
   * @param rw     If true, open read/write, false read-only
   * @param name   Name of the table 
   * @param mode   Default mode for new table, or override mode 
   *               for existing tables if not default
   * 
   * @return Table or (if rw=true) possibly newly created table
   */
  Table* GetTableFromFile(TFile* file, Bool_t rw, 
			  const TString& name,
			  const TString& mode) const;

  /** 
   * Get a table (TTree) from a file and attach
   * 
   * @param file   File to look in 
   * @param rw     If true, open read/write, false read-only
   * @param name   Name of the table 
   * @param mode   Default mode for new table, or override mode 
   *               for existing tables if not default
   * @param fallback Enable fall-backs
   * @param verb If true, be verbose 
   * 
   * @return Table or (if rw=true) possibly newly created table
   */
  void OpenTable(TFile* file, Bool_t rw, const TString& name, 
		 const TString& mode, Bool_t verb, Bool_t fallback);
  /** 
   * Helper function to append to query string 
   * 
   * @param q         String to attach to
   * @param s         What to attach 
   * @param andNotOr  If true, assume @b and, otherwise @b or
   */
  static void AppendToQuery(TString& q, const TString& s, Bool_t andNotOr=true);
  /** 
   * Helper function to build a query string 
   * 
   * @param sys    Collision system (1:pp, 2: PbPb, 3:pPb)
   * @param sNN    Collision energy in GeV
   * @param fld    L3-Field strength and polarity 
   * @param mc     For Monte-carlo
   * @param sat    For satelitte collisions 
   * 
   * @return Query string 
   */
  static TString Conditions(UShort_t       sys    = 0,
			    UShort_t       sNN    = 0, 
			    Short_t        fld    = kInvalidField,
			    Bool_t         mc     = false,
			    Bool_t         sat    = false);
  TMap fTables;
  ClassDef(AliOADBForward,0); // PWGLF/Forward OADB interface

public:
  // =================================================================
  /** 
   * @{ 
   * @name Tests
   */
  static void TestGet(AliOADBForward& t, 
		      const TString& table,
		      ULong_t  runNo      = 0,
		      ERunSelectMode mode = kNear,
		      UShort_t sys        = 2,
		      UShort_t sNN        = 2760, 
		      Short_t  fld        = -5,
		      Bool_t   mc         = false,
		      Bool_t   sat        = false);
  static void TestInsert(AliOADBForward& t, 
			 const TString&  table,
			 ULong_t         runNo  = 0,
			 UShort_t        sys    = 2,
			 UShort_t        sNN    = 2760, 
			 Short_t         fld    = -5,
			 Bool_t          mc     = false,
			 Bool_t          sat    = false);
  static void Test();
  /* @} */
};

#endif
/* Local Variables:
 *  mode: C++
 * End:
 */
