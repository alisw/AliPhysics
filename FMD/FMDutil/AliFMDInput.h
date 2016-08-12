#ifndef AliFMDInput_H
#define AliFMDInput_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights
 * reserved. 
 *
 * See cxx source for full Copyright notice                               
 */
//___________________________________________________________________
//
// The classes defined here, are utility classes for reading in data
// for the FMD.  They are  put in a seperate library to not polute the
// normal libraries.  The classes are intended to be used as base
// classes for customized class that do some sort of analysis on the
// various types of data produced by the FMD. 
/** @file    AliFMDInput.h
    @author  Christian Holm Christensen <cholm@nbi.dk>
    @date    Mon Mar 27 12:42:40 2006
    @brief   FMD utility classes for reading FMD data
*/
//___________________________________________________________________
/** @defgroup FMD_util Utility classes. 

    The classes defined here, are utility classes for reading in data
    for the FMD.  They are put in a seperate library to not polute the
    normal libraries.  The classes are intended to be used as base
    classes for customized class that do some sort of analysis on the
    various types of data produced by the FMD.
    
    @ingroup FMD
*/
#include <TNamed.h>
#ifndef ROOT_TString
# include <TString.h>
#endif
#ifndef ROOT_TArrayF
# include <TArrayF.h>
#endif
class AliTrackReference;
class AliRunLoader;
class AliLoader;
class AliStack;
class AliRun;
class AliRawReader;
class AliFMDRawReader;
class AliFMD;
class AliFMDHit;
class AliFMDDigit;
class AliFMDSDigit;
class AliFMDRecPoint;
class AliESDEvent;
class AliESDFMD;
class AliHeader;
class TString;
class TClonesArray;
class TTree;
class TGeoManager;
class TParticle;
class TChain;
class TSystemDirectory;

//___________________________________________________________________
/** @class AliFMDInput 
    @brief Base class for reading in various FMD data. 
    The class loops over all found events.  For each event the
    specified data is read in.   The class then loops over all
    elements of the read data, and process these with user defined
    code. 
    @code 
    struct DigitInput : public AliFMDInput 
    { 
      DigitInput() 
      { 
         // Load digits 
         AddLoad(kDigits);
	 // Make a histogram 
         fHist = new TH1F("adc", "ADC spectra", 1024, -.5, 1023.5);
      }
      // Process one digit. 
      Bool_t ProcessDigit(AliFMDDigit* d)
      { 
         fHist->Fill(d->Counts());
	 return kTRUE;
      }
      // After processing all events, display spectrum
      Bool_t Finish()
      {
        fHist->Draw();
      }
      TH1F* fHist;
    };

    void AdcSpectrum()
    { 
       DigitInput di;
       di.Run();
    }
    @endcode
    This class allows for writing small scripts, that can be compiled
    with AcLIC, to do all sorts of tests, quick prototyping, and so
    on.  It has proven to be quiet useful.   One can load more than
    one type of data in one derived class, to for example to make
    comparisons between hits and reconstructed points.   See also the
    various scripts in @c FMD/scripts. 
    @ingroup FMD_util
 */
class AliFMDInput : public TNamed
{
public:
  /** The kinds of data that can be read in. */
  enum ETrees {
    kHits       = 1,  // Hits
    kKinematics,      // Kinematics (from sim)
    kDigits,          // Digits
    kSDigits,         // Summable digits 
    kHeader,          // Header information 
    kRecPoints,       // Reconstructed points
    kESD,             // Load ESD's
    kRaw,             // Read raw data 
    kGeometry,        // Not really a tree 
    kTrackRefs,       // Track references - also for BG study
    kRawCalib,        // Read raws and calibrate them
    kUser
  };
  /** CTOR  */
  AliFMDInput();
  /** CTOR
      @param gAliceFile galice file  */
  AliFMDInput(const char* gAliceFile);
  /** DTOR */
  virtual ~AliFMDInput() {}

  /** Add a data type to load 
      @param tree Data to load */
  virtual void   AddLoad(ETrees tree)     { SETBIT(fTreeMask, tree); }
  /** Remove a data type to load 
      @param tree Data to @e not load */
  virtual void   RemoveLoad(ETrees tree)  { CLRBIT(fTreeMask, tree); }
  /** @return # of available events */
  virtual Int_t  NEvents() const;
  /** @return true if passed tree is loaded */
  virtual Bool_t IsLoaded(ETrees tree)const { return TESTBIT(fTreeMask, tree); }
  /** 
   * Set the trees to load.  
   * 
   * @param mask Bit mask of trees to load.  Should be constructed
   * like for example 
   *
   * @code 
   * UInt_t mask = ((1 << AliFMDInput::kHits) | 
   *                (1 << AliFMDInput::kDigits) | 
   *                (1 << AliFMDInput::kSDigits));
   * @endcode
   */
  virtual void SetLoads(UInt_t mask);
  /** 
   * Set the trees to load. 
   * 
   * @param mask A comma or space separated list of trees to load.
   * The case is not important, and a short from of the tree name can
   * be used.  
   */
  virtual void SetLoads(const char* mask);
  /** 
   * Get a string that represents the loaded trees 
   * 
   * @param dataOnly If true, then only show data 
   * 
   * @return String representing the loaded trees. 
   */
  virtual const char* LoadedString(Bool_t dataOnly=false) const;

  /** Initialize the class.  If a user class overloads this member
      function, then this @e must be explicitly called
      @return @c false on error */
  virtual Bool_t Init();
  /** Callled at the beginning of each event. If a user class
      overloads this member  function, then this @e must be explicitly
      called. 
      @param event Event number
      @return @c false on error */
  virtual Bool_t Begin(Int_t event);
  /** Process one event.  This loops over all the loaded data.   Users
      can overload this member function, but then it's @e strongly
      recommended to explicitly call this classes version. 
      @return @c false on error  */
  virtual Bool_t Event();
  /** Called at the end of each event. 
      @return @c false on error  */
  virtual Bool_t End();
  /** Called at the end of the run.
      @return  @c false on error */
  virtual Bool_t Finish() { return kTRUE; }
  /** Run a full job. 
      @return @c false on error  */
  virtual Bool_t Run(UInt_t maxEvents=0);

  /** Loop over all hits, and call ProcessHit with that hit, and
      optionally the corresponding kinematics track. 
      @return @c false on error  */
  virtual Bool_t ProcessHits();
  /** Loop over all track refs, and call ProcessTrackRef with that hit, and
      optionally the corresponding kinematics track. 
      @return @c false on error  */
  virtual Bool_t ProcessTrackRefs();
  /** Loop over all tracks, and call ProcessTrack with each hit for
      that track
      @return @c false on error  */
  virtual Bool_t ProcessTracks();
  /** Loop over all tracks, and call ProcessTrack with each hit for
      that track
      @return @c false on error  */
  virtual Bool_t ProcessStack();
  /** Loop over all digits, and call ProcessDigit for each digit.
      @return @c false on error  */
  virtual Bool_t ProcessDigits();
  /** Loop over all summable digits, and call ProcessSDigit for each
      digit. 
      @return @c false on error  */
  virtual Bool_t ProcessSDigits();
  /** Loop over all digits read from raw data files, and call
      ProcessRawDigit for each digit. 
      @return @c false on error  */
  virtual Bool_t ProcessRawDigits();
  /** Loop over all digits read from raw data files, and call
      ProcessRawDigit for each digit. 
      @return @c false on error  */
  virtual Bool_t ProcessRawCalibDigits();
  /** Loop over all reconstructed points, and call ProcessRecPoint for
      each reconstructed point. 
      @return @c false on error  */
  virtual Bool_t ProcessRecPoints();
  /** Loop over all ESD data, and call ProcessESD for each entry.
      @return  @c false on error  */
  virtual Bool_t ProcessESDs();
  /** Loop over all strips and ask user routine to supply the data.
      @return  @c false on error  */
  virtual Bool_t ProcessUsers();

  /** Process one hit, and optionally it's corresponding kinematics
      track.  Users should over this to process each hit. 
      @param h Hit 
      @param p Associated track
      @return  @c false on error   */
  virtual Bool_t ProcessHit(AliFMDHit* h, TParticle* p);
  /** Process one track reference, and optionally it's corresponding kinematics
      track.  Users should overload this to process each track reference. 
      @param trackRef Track Reference 
      @param track Associated track
      @return  @c false on error   */
  virtual Bool_t ProcessTrackRef(AliTrackReference* trackRef, TParticle* track);
  /** Process one hit per track. Users should over this to process
      each hit. 
      @param i Track number 
      @param p Track  
      @param h Associated Hit
      @return  @c false on error   */
  virtual Bool_t ProcessTrack(Int_t i, TParticle* p, AliFMDHit* h);
  /** Process stack particle 
      @param i Track number 
      @param p Track  
      @return  @c false on error   */
  virtual Bool_t ProcessParticle(Int_t i , TParticle* p);
  /** Process one digit.  Users should over this to process each
      digit. 
      @param digit Digit
      @return  @c false on error   */
  virtual Bool_t ProcessDigit(AliFMDDigit* digit);
  /** Process one summable digit.  Users should over this to process
      each summable digit.  
      @param sdigit Summable digit
      @return  @c false on error   */
  virtual Bool_t ProcessSDigit(AliFMDSDigit* sdigit);
  /** Process one digit from raw data files.  Users should over this
      to process each raw digit.  
      @param digit Raw digit
      @return  @c false on error   */
  virtual Bool_t ProcessRawDigit(AliFMDDigit* digit);
  /** Process one digit from raw data files.  Users should over this
      to process each raw digit.  
      @param digit Raw digit
      @return  @c false on error   */
  virtual Bool_t ProcessRawCalibDigit(AliFMDDigit* digit);
  /** Process one reconstructed point.  Users should over this to
      process each reconstructed point.  
      @param point Reconstructed point 
      @return  @c false on error   */
  virtual Bool_t ProcessRecPoint(AliFMDRecPoint* point);
  /** Process ESD data for the FMD.  Users should overload this to
      deal with ESD data. 
      @param d    Detector number (1-3)
      @param r    Ring identifier ('I' or 'O')
      @param s    Sector number (0-19, or 0-39)
      @param t    Strip number (0-511, or 0-255)
      @param eta  Psuedo-rapidity 
      @param mult Psuedo-multiplicity 
      @return  @c false on error  */
  virtual Bool_t ProcessESD(UShort_t d, Char_t r, UShort_t s, UShort_t t, 
			    Float_t eta, Float_t mult);
  /** Process User data for the FMD.  Users should overload this to
      deal with ESD data. 
      @param d    Detector number (1-3)
      @param r    Ring identifier ('I' or 'O')
      @param s    Sector number (0-19, or 0-39)
      @param t    Strip number (0-511, or 0-255)
      @param v    Value
      @return  @c false on error  */
  virtual Bool_t ProcessUser(UShort_t d, Char_t r, UShort_t s, UShort_t t, 
			     Float_t v);
  /** Service function to make a logarithmic axis. 
      @param n    Number of bins 
      @param min  Minimum of axis 
      @param max  Maximum of axis. 
      @return An array with the bin boundaries. */
  static TArrayF MakeLogScale(Int_t n, Double_t min, Double_t max);

  /** Set the raw data input 
      @param file File name - if empty, assume simulated raw. */
  void SetRawFile(const char* file) { if (file) fRawFile = file; }
  void SetInputDir(const char* dir) { fInputDir = (dir && dir[0] != '\0') 
      ? dir : "."; }
  /** 
   * Parse a string as a load option
   * 
   * @param what String to pass
   * 
   * @return Load option value, or 0 in case of error
   */
  static ETrees ParseLoad(const char* what);     
protected:
  /** Copy ctor 
      @param o Object to copy from  */
  AliFMDInput(const AliFMDInput& o) 
    : TNamed(o),
      fGAliceFile(""),
      fLoader(0),
      fRun(0),
      fStack(0),
      fFMDLoader(0),
      fReader(0),
      fFMDReader(0),
      fFMD(0),
      fESD(0),
      fESDEvent(0),
      fTreeE(0),
      fTreeH(0),
      fTreeTR(0),
      fTreeD(0),
      fTreeS(0),
      fTreeR(0),
      fTreeA(0),
      fChainE(0),
      fArrayE(0),
      fArrayH(0),
      fArrayTR(0),
      fArrayD(0),
      fArrayS(0),
      fArrayR(0),
      fArrayA(0),
      fHeader(0),
      fGeoManager(0),
      fTreeMask(0),
      fRawFile(""),      
      fInputDir("."),
      fIsInit(kFALSE),
      fEventCount(0),
      fNEvents(-1)
  {}
  /** Assignement operator 
      @return  REference to this */
  AliFMDInput& operator=(const AliFMDInput&) { return *this; }
  /** 
   * Get user supplued data
   * 
   * @param d Detector
   * @param r Ring
   * @param s Sector
   * @param t Strip
   * 
   * @return Value
   */
  virtual Float_t GetSignal(UShort_t d, Char_t r, UShort_t s, UShort_t t);

  static const char* TreeName(ETrees tree, bool shortest=false);

  /** 
   * Make a chain of specified data 
   * 
   * @param what       What data to chain.  Possible values are 
   *                   - ESD Event summary data (AliESD)
   *                   - MC  Simulation data (galice)
   * @param datadir    Data directory to scan 
   * @param recursive  Whether to recurse into sub-directories 
   * 
   * @return Pointer to newly create chain, or null
   */
  static TChain* MakeChain(const char* what, const char* datadir, 
			   bool recursive=false);
  /** 
   * Scan a directory (optionally recursive) for data files to add to
   * the chain.  Only ROOT files, and files which name contain the
   * passed pattern are considered.
   * 
   * @param dir        Directory to scan
   * @param chain      Chain to add data to 
   * @param pattern    Pattern that the file name must contain
   * @param recursive  Whether to scan recursively 
   */
  static void ScanDirectory(TSystemDirectory* dir, 
			    const TString& olddir, 
			    TChain* chain, 
			    const char* pattern, bool recursive);

  TString          fGAliceFile; // File name of gAlice file
  AliRunLoader*    fLoader;     // Loader of FMD data 
  AliRun*          fRun;        // Run information
  AliStack*        fStack;      // Stack of particles 
  AliLoader*       fFMDLoader;  // Loader of FMD data 
  AliRawReader*    fReader;     // Raw data reader 
  AliFMDRawReader* fFMDReader;  // FMD raw reader
  AliFMD*          fFMD;        // FMD object
  AliESDFMD*       fESD;        // FMD ESD data  
  AliESDEvent*     fESDEvent;   // ESD Event object. 
  TTree*           fTreeE;      // Header tree 
  TTree*           fTreeH;      // Hits tree
  TTree*           fTreeTR;     // Track Reference tree
  TTree*           fTreeD;      // Digit tree 
  TTree*           fTreeS;      // SDigit tree 
  TTree*           fTreeR;      // RecPoint tree
  TTree*           fTreeA;      // Raw data tree
  TChain*          fChainE;     // Chain of ESD's
  TClonesArray*    fArrayE;     // Event info array
  TClonesArray*    fArrayH;     // Hit info array
  TClonesArray*    fArrayTR;    // Hit info array
  TClonesArray*    fArrayD;     // Digit info array
  TClonesArray*    fArrayS;     // SDigit info array
  TClonesArray*    fArrayR;     // Rec points info array
  TClonesArray*    fArrayA;     // Raw data (digits) info array
  AliHeader*       fHeader;     // Header 
  TGeoManager*     fGeoManager; // Geometry manager
  Int_t            fTreeMask;   // Which tree's to load
  TString          fRawFile;    // Raw input file
  TString          fInputDir;   // Input directory
  Bool_t           fIsInit;     // Have we been initialized 
  Int_t            fEventCount; // Event counter 
  Int_t            fNEvents;    // The maximum number of events
  static const ETrees fgkAllLoads[kUser+1]; // List of all possible loads
  ClassDef(AliFMDInput,0)  //Hits for detector FMD
};

inline Bool_t AliFMDInput::ProcessHit(AliFMDHit*,TParticle*) { return kTRUE; }
inline Bool_t AliFMDInput::ProcessTrackRef(AliTrackReference*, 
					   TParticle*) { return kTRUE; }
inline Bool_t AliFMDInput::ProcessTrack(Int_t,TParticle*,
					AliFMDHit*) { return kTRUE; }
inline Bool_t AliFMDInput::ProcessParticle(Int_t,TParticle*) { return kTRUE; }
inline Bool_t AliFMDInput::ProcessDigit(AliFMDDigit*) { return kTRUE; }
inline Bool_t AliFMDInput::ProcessSDigit(AliFMDSDigit*) { return kTRUE; }
inline Bool_t AliFMDInput::ProcessRawDigit(AliFMDDigit*) { return kTRUE; }
inline Bool_t AliFMDInput::ProcessRawCalibDigit(AliFMDDigit*) { return kTRUE; }
inline Bool_t AliFMDInput::ProcessRecPoint(AliFMDRecPoint*) { return kTRUE; }
inline Bool_t AliFMDInput::ProcessESD(UShort_t,Char_t,UShort_t,UShort_t,
				      Float_t,Float_t) { return kTRUE; }
inline Bool_t AliFMDInput::ProcessUser(UShort_t,Char_t,UShort_t,UShort_t,
				       Float_t) { return kTRUE; }
inline Float_t AliFMDInput::GetSignal(UShort_t, Char_t, UShort_t, UShort_t) { 
  return 0.; }


#endif
//____________________________________________________________________
//
// Local Variables:
//   mode: C++
// End:
//
// EOF
//
