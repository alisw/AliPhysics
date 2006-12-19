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
*/
#include <TObject.h>
#ifndef ROOT_TString
# include <TString.h>
#endif
class AliRunLoader;
class AliLoader;
class AliStack;
class AliRun;
class AliRawReader;
class AliFMD;
class AliFMDHit;
class AliFMDDigit;
class AliFMDSDigit;
class AliFMDRecPoint;
class AliESD;
class AliESDFMD;
class TString;
class TClonesArray;
class TTree;
class TGeoManager;
class TParticle;
class TChain;

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
class AliFMDInput : public TObject
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
    kGeometry         // Not really a tree 
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
  virtual Bool_t Run();

  /** Loop over all hits, and call ProcessHit with that hit, and
      optionally the corresponding kinematics track. 
      @return @c false on error  */
  virtual Bool_t ProcessHits();
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
  /** Loop over all reconstructed points, and call ProcessRecPoint for
      each reconstructed point. 
      @return @c false on error  */
  virtual Bool_t ProcessRecPoints();
  /** Loop over all ESD data, and call ProcessESD for each entry.
      @return  @c false on error  */
  virtual Bool_t ProcessESDs();

  /** Process one hit, and optionally it's corresponding kinematics
      track.  Users should over this to process each hit. 
      @param h Hit 
      @param p Associated track
      @return  @c false on error   */
  virtual Bool_t ProcessHit(AliFMDHit* h, TParticle* p);
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
  virtual Bool_t ProcessESD(UShort_t, Char_t, UShort_t, UShort_t, 
			    Float_t, Float_t);
  
protected:
  /** Copy ctor 
      @param o Object to copy from  */
  AliFMDInput(const AliFMDInput& o) 
    : TObject(o),
      fGAliceFile(""),
      fLoader(0),
      fRun(0),
      fStack(0),
      fFMDLoader(0),
      fReader(0),
      fFMD(0),
      fMainESD(0),
      fESD(0),
      fTreeE(0),
      fTreeH(0),
      fTreeD(0),
      fTreeS(0),
      fTreeR(0),
      fTreeA(0),
      fChainE(0),
      fArrayE(0),
      fArrayH(0),
      fArrayD(0),
      fArrayS(0),
      fArrayR(0),
      fArrayA(0),
      fGeoManager(0),
      fTreeMask(0),
      fIsInit(kFALSE)
  {}
  /** Assignement operator 
      @return  REference to this */
  AliFMDInput& operator=(const AliFMDInput&) { return *this; }

  TString       fGAliceFile; // File name of gAlice file
  AliRunLoader* fLoader;     // Loader of FMD data 
  AliRun*       fRun;        // Run information
  AliStack*     fStack;      // Stack of particles 
  AliLoader*    fFMDLoader;  // Loader of FMD data 
  AliRawReader* fReader;     // Raw data reader 
  AliFMD*       fFMD;        // FMD object
  AliESD*       fMainESD;    // ESD Object
  AliESDFMD*    fESD;        // FMD ESD data  
  TTree*        fTreeE;      // Header tree 
  TTree*        fTreeH;      // Hits tree
  TTree*        fTreeD;      // Digit tree 
  TTree*        fTreeS;      // SDigit tree 
  TTree*        fTreeR;      // RecPoint tree
  TTree*        fTreeA;      // Raw data tree
  TChain*       fChainE;     // Chain of ESD's
  TClonesArray* fArrayE;     // Event info array
  TClonesArray* fArrayH;     // Hit info array
  TClonesArray* fArrayD;     // Digit info array
  TClonesArray* fArrayS;     // SDigit info array
  TClonesArray* fArrayR;     // Rec points info array
  TClonesArray* fArrayA;     // Raw data (digits) info array
  TGeoManager*  fGeoManager; // Geometry manager
  Int_t         fTreeMask;   // Which tree's to load
  Bool_t        fIsInit;     // Have we been initialized 
  ClassDef(AliFMDInput,0)  //Hits for detector FMD
};

inline Bool_t AliFMDInput::ProcessHit(AliFMDHit*,TParticle*) { return kTRUE; }
inline Bool_t AliFMDInput::ProcessDigit(AliFMDDigit*) { return kTRUE; }
inline Bool_t AliFMDInput::ProcessSDigit(AliFMDSDigit*) { return kTRUE; }
inline Bool_t AliFMDInput::ProcessRawDigit(AliFMDDigit*) { return kTRUE; }
inline Bool_t AliFMDInput::ProcessRecPoint(AliFMDRecPoint*) { return kTRUE; }
inline Bool_t AliFMDInput::ProcessESD(UShort_t,Char_t,UShort_t,UShort_t,
				      Float_t,Float_t) { return kTRUE; }


#endif
//____________________________________________________________________
//
// Local Variables:
//   mode: C++
// End:
//
// EOF
//
