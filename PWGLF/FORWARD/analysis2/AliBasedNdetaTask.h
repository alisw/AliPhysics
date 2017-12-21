//
// Task to analyse the AOD for for dN/deta in the base regions 
//
#ifndef ALIBASEDNDETATASK_H
#define ALIBASEDNDETATASK_H
/**
 * @file   AliBasedNdetaTask.h
 * @author Christian Holm Christensen <cholm@dalsgaard.hehi.nbi.dk>
 * @date   Wed Mar 23 13:58:12 2011
 * 
 * @brief  
 * 
 * @ingroup pwglf_forward_dndeta
 * 
 */
#include "AliBaseAODTask.h"
#include <AliAnalysisUtils.h>
class TAxis;
class TList;
class TH2D;
class TH2F;
class TH1D;
class TH1I;
class TF1;
class TProfile;
class AliAODEvent;
class AliAODForwardMult;
class TObjArray;

/** 
 * @defgroup pwglf_forward_tasks_dndeta dN/deta tasks 
 *
 * Code to produce @f$ dN/d\eta@f$
 *
 * @ingroup pwglf_forward_tasks 
 */
/**
 * @defgroup pwglf_forward_dndeta dN/deta
 *
 * @f$ dN/d\eta@f$ code 
 *
 * @ingroup pwglf_forward_topical
 */
/**
 * Base class for tasks to determine @f$ dN/d\eta@f$ 
 *
 * @ingroup pwglf_forward_dndeta
 */
class AliBasedNdetaTask : public AliBaseAODTask
{
public:
  /** 
   * Bit mask values of the normalisation scheme 
   */
  enum {
    /** Only normalize to accepted events */
    kNone = 0,
    /** 
     * Do the full normalisation 
     * @f[ 
     *   N = (N_A-N_A/N_V(N_T-N_V)) = \frac{1}{\epsilon_V}N_A
     * @f]
     */
    kEventLevel = 0x1,
    /** 
     * Dummy 
     */
    kDummy = 0x2, 
    /** 
     * Correct for background events (A+C-E)
     */
    kBackground = 0x4,
    /**
     * Correct for the trigger efficiency
     *
     * @f[
     *   N = N_A \frac{1}{\epsilon_X}
     * @f]
     */
    kTriggerEfficiency = 0x8,
    /** 
     * Correct using zero-bin efficiency only 
     */
    kZeroBin = 0x10,
    /**
     * Do the full correction
     */
    kFull = kEventLevel | kBackground | kTriggerEfficiency,
  };
  /**
   * Mask for selecting pile-up 
   */
  enum EPileupMask { 
    /**
     * Use the flag from AOD 
     */
    kPileupNormal = 0, 
    /** 
     * Check the pile-up flag from SPD 
     */
    kPileupSPD    = 0x1, 
    /** 
     * Check the pileup flag from tracks 
     */
    kPileupTrk    = 0x2, 
    /**
     * Check the out-of-bunch pileup flag 
     */
    kPileupBC     = 0x4, 
    /**
     * Use the flag from AOD 
     */
    kPileupFull   = 0x8,
    /** 
     * Also accept pileup 
     */
    kPileupIgnore = 0x10,
    /** 
     * Use analysis utility class 
     */
    kPileupUtil = 0x20,
    /** 
     * Also accept pileup 
     */
    kPileupBins = 0x40,
  };

  enum ECentralityEstimator { 
    kCentNone,
    kCentDefault,   // What ever stored in AOD 
    kCentV0M,       // VZERO multiplcity 
    kCentV0A,       // VZERO-A 
    kCentV0A123,    // VZERO-A 
    kCentV0C,       // VZERO-C
    kCentFMD,       // FMD
    kCentTrk,       // Tracks
    kCentTkl,       // Tracklets
    kCentCL0,       // Clusters in SPD-0
    kCentCL1,       // Clusters in SPD-1
    kCentCND,       // Candle - tracks+tracklets
    kCentZNA,       // ZDC neutrons A-side
    kCentZNC,       // ZDC neutrons C-side
    kCentZPA,       // ZDC protons A-side
    kCentZPC,       // ZDC protons C-side
    kCentNPA,       // ?
    kCentV0MvsFMD,  // V0M vs FMD
    kCentTklvsV0M,  // Tracks vs V0M
    kCentZEMvsZDC,  // ZDC veto vs neutrons
    kMult,          // Reference multiplicity in @f$|\eta|<0.8@f$ 
    kMultV0A,       // HMTF V0A
    kMultV0M,       // HMTF V0M
    kMultV0C,       // HMTF V0C 
    kCentTrue = 0x100, 
    kCentEq   = 0x200
  };
  /** 
   * Constructor 
   * 
   */
  AliBasedNdetaTask();
  /** 
   * Constructor
   * 
   * @param name    Name of task 
   */
  AliBasedNdetaTask(const char* name);
  /**
   * Destructor
   * 
   */
  virtual ~AliBasedNdetaTask();
  
  /** 
   * @{ 
   * @name Task configuration 
   */
  /** 
   * Set the debug level 
   * 
   * @param level Debug level
   */
  virtual void SetDebugLevel(Int_t level);
  /** 
   * Set whether to correct for empty bins when projecting on the X axis. 
   * 
   * @param use Whether to correct for empty bins 
   */
  void SetCorrEmpty(Bool_t use) { fCorrEmpty = use; }
  /** 
   * Set whether to use the ROOT TH2::ProjectionX method when
   * projecting on the X axis.
   * 
   * @param use Whether to use TH2::ProjectionX
   */
  void SetUseROOTProjectX(Bool_t use) { fUseROOTProj = use; }
  void SetUseUtilPileup(Bool_t use) { fUseUtilPileup = use; }
  /** 
   * Trigger efficiency for selected trigger(s)
   * 
   * @param e Trigger efficiency 
   */
  void SetTriggerEff(Double_t e) { fTriggerEff = e; } 
  /** 
   * Trigger efficiency for 0-bin for selected trigger(s)
   * 
   * @param e Trigger efficiency for 0-bin
   */
  void SetTriggerEff0(Double_t e) { fTriggerEff0 = e; } 
  /**
   * Set satellite vertex flag
   *
   * @param satVtx
   */
  void SetSatelliteVertices(Bool_t satVtx) { fSatelliteVertices = satVtx; }
  /** 
   * Set which centrality estimator to use - if not set, use the one
   * from the Forward AOD object.  Note, the string is diagnosed, and
   * if found not to be valid, then a the program terminates via a
   * Fatal.
   *
   * @param method String definining centrality method (case insensitive)
   * Typical values are 
   * - V0M (e.g., PbPb)
   * - V0A 
   * - V0C 
   * - FMD 
   * - ZNA (e.g., pPb) 
   * - ZNC (e.g., Pbp)
   * - ZPA 
   * - ZPC 
   * - ZEMvsZDC
   *
   * @return true if @a method is valid estimator 
   */
  Bool_t SetCentralityMethod(const TString& method);
  /**
   * Get reference to the analysis utility 
   *
   * @return reference to the AliAnalysisUtils object
   */
  AliAnalysisUtils& GetAnalysisUtils() { return fAnaUtil; }
  /** 
   * Get a string representing the normalization scheme 
   * 
   * @param scheme Normalization scheme bits 
   * 
   * @return String representation 
   */
  static const Char_t* NormalizationSchemeString(UShort_t scheme);
  /** 
   * Parse a string representing the normalization scheme 
   * 
   * @param what String of the normalization scheme 
   * 
   * @return normalization scheme bits
   */
  static UShort_t ParseNormalizationScheme(const Char_t* what);
  /** 
   * Setthe normalisation scheme to use 
   * 
   * @param scheme Normalisation scheme 
   */
  void SetNormalizationScheme(UShort_t scheme);
  /** 
   * Space, pipe, or comma separated list of options
   * 
   * @param what List of options 
   */
  void SetNormalizationScheme(const char* what);
  /** 
   * Filename of final MC correction
   * 
   * @param filename filename
   */
  void SetMCFinalCorrFilename(const char* filename) { 
    fFinalMCCorrFile.Clear();
    fFinalMCCorrFile.Append(filename); 
  }
  /** @} */
  /** 
   * Print information 
   * 
   * @param option If it contains R, print recursive
   */
  void Print(Option_t* option="R") const;
  /** @{ 
   *  @name Task interface 
   */
  /** 
   * Create output objects.  
   *
   * This is called once per slave process 
   *
   * @return true on success
   */
  virtual Bool_t Book();
  /** 
   * Called before processing a single event - should not do anything
   * but clear data, etc.
   * 
   * @return true on success
   */
  virtual Bool_t PreEvent() {
    fCacheCent = -10; fCacheQual = 0xFFFF; return true; } 
  /** 
   * Process a single event 
   * 
   * @return true on success
   */
  virtual Bool_t Event(AliAODEvent& aod);
  /** 
   * Called at end of event processing.
   *
   * This is called once in the master 
   * 
   * @return true on success
   */
  virtual Bool_t Finalize();
  /* @} */

  /** 
   * @{ 
   * @name Services member functions 
   */
  /** 
   * Project onto the X axis 
   * 
   * @param h         2D histogram 
   * @param name      New name 
   * @param firstbin  First bin to use 
   * @param lastbin   Last bin to use
   * @param useROOT   Use TH2::ProjectionX instead of custom code 
   * @param corr      Whether to do corrections or not 
   * @param error     Whether to calculate errors
   * 
   * @return Newly created histogram or null
   */
  static TH1D* ProjectX(const TH2D* h, 
			const char* name,
			Int_t firstbin, 
			Int_t lastbin, 
			bool  useROOT=false,
			bool  corr=true,
			bool  error=true);
  /** 
   * Scale the copy of the 2D histogram by coverage in supplied 1D histogram
   *  
   * @param copy Data to scale 
   * @param norm Coverage histogram 
   */
  static void ScaleToCoverage(TH2D* copy, const TH1D* norm);
  /** 
   * Scale the copy of the 1D histogram by coverage in supplied 1D histogram
   *  
   * @param copy Data to scale 
   * @param norm Coverage histogram 
   */
  static void ScaleToCoverage(TH1D* copy, const TH1D* norm);
  /** 
   * Set histogram graphical options, etc. 
   * 
   * @param h       Histogram to modify
   * @param colour  Marker color 
   * @param marker  Marker style
   * @param title   Title of histogram
   * @param ytitle  Title on y-axis. 
   */
  static void SetHistogramAttributes(TH1D* h, Int_t colour, Int_t marker, 
				     const char* title, 
				     const char* ytitle=0);
  /** @} */

  /**
   * Marker styles 
   */
  enum { 
    kSolid        = 0x000, 
    kHollow       = 0x001, 
    kCircle       = 0x002,
    kSquare       = 0x004, 
    kUpTriangle   = 0x006, 
    kDownTriangle = 0x008, 
    kDiamond      = 0x00a,
    kCross        = 0x00c,
    kStar         = 0x00e
  };
  /** 
   * Get the marker style from option bits
   * 
   * @param bits Option bits 
   * 
   * @return Marker style 
   */
  static Int_t GetMarkerStyle(UShort_t bits);
  /** 
   * Get the marker option bits from a style 
   * 
   * @param style Style
   * 
   * @return option bits
   */
  static UShort_t GetMarkerBits(Int_t style);
  /** 
   * Flip an option bit 
   * 
   * @param style Style parameter
   * 
   * @return New style 
   */
  static Int_t FlipHollowStyle(Int_t style);
  /**
   * Setter of empirical correction
   *
   * @param h 2D histogram of ratio of nominal @f$ 1/N
   * dN_{ch}/d\eta@f$ to satellite @f$ 1/N dN_{ch}/d\eta@f$ in PbPb
   * collisions as a function of @f$\eta@f$ and interaction point
   * Z-coordinate @f$ IP_{z}@f$
   */
  void SetGlobalEmpiricalcorrection(TH2D* h){fEmpiricalCorrection=h;}
  /** 
   * Set re-weighing function to correct for different primary vertex
   * distribution in analysed data from what was used in the empirical
   * correction.  To use, put code like the following into the
   * configuration file (default <tt>dNdetaConfig.C</tt>):
   * @code 
   TF1* f = new TF1("reweight",
                    "TMath::Gaus(x,[0],[1],true)/TMath::Gaus(x,[2],[3],true",
                    -10,10,4);
   f->SetParNames("#mu_{1},#sigma_{1},#mu_{2},#sigma_{2}");
   f->SetParameters(0.592,6.836,muIpz,sigmaIpz);
   f->SetParErrors(0.023,0.029,eMuIpz,eSigmaIpz);
   task->SetIpzReweight(f);
   @endcode 
   *
   * where @a muIpz and @a sigmaIpz are the mean and sigma of the data
   * being analysed.  The values 0.592 (+/- 0.023) and 6.836 (+/-
   * 0.029) are the best-fit (@f$ \chi^2/\nu=0.76@f$) values from the
   * empirical base run (138190).  Currently, the errors on the
   * parameters are not used.
   * 
   * @param f Function to use 
   */
  void SetIpzReweight(TF1* f) { fIpzReweight = f; }
protected:
  /** 
   * Copy contructor - not defined
   */
  AliBasedNdetaTask(const AliBasedNdetaTask&){;}
  /** 
   * Assignment operator - not defined
   * 
   * 
   * @return 
   */
  AliBasedNdetaTask& operator=(const AliBasedNdetaTask&);
  // Forward declaration 
  class CentralityBin;
  /** 
   * Check if the event corresponds to the selected trigger(s),
   * vertex, and centrality.  Derived classes can overload this to
   * enable event processing - even if the event is not within cuts.
   * 
   * @param forward Forward object
   * 
   * @return true if the event is within the cuts. 
   */
  virtual Bool_t CheckEvent(const AliAODForwardMult& forward);
  /** 
   * Create the CentralityBin objects if not already done.
   * 
   */
  virtual void InitializeCentBins();
  /** 
   * Retrieve the histogram 
   * 
   * @param aod AOD event 
   * @param mc  Whether to get the MC histogram or not
   * 
   * @return Retrieved histogram or null
   */
  virtual TH2D* GetHistogram(const AliAODEvent& aod, Bool_t mc=false) = 0;
  /** 
   * Get the colour to use for markers (only pp - in PbPb we use a rainbow)
   * 
   * @return Marker colour 
   */
  virtual Int_t GetColor() const { return kBlack; }
  /** 
   * Get the marker style 
   * 
   * @return Marker style 
   */
  virtual Int_t GetMarker() const { return GetMarkerStyle(kCircle); }
  /** 
   * Massage data histograms if needed
   * 
   * @param vtx 
   * @param data 
   * @param mcData 
   */
  virtual void CheckEventData(Double_t vtx, 
			      TH2*     data, 
			      TH2*     mcData);
  /** 
   * Add a centrality bin 
   * 
   * @param at   Where in the list to add this bin 
   * @param low  Low cut
   * @param high High cut
   */
  void AddCentralityBin(UShort_t at, Float_t low, Float_t high);
  /** 
   * Make a centrality bin 
   * 
   * @param name  Name used for histograms
   * @param low   Low cut in percent
   * @param high  High cut in percent
   * 
   * @return A newly created centrality bin 
   */
  virtual CentralityBin* MakeCentralityBin(const char* name, Float_t low, 
					   Float_t high) const;
  
  // function which applies empirical correction to the AOD object 
  Bool_t ApplyEmpiricalCorrection(const AliAODForwardMult* aod,TH2D* data);
  /** 
   * Static function toe calculate the asymmetry of the input
   * histogram around 0.
   * 
   * @param h Input histogram
   * 
   * @return Newly created histogram or null
   */
  static TH1* Asymmetry(TH1* h);
  /** 
   * Static member to get a simple ID from a centrality esitmator
   * string.  Used to check the validity of the centrality estimator 
   * 
   * @param meth Method string 
   * 
   * @return Id or negative value in case of invalid estimator 
   */
  static Int_t GetCentMethodID(const TString& meth);
  /** 
   * Get the estimator string from simple ID
   * 
   * @param id Simple ID
   * 
   * @return Estimator string 
   */
  static const char* GetCentMethod(UShort_t id);
  /** 
   * Get the centrality.  The trigger mask of the forward object is
   * not modified
   * 
   * @param event    Our event 
   * @param forward  Our FMD event 
   * @param qual     On return, the quality flag 
   *  
   * @return The centrality percentage 
   */
  Double_t GetCentrality(AliAODEvent&       event,
			 AliAODForwardMult* forward,
			 Int_t&             qual);
  /**
   * Get the centrality.  If the quality is bad, set the corresponding
   * bit on the forward object.
   *
   * @param event    Our event 
   * @param forward  Our FMD event 
   */
  Double_t GetCentrality(AliAODEvent&       event,
			 AliAODForwardMult* forward);
  //==================================================================
  /**
   * Class that holds the sum of the data - possibly split into 0 or
   * non-zero bins 
   * 
   */
  struct Sum : public TNamed
  {
    TH2D* fSum;     // Sum of non-zero events
    TH2D* fSum0;    // Sum of zero events 
    TH1I* fEvents;  // Distribution of events 
    Int_t fDebug;   // Debug level
    /** 
     * I/O Constructor - do not use
     */    
    Sum() : fSum(0), fSum0(0), fEvents(0), fDebug(0) {}
    /** 
     * Constructor 
     * 
     * @param name      Name
     * @param postfix   Possible post-fix 
     */
    Sum(const char* name, const char* postfix) 
      : TNamed(name,postfix), 
	fSum(0), 
	fSum0(0), 
	fEvents(0), 
	fDebug(0) 
    {}
    /** 
     * Copy constructor
     * 
     * @param o Object to copy from 
     */
    Sum(const Sum& o) 
      : TNamed(o), 
	fSum(o.fSum), 
	fSum0(o.fSum0), 
	fEvents(o.fEvents), 
	fDebug(o.fDebug) 
    {}
    /** 
     * Assignment operator 
     * 
     * @param o Object to assign from 
     * 
     * @return Reference to this object 
     */
    Sum& operator=(const Sum& o) 
    {
      if (&o == this) return *this;
      SetName(o.GetName()); fSum = o.fSum; fSum0 = o.fSum0; fEvents=o.fEvents;
      return *this;
    }
    /** 
     * Destructor 
     */
    ~Sum() {}
    /** 
     * Initialise this object.  
     * 
     * @param list  List to add histograms to
     * @param data  Format of data to be cloned here
     * @param col   Color 
     */
    void Init(TList* list, const TH2D* data, Int_t col);
    /** 
     * Add an event 
     * 
     * @param data    Data to add
     * @param isZero  If this is zero event
     * @param weight  Event weight 
     */
    void Add(const TH2D* data, Bool_t isZero, Double_t weight);
    /** 
     * Get the histogram name 
     * 
     * @param name Base name 
     * @param what Which one 
     * @param post Possible postfix
     * 
     * @return Name 
     */
    static TString GetHistName(const char* name, Int_t what=0, 
			       const char* post=0);
    /** 
     * Get the histogram name 
     * 
     * @param what Which one 
     * 
     * @return Name 
     */
    TString GetHistName(Int_t what=0) const;
    /** 
     * Get the sum 
     * 
     * @param o          Output list
     * @param ntotal     On return, the total number of events
     * @param zeroEff    Zero-bin efficiency
     * @param otherEff   Non-zero-bin efficiency 
     * @param marker     Marker to use 
     * @param rootXproj  Whether to use TH2::ProjectionX
     * @param corrEmpty  Correct for empty bins 
     * 
     * @return The total sum histogram 
     */
    TH2D* CalcSum(TList* o, Double_t& ntotal,
		  Double_t zeroEff, Double_t otherEff=1, Int_t marker=20,
		  Bool_t rootXproj=false, Bool_t corrEmpty=true) const;

    /** 
     * Print this sum container 
     * 
     * @param option Not used
     */
    virtual void Print(Option_t* option="") const;
    // ClassDef(Sum,2); // Summed histograms
  };
    
  //==================================================================
  /**
   * Calculations done per centrality.  These objects are only used
   * internally and are never streamed.  We do not make dictionaries
   * for this (and derived) classes as they are constructed on the
   * fly.
   * 
   */
  class CentralityBin : public TNamed
  {
  public:
    /** 
     * Constructor 
     */
    CentralityBin();
    /** 
     * Constructor 
     * 
     * @param name Name used for histograms (e.g., Forward)
     * @param low  Lower centrality cut in percent 
     * @param high Upper centrality cut in percent 
     */
    CentralityBin(const char* name, Float_t low, Float_t high);
    /** 
     * Copy constructor 
     * 
     * @param other Object to copy from 
     */
    CentralityBin(const CentralityBin& other);
    /** 
     * Destructor 
     */
    virtual ~CentralityBin();
    /** 
     * Assignment operator 
     * 
     * @param other Object to assign from 
     * 
     * @return Reference to this 
     */
    CentralityBin& operator=(const CentralityBin& other);
    /** 
     * Check if this is the 'all' bin 
     * 
     * @return true if low and high cuts are both zero
     */    
    Bool_t IsAllBin() const { return fLow == 0 && fHigh == 0; }
    Bool_t IsInclusiveBin() const { return fLow == 0 && fHigh == 101; }
    /** 
     * Get the list name 
     * 
     * @return List Name 
     */
    const char* GetListName() const;
    /** 
     * Create output objects 
     * 
     * @param dir   Parent list
     * @param mask  Trigger mask
     */
    virtual void CreateOutputObjects(TList* dir, Int_t mask);
    /** 
     * Process an event
     * 
     * @param forward     Forward data (for trigger, vertex, & centrality)
     * @param triggerMask Trigger mask 
     * @param isZero      True if this is a zero bin event 
     * @param vzMin       Minimum IP z coordinate
     * @param vzMax       Maximum IP z coordinate
     * @param data        Data histogram 
     * @param mc          MC histogram
     * @param filter      Mask of trigger bits to filter out
     * @param weight      Event weight
     *
     * @return true if the event was selected
     */
    virtual Bool_t ProcessEvent(const AliAODForwardMult* forward, 
				UInt_t                   triggerMask,
				Bool_t                   isZero,
				Double_t                 vzMin, 
				Double_t                 vzMax, 
				const TH2D*              data, 
				const TH2D*              mc,
				UInt_t                   filter,
				Double_t                 weight);
    /** 
     * Calculate the Event-Level normalization. 
     * 
     * The full event level normalization for trigger @f$X@f$ is given by 
     * @f{eqnarray*}{
     *    N &=& \frac{1}{\epsilon_X}
     *          \left(N_A+\frac{N_A}{N_V}(N_{-V}-\beta)\right)\\
     *      &=& \frac{1}{\epsilon_X}N_A
     *          \left(1+\frac{1}{N_V}(N_T-N_V-\beta)\right)\\
     *      &=& \frac{1}{\epsilon_X}N_A
     *          \left(1+\frac{N_T}{N_V}-1-\frac{\beta}{N_V}\right)\\
     *      &=& \frac{1}{\epsilon_X}N_A
     *          \left(\frac{1}{\epsilon_V}-\frac{\beta}{N_V}\right)
     * @f}
     * where 
     *
     * - @f$\epsilon_X=\frac{N_{T,X}}{N_X}@f$ is the trigger
     *   efficiency evaluated in simulation.
     * - @f$\epsilon_V=\frac{N_V}{N_T}@f$ is the vertex efficiency 
     *   evaluated from the data 
     * - @f$N_X@f$ is the Monte-Carlo truth number of events of type 
     *   @f$X@f$. 
     * - @f$N_{T,X}@f$ is the Monte-Carlo truth number of events of type 
     *   @f$X@f$ which was also triggered as such. 
     * - @f$N_T@f$ is the number of data events that where triggered 
     *   as type @f$X@f$ and had a collision trigger (CINT1B)
     * - @f$N_V@f$ is the number of data events that where triggered
     *   as type @f$X@f$, had a collision trigger (CINT1B), and had 
     *   a vertex. 
     * - @f$N_{-V}@f$ is the number of data events that where triggered
     *   as type @f$X@f$, had a collision trigger (CINT1B), but no
     *   vertex. 
     * - @f$N_A@f$ is the number of data events that where triggered
     *   as type @f$X@f$, had a collision trigger (CINT1B), and had 
     *   a vertex in the selected range. 
     * - @f$\beta=N_a+N_c-N_e@f$ is the number of control triggers that 
     *   were also triggered as type @f$X@f$. 
     * - @f$N_a@f$ Number of beam-empty events also triggered as type 
     *   @f$X@f$ events (CINT1-A or CINT1-AC). 
     * - @f$N_c@f$ Number of empty-beam events also triggered as type 
     *   @f$X@f$ events (CINT1-C). 
     * - @f$N_e@f$ Number of empty-empty events also triggered as type 
     *   @f$X@f$ events (CINT1-E). 
     * 
     * Note, that if @f$ \beta \ll N_A@f$ the last term can be ignored, and 
     * the expression simplyfies to  
     * @f[
     *  N = \frac{1}{\epsilon_X}\frac{1}{\epsilon_V}N_A
     * @f]
     *
     * @param t       Histogram of triggers 
     * @param scheme  Normalisation scheme 
     * @param trgEff  Trigger efficiency 
     * @param ntotal  On return, the total number of events to normalise to.
     * @param text    If non-null, fill with normalization calculation
     * 
     * @return @f$N_A/N@f$ or negative number in case of errors. 
     */
    virtual Double_t Normalization(const TH1I& t, 
				   UShort_t    scheme,
				   Double_t    trgEff,
				   Double_t&   ntotal,
				   TString*    text) const;
    /** 
     * Generate the dN/deta result from input 
     * 
     * @param sum        Sum of 2D hists 
     * @param postfix    Post fix on names
     * @param rootProj   Whether to use ROOT TH2::ProjectionX
     * @param corrEmpty  Correct for empty bins 
     * @param scaler     Event-level normalization scaler  
     * @param marker     Marker style 
     * @param color       Color of markers 
     * @param mclist      List of MC data 
     * @param truthlist   List of MC truth data 
     */
    virtual void MakeResult(const TH2D* sum,  
			    const char* postfix, 
			    bool        rootProj, 
			    bool        corrEmpty,
			    Double_t    scaler,
			    Int_t       marker,
			    Int_t       color, 
			    TList*      mclist,
			    TList*      truthlist);
    /** 
     * End of processing 
     * 
     * @param sums        List of sums
     * @param results     Output list of results
     * @param scheme      Normalisation scheme options
     * @param trigEff     Trigger efficiency 
     * @param trigEff0    0-bin trigger efficiency 
     * @param rootProj    If true, use TH2::ProjectionX
     * @param corrEmpty   Whether to correct for empty bins
     * @param triggerMask Trigger mask 
     * @param marker      Marker style 
     * @param color       Color of markers 
     * @param mclist      List of MC data 
     * @param truthlist   List of MC truth data 
     *
     * @return true on success
     */
    virtual bool End(TList*      sums, 
		     TList*      results,
		     UShort_t    scheme,
		     Double_t    trigEff,
		     Double_t    trigEff0,
		     Bool_t      rootProj,
		     Bool_t      corrEmpty, 
		     Int_t       triggerMask,
		     Int_t       marker,
		     Int_t       color,
		     TList*      mclist,
		     TList*      truthlist);
    /**
     * @{
     * @name Access histograms
     */
    /** 
     * Get sum histogram 
     * 
     * @param mc If true, return MC histogram 
     * 
     * @return Sum histogram
     */
    const Sum* GetSum(Bool_t mc=false) const { return mc ? fSumMC : fSum; }
    /** 
     * Get sum histogram 
     * 
     * @param mc If true, return MC histogram 
     * 
     * @return Sum histogram
     */
    Sum* GetSum(Bool_t mc=false) { return mc ? fSumMC : fSum; }
    /** 
     * Get trigger histogram
     * 
     * @return Trigger histogram
     */
    const TH1I* GetTriggers() const { return fTriggers; } 
    /** 
     * Get trigger histogram
     * 
     * @return Trigger histogram 
     */
    TH1I* GetTriggers() { return fTriggers; }
    /** 
     * Get trigger histogram
     * 
     * @return Trigger histogram
     */
    const TH1I* GetStatus() const { return fStatus; } 
    /** 
     * Get trigger histogram
     * 
     * @return Trigger histogram 
     */
    TH1I* GetStatus() { return fStatus; }
    /** @} */

    /** 
     * Get the color of the markers
     *
     * @param fallback Fall-back color 
     *
     * @return Color for this centrality bin 
     */
    Int_t GetColor(Int_t fallback=kRed+2) const;
    void SetColor(Color_t colour) { fColor = colour; }
    /** 
     * Get list of results 
     * 
     * @return List of results
     */
    TList* GetResults() const { return fOutput; }
    /** 
     * Get name of result histogram. Note, the returned pointer points
     * to static memory and should be copied/used immediately.
     * 
     * @param postfix  Possible postfix (e.g., "MC")
     * 
     * @return 
     */
    const char* GetResultName(const char* postfix="") const;
    /** 
     * Get a result 
     * 
     * @param postfix  Possible postfix (e.g., "MC")
     * @param verbose  If true, complain about missing histogram
     * 
     * @return Pointer to histogram or null
     */
    TH1* GetResult(const char* postfix="",
		   Bool_t      verbose=true) const;
    /** 
     * Set the debug level
     * 
     * @param lvl Debug level
     */
    void SetDebugLevel(Int_t lvl);
    /**
     * Set satellite vertex flag
     *
     * @param satVtx
     */
    void SetSatelliteVertices(Bool_t satVtx) { fSatelliteVertices = satVtx; }
    /** 
     * Print this centrality bin. 
     * 
     * @param option Options. 
     *
     * - R  Also print sums 
     */
    virtual void Print(Option_t* option="") const;
  protected:
    /** 
     * Read in sum hisotgram from list 
     * 
     * @param list List to read from 
     * @param mc   True for MC input 
     * 
     * @return true if sum histogram is found
     */
    virtual Bool_t ReadSum(TList* list, bool mc=false);
    /** 
     * Create sum histogram 
     * 
     * @param data  Data histogram to clone 
     * @param mc    (optional) MC histogram to clone 
     */
    virtual void CreateSums(const TH2D* data, const TH2D* mc);
    /** 
     * Check the trigger, vertex, and centrality
     * 
     * @param forward      Event input 
     * @param triggerMask  The used trigger mask 
     * @param vzMin        Least @f$ v_z@f$
     * @param vzMax        Largest @f$ v_z@f$
     * @param filter       Mask of trigger bits to filter out
     * 
     * @return true if the event is to be used 
     */
    virtual Bool_t CheckEvent(const AliAODForwardMult* forward, 
			      Int_t                    triggerMask,
			      Double_t                 vzMin, 
			      Double_t                 vzMax,
			      Int_t                    filter);
    TList*   fSums;      // Output list 
    TList*   fOutput;    // Output list 
    Sum*     fSum;       // Sum histogram
    Sum*     fSumMC;     // MC sum histogram
    TH1I*    fTriggers;  // Trigger histogram 
    TH1I*    fStatus;    // Trigger histogram 
    Float_t  fLow;       // Lower limit (inclusive)
    Float_t  fHigh;      // Upper limit (exclusive)
    Bool_t   fDoFinalMCCorrection; //Do final MC correction
    Bool_t   fSatelliteVertices; // Satellite vertex flag
    Int_t    fDebug;     // Debug level
    Color_t  fColor;     // Colour used 

    // ClassDef(CentralityBin,5); // A centrality bin 
  };
  Bool_t          fCorrEmpty;    // Correct for empty bins 
  Bool_t          fUseROOTProj;  // Whether to use ROOT's ProjectionX
  Double_t        fTriggerEff;   // Trigger efficiency for selected trigger(s)
  Double_t        fTriggerEff0;  // Bin-0 Trigger efficiency for sel trigger(s)
  TObjArray*      fListOfCentralities; // Centrality bins 
  UShort_t        fNormalizationScheme; // Normalization scheme
  TString         fFinalMCCorrFile; //Filename for final MC corr
  Bool_t          fSatelliteVertices; // satellite vertex flag
  TH2D*           fEmpiricalCorrection; // Empirical correction 
  TProfile* 	  fMeanVsC;         //mean signal per event vs cent
  TH1D*           fSeenCent;      // The seen centrality 
  TH1D*           fTakenCent;     // The taken centrality 
  TString         fCentMethod;    // Centrality estimator 
  AliAnalysisUtils fAnaUtil;      // Analysis utility 
  Bool_t          fUseUtilPileup; // Check for SPD outliers
  TF1*            fIpzReweight;   // Re-weighing function
  Double_t        fCacheCent;     // Stored centrality
  Int_t           fCacheQual;     // Stored centrality quality flag 
  ClassDef(AliBasedNdetaTask,19); // Determine charged particle density
};

#endif
//
// Local Variables:
//  mode: C++
// End:
//
