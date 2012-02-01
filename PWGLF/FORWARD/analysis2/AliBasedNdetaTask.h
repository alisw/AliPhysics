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
 * @ingroup pwg2_forward_dndeta
 * 
 */
#include <AliAnalysisTaskSE.h>
class TAxis;
class TList;
class TH2D;
class TH1D;
class TH1I;
class AliAODEvent;
class AliAODForwardMult;
class TObjArray;

/** 
 * @defgroup pwg2_forward_tasks_dndeta dN/deta tasks 
 * @ingroup pwg2_forward_tasks 
 */
/**
 * @defgroup pwg2_forward_dndeta dN/deta
 *
 * @f$ dN/d\eta@f$ code 
 *
 * @ingroup pwg2_forward_topical
 */
/**
 * Base class for tasks to determine @f$ dN/d\eta@f$ 
 *
 * @ingroup pwg2_forward_tasks_dndeta
 * @ingroup pwg2_forward_dndeta
 */
class AliBasedNdetaTask : public AliAnalysisTaskSE
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
     *   N = \frac{1}{\epsilon_X}(N_A-N_A/N_V(N_T-N_V)) = 
     *       \frac{1}{\epsilon_X}\frac{1}{\epsilon_V}N_A
     * @f]
     */
    kEventLevel = 0x1,
    /** 
     * Do the shape correction
     */
    kShape = 0x2, 
    /** 
     * Correct for background events (A+C-E). Not implemented yet
     */
    kBackground = 0x4,
    /**
     * Correct for the trigger efficiency from MC 
     */
    kTriggerEfficiency = 0x8,
    /** 
     * Correct using zero-bin efficiency only 
     */
    kZeroBin = 0x10,
    /**
     * Do the full correction
     */
    kFull = kEventLevel | kShape | kBackground | kTriggerEfficiency,
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
   * Set the vertex range to use 
   * 
   * @param min Minimum (in centermeter)
   * @param max Maximum (in centermeter)
   */  
  void SetVertexRange(Double_t min, Double_t max) { fVtxMin=min; fVtxMax=max; }
  /** 
   * Set the rebinning factor 
   * 
   * @param rebin Rebinning factor 
   */
  void SetRebinning(Int_t rebin) { fRebin = rebin; }
  /** 
   * Set the trigger maskl 
   * 
   * @param mask Trigger mask
   */
  void SetTriggerMask(UShort_t mask);
  /** 
   * Set the trigger mask 
   * 
   * @param mask trigger mask 
   */
  void SetTriggerMask(const char* mask);
  /** 
   * Set the centrality bins to use. 
   * 
   * @code 
   *   UShort_t bins[] = { 0, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100 };
   *   task->SetCentralityBins(11, bins);
   * @endcode 
   * 
   * @param n     Number of bins (elements in @a bins minus 1)
   * @param bins  Bin limits 
   */
  void SetCentralityAxis(UShort_t n, Short_t* bins);
  /** 
   * Whether to cut edges when merging 
   * 
   * @param cut If true, cut edges 
   */
  void SetCutEdges(Bool_t cut) {fCutEdges = cut;}
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
  /** 
   * Trigger efficiency for selected trigger(s)
   * 
   * @param e Trigger efficiency 
   */
  void SetTriggerEff(Double_t e) { fTriggerEff = e; } 
  /** 
   * Set the shape correction (a.k.a., track correction) for selected
   * trigger(s)
   * 
   * @param h Correction
   */
  void SetShapeCorrection(const TH1* h);
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
  /** 
   * Load the normalization data - done automatically if not set from outside
   * 
   * @param sys system
   * @param energy energy
   */
  void LoadNormalizationData(UShort_t sys, UShort_t energy);  
  /** @} */
  /** 
   * Print information 
   * 
   * @param option Not used
   */
  void Print(Option_t* option="") const;
  /** @{ 
   *  @name Task interface 
   */
  /** 
   * Initialise on master - does nothing
   * 
   */
  virtual void   Init() {}
  /** 
   * Create output objects.  
   *
   * This is called once per slave process 
   */
  virtual void UserCreateOutputObjects();
  /** 
   * Process a single event 
   * 
   * @param option Not used
   */
  virtual void UserExec(Option_t* option);
  /** 
   * Called at end of event processing.
   *
   * This is called once in the master 
   * 
   * @param option Not used 
   */
  virtual void Terminate(Option_t* option);
  /* @} */

  /** 
   * @{ 
   * @name Services member functions 
   */
  /** 
   * Make a copy of the input histogram and rebin that histogram
   * 
   * @param h         Histogram to rebin
   * @param rebin     Rebinning factor 
   * @param cutEdges  Whether to cut edges when rebinning
   * 
   * @return New (rebinned) histogram
   */
  static TH1D* Rebin(const TH1D* h, Int_t rebin, Bool_t cutEdges=false);
  /** 
   * Make an extension of @a h to make it symmetric about 0 
   * 
   * @param h Histogram to symmertrice 
   * 
   * @return Symmetric extension of @a h 
   */
  static TH1* Symmetrice(const TH1* h);
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
				     const char* ytitle="#frac{1}{N} #frac{dN_{ch}}{d#eta}");
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
  static Int_t GetMarkerStyle(UShort_t bits);
  static UShort_t GetMarkerBits(Int_t style);
  static Int_t FlipHollowStyle(Int_t style);
protected:
  /** 
   * Copy contructor
   */
  AliBasedNdetaTask(const AliBasedNdetaTask&);
  /** 
   * Assignment operator 
   * 
   * 
   * @return 
   */
  AliBasedNdetaTask& operator=(const AliBasedNdetaTask&) { return *this; }
  // Forward declaration 
  class CentralityBin;

  /** 
   * Retrieve the histogram 
   * 
   * @param aod AOD event 
   * @param mc  Whether to get the MC histogram or not
   * 
   * @return Retrieved histogram or null
   */
  virtual TH2D* GetHistogram(const AliAODEvent* aod, Bool_t mc=false) = 0;
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
   * Add a centrality bin 
   * 
   * @param at   Where in the list to add this bin 
   * @param low  Low cut
   * @param high High cut
   */
  void AddCentralityBin(UShort_t at, Short_t low, Short_t high);
  /** 
   * Make a centrality bin 
   * 
   * @param name  Name used for histograms
   * @param low   Low cut in percent
   * @param high  High cut in percent
   * 
   * @return A newly created centrality bin 
   */
  virtual CentralityBin* MakeCentralityBin(const char* name, Short_t low, 
					   Short_t high) const;
  
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
    /** 
     * I/O Constructor - do not use
     */    
    Sum() : fSum(0), fSum0(0), fEvents(0) {}
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
	fEvents(0) 
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
	fEvents(o.fEvents) 
    {}
    /** 
     * Assignment operator 
     * 
     * @param o Object to assign from 
     * 
     * @return Reference to this object 
     */
    Sum& operator=(const Sum& o) {
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
     */
    void Add(const TH2D* data, Bool_t isZero=false);
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
     * @param l          List to get histograms from 
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
    TH2D* GetSum(const TList* l, TList* o, Double_t& ntotal,
		 Double_t zeroEff, Double_t otherEff=1, Int_t marker=20,
		 Bool_t rootXproj=false, Bool_t corrEmpty=true) const;
  };
    
  //==================================================================
  /**
   * Calculations done per centrality 
   * 
   */
  class CentralityBin : public TNamed
  {
  public:
    /** dN
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
    CentralityBin(const char* name, Short_t low, Short_t high);
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
     */
    virtual void CreateOutputObjects(TList* dir);
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
     */
    virtual void ProcessEvent(const AliAODForwardMult* forward, 
			      Int_t                    triggerMask,
			      Bool_t                   isZero,
			      Double_t                 vzMin, 
			      Double_t                 vzMax, 
			      const TH2D*              data, 
			      const TH2D*              mc);
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
     * 
     * @return @f$N_A/N@f$ or negative number in case of errors. 
     */
    virtual Double_t Normalization(const TH1I& t, 
				   UShort_t    scheme,
				   Double_t    trgEff,
				   Double_t&   ntotal) const;
    /** 
     * Generate the dN/deta result from input 
     * 
     * @param sum        Sum of 2D hists 
     * @param postfix    Post fix on names
     * @param rootProj   Whether to use ROOT TH2::ProjectionX
     * @param corrEmpty  Correct for empty bins 
     * @param shapeCorr  Shape correction to use 
     * @param scaler     Event-level normalization scaler  
     * @param symmetrice Whether to make symmetric extensions 
     * @param rebin      Whether to rebin
     * @param cutEdges   Whether to cut edges when rebinning 
     * @param marker     Marker style 
     * @param color       Color of markers 
     * @param mclist      List of MC data 
     * @param truthlist   List of MC truth data 
     */
    virtual void MakeResult(const TH2D* sum,  
			    const char* postfix, 
			    bool        rootProj, 
			    bool        corrEmpty,
			    const TH1*  shapeCorr,
			    Double_t    scaler,
			    bool        symmetrice, 
			    Int_t       rebin, 
			    bool        cutEdges, 
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
     * @param shapeCorr   Shape correction or nil
     * @param trigEff     Trigger efficiency 
     * @param symmetrice  Whether to symmetrice the results
     * @param rebin       Whether to rebin the results
     * @param rootProj    If true, use TH2::ProjectionX
     * @param corrEmpty   Whether to correct for empty bins
     * @param cutEdges    Whether to cut edges when rebinning
     * @param triggerMask Trigger mask 
     * @param marker      Marker style 
     * @param color       Color of markers 
     * @param mclist      List of MC data 
     * @param truthlist   List of MC truth data 
     */
    virtual void End(TList*      sums, 
		     TList*      results,
		     UShort_t    scheme,
		     const TH1*  shapeCorr, 
		     Double_t    trigEff,
		     Bool_t      symmetrice,
		     Int_t       rebin, 
		     Bool_t      rootProj,
		     Bool_t      corrEmpty, 
		     Bool_t      cutEdges, 
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
    TH1I* GetTrigggers() { return fTriggers; }
    /** @} */

    /** 
     * Get the color of the markers
     *
     * @return 
     */
    Int_t GetColor(Int_t fallback=kRed+2) const;
    /** 
     * Get list of results 
     * 
     * 
     * @return 
     */
    TList* GetResults() const { return fOutput; }
    /** 
     * Get name of result histogram 
     * 
     * @param rebin 
     * @param sym 
     * @param postfix 
     * 
     * @return 
     */
    const char* GetResultName(Int_t rebin, Bool_t sym, 
			      const char* postfix="") const;
    /** 
     * Get a result 
     * 
     * @param rebin 
     * @param sym 
     * @param postfix 
     * 
     * @return 
     */
    TH1* GetResult(Int_t rebin, Bool_t sym, 
		   const char* postfix="") const;

  protected:
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
     * @param forward Event input 
     * @param triggerMask  The used trigger mask 
     * @param vzMin        Least @f$ v_z@f$
     * @param vzMax        Largest @f$ v_z@f$
     * 
     * @return true if the event is to be used 
     */
    virtual Bool_t CheckEvent(const AliAODForwardMult* forward, 
			      Int_t                    triggerMask,
			      Double_t                 vzMin, 
			      Double_t vzMax);
    TList*   fSums;      // Output list 
    TList*   fOutput;    // Output list 
    Sum*     fSum;       // Sum histogram
    Sum*     fSumMC;     // MC sum histogram
    TH1I*    fTriggers;  // Trigger histogram 
    UShort_t fLow;       // Lower limit (inclusive)
    UShort_t fHigh;      // Upper limit (exclusive)
    Bool_t   fDoFinalMCCorrection; //Do final MC correction
    
    ClassDef(CentralityBin,1); // A centrality bin 
  };
  TList*          fSums;         // Container of sums 
  TList*          fOutput;       // Container of outputs 
  Double_t        fVtxMin;       // Minimum v_z
  Double_t        fVtxMax;       // Maximum v_z
  Int_t           fTriggerMask;  // Trigger mask 
  Int_t           fRebin;        // Rebinning factor 
  Bool_t          fCutEdges;     // Whether to cut edges when rebinning
  Bool_t          fSymmetrice;   // Whether to symmetrice data 
  Bool_t          fCorrEmpty;    // Correct for empty bins 
  Bool_t          fUseROOTProj;  // Whether to use ROOT's ProjectionX
  Double_t        fTriggerEff;   // Trigger efficiency for selected trigger(s)
  TH1*            fShapeCorr;    // Shape correction 
  TObjArray*      fListOfCentralities; // Centrality bins 
  TNamed*         fSNNString;    // sqrt(s_NN) string 
  TNamed*         fSysString;    // Collision system string 
  TH1D*           fCent;         // Centrality distribution 
  TAxis*          fCentAxis;     // Centrality axis
  UShort_t        fNormalizationScheme; // Normalization scheme
  TNamed*         fSchemeString;     
  TNamed*         fTriggerString; 
  TString         fFinalMCCorrFile; //Filename for final MC corr
  
  ClassDef(AliBasedNdetaTask,4); // Determine multiplicity in base area
};

#endif
//
// Local Variables:
//  mode: C++
// End:
//
