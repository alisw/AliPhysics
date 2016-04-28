/**
 * @file   AliTrackletAODdNdeta.C
 * @author Christian Holm Christensen <cholm@nbi.dk>
 * @date   Wed Apr 27 16:46:04 2016
 * 
 * @brief  AOD tasks to do final dN/deta in midrapidity
 * 
 * @ingroup pwglf_forward_tracklets
 * 
 */

#include <AliAnalysisTaskSE.h>
#include <AliTrackletAODUtils.C>
#include <AliAODTracklet.C>
#ifndef __CINT__
#include <cctype>
#include "AliAODSimpleHeader.C"
#include <AliVVertex.h>
#include <AliVertex.h>
#include <AliAnalysisManager.h>
#include <AliVEventHandler.h>
#include <AliInputEventHandler.h>
#include <AliMultSelection.h>
#include <AliCentrality.h>
#include <AliLog.h>
#include <TClonesArray.h>
#include "AliTrackletWeights.C"
#include <TUrl.h>
#include <TFile.h>
#else
class AliAODTracklet;
class AliTrackletWeights;
class AliVEvent;
class AliMultSelection;  // Auto-load 
class TClonesArray;
#endif

//====================================================================
/**
 * Task to analyse AOD tracklets for dNch/deta 
 * 
 */
class AliTrackletAODdNdeta : public AliAnalysisTaskSE,
			     public AliTrackletAODUtils
{
public:
  /** Type of containers */
  typedef TList Container;
  /** 
   * Default constructor - for ROOT I/O only 
   */
  AliTrackletAODdNdeta();
  /** 
   * Named - user - constructor
   */
  AliTrackletAODdNdeta(const char* name);
  /**
   * Copy constructor 
   *
   * @param o Object to copy from 
   */
  AliTrackletAODdNdeta(const AliTrackletAODdNdeta& o);
  /**
   * Destructor 
   */
  virtual ~AliTrackletAODdNdeta() {}
  /** 
   * Assignment operator 
   * 
   * @param o Object to assign from 
   * 
   * @return Reference to this object 
   */
  AliTrackletAODdNdeta& operator=(const AliTrackletAODdNdeta& o);
  /** 
   * @{ 
   * @name Auxiliary methods 
   */
  /** 
   * Print information to standard output 
   * 
   * @param option Ignored 
   */
  virtual void Print(Option_t* option="") const;
  /** 
   * Connect this task to analysis manager 
   * 
   * @param sumFile (Optional) name of sum file 
   * @param resFile (Optional) name of result file 
   */
  Bool_t Connect(const char* sumFile=0,const char* resFile=0);
  /** 
   * Create object of this (or derived) class 
   * 
   * @param mc      Whether this is for MC or not
   * @param sumFile (Optional) name of sum file 
   * @param resFile (Optional) name of result file 
   * @param weights (Optional) name of weights file 
   * 
   * @return Newly allocated task or null
   */
  static AliTrackletAODdNdeta* Create(Bool_t      mc=false,
				      const char* weights=0,
				      const char* sumFile=0,
				      const char* resFile=0);
  /* @} */

  // -----------------------------------------------------------------
  /** 
   * @{ 
   * @name Task interface 
   */
  /**
   * Delegate worker initialization 
   * 
   */
  void UserCreateOutputObjects();
  /**
   * Event processing 
   */
  void UserExec(Option_t*); 
  /** 
   * Called at end of worker job. 
   * 
   */
  void FinishTaskOutput() { /*WorkerFinalize();*/ }
  /** 
   * Called at end of master job on merged results. 
   * 
   */
  void Terminate(Option_t*);
  /* @} */

  // -----------------------------------------------------------------
  /** 
   * @{ 
   * @name Set parameters on the task 
   */
  /** 
   * Set the centrality method to use 
   * 
   * @param name Name of centrality method 
   */
  void SetCentralityMethod(const TString& name) { fCentMethod = name; }
  /** 
   * Set the centrality axis 
   * 
   * @param n    Number of bins 
   * @param bins Bin borders (n+1 entries)
   */
  void SetCentralityAxis(Int_t n, Double_t* bins)
  {
    SetAxis(fCentAxis,n,bins);
  }
  /** 
   * Set the centraliy axis 
   * 
   * @param spec bin specification  
   */
  void SetCentralityAxis(const TString& spec)
  {
    SetAxis(fCentAxis, spec, "-:,");
  }
  /** 
   * Set the @f$ \mathrm{IP}_z@f$ axis 
   * 
   * @param n   Number of bins
   * @param min Least value 
   * @param max Largest value 
   */
  void SetIPzAxis(Int_t n, Double_t min, Double_t max)
  {
    SetAxis(fIPzAxis, n, min, max);
  }
  /** 
   * Set the @f$ \mathrm{IP}_z@f$ axis 
   * 
   * @param n   Number of bins
   * @param max Largest absolute value 
   */
  void SetIPzAxis(Int_t n, Double_t max)
  {
    SetAxis(fIPzAxis, n, max);
  }
  /** 
   * Set the interaction point Z axis 
   * 
   * @param spec bin specification  
   */
  void SetIPzAxis(const TString& spec)
  {
    SetAxis(fIPzAxis, spec);
  }
  /** 
   * Set the @f$ \varphi@f$ axis 
   * 
   * @param n   Number of bins
   * @param min Least value 
   * @param max Largest value 
   */
  void SetPhiAxis(Int_t n, Double_t min, Double_t max)
  {
    SetAxis(fPhiAxis, n, min, max);
  }
  /** 
   * Set the @f$ \varphi@f$ axis 
   * 
   * @param n   Number of bins
   * @param max Largest absolute value 
   */
  void SetPhiAxis(Int_t n, Double_t max)
  {
    SetAxis(fPhiAxis, n, max);
  }
  /** 
   * Set the @f$ \vareta@f$ axis 
   * 
   * @param n   Number of bins
   * @param min Least value 
   * @param max Largest value 
   */
  void SetEtaAxis(Int_t n, Double_t min, Double_t max)
  {
    SetAxis(fEtaAxis, n, min, max);
  }
  /** 
   * Set the @f$ \vareta@f$ axis 
   * 
   * @param n   Number of bins
   * @param max Largest absolute value 
   */
  void SetEtaAxis(Int_t n, Double_t max)
  {
    SetAxis(fEtaAxis, n, max);
  }
  /** 
   * Set the pseudorapidity axis 
   * 
   * @param spec bin specification  
   */
  void SetEtaAxis(const TString& spec)
  {
    SetAxis(fEtaAxis, spec);
  }
  /**
   * Set @f$\delta_{\phi}@f$
   *
   * @param x Shift of @f$\Delta\phi@f$ 
   */
  void SetDPhiShift(Double_t x=0.0045) { fDPhiShift = x; }
  /**
   * Set cut on @f$\Delta\phi-\delta_{\phi}@f$ 
   *
   * @param x Cut on shifted @f$\Delta\phi@f$ 
   */
  void SetShiftedDPhiCut(Double_t x=-1) { fShiftedDPhiCut = x; }
  /**
    * Set upper cut on @f$\Delta@f$ for signals 
    *
    * @param x Value 
    */
  void SetDeltaCut(Double_t x=1.5) { fDeltaCut = x; }
  /**
   * Set Maximum @f$ \Delta@f$ to consider 
   *
   * @param x Value 
   */
  void SetMaxDelta(Double_t x=25) { fMaxDelta = x; }
  /**
   * Set lower cut on tail of @f$\Delta@f$ distributions 
   *
   * @param x Value 
   */
  void SetTailDelta(Double_t x=5) { fTailDelta = x; }
  /** 
   * Set the very least centrality to consider.  This is for cases
   * where the centrality calibration of simulated data doesn't really
   * match the real-data one, but we want to keep the original
   * centrality bins.  E.g., for LHC15k1b1, the calibration of the
   * DPMJet centrality does not give the same mean number of tracklets
   * for the very most central bin.  We therefor need to rule out the
   * event with a very large number of tracklets from the sample by
   * setting this parameter to for example 0.1.  If this parameter is
   * set to a negative value (default) it is not considered.
   * 
   * @param x Absolute lowest centrality to consider 
   */
  void SetAbsMinCent(Double_t x=-1) { fAbsMinCent = x; }
  /** 
   * Set weights.  This implementatio does nothing 
   * 
   * @param w Weights to use 
   */
  virtual void SetWeights(AliTrackletWeights* w) {};
  /* @} */
  //__________________________________________________________________
  /**
   * Base class for sub-components 
   */
  struct Sub : public TObject
  {
    /** 
     * Constructor 
     * 
     */
    Sub(const char* name="") : fName(name), fContainer(0), fDebug(0) {}
    /** 
     * Copy constructor 
     * 
     * @param o Object to copy from 
     */
    Sub(const Sub& o) : TObject(o), fName(o.fName), fContainer(0), fDebug(0) {}
    /**
     * Destructor 
     */
    virtual ~Sub() {}
    /** 
     * Assignment operator 
     *
     * @return reference to this 
     */
    Sub& operator=(const Sub&) { return *this; }
    /** 
     * @return The name 
     */
    const char* GetName() const { return fName.Data(); }
    /** 
     * Initialize the bin 
     * 
     * @param parent    Parent container 
     * @param etaAxis   pseudorapidity axis to use 
     * @param ipzAxis   Interaction point Z coordinate axis 
     * @param deltaAxis @f$\Delta@f$ axis to use 
     * 
     * @return true on success 
     */
    virtual Bool_t WorkerInit(Container*   parent,
			      const TAxis& etaAxis,
			      const TAxis& ipzAxis,
			      const TAxis& deltaAxis)
    {
      fContainer = new Container;
      fContainer->SetName(fName);
      fContainer->SetOwner();
      if (parent) parent->Add(fContainer);
      return true;
    }
    /** 
     * Process a single tracklet 
     * 
     * @param tracklet The tracklet 
     * @param ipz      Z-coordinate of the IP
     * @param signal   True if a signal 
     * @param weight   Weight of tracklet 
     * 
     * @return true on success 
     */
    virtual Bool_t ProcessTracklet(AliAODTracklet* tracklet,
				   Double_t        ipz,
				   UShort_t        signal,
				   Double_t        weight) = 0;
    /** 
     * Initialize this sub-component at the time of finalizing the
     * job.  Should find sum container in @a parent and extract data
     * from that container.
     * 
     * @param parent Parent container of sum data 
     *  
     * @return true on success
     */
    virtual Bool_t FinalizeInit(Container* parent) = 0;
    /** 
     * Called on master when terminating 
     * 
     * @param parent  Parent container 
     * @param ipz     Distribution of interaction point Z coordinate
     * @param tailCut Cut on tails 
     * 
     * @return 
     */
    virtual Bool_t MasterFinalize(Container* parent,
				  TH1*       ipz,
				  Double_t   tailCut) = 0;
    /** 
     * Set debug flag 
     *
     * @param lvl Debug level 
     */
    virtual void SetDebug(UShort_t lvl) { fDebug = lvl; }
  protected:
    /** The name of the sub-component */
    TString fName;
    /** The sum container of the sub-component */
    Container* fContainer;
    /** Debug flag */
    UShort_t fDebug;
    ClassDef(Sub,1); 
  };    
  
  //__________________________________________________________________
  /**
   * A set of histograms. Contains 
   *
   * - @f$\eta,\mathrm{IP}_z@f$ distribution @f$ X_{\eta,z}@f$ 
   *
   * and optionally
   *
   * - @f$\eta,\Delta@f$ distribution @f$ f_X(\Delta)@f$ 
   */
  struct Histos : public Sub
  {
    /** 
     * Constructor 
     * 
     */
    Histos(const char* name="", Color_t color=kBlack, Style_t style=1,
	   UChar_t mask=0, UChar_t veto=0)
      : Sub(name),
	fMask(mask),
	fVeto(veto),
	fEtaIPz(0),
	fEtaDeltaIPz(0),
	fColor(color),
	fStyle(style)
    {}
    /** 
     * Copy constructor 
     * 
     * @param o Object to copy from 
     */
    Histos(const Histos& o)
      : Sub(o),
	fMask(o.fMask),
	fVeto(o.fVeto),
	fEtaIPz(0),
	fEtaDeltaIPz(0)
    {}
    /**
     * Destructor 
     */
    virtual ~Histos() {}
    /** 
     * Assignment operator 
     *
     * @return reference to this 
     */
    Histos& operator=(const Histos&) { return *this; }
    /** 
     * @return Selection corresponds to observed tracklets
     */    
    Bool_t IsMeasured() const
    {
      return fMask == kMeasuredMask && fVeto == kMeasuredVeto;
    }
    /** 
     * @return Selection corresponds to tracklets from injection events
     */
    Bool_t IsInjected() const
    {
      return fMask == kInjectedMask && fVeto == kInjectedVeto;
    }
    /** 
     * @return Selection corresponds to tracklets from two particles
     */
    Bool_t IsCombinatoric() const
    {
      return fMask == kCombinatoricMask && fVeto == kCombinatoricVeto;
    }
    /** 
     * @return Selection corresponds to primaries particles 
     */
    Bool_t IsPrimary() const
    {
      return fMask == kPrimaryMask && fVeto == kPrimaryVeto;
    }
    /** 
     * @return Selection corresponds to secondary particles 
     */
    Bool_t IsSecondary() const
    {
      return fMask == kSecondaryMask && fVeto == kSecondaryVeto;
    }
    /** 
     * @return Selection corresponds to generated particles 
     */
    Bool_t IsGenerated() const
    {
      return fMask == kGeneratedMask && fVeto == kGeneratedVeto;
    }
    /** 
     * Initialize the bin 
     * 
     * @param parent    Parent container 
     * @param etaAxis   pseudorapidity axis to use 
     * @param ipzAxis   Interaction point Z coordinate axis 
     * @param deltaAxis @f$\Delta@f$ axis to use 
     * 
     * @return true on success 
     */
    Bool_t WorkerInit(Container*   parent,
		      const TAxis& etaAxis,
		      const TAxis& ipzAxis,
		      const TAxis& deltaAxis);
    /** 
     * Set Attributes on histograms 
     * 
     * @param c Marker/line/fill color 
     * @param m Marker style 
     */
    void SetAttr(Color_t c, Style_t m);
    /** 
     * Process a single tracklet 
     * 
     * @param tracklet The tracklet 
     * @param ipz      Z-coordinate of the IP
     * @param signal   True if a signal 
     * @param weight   Weight of tracklet 
     * 
     * @return true on success 
     */
    Bool_t ProcessTracklet(AliAODTracklet* tracklet,
			   Double_t        ipz,
			   UShort_t        signal,
			   Double_t        weight);
    /** 
     * Initialize this sub-component at the time of finalizing the
     * job.  Should find sum container in @a parent and extract data
     * from that container.
     * 
     * @param parent Parent container of sum data 
     *  
     * @return true on success
     */
    Bool_t FinalizeInit(Container* parent);
    /** 
     * Called on master when terminating 
     * 
     * @param parent  Parent container 
     * @param ipz     Interaction point Z coordinate distribution
     * @param tailCut Cut on tails 
     * 
     * @return true on success
     */
    Bool_t MasterFinalize(Container* parent,
			  TH1*       ipz,
			  Double_t   tailCut);

    UChar_t GetMask() const { return fMask; }
    UChar_t GetVeto() const { return fVeto; }
    /** 
     * Print information to standard output 
     * 
     * @param option Ignored 
     */
    void Print(Option_t* option="") const;
  protected:
    UChar_t fMask;
    UChar_t fVeto;
    TH2*    fEtaIPz;    //!
    TH3*    fEtaDeltaIPz;  //!
  public:
    Color_t fColor;
    Style_t fStyle;
    
    ClassDef(Histos,1); 
  };    


  //__________________________________________________________________
  /**
   * A centrality bin.  Here, we need 
   *
   * - @f$\eta,\Delta@f$ distribution of measured @f$ f_M(\Delta)@f$ 
   * - @f$\eta,\Delta@f$ distribution of injected @f$ f_I(\Delta)@f$ 
   * - @f$\eta,\mathrm{IP}_z@f$ distribution of measured @f$ M_{\eta,z}@f$ 
   * - @f$\eta,\mathrm{IP}_z@f$ distribution of injected @f$ I_{\eta,z}@f$  
   *
   * For MC we also need 
   *
   * - @f$\eta,\mathrm{IP}_z@f$ distribution of combi @f$ C_{\eta,z}@f$ 
   * - @f$\eta,\mathrm{IP}_z@f$ distribution of gen @f$ P_{\eta,z}@f$ 
   * 
   * 
   */
  struct CentBin : public Sub
  {
    /** 
     * Default constructor - for ROOT I/O only  
     */
    CentBin()
      : Sub(""),
	fSubs(0),
	fLow(0),
	fHigh(0),
	fIPz(0),
	fCent(0),
	fMeasured(0),
	fInjection(0)
    {
    } 
    /** 
     * User constructor 
     * 
     * @param c1  Lower bound on centrality 
     * @param c2  Upper bound on centrality
     */
    CentBin(Double_t c1, Double_t c2);
    /** 
     * Copy constructor 
     * 
     * @param o Object to copy from 
     */
    CentBin(const CentBin& o)
      : Sub(o),
	fSubs(0),
	fLow(o.fLow),
	fIPz(0),
	fCent(0),
	fHigh(o.fHigh),
	fMeasured(0),
	fInjection(0)	
    {}
    /**
     * Destructor 
     */
    virtual ~CentBin() {}
    /** 
     * Assignment operator 
     *
     * @return reference to this 
     */
    CentBin& operator=(const CentBin&) { return *this; }
    /** 
     * Initialize the bin 
     * 
     * @param etaAxis  pseudorapidity axis to use 
     * @param ipzAxis  Interaction point Z coordinate axis 
     * @param deltaMax Largest @f$\Delta@f$ to consider 
     * 
     * @return true on success 
     */
    Bool_t WorkerInit(Container*   parent,
		      const TAxis& etaAxis,
		      const TAxis& ipzAxis,
		      const TAxis& deltaAxis);
    /** 
     * Check if we should process this event 
     * 
     * @param cent Event centrality 
     * @param ipz  Event Z-coordinate of the interaction 
     * 
     * @return true if we should process the event 
     */
    Bool_t Accept(Double_t cent, Double_t ipz);
    /** 
     * Process a single tracklet 
     * 
     * @param tracklet The tracklet 
     * @param ipz      Z-coordinate of the IP
     * @param signal   True if a signal 
     * @param weight   Weight of tracklet 
     * 
     * @return true on success 
     */
    Bool_t ProcessTracklet(AliAODTracklet* tracklet,
			   Double_t        ipz,
			   UShort_t        signal,
			   Double_t        weight);
    /** 
     * Initialize this sub-component at the time of finalizing the
     * job.  Should find sum container in @a parent and extract data
     * from that container.
     * 
     * @param parent Parent container of sum data 
     *  
     * @return true on success
     */
    Bool_t FinalizeInit(Container* parent);
    /** 
     * Called on master when terminating 
     * 
     * @param parent  Parent container 
     * @param ipz     Z-coordinate of the IP
     * @param tailCut Cut on tails 
     * 
     * @return 
     */
    Bool_t MasterFinalize(Container* parent,
			  TH1*       ipz,
			  Double_t   tailCut);
    /** 
     * Estimate the background a given histogram set 
     * 
     * @param result   Output container 
     * @param measCont The measured results 
     * @param genCont  The generator results (if applicable)
     * @param h        The histogram container 
     * @param tailCut  Cut on the tail distribution 
     * 
     * @return true on success 
     */
    Bool_t EstimateBackground(Container* result,
			      Container* measCont,
			      Container* genCont,
			      Histos*    h,
			      Double_t   tailCut);
    /** 
     * Print information to standard output 
     * 
     * @param option Ignored 
     */
    void Print(Option_t* option="") const;
    /** 
     * Set debug flag 
     *
     * @param lvl Debug level 
     */
    virtual void SetDebug(UShort_t lvl);
  protected:
    Container* fSubs;
    Double_t   fLow;
    Double_t   fHigh;
    TH1*       fIPz;  //! 
    TH1*       fCent; //! 
    Histos*    fMeasured; 
    Histos*    fInjection;

    ClassDef(CentBin,1);
  };

protected:
  // -----------------------------------------------------------------
  /** 
   * @{ 
   * @name Initialization 
   */
  /** 
   * Initialize the task on worker 
   * 
   * @return true on success 
   */
  virtual Bool_t WorkerInit();
  /** 
   * Initialize our centrality bins.  
   * 
   * @param existing If null, also create new sub-sets.  If non-null,
   * read sub-sets from the container passed (using in Terminate).
   */
  Bool_t InitCentBins(Container* existing);
  /** 
   * Make a centrality bin 
   * 
   * @param c1 Low edge 
   * @param c2 High edge 
   * 
   * @return 
   */
  virtual CentBin* MakeCentBin(Float_t c1, Float_t c2)
  {
    return new CentBin(c1, c2);
  }
  /* @} */
  
  // -----------------------------------------------------------------
  /** 
   * @{ 
   * @name Event inspection 
   */
  virtual const char* GetBranchName() const { return "AliAODTracklets"; }
  /** 
   * Check event 
   * 
   * @param cent        On return, the event centrality 
   * @param ip          On return, the event interaction point 
   * @param tracklets   On return, the list of tracklets 
   * 
   * @return true if all selections pass x5
   */
  Bool_t CheckEvent(Double_t&          cent,
		    const AliVVertex*& ip,
		    TClonesArray*&     tracklets);
  /** 
   * Find the tracklet list in the event 
   * 
   * @param event Event 
   * 
   * @return List of tracklets or null
   */
  TClonesArray* FindTracklets(AliVEvent* event);
  /** 
   * See if we can find an IP position 
   * 
   * @param event         Event 
   * @param maxDispersion Max uncertainty 
   * @param maxZError     Max uncertainty in Z 
   * 
   * @return Pointer to IP position or null
   */
  const AliVVertex* FindSimpleIP(AliVEvent* event,
				 Double_t   maxDispersion=0.04,
				 Double_t   maxZError=0.25);
  /** 
   * See if we can find an IP position 
   * 
   * @param event         Event 
   * @param maxDispersion Max uncertainty 
   * @param maxZError     Max uncertainty in Z 
   * 
   * @return Pointer to IP position or null
   */
  const AliVVertex* FindRealIP(AliVEvent* event,
			       Double_t   maxDispersion=0.04,
			       Double_t   maxZError=0.25);
  /** 
   * Check if an IP position is OK
   * 
   * @param ip            IP position 
   * @param maxDispersion Max uncertainty 
   * @param maxZError     Max uncertainty in Z 
   * 
   * @return Pointer to IP position or null
   */
  const AliVVertex* CheckIP(const AliVVertex* ip,
			    Double_t   maxDispersion=0.04,
			    Double_t   maxZError=0.25);
				 
  /** 
   * Find the interaction point location 
   * 
   * @param event         Event 
   * @param maxDispersion Max uncertainty 
   * @param maxZError     Max uncertainty in Z 
   * 
   * @return Pointer to vertex, or null in case of problems 
   */
  const AliVVertex* FindIP(AliVEvent* event,
			   Double_t   maxDispersion=0.04,
			   Double_t   maxZError=0.25);
  /** 
   * Find the centrality of the event. This looks for the
   * AliMultSelection structure.
   * 
   * @param event Event 
   * @param nTracklets Number of tracklets for benchmarking centrality 
   * 
   * @return Centrality percentile, or negative number in case of
   * problems
   */
  Double_t FindMultCentrality(AliVEvent* event, Int_t& nTracklets);
  /** 
   * Find the centrality of the event. This looks for the
   * AliCentrality structure.
   * 
   * @param event Event 
   * 
   * @return Centrality percentile, or negative number in case of
   * problems
   */
  Double_t FindCompatCentrality(AliVEvent* event);  
  /** 
   * Find the centrality of the event.  If an AliMultSelection
   * structure isn't found, fall back to AliCentrality, and that isn't
   * there either, we fail.
   * 
   * @param event      Event 
   * @param nTracklets Number of tracklets for benchmarking centrality 
   * 
   * @return Centrality percentile, or negative number in case of
   * problems
   */
  Double_t FindCentrality(AliVEvent* event, Int_t& nTracklets);
  /** 
   * Check if we got a selected trigger 
   * 
   * @return true if the event was triggered 
   */  
  Bool_t FindTrigger();
  /* @} */

  // -----------------------------------------------------------------
  /** 
   * @{ 
   * @name Tracklet selection and inspection 
   */
  /** 
   * Find the tracklet weight 
   * 
   * @param tracklet Tracklet 
   * @param cent     Centrality 
   * 
   * @return The weight - in this class always 1
   */
  virtual Double_t LookupWeight(AliAODTracklet* tracklet, Double_t cent);
  /** 
   * Check if the tracklet passed is a signal tracklet 
   * 
   * @param tracklet Tracklet to inspect 
   * 
   * @return Bit mask of status 
   * - 0x1 if tracklet @f$\Delta<\Delta_{\mathrm{cut}}@f$ 
   * - 0x2 if tracklet @f$d\phi-\delta\phi <C@f$
   */
  virtual UShort_t CheckTracklet(AliAODTracklet* tracklet) const;
  /* @} */

  // -----------------------------------------------------------------
  /** 
   * @{ 
   * @name Event processing 
   */
  /** 
   * Process tracklets of an event 
   * 
   * @param cent        Event centrality 
   * @param ip          Event interaction point 
   * @param tracklets   List of tracklets 
   */
  void ProcessEvent(Double_t cent,const AliVVertex* ip,TClonesArray* tracklets);
  /* @} */

  // -----------------------------------------------------------------
  /** 
   * @{ 
   * @name Finalize the job 
   */
  /** 
   * Finalize the job on the master 
   * 
   * @param results The container to add the results to 
   *
   * @return true on success 
   */
  virtual Bool_t MasterFinalize(Container* results);
  /* @} */

  // -----------------------------------------------------------------
  /** Status of task */
  enum {
    kAll=1,          // Count all events
    kEvent,          // Have event
    kTracklets,      // Have tracklets
    kTrigger,        // Have trigger 
    kIP,             // Have IP
    kCentrality,     // Have centrality
    kCompleted       // Have completed
  };
  
  // -----------------------------------------------------------------
  /** Container */
  Container* fContainer;
  /** List of centrality bins */
  Container* fCentBins;
  /** Histogram of Interaction point z coordinate */
  TH1* fIPz;       //! 
  /** Histogram of centrality */
  TH1* fCent;      //! 
  /** Histogram of centrality */
  TH1* fStatus;    //! 
  /** Histogram of all eta phi */
  TH2* fEtaPhi;    //!
  /** Histogram of centrality, nTracklets correlation */
  TProfile* fCentTracklets; //! 
  /** Centrality method to use */
  TString    fCentMethod;
  /** Centrality axis */
  TAxis      fCentAxis;
  /** Interaction point Z axis */
  TAxis      fIPzAxis;
  /** Pseudorapidity axis */
  TAxis      fEtaAxis;
  /** Azimuthal angle axis */
  TAxis      fPhiAxis;
  /** Maximum @f$\Delta@f$ to consider */
  Double_t   fMaxDelta;
  /** Least value of @f$\Delta@f$ considered background tail */
  Double_t   fTailDelta;
  /** Shift @f$\delta_{\phi}@f$ of @f$\Delta\phi@f$ */
  Double_t   fDPhiShift;
  /** Signal cut on @f$\Delta\phi-\delta_{\phi}@f$ */
  Double_t   fShiftedDPhiCut;
  /** Signal cut on @f$\Delta@f$ */
  Double_t   fDeltaCut;
  /** The absolute minimum centrality to consider  - for MC with poor match*/
  Double_t   fAbsMinCent;
  
  ClassDef(AliTrackletAODdNdeta,1); 
};
//====================================================================
/**
 * Task to analyse AOD tracklets for dNch/deta 
 * 
 * @ingroup pwglf_forward_tracklets
 */
class AliTrackletAODMCdNdeta : public AliTrackletAODdNdeta
{
public:
  /** Type of containers */
  typedef TList Container;
  /** 
   * Default constructor - for ROOT I/O only 
   */
  AliTrackletAODMCdNdeta()
    : AliTrackletAODdNdeta()
  {}
  /** 
   * Named - user - constructor
   */
  AliTrackletAODMCdNdeta(const char* name)
    : AliTrackletAODdNdeta(name)
  {}
  /**
   * Copy constructor 
   *
   * @param o Object to copy from 
   */
  AliTrackletAODMCdNdeta(const AliTrackletAODdNdeta& o)
    : AliTrackletAODdNdeta(o)
  {}
  /**
   * Destructor 
   */
  virtual ~AliTrackletAODMCdNdeta() {}
  //__________________________________________________________________
  /**
   * A centrality bin.  Here, we need 
   *
   * - @f$\eta,\Delta@f$ distribution of measured @f$ f_M(\Delta)@f$ 
   * - @f$\eta,\Delta@f$ distribution of injected @f$ f_I(\Delta)@f$ 
   * - @f$\eta,\mathrm{IP}_z@f$ distribution of measured @f$ M_{\eta,z}@f$ 
   * - @f$\eta,\mathrm{IP}_z@f$ distribution of injected @f$ I_{\eta,z}@f$  
   *
   * For MC we also need 
   *
   * - @f$\eta,\mathrm{IP}_z@f$ distribution of combi @f$ C_{\eta,z}@f$ 
   * - @f$\eta,\mathrm{IP}_z@f$ distribution of gen @f$ P_{\eta,z}@f$ 
   * 
   * 
   */
  struct CentBin : public AliTrackletAODdNdeta::CentBin
  {
    /** 
     * Default constructor - for ROOT I/O only  
     */
    CentBin()
      : AliTrackletAODdNdeta::CentBin(),
	fCombinatorics(0),
	fPrimaries(0),
	fSecondaries(0),
	fGenerated(0)
    {
    } 
    /** 
     * User constructor 
     * 
     * @param c1 Lower bound on centrality 
     * @param c2 Upper bound on centrality
     */
    CentBin(Double_t c1, Double_t c2);
    /** 
     * Copy constructor 
     * 
     * @param o Object to copy from 
     */
    CentBin(const CentBin& o)
      : AliTrackletAODdNdeta::CentBin(o),
	fCombinatorics(0),
	fPrimaries(0),
	fSecondaries(0),
	fGenerated(0)
    {}
    /**
     * Destructor 
     */
    virtual ~CentBin() {}
    /** 
     * Assignment operator 
     *
     * @return reference to this 
     */
    CentBin& operator=(const CentBin&) { return *this; }
  protected:
    Histos*    fCombinatorics;
    Histos*    fPrimaries;
    Histos*    fSecondaries;
    Histos*    fGenerated;

    ClassDef(CentBin,1);
  };

protected:
  // -----------------------------------------------------------------
  /** 
   * @{ 
   * @name Initialization 
   */
  /** 
   * Make a centrality bin 
   * 
   * @param c1 Low edge 
   * @param c2 High edge 
   * 
   * @return 
   */
  virtual AliTrackletAODdNdeta::CentBin* MakeCentBin(Float_t c1, Float_t c2)
  {
    return new CentBin(c1, c2);
  }
  /* @} */
  
  // -----------------------------------------------------------------
  /** 
   * @{ 
   * @name Event inspection 
   */
  virtual const char* GetBranchName() const { return "AliAODMCTracklets"; }
  /* @} */

  ClassDef(AliTrackletAODMCdNdeta,1); 
};
//====================================================================
/**
 * Task to analyse AOD tracklets for dNch/deta 
 * 
 * @ingroup pwglf_forward_tracklets
 */
class AliTrackletAODWeightedMCdNdeta : public AliTrackletAODMCdNdeta
{
public:
  /** Type of containers */
  typedef TList Container;
  /** 
   * Default constructor - for ROOT I/O only 
   */
  AliTrackletAODWeightedMCdNdeta()
    : AliTrackletAODMCdNdeta(),
      fWeights(0),
      fEtaWeight(0)
  {}
  /** 
   * Named - user - constructor
   */
  AliTrackletAODWeightedMCdNdeta(const char* name)
    : AliTrackletAODMCdNdeta(name),
      fWeights(0),
      fEtaWeight(0)
  {}
  /**
   * Copy constructor 
   *
   * @param o Object to copy from 
   */
  AliTrackletAODWeightedMCdNdeta(const AliTrackletAODWeightedMCdNdeta& o)
    : AliTrackletAODMCdNdeta(o),
      fWeights(0),
      fEtaWeight(0)
  {}
  /**
   * Destructor 
   */
  virtual ~AliTrackletAODWeightedMCdNdeta() {}
  /** 
   * Assignment operator 
   * 
   * @param o Object to assign from 
   * 
   * @return Reference to this object 
   */
  AliTrackletAODWeightedMCdNdeta&
  operator=(const AliTrackletAODWeightedMCdNdeta& o);
  /** 
   * Print information to standard output 
   * 
   * @param option Ignored 
   */
  void Print(Option_t* option="") const;

  // -----------------------------------------------------------------
  /** 
   * @{ 
   * @name Set parameters on the task 
   */
  /** 
   * Set the weights to use 
   * 
   * @param w 
   */
  void SetWeights(AliTrackletWeights* w) { fWeights = w; }
  /* @} */
protected:
  // -----------------------------------------------------------------
  /** 
   * @{ 
   * @name Tracklet selection and inspection 
   */
  /** 
   * Initialize the task on worker 
   * 
   * @return true on success 
   */
  Bool_t WorkerInit();
  /** 
   * Find the tracklet weight 
   * 
   * @param tracklet Tracklet 
   * @param cent     Centrality 
   * 
   * @return The weight - in this class always 1
   */
  virtual Double_t LookupWeight(AliAODTracklet* tracklet, Double_t cent);
  /* @} */
  // -----------------------------------------------------------------
  /** 
   * @{ 
   * @name Finalize the job 
   */
  /** 
   * Finalize the job on the master 
   * 
   * @param results The container to add the results to 
   *
   * @return true on success 
   */
  Bool_t MasterFinalize(Container* results);
  /* @} */
  // -----------------------------------------------------------------
  AliTrackletWeights* fWeights;
  TProfile2D*         fEtaWeight;
  
  ClassDef(AliTrackletAODWeightedMCdNdeta,1); 
};

//____________________________________________________________________
AliTrackletAODdNdeta::AliTrackletAODdNdeta()
  : AliAnalysisTaskSE(),
    fContainer(0),
    fCentBins(0),
    fIPz(0),
    fCent(0),
    fStatus(0),
    fEtaPhi(0),
    fCentTracklets(0),
    fCentMethod(""),
    fCentAxis(1,0,0),
    fIPzAxis(1,0,0),
    fEtaAxis(1,0,0),
    fPhiAxis(1,0,0),
    fMaxDelta(0),
    fTailDelta(0),
    fDPhiShift(0),
    fShiftedDPhiCut(0),
    fDeltaCut(0),
    fAbsMinCent(-1)
{}
//____________________________________________________________________
AliTrackletAODdNdeta::AliTrackletAODdNdeta(const char* name)
  : AliAnalysisTaskSE(name),
    fContainer(0),
    fCentBins(0),
    fIPz(0),
    fCent(0),
    fStatus(0),
    fEtaPhi(0),
    fCentTracklets(0),
    fCentMethod("V0M"),
    fCentAxis(10,0,100),
    fIPzAxis(30,-15,+15),
    fEtaAxis(16,-2,+2),
    fPhiAxis(100,0,TMath::TwoPi()),
    fMaxDelta(25),
    fTailDelta(5),
    fDPhiShift(0.0045),
    fShiftedDPhiCut(-1),
    fDeltaCut(1.5),
    fAbsMinCent(-1)
{
  FixAxis(fCentAxis, "Centrality [%]");
  FixAxis(fIPzAxis,  "IP_{#it{z}} [cm]");
  FixAxis(fEtaAxis,  "#eta");
  FixAxis(fPhiAxis,  "#phi");

  DefineOutput(1, Container::Class());
  DefineOutput(2, Container::Class());
}
//____________________________________________________________________
AliTrackletAODdNdeta::AliTrackletAODdNdeta(const AliTrackletAODdNdeta& o)
  : AliAnalysisTaskSE(o),
    fContainer(0),
    fCentBins(0),
    fIPz(0),
    fCent(0),
    fStatus(0),
    fEtaPhi(0),
    fCentTracklets(0),
    fCentMethod(o.fCentMethod),
    fCentAxis(o.fCentAxis),
    fIPzAxis(o.fIPzAxis),
    fEtaAxis(o.fEtaAxis),
    fPhiAxis(o.fPhiAxis),
    fMaxDelta(o.fMaxDelta),
    fTailDelta(o.fTailDelta),
    fDPhiShift(o.fDPhiShift),
    fShiftedDPhiCut(o.fShiftedDPhiCut),
    fDeltaCut(o.fDeltaCut),
    fAbsMinCent(o.fAbsMinCent)
{}
//____________________________________________________________________
AliTrackletAODdNdeta&
AliTrackletAODdNdeta::operator=(const AliTrackletAODdNdeta& o)
{
  if (&o == this) return *this;
  if (fContainer) {
    delete fContainer;
    fContainer = 0;
  }
  if (fCentBins) {
    delete fCentBins;
    fCentBins = 0;
  }
  fIPz            = 0;
  fCent           = 0;
  fCentMethod     = o.fCentMethod;
  fCentAxis       = o.fCentAxis;
  fIPzAxis        = o.fIPzAxis;
  fEtaAxis        = o.fEtaAxis;
  fPhiAxis        = o.fPhiAxis;
  fMaxDelta       = o.fMaxDelta;
  fTailDelta      = o.fTailDelta;
  fDPhiShift      = o.fDPhiShift;
  fShiftedDPhiCut = o.fShiftedDPhiCut;
  fDeltaCut       = o.fDeltaCut;
  fAbsMinCent     = o.fAbsMinCent;
  return *this;
}
//____________________________________________________________________
AliTrackletAODWeightedMCdNdeta&
AliTrackletAODWeightedMCdNdeta::operator=(const AliTrackletAODWeightedMCdNdeta&
					  o)
{
  if (&o == this) return *this;
  AliTrackletAODdNdeta::operator=(o);
  return *this;
}
//____________________________________________________________________
Bool_t AliTrackletAODdNdeta::Connect(const char* sumFile,
				     const char* resFile)
{
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    AliError("No analysis manager to connect to.");
    return false;
  }   

  // Add to the manager 
  mgr->AddTask(this);
  
  // Create and connect output containers 
  TString sumOut;
  TString resOut;
  if      (sumFile && sumFile[0] != '\0') sumOut = sumFile;
  if      (resFile && resFile[0] != '\0') resOut = resFile;
  else if (sumFile && sumFile[0] != '\0') resOut = sumFile;

  // If the string is null or 'default' connect to standard output file 
  if (sumOut.IsNull() || sumOut.EqualTo("default", TString::kIgnoreCase)) 
    sumOut = AliAnalysisManager::GetCommonFileName();
  // If the string is null or 'default' connect to standard output file 
  if (resOut.IsNull() || resOut.EqualTo("default", TString::kIgnoreCase)) 
    resOut = AliAnalysisManager::GetCommonFileName();

  // Always connect input 
  mgr->ConnectInput(this, 0, mgr->GetCommonInputContainer());

  // Connect sum list unless the output 'none' is specified
  if (!sumOut.EqualTo("none", TString::kIgnoreCase)) {
    TString sumName(Form("%sSums", GetName()));      
    AliAnalysisDataContainer* sumCon = 
      mgr->CreateContainer(sumName, TList::Class(), 
			   AliAnalysisManager::kOutputContainer, sumOut);
    mgr->ConnectOutput(this, 1, sumCon);
  }

  // Connect the result list unless the output 'none' is specified
  if (!resOut.EqualTo("none", TString::kIgnoreCase)) {
    TString resName(Form("%sResults", GetName()));
    AliAnalysisDataContainer* resCon = 
      mgr->CreateContainer(resName, TList::Class(), 
			   AliAnalysisManager::kParamContainer, resOut);
    mgr->ConnectOutput(this, 2, resCon);
  }
  return true;
}
//____________________________________________________________________
void AliTrackletAODdNdeta::Print(Option_t* option) const
{
  Double_t shiftedDPhiCut = fShiftedDPhiCut;
  if (shiftedDPhiCut < 0) shiftedDPhiCut = TMath::Sqrt(fDeltaCut)*0.06;
  
  Printf("%s: %s", ClassName(), GetName());
  Printf(" %22s: 0x%08x", "Off-line trigger mask", fOfflineTriggerMask);
  Printf(" %22s: %f",   "Delta phi shift",	   fDPhiShift);
  Printf(" %22s: %f",   "Shifted Delta phi cut",   shiftedDPhiCut);
  Printf(" %22s: %f",   "Delta cut",	           fDeltaCut);
  Printf(" %22s: %f",   "max Delta",	           fMaxDelta);
  Printf(" %22s: %f",   "tail Delta",	           fTailDelta);
  Printf(" %22s: %f%%", "Absolute least c",     fAbsMinCent);
  PrintAxis(fEtaAxis);
  PrintAxis(fPhiAxis);
  PrintAxis(fIPzAxis,1,"IPz");
  PrintAxis(fCentAxis,0);

  if (!fCentBins) return;

  Printf("--- Centrality bins");
  TIter next(fCentBins);
  CentBin* bin = 0;
  while ((bin = static_cast<CentBin*>(next()))) {
    bin->Print(option);
  }
}
//____________________________________________________________________
void AliTrackletAODWeightedMCdNdeta::Print(Option_t* option) const
{
  AliTrackletAODMCdNdeta::Print(option);
  if (!fWeights) return;
  fWeights->Print(option);
}

//____________________________________________________________________
void AliTrackletAODdNdeta::CentBin::Print(Option_t* option) const
{
  Printf(" Centrality bin: %s", fName.Data());
  Printf("  Low cut:       %5.1f", fLow);
  Printf("  High cut:      %5.1f", fHigh);

  if (!fSubs) return;
  Printf("  --- Histogram sets");
  TIter next(fSubs);
  Histos* h = 0;
  while ((h = static_cast<Histos*>(next()))) {
    h->Print(option);
  }
}

//____________________________________________________________________
void AliTrackletAODdNdeta::CentBin::SetDebug(UShort_t lvl)
{
  Sub::SetDebug(lvl);
  TIter next(fSubs);
  Histos* h = 0;
  while ((h = static_cast<Histos*>(next()))) {
    h->SetDebug(lvl);
  }
}  
namespace {
  void Bits2String(UChar_t m, char out[7])
  {
    if (m & AliAODTracklet::kInjection)     out[0] = 'I'; else out[0] = '-';
    if (m & AliAODTracklet::kCombinatorics) out[1] = 'C'; else out[1] = '-';
    if (m & AliAODTracklet::kSecondary)     out[2] = 'S'; else out[2] = '-';
    if (m & AliAODTracklet::kDistinct)      out[3] = 'D'; else out[3] = '-';
    if (m & AliAODTracklet::kSimulated)     out[4] = 'X'; else out[4] = '-';
    if (m & AliAODTracklet::kGenerated)     out[5] = 'G'; else out[5] = '-';
    out[6] = '\0';    
  }
}

//____________________________________________________________________
void AliTrackletAODdNdeta::Histos::Print(Option_t*) const
{
  char cMask[7]; Bits2String(fMask, cMask);
  char cVeto[7]; Bits2String(fVeto, cVeto);  
  Printf("  Histograms: %s", fName.Data());
  Printf("   Mask:         0x%02x (%s)", fMask, cMask);
  Printf("   Veto:         0x%02x (%s)", fVeto, cVeto);
  Printf("   Delta:        %s", fEtaDeltaIPz ? "yes" : "no");
}

  
//____________________________________________________________________
void 
AliTrackletAODdNdeta::UserCreateOutputObjects()
{
  if (!WorkerInit()) {
    AliWarning("Failed to initialize on worker");
    return;
  }
  PostData(1,fContainer);
}

//____________________________________________________________________
Bool_t AliTrackletAODdNdeta::WorkerInit()
{
  if (DebugLevel() > 1) Printf("Initialising on worker");
  if (fShiftedDPhiCut < 0) fShiftedDPhiCut = TMath::Sqrt(fDeltaCut)*0.06;
  fContainer = new Container;
  fContainer->SetName(Form("%sSums", GetName()));
  fContainer->SetOwner();

  fIPz    = Make1D(fContainer, "ipz",  "", kMagenta+2, 20, fIPzAxis);
  fCent   = Make1D(fContainer, "cent", "", kMagenta+2, 20, fCentAxis);
  fEtaPhi = Make2D(fContainer, "etaPhi","",kMagenta+2, 20, fEtaAxis,fPhiAxis);
  fCentTracklets = new TProfile("centTracklets",
				"Mean number of tracklets per centrality",
				1050, 0, 105);
  fCentTracklets->SetXTitle("Centrality [%]");
  fCentTracklets->SetYTitle("#LTtracklets#GT");
  fCentTracklets->SetDirectory(0);
  fCentTracklets->SetMarkerStyle(20);
  fCentTracklets->SetMarkerColor(kMagenta+2);
  fCentTracklets->SetLineColor(kMagenta+2);
  fContainer->Add(fCentTracklets);

  fStatus = new TH1F("status", "Status of task",
		     kCompleted, .5, kCompleted+.5);
  fStatus->SetMarkerSize(2);
  fStatus->SetMarkerColor(kMagenta+2);
  fStatus->SetLineColor(kMagenta+2);
  fStatus->SetFillColor(kMagenta+2);
  fStatus->SetFillStyle(1001);
  fStatus->SetBarOffset(0.1);
  fStatus->SetBarWidth(0.4);
  fStatus->SetDirectory(0);
  fStatus->SetStats(0);
  fStatus->SetXTitle("Event have");
  fStatus->SetYTitle("# Events");
  fStatus->GetXaxis()->SetBinLabel(kAll,            "Been seen");
  fStatus->GetXaxis()->SetBinLabel(kEvent,          "Event data");
  fStatus->GetXaxis()->SetBinLabel(kTracklets,      "Tracklets");
  fStatus->GetXaxis()->SetBinLabel(kTrigger,        "Trigger");
  fStatus->GetXaxis()->SetBinLabel(kIP,             "IP");
  fStatus->GetXaxis()->SetBinLabel(kCentrality,     "Centrality");
  fStatus->GetXaxis()->SetBinLabel(kCompleted,      "Completed");
  fContainer->Add(fStatus);

  typedef TParameter<double> DP;
  typedef TParameter<bool>   BP;
  typedef TParameter<int>    IP;
  Container* params = new Container;
  params->SetName("parameters");
  params->SetOwner();
  fContainer->Add(params);
  params->Add(new DP("DPhiShift",      fDPhiShift,      'f'));
  params->Add(new DP("ShiftedDPhiCut", fShiftedDPhiCut, 'f'));
  params->Add(new DP("DeltaCut",       fDeltaCut,       'f'));
  params->Add(new DP("MaxDelta",       fMaxDelta,       'f'));
  params->Add(new DP("TailDelta",      fTailDelta,      'f'));
  params->Add(new DP("AbsMinCent",     fAbsMinCent,     'f'));
  // Create our centrality bins 
  if (!InitCentBins(0)) {
    AliWarning("Failed to initialize centrality bins");
    return false;
  }

  // Print information to log
  Print();
  return true;
}
//____________________________________________________________________
Bool_t AliTrackletAODWeightedMCdNdeta::WorkerInit()
{
  if (!AliTrackletAODMCdNdeta::WorkerInit()) return false;
  fEtaWeight = 0;
  if (!fWeights) {
    AliFatal("No weights set!");
    return false;
  }
  
  fEtaWeight = Make2P(fContainer, "etaWeight", "#LTw#GT", kYellow+2, 24,
		      fEtaAxis, fCentAxis);
  fWeights->Store(fContainer);
  return true;
}
//____________________________________________________________________
Bool_t AliTrackletAODdNdeta::InitCentBins(Container* existing)
{
  if (DebugLevel() > 1)
    Printf("Initialising on centrality bins on %s",
	   existing ? "master" : "worker");
  if (fCentBins) return true;
  fCentBins = new Container;
  fCentBins->SetName("centralityBins");
  fCentBins->SetOwner();

  TAxis    deltaAxis (Int_t(5*fMaxDelta),0,         fMaxDelta);
  FixAxis(deltaAxis,
	  "#Delta=[(#Delta#phi-#delta#phi)/#sigma_{#phi}]^{2}+"
	  "[#Delta#thetasin^{-2}(#theta)/#sigma_{#theta}]^{2}");

  // Add min-bias bin
  Bool_t   ret  = true;
  CentBin* bin  = MakeCentBin(0, 100);
  if (!existing) 
    ret = bin->WorkerInit(fContainer,fEtaAxis,fIPzAxis,deltaAxis);
  else
    ret = bin->FinalizeInit(existing);
  if (!ret) {
    AliWarningF("Failed to initialize bin %s", bin->GetName());
    return false;
  }
  bin->SetDebug(DebugLevel());
  fCentBins->AddAt(bin, 0);

  // Add other bins
  Int_t nCentBins = fCentAxis.GetNbins();
  for (Int_t i = 1; i <= nCentBins; i++) {
    Float_t  c1 = fCentAxis.GetBinLowEdge(i);
    Float_t  c2 = fCentAxis.GetBinUpEdge(i);
    bin         = MakeCentBin(c1, c2);
    if (!existing) 
      ret = bin->WorkerInit(fContainer,fEtaAxis,fIPzAxis,deltaAxis);
    else
      ret = bin->FinalizeInit(existing);
    if (!ret) {
      AliWarningF("Failed to initialize %s", bin->GetName());
      return false;
    }
    bin->SetDebug(DebugLevel());
    fCentBins->AddAt(bin, i);
  }
  return true;
}

//____________________________________________________________________
AliTrackletAODdNdeta::CentBin::CentBin(Double_t c1, Double_t c2)
  : Sub(""),
    fSubs(0),
    fLow(c1),
    fHigh(c2),
    fIPz(0),
    fCent(0),
    fMeasured(0),
    fInjection(0)
{
  fName.Form("cent%03dd%02d_%03dd%02d",
	     Int_t(fLow), Int_t(fLow*100)%100,
	     Int_t(fHigh), Int_t(fHigh*100)%100);
  fMeasured  = new Histos("measured", kRed+2, 20,
			  0x00, // No requirements, just veto 
			  AliAODTracklet::kInjection|
			  AliAODTracklet::kGenerated);
  fInjection = new Histos("injected", kOrange+2, 21,
			  AliAODTracklet::kInjection,
			  AliAODTracklet::kGenerated);
  fSubs = new Container;
  fSubs->SetOwner(true);
  fSubs->Add(fMeasured);
  fSubs->Add(fInjection);
}

//____________________________________________________________________
AliTrackletAODMCdNdeta::CentBin::CentBin(Double_t c1, Double_t c2)
  : AliTrackletAODdNdeta::CentBin(c1, c2),
    fCombinatorics(0),
    fPrimaries(0),
    fSecondaries(0),
    fGenerated(0)
{
  // Combinatorics is everything that has the combinatorics bit
  // Primaries all with simulated bit, but not secondary or
  // combinatorics bit.
  // Secondaries are all those with the secondary bit set 
  fCombinatorics  = new Histos("combinatorics", kMagenta+2, 30, 
			       AliAODTracklet::kCombinatorics, 0x00);
  fPrimaries  = new Histos("primaries", kGreen+2, 26, 
			   0,
			   AliAODTracklet::kInjection|
			   AliAODTracklet::kCombinatorics|
			   AliAODTracklet::kSecondary|
			   AliAODTracklet::kGenerated);
  fSecondaries = new Histos("secondaries", kBlue+2, 32, 
			    AliAODTracklet::kSecondary, 0x00);
  fGenerated = new Histos("generated", kGray+1, 28, 
			  AliAODTracklet::kGenerated, 0x00);
  fMeasured->fStyle = 24;
  fInjection->fStyle = 25;
  fSubs->Add(fCombinatorics);
  fSubs->Add(fPrimaries);
  fSubs->Add(fSecondaries);
  fSubs->AddAfter(fMeasured, fGenerated);
}

//____________________________________________________________________
Bool_t AliTrackletAODdNdeta::CentBin::WorkerInit(Container* parent,
						 const TAxis& etaAxis,
						 const TAxis& ipzAxis,
						 const TAxis& deltaAxis)
{
  if (fDebug > 0)
    Printf("Initializing centrality bin %s", fName.Data());
  if (!Sub::WorkerInit(parent, etaAxis, ipzAxis, deltaAxis)) return false;
  
  TAxis centAxis(20, fLow, fHigh);
  FixAxis(centAxis, "Centrality [%]");
  fCent      = Make1D(fContainer,"cent","Centrality [%]",
		      kMagenta+2,20,centAxis);
  fIPz       = Make1D(fContainer,"ipz","IP_{#it{z}} [cm]",kRed+2,20,ipzAxis);

  TIter   next(fSubs);
  Histos* h = 0;
  while ((h = static_cast<Histos*>(next()))) {
    if (!h ->WorkerInit(fContainer, etaAxis,ipzAxis,deltaAxis))
      return false;
  }
      
  return true;
}
//____________________________________________________________________
Bool_t AliTrackletAODdNdeta::Histos::WorkerInit(Container* parent,
						const TAxis& etaAxis,
						const TAxis& ipzAxis,
						const TAxis& deltaAxis)
{
  if (!Sub::WorkerInit(parent, etaAxis, ipzAxis, deltaAxis)) return false;
  TString shrt(Form("%c", GetName()[0]));
  shrt.ToUpper();
  if (GetC(parent, "generated", false) != 0) shrt.Append("'");
  
  // Do not make eta vs IPz for secondaries and primaries 
  if (fMask != AliAODTracklet::kSecondary &&
      (fVeto != (AliAODTracklet::kSecondary     |
		 AliAODTracklet::kInjection     |
		 AliAODTracklet::kCombinatorics |
		 AliAODTracklet::kGenerated)))
    fEtaIPz   = Make2D(fContainer, "etaIPz", shrt.Data(),
		       kRed+2, 20, etaAxis, ipzAxis);
  // Always make eta vs Delta distribution
  if (fMask != AliAODTracklet::kGenerated) 
    fEtaDeltaIPz = Make3D(fContainer, "etaDeltaIPz",
			  Form("#Delta_{%s}",shrt.Data()), 
			  kBlue+2, 21, etaAxis, deltaAxis, ipzAxis);
  SetAttr(fColor, fStyle);
  return true;
}
//____________________________________________________________________
void AliTrackletAODdNdeta::Histos::SetAttr(Color_t c, Style_t m)
{
  if (fEtaIPz) {
    fEtaIPz->SetMarkerStyle(m);
    fEtaIPz->SetMarkerColor(c);
    fEtaIPz->SetLineColor(c);
    fEtaIPz->SetFillColor(c);
  }
  if (fEtaDeltaIPz) {
    fEtaDeltaIPz->SetMarkerStyle(m);
    fEtaDeltaIPz->SetMarkerColor(c);
    fEtaDeltaIPz->SetLineColor(c);
    fEtaDeltaIPz->SetFillColor(c);
  }
}

//____________________________________________________________________
void 
AliTrackletAODdNdeta::UserExec(Option_t*)
{
  if (DebugLevel() > 0) Printf("In user exec");
  Double_t          cent      = -1;
  const AliVVertex* ip        = 0;
  TClonesArray*     tracklets = 0;
  if (!CheckEvent(cent, ip, tracklets)) {
    AliWarningF("Event didn't pass %f, %p, %p", cent, ip, tracklets);
    Printf("Argh, check data failed %f, %p, %p", cent, ip, tracklets);
    return;
  }
  if (DebugLevel() > 0) Printf("Got centrality=%f ipZ=%f %d tracklets",
			       cent, ip->GetZ(), tracklets->GetEntriesFast());
  ProcessEvent(cent, ip, tracklets);

  PostData(1,fContainer);
  fStatus->Fill(kCompleted);
}

//____________________________________________________________________
Bool_t AliTrackletAODdNdeta::CheckEvent(Double_t&          cent,
					const AliVVertex*& ip,
					TClonesArray*&     tracklets)
{
  // Count all events 
  fStatus->Fill(kAll);

  // Check for event 
  AliVEvent* event = InputEvent();
  if (!event) {
    AliWarning("No event");
    return false;
  }
  fStatus->Fill(kEvent);

  // Check if we have the tracklets 
  tracklets = FindTracklets(event);
  if (!tracklets) return false;
  fStatus->Fill(kTracklets);
    
  // Check if event was triggered 
  Bool_t trg = FindTrigger();
  if (!trg) return false;
  fStatus->Fill(kTrigger);
    
  // Check the interaction point 
  ip = FindIP(event);
  if (!ip) return false;
  fStatus->Fill(kIP);

  // Check the centrality
  Int_t nTracklets = 0;
  cent = FindCentrality(event, nTracklets);
  if (cent < 0) return false;
  fStatus->Fill(kCentrality);

  fIPz->Fill(ip->GetZ());
  fCent->Fill(cent);
  fCentTracklets->Fill(cent, nTracklets); 
  
  return true;
}

//____________________________________________________________________
TClonesArray* AliTrackletAODdNdeta::FindTracklets(AliVEvent* event)
{
  // Check the multiplicity
  TObject* obj = event->FindListObject(GetBranchName());
  if (!obj) {
    AliWarningF("Couldn't get object %s", GetBranchName());
    // event->GetList()->Print();
    return 0;
  }
  if (!obj->IsA()->InheritsFrom(TClonesArray::Class())) {
    AliWarningF("Object %s is not a TClonesArray but a %s",
		obj->GetName(), obj->ClassName());
    return 0;
  }
  return static_cast<TClonesArray*>(obj);
}  

//____________________________________________________________________
const AliVVertex* AliTrackletAODdNdeta::CheckIP(const AliVVertex* ip,
						Double_t maxDispersion,
						Double_t maxZError)
{
  if (!ip) return 0;
  if (ip->GetNContributors() <= 0) {
    AliWarning("Not enough contributors for IP");
    return 0;
  }   
  // If this is from the Z vertexer, do some checks 
  if (ip->IsFromVertexerZ()) {
    // Get covariance matrix
    Double_t covar[6];
    ip->GetCovarianceMatrix(covar);
    Double_t sigmaZ = TMath::Sqrt(covar[5]);
    if (sigmaZ >= maxZError) {
      AliWarningF("IPz resolution = %f >= %f", sigmaZ, maxZError);
      return 0;
    }
      
    // If this IP doesn not derive from AliVertex, don't check dispersion. 
    if (ip->IsA()->InheritsFrom(AliVertex::Class())) {
      const AliVertex* ipv = static_cast<const AliVertex*>(ip);
      // Dispersion is the parameter used by the vertexer for finding the IP. 
      if (ipv->GetDispersion() >= maxDispersion) {
	AliWarningF("IP dispersion = %f >= %f",
		    ipv->GetDispersion(), maxDispersion);
	return 0;
      }
    }
  }
    
  // If we get here, we either have a full 3D vertex or track
  // vertex, and we should check if it is in range
  if (ip->GetZ() < fIPzAxis.GetXmin() || ip->GetZ() > fIPzAxis.GetXmax()) {
    AliWarningF("IPz = %fcm out of range [%f,%f]cm",
		ip->GetZ(), fIPzAxis.GetXmin(), fIPzAxis.GetXmax());
    return 0;
  }
  return ip;
}

//____________________________________________________________________
const AliVVertex* AliTrackletAODdNdeta::FindRealIP(AliVEvent* event,
						   Double_t maxDispersion,
						   Double_t maxZError)
{
  const AliVVertex* ip   = event->GetPrimaryVertex();
  if (!ip) {
    if (DebugLevel() > 1) AliWarning("No real IP for this event found!");
    return 0;
  }
  return CheckIP(ip, maxDispersion, maxZError);
}
						   
//____________________________________________________________________
const AliVVertex* AliTrackletAODdNdeta::FindSimpleIP(AliVEvent* event,
						     Double_t maxDispersion,
						     Double_t maxZError)   
{
  static AliVertex* ip;
  AliAODSimpleHeader* head =
    static_cast<AliAODSimpleHeader*>(event
				     ->FindListObject("AliAODSimpleHeader"));
  if (!head) {
    if (DebugLevel() > 1) AliWarning("No simple header");
    return 0;
  }
  if (!ip)
    ip = new AliVertex;
  ip->SetXv(head->fRecIP.X());
  ip->SetYv(head->fRecIP.Y());
  ip->SetZv(head->fRecIP.Z());
  ip->SetNContributors(10);
  ip->SetTitle("");
  return CheckIP(ip, maxDispersion, maxZError);
}
  
//____________________________________________________________________
const AliVVertex* AliTrackletAODdNdeta::FindIP(AliVEvent* event,
					       Double_t maxDispersion,
					       Double_t maxZError)
{
  const AliVVertex* ip   = FindRealIP(event,maxDispersion,maxZError);  
  if (!ip) {
    ip = FindSimpleIP(event,maxDispersion,maxZError);
    if (!ip) {
      if (DebugLevel() > 1) AliWarning("No IP for this event found!");
      return 0;
    }
  }
  // Good vertex, return it
  return ip;
}
//____________________________________________________________________
Bool_t AliTrackletAODdNdeta::FindTrigger()
{
  UInt_t evBits = fInputHandler->IsEventSelected();
  Bool_t trgOK  = (evBits & fOfflineTriggerMask || fOfflineTriggerMask == 0);
  if (DebugLevel())
    Printf("Trigger bits=0x%08x  mask=0x%08x  masked=0x%08x -> %s",
	   evBits, fOfflineTriggerMask, evBits & fOfflineTriggerMask,
	   trgOK ? "selected" : "rejected");
  return trgOK;
}
//____________________________________________________________________
Double_t AliTrackletAODdNdeta::FindCompatCentrality(AliVEvent* event)
{
  static TString centMeth;
  if (centMeth.IsNull()) {
    centMeth = fCentMethod(3,fCentMethod.Length()-3);
  }
  AliCentrality* cent   = event->GetCentrality();
  if (!cent) return -1;

  Double_t centPer = cent->GetCentralityPercentileUnchecked(centMeth);
  return centPer;
}
//____________________________________________________________________
Double_t AliTrackletAODdNdeta::FindMultCentrality(AliVEvent* event,
						  Int_t& nTracklets)
{
  AliMultSelection* cent =
    static_cast<AliMultSelection*>(event->FindListObject("MultSelection"));
  if (!cent) {
    AliWarning("No centrality in event");
    event->GetList()->Print();
    return -1;
  }
  AliMultEstimator* estTracklets = cent->GetEstimator("SPDTracklets");
  if (estTracklets)    
    nTracklets = estTracklets->GetValue();

  return cent->GetMultiplicityPercentile(fCentMethod);
}
  
//____________________________________________________________________
Double_t AliTrackletAODdNdeta::FindCentrality(AliVEvent* event,
					      Int_t& nTracklets)
{
  if (fCentMethod.EqualTo("MB", TString::kIgnoreCase)) {
    Printf("MB centrality - not checked");
    return 0;
  }
  Double_t centPer = -1;
  if (fCentMethod.BeginsWith("OLD"))
    centPer = FindCompatCentrality(event);
  else
    centPer = FindMultCentrality(event, nTracklets);
  const Double_t safety = 1e-3;
  if (DebugLevel() > 1) Printf("Read centrality: %f%%", centPer);
  if      (centPer < -safety)    return -2;
  if      (centPer < +safety)    centPer = safety;
  else if (centPer > 100-safety) centPer = 100-safety;

  if (fAbsMinCent >= 0 && centPer < fAbsMinCent) {
    AliWarningF("Centrality = %f lower than absolute minimum [%f]",
		centPer, fAbsMinCent);
    return -3;
  }
  if (centPer < fCentAxis.GetXmin() || centPer > fCentAxis.GetXmax()) {
    AliWarningF("Centrality = %f out of range [%f,%f]",
		centPer, fCentAxis.GetXmin(), fCentAxis.GetXmax());
    return -3;
  }
  return centPer;    
}

//____________________________________________________________________
Double_t AliTrackletAODdNdeta::LookupWeight(AliAODTracklet* tracklet,
					    Double_t        cent)
{
  return 1;
}
//____________________________________________________________________
Double_t AliTrackletAODWeightedMCdNdeta::LookupWeight(AliAODTracklet* tracklet,
					      Double_t        cent)
{
  // We don't check for weights, as we must have them to come this far 
  // if (!fWeights) {
  // AliWarning("No weights defined");
  // return 1;
  // }
  Double_t w = fWeights->LookupWeight(tracklet, cent);
  fEtaWeight->Fill(tracklet->GetEta(), cent, w);
  // printf("Looking up weight of tracklet -> %f ", w);
  // tracklet->Print();
  return w;
}
//____________________________________________________________________
UShort_t AliTrackletAODdNdeta::CheckTracklet(AliAODTracklet* tracklet) const
{
  Double_t dPhiS   = (tracklet->GetDPhi() -
		      TMath::Sign(fDPhiShift,Double_t(tracklet->GetDPhi())));
  UShort_t ret = 0;
  ret |= ((tracklet->GetDelta() <= fDeltaCut) ? 0x1 : 0);
  ret |= ((dPhiS < fShiftedDPhiCut) ? 0x2 : 0);  
  // Printf("dPhiS=%f (%f) Delta=%f (%f) -> 0x%x",
  // 	 dPhiS, fShiftedDPhiCut,
  // 	 tracklet->GetDelta(), fDeltaCut,
  // 	 ret);
  return ret;  
}

    
//____________________________________________________________________
void AliTrackletAODdNdeta::ProcessEvent(Double_t          cent,
					const AliVVertex* ip,
					TClonesArray*     tracklets)
{
  // Figure out which centrality bins to fill 
  Int_t    nAcc = 0;
  TIter    nextAcc(fCentBins);
  CentBin* bin = 0;
  TList    toRun;
  while ((bin = static_cast<CentBin*>(nextAcc()))) {
    if (!bin->Accept(cent, ip->GetZ())) continue; // Not in range for this bin
    toRun.Add(bin);
    nAcc++;
  }
  // If we have no centrality bins  to fill, we return immediately 
  if (nAcc <= 0) return;

  AliAODTracklet* tracklet = 0;
  TIter           nextTracklet(tracklets);
  while ((tracklet = static_cast<AliAODTracklet*>(nextTracklet()))) {
    Double_t weight = LookupWeight(tracklet, cent);
    UShort_t signal = CheckTracklet(tracklet);
    if (signal) fEtaPhi->Fill(tracklet->GetEta(), tracklet->GetPhi());
    TIter nextBin(&toRun);
    while ((bin = static_cast<CentBin*>(nextBin()))) {
      bin->ProcessTracklet(tracklet, ip->GetZ(), signal, weight);
    }    
  }
}    

//____________________________________________________________________
Bool_t AliTrackletAODdNdeta::CentBin::Accept(Double_t cent, Double_t ipz)
{
  if (cent < fLow || cent >= fHigh) return false;
  fCent->Fill(cent);
  fIPz ->Fill(ipz);
  return true;
}

//____________________________________________________________________
Bool_t AliTrackletAODdNdeta::CentBin::ProcessTracklet(AliAODTracklet* tracklet,
						      Double_t        ipZ,
						      UShort_t        signal,
						      Double_t        weight)
{
  TIter   next(fSubs);
  Histos* h = 0;
  if (fDebug > 3) tracklet->Print();
  while ((h = static_cast<Histos*>(next()))) 
    h->ProcessTracklet(tracklet, ipZ, signal, weight);
  
  return true;
}
//____________________________________________________________________
Bool_t AliTrackletAODdNdeta::Histos::ProcessTracklet(AliAODTracklet* tracklet,
						     Double_t        ipZ,
						     UShort_t        signal,
						     Double_t        weight)
{
  if (!fEtaIPz && !fEtaDeltaIPz) return true;
  char m[7];
  if (fDebug > 3) Bits2String(tracklet->GetFlags(), m);
  if (fMask != 0 && (tracklet->GetFlags() & fMask) == 0) {
    if (fDebug > 3)
      Printf("%14s (0x%02x,----) %6s %7s (0x%02x) ",
	     GetName(), fMask, "reject", m, tracklet->GetFlags());
    return false;
  }
  if (fVeto != 0 && (tracklet->GetFlags() & fVeto) != 0) {
    if (fDebug > 3) 
      Printf("%14s (----,0x%02x) %6s %7s (0x%02x) ",
	     GetName(), fVeto, "veto", m, tracklet->GetFlags());
    return false;
  }
  if (fDebug > 3)
    Printf("%14s (0x%02x,0x%02x) %6s %7s (0x%02x) ",
	   GetName(), fMask, fVeto, "accept", m, tracklet->GetFlags());
  if (fEtaIPz && (signal == 0x3)) // both reguirements 
    fEtaIPz->Fill(tracklet->GetEta(), ipZ, weight);
  if (fEtaDeltaIPz && (signal & 0x2)) // just check dPhi
    fEtaDeltaIPz->Fill(tracklet->GetEta(), tracklet->GetDelta(), ipZ, weight);
  return true;
}


//____________________________________________________________________
void 
AliTrackletAODdNdeta::Terminate(Option_t*)
{
  Container* results = new Container;
  results->SetName(Form("%sResults",GetName()));
  results->SetOwner();

  Print("");
  fContainer = static_cast<Container*>(GetOutputData(1));
  if (!fContainer) {
    AliWarning("No sum container found!");
    return;
  }
    
  if (!InitCentBins(fContainer)) {
    AliWarningF("Failed to initialize centrality bins from %s",
		fContainer->GetName());
    return;
  }
    
  if (!MasterFinalize(results)) {
    AliWarning("Failed to finalize results");
    return;
  }
    
  PostData(2, results);
}


//____________________________________________________________________
Bool_t AliTrackletAODdNdeta::CentBin::FinalizeInit(Container* parent)
{
  fContainer     = GetC(parent, fName);
  fCent          = GetH1(fContainer, "cent");
  fIPz           = GetH1(fContainer, "ipz");
  if (!fContainer || !fCent || !fIPz) return false;
  TIter next(fSubs);
  Histos* h = 0;
  while ((h = static_cast<Histos*>(next()))) 
    if (!h->FinalizeInit(fContainer)) return false;
  return true;
}
//____________________________________________________________________
Bool_t AliTrackletAODdNdeta::Histos::FinalizeInit(Container* parent)
{
  fContainer   = GetC(parent, fName);
  fEtaIPz      = GetH2(fContainer, "etaIPz", false); // No complaints
  fEtaDeltaIPz = GetH3(fContainer, "etaDeltaIPz", false);
  if (GetName()[0] == 'm' && GetC(parent,"generated")) {
    // Fix up titles 
    if (fEtaIPz)
      fEtaIPz->SetTitle(Form("%s'", fEtaIPz->GetTitle()));
    if (fEtaDeltaIPz)
      fEtaDeltaIPz->SetTitle(Form("#Delta_{%s}", fEtaIPz->GetTitle()));
  }
  return (fContainer != 0); //  && fEtaIPz != 0); 
}

//____________________________________________________________________
Bool_t AliTrackletAODdNdeta::MasterFinalize(Container* results)
{
  // Make copies of histograms and store
  fIPz    = static_cast<TH1*>(CloneAndAdd(results,GetH1(fContainer,"ipz")));
  fCent   = static_cast<TH1*>(CloneAndAdd(results,GetH1(fContainer,"cent")));  
  fStatus = static_cast<TH1*>(CloneAndAdd(results,GetH1(fContainer,"status")));
  fEtaPhi = static_cast<TH2*>(CloneAndAdd(results,GetH2(fContainer,"etaPhi")));
  fCentTracklets =
    static_cast<TProfile*>(CloneAndAdd(results,GetP(fContainer,
						    "centTracklets")));
  Double_t nEvents = fIPz->GetEntries();
  Printf("Event summary:");
  for (Int_t i = 1; i <= fStatus->GetNbinsX(); i++) 
    Printf("  %10d %s",
	   Int_t(fStatus->GetBinContent(i)),
	   fStatus->GetXaxis()->GetBinLabel(i));
  for (Int_t i = 1; i <= fCent->GetNbinsX(); i++) 
    Printf("  %6.2f-%6.2f%%: %d",
	   fCent->GetXaxis()->GetBinLowEdge(i),
	   fCent->GetXaxis()->GetBinUpEdge(i),
	   Int_t(fCent->GetBinContent(i)));

  
  fIPz   ->Scale(1./nEvents);
  fCent  ->Scale(1./fCent->GetEntries());
  fStatus->Scale(1./fStatus->GetBinContent(1));
  fEtaPhi->Scale(1./nEvents);

  TIter    next(fCentBins);
  CentBin* bin = 0;
  while ((bin = static_cast<CentBin*>(next()))) {
    if (!bin->MasterFinalize(results, 0, fTailDelta)) {
      AliWarningF("Failed to finalize %s", bin->GetName());
      return false;
    }
  }
  return true;
}
//____________________________________________________________________
Bool_t AliTrackletAODWeightedMCdNdeta::MasterFinalize(Container* results)
{
  if (!AliTrackletAODMCdNdeta::MasterFinalize(results)) return false;

  TObject* o = fContainer->FindObject("etaWeight");
  if (o && o->IsA()->InheritsFrom(TH1::Class())) {
    TH1* etaWeight = static_cast<TH1*>(o->Clone());
    etaWeight->SetDirectory(0);
    results->Add(etaWeight);
  }
  else {
    AliWarningF("Object %p (etaWeight) is not a TH1 or not found",o);
  }
  if (!fWeights) return true;

  if (!fWeights->Retrieve(fContainer)) return false;
  fWeights->Store(results);

  return true;
}

  
//____________________________________________________________________
Bool_t AliTrackletAODdNdeta::CentBin::MasterFinalize(Container* parent,
						     TH1*       ,
						     Double_t   tailDelta)
{
  Container* result = new Container;
  result->SetName(fName);
  result->SetOwner(true);
  parent->Add(result);

  Double_t   nEvents = fIPz->GetEntries();
  // Copy ipZ histogram and scale by number of events 
  TH1* ipZ = static_cast<TH1*>(CloneAndAdd(result, fIPz));
  ipZ->Scale(1./nEvents);

  TH1* cent = static_cast<TH1*>(CloneAndAdd(result, fCent));
  cent->Scale(1./nEvents);

  Container* measCont = 0;
  Container* genCont  = 0;
  TIter      next(fSubs);
  Histos*    h = 0;
  while ((h = static_cast<Histos*>(next()))) {
    if (h->GetMask() != AliAODTracklet::kInjection &&
	h->GetMask() != AliAODTracklet::kCombinatorics) {
      // For the anything but injection or MC-labels, we just finalize
      if (!h->MasterFinalize(result, fIPz, tailDelta)) {
	AliWarningF("Failed to finalize %s/%s", GetName(), h->GetName());
	return false;
      }
      if (h->GetMask() == AliAODTracklet::kGenerated)
	// Set MC truth container 
	genCont = GetC(result, h->GetName());
      if (h == fMeasured)
	measCont = GetC(result, h->GetName());	
      continue;
    }
    if (!EstimateBackground(result, measCont, genCont, h, tailDelta)) {
      AliWarningF("Failed to estimate Bg in %s/%s", GetName(), h->GetName());
      return false;
    }      
  }
  return true;
}

//____________________________________________________________________
Bool_t 
AliTrackletAODdNdeta::Histos::MasterFinalize(Container* parent,
					     TH1*       ipz,
					     Double_t   tailDelta)
{
  Container* result = new Container;
  result->SetName(fName);
  result->SetOwner(true);
  parent->Add(result);

  // Get the number of events
  Double_t nEvents = ipz->GetEntries();
  
  // Scale each vertex range by number of events in that range
  TH2* etaIPz = 0;
  if (fEtaIPz) {
    etaIPz = ScaleToIPz(fEtaIPz, ipz);
    result->Add(etaIPz);
  }

  // If we do not have eta vs Delta, just return 
  if (!fEtaDeltaIPz) return true;

  // Normalize delta distribution to integral number of events
  // static_cast<TH3*>(CloneAndAdd(result, fEtaDeltaIPz));
  TH3* etaDeltaIPz = ScaleToIPz(fEtaDeltaIPz, ipz,
				ipz->GetEntries()>1000);
  result->Add(etaDeltaIPz);
  
  // Make 2D projection to eta,Delta
  TH2* etaDelta = ProjectEtaDelta(etaDeltaIPz);
  result->Add(etaDelta);
  
  // Make projection of delta 
  TH1* delta = ProjectDelta(etaDelta);
  result->Add(delta);

  // PArameters of integrals
  Double_t maxDelta = etaDeltaIPz->GetYaxis()->GetXmax();
  Int_t    lowBin   = etaDeltaIPz->GetYaxis()->FindBin(tailDelta);
  Int_t    highBin  = etaDeltaIPz->GetYaxis()->GetNbins();  

  
  TH1* etaDeltaTail    = etaDelta->ProjectionX("etaDeltaTail");
  etaDeltaTail   ->SetDirectory(0);
  etaDeltaTail   ->Reset();
  etaDeltaTail   ->SetTitle(Form("#scale[.7]{#int}_{%3.1f}^{%4.1f}"
				 "d%s d#it{N}/d%s",
				 tailDelta,maxDelta,
				 etaDeltaIPz->GetTitle(), 
				 etaDeltaIPz->GetTitle()));
  etaDeltaTail   ->SetYTitle("#scale[.7]{#int}_{tail}d#Delta d#it{N}/d#Delta");

  TH2* etaIPzDeltaTail = static_cast<TH2*>(etaDeltaIPz->Project3D("zx e"));
  etaIPzDeltaTail->SetName("etaIPzDeltaTail");
  etaIPzDeltaTail->SetDirectory(0);
  etaIPzDeltaTail->Reset();
  etaIPzDeltaTail->SetTitle(etaDeltaTail->GetTitle());
  etaIPzDeltaTail->SetZTitle(etaDelta->GetYaxis()->GetTitle());
  // Loop over eta
  Double_t intg = 0, eintg = 0;
  for (Int_t i = 1; i <= etaDeltaTail->GetNbinsX(); i++) {
    // Integrate over Delta 
    intg        = etaDelta->IntegralAndError(i, i, lowBin, highBin, eintg);    
    etaDeltaTail->SetBinContent(i, intg);
    etaDeltaTail->SetBinError  (i, eintg);
    // Loop over IPz
    for (Int_t j = 1; j <= etaIPzDeltaTail->GetNbinsY(); j++) {
      // Integrate over Delta 
      intg = etaDeltaIPz->IntegralAndError(i,i,lowBin,highBin,j,j,eintg);
      etaIPzDeltaTail->SetBinContent(i, j, intg);
      etaIPzDeltaTail->SetBinError  (i, j, eintg);
    }
  }
  result->Add(etaIPzDeltaTail);
  result->Add(etaDeltaTail);

  // Integrate full tail
  intg = etaDeltaIPz->IntegralAndError(1,etaDeltaIPz->GetNbinsX(),
				       lowBin, highBin,
				       1,etaDeltaIPz->GetNbinsZ(),
				       eintg);
  result->Add(new TParameter<double>("deltaTail",      intg));
  result->Add(new TParameter<double>("deltaTailError", eintg));

  // Some consistency checks:
  Printf("%10s: Integral over eta,IPz: %9.4f +/- %9.4f",
	 GetName(), intg, eintg);
  intg = etaDelta->IntegralAndError(1,etaDeltaIPz->GetNbinsX(),
				    lowBin, highBin,
				    eintg);
  Printf("%10s: Integral over eta:     %9.4f +/- %9.4f",
	 GetName(), intg, eintg);
  intg = delta->IntegralAndError(lowBin, highBin, eintg);
  Printf("%10s: Integral:              %9.4f +/- %9.4f",
	 GetName(), intg, eintg);
  
				    
  return true;
}

//____________________________________________________________________
Bool_t AliTrackletAODdNdeta::CentBin::EstimateBackground(Container* result,
							 Container* measCont,
							 Container* genCont,
							 Histos*    h,
							 Double_t   tailCut)
{
  if (!h || !measCont) {
    AliWarningF("No sub-histos or measured container in %s", GetName());
    return false;
  }

  if (!h->MasterFinalize(result, fIPz, tailCut)) { 
    AliWarningF("Failed to finalize %s/%s", GetName(), h->GetName());
    return false;
  }

  Container* bgCont = GetC(result, h->GetName());
  if (!bgCont) {
    AliWarningF("%s/%s didn't put a container on output",
		GetName(), h->GetName());
    return false;
  }
  const char* sub = h->GetName();
  TString shrt(Form("%c", sub[0]));
  shrt.ToUpper();
  if (genCont) shrt.Append("'");
  
  TH2* background    = 0;
  TH2* etaIPzScale   = 0;
  if (h->GetMask() == AliAODTracklet::kInjection) {
    // Get the tail ratio per eta, ipZ
    // We get the integral of the measured eta,IPz,Delta dist
    // We get the integral of the injected eta,IPz,Delta dist
    // Then we take the ratio 
    etaIPzScale = CopyH2(measCont, "etaIPzDeltaTail",
			 "etaIPzScale");
    etaIPzScale->Divide(GetH2(bgCont, "etaIPzDeltaTail"));
    etaIPzScale->SetZTitle("k_{#eta,IP_{#it{z}}}");
    etaIPzScale->SetTitle(Form("k_{%s,#eta,IP_{#it{z}}}",shrt.Data()));
    bgCont->Add(etaIPzScale);  

    // Get the tail ratio per eta
    // We get the integral of the measured eta,Delta dist
    // We get the integral of the injected eta,Detla dist
    // Then we take the ratio  
    TH1* etaScale = CopyH1(measCont, "etaDeltaTail",
			   "etaScale");
    etaScale->Divide(GetH1(bgCont, "etaDeltaTail"));
    etaScale->SetYTitle("k_{#eta}");
    etaScale->SetTitle(Form("k_{%s,#eta}", shrt.Data()));
    bgCont->Add(etaScale);

    // Get the integrated tail ratio.
    // We get the integrated tail of the measured delta dist
    // We get the integrated tail of the ianjection delta dist
    // We then get the ratio.   
    Double_t measIntg = GetD(measCont, "deltaTail",      -1);
    Double_t measIntE = GetD(measCont, "deltaTailError", -1);
    Double_t bgIntg   = GetD(bgCont,   "deltaTail",      -1);
    Double_t bgIntE   = GetD(bgCont,   "deltaTailError", -1);
    Double_t scaleE   = 0;
    Double_t scale    = RatioE(measIntg, measIntE, bgIntg, bgIntE, scaleE);
    bgCont->Add(new TParameter<double>("scale",      scale));
    bgCont->Add(new TParameter<double>("scaleError", scaleE));

    // Get the fully differential Delta distribution and scale by the
    // eta,IPz scalar.
    TH3* scaledEtaDeltaIPz = ScaleDelta(CopyH3(bgCont, "etaDeltaIPz",
					       "scaleEtaDeltaIPz"),
					etaIPzScale);
    scaledEtaDeltaIPz->SetTitle(Form("%5.3f#times%s",
				     scale, scaledEtaDeltaIPz->GetTitle()));
    scaledEtaDeltaIPz->SetDirectory(0);
    scaledEtaDeltaIPz->SetYTitle("k#timesd^{3}#it{N}/"
				 "(d#Deltad#etadIP_{#it{z}})");
    bgCont->Add(scaledEtaDeltaIPz);
#if 0
    // scale by derived scalars, rather than by taking the scaled full
    // distribution.
    TH2* scaledEtaDelta = CopyH2(bgCont, "etaDelta", "scaledEtaDelta");
    scaledEtaDelta->SetTitle(scaledEtaDeltaIPz->GetTitle());
    scaledEtaDelta->SetZTitle("k#timesd^{2}#it{N}/(d#Deltad#eta)");
    Scale(scaledEtaDelta, etaScale);
    bgCont->Add(scaledEtaDelta);
    
    TH1* scaledDelta = CopyH1(bgCont, "delta", "scaledDelta");
    scaledDelta->SetTitle(scaledEtaDeltaIPz->GetTitle());
    scaledDelta->SetYTitle("k#timesd#it{N}/d#Delta");
    Scale(scaledDelta,scale,scaleE);
    bgCont->Add(scaledDelta);
#else 
    // Make 2D projection to eta,Delta
    TH2* scaledEtaDelta = ProjectEtaDelta(scaledEtaDeltaIPz);
    scaledEtaDelta->SetName("scaledEtaDelta");
    scaledEtaDelta->SetTitle(scaledEtaDeltaIPz->GetTitle());
    scaledEtaDelta->SetYTitle("k#timesd^{2}#it{N}/(d#Deltad#eta)");
    bgCont->Add(scaledEtaDelta);
  
    // Make projection of delta 
    TH1* scaledDelta = ProjectDelta(scaledEtaDelta);
    scaledDelta->SetName("scaledDelta");
    scaledDelta->SetTitle(scaledEtaDeltaIPz->GetTitle());
    scaledDelta->SetYTitle("k#timesd#it{N}/d#Delta");
    bgCont->Add(scaledDelta);
#endif
    // Make background scaled by full tail 
    background = CopyH2(bgCont,  "etaIPz", "background");
    if (!background) AliWarningF("Didn't get background in %s", sub);
    else             background->Multiply(etaIPzScale);
    // else             Scale(background, scale, scaleE);  
  }
  else {
    // For non-injection sets (i.e., combinatorics) we cannot form the
    // scaled background until we have the real data, so we calculate
    // beta instead. 
    background = CopyH2(bgCont, "etaIPz", "background");       
    TH2* beta  = CopyH2(bgCont, "etaIPz", "beta");
    if (!background || !beta)
      AliWarningF("Didn't get background or beta in %s", sub);
    else {
      beta->Divide(GetH1(measCont, "etaIPz"));
      beta->SetTitle(Form("#beta_{%s}", shrt.Data()));
      bgCont->Add(beta);
    }
  }
  if (!background) {
    AliWarningF("Didn't get background in %s", sub);
    return false;
  }
  background->SetTitle(Form("#it{B}_{%s}", shrt.Data()));
  bgCont->Add(background);

  TH2* signal = CopyH2(measCont, "etaIPz", "signal");
  if (!signal) {
    AliWarningF("Didn't get signal in %s", sub);
    return false;
  }
  else {
    signal->SetTitle(Form("#it{S}_{%s}", shrt.Data()));
    signal->Add(background,-1);
    // Zero small bins 
    for (Int_t i = 1; i <= signal->GetNbinsX(); i++) {
      for (Int_t j = 1; j <= signal->GetNbinsX(); j++) {
	if (signal->GetBinContent(i,j)<1e-6) {
	  signal->SetBinContent(i,j,0);
	  signal->SetBinError  (i,j,0);
	}
      }
    }
    CopyAttr(background, signal);
    bgCont->Add(signal);
  }

  TH1* alpha = 0;
  if (genCont) {
    alpha = CopyH2(genCont, "etaIPz", "alpha");
    if (alpha && signal) {
      alpha->Divide(signal);
      alpha->SetTitle(Form("#alpha_{%s}", shrt.Data()));
      CopyAttr(signal, alpha);
      bgCont->Add(alpha);
    }
  }

  return true;
}

//====================================================================
AliTrackletAODdNdeta*
AliTrackletAODdNdeta::Create(Bool_t      mc,
			     const char* weights, 
			     const char* sumFile,
			     const char* resFile)
{
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("Create","No analysis manager to connect to.");
    return 0;
  }   
  AliTrackletAODdNdeta* ret = 0;
  if (mc)  {
    if (weights && weights[0] != '\0') {
      AliTrackletAODWeightedMCdNdeta* wret =
	new AliTrackletAODWeightedMCdNdeta("MidRapidityMC");
      TUrl   wurl(weights);
      TFile* wfile = TFile::Open(wurl.GetFile());
      if (!wfile) {
	::Warning("Create", "Failed to open weights file: %s",
		  wurl.GetUrl());
	return 0;
      }
      TString wnam(wurl.GetAnchor());
      if (wnam.IsNull()) wnam = "weights";
      TObject* wobj = wfile->Get(wnam);
      if (!wobj) {
	::Warning("Create", "Failed to get weights %s from file %s",
		  wnam.Data(), wfile->GetName());
	return 0;
      }
      if (!wobj->IsA()->InheritsFrom(AliTrackletWeights::Class())) {
	::Warning("Create", "Object %s from file %s not an "
		  "AliTrackletWeights but a %s",
		  wnam.Data(), wfile->GetName(), wobj->ClassName());
	return 0;
      }
      wret->SetWeights(static_cast<AliTrackletWeights*>(wobj));
      ret = wret;
    }
    else 
      ret = new AliTrackletAODMCdNdeta("MidRapidityMC");
  }
  else           ret = new AliTrackletAODdNdeta("MidRapidity");
  if (ret)       ret->Connect();

  return ret;  

}




//____________________________________________________________________
