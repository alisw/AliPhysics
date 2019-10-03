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
#include <TGraphErrors.h>
#else
// class AliAODTracklet;
class AliVEvent;
class AliTrackletBaseWeights;
class AliMultSelection;  // Auto-load 
class TClonesArray;
#endif

//====================================================================
/**
 * Task to analyse AOD tracklets for dNch/deta 
 * 
 * @ingroup pwglf_forward_tracklets
 */
class AliTrackletAODdNdeta : public AliAnalysisTaskSE,
			     public AliTrackletAODUtils
{
public:
  /**
   * Masks and vetos used by histogram sets 
   */
  enum {
    kMeasuredMask     = 0x0,
    kMeasuredVeto     = AliAODTracklet::kInjection|AliAODTracklet::kGenerated,
    kInjectedMask     = AliAODTracklet::kInjection,
    kInjectedVeto     = AliAODTracklet::kGenerated,
    kCombinatoricMask = AliAODTracklet::kCombinatorics,
    kCombinatoricVeto = 0x0,
    kDistinctMask     = (AliAODTracklet::kCombinatorics|
                         AliAODTracklet::kDistinct),
    kDistinctVeto     = 0x0,
    kPrimaryMask      = 0x0,
    kPrimaryVeto      = (AliAODTracklet::kInjection|
			 AliAODTracklet::kCombinatorics|
			 AliAODTracklet::kSecondary|
			 AliAODTracklet::kGenerated),
    kSecondaryMask    = AliAODTracklet::kSecondary,
    kSecondaryVeto    = 0x0,
    kGeneratedMask    = AliAODTracklet::kGenerated,
    kGeneratedVeto    = AliAODTracklet::kNeutral|AliAODTracklet::kSuppressed
  };
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
  enum EStatus {
    kOKEvent       = ((1<<(kAll      -1))|
		      (1<<(kEvent    -1))|
		      (1<<(kTracklets-1))),
    kOKTrigger     = (kOKEvent|(1<<(kTrigger  -1))),
    kOKIPz         = (kOKTrigger|(1<<kIP-1)),
    kOkCentrality  = (kOKIPz|(1<<(kCentrality-1)))
  };
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
   * @{ 
   * @name Centrality 
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
   * Set maximum number of tracklets for fake centrality.  Note, in
   * case of reweighting, the total number of tracklets can be
   * non-integer.  Hence, we pass a double here. 
   *
   * @param maxN Set maximum 
   */
  void SetMaxNTracklet(Double_t maxN) { fMaxNTracklet = maxN; }
  /* @} */
  //------------------------------------------------------------------
  /** 
   * @{ 
   * @name Interaction point 
   */
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
  /* @} */
  //------------------------------------------------------------------
  /** 
   * @{ 
   * @name Azimuth 
   */
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
  /* @} */
  //------------------------------------------------------------------
  /** 
   * @{ 
   * @name Pseudorapidity 
   */
  /** 
   * Set the @f$ \eta@f$ axis 
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
   * Set the @f$ \eta@f$ axis 
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
  /* @} */
  //------------------------------------------------------------------
  /** 
   * @{ 
   * @name Trackletting
   */
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
   * Set lower cut on tail of @f$\Delta@f$ distributions 
   *
   * @param x Value 
   */
  void SetTailMaximum(Double_t x=-1) { fTailMax = x; }
  /* @} */
  //------------------------------------------------------------------
  /** 
   * @{ 
   * @name Reweighting (deprecated)
   */
  /** 
   * Set weights.  This implementatio does nothing 
   * 
   * @param w Weights to use 
   */
  virtual void SetWeights(AliTrackletBaseWeights* w) {};
  /** 
   * Whether to use square root of square (product) weights 
   * 
   * @param mode If true, take square root of square (product) weights
   */
  virtual void SetWeightCalc(UChar_t mode=0) {}
  /** 
   * Set the tracklet type mask to use on weights 
   * 
   * @param mask Tracklet type mask 
   */
  virtual void SetWeightMask(UChar_t mask=0xFF) {}
  /** 
   * Inverse the weights calculated.  That is, if this option is
   * enabled, then the weight used is @f$1/w@f$ where @f$w@f$ is the
   * normal weight.
   * 
   * @param inv If true, inverse weights 
   */
  virtual void SetWeightInverse(Bool_t inv) {}
  /** 
   * Set the tracklet type veto to use on weights 
   * 
   * @param veto Tracklet type veto 
   */
  virtual void SetWeightVeto(UChar_t veto=0xFF) {}
  /* @} */
  virtual void SetTriggerEfficiency(Double_t eff) { fTriggerEff = eff; }
  virtual void SetMinEta1(Int_t n) { fMinEta1 = n; }
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
     * @param tailMax Maximum to integrate tail to
     * 
     * @return true on success
     */
    virtual Bool_t MasterFinalize(Container* parent,
				  TH1*       ipz,
				  Double_t   tailCut,
				  Double_t   tailMax) = 0;
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
	fEtaDeltaPdg(0),
	fEtaPdgIPz(0),
	fEtaPdg(0),
	fEtaPt(0),
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
	fEtaDeltaIPz(0),
	fEtaDeltaPdg(0),
	fEtaPdgIPz(0),
	fEtaPdg(0),
	fEtaPt(0),
	fColor(o.fColor),
	fStyle(o.fStyle)
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
     * @return Selection corresponds to tracklets from two particles,
     * but distinct mothers.
     */
    Bool_t IsDistinct() const
    {
      return fMask == kDistinctMask && fVeto == kDistinctVeto;
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
     * @param tailMax Maximum to integrate tail to
     * 
     * @return true on success
     */
    Bool_t MasterFinalize(Container* parent,
			  TH1*       ipz,
			  Double_t   tailCut,
			  Double_t   tailMax);

    UChar_t GetMask() const { return fMask; }
    UChar_t GetVeto() const { return fVeto; }
    /** 
     * Print information to standard output 
     * 
     * @param option Ignored 
     */
    void Print(Option_t* option="") const;
  protected:
    Bool_t ProjectEtaPdg(Container* result, Int_t nEvents);
    Bool_t ProjectEtaDeltaPdgPart(Container*     result,
				  Int_t          nEvents,
				  Double_t       tailDelta,
				  Double_t       tailMax,
				  const TString& pre,
				  const TString& tit);
    Bool_t ProjectEtaDeltaPdg(Container* result,
			      Int_t      nEvents,
			      Double_t   tailCut,
			      Double_t   tailMax);
    Bool_t ProjectEtaPdgIPz(Container*     result,
			    TH1*           ipz,
			    const TString& shn);
    const UChar_t fMask;
    const UChar_t fVeto;
    TH2*          fEtaIPz;       //!
    TH3*          fEtaDeltaIPz;  //!
    TH3*          fEtaDeltaPdg;  //!
    TH3*          fEtaPdgIPz;    //! 
    TH2*          fEtaPdg;       //! 
    TH2*          fEtaPt;        //! 
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
	fStatus(0),
	fIPz(0),
	fCent(0),
	fCentIPz(0),
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
	fStatus(0),
	fIPz(0),
	fCent(0),
	fCentIPz(0),
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
     * Create a histogram set 
     * 
     * @param name  Name of histogram set 
     * @param color Color used by histograms 
     * @param style Style used by histograms 
     * @param mask  Tracklet selection mask  
     * @param veto  Tracklet veto mask 
     * 
     * @return Newly allocated histogram set 
     */
    Histos* MakeHistos(const char* name,
		       Color_t     color,
		       Style_t     style,
		       UShort_t    mask,
		       UShort_t    veto);
    /** 
     * Initialize the bin 
     * 
     * @param parent   Parent container 
     * @param etaAxis  pseudorapidity axis to use 
     * @param ipzAxis  Interaction point Z coordinate axis 
     * @param deltaAxis @f$\Delta@f$ axis to use 
     * 
     * @return true on success 
     */
    Bool_t WorkerInit(Container*   parent,
		      const TAxis& etaAxis,
		      const TAxis& ipzAxis,
		      const TAxis& deltaAxis);
    /** 
     * Check if this is the MB "centrality" bin
     * 
     * @return True if MB bin 
     */
    Bool_t IsAllBin() const;
    /** 
     * Check if we should process this event 
     * 
     * @param status Event status 
     * @param cent   Event centrality 
     * @param ipz    Event Z-coordinate of the interaction 
     * 
     * @return true if we should process the event 
     */
    Bool_t Accept(UInt_t status, Double_t cent, Double_t ipz);
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
     * Tell bin we're done with the processing. 
     * 
     */
    void Completed();
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
     * @param tailMax Maximum to integrate tail to
     * 
     * @return true on success
     */
    Bool_t MasterFinalize(Container* parent,
			  TH1*       ipz,
			  Double_t   tailCut,
			  Double_t   tailMax);
    /** 
     * Estimate the background a given histogram set 
     * 
     * @param result   Output container 
     * @param measCont The measured results 
     * @param genCont  The generator results (if applicable)
     * @param h        The histogram container 
     * @param tailCut  Cut on the tail distribution 
     * @param tailMax  Maximum to integrate tail to
     * 
     * @return true on success 
     */
    Bool_t EstimateBackground(Container* result,
			      Container* measCont,
			      Container* genCont,
			      Histos*    h,
			      Double_t   tailCut,
			      Double_t   tailMax);
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
    Container*  fSubs;
    Double_t    fLow;
    Double_t    fHigh;
    TH1*        fStatus;  //! 
    TH1*        fIPz;     //! 
    TH1*        fCent;    //!
    TProfile*   fCentIPz; //! 
    Histos*     fMeasured; 
    Histos*     fInjection;

    ClassDef(CentBin,2);
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
  /** 
   * Create and initialize a centrality bin 
   * 
   * @param bin        Bin number 
   * @param c1         Least centrality 
   * @param c2         Largest centrality 
   * @param existing   Possible existing container to initialize from 
   * @param deltaAxis  Axis for @f$ \Delta @f$ 
   *
   * @return Pointer to bin object or null
   */
  virtual CentBin* InitCentBin(Int_t        bin,
			       Float_t      c1,
			       Float_t      c2,
			       Container*   existing,
			       const TAxis& deltaAxis);
  /* @} */
  
  // -----------------------------------------------------------------
  /** 
   * @{ 
   * @name Event inspection 
   */
  virtual const char* GetBranchName() const { return "AliAODTracklets"; }
  /** 
   * Calculate bit flag corresponding to the bin number 
   * 
   * @param bin Bin number 
   * 
   * @return Bit flag 
   */
  static UInt_t Bin2Flag(UInt_t bin);
  /** 
   * Check event 
   * 
   * @param cent        On return, the event centrality 
   * @param ip          On return, the event interaction point 
   * @param tracklets   On return, the list of tracklets 
   * 
   * @return Bit mask of completed requirements
   */
  UInt_t CheckEvent(Double_t&          cent,
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
   * @param event      Event 
   * @param value      On return, estimator value 
   * @param nTracklets Number of tracklets for benchmarking centrality 
   * @param nCl0       On return, number of clusters on layer 0
   * @param nCl1       On return, number of clusters on layer 1
   * 
   * @return Centrality percentile, or negative number in case of
   * problems
   */
  Double_t FindMultCentrality(AliVEvent* event,
			      Double_t&  value,
			      Int_t&     nTracklets,
			      Int_t&     nCl0,
			      Int_t&     nCl1);
  Double_t FindMultCentrality(AliVEvent* event,
			      Int_t&     nTracklets);
  /** 
   * Find the centrality of the event. This looks for the
   * AliCentrality structure.
   * 
   * @param event Event 
   * @param nTracklets Number of tracklets for benchmarking centrality 
   * 
   * @return Centrality percentile, or negative number in case of
   * problems
   */
  Double_t FindFakeCentrality(AliVEvent* event, Int_t& nTracklets);  
  /** 
   * Find the fake centrality of the event. This look at the number of
   * tracklets
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
   * @param ipz      Z coordinate of IP
   * 
   * @return The weight - in this class always 1
   */
  virtual Double_t LookupWeight(AliAODTracklet* tracklet,
				Double_t        cent,
				Double_t        ipz);
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
   * @param status      The event status flags 
   * @param cent        Event centrality 
   * @param ip          Event interaction point 
   * @param tracklets   List of tracklets 
   */
  void ProcessEvent(UInt_t            status,
		    Double_t          cent,
		    const AliVVertex* ip,
		    TClonesArray*     tracklets);
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

  /** 
   * Make status histogram 
   * 
   * @param container Container to add to 
   * 
   * @return The status histogram created 
   */
  static TH1* MakeStatus(Container* container);
  
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
  /** Average good number of tracklets (weighed) vs total (unweighed) */
  TProfile* fNBareVsGood;
  /** Average fake number of tracklets (weighed) vs total (unweighed) */
  TProfile* fNBareVsFake;
  /** Average good number of tracklets (weighed) vs total (weighed) */
  TProfile* fNTrackletVsGood;
  /** Average fake number of tracklets (weighed) vs total (weighed) */
  TProfile* fNTrackletVsFake;
  /** Average good number of tracklets (weighed) vs generated tracks */
  TProfile* fNGeneratedVsGood;
  /** Average fake number of tracklets (weighed) vs generated tracks */
  TProfile* fNGeneratedVsFake;  
  /** Histogram of centrality, nTracklets correlation */
  TProfile* fCentTracklets; //!
  /** Histogram of centrality vs average mult estimator */
  TProfile* fCentEst; //!
  /** Histogram of correlation of clusters and tracklets */
  TProfile2D* fCl0Cl1Tracklets;
  /** Centrality method to use */
  TString    fCentMethod;
  /** Cached index of centrality estimator */
  Int_t      fCentIdx;
  /** Cached index of tracklet estimator */
  Int_t      fTrkIdx;
  /** Cached index of layer 0 clusters estimator */
  Int_t      fCl0Idx;
  /** Cached index of layer 1 clusters estimator */
  Int_t      fCl1Idx;
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
  /** Largest value of @f$\Delta@f$ considered background tail */
  Double_t   fTailMax;
  /** Shift @f$\delta_{\phi}@f$ of @f$\Delta\phi@f$ */
  Double_t   fDPhiShift;
  /** Signal cut on @f$\Delta\phi-\delta_{\phi}@f$ */
  Double_t   fShiftedDPhiCut;
  /** Signal cut on @f$\Delta@f$ */
  Double_t   fDeltaCut;
  /** The absolute minimum centrality to consider  - for MC with poor match*/
  Double_t   fAbsMinCent;
  /** Maximum number of tracklets for fake centrality */
  Double_t   fMaxNTracklet;
  /** External set trigger efficiency */
  Double_t   fTriggerEff;
  /** Least number of tracklets in |eta|<1 */
  Int_t      fMinEta1;
  
  ClassDef(AliTrackletAODdNdeta,2); 
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
	fDistinct(0),
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
    Histos*    fDistinct;
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
      fEtaWeight(0),
      fWeightCorr(0)
  {}
  /** 
   * Named - user - constructor
   */
  AliTrackletAODWeightedMCdNdeta(const char* name)
    : AliTrackletAODMCdNdeta(name),
      fWeights(0),
      fEtaWeight(0),
      fWeightCorr(0)
  {}
  /**
   * Copy constructor 
   *
   * @param o Object to copy from 
   */
  AliTrackletAODWeightedMCdNdeta(const AliTrackletAODWeightedMCdNdeta& o)
    : AliTrackletAODMCdNdeta(o),
      fWeights(0),
      fEtaWeight(0),
      fWeightCorr(0)
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
  void SetWeights(AliTrackletBaseWeights* w) { fWeights = w; }
  /** 
   * Whether to use square root of square (product) weights 
   * 
   * @param mode If true, take square root of square (product) weights
   */
  void SetWeightCalc(UChar_t mode=0) {
    Info("SetWeightCalc", "Use=%d", mode);
    if (fWeights) fWeights->SetCalc(mode);}
  /** 
   * Set the tracklet type mask to use on weights 
   * 
   * @param mask Tracklet type mask 
   */
  void SetWeightMask(UChar_t mask=0xFF) {
    Info("SetWeightMask", "mask=0x%x", mask);
    if (fWeights) fWeights->SetMask(mask); }
  /** 
   * Set the tracklet type veto to use on weights 
   * 
   * @param veto Tracklet type veto 
   */
  void SetWeightVeto(UChar_t veto=0x0) {
    Info("SetWeightVeto", "veto=0x%x", veto);
    if (fWeights) fWeights->SetVeto(veto); }
  /** 
   * Inverse the weights calculated.  That is, if this option is
   * enabled, then the weight used is @f$1/w@f$ where @f$w@f$ is the
   * normal weight.
   * 
   * @param inv If true, inverse weights 
   */
  void SetWeightInverse(Bool_t inv) { if (fWeights) fWeights->SetInverse(inv); }
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
   * @param ipz      Z coordinate of IP
   * 
   * @return The weight - in this class always 1
   */
  virtual Double_t LookupWeight(AliAODTracklet* tracklet,
				Double_t        cent,
				Double_t        ipz);
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
  AliTrackletBaseWeights* fWeights;
  TProfile2D*             fEtaWeight;  //! 
  TH2*                    fWeightCorr;    //! 
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
    fNBareVsGood(0),
    fNBareVsFake(0),
    fNTrackletVsGood(0),
    fNTrackletVsFake(0),
    fNGeneratedVsGood(0),
    fNGeneratedVsFake(0),    
    fCentTracklets(0),
    fCentEst(0),
    fCl0Cl1Tracklets(0),
    fCentMethod(""),
    fCentIdx(-1),
    fTrkIdx(-1),
    fCl0Idx(-1),
    fCl1Idx(-1),
    fCentAxis(1,0,0),
    fIPzAxis(1,0,0),
    fEtaAxis(1,0,0),
    fPhiAxis(1,0,0),
    fMaxDelta(0),
    fTailDelta(0),
    fTailMax(-1),
    fDPhiShift(0),
    fShiftedDPhiCut(0),
    fDeltaCut(0),
    fAbsMinCent(-1),
    fMaxNTracklet(6000),
    fTriggerEff(0),
    fMinEta1(0)
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
    fNBareVsGood(0),
    fNBareVsFake(0),
    fNTrackletVsGood(0),
    fNTrackletVsFake(0),
    fNGeneratedVsGood(0),
    fNGeneratedVsFake(0),    
    fCentTracklets(0),
    fCentEst(0),
    fCl0Cl1Tracklets(0),
    fCentMethod("V0M"),
    fCentIdx(-1),
    fTrkIdx(-1),
    fCl0Idx(-1),
    fCl1Idx(-1),
    fCentAxis(10,0,100),
    fIPzAxis(30,-15,+15),
    fEtaAxis(16,-2,+2),
    fPhiAxis(100,0,TMath::TwoPi()),
    fMaxDelta(25),
    fTailDelta(5),
    fTailMax(-1),
    fDPhiShift(0.0045),
    fShiftedDPhiCut(-1),
    fDeltaCut(1.5),
    fAbsMinCent(-1),
    fMaxNTracklet(6000),
    fTriggerEff(0),
    fMinEta1(0)
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
    fNBareVsGood(0),
    fNBareVsFake(0),
    fNTrackletVsGood(0),
    fNTrackletVsFake(0),
    fNGeneratedVsGood(0),
    fNGeneratedVsFake(0),    
    fCentTracklets(0),
    fCentEst(0),
    fCl0Cl1Tracklets(0),
    fCentMethod(o.fCentMethod),
    fCentIdx(o.fCentIdx),
    fTrkIdx(o.fTrkIdx),
    fCl0Idx(o.fCl0Idx),
    fCl1Idx(o.fCl1Idx),
    fCentAxis(o.fCentAxis),
    fIPzAxis(o.fIPzAxis),
    fEtaAxis(o.fEtaAxis),
    fPhiAxis(o.fPhiAxis),
    fMaxDelta(o.fMaxDelta),
    fTailDelta(o.fTailDelta),
    fTailMax(o.fTailMax),
    fDPhiShift(o.fDPhiShift),
    fShiftedDPhiCut(o.fShiftedDPhiCut),
    fDeltaCut(o.fDeltaCut),
    fAbsMinCent(o.fAbsMinCent),
    fMaxNTracklet(o.fMaxNTracklet),
    fTriggerEff(o.fTriggerEff),
    fMinEta1(o.fMinEta1)
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
  fCentIdx        = o.fCentIdx;
  fTrkIdx         = o.fTrkIdx;
  fCl0Idx         = o.fCl0Idx;
  fCl1Idx         = o.fCl1Idx;
  fCentAxis       = o.fCentAxis;
  fIPzAxis        = o.fIPzAxis;
  fEtaAxis        = o.fEtaAxis;
  fPhiAxis        = o.fPhiAxis;
  fMaxDelta       = o.fMaxDelta;
  fTailDelta      = o.fTailDelta;
  fTailMax        = o.fTailMax;
  fDPhiShift      = o.fDPhiShift;
  fShiftedDPhiCut = o.fShiftedDPhiCut;
  fDeltaCut       = o.fDeltaCut;
  fAbsMinCent     = o.fAbsMinCent;
  fMaxNTracklet   = o.fMaxNTracklet;
  fTriggerEff     = o.fTriggerEff;
  fMinEta1        = o.fMinEta1;
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
  Double_t tailMax = fTailMax;
  if (tailMax < 0) tailMax = fMaxDelta;
  
  Printf("%s: %s", ClassName(), GetName());
  Printf(" %22s: 0x%08x", "Off-line trigger mask", fOfflineTriggerMask);
  Printf(" %22s: %f",   "Delta phi shift",	   fDPhiShift);
  Printf(" %22s: %f",   "Shifted Delta phi cut",   shiftedDPhiCut);
  Printf(" %22s: %f",   "Delta cut",	           fDeltaCut);
  Printf(" %22s: %f",   "max Delta",	           fMaxDelta);
  Printf(" %22s: %f",   "tail Delta",	           fTailDelta);
  Printf(" %22s: %f",   "tail maximum",	           tailMax);
  Printf(" %22s: %f%%", "Absolute least c",        fAbsMinCent);
  Printf(" %22s: %f",   "Max(Ntracklet)",          fMaxNTracklet);
  Printf(" %22s: %f",   "Trigger eff)",            fTriggerEff);
  Printf(" %22s: %d",   "Min(Ntracklet)",          fMinEta1);
  Printf(" %22s: %s (%d)", "Centrality estimator", fCentMethod.Data(),fCentIdx);
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
  Printf("   Mask:          0x%02x (%s)", fMask, cMask);
  Printf("   Veto:          0x%02x (%s)", fVeto, cVeto);
  Printf("   Delta:         %s", fEtaDeltaIPz ? "yes" : "no");
  Printf("   Delta per PDG: %s", fEtaDeltaPdg ? "yes" : "no");
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
  if (fTailMax        < 0) fTailMax        = fMaxDelta;
  
  fContainer = new Container;
  fContainer->SetName(Form("%sSums", GetName()));
  fContainer->SetOwner();

  Bool_t save = TH1::GetDefaultSumw2();
  TH1::SetDefaultSumw2();
  fIPz    = Make1D(fContainer, "ipz",  "", kMagenta+2, 20, fIPzAxis);
  fCent   = Make1D(fContainer, "cent", "", kMagenta+2, 20, fCentAxis);
  fEtaPhi = Make2D(fContainer, "etaPhi","",kMagenta+2, 20, fEtaAxis,fPhiAxis);

  TAxis trackletAxis(1000, 0, 10000);
  FixAxis(trackletAxis, "#it{N}_{tracklets}");
  fNBareVsGood  = Make1P(fContainer, "nBareVsGood", "Good",
			 kGreen+2, 20, trackletAxis);
  fNBareVsFake  = Make1P(fContainer, "nBareVsFake", "Fake",
			 kRed+2, 20, trackletAxis);
  fNBareVsGood->SetYTitle("#LT#it{N}_{tracklet}#GT");
  fNBareVsFake->SetYTitle("#LT#it{N}_{tracklet}#GT");
  fNBareVsGood->SetStats(0);
  fNBareVsFake->SetStats(0);

  fNTrackletVsGood  = Make1P(fContainer, "nTrackletVsGood", "Good",
			     kGreen+2, 20, trackletAxis);
  fNTrackletVsFake  = Make1P(fContainer, "nTrackletVsFake", "Fake",
			     kRed+2, 20, trackletAxis);
  fNTrackletVsGood->SetYTitle("#LT#it{N}_{tracklet}#GT");
  fNTrackletVsFake->SetYTitle("#LT#it{N}_{tracklet}#GT");
  fNTrackletVsGood->SetStats(0);
  fNTrackletVsFake->SetStats(0);
  
  TAxis generatedAxis(1000, 0, 15000);
  FixAxis(generatedAxis, "#it{N}_{generated,|#eta|<2}");
  fNGeneratedVsGood = Make1P(fContainer, "nGeneratedVsGood", "Good",
			     kGreen+2, 24, generatedAxis);
  fNGeneratedVsFake = Make1P(fContainer, "nGeneratedVsFake", "Fake",
			     kRed+2, 24, generatedAxis);
  fNGeneratedVsGood->SetYTitle("#LT#it{N}_{tracklet}#GT");
  fNGeneratedVsFake->SetYTitle("#LT#it{N}_{tracklet}#GT");
  fNGeneratedVsGood->SetStats(0);
  fNGeneratedVsFake->SetStats(0);

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
  fCentEst = Make1P(fContainer,"centEstimator","",kMagenta+2, 20, fCentAxis);
  fCentEst->SetYTitle(Form("#LT%s#GT",fCentMethod.Data()));

  TAxis clAxis(150, 0, 15000);
  FixAxis(clAxis, "#it{N}_{cluster}");
  fCl0Cl1Tracklets = Make2P(fContainer, "cl0cl1Tracklets","",kMagenta+2, 20,
			    clAxis, clAxis);
  fCl0Cl1Tracklets->SetXTitle(Form("Layer 0 %s", clAxis.GetTitle()));
  fCl0Cl1Tracklets->SetYTitle(Form("Layer 1 %s", clAxis.GetTitle()));
  fCl0Cl1Tracklets->SetZTitle("#LT#it{N}_{tracklet}#GT");

  TH1::SetDefaultSumw2(save);
  
  fStatus = MakeStatus(fContainer);

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
  params->Add(new DP("TailMax",        fTailMax,        'f'));
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
  fWeights->SetDebug(fDebug);

  TAxis wAxis(100,0,10);
  FixAxis(wAxis,"#it{w}_{i}");
  fEtaWeight = Make2P(fContainer, "etaWeight", "#LTw#GT", kYellow+2, 24,
		      fEtaAxis, fCentAxis);
  fWeightCorr = Make2D(fContainer, "weightCorr", "w_{1} vs w_{2}",
		       kCyan+2, 24, wAxis, wAxis);
  fWeights->Store(fContainer);
  return true;
}
//____________________________________________________________________
AliTrackletAODdNdeta::CentBin*
AliTrackletAODdNdeta::InitCentBin(Int_t        bin,
				  Float_t      c1,
				  Float_t      c2,
				  Container*   existing,
				  const TAxis& deltaAxis)
{
  CentBin* ret  = MakeCentBin(c1, c2);
  Bool_t   ok   = true;
  if (!existing) 
    ok = ret->WorkerInit(fContainer,fEtaAxis,fIPzAxis,deltaAxis);
  else
    ok = ret->FinalizeInit(existing);
  if (!ok) {
    AliWarningF("Failed to initialize bin %s", ret->GetName());
    delete ret;
    return 0;
  }
  ret->SetDebug(DebugLevel());
  fCentBins->AddAt(ret, bin);
  return ret;

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

  TAxis    deltaAxis (Int_t(5*fMaxDelta),0,fMaxDelta);
  FixAxis(deltaAxis,
	  "#Delta=[(#Delta#phi-#delta#phi)/#sigma_{#phi}]^{2}+"
	  "[#Delta#thetasin^{-2}(#theta)/#sigma_{#theta}]^{2}");

  // Add min-bias bin
  if (!InitCentBin(0,0,0,existing,deltaAxis)) return false;

  // Add other bins
  Int_t nCentBins = fCentAxis.GetNbins();
  for (Int_t i = 1; i <= nCentBins; i++) {
    Float_t  c1 = fCentAxis.GetBinLowEdge(i);
    Float_t  c2 = fCentAxis.GetBinUpEdge(i);
    if (!InitCentBin(i, c1, c2, existing, deltaAxis)) return false;
  }
  if (!InitCentBin(nCentBins+1, 0, 100, existing, deltaAxis)) return false;
  return true;
}

//____________________________________________________________________
AliTrackletAODdNdeta::CentBin::CentBin(Double_t c1, Double_t c2)
  : Sub(""),
    fSubs(0),
    fLow(c1),
    fHigh(c2),
    fStatus(0),
    fIPz(0),
    fCent(0),
    fCentIPz(0),
    fMeasured(0),
    fInjection(0)
{
  if (c1 >= c2)
    fName = "all";
  else 
    fName = AliTrackletAODUtils::CentName(c1, c2);
  fMeasured  = MakeHistos("measured", kRed+2, 20,
			  kMeasuredMask, // No requirements, just veto
			  kMeasuredVeto);
  fInjection = MakeHistos("injected", kOrange+2, 21,
			  kInjectedMask,
			  kInjectedVeto);
  fSubs = new Container;
  fSubs->SetOwner(true);
  fSubs->Add(fMeasured);
  fSubs->Add(fInjection);
}

//____________________________________________________________________
AliTrackletAODMCdNdeta::CentBin::CentBin(Double_t c1, Double_t c2)
  : AliTrackletAODdNdeta::CentBin(c1, c2),
    fCombinatorics(0),
    fDistinct(0),
    fPrimaries(0),
    fSecondaries(0),
    fGenerated(0)
{
  // Combinatorics is everything that has the combinatorics bit
  // Primaries all with simulated bit, but not secondary or
  // combinatorics bit.
  // Secondaries are all those with the secondary bit set 
  fCombinatorics  = MakeHistos("combinatorics", kMagenta+2, 30,
			       kCombinatoricMask,
			       kCombinatoricVeto);
  fDistinct       = MakeHistos("distinct", kCyan+2, 27,
			       kDistinctMask,
			       kDistinctVeto);
  fPrimaries      = MakeHistos("primaries", kGreen+2, 26, 
			       kPrimaryMask,
			       kPrimaryVeto);
  fSecondaries    = MakeHistos("secondaries", kBlue+2, 32, 
			       kSecondaryMask,
			       kSecondaryVeto);
  fGenerated      = MakeHistos("generated", kGray+1, 28, 
			       kGeneratedMask,
			       kGeneratedVeto);
  fMeasured->fStyle = 24;
  fInjection->fStyle = 25;
  fSubs->Add(fCombinatorics);
  fSubs->Add(fDistinct);
  fSubs->Add(fPrimaries);
  fSubs->Add(fSecondaries);
  fSubs->AddAfter(fMeasured, fGenerated);
}

//____________________________________________________________________
AliTrackletAODdNdeta::Histos*
AliTrackletAODdNdeta::CentBin::MakeHistos(const char* name,
					  Color_t     color,
					  Style_t     style,
					  UShort_t    mask,
					  UShort_t    veto)
{
  return new Histos(name, color, style, mask, veto);
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
  fStatus    = MakeStatus(fContainer);
  fCent      = Make1D(fContainer,"cent","Centrality [%]",
		      kMagenta+2,20,centAxis);
  fIPz       = Make1D(fContainer,"ipz","IP_{#it{z}} [cm]",kRed+2,20,ipzAxis);
  fCentIPz   = Make1P(fContainer,"centIpz","#LTc#GT vs IP_{#it{z}}",kPink+2,
		      20,ipzAxis);
  fCentIPz->SetYTitle("#LTc#GT");
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
  if (!IsPrimary() && !IsSecondary()) 
    fEtaIPz   = Make2D(fContainer, "etaIPz", shrt.Data(),
		       kRed+2, 20, etaAxis, ipzAxis);
  // Always make eta vs Delta distribution, except for MC truth 
  if (!IsGenerated()) 
    fEtaDeltaIPz = Make3D(fContainer, "etaDeltaIPz",
			  Form("#Delta_{%s}",shrt.Data()), 
			  kBlue+2, 21, etaAxis, deltaAxis, ipzAxis);
  if (!IsGenerated() &&
      // !IsPrimary() &&
      // !IsSecondary() &&
      !IsInjected() &&
      fStyle != 20) // Last condition to not make this for real data
    fEtaPdgIPz = Make3D(fContainer, "etaPdgIPz",
			"Parent particle type",
			kGreen+2, 22, etaAxis, PdgAxis(), ipzAxis);

  if (IsGenerated()) {
    fEtaPdg = Make2D(fContainer, "etaPdg",
		     "Primary particle type",
		     kYellow+2, 30, etaAxis, PdgAxis());
    TAxis ptAxis(100, 0, 5);
    ptAxis.SetTitle("#it{p}_{T}");
    FixAxis(ptAxis);
    fEtaPt = Make2D(fContainer, "etaPt",
		    "Primary transverse momentum",
		    kCyan+2, 30, etaAxis, ptAxis);
  }
  if (IsPrimary() || IsSecondary() || IsCombinatoric()) {
    fEtaDeltaPdg = Make3D(fContainer, "etaDeltaPdg",
			  "#Delta by primary particle type",
			  kMagenta, 22, etaAxis, deltaAxis, PdgAxis());
  }
  
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
  if (fEtaDeltaPdg) {
    fEtaDeltaPdg->SetMarkerStyle(m);
    fEtaDeltaPdg->SetMarkerColor(c);
    fEtaDeltaPdg->SetLineColor(c);
    fEtaDeltaPdg->SetFillColor(c);
  }
  if (fEtaPdgIPz) {
    fEtaPdgIPz->SetMarkerStyle(m);
    fEtaPdgIPz->SetMarkerColor(c);
    fEtaPdgIPz->SetLineColor(c);
    fEtaPdgIPz->SetFillColor(c);
  }    
  if (fEtaPdg) {
    fEtaPdg->SetMarkerStyle(m);
    fEtaPdg->SetMarkerColor(c);
    fEtaPdg->SetLineColor(c);
    fEtaPdg->SetFillColor(c);
  }    
  if (fEtaPt) {
    fEtaPt->SetMarkerStyle(m);
    fEtaPt->SetMarkerColor(c);
    fEtaPt->SetLineColor(c);
    fEtaPt->SetFillColor(c);
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
  UInt_t            status    = CheckEvent(cent, ip, tracklets);
  if ((status & kOKEvent) != kOKEvent) {
    AliWarningF("Event didn't pass %f, %p, %p", cent, ip, tracklets);
    return;
  }
  if (DebugLevel() > 0) Printf("Got centrality=%f ipZ=%f %d tracklets",
			       cent,
			       (ip ? ip->GetZ() : -1000),
			       tracklets->GetEntriesFast());
  ProcessEvent(status, cent, ip, tracklets);

  PostData(1,fContainer);
  fStatus->Fill(kCompleted);
}

//____________________________________________________________________
UInt_t AliTrackletAODdNdeta::Bin2Flag(UInt_t bin)
{
  return (1 << (bin-1));
}

//____________________________________________________________________
TH1* AliTrackletAODdNdeta::MakeStatus(Container* container)
{
  TH1* status = new TH1F("status", "Status of task",
			 kCompleted, .5, kCompleted+.5);
  status->SetMarkerSize(2);
  status->SetMarkerColor(kMagenta+2);
  status->SetLineColor(kMagenta+2);
  status->SetFillColor(kMagenta+2);
  status->SetFillStyle(1001);
  status->SetBarOffset(0.1);
  status->SetBarWidth(0.4);
  status->SetDirectory(0);
  status->SetStats(0);
  status->SetXTitle("Event have");
  status->SetYTitle("# Events");
  status->GetXaxis()->SetBinLabel(kAll,            "Been seen");
  status->GetXaxis()->SetBinLabel(kEvent,          "Event data");
  status->GetXaxis()->SetBinLabel(kTracklets,      "Tracklets");
  status->GetXaxis()->SetBinLabel(kTrigger,        "Trigger");
  status->GetXaxis()->SetBinLabel(kIP,             "IP");
  status->GetXaxis()->SetBinLabel(kCentrality,     "Centrality");
  status->GetXaxis()->SetBinLabel(kCompleted,      "Completed");
  container->Add(status);
  return status;
}

//____________________________________________________________________
UInt_t AliTrackletAODdNdeta::CheckEvent(Double_t&          cent,
					const AliVVertex*& ip,
					TClonesArray*&     tracklets)
{
  UInt_t ret = Bin2Flag(kAll);
  // Count all events 
  fStatus->Fill(kAll);

  // Check for event 
  AliVEvent* event = InputEvent();
  if (!event) {
    AliWarning("No event");
    return ret;
  }
  ret |= Bin2Flag(kEvent);
  fStatus->Fill(kEvent);

  // Check if we have the tracklets 
  tracklets = FindTracklets(event);
  if (!tracklets) return ret;
  ret |= Bin2Flag(kTracklets);
  fStatus->Fill(kTracklets);

  // Check for least number of tracklets in |Eta|<1
  if (fMinEta1 > 0) {
    Int_t nOK = 0;
    AliAODTracklet* tracklet   = 0;
    TIter           nextTracklet(tracklets);
    while ((tracklet = static_cast<AliAODTracklet*>(nextTracklet()))) {
      UShort_t signal = CheckTracklet(tracklet);
      if (signal && tracklet->IsMeasured() && TMath::Abs(tracklet->GetEta()) < 1) nOK++;
    }
    if (nOK < fMinEta1) return ret;    
  }
  // Check if event was triggered 
  Bool_t trg = FindTrigger();
  if (!trg) return ret;
  ret |= Bin2Flag(kTrigger);
  fStatus->Fill(kTrigger);
    
  // Check the interaction point 
  ip = FindIP(event);
  if (!ip) return ret;
  ret |= Bin2Flag(kIP);
  fStatus->Fill(kIP);

  // Check the centrality
  Int_t nTracklets = 0;
  cent = FindCentrality(event, nTracklets);
  // Do not fail on missing centrality 
  if (cent >= 0) {
    ret |= Bin2Flag(kCentrality);
    fStatus->Fill(kCentrality);
  }
  
  fIPz->Fill(ip->GetZ());
  fCent->Fill(cent);
  fCentTracklets->Fill(cent, nTracklets); 
  
  return ret;
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
						  Double_t&  value,
						  Int_t&     nTracklets,
						  Int_t&     nCl0,
						  Int_t&     nCl1)
{
  AliMultSelection* cent =
    static_cast<AliMultSelection*>(event->FindListObject("MultSelection"));
  if (!cent) {
    AliWarning("No centrality in event");
    event->GetList()->Print();
    return -1;
  }
  if (fCentIdx < 0) {
    TString estName(fCentMethod);
    TString trkName("SPDTracklets");
    TString cl0Name("CL0");
    TString cl1Name("CL1");
    if (estName.EqualTo("fake", TString::kIgnoreCase)) estName = "V0M";
    for (Int_t i = 0; i < cent->GetNEstimators(); i++) {
      AliMultEstimator* e = cent->GetEstimator(i);
      if      (!e)                                continue;
      if      (estName.EqualTo(e->GetName()))     fCentIdx = i;
      if      (trkName.EqualTo(e->GetName()))     fTrkIdx  = i;
      else if (cl0Name.EqualTo(e->GetName()))     fCl0Idx  = i;
      else if (cl1Name.EqualTo(e->GetName()))     fCl1Idx  = i;
    }
    if (fTrkIdx < 0) 
      AliWarning("Index for tracklet estimator not found!");
    if (fCl0Idx < 0) 
      AliWarning("Index for layer 0 cluster estimator not found!"); 
    if (fCl1Idx < 0) 
      AliWarning("Index for layer 1 cluster estimator not found!");
 }
  if (fCentIdx < 0) {
    AliWarningF("Centrality estimator %s not found", fCentMethod.Data());
    return -1;
  }

  // Use cached index to look up the estimator 
  AliMultEstimator* estCent = cent->GetEstimator(fCentIdx);
  if (!estCent) {
    AliWarningF("Centrality estimator %s not available for this event",
		fCentMethod.Data());
    return -1;
  }
  value = estCent->GetValue();

  // Now look up tracklet value using cached index 
  if (fTrkIdx >= 0) {
    AliMultEstimator* est        = cent->GetEstimator(fTrkIdx);
    if (est)          nTracklets = est->GetValue();
  }
  if (fCl0Idx >= 0) {
    AliMultEstimator* est  = cent->GetEstimator(fCl0Idx);
    if (est)          nCl0 = est->GetValue();
  }
  if (fCl1Idx >= 0) {
    AliMultEstimator* est  = cent->GetEstimator(fCl1Idx);
    if (est)          nCl1 = est->GetValue();
  }
  
  return estCent->GetPercentile();
}
//____________________________________________________________________
Double_t AliTrackletAODdNdeta::FindMultCentrality(AliVEvent* event,
						  Int_t& nTracklets)
{  
  Int_t nCl0 = 0, nCl1 = 0;
  Double_t value = 0;
  Double_t cent = FindMultCentrality(event, value, nTracklets, nCl0, nCl1);
  fCentEst        ->Fill(cent, value);
  fCl0Cl1Tracklets->Fill(nCl0, nCl1, nTracklets);
  return cent;
}
//____________________________________________________________________
Double_t AliTrackletAODdNdeta::FindFakeCentrality(AliVEvent* event,
						  Int_t& nTracklets)

{
  Int_t nCl0 = 0, nCl1 = 0;
  Double_t value = 0;
  FindMultCentrality(event, value, nTracklets, nCl0, nCl1);

  AliAODTracklet* tracklet   = 0;
  TClonesArray*   tracklets  = FindTracklets(event);
  TIter           nextTracklet(tracklets);
  Int_t           nTracklet  = 0;
  while ((tracklet = static_cast<AliAODTracklet*>(nextTracklet()))) {
    if (tracklet->IsMeasured())
      nTracklet += LookupWeight(tracklet,0,0);
  }
  nTracklets = nTracklet;
  fCl0Cl1Tracklets->Fill(nCl0, nCl1, nTracklets);

  if (nTracklet > fMaxNTracklet) return -1;

  Double_t c = 100*(1-Float_t(nTracklet)/fMaxNTracklet);

  fCentEst->Fill(c, nTracklet);
  return c;
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
  
  if (fCentMethod.EqualTo("fake",TString::kIgnoreCase)) 
    centPer = FindFakeCentrality(event, nTracklets);  
  else if (fCentMethod.BeginsWith("OLD"))
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
					    Double_t        cent,
					    Double_t        ipz)
{
  return 1;
}
//____________________________________________________________________
Double_t AliTrackletAODWeightedMCdNdeta::LookupWeight(AliAODTracklet* tracklet,
						      Double_t        cent,
						      Double_t        ipz)
{
  // We don't check for weights, as we must have them to come this far 
  // if (!fWeights) {
  // AliWarning("No weights defined");
  // return 1;
  // }
  Double_t w = fWeights->LookupWeight(tracklet, cent, ipz, fWeightCorr);
  if (tracklet->IsMeasured()) 
    fEtaWeight->Fill(tracklet->GetEta(), cent, w);

  if (fDebug > 4) {
    Bool_t ok=true;
    switch (TMath::Abs(tracklet->GetParentPdg())) {
    case 310:  printf("K0s   "); break;
    case 321:  printf("K+-   "); break;
    case 3122: printf("Lambda"); break;
    case 3112: printf("Sigma-"); break;
      // case 3212: printf("Sigma0"); break;// Old, wrong code
    case 3222: printf("Sigma+"); break;
    case 3312: printf("Xi-   "); break; 
      // case 3322: printf("Xi0   "); break;// Old, wrong code 
    default:   ok = false;        break;
    }
    if (ok) {
      switch (TMath::Abs(tracklet->GetParentPdg(true))) {
      case 0:    printf("       "); break;
      case 310:  printf("-K0s   "); break;
      case 321:  printf("-K+-   "); break;
      case 3122: printf("-Lambda"); break;
      case 3112: printf("-Sigma-"); break;
	// case 3212: printf("-Sigma0"); break;// Old, wrong code
      case 3222: printf("-Sigma+"); break;
      case 3312: printf("-Xi-   "); break; 
	// case 3322: printf("-Xi0   "); break;// Old, wrong code
      default:   printf("-Other ");
      }
      printf(": %f", w); tracklet->Print();
    }
  }
    
  if (fDebug > 3) {
    printf("Looking up weight of tracklet -> %f ", w);
    tracklet->Print();
  }
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
void AliTrackletAODdNdeta::ProcessEvent(UInt_t            status,
					Double_t          cent,
					const AliVVertex* ip,
					TClonesArray*     tracklets)
{
  // Figure out which centrality bins to fill 
  Double_t ipz        = (ip ? ip->GetZ() : -1000);
  Int_t    nAcc = 0;
  TIter    nextAcc(fCentBins);
  CentBin* bin = 0;
  TList    toRun;
  while ((bin = static_cast<CentBin*>(nextAcc()))) {
    // Not in range for this bin
    if (!bin->Accept(status, cent, ipz)) continue; 
    toRun.Add(bin);
    nAcc++;
  }
  // If we have no centrality bins  to fill, we return immediately 
  if (nAcc <= 0) return;

  Double_t        nBare      = 0;
  Double_t        nMeasured  = 0;
  Double_t        nGenerated = 0;
  Double_t        nGood      = 0;
  Double_t        nFake      = 0;
  AliAODTracklet* tracklet   = 0;
  TIter           nextTracklet(tracklets);
  while ((tracklet = static_cast<AliAODTracklet*>(nextTracklet()))) {
    Double_t weight = LookupWeight(tracklet, cent, ipz);
    UShort_t signal = CheckTracklet(tracklet);
    if (signal) fEtaPhi->Fill(tracklet->GetEta(), tracklet->GetPhi());
    if (tracklet->IsMeasured() && TMath::Abs(tracklet->GetEta()) < 0.7) {
      nMeasured += weight;
      nBare     ++;
      if      (tracklet->IsCombinatorics()) nFake += weight;
      else                                  nGood += weight;
    }
    else if (tracklet->IsGenerated() &&
	     !tracklet->IsNeutral() &&
	     TMath::Abs(tracklet->GetEta()) <= 1) nGenerated ++; // += weight;
    TIter nextBin(&toRun);
    while ((bin = static_cast<CentBin*>(nextBin()))) {
      bin->ProcessTracklet(tracklet, ip->GetZ(), signal, weight);
    }    
  }
  TIter nextBin(&toRun);
  while ((bin = static_cast<CentBin*>(nextBin()))) {
    bin->Completed();
  }    
  fNBareVsFake     ->Fill(nBare,      nFake);
  fNBareVsGood     ->Fill(nBare,      nGood);
  fNTrackletVsFake ->Fill(nMeasured,  nFake);
  fNTrackletVsGood ->Fill(nMeasured,  nGood);
  fNGeneratedVsFake->Fill(nGenerated, nFake);
  fNGeneratedVsGood->Fill(nGenerated, nGood);
}    

//____________________________________________________________________
Bool_t AliTrackletAODdNdeta::CentBin::IsAllBin() const
{
  return fLow >= fHigh;
}
//____________________________________________________________________
Bool_t AliTrackletAODdNdeta::CentBin::Accept(UInt_t   status,
					     Double_t cent,
					     Double_t ipz)
{
  Bool_t centOK = (IsAllBin() ||
		   (status & Bin2Flag(kCentrality) &&
		    cent >= fLow && cent < fHigh));
  // We get out here already as we're not in this centrality bin 
  if (!centOK) return false;
  fStatus->Fill(kAll);
  if (status & Bin2Flag(kEvent))     fStatus->Fill(kEvent);
  if (status & Bin2Flag(kTracklets)) fStatus->Fill(kTracklets);
  if (status & Bin2Flag(kTrigger))
    fStatus->Fill(kTrigger);
  else
    return false; // In case we have no trigger
  if (status & Bin2Flag(kCentrality)) {
    fStatus->Fill(kCentrality);  
    fCent->Fill(cent);
  }
  // We explicity check the IP status here, since we need that to
  // calculate the per-bin vertex efficiency
  if (status & Bin2Flag(kIP)) {
    fStatus ->Fill(kIP);
    fIPz    ->Fill(ipz);
    fCentIPz->Fill(ipz, cent);
    return true;
  }
  return false;
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
  // Get tracklet info 
  Double_t eta     = tracklet->GetEta();
  Double_t delta   = tracklet->GetDelta();
  UChar_t  flags   = tracklet->GetFlags();
  Int_t    pdgBin  = -1;
  Int_t    pdgBin2 = -1;

  // For debugging 
  char m[7];
  if (fDebug > 3) Bits2String(flags, m);

  // Check the filter mask 
  if (fMask != 0 && (flags & fMask) == 0) {
    if (fDebug > 3) Printf("%14s (0x%02x,----) %6s %7s (0x%02x) ",
			   GetName(),fMask,"reject",m,flags);
    return false;
  }

  // If we need the PDG code, get that here
  if (fEtaPdgIPz || fEtaPdg || fEtaDeltaPdg) {
    pdgBin  = PdgBin(tracklet->GetParentPdg());
    if (tracklet->IsCombinatorics())
      pdgBin2 = PdgBin(tracklet->GetParentPdg(true));
  }
  
  // If we have the eta-vs-pdg histogram, we should fill before the
  // veto, which filters out the neutral particles.  We do however,
  // check if the particle was suppressed
  if (fEtaPdg && !(tracklet->IsSuppressed()))
    fEtaPdg->Fill(eta, pdgBin, weight);

  // Check the veto mask 
  if (fVeto != 0 && (flags & fVeto) != 0) {
    if (fDebug > 3) Printf("%14s (----,0x%02x) %6s %7s (0x%02x) ",
			   GetName(), fVeto, "veto", m, flags);
    return false;
  }
  // Debug message that we accepted tracklet 
  if (fDebug > 3) Printf("%14s (0x%02x,0x%02x) %6s %7s (0x%02x) ",
			 GetName(), fMask, fVeto, "accept", m, flags);
  
  // Fill PDG,eta dependent Delta distributions
  if (fEtaDeltaPdg) {
    fEtaDeltaPdg                ->Fill(eta, delta, pdgBin,  weight);
    // Fill both particles 
    if (pdgBin2 >= 0) fEtaPdgIPz->Fill(eta, delta, pdgBin2, weight);
  }
  
  // Both reguirements (Delta < cut, dPhi < cut)
  if (fEtaIPz && (signal == 0x3))     fEtaIPz->Fill(eta, ipZ, weight);
  // Just check dPhi
  if (fEtaDeltaIPz && (signal & 0x2)) fEtaDeltaIPz->Fill(eta,delta,ipZ,weight);

  // If we do not need the eta-PDG-IPz or eta-Pt, return now 
  if (!fEtaPdgIPz && !fEtaPt) return true;

  if (fEtaPdgIPz) {
    fEtaPdgIPz                  ->Fill(eta, pdgBin,  ipZ, weight);
    // Fill both particles 
    if (pdgBin2 >= 0) fEtaPdgIPz->Fill(eta, pdgBin2, ipZ, weight);
  }

  if (fEtaPt) fEtaPt ->Fill(eta, tracklet->GetParentPt(), weight);

  return true;
}

//____________________________________________________________________
void AliTrackletAODdNdeta::CentBin::Completed()
{
  fStatus->Fill(kCompleted);
}

//____________________________________________________________________
void 
AliTrackletAODdNdeta::Terminate(Option_t*)
{
  Bool_t save = TH1::GetDefaultSumw2();
  TH1::SetDefaultSumw2();
  
  TString resName; resName.Form("%sResults",GetName());
  if (GetOutputData(2)) {
    Warning("Terminate", "Already have a result container, making a new one");
    resName.Append("New");
  }
  Container* results = new Container;
  results->SetName(resName);
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
  TH1::SetDefaultSumw2(save);
}


//____________________________________________________________________
Bool_t AliTrackletAODdNdeta::CentBin::FinalizeInit(Container* parent)
{
  fContainer     = GetC(parent, fName);
  fStatus        = GetH1(fContainer, "status");
  fCent          = GetH1(fContainer, "cent");
  fIPz           = GetH1(fContainer, "ipz");
  fCentIPz       = GetP1(fContainer, "centIpz");
  if (!fContainer || !fStatus || !fCent || !fIPz) return false;
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
  fEtaIPz      = GetH2(fContainer, "etaIPz",      false); // No complaints
  fEtaDeltaIPz = GetH3(fContainer, "etaDeltaIPz", false);
  fEtaDeltaPdg = GetH3(fContainer, "etaDeltaPdg", false);
  fEtaPdgIPz   = GetH3(fContainer, "etaPdgIPz",   false);
  fEtaPdg      = GetH2(fContainer, "etaPdg",      false);
  fEtaPt       = GetH2(fContainer, "etaPt",       false);
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
  fCentEst =
    static_cast<TProfile*>(CloneAndAdd(results,GetP(fContainer,
						    "centEstimator")));
  typedef TParameter<double> DP;
  results->Add(new DP("triggerEfficiency", fTriggerEff));

  TString tstr("UNKNOWN");
  if      (fOfflineTriggerMask & AliVEvent::kMB) {
    tstr = "MBOR ";
    if (fTriggerEff > 0 && fTriggerEff < 1) {
      tstr = "INEL";
      if (fMinEta1 > 0) tstr.Append(Form(">%d",fMinEta1-1));
    }
  }
  else if (fOfflineTriggerMask & AliVEvent::kINT7) {
    tstr = "V0AND ";
    if (fTriggerEff > 0 && fTriggerEff < 1) tstr = "NSD";
  }
  else if (fOfflineTriggerMask & AliVEvent::kINT5)
    tstr = "MBAND";
  else if (fOfflineTriggerMask == AliVEvent::kAny)
    tstr = "OFFLINE";
  results->Add(new TNamed("trigger", tstr.Data()));
  
  
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

  if (fTailMax < 0) fTailMax = fMaxDelta;
  TIter    next(fCentBins);
  CentBin* bin = 0;
  while ((bin = static_cast<CentBin*>(next()))) {
    if (!bin->MasterFinalize(results, 0, fTailDelta, fTailMax)) {
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
    TH1* h = static_cast<TH1*>(o->Clone());
    h->SetDirectory(0);
    results->Add(h);
  }
  else {
    AliWarningF("Object %p (etaWeight) is not a TH1 or not found",o);
  }
  o = fContainer->FindObject("weightCorr");
  if (o && o->IsA()->InheritsFrom(TH1::Class())) {
    TH1* h = static_cast<TH1*>(o->Clone());
    h->SetDirectory(0);
    results->Add(h);
  }
  else {
    AliWarningF("Object %p (weightCorr) is not a TH1 or not found",o);
  }
  if (!fWeights) return true;

  if (!fWeights->Retrieve(fContainer)) return false;
  fWeights->Store(results);

  return true;
}

  
//____________________________________________________________________
Bool_t AliTrackletAODdNdeta::CentBin::MasterFinalize(Container* parent,
						     TH1*       ,
						     Double_t   tailDelta,
						     Double_t   tailMax)
{
  Container* result = new Container;
  result->SetName(fName);
  result->SetOwner(true);
  parent->Add(result);

  TH1* status = static_cast<TH1*>(CloneAndAdd(result, fStatus));
  status->Scale(1./status->GetBinContent(kAll));

  Double_t   trigEff = status->GetBinContent(kTrigger);
  Double_t   vtxEff  = status->GetBinContent(kIP) / trigEff;
  typedef TParameter<double> DP;

  result->Add(new DP("triggerEfficiency", trigEff));
  result->Add(new DP("ipEfficiency",      vtxEff));
  
  Double_t   nEvents = fIPz->GetEntries();
  // Copy ipZ histogram and scale by number of events 
  TH1* ipZ = static_cast<TH1*>(CloneAndAdd(result, fIPz));
  ipZ->Scale(1./nEvents);

  TH1* cent = static_cast<TH1*>(CloneAndAdd(result, fCent));
  cent->Scale(1./nEvents);

				  
  CloneAndAdd(result, fCentIPz);

  Container* measCont = 0;
  Container* genCont  = 0;
  TIter      next(fSubs);
  Histos*    h = 0;
  while ((h = static_cast<Histos*>(next()))) {
    if (h->GetMask() != AliAODTracklet::kInjection &&
	h->GetMask() != AliAODTracklet::kCombinatorics) {
      // For the anything but injection or MC-labels, we just finalize
      if (!h->MasterFinalize(result, fIPz, tailDelta, tailMax)) {
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
    if (!EstimateBackground(result, measCont, genCont, h, tailDelta, tailMax)) {
      AliWarningF("Failed to estimate Bg in %s/%s", GetName(), h->GetName());
      return false;
    }      
  }
  return true;
}

//____________________________________________________________________
Bool_t 
AliTrackletAODdNdeta::Histos::ProjectEtaPdg(Container* result,
					    Int_t      nEvents)
{
  // Scale distribution of PDGs to number of events and bin width
  if (!fEtaPdg) return true;

  
  TH2* etaPdg = static_cast<TH2*>(fEtaPdg->Clone());
  etaPdg->SetDirectory(0);
  etaPdg->Scale(1. / nEvents, "width");
  result->Add(etaPdg);
  
  // Loop over PDG types and create 2D and 1D distributions 
  TAxis*     yaxis  = etaPdg->GetYaxis();
  Int_t      first  = yaxis->GetFirst();
  Int_t      last   = yaxis->GetLast();
  THStack*   pdgs   = new THStack("all","");
  THStack*   toPion = new THStack("toPion", "");
  THStack*   toAll  = new THStack("toAll", "");
  TH1*       pion   = 0;
  Container* pdgOut = new Container();
  pdgOut->SetName("mix");
  result->Add(pdgOut);
  pdgOut->Add(pdgs);
  pdgOut->Add(toPion);
  pdgOut->Add(toAll);

  TH1* all = static_cast<TH1*>(etaPdg->ProjectionX("total",
						   1,etaPdg->GetNbinsY()));
  all->SetDirectory(0);
  all->SetName("total");
  all->SetTitle("All");
  all->SetYTitle(Form("d#it{N}_{%s}/d#eta", "all"));
  all->SetFillColor(kBlack);
  all->SetMarkerColor(kBlack);
  all->SetLineColor(kBlack);
  all->SetMarkerStyle(20);
  all->SetFillStyle(0);
  all->SetFillColor(0);
  all->Reset();
  pdgs->Add(all);
  for (Int_t i = 1; i <= etaPdg->GetNbinsY(); i++) {
    Int_t   pdg = TString(yaxis->GetBinLabel(i)).Atoi();
    TString nme;
    Style_t sty;
    Color_t col;
    PdgAttr(pdg, nme, col, sty);
    if (pdg < 0) pdg = 0;
    if (pdg == 22) continue; // Ignore photons 

    TH1* h1 = static_cast<TH1*>(etaPdg->ProjectionX(Form("h%d", pdg),i,i));
    if (h1->GetEntries() <= 0) continue; // Do not store if empty
    h1->SetDirectory(0);
    h1->SetName(Form("eta_%d", pdg));
    h1->SetTitle(nme);
    h1->SetYTitle(Form("d#it{N}_{%s}/d#eta", nme.Data()));
    h1->SetFillColor(col);
    h1->SetMarkerColor(col);
    h1->SetLineColor(col);
    h1->SetMarkerStyle(sty);
    h1->SetFillStyle(0);
    h1->SetFillColor(0);
    h1->SetBinContent(0,0);
    all->Add(h1);
    pdgs->Add(h1);
    switch (pdg) {
    case 321:  h1->SetBinContent(0, 0.15);   break; // PRC88,044910
    case 2212: h1->SetBinContent(0, 0.05);   break; // PRC88,044910
    case 310:  h1->SetBinContent(0, 0.075);  break; // PRL111,222301
    case 3122: h1->SetBinContent(0, 0.018);  break; // PRL111,222301
    case 3212: h1->SetBinContent(0, 0.0055); break; // NPA904,539
    case 3322: h1->SetBinContent(0, 0.005);  break; // PLB734,409
    case 211:  h1->SetBinContent(0, 1);      break; // it self 
    default:   h1->SetBinContent(0, -1);     break; // Unknown
    }
    pdgOut->Add(h1);
      
    if (pdg == 211) pion = h1;
  }
  if (!pdgs->GetHists()) return true;
  
  TIter    next(pdgs->GetHists());
  TH1*     tmp  = 0;
  Double_t rmin = +1e9;
  Double_t rmax = -1e9;
  while ((tmp = static_cast<TH1*>(next()))) {
    if (tmp == all)  continue;
    // Calculate ratio to all
    TH1* rat = static_cast<TH1*>(tmp->Clone());
    rat->Divide(all);
    rat->SetDirectory(0);
    rat->SetTitle(Form("%s / all", tmp->GetTitle()));
    toAll->Add(rat);
	
    if (tmp == pion) continue;
    Double_t r276 = tmp->GetBinContent(0);
    if (r276 < 0 || r276 >= 1) continue;
	
    // Calulate ratio to pions 
    rat = static_cast<TH1*>(tmp->Clone());
    rat->Divide(pion);
    rat->SetTitle(Form("%s / %s", tmp->GetTitle(), pion->GetTitle()));
    rat->SetDirectory(0);

    TGraphErrors* g = new TGraphErrors(1);
    g->SetName(Form("%s_2760", rat->GetName()));
    g->SetTitle(Form("%s in #sqrt{s_{NN}}=2.76TeV", rat->GetTitle()));
    g->SetPoint(0,0,r276);
    g->SetPointError(0,.5,0);
    g->SetLineColor(rat->GetLineColor());
    g->SetLineStyle(rat->GetLineStyle());
    g->SetMarkerColor(rat->GetMarkerColor());
    g->SetMarkerStyle(rat->GetMarkerStyle());
    g->SetMarkerSize(1.5*rat->GetMarkerSize());
    rat->GetListOfFunctions()->Add(g,"p");
    rat->SetMaximum(TMath::Max(rat->GetMaximum(),r276));
    rat->SetMinimum(TMath::Min(rat->GetMinimum(),r276));
    rmin = TMath::Min(rat->GetMinimum(),rmin);
    rmax = TMath::Max(rat->GetMaximum(),rmax);
	
    toPion->Add(rat);
  }
  // toPion->SetMinimum(rmin);
  toPion->SetMaximum(1.1*rmax);

  return true;
}

//____________________________________________________________________
Bool_t 
AliTrackletAODdNdeta::Histos::ProjectEtaDeltaPdgPart(Container*     result,
						     Int_t          nEvents,
						     Double_t       tailDelta,
						     Double_t       tailMax,
						     const TString& pre,
						     const TString& tit)
{
  Container* pdgOut = new Container();
  pdgOut->SetName(pre);
  result->Add(pdgOut);
  TAxis*     zaxis  = fEtaDeltaPdg->GetZaxis();
  Int_t      first  = zaxis->GetFirst();
  Int_t      last   = zaxis->GetLast();
  Bool_t     isMid  = TString(pre).EqualTo("mid");
  Int_t      nX     = fEtaDeltaPdg->GetNbinsX();
  Int_t      nY     = fEtaDeltaPdg->GetNbinsY();
  Int_t      nZ     = fEtaDeltaPdg->GetNbinsZ();
  Int_t      l1     = fEtaDeltaPdg->GetXaxis()->FindBin(-1);
  Int_t      r1     = fEtaDeltaPdg->GetXaxis()->FindBin(+1);
  Int_t      b1     = (isMid ? l1 : 1);
  Int_t      b2     = (isMid ? r1 : l1-1);
  Int_t      b3     = (isMid ? -1 : r1+1);
  Int_t      b4     = (isMid ? -1 : nX);
  THStack*   stack  = new THStack("all", tit);
  pdgOut->Add(stack);

  TH1* total = fEtaDeltaPdg->ProjectionY("totalMid",1, 1,1, nZ);
  total->SetDirectory(0);
  total->SetName("total");
  total->SetTitle("Total");
  total->SetYTitle(Form("d#it{N}_{%s}/d#Delta", "total"));
  total->SetFillColor(kBlack);
  total->SetMarkerColor(kBlack);
  total->SetLineColor(kBlack);
  total->SetMarkerStyle(24);
  total->SetFillStyle(0);
  total->SetFillColor(0);
  total->Reset();
  stack->Add(total);

  THStack* ratios = new THStack("ratios","");
  TH1* ratioSig = Make1D(0, "ratioSig", "Ratio to all", kGreen+1, 24, *zaxis);
  TH1* ratioBg  = Make1D(0, "ratioBg",  "Ratio to all", kRed+1,   25, *zaxis);
  ratioSig->GetXaxis()->LabelsOption("v");
  ratioBg ->GetXaxis()->LabelsOption("v");
  ratioSig->SetFillColor(kGreen+1);
  ratioBg ->SetFillColor(kRed  +1);
  ratioSig->SetFillStyle(1001);
  ratioBg ->SetFillStyle(1001);
  ratioSig->SetBarWidth(0.4);
  ratioBg ->SetBarWidth(0.4);
  ratioSig->SetBarOffset(0.05);
  ratioBg ->SetBarOffset(0.55);
  ratios->Add(ratioSig, "bar0 e text30");
  ratios->Add(ratioBg,  "bar0 e text30");
  pdgOut->Add(ratios);

  THStack* rows = new THStack("rows","");
  TH1* rowSig = new TH1D("rowSig","row",6,0,6);
  rowSig->SetDirectory(0);
  rowSig->SetYTitle("Fraction of signal");
  rowSig->GetXaxis()->SetBinLabel(1, "K^{0}_{S}");
  rowSig->GetXaxis()->SetBinLabel(2, "K^{#pm}");
  rowSig->GetXaxis()->SetBinLabel(3, "#Lambda");
  rowSig->GetXaxis()->SetBinLabel(4, "#Xi");
  rowSig->GetXaxis()->SetBinLabel(5, "#Sigma");
  rowSig->GetXaxis()->SetBinLabel(6, "Other");
  rowSig->GetXaxis()->LabelsOption("v");
  rowSig->SetMaximum(100);
  rowSig->SetMinimum(0);
  rowSig->SetLineColor(kGreen+1);
  rowSig->SetMarkerColor(kGreen+1);
  rowSig->SetMarkerStyle(24);
  TH1* rowBg = static_cast<TH1*>(rowSig->Clone("rowBg"));
  rowBg->SetDirectory(0);
  rowBg->SetYTitle("Fraction of background");
  rowBg->SetLineColor(kRed+1);
  rowBg->SetMarkerColor(kRed+1);
  rowBg->SetMarkerStyle(25);
  rowSig->SetFillStyle(1001);
  rowBg ->SetFillColor(1001);
  rowSig->SetFillColor(kGreen+1);
  rowBg ->SetFillColor(kRed  +1);
  rowSig->SetBarWidth(0.4);
  rowBg ->SetBarWidth(0.4);
  rowSig->SetBarOffset(0.05);
  rowBg ->SetBarOffset(0.55);
  rows->Add(rowSig, "bar0 e text30");
  rows->Add(rowBg,  "bar0 e text30");
  pdgOut->Add(rows);

  TH1* hK0S    = 0;
  TH1* hKpm    = 0;
  TH1* hLambda = 0;
  TH1* hXi     = 0;
  TH1* hSigma  = 0;
  TH1* hOther  = 0;
  for (Int_t i = 1; i <= zaxis->GetNbins(); i++) {
    Int_t   pdg = TString(zaxis->GetBinLabel(i)).Atoi();
    TString nme;
    Style_t sty;
    Color_t col;
    Int_t   ipdg = pdg;    
    PdgAttr(pdg, nme, col, sty);
    if (pdg < 0) pdg = 0;
    ratioSig->GetXaxis()->SetBinLabel(i, nme);
    ratioBg ->GetXaxis()->SetBinLabel(i, nme);
    if (pdg == 22) continue; // Ignore photons 

    
    TH1* h1 = fEtaDeltaPdg->ProjectionY(Form("%d", pdg), b1, b2, i,i);
    h1->SetDirectory(0);
    if (b3 < b4) {
      TH1* h2 = fEtaDeltaPdg->ProjectionY(Form("t%d", pdg), b3, b4, i,i);
      h2->SetDirectory(0);
      h1->Add(h2);
      delete h2;
    }
    h1->SetUniqueID(i); // Store bin number 
    if (h1->GetEntries() <= 0) continue; // Do not store if empty
    h1->SetName(Form("%d", pdg));
    h1->SetTitle(nme);
    h1->SetYTitle(Form("d#it{N}_{%s}/d#Delta", nme.Data()));
    h1->SetFillColor(col);
    h1->SetMarkerColor(col);
    h1->SetLineColor(col);
    h1->SetMarkerStyle(sty);
    h1->SetFillStyle(0);
    h1->SetFillColor(0);
    h1->SetBinContent(0,ipdg);
    total->Add(h1);
    stack->Add(h1);

    TH1** ptr = 0;
    switch (pdg) {
    case 321:  ptr = &hKpm;   break;  // K^{+}
    case 310:  ptr = &hK0S;   break;  // K^{0}_{S}
      // case 3212:                   // #Sigma^{0} // Old, wrong code
    case 3112:                        // #Sigma^{-}
    case 3222: ptr = &hSigma; break;  // #Sigma^{+} // New, right code 	 
      
      // case 3322:                    // #Xi^{0}   // Old, wrong code 
    case 3312: ptr = &hXi;    break;   // #Xi^{-}   // New, right code 
    case 3122: ptr = &hLambda;break;   // #Lambda  
    default:   ptr = &hOther; break;
    }
    if (!*ptr) {
      *ptr = static_cast<TH1*>(h1->Clone());
      (*ptr)->Reset();
      (*ptr)->SetDirectory(0);
      if (ptr == &hOther) {
	(*ptr)->SetTitle("Other");
	(*ptr)->SetName("0");
      }
    }
    // Printf("Adding %s to %s", h1->GetName(), (*ptr)->GetName());
    (*ptr)->Add(h1);
  }
  total->SetBinContent(0,0);
  total->Scale(1./nEvents);
  Double_t deltaCut   = 1.5;
  Double_t eTot, iTot = Integrate(total, 0,         tailMax,  eTot);
  Double_t eSig, iSig = Integrate(total, 0,         deltaCut, eSig);
  Double_t eBg,  iBg  = Integrate(total, tailDelta, tailMax,  eBg);

  TH1* totInt = new TH1D("totalIntegrals", "Integrals of total:", 3, 0, 3);
  totInt->SetDirectory(0);
  totInt->GetXaxis()->SetBinLabel(1,Form("0#minus%4.1f", tailMax));
  totInt->GetXaxis()->SetBinLabel(2,Form("0#minus%3.1f",deltaCut));
  totInt->GetXaxis()->SetBinLabel(3,Form("%3.1f#minus%4.1f",tailDelta,tailMax));
  totInt->SetBinContent(1,iTot); totInt->SetBinError(1, eTot);
  totInt->SetBinContent(2,iSig); totInt->SetBinError(2, eSig);
  totInt->SetBinContent(3,iBg);  totInt->SetBinError(3, eBg);
  pdgOut->Add(totInt);
  
  if (!stack->GetHists()) {
    AliWarningF("No histograms in the stack %s", stack->GetName());    
    return true;
  }
  TH1*  spec[] = { hK0S, hKpm, hLambda, hXi, hSigma, hOther, 0 };
  TH1** ptr    = spec;
  Int_t pbin   = 1;
  for (Int_t pbin=1; pbin<=6; pbin++) {
    TH1*     h            = spec[pbin-1];
    if (!h) continue;
    h->Scale(1./  nEvents);
    Double_t eiSig, iiSig = Integrate(h, 0,         deltaCut, eiSig);
    Double_t eiBg,  iiBg  = Integrate(h, tailDelta, tailMax,  eiBg);
    Double_t erS,  rS     = RatioE(iiSig, eiSig, iSig, eSig,  erS);
    Double_t erB,  rB     = RatioE(iiBg,  eiBg,  iBg,  eBg,   erB);
#if 0
    Printf("%10s (%d) signal=%6.4f+/-%6.4f  background=%6.4f+/-%6.4f  "
	   "ratios %6.4f+/-%6.4f  %6.4f+/-%6.4f",
	   h->GetTitle(), pbin, iiSig, eiSig, iiBg, eiBg, rS, erS, rB, erB);
#endif 
    rowSig->SetBinContent(pbin,100*rS);
    rowSig->SetBinError  (pbin,100*TMath::Sqrt(erS));
    rowBg ->SetBinContent(pbin,100*rB);
    rowBg ->SetBinError  (pbin,100*TMath::Sqrt(erB));
    pdgOut->Add(h);
  }
  TIter    next(stack->GetHists());
  TH1*     tmp  = 0;
  while ((tmp = static_cast<TH1*>(next()))) {
    // Calculate ratio to all
    // Scale(tmp, iTot, eTot);
    if (tmp == total) continue;
    Int_t    pdg          = tmp->GetBinContent(0); tmp->SetBinContent(0,0);
    tmp->Scale(1./  nEvents);
    Double_t eiSig, iiSig = Integrate(tmp, 0,         deltaCut, eiSig);
    Double_t eiBg,  iiBg  = Integrate(tmp, tailDelta, tailMax,  eiBg);
    Double_t erS,  rS     = RatioE(iiSig, eiSig, iSig, eSig,  erS);
    Double_t erB,  rB     = RatioE(iiBg,  eiBg,  iBg,  eBg,   erB);
    Int_t    bin          = tmp->GetUniqueID();
#if 0
    Printf(" %10s(%4d,%3d) "
	   "signal=%6.1f/%6.1f=%5.3f+/-%5.3f bg=%6.1f/%6.1f=%5.3f+/-%5.3f",
	   tmp->GetTitle(), pdg, bin, iiSig, iSig, rS, erS, iiBg, iBg, rB, erB); 
#endif 
    ratioSig->SetBinContent(bin, rS);
    ratioSig->SetBinError  (bin, erS);
    ratioBg ->SetBinContent(bin, rB);
    ratioBg ->SetBinError  (bin, erB);
  } // while(tmp)

  return true;
}

//____________________________________________________________________
Bool_t 
AliTrackletAODdNdeta::Histos::ProjectEtaDeltaPdg(Container* result,
						 Int_t      nEvents,
						 Double_t   tailDelta,
						 Double_t   tailMax)
{
  // Scale distribution of PDGs to number of events and bin width
  if (!fEtaDeltaPdg) return true;

  Container* pdgOut = new Container();
  pdgOut->SetName("specie");
  result->Add(pdgOut);
  
  ProjectEtaDeltaPdgPart(pdgOut,
			 nEvents,
			 tailDelta,
			 tailMax,
			 "mid",
			 "|#eta|<1");
  ProjectEtaDeltaPdgPart(pdgOut,
			 nEvents,
			 tailDelta,
			 tailMax,
			 "fwd",
			 "|#eta|>1");
}

//____________________________________________________________________
Bool_t 
AliTrackletAODdNdeta::Histos::ProjectEtaPdgIPz(Container*     result,
					       TH1*           ipz,
					       const TString& shn)
{
  // If we have the PDG distributions, we project on eta,IPz, and then
  // average over IPz to get dNpdg/deta.
  if (!fEtaPdgIPz) return true;
  
  // Scale each vertex range by number of events in that range
  TH3* etaPdgIPz = ScaleToIPz(fEtaPdgIPz, ipz, false);
  result->Add(etaPdgIPz);
  
  // Loop over PDG types and create 2D and 1D distributions 
  TAxis*     yaxis  = etaPdgIPz->GetYaxis();
  Int_t      first  = yaxis->GetFirst();
  Int_t      last   = yaxis->GetLast();
  THStack*   pdgs   = new THStack("all","");
  THStack*   ratios = new THStack("toPion", "");
  TH1*       pion   = 0;
  TH2*       dtfs   = 0;
  Container* pdgOut = new Container();
  pdgOut->SetName("types");
  result->Add(pdgOut);
  pdgOut->Add(pdgs);
  pdgOut->Add(ratios);
  for (Int_t i = 1; i <= etaPdgIPz->GetNbinsY(); i++) {
    yaxis->SetRange(i,i);
    
    Int_t   pdg = TString(yaxis->GetBinLabel(i)).Atoi();
    TString nme;
    Style_t sty;
    Color_t col;
    PdgAttr(pdg, nme, col, sty);
    if (pdg < 0) pdg = 0;
    
    TH2*   h2    = static_cast<TH2*>(etaPdgIPz->Project3D("zx e"));
    if (h2->GetEntries() <= 0) continue; // Do not store if empty
    h2->SetDirectory(0);
    h2->SetName(Form("etaIPz_%d", pdg));
    h2->SetTitle(Form("%s#rightarrowX_{%s}", nme.Data(), shn.Data()));
    h2->SetYTitle(Form("d#it{N}^{2}_{%s#rightarrow%s}/"
		       "(d#etadIP_{#it{z}})",
		       nme.Data(), shn.Data()));
    h2->SetFillColor(col);
    h2->SetMarkerColor(col);
    h2->SetLineColor(col);
    pdgOut->Add(h2);
    
    TH1* h1 = AverageOverIPz(h2, Form("eta_%d",pdg), 1, ipz, 0, 0, false);
    h1->SetDirectory(0);
    h1->SetYTitle(Form("d#it{N}_{%s#rightarrow%s}/d#eta",
		       nme.Data(), shn.Data()));
    h1->SetMarkerStyle(sty);
    pdgs->Add(h1);
    
    if (pdg == 211) pion = h1;
    TH2* tmp = 0;
    switch (pdg) {
    case 321: 	  // Strange meson K^{+} 	     break;
    case 323: 	  // Strange meson K^{*+} 	     break;
    case 310: 	  // Strange meson K^{0}_{S} 	     break;
    case 130: 	  // Strange meson K^{0}_{L} 	     break;
    case 311: 	  // Strange meson K^{0} 	     break;
    case 313: 	  // Strange meson K^{*} 	     break;
    case 221: 	  // Strange meson #eta 	     break;
    case 333: 	  // Strange meson #varphi 	     break;
    case 331: 	  // Strange meson #eta' 	     break;
    case 3112: 	  // Strange baryon #Sigma^{-}       break;
    case 3222: 	  // Strange baryon #Sigma^{+}       break;
    case 3114: 	  // Strange baryon #Sigma^{*-}      break;
    case 3224: 	  // Strange baryon #Sigma^{*+}      break;
    case 3312: 	  // Strange baryon #Xi^{-} 	     break;
    case 3314: 	  // Strange baryon #Xi^{*-} 	     break;
    case 3122: 	  // Strange baryon #Lambda 	     break;
    case 3212: 	  // Strange baryon #Sigma^{0}       break;
    case 3214: 	  // Strange baryon #Sigma^{*0}      break;
    case 3322: 	  // Strange baryon #Xi^{0} 	     break;
    case 3324: 	  // Strange baryon #Xi^{*0} 	     break;
      tmp = h2;
      break;
    default: break;
    }
    if (tmp) {
      if (!dtfs) {
	dtfs = static_cast<TH2*>(tmp->Clone("dtfs"));
	dtfs->SetDirectory(0);
	dtfs->SetTitle("X_{s}#rightarrowX");
	dtfs->SetYTitle("IP_{#it{z}} [cm]");
	dtfs->Reset();
	result->Add(dtfs);
      }
      // Add up contributions from strange particles 
      dtfs->Add(tmp); 
    }	
  }
  if (!pdgs->GetHists()) {
    yaxis->SetRange(first, last);    
    return true;
  }
  
  TIter    next(pdgs->GetHists());
  TH1*     tmp = 0;
  while ((tmp = static_cast<TH1*>(next()))) {
    if (tmp == pion) continue;
    TH1* rat = static_cast<TH1*>(tmp->Clone());
    rat->Divide(pion);
    rat->SetDirectory(0);
      ratios->Add(rat);
  }
  
  yaxis->SetRange(first, last);    
}

//____________________________________________________________________
Bool_t 
AliTrackletAODdNdeta::Histos::MasterFinalize(Container* parent,
					     TH1*       ipz,
					     Double_t   tailDelta,
					     Double_t   tailMax)
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
  // Scale distribution of Pt to number of events and bin width
  if (fEtaPt) {
    TH2* etaPt = static_cast<TH2*>(fEtaPt->Clone());
    etaPt->SetDirectory(0);
    etaPt->Scale(1. / nEvents, "width");
    result->Add(etaPt);
  }
  // Short-hand-name
  TString shn(etaIPz ? etaIPz->GetTitle() : "X"); 
  
  ProjectEtaPdg     (result, nEvents);
  ProjectEtaDeltaPdg(result, nEvents, tailDelta, tailMax);
  ProjectEtaPdgIPz  (result, ipz, shn);
  
  // If we do not have eta vs Delta, just return 
  if (!fEtaDeltaIPz) return true;

  // Normalize delta distribution to integral number of events
  // static_cast<TH3*>(CloneAndAdd(result, fEtaDeltaIPz));
  TH3* etaDeltaIPz = ScaleToIPz(fEtaDeltaIPz, ipz, false); 
  // ipz->GetEntries()>1000);
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
  Int_t    sigBin   = etaDeltaIPz->GetYaxis()->FindBin(1.5);
  Int_t    highBin  = TMath::Min(etaDeltaIPz->GetYaxis()->FindBin(tailMax),
				 etaDeltaIPz->GetYaxis()->GetNbins());  

  
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
  TH2* etaIPzDeltaMain = static_cast<TH2*>(etaDeltaIPz->Project3D("zx e"));
  etaIPzDeltaMain->SetName("etaIPzDeltaMain");
  etaIPzDeltaMain->SetDirectory(0);
  etaIPzDeltaMain->Reset();
  etaIPzDeltaMain->SetTitle(etaDeltaTail->GetTitle());
  etaIPzDeltaMain->SetZTitle(etaDelta->GetYaxis()->GetTitle());
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

      intg = etaDeltaIPz->IntegralAndError(i,i,1,sigBin,j,j,eintg);
      etaIPzDeltaMain->SetBinContent(i, j, intg);
      etaIPzDeltaMain->SetBinError  (i, j, eintg);
    }
  }
  result->Add(etaIPzDeltaTail);
  result->Add(etaIPzDeltaMain);
  result->Add(etaDeltaTail);

  // Integrate full tail
  intg = etaDeltaIPz->IntegralAndError(1,etaDeltaIPz->GetNbinsX(),
				       lowBin, highBin,
				       1,etaDeltaIPz->GetNbinsZ(),
				       eintg);
  result->Add(new TParameter<double>("deltaTail",      intg));
  result->Add(new TParameter<double>("deltaTailError", eintg));

  // Some consistency checks:
  if (fDebug > 1) {
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
  }
				    
  return true;
}

//____________________________________________________________________
Bool_t AliTrackletAODdNdeta::CentBin::EstimateBackground(Container* result,
							 Container* measCont,
							 Container* genCont,
							 Histos*    h,
							 Double_t   tailCut,
							 Double_t   tailMax)
{
  if (!h || !measCont) {
    AliWarningF("No sub-histos or measured container in %s", GetName());
    return false;
  }

  if (!h->MasterFinalize(result, fIPz, tailCut,tailMax)) { 
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
      if (!wobj->IsA()->InheritsFrom(AliTrackletBaseWeights::Class())) {
	::Warning("Create", "Object %s from file %s not an "
		  "AliTrackletBaseWeights but a %s",
		  wnam.Data(), wfile->GetName(), wobj->ClassName());
	return 0;
      }
      wret->SetWeights(static_cast<AliTrackletBaseWeights*>(wobj));
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
