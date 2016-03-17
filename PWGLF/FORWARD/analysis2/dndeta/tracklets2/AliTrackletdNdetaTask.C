#ifndef ALITRACKLETDNDETATASK_H
#define ALITRACKLETDNDETATASK_H
#include <AliAnalysisTaskSE.h>
#include <AliTrackletdNdetaUtils.C>
#ifndef __CINT__
#include <TGeoGlobalMagField.h>
#include <TProfile.h>
#include <AliVEvent.h>
#include <AliVVertex.h>
#include <AliVertex.h>
#include <AliVMultiplicity.h>
#include <AliMultiplicity.h>
#include <AliAnalysisManager.h>
#include <AliVEventHandler.h>
#include <AliESDInputHandlerRP.h>
#include <AliMultSelection.h>
#include <AliITSMultRecBg.h>
#include <AliCDBManager.h>
#include <AliCDBEntry.h>
#include <AliCDBPath.h>
#include <AliCDBId.h>
#include <AliGeomManager.h>
#include <TGeoManager.h>
#else
class AliVEvent;
class AliVVertex;
class AliMultiplicity;
class AliAnalysisDataContainer;
class AliITSMultRecBg;          
class AliMultSelection;         // Auto-load
class TGeoGlobalMagField;       // Auto-load
class TGeoManager;              // Auto-load 
class AliGeomManager;           // Auto-load
class AliCDBManager;            // Auto-load 
#endif

namespace {
  Int_t lDebug = 0;
}
/**
 * A task to analyse SPD tracklet data for the charged particle
 * pseudorapidity density - real data case 
 * 
 */
class AliTrackletdNdetaTask : public AliAnalysisTaskSE,
			      public AliTrackletdNdetaUtils
{
public:
  //==================================================================
  /** 
   * Flags for what to reconstruct 
   */
  enum {
    kNormal    = 0x1,
    kInjection = 0x2,
    kRotation  = 0x4
  };
  enum {
    kAll=1,          // Count all events
    kEvent,          // Have event
    kMultiplicity,   // Have multiplicity
    kClusters,       // Have clusters
    kTrigger,        // Have trigger
    kIP,             // Have IP
    kCentrality,     // Have Centrality
    kReconstructed,  // Have reconstructed
    kInjected,       // Have injected
    kRotated,        // Have rotated
    kCompleted       // Have completed
  };
  //==================================================================
  /** 
   * Base class for sub components 
   */
  struct SubBase : public TObject
  {
    typedef AliTrackletdNdetaTask::Container Container;
    /** Container of histograms for this sub-component */
    Container* fContainer;
    /** 
     * Constructor
     */    
    SubBase() : fContainer(0) {}
    /**
     * Copy constructor 
     */
    SubBase(const SubBase& o) : TObject(o), fContainer(0) {}
    /**
     * Destructor 
     */
    virtual ~SubBase() {}
    /** 
     * Assignment operator 
     * 
     * @return Reference to this object 
     */
    SubBase& operator=(SubBase& o) { return *this; }
    /** 
     * Get the name - must be overloaded
     * 
     * @return Name of sub 
     */
    virtual const char* Name() const = 0;
    /** 
     * Create our container and add to parent 
     * 
     * @param parent  Parent container 
     */
    virtual Container* CreateContainer(Container& parent);
    /** 
     * Initialize this sub-component 
     * 
     * @param parent      Parent container of histograms 
     * @param etaAxis     Psuedorapidity axis 
     * @param ipzAxis     IP's z-coordinate axis 
     * @param deltaAxis   Axis of @f$ \Delta@f$ 
     * @param dThetaAxis  Axis of @f$ d\theta@f$ 
     * @param dPhiAxis    Axis of @f$ d\phi@f$ 
     */
    virtual void WorkerInit(Container&     parent,
			    const TAxis&   etaAxis,
			    const TAxis&   ipzAxis,
			    const TAxis&   deltaAxis,
			    const TAxis&   dThetaAxis,
			    const TAxis&   dPhiAxis) = 0;
    /** 
     * Recreate from container 
     * 
     * @param parent Parent container 
     */
    virtual void FinalizeInit(Container* parent) = 0;
    /** 
     * Called on master on merged result. 
     * 
     * @param results   Output container for results 
     * @param ipz       Histogram of ipZ for normalization 
     * @param deltaCut  Upper cut on @f$\Delta@f$ signal
     * @param deltaTail Lower cut on @f$\Delta@f$ tail 
     *
     * @return Result container 
     */
    virtual Container* MasterFinalize(Container* parent,
				      TH1*       ipz,
				      Double_t   deltaCut, 
				      Double_t   deltaTail) = 0;

    ClassDef(SubBase,1); // Base class for sub structures
  };
  //==================================================================
  struct HistoSet : public SubBase
  {
    TString fName;
    TH2*    fEtaVsDelta;    // Eta vs delta
    TH2*    fDPhiVsdThetaX; // Normalized dphi vs dtheta*x
    TH2*    fDPhiVsdTheta;  // dphi vs dtheta
    TH1*    fDelta;         // Distance
    TH2*    fEtaVsIPz;      // Eta vs IPz
    Color_t fColor;         // Color for this histogram set
    Style_t fStyle; 
    /**
     * Constructor 
     */
    HistoSet(const char* name="",
	     Color_t color=kBlack,
	     Style_t style=20)
      : SubBase(),
	fName(name),
	fEtaVsDelta(0),
	fDPhiVsdTheta(0),
	fDPhiVsdThetaX(0),
	fDelta(0),
	fEtaVsIPz(0),
	fColor(color),
	fStyle(style)
    {}
    HistoSet(const HistoSet& o) 
      : SubBase(o),
	fName(o.fName),
	fEtaVsDelta(0),
	fDPhiVsdTheta(0),
	fDPhiVsdThetaX(0),
	fDelta(0),
	fEtaVsIPz(0),
	fColor(o.fColor),
	fStyle(o.fStyle)
    {}
    /** 
     * Assignment operator 
     * 
     * @return Reference to this object 
     */
    HistoSet& operator=(const HistoSet&) { return *this; }    
    /**
     * Destructor
     */
    virtual ~HistoSet() {}
    /** 
     * Get the name 
     * 
     * @return The name 
     */
    virtual const char* Name() const { return fName.Data(); }
    /** 
     * Initialize this object 
     * 
     * @param parent      Parent container of histograms 
     * @param etaAxis     Psuedorapidity axis 
     * @param ipzAxis     IP's z-coordinate axis 
     * @param deltaAxis   Axis of @f$ \Delta@f$ 
     * @param dThetaAxis  Axis of @f$ d\theta@f$ 
     * @param dPhiAxis    Axis of @f$ d\phi@f$ 
     */
    virtual void WorkerInit(Container&     parent,
			    const TAxis&   etaAxis,
			    const TAxis&   ipzAxis,
			    const TAxis&   deltaAxis,
			    const TAxis&   dThetaAxis,
			    const TAxis&   dPhiAxis);
    /** 
     * Recreate from container 
     * 
     * @param parent Parent container 
     */
    virtual void FinalizeInit(Container* parent);
    /** 
     * Called on master on merged result. 
     * 
     * @param parent    Output container for results 
     * @param ipz       Histogram of ipZ for normalization 
     * @param deltaCut  Upper cut on @f$\Delta@f$ signal
     * @param deltaTail Lower cut on @f$\Delta@f$ tail 
     *
     * @return Result container 
     */
    virtual Container* MasterFinalize(Container* parent,
				      TH1*       ipz,
				      Double_t   deltaCut,
				      Double_t   deltaTail);
    /** 
     * Fill into histograms 
     *  
     * @param ipZ      IP's z-coordiante 
     * @param eta      Pseudorapidity 
     * @param dphi     @f$ d\phi@f$ 
     * @param dphis    Shifted @f$ d\phi@f$ 
     * @param dtheta   @f$ d\theta@f$ 
     * @param dthetax  Scaled @f$ d\theta@f$ 
     * @param delta    Distance @f$ \Delta@f$ 
     * @param weight   Particle weight (for MC)
     */
    virtual void Fill(Double_t ipZ,
		      Double_t eta,
		      Double_t dphi,
		      Double_t dphis,
		      Double_t dtheta,
		      Double_t dthetax,
		      Double_t delta,
		      Double_t weight=1);
    ClassDef(HistoSet,1); // Class to hold histograms 
  };    
  //==================================================================
  /** 
   * A centrality bin  
   */
  struct CentBin : public SubBase
  {
    /** 
     * Constructor 
     * 
     * @param cmin     Least centrality (inclusive)
     * @param cmax     Largest centrality (exclusive)
     * @param recFlags Reconstruction flags 
     * @param colOff   Color offset for histograms 
     */
    CentBin(Float_t  cmin,
	    Float_t  cmax,
	    UShort_t recFlags,
	    Int_t    colOff=1,
	    Int_t    styOff=0);
    /** 
     * Destructor 
     */
    virtual ~CentBin() {}
    /** 
     * Get the name of the bin.  Return value must be copied by user,
     * as it is formated in a circular buffer.
     * 
     * @return Name of bin 
     */
    virtual const char* Name() const;
    /** 
     * Initialize this object on worker - e.g., called at start of
     * (slave) job
     * 
     * @param parent  Parent container 
     * @param etaAxis Eta axis used 
     * @param ipzAxis IPz axis used 
     */
    virtual void WorkerInit(Container&       parent,
			    const TAxis&     etaAxis,
			    const TAxis&     ipzAxis,
			    const TAxis&     deltaAxis,
			    const TAxis&     dThetaAxis,
			    const TAxis&     dPhiAxis);
    /** 
     * Recreate from container 
     * 
     * @param parent Parent container 
     */
    virtual void FinalizeInit(Container* parent);
    /** 
     * Calculate background delta distribution 
     * 
     * @param dataCon   Data container 
     * @param set       Background histogram set 
     * @param result    Where to store the results 
     * @param deltaCut  Upper cut on @f$\Delta@f$ signal
     * @param deltaTail Lower cut on @f$\Delta@f$ tail 
     */
    virtual Container* EstimateBackground(Container* dataCon,
					  HistoSet*  set,
					  Container* result,
					  Double_t   deltaCut,
					  Double_t   deltaTail);
    /** 
     * Called on master on merged result. 
     * 
     * @param parent    Output container for results 
     * @param ipz       Histogram of ipZ for normalization 
     * @param deltaCut  Upper cut on @f$\Delta@f$ signal
     * @param deltaTail Lower cut on @f$\Delta@f$ tail 
     *
     * @return Result container 
     */
    virtual Container* MasterFinalize(Container* parent,
				      TH1*       ipz,
				      Double_t   deltaCut,
				      Double_t   deltaTail);
    /** 
     * @{ 
     * @name 
     */
    /** 
     * Check if the value belong in this bin 
     * 
     * @param value Centrality 
     * 
     * @return true if value is within this bin 
     */
    virtual Bool_t Accept(Float_t value) const
    {
      return (fMin <= value && value < fMax);
    }
    /** 
     * Fill IP z coordinate 
     * 
     * @param ipZ IP z coordinate 
     */
    virtual void FillIPz(Double_t ipZ) { fIPz->Fill(ipZ); }
    /** 
     * Fill information into histogram set. This is delegated so that
     * sub-classes my override to do more (e.g., for MC).
     * 
     * @param set             Histogram set to fill into 
     * @param isSignal        True if signal survived cuts 
     * @param ipZ             Current event IP's z coordinate 
     * @param eta             Pseudorapidity 
     * @param dPhi            Opening angle in azimuth
     * @param dPhiS           Shifted opening angle in azimuth
     * @param dTheta          Opening polar angle
     * @param dThetaX         Scaled opening polar angle 
     * @param delta           Tracklet @f$\chi^2@f$ 
     * @param weight          Tracklet weight (for MC) 
     * 
     * @return true on success 
     */
    virtual Bool_t FillSet(UShort_t                set,
			   Bool_t                  isSignal,
			   Double_t                ipZ,
			   Double_t                eta,
			   Double_t                dPhi,
			   Double_t                dPhiS,
			   Double_t                dTheta,
			   Double_t                dThetaX,
			   Double_t                delta,
			   Double_t                weight);
    /* @} */
  protected:
    /** 
     * @{ 
     * @name Construction, etc. 
     */
    /** 
     * Constructor 
     * 
     * @param cmin Least centrality (inclusive)
     * @param cmax Largest centrality (exclusive)
     */
    CentBin()
      : SubBase(), 
	fMin(0),
	fMax(0),
	fHistoSets(0),
	fDataSet(0),
	fInjSet(0),
	fRotSet(0),
	fIPz(0)
    {}					  
    /** 
     * Copy Constructor
     */
    CentBin(const CentBin& o)
      : SubBase(o),
	fMin(o.fMin),
	fMax(o.fMax),
	fHistoSets(0),
	fDataSet(0),
	fInjSet(0),
	fRotSet(0),
	fIPz(0)
    {
    }
    /** 
     * Assignment operator 
     *
     * @return reference to this object 
     */
    CentBin& operator=(const CentBin&) { return *this; }
    /* @} */
			 
    /** 
     * @{ 
     * @name 
     */
    Float_t    fMin;          // Least value 
    Float_t    fMax;          // Largest value
    TList*     fHistoSets;    // List of histogram sets
    HistoSet*  fDataSet;      // Normal reconstruction histograms 
    HistoSet*  fInjSet;       // Injected histograms 
    HistoSet*  fRotSet;       // Rotation histograms
    TH1*       fIPz;          // Vertex distribution
    /* @} */
    ClassDef(CentBin,1); // Class to centrality dependent stuff
  };
  //==================================================================
  /** 
   * @{ 
   * @name Construction, copy, assign
   */
  /** 
   * Constructor
   *
   * @param name Name of task 
   */
  AliTrackletdNdetaTask(const char* name);
  /** 
   * Default constructor - for ROOT I/O only 
   */
  AliTrackletdNdetaTask()
    : AliAnalysisTaskSE(),
      AliTrackletdNdetaUtils(),
      fContainer(0),
      fCentBins(0),
      fCentMethod("V0M"),
      fCentAxis(1,0,100),
      fIPzAxis(20,-10,10),
      fEtaAxis(80,-2,2),
      fPhiAxis(200,0,TMath::TwoPi()),
      fCent(0),
      fIPz(0),
      fEtaPhi(0),
      fUsedClusters0(0),
      fUsedClusters1(0),
      fAllClusters0(0),
      fAllClusters1(0),
      fStatus(0),
      fRecMode(kNormal|kInjection),
      fScaleDTheta(false),
      fDPhiShift(0.0045),
      fShiftedDPhiCut(-1),
      fScaledDThetaCut(-1),
      fDeltaCut(1.5),
      fMaxDelta(25),
      fTailDelta(5),
      fDThetaWindow(0.025),
      fDPhiWindow(0.06),
      fPhiOverlapCut(0.005),
      fZEtaOverlapCut(0.05),
      fPhiRotation(TMath::Pi())
  {} 
  /** 
   * Destructor
   */
  virtual ~AliTrackletdNdetaTask() {}
  //__________________________________________________________________
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
  void FinishTaskOutput() { WorkerFinalize(); }
  /** 
   * Called at end of master job on merged results. 
   * 
   */
  void Terminate(Option_t*);
  /** 
   * Connect this task to the train
   * 
   * @param sumFile (optional) event sum file 
   * @param resFile (optional) result file 
   * 
   * @return true on success 
   */
  virtual Bool_t Connect(const char* sumFile=0,
			 const char* resFile=0);
  /** 
   * @{ 
   * @name Setters 
   */
  void SetReconstructionMode(const TString& mode)
  {
    UShort_t flags = 0;
    TString m(mode); m.ToLower();
    if (m.Contains("nor"))  flags |= kNormal;
    if (m.Contains("inj"))  flags |= kInjection;
    if (m.Contains("rot"))  flags |= kRotation;
    SetReconstructionMode(flags);
  }
  /** 
   * Set the reconstruction mode to use 
   * 
   * @param flags Bit mask of flags 
   */
  void SetReconstructionMode(UShort_t flags) { fRecMode = flags; }
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
  void SetEtaAxis(const TString& spec)
  {
    SetAxis(fEtaAxis, spec);
  }
  /**
   * Set wether to scale @f$\Delta\theta@f$ by @f$\sin^2\theta@f$ 
   *
   * @param x If true, scale 
   */
  void SetScaleDTheta(Bool_t x=false) { fScaleDTheta = x; }
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
   * Set (optional) cut on @f$\Delta\theta/\sin^2\theta@f$ 
   * 
   * @param x Cut on scaled @f$\Delta\theta@f$ 
   */
  void SetScaleDThetaCut(Double_t x=-1) { fScaledDThetaCut = x; }
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
   * Set DThetaWindow
   *
   * @param x Value 
   */
  void SetDThetaWindow(Double_t x=0.025) { fDThetaWindow = x; }
  /**
   * Set @f$ d\phi@f$ window - used for reconstruction only
   *
   * @param x Value 
   */
  void SetDPhiWindow(Double_t x=0.06) { fDPhiWindow = x; }
  /**
   * Set PhiOverlapCut - used for reconstruction only
   *
   * @param x Value
   */
  void SetPhiOverlapCut(Double_t x=0.005) { fPhiOverlapCut = x; }
  /**
   * Set ZEtaOverlapCut - used for reconstruction only
   *
   * @param x Value
   */
  void SetZEtaOverlapCut(Double_t x=0.05) { fZEtaOverlapCut = x; }
  /**
   * Set PhiRotation - used for rotated reconstruction only
   *
   * @param x Value
   */
  void SetPhiRotation(Double_t x=TMath::Pi()) { fPhiRotation = x; }
    /* @} */
  /** 
   * @{
   * @name Other functions 
   */
  virtual void Print(Option_t*) const;
protected:
  //__________________________________________________________________
  /** 
   * @{ 
   * @name Construction, copy, assign
   */
  /** 
   * Copy constructor 
   * 
   * @param o Object to copy from 
   */
  AliTrackletdNdetaTask(const AliTrackletdNdetaTask& o)
    : AliAnalysisTaskSE(o),
      AliTrackletdNdetaUtils(o),
      fContainer(0),
      fCentBins(0),
      fCentMethod("V0M"),
      fCentAxis(1,0,100),
      fIPzAxis(20,-10,10),
      fEtaAxis(80,-2,2),
      fPhiAxis(200,0,TMath::TwoPi()),
      fCent(0),
      fIPz(0),
      fEtaPhi(0),
      fUsedClusters0(0),
      fUsedClusters1(0),
      fAllClusters0(0),
      fAllClusters1(0),
      fStatus(0),
      fRecMode(kNormal|kInjection),
      fScaleDTheta(false),
      fDPhiShift(0.0045),
      fShiftedDPhiCut(-1),
      fScaledDThetaCut(-1),      
      fDeltaCut(1.5),
      fMaxDelta(25),
      fTailDelta(5),
      fDThetaWindow(0.025),
      fDPhiWindow(0),
      fPhiOverlapCut(0.005),
      fZEtaOverlapCut(0.05),
      fPhiRotation(TMath::Pi())
  {}
  /** 
   * Assignment operator 
   * 
   * @param o Object to assign from 
   * 
   * @return Reference to this object 
   */
  AliTrackletdNdetaTask& operator=(const AliTrackletdNdetaTask& o)
  {
    if (&o == this) return *this;
    return *this;
  }
  /** @} */
  //__________________________________________________________________
  /** 
   * Initialize our centrality bins.  
   * 
   * @param existing If null, also create new sub-sets.  If non-null,
   * read sub-sets from the container passed (using in Terminate).
   */
  virtual void InitCentBins(Container* existing);
  /** 
   * Make sure CDB is initialized 
   * 
   * @return true on success
   */
  virtual Bool_t InitCDB();
  /** 
   * Check that we have an initialized geometry if needed
   * 
   * @return true if all is good 
   */
  virtual Bool_t InitGeometry();
  /** 
   * Get the CDB reference run number 
   * 
   * @return A run from LHC10h
   */
  virtual Int_t GetCDBReferenceRun() const { return 137161; }
  /** 
   * Get the CDB reference URL 
   * 
   * @return A fixed string pointing to 2010 
   */
  virtual const char* GetCDBReferenceURL() const
  {
    return "alien://Folder=/alice/data/2010/OCDB";
  }
  /** 
   * @{ 
   * @name Interface 
   */
  virtual void WorkerInit();
  /** 
   * Called at end of worker job 
   * 
   */
  virtual void WorkerFinalize();
  /** 
   * Called on master on merged result. 
   * 
   * @param results Output container for results 
   */
  virtual void MasterFinalize(Container* results);
  /* @} */
  /** 
   * @{ 
   * @name Various worker functions 
   */
  /** 
   * Find cluster (rec.point) tree - if needed
   * 
   * @return Pointer to tree or null
   */
  virtual TTree* FindClusters(Bool_t needed);
  /** 
   * Find the centrality of the event 
   * 
   * @param event Event 
   * 
   * @return Centrality percentile, or negative number in case of
   * problems
   */
  virtual Double_t FindCentrality(AliVEvent* event);
  /** 
   * Find the interaction point location 
   * 
   * @param event Event 
   * 
   * @return Pointer to vertex, or null in case of problems 
   */
  virtual const AliVVertex* FindIP(AliVEvent* event,
				   Double_t maxDispersion=0.04,
				   Double_t maxZError=0.25);
  /** 
   * Check if we got a selected trigger 
   * 
   * @return true if the event was triggered 
   */
  virtual Bool_t FindTrigger();
  /** 
   * Get weight of a tracklet (for MC)
   * 
   * @param trackletNumber Tracklet number
   * @param mult           Tracklet container
   * 
   * @return The weight 
   */
  virtual Double_t TrackletWeight(Int_t                  trackletNumber,
				  const AliMultiplicity* mult) const
  {
    return 1;
  }
  /** 
   * Process tracklets from a reconstructon. 
   * 
   * @param toRun Centrality bins to process 
   * @param mode  Which reconstruction 
   * @param cent  The centrality 
   * @param ipZ   The Z coordinate of the vertex 
   * @param mult  The tracklet container  
   * 
   * @return true if any bin processed. 
   */
  virtual Bool_t ProcessTracklets(TList&            toRun,
				  AliITSMultRecBg*  reco,
				  Double_t          cent,
				  Double_t          ipZ,
				  AliMultiplicity*  mult);
  /** 
   * Find the reconstruction mode associated with reco object 
   * 
   * @param reco Reconstruction object 
   * 
   * @return Mode, or large number if invalid 
   */
  UShort_t FindMode(AliITSMultRecBg* reco) const;
  /** 
   * Fill selected centrality bins with tracklet information 
   * 
   * @param toRun           Selected centrality bins 
   * @param mode            What we're doing 
   * @param mult            Container of tracklets 
   * @param trackletNumber  Tracklet number
   * @param isSignal        If this is considered a signal 
   * @param ipz             Reconstructed IP z coordinate 
   * @param eta             Pseudorapidity of tracklet 
   * @param dPhi            Azimuthal opening angle  
   * @param dPhiS           Shifted azimuthal opening angle 
   * @param dTheta          Polar opening angle 
   * @param dThetaX         Scaled polar opening angle 
   * @param delta           Tracklet @f$\chi^2@f$ 
   * @param weight          Tracklet weight (for MC)
   */
  virtual void FillBins(TList&           toRun,
			AliITSMultRecBg* reco,
			AliMultiplicity* mult,
			Int_t            trackletNumber,
			Bool_t           isSignal,
			Double_t         ipz,
			Double_t         eta,
			Double_t         dPhi,
			Double_t         dPhiS,
			Double_t         dTheta,
			Double_t         dThetaX,
			Double_t         delta,
			Double_t         weight) const;
  /** 
   * (Re-)Run reconstruction from clusters.  
   * 
   * @param toRun     Centrality bins to process 
   * @param mode      Reconstruction mode 
   * @param clusters  Cluster tree
   * @param cent      Centrality 
   * @param ip        Coordinates of IP
   */
  virtual Bool_t Reconstruct(TList&            toRun,
			     UShort_t          mode,
			     TTree*            clusters,
			     Double_t          cent,
			     const AliVVertex* ip);
  /** 
   * Fill histograms of used and all clusters in both layers using
   * information from reconstructor directly.
   * 
   * @param reco Reconstructor 
   */
  virtual void FillClusters(AliITSMultRecBg* reco);
  /** 
   * Fill histograms of used and all clusters in both layers using
   * information from tracklets directly.
   * 
   * @param mult Tracklet container 
   */
  virtual void FillClusters(AliMultiplicity* mult, Double_t ipZ);
  /** 
   * Check the event.  Results are returned with no side-effects other
   * than to fill diagnostics histograms.
   * 
   * @param cent      On return, the centrality 
   * @param ip        On return, the IP coordinates 
   * @param mult      On return, the ESD multiplicity 
   * @param clusters  On return, the clusters (if needed)
   * 
   * @return true if event is OK 
   */
  virtual Bool_t CheckEvent(Double_t&          cent,
			    const AliVVertex*& ip,
			    AliMultiplicity*&  mult,
			    TTree*&            clusters);  
  /** 
   * Do the actual event processing.  If we get here, then the event
   * passed the physics and centrality selections.
   * 
   * @param cent      Centrality 
   * @param ip        Interaction point 
   * @param mult      Reconstructed multiplicity (tracklets)
   * @param clusters  Clusters (if needed)
   */
  virtual void ProcessEvent(Double_t          cent,
			    const AliVVertex* ip,
			    AliMultiplicity*  mult,
			    TTree*            clusters);
  /** 
   * Make a centrality bin 
   * 
   * @param c1 Low edge 
   * @param c2 High edge 
   * 
   * @return 
   */
  virtual CentBin* MakeCentBin(Float_t c1, Float_t c2, UShort_t recFlags)
  {
    return new CentBin(c1, c2, recFlags);
  }
  /* @} */
  //__________________________________________________________________
  /** 
   * @{ 
   * @name 
   */
  Container* fContainer;
  Container* fCentBins;
  TString    fCentMethod;
  TAxis      fCentAxis;
  TAxis      fIPzAxis;
  TAxis      fEtaAxis;
  TAxis      fPhiAxis;
  TH1*       fCent;
  TH1*       fIPz;
  TH2*       fEtaPhi;
  TH2*       fUsedClusters0;
  TH2*       fUsedClusters1;
  TH2*       fAllClusters0;
  TH2*       fAllClusters1;
  TH1*       fStatus;
  UShort_t   fRecMode;
  Bool_t     fScaleDTheta;
  Double_t   fMaxDelta;
  Double_t   fTailDelta;
  Double_t   fDPhiShift;
  Double_t   fShiftedDPhiCut;
  Double_t   fScaledDThetaCut;
  Double_t   fDeltaCut;
  Double_t   fDThetaWindow;
  Double_t   fDPhiWindow;
  Double_t   fPhiOverlapCut;
  Double_t   fZEtaOverlapCut;
  Double_t   fPhiRotation;
  /* @} */

  ClassDef(AliTrackletdNdetaTask,1); // Base class for dN/deta u/tracklets
};


//====================================================================
// Sub base member functions
//____________________________________________________________________
AliTrackletdNdetaTask::Container*
AliTrackletdNdetaTask::SubBase::CreateContainer(Container& parent)
{
  Container* c = new Container;
  c->SetName(Name());
  c->SetOwner();
  parent.Add(c);
  return c;
}

//====================================================================
// Histogram set member functions
//____________________________________________________________________
void AliTrackletdNdetaTask::HistoSet::WorkerInit(Container&     parent,
						 const TAxis&   etaAxis,
						 const TAxis&   ipzAxis,
						 const TAxis&   deltaAxis,
						 const TAxis&   dThetaAxis,
						 const TAxis&   dPhiAxis)
{
  fContainer   = CreateContainer(parent);
  Container* c = fContainer;
      
  fEtaVsIPz      = Make2D(*c, "etaVsIPz",
			  Form("#eta vs IP_{#it{z}} %s",Name()),
			  fColor,fStyle, etaAxis, ipzAxis);
  fEtaVsDelta    = Make2D(*c, "etaVsDelta","",fColor,fStyle,etaAxis,deltaAxis);
  fDPhiVsdTheta  = Make2D(*c, "dphiVsdtheta","",fColor,fStyle,
			  dPhiAxis,dThetaAxis);	
  fDPhiVsdThetaX = Make2D(*c, "dphiVsdthetaX","",fColor,fStyle,
			  dPhiAxis,dThetaAxis);
  fDPhiVsdThetaX->SetXTitle("#Delta#phi-#delta#phi");
  fDPhiVsdThetaX->SetYTitle("#Delta#theta/sin^{2}(#theta)");      
  fDelta         = Make1D(*c, "delta", Form("#Delta %s", Name()),
			  fColor, fStyle,deltaAxis);
}
//____________________________________________________________________
void AliTrackletdNdetaTask::HistoSet::FinalizeInit(Container* parent)
{
  fContainer     = GetC(parent, Name());
  fEtaVsDelta    = GetH2(fContainer, "etaVsDelta");
  fDPhiVsdThetaX = GetH2(fContainer, "dphiVsdthetaX");
  fDPhiVsdTheta  = GetH2(fContainer, "dphiVsdtheta");
  fEtaVsIPz      = GetH2(fContainer, "etaVsIPz");
  fDelta         = GetH1(fContainer, "delta");
}
//____________________________________________________________________
AliTrackletdNdetaTask::Container*
AliTrackletdNdetaTask::HistoSet::MasterFinalize(Container* parent,
						TH1*       ipz,
						Double_t   deltaCut,
						Double_t   tailDelta)
{
  Container* result = CreateContainer(*parent);

  // Get the number of events
  Double_t nEvents = ipz->GetEntries();
  
  // Scale each vertex range by number of events in that range
  TH2* etaVsIPz = ScaleToIPz(fEtaVsIPz, ipz);
  result->Add(etaVsIPz);

  // Normalize delta distribution to integral over 5 to 25 
  TH1* delta = static_cast<TH1*>(CloneAndAdd(result, fDelta));
  delta->Scale(1./nEvents);
  Double_t maxDelta = fDelta->GetXaxis()->GetXmax();
  Double_t eintg;
  Double_t intg = Integrate(delta, tailDelta, maxDelta,  eintg);
      
  result->Add(new TParameter<double>("deltaTailIntg", intg));
  result->Add(new TParameter<double>("deltaTailIntE", eintg));

  TH2* etaVsDelta    = static_cast<TH2*>(CloneAndAdd(result,fEtaVsDelta));
  TH2* dPhiVsdThetaX = static_cast<TH2*>(CloneAndAdd(result,fDPhiVsdThetaX));
  TH2* dPhiVsdTheta  = static_cast<TH2*>(CloneAndAdd(result,fDPhiVsdTheta));
  etaVsDelta   ->Scale(1./nEvents);
  dPhiVsdThetaX->Scale(1./nEvents);
  dPhiVsdTheta ->Scale(1./nEvents);
  
  TH1* deltaIntg = etaVsDelta->ProjectionX("deltaTailInt");
  deltaIntg->SetDirectory(0);
  deltaIntg->Reset();
  deltaIntg->SetYTitle(Form("#int_{%3.1f}^{%4.1f}d#Delta",
			    tailDelta, maxDelta));
  for (Int_t i = 1; i <= deltaIntg->GetNbinsX(); i++) {
    TH1* tmp = etaVsDelta->ProjectionY("tmp",i,i);
    intg     = Integrate(tmp, tailDelta, maxDelta, eintg);
    deltaIntg->SetBinContent(i, intg);
    deltaIntg->SetBinError  (i, eintg);
    delete tmp;
  }
  result->Add(deltaIntg);

  return result;
}
//____________________________________________________________________
void AliTrackletdNdetaTask::HistoSet::Fill(Double_t ipZ,
					   Double_t eta,
					   Double_t dphi,
					   Double_t dphis,
					   Double_t dtheta,
					   Double_t dthetax,
					   Double_t delta,
					   Double_t weight)
{
  if (delta > fEtaVsDelta->GetYaxis()->GetXmax())   return;
  fEtaVsDelta               ->Fill(eta,   delta,   weight);
  fDPhiVsdThetaX            ->Fill(dphis, dthetax, weight);
  fDPhiVsdTheta             ->Fill(dphi,  dtheta,  weight);
  fDelta                    ->Fill(delta,          weight);
  if (ipZ > -999.) fEtaVsIPz->Fill(eta,   ipZ,     weight);
}
//====================================================================
// Centrality bin member functions
//____________________________________________________________________
AliTrackletdNdetaTask::CentBin::CentBin(Float_t  cmin,
					Float_t  cmax,
					UShort_t recFlags,
					Int_t    colOff,
					Int_t    styOff)
  : SubBase(), 
    fMin(cmin),
    fMax(cmax),
    fIPz(0),
    fHistoSets(0),
    fDataSet(0),
    fInjSet(0),
    fRotSet(0)
{
  fHistoSets = new TList;
  fHistoSets->SetOwner(true);

  fDataSet = new HistoSet("data",       kRed+colOff, 20+styOff);
  fHistoSets->Add(fDataSet);
  if (recFlags & kInjection) {
    fInjSet = new HistoSet("injection", kGreen+colOff, 21+styOff);
    fHistoSets->Add(fInjSet);
  }
  if (recFlags & kRotation) {	
    fRotSet = new HistoSet("rotation",  kBlue+colOff, 22+styOff);
    fHistoSets->Add(fRotSet);
  }
}					  
//____________________________________________________________________
const char* AliTrackletdNdetaTask::CentBin::Name() const
{
  return Form("%s%03dd%02d_%03dd%02d", "cent",
	      Int_t(fMin), Int_t(fMin*100) % 100, 
	      Int_t(fMax), Int_t(fMax*100) % 100);
}
//____________________________________________________________________
void AliTrackletdNdetaTask::CentBin::WorkerInit(Container&       parent,
						const TAxis&     etaAxis,
						const TAxis&     ipzAxis,
						const TAxis&     deltaAxis,
						const TAxis&     dThetaAxis,
						const TAxis&     dPhiAxis)
{
  fContainer = CreateContainer(parent);
  fIPz = Make1D(*fContainer, "ipz",  "", kMagenta+2, 20, ipzAxis);
  TIter next(fHistoSets);
  HistoSet* set = 0;
  while ((set = static_cast<HistoSet*>(next())))
    set->WorkerInit(*fContainer,etaAxis,ipzAxis,deltaAxis,dThetaAxis,dPhiAxis);
}
//____________________________________________________________________
void AliTrackletdNdetaTask::CentBin::FinalizeInit(Container* parent)
{
  fContainer     = GetC(parent, Name());
  fIPz           = GetH1(fContainer, "ipz");
  TIter next(fHistoSets);
  HistoSet* set = 0;
  while ((set = static_cast<HistoSet*>(next())))
    set->FinalizeInit(fContainer);
}
//____________________________________________________________________
AliTrackletdNdetaTask::Container*
AliTrackletdNdetaTask::CentBin::EstimateBackground(Container* dataCon,
						   HistoSet*  set,
						   Container* result,
						   Double_t   deltaCut,
						   Double_t   deltaTail)
{
  if (!set || !dataCon) return 0;
  
  Container* bgCon = set->MasterFinalize(result, fIPz, deltaCut, deltaTail);
  if (!bgCon) return 0;

  Bool_t isComb = set->fName.Contains("combinatorics", TString::kIgnoreCase);
  // Info("EstimateBackground", "%s - combinatorics? %s",
  //      set->fName.Data(), isComb ? "yes" : "no");

  Double_t dataIntg = GetD(dataCon, "deltaTailIntg", -1);
  Double_t dataIntE = GetD(dataCon, "deltaTailIntE", -1);
  Double_t bgIntg   = GetD(bgCon,   "deltaTailIntg", -1);
  Double_t bgIntE   = GetD(bgCon,   "deltaTailIntE", -1);
  if (dataIntg <= 0 || bgIntg <= 0) return bgCon;
  if (lDebug > 2)  {
    AliInfoF("int(measured) Delta tail: %f +/- %f", dataIntg, dataIntE);
    AliInfoF("int(%s) Delta tail: %f +/- %f", bgCon->GetName(),
	     bgIntg, bgIntE);
  }
  Double_t scaleE   = 0;
  Double_t scale    = RatioE(dataIntg, dataIntE, bgIntg, bgIntE, scaleE);

  TH1* bgDelta = GetH1(bgCon, "delta");
  if (!bgDelta) return bgCon;

  TH1* dataDelta = GetH1(dataCon, "delta");
  if (!dataDelta) {
    AliWarningF("No delta distribution in %s", dataCon->GetName());
  }

  TH1* deltaScaled = CopyH1(bgCon, "delta", "deltaScaled");
  bgCon->Add(deltaScaled);
  deltaScaled->SetLineStyle(7);
  deltaScaled->SetTitle(Form("%s#times%6.3f", bgDelta->GetTitle(), scale));
  if (dataDelta) deltaScaled->SetMarkerStyle(dataDelta->GetMarkerStyle());
  else           deltaScaled->SetMarkerStyle(20);

  if (lDebug > 2)
    AliInfoF("Scaling %s Delta dist by %f +/- %f",
	     bgCon->GetName(), scale, scaleE);
  Scale(deltaScaled, scale, scaleE);

  TH2* backgroundEst = CopyH2(bgCon,  "etaVsIPz", "backgroundEst");
  backgroundEst->SetTitle("Background");
  if (!isComb) Scale(backgroundEst, scale, scaleE); // ->Scale(scale);
  // else         Info("EstimateBackground", "Combinators, no scaling of BG");
  bgCon->Add(backgroundEst);

  TH2* signalEst = CopyH2(dataCon, "etaVsIPz","signalEst");
  signalEst->SetMarkerStyle(backgroundEst->GetMarkerStyle());
  signalEst->SetMarkerColor(backgroundEst->GetMarkerColor());
  signalEst->SetMarkerSize (backgroundEst->GetMarkerSize());
  signalEst->SetLineStyle  (backgroundEst->GetLineStyle());
  signalEst->SetLineColor  (backgroundEst->GetLineColor());
  signalEst->SetLineWidth  (backgroundEst->GetLineWidth());
  signalEst->SetFillStyle  (backgroundEst->GetFillStyle());
  signalEst->SetFillColor  (backgroundEst->GetFillColor());
  signalEst->SetTitle("Signal");
  signalEst->Add(backgroundEst,-1);
  bgCon->Add(signalEst);

  TH2* measured  = GetH2(dataCon, "etaVsIPz");
  TH2* beta      = static_cast<TH2*>(backgroundEst->Clone("beta"));
  beta->SetTitle("#beta");
  beta->SetDirectory(0);
  beta->Divide(measured);
  bgCon->Add(beta);
    
  TH2* bg1MBeta   = static_cast<TH2*>(beta->Clone("oneMinusBeta"));
  bg1MBeta->SetTitle("1-#beta");
  bg1MBeta->SetDirectory(0);
  bg1MBeta->Reset();
  bgCon->Add(bg1MBeta);  
  for (Int_t ix = 1; ix <= bg1MBeta->GetNbinsX(); ix++) {
    for (Int_t iy = 1; iy <= bg1MBeta->GetNbinsY(); iy++) {
      bg1MBeta->SetBinContent(ix,iy,1);
      bg1MBeta->SetBinError  (ix,iy,0);
    }
  }
  bg1MBeta->Add(beta,-1);

  // if (lDebug > 2) bgCon->ls();
  return bgCon;
}
//____________________________________________________________________
AliTrackletdNdetaTask::Container*
AliTrackletdNdetaTask::CentBin::MasterFinalize(Container* parent,
					       TH1*       ipz,
					       Double_t   deltaCut,
					       Double_t   deltaTail)
{
  Double_t   nEvents = fIPz->GetEntries();
  Container* result  = CreateContainer(*parent);
  Printf("Bin %5.1f-%5.1f%%: %10d events", fMin, fMax, Int_t(nEvents));

  // Copy ipZ histogram and scale by number of events 
  TH1* ipZ = static_cast<TH1*>(CloneAndAdd(result, fIPz));
  ipZ->Scale(1./nEvents);
  
  Container* dataRes = 0;
  TIter      next(fHistoSets);
  HistoSet*  set = 0;
  while ((set = static_cast<HistoSet*>(next()))) {
    if (set == fDataSet) {
      dataRes = set->MasterFinalize(result, fIPz, deltaCut, deltaTail);
      continue;
    }
    EstimateBackground(dataRes, set, result, deltaCut, deltaTail);
  }

  if (lDebug > 5) {
    result->ls();
  }
  return result;
}
//____________________________________________________________________
Bool_t AliTrackletdNdetaTask::CentBin::FillSet(UShort_t set,
					       Bool_t   isSignal,
					       Double_t ipZ,
					       Double_t eta,
					       Double_t dPhi,
					       Double_t dPhiS,
					       Double_t dTheta,
					       Double_t dThetaX,
					       Double_t delta,
					       Double_t weight)
{
  HistoSet* hset = 0;
  if      (set == kNormal)    hset = fDataSet;
  else if (set == kInjection) hset = fInjSet;
  else if (set == kRotation)  hset = fRotSet;

  if (!hset) return false;
      
  hset->Fill(isSignal ? ipZ : -1000,eta,dPhi,dPhiS,dTheta,dThetaX,delta,weight);

  return true;
}

//====================================================================
// Task member functions
//____________________________________________________________________
AliTrackletdNdetaTask::AliTrackletdNdetaTask(const char* name)
  : AliAnalysisTaskSE(name),
    AliTrackletdNdetaUtils(),
    fContainer(0),
    fCentBins(0),
    fCentMethod("V0M"),
    fCentAxis(1,0,100),
    fIPzAxis(20,-10,10),
    fEtaAxis(80,-2,2),
    fPhiAxis(200,0,TMath::TwoPi()),
    fCent(0),
    fIPz(0),
    fEtaPhi(0),
    fUsedClusters0(0),
    fUsedClusters1(0),
    fAllClusters0(0),
    fAllClusters1(0),
    fStatus(0),
    fRecMode(kNormal|kInjection),
    fScaleDTheta(false),
    fDPhiShift(0.0045),
    fShiftedDPhiCut(-1),
    fScaledDThetaCut(-1),
    fDeltaCut(1.5),
    fMaxDelta(25),
    fTailDelta(5),
    fDThetaWindow(0.025),
    fDPhiWindow(0.06),
    fPhiOverlapCut(0.005),
    fZEtaOverlapCut(0.05),
    fPhiRotation(TMath::Pi())
{
  FixAxis(fCentAxis, "Centrality [%]");
  FixAxis(fIPzAxis,  "IP_{#it{z}} [cm]");
  FixAxis(fEtaAxis,  "#eta");
  FixAxis(fPhiAxis,  "#phi");

  // Set branches we want to look for 
  fBranchNames = 
    "ESD:AliESDRun.,AliESDHeader.,AliMultiplicity.,"
    "SPDVertex.,PrimaryVertex.";

  DefineOutput(1, Container::Class());
  DefineOutput(2, Container::Class());
}
//____________________________________________________________________
void AliTrackletdNdetaTask::UserCreateOutputObjects()
{
  WorkerInit();
  lDebug = DebugLevel();
  PostData(1, fContainer);
}
//____________________________________________________________________
void AliTrackletdNdetaTask::UserExec(Option_t*)
{
  Double_t          cent     = -1;
  const AliVVertex* ip       = 0;
  AliMultiplicity*  mult     = 0;
  TTree*            clusters = 0;
  if (!CheckEvent(cent, ip, mult, clusters)) {
    AliWarningF("Event didn't pass %f, %p, %p, %p",
		cent, ip, mult, clusters);
    return;
  }
    
  ProcessEvent(cent, ip, mult, clusters);

  PostData(1,fContainer);

  fStatus->Fill(kCompleted);
}
//____________________________________________________________________
void AliTrackletdNdetaTask::Terminate(Option_t*)
{
  lDebug = DebugLevel();

  Container* results = new Container;
  results->SetName(Form("%sResults",GetName()));
  results->SetOwner();

  Print("");
  fContainer = static_cast<Container*>(GetOutputData(1));
  if (!fContainer) {
    AliWarning("No sum container found!");
    return;
  }
    
  InitCentBins(fContainer);
    
  MasterFinalize(results);
    
  PostData(2, results);
}
//____________________________________________________________________
Bool_t AliTrackletdNdetaTask::Connect(const char* sumFile,
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
void AliTrackletdNdetaTask::Print(Option_t*) const
{
  Double_t shiftedDPhiCut = fShiftedDPhiCut;
  if (shiftedDPhiCut < 0)
    shiftedDPhiCut = TMath::Sqrt(fDeltaCut)*fDPhiWindow;
  
  Printf("%s: %s", ClassName(), GetName());
  Printf(" %22s: 0x%x", "Reconstruction mode",     fRecMode);
  Printf(" %22s: %d",   "Scale by sin^2(theta)",   fScaleDTheta);
  Printf(" %22s: %f",   "Delta phi shift",	   fDPhiShift);
  Printf(" %22s: %f",   "Shifted Delta phi cut",   shiftedDPhiCut);
  Printf(" %22s: %f",   "Scaled Delta Theta cut",  fScaledDThetaCut);
  Printf(" %22s: %f",   "Delta cut",	           fDeltaCut);
  Printf(" %22s: %f",   "max Delta",	           fMaxDelta);
  Printf(" %22s: %f",   "tail Delta",	           fTailDelta);
  Printf(" %22s: %f",   "Delta theta window",	   fDThetaWindow);
  Printf(" %22s: %f",   "Delta phi window",	   fDPhiWindow);
  Printf(" %22s: %f",   "phi overlap cut",	   fPhiOverlapCut);
  Printf(" %22s: %f",   "z-eta overlap cut",	   fZEtaOverlapCut);
  Printf(" %22s: %f",   "phi rotation",            fPhiRotation);
  PrintAxis(fEtaAxis);
  PrintAxis(fPhiAxis);
  PrintAxis(fIPzAxis,1,"IPz");
  PrintAxis(fCentAxis,0);
}
//____________________________________________________________________
void AliTrackletdNdetaTask::InitCentBins(Container* existing)
{
  if (fCentBins) return;
  fCentBins = new Container;
  fCentBins->SetName("centralityBins");
  fCentBins->SetOwner();

  Double_t maxdPhi   = fDPhiWindow  *TMath::Sqrt(fMaxDelta);
  Double_t maxdTheta = fDThetaWindow*TMath::Sqrt(fMaxDelta);
  TAxis    deltaAxis (Int_t(5*fMaxDelta),0,         fMaxDelta);
  TAxis    dThetaAxis(100,               -maxdTheta,+maxdTheta);
  TAxis    dPhiAxis  (100,               -maxdPhi,  +maxdPhi);
  FixAxis(deltaAxis,
	  Form("#Delta=[(#Delta#phi-#delta#phi)/#sigma_{#phi}]^{2}+"
	       "[#Delta#theta%s/#sigma_{#theta}]^{2}",
	       fScaleDTheta ? "sin^{-2}(#theta)" : ""));
  FixAxis(dThetaAxis, "#Delta#theta");
  FixAxis(dPhiAxis,   "#Delta#phi");

  // Add min-bias bin 
  CentBin* bin  = MakeCentBin(0, 100, fRecMode);
  if (!existing) 
    bin->WorkerInit(*fContainer,fEtaAxis,fIPzAxis,
		    deltaAxis, dThetaAxis, dPhiAxis);
  else
    bin->FinalizeInit(existing);
  fCentBins->AddAt(bin, 0);

  // Add other bins
  Int_t nCentBins = fCentAxis.GetNbins();
  for (Int_t i = 1; i <= nCentBins; i++) {
    Float_t  c1 = fCentAxis.GetBinLowEdge(i);
    Float_t  c2 = fCentAxis.GetBinUpEdge(i);
    bin         = MakeCentBin(c1, c2, fRecMode);
    if (!existing) 
      bin->WorkerInit(*fContainer,fEtaAxis,fIPzAxis,
		      deltaAxis, dThetaAxis, dPhiAxis);
    else
      bin->FinalizeInit(existing);
    fCentBins->AddAt(bin, i);
  }
}
//____________________________________________________________________
Bool_t AliTrackletdNdetaTask::InitCDB()
{
  Printf("Now initialising CDB");
  AliAnalysisManager* anaMgr = AliAnalysisManager::GetAnalysisManager();
  if (!anaMgr) {
    AliError("No manager defined!");
    return false;
  }
  // Check if we have the CDB connect task, and if so, do nothing as
  // we rely on that task set up things properly.
  const char*  cdbNames[] = {"CDBconnect", "cdb", 0 };
  const char** ptr        = cdbNames;
  while (*ptr) { 
    AliAnalysisTask* cdbConnect = anaMgr->GetTask(*ptr);
    if (cdbConnect && cdbConnect->IsA()->InheritsFrom("AliTaskCDBconnect")) {
      AliInfoF("CDB-connect task (%s: %s) present, do nothing",
	       cdbConnect->ClassName(), *ptr);
      return true;
    }
    ptr++;
  }
  // Otherwise, we need to do stuff ourselves
  Printf("Get the CDB manager");
  AliCDBManager* cdbMgr = AliCDBManager::Instance();
  if (!cdbMgr) {
    AliError("Failed to get instance of CDB manager");
    return false;
  }
  Int_t   refRun = GetCDBReferenceRun();
  TString refUrl = GetCDBReferenceURL();
  AliWarningF("Using reference CDB storage \"%s\" and run \"%d\"",
	      refUrl.Data(), refRun);
  // Set a very particular default storage. Perhaps we can do this
  // with specific storages instead!  Depends on whether the
  // reconstruction also uses CDB - probably does - in which case we
  // need specific storages for that too.
  cdbMgr->SetDefaultStorage(refUrl);
  // Now load our geometry - from a LHC10h run 
  Printf("Get Geometry entry");
  AliCDBEntry* cdbEnt = cdbMgr->Get("GRP/Geometry/Data", refRun);
  if (!cdbEnt) {
    AliErrorF("No geometry found from %d", refRun);
    return false;
  }
  // Initialize the geometry manager
  Printf("Set Geometry");
  AliGeomManager::SetGeometry(static_cast<TGeoManager*>(cdbEnt->GetObject()));
  // Now perform mis-alignment - again based on an LHC10h run!
  Printf("Misalign geometry");
  if (!AliGeomManager::ApplyAlignObjsToGeom("ITS",refRun,-1,-1)) {
    AliErrorF("Failed to misalign geometry from %d", refRun);
    return false;
  }
  return true;
}
//____________________________________________________________________
Bool_t AliTrackletdNdetaTask::InitGeometry()
{
  if (fRecMode == 0) return true;

  if (!AliGeomManager::GetGeometry()) {
    AliError("No geometry loaded, needed for reconstruction");
    AliError("Add the AliTaskCDBconnect to the train");
    return false;
  }
  return true;
}
//____________________________________________________________________
void AliTrackletdNdetaTask::WorkerInit()
{
  if (!InitCDB()) {
    // Argh! - make us a a zombie
    SetZombie(true);
    return;
  }
  
  if (fShiftedDPhiCut < 0)
    fShiftedDPhiCut = TMath::Sqrt(fDeltaCut)*fDPhiWindow;
    
  fContainer = new Container;
  fContainer->SetName(Form("%sSums",GetName()));
  fContainer->SetOwner();
    
  fIPz    = Make1D(*fContainer, "ipz",  "", kMagenta+2, 20, fIPzAxis);
  fCent   = Make1D(*fContainer, "cent", "", kMagenta+2, 20, fCentAxis);
  fEtaPhi = Make2D(*fContainer, "etaPhi","",kMagenta+2, 20, fEtaAxis,fPhiAxis);
  fIPz   ->SetFillStyle(0);
  fCent  ->SetFillStyle(0);

  TAxis zAxis(64,-16,16);     FixAxis(zAxis, "#it{z} [cm]");
  TAxis phiAxis(80, 0,TMath::TwoPi()); FixAxis(phiAxis,"#phi");
  TAxis ph2Axis(160,0,TMath::TwoPi()); FixAxis(ph2Axis,"#phi");

  fUsedClusters0 = Make2D(*fContainer, "usedClusters0",
			  "Used clusters on layer 0",
			  kBlack, 1, zAxis, phiAxis);
  fUsedClusters1 = Make2D(*fContainer, "usedClusters1",
			  "Used clusters on layer 1",
			  kBlack, 1, zAxis, ph2Axis);
  fAllClusters0 = Make2D(*fContainer, "allClusters0",
			 "All clusters on layer 0",
			 kBlack, 1, zAxis, phiAxis);
  fAllClusters1 = Make2D(*fContainer, "allClusters1",
			 "All clusters on layer 1",
			 kBlack, 1, zAxis, ph2Axis);
    
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
  fStatus->GetXaxis()->SetBinLabel(kMultiplicity,   "Multiplicity");
  fStatus->GetXaxis()->SetBinLabel(kClusters,       "Clusters");
  fStatus->GetXaxis()->SetBinLabel(kTrigger,        "Trigger");
  fStatus->GetXaxis()->SetBinLabel(kIP,             "IP");
  fStatus->GetXaxis()->SetBinLabel(kCentrality,     "Centrality");
  fStatus->GetXaxis()->SetBinLabel(kReconstructed,  "Been reconstructed");
  fStatus->GetXaxis()->SetBinLabel(kInjected,       "Done Injection");
  fStatus->GetXaxis()->SetBinLabel(kRotated,        "Been rotated");
  fStatus->GetXaxis()->SetBinLabel(kCompleted,      "Completed");
  fContainer->Add(fStatus);

  typedef TParameter<double> DP;
  typedef TParameter<bool>   BP;
  typedef TParameter<int>    IP;
  Container* params = new Container;
  params->SetName("parameters");
  params->SetOwner();
  fContainer->Add(params);
  params->Add(new IP("RecMode",        fRecMode,        'f'));
  params->Add(new BP("ScaleDTheta",    fScaleDTheta,    'f'));
  params->Add(new DP("DPhiShift",      fDPhiShift,      'f'));
  params->Add(new DP("ShiftedDPhiCut", fShiftedDPhiCut, 'f'));
  params->Add(new DP("ScaledDThetaCut",fScaledDThetaCut,'f'));
  params->Add(new DP("DeltaCut",       fDeltaCut,       'f'));
  params->Add(new DP("MaxDelta",       fMaxDelta,       'f'));
  params->Add(new DP("TailDelta",      fTailDelta,      'f'));
  params->Add(new DP("DThetaWindow",   fDThetaWindow,   'f'));
  params->Add(new DP("DPhiWindow",     fDPhiWindow,     'f'));
  params->Add(new DP("PhiOverlapCut",  fPhiOverlapCut,  'f'));
  params->Add(new DP("ZEtaOverlapCut" ,fZEtaOverlapCut, 'f'));
  params->Add(new DP("PhiRotation",    fPhiRotation,    'f'));

  // Fix some fill styles 

  InitCentBins(0);
}
//____________________________________________________________________
void AliTrackletdNdetaTask::WorkerFinalize()
{
  Printf("This worker had %f events within cuts",
	 fIPz->GetEntries());
}
//____________________________________________________________________
void AliTrackletdNdetaTask::MasterFinalize(Container* results)
{
  // Make copies of histograms and store
  fIPz    = static_cast<TH1*>(CloneAndAdd(results,GetH1(fContainer,"ipz")));
  fStatus = static_cast<TH1*>(CloneAndAdd(results,GetH1(fContainer,"status")));
  fCent   = static_cast<TH1*>(CloneAndAdd(results,GetH1(fContainer,"cent")));  

  Double_t nEvents = fIPz->GetEntries();
  Printf("Event summary:");
  for (Int_t i = 1; i <= fStatus->GetNbinsX(); i++) 
    Printf("  %10d %s", Int_t(fStatus->GetBinContent(i)),
	   fStatus->GetXaxis()->GetBinLabel(i));

  fIPz   ->Scale(1./nEvents);
  fStatus->Scale(1./fStatus->GetBinContent(1));
  fCent  ->Scale(1./fCent->GetEntries());

  fUsedClusters0 = static_cast<TH2*>(CloneAndAdd(results,
						 GetH2(fContainer,
						       "usedClusters0")));
  fUsedClusters1 = static_cast<TH2*>(CloneAndAdd(results,
						 GetH2(fContainer,
						       "usedClusters1")));
  fAllClusters0 = static_cast<TH2*>(CloneAndAdd(results,
						 GetH2(fContainer,
						       "allClusters0")));
  fAllClusters1 = static_cast<TH2*>(CloneAndAdd(results,
						 GetH2(fContainer,
						       "allClusters1")));
  fUsedClusters0->Scale(1/nEvents);
  fUsedClusters1->Scale(1/nEvents);
  fAllClusters0 ->Scale(1/nEvents);
  fAllClusters1 ->Scale(1/nEvents);

  TIter    next(fCentBins);
  CentBin* bin = 0;
  while ((bin = static_cast<CentBin*>(next()))) {
    bin->MasterFinalize(results, 0, fDeltaCut, fTailDelta);
  }
}
//____________________________________________________________________
TTree* AliTrackletdNdetaTask::FindClusters(Bool_t needed)
{
  if (!needed) return 0;
  AliAnalysisManager* mgr = AliAnalysisManager::GetAnalysisManager();
  AliVEventHandler*   inh = mgr->GetInputEventHandler();
  if (!inh->IsA()->InheritsFrom(AliESDInputHandlerRP::Class())) {
    AliErrorF("Not the right kind of input handler: %s",
	      inh->ClassName());
    return 0;
  }
  AliESDInputHandlerRP* rph = static_cast<AliESDInputHandlerRP*>(inh);
  TTree* tree = rph->GetTreeR("ITS");
  if (!tree) {
    AliError("Tree of clusters (rec.points) not found");
    return 0;
  }
  return tree;
}
//____________________________________________________________________
Double_t AliTrackletdNdetaTask::FindCentrality(AliVEvent* event)
{
  if (fCentMethod.EqualTo("MB", TString::kIgnoreCase)) return 0;
  AliMultSelection* cent =
    static_cast<AliMultSelection*>(event->FindListObject("MultSelection"));
  if (!cent) {
    AliWarning("No centrality in event");
    return -1;
  }
  const Double_t safety = 1e-3;
  Double_t centPer = cent->GetMultiplicityPercentile(fCentMethod);
  if      (centPer < -safety)    return -2;
  if      (centPer < +safety)    centPer = safety;
  else if (centPer > 100-safety) centPer = 100-safety;

  if (centPer < fCentAxis.GetXmin() || centPer > fCentAxis.GetXmax()) {
    AliWarningF("Centrality = %f out of range [%f,%f]",
		centPer, fCentAxis.GetXmin(), fCentAxis.GetXmax());
    return -3;
  }
  return centPer;    
}
//____________________________________________________________________
const AliVVertex* AliTrackletdNdetaTask::FindIP(AliVEvent* event,
						Double_t maxDispersion,
						Double_t maxZError)
{
  const AliVVertex* ip   = event->GetPrimaryVertex();
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
  // Good vertex, return it
  return ip;
}
//____________________________________________________________________
Bool_t AliTrackletdNdetaTask::FindTrigger()
{
  UInt_t evBits = fInputHandler->IsEventSelected();
  Bool_t trgOK  = evBits & fOfflineTriggerMask;

  return trgOK;
}
//____________________________________________________________________
Bool_t AliTrackletdNdetaTask::ProcessTracklets(TList&            toRun,
					       AliITSMultRecBg*  reco,
					       Double_t          cent,
					       Double_t          ipZ,
					       AliMultiplicity*  mult)
{
  // This can probably be optimized by using the reco object
  // directly. For now, we leave as is and accept the extra
  // processing time it implies.
  if      (!reco)                     FillClusters(mult, ipZ);
  else if (FindMode(reco) == kNormal) FillClusters(reco);

  // Now loop over all tracklets.  Since this method is only called
  // once per "reconstruction" pass, we save a lot of computing time
  // since we loop over the tracklets of each reconstruction exactly
  // once.  It also means that we implement the calculations and
  // cuts in exacly one place (except caveat on filling clusters above). 
  Int_t nTracklets = mult->GetNumberOfTracklets();
  for (Int_t trackletNumber = 0; trackletNumber < nTracklets;
       trackletNumber++) {
    Double_t weight  = TrackletWeight(trackletNumber, mult);
    Double_t theta   = mult->GetTheta(trackletNumber);
    Double_t phi     = mult->GetPhi(trackletNumber);
    Double_t dTheta  = mult->GetDeltaTheta(trackletNumber);
    Double_t dPhi    = mult->GetDeltaPhi(trackletNumber);
    Double_t delta   = mult->CalcDist(trackletNumber);
    Double_t eta     = -TMath::Log(TMath::Tan(theta/2));
    Double_t dThetaX = dTheta;
    Double_t dPhiS   = dPhi - TMath::Sign(fDPhiShift,dPhi); //+ if dPhi<0
    if (fScaleDTheta) dThetaX /= TMath::Power(TMath::Sin(theta),2);

    // Optional cut on scaled DTheta
    if (fScaledDThetaCut >= 0 &&
	TMath::Abs(dThetaX) > fScaledDThetaCut) continue;
    // If not within the bins we're looking at, continue 
    if (eta < fEtaAxis.GetXmin() || eta > fEtaAxis.GetXmax())  continue;
    // if (phi < fPhiAxis.GetXmin() || phi > fPhiAxis.GetXmax())  continue;

    // Assume a signal 
    Bool_t   isSignal = true;
    if (delta >  fDeltaCut)       isSignal = false;
    if (dPhiS >= fShiftedDPhiCut) isSignal = false;

    // Now fill all found centrality bins 
    FillBins(toRun, reco, mult, trackletNumber, isSignal, ipZ,
	     eta, dPhi, dPhiS, dTheta, dThetaX, delta, weight);

    // If this is a signal, fill into eta,phi map 
    if (isSignal) fEtaPhi->Fill(eta,phi);
  }
  return true; // nAcc > 0;
}
//____________________________________________________________________
UShort_t AliTrackletdNdetaTask::FindMode(AliITSMultRecBg* reco) const
{
  if (!reco) return kNormal;
  switch (reco->GetRecType()) {
  case AliITSMultRecBg::kData:  return kNormal;
  case AliITSMultRecBg::kBgInj: return kInjection;
  case AliITSMultRecBg::kBgRot: return kRotation; 
  default:
    AliWarningF("Unsupported reconstruction mode: %d", reco->GetRecType());
    break;
  }
  return 999;
}
//____________________________________________________________________
void AliTrackletdNdetaTask::FillBins(TList&           toRun,
				     AliITSMultRecBg* reco,
				     AliMultiplicity* mult,
				     Int_t            trackletNumber,
				     Bool_t           isSignal,
				     Double_t         ipz,
				     Double_t         eta,
				     Double_t         dPhi,
				     Double_t         dPhiS,
				     Double_t         dTheta,
				     Double_t         dThetaX,
				     Double_t         delta,
				     Double_t         weight) const
{
  UShort_t mode = FindMode(reco);
  if (mode >= 999) return;
    
  TIter nextB(&toRun);
  CentBin* bin = 0;
  while ((bin = static_cast<CentBin*>(nextB())))
    bin->FillSet(mode, isSignal, ipz, eta,
		 dPhi, dPhiS, dTheta, dThetaX, delta, weight);
}    
//____________________________________________________________________
Bool_t AliTrackletdNdetaTask::Reconstruct(TList&            toRun,
					  UShort_t          mode,
					  TTree*            clusters,
					  Double_t          cent,
					  const AliVVertex* ip)
{
  UShort_t firstMode = mode;
  AliITSMultRecBg reco;
  reco.SetCreateClustersCopy        (true);
  reco.SetScaleDThetaBySin2T        (fScaleDTheta);
  reco.SetNStdDev                   (fMaxDelta);
  reco.SetPhiWindow                 (fDPhiWindow);
  reco.SetThetaWindow               (fDThetaWindow);
  reco.SetPhiShift                  (fDPhiShift);
  reco.SetRemoveClustersFromOverlaps(fPhiOverlapCut > 0 ||
				     fZEtaOverlapCut   > 0);
  reco.SetPhiOverlapCut             (fPhiOverlapCut);
  reco.SetZetaOverlapCut            (fZEtaOverlapCut);
  reco.SetHistOn                    (false);
  switch (mode) {
  case kNormal:
  case kInjection:
    // If we're to do injection, first run normally, since the
    // injection procedure relies on this, and there's no reason to
    // do it twice
    firstMode = kNormal;
    reco.SetRecType(AliITSMultRecBg::kData);
    break;
  case kRotation:
    reco.SetRecType(AliITSMultRecBg::kBgRot);
    reco.SetPhiRotationAngle(fPhiRotation);
    break;
  }      
    
  // Now run the reconstruction
  Float_t ipv[3];
  ipv[0] = ip->GetX();
  ipv[1] = ip->GetY();
  ipv[2] = ip->GetZ();
  reco.Run(clusters, ipv);

  // And fill results into centrality bins
  Bool_t ret = ProcessTracklets(toRun, &reco,  cent, ip->GetZ(),
				reco.GetMultiplicity());
   
  // If we're not doing injection, return 
  if (mode != kInjection) return ret;

  // If we do injection, run again, but in injection mode
  reco.SetRecType(AliITSMultRecBg::kBgInj);
  reco.Run(clusters, ipv);

  // And fill results into centrality bins
  ret = ProcessTracklets(toRun, &reco,  cent, ip->GetZ(),
			 reco.GetMultiplicity());

  return ret;
}
//____________________________________________________________________
void AliTrackletdNdetaTask::FillClusters(AliITSMultRecBg* reco)
{
  if (!reco) return;

  TH2*  usedClusters[] = { fUsedClusters0, fUsedClusters1 };
  TH2*  allClusters[]  = { fAllClusters0,  fAllClusters1 };

  // Loop over tracklets
  Int_t nTracklets     = reco->GetNTracklets();
  for (Int_t trackletNo = nTracklets; trackletNo--; ) {
    // Get tracklet parameters
    Float_t* tracklet = reco->GetTracklet(trackletNo);

    Float_t dPhi      = tracklet[AliITSMultReconstructor::kTrDPhi];
    Float_t dTheta    = tracklet[AliITSMultReconstructor::kTrDTheta];
    Float_t theta     = tracklet[AliITSMultReconstructor::kTrTheta];
    if (TMath::Abs(TMath::Abs(dPhi)-fDPhiShift) > fShiftedDPhiCut)
      continue;

    Float_t delta     = reco->CalcDist(dPhi, dTheta, theta);
    if (delta > fDeltaCut) continue;

    Int_t clusterID[] = { Int_t(tracklet[AliITSMultReconstructor::kClID1]),
			  Int_t(tracklet[AliITSMultReconstructor::kClID2]) };
    for (Int_t layer = 0; layer < 2; layer++) {
      // Get cluster information
      Float_t* cluster = reco->GetClusterOfLayer(layer, clusterID[layer]);
      // Fill used cluster z,phi
      usedClusters[layer]->Fill(cluster[AliITSMultReconstructor::kClZ],
				cluster[AliITSMultReconstructor::kClPh]);
    }
  }
  // Loop over all clusters and fill information
  for (Int_t layer = 0; layer < 2; layer++) {
    for (Int_t clusterNo = reco->GetNClustersLayer(layer); clusterNo--; ) {
      Float_t* cluster = reco->GetClusterOfLayer(layer, clusterNo);
      allClusters[layer]->Fill(cluster[AliITSMultReconstructor::kClZ],
			       cluster[AliITSMultReconstructor::kClPh]);
    }
  }
}
//____________________________________________________________________
void AliTrackletdNdetaTask::FillClusters(AliMultiplicity* mult, Double_t ipZ)
{
  if (!mult) return;
  const double kRSPD2 = 3.9;

  // Loop over tracklets 
  Int_t nTracklets = mult->GetNumberOfTracklets();
  for (Int_t trackletNo = nTracklets; trackletNo--; ) {
    Double_t phi  = mult->GetPhi(trackletNo);
    Double_t z    = kRSPD2 / TMath::Tan(mult->GetTheta(trackletNo)) + ipZ;
    Double_t dPhi = mult->GetDeltaPhi(trackletNo);
    if ((TMath::Abs(dPhi-fDPhiShift) <= fShiftedDPhiCut) &&
	(mult->CalcDist(trackletNo  )<= fDeltaCut))
      fUsedClusters0->Fill(z, phi);
    fAllClusters0->Fill(z, phi);
  }
  // Loop over free clusters
  Int_t nClusters = mult->GetNumberOfSingleClusters();
  for (Int_t clusterNo = nClusters; clusterNo--; ) {
    Double_t phi = mult->GetPhiSingle(clusterNo);
    Double_t z   = kRSPD2 / TMath::Tan(mult->GetThetaSingle(clusterNo)) + ipZ;
    fAllClusters0->Fill(z, phi);
  }
}
//____________________________________________________________________
Bool_t AliTrackletdNdetaTask::CheckEvent(Double_t&          cent,
					 const AliVVertex*& ip,
					 AliMultiplicity*&  mult,
					 TTree*&            clusters)
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

  // Check for geometry
  if (!InitGeometry()) {
    // Argh! - make us a a zombie
    SetZombie(true);
    return false;
  }

  // Check magnetic field 
  if (!TGeoGlobalMagField::Instance()->GetField() &&
      !event->InitMagneticField()) {
    AliWarning("Failed to initialize magnetic field");
    return false;
  }
    
  // Check the multiplicity 
  AliVMultiplicity* vmult = event->GetMultiplicity();
  if (!vmult) return false;
  if (!vmult->IsA()->InheritsFrom(AliMultiplicity::Class())) {
    AliWarningF("Multiplicity is a %s object", vmult->ClassName());
    return false;
  }
  mult = static_cast<AliMultiplicity*>(vmult);
  fStatus->Fill(kMultiplicity);
    
  // Check that we have clusters - if so needed 
  Bool_t needRP = (fRecMode != 0);
  clusters      = FindClusters(needRP);
  if (needRP && !clusters) return false;
  if (needRP) fStatus->Fill(kClusters);
    
  // Check if event was triggered 
  Bool_t trg = FindTrigger();
  if (!trg) return false;
  fStatus->Fill(kTrigger);
    
  // Check the interaction point 
  ip = FindIP(event);
  if (!ip) return false;
  fStatus->Fill(kIP);

  // Check the centrality 
  cent = FindCentrality(event);
  if (cent < 0) return false;
  fStatus->Fill(kCentrality);

  fIPz->Fill(ip->GetZ());
  fCent->Fill(cent);
  return true;
}
  
//____________________________________________________________________
void AliTrackletdNdetaTask::ProcessEvent(Double_t          cent,
					 const AliVVertex* ip,
					 AliMultiplicity*  mult,
					 TTree*            clusters)
{
  // Figure out which centrality bins to fill 
  Int_t    nAcc = 0;
  TIter    next(fCentBins);
  CentBin* bin = 0;
  TList    toRun;
  while ((bin = static_cast<CentBin*>(next()))) {
    if (!bin->Accept(cent)) continue; // Not in range for this bin
    toRun.Add(bin);
    bin->FillIPz(ip->GetZ());
    nAcc++;
  }
  // If we have no centrality bins  to fill, we return immediately 
  if (nAcc <= 0) return;

  if (!(fRecMode & (kNormal | kInjection))) {
    // In case we do not do normal reconstruction nor injection, we
    // should fill our data histograms from the read (ESD)
    // multiplicity.
    ProcessTracklets(toRun, 0, cent, ip->GetZ(), mult);
  }
  if ((fRecMode & kNormal) && !(fRecMode & kInjection)) {
    // Only do explicit normal reconstruction if we're not
    // injecting.  In case of injection, we will also run the normal
    // mode with the same background estimator object.
    Reconstruct(toRun, kNormal, clusters, cent, ip);
    fStatus->Fill(kReconstructed);
  }
  if (fRecMode & kInjection) {
    Reconstruct(toRun, kInjection, clusters, cent, ip);
    fStatus->Fill(kReconstructed);
    fStatus->Fill(kInjected);
  }
  if (fRecMode & kRotation)  {
    Reconstruct(toRun, kRotation, clusters, cent, ip);
    fStatus->Fill(kRotated);
  }
}    


/*********************************************************************
 *
 * Code for processing simulations 
 */
#ifndef __CINT__
#include <AliStack.h>
#include <AliMCEvent.h>
#include <AliGenEventHeader.h>
#include <TParticle.h>
#include <TParticlePDG.h>
#include <TDatabasePDG.h>
#include <TGraphAsymmErrors.h>
#else
class TParticlePDG;             // Auto-load
class TGraphAsymmErrors;        // Auto-load 
class TParticle;
class AliStack;
#endif 

//====================================================================
/**
 * A task to analyse SPD tracklet data for the charged particle
 * pseudorapidity density - simulated data case 
 * 
 */
class AliTrackletdNdetaMCTask : public AliTrackletdNdetaTask
{
public:
  enum {
    kPrimaries = 0x8,
    kSecondaries = 0x10,
    kCombinatorics = 0x20,
    kUncorrCombinatorics = 0x40
  };
  
  //==================================================================
  struct CentBin : public AliTrackletdNdetaTask::CentBin
  {
    CentBin()
      : AliTrackletdNdetaTask::CentBin(),
	fPrimSet(0),
	fSecSet(0),
	fCombSet(0),
	fCombUSet(0),
	fEtaVsIPzMC(0),
	fEtaVsIPzMCSel(0),
	fPdgMC(0),
	fPrimaryPdg(0),
	fSecondaryPdg(0),
	fPrimaryParentPdg(0),
	fSecondaryParentPdg(0),
	fIPzVsMC(0),
	fIPzGen(0),
	fIPzSel(0)
    {}
    CentBin(Float_t c1, Float_t c2, UShort_t recFlags);
    /** 
     * Initialize this object on worker - e.g., called at start of
     * (slave) job
     * 
     * @param parent  Parent container 
     * @param etaAxis Eta axis used 
     * @param ipzAxis IPz axis used 
     */
    virtual void WorkerInit(Container&       parent,
			    const TAxis&     etaAxis,
			    const TAxis&     ipzAxis,
			    const TAxis&     deltaAxis,
			    const TAxis&     dThetaAxis,
			    const TAxis&     dPhiAxis);
    /** 
     * Fill information into histogram set. This is delegated so that
     * sub-classes my override to do more (e.g., for MC).
     * 
     * @param set             Histogram set to fill into 
     * @param isSignal        True if signal survived cuts 
     * @param ipZ             Current event IP's z coordinate 
     * @param eta             Pseudorapidity 
     * @param dPhi            Opening angle in azimuth
     * @param dPhiS           Shifted opening angle in azimuth
     * @param dTheta          Opening polar angle
     * @param dThetaX         Scaled opening polar angle 
     * @param delta           Tracklet @f$\chi^2@f$ 
     * @param weight          Tracklet weight (for MC) 
     * 
     * @return true on success 
     */
    virtual Bool_t FillSet(UShort_t                set,
			   Bool_t                  isSignal,
			   Double_t                ipZ,
			   Double_t                eta,
			   Double_t                dPhi,
			   Double_t                dPhiS,
			   Double_t                dTheta,
			   Double_t                dThetaX,
			   Double_t                delta,
			   Double_t                weight);
    /** 
     * Fill specie information with delta.  Later we will integrate
     * over these to get the before/after histograms.
     * 
     * @param primary    Whether this is a primary 
     * @param pdg        PDG of this particle 
     * @param parentPdg  PDG of parent particle 
     * @param delta      Delta of this particle 
     */
    virtual void FillSpecie(Bool_t   primary,
			    Int_t    pdg,
			    Int_t    parentPdg,
			    Double_t delta,
			    Double_t weight);
    /** 
     * Fill primary interaction point resolution 
     * 
     * @param ipZ     Reconstructed IPz (large negative if invalid) 
     * @param genIPz  Generated IPz 
     */
    virtual void FillMCIPz(Double_t ipZ, Double_t genIPz);
    /** 
     * Fill primary information
     * 
     * @param ipZ      Reconstructed IPz (large negative if invalid) 
     * @param genIPz   Generated IPz 				    
     * @param eta      Pseudorapidity of primary 
     * @param charged  True if primary is charged 
     * @param pdg      Particle type (absolute value)
     * @param weight   Particle weight 
     */
    virtual void FillPrimary(Double_t ipZ, Double_t genIPz,
			     Double_t eta, Bool_t   charged,
			     Int_t    pdg, Double_t weight);
    /** 
     * Calculate alpha 
     * 
     * @param name      Name of generated histogram
     * @param title     Title of generated histogram 
     * @param result    Container to add result to 
     * @param primaries Histogram of primaries 
     * @param data      Histogram of reconstructed tracklets 
     * 
     * @return Histogram containing alpha 
     */
    virtual TH2* MakeAlpha(const char* name,
			   const char* title,
			   Container*  result,
			   TH2*        primaries) const;
    /** 
     * Create a mask with fiducial cuts 
     * 
     * @param a1 One alpha
     * @param a2  Another alpha 
     * @param min Least value to consider 
     * @param max Largest value to consider 
     * 
     * @return Mask 
     */
    virtual TH2* MakeAlphaMask(TH2* a1, TH2* a2,
			       Double_t min=0, Double_t max=2.5) const;
    /** 
     * Make projection of PDG histogram of @f$\Delta@f$ bins within cut. 
     * 
     * @param name   name of projection 
     * @param title  title of projection 
     * @param orig   Original PDG vs @f$\Delta@f$ 
     * @param cut    Upper cut on @f$\Delta@f$ - if negative, project all. 
     * 
     * @return Newly created projection 
     */
    virtual TH1* MakePdgProj(const char* name,
			     const char* title,
			     TH2*        orig,
			     Double_t    cut=-1) const;
    /**
     * Called on master on merged result. 
     * 
     * @param parent    Output container for results 
     * @param ipz       Histogram of ipZ for normalization 
     * @param deltaTail Lower cut on @f$\Delta@f$ tail 
     *
     * @return Result container 
     */
    virtual Container* MasterFinalize(Container* parent,
				      TH1*       ipz,
				      Double_t   deltaCut,
				      Double_t   deltaTail);
  protected:
    /** 
     * Copy constructor 
     * 
     * @param o Object to copy from 
     */
    CentBin(const CentBin& o)
      : AliTrackletdNdetaTask::CentBin(o),
	fPrimSet(0),
	fSecSet(0),
	fCombSet(0),
	fCombUSet(0),
	fEtaVsIPzMC(0),
	fEtaVsIPzMCSel(0),
	fPdgMC(0),
	fPrimaryPdg(0),
	fSecondaryPdg(0),
	fPrimaryParentPdg(0),
	fSecondaryParentPdg(0),
	fIPzVsMC(0),
	fIPzGen(0),
	fIPzSel(0)
    {}
    /** 
     * Assignment operator 
     * 
     * 
     * @return Reference to this object 
     */
    CentBin& operator=(const CentBin&) { return *this; }

    HistoSet* fPrimSet;
    HistoSet* fSecSet;
    HistoSet* fCombSet;
    HistoSet* fCombUSet;
    TH2*      fEtaVsIPzMC;
    TH2*      fEtaVsIPzMCSel;
    TH2*      fIPzVsMC;
    TH1*      fIPzGen;
    TH1*      fIPzSel;
    TH1*      fPdgMC;
    TH2*      fPrimaryPdg;
    TH2*      fSecondaryPdg;
    TH2*      fPrimaryParentPdg;
    TH2*      fSecondaryParentPdg;
    ClassDef(CentBin,1); // Centrality bin 
  };
  //==================================================================
  /** 
   * Constructor 
   */
  AliTrackletdNdetaMCTask()
    : AliTrackletdNdetaTask(),
      fMCStatus(0),
      fReweigh(false)
  {}
  /** 
   * Constructor 
   */
  AliTrackletdNdetaMCTask(const char* name, Bool_t reweigh=false)
    : AliTrackletdNdetaTask(name),
      fMCStatus(0),
      fReweigh(reweigh)
  {}
  /** 
   * Make a centrality bin 
   * 
   * @param c1 Low edge 
   * @param c2 High edge 
   * 
   * @return 
   */
  virtual AliTrackletdNdetaTask::CentBin*
  MakeCentBin(Float_t c1, Float_t c2, UShort_t recFlags)
  {
    return new AliTrackletdNdetaMCTask::CentBin(c1, c2, recFlags);
  }
protected:

  /** 
   * @{ 
   * @name 
   */
  /** 
   * Constructor 
   */
  AliTrackletdNdetaMCTask(const AliTrackletdNdetaMCTask& o)
    : AliTrackletdNdetaTask(o),
      fMCStatus(0),
      fReweigh(false)
  {}
  /** 
   * Assingment operator 
   *
   * @return Reference to this object 
   */
  AliTrackletdNdetaMCTask& operator=(const AliTrackletdNdetaMCTask& o)
  {
    return *this; 
  }
  /** 
   * @{
   * @name Other functions 
   */
  virtual void Print(Option_t* option) const
  {
    AliTrackletdNdetaTask::Print(option);
    
    Printf(" %22s: %s", "Reweigh",   fReweigh ? "yes" : "no");
  }
  /* @} */
  /** 
   * @{ 
   * @name Interface 
   */
  enum {
    kAllMC = 1,
    kMCEvent,
    kStack,
    kGenHeader,
    kDataOK,
    kMCTruth,
    kReweighed,
    kMCCompleted
  };
  /** 
   * Get the CDB reference URL 
   * 
   * @return A fixed string pointing to 2008
   */
  virtual const char* GetCDBReferenceURL() const
  {
    return "alien://Folder=/alice/simulation/2008/v4-15-Release/Residual";
  }
  /** 
   * Initialize histograms.  Called on worker 
   * 
   */
  virtual void WorkerInit();
  /** 
   * Check the event.  Side effect is to set the internal pointer to
   * the MC stack.
   * 
   * @param cent     On return, the centrality 
   * @param ip       On return, the IP coordiantes
   * @param mult     On return, the container of tracklets 
   * @param clusters On return, the clusters - if requested 
   * 
   * @return true if the event is selected 
   */
  virtual Bool_t CheckEvent(Double_t&          cent,
			    const AliVVertex*& ip,
			    AliMultiplicity*&  mult,
			    TTree*&            clusters);
  /** 
   * Called on master on merged result. 
   * 
   * @param results Output container for results 
   */
  virtual void MasterFinalize(Container* results);
  /** 
   * Get weight of a tracklet (for MC)
   * 
   * @param trackletNumber Tracklet number
   * @param mult           Tracklet container
   * 
   * @return The weight 
   */
  virtual Double_t TrackletWeight(Int_t                  trackletNumber,
				  const AliMultiplicity* mult) const;
  /** 
   * Reweigh the input stack.  We do this by assigning weights to each
   * particle according to the primary mother particle of each
   * particle.
   * 
   * @param cent Event centrality 
   */
  void ReweighStack(Double_t cent) const;
  /** 
   * Look-up weight of a particle (secondary or primary).  This is
   * done by recursively calling this function until we hit the
   * primary mother particle. 
   * 
   * @param trackNo   Track number 
   * @param particle  Particle information 
   * @param cent      Event centrality 
   * 
   * @return Weight of particle 
   */
  Double_t LookupWeight(Int_t      trackNo,
			TParticle* particle,
			Double_t   cent) const;
  /** 
   * Function to return the weight of a primary particle.  Here, it
   * always returns 1, but it can be overloaded to do more.
   * 
   * @param particle Primary particle to look-up weight for 
   * @param cent     Event centrality 
   *
   * @return The weight 
   * @todo Implement the actual look-up  
   */
  virtual Double_t LookupWeight(TParticle* particle, Double_t cent) const
  {
    return 1;
  }
  /** 
   * Fill primary information 
   * 
   * @param cent    Centrality 
   * @param ipZ     Reconstructed IP Z coordinate 
   * @param genIPz  Generated IP Z coordinate 
   * @param dataOK  True if event is selected 
   */
  virtual Bool_t FillPrimaries(Double_t cent,
			       Double_t ipZ,
			       Double_t genIPz,
			       Bool_t   dataOK);
  /** 
   * Check if two particle labels eventually leads back to the same
   * mother particle.
   * 
   * @param labels0 First label 
   * @param labels1 Second label
   * 
   * @return true if both labels lead back to the same mother particle 
   */
  Bool_t HaveCommonParent(const Float_t* labels0,
			  const Float_t* labels1) const;
  /** 
   * Fill selected centrality bins with tracklet information 
   * 
   * @param toRun           Selected centrality bins 
   * @param mode            What we're doing 
   * @param mult            Container of tracklets 
   * @param trackletNumber  Tracklet number
   * @param isSignal        If this is considered a signal 
   * @param ipz             Reconstructed IP z coordinate 
   * @param eta             Pseudorapidity of tracklet 
   * @param dPhi            Azimuthal opening angle  
   * @param dPhiS           Shifted azimuthal opening angle 
   * @param dTheta          Polar opening angle 
   * @param dThetaX         Scaled polar opening angle 
   * @param delta           Tracklet @f$\chi^2@f$ 
   * @param weight          Tracklet weight (for MC)
   */
  virtual void FillBins(TList&           toRun,
			AliITSMultRecBg* reco,
			AliMultiplicity* mult,
			Int_t            trackletNumber,
			Bool_t           isSignal,
			Double_t         ipz,
			Double_t         eta,
			Double_t         dPhi,
			Double_t         dPhiS,
			Double_t         dTheta,
			Double_t         dThetaX,
			Double_t         delta,
			Double_t         weight) const;
  /* @} */
  /** 
   * @{ 
   * @name 
   */
  /** 
   * Return an array of interesting PDG codes. 
   * 
   * @param size On return, the size of the array 
   * 
   * @return Pointer to static array 
   */
  static Int_t* PdgArray(Int_t& size);
  /** 
   * Return pointer to a sorted array of interesting PDG codes 
   * 
   * @param size On return, the size of the array 
   * 
   * @return Pointer to static, sorted array of PDGs 
   */
  static Int_t* SortedPdgArray(Int_t& size);
  /** 
   * Find the index of a PDG in the sorted array of PDGs.  Note, if
   * the PDG is 1B or larger, it will be considered valid.  If the
   * exact PDG isn't found, then the size of the array will be
   * returned.
   * 
   * @param pdg PDG code to look up
   * 
   * @return Index of PDG code, or the size of the array 
   */
  static Int_t FindPdg(Int_t pdg);
  /** 
   * Make an axis of PDG codes.  Note, the axis is sorted by PDG
   * number - not abundance.
   * 
   * @return Pointer to static axis object 
   */
  static TAxis* MakePdgAxis();
  /* @} */
  TH1* fMCStatus; // Status of task
  Bool_t fReweigh;
  ClassDef(AliTrackletdNdetaMCTask,1); // MC class for dN/deta u/tracklets
};

//====================================================================
AliTrackletdNdetaMCTask::CentBin::CentBin(Float_t c1, Float_t c2,
					  UShort_t recFlags)
  : AliTrackletdNdetaTask::CentBin(c1,c2,recFlags, +2, 4),
    fPrimSet(0),
    fSecSet(0),
    fCombSet(0),
    fCombUSet(0),
    fEtaVsIPzMC(0),
    fEtaVsIPzMCSel(0),
    fPdgMC(0),
    fPrimaryPdg(0),
    fSecondaryPdg(0),
    fPrimaryParentPdg(0),
    fSecondaryParentPdg(0),
    fIPzVsMC(0),
    fIPzGen(0),
    fIPzSel(0)
{
  fPrimSet  = new HistoSet("primaries", kCyan+2, 34);
  fSecSet   = new HistoSet("secondaries", kMagenta+2, 28);
  fCombSet  = new HistoSet("combinatorics", kYellow+2, 30);
  fCombUSet = new HistoSet("uncorrelatedCombinatorics", kPink+2, 27);
  fHistoSets->Add(fPrimSet);
  fHistoSets->Add(fSecSet);
  fHistoSets->Add(fCombSet);
  fHistoSets->Add(fCombUSet);
}
//____________________________________________________________________
void AliTrackletdNdetaMCTask::CentBin::WorkerInit(Container&       parent,
						  const TAxis&     etaAxis,
						  const TAxis&     ipzAxis,
						  const TAxis&     deltaAxis,
						  const TAxis&     dThetaAxis,
						  const TAxis&     dPhiAxis)
{
  AliTrackletdNdetaTask::CentBin::WorkerInit(parent, etaAxis,
					     ipzAxis, deltaAxis,
					     dThetaAxis, dPhiAxis);
  Container* c        = fContainer;
  TAxis& pdg          = *MakePdgAxis();
  fEtaVsIPzMC         = Make2D(*c, "etaVsIPzMC",
			       "#eta vs IP_{#it{z}} signal",
			       kCyan+2,24, etaAxis, ipzAxis);
  fEtaVsIPzMC->SetYTitle("IP_{#it{z},gen}");
  fEtaVsIPzMCSel      = Make2D(*c, "etaVsIPzMCSel",
			       "#eta vs IP_{#it{z}} selected signal",
			       kMagenta+2,24, etaAxis, ipzAxis);
  fEtaVsIPzMCSel->SetYTitle("IP_{#it{z},rec}");
  fPdgMC              = Make1D(*c, "pdgMC", "Particle types",
			       kRed+2,20, pdg);
  fPdgMC->GetXaxis()->LabelsOption("v");
  fPrimaryPdg         = Make2D(*c, "primaryPdg", "Primary specie",
				   kBlue+2, 20, pdg, deltaAxis);
  fSecondaryPdg       = Make2D(*c, "secondaryPdg", "Secondary specie",
			       kGreen+2, 21, pdg, deltaAxis);
  fPrimaryParentPdg   = Make2D(*c, "primaryParentPdg",
			       "Parent of primary specie",
			       kBlue+2, 20, pdg, deltaAxis);
  fSecondaryParentPdg = Make2D(*c, "secondaryParentPdg",
			       "Parent of secondary specie",
			       kGreen+2, 21, pdg, deltaAxis);	  
  fPrimaryPdg->GetXaxis()->LabelsOption("v");
  fSecondaryPdg->GetXaxis()->LabelsOption("v");
  fPrimaryParentPdg->GetXaxis()->LabelsOption("v");
  fSecondaryParentPdg->GetXaxis()->LabelsOption("v");
  TAxis dIPz(ipzAxis);
  ScaleAxis(dIPz, .1);
  dIPz.SetTitle("IP_{#it{z}} - IP_{#it{z},gen}");      
  fIPzVsMC            = Make2D(*c,"ipzVsGenIPz", "Resolution vs IP_{#it{z}}",
			       kCyan+2, 24, ipzAxis, dIPz);
  fIPzGen             = Make1D(*c,"ipzGen", "Generated IP_{z}",
			       kCyan+2, 25, ipzAxis);
  fIPzSel             = Make1D(*c,"ipzSel", "Selected, generated IP_{#it{z}}",
			       kCyan+2, 21, ipzAxis);
  fIPzGen->SetMarkerSize(1.2*fIPzGen->GetMarkerSize());
  fIPz->SetMarkerStyle(24);
  fIPz->SetMarkerColor(kCyan+2);
  fIPz->SetMarkerSize(1.4*fIPz->GetMarkerSize());
  fIPz->SetLineColor(kCyan+2);
}
//____________________________________________________________________
Bool_t AliTrackletdNdetaMCTask::CentBin::FillSet(UShort_t set,
						 Bool_t   isSignal,
						 Double_t ipZ,
						 Double_t eta,
						 Double_t dPhi,
						 Double_t dPhiS,
						 Double_t dTheta,
						 Double_t dThetaX,
						 Double_t delta,
						 Double_t weight)
{
  if (AliTrackletdNdetaTask::CentBin::FillSet(set, isSignal, ipZ, eta,
					      dPhi, dPhiS, dTheta, dThetaX,
					      delta, weight)) return true;
  HistoSet* hset = 0;
  if      (set == kPrimaries)            hset = fPrimSet;
  else if (set == kSecondaries)          hset = fSecSet;
  else if (set == kCombinatorics)        hset = fCombSet;
  else if (set == kUncorrCombinatorics)  hset = fCombUSet;

  if (!hset) return false;
      
  hset->Fill(isSignal ? ipZ : -1000,
	     eta, dPhi, dPhiS, dTheta, dThetaX, delta, weight);

  return true;
}
//____________________________________________________________________
void AliTrackletdNdetaMCTask::CentBin::FillSpecie(Bool_t   primary,
						  Int_t    pdg,
						  Int_t    parentPdg,
						  Double_t delta,
						  Double_t weight)
{
  Int_t bin       = FindPdg(pdg);
  Int_t parentBin = FindPdg(parentPdg);
  TH2*  pdgHist   = (primary ? fPrimaryPdg       : fSecondaryPdg);
  TH2*  parHist   = (primary ? fPrimaryParentPdg : fSecondaryParentPdg);
  if (pdg >= 0)       pdgHist->Fill(bin,       delta, weight);
  if (parentPdg >= 0) parHist->Fill(parentBin, delta, weight);
}
//____________________________________________________________________
void AliTrackletdNdetaMCTask::CentBin::FillMCIPz(Double_t ipZ, Double_t genIPz)
{
  fIPzGen->Fill(genIPz);
  if (ipZ < -999) return;
  fIPzSel->Fill(genIPz);
  fIPzVsMC->Fill(ipZ, ipZ-genIPz);
}
//____________________________________________________________________
void AliTrackletdNdetaMCTask::CentBin::FillPrimary(Double_t ipZ,
						   Double_t genIPz,
						   Double_t eta,
						   Bool_t   charged,
						   Int_t    pdg,
						   Double_t weight)
{
  // Fill PDG histogram here
  fPdgMC->Fill(FindPdg(pdg));
  if (!charged) return;
  fEtaVsIPzMC->Fill(eta, genIPz, weight);
  if (ipZ > -999) fEtaVsIPzMCSel->Fill(eta, ipZ, weight);    
}

//____________________________________________________________________
TH2*
AliTrackletdNdetaMCTask::CentBin::MakeAlpha(const char* name,
					    const char* title,
					    Container*  result,
					    TH2*        primaries) const
{
  if (!primaries) {
    AliWarning("No primaries supplied for alpha calculation");
    return 0;
  }
  if (!result) {
    AliWarning("No data supplied for alpha calculations");
    return 0;
  }
  TH2* data  = GetH2(result, "signalEst");
  TH2* alpha = static_cast<TH2*>(CloneAndAdd(result,primaries));
  alpha->SetName(name);
  alpha->SetTitle(title);
  alpha->Divide(data);
  return alpha;
}

//____________________________________________________________________
TH2*
AliTrackletdNdetaMCTask::CentBin::MakeAlphaMask(TH2*     a1,
						TH2*     a2,
						Double_t min,
						Double_t max) const
{
  if (!a1 || !a2) return 0;
  TH2* mst = (a1 ? a1 : a2);
  TH2* ret = static_cast<TH2*>(mst->Clone("fiducial"));
  ret->SetTitle("Fiducial cut");
  ret->SetDirectory(0);
  ret->SetMarkerColor(kBlack);
  ret->SetMarkerSize(1);
  ret->SetMarkerStyle(1);
  ret->SetFillColor(kBlack);
  ret->SetFillStyle(1001);
  ret->SetLineColor(kBlack);
  ret->Reset();
  
  for (Int_t ix = 1; ix <= mst->GetNbinsX(); ix++) {
    for (Int_t iy = 1; iy <= mst->GetNbinsY(); iy++) {
      Double_t c1 = (a1 ? a1->GetBinContent(ix,iy) : min);
      Double_t c2 = (a2 ? a2->GetBinContent(ix,iy) : min);
      if ((c1 > min && c1 <= max) ||
	  (c2 > min && c2 <= max))
	ret->SetBinContent(ix,iy,1);
    }
  }
  return ret;  
}

//____________________________________________________________________
TH1*
AliTrackletdNdetaMCTask::CentBin::MakePdgProj(const char* name,
					      const char* title,
					      TH2*        orig,
					      Double_t    cut) const
{
  if (!orig) {
    AliWarningF("Nothing to make %s from", name);
    return 0;
  }
  Int_t first = 1;
  Int_t last  = (cut < 0 ?
		 orig->GetYaxis()->GetNbins() :
		 orig->GetYaxis()->FindBin(cut-1e-6));
  Double_t intg = orig->Integral();
  TH1*     proj = orig->ProjectionX(name,first,last,"e");
  if (!proj) {
    AliWarningF("Failed to project %s between 0 and %f",
		orig->GetName(), cut);
    return 0;
  }
  if (cut > 0) {
    proj->SetMarkerStyle(proj->GetMarkerStyle()+4);
    proj->SetMarkerSize(1.2*proj->GetMarkerSize());
  }
  proj->SetFillStyle(0);
  proj->SetFillColor(0);
  proj->SetTitle(title);
  proj->SetDirectory(0);
  // Scale to integral of whole thing.  In this way, we can see how
  // many particles of each kind (primary, secondary) we cut away with
  // our delta cut, or in case of parents, which kind of parents are
  // most likely to result in ignored particles.  This differs a
  // little from how it was done in Robertos script, where the
  // normalisation is done to the total within the cut.  In that case,
  // the plot shows merely the relative abundance of particles before
  // and after applying the cut - not which were cut away.  On the
  // other hand, in Robertos task the histogram filled was PDG vs
  // centrality - not PDG vs Delta - as is done here and what seems to
  // be assumed in Roberto's script.
  if (intg > 0) proj->Scale(1./intg);
  return proj;
}


//____________________________________________________________________
AliTrackletdNdetaTask::Container*
AliTrackletdNdetaMCTask::CentBin::MasterFinalize(Container* parent,
						 TH1*       ipz,
						 Double_t   deltaCut,
						 Double_t   deltaTail)
{
  Container* result  =
    AliTrackletdNdetaTask::CentBin::MasterFinalize(parent, ipz,
						   deltaCut, deltaTail);
  TH1* ipzRec    = GetH1(fContainer, "ipz"); // Original - not scaled
  TH1* ipzGen    = GetH1(fContainer, "ipzGen");

  fEtaVsIPzMC    = ScaleToIPz(GetH2(fContainer,"etaVsIPzMC"),   ipzGen);
  fEtaVsIPzMCSel = ScaleToIPz(GetH2(fContainer,"etaVsIPzMCSel"),ipzRec);
  fIPzVsMC       = GetH2(fContainer, "ipzVsGenIPz");

  TProfile* ipzVsMC =  fIPzVsMC->ProfileX(fIPzVsMC->GetName());
  ipzVsMC->SetYTitle(Form("#LT%s#GT",
			  fIPzVsMC->GetYaxis()->GetTitle()));
  result->Add(ipzVsMC);
  result->Add(fEtaVsIPzMC);
  result->Add(fEtaVsIPzMCSel);

  TH1* dNdetaAll = AverageOverIPz(fEtaVsIPzMC,    "truedNdeta",    1, ipzGen);
  TH1* dNdetaSel = AverageOverIPz(fEtaVsIPzMCSel, "truedNdetaSel", 1, ipzRec);
  dNdetaAll->SetMarkerSize(1.4*dNdetaAll->GetMarkerSize());
  dNdetaAll->SetTitle("All events");
  dNdetaSel->SetTitle("Selected events");
  dNdetaAll->SetYTitle("d#it{N}_{ch}/d#eta");
  dNdetaSel->SetYTitle("d#it{N}_{ch}/d#eta");
  result->Add(dNdetaAll);
  result->Add(dNdetaSel);
  
  result->Add(MakePdgProj("allPrimaryPdg",  "All primary particle PDGs",
			  GetH2(fContainer, "primaryPdg")));
  result->Add(MakePdgProj("allSecondaryPdg","All secondary particle PDGs",
			  GetH2(fContainer, "secondaryPdg")));
  result->Add(MakePdgProj("cutPrimaryPdg",  "Selected primary particle PDGs",
			  GetH2(fContainer, "primaryPdg"),deltaCut));
  result->Add(MakePdgProj("cutSecondaryPdg","Selected secondary particle PDGs",
			  GetH2(fContainer, "secondaryPdg"),deltaCut));
  result->Add(MakePdgProj("allPrimaryParentPdg",
			  "All primary particle parents PDGs",
			  GetH2(fContainer, "primaryParentPdg")));
  result->Add(MakePdgProj("allSecondaryParentPdg",
			  "All secondary particle parents PDGs",
			  GetH2(fContainer, "secondaryParentPdg")));
  result->Add(MakePdgProj("cutPrimaryParentPdg",
			  "Selected primary particle parents PDGs",
			  GetH2(fContainer, "primaryParentPdg"),deltaCut));
  result->Add(MakePdgProj("cutSecondaryParentPdg",
			  "Selected secondary particle parents PDGs",
			  GetH2(fContainer, "secondaryParentPdg"),deltaCut));
    

  
  fIPzGen = static_cast<TH1*>(CloneAndAdd(result,ipzGen));
  fIPzSel = static_cast<TH1*>(CloneAndAdd(result,GetH1(fContainer, "ipzSel")));
  TGraphAsymmErrors* ipzEff = new TGraphAsymmErrors;
  ipzEff->SetName("ipzEff");
  ipzEff->SetTitle("IP_{z} efficiency");
  ipzEff->BayesDivide(fIPzSel, fIPzGen);
  ipzEff->SetMarkerStyle(20);
  ipzEff->SetMarkerColor(kRed+2);
  ipzEff->SetLineColor(kRed+2);
  ipzEff->SetFillColor(kRed-2);
  ipzEff->SetFillStyle(1001);
  fIPzGen->Scale(1./fIPzGen->GetEntries());
  fIPzSel->Scale(1./fIPzGen->GetEntries());
  result->Add(ipzEff);

  const char* subs[] = { "injection", 
			 "combinatorics",
			 "uncorrelatedCombinatorics",
			 0 };
  const char** ptr = subs;
  while (*ptr) {
    Container* sub = GetC(result, *ptr);
    if (!sub) {
      AliWarningF("Container %s not found in %s", *ptr, result->GetName());
      ptr++;
    }
    TH2* a1 = MakeAlpha("alpha", "#alpha all MC events",
			sub, fEtaVsIPzMC);
    TH2* a2 = MakeAlpha("alphaSel","#alpha selected MC events",sub,
			fEtaVsIPzMCSel);
    TH2* am = MakeAlphaMask(a1, a2);
    sub->Add(am);
    ptr++;
  }
  if (lDebug > 5) {
    result->ls();
  }
  return result;
}

//====================================================================
void AliTrackletdNdetaMCTask::WorkerInit()
{
  AliTrackletdNdetaTask::WorkerInit();
    
  fMCStatus = new TH1F("statusMC", "Status of task (MC part)",
		       kMCCompleted, .5, kMCCompleted+.5);
  fMCStatus->SetMarkerColor(kCyan+2);
  fMCStatus->SetMarkerSize(2);
  fMCStatus->SetLineColor(kCyan+2);
  fMCStatus->SetFillColor(kCyan+2);
  fMCStatus->SetFillStyle(1001);
  fMCStatus->SetDirectory(0);
  fMCStatus->SetStats(0);
  fMCStatus->SetXTitle("Event have");
  fMCStatus->SetYTitle("# Events");
  fMCStatus->GetXaxis()->SetBinLabel(kAllMC, "Been seen");
  fMCStatus->GetXaxis()->SetBinLabel(kMCEvent, "MC Event");
  fMCStatus->GetXaxis()->SetBinLabel(kStack, "Kinematics stack");
  fMCStatus->GetXaxis()->SetBinLabel(kGenHeader, "Generator header");
  fMCStatus->GetXaxis()->SetBinLabel(kDataOK, "Been selected");
  fMCStatus->GetXaxis()->SetBinLabel(kMCTruth, "MC 'Truth'");
  fMCStatus->GetXaxis()->SetBinLabel(kReweighed, "Been reweighed");
  fMCStatus->GetXaxis()->SetBinLabel(kMCCompleted, "Completed for MC");
  fContainer->Add(fMCStatus);

  fStatus->SetFillColor(kCyan+2);
  fStatus->SetLineColor(kCyan+2);
  fStatus->SetFillStyle(0);
  fStatus->SetBarOffset(0.5);

  fIPz   ->SetMarkerSize(1.2*fIPz->GetMarkerSize());
  fIPz   ->SetMarkerStyle(24);
  fIPz   ->SetMarkerColor(kCyan+2);
  fIPz   ->SetLineColor(kCyan+2);

  fCent  ->SetMarkerSize(1.2*fIPz->GetMarkerSize());
  fCent  ->SetMarkerStyle(24);
  fCent  ->SetMarkerColor(kCyan+2);
  fCent  ->SetLineColor(kCyan+2);
}
//____________________________________________________________________
Bool_t AliTrackletdNdetaMCTask::CheckEvent(Double_t&          cent,
					   const AliVVertex*& ip,
					   AliMultiplicity*&  mult,
					   TTree*&            clusters)
{
  // Count all seen 
  fMCStatus->Fill(kAllMC);

  // Get the MC event 
  AliMCEvent* mcEvent = MCEvent();
  if (!mcEvent) {
    AliWarning("No MC event found");
    return false;
  }
  fMCStatus->Fill(kMCEvent);
    
  // Get the stack 
  if (!mcEvent->Stack()) {
    AliWarning("No particle stack in MC event");
    return false;
  }
  fMCStatus->Fill(kStack);
    
  // Get the header and generator IP
  TArrayF            genIP(3);
  AliGenEventHeader* genHeader = mcEvent->GenEventHeader();
  if (!genHeader) {
    AliWarning("No generator header found in MC event");
  }
  else {
    genHeader->PrimaryVertex(genIP);
    fMCStatus->Fill(kGenHeader);
  }

  // Now do the same event selection as for "real" data. 
  Bool_t dataOK = AliTrackletdNdetaTask::CheckEvent(cent,ip,mult,clusters);
  if (dataOK) fMCStatus->Fill(kDataOK);
    
  // We reweigh the MC particles
  if (fReweigh) {
    ReweighStack(cent);
    fMCStatus->Fill(kReweighed);
  }
    
  // Fill information on primaries 
  if (FillPrimaries(cent, ip ? ip->GetZ() : -1000, genIP[2], dataOK))
    fMCStatus->Fill(kMCTruth);

  fMCStatus->Fill(kMCCompleted);
  return dataOK;
}
//____________________________________________________________________
void AliTrackletdNdetaMCTask::MasterFinalize(Container* results)
{
  AliTrackletdNdetaTask::MasterFinalize(results);
  // Make copies of histograms and store 
  fMCStatus =
    static_cast<TH1*>(CloneAndAdd(results,GetH1(fContainer,"statusMC")));
  fMCStatus->Scale(1./fMCStatus->GetBinContent(1));
}

//____________________________________________________________________
Double_t
AliTrackletdNdetaMCTask::TrackletWeight(Int_t        trackletNumber,
					const AliMultiplicity* mult) const
{
  if (!fReweigh) return 1;

  Double_t weight  = 1;
  Int_t    label[] = { mult->GetLabel(trackletNumber, 0),
		       mult->GetLabel(trackletNumber, 1) };
  Int_t    nLayer   = 2;
  if (label[0] == label[1]) nLayer = 1;  // from the same track
  for (Int_t layer = 0; layer < nLayer; layer++) {
    if (label[layer] < 0)
      // If invalid label, give this tracklet 0 weight
      return 0;
    TParticle* particle = MCEvent()->Stack()->Particle(label[layer]);
    if (!particle)
      // If particle isn't found, give this tracklet 0 weight
      return 0;
    // Assume factorization of weights 
    weight *= particle->GetWeight();
  }
    
  return weight;
}
//____________________________________________________________________
void AliTrackletdNdetaMCTask::ReweighStack(Double_t cent) const
{
  if (!fReweigh) return;
    
  AliStack* stack   = MCEvent()->Stack();
  Int_t     nTracks = stack->GetNtrack();
  for (Int_t trackNo = nTracks; trackNo--;) {
    TParticle* particle = stack->Particle(trackNo);
    if (!particle) continue;
    particle->SetWeight(LookupWeight(trackNo, particle, cent));
  }
}
//____________________________________________________________________
Double_t AliTrackletdNdetaMCTask::LookupWeight(Int_t      trackNo,
					       TParticle* particle,
					       Double_t   cent) const
{
  if (!particle)
    // If we do not have a particle, give this tracklet 0 weight!
    return 0;
  AliStack* stack   = MCEvent()->Stack();
  if (!stack->IsPhysicalPrimary(trackNo)) {
    // If this is not a primary, recursively call this member
    // function until we do find the primary mother.
    Int_t mother = particle->GetFirstMother();
    return LookupWeight(mother, stack->Particle(mother), cent);
  }
  // Now, we can do a look-up on the primary particle properties to
  // get the associated weight.
  return LookupWeight(particle, cent);
}

//____________________________________________________________________
Bool_t AliTrackletdNdetaMCTask::FillPrimaries(Double_t cent,
					      Double_t ipZ,
					      Double_t genIPz,
					      Bool_t   dataOK)
{
  AliInfoF("Fill primaries for ipZ=%f, genIPz=%f centrality=%5.1f%% (%s)",
	   ipZ, genIPz, cent, (dataOK ? "selected" : "ignored"));
  Int_t    nAcc = 0;
  TIter    next(fCentBins);
  CentBin* bin = 0;
  TList    toRun;
  while ((bin = static_cast<CentBin*>(next()))) {
    if (!bin->Accept(cent)) continue; // Not in range for this bin
    toRun.Add(bin);
    bin->FillMCIPz(dataOK ? ipZ : -1000, genIPz);
    nAcc++;
  }
  // If there's no bins to fill, give up now 
  if (nAcc <= 0) {
    AliWarningF("No centrality bins found for %f%%", cent);
    return false;
  }

  AliStack* stack   = MCEvent()->Stack();
  Int_t     nTracks = stack->GetNtrack();
  for (Int_t trackNo = nTracks; trackNo--; ) {
    if (!stack->IsPhysicalPrimary(trackNo)) {
      // AliWarningF("Track # %6d not a primary", trackNo);
      // Not a primary, go on
      continue;
    }
    // Get the particle 
    TParticle* particle = stack->Particle(trackNo);
    if (!particle) {
      AliWarningF("No particle found for track # %d", trackNo);
      continue;
    }
    // Get theta
    Double_t theta = particle->Theta();
    // Check for beam-like particle 
    if (theta < 1e-6 || TMath::Abs(theta-TMath::Pi()) < 1e-6) {
      AliWarningF("Track # %6d is beam-like (%f)", trackNo,
		  TMath::RadToDeg()*theta);    
      continue;
    }
    // Get pseudorapidity, transverse momentum, weight, and type
    Double_t      eta    = particle->Eta();
    Double_t      pT     = particle->Pt();
    Double_t      weight = particle->GetWeight();
    Int_t         pdgID  = TMath::Abs(particle->GetPdgCode());
    TParticlePDG* pdg    = particle->GetPDG();
    // @todo Fill pT spectra of charged here
      
    TIter nextB(&toRun);
    while ((bin = static_cast<CentBin*>(nextB()))) {
      // AliInfoF("Filling track %6d into %s", trackNo, bin->Name());
      bin->FillPrimary(dataOK ? ipZ : -1000, genIPz,
		       eta, pdg->Charge()!=0,pdgID, weight);
    }
  }
  return true;
}

//____________________________________________________________________
Bool_t AliTrackletdNdetaMCTask::HaveCommonParent(const Float_t* labels0,
						 const Float_t* labels1) const
{
  // Max search depth 
  const Int_t kMaxParents = 50;

  // Arrays and counters of parents
  Int_t parents[2][kMaxParents];
  Int_t nParents[] = { 0, 0 };
    
  // A convinience (double) array
  const Float_t* labels[2] = { labels0, labels1 };

  // Total number of tracks
  Int_t nTracks = MCEvent()->Stack()->GetNtrack();
    
  // Loop over layers
  for (Int_t layer = 0; layer < 2; layer++) {
    // Loop over labels - never more than 3 
    for (Int_t labelNo = 0; labelNo < 3; labelNo++) {
      Int_t label = Int_t(labels[layer][labelNo]);
      // Check validity
      if (label < 0 || label >= nTracks) continue;

      while (nParents[layer] < kMaxParents) {
	// Set parent particle
	parents[layer][nParents[layer]++] = label;
	// Get the particle
	TParticle* particle = MCEvent()->Stack()->Particle(label);
	// if there's no particle, break out
	if (!particle) break;
	// Look at mother of this particle
	label = particle->GetFirstMother();
	// If there's no mother, break out
	if (label < 1 || label > nTracks) break;
      }
    }
  }
  // Now compare the two arrays of parents to see if we can find a
  // match
  for (Int_t label0No = nParents[0]; label0No--; ) {
    for (Int_t label1No = nParents[1]; label1No--; ) {
      if (parents[0][label0No] == parents[1][label1No])
	// Found common parent
	return true;
    }
  }
  return false;
}

//____________________________________________________________________
void AliTrackletdNdetaMCTask::FillBins(TList&           toRun,
				       AliITSMultRecBg* reco,
				       AliMultiplicity* mult,
				       Int_t            trackletNumber,
				       Bool_t           isSignal,
				       Double_t         ipz,
				       Double_t         eta,
				       Double_t         dPhi,
				       Double_t         dPhiS,
				       Double_t         dTheta,
				       Double_t         dThetaX,
				       Double_t         delta,
				       Double_t         weight) const
{
  UShort_t mode = FindMode(reco);
  if (mode >= 999) return;

  UShort_t which     = 999;
  Bool_t   uncorr    = false;
  Int_t    pdg       = -1;
  Int_t    parentPdg = -1;
  if (mode == kNormal) { // Perhaps need combinatorics from ESD
    // Get MC labels, and derive the type off the labels
    // - Labels are equal, check if primary (0) of secondary (1)
    // - Labels are not equal - combinatorial background (2)
    // Do this when looking at normal reconstruction data 
    Int_t label0 = mult->GetLabel(trackletNumber, 0);
    Int_t label1 = mult->GetLabel(trackletNumber, 1);
    which        = kCombinatorics; // Assume combinatorial background
    if (label0 == label1)
      which = (MCEvent()->Stack()->IsPhysicalPrimary(label0) ?
	       kPrimaries: kSecondaries);

    if (reco && which == kCombinatorics) {
      // Get tracklet parameters from reconstructor
      Float_t* trl = reco->GetTracklet(trackletNumber);
      // Get cluster identifiers
      Int_t cl0    = Int_t(trl[AliITSMultReconstructor::kClID1]);
      Int_t cl1    = Int_t(trl[AliITSMultReconstructor::kClID2]);
      // Get array of labels
      Float_t* labels0 = (reco->GetClusterOfLayer(0,cl0)+
			  AliITSMultRecBg::kClMC0);
      Float_t* labels1 = (reco->GetClusterOfLayer(1,cl1)+
			  AliITSMultRecBg::kClMC0);
      uncorr = !HaveCommonParent(labels0, labels1);
    }
    if (which != kCombinatorics) {
      // Track to mother particle
      AliStack*   stack    = MCEvent()->Stack();
      Int_t       nTracks  = stack->GetNtrack();
      TParticle*  particle = stack->Particle(label0);
      pdg                  = TMath::Abs(particle->GetPdgCode());
      parentPdg            = -1;
      Int_t       parentID = particle->GetFirstMother();
      while (parentID >= 0 && parentID < nTracks) {
	particle  = stack->Particle(parentID);
	parentPdg = TMath::Abs(particle->GetPdgCode());
	parentID  = particle->GetFirstMother();
      }
    }
  }

  // Re-implemented from normal case - except part with the
  // uncorrelated background.
  TIter nextB(&toRun);
  CentBin* bin = 0;
  while ((bin = static_cast<CentBin*>(nextB()))) {
    // Normal operation - as per real data 
    bin->FillSet(mode, isSignal, ipz, eta,
		 dPhi, dPhiS, dTheta, dThetaX, delta, weight);
    if (which >= 999) continue;
      
    // Additional stuff for MC - i.e., fill primary, secondary, or
    // combinatorial distributions
    bin->FillSet(which, isSignal, ipz, eta, dPhi, dPhiS,
		 dTheta, dThetaX, delta, weight);

    if (which != kCombinatorics)
      bin->FillSpecie(which == kPrimaries, pdg, parentPdg, delta, weight);   
    if (!uncorr) continue;
    bin->FillSet(kUncorrCombinatorics, isSignal, ipz, eta, dPhi, dPhiS,
		 dTheta, dThetaX, delta, weight);
  }
}    

//====================================================================
Int_t* AliTrackletdNdetaMCTask::PdgArray(Int_t& size)
{
  static Int_t codes[] = {
    211,        // #pi^{+}			 
    2212, 	  // p				 
    321, 	  // K^{+}			 
    323, 	  // K^{*+}			 
    11, 	  // e^{-}			 
    13, 	  // #mu^{-}			 
    213, 	  // #rho^{+}			 
    411, 	  // D^{+}			 
    413, 	  // D^{*+}			 
    431, 	  // D_{s}^{+}			 
    433, 	  // D_{s}^{*+}			 
    1114, 	  // #Delta^{-}			 
    2214, 	  // #Delta^{+}			 
    2224, 	  // #Delta^{++}			 
    3112, 	  // #Sigma^{-}			 
    3222, 	  // #Sigma^{+}			 
    3114, 	  // #Sigma^{*-}			 
    3224, 	  // #Sigma^{*+}			 
    4214, 	  // #Sigma^{*+}_{c}		 
    4224, 	  // #Sigma^{*++}_{c}		 
    3312, 	  // #Xi^{-}			 
    3314, 	  // #Xi^{*-}			 
    4122, 	  // #Lambda^{+}_{c}		 
    2112, 	  // n				 
    2114, 	  // #Delta^{0}			 
    22, 	  // #gamma			 
    310, 	  // K^{0}_{S}			 
    130, 	  // K^{0}_{L}			 
    311, 	  // K^{0}			 
    313, 	  // K^{*}			 
    221, 	  // #eta				 
    111, 	  // #pi^{0}			 
    113, 	  // #rho^{0}			 
    333, 	  // #varphi			 
    331, 	  // #eta'			 
    223, 	  // #omega			 
    3122, 	  // #Lambda			 
    3212, 	  // #Sigma^{0}			 
    4114, 	  // #Sigma^{*0}_{c}		 
    3214, 	  // #Sigma^{*0}			 
    421, 	  // D^{0}			 
    423, 	  // D^{*0}			 
    3322, 	  // #Xi_{0}			 
    3324, 	  // #Xi^{*0}			 
    4132, 	  // #Xi^{0}_{c}			 
    4314,	  // #Xi^{*0}_{c}
    1000000000  // Nuclei			 
  };
  size = sizeof(codes) / sizeof(Int_t);
  return codes; // Others                         
}

//____________________________________________________________________
Int_t* AliTrackletdNdetaMCTask::SortedPdgArray(Int_t& size)
{
  static Int_t  nCodes = 0;
  static Int_t* codes  = 0;
  static Int_t* sorted = 0;
  if (sorted) {
      size = nCodes;
      return sorted;
  }
  
  codes      = PdgArray(nCodes);
  size       = nCodes;
  sorted     = new Int_t[nCodes];
  Int_t* idx = new Int_t[nCodes];
  TMath::Sort(nCodes, codes, idx, false);
  for (Int_t i = 0; i < nCodes; i++) {
    sorted[i] = codes[idx[i]];
  }
  delete [] idx;
  return sorted;
}

//____________________________________________________________________
Int_t AliTrackletdNdetaMCTask::FindPdg(Int_t pdg)
{
  Int_t  size  = 0;
  Int_t* array = SortedPdgArray(size);
  Int_t  idx   = TMath::BinarySearch(size, array, pdg);
  // Printf("Lookup of PDG %7d -> %d", pdg, idx);
  if (idx == size-1)     return idx;
  if (array[idx] != pdg) return size;
  return idx;
}

//____________________________________________________________________
TAxis* AliTrackletdNdetaMCTask::MakePdgAxis()
{
  static TAxis* axis = 0;
  if (axis) return axis;
  
  Int_t  size  = 0;
  Int_t* array = SortedPdgArray(size);
  // Printf("Got array of size %d", size);
  axis         = new TAxis(size+1, -.5, size+.5);
  // axis->LabelOption("v");
  FixAxis(*axis, "");
  TDatabasePDG* pdgDb = TDatabasePDG::Instance();
  for (Int_t i = 1; i < size; i++) {
    TString       name = "?";
    Int_t         pdg  = array[i-1];
    TParticlePDG* pdgP = pdgDb->GetParticle(pdg);
    if (pdgP)     name = pdgP->GetName();
    // Printf("%3d/%3d: %6d - %s", i-1, size, pdg, name.Data());
    axis->SetBinLabel(i, name.Data());
  }
  axis->SetBinLabel(size,   "Nuclei");
  axis->SetBinLabel(size+1, "Unknown");
  return axis;
}

#endif
/**
 * Post-processing notes - for bg = Comb
 *
 * MC:
 * - Find eta vs Z from MC data 
 *   (mc_TrData_ZvEtaCutT -> mc_RawWithCut)
 * - Clone (mc_RawWithCut -> mc_SignalWithCut)
 * - Clone for lables (mc_RawWithCut -> SignalWithCut_bgMCLabels)
 * - Find eta vs Z from MC combinatorial background processing 
 *   (mc_TrComb_ZvEtaCutT -> mc_BgEst_MCLB)
 * - Make 1-beta from MC comb
 *   (mc_BgEst_MCLB -> mc_1mBeta)
 * - Make 1-beta from MC comb 
 *   (mc_TrComb_ZvEtaCutT -> mc_BgMC)
 *   - Clone (mc_BgMC -> mc_h1mBetaMC)
 *   - Scale mc_h1mBetaMC by mc_RawWithCut
 *   - Clone scaled (mc_BgMC -> mc_h1mBetaMC_scl)
 *   - Calculate 1-beta (mc_h1mBetaMC and mc_h1mBetaMC_scl)
 * - Multiply mc_SignalWithCut_bgMCLabels with mc_h1mBetaMC
 * 
 * - Obtain delta from MC data processing
 *    (mc_TrData_WDist -> mc_DistRawData)
 *   - Calculate integral of tail and normalize 
 *     (mc_DistRawData)
 * - Obtain delta from MC combinatorics processing 
 *     (mc_TrComb_WDist -> mc_DistBgMC)
 *   - Scale mc_DistBgMC by tail of mc_DistRawData
 * - Clone delta from MC combinatorics 
 *   (mc_DistBgMC -> mc_DistRawGenBgNorm)
 *   - Normalize mc_DistRawGenBgNorm relative to mc_DistRawData 
 *
 * - Scale mc_BgEst_MCLB by tail of mc_DistRawData
 * - Multiply mc_SignalWithCut by mc_h1mBetaMC
 * - Calculate mc_1mBeta as ratio of mc_BgEst_MCLB to mc_RawWithCut
 * 
 * Data:
 * - Obtain eta vs vz 
 *   (dt_TrData_ZvEtaCutT -> dt_RawWithCut)
 * - Clone (dt_RawWithCut -> dt_SignalWithCut)
 * - Find eta vs vz from MC comb processing 
 *   (mc_TrComb_ZvEtaCutT -> dt_BgEst_MCLB)
 * - Make 1-beta 
 *   (dt_BgEst_MCLB -> dt_1mBeta)
 * 
 * - Obtain delta from real data processing 
 *   (dt_TrData_WDist -> dt_DistRawData)
 *   - Integrate tail of dt_DistRawData and scale 
 * - Get previously calculated mc_DistBgMC 
 * - Clone delta from MC combinatorics
 *   (mc_DistBgMC -> dt_DistRawGenBgNorm)
 *   - Normalize tail of dt_DistRawGenBgNorm relative to dt_DistRawData 
 *
 * - scale dt_BgEst_MCLB by integral of tail 
 * - Multiply dt_SignalWithCut by mc_h1mBetaMC_scl
 * - Calculate dt_1mBeta 1 - ratio of dt_BgEst_MCLB to dt_RawWithCut
 * 
 * 
 * Both:
 * - Get eta vs vz from MC primaries 
 *   (zvEtaPrimMC -> Alpha)
 *   - Scale by mc_1mBeta and mc_RawWithCut
 * - Get eta vs vz from MC primaries 
 *   (zvEtaPrimMC -> AlphaMC)
 *   - Scale mc_h1mBetaMC and mc_RawWithCut
 * - Get eta vs vz minus comb. 
 *   (dt_SignalWithCut -> SignalEstCorr)
 *   - Multiply by AlphaMC -> DataCorrSignalX
 *   - Scale DataCorrSignalX by bin width
 */
/* 
 * Result calculation in case of use of MC label combinatorical
 * background (Delta distributions ignored).
 *
 *           P           C        P      C
 *    R = ---------- (1-k--)M = ---- (1-k--)M 
 *        (1-C/M*)M*     M*     M*-C     M*
 *
 * where 
 * 
 *  - P  is the primary particle distribution 
 *  - C  is the distribution of tracklets with two mothers 
 *  - M* is the observed distribution in simulated data 
 *  - M  is the observed distribution in real data 
 *  - k  is a scaling constant (1.3?)
 *
 * In case k=1, then we have 
 *
 *         P 
 *     R = -- M
 *         M*
 *
 * and we see that P/M* is an effective correction for secondaries and
 * other algorithmic effects.
 */
/**
 * Post processing notes for bg = Inj 
 *
 * For MC 
 *
 * - Get TrData_ZvEtaCutT -> mc_RawWithCut
 * - Get vertex dist. and scale by integral 
 * - Copy mc_RawWithCut -> mc_SignalWithCut
 * - Copy mc_RawWithCut -> mc_SignalWithCut_bgMCLabels
 * - Get TrInj_ZvEtaCutT -> mc_BgEst
 * - Copy mc_BgEst -> mc_1mBeta
 * - Comb stuff 
 *   - TrComb_ZvEtaCutT -> mc_BgMC
 *   - Copy mc_BgMC -> mc_h1mBetaMC
 *   - Scale mc_h1mBetaMC by mc_RawWithCut
 *   - Copy mc_h1mBetaMC -> mc_h1mBetaMC_scl
 *   - Calculate 1-beta and 1-k*beta 
 *   - Scale mc_SignalWithCut_bgMCLabels by mc_h1mBetaMC
 * - Get TrData_WDist -> mc_DistRawData
 * - Integrate mc_DistRawData oer tail 
 * - Scale mc_DistRawData by tail 
 * - Get TrInj_WDist -> mc_DistRawGenBgNorm
 * - Scale mc_DistRawGenBgNorm by tail of mc_DistRawData
 * - Copy mc_DistRawGenBgNorm -> mc_DistRawGenBgNorm_bgNorm
 * - Scale by mc_DistRawGenBgNorm_bgNorm by ratio of tail 
 *   from mc_DistRawDat and mc_DistRawGenBgNorm and return scale
 * - Scale mc_BgEst by ratios of tail 
 * - if use lables 
 *     - Scale  mc_SignalWithCut by mc_h1mBetaMC
 * - else 
 *     - Subtract mc_BgEst from mc_SignalWithCut
 * - Calculate mc_1mBeta as 1 minus mc_BgEst over mc_RawWithCut
 *
 * Real data 
 *
 * - Get TrData_ZvEtaCutT -> dt_RawWithCut
 * - Get vertex dist. and scale by integral 
 * - Copy dt_RawWithCut -> dt_SignalWithCut
 * - Get TrInj_ZvEtaCutT -> dt_BgEst
 * - Copy dt_BgEst -> dt_1mBeta
 * - Get TrData_WDist -> dt_DistRawData
 * - Integrate dt_DistRawData over tail 
 * - Scale dt_DistRawData by tail 
 * - Get TrInj_WDist -> dt_DistRawGenBgNorm
 * - Scale dt_DistRawGenBgNorm by tail of mc_DistRawData
 * - Copy dt_DistRawGenBgNorm -> dt_DistRawGenBgNorm_bgNorm
 * - Scale by mc_DistRawGenBgNorm_bgNorm by ratio of tail 
 *   from dt_DistRawDat and dt_DistRawGenBgNorm and return scale
 * - Scaling dt_BgEst by ratios of tail 
 * - if labels
 *   - Scale dt_SignalWithCut by mc_h1mBetaMC_scl
 * - else 
 *   - Subtract dt_BgEst from dt_SignalWithCut
 * - Calculate dt_1mBeta as 1 minus dt_BgEst over dt_RawWithCut
 *
 * Both 
 *
 * - Get zvEtaPrimMC -> Alpha
 *   - Scale by mc_1mBeta * mc_RawWithCut
 * - Get zvEtaPrimMC -> AlphaMC
 *   - Scale by mc_h1mBetaMC * mc_RawWithCut
 * - Copy dt_SignalWithCut -> SignalEstCorr
 * - if labels 
 *   - Scale SignalEstCorr by AlphaMC
 * - else 
 *   - Scale SignalEstCorr by Alpha
 * - Copy SignalEstCorr to SignalEstCorrX
 * - Do normal mean over IPz of SignalEstCorrX -> DataCorrSignalX
 * - Scale DataCorrSignalX by bin width
 * 
 * 
 */
// Local Variables:
//  mode: C++
// End:
//
// EOF
// 
