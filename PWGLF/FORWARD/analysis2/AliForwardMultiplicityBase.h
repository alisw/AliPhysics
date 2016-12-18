// 
// Base class for classes that calculate the multiplicity in the
// forward regions event-by-event
// 
#ifndef ALIFORWARDMULTIPLICITYBASE_H
#define ALIFORWARDMULTIPLICITYBASE_H
/**
 * @file   AliForwardMultiplicityBase.h
 * @author Christian Holm Christensen <cholm@dalsgaard.hehi.nbi.dk>
 * @date   Wed Mar 23 14:06:29 2011
 * 
 * @brief  
 * 
 * 
 * @ingroup pwglf_forward_aod
 */
#include "AliBaseESDTask.h"
#include "AliForwardUtil.h"
#include "AliAODForwardMult.h"
#include "AliAODForwardEP.h"
#include "AliAODMultEventClass.h"
// class AliFMDEnergyFitter;
class AliFMDESDFixer;
class AliFMDSharingFilter;
class AliFMDDensityCalculator;
class AliFMDCorrector;
class AliFMDHistCollector;
class AliFMDEventPlaneFinder;
// class AliMultEventClassifier;
class AliAODHandler;
class AliESDEvent;
class TH2D;
class TList;
class TTree;
class TAxis;
class TProfile;

/** 
 * @defgroup pwglf_forward_aod AOD
 * 
 * Code to do with AOD production 
 *
 * @ingroup pwglf_forward_topical
 */
/** 
 * Base class for classes that calculate the multiplicity in the
 * forward regions event-by-event
 * 
 * @par Inputs: 
 *   - AliESDEvent 
 *
 * @par Outputs: 
 *   - AliAODForwardMult 
 * 
 * @par Histograms 
 *   
 * @par Corrections used 
 * 
 * @ingroup pwglf_forward_aod
 * 
 */
class AliForwardMultiplicityBase : public AliBaseESDTask
{
public:
  /** 
   * Bins in timing histogram 
   */
  enum { 
    kTimingEventInspector    = 1,
    kTimingSharingFilter     = 2, 
    kTimingDensityCalculator = 3, 
    kTimingCorrections       = 4, 
    kTimingHistCollector     = 5, 
    kTimingEventPlaneFinder  = 6, 
    kTimingTotal             = 7
  };
  /** 
   * Bins in status histogram 
   */
  enum {
    kStatusNoEvent = 1,
    kStatusNoTrigger,
    kStatusNoSPD,
    kStatusNoFMD,
    kStatusNoVertex,
    kStatusPileup,
    kStatusSPDOutlier,
    kStatusIPzOutOfRange,
    kStatusFailSharing,
    kStatusFailDensity,
    kStatusFailEventPlane,
    kStatusOutlier,
    kStatusFailCorrector,
    kStatusFailCollector,
    kStatusNotAdded,
    kStatusAllThrough       
  };
  /** 
   * @{ 
   * @name Interface methods 
   */
  /** 
   * Create output objects 
   *
   * @return true on success
   */
  virtual Bool_t Book();
  /** 
   * Initialise the sub objects and stuff.  Called on first event
   *
   * @param vertex Vertex axis to use 
   * @param eta    Eta axis to use 
   *
   * @return false on errors 
   */
  virtual Bool_t PreData(const TAxis& vertex, const TAxis& eta);
  /** 
   * Called after processing a single event - should not do anything
   * but clear data, etc.
   * 
   * @return true on success
   */
  virtual Bool_t PostEvent();
  /** 
   * End of job
   * 
   * @return true on success
   */
  virtual Bool_t Finalize();
  /** 
   * @} 
   */
  /** 
   * Print information 
   * 
   * @param option Not used
   */
  virtual void Print(Option_t* option="") const;
  /** 
   * Set whether to make a timing histogram 
   * 
   * @param enable 
   */
  virtual void SetDoTiming(Bool_t enable=true) { fDoTiming = enable; }
  /** 
   * Whether to enable low-flux code 
   * 
   * @param use IF true, enable low-flux code 
   */
  virtual void SetEnableLowFlux(Bool_t use=true) { fEnableLowFlux = use; }
  /** 
   * @{ 
   * @name Access to sub-algorithms 
   */
  /** 
   * Get a reference to the event inspector. User must override this
   * to return proper object
   * 
   * @return Reference to the event inspector 
   */
  // virtual AliMultEventClassifier& GetMultEventClassifier() = 0;
  /** 
   * Get a reference to the event inspector. User must override this
   * to return proper object
   * 
   * @return Reference to the event inspector 
   */
  // virtual const AliMultEventClassifier& GetMultEventClassifier() const = 0;
  /**
   * Get reference to the ESDFixer algorithm 
   * 
   * @return Reference to AliFMDESDFixer object 
   */
  virtual AliFMDESDFixer& GetESDFixer() = 0;
  /**
   * Get reference to the SharingFilter algorithm 
   * 
   * @return Reference to AliFMDSharingFilter object 
   */
  virtual AliFMDSharingFilter& GetSharingFilter() = 0;
  /**
   * Get reference to the DensityCalculator algorithm 
   * 
   * @return Reference to AliFMDDensityCalculator object 
   */
  virtual AliFMDDensityCalculator& GetDensityCalculator() = 0;
  /**
   * Get reference to the Corrections algorithm 
   * 
   * @return Reference to AliFMDCorrector object 
   */
  virtual AliFMDCorrector& GetCorrections() = 0;
  /**
   * Get reference to the HistCollector algorithm 
   * 
   * @return Reference to AliFMDHistCollector object 
   */
  virtual AliFMDHistCollector& GetHistCollector() = 0;
  /**
   * Get reference to the ESDFixer algorithm 
   * 
   * @return Reference to AliFMDESDFixer object 
   */
  virtual const AliFMDESDFixer& GetESDFixer() const = 0;
  /**
   * Get reference to the SharingFilter algorithm 
   * 
   * @return Reference to AliFMDSharingFilter object 
   */
  virtual const AliFMDSharingFilter& GetSharingFilter() const = 0;
  /**
   * Get reference to the DensityCalculator algorithm 
   * 
   * @return Reference to AliFMDDensityCalculator object 
   */
  virtual const AliFMDDensityCalculator& GetDensityCalculator() const = 0;
  /**
   * Get reference to the Corrections algorithm 
   * 
   * @return Reference to AliFMDCorrector object 
   */
  virtual const AliFMDCorrector& GetCorrections() const = 0;
  /**
   * Get reference to the HistCollector algorithm 
   * 
   * @return Reference to AliFMDHistCollector object 
   */
  virtual const AliFMDHistCollector& GetHistCollector() const = 0;
   /**
   * Get reference to the EventPlaneFinder algorithm 
   * 
   * @return Reference to AliFMDEventPlaneFinder object 
   */
  virtual AliFMDEventPlaneFinder& GetEventPlaneFinder() = 0;
  /**
   * Get reference to the EventPlaneFinder algorithm 
   * 
   * @return Reference to AliFMDEventPlaneFinder object 
   */
  virtual const AliFMDEventPlaneFinder& GetEventPlaneFinder() const = 0;
  /* @} */

  /** 
   * Set the debug level 
   * 
   * @param dbg 
   */
  virtual void SetDebug(Int_t dbg);
  /** 
   * Set whether to make separate branches for each ring.  If enabled
   * there will be 5 additional branches on the AOD tree - each
   * holding a TH2D object of the charged particle multiplicity in
   * @f$(\eta,\varphi)@f$ bins for that event.
   * 
   * @param use If true, make separate branches for each ring. 
   */
  void SetStorePerRing(Bool_t use) { fStorePerRing = use; }
  /** 
   * For which triggers to add internally
   * 
   * @param mask Trigger mask as defined in AliAODForwardMult
   */
  void SetAddMask(UInt_t mask)  { fAddMask = mask; }
protected: 
  /** 
   * Constructor 
   * 
   * @param name Name of task 
   */
  AliForwardMultiplicityBase(const char* name); 
  /** 
   * Constructor
   */
  AliForwardMultiplicityBase() 
    : AliBaseESDTask(), 
      fEnableLowFlux(true), 
      fStorePerRing(false),
      fHData(0),
      fHistos(),
      fAODFMD(),
      fAODEP(),
      // fAODRef(),
      fRingSums(),
      fDoTiming(false), 
      fHTiming(0),
      fHStatus(0),
      fAddMask(AliAODForwardMult::kInel)
  {}
  /** 
   * Copy constructor 
   * 
   * @param o Object to copy from 
   */
  AliForwardMultiplicityBase(const AliForwardMultiplicityBase& o);
  /** 
   * Assignment operator 
   * 
   * @param o Object to assign from 
   * 
   * @return Reference to this object 
   */
  AliForwardMultiplicityBase& operator=(const AliForwardMultiplicityBase& o);
  /** 
   * Initialize members based on eta and vertex axis - only available
   * after first event - called from SetupForData.
   * 
   * @param pe @f$\eta@f$ axis
   * @param pv Interaction point Z-coordinate axis 
   */
  virtual void InitMembers(const TAxis& pe, const TAxis& pv);
  /**
   * Create output branches - called from UserCreateOutputObjects
   */
  virtual void CreateBranches(AliAODHandler* ah);
  /** 
   * Do estimates of @f$dN/d\eta@f$  - called at Terminate
   * 
   * @param input  Input list
   * @param output Output list
   */
  virtual void EstimatedNdeta(const TList* input, TList* output) const;
  /** 
   * Calculate a simple dN/deta from all accepted events 
   * 
   * @param input  Input list
   * @param output Output list
   * @param nTr    On return, number of triggers
   * @param nTrVtx On return, number of trigger+vertex events
   * @param nAcc   On return, number of accepted events
   * 
   * @return true on success 
   */
  virtual Bool_t MakeSimpledNdeta(const TList* input, 
				  TList*       output,
				  Double_t&    nTr, 
				  Double_t&    nTrVtx, 
				  Double_t&    nAcc);
  /** 
   * Make Ring @f$ dN/d\eta @f$ histogram and a stack 
   * 
   * @param input      List with summed signals 
   * @param output     Output list 
   * @param inName     Input name 
   * @param outName    Output name
   * @param style      Style 
   */
  virtual void MakeRingdNdeta(const TList* input, 
			      const char*  inName,
			      TList*       output,
			      const char*  outName,
			      Int_t        style=20) const;
  TAxis* DefaultEtaAxis() const { return new TAxis(200,-4,6); }
  TAxis* DefaultVertexAxis() const { return new TAxis(10,-10,10); }
  Bool_t                 fEnableLowFlux;// Whether to use low-flux specific code
  Bool_t                 fStorePerRing; // Store each ring on separate branch
  TH2D*                  fHData;        // Summed 1/Nd^2N_{ch}/dphideta
  AliForwardUtil::Histos fHistos;       // Cache histograms 
  AliAODForwardMult      fAODFMD;       // Output object
  AliAODForwardEP        fAODEP;        // Output object
  // AliAODMultEventClass   fAODRef;       // Reference multiplicity
  AliForwardUtil::Histos fRingSums;     // Cache histograms 
  Bool_t                 fDoTiming;     // Whether to do timing or not
  TProfile*              fHTiming;      // Timing histogram 
  TH1*                   fHStatus;      // Status histogram
  UInt_t                 fAddMask;      // For which triggers to add internally
  ClassDef(AliForwardMultiplicityBase,7) // Forward multiplicity class
};

#endif

// Local Variables:
//  mode: C++
// End:

