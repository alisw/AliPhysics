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
#include <AliAnalysisTaskSE.h>
#include "AliForwardUtil.h"
#include "AliAODForwardMult.h"
#include "AliAODForwardEP.h"
class AliFMDEventInspector;
class AliFMDEnergyFitter;
class AliFMDSharingFilter;
class AliFMDDensityCalculator;
class AliFMDCorrector;
class AliFMDHistCollector;
class AliAODHandler;
class AliForwardCorrectionManager;
class AliFMDEventPlaneFinder;
class AliESDEvent;
class TH2D;
class TList;
class TTree;
class TAxis;

/** 
 * @defgroup pwglf_forward PWGLF Forward analysis
 *
 * Code to do the multiplicity analysis in the forward psuedo-rapidity
 * regions
 *
 */
/** 
 * @defgroup pwglf_forward_tasks Tasks
 *
 * Code to do the multiplicity analysis in the forward psuedo-rapidity
 * regions
 *
 * @ingroup pwglf_forward 
 */
/** 
 * @defgroup pwglf_forward_topical Topical
 *
 * The code divided according to topic
 */
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
 * @ingroup pwglf_forward_tasks
 * @ingroup pwglf_forward_aod
 * 
 */
class AliForwardMultiplicityBase : public AliAnalysisTaskSE
{
public:
  /** 
   * @{ 
   * @name Interface methods 
   */
  /** 
   * Initialize the task 
   * 
   */
  virtual void Init() { fFirstEvent = true; }
  /** 
   * Create output objects 
   * 
   */
  virtual void UserCreateOutputObjects();
  /** 
   * Process each event 
   *
   * @param option Not used
   */  
  virtual void UserExec(Option_t* option) = 0;
  /** 
   * End of job
   * 
   * @param option Not used 
   */
  virtual void Terminate(Option_t* option);
  /** 
   * @} 
   */
  /** 
   * Configure this task via a macro 
   * 
   * @param macro Macro to configure va 
   * 
   * @return true on success, false otherwise
   */
  virtual Bool_t Configure(const char* macro="ForwardAODConfig.C");
  /** 
   * Print information 
   * 
   * @param option Not used
   */
  virtual void Print(Option_t* option="") const;
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
   * Get reference to the EventInspector algorithm 
   * 
   * @return Reference to AliFMDEventInspector object 
   */
  virtual AliFMDEventInspector& GetEventInspector() = 0;
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
   * Get reference to the EventInspector algorithm 
   * 
   * @return Reference to AliFMDEventInspector object 
   */
  virtual const AliFMDEventInspector& GetEventInspector() const = 0;
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
   * Overload super class method for setting debug level to call our
   * SetDebug member function.
   * 
   * @param dbg Debug level (0: no output, 1: essentials, 3: a whole lot)
   */
  virtual void SetDebugLevel(Int_t dbg) 
  { 
    AliAnalysisTaskSE::SetDebugLevel(dbg); 
    SetDebug(dbg);
  }
  /** 
   * Set whether to make separate branches for each ring.  If enabled
   * there will be 5 additional branches on the AOD tree - each
   * holding a TH2D object of the charged particle multiplicity in
   * @f$(\eta,\varphi)@f$ bins for that event.
   * 
   * @param use If true, make separate branches for each ring. 
   */
  void SetStorePerRing(Bool_t use) { fStorePerRing = use; }
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
  : AliAnalysisTaskSE(), 
    fEnableLowFlux(true), 
    fFirstEvent(true),
    fStorePerRing(false),
    fList(0),
    fHData(0),
    fHistos(),
    fAODFMD(),
    fAODEP(),
    fRingSums(),
    fCorrManager(0)
  {}
  /** 
   * Copy constructor 
   * 
   * @param o Object to copy from 
   */
  AliForwardMultiplicityBase(const AliForwardMultiplicityBase& o)
    : AliAnalysisTaskSE(o),
      fEnableLowFlux(o.fEnableLowFlux), 
      fFirstEvent(o.fFirstEvent),
      fStorePerRing(o.fStorePerRing),
      fList(o.fList),
      fHData(o.fHData),
      fHistos(o.fHistos),
      fAODFMD(o.fAODFMD),
      fAODEP(o.fAODEP),
      fRingSums(o.fRingSums),
      fCorrManager(o.fCorrManager)
  {
    DefineOutput(1, TList::Class());
    DefineOutput(2, TList::Class());
  }
  /** 
   * Assignment operator 
   * 
   * @param o Object to assign from 
   * 
   * @return Reference to this object 
   */
  AliForwardMultiplicityBase& operator=(const AliForwardMultiplicityBase& o);
  /** 
   * Check if all needed corrections are there and accounted for.  If not,
   * do a Fatal exit 
   * 
   * @param what Which corrections is needed
   * 
   * @return true if all present, false otherwise
   */  
  Bool_t CheckCorrections(UInt_t what) const;
  /** 
   * Read corrections
   * 
   * 
   * @param pe  On return, the eta axis
   * @param pv  On return ,the vertex axis 
   * @param mc  True assume MC input
   * 
   * @return true ons succcss
   */
  virtual Bool_t ReadCorrections(const TAxis*& pe, 
				 const TAxis*& pv,
				 Bool_t mc=false,
				 Bool_t sat=false);
  /**
   * Get the ESD event. IF this is the first event, initialise
   *
   * @return Pointer to ESD event structore 
   */
  virtual AliESDEvent* GetESDEvent();
  /** 
   * Initialise the sub objects and stuff.  Called on first event
   *
   * @return false on errors 
   */
  virtual Bool_t SetupForData();
  /** 
   * Initialize members based on eta and vertex axis - only available
   * after first event - called from SetupForData.
   * 
   * @param pe @f$\eta@f$ axis
   * @param pv Interaction point Z-coordinate axis 
   */
  virtual void InitMembers(const TAxis* pe, const TAxis* pv);
  /**
   * Create output branches - called from UserCreateOutputObjects
   */
  virtual void CreateBranches(AliAODHandler* ah);
  /**
   * Mark this event as one to store in the AOD 
   * 
   */
  virtual void MarkEventForStore() const;
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
  Bool_t                 fEnableLowFlux;// Whether to use low-flux specific code
  Bool_t                 fFirstEvent;   // Whether the event is the first seen 
  Bool_t                 fStorePerRing; // Store each ring on separate branch
  TList*                 fList;         // Output list
  TH2D*                  fHData;        // Summed 1/Nd^2N_{ch}/dphideta
  AliForwardUtil::Histos fHistos;       // Cache histograms 
  AliAODForwardMult      fAODFMD;       // Output object
  AliAODForwardEP        fAODEP;        // Output object
  AliForwardUtil::Histos fRingSums;     // Cache histograms 
private:
  /**
   * A pointer to the corrections manager.  This is here to make the
   * corrections manager persistent - that is, when we write the
   * analysis train to a file (as done in PROOF) we should also write
   * down the corrections mananger.   This pointer ensures that. 
   * 
   */
  AliForwardCorrectionManager* fCorrManager; // Pointer to corrections manager

  ClassDef(AliForwardMultiplicityBase,4) // Forward multiplicity class
};

#endif

// Local Variables:
//  mode: C++
// End:

