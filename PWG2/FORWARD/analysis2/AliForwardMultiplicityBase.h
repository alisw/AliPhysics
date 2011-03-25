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
 * @ingroup pwg2_forward_aod
 */
#include <AliAnalysisTaskSE.h>
class AliFMDEventInspector;
class AliFMDEnergyFitter;
class AliFMDSharingFilter;
class AliFMDDensityCalculator;
class AliFMDCorrector;
class AliFMDHistCollector;
class AliForwardCorrectionManager;
class AliESDEvent;
class TH2D;
class TList;
class TTree;
class TAxis;

/** 
 * @mainpage ALICE PWG2 Forward Multiplcity Analysis 
 * 
 * This is the analysis code for analysis of the Forward data. 
 * 
 * @par Code overview 
 * 
 * See the <a href="modules.html">Modules</a> page 
 * 
 * @par Run.sh script 
 * 
 * @verbatim 
 * Usage: Run.sh [OPTIONS]
 * 
 * Do Pass1 and Pass2 on ESD files in current directory.  
 * 
 * Options:
 * 	-h,--help		This help                  
 * 	-n,--events N		Number of events            (-1)
 * 	-1,--pass1 		Run pass 1, only AOD        (0)
 * 	-2,--pass2		Run pass 2, only Hists      (0)
 * 	-3,--pass3		Draw results                (0)
 * 	-v,--vz-min CM          Minimum value of vz         (-10)
 * 	-V,--vz-max CM          Maximum value of vz         (10)
 * 	-t,--trigger TYPE       Select trigger TYPE         (INEL)
 * 	-b,--batch              Do batch processing         (0)
 * 	-P,--proof NWORKERS	Run in PROOF(Lite) mode     (0)
 * 	-M,--mc			Run over MC data            (0)
 * 	-g,--gdb		Run in GDB mode    	    (0)
 * 	-E,--eloss		Run energy loss script      
 *      -r,--rebin              Rebin factor                (1)
 *      -C,--use-centrality     Run centrality task         (0)
 * 	-O,--show-older		Show older data	            (0)
 * 	-J,--show-published	Show ALICE published data   (1)
 * 	-R,--show-ratios	Show ratios to other data   (1)
 * 	-Z,--show-asymmetry	Show asymmetry 		    (1)
 * 	-S,--scheme SCHEME	Normalisation scheme	    (full)
 * 	-T,--title STRING       Title on plots              ()
 * 
 * TYPE is a comma or space separated list of 
 *  
 *   INEL	      Inelastic triggers (V0A|V0C|SPD)
 *   INEL>0      As above + N_ch > 0 in -0.5<eta<+0.5
 *   NSD         Non-single diffractive ((VOA&VOC)|N_ch > 5 -1.9<eta<+1.9)
 * 
 * SCHEME is a comma or space separated list of 
 * 
 *   NONE          No event-level normalization except trivial one 
 *   EVENTLEVEL    Event-level normalization 
 *   ALTEVENTLEVEL Event-level normalization (alternative version)
 *   BACKGROUND    Not implemented yet 
 *   SHAPE         Shape correction 
 *   FULL          Same as EVENTLEVEL,BACKGROUND,SHAPE
 *   ALTFULL       Same as ALTEVENTLEVEL,BACKGROUND,SHAPE
 * 
 * If NWORKERS is 0, then the analysis will be run in local mode.
 * @endverbatim
 */
/** 
 * @defgroup pwg2_forward PWG2 Forward analysis
 *
 * Code to do the multiplicity analysis in the forward psuedo-rapidity
 * regions
 *
 */
/** 
 * @defgroup pwg2_forward_tasks Tasks
 *
 * Code to do the multiplicity analysis in the forward psuedo-rapidity
 * regions
 *
 * @ingroup pwg2_forward 
 */
/** 
 * @defgroup pwg2_forward_topical Topical
 */
/** 
 * @defgroup pwg2_forward_aod AOD
 * @ingroup pwg2_forward_topical
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
 * @ingroup pwg2_forward_tasks
 * @ingroup pwg2_forward_aod
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
  virtual void UserCreateOutputObjects() = 0;
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
  virtual void Terminate(Option_t* option) = 0;
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
   * @} 
   */
  virtual void SetDebug(Int_t dbg) = 0;
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
      fCorrManager(o.fCorrManager)
  {}
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
   */
  virtual Bool_t ReadCorrections(const TAxis*& pe, 
				 const TAxis*& pv,
				 Bool_t mc=false);
  /**
   * Get the ESD event. IF this is the first event, initialise
   */
  virtual AliESDEvent* GetESDEvent();
  /** 
   * Initialise the sub objects and stuff.  Called on first event
   *
   */
  virtual void InitializeSubs() = 0;
  /**
   * Mark this event as one to store in the AOD 
   * 
   */
  virtual void MarkEventForStore() const;

  Bool_t                 fEnableLowFlux;// Whether to use low-flux specific code
  Bool_t                 fFirstEvent;   // Whether the event is the first seen 
private:
  /**
   * A pointer to the corrections manager.  This is here to make the
   * corrections manager persistent - that is, when we write the
   * analysis train to a file (as done in PROOF) we should also write
   * down the corrections mananger.   This pointer ensures that. 
   * 
   */
  AliForwardCorrectionManager* fCorrManager; // Pointer to corrections manager

  ClassDef(AliForwardMultiplicityBase,2) // Forward multiplicity class
};

#endif
// Local Variables:
//  mode: C++
// End:

