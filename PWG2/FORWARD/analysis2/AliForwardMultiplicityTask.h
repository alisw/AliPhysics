#ifndef ALIROOT_PWG2_FORWARD_ALIFORWARDMULTIPLICITY_H
#define ALIROOT_PWG2_FORWARD_ALIFORWARDMULTIPLICITY_H
#include <AliAnalysisTaskSE.h>
#include "AliForwardUtil.h"
#include "AliFMDEventInspector.h"
#include "AliFMDEnergyFitter.h"
#include "AliFMDSharingFilter.h"
#include "AliFMDDensityCalculator.h"
#include "AliFMDCorrections.h"
#include "AliFMDHistCollector.h"
#include "AliAODForwardMult.h"
#include "AliFMDEnergyFitter.h"
#include <AliESDFMD.h>
#include <TH1I.h>
class AliFMDAnaParameters;
class AliESDEvent;
class TH2D;
class TList;
class TTree;


/** 
 * @mainpage ALICE PWG2 Forward Multiplcity Analysis 
 */
/** 
 * @defgroup pwg2_forward_analysis PWG2 Forward analysis
 *
 * Code to do the multiplicity analysis in the forward psuedo-rapidity
 * regions
 *
 */
/** 
 * Calculate the multiplicity in the forward regions event-by-event 
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
 * @ingroup pwg2_forward_analysis 
 * 
 */
class AliForwardMultiplicityTask : public AliAnalysisTaskSE
{
public:
  /** 
   * Constructor 
   * 
   * @param name Name of task 
   */
  AliForwardMultiplicityTask(const char* name);
  /** 
   * Constructor
   */
  AliForwardMultiplicityTask();
  /** 
   * Copy constructor 
   * 
   * @param o Object to copy from 
   */
  AliForwardMultiplicityTask(const AliForwardMultiplicityTask& o);
  /** 
   * Assignment operator 
   * 
   * @param o Object to assign from 
   * 
   * @return Reference to this object 
   */
  AliForwardMultiplicityTask& operator=(const AliForwardMultiplicityTask& o);
  /** 
   * @{ 
   * @name Interface methods 
   */
  /** 
   * Initialize the task 
   * 
   */
  virtual void Init();
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
  virtual void UserExec(Option_t* option);
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
   * Print information 
   * 
   * @param option Not used
   */
  void Print(Option_t* option="") const;
  /** 
   * Whether to enable low-flux code 
   * 
   * @param use IF true, enable low-flux code 
   */
  void SetEnableLowFlux(Bool_t use=true) { fEnableLowFlux = use; }
  /** 
   * @{ 
   * @name Access to sub-algorithms 
   */
  /**
   * Get reference to the EventInspector algorithm 
   * 
   * @return Reference to AliFMDEventInspector object 
   */
  AliFMDEventInspector& GetEventInspector() { return fEventInspector; }
  /**
   * Get reference to the EnergyFitter algorithm 
   * 
   * @return Reference to AliFMDEnergyFitter object 
   */
  AliFMDEnergyFitter& GetEnergyFitter() { return fEnergyFitter; }
  /**
   * Get reference to the SharingFilter algorithm 
   * 
   * @return Reference to AliFMDSharingFilter object 
   */
  AliFMDSharingFilter& GetSharingFilter() { return fSharingFilter; }
  /**
   * Get reference to the DensityCalculator algorithm 
   * 
   * @return Reference to AliFMDDensityCalculator object 
   */
  AliFMDDensityCalculator& GetDensityCalculator() { return fDensityCalculator; }
  /**
   * Get reference to the Corrections algorithm 
   * 
   * @return Reference to AliFMDCorrections object 
   */
  AliFMDCorrections& GetCorrections() { return fCorrections; }
  /**
   * Get reference to the HistCollector algorithm 
   * 
   * @return Reference to AliFMDHistCollector object 
   */
  AliFMDHistCollector& GetHistCollector() { return fHistCollector; }
  /** 
   * @} 
   */
  void SetDebug(Int_t dbg);
protected: 
  /** 
   * Initialise the sub objects and stuff.  Called on first event 
   * 
   */
  virtual void   InitializeSubs();
  /** 
   * Mark this event as one to store in the AOD 
   * 
   */
  virtual void MarkEventForStore() const;

  Bool_t                 fEnableLowFlux;// Whether to use low-flux specific code
  TH2D*                  fHData;        // Summed 1/Nd^2N_{ch}/dphideta
  Bool_t                 fFirstEvent;   // Whether the event is the first seen 
  AliESDFMD              fESDFMD;       // Sharing corrected ESD object
  AliForwardUtil::Histos fHistos;       // Cache histograms 
  AliAODForwardMult      fAODFMD;       // Output object

  AliFMDEventInspector    fEventInspector;    // Algorithm
  AliFMDEnergyFitter      fEnergyFitter;      // Algorithm
  AliFMDSharingFilter     fSharingFilter;     // Algorithm
  AliFMDDensityCalculator fDensityCalculator; // Algorithm
  AliFMDCorrections       fCorrections;       // Algorithm
  AliFMDHistCollector     fHistCollector;     // Algorithm

  TList* fList; // Output list 

  ClassDef(AliForwardMultiplicityTask,1) // Forward multiplicity class
};

#endif
// Local Variables:
//  mode: C++
// End:

