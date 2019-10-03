// 
// Calculate the multiplicity in the forward regions event-by-event 
// 
#ifndef ALIFORWARDMULTIPLICITYTASK_H
#define ALIFORWARDMULTIPLICITYTASK_H
/**
 * @file   AliForwardMultiplicityTask.h
 * @author Christian Holm Christensen <cholm@dalsgaard.hehi.nbi.dk>
 * @date   Wed Mar 23 14:06:42 2011
 * 
 * @brief  
 * 
 * 
 * @ingroup pwglf_forward_aod
 */
#include "AliForwardMultiplicityBase.h"
#include "AliForwardUtil.h"
#include "AliFMDEventInspector.h"
#include "AliMultEventClassifier.h"
#include "AliFMDESDFixer.h"
#include "AliFMDSharingFilter.h"
#include "AliFMDDensityCalculator.h"
#include "AliFMDCorrector.h"
#include "AliFMDHistCollector.h"
// #include "AliFMDEnergyFitter.h"
#include "AliFMDEventPlaneFinder.h"
#include <AliESDFMD.h>
class AliESDEvent;
class TH2D;
class TList;
class TH3D;	

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
 * @ingroup pwglf_forward_tasks
 * @ingroup pwglf_forward_aod
 * 
 */
class AliForwardMultiplicityTask : public AliForwardMultiplicityBase
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
   * @{ 
   * @name Interface methods 
   */
  /** 
   * Called on first event _before_ reading corrections.  Here, the
   * user class can do additional checking to see if the some (more or
   * less) corrections are needed.
   * 
   * @param esd Event 
   */
  virtual void PreCorrections(const AliESDEvent* esd);
  /** 
   * Called before processing a single event - should not do anything
   * but clear data, etc.
   * 
   * @return true on success
   */
  virtual Bool_t PreEvent();
  /** 
   * Process each event 
   *
   * @param esd Event
   * 
   * @return true on success
   */  
  virtual Bool_t Event(AliESDEvent& esd);
  /** 
   * @} 
   */
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
   * Get a reference to the event inspector. User must override this
   * to return proper object
   * 
   * @return Reference to the event inspector 
   */
  // AliMultEventClassifier& GetMultEventClassifier() { return fMultEventClassifier; }
  /**
   * Get reference to the ESDFixer algorithm 
   * 
   * @return Reference to AliFMDESDFixer object 
   */
  AliFMDESDFixer& GetESDFixer() { return fESDFixer; }
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
   * @return Reference to AliFMDCorrector object 
   */
  AliFMDCorrector& GetCorrections() { return fCorrections; }
  /**
   * Get reference to the HistCollector algorithm 
   * 
   * @return Reference to AliFMDHistCollector object 
   */
  AliFMDHistCollector& GetHistCollector() { return fHistCollector; }
  /**
   * Get reference to the EventInspector algorithm 
   * 
   * @return Reference to AliFMDEventInspector object 
   */
  const AliFMDEventInspector& GetEventInspector() const { return fEventInspector; }
  /** 
   * Get a reference to the event inspector. User must override this
   * to return proper object
   * 
   * @return Reference to the event inspector 
   */
  // const AliMultEventClassifier& GetMultEventClassifier() const { return fMultEventClassifier; }
  /**
   * Get reference to the ESDFixer algorithm 
   * 
   * @return Reference to AliFMDESDFixer object 
   */
  const AliFMDESDFixer& GetESDFixer() const { return fESDFixer; }
  /**
   * Get reference to the SharingFilter algorithm 
   * 
   * @return Reference to AliFMDSharingFilter object 
   */
  const AliFMDSharingFilter& GetSharingFilter() const { return fSharingFilter; }
  /**
   * Get reference to the DensityCalculator algorithm 
   * 
   * @return Reference to AliFMDDensityCalculator object 
   */
  const AliFMDDensityCalculator& GetDensityCalculator() const { return fDensityCalculator; }
  /**
   * Get reference to the Corrections algorithm 
   * 
   * @return Reference to AliFMDCorrector object 
   */
  const AliFMDCorrector& GetCorrections() const { return fCorrections; }
  /**
   * Get reference to the HistCollector algorithm 
   * 
   * @return Reference to AliFMDHistCollector object 
   */
  const AliFMDHistCollector& GetHistCollector() const { return fHistCollector; }
 /**
   * Get reference to the EventPlaneFinder algorithm 
   * 
   * @return Reference to AliFMDEventPlaneFinder object 
   */
  AliFMDEventPlaneFinder& GetEventPlaneFinder() { return fEventPlaneFinder; }
  /**
   * Get reference to the EventPlaneFinder algorithm 
   * 
   * @return Reference to AliFMDEventPlaneFinder object 
   */
  const AliFMDEventPlaneFinder& GetEventPlaneFinder() const { return fEventPlaneFinder; }
  /** 
   * @} 
   */
  /** 
   * Set whether to make a timing histogram 
   * 
   * @param enable 
   */
  virtual void SetDoTiming(Bool_t enable=true);
protected: 
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

  AliESDFMD               fESDFMD;            // Sharing corrected ESD object
  AliFMDEventInspector    fEventInspector;    // Algorithm
  // AliMultEventClassifier  fMultEventClassifier;//Event class
  AliFMDESDFixer          fESDFixer;          // Algorithm
  AliFMDSharingFilter     fSharingFilter;     // Algorithm
  AliFMDDensityCalculator fDensityCalculator; // Algorithm
  AliFMDCorrector         fCorrections;       // Algorithm
  AliFMDHistCollector     fHistCollector;     // Algorithm
  AliFMDEventPlaneFinder  fEventPlaneFinder;  // Algorithm

  ClassDef(AliForwardMultiplicityTask,8) // Forward multiplicity class
};

#endif
// Local Variables:
//  mode: C++
// End:

