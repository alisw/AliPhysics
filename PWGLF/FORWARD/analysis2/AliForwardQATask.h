// 
// Calculate the qa in the forward regions event-by-event 
// 
#ifndef ALIFORWARDQATASK_H
#define ALIFORWARDQATASK_H
/**
 * @file   AliForwardQATask.h
 * @author Christian Holm Christensen <cholm@dalsgaard.hehi.nbi.dk>
 * @date   Wed Mar 23 14:06:42 2011
 * 
 * @brief  
 * 
 * 
 * @ingroup pwglf_forward_aod
 */
#include "AliBaseESDTask.h"
#include "AliFMDEventInspector.h"
#include "AliFMDESDFixer.h"
#include "AliFMDSharingFilter.h"
#include "AliFMDDensityCalculator.h"
#include "AliFMDEnergyFitter.h"
#include <AliESDFMD.h>
class AliESDEvent;
class TH2D;
class TAxis;

/** 
 * Calculate the QA in the forward regions
 * 
 * @par Inputs: 
 *   - AliESDEvent 
 *
 * @par Outputs: 
 *   - Histograms 
 * 
 * @par Histograms 
 *   
 * @par Corrections used 
 * 
 * @ingroup pwglf_forward_tasks
 * 
 */
class AliForwardQATask : public AliBaseESDTask
{
public:
  /** 
   * Constructor 
   * 
   * @param name Name of task 
   */
  AliForwardQATask(const char* name);
  /** 
   * Constructor
   */
  AliForwardQATask();
  /** 
   * @{ 
   * @name Interface methods 
   */
  /** 
   * Called when initialising the train. 
   * 
   * @return Always true
   */
  virtual Bool_t Setup();
  /** 
   * Book output objects. Derived class should define this to book
   * output objects on the processing output list @c fList before the
   * actual event processing.  This is called on the master and on
   * each slave.
   * 
   * If this member function returns false, the execution is stopped
   * with a fatal signal.
   *
   * @return true on success. 
   */
  virtual Bool_t Book();
  /** 
   * Called after reading in the first event. Here we can setup stuff
   * depending on the conditions we're running under.
   * 
   * @return true on success.  If this returns false, then we turn the
   * task into a zombie and we do no more processing.
   */
  virtual Bool_t PreData(const TAxis& vertex, const TAxis& eta);
  
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
   * End of job
   * 
   * @return true on success
   */
  virtual Bool_t Finalize();
  /** 
   * @} 
   */
  /** 
   * @{ 
   * @name Default axes 
   */
  /** 
   * Set the default eta axis to use in case we didn't get one from
   * the read-in corretions.  Override this if the sub class should go
   * on even without a valid eta axis from the corrections (e.g. QA
   * task)
   * 
   * @return null
   */
  virtual TAxis* DefaultEtaAxis() const;
  /** 
   * Set the default eta axis to use in case we didn't get one from
   * the read-in corretions.  Override this if the sub class should go
   * on even without a valid eta axis from the corrections (e.g. QA
   * task)
   * 
   * @return null
   */
  virtual TAxis* DefaultVertexAxis() const;
  /* @} */
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
   * Get reference to the ESDFixer algorithm 
   * 
   * @return Reference to AliFMDESDFixer object 
   */
  AliFMDESDFixer& GetESDFixer() { return fESDFixer; }
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
   * Get reference to the EventInspector algorithm 
   * 
   * @return Reference to AliFMDEventInspector object 
   */
  const AliFMDEventInspector& GetEventInspector() const { return fEventInspector; }
  /**
   * Get reference to the ESDFixer algorithm 
   * 
   * @return Reference to AliFMDESDFixer object 
   */
  const AliFMDESDFixer& GetESDFixer() const { return fESDFixer; }
  /**
   * Get reference to the EnergyFitter algorithm 
   * 
   * @return Reference to AliFMDEnergyFitter object 
   */
  const AliFMDEnergyFitter& GetEnergyFitter() const { return fEnergyFitter; }
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
   * @} 
   */
  /** 
   * Set debug level 
   * 
   * @param dbg Debug level
   */
  void SetDebug(Int_t dbg);
  /** 
   * Print information 
   * 
   * @param option Not used
   */
  void Print(Option_t* option="") const;
  /** 
   * Check if we're running over MC data
   * 
   * @return true if the event inspector thinks it's MC
   */
  Bool_t IsMC() const { return GetEventInspector().IsMC(); }
protected: 
  /** 
   * Copy constructor 
   * 
   * @param o Object to copy from 
   */
  AliForwardQATask(const AliForwardQATask& o);
  /** 
   * Assignment operator 
   * 
   * @param o Object to assign from 
   * 
   * @return Reference to this object 
   */
  AliForwardQATask& operator=(const AliForwardQATask& o);

  Bool_t                  fEnableLowFlux;// Whether to use low-flux code
  AliESDFMD               fESDFMD;       // Sharing corrected ESD object
  AliForwardUtil::Histos  fHistos;       // Cache histograms 
  AliFMDEventInspector    fEventInspector;    // Algorithm
  AliFMDESDFixer          fESDFixer;          // Algorithm
  AliFMDEnergyFitter      fEnergyFitter;      // Algorithm
  AliFMDSharingFilter     fSharingFilter;     // Algorithm
  AliFMDDensityCalculator fDensityCalculator; // Algorithm

  ClassDef(AliForwardQATask,4) // Forward QA class
};

#endif
// Local Variables:
//  mode: C++
// End:

