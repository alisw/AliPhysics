//
// Calculate the event plane in the forward regions using the FMD 
//
#ifndef ALIFMDEVENTPLANETASK_H
#define ALIFMDEVENTPLANETASK_H
/**
 * @file AliFMDEventPlaneTask.h
 * @author Alexander Hansen
 * @date   Tue Feb 21 2012
 * 
 * @brief
 * 
 * 
 * @ingroup pwglf_forward_flow
 */
#include "AliBaseAODTask.h"
#include "AliFMDEventPlaneFinder.h"
class AliAODForwardMult;
class TH1D;
class TH2D;

/**
 * Calculate the event plane in the forward regions using the FMD
 *
 * @par Inputs:
 *   - AliAODEvent
 *
 * Outputs:
 *   - AnalysisResults.root
 *
 * @ingroup pwglf_forward_tasks_flow
 * @ingroup pwglf_forward_flow
 *
 *
 */
class AliFMDEventPlaneTask : public AliBaseAODTask
{
public:
  /** 
   * Constructor 
   */
  AliFMDEventPlaneTask();
  /** 
   * Constructor
   * 
   * @param name Name of task 
   */
  AliFMDEventPlaneTask(const char* name);
  /**
   * Destructor
   */
  virtual ~AliFMDEventPlaneTask() {}
  /** 
   * @{ 
   * @name Task interface methods 
   */
  /** 
   * Create output objects 
   * 
   * @return true on success
   */
  virtual Bool_t Book();
  /** 
   * Process each event 
   *
   * @param aod AOD Event
   * 
   * @return true on success
   */  
  virtual Bool_t Event(AliAODEvent& aod);
  /** 
   * End of job
   * 
   * @return true on success
   */
  virtual Bool_t Finalize();
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
protected:
  /** 
   * Copy constructor 
   * 
   * @param o Object to copy from 
   */
  AliFMDEventPlaneTask(const AliFMDEventPlaneTask& o);
  /** 
   * Assignment operator 
   * 
   * @return Reference to this object 
   */
  AliFMDEventPlaneTask& operator=(const AliFMDEventPlaneTask&);

  AliFMDEventPlaneFinder  fEventPlaneFinder; //  Eventplane finder for the FMD
  TH1D*                   fHistVertexSel;    //  Diagnostics histogram

  ClassDef(AliFMDEventPlaneTask, 3); // Analysis task for FMD analysis
};
 
#endif
// Local Variables:
//   mode: C++ 
// End:
