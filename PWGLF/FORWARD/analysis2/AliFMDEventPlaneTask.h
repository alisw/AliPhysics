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
#include "AliAnalysisTaskSE.h"
#include "AliFMDEventPlaneFinder.h"
class AliAODForwardMult;
class TH1D;
class TH2D;

 /**
 * @defgroup pwg2_forward_tasks_flow Flow tasks 
 * @ingroup pwg2_forward_tasks
 */
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
class AliFMDEventPlaneTask : public AliAnalysisTaskSE
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
   */
  virtual void UserCreateOutputObjects();
  /**
   * Initialize the task
   */
  virtual void Init() {} 
  /** 
   * Process each event 
   *
   * @param option Not used
   */  
  virtual void UserExec(Option_t *option);
  /** 
   * End of job
   * 
   * @param option Not used 
   */
  virtual void Terminate(Option_t *option);
 /**
  * Check AODForwardMult object for trigger, vertex and centrality
  * returns true if event is OK
  * 
  * @param const aodfm
  * 
  * @return Bool_t 
  */
  Bool_t AODCheck(const AliAODForwardMult* aodfm);
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
  /*
   * Set MC input flag - currently does nothing special
   *
   * @ param mc MC input flag
   */
  void SetMCInput(Bool_t mc = true) { fMC = mc; }

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

  TList*                  fSumList;          //  Sum list
  TList*                  fOutputList;       //  Output list
  AliAODEvent*            fAOD;              //  AOD event
  Bool_t                  fMC;               //  MC input?
  AliFMDEventPlaneFinder  fEventPlaneFinder; //  Eventplane finder for the FMD
  Float_t                 fZvertex;	     //  Z vertex
  Double_t                fCent;             //  Centrality
  TH1D*                   fHistCent;         //  Diagnostics histogram
  TH1D*                   fHistVertexSel;    //  Diagnostics histogram
  TH1D*                   fHistVertexAll;    //  Diagnostics histogram

  ClassDef(AliFMDEventPlaneTask, 1); // Analysis task for FMD analysis
};
 
#endif
// Local Variables:
//   mode: C++ 
// End:
