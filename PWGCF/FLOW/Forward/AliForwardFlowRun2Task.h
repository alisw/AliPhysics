//
// Calculate the flow in the forward regions using the Q cumulants method
//
#ifndef AliForwardFlowRun2Task_H
#define AliForwardFlowRun2Task_H
/**
 * @file AliForwardFlowRun2Task.h
 * @author Freja Thoresen <freja.thoresen@cern.ch>
 *
 * @brief
 *
 * @ingroup pwgcf_forward_flow
 */
#include "AliAnalysisTaskSE.h"
#include "TString.h"
#include <TH2D.h>
#include "TRandom.h"
#include "AliForwardSettings.h"
#include "AliForwardTaskValidation.h"
#include "AliForwardFlowUtil.h"
#include "AliForwardFlowResultStorage.h"
#include "AliForwardGenericFramework.h"
#include "AliForwardWeights.h"

class TH2D;
class THn;

/**
 * @defgroup pwglf_forward_tasks_flow Flow tasks
 *
 * Code to do flow
 *
 * @ingroup pwglf_forward_tasks
 */
/**
 * Calculate the flow in the forward regions using the Q cumulants method
 *
 * @par Inputs:
 *   - AliAODEvent
 *
 * Outputs:
 *   - forward_flow.root
 *
 * @ingroup pwglf_forward_flow
 *
 */


class AliForwardFlowRun2Task : public AliAnalysisTaskSE
{
public:
  /**
   * Constructor
   */
  AliForwardFlowRun2Task();
  /**
   * Constructor
   *
   * @param name Name of task
   */
  AliForwardFlowRun2Task(const char* name);
  /**
   * Destructor
   */
  virtual ~AliForwardFlowRun2Task() {}

  /**
   * Copy constructor
   *
   * @param o Object to copy from
   */
  AliForwardFlowRun2Task(const AliForwardFlowRun2Task& o);

  /**
   * @{
   * @name Task interface methods
   */

  /**
   * Create output objects
   */
  virtual void UserCreateOutputObjects();

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

  //private:
  AliVEvent* fAOD;      //! input event
  TList* fOutputList;   //! output list
  TList* fAnalysisList; //!
  TList* fReferenceList; //!
  TList* fStandardList; //!
  TList* fMixedList; //!
  TList* fEventList;    //!

  TRandom fRandom; //!

  TH2D*   centralDist; //!
  TH2D*   refDist;     //!
  TH2D*   forwardDist; //!

  AliForwardFlowResultStorage* fStorage; //!

  // A class combining all the settings for this analysis
  AliForwardSettings fSettings;

  // Utility class for filling histograms
  AliForwardFlowUtil fUtil; 
  
  // Class for flow calculations using the Generic Framework
  AliForwardGenericFramework fCalculator;

  AliForwardWeights fWeights;
  ClassDef(AliForwardFlowRun2Task, 1); // Analysis task for flow analysis
};

#endif
// Local Variables:
//   mode: C++
// End:
