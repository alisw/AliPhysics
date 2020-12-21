//
// Calculate the flow in the forward regions using the Q cumulants method
//
#ifndef AliForwardNUETask_H
#define AliForwardNUETask_H
/**
 * @file AliForwardNUETask.h
 * @author Freja Thoresen <freja.thoresen@cern.ch>
 *
 * @brief
 *
 * @ingroup pwgcf_forward_flow
 */
#include "AliAnalysisTaskSE.h"
#include "TString.h"
#include <TH2D.h>
#include <TH3F.h>
#include "TRandom.h"
#include "AliForwardSettings.h"
#include "AliEventCuts.h"
#include <TF1.h>
#include "AliForwardFlowUtil.h"
#include <AliForwardTaskValidation.h>

class AliAODForwardMult;
class TH2D;
class AliESDEvent;
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
class AliForwardNUETask : public AliAnalysisTaskSE
{
public:
  /**
   * Constructor
   */
  AliForwardNUETask();
  /**
   * Constructor
   *
   * @param name Name of task
   */
  AliForwardNUETask(const char* name);
  /**
   * Destructor
   */
  virtual ~AliForwardNUETask() {}

  /**
   * Copy constructor
   *
   * @param o Object to copy from
   */
  AliForwardNUETask(const AliForwardNUETask& o);

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
  AliAODEvent*            fAOD;           //! input event
  TList*                  fOutputList;    //! output list
  TList* fEventList; //!

  // A class combining all the settings for this analysis
  AliForwardSettings fSettings;
  AliForwardFlowUtil fUtil;

Bool_t nua_mode;

  enum {
    kTPCOnly = 128, // TPC only tracks
    kHybrid = 768, // Hybrid tracks
    kGlobalOnly = 32, // Global tracks
    kGlobal = 96, // Global tracks
    kphiAcceptanceBin = 21 // phi acceptance bin in the FMD histogram (dNdetadphi)
  };

  ClassDef(AliForwardNUETask, 1); // Analysis task for flow analysis
};

#endif
// Local Variables:
//   mode: C++
// End:
