//
// Calculate the flow in the forward regions using the Q cumulants method
//
#ifndef AliForwardNUATask_H
#define AliForwardNUATask_H
/**
 * @file AliForwardNUATask.h
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
#include "AliForwardFlowRun2Settings.h"
#include "AliEventCuts.h"
#include <TF1.h>

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
class AliForwardNUATask : public AliAnalysisTaskSE
{
public:
  /** 
   * Constructor 
   */
  AliForwardNUATask();
  /** 
   * Constructor
   * 
   * @param name Name of task 
   */
  AliForwardNUATask(const char* name);
  /**
   * Destructor
   */
  virtual ~AliForwardNUATask() {}

  /** 
   * Copy constructor 
   * 
   * @param o Object to copy from 
   */
  AliForwardNUATask(const AliForwardNUATask& o);

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
  AliForwardFlowRun2Settings fSettings;

AliEventCuts fEventCuts;
TF1 *fMultTOFLowCut; //!
TF1 *fMultTOFHighCut; //!
TF1 *fMultCentLowCut; //!

  enum {
    kTPCOnly = 128, // TPC only tracks
    kHybrid = 768, // TPC only tracks
    kphiAcceptanceBin = 21 // phi acceptance bin in the FMD histogram (dNdetadphi)
  };

  ClassDef(AliForwardNUATask, 1); // Analysis task for flow analysis
};

#endif
// Local Variables:
//   mode: C++ 
// End:
