//
// Calculate the flow in the forward regions using the Q cumulants method
//
#ifndef AliForwardMCClosure_H
#define AliForwardMCClosure_H
/**
 * @file AliForwardMCClosure.h
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
class AliMCParticle;
class THn;
class AliTrackReference;
class TParticle;
class TCutG;
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
class AliForwardMCClosure : public AliAnalysisTaskSE
{
public:
  /** 
   * Constructor 
   */
  AliForwardMCClosure();
  /** 
   * Constructor
   * 
   * @param name Name of task 
   */
  AliForwardMCClosure(const char* name);
  /**
   * Destructor
   */
  virtual ~AliForwardMCClosure() {}

  /** 
   * Copy constructor 
   * 
   * @param o Object to copy from 
   */
  AliForwardMCClosure(const AliForwardMCClosure& o);

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


// Check if a given particle itself hit the FMD. If so, return the
  // (first) track reference of such a hit
  AliTrackReference* IsHitFMD(AliMCParticle* p);

// Check if a given particle itself hit the FMD. If so, return the
  // (first) track reference of such a hit
  AliTrackReference* IsHitTPC(AliMCParticle* p);


  /** 
   * End of job
   * 
   * @param option Not used 
   */
  virtual void Terminate(Option_t *option);

  //private:
  AliAODEvent*            fAOD;           //! input event
  TList*                  fOutputList;    //! output list
  TList*    fStdQCList; //! 
  TList*    fGFList; //!
  //TList* fRefList; //!
  //TList* fDiffList; //! 
  TList* fEventList; //! 
  TRandom fRandom;
  
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

  ClassDef(AliForwardMCClosure, 1); // Analysis task for flow analysis
};

#endif
// Local Variables:
//   mode: C++ 
// End:
