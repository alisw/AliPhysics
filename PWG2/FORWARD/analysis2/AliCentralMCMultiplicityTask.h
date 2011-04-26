// 
// Base class for classes that calculate the multiplicity in the
// SPD clusters event-by-event
// 
#ifndef ALICENTRALMCMULTIPLICITYTASK_H
#define ALICENTRALMCMULTIPLICITYTASK_H
/**
 * @file   AliCentralMCMultiplicityTask.h
 * @author Hans Hjersing Dalsgaard
 * @date   Wed Mar 23 14:00:03 2011
 * 
 * @brief  
 * 
 * @ingroup pwg2_forward_aod
 * 
 */
#include "AliCentralMultiplicityTask.h"
#include "AliSPDMCTrackDensity.h"
//class AliForwardCorrectionManager;
class AliESDEvent;
class AliMCEvent;

/** 
 * Class that calculates the multiplicity in the
 * central region event-by-event
 * 
 * @par Inputs: 
 *   - AliESDEvent 
 *
 * @par Outputs: 
 *   - 2 AliAODCentralMult (one from data and one from MC)
 * 
 * @par Histograms 
 *   
 * @par Corrections used 
 * 
 * @ingroup pwg2_forward_tasks
 * @ingroup pwg2_forward_aod
 * 
 */
class AliCentralMCMultiplicityTask : public AliCentralMultiplicityTask
{
public:
  /** 
   * @{ 
   * @name Interface methods 
   */
   /** 
   * Constructor 
   * 
   * @param name Name of task 
   */
  AliCentralMCMultiplicityTask(const char* name); 
  /** 
   * Constructor 
   *
   * Reserved for ROOT's I/O system - do not use
   */
  AliCentralMCMultiplicityTask();
  /** 
   * Copy constructor 
   * 
   * @param o Object to copy from 
   */
  AliCentralMCMultiplicityTask(const AliCentralMCMultiplicityTask& o);
  /** 
   * Assignment operator 
   * 
   * @param o Object to assign from 
   * 
   * @return Reference to this object 
   */
  AliCentralMCMultiplicityTask& operator=(const AliCentralMCMultiplicityTask&o);
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
   * Print information 
   * 
   * @param option Not used
   */
  virtual void Print(Option_t* option="") const;

protected: 
  AliSPDMCTrackDensity   fTrackDensity;     // Calculate N_ch,incl
					    // from MC
  AliAODCentralMult      fAODMCCentral;     // Output object
  ClassDef(AliCentralMCMultiplicityTask,1)  // Forward multiplicity class
};

#endif
// Local Variables:
//  mode: C++
// End:

