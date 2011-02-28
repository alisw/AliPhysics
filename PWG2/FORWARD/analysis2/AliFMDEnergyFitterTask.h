// Histogram and fit the energy loss distributions for the FMD
// 
// Wraps AliFMDEnergyFitter 
#ifndef ALIFMDENERGYFITTERTASK_H
#define ALIFMDENERGYFITTERTASK_H
#include <AliAnalysisTaskSE.h>
#include "AliForwardUtil.h"
#include "AliFMDEventInspector.h"
#include "AliFMDEnergyFitter.h"
#include <AliESDFMD.h>
#include <TH1I.h>
class AliESDEvent;
class TH2D;
class TList;
class TTree;


/** 
 * Histogram and fit the energy loss distributions for the FMD
 * 
 * @par Inputs: 
 *   - AliESDEvent 
 *
 * @par Outputs: 
 *   - None
 * 
 * @par Histograms 
 *   
 * @par Corrections used 
 *   - None
 * 
 * @ingroup pwg2_forward_tasks
 * 
 */
class AliFMDEnergyFitterTask : public AliAnalysisTaskSE
{
public:
  /** 
   * Constructor 
   * 
   * @param name Name of task 
   */
  AliFMDEnergyFitterTask(const char* name);
  /** 
   * Constructor
   */
  AliFMDEnergyFitterTask();
  /** 
   * Copy constructor 
   * 
   * @param o Object to copy from 
   */
  AliFMDEnergyFitterTask(const AliFMDEnergyFitterTask& o);
  /** 
   * Assignment operator 
   * 
   * @param o Object to assign from 
   * 
   * @return Reference to this object 
   */
  AliFMDEnergyFitterTask& operator=(const AliFMDEnergyFitterTask& o);
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
   * @} 
   */
  /** 
   * Set the debug level 
   * 
   * @param dbg Debug level
   */
  void SetDebug(Int_t dbg);

  void SetBLow(Float_t b) {fbLow = b;}
  void SetBHigh(Float_t b) {fbHigh = b;}

protected: 
  /** 
   * Initialise the sub objects and stuff.  Called on first event 
   * 
   */
  virtual void   InitializeSubs();

  Bool_t               fFirstEvent;     // Whether the event is the first seen 
  AliFMDEventInspector fEventInspector; // Algorithm
  AliFMDEnergyFitter   fEnergyFitter;   // Algorithm
  TList*               fList;           // Output list 
  Float_t              fbLow;
  Float_t              fbHigh;
  
  ClassDef(AliFMDEnergyFitterTask,1) // Forward multiplicity class
};

#endif
// Local Variables:
//  mode: C++
// End:

