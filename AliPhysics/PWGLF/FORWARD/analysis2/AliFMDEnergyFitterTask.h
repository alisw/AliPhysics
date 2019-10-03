// Histogram and fit the energy loss distributions for the FMD
// 
// Wraps AliFMDEnergyFitter 
#ifndef ALIFMDENERGYFITTERTASK_H
#define ALIFMDENERGYFITTERTASK_H
/**
 * @file   AliFMDEnergyFitterTask.h
 * @author Christian Holm Christensen <cholm@dalsgaard.hehi.nbi.dk>
 * @date   Wed Mar 23 14:02:39 2011
 * 
 * @brief  
 * 
 * @ingroup pwglf_forward_eloss
 * 
 */
#include "AliBaseESDTask.h"
#include "AliFMDEventInspector.h"
#include "AliFMDEnergyFitter.h"
#include "AliFMDESDFixer.h"
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
 * @ingroup pwglf_forward_tasks
 * @ingroup pwglf_forward_eloss
 * 
 */
class AliFMDEnergyFitterTask : public AliBaseESDTask
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
   * @{ 
   * @name Interface methods 
   */
  /** 
   * Called on master when setting up the train. 
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
   * Called on first event _before_ reading corrections.  Here, the
   * user class can do additional checking to see if the some (more or
   * less) corrections are needed.
   * 
   * @param esd Event 
   */
  virtual void PreCorrections(const AliESDEvent* esd);
  /** 
   * Called after reading in the first event. Here we can setup stuff
   * depending on the conditions we're running under.
   * 
   * @return true on success.  If this returns false, then we turn the
   * task into a zombie and we do no more processing.
   */
  virtual Bool_t PreData(const TAxis& vertex, const TAxis& eta);
  /** 
   * Process each event 
   *
   * @param esd Event to analyse
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
   * Get reference to the EventInspector algorithm 
   * 
   * @return Reference to AliFMDEventInspector object 
   */
  const AliFMDEventInspector& GetEventInspector() const{return fEventInspector;}
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
   * @} 
   */
  /** 
   * @{ 
   * @name Settings 
   */
  /** 
   * Set the debug level 
   * 
   * @param dbg Debug level
   */
  void SetDebug(Int_t dbg);
  /** 
   * Set whether to only look at MB (INEL) data, so as to avoid 
   * bias from different trigger scalars. 
   * 
   * @param onlyMB if true, only analyse MB events
   */
  void SetOnlyMB(Bool_t onlyMB) { fOnlyMB = onlyMB; }
  /* @} */
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
protected: 
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

  virtual const char* DefaultConfig() const { return "elossFitConfig.C"; }

  AliFMDEventInspector fEventInspector; // Algorithm
  AliFMDESDFixer       fESDFixer;       // Algorithm
  AliFMDEnergyFitter   fEnergyFitter;   // Algorithm
  Bool_t               fOnlyMB;         // Only MB flag

  ClassDef(AliFMDEnergyFitterTask,4) // Forward multiplicity class
};

#endif
// Local Variables:
//  mode: C++
// End:

