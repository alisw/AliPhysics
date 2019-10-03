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
 * @ingroup pwglf_forward_aod
 * 
 */
#include "AliCentralMultiplicityTask.h"
#include "AliSPDMCTrackDensity.h"
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
 * @ingroup pwglf_forward_tasks
 * @ingroup pwglf_forward_aod
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
   * Create output objects 
   * 
   * 
   * @return true on success
   */
  virtual Bool_t Book();
  /** 
   * Creatre ouput objects
   * 
   * @param ah AOD output handler 
   */
  virtual void CreateBranches(AliAODHandler* ah);
  /** 
   * Set-up for data, called before first event 
   * 
   * @param v Vertex axis
   * @param e @f$\eta@f$ axis 
   * 
   * @return true on success
   */
  virtual Bool_t PreData(const TAxis& v, const TAxis& e);
  /** 
   * End of job
   * 
   * 
   * @return true on success
   */
  virtual Bool_t PreEvent();
  /** 
   * Process each event 
   *
   * @param esd ESD event
   *
   * @return true on success
   */  
  virtual Bool_t Event(AliESDEvent& esd);
  /** 
   * Print information 
   * 
   * @param option Not used
   */
  virtual void Print(Option_t* option="") const;
  /** 
   * Return the track density calculator 
   * 
   * @return Track density calculator 
   */
  const AliSPDMCTrackDensity& GetTrackDensity() const { return fTrackDensity; }
  /** 
   * Return the track density calculator 
   * 
   * @return Track density calculator 
   */
  AliSPDMCTrackDensity& GetTrackDensity() { return fTrackDensity; }

protected: 
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
  AliSPDMCTrackDensity   fTrackDensity;     // Calculate N_ch,incl from MC
  AliAODCentralMult      fAODMCCentral;     // Output object
  ClassDef(AliCentralMCMultiplicityTask,2)  // Forward multiplicity class
};

#endif
// Local Variables:
//  mode: C++
// End:

