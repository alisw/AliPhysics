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
//class AliForwardCorrectionManager;
class AliESDEvent;
class AliMCEvent;
class TH1D;
class TH2D;
class TH3D;
class TList;
class TTree;

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
  AliCentralMCMultiplicityTask& operator=(const AliCentralMCMultiplicityTask& o);
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
  /** 
   * Process the MC event
   * 
   * @param hist   Histogram to fill 
   * @param event  Event structure 
   */
  void ProcessMC(TH2D& hist, const AliMCEvent* event) const;

  AliAODCentralMult      fAODMCCentral;     // Output object
  Double_t               fMinR;             // Min radius 
  Double_t               fMaxR;             // Max radius 
  Double_t               fMinZ;             // Min z
  Double_t               fMaxZ;             // Max z
  TH2D*                  fRZ;               // Location in (r,z)
  TH3D*                  fXYZ;              // Location in (x,y,z)
  TH1D*                  fNRefs;            // Refs per track 
  ClassDef(AliCentralMCMultiplicityTask,1)  // Forward multiplicity class
};

#endif
// Local Variables:
//  mode: C++
// End:

