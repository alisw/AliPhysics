// 
// Calculate the qa in the forward regions event-by-event 
// 
#ifndef ALIFORWARDQATASK_H
#define ALIFORWARDQATASK_H
/**
 * @file   AliForwardQATask.h
 * @author Christian Holm Christensen <cholm@dalsgaard.hehi.nbi.dk>
 * @date   Wed Mar 23 14:06:42 2011
 * 
 * @brief  
 * 
 * 
 * @ingroup pwg2_forward_aod
 */
#include <AliAnalysisTaskSE.h>
#include "AliFMDEventInspector.h"
#include "AliFMDSharingFilter.h"
#include "AliFMDDensityCalculator.h"
#include "AliFMDEnergyFitter.h"
#include <AliESDFMD.h>
class AliForwardCorrectionManager;
class AliESDEvent;
class TH2D;
class TList;
class TAxis;

/** 
 * Calculate the QA in the forward regions
 * 
 * @par Inputs: 
 *   - AliESDEvent 
 *
 * @par Outputs: 
 *   - Histograms 
 * 
 * @par Histograms 
 *   
 * @par Corrections used 
 * 
 * @ingroup pwg2_forward_tasks
 * 
 */
class AliForwardQATask : public AliAnalysisTaskSE
{
public:
  /** 
   * Constructor 
   * 
   * @param name Name of task 
   */
  AliForwardQATask(const char* name);
  /** 
   * Constructor
   */
  AliForwardQATask();
  /** 
   * Copy constructor 
   * 
   * @param o Object to copy from 
   */
  AliForwardQATask(const AliForwardQATask& o);
  /** 
   * Assignment operator 
   * 
   * @param o Object to assign from 
   * 
   * @return Reference to this object 
   */
  AliForwardQATask& operator=(const AliForwardQATask& o);
  /** 
   * @{ 
   * @name Interface methods 
   */
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
   * Get reference to the SharingFilter algorithm 
   * 
   * @return Reference to AliFMDSharingFilter object 
   */
  AliFMDSharingFilter& GetSharingFilter() { return fSharingFilter; }
  /**
   * Get reference to the DensityCalculator algorithm 
   * 
   * @return Reference to AliFMDDensityCalculator object 
   */
  AliFMDDensityCalculator& GetDensityCalculator() { return fDensityCalculator; }
  /**
   * Get reference to the EventInspector algorithm 
   * 
   * @return Reference to AliFMDEventInspector object 
   */
  const AliFMDEventInspector& GetEventInspector() const { return fEventInspector; }
  /**
   * Get reference to the EnergyFitter algorithm 
   * 
   * @return Reference to AliFMDEnergyFitter object 
   */
  const AliFMDEnergyFitter& GetEnergyFitter() const { return fEnergyFitter; }
  /**
   * Get reference to the SharingFilter algorithm 
   * 
   * @return Reference to AliFMDSharingFilter object 
   */
  const AliFMDSharingFilter& GetSharingFilter() const { return fSharingFilter; }
  /**
   * Get reference to the DensityCalculator algorithm 
   * 
   * @return Reference to AliFMDDensityCalculator object 
   */
  const AliFMDDensityCalculator& GetDensityCalculator() const { return fDensityCalculator; }
  /** 
   * @} 
   */
  /** 
   * Set debug level 
   * 
   * @param dbg Debug level
   */
  void SetDebug(Int_t dbg);
  void Print(Option_t* option="") const;
protected: 
  /** 
   * Check if all needed corrections are there and accounted for.  If not,
   * do a Fatal exit 
   * 
   * @param what Which corrections is needed
   * 
   * @return true if all present, false otherwise
   */  
  Bool_t CheckCorrections(UInt_t what) const;
  /**
   * Read corrections
   *
   */
  virtual Bool_t ReadCorrections(const TAxis*& pe, 
				 const TAxis*& pv,
				 Bool_t mc=false);
  /**
   * Get the ESD event. IF this is the first event, initialise.  If
   * initialisation failes, return a null pointer. 
   */
  virtual AliESDEvent* GetESDEvent();
  /** 
   * Initialise the sub objects and stuff.  Called on first event 
   * @return false on error. 
   */
  virtual Bool_t  InitializeSubs();

  Bool_t fEnableLowFlux;// Whether to use low-flux specific code
  Bool_t fFirstEvent;   // Whether the event is the first seen 
  /**
   * A pointer to the corrections manager.  This is here to make the
   * corrections manager persistent - that is, when we write the
   * analysis train to a file (as done in PROOF) we should also write
   * down the corrections mananger.   This pointer ensures that. 
   */
  AliForwardCorrectionManager* fCorrManager; // Pointer to corrections manager

  AliESDFMD              fESDFMD;       // Sharing corrected ESD object
  AliForwardUtil::Histos fHistos;       // Cache histograms 

  AliFMDEventInspector    fEventInspector;    // Algorithm
  AliFMDEnergyFitter      fEnergyFitter;      // Algorithm
  AliFMDSharingFilter     fSharingFilter;     // Algorithm
  AliFMDDensityCalculator fDensityCalculator; // Algorithm

  TList* fList; // Output list 
  Int_t fDebug; // Debug flag

  ClassDef(AliForwardQATask,1) // Forward QA class
};

#endif
// Local Variables:
//  mode: C++
// End:

