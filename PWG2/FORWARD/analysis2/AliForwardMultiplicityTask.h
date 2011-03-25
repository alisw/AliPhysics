// 
// Calculate the multiplicity in the forward regions event-by-event 
// 
#ifndef ALIFORWARDMULTIPLICITYTASK_H
#define ALIFORWARDMULTIPLICITYTASK_H
/**
 * @file   AliForwardMultiplicityTask.h
 * @author Christian Holm Christensen <cholm@dalsgaard.hehi.nbi.dk>
 * @date   Wed Mar 23 14:06:42 2011
 * 
 * @brief  
 * 
 * 
 * @ingroup pwg2_forward_aod
 */
#include "AliForwardMultiplicityBase.h"
#include "AliForwardUtil.h"
#include "AliFMDEventInspector.h"
#include "AliFMDSharingFilter.h"
#include "AliFMDDensityCalculator.h"
#include "AliFMDCorrector.h"
#include "AliFMDHistCollector.h"
#include "AliAODForwardMult.h"
#include "AliFMDEnergyFitter.h"
#include <AliESDFMD.h>
class AliESDEvent;
class TH2D;
class TList;

/** 
 * Calculate the multiplicity in the forward regions event-by-event 
 * 
 * @par Inputs: 
 *   - AliESDEvent 
 *
 * @par Outputs: 
 *   - AliAODForwardMult 
 * 
 * @par Histograms 
 *   
 * @par Corrections used 
 * 
 * @ingroup pwg2_forward_tasks
 * @ingroup pwg2_forward_aod
 * 
 */
class AliForwardMultiplicityTask : public AliForwardMultiplicityBase
{
public:
  /** 
   * Constructor 
   * 
   * @param name Name of task 
   */
  AliForwardMultiplicityTask(const char* name);
  /** 
   * Constructor
   */
  AliForwardMultiplicityTask();
  /** 
   * Copy constructor 
   * 
   * @param o Object to copy from 
   */
  AliForwardMultiplicityTask(const AliForwardMultiplicityTask& o);
  /** 
   * Assignment operator 
   * 
   * @param o Object to assign from 
   * 
   * @return Reference to this object 
   */
  AliForwardMultiplicityTask& operator=(const AliForwardMultiplicityTask& o);
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
   * Get reference to the Corrections algorithm 
   * 
   * @return Reference to AliFMDCorrector object 
   */
  AliFMDCorrector& GetCorrections() { return fCorrections; }
  /**
   * Get reference to the HistCollector algorithm 
   * 
   * @return Reference to AliFMDHistCollector object 
   */
  AliFMDHistCollector& GetHistCollector() { return fHistCollector; }
  /**
   * Get reference to the EventInspector algorithm 
   * 
   * @return Reference to AliFMDEventInspector object 
   */
  const AliFMDEventInspector& GetEventInspector() const { return fEventInspector; }
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
   * Get reference to the Corrections algorithm 
   * 
   * @return Reference to AliFMDCorrector object 
   */
  const AliFMDCorrector& GetCorrections() const { return fCorrections; }
  /**
   * Get reference to the HistCollector algorithm 
   * 
   * @return Reference to AliFMDHistCollector object 
   */
  const AliFMDHistCollector& GetHistCollector() const { return fHistCollector; }
  /** 
   * @} 
   */
  /** 
   * Set debug level 
   * 
   * @param dbg Debug level
   */
  void SetDebug(Int_t dbg);
protected: 
  /** 
   * Initialise the sub objects and stuff.  Called on first event 
   * 
   */
  virtual void   InitializeSubs();

  TH2D*                  fHData;        // Summed 1/Nd^2N_{ch}/dphideta
  AliESDFMD              fESDFMD;       // Sharing corrected ESD object
  AliForwardUtil::Histos fHistos;       // Cache histograms 
  AliAODForwardMult      fAODFMD;       // Output object

  AliFMDEventInspector    fEventInspector;    // Algorithm
  AliFMDSharingFilter     fSharingFilter;     // Algorithm
  AliFMDDensityCalculator fDensityCalculator; // Algorithm
  AliFMDCorrector         fCorrections;       // Algorithm
  AliFMDHistCollector     fHistCollector;     // Algorithm

  TList* fList; // Output list 

  ClassDef(AliForwardMultiplicityTask,1) // Forward multiplicity class
};

#endif
// Local Variables:
//  mode: C++
// End:

