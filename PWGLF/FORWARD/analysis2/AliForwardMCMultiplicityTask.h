// 
// Calculate the multiplicity in the forward regions event-by-event 
// 
#ifndef ALIFORWARDMCMULTIPLICITYTASK_H
#define ALIFORWARDMCMULTIPLICITYTASK_H
/**
 * @file   AliForwardMCMultiplicityTask.h
 * @author Christian Holm Christensen <cholm@dalsgaard.hehi.nbi.dk>
 * @date   Wed Mar 23 14:06:13 2011
 * 
 * @brief  
 * 
 * 
 * @ingroup pwglf_forward_aod
 */
#include "AliForwardMultiplicityBase.h"
#include "AliForwardUtil.h"
#include "AliFMDMCEventInspector.h"
#include "AliFMDMCSharingFilter.h"
#include "AliFMDMCDensityCalculator.h"
#include "AliFMDMCCorrector.h"
#include "AliFMDHistCollector.h"
#include "AliAODForwardMult.h"
#include "AliAODForwardEP.h"
#include "AliFMDEnergyFitter.h"
#include "AliFMDEventPlaneFinder.h"
#include <AliESDFMD.h>
class AliESDEvent;
class TH2D;
class TList;

/** 
 * Calculate the multiplicity in the forward regions event-by-event 
 * 
 * @par Inputs: 
 *   - AliESDEvent 
 *   - Kinematics
 *   - Track references
 *
 * @par Outputs: 
 *   - AliAODForwardMult 
 * 
 * @par Histograms 
 *   
 * @par Corrections used 
 * 
 * @ingroup pwglf_forward_tasks
 * @ingroup pwglf_forward_mc
 * @ingroup pwglf_forward_aod
 * 
 */
class AliForwardMCMultiplicityTask : public AliForwardMultiplicityBase
{
public:
  /** 
   * Constructor 
   * 
   * @param name Name of task 
   */
  AliForwardMCMultiplicityTask(const char* name);
  /** 
   * Constructor
   */
  AliForwardMCMultiplicityTask();
  /** 
   * Copy constructor 
   * 
   * @param o Object to copy from 
   */
  AliForwardMCMultiplicityTask(const AliForwardMCMultiplicityTask& o);
  /** 
   * Assignment operator 
   * 
   * @param o Object to assign from 
   * 
   * @return Reference to this object 
   */
  AliForwardMCMultiplicityTask& 
  operator=(const AliForwardMCMultiplicityTask& o);
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
   * Process only primary MC tracks 
   * 
   * @param use 
   */
  void SetOnlyPrimary(Bool_t use);
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
   * Get reference to the EventPlaneFinder algorithm 
   * 
   * @return Reference to AliFMDEventPlaneFinder object 
   */
  AliFMDEventPlaneFinder& GetEventPlaneFinder() { return fEventPlaneFinder; }
  /**
   * Get reference to the EventPlaneFinder algorithm 
   * 
   * @return Reference to AliFMDEventPlaneFinder object 
   */
  const AliFMDEventPlaneFinder& GetEventPlaneFinder() const { return fEventPlaneFinder; }
  /** 
   * @} 
   */
  /** 
   * Set debug level 
   * 
   * @param dbg debug level
   */
  void SetDebug(Int_t dbg);
protected: 
  /** 
   * Initialise the sub objects and stuff.  Called on first event 
   * 
   * @return false on errors
   */
  virtual Bool_t SetupForData();

  TH2D*                  fHData;        // Summed 1/Nd^2N_{ch}/dphideta
  AliESDFMD              fESDFMD;       // Sharing corrected ESD object
  AliForwardUtil::Histos fHistos;       // Cache histograms 
  AliAODForwardMult      fAODFMD;       // Output object
  AliAODForwardEP        fAODEP;       // Output object
  AliESDFMD              fMCESDFMD;     // MC 'Sharing corrected' ESD object
  AliForwardUtil::Histos fMCHistos;     // MC Cache histograms 
  AliAODForwardMult      fMCAODFMD;     // MC Output object
  AliForwardUtil::Histos fRingSums;     // Cache histograms 
  AliForwardUtil::Histos fMCRingSums;   // Cache histograms 
  TH2D*                  fPrimary;      // Per event primary particles 

  AliFMDMCEventInspector    fEventInspector;    // Algorithm
  AliFMDMCSharingFilter     fSharingFilter;     // Algorithm
  AliFMDMCDensityCalculator fDensityCalculator; // Algorithm
  AliFMDMCCorrector         fCorrections;       // Algorithm
  AliFMDHistCollector       fHistCollector;     // Algorithm
  AliFMDEventPlaneFinder    fEventPlaneFinder;  // Algorithm

  TList* fList; // Output list 
  TList* fListVertexBins; // list of the signal  in vertex bin	

  ClassDef(AliForwardMCMultiplicityTask,3) // Forward multiplicity class
};

#endif
// Local Variables:
//  mode: C++
// End:

