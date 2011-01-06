// 
// Calculate the multiplicity in the forward regions event-by-event 
// 
#ifndef ALIFORWARDMCMULTIPLICITYTASK_H
#define ALIFORWARDMCMULTIPLICITYTASK_H
#include "AliForwardMultiplicityBase.h"
#include "AliForwardUtil.h"
#include "AliFMDEventInspector.h"
#include "AliFMDEnergyFitter.h"
#include "AliFMDMCSharingFilter.h"
#include "AliFMDMCDensityCalculator.h"
#include "AliFMDMCCorrections.h"
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
 * @ingroup pwg2_forward_tasks
 * @ingroup pwg2_forward_mc
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
   * @return Reference to AliFMDCorrections object 
   */
  AliFMDCorrections& GetCorrections() { return fCorrections; }
  /**
   * Get reference to the HistCollector algorithm 
   * 
   * @return Reference to AliFMDHistCollector object 
   */
  AliFMDHistCollector& GetHistCollector() { return fHistCollector; }
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
   */
  virtual void   InitializeSubs();

  TH2D*                  fHData;        // Summed 1/Nd^2N_{ch}/dphideta
  AliESDFMD              fESDFMD;       // Sharing corrected ESD object
  AliForwardUtil::Histos fHistos;       // Cache histograms 
  AliAODForwardMult      fAODFMD;       // Output object
  AliESDFMD              fMCESDFMD;     // MC 'Sharing corrected' ESD object
  AliForwardUtil::Histos fMCHistos;     // MC Cache histograms 
  AliAODForwardMult      fMCAODFMD;     // MC Output object

  AliFMDEventInspector      fEventInspector;    // Algorithm
  AliFMDEnergyFitter        fEnergyFitter;      // Algorithm
  AliFMDMCSharingFilter     fSharingFilter;     // Algorithm
  AliFMDMCDensityCalculator fDensityCalculator; // Algorithm
  AliFMDMCCorrections       fCorrections;       // Algorithm
  AliFMDHistCollector       fHistCollector;     // Algorithm

  TList* fList; // Output list 

  ClassDef(AliForwardMCMultiplicityTask,1) // Forward multiplicity class
};

#endif
// Local Variables:
//  mode: C++
// End:

