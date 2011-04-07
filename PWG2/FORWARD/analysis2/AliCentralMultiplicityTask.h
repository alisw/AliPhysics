// 
// Base class for classes that calculate the multiplicity in the
// SPD clusters event-by-event
// 
#ifndef ALICENTRALMULTIPLICITYTASK_H
#define ALICENTRALMULTIPLICITYTASK_H
/**
 * @file   AliCentralMultiplicityTask.h
 * @author Hans Hjersing Dalsgaard
 * @date   Wed Mar 23 14:00:03 2011
 * 
 * @brief  
 * 
 * @ingroup pwg2_forward_aod
 * 
 */
#include <AliAnalysisTaskSE.h>
#include "AliFMDEventInspector.h"
#include "AliAODCentralMult.h"
#include "AliCentralCorrAcceptance.h"
#include "AliCentralCorrSecondaryMap.h"
//class AliForwardCorrectionManager;
class AliESDEvent;
class AliMultiplicity;
class TH2D;
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
 *   - AliAODCentralMult 
 * 
 * @par Histograms 
 *   
 * @par Corrections used 
 * 
 * @ingroup pwg2_forward_tasks
 * @ingroup pwg2_forward_aod
 * 
 */
class AliCentralMultiplicityTask : public AliAnalysisTaskSE
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
  AliCentralMultiplicityTask(const char* name); 
  /** 
   * Constructor 
   *
   * Reserved for ROOT's I/O system - do not use
   */
  AliCentralMultiplicityTask();
  /** 
   * Copy constructor 
   * 
   * @param o Object to copy from 
   */
  AliCentralMultiplicityTask(const AliCentralMultiplicityTask& o);
  /** 
   * Assignment operator 
   * 
   * @param o Object to assign from 
   * 
   * @return Reference to this object 
   */
  AliCentralMultiplicityTask& operator=(const AliCentralMultiplicityTask& o);
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
  /** 
   * Set whether to use the secondary corrections 
   * 
   * @param use Whether to use secondary corrections 
   */
  virtual void SetUseSecondary(Bool_t use) { fUseSecondary = use; }
  /** 
   * Set whether to use the acceptance corrections 
   * 
   * @param use Whether to use acceptance corrections 
   */
  virtual void SetUseAcceptance(Bool_t use) { fUseAcceptance = use; }

  AliFMDEventInspector& GetInspector() { return fInspector; }
  const AliFMDEventInspector& GetInspector() const { return fInspector; }

  //__________________________________________________________________
  /**
   * Manager of corrections 
   *
   * This is a small class to fetch corrections for secondaries and
   * dead channels.
   * 
   */
  class Manager 
  {
  public:
    /** 
     * Constructor
     * 
     */
    Manager();
    /** 
     * Copy constructor 
     * 
     * @param o 
     */
    Manager(const Manager& o);
    /** 
     * Destructor
     */
    virtual ~Manager() {}
    /** 
     * Assignment operator 
     * 
     * @param o Object to assign from 
     * 
     * @return Reference to this object 
     */
    Manager& operator=(const Manager& o);
    
    /** 
     * Initialize 
     * 
     * @param sys    Collision system (1: pp, 2: PbPb)
     * @param sNN    Center of mass energy per nucleon pair [GeV]
     * @param field  Magnetic field [kG]
     */
    void Init(UShort_t sys, UShort_t sNN, Short_t field);
    
    /** 
     * Is initialized 
     * 
     */
    Bool_t IsInit() { return fIsInit; }
    
    
    /** 
     * Get the acceptance path
     * 
     * @return 
     */
    const char* GetAcceptancePath() const {return fAcceptancePath.Data(); }
    /** 
     * Get the secondary path 
     * 
     * @return 
     */
    const char* GetSecMapPath() const {return fSecMapPath.Data(); }
    /** 
     * Set the path to the acceptance maps 
     * 
     * @param path PAth to object file 
     */
    void SetAcceptancePath(const char* path) {fAcceptancePath=path; }
    /** 
     * Set the path to the secondary maps 
     * 
     * @param path Path to object files 
     */
    void  SetSecMapPath(const char* path) {fSecMapPath=path; }
    /** 
     * Get full path name to object file 
     * 
     * @param what   What to get 
     * @param sys    Collision system
     * @param sNN    Center of mass energy 
     * @param field  Magnetic field 
     * 
     * @return 
     */
    const char* GetFullFileName(UShort_t what, UShort_t sys, UShort_t sNN, 
				Short_t  field) const;
    /** 
     * Get the acceptance object name 
     * 
     * @return 
     */
    const char* GetAcceptanceName() const {return fAcceptanceName.Data(); }
    /** 
     * Get the secondary object name 
     * 
     * @return 
     */
    const char* GetSecMapName() const {return fSecMapName.Data(); }
    
    /** 
     * Get the secondary map
     * 
     * @param vtxbin 
     * 
     * @return 
     */
    TH2D* GetSecMapCorrection(UShort_t vtxbin) const;
    /** 
     * Get the acceptance correction 
     * 
     * @param vtxbin 
     * 
     * @return 
     */
    TH1D* GetAcceptanceCorrection(UShort_t vtxbin) const;
    /** 
     * Get the secondary correction map object 
     */
    AliCentralCorrSecondaryMap* GetSecMap() const { return fSecmap; }

    void Print(Option_t* option="") const;
  private:
    /** 
     * Get the full path name 
     * 
     * @param what   What to get
     * @param sys    Collision system
     * @param sNN    Center of mass energy 
     * @param field  Magnetic field 
     * 
     * @return 
     */
    const char* GetFileName(UShort_t what, UShort_t sys, UShort_t sNN,
			    Short_t field) const;

    
    TString                     fAcceptancePath; // Path to acceptance 
    TString                     fSecMapPath;     // Path to secondary map
    AliCentralCorrAcceptance*   fAcceptance;     // Acceptance 
    AliCentralCorrSecondaryMap* fSecmap;         // Secindary map
    TString                     fAcceptanceName; // Acceptance name
    TString                     fSecMapName;     // Secindary name
    Bool_t                      fIsInit;         // Are we init

    ClassDef(Manager,1); // Manager of data 
  };
  /** 
   * Get the ESD event and initialise manager on first event if not
   * done already
   * 
   * @return Pointer to valid ESD event object 
   */
  virtual AliESDEvent* GetESDEvent();
  /** 
   * Mark this event for storage in AOD output
   * 
   */
  virtual void MarkEventForStore() const;
  /** 
   * Process the ESD SPD information 
   * 
   * @param hist    Histogram to fill
   * @param spdmult SPD multiplicity object
   */
  virtual void ProcessESD(TH2D& hist, const AliMultiplicity* spdmult) const;
  /** 
   * Corret the data 
   * 
   * @param hist    Histogram to correct
   * @param vtxbin  Vertex bin 
   */
  virtual void CorrectData(TH2D& hist, UShort_t vtxbin) const;
  /** 
   * Get a reference to the manager 
   * 
   * @return Reference to corrections manager 
   */
  Manager& GetManager() { return fManager; }
  /** 
   * Get a reference to the manager 
   * 
   * @return Reference to corrections manager 
   */
  const Manager& GetManager() const { return fManager; }


protected: 
  AliFMDEventInspector   fInspector;      // Inspect events 
  TH2D*                  fData;           // sum histogram if needed
  TList*                 fList;           // Output List for diagnostics
  AliAODCentralMult      fAODCentral;     // Output object
  Manager                fManager;        // Manager object for corrections
  Bool_t                 fUseSecondary;   // Whether to secondary map
  Bool_t                 fUseAcceptance;  // Whether to use acceptance corr.
  Bool_t                 fFirstEventSeen; // Have we seen first event     
  Int_t                  fIvz;            // Event's vertex bin 
  ClassDef(AliCentralMultiplicityTask,2)  // Forward multiplicity class
};

#endif
// Local Variables:
//  mode: C++
// End:

