// 
// Base class for classes that calculate the multiplicity in the
// SPD clusters event-by-event
// 
#ifndef ALICENTRALMULTIPLICITYTASK_H
#define ALICENTRALMULTIPLICITYTASK_H
#include <AliAnalysisTaskSE.h>
#include "AliForwardUtil.h"
#include "AliAODCentralMult.h"
#include "AliCentralCorrAcceptance.h"
#include "AliCentralCorrSecondaryMap.h"
//class AliForwardCorrectionManager;
class AliESDEvent;
class TH2D;
class TList;
class TTree;


/** 
 * @mainpage ALICE PWG2 Forward Multiplcity Analysis 
 */
/** 
 * @defgroup pwg2_forward PWG2 Forward analysis
 *
 * Code to do the multiplicity analysis in the central pseudo-rapidity
 * regions
 *
 */
/** 
 * @defgroup pwg2_forward_tasks Tasks
 *
 * Code to do the multiplicity analysis in the central pseudo-rapidity
 * regions
 *
 * @ingroup pwg2_forward 
 */
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

    ClassDef(Manager,1); // Manager of data 
  };

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
  
  TH2D*                  fData;          //sum histogram if needed
  TList*                 fList;          //Output List for diagnostics
  AliAODCentralMult      fAODCentral;    // Output object
  Manager                fManager;       //Manager object for corrections
  ClassDef(AliCentralMultiplicityTask,1) // Forward multiplicity class
};

#endif
// Local Variables:
//  mode: C++
// End:

