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
   */
  AliCentralMultiplicityTask() 
    : AliAnalysisTaskSE(),
      fData(0),
      fList(0),
      fAODCentral(),
      fManager()
  {
    DefineOutput(1, TList::Class());
  }
  /** 
   * Copy constructor 
   * 
   * @param o Object to copy from 
   */
  AliCentralMultiplicityTask(const AliCentralMultiplicityTask& o)
    : AliAnalysisTaskSE(o),
      fData(o.fData),
      fList(o.fList),
      fAODCentral(o.fAODCentral),
      fManager(o.fManager)
  {
    DefineOutput(1, TList::Class());
  }
  /** 
   * Assignment operator 
   * 
   * @param o Object to assign from 
   * 
   * @return Reference to this object 
   */
  AliCentralMultiplicityTask& operator=(const AliCentralMultiplicityTask& o)
  {
    fData       = o.fData;
    fList       = o.fList;
    fAODCentral = o.fAODCentral;
    fManager    = o.fManager;
    
    
    DefineOutput(1, TList::Class());
    
    return *this;
  }
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
   * Init the task and the manager  
   * 
   * @param option Not used
   */
  void InitManager(UShort_t sys, 
		   UShort_t  sNN,
		   Short_t   field) {fManager.Init(sys, sNN, field);}
  /** 
   * @} 
   */
  /** 
   * Print information 
   * 
   * @param option Not used
   */
  virtual void Print(Option_t* option="") const;
  /** 
   * Set Path for acceptance 
   * 
   * @param path
   */
  void          SetAcceptancePath(const char* path) {fManager.SetAcceptancePath(path); }
  /** 
   * Set Path for Secondary Map 
   * 
   * @param path
   */
  
  void          SetSecMapPath(const char* path) {fManager.SetSecMapPath(path); }
  
  char*         GetFullFileName(UShort_t  what ,
				UShort_t  sys, 
				UShort_t  sNN,
				Short_t   field) {return fManager.GetFullFileName(what ,sys, sNN, field); }
  
  const char*   GetAcceptanceName() {return fManager.GetAcceptanceName(); }
  const char*   GetSecMapName() {return fManager.GetSecMapName(); }
  
  
  class Manager {
    
    // This is a small class to fetch corrections for secondaries and dead
    // channels.
    
  public:
    Manager();
    Manager(const Manager& o) :
      fAcceptancePath(o.fAcceptancePath),
      fSecMapPath(o.fSecMapPath),
      fAcceptance(o.fAcceptance),
      fSecmap(o.fSecmap),
      fAcceptanceName(o.fAcceptanceName),
      fSecMapName(o.fSecMapName) {}
    
  /** 
   * Assignment operator 
   * 
   * @param o Object to assign from 
   * 
   * @return Reference to this object 
   */
  Manager& operator=(const Manager& o)
  {
    fAcceptancePath = o.fAcceptancePath;
    fSecMapPath     = o.fSecMapPath;
    fAcceptance     = o.fAcceptance;
    fSecmap         = o.fSecmap;
    fAcceptanceName = o.fAcceptanceName;
    fSecMapName     = o.fSecMapName;
    return *this;
  }

    void          Init(UShort_t  sys, 
		       UShort_t  sNN,
		       Short_t   field);
    const char*   GetAcceptancePath() {return fAcceptancePath.Data(); }
    const char*   GetSecMapPath() {return fSecMapPath.Data(); }
    void          SetAcceptancePath(const char* path) {fAcceptancePath=path; }
    void          SetSecMapPath(const char* path) {fSecMapPath=path; }
    char*         GetFullFileName(UShort_t  what ,
				  UShort_t  sys, 
				  UShort_t  sNN,
				  Short_t   field) {return Form("%s/%s",
what == 0 ? GetSecMapPath() : GetAcceptancePath(), GetFileName(what, sys, sNN, field));}
    const char*   GetAcceptanceName() {return fAcceptanceName.Data(); }
    const char*   GetSecMapName() {return fSecMapName.Data(); }
    
    TH2D* GetSecMapCorrection(UShort_t vtxbin) {return fSecmap->GetCorrection(vtxbin);}
    TH1D* GetAcceptanceCorrection(UShort_t vtxbin) {return fAcceptance->GetCorrection(vtxbin);}
    
  private:
    
    
    const char*   GetFileName(UShort_t  what ,
			      UShort_t  sys, 
			      UShort_t  sNN,
			      Short_t   field);
    
    
    TString fAcceptancePath;
    TString fSecMapPath;
    AliCentralCorrAcceptance* fAcceptance;
    AliCentralCorrSecondaryMap*     fSecmap;
    TString fAcceptanceName;
    TString fSecMapName;

  };

protected: 
 
 
private:
  
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

