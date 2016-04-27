 /**
 * @file   AliForwardTriggerBiasCorrection.h
 * @author Valentina Zaccolo
 * @date   Mon Feb  3 12:31:02 2013
 * 
 * @brief  
 * 
 * @ingroup pwglf_forward_multdist
 */
#ifndef ALIFORWARDTRIGGERBIASCORRECTION_H
#define ALIFORWARDTRIGGERBIASCORRECTION_H  
#include <AliBaseMultTask.h>
/**
 * Task to make the reponse matrices used by the multiplicity
 * distibution analysis
 * 
 * @ingroup pwglf_forward Tasks
 * @ingroup pwglf_forward_multdist
 */
class AliForwardTriggerBiasCorrection : public AliBaseMultTask
{
public:
  /** 
   * 
   * Default Constructor 
   * 
   */
  AliForwardTriggerBiasCorrection()
    : AliBaseMultTask()
  {}    
  /**
   *  Constructor
   */
  AliForwardTriggerBiasCorrection(const char* name)
    : AliBaseMultTask(name)
  {}
  /** 
   * Destructor
   */
  virtual ~AliForwardTriggerBiasCorrection(){}
  /**
   *  Embedded Class begins here
   */
  struct Bin : public AliBaseMultTask::Bin
  {
    /**
     * Default Constructor
     */
    Bin()
      : AliBaseMultTask::Bin(),
	fESDClass(0),
	fMCClass(0),
	fMCESDClass(0)
    {}
    /**
     * Constructor
     */
    Bin(Double_t etaLow, Double_t etaHigh)
      : AliBaseMultTask::Bin(etaLow, etaHigh),
	fESDClass(0),
	fMCClass(0),
	fMCESDClass(0)
    {}
    /**
     * Copy Constructor
     */ 
    Bin(const Bin& o)
      : AliBaseMultTask::Bin(o),
	fESDClass(0),
	fMCClass(0),
	fMCESDClass(0)
    {}
    /**
     * Assignment Operator
     */
    Bin& operator=(const Bin&) { return *this; }
    /**
     * Destructor
     */
    ~Bin() {}
    /**
     * Define outputs of a single eta bin
     */
    virtual void CreateOutputObjects(TList* cont,  Int_t max);
    /** 
     * Process a single eta bin
     * 
     * @param dndetaForward   Forward observations 
     * @param dndetaCentral   Central observations 
     * @param normForward     Acceptance 
     * @param normCentral     Acceptance 
     * @param dndetaMC        MC-truth distribution 
     * @param ipZ             Interaction point 
     * @param pileup          True if flagged as pile-up
     * @param selectedTrigger Is event selected
     * @param isMCClass       Is event MC NSD 
     * @param isESDClass      Is event real NSD 
     * @param aodevent        Full event 
     * @param minIPz          Least Z coordinate of IP 
     * @param maxIPz          Largest Z coordinate of IP 
     */
    virtual void Process(TH1D*              dndetaForward,
			 TH1D*              dndetaCentral,
			 TH1D*              normForward,
			 TH1D*              normCentral,
			 TH1D*              dndetaMC,
			 Double_t           ipZ,
			 Bool_t             pileup,
			 Bool_t             selectedTrigger,
			 Bool_t             isMCClass,
			 Bool_t             isESDClass,
			 const AliAODEvent& aodevent,
			 Double_t           minIPz,
			 Double_t           maxIPz);
    TH1D* fESDClass;   //! number of events found as ev. class
		       //! selected by the analysis vs. multiplicity
    TH1D* fMCClass;    //! number of events found as ev. class
		       //! selected by the MC truth vs. multiplicity
    TH1D* fMCESDClass; //! number of events found as ev. class
		       //! selected by both analysis and MC truth
		       //! vs. multiplicity
    ClassDef(Bin,2); // Manager of data 
  };
  /** 
   * Create a bin
   */
  AliBaseMultTask::Bin* MakeBin(Double_t etaLow, Double_t etaHigh);
protected:
  /**
   * Copy Constructor
   *
   */
  AliForwardTriggerBiasCorrection(const AliForwardTriggerBiasCorrection& o)
    : AliBaseMultTask(o)
  {}
  /**
   * Assignment operator
   *
   */
  AliForwardTriggerBiasCorrection& 
  operator=(const AliForwardTriggerBiasCorrection&) { return *this; }
  /** 
   * Check if this is OK class for reconstruction 
   * 
   * @return true if OK
   */
  Bool_t IsESDClass(AliAODForwardMult*) const;
  ClassDef(AliForwardTriggerBiasCorrection, 5); 
};

#endif
// Local Variables:
//   mode: C++
// End:
