/**
 * @file   AliForwardMultiplicityDistribution.h
 * @author Christian Holm Christensen <cholm@master.hehi.nbi.dk>
 * @date   Thu Feb  7 01:02:42 2013
 * 
 * @brief  
 * 
 * 
 * @ingroup pwglf_forward_multdist
 */
#ifndef ALIFORWARDMULTIPLICITYDISTRIBUTION_H
#define ALIFORWARDMULTIPLICITYDISTRIBUTION_H
#include "AliBaseMultTask.h"

/**
 * Task to do the multiplicity distibution
 * 
 * @ingroup pwglf_forward Tasks
 * @ingroup pwglf_forward_multdist
 */
class AliForwardMultiplicityDistribution : public AliBaseMultTask
{
public:
  /**
   * Default Constructor
   */
  AliForwardMultiplicityDistribution() : AliBaseMultTask() {}
  /**
   * Constructor
   */
  AliForwardMultiplicityDistribution(const char* name)
    : AliBaseMultTask(name)
  {}
  /**
   * Destructor
   */
  virtual ~AliForwardMultiplicityDistribution(){}
  /**
   * Embedded Class begins here
   */
  struct Bin : public AliBaseMultTask::Bin
  {
    /**
     * Default Constructor
     */
    Bin() : AliBaseMultTask::Bin(), fHistPileUp(0) {}
    /**
     * Constructor
     */
    Bin(Double_t etaLow, Double_t etaHigh)
      : AliBaseMultTask::Bin(etaLow, etaHigh), fHistPileUp(0)
    {}    
    /**
     * Copy Constructor
     */    
    Bin(const Bin& o)
       : AliBaseMultTask::Bin(o), fHistPileUp(0)
    {}
    /**
     * Assignment operator
     */
    Bin& operator=(const Bin&){return*this;}
    /**
     * Destructor
     */    
    virtual ~Bin(){}
    /**
     *  Define outputs of a single eta bin
     */
    virtual void CreateOutputObjects(TList* cont,  Int_t max);
    /** 
     * Process a single eta bin
     * 
     * @param dndetaForward   Forward observations 
     * @param dndetaCentral   Central observations 
     * @param normForward     Acceptance 
     * @param normCentral     Acceptance 
     * @param mc              Primary "observations"
     * @param ipZ             Interaction point 
     * @param pileup          True if flagged as pile-up
     * @param selectedTrigger Is event selected
     * @param isMCNSDm        Is event MC NSD 
     * @param isESDNSD        Is event real NSD 
     * @param aodevent        Full event 
     * @param minIPz          Least Z coordinate of IP 
     * @param maxIPz          Largest Z coordinate of IP 
     */
    virtual void Process(TH1D*              dndetaForward,
			 TH1D*              dndetaCentral,
			 TH1D*              normForward,
			 TH1D*              normCentral,
			 TH1D*              mc,
			 Double_t           ipZ,
			 Bool_t             pileup,
			 Bool_t             selectedTrigger,
			 Bool_t             isMCNSDm,
			 Bool_t             isESDNSD,
			 const AliAODEvent& aodevent,
			 Double_t           minIPz,
			 Double_t           maxIPz);
    TH1D*    fHistPileUp;       // multiplicity distribution hist 

    
    ClassDef(Bin,4);  // Manager of data 
  };
  /** 
   * Add another eta bin to the task
   */
  AliBaseMultTask::Bin* MakeBin(Double_t etaLow, Double_t etaHigh);
protected:
  /**
   * Copy Constructor
   */ 
  AliForwardMultiplicityDistribution(const AliForwardMultiplicityDistribution&);
  /**
   * Assignment Operator
   */
  AliForwardMultiplicityDistribution& 
  operator=(const AliForwardMultiplicityDistribution&);
  /** 
   * Check the event
   * 
   * @param fwd Forwarddata 
   * 
   * @return true on success
   */
  virtual Bool_t CheckEvent(const AliAODForwardMult& fwd);
  ClassDef(AliForwardMultiplicityDistribution, 4); 
};

#endif
// Local Variables:
//   mode: C++
// End:
