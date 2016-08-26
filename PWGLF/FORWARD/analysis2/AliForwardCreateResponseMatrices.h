/**
 * @file   AliForwardCreateResponseMatrices.h
 * @author Valentina Zaccolo <Valentina.Zaccolo@cern.ch>
 * @date   Thu Feb  7 00:56:02 2013
 * 
 * @brief Task to make the reponse matrices used by the multiplicity
 * distibution analysis
 * 
 * @ingroup pwglf_forward_multdist
 */
/** 
 * @defgroup pwglf_forward_multdist Multiplicity Distributions
 * 
 * Code to do with @f$P(N_{ch})@f$ analysis - based on Valentina's code
 *
 * @ingroup pwglf_forward_topical
 */
#ifndef ALIFORWARDCREATERESPONSEMATRICES_H
#define ALIFORWARDCREATERESPONSEMATRICES_H
#include "AliBaseMultTask.h"

/**
 * Task to make the reponse matrices used by the multiplicity
 * distibution analysis
 * 
 * @ingroup pwglf_forward Tasks
 * @ingroup pwglf_forward_multdist
 */
class AliForwardCreateResponseMatrices : public AliBaseMultTask
{
public:
  /** 
   * 
   * Default Constructor 
   * 
   */
  AliForwardCreateResponseMatrices() : AliBaseMultTask() {}
  /**
   *  Constructor
   *
   */
  AliForwardCreateResponseMatrices(const char* name)
    : AliBaseMultTask(name)
  {}
  /** 
   * 
   * Destructor
   * 
   */
  virtual ~AliForwardCreateResponseMatrices(){}
  /**
   *  Embedded Class begins here
   */
  struct Bin : public AliBaseMultTask::Bin
  {
    /**
     * Default Constructor
     */
    Bin() : AliBaseMultTask::Bin(), fResponseMatrix(0)  {}
    /**
     * Constructor
     */
    Bin(Double_t etaLow, Double_t etaHigh)
      : AliBaseMultTask::Bin(etaLow, etaHigh),
	fResponseMatrix(0)
    {}
    /**
     * Copy Constructor
     */ 
    Bin(const Bin& o)
      : AliBaseMultTask::Bin(o),
	fResponseMatrix(0)
    {}
    /**
     * Assignment Operator
     */
    Bin& operator=(const Bin&) { return *this; }
    /**
     * Destructor
     */
    ~Bin(){}
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
    TH2D* fResponseMatrix;//! Response matrix (MC truth vs. ana multiplicity)
    ClassDef(Bin,3); // Manager of data 
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
  AliForwardCreateResponseMatrices(const AliForwardCreateResponseMatrices& o)
    : AliBaseMultTask(o)
  {}
  /**
   * Assignment operator
   *
   */
  AliForwardCreateResponseMatrices& 
  operator=(const AliForwardCreateResponseMatrices&) { return *this; }
  ClassDef(AliForwardCreateResponseMatrices, 5); 
};

#endif
// Local Variables:
//   mode: C++
// End:
