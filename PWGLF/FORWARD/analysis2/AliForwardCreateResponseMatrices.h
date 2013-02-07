/**
 * @file   AliForwardCreateResponseMatrices.h
 * @author Christian Holm Christensen <cholm@master.hehi.nbi.dk>
 * @date   Thu Feb  7 00:56:02 2013
 * 
 * @brief  
 * 
 * @ingroup pwglf_forward_multdist
 */
/** 
 * @defgroup pwglf_forward_multdist Multiplicity Distributions
 * 
 * Code to do with @f$P(N_{ch})@f$ analysis
 *
 * @ingroup pwglf_forward_topical
 */
#ifndef ALIFORWARDCREATERESPONSEMATRICES_H
#define ALIFORWARDCREATERESPONSEMATRICES_H
  
#include "AliAnalysisTaskSE.h"
#include "AliBasedNdetaTask.h"
#include <TList.h>

class TH2D;

/**
 * Task to make the reponse matrices used by the multiplicity
 * distibution analysis
 * 
 * @ingroup pwglf_forward Tasks
 * @ingroup pwglf_forward_multdist
 * @todo Should not inherit from AliBasedNdetaTask 
 */
class AliForwardCreateResponseMatrices : public AliBasedNdetaTask
{
public:
  /** 
   * 
   * Default Constructor 
   * 
   */
  AliForwardCreateResponseMatrices();
  /**
   *  Constructor
   *
   */
  AliForwardCreateResponseMatrices(const char* name);
  /**
   * Copy Constructor
   *
   */
  AliForwardCreateResponseMatrices(const AliForwardCreateResponseMatrices& o) : AliBasedNdetaTask(o), fTrigger(),fBins(), fOutput() { }
  /**
   * Assignment operator
   *
   */
  AliForwardCreateResponseMatrices& operator=(const AliForwardCreateResponseMatrices&){return *this;}
  /** 
   * 
   * Destructor
   * 
   */
  virtual ~AliForwardCreateResponseMatrices(){}
  /**
   *  Embedded Class begins here
   */
  struct Bin : public TNamed
  {
    /**
     * Default Constructor
     */
    Bin();
    /**
     * Constructor
     */
    Bin(Double_t etaLow, Double_t etaHigh);
    /**
     * Copy Constructor
     */ 
    Bin(const Bin&);
    /**
     * Assignment Operator
     */
    Bin&operator=(const Bin&){return*this;}
    /**
     * Destructor
     */
    ~Bin(){}
    /**
     * Define outputs of a single eta bin
     */
    virtual void CreateOutputObjectss(TList* cont,  Int_t max);
    /**
     * Process a single eta bin
     */    
    virtual void Process(TH1D* dndetaForward, TH1D* dndetaCentral, TH1D* normForward,   TH1D* normCentral, TH1D* dndetaMC, Double_t VtxZ, Bool_t selectedTrigger,  Bool_t isMCNSDm, Bool_t isESDNSD, AliAODEvent* aodevent);
    Double_t fEtaLow;                  // low eta limit 
    Double_t fEtaHigh;                 // high eta limit 
    TH1D*    fHist;                    // multiplicity histogram 
    TH1D*    fHistMC;                  // multiplicity histogram MC truth primaries
    TH2D*    fAcceptance;              // histogram showing the 'holes' in acceptance. BinContent of 1 shows a hole, and BinContent of 10 shows data coverage
    TH2D*    fVtxZvsNdataBins;         // VtxZ vs. number of data acceptance bins (normalised to the eta range) 
    TH2D*    fResponseMatrix;          //Response matrix (MC truth vs. analysed multiplicity)
    TH2D*    fResponseMatrixPlus05;    //Response matrix with analysed multiplicity scaled up by 5%
    TH2D*    fResponseMatrixPlus075;   //Response matrix  with analysed multiplicity scaled up by 7.5%
    TH2D*    fResponseMatrixPlus10;    //Response matrix with analysed multiplicity scaled up by 10%
    TH2D*    fResponseMatrixMinus05;   //Response matrix with analysed multiplicity scaled down by 5%
    TH2D*    fResponseMatrixMinus075;  //Response matrix with analysed multiplicity scaled down by 7.55%
    TH2D*    fResponseMatrixMinus10;   //Response matrix with analysed multiplicity scaled down by 10%
    TH2D*    fResponseMatrixMinusSys;  //Response matrix with analysed multiplicity scaled up by event mult uncertainty
    TH2D*    fResponseMatrixPlusSys;   //Response matrix with analysed multiplicity scaled down by event mult uncertainty
    TH1D*    fESDNSD;                  //number of events found as NSD by the analysis vs. multiplicity
    TH1D*    fMCNSD;                   //number of events found as NSD by the MC truth vs. multiplicity
    TH1D*    fMCESDNSD;                //number of events found as NSD by both analysis and MC truth vs. multiplicity
    TH1D*    fTriggerBias;             // histogram for trigger vertex bias correction
   ClassDef(Bin,2); // Manager of data 
  };
  /**
   * Create Output Objects
   */
  virtual void UserCreateOutputObjects();
  /**
   * User Exec
   */
  void UserExec(Option_t *option);
  /**
   * Terminate
   */
  void Terminate(Option_t *option);
  /** 
   * implementation of pure virtual function, always returning 0
   */
  virtual TH2D* GetHistogram(const AliAODEvent* aod, Bool_t mc);
  /** 
   * Get single event forward and central @f$d^2N/d\eta d\phi@f$
   * histograms
   * 
   */  
  virtual void GetHistograms(const AliAODEvent* aod, TH2D& forward, TH2D& central); 
  /** 
   * Add another eta bin to the task
   */
  void AddBin(Double_t etaLow, Double_t etaHigh){fBins.Add(new Bin(etaLow, etaHigh)); }
 protected:
  TH1I* fTrigger;  //Trigger histogram
  TList  fBins;    // List of eta bins
  TList*  fOutput; // Output list
  ClassDef(AliForwardCreateResponseMatrices, 2); 
};

#endif
// Local Variables:
//   mode: C++
// End:
