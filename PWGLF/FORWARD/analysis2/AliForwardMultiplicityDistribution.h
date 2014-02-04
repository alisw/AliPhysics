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
#include "AliBaseAODTask.h"

class TH2D;

/**
 * Task to do the multiplicity distibution
 * 
 * @ingroup pwglf_forward Tasks
 * @ingroup pwglf_forward_multdist
 * @todo Should not inherit from AliBasedNdetaTask 
 */
class AliForwardMultiplicityDistribution : public AliBaseAODTask
{
public:
  /**
   * Default Constructor
   */
  AliForwardMultiplicityDistribution();
  /**
   * Constructor
   */
  AliForwardMultiplicityDistribution(const char* name);
  /**
   * Destructor
   */
  virtual ~AliForwardMultiplicityDistribution(){}
  /**
   * Embedded Class begins here
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
     * Assignment operator
     */
    Bin& operator=(const Bin&){return*this;}
    /**
     * Destructor
     */    
    ~Bin(){}
    /**
     *  Define outputs of a single eta bin
     */
    virtual void CreateOutputObjectss(TList* cont,  Int_t max);
    /**
     * Process a single eta bin
     */
    virtual void Process(TH1D* dndetaForward, TH1D* dndetaCentral, TH1D* normForward,   TH1D* normCentral, Double_t VtxZ);
    Double_t fEtaLow;           // low eta limit 
    Double_t fEtaHigh;          // high eta limit 
    TH1D*    fHist;             // multiplicity distribution hist 
    TH1D*    fHistPlus05;       // mult. dist. hist scaled up with 5%
    TH1D*    fHistPlus075;      // mult. dist. hist scaled up with 7.5%
    TH1D*    fHistPlus10;       // mult. dist. hist scaled up with 10%
    TH1D*    fHistMinus05;      // mult. dist. hist scaled down with 5%
    TH1D*    fHistMinus075;     // mult. dist. hist scaled down with 7.5%
    TH1D*    fHistMinus10;      // mult. dist. hist scaled down with 10%
    TH1D*    fHistPlusSys;      // mult. dist. hist scaled up with the event uncertainty
    TH1D*    fHistMinusSys;     // mult. dist, hist scaled down with the event uncertainty
    TH2D*    fAcceptance;       // histogram showing the 'holes' in acceptance. 
                                // BinContent of 1 shows a hole, and BinContent of 10 shows data coverage
    TH2D*    fVtxZvsNdataBins;  // VtxZ vs. number of data acceptance bins (normalised to the eta range) 
    
    ClassDef(Bin,3);  // Manager of data 
  };
  /**
   * Create Output Objects
   */
  virtual Bool_t Book();
  /**
   * User Exec
   */
  Bool_t Event(AliAODEvent& aod);
  /**
   * Terminate
   */
  Bool_t Finalize() { return true; }
  /**
   * Set Centrality
   */
  void SetCentrality(Int_t low, Int_t high) { SetCentralityAxis(low, high); }
  /**
   * Set fNBins, multiplicity histos run from 0 to fNBins
   */
  void SetNBins(Int_t n){fNBins= n;}
  /** 
   * Add another eta bin to the task
   */
  void AddBin(Double_t etaLow, Double_t etaHigh){fBins.Add(new Bin(etaLow, etaHigh)); }
  /**
   *  Form name of eta bin
   */
  static const Char_t* FormBinName(Double_t etaLow, Double_t etaHigh);
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

  TList  fBins;      // eta bin list
  Int_t  fNBins;     // multiplicity axis' runs from 0 to fNbins
  ClassDef(AliForwardMultiplicityDistribution, 3); 
};

#endif
// Local Variables:
//   mode: C++
// End:
