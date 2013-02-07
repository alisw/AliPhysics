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
#include "AliAnalysisTaskSE.h"
#include "AliBasedNdetaTask.h"
#include <TList.h>

class TH2D;

/**
 * Task to do the multiplicity distibution
 * 
 * @ingroup pwglf_forward Tasks
 * @ingroup pwglf_forward_multdist
 * @todo Should not inherit from AliBasedNdetaTask 
 */
class AliForwardMultiplicityDistribution : public AliBasedNdetaTask
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
   * Copy Constructor
   */ 
  AliForwardMultiplicityDistribution(const AliForwardMultiplicityDistribution& o) : AliBasedNdetaTask(o), fTrigger(o.fTrigger),fBins(), fOutput(o.fOutput), fLowCent(o.fLowCent), fHighCent(o.fHighCent),fNBins(o.fNBins), fCent(o.fCent){ }
  /**
   * Assignment Operator
   */
  AliForwardMultiplicityDistribution& operator=(const AliForwardMultiplicityDistribution&){return *this;}
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
    Bin&operator=(const Bin&){return*this;}
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
    TH1D*    fHistPlus05;       // multiplicity distribution hist scaled up with 5%
    TH1D*    fHistPlus075;      // multiplicity distribution hist scaled up with 7.5%
    TH1D*    fHistPlus10;       // multiplicity distribution hist scaled up with 10%
    TH1D*    fHistMinus05;      // multiplicity distribution hist scaled down with 5%
    TH1D*    fHistMinus075;     // multiplicity distribution hist scaled down with 7.5%
    TH1D*    fHistMinus10;      // multiplicity distribution hist scaled down with 10%
    TH1D*    fHistPlusSys;      // multiplicity distribution hist scaled up with the event uncertainty
    TH1D*    fHistMinusSys;     // multiplicity distribution hist scaled down with the event uncertainty
    TH2D*    fAcceptance;       // histogram showing the 'holes' in acceptance. 
                                // BinContent of 1 shows a hole, and BinContent of 10 shows data coverage
    TH2D*    fVtxZvsNdataBins;  // VtxZ vs. number of data acceptance bins (normalised to the eta range) 
    
    ClassDef(Bin,2);  // Manager of data 
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
   * Set Centrality
   */
  void SetCentrality(Int_t lowCent, Int_t highCent){fLowCent= lowCent; fHighCent= highCent;}
  /**
   * Set fNBins, multiplicity histos run from 0 to fNBins
   */
  void SetNBins(Int_t n){fNBins= n;}
  /** 
   * implementation of pure virtual function, always returning 0
   */
  virtual TH2D* GetHistogram(const AliAODEvent* aod, Bool_t mc);
  /**
   * Get single event forward and central dN²/dEta dPhi histograms 
   */
  virtual void GetHistograms(const AliAODEvent* aod, TH2D& forward, TH2D& central , Bool_t mc=false);
  /** 
   * Add another eta bin to the task
   */
  void AddBin(Double_t etaLow, Double_t etaHigh){fBins.Add(new Bin(etaLow, etaHigh)); }
  /**
   *  Form name of eta bin
   */
  static const Char_t* FormBinName(Double_t etaLow, Double_t etaHigh);
protected:
  TH1I*  fTrigger;   // trigger histogram
  TList  fBins;      // eta bin list
  TList* fOutput;    // output list
  Int_t  fLowCent;   // lower centrality limit
  Int_t  fHighCent;  // upper centrality limit
  Int_t  fNBins;     // multiplicity axis' runs from 0 to fNbins
  TH1D*  fCent;      // centrality
  ClassDef(AliForwardMultiplicityDistribution, 2); 
};

#endif
// Local Variables:
//   mode: C++
// End:
