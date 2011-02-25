//
// Calculate the flow in the forward regions using the Q cumulants method
//
#ifndef ALIFORWARDFLOWTASKQC_H
#define ALIFORWARDFLOWTASKQC_H
#include "AliAnalysisTaskSE.h"
#include "AliAODEvent.h"

/**
 * Calculate the flow in the forward regions using the Q cumulants method
 *
 * @par Inputs:
 *   - AliAODEvent
 *
 * Outputs:
 *   - AnalysisResults.root
 *
 * @ingroup pwg2_forward_tasks
 *
 */
class AliForwardFlowTaskQC : public AliAnalysisTaskSE
{
public:
  /** 
   * Constructor 
   */
   AliForwardFlowTaskQC();
  /** 
   * Constructor
   * 
   * @param name Name of task 
   */
   AliForwardFlowTaskQC(const char* name);
  /**
   * Destructor
   */
   virtual ~AliForwardFlowTaskQC() {;}
  /** 
   * Copy constructor 
   * 
   * @param o Object to copy from 
   */
   AliForwardFlowTaskQC(const AliForwardFlowTaskQC& o);
  /** 
   * Assignment operator 
   * 
   * @param o Object to assign from 
   * 
   * @return Reference to this object 
   */
   AliForwardFlowTaskQC& operator=(const AliForwardFlowTaskQC&) { return *this; }
   /** 
   * @{ 
   * @name Interface methods 
   */
  /** 
   * Create output objects 
   * 
   */
   virtual void CreateOutputObjects();
  /**
   * Initialize the task
   *
   */
   virtual void Init() {}
  /** 
   * Process each event 
   *
   * @param option Not used
   */  
   virtual void UserExec(Option_t *option);
  /**
   * Calculate Q cumulant
   * Parameters:
   *  type: Determines which histograms should be used
   *        - "" = data histograms
   *        - "TrRef" = track reference histograms
   *        - "MC" = MC truth histograms
   *  harmonic: Which harmonic to calculate
   */
   void CumulantsMethod(TString type, Int_t harmonic);
  /** 
   * End of job
   * 
   * @param option Not used 
   */
   virtual void Terminate(Option_t *option);
  /*
   * if MC information is available do analysis on Monte Carlo truth and track references
   *
   */
   void ProcessPrimary();
  /*
   * Returns the outputlist
   *
   */
   TList* GetOutputList() { return fOutputList; }
  /* 
   * Set Number of eta bins to be used in flow analysis
   *
   */
   void SetUseNEtaBins(Int_t nbins) { fEtaBins = nbins; }
  /* 
   * Set which harmonics to calculate
   * v_{1} to v_{4} is available and calculated as default
   *
   */
   void SetDoHarmonics(Bool_t v1 = kTRUE, Bool_t v2 = kTRUE, Bool_t v3 = kTRUE, Bool_t v4 = kTRUE)
    { fv[1] = v1; fv[2] = v2; fv[3] = v3; fv[4] = v4; }
  /** 
   * @} 
   */
   virtual void SetDebugLevel(Int_t level) { fDebug = level; }
  
private:
  /*
   * Caclulate the variance of x squared - used to finalize calculations in Terminate()
   *
   */
   Double_t VarSQ(Double_t wxx2, Double_t x, Double_t wx, Double_t wxx, Double_t sqrtwx2);
  /*
   * Caclulate the covariance between x and y - used to finalize calculations in Terminate()
   *
   */
   Double_t CovXY(Double_t wxwyxy, Double_t wxwy, Double_t XY, Double_t wx, Double_t wy);

   Int_t          fDebug;        //  Debug flag
   TList*         fOutputList;   //  Output list
   AliAODEvent*   fAOD;          //  AOD event
   Bool_t         fMC;           //  Is MC flags
   Int_t          fEtaBins;      //  Number of eta bins in histograms
   Bool_t         fv[5];         //  Calculate v_{n} flag
   
   ClassDef(AliForwardFlowTaskQC, 1); // Analysis task for FMD analysis
};
 
#endif
// Local Variables:
//   mode: C++ 
// End:
