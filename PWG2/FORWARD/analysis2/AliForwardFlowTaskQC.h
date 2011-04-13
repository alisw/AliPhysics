//
// Calculate the flow in the forward regions using the Q cumulants method
//
#ifndef ALIFORWARDFLOWTASKQC_H
#define ALIFORWARDFLOWTASKQC_H
/**
 * @file   AliForwardFlowTaskQC.h
 * @author Alexander Hansen alexander.hansen@cern.ch
 * @date   Fri Mar 25 13:53:00 2011
 * 
 * @brief  
 * 
 * 
 * @ingroup pwg2_forward_flow
 */
#include "AliAnalysisTaskSE.h"
#include "AliForwardFlowUtil.h"
class AliAODEvent;

 /**
 * @defgroup pwg2_forward_tasks_flow Flow tasks 
 * @ingroup pwg2_forward_tasks
 */
/**
 * Calculate the flow in the forward regions using the Q cumulants method
 *
 * @par Inputs:
 *   - AliAODEvent
 *
 * Outputs:
 *   - AnalysisResults.root
 *
 * @ingroup pwg2_forward_tasks_flow
 * @ingroup pwg2_forward_flow
 *
 * @todo Add centrality stuff
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
  virtual ~AliForwardFlowTaskQC() {}
  /** 
   * @{ 
   * @name Task interface methods 
   */
  /** 
   * Create output objects 
   * 
   */
  virtual void UserCreateOutputObjects();
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
   * End of job
   * 
   * @param option Not used 
   */
  virtual void Terminate(Option_t *option);
  /* @} */
  /*
   * Returns the outputlist
   *
   */
  TList* GetOutputList() { return fOutputList; }
  /* 
   * Set Number of @f$ \eta@f$ bins to be used in flow analysis
   *
   */
  void SetUseNEtaBins(Int_t nbins) { fEtaBins = nbins; }
  /*
   * Set which harmonics to calculate. @f$ v_{1}@f$ to @f$ v_{4}@f$ is
   * available and calculated as default
   *
   * @param v1  Do @f$ v_{1}@f$ 
   * @param v2  Do @f$ v_{2}@f$ 
   * @param v3  Do @f$ v_{3}@f$ 
   * @param v4  Do @f$ v_{4}@f$ 
   */
  void SetDoHarmonics(Bool_t v1 = kTRUE, Bool_t v2 = kTRUE, 
		      Bool_t v3 = kTRUE, Bool_t v4 = kTRUE) { 
    fv[1] = v1; fv[2] = v2; fv[3] = v3; fv[4] = v4; }
  /*
   * Set string to add flow to MC truth particles
   *
   * @param type String
   */
  void AddFlow(TString type = "") { fAddFlow = type; }
  /*
   * Set which function fAddFlow should use
   *
   * @param type of AddFlow 
   */
  void AddFlowType(Int_t number = 0) { fAddType = number; }
  /*
   * Set which order of flow to add
   *
   * @param order Flow order 
   */
  void AddFlowOrder(Int_t order = 2) { fAddOrder = order; }
protected:
  /** 
   * Copy constructor 
   * 
   * @param o Object to copy from 
   */
  AliForwardFlowTaskQC(const AliForwardFlowTaskQC& o);
  /** 
   * Assignment operator 
   * 
   * @return Reference to this object 
   */
  AliForwardFlowTaskQC& operator=(const AliForwardFlowTaskQC&) { return *this; }
  /*
   * if MC information is available do analysis on Monte Carlo truth
   * and track references
   *
   */
  void ProcessPrimary();
  /**
   * Calculate Q cumulant
   * 
   * @param type     Determines which histograms should be used
   *                 - "FMD" = FMD data histograms
   *                 - "SPD" = SPD data histograms
   *                 - "TrRef" = track reference histograms
   *                 - "MC" = MC truth histograms
   * @param harmonic Which harmonic to calculate
   */
  void CumulantsMethod(TString type, Int_t harmonic);
  /**
   * Caclulate the variance of x squared - used to finalize
   * calculations in Terminate()
   *
   * @param wxx2      Weight of @f$ x^2@f$
   * @param x         @f$ x@f$
   * @param wx        Weight of @f$ x@f$
   * @param wxx       Weight of @f$ x^2@f$
   * @param sqrtwx2   @f$ \sqrt{wx^2}@f$ 
   * 
   * @return Variance squared 
   */
  Double_t VarSQ(Double_t wxx2, Double_t x, Double_t wx, 
		 Double_t wxx, Double_t sqrtwx2) const ;
  /**
   * Caclulate the covariance between x and y - used to finalize
   * calculations in Terminate()
   *
   * @param wxwyxy @f$ w_x w_y x y@f$
   * @param wxwy   @f$ w_x w_y@f$
   * @param XY     @f$ xy@f$ 
   * @param wx     @f$ w_x@f$
   * @param wy     @f$ w_y@f$
   * 
   * @return 
   */
  Double_t CovXY(Double_t wxwyxy, Double_t wxwy, Double_t XY, 
		 Double_t wx, Double_t wy) const;

  TList*         fOutputList;   //  Output list
  AliForwardFlowUtil* fFlowUtil;//  AliForwardFlowUtil
  AliAODEvent*   fAOD;          //  AOD event
  Bool_t         fMC;           //  Is MC flags
  Int_t          fEtaBins;      //  Number of eta bins in histograms
  Bool_t         fv[5];         //  Calculate v_{n} flag
  TString        fAddFlow;	//  Add flow string
  Int_t          fAddType;	//  Add flow type #
  Int_t          fAddOrder;	//  Add flow order
  Float_t  	 fZvertex;	//  Z vertex bin
  Double_t       fCent;         //  Centrality
  

  ClassDef(AliForwardFlowTaskQC, 2); // Analysis task for FMD analysis
};
 
#endif
// Local Variables:
//   mode: C++ 
// End:
