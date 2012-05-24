//
// Calculate the flow in the forward regions using the Q cumulants method
//
#ifndef ALIFORWARDMCFLOWTASKQC_H
#define ALIFORWARDMCFLOWTASKQC_H
/**
 * @file   AliForwardMCFlowTaskQC.h
 * @author Alexander Hansen alexander.hansen@cern.ch
 * @date   Tue Feb 14 2012
 * 
 * @brief  
 * 
 * 
 * @ingroup pwglf_forward_flow
 */
#include "AliForwardFlowTaskQC.h"
#include "AliForwardFlowWeights.h"
#include <TH2D.h>
class TGraph;

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
 *
 */
class AliForwardMCFlowTaskQC : public AliForwardFlowTaskQC
{
public:
  /**
   * Constructor
   */
  AliForwardMCFlowTaskQC();
  /*
   * Constructor
   *
   * @param name Name of task
   */
  AliForwardMCFlowTaskQC(const char* name);
  /**
   * Destructor 
   */
  virtual ~AliForwardMCFlowTaskQC() {}
  /**
   * Check for centrality in AliAODForwardMult object, 
   * if present return true - also sets fCent value
   * can be used to get centrality from impact parameter
   *
   * @param aodfm AliAODForwardMultObject
   * 
   * @return Bool_t 
   */
  virtual Bool_t GetCentrality(const AliAODForwardMult* aodfm);
  /*
   * Set use parametrization from impact parameter for centrality
   *
   * @param use Use impact par
   */
  void SetUseImpactParameter(Bool_t use) { fUseImpactPar = use; }
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
  /*
   * Copy constructor
   *
   * @param o Object to copy from
   */
  AliForwardMCFlowTaskQC(const AliForwardMCFlowTaskQC& o);
  /** 
   * Assignment operator 
   * 
   * @return Reference to this object 
   */
  AliForwardMCFlowTaskQC& operator=(const AliForwardMCFlowTaskQC& o);
  /*
   * Initiate vertex bin objects
   */
  void InitVertexBins();
   /*
   * Initiate diagnostics histograms
   */
  void InitHists();
  /*
   * Analyze event
   */
  Bool_t Analyze();
  /*
   * Finalize analysis
   */
  void Finalize();
  /**
   * Loop over AliAODMCParticle branch object and fill d^2N/detadphi histograms
   * add flow if arguments are set
   * 
   * @return true on success
   */
  Bool_t LoopAODMC();
  /**
   * Get centrality form MC impact parameter
   */
  Double_t GetCentFromB() const;
  
  TList         fBinsFMDTR;         //  List with FMDTR VertexBin objects
  TList         fBinsSPDTR;         //  List with SPDTR VertexBin objects
  TList         fBinsMC;            //  List with MC VertexBin objects
  TH2D          fdNdedpMC;          //  d^2N/detadphi MC truth histogram
  AliForwardFlowWeights fWeights;   //  Flow after burner 
  TGraph*       fImpactParToCent;   //  Parametrization of b to centrality
  Bool_t        fUseImpactPar;      //  Flag to use impact parameter for cent
  TString       fAddFlow;           //  Add flow string
  Int_t         fAddType;           //  Add flow type #
  Int_t         fAddOrder;          //  Add flow order

  ClassDef(AliForwardMCFlowTaskQC, 2); // FMD MC analysis task 
};
 
#endif
// Local Variables:
//   mode: C++ 
// End:
