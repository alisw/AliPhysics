//
// Calculate the flow in the forward regions using the Q cumulants method
//
#ifndef ALIFORWARDMCFLOWTASKQC_H
#define ALIFORWARDMCFLOWTASKQC_H
/**
 * @file   AliForwardMCFlowTaskQC.h
 * @author Alexander Hansen alexander.hansen@cern.ch
 * 
 * @brief  
 * 
 * @ingroup pwglf_forward_flow
 */
#include "AliForwardFlowTaskQC.h"
#include "AliForwardFlowWeights.h"
class TGraph;
class AliAODMCHeader;

/**
 * Calculate the flow in the forward regions using the Q cumulants method
 *
 * @par Inputs:
 *   - AliAODEvent
 *
 * Outputs:
 *   - forward_flow.root
 *
 * @ingroup pwglf_forward_tasks_flow
 * @ingroup pwglf_forward_flow
 *
 */
class AliForwardMCFlowTaskQC : public AliForwardFlowTaskQC
{
public:
  /**
   * Constructor
   */
  AliForwardMCFlowTaskQC();
  /**
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
   * Set use parametrization from impact parameter for centrality
   *
   * @param use Use impact par
   */
  void SetUseImpactParameter(Bool_t use = kTRUE) { fUseImpactPar = use; }
  /**
   * Set to get vertex from MC header
   *
   * @param use Get from MC header
   */
  void SetUseMCHeaderVertex(Bool_t use = kTRUE) { fUseMCVertex = use; }
  /**
   * Add flow to MC particles
   */
  void SetUseFlowWeights(Bool_t use = kTRUE) { fUseFlowWeights = use; }
   
protected:
  /**
   * Copy constructor
   *
   * @param o Object to copy from
   */
  AliForwardMCFlowTaskQC(const AliForwardMCFlowTaskQC& o);
  /** 
   * Assignment operator 
   * 
   * @param o Object to assing from 
   *
   * @return Reference to this object 
   */
  AliForwardMCFlowTaskQC& operator=(const AliForwardMCFlowTaskQC& o);
  /**
   * Initiate vertex bin objects
   */
  void InitVertexBins();
  /**
   * Initiate diagnostics histograms
   */
  void InitHists();
  /**
   * Analyze event
   *
   * @return true on success 
   */
  Bool_t Analyze();
  /**
   * Finalize analysis
   */
  void Finalize();
    /**
   * Check AODevent object for trigger, vertex and centrality
   * uses aod header if object is null
   * returns true if event is OK
   *
   * @param aodfm AliAODForwardMult object
   * 
   * @return Bool_t 
   */
  virtual Bool_t CheckEvent(const AliAODForwardMult* aodfm);
  /**
   * Check trigger from AODForwardMult object
   * returns true if B trigger is present
   *
   * @param aodfm AliAODForwardMultObject
   * 
   * @return Bool_t 
   */
  virtual Bool_t CheckTrigger(const AliAODForwardMult* aodfm) const;
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
  /**
   * Check for vertex in MCHeader
   * returns true if in range of fVtxAXis, also sets fVtx value
   *
   * @param aodfm Not used
   * 
   * @return Bool_t 
   */
  virtual Bool_t GetVertex(const AliAODForwardMult* aodfm);
  /**
   * Loop over AliAODMCParticle branch object and fill d^2N/detadphi histograms
   * add flow if arguments are set
   * 
   * @return true on success
   */
  Bool_t FillMCHist();
  /**
   * Get centrality form MC impact parameter
   *
   * @return Centrality
   */
  Double_t GetCentFromB() const;
  
  TList                 fBinsForwardTR;   //  List with FMDTR VertexBin objects
  TList                 fBinsCentralTR;   //  List with SPDTR VertexBin objects
  TList                 fBinsMC;          //  List with MC VertexBin objects
  AliAODMCHeader*       fAODMCHeader;     //  MC header object
  TH2D                  fHistdNdedpMC;    //  d^2N/detadphi MC particles histogram
  TH2D*                 fHistFMDMCCorr;   //  Diagnostics for mult. corr. between FMD and MC
  TH2D*                 fHistSPDMCCorr;   //  Diagnostics for mult. corr. between SPD and MC
  AliForwardFlowWeights* fWeights;         //  Flow after burner 
  TGraph*               fImpactParToCent; //  Parametrization of b to centrality
  Bool_t                fUseImpactPar;    //  Flag to use impact parameter for cent
  Bool_t                fUseMCVertex;     //  Get vertex from MC header
  Bool_t                fUseFlowWeights;  //  Add flow

  ClassDef(AliForwardMCFlowTaskQC, 6); // FMD MC analysis task 
};
 
#endif
// Local Variables:
//   mode: C++ 
// End:
