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
   * @{ 
   * @name Task interface methods 
   */
  /**
   * Loop over AliAODMCParticle branch object and fill d^2N/detadphi histograms
   * add flow if arguments are set
   * 
   * @return true on success
   */
  Bool_t LoopAODMC();
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
   * Add pt dependent flow factor
   *
   * @param Pt   @f$ p_T@f$
   * @param type Type of flow 
   */
  Double_t AddptFlow(Double_t pt) const;
  /**
   * Add pid dependent flow factor
   *
   * @param ID   Particle ID 
   * @param type Type of flow
   */
  Double_t AddpidFlow(Int_t id) const;
  /**
   * Add eta dependent flow factor
   * 
   * @param Eta  @f$\eta@f$ 
   * @param type Type of flow 
   */
  Double_t AddetaFlow(Double_t eta) const;
  /**
   * Get centrality form MC impact parameter
   */
  Double_t GetCentFromB() const;
  
  TList         fBinsFMDTR;         //  List with FMDTR VertexBin objects
  TList         fBinsSPDTR;         //  List with SPDTR VertexBin objects
  TList         fBinsMC;            //  List with MC VertexBin objects
  TH2D          fdNdedpMC;          //  d^2N/detadphi MC truth histogram
  TGraph*       fAliceCent4th;      //  Parametrization of ALICE QC4 vs. cent. data
  TGraph*       fAlicePt2nd4050;    //  Parametrization of ALICE QC2 vs. pT data
  TGraph*       fAlicePt4th3040;    //  Parametrization of ALICE QC4 vs. pT data
  TGraph*       fAlicePt4th4050;    //  Parametrization of ALICE QC4 vs. pT data
  TGraph*       fImpactParToCent;   //  Parametrization of b to centrality datapoints
  TString       fAddFlow;           //  Add flow string
  Int_t         fAddType;           //  Add flow type #
  Int_t         fAddOrder;          //  Add flow order

  ClassDef(AliForwardMCFlowTaskQC, 1); // FMD MC analysis task 
};
 
#endif
// Local Variables:
//   mode: C++ 
// End:
