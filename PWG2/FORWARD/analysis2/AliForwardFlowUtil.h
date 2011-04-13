//
// Class used to handle the input from AODs and put it into histograms
// the Forward Flow tasks can run on
//
#ifndef ALIFORWARDFLOWUTIL_H
#define ALIFORWARDFLOWUTIL_H
/**
 * @file   AliForwardFlowUtil.h
 * @author Alexander Hansen alexander.hansen@cern.ch
 * @date   Fri Mar 25 13:15:40 2011
 * 
 * @brief  
 * 
 * 
 * @ingroup pwg2_forward_flow
 */
#include "TNamed.h"
class AliAODForwardMult;
class AliAODEvent;
class TList;

/**
 * 
 * Class used to handle the input from AODs and put it into histograms
 * the Forward Flow tasks can run on.
 *
 * @ingroup pwg2_forward_tasks_flow
 * @ingroup pwg2_forward_flow
 */
class AliForwardFlowUtil : public TNamed
{
public:
  /**
   * Constructor
   */
  AliForwardFlowUtil();
  /*
   * Constructor
   *
   * @param l list of histograms for flow analysis
   */
  AliForwardFlowUtil(TList* l);
  /**
   * Check that AOD event meet trigger requirements
   * 
   * @param aodfm Forward multplicity AOD event structure 
   * 
   * @return true on success 
   */
  Bool_t AODCheck(const AliAODForwardMult* aodfm) const;
  /**
   * Loop over AliAODForwardMult object and fill flow histograms
   * 
   * @param AODevent AOD event structure 
   * 
   * @return true on success
   */
  Bool_t LoopAODFMD(const AliAODEvent* AODevent);
  /*
   * Loop over AliAODCentralMult object and fill flow histograms
   * 
   * @param AODevent AOD event structure 
   * 
   * @return true on success
   */
  Bool_t LoopAODSPD(const AliAODEvent* AODevent) const;
  /**
   * Loop over AliAODForwardMult object and fill flow histograms from
   * track refs
   * 
   * @param AODevent AOD event structure 
   * 
   * @return true on success
   */
  Bool_t LoopAODFMDtrrefHits(const AliAODEvent* AODevent) const;
  /**
   * Loop over AliAODCentralMult object and fill flow histograms from
   * track refs
   * 
   * @param AODevent AOD event structure 
   * 
   * @return true on success
   */
  Bool_t LoopAODSPDtrrefHits(const AliAODEvent* AODevent) const;
  /**
   * Loop over AliAODMCParticle branch object and fill flow histograms
   * add flow if arguments are set
   * 
   * @param AODevent AOD event structure 
   * @param addFlow  What to add flow to 
   * @param type     Type of flow 
   * @param order    Order of added flow 
   * 
   * @return true on success
   */
  Bool_t LoopAODmc(const AliAODEvent* AODevent, TString addFlow, 
                  Int_t type, Int_t order) const;
 /**
   * Get centrality from newest processed event
   */
  Double_t GetCentrality() const { return fCent; }
 /**
   * Get z vertex coordinate from newest processed event
   */
  Float_t GetVertex() const { return fVertex; }
   
protected:
  /*
   * Copy constructor
   *
   * @param o Object to copy from
   */
  AliForwardFlowUtil(const AliForwardFlowUtil& o) : TNamed(),
						    fList(o.fList),
						    fCent(o.fCent),
						    fVertex(o.fVertex) {}
  /** 
   * Assignment operator 
   * 
   * @return Reference to this object 
   */
  AliForwardFlowUtil& operator=(const AliForwardFlowUtil&) { return *this; }
  /**
   * Add pt dependent flow factor
   *
   * @param Pt   @f$ p_T@f$
   * @param type Type of flow 
   */
  Double_t AddptFlow(Double_t Pt, Int_t type) const;
  /**
   * Add pid dependent flow factor
   *
   * @param ID   Particle ID 
   * @param type Type of flow
   */
  Double_t AddpidFlow(Int_t ID, Int_t type) const;
  /**
   * Add eta dependent flow factor
   * 
   * @param Eta  @f$\eta@f$ 
   * @param type Type of flow 
   */
  Double_t AddetaFlow(Double_t Eta, Int_t type) const;
  /**
   * Get centrality form MC impact parameter
   * 
   * @param AODevent AOD event structure 
   */
  Double_t GetCentFromMC(const AliAODEvent* AODevent) const;

  TList* 	fList;   // List of flow histograms
  Double_t	fCent;   // centrality
  Float_t         fVertex; // z vertex coordinate

  ClassDef(AliForwardFlowUtil, 1); 
};
 
#endif
// Local Variables:
//   mode: C++ 
// End:
