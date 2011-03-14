//
// Class used to handle the input from AODs and put it into histograms
// the Forward Flow tasks can run on
//
#ifndef ALIFORWARDFLOWUTIL_H
#define ALIFORWARDFLOWUTIL_H
#include "TNamed.h"
class AliAODForwardMult;
class AliAODEvent;
class TList;

/**
 * 
 * Class used to handle the input from AODs and put it into histograms
 * the Forward Flow tasks can run on.
 *
 * @ingroup pwg2_forward_tasks
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
   * @param fList list of histograms for flow analysis
   */
  AliForwardFlowUtil(TList* fList);
  /*
   * Copy constructor
   *
   * @param o Object to copy from
   */
  AliForwardFlowUtil(const AliForwardFlowUtil& o) : TNamed(),
						    fList(o.fList),
						    fZvertex(o.fZvertex) {}
  /** 
   * Assignment operator 
   * 
   * @param o Object to assign from 
   * 
   * @return Reference to this object 
   */
  AliForwardFlowUtil& operator=(const AliForwardFlowUtil&) { return *this; }
  /**
   * Check that AOD event meet trigger requirements
   */
  Bool_t AODCheck(const AliAODForwardMult* aodfm) const;
  /**
   * Loop over AliAODForwardMult object and fill flow histograms
   */
  Bool_t LoopAODFMD(const AliAODEvent* AODevent) const;
  /*
   * Loop over AliAODForwardCentral object and fill flow histograms
   */
  Bool_t LoopAODSPD(const AliAODEvent* AODevent) const;
  /**
   * Loop over AliAODForwardMult object and fill flow histograms from
   * track refs
   */
  Bool_t LoopAODtrrefHits(const AliAODEvent* AODevent) const;
 /**
   * Loop over AliAODMCParticle branch object and fill flow histograms
   * add flow if arguments are set
   */
  Bool_t LoopAODmc(const AliAODEvent* AODevent, TString addFlow, Int_t type, Int_t order) const;
 /**
   * Set Z vertex range - Used by flow task
   */
  void SetVertexRange(Int_t vertex = 2) { fZvertex = vertex; }
   
protected:
  /**
   * Add pt dependent flow factor
   */
  Double_t AddptFlow(Double_t Pt, Int_t type) const;
  /**
   * Add pid dependent flow factor
   */
  Double_t AddpidFlow(Int_t ID, Int_t type) const;
  /**
   * Add eta dependent flow factor
   */
  Double_t AddetaFlow(Double_t Eta, Int_t type) const;
  /**
   * Get centrality form MC impact parameter
   */
  Double_t GetCentFromMC(const AliAODEvent* AODevent) const;

  TList* 	fList;  // List of flow histograms
  Int_t		fZvertex; // Z vertex range

  ClassDef(AliForwardFlowUtil, 1); 
};
 
#endif
// Local Variables:
//   mode: C++ 
// End:
