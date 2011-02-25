//
// Class used to handle the input from AODs and put it into histograms
// the Forward Flow tasks can run on
//
#ifndef ALIFORWARDFLOWBASE_H
#define ALIFORWARDFLOWBASE_H
#include "AliAODEvent.h"
#include "AliAODForwardMult.h"

/**
 * 
 * Class used to handle the input from AODs and put it into histograms
 * the Forward Flow tasks can run on.
 * 
 * @ingroup pwg2_forward_tasks
 */
class AliForwardFlowBase : public TNamed
{
public:
  /**
   * Constructor
   */
   AliForwardFlowBase();
  /*
   * Constructor
   *
   * @param fList list of histograms for flow analysis
   */
   AliForwardFlowBase(TList* fList);
  /*
   * Copy constructor
   *
   * @param o Object to copy from
   */
    AliForwardFlowBase(const AliForwardFlowBase& o) : TNamed(),
  	fList(o.fList) {}
 /** 
   * Assignment operator 
   * 
   * @param o Object to assign from 
   * 
   * @return Reference to this object 
   */
  AliForwardFlowBase& operator=(const AliForwardFlowBase&) { return *this; }
  /*
   * Check that AOD event meet trigger requirements
   */
   Bool_t  AODCheck(AliAODForwardMult* aodfm);
  /*
   * Loop over AliAODForwardMult object and fill flow histograms
   */
   Bool_t  LoopAODFMD(AliAODEvent* AODevent);
  /*
   * Loop over AliAODForwardCentral object and fill flow histograms
   */
   Bool_t  LoopAODSPD(AliAODEvent* AODevent);
  /*
   * Loop over AliAODForwardMult and Central object and fill flow histograms
   */
   Bool_t  LoopAODFMDandSPD(AliAODEvent* AODevent);
  /*
   * Loop over AliAODMCParticle branch object and fill flow histograms
   */
   Bool_t  LoopAODmc(AliAODEvent* AODevent);
  /*
   * Loop over AliAODForwardMult object and fill flow histograms from track refs
   */
   Bool_t  LoopAODtrrefHits(AliAODEvent* AODevent);
  /*
   * Loop over AliAODMCParticle branch object, add pt dep. flow and fill flow histograms
   */
   Bool_t  LoopMCaddptFlow(AliAODEvent* AODevent);
  /*
   * Loop over AliAODMCParticle branch object, add pit dep. flow and fill flow histograms
   */
   Bool_t  LoopMCaddpdgFlow(AliAODEvent* AODevent);
  /*
   * Loop over AliAODMCParticle branch object, add eta dep flow and fill flow histograms
   */
   Bool_t  LoopMCaddetaFlow(AliAODEvent* AODevent);

private:  
  TList* 	fList;  // List of flow histograms
  
  ClassDef(AliForwardFlowBase, 0); 
};
 
#endif
// Local Variables:
//   mode: C++ 
// End:
