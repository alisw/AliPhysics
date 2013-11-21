//
// Class AliRsnCutEventUtils
//
// This cut implementation checks the quality of event primary vertex.
// It currently works only with ESD events (not AOD).
//
// authors: Martin Vala (martin.vala@cern.ch)
//          Alberto Pulvirenti (alberto.pulvirenti@ct.infn.it)
//

#include "AliRsnCutEventUtils.h"
#include "AliAnalysisUtils.h"

ClassImp(AliRsnCutEventUtils)

//_________________________________________________________________________________________________
AliRsnCutEventUtils::AliRsnCutEventUtils(const char *name, Bool_t rmFirstEvInChunck, Bool_t checkPileUppA2013) :
AliRsnCut(name, AliRsnCut::kEvent),
  fIsRmFirstEvInChunck(rmFirstEvInChunck),
  fCheckPileUppA2013(checkPileUppA2013),
  fUseMVPlpSelection(kFALSE),
  fMinPlpContribMV(5),
  fMinPlpContribSPD(5),
  fUseVertexSelection2013pA(kFALSE),
  fUtils(0x0)
{
  //
  // Main constructor.
  //This class is mainly used for pPb 2013 
  //If the rmFirstEvInChunck flag is kTRUE, it removed the first event in chunk
  //
  //If the checkPileUp flag is kTRUE, it removes the events from pile-up
  //- to be used for pA 2013, for other periods the rejection of pileup is implemented 
  //ad an additional cut in the primary vertex
  
}

//_________________________________________________________________________________________________
AliRsnCutEventUtils::AliRsnCutEventUtils(const AliRsnCutEventUtils &copy) :
  AliRsnCut(copy),
  fIsRmFirstEvInChunck(copy.fIsRmFirstEvInChunck),
  fCheckPileUppA2013(copy.fCheckPileUppA2013),
  fUseMVPlpSelection(copy.fUseMVPlpSelection),
  fMinPlpContribMV(copy.fMinPlpContribMV),
  fMinPlpContribSPD(copy.fMinPlpContribSPD),
  fUseVertexSelection2013pA(copy.fUseVertexSelection2013pA),
  fUtils(copy.fUtils)
{
  //
  // Copy constructor.
  //
}

//-------------------------------------------------------------------------------------------------
AliRsnCutEventUtils &AliRsnCutEventUtils::operator=(const AliRsnCutEventUtils &copy)
{
  //
  // Assignment operator.
  // Works like copy constructor.
  //
  AliRsnCut::operator=(copy);
  if (this == &copy)
    return *this;
  
  fIsRmFirstEvInChunck=copy.fIsRmFirstEvInChunck;
  fCheckPileUppA2013=copy.fCheckPileUppA2013;
  fUseMVPlpSelection=copy.fUseMVPlpSelection;
  fMinPlpContribMV=copy.fMinPlpContribMV;
  fMinPlpContribSPD=copy.fMinPlpContribSPD;
  fUseVertexSelection2013pA=copy.fUseVertexSelection2013pA;

  fUtils=copy.fUtils;
	
    return (*this);
}
//_________________________________________________________________________________________________
Bool_t AliRsnCutEventUtils::IsSelected(TObject *object)
{
//
// Cut checker
//
   // coherence check
   // which also fills data member objects
   if (!TargetOK(object)) return kFALSE;

   // retrieve event
   AliVEvent *vevt = dynamic_cast<AliVEvent *>(fEvent->GetRef());
   fUtils = new AliAnalysisUtils();
   fUtils->SetUseMVPlpSelection(fUseMVPlpSelection);
   fUtils->SetMinPlpContribMV(fMinPlpContribMV);
   fUtils->SetMinPlpContribSPD(fMinPlpContribSPD);

   // pile-up check
   if ((fCheckPileUppA2013) && (fUtils->IsPileUpEvent(vevt)))
     return kFALSE;
  
   //remove first event in chunk 
   if ((fIsRmFirstEvInChunck) && (fUtils->IsFirstEventInChunk(vevt)))
     return kFALSE;
   
   //apply vertex selection - for 2013 pPb data
   if((fUseVertexSelection2013pA) && (!fUtils->IsVertexSelected2013pA(vevt)))
      return kFALSE;
 
   return kTRUE;
}
