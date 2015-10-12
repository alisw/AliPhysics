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
#include <AliHeader.h>
#include <AliAODMCHeader.h>
#include <AliGenDPMjetEventHeader.h>

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
  fMaxVtxZ(10.0),
  fFilterNSDeventsDPMJETpA2013(kFALSE),
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
  fMaxVtxZ(copy.fMaxVtxZ),
  fFilterNSDeventsDPMJETpA2013(copy.fFilterNSDeventsDPMJETpA2013),
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
  fMaxVtxZ=copy.fMaxVtxZ;
  fFilterNSDeventsDPMJETpA2013=copy.fFilterNSDeventsDPMJETpA2013;
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
   //set cuts in analysis utils
   fUtils = new AliAnalysisUtils();
   fUtils->SetUseMVPlpSelection(fUseMVPlpSelection);
   fUtils->SetMinPlpContribMV(fMinPlpContribMV);
   fUtils->SetMinPlpContribSPD(fMinPlpContribSPD);
   fUtils->SetMaxVtxZ(fMaxVtxZ);

   // pile-up check
   if ((fCheckPileUppA2013) && (fUtils->IsPileUpEvent(vevt)))
     return kFALSE;
  
   //remove first event in chunk 
   if ((fIsRmFirstEvInChunck) && (fUtils->IsFirstEventInChunk(vevt)))
     return kFALSE;
   
   //apply vertex selection - for 2013 pPb data
   if((fUseVertexSelection2013pA) && (!fUtils->IsVertexSelected2013pA(vevt)))
      return kFALSE;

   //apply filter for NSD events in DPMJET MC for pA
   if (fFilterNSDeventsDPMJETpA2013){
     //Adaptation of snippet received from David D. Chinellato on 07/09/2015 +
     // this: http://svnweb.cern.ch/world/wsvn/AliRoot/trunk/PWGLF/SPECTRA/ChargedHadrons/dNdPt/AlidNdPtHelper.cxx?rev=61295
     AliGenDPMjetEventHeader* dpmHeader;
     
     if (fEvent->IsAOD()){
       AliAODEvent *aodEvt = (AliAODEvent*)(fEvent->GetRef());
       AliAODMCHeader *aodMCheader = (AliAODMCHeader*)aodEvt->GetList()->FindObject(AliAODMCHeader::StdBranchName());
       if (!aodMCheader) {
	 AliError("MC header branch in AOD not found!\n");
	 return kFALSE;
       }
       dpmHeader = dynamic_cast<AliGenDPMjetEventHeader*>(aodMCheader->GetCocktailHeader(0));
     } 
     else if (fEvent->IsESD()) {       
       AliMCEvent *mcEvt = (AliMCEvent*)(fEvent->GetRefMC());
       if (!mcEvt) {
	 AliError("MC header branch in AOD not found!\n");
	 return kFALSE;
       }
       AliHeader * header = mcEvt->Header();
       dpmHeader = dynamic_cast<AliGenDPMjetEventHeader*>(header->GenEventHeader());
       if (!dpmHeader) {
	 AliError("Invalid DPMJET event header.");
	 return kFALSE;
       }
     } else 
       return kFALSE;
     
     Int_t nsd1 = 0;
     Int_t nsd2 = 0;
     Int_t ndd  = 0;
     dpmHeader->GetNDiffractive(nsd1, nsd2, ndd);
     if ( ((dpmHeader->ProjectileParticipants()==nsd1) && (ndd==0)) || 
	  ((dpmHeader->ProjectileParticipants()==nsd2) && (ndd==0))  ) {
       Printf("fFilterNSDeventsDPMJETpA2013 just rejected a non-NSD event!");
       return kFALSE;  
     }
   }
   
   return kTRUE;
}
